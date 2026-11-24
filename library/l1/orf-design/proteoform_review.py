"""Review every L1 reference that differs from its UniProt canonical: evidence for each form.

Becca asked on 2026-11-22 (comments on the brief): "are you sure it is just these 18 [named
exceptions], how many have been thoroughly checked", "have you checked all of the proteins with
major differences vs UniProt?" and "how sure are we that these differences are correct and the
better protein form?". The proteoform audit (``audit_l1_proteoforms.py``) flags canonical
residues the reference *lacks*; it never flags a reference that is much *longer* than the
canonical, and PeptideAtlas had been queried only for the MANE/APPRIS disagreements. This
script closes both gaps for every gene whose reference differs from the UniProt canonical.

Per gene (``changes_current_reference`` in the current reference record, the record
``library_pipeline.REFERENCE_RECORD`` names):

* the two proteins: ours (the designed CDS translated, ``library_design.CDS_KEY``) and the
  UniProt canonical (``uniprot_human_reviewed_isoforms.tsv.gz``, release pinned in the
  pipeline inputs), their lengths, the segments unique to each (``difflib`` opcodes of at least
  ``MIN_REGION`` residues, as ``select_reference.discriminating_regions``);
* DeepTMHMM v1 for each (signal peptide, TM count, N-terminal side), by exact sequence;
* PeptideAtlas build 607 peptides observed for *any* reviewed accession of the gene (canonical
  and every isoform: PeptideAtlas lists a peptide under the biosequence it maps to, so an
  isoform-only peptide never appears under the canonical), classed as distinguishing for one
  form when it is contained in that form and not the other. Unlike
  ``select_reference.peptides_supporting`` this counts peptides that span a deletion junction,
  so a form that is a pure deletion of the other can score. Proteotypic only (in no other
  gene's reviewed proteoform). The cache is the pinned one extended with the accessions it
  lacked, published under new dated keys (``ATLAS_KEY``);
* GENCODE Ribo-seq ORF translations overlapping each side's unique regions;
* the Ensembl 115 MANE Select transcript read from the GTF tag (``mane_select``), and whether
  the reference record saw it: the ladder reads MANE from the APPRIS file, which omits some;
* an evidence class: ``ours_supported`` / ``canonical_supported`` / ``both`` / ``neither``;
  ``canonical_strong_ours_none`` when a canonical-only peptide is seen in at least
  ``MIN_SAMPLES`` samples and ours has none, and ``canonical_strong_by_ladder_rule`` when it
  meets the ladder's own bar (``STRONG_EVIDENCE`` such peptides).

Writes ``proteoform_review_2026-11-22.tsv``, ``proteoform_review_2026-11-22_counts.json`` and
``proteoform_review_2026-11-22_mane_gap.tsv`` (every library gene whose MANE Select the record
missed).
The per-gene recommendations for the length-divergent genes and the flagged genes are in
``PROTEOFORM_REVIEW_2026-11-22.md``; this script records evidence, not decisions.

Implemented by Claude Code session l1-review-export-scope (subagent), 2026-11-22.

    uv run --locked python proteoform_review.py --fetch     # fetch missing accessions, publish
    uv run --locked python proteoform_review.py             # offline from the published cache
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import library_pipeline as lp  # noqa: E402

sys.path.insert(0, str(lp.AP))
import select_reference as sr  # noqa: E402
from select_principal import load_uniprot_proteoforms  # noqa: E402

from deliverome_analysis import library_design as ld  # noqa: E402
from deliverome_analysis import seqtools, topology  # noqa: E402

OUT = HERE / "proteoform_review_2026-11-22.tsv"
COUNTS = HERE / "proteoform_review_2026-11-22_counts.json"
MANE_GAP = HERE / "proteoform_review_2026-11-22_mane_gap.tsv"
#: The cache this review extends, and where the extension is published (append-only). The
#: 2026-11-22 key added the canonical accessions; ``-isoforms`` then added every other reviewed
#: isoform of the same genes, because PeptideAtlas returns a peptide under the biosequence it
#: maps to: a peptide unique to our isoform is listed under that isoform's accession, never
#: under the canonical's, so a canonical-only query cannot find evidence for our form.
ATLAS_BASE_KEY = (
    "peptideatlas/build-607/queries-2026-11-22-isoforms/peptideatlas_build607_pools.json"
)
ATLAS_BASE_SHA = "ae56bae3aa2df203cd63381193f636fdbdbe72fb89379aca668b02c8cc64516b"
#: ``-round2`` adds the accessions of the genes the 2026-11-22 round-2 additions brought in.
ATLAS_KEY = "peptideatlas/build-607/queries-2026-11-22-round2/peptideatlas_build607_pools.json"
#: The record whose differing genes were "the 206" in the brief.
RECORD_206 = "library-design/reference-selection/2026-11-12-reference-fixes-2/"
#: A reference this much longer or shorter than the canonical is a "major" length difference.
LENGTH_FRACTION = 0.10

FIELDS = [
    "gene_symbol",
    "in_206",
    "decision",
    "selected_source",
    "canonical_acc",
    "ours_uniprot_equivalent",
    "mane_enst",
    "appris_enst",
    "mane_select_gtf",
    "mane_select_is",
    "mane_missing_from_record",
    "ours_aa",
    "canonical_aa",
    "length_diff",
    "length_diff_fraction",
    "major_length_difference",
    "ours_unique_residues",
    "canonical_unique_residues",
    "ours_unique_segments",
    "canonical_unique_segments",
    "ours_v1",
    "canonical_v1",
    "topology_differs",
    "atlas_peptides_shared",
    "atlas_ours_only",
    "atlas_ours_best_samples",
    "atlas_ours_ge_min_samples",
    "atlas_ours_examples",
    "atlas_canonical_only",
    "atlas_canonical_best_samples",
    "atlas_canonical_ge_min_samples",
    "atlas_canonical_examples",
    "riboseq_ours",
    "riboseq_canonical",
    "evidence_class",
    "canonical_strong_ours_none",
    "canonical_strong_by_ladder_rule",
    "recommendation",
    "recommendation_basis",
]

#: The 11 references more than 10% longer or shorter than the canonical that no earlier pass
#: reviewed, and the genes the PeptideAtlas pass flags. A recommendation, not a decision:
#: every "switch" changes the library and is Becca's call. Others take ``default_recommendation``.
UNREVIEWED_11 = "KCNN2 PLPPR2 PCDH7 PCDH1 HRH2 KCNC1 GPM6B TSPAN17 TMEM145 LY6H GFRA4".split()
RECOMMEND = {
    "KCNN2": (
        "keep",
        "MANE and APPRIS agree (ENST00000673685, no UniProt isoform). The 212-aa N-terminal extension is cytoplasmic (v1: N-in, 6 TM in both); every extracellular residue is shared. No distinguishing peptide either way. The extended form is a described, functional, surface-expressed SK2 isoform (SK2-L, mouse; Strassmaier et al. 2005, doi:10.1074/jbc.M413125200)",
    ),
    "PLPPR2": (
        "keep",
        "MANE and APPRIS agree (ENST00000688289). Our 452-aa form has its own C-terminal tail (cytoplasmic; v1 6 TM N-in in both) and 9 distinguishing peptides (best 15 samples) vs 3 for the canonical (best 8)",
    ),
    "PCDH7": (
        "keep",
        "MANE (ENST00000695919; APPRIS calls it ALTERNATIVE:1 and prefers ENST00000511884, an Ensembl-only 1,247-aa form). Differs only in the cytoplasmic C-terminal tail; extracellular cadherin repeats identical. One canonical-only peptide in 5 samples, none for ours: below the ladder's bar, so the MANE default stands",
    ),
    "PCDH1": (
        "keep",
        "MANE and APPRIS agree (ENST00000287008 = Q08174-2). 93 distinguishing peptides for our longer cytoplasmic tail (best 377 samples) vs 14 for the canonical (best 17)",
    ),
    "HRH2": (
        "keep, low confidence",
        "MANE and APPRIS agree (ENST00000636584, no UniProt isoform). Ours adds 63 aa at the cytoplasmic C-terminus of the 359-aa receptor; v1 7 TM N-out in both, ligand-binding and extracellular loops identical. No distinguishing peptide and no Ribo-seq ORF for either form; the extension rests on the transcript model alone",
    ),
    "KCNC1": (
        "keep",
        "MANE and APPRIS agree (ENST00000265969 = P48547-2). Ours carries the longer cytoplasmic C-terminus; 10 distinguishing peptides (best 59 samples) vs none for the canonical",
    ),
    "GPM6B": (
        "keep",
        "MANE and APPRIS agree (ENST00000316715 = Q13491-4). The 40-aa N-terminal extension is cytoplasmic (v1 4 TM N-in in both). 47 distinguishing peptides (best 247 samples) vs 6 for the canonical (best 58): both forms exist, ours is the better observed",
    ),
    "TSPAN17": (
        "keep",
        "MANE and APPRIS agree (ENST00000508164 = Q96FV3-4). The 59-aa extension is the cytoplasmic C-terminus; 3 distinguishing peptides (best 3 samples), none for the canonical",
    ),
    "TMEM145": (
        "switch to Q8NBT3 (canonical), ENST00000301204",
        "The ladder called this APPRIS-only because the APPRIS file does not list the MANE Select transcript; Ensembl 115 tags ENST00000301204 MANE Select and it encodes the 493-aa canonical exactly. MANE and APPRIS therefore disagree, no peptide or Ribo-seq evidence either way, and the policy default for a disagreement is MANE. Ours (ENST00000673187, 563 aa) adds 70 aa in the cytoplasmic C-terminus and a 14-aa insertion; v1 7 TM N-out in both",
    ),
    "LY6H": (
        "switch to O94772 (canonical), ENST00000430474",
        "Ours (MANE ENST00000342752 = O94772-2, APPRIS ALTERNATIVE:1) starts 21 aa upstream of the signal peptide; the canonical, which is the APPRIS principal, starts at it. The mature GPI-anchored protein is identical. Same case as RGMB (upstream_extension, switched 2026-11-12): an extension ahead of the signal peptide risks targeting and buys nothing. No distinguishing peptide either way",
    ),
    "GFRA4": (
        "keep, low confidence",
        "MANE and APPRIS agree (ENST00000290417 = Q9GZZ7-2). Differs from the canonical inside the second GDNF-binding domain (canonical residues 132-196); both ends, the signal peptide and the C-terminal GPI signal are shared (v1 SP, 0 TM in both). One canonical-only peptide in 1 sample; no evidence for ours. 10.03% shorter, just over the threshold",
    ),
    "ACHE": (
        "keep",
        "Named GPI form (P22303-2, H subunit; reference_overrides.tsv). Canonical evidence is expected: the canonical T subunit is the abundant secreted form",
    ),
    "BCAN": (
        "keep",
        "Named GPI form (Q96GW7-2; reference_overrides.tsv). The canonical is secreted brevican, the abundant form",
    ),
    "NCAM1": ("keep", "Named GPI form (P13591-3, NCAM-120; reference_overrides.tsv)"),
    "TFPI": ("keep", "Named GPI form (P10646-2, TFPI-beta; reference_overrides.tsv)"),
    "MUC1": (
        "keep, flag to Becca",
        "APPRIS principal ENST00000337604 (P15941-8, 273 aa) lacks the canonical's 980-residue tandem-repeat region; the canonical has 29 distinguishing peptides (best 97 samples). The repeat array cannot be synthesized as a pooled ORF, so the canonical is not a buildable alternative; the short form keeps the SEA domain, TM and cytoplasmic tail. Ensembl 115 tags no MANE Select for MUC1",
    ),
    "PTPRM": (
        "switch to P28827 (canonical), ENST00000332175",
        "MANE ENST00000580170 (P28827-2, 1,465 aa) and the APPRIS principal ENST00000332175 (P28827, 1,452 aa) disagree. The ladder found no evidence because the APPRIS form is a pure deletion and its test only looks inside residues unique to one side; counting junction-spanning peptides, the canonical has 7 distinguishing peptides (4 in >= 10 samples, best 28) and ours none, which meets the ladder's own bar for overruling MANE",
    ),
    "SLC45A4": (
        "revisit (topology differs)",
        "Ladder missed MANE Select ENST00000517878 (808 aa, no UniProt isoform; neither ours nor the canonical). Ours (Q5BKX6-3, 798 aa) has 0 distinguishing peptides, the canonical 10 (5 in >= 10 samples), and v1 places their N-termini on opposite sides (ours N-out, canonical N-in, 12 TM each). Needs a three-way look at MANE, ours and the canonical",
    ),
    "GPNMB": (
        "revisit",
        "Ladder missed MANE Select ENST00000258733 (Q14956-2, 560 aa). Ours (588 aa) has no distinguishing peptide; the canonical (572 aa) has 2 (best 16 samples). Under the policy a MANE/APPRIS disagreement with weak evidence defaults to MANE, i.e. Q14956-2",
    ),
}


def default_recommendation(r: dict) -> tuple[str, str]:
    """Rule for genes outside :data:`RECOMMEND`, applying the ladder's own policy."""
    if r["selected_source"] == "named":
        return "keep", "named reference (reference_overrides.tsv)"
    if r["mane_missing_from_record"] and r["mane_select_is"] != "ours":
        tgt = f"MANE Select {r['mane_select_gtf']}"
        if r["mane_select_is"] == "canonical":
            tgt += f" (= canonical {r['canonical_acc']})"
        verb = "switch to" if r["mane_select_is"] == "canonical" else "revisit: policy points to"
        if r["evidence_class"] in ("ours_supported", "both") and r["atlas_ours_ge_min_samples"]:
            verb = "revisit (ours also observed): policy points to"
        return f"{verb} {tgt}", (
            "the ladder read MANE from the APPRIS file, which omits this MANE Select transcript, "
            "so it called the gene APPRIS-only; MANE and APPRIS disagree and the policy default "
            f"is MANE unless peptides strongly favor APPRIS (ours {r['atlas_ours_only']} "
            f"distinguishing, {r['atlas_ours_ge_min_samples']} in >= 10 samples; canonical "
            f"{r['atlas_canonical_only']}, {r['atlas_canonical_ge_min_samples']})"
        )
    if r["canonical_strong_by_ladder_rule"]:
        return "revisit", "canonical meets the ladder's evidence bar and ours has none"
    return "keep", f"no evidence against the reference ({r['evidence_class']})"


def mane_select(gtf: Path, pep_fa: Path) -> dict[str, tuple[str, str]]:
    """``{gene_name: (ENST, protein)}`` for every Ensembl 115 transcript tagged MANE_Select.

    The reference ladder reads MANE from the APPRIS file's MANE column, and APPRIS e115v50
    does not list every MANE Select transcript, so a gene can look "APPRIS only" while
    Ensembl tags a MANE Select transcript for it. This reads the tag from the GTF itself.
    """
    import gzip

    tx: dict[str, str] = {}
    with gzip.open(gtf, "rt") as fh:
        for line in fh:
            if "\ttranscript\t" not in line or '"MANE_Select"' not in line:
                continue
            name = line.split('gene_name "')[1].split('"')[0] if "gene_name" in line else ""
            tx[name] = line.split('transcript_id "')[1].split('"')[0]
    seqs: dict[str, list[str]] = {}
    cur = ""
    with gzip.open(pep_fa, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                cur = next(f for f in line.split() if f.startswith("transcript:"))
                cur = cur.split(":")[1].split(".")[0]
                seqs[cur] = []
            else:
                seqs[cur].append(line.strip())
    return {g: (t, "".join(seqs.get(t, [])).rstrip("*")) for g, t in tx.items()}


def v1(seq: str) -> dict | None:
    hit = topology.for_sequence(seq)
    if hit is None:
        return None
    topo = hit["per_residue_topology"]
    sp = topo.count("S")
    first = next((c for c in topo if c in "IO"), "")
    return {
        "label": hit["deeptmhmm_label"],
        "sp": sp,
        "tm": int(hit["tm_helix_count"] or 0),
        "nterm": {"O": "out", "I": "in"}.get(first, "?"),
    }


def v1_text(d: dict | None) -> str:
    if d is None:
        return "no v1 row"
    return f"{d['label']} SP{d['sp']} TM{d['tm']} N-{d['nterm']}"


def distinguishing(
    one: str, other: str, observed: dict[str, tuple[int, int]], proteotypic
) -> list[tuple[str, int, int]]:
    """Observed proteotypic peptides contained in ``one`` and absent from ``other``.

    Unlike ``select_reference.peptides_supporting``, which tests only substrings of the
    regions unique to one form, this tests every observed peptide, so a peptide spanning
    the junction of a deletion counts for the shorter form. Without that a form that is a
    pure deletion of the other (no residues of its own) could never score any evidence.
    """
    hits = [
        (pep, v[0], v[1])
        for pep, v in observed.items()
        if sr.MIN_PEP <= len(pep) <= sr.MAX_PEP and pep in one and pep not in other
    ]
    hits = [h for h in hits if proteotypic(h[0])]
    return sorted(hits, key=lambda x: (-x[2], -x[1], -len(x[0]), x[0]))


def segments(regions: list[tuple[int, str]]) -> str:
    return ";".join(f"{s}-{s + len(r) - 1}" for s, r in regions)


def _review_row(
    r: dict, ours: str, can: str, obs: dict, prot, in_206: set, orf_seqs: set, mane_pair: tuple
) -> dict:
    """One reviewed gene's row: our reference against the canonical, with its evidence."""
    g = r["symbol"]
    only_ours, only_can = sr.discriminating_regions(ours, can)
    sup_o = distinguishing(ours, can, obs, prot)
    sup_c = distinguishing(can, ours, obs, prot)
    shared = sum(1 for pep in obs if pep in ours and pep in can)
    v_o, v_c = v1(ours), v1(can)
    m_enst, m_seq = mane_pair
    diff = len(ours) - len(can)
    frac = abs(diff) / len(can)
    best_o = max((x[2] for x in sup_o), default=0)
    best_c = max((x[2] for x in sup_c), default=0)
    cls = (
        "both"
        if sup_o and sup_c
        else "ours_supported"
        if sup_o
        else "canonical_supported"
        if sup_c
        else "neither"
    )
    return {
        "gene_symbol": g,
        "in_206": int(g in in_206),
        "decision": r["decision"],
        "selected_source": r["selected_source"],
        "canonical_acc": r["current_accession"],
        "ours_uniprot_equivalent": r["selected_uniprot_equivalent"],
        "mane_enst": r["mane_enst"],
        "appris_enst": r["appris_enst"],
        "mane_select_gtf": m_enst,
        "mane_select_is": (
            ""
            if not m_enst
            else "ours"
            if m_seq == ours
            else "canonical"
            if m_seq == can
            else "neither"
        ),
        "mane_missing_from_record": int(bool(m_enst) and not r["mane_enst"]),
        "ours_aa": len(ours),
        "canonical_aa": len(can),
        "length_diff": diff,
        "length_diff_fraction": round(frac, 3),
        "major_length_difference": int(frac > LENGTH_FRACTION),
        "ours_unique_residues": sum(len(s) for _, s in only_ours),
        "canonical_unique_residues": sum(len(s) for _, s in only_can),
        "ours_unique_segments": segments(only_ours),
        "canonical_unique_segments": segments(only_can),
        "ours_v1": v1_text(v_o),
        "canonical_v1": v1_text(v_c),
        "topology_differs": int(
            v_o is not None
            and v_c is not None
            and (v_o["sp"] > 0, v_o["tm"], v_o["nterm"]) != (v_c["sp"] > 0, v_c["tm"], v_c["nterm"])
        ),
        "atlas_peptides_shared": shared,
        "atlas_ours_only": len(sup_o),
        "atlas_ours_best_samples": best_o,
        "atlas_ours_ge_min_samples": sum(1 for x in sup_o if x[2] >= sr.MIN_SAMPLES),
        "atlas_ours_examples": ";".join(x[0] for x in sup_o[:3]),
        "atlas_canonical_only": len(sup_c),
        "atlas_canonical_best_samples": best_c,
        "atlas_canonical_ge_min_samples": sum(1 for x in sup_c if x[2] >= sr.MIN_SAMPLES),
        "atlas_canonical_examples": ";".join(x[0] for x in sup_c[:3]),
        "riboseq_ours": len(sr.riboseq_overlap(only_ours, orf_seqs)),
        "riboseq_canonical": len(sr.riboseq_overlap(only_can, orf_seqs)),
        "evidence_class": cls,
        "canonical_strong_ours_none": int(best_c >= sr.MIN_SAMPLES and not sup_o),
        # The ladder's own bar for overruling MANE: >= STRONG_EVIDENCE peptides seen in
        # >= MIN_SAMPLES samples on one side and none on the other.
        "canonical_strong_by_ladder_rule": int(
            sum(1 for x in sup_c if x[2] >= sr.MIN_SAMPLES) >= sr.STRONG_EVIDENCE and not sup_o
        ),
    }


def _mane_gap_rows(refsel: list[dict], mane: dict, cds: dict, uniprot: dict) -> list[dict]:
    """Library genes whose reference record lacks the MANE Select the GTF names."""
    gap = []
    for r in refsel:
        g = r["symbol"]
        if r["mane_enst"] or g not in mane or f"lib:{g}" not in cds:
            continue
        ours = seqtools.translate(cds[f"lib:{g}"]["cds"]).rstrip("*")
        m_enst, m_seq = mane[g]
        can = next(
            (p["seq"] for p in uniprot.get(g, []) if p["accession"] == r["current_accession"]), ""
        )
        gap.append(
            {
                "gene_symbol": g,
                "decision": r["decision"],
                "reference_rung": r["reference_rung"],
                "selected_source": r["selected_source"],
                "appris_enst": r["appris_enst"],
                "mane_select_gtf": m_enst,
                "mane_aa": len(m_seq),
                "ours_aa": len(ours),
                "canonical_acc": r["current_accession"],
                "canonical_aa": len(can),
                "ours_is_mane_protein": int(m_seq == ours),
                "mane_is_canonical_protein": int(bool(can) and m_seq == can),
            }
        )
    return gap


def _counts(
    changed: list, no_cds: list, rows: list[dict], accessions: list, gap: list[dict]
) -> dict:
    """The summary counts written beside the review table."""
    from collections import Counter

    counts = {
        "reference_record": lp.REFERENCE_RECORD,
        "atlas_cache": f"external://{ATLAS_KEY}",
        "genes_differing_from_canonical": len(changed),
        "not_in_library_no_designed_cds": sorted(no_cds),
        "reviewed": len(rows),
        "reviewed_in_206": sum(r["in_206"] for r in rows),
        "accessions_queried": len(accessions),
        "genes_with_any_atlas_peptide": sum(
            1
            for r in rows
            if r["atlas_peptides_shared"] + r["atlas_ours_only"] + r["atlas_canonical_only"] > 0
        ),
        "evidence_class": dict(Counter(r["evidence_class"] for r in rows)),
        "evidence_class_in_206": dict(Counter(r["evidence_class"] for r in rows if r["in_206"])),
        "canonical_strong_ours_none": sum(r["canonical_strong_ours_none"] for r in rows),
        "canonical_strong_by_ladder_rule": sum(r["canonical_strong_by_ladder_rule"] for r in rows),
        "canonical_strong_by_ladder_rule_genes": [
            r["gene_symbol"] for r in rows if r["canonical_strong_by_ladder_rule"]
        ],
        "canonical_strong_ours_none_genes": [
            r["gene_symbol"] for r in rows if r["canonical_strong_ours_none"]
        ],
        "major_length_difference": sum(r["major_length_difference"] for r in rows),
        "major_length_difference_genes": [
            r["gene_symbol"] for r in rows if r["major_length_difference"]
        ],
        "library_genes_mane_select_missing_from_record": len(gap),
        "of_which_ours_differs_from_mane_protein": sum(
            1 for x in gap if not x["ours_is_mane_protein"]
        ),
        "of_which_mane_is_canonical_protein": sum(
            1 for x in gap if not x["ours_is_mane_protein"] and x["mane_is_canonical_protein"]
        ),
        "recommendation": dict(Counter(r["recommendation"].split()[0].rstrip(",:") for r in rows)),
        "switch_genes": sorted(
            r["gene_symbol"] for r in rows if r["recommendation"].startswith("switch")
        ),
        "revisit_genes": sorted(
            r["gene_symbol"] for r in rows if r["recommendation"].startswith("revisit")
        ),
        "unreviewed_11": UNREVIEWED_11,
        "topology_differs": sum(r["topology_differs"] for r in rows),
        "topology_differs_genes": [r["gene_symbol"] for r in rows if r["topology_differs"]],
    }
    return counts


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--fetch", action="store_true", help="query missing accessions and publish")
    a = ap.parse_args(argv)
    from deliverome_analysis.s3 import external, processed

    work = HERE / "data" / "proteoform_review"  # gitignored
    work.mkdir(parents=True, exist_ok=True)
    refsel = list(
        csv.DictReader(
            open(processed.download(lp.PREFIX + lp.REFERENCE_RECORD + "reference_selection.tsv")),
            delimiter="\t",
        )
    )
    changed = [r for r in refsel if r["changes_current_reference"] == "True"]
    in_206 = {
        r["symbol"]
        for r in csv.DictReader(
            open(processed.download(lp.PREFIX + RECORD_206 + "reference_selection.tsv")),
            delimiter="\t",
        )
        if r["changes_current_reference"] == "True"
    }
    no_cds = []
    cds = json.loads(Path(processed.download(ld.CDS_KEY)).read_text())
    got = lp.stage_inputs(
        [
            "uniprot_human_reviewed_isoforms.tsv.gz",
            "Homo_sapiens.GRCh38.115.gtf.gz",
            "Homo_sapiens.GRCh38.pep.all.fa.gz",
            "Ribo-seq_ORFs.all.bed",
            "Ribo-seq_ORFs.all.faa",
            "Ribo-seq_ORFs.primary.bed",
            "Ribo-seq_ORFs.primary.faa",
            "Ribo-seq_ORFs.comprehensive.bed",
            "Ribo-seq_ORFs.comprehensive.faa",
        ],
        work / "in",
    )
    mane = mane_select(
        got["Homo_sapiens.GRCh38.115.gtf.gz"], got["Homo_sapiens.GRCh38.pep.all.fa.gz"]
    )
    uniprot = load_uniprot_proteoforms(got["uniprot_human_reviewed_isoforms.tsv.gz"])
    ribo_dir = work / "ribo"
    ribo_dir.mkdir(exist_ok=True)
    for name, path in got.items():
        if name.startswith("Ribo-seq"):
            shutil.copy(path, ribo_dir / name)
    _types, ribo_seqs = sr.load_riboseq(ribo_dir)

    # PeptideAtlas: the published cache, extended once with the accessions it lacks.
    cache = work / "atlas"
    cache.mkdir(exist_ok=True)
    store = cache / "peptideatlas_build607_pools.json"
    if a.fetch:
        base = Path(external.download(ATLAS_BASE_KEY))
        assert hashlib.sha256(base.read_bytes()).hexdigest() == ATLAS_BASE_SHA
        shutil.copy(base, store)
        old = json.loads(store.read_text())
    else:
        shutil.copy(external.download(ATLAS_KEY), store)
    gene_accs = {
        r["symbol"]: sorted(
            {r["current_accession"]} | {p["accession"] for p in uniprot.get(r["symbol"], [])}
        )
        for r in changed
        if r["current_accession"]
    }
    accessions = sorted({acc for accs in gene_accs.values() for acc in accs})
    observed = sr.fetch_peptides(accessions, cache, offline=not a.fetch)
    if a.fetch:
        new = json.loads(store.read_text())
        assert all(new["pools"][k] == v for k, v in old["pools"].items()), "a pinned pool changed"
        assert all(new["index"][k] == v for k, v in old["index"].items()), "a pinned index moved"
        if external.exists(ATLAS_KEY):
            raise SystemExit(f"external://{ATLAS_KEY} exists; bump the date")
        external.upload(store, ATLAS_KEY)
        print(f"published external://{ATLAS_KEY} ({len(old['index'])} -> {len(new['index'])})")

    all_proteoforms = [(sym, p["seq"]) for sym, recs in uniprot.items() for p in recs]

    def proteotypic_for(symbol: str):
        return lambda sub: not any(sub in s for sym, s in all_proteoforms if sym != symbol)

    rows = []
    for r in sorted(changed, key=lambda x: x["symbol"]):
        g = r["symbol"]
        key = f"lib:{g}"
        if key not in cds:
            no_cds.append(g)
            continue
        ours = seqtools.translate(cds[key]["cds"]).rstrip("*")
        can = next(
            (p["seq"] for p in uniprot.get(g, []) if p["accession"] == r["current_accession"]), ""
        )
        if not can:
            continue
        obs: dict[str, tuple[int, int]] = {}
        for acc in gene_accs[g]:
            for pep, v in observed.get(acc, {}).items():
                if pep not in obs or v[1] > obs[pep][1]:
                    obs[pep] = v
        rows.append(
            _review_row(
                r,
                ours,
                can,
                obs,
                proteotypic_for(g),
                in_206,
                ribo_seqs.get(g, set()),
                mane.get(g, ("", "")),
            )
        )
    for r in rows:
        r["recommendation"], r["recommendation_basis"] = RECOMMEND.get(
            r["gene_symbol"], default_recommendation(r)
        )
    assert all(g in {r["gene_symbol"] for r in rows} for g in RECOMMEND), "a reviewed gene is gone"
    with OUT.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    gap = _mane_gap_rows(refsel, mane, cds, uniprot)
    with MANE_GAP.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(gap[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(gap)

    counts = _counts(changed, no_cds, rows, accessions, gap)
    COUNTS.write_text(json.dumps(counts, indent=1) + "\n")
    print(json.dumps(counts, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
