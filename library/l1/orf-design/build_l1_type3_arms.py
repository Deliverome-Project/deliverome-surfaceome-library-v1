"""Build the tagged inserts for the L1 "N at residue 1, single-pass (type III)" rung.

Every member of the rung gets one of three designs, assigned in ``design_for`` from the
subgroups in ``type3_nterm_calls.tsv`` (see TYPE3_NTERM.md):

* ``igk_tag`` -- Igk leader + Asp + the member's Met1 + terminal cassette + native 2..end.
* ``tag_alone`` -- Met1 + terminal cassette + native 2..end.
* ``held`` -- not built; the reason is recorded.

**Decisions and provenance** (also written into every output row):

| Decision | Decided by | Implemented by |
|---|---|---|
| Option B: Igk + tag for the 18 members with an identified risk, tag alone for the other 33 genuine type III members | Becca, 2026-11-11 ("go with B") | Claude Code session l1-nterm-type3 (PR #282) |
| Igk + tag for the five zero-TM GPI members (BST1, CDH13, OPCML, LYPD5 isoforms; cyno ACHE) | Becca, 2026-11-10 ("fix these 5") | same |
| Full 51-aa cassette at residue 1, leading GGGGS kept | Becca, 2026-11-11 ("keep full") | same |
| Second review applied: OPALIN and SIGIRR to tag alone; six isoforms held for redraw; CD200R1L, LRRC37B and CACHD1 held for a swap to their signal-peptide forms; cyno CEACAM8 held for gene-model repair | Becca, 2026-11-11 ("Apply all"); review by session type3-risk-review | same |
| The five GPI members sent for redraw (superseding Igk + tag) | Becca, 2026-11-11 ("Redraw them too") | same |
| Held: LMLN, the three members whose annotated signal peptide is uncalled, and cyno TNFRSF1B / TIE1 pending CDS repair | no decision taken (open questions in TYPE3_NTERM.md) | same |

**Cassette.** The full L1 cassette, GGGGS | SpyTag003 | GGGGSGGGGS | P-ALFA-P | GGGGS (51 aa),
goes after Met1, with p003's SpyTag003 and P-ALFA-P codons and this arm's own GGGGS encoding
(``l1_arm_cassette_encodings.tsv``; Becca, 2026-11-21: "we still don't want shared linkers in the
middle"). Dropping the leading GGGGS
was tested and rejected. With SpyTag003 directly after Met1, DeepTMHMM v1 called a cleavable
signal peptide ending inside SpyTag003 for 12 of the 62 tag-alone arms. The calls fell with
spacer length: S 7, G 7, GS 3, GG 2, GGS 0, GGGGS 0 (TYPE3_NTERM.md).

**Option B membership** (Becca, 2026-11-11). Igk + tag where a risk was identified:

* the ordered form lacks the gene's signal peptide (subgroup ``sp_missing``), except the two
  members held for CDS repair;
* a native N-tail of <= 20 residues (``short_tail``: RELL2, TEX38);
* KLB, whose tag-alone arm v1 calls with the cassette cytosolic;
* mid-tail members with a native flank-charge difference Delta <= 1 (OPALIN, TEX29, SIGIRR,
  TNFRSF17 human and cyno). The <= 1 cut is the analysing agent's choice, not a published
  threshold.

Tag alone: ``long_tail`` (except KLB), ``dual_topology`` (MRAP, MRAP2), and ``mid_tail``
with Delta >= 2.

**Layout** (the L1 Esp3I wrapper, as ``build_l1_cterm_arms.py`` and the GPCR arms)::

    CGTCTCA | Kozak | [Igk leader (60 nt) | GAC] | ATG | cassette (51 aa) | ORF codons 2..end
            incl. native stop | CC | T7 | 5' handle | barcode (20 N) | Feldman 3' | 10x CS1 | AGA GAGACG

Per member the build asserts the encoded protein exactly, one stop at the end, and that the
insert carries only the two wrapper Esp3I sites (no other BbsI/Esp3I/BsaI), alone and in the
p0039 backbone.

    uv run --locked python build_l1_type3_arms.py
    uv run --locked python build_l1_type3_arms.py --check
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path

from build_constructs import CASSETTE_AA
from diversify_arm_cassettes import arm_cassette

from deliverome_analysis.cloning import library_design as ld
from deliverome_analysis.cloning import plasmid_qc as qc
from deliverome_analysis.cloning import seqtools

LEADER = qc.IGK_LEADER_NT
LEADER_AA = qc.IGK_LEADER_AA
PLUS1 = qc.IGK_PLUS1_NT
PLUS1_AA = qc.IGK_PLUS1_AA

HERE = Path(__file__).parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO / "library/l1/proteoform-selection"))

MANIFEST = HERE / "candidate_manifest_tranches.tsv"
TAG_PATHS = HERE / "l1_tag_paths.tsv"
CALLS = HERE / "type3_nterm_calls.tsv"
OUT_TSV = HERE / "constructs_l1_type3.tsv"
OUT_FASTA = HERE / "constructs_l1_type3.fasta"
#: Protein-level view of the built designs, for refolding: the translated precursor, the
#: mature chain, the native sequence (whose MSA a refold can reuse with gap columns for the
#: inserted residues), and where the leader, the cassette and native TM1 sit.
OUT_PROTEIN_TSV = HERE / "constructs_l1_type3_proteins.tsv"
OUT_PROTEIN_FASTA = HERE / "constructs_l1_type3_proteins.fasta"
#: One row per member: the design, who decided it, why, and the alternative rejected.
OUT_DECISIONS = HERE / "type3_decisions.tsv"
FEATURES = HERE / "type3_nterm_features.tsv"
#: Dispositions for refold flags, decided under the delegation (design_id -> (disposition,
#: reason, rejected alternative)). A flag that the evidence shows is not a topology failure
#: keeps its design; one that is would change the design here.
REFOLD_DISPOSITIONS = {
    pid: ("keep design",
          "refold flag 'confident SpyTag003 helix in the cassette' (AF2 model_1_ptm): the helix makes no contact with the native protein (cassette_contact_res 0), TM1 stays one confident helix, the native domain is intact (lDDT >= 0.96 vs the native model), DeepTMHMM v1 calls no TM or signal peptide in the cassette, and passing members carry SpyTag003 helices of the same length (10-11). Free SpyTag peptides are disordered and bind SpyCatcher as a strand; a predicted helix in the unbound tag is not evidence of a topology failure. KCNE3 is a pilot pair, which tests display directly",
          "Igk + tag: rejected, the cassette follows Met1 in both designs so the leader does not change the tag's local conformation")
    for pid in ("l1-p1002", "l1-p3087", "l1-p1003", "l1-p2053")  # KCNE3 human, mouse; KCNE4; SMAGP
}  # fmt: skip
REFOLD_DISPOSITIONS["l1-p0906"] = (
    "keep design",
    "IGSF5 refold flag 'native TM1 not a single helix (longest helical run 9)': the TM1 sequence VGLGLAGTMLLTPTCTLTI carries a mid-helix Pro and is Gly/Thr-rich (Kyte-Doolittle 1.45, the weakest on the rung), the drop in confidence follows that sequence, the cassette sits 266 residues upstream, the native domain matches its model (lDDT 0.98), and DeepTMHMM v1 calls the TM in both native and tagged forms; the flag reflects the native TM in a one-model, three-recycle fold. IGSF5's N-terminus (no signal peptide by SignalP 6 or v1, a 38-aa SP by v2, a Met-led hydrophobic stretch near residue 70) is referred to the proteoform audit",
    "Igk + tag: rejected, the leader does not change the TM sequence; dropping the member: rejected, no evidence the TM fails in cells",
)

#: Refold verdicts from session l1-tag-refold (AF2 model_1_ptm, 3 recycles, no relax; tag_alone
#: folds the precursor, igk_tag the mature chain), keyed by design_id. Read when present.
REFOLD_KEY = (
    "sprint4/s4e9-orf-design-fg-panel/predicted-structures/tagged/2026-11-12/verdicts_type3.tsv"
)
#: The domesticated ORF record, versioned with the CDS record (library_design.CONSTRUCTS_KEY).
CONSTRUCTS_KEY = ld.CONSTRUCTS_KEY
RUNG = "N at residue 1, single-pass (type III)"
BARCODE_SLOT = qc.BARCODE_SLOT
ENZYMES = list(qc.DOMESTICATION_ENZYMES)


IMPLEMENTED_BY = "Claude Code session l1-nterm-type3 (PR #282)"
B_DECIDED = "Becca, 2026-11-11 (option B)"
GPI_DECIDED = "Becca, 2026-11-10"
CASSETTE_DECIDED = "Becca, 2026-11-11 (keep the full 51-aa cassette)"
NO_DECISION = "no decision taken"

#: Zero-TM GPI members built with the leader (Becca, 2026-11-10).
GPI_IGK = {"l1-p2401", "l1-p2417", "l1-p2476", "l1-p2727", "l1-p2531"}
#: Held for CDS repair: cyno TNFRSF1B (Met110 begins the conserved signal peptide) and cyno
#: TIE1 (no ATG start; library_cds_issues.tsv CDS101).
CDS_REPAIR = {"l1-p2623": "CDS starts 109 codons upstream of Met110, which begins the conserved signal peptide",
              "l1-p2618": "ORF has no ATG start (library_cds_issues.tsv CDS101)"}  # fmt: skip

DELEGATED_DECIDED = 'l1-nterm-type3 (delegated by Becca 2026-11-12, "use your best judgment")'
#: Protein-level calls this session made under Becca's 2026-11-12 delegation, keyed by
#: design_id: (design, reason with evidence, rejected alternative). Each is reversible here.
#: SignalP 6 calls come from the signalp6-runs Modal volume, run sp6_2026_09_28
#: (prediction_results.txt, by UniProt accession), whose input was the exact sequence ordered.
DELEGATED = {
    "l1-p0181": ("tag_alone",
                 "BMPR1B: the ordered form is canonical O00238; SignalP 6 calls no signal peptide (OTHER, P(SP) 0.000), so UniProt's 1-13 annotation is not supported and TM1 is the signal-anchor. Mid N-tail 126, Delta +3: option B gives the tag alone",
                 "N-after-signal-peptide: rejected, no signal peptide predicted"),
    "l1-p0201": ("tag_alone",
                 "C11orf87: the ordered form is canonical Q6NUJ2; SignalP 6 calls no signal peptide (OTHER, P(SP) 0.13). Mid N-tail 61, Delta +8: option B gives the tag alone",
                 "N-after-signal-peptide: rejected, no signal peptide predicted"),
    "l1-p2641": ("tag_alone",
                 "mouse C11orf87: no SignalP 6 record; its first 24 residues, the region annotated as signal peptide in human, are identical to human C11orf87, which SignalP 6 calls OTHER (P(SP) 0.13). N-tail 60, Delta +8: tag alone, matching human",
                 "N-after-signal-peptide: rejected by the human call; a mouse-specific SignalP run would firm it up"),
    "l1-p3048": ("tag_alone",
                 "mouse SIGIRR: the Delta <= 1 rule gives Igk + tag, but human SIGIRR (l1-p1809) went to tag alone because its type III topology is established (UniProt Q6IA17); the orthologs share the topology",
                 "Igk + tag: rejected, the risk flag is the weak-Delta rule the second review showed does not mark cassette risk"),
    "l1-p0906": ("igk_tag",
                 "IGSF5: human IGSF5 has no signal peptide (SignalP 6 P(SP) 0.000; the proteoform audit confirms Met1 and finds the loss is divergence, while mouse and cyno IGSF5 keep one), so as ordered a folded 266-residue Ig ectodomain would have to be translocated by TM1 acting as a signal-anchor, which folded N-domains block (Denzer et al. 1995), and TM1 is the weakest on the rung (Kyte-Doolittle 1.45, mid-helix Pro; it folds as a broken helix in the tag-alone refold). The leader gives the ectodomain co-translational entry, as the orthologs' own signal peptides do; v1 keeps the cassette outside and TM1 as the single TM in the Igk arm",
                 "tag alone (option B's long-tail call): rejected, it rests on unassisted translocation of a folded Ig domain; replacing residues 1-34 with the leader (mature start at S35, the ortholog-shared N-terminus): rejected, it deletes human sequence and departs from the shared Igk layout"),
    "l1-p1106": ("held",
                 "LMLN: the L1 reference (ENST00000482695) starts at canonical residue 233, past the signal peptide that SignalP 6 calls on canonical Q96KR4 (SP, P 0.83, cleavage after 43); as ordered it never enters the ER. Swap the reference to canonical Q96KR4 (ENST00000330198), which keeps that signal peptide and the GPI anchor, on the N-after-SP rung; carried out by the manifest writer on PR #283, which preferred the canonical over Q96KR4-3 (l2-p0894, a 37-aa insertion with no supporting evidence) and keeps Q96KR4-3 only as a fallback",
                 "Igk + tag on the truncated form: rejected, it would display a fragment no cell makes, as for the redrawn GPI isoforms"),
}  # fmt: skip

GPI_REDRAW_DECIDED = (
    'Becca, 2026-11-11 ("Redraw them too"), superseding the 2026-11-10 Igk + tag decision'
)
#: The five zero-TM GPI members. Built with Igk + tag on 2026-11-10, then sent for redraw on
#: 2026-11-11 by the second review's surface-form test: each ordered form lacks its gene's
#: signal peptide, so with a leader it would be a GPI protein the cell may never make, and each
#: gene's canonical form is already in L1 on the N-after-SP rung.
GPI_REDRAW = {
    "l1-p2401": "redraw: BST1 isoform Q10588-2 starts at canonical residue 59, past the signal peptide; canonical in L1 (l1-p0187)",
    "l1-p2417": "redraw: CDH13 isoform P55290-4 has a divergent N-terminus with no signal peptide; canonical in L1 (l1-p0336)",
    "l1-p2476": "redraw: OPCML isoform Q14982-3 starts at canonical residue 42, past the signal peptide; canonical in L1 (l1-p1347)",
    "l1-p2727": "redraw: LYPD5 isoform Q6UWN5-2 starts at canonical residue 44, past the signal peptide; canonical in L1 (l1-p1180)",
    "l1-p2531": "redraw (ortholog panel): cyno ACHE gene model starts at canonical residue 116, past the signal peptide; human canonical in L1 (l1-p0017)",
}  # fmt: skip
WEAK_DELTA = 1

REVIEW_DECIDED = (
    'Becca, 2026-11-11 (applied the second review: "Apply all"; scope "tIII only for now")'
)
REVIEWED_BY = "Claude Code session type3-risk-review (read-only second review)"
#: Overrides from the independent second review of option B's 18 risk members, applied by
#: Becca's decision ("Apply all", 2026-11-11). Full evidence: TYPE3_NTERM.md, "Second review".
#: candidate_id -> (design, reason). A held design names the action that replaces building.
REVIEW = {
    "l1-p1346": ("tag_alone", "OPALIN: native N-out established (glycosylated N-terminal domain, PMID 18490449); weak delta does not flag cassette risk"),
    "l1-p1809": ("tag_alone", "SIGIRR: native type III already translocating a 118-residue Ig domain; weak delta does not flag cassette risk"),
    "l1-p2502": ("held", "redraw: SIGLEC15 isoform 2 is an N-terminal truncation from canonical residue 155 with no signal peptide; canonical in L1 (l1-p1814)"),
    "l1-p2484": ("held", "redraw: PDPN divergent first exon replaces the signal peptide; with a leader it nearly duplicates the canonical reference (l1-p1636)"),
    "l1-p2402": ("held", "redraw: BTN3A3 isoform 2 starts inside the IgV domain at Met43; canonical in L1 (l1-p0196)"),
    "l1-p2451": ("held", "redraw: IL9R isoform 3 has an internal signal peptide, a disrupted cytokine-binding module and no Box1; canonical in L1 (l1-p0948)"),
    "l1-p2529": ("held", "redraw: ZP3 isoform 3 starts at canonical Met52 with no signal peptide, and ZP3 is furin-processed and secreted (PMID 11804956); canonical in L1 (l1-p2359)"),
    "l1-p2424": ("held", "redraw: CPD isoform 2 starts at canonical residue 248, inside the first carboxypeptidase domain; canonical in L1 (l1-p0468)"),
    "l1-p0272": ("held", "swap reference: ordered CD200R1L starts at canonical Met22; use the signal-peptide form Q6Q8B3-1 (l2-p0276) on the N-after-SP rung"),
    "l1-p1136": ("held", "swap reference: ordered LRRC37B starts at canonical Met83; use the signal-peptide form Q96QE4-1 (l2-p0908) on the N-after-SP rung"),
    "l1-p0215": ("held", "swap reference: ordered CACHD1 starts at canonical Met135; use the signal-peptide form Q5VU97-1 (l2-p0248) on the N-after-SP rung"),
    "l1-p2549": ("held", "repair gene model or drop: cyno CEACAM8 is GPI-anchored (v1's 11-residue TM is the GPI signal) and its N-terminus is a gene-model error; human CEACAM8 is in L1 (l1-p0368)"),
}  # fmt: skip


LATE_DECIDED = (
    "Claude Code session l1-review-export-scope (subagent), under Becca's 2026-11-22 isoform redraw"
)
#: Members added 2026-11-21/22 whose subgroup needs a call, decided as the rung decided the same
#: gene's existing member: keyed by design_id, (design, reason with evidence).
LATE = {
    "l1-p3228": ("tag_alone",
                 "BMPR1B isoform (redraw 2026-11-22), first 127 residues identical to l1-p0181: "
                 "SignalP 6 calls no signal peptide (OTHER, P(SP) 0.000; run type3_sp6_2026_11_22, "
                 "processed://sprint4/s4e9-orf-design-fg-panel/type3-nterm/signalp6-runs/), so "
                 "UniProt's 1-13 annotation is not supported and TM1 is the signal-anchor, as for "
                 "l1-p0181: tag alone"),
}  # fmt: skip


def design_for(row: dict) -> tuple[str, str, str]:
    """(design, decided_by, reason) for one member of the rung."""
    # Decisions are keyed by design_id, which survives member renumbering (candidate_id is
    # plate order and was renumbered by PR #272 and again by PR #283).
    cid, sub = row["design_id"], row["subgroup"]
    if cid in REVIEW:
        design, reason = REVIEW[cid]
        return design, REVIEW_DECIDED, reason
    if cid in DELEGATED:
        design, reason, _alt = DELEGATED[cid]
        return design, DELEGATED_DECIDED, reason
    if cid in LATE:
        design, reason = LATE[cid]
        return design, LATE_DECIDED, reason
    if cid in GPI_REDRAW:
        return "held", GPI_REDRAW_DECIDED, GPI_REDRAW[cid]
    if cid in GPI_IGK:
        return ("igk_tag", GPI_DECIDED,
                "zero-TM GPI form without a signal peptide: never enters the ER as ordered; "
                "the leader supplies ER entry and the C-terminal GPI signal is intact")  # fmt: skip
    if sub == "no_tm":
        return "held", NO_DECISION, "zero-TM member with no signal peptide; undecided (LMLN)"
    if cid in CDS_REPAIR:
        return "held", NO_DECISION, f"CDS repair first: {CDS_REPAIR[cid]}"
    if sub == "sp_in_member":
        return ("held", NO_DECISION,
                "carries its gene's annotated signal peptide, uncalled by v1/v2; SignalP 6 decides "
                "between this rung and N-after-signal-peptide")  # fmt: skip
    if sub == "sp_missing":
        return ("igk_tag", B_DECIDED,
                "ordered form lacks the gene's signal peptide; the leader replaces it")  # fmt: skip
    if sub == "short_tail":
        return ("igk_tag", B_DECIDED,
                "native N-tail <= 20 residues; the cassette moves it past the ~20-residue point where "
                "signal-anchor insertion changes (Kocik 2012)")  # fmt: skip
    if row["gene_symbol"] == "KLB":
        return ("igk_tag", B_DECIDED,
                "internal signal-like segment; v1 puts the cassette in the cytosol with the tag alone")  # fmt: skip
    if sub == "mid_tail" and int(row["delta_charge_native"]) <= WEAK_DELTA:
        return ("igk_tag", B_DECIDED,
                f"flank-charge difference {row['delta_charge_native']} <= {WEAK_DELTA}: weak "
                "positive-inside drive for N-out")  # fmt: skip
    if sub in {"mid_tail", "long_tail", "dual_topology"}:
        return "tag_alone", B_DECIDED, f"genuine type III ({sub}) with no identified risk"
    raise ValueError(f"{cid}: unhandled subgroup {sub}")


def read_table(path: Path) -> list[dict]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def build() -> list[dict]:
    from deliverome_analysis.core.s3 import processed

    paths = {r["candidate_id"]: r for r in read_table(TAG_PATHS)}
    calls = read_table(CALLS)
    rung = sorted(c for c, r in paths.items() if r["placement"] == RUNG)
    assert sorted(r["candidate_id"] for r in calls) == rung, "calls table is not the current rung"
    manifest = {r["candidate_id"]: r for r in read_table(MANIFEST) if ld.is_active(r)}
    proteins = ld.designed_proteins(MANIFEST, library=ld.LIBRARY_L1)
    constructs = json.loads(Path(processed.download(CONSTRUCTS_KEY)).read_text())
    assert seqtools.translate(qc.P003_CASSETTE_NT) == CASSETTE_AA and len(CASSETTE_AA) == 51
    assert seqtools.translate(LEADER) == LEADER_AA and seqtools.translate(PLUS1) == PLUS1_AA

    on_rung = {c["design_id"] for c in calls}
    for name, keyed in (("REVIEW", REVIEW), ("GPI_REDRAW", GPI_REDRAW), ("CDS_REPAIR", CDS_REPAIR)):
        gone = sorted(set(keyed) - on_rung)
        if gone:
            print(f"{name}: {len(gone)} decided member(s) no longer on the rung: {gone}")
    rows = []
    for c in calls:
        cid = c["candidate_id"]
        design, decided_by, reason = design_for(c)
        p = paths[cid]
        native = proteins[cid]["sequence"]
        assert p["sequence_sha256"] == topology_sha(native), f"{cid}: sequence moved"
        row = {
            "candidate_id": cid,
            "design_id": manifest[cid].get("design_id", ""),
            "gene_symbol": c["gene_symbol"],
            "candidate_set": c["candidate_set"],
            "species": c["species"],
            "subgroup": c["subgroup"],
            "design": design,
            "native_aa": len(native),
            "protein_aa": "",
            "mature_aa": "",
            "coding_nt": "",
            "insert_nt": "",
            "cassette": "",
            "type_iis_sites": "",
            "decided_by": decided_by,
            "cassette_decided_by": "",
            "reviewed_by": REVIEWED_BY if c["design_id"] in REVIEW else "",
            "rejected_alternative": DELEGATED[c["design_id"]][2]
            if c["design_id"] in DELEGATED
            else "",
            "implemented_by": IMPLEMENTED_BY,
            "reason": reason,
            "protein_sha256": "",
            "insert": "",
        }
        if design != "held":
            orf = constructs[manifest[cid]["key"]]["construct"]
            assert seqtools.translate(orf) == native + "*", f"{cid}: ORF does not encode protein"
            assert orf[:3] == "ATG", f"{cid}: ORF does not start with ATG"
            head = (LEADER + PLUS1) if design == "igk_tag" else ""
            coding = head + orf[:3] + arm_cassette(cid) + orf[3:]
            mature = native[0] + CASSETTE_AA + native[1:]
            if design == "igk_tag":
                mature = PLUS1_AA + mature
            protein = (LEADER_AA if design == "igk_tag" else "") + mature
            got = seqtools.translate(coding)
            assert got == protein + "*" and got.count("*") == 1, f"{cid}: protein mismatch"
            insert = qc.wrap_insert(coding)
            assert qc.wrapper_sites_only(insert), f"{cid}: {seqtools.find_sites(insert, ENZYMES)}"
            row.update(
                protein_aa=len(protein),
                mature_aa=len(mature),
                coding_nt=len(coding),
                insert_nt=len(insert),
                cassette="full 51 aa; p003 tag codons, own linker encoding",
                type_iis_sites="Esp3I x2 (wrapper)",
                cassette_decided_by=CASSETTE_DECIDED,
                protein_sha256=hashlib.sha256(protein.encode()).hexdigest(),
                insert=insert,
                _protein=protein,
                _mature=mature,
                _native=native,
            )
        rows.append(row)
    feats = {r["candidate_id"]: r for r in read_table(FEATURES)}
    for row in rows:
        f = feats[row["candidate_id"]]
        row["_tm1"] = (int(f["tm1_start"]), int(f["tm1_end"])) if f.get("tm1_start") else None
    attach_refold(rows)
    vector_check([r for r in rows if r["insert"]])
    return rows


def topology_sha(aa: str) -> str:
    from deliverome_analysis.surfaceome import topology

    return topology.sequence_sha256(aa)


def vector_check(rows: list[dict]) -> None:
    """Every insert in the p0039 backbone: backbone sites plus the two wrapper Esp3I, no more."""
    import candidate_inputs
    from Bio import SeqIO

    bb = str(SeqIO.read(candidate_inputs.backbone(), "genbank").seq).upper()
    qc.check_inserts_in_backbone([(r["candidate_id"], r["insert"]) for r in rows], bb)


def render(rows: list[dict]) -> tuple[str, str]:
    cols = [c for c in rows[0] if c != "insert" and not c.startswith("_")]
    tsv = "\t".join(cols) + "\n"
    tsv += "".join("\t".join(str(r[c]) for c in cols) + "\n" for r in rows)
    fasta = "".join(
        f">{r['candidate_id']} {r['gene_symbol']} {r['species']} {r['design']}\n{r['insert']}\n"
        for r in rows
        if r["insert"]
    )
    return tsv, fasta


def render_proteins(rows: list[dict]) -> tuple[str, str]:
    """Protein table + FASTA for the built designs, coordinates 1-based in the precursor."""
    cols = [
        "candidate_id", "gene_symbol", "species", "design", "precursor_aa", "leader_span",
        "cassette_span", "native_tm1_span", "inserted_after_native_residue", "precursor_sha256",
        "native_sha256", "precursor", "mature", "native",
    ]  # fmt: skip
    out = []
    for r in rows:
        if r["design"] == "held":
            continue
        lead = len(LEADER_AA) + len(PLUS1_AA) if r["design"] == "igk_tag" else 0
        cas0 = lead + 1  # Met1, then the cassette
        tm = r["_tm1"]
        shift = lead + len(CASSETTE_AA)
        out.append(
            {
                "candidate_id": r["candidate_id"],
                "gene_symbol": r["gene_symbol"],
                "species": r["species"],
                "design": r["design"],
                "precursor_aa": len(r["_protein"]),
                "leader_span": f"1-{len(LEADER_AA)} (cleaved; +1 Asp)" if lead else "",
                "cassette_span": f"{cas0 + 1}-{cas0 + len(CASSETTE_AA)}",
                "native_tm1_span": f"{tm[0] + shift}-{tm[1] + shift}" if tm else "",
                "inserted_after_native_residue": 1,
                "precursor_sha256": hashlib.sha256(r["_protein"].encode()).hexdigest(),
                "native_sha256": topology_sha(r["_native"]),
                "precursor": r["_protein"],
                "mature": r["_mature"],
                "native": r["_native"],
            }
        )
    for o in out:
        p = o["precursor"]
        a, b = (int(x) for x in o["cassette_span"].split("-"))
        assert p[a - 1 : b] == CASSETTE_AA, o["candidate_id"]
    tsv = "\t".join(cols) + "\n" + "".join("\t".join(str(o[c]) for c in cols) + "\n" for o in out)
    fasta = "".join(
        f">{o['candidate_id']}|{o['gene_symbol']}|{o['species']}|{o['design']}\n{o['precursor']}\n"
        for o in out
    )
    return tsv, fasta


def _first(v: dict, *names: str) -> str:
    return next((v[n] for n in names if v.get(n, "") != ""), "")


def attach_refold(rows: list[dict]) -> None:
    """Copy each built member's refold verdict onto its row (blank until the file exists)."""
    import pandas as pd

    from deliverome_analysis.core.s3 import processed

    # Matched on the folded sequence's sha256 (tagged_sha256): the precursor for tag_alone,
    # the mature chain for igk_tag. Member ids drift with renumbering; sequences do not.
    by_sha: dict[str, dict] = {}
    local_verdicts = HERE / "tag_refold_verdicts.tsv"
    if local_verdicts.exists():
        df_loc = pd.read_csv(local_verdicts, sep="\t", dtype=str).fillna("")
        by_sha.update({r["tagged_sha256"]: r for r in df_loc.to_dict("records")})
    if processed.exists(REFOLD_KEY):
        df = pd.read_csv(processed.download(REFOLD_KEY), sep="\t", dtype=str).fillna("")
        by_sha.update({r["tagged_sha256"]: r for r in df.to_dict("records")})
    for row in rows:
        v: dict = {}
        row["refold_folded_form"] = row["refold_tagged_sha256"] = ""
        if row["design"] != "held":
            for form, seq in (("precursor", row["_protein"]), ("mature", row["_mature"])):
                sha = hashlib.sha256(seq.encode()).hexdigest()
                v = by_sha.get(sha, {})
                if v:
                    # Which sequence the verdict is for: tag_alone folds the precursor, igk_tag
                    # the mature (leader-cleaved) chain, so igk rows will not match protein_sha256.
                    row["refold_folded_form"], row["refold_tagged_sha256"] = form, sha
                    break
        verdict = _first(v, "verdict", "status")
        row["refold_verdict"] = (
            "pending" if (row["design"] != "held" and verdict in ("", "not_folded")) else verdict
        )
        row["refold_reason"] = _first(v, "verdict_reasons", "reason")
        row["refold_model_key"] = _first(v, "pdb_key", "model_key")
        row["refold_spytag003_longest_helix"] = _first(v, "spytag003_longest_helix")
        row["refold_tm1_longest_helix"] = _first(v, "tm1_longest_helix")
        row["refold_lddt_confident"] = _first(v, "lddt_confident")
        disp = REFOLD_DISPOSITIONS.get(row["design_id"])
        if disp and verdict == "fail":
            row["refold_disposition"] = (
                f"{disp[0]} ({DELEGATED_DECIDED}): {disp[1]}; rejected: {disp[2]}"
            )
        else:
            row["refold_disposition"] = "fail: undecided" if verdict == "fail" else ""


def _date(attribution: str) -> str:
    """The ISO date in an attribution string (every decision constant carries one)."""
    import re

    m = re.search(r"\d{4}-\d{2}-\d{2}", attribution)
    return m.group(0) if m else ""


def render_decisions(rows: list[dict]) -> str:
    for r in rows:
        r["decided_on"] = _date(r.get("decided_by", ""))
        assert r["decided_on"], f"{r['design_id']}: decision without a date"
    cols = ["design_id", "candidate_id", "gene_symbol", "species", "candidate_set", "subgroup",
            "design", "refold_verdict", "refold_folded_form", "refold_tagged_sha256", "refold_reason", "refold_disposition", "decided_on", "decided_by", "reviewed_by", "implemented_by", "reason", "rejected_alternative"]  # fmt: skip
    return (
        "\t".join(cols)
        + "\n"
        + "".join("\t".join(str(r.get(c, "")) for c in cols) + "\n" for r in rows)
    )


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--check", action="store_true", help="rebuild and compare with the record")
    args = ap.parse_args(argv)
    rows = build()
    tsv, fasta = render(rows)
    ptsv, pfasta = render_proteins(rows)
    decisions = render_decisions(rows)
    if args.check:
        same = (
            OUT_TSV.read_text() == tsv
            and OUT_FASTA.read_text() == fasta
            and OUT_PROTEIN_TSV.read_text() == ptsv
            and OUT_PROTEIN_FASTA.read_text() == pfasta
            and OUT_DECISIONS.read_text() == decisions
        )
        print("record matches rebuild" if same else "record DIFFERS from rebuild")
        return 0 if same else 1
    OUT_TSV.write_text(tsv)
    OUT_FASTA.write_text(fasta)
    OUT_PROTEIN_TSV.write_text(ptsv)
    OUT_PROTEIN_FASTA.write_text(pfasta)
    OUT_DECISIONS.write_text(decisions)
    counts: dict[str, int] = {}
    for r in rows:
        counts[r["design"]] = counts.get(r["design"], 0) + 1
    print(f"wrote {OUT_TSV.name} and {OUT_FASTA.name}: {counts}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
