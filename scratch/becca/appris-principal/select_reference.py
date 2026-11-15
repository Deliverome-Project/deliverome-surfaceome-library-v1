"""Choose one reference proteoform per library member: MANE + APPRIS where they agree,
adjudicated on observed evidence where they do not.

The policy, in order:

1. **MANE Select and the APPRIS principal agree** -> take it. This is the great majority,
   and it is the case where two independent methods -- RefSeq/Ensembl consensus on one side,
   structural and conservation inference on the other -- land on the same protein.
2. **Only one of the two exists** -> take that one, and say which.
3. **They disagree** -> do not pick on annotation alone. Work out the residues that actually
   distinguish the two candidates, and ask whether anything was ever *observed* there:
   peptides in PeptideAtlas, translated ORFs in the GENCODE Ribo-seq catalogue. Only the
   peptides decide: a gene with strong peptide evidence for one side and none for the other
   takes that side (``resolved_mane`` / ``resolved_appris``), including off MANE; weak,
   two-sided or absent evidence leaves it on MANE. Ribo-seq ORFs overlapping either side are
   counted and reported but never decide.

Agreement is tested on the **protein sequence**, not the transcript ID: a gene commonly has
several transcripts encoding one protein, and calling that a disagreement would be wrong.

Why evidence rather than a rule for the disagreements: 74 of the 111 are APPRIS
PRINCIPAL:4, APPRIS's own low-confidence tier, and only one is PRINCIPAL:1. The
disagreements are concentrated exactly where APPRIS is least sure, so a blanket
"APPRIS wins" would be backwards and a blanket "MANE wins" would discard the cases where
APPRIS is confident and right.

This also fixes the NPFFR2 class for free. There, MANE and APPRIS agree on the 420-residue
form and our UniProt canonical is the odd one out, so the policy picks a transcript-backed
sequence without anyone having to notice.

    uv run --locked python select_reference.py --surfaceome-repo ~/Git/accessible-surfaceome
"""

from __future__ import annotations

import argparse
import csv
import difflib
import gzip
import hashlib
import json
import re
import time
from collections import Counter, defaultdict
from pathlib import Path

import requests

from gpi import gpi_genes as gpi_gene_set
from select_principal import (
    APPRIS_RELEASE,
    UNIPROT_HUMAN,
    UNIPROT_RELEASE,
    load_uniprot_proteoforms,
    parse_appris,
    parse_ensembl_pep,
    parse_uniprot_fasta,
    select_principal,
)

ATLAS_BUILD_ID = 607
ATLAS_BUILD = "Human 2026-01"
ATLAS_CGI = "https://db.systemsbiology.net/sbeams/cgi/PeptideAtlas/GetPeptides"
ATLAS_BATCH = 20  # accessions per request; PeptideAtlas asks clients not to hammer the CGI
ATLAS_PAUSE = 3.0  # seconds between requests

# A peptide shorter than this is not diagnostic; longer than this is not routinely observed.
MIN_PEP, MAX_PEP = 7, 30

# A difference shorter than this is not worth calling a discriminating region -- a one or two
# residue change gives no peptide that is unique to it in practice.
MIN_REGION = 3

# Genes where both candidates are genuinely expressed and the evidence sources disagree, so
# neither is "the" reference. The peptide-favoured form stays the reference and the other is
# carried as an explicit second lane rather than discarded. CD44: 49 diagnostic peptides and
# MS 229 vs 158 favour the long variant-exon form, while OpenProt ribosome profiling favours
# the short standard form 17 to 5 -- CD44s and CD44v are both real proteins. ATP2B2:
# peptides call MANE and OpenProt is tied.
BOTH_LANES = {"CD44", "ATP2B2"}

# Default source precedence for the reference, strongest first. MANE is a RefSeq/Ensembl
# agreement and transcript-backed by construction; APPRIS is structural and conservation
# inference; the UniProt canonical is the last resort, since its selection rule permits a
# pure length tiebreak with no transcript support. The one exception is a MANE/APPRIS
# disagreement that PeptideAtlas peptides resolve strongly for APPRIS: there the evidence
# rung goes first (precedence_options). OpenProt enters only through BOTH_LANES, and
# GENCODE Ribo-seq ORFs only annotate.
PRECEDENCE = ("mane", "appris", "evidence", "uniprot_canonical")

# A reference must not be a soluble form of a surface protein. PILRA and LMLN were both
# moved onto one by taking APPRIS unconditionally, so each rung is checked for membrane
# topology before it is accepted and skipped if it is explicitly non-surface. "Explicitly"
# matters: a form with no topology call is not rejected, only a form called SP or GLOB that
# carries no GPI anchor.
SURFACE_LABELS = {"TM", "SP+TM", "BETA"}
# The two explicitly soluble calls. A form carrying one of these and no GPI anchor is
# not reachable from outside the cell, so it cannot serve as a library reference.
NON_SURFACE_LABELS = {"SP", "GLOB"}

# Peptides needed on one side, with none on the other, before the evidence decides it. One
# short proteotypic peptide is a real observation but a thin basis for changing a reference,
# and the two sides are not symmetric in practice: MANE-supporting hits run to dozens of
# long peptides while APPRIS-supporting hits are often a single 7-mer. Below this, the gene
# is flagged for a person rather than resolved.
STRONG_EVIDENCE = 3

# ...and its best peptide has to have been seen in at least this many distinct samples.
# Peptide count alone let two shallow cases through: PTPRG passed on three peptides whose
# best was seen in one sample. The data separates cleanly -- the cases that pass both tests
# run 17 to 137 samples, every case that fails runs 1 to 8 -- so the cut sits in the gap.
MIN_SAMPLES = 10

RIBOSEQ_SETS = ("Ribo-seq_ORFs.all.bed", "Ribo-seq_ORFs.primary.bed",
                "Ribo-seq_ORFs.comprehensive.bed")


def precedence_options(
    decision: str, mane: str | None, appris: str | None, current: str | None
) -> list[tuple[str, str | None]]:
    """The rungs to try for one gene, strongest first, given its final decision.

    The default order is :data:`PRECEDENCE`: MANE, APPRIS, the evidence pick, the UniProt
    canonical. Where MANE and APPRIS disagree, the policy is that observed evidence decides,
    so a gene resolved for APPRIS (``resolved_appris``: at least :data:`STRONG_EVIDENCE`
    proteotypic PeptideAtlas peptides on APPRIS's discriminating residues, the best seen in
    at least :data:`MIN_SAMPLES` samples, and none on MANE's) puts the evidence rung first.
    A gene resolved for MANE needs no reordering: MANE already leads. Weak, two-sided or
    absent evidence leaves MANE first, which is the stated default for disagreements.

    Ribo-seq plays no part here. GENCODE Ribo-seq ORFs overlapping either side's
    discriminating residues are counted (``riboseq_supports_mane`` / ``_appris``) and
    reported, but they never set a decision: only PeptideAtlas peptides do. Making
    translation evidence decide is a policy change for Becca, not a code fix.
    """
    evidence = {
        "resolved_appris": appris,
        "resolved_mane": mane,
        "weak_evidence_mane": mane,
    }.get(decision)
    options = [
        ("mane", mane),
        ("appris", appris),
        ("evidence", evidence),
        ("uniprot_canonical", current),
    ]
    if decision == "resolved_appris":
        options.insert(0, options.pop(2))
    return options


def choose_rung(
    options: list[tuple[str, str | None]], surface_ok
) -> tuple[str, str | None, list[str], bool]:
    """(rung, sequence, rungs skipped as non-surface, every candidate non-surface).

    Takes the first available rung whose form is not explicitly soluble. If every available
    form is, keeps the strongest one and says so rather than leaving the gene empty.
    """
    skipped = []
    for rung, seq in options:
        if seq is None:
            continue
        if surface_ok(seq) is False:
            skipped.append(rung)
            continue
        return rung, seq, skipped, False
    for rung, seq in options:
        if seq is not None:
            return rung, seq, skipped, True
    return "", None, skipped, False


def translations_by_gene(path: Path) -> dict[str, set[str]]:
    """{ENSG: {protein sequence}} over *every* translation in the Ensembl peptide FASTA.

    The transcripts APPRIS lists are not the whole gene -- it annotates a subset -- so the
    transcript-support test has to run against the full file or it reports sequences as
    unsupported when a transcript does encode them. SUSD6 is the case that caught this.
    """
    out: dict[str, set[str]] = defaultdict(set)
    cur, buf, hdr = None, [], ""
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if cur:
                    m = re.search(r"\bgene:(\S+)", hdr)
                    if m:
                        out[m.group(1).split(".")[0]].add("".join(buf).rstrip("*"))
                hdr = line
                cur = True
                buf = []
            elif cur:
                buf.append(line.strip())
    if cur:
        m = re.search(r"\bgene:(\S+)", hdr)
        if m:
            out[m.group(1).split(".")[0]].add("".join(buf).rstrip("*"))
    return dict(out)


# --------------------------------------------------------------------- evidence sources


def fetch_peptides(
    accessions: list[str], cache: Path, offline: bool = False
) -> dict[str, dict[str, tuple[int, int]]]:
    """{accession: {peptide: (n_observations, n_samples)}} from PeptideAtlas, cached on disk.

    ``n_samples`` is how many distinct samples the peptide was seen in, and it is the
    quantity that separates a real observation from a one-off: the same build reports
    peptides seen in a single sample and peptides seen in 1,852, and counting those equally
    is what made the first pass over-trust single-peptide hits.

    Queried in batches through ``biosequence_name_constraint``, which accepts a
    semicolon-separated list. The response does not say which protein each peptide came
    from, which does not matter here: every test below is containment in one candidate
    sequence and absence from the other, so a peptide from an unrelated protein in the same
    batch cannot create a false hit unless it is also a substring of this protein, and that
    case is excluded explicitly.
    """
    # Stored as one peptide table per batch plus an accession->batch index. Keeping a full
    # copy of the pooled result under every accession in the batch duplicated it ~20x and
    # ran to hundreds of megabytes.
    store = cache / f"peptideatlas_build{ATLAS_BUILD_ID}_pools.json"
    legacy = cache / f"peptideatlas_build{ATLAS_BUILD_ID}_meta.json"
    if store.exists():
        blob = json.loads(store.read_text())
    elif legacy.exists():
        pools, index = {}, {}
        for acc, peps in json.loads(legacy.read_text()).items():
            key = str(hash(tuple(sorted(peps))) & 0xFFFFFFFF)
            pools.setdefault(key, peps)
            index[acc] = key
        blob = {"pools": pools, "index": index}
        store.write_text(json.dumps(blob))
        legacy.unlink()
        print(f"  migrated peptide cache to {store.name}")
    else:
        blob = {"pools": {}, "index": {}}
    pools: dict[str, dict] = blob["pools"]
    index: dict[str, str] = blob["index"]
    seen = index
    todo = [a for a in accessions if a not in seen]
    if offline and todo:
        # Offline means "read the pinned cache", never "pretend nothing was observed": the
        # old offline branch returned no peptides at all, which is how a published
        # reference_selection.tsv came to carry zero evidence for all 111 disagreements.
        raise SystemExit(
            f"--offline: {len(todo)} accession(s) are not in {store.name}, e.g. {todo[:3]}"
        )
    for i in range(0, len(todo), ATLAS_BATCH):
        batch = todo[i : i + ATLAS_BATCH]
        resp = requests.get(
            ATLAS_CGI,
            params={
                "atlas_build_id": ATLAS_BUILD_ID,
                "biosequence_name_constraint": ";".join(batch),
                "apply_action": "QUERY",
                "output_mode": "tsv",
            },
            timeout=600,
        )
        resp.raise_for_status()
        rows = list(csv.DictReader(resp.text.splitlines(), delimiter="\t"))
        peptides = {}
        for r in rows:
            seq = r.get("peptide_sequence")
            if not seq or r.get("organism_full_name") not in (None, "", "Homo sapiens"):
                continue
            peptides[seq] = (
                int(r.get("n_observations") or 0),
                int(r.get("n_samples") or 0),
            )
        # The batch response is pooled; each accession in the batch gets the pool, and the
        # containment tests separate them.
        key = str(i // ATLAS_BATCH) + "_" + str(len(pools))
        pools[key] = peptides
        for acc in batch:
            index[acc] = key
        print(f"  PeptideAtlas {i + len(batch)}/{len(todo)}: {len(peptides)} peptides")
        time.sleep(ATLAS_PAUSE)
    store.write_text(json.dumps({"pools": pools, "index": index}))
    return {
        a: {k: tuple(v) for k, v in pools.get(index.get(a, ""), {}).items()}
        for a in accessions
    }


def load_riboseq(cache: Path) -> tuple[dict[str, list[str]], dict[str, set[str]]]:
    """Per gene: the ORF types, and the ORF *translations*.

    The types alone cannot adjudicate two candidate CDSs -- none of them describes an
    in-frame alternative CDS -- so the translations are loaded too and tested for overlap
    with the discriminating regions. That is the version of this test that can actually
    produce a hit.
    """
    types: dict[str, list[str]] = defaultdict(list)
    for fname in RIBOSEQ_SETS:
        path = cache / fname
        if not path.exists():
            continue
        for line in path.read_text().splitlines():
            if not line:
                continue
            parts = line.split("\t")[3].split(":")
            if len(parts) >= 3:
                types[parts[1]].append(parts[-1])

    seqs: dict[str, set[str]] = defaultdict(set)
    for fname in RIBOSEQ_SETS:
        path = cache / fname.replace(".bed", ".faa")
        if not path.exists():
            continue
        header, buf = None, []
        for line in path.read_text().splitlines():
            if line.startswith(">"):
                if header:
                    seqs[header].add("".join(buf))
                parts = line[1:].split()[0].split(":")
                header = parts[1] if len(parts) >= 3 else None
                buf = []
            elif header:
                buf.append(line.strip())
        if header:
            seqs[header].add("".join(buf))
    return dict(types), dict(seqs)


def riboseq_overlap(regions: list[tuple[int, str]], orf_seqs: set[str]) -> list[str]:
    """Subsequences of a discriminating region that appear in a Ribo-seq ORF translation."""
    hits = set()
    for _start, region in regions:
        for i in range(len(region)):
            for ln in range(MIN_PEP, min(MAX_PEP, len(region) - i) + 1):
                sub = region[i : i + ln]
                if any(sub in orf for orf in orf_seqs):
                    hits.add(sub)
    return sorted(hits, key=len, reverse=True)


# --------------------------------------------------------------------- adjudication


def discriminating_regions(a: str, b: str) -> tuple[list[tuple[int, str]], list[tuple[int, str]]]:
    """Stretches present in `a` but not `b`, and vice versa, as (1-based start, subsequence)."""
    sm = difflib.SequenceMatcher(None, a, b, autojunk=False)
    only_a, only_b = [], []
    for tag, i1, i2, j1, j2 in sm.get_opcodes():
        if tag == "equal":
            continue
        if i2 - i1 >= MIN_REGION:
            only_a.append((i1 + 1, a[i1:i2]))
        if j2 - j1 >= MIN_REGION:
            only_b.append((j1 + 1, b[j1:j2]))
    return only_a, only_b


def peptides_supporting(
    regions: list[tuple[int, str]],
    observed: dict[str, tuple[int, int]],
    other: str,
    proteotypic: callable,
) -> list[tuple[str, int, int]]:
    """Observed peptides lying inside a region unique to one candidate and absent from the other.

    Every substring of the region in the detectable length window is tested, so this finds
    tryptic, semi-tryptic and missed-cleavage peptides alike without assuming a protease.

    A hit also has to be *proteotypic* -- found in no other gene's proteoforms. Without that,
    a short sequence shared with an unrelated protein would be read as evidence for this
    gene's splice form, which is the one way this test could manufacture a false answer.
    """
    hits: dict[str, tuple[int, int]] = {}
    for _start, region in regions:
        for i in range(len(region)):
            for ln in range(MIN_PEP, min(MAX_PEP, len(region) - i) + 1):
                sub = region[i : i + ln]
                if sub in observed and sub not in other and proteotypic(sub):
                    hits[sub] = observed[sub]
    # Deepest observation first -- n_samples, then n_observations, then length.
    return sorted(
        ((k, v[0], v[1]) for k, v in hits.items()),
        key=lambda x: (x[2], x[1], len(x[0])),
        reverse=True,
    )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--surfaceome-repo", default="~/Git/accessible-surfaceome")
    ap.add_argument("--cache-dir", default="~/.cache/deliverome/appris")
    ap.add_argument("--riboseq-dir", default="~/.cache/deliverome/riboseq")
    ap.add_argument("--library", default="../tag-site-architecture/library_members.json")
    ap.add_argument("--out-dir", default=".")
    ap.add_argument(
        "--overrides", default="reference_overrides.tsv",
        help="genes whose reference is named outright. The precedence ranks forms by how "
             "well sourced they are, which cannot see whether the best-sourced form is "
             "anchored at all -- not that a GPI gene's carries no omega site, nor that "
             "MICA's keeps no transmembrane helix.",
    )
    ap.add_argument("--offline", action="store_true",
                    help="read PeptideAtlas only from the cached build; fail if any is missing")
    args = ap.parse_args()

    repo = Path(args.surfaceome_repo).expanduser()
    cache = Path(args.cache_dir).expanduser()
    out_dir = Path(args.out_dir)
    seqdir = repo / "data/external/deeptmhmm_surfaceome_inputs"

    library = sorted(set(json.loads(Path(args.library).read_text())))
    cohort = parse_uniprot_fasta(seqdir / "human_canonical_non_hla.fasta")
    cohort_by_symbol = {r["symbol"]: r for r in cohort.values() if r["symbol"]}
    uniprot = load_uniprot_proteoforms(cache / UNIPROT_HUMAN)
    appris = parse_appris(cache / f"appris_homo_sapiens_{APPRIS_RELEASE}.txt")
    pep = parse_ensembl_pep(cache / "Homo_sapiens.GRCh38.pep.all.fa.gz")
    gene_translations_all = translations_by_gene(
        cache / "Homo_sapiens.GRCh38.pep.all.fa.gz"
    )
    riboseq, riboseq_seqs = load_riboseq(Path(args.riboseq_dir).expanduser())
    print(
        f"library {len(library)}; APPRIS genes {len(appris)}; Ensembl peptides {len(pep)}; "
        f"Ribo-seq genes {len(riboseq)} ({len(riboseq_seqs)} with translations)"
    )

    appris_by_symbol: dict[str, list[dict]] = defaultdict(list)
    for rows in appris.values():
        for r in rows:
            appris_by_symbol[r["symbol"]].append(r)

    ensg_for_acc: dict[str, set[str]] = defaultdict(set)
    with (repo / "data/external/uniprot_ensembl_xrefs/ensg_to_uniprot.tsv").open() as fh:
        next(fh)
        for line in fh:
            g, a = line.split("\t")[:2]
            ensg_for_acc[a].add(g)

    # ---- pass 1: MANE vs APPRIS per gene --------------------------------------
    records = []
    for symbol in library:
        stored = cohort_by_symbol.get(symbol)
        pf = uniprot.get(symbol)
        current = stored or (pf[0] if pf and not pf[0]["is_isoform"] else None)
        acc = current["accession"] if current else ""

        rows = [r for g in sorted(ensg_for_acc.get(acc, ())) for r in appris.get(g, [])]
        if not rows:
            rows = appris_by_symbol.get(symbol, [])

        # Is the sequence we currently build from encoded by any transcript at all? NPFFR2
        # showed this is not a given: UniProt may display a sequence that no Ensembl, RefSeq
        # or CCDS entry carries, and nothing downstream would notice. Checked against every
        # translation of the gene, not only the ones APPRIS lists.
        gene_ensgs = sorted(ensg_for_acc.get(acc, ())) or sorted(
            {r["ensg"] for r in rows}
        )
        gene_translations = {
            seq for g in gene_ensgs for seq in gene_translations_all.get(g, ())
        }

        mane_rows = [r for r in rows if r["mane"] and r["enst"] in pep]
        mane_seq = pep[mane_rows[0]["enst"]] if mane_rows else None
        sel = select_principal(rows, pep) if rows else None
        appris_seq = sel["seq"] if sel else None

        rec = {
            "symbol": symbol,
            "current_accession": acc,
            "current_length": len(current["seq"]) if current else "",
            "mane_enst": mane_rows[0]["enst"] if mane_rows else "",
            # What APPRIS itself calls the MANE transcript. ALTERNATIVE:1 means the model is
            # conserved in at least three tested non-primate species, ALTERNATIVE:2 fewer --
            # so this says how hard APPRIS is arguing against MANE, not merely that it is.
            "mane_appris_label": mane_rows[0]["label"] if mane_rows else "",
            "mane_length": len(mane_seq) if mane_seq else "",
            "appris_class": sel["label"] if sel else "",
            "appris_enst": sel["enst"] if sel else "",
            "appris_length": len(appris_seq) if appris_seq else "",
            "_mane_seq": mane_seq,
            "_appris_seq": appris_seq,
            "_current_seq": current["seq"] if current else None,
            # "" means the question could not be asked -- no translations found for the
            # gene -- which is not the same as the sequence being unsupported.
            "current_has_ensembl_transcript": (
                str(current["seq"] in gene_translations)
                if current and gene_translations
                else ""
            ),
            "n_ensembl_translations": len(gene_translations),
        }
        # The decision records how the two annotations relate; the rung chosen further down
        # records which source the reference actually came from. They are different
        # questions and both are kept.
        if mane_seq and appris_seq:
            rec["decision"] = "mane_appris_agree" if mane_seq == appris_seq else "disagree"
        elif mane_seq:
            rec["decision"] = "mane_only"
        elif appris_seq:
            rec["decision"] = "appris_only"
        else:
            # No MANE and no APPRIS principal: the KIR cluster and a few olfactory/taste
            # receptors, which sit on alternate haplotypes rather than the primary assembly.
            rec["decision"] = "no_reference_uniprot_canonical"
        rec["_candidates"] = {"mane": mane_seq, "appris": appris_seq}
        records.append(rec)

    # ---- topology, for the soluble guard ---------------------------------------
    # The processed cohort TSVs carry no sequence column, only accessions, so the label has
    # to be joined back to a sequence through UniProt. Keying on a "sequence" field that does
    # not exist silently produced an empty map, and PILRA kept a soluble reference because
    # the guard had no call to reject it.
    seq_for_acc: dict[str, str] = {}
    for forms in uniprot.values():
        for pf in forms:
            seq_for_acc[pf["accession"]] = pf["seq"]
    for rec_ in cohort.values():
        seq_for_acc.setdefault(rec_["accession"], rec_["seq"])

    label_for_seq: dict[str, str] = {}
    unmapped = 0
    for name in ("deeptmhmm_human_canonical.tsv", "deeptmhmm_human_isoforms.tsv"):
        path = repo / "data/processed/deeptmhmm" / name
        if not path.exists():
            continue
        with path.open() as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                acc = row.get("uniprot_accession_full") or row.get("uniprot_accession")
                seq = seq_for_acc.get(acc) or seq_for_acc.get((acc or "").split("-")[0])
                if seq:
                    label_for_seq.setdefault(seq, row["deeptmhmm_label"])
                else:
                    unmapped += 1
    print(f"topology labels keyed to a sequence: {len(label_for_seq)} ({unmapped} unmapped)")
    for extra in (Path("library_reference_missing_topology.3line"), Path("topology_shortlist.3line")):
        if not extra.exists():
            continue
        lines = extra.read_text().splitlines()
        for i in range(len(lines) - 2):
            if lines[i].startswith(">"):
                _n, lab = lines[i][1:].rsplit("|", 1)
                label_for_seq.setdefault(lines[i + 1].strip(), lab.strip())
    gpi_genes = gpi_gene_set()

    def surface_ok(seq: str | None, gene: str) -> bool | None:
        """True / False / None when no topology call exists. Only an explicit SP or GLOB,
        on a gene the architecture attribute does not call GPI-anchored, is a rejection."""
        if not seq:
            return None
        lab = label_for_seq.get(seq)
        if lab is None:
            return None
        return lab in SURFACE_LABELS or gene in gpi_genes

    tally = Counter(r["decision"] for r in records)
    print("\npass 1 (annotation only):", dict(tally.most_common()))

    # ---- pass 2: evidence for the disagreements --------------------------------
    disagree = [r for r in records if r["decision"] == "disagree"]
    accessions = sorted({r["current_accession"] for r in disagree if r["current_accession"]})
    print(f"\nadjudicating {len(disagree)} disagreements over {len(accessions)} accessions")
    observed = (
        fetch_peptides(accessions, cache, offline=args.offline)
    )

    # Proteotypic test: the peptide must occur in no gene other than this one. Checked
    # against every reviewed human proteoform, canonical and isoform.
    all_proteoforms = [(sym, p["seq"]) for sym, recs in uniprot.items() for p in recs]

    def make_proteotypic(symbol: str):
        def is_proteotypic(sub: str) -> bool:
            return not any(sub in seq for sym, seq in all_proteoforms if sym != symbol)

        return is_proteotypic

    for r in disagree:
        m, a = r["_mane_seq"], r["_appris_seq"]
        only_m, only_a = discriminating_regions(m, a)
        obs = observed.get(r["current_accession"], {})
        proteotypic = make_proteotypic(r["symbol"])
        sup_m = peptides_supporting(only_m, obs, a, proteotypic)
        sup_a = peptides_supporting(only_a, obs, m, proteotypic)

        orfs = riboseq.get(r["symbol"], [])
        orf_seqs = riboseq_seqs.get(r["symbol"], set())
        ribo_m = riboseq_overlap(only_m, orf_seqs)
        ribo_a = riboseq_overlap(only_a, orf_seqs)
        r |= {
            "mane_unique_residues": sum(len(s) for _, s in only_m),
            "appris_unique_residues": sum(len(s) for _, s in only_a),
            "peptides_supporting_mane": len(sup_m),
            "peptides_supporting_appris": len(sup_a),
            "peptide_examples_mane": ";".join(x[0] for x in sup_m[:3]),
            "peptide_examples_appris": ";".join(x[0] for x in sup_a[:3]),
            "mane_best_n_samples": max((x[2] for x in sup_m), default=0),
            "appris_best_n_samples": max((x[2] for x in sup_a), default=0),
            "mane_best_n_observations": max((x[1] for x in sup_m), default=0),
            "appris_best_n_observations": max((x[1] for x in sup_a), default=0),
            "riboseq_orfs": ";".join(sorted(set(orfs))),
            "riboseq_n_orfs": len(orf_seqs),
            "riboseq_supports_mane": len(ribo_m),
            "riboseq_supports_appris": len(ribo_a),
            "n_peptides_for_protein": len(obs),
        }
        def strong(hits: list[tuple[str, int, int]]) -> bool:
            return len(hits) >= STRONG_EVIDENCE and max(h[2] for h in hits) >= MIN_SAMPLES

        if sup_m and not sup_a:
            r["decision"] = "resolved_mane" if strong(sup_m) else "weak_evidence_mane"
        elif sup_a and not sup_m:
            r["decision"] = "resolved_appris" if strong(sup_a) else "weak_evidence_appris"
        elif sup_m and sup_a:
            r["decision"] = "both_observed"
        else:
            # No observation either way. MANE is the default: it is a RefSeq/Ensembl
            # agreement and so transcript-backed by construction, whereas most of these
            # disagreements are APPRIS PRINCIPAL:4, its own low-confidence tier.
            r["decision"] = "unresolved_default_mane"

    # ---- apply the precedence ---------------------------------------------------
    # After pass 2, so that an evidence decision can move a gene. Run before it (as it was
    # until 2026-11-14), every gene was picked on its pass-1 decision, ``resolved_appris``
    # did not exist yet, and observed evidence could never take a gene off MANE.
    for r in records:
        cands = r.pop("_candidates")
        options = precedence_options(r["decision"], cands["mane"], cands["appris"],
                                     r["_current_seq"])
        rung, seq, skipped, all_non_surface = choose_rung(
            options, lambda s, g=r["symbol"]: surface_ok(s, g)
        )
        if all_non_surface:
            r["all_candidates_non_surface"] = True
        r["reference_rung"] = rung
        r["soluble_guard_skipped"] = ";".join(skipped)
        r["_picked_seq"] = seq

    # ---- named references ------------------------------------------------------
    override: dict[str, str] = {}
    op = Path(args.overrides)
    if op.exists():
        with op.open() as fh:
            override = {row["gene"]: row["form"] for row in csv.DictReader(fh, delimiter="\t")}
    else:
        # Silence here would be the dangerous outcome: the run would simply pick the
        # best-sourced form for every gene, which for the named genes is the one
        # without the anchor, and nothing would say so.
        print(
            f"  WARNING: no named-reference file at {op}. Seven genes are expected to name "
            "their reference outright -- six GPI genes plus MICA, whose precedence-ranked "
            "transcript keeps no transmembrane helix; without it they fall back to the form "
            "the precedence ranks first, which carries no anchor."
        )
    applied = 0
    for r in records:
        acc = override.get(r["symbol"])
        if not acc:
            continue
        seq = next((p["seq"] for p in uniprot.get(r["symbol"], []) if p["accession"] == acc), None)
        if seq is None:
            print(f"  override for {r['symbol']} names {acc}, which no stored form matches")
            continue
        r["_picked_seq"] = seq
        r["reference_rung"] = "named"
        applied += 1
    if applied:
        print(f"named references applied: {applied}")

    # ---- selected sequence and what it changes ---------------------------------
    for r in records:
        pick = r.get("_picked_seq")
        # The rung IS the source now; the precedence decides it.
        r["selected_source"] = r.get("reference_rung") or ""
        # Anything a person still has to look at, called out rather than buried in a code.
        if r["symbol"] in BOTH_LANES:
            r["second_lane"] = "appris"
            r["second_lane_length"] = (
                len(r["_appris_seq"]) if r.get("_appris_seq") else ""
            )
        r["needs_review"] = r["decision"] in (
            "weak_evidence_mane", "weak_evidence_appris", "both_observed"
        )
        r["selected_length"] = len(pick) if pick else ""
        lab = label_for_seq.get(pick, "") if pick else ""
        r["selected_topology_label"] = lab
        r["selected_gene_is_gpi"] = r["symbol"] in gpi_genes
        # The restriction the guard alone cannot make: the guard only declines to MOVE a
        # reference onto a soluble form, so a gene whose every candidate is soluble keeps
        # one. This flags those for removal from the library instead.
        r["restrict_soluble_not_gpi"] = bool(
            lab in NON_SURFACE_LABELS and r["symbol"] not in gpi_genes
        )
        pf = uniprot.get(r["symbol"], [])
        r["selected_uniprot_equivalent"] = (
            next((p["accession"] for p in pf if p["seq"] == pick), "") if pick else ""
        )
        r["changes_current_reference"] = (
            "" if pick is None or r["_current_seq"] is None
            else str(pick != r["_current_seq"])
        )
        if pick is not None and r["_current_seq"] is not None and pick != r["_current_seq"]:
            r["current_minus_selected"] = len(r["_current_seq"]) - len(pick)

    columns = [
        "symbol", "decision", "reference_rung", "soluble_guard_skipped",
        "all_candidates_non_surface", "needs_review", "selected_source",
        "selected_topology_label", "selected_gene_is_gpi", "restrict_soluble_not_gpi",
        "second_lane", "second_lane_length", "selected_length",
        "selected_uniprot_equivalent", "changes_current_reference",
        "current_accession", "current_length", "current_minus_selected",
        "current_has_ensembl_transcript", "n_ensembl_translations",
        "mane_enst", "mane_length", "mane_appris_label",
        "appris_class", "appris_enst", "appris_length",
        "mane_unique_residues", "appris_unique_residues",
        "peptides_supporting_mane", "peptides_supporting_appris",
        "mane_best_n_samples", "appris_best_n_samples",
        "mane_best_n_observations", "appris_best_n_observations",
        "peptide_examples_mane", "peptide_examples_appris",
        "n_peptides_for_protein",
        "riboseq_orfs", "riboseq_n_orfs", "riboseq_supports_mane", "riboseq_supports_appris",
    ]
    with (out_dir / "reference_selection.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(records)
    print(f"\nwrote reference_selection.tsv ({len(records)} genes)")

    with (out_dir / "reference_disagreements.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(sorted(disagree, key=lambda r: r["symbol"]))
    print(f"wrote reference_disagreements.tsv ({len(disagree)} genes)")

    # The chosen protein itself, by hash. reference_selection.tsv carries the source and the
    # length but not the sequence, so without this a change of protein at the same source
    # and length would pass a byte-for-byte check of that table unnoticed.
    with (out_dir / "reference_selection_sequences.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["symbol", "selected_source", "selected_length", "selected_sequence_sha256"])
        for r in records:
            pick = r.get("_picked_seq")
            w.writerow(
                [
                    r["symbol"],
                    r["selected_source"],
                    r["selected_length"],
                    hashlib.sha256(pick.encode()).hexdigest() if pick else "",
                ]
            )
    print("wrote reference_selection_sequences.tsv")

    changed = [r for r in records if r.get("changes_current_reference") == "True"]
    summary = {
        "appris_release": APPRIS_RELEASE,
        "uniprot_release": UNIPROT_RELEASE,
        "peptideatlas_build": ATLAS_BUILD,
        "library_members": len(library),
        "decision": dict(Counter(r["decision"] for r in records).most_common()),
        "selected_source": dict(
            Counter(r["selected_source"] for r in records if r["selected_source"]).most_common()
        ),
        "changes_current_reference": len(changed),
        "needs_review": sum(1 for r in records if r.get("needs_review")),
        "both_lanes": sorted(BOTH_LANES),
        "current_reference_unsupported_by_any_transcript": sum(
            1 for r in records if r.get("current_has_ensembl_transcript") == "False"
        ),
        "transcript_support_unknown": sum(
            1 for r in records if r.get("current_has_ensembl_transcript") == ""
        ),
        "disagreements": {
            "n": len(disagree),
            "appris_class": dict(Counter(r["appris_class"] for r in disagree).most_common()),
            "with_any_peptide_data": sum(
                1 for r in disagree if r.get("n_peptides_for_protein")
            ),
            "with_riboseq_orf": sum(1 for r in disagree if r.get("riboseq_orfs")),
            "riboseq_overlapping_a_discriminating_region": sum(
                1
                for r in disagree
                if r.get("riboseq_supports_mane") or r.get("riboseq_supports_appris")
            ),
            "mane_appris_label": dict(
                Counter(r.get("mane_appris_label", "") for r in disagree).most_common()
            ),
            "strong_evidence_threshold": STRONG_EVIDENCE,
            "min_samples_threshold": MIN_SAMPLES,
        },
    }
    (out_dir / "reference_selection_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
