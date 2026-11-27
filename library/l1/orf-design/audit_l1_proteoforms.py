"""Audit the designed proteoform of every active L1 member, one row per member, with flag classes.

Requested by Becca on 2026-11-11 ("systematically sort through the library proteoform design ...
review each individual sequence and record overrides or modifications to mane or appris").
Implemented by the Claude Code session l1-proteoform-audit (branch
``claude/l1-proteoform-audit-c4b011``). Report only: this script changes no manifest, panel,
CDS record or S3 key. Decisions taken on its flags live in ``l1_proteoform_decisions.tsv``.

What each row says about one member (the member's protein is its designed CDS translated,
``library_design``; every comparison below is by sequence, never by accession or length):

* identity: ``design_id`` (the design key), ``candidate_id`` (the manifest row id; differs from
  ``design_id`` for the forms moved in from L2), ``member_id`` (``l1_member_renumbering.tsv``,
  blank for members drawn after that table was written), set, species, gene, transcript;
* build source: the reference ladder's decision for the gene (``reference_selection.tsv``) and
  which of the gene's MANE Select, APPRIS principal and UniProt canonical proteins the member
  equals; whether MANE and APPRIS agree on the protein sequence;
* CDS: start codon, native or replaced (``start_codon_replaced`` in the CDS record), internal
  stops, frame, and whether the translation equals the stated protein (Ensembl 115 peptide
  for human transcripts, Ensembl 116 for mouse / cyno, the rescue record's ``protein`` for
  expression-rescue arms);
* length and alignment against the comparator: the gene's UniProt canonical for library
  references and rescue arms; the gene's L1 library reference for isoforms and orthologs
  (falling back to the canonical when the gene has none). Global alignment, BLOSUM62, gap
  11/1, free end gaps (as ``review_ortholog_correspondence.py``);
* signal peptide: UniProt (``reference_signal_peptides_uniprot.tsv``, release 2026_04) and
  DeepTMHMM v1 by exact sequence for member and comparator (``topology.TOOL_V1``, pinned
  snapshot ``topology.SNAPSHOT_KEY``), TM count and N-terminal orientation;
* Pfam 39.0 (gathering tier) families on member and comparator, and the families lost;
* evidence: APPRIS label and MANE flag of the member's transcript (APPRIS e115v50),
  PeptideAtlas build 607 peptides that distinguish member from canonical (genes in the pinned
  pool only), GENCODE Ribo-seq ORFs contained in the member, the isoform cohort tier;
* orthologs: identity and coverage against the human L1 reference, N-terminal extension
  beyond the human protein, signal peptide present in human but absent in the ortholog;
* structure: the model accession and whether it is our protein (``proteoform_match.tsv``,
  joined on ``candidate_id``: ``identical``, ``contained`` -- the model is longer -- or
  ``different``; no row means no model was matched);
* duplicates: other active L1 (and L2) members with the identical protein.

Flags fall in two classes (see ``FLAGS``). ``apply`` flags are the objective defects the
overnight plan lets a session apply without Becca (2026-11-11, OVERNIGHT-L1-2026-11-12.md
"Hard rules"): a non-ATG start, an internal stop, a CDS whose translation differs from the
stated protein, a protein that duplicates another active member, or a failed isoform cohort
gate. Every other flag is ``proposed``: it marks a member for hand review, and the thresholds
below are screening thresholds, not decisions. ``info`` flags are reported but route nowhere.

    uv run --locked --group domains python \
        library/l1/orf-design/audit_l1_proteoforms.py [--cache DIR]
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(HERE))

import library_pipeline as lp  # noqa: E402
import pandas as pd  # noqa: E402
import review_ortholog_correspondence as roc  # noqa: E402

from deliverome_analysis.cloning import library_design as ld  # noqa: E402
from deliverome_analysis.cloning import seqtools  # noqa: E402
from deliverome_analysis.core import s3  # noqa: E402
from deliverome_analysis.surfaceome import topology  # noqa: E402

MANIFEST = HERE / "candidate_manifest_tranches.tsv"
RENUMBERING = HERE / "l1_member_renumbering.tsv"
STRUCTURES = HERE / "proteoform_match.tsv"
GATE = HERE / "isoform_gate.tsv"
UNIPROT_SP = HERE / "reference_signal_peptides_uniprot.tsv"
RESCUE_CDS = REPO / "library/l1/proteoform-selection/rescue_cds.json"
#: The 25 expression-rescue arms are ordered as an Igk leader (+/- recoded receptor), so their
#: ordered protein is the translated ``construct`` here, not the bare receptor CDS.
RESCUE_CONSTRUCTS = REPO / "library/l1/proteoform-selection/rescue_constructs.tsv"
#: Declared in-frame N-terminal trims (5' CDS trims, no codon edited): the reference-fix
#: ledger for human references, the ortholog restart ledger for cyno orthologs.
ORTHOLOG_START_TRIMS = HERE / "ortholog_start_trims.tsv"
#: The current reference-ladder record (library_pipeline.REFERENCE_RECORD), not the
#: pre-fix ``library-design/reference_selection.tsv``. Both hashes (it and ``ld.CDS_KEY``) are
#: pinned in mane_gtf_changes_2026-11-22_summary.json, the summary of the run that published
#: the 2026-11-22-mane-gtf records (``apply_mane_gtf.py``). Before that rebuild they came from
#: reference_changes_summary.json and ortholog_start_trims_summary.json.
REFSEL_KEY = f"sprint4/s4e9-orf-design-fg-panel/{lp.REFERENCE_RECORD}reference_selection.tsv"
#: Records that explain a reference built from something other than its ladder rung.
REFERENCE_OVERRIDES = REPO / "library/l1/proteoform-selection/reference_overrides.tsv"
REFERENCE_CHANGES = HERE / "reference_changes.tsv"
REFSEL_SHA_FROM = HERE / "cut_2026_11_23_membership_changes.json"
CDS_SHA_FROM = HERE / "cut_2026_11_23_membership_changes.json"

OUT = HERE / "l1_proteoform_audit.tsv"
SUMMARY = HERE / "l1_proteoform_audit_summary.json"
INPUTS_OUT = HERE / "l1_proteoform_audit_inputs.json"

PIPELINE_INPUTS = [
    "appris_homo_sapiens_e115v50.txt",
    "Homo_sapiens.GRCh38.pep.all.fa.gz",
    "uniprot_human_reviewed_isoforms.tsv.gz",
    "Ribo-seq_ORFs.all.faa",
    "peptideatlas_build607_pools.json",
]
#: Ensembl 116 ortholog peptides, pinned in compara116_cds_inputs.json.
ORTHOLOG_PEP = {
    "mouse": "Mus_musculus.GRCm39.pep.all.fa.gz",
    "cyno": "Macaca_fascicularis.Macaca_fascicularis_6.0.pep.all.fa.gz",
}

# Screening thresholds. They route a member to hand review; none of them decides anything.
FRAGMENT_FRACTION = 0.5  # member shorter than half its comparator
FRAGMENT_MIN_AA = 100  # or shorter than 100 aa when the comparator is >= 200 aa
TRUNCATION_AA = 20  # comparator residues missing at one end
REFERENCE_RETAINED = 0.9  # library reference keeping < 90% of its canonical
#: A member more than this fraction longer *or* shorter than its UniProt canonical. The
#: truncation and retention flags only see canonical residues the member lacks, so a member
#: much longer than the canonical (KCNN2 791 vs 579 aa) passed unflagged; this column
#: records the length difference in both directions. Output only: no flag or class changes.
LENGTH_DIFF_FRACTION = 0.10
ORTHOLOG_IDENTITY = 0.6  # aligned identity to the human reference
ORTHOLOG_COVERAGE = 0.8  # fraction of the human reference aligned
ORTHOLOG_N_EXTENSION = 10  # ortholog residues before the human Met1
RIBOSEQ_MIN_AA = 8

#: flag -> (class, sets it is evaluated on, meaning)
REFS = ("library_reference", "expression_rescue")
FLAGS = {
    "non_atg_start": ("apply", "all", "CDS does not start with ATG"),
    "internal_stop": ("apply", "all", "stop codon inside the CDS"),
    "cds_frame": ("apply", "all", "CDS length is not a multiple of 3"),
    "translation_mismatch": (
        "apply",
        "all",
        "translated CDS differs from the stated protein beyond declared edits",
    ),
    "duplicate_protein_l1": ("apply", "all", "identical protein to another active L1 member"),
    "cohort_gate_fail": ("apply", ("isoform",), "isoform fails the 2026-11-12 cohort gate"),
    "fragment": ("proposed", "all", "member is a fragment of its comparator"),
    "sp_lost_vs_canonical": (
        "proposed",
        (*REFS, "isoform"),
        "UniProt canonical has a signal peptide the member lacks",
    ),
    "n_terminal_truncation": ("proposed", REFS, "comparator N-terminal residues missing"),
    "c_terminal_truncation": ("proposed", REFS, "comparator C-terminal residues missing"),
    "reference_short": ("proposed", REFS, "library reference keeps < 90% of its canonical"),
    "tm_count_change": ("proposed", (*REFS, "ortholog"), "v1 TM count differs from comparator"),
    "orientation_change": ("proposed", "all", "v1 N-terminal orientation differs from comparator"),
    "tm_all_lost": ("proposed", ("isoform",), "comparator has TM helices, the isoform none"),
    "pfam_domain_lost": (
        "proposed",
        (*REFS, "ortholog"),
        "a Pfam family on the comparator is absent",
    ),
    "pfam_all_lost": ("proposed", ("isoform",), "no comparator Pfam family left on the isoform"),
    "ortholog_divergent": ("proposed", ("ortholog",), "low identity or coverage to human"),
    "ortholog_n_extension": ("proposed", ("ortholog",), "N-terminus extends beyond human Met1"),
    "ortholog_sp_missing": (
        "proposed",
        ("ortholog",),
        "human reference has a v1 signal peptide, the ortholog none",
    ),
    "ladder_source_unrecorded": (
        "proposed",
        ("library_reference",),
        "built protein is not the ladder's selected source and no override or fix records why",
    ),
    "no_v1_topology": ("info", "all", "no DeepTMHMM v1 row for this exact sequence"),
    "structure_not_ours": ("info", "all", "no model of this exact protein"),
    "duplicate_protein_l2": ("info", "all", "identical protein to an active L2 member"),
    "reference_not_mane_appris_canonical": (
        "info",
        ("library_reference",),
        "member equals none of MANE, APPRIS principal, UniProt canonical",
    ),
}


def sha(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


def file_sha(path: Path) -> str:
    return roc.sha256_file(Path(path))


def read_tsv(path: Path) -> list[dict]:
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_fasta_gz(path: Path) -> dict[str, str]:
    """Ensembl peptide FASTA -> {unversioned transcript id: protein}."""
    out: dict[str, str] = {}
    tx, buf = None, []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if tx:
                    out[tx] = "".join(buf)
                m = re.search(r"transcript:(\S+)", line)
                tx = m.group(1).split(".")[0] if m else None
                buf = []
            else:
                buf.append(line.strip())
    if tx:
        out[tx] = "".join(buf)
    return out


def riboseq_by_gene(path: Path) -> dict[str, list[str]]:
    out: dict[str, list[str]] = defaultdict(list)
    name, buf = None, []

    def flush():
        if name:
            seq = "".join(buf).rstrip("*")
            if len(seq) >= RIBOSEQ_MIN_AA:
                out[name].append(seq)

    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                flush()
                parts = line[1:].split()[0].split(":")
                name = parts[1] if len(parts) > 1 else None
                buf = []
            else:
                buf.append(line.strip())
    flush()
    return out


def v1_summary(row: dict | None) -> dict:
    """SP length, TM count and N-terminal orientation from a v1 per-residue string."""
    if not row or not row.get("found"):
        return {"found": False, "sp": "", "tm": "", "nterm": "", "label": ""}
    s = row.get(topology.TOOL_COLUMN[topology.TOOL_V1]) or ""
    sp = len(s) - len(s.lstrip("S"))
    tm = len(re.findall(r"M+", s)) + len(re.findall(r"B+", s))
    nterm = next((c for c in s[sp:] if c in "IO"), "")
    if sp and nterm == "O":
        nterm = "O(after SP)"
    return {"found": True, "sp": sp, "tm": tm, "nterm": nterm, "label": row.get("label", "")}


def length_vs_canonical(member_aa: int, canonical_aa: int) -> tuple[float | str, int]:
    """``(|member - canonical| / canonical, 1 if that exceeds LENGTH_DIFF_FRACTION)``.

    Symmetric: a member 15% longer and one 15% shorter are both flagged. Without a canonical
    length the fraction is blank and the flag 0.
    """
    if not member_aa or not canonical_aa:
        return "", 0
    frac = abs(member_aa - canonical_aa) / canonical_aa
    return round(frac, 3), int(frac > LENGTH_DIFF_FRACTION)


def alignment(aligner, ref: str, member: str) -> dict:
    """Coverage, identity and end gaps of member against comparator ``ref``."""
    if not ref or not member:
        return {}
    if ref == member:
        return {
            "retained": 1.0,
            "identity": 1.0,
            "ref_n_missing": 0,
            "ref_c_missing": 0,
            "member_n_extra": 0,
            "start_in_ref": 1,
        }
    mapping = roc.align_map(aligner, ref, member)
    aligned = [(i, j) for i, j in enumerate(mapping) if j is not None]
    if not aligned:
        return {"retained": 0.0, "identity": 0.0}
    same = sum(ref[i] == member[j] for i, j in aligned)
    first_i, first_j = aligned[0]
    last_i, _ = aligned[-1]
    return {
        "retained": round(len(aligned) / len(ref), 3),
        "identity": round(same / len(aligned), 3),
        "ref_n_missing": first_i,
        "ref_c_missing": len(ref) - 1 - last_i,
        "member_n_extra": first_j,
        "start_in_ref": first_i + 1 - first_j,
    }


def pfam_families(hits: list[dict]) -> Counter:
    return Counter(d["pfam"] for d in hits if d["tier"] == "gathering")


def fmt_pfam(hits: list[dict]) -> str:
    return ";".join(
        f"{d['pfam']}:{d['name']}:{d['start']}-{d['end']}" for d in hits if d["tier"] == "gathering"
    )


def stage(cache: Path) -> dict[str, Path]:
    """Every pinned input, checksum-verified. Returns name -> local path."""
    cache.mkdir(parents=True, exist_ok=True)
    paths = lp.stage_inputs(PIPELINE_INPUTS, cache / "pipeline")
    spec = json.loads((HERE / "compara116_cds_inputs.json").read_text())["files"]
    for name in ORTHOLOG_PEP.values():
        entry = spec[name]
        bucket, _, key = entry["source"].partition(":")
        dest = cache / "e116" / name
        dest.parent.mkdir(parents=True, exist_ok=True)
        if not dest.exists() or file_sha(dest) != entry["sha256"]:
            getattr(s3, bucket).download(key, dest)
        if file_sha(dest) != entry["sha256"]:
            raise SystemExit(f"{name}: checksum does not match compara116_cds_inputs.json")
        paths[name] = dest
    pinned = (
        ("reference_selection.tsv", REFSEL_KEY, REFSEL_SHA_FROM, "reference_selection.tsv"),
        ("candidate_cds.json", ld.CDS_KEY, CDS_SHA_FROM, "cds.json"),
    )
    for name, key, sha_from, local in pinned:
        want = json.loads(sha_from.read_text())["record_sha256"][name]
        got = Path(s3.processed.download(key, cache / local))
        if file_sha(got) != want:
            raise SystemExit(f"{key}: checksum does not match {sha_from.name}")
        paths[name] = got
    topo = Path(s3.external.download(topology.SNAPSHOT_KEY, cache / "topology.parquet"))
    if file_sha(topo) != topology.SNAPSHOT_SHA256:
        raise SystemExit(
            f"{topology.SNAPSHOT_KEY}: checksum does not match topology.SNAPSHOT_SHA256"
        )
    paths["topology.parquet"] = topo
    paths["pfam"] = roc.fetch("pfam", cache / "pfam")
    paths["pfam_clans"] = roc.fetch("pfam_clans", cache / "pfam")
    return paths


def run_pfam(sequences: dict[str, str], paths: dict, cache: Path) -> dict[str, list[dict]]:
    """Pfam hits per sequence hash, cached on disk by hash so reruns only scan new sequences."""
    store = cache / "pfam_hits_by_sha.json"
    known = json.loads(store.read_text()) if store.exists() else {}
    todo = {h: s for h, s in sequences.items() if h not in known}
    if todo:
        clans = roc.load_clans(paths["pfam_clans"])
        print(f"pfam: scanning {len(todo)} sequences", flush=True)
        known |= roc.pfam_scan(todo, paths["pfam"], clans)
        store.write_text(json.dumps(known))
    return {h: known[h] for h in sequences}


def declared_trims() -> dict[str, int]:
    """Key -> residues trimmed from the stated protein's N-terminus, for applied trims only."""
    import apply_reference_fixes as fixes

    out = {k: int(v["codons"]) for k, v in fixes.TRIMS.items()}
    if ORTHOLOG_START_TRIMS.exists():
        with ORTHOLOG_START_TRIMS.open() as fh:
            for t in csv.DictReader(fh, delimiter="\t"):
                if t["status"] == "applied":
                    out[t["key"]] = int(t["restart_codon"]) - 1
    return out


def _write_inputs(paths: dict, cache: Path) -> None:
    """Record every input the audit read, with its SHA-256, in INPUTS_OUT."""
    inputs = {
        name: {"path_in_cache": str(p.relative_to(cache)), "sha256": file_sha(p)}
        for name, p in paths.items()
    }
    inputs["topology_snapshot"] = {"key": f"external://{topology.SNAPSHOT_KEY}"}
    inputs["cds_record"] = {
        "key": f"processed://{ld.CDS_KEY}",
        "sha256": inputs["candidate_cds.json"]["sha256"],
    }
    for local in (
        MANIFEST,
        RENUMBERING,
        STRUCTURES,
        GATE,
        UNIPROT_SP,
        RESCUE_CDS,
        RESCUE_CONSTRUCTS,
    ):
        inputs[str(local.relative_to(REPO))] = {"sha256": file_sha(local)}
    INPUTS_OUT.write_text(json.dumps(inputs, indent=1) + "\n")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cache", type=Path, default=Path.home() / ".cache/s4e9-proteoform-audit")
    args = ap.parse_args(argv)
    paths = stage(args.cache)

    # ---- members --------------------------------------------------------------------
    manifest = read_tsv(MANIFEST)
    active = [r for r in manifest if ld.is_active(r) and r["library"] in ("l1", "l2")]
    l1 = [r for r in active if r["library"] == "l1"]
    cds_rec = json.loads(paths["candidate_cds.json"].read_text())
    rescue = json.loads(RESCUE_CDS.read_text())
    rescue_arms = {r["symbol"]: r for r in read_tsv(RESCUE_CONSTRUCTS)}
    member_ids = {r["design_id"]: r["candidate_id"] for r in read_tsv(RENUMBERING)}

    def cds_of(row) -> tuple[str, dict, str]:
        e = cds_rec.get(row["key"])
        if e:
            return e["cds"], e, "candidate_cds"
        if row["key"].startswith("rescue:"):
            gene = ld.gene_from_key(row["key"])
            arm = rescue_arms[gene]
            return arm["construct"], {**rescue[gene], "arm": arm["arm"]}, "rescue_constructs"
        return "", {}, "none"

    proteins: dict[str, str] = {}
    for r in active:
        dna, _, _ = cds_of(r)
        if dna:
            proteins[r["design_id"]] = seqtools.translate(dna).rstrip("*")

    # ---- reference data ---------------------------------------------------------------
    refsel = {r["symbol"]: r for r in read_tsv(paths["reference_selection.tsv"])}
    e115 = read_fasta_gz(paths["Homo_sapiens.GRCh38.pep.all.fa.gz"])
    e116 = {sp: read_fasta_gz(paths[f]) for sp, f in ORTHOLOG_PEP.items()}
    with gzip.open(paths["uniprot_human_reviewed_isoforms.tsv.gz"], "rt") as fh:
        uniprot = {r["Entry"]: r["Sequence"] for r in csv.DictReader(fh, delimiter="\t")}
    uniprot_sp = {r["accession"]: int(r["sp_end"]) for r in read_tsv(UNIPROT_SP)}
    appris: dict[str, tuple[str, str]] = {}
    with open(paths["appris_homo_sapiens_e115v50.txt"]) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            appris[r["Transcript ID"]] = (r["APPRIS Annotation"], r["MANE"])
    atlas = json.loads(paths["peptideatlas_build607_pools.json"].read_text())
    ribo = riboseq_by_gene(paths["Ribo-seq_ORFs.all.faa"])
    gate = {r["transcript"]: r for r in read_tsv(GATE) if r["transcript"]}
    structures = {r["candidate_id"]: r for r in read_tsv(STRUCTURES)}

    def pep(tx: str) -> str:
        return e115.get(tx.split(".")[0], "")

    def canonical_of(gene: str) -> tuple[str, str]:
        acc = refsel.get(gene, {}).get("current_accession", "")
        return acc, uniprot.get(acc, "")

    # The human L1 library reference per gene: comparator for isoforms and orthologs.
    l1_ref = {
        ld.gene_from_key(r["key"]): r["design_id"]
        for r in l1
        if r["set"] == "library_reference" and r["design_id"] in proteins
    }

    # ---- duplicates -----------------------------------------------------------------
    by_protein: dict[str, list[dict]] = defaultdict(list)
    for r in active:
        if r["design_id"] in proteins:
            by_protein[proteins[r["design_id"]]].append(r)

    # ---- topology (v1 by exact sequence, per species) --------------------------------
    comparators: dict[str, tuple[str, str, str]] = {}  # design_id -> (label, seq, species)
    for r in l1:
        gene = ld.gene_from_key(r["key"])
        acc, can = canonical_of(gene)
        if r["set"] in REFS:
            comparators[r["design_id"]] = (f"uniprot:{acc}", can, "human")
        elif gene in l1_ref:
            comparators[r["design_id"]] = (
                f"l1_reference:{l1_ref[gene]}",
                proteins[l1_ref[gene]],
                "human",
            )
        else:
            comparators[r["design_id"]] = (f"uniprot:{acc}", can, "human")

    topo: dict[str, dict] = {}
    by_species: dict[str, dict[str, str]] = defaultdict(dict)
    for r in l1:
        if r["design_id"] in proteins:
            by_species[r["species"]][sha(proteins[r["design_id"]])] = proteins[r["design_id"]]
    for _, seq, sp in comparators.values():
        if seq:
            by_species[sp][sha(seq)] = seq
    snapshot = pd.read_parquet(paths["topology.parquet"])
    for sp, seqs in by_species.items():
        rep = topology.resolve(
            seqs, tool=topology.TOOL_V1, cohorts=topology.cohorts_for_species(sp), snapshot=snapshot
        )
        for rec in rep.to_dict("records"):
            topo[(sp, rec["key"])] = rec

    # ---- Pfam -------------------------------------------------------------------------
    pf_seqs = {sha(s): s for s in proteins.values()}
    pf_seqs |= {sha(s): s for _, s, _ in comparators.values() if s}
    pfam = run_pfam(pf_seqs, paths, args.cache)

    aligner = roc.make_aligner()
    trims = declared_trims()
    overridden = {o["gene"] for o in read_tsv(REFERENCE_OVERRIDES)}
    fixed = {c["key"] for c in read_tsv(REFERENCE_CHANGES)}
    rows = []
    for r in l1:
        did = r["design_id"]
        gene = ld.gene_from_key(r["key"])
        sp = r["species"]
        dna, entry, cds_source = cds_of(r)
        prot = proteins.get(did, "")
        rs = refsel.get(gene, {})
        acc, can = canonical_of(gene)
        mane = pep(rs.get("mane_enst", "")) if rs.get("mane_enst") else ""
        appris_txs = [t for t in (rs.get("appris_enst") or "").split(";") if t]
        appris_seqs = {pep(t) for t in appris_txs} - {""}
        out = {
            "design_id": did,
            "candidate_id": r["candidate_id"],
            "member_id": member_ids.get(did, ""),
            "key": r["key"],
            "candidate_set": r["set"],
            "species": sp,
            "gene_symbol": gene,
            "transcript": r["transcript"] or entry.get("transcript_id", ""),
            "route": r["route"],
            "panel": r["panel"],
            "reference_decision": rs.get("decision", ""),
            "ladder_selected_source": rs.get("selected_source", ""),
            "mane_enst": rs.get("mane_enst", ""),
            "appris_enst": rs.get("appris_enst", ""),
            "mane_appris_same_protein": (
                str(bool(mane) and mane in appris_seqs) if mane and appris_seqs else ""
            ),
            "member_equals": ";".join(
                n
                for n, ok in (
                    ("MANE", bool(mane) and prot == mane),
                    ("APPRIS", prot in appris_seqs),
                    ("UniProt_canonical", bool(can) and prot == can),
                )
                if ok
            ),
            "cds_source": cds_source,
            "cds_nt": len(dna),
        }
        # CDS checks.
        codons = [dna[i : i + 3] for i in range(0, len(dna) - len(dna) % 3, 3)]
        out["start_codon"] = codons[0] if codons else ""
        out["start_codon_replaced"] = entry.get("start_codon_replaced", "")
        out["stop_added"] = entry.get("stop_added", "")
        out["terminal_stop"] = bool(codons) and codons[-1] in ("TAA", "TAG", "TGA")
        out["internal_stops"] = sum(c in ("TAA", "TAG", "TGA") for c in codons[:-1])
        out["frame_ok"] = len(dna) % 3 == 0
        if r["set"] == "expression_rescue":
            stated, stated_src = entry.get("protein", ""), "rescue_cds protein"
        elif sp == "human":
            stated, stated_src = pep(out["transcript"]), "Ensembl 115 pep"
        else:
            stated, stated_src = (
                e116.get(sp, {}).get(out["transcript"].split(".")[0], ""),
                "Ensembl 116 pep",
            )
        stated = stated.rstrip("*")
        out["stated_protein_source"] = stated_src
        if not stated:
            tcheck = "transcript_not_in_release"
        elif prot == stated:
            tcheck = "exact"
        elif r["set"] == "expression_rescue" and prot.endswith(stated[1:]):
            tcheck = f"exact_after_{entry.get('arm', 'leader').replace(' ', '')}"
        elif r["key"] in trims and prot == stated[trims[r["key"]] :]:
            tcheck = f"exact_after_declared_trim_{trims[r['key']]}"
        elif len(prot) == len(stated):
            diffs = sum(a != b for a, b in zip(prot, stated, strict=True))
            declared = int(r.get("variant_residues") or 0) + (
                1 if out["start_codon_replaced"] else 0
            )
            tcheck = f"{diffs}_substitutions" + ("_declared" if diffs <= declared else "")
        else:
            tcheck = f"length_{len(prot)}_vs_{len(stated)}"
        out["translation_check"] = tcheck
        # Lengths and alignment.
        clabel, cseq, _ = comparators[did]
        out["protein_aa"] = len(prot)
        out["sequence_sha256"] = sha(prot) if prot else ""
        out["canonical_acc"] = acc
        out["canonical_aa"] = len(can) if can else ""
        out["mane_aa"] = len(mane) if mane else ""
        (
            out["canonical_length_diff_fraction"],
            out["length_differs_from_canonical"],
        ) = length_vs_canonical(len(prot), len(can) if can else 0)
        out["comparator"] = clabel
        out["comparator_aa"] = len(cseq) if cseq else ""
        al = alignment(aligner, cseq, prot)
        out["comparator_retained"] = al.get("retained", "")
        out["comparator_identity"] = al.get("identity", "")
        out["comparator_n_missing"] = al.get("ref_n_missing", "")
        out["comparator_c_missing"] = al.get("ref_c_missing", "")
        out["member_n_extra"] = al.get("member_n_extra", "")
        out["start_in_comparator"] = al.get("start_in_ref", "")
        if mane and prot and clabel.startswith("uniprot") and mane != cseq:
            out["mane_retained"] = alignment(aligner, mane, prot).get("retained", "")
        else:
            out["mane_retained"] = ""
        # Signal peptide and topology.
        sp_end = uniprot_sp.get(acc)
        out["uniprot_sp_end"] = sp_end or ""
        out["member_has_canonical_sp"] = (
            "" if not (sp_end and can and sp == "human") else str(can[1:sp_end] in prot)
        )
        mv = v1_summary(topo.get((sp, sha(prot)))) if prot else v1_summary(None)
        cv = v1_summary(topo.get(("human", sha(cseq)))) if cseq else v1_summary(None)
        for p, v in (("member", mv), ("comparator", cv)):
            out[f"{p}_v1_found"] = v["found"]
            out[f"{p}_v1_sp"] = v["sp"]
            out[f"{p}_v1_tm"] = v["tm"]
            out[f"{p}_v1_nterm"] = v["nterm"]
        if mv["found"]:
            hit = topo[(sp, sha(prot))]
            out["member_v1_version"] = f"{hit.get('cohort', '')}/{hit.get('topology_version', '')}"
        else:
            out["member_v1_version"] = ""
        # Pfam.
        mp = pfam.get(sha(prot), []) if prot else []
        cp = pfam.get(sha(cseq), []) if cseq else []
        lost = pfam_families(cp) - pfam_families(mp)
        out["member_pfam"] = fmt_pfam(mp)
        out["comparator_pfam"] = fmt_pfam(cp)
        out["pfam_lost"] = ";".join(f"{k}x{n}" if n > 1 else k for k, n in sorted(lost.items()))
        out["pfam_comparator_families"] = len(pfam_families(cp))
        out["pfam_families_kept"] = len(pfam_families(cp)) - sum(
            1 for k in pfam_families(cp) if k not in pfam_families(mp)
        )
        # Evidence.
        tx = out["transcript"].split(".")[0]
        label, is_mane = appris.get(tx, ("", ""))
        out["appris_label"] = label
        out["transcript_is_mane"] = is_mane
        pool = atlas["index"].get(acc)
        if pool and sp == "human" and can and prot:
            peps = atlas["pools"][pool]
            in_m = {p for p in peps if p in prot}
            in_c = {p for p in peps if p in can}
            out["atlas_peptides_member"] = len(in_m)
            out["atlas_peptides_member_only"] = len(in_m - in_c)
            out["atlas_peptides_canonical_only"] = len(in_c - in_m)
        else:
            out["atlas_peptides_member"] = out["atlas_peptides_member_only"] = ""
            out["atlas_peptides_canonical_only"] = "not in pinned pool"
        out["riboseq_orfs_in_member"] = (
            sum(o in prot for o in ribo.get(gene, [])) if sp == "human" and prot else ""
        )
        g = gate.get(tx)
        out["isoform_cohort_tier"] = g["cohort_tier"] if g else ""
        out["cohort_gate_verdict"] = g["verdict"] if g and r["set"] == "isoform" else ""
        out["cohort_gate_reason"] = g["reason"] if g and r["set"] == "isoform" else ""
        # Ortholog checks.
        if r["set"] == "ortholog":
            out["ortholog_identity_to_human"] = al.get("identity", "")
            out["ortholog_human_coverage"] = al.get("retained", "")
            out["ortholog_n_extension"] = al.get("member_n_extra", "")
        else:
            out["ortholog_identity_to_human"] = out["ortholog_human_coverage"] = ""
            out["ortholog_n_extension"] = ""
        # Structure.
        st = structures.get(r["candidate_id"])
        out["structure_model"] = st["accession"] if st else ""
        out["structure_match"] = st["verdict"] if st else "no model matched"
        out["structure_is_our_protein"] = str(bool(st) and st["verdict"] == "identical")
        # Duplicates.
        same = [x for x in by_protein.get(prot, []) if x["design_id"] != did]
        out["duplicate_l1"] = ";".join(x["design_id"] for x in same if x["library"] == "l1")
        out["duplicate_l2"] = ";".join(x["design_id"] for x in same if x["library"] == "l2")

        # ---- flags ----------------------------------------------------------------------
        flags = []

        def flag(name: str, cond: bool, flags=flags, member_set=r["set"]) -> None:
            _cls, sets, _ = FLAGS[name]
            if cond and (sets == "all" or member_set in sets):
                flags.append(name)

        flag("non_atg_start", out["start_codon"] != "ATG")
        flag("internal_stop", out["internal_stops"] > 0)
        flag("cds_frame", not out["frame_ok"])
        flag(
            "translation_mismatch",
            tcheck not in ("exact", "transcript_not_in_release")
            and not tcheck.endswith("_declared")
            and not tcheck.startswith("exact_after_"),
        )
        flag("duplicate_protein_l1", bool(out["duplicate_l1"]))
        flag("cohort_gate_fail", out["cohort_gate_verdict"] == "fail")
        clen = len(cseq) if cseq else 0
        flag(
            "fragment",
            bool(clen)
            and (
                len(prot) < FRAGMENT_FRACTION * clen
                or (clen >= 200 and len(prot) < FRAGMENT_MIN_AA)
            ),
        )
        flag("sp_lost_vs_canonical", out["member_has_canonical_sp"] == "False")
        flag("n_terminal_truncation", (al.get("ref_n_missing") or 0) >= TRUNCATION_AA)
        flag("c_terminal_truncation", (al.get("ref_c_missing") or 0) >= TRUNCATION_AA)
        flag("reference_short", al.get("retained", 1) < REFERENCE_RETAINED)
        both = mv["found"] and cv["found"] and prot != cseq
        flag("tm_count_change", both and mv["tm"] != cv["tm"])
        flag(
            "orientation_change",
            both and mv["nterm"][:1] != cv["nterm"][:1] and bool(mv["nterm"] and cv["nterm"]),
        )
        flag("tm_all_lost", both and cv["tm"] > 0 and mv["tm"] == 0)
        flag("pfam_domain_lost", bool(lost) and prot != cseq)
        flag(
            "pfam_all_lost",
            out["pfam_comparator_families"] > 0 and out["pfam_families_kept"] == 0,
        )
        flag(
            "ortholog_divergent",
            bool(al) and (al["identity"] < ORTHOLOG_IDENTITY or al["retained"] < ORTHOLOG_COVERAGE),
        )
        flag("ortholog_n_extension", (al.get("member_n_extra") or 0) >= ORTHOLOG_N_EXTENSION)
        flag(
            "ortholog_sp_missing",
            mv["found"] and cv["found"] and cv["sp"] > 0 and mv["sp"] == 0,
        )
        flag("no_v1_topology", not mv["found"])
        flag("structure_not_ours", out["structure_is_our_protein"] != "True")
        flag("duplicate_protein_l2", bool(out["duplicate_l2"]))
        flag("reference_not_mane_appris_canonical", not out["member_equals"])
        src, eq = out["ladder_selected_source"], out["member_equals"].split(";")
        built_as_ladder = {
            "mane": "MANE" in eq,
            "appris": "APPRIS" in eq,
            "uniprot_canonical": "UniProt_canonical" in eq,
            "named": gene in overridden,
        }.get(src, False)
        flag("ladder_source_unrecorded", not built_as_ladder and r["key"] not in fixed)

        out["flags"] = ";".join(flags)
        classes = {FLAGS[f][0] for f in flags}
        out["flag_class"] = (
            "apply"
            if "apply" in classes
            else "proposed"
            if "proposed" in classes
            else "info"
            if classes
            else "clean"
        )
        out["decided_by"] = "no decision taken (audit only)"
        out["implemented_by"] = "Claude Code session l1-proteoform-audit"
        out["reason"] = (
            "Becca 2026-11-11: review each L1 sequence and record overrides or "
            "modifications to MANE or APPRIS"
        )
        rows.append(out)

    rows.sort(key=lambda x: x["design_id"])
    with OUT.open("w", newline="") as fh:
        w = csv.DictWriter(fh, list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    per_flag = Counter(f for x in rows for f in x["flags"].split(";") if f)
    summary = {
        "members": len(rows),
        "by_set": dict(Counter(x["candidate_set"] for x in rows)),
        "by_flag_class": dict(Counter(x["flag_class"] for x in rows)),
        "by_flag": {f: per_flag.get(f, 0) for f in FLAGS},
        "flag_definitions": {
            f: {"class": c, "sets": s, "meaning": m} for f, (c, s, m) in FLAGS.items()
        },
        "length_differs_from_canonical": sum(x["length_differs_from_canonical"] for x in rows),
        "by_flag_and_set": {
            f: dict(Counter(x["candidate_set"] for x in rows if f in x["flags"].split(";")))
            for f in FLAGS
        },
        "thresholds": {
            "fragment_fraction": FRAGMENT_FRACTION,
            "fragment_min_aa": FRAGMENT_MIN_AA,
            "truncation_aa": TRUNCATION_AA,
            "reference_retained": REFERENCE_RETAINED,
            "length_diff_fraction": LENGTH_DIFF_FRACTION,
            "ortholog_identity": ORTHOLOG_IDENTITY,
            "ortholog_coverage": ORTHOLOG_COVERAGE,
            "ortholog_n_extension": ORTHOLOG_N_EXTENSION,
        },
        "topology_tool": topology.TOOL_V1,
        "topology_snapshot": topology.SNAPSHOT_KEY,
        "cds_key": ld.CDS_KEY,
    }
    SUMMARY.write_text(json.dumps(summary, indent=1) + "\n")
    _write_inputs(paths, args.cache)
    print(json.dumps({k: summary[k] for k in ("members", "by_flag_class", "by_flag")}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
