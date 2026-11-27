"""Build a single flat reviewer construct catalog for the L1 order of record.

Joins the order plan (`l1_order_plan.tsv`), synthesis compliance table
(`synthesis_compliance.tsv`), five rung build tables (`constructs_l1_*.tsv`),
topology/shedding table (`l1_tag_paths.tsv`), candidate manifest
(`candidate_manifest_tranches.tsv`), ColabFold refold index
(`tag_refold_structures.tsv`), and the synthesis-recoded order FASTA
(`l1_order_inserts_recoded.fasta`) into one self-contained, plate-ordered TSV
(`l1_order_catalog.tsv.gz`, 3,352 rows).

Every recoded DNA insert is checked against `synthesis_compliance.tsv`'s
`recoded_insert_sha256` and translated to verify its tagged protein sequence.

Usage:
    uv run --locked python build_l1_order_catalog.py
    uv run --locked python build_l1_order_catalog.py --plain /tmp/l1_order_catalog.tsv
    uv run --locked python build_l1_order_catalog.py --check
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path

from deliverome_analysis.cloning import plasmid_qc as qc
from deliverome_analysis.core import s3

RP = Path(__file__).resolve().parent
S4E9 = RP.parent
DEFAULT_OUT = RP / "l1_order_catalog.tsv.gz"

COLUMNS = [
    "insert_id",
    "plate",
    "well",
    "arm_id",
    "design_id",
    "gene_symbol",
    "candidate_set",
    "species",
    "source_key",
    "transcript_or_accession",
    "also_encodes",
    "rung",
    "placement",
    "insert_after_residue",
    "n_tm",
    "signal_peptide_aa",
    "native_aa",
    "tagged_aa",
    "insert_nt",
    "vector_bp",
    "synthesis_category",
    "codons_recoded",
    "refold_verdict",
    "refold_mean_plddt",
    "shed_risk",
    "decided_by",
    "reason",
    "insert_sha256",
    "tagged_protein_aa",
    "insert_dna",
]


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _read_fasta_by_arm(path: Path) -> dict[str, str]:
    seqs: dict[str, list[str]] = {}
    cur: str | None = None
    with path.open(encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                cur = line[1:].split()[0]
                seqs[cur] = []
            elif cur is not None:
                seqs[cur].append(line.upper())
    return {k: "".join(v) for k, v in seqs.items()}


def _load_rung_rows() -> dict[tuple[str, str], dict[str, str]]:
    """Index per-arm rung metadata by (rung, arm_id) or (rung, design_id)."""
    by_key: dict[tuple[str, str], dict[str, str]] = {}
    for rung, fname in (
        ("nterm", "constructs_l1_nterm.tsv"),
        ("multipass", "constructs_l1_multipass.tsv"),
        ("type3", "constructs_l1_type3.tsv"),
        ("cterm", "constructs_l1_cterm.tsv"),
        ("internal", "constructs_l1_internal.tsv"),
    ):
        for row in _read_tsv(S4E9 / fname):
            arm = row.get("arm_id") or f"{row['design_id']}-a1"
            by_key[(rung, arm)] = row
            by_key.setdefault((rung, row["design_id"]), row)
    return by_key


def _placement_and_site(rung: str, rrow: dict[str, str]) -> tuple[str, str]:
    if rung == "nterm":
        p = rrow.get("placement", "n_after_sp")
        site = rrow.get("native_after_residue") or rrow.get("cassette_after_residue") or ""
        return p, site
    if rung == "multipass":
        d = rrow.get("design", "igk_tag")
        return f"multipass_{d}", "1"
    if rung == "type3":
        d = rrow.get("design", "tag_alone")
        return f"type3_{d}", "1"
    if rung == "cterm":
        p = rrow.get("placement", "cterm")
        site = rrow.get("cassette_after_residue") or rrow.get("native_aa") or ""
        return p, site
    if rung == "internal":
        role = rrow.get("arm_role", "")
        p = rrow.get("placement", "internal")
        if role == "snorkel" or p == "snorkel":
            return "snorkel", "C-term_snorkel"
        return "internal_loop", rrow.get("cassette_after_residue", "")
    return rung, ""


def build_tsv_text() -> str:
    summary = json.loads((S4E9 / "synthesis_compliance_summary.json").read_text(encoding="utf-8"))
    fasta_key = summary["recoded_fasta"].removeprefix("processed://")
    fasta_path = Path(s3.processed.download(fasta_key))
    fasta_sha = hashlib.sha256(fasta_path.read_bytes()).hexdigest()
    if fasta_sha != summary["recoded_fasta_sha256"]:
        raise ValueError(f"FASTA SHA-256 {fasta_sha} != summary {summary['recoded_fasta_sha256']}")
    fasta = _read_fasta_by_arm(fasta_path)

    plan = _read_tsv(RP / "l1_order_plan.tsv")
    synth = {r["arm_id"]: r for r in _read_tsv(S4E9 / "synthesis_compliance.tsv")}
    manifest = {r["design_id"]: r for r in _read_tsv(S4E9 / "candidate_manifest_tranches.tsv")}
    tag_paths = {r["candidate_id"]: r for r in _read_tsv(S4E9 / "l1_tag_paths.tsv")}
    rung_rows = _load_rung_rows()

    # Index refold metrics by (design_id, tagged_sha256), then tagged_sha256, then design_id
    refold_by_sha: dict[tuple[str, str], dict[str, str]] = {}
    refold_by_prot: dict[str, dict[str, str]] = {}
    refold_by_design: dict[str, dict[str, str]] = {}
    for row in _read_tsv(S4E9 / "tag_refold_structures.tsv"):
        if row.get("kind") != "tagged":
            continue
        refold_by_sha[(row["design_id"], row["tagged_sha256"])] = row
        refold_by_prot.setdefault(row["tagged_sha256"], row)
        refold_by_design.setdefault(row["design_id"], row)

    buf = io.StringIO()
    writer = csv.DictWriter(buf, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
    writer.writeheader()

    for prow in plan:
        arm_id = prow["arm_id"]
        design_id = prow["design_id"]
        rung = prow["rung"]
        srow = synth[arm_id]
        mrow = manifest[design_id]
        trow = tag_paths.get(mrow["candidate_id"], {})
        rrow = rung_rows.get((rung, arm_id)) or rung_rows.get((rung, design_id), {})

        dna = fasta[arm_id]
        sha = hashlib.sha256(dna.encode("ascii")).hexdigest()
        if sha != srow["recoded_insert_sha256"]:
            raise ValueError(
                f"{arm_id}: FASTA SHA-256 {sha} != synthesis_compliance.tsv "
                f"{srow['recoded_insert_sha256']}"
            )
        if len(dna) != int(prow["insert_nt"]):
            raise ValueError(f"{arm_id}: FASTA length {len(dna)} != plan {prow['insert_nt']}")

        prot = qc.translate_insert(dna)
        prot_sha = hashlib.sha256(prot.encode("ascii")).hexdigest()
        rf = (
            refold_by_sha.get((design_id, prot_sha))
            or refold_by_prot.get(prot_sha)
            or refold_by_design.get(design_id, {})
        )

        placement, site = _placement_and_site(rung, rrow)
        refold_verdict = (
            rrow.get("refold_verdict") or rrow.get("order_verdict") or rf.get("verdict") or ""
        )
        writer.writerow(
            {
                "insert_id": prow["insert_id"],
                "plate": prow["plate"],
                "well": prow["well"],
                "arm_id": arm_id,
                "design_id": design_id,
                "gene_symbol": srow["gene_symbol"],
                "candidate_set": prow["set"],
                "species": mrow["species"],
                "source_key": prow["key"],
                "transcript_or_accession": mrow.get("transcript", ""),
                "also_encodes": mrow.get("also_encodes", ""),
                "rung": rung,
                "placement": placement,
                "insert_after_residue": site,
                "n_tm": trow.get("n_tm", ""),
                "signal_peptide_aa": trow.get("signal_peptide_aa", ""),
                "native_aa": rrow.get("native_aa") or trow.get("length_aa", ""),
                "tagged_aa": str(len(prot)),
                "insert_nt": prow["insert_nt"],
                "vector_bp": prow["vector_bp"],
                "synthesis_category": srow.get("after_category", ""),
                "codons_recoded": srow.get("codons_recoded", ""),
                "refold_verdict": refold_verdict,
                "refold_mean_plddt": rf.get("mean_plddt", ""),
                "shed_risk": trow.get("shed_risk", ""),
                "decided_by": rrow.get("decided_by") or srow.get("decided_by", ""),
                "reason": rrow.get("reason", ""),
                "insert_sha256": sha,
                "tagged_protein_aa": prot,
                "insert_dna": dna,
            }
        )

    return buf.getvalue()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT, help="gzipped TSV output path")
    ap.add_argument("--plain", type=Path, help="optional uncompressed TSV output path")
    ap.add_argument(
        "--check",
        action="store_true",
        help="verify committed gzipped catalog matches regenerated contents",
    )
    args = ap.parse_args()

    text = build_tsv_text()
    rows = text.count("\n") - 1
    if args.check:
        if not args.out.is_file():
            raise SystemExit(f"missing {args.out}")
        with gzip.open(args.out, "rt", encoding="utf-8") as fh:
            got = fh.read()
        if got != text:
            raise SystemExit(f"{args.out} is out of date; rerun build_l1_order_catalog.py")
        print(f"OK {args.out.name}: {rows} constructs match order of record")
        return 0

    raw = text.encode("utf-8")
    # mtime=0 makes gzip output deterministic across runs
    with (
        args.out.open("wb") as raw_fh,
        gzip.GzipFile(filename="l1_order_catalog.tsv", mode="wb", fileobj=raw_fh, mtime=0) as gz,
    ):
        gz.write(raw)
    if args.plain:
        args.plain.write_bytes(raw)
        print(f"wrote {args.plain}: {rows} rows ({len(raw) / 1e6:.2f} MB)")
    print(f"wrote {args.out}: {rows} rows ({args.out.stat().st_size / 1e6:.2f} MB gz)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
