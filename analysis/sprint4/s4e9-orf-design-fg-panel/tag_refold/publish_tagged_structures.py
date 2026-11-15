"""Publish the tagged refolds, their metrics and per-rung verdicts under a dated prefix.

    processed://sprint4/s4e9-orf-design-fg-panel/predicted-structures/tagged/<date>/
        <fold_id>.pdb               rank-1 model, pLDDT in the B-factor column (write-once)
        <fold_id>.scores.json.gz    ColabFold scores: pLDDT, PAE, pTM (write-once)
        structures.tsv              one row per published model (= tag_refold_structures.tsv)
        tagged_fold_metrics.tsv     one row per fold, every metric (= tag_refold_metrics.tsv)
        verdicts.tsv                one row per rung construct (= tag_refold_verdicts.tsv)
        verdicts_<rung>.tsv         the same, one file per rung
        inputs/<sha12>-<file>       the rung inputs (written by build_refold_inputs.py)

Models are uploaded once each with ``processed.upload`` (provenance-stamped) and never
overwritten; the tables are the index and are re-published as the run grows.
Implemented by Claude Code session l1-tag-refold.

    uv run --locked python analysis/sprint4/s4e9-orf-design-fg-panel/tag_refold/publish_tagged_structures.py \\
        --folds ~/.cache/deliverome-analysis/tagged-folds/tagged-2026-11-12 --date 2026-11-12
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
EXP = HERE.parent
sys.path.insert(0, str(HERE))

from build_refold_inputs import SNAPSHOT  # noqa: E402
from refold_queue import read_queue  # noqa: E402
from score_tagged_folds import COLABFOLD_VERSION, MSA_SOURCE, SETTINGS  # noqa: E402

PREFIX = "sprint4/s4e9-orf-design-fg-panel/predicted-structures/tagged"

#: Columns carried into the verdict tables, after the rung's own ids.
VERDICT_COLS = [
    "verdict", "verdict_reasons", "lddt_confident", "lddt_site", "tm_native", "tm_confident",
    "rmsd_confident", "dplddt_native", "dplddt_site", "mean_plddt", "cassette_plddt",
    "spytag003_plddt", "alfa_plddt", "spytag003_rsa", "alfa_rsa", "cassette_contact_res",
    "epitope_confident_contact_res", "contacted_native_residues",
    "cassette_longest_helix_excl_alfa", "spytag003_longest_helix", "tm1_plddt",
    "tm1_longest_helix", "tm1_longest_helix_native", "sp_kind", "sp_state", "sp_longest_helix", "sp_plddt", "sp_note",
    "native_model", "decided_by",
]  # fmt: skip


def read(path: Path) -> list[dict[str, str]]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write(path: Path, rows: list[dict], cols: list[str]) -> None:
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(
            fh, fieldnames=cols, delimiter="\t", lineterminator="\n", extrasaction="ignore"
        )
        w.writeheader()
        w.writerows(rows)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--folds", required=True)
    ap.add_argument("--date", required=True)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--skip-msas", action="store_true", help="publish tables and models only")
    ap.add_argument(
        "--msa-dir", default=str(Path.home() / ".cache/deliverome-analysis/tagged-folds/msa")
    )
    ap.add_argument(
        "--allow-shrink",
        default="",
        help="publish even if a verdict table loses rows; give the reason (recorded in the log)",
    )
    args = ap.parse_args()

    from deliverome_analysis.s3 import processed

    folds = Path(args.folds).expanduser()
    base = f"{PREFIX}/{args.date}"
    queue = {r["fold_id"]: r for r in read_queue(EXP / "tag_refold_queue.tsv")}
    metrics = {
        r["tagged_sha256"]: r
        for r in read(EXP / "tag_refold_metrics.tsv")
        if r.get("tagged_sha256")
    }
    existing = set(processed.list_keys(base + "/"))
    s_prev = EXP / "tag_refold_structures.tsv"
    prev_pdb_sha = {r["fold_id"]: r["pdb_sha256"] for r in read(s_prev)} if s_prev.exists() else {}
    prev_pdb_sha.pop("l1-m1153.nterm_built", None)  # recorded from the refold; see README

    structures, uploaded = [], 0
    for fid, q in sorted(queue.items()):
        pdb, scores = folds / f"{fid}.pdb", folds / f"{fid}.scores.json.gz"
        if not pdb.exists() or not scores.exists():  # a collect may be mid-way; next run gets it
            continue
        assert hashlib.sha256(q["tagged_seq"].encode()).hexdigest() == q["tagged_sha256"]
        for f in (pdb, scores):
            key = f"{base}/{f.name}"
            if key not in existing and not args.dry_run:
                processed.upload(f, key)
                uploaded += 1
        # Write-once: a local model that differs from the copy already published means it was
        # refolded after publishing (IL1R2, 2026-11-12). The published object wins; stop here.
        prev_sha = prev_pdb_sha.get(fid)
        if prev_sha and prev_sha != hashlib.sha256(pdb.read_bytes()).hexdigest():
            raise SystemExit(f"{fid}: local PDB differs from the published one; restore it from S3")
        m = metrics.get(q["tagged_sha256"], {})
        c0 = q.get("cassette_start", "")
        structures.append({
            "fold_id": fid, "kind": q.get("kind", ""), "rung": q["rung"], "design_id": q.get("design_id", ""),
            "gene_symbol": q["gene_symbol"], "tagged_sha256": q["tagged_sha256"],
            "native_sha256": q["native_sha256"],
            "cassette_start_1based": int(c0) + 1 if c0 not in ("", None) else "",
            "cassette_end_1based": q.get("cassette_end", ""),
            "insert_after_native_residue": q.get("insert_after_native_residue", ""),
            "mean_plddt": m.get("mean_plddt", ""), "verdict": m.get("verdict", ""),
            "pdb_key": f"processed://{base}/{pdb.name}", "scores_key": f"processed://{base}/{scores.name}",
            "pdb_sha256": hashlib.sha256(pdb.read_bytes()).hexdigest(),
            "colabfold_version": COLABFOLD_VERSION, "settings": SETTINGS, "msa_source": MSA_SOURCE,
            "implemented_by": "Claude Code session l1-tag-refold",
        })  # fmt: skip
    s_out = EXP / "tag_refold_structures.tsv"
    write(s_out, structures, list(structures[0]))

    # one verdict row per rung construct, keyed as the rung keys it
    pdb_key = {s["tagged_sha256"]: s["pdb_key"] for s in structures}
    sources = read(EXP / "tag_refold_sources.tsv")
    verdicts = []
    for src in sources:
        m = metrics.get(src["tagged_sha256"], {})
        v = {**src, "pdb_key": pdb_key.get(src["tagged_sha256"], "")}
        v.update({c: m.get(c, "") for c in VERDICT_COLS})
        if m.get("verdict") in ("", None, "not_folded"):
            v["verdict"] = "not_folded"
        verdicts.append(v)
    cols = [*sources[0].keys(), "pdb_key", *VERDICT_COLS]
    v_out = EXP / "tag_refold_verdicts.tsv"
    write(v_out, verdicts, cols)
    tmp = SNAPSHOT / "verdicts"
    tmp.mkdir(parents=True, exist_ok=True)
    per_rung: dict[str, list] = {}
    for v in verdicts:
        per_rung.setdefault(v["rung"], []).append(v)
    for rung, rows in per_rung.items():
        write(tmp / f"verdicts_{rung}.tsv", rows, cols)

    # Viewer manifest: one row per built construct with both models, its MSA and its verdict.
    index = {r["native_sha256"]: r["ref"] for r in read(EXP / "tag_refold_native_models.tsv")}
    own = {s["native_sha256"]: s["pdb_key"] for s in structures if s["kind"] == "native"}
    msa_dir = Path(args.msa_dir).expanduser()
    native_of = {q["tagged_sha256"]: q["native_sha256"] for q in queue.values()}
    # the protein as ordered (precursor) for each folded protein, from the order-of-record
    # check: an Igk design is folded as its mature chain, so a reader joining on the ordered
    # protein's SHA-256 needs this column (l1-viewer, 2026-11-13)
    ordered: dict[str, set[str]] = {}
    oc = EXP / "tag_refold_order_check.tsv"
    fold_sha = {fid: q["tagged_sha256"] for fid, q in queue.items()}
    for r in read(oc) if oc.exists() else []:
        if r["fold_id"] in fold_sha:
            ordered.setdefault(fold_sha[r["fold_id"]], set()).add(r["protein_sha256"])
    for v in verdicts:
        h = native_of.get(v["tagged_sha256"], "")
        v["native_sha256"] = h
        v["ordered_protein_sha256"] = ";".join(sorted(ordered.get(v["tagged_sha256"], ())))
        v["native_model_path"] = own.get(h) or index.get(h, "")
        v["msa_key"] = (
            f"processed://{base}/msa/{h}.a3m.gz" if (msa_dir / f"{h}.a3m.gz").exists() else ""
        )
    m_cols = ["rung", "source_id", "design_id", "gene_symbol", "tier", "tagged_sha256",
              "ordered_protein_sha256", "pdb_key",
              "native_sha256", "native_model_path", "msa_key", "mean_plddt", "cassette_plddt",
              "verdict", "verdict_reasons", "lddt_confident", "lddt_site", "tm_confident",
              "rmsd_confident", "spytag003_rsa", "alfa_rsa", "epitope_confident_contact_res",
              "sp_state", "tm1_longest_helix", "decided_by"]  # fmt: skip
    manifest_out = EXP / "tag_refold_viewer_manifest.tsv"
    write(manifest_out, verdicts, m_cols)

    # Published verdict tables only grow: every (tagged_sha256, source_id) a previous publish
    # carried must still be here (2026-11-12, a re-ingest of an emptied request file dropped
    # 52 loop rows from verdicts_loop.tsv).
    tables = {
        "verdicts.tsv": verdicts,
        **{f"verdicts_{r}.tsv": rows for r, rows in per_rung.items()},
    }
    lost_total = 0
    for name, rows_now in tables.items():
        key = f"{base}/{name}"
        if key not in existing:
            continue
        before = read(Path(processed.download(key)))
        now_keys = {(r["tagged_sha256"], r["source_id"]) for r in rows_now}
        lost = [r for r in before if (r["tagged_sha256"], r["source_id"]) not in now_keys]
        if lost or len(rows_now) < len(before):
            lost_total += len(lost)
            print(f"  REGRESSION {name}: {len(before)} -> {len(rows_now)} rows, {len(lost)} keys lost, "
                  f"e.g. {[r['source_id'] for r in lost[:5]]}", file=sys.stderr)  # fmt: skip
    if lost_total and not args.allow_shrink:
        raise SystemExit(
            "refusing to publish verdict tables that lose rows (pass --allow-shrink REASON)"
        )
    if lost_total:
        print(f"  publishing despite {lost_total} lost rows: {args.allow_shrink}", file=sys.stderr)

    if not args.dry_run:
        processed.upload(manifest_out, f"{base}/viewer_manifest.tsv")
        processed.upload(s_out, f"{base}/structures.tsv")
        processed.upload(EXP / "tag_refold_metrics.tsv", f"{base}/tagged_fold_metrics.tsv")
        processed.upload(v_out, f"{base}/verdicts.tsv")
        for rung in per_rung:
            processed.upload(tmp / f"verdicts_{rung}.tsv", f"{base}/verdicts_{rung}.tsv")
    # MSAs last: they are large and slow to upload, and the tables above must not wait on them.
    if not args.dry_run and not args.skip_msas:
        for f in sorted(msa_dir.glob("*.a3m.gz")):
            key = f"{base}/msa/{f.name}"
            if key not in existing:
                processed.upload(f, key)
                uploaded += 1
    from collections import Counter

    print(f"{len(structures)} models listed, {uploaded} new objects under processed://{base}/")
    for (rung, v), n in sorted(Counter((v["rung"], v["verdict"]) for v in verdicts).items()):
        print(f"  {rung:<10} {v:<16} {n}")


if __name__ == "__main__":
    main()
