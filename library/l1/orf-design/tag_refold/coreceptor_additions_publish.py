"""Publish the co-receptor additions' refolds and append them to every refold table.

``publish_tagged_structures.py`` rebuilds its tables from every fold held locally, so this
appends instead: the new folds' metrics (scored by ``score_tagged_folds.py`` against the queue
``tag_refold_queue_coreceptor_additions.tsv``), structures, verdicts and viewer rows, with the
same columns and the same rules as that script. Models, scores and MSAs go to a fresh, dated
prefix (write-once). Every existing row is kept byte for byte (asserted).

Implemented by Claude Code session l1-review-export-scope (subagent), 2026-11-21, under
Becca's "Yes, rebuild now".

    uv run --locked python library/l1/orf-design/tag_refold/coreceptor_additions_publish.py \\
        --folds ~/.cache/deliverome-analysis/tagged-folds/tagged-2026-11-21-coreceptor \\
        --metrics <scored metrics tsv> --msa-dir ~/.cache/deliverome-analysis/tagged-folds/msa-coreceptor

``--queue`` and ``--prefix`` publish another batch the same way, e.g. the 2026-11-22 all-rungs
batch (``tag_refold_queue_2026_11_22_all_rungs.tsv``, prefix ``…/tagged/2026-11-22-all-rungs``).
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

from publish_tagged_structures import VERDICT_COLS  # noqa: E402
from score_tagged_folds import COLABFOLD_VERSION, MSA_SOURCE, SETTINGS  # noqa: E402

PREFIX = "sprint4/s4e9-orf-design-fg-panel/predicted-structures/tagged/2026-11-21-coreceptor"
QUEUE = EXP / "tag_refold_queue_coreceptor_additions.tsv"
IMPLEMENTED = "Claude Code session l1-review-export-scope (subagent), 2026-11-21"


def tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open() as fh:
        r = csv.DictReader(fh, delimiter="\t")
        return list(r.fieldnames or []), list(r)


def append(path: Path, rows: list[dict], cols: list[str] | None = None) -> int:
    have_cols, _ = tsv(path)
    # By column name into the committed header: the scorer writes only the columns its rows
    # use (signal-peptide and TM1 metrics appear only where they apply).
    assert set(cols or have_cols) <= set(have_cols), f"{path.name}: unknown columns"
    cols = have_cols
    before = path.read_bytes()
    with path.open("a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n",
                           extrasaction="ignore")  # fmt: skip
        w.writerows(rows)
    assert path.read_bytes().startswith(before), f"{path.name}: an existing row changed"
    return len(rows)


def main() -> None:
    global PREFIX
    ap = argparse.ArgumentParser()
    ap.add_argument("--folds", required=True)
    ap.add_argument("--metrics", required=True)
    ap.add_argument("--msa-dir", required=True)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--queue", default=str(QUEUE), help="the fold queue these folds came from")
    ap.add_argument(
        "--prefix", default=PREFIX, help="fresh processed prefix for this batch (no trailing /)"
    )
    a = ap.parse_args()
    PREFIX = a.prefix
    from deliverome_analysis.s3 import processed

    folds = Path(a.folds).expanduser()
    queue = {r["fold_id"]: r for r in tsv(Path(a.queue))[1]}
    mcols, new_metrics = tsv(Path(a.metrics))
    have = {r["fold_id"] for r in tsv(EXP / "tag_refold_metrics.tsv")[1]}
    new_metrics = [r for r in new_metrics if r["fold_id"] not in have]
    metrics = {r["tagged_sha256"]: r for r in new_metrics}
    existing = set(processed.list_keys(PREFIX + "/"))

    structures = []
    for fid, q in sorted(queue.items()):
        pdb, scores = folds / f"{fid}.pdb", folds / f"{fid}.scores.json.gz"
        assert hashlib.sha256(q["tagged_seq"].encode()).hexdigest() == q["tagged_sha256"]
        for f in (pdb, scores):
            key = f"{PREFIX}/{f.name}"
            if key not in existing and not a.dry_run:
                processed.upload(f, key)
        m = metrics.get(q["tagged_sha256"], {})
        c0 = q.get("cassette_start", "")
        structures.append({
            "fold_id": fid, "kind": q["kind"], "rung": q["rung"], "design_id": q["design_id"],
            "gene_symbol": q["gene_symbol"], "tagged_sha256": q["tagged_sha256"],
            "native_sha256": q["native_sha256"],
            "cassette_start_1based": int(c0) + 1 if c0 not in ("", None) else "",
            "cassette_end_1based": q.get("cassette_end", ""),
            "insert_after_native_residue": q.get("insert_after_native_residue", ""),
            "mean_plddt": m.get("mean_plddt", ""), "verdict": m.get("verdict", ""),
            "pdb_key": f"processed://{PREFIX}/{pdb.name}",
            "scores_key": f"processed://{PREFIX}/{scores.name}",
            "pdb_sha256": hashlib.sha256(pdb.read_bytes()).hexdigest(),
            "colabfold_version": COLABFOLD_VERSION, "settings": SETTINGS, "msa_source": MSA_SOURCE,
            "implemented_by": IMPLEMENTED,
        })  # fmt: skip
    pdb_key = {s["tagged_sha256"]: s["pdb_key"] for s in structures}
    sources = [r for r in tsv(EXP / "tag_refold_sources.tsv")[1] if r["tagged_sha256"] in metrics]
    verdicts = []
    for src in sources:
        m = metrics[src["tagged_sha256"]]
        v = {**src, "pdb_key": pdb_key.get(src["tagged_sha256"], "")}
        v.update({c: m.get(c, "") for c in VERDICT_COLS})
        verdicts.append(v)
    # viewer manifest
    index = {r["native_sha256"]: r["ref"] for r in tsv(EXP / "tag_refold_native_models.tsv")[1]}
    own = {s["native_sha256"]: s["pdb_key"] for s in structures if s["kind"] == "native"}
    msa_dir = Path(a.msa_dir).expanduser()
    native_of = {q["tagged_sha256"]: q["native_sha256"] for q in queue.values()}
    viewer = []
    for v in verdicts:
        h = native_of[v["tagged_sha256"]]
        viewer.append({**v, "native_sha256": h, "ordered_protein_sha256": "",
                       "native_model_path": own.get(h) or index.get(h, ""),
                       "msa_key": f"processed://{PREFIX}/msa/{h}.a3m.gz"
                       if (msa_dir / f"{h}.a3m.gz").exists() else ""})  # fmt: skip
    if a.dry_run:
        print(f"dry run: {len(new_metrics)} metrics, {len(structures)} structures, "
              f"{len(verdicts)} verdicts")  # fmt: skip
        return
    print("metrics", append(EXP / "tag_refold_metrics.tsv", new_metrics, mcols))
    print("structures", append(EXP / "tag_refold_structures.tsv", structures))
    print("verdicts", append(EXP / "tag_refold_verdicts.tsv", verdicts))
    print("viewer", append(EXP / "tag_refold_viewer_manifest.tsv", viewer))
    tmp = EXP / "data" / "coreceptor-refold-tables" / Path(a.queue).stem
    tmp.mkdir(parents=True, exist_ok=True)
    for name, rows in (("structures.tsv", structures), ("verdicts.tsv", verdicts),
                       ("tagged_fold_metrics.tsv", new_metrics)):  # fmt: skip
        with (tmp / name).open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(rows)
        suffix = "" if Path(a.queue) == QUEUE else "_" + Path(a.queue).stem.rsplit("_", 1)[-1]
        processed.upload(tmp / name, f"{PREFIX}/{Path(name).stem}{suffix}.tsv")
    for f in sorted(msa_dir.glob("*.a3m.gz")):
        key = f"{PREFIX}/msa/{f.name}"
        if key not in existing:
            processed.upload(f, key)
    print(f"published under processed://{PREFIX}/")


if __name__ == "__main__":
    main()
