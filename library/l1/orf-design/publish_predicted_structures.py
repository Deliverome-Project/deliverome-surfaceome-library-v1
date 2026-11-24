"""Publish the ColabFold predictions to S3, so the tag-site lanes have a model to read.

`fold_missing_structures_modal.py` writes PDBs and per-model score JSONs locally; this
uploads them under the experiment prefix with the usual provenance stamp (commit,
dirty flag, script, lockfile hashes), together with the manifest saying which candidate
each accession belongs to and what its mean pLDDT was.

Predictions are NOT models from AlphaFold DB and must not be mistaken for them: they are
ours, made with ColabFold's MSA server rather than the DeepMind pipeline, so the key names
them `colabfold/` and the summary carries the version that produced them.

    uv run --locked python publish_predicted_structures.py
    uv run --locked python publish_predicted_structures.py --dry-run
"""

from __future__ import annotations

import argparse
import csv
from datetime import UTC, datetime
from pathlib import Path

from deliverome_analysis.s3 import processed

HERE = Path(__file__).parent
#: Outside the experiment folder on purpose: upload-experiment-outputs.sh walks every
#: file under an experiment, so coordinates kept there get published a second time.
LOCAL = Path.home() / ".cache" / "deliverome-analysis" / "predicted_structures"
PREFIX = "sprint4/s4e9-orf-design-fg-panel/predicted-structures/colabfold"
COLABFOLD_VERSION = "1.6.2-cuda12"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--local", type=Path, default=LOCAL)
    ap.add_argument("--date", default=datetime.now(UTC).strftime("%Y-%m-%d"))
    args = ap.parse_args(argv)

    local: Path = args.local
    if not local.exists():
        raise SystemExit(f"{local} does not exist — run the Modal job first")

    summary = local / "fold_summary.tsv"
    manifest = HERE / "fold_inputs.tsv"
    rows = (
        {r["accession"]: r for r in csv.DictReader(manifest.open(), delimiter="\t")}
        if manifest.exists()
        else {}
    )

    files = sorted(
        [
            p
            for p in local.iterdir()
            if p.suffix in (".pdb", ".json") or p.name == "fold_summary.tsv"
        ]
    )
    if not files:
        raise SystemExit(f"nothing to publish in {local}")

    base = f"{PREFIX}/{args.date}"
    print(f"publishing {len(files)} file(s) to processed://{base}/")
    for path in files:
        key = f"{base}/{path.name}"
        meta = rows.get(path.name.split(".")[0], {})
        label = f"  {path.name:<34} {path.stat().st_size / 1024:>8.1f} KiB"
        if meta:
            label += f"  {meta.get('symbol', '')}/{meta.get('species', '')}/{meta.get('set', '')}"
        print(label)
        if not args.dry_run:
            processed.upload(path, key)

    if summary.exists() and not args.dry_run:
        print(f"\ncolabfold {COLABFOLD_VERSION}; summary at processed://{base}/{summary.name}")
    if args.dry_run:
        print("\ndry run — nothing uploaded")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
