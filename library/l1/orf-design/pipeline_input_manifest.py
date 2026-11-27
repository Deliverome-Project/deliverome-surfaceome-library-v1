"""Saved inputs, not just pins: every S3 object the L1 pipeline reads, with its SHA-256.

``library_pipeline.py verify`` rebuilds L1 offline from saved inputs. This file says exactly
which inputs those are, so a rerun years from now can show each one is still there and
unchanged. The list is measured, not hand-written: run the full ``verify`` from a fresh
clone with an empty S3 cache (``DELIVEROME_ANALYSIS_S3_CACHE`` pointing at a new
directory). The cache then holds exactly the objects the 23 stages downloaded, one file per
object, named ``{bucket}/{key}.{etag}``. ``--from-cache`` turns that directory into
``library_pipeline_inputs_saved.tsv``: bucket, key, size, ETag and SHA-256 per object.

``--check`` (pipeline stage ``saved_inputs``) fails if any listed object is missing from S3,
or its SHA-256 no longer matches. ``external://`` objects are append-only by policy; for
``processed://`` objects the check is what guarantees they are unchanged.

What the manifest cannot hold is recorded beside it (``library_pipeline_inputs_saved.json``):

* tool versions and their identity records (DeepTMHMM 1.0 academic runtime fingerprint, the
  Python and R lockfiles' SHA-256);
* the git commits the stages replay at (``synthesis_compliance`` reads the rung builds by
  commit; ``isoform_universe`` and ``ortholog_repick`` rerun at their producing commits);
* live web services the design used once and that publish no downloadable snapshot of the
  response. Their responses were saved at the time and are among the listed objects
  (MyGene.info, the UniProt REST features, the PeptideAtlas queries); the services
  themselves cannot be saved.

Measured on 2026-11-13 by the Claude Code session l1-reviewer-package (from a clean clone
at c315c44b with an empty cache: 99 objects, 23 of 23 stages reproduced, and the order plan
byte for byte); requested by the coordinator under Becca's 2026-11-12 delegation.

    uv run --locked python pipeline_input_manifest.py --from-cache <empty-cache-dir-after-verify>
    uv run --locked python pipeline_input_manifest.py --check
"""

from __future__ import annotations

import argparse
import csv
import json
import platform
import sys
from pathlib import Path

from deliverome_analysis.core.provenance import sha256_file, verify_s3_objects

HERE = Path(__file__).resolve().parent
#: How the current manifest was measured.
MEASURED = (
    "2026-11-14, clean clone at f5222462 (main with the final-c order of record, #340), empty S3 "
    "and structure-model caches, full verify: 30 of 30 stages reproduced and the order plan "
    "tables byte for byte"
)
REPO = HERE.parents[2]
MANIFEST = HERE / "library_pipeline_inputs_saved.tsv"
META = HERE / "library_pipeline_inputs_saved.json"
BUCKETS = {"deliverome-external": "external", "deliverome-processed": "processed"}
COLS = ("bucket", "key", "bytes", "etag", "sha256")


def from_cache(cache: Path) -> int:
    rows = []
    for path in sorted(p for p in cache.rglob("*") if p.is_file()):
        rel = path.relative_to(cache).as_posix()
        bucket, rest = rel.split("/", 1)
        if bucket not in BUCKETS:
            raise SystemExit(f"unexpected bucket in cache: {bucket}")
        key, etag = rest.rsplit(".", 1)
        rows.append(
            {
                "bucket": BUCKETS[bucket],
                "key": key,
                "bytes": path.stat().st_size,
                "etag": etag,
                "sha256": sha256_file(path),
            }
        )
    with MANIFEST.open("w", newline="") as fh:
        w = csv.DictWriter(fh, COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    sys.path.insert(0, str(HERE))
    import synthesis_compliance as sc

    runtime = HERE / "topology_runs/2026-11-12_cyno_restarts.3line.runtime.json"
    meta = {
        "measured": MEASURED,
        "objects": len(rows),
        "bytes": sum(r["bytes"] for r in rows),
        "by_bucket": {b: sum(r["bucket"] == b for r in rows) for b in BUCKETS.values()},
        # the interpreter of the measuring run; .python-version pins only the minor version
        "python_version": platform.python_version(),
        "lockfiles": {
            name: sha256_file(REPO / name)
            for name in ("uv.lock", "renv.lock")
            if (REPO / name).exists()
        },
        "git_refs_replayed": {
            "synthesis_compliance": sc.RUNG_REFS,
            "note": "isoform_universe and ortholog_repick rerun at the producing commits named "
            "in library_pipeline.py; a clone needs the full history, not a shallow one",
        },
        "tools": {
            "deeptmhmm_v1": {
                "version": "DeepTMHMM 1.0 academic (deeptmhmm-1.0-academic / 1.0.24)",
                "where": "Modal volume deeptmhmm-v1 (license-gated, not redistributable)",
                "identity": json.loads(runtime.read_text()).get("package_files_sha256", {}),
            },
            "signalp6": {
                "version": "SignalP 6.0, package 6.0+h (slow_sequential models)",
                "install_source": "DTU Health Tech academic download "
                "signalp-6.0i.slow_sequential.tar.gz (licensed, not redistributable); held on "
                "Modal volumes signalp6-package and signalp6-models",
                "identity": "reviewer-package/nterm/signalp6_identity.json: SHA-256 of every "
                "package and checkpoint file; tarball SHA-256 3193855d… in signalp6_modal.py",
            },
            "netgpi": {
                "version": "NetGPI 1.1",
                "install_source": "the DTU Health Tech web service (services.healthtech.dtu.dk); "
                "no local install. Each run's raw HTTP exchange and tool_version are kept with "
                "its outputs (run_ortholog_netgpi.py, run_engineered_netgpi.py)",
            },
            "colabfold": "l1-tag-refold's run records (PR #301); folds are outputs, not pipeline inputs",
            "python_packages": "uv.lock (SHA-256 above)",
            "r_packages": "renv.lock (SHA-256 above)",
        },
        "cannot_be_saved": [
            "the live services behind saved responses: MyGene.info, the UniProt REST API, the "
            "PeptideAtlas query endpoint; their responses as used are listed objects",
            "the licensed DeepTMHMM 1.0 package itself, which may not be redistributed; its file "
            "hashes are recorded above so a licensed copy can be checked",
            "vendor-side synthesis rules (Twist); only their report is saved",
        ],
    }
    META.write_text(json.dumps(meta, indent=2) + "\n")
    print(f"{len(rows)} objects ({meta['bytes'] / 1e9:.2f} GB) -> {MANIFEST.name}, {META.name}")
    return 0


def check() -> str:
    with MANIFEST.open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    items = [(r["bucket"], r["key"], r["sha256"]) for r in rows]
    missing, changed = verify_s3_objects(items, max_workers=12)
    if missing or changed:
        raise SystemExit(
            f"saved inputs: missing {missing[:5]} ({len(missing)}); changed {changed[:5]} ({len(changed)})"
        )
    total = sum(int(r["bytes"]) for r in rows)
    return f"all {len(rows)} saved pipeline inputs present on S3 with their recorded SHA-256 ({total / 1e9:.2f} GB)"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--from-cache", type=Path)
    ap.add_argument("--check", action="store_true")
    args = ap.parse_args(argv)
    if args.from_cache:
        return from_cache(args.from_cache)
    if args.check:
        print(check())
        return 0
    ap.print_help()
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
