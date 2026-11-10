"""Rerun the L1 library's decision steps from the frozen bundle and check them against the record.

For an outside reviewer: no Deliverome credentials and no network are needed. The bundle is
a folder laid out as ``s3/<external|processed>/<key>``, a copy of exactly the objects the
library pipeline reads, with a SHA-256 per file in ``bundle_manifest.tsv``. This script

1. checks every bundle file against the manifest, and refuses to run if any differs;
2. points the repo's S3 helper at the bundle (``DELIVEROME_S3_MIRROR``), and makes any
   attempt to reach real S3 fail, so a pass cannot have quietly used live data;
3. reruns each stage of ``library_pipeline.py`` -- the same verifiers the lab runs, unchanged --
   and records whether it reproduces its published record;
4. recomputes every number quoted in ``BRIEF.md`` (``brief_numbers.py``) and compares it with
   the committed ``brief_numbers.json``.

    uv run --locked python run_review.py --bundle /path/to/l1-reviewer-bundle
    uv run --locked python run_review.py --bundle ... --stage reference --stage tag_paths

Two stages, ``isoform_universe`` and ``ortholog_repick``, rerun scripts at the commits that
produced their records. They read them from this repository's git history or, in an exported
review repo, from the archives in ``historical/`` (see ``run_stage.py``); with neither, they
are reported as skipped, not passed.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path

HERE = Path(__file__).parent
S4E9 = HERE.parent
REPO = S4E9.parents[2]
MANIFEST = HERE / "bundle_manifest.tsv"
NEEDS_HISTORY = {"isoform_universe": "88786350", "ortholog_repick": "c4dbe752"}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def check_bundle(bundle: Path) -> int:
    """Every file the replay reads present with its recorded hash; returns the file count.

    Rows marked ``read_by_replay = no`` are the rest of the pinned releases, shipped as a
    separate supplement tarball; they are checked when present and not required.
    """
    with MANIFEST.open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    bad = []
    for r in rows:
        path = bundle / "s3" / r["bucket"] / r["key"]
        if not path.is_file():
            if r.get("read_by_replay", "yes") == "yes":
                bad.append(f"missing  {r['bucket']}://{r['key']}")
        elif sha256(path) != r["sha256"]:
            bad.append(f"changed  {r['bucket']}://{r['key']}")
    if bad:
        raise SystemExit("bundle does not match bundle_manifest.tsv:\n  " + "\n  ".join(bad))
    return len(rows)


def offline_env(bundle: Path, home: Path) -> dict[str, str]:
    """Environment for the replay: bundle reads only, real S3 unreachable, fixed hash seed."""
    env = dict(os.environ)
    env.update(
        {
            "DELIVEROME_S3_MIRROR": str(bundle / "s3"),
            "DELIVEROME_S3_MIRROR_MODE": "read",
            # If anything bypasses the mirror it meets a closed port and fails loudly.
            "S3_ENDPOINT_URL": "http://127.0.0.1:9",
            "S3_ACCESS_KEY_ID": "offline-replay",
            "S3_SECRET_ACCESS_KEY": "offline-replay",
            # A fresh home: no developer cache can stand in for a bundle file.
            "HOME": str(home),
            "DELIVEROME_ANALYSIS_S3_CACHE": str(home / "s3-cache"),
            "PYTHONHASHSEED": "0",
        }
    )
    return env


def has_history(commit: str) -> bool:
    if (REPO / "historical" / f"{commit}.tar").is_file():
        return True
    proc = subprocess.run(
        ["git", "-C", str(REPO), "cat-file", "-e", f"{commit}^{{commit}}"], capture_output=True
    )
    return proc.returncode == 0


def run_stage(stage: str, env: dict[str, str]) -> tuple[str, str, float]:
    if stage in NEEDS_HISTORY and not has_history(NEEDS_HISTORY[stage]):
        return "SKIP", f"needs git history (commit {NEEDS_HISTORY[stage]} not in this checkout)", 0
    t0 = time.time()
    proc = subprocess.run(
        [sys.executable, str(HERE / "run_stage.py"), stage],
        cwd=S4E9,
        env=env,
        capture_output=True,
        text=True,
    )
    lines = [ln for ln in proc.stdout.splitlines() if ln.startswith(("OK ", "FAIL "))]
    line = lines[-1] if lines else (proc.stderr.strip().splitlines() or ["no output"])[-1]
    status = "OK" if proc.returncode == 0 and line.startswith("OK ") else "FAIL"
    detail = line.split(None, 2)[-1] if line.startswith(("OK ", "FAIL ")) else line
    return status, detail, time.time() - t0


def check_numbers(env: dict[str, str], work: Path) -> tuple[str, str]:
    out = work / "brief_numbers.json"
    proc = subprocess.run(
        [sys.executable, str(HERE / "brief_numbers.py"), "--no-figures", "--out", str(out)],
        cwd=HERE,
        env=env,
        capture_output=True,
        text=True,
    )
    if proc.returncode:
        return "FAIL", proc.stderr.strip().splitlines()[-1]
    new = json.loads(out.read_text())
    old = json.loads((HERE / "brief_numbers.json").read_text())
    differ = sorted(k for k in old if new.get(k) != old[k])
    if differ:
        return "FAIL", f"sections differ from brief_numbers.json: {differ}"
    return "OK", f"all {len(old)} sections of brief_numbers.json reproduce exactly"


def main(argv=None) -> int:
    sys.path.insert(0, str(S4E9))
    import library_pipeline as lp

    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--bundle", type=Path, required=True, help="the extracted bundle folder")
    ap.add_argument("--stage", action="append", choices=[s.id for s in lp.STAGES])
    ap.add_argument("--report", type=Path, default=Path("review_report.json"))
    args = ap.parse_args(argv)

    n = check_bundle(args.bundle)
    print(f"bundle: {n} files, every SHA-256 matches bundle_manifest.tsv\n")
    home = Path(tempfile.mkdtemp(prefix="l1-review-home-"))
    env = offline_env(args.bundle.resolve(), home)

    results = []
    for i, s in enumerate(lp.STAGES, 1):
        if args.stage and s.id not in args.stage:
            continue
        print(f"{i:>2}. {s.id:<17} decides: {s.decides}")
        status, detail, secs = run_stage(s.id, env)
        print(f"    {status:<4} {detail}  ({secs:.0f} s)\n")
        results.append({"stage": s.id, "status": status, "detail": detail, "seconds": round(secs)})

    status, detail = check_numbers(env, home)
    print(f"brief numbers: {status} {detail}")
    results.append({"stage": "brief_numbers", "status": status, "detail": detail})

    args.report.write_text(json.dumps({"bundle_files": n, "results": results}, indent=1) + "\n")
    failed = [r["stage"] for r in results if r["status"] == "FAIL"]
    skipped = [r["stage"] for r in results if r["status"] == "SKIP"]
    print(
        f"\n{sum(r['status'] == 'OK' for r in results)} reproduced, {len(failed)} failed"
        f"{', ' + str(len(skipped)) + ' skipped' if skipped else ''}; report: {args.report}"
    )
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
