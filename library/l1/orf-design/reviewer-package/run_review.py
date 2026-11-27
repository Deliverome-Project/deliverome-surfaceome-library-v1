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
import json
import os
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).parent
S4E9 = HERE.parent
REPO = S4E9.parents[2]
MANIFEST = HERE / "bundle_manifest.tsv"
NEEDS_HISTORY = {"isoform_universe": "88786350", "ortholog_repick": "c4dbe752"}


def sha256(path: Path) -> str:
    from deliverome_analysis.core.provenance import sha256_file

    return sha256_file(Path(path))


def verify_bundle(bundle: Path) -> tuple[int, int]:
    """Verify bundle files against ``bundle_manifest.tsv``.

    Returns ``(checked_count, skipped_optional_count)``. Rows marked ``read_by_replay = no`` are
    the rest of the pinned releases, shipped as a separate supplement tarball; they are checked
    when present on disk and counted as skipped (not checked) when absent.
    """
    with MANIFEST.open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    def _check_one(r: dict[str, str]) -> tuple[str, str | None]:
        path = bundle / "s3" / r["bucket"] / r["key"]
        if not path.is_file():
            if r.get("read_by_replay", "yes") == "yes":
                return "bad", f"missing  {r['bucket']}://{r['key']}"
            return "skipped", None
        if sha256(path) != r["sha256"]:
            return "bad", f"changed  {r['bucket']}://{r['key']}"
        return "ok", None

    with ThreadPoolExecutor(max_workers=min(8, os.cpu_count() or 4)) as pool:
        outcomes = list(pool.map(_check_one, rows))
    bad = [err for status, err in outcomes if status == "bad" and err is not None]
    if bad:
        raise SystemExit("bundle does not match bundle_manifest.tsv:\n  " + "\n  ".join(bad))
    checked = sum(1 for status, _ in outcomes if status == "ok")
    skipped = sum(1 for status, _ in outcomes if status == "skipped")
    return checked, skipped


def check_bundle(bundle: Path) -> int:
    """Every present bundle file matches its recorded SHA-256; returns count actually verified."""
    checked, _ = verify_bundle(bundle)
    return checked


def offline_env(bundle: Path, home: Path) -> dict[str, str]:
    """Environment for the replay: bundle reads only, S3 and HTTP(S) unreachable, fixed hash seed.

    Sets ``S3_ENDPOINT_URL`` and standard HTTP(S) proxy variables to ``http://127.0.0.1:9``
    (discard/closed port) so accidental live S3 or HTTP requests via ``boto3``, ``requests``,
    ``urllib`` or ``httpx`` fail immediately. For full kernel-level network isolation (blocking
    raw sockets / DNS as well), run inside the container with ``docker run --network=none``.
    """
    env = dict(os.environ)
    env.update(
        {
            "DELIVEROME_S3_MIRROR": str(bundle / "s3"),
            "DELIVEROME_S3_MIRROR_MODE": "read",
            # If anything bypasses the mirror or makes an HTTP(S) request, it meets a closed
            # local port and fails loudly.
            "S3_ENDPOINT_URL": "http://127.0.0.1:9",
            "S3_ACCESS_KEY_ID": "offline-replay",
            "S3_SECRET_ACCESS_KEY": "offline-replay",
            "HTTP_PROXY": "http://127.0.0.1:9",
            "HTTPS_PROXY": "http://127.0.0.1:9",
            "ALL_PROXY": "http://127.0.0.1:9",
            "NO_PROXY": "",
            "http_proxy": "http://127.0.0.1:9",
            "https_proxy": "http://127.0.0.1:9",
            "all_proxy": "http://127.0.0.1:9",
            "no_proxy": "",
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
    try:
        proc = subprocess.run(
            ["git", "-C", str(REPO), "cat-file", "-e", f"{commit}^{{commit}}"], capture_output=True
        )
    except FileNotFoundError:
        return False
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
    ap.add_argument(
        "-j",
        "--jobs",
        type=int,
        default=min(4, os.cpu_count() or 2),
        help="max stages to run in parallel (default: min(4, cpu_count))",
    )
    args = ap.parse_args(argv)

    checked, skipped_optional = verify_bundle(args.bundle)
    note = (
        f" ({skipped_optional} optional supplement files not unpacked)" if skipped_optional else ""
    )
    print(f"bundle: {checked} files verified against bundle_manifest.tsv{note}\n")
    home = Path(tempfile.mkdtemp(prefix="l1-review-home-"))
    env = offline_env(args.bundle.resolve(), home)

    selected = [(i, s) for i, s in enumerate(lp.STAGES, 1) if not args.stage or s.id in args.stage]
    jobs = max(1, min(args.jobs, len(selected)))
    if jobs == 1:
        stage_runs = [(i, s, *run_stage(s.id, env)) for i, s in selected]
    else:
        with ThreadPoolExecutor(max_workers=jobs) as pool:
            stage_runs = list(
                pool.map(lambda item: (item[0], item[1], *run_stage(item[1].id, env)), selected)
            )

    results = []
    for i, s, status, detail, secs in stage_runs:
        print(f"{i:>2}. {s.id:<17} decides: {s.decides}")
        print(f"    {status:<4} {detail}  ({secs:.0f} s)\n")
        results.append({"stage": s.id, "status": status, "detail": detail, "seconds": round(secs)})

    status, detail = check_numbers(env, home)
    print(f"brief numbers: {status} {detail}")
    results.append({"stage": "brief_numbers", "status": status, "detail": detail})

    args.report.write_text(
        json.dumps(
            {
                "bundle_files_verified": checked,
                "bundle_files_optional_absent": skipped_optional,
                "results": results,
            },
            indent=1,
        )
        + "\n"
    )
    failed = [r["stage"] for r in results if r["status"] == "FAIL"]
    skipped = [r["stage"] for r in results if r["status"] == "SKIP"]
    print(
        f"\n{sum(r['status'] == 'OK' for r in results)} reproduced, {len(failed)} failed"
        f"{', ' + str(len(skipped)) + ' skipped' if skipped else ''}; report: {args.report}"
    )
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
