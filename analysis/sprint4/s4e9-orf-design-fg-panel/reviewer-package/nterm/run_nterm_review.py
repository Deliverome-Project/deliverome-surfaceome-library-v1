"""Reproduce N-terminal tag placement for the L1 library from the frozen bundle, offline.

For an outside reviewer: no Deliverome credentials, no network, and no SignalP licence are
needed. The bundle is a folder laid out as ``s3/<external|processed>/<key>`` with a SHA-256
per file in ``nterm_bundle_manifest.tsv``. This script

1. checks every bundle file against the manifest and refuses to run if any differs;
2. points the repo's S3 helper at the bundle and makes real S3 unreachable (the parent
   package's ``offline_env``), so a pass cannot have used live data;
3. ``rungs``: reruns ``l1_tag_paths.py``, which assigns every L1 member its tag rung from
   DeepTMHMM v1 on its own designed protein, and compares the table byte for byte;
4. ``placement``: reruns ``nterm_placement.py`` (signal peptide, admissible interval,
   candidates, re-prediction verdicts and frame checks for every "N after SP" member) and
   compares its three outputs byte for byte with the committed ones, then prints what they
   say: agreement with the interval census, agreement with the sweep's verdicts, and any
   frame-check failure.

    uv run --locked python run_nterm_review.py --bundle /path/to/nterm-bundle

SignalP 6.0 and the disorder predictors are not rerun: their outputs are pinned in the
bundle, keyed by the SHA-256 of the sequence predicted. A licensed user can rerun SignalP
and compare (``signalp_rerun.py``).
"""

from __future__ import annotations

import argparse
import csv
import filecmp
import json
import subprocess
import sys
import tempfile
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
S4E9 = HERE.parents[1]
MANIFEST = HERE / "nterm_bundle_manifest.tsv"
OUTPUTS = ("nterm_members.tsv", "nterm_positions.tsv", "nterm_numbers.json")

sys.path.insert(0, str(HERE.parent))
from run_review import offline_env, sha256  # noqa: E402  (the parent package's runner)


def check_bundle(bundle: Path) -> int:
    with MANIFEST.open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    bad = []
    for r in rows:
        path = bundle / "s3" / r["bucket"] / r["key"]
        if not path.is_file():
            bad.append(f"missing  {r['bucket']}://{r['key']}")
        elif sha256(path) != r["sha256"]:
            bad.append(f"changed  {r['bucket']}://{r['key']}")
    if bad:
        raise SystemExit(f"bundle does not match {MANIFEST.name}:\n  " + "\n  ".join(bad))
    return len(rows)


def run(cmd: list[str], cwd: Path, env: dict) -> tuple[bool, str, float]:
    t0 = time.time()
    p = subprocess.run(cmd, cwd=cwd, env=env, capture_output=True, text=True)
    tail = (p.stderr.strip() or p.stdout.strip()).splitlines()[-1:] or [""]
    return p.returncode == 0, tail[0], time.time() - t0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--bundle", type=Path, required=True)
    ap.add_argument("--report", type=Path, default=Path("nterm_review_report.json"))
    args = ap.parse_args(argv)

    n = check_bundle(args.bundle)
    print(f"bundle: {n} files, every SHA-256 matches {MANIFEST.name}\n")
    home = Path(tempfile.mkdtemp(prefix="nterm-review-home-"))
    env = offline_env(args.bundle.resolve(), home)
    work = home / "out"
    work.mkdir()
    results = []

    ok, tail, secs = run(
        [sys.executable, "l1_tag_paths.py", "--out", str(work / "l1_tag_paths.tsv")], S4E9, env
    )
    same = ok and filecmp.cmp(work / "l1_tag_paths.tsv", S4E9 / "l1_tag_paths.tsv", shallow=False)
    status = "OK" if same else "FAIL"
    detail = (
        "l1_tag_paths.tsv reproduces byte for byte"
        if same
        else (tail if not ok else "l1_tag_paths.tsv differs")
    )
    print(f"1. rungs      {status}  {detail}  ({secs:.0f} s)")
    results.append({"stage": "rungs", "status": status, "detail": detail})

    ok, tail, secs = run(
        [sys.executable, str(HERE / "nterm_placement.py"), "--out", str(work)], HERE, env
    )
    differ = (
        [] if not ok else [f for f in OUTPUTS if not filecmp.cmp(work / f, HERE / f, shallow=False)]
    )
    status = "OK" if ok and not differ else "FAIL"
    detail = (
        tail
        if not ok
        else (f"differs: {differ}" if differ else "all three outputs reproduce byte for byte")
    )
    print(f"2. placement  {status}  {detail}  ({secs:.0f} s)")
    results.append({"stage": "placement", "status": status, "detail": detail})

    if ok:
        num = json.loads((work / "nterm_numbers.json").read_text())
        c, r, f = num["census"], num["reprediction"], num["frame"]
        print(
            f"\n   L1 'N after SP' members: {num['nafter_sp_members']}  classes: {num['classes']}"
            f"\n   interval census: {c['agree_all_fields']} of {c['members_in_census']} shared members agree on"
            f" class, cleavage, onsets, interval and candidates"
            f"\n   re-prediction: {r['positions_repredicted']} tagged designs found by SHA-256;"
            f" {r['verdict_agrees_with_sweep']} verdicts recomputed equal to the sweep's;"
            f" {r['members_not_covered']} members with candidates not re-predicted"
            f"\n   frame checks: {f['positions_checked'] - f['positions_failing']} of {f['positions_checked']}"
            f" candidate sites pass; failing members: {f['members_failing'] or 'none'}"
        )
    args.report.write_text(json.dumps({"bundle_files": n, "results": results}, indent=1) + "\n")
    failed = [x["stage"] for x in results if x["status"] != "OK"]
    print(f"\n{len(results) - len(failed)} reproduced, {len(failed)} failed; report: {args.report}")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
