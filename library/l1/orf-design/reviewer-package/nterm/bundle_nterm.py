"""Build the frozen input bundle that ``run_nterm_review.py`` replays. Lab-side: needs S3.

Same approach as the parent package's ``bundle.py``, on a much smaller input set. ``collect``
reruns the two steps the reviewer runner replays (``l1_tag_paths.py`` and
``nterm_placement.py``) with the S3 helper in record mode and a fresh home directory, so the
bundle holds exactly the objects they read. ``manifest`` hashes each file and labels where it
comes from, which decides what may be shared. ``pack`` writes a tarball.

    uv run --locked python bundle_nterm.py collect --out ~/nterm-bundle
    uv run --locked python bundle_nterm.py manifest --bundle ~/nterm-bundle
    uv run --locked python bundle_nterm.py pack --bundle ~/nterm-bundle --out nterm-reviewer-bundle.tar

Every file here is also needed by, or identical in kind to, the parent bundle's; a reviewer
who already holds the parent bundle can unpack this one over it.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
S4E9 = HERE.parents[1]
MANIFEST = HERE / "nterm_bundle_manifest.tsv"
P = r"processed/sprint4/s4e9-orf-design-fg-panel/"

#: Where each bundled object comes from (first match wins). The class decides sharing, with
#: the parent package's meanings; ``licensed_tool_output`` is new here and is the class to
#: ask about before outside reviewers receive it.
SOURCES: list[tuple[str, str, str]] = [
    (
        P + r"reviewer-package/nterm/[^/]+/predictor_rows",
        "licensed_tool_output",
        "SignalP 6.0 (DTU licence) and NetSurfP-3.0 / metapredict / AlphaFold-disorder rows "
        "from the accessible-surfaceome public D1 database, frozen 2026-11-04",
    ),
    (
        P + r"reviewer-package/nterm/[^/]+/canary_",
        "licensed_tool_output",
        "SignalP 6.0 panel re-prediction, canary_rp_1792371023 (Modal volume copy)",
    ),
    (
        P + r"reprediction/",
        "licensed_tool_output",
        "SignalP 6.0 / metapredict / NetSurfP-3.0 re-prediction runs (Deliverome Modal)",
    ),
    (
        P + r"triage/l1_interval_census",
        "deliverome_record",
        "L1 interval census, 2026-10-18 (PR 191)",
    ),
    (
        r"external/accessible-surfaceome/topology/",
        "deliverome_curation",
        "DeepTMHMM 1.0.24 predictions run by Deliverome (topology snapshot)",
    ),
    (P + r"library-candidates/", "deliverome_record", "published L1 library record"),
]


def classify(rel: str) -> tuple[str, str]:
    for pattern, cls, origin in SOURCES:
        if re.match(pattern, rel):
            return cls, origin
    raise SystemExit(f"no SOURCES rule for {rel}; add one before sharing the bundle")


def sha256(path: Path) -> str:
    import hashlib

    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def collect(out: Path) -> int:
    mirror = out / "s3"
    mirror.mkdir(parents=True, exist_ok=True)
    env = {
        **os.environ,
        "DELIVEROME_S3_MIRROR": str(mirror),
        "DELIVEROME_S3_MIRROR_MODE": "record",
        "HOME": tempfile.mkdtemp(prefix="nterm-bundle-home-"),
        "PYTHONHASHSEED": "0",
    }
    for var, name in (
        ("AWS_CONFIG_FILE", "config"),
        ("AWS_SHARED_CREDENTIALS_FILE", "credentials"),
    ):
        real = Path.home() / ".aws" / name
        if var not in env and real.is_file():
            env[var] = str(real)
    work = Path(tempfile.mkdtemp(prefix="nterm-collect-"))
    steps = [
        ([sys.executable, "l1_tag_paths.py", "--out", str(work / "l1_tag_paths.tsv")], S4E9),
        ([sys.executable, str(HERE / "nterm_placement.py"), "--out", str(work)], HERE),
    ]
    for cmd, cwd in steps:
        rc = subprocess.run(cmd, cwd=cwd, env=env, capture_output=True, text=True)
        if rc.returncode:
            print(rc.stdout[-2000:], rc.stderr[-2000:])
            return rc.returncode
    print(f"collected into {mirror}")
    return 0


def manifest(bundle: Path) -> int:
    root = bundle / "s3"
    rows = []
    for path in sorted(p for p in root.rglob("*") if p.is_file()):
        rel = path.relative_to(root).as_posix()
        bucket, _, key = rel.partition("/")
        cls, origin = classify(rel)
        rows.append(
            {
                "bucket": bucket,
                "key": key,
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
                "class": cls,
                "origin": origin,
            }
        )
    if not rows:
        raise SystemExit(f"{root} holds no files; refusing to overwrite {MANIFEST.name}")
    with MANIFEST.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    for r in rows:
        print(f"{r['class']:<22} {int(r['bytes']) / 1e6:>7.1f} MB  {r['bucket']}://{r['key']}")
    print(f"wrote {MANIFEST.name}: {len(rows)} files")
    return 0


def pack(bundle: Path, out: Path, exclude: list[str]) -> int:
    with MANIFEST.open() as fh:
        rows = [r for r in csv.DictReader(fh, delimiter="\t") if r["class"] not in exclude]
    with tarfile.open(out, "w") as tar:
        tar.add(MANIFEST, arcname=MANIFEST.name)
        for r in rows:
            tar.add(bundle / "s3" / r["bucket"] / r["key"], arcname=f"s3/{r['bucket']}/{r['key']}")
    print(f"wrote {out}: {len(rows)} files" + (f" (excluded {exclude})" if exclude else ""))
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("collect").add_argument("--out", type=Path, required=True)
    sub.add_parser("manifest").add_argument("--bundle", type=Path, required=True)
    p = sub.add_parser("pack")
    p.add_argument("--bundle", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--exclude-class", action="append", default=[])
    a = ap.parse_args(argv)
    if a.cmd == "collect":
        return collect(a.out.expanduser())
    if a.cmd == "manifest":
        return manifest(a.bundle.expanduser())
    return pack(a.bundle.expanduser(), a.out, a.exclude_class)


if __name__ == "__main__":
    raise SystemExit(main())
