"""Save the five rung builds at every commit an order of record was made from, so no order
depends on a git branch staying alive.

``synthesis_compliance.py`` reads each rung's ``constructs_l1_<rung>.fasta`` / ``.tsv`` with
``git show`` at a pinned commit. The final order's commit is on ``main``, but two earlier
orders and one refold proof were made from the draft integration branch
(``claude/l1-draft-integration-b6c53e``), whose commits are not on ``main``; deleting that
branch would destroy their inputs. This script packs the ten files at each commit into one
tarball per commit under ``BUNDLE_PREFIX`` (new dated keys, never overwritten) and records each
file's SHA-256 in ``rung_builds_saved.json``. ``check`` re-downloads every bundle and compares.

Requested by the coordinator (2026-11-14, independent audit item B1); implemented by the
Claude Code session l1-reviewer-package.

    uv run --locked python save_rung_builds.py save
    uv run --locked python save_rung_builds.py check
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
REL = HERE.relative_to(REPO).as_posix()
RECORD = HERE / "rung_builds_saved.json"
BUNDLE_PREFIX = "sprint4/s4e9-orf-design-fg-panel/pinned-inputs/rung-builds/"
RUNGS = ("nterm", "multipass", "type3", "internal", "cterm")
#: commit -> what was made from it
COMMITS = {
    "541c6c18d50c4da2c77a16122ba63aa22e9abf45": "order of record 2026-11-15-final-e (cut in its PR)",
    "4532a4d6374c1af156fa0ad7e69af208cb48bd7c": "order of record 2026-11-14-final-d (cut in its PR)",
    "903f913595d5ad34d8673dadf5e88d26740d5768": "order of record 2026-11-14-final-c (main)",
    "d6ab509371ddd9b84822afc171fb0cc2640746fe": "order of record 2026-11-14-final-b (main)",
    "cdb780686e700831271d68adbb9a3028788efa3d": "final order of record 2026-11-14-final (main)",
    "496d34825df6362b9abedb9c9a361e588c71e75d": (
        "superseded order 2026-11-13-one-arm-renamed (draft integration branch, not on main)"
    ),
    "5cfd91c02d2d7e889097e55b48d874277ff1306e": (
        "l1-tag-refold's ordered-construct proof, tag_refold/verify_ordered_constructs.py "
        "(draft integration branch, not on main)"
    ),
}


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def files_at(commit: str) -> dict[str, bytes]:
    out = {}
    for rung in RUNGS:
        for ext in ("fasta", "tsv"):
            name = f"constructs_l1_{rung}.{ext}"
            out[name] = subprocess.run(["git", "show", f"{commit}:{REL}/{name}"], cwd=REPO,
                                       check=True, capture_output=True).stdout  # fmt: skip
    return out


def pack(files: dict[str, bytes]) -> bytes:
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w", format=tarfile.PAX_FORMAT) as tar:
        for name in sorted(files):
            info = tarfile.TarInfo(name)
            info.size, info.mtime = len(files[name]), 0
            tar.addfile(info, io.BytesIO(files[name]))
    return buf.getvalue()


def save() -> int:
    from deliverome_analysis.s3 import processed

    record = {"bundles": {}}
    for commit, what in COMMITS.items():
        files = files_at(commit)
        data = pack(files)
        key = f"{BUNDLE_PREFIX}{commit[:8]}/rung_builds.tar"
        if not processed.exists(key):
            with tempfile.TemporaryDirectory() as tmp:
                path = Path(tmp) / "rung_builds.tar"
                path.write_bytes(data)
                processed.upload(path, key)
        record["bundles"][commit] = {
            "made": what,
            "key": f"processed://{key}",
            "sha256": sha256(data),
            "files": {n: sha256(b) for n, b in sorted(files.items())},
        }
    record["implemented_by"] = "Claude Code session l1-reviewer-package, 2026-11-14"
    RECORD.write_text(json.dumps(record, indent=2) + "\n")
    print(f"saved {len(COMMITS)} rung-build bundles")
    return 0


def check() -> str:
    from deliverome_analysis.s3 import processed

    record = json.loads(RECORD.read_text())["bundles"]
    for rec in record.values():
        data = Path(processed.download(rec["key"].removeprefix("processed://"))).read_bytes()
        if sha256(data) != rec["sha256"]:
            raise SystemExit(f"{rec['key']}: SHA-256 changed")
        with tarfile.open(fileobj=io.BytesIO(data)) as tar:
            got = {m.name: sha256(tar.extractfile(m).read()) for m in tar.getmembers()}
        if got != rec["files"]:
            raise SystemExit(f"{rec['key']}: member files differ from the record")
    # the final order's inputs: the files synthesis_compliance read must be the saved ones
    synth = json.loads((HERE / "synthesis_compliance_summary.json").read_text())
    ref = next(iter(set(synth["input_refs"].values())))
    if record.get(ref, {}).get("files") != synth["input_sha256"]:
        raise SystemExit(f"the order's rung builds at {ref[:8]} are not the saved ones")
    return f"{len(record)} rung-build bundles saved and unchanged; the order's inputs match"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("command", choices=("save", "check"))
    args = ap.parse_args(argv)
    if args.command == "save":
        return save()
    print(check())
    return 0


if __name__ == "__main__":
    sys.exit(main())
