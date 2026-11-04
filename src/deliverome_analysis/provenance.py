"""Stamp outputs with the code and package versions that made them.

``save_figure`` writes this into every PNG/PDF and ``processed.upload`` / ``external.upload``
attach it to every S3 object, so any figure or table answers "which commit, and which
locked package versions, produced me?":

=====================  ======================================================================
``git-commit``         the commit checked out when the output was made
``git-dirty``          ``true`` if there were uncommitted changes (then the commit alone does
                       not reproduce the output)
``script``             the script that ran, relative to the repo root
``uv-lock-sha256``     first 16 hex digits of ``uv.lock``'s SHA-256 (Python package versions)
``renv-lock-sha256``   same for ``renv.lock`` (R package versions)
``created``            UTC timestamp
=====================  ======================================================================

To reproduce: ``git checkout <git-commit>``, ``uv sync --frozen`` (and
``scripts/install-r-deps.sh`` for R). The lock hashes confirm the environment matches.

Read it back::

    uv run python -m deliverome_analysis.provenance figures/s1e3_main.png
    uv run python -m deliverome_analysis.provenance processed://sprint1/s1e3-pilot-screen/hits.tsv
"""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from datetime import UTC, datetime
from pathlib import Path

PNG_KEY = "Deliverome-Provenance"
_PKG_DIR = Path(__file__).resolve().parent


def _git(*args: str, cwd: Path) -> str | None:
    try:
        proc = subprocess.run(
            ["git", *args], cwd=cwd, capture_output=True, text=True, timeout=10, check=False
        )
    except OSError, subprocess.TimeoutExpired:
        return None
    return proc.stdout.strip() if proc.returncode == 0 else None


def _sha16(path: Path) -> str | None:
    if not path.is_file():
        return None
    return hashlib.sha256(path.read_bytes()).hexdigest()[:16]


def provenance() -> dict[str, str]:
    """Current code/environment stamp; fields that cannot be determined are left out."""
    root_s = _git("rev-parse", "--show-toplevel", cwd=_PKG_DIR)
    out: dict[str, str] = {"created": datetime.now(tz=UTC).strftime("%Y-%m-%dT%H:%M:%SZ")}
    if not root_s:
        return out  # installed outside a git checkout
    root = Path(root_s)
    commit = _git("rev-parse", "HEAD", cwd=root)
    if commit:
        out["git-commit"] = commit
        status = _git("status", "--porcelain", "--untracked-files=no", cwd=root)
        out["git-dirty"] = "true" if status else "false"
    script = Path(sys.argv[0]).resolve() if sys.argv and sys.argv[0] else None
    if script and script.is_file() and script.is_relative_to(root):
        out["script"] = script.relative_to(root).as_posix()
    for name in ("uv.lock", "renv.lock"):
        digest = _sha16(root / name)
        if digest:
            out[f"{name.replace('.', '-')}-sha256"] = digest
    return out


def warn_if_dirty(stamp: dict[str, str], what: str) -> None:
    if stamp.get("git-dirty") == "true":
        print(
            f"  note: {what} was made with uncommitted changes; commit first if it must be "
            f"reproducible from git-commit {stamp['git-commit'][:12]}."
        )


def read(target: str | Path) -> dict[str, str]:
    """Provenance of a local PNG/PDF or of an S3 object (``processed://key``/``external://key``)."""
    t = str(target)
    if "://" in t:
        from deliverome_analysis.s3 import client, external, processed

        scheme, key = t.split("://", 1)
        bucket = {"processed": processed, "external": external}[scheme]._name
        return dict(client().head_object(Bucket=bucket, Key=key).get("Metadata", {}))
    path = Path(t)
    if path.suffix.lower() == ".png":
        from PIL import Image

        with Image.open(path) as im:
            raw = im.text.get(PNG_KEY)
        return json.loads(raw) if raw else {}
    if path.suffix.lower() == ".pdf":
        from pypdf import PdfReader

        raw = (PdfReader(path).metadata or {}).get("/Keywords", "")
        return json.loads(raw) if raw.startswith("{") else {}
    raise ValueError(f"no provenance reader for {path.suffix or t}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(
            "usage: python -m deliverome_analysis.provenance <file.png|file.pdf|processed://key>"
        )
    stamp = read(sys.argv[1])
    print(json.dumps(stamp, indent=1) if stamp else "no provenance recorded")
