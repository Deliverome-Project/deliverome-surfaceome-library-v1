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
from collections.abc import Mapping, Sequence
from concurrent.futures import ThreadPoolExecutor
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, cast

PNG_KEY = "Deliverome-Provenance"
_PKG_DIR = Path(__file__).resolve().parent
PINNED_SNAPSHOT_DIR = "historical/git-snapshots"
PINNED_MANIFEST_NAME = "SHA256SUMS.json"


def sha256_file(path: Path, *, chunk_bytes: int = 8 << 20) -> str:
    """Streaming SHA-256 hex digest of a file on disk."""
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(chunk_bytes), b""):
            h.update(block)
    return h.hexdigest()


def verify_s3_objects(
    objects: Sequence[tuple[str, str, str]],
    *,
    max_workers: int = 12,
) -> tuple[list[str], list[str]]:
    """Verify ``(bucket_name, key, expected_sha256)`` tuples against S3 in parallel.

    Returns ``(missing_or_error_labels, changed_labels)`` in input order.
    """
    from deliverome_analysis.s3 import external, processed

    stores = {"external": external, "processed": processed}

    def _check_one(item: tuple[str, str, str]) -> tuple[str, str] | None:
        bucket, key, expected = item
        label = f"{bucket}://{key}"
        try:
            local = Path(stores[bucket].download(key))
        except Exception as exc:
            return ("missing", f"{label}: {exc}")
        if sha256_file(local) != expected:
            return ("changed", label)
        return None

    with ThreadPoolExecutor(max_workers=max(1, max_workers)) as pool:
        results = [res for res in pool.map(_check_one, objects) if res is not None]
    missing = [msg for kind, msg in results if kind == "missing"]
    changed = [msg for kind, msg in results if kind == "changed"]
    return missing, changed


def _repo_root(repo: Path | None) -> Path:
    if repo is not None:
        return repo
    for cand in (_PKG_DIR, *_PKG_DIR.parents):
        if (cand / PINNED_SNAPSHOT_DIR).is_dir() or (cand / "pyproject.toml").is_file():
            return cand
    return _PKG_DIR.parents[1]


def read_pinned_bytes(
    commit: str,
    rel_path: str,
    *,
    repo: Path | None = None,
    expected_sha256: str | None = None,
) -> bytes:
    """Read a file pinned at ``commit`` without requiring git history in exported archives.

    Resolution order:
    1. Materialized snapshot at ``<repo>/historical/git-snapshots/<commit>/<rel_path>``,
       verified against ``SHA256SUMS.json`` when present.
    2. Fallback to ``git -C <repo> show <commit>:<rel_path>`` inside a live git checkout.
    """
    root = _repo_root(repo)
    norm_rel = rel_path.lstrip("/")
    snap_root = root / PINNED_SNAPSHOT_DIR
    snap_file = snap_root / commit / norm_rel
    if snap_file.is_file():
        data = snap_file.read_bytes()
        manifest_path = snap_root / PINNED_MANIFEST_NAME
        if manifest_path.is_file():
            manifest = json.loads(manifest_path.read_text())
            key = f"{commit}/{norm_rel}"
            pinned = manifest.get(key)
            actual = hashlib.sha256(data).hexdigest()
            if pinned is None:
                raise ValueError(f"pinned snapshot missing from {PINNED_MANIFEST_NAME}: {key}")
            if actual != pinned:
                raise ValueError(
                    f"pinned snapshot checksum mismatch for {key}: expected {pinned}, got {actual}"
                )
    else:
        proc = subprocess.run(
            ["git", "-C", str(root), "show", f"{commit}:{norm_rel}"],
            capture_output=True,
            check=False,
        )
        if proc.returncode != 0:
            err = proc.stderr.decode("utf-8", errors="replace").strip()
            raise FileNotFoundError(
                f"pinned snapshot {commit}:{norm_rel} not found at {snap_file} "
                f"and git show failed: {err}"
            )
        data = proc.stdout

    if expected_sha256 is not None:
        actual = hashlib.sha256(data).hexdigest()
        if actual != expected_sha256:
            raise ValueError(
                f"checksum mismatch for {commit}:{norm_rel}: "
                f"expected {expected_sha256}, got {actual}"
            )
    return data


def read_pinned_text(
    commit: str,
    rel_path: str,
    *,
    repo: Path | None = None,
    expected_sha256: str | None = None,
) -> str:
    """UTF-8 text wrapper around :func:`read_pinned_bytes`."""
    return read_pinned_bytes(commit, rel_path, repo=repo, expected_sha256=expected_sha256).decode(
        "utf-8"
    )


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
        from PIL.PngImagePlugin import PngImageFile

        with Image.open(path) as im:
            raw = cast(PngImageFile, im).text.get(PNG_KEY)
        return json.loads(raw) if raw else {}
    if path.suffix.lower() == ".pdf":
        from pypdf import PdfReader

        info: Mapping[str, Any] = PdfReader(path).metadata or {}
        raw = info.get("/Keywords", "")
        return json.loads(raw) if raw.startswith("{") else {}
    raise ValueError(f"no provenance reader for {path.suffix or t}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(
            "usage: python -m deliverome_analysis.provenance <file.png|file.pdf|processed://key>"
        )
    stamp = read(sys.argv[1])
    print(json.dumps(stamp, indent=1) if stamp else "no provenance recorded")
