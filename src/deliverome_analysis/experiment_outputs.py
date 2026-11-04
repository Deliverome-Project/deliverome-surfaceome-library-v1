"""Which files in an experiment folder are outputs, and whether S3 already holds them.

Shared by ``scripts/upload-experiment-outputs.sh`` (what to upload) and the pre-push gate
``scripts/check_experiment_outputs_uploaded.py`` (what must be in S3 before a push), so the
two cannot disagree about what an output is. They used to keep separate copies of the rules.

**What counts as an output.** Everything an experiment produces: tables, JSON, figures,
archives, logs, and any extension not listed below. Code, documentation and configuration are
not outputs; git is their durable home. The list is a denylist on purpose, so a new data
format is covered without anyone remembering to add it. Before this module, every file except
``README.md`` and ``experiment.py`` counted, which made each experiment's scripts and notes a
standing "missing from S3" backlog.

**What counts as uploaded.** An object with the same content anywhere under the experiment's
prefix, not only at the mirror path. Experiments publish dated or versioned snapshots
(``partner-audit-2026-09-10/``, Agentflow run folders), and a byte-identical copy there is just
as durable. Content is compared by ETag: the MD5 for single-part uploads, and for multipart
uploads the same digest boto3 computes with its default 8 MiB parts. A multipart object made
with some other part size cannot be verified, so it counts only at the mirror path.
"""

from __future__ import annotations

import hashlib
from collections.abc import Iterable, Mapping
from pathlib import Path

# Code, documentation and configuration: committed to git, never uploaded as outputs.
NOT_OUTPUT_SUFFIXES = frozenset(
    {
        ".py", ".r", ".sh", ".ipynb", ".sql",  # code
        ".md", ".rst", ".qmd", ".rmd",  # documentation
        ".yaml", ".yml", ".toml", ".cfg", ".ini",  # configuration and gate recipes
    }
)  # fmt: skip
# Folders holding inputs or caches rather than outputs.
SKIP_PARTS = frozenset({".venv", "__pycache__", ".ipynb_checkpoints", "data"})

# boto3's defaults (TransferConfig.multipart_threshold and multipart_chunksize), which
# ``deliverome_analysis.s3`` uploads with.
MULTIPART_THRESHOLD = 8 * 1024 * 1024
MULTIPART_CHUNKSIZE = 8 * 1024 * 1024


def is_output(path: Path) -> bool:
    """True when ``path`` (relative to its experiment folder) is an output that belongs in S3."""
    if path.name.startswith("."):
        return False
    if any(part in SKIP_PARTS for part in path.parts[:-1]):
        return False
    return path.suffix.lower() not in NOT_OUTPUT_SUFFIXES


def outputs(folder: Path) -> list[Path]:
    """Every output file under ``folder``, sorted."""
    return [
        path
        for path in sorted(folder.rglob("*"))
        if path.is_file() and is_output(path.relative_to(folder))
    ]


def local_etag(path: Path) -> str:
    """The ETag S3 gives ``path`` when ``deliverome_analysis.s3`` uploads it."""
    size = path.stat().st_size
    with path.open("rb") as handle:
        if size < MULTIPART_THRESHOLD:
            return hashlib.md5(handle.read()).hexdigest()
        digests = [
            hashlib.md5(chunk).digest()
            for chunk in iter(lambda: handle.read(MULTIPART_CHUNKSIZE), b"")
        ]
    return f"{hashlib.md5(b''.join(digests)).hexdigest()}-{len(digests)}"


def is_uploaded(path: Path, mirror_key: str, etags: Mapping[str, str]) -> bool:
    """True when S3 already holds ``path``'s content.

    ``etags`` maps every key under the experiment's prefix to its ETag
    (``Bucket.list_etags``). Any key with the same content counts. A multipart
    object whose ETag cannot be reproduced locally counts only at ``mirror_key``.
    """
    etag = local_etag(path)
    if etag in set(etags.values()):
        return True
    remote = etags.get(mirror_key)
    if remote is None or "-" not in remote or "-" not in etag:
        return False
    # Both multipart. The same number of parts means the same part size, so a different
    # ETag is different content. A different number of parts means another tool chose the
    # part size, and the content cannot be compared without downloading it.
    return remote.rsplit("-", 1)[1] != etag.rsplit("-", 1)[1]


def missing(folder: Path, s3_prefix: str, etags: Mapping[str, str]) -> list[Path]:
    """Outputs under ``folder`` whose content is not in S3 under ``s3_prefix``."""
    return [
        path
        for path in outputs(folder)
        if not is_uploaded(path, f"{s3_prefix}/{path.relative_to(folder).as_posix()}", etags)
    ]


def must_block(missing_outputs: Iterable[Path], changed: set[str], tracked: set[str]) -> list[Path]:
    """The missing outputs a push may not leave behind.

    An output blocks when this branch adds or changes it, or when git does not track it at
    all (a large or gitignored output whose only copy is on this machine). A tracked output
    the branch leaves alone is someone's older backlog: reported, not blocking.
    """
    return [
        path
        for path in missing_outputs
        if path.as_posix() in changed or path.as_posix() not in tracked
    ]
