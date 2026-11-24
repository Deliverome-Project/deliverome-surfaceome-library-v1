"""Read ``tag_refold_queue.tsv`` with its sequences joined back in.

The queue keeps one row per fold but not the sequences: native sequences live once each in
``tag_refold_natives.tsv`` and tagged sequences in ``tag_refold_tagged.tsv``, both keyed by
SHA-256, so every committed table stays under the repository's 5 MB limit (the queue passed it
on 2026-11-15 at about 5,400 folds). A queue written before the split still carries the
columns and is read as is.

Implemented by Claude Code session l1-tag-refold.
"""

from __future__ import annotations

import csv
from pathlib import Path


def _sequences(path: Path, key: str, col: str) -> dict[str, str]:
    if not path.exists():
        return {}
    with path.open() as fh:
        return {r[key]: r[col] for r in csv.DictReader(fh, delimiter="\t")}


def read_queue(path: Path) -> list[dict[str, str]]:
    """Queue rows, each with ``native_seq`` and ``tagged_seq``; stops if a sequence is missing."""
    path = Path(path)
    with path.open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    natives = _sequences(path.with_name("tag_refold_natives.tsv"), "native_sha256", "native_seq")
    tagged = _sequences(path.with_name("tag_refold_tagged.tsv"), "tagged_sha256", "tagged_seq")
    for r in rows:
        if not r.get("native_seq"):
            r["native_seq"] = natives[r["native_sha256"]]
        if not r.get("tagged_seq"):
            r["tagged_seq"] = tagged[r["tagged_sha256"]]
    return rows
