"""Current library membership, read from the one list that defines it.

Several tables here are snapshots of upstream sources -- the APPRIS release, the D1
attribute dump, the isoform universe -- and they keep every gene they were built with. A
gene removed from the library therefore lingers in them, and any script that inherits
membership from a snapshot instead of checking it will keep emitting that gene. EPPIN was
removed on 2026-10-11 and reappeared in four outputs for exactly this reason.

``library_exclusions.tsv`` beside the list records why anything was removed, so a
regeneration cannot silently put it back.
"""

from __future__ import annotations

import json
from pathlib import Path

LIBRARY_FILE = Path(__file__).parent.parent / "tag-site-architecture" / "library_members.json"


def library_members(path: str | Path = LIBRARY_FILE) -> set[str]:
    """Gene symbols currently in the library."""
    p = Path(path)
    return set(json.loads(p.read_text())) if p.exists() else set()
