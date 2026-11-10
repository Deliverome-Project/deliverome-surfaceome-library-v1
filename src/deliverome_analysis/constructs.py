"""Construct (plasmid) maps from ``processed://constructs/``, looked up by name.

Maps are filed by status (``ordered/``, ``designed/``, ``reference/<source>/``), and a
map moves between those folders when its status changes. An analysis should not
hard-code which folder a map is in today, so it asks by name instead::

    from deliverome_analysis import constructs

    p008 = constructs.fetch("p008")                # constructs/ordered/p008.gb
    cmd = constructs.fetch("addgene/235285")       # constructs/reference/addgene/235285.gbk

A name matches a key when it equals the key's trailing path components with the
extension removed, ignoring case: ``"p008"`` matches ``constructs/ordered/p008.gb`` but
not ``constructs/ordered/p0080.gb``. Exactly one key must match.

S3 is the record. A committed copy of a map is allowed for review or offline runs
(``fallback=``), and is used only when S3 cannot be reached at all; a map that S3
answers for but does not hold is an error, not a reason to read the local copy, so a
missing upload fails loudly instead of silently diverging.
"""

from __future__ import annotations

import sys
from pathlib import Path

PREFIX = "constructs/"
EXTENSIONS = (".gb", ".gbk", ".genbank")


class ConstructNotFound(LookupError):
    """No key, or more than one key, under ``constructs/`` matches the name."""


def _stem(key: str) -> str:
    lowered = key.lower()
    for ext in EXTENSIONS:
        if lowered.endswith(ext):
            return lowered[: -len(ext)]
    return ""


def find_key(name: str, keys: list[str]) -> str:
    """The one key in ``keys`` that ``name`` identifies (see the module docstring)."""
    wanted = name.strip("/").lower()
    matches = [k for k in keys if (s := _stem(k)) and (s == wanted or s.endswith("/" + wanted))]
    if len(matches) == 1:
        return matches[0]
    if not matches:
        raise ConstructNotFound(
            f"no map named {name!r} under processed://{PREFIX}. Upload it first, filed by "
            f"status: {PREFIX}ordered/, {PREFIX}designed/ or {PREFIX}reference/<source>/ "
            f"(see docs/s3-data-access.md)."
        )
    raise ConstructNotFound(f"{name!r} is ambiguous under processed://{PREFIX}: {matches}")


def key_for(name: str, *, bucket=None) -> str:
    """The ``processed://`` key of the map called ``name``."""
    if bucket is None:
        from deliverome_analysis import s3

        bucket = s3.processed
    return find_key(name, bucket.list_keys(prefix=PREFIX))


def map_key(name: str, *, new_status: str = "designed", bucket=None) -> str:
    """The key to read or write the map called ``name``, wherever it is filed now.

    An existing map is found in its current status folder, so a script that edits a map
    and writes it back writes it to the same place even after the map moved from
    ``designed/`` to ``ordered/``. A name with no map yet gets
    ``constructs/<new_status>/<name>.gb``, the home of a new design. An ambiguous name
    is still an error.
    """
    if bucket is None:
        from deliverome_analysis import s3

        bucket = s3.processed
    try:
        return find_key(name, bucket.list_keys(prefix=PREFIX))
    except ConstructNotFound as error:
        if "ambiguous" in str(error):
            raise
        return f"{PREFIX}{new_status}/{name.strip('/')}.gb"


def fetch(name: str, *, fallback: Path | None = None, bucket=None) -> Path:
    """Download the map called ``name`` and return its local path.

    ``fallback`` is a committed review copy to use when S3 is unreachable (no
    credentials, no network). It is never used when S3 answers.
    """
    try:
        if bucket is None:
            from deliverome_analysis import s3

            bucket = s3.processed
        keys = bucket.list_keys(prefix=PREFIX)
    except Exception as error:
        if fallback is None:
            raise
        print(
            f"note: S3 unreachable ({error.__class__.__name__}); using the committed "
            f"review copy {fallback.name} for {name!r}",
            file=sys.stderr,
        )
        return fallback
    return bucket.download(find_key(name, keys))
