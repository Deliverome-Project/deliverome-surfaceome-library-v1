"""Was a search run on the same raw file that is archived on ``raw://``?

A search records the SHA-256 and size of the raw file it read (deliverome-internal's runner,
``benchmark.json`` ``inputs[0]``). The archived copy on ``raw://`` carries its whole-file
SHA-256 as ``x-amz-meta-sha256`` (our uploader, since 2026-11-18) and its size. Comparing the
two catches a search of a partial or different file: on 2026-11-17/18 two S3E12 runs (AC1_01,
BC1_01) were searched on ~1 GB of a ~19 GB file and nothing downstream noticed.

Statuses (:data:`STATUSES`), from best to worst evidence:

``match``          sizes and SHA-256 agree.
``unverified``     sizes agree, but one side has no SHA-256 (older raw uploads carry only an
                   ETag, which for multipart uploads is not a content hash).
``mismatch``       sizes agree but SHA-256 differs.
``size_mismatch``  sizes differ (a partial copy, typically): the search did not see the file.
``raw_missing``    no raw object where the search's file should be archived.
``unknown``        could not check: no mapping from the search to a raw key, an ambiguous
                   mapping, a failed lookup (permissions, network), or no size/SHA-256 recorded
                   by the search.

S3's own ``ChecksumSHA256`` is deliberately not used: for a multipart upload it is a composite
checksum of the parts (``...-284``), not the SHA-256 of the file.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from typing import Any, Final, Literal, Protocol

Status = Literal["match", "mismatch", "size_mismatch", "unverified", "raw_missing", "unknown"]
STATUSES: Final[tuple[Status, ...]] = (
    "match",
    "unverified",
    "mismatch",
    "size_mismatch",
    "raw_missing",
    "unknown",
)
#: Statuses that mean the numbers may not describe the archived acquisition.
PROBLEM_STATUSES: Final[frozenset[Status]] = frozenset({"mismatch", "size_mismatch", "raw_missing"})
#: Statuses that cannot change for an immutable raw object, so a cached result is final.
FINAL_STATUSES: Final[frozenset[Status]] = frozenset(
    {"match", "mismatch", "size_mismatch", "unverified"}
)
SHA256_METADATA: Final = "sha256"  # x-amz-meta-sha256, written by the raw uploader
# S3 error codes that mean "this credential cannot look": stop asking for the rest of a build.
_DENIED: Final = frozenset({"403", "AccessDenied", "Forbidden", "InvalidAccessKeyId"})


@dataclass(frozen=True)
class RawObject:
    """What ``raw://`` says about one key.

    ``exists`` False means a HEAD returned 404 (``size``/``sha256``/``etag`` are then None).
    ``sha256`` is the uploader's whole-file ``x-amz-meta-sha256``; None for older uploads.
    """

    key: str
    exists: bool
    size: int | None = None
    sha256: str | None = None
    etag: str | None = None

    @classmethod
    def from_head(cls, key: str, head: Mapping[str, Any] | None) -> RawObject:
        """Build from an S3 HEAD response (None = 404)."""
        if head is None:
            return cls(key=key, exists=False)
        meta = head.get("Metadata") or {}
        sha = meta.get(SHA256_METADATA)
        size = head.get("ContentLength")
        etag = head.get("ETag")
        return cls(
            key=key,
            exists=True,
            size=int(size) if size is not None else None,
            sha256=str(sha).lower() if sha else None,
            etag=str(etag).strip('"') if etag else None,
        )


def classify(
    searched_sha256: str | None, searched_bytes: int | None, obj: RawObject | None
) -> Status:
    """Provenance of a search's input against the archived raw object.

    ``obj`` None means the raw object could not be determined or looked up (``unknown``).
    A size difference wins over a SHA-256 difference because it says more (a partial copy).
    """
    if obj is None:
        return "unknown"
    if not obj.exists:
        return "raw_missing"
    if searched_bytes is not None and obj.size is not None and searched_bytes != obj.size:
        return "size_mismatch"
    if searched_sha256 and obj.sha256:
        return "match" if searched_sha256.lower() == obj.sha256 else "mismatch"
    if searched_bytes is not None and obj.size is not None:
        return "unverified"  # sizes equal, no checksum on one side
    return "unknown"


def raw_prefix_candidates(folder: str) -> list[str]:
    """Raw prefixes that a published alphaDIA run folder can have come from, most likely first.

    ``scripts/publish_alphadia_runs.py`` publishes a batch whose raw files are in
    ``<raw prefix>`` to ``scratch/<scientist>/<raw prefix without a <scientist>/ segment>/
    <method>/<run>/``. Inverting it is ambiguous only in where the scientist segment was, so
    the candidates are the prefix as-is and with ``<scientist>/`` inserted before each segment
    but the batch folder. Folders outside ``scratch/`` (experiment folders) have no known
    mapping and give no candidates.
    """
    parts = folder.strip("/").split("/")
    if len(parts) < 5 or parts[0] != "scratch":
        return []
    scientist, batch = parts[1], parts[2:-2]  # drop <method>/<run>
    out = ["/".join(batch) + "/"]
    for i in range(len(batch)):  # the batch folder itself (batch[-1]) stays last
        out.append("/".join([*batch[:i], scientist, *batch[i:]]) + "/")
    return list(dict.fromkeys(out))


class RawBucket(Protocol):
    """The slice of :data:`deliverome_analysis.s3.raw` used here."""

    def head(self, key: str) -> Mapping[str, Any] | None: ...

    def list_keys(self, prefix: str = "", *, max_keys: int | None = None) -> list[str]: ...


class LookupFailed(Exception):
    """A raw:// lookup could not be answered (as opposed to: the object is absent)."""


def _error_code(exc: BaseException) -> str:
    response = getattr(exc, "response", None)
    if isinstance(response, Mapping):
        return str((response.get("Error") or {}).get("Code", ""))
    return ""


class RawLookup:
    """Resolves search folders to raw objects with as few S3 calls as possible.

    * ``cache``: raw objects known to exist from earlier builds; they are immutable, so they
      are never looked up again. Absent objects are always re-checked (an upload may be late).
    * Each batch's raw prefix is resolved once per lookup object (one LIST per candidate).
    * After the first permission error every further lookup is answered ``None`` without
      calling S3, so a credential without raw:// access costs one failed call per build.
    """

    def __init__(self, bucket: RawBucket | None, cache: Iterable[RawObject] = ()) -> None:
        self.bucket = bucket
        self.objects: dict[str, RawObject] = {o.key: o for o in cache if o.exists}
        self._prefixes: dict[str, str | None] = {}
        self.denied: str | None = None  # the first permission error, if any
        self.calls = 0

    def _guard(self, exc: Exception) -> LookupFailed:
        if _error_code(exc) in _DENIED or type(exc).__name__ in {
            "NoCredentialsError",
            "PartialCredentialsError",
        }:
            self.denied = f"{type(exc).__name__}: {exc}"
        return LookupFailed(f"{type(exc).__name__}: {exc}")

    def _check(self) -> RawBucket:
        if self.bucket is None:
            raise LookupFailed("no raw bucket")
        if self.denied:
            raise LookupFailed(f"raw:// not readable ({self.denied})")
        return self.bucket

    def prefix(self, folder: str) -> str | None:
        """The raw prefix of a published run folder; None when no candidate exists on raw://.

        Raises :class:`LookupFailed` when there is no mapping, the mapping is ambiguous (two
        candidates exist), or S3 could not be asked.
        """
        candidates = raw_prefix_candidates(folder)
        if not candidates:
            raise LookupFailed(f"no raw:// mapping for {folder!r}")
        batch_key = candidates[0]
        if batch_key in self._prefixes:
            return self._prefixes[batch_key]
        bucket = self._check()
        found = []
        for c in candidates:
            try:
                self.calls += 1
                if bucket.list_keys(c, max_keys=1):
                    found.append(c)
            except Exception as exc:
                raise self._guard(exc) from exc
        if len(found) > 1:
            raise LookupFailed(f"ambiguous raw:// prefix for {folder!r}: {found}")
        self._prefixes[batch_key] = found[0] if found else None
        return self._prefixes[batch_key]

    def object(self, key: str) -> RawObject:
        """HEAD one raw key (cached when it exists). Raises :class:`LookupFailed` on errors."""
        if key in self.objects:
            return self.objects[key]
        bucket = self._check()
        try:
            self.calls += 1
            obj = RawObject.from_head(key, bucket.head(key))
        except Exception as exc:
            raise self._guard(exc) from exc
        if obj.exists:
            self.objects[key] = obj
        return obj

    def lookup(self, folder: str, raw_file: str) -> tuple[str | None, RawObject | None]:
        """(raw key, raw object) for a run folder and its raw file name.

        The object is None when the lookup failed (``unknown``). A missing raw folder gives
        the most likely key with an absent object (``raw_missing``).
        """
        try:
            prefix = self.prefix(folder)
            if prefix is None:
                key = raw_prefix_candidates(folder)[0] + raw_file
                return key, RawObject(key=key, exists=False)
            key = prefix + raw_file
            return key, self.object(key)
        except LookupFailed:
            return None, None


@dataclass(frozen=True)
class SearchRaw:
    """What the search index knows about one search's raw input and its archived copy."""

    run_uid: str | None
    raw_sha256: str | None
    raw_s3_key: str | None
    raw_s3_bytes: int | None
    raw_s3_sha256: str | None
    raw_provenance: str | None


def archived_object(
    searched_sha256: str | None, run_uid: str | None, searches: Sequence[SearchRaw]
) -> RawObject | None:
    """The archived raw object behind a QC run, from the search index; None if not unique.

    A QC run is matched to the searches of the same searched file (its search's SHA-256),
    or failing that to the searches of the same acquisition (``run_uid``). The matches must
    all point at one raw key that the search index has looked up (status not ``unknown``).
    """
    known = [s for s in searches if s.raw_s3_key and s.raw_provenance not in (None, "unknown")]
    for hits in (
        [s for s in known if searched_sha256 and s.raw_sha256 == searched_sha256],
        [s for s in known if run_uid and s.run_uid == run_uid],
    ):
        if len({s.raw_s3_key for s in hits}) == 1:
            s = hits[0]
            assert s.raw_s3_key is not None
            if s.raw_provenance == "raw_missing":
                return RawObject(key=s.raw_s3_key, exists=False)
            return RawObject(
                key=s.raw_s3_key, exists=True, size=s.raw_s3_bytes, sha256=s.raw_s3_sha256
            )
    return None
