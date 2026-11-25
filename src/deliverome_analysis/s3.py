"""Thin AWS S3 client for scientists.

Three buckets — ``raw``, ``processed``, ``external`` — each accessed
through a module-level helper (``s3.raw``, ``s3.processed``,
``s3.external``). The bucket *is* the security boundary; we don't need
prefix scoping or a bucket policy.

The default endpoint is AWS S3 (the boto3 default for the configured
region). Set ``S3_ENDPOINT_URL`` in ``.env`` to point at any S3-
compatible provider (Cloudflare R2, DigitalOcean Spaces, MinIO, Wasabi,
…) — the rest of this module is provider-agnostic.

==============================================================================
Bucket layout (raw mirrors processed exactly)
==============================================================================

::

    raw         sprint{N}/s{N}e{M}-{slug}/...          read-only, mirrors processed
    processed   sprint{N}/s{N}e{M}-{slug}/...          read + write + delete
    external    {source}/{version}/...                  read + write, NO delete (append-only)
                experiments/sprint{N}/s{N}e{M}[-{slug}]/  offsite / off-instrument raw data

==============================================================================
IAM model (AWS)
==============================================================================

One IAM user per scientist, one access-key pair (``S3_ACCESS_KEY_ID`` /
``S3_SECRET_ACCESS_KEY``). The IAM policy is least-privilege per bucket:

* **raw** — ``s3:GetObject`` + ``s3:ListBucket`` only. The immutable
  instrument-output bucket; the single key cannot write to it, and a bucket
  policy reserves deletion for the account root.
* **external** — ``s3:GetObject``, ``s3:ListBucket``, ``s3:PutObject``, and
  the multipart actions, with an **explicit Deny** on ``s3:DeleteObject`` /
  ``s3:DeleteObjectVersion`` / ``s3:PutBucketVersioning``. Scientists can
  publish new reference snapshots but can never destroy one. The bucket is
  **versioned**, so an overwrite of an existing key supersedes rather than
  destroys — the prior version survives and cannot be purged with a
  scientist key. Versioning is what makes "write but not delete" mean
  something here; the Deny on ``PutBucketVersioning`` stops anyone turning
  it back off.
* **processed** — full ``s3:GetObject``, ``s3:PutObject``,
  ``s3:DeleteObject``, ``s3:ListBucket``. Read and write on the mutable
  bucket. ``s3:DeleteObjectVersion`` is *not* granted and the bucket is
  versioned, so a scientist delete leaves a recoverable delete marker rather
  than destroying the object.

The same ``client()`` object backs all three accessors; per-bucket
enforcement is done server-side by the IAM policy, not by different keys.

See ``docs/s3-data-access.md`` for the policy JSON and one-time IAM setup.

``raw`` doesn't expose an ``upload`` method, so ``s3.raw.upload(...)`` is an
``AttributeError`` at the Python layer. ``s3.processed.upload()`` validates
the key against the experiment-folder layout (see
:data:`_PROCESSED_KEY_PATTERN`); ``s3.external.upload()`` validates the
``{source}/{version}/...`` layout instead (see
:data:`_EXTERNAL_KEY_PATTERN`). Neither read bucket exposes a delete, and
nothing in this module deletes from ``external`` — the IAM Deny is the
server-side backstop.

==============================================================================
Quick start
==============================================================================

::

    from deliverome_analysis.s3 import raw, processed, external

    # Browse
    processed.list_keys(prefix="sprint1/")
    raw.list_keys(prefix="sprint1/")

    # Read
    df = processed.read_tsv("sprint1/s1e3-pilot-screen/hits.tsv")
    df = processed.read_parquet("sprint1/s1e3-pilot-screen/expression.parquet")
    fastq_path = raw.download("sprint1/s1e3-pilot-screen/sample01.fastq.gz")
    df = external.read_tsv("uniprot/2026_03/human.tsv")

    # Write (processed only, with key validation)
    processed.upload("/tmp/hits.tsv", "sprint1/s1e3-pilot-screen/hits.tsv")

All reads are cached at ``~/.cache/deliverome-analysis/s3/{bucket}/{key}.{etag}``
so repeats in a notebook are near-instant.

S3 bucket names are GLOBALLY unique across all of AWS, so the conceptual
names (``raw``, ``processed``, ``external``) almost certainly collide
with someone else's. Set ``S3_BUCKET_RAW`` / ``S3_BUCKET_PROCESSED`` /
``S3_BUCKET_EXTERNAL`` in ``.env`` to the actual bucket names you
provisioned (e.g. ``deliverome-raw``).
"""

from __future__ import annotations

import io
import os
import re
import shutil
import tempfile
from collections.abc import Callable
from functools import lru_cache
from pathlib import Path
from typing import TYPE_CHECKING, Any, BinaryIO, cast

from .env import load_env

if TYPE_CHECKING:
    import pandas as pd

# Bucket helpers are constructed at module import time, so load the repository
# configuration before resolving their names. Without this, direct imports such
# as ``from deliverome_analysis.s3 import raw`` silently bind to the conceptual
# defaults (``raw``, ``processed``, ``external``) instead of the globally unique
# AWS bucket names configured in ``.env``.
load_env()

# A key under processed must live inside an experiment folder:
#   sprint{N}/s{N}e{M}-{slug}/<anything>
# Mirrors the analysis/ on-disk layout, the raw/ layout, and the experiment
# naming convention enforced by .github/workflows/enforce-experiment-naming.yml.
_PROCESSED_KEY_PATTERN = re.compile(r"^sprint\d+/s\d+e\d+-[a-z0-9]+(-[a-z0-9]+)*/.+$")

# Allowlisted prefixes that bypass the experiment-folder pattern check:
#   _ci_smoke/   — transient probe objects written by CI smoke tests
#   scratch/     — a throwaway area with per-user subfolders (e.g. scratch/becca/)
#                  for onboarding + quick exploration, not tied to an experiment ID
#   constructs/  — the machine-readable mirror of the plasmid registry. These are
#                  inputs shared across experiments, not the output of any one, so
#                  the experiment-folder shape does not fit. Benchling stays the
#                  canonical source of truth; this is the copy analyses (and
#                  agents) read without anything living on a laptop. It is in the
#                  processed bucket on purpose: processed allows overwrite, so a
#                  re-export can correct a map, whereas external is append-only
#                  and a mis-named key there is permanent. Layout is by status,
#                  read from the Plasmid inventory: constructs/ordered/p003.gb,
#                  constructs/designed/p0076.gb, and third-party maps under
#                  constructs/reference/{source}/, e.g. reference/addgene/52961.gbk.
# Anything matched here is still confined to the processed bucket — the
# exemption is only about the in-bucket key shape, not access scope.
# db/: the rebuilt SQLite indexes (deliverome_analysis.qcdb). Only the build writes there, and
# every row in them is derived from files inside experiment folders.
_PROCESSED_KEY_EXEMPT_PREFIXES: tuple[str, ...] = ("_ci_smoke/", "scratch/", "constructs/", "db/")

# A key under external must be source- and version-scoped:
#   {source}/{version}/<anything>          e.g. uniprot/2026_03/human.tsv
# The version segment is what makes external append-only in practice: a new
# snapshot goes to a NEW version prefix rather than overwriting the old one.
# Scientists have no delete on this bucket (explicit IAM Deny), so writing to
# an existing key would only pile up object versions nobody can clean —
# publish under a fresh {version} instead.
_EXTERNAL_KEY_PATTERN = re.compile(r"^[a-z0-9][a-z0-9._-]*/[a-z0-9][a-z0-9._-]*/.+$")

_CACHE_ROOT = Path(
    os.environ.get(
        "DELIVEROME_ANALYSIS_S3_CACHE",
        Path.home() / ".cache" / "deliverome-analysis" / "s3",
    )
)

# Offline mirror: a local folder laid out as ``<mirror>/<raw|processed|external>/<key>``.
# With DELIVEROME_S3_MIRROR set, reads come from the folder and never reach S3, so a
# frozen input bundle can be replayed without credentials (the s4e9 reviewer package).
# DELIVEROME_S3_MIRROR_MODE=record instead reads from S3 as usual and copies each object
# it reads into the folder, which is how such a bundle is collected. Writes always go to
# S3 and are refused in read mode, so a replay cannot publish anything.
_MIRROR_ENV = "DELIVEROME_S3_MIRROR"
_MIRROR_MODE_ENV = "DELIVEROME_S3_MIRROR_MODE"


def _mirror() -> tuple[Path, str] | None:
    """(mirror folder, "read" | "record") when an offline mirror is configured."""
    root = os.environ.get("DELIVEROME_S3_MIRROR")
    if not root:
        return None
    mode = os.environ.get("DELIVEROME_S3_MIRROR_MODE") or "read"
    if mode not in ("read", "record"):
        raise ValueError(f"{_MIRROR_MODE_ENV} must be 'read' or 'record', not {mode!r}")
    return Path(root), mode


# ---------------------------------------------------------------------------
# Account-level plumbing — region + optional endpoint URL + client factory.
# AWS S3: boto3 auto-builds the endpoint from the region. S3-compatible
# providers (R2, Spaces, MinIO, …): set S3_ENDPOINT_URL explicitly.
# ---------------------------------------------------------------------------


def _region() -> str:
    """AWS region. Default us-east-1; override with S3_REGION."""
    load_env()
    return os.environ.get("S3_REGION") or os.environ.get("AWS_REGION") or "us-east-1"


def _endpoint_url() -> str | None:
    """Optional endpoint override for S3-compatible providers.

    Returns None for plain AWS S3 (boto3 builds the endpoint from the
    region). Returns the override for R2 / Spaces / MinIO / etc.
    """
    load_env()
    return os.environ.get("S3_ENDPOINT_URL") or None


def _build_client(access_key: str | None = None, secret_key: str | None = None) -> Any:
    """Build an S3 client.

    With explicit keys, use them. Without, omit the credential kwargs entirely so
    boto3 falls back to its default chain -- which is what lets CI authenticate by
    assuming a role via GitHub OIDC instead of holding a long-lived key.
    """
    import boto3
    from botocore.config import Config

    kwargs: dict[str, Any] = {
        "service_name": "s3",
        "region_name": _region(),
        "config": Config(
            signature_version="s3v4",
            retries={"max_attempts": 5, "mode": "standard"},
            s3={"addressing_style": "virtual"},
        ),
    }
    if access_key and secret_key:
        kwargs["aws_access_key_id"] = access_key
        kwargs["aws_secret_access_key"] = secret_key

    endpoint = _endpoint_url()
    if endpoint:
        kwargs["endpoint_url"] = endpoint

    return boto3.client(**kwargs)


def credentials_available() -> bool:
    """True if boto3 can resolve credentials by any means.

    Covers the explicit ``S3_ACCESS_KEY_ID`` pair, standard ``AWS_*`` variables, an
    instance profile, and an assumed web-identity role. Tests gate on this rather than
    on the presence of two specific variables, so migrating CI to OIDC cannot silently
    turn the live suite into a no-op.
    """
    import boto3

    load_env()
    if os.environ.get("S3_ACCESS_KEY_ID") and os.environ.get("S3_SECRET_ACCESS_KEY"):
        return True
    try:
        return boto3.Session().get_credentials() is not None
    except Exception:
        # Unresolvable credentials are the answer to this question, not an error.
        return False


@lru_cache(maxsize=1)
def client() -> Any:
    """Boto3 S3 client, preferring ``S3_ACCESS_KEY_ID`` / ``S3_SECRET_ACCESS_KEY``.

    One credential per scientist, least-privilege per bucket: read-only on ``raw``,
    read+write-never-delete on ``external``, read+write on ``processed``. So ``raw``
    and ``external`` stay immutable server-side even though a single key backs all
    three accessors.

    When those variables are absent, fall back to boto3's default credential chain.
    That is how CI authenticates -- it assumes a role via GitHub OIDC and holds no
    long-lived key. Scientists are unaffected: an explicit pair always wins.
    """
    load_env()
    access_key = os.environ.get("S3_ACCESS_KEY_ID")
    secret_key = os.environ.get("S3_SECRET_ACCESS_KEY")
    if access_key and secret_key:
        return _build_client(access_key, secret_key)

    if not credentials_available():
        raise RuntimeError(
            "No AWS credentials found. Set S3_ACCESS_KEY_ID and S3_SECRET_ACCESS_KEY "
            "in .env (see .env.example), or provide credentials another way boto3 "
            "understands -- an assumed role, an instance profile, or AWS_* variables."
        )
    return _build_client()


# ---------------------------------------------------------------------------
# Bucket accessors
# ---------------------------------------------------------------------------


def _is_not_found(exc: Exception) -> bool:
    resp = getattr(exc, "response", None)
    if not isinstance(resp, dict):
        return False
    code = str(resp.get("Error", {}).get("Code", ""))
    status = resp.get("ResponseMetadata", {}).get("HTTPStatusCode")
    return code in {"404", "NoSuchKey", "NotFound"} or status == 404


def _refuse_write_in_mirror(key: str) -> None:
    m = _mirror()
    if m is not None and m[1] == "read":
        raise RuntimeError(
            f"refusing to upload {key!r}: {_MIRROR_ENV} is set, so this is an offline replay"
        )


def _provenance_metadata(key: str) -> dict[str, str]:
    """S3 user metadata: the commit, dirty flag, script and lock hashes that made the file."""
    from deliverome_analysis.provenance import provenance, warn_if_dirty

    stamp = provenance()
    warn_if_dirty(stamp, key)
    return stamp


class _ReadBucket:
    """A read-only accessor for one S3 bucket.

    Exposes list / download / read / open / exists. No upload. The bound
    ``client_factory`` callable returns the boto3 client to use for reads on
    this bucket. All buckets share the single ``client()`` factory (one
    credential per scientist); the param is kept for flexibility but there's
    one identity now.
    """

    def __init__(
        self, name: str, client_factory: Callable[[], Any], logical: str | None = None
    ) -> None:
        self._name = name
        self._client_factory = client_factory
        # The conceptual name (raw / processed / external): the mirror's folder name,
        # independent of the deployment's globally unique bucket names.
        self._logical = logical or name

    def _mirror_path(self, key: str) -> Path | None:
        """The object's path in a read-mode mirror, or None when reads go to S3."""
        m = _mirror()
        if m is None or m[1] != "read":
            return None
        return m[0] / self._logical / key

    def _record(self, key: str, path: Path) -> None:
        m = _mirror()
        if m is not None and m[1] == "record":
            target = m[0] / self._logical / key
            if not target.exists():
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(path, target)

    # Public so a power-user can drop down to raw boto3 if needed — but file an
    # issue first; we'd rather add the helper.
    @property
    def name(self) -> str:
        return self._name

    def list_keys(self, prefix: str = "", *, max_keys: int | None = None) -> list[str]:
        """List object keys under ``prefix`` (empty = whole bucket)."""
        local = self._mirror_path("")
        if local is not None:
            found = sorted(p.relative_to(local).as_posix() for p in local.rglob("*") if p.is_file())
            found = [k for k in found if k.startswith(prefix)]
            return found[:max_keys] if max_keys is not None else found
        s3 = self._client_factory()
        keys: list[str] = []
        continuation: str | None = None
        while True:
            kwargs: dict[str, Any] = {"Bucket": self._name, "Prefix": prefix}
            if continuation:
                kwargs["ContinuationToken"] = continuation
            resp = s3.list_objects_v2(**kwargs)
            for obj in resp.get("Contents", []):
                keys.append(obj["Key"])
                if max_keys is not None and len(keys) >= max_keys:
                    return keys
            if not resp.get("IsTruncated"):
                return keys
            continuation = resp.get("NextContinuationToken")

    def exists(self, key: str) -> bool:
        local = self._mirror_path(key)
        if local is not None:
            return local.is_file()
        try:
            self._head(key)
            return True
        except Exception as exc:
            if _is_not_found(exc):
                return False
            raise

    def download(self, key: str, dest: str | Path | None = None) -> Path:
        """Download an object to local disk. Returns the path.

        Caches by ETag at ``~/.cache/deliverome-analysis/s3/<bucket>/<key>.<etag>``.
        In a read-mode offline mirror the mirror's own file is the cache.
        """
        local = self._mirror_path(key)
        if local is not None:
            if not local.is_file():
                raise FileNotFoundError(
                    f"{self._logical}://{key} is not in the offline mirror {local.parents[0]}"
                    f" ({_MIRROR_ENV}); the bundle does not carry this input"
                )
            if dest is None:
                return local
            dest_path = Path(dest)
            dest_path.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(local, dest_path)
            return dest_path
        head = self._head(key)
        etag = head.get("ETag")
        cached = self._cache_path(key, etag)
        if not cached.exists():
            cached.parent.mkdir(parents=True, exist_ok=True)
            fd, tmp_name = tempfile.mkstemp(
                dir=cached.parent, prefix=f".{cached.name}.", suffix=".part"
            )
            os.close(fd)
            tmp = Path(tmp_name)
            try:
                self._client_factory().download_file(self._name, key, str(tmp))
                os.replace(tmp, cached)
            except BaseException:
                tmp.unlink(missing_ok=True)
                raise
        self._record(key, cached)

        if dest is None:
            return cached
        dest_path = Path(dest)
        dest_path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(cached, dest_path)
        return dest_path

    def open_raw(self, key: str) -> BinaryIO:
        """Return a binary file handle for a large object (FASTQ, BAM, FCS).

        Streams via the local cache, so the handle survives after the call.
        """
        return self.download(key).open("rb")

    def read_tsv(self, key: str, **kwargs: Any) -> pd.DataFrame:
        """Read a tab-separated table into a pandas DataFrame.

        kwargs are forwarded to ``pandas.read_csv``.
        """
        import pandas as pd

        kwargs.setdefault("sep", "\t")
        return pd.read_csv(io.BytesIO(self._read_bytes(key)), **kwargs)

    def read_csv(self, key: str, **kwargs: Any) -> pd.DataFrame:
        import pandas as pd

        return pd.read_csv(io.BytesIO(self._read_bytes(key)), **kwargs)

    def read_parquet(self, key: str, **kwargs: Any) -> pd.DataFrame:
        import pandas as pd

        return pd.read_parquet(io.BytesIO(self._read_bytes(key)), **kwargs)

    # ---- internal -----------------------------------------------------

    def list_etags(self, prefix: str = "") -> dict[str, str]:
        """``{key: ETag}`` for every object under ``prefix``, from one paginated listing.

        ETags change whenever an object is rewritten, so they fingerprint content without
        downloading it (used by the incremental index build in ``deliverome_analysis.qcdb``).
        """
        if self._mirror_path("") is not None:
            raise RuntimeError("list_etags is not available in an offline mirror")
        out: dict[str, str] = {}
        pages = self._client_factory().get_paginator("list_objects_v2")
        for page in pages.paginate(Bucket=self._name, Prefix=prefix):
            for obj in page.get("Contents", []):
                out[obj["Key"]] = obj["ETag"].strip('"')
        return out

    def etag(self, key: str) -> str | None:
        """The object's current ETag, or None if it does not exist."""
        if self._mirror_path(key) is not None:
            raise RuntimeError("etag is not available in an offline mirror")
        try:
            return cast(str, self._head(key)["ETag"]).strip('"')
        except Exception as exc:
            if _is_not_found(exc):
                return None
            raise

    def head(self, key: str) -> dict[str, Any] | None:
        """The object's S3 HEAD response (``ContentLength``, ``ETag``, ``Metadata``...), or None.

        None means the object does not exist (404). Any other failure (no permission,
        network) raises, so a caller can tell "absent" from "could not check". Not available
        in an offline mirror, which has no S3 metadata.
        """
        if self._mirror_path(key) is not None:
            raise RuntimeError("head is not available in an offline mirror")
        try:
            return self._head(key)
        except Exception as exc:
            if _is_not_found(exc):
                return None
            raise

    def _head(self, key: str) -> dict[str, Any]:
        return cast(dict[str, Any], self._client_factory().head_object(Bucket=self._name, Key=key))

    def _read_bytes(self, key: str) -> bytes:
        return self.download(key).read_bytes()

    def _cache_path(self, key: str, etag: str | None) -> Path:
        suffix = (etag or "noetag").strip('"').replace("/", "_")
        return _CACHE_ROOT / self._name / f"{key}.{suffix}"

    def __repr__(self) -> str:
        return f"<S3 read-only bucket {self._name!r}>"


class _WritableBucket(_ReadBucket):
    """A read+write accessor for one S3 bucket.

    Inherits everything from ``_ReadBucket`` and adds ``upload``. Backed by
    ``client()``, whose IAM policy grants write access to this (``processed``)
    bucket only — writes to other buckets fail server-side at S3.
    """

    def upload(
        self,
        local_path: str | Path,
        key: str,
        *,
        content_type: str | None = None,
    ) -> str:
        """Upload a local file to the bucket at ``key``. Returns the key.

        Refuses (raises ``ValueError``) if ``key`` doesn't fit the
        experiment-folder shape ``sprint{N}/s{N}e{M}-{slug}/<anything>``.
        That layout mirrors the on-disk ``analysis/`` and ``raw://``
        conventions, so processed outputs always trace back to one experiment
        without having to inspect the file.

        Bypassed only by the prefixes in
        :data:`_PROCESSED_KEY_EXEMPT_PREFIXES` (CI smoke probes).
        """
        if not (
            key.startswith(_PROCESSED_KEY_EXEMPT_PREFIXES) or _PROCESSED_KEY_PATTERN.match(key)
        ):
            raise ValueError(
                f"refusing to upload to processed://{key!r}: keys must match "
                f"the experiment-folder layout sprint{{N}}/s{{N}}e{{M}}-{{slug}}/... "
                f"(e.g. 'sprint1/s1e3-pilot-screen/hits.tsv'). This mirrors the "
                f"analysis/ and raw:// layouts so every processed output traces "
                f"back to a single experiment. If you genuinely need a non-"
                f"experiment key, allow-list it in s3._PROCESSED_KEY_EXEMPT_PREFIXES."
            )

        _refuse_write_in_mirror(key)
        local = Path(local_path)
        if not local.is_file():
            raise FileNotFoundError(f"upload source does not exist: {local}")

        extra: dict[str, Any] = {"Metadata": _provenance_metadata(key)}
        if content_type:
            extra["ContentType"] = content_type

        client().upload_file(
            str(local),
            self._name,
            key,
            ExtraArgs=extra,
        )
        return key

    def upload_if(
        self,
        local_path: str | Path,
        key: str,
        *,
        if_match: str | None,
        content_type: str | None = None,
    ) -> str:
        """Upload only if ``key`` is unchanged since its ETag was ``if_match`` (None: absent).

        S3 conditional write: if another writer replaced (or created) the object in between,
        S3 refuses with ``PreconditionFailed`` / ``ConditionalRequestConflict`` and nothing is
        overwritten. Single PUT, so for files up to 5 GB. Same key rules as :meth:`upload`.
        """
        if not (
            key.startswith(_PROCESSED_KEY_EXEMPT_PREFIXES) or _PROCESSED_KEY_PATTERN.match(key)
        ):
            raise ValueError(f"refusing to upload to processed://{key!r}: not an allowed key")
        _refuse_write_in_mirror(key)
        local = Path(local_path)
        if not local.is_file():
            raise FileNotFoundError(f"upload source does not exist: {local}")
        extra: dict[str, Any] = {"Metadata": _provenance_metadata(key)}
        if content_type:
            extra["ContentType"] = content_type
        extra |= {"IfMatch": if_match} if if_match else {"IfNoneMatch": "*"}
        with local.open("rb") as body:
            client().put_object(Bucket=self._name, Key=key, Body=body, **extra)
        return key

    def _check_key(self, key: str, verb: str) -> None:
        if not (
            key.startswith(_PROCESSED_KEY_EXEMPT_PREFIXES) or _PROCESSED_KEY_PATTERN.match(key)
        ):
            raise ValueError(f"refusing to {verb} processed://{key!r}: not an allowed key")
        _refuse_write_in_mirror(key)

    def copy(self, src_key: str, dst_key: str) -> str:
        """Server-side copy of ``src_key`` to ``dst_key`` in this bucket (any size; metadata,
        including the provenance stamp, is kept). Returns ``dst_key``."""
        self._check_key(dst_key, "copy to")
        client().copy(
            {"Bucket": self._name, "Key": src_key},
            self._name,
            dst_key,
            ExtraArgs={"MetadataDirective": "COPY"},
        )
        return dst_key

    def delete(self, key: str) -> None:
        """Delete ``key``. processed:// is versioned, so the previous version stays
        recoverable; prefer moving (copy, then delete) over deleting results outright."""
        self._check_key(key, "delete")
        client().delete_object(Bucket=self._name, Key=key)

    def __repr__(self) -> str:
        return f"<S3 read+write bucket {self._name!r}>"


class _AppendOnlyBucket(_ReadBucket):
    """A read + write (never delete) accessor for the ``external`` bucket.

    Adds ``upload`` but no delete of any kind. That is a convenience mirror of
    the server-side rule, not the rule itself: the IAM policy carries an
    explicit ``Deny`` on ``s3:DeleteObject`` / ``s3:DeleteObjectVersion`` for
    this bucket, so a delete fails at S3 even via raw boto3.

    The bucket is versioned, so re-writing an existing key supersedes it
    rather than destroying it — but prefer a fresh ``{version}`` prefix, since
    nobody holding a scientist key can clean up superseded versions.
    """

    def upload(
        self,
        local_path: str | Path,
        key: str,
        *,
        content_type: str | None = None,
    ) -> str:
        """Upload a local file to the external bucket at ``key``.

        Refuses (raises ``ValueError``) if ``key`` doesn't fit the
        source/version shape ``{source}/{version}/<anything>`` — e.g.
        ``uniprot/2026_03/human.tsv``. Third-party reference data is only
        useful when you can tell which snapshot you joined against.
        """
        if not _EXTERNAL_KEY_PATTERN.match(key):
            raise ValueError(
                f"refusing to upload to external://{key!r}: keys must match "
                f"the source/version layout {{source}}/{{version}}/... "
                f"(e.g. 'uniprot/2026_03/human.tsv'). Reference data without a "
                f"version prefix can't be traced to a snapshot, and you cannot "
                f"delete it afterwards to fix the name — scientists have no "
                f"delete on this bucket."
            )

        _refuse_write_in_mirror(key)
        local = Path(local_path)
        if not local.is_file():
            raise FileNotFoundError(f"upload source does not exist: {local}")

        extra: dict[str, Any] = {"Metadata": _provenance_metadata(key)}
        if content_type:
            extra["ContentType"] = content_type

        client().upload_file(
            str(local),
            self._name,
            key,
            ExtraArgs=extra,
        )
        return key

    def __repr__(self) -> str:
        return f"<S3 append-only bucket {self._name!r}>"


# Module-level accessors. Bucket names default to the conceptual names but
# S3 bucket names are GLOBALLY unique across AWS, so the user will almost
# always override these to something prefixed (e.g. `deliverome-raw`).
#
# Load .env BEFORE resolving the names below. Each bucket name is read once,
# here at import time, so the S3_BUCKET_* values must already be in
# os.environ. Without this, importing `s3` before calling load_env() (as the
# live tests and most scripts naturally do) silently falls back to the bare
# default names and every call hits a bucket we don't own. load_env() is
# idempotent and never overrides shell-exported vars, so this is safe.
load_env()

# All three accessors share the single client() factory. Per-bucket access
# control is enforced server-side by the IAM policy attached to the one key
# pair (S3_ACCESS_KEY_ID / S3_SECRET_ACCESS_KEY): read-only on raw,
# read+write-no-delete on external, read+write+delete on processed.
raw = _ReadBucket(os.environ.get("S3_BUCKET_RAW", "raw"), client_factory=client, logical="raw")
processed = _WritableBucket(
    os.environ.get("S3_BUCKET_PROCESSED", "processed"),
    client_factory=client,
    logical="processed",
)
external = _AppendOnlyBucket(
    os.environ.get("S3_BUCKET_EXTERNAL", "external"),
    client_factory=client,
    logical="external",
)


# ---------------------------------------------------------------------------
# Module-level cache utility
# ---------------------------------------------------------------------------


def clear_cache() -> None:
    """Wipe the local S3 cache. Useful in CI or after schema changes."""
    if _CACHE_ROOT.exists():
        shutil.rmtree(_CACHE_ROOT)
