"""Where every file of the reviewer bundle lives for a reviewer, and the public copies.

The bundle (``bundle_manifest.tsv``) has three homes, decided per file by the rules below and
written to ``data_locations.tsv``:

* ``public_r2``: a public Cloudflare R2 bucket, for files from public releases whose terms
  allow redistribution (each row names the license). Anyone can download them, no account.
* ``review_repo``: small tables (≤ 5 MB, tabular), committed to the review repository under
  ``data/<bucket>/<key>`` by ``export_review_repo.py --bundle``.
* ``public_source``: a public file too large to upload through the API (over ~300 MB), whose
  original source serves the identical bytes (same size and SHA-256, checked when it was added).
* ``private_release``: everything else (our own predictions and records, and public sources
  whose redistribution terms are unclear or non-commercial), only in the private review
  repository's release tarballs.

``tests/test_reviewer_data_locations.py`` fails unless every manifest row has exactly one home
with a matching size and SHA-256, so no file can be dropped between the three.

    uv run --locked python public_data.py plan                       # write data_locations.tsv
    uv run --locked python public_data.py upload --bundle ~/l1-bundle  # CLOUDFLARE_* in .env
    uv run --locked python public_data.py fetch --out ../l1-reviewer-bundle   # no credentials
    uv run --locked python public_data.py check                      # every public object online

``upload`` reads ``CLOUDFLARE_ACCOUNT_ID`` and ``CLOUDFLARE_API_TOKEN`` (a token with R2 edit
rights) from the repository's ``.env``; ``fetch`` and ``check`` need no credentials.
"""

from __future__ import annotations

import argparse
import csv
import datetime
import email.utils
import hashlib
import os
import re
import sys
import tempfile
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any

HERE = Path(__file__).parent
MANIFEST = HERE / "bundle_manifest.tsv"
LOCATIONS = HERE / "data_locations.tsv"

BUCKET = "l1-review-public-data"
PUBLIC_BASE = "https://pub-e53f37f6fe6b4745829c2066350fe26f.r2.dev"
REPO_MAX_BYTES = 5_000_000
#: r2.dev refuses the default Python-urllib User-Agent (HTTP 403); a plain client name is enough.
USER_AGENT = "l1-review-data/1.0"
TABULAR = re.compile(r"\.(tsv|csv|json|jsonl|txt|fasta|fa)(\.gz)?$")

#: Public sources whose terms allow redistribution of unmodified files, checked 2026-11-21.
#: (pattern on "<bucket>/<key>", license as stated by the provider). Anything not matched
#: here stays private: PeptideAtlas is CC BY-NC 4.0. APPRIS states CC BY-NC-ND 4.0
#: (github.com/appris/appris LICENSE.md, read 2026-11-22), which permits unmodified
#: non-commercial redistribution with attribution; Becca decided on 2026-11-22 to publish the
#: unmodified APPRIS data files ("yes put appris on the public bucket"), so they are matched
#: below and our own files beside them (provenance, availability) are not. OpenProt, the
#: GENCODE Ribo-seq ORF files and the Tedman et al. supplementary tables had no license we
#: could confirm.
PUBLIC_SOURCES: list[tuple[str, str]] = [
    (
        r"^external/alphafold-db/",
        "CC BY 4.0 (AlphaFold DB; cite Jumper et al. 2021, Fleming et al. 2025)",
    ),
    (
        r"^processed/sprint4/s4e9-orf-design-fg-panel/pinned-inputs/afdb-v6/",
        "CC BY 4.0 (AlphaFold DB models, repackaged unmodified)",
    ),
    (r"^external/hpa/", "CC BY 4.0 (Human Protein Atlas; cite the HPA and proteinatlas.org)"),
    (
        r"^external/human-protein-atlas/",
        "CC BY 4.0 (Human Protein Atlas; cite the HPA and proteinatlas.org)",
    ),
    (
        r"^external/ensembl/",
        "EMBL-EBI terms: no additional restrictions on redistribution (Ensembl)",
    ),
    (
        r"^external/uniprot/",
        "EMBL-EBI terms: no additional restrictions on redistribution (UniProt)",
    ),
    (r"^external/pfam/", "EMBL-EBI terms: no additional restrictions on redistribution (Pfam)"),
    (r"^external/accessible-surfaceome/release-", "CC0 (accessible-surfaceome release, Zenodo)"),
    (
        r"^external/appris/.*appris_[A-Za-z_.0-9]+\.txt$",
        "CC BY-NC-ND 4.0 (APPRIS, unmodified; non-commercial use; cite Rodriguez et al. 2022, "
        "appris.bioinfo.cnio.es)",
    ),
]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def manifest_rows() -> list[dict[str, str]]:
    with MANIFEST.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


#: Public files served unchanged at their original source (bundle key -> URL), each checked
#: byte for byte against the bundle's SHA-256 when added. Pfam-A is over the API upload size.
PUBLIC_SOURCE_URLS: dict[str, str] = {
    "external/pfam/39.0/Pfam-A.hmm.gz": (
        "https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam39.0/Pfam-A.hmm.gz"
    ),
}


def home_of(row: dict[str, str]) -> tuple[str, str]:
    """(home, license) for one manifest row."""
    loc = f"{row['bucket']}/{row['key']}"
    for pattern, license_ in PUBLIC_SOURCES:
        if re.match(pattern, loc):
            return ("public_source" if loc in PUBLIC_SOURCE_URLS else "public_r2"), license_
    if TABULAR.search(row["key"]) and int(row["bytes"]) <= REPO_MAX_BYTES:
        return "review_repo", ""
    return "private_release", ""


def location(row: dict[str, str], home: str) -> str:
    if home == "public_r2":
        return f"{PUBLIC_BASE}/{row['bucket']}/{row['key']}"
    if home == "public_source":
        return PUBLIC_SOURCE_URLS[f"{row['bucket']}/{row['key']}"]
    if home == "review_repo":
        return f"data/{row['bucket']}/{row['key']}"
    return f"release asset: s3/{row['bucket']}/{row['key']}"


FIELDS = ["bucket", "key", "bytes", "sha256", "home", "license", "location"]


def plan() -> int:
    out = []
    for r in manifest_rows():
        home, license_ = home_of(r)
        out.append(
            {
                "bucket": r["bucket"],
                "key": r["key"],
                "bytes": r["bytes"],
                "sha256": r["sha256"],
                "home": home,
                "license": license_,
                "location": location(r, home),
            }
        )
    with LOCATIONS.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out)
    counts: dict[str, list[int]] = {}
    for r in out:
        c = counts.setdefault(r["home"], [0, 0])
        c[0] += 1
        c[1] += int(r["bytes"])
    for home, (n, b) in sorted(counts.items()):
        print(f"{home:16s} {n:5d} files {b / 1e9:6.2f} GB")
    return 0


def locations(*homes: str) -> list[dict[str, str]]:
    with LOCATIONS.open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return [r for r in rows if not homes or r["home"] in homes]


def _is_retryable_status(code: int) -> bool:
    return code == 429 or 500 <= code < 600


def _retry_delay(attempt: int, headers: Any = None) -> float:
    """Seconds to wait before retry, honouring numeric or HTTP-date Retry-After."""
    raw = headers.get("Retry-After") if headers is not None else None
    if raw is not None:
        text = str(raw).strip()
        try:
            return max(0.0, float(text))
        except ValueError:
            try:
                dt = email.utils.parsedate_to_datetime(text)
                if dt.tzinfo is None:
                    dt = dt.replace(tzinfo=datetime.UTC)
                return max(0.0, dt.timestamp() - time.time())
            except TypeError, ValueError, OverflowError:
                pass
    return float(2**attempt)


def _head(url: str) -> int | None:
    """Content-Length of ``url``, or None if absent; backs off on rate limits (r2.dev returns 429)."""
    for attempt in range(8):
        req = urllib.request.Request(url, method="HEAD", headers={"User-Agent": USER_AGENT})
        try:
            with urllib.request.urlopen(req, timeout=60) as resp:
                return int(resp.headers.get("Content-Length", "-1"))
        except urllib.error.HTTPError as exc:
            if exc.code == 404:
                return None
            if not _is_retryable_status(exc.code) or attempt == 7:
                raise
            time.sleep(_retry_delay(attempt, exc.headers))
            continue
        time.sleep(2**attempt)
    return None


def _put(path: Path, key: str, account: str, token: str) -> None:
    url = (
        f"https://api.cloudflare.com/client/v4/accounts/{account}/r2/buckets/{BUCKET}/objects/{key}"
    )
    data = path.read_bytes()
    for attempt in range(8):
        req = urllib.request.Request(
            url, data=data, method="PUT", headers={"Authorization": f"Bearer {token}"}
        )
        try:
            with urllib.request.urlopen(req, timeout=900) as resp:
                if resp.status == 200:
                    return
        except urllib.error.HTTPError as exc:
            if not _is_retryable_status(exc.code) or attempt == 7:
                raise
            time.sleep(_retry_delay(attempt, exc.headers))
            continue
        time.sleep(2**attempt)


def _download_file(url: str, dest: Path, *, max_attempts: int = 8) -> None:
    """Download ``url`` to ``dest`` atomically, retrying HTTP 429 and 5xx with backoff."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    for attempt in range(max_attempts):
        req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
        tmp: Path | None = None
        try:
            with urllib.request.urlopen(req, timeout=900) as resp:
                fd, tmp_name = tempfile.mkstemp(
                    dir=dest.parent, prefix=f".{dest.name}.", suffix=".part"
                )
                tmp = Path(tmp_name)
                with os.fdopen(fd, "wb") as fh:
                    for chunk in iter(lambda: resp.read(1 << 20), b""):
                        fh.write(chunk)
            os.replace(tmp, dest)
            return
        except urllib.error.HTTPError as exc:
            if tmp is not None:
                tmp.unlink(missing_ok=True)
            if not _is_retryable_status(exc.code) or attempt == max_attempts - 1:
                raise
            time.sleep(_retry_delay(attempt, exc.headers))
        except BaseException:
            if tmp is not None:
                tmp.unlink(missing_ok=True)
            raise


def upload(bundle: Path, workers: int) -> int:
    account, token = os.environ["CLOUDFLARE_ACCOUNT_ID"], os.environ["CLOUDFLARE_API_TOKEN"]
    rows = locations("public_r2")
    failed: list[str] = []

    def one(r: dict[str, str]) -> str:
        key = f"{r['bucket']}/{r['key']}"
        if _head(r["location"]) == int(r["bytes"]):
            return "skip"
        path = bundle / "s3" / r["bucket"] / r["key"]
        if sha256(path) != r["sha256"]:
            raise SystemExit(f"{path}: SHA-256 differs from the manifest")
        _put(path, key, account, token)
        return "put"

    with ThreadPoolExecutor(workers) as pool:
        futures = {pool.submit(one, r): r for r in rows}
        done = {"put": 0, "skip": 0}
        for fut, r in futures.items():
            try:
                done[fut.result()] += 1
            except Exception as exc:  # report every failure, then fail
                failed.append(f"{r['bucket']}/{r['key']}: {type(exc).__name__}: {exc}")
    print(f"uploaded {done['put']}, already present {done['skip']}, failed {len(failed)}")
    for line in failed:
        print("  " + line, file=sys.stderr)
    return 1 if failed else 0


def fetch(out: Path, workers: int) -> int:
    """Download every public file into ``out/s3/<bucket>/<key>``, verifying SHA-256."""
    rows = locations("public_r2", "public_source")
    bad: list[str] = []

    def one(r: dict[str, str]) -> None:
        dest = out / "s3" / r["bucket"] / r["key"]
        if dest.is_file() and sha256(dest) == r["sha256"]:
            return
        _download_file(r["location"], dest)
        if sha256(dest) != r["sha256"]:
            bad.append(f"{r['location']}: SHA-256 differs")

    with ThreadPoolExecutor(workers) as pool:
        list(pool.map(one, rows))
    print(f"fetched {len(rows)} public files into {out}; {len(bad)} failed verification")
    for line in bad:
        print("  " + line, file=sys.stderr)
    return 1 if bad else 0


def _list_r2_sizes(account: str, token: str) -> dict[str, int]:
    """List all object sizes in ``BUCKET`` via the Cloudflare R2 API (fast, no r2.dev 429s)."""
    import json

    have: dict[str, int] = {}
    cursor: str | None = None
    while True:
        url = (
            f"https://api.cloudflare.com/client/v4/accounts/{account}"
            f"/r2/buckets/{BUCKET}/objects?per_page=1000"
        )
        if cursor:
            url += f"&cursor={cursor}"
        req = urllib.request.Request(url, headers={"Authorization": f"Bearer {token}"})
        data = json.loads(urllib.request.urlopen(req, timeout=30).read())
        if not data.get("success"):
            raise RuntimeError(f"Cloudflare R2 list failed: {data}")
        for obj in data["result"]:
            have[obj["key"]] = int(obj["size"])
        ri = data.get("result_info") or {}
        if not ri.get("is_truncated"):
            break
        cursor = ri.get("cursor")
    return have


def check(workers: int) -> int:
    """Every public file is online, at its bucket or source URL, with the manifest's size."""
    rows = locations("public_r2", "public_source")
    account = os.environ.get("CLOUDFLARE_ACCOUNT_ID")
    token = os.environ.get("CLOUDFLARE_API_TOKEN")
    bad: list[str] = []
    if account and token:
        r2_sizes = _list_r2_sizes(account, token)
        for r in rows:
            if r["home"] == "public_r2":
                if r2_sizes.get(f"{r['bucket']}/{r['key']}") != int(r["bytes"]):
                    bad.append(r["location"])
            elif _head(r["location"]) != int(r["bytes"]):
                bad.append(r["location"])
    else:
        with ThreadPoolExecutor(workers) as pool:
            sizes = list(pool.map(lambda r: _head(r["location"]), rows))
        bad = [r["location"] for r, n in zip(rows, sizes, strict=True) if n != int(r["bytes"])]
    print(f"{len(rows) - len(bad)} of {len(rows)} public objects online at the manifest size")
    for url in bad[:20]:
        print("  missing or wrong size: " + url, file=sys.stderr)
    return 1 if bad else 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("plan")
    up = sub.add_parser("upload")
    up.add_argument("--bundle", type=Path, required=True)
    up.add_argument("--workers", type=int, default=6)
    fe = sub.add_parser("fetch")
    fe.add_argument("--out", type=Path, required=True)
    fe.add_argument("--workers", type=int, default=8)
    ch = sub.add_parser("check")
    ch.add_argument("--workers", type=int, default=4)
    a = ap.parse_args(argv)
    if a.cmd == "plan":
        return plan()
    if a.cmd == "upload":
        from deliverome_analysis.core.env import load_env

        load_env()
        return upload(a.bundle.expanduser(), a.workers)
    if a.cmd == "fetch":
        return fetch(a.out.expanduser(), a.workers)
    return check(a.workers)


if __name__ == "__main__":
    raise SystemExit(main())
