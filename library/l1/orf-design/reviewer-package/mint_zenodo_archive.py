"""Stage, verify, and mint an immutable 5-year archival release of the L1 package to Zenodo.

Bundles the complete offline replay artifacts into one checksum-verified folder:

1. ``deliverome-surfaceome-library-v1.bundle`` — full git bundle of the exported review repository;
2. ``l1-reviewer-bundle.tar.gz`` — frozen S3 inputs required by ``run_review.py``;
3. ``l1-pinned-releases.tar`` (optional) — remaining pinned upstream releases;
4. ``deliverome-surfaceome-library-v1-image.tar.gz`` (optional) — frozen OCI container built from ``Dockerfile``;
5. ``RERUN.md`` — three-command offline replay instructions for 2031+;
6. ``zenodo_metadata.json`` — Zenodo deposition metadata;
7. ``SHA256SUMS`` — SHA-256 of every file above.

Upload creates a **draft** deposition by default (printing its URL and reserved DOI) so a human
can inspect the record before publishing; pass ``--publish`` to mint the permanent DOI, or
``--sandbox`` to test against ``https://sandbox.zenodo.org``.

    uv run --locked python mint_zenodo_archive.py stage \
        --review-repo ~/git/deliverome-surfaceome-library-v1 \
        --bundle ~/l1-reviewer-bundle.tar.gz \
        --out ~/l1-zenodo-archive-2026-11-22 \
        --version 1.0.0
    uv run --locked python mint_zenodo_archive.py verify --archive ~/l1-zenodo-archive-2026-11-22
    uv run --locked python mint_zenodo_archive.py upload --archive ~/l1-zenodo-archive-2026-11-22 --sandbox
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
from datetime import date
from pathlib import Path
from typing import Any

import requests

from deliverome_analysis.core.env import load_env

ZENODO_API = "https://zenodo.org/api"
ZENODO_SANDBOX_API = "https://sandbox.zenodo.org/api"
SUMS_FILE = "SHA256SUMS"
META_FILE = "zenodo_metadata.json"
RERUN_FILE = "RERUN.md"

RERUN_TEMPLATE = """# Offline Replay Instructions (L1 Surface-Protein ORF Library)

Version: {version}
Archived: {archived_on}

## 1. Verify archive integrity

```bash
sha256sum -c SHA256SUMS
```

## 2. Unpack repository and frozen input bundle

```bash
git clone deliverome-surfaceome-library-v1.bundle deliverome-surfaceome-library-v1
mkdir -p l1-reviewer-bundle && tar -xzf l1-reviewer-bundle.tar.gz -C l1-reviewer-bundle
```

## 3. Replay with zero network access

With the archived OCI container image (recommended for 5+ year reproducibility):

```bash
docker load -i deliverome-surfaceome-library-v1-image.tar.gz
docker run --rm --network=none \\
  -v "$(pwd)/l1-reviewer-bundle:/bundle:ro" \\
  deliverome-surfaceome-library-v1:frozen
```

Or directly on a host with Python 3.14.10 and `uv >= 0.13.4`:

```bash
cd deliverome-surfaceome-library-v1
uv sync --locked --no-group flow --no-group flow-gui --no-group modal
uv run --locked python library/l1/orf-design/reviewer-package/run_review.py \\
  --bundle ../l1-reviewer-bundle
```
"""


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def default_metadata(version: str, *, access_right: str = "closed") -> dict[str, Any]:
    """Zenodo deposition metadata block (defaults to closed/restricted access for private review)."""
    return {
        "metadata": {
            "title": f"Deliverome L1 Surface-Protein ORF Library — Frozen Review Archive (v{version})",
            "upload_type": "dataset",
            "version": version,
            "publication_date": date.today().isoformat(),
            "description": (
                "Self-contained offline reproducibility archive for the Deliverome L1 "
                "surface-protein ORF library design: exported review repository git bundle, "
                "SHA-256-pinned input bundle, optional frozen OCI container image, and "
                "checksum manifest."
            ),
            "creators": [{"name": "Deliverome Project"}],
            "access_right": access_right,
            "keywords": ["surfaceome", "ORF library", "functional genomics", "reproducibility"],
        }
    }


def _copy_or_link(src: Path, dest: Path, *, hardlink: bool) -> None:
    if dest.exists():
        dest.unlink()
    if hardlink:
        try:
            os.link(src, dest)
            return
        except OSError:
            pass
    shutil.copy2(src, dest)


def stage_archive(
    *,
    review_repo: Path,
    bundle_tar: Path,
    out_dir: Path,
    version: str,
    pinned_releases_tar: Path | None = None,
    image_tar: Path | None = None,
    access_right: str = "closed",
    hardlink: bool = False,
) -> dict[str, str]:
    """Assemble the archival directory and write ``SHA256SUMS``; returns ``{filename: sha256}``."""
    out_dir.mkdir(parents=True, exist_ok=True)
    staged: list[Path] = []

    git_bundle = out_dir / "deliverome-surfaceome-library-v1.bundle"
    subprocess.run(
        ["git", "-C", str(review_repo), "bundle", "create", str(git_bundle), "--all"],
        check=True,
        capture_output=True,
    )
    staged.append(git_bundle)

    for src, name in (
        (bundle_tar, "l1-reviewer-bundle.tar.gz"),
        (pinned_releases_tar, "l1-pinned-releases.tar"),
        (image_tar, "deliverome-surfaceome-library-v1-image.tar.gz"),
    ):
        if src is None:
            continue
        if not src.is_file():
            raise FileNotFoundError(f"missing archive input: {src}")
        dest = out_dir / name
        _copy_or_link(src, dest, hardlink=hardlink)
        staged.append(dest)

    rerun = out_dir / RERUN_FILE
    rerun.write_text(RERUN_TEMPLATE.format(version=version, archived_on=date.today().isoformat()))
    staged.append(rerun)

    meta = out_dir / META_FILE
    meta.write_text(
        json.dumps(default_metadata(version, access_right=access_right), indent=2) + "\n"
    )
    staged.append(meta)

    sums = {p.name: sha256_file(p) for p in sorted(staged, key=lambda x: x.name)}
    lines = [f"{digest}  {name}" for name, digest in sums.items()]
    (out_dir / SUMS_FILE).write_text("\n".join(lines) + "\n")
    return sums


def verify_archive(archive_dir: Path) -> dict[str, str]:
    """Verify every file recorded in ``SHA256SUMS``; raises ``ValueError`` on any mismatch."""
    sums_path = archive_dir / SUMS_FILE
    if not sums_path.is_file():
        raise FileNotFoundError(f"missing {sums_path}")
    recorded: dict[str, str] = {}
    problems: list[str] = []
    for line in sums_path.read_text().splitlines():
        if not line.strip():
            continue
        digest, name = line.split(None, 1)
        target = archive_dir / name.strip()
        if not target.is_file():
            problems.append(f"missing: {name}")
        elif sha256_file(target) != digest:
            problems.append(f"checksum mismatch: {name}")
        recorded[name.strip()] = digest
    if problems:
        raise ValueError("archive verification failed:\n  " + "\n  ".join(problems))
    return recorded


def upload_archive(
    archive_dir: Path,
    *,
    token: str,
    sandbox: bool = False,
    publish: bool = False,
    session: Any | None = None,
) -> dict[str, Any]:
    """Upload a verified archive to Zenodo as a draft (or publish with ``publish=True``)."""
    recorded = verify_archive(archive_dir)
    meta = json.loads((archive_dir / META_FILE).read_text())
    api = ZENODO_SANDBOX_API if sandbox else ZENODO_API
    http = session or requests.Session()
    headers = {"Authorization": f"Bearer {token}"}

    resp = http.post(f"{api}/deposit/depositions", json={}, headers=headers, timeout=60)
    resp.raise_for_status()
    dep = resp.json()
    dep_id = dep["id"]
    bucket_url = dep["links"]["bucket"]

    for name in [*sorted(recorded), SUMS_FILE]:
        path = archive_dir / name
        with path.open("rb") as fh:
            up = http.put(f"{bucket_url}/{name}", data=fh, headers=headers, timeout=3600)
            up.raise_for_status()

    put_meta = http.put(
        f"{api}/deposit/depositions/{dep_id}", json=meta, headers=headers, timeout=60
    )
    put_meta.raise_for_status()
    result = put_meta.json()

    if publish:
        pub = http.post(
            f"{api}/deposit/depositions/{dep_id}/actions/publish", headers=headers, timeout=120
        )
        pub.raise_for_status()
        result = pub.json()
    return result


def main(argv: list[str] | None = None) -> int:
    load_env()
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    st = sub.add_parser("stage", help="assemble archive folder and SHA256SUMS")
    st.add_argument("--review-repo", type=Path, required=True)
    st.add_argument("--bundle", type=Path, required=True)
    st.add_argument("--out", type=Path, required=True)
    st.add_argument("--version", required=True)
    st.add_argument("--pinned-releases", type=Path, default=None)
    st.add_argument("--image", type=Path, default=None)
    st.add_argument("--access-right", default="closed", choices=["closed", "restricted", "open"])
    st.add_argument(
        "--hardlink", action="store_true", help="hardlink large tarballs instead of copying"
    )

    vf = sub.add_parser("verify", help="verify SHA256SUMS in a staged archive")
    vf.add_argument("--archive", type=Path, required=True)

    up = sub.add_parser("upload", help="upload a staged archive to Zenodo (draft by default)")
    up.add_argument("--archive", type=Path, required=True)
    up.add_argument("--sandbox", action="store_true", help="target sandbox.zenodo.org")
    up.add_argument("--publish", action="store_true", help="publish immediately and mint DOI")

    args = ap.parse_args(argv)
    if args.cmd == "stage":
        sums = stage_archive(
            review_repo=args.review_repo.expanduser(),
            bundle_tar=args.bundle.expanduser(),
            out_dir=args.out.expanduser(),
            version=args.version,
            pinned_releases_tar=args.pinned_releases.expanduser() if args.pinned_releases else None,
            image_tar=args.image.expanduser() if args.image else None,
            access_right=args.access_right,
            hardlink=args.hardlink,
        )
        print(f"staged {len(sums)} files to {args.out}:")
        for name, digest in sums.items():
            print(f"  {digest[:12]}  {name}")
        return 0
    if args.cmd == "verify":
        sums = verify_archive(args.archive.expanduser())
        print(f"verified {len(sums)} files in {args.archive}")
        return 0
    token = os.environ.get("ZENODO_TOKEN", "").strip()
    if not token:
        raise SystemExit("ZENODO_TOKEN is not set (add it to .env or export it)")
    dep = upload_archive(
        args.archive.expanduser(), token=token, sandbox=args.sandbox, publish=args.publish
    )
    state = "published" if args.publish else "draft"
    print(
        f"Zenodo {state} deposition id={dep.get('id')} url={dep.get('links', {}).get('html', '')}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
