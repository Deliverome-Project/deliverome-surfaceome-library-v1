"""Published spectral libraries for DIA searches (alphaDIA), stored like the search FASTA.

A search is only reproducible if the library it ran against can be fetched again. Libraries are
large (the full predicted human library is 2.3 GB) and live in the append-only external bucket:

    external://spectral-library/<name>/speclib.hdf
    external://spectral-library/<name>/manifest.json
    external://spectral-library/<name>/sources/...   # configs, logs, build records

The manifest records the library's SHA-256 and size, the published search FASTA it was predicted
from (by name and SHA-256), how it was made (prediction settings and package versions, or for a
subset the source library, runs and thresholds), and its precursor counts. Search records name
their library by SHA-256 only; `find(sha256)` maps that back to a published library. A published
library is never overwritten. See docs/spectral-library/README.md.

    uv run python -m deliverome_analysis.spectral_library publish speclib.hdf \\
        --name alphadia-2.4.0_human-sp-iso_cont_irt_2026_03_predicted --kind predicted \\
        --fasta human-sp-iso_cont-universal-and-cell-culture_irt/2026_03 \\
        --build build.json --source frozen_config.yaml --source log.txt
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import tempfile
from datetime import UTC, datetime
from pathlib import Path

from deliverome_analysis.search_fasta import sha256

PREFIX = "spectral-library/"
LIBRARY_FILE = "speclib.hdf"
KINDS = ("predicted", "subset")
_NAME = re.compile(r"^[a-z0-9][a-z0-9._-]*$")


def build_manifest(
    library: str | os.PathLike,
    *,
    name: str,
    kind: str,
    fasta_name: str,
    fasta_sha256: str,
    build: dict,
    sources: dict[str, str],
) -> dict:
    """Manifest for one library: identity, FASTA, how it was made, archived sources."""
    if not _NAME.match(name):
        raise ValueError(f"library name {name!r} must be lower-case [a-z0-9._-]")
    if kind not in KINDS:
        raise ValueError(f"kind must be one of {KINDS}")
    lib = Path(library)
    return {
        "name": name,
        "kind": kind,
        "library": {"file": LIBRARY_FILE, "sha256": sha256(lib), "bytes": lib.stat().st_size},
        "search_fasta": {"name": fasta_name, "sha256": fasta_sha256},
        "build": build,
        "sources": sources,
        "published_utc": datetime.now(UTC).isoformat(timespec="seconds"),
    }


def publish(
    library: str | os.PathLike,
    *,
    name: str,
    kind: str,
    fasta_name: str,
    build: dict,
    sources: list[str | os.PathLike] = (),
    dry_run: bool = False,
) -> dict:
    """Upload a library to ``external://spectral-library/<name>/``; refuse if it exists.

    The FASTA must be a published search database: its SHA-256 is read from that snapshot's
    manifest, and if the build record names a FASTA SHA-256 the two must agree.
    """
    from deliverome_analysis.s3 import external

    prefix = f"{PREFIX}{name}/"
    if external.list_keys(prefix=prefix):
        raise FileExistsError(f"external://{prefix} already exists; libraries are immutable")
    fasta_manifest = json.loads(
        Path(external.download(f"search-fasta/{fasta_name}/manifest.json")).read_text()
    )
    fasta_sha = fasta_manifest["search_fasta"]["sha256"]
    if build.get("fasta_sha256") not in (None, fasta_sha):
        raise ValueError(f"build record FASTA {build['fasta_sha256'][:12]} is not {fasta_name}")
    src = {f"sources/{Path(s).name}": Path(s) for s in sources}
    if len(src) != len(list(sources)):
        raise ValueError("source files must have distinct names")
    manifest = build_manifest(
        library,
        name=name,
        kind=kind,
        fasta_name=fasta_name,
        fasta_sha256=fasta_sha,
        build=build,
        sources={k: sha256(p) for k, p in src.items()},
    )
    if dry_run:
        return manifest
    for key, path in src.items():
        external.upload(path, prefix + key)
    external.upload(library, prefix + LIBRARY_FILE)
    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d) / "manifest.json"
        tmp.write_text(json.dumps(manifest, indent=2) + "\n")
        external.upload(tmp, prefix + "manifest.json", content_type="application/json")
    return manifest


def manifest(name: str) -> dict:
    from deliverome_analysis.s3 import external

    return json.loads(Path(external.download(f"{PREFIX}{name}/manifest.json")).read_text())


def download(name: str, dest: str | os.PathLike | None = None) -> Path:
    """Fetch a published library, checking it against its manifest's SHA-256."""
    from deliverome_analysis.s3 import external

    man = manifest(name)
    path = Path(external.download(f"{PREFIX}{name}/{LIBRARY_FILE}"))
    if sha256(path) != man["library"]["sha256"]:
        raise ValueError(f"{LIBRARY_FILE} for {name} does not match its manifest")
    if dest is not None:
        shutil.copyfile(path, dest)
        return Path(dest)
    return path


def published() -> dict[str, dict]:
    """name -> manifest for every published library."""
    from deliverome_analysis.s3 import external

    names = {k[len(PREFIX) :].split("/")[0] for k in external.list_keys(prefix=PREFIX)}
    return {n: manifest(n) for n in sorted(names)}


def find(library_sha256: str) -> str | None:
    """Name of the published library with this SHA-256 (as search records store it)."""
    for name, man in published().items():
        if man["library"]["sha256"] == library_sha256:
            return name
    return None


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("publish")
    p.add_argument("library")
    p.add_argument("--name", required=True)
    p.add_argument("--kind", required=True, choices=KINDS)
    p.add_argument("--fasta", required=True, help="published search-fasta name")
    p.add_argument("--build", required=True, help="JSON file: how the library was made")
    p.add_argument("--source", action="append", default=[], help="file archived under sources/")
    p.add_argument("--dry-run", action="store_true")
    sub.add_parser("list")
    f = sub.add_parser("find")
    f.add_argument("sha256")
    a = ap.parse_args()
    if a.cmd == "publish":
        m = publish(
            a.library,
            name=a.name,
            kind=a.kind,
            fasta_name=a.fasta,
            build=json.loads(Path(a.build).read_text()),
            sources=a.source,
            dry_run=a.dry_run,
        )
        print(json.dumps(m, indent=2))
    elif a.cmd == "list":
        for n, m in published().items():
            print(n, m["kind"], m["library"]["sha256"][:12], m["library"]["bytes"])
    else:
        print(find(a.sha256))
