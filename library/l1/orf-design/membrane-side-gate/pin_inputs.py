"""Pin every external input the membrane-side gate reads, so a clean clone can rerun it.

The gate and its calibration read four kinds of input that are not in Git:

| Input | Where it comes from | Pinned how |
|---|---|---|
| S3 objects (topology snapshot, designed-CDS record, ColabFold models, gold standard) | `external://` and `processed://` | SHA-256 recorded; the topology snapshot is also checked against a pin in code |
| AlphaFold DB v6 models | alphafold.ebi.ac.uk (public, re-issued between versions) | copy bundled, SHA-256 recorded |
| TmAlphaFold TMDET entries and transformed models | tmalphafold.ttk.hu (live service) | copy bundled, SHA-256 recorded |
| Tag-site API responses, UniProt sequences and isoform lists, EBI features | live web services | copy bundled, SHA-256 recorded |

`pin` writes ``input_manifest.tsv`` (one row per file: kind, file, size, SHA-256, and the
object key or URL it came from) and ``inputs_bundle.tar.gz`` of every non-S3 file. The bundle
is uploaded once to ``processed://`` under a new dated key. `restore` downloads that bundle,
checks its SHA-256 against the manifest's, and unpacks it into the gate's cache folder, after
which the gate and calibration run with no live service.

    uv run --locked python pin_inputs.py pin --upload-to <processed key prefix>
    uv run --locked python pin_inputs.py restore
"""

from __future__ import annotations

import argparse
import hashlib
import json
import tarfile
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
CACHE = Path.home() / ".cache" / "deliverome-analysis" / "membrane-side-gate"
MANIFEST = HERE / "input_manifest.tsv"
BUNDLE_POINTER = HERE / "input_bundle.json"
#: Cache folders holding copies of live-service responses and public model files.
CACHED_KINDS = {
    "models": "AlphaFold DB model or ColabFold model copy",
    "tmaf": "TmAlphaFold (TMDET) entry or transformed model",
    "api": "accessible-surfaceome tag-site API response",
    "uniprot": "UniProt FASTA or isoform list",
    "ebi": "EBI Proteins API features",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _source(kind: str, name: str) -> str:
    if kind == "models" and name.startswith("AF-"):
        return f"https://alphafold.ebi.ac.uk/files/{name}"
    if kind == "models" and name.startswith("colabfold_"):
        return "processed:// ColabFold fold (key in the member table's model_ref)"
    if kind == "api":
        return f"https://api.deliverome.org/surfaceome/v1/tag-sites/{name.removesuffix('.json')}"
    if kind == "tmaf":
        return "https://tmalphafold.ttk.hu/"
    if kind == "uniprot":
        return "https://rest.uniprot.org/uniprotkb/"
    return "https://www.ebi.ac.uk/proteins/api/features/"


def s3_inputs() -> list[dict]:
    """The S3 objects the gate reads, with the SHA-256 of the copy read."""
    import gate_inputs as gi

    from deliverome_analysis.s3 import external, processed

    rows = []
    for bucket, key in [
        (external, gi.TOPOLOGY_SNAPSHOT_KEY),
        (processed, gi.CDS_KEY),
        # The pins of record on main since 2026-11-12, read by the 2026-11-14 snapshot check.
        (external, "accessible-surfaceome/topology/2026-11-12-cyno-restarts/topology.parquet"),
        (
            processed,
            "sprint4/s4e9-orf-design-fg-panel/library-candidates/"
            "2026-11-12-ortholog-start-trims-2/candidate_cds.json",
        ),
        (processed, "library/l1/tag-site-corpus-audit/gold_standard.tsv"),
        (processed, "library/l1/tag-site-corpus-audit/failed_sites.tsv"),
    ]:
        path = Path(bucket.download(key))
        rows.append(
            dict(
                kind="s3",
                file=key,
                size=path.stat().st_size,
                sha256=sha256(path),
                source=f"{'external' if bucket is external else 'processed'}://{key}",
            )
        )
    return rows


def pin(upload_to: str | None) -> int:
    rows = s3_inputs()
    files = []
    for kind in CACHED_KINDS:
        d = CACHE / kind
        if not d.is_dir():
            continue
        for f in sorted(d.iterdir()):
            if f.is_file() and f.stat().st_size > 0:
                files.append(f)
                rows.append(
                    dict(
                        kind=kind,
                        file=f"{kind}/{f.name}",
                        size=f.stat().st_size,
                        sha256=sha256(f),
                        source=_source(kind, f.name),
                    )
                )
    pd.DataFrame(rows).to_csv(MANIFEST, sep="\t", index=False)
    bundle = CACHE / "inputs_bundle.tar.gz"
    with tarfile.open(bundle, "w:gz") as tar:
        for f in files:
            tar.add(f, arcname=str(f.relative_to(CACHE)))
    digest = sha256(bundle)
    pointer = dict(bundle_sha256=digest, files=len(files), bytes=bundle.stat().st_size)
    if upload_to:
        from deliverome_analysis.s3 import processed

        key = upload_to.rstrip("/") + "/inputs_bundle.tar.gz"
        if processed.list_keys(prefix=key):
            raise SystemExit(f"{key} exists; pin under a new dated prefix")
        processed.upload(bundle, key)
        pointer["key"] = key
    BUNDLE_POINTER.write_text(json.dumps(pointer, indent=1) + "\n")
    print(json.dumps(pointer, indent=1), f"\nmanifest: {len(rows)} rows")
    return 0


def restore() -> int:
    from deliverome_analysis.s3 import processed

    pointer = json.loads(BUNDLE_POINTER.read_text())
    CACHE.mkdir(parents=True, exist_ok=True)
    bundle = CACHE / "inputs_bundle.tar.gz"
    processed.download(pointer["key"], bundle)
    if sha256(bundle) != pointer["bundle_sha256"]:
        raise SystemExit("bundle SHA-256 does not match input_bundle.json")
    with tarfile.open(bundle) as tar:
        tar.extractall(CACHE, filter="data")
    manifest = pd.read_csv(MANIFEST, sep="\t")
    bad = [
        r.file
        for r in manifest[manifest["kind"] != "s3"].itertuples()
        if sha256(CACHE / r.file) != r.sha256
    ]
    if bad:
        raise SystemExit(f"{len(bad)} restored files differ from the manifest: {bad[:5]}")
    print(f"restored {len(manifest[manifest['kind'] != 's3'])} files into {CACHE}")
    return 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("pin")
    p.add_argument("--upload-to", help="processed:// key prefix (new, dated) for the bundle")
    sub.add_parser("restore")
    args = ap.parse_args(argv)
    return pin(args.upload_to) if args.cmd == "pin" else restore()


if __name__ == "__main__":
    raise SystemExit(main())
