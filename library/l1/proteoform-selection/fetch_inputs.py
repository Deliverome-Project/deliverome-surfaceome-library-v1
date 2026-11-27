"""Download the APPRIS principal-isoform tables and the Ensembl peptide FASTAs they index.

APPRIS release e115v50 is built on Ensembl 115, so the peptide sequences must come from
Ensembl 115 and not from the live REST service (currently serving 116) — a transcript
APPRIS labels can be retired or re-versioned between releases, and a sequence fetched from
the wrong release is silently the wrong answer.

Files land in a local cache (gitignored, ~30 MB) and the two small APPRIS tables are
published to ``external://appris/e115v50/`` so the selection is reproducible from S3.

    uv run --locked python fetch_inputs.py --cache-dir ~/.cache/deliverome/appris
    uv run --locked python fetch_inputs.py --cache-dir ~/.cache/deliverome/appris --publish
"""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import requests

APPRIS_RELEASE = "e115v50"
ENSEMBL_RELEASE = "115"

# Pin the dated release directory, not ``current_release/``, which is a moving symlink that
# serves different content after every APPRIS rebuild. 2026_06.v50 is the first named
# release carrying the e115v50 dataset; its files are byte-identical to the current_release
# copies this analysis was first run on (sha256 in provenance.json).
APPRIS_SNAPSHOT = "2026_06.v50"
APPRIS_BASE = f"https://apprisws.bioinfo.cnio.es/pub/releases/{APPRIS_SNAPSHOT}/datafiles"
APPRIS_FILES = {
    "homo_sapiens": "appris_data.principal.txt",
    "mus_musculus": "appris_data.principal.txt",
}
ENSEMBL_PEP = {
    "homo_sapiens": "Homo_sapiens.GRCh38.pep.all.fa.gz",
    "mus_musculus": "Mus_musculus.GRCm39.pep.all.fa.gz",
}

# APPRIS covers 11 species at e115v50; Macaca fascicularis is not one of them. Our cyno
# orthologs are M. fascicularis (taxon 9541), so no principal-isoform call is available
# for that cohort and none is substituted — see README.
SPECIES_WITHOUT_APPRIS = {"macaca_fascicularis": 9541}


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def _download(url: str, dest: Path) -> dict:
    last_modified = ""
    try:
        last_modified = requests.head(url, timeout=60).headers.get("Last-Modified", "")
    except Exception:
        pass
    if not dest.exists():
        print(f"  fetching {url}")
        with requests.get(url, stream=True, timeout=600) as resp:
            resp.raise_for_status()
            with dest.open("wb") as out:
                for block in resp.iter_content(1 << 20):
                    out.write(block)
    else:
        print(f"  cached   {dest.name}")
    return {
        "source_url": url,
        "source_last_modified": last_modified,
        "local_path": str(dest),
        "size_bytes": dest.stat().st_size,
        "sha256": _sha256(dest),
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cache-dir", default="~/.cache/deliverome/appris")
    ap.add_argument(
        "--publish",
        action="store_true",
        help="also upload the APPRIS tables to external://appris/e115v50/",
    )
    args = ap.parse_args()

    cache = Path(args.cache_dir).expanduser()
    cache.mkdir(parents=True, exist_ok=True)

    records = []
    for species, fname in APPRIS_FILES.items():
        url = f"{APPRIS_BASE}/{species}/{APPRIS_RELEASE}/{fname}"
        dest = cache / f"appris_{species}_{APPRIS_RELEASE}.txt"
        rec = _download(url, dest)
        rec |= {
            "dataset": "appris",
            "species": species,
            "release": APPRIS_RELEASE,
            "snapshot": APPRIS_SNAPSHOT,
        }
        records.append(rec)

    for species, fname in ENSEMBL_PEP.items():
        url = f"https://ftp.ensembl.org/pub/release-{ENSEMBL_RELEASE}/fasta/{species}/pep/{fname}"
        rec = _download(url, cache / fname)
        rec |= {
            "dataset": "ensembl_pep",
            "species": species,
            "release": ENSEMBL_RELEASE,
        }
        records.append(rec)

    prov = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "appris_release": APPRIS_RELEASE,
        "appris_snapshot": APPRIS_SNAPSHOT,
        "ensembl_release": ENSEMBL_RELEASE,
        "species_without_appris": SPECIES_WITHOUT_APPRIS,
        "records": records,
    }
    out = cache / "provenance.json"
    out.write_text(json.dumps(prov, indent=2) + "\n")
    print(f"wrote {out}")

    if args.publish:
        from deliverome_analysis.core import s3

        for rec in records:
            if rec["dataset"] != "appris":
                continue  # the Ensembl FASTAs are large and already versioned upstream
            key = f"appris/{APPRIS_RELEASE}/{Path(rec['local_path']).name}"
            s3.external.upload(rec["local_path"], key)
            print(f"  published external://{key}")
        s3.external.upload(str(out), f"appris/{APPRIS_RELEASE}/provenance.json")
        print(f"  published external://appris/{APPRIS_RELEASE}/provenance.json")


if __name__ == "__main__":
    main()
