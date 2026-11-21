"""Per-protein publication counts for the surface-proteomics challenge set, from pinned inputs.

The challenge set (section 2 of the reviewer brief) is built from the cross-method survey's
published snapshots on S3. ``surface_proteomics_inputs.json`` pins, for every snapshot the brief
has cited, the S3 key and SHA-256 of each file the counting reads; ``current`` names the
snapshot the brief uses. This script fetches those keys, refuses any whose bytes differ from the
pin, and runs the survey's own counting code from ``scratch/becca/s3e8-library-gap/`` unchanged:

1. ``abundance_roles.py`` (abundance snapshots only) overlays the abundance screen's review roles
   onto the sweep's detections, which the abundance snapshot copies through unchanged.
2. ``surface_specific.py`` counts, per protein, the surface-labeling, membrane-fractionation and
   whole-lysate glycoproteome publications that detect it, and joins library membership on
   Ensembl gene via ``gene_resolver`` against ``candidate_manifest_tranches.tsv`` (its SHA-256 is
   pinned too, so a changed library manifest stops the run).

It writes ``surface_proteomics_counts.tsv``: one row per (snapshot, protein) for proteins with at
least 20 surface publications (no protein below that can pass rule B or any variant), with the
number of publications of each kind in that snapshot as ``N_*`` columns. The rule and its
variants are applied downstream by ``surface_proteomics_logodds_compare.py`` and
``surface_proteomics_gap_candidates.py``.

    uv run --locked python analysis/sprint4/s4e9-orf-design-fg-panel/surface_proteomics_snapshot.py
    # re-pin after deliberately adding a snapshot to the record (writes the SHA-256s):
    uv run --locked python analysis/sprint4/s4e9-orf-design-fg-panel/surface_proteomics_snapshot.py --pin
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
GAP = REPO / "scratch/becca/s3e8-library-gap"
INPUTS = HERE / "surface_proteomics_inputs.json"
COUNTS = HERE / "surface_proteomics_counts.tsv"
MIN_SURFACE = 20

#: S3 file name -> path the survey scripts expect inside a sweep folder.
SWEEP_LAYOUT = {
    "detections_long.tsv": "detections_long.tsv",
    "protein_method_matrix.tsv": "protein_method_matrix.tsv",
    "novel_candidates.tsv": "novel_candidates.tsv",
    "clustering__candidate_review.tsv": "clustering/candidate_review.tsv",
}
ABUNDANCE_FILES = ("detections_long.tsv", "manifest.json")

csv.field_size_limit(10**9)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def fetch(entry: dict, dest: Path, pin: bool) -> Path:
    from deliverome_analysis import s3

    bucket, key = entry["key"].split("://", 1)
    dest.parent.mkdir(parents=True, exist_ok=True)
    getattr(s3, bucket).download(key, dest)
    got = sha256(dest)
    if pin:
        entry["sha256"] = got
    elif got != entry["sha256"]:
        raise SystemExit(
            f"{entry['key']}: SHA-256 {got[:12]} is not the pinned {entry['sha256'][:12]}"
        )
    return dest


def publications_per_kind(detections: Path) -> dict[str, int]:
    """Exactly rule_variants.py's classification of publications into the three kinds."""
    pubs = {"surface": set(), "membrane": set(), "global": set()}
    with open(detections, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["method_group"] == "membrane-fractionation":
                pubs["membrane"].add(r["publication_id"])
            elif r["review_role"] == "complementary":
                pubs["global"].add(r["publication_id"])
            elif r["admitted"] == "True":
                pubs["surface"].add(r["publication_id"])
    return {k: len(v) for k, v in pubs.items()}


def run(script: str, *args: Path) -> None:
    proc = subprocess.run(
        [sys.executable, script, *map(str, args)], cwd=GAP, capture_output=True, text=True
    )
    if proc.returncode:
        raise SystemExit(f"{script} failed:\n{proc.stderr[-2000:]}")


def build(name: str, snap: dict, library: Path, work: Path, pin: bool) -> pd.DataFrame:
    sweep = work / name / "sweep"
    for fname, rel in SWEEP_LAYOUT.items():
        fetch(snap["sweep"][fname], sweep / rel, pin)
    data = sweep
    if "abundance" in snap:
        abundance = work / name / "abundance"
        for fname in ABUNDANCE_FILES:
            fetch(snap["abundance"][fname], abundance / fname, pin)
        data = work / name / "overlaid"
        run("abundance_roles.py", sweep, abundance, data)
    d = work / name / "counts"
    d.mkdir(parents=True)
    shutil.copy(library, d / "manifest.tsv")
    run("surface_specific.py", d, data)
    n = publications_per_kind(data / "detections_long.tsv")
    df = pd.read_csv(d / "surface_counts_all_proteins.tsv", sep="\t", keep_default_na=False)
    df = df[df.n_surface_pubs >= MIN_SURFACE].copy()
    df.insert(0, "snapshot", name)
    for kind, count in n.items():
        df[f"N_{kind}"] = count
    print(
        f"{name}: publications {n}; {len(df)} proteins with >= {MIN_SURFACE} surface publications"
    )
    return df[
        [
            "snapshot",
            "uniprot_accession",
            "gene_symbol",
            "catalog_tier",
            "n_surface_pubs",
            "n_membrane_pubs",
            "n_global_pubs",
            "N_surface",
            "N_membrane",
            "N_global",
            "library",
        ]
    ]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pin", action="store_true", help="record SHA-256s instead of checking them")
    args = ap.parse_args()
    spec = json.loads(INPUTS.read_text())
    library = REPO / spec["library_manifest"]["path"]
    got = sha256(library)
    if args.pin:
        spec["library_manifest"]["sha256"] = got
    elif got != spec["library_manifest"]["sha256"]:
        raise SystemExit(f"{library.name} has changed since the pin ({got[:12]})")
    work = Path(tempfile.mkdtemp(prefix="surface_proteomics_"))
    frames = [
        build(name, snap, library, work, args.pin) for name, snap in spec["snapshots"].items()
    ]
    pd.concat(frames).to_csv(COUNTS, sep="\t", index=False)
    if args.pin:
        INPUTS.write_text(json.dumps(spec, indent=2) + "\n")
        print(f"pinned {INPUTS.name}")
    print(f"wrote {COUNTS.name}; current snapshot: {spec['current']}")


if __name__ == "__main__":
    main()
