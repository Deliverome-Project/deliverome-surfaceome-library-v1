"""Build the frozen input bundle that ``run_review.py`` replays. Lab-side: needs S3 credentials.

The bundle is not a hand-picked list. ``collect`` reruns every ``library_pipeline.py``
stage (and ``brief_numbers.py``) with the S3 helper in record mode
(``DELIVEROME_S3_MIRROR_MODE=record``) and a fresh home directory, so exactly the objects
the pipeline reads are copied, and no warm developer cache can hide one. ``manifest``
hashes each file and labels where it comes from (``SOURCES`` below), which is what decides
what can be shared outside the lab. ``pack`` writes a tarball; ``stage`` uploads it to this
experiment's own ``processed://`` prefix and nowhere else.

    uv run --locked python bundle.py collect --out ~/l1-bundle        # ~4 min, ~800 MB
    uv run --locked python bundle.py manifest --bundle ~/l1-bundle    # -> bundle_manifest.tsv
    uv run --locked python bundle.py pack --bundle ~/l1-bundle --out l1-reviewer-bundle.tar
    uv run --locked python bundle.py pack ... --exclude-class deliverome_curation
    uv run --locked python bundle.py stage --tar l1-reviewer-bundle.tar --date 2026-11-02
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
import re
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path

HERE = Path(__file__).parent
S4E9 = HERE.parent
MANIFEST = HERE / "bundle_manifest.tsv"
STAGE_PREFIX = "sprint4/s4e9-orf-design-fg-panel/reviewer-package/"
SNAP = r"sprint4/s4e9-orf-design-fg-panel/library-design/(pipeline|panel|universe)-inputs/[^/]+/"

#: Where each bundled object comes from, first match wins. The class decides sharing:
#:   public_release      byte-identical to a named public release; anyone can fetch it
#:   public_derived      public data, but queried or assembled by Deliverome
#:   surfaceome_repo     a file of the accessible-surfaceome repository at a pinned commit
#:   deliverome_curation an internal input: curation calls, internal tables, predictions
#:   deliverome_record   a published library record that the replay is checked against
SOURCES: list[tuple[str, str, str]] = [
    (r"external/ensembl/release-115/", "public_release", "Ensembl release 115 FTP"),
    (r"external/uniprot/2026-10-11/", "public_release", "UniProt 2026_03 table download"),
    (r"external/uniprot/2026_03/", "public_derived", "UniProt 2026_03 REST query"),
    (r"external/appris/", "public_release", "APPRIS e115v50 / 2026_06.v50"),
    (r"external/openprot/2.2/", "public_release", "OpenProt 2.2"),
    (
        r"external/accessible-surfaceome/release-1.3.0/",
        "public_release",
        "accessible-surfaceome 1.3.0, Zenodo doi:10.5281/zenodo.23445521",
    ),
    (
        r"external/tedman-gpcr-surface-screen/",
        "public_release",
        "Tedman 2025 GPCR screen supplementary table",
    ),
    (r"external/peptideatlas/", "public_derived", "PeptideAtlas build 607, replayed queries"),
    (rf"processed/{SNAP}Ribo-seq_ORFs\.", "public_release", "GENCODE Ribo-seq ORFs"),
    (rf"processed/{SNAP}peptideatlas_", "public_derived", "PeptideAtlas build 607 queries"),
    (rf"processed/{SNAP}s1e3/mygene_", "public_derived", "MyGene.info responses, 2026-11-02"),
    (rf"processed/{SNAP}uniprot_human_gpi", "public_derived", "UniProt KW-0336 query"),
    (
        r"external/accessible-surfaceome/f9c73691/",
        "surfaceome_repo",
        "accessible-surfaceome commit f9c73691",
    ),
    (
        rf"processed/{SNAP}(deeptmhmm_|human_canonical_non_hla|ensg_to_uniprot)",
        "surfaceome_repo",
        "accessible-surfaceome commit f9c73691 (snapshot copy)",
    ),
    (
        r"external/accessible-surfaceome/topology/",
        "deliverome_curation",
        "DeepTMHMM 1.0.24 predictions run by Deliverome (topology snapshot)",
    ),
    (
        r"external/expression-tpm/",
        "deliverome_curation",
        "expression TPM tables (GTEx v10 / HPA / TCGA, Deliverome-compiled)",
    ),
    (
        r"processed/sprint4/s4e9-orf-design-fg-panel/p0039_",
        "deliverome_curation",
        "p0039 backbone map",
    ),
    (rf"processed/{SNAP}", "deliverome_curation", "Deliverome pipeline-input snapshot"),
    (
        r"processed/sprint4/s4e9-orf-design-fg-panel/library-(design|candidates)/",
        "deliverome_record",
        "published L1 library record",
    ),
    # Reviewer records the brief cites; brief_numbers.py reads them, so collect records them.
    (
        r"processed/scratch/becca/tag-site-corpus-audit/audit\.tsv",
        "deliverome_curation",
        "audit of the 176 S5E13 literature tag-site positives against their papers (PR 236)",
    ),
    (
        r"processed/scratch/becca/tag-site-corpus-audit/gold_standard\.tsv",
        "deliverome_curation",
        "gold standard: the audit's confirmed sites, one per (accession, residue) (PR 236)",
    ),
    (
        r"processed/scratch/becca/tag-site-corpus-audit/(failed_sites|added_sites|region_labels)\.tsv",
        "deliverome_curation",
        "tag-site audit: failed sites, sites recovered from wrong-gene papers, region labels (PR 236)",
    ),
    (
        r"processed/sprint5/s5e13-tag-site-benchmarking/corpus_audit_corpus\.tsv",
        "deliverome_curation",
        "S5E13 training corpus with the audit's corrections (PR 232)",
    ),
]


def classify(rel: str) -> tuple[str, str]:
    for pattern, cls, origin in SOURCES:
        if re.match(pattern, rel):
            return cls, origin
    raise SystemExit(f"no SOURCES rule for {rel}; add one before sharing the bundle")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def collect(out: Path) -> int:
    mirror = out / "s3"
    mirror.mkdir(parents=True, exist_ok=True)
    env = {
        **os.environ,
        "DELIVEROME_S3_MIRROR": str(mirror),
        "DELIVEROME_S3_MIRROR_MODE": "record",
        "HOME": tempfile.mkdtemp(prefix="l1-bundle-home-"),
        "PYTHONHASHSEED": "0",
    }
    # The fresh HOME hides developer caches; it must not also hide the AWS profile that
    # the credentials come from (AWS_PROFILE would then name a profile that does not exist).
    for var, name in (
        ("AWS_CONFIG_FILE", "config"),
        ("AWS_SHARED_CREDENTIALS_FILE", "credentials"),
    ):
        real = Path.home() / ".aws" / name
        if var not in env and real.is_file():
            env[var] = str(real)
    rc = subprocess.run([sys.executable, "library_pipeline.py", "verify"], cwd=S4E9, env=env)
    work = Path(tempfile.mkdtemp())
    rc2 = subprocess.run(
        [sys.executable, "brief_numbers.py", "--no-figures", "--out", str(work / "n.json")],
        cwd=HERE,
        env=env,
    )
    return rc.returncode or rc2.returncode


def manifest(bundle: Path) -> int:
    root = bundle / "s3"
    rows = []
    for path in sorted(p for p in root.rglob("*") if p.is_file()):
        rel = path.relative_to(root).as_posix()
        bucket, _, key = rel.partition("/")
        cls, origin = classify(rel)
        rows.append(
            {
                "bucket": bucket,
                "key": key,
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
                "class": cls,
                "origin": origin,
            }
        )
    if not rows:
        raise SystemExit(f"{root} holds no files; refusing to overwrite {MANIFEST.name}")
    with MANIFEST.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    by: dict[str, list[int]] = {}
    for r in rows:
        by.setdefault(r["class"], []).append(r["bytes"])
    for cls, sizes in sorted(by.items()):
        print(f"{cls:<20} {len(sizes):>3} files {sum(sizes) / 1e6:>8.1f} MB")
    print(f"wrote {MANIFEST.name}: {len(rows)} files")
    return 0


def pack(bundle: Path, out: Path, exclude: list[str]) -> int:
    with MANIFEST.open() as fh:
        rows = [r for r in csv.DictReader(fh, delimiter="\t") if r["class"] not in exclude]
    with tarfile.open(out, "w") as tar:
        tar.add(MANIFEST, arcname="bundle_manifest.tsv")
        for r in rows:
            tar.add(bundle / "s3" / r["bucket"] / r["key"], arcname=f"s3/{r['bucket']}/{r['key']}")
    print(f"wrote {out}: {len(rows)} files" + (f" (excluded {exclude})" if exclude else ""))
    return 0


def stage(tar: Path, date: str) -> int:
    """Upload the tarball under this experiment's prefix only; never anywhere public."""
    from deliverome_analysis.s3 import processed

    key = f"{STAGE_PREFIX}{date}/{tar.name}"
    if processed.exists(key):
        raise SystemExit(f"processed://{key} exists; stage under a new date")
    processed.upload(tar, key)
    print(f"staged processed://{key}")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("collect").add_argument("--out", type=Path, required=True)
    sub.add_parser("manifest").add_argument("--bundle", type=Path, required=True)
    p = sub.add_parser("pack")
    p.add_argument("--bundle", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--exclude-class", action="append", default=[])
    s = sub.add_parser("stage")
    s.add_argument("--tar", type=Path, required=True)
    s.add_argument("--date", required=True)
    a = ap.parse_args(argv)
    if a.cmd == "collect":
        return collect(a.out.expanduser())
    if a.cmd == "manifest":
        return manifest(a.bundle.expanduser())
    if a.cmd == "pack":
        return pack(a.bundle.expanduser(), a.out, a.exclude_class)
    return stage(a.tar, a.date)


if __name__ == "__main__":
    raise SystemExit(main())
