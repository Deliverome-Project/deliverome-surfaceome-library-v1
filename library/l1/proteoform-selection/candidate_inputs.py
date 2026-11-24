"""Every external input this pipeline needs, pinned to an S3 key.

The scripts used to read a developer's ``~/.cache/deliverome`` directory, which
meant a clean clone could not reproduce the run at all. Each input is now named
by the exact S3 key it comes from and downloaded on demand into a local cache,
so the only thing that varies between machines is whether the cache is already
warm.

Pins are deliberate. Bump a version here and the whole pipeline moves together.
"""

from __future__ import annotations

from pathlib import Path

from deliverome_analysis import s3

CACHE = Path("~/.cache/deliverome/library-candidates").expanduser()

ENSEMBL_RELEASE = "115"
UNIPROT_SNAPSHOT = "2026-10-11"

ASSEMBLIES = {
    "human": "Homo_sapiens.GRCh38",
    "mouse": "Mus_musculus.GRCm39",
    "cyno": "Macaca_fascicularis.Macaca_fascicularis_6.0",
}
UNIPROT_FILES = {
    "human": ["uniprot_human_reviewed.tsv.gz", "uniprot_human_reviewed_isoforms.tsv.gz"],
    "mouse": ["uniprot_mouse_reviewed_isoforms.tsv.gz"],
    "cyno": ["uniprot_cyno_all_isoforms.tsv.gz"],
}
BACKBONE_KEY = "sprint4/s4e9-orf-design-fg-panel/p0039_EFS_IRES_puro_backbone.gb"
REFERENCE_VECTOR_KEY = "sprint4/s4e9-orf-design-fg-panel/p0039_AXL_vector.gb"


def _fetch(bucket, key: str) -> Path:
    dest = CACHE / key.replace("/", "__")
    if not dest.exists():
        dest.parent.mkdir(parents=True, exist_ok=True)
        print(f"    fetching {key}")
        bucket.download(key, dest)
    return dest


def cds_fasta(species: str) -> Path:
    stem = ASSEMBLIES[species]
    return _fetch(s3.external, f"ensembl/release-{ENSEMBL_RELEASE}/cds/{stem}.cds.all.fa.gz")


def gtf(species: str = "human") -> Path:
    """The pinned Ensembl GTF, for the canonical-transcript tag (``Ensembl_canonical``)."""
    stem = ASSEMBLIES[species]
    return _fetch(
        s3.external, f"ensembl/release-{ENSEMBL_RELEASE}/gtf/{stem}.{ENSEMBL_RELEASE}.gtf.gz"
    )


def pep_fasta(species: str) -> Path:
    stem = ASSEMBLIES[species]
    return _fetch(s3.external, f"ensembl/release-{ENSEMBL_RELEASE}/pep/{stem}.pep.all.fa.gz")


def uniprot_tables(species: str) -> list[Path]:
    return [_fetch(s3.external, f"uniprot/{UNIPROT_SNAPSHOT}/{n}")
            for n in UNIPROT_FILES.get(species, [])]


def backbone() -> Path:
    return _fetch(s3.processed, BACKBONE_KEY)


def reference_vector() -> Path:
    """p0039_AXL: the vector whose Esp3I wrapper grammar this build copies."""
    return _fetch(s3.processed, REFERENCE_VECTOR_KEY)


def provenance() -> dict[str, str]:
    """What a run should record so the inputs can be identified later."""
    out = {
        "ensembl_release": ENSEMBL_RELEASE,
        "uniprot_snapshot": UNIPROT_SNAPSHOT,
        "backbone": f"processed://{BACKBONE_KEY}",
        "reference_vector": f"processed://{REFERENCE_VECTOR_KEY}",
    }
    for sp, stem in ASSEMBLIES.items():
        out[f"cds_{sp}"] = f"external://ensembl/release-{ENSEMBL_RELEASE}/cds/{stem}.cds.all.fa.gz"
        out[f"pep_{sp}"] = f"external://ensembl/release-{ENSEMBL_RELEASE}/pep/{stem}.pep.all.fa.gz"
    for sp, names in UNIPROT_FILES.items():
        out[f"uniprot_{sp}"] = ", ".join(
            f"external://uniprot/{UNIPROT_SNAPSHOT}/{n}" for n in names)
    return out
