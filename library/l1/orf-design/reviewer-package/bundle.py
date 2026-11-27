"""Build the frozen input bundle that ``run_review.py`` replays. Lab-side: needs S3 credentials.

The bundle is not a hand-picked list. ``collect`` reruns every ``library_pipeline.py``
stage (and ``brief_numbers.py``) with the S3 helper in record mode
(``DELIVEROME_S3_MIRROR_MODE=record``) and a fresh home directory, so exactly the objects
the pipeline reads are copied, and no warm developer cache can hide one. It then adds every
pinned release the design rests on (``PINNED`` below), whole, so a reviewer can inspect a
release and not only the objects the replay happens to read; each file pinned by SHA-256 in
an input manifest (``PIN_MANIFESTS``) is checked against that hash. ``audit`` fails if any
pinned input is missing from the bundle manifest. ``manifest``
hashes each file and labels where it comes from (``SOURCES`` below), which is what decides
what can be shared outside the lab. ``pack`` writes a tarball; ``stage`` uploads it to this
experiment's own ``processed://`` prefix and nowhere else.

    uv run --locked python bundle.py collect --out ~/l1-bundle        # ~4 min, ~800 MB
    uv run --locked python bundle.py manifest --bundle ~/l1-bundle    # -> bundle_manifest.tsv
    uv run --locked python bundle.py audit                            # pinned but not bundled
    uv run --locked python bundle.py pack --bundle ~/l1-bundle --out l1-reviewer-bundle.tar
    uv run --locked python bundle.py pack ... --part replay       # or supplement (2 GiB asset cap)
    uv run --locked python bundle.py pack ... --exclude-class deliverome_curation
    uv run --locked python bundle.py stage --tar l1-reviewer-bundle.tar --date 2026-11-02
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).parent
S4E9 = HERE.parent
MANIFEST = HERE / "bundle_manifest.tsv"
STAGE_PREFIX = "sprint4/s4e9-orf-design-fg-panel/reviewer-package/"
SNAP = r"sprint4/s4e9-orf-design-fg-panel/library-design/(pipeline|panel|universe)-inputs/[^/]+/"
REPO = S4E9.parents[2]
E = "sprint4/s4e9-orf-design-fg-panel/"
ASSET_LIMIT = 2 * 1024**3  # GitHub release assets

#: Every pinned release the L1 design rests on, bundled whole:
#: (bucket, key prefix, release, version, source URL).
PINNED: list[tuple[str, str, str, str, str]] = [
    (
        "external",
        "archs4/v2.5/cell-line-panel-consensus/2026-09-04/",
        "ARCHS4 (cell-line panel consensus)",
        "2.5",
        "https://archs4.org/",
    ),
    (
        "external",
        "human-protein-atlas/2023-05/",
        "Human Protein Atlas",
        "23.0",
        "https://v23.proteinatlas.org/about/download",
    ),
    (
        "external",
        "hpa/v23/",
        "Human Protein Atlas",
        "23.0",
        "https://v23.proteinatlas.org/about/download",
    ),
    (
        "external",
        "ensembl/release-115/",
        "Ensembl",
        "115",
        "https://ftp.ensembl.org/pub/release-115/",
    ),
    (
        "external",
        "ensembl/release-116/",
        "Ensembl",
        "116",
        "https://ftp.ensembl.org/pub/release-116/",
    ),
    (
        "external",
        "appris/e115v50/",
        "APPRIS principal isoforms",
        "e115v50 (snapshot 2026_06.v50)",
        "https://apprisws.bioinfo.cnio.es/pub/releases/2026_06.v50/",
    ),
    (
        "external",
        "appris/2026_06.v50/",
        "APPRIS full annotations, incl. TRIFID scores",
        "2026_06.v50 (e115v50)",
        "https://apprisws.bioinfo.cnio.es/pub/releases/2026_06.v50/",
    ),
    (
        "external",
        "appris/2026_08.v50/",
        "APPRIS principal isoforms, mouse",
        "2026_08.v50 (e116v50)",
        "https://apprisws.bioinfo.cnio.es/pub/releases/2026_08.v50/",
    ),
    (
        "external",
        "uniprot/2026-10-11/",
        "UniProtKB",
        "2026_03 (downloaded 2026-10-11)",
        "https://rest.uniprot.org/",
    ),
    (
        "external",
        "uniprot/2026_03/",
        "UniProtKB REST queries",
        "2026_03",
        "https://rest.uniprot.org/",
    ),
    (
        "external",
        "uniprot/2026_04/",
        "UniProtKB REST queries",
        "2026_04",
        "https://rest.uniprot.org/",
    ),
    ("external", "openprot/2.2/", "OpenProt", "2.2", "https://www.openprot.org/"),
    (
        "external",
        "peptideatlas/build-607/",
        "PeptideAtlas human build, replayed queries",
        "build 607",
        "https://peptideatlas.org/",
    ),
    (
        "external",
        "pfam/39.0/",
        "Pfam",
        "39.0",
        "https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam39.0/",
    ),
    (
        "external",
        "accessible-surfaceome/release-1.3.0/",
        "accessible-surfaceome catalog",
        "1.3.0",
        "https://doi.org/10.5281/zenodo.23445521",
    ),
    (
        "external",
        "accessible-surfaceome/f9c73691/",
        "accessible-surfaceome repository files",
        "commit f9c73691",
        "https://github.com/Deliverome-Project/accessible-surfaceome (private)",
    ),
    (
        "external",
        "accessible-surfaceome/topology/",
        "DeepTMHMM topology snapshots (Deliverome runs)",
        "see key date",
        "external://accessible-surfaceome/topology/",
    ),
    (
        "external",
        "expression-tpm/",
        "Expression TPM, GTEx v10 / HPA / TCGA (Deliverome-compiled)",
        "see key",
        "external://expression-tpm/",
    ),
    (
        "external",
        "tedman-gpcr-surface-screen/",
        "Tedman et al. 2026 GPCR surface screen, supplementary table",
        "2025-09-19 download",
        "https://doi.org/10.1038/s41467-026-76564-7",
    ),
    (
        "processed",
        E + "library-design/ortholog-cds-eligibility/2026-11-05/",
        "Compara 116 ortholog CDS-eligibility draw (final-v2) and its verification",
        "2026-11-05",
        "processed://" + E + "library-design/ortholog-cds-eligibility/2026-11-05/",
    ),
    (
        "processed",
        E + "ortholog-cds-policy-predictions-2026-11-05/",
        "DeepTMHMM v1 and NetGPI 1.1 predictions for the ortholog draw",
        "2026-11-05",
        "processed://" + E + "ortholog-cds-policy-predictions-2026-11-05/",
    ),
    (
        "processed",
        E + "library-design/ortholog-rebuild-inputs/2026-11-05/",
        "Ortholog rebuild input snapshot",
        "2026-11-05",
        "processed://" + E + "library-design/ortholog-rebuild-inputs/2026-11-05/",
    ),
    (
        "processed",
        E + "ortholog-correspondence-review/2026-11-07/",
        "Ortholog correspondence review",
        "2026-11-07",
        "processed://" + E + "ortholog-correspondence-review/2026-11-07/",
    ),
]
#: Input manifests whose SHA-256 pins every bundled copy must match. Paths from the repo root.
PIN_MANIFESTS = [
    "library/l1/orf-design/library_pipeline_inputs.json",
    "library/l1/orf-design/reproduce_isoform_universe_inputs.json",
    "library/l1/orf-design/ortholog_rebuild_inputs.json",
    "library/l1/orf-design/ortholog_cds_eligibility_v2_inputs.json",
    "library/l1/orf-design/ortholog_cds_policy_final_run.json",
    "library/l1/orf-design/ortholog_accession_replay_inputs.json",
    "library/l1/proteoform-selection/panel_inputs.json",
    "library/l1/proteoform-selection/candidate_inputs_provenance.json",
]

#: Where each bundled object comes from, first match wins. The class decides sharing:
#:   public_release      byte-identical to a named public release; anyone can fetch it
#:   public_derived      public data, but queried or assembled by Deliverome
#:   surfaceome_repo     a file of the accessible-surfaceome repository at a pinned commit
#:   deliverome_curation an internal input: curation calls, internal tables, predictions
#:   deliverome_record   a published library record that the replay is checked against
SOURCES: list[tuple[str, str, str]] = [
    (r"external/ensembl/release-11[56]/", "public_release", "Ensembl FTP"),
    (r"external/pfam/39\.0/", "public_release", "Pfam 39.0 (EBI release copy)"),
    (r"external/uniprot/2026_04/", "public_derived", "UniProt 2026_04 REST query"),
    (
        r"external/archs4/v2\.5/",
        "public_derived",
        "ARCHS4 v2.5 consensus over the S2E18 panel lines",
    ),
    (
        r"external/human-protein-atlas/2023-05/",
        "public_release",
        "HPA v23 (same bytes as hpa/v23/)",
    ),
    (
        r"external/topdb/v3\.3/",
        "public_release",
        "TOPDB 3.3 full download (topology_vs_experiment.py)",
    ),
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
    (
        r"processed/sprint4/s4e9-orf-design-fg-panel/(library-design/ortholog-rebuild-inputs/"
        r"|ortholog-cds-policy-predictions-)",
        "deliverome_curation",
        "ortholog draw inputs and Deliverome predictions (DeepTMHMM v1, NetGPI 1.1)",
    ),
    (
        r"processed/sprint4/s4e9-orf-design-fg-panel/ortholog-correspondence-review/",
        "deliverome_record",
        "ortholog correspondence review (PR 273)",
    ),
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
        r"external/accessible-surfaceome/d1/",
        "deliverome_curation",
        "surfaceome D1 tables saved whole by S5E13 d1_snapshot.py (dated, append-only)",
    ),
    (
        r"external/expression-tpm/",
        "deliverome_curation",
        "expression TPM tables (GTEx v10 / HPA / TCGA, Deliverome-compiled)",
    ),
    (
        r"external/hpa/v23/",
        "public_release",
        "Human Protein Atlas 23.0 cell-line RNA (rna_celline.tsv), the source of the co-receptor medians",
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
        r"processed/scratch/becca/tag-site-corpus-audit/",
        "deliverome_curation",
        "tag-site literature audit: other files in the audit folder the brief cites",
    ),
    (
        r"processed/sprint5/s5e13-tag-site-benchmarking/corpus_audit_corpus\.tsv",
        "deliverome_curation",
        "S5E13 training corpus with the audit's corrections (PR 232)",
    ),
    (
        r"processed/sprint4/s4e9-orf-design-fg-panel/synthesis/",
        "deliverome_record",
        "recoded synthesis order file of record, cited by the brief",
    ),
    (
        r"processed/sprint4/s4e9-orf-design-fg-panel/type3-nterm/",
        "deliverome_curation",
        "type III N-terminal spacer test, raw DeepTMHMM v1 calls, cited by the brief",
    ),
    # Objects the 2026-11-14 inputs added (AFDB archive, saved UniProt snapshots, rung builds).
    (
        r"external/alphafold-db/v6/",
        "public_release",
        "AlphaFold DB model_v6 files, archived write-once",
    ),
    (
        r"external/gencode-riboseq-orfs/",
        "public_release",
        "GENCODE Ribo-seq ORFs, fetched 2026-11-01",
    ),
    (r"external/mygene/", "public_derived", "MyGene.info responses, 2026-11-02"),
    (
        r"external/uniprot/rest-2026-11-[12][0-9]/",
        "public_derived",
        "UniProt REST snapshots, C-terminal review",
    ),
    (
        r"processed/sprint4/s4e9-orf-design-fg-panel/reviewer-package/nterm/uniprot-",
        "public_derived",
        "UniProt 2026_04 feature queries (N-terminal review)",
    ),
    (
        r"processed/sprint4/s4e9-orf-design-fg-panel/pinned-inputs/afdb-v6/",
        "public_derived",
        "AlphaFold DB model_v6 files, tarred by Deliverome",
    ),
    (
        r"processed/sprint4/s4e9-orf-design-fg-panel/pinned-inputs/rung-builds/",
        "deliverome_record",
        "rung arm tables at each cut's pinned commit",
    ),
    (
        r"processed/sprint4/s4e9-orf-design-fg-panel/(predicted-structures|membrane-side-gate|"
        r"reprediction|reviewer-package/nterm|triage)/",
        "deliverome_curation",
        "Deliverome predictions and gate inputs (ColabFold, SignalP re-prediction, gate bundle)",
    ),
    (
        r"processed/sprint5/s5e13-tag-site-benchmarking/(l1_snorkel_designs|v5/)",
        "deliverome_curation",
        "S5E13 snorkel designs and v5 training table",
    ),
]


def classify(rel: str) -> tuple[str, str]:
    for pattern, cls, origin in SOURCES:
        if re.match(pattern, rel):
            return cls, origin
    raise SystemExit(f"no SOURCES rule for {rel}; add one before sharing the bundle")


def sha256(path: Path) -> str:
    from deliverome_analysis.core.provenance import sha256_file

    return sha256_file(path)


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
    return rc.returncode or rc2.returncode or add_pinned(out)


def _norm(loc: str) -> tuple[str, str] | None:
    m = re.match(r"(?:s3://deliverome-)?(external|processed|raw)(?::/?/?|/)(\S+)", loc.strip())
    return (m.group(1), m.group(2)) if m else None


def manifest_pins() -> dict[tuple[str, str], dict]:
    """Every (bucket, key) an input manifest pins, with its SHA-256 where it gives one."""
    pins: dict[tuple[str, str], dict] = {}

    def add(loc: str, sha: str, src: str) -> None:
        for part in loc.split(","):
            bk = _norm(part.split(" (")[0])
            if bk:
                p = pins.setdefault(bk, {"sha256": "", "pinned_by": set()})
                p["pinned_by"].add(src)
                p["sha256"] = p["sha256"] or sha

    def walk(o, src: str, prefix: str | None) -> None:
        if isinstance(o, dict):
            prefix = o.get("prefix", prefix) if isinstance(o.get("prefix"), str) else prefix
            for name, v in o.items():
                if isinstance(v, dict) and "sha256" in v:
                    loc = next(
                        (
                            v[k]
                            for k in ("source", "key", "s3_key", "uri")
                            if isinstance(v.get(k), str)
                        ),
                        None,
                    )
                    if loc is None and prefix and "/" not in name:
                        loc = f"processed:{prefix}{name}"
                    if loc:
                        add(loc, v.get("sha256") or "", src)
                if name not in ("origin", "note", "notes", "description"):
                    walk(v, src, prefix)
        elif isinstance(o, list):
            for v in o:
                walk(v, src, prefix)
        elif isinstance(o, str) and re.match(r"(external|processed)://", o) and "..." not in o:
            add(o, "", src)

    for rel in PIN_MANIFESTS:
        walk(json.loads((REPO / rel).read_text()), Path(rel).name, None)
    return pins


#: Documents a reviewer reads; every S3 file they cite is bundled, so it opens offline.
CITED_DOCS = ("DESIGN.md", "BRIEF.md")
#: An S3 location cited in prose, e.g. ``processed://sprint5/.../corpus_audit_corpus.tsv``.
S3_CITATION = re.compile(r"\b(raw|processed|external)://([^\s`'\")<>|,]+)")
#: The same location written as its path inside the bundle, ``s3/<bucket>/<key>``.
BUNDLE_CITATION = re.compile(r"(?<![\w/])s3/(raw|processed|external)/([^\s`'\")<>|,]+)")


def cited_locations(text: str) -> set[tuple[str, str]]:
    """(bucket, key) for each S3 location ``text`` cites; a key ending in ``/`` is a folder.

    Skipped: patterns and elisions (``…``, ``{``, ``<``), and this package's own staging
    prefix, which names the bundle tarballs rather than files in them.
    """
    out = set()
    for bucket, key in S3_CITATION.findall(text) + BUNDLE_CITATION.findall(text):
        key = key.rstrip(".;:")
        if "…" in key or "{" in key or "<" in key or key.startswith(STAGE_PREFIX):
            continue
        out.add((bucket, key))
    return out


def doc_citations() -> set[tuple[str, str]]:
    return set().union(*(cited_locations((HERE / d).read_text()) for d in CITED_DOCS))


def release_of(bucket: str, key: str) -> tuple[str, str, str]:
    for b, prefix, release, version, url in PINNED:
        if b == bucket and key.startswith(prefix):
            if version.startswith("see key"):
                version = key[len(prefix) :].split("/")[0]
            return release, version, url
    if re.match(SNAP, key):
        return "Deliverome pipeline-input snapshot", key.split("/")[3], f"processed://{key}"
    return "Deliverome library record", "", f"{bucket}://{key}"


def add_pinned(out: Path) -> int:
    """Copy every pinned release and manifest-pinned file into the mirror; check SHA-256 pins."""
    from deliverome_analysis.core import s3

    mirror = out / "s3"
    # What the replay read is recorded once, right after collect's replay; a later add-pinned
    # (to add a newly pinned release) must not count the pinned files as read by the replay.
    keys = out / "replay_keys.txt"
    if keys.exists():
        replay = keys.read_text().split()
    else:
        replay = sorted(p.relative_to(mirror).as_posix() for p in mirror.rglob("*") if p.is_file())
        keys.write_text("\n".join(replay) + "\n")
    wanted: dict[tuple[str, str], str] = {}
    for b, prefix, *_ in PINNED:
        for key in getattr(s3, b).list_keys(prefix=prefix):
            wanted[(b, key)] = ""
    for bk, pin in manifest_pins().items():
        wanted[bk] = pin["sha256"]
    # Files the brief and cover note cite but no stage reads (the S5E13 corpus, the order
    # FASTA): a reviewer working offline from the bundle must still be able to open them.
    for b, key in sorted(doc_citations()):
        keys = getattr(s3, b).list_keys(prefix=key) if key.endswith("/") else [key]
        for k in keys:
            wanted.setdefault((b, k), "")

    def _fetch_one(item: tuple[tuple[str, str], str]) -> str | None:
        (b, key), sha = item
        dest = mirror / b / key
        if not dest.is_file():
            bucket = getattr(s3, b)
            if not bucket.exists(key):
                return f"{b}://{key}: pinned but absent from S3"
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(bucket.download(key), dest)
        if sha and sha256(dest) != sha:
            return f"{b}://{key}: SHA-256 differs from its pin"
        return None

    items = sorted(wanted.items())
    with ThreadPoolExecutor(max_workers=8) as pool:
        bad = [err for err in pool.map(_fetch_one, items) if err is not None]
    for line in bad:
        print(line, file=sys.stderr)
    print(f"pinned objects: {len(wanted)}; replay read {len(replay)}")
    return 1 if bad else 0


def active_record_pins() -> dict[tuple[str, str], dict]:
    """Active code-pinned S3 records read by ``library_pipeline.py`` and ``brief_numbers.py``."""
    if str(S4E9) not in sys.path:
        sys.path.insert(0, str(S4E9))
    import library_pipeline as lp
    import synthesis_compliance as sc

    from deliverome_analysis.cloning import library_design as ld
    from deliverome_analysis.surfaceome import topology

    sc_summary = json.loads((S4E9 / "synthesis_compliance_summary.json").read_text())
    return {
        ("external", topology.SNAPSHOT_KEY): {
            "sha256": topology.SNAPSHOT_SHA256,
            "pinned_by": {"topology.SNAPSHOT_KEY"},
        },
        ("processed", ld.CDS_KEY): {
            "sha256": "",
            "pinned_by": {"library_design.CDS_KEY"},
        },
        ("processed", ld.CONSTRUCTS_KEY): {
            "sha256": "",
            "pinned_by": {"library_design.CONSTRUCTS_KEY"},
        },
        ("processed", sc.FASTA_KEY): {
            "sha256": sc_summary["recoded_fasta_sha256"],
            "pinned_by": {"synthesis_compliance.FASTA_KEY"},
        },
        ("processed", lp.PREFIX + lp.REFERENCE_RECORD + "reference_selection.tsv"): {
            "sha256": "",
            "pinned_by": {"library_pipeline.REFERENCE_RECORD"},
        },
        ("processed", lp.PREFIX + lp.REFERENCE_RECORD + "reference_disagreements.tsv"): {
            "sha256": "",
            "pinned_by": {"library_pipeline.REFERENCE_RECORD"},
        },
    }


def audit() -> int:
    """Pinned inputs (manifests, active records, and PINNED releases) absent from ``bundle_manifest.tsv``."""
    with MANIFEST.open() as fh:
        have = {(r["bucket"], r["key"]): r for r in csv.DictReader(fh, delimiter="\t")}
    pins = manifest_pins()
    for bk, info in active_record_pins().items():
        p = pins.setdefault(bk, {"sha256": "", "pinned_by": set()})
        p["pinned_by"].update(info["pinned_by"])
        p["sha256"] = p["sha256"] or info["sha256"]
    missing, wrong = [], []
    for bk, pin in sorted(pins.items()):
        if bk not in have:
            missing.append(f"{bk[0]}://{bk[1]}  (pinned by {', '.join(sorted(pin['pinned_by']))})")
        elif pin["sha256"] and have[bk]["sha256"] != pin["sha256"]:
            wrong.append(f"{bk[0]}://{bk[1]}")
    prefixes = {(b, p) for b, p, *_ in PINNED}
    empty = [
        f"{b}://{p}"
        for b, p in sorted(prefixes)
        if not any(k[0] == b and k[1].startswith(p) for k in have)
    ]
    for label, rows in (
        ("pinned but not bundled", missing),
        ("SHA-256 differs", wrong),
        ("release prefix with no bundled file", empty),
    ):
        print(f"{label}: {len(rows)}")
        for r in rows:
            print(f"  {r}")
    return 1 if missing or wrong or empty else 0


def releases(out: Path) -> int:
    """Write the "Pinned releases" table (Markdown) from bundle_manifest.tsv.

    One row per release prefix in ``PINNED``: release, version, the date the copy was pinned
    on S3 (earliest upload of its files), file count and size, and where it sits in the bundle.
    """
    from deliverome_analysis.core import s3

    with MANIFEST.open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    lines = [
        "| Release | Version | Pinned on S3 | Files | Size | In the bundle | Source |",
        "|---|---|---|---:|---:|---|---|",
    ]
    for b, prefix, release, version, url in PINNED:
        got = [r for r in rows if r["bucket"] == b and r["key"].startswith(prefix)]
        if not got:
            raise SystemExit(f"{b}://{prefix} has no bundled file; run audit")
        dates = sorted(
            getattr(s3, b)._head(r["key"])["LastModified"].date().isoformat() for r in got
        )
        if version.startswith("see key"):
            version = ", ".join(sorted({r["key"][len(prefix) :].split("/")[0] for r in got}))
        mb = sum(int(r["bytes"]) for r in got) / 1e6
        size = f"{mb:,.0f} MB" if mb >= 1 else f"{mb * 1000:,.0f} kB"
        lines.append(
            f"| {release} | {version} | {dates[0]} | {len(got)} | {size} | `s3/{b}/{prefix}` | {url} |"
        )
    out.write_text("\n".join(lines) + "\n")
    print(f"wrote {out}")
    return 0


def manifest(bundle: Path) -> int:
    root = bundle / "s3"
    replay_file = bundle / "replay_keys.txt"
    replay = set(replay_file.read_text().split()) if replay_file.exists() else set()
    paths = sorted(p for p in root.rglob("*") if p.is_file())
    with ThreadPoolExecutor(max_workers=min(8, os.cpu_count() or 4)) as pool:
        digests = list(pool.map(sha256, paths))
    rows = []
    for path, digest in zip(paths, digests, strict=True):
        rel = path.relative_to(root).as_posix()
        bucket, _, key = rel.partition("/")
        cls, origin = classify(rel)
        release, version, url = release_of(bucket, key)
        rows.append(
            {
                "bucket": bucket,
                "key": key,
                "bytes": path.stat().st_size,
                "sha256": digest,
                "class": cls,
                "origin": origin,
                "release": release,
                "version": version,
                "source_url": url,
                "read_by_replay": "yes" if rel in replay else "no",
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


def pack(bundle: Path, out: Path, exclude: list[str], part: str = "all") -> int:
    """Write a tarball. ``part`` splits it for GitHub's 2 GiB asset limit: ``replay`` holds
    what ``run_review.py`` reads, ``supplement`` the rest of every pinned release."""
    keep = {"all": ("yes", "no"), "replay": ("yes",), "supplement": ("no",)}[part]
    with MANIFEST.open() as fh:
        rows = [
            r
            for r in csv.DictReader(fh, delimiter="\t")
            if r["class"] not in exclude and r.get("read_by_replay", "yes") in keep
        ]
    # ``.tar.gz`` compresses: the replay part holds ~2,700 AlphaFold DB PDB text files, which
    # left it at 2.76 GB uncompressed on 2026-11-14, over the asset limit. ``tar -xf`` reads both.
    # Use parallel pigz when available, otherwise gzip level 6 (~3x faster than default level 9).
    if out.name.endswith(".tar.gz") and shutil.which("pigz"):
        with out.open("wb") as fh:
            proc = subprocess.Popen(["pigz", "-c"], stdin=subprocess.PIPE, stdout=fh)
            assert proc.stdin is not None
            with tarfile.open(fileobj=proc.stdin, mode="w|") as tar:
                tar.add(MANIFEST, arcname="bundle_manifest.tsv")
                for r in rows:
                    tar.add(
                        bundle / "s3" / r["bucket"] / r["key"],
                        arcname=f"s3/{r['bucket']}/{r['key']}",
                    )
            proc.stdin.close()
            if proc.wait() != 0:
                raise SystemExit(f"pigz failed with exit code {proc.returncode}")
    else:
        open_kw = {"compresslevel": 6} if out.name.endswith(".tar.gz") else {}
        mode = "w:gz" if out.name.endswith(".tar.gz") else "w"
        with tarfile.open(out, mode, **open_kw) as tar:
            tar.add(MANIFEST, arcname="bundle_manifest.tsv")
            for r in rows:
                tar.add(
                    bundle / "s3" / r["bucket"] / r["key"],
                    arcname=f"s3/{r['bucket']}/{r['key']}",
                )
    size = out.stat().st_size
    print(f"wrote {out}: {len(rows)} files, {size / 1e9:.2f} GB"
          + (f" (excluded {exclude})" if exclude else ""))  # fmt: skip
    if size >= ASSET_LIMIT:
        raise SystemExit(f"{out.name} is over GitHub's 2 GiB release-asset limit; write .tar.gz")
    return 0


def stage(tar: Path, date: str) -> int:
    """Upload the tarball under this experiment's prefix only; never anywhere public."""
    from deliverome_analysis.core.s3 import processed

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
    sub.add_parser("add-pinned").add_argument("--out", type=Path, required=True)
    sub.add_parser("audit")
    sub.add_parser("releases").add_argument("--out", type=Path, default=HERE / "pinned_releases.md")
    p = sub.add_parser("pack")
    p.add_argument("--bundle", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--exclude-class", action="append", default=[])
    p.add_argument("--part", choices=("all", "replay", "supplement"), default="all")
    s = sub.add_parser("stage")
    s.add_argument("--tar", type=Path, required=True)
    s.add_argument("--date", required=True)
    a = ap.parse_args(argv)
    if a.cmd == "collect":
        return collect(a.out.expanduser())
    if a.cmd == "manifest":
        return manifest(a.bundle.expanduser())
    if a.cmd == "add-pinned":
        return add_pinned(a.out.expanduser())
    if a.cmd == "audit":
        return audit()
    if a.cmd == "releases":
        return releases(a.out)
    if a.cmd == "pack":
        return pack(a.bundle.expanduser(), a.out, a.exclude_class, a.part)
    return stage(a.tar, a.date)


if __name__ == "__main__":
    raise SystemExit(main())
