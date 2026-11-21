# L1 surface-protein ORF library: design review

This repository holds the design of Deliverome's L1 ORF library for external review, and
the code to rerun it. Start with the brief,
[`BRIEF.md`](analysis/sprint4/s4e9-orf-design-fg-panel/reviewer-package/BRIEF.md). The issues
still being worked on before the order are tracked in the
[cover note](analysis/sprint4/s4e9-orf-design-fg-panel/reviewer-package/COVER_NOTE.md).

## Rerun the design

You need [uv](https://docs.astral.sh/uv/) and the frozen input bundle, attached to this
repository's release
[`bundle-2026-11-21`](https://github.com/Deliverome-Project/l1-library-review/releases/tag/bundle-2026-11-21)
as two tarballs. `l1-reviewer-bundle-2026-11-21.tar.gz` holds the 3,123 files the rerun reads
(1.38 GB, SHA-256 `9b4356cdd6bab64a22e6dc4ae22c1451894ea0d10baafc44bb01bc84c32d1695`).
`l1-pinned-releases-2026-11-21.tar` holds the rest of every pinned release, for inspection only
(248 files, 1.65 GB, SHA-256 `c1d36285f157e24a8a1f3a7fe0825184adf9cf97f600372e317fde256133a37e`).
After the download, no credentials or network access are used.

**Where the data lives.** Every bundle file has one listed home in
[`data_locations.tsv`](analysis/sprint4/s4e9-orf-design-fg-panel/reviewer-package/data_locations.tsv),
and a test fails if any file of the bundle is missing from it:

| Home | What | How to get it |
|---|---|---|
| public bucket | the public releases whose terms allow redistribution: Ensembl 115/116 and Compara, the AlphaFold DB models, UniProt, Pfam, HPA v23 and the accessible-surfaceome release (2,775 files, 3.2 GB; each row names its license) | anyone, no account: <https://pub-e53f37f6fe6b4745829c2066350fe26f.r2.dev/> + the `location` column (Pfam-A, 658 MB, at EBI's identical copy instead), or `uv run --locked python public_data.py fetch --out ../l1-reviewer-bundle` (checks every SHA-256) |
| this repository | small tables (≤ 5 MB), 188 files | already here, under `data/<bucket>/<key>` |
| release | everything else: our own predictions and records, and public sources whose redistribution terms are unclear or non-commercial (PeptideAtlas, OpenProt, APPRIS, GENCODE Ribo-seq ORFs, the Tedman et al. tables) | the release tarballs above, which also hold the complete bundle |

The release tarballs alone are enough to rerun everything; the public bucket and `data/` are
for reading or citing individual files without them.

```bash
gh release download bundle-2026-11-21 --repo Deliverome-Project/l1-library-review --dir ..
```

```bash
mkdir ../l1-reviewer-bundle && tar -xf ../l1-reviewer-bundle-2026-11-21.tar.gz -C ../l1-reviewer-bundle
```

To inspect the full releases as well, unpack the second tarball into the same folder:

```bash
tar -xf ../l1-pinned-releases-2026-11-21.tar -C ../l1-reviewer-bundle
```

```bash
uv sync --locked
```

```bash
cd analysis/sprint4/s4e9-orf-design-fg-panel/reviewer-package
uv run --locked python run_review.py --bundle ../../../../../l1-reviewer-bundle
```

The run checks the bundle, reruns the library's 16 decision stages and reports whether each
reproduces its published record, then recomputes every number in the brief. It takes
about three minutes. The package README beside the brief describes each stage and the
bundle's contents.

## Pinned releases

Every release the design rests on is in the bundle whole. Each file's release, version,
source URL and SHA-256 are in `bundle_manifest.tsv` beside the brief.

| Release | Version | Pinned on S3 | Files | Size | In the bundle | Source |
|---|---|---|---:|---:|---|---|
| Ensembl | 115 | 2026-10-11 | 17 | 490 MB | `s3/external/ensembl/release-115/` | https://ftp.ensembl.org/pub/release-115/ |
| Ensembl | 116 | 2026-11-05 | 14 | 531 MB | `s3/external/ensembl/release-116/` | https://ftp.ensembl.org/pub/release-116/ |
| APPRIS principal isoforms | e115v50 (snapshot 2026_06.v50) | 2026-10-08 | 3 | 8 MB | `s3/external/appris/e115v50/` | https://apprisws.bioinfo.cnio.es/pub/releases/2026_06.v50/ |
| APPRIS full annotations, incl. TRIFID scores | 2026_06.v50 (e115v50) | 2026-11-02 | 1 | 41 MB | `s3/external/appris/2026_06.v50/` | https://apprisws.bioinfo.cnio.es/pub/releases/2026_06.v50/ |
| APPRIS principal isoforms, mouse | 2026_08.v50 (e116v50) | 2026-11-06 | 2 | 6 MB | `s3/external/appris/2026_08.v50/` | https://apprisws.bioinfo.cnio.es/pub/releases/2026_08.v50/ |
| UniProtKB | 2026_03 (downloaded 2026-10-11) | 2026-10-11 | 4 | 41 MB | `s3/external/uniprot/2026-10-11/` | https://rest.uniprot.org/ |
| UniProtKB REST queries | 2026_03 | 2026-11-02 | 1 | 759 kB | `s3/external/uniprot/2026_03/` | https://rest.uniprot.org/ |
| UniProtKB REST queries | 2026_04 | 2026-11-05 | 36 | 4 MB | `s3/external/uniprot/2026_04/` | https://rest.uniprot.org/ |
| OpenProt | 2.2 | 2026-11-02 | 1 | 231 MB | `s3/external/openprot/2.2/` | https://www.openprot.org/ |
| PeptideAtlas human build, replayed queries | build 607 | 2026-11-02 | 2 | 17 MB | `s3/external/peptideatlas/build-607/` | https://peptideatlas.org/ |
| Pfam | 39.0 | 2026-11-07 | 4 | 659 MB | `s3/external/pfam/39.0/` | https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam39.0/ |
| accessible-surfaceome catalog | 1.3.0 | 2026-11-02 | 3 | 221 MB | `s3/external/accessible-surfaceome/release-1.3.0/` | https://doi.org/10.5281/zenodo.23445521 |
| accessible-surfaceome repository files | commit f9c73691 | 2026-11-02 | 10 | 10 MB | `s3/external/accessible-surfaceome/f9c73691/` | https://github.com/Deliverome-Project/accessible-surfaceome (private) |
| DeepTMHMM topology snapshots (Deliverome runs) | 2026-10-14, 2026-10-14r2, 2026-11-01, 2026-11-05-l1-compara116, 2026-11-06-ortholog-cds-policy | 2026-10-15 | 5 | 47 MB | `s3/external/accessible-surfaceome/topology/` | external://accessible-surfaceome/topology/ |
| Expression TPM, GTEx v10 / HPA / TCGA (Deliverome-compiled) | deliverome-internal-5a1416ea, rederived-2026-08-04 | 2026-11-02 | 3 | 7 MB | `s3/external/expression-tpm/` | external://expression-tpm/ |
| Tedman et al. 2026 GPCR surface screen, supplementary table | 2025-09-19 download | 2026-08-14 | 2 | 909 kB | `s3/external/tedman-gpcr-surface-screen/` | https://doi.org/10.1038/s41467-026-76564-7 |
| Compara 116 ortholog CDS-eligibility draw (final-v2) and its verification | 2026-11-05 | 2026-11-06 | 123 | 123 MB | `s3/processed/sprint4/s4e9-orf-design-fg-panel/library-design/ortholog-cds-eligibility/2026-11-05/` | processed://sprint4/s4e9-orf-design-fg-panel/library-design/ortholog-cds-eligibility/2026-11-05/ |
| DeepTMHMM v1 and NetGPI 1.1 predictions for the ortholog draw | 2026-11-05 | 2026-11-06 | 43 | 21 MB | `s3/processed/sprint4/s4e9-orf-design-fg-panel/ortholog-cds-policy-predictions-2026-11-05/` | processed://sprint4/s4e9-orf-design-fg-panel/ortholog-cds-policy-predictions-2026-11-05/ |
| Ortholog rebuild input snapshot | 2026-11-05 | 2026-11-05 | 4 | 15 MB | `s3/processed/sprint4/s4e9-orf-design-fg-panel/library-design/ortholog-rebuild-inputs/2026-11-05/` | processed://sprint4/s4e9-orf-design-fg-panel/library-design/ortholog-rebuild-inputs/2026-11-05/ |
| Ortholog correspondence review | 2026-11-07 | 2026-11-07 | 6 | 534 kB | `s3/processed/sprint4/s4e9-orf-design-fg-panel/ortholog-correspondence-review/2026-11-07/` | processed://sprint4/s4e9-orf-design-fg-panel/ortholog-correspondence-review/2026-11-07/ |

## Layout

The folders keep the paths they have in the lab's analysis repository, so the scripts and
the brief's links work unchanged.

| Path | What it is |
|---|---|
| `analysis/sprint4/s4e9-orf-design-fg-panel/` | the library pipeline, its committed tables, and `reviewer-package/` (brief, runner, bundle tools) |
| `analysis/sprint1/s1e3-surface-library-filters/` | the gene filter |
| `scratch/becca/appris-principal/` | proteoform selection, panels, CDS resolution and construct assembly |
| `src/deliverome_analysis/` | the shared utilities this code imports (S3 access and the offline mirror, sequence tools, topology, plotting); the rest of the lab's package, unrelated to L1, is not exported |
| `historical/` | script versions two stages rerun as they were when their records were made |
| `EXPORT_MANIFEST.json` | which analysis-repo commit each managed file came from, with its SHA-256 |

## Adding to this repository

Add files freely. This repository is exported from Deliverome-Project/deliverome-analysis
(`export_review_repo.py`), and a re-export only updates the files listed in
`EXPORT_MANIFEST.json`. If you edit one of those, the next export leaves your version in
place and lists it, so the change can be carried back to the analysis repository by hand.
Files you add are never touched.

New utilities go in `src/deliverome_analysis/`, with dependencies declared in
`pyproject.toml` and locked in `uv.lock` (`uv add <package>`), never installed ad hoc.
