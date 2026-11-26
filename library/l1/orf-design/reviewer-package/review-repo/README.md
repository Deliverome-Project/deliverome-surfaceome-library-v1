# Deliverome surfaceome library v1

The Deliverome surfaceome library v1 (L1 in the code) is a pooled, barcoded lentiviral library of
human cell-surface proteins, each expressed from its own open reading frame (ORF) with an
extracellular epitope tag: SpyTag003 and ALFA in one 51-amino-acid cassette. It is built to screen
which surface proteins carry a bound binder into the cell, and how far, in one pooled experiment:
every cell carries one construct, the tag reports display and uptake, and a 20-nucleotide barcode
says which protein the cell expressed.

This repository holds the library's design for external review, and the code and frozen inputs
to rerun every design decision offline.

## Where to start

| Read | What it gives you |
|---|---|
| [Design brief](library/l1/orf-design/reviewer-package/BRIEF.md) | the design as one argument: which genes, which protein form of each, where the tag goes and why, the construct and vector, barcodes and readout, and the questions we would like reviewers to weigh in on (section 12) |
| [Cover note](library/l1/orf-design/reviewer-package/COVER_NOTE.md) | what is still being worked on before the order, with links to the open issues |
| [Library layout map](library/l1/orf-design/l1_example_AXL_ires_puro.gb) | one complete construct (AXL as the example ORF) in the library vector, EFS – ORF – IRES2 – PuroR, annotated element by element (GenBank) |
| [Alternative layout map](library/l1/orf-design/l1_example_AXL_sv40_puro.gb) | the same insert with PuroR from its own SV40 promoter, the layout we did not use (GenBank) |
| [Decision log](library/l1/orf-design/reviewer-package/DECISION_LOG.md) | every design decision, dated, with who made it and why |

Comments are welcome as issues on this repository, or as comments on the brief.

## Rerun the design

You need [uv](https://docs.astral.sh/uv/) and the frozen input bundle, attached to this
repository's release
[`bundle-2026-11-26`](https://github.com/Deliverome-Project/deliverome-surfaceome-library-v1/releases/tag/bundle-2026-11-26)
as two tarballs matching [`bundle_manifest.tsv`](library/l1/orf-design/reviewer-package/bundle_manifest.tsv),
which indexes **3,570 files**: **3,286** `read_by_replay=yes` files [3.00 GB uncompressed]
and **284** `read_by_replay=no` supplement files [1.96 GB uncompressed].
For `bundle-2026-11-26`: `l1-reviewer-bundle-2026-11-26.tar.gz` holds all 3,286 replay files
(1.45 GB compressed, SHA-256 `281b15519be77a1732a3b9770c342cf2cbc991f1231f381aaefdd63f77c84c90`)
and `l1-pinned-releases-2026-11-26.tar` holds the 284 supplement files for inspection only
(1.96 GB, SHA-256 `cd3685894209a9d559091b64450729e9de9f3de65c91266305c85f54b004cc8f`).
After the download, no credentials or network access are used.

**Where the data lives.** Every bundle file has one listed home in
[`data_locations.tsv`](library/l1/orf-design/reviewer-package/data_locations.tsv),
and a test fails if any file of the bundle is missing from it:

| Home | What | How to get it |
|---|---|---|
| public bucket / upstream source | the public releases whose terms allow redistribution: Ensembl 115/116 and Compara, the AlphaFold DB models, UniProt, APPRIS (CC BY-NC-ND 4.0 unmodified, non-commercial), Pfam, HPA v23 and the accessible-surfaceome release (**2,932** files on public R2 [2.86 GB] + **1** Pfam-A file [0.66 GB] at EBI's identical copy = **2,933** public files, 3.52 GB; each row names its license) | anyone, no account: <https://pub-e53f37f6fe6b4745829c2066350fe26f.r2.dev/> + the `location` column (Pfam-A, 658 MB, at EBI's identical copy instead), or `uv run --locked python public_data.py fetch --out ../l1-reviewer-bundle` (checks every SHA-256) |
| this repository | small tables (≤ 5 MB), **210** files (0.09 GB) | already here, under `data/<bucket>/<key>` |
| release | everything else (**427** files, 1.34 GB): our own predictions and records, and public sources whose redistribution terms are unclear or restricted (PeptideAtlas, OpenProt, ARCHS4, the Tedman et al. tables, GTEx v10 summary TPMs, and private `accessible-surfaceome` repo snapshots) | the release tarballs above, which also hold the complete bundle |

The release tarballs alone are enough to rerun everything; the public bucket and `data/` are
for reading or citing individual files without them.

```bash
gh release download bundle-2026-11-26 --repo Deliverome-Project/deliverome-surfaceome-library-v1 --dir ..
```

```bash
mkdir -p ../l1-reviewer-bundle && tar -xf ../l1-reviewer-bundle-2026-11-26.tar.gz -C ../l1-reviewer-bundle
```

To inspect the full releases as well, unpack the second tarball into the same folder:

```bash
tar -xf ../l1-pinned-releases-2026-11-26.tar -C ../l1-reviewer-bundle
```

```bash
uv sync --locked --no-group flow --no-group flow-gui --no-group modal
```

```bash
cd library/l1/orf-design/reviewer-package
uv run --locked python run_review.py --bundle ../../../../../l1-reviewer-bundle
```

The run checks every unpacked bundle file against `bundle_manifest.tsv`, reruns all **30** pipeline
stages in `library_pipeline.py` offline and reports whether each reproduces its published record,
then recomputes every number in the brief. It takes about three minutes. The package README beside
the brief describes each stage and the bundle's contents.

## Pinned releases

Every release the design rests on is in the bundle whole. Each file's release, version,
source URL and SHA-256 are in `bundle_manifest.tsv` beside the brief.

| Release | Version | Files | Size | In the bundle | Source |
|---|---|---:|---:|---|---|
| ARCHS4 (cell-line panel consensus) | 2.5 | 3 | 9 MB | `s3/external/archs4/v2.5/cell-line-panel-consensus/2026-09-04/` | https://archs4.org/ |
| Human Protein Atlas | 23.0 | 1 | 206 MB | `s3/external/human-protein-atlas/2023-05/` | https://v23.proteinatlas.org/about/download |
| Human Protein Atlas | 23.0 | 2 | 206 MB | `s3/external/hpa/v23/` | https://v23.proteinatlas.org/about/download |
| Ensembl | 115 | 18 | 491 MB | `s3/external/ensembl/release-115/` | https://ftp.ensembl.org/pub/release-115/ |
| Ensembl | 116 | 14 | 531 MB | `s3/external/ensembl/release-116/` | https://ftp.ensembl.org/pub/release-116/ |
| APPRIS principal isoforms | e115v50 (snapshot 2026_06.v50) | 3 | 8 MB | `s3/external/appris/e115v50/` | https://apprisws.bioinfo.cnio.es/pub/releases/2026_06.v50/ |
| APPRIS full annotations, incl. TRIFID scores | 2026_06.v50 (e115v50) | 1 | 41 MB | `s3/external/appris/2026_06.v50/` | https://apprisws.bioinfo.cnio.es/pub/releases/2026_06.v50/ |
| APPRIS principal isoforms, mouse | 2026_08.v50 (e116v50) | 2 | 6 MB | `s3/external/appris/2026_08.v50/` | https://apprisws.bioinfo.cnio.es/pub/releases/2026_08.v50/ |
| UniProtKB | 2026_03 (downloaded 2026-10-11) | 4 | 41 MB | `s3/external/uniprot/2026-10-11/` | https://rest.uniprot.org/ |
| UniProtKB REST queries | 2026_03 | 2 | 761 kB | `s3/external/uniprot/2026_03/` | https://rest.uniprot.org/ |
| UniProtKB REST queries | 2026_04 | 40 | 9 MB | `s3/external/uniprot/2026_04/` | https://rest.uniprot.org/ |
| OpenProt | 2.2 | 1 | 231 MB | `s3/external/openprot/2.2/` | https://www.openprot.org/ |
| PeptideAtlas human build, replayed queries | build 607 | 10 | 33 MB | `s3/external/peptideatlas/build-607/` | https://peptideatlas.org/ |
| Pfam | 39.0 | 4 | 659 MB | `s3/external/pfam/39.0/` | https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam39.0/ |
| accessible-surfaceome catalog | 1.3.0 | 3 | 221 MB | `s3/external/accessible-surfaceome/release-1.3.0/` | https://doi.org/10.5281/zenodo.23445521 |
| accessible-surfaceome repository files | commit f9c73691 | 10 | 10 MB | `s3/external/accessible-surfaceome/f9c73691/` | https://github.com/Deliverome-Project/accessible-surfaceome (private) |
| DeepTMHMM topology snapshots (Deliverome runs) | 14 dated snapshots (`2026-10-14` through `2026-11-23-pre-freeze`) | 14 | 154 MB | `s3/external/accessible-surfaceome/topology/` | external://accessible-surfaceome/topology/ |
| Expression TPM, GTEx v10 / HPA / TCGA (Deliverome-compiled) | deliverome-internal-5a1416ea, rederived-2026-08-04 | 3 | 7 MB | `s3/external/expression-tpm/` | external://expression-tpm/ |
| Tedman et al. 2026 GPCR surface screen, supplementary table | 2025-09-19 download | 2 | 909 kB | `s3/external/tedman-gpcr-surface-screen/` | https://doi.org/10.1038/s41467-026-76564-7 |
| Compara 116 ortholog CDS-eligibility draw (final-v2) and its verification | 2026-11-05 | 123 | 123 MB | `s3/processed/sprint4/s4e9-orf-design-fg-panel/library-design/ortholog-cds-eligibility/2026-11-05/` | processed://sprint4/s4e9-orf-design-fg-panel/library-design/ortholog-cds-eligibility/2026-11-05/ |
| DeepTMHMM v1 and NetGPI 1.1 predictions for the ortholog draw | 2026-11-05 | 43 | 21 MB | `s3/processed/sprint4/s4e9-orf-design-fg-panel/ortholog-cds-policy-predictions-2026-11-05/` | processed://sprint4/s4e9-orf-design-fg-panel/ortholog-cds-policy-predictions-2026-11-05/ |
| Ortholog rebuild input snapshot | 2026-11-05 | 4 | 15 MB | `s3/processed/sprint4/s4e9-orf-design-fg-panel/library-design/ortholog-rebuild-inputs/2026-11-05/` | processed://sprint4/s4e9-orf-design-fg-panel/library-design/ortholog-rebuild-inputs/2026-11-05/ |
| Ortholog correspondence review | 2026-11-07 | 6 | 534 kB | `s3/processed/sprint4/s4e9-orf-design-fg-panel/ortholog-correspondence-review/2026-11-07/` | processed://sprint4/s4e9-orf-design-fg-panel/ortholog-correspondence-review/2026-11-07/ |

## Layout

The folders keep the paths they have in the lab's analysis repository, so the scripts and
the brief's links work unchanged.

| Path | What it is |
|---|---|
| `library/l1/orf-design/` | the library pipeline, its committed tables, and `reviewer-package/` (brief, runner, bundle tools) |
| `library/l1/surface-filters/` | the gene filter |
| `library/l1/proteoform-selection/` | proteoform selection, panels, CDS resolution and construct assembly |
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

## License and citation

- **Code:** [MIT License](LICENSE).
- **Derived library designs, construct tables, brief, and figures:** dedicated to the public domain under [CC0 1.0 Universal](https://creativecommons.org/publicdomain/zero/1.0/) — see [`LICENSING.md`](LICENSING.md) for full terms and upstream source licenses.
- **Citation:** see [`CITATION.cff`](CITATION.cff) (also available via GitHub's *Cite this repository* button).
