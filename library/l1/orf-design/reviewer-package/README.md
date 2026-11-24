# L1 library reviewer package

A self-contained way to rerun the L1 ORF library's design decisions without Deliverome
credentials. Read [`BRIEF.md`](BRIEF.md) first: it is the design argument, and every number in it is
checked by this package.

N-terminal tag placement behind a signal peptide has its own companion package and brief,
in [`nterm/`](nterm/README.md).

## What you need

| Item | Notes |
|---|---|
| This repository, with its git history | a clone, or a `git bundle` of the branch. History is needed only for two stages (below). |
| The frozen input bundle | a folder (or tarball) holding `s3/external/...` and `s3/processed/...`: 98 files the replay reads, ~1.14 GB. A second tarball adds the rest of every pinned release (238 files, ~1.40 GB); the replay does not need it |
| [uv](https://docs.astral.sh/uv/) | installs the locked Python environment from `uv.lock` |
| ~10 GB free disk, ~5 minutes | the whole run is offline |

No network access, AWS account or `.env` file is used. The runner points the repository's
S3 helper at the bundle (`DELIVEROME_S3_MIRROR`) and sends any other S3 call to a closed
local port, so a stage that tried to read live data would fail rather than pass.

## Run it

```bash
uv sync --locked
```

```bash
cd library/l1/orf-design/reviewer-package
uv run --locked python run_review.py --bundle /path/to/l1-reviewer-bundle
```

It first checks every bundle file against [`bundle_manifest.tsv`](bundle_manifest.tsv) (SHA-256), then
reruns each of the 23 pipeline stages and prints what each decides and whether it reproduces
its published record, then recomputes the brief's numbers and compares them with
[`brief_numbers.json`](brief_numbers.json). A summary is written to `review_report.json`. One stage:

```bash
uv run --locked python run_review.py --bundle /path/to/bundle --stage reference
```

## What each stage reproduces

| # | Stage | Decides | Checked against |
|---:|---|---|---|
| 1 | gene_list | which genes: the cohort filter over surfaceome release 1.3.0 | 2,370 pass; minus the recorded EPPIN exclusion = the 2,369 library genes |
| 2 | gpi_set | which forms are GPI-anchored | `gpi_set.tsv`, byte for byte |
| 3 | reference | each gene's reference proteoform (the ladder) | `reference_selection.tsv`, byte for byte |
| 4 | isoform_universe | every alternate form and its evidence tier | summaries byte for byte; cohort table row for row (see the hash-seed limit in the brief) |
| 5 | ortholog_repick | each gene's mouse and cyno proteoform | byte for byte |
| 5b | universe_coverage | how many library genes the isoform and ortholog pools could draw from (Compara 115) | `universe_coverage.tsv` and its summary, byte for byte |
| 6 | panels | the 150-isoform and 200-ortholog draws | byte for byte |
| 7 | rescue_arms | 25 expression-rescue arms | byte for byte |
| 8 | giants | truncations of oversized genes (L2) | byte for byte |
| 9 | candidate_cds | each member's designed CDS | all 7,017 recorded CDS |
| 10 | constructs | domesticated ORFs and vectors | all 7,017 ORFs and 7,066 vectors |
| 11 | tranches | L1/L2 membership, duplicates retired | no two active members share a protein |
| 12 | topology | DeepTMHMM v1 per designed protein | every active member has one |
| 13 | tag_paths | each member's placement rung | `l1_tag_paths.tsv`, `l2_tag_paths.tsv`, byte for byte |
| 14 | cassettes | each C-terminal member's linker encoding | every issued batch replays |

Stages 4 and 5 rerun scripts at the commits that produced their records, so they read the
repository's git history; without it the runner reports them as skipped, never as passed.

DeepTMHMM and NetGPI are not rerun: their predictions are inputs, pinned as committed
outputs and in the bundle's topology snapshot. Stage 12 checks that every member has a
prediction for its exact sequence, not that the predictor would return the same call.

## What is in the bundle, and what is public

Each bundle file is labeled in `bundle_manifest.tsv` with its sharing class, its release,
version and source URL, and whether the replay reads it (`read_by_replay`). Every release the
design rests on is bundled whole, not only the files the replay opens, so a reviewer can
inspect it; [`pinned_releases.md`](pinned_releases.md) lists them (release, version, date
pinned, size, place in the bundle), and `bundle.py audit` fails if anything an input manifest
pins is missing.

| Class | Read by the replay | Rest of the pinned releases | Meaning |
|---|---:|---:|---|
| `public_release` | 26 | 31 | byte-identical to a named public release (Ensembl 115 and 116 including the Compara homology dumps, UniProt, APPRIS, OpenProt 2.2, Pfam 39.0, surfaceome 1.3.0 on Zenodo, GENCODE Ribo-seq ORFs, the published GPCR screen table) |
| `public_derived` | 7 | 36 | public data, as queried by us (PeptideAtlas build 607, MyGene.info, UniProt keyword and accession queries) |
| `surfaceome_repo` | 15 | 1 | files of the accessible-surfaceome repository at commit f9c73691 |
| `deliverome_curation` | 35 | 48 | internal inputs: curation calls, the compiled expression table, the vector backbone map, DeepTMHMM and NetGPI predictions, the ortholog draw inputs, and the audited S5E13 tag-site positives the brief cites |
| `deliverome_record` | 15 | 122 | published library records: what the replay is compared with, the Compara 116 ortholog draw and its verification, and the ortholog correspondence review |

Reviewers receive the full bundle, all five classes (Becca, 2026-11-04). A bundle can still
be packed without a class (`bundle.py pack --exclude-class ...`); the stages
that read the missing files then fail with a message naming the absent key.

## Building the bundle (lab side)

```bash
uv run --locked python bundle.py collect --out ~/l1-bundle
```

```bash
uv run --locked python bundle.py manifest --bundle ~/l1-bundle
```

```bash
uv run --locked python bundle.py audit
```

```bash
uv run --locked python bundle.py pack --bundle ~/l1-bundle --out l1-reviewer-bundle.tar --part replay
```

```bash
uv run --locked python bundle.py pack --bundle ~/l1-bundle --out l1-pinned-releases.tar --part supplement
```

`collect` reruns every stage with S3 in record mode and a fresh home directory, so the
bundle holds everything the pipeline reads, then adds every pinned release whole and checks
each file an input manifest pins against its SHA-256. The two parts are packed separately to
stay under GitHub's 2 GiB limit per release asset. `stage` uploads a tarball under this
experiment's `processed://` prefix; nothing in the package publishes anywhere else.

## The review repository

Reviewers get the code from a separate private repository,
[Deliverome-Project/deliverome-surfaceome-library-v1](https://github.com/Deliverome-Project/deliverome-surfaceome-library-v1),
generated from this one by `export_review_repo.py`. It contains:

- the library's folders whole (s4e9, the S1E3 gene filter, `library/l1/proteoform-selection`,
  `library/l1/tag-site-architecture`);
- every other file the replay opens (listed in `review_repo_files.txt`, found by tracing a
  replay);
- the complete `deliverome_analysis` package, `pyproject.toml` and `uv.lock`;
- in `historical/`, the old script versions that stages 4 and 5 rerun.

The input bundle is attached to that repository's release
[`bundle-2026-11-20`](https://github.com/Deliverome-Project/deliverome-surfaceome-library-v1/releases/tag/bundle-2026-11-20)
as two tarballs: `l1-reviewer-bundle-2026-11-20.tar.gz` (the 3,120 files the replay reads, 1.36 GB,
SHA-256 `8bdb0006b739f28a860017aef7dd480d24e50bb2c2fa6d68120276a00a3be165`; gzipped to stay under
GitHub's 2 GiB asset limit) and `l1-pinned-releases-2026-11-20.tar` (the rest of every pinned
release and every S3 file the brief cites, 246 files, 1.65 GB, SHA-256
`56bd5618245f7ef0ae7e7ab85ee2445c31f39ec501af771ae6287c68adb73d13`). The same files are staged at
`processed://sprint4/s4e9-orf-design-fg-panel/reviewer-package/2026-11-20/`. `tar -xf` reads both.
A rebuilt bundle goes on a new release with a new date, never over an existing one; earlier
releases are `bundle-2026-11-02` (before l1-p2740..2742), `bundle-2026-11-05` (before the
audited tag-site positives were added), `bundle-2026-11-05-r2` (before the stop-codon fix of 13
designed CDS and the extended audit), `bundle-2026-11-05-r3` (before the tag-site gold
standard), `bundle-2026-11-05-r4` (before the final audit refresh), `bundle-2026-11-05-r5`
(before the universe-coverage stage and its Compara dumps) and `bundle-2026-11-05-r6` (before
the Compara 116 ortholog panel was adopted and every pinned release was bundled whole), and
`bundle-2026-11-10` (before the 2026-11-14 order of record and the files the brief cites).

`EXPORT_MANIFEST.json` there records the source commit and a SHA-256 per managed file. To
update it after changes here:

```bash
uv run --locked python export_review_repo.py --dest ~/git/deliverome-surfaceome-library-v1
```

A re-export never touches files added in the review repo. A managed file that was edited there is left in place and listed,
so the edit can be carried back here by hand. If a new script is needed by the replay,
rebuild the list from a traced run (`--trace`) before exporting.
