# L1 library reviewer package

A self-contained way to rerun the L1 ORF library's design decisions without Deliverome
credentials. Read [`BRIEF.md`](BRIEF.md) first: it is the design argument, and every number in it is
checked by this package.

N-terminal tag placement behind a signal peptide has its own companion package and brief,
in [`nterm/`](nterm/README.md).

## What you need

| Item | Notes |
|---|---|
| This repository, with its git history (or an exported review repo with `historical/`) | a clone, or a `git bundle` of the branch. Exported review repos ship `historical/` tarballs and `historical/git-snapshots/` so all 30 stages run with zero git history. |
| The frozen input bundle | a folder (or tarballs) holding `s3/external/...` and `s3/processed/...`: [`bundle_manifest.tsv`](bundle_manifest.tsv) currently indexes **3,591** files (**3,310** `read_by_replay=yes`, ~3.04 GB uncompressed; plus **281** `read_by_replay=no` supplement files, ~1.96 GB uncompressed, for inspection only) |
| [uv](https://docs.astral.sh/uv/) | installs the locked Python environment from `uv.lock` |
| ~10 GB free disk, ~5 minutes | the whole run is offline |

No network access, AWS account or `.env` file is used. The runner points the repository's
S3 helper at the bundle (`DELIVEROME_S3_MIRROR`) and sends any other S3 or HTTP(S) call to a closed
local port (`127.0.0.1:9`), so a stage that tried to read live data fails loudly rather than passing
(for kernel-level isolation, run via `docker run --network=none`).

## Run it

```bash
uv sync --locked --no-group flow --no-group flow-gui --no-group modal
```

```bash
cd library/l1/orf-design/reviewer-package
uv run --locked python run_review.py --bundle /path/to/l1-reviewer-bundle
```

It first checks every unpacked bundle file against [`bundle_manifest.tsv`](bundle_manifest.tsv) (SHA-256), then
reruns all **30** pipeline stages in `library_pipeline.py` and prints what each decides and whether it reproduces
its published record, then recomputes the brief's numbers and compares them with
[`brief_numbers.json`](brief_numbers.json). A summary is written to `review_report.json`. One stage:

```bash
uv run --locked python run_review.py --bundle /path/to/bundle --stage reference
```

## What each stage reproduces

| # | Stage | Decides | Checked against |
|---:|---|---|---|
| 1 | `gene_list` | which genes: S1E3 filter over surfaceome release 1.3.0 + restated co-receptor & gap additions, minus exclusions and withdrawals | the 2,399 library genes exactly |
| 2 | `gpi_set` | which forms are GPI-anchored | `gpi_set.tsv`, byte for byte |
| 3 | `reference` | each gene's reference proteoform (the ladder) | `reference_selection.tsv`, `reference_disagreements.tsv`, `reference_selection_sequences.tsv` |
| 4 | `isoform_universe` | every alternate form and its evidence tier | summaries byte for byte; cohort table row for row (see the hash-seed limit in the brief) |
| 5 | `ortholog_repick` | each gene's mouse and cyno proteoform (April Compara 115 pool) | byte for byte |
| 6 | `universe_coverage` | how many library genes the isoform and ortholog pools could draw from (Compara 115) | `universe_coverage.tsv` and its summary, byte for byte |
| 7 | `isoform_gate` | which isoform-cohort forms may be drawn (native ATG, CDS, N-terminal architecture) | `isoform_gate.tsv` and pinned S3 key |
| 8 | `panels` | original and 2026-11-12 gated panel draws, plus 1:1 check of current `isoform_panel.tsv` | `isoform_panel_2026-11-12.tsv`, `ortholog_panel.tsv`, active manifest |
| 9 | `ortholog_adoption` | L1 Compara 116 CDS-eligible ortholog panel (final-v2) + correspondence review + 2026-11-22 restorations | `ortholog_panel_adopted.tsv`, manifest, S3 records |
| 10 | `rescue_arms` | 25 expression-rescue arms (retired record) | byte for byte |
| 11 | `giants` | truncations of oversized genes (L2) | byte for byte |
| 12 | `candidate_cds` | each member's designed CDS | `candidate_cds.json` (all entries end in a stop codon) |
| 13 | `reference_fixes` | named reference fixes, ATG start replacements, EPHA6 trim | `reference_changes.tsv` and S3 records |
| 14 | `ortholog_start_trims` | 4 cyno 5' Met restarts | `ortholog_start_trims.tsv` and S3 records |
| 15 | `constructs` | domesticated ORFs (BbsI/Esp3I/BsaI) and p0039 vectors | `candidate_constructs.json` and manifest vectors |
| 16 | `tranches` | L1/L2 membership, duplicates retired | no two active members share a protein |
| 17 | `topology` | DeepTMHMM v1 per designed protein | every active L1 member has a v1 row in the pinned snapshot |
| 18 | `tag_paths` | each member's placement rung | `l1_tag_paths.tsv`, `l2_tag_paths.tsv`, byte for byte |
| 19 | `cassettes` | each C-terminal member's linker encoding | every issued batch replays |
| 20 | `renumber` | plate-ordered member ids (`l1-mNNNN`) and design-id keying | `l1_member_renumbering_summary.json` |
| 21 | `cterm_review` | C-terminal proteolysis/burial review and decisions | `l1_cterm_review.tsv` and `l1_cterm_decisions.tsv` |
| 22 | `cterm_arms` | C-terminal and snorkel inserts + diversified snorkel modules | `constructs_l1_cterm.tsv/.fasta` |
| 23 | `rung_arms` | N-after-SP, multipass, type III and internal arms + diversified linkers | all four `constructs_l1_*.tsv/.fasta` tables |
| 24 | `consistency` | protein, cassette, layout, fold, and completeness across all built arms | `l1_consistency.py --check`, `l1_completeness.py --check` |
| 25 | `synthesis_compliance` | vendor synthesis check and minimal synonymous recodes | `synthesis_compliance.py --check` |
| 26 | `saved_inputs` | SHA-256 of every S3 object recorded in `library_pipeline_inputs_saved.tsv` | `pipeline_input_manifest.py` |
| 27 | `pinned_refs` | every pinned git ref is an ancestor of `HEAD` (or materialized in `historical/git-snapshots/` in exported repos) | `pinned_refs.py` |
| 28 | `structure_models` | every AlphaFold DB v6 model read by the design is archived with its SHA-256/size | `pin_structure_models.py` |
| 29 | `rung_builds` | five rung builds per cut saved on S3 | `save_rung_builds.py` |
| 30 | `external_inputs` | every external source has a saved copy or recorded tool identity | `external_inputs.py` |

Stages 4 and 5 rerun scripts at the commits that produced their records, reading from this
repository's git history or from `historical/<commit>.tar` in an exported review repo; without
either, the runner reports them as skipped, never as passed.

DeepTMHMM and NetGPI are not rerun: their predictions are inputs, pinned as committed
outputs and in the bundle's topology snapshot. Stage 17 checks that every member has a
prediction for its exact sequence, not that the predictor would return the same call.

## What is in the bundle, and what is public

Each bundle file is labeled in `bundle_manifest.tsv` with its sharing class, its release,
version and source URL, and whether the replay reads it (`read_by_replay`). Every release the
design rests on is bundled whole, not only the files the replay opens, so a reviewer can
inspect it; [`pinned_releases.md`](pinned_releases.md) lists them (release, version, date
pinned, size, place in the bundle), and `bundle.py audit` fails if anything an input manifest
or active record constant pins is missing.

| Class | Read by the replay | Rest of the pinned releases | Meaning |
|---|---:|---:|---|
| `public_release` | 2,875 | 36 | byte-identical to a named public release (AlphaFold DB v6 archive, Ensembl 115 and 116 including Compara homology dumps, UniProt, APPRIS, OpenProt 2.2, Pfam 39.0, HPA v23, GTEx v10, surfaceome 1.3.0 on Zenodo, GENCODE Ribo-seq ORFs, the published GPCR screen table) |
| `public_derived` | 21 | 47 | public data, as queried by us (PeptideAtlas build 607, ARCHS4 v2.5, MyGene.info, UniProt keyword and accession queries) |
| `surfaceome_repo` | 15 | 1 | files of the accessible-surfaceome repository at commit f9c73691 |
| `deliverome_curation` | 345 | 54 | internal inputs: curation calls, the compiled expression table, the vector backbone map, DeepTMHMM and NetGPI predictions, ColabFold structures, the ortholog draw inputs, and the audited S5E13 tag-site positives the brief cites |
| `deliverome_record` | 54 | 143 | published library records: what the replay is compared with, the Compara 116 ortholog draw and its verification, the ortholog correspondence review, and the Section 8 confirmed-site ranker benchmark tables |

Data homes in [`data_locations.tsv`](data_locations.tsv) mirror the same 3,591 files: **2,946**
public files on R2 (including APPRIS under CC BY-NC-ND 4.0 unmodified) + **1** upstream Pfam-A file
at EBI (3.53 GB total), **213** small tables committed in the review repo (`data/`, 0.09 GB), and
**431** files in the release tarballs (1.38 GB). Reviewers receive the full bundle, all five
classes (Becca, 2026-11-04). A bundle can still be packed without a class
(`bundle.py pack --exclude-class ...`); the stages that read the missing files then fail with a
message naming the absent key.

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

Reviewers get the code from the public review repository,
[Deliverome-Project/deliverome-surfaceome-library-v1](https://github.com/Deliverome-Project/deliverome-surfaceome-library-v1),
generated from this one by `export_review_repo.py`. It contains:

- the library's folders whole (s4e9, the S1E3 gene filter, `library/l1/proteoform-selection`,
  `library/l1/tag-site-architecture`);
- every other file the replay opens (listed in `review_repo_files.txt`, found by tracing a
  replay);
- the complete `deliverome_analysis` package, `pyproject.toml` and `uv.lock`;
- in `historical/`, the old script versions that stages 4 and 5 rerun.

The latest published GitHub release on that repository is
[`bundle-2026-11-26b`](https://github.com/Deliverome-Project/deliverome-surfaceome-library-v1/releases/tag/bundle-2026-11-26b)
as two tarballs: `l1-reviewer-bundle-2026-11-26b.tar.gz` (3,310 replay files, 1.46 GB compressed,
SHA-256 `a18660e48f6accf9ae10a2530abad6c8cbd07c4b0172707a9d0818dcdda3a3ae`; gzipped to stay under
GitHub's 2 GiB asset limit) and `l1-pinned-releases-2026-11-26b.tar` (the rest of every pinned
release and every S3 file the brief cites, 281 files, 1.96 GB, SHA-256
`a86c270e2fc2b53315d4f7ceaa2a63dd94bc600bdd4a2f6e0b785818af9cde80`). The working
[`bundle_manifest.tsv`](bundle_manifest.tsv) indexes those exact **3,591** files (**3,310** replay /
**281** supplement) matching order of record `2026-11-25-pre-freeze-e`.
`tar -xf` reads both.
A rebuilt bundle goes on a new release with a new date, never over an existing one; earlier
releases are `bundle-2026-11-02` (before l1-p2740..2742), `bundle-2026-11-05` (before the
audited tag-site positives were added), `bundle-2026-11-05-r2` (before the stop-codon fix of 13
designed CDS and the extended audit), `bundle-2026-11-05-r3` (before the tag-site gold
standard), `bundle-2026-11-05-r4` (before the final audit refresh), `bundle-2026-11-05-r5`
(before the universe-coverage stage and its Compara dumps), `bundle-2026-11-05-r6` (before
the Compara 116 ortholog panel was adopted and every pinned release was bundled whole),
`bundle-2026-11-10` (before the 2026-11-14 order of record and the files the brief cites),
`bundle-2026-11-20` (3,120 replay / 246 supplement files),
`bundle-2026-11-21b` (3,125 replay / 251 supplement files), and
`bundle-2026-11-26` (3,286 replay / 284 supplement files, prior to 30/30 offline replay hardening).

`EXPORT_MANIFEST.json` there records the source commit and a SHA-256 per managed file. To
update it after changes here:

```bash
uv run --locked python export_review_repo.py --dest ~/git/deliverome-surfaceome-library-v1
```

A re-export never touches files added in the review repo. A managed file that was edited there is left in place and listed,
so the edit can be carried back here by hand. If a new script is needed by the replay,
rebuild the list from a traced run (`--trace`) before exporting.
