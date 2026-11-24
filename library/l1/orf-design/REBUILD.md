# Rebuilding the L1 order from zero

This page is the one procedure for rebuilding the L1 library order from a clean machine:
which tools are needed, which steps run in which order, and what each step proves. The
detailed procedures it points to are listed at the end. It is current as of the final
regeneration (2026-11-14, order of record `synthesis/2026-11-14-final/`).

## What "rebuild" means here

Every L1 design decision is a committed table, and every input is saved with its SHA-256.
A rebuild reruns each step from those saved inputs and must reproduce the committed tables
byte for byte. It does not re-ask today's databases or re-run the predictors: the library
is defined by what was decided, so a step is checked against its record, never redone from
current sources. Re-predicting (SignalP 6, DeepTMHMM, ColabFold) is possible but is a new
analysis, not a rebuild.

## 1. Set up

| Need | How | Notes |
|---|---|---|
| The repository with its full history | `git clone https://github.com/Deliverome-Project/deliverome-analysis.git` (not `--depth 1`) | stages replay pinned commits; `pinned_refs` refuses a shallow clone |
| Python 3.14 and the locked environment | `uv sync --locked` | the clean rebuild of 2026-11-14 ran Python 3.14.10; `uv.lock` pins every package (its SHA-256 is in `library_pipeline_inputs_saved.json`) |
| R 4.6 and Bioconductor 3.23 (optional) | `scripts/install-r-deps.sh` | not needed by any L1 stage |
| Read access to the `processed` and `external` S3 buckets | credentials in `.env` (see `.env.example` and the repo's getting-started guide) | everything a stage reads is listed in `library_pipeline_inputs_saved.tsv` (3,104 objects) and `external_inputs.tsv` |
| Empty caches, to prove nothing comes from a local copy | `export DELIVEROME_ANALYSIS_S3_CACHE=$(mktemp -d) DELIVEROME_CACHE_DIR=$(mktemp -d)` | optional for a routine check |

Licensed tools (DeepTMHMM 1.0.24, SignalP 6.0) are **not** needed to rebuild: their
predictions are saved. They are needed only to re-predict, and their identity (version,
install source, SHA-256 of every file) is in `library_pipeline_inputs_saved.json` and
`reviewer-package/nterm/signalp6_identity.json`, so a licensed copy can be checked.

## 2. Rebuild and check everything the order depends on

```bash
cd library/l1/orf-design
uv run --locked python library_pipeline.py plan      # the 30 stages, what each decides
uv run --locked python library_pipeline.py verify    # rerun every stage from saved inputs
```

The stages run in this order. Each prints `OK` with what it reproduced, or `FAIL` with why.

| Group | Stages | What they reproduce |
|---|---|---|
| Genes and proteoforms | `gene_list`, `gpi_set`, `reference`, `isoform_universe`, `ortholog_repick`, `universe_coverage`, `isoform_gate`, `panels`, `ortholog_adoption`, `rescue_arms`, `giants`, `candidate_cds`, `reference_fixes`, `ortholog_start_trims`, `constructs`, `tranches`, `topology`, `tag_paths` | which genes, which protein per member, its designed CDS and its DeepTMHMM v1 rung |
| Tagged arms | `cassettes`, `renumber`, `cterm_review`, `cterm_arms`, `rung_arms`, `consistency` | the cassette encodings, member ids and plates, every rung's tagged arm, and that every built arm is consistent and every refold failure has a disposition |
| The order | `synthesis_compliance`, `rung_builds` | the synthesis recodes of every insert; the rung builds the order was made from, saved on S3 |
| The inputs themselves | `saved_inputs`, `pinned_refs`, `structure_models`, `external_inputs` | every saved input present and unchanged; every pinned commit on `main`; every AlphaFold DB model the design read archived under `external://alphafold-db/v6/`; every external source listed with its saved copy |

Two inputs come from steps outside this pipeline, each with its own record:

| Step | Rebuild | Record |
|---|---|---|
| Internal-site decisions (S5E13 ranker v5, contact veto, refold verdicts) | `cd library/l1/tag-site-benchmarking && uv run --locked python l1_internal_decisions.py --restore && uv run --locked python l1_internal_decisions.py --check` | `l1_internal_decisions.tsv` / `.json`; the D1 rows it reads are saved (`d1_snapshot_2026-11-14.json`) |
| Membrane-side gate | `cd membrane-side-gate && uv run --locked python pin_inputs.py restore`, then the gate's rerun in its README | `membrane-side-gate/*.tsv`; 4,402 inputs in one hashed bundle |
| Tagged refolds | ColabFold on Modal (`scripts/fold_tagged_modal.py`); a new run, not a rebuild | every model saved with its SHA-256 in `tag_refold_structures.tsv`; `reviewer-package/fold_order_check.py` re-checks each against the order |

## 3. Re-make the order of record

Only after `verify` passes, and only into a new dated key (published files are never
overwritten):

```bash
uv run --locked python synthesis_compliance.py --check            # tables reproduce
uv run --locked python reviewer-package/fold_order_check.py       # every insert = its fold
uv run --locked python l1_consistency.py --check                  # dispositions complete
(cd reviewer-package && uv run --locked python order_plan.py)     # plates and pooling
```

To publish a changed order: bump `FASTA_KEY` in `synthesis_compliance.py` to a new date,
commit, then `synthesis_compliance.py --publish` from the clean tree (the S3 object is
stamped with the commit).

## 4. Regenerate the brief

```bash
cd reviewer-package
uv run --locked python brief_numbers.py        # every count and figure, from the tables
uv run --locked pytest -q ../tests/test_brief_prose_numbers.py   # prose = tables
uv run --locked python render_brief.py --dest <folder>           # PDF, Word and Markdown
```

## Last full rebuild

2026-11-14, a fresh clone of `main` at `f5222462` (the final-c order of record merged) with
empty caches: 30 of 30 stages reproduced, and the order plan regenerated its tables byte for
byte. The 73 AlphaFold DB models not yet archived are reported by `structure_models` as a
listed exception (`tag_refold_afdb_unarchived.tsv`).

## Detailed procedures

| Topic | Where |
|---|---|
| Every stage, script and table | [`CODE_MAP.md`](CODE_MAP.md) |
| External sources and their saved copies | [`external_inputs.tsv`](external_inputs.tsv) |
| The design and every decision | [`reviewer-package/BRIEF.md`](reviewer-package/BRIEF.md) |
| N-after-SP sites and SignalP runs | [`reviewer-package/nterm/README.md`](reviewer-package/nterm/README.md) |
| Tagged refolds | [`tag_refold/README.md`](tag_refold/README.md) |
| Membrane-side gate | [`membrane-side-gate/README.md`](membrane-side-gate/README.md) |
| Internal-site ranker | `library/l1/tag-site-benchmarking/README.md`, "Reproducing v5 and the internal decisions" |
