# s1e3 — surface-library-filters

A **tractability cut** over the accessible-surfaceome catalog: which surface proteins,
when expressed in a cell line, reach the plasma membrane at reasonable levels *without*
requiring strong cofactors — plus which **isoforms** and **orthologs** are worth
including, on an explicit internalization / ICD-dynamics axis.

This is orthogonal to the catalog's `canonical` preset. That preset certifies *"is this a
surface protein?"*; this one asks *"will it traffic to the surface on its own when I
transfect it?"*

Promoted from `scratch/becca/surface-library-filters` (unreviewed) into reviewed
`analysis/` — see the Decision log in [DESIGN.md](DESIGN.md) for the filter-design history.

## Notion page

https://www.notion.so/3b8ab86c768381bc8a73dc79ec72d832
(ELN: **S1E3** in the Experiments database)

## Inputs

All inputs are **fetched, not committed** — nothing here depends on a local file that
only exists on one machine:

| Input | Source | Fetched by |
|---|---|---|
| 5,130-gene deep-dive cohort (`cohort.jsonl`) | `https://api.deliverome.org/surfaceome/v1/` — public, no auth, edge-cached, D1-backed | [`fetch_cohort.py`](fetch_cohort.py) (resumable; re-runs only fill gaps) |
| Expression TPM reference | HPA / consensus TPM | [`pipeline/build_expression_reference.py`](pipeline/build_expression_reference.py) — see [`external/expression_tpm/PROVENANCE.md`](external/expression_tpm/PROVENANCE.md) |
| Willow GPCR screen | third-party screen | see [`external/willow_gpcr_screen/PROVENANCE.md`](external/willow_gpcr_screen/PROVENANCE.md) |
| Symbol → Ensembl canonicalization | MyGene, via `deliverome_analysis.gene_resolver` | resolver cache (keyed by input-symbol-set hash) |

**Preset pinning:** `presets_pr130.py` / `presets_pr137.py` copy the `passes_canonical` /
`passes_likely` predicates verbatim from the corresponding accessible-surfaceome PRs
rather than importing from its `main`. Re-pin if those presets change.

## Method

- Fetch the deep-dive cohort once; cache to `cohort.jsonl` so re-runs are free.
- Layer a `displayable` tier cut on top of `canonical` / `likely`, using fields those
  presets omit (co-receptor dependency, membrane architecture, state dependence).
- Flag rather than exclude, except for hard obligate-cofactor cases.
- Score isoforms and orthologs on an internalization / ICD axis.
- Co-receptor **ubiquity** replaces the earlier hand-curated housekeeping heuristic — see
  [COReceptor_UBIQUITY.md](COReceptor_UBIQUITY.md).

## Outputs

Figures land in [`figures/`](figures/) (via `save_figure()`, per project plotting style):

| Figure | Produced by |
|---|---|
| `figures/tier_gene_isoform_ortholog.png` | [`tier_plot.py`](tier_plot.py) |
| `figures/tier_variant_coverage.png` | [`tier_plot.py`](tier_plot.py) |
| `figures/category_by_tier_with_total.png` | [`tier_plot.py`](tier_plot.py) |
| `figures/isoform_outcomes.png` | [`plots.py`](plots.py) |
| `figures/ortholog_outcomes.png` | [`plots.py`](plots.py) |
| `figures/reference_set_coverage.png` | [`reference_coverage.py`](reference_coverage.py) |
| `figures/surface_library_funnel.png` | [`build_library.py`](build_library.py) |
| `figures/category_totals_by_tier.png` | [`tier_plot.py`](tier_plot.py) |

Tables / notes: `tier_counts.json`, [`counts.md`](counts.md),
[`willow_validation.md`](willow_validation.md),
[`list_comparison_surfy_tpm.md`](list_comparison_surfy_tpm.md).

Library export to S3 + viewer: [`export_fg_library.py`](export_fg_library.py).

## Status

2026-09-10 partner-rule amendment: **2,340 -> 2,344 genes**, adding ATP1B3/CD298,
CDC50A/TMEM30A, HFE and P2RX6. See [CORECEPTOR_AUDIT.md](CORECEPTOR_AUDIT.md)
for the six reviewed partner routes, validation, and the local frozen input/output bundle.
The amendment is pending dev review; canonical published exports remain at 2,340.

Promoted to reviewed. Filter design decisions are locked (graded tiers, flag-not-exclude
architecture, hard obligate-cofactor exclusion, ICD/internalization axis for isoforms and
orthologs); see [DESIGN.md](DESIGN.md) §1–§3 first if picking this up cold.

## Signal-peptide review snapshot

[Complete SP prioritization and reproduction](SIGNAL_PEPTIDE_REVIEW.md): exact
2,344 main-protein classifications, including the 519 precedent-informed and
1,667 broad candidates, archived inputs/outputs, and Notion downloads.
