# Genes at each L1 gene-list gate

One table per gate of the cohort filter and one per entry tier, written by
[`gate_lists.py`](../gate_lists.py) from an offline rerun of the filter over
accessible-surfaceome release 1.3.0 (the pipeline's `gene_list` stage). Each gate
table lists every gene the gate was applied to, `fail` rows first, with the value
the gate tested. Gates apply in order and stop at a gene's first failure, so a gene
appears in the table of every gate up to and including the one that excluded it,
and in no later one. The per-gene summary across all gates is
[`gene_gate_outcomes.tsv`](../../gene_gate_outcomes.tsv); the script asserts every
pass/fail here agrees with it. Topology values are DeepTMHMM v1 (1.0.24) on the
UniProt canonical sequence, as recorded in the release.

| File | Gate or tier | Genes tested | Failed (brief: "Excluded here") |
|---|---|---:|---:|
| [`01_surface_evidence.tsv`](01_surface_evidence.tsv) | surface evidence | 5,332 | 2,297 |
| [`02_accessibility.tsv`](02_accessibility.tsv) | accessibility | 3,035 | 0 |
| [`03_surface_fraction.tsv`](03_surface_fraction.tsv) | surface fraction | 3,035 | 362 |
| [`04_anchored.tsv`](04_anchored.tsv) | anchored | 2,673 | 141 |
| [`05_coreceptor.tsv`](05_coreceptor.tsv) | co-receptor | 2,532 | 160 |
| [`06_expression.tsv`](06_expression.tsv) | expression | 2,372 | 81 |
| [`07_localization.tsv`](07_localization.tsv) | localization | 2,291 | 7 |
| [`08_orf_size.tsv`](08_orf_size.tsv) | ORF size | 2,284 | 100 |
| [`09_ectodomain.tsv`](09_ectodomain.tsv) | ectodomain | 2,184 | 50 |

| File | Tier | Genes |
|---|---|---:|
| [`10_tier_core.tsv`](10_tier_core.tsv) | passes every gate (tiers T1-T3 by architecture risk) | 2,134 |
| [`11_tier_low_literature_uniprot.tsv`](11_tier_low_literature_uniprot.tsv) | low-literature UniProt surface tier | 139 |
| [`12_tier_screen_validated_gpcr.tsv`](12_tier_screen_validated_gpcr.tsv) | screen-validated GPCR tier (Tedman et al., 2026) | 85 |
| [`13_tier_overexpression_surface.tsv`](13_tier_overexpression_surface.tsv) | overexpression-surface tier | 12 |

The filter passes 2,370 genes (core 2,134 plus the three tiers); 1 removed afterwards ([`gene_list_exclusions.tsv`](../../gene_list_exclusions.tsv): EPPIN) leaves 2,369, the library's gene count. Cohort: 5,332 genes.
