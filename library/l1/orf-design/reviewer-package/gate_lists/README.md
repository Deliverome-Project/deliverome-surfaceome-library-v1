# Genes at each L1 gene-list gate

Per-gate and per-tier gene tables written by
[`gate_lists.py`](../gate_lists.py) from an offline rerun of the filter over
accessible-surfaceome release 1.3.0 (the pipeline's `gene_list` stage). For each gate,
`*_passed.tsv` holds only the genes remaining after the gate, `*_excluded.tsv` holds only
the genes excluded at that gate (so each file's row count matches the count shown in the
brief), and the combined `XX_<gate>.tsv` holds all genes tested at that gate (`fail` rows
first). Gates apply in order and stop at a gene's first failure, so a gene appears up to
and including the gate that excluded it, and in no later one. The per-gene summary across
all gates is [`gene_gate_outcomes.tsv`](../../gene_gate_outcomes.tsv); the script asserts
every pass/fail here agrees with it. Topology values are DeepTMHMM v1 (1.0.24) on the
UniProt canonical sequence, as recorded in the release.

| Gate | Genes tested (combined) | Genes remaining (`*_passed.tsv`) | Genes excluded (`*_excluded.tsv`) |
|---|---:|---:|---:|
| surface evidence | [5,332](01_surface_evidence.tsv) | [3,035](01_surface_evidence_passed.tsv) | [2,297](01_surface_evidence_excluded.tsv) |
| accessibility | [3,035](02_accessibility.tsv) | [3,035](02_accessibility_passed.tsv) | [0](02_accessibility_excluded.tsv) |
| surface fraction | [3,035](03_surface_fraction.tsv) | [2,673](03_surface_fraction_passed.tsv) | [362](03_surface_fraction_excluded.tsv) |
| anchored | [2,673](04_anchored.tsv) | [2,532](04_anchored_passed.tsv) | [141](04_anchored_excluded.tsv) |
| co-receptor | [2,532](05_coreceptor.tsv) | [2,372](05_coreceptor_passed.tsv) | [160](05_coreceptor_excluded.tsv) |
| expression | [2,372](06_expression.tsv) | [2,291](06_expression_passed.tsv) | [81](06_expression_excluded.tsv) |
| localization | [2,291](07_localization.tsv) | [2,284](07_localization_passed.tsv) | [7](07_localization_excluded.tsv) |
| ORF size | [2,284](08_orf_size.tsv) | [2,184](08_orf_size_passed.tsv) | [100](08_orf_size_excluded.tsv) |
| ectodomain | [2,184](09_ectodomain.tsv) | [2,134](09_ectodomain_passed.tsv) | [50](09_ectodomain_excluded.tsv) |

| File | Tier | Genes |
|---|---|---:|
| [`10_tier_core.tsv`](10_tier_core.tsv) | passes every gate (tiers T1-T3 by architecture risk) | 2,134 |
| [`11_tier_low_literature_uniprot.tsv`](11_tier_low_literature_uniprot.tsv) | low-literature UniProt surface tier | 119 |
| [`12_tier_screen_validated_gpcr.tsv`](12_tier_screen_validated_gpcr.tsv) | screen-validated GPCR tier (Tedman et al., 2026) | 85 |
| [`13_tier_overexpression_surface.tsv`](13_tier_overexpression_surface.tsv) | overexpression-surface tier | 12 |

The filter passes 2,350 genes (core 2,134 plus the three tiers); 1 removed afterwards ([`gene_list_exclusions.tsv`](../../gene_list_exclusions.tsv): EPPIN) leaves 2,349, the library's gene count. Cohort: 5,332 genes.
[`14_library_genes.tsv`](14_library_genes.tsv) lists all 2,399 library genes: the 2,349 above and 50 dated additions ([`gene_list_additions.tsv`](../../gene_list_additions.tsv)).
