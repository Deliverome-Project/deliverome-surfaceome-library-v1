# Current library vs the original "Surfy + TPM" list

**Question:** how different is the current filter-based library from the original Deliverome candidate list (Surfy surfaceome membership, expression-filtered) on deliverome-internal, ~2,200 genes?

**Reconstruction** (approximate — rebuilt in-repo from `catalog_rows.json` Surfy `db` bit + the full re-derived max-TPM reference; the exact deliverome-internal list may differ slightly in the Surfy release / TPM source):
- Surfy-positive genes in our 19,324 universe: **2,733**
- **Surfy ∩ max-TPM ≥ 10 = 2,303** ← matches the "~2,200" recollection
- Current filter-based library: **2,346**

## They're the same size, but only 62% the same genes

| | count |
|---|---|
| in **both** | 1,787 |
| **current only** (not in Surfy+TPM) | 559 |
| **Surfy+TPM only** (dropped by our filter) | 516 |
| Jaccard overlap | **62%** |

Near-identical totals (2,346 vs 2,303) hide ~¼ turnover in each direction. The current filter is **not** "Surfy + expression" — it trades genes.

## Why 516 Surfy+TPM genes are DROPPED (primary failing gate)
Our filter is stricter than Surfy∩TPM on evidence, cofactor-dependence, localization, and size:

| dropped | primary reason |
|---|---|
| 232 | **below-Likely** — weak/conflicting deep-dive evidence grade (Surfy's SVM called them surface; our evidence gate doesn't) |
| 103 | **obligate co-receptor** — heterodimer whose partner isn't ubiquitous (no co-delivery) |
| 93 | **oversized > 5 kb** — on the truncation side-list, not truly lost |
| 75 | **mostly_intracellular** — deep-dive says predominantly internal |
| 6 / 4 / 3 | not anchored / non-surface reason / ECD < 30 |

## What 559 current-only genes ADD
- **399 are NOT in Surfy at all** — surface proteins Surfy's SVM missed, caught by the deep-dive (GO/CSPA/HPA/UniProt evidence + DeepTMHMM topology). Heavily **multi-pass**: T2 194, T1 106, T3 55, LL 42. Examples: ABCC6, AQP3/7, ATP2B1, ATP11A/C, ANO4, ASIC2/3, BEST1/3, BMPR1B — the channels / transporters / P-type ATPases that Surfy's feature model under-calls.
- **160 are Surfy-positive but TPM < 10** — kept via the `low_endogenous_expression = False` path or the SV/OE/co-receptor rescues (low endogenous but real targets): FSHR, GALR3, GPR52, HRH4, GFRAL, the KIR family, GLRA1.

## Takeaway
Same ballpark size as the original Surfy+TPM cut, but a **different list** (62% overlap). The filter removes ~516 Surfy genes it judges non-surface / non-targetable / cofactor-dependent / oversized, and adds ~559 (mostly multi-pass channels & transporters Surfy missed, plus low-but-real-TPM receptors rescued on evidence). Per-gene detail in `list_comparison_surfy_tpm.tsv`.
