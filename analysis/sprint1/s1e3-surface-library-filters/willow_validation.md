# Validation vs Willow's GPCR surface screen (Tedman et al. 2025)

Cross-check of the s1e3 surface-library filter against direct experimental surface-expression
data for GPCRs (and their isoforms). The screen FACS-sorts GPCR ORFs on cell-surface
immunostaining in a cultured line **without olfactory RTP chaperones** — so a high value means
the receptor reaches the surface autonomously when expressed. This is exactly the question our
filter tries to answer with proxies, so it's a clean ground-truth test.

- **Dataset:** `s3://deliverome-external/tedman-gpcr-surface-screen/2025-09-19/` (see its
  `PROVENANCE.md`). bioRxiv 2025.09.19.677468.
- **Surface-positive threshold:** anti-HA immunostaining **≥ 1000** (`SCREEN_MIN`).
- **Scripts:** `analyze_willow_gpcr.py` (concordance), `build_isoform_map.py` (ENST→UniProt-
  isoform), `willow_rescue.py` (screen-based SV inclusion). Outputs: `willow_*.csv/json`.

> **Refreshed 2026-08-15 on the fresh cohort + live gate stack.** `analyze_willow_gpcr.py` now
> imports `final_library`'s gate functions directly, so this validation can no longer drift from
> the real filter. The live gates are: membership (canonical ∪ likely), `surface_accessibility ≠
> no`, `surface_specificity ∈ {surface_dominant, mixed}`, anchored, co-receptor (INCLUDE_ALONE
> allowlist + ubiquitous-partner rule), expression (`low_endog=False` OR max TPM ≥ 5 on the full
> re-derived reference), non-nonsurface reason, ORF ≤ 5 kb, ECD ≥ 30. (Supersedes the earlier
> TPM≥10 / no-membership / no-ECD run.)

## Canonical GPCRs through the CURRENT core filter (screen-independent)
763 screened GPCRs, 758 in our cohort. This is the honest test of the proxy filter **without**
the screen-based SV rescue:

| class | in-cohort | core-pass | rate |
|---|---|---|---|
| non-olfactory (A/B/C/F) | 338 | 314 | **93%** |
| olfactory / taste / vomeronasal | 420 | 72 | 17% |
| **all** | 758 | 386 | **51%** |

Core passers are T2 (287) + T3 (99), **zero T1** — correct, every GPCR is 7-TM (never single-pass
"low-risk"). Of the 372 in-cohort GPCRs the core filter excluded, gate-failure counts are
**expr 284, membership 273, spec 37, coreceptor 14** — exclusions are driven by the expression
proxy + membership, not topology. The 93% non-olfactory rate is the meaningful number; olfactory
receptors (correctly low-expression / below-Likely) drag the overall figure down.

## Library-inclusive coverage (core ∪ SV ∪ LL)
The **SV (screen-validated)** tier adds back GPCRs the core filter excluded but the screen shows
at the surface (immunostaining ≥ 1000). Counting the full library:

| class | in-cohort | in-library | rate |
|---|---|---|---|
| non-olfactory | 338 | 323 | **96%** |
| olfactory / taste | 420 | 199 | 47% |
| **all** | 758 | 522 | **69%** |

**Well-surface-expressed screen GPCRs (immunostaining ≥ 1000): 275. The library captures 273 —
only 2 misses:** GPR182 (not in our cohort at all — a catalog gap, flag upstream) and one
olfactory boundary case. So of the receptors the screen shows *clearly* at the surface, coverage
is ~99%.

## SV tier + totals
85 GPCRs (screen ≥ 1000, core-excluded) included as the flagged SV tier, plus 13 native
screen-validated isoform variants. Current library total **2,330** (2,119 core + 126 LL + 85 SV;
incl. ITGB1).

## Isoform concordance (ENST→UniProt-isoform map)
174 screened isoform transcripts, 105 genes; 163/174 parents pass the core filter. **13 "make it
through"** (parent passes AND our topology calls the matched isoform surface). The dominant limit
is annotation coverage: many screened transcripts are **Ensembl-only** (including the `_SP_HA`
constructs) with no UniProt isoform, so they cannot appear in our UniProt-isoform-keyed catalog —
an annotation-source limit, not a filter decision. The screen confirms **133/174** isoforms lose
surface expression vs their canonical (% change < 0), supporting the topology-first isoform rule
(judge surface-attachment, not % identity).

## Filter issues — status vs the original run
1. **Expression gate over-filtered GPCRs** — softened (TPM ≥ 5 on the full re-derived reference)
   and the SV tier now overrides the proxy for screen-positive receptors. *(addressed)*
2. **`surface_specificity = mostly_intracellular` false negatives** — OXTR (5387), ACKR2 (6397),
   TAAR1 are now captured by the SV tier (experimental surface evidence overrides the LLM
   specificity call). *(addressed)*
3. **Obligate-co-receptor direction (GABBR2 traffics alone)** — now handled by the INCLUDE_ALONE
   allowlist built from the co-receptor subagent calls (`coreceptor_calls.tsv`). *(addressed)*
4. **Olfactory blanket T3 demotion partly contradicted** — screen-positive ORs now enter as SV
   rather than being discarded. *(addressed)*
5. **Isoform annotation coverage** remains the isoform bottleneck (Ensembl-only transcripts have
   no UniProt isoform). *(open — upstream data source, not a filter change)*
