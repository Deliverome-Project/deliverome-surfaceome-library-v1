# Surface-protein library — filter design

**Status:** Draft for review
**Date:** 2026-08-04
**Owner:** Becca Carlson
**Home:** `library/l1/surface-filters/` (promoted from `scratch/becca/surface-library-filters`)
**Related:** accessible-surfaceome PR130 (`feat(surface-proteome): gate canonical on overall confidence`); s1e2 fragmid vector design (the 2,243-ORF synthesis set == canonical count)

---

## 1. Purpose & scope

Produce an **initial list of surface proteins for the Deliverome library**: proteins that, **when expressed in a cell line, reach the plasma membrane at reasonable levels without requiring strong cofactors**. The list is a *tractability* cut over the accessible-surfaceome catalog, plus an assessment of which **isoforms** and **orthologs** are worth including — with an explicit **internalization / ICD-dynamics** axis, since delivery (receptor-mediated uptake) is the point of the library.

Out of scope for this pass: binder/epitope design, dosing, construct-level cloning design.

## 2. Background — why this is a *new* cut, not a tweak to `canonical`

The catalog's `canonical` preset (2,243 genes after PR130) certifies **"is this a surface protein?"**. It *deliberately* ignores heterologous-expression tractability — its docstring keeps `state_dependence='high'` genes and never looks at co-receptor dependency or membrane architecture, treating "*when* it is surface" as a facet, not a gate.

Our requirement is orthogonal: **"will it traffic to the surface on its own when I transfect it?"** So this is a new `displayable` cut **layered on top of** `canonical`/`likely`, using fields those presets omit.

## 3. Data source & schema

- **API:** `https://api.deliverome.org/surfaceome/v1/` (public, no auth, edge-cached).
  - `/v1/genes` → 5,130 deep-dive cohort symbols (summary fields only).
  - `/v1/genes/{SYMBOL}` → full `SurfaceomeRecord` (JSON), carrying `filters` **and** `deterministic_features` (topology + isoform + ortholog **sequences + per-residue topology**).
- **Schema:** `viewer/lib/surfaceome-types.ts` (v1.0.0 / served as 1.1.x). Preset predicates: `src/accessible_surfaceome/release/catalog_presets.py` (mirror: `viewer/lib/catalog-presets.ts`).
- **Preset version — pin to PR130.** PR130 is still **OPEN**; `main` still ships the *old* canonical (`evidence_grade ∈ {direct_multi_method, direct_single_method}` → **1,452** genes). This analysis copies `passes_canonical` / `passes_likely` **verbatim from the PR130 branch** (`claude/canonical-constitutive-gate`) — the confidence-gated recast (**1,452 → 2,243**, +791 incl. ICAM1-class) — never importing from `main`. Re-pin if PR130 changes before merge. PR130 is a classifier change, not a re-annotation, so the live API `filters` fields are unaffected: live data + PR130 predicate reproduces the 2,243.
- **One fetch pass** over the 5,130 full records supplies everything Components A–D need; cache to `cohort.jsonl` (compact — keep only `filters` + `deterministic_features.{canonical_topology, isoform_topologies, orthologs}`) so re-runs are free.

## 4. Requirements (locked with Becca)

| # | Decision |
|---|---|
| R1 | **Graded tiers** (T1 core / T2 extended / T3 full-flagged), not one set. **Membership base = `canonical` ∪ `likely`** (both); membership provenance (`canonical` vs `likely_only`) is carried as a flag and is the top confidence axis. The additional red-line filters + flags apply across the whole union. |
| R2 | **All architectures kept, risk-flagged** — architecture is a tiering axis, not a hard exclusion. |
| R3 | **Obligate cofactor = hard red line**: drop `co_receptor_dependency='required'` OR `surface_call_reason='stable_complex_partner'`. |
| R4 | **Secreted-only excluded** — won't attach on its own (`surface_call_reason='secreted_only'`). |
| R5 | **Ignore** `overexpression_surface_localization_observed`. |
| R6 | **Soft** filter on `expression_level` + `expression_breadth` (rank/down-weight only, never exclude). |
| R7 | **Isoforms**: inclusion is **topology-only** — keep surface-attached topologies, drop secreted/intracellular ones; do **not** bin on magnitude of sequence change (Tedman: divergence doesn't predict expression). ECD + full-length %-identity differences are **annotated for review, not gates**. GPCR measured data = soft estimate only. |
| R8 | **Orthologs**: IN, evaluated by the **same topology-first logic as isoforms** — keep surface-attached, flag topology changes; annotate ECD + full-length %id vs human canonical to surface which orthologs diverge meaningfully. **No internalization / ICD-motif axis** (descoped). High-identity surrogates already live on the site; this pass adds the topology/divergence evaluation. |

## 5. Design

### 5A. Per-gene "displayable" tiers

**Red lines (exclude from all tiers):**
1. Fails membership — neither `passes_canonical(filters)` nor `passes_likely(filters)` is `True`. Carry `membership_tier ∈ {canonical, likely_only}` as a flag (`canonical` if it passes canonical; `likely_only` if it passes only likely) — this is the top confidence axis.
2. `co_receptor_dependency == "required"` (obligate partner for surface expression).
3. `surface_call_reason ∈ {stable_complex_partner, secreted_only, inner_leaflet_anchored, cytoplasmic, nuclear, mitochondrial_internal, endomembrane_resident, nuclear_envelope, pmhc_only_intracellular}`.

**Flags (computed on survivors — flag, don't exclude):**

```python
# arch_risk — trafficking autonomy by architecture
low  = subcategory in {single_pass_T1, single_pass_T2, GPI_anchored, tetraspanin}
med  = subcategory in {multi_pass, GPCR} and not is_olfactory_or_taste(uniprot_family, hgnc_gene_groups)
high = (subcategory == GPCR and is_olfactory_or_taste(...))   # RTP-chaperone-dependent; barely express in HEK293T
       or subcategory == "other"

# conditionality — does it present without a state change?
constitutive = (state_dependence in {low, moderate}
                and low_endogenous_expression is False
                and surface_call_reason not in {cell_state_induced, lysosomal_exocytosis,
                                                dual_localization, tissue_restricted_surface}
                and induction_trigger == "none")
conditional  = not constitutive

# retention — does it stay on the surface?
retained = (surface_specificity == "surface_dominant"
            and not has_shed_form and not has_secreted_form)
leaky    = not retained

# expression_prior — SOFT, ranking only (R6)
expr_prior = {high:2, moderate:1, low:0, absent:0}[expression_level] \
           + (1 if expression_breadth in {pan_tissue, broad} else 0)
```

**Tiers:**

| Tier | Predicate | Intent |
|---|---|---|
Core acceptance (shared by T1+T2): `passes_canonical` · `constitutive` · `surface_specificity ∈ {surface_dominant, mixed}` (already implied by canonical) · `co_receptor_dependency ∈ {none, modulatory}`. **T1 vs T2 differ only on architecture.**

| Tier | Predicate | Intent |
|---|---|---|
| **T1 — Core "will present"** | core acceptance · `subcategory ∈ {single_pass_T1, single_pass_T2, GPI_anchored, tetraspanin}` (`arch_risk=low`) | Autonomous-trafficking architectures; highest odds it presents when transfected |
| **T2 — Extended** | core acceptance · `subcategory ∈ {multi_pass, GPCR}` & not olfactory/taste (`arch_risk=med`) | Adds multi-pass + non-olfactory GPCRs (less reliable trafficking) |
| **T3 — Full flagged** | everything else in `canonical ∪ likely` clearing the pool red lines — incl. `passes_likely`-only, `arch_risk=high` (olfactory/taste/other), `conditional` (`low_endogenous_expression=True` or induced `surface_call_reason`), `co_receptor='unknown'` · every flag exposed | The "everything, annotated" universe for re-cutting |

First-pass counts (2026-08-04, live cohort=5,130): pool=2,582 (canonical 2,010 + likely-only 572) → **T1=731 · T2=622 · T3=1,229**. Validation: raw `passes_canonical`=2,243 (== PR130).

Rank within tier by `expr_prior`, then `evidence_density`, then `ecd_accessibility_class`.

### 5B. Isoform candidates (R7)

**Inclusion is topology-only — do NOT bin on the size/magnitude of sequence change.** Tedman et al. show divergence magnitude doesn't predict surface expression, so an ECD/ICD %-identity window is the wrong gate. Judge each isoform in `deterministic_features.isoform_topologies[]` purely on whether its topology keeps it a **surface-attached membrane protein**:

- **Keep** if the isoform retains a plasma-membrane topology with an extracellular-exposed region (≥1 TM / signal-anchor **plus** an `O` (outside) segment in `per_residue_topology`).
- **Exclude** if the topology makes it **not surface-attached**: a **secreted isoform** (no TM/anchor → soluble — not useful, it won't attach on its own) or a fully intracellular / inner-leaflet isoform.
- **Flag, don't drop,** isoforms that keep an anchor but change TM count or ECD extent (`topology_changed`) — informative (ER-retention risk per Tedman), surfaced with the flag rather than gated on magnitude.

**Annotate (never gate):** ECD %id and full-length %id vs canonical (the gap between them flags divergence *outside* the ECD, e.g. in the ICD) — reported so we can see which *kept* isoforms diverge, but they do not decide inclusion.

**GPCR estimate (soft, R7):** if the gene is a GPCR and the isoform appears in Tedman et al. Doc. S1, attach its **measured** surface expression as a prior (`gpcr_measured ∈ {expressed, low, null}`) — an estimate of importance, never a gate.

Output columns: `gene, isoform_id, topology_status ∈ {surface, secreted, intracellular}, tm_delta, ecd_pct_id, full_length_pct_id, gpcr_measured, recommendation ∈ {include, flag}`.

### 5C. Ortholog candidates (R8)

Evaluate each gene's mouse / cyno ortholog (`deterministic_features.orthologs.{mouse,cynomolgus}[]`, which carry sequence + projected per-residue topology) by the **same topology-first logic as isoforms** (§5B) — no internalization / ICD-motif axis:

- **Keep** if the ortholog preserves a surface-attached topology (≥1 TM/anchor + an extracellular `O` segment); **flag** topology changes (`topology_changed` — e.g. a human TM helix `tm_absent_from_model` in a truncated ortholog); **drop** secreted/intracellular.
- **Annotate (never gate):** `ecd_pct_identity_to_human_canonical` + full-length identity + `tm_delta` — surfaces which orthologs are faithful surface surrogates vs which **diverge meaningfully and warrant separate evaluation**. (The high-identity surrogates already on the site fall out as the high-ECD-identity, topology-matched rows.)

Output columns: `gene, species, topology_status ∈ {surface, secreted, intracellular}, tm_delta, ecd_pct_id, full_length_pct_id, recommendation ∈ {surrogate, evaluate, drop}`.

### 5D. Shared topology parse (engine for 5B & 5C)

`per_residue_topology` (`M`/`O`/`I`/`S`/`B` per residue, aligned 1:1 with `sequence`; isoforms/orthologs also carry `per_residue_topology_canonical_frame`) → derive, per variant: TM-helix count, presence of an extracellular (`O`) segment, signal-peptide/anchor presence, ECD extent, and `topology_status ∈ {surface, secreted, intracellular}`. Compare against the human canonical to set `tm_delta` and the topology flags used in §5B / §5C.

> **Removed:** the earlier ICD / internalization-motif scan (endocytic-sorting-motif detection). Internalization is out of scope for a surface-presentation filter; isoform/ortholog divergence is captured by the ECD + full-length %-identity annotations, not motif scanning.

## 6. Calibration — Tedman et al. 2025 (deep receptor scanning)

*Deep Receptor Scanning Reveals General Sequence Constraints on GPCR Biosynthesis*, bioRxiv 10.1101/2025.09.19.677468 (Coyote-Maestas et al.). Measured plasma-membrane expression of 766 GPCRs + 174 validated splice variants in HEK293T. Load-bearing findings:

- **133/174 (~76%) GPCR isoforms show little/no surface expression**; 65 match their canonical, 20 exceed it by ≥20%.
- **Divergence magnitude does not predict expression** — no indel-size relationship; CXCR1 vs CXCR2 (77% identical) differ **22×** in surface expression. ⇒ **no "% identity = too different" gate**; key on topology.
- **Structural predictors** (RF AUROC 0.81): intact TM topology (truncations fail — ER-retained/dominant-negative), total SASA > ~1900 Å², hydrophobic TM helices (esp. TMD6).
- **Olfactory/taste/vomeronasal barely express** in HEK293T (no RTP chaperones) → supports `arch_risk=high` for those.
- Even canonical GPCR surface expression spans **4 orders of magnitude**; ~64% only modest.

**Design consequence:** isoform inclusion is **topology-only** — keep surface-attached topologies, drop secreted/intracellular ones; never gate on magnitude of sequence change (it doesn't predict expression). Sequence/ICD differences are annotated for review, not used to include or exclude. GPCR isoforms carry the measured Doc. S1 estimate as a soft prior.

## 7. Outputs

- `displayable_tiers.tsv` — one row per gene: tier, all flags, raw facets.
- `isoform_candidates.tsv` — §5B columns.
- `ortholog_candidates.tsv` — §5C columns.
- `counts.md` — per-gate / per-tier funnel counts (with any silent-drop notes).
- `funnel.pdf` / `funnel.png` — cascade figure via `setup_plotting_style()` + `save_figure()` (project style; no titles).
- `README.md` — summary + exact reproduce command.

## 8. Implementation plan

1. Scaffold `scratch/becca/surface-library-filters/` (uv env; `experiment.py` importing `setup_plotting_style`, `save_figure`).
2. `fetch_cohort.py` — `/v1/genes` → symbols; fetch full records concurrently (bounded ~8–12, backoff on 429), extract compact record → `cohort.jsonl`.
3. `displayable.py` — Component A predicates → tiers + counts + funnel figure.
4. `topology.py` (5D — shared parse of `per_residue_topology` → TM count / ECD / orientation / `topology_status`) → `isoforms.py` (5B) + `orthologs.py` (5C), both topology-first.
5. *(Optional, ask first)* pull Tedman Doc. S1 for the GPCR isoform estimate (external file download).
6. Review counts with Becca; iterate thresholds; **promote to a numbered experiment** (`s1e{N}`) once the cut looks right.

## 9. Status / open questions

- **Membership base:** RESOLVED — `canonical ∪ likely` (both), pinned to the PR130 predicate (§3).
- **Home:** RESOLVED — `scratch/becca/surface-library-filters/`.
- **Internalization:** REMOVED — out of scope. The ICD-motif module is gone; §5D is now just the shared topology parse. Isoform/ortholog divergence is captured by ECD + full-length %-identity annotations.
- **Orthologs (§5C):** RESOLVED — IN, evaluated topology-first (same logic as isoforms); no internalization axis.
- **Pull Tedman Doc. S1?** deferred to a later pass (soft, GPCR-only estimate; separate download).

## 10. References

- Tedman et al. 2025, bioRxiv 10.1101/2025.09.19.677468 — GPCR isoform surface expression (calibration).
- accessible-surfaceome PR130 — canonical preset rework; `catalog_presets.py`, `surfaceome-types.ts`.
- Surfaceome public API skill — `viewer/public/surfaceome-api.skill.md`.
