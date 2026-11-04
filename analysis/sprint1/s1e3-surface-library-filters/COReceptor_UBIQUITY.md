# Co-receptor ubiquity rescue — implementation

## Current amendment (2026-09-10)

The operative threshold in `final_library.py` is **median HPA cell-line TPM >= 5**.
The max-TPM implementation and counts below are historical, not the current definition.
The six reviewed routes in `coreceptor_rules.py` distinguish alternative partners from
simultaneously required roles. Four integral surface complexes also receive a reviewed
exception to the stable-complex reason veto. This changes the reproduced 2,340-gene
library to **2,344**. See [CORECEPTOR_AUDIT.md](CORECEPTOR_AUDIT.md) for all affected
proteins, evidence, deferred cases, and exact reproduction instructions.

## Problem
The catalog flags **333** cohort genes `co_receptor_dependency == 'required'` (need a partner
to reach the surface). Hard-excluding all 333 is too blunt: if the required partner is itself
broadly expressed, a standard cell line already supplies it and the gene traffics fine. Only
genes whose partner is *rare / restricted* truly need a special cell context.

## Inputs
- **Required genes**: `filters.co_receptor_dependency == 'required'` (333 of the 5,130 cohort).
- **Partners per gene**: `accessibility_risks.co_receptor_requirements.partners` — a list of
  partner gene symbols (e.g. `CD79B ← [CD79A, IgM, IgG]`). Fetched per required gene from
  `/v1/genes/{SYMBOL}` → `coreceptor_partners.jsonl`.
- **Partner expression**: deterministic HPA cell-line TPM (`max_tpm_hpa_celline`) from the
  genome-wide TPM table, joined **via the gene resolver on Ensembl id**
  (see `external/expression_tpm/PROVENANCE.md`).

## Historical ubiquity call (max-based prototype)
A partner is **common** iff its **HPA cell-line max TPM ≥ 10** (present in cancer cell lines).
This *replaces* the earlier heuristic (a hand-curated housekeeping set + the LLM
`expression_breadth` label) with a measured, reproducible signal.

## Rescue rule
```
coreceptor_ok(gene):
    if co_receptor_dependency != 'required':          return True
    partners = co_receptor_requirements.partners
    known = [partner_common(p) for p in partners if celline_tpm(p) is not None]
    return bool(known) and all(known)     # rescue iff EVERY expression-known partner is common
```
- all partners common → eligible for T1/T2.
- ≥1 rare partner (or no expression-known partner) → stays excluded / flagged.

## Result (max cell-line ≥ 10)
- **Rescued: 196** — flippases (`ATP11A/C ← CDC50A/TMEM30A`), Na/K-ATPase (`ATP1A/B`), `ABCC8/9`.
- **Kept out: ~137** — `ABCG5↔ABCG8` (gut/liver), `CD79B ← CD79A` (B-cell), TCR `← CD3`.

## Known limitations / knobs
- **max vs median across cell lines** — "max ≥10" = expressed in ≥1 cell line (permissive).
  A **median ≥10** (expressed in ≥half of cell lines) is a stronger "ubiquitous" test and would
  tighten rescues to genuinely-broad partners. Needs the full HPA `rna_celline` matrix
  (per-cell-line), not just the `max` column we currently carry.
- **alternatives vs obligate complex** — the `partners` list mixes "need any one" (Na/K-ATPase
  β subunits) with "need all" (complexes). `all(common)` is correct for complexes and
  conservative for alternatives; here it still rescues the alternative cases because every
  alternative subunit is itself cell-line-expressed.
- **partner-phrase noise** — some partners are LLM free-text
  (`"ATP4B-like beta subunit (ATP1B1/FXYD family)"`) that don't resolve; treated as
  expression-unknown (excluded from the `all()`).
- **threshold** — 10 TPM is the team-standard detection floor; raise for a stricter bar.
