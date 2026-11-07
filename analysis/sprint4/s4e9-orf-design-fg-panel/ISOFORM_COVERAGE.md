# Isoform cohort coverage

All 451 current library genes with a non-canonical reviewed UniProt isoform but no
selected isoform cohort member occur in the scored isoform universe. They are
excluded by the existing evidence and surface-anchoring rules. April query-set
membership is not the boundary of this universe; the 92 genes outside that set
also occur in the scored universe. The current panels and constructs are unchanged.

## Selection and accounting

The universe is the set of distinct UniProt and Ensembl protein sequences offered
for each gene in the committed library list. The cohort is the subset of alternate
forms passing the evidence and anchoring rules. The audit applies the historical
cohort's filter precedence and records the input evidence for every form, including
reference forms, so no gene is assigned a single reason that hides other exclusions.

| First applicable decision | Forms among the 451 genes | Historical rule |
|---|---:|---|
| Reference | 450 | The reference is not an alternate cohort candidate |
| Insufficient evidence | 2,400 | No diagnostic peptide, no positive OpenProt translation score, and no APPRIS TRIFID score above 0.5 |
| Loses all reference transmembrane helices | 30 | Historical relative topology filter |
| Lost surface anchoring | 16 | Explicit historical topology verdict |
| Soluble restriction | 15 | Evidence-passing form with a non-surface topology label and no retained GPI anchor |
| Unexplained | 0 | An unexplained eligible form makes the audit fail |

399 genes have only insufficient-evidence alternate forms. The remaining 52 have
one or more topology exclusions, sometimes alongside insufficient-evidence forms.
These counts explain selection; insufficient evidence does not establish that a
protein isoform is absent or nonfunctional. The audit preserves the existing
criterion, including its prediction-based TRIFID cutoff, rather than choosing a
new evidence threshold.

`isoform_coverage_genes.tsv` records every gene and counts all applicable form
reasons. `isoform_coverage_forms.tsv` contains all 2,911 scored rows for those genes,
their evidence columns, and the historical verdict/restriction supporting each
first applicable decision. `isoform_coverage_summary.json` records counts,
source-manifest hashes, script versions, output hashes and replay limits.

## Reference enumeration limitation

LILRB1 has 17 scored forms but no row marked as its reference. The frozen April
cross-reference maps Q8NHL6 only to ENSG00000274669, ENSG00000277134 and
ENSG00000277807. APPRIS 115 instead supplies primary gene ENSG00000104972 and MANE
transcript ENST00000324602, encoding 652 residues. The historical builder uses the
nonempty cross-reference set to enumerate Ensembl forms, suppressing its
APPRIS-by-symbol fallback; reference selection does fall back by symbol, so the
chosen MANE sequence is absent from the enumerated forms. UniProt Q8NHL6-3 also has
652 residues but differs from that MANE sequence at four positions.

All 17 enumerated LILRB1 forms fail the existing evidence criterion. Their cohort
absence is accounted for, but reference enumeration across releases remains an
upstream followup. The audit flags genes without reference rows and does not repair
this mapping or silently add forms to the cohort.

## Historical topology limitation

The producing scripts join some verdicts by gene, length and accession and resolve
some Ensembl-only sequence ties by unrecorded hash order. The frozen replay retains
the committed topology shortlist. In this run, 94 replayed shortlist sequences are
a same-length sibling of the committed sequence; all 1,894 cohort rows and their
gene order reproduce. The ledger explains historical decisions and does not
establish exact-sequence topology for a new construct. A future rebuild must use
sequence hashes before applying those predictions to new candidates.

## Reproduction and provenance

```bash
uv sync --locked
uv run --locked python analysis/sprint4/s4e9-orf-design-fg-panel/reproduce_isoform_universe.py \
  --step isoforms --work /tmp/l1-isoforms
uv run --locked python analysis/sprint4/s4e9-orf-design-fg-panel/audit_isoform_coverage.py \
  --replay /tmp/l1-isoforms/iso --out /tmp/l1-isoform-coverage
uv run --locked pytest tests/test_isoform_coverage.py -q
```

| Element | Source |
|---|---|
| Audit code and tests | [Audit script](https://github.com/Deliverome-Project/deliverome-analysis/blob/f008eb9175e33c80a38bb3847d7376764c56c11f/analysis/sprint4/s4e9-orf-design-fg-panel/audit_isoform_coverage.py), [filter precedence and failure tests](https://github.com/Deliverome-Project/deliverome-analysis/blob/f008eb9175e33c80a38bb3847d7376764c56c11f/tests/test_isoform_coverage.py) |
| Historical scripts and immutable inputs | [Replay script and producing commits](https://github.com/Deliverome-Project/deliverome-analysis/blob/f008eb9175e33c80a38bb3847d7376764c56c11f/analysis/sprint4/s4e9-orf-design-fg-panel/reproduce_isoform_universe.py), [versioned S3 keys and SHA-256 manifest](https://github.com/Deliverome-Project/deliverome-analysis/blob/f008eb9175e33c80a38bb3847d7376764c56c11f/analysis/sprint4/s4e9-orf-design-fg-panel/reproduce_isoform_universe_inputs.json) |
| Current gene list and recorded cohort | [Pinned library input manifest](https://github.com/Deliverome-Project/deliverome-analysis/blob/f008eb9175e33c80a38bb3847d7376764c56c11f/analysis/sprint4/s4e9-orf-design-fg-panel/library_pipeline_inputs.json), plus committed `gene_list_additions.tsv` |
| Durable audit bundle | `processed://sprint4/s4e9-orf-design-fg-panel/library-design/isoform-coverage-audit/2026-11-05-f008eb91/` (create-only; code, tests, manifests and derived ledgers) |
| Pull request | [Isoform coverage audit](https://github.com/Deliverome-Project/deliverome-analysis/pull/257), dependent on [ortholog rebuild](https://github.com/Deliverome-Project/deliverome-analysis/pull/234) |
| Review question | [Library review issue 18](https://github.com/Deliverome-Project/l1-library-review/issues/18) |

## Decision log

| Date | Was | Now | Why |
|---|---|---|---|
| 2026-11-05 | Cohort absence described as universe absence; 92 genes attributed to April query-set exclusion | All 451 genes traced through their original universe rows and cohort decisions | The builder iterates its committed library list. April membership supplies topology/cohort annotations, not a gene-eligibility gate. No cohort or panel change follows from this audit. |
