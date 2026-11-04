# Deterministic expression TPM — provenance & processing

## What this is
Per-gene **maximum expression (TPM/nTPM)** across four human reference datasets,
used as a *deterministic* expression filter for the surface-protein library
(alongside the accessible-surfaceome LLM `low_endogenous_expression` flag — used
**in addition to**, per Becca's instruction).

- File: `genomewide_expression_tpm.csv`
- Columns: `gene_symbol_canonical, ensembl_bases, max_tpm_gtex,
  max_ntpm_hpa_tissue, max_tpm_hpa_celline, max_tpm_tcga, max_expression_any,
  passes_expression_threshold, in_surfaceome_expressed, is_protein_coding_hgnc`
- Coverage: 18,457 HGNC-approved protein-coding genes reaching ≥10 TPM/nTPM in ≥1 dataset.

## Re-derived in-repo (2026-08-04) — `genomewide_expression_tpm_rederived.csv`
A **fresh, independent re-derivation** from raw sources now exists alongside the reused
table, produced by `pipeline/build_expression_reference.py` (this repo). It adds a
**`median_tpm_hpa_celline`** column (median TPM across all 1,206 HPA cell lines — the
ubiquity signal for the co-receptor rescue). Validation: **65,098 genes evaluated,
35,105 reach ≥10 TPM/nTPM** — matching deliverome-internal's documented counts exactly,
confirming the port reproduces their pipeline. Per-source download traceability (URL,
sha256, bytes, timestamp) is in `download_traceability.json`. Raw source URLs (all
verified 200): GTEx `adult-gtex/bulk-gex/v10/rna-seq/…gene_median_tpm.gct.gz`;
HPA `proteinatlas.org/download/tsv/{rna_tissue_consensus,rna_celline}.tsv.zip`;
TCGA `toil-xena-hub…/TcgaTargetGtex_rsem_gene_tpm.gz`; HGNC `…/hgnc_complete_set.txt`.

## Source of the reused file — REUSE, not re-derivation (documented chain)
Fetched verbatim from the Deliverome team's version-pinned pipeline in
**`Deliverome-Project/deliverome-internal`** (avoids version drift vs. re-downloading raw):
- Output reused: `data/analysis/genomewide_expression_filtering/genomewide_expressed_genes_protein_coding_collapsed.csv`
- Build script: `data/analysis/genomewide_expression_filtering/build_genomewide_expression_list.py`
- Repo commit (pinned): `5a1416eab0e1e16f3e7fa66ae539f7b2fe693730` (branch main)
- Retrieved via Git LFS on **2026-08-04**.

To instead establish an *independent* copy in this repo's `deliverome-external` S3
bucket, re-download the raw sources below and re-run the build script — this needs
write access to `deliverome-external` (admin step; not the scientist read-only role).

## Raw sources + versions (as consumed by the build script)
| Dataset | File | Content | Per-gene aggregation |
|---|---|---|---|
| **GTEx v10** | `GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct` | gene **median TPM**, 59,033 genes × 68 normal tissues | max across tissues → `max_tpm_gtex` |
| **HPA Tissue Consensus** | `hpa_rna_tissue_consensus.tsv` | consensus (HPA/GTEx/FANTOM5) **nTPM**, 20,151 genes × 51 tissues | max nTPM per gene → `max_ntpm_hpa_tissue` |
| **HPA Cell Lines** | `hpa_rna_celline.tsv` | **TPM**, 20,151 genes × 1,206 cancer cell lines | max TPM per gene → `max_tpm_hpa_celline` |
| **TCGA primary tumors** | `TcgaTargetGtex_rsem_gene_tpm.gz` (UCSC Xena "TCGA TARGET GTEx", RSEM gene TPM, values are `log2(TPM+0.001)`) | TCGA-* samples only | max over TCGA samples, then `2^max_log − 0.001` → `max_tpm_tcga` |
| **HGNC** | `hgnc_complete_set.tsv` | approved protein-coding genes + alias/prev-symbol + Ensembl | protein-coding filter, Ensembl-base collapse, symbol canonicalization |

Portals (versions as above): GTEx v10 — gtexportal.org/home/downloads ; HPA —
proteinatlas.org/about/download (`rna_tissue_consensus.tsv.zip`, `rna_celline.tsv.zip`) ;
TCGA — xenabrowser.net (TcgaTargetGtex) ; HGNC — genenames.org.

## Derived field + threshold
```
max_expression_any = max(max_tpm_gtex, max_ntpm_hpa_tissue, max_tpm_hpa_celline, max_tpm_tcga)
passes_expression_threshold = (max_expression_any >= 10)   # TPM/nTPM; deliverome-internal standard
```

## Units & caveats (read before using)
- **Mixed units**: GTEx / HPA-cellline / TCGA are TPM; HPA tissue is nTPM (normalized).
  `max_expression_any` takes the max across these — not strictly unit-homogeneous, but
  matches the established team pipeline. Documented, not silently mixed.
- **≥10-only file**: this collapsed table contains only genes reaching ≥10 in ≥1 dataset;
  a gene absent from it is <10 everywhere (fails the ≥10 filter). Exact sub-10 values not retained.
- **"Expressed somewhere" is permissive**: max over 68 tissues + 1,206 cell lines + TCGA
  → ~96% of protein-coding genes qualify. Read as **"not truly absent,"** NOT "ubiquitous."
  For a "ubiquitous co-receptor" call, use a *higher* threshold and/or the cell-line column,
  not bare ≥10.

## How it's joined here
Cohort genes carry accessible-surfaceome HGNC symbols (which can be aliases, e.g.
`CDC50A` vs approved `TMEM30A`). We map cohort symbols → HGNC-approved symbols via
`hgnc_complete_set.tsv` (alias_symbol + prev_symbol) before joining to
`gene_symbol_canonical`, so alias drift doesn't cause false "low-expression" calls.
