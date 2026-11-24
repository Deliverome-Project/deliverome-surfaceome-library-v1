# Willow's GPCR surface-expression screen (Tedman et al. 2025)

External dataset — experimental **ground truth** for GPCR (and GPCR-isoform) cell-surface
expression, used to validate the surface-library filter (s1e3).

## Source
- **Paper:** Tedman, A. et al. "Deep receptor scanning..." bioRxiv 2025.09.19.677468 (Coyote-Maestas lab).
  https://www.biorxiv.org/content/10.1101/2025.09.19.677468v1.full
  DOI: 10.1101/2025.09.19.677468
- **File:** Supplementary `media-2.xlsx` (downloaded by Becca; local copy
  `tedman2025_gpcr_screen_media-2.xlsx`).
- **Accessed:** 2026-08-14. sha256 `c0fbcbd1df5ab488e026239fb55aa92a7a8d9171154bec20bc837f6fb9685286`.

## Assay (what the numbers mean)
A pooled surface-display screen of GPCR ORFs in a cultured cell line: constructs are
FACS-sorted on cell-surface immunostaining and read out by RNAseq/riboseq. The key
per-construct readout is **surface immunostaining intensity** (higher = more protein at the
plasma membrane). Screen is done **without** olfactory RTP chaperones, so a high value for an
olfactory receptor means it reaches the surface autonomously.

## Sheets
- **Canonical** (766 constructs, 763 unique `gene_id`) — one row per canonical GPCR ORF.
  Key cols: `Receptor Name` (`GENE_ENST…[_SP_HA]`), `gene_id` (symbol), `GPCR Class`
  (A / A (Olfactory) / B / C / F / Taste / Vomeronasal), **`Immunostaining Intensity`** (+ SD)
  = surface readout, `Low GFP Enrichment`, `gfp z score`, RNAseq/riboseq/TE, UMAP coords.
- **Isoforms** (174 transcripts, 105 genes) — one row per non-canonical isoform.
  Key cols: `Isoform Transcript` (`GENE_ENST…`), **`Surface Expression`** (+ SD),
  `Canonical Surface Expression` (+ SD), `length (AA)`, `canonical length`,
  `Difference in Isoform PME`, `% change in PME` (PME = plasma-membrane expression).
- **MACHINE_LEARN_PROTEIN_FEATURE** (779 transcripts) — per-transcript `GPCR_PME` + ~75
  sequence/codon/structure features used for the paper's ML model.

## ID note (matters for joining)
Screen uses **Ensembl transcript IDs** (`ENST…`); our catalog keys isoforms by **UniProt
isoform IDs** (`P35372-10`). Canonical rows join to our cohort cleanly on **gene symbol**;
per-isoform matching is approximate (gene + protein length) unless an ENST↔UniProt map is built.

## S3 location (canonical copy)
Uploaded to the immutable third-party bucket following the `{source}/{version}/` convention:

- **`s3://deliverome-external/tedman-gpcr-surface-screen/2025-09-19/`**
  - `tedman2025_gpcr_screen_media-2.xlsx` — the source supplementary data
  - `PROVENANCE.md` — this file

Version `2025-09-19` = the bioRxiv **v1** posting date. Uploaded 2026-08-14 via the
`becca-admin` IAM profile (scientist creds are read-only on `deliverome-external`).

Read it back through the helper (read-only on external):

```python
from deliverome_analysis.s3 import external
p = external.download("tedman-gpcr-surface-screen/2025-09-19/tedman2025_gpcr_screen_media-2.xlsx")
```

The copy under `scratch/becca/surface-library-filters/external/willow_gpcr_screen/` is a local
working mirror (not committed — data-file policy); S3 is the canonical source.
