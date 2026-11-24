# Licensing and redistribution

This repository holds the design pipeline, construct tables, verification tests,
and external review brief for the Deliverome v1 surface-protein ORF library (L1
in code).

- **Code** (Python scripts, library utilities, tests, and workflow definitions):
  licensed under the **MIT License** (see [`LICENSE`](LICENSE)).
- **Derived library designs, tables, sequences, brief, and figures** created by
  Deliverome Bio, Inc.: dedicated to the public domain under
  **[CC0 1.0 Universal](https://creativecommons.org/publicdomain/zero/1.0/)**,
  to the extent we hold rights in them. This covers our gene-gate tables,
  proteoform/panel selections, tag-placement calls, construct and barcode
  tables, GenBank maps, design brief, and figures, so both academic and
  commercial users can reuse them without legal friction. Please still cite the
  repository using [`CITATION.cff`](CITATION.cff) when using these designs in
  published work.

## Upstream data sources and redistribution

The library design integrates public databases, published datasets, and
third-party structure/topology predictors. Every input file's source, home, and
per-file license is recorded in
[`data_locations.tsv`](library/l1/orf-design/reviewer-package/data_locations.tsv).
Upstream terms apply to any original upstream file redistributed for offline
reproducibility:

| Source | Used for | Upstream license / terms | Where it lives | Attribution |
|---|---|---|---|---|
| **accessible-surfaceome v1.3.0** | Surfaceome gene universe, evidence tiers, deep-dive flags | CC0 1.0 (derived data, Zenodo `10.5281/zenodo.23445521`); MIT (code) | Public R2 bucket & `data/` | Carlson, Hollingsworth & Smallegan, 2026 (bioRxiv `10.64898/2026.10.18.759859`) |
| **UniProtKB / Swiss-Prot** | Reference sequences, isoforms, topology & binding features | [CC BY 4.0](https://www.uniprot.org/help/license) / EMBL-EBI terms | Public R2 bucket & derived tables | UniProt Consortium |
| **Ensembl (115 / 116) & Compara** | Coding sequences, MANE Select GTF, mouse & cynomolgus orthologs | Open ([EMBL-EBI terms](https://www.ensembl.org/info/about/legal/disclaimer.html), no redistribution restrictions; citation requested) | Public R2 bucket & derived tables | EMBL-EBI / Ensembl |
| **Pfam 39.0** | Domain annotations for ortholog correspondence | CC0 / EMBL-EBI terms | EBI FTP / public source | Paysan-Lafosse et al. / InterPro–Pfam |
| **AlphaFold DB** | 3D protein structures for ranking, orientation, and refold comparisons | [CC BY 4.0](https://alphafold.ebi.ac.uk/) | Public R2 bucket (unmodified models) | Jumper et al., 2021; Varadi et al., 2022; Fleming et al., 2025 |
| **Human Protein Atlas (HPA)** | Expression & subcellular localization | [CC BY 4.0](https://www.proteinatlas.org/about/licence) | Public R2 bucket & derived flags | Uhlén et al., 2015; Thul et al., 2017 |
| **APPRIS (`e115v50` / `e116v50`)** | Principal isoform annotations & TRIFID scores | [CC BY-NC-ND 4.0](https://appris.bioinfo.cnio.es/) (non-commercial) | Release tarball / public mirror unmodified; derived calls only in library tables | Rodriguez et al., 2022; Pozo et al., 2021 |
| **DeepTMHMM v1, SignalP 6, NetGPI 1.1** | Transmembrane topology, signal peptide cleavage, GPI anchor predictions | DTU Health Tech academic / licensed software | Derived prediction outputs only (tools themselves are not redistributed) | Hallgren et al., 2022; Teufel et al., 2022; Gíslason et al., 2021 |
| **TmAlphaFold / TMDET** | Structure-based membrane-plane placement | Open academic resource ([CC BY 4.0](https://tmalphafold.ttk.hu/)) | Derived membrane-side flags only | Tusnády et al., 2005; Dobson et al., 2023 |
| **PeptideAtlas, OpenProt, GENCODE Ribo-seq, Tedman et al. 2026** | Proteoform evidence and GPCR display benchmarks | Source-specific / unclear redistribution terms | Release tarball only (not in the public R2 bucket); derived decisions only in library tables | Respective authors (cited in `BRIEF.md`) |

## Note on third-party redistribution

- If you redistribute only **our code and derived library designs** (the
  scripts, construct catalog, and brief), MIT (code) and CC0 1.0 (derived
  designs) apply.
- If you redistribute **upstream source bundles** from the public bucket or
  release tarballs, comply with the per-source licenses listed above and in
  [`data_locations.tsv`](library/l1/orf-design/reviewer-package/data_locations.tsv)
  (in particular, APPRIS carries a non-commercial CC BY-NC-ND 4.0 term, and
  AlphaFold DB / UniProt / HPA require attribution under CC BY 4.0).
