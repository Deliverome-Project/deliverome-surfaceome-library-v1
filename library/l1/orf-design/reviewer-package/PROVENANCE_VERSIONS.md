# Software, Predictor, and Database Provenance Versions

This file consolidates the exact runtime environment, non-replayed upstream predictor builds,
model weight checksums, and reference database releases used to design and verify the
Deliverome v1 surface-protein ORF library (order of record `2026-11-25-pre-freeze-e`).

All **3,591** frozen input files used by the offline replay pipeline (`run_review.py`, 30/30
stages + `brief_numbers.json`) are SHA-256-pinned in [`bundle_manifest.tsv`](bundle_manifest.tsv).
The licensed/GPU predictors below are **not** re-executed during offline replay — their outputs
are frozen in the bundle — and are documented here so future researchers can reproduce or extend
upstream predictions on new targets.

---

## 1. Offline Replay Runtime (Container & Python)

| Component | Exact Pinned Version / Digest | Source |
|---|---|---|
| Base OS / Python image | `python:3.14.10-slim-bookworm@sha256:3fe1a0b6c2d6b2ebb572906cdabeca75aeb558c0c17db5ffefbcf8094871c32b` | `Dockerfile` |
| `uv` package manager | `ghcr.io/astral-sh/uv:0.13.6@sha256:af22a52f5eafef3d78e551599d405c02df61919d18028f0cea8eda93bc9986d0` | `Dockerfile` |
| Python runtime | `3.14.10` (`PYTHONHASHSEED=0`) | `.python-version` |
| Python dependency lockfile | `uv.lock` (all runtime packages pinned by hash; `dev`/`flow`/`modal` groups excluded from replay container) | `uv.lock` |

---

## 2. Structure, Topology, Signal Peptide, and Disorder Predictors

| Tool | Exact Version & Runtime | Settings / Notes | Weight / Package Provenance |
|---|---|---|---|
| **ColabFold** | `1.6.2` (`ghcr.io/sokrypton/colabfold@sha256:22f1243a8feca427b6942889d4dc3b549a3b75b23df9cf85c5f512a85e4f7f77`), `alphafold-colabfold 2.3.18`, `jax 0.10.2`, `jaxlib 0.10.2`, `dm-haiku 0.0.16`, `numpy 2.5.1` | `alphafold2_ptm`, `model_1`, `num_recycle=3`, unrelaxed, no templates; MSAs from `api.colabfold.com` (`mmseqs2_uniref_env`, `use_env=True`, `use_filter=True`), cached per native protein SHA-256 with cassette gap columns | All 12 `alphafold_params_2021-07-14` weight files SHA-256-pinned in [`tag_refold/colabfold_provenance.json`](../tag_refold/colabfold_provenance.json); all 3,352 tagged PDBs hashed in [`tag_refold_structures.tsv`](../tag_refold_structures.tsv) |
| **DeepTMHMM** | `1.0.24` (`deeptmhmm-1.0-academic`, DTU, license-gated) | Executed on exact built protein sequences on Modal GPU volume | Model checkpoints (`deeptmhmm_cv_0..4.model`) SHA-256-pinned in [`library_pipeline_inputs_saved.json`](../library_pipeline_inputs_saved.json); predictions frozen in pinned topology parquet |
| **SignalP 6** | `6.0+h` (`torch 1.13.1+cu117`, DTU, license-gated) | Slow/accurate eukaryotic mode (`--organism eukarya`); run on native and candidate tagged N-termini | All CPU/GPU `.pt` checkpoints and package files SHA-256-pinned in [`nterm/signalp6_identity.json`](nterm/signalp6_identity.json); outputs frozen in [`nterm/nterm_bundle_manifest.tsv`](nterm/nterm_bundle_manifest.tsv) |
| **NetGPI** | `1.1` (DTU web service) | Executed on candidate human, isoform, and ortholog sequences | Raw HTTP exchanges and parsed outputs archived with SHA-256 in [`library_pipeline_inputs_saved.tsv`](../library_pipeline_inputs_saved.tsv) |
| **metapredict** | `3.0.2` | Per-residue disorder tracks for human/ortholog internal-loop ranking and N-after-SP onset | Outputs frozen in committed feature tables & D1 snapshot |
| **NetSurfP-3.0** | Snapshot `dis_2026_09_28` (via `surfaceome_public` D1) | Per-residue RSA and disorder alongside `metapredict` and AlphaFold `1 - pLDDT` | Frozen in `surfaceome_public_s5e13.sqlite.gz` and `nterm_bundle_manifest.tsv` |
| **ProDomino** | `v8_prodomino` internal-site ranker | Combined with AlphaFold `-plddt_w5`, empirical Bayes penalties, and contact/disulfide/glycan vetoes | Per-residue predictions and benchmark tables bundled in `bundle_manifest.tsv` |
| **TmAlphaFold / TMDET** | Fetched `2026-11-14` | Membrane plane placement for extracellular membrane-side gate | All fetched models & placements bundled in `membrane-side-gate/inputs/2026-11-14/inputs_bundle.tar.gz` (SHA-256 `90215912...`) |

---

## 3. Pinned External Reference Databases

| Database / Resource | Release / Version | Role in v1 Design | Bundle / Archive Status |
|---|---|---|---|
| **accessible-surfaceome** | `v1.3.0` (commit `f9c73691`, Zenodo `10.5281/zenodo...`) | Candidate surfaceome catalog, evidence tiers, topology & partner annotations | Bundled whole (`bundle_manifest.tsv`) & public R2 |
| **Ensembl / Compara** | Releases `115` and `116` | MANE Select GTF, human CDS/isoforms, mouse & cyno 1:1 orthologs | Bundled whole & public R2 |
| **UniProtKB** | `2026_03` (catalog/giants) & `2026_04` (C-term/propeptide/topology reviews) | Canonical sequences, signal/propeptide/GPI/shedding features, LoF restorations (`OR52P1` `Q8NH57`) | Bundled whole & public R2 |
| **AlphaFold DB** | `v6` (`model_v6` PDBs) | Native reference structures (`2,847` models archived) | Write-once `external://alphafold-db/v6/` + public R2 |
| **APPRIS** | `2026_08.v50` (Ensembl 115) | Principal isoform selection ladder (CC BY-NC-ND 4.0, redistributed unmodified) | Bundled whole & public R2 |
| **PeptideAtlas** | Build `607` (`queries-2026-11-26-pre-freeze-e`) | Observed peptide support for reference proteoforms and isoforms | Bundled whole & repo `data/` |
| **OpenProt** | `2.2` (`min_2_pep-uniprot2022_06_01`) | Isoform proteomics support | Bundled in release tarballs |
| **Pfam** | `39.0` (`Pfam-A.hmm.gz`) | Domain boundaries & structural features | Identical copy at EBI FTP (`public_source`) + bundled in release tarballs |
| **Human Protein Atlas** | `v23` | Cell-line and tissue RNA consensus TPM | Bundled whole & public R2 |
| **GTEx** | `v10` (`RNASeQCv2.4.2` gene median TPM) | Tissue expression filter & compiled expression join | Bundled whole & lab S3 |
| **TCGA / TARGET / GTEx Toil** | UCSC Xena `TcgaTargetGtex_rsem_gene_tpm.gz` (`1.32 GB`, SHA-256 `a8c36cb1...`) | Compiled expression join (`genomewide_expression_tpm.csv`) | Compiled table bundled & committed in repo `data/`; raw 1.32 GB matrix archived on lab S3 |
| **ARCHS4** | `v2.5` (`2026-09-04` 42-line consensus) | Screening-line co-receptor partner expression rule | Bundled & committed in repo `data/` |
| **TOPDB** | `v3.3` (fetched `2026-11-20`) | Experimental membrane topology validation | Bundled whole |
| **Tedman et al. (2026)** | Published GPCR surface-display screen table | Independent display benchmark for OR/GPCR and accessory-factor families | Bundled in release tarballs |

Full per-script external input ledger (all 46 external sources with SHA-256 checksums): [`external_inputs.tsv`](../external_inputs.tsv).
