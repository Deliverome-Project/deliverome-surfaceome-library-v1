# Scoping: which surface isoforms to include (expression + features)

**Problem.** The library today counts a surface *isoform* by **topology alone** (`final_library.n_surface_isoforms`: is the spliced variant still surface-attached?). That over-counts what we'd actually synthesize — many surface-topology-competent isoforms are minor/non-expressed or redundant with the canonical. We need a **selection/prioritization** layer.

## Features to prioritize on

**A. Expression (is the isoform actually made?)**
- **GTEx v10 transcript-level TPM** — per-`ENST` quantifications across ~50 tissues. Isoform-level TPM is noisier than gene-level, so use with a dominance ratio (isoform TPM / gene TPM) rather than an absolute cut.
- **APPRIS principal isoform** (`PRINCIPAL:1…5`) — the computationally-supported dominant transcript per gene; cheap, high-signal.
- **MANE Select** (`ENST` per gene) — the agreed canonical clinical transcript; usually = our canonical, so it's the *anchor*, not a candidate to add.
- **GENCODE / Ensembl Transcript Support Level (TSL)** — deprioritize TSL 4–5 / unsupported models.
- **HPA transcript data** where available; **recount3** as an alternative uniform GTEx/TCGA transcript quantification source.

**B. Non-expression features (already in our catalog per isoform)**
- **ECD residues** (`per_residue_topology` "O" count) — does the isoform keep a targetable ectodomain, or splice it out?
- **ECD % identity to canonical** (`ecd_pct_identity_to_canonical`) — high identity ⇒ redundant epitope (don't bother); low ⇒ a **distinct targetable surface** worth adding.
- **Signal peptide / TM count** changes vs canonical (topology-flipping splices).
- **Tissue / tumor specificity** of the isoform (an isoform restricted to tumor is more valuable than a ubiquitous one).

**C. GPCRs specifically — direct experimental data (Tedman et al. 2025)**
The Willow / Tedman GPCR surface screen (`s3://deliverome-external/tedman-gpcr-surface-screen/2025-09-19/`, bioRxiv 2025.09.19.677468) measured **per-isoform cell-surface expression** (anti-HA immunostaining) for 174 GPCR isoform transcripts, plus **"% change in PME" vs the canonical**. This is a *direct* surface readout — the gold standard for the GPCR isoform class, and it already showed **76% of GPCR isoforms lose surface** and that topology (not % identity) predicts it. We have the plumbing: `willow_enst_uniprot_isoform_map.json` (ENST→UniProt-isoform) and `willow_isoform_concordance.csv`. For GPCRs, use the screen's surface value directly instead of an expression proxy.

## The mapping problem (do this first)
Our catalog isoforms are keyed by **UniProt isoform ID**; expression is keyed by **Ensembl `ENST`**. Need a robust `ENST ↔ UniProt-isoform` map (extend the Willow one to all genes, via UniProt's `Ensembl_TRS` xrefs / MANE). **~half of screened isoforms were Ensembl-only** (no UniProt isoform) — a real annotation gap: if isoform breadth matters, source isoforms from Ensembl/GENCODE + run per-transcript DeepTMHMM rather than relying on UniProt isoforms only.

## Concrete data pulls
1. **APPRIS principal-isoform table** (human, GRCh38) — small, static; gives the dominant transcript per gene immediately.
2. **GTEx v10 transcript TPM matrix** (or recount3) — for the isoform TPM / dominance ratio. Large; pull per-gene for library genes rather than genome-wide.
3. **MANE Select + Ensembl TSL** (GENCODE) — transcript-quality metadata.
4. **UniProt↔Ensembl transcript xref** — to bridge our UniProt-isoform catalog to the expression `ENST`s.
5. **(GPCRs) already in hand** — Tedman Isoforms sheet; just join through the existing ENST↔UniProt map.

## Proposed rule (draft, to iterate)
Include an isoform if: surface-topology-competent (current gate) **AND** (isoform TPM / gene TPM ≥ threshold **OR** APPRIS principal **OR** — for GPCRs — Tedman surface ≥ 1000) **AND** it adds a distinct targetable surface (ECD present; ECD %-identity to canonical below a redundancy cut). Tune thresholds against the Tedman GPCR ground truth first, then generalize.
