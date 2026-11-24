# How good is the processed-form flag, and what else could supply it

`check_processed_forms.py` flags library members whose full-length ORF is not
the species on the cell. Source is UniProt, via one streamed bulk query for
`ft_chain`, `ft_site` and `ft_propep` over all reviewed human entries.

## It is 41% measured, 59% inferred

Evidence codes on the Chain and Site annotations backing the 185 flagged
library members:

| evidence | share |
|---|---|
| experimental (ECO:0000269) | 40% |
| by similarity (ECO:0000250) | 26% |
| sequence analysis, i.e. predicted (ECO:0000255) | 17% |
| curator inference (ECO:0000305) | 16% |
| large-scale / combinatorial | 1% |

So a flag is a reason to look, not a finding. Before a gene is truncated on
this basis its own annotation should be read — an ECO:0000269 cleavage with a
named protease is a different proposition from one inherited by similarity.

It also **under-calls**: PTPRF is furin-processed into E- and P-subunits with no
annotated cleavage Site, and MUC16's Golgi cleavage is not annotated as one. A
negative means "not annotated", never "not processed".

Not every miss is an error. CELSR1 comes back clean because it genuinely is
cleavage-deficient — its catalytic threonine is an alanine.

## A bare furin motif does not work — measured, not assumed

The obvious cheap addition is to scan for the proprotein-convertase consensus
R-X-K/R-R. Tested against the whole library:

| | |
|---|---|
| genes matching the motif | 570 (24% of the library) |
| genes UniProt calls processed | 185 |
| both | 62 |
| **precision of the bare motif** | **10.9%** |

A four-residue motif with two wildcard-ish positions occurs by chance in a
quarter of all surface proteins. It would add 508 false flags to find 62 real
ones. **Do not use a bare motif scan.** Anything sequence-based has to score
flanking context, not just match the consensus.

## What would actually improve on UniProt

1. **Experimental N-terminomics** is the right answer, because it is evidence
   rather than prediction: a protein whose observed N-terminus differs from its
   annotated start *is* processed. TopFIND aggregates these datasets.
2. **Our own PeptideAtlas data**, already in this pipeline for the
   reference-selection evidence rung. Semi-tryptic N-terminal peptides at an
   internal position are the same signal, on data we already hold.
3. **Context-scoring predictors** — ProP for proprotein convertases, and the
   ML protease-site tools (PROSPER, DeepCleave, Procleave). These weigh flanking
   residues rather than matching a consensus, which is the failure the motif test
   above demonstrates.
4. **MEROPS** for curated protease–substrate pairs where the protease is known.

Of these, (2) is the cheapest by a wide margin and uses data already on S3.

## Should the short ones be expressed full-length and left to the cell?

Yes, by default. The cell's furin, GPS and SEA machinery knows the boundary
exactly and we would be guessing it. 173 of the 185 processed members fit the
packaging budget, so this is the common case, not the exception.

Actively removing a cleavage site is worse than leaving it. For polycystin-1 a
cleavage-dead mutant **never reaches the surface at all**, even with PC2, and
MUC16's SEA cleavage depends on flanking context that over-trimming destroys.
Cleavage is often part of the trafficking route, not an obstacle to it.

But "let the cell cleave it" has one consequence that has to be designed for,
and it is a tag-placement consequence rather than a construct one.

### The cut-off piece does not always stay

Of the 173 processed members we would express full-length:

| what happens to the cut-off piece | n |
|---|---|
| protease not named | 102 |
| other | 30 |
| stays, non-covalently (furin, autolysis) | 29 |
| **leaves (sheddase)** | **12** |

For the furin and autolysis group the two pieces remain associated — the PTPR
E- and P-subunits, the NOTCH S1 heterodimer, the adhesion-GPCR NTF and CTF — so
a tag on the distal piece is still on the cell, if non-covalently.

For the sheddase group it leaves. CHL1 (ADAM8), PTK7 (MMP14), EPHB2, DAG1
(MMP9), CDH1, APP, IL6R (ADAM10/17), FASLG, FCGR3A (ADAM17), TNF. These are the
textbook shed ectodomains — soluble IL6R, soluble TNF, sAPPα, soluble
E-cadherin — which is a good sign the classification is working.

Twelve is a floor, not a count: 102 have no named protease.

### Three concerns, in order

1. **An N-terminal tag on a sheddase substrate leaves with the ectodomain.** The
   surface stain reads absent while the protein is present and behaving
   correctly. Indistinguishable from a failed construct unless designed for.
2. **Shedding is regulated, so the readout moves with cell state.** ADAM17
   responds to signalling; EGFR activity drives LAR shedding; TNFα drives MUC17.
   The same construct reads differently across lines — a confound that looks
   like biology in a pooled screen.
3. **Cleavage is not quantitative**, so a single cell carries both cleaved and
   uncleaved species. Heterogeneity nobody designed in.

### What follows

Do not remove cleavage sites. Make tag placement cleavage-aware instead: on a
sheddase substrate the tag belongs membrane-proximal to the cut or on the
cytoplasmic side, never on the distal ectodomain. That is a constraint on the
tag-placement decision that is still open, and it applies to roughly 12 members
known and up to 102 more unresolved — so the protease for those 102 is worth
settling before tag positions are fixed.
