# Audit of the literature tag-site positives behind the S5E13 ranker

The S5E13 tag-site ranker (draft PR [#180](https://github.com/Deliverome-Project/deliverome-analysis/pull/180))
is trained on 176 published internal extracellular tag-insertion sites across 127 human surface proteins.
An LLM literature agent (prompt 1.4.0) retrieved them and published them to the public tag-site API as
`literature_retrieved` + `internal` rows; the checks on each row were the agent's own. This audit reads
every cited paper and asks whether it places a tag at the row's residue in UniProt canonical numbering of the
human protein. Every one of the 176 has a verdict; no paper was inaccessible, although 17 sites rest on an
abstract only (below).

A row says `insert_after_residue = N`: the tag sits between residues N and N+1. Every training row's
`residue_before` / `residue_after` matches the human canonical sequence, so the agent's flank check cannot
catch a wrong paper, a wrong gene or a mouse-numbered construct; that is what this audit adds.

## Tally

| Verdict | Meaning | Sites |
|---|---|---|
| confirmed_exact | paper states this position or pins it by flanking residues; human canonical numbering, or an ortholog with zero offset | 125 |
| confirmed_shifted | paper supports a specific position, offset from the row (species, isoform, mature numbering) | 23 |
| region_only | paper names a loop, domain or landmark, not a junction; the row pinned one anyway | 11 |
| not_supported | the paper tags a different gene, or shows the site is not extracellular | 7 |
| wrong_species_numbering | non-human (or unstated) numbering with no determinable human equivalent | 5 |
| terminal_misfiled | the tag is the N- or C-terminus of a truncated construct | 5 |
| paper_inaccessible | — | 0 |

By the agent's own evidence tier (its `evidence_type`):

| Tier | Sites | confirmed_exact | confirmed_shifted | region_only | not_supported | wrong_species | terminal |
|---|---|---|---|---|---|---|---|
| exact-site ("published tag insertion at this exact site") | 111 | 92 | 8 | 5 | 3 | 0 | 3 |
| same-loop ("in the same loop or domain") | 42 | 17 | 11 | 6 | 4 | 2 | 2 |
| tolerance ("transposon, FP fusion") | 23 | 16 | 4 | 0 | 0 | 3 | 0 |

By `position_evidence`: validated 125 sites, 102 confirmed_exact (82%); inferred 49, 23 confirmed_exact (47%);
topology_derived 2, neither confirmed. By `validation_level`: surface_and_function 62, 47 confirmed_exact;
surface_only 89, 66 confirmed_exact. Two of these counts moved after a second read of a cited construct paper and a preprint's oligo table (MMP14 111 and TMEM123 33, decision log below).

Species of the tagged construct: human 94, mouse 36, rat 19, unstated 20, pig 3, chicken 2, zebrafish 1,
*Xenopus* 1. Of the 125 confirmed_exact sites, 46 are mouse, rat or unstated constructs whose local alignment
to human has zero offset; 27 are loop replacements, residue substitutions or amber-codon sites rather than
pure insertions (the position is right, the construct is not an insertion). The 23 shifted sites are off by
−16 to +34 residues; the causes are mouse numbering (15), mature-protein numbering (3, ITGB7), an isoform
(1, NRG1 type III), rat (1), chicken (1), a junction read on the wrong side of a motif (1) and a tag placed
before rather than after the residue the text names (1, TMEM123).

The biology claims are overstated more often than the residue. Of 151 rows labeled `surface_only` or
`surface_and_function`, 36 have no non-permeabilized surface assay of the tagged construct in the cited paper
(10 show the tag is not at the surface or the construct fails to traffic; 26 report no surface readout). Of
67 rows claiming function, 18 have no functional measurement of the tagged construct (function was measured
on an untagged construct, on a different protein, or cited from elsewhere).

## What this means for the 176

| Use | Sites | Composition |
|---|---|---|
| Residue-level positives that survive | 148 | 125 confirmed_exact + 23 confirmed_shifted at the corrected residue |
| of which tolerated (gold standard) | 143 | minus the 5 failed sites below |
| of which exact-site tier | 100 | 92 exact + 8 shifted |
| Clean exact-site positives (no correction needed) | 92 | exact-site tier, confirmed_exact |
| Should become region labels | 11 | region_only; the loop or domain is supported, the residue is not; windows in `region_labels.tsv` |
| Should be dropped | 17 | 7 not_supported + 5 wrong_species_numbering + 5 terminal_misfiled; 8 of them are labeled surface_and_function |

So 28 of the 176 positives (16%) are wrong at the residue and 17 (10%) have no support for the gene at all.
The ranker's folds and the sibling evidence-tier retraining (PR
[#232](https://github.com/Deliverome-Project/deliverome-analysis/pull/232)) both use the agent's tiers; the
table above shows the "validated" tier is 82% clean and the "inferred" tier 45%, so the tier split carries
real signal but is not a clean split. The corrected label set is `gold_standard.tsv` (below); `region_labels.tsv` carries the loop-level sites.

## The three most consequential corrections

1. **All five ABCB1 positives are mouse Mdr3 numbering, +4 on human P-gp.** Kast et al. 1995 and 1996
   ([PMID 7535563](https://pubmed.ncbi.nlm.nih.gov/7535563/), [PMID 8621583](https://pubmed.ncbi.nlm.nih.gov/8621583/))
   inserted HA "at six different positions of the Mdr3 protein (101, 161, 206, 244, 320, and 376)". Mouse 101
   (..KRAMFA|KLEEEM..) aligns to human 105 (..DTGFFM|NLEEDM..); the same +4 holds at 320→324, 736→740, 741→745
   and 961→965. The rows carry `surface_and_function`, which the papers do support, on the mouse protein.
2. **Six positives cite a paper about a different gene.** ABCC2-163 and -574 cite the MRP1 (ABCC1) topology
   paper; TMEM225-125 is the PMP22 row duplicated under another gene; KCNK5-214 is a TASK-1 (KCNK3) construct;
   TMPRSS9-270 is secreted polyserase-3; GABRP-27 cites a paper that states "We excluded GABA A receptor
   subunit variants derived from 3 GABRs (GABRE, GABRP, GABRR2)". In every case the human flank matched by
   coincidence, which is why the agent's sequence check passed. A seventh not_supported site, TLCD1-201, is a
   chicken construct whose paper shows that loop is intracellular ("without detergent ... no signal is apparent").
3. **The validation tiers are inflated.** 36 of 151 surface-claimed rows lack a surface assay of the tagged
   construct and 18 of 67 function-claimed rows lack a functional one. Examples: LRRC8A-76 is a BRIL fusion
   expressed in Sf9 cells for cryo-EM; HCN2-311's HA tag is never used in any experiment; SLC6A3-288 and -387
   "failed to mature and traffic to the cell surface"; both OTOP1 loop replacements "still produced no
   current"; GLP1R-376 ("completely inactive") and -377 are amber substitutions that killed or crippled the
   receptor. A model trained to find tolerated sites is learning from several that were not tolerated.

## Every non-confirmed site

Sorted by verdict. Quotes are verbatim from the paper (trimmed); the full quote, numbering basis, species,
surface and function columns and notes are in [`audit.tsv`](audit.tsv).

| Site | Tier | Validation claimed | Verdict | Corrected residue | Paper quote | What the paper actually supports |
|---|---|---|---|---|---|---|
| `ABCB1-internal-101-lit` | exact-site | surface_and_function | confirmed_shifted | 105 | “a small antigenic peptide epitope (YPYDVPDYAIEGR) ... was inserted at six different positions of the Mdr3...” | mouse Mdr3 (Abcb1a) numbering, +4 on human |
| `ABCB1-internal-320-lit` | exact-site | surface_and_function | confirmed_shifted | 324 | “the tags at positions 101 and 320 were located on the extracellular side of the membrane” | mouse Mdr3 numbering, +4 |
| `ABCB1-internal-736-lit` | exact-site | surface_and_function | confirmed_shifted | 740 | “The epitope tags inserted at positions 736, 849, and 961 were localized extracellularly, whereas tags at...” | mouse Mdr3 numbering, +4 (+2 if mdr1b) |
| `ABCB1-internal-741-lit` | same-loop | function_only | confirmed_shifted | 745 | “Insertions of epitopes at positions 235, 736, 741, 849, 887, 961, and 1024 resulted in functional proteins” | mouse Mdr3 numbering, +4 (+2 if mdr1b) |
| `ABCB1-internal-961-lit` | same-loop | surface_and_function | confirmed_shifted | 965 | “The epitope tags inserted at positions 736, 849, and 961 were localized extracellularly, whereas tags at...” | mouse Mdr3 numbering, +4 (+2 if mdr1b) |
| `ADGRL2-internal-825-lit` | tolerance | surface_only | confirmed_shifted | 824 | “Enterokinase-mediated cleavage following DDDDK/ was designed to expose the native TA region for the aGPCR:...” | tag precedes T825 of the tethered agonist; NTF (1-824) replaced, not an insertion |
| `ANO1-internal-824-lit` | exact-site | surface_only | confirmed_shifted | 850 | “In non-permeabilized cells, HA epitopes were accessible to extracellular anti-HA antibody only at...” | mouse Ano1 numbering, +26 |
| `CACHD1-internal-52-lit` | same-loop | surface_only | confirmed_shifted | 38 | “the 2HA tag followed by a TEV cleavage site was introduced into the pβA-Cachd1 between the third and...” | mouse Cachd1 numbering (signal peptide 49), -14 |
| `GPIHBP1-internal-31-lit` | exact-site | surface_only | confirmed_shifted | 34 | “A vector for S-GPIHBP1 was generated by replacing sequences for GPIHBP1 residues 32 to 48 with sequences...” | mouse residues 32-48 replaced by S-tag, +3 (low-complexity, +/-2) |
| `ITGB7-internal-29-lit` | same-loop | surface_and_function | confirmed_shifted | 48 | “PA tag into 4 locations (PAins 1: 23/24, PAins 2: 29/30, PAins 3: 427/428, PAins 4: 431/432) were...” | mouse mature-protein numbering; primers give precursor 48 |
| `ITGB7-internal-427-lit` | same-loop | surface_only | confirmed_shifted | 446 | “PA tag into 4 locations (PAins 1: 23/24, PAins 2: 29/30, PAins 3: 427/428, PAins 4: 431/432) were...” | mouse mature-protein numbering; primers give precursor 446 |
| `ITGB7-internal-431-lit` | same-loop | surface_only | confirmed_shifted | 450 | “PA tag into 4 locations (PAins 1: 23/24, PAins 2: 29/30, PAins 3: 427/428, PAins 4: 431/432) were...” | mouse mature-protein numbering; primers give precursor 450 |
| `MXRA8-internal-98-lit` | same-loop | surface_only | confirmed_shifted | 99 | “The mouse Mxra8 + moo was generated by inserting residues GEQRVGEQRLGEQRV after amino acid residue 98.” | mouse Mxra8 numbering, +1 |
| `NRG1-internal-173-lit` | same-loop | detected_only | confirmed_shifted | 178 | “HA-tag (amino-acid residues: YPYDVPDYA ...) was inserted into amino-acid residue 233-234 of Nrg1 type III” | type III isoform 233|234 = canonical H178|L179 |
| `NRG3-internal-349-lit` | tolerance | detected_only | confirmed_shifted | 333 | “carboxyl-terminal to L 362 [no deletion], P 358 [deltaTDHL], S 353 [deltaILSDP], or P 349 [deltaKTDS]) ......” | mouse proNRG3; deletion motif KTDS pins mouse P335 = human P333 |
| `NRG3-internal-353-lit` | tolerance | detected_only | confirmed_shifted | 337 | “carboxyl-terminal to L 362 [no deletion], P 358 [deltaTDHL], S 353 [deltaILSDP], or P 349 [deltaKTDS]) ......” | mouse proNRG3; motif ILSDP pins mouse S339 = human S337 |
| `NRG3-internal-358-lit` | tolerance | detected_only | confirmed_shifted | 342 | “carboxyl-terminal to L 362 [no deletion], P 358 [deltaTDHL], S 353 [deltaILSDP], or P 349 [deltaKTDS]) ......” | mouse proNRG3; motif TDHL pins mouse P344 = human P342 |
| `SDC4-internal-32-lit` | exact-site | surface_only | confirmed_shifted | 27 | “Syndecan-4-HA: The HA-tag was fused into the extracellular domain after aa 32, and was a gift from John...” | rat syndecan-4 (source paper Gopal 2010), -5 |
| `SEMA7A-internal-352-lit` | exact-site | surface_only | confirmed_shifted | 353 | “an HA tag was inserted using the Q5 mutagenesis kit at an exposed loop at amino acid 352, removing only 2...” | mouse Sema7A, two residues replaced; +1 |
| `SMO-internal-153-lit` | same-loop | detected_only | confirmed_shifted | 149 | “The sequence coding for the 3-HA epitope was inserted at position 153 of mouse Smo cDNA (exon 1)” | mouse Smo numbering, -4; no non-permeabilized staining |
| `SORCS1-internal-110-lit` | same-loop | surface_only | confirmed_shifted | 144 | “using CRISPR/Cas9 to insert a hemagglutinin (HA) epitope tag in the Sorcs1 locus after the propeptide...” | donor DNA places HA after mouse 144 = human 144, not at the furin site |
| `TLCD1-internal-74-lit` | same-loop | surface_only | confirmed_shifted | 73 | “a myc-tag was introduced after aminoacids 3, 42, 74, 107, 163, 201 and 252 of Calfacilitin” | chicken Calfacilitin numbering, -1 |
| `AQP1-internal-120-lit` | exact-site | surface_only | region_only | — | “The E1 epitope at Thr-120 was protease-sensitive in intact oocytes, indicating that loop C is exofacial.” | abstract names Thr-120 only; full text behind a bot check |
| `ATP1A3-internal-109-lit` | same-loop | surface_and_function | region_only | — | “double tags at positions 119 and 318 were located extracellularly” | paper tags rat alpha1 only; alpha3 position is a paralog alignment |
| `BACE1-internal-45-lit` | same-loop | surface_only | region_only | — | “A single HA tag was inserted after the propeptide cleavage site and a stop codon was introduced after...” | "after the propeptide cleavage site"; residue from UniProt annotation |
| `BST2-internal-146-lit` | exact-site | detected_only | region_only | — | “CD317 with a Flag tag at residue 146 (CD317-Flag) was generated by PCR amplification and inserted into the...” | "at residue 146"; construct from another paper, surface not shown |
| `CD1D-internal-301-lit` | same-loop | surface_only | region_only | — | “inserting the SortA peptide recognition motif LPETG together with a glycine-/serine-based linker between...” | "between the extracellular and transmembrane domain" of mouse Cd1d1 |
| `CFTR-internal-897-lit` | exact-site | surface_only | region_only | — | “The TMD2-3HA construct was generated by Gibson assembly from pBI-CMV2-CFTR and a gene block consisting of...” | "at position S898" in an isolated TMD2 fragment; no full-length construct |
| `CLCA4-internal-594-lit` | same-loop | function_only | region_only | — | “A Flag tag was inserted into CLCA4 sequence at nucleotide position1782 after the start codon.” | nucleotide 1782 only (codon 594 or 595); western only |
| `GRIA1-internal-396-lit` | exact-site | surface_only | region_only | — | “a SNAP-tag was either inserted at the N-terminus (SNAP::GluA1; no FRET expected) or within the ATD-LBD...” | "at position 396 aa" in a Tn5 transposon construct; junction not given |
| `KCNJ1-internal-113-lit` | same-loop | surface_and_function | region_only | — | “The rat ROMK-K80M construct ... contained a single HA tag (GGYPYDVPDYAGG) inserted at Y113 in the...” | rat ROMK "inserted at Y113"; before/after unresolved |
| `MMP14-internal-111-lit` | same-loop | detected_only | region_only | — | “WT or catalytically inactive human MT1-MMP constructs with an N-terminal Flag tag between the propeptide...” | "between the propeptide and catalytic domain", citing another paper |
| `SLC46A2-internal-57-lit` | same-loop | surface_only | region_only | — | “the codon optimized CDS of SLC46A2 (pTwist-SLC46A2) was linearized with AfeI and an oligonucleotide...” | AfeI site in an undisclosed codon-optimized CDS; residue undetermined |
| `SLC6A3-internal-199-lit` | exact-site | detected_only | region_only | — | “pHluorin was inserted at multiple sites, T173/E174, P194/G195, G199/S201, G209/T210, R219/G220 outside the...” | junction written G199/S201; "data not shown" screen only |
| `TMEM123-internal-34-lit` | exact-site | surface_only | region_only | — | “we added the Alfa tag to the N-terminal region following the signal peptide (Alanine 34)” | "following the signal peptide (Alanine 34)"; before/after not stated |
| `ABCC2-internal-163-lit` | exact-site | surface_and_function | not_supported | — | “Insertion of epitopes at positions 4, 163, 271, 574, and 653 produced functional proteins ... Epitopes...” | paper is MRP1 (ABCC1); ABCC2 never tagged |
| `ABCC2-internal-574-lit` | same-loop | surface_and_function | not_supported | — | “Insertion of epitopes at positions 4, 163, 271, 574, and 653 produced functional proteins ... Epitopes...” | paper is MRP1 (ABCC1); ABCC2 never tagged |
| `GABRP-internal-27-lit` | same-loop | surface_and_function | not_supported | — | “We excluded GABA A receptor subunit variants derived from 3 GABRs (GABRE, GABRP, GABRR2), whose...” | paper explicitly excludes GABRP; no pi-subunit construct |
| `KCNK5-internal-214-lit` | exact-site | surface_only | not_supported | — | “was inserted in the extracellular P2-M4 loop of human TASK-Δi20 at amino acid position 214” | construct is TASK-1 (KCNK3) delta-i20; KCNK5 never tagged |
| `TLCD1-internal-201-lit` | same-loop | surface_only | not_supported | — | “cells stained with Calfacilitin-myc6 without detergent seen by phase contrast (c) and fluorescence (d) -...” | chicken construct; paper shows this loop is intracellular (no staining without detergent) |
| `TMEM225-internal-125-lit` | exact-site | surface_and_function | not_supported | — | “A DNA fragment encoding an intron-less human PMP22 ORF ( NM_153321 ) with an internal myc epitope tag...” | paper tags PMP22; duplicate of PMP22-internal-125 assigned to the wrong gene |
| `TMPRSS9-internal-270-lit` | same-loop | detected_only | not_supported | — | “A Bss HII site was created in the coding sequence (positions 809 to 814) which allowed us to introduce a...” | construct is secreted polyserase-3, a different gene; no surface staining |
| `ADGRG6-internal-550-lit` | tolerance | surface_and_function | wrong_species_numbering | — | “Zebrafish ADGRG6 was cloned in the pcDNA3.1(+) vector ... All constructs with an amber codon (TAG)...” | zebrafish Y551, no human residue-level equivalent |
| `GABRA2-internal-212-lit` | tolerance | surface_and_function | wrong_species_numbering | — | “Amber mutants S181TAG and S201TAG showed the most efficient labeling demonstrated by a continuous...” | S181 is Ser in no known numbering of human alpha2; +32 offset is the agent's |
| `GABRA2-internal-232-lit` | tolerance | surface_only | wrong_species_numbering | — | “Amber mutants S181TAG and S201TAG showed the most efficient labeling demonstrated by a continuous...” | S201 as above |
| `TRIL-internal-499-lit` | same-loop | surface_and_function | wrong_species_numbering | — | “Sequence encoding an HA epitope tag was inserted into the ECD of Xenopus tril.S between amino acids 499...” | Xenopus tril.S 499|500 maps to human ~473-475, unpinnable |
| `ZP4-internal-468-lit` | same-loop | detected_only | wrong_species_numbering | — | “and a V5 tag to ZP4 (469-482 aa)” | porcine ZP4, secreted; V5 at construct coordinates 469-482 |
| `CDH13-internal-138-lit` | exact-site | surface_only | terminal_misfiled | — | “The mature form of CDH13 cDNA lacking the prodomain (amino acids 1 to 138) was subcloned into the pCAG-Ble...” | N-terminal tag on a construct deleted for 1-138 |
| `CFC1-internal-191-lit` | exact-site | surface_only | terminal_misfiled | — | “We also generated the V5-His tagged form of CR-CFC L191stop (CR-CFC L191-V5)” | C-terminal V5-His on a L191stop truncation; tag not at the surface |
| `TEK-internal-442-lit` | exact-site | surface_and_function | terminal_misfiled | — | “cDNA encoding Tie2 ectodomain (residues 1-442) followed by a FLAG-tag and residues 514-562 of...” | C-terminal FLAG on ectodomain 1-442 fused to PDGFRB TM |
| `ZP2-internal-641-lit` | same-loop | detected_only | terminal_misfiled | — | “A histidine tag (6x) was added to the N-terminus of ZP2N (38-43 aa), ZP3 (41-46 aa) and ZP4 (22-27 aa),...” | porcine secreted ZP2N; FLAG is its C-terminus |
| `ZP2-internal-643-lit` | same-loop | detected_only | terminal_misfiled | — | “A histidine tag (6x) was added to the N-terminus of ZP2N (38-43 aa), ZP3 (41-46 aa) and ZP4 (22-27 aa),...” | porcine secreted ZP2C; 6xHis is its C-terminus |

## The 15 live rows outside the training set

The API holds 15 more literature internal rows for these genes that the ranker did not train on: 3 with a
residue and 12 region rows with no residue (`insert_after_residue` null). They were read the same way and
are in `audit.tsv` with `in_training = 0`. Six of the twelve region rows can be pinned to a residue from
the paper or its primers, which the agent did not do.

| Site | Row | Verdict | Residue | What pins it |
|---|---|---|---|---|
| `CFTR-internal-None-lit` | ECL4 | confirmed_exact | 891 | validated HiBiT knock-in clone B38; donor ssODN homology arm ends at codon 891 (D891|K892); endogenous CFTR, live-cell extracellular assay and Ussing currents |
| `NRXN1-internal-None-lit` | stalk | confirmed_exact | 371 | Fig 1B of the knock-in paper shows 2xHA between mouse ...GGLA|NPTRVGG; human A371|N372, offset 0 |
| `TRPV1-internal-None-lit` | turret | confirmed_exact | 606 | Gibson primers translate to rat S606|cpHalo|L607; rat S606 = human S606 |
| `NRG3-internal-None-lit` | upstream of EGF | confirmed_shifted | 275 | "between amino acids 277 (H) and 278 (T)" is mouse numbering; human H275|T276 |
| `SLC12A5-internal-None-lit` | EL2 | confirmed_shifted | 248 | "between Ala225 and Met226" of mouse KCC2b; human canonical A248|M249 (+23 isoform offset, as for SLC12A5-406) |
| `SLC4A2-internal-None-lit` | ECL3 | confirmed_shifted | 862 | mouse AE2 E858|A859 aligns to human E862|V863 (+4; low-complexity loop, ±2) |
| `ITGB7-internal-23-lit` | 23 | confirmed_shifted | 42 | PAins1 primer translates to after mouse precursor Q42 (mature 23); same +19 mature-numbering shift as the three training ITGB7 rows |
| `ABCC1-internal-4-lit` | 4 | confirmed_exact | 4 | abstract names position 4 as extracellular and functional; ±1, abstract only |
| `CNR1-internal-25-lit` | 25 | not_supported | — | the paper's construct is an N-terminal GLuc fusion; "between residues 25 and 26" describes another group's work (its ref 27) |
| `ABCB1`, `ATP1A3`, `CLCA4`, `NRG1`, `PMP22`, `SLC6A3` region rows | — | region_only | — | the paper names only the loop or domain (PMP22's source, Tobler 1999, says "third loop" with no residue; ATP1A3 is again the rat alpha1 paralog) |

Validation claims are inflated here too: NRG3 claims surface and function with neither measured on the
tagged construct, CLCA4's function was assayed untagged, PMP22's function is cited rather than measured,
and SLC12A5 has no surface assay.

## Sites dropped as not tolerated

Five confirmed positions are excluded from the positive set because the paper reports that the tagged
construct failed: the position is right, the tag was not tolerated there. All five are substitutions, not
insertions. They are written to [`failed_sites.tsv`](failed_sites.tsv) as candidate negatives.

| Site | What the paper reports |
|---|---|
| `SLC6A3-internal-288-lit`, `SLC6A3-internal-387-lit` | tetracysteine EL3 and EL4 constructs "failed to mature and traffic to the cell surface"; uptake severely compromised |
| `NFASC-internal-533-lit`, `NFASC-internal-570-lit` | amber-codon mutants excluded by the authors for poor expression and weak or absent labeling |
| `TRPV1-internal-464-lit` | the ncAA mutant "was not efficiently labelled" |

The four other surface-negative confirmed sites stay in: PTH1R-93 has wild-type cAMP and ERK signaling
(folding shown by function), and SHISA9-29, SMO-153 and LRRC8A-76 were simply never assayed at the surface.

## Gold standard

[`gold_standard.tsv`](gold_standard.tsv), from [`gold_standard.py`](gold_standard.py): every confirmed row
from `audit.tsv` plus the added sites below, minus the five failed sites, deduplicated on (accession,
residue). Residue is `corrected_residue` where the audit set one, in human canonical `insert_after`
numbering.

| | Sites | Proteins |
|---|---|---|
| Gold standard | 156 | 113 |
| with surface shown | 120 | 93 |
| surface negative but kept (see above) | 4 | 4 |
| with function shown | 75 | — |
| of the 127 training proteins | — | 108 |

Three sites merge two readings: SDC4 27 (two papers), GABRA4 39 (the training row and the added-sites
reread of the same paper) and GABRG2 43 (HA and FLAG constructs).

## Added sites: real insertions behind wrong-gene rows

Four not_supported rows cited a paper that tags a different gene. Rereading those papers for the gene they
actually tagged gives six new sites, in [`added_sites.tsv`](added_sites.tsv); they are not live API rows.

| Gene | Residue | Flank | Tag | Source | Evidence |
|---|---|---|---|---|---|
| KCNK3 | 214 | L\|Q | HA | the KCNK5-214 row's paper ([PMID 26794006](https://pubmed.ncbi.nlm.nih.gov/26794006/)) | anti-HA on intact oocytes; numbering taken as full-length because 214 sits in the P2-M4 loop (205-222) and the Δi20 deletion is defined only in a paywalled paper |
| GABRA5 | 35 | T\|S | HA | the GABRP-27 row's paper ([PMID 27622563](https://pubmed.ncbi.nlm.nih.gov/27622563/)) | non-permeabilized flow cytometry; "between amino acids 4 and 5 of the mature" subunit, plus the UniProt signal peptide |
| GABRA6 | 23 | V\|E | HA | same | same |
| GABRG1 | 39 | A\|D | HA | same | same |
| GABRG2 | 43 | D\|D | HA and FLAG | same | same; the γ2L insert is downstream (377-384) so 43 holds for both isoforms |
| GABRG3 | 21 | K\|V | HA | same | same |

GABRA4 39 from the same paper was already in training and confirms the signal-peptide arithmetic. Function
is not attributed to the tagged subunits because the paper never states the tagged subunit was in the
recorded cells. Nothing is recoverable for CNR1-25 (its source, McDonald 2007, is subscription-only and its
abstract describes GFP and pHluorin fusions, not an HA site) or TMPRSS9-270 (polyserase-3 is PRSS53, secreted).

## Region-level labels

[`region_labels.tsv`](region_labels.tsv), from [`region_labels.py`](region_labels.py): for every
region_only or wrong_species_numbering row, the window of human canonical `insert_after` values the paper's
wording supports, so the ranker can use these sites as loop-level positives. 22 sites, 19 with a window on
16 proteins; 8 training proteins are covered only this way (ATP1A3, AQP1, BST2, CD1D, CLCA4, GRIA1, KCNJ1,
SLC46A2 among them, plus TRIL and ADGRG6 by alignment). Window widths run from 1 (BACE1, MMP14: the junction
after an annotated propeptide) through 2 ("at residue N") to 64 (a whole loop). GABRA2's two sites and ZP4
have no defensible window and are written with an empty one.

| Basis | Sites | Window |
|---|---|---|
| paper says "at residue N" | 10 | N-1..N |
| UniProt loop, domain linker or propeptide end named by the paper | 7 | the feature's interior |
| ortholog position mapped into a poorly conserved stretch | 2 | mapped position ±3 |
| none (unstated numbering frame; secreted construct) | 3 | empty |

## Recovery pass on the unpinned sites

Eight region or wrong-species sites were reread from the source the first paper points to. Two were pinned:
MMP14 111 (its ref 33, Weaver 2014, "A FLAG tag was inserted at the end of the propeptide (between Arg111 and
Tyr112)", with surface and collagen-degradation data) and TMEM123 33 (the preprint's oligo table labels the
site "A33" and the donor translates to residues 14-33, GS-ALFA-GS, then A34; "Alanine 34" in the text is the
first residue after the insert). Six stay where they were because the source is paywalled or silent: KCNJ1
(Yoo 2003, no open copy), GABRA2 (the Addgene #49169 sequence is behind a login), CD1D and CLCA4 (no
primers in text or supplement), SLC46A2 (the codon-optimized CDS is not published), BACE1 (supplement has no
construct detail). Of the 127 training proteins, 108 are in the gold standard and 8 more have a region
window; the remaining 11 (GABRP, KCNK5, TMEM225, TMPRSS9, ABCC2, CDH13, CFC1, TEK, ZP2, ZP4, GABRA2) need
papers the corpus does not have.

## Caveats on the audit itself

| Limitation | Sites affected | What would resolve it |
|---|---|---|
| Abstract-only evidence: the publisher page returned HTTP 403 or a bot check and there is no PMC copy (pre-2000 J Biol Chem and Biochemistry topology papers, plus one 2004 J Virol). The abstracts name the positions, so these are not `paper_inaccessible`, but "at position N" leaves before/after/replacing open. | 17: ABCB1 ×5, ABCC1 ×4, ATP1A1 ×2, SLC19A1 ×2, ERVW-1 ×3, AQP1-120 | Reading the methods through institutional access |
| "At residue N" wording with no before/after; graded confirmed when the row's choice is one of the two readings | 20 flagged in notes | Construct sequences or primers where the paper gives them |
| One reader per site. Each verdict was made by one agent reading the paper, with the quote recorded so a second reader can check it; no site was read twice. | all 191 | Spot-check the 52 non-confirmed rows and a sample of the 124 confirmed ones |
| Secondary sources were not followed. Where a paper cites an earlier paper for the construct (KCNJ1 → Yoo 2003, SDC4 → Gopal 2010, BST2, MMP14, GRIA1), the earlier paper was read only when it was open access. | 5 | Reading the cited construct papers |

The API and the D1 records were not edited; this folder is the proposed correction, not the applied one.

## Provenance

| Element | Source | Notes |
|---|---|---|
| Training positives | `library/l1/tag-site-benchmarking/site_features.tsv` at `f07042b` on `claude/protein-tag-location-cv-83085b` (PR [#180](https://github.com/Deliverome-Project/deliverome-analysis/pull/180)) | 176 `label == 1` rows, 127 accessions; gitignored on `main` |
| Live site records | `https://api.deliverome.org/surfaceome/v1/tag-sites/{SYMBOL}`, fetched 2026-10-29 | 191 `literature_retrieved` internal rows for the 127 genes; all 176 training ids present |
| Human canonical sequences | EBI Proteins API, UniProt REST as fallback, 2026-10-29 | UniProt REST returned HTTP 500 for a subset of accessions that day |
| Paper resolution and open-access text | Europe PMC REST (search, `fullTextXML`), 2026-10-29 | 139 cited papers (134 as primary source); 135 PMIDs resolved there, 4 bioRxiv DOIs read on bioRxiv; 124 full texts fetched here; the rest were read through PubMed Central, bioRxiv or the publisher page where one was reachable, and 17 sites rest on an abstract |
| Ortholog sequences | EBI Proteins API (mouse, rat, chicken, zebrafish, pig, *Xenopus* entries named in `notes`) | aligned locally around each site |
| [`fetch_sites.py`](fetch_sites.py) | this folder | API rows, sequences, `data/sites.tsv` with the training flag and flank check |
| [`fetch_papers.py`](fetch_papers.py) | this folder | `data/papers.tsv`, abstracts, full-text XML |
| [`merge_audit.py`](merge_audit.py) | this folder | merges per-batch verdict TSVs into `audit.tsv`, fails if any training site lacks a verdict, prints the tallies above |
| [`audit.tsv`](audit.tsv) | this folder, committed with `git add -f` | one row per live literature internal row, 191 rows (176 `in_training = 1`), 19 columns |
| [`gold_standard.py`](gold_standard.py), [`gold_standard.tsv`](gold_standard.tsv), [`failed_sites.tsv`](failed_sites.tsv) | this folder; also under `processed://library/l1/tag-site-corpus-audit/` | 156 sites on 113 proteins (see Gold standard); 5 failed sites as candidate negatives |
| [`added_sites.tsv`](added_sites.tsv) | this folder; also in S3 | 6 real sites recovered from wrong-gene papers, 2 not_supported rows kept for the record |
| [`region_labels.py`](region_labels.py), [`region_labels.tsv`](region_labels.tsv) | this folder; also in S3; UniProt features cached in gitignored `data/uniprot_features/` | 22 region-level sites, 19 with a window |
| `data/` | gitignored | downloaded API rows, sequences, papers; regenerate with the two fetch scripts |
| Pull request | [PR 236](https://github.com/Deliverome-Project/deliverome-analysis/pull/236) | draft, against `main` |

## Decision log

| Date | Was | Now | Why |
|---|---|---|---|
| 2026-10-29 | 176 positives taken at the agent's word | 146 residue-level positives (22 at corrected residues), 13 region labels, 17 to drop | this audit |
| 2026-11-05 | the 15 live rows outside training unaudited | 6 of 12 region rows pinned to a residue, ITGB7-23 shifted to 42, CNR1-25 unsupported | Becca asked whether any were recoverable |
| 2026-11-05 | 9 surface-negative confirmed sites kept as positives | 5 failed sites dropped to `failed_sites.tsv`; 4 unassayed kept | the paper shows those 5 constructs did not traffic or label |
| 2026-11-05 | MMP14-111 region_only, TMEM123-34 region_only | MMP14 confirmed_exact 111, TMEM123 confirmed_shifted 33 | second read of the cited construct paper and the preprint's oligo table |
| 2026-11-05 | wrong-gene rows discarded | 6 real sites on KCNK3 and five GABA-A subunits added | Becca: "if it's a different gene but real let's keep that" |
