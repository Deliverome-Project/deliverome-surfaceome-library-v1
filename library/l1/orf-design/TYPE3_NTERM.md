# N-terminal tag design for the type III rung

The L1 members on the "N at residue 1, single-pass (type III)" rung of
[TAG_PATHS.md](TAG_PATHS.md) take the tag cassette directly after Met1, ahead of a
transmembrane helix (TM) that must translocate the N-terminus without a signal peptide.
This page tests whether an Igκ leader ahead of the cassette is needed for them.

The analysis was run on the 62 members the rung held on 2026-11-10. They turned out to
hold seven populations, not one. Becca chose option B: the leader only where a specific
risk was found. She then applied an independent second review of the 18 risk members and
extended its surface-form test to the five GPI members. That sent 15 members for a library
change rather than a tag design: a redraw, a reference swap or a gene-model repair. Those
changes, the Compara 116 ortholog adoption and the isoform cohort gate leave the
rung at 45 members. Becca then delegated the remaining protein-level calls ("use your best
judgment"). A second membership pass moved LMLN to its signal-peptide form, off this rung,
and swapped one SUSD3 isoform for another. All 44 members now on the rung are built: 8 with
Igκ + tag and 36 with the tag alone. Option B rests on the tag alone being tolerated by the genuine
type III members that carry it, and a five-pair pilot tests that.

## Terms

| Term | Meaning on this page |
|---|---|
| Cassette | The L1 tag cassette, GGGGS–SpyTag003–GGGGSGGGGS–P-ALFA-P–GGGGS: 51 aa, net charge +3 (K+R − D−E; one His left uncharged). Asserted against every issued encoding in `l1_cassette_encodings.tsv`. The leading GGGGS stays even at residue 1 (see "Why the full cassette stays at residue 1"). |
| Tag alone | The current rung design: Met1, cassette, native residues 2…end. |
| Igκ + tag | Murine Igκ V-J2-C leader METDTLLLWVLLLWVPGSTG, the +1 Asp every Igκ vector keeps, then the tag-alone protein (Met1, cassette, native residues 2…end). The mature protein is the tag-alone protein plus one Asp, so the two designs differ by the leader alone. This is the layout of the ordered GPCR rescue arms p0041, p0043 and p0045 (`build_gpcr_rescue_arms.py`). |
| N-tail | Residues before TM1, Met1 included. The cassette adds 51 to it. |
| Δcharge | Charge of the 15 residues after TM1 minus the 15 before it. Positive favors an extracellular N-terminus ([Hartmann, Rapoport & Lodish, 1989](https://pubmed.ncbi.nlm.nih.gov/2762295/)). |
| TM hydropathy | Mean Kyte-Doolittle score over the v1 TM1 span ([Kyte & Doolittle, 1982](https://pubmed.ncbi.nlm.nih.gov/7108955/)). |
| v1 / v2 | DeepTMHMM 1.0.24 (topology of record) and DeepTMHMM2, matched to each member by sequence hash. |

## How the rung is defined, and what it actually contains

`l1_tag_paths.classify` puts a member on "N at residue 1" when its v1 string starts
extracellular with no signal peptide, then splits multipass (two or more TM helices) from
the rest. "The rest" includes members with **zero** TM helices, so the rung is "N-terminus
extracellular, no signal peptide, fewer than two TMs", which is broader than type III.

Each member was read as its own ordered protein, and its gene's reviewed human UniProt record
was aligned to it to see whether the ordered form starts at, before or after the gene's
signal peptide (`type3_nterm.py`). That splits the 62 into seven subgroups:

| Subgroup | Members | Defining test | Who |
|---|---:|---|---|
| No TM | 6 | v1 calls no TM and no signal peptide (globular) | LMLN; BST1, CDH13, OPCML, LYPD5 isoforms; cyno ACHE |
| Gene signal peptide absent from the ordered form | 12 | The UniProt record has a signal peptide; the ordered form starts downstream of it or carries a divergent N-terminus | CACHD1, CD200R1L, LRRC37B; BTN3A3, CPD, IL9R, PDPN, SIGLEC15, ZP3 isoforms; cyno CEACAM8, TIE1, TNFRSF1B |
| Gene signal peptide present but uncalled | 3 | The ordered form starts at UniProt residue 1, so it carries the annotated signal peptide; neither v1 nor v2 calls it | BMPR1B; C11orf87 (human, mouse) |
| Dual topology | 2 | The protein is known to insert in both orientations | MRAP, MRAP2 |
| Type III, short N-tail | 2 | N-tail ≤ 20 | RELL2, TEX38 (14 each) |
| Type III, mid N-tail | 26 | N-tail 21–129 (median 54.5) | KCNE2–5, GYPC ×3, TNFRSF17 ×2, TNFRSF13C, LAX1, SMAGP, PMEPA1, PRRG1 ×2, PRRG3, SUSD3 ×3, SIGIRR, OPALIN, TEX29, SERTM2, C16orf54, NRG4, KCNN2 isoform |
| Type III, long N-tail | 11 | N-tail ≥ 130 (139–996) | NRG1 ×5, KLB, TNFRSF13B, EDA2R, IGSF5, IGSF23, cyno TMEM132B |

The two N-tail cutoffs come from the literature below. At ≤ 20 residues the cassette moves a
member across the ~20-residue boundary where Kocik et al. saw insertion behavior change. At
≥ 130 the cassette is under 40% of an N-domain the TM already translocates natively.

By candidate set the 62 are 35 library references, 16 isoforms and 11 orthologs. Isoforms and
orthologs are over-represented in the first two subgroups (14 of 18), because an alternative
start or an ortholog gene model is the usual way a signal-peptide gene loses its signal
peptide.

## What decides signal-anchor orientation

A type III protein depends on its one TM acting as a signal-anchor that translocates the
N-terminus. The published determinants, each verified against its PubMed record:

| Determinant | Finding | Source |
|---|---|---|
| Flanking charge | The more positive flank stays cytosolic. Adding four Lys to the N-terminus flipped a bacterial protein from N-out to N-in. | [von Heijne, 1989](https://pubmed.ncbi.nlm.nih.gov/2677744/) |
| Which residues count | The charge difference over the 15 residues on each side of the first signal-anchor predicts its orientation. | [Hartmann, Rapoport & Lodish, 1989](https://pubmed.ncbi.nlm.nih.gov/2762295/) |
| Charge against hydrophobicity | Whether a hydrophobic segment acts as a signal or a signal-anchor depends on the balance between N-terminal charge and segment length. | [Sakaguchi et al., 1992](https://pubmed.ncbi.nlm.nih.gov/1729684/) |
| Hydrophobic core | Longer, more hydrophobic signal-anchors favor N-terminal translocation; core, charge and N-domain effects add. | [Wahlberg & Spiess, 1997](https://pubmed.ncbi.nlm.nih.gov/9151664/) |
| Mechanism | N-terminal signals insert head-on (N-out) and then invert; more positive N-terminal charge speeds inversion and greater hydrophobicity slows it. | [Goder & Spiess, 2003](https://pubmed.ncbi.nlm.nih.gov/12853479/); review [Higy, Junne & Spiess, 2004](https://pubmed.ncbi.nlm.nih.gov/15461443/) |
| N-domain folding | Stably folded N-domains hinder or block N-translocation. Within 40–237 residues, the length of an unfolded N-domain did not change translocation efficiency. | [Denzer, Nabholz & Spiess, 1995](https://pubmed.ncbi.nlm.nih.gov/8557050/) |
| N-domain length | Once the hydrophilic N-domain exceeds about 20 residues, apparent signal inversion and the effect of hydrophobicity are largely lost, and longer N-domains sterically hinder N-translocation in a length-dependent manner. | [Kocik, Junne & Spiess, 2012](https://pubmed.ncbi.nlm.nih.gov/23084973/) |
| Two signals in one chain | Two signals more than 60 residues apart insert linearly: the second acts as a stop-transfer. | [Goder, Bieri & Spiess, 1999](https://pubmed.ncbi.nlm.nih.gov/10525533/) |

The two N-domain rows disagree on length. Denzer et al. varied a domain fused ahead of a
signal-anchor and saw no length effect; Kocik et al. varied the hydrophilic segment directly
preceding the signal core and did. The cassette is the second case, but neither study used a
human type III protein, so the size of the effect on these members is unknown.

## What the cassette changes

The cassette does not change Δcharge for any of the 62 members. All of its charges sit more
than 15 residues upstream of TM1, so the Hartmann window is untouched. What it does change:

- **N-tail length**, by 51 residues for every member. For RELL2 and TEX38 (14 → 65) this
  crosses Kocik's ~20-residue boundary; for the mid-tail members it adds 43–232% to the
  native N-domain.
- **Total N-tail charge**, by +3. This is the effect Goder & Spiess tie to faster inversion
  toward N-in, but it was measured on N-terminal signals, and with the cassette in place none
  of these members has one.
- **Any hydrophobic segment upstream of TM1**, which the cassette pushes further from the
  N-terminus. That matters where the N-region already holds a signal-peptide-like stretch,
  which is how five members flip in silico (next section).

## Re-prediction of each design

Every member was re-predicted with v1 in three forms: native, tag alone, and Igκ + tag
(`type3_arms.fasta`, 186 sequences; results in `type3_nterm_calls.tsv`). The native form
reproduces the stored string's single TM for all 56 members that have one.

| Arm | Cassette extracellular | Cassette cytosolic | Signal peptide called |
|---|---:|---:|---|
| Tag alone | 57 | **5** | none |
| Igκ + tag | **62** | 0 | 20 aa for 61; 17 aa for KLB |

The five that flip with the tag alone gain a second TM and put the cassette in the cytosol.
Four of them have a signal-peptide-like stretch in the N-region that the cassette pushes
inward, where v1 reads it as a TM:

| Member | Subgroup | Hydrophobic stretch (Kyte-Doolittle, 12-residue window) | v2 on the native form |
|---|---|---|---|
| KLB `l1-m2478` | long N-tail | residues 39–50, 2.57 | signal peptide, 52 aa |
| IL9R isoform `l1-p2451, retired` | gene SP absent | 74–85, 2.87 (the canonical signal peptide, now internal) | 2 TM, N-terminus cytoplasmic |
| CEACAM8 cyno `l1-p2549, retired` | gene SP absent | 43–54, 2.59 (after a 32-residue N-terminal extension) | signal peptide, 67 aa |
| TNFRSF1B cyno `l1-p2623, retired` | gene SP absent | 110–121, 1.87, beginning MAPAAVWAALAVGLEL, the human signal peptide's start | TM |
| CDH13 isoform `l1-p2417, retired` | no TM | the tag alone creates a TM call | globular |

With Igκ + tag, all five come back with the cassette extracellular and one TM. The
predictor does not register the length effect for any member whose TM1 is its only
hydrophobic segment: all 56 such members keep N-out with the tag alone. It was trained on
natural proteins, so a "no change" call is uninformative about Kocik's mechanism, and the
pilot below is what tests it.

v2 disagrees with v1 on 7 of the 62 natives, in two ways. It calls a signal peptide on
BMPR1B, CD200R1L, IGSF5, KLB, TNFRSF13B and cyno CEACAM8 (UniProt annotates none on
IGSF5, KLB or TNFRSF13B). It calls the IL9R isoform N-in with two TMs. v2 agrees with v1 on the other 54 that have a v2 row,
including "globular" for the six no-TM members; SUSD3 isoform `l1-m0168` has none. The IL9R isoform is the one member
where the "really type II" concern has support: two predictors and its own hydrophobic
stretch disagree about which end is out.

## Refold check

Each built construct was refolded on 2026-11-12 by session `l1-tag-refold`: AlphaFold2
model_1_ptm, 3 recycles, no relax, reusing the native MSA. Tag-alone constructs fold the
precursor and Igκ constructs the mature chain. The verdicts
(`processed://sprint4/s4e9-orf-design-fg-panel/predicted-structures/tagged/2026-11-12/verdicts_type3.tsv`)
are matched to the constructs by the folded sequence's sha256 and copied onto each row.
`refold_verdict`, `refold_reason` and `refold_disposition` hold the verdict and what was done
about it. `refold_folded_form` and `refold_tagged_sha256` say which sequence the verdict is
for: the precursor for tag-alone rows, the mature chain for Igκ rows. An Igκ row's verdict
therefore never matches its `protein_sha256`, which hashes the precursor.

All 44 built constructs are folded: 40 pass, including all 8 Igκ constructs, and 4 are
flagged. The four flags are tag-alone constructs with the same reason: AlphaFold models 10–11 residues of SpyTag003
as a confident helix. They are KCNE3 (human and mouse), KCNE4 and SMAGP.

These four flags are not topology failures, so the four keep their design:
- **The helix touches nothing.** It makes no contact with the native protein, so it does
  not pack as a second membrane segment.
- **TM1 stays one confident helix,** with pLDDT 87–94.
- **The native domain is intact,** with lDDT 0.96–1.00 against the native model.
- **DeepTMHMM v1 calls no TM or signal peptide in the cassette** of these arms.
- **Passing members show the same helix.** Five passing constructs carry SpyTag003
  helices of the same length that were not called confident.
- **A helix in the free tag means little.** Free SpyTag peptides are disordered and bind
  SpyCatcher as a strand.

The Igκ leader was rejected as a remedy, because the cassette follows Met1 in both designs.
KCNE3 is a pilot pair, which measures display directly.

The fifth flag is IGSF5. Its native TM1 folds as a 9-residue helix in the tag-alone
construct. The sequence explains the break: VGLGLAGTMLLTPTCTLTI has a Pro mid-helix and is
Gly- and Thr-rich; at Kyte-Doolittle 1.45 it is the weakest on the rung.

The flag prompted a closer look at IGSF5's N-terminus, which changed its design to Igκ + tag:
- **Human IGSF5 has no signal peptide.** SignalP 6 gives P(SP) 0.000. The proteoform audit
  confirmed that Met1 is the true start and that the loss is sequence divergence.
- **The orthologs keep one.** Mouse and cyno IGSF5 both carry a signal peptide.
- **Tag alone would be a hard translocation.** It would need the weak TM1, acting as a
  signal-anchor, to translocate a folded 266-residue Ig ectodomain. Folded N-domains block
  that ([Denzer et al., 1995](https://pubmed.ncbi.nlm.nih.gov/8557050/)).
- **The leader gives co-translational entry,** as the orthologs' own signal peptides do,
  and v1 keeps the cassette outside and TM1 as the single TM in the Igκ arm.

Two alternatives were rejected. Replacing residues 1–34 with the leader would start the
mature protein at the ortholog-shared S35, but it deletes human sequence and departs from
the shared Igκ layout. Keeping the tag alone rests on unassisted translocation of a folded
Ig domain.

The IGSF5 Igκ construct passes: the native domain matches (lDDT 0.984) and the cassette
adds no helix. Its TM1 again folds as a 9-residue helix at pLDDT about 62, but IGSF5's own
AlphaFold DB model has no helix there at all (pLDDT about 50). The broken TM is therefore
an AlphaFold gap for this sequence in both models, not an effect of the tag. That also
explains the earlier tag-alone flag. The refold check now fails a TM only when the native
model has a TM1 helix of at least 15 residues and the tagged fold loses it.

## Why the full cassette stays at residue 1

The cassette's leading GGGGS looks dispensable at a terminus, because nothing but Met1
precedes SpyTag003 there. Dropping it was tested and rejected. Without it, SpyTag003's
hydrophobic VPHIVMV sits right behind Met-Arg, a basic residue followed by a short
hydrophobic core: the shape of a signal peptide's n- and h-regions. DeepTMHMM v1 then called
a cleavable signal peptide on 12 of the 62 tag-alone arms, cleaved at residue 12 (17 for one
arm), inside SpyTag003. A cleavage there would remove the SpyTag003 epitope and leave ALFA,
so a protein built that way could show ALFA but fail SpyCatcher capture. The calls fall as
the spacer after Met1 grows:

| Residues between Met1 and SpyTag003 | Tag-alone arms called a signal peptide (of 62) |
|---|---:|
| none | 12 |
| S | 7 |
| G | 7 |
| GS | 3 |
| GG | 2 |
| GGS | 0 |
| GGGGS (full cassette) | 0 |

Among members built with the tag alone, the trimmed cassette was called for eight: NRG1 (all
five forms), PRRG1 (two) and the KCNN2 isoform. With the Igκ leader in front, the trimmed
cassette drew no call, because SpyTag003 then sits behind the leader's own cleavage site.
The full cassette keeps one design for every residue-1 member. These are predictor calls,
not SignalP 6, so they flag a risk rather than prove a cleavage. The spacer arms were run
with the same DeepTMHMM v1 runner (`spacer_test`, below).

## Published tagged overexpression of these proteins

Seven of the 62 members' proteins (or a close family member) have been overexpressed in
cells with a short tag at or just after Met1 and no leader. All seven reached the surface
N-out. No publication added an Igκ or other cleavable leader to any of them as a membrane
protein, and none tested a tag as long as the 51-aa cassette.

| Protein | Construct | Readout | Source |
|---|---|---|---|
| KCNE1, KCNE2 | HA at the N-terminus, no leader | Non-permeabilized staining in HEK293; a trafficking mutant fell to 6%, so the assay detects failure | [Hu et al., 2019](https://pubmed.ncbi.nlm.nih.gov/31679457/) |
| KCNE5 | HA in the extracellular N-terminus | Surface staining in HEK293 | [David et al., 2015](https://pubmed.ncbi.nlm.nih.gov/26242757/) |
| KLB | FLAG at the N-terminus | Surface antibody binding in COS7 | [Xu et al., 2017](https://pubmed.ncbi.nlm.nih.gov/28754744/) |
| OPALIN | HA at the N-terminus | Seen with and without permeabilization; a C-terminal HA only with | [Teng et al., 2024](https://pubmed.ncbi.nlm.nih.gov/39083419/) |
| SIGIRR | Myc or FLAG at the N-terminus | Non-permeabilized staining | [Ueno-Shuto et al., 2022](https://pubmed.ncbi.nlm.nih.gov/35887095/) |
| MRAP | 15-aa biotin-acceptor peptide right after Met1 | Native dual topology and ACTH signaling kept | [Maben et al., 2016](https://pubmed.ncbi.nlm.nih.gov/27486435/) |

These tags were 8–15 residues, so the record supports tagging this class after Met1 but
does not cover the cassette's length, which is what option B rests on.

Caveats that bear on the design:

- **TEX38.** Mouse TEX38 alone stays in the Golgi of COS-7 cells and reaches the plasma
  membrane only with ZDHHC19 ([Kaneda et al., 2025](https://pubmed.ncbi.nlm.nih.gov/40030029/)).
  A low surface signal for TEX38 here is therefore uninformative about the tag.
- **KCNE3.** An N-terminal HA altered KCNQ1–KCNE3 gating in oocytes
  ([Gage & Kobertz, 2004](https://pubmed.ncbi.nlm.nih.gov/15572349/)).
- **KCNE1, KCNE2.** Their surface level depends on co-expressed KCNQ1 (Hu et al., 2019).
- **BCMA (TNFRSF17).** The one published leader addition was on soluble BCMA. Its
  ectodomain alone stayed intracellular until an interferon-α2 signal peptide was added
  ([Nguyen et al., 2026](https://pubmed.ncbi.nlm.nih.gov/42528274/)). That shows the
  ectodomain cannot enter the secretory pathway on its own; it says nothing about the
  membrane-anchored protein.

## Risks of the Igκ leader

| Risk | Assessment |
|---|---|
| Leader poorly cleaved | SignalP 6 reads the first 70 residues. Leader + Asp + Met1 + cassette is 73, so every Igκ + tag construct here and the three ordered GPCR Igκ arms (p0041, p0043, p0045) present SignalP the identical input. That input has not been predicted. v1 places the cleavage site at 20 for 61 of 62 Igκ arms and at 17 for KLB. A related Igκ construct, VSIG10 with the cassette downstream, scored a cleavage-site probability of 0.982 (`build_sp_swap_arms.py`), but its first 70 residues differ. If the leader is not cleaved, it would itself act as an N-in signal-anchor. The cassette would still be translocated ahead of the native TM, so the expected failure is an extra N-terminal anchor, not a lost tag. |
| Native TM read as a second signal | The leader's hydrophobic core and TM1 are separated by at least 66 residues (Asp, Met1, the cassette and native residues 2–14 of the shortest N-tail), past the 60 at which Goder et al. saw linear insertion with the second signal as a stop-transfer. v1 keeps TM1 as the single TM in all 56 Igκ arms that have one. TM1 needs to be hydrophobic enough to stop transfer. The weakest here are IGSF5 (1.45 over 19 residues), the KCNN2 isoform (1.58) and the BTN3A3 isoform (1.63), all still called TM by v1 in the Igκ arm. |
| Member is not type III | A natively N-in protein would be forced N-out by the leader. With the tag alone, the same protein would have put the cassette in the cytosol, so neither design displays such a member correctly; the leader at least makes the failure uniform. The IL9R isoform is the one candidate. |
| Dual topology | MRAP forms antiparallel homodimers, and both termini reach the surface ([Sebag & Hinkle, 2007](https://pubmed.ncbi.nlm.nih.gov/18077336/)). Its residues 31–37, just before the TM, are required for dual topology and MC2R function, while most of the N-terminus can be deleted without losing it ([Sebag & Hinkle, 2009](https://pubmed.ncbi.nlm.nih.gov/18981183/)). MRAP2 is reported as parallel or antiparallel dimers ([Xu et al., 2022](https://pubmed.ncbi.nlm.nih.gov/35741395/)). A leader would force the uniform N-out population these proteins do not have. |
| Construct length | The leader adds 63 nt to each insert. |

## Recommendation by subgroup

| Subgroup | Members | Recommendation | Basis |
|---|---:|---|---|
| No TM | 6 | **Do not build on this rung.** Redraw the five GPI slots; LMLN is undecided | A GPI anchor is added in the ER lumen, so the C-terminal GPI signal only works on a protein that a signal peptide has already brought into the ER. Each of these ordered forms starts downstream of where its gene's signal peptide would be (canonical residue 42 to 249), and none has a hydrophobic N-terminal stretch (best 12-residue Kyte-Doolittle 0.48–1.72). As ordered they stay cytosolic. A leader would make them GPI-anchored surface proteins the cell may never make, and each gene's canonical form is already in L1 on the N-after-signal-peptide rung. LMLN is a library reference whose ordered form starts at canonical residue 233; UniProt annotates no signal peptide for it. |
| Gene signal peptide absent | 12 | **Do not build on this rung.** Redraw 6 isoform slots, swap 3 L1 references to their signal-peptide forms, repair or drop cyno CEACAM8, and repair the cyno TNFRSF1B and TIE1 CDS | Each ordered form lacks a signal peptide its gene has. A leader would display an ectodomain the cell never makes, or duplicate a reference already in the library. With the tag alone, three flip in silico. Detail per member is in the second review. |
| Gene signal peptide present, uncalled | 3 | Hold; run SignalP 6 on the 3 natives | If the signal peptide is real, a residue-1 cassette leaves with it, so these belong on the N-after-signal-peptide rung. If it is not, they follow the mid- or long-tail row. BMPR1B's annotated 1–13 is not hydrophobic (best 8-residue Kyte-Doolittle 0.45), while v2 calls 15 aa; C11orf87's 1–23 reaches 2.04 and neither predictor calls it. |
| Dual topology | 2 | **Tag alone** | Keeps residues 31–37 and the native mixed orientation; only the N-out half of the population displays the cassette. |
| Short N-tail | 2 | **Igκ + tag** | The cassette crosses Kocik's ~20-residue boundary (14 → 65); this is the strongest literature case for risk. |
| Mid N-tail | 26 | **Igκ + tag** for TEX29 and TNFRSF17 (human, cyno); **tag alone** for the other 23 | Option B first gave the leader to the five with Δcharge ≤ 1. The second review moved OPALIN and SIGIRR to the tag alone, because their N-out topology is established, and kept the leader for the two with Cys-rich N-domains and no experimental topology or tagged construct. The ≤ 1 cut is the analyzing agent's choice: the Hartmann rule uses only Δcharge's sign, and no published threshold exists. The closest calls among the 21 are the KCNN2 isoform (22-residue tail, the weakest TM at 1.58) and KCNE4 (TM 1.94). |
| Long N-tail | 11 | **Tag alone**, except KLB → Igκ + tag | The TM already translocates 139–996 native residues, and the cassette adds under 40% of that. KLB flips with the tag alone in silico and has a predictor disagreement over its N-terminal signal. |

Three ways to apply this:

| Option | Igκ + tag | Tag alone | Held or moved | Trade-off |
|---|---:|---:|---:|---|
| A. Leader wherever the native signal-anchor is untested (recommended) | 39 (gene SP absent 10, short 2, mid 26, KLB) | 12 (long tail 10, MRAP, MRAP2) | 11 (no TM 6, SP uncalled 3, 2 CDS repairs) | Fewest silent failures. One constant element across 39 constructs, so a cleavage problem would be shared, which the pilot and one SignalP run check first. |
| B. Leader only where a risk is identified | 18 (gene SP absent 10, short 2, KLB, weakest-Δcharge mid 5) | 33 | 11 | Keeps native biogenesis for most type III proteins. Rests the mid-tail members on the length effect being small. |
| C. Tag alone everywhere, pilot first | 0 | 51 | 11 | Changes nothing until data exist, but the five in-silico flips are built knowing they likely fail. |

Becca chose **option B** on 2026-11-11, then applied the second review below. That changes the rows above for the signal-peptide-absent members and for OPALIN and SIGIRR. The final designs are in "Built constructs".

## Built constructs

[build_l1_type3_arms.py](build_l1_type3_arms.py) assigns each member now on the rung a design
from its subgroup and builds the insert (`constructs_l1_type3.tsv`, `constructs_l1_type3.fasta`;
protein view `constructs_l1_type3_proteins.tsv`). Decisions are keyed by `design_id`, which
survives member renumbering. Every row records who decided its design and when (`decided_on`),
who decided the cassette, the reviewing and implementing sessions, the reason and, for
calls made under delegation, the rejected alternative. A held row names what
it waits on.

| Subgroup | Igκ + tag | Tag alone | Held |
|---|---:|---:|---:|
| Gene signal peptide annotated, not predicted | | 3 (BMPR1B; C11orf87 human, mouse) | |
| Dual topology | | 3 (MRAP, MRAP2, cyno MRAP) | |
| Short N-tail | 2 (RELL2, TEX38) | | |
| Mid N-tail | 3 (TEX29 human, cyno; TNFRSF17) | 22 | |
| Long N-tail | 3 (KLB human, cyno; IGSF5) | 8 | |
| **Total (44)** | **8** | **36** | **0** |

**Calls made under delegation** (2026-11-12, recorded in `type3_decisions.tsv` with the
evidence and the rejected alternative):

| Member | Call | Evidence | Rejected alternative |
|---|---|---|---|
| BMPR1B | Tag alone | SignalP 6 on the ordered canonical O00238: no signal peptide (P(SP) 0.000), so UniProt's 1–13 annotation is not supported. Mid N-tail 126, Δcharge +3 | N-after-signal-peptide |
| C11orf87 human | Tag alone | SignalP 6 on canonical Q6NUJ2: no signal peptide (P(SP) 0.13). N-tail 61, Δcharge +8 | N-after-signal-peptide |
| C11orf87 mouse | Tag alone | No SignalP record. Its first 24 residues, the annotated signal-peptide region, are identical to human | N-after-signal-peptide; a mouse SignalP run would firm it up |
| SIGIRR mouse | Tag alone | Human SIGIRR's type III topology is established, and Becca moved it to tag alone; the weak-Δcharge flag does not mark cassette risk | Igκ + tag (the charge rule's call) |
| IGSF5 | Igκ + tag | Human IGSF5 lacks the signal peptide its mouse and cyno orthologs keep (audit-confirmed). Tag alone would need its weak, Pro-containing TM1 to translocate a folded 266-residue Ig domain | Tag alone; replacing residues 1–34 with the leader |
| LMLN | Held: swap the L1 reference to canonical Q96KR4 (N-after-signal-peptide), done on the manifest side; Q96KR4-3 (`l2-p0894`) is the fallback | The ordered form starts at canonical residue 233, past the signal peptide SignalP 6 calls on Q96KR4 (P 0.83, cleavage after 43); as ordered it never enters the ER | Igκ + tag on the truncated form |

The SignalP calls come from the 2026-09-28 sweep on the `signalp6-runs` Modal volume.
SignalP reads only the first 70 residues, and that sweep's input was the exact sequence
ordered.

**Who left the rung, and why.** Of the 62 analyzed, 23 are no longer on it:
- **Redrawn under the isoform cohort gate (10):** the six signal-peptide-less isoforms the
  second review named (SIGLEC15, PDPN, BTN3A3, IL9R, ZP3, CPD) and the four GPI isoforms
  (BST1, CDH13, OPCML, LYPD5). Each failed the gate's signal-peptide check.
- **Moved to N-after-signal-peptide (3):** CD200R1L, LRRC37B and CACHD1, whose L1
  references were swapped to their signal-peptide forms.
- **Retired by the Compara 116 ortholog adoption (9):** cyno ACHE, CEACAM8, SUSD3, TIE1,
  TMEM132B, TNFRSF17 and TNFRSF1B, and mouse GYPC and NRG1. This retires the cyno TNFRSF1B
  and TIE1 CDS repairs too.
- **Dropped by the cohort gate (1):** the KCNN2 isoform.

**Who joined.** Six new Compara 116 orthologs: cyno MRAP, TEX29, TNFRSF13B and KLB, and
mouse SIGIRR and KCNE3. They were re-predicted with DeepTMHMM v1 in all three designs.
The second pass then moved LMLN off the rung (reference swapped to canonical Q96KR4, with
its signal peptide) and replaced the SUSD3 isoform with Q96L08-2, built with the tag alone.

**Insert layout.** Each insert is the standard L1 Esp3I wrapper: Kozak, then [Igκ leader +
GAC], ATG, the 51-aa cassette, the member's domesticated ORF from codon 2 with its
native stop, then the barcode amplicon.

**Build checks.** The build asserts each encoded protein exactly. Every insert carries only
the two wrapper Esp3I sites, alone and in the p0039 backbone, and no domestication edit was
needed.

**Cassette codons.** The cassette uses p003's codons. Linker recoding exists to stop
C-terminal members sharing the stretch next to the barcode, and an N-terminal cassette sits
upstream of the member's own ORF.

**Still open.** Barcodes are placeholders until reserved. The encoded proteins are exactly
the tag-alone and Igκ + tag arms re-predicted above. Neither SignalP 6 nor NetGPI has been
run on them.

## Pilot

Option B assumes the 51-aa cassette at residue 1 does not reduce surface display of a
genuine type III protein with decent charge balance. A paired design tests that assumption
directly: the same protein built with the tag alone and with Igκ + tag.

| Member | Why | Arms |
|---|---|---|
| RELL2 `l1-m0482` | Short N-tail (14), Δcharge 9 | tag alone, Igκ + tag |
| TNFRSF17 `l1-m0104` | Mid N-tail (53), Δcharge −1; well-characterized type III, with ectodomain antibodies that read orientation independently of the tag | tag alone, Igκ + tag |
| KCNE3 `l1-m0009` | Mid N-tail (57), Δcharge 4; built with the tag alone under B | tag alone, Igκ + tag |
| TNFRSF13C `l1-m0103` | Mid N-tail (80), Δcharge 7; built with the tag alone under B | tag alone, Igκ + tag |
| KLB `l1-m2478` | Long N-tail with an internal signal-like stretch; predicted to flip with the tag alone | tag alone, Igκ + tag |

RELL2 replaces TEX38 as the short-tail member, because TEX38 needs ZDHHC19 to reach the
plasma membrane. In each pair, one arm is the library construct and the other is built for
the pilot.

**Readout.** Live-cell (non-permeabilized) anti-ALFA or SpyCatcher staining by flow,
normalized to permeabilized total. A failed orientation then reads as total signal without
surface signal, not as non-expression. For TNFRSF17, add an ectodomain antibody to check that
the native domain is outside as well, and a γ-secretase-inhibitor condition, because
γ-secretase sheds BCMA together with the tag (Laurent et al., 2015).

**Interpretation.** If tag alone matches Igκ + tag for KCNE3 and TNFRSF13C, option B stands
for the 23 tag-alone mid-tail members. If it falls short, those members move to Igκ + tag.

## Second review of the 18 risk members

An independent review of the 18 option-B leader members was carried out on 2026-11-11 by a
separate Claude Code session (`type3-risk-review`, read-only). Its full table, with UniProt and
PubMed evidence, is `type3-risk-review.md` in the L1 session notes. It agrees with Igκ + tag
for 6 and recommends changing the other 12. Becca applied all of it on 2026-11-11, and the
built designs above include it.

| Recommendation | Members | Reason given |
|---|---|---|
| Keep Igκ + tag | RELL2, TEX38, KLB, TEX29, TNFRSF17 (human, cyno) | Short or Cys-rich N-domains with no published tagged construct; KLB flips in silico |
| Tag alone | OPALIN, SIGIRR | Native N-out is established (OPALIN's glycosylated N-terminus, [Yoshikawa et al., 2008](https://pubmed.ncbi.nlm.nih.gov/18490449/); SIGIRR already translocates a 118-residue Ig domain), so a weak Δcharge does not flag cassette risk |
| Redraw the slot | SIGLEC15, PDPN, BTN3A3, IL9R, ZP3, CPD isoforms | Each is a truncation or divergent first exon of a type I protein whose canonical form is already in L1; a leader would display an ectodomain fragment no cell makes |
| Swap the L1 reference to its signal-peptide form (already in L2) | CD200R1L (`l2-p0276`), LRRC37B (`l2-p0908`), CACHD1 (`l2-p0248`) | The ordered references start downstream of their gene's signal peptide; the SP-bearing forms would go on the N-after-signal-peptide rung |
| Repair the gene model or drop | cyno CEACAM8 | GPI-anchored; v1's 11-residue "TM" is the GPI signal and the N-terminus is a gene-model error |

The review did not cover the five GPI members, but its surface-form test applies to them
too. Becca extended it to them on 2026-11-11, so they are held for redraw. The review also recommends a γ-secretase
inhibitor condition for TNFRSF17 in the pilot, because BCMA is shed directly by γ-secretase
([Laurent et al., 2015](https://pubmed.ncbi.nlm.nih.gov/26065893/)).

## Open questions

| Question | Why it matters | What would resolve it |
|---|---|---|
| Is the Igκ leader cleaved in front of the cassette? | Shared by all 23 Igκ + tag constructs and the ordered GPCR Igκ arms | One SignalP 6 prediction on the shared 70-residue window; surface data from p0041/p0043/p0045 |
| Should the no-TM members be on this rung at all? | `classify` admits zero-TM proteins; a residue-1 tag on them is cytosolic | A rung rule that requires exactly one TM, and a panel decision for the six |
| Are the five no-TM GPI forms real surface isoforms? | Igκ + tag makes them GPI-anchored whether or not the native transcript ever reaches the ER | Evidence that the SP-less transcripts are expressed; otherwise redraw the slots. The C-terminal GPI signal is checked by sequence only, not by re-running NetGPI on the ordered forms |
| Does the 51-aa cassette at residue 1 reduce display of a genuine type III protein? | Option B builds 33 members with the tag alone on that assumption; published tags after Met1 were 8–15 aa | The pilot |

## Next steps

1. Run SignalP 6 on the shared Igκ + cassette window.
2. Run the five-pair pilot to confirm option B for the tag-alone members.
4. Repair the cyno TNFRSF1B CDS start (Met110) and the cyno TIE1 non-ATG start, then re-place both.
5. Reserve barcodes for the 56 built inserts.
6. Carry out the library changes for the 10 held members: redraw six isoform slots, swap three L1 references to their L2 signal-peptide forms, and repair or drop cyno CEACAM8.

## Materials and provenance

| Element | Source |
|---|---|
| Rung membership | `l1_tag_paths.tsv` (62 rows with placement "N at residue 1, single-pass (type III)") |
| Library membership | The manifest, tag paths, CDS record (`2026-11-12-reference-fixes`) and member numbering of [PR #283](https://github.com/Deliverome-Project/deliverome-analysis/pull/283) at 9928f4b1, merged into this branch: Compara 116 orthologs, the isoform cohort gate and redraw, and the reference swaps |
| Ortholog panel (superseded) | The panel on `origin/main` at 514eb607, before the Compara 116 ortholog panel is merged. The 11 ortholog members and every count on this page reflect that panel; both scripts read membership from `l1_tag_paths.tsv` and the manifest, so rerunning them after the Compara 116 panel lands (and the arm re-prediction for any new ortholog sequences) updates the page. |
| Ordered proteins | `processed://sprint4/s4e9-orf-design-fg-panel/library-candidates/candidate_cds.json`, via `library_design.designed_proteins` |
| Topology | DeepTMHMM v1 (`deeptmhmm-1.0.24`, `topo_2026_*` rows, asserted per member) and v2, by sequence hash through `deliverome_analysis.topology` |
| Gene annotation | Reviewed human UniProt record per HGNC symbol, gene name checked, cached in `type3_uniprot_annotations.json` |
| Arm re-prediction | Current: `processed://sprint4/s4e9-orf-design-fg-panel/type3-nterm/type3_arms_rekey_6bf341f2.v1.3line`. It holds the 2026-11-11 run on the 62-member rung plus 18 sequences for the six orthologs that joined and 3 for the replacement SUSD3 isoform; records are matched by sequence, so renumbering does not affect them. DeepTMHMM 1.0.24 on Modal (accessible-surfaceome `modal/deeptmhmm_v1_app.py` at f66c1ba7), 0 failures. Built designs: `processed://sprint4/s4e9-orf-design-fg-panel/type3-nterm/type3_arms_met1.v1.3line` (186 sequences). Spacer test: `type3-nterm/type3_spacer_test.v1.3line` (310 sequences: the five spacer variants per member) and `type3-nterm/type3_arms_terminal46.v1.3line` (the arms with no spacer). The first run, with Met1 after the cassette in the Igκ arm, is kept as `type3_arms.v1.3line` |
| Igκ leader | `processed://sprint4/s4e9-orf-design-fg-panel/library-design/igk_leader.tsv`, as used by `build_gpcr_rescue_arms.py` |
| Refold verdicts | `processed://sprint4/s4e9-orf-design-fg-panel/predicted-structures/tagged/2026-11-12/verdicts_type3.tsv`, from session `l1-tag-refold` (AlphaFold2 model_1_ptm, 3 recycles, no relax, native MSA reused); matched by folded-sequence sha256 |
| Code | [type3_nterm.py](type3_nterm.py) (features, UniProt alignment, arm FASTA); [type3_arms_topology.py](type3_arms_topology.py) (arm calls, subgroups); [build_l1_type3_arms.py](build_l1_type3_arms.py) (design per member and the inserts) |
| Built inserts | `constructs_l1_type3.tsv` / `.fasta`; ORFs from `processed://sprint4/s4e9-orf-design-fg-panel/library-candidates/candidate_constructs.json` |

    uv run --locked python type3_nterm.py
    uv run --locked python type3_arms_topology.py
    uv run --locked python build_l1_type3_arms.py --check

## Decision log

| Date | Was | Now | Why | Decided by / implemented by |
|---|---|---|---|---|
| 2026-11-10 | One rung, "topogenesis, untested", all 62 tagged after Met1 | Seven subgroups with a per-subgroup recommendation; Igκ + tag proposed for 39, nothing built | The rung admits zero-TM members and members missing their gene's signal peptide, and the cassette flips five members in silico | Analysis by Claude Code session `l1-nterm-type3`; no decision taken |
| 2026-11-10 | No-TM members: decide between Igκ + tag and redraw | Igκ + tag for the five GPI genes; LMLN separately | Their C-terminal GPI signals are intact but they lack a signal peptide, so as ordered they never enter the ER; the leader makes them standard GPI constructs | Becca ("fix these 5") / Claude Code session `l1-nterm-type3` |
| 2026-11-10 | Igκ arm laid out as leader, Asp, cassette, native from Met1 | Leader, Asp, Met1, cassette, native 2…end; Igκ arms re-predicted with DeepTMHMM v1 | Matches the ordered GPCR Igκ arms, so the leader is the only difference from the tag-alone arm and the first 73 residues are shared with p0041/p0043/p0045 | Claude Code session `l1-nterm-type3`, on reading `build_gpcr_rescue_arms.py` |
| 2026-11-10 | The five GPI members unbuilt | Inserts built in `constructs_l1_type3_gpi_igk.tsv`, p003 cassette codons, barcode placeholder | Executes the decision above; p003 codons because linker recoding only matters next to the barcode, and an N-terminal cassette is upstream of the member's own ORF | Becca / Claude Code session `l1-nterm-type3` |
| 2026-11-11 | Igκ + tag proposed for 39 (option A), nothing else built | Option B: Igκ + tag for the 18 risk members, tag alone for 33; all 56 buildable members built in `constructs_l1_type3.tsv` | Becca chose B after the A/B comparison and the literature check; B keeps native biogenesis where no risk was found | Becca / Claude Code session `l1-nterm-type3` |
| 2026-11-11 | 51-aa cassette after Met1 | Dropping the leading GGGGS proposed for terminal tags; tested, then rejected: the full 51-aa cassette stays | v1 called a cleavable signal peptide ending inside SpyTag003 for 12 of 62 tag-alone arms without the GGGGS, and for 0 of 62 with it (spacer test) | Proposal and decision: Becca ("keep full"); test: Claude Code session `l1-nterm-type3` |
| 2026-11-11 | Five GPI members built by their own script | Folded into `build_l1_type3_arms.py`; the GPI-only outputs are retired | One builder for the rung, one design table | Claude Code session `l1-nterm-type3` |
| 2026-11-11 | TEX38 as the short-tail pilot member | RELL2 | TEX38 reaches the plasma membrane only with ZDHHC19 (Kaneda et al., 2025), so its surface signal would not measure the tag | Claude Code session `l1-nterm-type3` (agent's choice) |
| 2026-11-11 | 18 risk calls by one analysis | Independent review of the 18: keep 6, tag alone 2, redraw 6, swap 3 references, repair 1; not yet decided | Becca asked for a second agent to review the risk calls | Requested by Becca; review by Claude Code session `type3-risk-review`; recorded by `l1-nterm-type3` |
| 2026-11-11 | Second review recorded, not applied | Applied: OPALIN, SIGIRR to tag alone; six isoforms held for redraw; three references held for a swap to their L2 signal-peptide forms; cyno CEACAM8 held for gene-model repair | The review found that 10 of the 18 are not natural surface forms, or have a signal-peptide form already in the library, and that OPALIN and SIGIRR have established N-out topology | Becca ("Apply all"; scope "tIII only for now": nothing carried to other rungs) / Claude Code session `l1-nterm-type3` |
| 2026-11-11 | Five zero-TM GPI members built with Igκ + tag | Held for redraw (four isoform slots, cyno ACHE from the ortholog panel) | The second review's surface-form test applies to them: each lacks its gene's signal peptide, and each gene's canonical form is already in L1 | Becca ("Redraw them too") / Claude Code session `l1-nterm-type3` |
| 2026-11-12 | Rung of 62 on main's membership | Rung of 45 on PR #283's membership (9928f4b1): 23 members left, 6 Compara 116 orthologs joined; decisions keyed by design_id; rebuilt as 7 Igκ + tag, 34 tag alone, 4 held | PR #283 carried out the redraws, swaps and ortholog adoption | Membership: Becca, executed by session l1-reviewer-package; rekey: Claude Code session `l1-nterm-type3` |
| 2026-11-12 | Mouse SIGIRR (new ortholog) would get Igκ + tag by the Δcharge rule | Built with the tag alone, marked proposed | Matches Becca's 2026-11-11 decision for human SIGIRR | Proposed by Claude Code session `l1-nterm-type3`; awaiting Becca |
| 2026-11-12 | BMPR1B, C11orf87 (human, mouse) held for SignalP; LMLN undecided; mouse SIGIRR proposed | BMPR1B and both C11orf87 tag alone; mouse SIGIRR tag alone; LMLN held for a swap to its signal-peptide form (canonical Q96KR4, chosen by the manifest writer over Q96KR4-3) | SignalP 6 calls no signal peptide on BMPR1B or C11orf87 and one on LMLN's canonical, which its ordered form lacks | l1-nterm-type3 (delegated by Becca 2026-11-12, "use your best judgment") / Claude Code session `l1-nterm-type3` |
| 2026-11-12 | Rung of 45 on the first redraw pass; LMLN held | Rung of 44 on the second pass (6bf341f2), all built: 7 Igκ + tag, 37 tag alone | LMLN's L1 reference was swapped to canonical Q96KR4 (signal peptide, N-after-SP rung). The isoform draw replaced SUSD3 `l1-p2733` with Q96L08-2 (`l1-p3152`), mid N-tail 92, Δcharge +6, tag alone | Membership: manifest writer session l1-reviewer-package; rebuild: Claude Code session `l1-nterm-type3` |
| 2026-11-12 | Refold check pending | All 44 folded: 39 pass, 5 flagged (a confident SpyTag003 helix in KCNE3 human and mouse, KCNE4 and SMAGP; a broken TM1 helix in IGSF5); the four SpyTag003 flags keep their design | The flagged helix makes no contact, TM1 and the native domain are intact, and passing members show the same helix | l1-nterm-type3 (delegated by Becca 2026-11-12, "use your best judgment") / Claude Code session `l1-nterm-type3` |
| 2026-11-12 | IGSF5 tag alone (long N-tail) | IGSF5 Igκ + tag | Human IGSF5 lacks the signal peptide its orthologs keep (audit-confirmed), so tag alone needs a weak TM to translocate a folded Ig domain | l1-nterm-type3 (delegated by Becca 2026-11-12, "use your best judgment") / Claude Code session `l1-nterm-type3` |
| 2026-11-12 | IGSF5 Igκ refold pending | Pass; refold check now 40 pass, 4 flagged (the SpyTag003 helix flags, designs kept) | The native AlphaFold model has no TM1 helix either, so the earlier TM flag was a model gap; the folding session's revised rule changes only IGSF5 | Verdict: session l1-tag-refold; recorded by Claude Code session `l1-nterm-type3` |
| 2026-11-13 | Verdicts matched to rows without saying which sequence was folded; decision rows undated | `refold_folded_form` and `refold_tagged_sha256` per row; `decided_on` on every decision row; IGSF5 added to the delegation table | Integration could not join the eight Igκ verdicts on `protein_sha256` (the precursor); close-out asked for dated decisions | Claude Code session `l1-nterm-type3` |

## References

| Key | Citation |
|---|---|
| vonHeijne1989 | von Heijne G. Control of topology and mode of assembly of a polytopic membrane protein by positively charged residues. *Nature* 341:456–458 (1989). [PMID 2677744](https://pubmed.ncbi.nlm.nih.gov/2677744/) |
| Hartmann1989 | Hartmann E, Rapoport TA, Lodish HF. Predicting the orientation of eukaryotic membrane-spanning proteins. *PNAS* 86:5786–5790 (1989). [PMID 2762295](https://pubmed.ncbi.nlm.nih.gov/2762295/) |
| KyteDoolittle1982 | Kyte J, Doolittle RF. A simple method for displaying the hydropathic character of a protein. *J Mol Biol* 157:105–132 (1982). [PMID 7108955](https://pubmed.ncbi.nlm.nih.gov/7108955/) |
| Sakaguchi1992 | Sakaguchi M et al. Functions of signal and signal-anchor sequences are determined by the balance between the hydrophobic segment and the N-terminal charge. *PNAS* 89:16–19 (1992). [PMID 1729684](https://pubmed.ncbi.nlm.nih.gov/1729684/) |
| Denzer1995 | Denzer AJ, Nabholz CE, Spiess M. Transmembrane orientation of signal-anchor proteins is affected by the folding state but not the size of the N-terminal domain. *EMBO J* 14:6311–6317 (1995). [PMID 8557050](https://pubmed.ncbi.nlm.nih.gov/8557050/) |
| Wahlberg1997 | Wahlberg JM, Spiess M. Multiple determinants direct the orientation of signal-anchor proteins: the topogenic role of the hydrophobic signal domain. *J Cell Biol* 137:555–562 (1997). [PMID 9151664](https://pubmed.ncbi.nlm.nih.gov/9151664/) |
| Goder1999 | Goder V, Bieri C, Spiess M. Glycosylation can influence topogenesis of membrane proteins and reveals dynamic reorientation of nascent polypeptides within the translocon. *J Cell Biol* 147:257–266 (1999). [PMID 10525533](https://pubmed.ncbi.nlm.nih.gov/10525533/) |
| Goder2003 | Goder V, Spiess M. Molecular mechanism of signal sequence orientation in the endoplasmic reticulum. *EMBO J* 22:3645–3653 (2003). [PMID 12853479](https://pubmed.ncbi.nlm.nih.gov/12853479/) |
| Higy2004 | Higy M, Junne T, Spiess M. Topogenesis of membrane proteins at the endoplasmic reticulum. *Biochemistry* 43:12716–12722 (2004). [PMID 15461443](https://pubmed.ncbi.nlm.nih.gov/15461443/) |
| Kocik2012 | Kocik L, Junne T, Spiess M. Orientation of internal signal-anchor sequences at the Sec61 translocon. *J Mol Biol* 424:368–378 (2012). [PMID 23084973](https://pubmed.ncbi.nlm.nih.gov/23084973/) |
| Sebag2007 | Sebag JA, Hinkle PM. Melanocortin-2 receptor accessory protein MRAP forms antiparallel homodimers. *PNAS* 104:20244–20249 (2007). [PMID 18077336](https://pubmed.ncbi.nlm.nih.gov/18077336/) |
| Sebag2009 | Sebag JA, Hinkle PM. Regions of melanocortin 2 (MC2) receptor accessory protein necessary for dual topology and MC2 receptor trafficking and signaling. *J Biol Chem* 284:610–618 (2009). [PMID 18981183](https://pubmed.ncbi.nlm.nih.gov/18981183/) |
| Xu2022 | Xu J et al. Reversion of MRAP2 protein sequence generates a functional novel pharmacological modulator for MC4R signaling. *Biology* 11:874 (2022). [PMID 35741395](https://pubmed.ncbi.nlm.nih.gov/35741395/) |
| Gage2004 | Gage SD, Kobertz WR. KCNE3 truncation mutants reveal a bipartite modulation of KCNQ1 K+ channels. *J Gen Physiol* 124:759–771 (2004). [PMID 15572349](https://pubmed.ncbi.nlm.nih.gov/15572349/) |
| David2015 | David JP, Stas JI, Schmitt N, Bocksteins E. Auxiliary KCNE subunits modulate both homotetrameric Kv2.1 and heterotetrameric Kv2.1/Kv6.4 channels. *Sci Rep* 5:12813 (2015). [PMID 26242757](https://pubmed.ncbi.nlm.nih.gov/26242757/) |
| Maben2016 | Maben ZJ, Malik S, Jiang LH, Hinkle PM. Dual topology of the melanocortin-2 receptor accessory protein is stable. *Front Endocrinol* 7:96 (2016). [PMID 27486435](https://pubmed.ncbi.nlm.nih.gov/27486435/) |
| Xu2017 | Xu C et al. KLB, encoding β-Klotho, is mutated in patients with congenital hypogonadotropic hypogonadism. *EMBO Mol Med* 9:1379–1397 (2017). [PMID 28754744](https://pubmed.ncbi.nlm.nih.gov/28754744/) |
| Hu2019 | Hu B et al. A conserved arginine/lysine-based motif promotes ER export of KCNE1 and KCNE2 to regulate KCNQ1 channel activity. *Channels* 13:483–497 (2019). [PMID 31679457](https://pubmed.ncbi.nlm.nih.gov/31679457/) |
| UenoShuto2022 | Ueno-Shuto K et al. A splice switch in SIGIRR causes a defect of IL-37-dependent anti-inflammatory activity in cystic fibrosis airway epithelial cells. *Int J Mol Sci* 23:7748 (2022). [PMID 35887095](https://pubmed.ncbi.nlm.nih.gov/35887095/) |
| Teng2024 | Teng XY et al. OPALIN is an LGI1 receptor promoting oligodendrocyte differentiation. *PNAS* 121:e2403652121 (2024). [PMID 39083419](https://pubmed.ncbi.nlm.nih.gov/39083419/) |
| Kaneda2025 | Kaneda Y et al. TEX38 localizes ZDHHC19 to the plasma membrane and regulates sperm head morphogenesis in mice. *PNAS* 122:e2417943122 (2025). [PMID 40030029](https://pubmed.ncbi.nlm.nih.gov/40030029/) |
| Nguyen2026 | Nguyen LPV et al. Engineering soluble recombinant BCMA for ide-cel labeling. *Eur J Immunol* 56:e70251 (2026). [PMID 42528274](https://pubmed.ncbi.nlm.nih.gov/42528274/) |
| Yoshikawa2008 | Yoshikawa F et al. Opalin, a transmembrane sialylglycoprotein located in the central nervous system myelin paranodal loop membrane. *J Biol Chem* 283:20830–20840 (2008). [PMID 18490449](https://pubmed.ncbi.nlm.nih.gov/18490449/) |
| Laurent2015 | Laurent SA et al. γ-Secretase directly sheds the survival receptor BCMA from plasma cells. *Nat Commun* 6:7333 (2015). [PMID 26065893](https://pubmed.ncbi.nlm.nih.gov/26065893/) |
