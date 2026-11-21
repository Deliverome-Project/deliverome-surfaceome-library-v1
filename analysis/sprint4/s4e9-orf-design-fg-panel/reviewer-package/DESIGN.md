# The L1 surface-protein ORF library: design and the case for it

L1 is a pooled, barcoded library of 2,710 cell-surface protein coding sequences (CDSs), each
expressed from the same lentiviral vector with an extracellular epitope-tag cassette. It has not
been ordered. This document sets out what L1 is meant to measure, the five decisions the design
rests on, the evidence behind each, what is still assumed, and what we propose to check before
and alongside the order. The [design brief](BRIEF.md) is the technical reference: every number
here is measured there from the library's own tables, and each section below links to the part
of it that carries the detail.

## 1. What L1 measures

Each cell carries one CDS. Its protein carries a 51-residue cassette with two epitopes,
SpyTag003 and ALFA, placed where the protein should present it outside the cell:

`GGGGS – SpyTag003 – GGGGSGGGGS – P-ALFA-P – GGGGS`

Two reagents read the cassette on live cells. SpyCatcher3 binds SpyTag003 covalently; a
nanobody binds ALFA, and in the internalization assay it carries a disulfide linker so that
label still on the surface can be stripped with a reducing agent while label taken into the
cell is protected. The two epitopes answer different questions:

| Readout | What it reports |
|---|---|
| SpyCatcher signal on live cells | the protein reaches the plasma membrane and its tag faces outside |
| anti-ALFA signal on live cells | the same, with a second, independent epitope |
| anti-ALFA retained after stripping | label that entered the cell and stayed (net internalization) |
| barcode count per sorted population | which CDS each sorted cell carries |

The first deliverable of L1 is a display map: which CDSs put the tag on the outer surface in the
screening line. The internalization screen the library was designed for builds on that map,
because a CDS that does not display cannot be scored for uptake. The screen's own design (host
line, coverage, sort scheme, hit calling) is planned separately and is not settled: the current
plan models a paired stripped-versus-unstripped readout at about 400 sorted cells per construct
per harvest, and names HEK293 as host, while recent assays use K562. This document treats K562
as the first screening line. What the library design must guarantee is narrower, and is what the
rest of this document defends: that a CDS's barcode names the protein in the cell, that the tag is extracellular when
the protein is correctly inserted, and that a CDS which fails to display fails for a reason
that can be traced.

**What a negative can mean.** A CDS whose barcode is present but whose tag is not seen may not
be transcribed, may be lost at puromycin selection (selection runs through the CDS's own
transcript), may be retained in the ER, may insert in the wrong orientation, may carry its tag on
a buried face, or may have had its signal peptide clipped inside the tag. The design reduces each
of these where it can, and section 5 lists the controls that would separate them. A negative is
therefore not yet evidence that a protein does not reach the surface.

## 2. The library

| | Count |
|---|---:|
| genes | 2,369 |
| CDSs: reference protein of each gene / alternative human isoform / mouse or cynomolgus ortholog | 2,360 / 150 / 200 |
| CDSs with a construct in the order | 2,673 |
| CDSs held (no construct passes the build checks) | 37 |
| constructs ordered (an internally tagged CDS can have two or three) | 3,255 |
| CDS length, median (range) | 1,401 nt (186–5,139; the shortest, CD52, is held) |
| full insert with tag and barcode, median (range) | 1,751 nt (527–5,423) |

The vector is `EFS – insert – IRES2 – PuroR – WPRE` between the LTRs. Each insert is Kozak, the
tagged ORF and its stop, a 20-nt barcode and a readout block for amplification and capture.
Puromycin resistance is translated from the same transcript as the CDS, so a cell survives
selection only if it transcribes the ORF.

## 3. What is measured and what is assumed

Each design claim, and the strongest evidence for it today:

| Claim | Status |
|---|---|
| every CDS encodes its designed protein, from a real transcript | measured for all 2,710 |
| every fold used in a check is of the protein ordered | measured for all 3,255 constructs |
| barcodes are distinct and screened | measured |
| the cassette displays both epitopes at the surface | shown on p003 (TfR1, C-terminal); 21 constructs at the other placements are being read in K562 on 2026-11-25 |
| an after-the-signal-peptide tag leaves cleavage intact | predicted by SignalP 6 on each tagged protein |
| the Igκ leader helps residue-1 multipass proteins insert | reported for two GPCRs without a tag; three matched pairs (CNR1, CCR5, OR51E2) are being read in K562 on 2026-11-25 |
| the internal-site ranker picks tolerated sites | 71–72% within 10 residues on held-out published sites (top three) |
| structure checks predict display | tested and not supported: working, failed and decoy sites all pass at about 9 in 10 |
| placement topology is right | agrees with experiment for 99.3% of the 1,083 CDSs that have experimental topology; untested for the rest |
| the pool keeps barcode and protein together | modeled; upper bound 5.5% of integrants report another protein |

## 4. The design in five decisions

### 4.1 Which genes

**Decision.** A gene enters if the accessible-surfaceome catalog calls it a surface protein and it
passes eight further gates: reachable from outside, mostly at the surface, membrane-anchored,
any obligate partner expressed, expressed somewhere, displayed on its own rather than only in a
complex, an ORF of at most 5,000 nt, and an ectodomain of at least 30 residues. Three named tiers
add genes the gates miss but that have direct evidence of surface display: little-studied
proteins UniProt places at the surface, GPCRs displayed in a published HEK293T screen, and
proteins with a surface pool when overexpressed. That gives 2,134 core genes and 236 tier genes;
EPPIN is removed by hand, leaving 2,369.

**Evidence.** Each gate is a field of the gene's catalog record, and each gate's genes are listed
with the value tested ([brief section 2](BRIEF.md#2-genes-from-the-surfaceome-universe-to-2369)).
An independent check uses published surface proteomics: of the 394 proteins that at least 20
intact-cell surface-labeling studies detect, and that are enriched over membrane fractionation,
all but 51 are in L1 or L2, and every one of the 51 fails a stated gate or has no catalog record.

**Where it is weakest: the host cell.** Two gates assume something about the cell the screen runs
in, and neither tests it:

- The co-receptor gate asks whether an obligate partner is expressed across all 1,206 HPA cell
  lines, not in the lines L1 will be screened in, and it requires every listed partner even where
  the list names alternatives. A screening-line rule (every partner ≥ 5 TPM in at least 2 lines of
  the screening shortlist, or in K562) keeps all 48 current genes and adds 29; over only the lines
  in house or on order it adds 9, plus RHCE and RHD. The T- and B-cell receptor chains stay
  excluded either way: their partners are present in only one shortlist line each.
- No gate asks whether the host carries the chaperones some families need. Odorant receptors
  need RTP1 or RTP2 to reach the surface of heterologous cells
  ([Saito et al., 2004](https://pubmed.ncbi.nlm.nih.gov/15550249/)), bitter-taste receptors are
  helped by RTP3 or RTP4, CHRNA7 needs NACHO and GABA-A α subunits need a β subunit. K562
  expresses none of these. That affects 213 CDSs (7.9%): 186 odorant receptors, 22 bitter-taste
  receptors, 4 GABA-A α subunits and CHRNA7. A published HEK293T screen, also without these
  factors, measured 204 of them: 85 displayed and 119 did not, so about 128 are expected to fail
  in K562 for want of a factor rather than for anything about the CDS
  ([brief section 2](BRIEF.md#cdss-that-need-a-host-accessory-factor)).

**Alternatives rejected.** A lower partner-expression threshold (it adds only integrin α chains
whose β partner K562 lacks); building every UniProt surface protein regardless of the catalog
(the catalog's evidence grading is what removes ER, Golgi and secreted proteins).

### 4.2 Which protein for each gene

**Decision.** Each gene is built as one reference protein: the protein MANE Select and the APPRIS
principal agree on (2,172 genes); the one that exists when only one does (77); UniProt
canonical when neither exists (9); MANE when they disagree, unless peptide evidence clearly favors
the other form (111, none of which has strong evidence for APPRIS). Eighteen genes take a named
exception where the chosen form would not be displayed as intended: a GPI anchor or signal
peptide missing from the chosen form, or a truncated model.

**Evidence.** The library orders a coding sequence, so the reference has to be a real transcript;
UniProt canonical is a protein, and for 44 library genes it is not the translation of any Ensembl
transcript. Every CDS translates to its designed protein, and 2,699 of 2,710 equal their pinned
Ensembl transcript base for base, the other 11 differing only by recorded edits
([brief section 3](BRIEF.md#3-proteoforms-one-reference-protein-per-gene)).

**Where it is weakest.** For 206 genes the built protein is not UniProt's canonical, so any
record looked up by accession, such as a structure model or a topology call, can describe a
different protein; everything downstream is therefore matched by sequence. Thirteen genes,
mostly olfactory receptors, are disrupted on the reference genome and are either built from an
intact copy on another scaffold or not built.

### 4.3 Where the tag goes

**Decision.** The tag goes where the protein's predicted topology puts an end or a loop outside
the cell, chosen in this order:

| Placement | Rule | CDSs |
|---|---|---:|
| after the signal peptide | the protein has a cleaved signal peptide | 1,285 |
| at residue 1, multipass | N-terminus outside, no signal peptide, more than one helix; a cleavable Igκ leader is added ahead of the tag | 565 |
| at residue 1, single-pass (type III) | the same with one helix; a leader only where a risk was found | 44 |
| C-terminal | C-terminus outside | 208 |
| internal loop | neither end outside; site chosen by a ranker and checked | 608 |

**Evidence.**

- *Topology.* Placement follows DeepTMHMM on each CDS's own protein. For the 1,083 CDSs (40%)
  whose tagged segment has experimental topology (UniProt features with experimental evidence,
  and TOPDB), the placement disagrees with experiment for 8 (0.7%; 95% CI 0.4–1.5%). Only 3 would
  put the tag on the wrong side, each resting on a single item of evidence; the other 5 are tag
  sites a few residues inside a helix end. The rate is likely optimistic for the 60% without
  evidence, which are less studied ([brief section 5](BRIEF.md#topology)).
- *Signal peptides.* After-the-signal-peptide sites are kept only where SignalP 6 still predicts
  cleavage upstream of the tag on the tagged protein; where the native signal peptide is predicted
  to fail, the Igκ leader replaces it (232 CDSs).
- *The Igκ leader at residue 1.* A cleavable signal peptide improves membrane insertion and
  surface expression of GPCRs that lack one ([Guan et al., 1992](https://pubmed.ncbi.nlm.nih.gov/1331042/);
  [Andersson et al., 2003](https://pubmed.ncbi.nlm.nih.gov/12920192/)). Neither study carried a
  51-residue tag. Three matched pairs with and without the leader (CNR1, CCR5, OR51E2) are
  transduced in K562 and are read on 2026-11-25.
- *GPI proteins* are tagged after the signal peptide (186 of 191), because the anchor replaces
  the C-terminus.
- *Internal sites.* The ranker puts one of its top three picks within 10 residues of a published,
  confirmed tag site for 71–72% of held-out proteins, against 16% for a random residue; it is not
  shown to beat low AlphaFold confidence alone. Each site must also sit outside the membrane in
  the CDS's own structure model, avoid known ligand contacts, and pass a refold check. Up to three
  sites are ordered per CDS to hedge the ranker's uncertainty.
- *The refold check.* Each tagged construct is folded with ColabFold and compared with the
  untagged protein; 3,119 of 3,255 pass. Tested against published insertions, the check does not
  discriminate: sites shown to work pass 92% of the time, sites shown to fail 89%, and decoy
  sites expected to break the protein 87–88% (AUROC 0.51–0.53). In proteins of 300 residues or
  fewer it fails half of the sites known to work. So a pass is not evidence a site works, and a
  fail is not evidence it does not ([brief section 6](BRIEF.md#how-well-the-refold-check-predicts-display)).

**Where it is weakest.** Three choices rest on arguments rather than data: tagging at residue 1
(609 CDSs, 565 of them behind the Igκ leader), replacing a native signal peptide with the leader
(232), and after-the-signal-peptide sites that fall inside the predicted fold (513). The refold
check that was meant to guard these choices does not predict display. The first and
third alone cover 1,122 CDSs, 41% of the library. Each class shares its failure mode, so a wrong
assumption would remove the whole class and look like biology. The running K562 pilot tests the
leader and after-the-signal-peptide tagging on a few proteins each; type III tagging and in-fold
sites are not yet in it (section 5).

### 4.4 How each construct is built

**Decision.** ORFs keep their native codons except where a synonymous change removes an
assembly enzyme site or makes a sequence synthesizable (790 of 3,255 constructs change, 5,312
codons in all, no protein change). The cassette's codons are those shown to display in our first
construct, p003 (TfR1 with the cassette, detected at the surface with both reagents). Every
construct sits in one vector with the EFS promoter and IRES-coupled puromycin selection.

**Evidence.** The vendor scored the first draft; local rules calibrated on its report reproduce
all 25 of its refusals, and after recoding none is refused and 23 remain high-complexity. Coupled
selection removes cells that survive puromycin without transcribing the ORF, which a separate
selection promoter allowed in an earlier construct ([brief section 8](BRIEF.md#8-the-vector-promoter-and-selection-marker)).

**Where it is weakest.** Coupling means a weakly transcribed CDS is lost at selection rather than
read as dark. EF1a gave 1.7–4.0-fold more transcript than EFS in three cell types but adds 967 bp
to every construct; K562 had no EFS well in that comparison. The promoter can be changed after
the order, because each construct carries its own transfer sites.

### 4.5 Keeping barcode and protein together in the pool

**Decision.** Barcodes sit 21–24 nt downstream of each ORF's stop, are screened against splice
donors, miRNA seeds, enzyme sites and RNA-binding motifs, and differ pairwise by at least three
substitutions. The pool is weighted by construct length to offset lower titer of longer genomes.

**Evidence.** Two genomes packaged in one particle can recombine. In a pool, a construct shares
a particle with a specific other construct rarely (median 1 in 11,000), so sequence shared by a
specific pair adds little: a median 0.001% of the cells carrying one construct's barcode carry
its sibling ([brief section 1](BRIEF.md#1-what-l1-is-for)).

**Where it is weakest.** Sequence shared by every construct matters more. Only the 310 C-terminal
and snorkel constructs have recoded linkers (SpyTag003 and ALFA stay on p003's codons in all
constructs). In the 2,942 N-terminally and internally tagged constructs the whole cassette,
linkers included, is on p003's codons, so it is identical in all of them, and a template
switch inside it would join one protein's N-terminal half to another's tag and C-terminal half.
The model puts that at up to 4.1% of integrants, an upper bound that assumes every particle is
mixed; recoding the cassette's linkers per construct would leave about 3.6%, and recoding every
codon about 1.8%. A fusion that borrows a working signal peptide could make a failing protein
look displayed. The length weighting (2- to 3.5-fold per kilobase) is not measured for this
vector, and across that range the longest construct is pooled at 13 to 100 times the median
construct's amount.

## 5. Before and alongside the order

Each check below changes what is ordered or how results are read.

| Check | What it decides | Proposed design |
|---|---|---|
| Arrayed pilot in K562 (running) | whether each placement displays, and whether the Igκ leader helps | 22 constructs in the library vector, transduced and puromycin-selected, read for SpyCatcher and anti-ALFA at 4 °C and 37 °C with and without the strip on 2026-11-24/25: five tags after a native signal peptide (ERBB2, AXL, VSIG10, NECTIN4, FOLR1), six loop sites (SLC2A4 ×2, MS4A1 ×2, CCR5, AVPR2), the C-terminal TfR1 benchmark, three residue-1 GPCRs with and without the Igκ leader (CNR1, CCR5, OR51E2), two more without (ACKR3, AVPR2), and VSIG10 with the leader in place of its own signal peptide ([S8E1](https://www.notion.so/41aab86c7683800fbb3de43755efb8e5)) |
| Extending the pilot | the classes it does not cover | type III proteins; more leader pairs, including ones whose own N-terminus is short; after-the-signal-peptide sites inside the predicted fold; constructs kept despite a failed refold; odorant receptors with and without RTP1S |
| Controls in the pool | whether a negative can be traced to its cause | per placement, positives with published display (p003 and gold-standard sites) and matched negatives (tag on the cytoplasmic side, signal peptide deleted, untagged ORF); barcodes counted from plasmid, before selection and after selection |
| Cassette encoding | the fusion term in 4.5 | recode the linkers (or every codon) before synthesis, or keep them and measure fusions by long-read sequencing of ORF and barcode from transduced cells |
| Titer per kilobase | the length weighting | four plasmids spanning the length range, titered in K562 with an EFS condition |
| Host chaperones | whether the 213 factor-dependent CDSs can be read at all | express RTP1S and RTP2 (and the other factors) in the K562 host, screen these CDSs as a separate sub-pool in a factor-expressing host, or analyze them apart and keep them out of the display denominator |
| Topology | whether any placement puts the tag inside | done for the 40% with experimental evidence (0.7% discordant); review the 4 internal sites inside a helix end (ABCC4, SLC12A5, SLC12A6, XPR1) and CATSPERH before the order |
| Refold check | whether 37 held CDSs should be held | tested: the check does not separate working from failing sites, and fails half of working sites in small proteins; we propose building each CDS held only for refold at its best site, flagged, and letting the screen decide (a decision in the cover note) |

## 6. Questions for reviewers

The nine questions we most want challenged are listed with their evidence in
[brief section 12](BRIEF.md#12-decisions-we-would-like-reviewers-to-weigh-in-on): the MANE
default, type III tagging, the genes disrupted on the reference genome, the isoform panel, the
length weighting, the promoter and selection coupling, the Igκ leader, the cassette encoding,
and the co-receptor host line.

Decisions that are ours to make before the order (the cassette encoding, the host for the
factor-dependent CDSs, the co-receptor rule, the gap genes, the isoform redraw) are listed with
what each changes in the [cover note](COVER_NOTE.md).
