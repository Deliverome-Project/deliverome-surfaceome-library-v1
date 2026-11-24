# Library design: from a surfaceome record to an orderable ORF

How a protein in the surfaceome database becomes a synthesis-ready insert. The intent is
that adding a member costs a database lookup and a decision-tree walk, not a bespoke
design session — so every branch below reads a field that already exists rather than a
judgement someone has to make again.

## Terms

| term | meaning |
|---|---|
| **cassette** | the constant tag block, `GGGGS-SpyTag003-GGGGSGGGGS-P-ALFA-P-GGGGS`, 51 aa |
| **rung** | which of the four placements a member takes; the decision tree assigns exactly one |
| **signal peptide (SP)** | N-terminal targeting sequence the signal peptidase removes |
| **ω-site** | the residue a GPI transamidase cuts at, replacing everything downstream with the anchor |
| **junction** | the bond the cassette is spliced into, written `after N` = between residues N and N+1 |
| **snorkel** | a transmembrane helix plus tag appended to a cytoplasmic C-terminus, displaying the tag extracellularly but structurally separate ([Brown 2013](https://pubmed.ncbi.nlm.nih.gov/24023844/)) |

## Inputs, and where each comes from

Everything the tree reads is already published, so no branch requires a new computation:

| field | source | used for |
|---|---|---|
| `canonical_topology.per_residue_topology` | `topology_public`, DeepTMHMM v1 | terminus accessibility, loop extents |
| `canonical_topology.signal_peptide_length` | same | SP end (DeepTMHMM call) |
| `cleavage_site` | `signalp_public`, SignalP 6 | SP end (second call) |
| `filters.subcategory == "GPI_anchored"` | surfaceome record | the GPI branch |
| `disorder_hex` | `disorder_public`, 3 predictors | post-SP junction selection |
| `plddt`, `rsa`, `ss` | tag-site signal records | internal-site gating |

**DeepTMHMM v1 is the topology of record.** Our own benchmark puts v1 at 99.5% signal-peptide
recall against v2's 92.9%, so v2 is used for membrane compartment and reentrant helices, not
for the calls this tree makes.

## The decision tree

Walk in order; the first branch that matches wins.

```
1. GPI-anchored?            filters.subcategory == "GPI_anchored"
   └─ yes → N-TERMINAL, post-SP (§ below). C-terminal is forbidden outright.
            A C-terminal cassette sits past the ω-site, so the transamidase either
            removes it with the GPI signal — anchored but untagged, reading as a
            non-expressor — or the 51-aa addition breaks the signal's length and
            hydrophobicity and the protein is secreted into the medium.
            BST2 is why this is branch 1 and not a sub-case of branch 2: it is
            GPI-anchored with NO cleaved SP, using an N-terminal signal anchor, so
            it reaches the C-terminal rung on topology alone.

2. Extracellular N-terminus behind a cleaved SP?    topology[sp_end+1] == 'O'
   └─ yes → N-TERMINAL, post-SP (§ below).  ~57% of the library.

3. Extracellular N-terminus at residue 1?           topology[1] == 'O', no SP
   ├─ multipass (≥2 TM) → N-TERMINAL at residue 1, bare.
   │                      Tedman's design; 834 of his 946 constructs, validated.
   └─ single-pass (1 TM) = TYPE III → INTERNAL, in the N-terminal ectodomain.
                           NOT bare N-terminal: one TM carries the whole topogenic
                           decision and the charge flanking it sets orientation, so a
                           cassette on the translocated N-tail risks inverting the
                           protein (§ Type III).

4. Extracellular C-terminus?                        topology[L] == 'O'
   └─ yes → C-TERMINAL.  ~6% of the library.
            Guarded: refuse unless the C-terminus really is extracellular. All eight
            of the first panel's C-terminal constructs failed this and were rebuilt.

5. Neither terminus accessible → INTERNAL, ranked by relative flexibility (§ below)
   └─ no candidate clears the gates → SNORKEL, or exclude (§ limits)
```

A member that reaches no branch is not designed. That is a result, not a failure: it
usually means the protein does not belong in a surface library (§ limits).

## N-terminal placement after a signal peptide

Anchor on disorder, not on the cleavage site. Full evidence in
[POST_SP_PLACEMENT.md](POST_SP_PLACEMENT.md); the rule:

1. SP end = the **later** of the DeepTMHMM v1 and SignalP 6 calls. The failure modes are
   asymmetric — a cassette inside the SP leaves with it and the member reads as a
   non-expressor, while one placed a few residues late is still presented.
2. Search `sp_end + 1 .. sp_end + 20`, the span Tedman's constructs occupy.
3. Take the residue with the highest mean within-protein disorder percentile across
   metapredict, NetSurfP-3.0 and AlphaFold-disorder.

This reproduces Tedman's geometry — his 108 SP-bearing insertions sit a median +4 past the
cleavage site, at roughly the 90th disorder percentile of their own receptor — and it is
robust to the two SP predictors disagreeing, which they do for 20 of 65 accessions by 1–4
residues, the same scale as the offset itself.

## Internal site selection

Candidates are scored from the persisted per-residue signals and must clear every clause:
pLDDT ≥ 70, junction-window RSA ≥ 0.4, ≥ 10 Å in 3D from any UniProt hazard feature
(disulfide, glycosylation, binding or active site), and a contiguous DSSP coil/turn run of
≥ 8 residues.

**Rank survivors by relative flexibility, not by which lane nominated them.** The pipeline
currently runs an ordered-loop lane and a disorder lane and lets the loop nomination win
when both name a residue. Neither of our two experimental datasets supports that:

| dataset | n | median disorder percentile at the junction |
|---|---|---|
| curated positive controls, internal sites only | 25 | 73rd–83rd, depending on predictor |
| Tedman screen, inserted-after-N | 112 | 88th–93rd |

Both say successful internal sites are the *most flexible residue available in that
protein*. Note the absolute scores are modest — median 0.25–0.44, mostly below the 0.5
disorder threshold — so the criterion is relative flexibility within an otherwise ordered
protein, not an intrinsically disordered region. The loop lane is already reaching for the
same thing when it ranks by indel tolerance; the two lanes are one criterion split in
half, and the preference between them is unjustified.

**The 8-residue loop requirement is a property of our cargo, not of the proteins.** It
exists because the cassette is a ~30–40 aa SpyTag+ALFA tandem needing a flexible loop to
protrude from. Of the members that reach branch 5 and fail, most have ample exposed
extracellular surface — median 30 residues clearing both pLDDT and RSA — and fail only
this clause, typically at 6 or 7 residues. A shorter tag would move most of them.

## Sequence sourcing, and what to do when the sources disagree

**Ensembl supplies the sequence; UniProt holds a veto.** `design_orfs.py` resolves the
symbol, takes the Ensembl canonical transcript, and that native CDS becomes the
synthesised ORF. The translation is then cross-checked against the reviewed UniProt entry.
Ensembl is primary because we want native codons, which requires a real transcript;
UniProt gives a sequence, not a CDS.

Tedman uses the same arrangement: genes defined from UniProt annotations, "the
corresponding canonical transcript sequences... collected from the Ensembl database", then
"cross checked their amino acid sequences against the Uniprot database"
([Tedman 2026](https://doi.org/10.1101/2025.09.19.677468)). He does not state what he does
on a mismatch.

We currently hard-stop, and that is stricter than the field. Of 268 ORFs attempted across
two batches, 23 were refused on an Ensembl/UniProt disagreement, including:

| gene | Ensembl | UniProt | the actual difference |
|---|---|---|---|
| SGCD | 290 aa | 289 aa | **one leading methionine** |
| SLC13A4 | 627 aa | 626 aa | one residue, internal |
| TMPRSS6 | 802 aa | 811 aa | 9 N-terminal residues, `MLLLFHSKR` |
| GRAMD1B | 877 aa | 738 aa | 139 residues — genuinely different proteins |

Two further facts bear on the policy. For the 17 refused genes that have a **MANE Select**,
MANE is the Ensembl canonical transcript *by definition* -- Ensembl adopts MANE as its
canonical -- so we are refusing genes where the Ensembl/RefSeq consensus is the sequence we
would have built, because UniProt dissents. And the field's own practice is to proceed with
a recorded confidence rather than refuse: the ORFeome Collaboration ships clones carrying
"clone confidence levels based on CCDS and RefSeq annotations"
([Wiemann 2016](https://pubmed.ncbi.nlm.nih.gov/27454283/)), and Rolfs selected clones "to
represent the best available GenBank reference sequence"
([Rolfs 2008](https://pubmed.ncbi.nlm.nih.gov/18213396/)). Both lean on **CCDS**, the NCBI/Ensembl agreed coding-sequence set -- but CCDS is not a
live option: its most recent release is **24, from October 2022**, and MANE identifiers were
added to CCDS report pages in November 2022 rather than the reverse. It should not be
adopted as an arbiter in 2026 on the strength of what those collections did in 2008 and 2016.

The benchmark literature points to MANE and APPRIS instead.
[Pozo 2022](https://pubmed.ncbi.nlm.nih.gov/36257961/) shows reference-transcript choice
materially changes clinical variant interpretation, that MANE Select and APPRIS principal
are the best references, and that they are "particularly powerful over the 94% of coding
genes in which they agree". The same study bears on our isoform question from the other
side: vanishingly few ClinVar pathogenic variants fall on alternative protein products, and
the alternatives that do matter are identifiable in advance by TRIFID score.

**Proposed policy, not yet implemented:** proceed on the Ensembl canonical, record the
UniProt disagreement and its size as a confidence field, and refuse only where MANE and
APPRIS *disagree* -- the 6% of genes Pozo identifies -- **and** the difference overlaps the
tag site. A one-residue start-codon
difference in a cytoplasmic tail cannot affect a C-terminal construct; 139 residues of
different protein is a real scientific question. The current gate cannot tell those apart.

## Other constraints worth matching

- **Synthesis length.** Tedman included "all receptors encoded by transcripts under 3 kb".
  We have no equivalent cutoff; our longest C-terminal ORF is TENM1 at 2,732 aa, roughly
  8.2 kb of CDS, well outside his bound.
- **Isoforms are coverage, not noise.** He added "174 experimentally-validated GPCR splice
  variants" deliberately, 940 constructs from 766 canonical receptors. We hold tag sites
  for isoforms but have built no isoform constructs, and have been treating
  canonical/isoform disagreement as a problem rather than as two orderable things.

## Numbering frames

A junction is meaningless without the frame it was measured in. Tedman's positions are in
the numbering of the transcript each plasmid was built from, which is not always the
canonical proteoform: `ENST00000579344` encodes CCR7 without its first six residues, so its
junction 23 is canonical residue 29. Comparing the integers directly made that look like a
signal-peptide disagreement -- the failure mode that silently places a tag inside the signal
peptide. `tedman.project_junction` resolves the frame by aligning the two proteins on
longest identical blocks and refuses to guess when the junction falls where they differ.

## Domestication and the orderable unit

Once a rung and junction are fixed, the insert is assembled as
`Kozak | ORF with cassette | stop | pad | primer | barcode | RT site | CS1`.

- ORF sequence: Ensembl canonical CDS, **native codons**, cross-checked against the
  reviewed UniProt entry. A length mismatch is a hard stop, not a warning.
- Type IIS sites (BbsI, Esp3I, BsaI) removed with the fewest possible edits; the protein is
  gated byte-identical before and after.
- The tag cassette gets a **per-gene synonymous encoding** so the constant region is not a
  shared template-switch substrate across the library.

## Limits, stated rather than discovered later

| limit | consequence |
|---|---|
| Ensembl canonical ≠ UniProt canonical | hard stop. 13 of 145 in the first C-terminal batch; 10 resolve to a named alternative transcript, 3 match no coding transcript at all |
| no extracellular surface | not a surface protein. ENO1, RAB5C, CAPNS1, DYNC2H1, VDAC1 and others reach branch 5 and emit nothing |
| beta-barrel topology | unhandled. VDAC1 is the library's only one and emits no sites |
| snorkel is unvalidated for internalization | Brown 2013 validated surface expression only; C-terminal lumenal tags are known to redirect trafficking, which is the readout |
| type III placement is extrapolated | the multipass arm is validated, the single-pass one is not (§ Type III) |

## Type III

26 members are single-pass, N-out, with no cleaved SP — a signal anchor instead. Adding
positive charge near the N-terminus of a type III signal sequence inverts 28–95% of the
protein ([Zhang 2004](https://pubmed.ncbi.nlm.nih.gov/15355969/)), and most documented type III
proteins carry a *negative* residue flanking the anchor, where single charge swaps flip
13–23% ([Parks & Lamb 1993](https://pubmed.ncbi.nlm.nih.gov/8103052/)). Our cassette leads with
SpyTag003, which is net positive. The N-terminal segment is also translocated
post-translationally, so its folding state and length matter
([Spiess 1995](https://pubmed.ncbi.nlm.nih.gov/7641889/)); 51 structured residues is not neutral.
An inverted protein traffics normally and presents nothing, so the failure is silent. These
members therefore take an internal site in their N-terminal ectodomain, which leaves the
charge context that sets orientation untouched.

## Decisions

| date | decision |
|---|---|
| 2026-10-08 | Post-SP placement anchors on disorder in a +1..+20 window, not the cleavage site. |
| 2026-10-08 | SP end is the later of DeepTMHMM v1 and SignalP 6. |
| 2026-10-08 | GPI is branch 1 and forbids C-terminal placement; BST2 shows it cannot be a sub-case of the SP branch. |
| 2026-10-08 | C-terminal placement is guarded on `c_terminal_orientation`; the first panel's eight C-terminal constructs were rebuilt. |
| 2026-10-08 | Type III members take an internal ectodomain site rather than a bare N-terminal tag. |
| 2026-10-08 | DeepTMHMM v1 is the topology of record (Becca). |
| 2026-10-08 | Internal sites rank by relative flexibility; the ordered-loop lane does not outrank the disorder lane (Becca, confirmed against both the positive controls and the Tedman screen). |
| 2026-10-11 | GPI membership reads the curated database attribute; UniProt KW-0336 and a cited supplement union in, and `gpi_exclusions.tsv` subtracts only where sources disagree with each other. |
| 2026-10-11 | The anchor is located per form, not per gene. NetGPI is a veto on confident negatives; source precedence still ranks. Six references move onto the anchored form. |
| 2026-10-11 | A gene with no anchored form is reviewed, never dropped — `CRIPTO`, `CFC1` and `GML` are predictor blind spots, not bad attributes. |
| 2026-10-11 | The ortholog surface guard extends from transmembrane helices to GPI anchors; six picks move onto a form that keeps the anchor. |
| 2026-10-11 | `EPPIN` leaves the library: its surface presence belongs to its complex partners, not to an anchor of its own (Becca). |
| 2026-10-11 | Downstream steps filter on `library_members.json` rather than inheriting membership from a snapshot. |
| 2026-10-11 | The ortholog panel is stratified by the human gene's class and by correspondence to the human reference. |

## Candidate isoform funnel

Alternate forms, pooled from UniProt 2026_03 and Ensembl 115 and deduplicated by protein
sequence, for the 2,369 human library members:

| stage | forms | lost |
|---|---|---|
| distinct alternate forms | **13,687** | — |
| fail the evidence/TRIFID bar | 2,110 | **−11,577** |
| lose every TM helix (already-called, s1e3 rule) | 2,037 | −73 |
| lose every TM helix (new DeepTMHMM v1 run) | **1,928** | −109 |

The bar is **observed evidence of expression** — a diagnostic peptide, or an OpenProt
translation score — **or TRIFID above 0.5**. The 11,577 that fail it are overwhelmingly
APPRIS MINOR with nothing behind them. **1,928 candidate isoforms survive**, of which 710 are
in the current cohort; forms that lose only *some* TM helices are kept, as a designer
judgement rather than a filter.

### Two evidence tiers

The funnel's 1,928 survivors are emitted at two tiers so a smaller order is a filter rather
than a rebuild:

| tier | forms | genes | already in cohort |
|---|---|---|---|
| **1 — full** | **1,916** | 957 | 633 |
| **2 — observed** | **816** | 560 | 437 |

| evidence behind each form | full | observed |
|---|---|---|
| protein, deep diagnostic peptides | 61 | 61 |
| protein, weak | 207 | 207 |
| translation only (OpenProt TE > 0) | 556 | 556 |
| **prediction only (TRIFID > 0.5)** | **1,104** | — |

Tier 2 is a **57% reduction**, and the line is deliberately *observed versus predicted*
rather than a threshold. The tighter cuts the data offers — TE ≥ 3 gives 280 forms, TE ≥ 5
gives 199 — are arbitrary points on a count of studies, and a study count measures how often
a protein has been looked at as much as whether the isoform is real. Observed-versus-predicted
is the boundary that means something, and it is also the one that can be defended to someone
who asks why a form was ordered.

### The absolute surface test, and why it was needed

The membrane filter was **relative** — does a form lose the TM helices its reference has? —
and that test cannot fire when the reference itself has no topology call. Twelve forms
reached the cohort that way, and they are real: `ADAM12`-S, soluble `PVR`, soluble `OSMR`,
shed `NRXN3`, `FCGR2C`, `PRSS8`. Genuine secreted isoforms of membrane proteins, and not
surface-displayable, which is what the library is for.

So s1e3's `topology_status` is now applied **absolutely as well as relatively**: a form whose
own DeepTMHMM label is `SP` or `GLOB`, and which carries no GPI anchor, is excluded whatever
its reference does. That removes 11 secreted and 1 intracellular form, taking the full tier
to **1,916** and the observed tier to **816**.

Twelve forms in 1,928 is a small error, and it would have stayed small under random
sampling. Stratification is what made it matter: a per-class floor gave the spurious
"no TM helix" class 8 slots per tier, turning 0.6% of the cohort into **8% of a 150-form
panel**. A stratified sampler amplifies any class that should not exist, so the sampler now
refuses to run if one appears.

## The balanced panel

150 forms — 100 from the observed tier, 50 from prediction-only — stratified by the
**topology class of the library gene the isoform belongs to**, not by the isoform's own
topology. An isoform of a type I receptor that drops its signal peptide computes as type III
on its own sequence, and binning on that is what produced a panel with 20 type III members.

### The library's own class distribution

Established only after running DeepTMHMM v1 over the **588 library references that had no
topology call** — a quarter of the library. The result moves materially:

| class | before the run | after |
|---|---|---|
| multipass | 807 (34.1%) | **1,223 (51.6%)** |
| type I | 760 (32.1%) | 836 (35.3%) |
| GPI | 121 | 135 (5.7%) |
| type II | 72 | **126 (5.3%)** |
| type III | **13** | **33 (1.4%)** |
| no topology call | **593** | **6** |

The missing quarter was disproportionately multipass, so the earlier figures understated the
library's dominant class by a third. **Type III went from 13 genes to 33**, which matters
because the panel's floor was set against the smaller number.

### Allocation

Floors of **2 for type III**, **5 for GPI** and **5 for type II**; everything else
proportional. Shares are
assigned by **largest remainder** rather than independent rounding — rounding each share
separately let the two big classes absorb the whole remainder and gave type II **zero**
slots despite 108 candidates.

| class | panel | library share |
|---|---|---|
| type I | 61 | 35.3% |
| multipass | 50 | 51.6% |
| **GPI** | **17** | 5.7% |
| **type II** | **17** | 5.3% |
| **type III** | **5** | 1.4% |

136 genes, 62 already in the current cohort. GPI, type II and type III all sit well above
their library share — together 39 of 150 against about 12% of the library — which is what
the floors are for. Multipass is correspondingly under its 51.6%, by design: it is the class
least likely to be under-sampled by accident.

GPI is tested **before** the TM rules, since a GPI protein's C-terminal signal is routinely
called as a transmembrane helix and these would otherwise be filed as type I. The KW-0336
keyword confirms all 113 genes in the existing GPI design set with no disagreement.

The seed is fixed and recorded, so a rerun returns the same 150.

### Source precedence, and the guard against soluble references

The reference is taken from the strongest source that offers a usable answer:

**MANE › APPRIS › observed evidence (PeptideAtlas peptides, OpenProt translation) › UniProt
canonical.**

MANE is a RefSeq/Ensembl agreement and transcript-backed by construction; APPRIS is
structural and conservation inference; observed evidence decides only where neither
annotation offers an answer; the UniProt canonical is last, because its selection rule
permits a pure length tiebreak with no transcript support. Outcome over the 2,368:

| rung | genes |
|---|---|
| MANE | 2,282 |
| APPRIS | 75 |
| UniProt canonical | 11 |

**A reference must also not be a soluble form of a surface protein.** Each rung is checked
for membrane topology before it is accepted, and skipped if the form is *explicitly* called
`SP` or `GLOB` with no GPI anchor — a form with no topology call is never rejected on
absence of evidence. The guard fires on **11 genes**:

- **`PILRA` (230 → 303), `LMLN` (452 → 647), `SUSD6` (114 → 303)** — APPRIS had moved each
  onto a soluble isoform and the guard falls through to the UniProt canonical. These are the
  concrete harm from taking APPRIS unconditionally, which is the same branch that agreed
  with Tedman 0 times in 6.
- **Eight genes where every candidate is non-surface**, flagged rather than silently kept.

The guard was initially keyed on a `sequence` column that the processed DeepTMHMM cohort
TSVs do not have, so the lookup was empty and `PILRA` kept its soluble reference anyway. The
labels are now joined back to sequences through UniProt accessions — 5,504 mapped, none
unmapped.

### The GPI set: a curated attribute, two annotation sources, and a cited supplement

A GPI-anchored protein has a signal peptide and no transmembrane helix, so DeepTMHMM calls
it `SP` — **correctly**. Without an anchor attribute those genes read as secreted and a
topology guard throws them out, so the attribute is load-bearing rather than bookkeeping.

The primary source is the architecture call already carried in the surfaceome database:
`surface_annotation.annotation_json → filters.subcategory == "GPI_anchored"`. It is
evidence-backed per gene and covers every library member, so it is **read, not predicted**.
Two further sources are unioned in because no one of them is a superset:

| source | genes | what it adds |
|---|---|---|
| database attribute | 156 | the primary call, with cited evidence per gene |
| UniProt KW-0336 | 144 | 7 the attribute misses, 3 of them in the library |
| `gpi_supplement.tsv` | 6 | literature the other two have not absorbed |
| `gpi_exclusions.tsv` | −2 | genes a source calls anchored that the evidence contradicts |

The three UniProt-only library genes are instructive rather than noise: `BST2` genuinely
carries both an N-terminal transmembrane helix and a C-terminal anchor, `NCAM1` is anchored
only as NCAM-120, and `TNFRSF10C` is a GPI-anchored decoy receptor whose anchor signal is
read as a helix.

The supplement exists because the keyword is incomplete. According to PubMed,
[Davies et al. 2010](https://doi.org/10.1073/pnas.0908735107) show GPI anchoring directly
for the voltage-gated calcium channel **alpha2delta-1, -2 and -3** subunits — release by
PI-PLC, a cross-reacting-determinant epitope, and loss of membrane localization when the
predicted anchor site is mutated — while UniProt still annotates them as single-pass type I.
An earlier version of this note had `CACNA2D3` the wrong way round, as the predictor erring
against UniProt, when the predictor was right and the annotation is stale.

Exclusion is deliberately narrow: a gene leaves only when **the sources disagree with each
other**, never on one predictor's say-so.

| gene | why it leaves | kept in the library by |
|---|---|---|
| `EPGN` | the attribute is alone — no KW-0336, and NetGPI rejects all 11 forms at 0.961–0.995 | its transmembrane helix |
| `RAET1G` | RAET1G1 encodes transmembrane and cytoplasmic domains "unique amongst NKG2D ligands" and RAET1G2 is secreted ([Eagle et al. 2009](https://doi.org/10.1371/journal.pone.0004503)) | its transmembrane helix |

### Locating the anchor on a form, not a gene

The attribute is a claim about the **gene**: that some form of it is anchored. It does not
say the form this library selected is that one, and for six genes it is not. NetGPI 1.1
([Gíslason et al. 2021](https://doi.org/10.1016/j.crbiot.2021.01.001)) is run over every
human form of every GPI gene — the selected reference, all 770 alternates rather than only
the evidence-filtered cohort, and the ortholog picks.

The published weights are not distributed: the GitHub repository is MIT-licensed training
code whose `.gitignore` excludes the `*.pt` ensemble, and DTU publishes only the dataset.
The service is therefore the only way to run the tool, and the sequences submitted are
public UniProt and Ensembl entries.

**A prediction is worth what its measured error rate on the data in hand says it is.** The
same service was run over 148 curated-GPI and 100 single-pass human references:

| | value |
|---|---|
| sensitivity | 89.9% raw, **95.0%** once form-level labelling is accounted for |
| specificity | 99.0% |
| MCC | 0.875, against the paper's 0.895 |

The raw figure is pessimistic because five of its fifteen misses are NetGPI correctly
finding the anchor on a different form, and three are bad attributes. Reading the calls
follows what that run measured:

- **Anchoring is a veto, not a sort key.** A form is excluded only on a confident negative;
  source precedence does the ranking, so a marginally better-anchored form cannot displace a
  better-sourced one.
- **A gene with no anchored form is reviewed, never dropped.** A confidence cutoff does not
  separate the failures: NetGPI calls all ten `CRIPTO` forms non-anchored at 0.745 though
  Cripto-1 is GPI-anchored ([Watanabe & Salomon 2010](https://doi.org/10.1016/j.bbrc.2010.10.128)),
  with `CFC1` missed at 0.93 and `GML` — whose HGNC name is "GPI anchored molecule like" —
  at 0.752. All three carry KW-0336, so UniProt and the attribute agree against the
  predictor. `netgpi_blind_spots.tsv` records them and they stay in the set.
- **A form under 50 residues cannot carry a signal that is 18–39 residues by itself**, so
  its negative is uninterpretable rather than evidence.

Six references move onto the form that carries the anchor:

| gene | reference, not anchored | anchor is on |
|---|---|---|
| `NCAM1` | `P13591` | `P13591-3` — NCAM-120 |
| `ACHE` | `P22303` | `P22303-2` — the H subunit |
| `TFPI` | `P10646` | `P10646-2` — TFPI-β |
| `BCAN` | `Q96GW7` | `Q96GW7-2` |
| `CD55` | the APPRIS-rung form | `P08174` — canonical DAF |
| `CD58` | `P19256` | `P19256-2` |

`CD55` is a selection fault rather than biology: DAF's canonical form is the anchored one,
and APPRIS had picked a transcript encoding the secreted form. `CD58` is named in
`gpi_form_overrides.tsv` rather than ranked — its two anchored candidates predict the same
omega site and are identical through it, so they give the same mature protein and differ
only in the cleaved signal, and the ranking preferred the Ensembl model purely because the
reviewed isoform carries no APPRIS label.

### The ortholog anchor, and why the guard was a no-op

The surface guard refused to move a reference onto a candidate that had lost every
transmembrane helix, but read `or human_symbol in gpi_human` — so on a GPI gene **every
candidate passed whatever its own sequence said**. The guard did nothing on exactly the
genes where the anchor is the only thing holding the protein to the surface.

The database records its ortholog topology as `projected_from_human_canonical`, so nothing
in hand measured it. The within-species C-terminal test used for isoforms cannot substitute:
the signal is conserved as architecture, not sequence. Human against mouse `NYX`, which
differ in anchoring ([O'Connor et al. 2005](https://doi.org/10.1093/hmg/ddi194)), scores
0.37 over the last 30 residues, while human against mouse `CD55`, which do not, scores 0.43.

NetGPI over all 761 candidates supplies the equivalent test. Six picks move onto a form that
keeps the anchor — `NCAM1` in both species, and `ACHE`, `CD55`, `CD58`, `EFNA4` in cyno —
each at a cost in correspondence, which is the same trade the helix guard already makes.
Where no candidate is anchored the existing fallback keeps the best-aligned one, so a
family-shaped blind spot cannot empty a gene.

### Library membership is a check, not an inheritance

`EPPIN` was removed from the library on 2026-10-11. Every source except the database
attribute says it carries no anchor of its own — no KW-0336, NetGPI rejects all five human
forms at 0.994–0.996, and DeepTMHMM gives the reference `SP` with no helix. It does reach
the sperm surface, but through the eppin protein complex with lactotransferrin, clusterin
and semenogelin ([Wang et al. 2007](https://doi.org/10.1095/biolreprod.107.060194)), so the
membrane attachment belongs to its partners. A surface-display library needs the protein to
carry its own anchor. `library_exclusions.tsv` records that beside the member list so a
regeneration cannot quietly restore it.

Removing it exposed a second problem. Several tables here are **snapshots** of upstream
sources and keep every gene they were built with, so a script that inherits membership from
a snapshot keeps emitting a gene the library has dropped — `EPPIN` reappeared in four
outputs that way. The cohort, form-call and slot-selection steps now read
`library_members.json` and filter on it.

## How the ortholog proteoform was chosen, and how it should be

The Compara one-to-one high-confidence mapping is **gene level** and is not in question.
What was weak is the proteoform taken inside each confirmed ortholog gene.

**It was not the canonical.** Measured against the data, the cohort took the **longest
symbol-matched entry, reviewed or not**: 654 of 1,971 mouse picks are unreviewed TrEMBL, and
for **543 of those a Swiss-Prot entry existed and was passed over**. Against the longest
reviewed proteoform the selection ran longer in 265 cases and shorter in 13. That is the same
length tiebreak that made UniProt's human canonical unreliable, without even the
reviewed-entry floor.

**What "correspondence" means here.** For a candidate sequence and the selected human
reference, `difflib.SequenceMatcher(autojunk=False)` is used to find the exactly-matching
contiguous blocks, and the score is **total matched residues divided by the longer of the
two lengths**. The candidate with the highest score wins. This follows the precedent in
`tag_sites/tedman.py::project_junction`, which uses the same approach and documents it as
*"the longest identical blocks between the two sequences rather than an alignment"*.

It is **not** sequence identity in the Needleman–Wunsch sense: there is no substitution
matrix, so a conservative I→V substitution counts as a mismatch; there are no gap penalties,
so blocks may be separated at no cost; and the block search is greedy rather than optimal.
The absolute numbers below are therefore systematically lower than a proper aligner would
report and **should not be quoted as percent identity**. For *ranking candidates within one
gene* against one reference the deflation is common to every candidate and cancels, which is
the only use it is put to.

Dividing by the **longer** length is what stops the metric preferring short candidates: a
fragment is capped by its own length over the larger denominator, so a 200-residue candidate
against a 600-residue reference cannot score above 0.33.

Four strategies, scored this way against the **selected human reference**, paired on the
1,296 orthologs where all four produce a candidate:

| strategy | median score | mean | ≥0.9 | better than current | worse |
|---|---|---|---|---|---|
| current (longest symbol match) | 0.8568 | 0.8317 | 447 | — | — |
| reviewed (Swiss-Prot canonical) | 0.8630 | 0.8381 | 477 | 233 | **54** |
| APPRIS (mouse principal) | 0.8637 | 0.8379 | 485 | 266 | **66** |
| **correspondence (best match to the human reference)** | **0.8656** | **0.8401** | **490** | **275** | **0** |

The medians separate by about a point, which is not the argument. The argument is the last
two columns: **reviewed and APPRIS each make the pick *worse* on 54 and 66 genes**, because
neither optimises the property that matters. Correspondence cannot, by construction — it
selects on resemblance to the human form the ortholog exists to partner, which also makes
"human principal and mouse principal are different splice forms" impossible rather than
merely detectable.

**APPRIS beats correspondence on 58 genes**, which is informative rather than contradictory:
APPRIS picks from Ensembl, and those candidates were not in the UniProt-derived pool scored
here. Widening the ortholog candidate pool to Ensembl should close that gap, exactly as it
did on the human side.

### The re-pick, as run

`repick_orthologs.py` applies correspondence with the candidate pool widened to Ensembl.

| | mouse | cynomolgus |
|---|---|---|
| orthologs re-picked | 1,452 | 1,429 |
| **pick changed** | **481** | **362** |
| median score, before → after | 0.854 → 0.864 | 0.947 → **0.961** |
| median length ratio to the human reference | 0.9961 → 0.9972 | 0.9957 → **1.0000** |
| gain above 5 points | 92 | **209** |
| **reviewed picks** | 949 → **1,227** | 53 → 58 |
| picks sourced from Ensembl | 116 | **326** |
| ever worse than the current pick | **no** | **no** |

Cynomolgus gains most, which follows from where it started: 96% of its picks were unreviewed
TrEMBL and it has no APPRIS, so correspondence is the only signal available and the baseline
was weakest. Its median length ratio reaching exactly **1.0000** means the typical re-picked
cyno ortholog is now the same length as the human form it partners.

Widening the pool to Ensembl was not cosmetic: **326 of the 362 changed cyno picks and 116
of the 481 mouse picks come from Ensembl translations** that no UniProt-only search would
have found. This is the same finding as on the human side, where Ensembl supplied 79% of
alternate forms.

565 mouse and 492 cyno orthologs are not re-picked. **None is a current library member** —
they are orthologs of genes that left the library, carried along in a cohort built against
an older membership list.

The re-pick is applied: **correspondence first, reviewed as a tiebreak among equally
correspondent candidates, APPRIS as a further tiebreak for mouse.** Cynomolgus gets
correspondence only — it has no APPRIS, no PeptideAtlas and no OpenProt, and 96% of its
current picks are unreviewed. 1,452 mouse and 1,429 cyno picks stand, every one carrying a
topology call.

Candidate order is sorted before scoring. The translations arrive as a set, so iterating
them followed the process hash order and the tiebreak between equally-scoring candidates
differed from run to run — `TMIGD2` flipped between two candidates both scoring 0.8489, and
seed 5 disagreed with seeds 1–4. The same latent instability sits in any other caller of
`translations_by_gene`.

### The ortholog panel

Orthologs need the same stratification the isoform panel needs, for the same reason, plus
one more: an ortholog's value for tag placement depends on how far it has drifted from the
human reference, and an unstratified draw returns mostly near-identical sequences that
exercise nothing.

So the draw is stratified twice — by the **human gene's class**, never the ortholog's own
topology, and by **correspondence**, split at the measured median into a near and a diverged
band. Only surface picks are eligible, and on a GPI gene only picks that keep the anchor, so
the panel cannot propose an ortholog the library would reject.

| | mouse | cyno |
|---|---|---|
| eligible | 1,448 | 1,413 |
| correspondence split | 0.864 | 0.962 |
| drawn | 150 | 150 |

300 picks over 288 human genes, 75 per band per species, with the same class floors as the
isoform panel. The first draw came from a third of the orthologs because most carry their
topology in the database rather than in any local run; `ortholog_topology_labels.tsv`
snapshots all 13,169 so the pool is the whole set.

## Membership tiers and construct sets are different things

These are the two layers readers most often conflate, because both get called "rescue".

A **tier** says why a gene is a library member. Every one of the 2,368 members carries
exactly one, and three of the six are rescues in the sense of admitting a gene the main
criteria would have excluded.

| Tier | Members | Admitted because |
|---|---|---|
| T1 | 967 | low architecture (GPI or single-pass) and constitutive expression |
| T2 | 880 | medium architecture (multipass, non-olfactory GPCR) and constitutive |
| T3 | 286 | high architecture (olfactory or taste), or conditional expression |
| LL | 139 | **low literature** — under ~100 papers, so no independent surface evidence; admitted on the optimized UniProt-surface cutoff alone |
| SV | 85 | **screen validated** — Tedman's GPCR screen placed it at the surface |
| OE | 12 | **overexpression evidence** only |

A **set** says what kind of construct a manifest row is: `library_reference`, `isoform`,
`ortholog`, `expression_rescue` or `giant_truncation`. A gene can appear in several.

So `expression_rescue` is **not** a tier and adds no members. All 25 of its genes already
hold a `library_reference` construct; the rescue row is a *second* construct for the same
gene, carrying an exogenous Igκ leader. Their tiers are T2 10, T3 10, LL 5 — and no SV,
which is structural rather than coincidental: SV means Tedman measured the gene *above* his
expressed bar, and the rescue arms are drawn from below it. One dataset, opposite ends.

## Expression-rescue arms

25 receptors get an Igκ leader and 10 of those are also fully recoded. The draw is entirely
Tedman-derived and therefore **GPCR-only**: there is no equivalent screen to threshold
against outside the GPCRome, so nothing else in the library has a rescue arm. That is a
scope limit worth stating, because the name suggests otherwise.

| Step | Count |
|---|---|
| Tedman expressed bar | 1025.2 a.u. |
| Library members below the bar | 276 |
| Eligible after topology filters | 244 |
| Leader arms drawn | 25 (10 also recoded) |

Two filters removed candidates that would have been mistakes: 30 already carry a signal
peptide, and one has a cytoplasmic N-terminus where a leader would invert the topology and
throw the sorting determinants outside. The draw is then balanced across GPCR classes
(A 7, olfactory 7, taste 7, B 2, vomeronasal 1, C 1) rather than ranked straight, because
olfactory receptors are 17 of the 30 lowest expressors and a straight ranking would return
an olfactory experiment that says nothing about the rest.

The leader itself is the murine Ig κ-chain V-J2-C signal peptide, read off the annotated
vector sequence in Thermo's pSecTag2/Hygro manual and archived. **Fragmid §4.1 specifies
the influenza HA signal peptide for N-out receptors and reserves Igκ for ectodomain-only
constructs**, so these 25 full-length arms sit on the far side of a line Fragmid drew the
other way. That is unresolved, and it is a decision rather than an oversight.

## Giant truncations

24 constructs. Members whose ORF exceeds what the vector carries are not dropped; a
truncation is specified instead, with the retained span chosen to keep the surface-exposed
region. They are kept out of the first synthesis tranche because a truncation is a design
judgement per gene rather than a rule, and the tranche is meant to be the part the tree
decided on its own.

## The orderable unit, and the two synthesis tranches

A manifest row is one orderable insert. Ids are `l1-pNNNN` and `l2-pNNNN`, assigned by
`split_library_tranches.py` from (set, key) with no randomness, so a re-run on the same
inputs is byte-identical and `prior_candidate_id` keeps the previous numbering recoverable.

| Tranche | Contents | Count |
|---|---|---|
| **L1** | every human reference, the 25 rescue arms, the isoform panel, the ortholog panel | **2,725** |
| **L2** | isoforms and orthologs outside the panels, plus the giant truncations | **4,341** |

**Renumbering is safe only before barcodes are reserved.** `tag_cassette` is `deferred` and
`barcode` is `not reserved` on all 7,066 rows, so no id is yet bound to a physical reagent.
Once a barcode is reserved the id becomes the key a sequencing readout joins on, and moving
it silently invalidates every downstream map.

### Joining a panel to the manifest

The manifest keys an isoform by **whichever identifier resolved it** — a UniProt accession
for some (`iso:ACHE:P22303-2`), an Ensembl transcript for others
(`iso:CACNA2D1:ENST00000957014`), because not every form has an accession. A join on
transcript alone recovers 62 of 150 panel rows; accepting either identifier recovers 144.
Ten rows still do not resolve and are listed in `library_tranche_summary.json` rather than
smoothed over: six have no construct in the manifest at all (LYPD1, NRG4, cyno ATP1B1,
mouse LYPD6B, PTPRJ, INSR) and four name a proteoform the manifest does not hold.

The ortholog join **refuses a (species, gene) fallback on purpose.** Several proteoforms of
one ortholog gene can be present, and taking the wrong one substitutes a construct the panel
never drew — the same silent wrong-protein failure the identifier discipline exists to
prevent. It reports the mismatch and the siblings it found instead.

## What is deliberately not resolved yet

| Held | Why |
|---|---|
| `tag_cassette: deferred` | Placement is unresolved at library scale, so the insert runs Kozak → ORF and the cassette slot is left empty. `MIN_LOOP_LEN`, `LOOP_SS` and the cargo size each change which site a member gets |
| `barcode: not reserved` | Reserving ~7,000 against a shared append-only registry before the inserts are final would bind them to sequences that are about to change |
| 8 references unbuildable | Pseudogenes with a premature stop (FCGR2C, OR52P1, OR8J2), no clean transcript (KIR2DS4), and four olfactory receptors whose transcripts sit a few residues from the selected form. They pass reference selection with `needs_review: False` and fail only at the CDS layer, so the selection table does not flag them |

## Reproducing an output

Every object published through `processed.upload()` carries the commit, a dirty flag, the
script and the lock hashes as S3 metadata. Read it rather than assume:

```
uv run python -m deliverome_analysis.provenance processed://<key>
```

`candidate_manifest_tranches.tsv` returns commit `086574d`, `dirty=false` and the script
that produced it, so `git checkout 086574d && uv sync --frozen` reproduces it. **Most older
keys under this prefix return "no provenance recorded"** — `candidate_manifest.tsv`,
`ortholog_panel.tsv` and even `constructs/designed/p0022.gb` — despite being written after
stamping landed, so something in their write path bypassed the helper. Treat an unstamped
key as untraceable until it is re-published.
