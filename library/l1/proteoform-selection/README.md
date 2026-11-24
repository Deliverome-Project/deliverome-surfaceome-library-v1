# Selecting a reference isoform with APPRIS, and what it disagrees with

> **Superseded. This page is the APPRIS-vs-UniProt comparison, not the selection policy.**
> Its recommendation below — keep UniProt as the build sequence — was **not** what we
> shipped. [`select_reference.py`](select_reference.py) is the policy of record, and it
> changed the reference for **218** members:
>
> 1. MANE Select and the APPRIS principal **agree** → take it (**2,171** members).
> 2. Only **one** of the two exists → take that one (**77**, all APPRIS-only).
> 3. They **disagree** → adjudicate on *observed* evidence (PeptideAtlas peptides,
>    GENCODE Ribo-seq ORFs), not on annotation (**111** members). In the event the
>    evidence discriminated in **0** of them — 51 had a Ribo-seq ORF but none overlapped
>    a region that distinguishes the candidates — so all 111 defaulted to MANE.
>
> Final build source: MANE **2,277**, APPRIS **74**, UniProt canonical **11**, named **6**.
> So UniProt canonical is the build sequence for 11 members, not the default.
> Agreement is tested on the **protein sequence**, never the transcript ID.
>
> The directory name says `appris-principal` for historical reasons; the ladder is
> MANE → APPRIS → experimental evidence → UniProt. Read `select_reference.py`'s module
> docstring before using anything on this page to answer "which proteoform did we pick".

**The recommendation: adopt APPRIS as an annotation, not as an override.** Its principal
isoform already *is* our build sequence for 2,085 of 2,369 library members (88%). For the
268 where it is not, switching would invalidate the topology, cleavage and disorder tracks
we computed on the UniProt sequence, in exchange for a different answer to a different
question. Keep UniProt as the build sequence, carry the APPRIS call as a column, and review
the 17 N-terminal disagreements by hand, because those are the ones that move a tag site.

**Do not apply APPRIS independently to the mouse orthologs.** Per-species selection is what
creates the failure mode this analysis was built to detect, and it fires: 363 mouse
principals are a UniProt isoform we do not store, and 27 human/mouse principal pairs differ
in length by more than 10%.

## Terms

| term | meaning |
|---|---|
| **APPRIS** | A per-transcript annotation of which splice form is the main one, inferred from protein structure, functional residues and cross-species conservation. Release **e115v50**, built on Ensembl 115. |
| **PRINCIPAL:1** | APPRIS's strongest call: one CCDS-supported transcript whose structural and functional evidence agree. Weaker classes **PRINCIPAL:2–5** and **PRINCIPAL:M** are assigned when the evidence does not single one transcript out. |
| **MANE Select** | One transcript per human gene on which RefSeq and Ensembl agree exactly. Human only — the mouse APPRIS table has no MANE column. |
| **principal** (here) | The transcript this analysis selects for a gene: the best APPRIS class present, collapsed by protein sequence. |
| **stored** | A UniProt sequence our pipeline actually holds — the DeepTMHMM cohort FASTAs in `accessible-surfaceome`. |

## The selection rule

1. Take the **best APPRIS class present for the gene**, never a mix of classes. 1,528 of
   2,369 library members reach PRINCIPAL:1; the rest are carried at the class they reach,
   recorded in `appris_class`, rather than dropped. Requiring PRINCIPAL:1 would silently
   discard 36% of the library.
2. Collapse the transcripts in that class **by protein sequence**. Several transcripts
   commonly encode one protein and differ only in UTR, which is not an ambiguity.
3. If distinct proteins remain, break the tie on MANE Select, then on length, and flag the
   gene. Only **CSF2RA** and **IL3RA** are ambiguous this way — both pseudoautosomal, both
   with two distinct PRINCIPAL:2 proteins.
4. Join on stable IDs throughout. The library symbol selects a stored accession; the
   accession carries to Ensembl through the UniProt cross-reference. 29 genes have no
   cross-reference and fall back to the APPRIS gene symbol, marked `ensg_source`.

## What APPRIS points at, per library member

| outcome | human (of 2,369) | mouse (of 1,971 orthologs) |
|---|---|---|
| the canonical we already build from | **2,085** | **1,383** |
| an isoform we already store | 139 | — (no mouse isoform cohort) |
| a UniProt proteoform we do not store | 34 | 363 |
| **no reviewed UniProt sequence for that gene** | **95** | **225** |
| no APPRIS entry for the gene | 9 | 0 |
| no UniProt canonical for the library symbol | 7 | — |

2,165 human principals are also MANE Select.

The 9 with no APPRIS entry are KIR and olfactory-receptor loci (`KIR2DL2`, `KIR2DL5A`,
`KIR2DS2`, `KIR2DS5`, `KIR3DS1`, `OR8U8`, `OR8U9`, `OR9G9`, `TAS2R45`) — genes that sit on
alternate haplotypes rather than the primary assembly APPRIS annotates. The 7 with no
UniProt canonical are stale symbols in the library list (`C14orf180`, `C19orf38`,
`C3orf52`, `C5orf15`, `OR2A1`, `OR2A42`, `OR4F16`); resolving the library through HGNC
would clear them.

## The principals that match nothing we store

`appris_principal_not_stored.tsv` holds all 717 (129 human, 588 mouse). The `status` column
separates the two fixes: `uniprot_proteoform_not_stored` is a cohort gap we can fill by
fetching the sequence, `no_uniprot_match` means APPRIS and UniProt disagree about what the
protein is. For the latter, `difference_kind` says where the difference falls, because that
decides whether it touches a design:

| difference | human | mouse | consequence for a design |
|---|---|---|---|
| `internal` — a splice difference | 65 | 105 | changes the ectodomain; may move an internal tag site |
| `n_terminal` | 17 | 15 | **moves the signal peptide, and every post-SP placement derived from it** |
| `c_terminal` | 1 | 2 | moves the GPI omega site and the C-terminal tag arm |
| `substitution_only` | 12 | 103 | none, if the substitutions are allelic |

**The substitution-only cases are mostly allele differences, not different proteins.**
UniProt's reference entry and the reference genome carry different alleles at common
polymorphisms, and UniProt annotates those positions. Five human cases (`F2RL1`, `GPR42`,
`KIR3DL3`, `OR4E2`, `PVR`) have *every* differing position annotated as a natural variant
whose alternative allele is exactly the Ensembl residue — PVR I340M and OR4E2 V118M/Q234R
were checked by hand against the UniProt feature table. The loci that dominate this bucket
(KIR, LILR, MICA, olfactory receptors) are the famously polymorphic ones, which is what a
reference-allele difference should look like.

The mouse column cannot be read the same way: **95 of the 103 mouse accessions carry no
annotated natural variant at all**, against 7 of 12 human accessions carrying at least one
(median 2.5). A mouse verdict of `none` therefore means UniProt does not annotate mouse
strain variation, not that the difference is real. Mouse substitution-only cases are
unresolved, not refuted.

## Three reference choices for the same gene: ours, Tedman's, APPRIS's

Every construct starts from a choice of which sequence the gene *is*, and three are on the
table: our UniProt canonical (what `build_constructs.py` splices into), the Ensembl
transcript Tedman built each screen-validated GPCR construct from, and the APPRIS principal.
`compare_reference_choices.py` puts all three on his 766 GPCRs — the only gene set where an
independent, experimentally validated choice exists.

**They are not very different.** Agreement is sequence identity, not similarity:

| comparison | same | differ | agreement |
|---|---|---|---|
| Tedman vs ours | 732 | 33 | **95.7%** |
| APPRIS vs Tedman | 723 | 40 | **94.8%** |
| APPRIS vs ours | 698 | 64 | **91.6%** |

Tedman's choice is closer to ours than APPRIS's is, and APPRIS is the outlier of the three.
None of his transcripts have been retired from Ensembl 115.

### Where they differ, by how much

| | Tedman vs ours | APPRIS vs ours | APPRIS vs Tedman |
|---|---|---|---|
| differ in length | 28 | 53 | 32 |
| differ at the **same** length (substitutions) | 5 | 11 | 8 |
| the alternative is **shorter** | 23 of 28 | 47 of 53 | 27 of 32 |
| median \|Δ\| | **14 aa** | **16 aa** | **16 aa** |
| 90th percentile \|Δ\| | 35 aa | 50 aa | 51 aa |
| largest \|Δ\| | 102 aa | 102 aa | 88 aa |
| within 5 aa | 3 | 7 | 8 |

The differences are not small adjustments: the distribution is bimodal, with a cluster at
5–35 aa and a tail past 100, and almost nothing in the 1–4 aa band that a numbering-offset
artefact would produce. **Our canonical is nearly always the longest of the three** — UniProt
selects canonical partly on completeness, APPRIS on structural evidence.

### And the difference is N-terminal, which is where the tag goes

| difference falls at | Tedman vs ours | APPRIS vs ours | APPRIS vs Tedman |
|---|---|---|---|
| **n_terminal** | **23** | **35** | 13 |
| internal | 4 | 16 | 17 |
| substitution_only | 5 | 11 | 8 |
| c_terminal | 1 | 2 | 2 |

These are GPCRs, so a terminal-N tag is the design. An N-terminal reference difference moves
the tag site by exactly that many residues — it is the worst place for the three sources to
disagree, and it is where they do.

**30 genes have both Tedman and APPRIS differing from us, and on 26 of those the two agree
with each other** — two independent sources against our canonical. 22 of the 26 are
N-terminal, and all but two are olfactory receptors, whose single coding exon makes the
start codon a prediction rather than an observation. Our canonical runs 4–47 aa longer than
both. **Nine are current library members**: `NPFFR2`, `OR11H1`, `OR13D1`, `OR2AK2`, `OR3A2`,
`OR3A3`, `OR4F5`, `OR51V1`, `OR8G5`.

The non-OR disagreements are few enough to name, and are real annotated alternative products
rather than annotation noise:

| gene | ours | APPRIS | Tedman | what APPRIS picks |
|---|---|---|---|---|
| `NPFFR2` | 522 | 420 | 420 | isoform `Q9Y5X5-2`, which we store — but at PRINCIPAL:4, a weak class |
| `GPR142` | 462 | 374 | = ours | isoform `Q7Z601-2`, which we do not store |
| `HRH2` | 359 | 422 | 363 | a C-terminally extended form no reviewed UniProt entry carries |
| `GPR17` | 367 | 339 | = ours | the short isoform `Q13304-2`, which we store |

For 12 of the 17 non-OR disagreements **Tedman sides with us**, which is the main reason to
treat APPRIS as an annotation rather than an override.

## Testing one disagreement against observed peptides (PeptideAtlas)

Sequence annotation cannot settle a start codon; observed mass spectra can speak to it.
**PeptideAtlas** build **Human 2026-01** (`atlas_build_id` 607, 5.3M distinct peptides) is a
uniform reprocessing of public MS data, so it is orthogonal to both UniProt and APPRIS.

`NPFFR2` is the test case: the one non-olfactory gene where Tedman and APPRIS both disagree
with our reference *and* agree with each other. We build 522 residues, both build 420,
differing by a 102-residue N-terminal extension. If that extension is not translated, our
terminal-N tag sits on 102 residues that do not exist.

PeptideAtlas calls `Q9Y5X5` **canonical** at ProteinProphet 0.999, with **8 distinct
peptides over 15 observations**, on the 522-residue sequence — so its coordinates are ours.
`peptideatlas_npffr2.py` asserts that rather than assuming it, and every peptide matches.

**None of the 8 observed peptides falls in residues 1–102.** The first observation is at
residue 310.

That is consistent with the shorter form, and it is **not enough to act on**:

| test | result |
|---|---|
| detectable tryptic peptides in the extension that were observed | **0 of 4** |
| same, from residue 310 on (where coverage exists) | 4 of 10 |
| Fisher exact, one-sided | **p = 0.21** |
| binomial: 0 of 8 peptides in a region that is 19.5% of the sequence | **p = 0.18** |

The extension would yield only four fully-tryptic peptides in the 7–30 aa detectable window
(`MNSFFGTPAASWCLLESDVSSAPDK`, `ALSVQQR`, `GGPAWSGSLEWSR`, `AWWILVPAADR`), and even the
well-covered C-terminal half gives up only 4 of its 10. Zero out of four in a protein seen
15 times total cannot be separated from the region simply not having been sampled. Two tests
on different assumptions agree on that.

**So PeptideAtlas is worth consulting and did not resolve this one.** It would resolve a
case with real depth of coverage; a low-abundance GPCR seen 15 times is not that case, and
the olfactory receptors that make up the other 21 N-terminal disagreements will be worse.
What would settle `NPFFR2` is an N-terminomics dataset (N-Cap/TAILS) that reports the actual
protein N-terminus, rather than inferring it from the absence of tryptic peptides.

## Prioritising alternate isoforms by evidence that they are expressed

Isoform membership today rests on UniProt having **annotated** the isoform, which says a
splice form was described, not that it is made. `isoform_expression.py` scores all **2,802
UniProt isoforms across the 1,223 library genes that have one** on observed protein
evidence: peptides in PeptideAtlas Human 2026-01 that are present in that isoform and in
**no other proteoform of the same gene**, weighted by how many distinct samples saw them.

| tier | definition | isoforms | of which in the current cohort |
|---|---|---|---|
| **expressed** | ≥2 diagnostic peptides, best seen in ≥10 samples | **142** | 117 |
| likely expressed | 1 such peptide, ≥10 samples | 12 | 11 |
| weak | diagnostic peptides, but all shallow | 342 | 290 |
| not detected | could have been seen, was not | 1,516 | 1,227 |
| **undetectable** | no peptide could ever distinguish it | **790** (28%) | 686 |

The top of the list is recognisable, which is the best available check that the method
works: `BSG`-2 (the two-Ig basigin/CD147 form) with **44 diagnostic peptides across 1,693
samples**, `PCDH1`-2 with 107, `NECTIN2`-2 with 62, `INSR`-2 — the IR-A form — with 29, plus
`ALCAM`, `L1CAM`, `NRCAM`, `PTK7`, `CLSTN1`.

### What to do with it

**26 isoforms scored `expressed` are not in the current cohort.** Two are `HLA-A`-2 and
`HLA-C`-2, which the cohort excludes on purpose (it is the `non_hla` set), so **24 are real
additions** — led by `SLC3A2`-1 (the long CD98hc form, 89 peptides in 238 samples),
`TNFSF12`-2 (80 peptides), `TMC4`-1, `SPINT1`-2, `LSR`-3 and `CD99L2`-5.

**1,227 cohort isoforms have no diagnostic peptide despite one being possible.** That is a
ranking signal, not a deletion list: PeptideAtlas is deep but not complete, and a surface
protein restricted to a tissue nobody profiled will look silent. Treat it as the bottom of
the priority order rather than as evidence of absence.

**790 isoforms (28%) cannot be evidenced by bottom-up MS at all** and no amount of data will
change that — see below.

### Two mistakes this took to get right

**Most isoforms have no unique residues, only a unique junction.** Alternate isoforms are
usually cassette-exon skips, so they are strict subsequences of the canonical — every
residue they contain is also in the longer form. CD44 has 17 of 19 isoforms like this, all
7 of `PTPRC`'s, both of `KCNN3`'s. What is novel is the **join** the skip creates. A
per-residue uniqueness test erases exactly that, because the residues either side of the
junction are individually present elsewhere; the test has to be **per window**. Fixing it
took `CD44` from 0 to 7 of 18 isoforms with something findable, `EGFR` to 3 of 3.

**PeptideAtlas indexes isoforms separately, and the first run did not ask for them.**
Querying `P16070` returns 913 peptide rows; `P16070-3` returns a different 865. Asking only
for canonical accessions hid precisely the isoform-specific peptides being looked for, and
the run found **3 expressed isoforms in 2,802**. Querying all 4,025 proteoforms gives 154.

A note on UniProt's format, since it looks wrong at first: the stream lists the displayed
sequence once under the bare accession and every **non-displayed** isoform as `-N`, never
duplicating — in 0 of 10,588 multi-isoform entries does the canonical reappear as a `-N`
row. So a `-1` row (365 entries) means isoform 1 is *not* the displayed sequence, which is
why `SLC3A2`-1 at 630 aa is a genuine alternate to a 529 aa canonical rather than a copy of it.

## APPRIS ranks the forms no experiment can reach

Observed evidence runs out fast: of 13,687 alternate forms, **11,631 have neither a usable
peptide nor an OpenProt translation score**, and 5,099 of those cannot be evidenced by
bottom-up MS even in principle. APPRIS annotates **every** transcript, not only the
principal, so it can order exactly that bucket. `isoform_appris.py` joins
`appris_data.appris.txt` onto the universe: **14,894 of 16,048 forms match**, 12,546 carry a
TRIFID score.

### It agrees with the evidence, which is why it is usable

TRIFID is APPRIS's machine-learning estimate of isoform functionality. Against our
independent observed evidence:

| evidence tier | forms | median TRIFID |
|---|---|---|
| expressed | 40 | **0.683** |
| likely expressed | 8 | 0.373 |
| weak | 189 | 0.281 |
| undetectable | 4,902 | 0.052 |
| not detected, low gene coverage | 241 | 0.036 |
| not detected | 7,166 | 0.033 |

Monotonic across every tier, median **0.252 for evidenced forms against 0.036 for the rest**
(Mann–Whitney p = 3e-120, **AUC 0.731**). A prediction that recovers the experimental
ordering on the forms we *can* measure is the only justification for trusting it on the ones
we cannot.

The label classes order the same way — PRINCIPAL 1.000, ALTERNATIVE:1 0.709,
ALTERNATIVE:2 0.237, **MINOR 0.025** — and **10,736 of the 13,687 alternates (78%) are
MINOR**, meaning APPRIS's core modules discarded them before tie-breaking even began. That
is the single biggest deprioritisation signal available.

Of the 11,631 forms experiment cannot reach, **1,119 score TRIFID above 0.5** and are worth
a look; the rest are predicted non-functional and can sit at the bottom of the order.

### A filter that matters more than the ranking — and a correction to it

An isoform that loses its membrane anchor cannot be surface-displayed at any score, so this
is an exclusion rather than a ranking. **It must not be built on APPRIS's THUMP column.**
Checked against known topologies, THUMP scores `ADRB2` — a seven-transmembrane GPCR — as
**0** helices, and returns 3 for `EGFR`, 5 for `CD44` and 5 for `TFRC` where each has one.
It is an APPRIS-internal scoring feature, not a topology predictor. DeepTMHMM, this
project's topology of record, returns 7, 1, 1 and 1 for the same four.

Rebuilt on DeepTMHMM, and splitting the two cases rather than lumping them:

| | forms |
|---|---|
| **lose every TM helix → exclude** | **372** |
| lose some but keep at least one → **keep** | 295 |
| forms with a DeepTMHMM call at all | 2,254 |

An earlier version of this note said 2,610 on THUMP, roughly seven times too many. THUMP
itself called >0 helices for 10 of the 372 DeepTMHMM says have none, so the two disagree in
both directions.

The 372 are recognisable as soluble forms, which is the check that matters: `FLT1`-2 at 687
residues against a 1,338-residue reference is **sFlt-1**, the secreted decoy VEGFR1, and the
list also holds `CEACAM1`-4, `ROR1`-3, `CADM1`-2, `EPHA10`-2 and `ANTXR1`-3. Losing *some*
helices is left as a designer's judgement: 295 forms, almost all of them undetected anyway.

**Coverage limit.** DeepTMHMM has been run on the UniProt cohort only, so only 2,254 of
13,687 forms can be tested — every Ensembl-only form is untested, and all 372 exclusions are
necessarily cohort members. Extending the membrane test to the Ensembl-only forms needs a
DeepTMHMM run over them.

### How the three sources divide the work

| | covers | forms |
|---|---|---|
| PeptideAtlas peptides | forms with a distinguishing peptide, deeply observed | 68 expressed + 175 weak |
| OpenProt TE | translation of a transcript, no peptide needed | 748 more |
| **APPRIS TRIFID** | **everything else, by prediction** | **11,631** |

Evidence first, prediction last, and the two kept in separate columns of
`isoform_universe_scored.tsv` rather than blended into one score.

## How well covered are the references, and can OpenProt pick between candidates?

### Coverage of the 2,361 selected references

| | references |
|---|---|
| **any peptide coverage at all** | **2,103 (89.1%)** |
| ≥10% of the sequence covered | 1,984 (84.0%) |
| ≥25% | 1,787 (75.7%) |
| ≥50% | 1,381 (58.5%) |
| ≥75% | 723 (30.6%) |
| **zero coverage** | **258 (10.9%)** |

Median coverage is **58.8%**. On the OpenProt side, 88.1% have a row, **59.7% score TE > 0**
(ribosome-profiling evidence of translation), 28.4% score TE = 0 and 11.9% have no row at
all — which, given OpenProt's 2022 build, may only mean the form is newer.

**Of the 258 references with no peptide coverage, only 8 have OpenProt translation
evidence.** So roughly 250 of the library's references — about one in ten — are proteins
that neither mass spectrometry nor ribosome profiling has ever seen. They are the expected
population: orphan and olfactory receptors, `C*orf` genes, taste receptors, tissue-restricted
channels. Nothing is wrong with building them, but no evidence-based method will ever rank
them, and they should not be penalised for it.

### OpenProt answers questions the peptide test structurally could not

The peptide test failed on 82 of the 111 MANE/APPRIS disagreements, mostly because where one
candidate is a subsequence of the other, **no peptide can tell them apart**. OpenProt does
not have that problem: it scores each form as its own row keyed on its transcript, so both
candidates get independent MS and TE scores even when one contains the other.

| OpenProt verdict on the 111 | genes |
|---|---|
| favours MANE (TE margin ≥2) | 32 |
| favours MANE on MS | 6 |
| favours APPRIS (TE margin ≥2) | 11 |
| favours APPRIS on MS | 4 |
| tied | 47 |
| only one candidate has a row / no data | 11 |

**It calls 30 of the 82 the peptides left unresolved** — 22 for MANE, 8 for APPRIS. Several
are decisive rather than marginal: `GYPC` TE 23 vs 1, `TMEM181` 21 vs 7, `EFNA4` 16 vs 0,
`PMEPA1` 16 vs 0, `TYRO3` 15 vs 1, `NEO1` 12 vs 0, `VLDLR` 12 vs 2. On the other side
`CD46` goes 4 vs 11 and `ROBO1` 6 vs 8, both for APPRIS.

**Where both methods spoke they mostly agree**: 18 concordant against 3 discordant. That is
a genuine cross-check, since the two rest on different data — PeptideAtlas peptides in one
case, aggregated translation studies in the other.

**And it partly closes the start-codon gap.** Five of the 13 N-terminal disagreements
flagged earlier as needing initiating-ribosome data — `GYPC`, `PMEPA1`, `TYRO3`, `EVI2A`,
`SLC1A2` — are now called by OpenProt, all in favour of MANE. The remaining eight still want
GWIPS-viz.

47 tied and 11 with no usable row means OpenProt leaves 58 of the 111 open, so this
supplements the default rather than replacing it.

## The full form universe: UniProt + Ensembl, with MS and Ribo-seq evidence

Tedman's isoform arm showed a UniProt-only list is incomplete — 85 of his 177 non-canonical
constructs match no UniProt proteoform. `isoform_universe.py` pools **UniProt and Ensembl
115**, deduplicates by protein sequence, and scores every form. `openprot_evidence.py` then
attaches OpenProt's two evidence axes.

| | forms |
|---|---|
| total distinct forms across 2,369 genes | **16,048** |
| of which the selected reference | 2,361 |
| **alternate forms** | **13,687** |
| …**Ensembl-only** | **10,880 (79%)** |
| …in both | 4,296 |
| …UniProt-only | 872 |

**Ensembl contributes seven times what UniProt does.** A UniProt-derived isoform list was
showing roughly a seventh of the candidate forms.

### "Not detected" now carries the gene's coverage

A form with no diagnostic peptide in a protein that is barely seen says nothing. Sequence
coverage of each gene's reference form is computed from the same observed peptides and
carried per row. **Median gene coverage is 58.8%**, and only **380 of 2,369 genes fall under
10%** — so for most genes an undetected form really is a negative, and only 294 forms are
set aside as `not_detected_low_gene_coverage`. That is a more reassuring answer than
expected, and it is why the tier is worth keeping separate rather than assuming the worst.

### OpenProt adds translation evidence where peptides cannot reach

OpenProt's **TE score** counts studies with a significant detection of *translation* —
ribosome profiling on annotated forms, which is what the GENCODE ORF catalogue could not
provide. **7,850 of 16,048 forms match an OpenProt row** (7,213 by Ensembl protein ID, 637
by UniProt accession) and **2,194 have TE > 0**.

The join only works because OpenProt carries **versioned** accessions (`ENSP00000307822.7`,
`NP_004876.3`) while our side is unversioned; keyed on the versioned spelling alone it
matched **zero** Ensembl proteins, which is how the first run reported 531 matches instead
of 7,850.

**678 alternate forms that the peptide tiers could not call have Ribo-seq evidence of
translation** — 350 of them tiered `undetectable`, meaning MS *cannot ever* distinguish
them, so this is the only evidence available. **300 of the 678 are Ensembl-only**, invisible
to a UniProt-based pipeline entirely. The strongest are `BSG`-2 (TE=47), `CD55` (47),
`APLP2` (46), `SLC6A6` (42, Ensembl-only), `APP` (41), `CD9` and `SLC3A2` (30, Ensembl-only).

### Combined evidence over the 13,687 alternate forms

| evidence | forms |
|---|---|
| protein evidence (diagnostic peptides, deep) | 68 |
| **translation evidence only** (Ribo-seq, no usable peptide) | **748** |
| weak protein evidence | 175 |
| gene too thinly covered to say | 289 |
| cannot be evidenced by MS at all | 5,099 |
| no evidence | 7,308 |

So **991 forms carry some positive evidence**, and three quarters of those would have been
missed on peptides alone.

### A flag on the references themselves

Of the 2,361 selected references, 2,081 have an OpenProt row and **671 of those score TE=0**
— no ribosome-profiling support for the form we intend to build. That is worth a look, but
not alarm: OpenProt 2.2 predates Ensembl 115 by several years and its TE score depends on
whether a tissue expressing the gene was profiled.

**Currency caveat, which applies to everything above.** OpenProt 2.2 is built on Ensembl
GRCh38.106 and UniProt 2022_06. 8,198 of our forms have no OpenProt row at all, and a form
absent from a 2022 build may simply postdate it. Absence of an OpenProt score is not
evidence of anything.

### neXtProt is not an option

neXtProt reached end of life and no longer serves data; its API redirects to an Expasy
archive page. The final release is September 2023, frozen on Zenodo
([10.5281/zenodo.14163587](https://doi.org/10.5281/zenodo.14163587)). It is not worth
importing: its MS evidence came from PeptideAtlas, which we query live and current, and its
curated protein-existence levels answer whether a *protein* exists rather than which
*isoform* is made.

### Checking the tiers against Tedman's isoform arm

Tedman built **177 deliberately non-canonical constructs** across 108 genes — an independent
opinion about which isoforms are worth testing. Mapping each to a UniProt proteoform:

| what his construct turned out to be | constructs |
|---|---|
| **matches no UniProt proteoform at all** | **85 (48%)** |
| maps to a scored library isoform | 79 |
| encodes the UniProt *canonical* despite his non-canonical flag | 10 |
| UniProt isoform, gene not in the library | 2 |
| no Ensembl 115 peptide | 1 |

**Nearly half the isoforms he thought worth building are not in UniProt's isoform universe.**
That is the sharpest limitation on prioritising isoforms from UniProt: the annotation is a
subset of what an experimentalist working from Ensembl would consider, so a UniProt-derived
isoform list is incomplete before any expression evidence is applied.

Of the 79 that do map:

| tier | constructs |
|---|---|
| expressed | **3** — `CCR2`-2, `ADGRE5`-2, `ADGRE5`-3 |
| weak | 12 |
| not detected | 48 |
| undetectable | 16 |

`ADGRE5` is CD97, whose EGF-domain isoforms are real and well described, so the two scoring
expressed is the right answer. But the broad picture is that **GPCR isoforms are almost
invisible to bottom-up MS** — 64 of 79 have no diagnostic peptide at all — which matches
what PeptideAtlas did for `NPFFR2` earlier. For this receptor class the peptide tiers rank
almost nothing, and `isoform_priority.tsv` should not be used to deprioritise GPCR isoforms.

`NPFFR2`-2 appears here as **undetectable**, which is consistent rather than contradictory:
it is the 420-residue form the reference policy now selects, and it is a strict subsequence
of the 522-residue UniProt canonical, so no peptide can be unique to it.

**The cross-validation that would have settled this is not available from our table.**
Tedman's measured surface expression is populated for all 769 canonical constructs but for
only **3 of the 177 non-canonical** ones, so his isoform arm cannot be compared against
these tiers without going back to his supplementary data. Worth doing if isoform
prioritisation becomes load-bearing.

## Does the new reference agree with Tedman better than the old one?

Tedman's 766 GPCRs remain the only gene set with an independently chosen, experimentally
exercised reference, so they are the check on whether the policy moved toward or away from
what somebody actually built and screened. `compare_to_tedman.py` puts four references side
by side on the **526 of his genes that are library members**.

| reference | agreement with Tedman |
|---|---|
| ours today (UniProt canonical) | 96.2% (507/527) |
| **the new policy** | **96.8%** (509/526) |
| MANE Select alone | **97.9%** (509/520) |

Counting only differences bigger than two residues — a real isoform choice rather than a
reference allele or a one-residue start call — the new policy is a **modest** improvement:
**14 real disagreements with Tedman become 12.**

### What moved, and the thing that actually matters

**13 genes moved toward Tedman, 12 away.** As a headline that is a wash, but the two groups
are not alike. The genes that moved toward him are **all nine library members previously
flagged as unresolved N-terminal start-codon disagreements** — `NPFFR2`, `OR11H1`, `OR13D1`,
`OR2AK2`, `OR3A2`, `OR3A3`, `OR4F5`, `OR51V1`, `OR8G5` — plus `OR8S1` and the three
reference-allele cases `F2RL1`, `GPR42`, `OR4E2`. **Every open question from the N-terminal
group is closed by the policy, in the direction Tedman's screen-validated constructs
used.** That was the group with no MS and no Ribo-seq coverage, and it resolved on
annotation consensus instead.

The twelve that moved away are mostly genes where **MANE and APPRIS both disagree with
Tedman** (10 of 12), so they are not a policy failure so much as two annotation sources
against one experimentalist's transcript pick — and a construct being screen-validated means
it expressed, not that it is the dominant form. The largest are `GPR142` (-88),
`HRH2` (+59), `OR5AU1` (-51) and `HTR7` (-47).

### The one branch worth changing

Agreement with Tedman by policy branch:

| branch | agreement |
|---|---|
| `mane_appris_agree` | 503/514 = **97.9%** |
| `unresolved_default_mane` | 6/6 = 100% |
| **`appris_only`** (no MANE Select exists) | **0/6 = 0%** |

Every gene where the policy falls back to APPRIS because the gene has no MANE Select
disagrees with Tedman. On a 97.9% baseline, 0 of 6 is `p = 9e-11` — but the number is small
and the differences are mostly not isoform-scale: four of the six are same-length or
one-residue, leaving `HTR7` (-47) and `OR12D1` (-5) as real. So this is a flag, not a
verdict.

It still matters library-wide, because `appris_only` covers **77 genes and changes 62 of
them**, median 16 aa, including `MUC1` (1,255 -> 273). 54 of the 77 are APPRIS PRINCIPAL:1,
so these are not weak APPRIS calls — they are genes where RefSeq and Ensembl never agreed
well enough for a MANE Select to exist, which is itself a sign of messy annotation.
**Recommendation: treat the 77 as review rather than as an automatic change**, and do not
let `MUC1` through without a decision on its VNTR length.

## And against Ribo-seq (GENCODE Ribo-seq ORFs)

Ribosome profiling reports translation directly, so the GENCODE Ribo-seq ORF catalogue is
the right instinct for a start codon. `riboseq_npffr2.py` checks all three sets — Phase 1
(7,264 ORFs, GENCODE v35), Phase 2 primary (10,127) and Phase 2 comprehensive (28,359, v45).

**`NPFFR2` has no Ribo-seq ORF in any of them**, under any alias (`NPFF2`, `GPR74`,
`NPGPR`). Nor do any of the other eight library genes with an N-terminal disagreement. The
catalogue reaches **804 of the 2,369 library members**, so this is not a resource that will
cover the library.

It also cannot answer this class of question by construction. GENCODE's own scope is ORFs
"not currently represented as protein-coding CDSs in reference gene annotation", and the
Phase 2 ORF-type vocabulary is `uORF`, `uoORF`, `intORF`, `dORF`, `doORF`, `lncRNA`, `PT`,
`CDS` — none of which is an in-frame N-terminal extension of a CDS, because that would be a
change *to* a CDS rather than a new ORF. A `uoORF` at `NPFFR2` would at least have shown
ribosomes initiating upstream; there is none of those either.

### But the lookup turned up the thing that does settle it

Checking what Ensembl encodes at `ENSG00000056291` answers the question the Ribo-seq
catalogue could not:

| transcript | length | CCDS | MANE | UniProt equivalent |
|---|---|---|---|---|
| ENST00000395999 | 423 | CCDS47072.1 | | `Q9Y5X5-3` |
| **ENST00000308744** | **420** | **CCDS3552.1** | **MANE Select** | `Q9Y5X5-2` |
| ENST00000358749 | 420 | CCDS3552.1 | | `Q9Y5X5-2` |
| ENST00000940437 | 420 | CCDS3552.1 | | `Q9Y5X5-2` |
| ENST00000344413 | 304 | — | | (none) |

**No Ensembl 115 transcript encodes our 522-residue reference.** The 420-residue form is an
exact suffix of ours, so the whole difference is the first 102 residues, and that form is
what MANE Select, both CCDS entries, APPRIS and Tedman all carry. Our 522-residue sequence
is UniProt's canonical choice (`Q9Y5X5-1`) and is supported by no transcript in the
annotation our ORFs are built against.

### UniProt's own cross-references say the same thing

The decisive evidence is inside the UniProt entry we take the sequence from. `Q9Y5X5` has
four proteoforms, and every transcript-level cross-reference points at `-2` or `-3`:

| database | entries | isoform |
|---|---|---|
| CCDS | CCDS3552.1 | `Q9Y5X5-2` (420) |
| CCDS | CCDS47072.1 | `Q9Y5X5-3` (423) |
| RefSeq | NP_004876.3, NP_444264.1 | `Q9Y5X5-2` |
| RefSeq | NP_001138228.1 | `Q9Y5X5-3` |
| Ensembl | ENST00000308744, ENST00000358749, ENST00000940437 | `Q9Y5X5-2` |
| Ensembl | ENST00000395999 | `Q9Y5X5-3` |
| **`Q9Y5X5-1` (522, the Displayed sequence)** | **—** | **no CCDS, no RefSeq, no Ensembl** |

The 522-residue form rests on EMBL cDNA clones alone (`AF119815`/`AAD22047.1` and five
others, one of which UniProt flags `ALT_FRAME`). "Displayed" is UniProt's presentation
choice, not a claim that this is the dominant form.

### What this actually means for us

It is not a defect in any script. It is that **our pipeline takes "the UniProt canonical" as
the reference proteoform for every library member and never checks that the choice is backed
by a transcript.** For `NPFFR2` that silently selects a sequence no transcript database
carries.

The consequence here is narrow but real. The tag *path* does not change — DeepTMHMM calls
both forms identically (7 TM helices, no signal peptide, N-out/C-in), so both take a bare
terminal-N tag. What changes is what sits behind the tag: the ORF would be
`tag – [102 aa] – [420 aa receptor]` rather than `tag – [420 aa receptor]`. That is 306 bp of
extra synthesis, and it puts 102 residues of doubtful sequence between the tag and the
receptor's real N-terminus — in the ectodomain of a peptide-binding GPCR, which is not a
neutral place to add length.

**Nothing is ordered wrong today.** `NPFFR2` is a library member with a tag-path assignment
(`N at residue 1`, length 522) but no designed ORF, so this is prospective.

**The fix worth making is the general one**: assert that each reference proteoform maps to at
least one Ensembl or RefSeq transcript. That check would have caught `NPFFR2` without anyone
looking, and would flag the rest of this class automatically.

## Human/mouse pair consistency

APPRIS selects within each species independently, so nothing makes the human principal and
the mouse principal the same splice form. Over the 1,455 one-to-one high-confidence
Compara pairs where both species have a principal, the relative length difference has a
median of **0.3%** and a 90th percentile of **3.4%** — so the selection is usually
consistent, and the tail is informative rather than noise. **27 pairs differ by more than
10%** (`appris_pair_consistency.tsv`, the tolerance is a column so it can be re-cut).

The largest, **CD276** (human 534 aa, mouse 316 aa, 41% apart), is real biology and not a
selection error: according to PubMed, mouse B7-H3 has only the two-Ig isoform while human
has both two-Ig and four-Ig forms, from a tandem exon duplication
([Sun et al., 2011](https://doi.org/10.1371/journal.pone.0024751), PMID 21931843). A mouse
CD276 construct is not the counterpart of the human one at any isoform choice, and that is
worth knowing before it becomes a cross-species result. The remaining 26 need the same
read: real species difference, or an artefact of choosing per species.

## Cynomolgus has no APPRIS call, and none is substituted

Our macaque orthologs are ***Macaca fascicularis*** (taxon 9541; verified on the stored
accessions). APPRIS e115v50 covers *Macaca mulatta* and not *M. fascicularis*, so there is
no principal-isoform call for this cohort. Rhesus is **not** used as a stand-in: according
to PubMed the two species differ at ~9.7M SNVs with 60 nonsynonymous variants completely
differentiating them, flagged as candidates for species-specific drug response
([Higashino et al., 2012](https://doi.org/10.1186/gb-2012-13-7-r58), PMID 22747675), and
cynomolgus genomes carry rhesus introgression
([Yan et al., 2011](https://doi.org/10.1038/nbt.1992), PMID 22002653). *M. fascicularis* is
also the species sequenced specifically for drug safety assessment
([Ebeling et al., 2011](https://doi.org/10.1101/gr.123117.111), PMID 21862625), which is
why the library is on it.

`cyno_no_appris.tsv` lists all 1,921 stored cyno orthologs with the human gene each was
found from; **1,432 pair to a current library member**. Those keep their existing
UniProt-selected sequence, with no isoform-level reference call available.

## A better Ribo-seq resource, and what it could actually reach

The GENCODE catalogue was the wrong *kind* of resource rather than a bad one: it catalogues
ORFs **not** represented as CDSs, so it will never hold two candidate CDSs to compare. What
this question needs is **read coverage over annotated CDS**, and for the start-codon cases,
**initiating-ribosome data**.

According to PubMed, the resource for that is **RiboSeq.Org**
([Tierney et al., 2025](https://doi.org/10.1093/nar/gkae1020), PMID 39540432): GWIPS-viz for
genome-wide coverage, Trips-Viz for transcript-level profiles, and the Ribo-seq Data Portal
holding **14,840 preprocessed libraries from 969 studies across 96 species**, with
standardised metadata and programmatic access. All three are live. For isoform-level
quantification specifically there is **ORQAS**
([Reixachs-Solé et al., 2020](https://doi.org/10.1038/s41467-020-15634-w), PMID 32286305),
which quantifies translation per transcript isoform from ribosome-protected fragments.

### But first, how much of our problem is reachable at all

| where the difference falls | genes | median unique region | still unresolved |
|---|---|---|---|
| internal | 94 | 18 aa | 67 |
| **n_terminal** | **13** | 23 aa | 11 |
| substitution_only | 3 | 3 aa | 3 |
| c_terminal | 1 | 0 aa | 1 |

- **13 N-terminal differences** — alternative start codons. Initiating-ribosome profiling
  (harringtonine/lactimidomycin, in GWIPS-viz's initiating tracks) reads start codons out
  directly, which is exactly this question. And **none of the 13 is an olfactory receptor**:
  they are `SLC1A2`, `GYPC`, `TYRO3`, `PMEPA1`, `SCNN1D`, `PVRIG`, `EVI2A`, `KCNN3`,
  `CRIPTO`, `LY6H`, `ANO7`, `CD200R1L`, `LRRC37B` — ordinary, well-expressed genes that
  public datasets will cover. This is the highest-value use.
- **37 genes differ by 30 aa or more**, likely a distinct exon, so isoform-level Ribo-seq
  (ORQAS, or Trips-Viz coverage across the junction) can speak to them.
- **33 differ by under 10 aa.** No translation data resolves those. The MANE default stands
  and should be left alone.

So a realistic ceiling is **roughly 40–50 of the 111**, not all of them, and the 13
N-terminal cases are where to start.

### One finding that bears on the whole exercise

ORQAS reports that **53% of expressed human genes have more than one translated isoform**,
with translation detectable for 40–50% of expressed isoforms. That cuts against the premise
of choosing a single reference: for many genes both the MANE and the APPRIS form are
genuinely made. "Which isoform is correct" is often the wrong question, and the right one is
which isoform we want to display — a design decision, not an annotation lookup.

### Correction: isoforms were already topology-filtered in s1e3

An earlier version of this note said no TM filter had ever been applied to isoforms. That is
wrong, and the prior work is in
[`library/l1/surface-filters`](../../../library/l1/surface-filters).
`build_library.py` classifies every isoform's `topology_status` from its DeepTMHMM label and
TM count, and computes the TM delta against the canonical. Its last build recorded:

| isoform topology_status | isoforms |
|---|---|
| surface | 2,679 |
| secreted | 639 |
| intracellular | 182 |

with recommendations **include 2,203, drop 821, flag 476**. The 821 dropped are the same
class this analysis rediscovered as "loses every TM helix".

`isoform_expression_scoping.md` in that folder had also already scoped most of what is built
here — APPRIS principal, MANE, TSL, ECD identity, and Tedman surface expression for GPCRs —
including the observation that *"~half of screened isoforms were Ensembl-only (no UniProt
isoform) ... source isoforms from Ensembl/GENCODE + run per-transcript DeepTMHMM rather than
relying on UniProt isoforms only."* That is precisely the gap this work closed, so it is an
execution of an existing plan rather than a new one.

What is genuinely added here: the forms are sourced from Ensembl as well as UniProt and
deduplicated by sequence; evidence is observed peptides and OpenProt translation rather than
the GTEx TPM proxy the earlier doc proposed; and the TM comparison is against the newly
selected MANE+APPRIS reference rather than the old UniProt canonical.

### The 9 genes with no MANE and no APPRIS principal

They now take the **reviewed UniProt canonical**, so every one of the 2,369 members has a
reference. They are not a random nine, and the evidence splits them cleanly:

| gene | MS coverage | peptides | verdict |
|---|---|---|---|
| `KIR2DL5A` | 60.5% | 21 | **keep** |
| `KIR3DS1` | 56.8% | 18 | **keep** |
| `KIR2DL2` | 54.0% | 23 | **keep** |
| `KIR2DS2` | 52.6% | 21 | **keep** |
| `KIR2DS5` | 26.0% | 10 | **keep** |
| `OR9G9` | 3.3% | 1 | candidate to drop |
| `OR8U9` | 2.3% | 1 | candidate to drop |
| `OR8U8` | 0% | 0 | candidate to drop |
| `TAS2R45` | 0% | 0 | candidate to drop |

**The five KIR genes are well-observed proteins** — half their sequence covered by peptides.
They lack MANE and APPRIS only because the KIR locus sits on alternate haplotypes rather
than the primary assembly, which is an annotation artefact, not a statement about the
protein. Keep them.

The four olfactory and taste receptors have essentially no evidence, and OpenProt adds
nothing (`TAS2R45` TE=0, the rest absent). Dropping them is a scope decision rather than an
evidence one: they are real genes that nothing has ever detected outside their native
epithelium, which is true of most of the library's olfactory arm. **Tedman settles the four.** Three of them are in his screen as canonical constructs and all
three displayed **above his all-construct median of 869**: `OR8U9` 1052.0, `OR9G9` 1048.8,
`OR8U8` 1044.5. (`TAS2R45` is a taste receptor and not in his GPCR panel.) Across his 393
olfactory receptors the median is 823, so these three are at the better end of that class.

That is a direct surface readout and it overrides the absence of mass-spectrometry evidence:
these proteins are not natively made outside olfactory epithelium, which is why no tissue
proteomics sees them, but **as constructs they display**. For a synthetic library the second
is the question that matters. **Keep all nine.**

### Why MANE is the default for the 82, and what OpenProt changed

The default rests on two things: MANE is a RefSeq/Ensembl agreement and therefore
transcript-backed by construction, and APPRIS labels the MANE transcript **ALTERNATIVE** in
all 111 — but 74 of the 111 reach only PRINCIPAL:4 or worse, APPRIS's own weakest tiebreak.
Choosing the weaker-supported candidate on no evidence would be the harder position to
defend.

That said, **OpenProt now calls 30 of the 82**, so for those the default is replaced by
evidence rather than retained — 22 going to MANE and 8 to APPRIS. 52 still rest on the
default alone.

### Do peptides and OpenProt agree on the 9 that peptides resolved?

**Seven agree, one has no call, and one disagrees.**

| gene | peptides | OpenProt TE (MANE vs APPRIS) | |
|---|---|---|---|
| `SPPL2B` | MANE | 22 vs 8 | agree |
| `ADAM15` | MANE | 13 vs 2 | agree |
| `DCBLD1` | MANE | 12 vs 1 | agree |
| `LRP8` | MANE | 11 vs 1 | agree |
| `SLC12A6` | MANE | 10 vs 7 | agree |
| `PTPRC` | MANE | 6 vs 0 | agree |
| `KCNN3` | MANE | 1 vs 0 (called on MS) | agree |
| `ATP2B2` | MANE | 1 vs 0 | tied, no call |
| **`CD44`** | **MANE** (49 peptides) | **5 vs 17** | **disagrees** |

**`CD44` and `ATP2B2` now carry both lanes.** The peptide-favoured form stays the reference
and the other is recorded as a `second_lane` to be built as an isoform, rather than one being
discarded: `CD44` 742 + 361, `ATP2B2` 1243 + 1198. `CD44` is the interesting one and probably
not an error in either method: the long
variant-exon form carries 49 diagnostic peptides and the higher MS score (229 vs 158), while
the short standard form has more ribosome-profiling studies behind it. Both `CD44s` and
`CD44v` are genuinely expressed in different tissues. This is the case ORQAS describes —
53% of expressed human genes translate more than one isoform — and it is a design choice
rather than a fact to look up.

## Where the human library stands

### One reference per gene — 2,369 members

| | genes |
|---|---|
| MANE and APPRIS agree | **2,172** |
| no observation either way, defaults to MANE | 82 |
| APPRIS principal only (gene has no MANE) | 77 |
| weak evidence, needs a look | 20 |
| resolved on peptides | 9 |
| no reference available | 9 |

**216 references change** against what we build from today, **44 of those were encoded by no
Ensembl transcript at all**, and **29 need a person**. MS coverage of the chosen references
is a median **58.8%**, with 89% covered at all and 258 (11%) never observed; 1,410 carry
OpenProt translation evidence.

### Alternate forms — 13,687 across UniProt and Ensembl

| evidence | forms |
|---|---|
| protein evidence (diagnostic peptides, deep) | 68 |
| weak protein evidence | 245 |
| translation evidence only (OpenProt TE) | 678 |
| **any positive evidence** | **991** |
| predicted only, TRIFID > 0.5 | 1,119 |
| APPRIS MINOR — deprioritise | 10,736 |

10,790 of the 13,687 are Ensembl-only, which is why the universe had to be widened.

### Membrane filter

**372 forms lose every TM helix and are excluded** (`FLT1`-2 = sFlt-1 among them); **295
lose some and are kept**, as a designer judgement. **11,433 forms have no topology call**,
because DeepTMHMM has only ever run over the UniProt cohort.

Running it over all of them is not worth the compute, and most of the question is already
settled: **s1e3 classifies every isoform in the UniProt cohort by `topology_status`** and
already drops the secreted and intracellular ones. Re-running those would re-decide what is
decided. What s1e3 never saw are the Ensembl-sourced forms — the gap its own
`isoform_expression_scoping.md` identified.

`shortlist_for_topology.py` scopes the run in three steps:

1. **Which forms could change a decision?** Those with observed evidence of expression, or
   TRIFID above 0.5 — **2,110 of 13,687**.
2. **Which of those already have a call?** 813 do. Of those, **752 are decidable right now**,
   because their reference has a call too, so the TM delta can be computed from data in
   hand. No compute needed.
3. **What is genuinely missing?** 1,297 alternates, **plus 230 references**.

That second number is the correction that matters. A TM count means nothing on its own —
it is only a comparison against the reference's count. The first version of this shortlist
queued forms whose *reference* had never been run, which would have produced 230 genes'
worth of isoform calls that could not be interpreted, and it also queued 61 forms that
already had calls.

| | |
|---|---|
| **run size** | **1,527 sequences over 719 genes** |
| alternates qualifying on evidence | 331 |
| on TRIFID alone | 883 |
| on both | 83 |
| references, needed for the comparison | 230 |
| total residues | 960k |

`topology_shortlist.fasta` and `topology_shortlist.tsv` are written and ready, and the
manifest marks which rows are references rather than candidates.

**Results must be read through s1e3's rule, not a new one** — `topology_status` from the
DeepTMHMM label (TM / SP+TM / BETA = surface; SP = secreted; GLOB = intracellular), so the
new forms are judged the same way the existing 2,238 were.

### The run, and what it decided

DeepTMHMM v1 ran on Modal over the 1,527 sequences: **17 shards, 0 failures**. The package
is license-gated and lives on a Modal volume rather than in Git; `modal/deeptmhmm_v1_app.py`
holds the app, and the canary is worth keeping — the first attempt failed on
`cannot mount volume on non-empty path: "/pkg"`, which costs nothing to find on 2 sequences
and would have wasted the whole sweep.

Judged by **s1e3's rule**, not a new one — DeepTMHMM label decides, TM / SP+TM / BETA are
surface, SP is secreted, GLOB is intracellular:

| verdict | forms |
|---|---|
| keep | 1,124 |
| **exclude — lost surface** | **109** |
| keep, lost some TM helices | 61 |
| undecidable, reference still has no call | 3 |

Of the 1,297 alternate forms: 1,137 surface, 116 secreted, 44 intracellular.

**104 of the 109 exclusions are Ensembl-only forms** — they exist only because the universe
was widened past UniProt, and they would otherwise have entered the library untested. Only
**6 of the 109 had any protein evidence**, so this is mostly removing predicted forms that
would have ranked on TRIFID alone.

The net effect of widening to Ensembl is **1,188 new forms that pass the membrane filter**
and 109 caught before they could be ordered.

### The current isoform cohort — 2,238 forms

| | forms |
|---|---|
| with any positive evidence | 507 |
| **flagged to exclude** (lose every TM helix) | **372** |
| expressed forms missing from the cohort | 25 |

### What is left to decide

1. **29 references need a person** — 20 weak-evidence calls plus 9 genes with no reference
   (the KIR cluster, three olfactory receptors, `TAS2R45`).
2. **`MUC1` 1,255 → 273** is a VNTR-length decision, not an annotation error. Do not let it
   through unreviewed.
3. **The 77 `appris_only` genes** change 62 references and agree with Tedman 0 for 6. Treat
   as review, not an automatic change.
4. ~~Run DeepTMHMM v1 over the 1,527 shortlisted sequences~~ — **done**: 109 excluded, 1,124
   kept, 61 kept having lost some helices. Publishing the calls to D1 is the remaining step,
   and needs a schema decision for the Ensembl-only forms, which have no UniProt accession.
5. **Eight N-terminal start-codon cases** remain open after OpenProt closed five; they want
   GWIPS-viz initiating-ribosome data.

## The selection policy, implemented

`select_reference.py` chooses one reference proteoform per library member: **MANE Select and
the APPRIS principal where they agree; observed evidence where they do not.** Agreement is
tested on the *protein sequence*, not the transcript ID — a gene commonly has several
transcripts encoding one protein, and calling that a disagreement would be wrong.

| outcome | genes |
|---|---|
| **MANE and APPRIS agree** → take it | **2,172** (91.7%) |
| they disagree → adjudicate on evidence | 111 |
| APPRIS principal only (gene has no MANE Select) | 77 |
| neither | 9 |

### What the APPRIS class actually means, and whether we are comparing like with like

The class is **the rule APPRIS used to pick**, not a quality score for the transcript.
PRINCIPAL:1 means the core modules — structure, functional residues, cross-species
conservation — chose a clear winner. For roughly a quarter of human genes they cannot, and a
tiebreak cascade runs instead: **P2** best TRIFID score by a wide margin, **P3** most
supporting proteomics evidence, **P4** best TRIFID score without a wide margin, **P5**
longest. APPRIS regards P1, P2 and P3 as reliable; P4 and P5 are the weak end.

Our 111 disagreements: **74 P4, 23 P2, 10 P3, 3 P:M, 1 P1.** So two thirds sit at APPRIS's
weakest tiebreak, but a third do not, and the P2/P3 cases are ones APPRIS itself calls
reliable.

**And yes, this is top transcript against top transcript.** The check that matters is what
APPRIS calls the MANE transcript, and in **all 111 it is ALTERNATIVE** — never a lower-ranked
principal, never a tie my own tiebreak had to break (zero cases). APPRIS has actively demoted
the MANE model to a candidate it did not choose:

| APPRIS's label for the MANE transcript | genes | meaning |
|---|---|---|
| ALTERNATIVE:2 | 64 | conserved in **fewer than three** tested non-primate species |
| ALTERNATIVE:1 | 45 | conserved in **at least three** non-primate species |
| ALTERNATIVE:M | 2 | the gene was manually curated afterwards |

That is a real argument against MANE on 64 of them, on conservation grounds. It does not
override MANE here — MANE is a RefSeq/Ensembl agreement and transcript-backed by
construction, and conservation is one axis — but it is why these cases go to evidence rather
than to a rule.

### Adjudicating the 111

For each, the script computes the residues that actually **distinguish** the two candidates
(`difflib` opcodes), then asks whether anything was ever observed there. A supporting
peptide must lie inside a region unique to one candidate, be **absent** from the other, and
be **proteotypic** — found in no other gene's reviewed proteoforms. Without that last test a
sequence shared with an unrelated protein would read as evidence for this gene's splice
form. All substrings of 7–30 aa are tested, so tryptic, semi-tryptic and missed-cleavage
peptides are all caught without assuming a protease.

A hit also has to be **deep**, not just present. PeptideAtlas reports `n_observations` and
`n_samples` per peptide, and the build spans peptides seen once and peptides seen in 1,852
samples. Counting those equally is what made the first pass over-trust single-peptide hits,
so a verdict now needs 3+ proteotypic peptides **and** a best peptide seen in **at least 10
distinct samples**. The data chooses that cut: the cases clearing both run 17–137 samples,
every case that fails runs 1–8, and the gap between is empty. It demoted `PTPRG` (3 peptides,
best seen once) and `FCRL6` (3 peptides, 3 samples) out of the resolved set.

| verdict | genes |
|---|---|
| **resolved for MANE** (≥3 proteotypic peptides, best in ≥10 samples) | 9 |
| weak evidence for MANE — review | 13 |
| weak evidence for APPRIS — review | 7 |

The 18 weak cases split 11 MANE / 7 APPRIS, and the direction often follows the geometry
again: four of the seven APPRIS-leaning genes (`ANO7`, `APLP2`, `SLC22A17`, `SPINT1`) have
**no MANE-unique residues**, so APPRIS was the only side the test could support. Peptide
quality varies more than the counts suggest — `NRG2` (`YIFFLEPTEQPLVFK`, 15 aa), `P2RX2`
(17 aa), `PCDH19` and `SLC24A1` (16 aa) are substantial, while `ART3` (`LEDHSEK`),
`TNFRSF18` (`SQCMWPR`) and `ASIC2` (`DMMLYCK`) rest on a single 7-mer. `ASIC2` deserves the
closest look: it is the only PRINCIPAL:1 disagreement, it has large unique regions on both
sides (141 and 92 aa), and one 7-mer decided it.
| no observation either way → default to MANE | 82 |

**Nothing resolved strongly for APPRIS — but that is mostly the test's geometry, not
evidence against APPRIS.** In **70 of the 111 the APPRIS form has no unique residues at
all**: it is a strict subsequence of the MANE form, so no peptide could ever support it. The
test only had an APPRIS target on 34 genes, against 70 for MANE, and the MANE-unique regions
total 3,381 residues against 1,041 — a 3.2x asymmetry in opportunity before any data is
looked at.

Conditioned on having a target, the two sides are much closer:

| | genes where a supporting peptide was geometrically possible (>=7 unique aa) | genes with a hit |
|---|---|---|
| MANE | 70 | 22 (**31%**) |
| APPRIS | 34 | 7 (**21%**) |

(Hit counts here are any supporting peptide, before the depth and count thresholds.)

So read `resolved_mane` as *"the MANE-unique region is translated"*, which is a sound reason
to build from MANE — it is the form we can confirm exists. It is **not** evidence that the
APPRIS form is absent; where one is contained in the other, both can be expressed and the
peptide data cannot separate them.

### Why the 11 resolved, and why the 82 did not

The nine are not genes where MANE is specially right — they are genes where the test had a
large target and deep coverage. Their median MANE-unique region is **66 aa against 11 aa across
all 111**, and their best supporting peptide was seen in a median of **44 samples** against
**2** for the weak cases, and they are the familiar large alternatively-spliced cassettes: `CD44`
(the variant exons), `PTPRC` (the CD45 RA/RB/RO exons), `KCNN3`, `DCBLD1`, `LRP8`, `ADAM15`,
`SPPL2B`. Where a cassette is big and genuinely expressed, mass spectrometry sees it.

The 82 failed for want of a target, not for want of data. **No protein among them had zero
peptides** — the median was 5,183 in its batch pool. What they lack is anywhere to look:
median MANE-unique region **8 aa**, median APPRIS-unique **0 aa**, and 32 of the 82 have
under 10 unique residues on *both* sides. An 8-residue region contains at most two distinct
7-mers, so the test was never going to fire. These are genes distinguished by a handful of
residues, and proteomics at this depth cannot speak to them.

The default for the 82 unresolved is **MANE**, because MANE is a RefSeq/Ensembl agreement
and therefore transcript-backed by construction, while two thirds of these are APPRIS
PRINCIPAL:4 — its weakest tiebreak. The `mane_appris_label` column carries APPRIS's
counter-argument per gene so the default can be revisited where it is ALTERNATIVE:2.

### Ribo-seq, tested on the sequences rather than the labels

The catalogue ships translated ORF sequences (`.faa`), not just types, so the real test is
whether any Ribo-seq ORF translation overlaps a discriminating region. All 35,169 ORF
translations across the three sets were tested against every 7–30 aa window of every
discriminating region.

**51 of the 111 genes have a Ribo-seq ORF; zero have one that overlaps a discriminating
region.** The lookup is not vacuous — translations loaded for all 51 — but there are only
**128 ORF translations across those 51 genes, a median of 2 each**, and the median
catalogued ORF is **21 aa**. Their types are what the catalogue is for: 40 `uORF`, 11
`uoORF`, 8 `intORF`, 1 `dORF`.

**It is not true, though, that the catalogue can never contain annotated CDS sequence**, and
an earlier version of this note said so. Testing every Phase 2 ORF against its own gene's
Ensembl translations, 569 of 28,359 are substrings of one. Almost all of that is an artifact
of length — the median matching ORF is **1 residue**, and shuffling the matched sequences
leaves 497 of 569 still matching, so they are chance hits of 1–4-mers. At a length where a
match means something, it is **4 ORFs in 28,359** (`TGFBR3L` 56 aa, `A1BG` 46 aa, `FSTL3`
39 aa, `ADAR` 24 aa), and five Phase 1 ORFs are identical to a whole annotated protein. So
the right statement is that in-frame CDS sequence is **possible but vanishingly rare
(~0.01%)**, not impossible — and the zero overlap on our 111 is a measured result, not a
foregone one.

**Why Phase 2 has fewer full-CDS matches than Phase 1, despite being four times larger.**
Phase 2 is **not a superset**: only **3,020 of Phase 1's 7,264 ORF ids survive into it**. It
is a re-derivation on GENCODE v45 with translation-signature scoring, not an extension of the
v35 set. And the catalogue is defined as ORFs *not* represented as protein-coding CDSs, so an
ORF adopted into the reference annotation leaves it — `ERCC6` and `CYP27B1` have Phase 1 ORFs
and **no Phase 2 entry at all**.

What the five Phase 1 full-protein matches actually matched is worth knowing before reading
too much into them:

| gene | ORF | all Ensembl 115 translation lengths for the gene |
|---|---|---|
| `TASOR2` | 43 aa | 8, 43, 71, 117, 304, 1255, … 2,691 |
| `A1BG` | 46 aa | 46, 305, 495, 517, 549 |
| `ERCC6` | 80 aa | 77, 80, 179, 224, 1061, … 1,493 |
| `CYP27B1` | 198 aa | 198, 347, 350, 508, 535 |
| **`ERVH48-1`** | **160 aa** | **160 — the gene's only translation** |

Four of the five match a **short alternative transcript** of a much larger gene, not the
gene's main CDS. Only `ERVH48-1` (suppressyn, an endogenous-retrovirus-derived protein) is a
case where the Ribo-seq ORF is the complete and only annotated CDS of its gene — which is
exactly the boundary case GENCODE describes, an ORF that has made it into the annotation.

The one Phase 2 ORF typed `CDS` is `PEG10`, a retrotransposon-derived gene known for -1
ribosomal frameshifting, so an ORF inside the CDS region in another frame is expected there
rather than anomalous.

**So the honest accounting is that Ribo-seq changed no decision, before or after the fix.**
It could have helped on 27 genes — unresolved by peptides and with at least 7 unique
residues somewhere to look — and helped on none of them. What the fix bought is the
difference between asserting that from the type vocabulary and testing it against all 35,169
ORF translations: the negative is now measured rather than argued. It also settles that this
resource will not help this class of question again, which is worth knowing before anyone
reaches for it a second time.

### A cross-check worth knowing about

PRINCIPAL:3 means APPRIS chose its form because it had **the most supporting proteomics
evidence** — the same axis this analysis uses. On those 10 genes our independent PeptideAtlas
test leans **MANE on 4** (`ALCAM`, `ART3`, `CLSTN1`, `SLC12A6`), **APPRIS on 2** (`SLC12A5`,
`SPINT1`) and finds nothing on 4. So APPRIS's proteomics argument and PeptideAtlas Human
2026-01 do not agree, and a P3 label should not be read as settling the question.

### What the policy changes

**216 of 2,369 references change**, median difference 16 aa, our current sequence longer in
137 of them — the signature of UniProt's length tiebreak.

Among those, **44 genes currently build from a sequence encoded by no Ensembl transcript at
all**, and the policy changes every one. `NPFFR2` was not special; it was one of 44. The
largest is `MUC1` (1,255 → 273), which is the VNTR tandem-repeat region and a genuine
design decision rather than an error — **look at that one by hand before accepting it.**

**No staged order gene is affected.** `HEG1`, `FOLH1`, `VSIG10`, `TFRC`, `AXL`, `ERBB2`,
`EGFR`, `MET`, `LDLR`, `CD19`, `CD22`, `NECTIN4` and `TMEM123` all keep their current
sequence.

**27 genes need a person**: the 18 weak-evidence cases plus the 9 with no reference at all
(the KIR cluster, three olfactory receptors and `TAS2R45`).

### A correction worth recording

The transcript-support count was **154** on the first run and is **44**. The first version
tested our sequence only against the transcripts APPRIS lists, which is a subset of the
gene's translations, so genes whose sequence *is* encoded — `SUSD6` among them — were
reported as unsupported. The check now scans every translation in the Ensembl peptide FASTA
by gene, and distinguishes "no translations found" (11 genes, unknown) from "translations
found, ours is not among them" (44, unsupported). Spot-checked against an independent
exhaustive scan for `MUC1`, `KCNN2`, `NPFFR2`, `SUSD6`, `EGFR` and `TFRC`.

## Running it

```bash
uv run --locked python fetch_inputs.py --cache-dir ~/.cache/deliverome/appris
uv run --locked python select_principal.py \
    --surfaceome-repo ~/Git/accessible-surfaceome \
    --cache-dir ~/.cache/deliverome/appris
uv run --locked python select_reference.py --surfaceome-repo ~/Git/accessible-surfaceome
```

`select_reference.py --offline` reruns the policy from the cached PeptideAtlas responses
without touching the network. Peptides are kept with their `n_observations`, `n_samples` and
`organism_full_name`; the species field is checked rather than assumed, though build 607 is
human-only. PeptideAtlas is queried through
`biosequence_name_constraint`, which accepts a semicolon-separated list, so the 111
accessions take six batched requests spaced three seconds apart rather than 111 — the CGI
asks clients not to hammer it.

`fetch_inputs.py --publish` uploads the two APPRIS tables and their provenance to
`external://appris/e115v50/`. The Ensembl 115 peptide FASTAs and the UniProt reviewed
streams stay in the local cache; they are large and already versioned upstream.

| input | version | why pinned |
|---|---|---|
| APPRIS principal tables | dataset **e115v50**, snapshot **2026_06.v50**; files built 2026-08-31 (human) and 2026-08-29 (mouse) | the release this selection is defined by. Fetched from the dated `releases/2026_06.v50/` directory, not `current_release/`, which is a moving symlink; the two were byte-identical when checked (sha256 in `provenance.json`) |
| Ensembl peptide FASTA | release 115 | APPRIS e115v50 indexes Ensembl 115 transcripts; the live REST service serves 116, where a transcript can be retired or re-versioned, so a sequence fetched from it is silently the wrong one |
| UniProt reviewed proteomes (+isoforms) | 2026_03 | the comparison set, and the source for the 587 library members absent from the DeepTMHMM cohort |
| DeepTMHMM cohort FASTAs, ortholog metadata | `accessible-surfaceome` | what the pipeline actually holds and builds from |

## Open

- The stored human cohort covers only **1,775** of 2,369 library members; the other 587 are
  compared against a fresh UniProt 2026_03 fetch, marked in `sequence_source`. The cohort
  needs regenerating against the current library before any of its precomputed topology or
  disorder tracks can be used library-wide.
- The 17 human and 15 mouse `n_terminal` disagreements are unreviewed. They are the only
  bucket that can move a post-SP tag placement, so they are the ones worth the time.
- Mouse substitution-only cases cannot be resolved against UniProt annotation. Comparing the
  Ensembl GRCm39 translation against the specific strain each UniProt entry was derived from
  would settle them; that strain is in the UniProt entry and is not read here.
- The 21 olfactory-receptor N-terminal start-codon disagreements (8 in the library) are
  unresolved. Neither MS nor Ribo-seq reaches them, so an N-terminomics dataset is the route
  worth trying next. `NPFFR2`, the 22nd, is resolved against us and needs fixing.
- Whether the other eight library members in that group are, like `NPFFR2`, unsupported by
  any Ensembl transcript has not been checked; `riboseq_npffr2.py` does that check for one
  gene and should be run across all nine.
- The 22 N-terminal start-codon disagreements (9 in the library) are unresolved. Tedman's
  shorter form is screen-validated, which is evidence our canonical is too long at the N
  terminus for those genes; a terminal-N tag there currently sits on 4-47 residues that may
  not be translated.
- The three-way comparison only covers Tedman's 766 GPCRs. There is no validated reference
  choice for the rest of the library, so the 91.6% APPRIS/ours agreement cannot be assumed
  to hold outside GPCRs.
- The pair-consistency check uses length only. Two isoforms of equal length that differ by a
  swapped cassette exon pass it.

## Decisions

| date | decision |
|---|---|
| 2026-10-08 | Import APPRIS P1 and select human isoforms by it (Becca). |
| 2026-10-08 | Selection takes the best APPRIS class present per gene, not PRINCIPAL:1 only — requiring P1 would drop 36% of the library. |
| 2026-10-08 | APPRIS is carried as an annotation; UniProt stays the build sequence, because the topology, cleavage and disorder tracks were computed on it. |
| 2026-10-08 | APPRIS is not applied independently to mouse orthologs; per-species selection is the failure mode, and 27 pairs fail the length check. |
| 2026-10-08 | Reference choice cross-checked against Tedman's 766 screen-validated GPCRs: all three agree on 92-96%, and where they differ it is N-terminal, so the disagreements land exactly on the tag site. |
| 2026-10-08 | NPFFR2's start codon tested against PeptideAtlas Human 2026-01; zero peptides in the disputed extension but p = 0.18-0.21, so the question stays open and the reference is unchanged. |
| 2026-10-08 | NPFFR2 reclassified from open question to correction: no Ensembl 115 transcript encodes our 522-residue reference, and MANE Select, both CCDS entries, APPRIS and Tedman all carry the 420-residue form. |
| 2026-10-08 | Reference policy: MANE + APPRIS where they agree (2,172 of 2,369); PeptideAtlas discriminating-region evidence where they do not; MANE as the default for the 82 with no observation (Becca). |
| 2026-10-08 | The policy closes all nine N-terminal start-codon questions in the direction Tedman's screen-validated constructs used, and cuts real isoform disagreements with him from 14 to 12. |
| 2026-10-08 | The `appris_only` branch (77 genes, no MANE Select) is flagged for review rather than applied automatically: it agrees with Tedman 0/6, and it changes 62 references including MUC1. |
| 2026-10-08 | "No strong APPRIS support" is a geometry artifact and is recorded as one: in 70 of 111 the APPRIS form is contained in the MANE form, so no peptide could support it. Hit rates conditioned on having a target are 31% MANE vs 21% APPRIS. |
| 2026-10-08 | The 9 genes with no MANE and no APPRIS principal take the reviewed UniProt canonical, so all 2,369 members have a reference. The five KIR genes are well-observed and stay; the four olfactory/taste receptors have no evidence and are a scope decision. |
| 2026-10-08 | APPRIS TRIFID ranks the 11,631 forms no experiment can reach; it recovers the observed ordering on the forms that can be measured (AUC 0.731), which is what justifies using it where they cannot. |
| 2026-10-08 | The membrane-loss filter uses DeepTMHMM, not APPRIS THUMP, which scores a 7TM GPCR as 0 helices: 372 forms lose every TM helix and are excluded, while 295 that lose only some are kept as a designer judgement. |
| 2026-10-08 | OpenProt is a second adjudicator for the reference choice, not only for isoforms: it calls 30 of the 82 disagreements the peptide test could not, agrees with it 18-3 where both speak, and closes 5 of the 13 N-terminal start-codon cases. |
| 2026-10-08 | ~250 references (1 in 10) have neither peptide coverage nor OpenProt translation evidence; they are orphan/olfactory/tissue-restricted and must not be penalised for being unrankable. |
| 2026-10-08 | The form universe is UniProt + Ensembl, not UniProt alone: Ensembl contributes 10,880 of 13,687 alternate forms. |
| 2026-10-08 | "Not detected" is qualified by the gene's own MS coverage; median coverage is 58.8% so most undetected forms are real negatives, and only 294 forms are set aside as too thinly covered to judge. |
| 2026-10-08 | OpenProt TE scores supply translation evidence for 678 forms the peptide tiers could not call, 350 of them undetectable by MS in principle. neXtProt is end-of-life and not used. |
| 2026-10-08 | Alternate isoforms are prioritised on observed diagnostic peptides rather than on UniProt annotating them; 142 score expressed, 24 of those are additions to the cohort, and 790 cannot be evidenced by bottom-up MS at all. |
| 2026-10-08 | 48% of Tedman's non-canonical constructs match no UniProt proteoform, so a UniProt-derived isoform list is incomplete before any evidence is applied; and GPCR isoforms are near-invisible to MS, so the peptide tiers must not be used to deprioritise them. |
| 2026-10-08 | "Not detected" is a ranking signal, not a deletion list -- PeptideAtlas depth is uneven and tissue-restricted proteins look silent. |
| 2026-10-08 | If Ribo-seq is revisited, use RiboSeq.Org (GWIPS-viz initiating tracks) on the 13 N-terminal cases, not the GENCODE ORF catalogue; ~33 of the 111 differ by under 10 aa and are out of reach of any translation data. |
| 2026-10-08 | Ribo-seq is tested on ORF translations, not ORF-type labels; it overlaps no discriminating region in any of the 111, so it contributes nothing here. |
| 2026-10-08 | Evidence has to clear 3 proteotypic peptides **and** a best peptide seen in 10+ samples; peptide count alone passed two one-sample cases. |
| 2026-10-08 | Evidence has to clear 3 proteotypic peptides to move a reference. Below that the gene is flagged for review, not flipped. |
| 2026-10-08 | Add a transcript-support assertion to reference selection: every reference proteoform must map to >=1 Ensembl or RefSeq transcript. NPFFR2 is the case that motivates it. |
| 2026-10-08 | No APPRIS proxy for cynomolgus. *M. mulatta* is not substituted for *M. fascicularis*. |
