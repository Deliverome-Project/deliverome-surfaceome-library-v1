# L1 references that differ from the UniProt canonical: evidence review, 2026-11-22

**Question.** For each L1 gene, is the protein we build a better choice than the UniProt canonical
where the two differ? Becca asked this in her comments on the brief: "have you checked all of the
proteins with major differences vs UniProt?" and "how sure are we that these differences are
correct and the better protein form?".

**Scope.** The review covers every gene where the reference record
(`library-design/reference-selection/2026-11-22-round2-additions/`) marks the reference as
different from the canonical. That is 214 genes. 212 of them are in the library, and 205 of
those belong to the 206 genes the brief cites. OR8J2 and SIGLEC16 have no designed CDS.

Two terms used throughout:

- **Ours** is the designed CDS translated.
- **Canonical** is the UniProt 2026_03 canonical isoform of the same gene.

Here "recommendation" means the review's proposal. The decisions taken on it are in the
decision log at the end.

## Method

| Evidence | How it is measured |
|---|---|
| Segments unique to each form | `difflib` regions of at least 3 residues, as in the reference ladder |
| Topology | DeepTMHMM v1 by exact sequence (`topology.for_sequence`): signal peptide, TM count and N-terminal side. "Differs" means one of those three differs. |
| Pfam | From the proteoform audit (`l1_proteoform_audit.tsv`) |
| PeptideAtlas build 607 | Peptides observed for **every** reviewed accession of the gene: the canonical plus all isoforms. A peptide **distinguishes** a form when that form contains it and the other does not, and when it is proteotypic (found in no other gene's reviewed proteoform). |
| GENCODE Ribo-seq ORFs | Translations overlapping each form's unique segments |
| MANE Select | The Ensembl 115 GTF `MANE_Select` tag, compared with ours and with the canonical |

Two choices in the PeptideAtlas test differ from the reference ladder's own peptide test, and each
changes the answer:

1. **All isoform accessions are queried.** PeptideAtlas lists a peptide under the protein it maps
   to, so a peptide unique to our isoform never appears under the canonical accession. A query on
   the canonical alone found no evidence for ours in any gene. Querying all isoforms finds evidence
   for ours in 99 genes.
2. **Peptides that span a deletion junction count.** The ladder's test looks only inside residues
   unique to one form, so a form that is a pure deletion of the other can never score. This
   matters for PTPRM, below.

"Strong" uses the ladder's own numbers (`select_reference.py`):

- a peptide is "seen in ≥ 10 samples" (`MIN_SAMPLES`), which separates real observations from
  single-sample hits;
- the ladder overrules MANE only when one side has at least three such peptides and the other
  side has none (`STRONG_EVIDENCE`).

The PeptideAtlas cache is extended and never overwritten. Each new key is a superset whose
earlier entries are asserted unchanged:

| Key (external bucket) | Adds |
|---|---|
| `peptideatlas/build-607/queries-2026-11-22/` | The canonical accessions of the 206 (113 → 289 accessions) |
| `peptideatlas/build-607/queries-2026-11-22-isoforms/` | Every other reviewed isoform of those genes (→ 793) |
| `peptideatlas/build-607/queries-2026-11-22-round2/` | The genes added in round 2 (→ 804); the review reads this one |

Code: [proteoform_review.py](proteoform_review.py). Tables:

- [proteoform_review_2026-11-22.tsv](proteoform_review_2026-11-22.tsv): one row per gene, with
  evidence, recommendation and basis;
- [proteoform_review_2026-11-22_mane_gap.tsv](proteoform_review_2026-11-22_mane_gap.tsv);
- counts in [proteoform_review_2026-11-22_counts.json](proteoform_review_2026-11-22_counts.json).

## Finding 1: the reference ladder misses 66 MANE Select transcripts

The ladder reads MANE from a column in the APPRIS e115v50 file. That file does not list every
Ensembl 115 MANE Select transcript. APP is an example: its MANE Select is ENST00000346798, the
770-aa protein, and the file omits it.

So 66 library genes look "APPRIS-only" even though Ensembl tags a MANE Select transcript for them.
In 59 of the 66, ours is not the MANE Select protein, and in 37 of those the MANE Select protein is
exactly the canonical. Under the policy of record, MANE and APPRIS disagree in each of these genes,
and a disagreement defaults to MANE unless peptides strongly favor APPRIS.

This one defect accounts for most of what follows:

- 15 of the 18 genes where the canonical has peptide evidence and ours has none;
- 28 of the 30 switch recommendations.

The fix belongs in the ladder: read MANE from the GTF and rerun. Correcting it gene by gene would
mean 59 hand overrides. Of the 59, 45 are in the review table. The other 14 differ from MANE but
not from the canonical, so the review table does not cover them; they are listed in the MANE-gap
table.

## Finding 2: PeptideAtlas classification of the 212

| Class | Genes | In the 206 |
|---|---|---|
| Ours supported only | 32 | 31 |
| Canonical supported only | 34 | 32 |
| Both supported | 67 | 63 |
| Neither | 79 | 79 |

18 genes have a canonical-only peptide seen in ≥ 10 samples and nothing for ours. 9 of them meet
the ladder's three-peptide bar: ADGRB2, CD99L2, GGT5, MUC1, MUC13, PTPRM, SLC39A14, SLC45A4 and
TMEM134. The decision for each of the 18:

| Decision | Genes | Basis |
|---|---|---|
| Keep | ACHE | Named GPI form: the canonical is the abundant non-anchored form, so evidence for it is expected |
| Keep, flagged to Becca | MUC1 | The canonical's 980-residue tandem-repeat array cannot be built as a pooled ORF. Our 273-aa APPRIS form keeps the SEA domain, the TM helix and the tail. Ensembl tags no MANE Select. |
| Switch to the canonical | ADGRB2, CD99L2, CDH3, GGT5, ILDR2, MUC13, NTRK3, OSMR, PCDH9, SLC39A14, TMEM134, TMEM94, TSPAN7 | Finding 1: the MANE Select protein is the canonical. Exact transcripts are in the table. |
| Switch to the canonical | PTPRM | MANE (P28827-2) and APPRIS (P28827, ENST00000332175) disagree. The APPRIS form is a pure deletion, so the ladder saw no evidence. Counting junction peptides, the canonical has 7 (4 in ≥ 10 samples, best 28) and ours has none, which meets the ladder's bar. |
| Revisit | GPNMB, SLC45A4 | The missed MANE Select is neither ours nor the canonical. For SLC45A4, v1 also puts the N-terminus on opposite sides (ours out, canonical in; 12 TM each). |

## Finding 3: the 11 references never reviewed

The proteoform audit flags canonical residues a reference lacks, so it missed these 11 references,
which are more than 10% longer or shorter than the canonical. The audit now records the length
difference in both directions (below).

In the table:

- **Peptides ours / canonical** is the count of distinguishing peptides, with the best sample
  count in parentheses.
- No gene has a GENCODE Ribo-seq ORF overlapping either form's unique segments.
- No gene loses a Pfam family.
- None changes signal peptide, TM count or orientation, except where the table says so.

| Gene | Ours vs canonical (aa) | Where they differ | Peptides ours / canonical | Recommendation |
|---|---|---|---|---|
| KCNN2 | 791 vs 579 | N-terminal 212 aa, cytoplasmic | 0 / 0 | **Keep.** MANE+APPRIS (ENST00000673685). The extended form is a described, functional SK2 isoform that reaches the surface (SK2-L; [Strassmaier et al., 2005](https://pubmed.ncbi.nlm.nih.gov/15797870/), in mouse). |
| PLPPR2 | 452 vs 343 | Cytoplasmic C-terminus | 9 (15) / 3 (8) | **Keep.** MANE+APPRIS (ENST00000688289). |
| PCDH7 | 1,255 vs 1,069 | Cytoplasmic C-terminus | 0 / 1 (5) | **Keep.** MANE (ENST00000695919). The evidence is below the ladder's bar. |
| PCDH1 | 1,237 vs 1,060 | Cytoplasmic C-terminus | 93 (377) / 14 (17) | **Keep.** MANE+APPRIS, Q08174-2. |
| HRH2 | 422 vs 359 | Cytoplasmic C-terminal 63 aa | 0 / 0 | **Keep, low confidence.** MANE+APPRIS (ENST00000636584). The extension rests on the transcript model alone. |
| KCNC1 | 585 vs 511 | Cytoplasmic C-terminus | 10 (59) / 0 | **Keep.** MANE+APPRIS, P48547-2. |
| GPM6B | 328 vs 265 | N-terminal 40 aa, cytoplasmic | 47 (247) / 6 (58) | **Keep.** MANE+APPRIS, Q13491-4. |
| TSPAN17 | 329 vs 270 | Cytoplasmic C-terminal 59 aa | 3 (3) / 0 | **Keep.** MANE+APPRIS, Q96FV3-4. |
| TMEM145 | 563 vs 493 | Cytoplasmic C-terminus and a 14-aa insertion | 0 / 0 | **Switch to Q8NBT3, ENST00000301204.** That transcript is the missed MANE Select (Finding 1). With no evidence either way, the policy default is MANE. |
| LY6H | 161 vs 140 | 21 aa upstream of the signal peptide | 0 / 0 | **Switch to O94772, ENST00000430474,** the APPRIS principal. The mature GPI protein is identical. This is the same case as RGMB (`upstream_extension`, switched 2026-11-12). |
| GFRA4 | 269 vs 299 | Inside the second GDNF-binding domain | 0 / 1 (1) | **Keep, low confidence.** MANE+APPRIS, Q9GZZ7-2. Both ends and the GPI signal are shared. |

Literature was checked for KCNN2 only. For the other ten, the recommendation rests on transcript
support, peptides and topology.

## Proteoform audit change

`audit_l1_proteoforms.py` gains two columns, `canonical_length_diff_fraction` and
`length_differs_from_canonical` (more than 10% in either direction), and a summary count. These
are output-only additions: no flag or flag class changes. GFRA4 is 10.03% shorter, just over the
threshold.

`tests/test_proteoform_length_flag.py` checks two things:

- the helper is symmetric and exact at the boundary;
- run on the committed audit table, it catches exactly these 11.

The committed `l1_proteoform_audit.tsv` was not regenerated. The new columns appear at the next
audit run.

## Open

- The 14 genes kept as open decisions (decision log) need a call between the current form and
  MANE Select: both are observed, or the current one only weakly.
- 6 isoform-panel slots, left by the retired duplicates, need refilling by the isoform redraw.
- 11 genes have a v1 topology difference between ours and the canonical, listed in the counts
  JSON. Apart from SLC45A4, this review does not adjudicate them.

## Decision log

**2026-11-22.** Becca approved rerunning the ladder with MANE Select read from the GTF in this
build ("Yes, in this build"), with the three-way calls delegated: MANE unless observed evidence
is strong for the other form, ambiguous genes kept and listed. Applied by
[apply_mane_gtf.py](apply_mane_gtf.py); the MANE table is
[derive_mane_select.py](derive_mane_select.py)'s, pinned in `library_pipeline_inputs.json`.

| Outcome | Genes |
|---|---|
| Switched to MANE Select | 42, including TMEM145 |
| Switched to a named non-MANE form | PTPRM (P28827), LY6H (O94772) |
| Kept: evidence favors the current form | ADGRL2, MS4A6A |
| Kept: open decision | ADGRG6, APP, ATP1A3, CSPG5, DPP10, DPP6, GRIN1, KCNMA1, KCNQ2, NFASC, ROBO2, SIRPA, SLC12A3, SLC4A7 |
| Unchanged (MANE already built, or named earlier) | the other 8 of the 66 |

The per-gene calls and their peptide counts are in
[mane_gtf_decisions_2026-11-22.tsv](mane_gtf_decisions_2026-11-22.tsv), and the 44 changed
members in [mane_gtf_changes_2026-11-22.tsv](mane_gtf_changes_2026-11-22.tsv). Six isoform-panel
members became identical to their gene's new reference and are retired as duplicates (ADGRB2,
LSR, MUC13, PCDH9, STEAP3, TMEM134); the isoform redraw refills those slots.
