# L1 proteoform audit: every active member's designed protein, checked against its gene

Each of the 2,735 active L1 members is ordered as one protein: its designed coding sequence (CDS)
translated. This audit asks, for every member, whether that protein is the one the library meant
to order. Is it the gene's reference form? Is it intact? Does it encode what its source record
says? For orthologs, is it a sound model of the human protein? The audit itself is a report.
Every member that fails a check gets a decision in a separate table, decided on the evidence
under Becca's delegation of protein-level calls (2026-11-12). The single writer for the
manifest applies them.

## Result

| Outcome | Members |
|---|---:|
| No flag that needs review | 2,651 |
| Objective defect, fixed by rule | 1 |
| Judgement call, decided on the evidence | 83 |
| of which: change (swap, repair or redraw) | 13 |
| of which: keep (flag explained, no change) | 70 |

The library is clean against every objective check except one start codon. No member has an
internal stop or a broken frame. Every member's CDS translates to its stated protein. No two
active members order the same protein, and every active isoform passes the cohort gate. The
substantive problems are 13 members whose protein is a fragment, lacks part of its gene's
architecture, or carries a gene-model error. Each gets a named replacement, repair or redraw.

The table above counts the members as they were decided, against the manifest before any change
was applied. `l1_proteoform_decisions_check.tsv` checks every decision against the current
manifest. In that manifest, the writer has applied the changes on PR 283.

| Decision | Outcome in the current manifest | Members |
|---|---|---:|
| Swap (TSPAN10, DISP2, LMLN) | Applied; each translates exactly to its new form | 3 |
| Repair (OR8G3, EPHA6, cyno CHRM3, CCR5, OR3A2, NIPAL4) | Applied. OR8G3's start substitution and the 5' trims (EPHA6 94 residues; cyno CHRM3 31, CCR5 27, OR3A2 32, NIPAL4 214) are declared, so each matches its source protein. The four cyno restarts are N-out under DeepTMHMM v1 (snapshot `2026-11-12-cyno-restarts`) | 6 |
| Redraw (CLEC7A, OPN3, SLC6A11, ABCC3, NECTIN2 isoforms) | Dropped. The panel rule redrew the isoform stratum and displaced 7 more isoforms. The 12 isoforms it added all audit clean | 5 |
| Keep | Unchanged | 70 |

After the apply, the audit flags exactly the 74 keeps: the 70 above, plus four that the
provenance review below found without a decision row. No objective defect remains, and every
change reached the manifest.

## What the audit checks

**The member's protein** is its designed CDS translated: the CDS record at
`processed://sprint4/s4e9-orf-design-fg-panel/library-candidates/2026-11-12-reference-fixes/candidate_cds.json`.
The 25 expression-rescue arms are the exception. They are ordered as an Igκ leader fused to the
receptor, so their protein is the translated construct. Every comparison is by exact sequence,
never by accession or by length.

**The comparator** is the protein each member is judged against.

| Member set | Count | Comparator |
|---|---:|---|
| Library reference | 2,360 | The gene's UniProt canonical |
| Expression-rescue arm | 25 | The gene's UniProt canonical |
| Isoform | 150 | The gene's L1 library reference |
| Ortholog (mouse 100, cyno 100) | 200 | The human L1 library reference |

The UniProt canonical is the comparator for references because it is the independent third
annotation. The reference ladder chooses between MANE Select and the APPRIS principal. A
difference from UniProt therefore marks a member to look at, not an error.

**Per-member columns** in `l1_proteoform_audit.tsv`:

| Group | What is recorded | Source |
|---|---|---|
| Build source | Ladder decision for the gene. Whether the member equals the MANE, APPRIS principal or UniProt canonical protein. Whether MANE and APPRIS agree on the protein sequence | `reference_selection.tsv`; Ensembl 115 peptides; UniProt 2026-10-11 isoforms |
| CDS | Start codon (and whether it was replaced), internal stops, frame. Whether the translation equals the stated protein | Ensembl 115 peptide (human); Ensembl 116 (mouse, cyno); rescue record |
| Length | Fraction of the comparator aligned. Comparator residues missing at each end. Member residues before the comparator's Met1 | Global alignment, BLOSUM62, gap 11/1, free end gaps |
| Signal peptide and topology | Whether the member contains the UniProt signal peptide. DeepTMHMM v1 signal peptide, TM count and N-terminal side, for member and comparator | UniProt 2026_04 signal peptides; DeepTMHMM v1 pinned snapshot, matched by sequence hash |
| Domains | Pfam families on member and comparator, and families lost | Pfam 39.0, gathering threshold |
| Evidence | APPRIS label and MANE flag of the member's transcript. PeptideAtlas peptides that distinguish member from canonical. GENCODE Ribo-seq ORFs inside the member. Isoform cohort tier | APPRIS e115v50; PeptideAtlas build 607 (genes in the pinned pool only); GENCODE Ribo-seq ORFs |
| Ortholog | Identity and coverage against the human reference. N-terminal extension beyond human Met1. Human signal peptide absent in the ortholog | Alignment above; DeepTMHMM v1 |
| Structure | Model accession, and whether the model is this exact protein | `proteoform_match.tsv`, joined on member id |
| Duplicates | Other active L1 or L2 members with the identical protein | Manifest |

## Flag classes

**Apply** flags are objective defects: a broken start codon, an internal stop, a CDS that does
not translate to the stated protein, a protein that duplicates another member, or a failed
cohort gate. The fix follows from the defect. **Proposed** flags mark a scientific question:
they route a member to hand review and decide nothing by themselves. **Info** flags are
reported but route nowhere.

The numeric cutoffs on proposed flags are screening thresholds, chosen to catch the cases Becca
had already found by hand (TFRC and SLC43A2 fragments, CD200R1L / LRRC37B / CACHD1 truncations,
cyno N-terminal gene-model errors) without flagging ordinary sequence variation. A member near a
cutoff is not therefore suspect.

| Flag | Class | Sets checked | Fires when | Members |
|---|---|---|---|---:|
| `non_atg_start` | apply | all | CDS does not start ATG | 1 |
| `internal_stop` | apply | all | Stop codon inside the CDS | 0 |
| `cds_frame` | apply | all | CDS length not a multiple of 3 | 0 |
| `translation_mismatch` | apply | all | Translation differs from the stated protein beyond declared edits | 0 |
| `duplicate_protein_l1` | apply | all | Identical protein to another active L1 member | 0 |
| `cohort_gate_fail` | apply | isoform | Isoform fails the 2026-11-12 cohort gate | 0 |
| `fragment` | proposed | all | Shorter than half the comparator, or < 100 aa against a comparator ≥ 200 aa | 6 |
| `sp_lost_vs_canonical` | proposed | reference, rescue, isoform | Member lacks the UniProt canonical's signal-peptide sequence | 9 |
| `n_terminal_truncation` | proposed | reference, rescue | ≥ 20 comparator residues missing at the N-terminus | 19 |
| `c_terminal_truncation` | proposed | reference, rescue | ≥ 20 comparator residues missing at the C-terminus | 16 |
| `reference_short` | proposed | reference, rescue | < 90% of the canonical aligned | 21 |
| `tm_count_change` | proposed | reference, rescue, ortholog | DeepTMHMM v1 TM count differs from comparator | 10 |
| `orientation_change` | proposed | all | DeepTMHMM v1 N-terminal side differs from comparator | 3 |
| `tm_all_lost` | proposed | isoform | Comparator has TM helices, the isoform none | 0 |
| `pfam_domain_lost` | proposed | reference, rescue, ortholog | A comparator Pfam family is absent | 15 |
| `pfam_all_lost` | proposed | isoform | No comparator Pfam family left | 2 |
| `ortholog_divergent` | proposed | ortholog | Identity < 60% or human coverage < 80% | 5 |
| `ortholog_n_extension` | proposed | ortholog | ≥ 10 residues before the residue aligned to human Met1 | 7 |
| `ortholog_sp_missing` | proposed | ortholog | Human has a v1 signal peptide, the ortholog none | 1 |
| `ladder_source_unrecorded` | proposed | reference | Built protein is not the ladder's selected source (MANE, APPRIS, UniProt canonical or a named override), and no reference fix records why | 4 |
| `no_v1_topology` | info | all | No DeepTMHMM v1 row for this exact sequence | 25 |
| `structure_not_ours` | info | all | No model of this exact protein | 346 |
| `reference_not_mane_appris_canonical` | info | reference | Member equals none of MANE, APPRIS principal, UniProt canonical | 8 |

A member can carry several flags, so the counts sum to more than the 84 flagged members.
Isoforms are not tested for truncation or domain loss: an isoform differs from its reference by
design, and the cohort gate already tests its start, signal peptide, orientation and first TM.
The `fragment` and `pfam_all_lost` checks are the isoform tests for "too little left to be a
surface form".

## Decisions

All 84 flagged members have a row in `l1_proteoform_decisions.tsv`. Each row records the
decision, target, `decided_by`, `implemented_by`, reason, sources and the alternative that was
rejected. Every row has status `apply`, and the manifest writer applies all of them.

On 2026-11-12 Becca delegated protein-level calls: "use your best judgment to complete a new
and improved evidence based l1 completely by the morning". The 83 judgement calls are therefore
decided here and attributed to `l1-proteoform-audit (delegated by Becca 2026-11-12)`. The
exceptions are the four GPI-form keeps, which stay attributed to Becca's own earlier overrides.
The rejected alternative on each row is what to restore if Becca overturns a call.

### Objective defect

| Member (design) | Gene | Defect | Fix |
|---|---|---|---|
| l1-m0550 (l1-p1525) | OR8G3 | Both Ensembl 115 transcripts start GTG, so the ordered protein begins Val | Codon 1 GTG → ATG. The protein becomes reviewed UniProt P0DMU2 exactly. This is the remedy Becca approved for SLC22A17 and TRPV6 on 2026-11-12 |

### Changes

| Member (design) | Gene | Set | Decision | Target | Why |
|---|---|---|---|---|---|
| l1-m0254 (l1-p2284) | TSPAN10 | reference | swap | l2-p1634 (Q9H1Z9) | APPRIS-only Ensembl model, no UniProt form. It loses TM4 (3 of 4 TM helices) and 123 C-terminal residues. The canonical has five Ensembl 115 transcripts and is already in L2 |
| l1-m2587 (l1-p0513) | DISP2 | reference | swap | l2-p0455 (A7MBM2) | Ensembl-only model. 217 canonical residues are missing internally, and 2 of 12 TM helices. The canonical is already in L2 |
| l1-m1324 (l1-p1106) | LMLN | reference | swap | Q96KR4 via ENST00000330198 (new CDS) | Starts at canonical residue 233, with no signal peptide and no TM. As ordered it cannot enter the secretory pathway |
| l1-m2549 (l1-p0578) | EPHA6 | reference | repair | Start at native Met95 → UniProt Q9UF33 exactly | MANE and APPRIS add 94 residues before the signal peptide. DeepTMHMM v1 then reads the signal peptide as a TM helix with the N-terminus inside. Rejected: keep the MANE form. The mouse and cyno models are also 1,130 aa, so the extension may be translated, but the library needs the form with a cleavable signal peptide at residue 1 |
| l1-m0003 (l1-p2420) | CLEC7A | isoform | redraw | Same stratum | Q9BXN2-9, 77 aa, no C-type lectin domain |
| l1-m0111 (l1-p2477) | OPN3 | isoform | redraw | Same stratum | 189 aa, 4 of 7 TM helices, no UniProt isoform |
| l1-m0160 (l1-p2511) | SLC6A11 | isoform | redraw | Same stratum | P48066-2, 208 aa, 3 of 12 TM helices |
| l1-m1723 (l1-p3119) | ABCC3 | isoform | redraw | Same stratum | O15438-5, 572 of 1,527 aa, no nucleotide-binding domain |
| l1-m1617 (l1-p2470) | NECTIN2 | isoform | redraw | Same stratum | Differs from the reference only by one residue deleted in the signal peptide. The mature protein is identical, so the slot adds nothing |
| l1-m1826 (l1-p2938) | CHRM3 | cyno ortholog | repair | Start at native Met32 | 31-residue extension upstream of the Met aligned to human Met1. It has no signal peptide, so it would sit on the extracellular N-terminus. The trimmed protein is exactly as long as human |
| l1-m1070 (l1-p2959) | CCR5 | cyno ortholog | repair | Start at native Met28 | As CHRM3, 27 residues. The trimmed protein is exactly as long as human |
| l1-m0887 (l1-p2963) | OR3A2 | cyno ortholog | repair | Start at native Met33 | As CHRM3, 32 residues. The trimmed protein is exactly as long as human |
| l1-m1820 (l1-p2962) | NIPAL4 | cyno ortholog | repair | Start at native Met215 | 214-residue Lys/Arg/Gly-rich extension, a gene-model error. The trimmed protein is exactly as long as human |

The four cyno repairs are explicit, recorded exceptions to the ortholog CDS policy
(`ORTHOLOG_CDS_POLICY.md`), which forbids implicit start repairs. Each gene has only one
Ensembl 116 model, so no alternative transcript exists to swap to. Each trim starts at a native
in-frame ATG and changes no codon. Rejected: a redraw from the Compara 116 pool, which would
remove these genes' cyno arms, including the therapeutic targets CCR5 and CHRM3.

### Keeps

Each group's flag is explained by a choice the library already made deliberately.

| Group | Members | Why the flag does not call for a change |
|---|---:|---|
| MANE and APPRIS agree; UniProt canonical differs | 29 | Rung 1 of the ladder. The difference is UniProt's, not a defect. Includes NPFFR2, CCR2, the olfactory receptors and AQP12B. For AQP12B, DeepTMHMM v1 calls the MANE form 7 TM N-out (aquaporins are 6 TM, N-in), so tag placement for it should not rest on v1. Also includes the two rescue arms (CCR2, GPR142) built on the MANE receptor |
| MANE kept over APPRIS on evidence | 4 | Rung 3 of the ladder (ANO7, SLC12A5, PVRIG, XG) |
| No MANE; APPRIS principal | 7 | Rung 2 (HTR7, DPP6, CD99L2, NFASC, SLC4A7, MUC1, LSR). The MUC1 form lacks the repeat array, which cannot be synthesized. LSR's canonical is in L2 if wanted |
| Signal-peptide sequence variant | 7 | DeepTMHMM v1 still calls a cleaved signal peptide on the member (ADGRB2, CADM2, CDH3, LILRA6, LILRB3, LINGO1, LYPD6B) |
| Becca's GPI-form overrides | 4 | ACHE, BCAN, NCAM1, TFPI. The "lost" domains belong to the non-GPI canonical |
| Named, membrane-anchored isoforms | 3 | HLA-G P17693-3, GRIA1 P42261-6, SCARA3 Q6AZY7-2 |
| Readthrough reference | 1 | SLCO1B3-SLCO1B7 orders UniProt Q9NPD5-2, a reviewed SLCO1B3 isoform distinct from the SLCO1B3 reference. Rejected: drop it |
| Ortholog extension before a signal peptide | 3 | Cyno FOLR2, LYPD2 and IL9R. The extension is cleaved with the signal peptide |
| Ortholog divergence | 5 | Mouse IL15RA, ICAM1, IZUMO2, MS4A7 and TMEM139. Coverage ≥ 97%, ordinary divergence |
| Ortholog Pfam call below threshold | 6 | Same family at lower copy number or score, or a clan-mate on the same segment |
| Built allele is not the ladder's rung | 4 | OR11H7: GRCh38's chr14 allele is loss-of-function (Q227 stop), so the build uses the intact patch-scaffold allele, identical to UniProt Q8NGC8, the KIR2DS4 / OR12D1 / OR5AL1 pattern without its override row. OR8U8, OR9G9, TAS2R45: alternative-haplotype Ensembl models that differ from UniProt canonical at 3, 1 and 1 declared residues; the build keeps the assembled allele rather than editing UniProt residues in |
| Human reference is the anomaly | 1 | Mouse OR10G6. The human reference carries 21 extra N-terminal residues that v1 calls a signal peptide |

## Calls Becca may want to revisit

Every call can be reversed from the `alternative_rejected` column. These are the ones where the
evidence was closest:

| Member | Call | Rejected alternative | Why it was close |
|---|---|---|---|
| EPHA6 | Repair to UniProt start (Met95) | Keep the 1,130-aa MANE form | MANE, APPRIS and both ortholog models carry the extension |
| Cyno CHRM3, CCR5, OR3A2, NIPAL4 | 5' trim to the conserved Met | Redraw from the pool | The trim is an exception to the ortholog CDS policy |
| LMLN | Swap to the canonical (new CDS) | Keep the APPRIS model | The canonical needs a CDS designed from ENST00000330198 |
| SLCO1B3-SLCO1B7 | Keep | Drop | Its protein is an isoform of another library gene |
| AQP12B | Keep MANE | Swap to the 295-aa canonical | DeepTMHMM v1 calls the MANE form in the wrong orientation for an aquaporin |

## Provenance review

Becca asked on 2026-11-12 that every in-silico design step be checked against its code and
every decision carry a rationale and date. For proteoform selection that is four checks.

**Every ordered ORF is the protein the manifest names.**
[check_l1_built_orfs.py](https://github.com/Deliverome-Project/deliverome-analysis/blob/main/analysis/sprint4/s4e9-orf-design-fg-panel/check_l1_built_orfs.py)
reads the domesticated ORF that is synthesized, not the CDS record, for all 2,735 members of the
order plan. Every one starts ATG, ends in a stop, has no internal stop and is a whole number of
codons. All 2,710 manifest members translate to the manifest's `sequence_sha256` and length,
and the 25 rescue arms equal their Igκ constructs. One downstream table is stale:
`reviewer-package/l1_order_plan.tsv` still carries the pre-trim `insert_nt` for cyno CHRM3,
NIPAL4 and OR3A2 (1,956, 1,947 and 1,134 nt, against 1,863, 1,305 and 1,038). Their plate
order and pooling volumes are therefore wrong until the order plan is regenerated.

**Every deviation from MANE / APPRIS has a record.** The 30 library references whose protein is
neither the MANE nor the APPRIS form are each explained:

| Explanation | Members | Record |
|---|---:|---|
| Named override (GPI forms, MICA, reference-allele LoF, signal-peptide and truncation fixes) | 16 | `reference_overrides.tsv` |
| No MANE or APPRIS (alternative-haplotype loci), so the ladder's UniProt-canonical rung | 9 | `reference_selection.tsv` rung; OR8U8, OR9G9 and TAS2R45 also differ from UniProt canonical at declared residues, which had no decision row and are decided here |
| APPRIS form is a fragment, so the ladder's UniProt-canonical rung (PILRA, SUSD6) | 2 | `reference_selection.tsv` rung |
| Declared fixes (EPHA6 trim, OR8G3 start) | 2 | `reference_changes.tsv` |
| Built from a different allele than the ladder's APPRIS pick, with no record (OR11H7) | 1 | decided here, under the delegation (Keeps above) |

The audit's new `ladder_source_unrecorded` flag makes the check permanent. The ten override rows written before the 2026-11-10 attribution
convention carry no date. Their dates of record are the commits that added them: 2026-10-11
for ACHE, BCAN, CD55, CD58, NCAM1 and TFPI, 2026-10-15 for MICA, and 2026-11-04 for KIR2DS4,
OR12D1 and OR5AL1.

**The isoform gate's drops of SLC9D1 and KCNN2 are right.** Neither dropped form can reach the
surface as the protein its gene makes:

| Form | Ordered | What it is | Call |
|---|---|---|---|
| SLC9D1 ENST00000622371 | 126 of 677 aa | An internal fragment: reference residues 232 onward, on a partial Ensembl model whose CDS opens ATA, with no signal peptide and no TM helix | Drop stands. SLC9D1 keeps no isoform in L1; its gate-passing Q6UWJ1-2 (524 aa) stays in L2, because naming it back would override the stratified draw for one gene with no evidence that it matters more |
| KCNN2 Q9H2S1-2 | 231 of 791 aa | Starts just before the pore loop's GYG selectivity filter (reference ~580), so it keeps the pore loop, S6 and the C-terminal tail but lacks S1–S5; it cannot form a channel | Drop stands |

The 931 gate failures marked "unverifiable" are a different matter (777 full-length-only and
154 observed forms, including KCNN2's 857-aa ENST00000512097). These forms were never tested;
they fail only because DeepTMHMM v1 has not been run on them. The isoform panel is therefore
drawn from forms that already had topology. Running v1 on them is a compute decision, not a
proteoform one, and is not taken here.

**Every input is on S3 with a hash, but not every loader checks it.** Each input the selection
reads exists on S3 and matches a recorded SHA-256:
- the Ensembl 115 human protein, CDS and GTF files, and the Ensembl 116 cyno and mouse protein
  and CDS files;
- APPRIS e115v50 and UniProt 2026-10-11;
- the PeptideAtlas build 607 pools and the six GENCODE Ribo-seq ORF files;
- Pfam 39.0, the DeepTMHMM snapshot, the CDS record and the reference record.

This audit now reads the current reference record (`2026-11-12-reference-fixes-2`, not the
pre-fix `library-design/reference_selection.tsv`). It also refuses any reference record, CDS
record or topology snapshot whose hash differs from the pinned one. Switching record changed only
the `ladder_selected_source` column, for 11 members now on the named rung. The gaps that remain
are in other code, and are listed for the input-manifest owner:

| Gap | Where |
|---|---|
| `topology.load_snapshot()` does not check `SNAPSHOT_SHA256` | `src/deliverome_analysis/topology.py` (checked only by the isoform gate and the pipeline) |
| The CDS record has no pinned hash in `library_design.py` | `CDS_KEY` |
| Live fetches by default: the PeptideAtlas query, the UniProt REST stream (no release pinned), and the APPRIS and Ensembl downloads | `select_reference.py` (without `--offline`), `audit_reference_signal_peptides.py`, `fetch_inputs.py` |
| No MANE release is recorded; MANE comes from the APPRIS file's column | input manifest |
| The Ribo-seq ORFs exist only in `processed://`, with no external copy or release label | input manifest |
| `candidate_inputs_provenance.json` lists URIs without hashes; the signal-peptide table's S3 copy is not in `bundle_manifest.tsv` | input manifest |

## Limitations

- **Thresholds are screens.** A member just inside a cutoff was not reviewed by hand. The full
  table carries every measurement, so a stricter screen can be rerun from it.
- **Topology is DeepTMHMM v1 only.** v1 is the topology of record. Where it disagrees with
  UniProt or biology (AQP12B, EPHA6, OR10G6 human), that is stated as a prediction, not
  evidence. No new predictions were run.
- **Evidence coverage is uneven.** The pinned PeptideAtlas pool covers only the genes the
  reference ladder adjudicated. GENCODE Ribo-seq ORFs are mostly non-canonical ORFs. Absence of
  peptides or ORFs is therefore not evidence against a form.
- **346 members have no model of their own protein.** 193 have a model of a different protein
  (188 references, mostly AlphaFold models of a UniProt canonical that differs from the MANE /
  APPRIS form ordered), 23 references have a longer model that contains their protein, and 130
  have no matched model (67 orthologs, 37 isoforms, 25 rescue arms, 1 reference). The audit
  reports these and decides nothing; folding them belongs to the tag-refold work.
- **Members are named by member id with the design id beside it.** The audit keys on
  `design_id`, as all design-stage code does; `candidate_id` is the `l1-mNNNN` member id.

## Provenance

| Item | Location |
|---|---|
| Audit script | [audit_l1_proteoforms.py](https://github.com/Deliverome-Project/deliverome-analysis/blob/main/analysis/sprint4/s4e9-orf-design-fg-panel/audit_l1_proteoforms.py) |
| Decision script | [l1_proteoform_decisions.py](https://github.com/Deliverome-Project/deliverome-analysis/blob/main/analysis/sprint4/s4e9-orf-design-fg-panel/l1_proteoform_decisions.py) |
| Built-ORF check | [check_l1_built_orfs.py](https://github.com/Deliverome-Project/deliverome-analysis/blob/main/analysis/sprint4/s4e9-orf-design-fg-panel/check_l1_built_orfs.py) → `l1_built_orf_check.tsv` |
| Outputs | `l1_proteoform_audit.tsv` and `l1_proteoform_audit_summary.json` (current manifest), `l1_proteoform_decisions.tsv` (the record, with member columns as decided), `l1_proteoform_decisions_check.tsv` (each decision against the current manifest) |
| Pinned inputs and checksums | `l1_proteoform_audit_inputs.json` (pipeline-input snapshot 2026-11-01, Ensembl 116 ortholog peptides, Pfam 39.0, CDS record, topology snapshot) |
| Pull request | [PR 300](https://github.com/Deliverome-Project/deliverome-analysis/pull/300), stacked on [PR 283](https://github.com/Deliverome-Project/deliverome-analysis/pull/283) (renumbering, reference fixes, isoform redraw) |

Rerun: `uv run --locked --group domains python analysis/sprint4/s4e9-orf-design-fg-panel/audit_l1_proteoforms.py`,
then `uv run --locked python analysis/sprint4/s4e9-orf-design-fg-panel/l1_proteoform_decisions.py`,
then `uv run --locked python analysis/sprint4/s4e9-orf-design-fg-panel/check_l1_built_orfs.py`.

## Decision log

| Date | Decision | Reason | Decided by / implemented by |
|---|---|---|---|
| 2026-11-11 | Audit every active L1 member's proteoform, and record overrides of MANE / APPRIS | Becca's overnight request | Becca Carlson / Claude Code session l1-proteoform-audit |
| 2026-11-12 | OR8G3 codon 1 GTG → ATG, sent to the manifest writer | Objective defect under the overnight rules; same remedy as SLC22A17 / TRPV6 | Rule (Becca, 2026-11-11) / Claude Code session l1-proteoform-audit |
| 2026-11-12 | 13 changes and 70 keeps proposed | Judgement calls went to Becca under the overnight rules | No decision taken (proposed) / Claude Code session l1-proteoform-audit |
| 2026-11-12 | All 83 judgement calls decided and sent to the manifest writer as apply; the cyno extensions changed from redraw to repair, and EPHA6 from "repair or keep" to repair | Becca delegated protein-level calls ("use your best judgment to complete a new and improved evidence based l1 completely by the morning") | l1-proteoform-audit (delegated by Becca 2026-11-12) / Claude Code session l1-proteoform-audit, PR 300 |
| 2026-11-12 | IGSF5 (l1-p0906, Q9NSI5, 407 aa) keeps its MANE / APPRIS / UniProt Met1; no start change sent to the writer. Referred by the type III N-tag session because the ordered form has no signal peptide ahead of a 266-residue Ig ectodomain | The start is not wrong: Met1 is the only annotated start (one human Ensembl 115 protein), and no human frame restores a signal peptide. Mouse (ENSMUSP00000109425, `MEGSWRDVLAVLVILAQLTASG`) and cyno (ENSMFAP00000002510, `…GKAVLTVLVVLAGLAGSGSG`) carry one, but human exon 2 encodes `SAAGLRWWQTAVVDG` in frame. No single-nucleotide deletion makes it hydrophobic, and the +1-frame Met aligned with cyno Met17 reads `MEVSGGSAMVANSSGG`. So the loss is sequence divergence, not a misannotated start. Rejected: a start at Met70 or Met76 (`GWKLIMWALSDMVVL`). That stretch lies inside the I-set domain (Pfam 39–140), is conserved with mouse `GWKLLMWTLNQMVVL` as an Ig strand, and starting there deletes 31–37 residues of the domain. DeepTMHMM v2's 38-aa signal peptide over residues 1–38 has no hydrophobic core and is set aside; SignalP 6 gives P(SP) 0.000 and DeepTMHMM v1 is the record. Adding a leader, e.g. Igκ replacing residues 1–34 so the ectodomain starts at S35 `SGSGNEVIEGPQN`, the mature N-terminus shared with mouse and cyno, is an engineered variant like the OR rescue arms, not a proteoform correction. It stays with the type III tag design | l1-proteoform-audit (delegated by Becca 2026-11-12) / Claude Code session l1-proteoform-audit, PR 300 |
| 2026-11-12 | Re-audited the manifest after the writer applied the 14 changes on PR 283 (head 5b97ddfc). 11 changes are applied, 3 cyno restarts are held for v1 topology, and 5 redraws are dropped. The audit now accepts a declared 5' start trim as an exact translation, and the decision record keeps its as-decided member columns, with a separate post-apply check | Verify that every decision reached the manifest without rewriting the record | l1-proteoform-audit (delegated by Becca 2026-11-12) / Claude Code session l1-proteoform-audit, PR 300 |
| 2026-11-12 | Re-audited on PR 283 head 13f3c1b8, after the writer restarted cyno CHRM3, OR3A2 and NIPAL4 (fa21edc1). All 14 changes are applied and verified. The audit flags only the 70 keeps | Close the post-apply check | l1-proteoform-audit (delegated by Becca 2026-11-12) / Claude Code session l1-proteoform-audit, PR 300 |
| 2026-11-12 | Provenance review on main. All 2,735 ordered ORFs check against the manifest. Four undocumented ladder deviations are decided as keeps (OR11H7, OR8U8, OR9G9, TAS2R45). The gate's drops of SLC9D1 ENST00000622371 and KCNN2 Q9H2S1-2 stand. The audit reads the current reference record with pinned hashes. Stale order-plan lengths, undated override rows and loader gaps are sent to the input-manifest owner | Becca's overnight request to check every design step against its code and record every decision with rationale and date | l1-proteoform-audit (delegated by Becca 2026-11-12) / Claude Code session l1-proteoform-audit |
