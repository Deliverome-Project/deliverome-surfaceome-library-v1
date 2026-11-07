# Coding-sequence integrity and expression variants

The L1 and L2 coding-sequence audit saves engineered expression variants separately
from the preserved input CDS/construct records, called native records here. Earlier
recorded edits, including appended stops, retain their source metadata. A separate Compara 116 candidate revision contains
the 200 selected orthologs; neither that revision nor the initiation edits activate
the canonical manifest. Tags and barcodes remain downstream of these pre-tag designs.

## Candidate scope and frozen audit

The ortholog revision retains 11 existing L1 proteins, promotes 139 exact-protein
L2 members, and introduces 50 new proteins. The 189 additions receive permanent
L1 IDs `l1-p2743`–`l1-p2931`. Replaced L1 rows remain as retired records; promoted
L2 rows retain their original IDs with `moved_to` reasons. The allocation ledger
reserves the new IDs even while the canonical manifest remains unchanged.

Every active member in each manifest is audited, including L2. Retired rows are
counted explicitly rather than silently disappearing from the record.

| Frozen snapshot | Active L1 | L1 non-ATG | Active L2 | L2 non-ATG | Retired rows excluded |
|---|---:|---:|---:|---:|---:|
| Current manifest | 2,735 | 7 | 4,279 | 53 | 68 |
| Compara 116 proposal | 2,735 | 7 | 4,140 | 51 | 396 |

Non-ATG means that the first DNA triplet differs from ATG. It is an initiation
review flag, not evidence that translation cannot occur. No other requested
sequence defect occurs in these frozen snapshots: the audited sequences have no
frame errors, ambiguous bases, internal stops, missing terminal stops, or unexpected
missing CDS/ORF records. Raw CDS and domesticated ORF proteins agree wherever both
records apply. The 25 expression-rescue and 24 giant-truncation arms instead have
separately sourced constructed ORFs; the audit records that evidence class and does
not label their absent transcript-CDS entries as missing designs.

| Audit rule | Test and evidence retained |
|---|---|
| `CDS001`, `CDS002`, `CDS301` | Construct/raw-CDS availability, with an explicit source requirement for construct-only rescue and giant arms |
| `CDS101` | First triplet, independently for raw CDS and domesticated ORF |
| `CDS102`, `CDS103` | Length divisible by three; DNA alphabet restricted to A/C/G/T |
| `CDS104`–`CDS106` | In-frame terminal stop, internal stop positions, standard nuclear-code translation |
| `CDS201` | Exact raw-versus-domesticated protein equality, excluding one terminal stop |

The ledger retains member ID, tranche, key, transcript, source, original codons,
nucleotide/protein hashes, explicit engineering metadata and rule hits. Summary
files record all input hashes, script hash and commit/dirty status, commands,
excluded-row counts and output hashes. Audit rule version is `library-cds-audit-v1`.

## Initiation evidence and proposed designs

The union of the baseline flags and additional proposed orthologs contains **61
distinct CDS keys**. Evidence classification uses pinned Ensembl transcript
annotations, with experimental initiation evidence where available. Absence of an
incomplete-CDS flag does not establish completeness.

| Evidence class across the 61 keys | Keys | Proposed expression design |
|---|---:|---|
| Annotated non-AUG start at the first CDS codon | 5 | Replace the first triplet with ATG; insert no extra residue |
| Experimentally supported non-AUG initiation, human and mouse TRPV6 | 2 | Replace the initiating ACG with ATG; insert no extra residue |
| Annotated CDS start not confirmed (`cds_start_NF`) | 11 | Held diagnostic ATG-prepend variant; retain the uncertainty flag |
| No annotated start; CDS begins at the transcript's 5′ boundary | 30 | Prepend ATG; retain unresolved completeness |
| No annotated start; CDS begins within the transcript | 12 | Prepend ATG; retain unknown-start status |
| Transcript annotation unavailable | 1 | Prepend ATG; retain unknown annotation status |

The `cds_start_NF` tag means that the coding start could not be confirmed; it does
not establish biological truncation ([GENCODE tag definitions](https://www.gencodegenes.org/pages/tags.html)).
The 11 tagged cases in this table cover the non-ATG subset, not every active library
member. Annotation review across all starts is a separate audit.

The diagnostic proposal set contains seven first-triplet replacements and 54 ATG-prepend
variants. The two experimentally supported and five annotation-supported replacements
have stronger initiation evidence. The 54 prepend proposals are not accepted repairs;
they remain excluded from canonical activation while transcript alternatives are reviewed. For TRPV6, experiments identify ACG initiation with
methionine; ordinary codon-table translation of that first ACG as threonine does
not describe the initiation event ([Fecher-Trost et al., 2013](https://pubmed.ncbi.nlm.nih.gov/23612980/)).

The 54 prepend decisions add an N-terminal methionine to the selected sequence;
they do not reconstruct missing upstream sequence or certify a full-length native
protein. Completeness labels remain attached to the design. Among the 60 baseline
cases, 26 have another usable ATG transcript from the same stable gene, but none
has an annotated ATG model encoding exactly methionine plus the selected protein.
Those alternatives remain review candidates rather than automatic isoform substitutions.

The decision engine (`cds-boundaries-1`) requires an original-CDS checksum, named
action, rationale and evidence. It rejects ambiguous/out-of-frame inputs and
internal stops instead of repairing them implicitly. Each edit records the original
and designed sequences and hashes, protein effect and unresolved completeness.
Native CDS and construct artifacts remain immutable; engineered records occupy a
separate design revision. A member ID or accession alone therefore does not specify
the protein across these revisions: the revision and exact sequence hash are required.

## Transcript alternatives remain review candidates

Of the 54 held forms, 25 have another transcript with annotated start and stop
features; 29 lack such an alternative in the pinned inputs. This stricter
start-and-stop-supported count covers the 54 held proposals; the 26-model count
above covers all 60 baseline cases under the broader sequence-integrity filter.
Only PRRG1 and NPFFR2
have strict N-terminal replacement candidates under the recorded alignment rule.
OPN3, ADCY2, OR8G3 and SLC9D1 remain blocked. These comparisons do not select a
replacement automatically: changing transcript can erase the intended isoform
contrast, and 12 cases have an alternative protein already present in the proposed
library. The duplicate index covers all 6,875 proposed active members, including
49 separately sourced constructed ORFs. Every held row retains
`curation_required` and canonical activation is disabled.

## Annotation review across all active members

The annotation audit includes ATG-starting sequences, main human references,
alternative isoforms and orthologs. Frozen baseline sequences are compared with
pinned Ensembl 115 annotations; this comparison does not assert that every original
sequence came from that release. The 50 newly resolved forms use their documented
release-116 inputs. Supplementary release-116 comparisons remain separate from
primary member counts, and annotation evidence is not transferred across a CDS mismatch.

| Snapshot | Active members | L1 members requiring review | L2 members requiring review | Construct-only arms |
|---|---:|---:|---:|---:|
| Current | 7,014 | 36 | 110 | 49 |
| Compara 116 proposal | 6,875 | 35 | 108 | 49 |

These review queues combine initiation, annotation availability and transcript-class
flags; they are not counts of broken proteins. In both snapshots, 13 forms have
`cds_start_NF`: one L1 isoform and 12 L2 forms. Two of the L2 forms begin with ATG:
KCNMA1 isoform ENST00000404857 and mouse FGFR2 ortholog ENSMUST00000119260.
No mapped transcript carries `cds_end_NF`; 26 transcripts lack a matching GTF record,
so their annotation status remains unknown. Seventeen mapped transcripts lack an
annotated stop feature even though the library sequences contain a terminal stop.
The audit keeps a recorded added stop distinct from native annotation evidence.

The baseline also contains 37 transcripts classified as nonsense-mediated decay and
two classified as protein-coding loss-of-function; the proposal contains 36 and two,
respectively. These labels trigger curation, not automatic exclusion. The two
loss-of-function annotations belong to the main human references CD300H and OR8G3.
Counts overlap across flags. Each row retains exact CDS identity, transcript version,
strand, assembly, gene, coding-feature length, phase, support level and spliced
start/stop coordinates. The 49 constructed arms retain their separate source class.

## Priority form recommendations

A first-codon replacement changes the existing triplet to ATG; a prepend inserts
an additional methionine. For the three human reference cases, replacement yields
the intended reviewed protein exactly, whereas a prepend does not. Evidence strength
still differs, so the recommendations below remain unapplied.

| Form | Evidence and recommendation |
|---|---|
| Human TRPV6 | Experimental non-AUG initiation supports replacing ACG with ATG while preserving the selected 765-residue form. |
| Human SLC22A17 | CTG-to-ATG matches the selected 631-residue Q8WUG5-2 protein exactly; annotated initiation and curator inference support review of this replacement. The 649-residue canonical isoform is a different form. |
| Human OR8G3 | GTG-to-ATG matches reviewed P0DMU2, but suggested GUG initiation and a loss-of-function transcript label conflict. Resolve this before accepting a reference-form rescue. |
| Cyno OPN3 | Available ATG models differ substantially from the selected counterpart. Seek independently supported full-length sequence or review the ortholog draw. |
| Cyno ADCY2 | Added methionine is an inference from alignment; available ATG models have internal differences. Validate initiation or review the selected ortholog form. |
| Cyno PRRG1 | Annotated alternative N termini merit form-repick review; they are distinct proteins, not accepted boundary repairs. |
| Human SLC9D1 alternative isoform | Preserve the screen's isoform distinction. The 127-residue form has unconfirmed coding and transcript starts; replacing it with the already represented 677-residue canonical protein would erase that distinction. |

The priority evidence retains raw reviewed UniProt records from release 2026_04,
including comments and evidence codes, alongside the reference-selection record and
exact sequence comparisons. The next curation queue starts with flagged main human
references and the three selected cyno orthologs, then the remaining L1 and L2 forms.
Prediction work remains held until actual sequences are selected.

## Construct and prediction consequences

The 61 engineered variants are domesticated and checked against their intended
proteins, start/stop boundaries and assembled restriction-site patterns. This
sequence verification does not establish expression, surface display or restoration
of the native N terminus. Reapplying all mechanical CDS rules to each engineered
snapshot gives zero flags in L1 and L2. No new within-tranche protein duplicate
involves an engineered key. These are mechanical checks, not biological validation.

Native D1 predictions remain historical records. A native accession, similar length
or inherited topology label does not authorize transfer to an engineered sequence.
Prediction reuse requires the identical protein hash and the intended tool/version;
topology requests specify DeepTMHMM v1 and retain its reported model version.
The build's pinned starting topology snapshot contains seven exact engineered-protein
matches and 54 gaps (29,416 residues); `missing_deeptmhmm.fasta` is the corresponding
worklist. New prediction results require their own exact-sequence provenance. Native,
engineered and curated evidence remain distinguishable.

## Reproduction and provenance

All stages write fresh review directories. The recorded commands and manifests,
rather than a live database query or an unrecorded manual filter, define the replay.

| Stage | Reproducible method | Review records |
|---|---|---|
| Ortholog candidate adoption | `run_compara116_adoption.py`; pinned `compara116_cds_inputs.json` | Run journal, member map, `compara116_id_allocations.tsv`, full candidate manifest/CDS/construct revision and 200 pre-tag vectors |
| All-member annotation audit | `audit_library_completeness.py`; `library_completeness_inputs.json` | Full member ledger, review queue, grouped counts and exact annotation/sequence provenance |
| Priority form review | `review_priority_cds_forms.py`; `priority_cds_forms_inputs.json` | Seven evidence-bound recommendations with native/intended protein hashes and raw source records |
| All-member sequence audit | `audit_library_cds.py`, once with each snapshot's manifest/CDS/construct triple | `library_cds_audit.tsv.gz`, `library_cds_issues.tsv`, `library_cds_transitions.tsv`, `library_cds_audit_summary.json` |
| Native-start evidence and alternative-model filters | `audit_library_non_atg.py`; `library_non_atg_inputs.json` and `compara116_non_atg_inputs.json` pin releases, S3 sources and checksums | `library_non_atg_evidence.json`, `library_non_atg_repair_candidates.tsv`, `library_cds_boundary_decisions.tsv` |
| Alternative adjudication | `audit_non_atg_alternatives.py`; `non_atg_alternatives_inputs.json` | Full alignment differences, start/stop evidence, exact topology matches and active-member duplicate checks |
| Checksum-bound edits | `repair_cds_boundaries.py`, using the explicit 61-row decision table | `cds_boundary_overrides.json`, `cds_boundary_edit_ledger.json` |
| Engineered construct verification | `build_cds_boundary_variants.py`, with pinned backbone and topology snapshot | `engineered_candidate_cds.json`, `engineered_candidate_constructs.json`, `cds_boundary_variants.tsv`, build summary and prediction FASTAs |
| Final code revision | Run journals pin each producing commit; [audit source](https://github.com/Deliverome-Project/deliverome-analysis/blob/c8c72328/analysis/sprint4/s4e9-orf-design-fg-panel/audit_library_cds.py) | Source, tests and exact rules used for the published run |
| Pull request | [Candidate adoption and CDS audit PR](https://github.com/Deliverome-Project/deliverome-analysis/pull/260) | Review of adoption, filters, decisions and remaining limitations |
| Durable review bundle | `processed://sprint4/s4e9-orf-design-fg-panel/library-design/cds-boundary-review/2026-11-05/` | Preserve native inputs, derived evidence, decisions, engineered outputs and run provenance together |

The annotation audit has 13 focused tests, including split codons on both strands, transcript-version/strand disagreement and coding-length conventions. Its clean replay uses [source revision 3b3e240f](https://github.com/Deliverome-Project/deliverome-analysis/tree/3b3e240f28609d26ff319cf62475f5bd50e61ff8), with the full report and priority evidence at `processed://sprint4/s4e9-orf-design-fg-panel/library-design/completeness-review/2026-11-05/`.

The sequence audit has 12 focused tests. Native-start evidence tests cover annotation
classification and transcript-model comparisons; boundary-edit tests check checksum
binding, protein effects and preservation of uncertainty. Final test/run results
belong to the published run record, not an assumption from the existence of outputs.

## Remaining decisions

| Question | Why it matters | What resolves it |
|---|---|---|
| Are the partial/unknown N termini suitable for the intended screen? | An added initiator does not restore missing domains, targeting signals or native regulation | Review retained evidence flags and, where needed, choose a documented model or validate the engineered arm experimentally |
| Are exact engineered-sequence predictions complete? | Native predictions may describe a different protein | Resolve which variants to keep before running the held worklist; retain exact-sequence prediction provenance |
| Which revision becomes the canonical member/construct record? | Candidate IDs span preserved native and engineered design revisions | Review the versioned bundle and explicitly activate the chosen revision; tag/barcode construction remains downstream |

## Decision log

| Date | Was | Now | Why |
|---|---|---|---|
| 2026-11-05 | Ortholog adoption and boundary integrity evaluated in separate candidate stages | Full L1/L2 audit accompanies the versioned 200-ortholog adoption and permanent ID ledger | Existing L2 CDSs need the same checks as new L1 picks; source records remain recoverable |
| 2026-11-05 | Non-ATG flags unresolved in native records | Seven explicit start replacements and 54 explicit ATG-prepend expression variants in a separate revision | Provide a reproducible expression design while preserving annotation uncertainty and native sequence/prediction provenance; canonical activation is separate |
| 2026-11-05 | Initiation review covers non-ATG forms | All active L1/L2 members receive annotation review; 61 diagnostic variants remain held | ATG does not establish a confirmed coding start, and intended isoforms must survive curation |
