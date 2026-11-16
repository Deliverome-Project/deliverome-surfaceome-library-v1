# L1 design review: cover note on work in progress

The [design brief](BRIEF.md) describes the L1 library as it stands. This note tracks what is
still being worked on before L1 is ordered, separately from the brief, so the brief can be read
as the design and this note updated as items close. The first table explains the six design
issues reviewers are most likely to ask about; the tracker below it lists every open item as an
issue, with its owner. Status is as of **2026-11-14**, the final regeneration of the order of
record; rows closed since 2026-11-05 keep their row with the date and the decision. Comments on any of them are welcome; they
are not the questions we are asking reviewers to settle.

Terms follow the brief: a *member* is one construct, a *rung* its tag position, and counts
refer to tranche L1, active members, DeepTMHMM v1 topology. Every count is measured by
[`brief_numbers.py`](brief_numbers.py).

| # | Issue | Why it matters | What resolves it | Status |
|---:|---|---|---|---|
| 1 | **Structure models of a different protein.** 198 members have an AlphaFold model of another proteoform, and 26 a model of a longer one that contains theirs. | Tag-site coordinates would be numbered against a sequence the construct does not contain, and nothing downstream flags it. | Fold each "different" member from its designed CDS; verify and apply the offset for the 26 "contained" ones. | closed 2026-11-14: every ordered insert is refolded from its own tagged protein, and an independent check matches all 2,629 folds to the order of record residue for residue (brief section 12) |
| 2 | **Internal-rung cargo**: ALFA alone (19 aa) or the 51-aa SpyTag003/ALFA tandem. | At 19 aa the cassette fits 16-residue loops (e.g. the SLC16A9/12/14 block); the tandem does not, so this decides whether whole paralog groups have a site. | Decide whether internal members need SpyTag003, then calibrate the loop-length floor at that cargo size. | closed: the full 51-aa cassette is used at every rung, internal included; each internal site is kept only if the tagged protein passes refold (brief section 5) |
| 3 | **Loop gates** `MIN_LOOP_LEN = 8` (minimum contiguous coil/turn run) and `LOOP_SS` (coil, turn, bend or 3₁₀ only). | They change the top-ranked site for 56% and 72% of internal members respectively. In a published-junction benchmark, 28% of large-cargo insertions sit in loops shorter than 8 (n = 18). | Wet data from a ladder of coil-run lengths; the benchmark corpus is 78% small epitope tags and cannot settle it. | superseded 2026-11-12: internal sites come from tag-site model v5, the membrane-side gate and refold, not from these gates (brief section 5) |
| 4 | **GPI single-pass members and the internal worklist.** The member-level ladder reroutes BST2 and RGMB to the internal rung, but the per-gene internal-site worklist and its site generator are topology-only and do not hold them. | A GPI single-pass protein needs an ectodomain site, not a terminal or inter-TM one; until the site lanes cover it, these members have a rung but no junction. | Build the GPI internal set from the anchored genes and add the GPI branch to the per-gene rung derivation. | closed 2026-11-12: BST2 is tagged in its ectodomain on the internal rung; RGMB, with no two-sided loop, is routed to N after SP |
| 5 | **Vendor manufacturability scan** (GC windows, long repeats, homopolymers, hairpins) not yet run on L1. | Synthesis houses quote and schedule on these features; finding them before the order is cheaper than a failed region. | Run the scan over the L1 vectors and separate backbone regions already synthesized from new ones. | closed 2026-11-13: Twist scored the first draft; local rules calibrated on its report recode 790 of the 3,255 ordered inserts, none left refused and 23 left High (protected signal-peptide codons). A vendor re-quote of the order of record is with Becca |
| 6 | **RAET1G GPI call.** Curation calls RAET1G GPI-anchored; NetGPI calls all five forms non-anchored with confidence (0.98–0.99); the set keeps it as GPI. CFC1, CRIPTO and GML carry the same conflict. | All four sit on the N-after-SP rung, so placement does not change; the conflict matters for any analysis split by anchor type, and for RAET1G it is unresolved which annotation is right. | A literature or experimental check of RAET1G's anchor, recorded as a curation decision. | open (flagged as a conflict in the GPI set) |

## Tracker: every open item, its issue and who is on it

Status as of **2026-11-14** (rows not rechecked since 2026-11-05 say so). Each item is an issue on
[Deliverome-Project/l1-library-review](https://github.com/Deliverome-Project/l1-library-review/issues),
labeled by area; "Needs Becca" marks a decision rather than work. An owner is the working
session (and its pull request) carrying the item; "unassigned" means nobody has picked it up.
The six issues above map to #11–#12 (structure models), #4 (cargo), #5–#6 (loop gates),
#7 (GPI single-pass), #13 (vendor scan) and #21 (RAET1G).

### Order-blocking

| Issue | Owner | Status | Needs Becca |
|---|---|---|---|
| [#1](https://github.com/Deliverome-Project/l1-library-review/issues/1) Reserve barcodes for the 2,710 active L1 members | Becca | a block of 3,200 reserved provisionally in the shared registry (2026-11-10, [`reserve_l1_block.py`](../reserve_l1_block.py)): each passes both candidate 3′ flanks, the in situ prefix screen and an RNA-binding-protein motif screen; assignment to members waits for the frozen order. 444 prefix-distinct barcodes remain in the pool for other constructs |  |
| Approve the membrane-side gate's changes to 45 failing internal tag sites (29 drops, 14 replacements, 2 re-ranks; [PR 291](https://github.com/Deliverome-Project/deliverome-analysis/pull/291), `l1_internal_pick_changes.tsv`) | Becca | closed 2026-11-12: the gate is part of the shipped internal-site policy (brief section 5) | |
| [#2](https://github.com/Deliverome-Project/l1-library-review/issues/2) Assemble the tag cassette into every L1 insert, per rung | the rung sessions | closed 2026-11-14: all five rungs built; order of record `synthesis/2026-11-15-final-e/` (3,255 inserts of 2,673 members, cut at 541c6c18; supersedes `final-d`) |  |
| [#3](https://github.com/Deliverome-Project/l1-library-review/issues/3) Choose a concrete internal site for each of the 581 internal-rung members (583 since ENTPD2 joined) | tag-site pipeline review ([#224](https://github.com/Deliverome-Project/deliverome-analysis/pull/224)) | closed 2026-11-14: 531 loop sites, 28 snorkels, 22 held, SLC44A1 orientation unconfirmed, SLC68A1 moved to C-terminal (brief section 5) |  |
| [#4](https://github.com/Deliverome-Project/l1-library-review/issues/4) Decide the internal-rung cargo: ALFA alone (19 aa) or the 51-aa SpyTag003/ALFA tandem | Becca | closed: the full tandem at every rung |  |
| [#8](https://github.com/Deliverome-Project/l1-library-review/issues/8) Gate v4 picks on the member's own structure placement (9 cytoplasmic, 13 membrane-core picks) | `topology-loop-qc` ([#246](https://github.com/Deliverome-Project/deliverome-analysis/pull/246)) | closed 2026-11-12: the membrane-side gate restricts every internal site |  |
| [#11](https://github.com/Deliverome-Project/l1-library-review/issues/11) Fold the 196 active L1 members whose AlphaFold model is a different protein, and the 61 new ortholog members with no model | unassigned | closed 2026-11-14: every ordered insert is folded as its own tagged protein |  |
| [#13](https://github.com/Deliverome-Project/l1-library-review/issues/13) Run the vendor manufacturability scan on the L1 vectors | unassigned | closed 2026-11-13: Twist report and local recodes (727 inserts) |  |
| [#14](https://github.com/Deliverome-Project/l1-library-review/issues/14) Review and merge the N-terminal placement package (1,300 N-after-SP members) | `nterm-reviewer-package` ([#254](https://github.com/Deliverome-Project/deliverome-analysis/pull/254)) | merged ([#254](https://github.com/Deliverome-Project/deliverome-analysis/pull/254)) |  |

### Tag sites

| Issue | Owner | Status | Needs Becca |
|---|---|---|---|
| [#5](https://github.com/Deliverome-Project/l1-library-review/issues/5) Calibrate MIN_LOOP_LEN at the chosen internal cargo | unassigned | superseded by v5, the gate and refold |  |
| [#6](https://github.com/Deliverome-Project/l1-library-review/issues/6) Test the LOOP_SS secondary-structure gate | unassigned | superseded by v5, the gate and refold |  |
| [#7](https://github.com/Deliverome-Project/l1-library-review/issues/7) Give GPI single-pass members (BST2, RGMB) an ectodomain site | ranker-pu-segment | closed: BST2 internal, RGMB N after SP |  |
| [#9](https://github.com/Deliverome-Project/l1-library-review/issues/9) Make site tables carry residue identity and the numbered sequence's SHA-256 | unassigned | proposed (not rechecked since 2026-11-05) |  |
| [#10](https://github.com/Deliverome-Project/l1-library-review/issues/10) Confirm the 25 v4 picks on members whose designed protein is not the modeled one are fixed at c7be01a4 | `topology-loop-qc` ([#246](https://github.com/Deliverome-Project/deliverome-analysis/pull/246)) | to verify (not rechecked since 2026-11-05) |  |
| [#12](https://github.com/Deliverome-Project/l1-library-review/issues/12) Apply verified offsets for the 26 'contained' members | unassigned | open (not rechecked since 2026-11-05) |  |
| [#15](https://github.com/Deliverome-Project/l1-library-review/issues/15) Choose which S5E13 ranker ships: v4, the audit-corrected retrain, or a PU variant | `tier-retrain` ([#232](https://github.com/Deliverome-Project/deliverome-analysis/pull/232)), `ranker-pu-segment` ([#245](https://github.com/Deliverome-Project/deliverome-analysis/pull/245)) | closed 2026-11-12: v5 (RF with ligand-inclusive hard negatives) plus a ±2 ligand-contact veto |  |
| [#16](https://github.com/Deliverome-Project/l1-library-review/issues/16) Decide what role the agentic tag-site proposer plays for L1 multipass members | `agentic-tag-sites` ([#233](https://github.com/Deliverome-Project/deliverome-analysis/pull/233)) | benchmark reported (not rechecked since 2026-11-05) | yes |
| [#22](https://github.com/Deliverome-Project/l1-library-review/issues/22) Apply the tag-site literature audit to the live API and D1 | `tag-site-live-fix` ([surfaceome #459](https://github.com/Deliverome-Project/accessible-surfaceome/pull/459)) | done 2026-11-05 (D1 written; PR not yet merged) |  |
| [#26](https://github.com/Deliverome-Project/l1-library-review/issues/26) Give the disorder lane a secondary-structure gate | unassigned | open (not rechecked since 2026-11-05) |  |

### Library composition

| Issue | Owner | Status | Needs Becca |
|---|---|---|---|
| [#17](https://github.com/Deliverome-Project/l1-library-review/issues/17) Ortholog pool built from an April gene list: 418 eligible library genes were never drawable | Becca | adopted 2026-11-10: the Compara 116 CDS-gated draw over all library genes, with 11 correspondence replacements (`adopt_ortholog_cds_policy.py`) | |
| [#18](https://github.com/Deliverome-Project/l1-library-review/issues/18) Isoform universe misses 451 library genes with non-canonical UniProt isoforms | ([#257](https://github.com/Deliverome-Project/deliverome-analysis/pull/257), merged) | resolved: all 451 explained by the cohort rules (399 below the evidence bar, 52 incl. topology) |  |
| [#19](https://github.com/Deliverome-Project/l1-library-review/issues/19) Make the ortholog gene-to-UniProt step replayable (pin UniProt 2026_04) | ([#258](https://github.com/Deliverome-Project/deliverome-analysis/pull/258), merged) | replayed against UniProt 2026_04: 1,492 of 3,892 picks confirmed; the original search responses were never saved (permanent limit) |  |
| [#20](https://github.com/Deliverome-Project/l1-library-review/issues/20) Five library genes without an intact Ensembl 115 copy: build the UniProt allele or exclude | Becca | open (not rechecked since 2026-11-05) | yes |
| [#21](https://github.com/Deliverome-Project/l1-library-review/issues/21) Resolve the RAET1G (and CFC1, CRIPTO, GML) GPI curation conflict | unassigned | open (not rechecked since 2026-11-05) |  |

### Records and provenance

| Issue | Owner | Status | Needs Becca |
|---|---|---|---|
| [#23](https://github.com/Deliverome-Project/l1-library-review/issues/23) Re-export the viewer's gene list (fg-library.json) from the library | unassigned | open (not rechecked since 2026-11-05) |  |
| [#24](https://github.com/Deliverome-Project/l1-library-review/issues/24) Re-publish unstamped s4e9 S3 keys with provenance | unassigned | open (not rechecked since 2026-11-05) |  |
| [#25](https://github.com/Deliverome-Project/l1-library-review/issues/25) Publish the ColabFold structures to R2 with a D1 index row | unassigned | open (not rechecked since 2026-11-05) |  |
| [#27](https://github.com/Deliverome-Project/l1-library-review/issues/27) Fix the stale tag label in constructs_loop.tsv | `l1-recode-dups` ([#234](https://github.com/Deliverome-Project/deliverome-analysis/pull/234)) | minor (not rechecked since 2026-11-05) |  |
| [#28](https://github.com/Deliverome-Project/l1-library-review/issues/28) Verify the Tedman 2026 Nat Commun DOI cited in the s4e9 README | unassigned | minor (not rechecked since 2026-11-05) |  |

### Reviewer questions and review process

| Issue | Owner | Status | Needs Becca |
|---|---|---|---|
| [#29](https://github.com/Deliverome-Project/l1-library-review/issues/29) Reviewer question: default to MANE for the 111 MANE/APPRIS disagreements | Becca | posed in the brief (section 11, question 1) |  |
| [#30](https://github.com/Deliverome-Project/l1-library-review/issues/30) Reviewer question: N-terminal tagging at residue 1 for type III single-pass proteins (44 members) | Becca | posed in the brief (section 11, question 2) |  |
| [#31](https://github.com/Deliverome-Project/l1-library-review/issues/31) Give outside reviewers access to this repository | Becca | not started (not rechecked since 2026-11-05) | yes |
| [#32](https://github.com/Deliverome-Project/l1-library-review/issues/32) Merge order: #234, then #243, then #254 | Becca | closed: all merged |  |

### L2

| Issue | Owner | Status | Needs Becca |
|---|---|---|---|
| [#33](https://github.com/Deliverome-Project/l1-library-review/issues/33) L2: fold the 1,155 members without a structure model (~85 GPU-hours) | Becca | measured, not folded (not rechecked since 2026-11-05) | yes |


## Open before the order, moved here from the brief

The brief describes the design; everything still open before L1 is ordered lives here, which
is the contract the brief states. These three tables were in the brief's section 12 until
2026-11-16 and are unchanged apart from renumbered cross-references.

### Still with Becca

| Item | Why only she can do it |
|---|---|
| **The plate map changed.** Plates and wells now belong to inserts: 3,255 ordered inserts of 2,673 built members on 34 plates, numbered `l1-i0001`–`l1-i3255` by ordered insert length with a permanent arm id each; members are named by design id and the `l1-m` member id is retired. Against `2026-11-14-final-d` every plate position changes | it changes what the vendor and the bench receive |
| Whether to keep the Igκ swaps that rest on weaker second sources: 18 whose only confirmation of a cleaved signal peptide is DeepTMHMM v2 (library references ATRNL1, C5orf15, CDH22, CNNM1, CNTN5, GLP2R, GRIK1, LYPD6B, MMP16, SUSD4 and VSIG10L; orthologs of GABRR3, IL9R, LYPD2, LYPD6B and PVR; two SEMA4C isoforms), and 7 re-anchored on a UniProt SIGNAL feature that is itself a sequence-analysis prediction (ECO:0000255: CACHD1, CACNA2D4, CDH24, DCBLD2, GRM1, TNFRSF10A, UNC5C; companion brief section 11) | each replaces native sequence on the strength of a prediction |
| Confirm 4 C-terminal members kept against UniProt: LYSMD3, LYSMD4, SLC24A1 and TMEM222, where UniProt predicts a cytoplasmic C-terminus and DeepTMHMM v1 and v2 both call it extracellular (`l1_cterm_decisions.tsv`). Each is ordered as built until she decides | the call rests on predictors that disagree; no experimental topology exists |
| Sign off the thresholds set by design sessions that decide what is built: the refold bars (native-region lDDT 0.85, site lDDT 0.70, epitope relative accessibility 0.25, 5 packed epitope residues, a 10-residue cassette helix, type III TM1 15 residues; `tag_refold/score_tagged_folds.py`), the N-after-SP cleavage softening tolerance 0.07 and the 0.01 borderline margin (`nterm_decide.py`), the membrane-side gate's 14.0 Å pass line and 2.3 Å review band (`membrane-side-gate/calibration.json`), the ±2 ligand-contact veto (`contact_veto.py`), and the stalk sites chosen for failed C-terminal folds (KLRF1, HPN; `build_l1_cterm_arms.py`) | each was set by a design session under the delegation and changes which arm, or whether any arm, is ordered |
| May the licensed DeepTMHMM 1.0.24 and SignalP 6.0 packages (and the NetGPI 1.1 install record) be kept in a private S3 prefix, so the predictors themselves are saved and not only their hashes? | a licence question |
| Let CI check out full history (`fetch-depth: 0` in the test workflow) so the `pinned_refs` test runs there; on today's one-commit checkout it skips and says so | a CI configuration change |
| Ask the vendor to re-quote the recoded inserts (section 7) | external send; the local rules approximate the vendor's |
| Promoter, EFS or EF1a (question 7); the barcode 3′ flank (question 6) | construct design, not a protein call |
| Deploy the structure viewer with tag sites and folds; publish the staged 2026-11-14 bundle as a release of the review repository | production and shared records |
| Whether to run DeepTMHMM v1 on the 931 isoform forms the gate cannot verify | a compute call |
| Whether the 0.85 native-region lDDT bar suits small 4-TM folds. It was calibrated on 17 untagged controls refolded with the same settings (0.897–1.000; `LDDT_PASS`, `tag_refold/score_tagged_folds.py`); on the internal rung 22 members are held with no arm passing refold (`l1_internal_decisions.tsv`, `hold_no_passing_arm`) | a scoring threshold that decides whether members ship |
| Whether to keep counting a verdict just under a bar as a pass. The N-after-SP rung counts six arms as "borderline, counted as pass": FOLR1, SLITRK2 and SMO at site lDDT 0.691–0.698 against 0.7; SNORC, LYVE1 and CIST1 at native-region lDDT 0.842–0.847 against 0.85 (rule `lddt_borderline_counted_as_pass`, `reviewer-package/nterm/nterm_decisions.tsv`) | a scoring rule, applied on one rung only |
| Whether to accept buried SpyTag003 on multipass arms: 53 of the 72 multipass refold failures list a buried SpyTag003 among their reasons and are kept as built (`multipass_refold_dispositions.tsv`) | SpyCatcher capture depends on the epitope being exposed; kept-as-built is a judgment on the model, not a measurement |
| `l1_rbp_motif_cassettes.tsv` carries member ids with no design id or other stable key, so it cannot be joined after a renumber (`tests/test_member_id_policy.py` names it) | whether to regenerate it before the order |
| The literature pass on ligand contacts covers 437 of the 583 internal members; the other 146 were finalized on deposited and transferred ligands only (`l1_internal_no_literature_pass.tsv`) | whether to finish it before the order or accept the ranker's calls |

### Still running

| Item | Expected |
|---|---|
| Four built loop arms are flagged because each is the member's only loop arm: DISP1, GJB2 and KCNK4 (two members) | an agent and literature review before the order |
| 7 internal arms split a predicted N-glycosylation sequon | review before the order |
| 14 arms kept on a soft literature contact (flagged under the contact veto) | the literature search, finished |
| TMEM81, MXRA8 and IL1RAPL2 keep their native weak signal peptides, flagged because the two cleavage calls differ by more than 5 residues | Becca's call on the cut site |
| 23 N-after-SP inserts stay High by the local synthesis rules: a 50-bp window of 90–98% GC inside protected signal-peptide codons | the vendor's quote decides whether they are made as is |
| The agent review of every built site and the literature search for ligand contacts are unfinished | before the order |
| Membrane-side gate r9 | L2 |
| Native folds of the restarted cyno CHRM3, NIPAL4 and OR3A2; until they land, these three have no chosen structure model | next refold refresh |
| 73 AlphaFold DB models the design read are not yet in the `external://alphafold-db/v6/` archive (`tag_refold_afdb_unarchived.tsv`); the `structure_models` stage names them | l1-tag-refold |

### Known risks

| Risk | Consequence or check |
|---|---|
| The synthesis rules are a local approximation of the vendor's | a vendor re-quote of the recoded inserts is the real test |
| Tag-site model v5 was trained with features misnumbered for 4 of 127 proteins (858 of 35,650 rows) | retraining deferred to L2; the 11 L1 members hit by the same gap are redesigned on corrected features |
| A 10-residue spacer between the cassette and native residue 2 rescued none of the 12 multipass arms first classed fold-disturbed (`verdicts_multipass_spacer.tsv`) | they are kept as built and flagged (`multipass_refold_dispositions.tsv`) |
| 9 of 21 multipass taste receptors fail the refold, with the cassette next to TM1 on a 6–8-residue N-terminal tail | kept as built and flagged; a useful bench check |
| 14 class A peptide and chemokine GPCRs have the epitope packed against the receptor in the model | kept, with ALFA as the display control |
| AarI/PaqCI (1,459 of 3,255 ordered inserts), BtgZI (1,387) and SapI (all, inside ALFA's codons) sites remain (`l1_domestication_barcode_audit.json`) | no cloning step uses them; domesticate if one is adopted |

## Keeping this note current

When an issue closes, change its status to the date and the decision, and keep the row:
a reviewer who read an earlier version should be able to see what changed. If a fix changes
a library table, rerun [`brief_numbers.py`](brief_numbers.py), update any count above and in
the brief, rebuild the input bundle, and publish it as a new dated release of the review
repository.
