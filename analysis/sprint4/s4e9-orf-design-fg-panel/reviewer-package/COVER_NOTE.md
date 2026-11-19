# L1 design review: what is still open before the order

The [design brief](BRIEF.md) is the design. This note is everything the design does not yet
settle: the decisions still to be made, the work still running, and the risks we are carrying
knowingly. Read the decisions first: they are the ones that cannot move without you.

Status is as of **2026-11-18**, against the order of record `2026-11-15-final-e`. Terms follow
the brief: a *member* is one construct, a *rung* its tag position, an *arm* one tag site on a
member. Counts refer to tranche L1, active members, DeepTMHMM v1 topology, and every one is
measured by [`brief_numbers.py`](brief_numbers.py). Comments on anything below are welcome; they
are not the questions we are asking reviewers to settle.

## What reviewers ask about first

Six design issues account for most of what reviewers raise. Each appears once, here.

| Issue | Why it matters | Where it stands |
|---|---|---|
| **Structure models of a different protein.** 198 members have an AlphaFold model of another proteoform, and 26 a model of a longer one that contains theirs. | Tag-site coordinates would be numbered against a sequence the construct does not contain, and nothing downstream flags it. | Closed 2026-11-14. Every ordered insert is refolded from its own tagged protein, and an independent check matches all 3,255 folds to the order of record residue for residue (brief section 12) |
| **Internal-rung cargo**: ALFA alone (19 aa) or the 51-aa SpyTag003/ALFA tandem. | At 19 aa the cassette fits 16-residue loops (e.g. the SLC16A9/12/14 block); the tandem does not, so this decides whether whole paralog groups have a site. | closed: the full 51-aa cassette is used at every rung, internal included; each internal site is kept only if the tagged protein passes refold (brief section 5) |
| **Loop gates** `MIN_LOOP_LEN = 8` (minimum contiguous coil/turn run) and `LOOP_SS` (coil, turn, bend or 3₁₀ only). | They change the top-ranked site for 56% and 72% of internal members respectively. In a published-junction benchmark, 28% of large-cargo insertions sit in loops shorter than 8 (n = 18). | superseded 2026-11-12: internal sites come from tag-site model v5, the membrane-side gate and refold, not from these gates (brief section 5) |
| **GPI single-pass members and the internal worklist.** The member-level ladder reroutes BST2 and RGMB to the internal rung, but the per-gene internal-site worklist and its site generator are topology-only and do not hold them. | A GPI single-pass protein needs an ectodomain site, not a terminal or inter-TM one; until the site lanes cover it, these members have a rung but no junction. | closed 2026-11-12: BST2 is tagged in its ectodomain on the internal rung; RGMB, with no two-sided loop, is routed to N after SP |
| **Vendor manufacturability scan** (GC windows, long repeats, homopolymers, hairpins) not yet run on L1. | Synthesis houses quote and schedule on these features; finding them before the order is cheaper than a failed region. | closed 2026-11-13: Twist scored the first draft; local rules calibrated on its report recode 790 of the 3,255 ordered inserts, none left refused and 23 left High (protected signal-peptide codons). A vendor re-quote of the order of record is outstanding |
| **RAET1G GPI call.** Curation calls RAET1G GPI-anchored; NetGPI calls all five forms non-anchored with confidence (0.98–0.99); the set keeps it as GPI. CFC1, CRIPTO and GML carry the same conflict. | All four sit on the N-after-SP rung, so placement does not change; the conflict matters for any analysis split by anchor type, and for RAET1G it is unresolved which annotation is right. | Open, and flagged as a conflict in the GPI set. A literature or experimental check of RAET1G's anchor, recorded as a curation decision, would settle it |

## Decisions still open

Each of these changes what is built or what the vendor receives, and none can be settled
by a rule.

| Item | Why it takes a decision rather than a rule |
|---|---|
| **The plate map changed.** Plates and wells now belong to inserts: 3,255 ordered inserts of 2,673 built members on 34 plates, numbered `l1-i0001`–`l1-i3255` by ordered insert length with a permanent arm id each; members are named by design id and the `l1-m` member id is retired. Against `2026-11-14-final-d` every plate position changes | it changes what the vendor and the bench receive |
| Whether to keep the Igκ swaps that rest on weaker second sources: 18 whose only confirmation of a cleaved signal peptide is DeepTMHMM v2 (library references ATRNL1, C5orf15, CDH22, CNNM1, CNTN5, GLP2R, GRIK1, LYPD6B, MMP16, SUSD4 and VSIG10L; orthologs of GABRR3, IL9R, LYPD2, LYPD6B and PVR; two SEMA4C isoforms), and 7 re-anchored on a UniProt SIGNAL feature that is itself a sequence-analysis prediction (ECO:0000255: CACHD1, CACNA2D4, CDH24, DCBLD2, GRM1, TNFRSF10A, UNC5C; companion brief section 11) | each replaces native sequence on the strength of a prediction |
| Confirm 4 C-terminal members kept against UniProt: LYSMD3, LYSMD4, SLC24A1 and TMEM222, where UniProt predicts a cytoplasmic C-terminus and DeepTMHMM v1 and v2 both call it extracellular (`l1_cterm_decisions.tsv`). Each is ordered as built until the call is made | the call rests on predictors that disagree; no experimental topology exists |
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
| **Nine internal members have no arm left that carries a tag site.** An *arm* is one tag site on a member; a *loop arm* sits in an extracellular loop and a *snorkel* reaches the tag out past the membrane from elsewhere, so a snorkel alone is not an internal site. Review round r2 rejects an arm on each of 28 members; 19 keep a loop arm and nine do not — CNGA3 `l1-p0443`, GJB2 `l1-p0711`, KCNA3 `l1-p0988`, KCNC3 `l1-p0997`, KCNJ2 `l1-p1019`, KCNK4 `l1-p1030`, SCNN1A `l1-p1772` and `l1-p2497`, TRPM5 `l1-p2275` — of which three keep only a snorkel. Each has a ranker replacement residue proposed for at least one rejected arm (the r2 internal-rung review, section 2) | whether to take the ranker's replacement, find another site, or order these nine without an internal tag |
| Whether an annotated risk 3 to 5 residues from the junction is a concern, as a table rule. The rubric says it is and no table implements it: the adjacency window stops at 2 residues, and the combination step silences an agent verdict whose driver belongs to a table rule at any distance. **60 built loop arms carry an annotated N-glycan in that band and are published `makes_sense`**; a reviewer raised 52 of them and was silenced. `band_gap.py` measures it and changes nothing, and the same sweep finds no disulfide cysteine in the band | it would apply to every built arm, not only the ones a reviewer happened to mention, so it is a scientific call rather than a bug fix |
| Which tier governs a contact labeled `literature_small_ligand; structure`, which carries two at once. It matters because the 5-residue interface concern applies to a structure-tier contact and not to a literature one. 85 arms across 51 members have such a contact, and 45 of those arms across 28 members visibly reason about which way to read it | r2 reports that ruling on it once would settle more arms than any other open question on the internal rung |
| Four members sit on the skip list, which holds a member out of review while the ranker revises its decision rows, with no reason recorded anywhere: `l1-p1543`, `l1-p1732`, `l1-p1894`, `l1-p1967`, one built loop arm each. Of the 31 on the list, 23 still have a `replace` row and are correctly held | whether to unskip and review them or write the reason down |

## Still running

| Item | Expected |
|---|---|
| Four built loop arms were flagged because each is the member's only loop arm: DISP1, GJB2 and KCNK4 (two members). r2 reaches three of the four. DISP1 `l1-p0512-a2` is a concern on a co-receptor contact at distance 0 and keeps its arm; GJB2 `l1-p0711-a1` (functional site) and KCNK4 `l1-p1030-a1` (glycan, distance 1) are rejected and join the nine above. The second KCNK4 member, `l1-p2457`, is on the skip list and unreviewed | DISP1's co-receptor concern settled, the two rejects handled with the nine, and `l1-p2457` either unskipped or its hold recorded |
| 7 internal arms split a predicted N-glycosylation sequon | review before the order |
| 14 arms kept on a soft literature contact (flagged under the contact veto) | the literature search, finished |
| TMEM81, MXRA8 and IL1RAPL2 keep their native weak signal peptides, flagged because the two cleavage calls differ by more than 5 residues | a decision on the cut site |
| 23 N-after-SP inserts stay High by the local synthesis rules: a 50-bp window of 90–98% GC inside protected signal-peptide codons | the vendor's quote decides whether they are made as is |
| The internal rung's agent review is finished as of r2 (2026-11-16): 1,133 arms over 567 of 586 members, with 41 arms on 19 skip-list members unreviewed. The terminal rung's review and the literature search for ligand contacts are unfinished | before the order |
| AREG `l1-p0135`, both arms published `makes_sense`. The reviewer wrote a concern on each from recollection rather than dossier evidence and was overruled by the rule that makes a literature-tier observation a note. The claim is that insertion after M186 sits at the ADAM17 shedding site, which would put the tag on the shed ectodomain and be a reject under the cleavage-site rule. It is a rescue-round member with no r8 region, so no mechanical rule applies either way | the claim checked against a source |
| Both rule tables measure the distance from a UniProt feature correctly but print the residue number in UniProt coordinates, so a label is in the wrong frame whenever the member's proteoform is not the UniProt canonical sequence. No verdict is affected, the dossier now states the frame, and `label_frame_check.py` fails if a label stops reproducing its own distance. Two rejected arms carry such a span: SCNN1A `l1-p2497-a1` and KCNK10 `l1-p1024-a3` | the one-line fix in `annotation_adjacency.py` and `pore_rule.py`, and the tables regenerated against the EBI Proteins API |
| Membrane-side gate r9 | L2 |
| Native folds of the restarted cyno CHRM3, NIPAL4 and OR3A2; until they land, these three have no chosen structure model | next refold refresh |
| 73 AlphaFold DB models the design read are not yet in the `external://alphafold-db/v6/` archive (`tag_refold_afdb_unarchived.tsv`); the `structure_models` stage names them | l1-tag-refold |

## Known risks

Carried knowingly: each is accepted for L1 with the check that would settle it named.

| Risk | Consequence or check |
|---|---|
| The synthesis rules are a local approximation of the vendor's | a vendor re-quote of the recoded inserts is the real test |
| Tag-site model v5 was trained with features misnumbered for 4 of 127 proteins (858 of 35,650 rows) | retraining deferred to L2; the 11 L1 members hit by the same gap are redesigned on corrected features |
| A 10-residue spacer between the cassette and native residue 2 rescued none of the 12 multipass arms first classed fold-disturbed (`verdicts_multipass_spacer.tsv`) | they are kept as built and flagged (`multipass_refold_dispositions.tsv`) |
| 9 of 21 multipass taste receptors fail the refold, with the cassette next to TM1 on a 6–8-residue N-terminal tail | kept as built and flagged; a useful bench check |
| 14 class A peptide and chemokine GPCRs have the epitope packed against the receptor in the model | kept, with ALFA as the display control |
| AarI/PaqCI (1,459 of 3,255 ordered inserts), BtgZI (1,387) and SapI (all, inside ALFA's codons) sites remain (`l1_domestication_barcode_audit.json`) | no cloning step uses them; domesticate if one is adopted |

## Issue tracker

Every open item, by area. Each is an issue on
[Deliverome-Project/l1-library-review](https://github.com/Deliverome-Project/l1-library-review/issues);
a row that has not been rechecked since 2026-11-05 says so. "Needs a decision" marks a call to make
rather than work to do, and an owner is the working session carrying the item, with its pull request,
or "unassigned" if nobody has picked it up. The six issues in "What reviewers ask about first"
are tracked here as issues 11 and 12 (structure models), 4 (cargo), 5 and 6 (loop gates),
7 (GPI single-pass), 13 (vendor scan) and 21 (RAET1G).

### Order-blocking
| Issue | Owner | Status | Needs a decision |
|---|---|---|---|
| [#1](https://github.com/Deliverome-Project/l1-library-review/issues/1) Reserve barcodes for the 2,710 active L1 members | Becca | a block of 3,200 reserved provisionally in the shared registry (2026-11-10, [`reserve_l1_block.py`](../reserve_l1_block.py)): each passes both candidate 3′ flanks, the in situ prefix screen and an RNA-binding-protein motif screen; assignment to members waits for the frozen order. 444 prefix-distinct barcodes remain in the pool for other constructs |  |
| Approve the membrane-side gate's changes to 45 failing internal tag sites (29 drops, 14 replacements, 2 re-ranks; [PR 291](https://github.com/Deliverome-Project/deliverome-analysis/pull/291), `l1_internal_pick_changes.tsv`) | Becca | closed 2026-11-12: the gate is part of the shipped internal-site policy (brief section 5) | |
| [#2](https://github.com/Deliverome-Project/l1-library-review/issues/2) Assemble the tag cassette into every L1 insert, per rung | the rung sessions | closed 2026-11-14: all five rungs built; order of record `synthesis/2026-11-15-final-e/` (3,255 inserts of 2,673 members, cut at 541c6c18; supersedes `final-d`) |  |
| [#3](https://github.com/Deliverome-Project/l1-library-review/issues/3) Choose a concrete internal site for each of the 581 internal-rung members (583 since ENTPD2 joined) | tag-site pipeline review ([#224](https://github.com/Deliverome-Project/deliverome-analysis/pull/224)) | closed 2026-11-14: 531 loop sites, 28 snorkels, 22 held, SLC44A1 orientation unconfirmed, SLC68A1 moved to C-terminal (brief section 5) |  |
| [#4](https://github.com/Deliverome-Project/l1-library-review/issues/4) Decide the internal-rung cargo: ALFA alone (19 aa) or the 51-aa SpyTag003/ALFA tandem | Becca | closed: the full tandem at every rung |  |
| [#8](https://github.com/Deliverome-Project/l1-library-review/issues/8) Gate v4 picks on the member's own structure placement (9 cytoplasmic, 13 membrane-core picks) | `topology-loop-qc` ([#246](https://github.com/Deliverome-Project/deliverome-analysis/pull/246)) | closed 2026-11-12: the membrane-side gate restricts every internal site |  |
| [#11](https://github.com/Deliverome-Project/l1-library-review/issues/11) Fold the 196 active L1 members whose AlphaFold model is a different protein, and the 61 new ortholog members with no model | unassigned | closed 2026-11-14: every ordered insert is folded as its own tagged protein |  |
| [#13](https://github.com/Deliverome-Project/l1-library-review/issues/13) Run the vendor manufacturability scan on the L1 vectors | unassigned | closed 2026-11-13: Twist report and local recodes (790 of the 3,255 ordered inserts, 23 left High; `synthesis_compliance_summary.json` against `2026-11-15-final-e`) |  |
| [#14](https://github.com/Deliverome-Project/l1-library-review/issues/14) Review and merge the N-terminal placement package (1,300 N-after-SP members) | `nterm-reviewer-package` ([#254](https://github.com/Deliverome-Project/deliverome-analysis/pull/254)) | merged ([#254](https://github.com/Deliverome-Project/deliverome-analysis/pull/254)) |  |

### Tag sites
| Issue | Owner | Status | Needs a decision |
|---|---|---|---|
| [#5](https://github.com/Deliverome-Project/l1-library-review/issues/5) Calibrate MIN_LOOP_LEN at the chosen internal cargo | unassigned | superseded by v5, the gate and refold |  |
| [#6](https://github.com/Deliverome-Project/l1-library-review/issues/6) Test the LOOP_SS secondary-structure gate | unassigned | superseded by v5, the gate and refold |  |
| [#7](https://github.com/Deliverome-Project/l1-library-review/issues/7) Give GPI single-pass members (BST2, RGMB) an ectodomain site | ranker-pu-segment | closed: BST2 internal, RGMB N after SP |  |
| [#9](https://github.com/Deliverome-Project/l1-library-review/issues/9) Make site tables carry residue identity and the numbered sequence's SHA-256 | unassigned | proposed (not rechecked since 2026-11-05) |  |
| [#10](https://github.com/Deliverome-Project/l1-library-review/issues/10) Confirm the 25 v4 picks on members whose designed protein is not the modeled one are fixed at c7be01a4 | `topology-loop-qc` ([#246](https://github.com/Deliverome-Project/deliverome-analysis/pull/246)) | to verify (not rechecked since 2026-11-05) |  |
| [#12](https://github.com/Deliverome-Project/l1-library-review/issues/12) Apply verified offsets for the 26 'contained' members | unassigned | open (not rechecked since 2026-11-05) |  |
| [#15](https://github.com/Deliverome-Project/l1-library-review/issues/15) Choose which S5E13 ranker ships: v4, the audit-corrected retrain, or a PU variant | `tier-retrain` ([#232](https://github.com/Deliverome-Project/deliverome-analysis/pull/232)), `ranker-pu-segment` ([#245](https://github.com/Deliverome-Project/deliverome-analysis/pull/245)) | closed 2026-11-12: v5 (RF with ligand-inclusive hard negatives) plus a ±2 ligand-contact veto |  |
| [#16](https://github.com/Deliverome-Project/l1-library-review/issues/16) Decide what role the agentic tag-site proposer plays for L1 multipass members | `agentic-tag-sites` ([#233](https://github.com/Deliverome-Project/deliverome-analysis/pull/233)) | r2 2026-11-16: the proposer reviews, it does not place. Every built internal arm is reviewed once against a published rubric and the verdict is advisory — 1,133 arms over 567 of 586 members, 780 makes-sense, 321 concern, 32 reject; 275 of the 353 verdicts above makes-sense are set by a table rule and 78 by the agent. Three decisions it raises are under "Decisions still open" | yes |
| [#22](https://github.com/Deliverome-Project/l1-library-review/issues/22) Apply the tag-site literature audit to the live API and D1 | `tag-site-live-fix` ([surfaceome #459](https://github.com/Deliverome-Project/accessible-surfaceome/pull/459)) | done 2026-11-05 (D1 written; PR not yet merged) |  |
| [#26](https://github.com/Deliverome-Project/l1-library-review/issues/26) Give the disorder lane a secondary-structure gate | unassigned | open (not rechecked since 2026-11-05) |  |

### Library composition
| Issue | Owner | Status | Needs a decision |
|---|---|---|---|
| [#17](https://github.com/Deliverome-Project/l1-library-review/issues/17) Ortholog pool built from an April gene list: 418 eligible library genes were never drawable | Becca | adopted 2026-11-10: the Compara 116 CDS-gated draw over all library genes, with 11 correspondence replacements (`adopt_ortholog_cds_policy.py`) | |
| [#18](https://github.com/Deliverome-Project/l1-library-review/issues/18) Isoform universe misses 451 library genes with non-canonical UniProt isoforms | ([#257](https://github.com/Deliverome-Project/deliverome-analysis/pull/257), merged) | resolved: all 451 explained by the cohort rules (399 below the evidence bar, 52 incl. topology) |  |
| [#19](https://github.com/Deliverome-Project/l1-library-review/issues/19) Make the ortholog gene-to-UniProt step replayable (pin UniProt 2026_04) | ([#258](https://github.com/Deliverome-Project/deliverome-analysis/pull/258), merged) | replayed against UniProt 2026_04: 1,492 of 3,892 picks confirmed; the original search responses were never saved (permanent limit) |  |
| [#20](https://github.com/Deliverome-Project/l1-library-review/issues/20) Five library genes without an intact Ensembl 115 copy: build the UniProt allele or exclude | Becca | open (not rechecked since 2026-11-05) | yes |
| [#21](https://github.com/Deliverome-Project/l1-library-review/issues/21) Resolve the RAET1G (and CFC1, CRIPTO, GML) GPI curation conflict | unassigned | open (not rechecked since 2026-11-05) |  |

### Records and provenance
| Issue | Owner | Status | Needs a decision |
|---|---|---|---|
| [#23](https://github.com/Deliverome-Project/l1-library-review/issues/23) Re-export the viewer's gene list (fg-library.json) from the library | unassigned | open (not rechecked since 2026-11-05) |  |
| [#24](https://github.com/Deliverome-Project/l1-library-review/issues/24) Re-publish unstamped s4e9 S3 keys with provenance | unassigned | open (not rechecked since 2026-11-05) |  |
| [#25](https://github.com/Deliverome-Project/l1-library-review/issues/25) Publish the ColabFold structures to R2 with a D1 index row | unassigned | open (not rechecked since 2026-11-05) |  |
| [#27](https://github.com/Deliverome-Project/l1-library-review/issues/27) Fix the stale tag label in constructs_loop.tsv | `l1-recode-dups` ([#234](https://github.com/Deliverome-Project/deliverome-analysis/pull/234)) | minor (not rechecked since 2026-11-05) |  |
| [#28](https://github.com/Deliverome-Project/l1-library-review/issues/28) Verify the Tedman 2026 Nat Commun DOI cited in the s4e9 README | unassigned | minor (not rechecked since 2026-11-05) |  |

### Reviewer questions and review process
| Issue | Owner | Status | Needs a decision |
|---|---|---|---|
| [#29](https://github.com/Deliverome-Project/l1-library-review/issues/29) Reviewer question: default to MANE for the 111 MANE/APPRIS disagreements | Becca | posed in the brief (section 11, question 1) |  |
| [#30](https://github.com/Deliverome-Project/l1-library-review/issues/30) Reviewer question: N-terminal tagging at residue 1 for type III single-pass proteins (44 members) | Becca | posed in the brief (section 11, question 2) |  |
| [#31](https://github.com/Deliverome-Project/l1-library-review/issues/31) Give outside reviewers access to this repository | Becca | not started (not rechecked since 2026-11-05) | yes |
| [#32](https://github.com/Deliverome-Project/l1-library-review/issues/32) Merge order: #234, then #243, then #254 | Becca | closed: all merged |  |

### L2
| Issue | Owner | Status | Needs a decision |
|---|---|---|---|
| [#33](https://github.com/Deliverome-Project/l1-library-review/issues/33) L2: fold the 1,155 members without a structure model (~85 GPU-hours) | Becca | measured, not folded (not rechecked since 2026-11-05) | yes |

## Decision log

| Date | Was | Now | Why |
|---|---|---|---|
| 2026-11-16 | Three decision tables sat in the brief's section 12 | Moved here, cross-references renumbered | The brief states the design; what is unsettled belongs in this note |
| 2026-11-16 | Issue 1 counted 2,629 refolded inserts | 3,255 | 2,629 was the superseded `final-c` cut; recounted against `2026-11-15-final-e` |
| 2026-11-18 | Issue 13 recorded 727 recoded inserts | 790, 23 left High | 727 predated `final-e`; `synthesis_compliance_summary.json` measures 790, which the headline table already said |
| 2026-11-18 | The note predated the internal rung's r2 review | r2 folded in | r2 landed 2026-11-16 at 22:00, nine hours after this note was last rendered |
| 2026-11-18 | Decisions sat below the issue tracker | Decisions first, tracker last | The reader of this note is the person who has to decide |

## Keeping this note current

When an issue closes, change its status to the date and the decision, and keep the row:
a reviewer who read an earlier version should be able to see what changed. If a fix changes
a library table, rerun [`brief_numbers.py`](brief_numbers.py), update any count above and in
the brief, rebuild the input bundle, and publish it as a new dated release of the review
repository.
