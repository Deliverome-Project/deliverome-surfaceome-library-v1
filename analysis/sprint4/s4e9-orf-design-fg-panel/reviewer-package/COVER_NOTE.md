# L1 design review: cover note on work in progress

The [design brief](BRIEF.md) describes the L1 library as it stands. This note tracks what is
still being worked on before L1 is ordered, separately from the brief, so the brief can be read
as the design and this note updated as items close. The first table explains the six design
issues reviewers are most likely to ask about; the tracker below it lists every open item as an
issue, with its owner. Status is as of **2026-11-05**. Comments on any of them are welcome; they
are not the questions we are asking reviewers to settle.

Terms follow the brief: a *member* is one construct, a *rung* its tag position, and counts
refer to tranche L1, active members, DeepTMHMM v1 topology. Every count is measured by
[`brief_numbers.py`](brief_numbers.py).

| # | Issue | Why it matters | What resolves it | Status |
|---:|---|---|---|---|
| 1 | **Structure models of a different protein.** 198 members have an AlphaFold model of another proteoform, and 26 a model of a longer one that contains theirs. | Tag-site coordinates would be numbered against a sequence the construct does not contain, and nothing downstream flags it. | Fold each "different" member from its designed CDS; verify and apply the offset for the 26 "contained" ones. | open |
| 2 | **Internal-rung cargo**: ALFA alone (19 aa) or the 51-aa SpyTag003/ALFA tandem. | At 19 aa the cassette fits 16-residue loops (e.g. the SLC16A9/12/14 block); the tandem does not, so this decides whether whole paralog groups have a site. | Decide whether internal members need SpyTag003, then calibrate the loop-length floor at that cargo size. | open |
| 3 | **Loop gates** `MIN_LOOP_LEN = 8` (minimum contiguous coil/turn run) and `LOOP_SS` (coil, turn, bend or 3₁₀ only). | They change the top-ranked site for 56% and 72% of internal members respectively. In a published-junction benchmark, 28% of large-cargo insertions sit in loops shorter than 8 (n = 18). | Wet data from a ladder of coil-run lengths; the benchmark corpus is 78% small epitope tags and cannot settle it. | open |
| 4 | **GPI single-pass members and the internal worklist.** The member-level ladder reroutes BST2 and RGMB to the internal rung, but the per-gene internal-site worklist and its site generator are topology-only and do not hold them. | A GPI single-pass protein needs an ectodomain site, not a terminal or inter-TM one; until the site lanes cover it, these members have a rung but no junction. | Build the GPI internal set from the anchored genes and add the GPI branch to the per-gene rung derivation. | open |
| 5 | **Vendor manufacturability scan** (GC windows, long repeats, homopolymers, hairpins) not yet run on L1. | Synthesis houses quote and schedule on these features; finding them before the order is cheaper than a failed region. | Run the scan over the L1 vectors and separate backbone regions already synthesized from new ones. | open |
| 6 | **RAET1G GPI call.** Curation calls RAET1G GPI-anchored; NetGPI calls all five forms non-anchored with confidence (0.98–0.99); the set keeps it as GPI. CFC1, CRIPTO and GML carry the same conflict. | All four sit on the N-after-SP rung, so placement does not change; the conflict matters for any analysis split by anchor type, and for RAET1G it is unresolved which annotation is right. | A literature or experimental check of RAET1G's anchor, recorded as a curation decision. | open |

## Tracker: every open item, its issue and who is on it

Status as of **2026-11-05**. Each item is an issue on
[Deliverome-Project/l1-library-review](https://github.com/Deliverome-Project/l1-library-review/issues),
labeled by area; "Needs Becca" marks a decision rather than work. An owner is the working
session (and its pull request) carrying the item; "unassigned" means nobody has picked it up.
The six issues above map to #11–#12 (structure models), #4 (cargo), #5–#6 (loop gates),
#7 (GPI single-pass), #13 (vendor scan) and #21 (RAET1G).

### Order-blocking

| Issue | Owner | Status | Needs Becca |
|---|---|---|---|
| [#1](https://github.com/Deliverome-Project/l1-library-review/issues/1) Reserve barcodes for the 2,735 active L1 members | unassigned | not started |  |
| [#2](https://github.com/Deliverome-Project/l1-library-review/issues/2) Assemble the tag cassette into every L1 insert, per rung | unassigned | blocked on site choices |  |
| [#3](https://github.com/Deliverome-Project/l1-library-review/issues/3) Choose a concrete internal site for each of the 567 internal-rung members | tag-site pipeline review ([#224](https://github.com/Deliverome-Project/deliverome-analysis/pull/224)) | in progress across the S5E13 PRs |  |
| [#4](https://github.com/Deliverome-Project/l1-library-review/issues/4) Decide the internal-rung cargo: ALFA alone (19 aa) or the 51-aa SpyTag003/ALFA tandem | Becca | open | yes |
| [#8](https://github.com/Deliverome-Project/l1-library-review/issues/8) Gate v4 picks on the member's own structure placement (9 cytoplasmic, 13 membrane-core picks) | `topology-loop-qc` ([#246](https://github.com/Deliverome-Project/deliverome-analysis/pull/246)) | finding reported; decision pending | yes |
| [#11](https://github.com/Deliverome-Project/l1-library-review/issues/11) Fold the 199 L1 members whose AlphaFold model is a different protein | unassigned | open |  |
| [#13](https://github.com/Deliverome-Project/l1-library-review/issues/13) Run the vendor manufacturability scan on the L1 vectors | unassigned | not run on L1 |  |
| [#14](https://github.com/Deliverome-Project/l1-library-review/issues/14) Review and merge the N-terminal placement package (1,300 N-after-SP members) | `nterm-reviewer-package` ([#254](https://github.com/Deliverome-Project/deliverome-analysis/pull/254)) | PR open (draft) |  |

### Tag sites

| Issue | Owner | Status | Needs Becca |
|---|---|---|---|
| [#5](https://github.com/Deliverome-Project/l1-library-review/issues/5) Calibrate MIN_LOOP_LEN at the chosen internal cargo | unassigned | blocked on the cargo decision |  |
| [#6](https://github.com/Deliverome-Project/l1-library-review/issues/6) Test the LOOP_SS secondary-structure gate | unassigned | open |  |
| [#7](https://github.com/Deliverome-Project/l1-library-review/issues/7) Give GPI single-pass members (BST2, RGMB) an ectodomain site | unassigned | open |  |
| [#9](https://github.com/Deliverome-Project/l1-library-review/issues/9) Make site tables carry residue identity and the numbered sequence's SHA-256 | unassigned | proposed |  |
| [#10](https://github.com/Deliverome-Project/l1-library-review/issues/10) Confirm the 25 v4 picks on members whose designed protein is not the modeled one are fixed at c7be01a4 | `topology-loop-qc` ([#246](https://github.com/Deliverome-Project/deliverome-analysis/pull/246)) | to verify |  |
| [#12](https://github.com/Deliverome-Project/l1-library-review/issues/12) Apply verified offsets for the 26 'contained' members | unassigned | open |  |
| [#15](https://github.com/Deliverome-Project/l1-library-review/issues/15) Choose which S5E13 ranker ships: v4, the audit-corrected retrain, or a PU variant | `tier-retrain` ([#232](https://github.com/Deliverome-Project/deliverome-analysis/pull/232)), `ranker-pu-segment` ([#245](https://github.com/Deliverome-Project/deliverome-analysis/pull/245)) | measurements done; decision pending | yes |
| [#16](https://github.com/Deliverome-Project/l1-library-review/issues/16) Decide what role the agentic tag-site proposer plays for L1 multipass members | `agentic-tag-sites` ([#233](https://github.com/Deliverome-Project/deliverome-analysis/pull/233)) | benchmark reported | yes |
| [#22](https://github.com/Deliverome-Project/l1-library-review/issues/22) Apply the tag-site literature audit to the live API and D1 | `tag-site-live-fix` ([surfaceome #459](https://github.com/Deliverome-Project/accessible-surfaceome/pull/459)) | PR open |  |
| [#26](https://github.com/Deliverome-Project/l1-library-review/issues/26) Give the disorder lane a secondary-structure gate | unassigned | open |  |

### Library composition

| Issue | Owner | Status | Needs Becca |
|---|---|---|---|
| [#17](https://github.com/Deliverome-Project/l1-library-review/issues/17) Ortholog pool built from an April gene list: 418 eligible library genes were never drawable | `l1-recode-dups` ([#234](https://github.com/Deliverome-Project/deliverome-analysis/pull/234)) | decision pending | yes |
| [#18](https://github.com/Deliverome-Project/l1-library-review/issues/18) Isoform universe misses 451 library genes with non-canonical UniProt isoforms | `l1-recode-dups` ([#234](https://github.com/Deliverome-Project/deliverome-analysis/pull/234)) | open |  |
| [#19](https://github.com/Deliverome-Project/l1-library-review/issues/19) Make the ortholog gene-to-UniProt step replayable (pin UniProt 2026_04) | `l1-recode-dups` ([#234](https://github.com/Deliverome-Project/deliverome-analysis/pull/234)) | in progress (UniProt download rate-limited) |  |
| [#20](https://github.com/Deliverome-Project/l1-library-review/issues/20) Five library genes without an intact Ensembl 115 copy: build the UniProt allele or exclude | Becca | open | yes |
| [#21](https://github.com/Deliverome-Project/l1-library-review/issues/21) Resolve the RAET1G (and CFC1, CRIPTO, GML) GPI curation conflict | unassigned | open |  |

### Records and provenance

| Issue | Owner | Status | Needs Becca |
|---|---|---|---|
| [#23](https://github.com/Deliverome-Project/l1-library-review/issues/23) Re-export the viewer's gene list (fg-library.json) from the library | unassigned | open |  |
| [#24](https://github.com/Deliverome-Project/l1-library-review/issues/24) Re-publish unstamped s4e9 S3 keys with provenance | unassigned | open |  |
| [#25](https://github.com/Deliverome-Project/l1-library-review/issues/25) Publish the ColabFold structures to R2 with a D1 index row | unassigned | open |  |
| [#27](https://github.com/Deliverome-Project/l1-library-review/issues/27) Fix the stale tag label in constructs_loop.tsv | `l1-recode-dups` ([#234](https://github.com/Deliverome-Project/deliverome-analysis/pull/234)) | minor |  |
| [#28](https://github.com/Deliverome-Project/l1-library-review/issues/28) Verify the Tedman 2026 Nat Commun DOI cited in the s4e9 README | unassigned | minor |  |

### Reviewer questions and review process

| Issue | Owner | Status | Needs Becca |
|---|---|---|---|
| [#29](https://github.com/Deliverome-Project/l1-library-review/issues/29) Reviewer question: default to MANE for the 111 MANE/APPRIS disagreements | Becca | posed in the brief (section 8, question 1) |  |
| [#30](https://github.com/Deliverome-Project/l1-library-review/issues/30) Reviewer question: N-terminal tagging at residue 1 for type III single-pass proteins (62 members) | Becca | posed in the brief (section 8, question 2) |  |
| [#31](https://github.com/Deliverome-Project/l1-library-review/issues/31) Give outside reviewers access to this repository | Becca | not started | yes |
| [#32](https://github.com/Deliverome-Project/l1-library-review/issues/32) Merge order: #234, then #243, then #254 | Becca | all three open as drafts | yes |

### L2

| Issue | Owner | Status | Needs Becca |
|---|---|---|---|
| [#33](https://github.com/Deliverome-Project/l1-library-review/issues/33) L2: fold the 1,155 members without a structure model (~85 GPU-hours) | Becca | measured, not folded | yes |

## Keeping this note current

When an issue closes, change its status to the date and the decision, and keep the row:
a reviewer who read an earlier version should be able to see what changed. If a fix changes
a library table, rerun [`brief_numbers.py`](brief_numbers.py), update any count above and in
the brief, rebuild the input bundle, and publish it as a new dated release of the review
repository.
