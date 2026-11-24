# Lane check: is "N at residue 1, multipass" right for each member on it?

**Scope.** The 561 L1 members that [`l1_tag_paths.tsv`](l1_tag_paths.tsv) places on the rung
"N at residue 1, multipass", where the tag cassette goes directly after Met1, on the library as
redrawn by [PR #283](https://github.com/Deliverome-Project/deliverome-analysis/pull/283)'s second
pass (isoform gate, Compara 116 ortholog panel, reference fixes and the lane overrides applied so
far). The rung, its GPCR and non-GPCR
composition and the leader question are described in
[NTERM_RES1_MULTIPASS.md](NTERM_RES1_MULTIPASS.md). This note asks a prior question: for each
member, is this the right lane at all?

**Identifiers.** Members are named by `design_id` (`l1-pNNNN`), which is stable across every
renumbering, and by `candidate_id` (`l1-mNNNN`) in #283's numbering. `design_id` is the first
column of every table and the key the review is joined on. Member ids quoted inside the review's
evidence text are as written on 2026-11-11, before #283's renumbering.

**Result.** For 560 members the lane holds; 1 stays unresolved (TSPO2). Every move or
replacement this check recommended is already applied on #283, so those members have left the rung
or been restarted (see "What PR #283's redraw settled" and "Decisions taken since"). Every verdict is a recommendation until Becca confirms it,
and none changes `l1_tag_paths.tsv`, the manifest or any placement record.

## The test

A member belongs on this lane only if its ordered protein has, at residue 1:

| Condition | Why the lane needs it |
|---|---|
| 1. An extracellular N-terminus | The cassette after Met1 is displayed only if residue 1 faces out. |
| 2. No cleavable signal peptide | A signal peptide would be cut off with the cassette, or carry it into the ER lumen and away; such members belong on "N after SP". |
| 3. At least two TM helices | With one helix the member is type III and goes on the single-pass lane, where that helix alone decides orientation. |

**DeepTMHMM v1 is the call of record** (Becca, 2026-11-11). The lane comes from the member's own
`deeptmhmm_v1` string, matched on the sha256 of its ordered sequence, and every column in the
output that carries the call says `deeptmhmm_v1`. Four other sources are read only as
cross-checks. A cross-check that disagrees with v1 on any condition sends the member to
literature review; it never replaces the v1 call.

| Cross-check | Conditions | When it is read |
|---|---|---|
| UniProt topological domains and SIGNAL feature, with evidence codes | 1, 2 (3 only for a member identical to the canonical) | Only when the member's first 20 residues are the canonical's, because the canonical's annotation says nothing about a different N-terminus |
| DeepTMHMM v2 on the member's own sequence | 1, 2, 3 | Whenever a v2 row exists (558 of 561; the three cyno orthologs restarted on #283 have none yet) |
| SignalP 6 | 2 | When an existing SignalP 6 run saw the member's first 70 residues (505 of 561). SignalP truncates every input to 70 residues, so a run on the same 70 residues is the member's own call. |
| The member's chosen AlphaFold DB or ColabFold model | 1, 3 | Whenever a model is chosen and published (549 of 561 fit a membrane) |

Three pattern checks also send a member to review: residue 1 is not Met; the member begins at or
after the canonical's first TM helix; and a GPCR of class B, C, F or adhesion, which normally
carry a signal peptide, has none.

### Reading v2 and the structure model

v2 writes its own state alphabet. Every v2 row in the pinned topology snapshot also stores v2
projected onto v1's letters, so the script derives the letter-by-letter mapping from all 33,015
v2 rows and refuses to run unless each letter maps to one state. Three v2 letters map to v1's
helix; only one of them reproduces the snapshot's own TM-helix count (100% of rows), and the
other two are reentrant and interfacial segments such as GPCR helix 8, which do not cross the
membrane.

The structure check fits a planar membrane to the model: the normal is the mean axis of the
helices that v1 or v2 places (their positions only, never their sides), from C-alpha atoms at
pLDDT ≥ 70. On that plane it counts **crossings**, the number of times the chain passes from one
face (more than 10 Å from the mid-plane, well inside the ~30 Å hydrocarbon core's edges) to the
other, ignoring residues below pLDDT 50, which are often modelled as extended ribbons through the
slab. It calls the **cytoplasmic face** by the positive-inside rule: the face with more Lys + Arg
within 10 residues of a membrane-core residue, if the difference is at least 3. It made that call
for 439 members; for the other 110 the flank charges were too even to call a side.

## Results

| Verdict | GPCR | Not a GPCR | Total |
|---|---:|---:|---:|
| Lane holds | 501 | 59 | **560** |
| Recommended elsewhere | 0 | 0 | **0** |
| Unresolved | 0 | 1 | **1** |

Unit: library member. 535 members had no cross-check disagreeing with v1; 26 did, and every one
has a literature review row. 36 members on the rung carry a review row: the 26, plus members the
review keeps on the lane that no cross-check now flags (an earlier UniProt parse error flagged six
TAS2Rs, SLC39A2 and POPDC1, see the decision log; OR8G3 and cyno NIPAL4, whose starts were recoded
or restarted).

| Lane holds, by confidence | Members | What it rests on |
|---|---:|---|
| high | 476 | At least three cross-checks agree with v1, one of them experimental: a UniProt ECO:0000269 topology, or the Tedman screen's surface display of an HA tag right after Met1 ([Tedman et al., 2026](https://pubmed.ncbi.nlm.nih.gov/42736269/)); or the review found structural or epitope evidence |
| medium | 81 | Two cross-checks agree, or the review rests on homology |
| low | 3 | Only one cross-check available: the adopted cyno OR51C1 and mouse POPDC3 orthologs (`l1-p2964`, `l1-p3098`; v2 only) and the restarted cyno OR3A2 (`l1-p2963`; the Tedman screen of the human receptor only) |

The Tedman intensities are gene-level. Two PTGDR members carry the same value, so a screen value
supports the orientation of a gene's canonical receptor, not the display of a member's isoform
or ortholog.

### Unresolved member

| Design id | Member | Gene | Set | Best call | What would settle it |
|---|---|---|---|---|---|
| `l1-p2302` | `l1-m0074` | TSPO2 | reference | stay | No TSPO2 topology data; TSPO1 evidence comes from mitochondrial membranes, and v1 and v2 disagree. Surface staining of the tagged member, or a C-terminal fallback if N-in. |

### What PR #283's redraw settled

Fourteen reviewed members are no longer on the rung, so their verdicts need no decision. Their
review rows are kept; retired ones are marked `retired` with #283's reason.

| Design ids | Members | What #283 did |
|---|---|---|
| `l1-p2400`, `l1-p2444`, `l1-p2486`, `l1-p2491`, `l1-p2516`, `l1-p2730`, `l1-p2477` | BEST1, GPM6B, PLPP2, PTCH1, SLC9D1, GPR161 and OPN3 isoforms | retired (isoform gate) |
| `l1-p2558`, `l1-p2606` | cyno CRHR2 and SLC22A17 | retired (Compara 116 ortholog policy) |
| `l1-p0067` | ADGRF3 | reference fixed to the signal-peptide isoform, which moves to "N after SP" |
| `l1-p1213`, `l1-p2466`, `l1-p0128` | MFSD2A ×2, AQP12B | moved to internal (lane overrides) |
| `l1-p0724` | GLP2R | moved to "N after SP" (lane override) |

#283 also changed the start of two reviewed members, and each review was redone on the new
sequence:

| Design id | Gene | Change on #283 | Verdict now |
|---|---|---|---|
| `l1-p1525` | OR8G3 | GTG start recoded to ATG | Lane holds, medium |
| `l1-p2962` | NIPAL4 (cyno) | Restarted at cyno Met215, the human NIPAL4 Met1 equivalent (404 aa); this removes the 214-residue Arg/Lys- and Gly-rich extension for which this check had recommended replacing the member | Lane holds, medium: v1 calls no signal peptide, a 55-residue extracellular tail and 9 helices, matching the human canonical's UniProt topology (residues 1–55 extracellular, by similarity); v2, SignalP 6 and a model do not yet cover the new sequence |

The restart also gave cyno CHRM3 and OR3A2 (`l1-p2938`, `l1-p2963`) new sequences with no v2
row or model yet; v1 still holds for both, at lower confidence.

### Decisions taken since

On 2026-11-12 Becca delegated the per-member calls to the multipass session
([NTERM_RES1_MULTIPASS.md](NTERM_RES1_MULTIPASS.md), "Decisions taken under Becca's delegation"),
and the library writer applies lane changes through `lane_overrides.tsv` on #283. Their calls
now agree with this check for every member:

| Member | Call | Agrees with this check |
|---|---|---|
| MFSD2A ×2, AQP12B, GLP2R | moved (applied on #283) | yes |
| OR8G3 | recoded to ATG, stays on this lane | yes |
| cyno NIPAL4 (`l1-p2962`) | restarted at Met215 (applied on #283) | yes |
| OPN3 isoform | retired | yes |
| TSPO2 (`l1-p2302`) | stays on this lane, Igκ + tag, low confidence | yes (unresolved here) |
| SLC61A1 (`l1-p1994`, `l1-p3105`) | stays on this lane, Igκ + tag | yes |

SLC61A1 stays on this lane on the topology evidence. Every cross-check is on v1's side: v2 calls
a 19-residue TM1 at residues 4–22 with residue 1 outside and 13 helices; the structure model
(pLDDT 89) crosses 13 times with residue 1 outside by flank K/R 19 against 6; SignalP 6 calls no
signal peptide. The rejected alternative was to move it to an internal site: UniProt (Q6N075)
annotates 12 helices with TM1 at residues 1–21, the major facilitator fold is usually N-in/C-in,
and an internal site chosen under the membrane-side gate would display either way. The writer
(l1-reviewer-package) kept it on this lane on 2026-11-13, under Becca's delegation, and
l1-integration confirmed.

### What the cross-checks resolved without a change

| Member | Gene | Flag | Resolution |
|---|---|---|---|
| `l1-p0475` | CRHR2 | UniProt SIGNAL 1–19; class B without a signal peptide | An uncleaved pseudo signal peptide ([Rutz et al., 2006](https://pubmed.ncbi.nlm.nih.gov/16766521/)); HA after Met1 displays at 2,385 a.u. |
| `l1-p0448`, `l1-p3071`, `l1-p0749`, `l1-p0758`, `l1-p0882`, `l1-p3123` | CNR1 (human, mouse), GPR135, GPR152, HTR7 (reference, isoform) | Model counts 8 crossings | A long, low-confidence N-tail passing through the fitted slab; CNR1 has a crystal structure ([Hua et al., 2016](https://pubmed.ncbi.nlm.nih.gov/27768894/)). |
| six TAS2Rs, OR4M2B, SLC39A2 | | UniProt N-terminus cytoplasmic | A parse error for the TAS2Rs and SLC39A2 (fixed); OR4M2B's unreviewed entry is inverted relative to every other olfactory receptor. TAS2R46 has a cryo-EM structure with an extracellular N-terminus ([Xu et al., 2022](https://pubmed.ncbi.nlm.nih.gov/36108005/)). |
| `l1-p2054`, `l1-p2518`, `l1-p2144` | SMIM7, TLCD1 | UniProt SIGNAL (sequence analysis) | SignalP 6 calls none; the review reads both as uncleaved N-out signal anchors. |
| `l1-p2946` | GPRC5A (cyno) | Class C GPCR without a signal peptide | The GPRC5 family lacks the class C ectodomain and signal peptide; the human receptor displays HA after Met1 at 6,447 a.u. |
| `l1-p2492` | PTGDR isoform | GPCR with 6 helices, C-terminus outside | Residue 1 is outside, so the lane is right for this sequence; whether a form missing helix 7 folds is a separate question about the isoform set. |

The per-member evidence and PMIDs are in [`nterm_res1_lane_review.tsv`](nterm_res1_lane_review.tsv).
Every PMID there was resolved in PubMed by the reviewer that cited it, and the ones cited on this
page were checked again against their PubMed records.

### Compared with the earlier non-GPCR review

The leader review in NTERM_RES1_MULTIPASS.md held three non-GPCRs. This check settles all three:
SLC61A1 (`l1-p1994`) passes every cross-check and stays on the lane; SLC9D1's Met-less isoform is retired by #283; and AQP12B moves to
internal. MFSD2A is off the rung, as before.

## Limits

- The cross-checks are not independent of each other in every case: v2 shares training lineage
  with v1, and the structure check fits its membrane on v1 and v2 helix positions.
- The structure check's thresholds (10 Å sides, pLDDT 50 and 70, a K/R margin of 3) are analysis
  choices, not calibrated against experimental topology here. The membrane-side calibration on
  PR [#291](https://github.com/Deliverome-Project/deliverome-analysis/pull/291) is the place to
  calibrate them once it merges; it uses the same helix-axis construction.
- 56 members (35 orthologs, 18 library references, 3 isoforms) have no SignalP 6 call: no stored
  run saw their first 70 residues. v2 still covers condition 2 for them.
- 12 members have no structure check: the three cyno orthologs restarted on #283 (CHRM3, NIPAL4,
  OR3A2) have no chosen model yet; eight adopted orthologs have one listed but not yet published in
  S3; and the SMIM7 isoform (`l1-p2518`) has too few confident helices to fit a membrane.
- A member no cross-check disagreed with passes on v1 plus agreement, not on literature. For the
  three low-confidence members only one cross-check was available.

## Methods

[`nterm_res1_lane_check.py`](nterm_res1_lane_check.py) takes the rung's members from
`l1_tag_paths.tsv` and the manifest (`candidate_manifest_tranches.tsv`), never from a list, so
rerunning it after a redraw checks the new members. It stops if a member's designed sequence no longer hashes to the value in
`l1_tag_paths.tsv`, if its v1 topology version differs, or if a reviewed member's ordered
sequence has changed since its review (each review row records the sequence hash it was read
on). The SignalP snapshot is joined by `design_id` and used only while the member's first 70
residues still hash to the value the run saw. `--rekey-review` brings the review table's `candidate_id` up to the manifest's numbering by
`design_id` and marks retired members with their reason. It reads GPCR status and Tedman
values from `nterm_res1_multipass.tsv`
([PR #280](https://github.com/Deliverome-Project/deliverome-analysis/pull/280)), and stops if
that table does not cover exactly the rung.

**Every input is pinned** (Becca, 2026-11-11), so a rerun is needed only when one moves. The
script runs offline from the snapshots committed beside it, and each row of
`nterm_res1_lane_check.tsv` carries the fingerprints of what it rests on:

| Input | Snapshot | Fingerprint columns |
|---|---|---|
| UniProt entries (517 accessions, release 2026_04) | [`nterm_res1_lane_uniprot.tsv`](nterm_res1_lane_uniprot.tsv) | `uniprot_release`, `uniprot_entry_version`, `uniprot_sequence_version`, `uniprot_sequence_sha256` |
| Topology | `topology.SNAPSHOT_KEY`, the pinned parquet | `topology_snapshot_key`, `deeptmhmm_v1_topology_version`, `deeptmhmm_v1_tool_version`, `deeptmhmm_v2_topology_version`, `deeptmhmm_v2_tool_version`, `sequence_sha256` |
| SignalP 6.0 | [`nterm_res1_signalp6.tsv`](nterm_res1_signalp6.tsv), from run `sp6_2026_09_28` on the `signalp6-runs` Modal volume | `signalp6_version`, `signalp6_source` |
| Structure models | `library_structure_payload.tsv`, the chosen model per member (AlphaFold DB v6 or a ColabFold release) | `model_version`, `model_ref`, `model_sha256` |

`--check-drift` re-fetches only UniProt's per-entry versions, recomputes each member's sequence
hash and v1 topology version, lists the members whose UniProt entry or sequence version,
sequence hash, topology version or rung membership changed, and exits nonzero if any did.

    uv run --locked python nterm_res1_lane_check.py                    # offline, from the snapshots
    uv run --locked python nterm_res1_lane_check.py --check-drift      # which members to rerun
    uv run --locked python nterm_res1_lane_check.py --refresh-uniprot  # re-pin UniProt, then rerun

Outputs: [`nterm_res1_lane_screen.tsv`](nterm_res1_lane_screen.tsv) (every cross-check value per
member) and [`nterm_res1_lane_check.tsv`](nterm_res1_lane_check.tsv) (one verdict per member,
keyed by `design_id` and `candidate_id`: `lane_ok`, `recommended_lane`, `confidence`, evidence, PMIDs, `decided_by`, `implemented_by`,
`reason`, and the fingerprints). Both are also at
`processed://sprint4/s4e9-orf-design-fg-panel/nterm-res1-lane-check/`.

## Decision log

| Date | Decision | Reason | Decided by / implemented by |
|---|---|---|---|
| 2026-11-11 | Check every member on the rung, not only those flagged by the leader review | The leader review found v1 wrong for MFSD2A and several isoforms | Becca / Claude Code session l1-nterm-lane-check |
| 2026-11-11 | DeepTMHMM v1 is the call of record; UniProt, v2, SignalP 6 and the structure model are cross-checks that only flag a disagreement and send the member to review | v1 is the repo's topology of record; replacing it member by member would make the lane depend on which source happened to be consulted | Becca / Claude Code session l1-nterm-lane-check |
| 2026-11-11 | Every input pinned with its version, fingerprints in every row, and a drift check | To know which members to rerun when UniProt or the topology snapshot moves | Becca / Claude Code session l1-nterm-lane-check |
| 2026-11-11 | Cross-check rules and thresholds as in the script's constants | Each constant's comment gives its reason | no decision taken: analysis choice, not yet reviewed by Becca / Claude Code session l1-nterm-lane-check |
| 2026-11-11 | UniProt parse fixed to read single-residue domains (`TOPO_DOM 1;`) and domains without an evidence code | The first parse skipped them and read the second domain as the N-terminus, flagging six TAS2Rs, SLC39A2 and POPDC1 in error; the review had already covered all eight and keeps them on the lane | no decision taken: bug fix / Claude Code session l1-nterm-lane-check, found by its literature-review subagents |
| 2026-11-11 | Verdicts for the 46 reviewed members | Literature review by four subagents (GPCR signal-peptide classes and starts; GPCR orientation and helix count; non-GPCR transporters; non-GPCR isoforms and small proteins) | recommendation, not yet confirmed by Becca / literature-review subagents of Claude Code session l1-nterm-lane-check |
| 2026-11-12 | `design_id` added as the first column of every output and input table this script writes or reads (verdict, screen, review and SignalP tables; the script stops if the review's `design_id` disagrees with the manifest) and quoted in the write-up; cyno SLC22A17 (`l1-p2606`) noted as retired by the Compara 116 adoption | `candidate_id` changes with PR #283's renumbering; `design_id` does not | Becca (2026-11-12, via the reviewer-package session) / Claude Code session l1-nterm-lane-check |
| 2026-11-12 | Rekeyed on PR #283's final redraw (9928f4b1) and #280's matching rung table (79c1c464): 566 members; review joined by `design_id`, with the reviewed sequence hash checked; 9 reviewed members retired or moved off the rung by #283; four newly flagged members reviewed (cyno GPRC5A, mouse CNR1 and the HTR7 isoform stay; cyno NIPAL4 to replace) | #283's isoform gate, Compara 116 adoption and reference fixes change the rung's membership and numbering | integration session's request, 2026-11-12 / Claude Code session l1-nterm-lane-check; new verdicts are recommendations, not yet confirmed by Becca |
| 2026-11-12 | Merged the multipass session's delegated decisions (671564f4); recorded where they differ from this check (SLC61A1, cyno NIPAL4, OPN3 isoform) and sent the evidence to that session | Two analyses of the same members should not disagree silently | no decision taken: record of a disagreement / Claude Code session l1-nterm-lane-check |
| 2026-11-12 | The multipass session reversed its three differing calls on this check's evidence (6070d861): SLC61A1 kept on the rung, cyno NIPAL4 and the OPN3 isoform replaced | v2, structure-model, charge and SignalP evidence outweigh the family-homology argument | multipass session under Becca's delegation / Claude Code session l1-nterm-lane-check (record only) |
| 2026-11-12 | SLC61A1 moved to internal after all by the writer and the multipass session (5f21e02a), on UniProt's 12-helix N-in annotation and the robustness of a gated internal site to a disputed orientation; this check's verdict is kept as the topology record and noted as their rejected alternative | The two disagree on how to act on disputed evidence, not on the evidence | l1-reviewer-package and the multipass session under Becca's delegation / Claude Code session l1-nterm-lane-check (record only) |
| 2026-11-12 | Rekeyed on PR #283's second pass (6bf341f2) and #280's matching rung (c9e67f3e): 561 members; MFSD2A ×2, AQP12B and GLP2R moved off the rung and the OPN3 isoform retired by #283; OR8G3 re-reviewed on its ATG-recoded sequence and kept on the lane | #283 applied the lane overrides and the ATG recode; the review table's sequence check stopped the run until OR8G3 was re-reviewed | integration session's request / Claude Code session l1-nterm-lane-check; OR8G3 verdict a recommendation, not yet confirmed by Becca |
| 2026-11-13 | Merged `main` after PR #280 and #283 landed (13f3c1b8 onwards); rerun on the merged rung: 561 members, 560 hold, TSPO2 unresolved. Cyno NIPAL4 re-reviewed on its Met215 restart and kept on the lane (was: replace) | #283 restarted cyno NIPAL4, CHRM3 and OR3A2; the review's sequence check stopped the run until NIPAL4 was re-reviewed | integration request / Claude Code session l1-nterm-lane-check; NIPAL4 verdict a recommendation, not yet confirmed by Becca |
| 2026-11-13 | SignalP snapshot joined by `design_id`, with its 70-residue prefix hash checked, instead of by `candidate_id` | Main's renumbering silently detached 109 members from their SignalP call; no verdict had used the detached state | no decision taken: bug fix / Claude Code session l1-nterm-lane-check |
| 2026-11-13 | SLC61A1 (human, mouse) recorded as staying on this lane, with the move to internal as the rejected alternative | The writer (l1-reviewer-package) kept it on the rung on 2026-11-13 and l1-integration confirmed | l1-reviewer-package under Becca's delegation / Claude Code session l1-nterm-lane-check (record only) |
