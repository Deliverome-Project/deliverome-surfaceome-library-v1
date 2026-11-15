# N-after-SP members: placement rule, Igκ leader, and what blocks the build

This note turns the N-terminal reviewer package ([`BRIEF.md`](BRIEF.md)) into build decisions
for the 1,300 L1 members on the "N after SP" rung. It is for Becca. Every count comes from
[`nterm_build_options.py`](nterm_build_options.py) → [`nterm_build_options.json`](nterm_build_options.json),
with one row per member in `nterm_build_options.tsv`. Topology is DeepTMHMM v1, the cleavage
site `c` is the later of the DeepTMHMM v1 and SignalP 6.0 calls, and "preserved" means SignalP
6.0 keeps the native cleavage site with the 51-aa p003 cassette inserted (sweep
`sweep_rp_1792381414`).

**The cassette and the site rule are settled.** Becca chose the 51-aa p003 tandem on
2026-11-10: `GGGGS`–SpyTag003–`GGGGSGGGGS`–P-ALFA-P–`GGGGS`, the cassette every re-prediction
here used. On 2026-11-11 she chose the re-prediction-gated site rule.
[`nterm_site_picks.py`](nterm_site_picks.py) applies it: one site per member, with a status.

| Status | L1 members |
|---|---:|
| Ready: a gated site before the predicted fold | 604 |
| Held: preserved sites only inside the fold (decision 3) | 354 |
| Held: cleavage moves but stays upstream of the cassette (not yet accepted) | 53 |
| Held: Igκ swap recommended (decision 4) | 111 |
| Refused: no admissible interval (decision 5) | 178 |

**The 604 ready arms are built** ([`../../build_l1_nterm_arms.py`](../../build_l1_nterm_arms.py) →
`constructs_l1_nterm.tsv`, `.fasta`, `_summary.json`). They use the C-terminal builder's layout
and checks (PR #284) and the p003 tandem's own DNA, one encoding shared by every N-terminal
member. Per-member encodings are not needed here. The barcode-swap model (`l1_swap_risk.py`)
and the integrity tests measure the stretch a member shares with another *ending at the
barcode*, and an N-terminal cassette sits upstream of the ORF's own C-terminus, so it never
enters that stretch. Each member's swap risk is therefore unchanged.

All 604 checks pass:
- the protein is exactly `native[:n] + cassette + native[n:]`, and its SHA-256 is the one
  SignalP scored;
- the cassette goes in after nucleotide `3n` of the domesticated ORF;
- there is one stop, at the end;
- no seam creates a BbsI, Esp3I or BsaI site, and in the p0039 backbone each vector has the
  backbone's 4 sites plus the 2 wrapper Esp3I.

Barcodes are a 20-N placeholder until the order is frozen, as for the C-terminal arms.

Barcodes also wait for the frozen order (cover note #1). No arms are built here.

**Ortholog panel.** These counts use the ortholog panel on `main` as of 2026-11-10, merge
152176f0. That panel has 103 N-after-SP ortholog members; 84 of them are refused for want of
a SignalP row. The Compara 116 ortholog panel Becca adopted on 2026-11-10 is not on `main`
yet. Nothing here lists members by hand: every script reads `l1_tag_paths.tsv` and the
tranche manifest. Once that panel lands, rerun in this order:
1. `snapshot_inputs.py` under a new date, to fetch predictor rows for the new members;
2. `nterm_placement.py`;
3. `nterm_build_options.py`.

New orthologs without a public SignalP or disorder row will be refused until they are
predicted.

## Three rules, compared on the 1,101 re-predicted members

Counts are those in [`nterm_build_options.json`](nterm_build_options.json) (`rule_verdicts_over_repredicted_members`, `reprediction_gated`), on the current membership; `nterm_build_options.py` regenerates both.

| Rule | Site | Preserved | Moved | Crosses the junction | Lost | Outside the interval |
|---|---|---:|---:|---:|---:|---:|
| Cleavage + 1 | the first mature residue | 757 | 57 | 107 | 180 | — |
| Disorder-anchored | most disordered residue in `c+1…c+20` (`tag_placement.post_sp_junction`) | 758 | 61 | 83 | 156 | 43 |
| Re-prediction-gated | a candidate whose tagged design is preserved: before the fold first, then softening ≤ 0.07, then highest tagged cleavage probability, then the smallest offset | 941 | — | — | — | — |

The gated rule finds a working site for 183 more members than either fixed rule, because
the fixed rules ignore the measurement already in hand. Its picks are mostly close to the
cleavage site: offset +1 for 287 members, +2 for 298, +3 for 90. For 176 of its picks the
tagged cleavage-site probability still falls by more than 0.07. The record found that moving
the junction does not recover a softened site, so those are flagged rather than avoided.

## Recommended placement, per subgroup

| Subgroup | L1 members | Recommendation | Why |
|---|---:|---|---|
| Preserved at a site before the fold | 602 | Re-prediction-gated site | The site is outside the predicted fold and SignalP keeps the cleavage site |
| Preserved only inside the fold (mostly "no room" members) | 353 | Re-prediction-gated site, built only if Becca accepts the fold risk; otherwise hold for p0028 (FOLR1) display data | SignalP cannot see fold insertion. FOLR1's arm is the only test of a cassette inside a predicted fold |
| Cleavage moves, but stays upstream of the cassette | 53 | Use the moved-but-upstream candidate, flagged | The cassette is still at the mature N-terminus. The record judged HEG1's +1 shift benign on this ground |
| Cleavage crosses the junction at every candidate | 21 | Igκ leader swap | All 21 keep the leader's cleavage site at every candidate once swapped (below) |
| Signal peptide lost at every candidate | 88 | Igκ leader swap | All 88 keep the leader's cleavage site at every candidate once swapped (below) |
| Added after the sweep, now re-predicted | 5 | PRSS8, KIR2DS4, ADGRE1: the gated site. CD58, PCDHGA7: Igκ swap, not yet re-predicted | PRSS8, KIR2DS4 and ADGRE1 keep their native site at 4, 2 and 9 candidates. CD58 and PCDHGA7 lose it at both of theirs (run `nterm_sp6_2026_11_10`) |
| Refused: no SignalP row | 121 | Admit the ortholog SignalP release, then rerun the interval and re-prediction | 81 of them have an exact-sequence row in `sp6_2026_09_28_ortho` (one-off database check, 2026-11-04) |
| Refused: DeepTMHMM v1 calls a signal peptide, SignalP 6.0 does not | 57 | **Decided 2026-11-14 (coordinator, under Becca's delegation): Igκ swap only where UniProt SIGNAL or DeepTMHMM v2 confirms a cleaved signal peptide within 3 residues; signal anchors routed; other sites re-anchored; unconfirmed held** (`sp_disagree_review.py`) | SignalP's OTHER call may mean a signal anchor, not a cleaved peptide, so deleting it could change topology |

The full exception lists are in `nterm_build_options.json` under `exceptions`, each
member given as its member ID (`l1-mNNNN`), then its design ID (`l1-pNNNN`), then its gene.

## Should an Igκ leader replace weak or uncertain native signal peptides?

**Not on native confidence alone. Use it for members whose native signal peptide fails
re-prediction once tagged.** That is the rule `build_sp_swap_arms.py` already states. The
data:

| Native cleavage-site probability | Members | Cleavage + 1 not preserved | No candidate preserved |
|---|---:|---:|---:|
| < 0.50 | 33 | 24 | 19 |
| 0.50–0.80 | 114 | 64 | 28 |
| 0.80–0.90 | 75 | 25 | 18 |
| 0.90–0.95 | 61 | 17 | 7 |
| ≥ 0.95 | 834 | 217 | 90 |

A weak native cleavage site does raise the risk: 56% of members below 0.5 have no preserved
candidate (18 of 32), against 11% at 0.95 or above (89 of 823). But a confidence threshold is
the wrong trigger. 96 of the 160 members with no preserved candidate have a native
cleavage-site probability of 0.9 or more, so a threshold misses most of them. And 14 of the 32
weak members keep their site at some candidate, so a threshold would swap signal peptides that work.

**What re-prediction says about the leader.** The pinned SignalP run `nterm_sp6_2026_11_10`
(2026-11-10; see "Reproducing the SignalP runs") swapped the Igκ leader plus Asp in for residues `1…c` of all 109 failing
members (88 lost, 21 crossing). It put the tandem after the same mature residue as each native
candidate.

| Arm | Designs | SignalP 6.0 call | Cleavage site | Cleavage-site probability |
|---|---:|---|---|---|
| Igκ-swapped, untagged | 109 | SP for all 109 | residue 20 (the leader's end) for all 109 | 0.813–0.984, median 0.981 |
| Igκ-swapped, tagged | 346 | SP for all 346 | residue 20 for all 346 | 0.980–0.983 |

Every one of the 346 tagged designs keeps the leader's cleavage site, so each of the 109
members has a usable site once swapped. This matches VSIG10's earlier canary (0.348 native,
lost when tagged, 0.982 swapped and tagged; `constructs_sp_swap.tsv`).

Read it for what it is. It is a prediction on a 70-residue window that is mostly the same
leader in every design, so its uniformity is expected rather than informative per member. It
says nothing about expression, folding or display, and the leader adds an Asp to each mature
N-terminus. No Igκ-swapped arm has expression data yet; VSIG10's pair (p0026 native, p0048
Igκ) will be the first.

The ALFA-only cargo is not measured, because Becca chose the tandem.

## Reproducing the SignalP runs

SignalP 6.0 runs on Modal through [`signalp6_modal.py`](signalp6_modal.py), to the standard of
the DeepTMHMM runner: the licensed package and checkpoints sit on Modal volumes and are
fingerprinted on every run. Each run writes `provenance.json` per shard and for the whole run.

| Piece | Pinned as |
|---|---|
| SignalP 6.0, version `6.0+h` | DTU tarball `signalp-6.0i.slow_sequential.tar.gz`, SHA-256 `3193855d628c78337c9df28a0db76d977588b6087e6d118d267bd461564ba024`; its 16 code files on volume `signalp6-package`, each checked against the tarball on every run |
| Model checkpoints | volume `signalp6-models`. The 7 shipped checkpoints (`cpu/`) are checked against the tarball on every run. The copies under `gpu/`, used to predict, are hashed and recorded, and are byte-identical to the tarball's. All hashes: [`signalp6_identity.json`](signalp6_identity.json) |
| Container | Python 3.10.17, torch 1.13.1+cu117 (CUDA 11.7, cuDNN 8500), numpy 1.26.4, matplotlib 3.8.4, tqdm 4.66.5; full `pip freeze` in each run's provenance |
| Hardware and flags | Tesla T4; `--mode slow-sequential --organism eukarya --bsize 64 --format none`, through the package's own `signalp6` entry point |
| Input order | FASTA order kept, shards of 1,024 at most, each shard's sequence order recorded (SignalP's numbers depend on batch composition) |

**Validation.** Rerunning shard `s0000` of the original sweep `sweep_rp_1792381414` (1,024
sequences, original order) reproduced every call and cleavage site. The largest difference in
any probability was 0.0 ([`signalp6_validation_sweep_s0000.json`](signalp6_validation_sweep_s0000.json)).

**Run `nterm_sp6_2026_11_10`** (from commit `9dd9d73a`, clean tree). It covers 486 sequences
from [`build_signalp_inputs.py`](build_signalp_inputs.py) (`signalp_inputs.tsv` maps each one
to its member, arm and site). Calls and raw outputs are at
`processed://sprint4/s4e9-orf-design-fg-panel/reviewer-package/nterm/signalp6-runs/nterm_sp6_2026_11_10/`,
and [`analyse_signalp_run.py`](analyse_signalp_run.py) reads them into
`signalp_run_verdicts.tsv` and `signalp_run_summary.json`.

Its one input that cannot be rebuilt from this repo is the package volume: it was uploaded
from Becca's licensed tarball on 2026-11-10. The tarball's hash and every file's hash are the
record. To restore it, upload the same tarball's `signalp-6-package/` code and its
`tar_member_sha256.json`; the runner refuses any file that differs.

## Decisions needed from Becca

| # | Decision | Options | Recommendation |
|---:|---|---|---|
| 1 | N-terminal cassette | — | **Decided 2026-11-10: the 51-aa p003 tandem** |
| 2 | Site rule | — | **Decided 2026-11-11: re-prediction-gated** |
| 3 | The 353 members preserved only inside the fold | Build now / hold for p0028 / route to another rung | Hold for p0028 unless the order date forces it |
| 4 | The 109 members that fail everywhere (and CD58, PCDHGA7) | Igκ swap / no N-terminal tag | Igκ swap: SignalP keeps the leader's site at all 346 tagged designs |
| 5 | The 178 refused members | Admit the ortholog SignalP release (81), then inspect the 57 | Yes, in that order |

## Every member decided

[`nterm_decide.py`](nterm_decide.py) writes one decided arm per member to `nterm_decisions.tsv`,
with the evidence, the reason, the rejected alternative, `decided_by`, `implemented_by` and the
refold rule; [`../../build_l1_nterm_arms.py`](../../build_l1_nterm_arms.py) builds every
native-site and Igκ arm, and `--check` rebuilds them byte for byte. Current counts are in
`nterm_decisions_summary.json` and the companion brief's section 7; each change to them is a
dated row below. At main after #301 and the second-source review (2026-11-14): 1,309 members,
1,089 native-site, 183 Igκ, 35 held, 1,272 built; OR10G6 `routed_multipass`; ATRN, NRG2 and both ENTPD2 members `routed_internal`, OR11H6 and TMEM63B `routed_multipass`.

## Decision log

| Date | Was | Now | Why | Decided by / implemented by |
|---|---|---|---|---|
| 2026-11-10 | Placement measured, no build rule | Three rules compared per member; re-prediction-gated recommended; the build waits on the cassette and the site rule | Becca split tag work by rung; this note is the N-after-SP rung's input | No decision taken; proposals by the Claude Code session nterm-reviewer-package (PR #254), on Becca's request relayed by the l1-recode-dups session |
| 2026-11-10 | Cassette open (51-aa tandem or ALFA-only) | The 51-aa p003 tandem | Becca: "we want the tandem" | Becca / recorded by nterm-reviewer-package (PR #254) |
| 2026-11-10 | 5 members never re-predicted; the Igκ swap measured only on VSIG10 | Both re-predicted with a pinned SignalP 6.0 run (`nterm_sp6_2026_11_10`); the runner reproduces sweep shard s0000 exactly | Close the coverage gap and answer the Igκ question with data | Becca (asked for the Modal run, with versions and hashes recorded as for DeepTMHMM) / nterm-reviewer-package (`signalp6_modal.py`, PR #254) |
| 2026-11-11 | Site rule open | Re-prediction-gated; one site per member in `nterm_site_picks.tsv` (604 ready) | Finds a site that keeps the signal peptide for 185 more members than either fixed rule | Becca / nterm-reviewer-package (PR #254) |
| 2026-11-11 | No N-terminal inserts | 604 arms built, p003 tandem, shared encoding | Becca asked this session to proceed; the swap model needs no per-member encoding for an N-terminal cassette | Becca / nterm-reviewer-package (`build_l1_nterm_arms.py`, PR #254) |
| 2026-11-11 | Trimming the leading GGGGS considered elsewhere | The full 51-aa cassette, leading GGGGS included, everywhere; the 604 arms already use it | With SpyTag003 directly after Met1, DeepTMHMM v1 called a spurious signal peptide inside the tag on 12 of 62 type III constructs, none with the full cassette (PR #282) | Becca ("keep full"), relayed by l1-recode-dups / recorded by nterm-reviewer-package (PR #254) |
| 2026-11-12 | Held groups awaiting Becca | Every member decided and built: 1,307 arms (1,115 native, 192 Igκ), 7 with no N-terminal arm | Becca: "use your best judgment to complete … by the morning" | nterm-reviewer-package (delegated by Becca 2026-11-12, "use your best judgment") / nterm-reviewer-package (PR #254) |
| 2026-11-12 | Refold verdicts not yet applied | First l1-tag-refold verdicts (98 folded: 89 pass, 8 fail, 1 no native model) applied. Only `fail` changes an arm; `not_folded` and `no_native_model` stay pending. Of 8 failed native arms, 4 move to their next SignalP-passing site (NRN1 29, LYPD6 26, NCR3 20, GABRA2 30; sent back for refold) and 4 with no other SignalP-passing form get no N-terminal arm (TMEM95, ICOS, FOLR1, LYVE1). Arms: 1,111 native, 192 Igκ, 11 none (1,303 built) | A failed refold means the cassette disturbs the native fold at that site; an Igκ swap keeps the same mature insertion residue, so it does not fix that | nterm-reviewer-package (delegated by Becca 2026-11-12, "use your best judgment") / nterm-reviewer-package (PR #254) |
| 2026-11-12 | First refold verdicts (98 folded) | Full l1-tag-refold verdicts applied (1,271 folded: 1,195 pass, 76 fail; 40 not yet folded stay pending). Every fail has a disposition: 61 members move to their next SignalP-passing form (56 native site, 5 Igκ; sent back for refold as `nterm_refold_requests_2026-11-12b`), 11 with no SignalP-passing form left get no N-terminal arm (NRN1, LYPD6, NCR3, GABRA2 after a second fail; TMEM95, ICOS, FOLR1, LYVE1, ITGB3, LRFN2). No current arm is a failed form. Arms: 1,104 native, 192 Igκ, 18 none (1,296 built) | A fold-check failure is evidence the cassette disturbs the native domain at that site; the next site that keeps SignalP cleavage is the least-changed design to try, and is refolded before it counts | nterm-reviewer-package (delegated by Becca 2026-11-12, "use your best judgment") / nterm-reviewer-package (PR #254) |
| 2026-11-13 | Site rule blind to pro-peptides | The rule knew the signal peptide and disorder only. Against UniProt 2026_04 PROPEP features (`propeptide_check.py`, snapshot hashed on S3), 63 built arms put the cassette inside an N-terminal pro-peptide that a convertase removes on the way to the surface (cadherins, desmogleins, desmocollins, ADAMs, MT-MMPs, BACE1/2, SORT1, SORCS3, PAM, SELPLG, PRRG2/4, TGFA, HBEGF, AREG). Those cassettes move to the mature N-terminus (after the pro-peptide, keeping the R-X-[KR]-R site whole): 61 keep the SignalP cleavage site (run `nterm_sp6_2026_11_12f`) and go to refold; CDH5 does not and has no N-terminal arm. Activation-only pro-peptides (F2R, F2RL1-3, PRSS8/21/41, MEP1B) and RGMA (fragments stay disulfide-linked) keep their site | A cassette in a constitutively removed pro-domain never reaches the surface; at an activation-only pro-peptide the moved cassette would cap the tethered ligand or the zymogen's new N-terminus. MMP14 and BACE1 junctions at the mature N-terminus are in the literature corpus. Classical cadherins need a free mature N-terminus (Trp2) for strand-swap adhesion, so these arms report surface display, not adhesion | nterm-reviewer-package (delegated by Becca 2026-11-12, "use your best judgment") / nterm-reviewer-package (PR #254) |
| 2026-11-13 | One alternative per refold round | 45 of 59 next-site arms failed again, most on epitope packing one residue from the failed site. Every remaining SignalP-passing alternative of a failing member (400 forms) is now sent to refold at once (`refold_requests_alternatives.tsv`), so one round settles each member | Adjacent sites share the failure mode; one residue per round would not converge overnight | nterm-reviewer-package (delegated by Becca 2026-11-12) / nterm-reviewer-package (PR #254) |
| 2026-11-13 | Furin motifs without annotation | 28 built arms sit upstream of an R-X-[KR]-R motif within 60 residues of cleavage, with no UniProt pro-peptide (e.g. ADAM8, ADAM9, CDH22, CDH26, whose UniProt CHAIN starts at the signal-peptide cleavage). Kept, and listed in `propeptide_check.tsv` | A motif alone is weak evidence (it occurs by chance and in mature domains); UniProt annotation is the single documented source the rule acts on | nterm-reviewer-package (delegated by Becca 2026-11-12) / nterm-reviewer-package (PR #254) |
| 2026-11-13 | Any refold fail overturns an arm | A fail whose only failing criteria are lDDT values within 0.01 of their bar (native region 0.85, site window 0.70) counts as a pass, labeled `borderline` in `refold_verdict`: SNORC, FOLR1, LYVE1, SMO, SLITRK2. Epitope packing, a buried tag and a helix inside the cassette stay decisive at any margin | The verdict rests on one ColabFold model (model_1_ptm, 3 recycles), whose lDDT moves by more than 0.01 between runs; FOLR1 at 0.698 is not evidence against its site | nterm-reviewer-package (delegated by Becca 2026-11-12) / nterm-reviewer-package (PR #254) |
| 2026-11-13 | Live refold verdicts read by key | `decision_inputs.py` pins every S3 input of the decisions and arms by SHA-256 (23 files) and freezes l1-tag-refold's verdicts under a content-addressed key; `nterm_decisions_summary.json` records the verdict SHA it used | The verdict file is republished as folds land; a decision must name the exact verdicts it rests on | nterm-reviewer-package (overnight review relayed from Becca) / nterm-reviewer-package (PR #254) |
| 2026-11-13 | Refold verdicts `7c8b5bcf` | Verdicts `ae63ff51` applied (1,229 pass, 47 fail, 25 replacements pending). 56 of 61 pro-peptide relocations pass; AREG, BACE1, CDH15 and ADAM32 fail refold at the mature N-terminus and are held, as is CDH5. Arms: 1,069 native, 190 Igκ, 55 held | Same rules, new evidence; a cassette left inside a removed pro-peptide is not an alternative | nterm-reviewer-package (delegated by Becca 2026-11-12) / nterm-reviewer-package (PR #254) |
| 2026-11-14 | Verdicts read from l1-tag-refold's live S3 file | Verdicts read from `tag_refold_verdicts.tsv` committed on main by #301 (3a7dfbd9), pinned by SHA-256 in `decision_inputs.json`. A replacement already folded and passing is taken first. Every built arm has a verdict: 1,068 native, 190 Igκ, 56 held (1,258 built). SLC39A6 failed at every SignalP-passing site and is held. CIST1 (native-region lDDT 0.847) joins the borderline passes | Finalize #254 for merge from verdicts that live in git | coordinator decision under Becca's delegation, 2026-11-14 / nterm-reviewer-package (PR #254) |
| 2026-11-14 | Borderline rule visible only in the verdict text | `refold_rule`, `refold_rule_metric`, `refold_rule_threshold`, `refold_rule_decided_by`, `refold_rule_date` on every decision row (`lddt_borderline_counted_as_pass`, `l1-tag-refold verdict`, or `hold_without_refold_verdict` under `--hold-pending`) | So a consistency check sees that this rung differs from #301 by rule, not by omission | coordinator request, 2026-11-14 / nterm-reviewer-package (PR #254) |
| 2026-11-14 | ATRN and NRG2 counted as N-after-SP holds | `routed_internal`: #245's proposal to move them onto this rung was evaluated here and rejected (no form keeps a SignalP cleavage site); both stay internal, where `l1_tag_paths.tsv` places and builds them. N after SP is 1,312 members: 1,258 built, 54 held | The member count matches the placement of record (l1-integration, 2026-11-14) | l1-integration request / nterm-reviewer-package (PR #254) |
| 2026-11-14 | Igκ swap applied to 54 members on DeepTMHMM v1's signal peptide alone | Second-source rule (`sp_disagree_review.py`): 42 confirmed and kept; OR11H6, TMEM63B, ENTPD2, CATSPERH are signal anchors or first TM helices, routed (residue-1 multipass, internal, C-terminal) and held here; 10 with a cleaved signal peptide elsewhere re-anchored there (UniProt first, else DeepTMHMM v2; SignalP run `nterm_sp6_2026_11_14g` keeps the leader's site for all ten), held until refold; FCAMR and MDGA2 held. Reason strings corrected; `sp_review_verdict` and `sp_sources_checked` per row | The written rule said "inspect before choosing; no Igκ by default"; an external audit (2026-11-14b, B2) found the code broader | coordinator under Becca's 2026-11-12 delegation / nterm-reviewer-package |
| 2026-11-14 | Round-2 candidates from DeepTMHMM v1's cleavage only; the softening bound restated in three modules | A round-2 site must follow the parent's SignalP cleavage (SMO ortholog 34→35; ADAM33 isoform past 29, then out of its pro-peptide), held until refold. One pick key (`pick_key`, `SOFTENING_BOUND`) for every caller; unit tests in `tests/test_nterm_site_rule.py`; no other arm changed | Audit Q1: two arms had the cassette at or before SignalP's cut | nterm-reviewer-package (audit follow-up) |
| 2026-11-14 | 34 native arms with tagged cleavage probability below 0.5 unlisted | Listed as flagged in the brief; no floor, no swap | Every Igκ form scores the leader's 0.97-0.98 regardless of member, so it is no evidence for swapping | coordinator default / nterm-reviewer-package |
| 2026-11-14 | Re-anchor at any independent site, UniProt first, else DeepTMHMM v2 | Re-anchor only on curated evidence: the UniProt SIGNAL site for 7 (evidence code recorded per row; all ECO:0000255 except PTPRG, none), and CNNM2 at v2's 65 on a verified primary report of a ~64-residue cleaved signal peptide (de Baaij et al. 2012, PMID 22399287). CDH3 and GRIA1 isoforms (v2 only; no alignment to a canonical UniProt SIGNAL) and FCAMR, MDGA2 held. OR11H6 and TMEM63B retired as `routed_multipass` (built by #337), ENTPD2 as `routed_internal`; CATSPERH held. Several v1 calls end at 60, an apparent cap. Now 1,309 on the rung: 1,066 native, 175 Igκ, 68 held, 1,241 built | A second predictor alone is not evidence enough to delete residues; rejected alternative: hold all 12 (loses seven curated signal peptides). Reversed for a member if its Igκ form fails refold | coordinator under Becca's 2026-11-12 delegation / nterm-reviewer-package |
| 2026-11-14 | Second-source rule covered the 58-member class only | Extended to 12 members with no census SignalP row and no native SignalP signal peptide (six scored natively for the first time, run `nterm_sp6_2026_11_14h`): 8 confirmed, cyno ENTPD2 deferred to the internal rung, CD84, IL17RC and CACHD1 held. Every Igκ row names its deletion evidence; a test enforces it. 1,309 on the rung: 1,066 native, 172 Igκ, 71 held, 1,238 built | The same DeepTMHMM-v1-only deletion existed outside the audited class | coordinator under Becca's 2026-11-12 delegation / nterm-reviewer-package |
| 2026-11-14 | Ortholog swaps decided on DeepTMHMM alone | The human ortholog's UniProt SIGNAL, aligned, is recorded as support; CD84, IL17RC and CACHD1 released at the aligned site (≥80% over ±5, no gap); the 49 second-source swaps split by source, DeepTMHMM-v2-only (16) stated as the weakest. 1,309 on the rung: 1,066 native, 174 Igκ, 69 held, 1,240 built | DeepTMHMM v2 confirming v1 is one tool twice | coordinator under Becca's 2026-11-12 delegation / nterm-reviewer-package |
| 2026-11-14 | Ten new forms held until refold | All ten pass (l1-tag-refold, verdicts SHA-256 0d2b1a23…): the eight re-anchored Igκ swaps and the re-sited SMO ortholog and ADAM33 isoform are built. The pro-peptide snapshot is read from the append-only `external://uniprot/2026_04/nterm-propeptides/` copy (byte-identical). 1,309 on the rung: 1,068 native, 182 Igκ, 59 held, 1,250 built; IL17RC awaits refold | Refold verdicts returned | nterm-reviewer-package |
| 2026-11-14 | Cyno ENTPD2 Igκ arm kept until its internal arm was built | Retired: `routed_internal` (internal arm after 458 passes gate and refold, PR #341). No Igκ deletion now lacks cleaved-signal-peptide evidence. 1,308 on the rung: 1,068 native, 181 Igκ, 59 held, 1,249 built | The deferral's condition was met | ranker session / nterm-reviewer-package |
| 2026-11-14 | IL17RC held until refold | Released: the Igκ swap re-anchored at the human ortholog's aligned site 20 passes refold (#348). Both ENTPD2 members leave the rung for the internal rung (#346). 1,308 on the rung: 1,068 native, 182 Igκ, 58 held, 1,250 built; nothing on the rung awaits a verdict | Refold verdict returned | nterm-reviewer-package |
| 2026-11-14 | Arms whose every form failed refold were held | The residue-1 multipass rung's disposition rule applied verbatim where the fold is intact: 23 members built and flagged on epitope packing or a buried SpyTag003 (native-region lDDT ≥ 0.85, site lDDT ≥ 0.70 where scored, ALFA RSA ≥ 0.25), at the form with the fewest packed residues; a helix inside the cassette disqualifies (EFNA1, PDCD1, FAM174A). `refold_disposition` per row; a test requires it on every built arm that fails refold. 1,308 on the rung: 1,090 native, 183 Igκ, 35 held, 1,273 built. CDH5 (Igκ form moves the cut) and VTCN1 (signal anchor by DeepTMHMM v1 and v2; SignalP OTHER) stay off | One rule for both rungs: 19 genes lost their only member to the inconsistency; a packed epitope in one model of a flexible cassette is weak evidence, and the library reads two epitopes. Would change: bench SpyCatcher003:ALFA capture ratio | coordinator under Becca's 2026-11-12 delegation and 2026-11-14 rescue request / nterm-reviewer-package |
| 2026-11-14 | Orthologs built whether or not the gene's human member was | Ortholog rule: an ortholog arm is built only if a human member of its gene is built on some rung. IL12RB1 cyno l1-p2947 held (human l1-p0913 fails refold on site lDDT). A test enforces the rule. 1,308 on the rung: 1,089 native, 183 Igκ, 36 held, 1,272 built | Becca: "Never substitute a mouse or cyno ortholog for a human and don't build orthologs for genes with no good human proteoform" | Becca / nterm-reviewer-package |
| 2026-11-14 | No floor on weak native signal peptides (34 flagged, kept) | Weak predicted SP (native or tagged cleavage-site probability < 0.5, 53 arms) moves to Igκ where the Igκ form keeps the leader's cleavage and passes refold; otherwise the native arm stays, flagged. Cut calls more than 5 residues apart are not swapped but flagged for Becca. Thresholds proposed here, adopted by the coordinator 2026-11-15, not yet confirmed by Becca; rejected: either measure alone, and the native SP probability | Becca: "we are relying more on igk now for ... weak predicted SP like VSIG10" | coordinator under Becca's delegation / nterm-reviewer-package |
| 2026-11-15 | OR10G6 held (no form keeps SignalP cleavage) | Retired as `routed_multipass`: the multipass rung builds it Igκ + tag at Met1 under its own lane override; this rung adds none. 1,307 on the rung, 35 held | Displayed leader-free in Tedman's screen; no signal peptide by SignalP 6.0 or UniProt | multipass session / nterm-reviewer-package |
| 2026-11-15 | 49 weak-SP members awaiting refold of their Igκ forms | All 49 swap to Igκ: 47 pass refold; ICAM4 and LYPD1 fail only on epitope packing with the fold intact and ALFA exposed, so they are built and flagged under the disposition rule. TMEM81, MXRA8, IL1RAPL2 stay native and flagged for Becca; KIAA0319 stays native. 1,307 on the rung: 1,040 native, 232 Igκ, 35 held, 1,272 built | Refold verdicts returned (l1-tag-refold) | coordinator under Becca's delegation / nterm-reviewer-package |
