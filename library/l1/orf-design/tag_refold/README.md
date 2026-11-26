# Refolding tagged L1 constructs

Every L1 member carries the 51-residue tag cassette, `GGGGS`–SpyTag003–`GGGGSGGGGS`–`P`–ALFA–`P`–`GGGGS`
(SpyTag003 = `RGVPHIVMVDAYKRYK`, ALFA = `SRLEEELRRRLTE`), placed by rung: N-terminal after the signal
peptide (SP), C-terminal, type III, residue-1 multipass behind an Igκ leader, or an internal loop.
This folder predicts each tagged protein's structure and asks two questions of it:

1. **Does the cassette disturb the native protein?** The native region of the tagged model is
   compared residue for residue with a model of the same native protein.
2. **Can the epitopes be reached?** SpyTag003 and ALFA should be solvent-exposed and not packed
   onto the protein.

Becca asked for the refolds on 2026-11-11 ("for proteins that have a clear position at the terminus
refold those to double check and for n term proteins run signal reprediction and disorder etc as we
had planned to select the site then also refold those"), within $230 of the night's $300 of Modal.
On 2026-11-12 she delegated the calls ("use your best judgment to complete a new and improved
evidence based l1 completely by the morning"), so each verdict here is a decision, not a proposal.

## What is folded, and in what order

The inputs are each rung's own final export after the PR #283 renumbering, snapshotted write-once to
S3 and listed with their SHA-256 in [`tag_refold_inputs.tsv`](../tag_refold_inputs.tsv). One protein
is folded per tagged SHA-256; [`tag_refold_sources.tsv`](../tag_refold_sources.tsv) maps every rung's
own id onto it. A night cannot fold everything, because the public MSA server returns only about two
native alignments a minute, so the queue is tiered by the value of the decision each fold informs.

| Rung (owner session, PR) | Tier 1: the refold can change a call | Tier 2: failure-rate sample | Tier 3 |
|---|---|---|---|
| N after SP (nterm-reviewer-package, #254), 1,311 arms | 591 Igκ arms and native-site arms inside the predicted fold, and 4 replacement arms for arms that failed | 100 before-fold arms | 616 |
| C-terminal (l1-cterm-build, #284), 216 | 44: structured terminus (native terminal pLDDT ≥ 85), M13 family, KLRF1, TNF ligands, ITM2B, named priorities | 100, including the 4 snorkel (TM-module) arms | 72 |
| type III (l1-nterm-type3, #282), 44 | all 44 | — | — |
| residue-1 multipass (l1-nterm-res1-multipass, #280), 559 mature chains, plus 12 spacer comparisons (`multipass_spacer`: GGGGSGGGGS between cassette and residue 2, for the 12 members whose native-region lDDT dropped) | the 12 spacer comparisons | 100: all 24 leader-free arms + 76 Igκ arms | 458 |
| internal loop (ranker-pu-segment, #245), 679 | both sites of the 100 members whose rank-2 loop is close to rank 1, and 20 decided sites that are their member's only site (one sent by message, recorded in [`tag_refold_extra_requests.tsv`](../tag_refold_extra_requests.tsv)) | 100 rank-1 sites | 359 |

Tier 2 samples take every k-th construct by tagged SHA-256, which is unrelated to anything about the
protein, so their failure rate estimates the rung's. A native protein with no model of exactly its
sequence ([`tag_refold_native_models.tsv`](../tag_refold_native_models.tsv): AFDB via main's structure
payload, or any ColabFold model on S3, matched by sequence) is folded untagged at its construct's tier.
Its MSA is fetched for the tagged form anyway, so the extra fold costs about $0.03.

## How each protein is folded

ColabFold 1.6.2 on Modal (`scripts/fold_tagged_modal.py`),
`alphafold2_ptm` model 1 only, 3 recycles, no Amber relaxation, no templates.

| Choice | Reason |
|---|---|
| One model, 3 recycles, no relax | Budget; the controls below measure what it costs in accuracy. |
| The **native** protein's MSA, with gap columns for the cassette or leader | Homologues carry no information about a synthetic tag, so a search on the tagged sequence returns the native family with the tag unaligned. One search serves every tagged form of a protein. |
| MSAs from the public ColabFold server, cached by native SHA-256 | Each native is sent once, in tickets of 20, three tickets at a time with the server's own back-off. |
| Igκ multipass and Igκ type III folded as mature chains | The leader is cleaved; N-after-SP Igκ arms are folded as precursors, so the leader's fate is visible. |

A fold costs about $0.03 on an A10G (benchmark: 20 folds, 10,929 residues, $0.68).

## Controls: what the cheap settings alone do

Seventeen L1 natives were refolded untagged with the same settings and compared with their AFDB
models: 10 spread over 150–1,200 residues, and 7 with signal peptides (`sp_control`).

| Measure (confident AFDB residues, pLDDT ≥ 70) | Range over the 17 controls |
|---|---|
| Native-region Cα lDDT | 0.897–1.000 |
| TM-score | 0.69–1.00 (multi-domain proteins lose TM on hinges; lDDT does not) |
| 30-residue window Cα lDDT, 1,471 windows | 1st percentile 0.70, median 0.99 |
| SP helix (≥ 7-residue helical run) kept, of 7 SP controls | 5 of 7; the other 2 become unresolved low-pLDDT ribbons |

## Verdict

A construct **passes** unless one of these holds; `verdict_reasons` names each failed test. The
thresholds were set on 2026-11-12 by l1-tag-refold under Becca's delegation, from the controls above;
the type III TM rule was made relative to the native model the same day, after IGSF5.

| Test | Fail when | Why this threshold |
|---|---|---|
| Native fold reproduced | native-region Cα lDDT < 0.85 | below every control (minimum 0.897) by a margin |
| Site neighborhood reproduced | lDDT of confident residues within 15 of the insertion < 0.70 | 1st percentile of control windows |
| SpyTag003 / ALFA exposed | mean relative SASA < 0.25 | the usual buried/exposed boundary on the Tien et al. 2013 maximum-ASA scale |
| Epitope not packed on the protein | ≥ 5 SpyTag003 or ALFA residues in confident contact | contact: heavy atoms within 4.5 Å of a native residue more than 4 positions from the junction, PAE < 10 Å both ways, tag residue placed (pLDDT ≥ 50). Five is about a third of either epitope. Linker contacts are reported, not counted. |
| No extra helix in the cassette | a confident (pLDDT ≥ 70) helix of ≥ 10 residues outside ALFA | ALFA is a helix by design; a 10-residue helix elsewhere is the "second signal or TM-like helix" type III and residue-1 designs must avoid |
| Type III TM stays one helix | TM1 is a helix of ≥ 15 residues in the native model but not in the tagged fold | 15 residues is the shortest helix that spans a membrane; AFDB leaves some TMs unmodelled, and that is not the tag |

lDDT carries the verdict rather than TM-score because it is superposition-free: a domain that swings
on a hinge costs only the contacts across the hinge. The **SP state** (`sp_state`: helix, unresolved,
placed_not_helical) is reported but never fails a construct, because the controls already lose 2 of 7
AFDB SP helices to the settings alone. Cleavage is SignalP's question, answered by
nterm-reviewer-package. TNF-family ligands are folded as monomers; their trimer interface is not
tested here.

## Results

One verdict per rung construct ([`tag_refold_verdicts.tsv`](../tag_refold_verdicts.tsv), 2026-11-13; 3,765 models).
The N-after-SP rung includes the 508 alternative sites nterm-reviewer-package sent for members whose
first arm failed, so its fail count is over sites tried, not members.

| Rung | Pass | Fail | Most common reason for failing |
|---|---:|---:|---|
| N after SP | 1,454 | 410 | an epitope packed onto the protein (169); site neighborhood (121); native fold (110) |
| C-terminal (incl. 47 snorkel and 6 stalk-site arms) | 219 | 20 | native fold (13) |
| type III | 40 | 4 | a confident helix inside the cassette (4) |
| residue-1 multipass | 489 | 72 | SpyTag003 buried (52): at the mature N-terminus it packs against the receptor |
| internal loop (incl. 52 alternatives, 4 last holds) | 641 | 101 | native fold (75); site neighborhood (60) |
| expression-rescue arms | 21 | 4 | site neighborhood (2) |
| multipass + 10-residue spacer (comparison) | 0 | 12 | native fold (12): the spacer rescues none of the 12 fold-disturbed arms |
| snorkel + 30-residue L1 linker (pre-registered test, 5 failing + 2 passing controls) | 0 | 7 | the longer linker rescues none of the 5 and breaks both controls (native fold, 6 of 7 lower) |

## The fold is the ordered protein

Each member's cassette DNA is diversified (member-specific linker codons in `l1_cassette_encodings.tsv`),
so [`verify_ordered_constructs.py`](verify_ordered_constructs.py) translates every built insert on the
draft-L1 integration branch (PR #299 at `5cfd91c0`) and checks it against what was folded
([`tag_refold_construct_check.tsv`](../tag_refold_construct_check.tsv)):

| Check | Result (2026-11-13) |
|---|---|
| translated protein SHA-256 equals the rung table's stated SHA-256 | 2,691 of 2,691 |
| cassette amino acids are exactly the p003 tandem, once | 2,691 of 2,691 |
| the insert carries the member's listed cassette encoding | 256 of 257 listed; TSPAN10 l1-p2284 moved to an internal site and carries a different encoding of the same amino acids |
| the protein (or the mature chain of an Igκ design) is a folded sequence | every construct, once the folds listed as unfolded in `tag_refold_ordered_unfolded.tsv` land |

The overnight review (PR #303 at `8d4b935e`) then recodes 4,655 codons in 751 inserts for Twist's
synthesis rules. [`tag_refold_recode_check.tsv`](../tag_refold_recode_check.tsv) applies each recode to its
#299 insert: all 4,655 are synonymous, each `from_codon` is where the table says, and all 751 recoded
inserts translate to the unchanged protein.

Neither codon diversification nor synthesis recoding changes a residue: the folded sequences are the
ordered proteins.

The **order of record** (`synthesis/2026-11-13-one-arm-renamed/l1_order_inserts_recoded.fasta`, SHA-256
`b426e433…`, 2,678 recoded inserts) is also a source of the queue: every ordered insert gets a live
verdict row, matched by protein (the precursor, or the mature chain of an Igκ design), so a construct
cannot drop out of the verdicts when the rung export that first carried it moves on. On 2026-11-13 all
2,678 ordered inserts carry the cassette amino acids exactly once and have a live fold (2,110 as the
precursor, 568 as an Igκ mature chain): **2,535 pass, 143 fail**. This closed a gap
found by l1-reviewer-package on 2026-11-13: 19 snorkel arms whose only fold had been marked superseded.
[`tag_refold_order_check.tsv`](../tag_refold_order_check.tsv) records, per ordered insert, the cassette
check, which fold carries it and its verdict. The 23 Igκ multipass references whose mature protein is
also a retired rescue arm's carry that fold, matched by sequence.

**IL1R2 (`l1-m1153.nterm_built`)**: on 2026-11-12 its model was refolded locally after being published,
so the table's `pdb_sha256` briefly described the refold. The published object (uploaded 2026-11-12
13:33Z) is the record; the local copy was restored from it and re-scored, and `publish_tagged_structures.py`
now stops if a local model differs from its published copy.

## Reproducibility

| What | Where |
|---|---|
| ColabFold, AlphaFold and JAX versions, SHA-256 of all 12 weight files, image registry digest | [`colabfold_provenance.json`](colabfold_provenance.json); the fold script pins the image by that digest (`ghcr.io/sokrypton/colabfold@sha256:22f1243a…`), not by the movable `1.6.2-cuda12` tag |
| Every AlphaFold DB model an s4e9 or s5e13 analysis reads: the 2,251 refold natives, every AFDB model in `tag_refold_native_models.tsv`, every `model_url` in `library_structure_payload.tsv`, the membrane-side gate's pinned models, the s5e13 corpus, and models other sessions request ([`tag_refold_afdb_requests.tsv`](../tag_refold_afdb_requests.tsv)) | `external://alphafold-db/v<N>/AF-<accession>-F1-model_v<N>.pdb`, written once and never overwritten. [`tag_refold_afdb_archive.tsv`](../tag_refold_afdb_archive.tsv) gives each model's AFDB URL, version, archive key, SHA-256, size, where the bytes came from and which readers use it (2,688 models, all v6; 70 of them are l1-reviewer-package's 2026-11-14 downloads, named in [`tag_refold_afdb_sources.tsv`](../tag_refold_afdb_sources.tsv); its independent downloads of 1,910 others match byte for byte). The refold scorer reads only these copies, checks the SHA-256, and stops if a model is missing; it never fetches AFDB live. [`tag_refold_afdb_unarchived.tsv`](../tag_refold_afdb_unarchived.tsv) lists the 73 models no local cache holds (54 payload, 19 native index); AFDB refuses this machine, so they are not archived yet |
| Licensed predictors (DeepTMHMM v1, SignalP 6.0) | [`licensed_predictors.json`](licensed_predictors.json): the Modal volumes that hold the bytes and their SHA-256s, hashed in place; the bytes are not copied anywhere |
| Every native MSA used | `processed://…/tagged/2026-11-12/msa/<native_sha256>.a3m.gz`; a tagged construct's MSA is that file re-indexed by `retarget()` in the fold script, deterministically |
| Every rung input as read | `…/tagged/2026-11-12/inputs/`, listed in `tag_refold_inputs.tsv`; request files are cumulative, so an input emptied upstream cannot drop a construct |
| Verdicts re-derived from the saved models | `score_tagged_folds.py --folds <dir> --check` re-scores into a temporary file and requires it to match `tag_refold_metrics.tsv` byte for byte |
| Published tables never shrink | `publish_tagged_structures.py` refuses to publish a verdict table that loses any (tagged SHA-256, source id) row a previous publish carried |

## Outputs

| File (s4e9 root) | One row per |
|---|---|
| [`tag_refold_queue.tsv`](../tag_refold_queue.tsv) | protein folded or to fold: tier, rung, SHA-256s, native↔tagged residue map, cassette span, provenance, attribution |
| [`tag_refold_natives.tsv`](../tag_refold_natives.tsv) | native protein sequence, by SHA-256 |
| [`tag_refold_sources.tsv`](../tag_refold_sources.tsv) | rung construct → tagged SHA-256 |
| [`tag_refold_metrics.tsv`](../tag_refold_metrics.tsv) | fold: every metric, verdict, reasons, settings |
| [`tag_refold_verdicts.tsv`](../tag_refold_verdicts.tsv) | rung construct: its id, model key, verdict and the metrics behind it |
| [`tag_refold_structures.tsv`](../tag_refold_structures.tsv) | published model: S3 keys, SHA-256s, cassette span, settings |
| [`tag_refold_viewer_manifest.tsv`](../tag_refold_viewer_manifest.tsv) | built construct, for the viewer: tagged and native model paths, MSA key, pLDDT, verdict and metrics |
| [`tag_refold_construct_check.tsv`](../tag_refold_construct_check.tsv) | built construct: translated-insert identity checks |
| [`tag_refold_recode_check.tsv`](../tag_refold_recode_check.tsv) | Twist-recoded construct: synonymous, found, protein unchanged |
| [`tag_refold_order_check.tsv`](../tag_refold_order_check.tsv) | insert in the order of record: cassette check, live fold, verdict |
| [`tag_refold_held_calibration.tsv`](../tag_refold_held_calibration.tsv) | held internal member: untagged control lDDT, its control-relative bar, best arm, arms passing under the rule |
| [`tag_refold_control_relative.tsv`](../tag_refold_control_relative.tsv) | tagged arm of a calibrated member: original verdict beside the control-relative one |
| [`figures/s4e9_tag_refold_lddt.png`](../figures/s4e9_tag_refold_lddt.png) | native-region and site lDDT by rung, beside the controls |
| [`figures/s4e9_tag_refold_epitopes.png`](../figures/s4e9_tag_refold_epitopes.png) | SpyTag003 vs ALFA exposure |

Models are under `processed://sprint4/s4e9-orf-design-fg-panel/predicted-structures/tagged/2026-11-12/`
(`<fold_id>.pdb`, pLDDT in the B-factor column; `<fold_id>.scores.json.gz` with PAE and pTM), with
`structures.tsv`, `verdicts.tsv`, `verdicts_<rung>.tsv` and `tagged_fold_metrics.tsv` beside them.

## Rerun

```bash
uv run --locked python library/l1/orf-design/tag_refold/native_models.py
```

```bash
uv run --locked python library/l1/orf-design/tag_refold/build_refold_inputs.py
```

```bash
uv run --locked --group modal modal run --detach scripts/fold_tagged_modal.py \
    --queue library/l1/orf-design/tag_refold_queue.tsv \
    --run tagged-2026-11-12 --tiers 1,2,3 --spawn --msa-only
```

```bash
uv run --locked --group modal modal run --detach scripts/fold_tagged_modal.py \
    --queue library/l1/orf-design/tag_refold_queue.tsv \
    --run tagged-2026-11-12 --tiers 1,2,3 --follow --budget-usd 200
```

```bash
uv run --locked --group modal modal run scripts/fold_tagged_modal.py --collect-only \
    --run tagged-2026-11-12 --out-dir ~/.cache/deliverome-analysis/tagged-folds/tagged-2026-11-12
```

```bash
uv run --locked python library/l1/orf-design/tag_refold/archive_afdb_models.py
```

```bash
uv run --locked python library/l1/orf-design/tag_refold/score_tagged_folds.py \
    --folds ~/.cache/deliverome-analysis/tagged-folds/tagged-2026-11-12
```

```bash
uv run --locked python library/l1/orf-design/tag_refold/plot_tagged_folds.py
```

```bash
uv run --locked python library/l1/orf-design/tag_refold/publish_tagged_structures.py \
    --folds ~/.cache/deliverome-analysis/tagged-folds/tagged-2026-11-12 --date 2026-11-12
```

```bash
uv run --locked python library/l1/orf-design/tag_refold/verify_ordered_constructs.py \
    --folds ~/.cache/deliverome-analysis/tagged-folds/tagged-2026-11-12
```

```bash
uv run --locked python library/l1/orf-design/tag_refold/score_tagged_folds.py \
    --folds ~/.cache/deliverome-analysis/tagged-folds/tagged-2026-11-12 --check
```

```bash
uv run --locked --group modal modal run scripts/fold_tagged_modal.py --collect-msas \
    --queue library/l1/orf-design/tag_refold_queue.tsv \
    --out-dir ~/.cache/deliverome-analysis/tagged-folds/msa
```

```bash
uv run --locked --group modal modal run scripts/fold_tagged_modal.py --provenance-only \
    --out-dir library/l1/orf-design/tag_refold/colabfold_provenance.json
```

```bash
uv run --locked --group modal modal run scripts/hash_predictor_volume.py --volume deeptmhmm-v1 \
    --out library/l1/orf-design/tag_refold/deeptmhmm_v1_volume_sha256.json
```

The builder reads its inputs from the S3 snapshots in `tag_refold_inputs.tsv`; a rebuild from them
is byte-identical. `--inputs-dir <folder>` ingests a new set of rung exports and snapshots them. The
control folds come from the benchmark run (`--run tagged-bench-2026-11-12`) and are copied into the
main folds folder before scoring.

## Decision log

| Date | Decision | Decided by | Implemented by |
|---|---|---|---|
| 2026-11-11 | Refold terminal constructs and the N-terminal site candidates; $230 of Modal for this work | Becca (request); coordinator (allocation) | Claude Code session l1-tag-refold |
| 2026-11-12 | One model, 3 recycles, no relax; reuse the native MSA with gap columns | coordinator (`OVERNIGHT-L1-2026-11-12.md`) | l1-tag-refold |
| 2026-11-12 | Verdict thresholds set from the 17 untagged controls; SP state reported, not a failure | l1-tag-refold (delegated by Becca 2026-11-12) | l1-tag-refold |
| 2026-11-12 | Inputs are the rungs' final post-#283 exports; one fold per tagged SHA-256 | l1-tag-refold, on the rung sessions' requests | l1-tag-refold |
| 2026-11-12 | Fold order by decision value (tiers above); tier 2 as 100-construct samples | l1-tag-refold (delegated by Becca 2026-11-12), following the coordinator's triage | l1-tag-refold |
| 2026-11-12 | Native models matched by sequence, reusing l1-reviewer-package's 66 ortholog folds; only natives with no exact model are refolded | coordinator (correction) | l1-tag-refold |
| 2026-11-12 | Type III TM1 fails only a helix the native model has and the tag loses (IGSF5's AFDB model leaves TM1 unmodelled) | l1-tag-refold (delegated by Becca 2026-11-12) | l1-tag-refold |
| 2026-11-12 | Request inputs cumulative; published verdict tables may only grow | l1-tag-refold, after a re-ingest dropped 52 loop rows | l1-tag-refold |
| 2026-11-13 | Folded sequences proven equal to the ordered (diversified-codon) inserts; MSAs, weights hashes and image digest saved; `--check` added | coordinator, for Becca ("save the structure with the diversified gggs but constant alfa and spy") | l1-tag-refold |
| 2026-11-13 | The order of record is a queue source; every ordered insert has a live verdict row | coordinator, after l1-reviewer-package found 19 superseded-only arms | l1-tag-refold |
| 2026-11-13 | IL1R2's published model is the record; the local refold was discarded | l1-tag-refold (delegated by Becca 2026-11-12) | l1-tag-refold |
| 2026-11-12 | Three MSA tickets of 20 in flight, server back-off honored; whatever has no MSA by the report is listed as not folded | l1-tag-refold (delegated by Becca 2026-11-12) | l1-tag-refold |
| 2026-11-14 | AFDB models archived write-once under `external://alphafold-db/v<N>/` (refold natives, then every s4e9/s5e13 reader) and the refold scorer reads only the archive; ColabFold image pinned by digest; licensed predictor bytes located and hashed in place, not copied | coordinator, after an independent reproducibility audit | l1-tag-refold |
| 2026-11-14 | Held internal members (22): an untagged control of each designed protein, folded and scored like the arms; where the control is below 0.85, the member's native-region bar is max(control − 0.05, 0.70), site and epitope rules unchanged; reported in a separate table, original verdicts kept (`calibrate_held_members.py`) | coordinator under Becca's delegation (rejected: lowering 0.85 for every member) | l1-tag-refold |
| 2026-11-15 | Rescue round (ranker, 74 internal ectodomain sites on held members; 116 first sent and folded), S5E13 order-rule arms (451 extra loops, 85 snorkels; approved to $26), 49 weak-SP Igk forms, SLC41A2 and OR10G6 residue-1 arms | coordinator for Becca ("find replacements where things were held"; "relying more on igk now for weak predicted SP") | l1-tag-refold |
