# Panel redraw of 2026-11-22

Draft README/brief section. Status: **applied to the manifest**
(`apply_panel_redraw_2026_11_22.py`, change log `panel_redraw_2026-11-22_changes.tsv`).
The panel tables of record (`isoform_panel.tsv`, `ortholog_panel_adopted.tsv`, the
`panel_inputs.json` pins) and the rung tables are not yet regenerated, and the new members
have no built arms yet.

## Isoform panel

The isoform panel is a seeded, class-stratified draw of 150 forms: 100 from the observed tier
(forms with peptide or Ribo-seq evidence) and 50 from the prediction-only tier. It uses seed
20261009 and floors of type III 2, GPI 5 and type II 5. The remaining slots in each tier go
to classes by largest remainder. Each gene's class comes from its selected reference protein,
using DeepTMHMM v1 matched by sequence.

The redraw changes two things:

| Change | What it does | Why |
|---|---|---|
| Rebuilt cohort and gate | Draws from `isoform_cohort_rebuild.tsv` (1,899 forms) gated by `isoform_gate_rebuild.tsv` (1,690 pass) | Always rebuild with the latest inputs (Becca, 2026-11-22) |
| Multipass split by N-terminal side | `multipass N-out` (reference N-terminus extracellular; 672 library genes) and `multipass N-in` (cytoplasmic; 551 genes) | An N-out and an N-in multipass protein take the tag on opposite sides of the membrane, so each orientation gets its own quota |

A form is replaced only when it loses a feature that decides display: the signal peptide, the
N-terminal side, the first TM helix or the membrane anchor. The isoform gate already tests
exactly these. It also requires a native ATG start and an exact Ensembl 115 CDS, because both
decide whether the form can be built. The rebuilt cohort already excludes forms that lose the
membrane or GPI anchor. So the narrower rule changes no isoform verdict.

As in `redraw_panels.py`, a pick is replaced within its stratum if no designed CDS encodes it or
if its protein repeats an active L1 member. Two picks were replaced this way: NCAM1
ENST00000401611 and KCNE3 ENST00000532569, both for having no designed CDS.

| Stratum | Observed | Prediction-only | Panel |
|---|---:|---:|---:|
| type I | 43 | 18 | 61 |
| multipass N-in | 21 | 11 | 32 |
| multipass N-out | 13 | 6 | 19 |
| type II | 10 | 7 | 17 |
| GPI | 10 | 6 | 16 |
| type III | 3 | 2 | 5 |
| **All** | **100** | **50** | **150 forms, 133 genes** |

Against the panel of record, 69 picks are kept, 81 are dropped and 81 are added:

| Class in the panel of record | Kept | Dropped |
|---|---:|---:|
| type I | 47 | 14 |
| multipass | 5 | 46 |
| type II | 8 | 9 |
| GPI | 6 | 10 |
| type III | 3 | 2 |

Most of the churn is in multipass, whose 51 slots are now two separate draws. The other
strata change too, because all strata draw in sequence from one seeded generator: any change
upstream shifts every later stratum's sample. No dropped pick fails the rebuilt gate. Each is
dropped only because the draw changed.

All 81 added forms were L2 members. Each moved to L1 under a new id, l1-p3215 to l1-p3295 in
table order, with its L2 row retired `moved_to:<new id>`. Four of them (ART3, EFNA1, PRSS8 and
SLC43A2) once held an L1 id that was retired `dropped_from_panel`. As in `redraw_panels.py`,
those old ids stay retired. The 81 dropped L1 members are retired
`dropped_from_panel:isoform_redraw_2026-11-22 (Becca 2026-11-22)`. Every added key is already
in the designed-CDS and construct records.

**Topology.** Every added form has a DeepTMHMM v1 row in the snapshot the gate read
(`2026-11-20-isoform-cohort-rebuild`). That snapshot is the gate's evidence. For 40 of the added
forms, the row (`topo_2026_11_20_isoform_cohort`) is missing from the shared default snapshot
(`2026-11-21-coreceptor-additions`). Those 40 rows are copied unchanged, under the
fresh pair `(topo_2026_11_22_isoform_redraw, human_isoforms)`, into
`external://accessible-surfaceome/topology/2026-11-22-isoform-redraw/topology.parquet`
(71,706 + 40 = 71,746 rows; SHA-256 `bfd4aeff821d49bd53455a9ae6b312a42b6b979acf1e94d75aa90ce810c8b9f5`;
`publish_isoform_redraw_topology.py`). Until `topology.SNAPSHOT_KEY` points at it,
`topology.for_sequence` with the default snapshot returns `None` for those 40.

## Ortholog feature rule

The correspondence review first replaced a drawn ortholog form if it lost any human feature
(signal peptide, any TM helix, N-terminal side, GPI anchor, or a Pfam domain across an indel).
Eleven of the 200 drawn forms were replaced. The rule now counts only the display-deciding
features listed above. This reuses the review's evidence unchanged.

| Drawn form | Lost under the earlier rule | Display feature lost | Under the narrow rule |
|---|---|---|---|
| cyno IL31RA | SP, Pfam | signal peptide | still replaced (CD84) |
| cyno IL6R | SP | signal peptide | still replaced (AMIGO3) |
| cyno ITGA10 | SP | signal peptide | still replaced (TNFRSF10B) |
| cyno NFASC | SP, N-side | signal peptide, N-terminal side | still replaced (NCR1) |
| cyno MS4A10 | TM1, N-side, Pfam | first TM helix, N-terminal side | still replaced (CCR5) |
| cyno ADGRB2 | one TSP_1 repeat | none | held on MRGPRX3 (l1-p2932): no CDS record or member encodes it yet |
| cyno FLT4 | Ig domain across an insertion | none | restored as l1-p3296; CSF2RB (l1-p2941) retired |
| cyno SLC8A1 | TM6–TM10 and the exchanger domain (C-terminal truncation) | none | held on MRGPRX4 (l1-p2970), pending Becca's decision |
| cyno TMEM67 | TM2 | none | restored as l1-p3297; GPR78 (l1-p2572) retired |
| mouse ADGRE5 | one EGF_CA domain | none | restored as l1-p3298; AVPR1B (l1-p3018) retired |
| mouse BTLA | Ig domain | none | restored as l1-p3299; IL15RA (l1-p3022) retired |

The ortholog panel diff is 196 kept, 4 dropped and 4 added. Each restored form is promoted
from the L2 member that already carries its exact protein and CDS (l2-p2115, l2-p2920,
l2-p3043, l2-p3126). Retired replacements carry `replaced_by_original:ortholog_rule_2026-11-22`.

Mouse OR10G6 is the one form kept on adoption that the narrow rule alone would replace, because
its human 21-residue N-terminus is called a signal peptide. Its recorded keep stands. The rule
does not override a recorded decision.

## Open questions

| Question | Why it matters |
|---|---|
| Restore cyno SLC8A1? | The narrow rule keeps it because its SP, N-side, TM1 and anchor are intact. It still lacks half its helices and the exchanger domain, so it displays but is not the human protein's fold. |
| Restore cyno ADGRB2? | The narrow rule keeps it, but no designed CDS or construct record encodes the form; restoring it needs a new CDS member published first. |

A per-stratum seeding scenario (`--per-stratum-seed`, each stratum drawn from
`random.Random("20261009|<tier>|<class>")`) was computed for comparison only. It keeps 21 of
the 150 record picks (129 change), against 69 for the shared stream, because every stratum,
not only multipass, gets a new stream. It had no unbuildable swaps.

## Provenance

| Element | Source |
|---|---|
| Isoform redraw | `redraw_isoform_panel_2026_11_22.py` → `isoform_panel_redraw_2026-11-22.tsv`, `_summary.json`; `--check` replays it |
| Ortholog rule | `ortholog_rule_2026_11_22.py` → `ortholog_rule_2026-11-22.tsv`, `_summary.json`; `--check` replays it |
| Manifest changes | `apply_panel_redraw_2026_11_22.py` → `panel_redraw_2026-11-22_changes.tsv`; `--check` fails unless applied |
| Topology merge | `publish_isoform_redraw_topology.py`; `--check` compares with the published key |
| Decision | Becca Carlson, 2026-11-22 |
| Implemented by | Claude Code session l1-review-export-scope (worker fork), 2026-11-22 |

## brief_numbers-ready counts

```
isoform_redraw_cohort_forms            1899
isoform_redraw_gate_pass               1690
isoform_redraw_panel_forms             150
isoform_redraw_panel_genes             133
isoform_redraw_multipass_n_out_genes   672
isoform_redraw_multipass_n_in_genes    551
isoform_redraw_multipass_n_out_picks   19   (13 observed, 6 prediction-only)
isoform_redraw_multipass_n_in_picks    32   (21 observed, 11 prediction-only)
isoform_redraw_kept                    69
isoform_redraw_dropped                 81
isoform_redraw_added                   81
isoform_redraw_replaced_in_draw        27
isoform_redraw_unbuildable_rejections  2
isoform_redraw_added_missing_shared_v1 40   (all 81 present in the gate snapshot)
ortholog_rule_replaced_earlier         11
ortholog_rule_still_replaced           5
ortholog_rule_narrow_rule_keeps        6
ortholog_rule_restored                 4
ortholog_rule_held                     2    (cyno SLC8A1 pending decision; cyno ADGRB2 no CDS record)
ortholog_rule_panel_kept               196
ortholog_rule_kept_by_decision_only    1    (mouse OR10G6)
```
