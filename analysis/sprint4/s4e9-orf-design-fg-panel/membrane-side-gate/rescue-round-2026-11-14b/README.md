# Membrane-side gate on the rescue-round sites (2026-11-14b)

**Question:** is each rescue-round tag site placed outside the membrane on the member's own model?
The rescue round proposes sites on held members that are not on the internal rung: mostly
single-pass and GPI-anchored ectodomains, plus six multipass members with a signal peptide.

**Inputs:** `processed://sprint5/s5e13-tag-site-benchmarking/rescue_round/2026-11-14b/`, with the
SHA-256 of each file in `rescue_round_gate_summary.json`. That is 74 sites on 39 members, the list
after the 23 members built on the N-after-SP rung (PR #353) were dropped. Produced by
`../gate_rescue_round.py`.

**Topology:** UniProt's. The member's extracellular residues come from the ranker's per-residue
table. Transmembrane and intramembrane spans come from the span-reference accession, mapped
onto the member by alignment (identical residues only). None of these members has a DeepTMHMM
v1 row, so no verdict here is a v1 call.

**Placement and the meaning of a verdict:** the `basis` column says which rule applies.

| basis | when | verdict |
|---|---|---|
| `multipass` | 2+ TMs; TMDET (TmAlphaFold model within 2 Å of ours) or the helix-axis fit | the gate's calibrated verdict, thresholds from `../calibration.json` |
| `single_pass_tmdet` | 1 TM and TmAlphaFold's model is within 2 Å of ours | at or above +14.0 Å is a pass; lower is `review_single_pass_geometry`, never a fail |
| `not_placeable` | GPI anchor (a lipid the model lacks), or 1 TM with no usable TMDET placement | no height; `reason` says why |

A single-pass site is never failed on its height. The thresholds were calibrated on multipass
proteins, and a single-pass ectodomain's angle to the membrane is not fixed by the model: a
stalk or hinge lets it swing. A low height can show where this model put the domain, but not
that the site sits in the membrane. Decided by the l1-membrane-side-gate session (Claude Code)
for the rescue round, under Becca's 2026-11-14 request.

`seq_distance_to_membrane` counts residues from the junction to the nearest TM or
intramembrane residue, or to the ω site for GPI members. Use it for the not-placeable sites.

**Result:**

| basis | pass | review | not placeable |
|---|---|---|---|
| multipass (GABRA2, OR10G6, SLC39A6, GRIA1, TRPV1) | 9 | 0 | 0 |
| single-pass, TMDET | 15 | 3 (`review_single_pass_geometry`) | 0 |
| not placeable | – | – | 47 |

There are no fails.

**The three reviews:**
- KREMEN2 after 362: −11.4 Å, 3 residues before the TM. This is the juxtamembrane stalk, so the
  low height is plausible rather than a pose artifact.
- MUC21 after 440 (−24.7 Å) and after 465 (−42.7 Å): both junctions are 15–40 residues from
  the TM, and the model places them below the mid-plane at pLDDT 58 and 31. That is this
  model's pose, not evidence that the sites are in the membrane.

**VTCN1 (l1-p2342 human, l1-p2628 cyno):** not placeable under any of the three readings of its
N terminus.
- UniProt TM 260–280 (predicted): one helix, and TmAlphaFold has no entry.
- GPI: there is no helix at all.
- DeepTMHMM signal anchor 7–29: still one helix.

The sites' distance from the membrane depends on which reading is right:

| Reading | After 40 | After 187 |
|---|---|---|
| C-terminal anchor (TM 260–280, or GPI near there) | 220 residues | 73 residues |
| N-terminal signal anchor (7–29) | 11 residues | 158 residues |

Only site 40 changes, from 220 residues away to 11.
