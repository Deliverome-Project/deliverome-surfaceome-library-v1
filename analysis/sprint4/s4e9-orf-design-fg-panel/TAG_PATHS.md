# Tag placement for the L1 library

> **2026-11-12:** RGMB's reference changed to the named form Q6NW40 (README, "Reference-form fixes"), whose DeepTMHMM v1 topology puts it on N after SP, so the GPI reroute below now touches BST2 only. The tables in this note are as measured before that change; current counts are in the reviewer brief, section 5.

> **2026-10-29 / 2026-11-01:** 7 L1 members whose protein duplicated another member's were retired, and the isoform and ortholog panels were redrawn by rule to their designed 150 and 200 (13 members added; see the README, "Panels are redrawn by rule, reproducibly"). `l1_tag_paths.tsv` places all 2,706 active members with a designed CDS; the three without a v1 topology (SLC43A2, SUSD3, cyno ECE1) were run on 2026-11-01 (`publish_topology_supplement.py`). C-terminal members are now 218 (SLC43A2 and cyno ECE1 added). **2026-11-02:** SLC67A1 was added as `l1-p2739` (internal rung), so 2,707 members are placed and internal-site members are 567. The per-rung totals below predate these changes.

Where the tag cassette goes for each of the 2,700 L1 members, decided from the topology of
the protein that member actually encodes.

## What makes this different from the earlier analyses

The unit is the **member**, not the canonical gene. We order 2,700 L1 members — 2,360
library references plus isoforms, orthologs and rescue arms — and a member's topology is
not its gene's. SLC9D1's ordered isoform has 3 TM helices where its canonical has 13;
TRPV1's has 3 against 6; KCNN2's has 1 against 6. A rung assigned from the canonical is
the wrong rung for those constructs.

Each member is matched to its topology by the **sha256 of its ordered protein sequence**,
through `deliverome_analysis.topology`. Accession lookups cannot do this: an
Ensembl-transcript ORF frequently corresponds to a UniProt isoform record, so the canonical
accession returns a different protein's topology.

Topology is DeepTMHMM **v1** throughout — the string the structure viewer colours by, and
the topology of record (v1 has 99.5% signal-peptide recall against v2's 92.9%).

## The ladder

Unchanged from `scratch/becca/tag-site-architecture/derive_tag_paths.py`. Changing the unit
without changing the rule is what makes the two comparable.

    ecto-N (at residue 1, or after the signal peptide) → ecto-C → internal

| Placement | Members | Share | Rests on |
|---|---:|---:|---|
| N-term, after signal peptide | 1,295 | 48.0% | predicted cleavage site |
| internal loop or snorkel | **566** | 21.0% | per-member topology |
| N-term at residue 1, multipass | 563 | 20.9% | Tedman GPCR screen |
| C-terminal | 216 | 8.0% | topology only |
| N-term at residue 1, single-pass (type III) | 60 | 2.2% | topogenesis, untested; per-subgroup design in [TYPE3_NTERM.md](TYPE3_NTERM.md) |

By candidate set:

| Placement | library_reference | isoform | ortholog |
|---|---:|---:|---:|
| N-term, after signal peptide | 1,121 | 74 | 100 |
| internal loop or snorkel | 517 | 23 | 26 |
| N-term at residue 1, multipass | 523 | 13 | 27 |
| C-terminal | 164 | 20 | 32 |
| N-term at residue 1, single-pass (type III) | 35 | 14 | 11 |

**The `library_reference` column reproduces the canonical-gene analysis** — 1,121 / 517 /
523 / 164 / 35 here against 1,129 / 515 / 525 / 166 / 32 there. The residual differences
are members whose ordered proteoform is not the canonical one, which is the entire point of
redoing it per member. That agreement is the check that the ladder was not altered.

## N-terminal tag at residue 1 on multipass members

The 563 members on "N at residue 1, multipass" are 495 GPCRs (489 in Tedman's screen) and 68
non-GPCRs. Whether they take an Igκ leader in front of the cassette, and which members need
different handling whatever is chosen, is worked through in
[NTERM_RES1_MULTIPASS.md](NTERM_RES1_MULTIPASS.md), with the per-member table in
`nterm_res1_multipass.tsv`. GPCRs take the Igκ leader by default (Becca, 2026-11-10); a
literature review recommends Igκ for 63 of the 68 non-GPCRs, pending her confirmation.

## Internal-site worklist

**566 members** need a loop insertion, against the 515 **canonical genes** the earlier
analysis reported. Those numbers measure different things and neither supersedes the other
arithmetically: 515 is genes, 566 is members, and 517 of the 566 are library references.

### Reconciling against the 515

The two counts differ by eight genes, and every one has a cause. Restricted to
`library_reference`, my internal set is 517 genes against the 515 table's 515:

    515 + 5 − 3 = 517

| | Gene | Why it differs |
|---|---|---|
| **+** | BST2 | GPI reroute. The 515 **excluded** GPI proteins as a stated known gap; here they are rerouted. |
| **+** | RGMB | GPI reroute, same reason. |
| **+** | EPHA6 | Its ordered proteoform is not the canonical one (folded as `library_designed_orf`, 2 TM). |
| **+** | SERINC5 | Ordered form is a UniProt isoform, 10 TM. |
| **+** | SLC17A2 | Ordered form is a UniProt isoform, 12 TM. |
| **−** | AQP12B | Ordered form has an **extracellular N-terminus**, so it has a terminal option the canonical did not. |
| **−** | LSR | Ordered form carries a 41-aa signal peptide → N-after-SP. |
| **−** | TSPAN10 | Ordered form has an **extracellular C-terminus** → ecto-C. |

So the 515 is not wrong. It is the right answer to "which canonical genes can only take a
loop insertion", and it reproduces here to within the two GPI proteins it deliberately set
aside and six genes whose ordered proteoform is not their canonical. The member-level answer
is **566**, and that is the number to plan constructs against, because the three genes in
the lower half of that table would otherwise be built with a loop insertion they do not
need, and the five in the upper half with a terminal tag that cannot be displayed.

Only **7 of the 566** have no inter-TM extracellular run of ≥ 8 aa — and for 2 of those 7
that 0 is an artefact rather than a finding.

### Inter-TM loop vs ectodomain: two different questions

`longest_inter_tm_ecto_run_aa` counts only runs with a TM on **both** sides. That is the
right measure for a multipass loop inventory, and it returns 0 for any single-pass protein
by construction — so reading a 0 as "nowhere to put a tag" is wrong for the GPI-rerouted
members, whose site goes somewhere in the ectodomain rather than between two helices.
`longest_ecto_run_aa` is the answer for those, with `longest_ecto_run_position` saying
whether the run sits before the first TM, between two, after the last, or in a protein with
no TM at all.

| Member | Gene | TM | inter-TM | longest ectodomain | Where |
|---|---|---:|---:|---:|---|
| `l1-p1731` | RGMB | 1 | 0 | **394 aa** from 85 | c_terminal |
| `l1-p0188` | BST2 | 1 | 0 | **137 aa** from 44 | c_terminal |
| `l1-p2196` | TMEM268 | 2 | 7 | 7 aa | inter_tm |
| `l1-p2202` | TMEM86A | 8 | 7 | 7 aa | inter_tm |
| `l1-p2203` | TMEM86B | 8 | 7 | 7 aa | inter_tm |
| `l1-p2172` | TMEM139 | 2 | 6 | 6 aa | inter_tm |
| `l1-p2178` | TMEM169 | 4 | 1 | 1 aa | inter_tm |

The split is clean, and the ectodomain measure confirms the inter-TM one rather than
contradicting it: for the five genuine multipass members the two agree, because their only
extracellular region **is** that short loop. Those five are the real no-site cases. The two
GPI members have 137 and 394 residues to work with. L2 behaves the same way — GPC1 539 aa
and SPRN 133 aa against TMEM268 7, SLC28A1 4, TMEM134 1. An N-terminal
ectodomain is not an inter-TM loop and is not counted — mistaking one for the other is the
error recorded in `scratch/becca/appris-principal/SHEDDING_AND_TAGS.md`, where ADGRB1's
918-residue N-terminal segment was once read as a loop candidate.

## GPI anchors change the C-terminal rung

A GPI omega site **deletes the C-terminus**, so an ecto-C tag on a GPI-anchored protein
leaves with the propeptide. `internal-tag-targets/README.md` names this as a known gap and
excludes GPI proteins; here they are rerouted to internal instead of dropped.

191 L1 members (147 genes) are GPI-anchored, but the reroute touches only **2** of them:

| Rung before the GPI check | Members |
|---|---:|
| N after SP | 182 |
| N at residue 1 | 7 |
| C-terminal → rerouted to internal | **2** (BST2, RGMB) |

That is not a small effect hiding a bug. A GPI protein with a usable extracellular
N-terminus never needed its C-terminus, so the omega site does not enter its placement
decision. Only a GPI protein that would otherwise have been tagged at the C-terminus is
affected. FOLR1, named alongside BST2 in the earlier note, lands on N-after-SP and so needs
no reroute.

The GPI set itself is built by `build_gpi_set.py` → `gpi_set.tsv` and is **157 genes**, not
the 130 a curated-annotation-only pull returns. See below.

## The C-terminal rung is checked member by member before it is built

A C-terminal rung call says the C-terminus is extracellular. It does not say the C-terminus
is still there at the cell surface. `review_l1_cterm_members.py` checks each of the 212
members (`l1_cterm_review.tsv`, member set of PR #283 at 13f3c1b8). It marks a member only when
the tag would be cut off, and records every other risk as a flag. One member, human ITM2B
(`l1-m0355`), loses a C-terminal tag (`cterm_tag = removed`): furin releases the C-terminal
Bri23 peptide (residues 244–266), which is the GPI rule above applied to a propeptide.
`build_l1_cterm_arms.py` builds the other 211 at the C-terminus, and builds ITM2B with the
cassette after residue 235, before the furin site (Becca, 2026-11-10). All 212 are
built.
The README section "C-terminal members carry their cassette" gives the flag counts, the
checks on each insert and the questions still open. The most relevant flags here are the
v1/v2 topology checks. Only one member's C-terminus is cytoplasmic by v2: the TRPV1 isoform
`l1-m1559`, which also has a 1-residue v1 extracellular tail. v2 calls no TM helix for three
members (BRICD5, TCTA, TNFSF13B) and has no row for three more (the new orthologs CLEC2L,
CORIN and SLC47A1).

## Every L1 member is anchored

A member with neither a transmembrane helix nor a GPI anchor has nothing holding it in the
membrane, so no tag site can display it — the protein is secreted.

| | Members |
|---|---:|
| L1 members | 2,700 |
| zero TM helices | 173 |
| &nbsp;&nbsp;of those, GPI-anchored | 173 |
| &nbsp;&nbsp;**neither TM nor GPI** | **0** |

MICA is the gene this test was written for. Its reference is `l1-p1216`, carrying the
383-aa form: a 22-aa signal peptide, one transmembrane helix at residues 307–328, and a
55-residue cytoplasmic tail.

The form the selection ladder originally picked, `ENST00000449934`, is 332 aa and keeps no
helix at all. The two are the same protein through `WQTFHVSAVAA` and then diverge, so what
the short form loses is the second half of the helix together with the `KKK` cluster that
sets its orientation:

| Form | Length | TM by v1 | TM by v2 | C-terminus from residue 307 |
|---|---:|---:|---:|---|
| `ENST00000400325` — the reference | 383 | 1 (res 307–328, 22 aa) | 1 (res 309–329, 21 aa) | `WQTFHVSAVAAAAIFVIIIFYVRCCKKK…` |
| `ENST00000449934` — superseded | 332 | 0 | 1 (res 315–324, 10 aa) | `WQTFHVSAVAAGCCYFCYYYFLCPLL` |

The decision does not rest on choosing between the two predictors. v1, the topology of
record, calls the short form `SP` — signal peptide and nothing else. v2 dissents, but what
it offers is a 10-residue segment, roughly half the length needed to cross a bilayer, with
8 cytoplasmic residues after it. Both tools agree the 383-aa form carries a full helix, so
it is the sound construct on either reading. MICA's shed risk is recorded **high** —
soluble MICA is a known immune-evasion route — which makes a secreted form the worse of the
two to build, not merely the less useful one.

Repointing cost no synthesis slot and changed no identifier. The replacement sequence was
already designed and domesticated as the isoform arm `iso:MICA:Q29983`, so the correction
was applied to the `lib:MICA` key in place; that arm, now an exact duplicate of the
reference, is retired. The rung is unchanged — both forms are `N after SP` — so every count
on this page is the same either side of the change. One number does improve: the reference
now matches a published v1 row (`topo_2026_05_16`) instead of the provisional
`library_designed_orf` run.

### Why the existing guard did not catch MICA

`scratch/becca/appris-principal/select_reference.py` already carries a soluble guard —
`SURFACE_LABELS = {TM, SP+TM, BETA}` against `NON_SURFACE_LABELS = {SP, GLOB}` — added
after PILRA and LMLN were both moved onto a soluble form by taking APPRIS unconditionally.
It is deliberately conservative: *"a form with no topology call is not rejected, only a form
called SP or GLOB that carries no GPI anchor."*

That conservatism is exactly why MICA passed. The guard reads
`deeptmhmm_human_canonical.tsv` and `deeptmhmm_human_isoforms.tsv` — the **UniProt**
proteoform cohorts. A library member's ORF is an Ensembl-transcript product, so for members
whose form is in neither cohort there was no call to reject. MICA's L1 reference is one of
those: its topology exists only as of the 2026-10-14 `library_designed_orf` run.

| Gene | Ordered form TM | Topology cohort | Outcome |
|---|---:|---|---|
| PILRA | 1 | `human_canonical` | guard had a call, works |
| LMLN | 0 | `library_designed_orf` | no call at selection time — saved only by being GPI |
| MICA | 0 | `library_designed_orf` | no call at selection time, and not GPI |

The guard was right and the topology was missing. Two things close that loop.
`scripts/check-construct-anchoring.py` re-applies the same test to **all 5,887** classified
members rather than only those that happen to be UniProt proteoforms, and CI runs it with
no exemptions. Upstream, MICA now names its reference outright in
`appris-principal/reference_overrides.tsv`, so a future run of the ladder picks Q29983
rather than re-deriving the secreted form — the same mechanism CD55 uses, where APPRIS had
likewise ranked a transcript encoding the secreted form first.

### Eleven L2 orthologs are unanchored for a different reason

Eleven L2 members still have neither a TM helix nor a GPI anchor: ASGR1, CLEC5A, CRB2,
ICAM1, ICAM5, LRFN4, LRRN4CL, NRXN1, PCDHB3, TP53I13 as cyno orthologs, and NAALADL1 as a
mouse one. They are not MICA's case, and the difference decides what can be done about
them.

**None has a TM-bearing form of its own species in the library.** Every alternative the
gene has is human or mouse, so swapping one in would substitute a different organism's
protein rather than a different proteoform of the same one — the ortholog arm would stop
being an ortholog arm. `check-construct-anchoring.py` keys its suggestion on
(gene, species) for this reason and reports these as not repointable.

Nor is there an obvious form to re-pick. In each species' reference proteome, 8 of the 11
genes have exactly one annotated proteoform — the one already chosen. The three with
alternatives (CLEC5A, CRB2, NRXN1) offer only shorter ones, which would have to gain a
helix the longest form lacks.

What the evidence points to instead is incomplete annotation outside human. Nine of the
eleven are shorter than their human reference, by 30 to 674 residues:

| Gene | Ortholog | Human reference | Difference |
|---|---:|---:|---:|
| NRXN1 | 833 | 1,507 | −674 |
| LRFN4 | 380 | 635 | −255 |
| CRB2 | 1,136 | 1,285 | −149 |
| ASGR1 | 183 | 291 | −108 |
| TP53I13 | 314 | 393 | −79 |
| ICAM5 | 863 | 924 | −61 |
| ICAM1 | 480 | 532 | −52 |
| PCDHB3 | 762 | 796 | −34 |
| LRRN4CL | 208 | 238 | −30 |
| CLEC5A | 188 | 188 | 0 |
| NAALADL1 | 745 | 740 | +5 |

A gene model missing its 3′ end is consistent with losing the exon that encodes the helix,
which is where every one of these nine would carry it. CLEC5A and NAALADL1 are the same
length as their human counterparts and so need their own explanation; neither has been
looked at residue by residue yet.

Resolving these is a design round, not a repoint: each needs a TM-bearing sequence sourced
from somewhere other than the current annotation, and then ORF design, domestication and a
topology run. Until that happens the honest options are to source those sequences or to
drop the eleven arms.

CI gates on L1 only, since that is the tranche we order from. These eleven are printed and
do not block, and will block once L2 is ordered (`--gate l1,l2`).

## Shedding is recorded and does not change the rung

695 genes carry a recorded shed form. The flag travels on every row
(`shed_risk`, `shed_evidence`) because an extracellular terminus on a shed ectodomain
leaves with the fragment.

It does not reroute anything, because there is nowhere to route to:
`SHEDDING_AND_TAGS.md` established that of the 122 members whose N-terminal tag is at
proteolytic risk, **zero** have a usable inter-TM loop — 108 are single-pass with a
one-piece ectodomain, 11 have no annotated TM span, and 3 are multi-pass with a largest
inter-TM loop under 30 aa. Rerouting them would invent an option that does not exist.

## Reproducing

```bash
uv run --locked python build_gpi_set.py        # -> gpi_set.tsv
uv run --locked python l1_tag_paths.py         # -> l1_tag_paths.tsv
```

`l1_tag_paths.py --library l2` runs the same analysis on L2, where 1,136 members still lack
a v1 topology for their exact sequence and are reported rather than classified.

| File | Role |
|---|---|
| `l1_tag_paths.py` | the analysis |
| `l1_tag_paths.tsv` | one row per member: placement, TM count, SP length, longest inter-TM loop, shedding, topology provenance |
| `build_gpi_set.py` | the GPI set, with scores |
| `gpi_set.tsv` | 159 genes scored, 157 anchored, with basis and likelihoods |
| `gpi_anchored_genes.tsv` | curated input (130 genes) — an **input**, not the answer |
| `shed_form_risk.tsv` | 695 genes with a recorded shed form |
| `repoint_mica_reference.py` | MICA's reference correction; idempotent, and refuses a replacement with no v1 helix |
| `build_library_topology_index.py` | per member, which topology row describes its ordered sequence |
| `library_topology_index.tsv` | that pointer: cohort, version, accession, ordered length, sequence hash |
| `build_library_structure_payload.py` | per member, which structure to render and which to hold back |
| `library_structure_payload.tsv` | one row per (member, structure), `is_chosen` on the one to show |
| `compute_residue_offsets.py` / `residue_offsets.tsv` | where the ordered protein sits inside a `contained` model |
| `build_library_manifest.py` | the merged per-member table, and the length assertions against the CDS |
| `library_manifest_l1.tsv` | **what a consumer should read**: one row per L1 member, identity + `topology_*` + `structure_*` |
| `../../../scripts/check-construct-anchoring.py` | the CI gate: every ordered member has a helix or an omega site |

### One table still cannot be re-run

`check_proteoform_match.py` treats any failure from the AlphaFold API as "no model" and
drops the member from its output, reporting only a count. Each run loses a different 25 to
26 members that way, so consecutive runs on identical inputs disagree about which rows
exist. Retrying a failed accession, and distinguishing a transient error from a genuine
absence, would make the table reproducible. Until then `proteoform_match.tsv` is edited in
place for a single member rather than regenerated.

`build_library_structure_payload.py` was in the same position and no longer is. It marked a
`contained` model usable only when its offset was zero, reading the placeholder `0` that
`proteoform_match.tsv` ships on every contained row as "no shift needed". It now reads the
measured offsets from `residue_offsets.tsv`, so all 26 contained members are chosen with
the span to render, and the script reproduces the committed table. A member whose offset was
never located is held back rather than rendered at a guessed one.

### A reproduction gap in the older GPI path

`scratch/becca/appris-principal/gpi.py::gpi_genes()` cannot be run: its three curated
inputs — `gpi_supplement.tsv`, `surface_subcategory.tsv`, `gpi_exclusions.tsv` — are
gitignored and are neither tracked nor present on disk, and it also reads a
`~/.cache` UniProt keyword pull. `build_gpi_set.py` deliberately avoids all of them,
taking NetGPI's committed raw output plus a curated list exported from
`triage_run_public`, so the set here rebuilds from tracked files alone. If `gpi.py` is
ever needed again, those three inputs have to be recovered or regenerated first.

Inputs not in this folder: the published topology snapshot
(`topology.SNAPSHOT_KEY` on `deliverome-external`), the designed CDS
(`processed://sprint4/s4e9-orf-design-fg-panel/library-candidates/candidate_cds.json`), and
NetGPI's raw output in `scratch/becca/appris-principal/netgpi_gpi_human_predictions.txt`.

## Decision log

**2026-11-10** — "N at residue 1, multipass" characterized per member in
`NTERM_RES1_MULTIPASS.md`: 495 GPCRs and 68 non-GPCRs, sorted into five subgroups by what the
Igκ-leader decision turns on. Two MFSD2A members are probably mis-placed (structure puts the
C-terminus, and so the N-terminus, inside). Placements in `l1_tag_paths.tsv` are unchanged.

**2026-10-15** — Redone per member rather than per canonical gene, on each member's own
ordered-proteoform v1 topology. Internal-site worklist is 566 members (517 of them library
references, reproducing the earlier 515 genes).

**2026-10-15** — GPI set rebuilt from NetGPI 1.1 over reference *and* alternative forms,
unioned with curation under the 0.55 confidence floor: 157 genes, up from 130. The first
pass of this applied the confidence floor to the strongest negative across all forms, which
read an alternative form's confident negative as the reference's and flipped ART4 (reference
0.534, alternatives 0.991) and FCGR3B (0.445 against 0.994) out of the set. The floor was
calibrated on reference sequences, so it is the reference form's negative that it tests.

**2026-10-15** — Shedding kept as a flag rather than a reroute, on the finding that no
shedding-risk member has a usable inter-TM loop.

**2026-10-15** — Fixed `library_design.gene_from_key`. The manifest key does not put the
gene in a constant field: `ort:cyno:ABCB11:D2CFL3` has the **species** in field 1, where
`lib:`, `iso:` and `rescue:` keys have the gene. Taking field 1 unconditionally returned
"cyno"/"mouse" for all 196 L1 ortholog members, which matched nothing in any gene-keyed
join, so GPI and shedding annotation was silently absent for every one of them. The
zero-TM audit is what surfaced it: it reported 26 members as neither-TM-nor-GPI, of which
25 were orthologs whose GPI status simply could not be looked up. The real count is 1.
Placements did not change — no ortholog was both GPI and ecto-C — but shed flags went from
600 to 637. Tests now pin every key shape.

**2026-10-15** — MICA's reference repointed from the 332-aa `ENST00000449934` to the 383-aa
`ENST00000400325`, and the CI exemption removed. The selection ladder is not at fault for
the original pick: it ranks transcripts by how well sourced they are, which cannot see
whether the encoded protein keeps its anchor. Where the ladder's answer and displayability
conflict, displayability wins, so the reference is named outright in
`reference_overrides.tsv` rather than left to be re-derived. Applied in place on the
`lib:MICA` key by `repoint_mica_reference.py`, which keeps every candidate id stable and
retires the duplicate `iso:MICA:Q29983` arm; safe to do because `tag_cassette` was
`deferred` and `barcode` `not reserved` on both rows, so no id was bound to a reagent.
Nothing was built or ordered under either form. L1 members remain 2,700 and no placement
count changes; L2 drops to 3,187.

**2026-10-15** — Corrected the claim that the 11 unanchored L2 orthologs each have a
TM-bearing form already in the library. They do not: every alternative is a different
species, and substituting one would change the organism rather than the proteoform. Eight
of the eleven have only one annotated proteoform in their species, so there is nothing to
re-pick, and nine are shorter than their human reference. `check-construct-anchoring.py`
now keys its suggestion on (gene, species); while it was gene-keyed it offered the human
member as available for all eleven, which is the substitution it exists to prevent.

**2026-10-15** — Fixed two bugs in `library_topology_index.tsv`, reported by the
extracellular-site-prediction session, which consumes it alongside
`library_structure_payload.tsv` and found the two disagreeing on species for 39 members.

The first was a species-blind match. `topology.resolve` matches on the sha256 of the
sequence, and a human protein is frequently byte-identical to its cyno or mouse ortholog,
so an unconstrained match could return the other organism's row. Measured against the
manifest rather than against the other file, **66 of 5,888** members were matched across
species (43 human→cyno, 15 cyno→human, 5 mouse→human, 2 human→mouse, 1 mouse→cyno) — more
than the cross-file check found, because that check cannot see a case where both files
inherited the same error. `topology.SPECIES_COHORTS` now constrains the match, and
`normalize_species` folds `cynomolgus`/`cyno` so a groupby cannot split one organism in
two.

The second was column clobbering. `topology_report` wrote the member's identity first and
then spread the matched topology row over it, and both sides have `species` and
`gene_symbol` — so a matched row replaced the member's with the matched protein's. That is
why the index read `cynomolgus` and `G7PPA4_MACFA` for a human library reference, and why
1,877 of 5,888 `gene_symbol` values were UniProt entry names. The matched row's own values
are now kept under `topology_*`, so the collision is unrepresentable rather than
order-dependent.

**The topology strings were never wrong.** A cross-species match only happened where the
sequences were byte-identical — that is what made the hashes collide — so the per-residue
string was the same either way. Regenerating changed **zero placements and zero TM counts**
across all 2,700 L1 members. What was wrong was the provenance: species, accession, entry
name and cohort. Eight members lost a match that should never have counted, so L2 goes
3,188 → 3,180 classified and the L2 v1 gap 1,129 → 1,136.

Both fixes are defect-seeded, and the first version of each test did **not** catch its
bug: the species test passed because the correct row won the version tiebreak anyway, and
the clobbering test passed because the `topology_*` rename makes spread order irrelevant.
The fixtures were strengthened until removing each fix fails a test.

**2026-10-15** — Added `longest_ecto_run_aa` / `_start` / `_position`, because the inter-TM
measure returns 0 for single-pass members by construction and that 0 was being read as "no
site available". It is not: BST2 has 137 ectodomain residues and RGMB 394.

**2026-10-15** — `INTERNAL_SITE_ONLY.md`, `internal_site_only.tsv` and
`internal_site_undetermined.tsv` deleted from `scratch/becca/appris-principal/`. They are
fully superseded by this page and `l1_tag_paths.tsv`, which carries per-member loop bounds,
and they were actively misleading: their 377/176 split reads UniProt's missing
topological-domain annotations as unknown orientation. The record is in git history at
commit `8eab499`, where they still carry the superseded notice explaining the arithmetic.

**2026-10-15** — Took the ordered protein's length off the sequence rather than off a codon
count, after the extracellular-site-prediction session found
`library_structure_payload.tsv` and `library_topology_index.tsv` disagreeing on it for four
members. The payload derived it as `orf_nt // 3 - 1`, which assumes every designed ORF ends
in a terminator codon. **Thirteen of 7,016 do not** — three in L1 (KIR2DL5A, KIR2DS2,
KIR2DS5) and ten in L2 — so those shipped a length one residue short of the protein. The
length is what pairs the structure row with the topology string, and an off-by-one pairing
puts every tag site one residue out against the model. `orf_aa` in `proteoform_match.tsv`
is measured on the translated CDS and agreed with the designed sequence for all 2,629
members that have both, so the formula was the only thing wrong.

The fourth member was MICA, and it was not a defect in either table but the two branches
being read together: the repoint had reached the index, which resolves through the designed
CDS on S3, and not the payload, which read the old `orf_nt` from the committed manifest.
Both are now on one branch and agree at 383 aa. It is worth naming the hazard rather than
only the instance — `designed_proteins` joins a **committed** manifest to an **S3** CDS
file, so overwriting that file silently changes what every committed derived table means,
with no version to compare. Today exactly one member's `orf_nt` disagrees with the CDS
length, and that was MICA before the merge.

Two guards added to `tests/test_library_structure_coverage.py`, both defect-seeded: the
payload's length must equal the measured `orf_aa`, and a member's chosen structure and its
topology string must agree on the length. Reinstating the formula fails both; moving only
the index's length fails only the second, which is the MICA shape. `build_library_topology_index.py`
also now writes LF — `csv` defaults to CRLF, so the index rewrote every line on each run
and git normalised it back, which is how a real diff would have hidden in a phantom one.

**2026-10-15** — Merged the topology index and the structure payload into one table,
`library_manifest_l1.tsv`, because two independent sessions consuming them made the same
two mistakes doing the join themselves — and both are mistakes the shape of the two files
invites. Both carry `species` and `gene_symbol`, so whose they are depends on which file
you read them from, and after the join every row says "human" whether or not it is: that
is how a human member's structure got paired with its cyno ortholog's topology invisibly.
And a cross-file comparison cannot see an error both files inherited, so had both read the
stale MICA reference, every check comparing them would have passed.

The merged table puts the matched topology under `topology_*` and the chosen structure
under `structure_*`, which makes a species conflict unrepresentable rather than
order-dependent, and it anchors every length to the **designed CDS** rather than to either
file. The CDS is the arbiter: a repoint that reaches one table and not the other is caught
without the two needing to disagree. 2,700 L1 members, 0 without a v1 topology, 0 matched
across species, 1 with no chosen structure (l1-p2256 TPO, which has no AlphaFold model at
all — a real absence, not a join artefact).

**The model-length rule is per match class, and this is the part a flat assertion gets
wrong.** `model_length_aa == ordered_length_aa` holds for `identical` (2,404) and
`ours_by_construction` (269) and fails for **all 26** `contained` rows, because the model
is longer there — that is what "contained" means. For those the invariant is the rendered
span, `display_last - display_first + 1 == ordered_length_aa`, which holds 26 of 26. A
flat comparison would have raised 26 false alarms on its first run, which is how a new
guard acquires an exemption list and then stops being read.

L2 is deliberately not in the table. `library_structure_payload.tsv` is L1 only, and a
blank structure column cannot distinguish "this table does not cover you" from "no model
exists" — so the 4,340 L2 members are counted and named instead.

Four guards, all defect-seeded: a duplicated row, a topology row from the wrong organism,
an `identical` model one residue off, and a `contained` span three residues short each fail
a test. The negative control matters as much — on the correct table all 26 contained models
differ in length from the ordered protein, so the flat rule really would have fired.

**2026-10-15** — Documented a wart in `topology_acc` rather than renaming it, after the
tag-sites-enrichment session hit it. It is the topology snapshot's join key, and for the
`library_designed_orf` cohort that key is the **candidate id** (`l1-p1037`), not a UniProt
accession — those predictions were published keyed on the member, because no UniProt record
exists for the form. 158 of 2,700 L1 rows. Renaming would break the join back to the
snapshot, so the column keeps the snapshot's meaning and `topology_cohort` is what says
which kind of key it holds; `build_library_manifest.py` now prints the count on every run.

The failure it causes is worth naming because it is silent in the worst way. A consumer
resolving `topology_acc` inside `cohort='human_canonical'` finds nothing for those 158, and
if it reads the miss as "junction out of range" it drops real sites while reporting
success. That happened to KCNN2, which has two L1 members — `l1-p1037` (791 aa, 6 TM, the
`library_designed_orf` form) and `l1-p2458` (a 231-aa single-pass isoform) — and whose two
literature sites, between residues 246/247 (S3–S4 loop) and 162/163 (S1–S2 loop), sit on
the 791-aa form. Resolve by `sequence_sha256`, which is on every row for this reason.

Also recorded from that thread, because it bears on how any per-site measure should be
scored: `site_kind` in `tag_site_public` does **not** discriminate terminal from internal
at low positions. The five literature sites at p 4–5 (ABCC1 R4, ABCC2 K4, SLC5A1 S4,
SLC5A2 T5, SLC5A4 T4) are topologically identical — no signal peptide, extracellular
N-terminus, junction inside the N-terminal extracellular run — yet two are labelled
`internal` and three `terminal_n`, with both paralogue pairs disagreeing with their own
partner. There is no defined boundary between the two for a tag a few residues into an
extracellular tail. A loop insertion should be derived from topology — an extracellular run
with a TM helix on **both** sides — which is `longest_inter_tm_extracellular_run`'s rule
applied per site instead of per member.

**2026-10-15** — Folded the species spelling at the source in `library_manifest_l1.tsv`,
after the extracellular-site-prediction session found the table still spelling the macaque
both ways: `species` is `cyno` for all 99 members while `topology_species` was `cynomolgus`
for 92 of them and `cyno` for 7. The cross-species count was right either way — there are
genuinely zero — but a consumer writing the obvious `species != topology_species` got **92
false positives**, which defeats the point of a table built to make a species conflict
visible. It is the same `cynomolgus`/`cyno` fold the matcher already applies internally,
surviving as a label in the output. Both columns now go through
`topology.normalize_species`, so plain equality is the correct check.

The test was the other half of the problem and is the more useful lesson. It normalised
both sides before comparing, so it passed while the table shipped the split — a test that
normalises hides exactly the defect the column exists to expose. It now compares with plain
equality, the check a consumer would write, and seeding the old spelling back into 99 rows
fails it.

**2026-10-15** — `gene_symbol` is **not unique** in any of these tables, and a
`{gene_symbol: member}` dict is lossy rather than redundant. Measured on
`l1_internal_members.tsv` (566 rows) after the tag-sites-enrichment session hit it:

| key | distinct | >1 member | disagree on the longest inter-TM loop |
|---|---:|---:|---:|
| `gene_symbol` | 518 | 45 | 14 |
| `(species, gene_symbol)` | 543 | 22 | **7** |

Nine of the 14 gene-keyed disagreements are a human protein against its mouse or cyno
ortholog, which *should* differ. The seven that matter are same-species
`library_reference` vs `isoform` pairs: LIM2 40 vs 82 aa, TMEM219 161 vs 160, KCNK4 93 vs
90, PLPP1 52 vs 48, SLC12A4 114 vs 110, SLC15A2 210 vs 206, and SLC4A7 83 vs 76 — the last
with 14 TM against 12, a difference large enough that the two members are not the same
design problem.

**Gene-keying also misclassifies, not just overcounts.** LIM2 and TMEM219 each have three
members spanning human and cyno, so a gene-keyed group reads as cross-species and its
genuine human-reference-vs-human-isoform conflict is attributed to the species difference.
LIM2's 40-vs-82 aa — both human, both 4 TM — is the sharpest same-species case in the set
and the one a gene-keyed analysis explains away. Key on `(species, gene_symbol)`, or better
on `candidate_id`, and reach for `sequence_sha256` when the question is about the protein.

Independently, that session's own implementation of the both-sides rule agrees with
`longest_inter_tm_ecto_run_aa` for **566 of 566** members — two implementations of the
rule matching on every row, which is a better check on it than this analysis had alone.

### Open question: LIM2's reference proteoform

**For Becca, not resolved here.** LIM2 is the one gene where the ortholog's length sides
with the human *isoform* rather than with our library reference, and the difference lands
squarely on the region a tag has to go in:

| member | species | set | length | longest inter-TM loop |
|---|---|---|---:|---:|
| `l1-p1099` | human | library_reference | 173 aa | **40 aa** |
| `l1-p2461` | human | isoform (P55344-2) | 215 aa | **82 aa** |
| `l1-p2577` | cyno | ortholog (G7PYH7) | 215 aa | 82 aa |

LIM2 is an internal-placement member either way, so the loop is the whole site inventory,
and 40 aa against 82 aa is a 2x difference in the only place the cassette can sit.

**The cyno evidence is weak on its own** and should not decide this. G7PYH7 is a TrEMBL
entry, and *M. fascicularis* gene models are frequently automated, so a length agreeing
with the human isoform may reflect the gene model rather than the dominant transcript. What
would settle it is human evidence for which form is expressed in the tissue we care about —
the selection ladder picked P55344 on MANE+APPRIS agreement, which ranks how well sourced a
transcript is and cannot see that one form has half the usable loop.

Scope of the check: of the 13 L1 genes that have a human reference, a human isoform **and**
a non-human ortholog, only LIM2 shows a substantive split. TMEM219's is 1 aa (initiator
methionine) and GYPC's ortholog is closer to the isoform without matching it. So this is one
gene to look at, not a systematic problem with the ladder.

## Designated MIN_LOOP_LEN floor probes (not built — no slot)

The floor at 10 aa is the largest single structural uncertainty in the library designs:
**185 extracellular regions of 7–9 residues sit on 140 of the 515** internal-placement
proteins, and if the floor is wrong those gain a candidate region. MS4A1 P79 (`p0032`)
probes it, but at a tight side of **0** and a clearance of **1** it tests an extreme —
only 2.8% of 176 published internal sites are that tight — so a null there is
attributable to the extreme rather than to length.

Two better probes, identified by the tag-location session and verified here against the
**L1 designed proteoform** (not the canonical) before being recorded. Neither is built;
both are ready if a slot appears.

| | member | cohort / acc | probe | run | clearance | control |
|---|---|---|---|---|---:|---|
| **TMEM26** | `l1-p2194` | `human_canonical` / Q6ZUK4 | **res 279** (F) | 277–284, **8 aa** | **3** | **res 105** (R), 86–142, 57 aa, clearance 20 |
| **TMEM268** | `l1-p2196` | `human_canonical` / Q5VZI3 | **res 126** (A) | 125–131, **7 aa** | 2 | *none — this is its only extracellular run* |

Both are `alphafold_db` / `proteoform_match = identical`, so the model scored is the
protein we would build and no pLDDT or RSA feature is imputed.

**TMEM26 is the better floor test.** It reproduces the MS4A1-162-against-79 structure with
each of that pair's weaknesses removed: the probe is 8 residues rather than 6, so it sits
where the floor actually bites instead of at the extreme; clearance 3 rather than 1 means a
null is not confounded with "too close to the membrane"; and the control at absolute 0.797
clears the median published region outright where MS4A1 162 at 0.552 does not. Its second
8-aa run (28–35, best residue 35) has clearance **1** and reintroduces the confound, so it
is not the first choice.

**TMEM268 answers a different question** — whether a protein with no region above the floor
is taggable at all. Its 7 residues are its *entire* extracellular surface, and its ceiling
of 0.693 is level with the median published region (0.695), so the model rates the stub as
well as a typical region where a tag demonstrably worked. Note its clearance is **2**, so
on that axis it sits between MS4A1 79 and TMEM26 279; its value is coverage, not clearance.
One imputed feature: `dist_to_glyc` is NaN across all seven residues because the protein
carries no sequon anywhere.

**DISP2 was considered and is out**, on both arms. The canonical A7MBM2 that scores well is
our **L2 isoform** at 1,401 aa; the L1 construct `l1-p0513` is 1,184 aa, a
`library_designed_orf` with no UniProt record, and residue 1015 there is **intracellular** —
a compartment flip, not an offset. Re-deriving on the L1 sequence would not rescue it
either: that member has `structure_feature_coverage 0.0`, so a score on its three 8-aa runs
would rest entirely on imputed pLDDT and RSA, which are the two features the model is built
from. The candidate generator now requires `human_canonical` on the same accession **and**
un-imputed structure features, which excludes 22 of the 185 regions across 19 genes.

**The rule this establishes:** score a floor-probe candidate on the L1 designed proteoform
matched by `sequence_sha256`, never on the canonical. 158 of the 2,700 L1 members are
`library_designed_orf` with no UniProt record and 821 resolve under `human_isoforms`, so for
a substantial minority a canonical-keyed score describes a different protein — and DISP2 is
the case where that flips a compartment call rather than shifting a number.
