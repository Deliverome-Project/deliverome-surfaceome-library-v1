# L1 design review: cover note on work in progress

The [design brief](BRIEF.md) describes the L1 library as it stands. Six issues are still
being worked on before L1 is ordered. They are listed here, separately from the brief, so the
brief can be read as the design and this note can be updated as each one closes. Each row
gives its status as of **2026-11-04**. Comments on any of them are welcome; they are not the
questions we are asking reviewers to settle.

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

## Keeping this note current

When an issue closes, change its status to the date and the decision, and keep the row:
a reviewer who read an earlier version should be able to see what changed. If a fix changes
a library table, rerun [`brief_numbers.py`](brief_numbers.py), update any count above and in
the brief, rebuild the input bundle, and publish it as a new dated release of the review
repository.
