
## The `hold` column in barcodes.tsv was never enforced (2026-10-18)

p0026 carried `hold = "HOLD 2026-10-18: ..."` from the moment its tagged SignalP call came
back as OTHER. No code read that column. `rebuild_order_maps.py` rebuilt the map and
published it to `processed://constructs/designed/p0026.gb` on every run, and
`export_order_batch.py` wrote it into the order folder, exactly as if nothing were flagged.
The column was documentation that read like a gate.

It is renamed `note`, which is what it is. A hold that must actually hold has to be read by
the publisher and the exporter, and the way to know it works is to set it on a construct and
watch that construct disappear from both outputs -- not to see the string in the file.

## Two things the re-prediction verdict does not tell you (2026-10-18)

Both found by the InterPro session's review of its own run, both confirmed here. Until the
upstream fields land, a reader of `comparison.verdict` has to recompute these.

**`preserved` carries no confidence, so a softening arm looks identical to a safe one.**
SLC9A3 reports `verdict: preserved, ok: True` while its cleavage-site probability falls
0.967 -> 0.650 and its p(SP) falls 0.867 -> 0.667 — the largest degradation in the set,
and the only arm that started below 0.9. The probabilities are on the records but not in
the verdict block, and the native's values live on a *different* record, so seeing this at
all requires a manual join across the two. The canary console prints p(SP) but not cleavage
probability, so it does not surface it either. **SLC9A3 is the arm with least margin in
this batch and nothing in the verdict says so.**

**`cleavage_moved` is not graded against the junction, so benign and fatal share a string.**
HEG1's 28 -> 29 at junction 32 is benign: the mature chain goes from 4 native residues
ahead of the cassette to 3. The fatal case is a shift to **>= the junction**, which would
put the cassette inside the signal peptide. HEG1 has 3 residues of margin. Confirmed by
arithmetic here rather than inherited: for any `cleavage_moved` arm, compare the new
cleavage site against the junction, not against the old site.

**And the two VSIG10 results are not equally sourced.** Junction 31 was run on the real
construct (`constructs/designed/p0026.gb`); junction 32 is a sequence assembled for that
test, residue-exact (sha256 `0df6fab8bb253e86`) but never built or registered. "The built
construct fails and the one alternative its interval allows also fails when modelled" is
the claim; "both positions tested" would overstate the second.

## The screen is UNTESTED, not refuted — and the batch never depended on it either way

SignalP reports two things and they fail independently. `sp_probability` says whether a
signal peptide exists; `cleavage_probability` says how confidently it knows where it ends.
Our two problem arms sit on **different** axes, and neither is low on the other:

| gene | p(SP) | cleavage-site p | |
|---|---|---|---|
| VSIG10 | 0.9998 | **0.348** | certain there IS one, no idea where it cuts |
| SLC9A3 | **0.867** | 0.967 | less certain there is one, confident where it cuts |
| AXL, ERBB2, FOLR1, HEG1, NECTIN4, PTPRC | ≥0.9997 | 0.964–0.982 | clean on both |

**NOT refuted — underpowered, and on the wrong population.** I first recorded this as
refuted on the strength of "three of the six losses had confident natives". That is an
anecdote, not a rate, and I committed it without asking for the denominators — the exact
error this work keeps turning up. Recomputed from the run itself (67 native/tagged pairs,
5 losses: ADGRD1, FZD6, FZD8, RXFP1, TSHR):

| native column, split at 0.9 | weak lost | strong lost | odds ratio | Fisher one-sided |
|---|---|---|---|---|
| `sp_probability` | 1/5 = 20.0% | 4/62 = 6.5% | 3.6 | p = 0.33 |
| `cleavage_probability` | 2/13 = 15.4% | 3/54 = 5.6% | 3.1 | p = 0.25 |

**The direction favours the screen on both columns.** Weak natives lost their signal peptide
about three times as often. Neither reaches significance — but 5 events over 67 pairs cannot
establish *or* refute a 3x effect, and "did not reach significance" is not "does not
predict". The honest status is **untested**.

**And the population is wrong for us anyway.** All 67 are GPCRs — Tedman's screen is
GPCR-only, and the five that lost are adhesion (ADGRD1), frizzled (FZD6/FZD8), relaxin
(RXFP1) and glycoprotein-hormone (TSHR) receptors. **Not one of our eight native-SP arms is
a GPCR**: HEG1, ERBB2, AXL, NECTIN4 and VSIG10 are single-pass type I, FOLR1 is GPI-anchored,
SLC9A3 an 11-TM SLC, PTPRC a receptor phosphatase. Transferring a GPCR result onto that set
is the scope error we have been catching all day, run in the other direction.

**What does not change either way: the batch never rested on the screen.** All eight arms
were re-predicted *with the cassette in place*; the six are fine because their **tagged**
predictions preserved the cleavage site, not because their natives looked clean. That was
true when I thought the screen worked and it is true now. The two-column completeness
argument was redundant, and a redundant argument is dangerous precisely because nothing
fails when it is wrong.

**And no screen is needed:** the full 7,101-sequence sweep is about **$4.87** of GPU. The
constraint is wall time unsharded, not money. Ranking which members to test solves a problem
that does not exist — test them all.

**A third verdict hazard, flagged by the InterPro session (2026-10-19):**
`cleavage_crosses_junction` returns `False` when the signal peptide is lost entirely, so a
filter on `crosses == False` collects VSIG10 as safe. Being changed to `None` upstream.
Nothing here depends on it — the numbers in this folder were recorded directly — but any
new consumer should treat "did not cross" and "there is nothing left to cross" as different
answers.

## "Move into the room you have" is not a general rule — AXL tests it and fails (2026-10-19)

SLC9A3's fix invited the obvious generalisation: four of our eight arms sit at `sp_end + 1`,
so move them all. Tested on the two that can move inside their own admissible interval
(`canary_rp_1792375661`), and the answer is no.

| gene | interval | junction | cleavage-site p |
|---|---|---|---|
| AXL | 32–36 | 33 (built) | 0.816 |
| AXL | 32–36 | 34 / 35 / 36 | 0.817 / 0.817 / 0.817 |
| NECTIN4 | 31–33 | 32 (built) | 0.973 |
| NECTIN4 | 31–33 | 33 | 0.979 |

**AXL is flat to three decimal places across its whole interval.** Its 0.978 → 0.816 is not
positional, and it is already at its best available position. NECTIN4 gains 0.006, which is
noise.

The reason is the window, and it is the same mechanism that explains SLC9A3 — read in the
other direction. What matters is how many **native mature residues** sit in SignalP's
70-residue window ahead of the cassette. SLC9A3's interval (25–43) lets it buy nine; AXL's
(32–36) caps it at four. The remedy needs *enough* room, not merely *some* room, and the
amount available is set by where the fold starts.

So the honest scope: of the four arms at +1, one was marginal (SLC9A3 at 0.650) and one is
improvable (SLC9A3, to 0.973). AXL at 0.816 and NECTIN4 at 0.973 are comfortable and
unimprovable by moving. This bounds the SLC9A3 result to the case it came from rather than
letting it become a rule nobody tested.

**Also: junction 46 was never admissible.** SLC9A3's interval is 25–43, so of the positions
swept, 34 and 40 are inside it and **46 is outside**. 46 scored marginally highest (0.974 vs
0.973) and was rejected on structural grounds before this was known — the right answer for a
reason that did not include the real one. Only 34 and 40 are selectable; 34 remains the
recommendation.

## My window mechanism was wrong — HEG1 falsifies it (2026-10-19)

I explained SLC9A3's recovery as SignalP's 70-residue window: at junction 26 it holds one
native mature residue and 44 of cassette, at 34 it holds nine, so the remedy is "buy more
native context" and it needs *enough* room rather than *some*. AXL being flat across its
narrow interval looked like confirmation.

**It was confounded, and the InterPro session spotted it.** Across those two arms, residue
count and *which axis was degraded* co-varied perfectly. SLC9A3 was the only arm with a
native p(SP) deficit (0.867) and its move recovered both axes. AXL's p(SP) was 1.000 native
and tagged everywhere — there was no deficit to recover, only a site, and the site did not
move. "Needs ≥9 native residues" and "only a p(SP) deficit is recoverable" predict both
results identically.

**HEG1 separates them and I lose.** No p(SP) deficit (0.99971 native), a site softened like
AXL's, and ample room:

| junction | native mature residues ahead of the cassette | cleavage-site p |
|---|---|---|
| native | — | 0.972 |
| 32 (built) | 4 | 0.840 |
| 40 | 12 | **0.789** |
| 48 | 20 | 0.859 |

Flat and **non-monotonic**. Twenty native residues is more than double what SLC9A3 bought,
and it recovers nothing — it is worse at 12 than at 4. The cleavage also stays shifted
28 → 29 at every junction.

**What the three arms actually support:** moving the junction recovers a **p(SP) deficit**,
and does not recover a softened **cleavage site**. SLC9A3 was the only one of the eight with
a p(SP) deficit, and the only one that improved. Residue count was incidental.

The SLC9A3 move stands — 0.650 → 0.973 and 0.667 → 1.000 are measurements on that construct,
not inferences. Only my explanation, and anything extrapolated from it, were wrong.

## FOLR1's cassette sits inside its predicted fold, and the SP check cannot see it

Flagged by the InterPro session's census and confirmed here against `disorder_public`.
FOLR1's admissible interval is **22–24**; p0028 is built at **27**, three residues past it.
Per-residue disorder at residue 27:

| predictor | score at 27 | calls it |
|---|---|---|
| NetSurfP-3.0 | 0.859 | disordered |
| metapredict | 0.369 | **ordered** |
| AlphaFold-disorder | 0.255 | **ordered** |

Median of three is 0.369, i.e. **ordered** — and metapredict and AlphaFold-disorder both put
the fold's start at residue 23, immediately after the signal peptide. FOLR1's SignalP
re-prediction came through untouched (0.982 → 0.981), which is the point: **the re-prediction
is blind to fold insertion.** A clean SP result says nothing about whether the cassette
landed in structure.

Not resolved, and not obviously fixable. InterPro puts FOLR1's first domain at 36, which
disagrees with the disorder consensus; and FOLR1 is GPI-anchored, so the N-terminus is the
only available end — there is no alternative terminus to move to, and the admissible
positions 23–24 are closer to the cleavage site, not further into open sequence. Recorded as
a known risk on p0028 rather than a pending change.

## HEG1's null result is also the only noise calibration we have (2026-10-19)

Moving HEG1 across its interval with nothing else in play varied its cleavage-site
probability over **0.789–0.859**, a spread of **0.070**. That is the only measurement
anyone has of how much this number moves for positional reasons alone, and it is the
yardstick a small drop has to be read against. Applying it to the eight arms:

| gene | native → tagged | drop | |
|---|---|---|---|
| SLC9A3 | 0.967 → 0.650 | **0.317** | clearly softened |
| AXL | 0.978 → 0.816 | **0.163** | clearly softened |
| HEG1 | 0.972 → 0.840 | **0.132** | clearly softened |
| ERBB2 | 0.965 → 0.905 | 0.060 | **inside the band — not distinguishable from noise** |
| NECTIN4 | 0.980 → 0.973 | 0.007 | inside the band |
| PTPRC | 0.981 → 0.979 | 0.002 | inside the band |
| FOLR1 | 0.982 → 0.981 | 0.001 | inside the band |
| VSIG10 | 0.348 → lost | — | signal peptide abolished |

**The 0.070 cutoff is not doing any work.** Sorted, the drops are 0.317, 0.163, 0.132,
then 0.060, 0.007, 0.002, 0.001 — an **empty interval 0.072 wide** between 0.060 and 0.132,
with the band sitting inside it. Any threshold from 0.061 to 0.131 returns the same three
arms. A cutoff that separates a clean gap is a different kind of claim from one that splits
a continuum, and only the second needs defending.

**But the band is one protein's positional variance, not a noise floor.** 0.070 is how much
HEG1's site probability moves when only its junction changes. Applying it to another protein
assumes comparable positional variance, which nothing has tested. It is the only calibration
available and it is fine for this batch, where the gap makes the exact value irrelevant — but
a future member landing at 0.08–0.09 would be **ungraded**, not passing, and the honest
answer there is to sweep a clean-p(SP) arm and widen the basis rather than lean on this one.

**Three arms clearly softened, not four.** ERBB2 was being counted as a softening and it
is not one — that was reading noise as signal, and nobody could have known until HEG1 gave
us a band. Spotted by the InterPro session, which struck it from its own "three of four
arms at +1 softened" line for the same reason.

## The library census, and why its shapes no longer imply a remedy

Earlier here: "low on both axes: none; clean on both: six" — true of the eight built arms
and the reason the batch has no hidden third case. **It does not generalise.** Across the
701 placeable library members, cut at 0.9 on both columns (InterPro session's census,
2026-10-19):

| | members | |
|---|---|---|
| clean on both | 539 | 77% |
| low p(SP) only | **7** | the SLC9A3 shape — the only shape with a positional remedy |
| low cleavage-site only | **110** | the AXL/HEG1 shape — no positional remedy |
| **low on BOTH** | **45** | **no representative among the eight built arms, nothing ever measured** |

Two things follow. The positional remedy, if real at all, reaches seven members — our one
success is not a general tool. And the 45 double-low members are the worst-looking set in
the cohort on the only two columns anyone has, with no measurement on any of them; VSIG10
was low on one axis and lost its signal peptide outright, so that is where the next one is
most likely to be.

## Barcode registry: 22 of 24 appended; p0034 and p0035 blocked on a stale P-number

Appended 2026-10-19 with Becca's go-ahead. The shared s1e2 registry now holds 36 rows and
carries 22 of the batch's 24 barcodes. `append_reservations.py` verified every appended row
is present exactly once.

**Two are held back, and the obstacle is a label rather than a barcode.** All 24 barcodes
are themselves free in the registry — none collides. But the registry still carries rows
from before the 2026-10-18 renumbering that attach those two P-numbers to different
constructs:

| registry row | barcode | where that barcode actually is now |
|---|---|---|
| `p0034` / TFRC / idx 21 | `ATCACGCGACGAGAGAAGCT` | **p0047**, AVPR2 leader-free |
| `p0035` / TFRC / idx 22 | `TGTCTCTCTCAGCTTCTCTC` | **p0037**, TFRC C-term |

So those barcodes are correctly reserved and nobody else will draw them; only the P-number
attribution is stale. Appending our `p0034` (CCR5 F263) and `p0035` (AVPR2 L302) would put
two rows under each number with different barcodes and genes, and a later reader asking
"what is p0034?" would get two answers.

The registry is append-only by convention, so the stale rows cannot be deleted. Three ways
out, and the choice is Becca's because it changes a shared team record:

1. **Append anyway with a comment line** recording that the earlier rows are superseded.
   Keeps our P-numbers, leaves the ambiguity visible rather than resolved.
2. **Renumber** CCR5 F263 and AVPR2 L302 to never-used numbers. Clean registry, but both
   are already registered in the Reagent Inventory as P-0034 and P-0035, which is the
   registration authority — so this makes the authority wrong instead.
3. **Leave the two unappended** and rely on the Reagent Inventory, which already records the
   current assignment correctly. Costs nothing today; the risk is that the barcodes stay
   invisible to anyone drawing from the pool against the registry alone.

Recommended: (1). The Reagent Inventory already says what p0034 and p0035 are, so the
registry only needs to not contradict it, and a supersession comment does that without
touching bytes anyone else wrote.
