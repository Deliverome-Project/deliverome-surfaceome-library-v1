# Where the cassette goes after a signal peptide

**The rule: anchor on disorder, not on the cleavage site.** For the 1,337 library members
whose only accessible terminus is an N-terminus behind a cleaved signal peptide, place the
cassette at the most disordered residue in a short window downstream of the predicted
cleavage site — not immediately at the cleavage site, which is what this repo currently
builds.

Terms used below: the **signal peptide (SP)** is the N-terminal targeting sequence the
signal peptidase removes; the **cleavage site** is the bond it cuts, after which the
**mature N-terminus** begins; the **junction** is the bond our cassette is inserted into,
written `after N` meaning between residues N and N+1 in UniProt canonical numbering.

## What we build today, and why it is wrong

`build_constructs.py` computes `cut = sp_end * 3` and splices the cassette there, so the
mature protein begins with the cassette. That was recorded as "Tedman's own rule for his
signal-peptide-bearing receptors". It is not his rule.

Tedman's 946 screen-validated HA insertions are in
[`tedman_gpcr_controls.tsv`](https://github.com/Deliverome-Project/accessible-surfaceome/blob/dev/data/tag_sites/tedman_gpcr_controls.tsv)
with the junction of each construct projected onto canonical UniProt numbering. Of those,
834 are bare N-terminal insertions before residue 1 — the GPCRs with no signal peptide —
and 112 are inserted after a residue. Restricting to the 108 whose accession has a
DeepTMHMM v1 signal peptide:

| offset of the junction from the SP cleavage site | constructs |
|---|---|
| at the cleavage site (0) | **0** |
| downstream (> 0) | **107** |
| inside the SP (< 0) | 1 |

Median offset **+4**, range −1 to +19. Not one construct does what we do.

## The junctions sit in disordered regions

Scored at the junction residue against all three predictors published to
`disorder_public`, over the 112 inserted-after-N constructs:

| predictor | score at the junction | percentile within that protein | above 0.5 |
|---|---|---|---|
| metapredict v3 | 0.61 | **88th** | 74% |
| NetSurfP-3.0 | 0.93 | **92nd** | 94% |
| AlphaFold-disorder (1 − pLDDT) | 0.63 | **93rd** | 86% |

Three methods that disagree about absolute disorder agree about rank: the junction lands
in roughly the top tenth most disordered residues of each receptor. A tag spliced into a
flexible stretch neither displaces structure nor gets buried, which is the mechanical
reason to prefer it over the cleavage site, where the mature N-terminus may begin a folded
domain immediately.

## Why anchoring on disorder matters more than getting the cleavage site right

The cleavage site itself is uncertain at the same scale as the offset. Over the 65
accessions in this set where both predictors call a signal peptide:

- identical cleavage site: **45**
- differ: **20**, by 1–4 residues, with SignalP 6 usually calling it later than
  DeepTMHMM v1 (e.g. O75084 28 vs 32, P32241 28 vs 30, Q03431 26 vs 22)

A rule written as "cleavage site + 4" inherits that uncertainty directly, and a tag placed
4 residues after a cleavage site that was called 4 residues early ends up at the cleavage
site after all. A rule written as "the most disordered residue in the window downstream"
does not: disorder is a property of the sequence, so the same residue wins whichever
predictor sets the window's start.

This is the same uncertainty that makes the 1,337 the library's largest single dependency.
Anchoring on disorder converts it from a correctness risk into a window-width choice.

## The rule to implement

1. Take the SP end as the **later** of the DeepTMHMM v1 and SignalP 6 calls. Erring late
   is the safe direction: a tag placed inside the SP leaves with it at cleavage and the
   member reads as a non-expressor, indistinguishable from a real biological negative,
   whereas a tag placed a few residues late is still presented.
2. Search residues `sp_end + 1` through `sp_end + 20`, the span Tedman's constructs
   occupy.
3. Score each candidate by the mean of the three predictors' disorder, each converted to
   its within-protein percentile so the methods are comparable, and take the maximum.
4. Record the chosen offset and its disorder percentile per gene, so a placement can be
   audited against this table rather than re-derived.

Window width is the one free parameter. 20 covers the observed range; a narrower window
would force some genes onto a less disordered residue, and a wider one leaves the region
the signal peptide's own geometry governs.

## Open

- Tedman's offsets are measured against *our* predictors' cleavage calls, not against his.
  He used TOPCONS, which we have not run. The +4 median is therefore our reading of his
  geometry and could shift by a residue or two under his own SP calls.
- Whether the rule should apply to the 382 multipass members with no signal peptide is
  untested. Those take a bare N-terminal tag in his hands (834 of 946 constructs) and that
  arm is already validated, so there is no reason to change it.

## Decisions

| date | decision |
|---|---|
| 2026-10-08 | Post-SP placement anchors on disorder in a +1..+20 window, not on the cleavage site (Becca, from Tedman's actual junctions). |
| 2026-10-08 | SP end is the later of DeepTMHMM v1 and SignalP 6, because the failure modes are asymmetric. |
| 2026-10-08 | DeepTMHMM v1 is the topology of record; our own benchmark puts v1 at 99.5% SP recall against v2's 92.9%. |
