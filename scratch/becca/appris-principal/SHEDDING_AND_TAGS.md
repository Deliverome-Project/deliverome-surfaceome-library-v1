# Can a loop or juxtamembrane tag rescue the shedding-risk genes?

Mostly no, and the reason is topology rather than design effort.

## There is no loop to use

Of the 122 members whose N-terminal tag is at risk from proteolysis:

| | n |
|---|---|
| single-pass — the ectodomain is one piece, no internal loop | 108 |
| no transmembrane span annotated | 11 |
| multi-pass, largest **inter-TM** loop < 30 aa | 3 |
| **multi-pass with a usable inter-TM loop** | **0** |

*Corrected.* An earlier pass counted ADGRB1 as a loop candidate on a 918-residue
extracellular segment. That segment is its N-terminal ectodomain, which sits
before the first transmembrane span, not an inter-TM loop — the filter excluded
segments starting at residue 1 but ADGRB1's ectodomain starts at 31. Its real
inter-TM loops are 7, 20 and 9 residues. **No at-risk gene has a usable
inter-TM loop**, so classic loop insertion is not available to any of them.

Classic loop insertion — the SLC2A4 and MS4A1 route from the panel — needs a
genuine inter-transmembrane loop. Exactly one of these genes has one. For the
108 single-pass members "loop insertion" would mean inserting into the middle of
a continuous ectodomain, which only helps if the insertion point survives the
cut.

## For single-pass, the question is the juxtamembrane window

A tag survives shedding only if it sits between the cut and the membrane. That
window is set by the **most membrane-proximal** cut, not by any cut — APP alone
carries three annotated cleavages and using the β-secretase site rather than the
θ-secretase one would overstate the window by nineteen residues.

Measured against a 25-residue budget (15-residue ALFA plus minimal Gly/Ser
linkers), for the nine genes with a mapped cut position:

| gene | TM starts | closest cut | window | fits |
|---|---|---|---|---|
| EPHB2 | 544 | 535 | 9 | no |
| CDH1 | 710 | 700 | 10 | no |
| IL6R | 366 | 355 | 11 | no |
| APP | 702 | 690 | 12 | no |
| FCGR3A | 209 | 195 | 14 | no |
| CLTRN | 142 | 125 | 17 | no |
| DAG1 | 750 | 715 | 35 | **yes** |
| CHL1 | 1083 | 1039 | 44 | **yes** |
| PTK7 | 705 | 621 | 84 | **yes** |

Six of nine have no room. The ADAM and metalloproteinase cuts cluster at **9–17
residues from the membrane**, which is the biology — these enzymes work close to
the bilayer — and it is simply smaller than a tag.

Only 9 of the 122 have a mapped cut position at all, so this is a sample, not a
survey. But it is a representative one: the six that fail are the textbook shed
ectodomains.

## What to do instead

1. **DAG1, CHL1, PTK7** — place the tag in the juxtamembrane window. Assigned in
   `tag_assignments_shed.tsv`: insert after residue 732, 1060 and 662
   respectively, each centred in its window with at least ten residues clear of
   both the cut and the membrane, so the tag is neither shed with the ectodomain
   nor crowded against the bilayer.
2. **The six with no window** — do not try to engineer one. Removing a cleavage
   site risks trafficking, and these cuts are physiological. Tag them
   cytoplasmically and read surface display by an antibody against the protein
   rather than against our tag, or accept that the tag reports the uncleaved
   pool only.
3. **The 110 with no named protease** — unresolved, and worth resolving before
   any of them is assigned a tag position. N-terminomics would settle both
   whether they are cut and where.

## The honest framing

Shedding is not a construct failure. A shed tag means the surface stain reads
absent while the protein is present and behaving normally — so for these genes
the readout has to distinguish "not expressed" from "expressed and shed", and no
tag position alone can do that. A cytoplasmic-side readout alongside the surface
stain is the control that separates them.
