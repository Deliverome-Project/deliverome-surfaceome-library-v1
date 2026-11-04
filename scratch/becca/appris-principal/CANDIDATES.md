# Library candidate constructs

Every library member, isoform, ortholog and expression-rescue arm, resolved to a
verified CDS, domesticated, and dropped into the p0039 backbone. IDs are
`l1-p0001` onwards — their own namespace, so provisional designs do not consume
P-numbers, which are for constructs that get built.

## Run order

```bash
uv run --locked python resolve_candidate_cds.py     # CDS per candidate, verified
uv run --locked python domesticate_candidates.py    # clear BbsI / Esp3I / BsaI
uv run --locked python assemble_candidates.py       # into p0039, Esp3I-wrapped
```

Every input is pinned in `candidate_inputs.py` and fetched from S3, so a clean
clone reproduces the run; the only thing that varies is whether the local cache
is already warm. `candidate_inputs_provenance.json` records the exact keys used.
The expression-rescue arms come from `build_rescue_constructs.py`, which needs
`fetch_rescue_cds.py` first.

## What is deliberately not here

- **Tag cassette.** Placement is unresolved at library scale, so the insert runs
  Kozak → ORF and the cassette slot is left for later.
- **Barcodes.** Not reserved. Reserving ~7,000 against a shared, append-only
  registry before the inserts are final would attach them to sequences that are
  about to change.
- **Giant surface proteins.** LRP1, NOTCH1, FAT1, PTPRF, USH2A and similar are
  absent, and not because of anything here: they are not in
  `library_members.json`. The largest member is NRXN2 at 1,712 aa.

## How a CDS is accepted

A gene's MANE or APPRIS transcript is not always the form the library selected —
the six named GPI references are exactly that case, and their canonical
transcript encodes the non-anchored protein. So the test is translation equality
against the selected protein, tried over four routes, closest match winning:

| route | what it matches on |
|---|---|
| named transcript | the MANE or APPRIS ENST the selection recorded |
| protein match | any transcript translating to the selected protein |
| allele variant | same length, within `VARIANT_TOLERANCE` residues |
| gene symbol | any transcript of the gene at the selected length |

`allele variant` exists because GRCh38 and UniProt disagree on a residue or two
for some genes. Those are reference alleles, not different isoforms, and the
distance is recorded per candidate in `variant_residues` rather than hidden.

A CDS is refused outright if it has a premature stop, an ambiguity code, or a
length that is not a multiple of three.

## Checks that run

- Domestication must not change the protein — asserted per construct.
- No construct may carry BbsI, Esp3I or BsaI; verified by a re-scan, not by
  trusting the counter that did the removing.
- Nothing may be created at the Kozak junction or, for the rescue arms, at the
  leader-to-receptor seam. The seam detector is itself checked against a seeded
  site before it is believed.
- Each assembled vector must carry exactly two Esp3I and no other site beyond
  the backbone's own four BsaI, and a simulated digest must release a fragment
  with p0039_AXL's exact shuttle grammar.
- A key collision anywhere aborts the run rather than silently dropping
  candidates, which it previously did for 438 of them.

## Known gaps

`candidate_cds_coverage.tsv` lists every unresolved candidate with its reason.
The library ones that remain are pseudogenes with a premature stop (FCGR2C,
OR52P1, OR8J2), a gene with no clean transcript at all (KIR2DS4), and four
olfactory receptors whose transcripts are a few residues from the selected form
(OR10AC1, OR12D1, OR4Q2, OR5AL1). These need a non-Ensembl source or a decision
to drop them, not a retry.

## Giant truncations

17 constructs (`l1-p7043`–`l1-p7059`) for genes the 5 kb ORF cap dropped, built
from the per-gene literature precedents in
`analysis/sprint1/s1e3-surface-library-filters/oversized_truncation_precedents.md`.

`giant_truncation_spec.tsv` holds one row per construct: the strategy, the
residues kept, the PMID it rests on, and a confidence grade —
`literature-exact` where the paper gives residue numbers, `literature-anchored`
where it gives a strategy and UniProt supplies the boundaries, `designed` where
no minimal construct has been published and the family recipe was followed.

Cuts are resolved against UniProt features, never against a residue number
copied from a paper, and two things are asserted for every construct: it still
contains its transmembrane span, and it starts with a signal sequence. Four
fragments begin mid-protein and so have no native leader — MUC16, STAB1, STAB2 —
and take the Igκ leader sourced for the rescue arms.

The ORFs are back-translated from protein, since that is how the literature
specifies them, then domesticated like every other candidate.

Each giant is a single human construct. No isoforms or orthologs are generated
for them: an alternative form of an already-truncated design is a second remove
from the physiological protein, and the point of these is to recover the
wild-type receptor's surface behaviour, not to survey its variants.

Carrying a transmembrane span is not the same as inserting into the membrane,
so `ecto_spacer_aa` records the residues between signal cleavage and the TM and
the build refuses anything under 30. PLXNB1 originally failed it — a signal
peptide butted straight against the TM leaves nothing for signal peptidase to
cut, and the TM is then read as a signal anchor with the topology inverted — so
it keeps a juxtamembrane stub.
