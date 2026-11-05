# s4e9 — orf-design-fg-panel

## Notion page

https://www.notion.so/3e9ab86c76838139a7acc9694de6e82f

## Question

For each receptor on Brenna's 15-member FG panel, where can the SpyTag003/ALFA cassette
go, what has to change in the ORF to make it assemblable, and which barcode does it
carry — without codon-optimizing, without inventing a signal peptide, and without
colliding with the running P-series registry?

## Inputs

All immutable; nothing reads a local path except the P003 map noted under Provenance.

- `external://tedman-gpcr-surface-screen/2025-09-19/tedman2025_gpcr_screen_media-2.xlsx`
  — surface immunostaining intensity per GPCR construct, and the `Receptor Name` field
  that records Tedman's tag placement. Used to identify which panel members already
  exist as his constructs and at what expression level.
- `external://accessible-surfaceome/api-catalog/2026-08-26/catalog.json` — pinned
  surfaceome catalog snapshot, universe `cu_2026_05_12_post_silentdrop_fix`. Used for
  the floor-control screen and for symbol-to-accession resolution.
- `processed://sprint1/s1e2-fragmid-vector-design/barcodes_20nt_30k.tsv` and
  `barcode_assignments.tsv` — the barcode pool and the P-series registry of what is
  already assigned.
- Ensembl canonical CDS and reviewed UniProt entries, fetched live per gene and
  cross-checked against each other (`annotate.py`, `design_orfs.py`).

## Method

1. Resolve each symbol through `deliverome_analysis.gene_resolver`, take the Ensembl
   canonical transcript, and confirm its translated protein matches the reviewed UniProt
   entry reached from the same stable ID. A mismatch is a hard stop.
2. Derive topology — signal peptide, both termini's accessibility, GPI anchor —
   from UniProt (`annotate.py`), never from a hand-written table.
3. Domesticate each ORF against BbsI, Esp3I and BsaI with minimal edits and no codon
   optimization (`design_orfs.py`).
4. Place the tag by topology: C-terminal where the C-terminus is extracellular,
   otherwise N-terminal after the signal-peptide cleavage site, otherwise an
   extracellular-loop insertion (`build_constructs.py`, `build_loop_inserts.py`).
5. Add inward-pointing Esp3I shuttle sites, with overhangs chosen on NEB ligation
   fidelity (`shuttle.py`).
6. Reserve barcodes from the registry with boundary-enzyme and splice-donor screens
   applied in the real junction context (`reserve_barcodes.py`).

Every step asserts its own correctness: protein unchanged, no Type IIS site created at
any junction, digest regenerates the insert, and the TFRC reference arm reproduces
P003's insert region byte-for-byte.

## Status

**In progress — nothing ordered.** The constructs are designed and verified but several
decisions are open; see Open questions on the Notion page. The tag-placement choice
(N- versus C-terminal) is untested and is the main risk.

## Scripts

| | |
|---|---|
| `design_orfs.py` | symbol → `gene_resolver` → Ensembl canonical CDS → UniProt sequence check → Golden Gate domestication. No codon optimization. |
| `build_constructs.py` | wraps a domesticated ORF in the p003 cassette format (Kozak, SpyTag003/ALFA cassette, stop, barcode amplicon). `--placement nterm\|cterm`. |
| `barcodes.tsv` | the barcodes, drawn from the s1e2 30k production pool by index. |
| `tag_arms.py` | the placement/order arms: SLC20A2 (N-term / C-term / C-term swapped) and TFRC (p003 reference / C-term swapped). |
| `tag_sites_provenance.json` | how the tag sites were generated, since the upstream pipeline's own store is not being written. |
| `puro_only_backbone.py` | builds and verifies the puro-only intermediate (IRES2→PuroR). No P-number: never ordered, released 2026-10-18. |
| `reserve_barcodes.py` | reserves barcodes from the running P-series registry on S3, with boundary-enzyme and splice-donor screens. |
| `shuttle.py` | wraps a finished insert in two inward-pointing Esp3I sites for one-pot transfer to another backbone (s1e2 §4.11). |

```bash
uv run --locked python design_orfs.py
uv run --locked python build_constructs.py --placement nterm
uv run --locked python shuttle.py                                 # terminal constructs
uv run --locked python shuttle.py --source constructs_loop.fasta  # loop constructs
```

## Esp3I shuttle sites

Added over p003 so the barcoded ORF can move to another backbone in one Golden Gate
reaction — Tedman's reason, and s1e2 §4.11's: if the barcode travels **with** the ORF,
the expensive long-read variant-barcode map is built once and amortised across every
future backbone.

```
…EFS… CGTCTC A [GACC] Kozak-ORF-cassette-stop-barcode-revprimer AG [TCAG] A GAGACG GGCGCGCC …IRES2…
      <-- cuts right                                                     cuts left -->
```

Both sites point **inward**, so both recognition sequences stay behind in the donor and
the transferred fragment carries no Esp3I site. Cost is **+16 nt** per insert.

### Overhangs are chosen on NEB ligation-fidelity data

NEBridge Ligase Fidelity Viewer, dataset **"Esp3I, 1x T4 DNA Ligase buffer, 37-16
cycling"**, read 2026-09-22 (<https://ligasefidelity.neb.com/>).

| Pair | Fidelity | Correct-pair counts | Mismatch counts |
|---|---|---|---|
| GACC + CGTC | **95%** | 217 / 232 | 8 (GACC:CGTC), 4 (GGTC:GACG) |
| **GACC + TCAG** | **100%** | **217 / 355** | **0** ← chosen |
| GACC + TCCA | 100% | 217 / 317 | 0 |
| GACC + TCAC | 100% | 217 / 239 | 0 |
| GACC + CACG | 100% | 217 / 177 | 0 |

The obvious overhang, `CGTC` (the natural last four bases of the reverse primer site),
**loses 5%**: `GACC` differs from `GACG` — CGTC's ligation partner — at a single
position, so the two ends cross-ligate. `TCAG` is clean and is also the strongest
ligator in the set at 355 counts.

`TCAG` costs a **2-nt pad (`AG`)** after the reverse primer site. It sits inside the
transferred fragment but in the non-coding gap before the IRES, so no coding sequence,
no primer site and no barcode is touched. `OH5 = GACC` is kept: it is the first four
bases of the Kozak block, so it adds nothing, and it pairs cleanly with TCAG.

Overhangs are distinct, non-palindromic, and not each other's complement, so the
fragment cannot self-ligate or insert backwards.

### Junctions, checked against the real p003 map

| Junction | Sequence | Result |
|---|---|---|
| EFS → Kozak | `…GGTTTGCCGCCAGAACACAG` `CGTCTCA` `GACC…` | 1 Esp3I (+), 0 BbsI, 0 BsaI |
| barcode → IRES | `…CACGTC` `AGAGAGACG` `GGCGCGCC` `…IRES2` | 1 Esp3I (−), 0 BbsI, 0 BsaI |

The AscI site (`GGCGCGCC`) before IRES2 is preserved. Because every ORF is domesticated
against Esp3I, these two are the **only** Esp3I sites in the construct, and p003, p004
and pRDA_722 are all Esp3I = 0 (s1e2 §4.9/§4.10).

**One-way transfer.** The recognition sequences stay with the donor, so the product of a
transfer carries no Esp3I sites — the destination vector has to supply its own if the
cassette is to move again. That is the normal Golden Gate shuttle arrangement and what
Tedman describes; say so if repeat-shuttleable is wanted instead.

⚠️ **This spends the Esp3I reservation.** s1e2 §1 reserves Esp3I/BsmBI for downstream
guide cloning; §4.11 records that reservation as contested and as the only objection.
Using it here means guide cloning needs another enzyme (AarI/PaqCI is the best
unreserved option, at the cost of one IRES2 edit).

## Signal peptides — none of the nine needs one added

Seven carry their own (`signal_peptide` column in `orf_design.tsv`): HEG1 1–29,
ERBB2 1–22, AXL 1–25, SLC9A3 1–25, VSIG10 1–30, NECTIN4 1–31, FOLR1 1–24. Keep them.

### InterPro check on the eight junctions (2026-10-18)

Reviewed by the InterPro site-placement session, which recovered all eight junctions
independently from the built nucleotides in `processed://constructs/designed/` and agreed
with every one. "Inside a domain" is restricted to InterPro entry type `domain`/`repeat`;
`family` entries routinely span a whole protein and would trivially contain any junction.

| gene | junction | verdict |
|---|---:|---|
| HEG1 | 29 | clear — first domain at 982 |
| PTPRC | 41 | clear — first domain at 391 |
| SLC9A3 | 26 | clear — first domain at 62 |
| FOLR1 | 23 | clear — first domain at 36 |
| AXL | 33 | profile over-extension only: IPR007110 starts at 27, five residues *inside* the signal peptide (cleavage 32). First genuine start is 34. |
| NECTIN4 | 32 | same pattern: IPR007110 at 30 against cleavage 31. First genuine start 34. |
| ERBB2 | 23 | flush with the Receptor L-domain superfamily's first residue, not inside a `domain` entry |
| **VSIG10** | **31** | **genuinely inside IPR013106 / PF07686 Ig V-set (31–105)** — the domain begins at the first mature residue |

**VSIG10 still does not move, because there is nowhere to move it to.** Tedman's
admissible interval (`accessible_surfaceome.tag_sites.sp_interval`) is 30–32. Excluding the
cleavage site itself — which none of Tedman's 108 signal-peptide-bearing constructs uses —
only 31 and 32 remain, and both are inside the V-set domain. AXL (interval 32–36) and
NECTIN4 (31–33) are in the same position.

The general result is the one worth carrying forward: **for compact single-pass receptors
the structured ectodomain begins essentially at the mature N-terminus**, so there is no
disordered linker after the signal peptide to tag into. Our disorder-consensus onset and
InterPro's first-domain start agree to 0–2 residues on exactly the four genes with no room
(VSIG10 31/31, ERBB2 23/24, AXL 34/35, NECTIN4 34/32) and diverge widely where there is
plenty (PTPRC 391/226, SLC9A3 62/42). Two methods corroborating on the decision that
matters.

**Caveat, carried verbatim from the review:** Tedman's own junctions have *not* been
measured against InterPro, so there is no evidence that sitting at an Ig V-set start is
safe — only that VSIG10 has no better option and that this geometry is common among the
receptors his screen validated. That measurement is a planned task, not a result.

On the Igκ premise, InterPro cannot answer: it carries no signal-peptide annotation at all
(signal peptides reach InterProScan output as a SignalP sequence feature, not an InterPro
entry). The review checked the predictors directly instead — SignalP 6 returns `OTHER` with
p_sp = 0.0000 for all five GPCRs, and DeepTMHMM v1 starts all five with `O` and no `S` in
the first 45 residues. With UniProt and Tedman's `_SP_HA` set, that is four independent
sources agreeing there is nothing for the leader to replace.

SLC2A4 and MS4A1 have **cytoplasmic N-termini**; an added SP would invert topology, the
same argument s1e2 §4.1 makes for TFRC. Neither gets one.

This matches Tedman, who adds no exogenous signal peptide anywhere and criticizes
PRESTO-Tango for doing so; for his own SP-bearing receptors he keeps the native SP and
puts the tag *after* its cleavage site.

### FOLR1's cassette sits inside its predicted fold — and p0028 is the only read on that

**p0028 is built at residue 27, three past the end of its admissible interval (22–24), and
two of three disorder predictors call that residue structured.** Measured against
`disorder_public` on 2026-10-19:

| predictor | score at residue 27 | calls it | puts the fold's start at |
|---|---|---|---|
| NetSurfP-3.0 | 0.859 | disordered | 31 |
| metapredict | 0.369 | **ordered** | 23 |
| AlphaFold-disorder | 0.255 | **ordered** | 23 |

The median of three is ordered, and the two that disagree with NetSurfP both put the fold
beginning immediately after the signal peptide. InterPro does not corroborate it — its first
FOLR1 domain is at 36 — so the two methods genuinely disagree about this protein.

**Its signal-peptide re-prediction came through untouched: 0.982 → 0.981.** That is the
point worth carrying rather than the FOLR1 result itself. A clean SignalP verdict says
nothing about whether the cassette landed in structure; the SP check and the fold check are
different checks, and across this batch one had been standing in for the other.

**So p0028's display result is the only read this batch has on whether a fold tolerates the
cassette inserted into it.** Nothing else here puts a cassette inside predicted structure.
That question is live for a whole class upstream — the members whose structured ectodomain
starts at the mature N-terminus, where the interval rule yields no strict candidate and
every available position is inside the fold (the InterPro site-placement session counts 419
such library members, and treats this arm as the test case). When FOLR1's expression and
display data come back, read them as that test and not only as one more panel member.

**Nobody has published a FOLR1 tag insertion.** All three sites `tag_site_public` records
for it are inference, not literature: `terminal-n` after 22, `terminal-c` after 257, and a
disorder site after 246. So there is no precedent to copy and our junction is not departing
from one. The corpus's own N-terminal suggestion is after **22**, the cleavage site itself,
with a 15-aa ALFA tag rather than our 51-aa cassette.

**And its `terminal-c` site cannot work — which is a library-wide defect, not a FOLR1 one.**
FOLR1's recorded C-terminal site sits after residue 257, inside the 235–257 GPI propeptide.

**The mechanism, corrected 2026-10-19.** This section first said the tag is *cleaved off
with the propeptide*. That is wrong. The GPI signal has to **be** the C-terminus to be
recognised by the transamidase, so appending anything past it **prevents the anchor being
added at all** — the protein is not anchored rather than anchored-and-untagged. Same
outcome for the construct, different reason, and the corrected version is the one to
propagate.

Counted against UniProt's own list (KW-0336, reviewed human, 144 proteins): **118 carry a
`terminal_c` site — 125 rows, 11.7% of all 1,068.** 116 inference-only, 9 literature-backed
(CA4, CD59, CNTN2, DPEP1, IGSF21, LYPD2, PLAUR, PRSS55, SPRN). The tag-sites session counts
115/10 with ISLR included, because it sources GPI membership from **NetGPI** rather than
UniProt; ISLR is NetGPI-called and not in KW-0336. Both splits are right within their own
source and neither of us labelled which we were using — the same denominator failure this
work keeps producing, mutual this time. Their own finding is that **UniProt should be read
first** and NetGPI used as corroboration, since NetGPI misses ITLN1, TECTA and TECTB.

**The nine literature rows are not counter-evidence.** All were pulled and read: every one
is an anchor-replacing or anchor-blocking fusion, not a tag at the mature terminus. IGSF21
states the rule outright — HA "inserted in frame 4 aa residues *upstream* of the putative
C-terminal GPI anchor signal sequence" (PMID 28864826); CA4 "replac[es] their C-terminal
GPI anchors **with** a hemagglutinin tag" (PMID 41175373). So the published constructs
avoid exactly the position the lane proposes.

**It is a missing input, not a missing check.** The rule is "topology string ends
extracellular → propose a C-terminal tag", and a GPI protein's topology string does end
extracellular. No column in `tag_site_public` can express the constraint: no GPI,
lipidation, anchor, propeptide or omega field among its 37 columns. Nothing was skipped;
there was nothing to check against.

**And it does not reach this batch or the L1 library.** Of the 120 genes with a refused
`terminal_c`, exactly **two are L1 members — BST2 (`l1-p0188`) and RGMB (`l1-p1731`)** — and
both are already flagged `gpi_rerouted` in `l1_internal_members.tsv`, rerouted because the
omega site removes the C-terminus. That column having only two entries was correct rather
than thin.

Nothing in this batch depends on it: FOLR1 is this panel's only GPI member and it is
N-terminal by design, for exactly this reason.

No change is proposed to p0028. FOLR1 is GPI-anchored, so the N-terminus is the only end it
has; the admissible positions 23–24 are *closer* to the cleavage site rather than further
into open sequence; and InterPro's 36 would permit 27. There is no move that is clearly
better, only different — so this is recorded as a known risk carried deliberately.

## Tag placement — the C-terminal cassette does not transfer

**TFRC is the only member of the whole panel with an extracellular C-terminus.** It is a
type II protein; the p003 cassette sits on its ectodomain. For these nine:

| | C-terminus | Terminal cassette |
|---|---|---|
| HEG1, ERBB2, AXL, SLC9A3, VSIG10, NECTIN4 | cytoplasmic | **N-terminal**, after the SP cleavage site |
| FOLR1 | GPI signal, residues 235–257 removed at Ser234 | **N-terminal** — a C-terminal tag is physically cleaved off |
| SLC2A4, MS4A1 | cytoplasmic (both termini in) | **impossible** — needs an extracellular-loop insertion |

`--placement nterm` is therefore the default and produces seven constructs.
`--placement cterm` reproduces p003 verbatim for all nine and is kept only for
comparison — for eight of them it puts SpyTag003 and ALFA inside the cell.

SLC2A4 and MS4A1 are the forced-insertion cohort (~19.8% of the surfaceome per
`analysis/sprint3/s3e13-flow-screen-scenarios/tag_sites.py`); the accessible-surfaceome
`tag_sites` pipeline already scores loop sites for ~800 genes and is the right source.

**Open:** the cassette is 51 aa against Tedman's 9-aa HA. Dropping 51 residues
immediately after a signal-peptide cleavage site is a much larger perturbation than he
validated. Tedman refines the exact insertion point per receptor with SignalP 6.0 +
NetSurfP-2.0, testing several positions between the cleavage site and the start of the
structured N-terminal domain. We place the cassette immediately after the cleavage site;
that refinement is **not yet done**, and splitting the two tags (e.g. ALFA N-terminal,
SpyTag elsewhere) has not been evaluated.

## Outputs

| File | Contents |
|---|---|
| `orf_design.tsv` | per gene: stable IDs, UniProt, topology, sites removed, edit count and positions, GC |
| `orf_native.fasta` / `orf_domesticated.fasta` | CDS before and after domestication |
| `constructs_nterm.tsv` / `.fasta` | the seven assembled inserts (default) |
| `constructs_cterm.tsv` / `.fasta` | p003-verbatim comparison, nine inserts |
| `skipped_nterm.json` | SLC2A4, MS4A1 and why |

## Two corrections to the s1e2 README

Both found by reading the GenBank maps rather than the doc.

1. **§3 row 11 calls `tagcactgagtcagatctgag` "likely a cloning scar, not the barcode."
   It *is* the barcode.** P010 (`…_NNNbc_…`) carries twenty `N`s at exactly that window;
   p003, P008, P009, P011, P012 and P013 all carry `TAGCACTGAGTCAGATCTGA` there.
2. **§3 row 12's `CCTTTGGGTAAGCACACGTC` is not the barcode** — it is part of the fixed
   reverse primer site, and all seven maps read `CCTTTGGGC AAGCACACGTC` (C, not T, at
   position 9).

Also worth noting: §4.5 warns that both stops are `TGA`. The p003 ORF stop is `TAA`
followed by the `GG` spacer — `TAAG`, already compliant with the `TAA[A,G]` policy.

## Library creation

The fifteen panel ORFs above are designed one at a time. The 2,369-member library is not:
adding a member has to cost a database lookup and a decision-tree walk, not a design
session. The tree, its inputs and the evidence behind each branch are in
[METHODS_library_design.md](../../../scratch/becca/tag-site-architecture/METHODS_library_design.md),
also published at `processed://sprint4/s4e9-orf-design-fg-panel/library-design/METHODS_library_design.md`.
This section records what the tree produces and where each output lives.

**DeepTMHMM v1 is the topology of record** — v1 reaches 99.5% signal-peptide recall against
v2's 92.9%, so v2 is used for membrane compartment and reentrant helices, not for placement.
Topology is read from public D1 `topology_public`, `topology_version='topo_2026_05_16'`.

### Placement is earmarked per member, separately from whether a site exists

Every member is assigned exactly one **rung** from topology alone (`derive_tag_paths.py`).
A rung is a commitment about *where* the cassette goes; it says nothing about whether a
concrete junction has been computed yet. Those are different facts and the library needs
both.

| Rung | Members | Rests on | Members with a site in `tag_site_public` |
|---|---|---|---|
| N-terminal, after signal-peptide cleavage | 1,129 | predicted cleavage site | 1,025 |
| N-terminal at residue 1, multipass | 525 | Tedman GPCR screen | 345 |
| N-terminal at residue 1, single-pass (type III) | 32 | topogenesis, untested | 13 |
| C-terminal | 166 | topology only | 104 |
| Internal loop or snorkel | 515 | AlphaFold model per gene | 431 insertion + 26 snorkel-only |

Of 2,367 members with canonical topology, 1,918 have at least one extracellular
deterministic site and 423 have none.

The internal rung is the hard one and is **canonical genes only** — both termini
cytoplasmic, so no terminal tag is displayed. Its TM counts are necessarily even, which
excludes every single-pass class and all 7-TM receptors. Space is rarely the limit: the
median longest contiguous extracellular run is 51 residues and only 5 of 515 fall under 8.
What blocks the thin members is `MIN_LOOP_LEN = 8` applied to the contiguous DSSP coil/turn
run *inside* those loops, which is 3–7 residues for them — KCNMB2 has a 129-residue
extracellular run and a best coil run of 7, so it is gate-blocked, not space-blocked.

A snorkel site is emitted for every member with no accessible terminus, as a documented
fallback ([Brown et al., 2013](https://pubmed.ncbi.nlm.nih.gov/24023844/)) at low
confidence. All 431 internal-rung members carrying a loop candidate also carry a snorkel
site, so its presence is not evidence that the loop lane failed; it is decisive only for
the 26 where no loop cleared the gates.

### GPI members are resolved per form, not per gene

GPI membership is **read** from the curated surfaceome architecture attribute, with UniProt
KW-0336 and a cited supplement unioned in — not predicted. NetGPI 1.1 is then used only to
locate the anchor *within* a gene, as a veto on confident negatives and never as a ranking
signal, because the attribute is a claim about the gene rather than about the form the
library picked. Calibrated against 148 curated-GPI and 100 single-pass references it reads
95.0% sensitivity and 99.0% specificity, MCC 0.875 against the published 0.895.

Six references move onto the form that carries the anchor: ACHE (P22303 → P22303-2, the H
subunit over the ColQ/PRiMA-assembling T subunit), BCAN, CD55, CD58, NCAM1 (NCAM-140 →
NCAM-120) and TFPI (TFPI-α → TFPI-β). Per-form verdicts across 147 genes are 132
`reference_is_anchored`, 6 `anchor_on_another_form`, 5 `no_anchored_form_found` and 4
`unresolved`. A gene the predictor alone rejects is reviewed, never dropped: CRIPTO, CFC1
and GML carry KW-0336 and are predictor blind spots.

A C-terminal cassette on a GPI member is forbidden outright: it sits past the ω-site, so
the transamidase either removes it with the GPI signal — anchored but untagged, reading as
a non-expressor — or the 51-residue addition breaks the signal's length and hydrophobicity
and the protein is secreted. BST2 is why this is the first branch of the tree rather than a
sub-case: it is GPI-anchored with no cleaved signal peptide, so on topology alone it reaches
the C-terminal rung.

### Isoform and ortholog panels

Library constructs are named `l1-pNNNN`. Two panels extend the canonical set, both drawn
class-balanced across type I, type II, type III, GPI and multipass rather than uniformly at
random, so a rare architecture is not lost to sampling.

The ortholog panel is 200 and the isoform panel 150 — 192 and 137 genes respectively.

| Panel | Members | Composition | Key |
|---|---|---|---|
| Isoforms | 150 | 100 observed + 50 prediction-only; type I 61, multipass 51, type II 17, GPI 16, type III 5 | `library-design/isoform_panel.tsv` |
| Orthologs | 200 | mouse 100 + cyno 100, split near/diverged at 50 per band, seed 20261011, minimum floors on type III, GPI and type II; type I 66, multipass 72, GPI 28, type II 26, type III 8 | `library-design/ortholog_panel.tsv` |

**The panels were drawn from April's gene set, not the library's.** Both universes were
built from accessible-surfaceome files frozen at commit f9c73691 (2026-04-27), whose human
query set is deliverome-internal's expression-filtered `surfaceome_expressed` list (2,359
proteins, 2,355 gene symbols). The library gene list came later, from the S1E3 filter over
surfaceome release 1.3.0, and only 1,775 of its 2,369 genes are in the April set. Later
Compara pulls (2026-06-01, 2026-09-20/21, and the 2026-10-12 consolidation onto Compara
115) went into the surfaceome D1 tables, not into the f9c73691 files these universes read,
so the other 594 library genes were never queried. `check_universe_coverage.py` measures it
from pinned inputs (Ensembl Compara 115 dumps, the f9c73691 query set and ortholog lists,
Ensembl 115, UniProt 2026_03) and writes `universe_coverage.tsv`, one row per library gene:

| | Mouse | Cyno |
|---|---:|---:|
| Library genes with a one-to-one high-confidence ortholog in Compara 115 | 1,870 | 1,848 |
| … in the ortholog list the panel drew from | 1,452 | 1,432 |
| … never queried (all outside the April set) | 418 | 416 |
| Recorded pairs that are not one-to-one high-confidence in Compara 115 | 0 of 1,972 | 0 of 1,927 |

The isoform universe has the same cause for 92 genes: 451 library genes with a
non-canonical UniProt isoform have no row in it, 92 of them outside the April set; the
other 359 were offered and excluded by the universe's own rules, which this does not
re-examine. The panels are a valid seeded draw from that smaller pool and are kept as they
are; the gap is recorded rather than repaired.

The 515 internal-rung members carry 607 surface isoforms across 287 genes and 818 surface
orthologs across 474 genes, so a build covering every form of that rung alone would be 1,940
constructs rather than 515.

**Walton swap-risk screening does not transfer to these panels.** The ≤100 bp shared-tract
bar is a property of the pool, and an isoform shares its parent's sequence by construction:
the median canonical-versus-isoform overlap is 357 aa (1,071 nt) across 2,150 pairs, with
zero pairs under the bar and the p10 still at 306 nt. Applying the tract bar to isoforms and
orthologs would remove essentially all of them, which defeats the reason they are included.
The cassette diversification in `diversify_cassettes.py` addresses the constant every member
carries, which is the substrate the bar is actually protecting against.

### Materials and provenance

| Element | Location | Status |
|---|---|---|
| Decision tree and evidence | `library-design/METHODS_library_design.md` | on S3 and in `scratch/becca/tag-site-architecture/` |
| Library membership (2,369) | `library-design/library_members.json` | on S3 |
| Exclusions, with reasons | `library-design/library_exclusions.tsv` | on S3 |
| GPI form calls and verdicts | `library-design/netgpi_*.tsv`, `gpi_selected_forms.tsv` | on S3, 8 NetGPI keys |
| Isoform panel | `library-design/isoform_panel.tsv` | on S3 |
| Ortholog panel and draw summary | `library-design/ortholog_panel.tsv`, `ortholog_panel_summary.json` | on S3 |
| Candidate constructs and CDS | `library-candidates/candidate_{manifest.tsv,constructs.json,cds.json}` | on S3, 31 keys |
| Tag-site constraints per member | `library-candidates/library_tag_constraints.tsv` | on S3 |
| Oversized and giant-truncation triage | `library-candidates/{OVERSIZED_SURVEY,giant_*}` | on S3 |
| Per-member rung assignment | `library-design/tag_path_assignments.tsv`, `tag_path_breakdown.tsv` | on S3, provenance-stamped; committed via `git add -f` |
| Internal-rung worklist (515) | `library-design/internal_tag_targets_515.tsv` | on S3, provenance-stamped |
| Per-member tag paths (all 2,707 active L1 with a designed CDS) | `l1_tag_paths.tsv`, `TAG_PATHS.md` | committed (force-added past the `*.tsv` ignore) and on S3; **the live answer** — 567 internal-site members |
| GPI set | `gpi_set.tsv`, `build_gpi_set.py` | committed and on S3; 157 genes from NetGPI over reference *and* alternative forms |

S3 paths are relative to `processed://sprint4/s4e9-orf-design-fg-panel/`; repo paths are from
the repository root. The repo-wide `*.tsv` ignore means a table reaches version control only
by `git add -f` — so a derived table is durable only if someone did that *and* published it.

**Provenance stamping is uneven and it is worth checking before trusting a key.**
`processed.upload()` records the commit, a dirty flag and the lock hashes as S3 metadata, and
the rung tables published here read back `git-commit 31784e1`, `git-dirty false`. Most older
keys under this prefix — `candidate_manifest.tsv`, `ortholog_panel.tsv`, even
`constructs/designed/p0022.gb` — return "no provenance recorded" despite being written after
stamping landed, so they cannot be traced to a commit. Read a key's stamp with
`uv run python -m deliverome_analysis.provenance processed://<key>` rather than assuming.

### Structural coverage: every member needs a model before it can be scored

The tag-site lanes and EpicTope both read pLDDT, relative solvent accessibility and
secondary structure off an AlphaFold model. A member without one is not mis-scored, it is
**unscoreable** — no lane can propose a site for it at all. Coverage is therefore a
prerequisite for the library, not a refinement of it.

`l1_structure_gap.py` measures it per tranche. Of L1's 2,725 candidates 2,654 already had
an AlphaFold DB model; 71 did not, one of which (TPO) was a gene-symbol lookup miss rather
than a real gap, and the remaining 70 are predicted here.

**Having a model is the weaker question, and answering only it overstates coverage.**
`check_proteoform_match.py` translates each candidate's designed CDS from
`library-candidates/candidate_cds.json` and compares it residue-for-residue against the
AlphaFold entry's own sequence. Over the 2,667 L1 members that have a model:

| Verdict | Members | Meaning |
|---|---:|---|
| identical | 2,417 (90.6%) | the model is of exactly this protein |
| contained | 26 (1.0%) | our protein is a subsequence of the modelled one; 25 of the 26 need a 6-134 residue offset the lanes do not apply |
| **different** | **199 (7.5%)** | the model is a different protein; its numbering cannot be carried onto the construct |
| not compared | 25 | no designed CDS, or no sequence returned |

The mismatch is concentrated in `library_reference` (191 of the 199), which is the one set
resolved by **gene symbol to a canonical accession** rather than by the accession its
sequence came from; the library does not always build the canonical proteoform. Neurexins
and MUC1 are the extremes — NRXN2 is modelled at 666 aa against an ORF of 1,712, MUC1 at
1,255 against 273, the latter because the canonical VNTR is unsynthesisable and a short
form was chosen deliberately.

**A length comparison is not sufficient and was tried first.** 23 of the 199 mismatches are
the *same length* as the model with a different sequence, so a length test reports them as
fine. The check is residue-exact for that reason.

This is worse than a missing model rather than better: a missing model is visible and gets
folded, while a wrong-proteoform model yields tag-site coordinates numbered against a
sequence the construct does not contain, and nothing downstream flags it. So **L1 is not
structurally complete**: on top of the 70 folded here, roughly 225 members either need
their own ORF folded or need a verified offset.

| Tranche | Candidates | Have a model | Gap |
|---|---:|---:|---:|
| L1 reference | 2,360 | 2,359 | 1, and that one resolves |
| L1 expression rescue | 25 | 25 | none |
| L1 isoform panel | 144 | 106 | 38 |
| L1 ortholog panel | 196 | 164 | 32 |
| **L1 total** | **2,725** | **2,654** | **70 folded, plus ~225 whose model is the wrong protein** |
| L2 | 4,341 | 3,186 | **1,155**, not folded |

**L2 is measured but not folded: 1,155 of its 4,341 candidates have no model** — 689
isoforms, 463 orthologs and 3 giant truncations. (An earlier extrapolation put it near
1,200, which was close, but the number above is counted rather than inferred.) Re-measure
with:

```bash
uv run --locked python l1_structure_gap.py --library l2 --refresh --out l2_structure_gap.tsv
```

Then `build_fold_inputs.py` fetches the sequences and
`scripts/fold_missing_structures_modal.py` folds them. The L1 set runs at roughly 4.3 GPU
minutes per protein, so 70 costs about 3 GPU-hours; L2 at ~1,200 would be nearer 85, which
is a spending decision rather than a rerun.

**These are ColabFold predictions, not AlphaFold DB entries.** The model weights are the
same AlphaFold2, but the multiple sequence alignment comes from MMseqs2 against
UniRef and environmental sets rather than DeepMind's JackHMMER and BFD pipeline, so the
two are comparable without being interchangeable. They publish under
`processed://sprint4/s4e9-orf-design-fg-panel/predicted-structures/colabfold/<date>/` and
the summary records the ColabFold version, so the distinction survives into the record.

Coordinates do not belong in D1, which caps a row at 2 MB. The viewer's convention, set by
`ingest_schweke_pdbs_to_r2.py` in accessible-surfaceome, is that coordinates live in the
public R2 bucket `surfaceome-structures` and the metadata index lives in D1. Publishing
these to R2 with a D1 index row is outstanding; S3 alone serves the analysis pipelines but
not the viewer.

### Cassette encoding: only the linkers vary

The tag cassette is GGGGS | SpyTag003 | GGGGSGGGGS | P-ALFA v1.5-P | GGGGS (51 aa). Every
member gets its own synonymous encoding so that members do not share a long identical
stretch before the barcode, which is where a reverse-transcription template switch moves
one member's barcode onto another's ORF unseen. The shared stretch is measured up to the
barcode and gated at 100 bp (5.0% predicted swapping, Walton's fit in
`seqtools.barcode_swap_rate`), 80 bp preferred.

**SpyTag003 and P-ALFA-P keep p003's codons in every new encoding; only the four GGGGS
blocks are recoded.** p003's tag codons are the ones validated in cells, and recoding the
tags buys no swap protection: what sets the stretch before the barcode in a C-terminal
member is its 3'-most GGGGS block, followed by the constant 24-nt tail (TAA, CC pad,
19-nt handle). `diversify_cassettes.py` and `l1_swap_risk.py` pass the tag codons to
`seqtools.diversify_library_linkers(fixed_codons=...)`.

Each C-terminal member must also end on a 3' GGGGS block of its own, because all of them
share the tail that follows. How many distinct ends exist depends on how many terminal
nucleotides must differ, after excluding runs of 5 or more identical bases:

| Ends required to differ in | Distinct encodings | Enough for the 218 C-terminal L1 members? |
|---|---:|---|
| last 9 nt (GGS) — the original rule | 60 | no: the encoder stopped at member 60 |
| last 12 nt (GGGS) | 180 | no |
| whole 15-nt GGGGS block — current rule | 540 | yes |

The rule is therefore that the whole terminal block is unique per member
(`LINKER_ONLY_JUNCTION_IDENTITY = 14`). Two C-terminal members can then share at most
24 + 14 = 38 nt before the barcode.

Measured over all 2,707 active L1 members (`l1_swap_risk.py`; per member in
`l1_swap_risk.tsv`, summary in `l1_swap_risk_summary.json`). The worst-partner tract is
the longest stretch a member shares, up to its barcode, with any other member:

| Scheme | All members: median worst-partner tract | C-terminal (218): median / max | Members over 100 bp |
|---|---:|---:|---:|
| one shared cassette (p003's) | 30 nt | 181 / 2,529 nt | 482 |
| every codon recoded, original 9-nt end rule | encoder stops after 60 of 218 | — | — |
| every codon recoded, 15-nt end rule (control) | 30 nt | 36 / 36 nt | 264 |
| **linkers only, 15-nt end rule (current)** | **30 nt** | **36 / 36 nt** | **264** |

The library median stays at 30 nt because 2,489 of the 2,707 members are tagged upstream:
their cassette is not between the ORF and the barcode, so they share only the 24-nt tail
plus whatever C-terminal codons two ORFs have in common. All 264 remaining members over the
bar are explained by the ORFs themselves, not by the cassette. 211 are same-gene
isoform/ortholog pairs that share a C-terminus, which the 2026-10-13 decision accepts. 53
are different genes with a shared 3' coding sequence: the clustered protocadherins (37,
which share constant exons; up to 486 nt) and close paralogs such as OR2T5/OR2T29
(898 nt), OR8U1/OR8U9, OR4F5/OR4F17 and NPY4R/NPY4R2. The cassette cannot shorten
those, because the shared sequence lies upstream of it.

One C-terminal member needs a junction edit when it is assembled. SLC41A2 (l1-p1957) ends
`…GGA GAC`, and any Gly codon after it completes an Esp3I site. The synonymous
last-codon change GAC→GAT (Asp) clears the site, and `l1_swap_risk.py` applies it.
The 218 linkers-only encodings are issued in `l1_cassette_encodings.tsv`. Once issued, an
encoding is pinned: a rerun draws only members that have none, exactly as
`diversify_cassettes.py` does for the panel. Without pinning, adding two C-terminal members
redrew 156 of the 216 encodings already issued, because members are served in risk order.

**The 25-construct order (p0022–p0049, and its p0051–p0075 copies in the P-0050 backbone)
is unchanged and keeps the earlier scheme.** 24 of its 25 constructs carry recoded tags,
with 23–24 distinct encodings of each cassette element. Only p0037 (TFRC) carries p003's
codons. The encoded peptide is identical in all 25, so protein-level readouts (SpyCatcher,
anti-ALFA) are unaffected. Any assay that targets the tag DNA, such as a primer or probe
on the SpyTag003 or ALFA nucleotides, will not read these constructs uniformly, and they
will not share tag DNA with L1. `cassette_encodings.tsv` pins all 23 issued encodings as
built and labels them `whole_cassette` in a new `scheme` column.

### Members whose protein duplicates another's are retired

A member's protein is its designed CDS translated. Comparing those proteins across the
library found 53 groups (108 members) in which two or more members encode the identical
protein, and each pair had a different gene symbol, accession or transcript.
`dedupe_library_members.py` keeps one member per protein, by precedence: library reference
before expression-rescue, isoform, ortholog or giant truncation; L1 before L2; then fewer
residues off the source record. It marks the others `retired_reason` in
`candidate_manifest_tranches.tsv`. A retired member keeps its `candidate_id`, since ids are
never renumbered, and `library_design` readers skip it. The full list with rationales is
in `library_duplicate_members.tsv`.

| Retired | Kept | Category | Root cause |
|---|---|---|---|
| l1-p2387 ACHE isoform P22303-2 | l1-p0017 ACHE | isoform is the reference | The GPI rule moved ACHE's reference onto P22303-2, but `isoform_universe.py` still took the MANE/canonical form as reference, so it offered P22303-2 as an alternate |
| l1-p2411 CD177 isoform Q8N6Q3 | l1-p0263 CD177 | allele variant of the reference | Q8N6Q3's 3 differences (V184G, D204N, M205R) are UniProt sequence conflicts from one cDNA. The resolver's 3-residue allele tolerance then built it from the reference's own transcript |
| l1-p2585 MMP16 cyno ortholog G7PC63 | l1-p1225 MMP16 | ortholog identical to human | The near band of the ortholog draw had no 100%-identity exclusion |
| l1-p1350 OPN1MW2, l1-p1351 OPN1MW3 | l1-p1349 OPN1MW | identical paralog | Membership is per gene symbol, and tandem copies in the opsin array encode one protein |
| l1-p1406 OR2A42 | l1-p1404 OR2A1 | identical paralog | as above (segmental duplication) |
| l1-p2318 UPK3BL2 | l1-p2317 UPK3BL1 | identical paralog | as above |

| Tranche | Active before | Retired as duplicates | Panel redraw | Active now |
|---|---:|---:|---:|---:|
| L1 (members with a designed CDS) | 2,700 | 7 | +13 | 2,706 |
| L1 (all rows) | 2,725 | 7 | +13 | 2,731 |
| L2 | 4,340 | 48 | −13 moved to L1 | 4,279 |

L2's 48 follow the same causes: 39 orthologs identical to human, 5 isoforms that are the
reference (BCAN, CD55 and TFPI are three of the other GPI-override genes; SUSD6 and PILRA),
2 members of another gene (an isoform of KIR2DS5 equal to the KIR2DL5A reference, and an
SLCO1B3 isoform equal to the SLCO1B3-SLCO1B7 read-through reference) and 2 duplicate
isoform keys (HLA-A, NRXN1). No retired member was C-terminal, so the 216 C-terminal
members are unchanged. OPN1MW2, OPN1MW3, OR2A42 and UPK3BL2 no longer have a member of
their own; their protein is still built under the kept symbol, which carries the others in
a new `also_encodes` manifest column (OPN1MW: OPN1MW2;OPN1MW3, and likewise OR2A1,
UPK3BL1, KIR2DL5A, SLCO1B3-SLCO1B7), so a readout for the kept member is read for every
symbol that encodes it.

The causes are fixed where the panels are built, in `scratch/becca/appris-principal/`.
`isoform_universe.py` now takes the library's selected form as the reference and treats a
same-length form within 3 residues of it as the reference; that stage is pinned (see "The
library pipeline"), so the fix applies to any future universe. Both panel samplers reject a
selected reference form or an ortholog identical to human after the draw, as described in
the next section. CDS resolution keeps its allele tolerance, so it reproduces its record;
a resolution that lands on another member's protein, as CD177's did, is caught downstream.
The guard of record is `dedupe_library_members.py --check`, which fails if any two active
members share a designed protein.

### Panels are redrawn by rule, reproducibly

The isoform and ortholog panels are seeded, class-stratified draws:
`sample_isoform_panel.py` (seed 20261009; 100 observed + 50 prediction-only, by gene class)
and `sample_ortholog_panel.py` (seed 20261011; 100 mouse + 100 cyno, by near/diverged band
and gene class), both in `scratch/becca/appris-principal/`. Their inputs were first read
from one laptop, and some were never committed. Every input is now a dated snapshot on S3,
`library-design/panel-inputs/2026-11-01/`, with a SHA-256 per file in `panel_inputs.json`;
`--snapshot` reads nothing else. With `--original`, both samplers reproduce the published
150 and 200 byte for byte.

A drawn form the library cannot use is replaced by one fixed rule (`panel_draw.py`). The
draw runs unchanged, so every valid pick keeps its place. Each rejected pick is then
replaced by continuing the same seeded draw within its own stratum (tier or band, and
class). A form is rejected when:

| Reason | Picks rejected | Detected by |
|---|---:|---|
| it is its gene's selected reference form | 1 (ACHE P22303-2) | the isoform sampler, from `reference_selection.tsv` |
| its protein is identical to human (correspondence 1.0) | 1 (cyno MMP16) | the ortholog sampler |
| its designed protein duplicates another member's | 1 (CD177 Q8N6Q3) | `redraw_panels.py`, on the translated CDS |
| no transcript encodes it, so it has no designed CDS | 10 (6 isoforms, 4 orthologs) | `redraw_panels.py`: no manifest key |

The last two are only knowable after the draw. `redraw_panels.py` writes them to
`panel_rejections.tsv` and reruns both samplers until a pass adds no rejection; two passes
reached that point. The ten no-CDS picks were never L1 members, which is why the panels
held 144 and 196 members instead of their designed 150 and 200. They fall under the same
rule, so the panels are now 150 and 200. The replacements, each in its rejected pick's
stratum, are logged in `isoform_panel_replacements.tsv` and
`ortholog_panel_replacements.tsv`. The manifest changes are in `panel_changes.tsv`:

| New member | Replaces | Stratum |
|---|---|---|
| l1-p2726 CD58 P19256-3 | ACHE P22303-2 | observed, GPI |
| l1-p2727 LYPD5 Q6UWN5-2 | CD177 Q8N6Q3 | observed, GPI |
| l1-p2728 PRSS8 Q16651-2 | LYPD1 Q8N2G4-3 | observed, GPI |
| l1-p2729 ADGRE1 Q14246-4 | SLC29A1 Q99808-2 | observed, multipass |
| l1-p2730 GPR161 Q8N6U8-4 | SCNN1A P37088-5 | observed, multipass |
| l1-p2731 SLC43A2 ENST00000572801 (62 aa) | SLC12A6 Q9UHW9-6 | observed, multipass |
| l1-p2732 PCDHGA7 Q9Y5G6-2 | LRRC37B ENST00000582117 | observed, type I |
| l1-p2733 SUSD3 ENST00000617293 | NRG4 ENST00000563910 | prediction-only, type III |
| l1-p2734 cyno ECE1 ENSMFAP00000023673 | cyno ATP1B1 Q4R4V5 | cyno near, type II |
| l1-p2735 cyno UNC5C G7P5Y1 | cyno MMP16 G7PC63 | cyno near, type I |
| l1-p2736 mouse EMB P21995 | mouse PTPRJ Q64455 | mouse diverged, type I |
| l1-p2737 mouse NT5E Q61503 | mouse LYPD6B Q9D7F2-1 | mouse diverged, GPI |
| l1-p2738 mouse CADM2 Q8BLQ9-2 | mouse INSR M4TKR7 | mouse near, type I |

Each new member was already an L2 candidate with a resolved CDS. It takes the next free L1
id, and its L2 row is retired `moved_to:<id>`; no id is renumbered. Rerunning
`split_library_tranches.py` would renumber the library, so it must not be used for this.
The paralog retirements are not refilled, because their protein is still in the library.

Three new members had no DeepTMHMM v1 topology. `publish_topology_supplement.py` closes
that gap the same way for any member. It writes the FASTA of every active member without a
v1 row and runs the licensed DeepTMHMM 1.0.24 on Modal (`modal/deeptmhmm_v1_app.py`). The
3-line output becomes rows under the fresh pair
`(topo_2026_11_01_l1_designed, library_designed_orf)`, appended to a new snapshot,
`external://accessible-surfaceome/topology/2026-11-01/topology.parquet` (71,313 → 71,316
rows). `topology.SNAPSHOT_KEY` now points at that snapshot. The run's FASTA, 3-line output
and record (Modal app `ap-mCYrvq2jb3oarooMuLQAwb`) are committed in `topology_runs/` and on
S3. These three rows are not yet in D1 `topology_public`, which the viewer reads.

| Member | v1 call | Placement |
|---|---|---|
| l1-p2731 SLC43A2 (62 aa) | 1 TM helix, N-terminus cytoplasmic, C-terminus extracellular | C-terminal |
| l1-p2733 SUSD3 | 1 TM helix, no signal peptide, N-terminus extracellular | N at residue 1, type III |
| l1-p2734 cyno ECE1 | 1 TM helix, N-terminus cytoplasmic, C-terminus extracellular | C-terminal |

Every active L1 member now has a placement. Cyno ECE1 has no AlphaFold DB model, so it was
folded from its designed CDS with ColabFold 1.6.2, as the other 71 were (754 aa, mean pLDDT
89.1, against 89.3 for human ECE1; under
`predicted-structures/colabfold/2026-11-02/`). The other twelve new members have an AlphaFold
model identical to the designed protein.
SLC43A2's drawn form is a 62-residue fragment; the cohort table admitted it, and the rule
does not second-guess the cohort.

### The library pipeline: gene list to ordered proteoforms

`library_pipeline.py` lists every step from the surfaceome release to a member's tag
cassette, in run order. For each step it gives what the step decides, the script that
decides it, and its inputs. `verify` reruns each step offline from pinned inputs and
compares the result with the published record. The approach is freeze-then-verify (Becca,
2026-11-01): a step must reproduce what the library actually used, so nothing is rebuilt
from today's sources.

Every input is either an immutable `external://` key (a named release of a public source)
or a file in a dated processed snapshot. `library_pipeline_inputs.json`,
`reproduce_isoform_universe_inputs.json` and `panel_inputs.json` record each input's source
key, SHA-256 and origin. `verify` refuses any input whose hash differs.

```bash
uv run --locked python library_pipeline.py plan                      # the stage table
uv run --locked python library_pipeline.py verify                    # rerun every stage
uv run --locked python library_pipeline.py verify --stage gene_list  # one stage
```

| # | Stage | Decides | How it is checked |
|---:|---|---|---|
| 1 | gene list | which genes: the S1E3 filter over surfaceome release 1.3.0, minus `gene_list_exclusions.tsv` | reruns to 2,370; minus EPPIN, it equals the 2,369 library genes |
| 2 | GPI set | which forms are GPI-anchored (NetGPI 1.1 per form ∪ curation) | byte for byte |
| 3 | reference | each gene's reference proteoform | byte for byte |
| 4 | isoform universe → cohorts | every alternate form, its evidence, and its cohort tier | every intermediate summary exactly; `isoform_cohorts.tsv` row for row (see below) |
| 5 | ortholog repick | each gene's mouse and cyno proteoform | byte for byte |
| 5b | universe coverage | which library genes the isoform and ortholog universes could draw from | every recorded ortholog pair is 1:1 high-confidence in Compara 115; coverage table byte for byte (see below) |
| 6 | panels | the 150-isoform and 200-ortholog draws | byte for byte (original draw and rule-based redraw) |
| 7 | rescue arms | 25 expression-rescue arms, their CDS and constructs | byte for byte |
| 8 | giants | truncations of oversized genes and their ORFs | byte for byte |
| 9 | candidate CDS | each candidate's designed CDS | all 7,017 recorded CDS |
| 10 | constructs | domesticated ORFs and p0039 vectors | all 7,017 ORFs and 7,066 vectors (about 30 s) |
| 11 | tranches | L1/L2 membership; duplicate proteins retired | no two active members share a protein |
| 12 | topology | DeepTMHMM v1 per designed protein | every active member has one |
| 13 | tag paths | each member's placement rung | byte for byte |
| 14 | cassettes | each C-terminal member's recoded linkers | every issued batch replays to the recorded encodings |

Each step runs from these pinned sources:

| Source | Version pinned | Used by |
|---|---|---|
| accessible-surfaceome | release 1.3.0 (Zenodo doi:10.5281/zenodo.23445521): deep-dive records, preset table | 1 |
| accessible-surfaceome | commit f9c73691: canonical and isoform FASTAs, DeepTMHMM tables, ortholog cohorts | 3, 4, 5, 5b, 7 |
| Ensembl | release 115: human, mouse and cyno cds and pep; human GTF | 3, 4, 5, 5b, 7, 9 |
| Ensembl Compara | release 115 protein homology dumps (human and mouse files, with Ensembl's CHECKSUMS) | 5b |
| UniProt | release 2026_03: reviewed and isoform tables (snapshot of 2026-10-11), giant features | 3, 4, 5, 8, 9 |
| APPRIS | e115v50 principal isoforms; 2026_06.v50 full scores (TRIFID) | 3, 4, 5 |
| OpenProt | 2.2 (refprots+altprots+isoforms, min 2 peptides, UniProt 2022_06_01) | 4 |
| PeptideAtlas | Human build 607, replayed in the original query batches | 3, 4 |
| GENCODE Ribo-seq ORFs | phase 1 all; phase 2 primary and comprehensive | 3 |
| Expression | deliverome-internal TPM table at 5a1416ea; S1E3 re-derivation 2026-08-04; MyGene responses of 2026-11-02 | 1 |
| DeepTMHMM | 1.0.24 (v1), topology snapshot 2026-11-01 | 4, 12, 13 |
| NetGPI | 1.1 predictions (committed) | 2, 5 |
| Tedman 2025 GPCR screen | supplementary table, 2025-09-19 | 7 |

Steps 4 and 5 run each script at the commit that produced the record
(`reproduce_isoform_universe.py`), because later code changes (allele tolerance, NetGPI as
the GPI set) would give a different table. Two inputs had to be reconstructed. Each
variant used is the only one that reproduced the record:

- The PeptideAtlas cache is replayed in its original four fill stages.
- The D1 GPI subcategory is taken as of 2026-10-11, before EPGN gained a row.

One limit is permanent. `isoform_cohorts.tsv` reproduces in every row and in gene order,
but 420–450 rows (varying with the run's hash seed) come back in a different order within their gene. The
1,527-form topology shortlist also picks a different same-length translation for 78–95
forms. Both orders came from Python's string-hash seed, which the 2026-10-08 run did not
record. The committed `topology_shortlist.fasta` and `.3line` therefore remain the record
of what DeepTMHMM was run on, and the verifier allows exactly these differences and no
others. Any new run of these steps should set `PYTHONHASHSEED`.

**The reference ladder.** Stage 3 tries each rung in turn:

1. MANE and APPRIS agree on the protein sequence.
2. Otherwise, take whichever of the two exists.
3. Where they disagree, use observed evidence: PeptideAtlas build 607 peptides seen in at
   least 10 samples, and GENCODE Ribo-seq ORF translations.
4. Without strong evidence, default to MANE.

Named overrides in `reference_overrides.tsv` take precedence; seven genes have one (six GPI
forms and MICA).

**GPI.** NetGPI 1.1 per form, unioned with curation under the 0.55 confidence rule
(`gpi_set.tsv`), is the GPI set of every stage (Becca, 2026-11-01). Before that, the
selection stages read a curated union (D1 subcategory ∪ UniProt KW-0336 ∪ supplement −
exclusions). The curated union is kept, as `gpi.curated_gpi_genes()`, only to choose which
genes' forms are sent to NetGPI. Within the library, the two sets differ only by GFRAL and
RAET1G. Neither is a NetGPI call: NetGPI finds no anchored form for either (RAET1G is
confidently negative on all 5 forms; GFRAL has no scored form), and both enter the new set
on curation through the 0.55 rule. Switching changed those two genes' `selected_gene_is_gpi`
flag and nothing else: no reference selection changed and no panel draw changed. CFC1,
CRIPTO and GML are the same conflict, curated GPI against a confident NetGPI negative. All
five are placed N-terminal after the signal peptide, so no placement depends on the call.

What freezing the pipeline turned up:

| Finding | Consequence |
|---|---|
| The published `reference_selection.tsv` carried zero PeptideAtlas evidence for all 111 MANE/APPRIS disagreements. `--offline` returned no peptides instead of reading the cache. | With build 607 loaded: 9 genes resolve to MANE on strong evidence, 13 on weak MANE evidence, 7 on weak APPRIS evidence, and 82 have none. No selection changes, because only strong evidence can move a gene off MANE. 20 genes are now flagged for review. The table is republished with the evidence, and `--offline` now reads the pinned cache and fails if any accession is missing. |
| MICA's designed CDS had been patched in place on S3 (`repoint_mica_reference.py`) | The override now reaches it through the pipeline: a fresh reference selection names Q29983, and CDS resolution produces the patched CDS unaided. The rerun also yields `iso:MICA:Q29983`, which is now MICA's reference and therefore a duplicate. |
| `gpi.gpi_genes()` returned an empty set without error when its files were missing | It now raises. |

**Adding a gene.** Members are appended, never renumbered, so a gene admitted later
cannot shift any existing member. List the gene in `gene_list_additions.tsv` with the
decision, the date and the surfaceome release that admitted it, then run
`add_library_genes.py`. It runs the reference, CDS and construct steps (3, 9 and 10) over
the list from the same pinned inputs. It fails unless every recorded reference row and CDS comes back unchanged.
The new member gets the next free L1 id, and its designed protein must not repeat an active
member's. Its designed protein also needs a v1 topology row; without one, the script stops
and names `publish_topology_supplement.py`. An added gene is a library reference only:
the isoform and ortholog panels draw from the pinned universe, which predates it. The
pre-tranche `candidate_manifest.tsv` and `candidate_vectors.fasta.gz` keep their positional
ids, a namespace the tranche split abandoned, so an addition's vector goes to
`library-candidates/added_vectors.fasta` under its final id.

```bash
uv run --locked python add_library_genes.py              # dry run
uv run --locked python add_library_genes.py --write      # append locally; then commit
uv run --locked python add_library_genes.py --publish    # upload the S3 records
uv run --locked python l1_tag_paths.py && uv run --locked python l1_swap_risk.py
```

SLC67A1 (formerly SLC22A18; UniProt Q96BI1) is the first addition (Becca, 2026-11-02). Surfaceome
release 1.3.0 calls its canonical 12 TM helices where the September cohort had 0, which
moves it from T3 to T2. It is the only gene the S1E3 filter adds over 1.3.0, and none is
removed. As `l1-p2739`:

| Property | Value |
|---|---|
| Reference | MANE and APPRIS agree on ENST00000649076, 424 aa, identical to Q96BI1 |
| Construct | 5 domestication edits; 8,842 bp vector |
| Topology | DeepTMHMM v1 (`topo_2026_05_16`, `human_canonical`): TM, 12 helices |
| Tag rung | internal; longest inter-TM extracellular run 16 aa at residue 262 |
| Swap risk | worst partner shares 27 nt (IGSF8) |
| Structure | AlphaFold model of Q96BI1, residue-identical |

None of the other 2,368 reference rows, 7,016 CDS or panels moved.

**Rebuilding a gene whose reference-genome copy is broken.** Eight library genes had no
member because GRCh38 carries a loss-of-function copy (Ensembl 115 biotype
`protein_coding_LoF`): none has a MANE transcript, the ladder's APPRIS pick was the broken
transcript, and CDS resolution refuses a transcript that does not translate to an intact
protein. Three of them have a transcript encoding the intact UniProt protein, residue for
residue, on an alternative-haplotype or patch scaffold of Ensembl 115. A named override with
basis `reference_allele_lof` in `reference_overrides.tsv` points each one's reference at that
protein; `add_library_genes.py` then builds it like an addition, and fails unless every other
recorded reference row and CDS reproduces. The overrides file is pinned as a new input
snapshot, `library-design/pipeline-inputs/2026-11-04/reference_overrides.tsv`.
Both halves of each override are recomputed from the pinned releases by
`scratch/becca/appris-principal/check_reference_allele_lof.py`, which fails unless every GRCh38
primary-assembly transcript of the gene is `protein_coding_LoF` with internal stops and the
transcript built from translates exactly to the named UniProt 2026_03 sequence; its per-transcript
table is `reference_allele_lof_evidence.tsv` (primary copies: KIR2DS4 3 internal stops,
OR12D1 8, OR5AL1 10).

| Member | Gene | Defect in GRCh38 | CDS from | Tag rung | Structure |
|---|---|---|---|---|---|
| `l1-p2740` | KIR2DS4 | 297-aa loss-of-function form | ENST00000617825 (alternative KIR haplotype), 304 aa = P43632 | N after SP | AlphaFold model of P43632, identical |
| `l1-p2741` | OR12D1 | frameshift | ENST00000703726 (alternative MHC haplotype), 320 aa = P0DN82 | N at residue 1, multipass | AlphaFold model of P0DN82, identical |
| `l1-p2742` | OR5AL1 | frameshift | ENST00000710330 (patch scaffold), 329 aa = P0C617 | N at residue 1, multipass | AlphaFold model of P0C617, identical |

All three have a published DeepTMHMM v1 row (`topo_2026_05_16`, `human_canonical`), need no
domestication edit and are not C-terminal, so no cassette encoding is issued or moved.
KIR2DS4 shares 125 nt of 3' coding sequence up to its barcode with KIR2DS2 (`l1-p1055`), a
paralog, which is over the 100-nt swap bar for the pair. The other five (FCGR2C and OR8J2,
with a premature stop; OR10AC1, OR52P1 and OR4Q2, frameshifted) have no intact Ensembl 115
transcript and remain without a member.

### Open questions

| Question | Why it matters | What would resolve it |
|---|---|---|
| The viewer's gene list is not the library's | The library is the 2,368-gene root plus SLC67A1. The viewer's `fg-library.json` holds 2,369 different genes: EPPIN, which was removed, and no SLC67A1. Anything computed from the viewer copy inherits both differences | Re-export `fg-library.json` from root + `gene_list_additions.tsv` |
| Cargo size on the internal rung | At 19 aa the ALFA-only loop cassette fits SLC16A9/12/14's 16-residue runs; the 51-aa tandem does not, so the cargo decision determines whether a whole paralog block is in the library | Decide whether internal-rung members carry SpyTag003, then calibrate `MIN_LOOP_LEN` at that size |
| `MIN_LOOP_LEN` and `LOOP_SS` | They change the top-ranked pick for 56% and 72% of members respectively, so they decide which site 515 members get | Wet data from a coil-run ladder; the benchmark corpus is 78% small epitope tags and cannot settle it |
| GPI members needing an internal site | A GPI single-pass protein needs a loop but reaches the C-terminal rung on topology alone, so it is absent from the 515 | Build the GPI internal set from the 132 `reference_is_anchored` genes |
| L2's 1,155 members with no model | Each is unscoreable until it is folded, and L2 has not been checked for proteoform mismatches at all | Fold the gap (~80 GPU-hours at the L1 rate), then run `check_proteoform_match.py` over it |
| 199 L1 members whose model is a different protein | Their tag sites would be numbered against a sequence the construct does not contain, and nothing flags it | Fold each one's designed CDS, or verify and apply the offset for the 26 contained cases |

### Next steps

1. Re-publish the unstamped keys under this prefix so each can be traced to a commit, or
   record why they cannot be.
2. Settle the ortholog panel size against `ortholog_panel_summary.json`.
3. Decide the internal-rung cargo (ALFA-only at 19 aa, or the 51-aa tandem), then calibrate
   `MIN_LOOP_LEN` at that size.
4. Fold the GPI branch into `derive_tag_paths.py`, which is topology-only today and so
   routes GPI single-pass members to the C-terminal rung.
5. Give the disorder lane a secondary-structure gate, or rank by one.
6. Measure L2's structural coverage exactly, then decide whether to fold its gap (~85 GPU-hours at the L1 rate).
7. Publish the predicted structures to R2 with a D1 index row, so the viewer can reach them as it reaches the Schweke models.
8. Fold the 199 L1 members whose AlphaFold model is a different proteoform, from their designed CDS rather than their canonical accession.

### Decision log

| Date | Was | Now | Why |
|---|---|---|---|
| 2026-10-09 | Isoform panel drawn uniformly at random | 100 observed + 50 prediction-only, balanced across type I, type II, GPI and multipass | A uniform draw loses the rare architectures, which are the ones whose placement is least certain |
| 2026-10-11 | Library constructs unnamed | `l1-pNNNN` | Readability against the `p00NN` panel series |
| 2026-10-11 | GPI called per gene from UniProt keyword | NetGPI per form, cross-checked against UniProt with documented blind spots | The anchored species is not always the reference form — ACHE's canonical is non-anchored at 0.993 while P22303-2 is anchored |
| 2026-10-12 | Internal-rung set reported as 377, from `internal_site_only.tsv` | 515, from both termini cytoplasmic in v1 topology | The 377 held back 176 members as topology-undetermined; all 176 have determined orientations in v1, 116 of them both-cytoplasmic, so nothing was waiting on DeepTMHMM2 |
| 2026-10-15 | Tag rung assigned per canonical gene (2,367) | Per L1 **member** (2,700), on each member's own ordered-proteoform v1 topology | We order members, not genes, and a member's topology is not its gene's — SLC9D1's ordered isoform has 3 TM against the canonical's 13. Internal-site worklist 515 genes → 566 members; restricting to `library_reference` reproduces the per-gene numbers |
| 2026-10-15 | `internal_site_only.tsv` / `internal_site_undetermined.tsv` / `INTERNAL_SITE_ONLY.md` and their generator | Deleted; `TAG_PATHS.md` + `l1_tag_paths.tsv` replace them | Their 377/176 split reads UniProt's missing topological-domain annotations as unknown orientation, so it misleads on sight. Record kept in git history at `8eab499` |
| 2026-10-14 | Structural coverage measured as "does a model exist" | Measured residue-exact against the designed CDS | Length agreement is necessary and not sufficient: 23 of 200 mismatched members are the same length as their model with a different sequence |
| 2026-10-14 | Members without an AlphaFold DB model treated as a reporting gap | Folded with ColabFold and published as our own predictions | A member with no model cannot be scored by any structural lane, so the gap blocks the library rather than annotating it |
| 2026-11-01 | Library steps run from local caches and hand-uploaded tables | One pipeline (`library_pipeline.py`): 7 stages verified from a checksummed S3 snapshot, 5 pinned with the reason | Reproducibility (Becca); freeze-then-verify, so the library that exists is the one reproduced |
| 2026-11-01 | Selection stages read the curated GPI union; tag placement read NetGPI | NetGPI per form (`gpi_set.tsv`) for every stage | One definition (Becca). Within the library it adds GFRAL and RAET1G, on curation rather than NetGPI, and changes no selection or panel |
| 2026-11-01 | `reference_selection.tsv` published with no PeptideAtlas evidence | Republished with build 607 evidence from the pinned cache | `--offline` had discarded the evidence; no selection changes, 20 genes flagged for review |
| 2026-10-29 | Every codon of the tag cassette recoded per member, ends unique in their last 9 nt | Only the GGGGS linkers recoded; SpyTag003 and P-ALFA-P on p003's codons; the terminal 15-nt GGGGS unique per member | Tags stay on the validated encoding (Becca); recoding them gives no swap benefit (36 nt either way), and the 9-nt end rule runs out at 60 of 216 C-terminal members. The 25-construct order is not changed |
| 2026-11-01 | Panel picks the library could not use left out (panels at 144 and 196) | Panels redrawn from a pinned S3 input snapshot; each rejected pick replaced by continuing the same seeded draw in its stratum; panels at 150 and 200 | Reproducibility (Becca): a hand-picked or separately seeded replacement cannot be regenerated by the panel code. The snapshot reproduces the original panels byte for byte |
| 2026-10-29 | L1 members accepted whenever gene, accession or transcript differed | One active member per translated designed protein; duplicates retired with ids kept | 7 L1 and 48 L2 members duplicated another member's protein (ACHE, CD177, MMP16, OPN1MW2/3, OR2A42, UPK3BL2 in L1) |
| 2026-11-02 | Gene list frozen at 2,368 | Frozen root plus dated additions (`gene_list_additions.tsv`), appended by `add_library_genes.py`; SLC67A1 added as `l1-p2739` | Becca: adopt the one gene surfaceome 1.3.0 adds. Appending, not re-deriving, leaves every existing member in place, and the script asserts it |
| 2026-11-02 | Constructs stage pinned (full rerun is hours) | Every construct and vector re-verified (about 30 s) | A full rerun takes 30 s. It found the S3 constructs record still held MICA's pre-repoint ORF, now republished |
| 2026-11-02 | Gene list, isoform universe, ortholog repick, rescue arms and giants pinned as unreproducible records | All rerun offline from named releases on external:// (surfaceome 1.3.0, Ensembl 115, UniProt 2026_03, APPRIS 2026_06.v50, OpenProt 2.2, PeptideAtlas build 607); the universe and repick at their producing commits | Reproducibility years out (Becca). OpenProt 2.2 and the APPRIS full scores were still downloadable, so not lost. EPPIN's removal is now a recorded exclusion, not a hand edit |
| 2026-11-02 | Linker encodings pinned as issued | Each encoding records its issue batch; verify replays the batches from scratch | Shows the 218 encodings follow from the seed and the issue order alone |
| 2026-11-02 | Cyno ECE1 (l1-p2734) listed as needing no fold | Folded from its designed CDS with ColabFold | The coverage baseline had counted its accession as resolved, so the fold-completeness test passed with the gap open |
| 2026-11-04 | KIR2DS4, OR12D1, OR5AL1 without a member (GRCh38 carries a loss-of-function copy) | Built as l1-p2740..2742 from their intact alternative-haplotype transcripts via `reference_allele_lof` overrides | Becca: the intact protein is in Ensembl 115; CDS resolution had searched only by the broken copy's length |
| 2026-11-03 | Snorkel cited to PMID 24391846; GFRAL and RAET1G described as NetGPI additions | Snorkel cited to Brown et al. 2013, PMID 24023844; both described as curated calls NetGPI does not support | The old PMID is an unrelated docking paper (Hsin et al. 2013); NetGPI finds no anchored form for either gene. Found by the reviewer-package review, with re-measured counts (191 GPI members, 199 different models) |
| 2026-11-04 | 13 designed CDSs ended on a sense codon: the Ensembl 115 record has no stop (l1-p1054..1056 KIR2DL5A, KIR2DS2, KIR2DS5; L2 l2-p0846 KIR2DS2 isoform and nine cyno orthologs) | `resolve_candidate_cds.py` appends TAA when the record ends without a stop and marks the entry `stop_added`; orf_nt, vector_bp and insert_nt grow by 3 for those rows. Proteins and domestication edits are unchanged | Becca. Without a stop the ORF reads through into the vector. TAA is the stop the C-terminal constructs already use; the native codon cannot be read from the pinned inputs (four KIR records sit on alternate-haplotype scaffolds the primary GTF omits). Previous S3 versions: candidate_cds.json `ftYKdsVZYDD8M49CEt4rDJMn.x.QmKkn`, candidate_constructs.json `Nkbr70dq6lfjWkb.ASYGWZpU8yRocj_Y`. Found by the N-terminal reviewer package's frame check |
| 2026-11-05 | Isoform and ortholog panels described as drawn from library genes with a form or ortholog | Recorded as drawn from the April 2026 `surfaceome_expressed` query set (1,775 of 2,369 library genes); coverage measured against Compara 115 by `check_universe_coverage.py` | Becca: record the gap and keep the panels. The 418 mouse / 416 cyno eligible genes were never queried because later Compara pulls went to D1, not to the f9c73691 files the universes read; every recorded pair is one-to-one high-confidence in Compara 115 |
| 2026-10-13 | Walton ≤100 bp tract bar applied library-wide | Bar applies to the pool and the shared cassette, not to isoforms and orthologs | Median canonical-isoform overlap is 1,071 nt with zero of 2,150 pairs under the bar; applying it would remove the panels entirely |

## Provenance

- p003 constants read from `processed://constructs/ordered/p003.gb`. This was a local-only
  file when the note was first written; the map was published to the construct registry on
  2026-10-06, so the gap it described — and the matching one in `s1e2/enzyme_scan.py` — is
  closed for p003.
- Tedman construct design: Nat Commun 2026, [doi:10.1038/s41467-026-76564-7](https://doi.org/10.1038/s41467-026-76564-7)
  (via PubMed, PMC13575135); screen values from
  `deliverome-external://tedman-gpcr-surface-screen/2025-09-19/`.
- Barcodes: `analysis/sprint1/s1e2-fragmid-vector-design/barcodes_20nt_30k.tsv`.


## TFRC cassette-order arm

`TFRC_cterm_swap` is the only arm built here. TFRC is type II with an extracellular
C-terminus and is **already built and working as p003**, so this is a one-variable
change against a known reference. The ORF is lifted verbatim from the assembled p003 map
(4831..7110 — already domesticated, 0 Type IIS sites), so the swap differs from p003 in
the cassette order and the barcode and in nothing else. `tag_arms.py` asserts the
unswapped arm reproduces p003's insert region byte-for-byte; it passes.

```
p003        …ORF | GGGGS SpyTag003 GGGGSGGGGS P-ALFA-P GGGGS |   ALFA at 33-45
cterm_swap  …ORF | GGGGS P-ALFA-P GGGGSGGGGS SpyTag003 GGGGS |   ALFA at  7-19
```

| Arm | Tag order | Insert | +shuttle | Barcode |
|---|---|---:|---:|---|
| `TFRC_cterm` | SpyTag003 → ALFA | 2529 nt | 2545 | `TAGCACTGAGTCAGATCTGA` — **p003, exists** |
| `TFRC_cterm_swap` | **ALFA → SpyTag003** | 2529 nt | 2545 | `CTCAAGTGCAGCGGCCATAT` — **the only new build** |

**What it is not.** It does not test whether ALFA needs more structure: the proline
flanks that isolate its helix (§4.7, Vicidomini 2024) are present in both orders, so
helicity is identical by construction. It also moves both tags at once, so a difference
cannot be attributed to one without follow-up.

**Nothing predicts a difference.** Both tags are already validated in the p003 order —
ALFA distal by **S2E9** (anti-ALFA nanobody on TfR1-Spy/ALFA: WT at the camera floor,
expressing wells 51–73× the WT signal span, clear membrane outlines) and SpyTag proximal
by **S2E10** (SpyCatcher–CF488A/CF647 conjugates labelling the same L-0004 line). So
this is a robustness check: if the arms match, cassette order stops being a
consideration for future terminal designs.

### ⚠️ Dropped, and still open

A five-arm study including **SLC20A2 N- versus C-terminal** placement was designed and
then cut (Becca, 2026-09-26). **The N-versus-C placement question is therefore
unresolved**, and it is the one the panel cannot answer for itself: TFRC is its only
member with an extracellular C-terminus, so the other fourteen are being built
N-terminally on topology alone, with no head-to-head measurement behind it.

If it is revived, SLC20A2 (PiT2) remains the screened choice — the only FG-library gene
besides SLC20A1 with both termini extracellular, monomeric, with observed overexpression
surface localisation and live-cell surface evidence, and SLC20A1's dual-terminus call is
uncorroborated. Note the N-terminal arm needs ALFA alone (15 aa) or an Igκ leader, not
the bare 51-aa cassette: SLC20A2 has no signal peptide and only 5 extracellular residues
before TM1, so TM1 is the signal-anchor, and the cassette takes the segment it must
translocate from 5 aa / net −2 to 56 aa / net +1 — against the positive-inside rule.

## Tag sites — recorded here, not written upstream

The accessible-surfaceome deterministic `tag_sites` pipeline is still in progress, and
its `viewer/public/tag-sites/*.json` store is the input to the public D1 sync. We did
**not** commit or publish our runs there: the generated JSON was read, the sites copied
into `tag_sites_provenance.json` with the source commit, exact command and input list,
and the untracked runs removed from that repo. Nothing was written to D1 (dry-run only).

⚠️ The upstream JSON is sorted **lexicographically by `site_id`**, which is not a
ranking — see the SLC2A4 note in `build_loop_inserts.py`.


## Barcodes and provisional P-numbers

Barcodes come from the running registry, **not** from the raw pool:

`processed://sprint1/s1e2-fragmid-vector-design/barcode_assignments.tsv`

Its header says to reserve from the pool and append there before assigning to a new
construct. At the time of reservation the registry held pool indices **1, 3, 4, 5, 6, 8**
(p0016–p0021) plus the pre-pool TFRC barcode, with **2 and 7 deliberately skipped**
("boundary-enzyme or splice-donor motif"). Highest P-number was **p0021**.

⚠️ **A first pass drew from pool index 1 upward and collided with seven of those.**
`reserve_barcodes.py` now reads the registry and refuses anything already listed.

Two screens are applied on top of "not already used", both reconstructed from why 2 and
7 were skipped, and both evaluated in the **real junction context** (forward primer site
+ barcode + reverse primer site + the Esp3I shuttle) rather than on the barcode alone:

- **boundary enzyme** — no BbsI / Esp3I / BsaI site created across a junction;
- **splice donor** — no `GT[AG]AG`. Pool index 7 is `GGTAAGACGCGCCTGAGTAG`, opening with
  the canonical `GGTAAG` 5′ splice site; in a lentiviral transcript a cryptic donor there
  would splice out everything downstream. The s1e2 generator's filter stack does **not**
  screen for this, so it has to be applied at reservation.

### The 30k pool is clean — the skips are context effects, not defects

The generator's filter stack screens BbsI / Esp3I / BsaI **inside the barcode, both
strands**, and the pool is clean against that. What it cannot do is screen a site created
**across a junction**, because the flanking sequence was unknown when the pool was built.
It also does not screen splice donors at all — that criterion is simply not in the stack.

Both registry skips are explained by those two gaps:

| Index | Inside the barcode | Why it was skipped |
|---|---|---|
| 2 `AAGGAACTCGACGCATAACG` | clean — no Type IIS, no splice donor | ends `…ATAACG`; a 3′ flank starting `TCTC` completes **`CGTCTC` (Esp3I)** |
| 7 `GGTAAGACGCGCCTGAGTAG` | contains **`GTAAG`** | intrinsic splice donor, *and* a `GAGAC` 5′ or `AAGAC` 3′ flank completes an Esp3I site |

So index 7 is bad everywhere — it carries the donor in its own sequence. **Index 2 is
only bad next to a `TCTC` flank**, which is the p0016–p0021 context, not ours: our 3′
flank is the reverse primer site, `GGATCTGG…`. That is why the screen here passes it, and
it is the correct result rather than a miss. The registry's skip is global and therefore
conservative — fine, and we honour it, but it is not evidence the pool is dirty.

The operative lesson is that **a barcode is only clean relative to the flanks it lands
next to**, so junction-context screening has to happen at reservation, per construct.

| Plasmid | Gene | Arm | Pool index | Barcode |
|---|---|---|---:|---|
| `p0022` | HEG1 | N-term | 9 | `AGCAATCTCCACGGCACTCG` |
| `p0023` | ERBB2 | N-term | 10 | `GCCTATAACGGAGACAGAAT` |
| `p0024` | AXL | N-term | 11 | `GGTGGTTCTTGTAGGAGCAT` |
| `p0025` | SLC9A3 | N-term | 12 | `AGGAATCAGTGCTAACAGAC` |
| `p0026` | VSIG10 | N-term | 13 | `AGATCAACCTACAGCAACCG` |
| `p0027` | NECTIN4 | N-term | 14 | `GAGGCACGGAACTAGCAATT` |
| `p0028` | FOLR1 | N-term | 15 | `CCAGACAAGAGGCGACCAGA` |
| `p0029` | PTPRC | N-term | 16 | `CTCCGCTCTTATCTCTGAAC` |
| `p0030` | SLC2A4 | loop E67 | 17 | `GCCGGAAGTTGCGTGGTATT` |
| `p0031` | SLC2A4 | loop V319 | 18 | `CTCAAGTGCAGCGGCCATAT` |
| `p0032` | MS4A1 | loop P79 | 19 | `ACTGGCAAGCAACGTGTCGA` |
| `p0033` | MS4A1 | loop I162 | 20 | `TACTTCTGTGTGCTTGTTGC` |
| `p0034` | **AVPR2** | N-term, leader-free, 51-aa tandem | 21 | `ATCACGCGACGAGAGAAGCT` |
| `p0035` | TFRC | C-term, p003 order | 22 | `TGTCTCTCTCAGCTTCTCTC` |
| — | — | *index 23 retired for a splice donor* | 23 | — |
| `p0040` | TFRC | outside this batch | 24 | `CAGAGCGTCAGGTAGAACAT` |
| `p0037` | **CNR1** | N-term, leader-free, 51-aa tandem | 25 | `ATAGTCAAGAGCGCACCAGA` |
| `p0038` | **CCR5** | N-term, leader-free, 19-aa ALFA | 26 | `ACAGTTAGTCGAGAAGCGCA` |

Three rows above were stale until 2026-10-18 and are corrected here: `p0031` read
`loop E193` and `p0032` `loop N171`, both retired on 2026-10-13, and `p0030`/`p0033`
carried the pre-2026-10-15 sites. `p0037`/`p0038` read pool index 23/24, which was
pre-skip numbering — `barcode_reservations.tsv` records 25/26 and notes that **index 23
is retired for a splice donor**. Index 24 belongs to `p0040`, so there was never a
double-assignment; the stale file only looked like one.

**`p0034`, `p0037` and `p0038` are recycled slots.** `p0034` was never built and has no
sequence anywhere, and the TFRC linker 3x/4x arms were cut from the batch and parked at
`processed://constructs/scratch/tfrc-linker-{3x,4x}.gb`. Their three barcodes now carry
the leader-free AVPR2, CNR1 and CCR5 arms — the GPCRs Brenna's panel prioritised, none of
which was in the batch. Assignments live in `barcodes_gpcr.tsv`; the parked rows in
`barcodes_tfrc.tsv` record the release, because the scratch maps still physically contain
those barcodes and reviving an arm therefore needs a fresh one.

**P-numbers are provisional** — assigned here in flight after p0021, not yet registered.
`barcode_reservations.tsv` holds the rows to append to the registry; **they have not been
appended**, so another design drawing from index 9 upward would collide until they are.


## Pool capacity under the two screens

Applied to all 30,000 (`reserve_barcodes.py` screens, our flanks):

| Screen | Fails | Scope |
|---|---:|---|
| intrinsic splice donor `GT[AG]AG` | **1,494 (4.98%)** | context-**independent** — bad in any construct |
| — stricter `GGT[AG]AGT` | 119 (0.40%) | |
| junction Type IIS, our flanks | 94 (0.31%) | flank-**dependent** — must be re-run per construct |
| **either** | **1,582 (5.27%)** | |
| **survivors in our context** | **28,418 (94.73%)** | |

Removing members never breaks the pool's pairwise Levenshtein-3 / Hamming-4 guarantee —
a subset of a distance-separated set is still distance-separated — so a filtered sub-pool
stays valid without regeneration.

**⚠️ This puts the pool below the library target.** s1e2 §1 sets library size at *up to
30,000* and the pool is exactly 30,000. Dropping the 1,494 intrinsic splice donors leaves
**28,506**, so the pool no longer covers a full-size library. Regenerating with a
splice-donor filter added to `generate_barcodes.py`'s stack is the fix; the stack
currently screens GC, homopolymer, ATG, polyA, hairpin, miRNA seed, Type IIS and
distance, but **not** splice sites.

### How skips are tracked today — and the gap

There is **no standing analysis**. The registry is a flat TSV whose skip record is a
single hand-written comment line naming indices 2 and 7. Nothing re-screens the pool, and
nothing recomputes when flanks change. The only systematic screening is what
`reserve_barcodes.py` does at reservation time, per construct, in that construct's own
junction context.

The two failure modes want different treatment:

- **splice donor is global** — screen once, remove from the pool (or regenerate with the
  filter in the stack). Re-checking it per construct is wasted work;
- **junction Type IIS is local** — cannot be baked into the pool, because it depends on
  flanks that differ per construct. Index 2 is the worked example: clean in itself, fatal
  only beside a `TCTC` 3′ flank.


## Backbone: puro only, no mScarlet3

All of these go into a **puro-only** backbone to mimic the library (Becca, 2026-09-26).

**Nothing in the insert designs changes.** mScarlet3 sits in the *backbone*, downstream
of the IRES (p003: `IRES2 → mScarlet3 → P2A → PuroR`). Our inserts stop at the barcode
amplicon, so all fifteen are already backbone-agnostic — this is a destination-vector
choice, not a redesign. Verified: the downstream Esp3I shuttle junction is clean in both
backbones (1 Esp3I on the minus strand, 0 BbsI, 0 BsaI), despite their different 3′
flanks — `GGCGCGCC…IRES2` in p003 versus `GCTTTAAGG…` in P008.

**`p0035` exists because of this.** p003 carries mScarlet3, so it is not backbone-matched
to a puro-only swap arm; comparing `p0034` against p003 would confound cassette order
with backbone. `p0035` is the p003 cassette order rebuilt in the puro-only backbone, and
`tag_arms.py` still asserts its coding region reproduces p003's byte-for-byte. So the
order comparison is `p0035` vs `p0034`, both puro-only, differing only in cassette order
and barcode.

### The puro-only intermediate: IRES2 → PuroR

**Decided (Becca, 2026-09-26): IRES2 → PuroR, matching p003 as closely as possible.**
P008 already provides a puro-only vector, but it drives PuroR from its own SV40 promoter
as a separate transcription unit, which **decouples selection from ORF expression** — a
cell can survive puromycin while transcribing the ORF poorly. p003 couples them through
the IRES, and this keeps that.

A single clean deletion:

```
p003    …IRES2 · linker·ATG | mScarlet3 · linker · P2A | PuroR · TAA · WPRE…
                              └───── delete 7969..8721 ─────┘
(this)  …IRES2 · linker·ATG | PuroR · TAA · WPRE…
```

| | |
|---|---|
| Deletion | **7969..8721**, 753 nt (mScarlet3 + linker + P2A) |
| Size | 10,759 → **10,006 bp** |
| New junction | `…GCCGCGCCATG\|ACCGAGTACAAGCCCACG…` |

**No initiator needs adding.** PuroR in p003 begins `ACC GAG TAC` (T-E-Y) — it has no ATG
of its own, because it is translated as a P2A downstream product. But the IRES linker
already carries an in-frame `ATG` at 7966..7968, which was mScarlet3's upstream Met.
Deleting exactly 7969..8721 joins that ATG straight onto PuroR's second codon and
regenerates the **native** PuroR protein, `MTEYKPTVRLATRDDVPRAVR…TWCMTRKPGA`.

Because nothing 5′ of 7969 is touched, the **IRES2-to-initiator spacing is bit-for-bit
what p003 already validated** — that spacing is what sets IRES initiation efficiency, so
leaving it alone is the point of cutting here rather than anywhere else.

Checks, all passing in `puro_only_backbone.py`: native PuroR N- and C-termini, single
clean reading frame with one stop, no Type IIS site created across the junction, and
everything 5′ of the deletion byte-identical to p003.


## Enzyme scope: BbsI + Esp3I + BsaI, on every ORF

**Three enzymes, deliberately.** BbsI is the Fragmid assembly enzyme, Esp3I is the
shuttle (and what Tedman used), BsaI is scrubbed because it was never cleared from the
backbone. `seqtools.ENZYME_SITES` holds exactly these three and every ORF here is
domesticated against all three, to **zero** sites.

`s1e2/enzyme_scan.py` scans a wider set — AarI/PaqCI and SapI as well — but that script's
job is surveying *backbones* for a free shuttle enzyme, not domesticating ORFs.
**AarI and SapI are out of scope for ORF domestication and are not screened here.** Do
not read their presence in `enzyme_scan.py` as a gap in this design.

### All fifteen ORFs, including the Tedman five

The Tedman receptors were initially listed as "reuse verbatim". **That was wrong for our
backbone.** Tedman silenced **BsmBI/Esp3I only**, so every one of his ORFs still carries
BbsI and/or BsaI sites — and BbsI is our assembly enzyme:

| Gene | BbsI | Esp3I | BsaI | Edits to clear |
|---|---:|---:|---:|---:|
| CNR1 | 1 | 0 | 0 | 1 |
| AVPR2 | 1 | 0 | 0 | 1 |
| ACKR3 | 3 | 1 | 2 | 7 |
| CCR5 | 2 | 1 | 0 | 4 |
| OR51E2 | 1 | 0 | 0 | 1 |

So they are designed here like any other ORF. **Tedman's design is not compromised by
this**: domestication is synonymous-only, so the protein and his tag placement are
unchanged, and he kept native codons outside his own BsmBI edits for the same reason we
do.

Full domestication table for all fifteen is in `orf_design.tsv`.
