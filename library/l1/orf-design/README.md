# s4e9 — orf-design-fg-panel

For continuation of the ortholog CDS work, start with the [agent handoff](AGENT_HANDOFF.md).
It identifies merged work, durable artifacts, remaining adoption decisions, and compute approvals.

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

**[`REBUILD.md`](REBUILD.md) is the procedure for rebuilding the L1 order from zero, and
[`CODE_MAP.md`](CODE_MAP.md) is the map of this folder.** It lists every pipeline stage with
the script it runs and the tables it writes, then every script with what it does and a
status: `stage`, `called` (another script runs it), `documented`, `test`, or `one-off` (ran
once; its output is the record and nothing reruns it). One-off scripts are kept, not
deleted, because each is the only record of how a committed or published output was made.
[`code_map.py`](code_map.py) regenerates the map, and `tests/test_code_map.py` fails when a
script is added or removed without it. The reviewer brief is rendered to PDF and Word by
[`reviewer-package/render_brief.py`](reviewer-package/render_brief.py).

Every external source a script reads (UniProt, Ensembl, MyGene, AlphaFold DB, TmAlphaFold,
NetGPI, D1, the licensed predictors and the Modal images) is listed in
[`external_inputs.tsv`](external_inputs.tsv) with its version, its saved copy and that copy's
SHA-256 ([`external_inputs.py`](external_inputs.py)). `tests/test_external_inputs.py` fails if a
script makes a network call that is not listed there.

**Tables of record for L1** are the `constructs_l1_<rung>.tsv` / `.fasta` arm tables and the order of record named in `synthesis_compliance.py`; `constructs_nterm.tsv`, `constructs_cterm.tsv` and `constructs_loop.tsv` are the nine-ORF FG panel's, and `reviewer-package/review-repo/` is an export copy for the private review repository, not a source. The table below covers the original nine-ORF FG panel scripts:

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
than thin. (Since 2026-11-12 RGMB's named reference, Q6NW40, sits on N after SP, so the reroute
touches BST2 only.)

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
[METHODS_library_design.md](../../../library/l1/tag-site-architecture/METHODS_library_design.md),
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

The ortholog panel is 200 and the isoform panel 150 (134 genes).

| Panel | Members | Composition | Key |
|---|---|---|---|
| Isoforms | 150 | 100 observed + 50 prediction-only; type I 61, multipass 51, type II 17, GPI 16, type III 5 | `library-design/isoform_panel.tsv` |
| Orthologs, L1 (adopted 2026-11-10) | 200 | the CDS-eligible Compara 116 draw with 11 correspondence replacements; see "The L1 ortholog panel" below | `ortholog_panel_adopted.tsv` |
| Orthologs, original draw (retired from L1) | 200 | mouse 100 + cyno 100, split near/diverged at 50 per band, seed 20261011, minimum floors on type III, GPI and type II; type I 66, multipass 72, GPI 28, type II 26, type III 8 | `library-design/ortholog_panel.tsv` |
| Ortholog replacement, Compara 116 | 200 | same seed, species/band sizes and class floors; type I 54, multipass 81, GPI 28, type II 28, type III 9 | `ortholog_panel_compara116.tsv` |

**The recorded ortholog panel draws from April's gene set.** Its pool uses
accessible-surfaceome files frozen at commit f9c73691 (2026-04-27), whose human query
set is deliverome-internal's expression-filtered `surfaceome_expressed` list (2,359
proteins, 2,355 gene symbols). Only 1,775 of the current library's 2,369 genes are in
that set. Later Compara pulls went into surfaceome D1 rather than these frozen
files, leaving 594 library genes unqueried for this ortholog pool.
`check_universe_coverage.py` measures the recorded ortholog pool from pinned
Compara 115 and April inputs and separately measures selected isoform cohort
membership; `universe_coverage.tsv` contains one row per current library gene.

| | Mouse | Cyno |
|---|---:|---:|
| Library genes with a one-to-one high-confidence ortholog in Compara 115 | 1,870 | 1,848 |
| … in the ortholog list the panel drew from | 1,452 | 1,432 |
| … never queried (all outside the April set) | 418 | 416 |
| Recorded pairs that are not one-to-one high-confidence in Compara 115 | 0 of 1,972 | 0 of 1,927 |

The isoform universe uses the committed library gene list. All 451 library genes
with a non-canonical UniProt isoform but no selected isoform cohort member have
scored universe rows, including the 92 outside April's set. Their absence is
explained by the existing cohort rules: 399 genes have only alternates below its
evidence criterion; the other 52 include topology exclusions. A form passes the
evidence criterion when it has a diagnostic peptide, positive OpenProt translation
score, or APPRIS TRIFID score above 0.5; evidence-passing forms can still be
excluded for lost surface anchoring. Insufficient evidence is not evidence that an isoform is absent.
`audit_isoform_coverage.py` records every form decision and fails on an unexplained
absence; [the audit note](ISOFORM_COVERAGE.md) gives the full accounting and
historical replay limits. The current panels, member IDs and constructs are unchanged.
The replacement ortholog panel draws over the full current library gene list;
its member-ID allocation and construct regeneration follow selection.

The [CDS eligibility review](ORTHOLOG_CDS_POLICY.md) gates the full transcript pool
before protein ranking and supplies a 200-member review panel with verified native
CDSs and no pending predictions. Its [panel-change ledger](ortholog_cds_panel_changes.tsv)
distinguishes the full redraw from sequence correction. The policy and prediction
records are saved. A [correspondence review](ORTHOLOG_CORRESPONDENCE_REVIEW.md) compares
every panel form with its human library protein (signal peptide, TM helices, topology, GPI,
Pfam domains): 183 forms keep the human architecture, 11 lack a human feature and 6 need a
decision. The reviewed panel is the L1 ortholog panel (next section).
A [proposed draw update](ORTHOLOG_DRAW_POLICY_UPDATE.md) applies the same rule to all 3,643
sampleable forms before the draw (195 fail, mostly cyno diverged) and adds a start-site rule;
it is not adopted.

### The L1 ortholog panel

The L1 ortholog panel is the CDS-eligible draw (`ortholog_panel_cds_eligible.tsv`) with
each form the correspondence review recommends replacing swapped for its proposed
replacement, the first architecture-conserving form after the drawn ones in the same
species, band and human topology class. The six forms the review asks about are kept as
drawn: none is shown to lose a human feature (Becca, 2026-11-10;
`ortholog_panel_adoption_judgements.tsv`).

| Replaced form | Missing human feature | Replacement |
|---|---|---|
| cyno ADGRB2 | one TSP_1 repeat (55-residue deletion) | cyno MRGPRX3 |
| cyno MS4A10 | TM1 (40-residue deletion) | cyno CCR5 |
| cyno SLC8A1 | TM6–TM10 and `Na_Ca_ex` (C-terminal truncation) | cyno MRGPRX4 |
| cyno TMEM67 | TM2 (33-residue deletion) | cyno GPR78 |
| cyno FLT4 | second Ig domain (47-residue insertion) | cyno CSF2RB |
| cyno IL31RA | signal peptide, `IL6Ra-bind`, first `fn3` | cyno CD84 |
| cyno IL6R | signal peptide (internal after an N-terminal extension) | cyno AMIGO3 |
| cyno ITGA10 | signal peptide (N-terminal extension) | cyno TNFRSF10B |
| cyno NFASC | signal peptide, N-terminal orientation | cyno NCR1 |
| mouse ADGRE5 | one `EGF_CA` domain | mouse AVPR1B |
| mouse BTLA | the Ig domain | mouse IL15RA |

Kept as drawn: cyno MARCO and mouse OR10G6, TMIGD1, IGSF1, TMPRSS7 and CACHD1.

Every form is built from its exact pinned CDS: the archived final-v2 `selected_cds.json.gz`
for a drawn form, and the final-v2 eligibility ledger's chosen encoding for a replacement.
Each CDS translates to the panel's protein hash. `adopt_ortholog_cds_policy.py` reconciles
members on the designed protein, with ids never renumbered:

| Form's protein | Members | Action |
|---|---:|---|
| an active L1 ortholog's, identical CDS | 18 | keeps its id |
| an active L2 ortholog's, identical CDS | 126 | new L1 id; L2 row retired `moved_to:<id>` |
| no member's | 56 | new L1 id, built from the pinned CDS |
| (old L1 ortholog not in the panel) | 182 | retired `dropped_from_panel:ortholog_cds_policy` |

The 182 new ids are `l1-p2932`–`l1-p3113`, above the 189 ids the never-activated first
Compara 116 panel reserved (`compara116_id_allocations.tsv`). An identical protein whose
recorded CDS differs would stop the run; none does. Eight moved members carry a different
transcript id in their record for the same DNA. L2 keeps its other 2,451 active orthologs:
the policy governs the L1 draw, and 42 of their CDSs would not pass its automatic gate
(39 need review, 3 fail source integrity) while 20 are absent from the Compara 116 ledger.

The designed-CDS and construct records are append-only: the new records
(`library-candidates/2026-11-10-ortholog-cds-policy/`, named by `library_design.CDS_KEY`
and `CONSTRUCTS_KEY`) hold every earlier entry unchanged plus the 56 new forms; their
p0039 vectors are in `added_vectors.fasta` there. Topology comes from the
`2026-11-06-ortholog-cds-policy` snapshot (`topology.SNAPSHOT_KEY`), which keeps all
71,316 rows of the 2026-11-01 snapshot. Every active L1 member has a v1 row; no existing
L1 member's match changed, and 166 L2 orthologs gained one. The 25 new C-terminal members
got linker encodings in issue batch 3; the 28 retired members that hold an encoding keep
it, and batch 3 is drawn around them.

| After adoption | Value |
|---|---:|
| Active L1 members (reference 2,360, isoform 150, ortholog 100 mouse + 100 cyno; the 25 rescue members retired 2026-11-13, one design per gene) | 2,710 |
| … with a designed CDS | 2,710 |
| Active L2 members (ortholog 2,451, isoform 1,678, giant truncation 24) | 4,153 |
| C-terminal members | 215 |
| Internal-site members | 571 |
| Members over the 100-nt swap bar | 255 |
| New members with an identical AlphaFold model | 116 of 182 |
| New members whose AlphaFold model is a different protein | 5 (same length, 1–few residues) |
| New members without an AlphaFold DB model | 61, folded with ColabFold 1.6.2 (2026-11-10) |

The 61 new members without an AlphaFold DB model were folded with ColabFold 1.6.2 on Modal
(`scripts/fold_missing_structures_modal.py`; Becca approved the spend on 2026-11-10; 61/61
succeeded, 5.7 GPU-hours), each from its member's own designed protein, translated from its
CDS. The models are at
`processed://sprint4/s4e9-orf-design-fg-panel/predicted-structures/colabfold/2026-11-10-l1-orthologs/`
and their rows are in `predicted_structures/fold_summary.tsv`.

The 5 adopted orthologs whose AlphaFold DB model is a different protein (cyno SLC47A1
`l1-p2969`, MC4R `l1-p3001`, SLC22A8 `l1-p3013`; mouse ICAM1 `l1-p3036`, TM7SF3 `l1-p3059`)
were folded the same way from their designed proteins, named by design id because the
accession's model is the wrong protein (Becca approved, 2026-11-10; implemented by the Claude
Code session l1-reviewer-package, PR #283). 5/5 succeeded, mean pLDDT 75.7–86.6; models at
`processed://sprint4/s4e9-orf-design-fg-panel/predicted-structures/colabfold/2026-11-10-l1-ortholog-mismatched/`,
summary `predicted_structures/fold_summary_mismatched.tsv`.

The replacement gene-level pool uses the union of the human, mouse and cynomolgus
Compara 116 protein homology dumps, read in both orientations. Ensembl documents each
species-specific dump as an arbitrary subset, so one species file alone does not prove
complete coverage. One-to-one high-confidence relationships remain the eligibility
rule. Genes missing an unambiguous primary-assembly mapping remain explicit accounting
rows. The selected reference transcript disambiguates pseudoautosomal X/Y copies and
renamed symbols. A pre-existing library candidate is not required for admission.

| Rebuilt pool, before proteoform and topology checks | Mouse | Cyno |
|---|---:|---:|
| One-to-one high-confidence pair with available Ensembl translations | 1,876 | 1,862 |
| Of these, absent from the April ortholog lists | 431 | 426 |
| No one-to-one high-confidence pair | 484 | 498 |
| Library gene lacks an unambiguous primary-assembly mapping | 9 | 9 |

The proteoform comparison uses each human reference's translated designed CDS, rather
than recovering a sequence by length. It scores 3,737 gene/species pairs against the
Compara-linked Ensembl translations, pinned UniProt 2026_03 proteoforms, and the
recorded cohort form where its ortholog gene still agrees. The correspondence metric
and reviewed/APPRIS/current-pick tie preferences are retained; a final accession sort
makes new ties deterministic. UniProt candidate lookup within an ortholog gene retains
the existing symbol rule, whose replayability is a separate open task.

The completed checks include 323 new DeepTMHMM v1 (1.0.24) predictions and 242 NetGPI
1.1 predictions. Predicting winners, then rescoring, avoids running every unselected
alternative. The usable surface pool has 1,870 mouse and 1,836 cynomolgus forms. Every
excluded or unresolved pair remains in the accounting tables. Cynomolgus OR4Q2 has a
qualifying pair but cannot be scored because the human reference lacks a designed CDS;
27 selected forms are known non-surface, and four anchor predictions are explicitly
uncertain. Neither missing predictions nor an old manifest entry restricts the draw.

The replacement panel contains 100 mouse and 100 cynomolgus forms, with 50 near and
50 diverged forms per species. It retains 11 of the existing panel's proteins and
replaces 189; 45 picks come from genes absent from the April lists. Every pick has an
exact Ensembl translation or an existing designed CDS, and no two picks share a protein
or duplicate an active non-ortholog L1 member. Fifty forms have no existing library
member ID. The new panel and its change/rejection tables are review records; the
current manifest and ordered constructs are unchanged. The isoform panel is unchanged.

Both private and public D1 contain the 323 new topology rows, 2,359 NetGPI records
(242 new and 2,117 historical), and 2,187 separate curation/source/derived records. NetGPI
records preserve the raw header and row, tool version, one-based prediction number,
submitted ID, exact likelihood string, omega position/residue, source checksum and
job ID when recorded. The 2,117 legacy predictions without a verified original input
remain `unresolved_legacy`; no sequence is inferred from accession or length.
The prior UniProt GPI-keyword and curated architecture inputs, supplements,
exclusions, blind-spot notes, manual form choices and derived gene calls
remain independent of predictor output. The human GPI set and its four recorded
curation conflicts are preserved, not re-decided by this ortholog run.

Methods and provenance: `rebuild_ortholog_pool.py`, `score_ortholog_pool.py` and
`draw_rebuilt_ortholog_panel.py`, `export_ortholog_predictions.py`, and
`ortholog_rebuild_inputs.json`; [Compara 116 homology dump documentation](https://ftp.ensembl.org/pub/release-116/tsv/ensembl-compara/homologies/homo_sapiens/README.gene_trees.tsv_dumps.txt).
The immutable run archive is `processed://sprint4/s4e9-orf-design-fg-panel/library-design/ortholog-rebuild/2026-11-05/`,
with `prediction_runs.json` recording input/output paths, external run IDs and the
curation policy. The D1 importer and independent GPI tables are in
[surfaceome PR 463](https://github.com/Deliverome-Project/accessible-surfaceome/pull/463).

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
| Decision tree and evidence | `library-design/METHODS_library_design.md` | on S3 and in `library/l1/tag-site-architecture/` |
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

### C-terminal members carry their cassette: 211 at the C-terminus, ITM2B before its furin site

The C-terminal rung is the 212 L1 members whose N-terminus is cytoplasmic and whose
C-terminus is extracellular by DeepTMHMM v1 (`l1_tag_paths.tsv`). For these members the
cassette goes straight onto the end of the protein: the native stop is removed, then the
member's issued cassette encoding follows, then TAA. Nothing about the cassette is chosen
at this step. The peptide (the 51-aa SpyTag003/P-ALFA-P tandem), the tag codons (p003's),
the per-member linker encodings (`l1_cassette_encodings.tsv`) and the layout (p003's,
`build_constructs.py`) are all fixed earlier in the record. `build_l1_cterm_arms.py`
assembles them into one synthesis insert per member, inside the same Esp3I wrapper that
`assemble_candidates.py` puts around every L1 ORF:

    CGTCTCA | Kozak | ORF, native stop removed | cassette | TAA | CC | 5' handle | barcode | RT site | 10x CS1 | AGA GAGACG

**Member ids.** Rows are keyed by `candidate_id` (`l1-mNNNN`, plate order) and also carry `design_id` (`l1-pNNNN`), which never changes. The UniProt snapshot and the ColabFold model files are keyed by `design_id`, so a renumbering does not orphan them.

**Member set.** These counts are for the member set of PR #283 at 13f3c1b8 (second pass plus the cyno start-site restarts): the Compara 116 ortholog adoption, the gated isoform redraw, and two rounds of reference fixes and isoform rejections. That gives 212 C-terminal members: 163 library references, 19 isoforms and 30 orthologs. Since the first build of 218, #283's decisions have removed 33 members, including the four partial-form isoforms, and TSPAN10 moved to the internal rung after its reference fix. 27 joined, two of them ECE1 isoforms. Members present in every build keep byte-identical inserts. Both scripts read their members from `l1_tag_paths.tsv` and `candidate_manifest_tranches.tsv` and hold no list of their own. The rerun after a membership change is: issue cassette encodings for any new C-terminal member (`l1_swap_risk.py`), rerun the review with `--fetch` to a new dated UniProt key, then rebuild.

**The barcode is a 20-N placeholder.** Barcodes are not reserved yet, and the boundary
screens in `reserve_barcodes.py` run when they are. Everything upstream of the barcode,
including the whole coding sequence and the stretch that sets swap risk, is final.

**A C-terminal topology call does not by itself show that a C-terminal tag survives.**
`review_l1_cterm_members.py` checks each member before the build and writes the result to
`l1_cterm_review.tsv`, one row per member. Each check uses the member's own designed
protein: the sequence for v1 and v2 topology, the last 15 residues for UniProt features, and
the member's own structure model. Only a member whose C-terminal tag would be cut off the
protein is marked `cterm_tag = removed`. Such a member gets no C-terminal tag, and it is built
only at a decided internal site; the one case here (human ITM2B) has one. Every other risk is recorded
as a flag and the member is built at the C-terminus.

| Check | Rule | Members | Outcome |
|---|---|---:|---|
| C-terminal peptide removed | A UniProt Peptide/Propeptide feature covers the last residue | 1: ITM2B human (`l1-m0355`); the mouse ortholog left with #283 | **Not tagged at the C-terminus.** Furin cuts at 243/244 and releases the 23-aa Bri23 peptide, so a C-terminal tag is secreted with it, which is the same mechanism as a GPI omega site. The cassette goes after residue 235 instead (below) |
| GPI | Gene in `gpi_set.tsv` | 0 | `l1_tag_paths.py` already moved BST2 (and, until its reference changed on 2026-11-12, RGMB) to the internal rung |
| Ectodomain released | UniProt records shedding, furin processing or a soluble or secreted form | 24 (74 with `shed_risk` included) | Built. On a type II protein the C-terminus sits on the released domain, so the tag goes with it. Any extracellular tag on these members has the same exposure, and no other rung exists for them |
| Two-chain enzyme | UniProt records more than one chain, named heavy/light or (non-)catalytic | 14 (GGT1/5/6/7, HPN, TMPRSS2/4/11D/11E/15) | Built. An activation cleavage puts the tag on the C-terminal chain, which stays on the cell only while it stays associated with the anchored chain (by a disulfide in the serine proteases); on a gel the two can run apart |
| Terminus inside a fold | Last-residue pLDDT ≥ 85 in the member's own model | 37 (198 of 212 have a model; 14 new or rebuilt members have none on S3) | Built. Two families account for most of these. The TNF-family ligands (TNF, FASLG, CD40LG, CD70, TNFSF11/12/13B/14/15) end on the last strand of the TNF-homology domain. The M13 metalloproteases (MME, MMEL1, ECE1, PHEX, KEL) end on a conserved …C(R/E/Q)(V/L)W motif. Whether a 51-aa extension is tolerated there is not established here |
| Partial form | Member under half its gene's UniProt length | 0. #283's isoform gate rejected SLC43A2, TFRC, CLEC7A and SLC6A11 | — |
| Short extracellular tail | v1 C-terminal extracellular run under 5 residues | 9 (SLC29A1/2, SLC47A1 ×2, SLC47A2, SLC38A8, RHBDL2, STEAP1B, TRPV1) | Built, except the TRPV1 isoform (held; see "C-out calls" below). The cassette starts at the membrane, behind its own GGGGS linker |
| v2 disagrees | DeepTMHMM v2 on the same sequence puts the C-terminus in the cytoplasm, or calls no TM | 2 cytoplasmic (TRPV1 isoform `l1-p2528`, which also has a 1-residue v1 tail, and SLC68A1, which is C-terminal by a lane override); 3 no TM (BRICD5, TCTA, TNFSF13B); 3 with no v2 row (orthologs CLEC2L, CORIN, SLC47A1) | Built on v1, the topology of record, unless an independent source contradicts the call. The TRPV1 isoform is held (see "C-out calls" below) |

**ITM2B gets its cassette after residue 235, before the furin site** (Becca, 2026-11-10).
Furin processing of BRI2 is constitutive ([Kim et al., 1999](https://doi.org/10.1038/14783)).
UniProt places the cut between R243 and E244, and a tag placed directly after R243 would leave
the whole cut motif (…IQKR243) in front of it. Residue 235 meets four constraints:

| Constraint | At residue 235 |
|---|---|
| Past the BRICHOS domain | the domain is residues 137–231, and its 164–223 disulfide is untouched |
| Past the confidently predicted fold | pLDDT is 97 at 230 and 86 at 234 in the member's own model |
| In the linker | pLDDT is 76 at 236 and 62 at 240 |
| Leaves the furin site's context intact | the eight residues before the cut (236–243: TIKGIQKR) stay contiguous, and the build asserts them |

After furin, the tag is on the membrane form. ADAM10 later sheds the BRICHOS-containing
ectodomain ([Martin et al., 2007](https://doi.org/10.1074/jbc.M706661200)), so ITM2B is in the
same ectodomain-released class as TNF or FASLG.

**Unknown:** whether a 51-aa insert eight residues before the motif changes how efficiently
furin cuts. That would show as uncut, Bri23-bearing protein on a blot. It cannot affect
whether the tag stays on the cell, because the cassette is upstream of either cut.

Each insert is checked as it is built, and `library_pipeline.py verify --stage cterm_arms`
rebuilds all 232 inserts byte for byte. 229 end ORF, cassette, stop (201 at the true C-terminus,
including SLC68A1, and 28 snorkels). In 3 the cassette is internal and native sequence follows it:
HPN and KLRF1 at their stalk sites, and ITM2B before its furin site. The `cassette_ends_orf`
column of `constructs_l1_cterm.tsv` records which is which, and every table is joined on design_id:

| Check | Result |
|---|---|
| Protein | exactly the designed protein followed by the 51-aa cassette (for ITM2B, the cassette after residue 235); one stop, at the end |
| Tag DNA | SpyTag003 and P-ALFA-P on p003's codons in every insert |
| Type IIS | 2 Esp3I (the wrapper) and no BbsI or BsaI in every insert; in the p0039 backbone, the backbone's 4 sites plus those 2 |
| Junction | none in the built arms. The rule (`l1_swap_risk.domesticate_junction`) recodes an ORF's last codon synonymously when it and the cassette's first Gly codon complete a Type IIS site; its one case, SLC41A2 (GAC→GAT), is held (see "C-out calls" below) |
| Swap risk, on the built DNA of all 2,710 active L1 members | all 232 arms: 36 nt median and 36 nt maximum; none over 80 nt |
| Vendor complexity features touching the cassette or its junctions | 46 inserts (12 C-terminal, 34 snorkel), each with one 50-nt window at 76–80% GC that spans the cassette or the ORF–cassette junction (band 25–75%); no repeats, homopolymers or hairpins |

Inserts are 596–4,778 nt (`constructs_l1_cterm.fasta`), and per-member build facts are in
`constructs_l1_cterm.tsv`. `candidate_manifest_tranches.tsv` still records
`tag_cassette = deferred` for every member. That manifest is the stage-10 record of the
untagged vectors, and this table supersedes it for the C-terminal rung.

**Calls made on the evidence** (Becca delegated protein-level calls on 2026-11-12: "use your
best judgment"). One row per member and decision is in `l1_cterm_decisions.tsv`, written by
`l1_cterm_decisions.py` from the review table. Each row names the rejected alternative, so
it can be reversed.

| Question | Call | Main evidence | What would change it |
|---|---|---|---|
| Furin-processed members (COL23A1 ×2, COL25A1 ×2, EDA ×3, TNFSF12, TNFSF13B): mutate the site or move the tag? | Neither. The cassette stays at the native C-terminus | The cut releases the whole functional ectodomain and leaves a short stalk stub on the membrane. A tag before the site would report the stub, and a mutated site would make a non-native, surface-retained protein. At the C-terminus the tag reports the native, uncut surface fraction. ITM2B differs because furin removes only a 23-aa peptide there and the folded domain stays | Evidence that the surface fraction is undetectable |
| The 36 C-terminally tagged members whose terminus ends inside a fold | The cassette stays at the C-terminus until each refold is read | The N-terminus is cytoplasmic, so the only alternative is a stalk site (`stalk_site_after` in the review). The most buried termini (MME, ECE1, KEL; 16–20 neighboring Cα within 10 Å) have no stalk of 8+ residues | A refold that disrupts the native fold. The member then moves to its stalk site, or is held where it has none |
| Stop context TAAC (pad `CC`) or TAAG (`GC`) | Keep TAAC | TAA is the strongest stop codon either way. CC is the pad of every construct designed and ordered since 2026-10-01, and the pad is shared by every insert on every rung | A readthrough measurement on a tagged construct |

**Still open on this rung:**

| Question | What would resolve it |
|---|---|
| ITM2B: does the cassette at residue 235 change furin processing? | A blot for uncut, Bri23-bearing ITM2B in the tagged line |
| Refold verdicts for the tagged proteins (sent to the refold session on 2026-11-11 and 2026-11-12) | Its verdicts. A failed fold triggers the stalk move or hold above, applied without further sign-off |
| The 46 inserts with a 76–80% GC window | The vendor scan (issue #13); if it flags them, the linker encodings would need re-issuing |

**Tagged refolds, and what changed because of them.** l1-tag-refold folded each C-terminal-rung
member with its cassette and compared it with the member's own native model (snapshot of its
verdicts: `l1_cterm_refold_verdicts.tsv`, joined on design_id, with the SHA-256 of each tagged
protein checked). Of the 212, 200 pass and 12 fail. GGT5, ECE1 (`l1-p0537`) and CEMIP2 had no verdict in the first
snapshot; their rows are refreshed from #301's final `tag_refold_verdicts.tsv` (same tagged-protein
SHA-256), and all three pass. Each structured-terminus row in `l1_cterm_decisions.tsv` now states
what its refold decided. One rule decides every failing member, and it is a function
(`refold_disposition` in `build_l1_cterm_arms.py`), not a list of members. The stalk site is the
midpoint of the stretch between the last TM helix and the folded domain, where that stretch is at
least 8 residues (`stalk_site_after` in the review).

| The C-terminally tagged protein | Stalk site | Outcome |
|---|---|---|
| Passes its refold | any | Cassette at the C-terminus |
| Fails (fold disturbed, epitope buried or packed on the native domain past the verdict's limit of 5 residues, or a helix inside the cassette) | The stalk-site design passes its own refold | Cassette moves to the stalk site |
| Fails | The stalk-site design fails too | Held |
| Fails on the fold (native-region lDDT under 0.85, or lDDT within 15 residues of the site under 0.7) | None usable | Held |
| Fails with the fold intact (epitope packing or a cassette helix alone) | None usable | Cassette stays at the C-terminus, flagged |
| Fails, and the ectodomain is natively released | any | Cassette stays at the C-terminus: a stalk site lies before the cut and would tag the membrane stub |

Where a stalk site exists, the trigger is the verdict, whatever its reason. KLRF1 moved on
epitope packing alone (fold intact, native-region lDDT 0.99) and stays at its stalk site, because
that arm passes its refold with no epitope residue packed. TNFSF18 has the same packing and stays
at the C-terminus because its stalk is 6 residues, under the 8 a site needs. The PLVAP isoform is
held: its failure is in the fold at the site (lDDT 0.39), and it has no native model from which a
stalk site could be designed. No member departs from the rule.

| Outcome | Members | Evidence |
|---|---|---|
| Moved to the stalk site, which passes | HPN (43), KLRF1 (84) | At the C-terminus: site lDDT 0.47 with 5 epitope residues packed (HPN); 6 packed (KLRF1). At the stalk site: native-region lDDT 0.915 and 0.992, no epitope packed |
| Held: the stalk site fails too | CD70 (46), CD209 (67), SCARA5 (86), BRICD5 (63) | Native-region lDDT at the C-terminus / at the stalk site: CD70 0.79 / 0.81, CD209 0.77 / 0.82, SCARA5 0.82 / 0.83, BRICD5 0.78 / 0.77 (site lDDT 0.32 / 0.33). Each is the gene's only L1 member, so the gene leaves L1 |
| Held: the fold fails and there is no usable stalk site | PLVAP isoform `l1-p2487` | Site lDDT 0.39; no native model, so no stalk site. The PLVAP reference and ortholog pass |
| Kept, flagged | EDA ×2, TNFSF18, MCEMP1, ITM2B | Fold intact (native-region lDDT ≥ 0.94); 3–6 epitope residues pack on the domain (MCEMP1: an 11-residue helix forms inside the cassette linker instead). EDA's stalk site lies before its furin cut, so it would tag the stub. TNFSF18 and MCEMP1 have no stalk. ITM2B keeps Becca's site |

l1-tag-refold then folded the six stalk-site designs (PR #301, 8af3b01d; snapshot
`l1_cterm_stalk_refold_verdicts.tsv`, joined on design_id with each built protein's SHA-256
checked). The rule recorded before either verdict existed holds a member whose stalk site also
fails, so CD70, CD209, SCARA5 and BRICD5 are held. The rejected alternative is to build the
stalk-site arm with a refold flag. That would keep the four genes in L1, and a weak arm would
show up as low display in the screen. Three of the four are oligomers in their native state (the CD70 TNF trimer, the
CD209 coiled-coil neck, the SCARA5 collagenous trimer), and BRICD5's reference is a ColabFold
model of a mismatched sequence. Each one's monomer prediction drops by a similar amount at both
sites. So the failures may reflect the reference as much as the tag. A passing site from the
ranker would reverse the hold. #301's verdicts for the other C-terminal-rung failures
(EDA ×2, TNFSF18, MCEMP1, ITM2B) are for the same tagged proteins with the same metrics, so
their calls stand. Its 22 snorkel failures are still under test by the ranker, and every snorkel
arm stays built until it decides. One row per decision, with the rejected alternative, is in
`l1_cterm_decisions.tsv`.

**C-out calls.** DeepTMHMM v1 is the topology of record, so a C-out call is kept unless an
independent source contradicts it. The review records three things from UniProt for each
member's own sequence: the side of its C-terminus with the evidence code (`uniprot_cterm_topology`;
ECO:0000269 is experimental, ECO:0000255 a sequence-model prediction), a cleaved signal peptide
(`uniprot_signal_peptide`) and the number of TM helices (`uniprot_tm_helices`). It also marks a
call weak when v1's extracellular tail is under 5 residues, or when DeepTMHMM v2 puts the
C-terminus in the cytoplasm, calls no TM helix or has no row. One rule (`cout_disposition` in
`build_l1_cterm_arms.py`) then decides:

| Evidence against C-out | Outcome |
|---|---|
| UniProt annotates the C-terminus cytoplasmic on experimental evidence | Held |
| UniProt describes a type I protein: a cleaved signal peptide, one TM helix, a cytoplasmic C-terminus. v1's C-out call then rests on reading the signal peptide as the only TM helix | Held |
| The call is weak, and UniProt annotates the C-terminus cytoplasmic on any evidence | Held |
| v2 puts the C-terminus in the cytoplasm, and neither a UniProt extracellular annotation at that C-terminus nor a lane override places it outside | Held |
| An ortholog with no UniProt entry whose human member is held under the first two rows | Held |
| UniProt's cytoplasmic annotation is a prediction, v1 and v2 agree, and the protein is not type I | Kept, flagged |
| The call is weak and nothing contradicts it | Kept |

Four members are held and five are kept with a flag. Each of these, and each of the 15 weak
calls, has one row in `l1_cterm_decisions.tsv`.

| Outcome | Members | Evidence |
|---|---|---|
| Held: cytoplasmic by experiment | SLC41A2 `l1-p1957` | The member is identical to UniProt Q96JW4, whose topological domains run from extracellular 1–162 to cytoplasmic 565–573, the last feature, all citing one study. Its abstract reports that heterologously expressed SLC41A2 has an N-terminus-outside, C-terminus-inside, 11-TM topology ([Sahni et al., 2007](https://pubmed.ncbi.nlm.nih.gov/16984228/)). v1 and v2 both call 11 helices in the reverse orientation. It is a candidate for the N at residue 1 rung |
| Held: type I by UniProt | VTCN1 `l1-p2342` and its cyno ortholog `l1-p2628` | UniProt Q7Z7D3 (identical to the human member): signal peptide 1–24, TM 260–280, cytoplasmic 281–282, all predictions. v1 and v2 call residues 7–29 the only TM helix and no helix at 260–280. The cyno ortholog has no UniProt entry, is 98.6% identical and has the same C-terminal segment. VTCN1 is not in `gpi_set.tsv`. Both are candidates for the N after SP rung |
| Held: v2 cytoplasmic, nothing places it outside | TRPV1 isoform `l1-p2528` | The 521-aa isoform is canonical TRPV1 (UniProt Q8NER1) residues 1–515 plus six residues of its own. It ends five residues into canonical helix S3 (511–532), which follows the cytoplasmic S2–S3 loop (498–510), and it lacks S4–S6 and the pore. v1's C-out call rests on a 10-residue third helix and one extracellular residue; v2 calls two helices and a cytoplasmic C-terminus. Canonical TRPV1's own C-terminus (684–839) is cytoplasmic. It is a candidate for the internal rung |
| Kept, flagged for Becca: predicted cytoplasmic, v1 and v2 agree C-out | LYSMD3, LYSMD4, SLC24A1, TMEM222 | LYSMD3 and LYSMD4: one TM helix (218–238), no signal peptide in UniProt, so UniProt's N-out and DeepTMHMM's N-in readings differ only in orientation, and neither is experimental. SLC24A1: v1 calls 13 helices, UniProt predicts 11, and UniProt's signal peptide is annotated as not cleaved. TMEM222: v1's call rests on an 11-residue helix (191–201) and a 7-residue tail; UniProt predicts three helices and a single cytoplasmic residue at 208. None is in `gpi_set.tsv` |
| Kept, weak: UniProt annotates the member's C-terminus extracellular | SLC29A1 and SLC68A1 (experimental); SLC47A1, SLC47A2, TCTA, TNFSF13B, CORIN mouse (prediction) | `uniprot_cterm_topology` in the review. SLC68A1 is C-terminal by a lane override |
| Kept, weak: v2 agrees with v1, and UniProt has no topology feature at the C-terminus | RHBDL2, SLC29A2, SLC38A8, STEAP1B | v2 C-terminus non-cytoplasmic |
| Kept, weak, on v1 alone | SLC47A1 cyno, CLEC2L mouse | No v2 row and no UniProt topology feature; nothing contradicts v1 |
| C-out call stands; held for its refolds | BRICD5 | v2 calls no TM helix |

**Unknown:** for SLC41A2 only the abstract was read, not the experiments behind the topology. For
VTCN1, LYSMD3, LYSMD4, SLC24A1 and TMEM222 no experimental topology is on record here, and an
experiment on any of them would settle its row.

**Snorkel arms for 28 internal-rung members.** These members have both termini cytoplasmic. The
ranker's list (`l1_snorkel_designs/2026-11-14b/l1_snorkel_designs.tsv`, PR #331, pinned
by SHA-256 in `l1_snorkel.py`) assigns them a snorkel: an extra TM helix after the native
C-terminus, so the cassette sits outside. Most are members whose every designed internal loop
site failed its refold (MS4A, TM4SF, KCNK, connexin, CD9 and others). NRG2 is a snorkel because
its loop site lies on the shed ectodomain. P2RX5's internal sites all lie within 2 residues of the
ATP-pocket ligand contact. P2RX2's one contact-free internal site fails its refold. GJA8, OTOP1 and SLC23A2 (added
2026-11-14) have every internal site vetoed by a verified literature ligand contact. The snorkel
arms of these four have not been refolded yet. The list in `constructs_l1_cterm.tsv` (placement
`snorkel`) is read from the pinned table, never typed here.

The first list had 45 members. When the snorkel arms themselves were refolded (PR #301), 21 of
them failed, and so did MS4A7's new arm (`l1-p3040`). A 30-residue spacer rescued none of the 5
the ranker tested. Those 22 members are held with no arm: CLRN3 (`l1-p0434`), EMP1–3, MS4A1, 2,
6A, 7 (both members), 10 and 13, PANX1, PLLP, PMP22, SLC1A4, SLC35G2, TM4SF1, TM4SF5, TMEM169,
TMEM174, TMEM204 and TSPAN12. SLC22A24 `l1-p1884` left on 2026-11-14, when its proteoform-fix
redesign gave it an internal site after residue 59; the other SLC22A24 member (`l1-p3149`) stays a
snorkel. The batch-6 encodings of all 22 stay issued and reserved
(`l1_snorkel.reserved_designs`): later batches are drawn around them, and batch 6 still replays
on their snorkel flanks.

| Step | What is done |
|---|---|
| Design | Read from `processed://sprint5/s5e13-tag-site-benchmarking/l1_snorkel_designs.tsv`, pinned by SHA-256 in `l1_snorkel.py`. The built protein must equal its `snorkel_protein` exactly |
| Module DNA | The 55 residues before the cassette are encoded per arm (`diversify_snorkel_modules.py`, `l1_snorkel_module_encodings.tsv`; Becca, 2026-11-21). The single encoding they replaced (`l1_snorkel.PREFIX_DNA`: DNA Chisel, human codon usage, no BbsI/Esp3I/BsaI, seed 20261112) stays pinned as the flank the cassette batches were issued against |
| Cassette | Each member gets its own linkers-only encoding from `l1_swap_risk.py`, which treats snorkel members as cassette-before-barcode: batch 6 for the first list, batch 7 for P2RX5 and P2RX2, batch 8 for SLC68A1 (C-terminal), batch 9 for GJA8, OTOP1 and SLC23A2. All 296 encodings replay from scratch in 9 batches, and none equals any encoding issued before, the retired rescue arms' included |
| Checks | All 28 pass the Type IIS check, alone and in p0039. The worst-partner tract is at most 36 nt. Each insert has one 76% GC window at the cassette's end, the same class as the 12 C-terminal ones above |

**Decisions on this rung.** Attribution follows the convention of 2026-11-10: every manual fix
names who decided it and who implemented it. "This session" is the Claude Code session on
branch `claude/l1-cterm-build-16eeeb`, which opened PR #284.

| Date | Decision | Decided by / implemented by | Why |
|---|---|---|---|
| 2026-11-10 | C-terminal-rung members carry the cassette at their true C-terminus | Becca (the brief for this build) / this session | The cassette and its encodings were already fixed, so only the junction remained |
| 2026-11-10 | Risks other than tag removal are flagged, not rerouted | Existing rule in TAG_PATHS ("Shedding is recorded and does not change the rung") / this session | No other rung exists for these members |
| 2026-11-10 | ITM2B (`l1-m0355`; the mouse ortholog has since left with #283): cassette after residue 235, before the furin site | Becca, choosing this over mutating the furin site / this session, which proposed residue 235 | A C-terminal tag leaves with Bri23. Directly after R243, the cut motif would sit in front of the tag |
| 2026-10-29 / 2026-11-10 | SLC41A2 (`l1-m1697`): last codon GAC→GAT | Rule `l1_swap_risk.domesticate_junction` (commit e545db2e, 2026-10-29) / reapplied to the domesticated ORF by this session | GAC followed by the cassette's first Gly codon completes an Esp3I site; the change is synonymous |
| 2026-11-10 | The four partial-form isoforms stay built and flagged | The L1 main session ("L1 library: stop recoding Spy/ALFA, fix dup bugs"), with the reject-or-length-floor choice put to Becca / no change made | The isoform panel belongs to `redraw_panels.py`, on the reviewer branch |
| 2026-11-12 | Furin-processed members keep the C-terminal tag (9 members); members ending inside a fold keep it pending the refold (36); stop context stays TAAC | l1-cterm-build (delegated by Becca 2026-11-12, "use your best judgment") / this session | See "Calls made on the evidence"; one row per member in `l1_cterm_decisions.tsv`, each with its rejected alternative |
| 2026-11-12 | Refold failures: 6 members moved to their stalk site, the PLVAP isoform held, 5 kept with a flag | l1-cterm-build (delegated by Becca 2026-11-12, "use your best judgment") / this session | The rule recorded before the verdicts; see "Tagged refolds" and `l1_cterm_decisions.tsv` |
| 2026-11-13 | Stalk-site refolds (PR #301): HPN and KLRF1 stay at their stalk sites; CD70, CD209, SCARA5 and BRICD5 held, so these four genes leave L1 | l1-cterm-build (delegated by Becca 2026-11-12, "use your best judgment") / l1-cterm-build on `claude/l1-cterm-stalk-refold-holds` | The same rule, recorded before the verdicts. Rejected: building the stalk-site arm with a flag |
| 2026-11-13 | Snorkel list moved to the ranker's 2026-11-13b list: 22 snorkel arms dropped (held with no arm, encodings reserved), P2RX5 and P2RX2 added as cassette batch 7 | ranker-pu-segment under the coordinator's hold rule (delegated by Becca 2026-11-12), PR #245 bf831be8 / l1-cterm-build on `claude/l1-cterm-stalk-refold-holds` | Each dropped snorkel arm fails its own refold, and the spacer test rescued none. An issued encoding is never redrawn or handed on |
| 2026-11-14 | Snorkel list moved to the ranker's 2026-11-14 list (29): GJA8, OTOP1 and SLC23A2 added as cassette batch 9; no removals | ranker-pu-segment (delegated by Becca 2026-11-12), PR #245 322a2bc0 / l1-cterm-build on `claude/l1-cterm-stalk-refold-holds` | Every internal site of the three is vetoed by a verified literature ligand contact |
| 2026-11-14 | Snorkel list moved to the ranker's 2026-11-14b list (28): SLC22A24 `l1-p1884` removed, its encoding reserved | ranker-pu-segment (delegated by Becca 2026-11-12), PR #331 / l1-cterm-build on `claude/l1-cterm-stalk-refold-holds` | The proteoform-fix redesign gives it an internal site after residue 59 |
| 2026-11-13 | SLC68A1 (`l1-p1995`) built as a C-terminal arm (cassette batch 8); the review takes its topology from the lane override, and its UniProt snapshot `2026-11-13-slc68a1` adds only Q14CX5 to the second-pass entries | l1-reviewer-package (lane override, delegated by Becca 2026-11-12, on ranker-pu-segment's proposal) / l1-cterm-build on `claude/l1-cterm-stalk-refold-holds` | Experimental topology (Q14CX5, ECO:0000269) puts both termini outside, the reverse of v1. Its terminus is a free tail (pLDDT 33) |
| 2026-11-14 | The C-terminal build joins its inputs on design_id, and the encoding ledger carries design_id | l1-cterm-build, on independent-audit finding 8 / l1-cterm-build | Member ids are regenerated by each renumbering; a stale one now fails the build instead of mis-joining |
| 2026-11-14 | One rule decides stalk moves and holds (`refold_disposition`), replacing the per-member lists; KLRF1 (`l1-p1070`) stays at its stalk site | coordinator under Becca's 2026-11-12 delegation and her 2026-11-14 instruction to manage the C-terminal rung / coordinator subagent, branch `claude/l1-cterm-audit-followups-2` | KLRF1's stalk arm passes its refold and its C-terminal form packs 6 epitope residues on the domain. Rejected: returning it to the C-terminus. An independent audit found the rule stated three ways. Would change it: evidence that the C-terminal epitopes are accessible |
| 2026-11-14 | The PLVAP isoform (`l1-p2487`) stays held, by a clause of the rule: a failed fold with no usable stalk site is held, and packing alone with no stalk site keeps the C-terminal tag, flagged | coordinator under Becca's 2026-11-12 delegation and her 2026-11-14 instruction to manage the C-terminal rung / coordinator subagent, branch `claude/l1-cterm-audit-followups-2` | A tag on a disturbed fold is not a usable arm, and a packed epitope on an intact fold may still be detected. Rejected: recording PLVAP as an exception to a rule that keeps every member with no stalk site. Would change it: a native model and a passing site |
| 2026-11-14 | TRPV1 isoform (`l1-p2528`) held: its arm is removed. Every other weak C-out call is kept, each with a decision row | coordinator under Becca's 2026-11-12 delegation and her 2026-11-14 instruction to manage the C-terminal rung / coordinator subagent, branch `claude/l1-cterm-audit-followups-2` | C-terminal orientation unconfirmed; canonical TRPV1 C-terminus is cytoplasmic; candidate for the internal rung. Rejected: building on the v1 call alone. Would change it: experimental evidence that this form's C-terminus is outside, or a passing internal site |
| 2026-11-14 | SLC41A2 (`l1-p1957`) held: its arm is removed | coordinator under Becca's 2026-11-12 delegation and her 2026-11-14 instruction to manage the C-terminal rung / coordinator subagent, branch `claude/l1-cterm-audit-followups-2` | UniProt's experimental topology puts the C-terminus in the cytoplasm, and the cited study's abstract says so ([Sahni et al., 2007](https://pubmed.ncbi.nlm.nih.gov/16984228/)). Rejected: building on the agreement of v1 and v2. Would change it: experimental evidence that the C-terminus is outside; a site on the N at residue 1 rung would move it there |
| 2026-11-14 | VTCN1 (`l1-p2342`) and its cyno ortholog (`l1-p2628`) held: their arms are removed | coordinator under Becca's 2026-11-12 delegation and her 2026-11-14 instruction to manage the C-terminal rung / coordinator subagent, branch `claude/l1-cterm-audit-followups-2` | The cyno ortholog has no UniProt entry and is held with its human member. UniProt describes a type I protein with a cleaved signal peptide and a two-residue cytoplasmic tail, so a C-terminal cassette would sit in the cytoplasm. Rejected: building on v1 and v2, which read the signal peptide as the only TM helix. Would change it: experimental topology; a site on the N after SP rung would move them there |
| 2026-11-14 | LYSMD3, LYSMD4, SLC24A1 and TMEM222 kept at the C-terminus, flagged for Becca | coordinator under Becca's 2026-11-12 delegation and her 2026-11-14 instruction to manage the C-terminal rung / coordinator subagent, branch `claude/l1-cterm-audit-followups-2` | v1 and v2 agree on C-out, UniProt's cytoplasmic annotation is a prediction, and none is type I by UniProt. Rejected: holding on a predicted annotation. Would change it: experimental topology; Becca's review |
| 2026-11-14 | EDA cyno (`l1-p2986`) left at the C-terminus | coordinator under Becca's 2026-11-12 delegation and her 2026-11-14 instruction to manage the C-terminal rung / coordinator subagent, branch `claude/l1-cterm-audit-followups-2` | The release call is transferred from the human gene, as the ortholog has no UniProt entry; its stalk site (126) was never refolded. Rejected: moving it to an unrefolded stalk site. Would change it: its own processing record, or a refold of the stalk site |
| 2026-11-14 | The review reads UniProt from the append-only copy `external://uniprot/rest-2026-11-13/cterm-review-slc68a1/uniprot_entries.json.gz`, pinned by SHA-256, and AlphaFold models only from the archive `external://alphafold-db/v6/`; an ortholog with no UniProt entry is flagged `uniprot_not_checked`; the released-ectodomain genes are derived from the review, not listed | l1-cterm-build, on the independent audit / l1-cterm-build on `claude/l1-cterm-audit-followups`, finished by the coordinator subagent on `claude/l1-cterm-audit-followups-2` | The external copy is byte-identical to the processed snapshot it replaces (same SHA-256). A model missing from the archive is recorded as "native model not archived" and the review continues; all 182 AlphaFold models this rung reads are archived |
| 2026-11-12 | Shared S3 keys `l1_cassette_encodings.tsv`, `l1_swap_risk.tsv` and `l1_swap_risk_summary.json` overwritten with #284's versions (#283's rows at 13f3c1b8 plus the 45 snorkel encodings, batch 6). Prior VersionIds: h5Yw20v3XM3yHFQv_rwLe4IJfTgBpFvo, LIVyOwv5z.TC0CosWaIAl8gNC3_J9xji, Nv5pne83aTM4ycfMhn73vqqzVN9tk4Lb (all from 13f3c1b8) | Becca, 2026-11-12 ("Yes, overwrite and push", in the l1-cterm-build session) / this session | The push hook requires the experiment-root copies; the bucket is versioned, so the prior versions can be restored |
| 2026-11-12 | Snorkel arms built for the 45 members of the ranker's final table at bb9cd99d (cassette batch 6, drawn once on the final list) | ranker-pu-segment (delegated by Becca 2026-11-12), PR #245, for the design / this session for the DNA, encodings and checks | Both termini are cytoplasmic; the snorkel's extra TM helix puts the cassette outside |
| 2026-11-12 | Member set rebuilt on PR #283 at 9928f4b1: 30 members out, 26 in | Becca, through #283's ortholog adoption, isoform gate and redraw / this session rebuilt from the tables, with no member chosen by hand | The arms follow the membership; this rung decides nothing about it |
| 2026-11-10 | Flag thresholds: pLDDT ≥ 85, tail < 5 aa, under half the reference length | No decision taken / this session | The thresholds set flags only, and nothing built depends on them |
| 2026-11-10 | `candidate_manifest_tranches.tsv` keeps `tag_cassette = deferred` | No decision taken / this session | The stage-10 verifier compares that column; `constructs_l1_cterm.tsv` carries the built inserts |

In `constructs_l1_cterm.tsv`, the `decided_by`, `implemented_by` and `reason` columns are filled
on the rows an adjustment touched: the ITM2B row (SLC41A2's junction edit left with its arm,
now held). They are
blank on rows built unchanged.

### Domestication and barcodes, audited on the final inserts

`audit_domestication_barcodes.py` recomputes, from the inserts themselves, whether domestication
and the barcode rules hold for everything L1 will order. It checks two sets: the 3,255 tagged
inserts (the five rung builds at 8432c010, from which the cut `2026-11-21-snorkel-modules` was made),
and the Twist-compliant recoded set made from that same build (`synthesis_compliance.py`; key
`synthesis/2026-11-21-snorkel-modules/`, sha 66332be4). Inserts are keyed by their permanent arm id.
The recoded set is the order of record, and the audit pins it by SHA-256. Its inserts are named by
current member ids, which renumber, so the audit matches each recoded insert to its build by the
permanent design id.

The 3,255 inserts are every built arm of the 2,673 built members; 37 of the 2,710 active members
are held with no insert ordered (13 on the N-after-SP rung, 20 internal, 4 C-terminal; brief
section 12).
The 25 expression-rescue
members (`l1-p2361`–`l1-p2385`) are retired, one arm per gene (Becca, 2026-11-13), so
`constructs_l1_rescue` is not read. Per-insert results are
in `l1_insert_audit.tsv`, every domestication edit is in `l1_domestication_edits.tsv`, and the
summary is in `l1_domestication_barcode_audit.json`. The tests in
`tests/test_domestication_barcode_audit.py` keep each result true.

**Domestication.** Every check holds on all 3,255 inserts, before and after recoding:

| Check | Result |
|---|---|
| Type IIS sites | exactly the two wrapper Esp3I sites; no BbsI, BsaI or other Esp3I site on either strand |
| Kozak and start | the designed Esp3I + p003 Kozak block + ATG in every insert |
| Protein | the coding sequence translates to the protein its arm table records, with one stop at the end |
| Edits | 8,201 edited codons in 2,389 members. Every one is synonymous, and each member's recorded `edits` count equals its logged nucleotide changes. The 25 retired rescue members are skipped |
| Unlogged changes | none. With the cassette cut out, every insert's coding sequence is the domesticated ORF. The Igκ-swap arms are the leader followed by the domesticated ORF from the first kept residue. SLC41A2, whose last codon (GAC→GAT) was the one logged exception, is held, so no ordered insert carries one |
| Barcode slot | exactly one 20-N slot per insert |

**Barcodes.** No member has a barcode yet. The audit re-derives the 3,200-barcode block reserved
for L1:

| Check | Result | Stated rule |
|---|---|---|
| Uniqueness | 3,200 distinct, and each appears exactly once in the shared registry | — |
| Composition | length 20, GC 40–60%, longest homopolymer 2, no ATG | the s1e2 pool's filters |
| Distance, full length | minimum Hamming 5, minimum edit distance 4 | pool: Levenshtein ≥ 3 and Hamming ≥ 4 |
| Distance, in situ prefix | minimum Hamming 3 over the first 10 nt | brief and code: 3 over 10 |
| Screens, both 3′ flanks | none fail (Type IIS, splice donor, junction miRNA seed, RBP motifs) | brief: the same screens, and the code's motif list matches the brief's |
| Inserts and backbone | no barcode occurs, on either strand, in any insert or in p0039 | not previously checked |
| Determinism | the draw takes the pool in index order, with no random number generator; replayed against the registry as it stood before the block, it returns the same 3,200 in the same order | — |

**Barcode reservation extended (2026-11-15).** With 2 or 3 arms per internal member, L1 grew to
about 3,265 inserts, more than the 3,200-barcode block. Becca approved extending the reservation
on 2026-11-15 ("Yes, extend it"). The extension is a second block appended to the registry, not a
regeneration: every pool barcode after the block's last index (20,892) that passes the same
screens under both 3′ flanks, with the in situ prefix screen run against the whole registry,
the 3,200 block included. That gives 444 barcodes (pool indices 20,910–29,909), named
`L1-res-3201`–`L1-res-3644` and recorded in `l1_barcode_block_extension_2026-11-15.tsv`
(`reserve_l1_block.py --extension`), for 3,644 reserved in all. They are every barcode left in
the pool that passes the screens (5,603 more fail only the prefix screen), so a library above
about 3,640 inserts would need a relaxed prefix rule or a second pool. The audit and the
synthesis check treat the block and the extension as one set, and the extension replays against
the registry rows that precede it. Decided by Becca Carlson; drawn and appended by the Claude
Code session l1-reviewer-package.

**Swap tract on the final DNA.** Across all 3,255 inserts, up to the barcode slot, against the
nearest insert of a different member: median 30 nt, maximum 3,852 nt, 288 inserts over 100 nt
(each insert's partner is in `l1_insert_audit.tsv`). Two inserts of the same member are measured
apart (`swap_tract_same_member`): 791 pairs on 373 members, shared tract median 366 nt, modelled
swap median 16.1% (brief section 1).

**Recoding.** The 4,502 codon changes across 727 inserts are all synonymous, and each insert's
count equals `synthesis_recodes.tsv`. No protein, cassette, readout block, barcode slot or wrapper
changes, and no BbsI, Esp3I or BsaI site or consensus splice donor (`GT[AG]AGT`) is added.

**Calls made on these results** (Becca delegated design calls on 2026-11-12: "use your best
judgment"):

| Finding | Call | Rejected alternative | Decided by / implemented by |
|---|---|---|---|
| 1,213 inserts carry an AarI/PaqCI site (`CACCTGC`), and BtgZI is in 1,092. Recodes add one of these sites to 9 inserts | No change. No cloning step uses AarI, PaqCI or BtgZI. The README names AarI/PaqCI as the likely enzyme for guide cloning; if it is adopted for a vector that carries these inserts, the ORFs must be domesticated against it then | Domesticate them all now, which would change those ORFs after their refolds and the Twist recode | l1-cterm-build (delegated by Becca 2026-11-12) / Claude Code session l1-cterm-build, this PR |
| Every insert carries a SapI site (`GAAGAGC`) inside p003's ALFA codons | No change. SapI is unused, and the tag codons are fixed to p003's | Recode the ALFA tag | same |
| Domestication edits create 471 new ATG trinucleotides (any frame), 51 `GT[AG]AG` and 2 poly(A) signals (KLRB1 codon 91, SLC46A1 codon 233) | No change. No stage claims these checks for ORFs. 808 of the 2,691 inserts already carry a native AATAAA/ATTAAA, and no recode adds the 7-mer donor the recode stage screens | Re-pick those codons, cascading into the constructs record and the recode | same |
| Each rescue arm shared 1.1–2.9 kb up to the barcode with its own library reference, a predicted 34–35% barcode swap within exactly the pairs that measured the Igκ leader | Superseded 2026-11-13: the 25 rescue members are retired, one arm per gene (Becca: "make a single choice for each L1"; PR #312), so no gene has two members and the finding no longer applies. Their 25 linker-only encodings are withdrawn | Recode each rescue arm's 3′ codons (PR #312's first fix), or pool each pair apart | finding: l1-cterm-build (delegated by Becca 2026-11-12); retirement: Becca, implemented by l1-nterm-res1-multipass (PR #312) |
| The first 2026-11-13 order of record (sha d1257f0a) held the 25 retired rescue inserts and 25 superseded multipass designs; the next (b426e433) was built from the draft integration branch, not `main` | Resolved: the order is regenerated from the rung builds (now `synthesis/2026-11-15-final-e/`, sha 20e8f275, 3,255 inserts, cut at 541c6c18) and re-pinned here; every insert matches its build | Filter the old order by hand | finding: l1-cterm-build / regeneration: l1-reviewer-package (PR #313) |
| Barcodes are screened under both candidate 3′ flanks | Unchanged. The choice of 3′ flank (today's RT site or the Feldman T9A padlock draft) remains Becca's | — | Becca (open) |

### Members whose protein duplicates another's are retired

A member's protein is its designed CDS translated. Comparing those proteins across the
library found 53 groups (108 members) in which two or more members encode the identical
protein, and each pair had a different gene symbol, accession or transcript.
`dedupe_library_members.py` keeps one member per protein, by precedence: library reference
before expression-rescue, isoform, ortholog or giant truncation; L1 before L2; then fewer
residues off the source record. It marks the others `retired_reason` in
`candidate_manifest_tranches.tsv`. A retired member keeps its design id and gets no member id
(see "Member ids follow plate order"), and `library_design` readers skip it. The full list with rationales is
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

The causes are fixed where the panels are built, in `library/l1/proteoform-selection/`.
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
and gene class), both in `library/l1/proteoform-selection/`. Their inputs were first read
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
`panel_rejections.tsv` and reruns the isoform sampler until a pass adds no rejection. It no
longer reruns the ortholog sampler: the L1 ortholog panel is the Compara 116 adoption (next
section), and the original ortholog panel's rejections are kept as recorded. The table below
is the first redraw, of 2026-11-01; the gated redraw of 2026-11-12 is described in "Isoform
forms are gated before the draw". The ten no-CDS picks were never L1 members, which is why the panels
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
design id, and its L2 row is retired `moved_to:<design id>`; no existing design id changes.
Rerunning `split_library_tranches.py` would reassign every design id, so it must not be used
for this.
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

### Isoform forms are gated before the draw

Every form in the isoform cohort (1,894) is tested before the panel draw by
`gate_isoform_cohort.py`, and the isoform sampler treats a failing form as unusable: a failing
pick is replaced within its stratum by the same seeded rule, and a failing form is never drawn
as a replacement. A form passes only if both criteria hold.

| Criterion | Test |
|---|---|
| Native start | its protein starts with M, and an Ensembl 115 transcript whose CDS starts ATG encodes exactly that protein (exact peptide match, no allele tolerance) |
| Reference N-terminus kept | by DeepTMHMM v1, against the gene's L1 reference protein: the form keeps the reference's signal peptide, its N-terminal orientation after any signal peptide, and, for a non-GPI reference, the reference's first TM helix sequence called membrane over at least half its residues |

Both proteins are matched to the pinned topology snapshot by sequence. A form with no v1 row
fails as unverifiable rather than being predicted here. The gate exists because a panel form
that starts at a non-native residue, or that loses the signal peptide or first helix the
reference uses for placement, tests a different architecture from the gene it stands for.

| Tier | Forms | Pass | Fail | Of the fails, unverifiable only |
|---|---:|---:|---:|---:|
| observed | 795 | 511 | 284 | |
| prediction-only | 1,099 | 291 | 808 | |
| all | 1,894 | 802 | 1,092 | 931 |

Most prediction-only forms have no v1 row in the snapshot. Many of them were run through
DeepTMHMM v1 for the isoform universe (`topology_shortlist.3line`), but those strings are not
in the topology snapshot, so the gate cannot read them. Publishing them to the snapshot would
make those forms drawable; until then the prediction-only tier draws from 291 forms.

The gate tests only the N-terminal architecture, so a form can pass it and still lack the
domains that make the protein what it is. Five such forms are rejected by name in
`panel_rejections.tsv`, which carries `decided_by`, `implemented_by`, `alternative` and `sources`
per row; the sampler replaces them in their stratum by the same seeded rule. Each decision comes
from the per-member proteoform audit (`proteoform_decisions_applied.tsv`, copied verbatim from
`l1_proteoform_decisions.tsv`).

| Form | Protein | Why it is not a surface isoform of its gene |
|---|---|---|
| ABCC3 O15438-5 | 572 of 1,527 aa | canonical 1–546 plus a 26-aa novel tail: TM helices 1–10 only, and neither nucleotide-binding domain nor the second transmembrane domain (UniProt O15438 features) |
| OPN3 ENST00000469376 | 189 of 402 aa | canonical 1–181 plus 8 divergent residues: TM1–4 only, without TM5–7 or the retinal-binding Lys299 |
| CLEC7A Q9BXN2-9 | 77 aa | cytoplasmic tail, TM and a stalk fragment; no C-type lectin domain |
| SLC6A11 P48066-2 | 208 aa | 3 of the transporter's 12 TM helices (33% of the reference) |
| NECTIN2 ENST00000917943 | 537 aa | one residue deleted inside the signal peptide, so the mature protein is the reference's |

No length floor is applied: a 50% floor would remove HLA-G's form (46% of its reference), which
is kept on evidence. SLC9D1's and KCNN2's short forms (19% and 29%) were considered too, and are
out of the panel on their own defects, not on length: the SLC9D1 form is a 126-aa internal
fragment with a non-ATG start and no TM helix, and the KCNN2 form lacks S1–S5 (l1-proteoform-audit,
PR #305).

Gating and these rejections have each been followed by one redraw (`panel_changes.tsv`, with
`decided_by` / `implemented_by` / `reason` per row). A form that passes the gate can still leave
the panel, because the replacement draw continues one seeded sequence in draw order: a new
rejection early in the order changes which forms later replacements land on. The current panel's
150 forms all pass the gate; every form drawn in the second redraw (`l1-p3141`–`l1-p3152`) is
at least 95% of its reference's length.

### Reference-form fixes

The reference ladder names a gene's form outright through `reference_overrides.tsv`; each row
carries `decided_by` and `implemented_by`. `apply_reference_fixes.py` reruns the ladder and the
CDS resolver from pinned inputs and writes the ledger `reference_changes.tsv`.

| Gene | Was | Now | Tag lane after |
|---|---|---|---|
| CD200R1L | ENST00000488794, 250 aa, no signal peptide | Q6Q8B3-1, ENST00000398214, 271 aa | N after SP (was N at residue 1, type III) |
| LRRC37B | ENST00000543378, 865 aa, no signal peptide | Q96QE4-1, ENST00000341671, 947 aa | N after SP (was type III) |
| CACHD1 | ENST00000650260, 1,140 aa, no signal peptide | Q5VU97-1, ENST00000651257, 1,274 aa | N after SP (was type III) |
| ADGRF3 | ENST00000651242, 1,004 aa, no signal peptide | Q8IZF5-1, ENST00000311519, 1,079 aa | N after SP (was N at residue 1) |
| RGMB | ENST00000308234, 478 aa, 41 residues upstream of the signal peptide | Q6NW40, ENST00000513185, 437 aa | N after SP (was internal, GPI-rerouted) |
| SLC22A17 | CTG start | the same CDS with ATG; no residue added | internal (unchanged) |
| TRPV6 | ACG start | the same CDS with ATG; no residue added | internal (unchanged) |
| OR8G3 | GTG start (both Ensembl 115 transcripts) | the same CDS with ATG; protein = P0DMU2, 310 aa | N at residue 1, multipass (unchanged) |
| EPHA6 | ENST00000389672, 1,130 aa: 94 residues upstream of the signal peptide, which v1 then reads as a TM helix | the same CDS from the codon of Met95; protein = Q9UF33, 1,036 aa | N after SP, 24 aa (was internal) |
| LMLN | ENST00000482695, 452 aa, starts at canonical residue 233 without the signal peptide | Q96KR4, ENST00000330198, 647 aa | N after SP, 43 aa (was type III) |
| DISP2 | ENST00000949524, 1,184 aa, 10 of 12 TM helices | A7MBM2, ENST00000267889, 1,401 aa | internal (unchanged) |
| TSPAN10 | ENST00000696047, 235 aa, 3 of 4 TM helices | Q9H1Z9, ENST00000611590, 355 aa | internal (was C-terminal) |

No Ensembl 115 transcript encodes Q9UF33, so EPHA6 is an in-frame 5′ trim of its MANE
transcript, with no sequence added. LMLN's canonical is preferred over Q96KR4-3 (already built,
684 aa): both keep the same 43-aa signal peptide and both are GPI-anchored by NetGPI, and
isoform 3 differs only by a 37-aa insertion after residue 403 of unknown function. The L2 forms
that the new CD200R1L, LRRC37B, CACHD1, DISP2 and TSPAN10 references duplicate are retired
`duplicate_of:` their reference. TSPAN10's linker encoding, issued in batch 1 while it was C-terminal, is unused (TSPAN10
is now an internal-loop arm) but stays reserved (`l1_swap_risk.RELEASED`), like the 21 snorkel
members held on 2026-11-13. It is reserved rather than released for two reasons. An issued
encoding is never handed to another member, so no two inserts in L1's history ever share one.
And batch 1 replays only with TSPAN10's flank in the draw. Later batches are drawn around it. The designed-CDS and construct records are versioned
under `library-candidates/2026-11-12-reference-fixes-2/`, and the reference record under
`library-design/reference-selection/2026-11-12-reference-fixes-2/`.

#### Cyno orthologs restarted at the human Met1

Four cyno Ensembl 116 models start upstream of the residue aligned to human Met1, at an
in-frame native ATG; the extension has no human counterpart, no signal peptide and would carry
the N-terminal tag. Each gene has only that one model, so the restart is a recorded exception to
the ortholog CDS policy's ban on start repairs (`ortholog_start_trims.tsv`;
`apply_ortholog_start_trims.py`). The restart drops the upstream codons and edits none. Each
restart's CDS and protein hash match the ones the ortholog review computed independently.

| Member | Restart | Protein | Status |
|---|---|---|---|
| cyno CCR5 `l1-p2959` | Met28 | 352 aa | applied |
| cyno CHRM3 `l1-p2938` | Met32 | 590 aa | applied 2026-11-12, after its v1 row was published (`topo_2026_11_12_cyno_restarts`) |
| cyno OR3A2 `l1-p2963` | Met33 | 315 aa | applied 2026-11-12, after its v1 row was published (`topo_2026_11_12_cyno_restarts`) |
| cyno NIPAL4 `l1-p2962` | Met215 | 404 aa | applied 2026-11-12, after its v1 row was published (`topo_2026_11_12_cyno_restarts`) |

A restart is applied only when its protein has a DeepTMHMM v1 row in the pinned topology
snapshot, because every active member needs one for its tag path. The other 196 active L1
orthologs are unchanged.

### Lane overrides

DeepTMHMM v1 is the topology of record, but for four members structures or homology contradict
its N-terminal call. `lane_overrides.tsv` (keyed on design id, bound to the protein's SHA-256)
sets their placement after classification in `l1_tag_paths.py`; `rung_before_gpi` keeps the v1
call and `lane_override` carries the reason.

| Member | v1 lane | Placement | Evidence |
|---|---|---|---|
| MFSD2A reference `l1-p1213` and isoform `l1-p2466` | N at residue 1 | internal | cryo-EM of mouse ([Wood et al., 2021](https://pubmed.ncbi.nlm.nih.gov/34349262/)) and chicken ([Cater et al., 2021](https://pubmed.ncbi.nlm.nih.gov/34135507/)) MFSD2A put the N-terminus in the cytoplasm and the TM5–TM6 loop outside; DeepTMHMM v2, the model and the K/R bias agree |
| AQP12B `l1-p0128` | N at residue 1 | internal | the aquaporin fold puts both termini in the cytoplasm; UniProt A6NM10 annotates 1–22 cytoplasmic; v1's extra helix at 8–21 is weak. AQP12 is intracellular, so surface display is doubtful on any lane |
| GLP2R `l1-p0724` | N at residue 1 | N after SP, 66 aa | mouse Glp2r's signal peptide aligns to human ~43–66 and v2 calls SP 1–66; an HA tag after Met1 gives ~0.12× the median surface signal in HEK293T ([Tedman et al., 2026](https://pubmed.ncbi.nlm.nih.gov/42736269/)), consistent with a cleaved N-terminus |

GLP2R is also an expression-rescue gene, so the contrast with its Igκ-leader arm
(`rescue:GLP2R`) is now between two constructs tagged after a leader. SLC61A1 (`l1-p1994`),
TSPO2 (`l1-p2302`), SLC44A1 (`l1-p1961`) and TMC2 (`l1-p2158`) keep their v1 lanes; the decision
log says why.

### Member ids follow plate order

Every L1 member has two ids:

| Id | Format (manifest column) | Assigned by | Changes when | Records keyed on it |
|---|---|---|---|---|
| Design id | `l1-pNNNN`, `l2-pNNNN` (`design_id`) | the design stage that created the row: `split_library_tranches.py`, `redraw_panels.py`, `add_library_genes.py`, panel adoption | never; a retired id is not reused | design-stage records: panel changes, rejections and replacements, retirements (`duplicate_of:`, `moved_to:`), adoption ledgers, completeness and CDS audits, fold inputs and their structure files, topology supplement rows |
| Member id | `l1-mNNNN` (`candidate_id`) | `renumber_library_members.py` | any membership change, until the order freeze | the current-library records listed below, and the plate order |

Member ids run in order of the *ordered* insert's length (`ordered_insert_nt`: tagged and
synthesis-recoded, as Twist makes it, from `synthesis_compliance.tsv`), shortest first, with
ties broken by design id (decided 2026-11-14 by the coordinator under Becca's delegation; until
then the untagged `insert_nt`, which is 59–359 nt shorter). That is the order the library is
ordered in: purified plasmid in 96-well plates, so each plate holds members of similar length
(reviewer brief, section 8). From the per-insert cut (`2026-11-15-final-e`) the plate position
belongs to the insert: insert N (`l1-iNNNN`, `l1_insert_numbering.tsv`) sits on plate ⌈N/96⌉,
wells column-major (A01, B01, …, H01, A02), sorted by ordered insert length then arm id. The
3,255 inserts fill 33 plates and 87 wells of a 34th; `l1-i0001` is CD24 (`l1-p0277-a1`, 527-nt
ordered insert) and `l1-i3255` is NRXN2 (`l1-p1330-a1`, 5,423 nt). Held members have no insert
and no plate. Members are named by their design id; the `l1-mNNNN` member id is retired and
survives only in dated records. The two formats differ on purpose, so an old id cannot be read
as a member id: `l1-m0026` and `l1-p0026` are different members, and a member id has no `p`
to be mistaken for a P-series construct number. Retired L1 rows and every L2 row have no
member id; their `candidate_id` is their design id.

`l1_member_renumbering.tsv` is the mapping: design id, member id, plate, well, `insert_nt`,
manifest key, and the SHA-256 of the member's designed CDS (for the 25 expression-rescue arms,
the IgK-leader + receptor construct in `rescue_constructs.tsv`, since their CDS is not in
`candidate_cds.json`). These records describe the current library and are keyed on member ids:

| Record | Id columns re-keyed |
|---|---|
| `candidate_manifest_tranches.tsv` | `candidate_id` (`design_id` added beside it) |
| `l1_tag_paths.tsv`, `l1_cassette_encodings.tsv` | `candidate_id` |
| `l1_swap_risk.tsv`, `l1_swap_risk_summary.json` | `candidate_id` and the three `partner_*` columns; the summary's member ids |
| `l1_structure_gap.tsv`, `proteoform_match.tsv`, `residue_offsets.tsv`, `structure_coverage_baseline.json` | `candidate_id`; the baseline's member lists |
| `library_topology_index.tsv`, `library_manifest_l1.tsv`, `library_structure_payload.tsv`, `l1_internal_members.tsv`, `l1_deeptmhmm_v1_topology_gap.tsv` | `candidate_id` |

Only whole id cells change. `topology_acc` is left alone: for a supplement row it is the
design id under which that protein sits in the topology snapshot, a pointer into another
record rather than a member id. Row counts are asserted equal, and every row that carries a
`sequence_sha256` is joined back to the designed protein of the member it now names. Designed
CDS and constructs on S3 (`candidate_cds.json`, `candidate_constructs.json`) are keyed on the
manifest `key`, not on an id, so they did not change. The re-keyed records are published
under `library-design/member-renumbering/2026-11-07/`; the older copies at the experiment
root on S3 carry design ids and are not overwritten.

**Membership is not frozen.** Five genes lack an intact Ensembl copy, and the order freeze
has not been set. Design stages therefore work only on design ids: they read the
manifest with `candidate_id` set to the design id and write it back unnumbered.
`renumber_library_members.py --write` then renumbers from the manifest and re-keys every
record above, including records still carrying member ids from the previous numbering, which
it maps back through the committed mapping. Until it runs, `library_pipeline.py verify`
fails at stage 15. The numbering that ships is the one taken at the order freeze: setting
`ORDER_FREEZE` in the script records that date, after which `--write` refuses to move any
assigned member id.

```bash
uv run --locked python renumber_library_members.py --check   # committed state is current
uv run --locked python renumber_library_members.py --write   # after any membership change
```

**Join on design id, never on member id** (decided 2026-11-13 by the coordinator under Becca's
2026-11-12 delegation). A member id is the plate position: it is authoritative only in
`l1_member_renumbering.tsv`, the order plan and the records in the table above. Those are the
only tables re-keyed when membership changes, along with the rung arm tables, which their
builders rebuild. The dated audit and screen tables below keep the member id current when
each was written, and are not re-keyed. Rewriting them would change checksums that verify
stages pin, for no change in content, and member ids keep moving until the order freeze. At
the freeze, one final re-key of everything is reasonable. Rejected: re-keying them all now;
freezing ids early with gaps. After the ordered-length renumber of 2026-11-14, 36 dated
tables carry member ids from an earlier numbering. They are listed by name in
`tests/test_member_id_policy.py` (`DATED`), every one carries `design_id`, and the test fails if
any other table names a member by an id that is not current.

Three tables carry a member id with no design id: `constructs_l1_type3_proteins.tsv` (join on
`precursor_sha256`), `membrane-side-gate/calibration_tmdet.tsv` (a calibration set, not a join input) and
`l1_rbp_motif_cassettes.tsv`, which has no stable key and must be regenerated keyed on design id
before it is joined. `tests/test_member_id_policy.py` fails if any other table carries member
ids without a design id, or names a design id that does not exist.

No L1 barcode is reserved yet (`barcode_reservations.tsv` holds only panel rows), so the
renumbering needed no barcode registry change.

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
| 5 | ortholog repick | each gene's mouse and cyno proteoform in the April pool (source of L2 orthologs) | byte for byte |
| 5b | universe coverage | which library genes the isoform and ortholog universes could draw from | every recorded ortholog pair is 1:1 high-confidence in Compara 115; coverage table byte for byte (see below) |
| 5c | isoform gate | which cohort forms the isoform draw may use | `isoform_gate.tsv` byte for byte, and equal to its pinned key |
| 6 | panels | the 150-isoform draw (through the gate) and the original 200-ortholog draw (retired from L1) | byte for byte (original draw and rule-based redraw) |
| 6b | ortholog adoption | the L1 ortholog panel: CDS-eligible draw, review swaps, member map | final-v2 draw replays to the panel; adopted panel, member map, manifest and appended CDS/construct records byte for byte |
| 7 | rescue arms | 25 expression-rescue arms, their CDS and constructs | byte for byte |
| 8 | giants | truncations of oversized genes and their ORFs | byte for byte |
| 9 | candidate CDS | each candidate's designed CDS | all 7,020 resolver-recorded CDS except the 5 named references, which equal the current record; 2 approved ATG starts; 56 appended |
| 9b | reference fixes | the named references, ATG start replacements and EPHA6's in-frame trim | ledger, CDS/construct and reference records byte for byte; the first round's record differs only in the second round's entries |
| 9c | ortholog start trims | cyno orthologs restarted at the native Met aligned to human Met1 | every restart's hashes; applied only with a v1 row; the other orthologs unchanged; records byte for byte |
| 10 | constructs | domesticated ORFs and p0039 vectors | all 7,076 ORFs and 7,125 vectors (about 30 s) |
| 11 | tranches | L1/L2 membership; duplicate proteins retired | no two active members share a protein |
| 12 | topology | DeepTMHMM v1 per designed protein | every active member has one |
| 13 | tag paths | each member's placement rung, with `lane_overrides.tsv` applied | byte for byte |
| 14 | cassettes | each C-terminal member's recoded linkers | every issued batch replays to the recorded encodings |
| 15 | member ids | each active L1 member's `l1-mNNNN`, in plate order | recomputed from the manifest; mapping and every re-keyed record current |
| 16 | C-terminal review | which C-terminal members keep a C-terminal tag on the cell (`cterm_tag = removed` when a C-terminal peptide is cut off) | `l1_cterm_review.tsv` byte for byte from the UniProt snapshot of 2026-11-12 |
| 17 | C-terminal arms | each C-terminal-rung member's insert: ORF, cassette, TAA, in the p003 layout (ITM2B: cassette after residue 235), plus the snorkel arms | all 256 inserts byte for byte; protein = native[:n] + cassette + native[n:]; no new Type IIS site |

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
| DeepTMHMM | 1.0.24 and 1.0 academic (v1), topology snapshot 2026-11-06-ortholog-cds-policy | 4, 12, 13 |
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

Named overrides in `reference_overrides.tsv` take precedence. Eighteen genes have one in the
current reference record (2026-11-12-reference-fixes-2): six GPI forms, MICA, three
reference-allele loss-of-function genes, five signal-peptide forms (CD200R1L, LRRC37B, CACHD1,
ADGRF3, RGMB) and the UniProt canonical for DISP2, TSPAN10 and LMLN. The brief's section 3
tabulates them with their reasons and dates.

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

**Adding a gene.** Members are appended under new design ids, so a gene admitted later
cannot shift any existing member's design id. It does shift member ids, so the renumber
stage reruns afterwards (see "Member ids follow plate order"). List the gene in `gene_list_additions.tsv` with the
decision, the date and the surfaceome release that admitted it, then run
`add_library_genes.py`. It runs the reference, CDS and construct steps (3, 9 and 10) over
the list from the same pinned inputs. It fails unless every recorded reference row and CDS comes back unchanged.
The new member gets the next free L1 design id, and its designed protein must not repeat an active
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
uv run --locked python export_internal_members.py   # internal members follow the tag paths
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
`library/l1/proteoform-selection/check_reference_allele_lof.py`, which fails unless every GRCh38
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

### Synthesis compliance: Twist's report and the recodes that clear it

Twist scored the first L1 draft (untagged ORFs in their Esp3I wrapper, no cassette, no
barcode): 2,290 Low, 316 Medium, 104 High and 25 Not Accepted. Its report is stored as vendor
output at `external://twist/l1-complexity/2026-11-12/l1_insert_complexity_analysis.csv`
(SHA-256 `b95071e0…`). [`twist_complexity_map.py`](twist_complexity_map.py) joins every row to
the current manifest on its design id → `twist_complexity_map.tsv`. Of the 25 refused, 24 are
still in L1 with the CDS Twist scored, and 1 is retired. Of the 104 High, 101 are unchanged and 3 are retired.

[`synthesis_compliance.py`](synthesis_compliance.py) (pipeline stage `synthesis_compliance`)
checks the inserts we will order: each rung's tagged arm as built (the five tag-placement
rungs at integration commit `355b524d`, and the 25 expression-rescue arms at `92a1dfa4`), in the
vendor's wrapper, with the barcode as `N`. Its rules are a local approximation of Twist's, calibrated on
Twist's calls for the rebuilt draft (`--calibrate`). They agree on all 25 refusals and err
strict elsewhere:

| Rule | Local check | Calibration |
|---|---|---|
| Tn5 mosaic end | edit distance ≤ 3 to the 19-nt ME, either strand | Twist's 13 "impermissible" hits are all at 3; every accepted insert is at ≥ 4 |
| 100-bp GC | refused above 90% | Twist refused windows of 91–98% |
| 50-bp GC | flagged above 80% | all 568 of Twist's flags, plus 134 it did not flag |
| Direct repeat | refused at ≥ 200 nt of ≥ 90% identity, tandem copies included | CLEC4M is the only refusal; Twist's rule is not reproducible from sequence (CD209 and MUC21 carry longer repeats and scored High) |

For each insert that misses a target it makes the fewest synonymous codon changes. The targets are:

- every 50-bp window at or below 80% GC, and 75% for windows that span the cassette's last base;
- ME distance ≥ 4;
- no exact repeat ≥ 40 nt and no ≥ 90% repeat ≥ 100 nt.

Only codons with human usage ≥ 10% are used. It never changes Met1, the cassette, one codon either side of it, an Igκ leader, or the signal peptide plus 3 codons. It rejects any edit that adds a BbsI, Esp3I or BsaI site, a consensus splice donor or an ME match. Protein, cassette DNA, wrapper and readout block are asserted unchanged.

| Inserts (3,255 rung arms at 8432c010) | Before | After |
|---|---:|---:|
| Not Accepted (local) | 46 | 0 |
| High (local) | 174 | 24 |
| Medium (local) | 676 | 255 |
| Low (local) | 2,359 | 2,976 |
| Cassette-end 50-bp window over 75% GC | 204 | 2 |

806 inserts change, by 5,319 codons in all: 4,668 for 50-bp GC, 458 for the cassette-end margin, 179 for repeats and 14 for ME matches. All but 24 of the 3,255 inserts pass every target after recoding; those 24, all native N-after-SP arms, keep a 50-bp window of 90–98% GC inside the protected signal-peptide codons and score High. Two inserts (l1-p0222-a1, l1-p1330-a1) keep a 76% window across the cassette's end, over the 75% margin and under the 80% target. No provisional L1 barcode occurs as a 20-mer in any of them. Medium is far rarer than in the previous order (255 against 938), consistent with the new linker encodings: p003's cassette repeats 17 nt of its own `GGGGS` codons, and no issued encoding repeats more than 9 (Medium is driven by exact repeats of 18–38 nt).

Every changed codon is in `synthesis_recodes.tsv`, one row per codon with its reason, `decided_by` and `implemented_by`. Per-insert metrics before and after are in `synthesis_compliance.tsv`. **The recoded inserts are the order of record** (superseded through subsequent pre-freeze cuts to the current order of record `processed://sprint4/s4e9-orf-design-fg-panel/synthesis/2026-11-25-pre-freeze-e/l1_order_inserts_recoded.fasta`, `3,352` inserts across `35` plates, SHA-256 `658682f8…` — see [`synthesis_compliance_summary.json`](synthesis_compliance_summary.json) and [`reviewer-package/BRIEF.md`](reviewer-package/BRIEF.md)). The historical `2026-11-21` cut documented in the table above was `processed://sprint4/s4e9-orf-design-fg-panel/synthesis/2026-11-21-snorkel-modules/l1_order_inserts_recoded.fasta` (SHA-256 `66332be4…`), with each insert named by its arm id, plate position and design id. It superseded nine earlier files, kept unchanged: `synthesis/2026-11-21-arm-linkers/` (8b1c4ed3…, the same inserts with the 109 snorkel arms on one module encoding; Becca, 2026-11-21: "Yes, recode them"), `synthesis/2026-11-15-final-e/` (20e8f275…, the same 3,255 inserts and proteins with the 2,942 N-terminal and internal arms on p003's one cassette encoding; Becca, 2026-11-21: "we still don't want shared linkers in the middle", so every arm now has its own linker encoding, `diversify_arm_cassettes.py`), `synthesis/2026-11-14-final-d/` (1ee73f40…, the interim one-insert-per-member cut, 2,621 inserts), `synthesis/2026-11-14-final-c/` (8c5fdc3a…, cut from `main` at 903f9135 before the C-terminal follow-up held SLC41A2, TRPV1 l1-p2528 and VTCN1 human and cyno; 2,625 inserts), `synthesis/2026-11-14-final-b/` (dd65ff76…, cut from `main` at d6ab5093 before mouse IL17RC was built), `synthesis/2026-11-14-final/` (b6ce7bc9…, cut from `main` at cdb78068 before the last rung changes), `synthesis/2026-11-13/` (d1257f0a…, 2,716 inserts including the 25 expression-rescue arms retired that day), `synthesis/2026-11-13-one-arm/` (1f29936b…, named by pre-renumber member ids) and `synthesis/2026-11-13-one-arm-renamed/` (b426e433…, 2,678 inserts built from the draft integration branch, not `main`). The recoder is deterministic, and `--check` reproduces every table byte for byte. The rung builds themselves are unchanged: the recode is a layer applied at export, so each rung's `--check` still holds. The barcodes, assigned later, are outside the recoded sequence. The swap-risk tract measure covers only the cassette encodings, which do not change. `l1_pool_qc.py` assembles every insert into p0039 and runs the plasmid checks on each (CI: `tests/test_l1_pool_qc.py`).

### RNA-binding-protein motifs: the shared elements are clean

Every active L1 member (2,735 at the scan, including the 25 expression-rescue members retired since) is assembled in silico into p0039 with the order-time insert
(Esp3I | Kozak | designed CDS + stop | pad | 5' handle | barcode | gap + RT | 10x CS1 | Esp3I)
and its EFS transcript, from the EFS end to the 3' LTR R end, is scanned for the motifs that
RNA-binding proteins read to lower expression: Pumilio PRE `TGTA[ACT]ATA`, CFIm `TGTA`
clusters, AREs, CPEs, cryptic poly(A) signals and splice sites. The motif panel, its flag
rules and the reason for the screen (four tiled PREs in the miR-L no-target arms) are in
the construct RBP-motif screen ([`rbp_motif_screen.py`](rbp_motif_screen.py)). Natural CDSs carry
some of these motifs, so hits are split by where they sit:

| Where | Result | Table |
|---|---|---|
| Shared: 5' leader, insert tail, p0039 backbone transcript | 11 hits, identical in every member, none flagged in designed sequence: 7 `TGTA` (3 IRES2, 2 WPRE, 1 3'UTR, 1 LTR), 2 `GT[AG]AG` (IRES2), 1 acceptor-heuristic hit (WPRE), and the LTR's own AATAAA. No PRE, ARE or CPE | `l1_rbp_motif_shared.tsv` |
| Junction: across the ORF boundary | 54 members, all one acceptor-heuristic hit: a pyrimidine-rich ORF end, the pad `CC` and the handle's leading `AG` read as tract + `CAG` (51 placed at the pad, 3 inside the ORF end) | `l1_rbp_motif_members.tsv`, `junction_hits` |
| ORF-intrinsic | 86 PREs in 85 members (SLC12A2 has 2), 49 ARE nonamers, 1,395 CPEs, 1,480 `AATAAA`/`ATTAAA`; 10 members carry a motif in identical repeats (e.g. CD209 and CLEC4M neck repeats, 6 copies of a donor) | `l1_rbp_motif_members.tsv`, `orf_*` |
| Issued tag cassettes (218) | no hit of any motif | `l1_rbp_motif_cassettes.tsv` |
| Barcode pool (s1e2, 28,506), in handle/RT context | no PRE, ARE, CPE or poly(A) signal; 1,409 carry `TGTA` (19 of them twice); 874 trip the acceptor heuristic, 755 of them with a tract that starts in the handle's pyrimidine-rich 3' end | `l1_rbp_motif_summary.json` |

The ORF-intrinsic hits are reported, not acted on: they are part of each protein's native
coding sequence, and the screen has no evidence that any of them lowers expression of these
constructs. Members are not yet assigned barcodes, so the barcode slot is masked in the member
scan and the pool is screened on its own. Regenerate with `rbp_motif_screen.py`.

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
9. Refold TSPAN10 (`l1-p2284`) and SLC22A17 (`l1-p1882`), whose ColabFold models are of earlier designed proteins (`stale_fold` in `library_structure_payload.tsv`).
10. Publish DeepTMHMM v1 rows for the restarted cyno CHRM3, OR3A2 and NIPAL4 proteins, then set their `ortholog_start_trims.tsv` status to `applied` and rerun `apply_ortholog_start_trims.py`.

### Decision log

| Date | Was | Now | Why | Decided by / implemented by |
|---|---|---|---|---|
| 2026-11-12 | Inserts synthesized as built; Twist refused 25 of the first draft | A synthesis-compliance layer recodes 751 tagged inserts synonymously, with every codon recorded, and all 2,691 meet the local targets | Twist's complexity report on the first draft (refusals for Tn5 mosaic-end matches, 100-bp GC and one repeat) | l1-reviewer-package (delegated by Becca 2026-11-12, "use your best judgment") / Claude Code session l1-reviewer-package, PR #283 |
| 2026-11-10 | Every L1 member's cassette deferred | Every C-terminal-rung member is checked one by one and carries its issued cassette in a built insert: at the C-terminus, except the ITM2B members, which get it after residue 235 | The cassette was fixed already, so only the junction remained. Other risks are flagged, not rerouted, because no other rung exists for these members | Becca (the brief for this build) / Claude Code session claude/l1-cterm-build-16eeeb, PR #284 |
| 2026-11-10 | ITM2B held from the C-terminal rung | Cassette after residue 235, before the furin site | Furin releases the C-terminal Bri23 peptide, so a C-terminal tag is lost. Directly after R243 the cut motif would sit in front of the tag. Residue 235 is past the BRICHOS fold and leaves the 8 residues before the cut contiguous | Becca, 2026-11-10 / Claude Code session claude/l1-cterm-build-16eeeb, PR #284 |
| 2026-10-09 | Isoform panel drawn uniformly at random | 100 observed + 50 prediction-only, balanced across type I, type II, GPI and multipass | A uniform draw loses the rare architectures, which are the ones whose placement is least certain | not recorded / implementer not recorded (before the 2026-11-10 attribution convention) |
| 2026-10-11 | Library constructs unnamed | `l1-pNNNN` | Readability against the `p00NN` panel series | not recorded / implementer not recorded (before the 2026-11-10 attribution convention) |
| 2026-10-11 | GPI called per gene from UniProt keyword | NetGPI per form, cross-checked against UniProt with documented blind spots | The anchored species is not always the reference form — ACHE's canonical is non-anchored at 0.993 while P22303-2 is anchored | not recorded / implementer not recorded (before the 2026-11-10 attribution convention) |
| 2026-10-12 | Internal-rung set reported as 377, from `internal_site_only.tsv` | 515, from both termini cytoplasmic in v1 topology | The 377 held back 176 members as topology-undetermined; all 176 have determined orientations in v1, 116 of them both-cytoplasmic, so nothing was waiting on DeepTMHMM2 | not recorded / implementer not recorded (before the 2026-11-10 attribution convention) |
| 2026-10-15 | Tag rung assigned per canonical gene (2,367) | Per L1 **member** (2,700), on each member's own ordered-proteoform v1 topology | We order members, not genes, and a member's topology is not its gene's — SLC9D1's ordered isoform has 3 TM against the canonical's 13. Internal-site worklist 515 genes → 566 members; restricting to `library_reference` reproduces the per-gene numbers | not recorded / implementer not recorded (before the 2026-11-10 attribution convention) |
| 2026-10-15 | `internal_site_only.tsv` / `internal_site_undetermined.tsv` / `INTERNAL_SITE_ONLY.md` and their generator | Deleted; `TAG_PATHS.md` + `l1_tag_paths.tsv` replace them | Their 377/176 split reads UniProt's missing topological-domain annotations as unknown orientation, so it misleads on sight. Record kept in git history at `8eab499` | not recorded / implementer not recorded (before the 2026-11-10 attribution convention) |
| 2026-10-14 | Structural coverage measured as "does a model exist" | Measured residue-exact against the designed CDS | Length agreement is necessary and not sufficient: 23 of 200 mismatched members are the same length as their model with a different sequence | not recorded / implementer not recorded (before the 2026-11-10 attribution convention) |
| 2026-10-14 | Members without an AlphaFold DB model treated as a reporting gap | Folded with ColabFold and published as our own predictions | A member with no model cannot be scored by any structural lane, so the gap blocks the library rather than annotating it | not recorded / implementer not recorded (before the 2026-11-10 attribution convention) |
| 2026-11-01 | Library steps run from local caches and hand-uploaded tables | One pipeline (`library_pipeline.py`): 7 stages verified from a checksummed S3 snapshot, 5 pinned with the reason | Reproducibility (Becca); freeze-then-verify, so the library that exists is the one reproduced | Becca / implementer not recorded (before the 2026-11-10 attribution convention) |
| 2026-11-01 | Selection stages read the curated GPI union; tag placement read NetGPI | NetGPI per form (`gpi_set.tsv`) for every stage | One definition (Becca). Within the library it adds GFRAL and RAET1G, on curation rather than NetGPI, and changes no selection or panel | Becca / implementer not recorded (before the 2026-11-10 attribution convention) |
| 2026-11-01 | `reference_selection.tsv` published with no PeptideAtlas evidence | Republished with build 607 evidence from the pinned cache | `--offline` had discarded the evidence; no selection changes, 20 genes flagged for review | not recorded / implementer not recorded (before the 2026-11-10 attribution convention) |
| 2026-10-29 | Every codon of the tag cassette recoded per member, ends unique in their last 9 nt | Only the GGGGS linkers recoded; SpyTag003 and P-ALFA-P on p003's codons; the terminal 15-nt GGGGS unique per member | Tags stay on the validated encoding (Becca); recoding them gives no swap benefit (36 nt either way), and the 9-nt end rule runs out at 60 of 216 C-terminal members. The 25-construct order is not changed | Becca / implementer not recorded (before the 2026-11-10 attribution convention) |
| 2026-11-01 | Panel picks the library could not use left out (panels at 144 and 196) | Panels redrawn from a pinned S3 input snapshot; each rejected pick replaced by continuing the same seeded draw in its stratum; panels at 150 and 200 | Reproducibility (Becca): a hand-picked or separately seeded replacement cannot be regenerated by the panel code. The snapshot reproduces the original panels byte for byte | Becca / implementer not recorded (before the 2026-11-10 attribution convention) |
| 2026-10-29 | L1 members accepted whenever gene, accession or transcript differed | One active member per translated designed protein; duplicates retired with ids kept | 7 L1 and 48 L2 members duplicated another member's protein (ACHE, CD177, MMP16, OPN1MW2/3, OR2A42, UPK3BL2 in L1) | not recorded / implementer not recorded (before the 2026-11-10 attribution convention) |
| 2026-11-02 | Gene list frozen at 2,368 | Frozen root plus dated additions (`gene_list_additions.tsv`), appended by `add_library_genes.py`; SLC67A1 added as `l1-p2739` | Becca: adopt the one gene surfaceome 1.3.0 adds. Appending, not re-deriving, leaves every existing member in place, and the script asserts it | Becca / implementer not recorded (before the 2026-11-10 attribution convention) |
| 2026-11-02 | Constructs stage pinned (full rerun is hours) | Every construct and vector re-verified (about 30 s) | A full rerun takes 30 s. It found the S3 constructs record still held MICA's pre-repoint ORF, now republished | not recorded / implementer not recorded (before the 2026-11-10 attribution convention) |
| 2026-11-02 | Gene list, isoform universe, ortholog repick, rescue arms and giants pinned as unreproducible records | All rerun offline from named releases on external:// (surfaceome 1.3.0, Ensembl 115, UniProt 2026_03, APPRIS 2026_06.v50, OpenProt 2.2, PeptideAtlas build 607); the universe and repick at their producing commits | Reproducibility years out (Becca). OpenProt 2.2 and the APPRIS full scores were still downloadable, so not lost. EPPIN's removal is now a recorded exclusion, not a hand edit | Becca / implementer not recorded (before the 2026-11-10 attribution convention) |
| 2026-11-02 | Linker encodings pinned as issued | Each encoding records its issue batch; verify replays the batches from scratch | Shows the 218 encodings follow from the seed and the issue order alone | not recorded / implementer not recorded (before the 2026-11-10 attribution convention) |
| 2026-11-02 | Cyno ECE1 (l1-p2734) listed as needing no fold | Folded from its designed CDS with ColabFold | The coverage baseline had counted its accession as resolved, so the fold-completeness test passed with the gap open | not recorded / implementer not recorded (before the 2026-11-10 attribution convention) |
| 2026-11-04 | KIR2DS4, OR12D1, OR5AL1 without a member (GRCh38 carries a loss-of-function copy) | Built as l1-p2740..2742 from their intact alternative-haplotype transcripts via `reference_allele_lof` overrides | Becca: the intact protein is in Ensembl 115; CDS resolution had searched only by the broken copy's length | Becca / implementer not recorded (before the 2026-11-10 attribution convention) |
| 2026-11-03 | Snorkel cited to PMID 24391846; GFRAL and RAET1G described as NetGPI additions | Snorkel cited to Brown et al. 2013, PMID 24023844; both described as curated calls NetGPI does not support | The old PMID is an unrelated docking paper (Hsin et al. 2013); NetGPI finds no anchored form for either gene. Found by the reviewer-package review, with re-measured counts (191 GPI members, 199 different models) | not recorded / implementer not recorded (before the 2026-11-10 attribution convention) |
| 2026-11-04 | 13 designed CDSs ended on a sense codon: the Ensembl 115 record has no stop (l1-p1054..1056 KIR2DL5A, KIR2DS2, KIR2DS5; L2 l2-p0846 KIR2DS2 isoform and nine cyno orthologs) | `resolve_candidate_cds.py` appends TAA when the record ends without a stop and marks the entry `stop_added`; orf_nt, vector_bp and insert_nt grow by 3 for those rows. Proteins and domestication edits are unchanged | Becca. Without a stop the ORF reads through into the vector. TAA is the stop the C-terminal constructs already use; the native codon cannot be read from the pinned inputs (four KIR records sit on alternate-haplotype scaffolds the primary GTF omits). Previous S3 versions: candidate_cds.json `ftYKdsVZYDD8M49CEt4rDJMn.x.QmKkn`, candidate_constructs.json `Nkbr70dq6lfjWkb.ASYGWZpU8yRocj_Y`. Found by the N-terminal reviewer package's frame check | Becca (2026-11-04) / Claude Code session nterm-reviewer-package, PR #234 bb393b72; per member in `stop_codon_fix.tsv` |
| 2026-11-05 | Ortholog panel described as drawn from library genes with an ortholog | Recorded as drawn from the April 2026 `surfaceome_expressed` query set (1,775 of 2,369 library genes); coverage measured against Compara 115 by `check_universe_coverage.py` | Becca: record the gap and keep the panels. The 418 mouse / 416 cyno eligible genes were never queried because later Compara pulls went to D1, not to the f9c73691 files the universes read; every recorded pair is one-to-one high-confidence in Compara 115 | Becca / implementer not recorded (before the 2026-11-10 attribution convention) |
| 2026-11-05 | Keep the current ortholog panel while documenting its restricted pool | Rebuild ortholog eligibility over all 2,369 current library genes using Compara 116, then redraw with the existing seed and strata after proteoform/topology checks | Becca clarifies that keeping the panel is temporary, and requests the latest Compara release even when newer than 115. The old panel is not the final choice. | Becca / implementer not recorded (before the 2026-11-10 attribution convention) |
| 2026-11-10 | No RNA-binding-protein motif screen | Every member's transcript scanned (`rbp_motif_screen.py`); shared elements reported once, ORF-intrinsic hits per member | Becca: the miR-L no-target arms carry four tiled PREs and underexpress in S6E5; the shared tail and backbone are clean of PRE, ARE and CPE. Decided by Becca (request of 2026-11-08); the cutoffs and the junction class by the Claude Code rbp-motif-screen session, recorded in `l1_rbp_motif_summary.json` `decisions`. Implemented by that session, PR #277 | Becca, 2026-11-08 / Claude Code session rbp-motif-screen, PRs #277 and #290 |
| 2026-11-07 | Member ids (`l1-pNNNN`) are never renumbered; the order would carry separate order numbers (`L1-O0001` …) | Active L1 members are renumbered `l1-mNNNN` in plate order (insert length, then design id); the old id is kept as `design_id` and still keys every design-stage record | Becca: one id per member that matches the plate it sits on. Nothing had been ordered and no L1 barcode was reserved, so no external record named the old ids as members. A new format, not reuse of `l1-pNNNN`, so an old id is never silently read as a different member | Becca, 2026-11-07 / Claude Code session l1-renumber, PR #272 |
| 2026-11-10 | L1 ortholog panel drawn from the April Compara 115 pool; CDS-eligible Compara 116 panel held for review | CDS-eligible panel adopted into L1 with the 11 correspondence replacements; 6 asks kept as drawn; 18 members kept, 126 moved from L2, 56 new (`l1-p2932`–`l1-p3113`), 182 retired; CDS/construct records and topology snapshot versioned | Becca: adopt the reviewed panel. A replaced form lacks a human feature the assay depends on; an ask shows no lost feature | Becca, 2026-11-10 / Claude Code session l1-reviewer-package (subagent), PR #283; correspondence review by l1-ortholog-corresp, PR #273 |
| 2026-10-13 | Walton ≤100 bp tract bar applied library-wide | Bar applies to the pool and the shared cassette, not to isoforms and orthologs | Median canonical-isoform overlap is 1,071 nt with zero of 2,150 pairs under the bar; applying it would remove the panels entirely | not recorded / implementer not recorded (before the 2026-11-10 attribution convention) |
| 2026-11-10 | 61 adopted orthologs with no AlphaFold DB model, unfolded | Folded with ColabFold 1.6.2 from each member's designed protein (61/61, 5.7 GPU-hours); member ids renumbered over the adopted panel by `renumber_library_members.py --write` | A member without a model cannot be scored by any structural lane; the renumber stage runs after every membership change | Becca, 2026-11-10 (approved the GPU spend) / Claude Code session l1-reviewer-package, PR #283 |
| 2026-11-11 | CD200R1L, LRRC37B and CACHD1 references lack the UniProt canonical's signal peptide (type III lane) | References are the signal-peptide forms already built in L2; the L2 rows retired as duplicates | The type III lane rests on an N-terminus the canonical protein does not have | Becca, 2026-11-11 ("Apply all"; "tIII only for now"; review by session type3-risk-review, recorded in PR #282 by session l1-nterm-type3) / Claude Code session l1-reviewer-package (subagent), PR #283 |
| 2026-11-11 | BST1 Q10588-2, CDH13 P55290-4, OPCML Q14982-3 and LYPD5 Q6UWN5-2 in the isoform panel without their gene's signal peptide | Redrawn; all four fail the cohort gate's signal-peptide test, so no explicit rejection was needed | A GPI isoform without the signal peptide cannot reach the surface the way its canonical does | Becca, 2026-11-11 ("Redraw them too"; recorded in PR #282 by session l1-nterm-type3) / Claude Code session l1-reviewer-package (subagent), PR #283 |
| 2026-11-12 | Isoform panel drawn from every cohort form | Cohort gated (native ATG start; reference signal peptide, orientation and first TM helix by DeepTMHMM v1), then redrawn once: 27 out, 27 in | Forms that start at a non-native residue or lose the reference's N-terminal architecture test a different protein from the gene | Becca, 2026-11-12 ("Gate cohort, redraw") / Claude Code session l1-reviewer-package (subagent), PR #283 |
| 2026-11-12 | ADGRF3 lacks its signal peptide; RGMB carries 41 residues upstream of it; SLC22A17 and TRPV6 start on CTG and ACG | ADGRF3 Q8IZF5-1 and RGMB Q6NW40 named; SLC22A17 and TRPV6 start on ATG with no residue added; EPHA6 unchanged (no exact Ensembl 115 encoding of Q9UF33) | Reference-form screen (PR #280, `REFERENCE_FORM_SCREEN.md`) | Becca, 2026-11-12 ("Send all fixes") / Claude Code session l1-reviewer-package (subagent), PR #283 |
| 2026-11-12 | ABCC3 O15438-5, OPN3 ENST00000469376, CLEC7A Q9BXN2-9, SLC6A11 P48066-2 and NECTIN2 ENST00000917943 in the isoform panel | Rejected by name and redrawn once: 12 out, 12 in (`l1-p3141`–`l1-p3152`, each ≥95% of its reference's length). Alternative: keep each form. A ≥50%-of-reference gate rule was considered and withdrawn | Fragments without their gene's defining domains, and one form whose mature protein is the reference's; a length floor would also drop isoforms kept on evidence (SLC9D1, KCNN2, HLA-G) | l1-proteoform-audit (delegated by Becca 2026-11-12, "use your best judgment"), PR #300 / Claude Code session l1-reviewer-package (subagent), PR #283 |
| 2026-11-12 | LMLN, DISP2 and TSPAN10 references are truncated APPRIS-only models; OR8G3 starts GTG; EPHA6 carries 94 residues upstream of its signal peptide | Q96KR4, A7MBM2, Q9H1Z9 named; OR8G3 GTG→ATG (= P0DMU2); EPHA6 trimmed in frame to Met95 (= Q9UF33). Alternatives: keep each model; LMLN Q96KR4-3 | `l1_proteoform_decisions.tsv` per-member audit; LMLN canonical over isoform 3 because isoform 3 adds only an uncharacterized 37-aa insertion | l1-proteoform-audit (delegated by Becca 2026-11-12, "use your best judgment"), PR #300 (LMLN reconciled by l1-reviewer-package with l1-nterm-type3); OR8G3 recorded there as an objective-defect rule / Claude Code session l1-reviewer-package (subagent), PR #283 |
| 2026-11-12 | Cyno CCR5, CHRM3, OR3A2 and NIPAL4 start upstream of the residue aligned to human Met1 | Restarted at the native in-frame Met (Met28, Met32, Met33, Met215), a recorded exception to the ortholog CDS policy; CCR5 applied, the other three held until their v1 rows are published. Alternative: redraw from the Compara 116 pool | The extension has no human counterpart and would carry the N-terminal tag; each gene has one Ensembl 116 model | l1-proteoform-audit and ortholog-correspondence-review (each delegated by Becca 2026-11-12), PRs #300 and 297 / Claude Code session l1-reviewer-package (subagent), PR #283 |
| 2026-11-12 | MFSD2A (reference and isoform), AQP12B and GLP2R on the N-at-residue-1 lane by DeepTMHMM v1 | `lane_overrides.tsv`: MFSD2A ×2 and AQP12B internal; GLP2R N after SP (66 aa). Alternative: the v1 lane | Structures, homology and the Tedman surface signal contradict v1's N-terminal call (lane check, PR #296) | l1-reviewer-package (delegated by Becca 2026-11-12, "use your best judgment") / Claude Code session l1-reviewer-package (subagent), PR #283 |
| 2026-11-12 | SLC61A1, TSPO2, SLC44A1 and TMC2 lanes questioned | No override: SLC61A1 stays N at residue 1 (v1, v2, SignalP 6 and the model agree; another session's case for internal, a 13-residue first helix too short to span the bilayer and UniProt's 12-TM annotation, is recorded here for Becca); TSPO2 stays N at residue 1 (v1 and v2 disagree, no experimental topology; fallback C-terminal if N-in); SLC44A1 and TMC2 stay internal (v1 and v2 both N-in/C-in; UniProt agrees for TMC2, and its C-out call for SLC44A1 has no experimental support) | The topology of record stands where the evidence does not overturn it | l1-reviewer-package (delegated by Becca 2026-11-12, "use your best judgment") / Claude Code session l1-reviewer-package (subagent), PR #283 |

## Provenance

- p003 constants read from `processed://constructs/ordered/p003.gb`. This was a local-only
  file when the note was first written; the map was published to the construct registry on
  2026-10-06, so the gap it described — and the matching one in `s1e2/enzyme_scan.py` — is
  closed for p003.
- Tedman construct design: [Tedman et al., 2026](https://pubmed.ncbi.nlm.nih.gov/42736269/),
  *Efficient experimental characterization of the GPCRome via deep receptor scanning*,
  Nat Commun. [doi:10.1038/s41467-026-76564-7](https://doi.org/10.1038/s41467-026-76564-7).
  (via PubMed, PMC13575135); screen values from
  `deliverome-external://tedman-gpcr-surface-screen/2025-09-19/`.
- Barcodes: `library/l1/vector-design/barcodes_20nt_30k.tsv`.


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
