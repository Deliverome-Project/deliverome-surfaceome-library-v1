# s1e2 — fragmid-vector-design

Design decisions for the modular (Fragmid-style) lentiviral expression + barcoding
construct for the delivery-receptor library. This is a **design doc**, not a data
analysis; the code here backs two quantitative decisions (barcode length, recombination
risk).

> **This document describes the CURRENT design state.** Decisions that were changed or
> retracted along the way are quarantined in the **[Decision log](#decision-log)** at the
> bottom, so everything above it can be read as live. If you are picking this up cold,
> read §1–§3 and then the Next steps.

## Notion page

https://app.notion.com/p/3a5ab86c768381fdbf5dd81d00ae8b19
(ELN: **S1E2** in the Experiments database — kept in sync with this file)

---

## 1. Project parameters

| | |
|---|---|
| First target ORF | **TFRC** (transferrin receptor 1, full length) |
| **Target library size** | **up to 30,000** |
| **Host** | **mammalian cells** — specific line not yet recorded (affects Kozak strength, IRES efficiency, and comparability to EndoNB, which used U2OS) |
| Golden Gate enzyme | **BbsI** (Fragmid); **Esp3I/BsmBI** reserved for downstream guide cloning — ⚠️ **this reservation is now contested**, see §4.11 |
| Parent — domesticated build | **pRDA_722 / Addgene #216077** |
| Parent — p006/p007 fast-iteration | **#236083** (deliberately non-domesticated) |

## Question

What is the canonical module architecture for our Fragmid-style delivery construct, and
what did we decide — and why — on the load-bearing sub-choices: signal peptide,
selection/reporter marker, barcode length and placement, codon usage, stop codons, and
domestication?

## 2. Construct architecture (current)

```
EFS promoter
  → Kozak block  gaccggttctagagccgccacc      (AgeI + XbaI + gccgccacc)
  → [optional SP — ONLY if the mature N-terminus is extracellular; NOT for TfR1]
  → TFRC                                      (native codon usage — see §4.2)
  → GGGGS → SpyTag003 → GGGGSGGGGS → P-ALFA(v1.5)-P → GGGGS
  → TAA[A,G]                                  (stop — see §4.5)
  → HT-recruit fwd primer site → 20-nt barcode → lentiGuide-BC 3′ site
  → IRES2  (Fragmid pF_AA333)
  → mScarlet3 → P2A → PuroR → TAA[A,G]
  → WPRE → 3′LTR / polyA
  → (miRNA control system — NOT yet integrated)
```

## 3. Element provenance

✅ = source recorded · ⚠️ = recorded but unconfirmed/discrepant · ❌ = no provenance

| # | Element | Source | Status |
|---|---|---|---|
| 1 | **EFS promoter** | Fragmid | ✅ |
| 2 | AgeI `accggt` / XbaI `tctaga` | part of the Kozak block | ✅ |
| 3 | **Kozak** `gaccggttctagagccgccacc` | *unstated* | ⚠️ confirm origin. BbsI/Esp3I-clean; −3 = `A` ✅, +4 = `T` (mildly suboptimal vs `gccRccATGG`). Separately `atgtac` = **M-Y**, not TfR1's native start (`MMDQARSAF…`) — open "issue 7" |
| 4 | **TFRC**, full length | Ensembl canonical; cross-checked vs Addgene #61505 (TfR–pHuji) | ✅ C-terminus verified `…DVWDIDNEF`; must be **native codon usage** (§4.2) |
| 5 | GGGGS linker | — | ❌ synthetic, no source; the repeat-diversification target (§5.1) |
| 6 | **SpyTag003** `RGVPHIVMVDAYKRYK` | Howarth ref — Addgene #133450, Keeble 2019 *PNAS* | ✅ AA identical; Twist mammalian codon-opt |
| 7 | GGGGSGGGGS + **P** | — | ❌ no source (the P is the ALFA N-flank) |
| 8 | **ALFA tag** `SRLEEELRRRLTE` | `pF_AA197` / Addgene #216049 | ✅ variant **v1.5**; proline-flanked in context (§4.7) |
| 9 | **P** + GGGGS + stop | — | ❌ no source. ⚠️ stop is currently `TGA` → must become `TAA[A,G]` (§4.5) |
| 10 | **HT-recruit handle** `ggtggATCCggaTCAggctcc` (UTR) | Tycko et al., *Cell* 2020 | ✅ also the barcode amplicon's 5′/fwd primer |
| 11 | 21 nt `tagcactgagtcagatctgag` (BglII) | — | ⚠️ **likely a cloning scar**, not the barcode. Confirm |
| 12 | `CCTTTGGGTAAGCACACGTC` (20 nt) | **designed in-house** | ✅ the barcode. ⚠️ clear the stray `n` nearby |
| 13 | **IRES2** | **Fragmid `pF_AA333`** | ✅ adopted (domesticated). ⚠️ functionally **unvalidated** (§5.2) |
| 14 | **mScarlet3** | FPbase (Gadella 2023) | ✅ |
| 15 | **P2A** | Fragmid `P2A_v1.1` | ✅ (feature map's "T2A" was a mis-annotation) |
| 16 | **PuroR** + stop | pOT_2 / #181971; identical to Fragmid `PuroR_v1.1` | ✅ two independent sources agree. ⚠️ stop → `TAA[A,G]` |
| — | **WPRE** | pRDA_722 / #216077 | ✅ BbsI-clean, but is the **native, non-X-inactivated** element (§5.3) |
| — | **Backbone / 3′LTR / polyA** | pRDA_722 / #216077 | ✅ verified: 2 BbsI (assembly junctions only) + 0 Esp3I |
| — | **snorkel** | — | ⚠️ undecided; would worsen recombination (§5.1) |
| — | **miRNA control system** | — | ⚠️ not integrated (§4.8) |

---

## 4. Key decisions

### 4.1 Signal peptide — only when the mature N-terminus is extracellular

- A cleavable N-terminal SP is only correct for **N-out** topologies (type I
  single-pass, secreted). When we add one, use the **influenza hemagglutinin (HA) SP**
  (+ FLAG) — the Manglik/Kobilka-lineage convention for boosting surface expression of
  N-out receptors like GPCRs.
- **Do NOT add an SP to a type II protein.** **TfR1 is type II** (cytoplasmic N-terminus
  ~1–67, TM ~68–88, ectodomain ~89–761); its TM is a **signal-anchor**, so ER targeting
  is already built in. An added SP would either invert the topology — throwing the
  cytoplasmic **YTRF/Tyr20 endocytosis motif** to the outside and breaking
  internalisation/recycling — or fight the signal-anchor and cause misfolding/ER
  retention. **TfR1 is cloned SP-free with its native N-terminus.**
- An SP is only needed for **ectodomain-only** constructs (secreted or re-anchored on a
  heterologous TM/GPI) — then HA SP or an Igκ leader.

### 4.2 Codon usage — do NOT codon-optimize native ORFs

Synonymous codon choice sets local translation elongation rate, and native codon usage —
including its rare/"slow" codons — encodes a **co-translational folding timing pattern**
that can be load-bearing, especially for multi-pass membrane proteins like TfR1. Bulk
optimization erases that pattern and can yield a protein with the **identical amino-acid
sequence but a different fold, stability, or function.**

- **Kimchi-Sarfaty 2007, *Science*** — a single **synonymous** SNP in the multi-pass
  transporter **MDR1/P-glycoprotein** changed substrate/inhibitor specificity; same
  mRNA/protein levels, altered conformation, attributed to a rare codon shifting
  **co-translational folding and membrane-insertion timing**.
- **Buhr 2016, *Mol Cell*** — synonymous variants of γB-crystallin gave **different final
  conformations** from the identical sequence ("synonymous codons provide a secondary code
  for protein folding").
- **Walsh & Clark 2020, *PNAS*** — synonymous substitutions perturbed the cotranslational
  folding pathway and impaired cell fitness.
- **Bartoszewski 2016** — review on **CFTR and MDR1/P-gp**, the same membrane-transporter
  class as TfR1.
- **Konczal 2019, *PLoS ONE*** — **re-introducing** non-optimal codons into optimized
  constructs *increased* soluble, properly-folded recovery.

**Scope:** native protein ORFs only — **not** the synthetic tags/linkers (SpyTag003,
ALFA, GS-linkers), which have no evolved folding pathway to preserve.

**Consequence:** BbsI/Esp3I sites still must be cleared from native ORFs for Golden Gate,
but via **minimal synonymous edits at only the offending codon(s)** — never whole-ORF
recoding. Non-coding elements (IRES, WPRE) can't use synonymous recoding at all, which is
why we take domesticated versions from Fragmid rather than re-deriving them (§5.2).

### 4.3 Marker cassette — IRES2 → mScarlet3 → P2A → PuroR

A bright fluorescent reporter plus a drug-selection marker in one cistron.

- **mScarlet3** (Gadella 2023) — cysteine-free monomeric RFP with fast/complete
  maturation, **record brightness (QY 75%)**, 4.0 ns lifetime, well-behaved as a fusion
  tag, no apparent cytotoxicity, and specifically strong "as a reporter in transient
  expression systems" — which is what we need for fast iteration.
- **P2A → PuroR** gives puromycin selection stoichiometric with mScarlet3, so you can
  drug-select the population **and** read per-cell expression by fluorescence. P2A leaves
  a scar on mScarlet3's C-terminus (harmless for an FP) and an N-terminal proline on
  PuroR (standard).
- Delivered via **IRES2**, which is frame-independent — so a C-terminal tag on the ORF
  that needs a free C-terminus is unaffected. Trade-off: the IRES cistron runs at
  ~10–50% of the cap-dependent ORF, so the reporter is dimmer than a 2A-linked design
  would be; mScarlet3's brightness offsets this.

**Trade-off accepted:** mScarlet3 **occupies the red channel**. Plan the imaging/FACS
panel accordingly, and note selection is now drug + fluorescence rather than
antibody/MACS. (This reverses an earlier constraint — see the [Decision log](#decision-log).)

#### Alternative: put the selection marker on its own promoter (MORF, REPROcode, RESOLUTE)

Three published barcoded-ORF libraries **decouple the marker from the ORF cistron** instead
of linking it. Worth recording because it is the natural fallback if the IRES A/B in §5.2
fails, and because their barcode choices bear on §4.4.

⚠️ **The table below is user-supplied and not independently verified** — treat the
per-library specifics as to-confirm against the primary sources.

| | scale | architecture | barcode | marker | regulation |
|---|---|---|---|---|---|
| **MORF** (Feng Zhang lab) | 3,548 TF isoforms / 1,836 genes | EF1α → TF ORF → [C-term V5] **plus SV40 → PuroR**, in pLX_TRC317 | unique **24-bp** DNA barcode per isoform; recoverable by bulk sequencing or dial-out from 10x libraries | puromycin | constitutive, 3rd-gen lentivirus |
| **REPROcode** | 408 immune-related TF constructs | SFFV → TF ORF → WPRE → barcode | unique **transcribed** barcode **downstream of WPRE** | **none listed** for the deposited SFFV-MCS vectors | constitutive lentivirus |
| **Human SLC Overexpression Library** (RESOLUTE) | 467 barcoded inserts, **codon-optimized** SLC-centered | TRE/TetO → SLC ORF **plus PGK → rtTA-VP16–2A–PuroR**, in pCW57.1 | ORF-associated barcode **+ UMI** for clone tracking | puromycin | **doxycycline-inducible**, 3rd-gen lentivirus |

**What the separate-promoter design buys — and costs.** MORF's `SV40 → PuroR` and RESOLUTE's
`PGK → rtTA-VP16–2A–PuroR` make marker expression **independent of the ORF cistron**:

- ✅ **It removes our single largest unvalidated dependency.** Our puromycin resistance
  *and* mScarlet3 both sit downstream of an IRES2 we have not functionally validated
  (§5.2). If that IRES under-performs, selection and readout weaken together. A separate
  promoter eliminates that coupling outright.
- ✅ Selection also survives a poorly-expressed or mildly toxic ORF.
- ❌ **You lose linkage as a quality filter — but the linkage is narrower than it sounds.**
  ⚠️ **Correction to an earlier version of this section, which overstated it.** IRES2
  initiates translation *internally and cap-independently*, so a frameshift or premature
  stop **in the ORF does not stop PuroR being made**. What IRES linkage actually certifies
  is that **the transcript exists and is intact through the IRES** — so it filters promoter
  silencing, 5′ truncation and large deletions, **but not the ORF point mutations and indels
  that dominate §5.4's error budget**. A separate promoter forfeits even that.
  **2A would give the strong version** of this filter (same ORF, so an upstream frameshift
  kills the marker) — but 2A is unavailable to us precisely because the ORF needs a free
  C-terminus for the C-terminal tag, which is why IRES2 was chosen in the first place.
- ➡️ **Net:** keep the IRES2-linked design as primary *because* of the linkage, and treat
  separate-promoter as the pre-planned fallback if §5.2's A/B fails. Recorded in §6.

**RESOLUTE's marker cistron is worth copying regardless.** `PGK → rtTA-VP16–2A–PuroR` puts
the transactivator and the resistance gene in **one cistron**, so puromycin selects for
cells with functional rtTA — exactly the thing dox induction depends on. If we ever go
inducible, that is the right way to wire it. Inducibility itself is worth considering if
constitutive overexpression of any delivery receptor proves toxic.

**Bearing on our barcode (§4.4):**

- **MORF's 24 bp** supports our 20 nt on convention grounds — prevailing practice is ≥20.
- **REPROcode places a transcribed barcode 3′ of WPRE**, which is precedent that the
  placement §4.4 leaves open does work, and it is the right choice for 3′ single-cell
  capture — but it maximises the variable-region→barcode span, which is the one thing §5.1
  says to minimise. The two goals genuinely conflict; REPROcode picked capture, we are
  leaning toward linkage.
- **RESOLUTE carries a UMI alongside the barcode.** We do not. Worth deciding explicitly.

⚠️ **RESOLUTE codon-optimized its ORFs**, which is a documented counterexample to §4.2's
"do not codon-optimize native ORFs" — and it is the *worst* case for that concern, since
SLCs are multi-pass TM proteins where co-translational folding matters most. Understand why
they chose it before treating §4.2 as settled.

### 4.4 Barcode — 20 nt, placed between two borrowed primer sites

**Length.** The collision math (`barcode_length_analysis.py`) shows **18 nt is already
more than enough**; we standardize on **20 nt** to match prevailing ORF-BC practice —
convention/interoperability, not because 18 was unsafe.

Precedent: **MORF uses 24 bp**, RESOLUTE and REPROcode both carry per-clone barcodes
(§4.3) — ≥20 is the norm, which is the whole argument for 20 over 18.

- 20 nt → 4²⁰ ≈ 1.1 × 10¹² raw; effective space after our constraints (GC window,
  homopolymer cap, restriction-site + AATAAA + start-ATG avoidance) ≈ **6.5 × 10¹¹**.
- **At the real library size (30,000): expected collision fraction ≈ 4.6 × 10⁻⁸**
  (~0.001 clones). Comfortably over-specified, which is fine.
- Generator: constraint-based rejection sampling in `deliverome_analysis.codon_opt`
  (**PR #10**), regenerated at 20 nt.

#### ⚠️ miRNA seed sites were NOT being excluded — now they are

**The gap.** The shipped constraint set was GC 40–60% · homopolymer ≤3 · no leading `ATG` ·
no `AATAAA` · no BbsI/Esp3I. **Nothing about miRNA seed sites.** That matters because
**our barcode is transcribed**: it sits between the ORF stop and IRES2, an intercistronic
region that behaves like a 3′UTR for cistron 1. A barcode that happens to present a miRNA
**7mer-m8** site can get its own transcript repressed — and the failure mode is a
**barcode-specific expression penalty that looks exactly like a weak surface expressor.**
In a display screen that is very hard to separate from real biology.

How often a random 20-mer carries a site, by how many seed families you exclude:

| distinct seed families | P(7-mer site) | P(8-mer site) |
|---|---|---|
| 108 (broadly conserved) | 7.4% | 1.8% |
| 219 (conserved+) | 15.7% | 3.9% |
| ~2,064 (all human) | ~72% | ~26% |

**Filtering is effectively free** — see `barcode_design.py`. Even the most aggressive option
leaves >10⁷× headroom over a 20,000-clone library, so there is no collision-budget reason
not to do it.

**Two mitigations that reduce but don't remove the risk**, worth knowing before over-engineering:
the barcode sits ~1.6 kb upstream of the polyA (IRES2 + marker + WPRE intervene) and miRNA
repression is strongest near the 3′-most end of a 3′UTR, so sites here are weaker than in a
true 3′UTR; and with ~13 clones per gene each carrying a *different* barcode, a per-barcode
penalty averages out **provided the unit of analysis is the gene, not the barcode.** Score
per-barcode and that protection disappears.

**The decisive argument for filtering now is §4.8:** the miRNA control system is recorded as
*not yet integrated*. Accidental seed sites in barcodes would directly confound a designed
miRNA control system later.

#### The published filter stack, applied to the deterministic barcode

`barcode_filter_stack.py` implements each published rule as a named predicate with its source,
and measures what survives. **Applies to the designed barcode only** — a degenerate window
inherits none of it (see the UMI discussion).

**Sources, and what we could actually verify:**

| Source | Rules | Verification |
|---|---|---|
| **Slutskin et al., Nat Commun 2018** (PMID [29410437](https://pubmed.ncbi.nlm.nih.gov/29410437/)) | 12 nt; differ "by 3 nt or more"; exclude seeds for **highly expressed** miRNAs; exclude `AATAAA` **+ point mutants** | ✅ **Verified verbatim** against the open-access copy, [PMC5802814](https://pmc.ncbi.nlm.nih.gov/articles/PMC5802814/) |
| **FREE — Hawkins et al., PNAS 2018** (PMID [29925596](https://pubmed.ncbi.nlm.nih.gov/29925596/)) | GC bounds; no homopolymer triple; no G→C; no 3-bp self-complementarity | ✅ **Verified against the reference implementation** (`finkelsteinlab/freebarcodes`, `generate.py`) |
| Shiroguchi 2012 | GC 40–60%, homopolymers ≤4 | Strictly weaker than FREE on homopolymers; FREE governs |
| Lyons 2017 | min Hamming distance, homopolymer cap, GC range, blacklists | Our min-Levenshtein-3 is stronger than Hamming |
| O'Brien 2026 preprint | FREE → quality filter → miRBase seed exclusion; 35 barcodes/construct | ⚠️ **Not verified** — preprint cited via ResearchGate. Recorded as an implementation precedent, not a basis for a rule |

**⚠️ Two corrections found by checking the sources:**

The summary was **accurate on Slutskin** but **wrong on FREE** — which is exactly why the
primary sources were worth chasing:

1. **FREE forbids every `GC` dinucleotide, not just `GGC`.** The code is
   `if prev_seq[-1] == 'G': bad_bases += 'C'` — the comment says GGC, the code implements GC.
   That is far stricter, and against a 40–60% GC window nearly self-contradictory. We measured
   both: `GGC` costs 18% of survivors, `GC` costs **75%**. Treated as optional, not adopted.
2. **"No 3-bp self-complementarity" is a hairpin rule, not a k-mer blacklist.** Implemented
   naively — "does any 3-mer's reverse complement appear anywhere" — it rejects `ACGT` outright,
   because revcomp(`ACG`) = `CGT` overlaps it at offset 1. That artifact alone made the filter
   look like it cost 79% of survivors; requiring **non-overlapping arms**, as FREE does, brings
   it to 59%. Our implementation enforces separation.

**Cumulative survival, 16-nt barcode** (Monte Carlo, 120,000 draws):

| filter | cumulative pass | effective space |
|---|---|---|
| GC 40–60% | 54.5% | 2.34 × 10⁹ |
| no homopolymer **triple** (FREE — stricter than our old run ≤3) | 28.9% | 1.24 × 10⁹ |
| no `ATG` anywhere | 21.9% | 9.42 × 10⁸ |
| no `AATAAA` + 18 point mutants | 21.6% | 9.29 × 10⁸ |
| no **BbsI / Esp3I / BsaI**, both strands | 21.1% | 9.06 × 10⁸ |
| no 3-bp self-complementarity, non-overlapping | 8.70% | 3.74 × 10⁸ |
| no `conserved+` 7mer-m8 seed | **7.60%** | **3.26 × 10⁸** |

Optional additions, as fraction of survivors lost: **6mer seed (Slutskin's definition) −35%**
— *since adopted, see below* · `GGC` −18% · `GC` dinucleotide −75% · all-human 7mer seed −62%.

**Self-complementarity is the most expensive rule, not the miRNA filter** — it costs 59% of
survivors where the seed filter costs 13%. Worth knowing, since the seed filter is the one that
feels expensive.

**✅ 16 nt still clears 30,000 barcodes at min Levenshtein 3 under the full stack** — including
Slutskin's stricter 6-mer seed definition. 15 nt reaches only 27,503 (23,318 with the 6-mer),
14 nt only 12,063. So **16 nt is the floor and the published stack costs us nothing at that
length.**

**Adopted (final):** GC 40–60% · no homopolymer triple · no `ATG` anywhere · `AATAAA` + 18
point mutants · **BbsI + Esp3I + BsaI** both strands · non-overlapping 3-bp
self-complementarity · **`conserved+` 6-mer seeds (Slutskin's definition)** · **min Levenshtein
3 AND min Hamming 4**.

**Not adopted:** the `GC`-dinucleotide rule (75% cost, a 2018 Illumina-chemistry concern) and
all-human seed tiers (62% cost, mostly miRNAs not expressed in any given line).

⚠️ **Two rules were tightened after the first build** — the seed width and the distance metric
— because a published method proved stricter than ours on each. The survival table above was
measured with the *7mer-m8* seed and is kept as the per-rule cost breakdown; the generator now
uses the 6-mer, which is why its constraint-rejection rate is 98.15% rather than the ~96.8%
the table implies.

**⚠️ BsaI is newly added here.** §4.10 shows BsaI was never cleared from any construct, so a
barcode carrying a BsaI site is a live hazard if a BsaI toolkit is ever used. It costs 0.5% of
survivors, so there is no reason not to exclude it now.

**On 6-mer vs 7-mer seeds.** Slutskin's Methods define the seed as "the MRE bases
complementary to bases 2–7 of the miRNA" — a 6-mer (verified, PMC5802814); we use the 7mer-m8
(bases 2–8). The 6-mer is strictly more aggressive and costs 35%
more of the space. Slutskin's own mitigation is instructive and matches §4.4's earlier
conclusion: they filtered only **highly expressed** miRNAs rather than all of them. That is the
right answer for us too, and it is blocked on the same thing — §5.11's unrecorded host line.

#### Two conservatism gaps closed after checking a third method

A published UMI generator — 10 nt, min **Hamming 4**, max homopolymer 2, GC 35–65%, common
restriction sites and ATG codons removed — exposed two places where we were the *looser* party.

**1. Our 7mer-m8 seed was less conservative than Slutskin's 6-mer.** Every 7mer-m8 target
*contains* its 6-mer core (219/219 verified), so the 6-mer filter rejects a strict **superset**.
We had been treating our choice as merely "different"; it was weaker. **Now switched to the
6-mer** (207 sites), matching the closest published precedent. We have ~10⁶× headroom at 20 nt,
so the extra ~35% cost buys strictly more caution for nothing.

**2. Levenshtein 3 alone permitted Hamming-3 pairs.** The metrics are not nested usefully: for
equal-length strings Levenshtein ≤ Hamming, so a **Hamming-4 code gives no indel guarantee at
all** (which is why FREE exists), while our **Levenshtein-3 code permitted Hamming-3 pairs** —
an exact census of the first build found **71 such pairs**. **Now enforcing both**, which is
stronger than either reference. Cost: 74 extra rejections out of 1.66 M draws.

| axis | that method | ours now |
|---|---|---|
| substitution separation | Hamming ≥4 | **Hamming ≥4** ✅ |
| indel separation | none | **Levenshtein ≥3** ✅ |
| GC window | 35–65% | **40–60%** (tighter) |
| homopolymer | max run 2 | max run 2 |
| ATG | "ATG codons" | **any position** (tighter) |
| restriction sites | "common" | 3 Type IIS named + all palindromes (below) |
| miRNA seed / polyA / hairpin | none | **all three** |

**Palindromic restriction sites are already excluded, for free.** Any palindromic 6-cutter `P`
satisfies `revcomp(P[0:3]) == P[3:6]` — its halves are non-overlapping reverse complements,
exactly what the hairpin rule rejects. A survey of 25–30 common 6-cutters found **zero** in the
set before any was named. So the self-complementarity rule, the most expensive filter in the
stack at 81% of survivors, is quietly buying exhaustive palindromic-site exclusion — better
than a hand-written blacklist. Only **non-palindromic** sites (our Type IIS enzymes) need
naming, which is why they are listed explicitly and nothing else needs to be.

#### Seed source: TargetScan vs miRBase

Adopted: **TargetScanHuman 8.0 `miR_Family_Info.txt`, human (Species ID 9606)**, reduced to
DNA 7mer-m8 target sites in [`mirna_seed_sites.tsv`](mirna_seed_sites.tsv) (2,064 sites,
tier-labelled). Reasoning:

| | TargetScan | miRBase |
|---|---|---|
| unit of data | **seed *families*** — exactly what we need to exclude | individual mature miRNAs; you derive seeds yourself from nt 2–8 |
| redundancy | already collapsed: many miRNAs share one seed | ~2,650 human matures collapse to far fewer distinct seeds — using the raw list **inflates the blacklist** |
| functional grading | ✅ **conservation tiers** (broadly conserved / conserved / poorly conserved), so you can exclude the ~108 families that actually matter | none; only a separate "high confidence" subset (~1,100) |
| site definitions | built around 7mer-m8 / 7mer-A1 / 8mer — the sites with *measured* repression | sequences only; site grammar is up to you |
| curation risk | a derived resource with its own inclusion calls | primary registry; includes many weakly-supported entries |
| practical cost here | **219 sites at `conserved+` → ~16% of barcode space** | ~2,064 sites → **~72% of barcode space, mostly for miRNAs not expressed in any given line** |

**Chosen default: `conserved+` (broadly conserved + conserved, 219 sites).** It captures the
families with real measured repression at ~16% of the barcode space. `all_human` is
available but over-filters on the strength of ~1,625 unassigned-tier entries.

⚠️ **Conservation is a proxy for the thing we actually care about, which is *expression in
our host line*.** A poorly-conserved miRNA can still be abundant in HEK293T. The correct
filter is "miRNAs expressed in the line we screen in" — which we cannot build until §5.11's
**unrecorded host line** is settled. Until then `conserved+` is the defensible default, and
this is a real dependency, not a nitpick.

⚠️ **Attribution:** `mirna_seed_sites.tsv` is derived from TargetScan (Agarwal / Bartel).
Confirm redistribution terms before publishing it outside the project.

#### ⚠️ Two further defects found in the shipped generator

1. **The `ATG` check is position-0 only.** It is `bc.startswith("ATG")`, so **~26% of
   barcodes that pass still contain an internal ATG**, despite the README having described
   this as "start-ATG avoidance." Consequence is mild here — the barcode is downstream of the
   ORF stop and IRES2 supplies its own start — but code and prose disagreed.
   `barcode_design.py` forbids ATG anywhere.
2. **There was no distance constraint at all.** Barcodes were independent random draws, so
   while collisions are vanishingly rare, **nothing guaranteed any error-correcting
   separation** — a single sequencing substitution could silently convert one valid barcode
   into another. `barcode_design.py` builds sets with a guaranteed **minimum Levenshtein
   distance** instead.
3. **Stale reference:** the docstring claimed these mirror "the production generator in
   `deliverome_analysis.codon_opt` (PR #10)." That module no longer exists — PR #10 renamed
   it to `seqtools` (commit `437aa34`), and neither is present. The local implementation is
   currently the *only* generator.

#### ✅ Shortening the barcode to free bases for a UMI — yes, 3–4 nt is affordable

`barcode_design.py` builds barcode sets with a **guaranteed minimum Levenshtein distance**
(the shipped generator had none). Achieved counts are **lower bounds** — the deletion-
neighbourhood filter over-rejects — so real capacity is at least this.

Target **20,000**, seed filter at `conserved+`:

| barcode nt | d=2 | **d=3** (corrects 1 edit) | d=4 | UMI bases freed | UMI codes |
|---|---|---|---|---|---|
| 14 | ✅ | 19,516 ✗ **at capacity** | 1,315 ✗ | 6 | 4,096 |
| 15 | ✅ | ✅ | 3,264 ✗ | 5 | 1,024 |
| **16** | ✅ | **✅** | 6,267 ✗ | **4** | **256** |
| **17** | ✅ | **✅** | 15,367 ✗ | **3** | **64** |
| 18 | ✅ | ✅ | ✅ | 2 | 16 |
| 20 (current) | ✅ | ✅ | ✅ | 0 | — |

**Also checked at the README's actual library size of 30,000** — 15, 16, 17 and 18 nt all
still clear it at d=3.

**How much margin?** The draw counts separate distance-driven rejection from
constraint-driven, which is what tells you whether a design point is comfortable or
saturated (at 30,000, expected draws with no distance rejection ≈ 81,700):

| barcode nt | draws for 30,000 | distance-driven | verdict |
|---|---|---|---|
| 15 | 244,094 | 2.99× | tight |
| **16** | 146,873 | **1.80×** | **workable margin** |
| **17** | 90,232 | **1.10×** | **distance barely binding — comfortable** |
| 14 (at 20k) | hit the 1.2M cap | 22.6× | **saturated — do not use** |

**Recommendation: 16 nt barcode + 4 nt UMI at min Levenshtein 3.** It clears 30,000 with
real margin and 256 UMIs is a genuinely useful dedup depth. Take **17 nt + 3 nt UMI** if you
want maximum slack and 64 UMIs is enough.

**⚠️ The hard trade: min Levenshtein 4 and a UMI are mutually exclusive.** d=4 needs 18 nt to
reach 20,000, which leaves 2 nt = **16 UMIs** — useless. So it is *either* single-edit
correction plus a real UMI, *or* distance 4 and no UMI. d=3 is the meaningful floor anyway,
since d=2 only **detects** one edit without correcting it.

#### Library scale, stated once: up to 30,000 ORFs at ≤500× representation

Earlier text in this doc read 30,000 as a *clone* count over 2,243 genes (≈13.4 clones per
gene). **That was wrong.** The library ceiling is **up to 30,000 distinct ORFs at 500× or
lower representation**, i.e. up to **15,000,000 transduced cells**. Consequences:

- **§5.4's "redundancy is the real buffer" holds, but for a different reason** — 500 cells
  per ORF, not 13. Synthesis error thins representation rather than losing ORFs, more
  strongly than stated.
- **The synthesis cost model is scoped to 2,243 ORFs and does not cover this ceiling.** At
  30,000 ORFs the Synplexity flat rate would be ~$1.5M against the **$221,446** modelled.
  Confirm whether 30,000 is a design ceiling or an order intent before treating §5.4's cost
  as the plan.
- **Barcode capacity is unaffected.** A 20-nt barcode admits ≥6.2 × 10⁸ codewords at min
  Levenshtein 3 (Gilbert-Varshamov), ~41× the 15M cells, so one 20-nt barcode can tag every
  cell clone unaided.
- **Clone-tracking arithmetic changes.** Cells sharing one barcode = 500 ÷ (barcodes per
  ORF in the plasmid pool). At 1 barcode/ORF, 500 cells are indistinguishable by *anything*
  encoded in the vector; at ~500 barcodes/ORF each cell is uniquely tagged. **Plasmid-pool
  complexity, not field width, decides whether expansion is measurable.** Sized in
  [S1E6](https://app.notion.com/p/3baab86c768381928004fea00ff5badc).

**⚠️ 3–4 nt cannot encode ORF identity.** Indexing 2,243 ORFs needs **≥6 nt raw** (4⁶ =
4,096), and more with any error correction — so if the intent was "read the ORF straight off
the barcode without a lookup table," that is a 6–9 nt budget, not 3–4. The freed bases can
only serve as a **molecular UMI** (dedup / counting independent integrations).

#### ⭐ Better use of the surplus: a SPLIT barcode that encodes ORF identity

Rather than spend freed bases on a UMI, spend them making the barcode *mean* something.
**Block-allocating barcode space to ORFs costs zero extra sequence space** — you need
2,243 × (clones per ORF) codewords either way; whether the assignment is structured or
arbitrary is just a convention. So the question is only whether synthesis can control which
barcode meets which ORF.

**The design that threads that constraint — split the barcode into two fields:**

```
… ORF … TAA[A,G] → [ORF-ID prefix, 13 nt, min Levenshtein 3, synthesized WITH the ORF]
                 → [clone suffix, 7 nt, degenerate, ligated from a pool] → …
```

| field | length | source | job |
|---|---|---|---|
| **ORF-ID prefix** | **13 nt** (min L3) | **rides along on each ORF's own synthesis order** — ~13 bp × 2,243 ≈ 29 kb total, negligible cost, no arrayed assembly | which ORF, error-corrected |
| **clone suffix** | **7 nt** degenerate | one bulk pool ligation after pooling | which clone |

Total **20 nt — exactly the current budget.**

**Prefix capacity (min Levenshtein 3, all our other constraints applied):**

| prefix nt | codewords built | covers 2,243? | margin |
|---|---|---|---|
| 11 | 1,151 | ✗ | — |
| **12** | 2,243 | ✅ | ⚠️ 12.6× distance-driven rejection — **at capacity** |
| **13** | 2,243 | ✅ | 3.4× — **headroom; use this** |

**12 nt is the floor and it is tight.** Use **13 nt**: §5.7 contemplates adding 632 isoforms
and 644 + 1,350 orthologs, which would take the set past 5,000 and 12 nt would not stretch.

**Clone-suffix width** — intra-ORF collision among ~13 clones of one ORF: 5 nt → 7.3% ·
**6 nt → 1.9%** · **7 nt → 0.47%** · 8 nt → 0.12%. 7 nt is the sweet spot inside a 20-nt
budget.

**⭐ The real prize is not convenience — it is a ~9× cut in wrong-ORF assignment.** §5.1's
~9% mispairing is a function of the **span between the variable region and the barcode**. Move
ORF identity to a prefix immediately 3′ of the stop and that span collapses:

| ORF-ID position | span | mispairing (5.3 × 10⁻⁴/bp) |
|---|---|---|
| **immediately 3′ of the ORF stop** | **~20 bp** | **1.05%** |
| prefix + short linker | ~60 bp | 3.13% |
| current barcode position | 177 bp | **8.96%** |
| barcode 3′ of WPRE (REPROcode-style) | 297 bp | 14.57% |

**This is the strongest argument for the split design**, and it attacks the problem §5.1
declared "not design-fixable" — correctly, for a *single distant* barcode, but not for ORF
identity specifically. ⚠️ The 20 bp figure assumes ORFs differ up to their 3′ end, so the
effective span is from the last ORF-distinguishing base to the prefix; verify per-ORF for
close paralogs and isoform pairs, which may be identical near the 3′ end.

**⚠️ Do not treat the designed prefix as ground truth.** Even at 1.05%, a mispaired clone
presents a *valid* prefix for the wrong ORF — it asserts a wrong answer confidently. The
empirical long-read variant↔barcode map (§5.1) is **still required**. What the prefix buys:
the error rate drops ~9×, a short-read amplicon alone reveals the *intended* ORF, and any
read whose prefix disagrees with the long-read-observed ORF is a **detected** chimera rather
than silent corruption.

**⚠️ Synthesis dependency — the "if my synthesis method allows" part.** This works because the
prefix is part of each ORF's own synthesis order, which per-gene clonal synthesis
(Synplexity, §5.4) supports directly. It does **not** work if the whole barcode arrives as a
degenerate oligo pool ligated onto already-pooled ORFs — then nothing controls which barcode
meets which ORF. The split is what makes it feasible: the *controlled* field is synthesized
per ORF, the *uncontrolled* field stays a pool and only ever confuses clone identity.

**Amplicon architecture** — the flanking primer sites are mix-and-matched:

```
[5′ fwd: Tycko et al. Cell 2020 (HT-recruit)] → 20-nt barcode → [3′: lentiGuide-BC (Feldman Cell 2019)]
```

Why not lentiGuide-BC's own 5′ site: (1) the sequence immediately 5′ of its barcode is
the **lower-Tm padlock probe site**, a poor PCR primer; and (2) its **NGS 5′ primer site
lands inside PuroR — i.e. coding sequence**, which shouldn't be used as an amplification
handle (couples amplicon design to a functional element and constrains future recoding).

**Placement vs WPRE** is still open — for 3′ single-cell capture the barcode wants to sit
3′ of WPRE near the polyA; for bulk amplicon it's flexible. Note this trades against
§5.1 (shorter variable-region→barcode span is better). **REPROcode chose 3′-of-WPRE and a
*transcribed* barcode** (§4.3), i.e. it optimised for capture at the cost of the longest
possible span — the opposite of our lean. That's the clearest precedent on either side of
this decision.

### 4.5 Stop codons — `TAA[A,G]`

*(Interacts with §4.11: a Type IIS shuttle site placed immediately 3′ of the stop is
compatible with this policy — and for BsaI the site's own first base supplies the +4 purine,
so `TAAGGTCTC…` costs no spacer at all.)*

**Policy: `TAA` followed by a purine** (`TAAA` or `TAAG`; `TAAR` in IUPAC). Write it as
**`TAA[A,G]`** in designs.

The stop signal is effectively **4 bases**: the +4 nucleotide is part of what release
factor reads and matters as much as the codon.

- **McCaughan, Brown, Dalphin, Berry & Tate 1995, *PNAS*** measured all 12 four-base
  signals: **+4 order `A ≈ G >> C ≈ U`**, with a **70-fold** efficiency range across
  UAA·N vs only ~8-fold for UGA·N and UAG·N. Strong signals are the ones most frequent at
  natural termination sites. → **`UAA` has the highest achievable termination strength**,
  but only with a favourable +4 base.
- **Tate 1995** (review) — the +4 base determines RF decoding efficiency in vivo; **weak
  signals cluster at readthrough/recoding sites, strong ones in highly expressed genes.**
- **Why avoid `TGA`:** UGA readthrough is decoded by near-cognate **Trp/Cys/Arg** tRNAs
  (vs Gln/Tyr/Lys at UAA/UAG — Blanchet 2014), and UGA is the canonical leaky stop.

**⚠️ Both stops in the current construct are `TGA`.** The ORF stop is the consequential
one: readthrough there extends the protein past the **ALFA tag** into the UTR/barcode
region, and since ALFA sits immediately before that stop, an appended extension risks
**occluding NbALFA access** — degrading the detection reagent's own target. It also
increases translation of the WPRE's X ORFs (§5.3). Tandem stops (`TAA TAA`) are the
belt-and-suspenders option.

### 4.6 P2A over T2A

**Kim et al. 2011, *PLoS ONE*** compared the four common 2A peptides (F2A/FMDV,
E2A/ERAV, **T2A**/*Thosea asigna*, **P2A**/porcine teschovirus-1) across three human cell
lines, zebrafish embryos and adult mice: **P2A had the highest cleavage efficiency in
every context examined.** Incomplete cleavage would leave an mScarlet3–PuroR **fusion**,
distorting the fluorescent readout and coupling reporter to selection marker. P2A also
matches Fragmid's `P2A_v1.1`.

### 4.7 ALFA v1.5 — single copy, proline-flanked

`pF_AA197` (Addgene #216049 — a Fragmid C-terminus fragment: ALFA + NpNLS + 3xHA)
supplied **5 synonymous codon variants** of the ALFA core `SRLEEELRRRLTE`, labeled
v1.1–v1.5. **We use a single copy of v1.5**, not a tandem array.

- **Why v1.5 over v1.1:** v1.1 uses the **same codon (`AGG`) for all 4 arginines**,
  creating a 9-nt run `AGGAGGAGG` across the R-R-R stretch — exactly the kind of local
  repeat that invites recombination/deletion. **v1.5 uses a distinct codon at all four
  arginine positions** (`CGC/AGG/CGG/AGA`), breaking that run. Deciding factor.
- **Proline flanks are present.** The construct's `…GGGGSGGGGS + P` / `P + GGGGS…`
  context gives **`P-SRLEEELRRRLTE-P`**. This matters: **Vicidomini et al. 2024** state
  the ALFA tag "is flanked by prolines to isolate its alpha-helical secondary structure,
  which is necessary for the recognition by the NbALFA." (Note the v1.x variants as
  supplied by `pF_AA197` are bare cores — the flanks come from our surrounding context.)

### 4.8 miRNA control system — not yet integrated

Placement depends on modality: a **Pol III shRNA/amiR cassette** is an independent
transcription unit (put it downstream of WPRE or in the backbone); **3′UTR miRNA target
sites** for de-targeting are fine alongside WPRE; an **intronic miRNA** must go in an
intron — never tacked after WPRE, where Drosha cleavage would sever the polyA. **Open.**

### 4.9 Domestication — recoded off pRDA_722

The build is recoded off **pRDA_722 (#216077)** to clear BbsI + Esp3I. **Verified
directly:**

| | **pRDA_722 (#216077)** | **#236083 (old parent)** |
|---|---|---|
| Total BbsI/Esp3I | **2** | **5** |
| BbsI | 2 — pos **2617**, **3042** | 2 — pos 6215, 6372 |
| Esp3I/BsmBI | **0** | 3 — pos 252, 306, 5445 |
| WPRE region | ✅ **CLEAN** (3057..3645) | ⚠️ 1 × BbsI (6215) |

pRDA_722's two BbsI sites **are the assembly junctions** — at 2617 and 3042, between the
cPPT/CTS (2440..2557) and the WPRE (3057..3645), exactly where the cassette drops in.
That's the correct signature for a Golden Gate destination vector: enzyme sites **only**
at the junctions, zero elsewhere.

**This also avoided two edits that would have been unattractive in #236083:** its BbsI at
**6372 sits inside the 3′LTR (ΔU3)** — an essential viral element, risky to hand-mutate —
and its Esp3I at **5445 sits inside its own `pac`/PuroR CDS**.

Other scanned elements: **Kozak block ✅ clean.**

#### The two hand-domesticated Esp3I edits in the #236083 backbone

Before pRDA_722 was adopted, **#236083's Esp3I sites 252 and 306 were knocked out by hand**
— these are the p003/p004 lineage (§4.10). Recording the exact edits, since they were not
written down anywhere and are not recoverable from the sequence alone:

| # | site | span | strand | edit | before → after |
|---|---|---|---|---|---|
| 1 | Esp3I **252** | 252..257 | + | **253 `g→A`** | `CGTCTC` → `CATCTC` |
| 2 | Esp3I **306** | 306..311 | − (`CGTCTC` on the bottom strand) | **310 `c→A`** | `GAGACG` → `GAGAAG` |

Both edits change **base 2 of the `CGTCTC` recognition sequence** on their respective
strands — the same position each time, so the pair is a consistent rule, not two ad-hoc
fixes. In context (edited bases capitalized):

```
#236083 242..313, before:  aagctgtgaccgtctccgggagctgcatgtgtcagaggttttcaccgtcatcaccgaaacgcgcgagacgaa
#236083 242..313, after:   aagctgtgaccAtctccgggagctgcatgtgtcagaggttttcaccgtcatcaccgaaacgcgcgagaAgaa
```

**Why these two were safe to hand-edit** — verified by `enzyme_scan.py`, not assumed:

- Both sit in **unannotated backbone**, between the pRS primer-bind region and the
  **`bla`/AmpR promoter (403..507)**. No annotated feature spans either site.
- **`pGEX` reverse primer site (257..279) intact** — edit 1 is 4 bp upstream of it.
- **`pBR322` forward primer site (317..335) intact** — edit 2 is 7 bp upstream of it.
- **`bla` promoter untouched** (starts 90 bp downstream of edit 2).
- **No new site created** by either edit — rescanned for BbsI, Esp3I/BsmBI, BsaI, AarI and
  SapI on both strands: zero hits in the edited window.

**The third site could not be fixed this way**, which is the whole argument for pRDA_722:
**Esp3I 5445 is inside `pac`/PuroR (5059..5658)**, so it needs a *synonymous* recode rather
than an arbitrary base change — and #236083's BbsI **6372 is inside the 3′LTR (ΔU3)**.

### 4.10 The construct series — p003/p004 (domesticated #236083) and p006/p007 (not)

There are **two distinct #236083-derived lineages**, and the difference is exactly whether
the domestication edits above were applied. Rescanned from the GenBanks by
`enzyme_scan.py` — both strands, treating each plasmid as circular:

| construct | bp | **BbsI** | **Esp3I** | BsaI | AarI | SapI |
|---|---|---|---|---|---|---|
| #236083 parent | 7,119 | 2 | **3** (252, 306, 5445) | 4 | 0 | 2 |
| **p004** wpre-only | 6,201 | **0** ✅ | **0** ✅ | 3 | 0 | 2 |
| **p003** assembled | 10,759 | **0** ✅ | **0** ✅ | 4 | 2 | 3 |
| **p003** cmv-empty-puro variant | 10,757 | **0** ✅ | **0** ✅ | 4 | 2 | 3 |
| p006 ctermtag | 10,027 | 4 | 3 (252, 306, **8353**) | 4 | 1 | 3 |
| p007 V291 | 10,027 | 4 | 3 (252, 306, **8353**) | 4 | 1 | 3 |
| pRDA_722 / #216077 | 7,965 | 2 (junctions) | **0** ✅ | 3 | 0 | 2 |

- **p004** — domesticated backbone, **WPRE only** (PacI/ClaI site).
- **p003** — domesticated, **fully assembled** EFS–TFRC–BC–ALFA–IRES(Fragmid)–WPRE–Puro.
  Two variants exist, 10,757 and 10,759 bp, differing by 2 bp.
- **p006 / p007** — **deliberately not** domesticated (below).

**p003/p004 reached BbsI = 0 and Esp3I = 0** — so the `pac`/PuroR site (5445 in the parent)
was also resolved in that lineage, presumably by taking Fragmid's PuroR rather than
point-editing the CDS.

⚠️ **"Domesticated" here means BbsI + Esp3I only — not Type IIS-free.** Worth stating
plainly so nobody over-reads it later:

- **BsaI was never cleared anywhere.** 1223 and 2854 persist in *every* construct including
  p003/p004, and pRDA_722 carries 3 of its own (836, 3771, 7112). Harmless — our documented
  pair is BbsI + Esp3I/BsmBI, and **BsaI (`GGTCTC`) is a different enzyme** from Esp3I
  (`CGTCTC`) — but it does mean a future switch to a BsaI-based toolkit needs fresh work.
- **SapI is present in all of them** (2–3 sites), including pRDA_722.
- **AarI appears only in the insert-bearing constructs** — 2 in p003, 1 in p006/p007, 0 in
  both parents and in p004. Consistent with §5.2's note that the Fragmid IRES2 carries
  1 AarI; the insert is where they come from, not the backbone.

**p006 and p007** are the fast-iteration pair: two quick constructs on a **CMV promoter,
non-domesticated background** (parental backbone #236083, not scrubbed of BbsI/Esp3I —
skips domestication for faster startup):

- **p006** — **C-terminal** tagged (SpyTag003/ALFA, as in §2).
- **p007** — **internally tagged at V291**, per the internal ALFA-tag insertion strategy
  in **Lenaerts et al., *EndoNB***, bioRxiv
  [10.1101/2025.06.08.658482](https://doi.org/10.1101/2025.06.08.658482), which ALFA-tags
  TFRC internally (extracellular loop insertion rather than N/C-terminal) for an
  endocytosis assay. ⚠️ **V291 not independently confirmed** against the paper's
  supplementary table — the exact residue wasn't recoverable from the rendered full text;
  taken as user-supplied pending a direct check.

Because these skip Golden Gate, their BbsI/Esp3I sites are **harmless and expected** —
and that makes them the ideal test bed for the IRES A/B in §5.2.

### 4.11 Planned: Type IIS shuttle sites flanking the insert — enzyme still TBD

*(Logged as an open decision in **§5.14** — this section holds the analysis.)*

**Goal:** flank the insert with a third Type IIS enzyme's sites so the whole barcoded ORF
collection can be **moved to another backbone in one Golden Gate reaction**, making these
vectors reusable by others. Precedent is **Tedman et al. 2025** (§References), who state it
explicitly:

> "we added BsmBI cleavage sites upstream of the UMI and at the immediate 3′ of the ORF to
> facilitate the transfer of this library from one plasmid backbone to another with a single
> golden gate assembly reaction"

— giving a library "backward-compatible with several current activity assays and
forward-compatible with future assays."

**⚠️ Correction to the framing: they used BsmBI (= Esp3I), not BsaI.** That matters,
because for *us* those two enzymes have opposite properties — see below.

**Why this is worth real effort for us specifically.** §5.1 requires an empirical
**long-read variant↔barcode map** because ~9% of clones are mispaired and can't be trusted
from the design. That map is the single most expensive QC artifact in the project. If the
barcode travels *with* the ORF through a shuttle reaction — which is exactly why Tedman put
the site **upstream of the UMI** rather than between UMI and ORF — then the map is built
**once** and amortised across every future backbone. Without shuttle sites, each new
backbone means re-deriving it.

**Which enzyme? The §4.10 scan already answers this, and it is not BsaI.**

| candidate | sites already present | verdict |
|---|---|---|
| **BbsI** | 2 (the assembly junctions) | ❌ already our assembly enzyme — can't also be the shuttle |
| **Esp3I/BsmBI** | **0** in p003, p004 *and* pRDA_722 ✅ | **Cleanest by far** — but §1 **reserves it for downstream guide cloning.** Also what Tedman used. The reservation is the only objection. |
| **AarI/PaqCI** | 0 in both parents and p004; 1–2 in insert-bearing constructs (from the Fragmid IRES2) | **Best unreserved option** — only the IRES2's single site to remove |
| **BsaI** | **2–4 everywhere** (1223, 2854 in all #236083-derived; pRDA_722 has 836, 3771, 7112) | ⚠️ **Worst option.** A BsaI shuttle would cut the backbone. Needs 2–3 knockouts first. |

So: **BsaI is the one enzyme we'd have to domesticate *into* cleanliness first**, precisely
because §4.10 shows it was never cleared. The decision is really *"is the guide-cloning
reservation on Esp3I worth more than a clean shuttle?"* — if we can move guide cloning to
another enzyme, Esp3I is free, zero extra domestication, and matches the published
precedent. Otherwise AarI/PaqCI at the cost of one IRES2 edit.

**✅ Compatible with the §4.5 stop policy either way.** A shuttle site can sit immediately
3′ of the stop without breaking `TAA[A,G]`, and for BsaI the site's own first base *is* the
purine, so no spacer is needed at all:

```
BsaI         TAA G GTCTC…   →  +4 = G ✅, site's first base doubles as the +4 purine
Esp3I/BsmBI  TAA G CGTCTC…  →  +4 = G ✅, costs one spacer base
AarI/PaqCI   TAA G CACCTGC… →  +4 = G ✅, costs one spacer base
```

⚠️ **Do not copy Tedman's stop handling verbatim.** They "replaced the native stop codon of
the coding sequence with a glycine codon, a stop codon, and an additional BsmBI cleavage
site" — i.e. **Gly + stop**. Our ORF's 3′ end already carries the C-terminal tag cassette
and a `TAA[A,G]` stop (§4.5); an extra glycine is a scar we don't need, and their `TGA`-era
assumptions don't apply since we've deliberately chosen `TAA`. Take the site placement, not
the stop recoding.

---

## 5. Risks and open questions

### 5.1 Recombination — ORF↔barcode decoupling (~9%, not design-fixable)

Lentivirus is pseudodiploid; reverse transcription template-switches between the two
co-packaged genomes, which can **mispair an ORF with another clone's barcode**. This
applies to whatever is *variable and must stay paired with the barcode* — internal or
C-terminal tags through to the barcode — but **not** the N-terminal region, which is
upstream of the informative element and invariant across the library.

**Walton, Qin & Blainey (CROPseq-multi)** model the rate directly:

> "lentiviral recombination frequency could be modeled by a fixed per-bp probability of
> template switching of **5.3×10⁻⁴ per bp**, within the range of previously reported
> measurements of **4.9–13.5×10⁻⁴ per bp**"

Their model is fit **as a function of the length of *homologous* intervening sequence** —
exactly our case: the linker/tag scaffold between the variable region and the barcode is
**constant across library members**, hence 100% homologous between two co-packaged
genomes, hence fully available as switch substrate.

`P(≥1 switch over L bp) = 1 − (1 − p)^L`. **The 5% threshold sits at ~97 bp**
(38–105 bp across the rate range).

| Span measured | bp | @4.9×10⁻⁴ | **@5.3×10⁻⁴** | @13.5×10⁻⁴ | of 30,000 clones |
|---|---|---|---|---|---|
| stop → barcode (**naive**) | **21** | 1.02% | **1.11%** | 2.80% | ~332 |
| linker start → barcode (**correct**) | **177** | 8.31% | **8.96%** | 21.27% | **~2,687** |
| + snorkel +60 bp | 237 | 10.97% | 11.81% | 27.40% | ~3,542 |
| + snorkel +120 bp | 297 | 13.55% | 14.57% | 33.05% | ~4,370 |

Figure with the rate range as error bars: `figures/s1e2_recombination_risk.svg`
(script `recombination_analysis.py`).

**The choice of measurement flips the verdict.** 21 bp → ~1.1%, comfortably under
threshold. 177 bp → ~9.0%, above it, and up to 21% at the high end. **Method (a)
underestimates by ~8.1×.** **177 bp is the correct measure**: the 21 bp figure implicitly
assumes everything upstream of the stop is invariant, but that upstream region is exactly
what varies. **A snorkel makes it worse** (12–15% at the point estimate).

**⚠️ Two distinct failure modes — don't conflate them; they have different fixes:**

**(1) Barcode↔variant uncoupling** — *inter*-genome switching in a heterozygous virion.
Driven by the **length of the constant homologous span**. Levers that work:
- **Shorten the span** — the only strong design-level lever.
- **Low-MOI packaging** — fewer heterozygous virions.
- **Long-read QC** — build the *empirical* variant↔barcode table; don't trust the design.

  **Codon-diversifying a *constant* linker does NOT reduce this risk** — every library
  member carries the same linker, so inter-genome homology stays 100%. (Giving each member
  a *different* synonymous linker would help, but is impractical for a shared scaffold.)

**(2) Internal deletion / repeat instability** — *intra*molecular recombination between
nearby direct repeats, during bacterial propagation or RT. Driven by **local repeats**.
Our scan found **9 exact 12-mer repeats**, all inside the **GS-linkers flanking SpyTag003
and ALFA** (Gly/Ser codons reuse a tiny codon set, so GS-linkers are intrinsically
repetitive at the DNA level). Same failure mode as the ALFA v1.1 `AGGAGGAGG` run.
- **Lever: codon-diversify the GS-linkers.** ← this is what the linker-diversification
  work item fixes, and it *is* important — for reason (2), not reason (1).

### 5.2 The domesticated IRES2 is unvalidated

Fragmid's `IRES2_v1` differs from the EMMA/#100673 version by **3 bp — all restriction-site
removals**, i.e. **domestication, not a sequence disagreement**:

| # | Position | Ours (#100673) | Fragmid `IRES2_v1` | Site destroyed |
|---|---|---|---|---|
| 1 | **9** | `T` | `G` | **BbsI** `GTCTTC` (at pos 6) |
| 2 | **108** | `A` | `C` | **BbsI** `GAAGAC` (at pos 106) |
| 3 | **352** | `G` | `A` | **ApaLI** `GTGCAC` (at pos 350) |

Full scan: ours has **2 × BbsI** + 1 ApaLI + 1 AarI; Fragmid's has **0 × BbsI**, 0 ApaLI,
1 AarI. Both free of BsaI/Esp3I/SapI. **Each toolkit domesticates for its own enzyme** —
EMMA doesn't use BbsI, Fragmid does. **Adopted `pF_AA333`**, which cleared a live bug: our
original #100673-derived IRES carried two BbsI sites that Golden Gate would have cut.

**But treat all three mutations as unvalidated:**

- **No functional data in the paper.** Searching McGee et al.: **"IRES" appears twice**,
  both generic architecture ("separated from the Cas protein by the use of 2A sites or
  internal ribosome entry site (IRES) elements"; "…fragments that contain IRES elements
  … and thus are **preceded by a stop codon**"), and **"IRES2" appears zero times.** No
  expression data, no selection-efficiency comparison, no characterisation.
- **No independent adoption.** Per an Addgene sequence search, **Fragmid is the only
  vector on Addgene carrying this particular domesticated IRES sequence** — so no other
  depositor's data can vouch for the 3 mutations either.
- **The riskiest change is pos 352** — only **6 nt from the 3′ end**, in the IRES 3′
  boundary / start-codon-proximal region, the most function-sensitive part of an
  EMCV-family IRES. (#1 and #2 sit ~349 and ~250 nt from the 3′ end.)

*(Fragmid cites PMID 31504667 = Yang & Wang 2019 for its IRES, but that is a **review
article**, not a sequence-of-record — a citation-of-convenience. Irrelevant to the
decision, since the 3 bases were changed on purpose, but the lineage is better traced
through the EMMA part #100673.)*

#### How to compare the two IRES variants — is fluorescence intensity enough?

**Short answer: raw fluorescence intensity is NOT enough — but *ratiometric* fluorescence
is, for a go/no-go.**

**Why absolute MFI fails.** In `ORF → IRES → mScarlet3-P2A-PuroR`, per-cell mScarlet3
signal is a product of IRES efficiency **and** nuisance variables that differ between
separately-transfected samples: transfection/transduction efficiency, copy number and
integration site, prep quality, cell state, expression duration, and day-to-day cytometer
drift. Comparing bare MFI **conflates IRES efficiency with delivery efficiency** — the
classic failure mode of this experiment.

**The fix: measure the IRES-dependent product *relative to* the cap-dependent product
from the same transcript, per cell.**

- **Option A — ratiometric flow with tags we already have (recommended first pass).**
  Cap-dependent = **TfR1 with the extracellular C-terminal ALFA tag**, live-stained with
  fluorescent **NbALFA** (no perm). IRES-dependent = **mScarlet3**. Read the **per-cell
  mScarlet3 ÷ ALFA ratio** and compare between variants. Valid because any inefficiency
  in NbALFA staining or TfR1 trafficking is **identical across the two constructs** (they
  differ only in the IRES), so it cancels. Requires same-day acquisition, matched
  settings, unstained + untransduced controls, **n ≥ 3** biological replicates,
  independent preps.
- **Option B — dual-luciferase** `cap-Renilla – IRES – Firefly`, report **Fluc/Rluc**.
  Wide dynamic range, ratiometric by construction. **Essential controls: promoterless and
  monocistronic versions** — apparent "IRES activity" that is really a **cryptic promoter
  or splicing artefact** is a well-documented pitfall; without these a positive result is
  uninterpretable.
- **Option C — puromycin dose–response** (cheap, functional, most decision-relevant). The
  IRES drives mScarlet3-P2A-PuroR, so IRES output sets puromycin resistance directly.
  Plate equal cell numbers, titrate puromycin, compare survival/IC50. Good orthogonal
  confirmation.

**Practical notes.** Run in the **p006/p007 non-domesticated background** — both variants
are clonable there. **Transient** is fastest but has large copy-number variance, so
ratiometric is mandatory; mScarlet3 was reported as strong in transient systems, so
transient + ratiometric is a reasonable first pass. Keep timepoints matched.
**If it fails, the fallback is architectural, not another IRES variant.** Move PuroR (and
optionally mScarlet3) onto **its own promoter** — MORF's `SV40 → PuroR` or RESOLUTE's
`PGK → rtTA-VP16–2A–PuroR` (§4.3). That removes the IRES from the selection path entirely.
The cost is losing linkage: puro resistance would no longer certify that the ORF cistron
was translated. Decide this *before* the A/B so the failure branch is pre-planned.

**Pre-specify the threshold:** within ~20% → equivalent, proceed with the domesticated
version; **≥2-fold reduction → a real problem**, and worth either re-domesticating with
different point mutations (avoiding the start-codon-proximal pos-352 change) or keeping
the undomesticated IRES and moving the assembly junction.

### 5.3 The WPRE is the native, non-X-inactivated element

The WPRE is a ~590 bp fragment of the **woodchuck hepatitis virus** genome that overlaps
the WHV **X protein (WHx)** gene, so the standard full-length element **carries an ORF
encoding a truncated WHx peptide**. WHx is a transactivator implicated in WHV-associated
hepatocellular carcinoma — hence the vector-safety concern.

Our 589 bp element (from pRDA_722) has ATG-initiated ORFs of **166, 162, 156, 140 aa** —
far exceeding **Schambach 2006**'s engineered-safe criterion of *no ORF > 25 aa*. And
pRDA_722's WPRE differs from #236083's by **exactly 1 bp** (pos 539 T→C = the BbsI
knockout), so **Fragmid domesticated BbsI but did not inactivate the X ORF.**

Two characterised fixes if we want it addressed:

- **`mut6`** — **Zanta-Boussif 2009, *Gene Therapy***: WPRE transcriptional activity was
  undetectable in a lentiviral vector **but "the element was capable of translating a
  polypeptide,"** and this **"was abrogated by mutating the WHX ORF translation start."**
  Critically, **"the native form or mutated derivatives functioned equivalently"** for
  expression — **so the swap is ~free.**
- **Schambach ORF-free PRE** — PRE mutants **"devoid of X protein open reading frames …
  as well as other ORFs exceeding 25 amino acids, without significant loss of RNA
  enhancement activity"**; the X promoter could also be deleted.

Both papers are careful that a causal oncogenic role for the PRE **"remained unproven"** —
this is precaution, not demonstrated hazard.

**Links to §4.5:** the WPRE sits in the **3′UTR**, so its ORFs translate only via
**readthrough or reinitiation** — which is exactly what Zanta-Boussif observed. Our leaky
`TGA` stops make WPRE-ORF translation *more* likely, so `TAA[A,G]` is a second,
independent reason for that change.

*Honest limit:* the native-vs-inactivated call is inferred from **ORF content and
length** — I could not match any ORF to a canonical WHx N-terminus. To settle it, align
our 589 bp against the WHV X reference and against `mut6`.

### 5.4 Synthesis error burden and cost — the 2,243-ORF library

**Gene set.** The **2,243 canonical** genes from accessible-surfaceome
**PR [#130](https://github.com/Deliverome-Project/accessible-surfaceome/pull/130)** —
obtained by running that PR's `passes_canonical` over the 5,130-gene deep-dive cohort in
the **read-only public D1 mirror** (`surface_annotation`). Count reconciles **exactly**
with the PR (1,452 → 2,243). Lengths are real `topology_public.protein_length` values,
**2,243/2,243 matched, 0 unmatched.**

| | aa | CDS bp |
|---|---|---|
| **mean** | 676.3 | **2,032** |
| median | 501 | 1,506 |
| p25 / p75 | — | 1,014 / 2,469 |
| p90 | — | 3,585 |
| **total** | — | **4,557,594 bp** |

#### Error model

Per-bp error rates as supplied: **Twist pools 1:3000** (≤1.8 kb), **Twist fragments
1:7500**, **Abstract Bio 1:20,000–1:40,000** (≤~1.2 kb; midpoint 1:30k used).

Substitution consequences use the **exact** single-nucleotide spectrum of the standard
code — enumerating all 549 changes across the 61 sense codons:

| | count | share |
|---|---|---|
| synonymous | 134 | **24.4%** |
| missense | 392 | **71.4%** |
| nonsense | 23 | **4.2%** |

→ **75.6% of substitutions are protein-altering.**

Two soft parameters, both flagged:
- **Indel fraction = 30%** (assumed). **Unverified for Twist** and the biggest lever on
  the severe-error column, since every indel frameshifts. *Ask Twist for the
  substitution/deletion/insertion breakdown.*
- **Fraction of missense that is severely damaging** — the scenario tables above use
  **32%** (AlphaMissense genome-wide). ✅ **Now measured on our own gene set: 0.427**,
  which raises the severe column by ~1.5 points (§5.5). The honest band is **10–16%**;
  see §5.5 for why this parameter is the *less* important of the two.

#### Cost and usable output

Figure: `figures/s1e2_synthesis_cost.svg` · table: `synthesis_cost_by_scenario.tsv` ·
script: `synthesis_cost_analysis.py`.

| Scenario | ORF only | + cassette | 2 variants/ORF | no protein-altering error | severe |
|---|---|---|---|---|---|
| **A** all fragments | $319,032 | $362,209 | $638,063 | 81.3% (1,824) | 11.8% |
| **B** pools ≤1.8k + frag | **$221,966** | $269,345 | $443,932 | 72.4% (1,624) | 17.6% |
| **C** AbstractBio ≤1.2k + frag | $276,335 | $338,153 | $552,669 | **83.8% (1,879)** | **10.3%** |
| **D** AbsBio ≤1.2k + pools 1.2–1.8k + frag | **$224,478** | $270,760 | $448,956 | 79.3% (1,778) | 13.3% |

**⚠️ Costs above are ORF-only — no tags, barcode, linkers or adapters.**

#### The dominant fact: 94.5% of cost is 870 long ORFs

Scenario D: **AbstractBio ≤1.2k = $7,116 (3.2%)** · **pools 1.2–1.8k = $5,281 (2.4%)** ·
**fragments >1.8k = 870 ORFs, 3.03 Mbp, $212,080 (94.5%)**.

So **the short-end vendor choice is ~3% of the budget** — it is a *quality* decision, not
a cost one. **D is the recommendation:** $2,512 more than B (+1.1%) for **+6.9 points of
usable library (~154 more ORFs)** ≈ **$16 per additional usable ORF**. C buys a further
4.5 points but costs **$52k more** (~$515/ORF) because its 1.2 kb cap pushes 549 extra
ORFs into fragments.

Length concentration (figure panel 3): **top 50 longest = 12.5% of all bp; top 100 =
19.1%; top 500 = 48.7%.**

#### ⚠️ The cassette is a $46k threshold trap — and our architecture already avoids it

Adding ~275 bp of cassette **to each synthesized ORF** costs **+$46,283 (+21%)** in
scenario D — mostly because **231 ORFs get bumped from cheap pools across the 1.8 kb
line into $0.07/bp fragments** (870 → 1,101 fragments), not because of the extra bases.

**In Fragmid/Golden Gate the tags, barcode and primer sites live in separate fragments
and the destination vector**, so the ORF is synthesized alone and assembled. Keep it that
way and this cost never appears. **Bundling the cassette into each ORF order would be a
~$46k mistake.**

#### Redundancy is the real buffer

30,000 barcoded clones ÷ 2,243 genes ≈ **13.4 clones per gene**. At scenario D's 79.3%
clean rate that's **~10.6 error-free clones per gene** — so synthesis error does not
*lose* genes, it just thins per-gene representation. Combined with the ~9% recombination
mispairing (§5.1), the library is well-buffered **provided** the long-read
variant↔barcode map is built so error-bearing clones can be *identified and excluded*
rather than silently scored.

#### ⭐ Vendor tiers beyond fragments — Twist clonal changes the recommendation

Scenarios A–D above all route long ORFs to **fragments ($0.07/bp, 1:7500)**, which is where
94.5% of the money and essentially all of the error burden sits. Two tiers exist above it,
and the metric that decides between them is **cost per *usable* construct**, not per bp:

| tier | $/bp | error rate | clean @ 2 kb | **effective $/bp per usable** |
|---|---|---|---|---|
| Twist fragments | $0.07 | 1:7,500 | 80.2% | $0.0873 |
| ⭐ **Twist clonal genes** | **$0.09** | **sequence-perfect** | **100%** | **$0.0900** |
| Elegen Infinia ⚠️ unverified | $0.14 | ~1:70,000 | 97.6% | $0.1434 |
| Elegen Infinia ⚠️ unverified | $0.20 | ~1:70,000 | 97.6% | $0.2049 |

**Clonal costs 1.29× fragments per base but is verified perfect**, so it wins outright
above the length where fragments drop below **q = 0.07/0.09 = 77.8% clean**: **2,273 bp**
(no protein-altering error) or **3,852 bp** (no severe error).

| strategy | cost | usable | severe |
|---|---|---|---|
| **D** — no clonal | $224,478 | 79.3% | 13.3% |
| F1 clonal >5 kb (103 giants) | $242,213 | 81.9% | 11.4% |
| F2 clonal >3.85 kb (187) | $249,460 | 83.3% | 10.5% |
| ⭐ **F3 clonal >2.27 kb (656)** | **$276,416** | **89.0%** | **7.0%** |
| F4 clonal all >1.8 kb (870) | $285,072 | 90.9% | 5.8% |

**F3 is the recommendation: +$51,938 over D (+23%) buys +9.7 points usable and nearly
halves the severe-error rate.** Two consequences:

- **⚠️ Clonal strictly dominates duplication — take 2× builds off the table.** For the
  870-ORF long tail: clonal **$272,675 @ 100% usable ($313/usable ORF)** vs duplicated
  fragments **$424,161 @ 89.5% ($545/usable ORF)**. Duplication was only ever a workaround
  for lacking a guaranteed-clean option; §5.8's duplication arithmetic is retained for the
  record but is no longer the plan.
- **Elegen is a fallback, not an upgrade.** It is cheaper per usable base than *fragments*
  above 6,269 bp ($0.14) or 9,496 bp ($0.20), but **clonal beats it at every length**. Its
  only role is **lengths clonal cannot reach** — so **Twist's clonal length ceiling is the
  gating unknown.** If clonal reaches 5 kb, F3 stands and Elegen isn't needed; capping at
  3.5 kb costs **+$72k** over F3 to cover the remainder with Elegen, for slightly *lower*
  usable output. Worth asking Elegen about turnaround, which may justify a premium
  independently of price. ⚠️ **Elegen's $0.14–0.20/bp and 1:70,000 are unverified.**

**MUC16 (43.5 kb) and MUC3B (40.4 kb) remain unorderable by any route** and still need
truncation or a short isoform. Full scenario tables, the giants audit and the
rejected-gene audit: [`REPORT_library_scoping.md`](REPORT_library_scoping.md) §5b–5c.

#### ⭐ Head-to-head — and Synplexity at $50/gene up to 3 kb changes the plan

Figure: `figures/s1e2_strategy_comparison.svg` · table: `strategy_comparison.tsv` · script:
`strategy_comparison.py`. Three families, each with one free parameter, so each traces a
cost-vs-quality frontier:

- **F-series** — scenario-D short end, then **Twist clonal** above a length threshold.
- **E-series** — ORFs >1.2 kb split into **AbstractBio pieces** and stitched. Parameter: the
  per-piece fee (**still unquoted**).
- **S-series** — ⭐ **Synplexity clonal genes, flat $50/gene, hard cap 3 kb.**

**What the 3 kb cap does.** It splits the library cleanly in two:

| | ORFs | share | bp | cost | quality |
|---|---|---|---|---|---|
| **≤3 kb — Synplexity covers it** | **1,897** | **84.6%** | 2,772,267 | **$94,850** flat | **sequence-perfect** |
| **>3 kb — must route elsewhere** | **346** | 15.4% | 1,785,327 | see below | the only remaining error source |

**So after Synplexity, 100% of the residual cost and error lives in 346 ORFs.** That tail
(mean 5,160 bp, max 43,524) costs **$124,973** as fragments, **$160,679** as Twist clonal, or
**$249,946** as Elegen.

| scenario | cost | usable | severe | $/usable |
|---|---|---|---|---|
| **D** old baseline, no clonal | $224,478 | 79.3% | 13.3% | $126 |
| **F3** *previous recommendation* | $276,416 | 89.0% | 7.0% | $138 |
| **E** stitch @ $44/piece | $222,896 | 94.7% | 3.2% | $105 |
| **S4** Synplexity only 1.8–3 kb + clonal tail | $199,277 | 90.9% | 5.8% | $98 |
| ⭐ **S7** AbstractBio ≤1.2 kb + **Synplexity 1.2–3 kb** + clonal tail | **$221,446** | **99.1%** | **0.5%** | **$100** |
| **S1** Synplexity **all** ≤3 kb + clonal tail | $255,529 | **100.0%** | **0.0%** | $114 |
| **S2** Synplexity all ≤3 kb + fragments tail | $219,823 | 93.7% | 4.2% | $105 |
| **S6** Synplexity all ≤3 kb + **stitched** tail @$44 | **$185,039** | 98.0% | 1.2% | **$84** |
| **S5** Synplexity all ≤3 kb + Elegen tail | $344,796 | 99.1% | 0.5% | $155 |

**⭐ S7 strictly dominates the old plan on both axes.** $221,446 vs the D baseline's
$224,478 — **slightly cheaper** — at **99.1% usable instead of 79.3%**, and severe errors
**13.3% → 0.5%**. It also beats F3 by **$55k while adding 10 points of usable library**.
**F3 is superseded as the recommendation.**

**Why S7 over S1.** S1 buys a literally perfect library (100.0% / 0.0%) for **+$34,083**.
The difference is entirely whether the 824 ORFs ≤1.2 kb go to AbstractBio at ~$0.01/bp
(~$12 each, 1:30,000) or to Synplexity at $50 each. Paying $50 to make a 600 bp gene perfect
when $6 makes it 99%+ clean is poor value — but **S1 is the right call if a single
error-bearing clone anywhere is unacceptable.** Both are defensible; S7 is the better trade.

**S4 is the trap.** Skipping Synplexity below 1.8 kb saves $22k but costs **8 points of
usable library** — the cheap short-end vendors are cheap precisely where per-bp pricing
favours them, but 1:3000 pools are the dirtiest thing in the plan. Don't optimise there.

**S6 is cheapest ($185,039, $84/usable) but not orderable today** — it depends on the same
unquoted Genyro/Sidewinder per-piece fee as scenario E. Pursue as the optimisation; S7 is
what can be ordered now.

⚠️ **Assumptions on Synplexity, all needing confirmation before this becomes the plan:**

1. **"Clonal" = sequence-verified, error rate 0.** Modelled as perfect, like Twist clonal.
   If it is only *mostly* clean, **every S-series quality number moves** — this is the single
   biggest dependency.
2. **$50 is flat at any length up to 3 kb**, with no per-bp component, no minimum order and
   no per-plate structure. A 3 kb gene at $50 is $0.0167/bp — 4× better than Twist clonal.
3. **The 3 kb cap is hard.** If it flexes to 3.5 kb, another 108 ORFs come inside
   (89.4% covered); at 5 kb it would be 95.4%.
4. **Turnaround unknown**, and 1,897 individually-clonal genes is a lot of separate
   deliverables. **Ask about delivery format** — individual tubes vs normalised plates vs
   pre-pooled — because at that count the logistics of pooling are real work, and a vendor
   that ships normalised plates or a pre-pooled library is worth a premium.

⚠️ **Even S1 does not fix the giants.** MUC16 (43,524 bp) and MUC3B (40,392 bp) are past
every vendor's ceiling including Twist clonal's — still truncation or a short isoform (§5.6).

*(Superseded: an earlier version of this section modelled a flat $50/gene with **no length
cap**, giving $55,897 / 90.9% and a "$102,267 at 99.9%" scenario. Those numbers assumed $50
would buy a 43 kb gene and are wrong. The 3 kb cap is what makes the real answer S7.)*

*(Model note: uses `MIS_SEV = 0.32` for comparability with §5.4's scenario table. §5.5's
measured 0.427 shifts every curve up ~1.5 points together, changing no ranking here.)*

### 5.5 AlphaMissense — severity weighting, and a licensing flag

**⚠️ LICENSE: AlphaMissense predictions are released under CC BY-NC-SA 4.0 —
non-commercial.** Confirmed from the Zenodo record (8208688). Using these outputs in
commercial R&D — or in anything published/shipped — likely needs legal review. Flagging
because Deliverome is a company and this is easy to miss. Permissively-licensed
alternatives if that's a blocker: **ESM-1v** (MIT), **ESM-2** derivatives, or classical
SIFT/PolyPhen-2.

**Headline parameter.** AlphaMissense classifies, across **all 71M possible** human
missense variants, **32% likely pathogenic / 57% likely benign / ~11% uncertain** (90%
precision cutoffs on ClinVar). This is the right reference class for *random synthesis
errors* — all possible substitutions, not a disease-enriched set. Adopted as
`MIS_SEV = 0.32`.

**Cross-checks from two orthogonal, experimental sources:**
- **VAMP-seq** (Matreyek et al., *Nat Genet* 2018) — 7,801 single-aa variants of PTEN +
  TPMT; **1,138 + 777 = 1,915 low-abundance ≈ 24.5%**.
- **Taipale — Lacoste et al., *Cell* 2024** — 3,448 missense across >1,000 genes:
  **mislocalization affects ~one-sixth of pathogenic missense**, and is *"primarily
  driven by effects on protein stability and membrane insertion rather than disruptions
  of trafficking signals."* **This is the most relevant endpoint we have**: for a
  surface-protein internalization screen, mislocalization *is* the phenotype — the
  protein never reaches the surface and the clone reads as a biological negative. Note
  the ~1/6 is conditional on *pathogenic* variants, so it is not directly substitutable.

#### ✅ RUN — and surface proteins are markedly LESS tolerant than the proteome average

Downloaded the full AlphaMissense per-variant table (1.2 GB) and filtered to our 2,242
accessions: **28,366,506 variants across 2,235 of 2,243 genes** (7 unmatched). Script:
`alphamissense_severity.py` · per-gene table: `alphamissense_by_gene.tsv` · figure:
`figures/s1e2_alphamissense_severity.svg`.

| | pathogenic | ambiguous | benign |
|---|---|---|---|
| **our surfaceome variants** | **12,486,217 (44.0%)** | 3,975,727 (14.0%) | 11,904,562 (42.0%) |
| AlphaMissense genome-wide | 32% | ~11% | 57% |

**Per-gene damaging fraction: mean 0.427, median 0.428, p10 0.246, p90 0.608.**

**So the flat 0.32 was an underestimate for this gene set — the real figure is ~0.43.**
Surface proteins are substantially more missense-intolerant than the proteome average,
which is mechanistically unsurprising: they are disulfide-rich, glycosylated, and must
fold *and* traffic through the secretory pathway.

**By topology class** — and the ordering is informative:

| DeepTMHMM class | n | mean frac. pathogenic |
|---|---|---|
| **TM** (multi-pass, no signal peptide) | 1,038 | **0.463** ← least tolerant |
| GLOB | 17 | 0.458 |
| SP (secreted / GPI) | 194 | 0.412 |
| **SP+TM** | 987 | **0.392** ← most tolerant |

**TM proteins are the least tolerant**, consistent with Taipale's finding that
mislocalization is driven primarily by **stability and membrane insertion** — multi-pass
TM proteins have the most insertion steps to get wrong.

⚠️ **But don't lean on that ordering.** AlphaMissense weights **evolutionary conservation
(MSA depth)** heavily, and TM helices are among the most conserved regions in the
proteome — so substitutions there score pathogenic partly *because* they're conserved.
The TM > SP+TM ordering is **consistent with** Taipale's mechanism but is not independent
evidence for it. See `REPORT_library_scoping.md` §5d.

**Effect on the model:** re-running scenario D with each gene's own damaging fraction
instead of the flat 0.32 raises the severe-error rate **13.3% → 14.8%**. Modest in
aggregate (the length term dominates), but it means every "severe" figure in §5.4/§5.8 is
**~1.5 points optimistic**, and the per-gene table lets you rank which genes most justify
paying for a clean (clonal) build.

Prior band was ~25–32% from the literature; **the measured value for our specific gene set
is 0.43**, so use the per-gene column where it matters.

⚠️ **Is 0.43 the right number? Not confidently — but it barely matters.** Two real
objections: the `am_class = pathogenic` cutoff was tuned for **90% precision on ClinVar**
(a clinical endpoint, and it leaves 14% "ambiguous" unassigned — count those and the
figure becomes **0.567**); and more importantly **clinical pathogenicity is the wrong
endpoint** — we need "does it still reach the surface", and chaining Taipale's ~1/6
mislocalization rate gives **0.427 × ~0.17 ≈ 0.07** instead. The saving grace is that a
**6× swing in this parameter moves the severe-error rate only 10.2% → 16.2%**, because the
**indel (30% of errors, all frameshifting) and nonsense terms dominate and don't depend on
it at all**. Sweep printed by `alphamissense_severity.py`; full argument in
`REPORT_library_scoping.md` §5d.

**Practical upshot:** keep **0.427** as the headline (measured on our actual genes, and
the conservative end), quote **10–16%** as the honest band on severe-error rate, and spend
the effort on the **Twist indel : substitution ratio** instead — that's the parameter that
actually drives the dominant term.

#### On BLOSUM — not the right tool here

BLOSUM/PAM measure **evolutionary exchangeability averaged over all proteins and
positions**. They are position-agnostic — blind to whether a residue sits in a TM helix,
a disulfide, or a disordered loop — and are strictly dominated for this purpose:

- **AlphaMissense publishes per-variant, per-position scores for the whole proteome**, so
  for a defined gene list you get the *actual* expected damaging fraction per gene rather
  than a flat rate (license permitting).
- **Taipale's dataset** gives the trafficking-specific endpoint.
- **Position weighting matters most for us**: TM domains and disulfide-bonded ectodomains
  are the least tolerant — and we already pulled `per_residue_topology`,
  `ecd_length_residues` and `tm_helix_count` for all 2,243 in the same D1 query, so a
  TM-vs-loop weighted estimate is a join away.

Also a framing correction: synthesis errors are **random single-nucleotide changes, not
designed codon recoding**. The right computation is: enumerate single-nt changes over the
**real codon usage of the actual sequences** → derive the true aa→aa substitution
spectrum → score each with an empirical severity source. BLOSUM only as a sanity
baseline.

### 5.6 Long-ORF QC — likely non-surface entries worth reviewing

Using DeepTMHMM topology pulled with the lengths. Of the **870 ORFs >1.8 kb**:

| label | all 2,243 | long only |
|---|---|---|
| SP+TM | 993 | 520 |
| TM | 1,038 | 291 |
| SP (no TM) | 195 | 54 |
| **GLOB** (no TM, no SP) | 17 | **5** |

**⚠️ 59 of 870 long ORFs have ZERO predicted TM helices** (54 SP-only + 5 GLOB).

**The 5 long GLOB entries — no TM, no signal peptide, ECD = 0 residues.** DeepTMHMM
predicts no membrane insertion *and* no extracellular residues, so these look like
mis-inclusions for a membrane-receptor library:

| gene | CDS bp | aa | TM | ECD |
|---|---|---|---|---|
| **CATSPERT** | 5,463 | 1,820 | 0 | 0 |
| **DIP2A** | 4,716 | 1,571 | 0 | 0 |
| **PITPNM3** | 2,925 | 974 | 0 | 0 |
| **ENOX1** | 1,932 | 643 | 0 | 0 |
| **ENOX2** | 1,833 | 610 | 0 | 0 |

*(DIP2A, PITPNM3 and ENOX2 are among PR #130's newly-admitted 791 — worth a second look
at whether the relaxed gate let these in.)*

**The 54 long SP-only entries are secreted / ECM, not integral membrane** — e.g.
**HSPG2** (perlecan, 13,176 bp), **PKHD1L1** (12,732), **CUBN** (cubilin, 10,872),
**LAMA3** (10,002), **COL6A3** (9,534). Legitimately surface-*associated*, but a
different category from an anchored receptor, and they're long/expensive.

**Extreme tail — beyond any practical synthesis limit.** **MUC16 = 43,524 bp (14,507 aa)**
and **MUC3B = 40,434 bp** cannot be ordered as fragments at all (Twist fragments top out
~5 kb). The giant mucins (MUC16, MUC3B, MUC12, MUC4, MUC17, MUC3A) are also mostly
**repetitive PTS/tandem-repeat domains** — terrible synthesis substrate *and* a
recombination hazard. For a binding/internalization screen the membrane-proximal region +
TM is almost certainly what's wanted, not 14,000 residues of mucin domain.

#### What trimming would save (scenario D baseline $224,478)

| Adjustment | ORFs dropped | New cost | Saves |
|---|---|---|---|
| drop 17 GLOB (no TM, no ECD) | 17 | $223,217 | $1,261 (0.6%) |
| drop 195 SP-only (secreted/ECM) | 195 | $208,163 | $16,315 (7.3%) |
| drop both → TM / SP+TM only | 212 | $206,902 | $17,576 (7.8%) |
| exclude CDS > 10 kb | 21 | $201,158 | $23,319 (10.4%) |
| exclude CDS > 7.5 kb | 41 | $189,139 | $35,339 (15.7%) |
| **exclude CDS > 5 kb** (practical fragment ceiling) | **103** | **$162,404** | **$62,074 (27.7%)** |
| **combined: TM/SP+TM only AND ≤5 kb** | **304** | **$151,561** | **$72,917 (32.5%)** |

**Excluding the 103 ORFs >5 kb saves 27.7% of the budget and removes constructs that
mostly cannot be synthesized as single fragments anyway** — the single highest-leverage
adjustment available.

### 5.7 Scope extensions — isoforms and orthologs

All from the same read-only public D1 mirror (`topology_public`, `compara_ortholog_ecd`).
**Two data hygiene points that materially changed the numbers:**

1. `topology_public` carries **multiple `topology_version` rows per sequence**
   (`topo_2026_05_16`, `topo_2026_05_25`, `topo_test_optA_full`) — a naive query
   **double-counts**. All numbers below dedupe to the latest version per
   `uniprot_acc_full`.
2. Orthologs must be joined via **`compara_ortholog_ecd.ortholog_uniprot_acc`**, not by
   upper-casing the ortholog gene symbol. Symbol matching found only 302 mouse
   orthologs; the proper key finds **1,777**.

Secreted-only (`deeptmhmm_label = SP`) sequences are **excluded throughout**, per scope.

#### Isoforms of our 2,243 genes

**2,513** non-secreted-only isoforms exist across **1,134** of our genes. Filtering to
"substantially different from the canonical" (both a % and an absolute-residue floor, so
short proteins aren't admitted on a small absolute change):

| Difference threshold | isoforms | genes | bp | **cost** |
|---|---|---|---|---|
| ≥10% **and** ≥50 aa | 1,076 | 662 | 1,518,168 | **$59,360** |
| **≥20% and ≥100 aa** | **632** | **410** | **783,129** | **$26,327** |
| ≥30% **and** ≥150 aa | 429 | 296 | 463,341 | **$11,984** |

The ≥20%/≥100 aa tier looks like the sensible default: **632 isoforms across 410 genes
for $26k**, ~12% on top of the $224k canonical build.

Worth noting the biggest-difference isoforms are mostly **domain-loss variants of the
giant genes** — e.g. MUC4 `Q99102-12` (5,412 → 1,125 aa, −79%), CDH23 `Q9H251-11`
(3,354 → 216 aa, −94%, and DeepTMHMM calls it **GLOB** so it has lost membrane
insertion entirely). Several of these are *cheaper and more tractable* than their
canonical parents — for the >5 kb genes flagged in §5.6, **a short isoform may be the
practical construct**.

#### Close mouse and cyno orthologs

"Close" = ECD percent identity from `compara_ortholog_ecd`:

| Threshold | mouse | | cyno | |
|---|---|---|---|---|
| | orthologs | cost | orthologs | cost |
| all (non-SP) | 1,777 | $198,570 | 1,722 | $179,764 |
| ≥70% ECD id | 1,451 | $180,690 | 1,631 | $173,583 |
| ≥80% ECD id | 1,119 | $147,584 | 1,560 | $169,824 |
| **≥90% ECD id** | **644** | **$92,292** | **1,350** | **$156,477** |
| ≥95% ECD id | 361 | $51,819 | 1,018 | $125,781 |

**Cyno is far more conserved than mouse**, as expected: at ≥90% ECD identity cyno retains
**1,350** orthologs vs mouse's **644**. So a stringent identity cut costs cyno relatively
little coverage but halves the mouse set.

#### Combined budget sketch (scenario D routing throughout)

| Component | cost |
|---|---|
| canonical human, 2,243 ORFs | **$224,478** |
| \+ isoforms (≥20% / ≥100 aa, 632) | +$26,327 |
| \+ mouse orthologs ≥90% ECD id (644) | +$92,292 |
| \+ cyno orthologs ≥90% ECD id (1,350) | +$156,477 |
| **total** | **≈$499,574** |

Orthologs are the expensive add — **$249k for the two species at ≥90%**, more than the
entire human canonical build. If cross-species is needed, consider a **prioritised
subset** (e.g. only the genes that actually score in the human screen) rather than
synthesising both full ortholog sets up front. Per the user's instruction, **AlphaMissense
was not run on isoforms or orthologs.**

### 5.8 Redundant builds — does synthesising long ORFs twice pay off?

**How the strategy works.** Order *N* independent builds of the same construct, sequence
them, and keep whichever is clean. Probability at least one of *N* is error-free is
`1 − (1 − q)^N` where `q = (1 − p_err)^L`. **This only works if you sequence and select** —
duplication without selection just doubles the number of unknown-quality molecules. Since
we're already doing long-read variant↔barcode mapping (§5.1), the selection step is free.

**Longest CDS still meeting the target, by number of builds:**

| criterion / target | 1 build | 2 | 3 | 4 |
|---|---|---|---|---|
| **≥95% no severe error** (fragments 1:7500) | **786** | 3,880 | 7,044 | 9,815 |
| ≥85% no severe error | **2,491** | 7,509 | 11,617 | 14,926 |
| **≥95% no non-synonymous error at all** | **464** | 2,289 | 4,156 | 5,792 |
| ≥85% no non-synonymous error | **1,470** | 4,431 | 6,855 | 8,807 |

*(Same table at the pool rate 1:3000 is ~2.5× shorter; at AbstractBio 1:30,000 ~4× longer.)*

**So the answer to "at what length": duplication starts paying at ~786 bp** for the
95%/no-severe target and **~2,491 bp** for 85%/no-severe. Below those a single build
already clears the bar; above them you need 2, then 3, then more.

**What it costs** (baseline scenario D, 1 build each = $224,478; duplication applied only
where needed):

| target | cap | cost | vs base | ORFs still short |
|---|---|---|---|---|
| ≥95% no severe | 2 | $441,839 | 1.97× | 383 |
| ≥95% no severe | 3 | $530,606 | 2.36× | 47 |
| **≥85% no severe** | **2** | **$394,235** | **1.76×** | **41** |
| ≥85% no severe | 3 | $429,574 | 1.91× | 16 |
| ≥85% no non-syn | 3 | $513,010 | 2.29× | 54 |
| ≥95% no non-syn | 3 | $627,629 | 2.80× | **247** |

**Two conclusions:**

1. **"95% of sequences with no non-synonymous error at all" is not a reachable target by
   duplication.** A single fragment build only clears it to 464 bp, and even 3 builds
   leaves **247 ORFs short at 2.80× cost ($628k)**. Drop this target or change approach.
2. **"85% with no severe error" is very reachable — $394k at 2 builds (1.76×), only 41
   ORFs short.** That's the cheapest credible quality bar.

**But stitching dominates duplication.** Duplication is brute force: you pay 2–3× to buy
more lottery tickets at the same per-base error rate. **Stitching from short verified
pieces attacks the error rate itself.** Compare single-build reach at ≥95% no-severe:
**786 bp** (Twist fragments 1:7500) vs **3,145 bp** (AbstractBio-grade 1:30,000) — a **4×
length improvement with no duplication at all**. Since 94.5% of our budget is the 870 long
ORFs, lowering their effective error rate is worth far more than duplicating them.

→ **Price Genyro / Sidewinder stitching before committing to a duplication strategy**
(§5.4). Duplication is the fallback if stitching isn't viable.

**⚠️ SUPERSEDED — duplication is off the table.** Once Twist clonal ($0.09/bp,
sequence-perfect) is in the picture, it **strictly dominates** duplication on the long tail
that duplication was meant to rescue: **$272,675 @ 100% usable ($313/usable ORF)** vs
**$424,161 @ 89.5% ($545/usable ORF)** for 2× fragments. Cheaper *and* perfect. The
thresholds above are retained because they answer the question as asked and still frame the
quality bars, but **the plan is F3 (clonal >~2.3 kb), not duplication** — see §5.4.
Duplication remains relevant only in the extreme tail, where nothing helps much anyway: a
second build takes MUC16 from 5.8% to 11.4% at 2× the cost.

### 5.9 Additional candidates — the "likely" tier with a membrane anchor

Ran PR #130's `passes_likely` over the same cohort: **625 genes are "likely" but not
canonical.** Filtering to those that could actually be screened:

- **≥1 TM helix OR GPI-anchored: 444 genes**, 714,942 bp, **$29,403**
  (433 with TM≥1, 12 GPI-anchored; subcategories: 146 multi_pass, 110 single_pass_T1,
  96 GPCR, 40 single_pass_T2, 16 tetraspanin).
- **After the same screen-ability filters** (no co-receptor requirement, ECD ≠ none,
  ≤5 kb, no OR/TAS2R): **330 genes, 499,710 bp, $19,931.**

**This is the best marginal value in the whole analysis: +330 screenable targets for
~$20k — +15% more genes for +9% cost**, because they're short (mean ~1,514 bp vs 2,032
for canonical). Cheapest additions are tiny single-pass and tetraspanin-family proteins:
SMAGP (294 bp), TCTA (312), NRG4 (348), VAMP5 (351), IFITM2/3 (399/402), TMEM100 (405),
RNASEK (414), SLC31A2 (432), plus GPI-anchored EPGN (465) and RNASE2 (486).

Worth considering, given canonical already excludes them only on evidence strength rather
than on being non-surface — and at this price the downside of a few false positives is
small relative to missing real targets.

### 5.10 How the field handles error rate in long ORFs

Standard mitigations, roughly in order of how well they fit our design:

1. **Don't chase a perfect library — sequence what you built.** For pooled libraries the
   accepted answer is **long-read variant↔barcode mapping** (PacBio/Nanopore full-length
   linkage), then analyse only clones whose sequence is confirmed. This is already §5.1's
   recommendation, and it simultaneously resolves recombination mispairing. **Best fit
   for us.**
2. **Exploit redundancy.** At ~13 clones/gene, per-gene representation survives a 20%
   error rate (~10.6 clean clones/gene). Increase coverage rather than perfect the
   synthesis.
3. **Clonal, sequence-verified genes** instead of pools/fragments — you buy a guaranteed
   perfect sequence. ✅ **Corrected: this is far cheaper than the 2–4× premium I originally
   assumed — Twist quotes $0.09/bp vs $0.07/bp for fragments, only 1.29×.** At that
   multiple it is **not prohibitive at 2,243 genes**; it is the recommended route above
   ~2.3 kb (§5.4). This is the single largest revision in the synthesis analysis.
4. **Build long ORFs from shorter verified fragments.** Sub-kb pieces have far lower
   per-molecule error; verify each piece cheaply, then Golden Gate/Gibson assemble.
   Aligns naturally with the Fragmid architecture we're already using.
5. **Enzymatic error correction** (mismatch-cleavage, e.g. ErrASE-class) — typically
   ~5–10× error reduction, applied post-synthesis.
6. **Design-level shortening** — ectodomain-only or truncated constructs, or splitting
   very large receptors. This is what §5.6's >5 kb exclusion amounts to.
7. **Accept and model it** — report the expected clean fraction and power the screen
   accordingly, rather than assuming a perfect library.

### 5.11 Other open gaps

- **What actually varies across the library is never stated** — and it is load-bearing
  for the §5.1 span. Define it.
- **`atgtac` → non-native M-Y** on TfR1's N-terminus (native `MMDQARSAF…`) — "issue 7".
- **Kozak block source** unstated; **GS-linker sequences** have no source.
- **Stray `n`** in the barcode region; confirm element 11 (21 nt, BglII) is a scar.
- **pRDA_722's dropout region (2617–3042, ~425 bp)** — confirm our cassette replaces
  exactly that span.
- **Packaging plan** (helper/envelope plasmids, titer, MOI) and **construct-verification
  plan** (full-plasmid long-read vs junction Sanger) are unrecorded.
- **Specific host line** unrecorded — the IRES A/B should run in the intended line.

### 5.12 OPEN DECISION — PuroR as a separate transcript, or IRES-linked?

**Current design:** IRES-linked (`IRES2 → mScarlet3 → P2A → PuroR`, §4.3).
**Alternative:** its own promoter — MORF's `SV40 → PuroR`, RESOLUTE's
`PGK → rtTA-VP16–2A–PuroR` (§4.3).

**Why this is a risk and not just a preference:** the whole marker cassette currently sits
downstream of an IRES2 that is **functionally unvalidated** (§5.2) and **not independently
adopted by anyone** — Fragmid is the only Addgene vector carrying this domesticated
sequence. If it under-performs we lose **selection and readout together**, and we find out
late, after the library is built.

**What each option actually buys — stated precisely, because the earlier version of §4.3
overstated the linkage:**

| | IRES-linked (current) | separate promoter |
|---|---|---|
| filters promoter silencing / 5′ truncation / large deletions | ✅ | ❌ |
| filters ORF point mutations + indels — **the §5.4 error budget** | ❌ (IRES initiates internally; upstream frameshift doesn't block PuroR) | ❌ |
| selection survives a weak or mildly toxic ORF | ❌ | ✅ |
| depends on the unvalidated IRES2 | ⚠️ **yes** | ✅ no |
| costs a second promoter (size, and possible promoter interference) | ✅ no | ⚠️ yes |

**So the strong "resistance proves the ORF is intact" argument does not hold for an IRES**
— only for 2A, which we can't use because the ORF needs a free C-terminus (§4.3). That
weakens the case for staying IRES-linked considerably: what we're really buying is a
transcript-integrity check, and paying for it with a dependency on an unvalidated element.

**Decide by:** the §5.2 ratiometric A/B in p006/p007. **Pre-commit the branch now** —
≥2-fold reduction → move PuroR to its own promoter rather than hunting another IRES variant.
⚠️ Do not defer this until after the result, or the decision gets made under schedule
pressure.

### 5.13 OPEN DECISION — separate UMI, or let the barcode serve both roles?

**RESOLUTE carries a barcode *and* a UMI; Tedman et al. carry a UMI** (§4.3, §4.11). We
currently have a 20-nt barcode and no UMI. The two do different jobs:

- **Barcode** → *which construct is this?*
- **UMI** → *which individual clone / molecule is this?* (distinguishes independent
  integration events, and PCR duplicates from real reads)

**✅ UPDATED — a UMI is now affordable without lengthening the construct.** §4.4's
edit-distance analysis shows **16 nt clears 30,000 barcodes at min Levenshtein 3**, so the
4 freed bases give **256 molecular UMIs at no change in total length** — and the barcode
gains single-edit error correction it never had. That reframes this decision: it is no longer
"pay bases for a UMI" but "reallocate 4 bases from unused barcode space." ⚠️ The cost is
giving up min distance 4, which would need 18 nt and leave only 16 UMIs.

**Original provisional answer (retained — still valid if we keep 20 nt of barcode): we don't
need a separate UMI, because our barcode space already dwarfs the library.** 20 nt gives ≈ **6.5 × 10¹¹** effective sequences after
constraints (§4.4) against a **30,000**-clone library. If barcodes are assigned **per clone**
rather than per ORF — which is what §5.1's empirical long-read variant↔barcode map assumes
anyway, since the designed pairing can't be trusted — then **the barcode already identifies
the clone and is doing the UMI's job.** A separate UMI is only necessary if we assign one
barcode per ORF.

**Why not add one anyway:** every extra base between the variable region and the barcode
**lengthens the §5.1 template-switching span** and so *increases* mispairing — the one
quantity we are trying to minimise. A UMI is not free.

**What would change the answer:** if we end up needing to count independent integrations
per barcode (e.g. to weight a hit by how many clones support it), or if PCR duplicate
collapse turns out to matter for the long-read map. Both are analysis-plan questions, not
construct questions — so **resolve the analysis plan first, then decide.**

**⚠️ Placement, if we do add one:** put it **inside** the shuttle sites and **upstream of
the barcode**, per Tedman (§4.11), so it travels with the ORF and the long-read map is
built once.

### 5.14 OPEN DECISION — surrounding Golden Gate sites so others can reuse the ORFs

**The decision:** flank the insert with a **third Type IIS enzyme's sites** so the whole
barcoded ORF collection can be lifted into a different backbone in **one Golden Gate
reaction** — by us, and by anyone we share it with. Full analysis, enzyme options and the
stop-codon compatibility check are in **§4.11**; this entry exists so it sits with the other
open decisions rather than only in the architecture section.

**Why it belongs here and not just in "nice to have":**

- **Reuse by others is the point.** A 2,243-ORF barcoded surfaceome is a reagent other labs
  will want in *their* backbone — different promoter, different marker, different assay.
  Without flanking sites, every recipient has to re-clone from scratch; with them it's one
  reaction. Tedman et al. explicitly designed for this, calling their library
  "backward-compatible with several current activity assays and forward-compatible with
  future assays" (§4.11).
- **It protects our own biggest sunk cost.** §5.1 forces an empirical long-read
  variant↔barcode map because ~9% of clones are mispaired. If the barcode travels with the
  ORF through the shuttle, that map is built **once**; without shuttle sites, *every* future
  backbone needs it re-derived. This is the self-interested argument, and it's stronger than
  the altruistic one.
- **It is only cheap if decided now.** Adding sites later means re-synthesising or re-cloning
  the whole library. Adding them at design time costs ~6–7 bp per flank.

**What's actually undecided:** not *whether* — **which enzyme**, and that reduces to a single
question: **is the Esp3I/BsmBI reservation for downstream guide cloning worth more than a
clean shuttle?**

- **Esp3I/BsmBI** is already at **0 sites** in p003, p004 *and* pRDA_722 (§4.10) and is what
  Tedman used — cleanest, best-precedented, zero extra domestication. Blocked only by our
  own reservation (§1).
- **AarI/PaqCI** is the best unreserved option — one IRES2 site to remove.
- **BsaI is the worst option**, despite being the obvious-sounding one: §4.10 shows **2–4
  BsaI sites in every construct**, so a BsaI shuttle would cut the backbone.

**Decide by:** before the TFRC CDS is re-derived (§4.2), since the flanks and the ORF's
terminal sequence are designed together. **Blocking for synthesis**, not for the p006/p007
pilots.

**Then verify:** re-run `enzyme_scan.py` and confirm the chosen enzyme cuts at **exactly the
two new flanks and nowhere else.**

---

## 6. Next steps — prioritised

**Blocking / before ordering synthesis**

- [ ] **Change both stops `TGA` → `TAA[A,G]`.** Highest priority at the ORF stop:
      readthrough extends past the ALFA tag (can occlude NbALFA) *and* raises WPRE X-ORF
      translation (§4.5, §5.3).
- [ ] **Define what varies across the library** (§5.4) — needed to fix the §5.1 span.
- [ ] **Decide snorkel in/out**, then **finalise barcode placement** to minimise the
      variable-region→barcode span (§4.4, §5.1).
- [ ] **Diversify GS-linker codons** — removes the 9× 12-mer repeat cluster (failure mode
      2). Cheap, no protein change (§5.1).
- [ ] **Clear the stray `n`**; confirm element 11 is a cloning scar (§3).
- [ ] **Re-derive the TFRC CDS with native codon usage** — minimal site-only synonymous
      edits for BbsI/Esp3I, not whole-ORF recoding (§4.2).

**Before committing the design**

- [ ] **Validate `IRES2_v1`** by ratiometric flow in p006/p007; pre-specify the
      ≤20% / ≥2-fold thresholds (§5.2).
- [ ] **Pre-commit to the IRES-failure branch:** if the A/B fails, move PuroR to its own
      promoter (MORF `SV40 → PuroR` / RESOLUTE `PGK → rtTA-VP16–2A–PuroR`, §4.3) rather
      than hunting another IRES variant. Accepts loss of ORF↔marker linkage — decide now,
      not after the result (§5.2).
- [ ] ⭐ **Decide on the SPLIT barcode: 13-nt ORF-ID prefix (min Levenshtein 3) + 7-nt
      degenerate clone suffix (§4.4).** Same 20-nt budget as now, but it encodes ORF identity
      with error correction **and cuts wrong-ORF mispairing from ~8.96% to ~1.05%** by moving
      ORF identity next to the stop — the largest available improvement to §5.1, which had
      been recorded as not design-fixable. Depends on per-ORF synthesis carrying the prefix
      (Synplexity does; a pooled degenerate barcode does not). **Verify the ~20 bp span
      assumption for close paralogs / isoform pairs identical near the 3′ end.**
- [x] ✅ **Production barcode set generated — `barcodes_20nt_30k.tsv`**, 30,000 × 20 nt under
      the full stack at min Levenshtein 3, independently verified (§4.4). Regenerate with
      `generate_barcodes.py`; deterministic from SEED=20260812.
- [ ] **Port the full published filter stack into the production barcode generator** (§4.4):
      GC 40–60%, no homopolymer triple, no `ATG`, `AATAAA` + point mutants, **BbsI + Esp3I +
      BsaI** both strands, non-overlapping 3-bp self-complementarity, `conserved+` 7mer-m8
      seeds, min Levenshtein 3. **16 nt clears 30,000 under the whole stack**, so it costs
      nothing at that length.
- [x] ✅ **Slutskin's barcode Methods confirmed verbatim** against [PMC5802814](https://pmc.ncbi.nlm.nih.gov/articles/PMC5802814/)
      (§4.4). Confirms 12 nt, "differ ... by 3 nt or more", highly-expressed-miRNA seeds, and
      `AATAAA` + point mutants — and confirms the seed is defined as bases **2–7**, a 6-mer.
- [x] ✅ **6-mer seed adopted** (Slutskin's definition, bases 2–7). Our 7mer-m8 was strictly
      less conservative — every 7mer target contains its 6-mer core — and at 20 nt the extra
      cost is free against ~10⁶× headroom (§4.4).
- [x] ✅ **Min Hamming 4 added alongside min Levenshtein 3.** Levenshtein alone permitted 71
      Hamming-3 pairs; Hamming alone gives no indel guarantee. Both now enforced (§4.4).
- [ ] **Fix the position-0-only ATG check** — ~25% of passing barcodes carry an internal ATG
      (§4.4).
- [ ] **Impose a minimum Levenshtein distance in the production generator** — there is
      currently none, so nothing guarantees error-correcting separation (§4.4).
- [ ] **Settle the host line (§5.11), then re-derive the seed blacklist from miRNAs actually
      expressed in it** — conservation tier is only a proxy (§4.4).
- [ ] **Confirm TargetScan redistribution terms** for `mirna_seed_sites.tsv` (§4.4).
- [ ] **Decide whether to add a UMI** — provisionally **no**: our 20-nt barcode space
      (6.5 × 10¹¹) already dwarfs a 30,000-clone library, so a per-clone barcode does the
      UMI's job, and extra bases lengthen the §5.1 mispairing span. **Settle the analysis
      plan first** (do we need to count independent integrations per barcode?), then confirm
      (§5.13).
- [ ] **Verify the MORF / REPROcode / RESOLUTE specifics in §4.3** against primary sources
      — currently user-supplied. In particular confirm REPROcode really ships no mammalian
      drug-resistance marker, and get the REPROcode citation.
- [ ] **Understand why RESOLUTE codon-optimized** its SLC ORFs before treating §4.2's
      "don't codon-optimize native ORFs" as settled — SLCs are the multi-pass TM case where
      co-translational folding should matter most, so this is a real counterexample (§4.3).
- [ ] **Decide the WPRE X-ORF question** — `mut6` is ~free (§5.3).
- [ ] ⚠️ **Add surrounding Golden Gate sites so others can reuse the ORFs, and pick the
      enzyme (§5.14, detail in §4.11).** Flanking Type IIS sites let the whole barcoded
      collection move to any backbone in one reaction — for other labs, and for us. Only
      cheap if decided **now**: retrofitting means re-cloning the library, adding it at
      design time costs ~6–7 bp per flank. **Blocking for synthesis** (settle it before the
      TFRC CDS is re-derived), not for the p006/p007 pilots. The choice reduces to one
      question: **is the Esp3I reservation worth more than a clean shuttle?** Esp3I/BsmBI is already at **0 sites** in p003/p004/pRDA_722 and is what
      Tedman et al. used, so it is the cleanest and best-precedented choice — *if* downstream
      guide cloning can move to another enzyme. Otherwise **AarI/PaqCI** (one IRES2 edit).
      **Not BsaI** — §4.10 shows 2–4 BsaI sites in every construct, so it is the only
      candidate needing new domestication work.
- [ ] **Then add the shuttle sites and re-run `enzyme_scan.py`** to confirm the chosen enzyme
      cuts *only* at the two new flanks (§4.11).
- [ ] **Confirm pRDA_722's dropout span (2617–3042)** matches our cassette (§5.4).
- [x] ✅ **Esp3I coverage — largely closed by direct rescan (§4.9, §4.10).** The worry was
      that the manual scan was described as "**BbsI and BsaI**" while our documented pair is
      **BbsI + Esp3I/BsmBI** (`CGTCTC`/`GAGACG`) — **BsaI (`GGTCTC`) being a different
      enzyme.** Rescanning every construct GenBank directly settles the outcome regardless
      of what was scanned: **p003 and p004 are BbsI = 0, Esp3I = 0.** ⚠️ Two residual
      caveats: (a) **BsaI was never cleared anywhere** (1223, 2854 in all constructs;
      pRDA_722 has 3 of its own) — harmless for our enzyme pair but don't read
      "domesticated" as "Type IIS-free"; (b) the rescan covers the **assembled plasmids**,
      so the **not-yet-derived TFRC CDS still needs its own BbsI/Esp3I scan** when it is
      re-derived with native codon usage.
- [ ] Confirm the **Kozak block's source**, and whether the non-native **M-Y** is intended.
- [ ] Confirm the exact ID behind *pAF0014* — likely **AA014 / #216057**; no plasmid
      literally named `pAF0014` exists in Addgene's public catalog.
- [ ] Definitively confirm the WPRE X-ORF status by aligning vs the WHV X reference and
      `mut6` (§5.3).
- [ ] Independently confirm **V291** for p007 against the EndoNB supplementary table
      (§4.10).
- [ ] **Check the construct GenBanks into version control** (or S3). §4.9/§4.10's enzyme
      counts come from local-only files in `~/Downloads`, so `enzyme_scan.py` is currently
      reproducible on one machine only. Plasmid maps are small and belong under git.
- [ ] **Confirm how Esp3I 5445 (`pac`/PuroR) was cleared in p003/p004** — presumed to be
      Fragmid's PuroR swapped in rather than a synonymous point edit, but not verified. If
      it *was* point-edited, confirm the change is synonymous (§4.10).

**Synthesis / library scope (§5.4–5.7)**

- [ ] ⚠️ **HIGH — fix the broken S3 pre-push gate.**
      `scripts/check-experiment-outputs-uploaded.sh` reads
      `os.environ.get("S3_ACCESS_KEY_ID")` **before** calling `load_env()`. Credentials
      live in `.env` and are not exported, so the guard short-circuits to `sys.exit(0)`
      and the hook reports **"Passed" without ever checking S3**. The gate has therefore
      **never actually run** for anyone using the documented `.env` workflow — which is
      exactly how s1e2's outputs sat un-uploaded while every push reported green. Fix:
      call `load_env()` before the guard. *Deliberately not changed in this PR* — turning
      a previously-inert gate on will start blocking pushes repo-wide, so it should be an
      announced change, not a side effect.
- [ ] ⚠️ **HIGHEST — confirm Synplexity's terms, then adopt S7** (§5.4). At $50/gene flat
      to 3 kb, **S7 strictly dominates the old plan: $221,446 (slightly under the $224,478
      baseline) at 99.1% usable vs 79.3%, severe errors 13.3% → 0.5%.** Four things to
      confirm, in order of how much they matter:
      **(a) is "clonal" sequence-verified (error 0)?** — the single biggest dependency; every
      S-series quality number moves if it is only *mostly* clean.
      **(b) is $50 truly flat to 3 kb**, no per-bp component, no minimum order, no per-plate
      structure?
      **(c) is the 3 kb cap hard?** At 3.5 kb another 108 ORFs come inside (89.4% covered);
      at 5 kb, 95.4%.
      **(d) turnaround and delivery format** — 1,897 individually-clonal genes is a lot of
      deliverables; normalised plates or a pre-pooled library would be worth a premium over
      individual tubes.
- [ ] ⚠️ **HIGH — get Twist's clonal-gene length ceiling.** It gates the whole routing
      decision: if clonal reaches 5 kb, **F3 is the answer and Elegen is unnecessary**; if
      it caps at ~3.5 kb, covering the remainder with Elegen costs **+$72k** for slightly
      *lower* usable output (§5.4).
- [ ] **Ask Twist for the indel : substitution ratio.** Assumed 30%; **the biggest single
      lever** on the severe-error estimate — it drives the dominant term, and matters far
      more than refining the missense-severity cutoff, which is worth only ±3 points
      (§5.4, §5.5).
- [ ] **Decide vendor routing** — recommend ⭐ **S7: AbstractBio ≤1.2 kb + Synplexity
      1.2–3 kb + Twist clonal >3 kb ($221,446, 99.1% usable, 0.5% severe).** **Supersedes
      F3**, which cost $55k more for 10 points less. Variants: **S1** (Synplexity for
      everything ≤3 kb) buys a literally perfect library for +$34k — take it if a single
      error-bearing clone is unacceptable; **S6** (stitched tail) is cheapest at $185,039 but
      needs the unquoted per-piece fee. ⚠️ All contingent on (a)–(d) above (§5.4).
- [ ] **Verify Elegen Infinia's specs before modelling on them** — $0.14–0.20/bp and
      1:70,000 are both **unverified**. Only relevant above the clonal ceiling; ask about
      turnaround too (§5.4).
- [ ] **Get Genyro / Sidewinder per-piece quotes.** Break-even vs scenario D is
      **$44/piece**; below that, AbstractBio-grade pieces + stitching beats D on cost *and*
      quality (94.7% vs 79.3% usable, 4× fewer severe errors) — the best long-run option,
      but it needs vendor qualification, so pursue it in parallel with clonal as the
      interim default (`REPORT_library_scoping.md` §2).
- [ ] **Manually review the 38 concerning rejects** — high/moderate-confidence,
      surface-accessible genes rejected only on `evidence_grade = weak`: **ANTXR2** (an
      established toxin-delivery receptor), **EPHB1**, **TNFRSF11A** (RANK), **EFNB3**,
      **EFNA2**, LCT, SI, CPM, ART4, XG, GP5, EDAR … These look like corpus-coverage
      failures, not biology failures (`REPORT_library_scoping.md` §1).
- [ ] **Do NOT bundle the cassette into the ORF orders** — it pushes 231 ORFs over the
      1.8 kb fragment threshold and costs **+$46k**. Fragmid keeps the cassette in
      separate fragments / the destination vector, so this is avoidable (§5.4).
- [ ] **Review the 5 GLOB long ORFs** (CATSPERT, DIP2A, PITPNM3, ENOX1, ENOX2 — no TM, no
      ECD) and the **54 long SP-only** secreted/ECM entries for membrane-library fit
      (§5.6).
- [ ] **Decide the >5 kb policy.** Excluding the 103 ORFs >5 kb saves **27.7%** and drops
      constructs that mostly cannot be made as single fragments anyway. **MUC16 (43.5 kb)
      and MUC3B (40.4 kb) are unsynthesisable as-is** — truncate to membrane-proximal +
      TM, or use a short isoform (§5.6, §5.7).
- [ ] **Decide isoform scope** — ≥20%/≥100 aa tier = 632 isoforms / 410 genes / **$26k**
      (§5.7).
- [ ] **Decide ortholog scope** — ≥90% ECD identity = mouse 644 (**$92k**) + cyno 1,350
      (**$156k**). At **$249k the pair costs more than the entire human build**; consider
      gating on genes that actually score in the human screen first (§5.7).
- [ ] ⚠️ **Resolve the AlphaMissense licence question before any commercial or published
      use** — predictions are **CC BY-NC-SA 4.0 (non-commercial)**. ESM-1v (MIT) is the
      permissive fallback (§5.5).

**Recombination QC (accept-and-measure)**

- [ ] **Assume ~9% mispairing** (~2,700 of 30,000) and build the **empirical
      variant↔barcode table by long-read sequencing** — don't trust the designed pairing.
      Use **low MOI** (§5.1).
- [ ] **Re-run the span + repeat scan on the final corrected insert.**
- [ ] **Place the §4.11 shuttle site upstream of the barcode, not between barcode and ORF**
      — so the barcode travels with the ORF and the expensive long-read variant↔barcode map
      is built **once** and reused across every future backbone. This is why Tedman et al.
      put theirs upstream of the UMI (§4.11).

**Housekeeping**

- [ ] **Pause/cancel the PE anti-CD271 ME20.4 order** (BioLegend #6428) — likely moot now
      that tNGFR is out (see Decision log).
- [ ] **Re-plan the fluorescence panel** — mScarlet3 occupies the red channel (§4.3).
- [ ] **Integrate the miRNA control system** and fix its placement (§4.8).
- [ ] `deliverome_analysis.codon_opt` (PR #10) needs a **"minimal site-avoidance recode"**
      mode; the whole-ORF recoder no longer applies to native ORFs (§4.2).
- [ ] Apply the native-codon, minimal-edit policy to future ORF payloads (§4.2).

---

## Inputs

- **Parent (domesticated build):** Addgene [#216077](https://www.addgene.org/216077/)
  (`pRDA_722`, Fragmid destination, 7,965 bp) — also the source of the **WPRE**.
- **Parent (p006/p007):** Addgene [#236083](https://www.addgene.org/236083/)
  (`pAG Lenti CMV Empty Puro`, 7,119 bp, dep. Melina Fan) — deliberately non-domesticated.
- **TFRC** canonical CDS via Ensembl; cross-checked against Addgene
  [#61505](https://www.addgene.org/61505/) (TfR–pHuji; Shen et al. 2014 *JCB*).
- Fragmid fragments: `pF_AA197` (#216049, ALFA), `pF_AA333`/`AA333` (#216067, IRES2 +
  PuroR), `AA014` (#216057, P2A + PuroR).
- `pOT_2` (#181971, P2A-Puro), `SpyTag003-MBP` (#133450), **FPbase** (mScarlet3).
- Local GenBanks used for the scans live in `~/Downloads` (Addgene downloads for #216077,
  #236083 among others). The construct-series scan in §4.10 used, specifically:
  `p003_pag-lenti_236083_esp3i_domesticated_efs_tfrc_bc_alfa_iresfragmid_wpre_puro_assembled.gb`
  (and its `pag-lenti-cmv-empty-puro` variant),
  `p004_pag-lenti_236083_paci_clai_site_wpre_only.gb`,
  `p006_pag_lenti_cmv_tfrc_ctermtag_ires_puro.gb`,
  `p007_pag_lenti_cmv_tfrc_v291_ires_puro.gb`.
  ⚠️ **These are local-only and not in git or S3** — the enzyme counts in §4.9/§4.10 are
  reproducible only on Becca's machine until the GenBanks are checked in somewhere. Worth
  fixing; plasmid maps are small and belong under version control.

## Method

Desk analysis, all reproducible from this folder:

1. **Barcode length** — Monte-Carlo estimate of the constraint-passing fraction, then
   birthday-collision math across candidate lengths and library sizes
   (`barcode_length_analysis.py`).
2. **Recombination risk** — Walton/CROPseq-multi per-bp rate applied to candidate spans,
   with the reported rate range propagated as error bars (`recombination_analysis.py`).
3. **Sequence/domestication checks** — translation, ORF scans, pairwise alignment and
   restriction-site scans of the GenBanks and supplied element sequences. The construct-
   series Type IIS scan and the two #236083 Esp3I edits are reproducible via
   `enzyme_scan.py` (§4.9, §4.10); the remaining one-off sequence findings are recorded
   inline above with positions so they can be re-derived.
4. **Synthesis cost / error model** — exact single-nt consequence spectrum over the 61
   sense codons, applied per ORF under each vendor-routing scenario
   (`synthesis_cost_analysis.py`).
5. **Gene set and scoping** — PR #130's `passes_canonical` / `passes_likely` re-run over
   the 5,130-gene deep-dive cohort in the read-only public D1 mirror; isoform, ortholog,
   rejected-gene and filter queries in the four `query_*` / `library_scoping_filters.py`
   scripts.
6. **Missense severity** — AlphaMissense (Zenodo 8208688) filtered to our 2,242 accessions
   and reduced to a per-gene damaging fraction, then swept across the defensible parameter
   range (`alphamissense_severity.py`; the download + filter recipe, and the two ways to
   get it wrong, are in that script's docstring).

## Outputs

- `barcode_length_collision.tsv` + `figures/s1e2_barcode_length.{svg,png}`
  (from `barcode_length_analysis.py`)
- `recombination_risk_by_span.tsv` + `figures/s1e2_recombination_risk.{svg,png}`
  (from `recombination_analysis.py`)
- `synthesis_cost_by_scenario.tsv` + `figures/s1e2_synthesis_cost.{svg,png}`
  (from `synthesis_cost_analysis.py`) — cost across the four vendor-routing scenarios,
  including the cassette penalty and 2-variants-per-ORF, plus the length-concentration
  curve.
- `enzyme_scan.py` — **no data outputs; prints the §4.9/§4.10 tables.** Rescans the
  construct series for BbsI / Esp3I / BsaI / AarI / SapI (both strands, circular) and
  re-derives the two #236083 Esp3I domestication edits from the before/after windows,
  including the feature-overlap and no-new-site checks. ⚠️ Reads **local-only GenBanks**
  from `~/Downloads` — it reports missing files rather than failing, so on another machine
  it degrades to a no-op until the plasmid maps are shared.
- ⭐ **`barcodes_20nt_30k.tsv` — the production barcode set.** 30,000 × 20 nt, full published
  filter stack, **min Levenshtein 3 AND min Hamming 4**, Slutskin's 6-mer seed definition,
  emitted deterministically by `generate_barcodes.py` (SEED=20260812). 1,657,359 draws:
  98.15% rejected by constraints, 0.03% by Levenshtein, 0.004% by Hamming. GC 0.40–0.60
  (mean 0.506). Independently re-verified: all unique, zero violations across every filter
  including a 25-enzyme common-6-cutter survey, and an **exact** Hamming census (pigeonhole
  over 4 blocks of 5, so all 4.5 × 10⁸ pairs are covered) gives **minimum Hamming = 4, zero
  violating pairs**.
- `barcode_filter_stack.tsv` (from `barcode_filter_stack.py`) — each published barcode rule as
  a named predicate with its source, cumulative survival, and the minimum length that still
  reaches 30,000 at min Levenshtein 3.
- `umi_width.tsv` (from `barcode_design.py`) — UMI width vs items sharing one barcode.
- `barcode_design_capacity.tsv` + `figures/s1e2_barcode_design.{svg,png}`
  (from `barcode_design.py`) — miRNA-seed exclusion cost by conservation tier, edit-distance
  barcode capacity by length, and the ORF-prefix / clone-suffix budget.
- `mirna_seed_sites.tsv` — 2,064 human miRNA 7mer-m8 target sites (DNA), tier-labelled,
  derived from TargetScanHuman 8.0 `miR_Family_Info.txt`. ⚠️ Confirm redistribution terms.
- `strategy_comparison.tsv` + `figures/s1e2_strategy_comparison.{svg,png}`
  (from `strategy_comparison.py`) — cost-vs-quality frontiers for fragments+clonal vs
  AbstractBio+stitching, plus the flat-per-gene-fee (H) scenarios and the fee sweep.
- `alphamissense_by_gene.tsv` + `figures/s1e2_alphamissense_severity.{svg,png}`
  (from `alphamissense_severity.py`) — per-gene damaging fraction over **28,366,506
  variants / 2,235 genes**, plus the severity-parameter sensitivity sweep printed to
  stdout. ⚠️ Derived from **CC BY-NC-SA 4.0 (non-commercial)** AlphaMissense predictions.
- [`REPORT_library_scoping.md`](REPORT_library_scoping.md) — long-form appendix:
  rejected-gene audit (§1), scenario E stitching (§2), duplication gains (§3), the four
  filter lists in full (§4), the giants (§5), Twist clonal (§5b), Elegen Infinia (§5c), the
  AlphaMissense cutoff challenge (§5d), consolidated recommendation (§6). Mirrored as the
  APPENDIX section of the Notion page.
- `canonical_topology.json` — the **2,243 canonical genes** (PR #130) with protein
  lengths + DeepTMHMM topology, pulled from the read-only public D1 mirror. This is the
  input to the cost/error model and the long-ORF QC.
- `canonical_accessions.{json,txt}` — gene → UniProt accession map / accession list; the
  `.txt` is the filter input for the AlphaMissense extraction.
- D1 query helpers: `query_rejected_anchored.py` (rejected-gene audit),
  `query_likely_tier.py` (the "likely" tier), `query_isoforms_orthologs.py` (isoform and
  mouse/cyno ortholog scope), `library_scoping_filters.py` (the four filter lists). Each
  carries a header documenting the two D1 gotchas — dedupe `topology_public` on
  `topology_version`, and join orthologs on `compara_ortholog_ecd.ortholog_uniprot_acc`,
  not gene symbol.

## Status

**design in progress** — architecture and all major element choices settled; blocking
sequence fixes (stops, linker diversification, barcode placement) outstanding before
synthesis.

---

## Decision log

Things that changed. Kept out of the body above so the body reads as current.

| Was | Now | Why |
|---|---|---|
| IRES → **tNGFR** (non-fluorescent, MACS-selectable, all channels free) | **IRES2 → mScarlet3 → P2A → PuroR** | Changed direction. Red channel now occupied; selection is drug + fluorescence, not antibody/MACS. **The PE anti-CD271 ME20.4 order (#6428) is likely moot.** |
| **DNA Chisel** whole-ORF codon optimization (TFRC done) | **No codon optimization of native ORFs**; minimal site-only edits | Co-translational folding risk (§4.2). **TFRC CDS needs re-deriving.** |
| Barcode **18 nt** | **20 nt** | ORF-BC convention; 18 was already sufficient. |
| IRES2 from **Addgene #100673** (EMMA) | **Fragmid `pF_AA333`** | #100673 carries **2 BbsI sites** → would be cut during assembly. |
| Parent **#236083** | **pRDA_722 / #216077** (domesticated build) | #236083 is undomesticated (2 BbsI + 3 Esp3I), incl. one **inside the WPRE** and one **inside the 3′LTR**. #236083 retained for p006/p007. |
| Feature map annotated **T2A** | **P2A** (Fragmid `P2A_v1.1`) | Mis-annotation; P2A also correct on the merits (§4.6). |
| "**WPRE3** (X-protein-inactivated)" | **`mut6`** / Schambach ORF-free PRE | **That label was wrong** — unverified, and WPRE3 appears to be a separate *shortened* variant, a different modification. |
| "Fragmid **empirically validated** `IRES2_v1`" | **Retracted** | **Unsupported inference.** McGee et al. contains no functional data ("IRES2" appears 0 times). |
| Recombination span **150 bp → 7.6%** | **177 bp → 8.96%** | Corrected span definition (linker start, not tag-cassette start). |
| ⭐ **F3** — Twist clonal >2.27 kb, $276,416 / 89.0% | ⭐ **S7** — AbstractBio ≤1.2 kb + **Synplexity $50/gene flat to 3 kb** + clonal >3 kb, **$221,446 / 99.1%** | Synplexity's flat fee covers **84.6% of the library (1,897 ORFs) for $94,850, sequence-perfect**. S7 is *cheaper than the original D baseline* at 99.1% usable vs 79.3%. Contingent on confirming the terms (§5.4). |
| Flat $50/gene modelled with **no length cap** ($55,897 / 90.9%; "$102,267 at 99.9%") | **3 kb cap applied** → S-series | My model assumed $50 would buy a 43 kb gene. With the real cap, 346 ORFs (15.4%) fall outside it and carry **all** the residual cost and error. |
| "Puro resistance proves the ORF cistron was **translated**" (§4.3, first version) | It proves the **transcript exists and is intact through the IRES** — nothing more | **My error.** IRES2 initiates internally and cap-independently, so an ORF frameshift does *not* block PuroR. Only **2A** gives the strong filter, and 2A is unavailable because the ORF needs a free C-terminus. Materially weakens the case for staying IRES-linked (§5.12). |
| "#236083 is the **non-domesticated** background" (implying it was never domesticated) | **Two #236083 lineages** — p003/p004 **were** hand-domesticated to BbsI 0 / Esp3I 0; p006/p007 were not | The doc only described p006/p007, so it read as though #236083 was left undomesticated wholesale. The p003/p004 line and its two Esp3I edits (252, 306) are now recorded in §4.9/§4.10. |
| Clonal genes cost "**2–4× fragments**", so prohibitive at scale | **1.29×** ($0.09 vs $0.07/bp) | My assumption was wrong. At 1.29× clonal is affordable, and above ~2.3 kb it is **cheaper per *usable* construct** than fragments. |
| Mitigate long-ORF errors by **synthesising twice (2× builds)** | **Twist clonal >~2.3 kb (F3)** | Clonal **strictly dominates** duplication on the 870-ORF long tail: $272,675 @ 100% usable vs $424,161 @ 89.5%. Cheaper *and* perfect. §5.8 retained for the record. |
| Missense-severity **flat 0.32** (genome-wide) | **0.427 per-gene**, measured on our 2,235 genes | AlphaMissense run over the actual gene set; our surfaceome is markedly less missense-tolerant. Severe rate 13.3% → 14.8%. |
| AlphaMissense **0.427 treated as settled** | **Quoted as a 10–16% band** | The `pathogenic` cutoff is tuned for ClinVar precision, and clinical pathogenicity ≠ surface display (Taipale-chained ≈ 0.07). Practically minor — the indel term dominates. `REPORT_library_scoping.md` §5d. |

### Superseded: tNGFR vs synthetic Fc–PDGFRβTM marker comparison

Retained because **if we ever need channel-free enrichment again, tNGFR is the documented
option.** The alternative considered was a synthetic display cassette: **Igκ leader +
human IgG1 Fc + PDGFRβ TM** (pDisplay-style anchor).

| | **tNGFR (ΔLNGFR / CD271)** | **Igκ-Fc-PDGFRβTM** |
|---|---|---|
| Topology / SP | Type I, native SP; live-surface stain | Type I; Igκ leader; live-surface stain |
| Detection | anti-CD271 **+ GMP MACS beads** | anti-human-IgG-Fc **+ Protein A/G**; Fc dimeric → bright, high-avidity |
| Bonus | none (inert by design) | Fc doubles as affinity/purification handle |
| **Nonspecific binding** | minimal — no Fc-receptor engagement | **binds FcγR / FcRn / Protein A** → background on immune/myeloid cells, elevated in high-IgG tissues; serum IgG competes |
| **In-vivo liability** | clean | **FcγR → ADCC/phagocytosis** unless effector-silenced (LALA-PG) |
| Clinical selection | established GMP CD271 MACS | less standard |
| Endogenous background | some cells express CD271 (neural crest, MSC, some tumours) | serum/tissue IgG + FcR⁺ cells |

Verdict at the time: **tNGFR won on cleanliness** for delivery/in-vivo work; Fc–PDGFRβTM
is excellent **in vitro** and a reasonable fallback *if* an effector-silenced Fc is used.

### tNGFR sequence notes (superseded, retained for reference)

The tNGFR CDS then in hand was self-contained: **Kozak + `ATG` + a functional SP** (`M-Y`
fused to the native NGFR hydrophobic core/cleavage site, res 11–28), cleaving to the
correct mature ectodomain — **no additional SP needed** (and don't stack one). Full NGFR
ectodomain + TM, truncated after the TM (`…KRWN`) → signaling-dead. Detection reagent was
PE anti-human **CD271, clone ME20.4** (BioLegend #6428), with Miltenyi CD271 MicroBeads
for channel-free MACS. tNGFR residues **11–276** — the same boundaries as pOT_2 (#181971).

---

## Notes

### Linker–SpyTag003–ALFA codon check

The tag cassette (GS–SpyTag003–GS–ALFA) was run through **Twist Bioscience's online codon
optimizer**.

- **SpyTag003:** amino-acid sequence matches the **Howarth reference**
  ([Addgene #133450](https://www.addgene.org/133450/), SpyTag003-MBP; Keeble et al. 2019
  *PNAS*). Our copy is the **same protein, codon-optimized for mammalian expression via
  Twist** (the reference construct is bacterial, pET28a/T7).
- **ALFA:** amino-acid sequence matches **Leo's reference** (`PSRLEEELRRRLTEP`); the only
  differences are **synonymous** — protein identical:

```
codon:   1    2    3    4    5    6    7    8    9   10   11   12   13   14   15
paper:  CCC  AGC  AGG  CTG  GAG  GAA  GAG  CTG  AGG  CGC  AGA  CTG  ACC  GAG  CCC
ours:   CCC  AGC  AGG  CTG  GAG  GAG  GAG  CTG  AGG  AGG  AGG  CTG  ACA  GAG  CCT
AA:      P    S    R    L    E    E    E    L    R    R    R    L    T    E    P
diff:                        6            10   11        13        15
```

Differing codons 6, 10, 11, 13, 15 are all synonymous (E: GAA/GAG · R: CGC,AGA/AGG ·
T: ACC/ACA · P: CCC/CCT) → ALFA epitope unchanged.

### Kozak context detail

`gaccggttctagagccgccacc` = **AgeI** (`accggt`, pos 2) + **XbaI** (`tctaga`, pos 8) +
**Kozak** (`gccgccacc`). Evaluated as `gccgccaccATGTAC`: **−3 = `A`** (purine ✅ — the
single most important Kozak position); **+4 = `T`** (⚠️ suboptimal; canonical strong Kozak
is `gccRccATGG`). Context `GCCACCATGT` vs ideal `GCCACCATGG`. Net: strong where it matters,
minor sub-optimality at +4 — worth noting rather than fixing, since a +4 `G` would change
the second residue (and native TfR1 `ATGATGGAT` gives +4 = `A`, also not ideal).

---

## References

**Fragmid / toolkit**

- McGee et al., *Modular vector assembly enables rapid assessment of emerging CRISPR
  technologies*, Cell Genomics 2024. doi:10.1016/j.xgen.2024.100519
- `pF_AA197` [#216049](https://www.addgene.org/216049/) · `AA333`/`pF_AA333`
  [#216067](https://www.addgene.org/216067/) · `AA014`
  [#216057](https://www.addgene.org/216057/) · `pRDA_722`
  [#216077](https://www.addgene.org/216077/)
- `pOT_2` [#181971](https://www.addgene.org/181971/) (Legut et al., *Nature* 2022) ·
  `pAG Lenti CMV Empty Puro` [#236083](https://www.addgene.org/236083/)
- **PR #10** — `deliverome_analysis.codon_opt` (barcode generator; the whole-ORF recoder is
  superseded for native ORFs): https://github.com/Deliverome-Project/deliverome-analysis/pull/10

**Tags / reporters**

- Keeble et al., PNAS 2019 (PMID 31822621) — SpyTag003/SpyCatcher003; Addgene #133450
- Vicidomini et al., Nat Commun 2024. doi:[10.1038/s41467-024-50462-2](https://doi.org/10.1038/s41467-024-50462-2)
  — ALFA proline flanks
- Lenaerts et al., *EndoNB*, bioRxiv 2025. doi:[10.1101/2025.06.08.658482](https://doi.org/10.1101/2025.06.08.658482)
  — internal ALFA tagging of TFRC (V291), basis for p007
- Gadella et al., *mScarlet3*, Nat Methods 2023;20(4):541-545.
  doi:[10.1038/s41592-023-01809-y](https://doi.org/10.1038/s41592-023-01809-y) ·
  [FPbase](https://www.fpbase.org/protein/mscarlet3/)
- Addgene [#61505](https://www.addgene.org/61505/) — TfR–pHuji (Shen et al. 2014 *JCB*)

**Barcoding / screens**

- Tycko et al., *High-Throughput Discovery and Characterization of Human Transcriptional
  Effectors*, Cell 2020;183(7):2020-2035.e16.
  doi:[10.1016/j.cell.2020.11.024](https://doi.org/10.1016/j.cell.2020.11.024) —
  **HT-recruit**; UTR handle + barcode 5′ primer
- Feldman et al., *Optical Pooled Screens in Human Cells*, Cell 2019;179(3):787-799.e17.
  doi:[10.1016/j.cell.2019.09.016](https://doi.org/10.1016/j.cell.2019.09.016) —
  **lentiGuide-BC**; barcode 3′ primer
- Feldman et al., Nat Protoc 2022;17(2):476-512.
  doi:[10.1038/s41596-021-00653-8](https://doi.org/10.1038/s41596-021-00653-8) — protocol
- Walton, Qin & Blainey, *CROPseq-multi*,
  [bioRxiv 2024.03.17.585235](https://www.biorxiv.org/content/10.1101/2024.03.17.585235v1.full)
  / [PMC10979941](https://pmc.ncbi.nlm.nih.gov/articles/PMC10979941/) —
  template-switching rate model

**Vector reuse / shuttle architecture (§4.11)**

- **Tedman, Goel, Shah, Howard, Chamness, Bonifasi, Adams, Gallagher, Penn, Nemec,
  McDonald, Corman, Robinson, Post, Clark, Babu, Manglik, Kuntz, Coyote-Maestas &
  Schlebach**, *Deep Receptor Scanning Reveals General Sequence Constraints on GPCR
  Biosynthesis*, bioRxiv 2025 (posted 21 Sep 2025).
  doi:[10.1101/2025.09.19.677468](https://doi.org/10.1101/2025.09.19.677468) ·
  [PMC12458215](https://pmc.ncbi.nlm.nih.gov/articles/PMC12458215/) — **BsmBI (not BsaI)**
  sites upstream of the UMI and immediately 3′ of the ORF for single-reaction
  backbone-to-backbone shuttling; 766 GPCRs + 174 splice variants surface-expression scan.
  ✅ Verified against the PMC full text.

**Barcoded-ORF libraries with a decoupled marker (§4.3)** — ⚠️ specifics user-supplied,
citations below are my attribution and need confirming:

- **MORF** — Joung et al., *A transcription factor atlas of directed differentiation*,
  Cell 2023 (Zhang lab). 3,548 TF isoforms / 1,836 genes; `EF1α → TF ORF → [V5]` plus
  `SV40 → PuroR` in Broad TRC vector **pLX_TRC317**; 24-bp barcode, dial-out-capable.
  ⚠️ Verify the citation and the vector ID.
- **REPROcode** — 408 immune-related TF constructs; `SFFV → TF ORF → WPRE → barcode`,
  transcribed barcode 3′ of WPRE; **no mammalian drug-resistance marker listed** for the
  deposited SFFV-MCS vectors. ⚠️ **No citation yet — needs one.**
- **Human SLC Overexpression Library** — **RESOLUTE** consortium (IMI-funded, CeMM /
  Superti-Furga). 467 barcoded, **codon-optimized** inserts; `TRE/TetO → SLC ORF` plus
  `PGK → rtTA-VP16–2A–PuroR` in **pCW57.1**; barcode + UMI; dox-inducible.
  ⚠️ Verify the primary citation.

**Translation: codons, stops, 2A**

- Kimchi-Sarfaty et al., Science 2007;315(5811):525-8. doi:[10.1126/science.1135308](https://doi.org/10.1126/science.1135308)
- Buhr et al., Mol Cell 2016;61(3):341-351. doi:[10.1016/j.molcel.2016.01.008](https://doi.org/10.1016/j.molcel.2016.01.008)
- Walsh, Bowman, Soto Santarriaga, Rodriguez & Clark, PNAS 2020;117(7):3528-3534. doi:[10.1073/pnas.1907126117](https://doi.org/10.1073/pnas.1907126117)
- Bartoszewski et al., Cell Mol Biol Lett 2016;21:23. doi:[10.1186/s11658-016-0025-x](https://doi.org/10.1186/s11658-016-0025-x)
- Konczal, Bower & Gray, PLoS ONE 2019;14(4):e0215892. doi:[10.1371/journal.pone.0215892](https://doi.org/10.1371/journal.pone.0215892)
- McCaughan, Brown, Dalphin, Berry & Tate, PNAS 1995;92(12):5431-5. doi:[10.1073/pnas.92.12.5431](https://doi.org/10.1073/pnas.92.12.5431)
- Tate et al., Biochem Cell Biol 1995;73(11-12):1095-103. doi:[10.1139/o95-118](https://doi.org/10.1139/o95-118)
- Blanchet et al., NAR 2014;42(15):10061-72. doi:[10.1093/nar/gku663](https://doi.org/10.1093/nar/gku663)
- Lombardi et al., Adv Protein Chem Struct Biol 2022;132:1-48. doi:[10.1016/bs.apcsb.2022.06.001](https://doi.org/10.1016/bs.apcsb.2022.06.001)
- Kim et al., PLoS ONE 2011;6(4):e18556. doi:[10.1371/journal.pone.0018556](https://doi.org/10.1371/journal.pone.0018556)
  — **P2A > T2A/E2A/F2A**

**WPRE**

- Zanta-Boussif et al., Gene Ther 2009;16(5):605-19. doi:[10.1038/gt.2009.3](https://doi.org/10.1038/gt.2009.3)
  — **`mut6`**; native and mutant "functioned equivalently"
- Schambach et al., Gene Ther 2006;13(7):641-5. doi:[10.1038/sj.gt.3302698](https://doi.org/10.1038/sj.gt.3302698)
  — ORF-free PRE (no ORF >25 aa)

**Other**

- DNA Chisel — [Nat Commun s41467-026-72227-9](https://www.nature.com/articles/s41467-026-72227-9#Sec30)
  *(no longer used for native ORFs — see Decision log)*
- Yang & Wang, J Mol Cell Biol 2019;11(10):911-919. doi:[10.1093/jmcb/mjz091](https://doi.org/10.1093/jmcb/mjz091)
  — cited by Fragmid for its IRES; **note this is a review, not a sequence-of-record**
- `YCe1760 HC_Kan_IRES2_p8b` [#100673](https://www.addgene.org/100673/) (Martella et al.,
  ACS Synth Biol 2017, EMMA) — the undomesticated IRES2
