# s1e2 Appendix — Library scoping, synthesis strategy and rejected-gene audit

Companion to [`README.md`](README.md) §5.4–5.11. All numbers from the read-only public
D1 mirror (accessible-surfaceome), gene set = the **2,243 canonical** genes from
[PR #130](https://github.com/Deliverome-Project/accessible-surfaceome/pull/130).
Vendor parameters as supplied: Twist pools **$0.00647/bp ≤1.8 kb, 1:3000**; Twist
fragments **$0.07/bp, 1:7500**; Abstract Bio **$0.01/bp ≤~1.2 kb, 1:20k–1:40k**.
Baseline throughout is **scenario D = $224,478** (AbstractBio ≤1.2 k + pools 1.2–1.8 k +
fragments >1.8 k), **ORF-only — no tags, barcode or linkers**.

---

## 1. Rejected genes with a membrane anchor

Of the **5,130** deep-dived genes: **2,243 canonical**, **625 likely-only**,
**2,262 rejected**.

**1,129 of the rejected genes carry ≥1 TM helix or a GPI anchor** (1,111 TM, 19 GPI).
That sounds alarming, but the rejection reasons are mostly sound:

| failing gate | n |
|---|---|
| `confidence = low` | 913 |
| `surface_accessibility = no` | 292 |
| `surface_accessibility = uncertain` | 33 |
| `evidence_grade = conflicting` | 27 |

By specificity: **528 `mostly_intracellular`**, 356 `surface_dominant`, 245 `mixed`.

### The 86 high-confidence rejections are correctly rejected

Every one I can identify is a **transmembrane protein in the wrong membrane** — which is
exactly what a *surface* catalogue should exclude:

- **Nuclear envelope / pore:** NUP210, NUP210L, NDC1, SYNE2
- **ER / ERGIC:** RTN2, RTN3, ESYT2, ESYT3, GRAMD1A/B/C, SSR1, TMX4, ATP13A1, HM13, WFS1
- **γ-secretase (ER/Golgi):** PSEN1, PSEN2, PSENEN, NCSTN
- **Vesicle / SNARE:** VAMP2, VAMP8, STX1A, SV2B, SYT6, SYT11
- **Lysosome / endosome:** NPC1, SLC11A1, SLC15A4, LMBRD1
- **Golgi / trans-Golgi:** ATP7B, CASD1, SLC33A1, B4GALNT2, FUT4
- **Ciliary / other:** TMEM231, PRCD, ZDHHC20, PIGT

This is a **good validation signal** — the pipeline is correctly distinguishing
"has a TM helix" from "is at the plasma membrane."

### ⚠️ The 38 that ARE concerning

`confidence ∈ {high, moderate}` **and** `surface_specificity ∈ {surface_dominant, mixed}`
**and** `surface_accessibility ∈ {high, moderate}` — i.e. the deep dive itself thinks
they're accessible surface proteins, yet they were rejected. **Every one fails on
`evidence_grade = weak` with mostly `evidence_density = low`** — a *corpus-coverage*
failure, not a biology failure.

The list contains well-established surface proteins:

| gene | why it stands out |
|---|---|
| **ANTXR2** | anthrax toxin receptor 2 — an **established toxin-delivery receptor** |
| **EPHB1**, **EFNB3**, **EFNA2** | canonical Eph/ephrin surface receptor-ligand pairs |
| **TNFRSF11A** | **RANK** — major TNF-receptor family member |
| **LCT**, **SI** | lactase, sucrase-isomaltase — classic intestinal brush-border surface enzymes |
| **CPM** | carboxypeptidase M, GPI-anchored |
| **ART4** | Dombrock blood group, GPI-anchored |
| **XG** | Xg blood group antigen |
| **GP5** | platelet glycoprotein V |
| **EDAR** | ectodysplasin A receptor |
| others | SLAMF8, CLMP, BTN1A1, PLXNA3, CNTN3, LRRC4B, CDH22, PCDHA1/A11/B15, GABRR3, CACNG4/6, KCNH8, ENPP6, ADAM33, ABCB5, EVI2A |

**Recommendation:** treat these 38 as a **manual-review shortlist**, not as rejects.
They look like genes the literature corpus under-covered rather than genes that aren't at
the surface. **ANTXR2, EPHB1, TNFRSF11A and the ephrins would be odd omissions from a
delivery-receptor library.** This is also indirect support for issue #131 (fixing
`evidence_grade` at source) — the same weakness that PR #130 worked around for canonical
is still suppressing these.

---

## 2. Scenario E — AbstractBio-grade pieces + Genyro/Sidewinder stitching

**Model:** split every ORF >1.2 kb into ≤1.2 kb pieces synthesised at Abstract Bio
($0.01/bp, **1:30,000**), then stitch. 4,854 pieces total across the 2,243 ORFs. The
assembly cost per piece is **unknown**, so it's swept.

| assembly cost | total cost | usable (no protein-altering error) | severe-error rate |
|---|---|---|---|
| **scenario D reference** | **$224,478** | **79.3%** | **13.3%** |
| E @ $0/piece | $45,576 | **94.7%** | **3.2%** |
| E @ $44/piece | ~$224,478 | **94.7%** | **3.2%** |
| E @ $50/piece | $247,076 | 94.7% | 3.2% |
| E @ $100/piece | $448,576 | 94.7% | 3.2% |
| E @ $200/piece | $851,576 | 94.7% | 3.2% |

**Break-even vs scenario D is $44 per piece.** Below that, **stitching is cheaper *and*
dramatically cleaner** — 94.7% usable vs 79.3%, and the severe-error rate drops
**13.3% → 3.2%** (a 4× reduction).

**This is the single highest-leverage open question in the whole plan.** The prize is
+15 points of usable library and a 4× cut in severe errors. **Get Genyro / Sidewinder
quotes on a per-piece / per-assembly basis** — that one number decides the entire
synthesis strategy. Even at $50–100/piece, E is arguably worth it for the quality alone.

⚠️ Caveats: this model assumes the assembled molecule inherits the AbstractBio per-bp
error rate and ignores (a) junction fidelity, (b) whether pieces are individually
sequence-verified before assembly (if they are, effective error is far lower still), and
(c) any per-project setup fee. Treat as an upper bound on quality / lower bound on cost
until quoted.

---

## 3. Duplication — the gain even when it doesn't reach the target

Duplication helps materially well below the 95% bar. **P(at least one of N builds has no
severe error)**, Twist fragments 1:7500:

| CDS | 1 build | 2 builds | 3 builds |
|---|---|---|---|
| 2,000 | 87.8% | **98.5%** | 99.8% |
| 3,000 | 82.2% | **96.8%** | 99.4% |
| 4,000 | 77.0% | **94.7%** | 98.8% |
| 6,000 | 67.6% | 89.5% | 96.6% |
| 8,000 | 59.3% | 83.5% | 93.3% |
| 10,000 | 52.1% | 77.0% | 89.0% |
| 15,000 | 37.6% | 61.0% | 75.7% |
| 20,000 | 27.1% | 46.9% | 61.3% |
| **43,524 (MUC16)** | **5.8%** | 11.4% | 16.5% |

Same, for **no non-synonymous error at all**:

| CDS | 1 build | 2 builds | 3 builds |
|---|---|---|---|
| 2,000 | 80.2% | **96.1%** | 99.2% |
| 4,000 | 64.3% | 87.2% | **95.4%** |
| 6,000 | 51.5% | 76.5% | 88.6% |
| 10,000 | 33.1% | 55.2% | 70.1% |
| 20,000 | 11.0% | 20.7% | 29.4% |
| 43,524 | 0.8% | 1.6% | 2.4% |

**Reading it:** a second build buys **+10 to +25 points** across the 2–10 kb range —
worth having even where 95% is out of reach. But it is **hopeless in the extreme tail**:
MUC16 goes 5.8% → 16.5% at 3× the cost. For anything >15 kb, duplication is not the
answer — truncation, a short isoform, or stitching is.

---

## 4. The filter lists in full

### 4a. The 79 secreted-only, non-GPI-anchored ORFs — $11,342

Families: 53 miscellaneous, 21 enzyme, 5 receptor. **This filter spares the 116
GPI-anchored proteins** that DeepTMHMM also labels `SP`.

`HSPG2(13176) PKHD1L1(12732) CUBN(10872) LAMA3(10002) COL6A3(9534) COL7A1(8835)
FN1(7434) TNC(6606) AGRN(6207) STRC(5328) C3(4992) PAPPA(4884) PRG4(4215)
ADAMTS16(3675) THBS1(3513) CACNA2D2(3453) EGFLAM(3054) SCUBE2(3000) SCUBE3(2982)
SCUBE1(2967) THBS4(2886) MMRN2(2850) ADAMTS5(2793) SULF1(2616) SULF2(2613) LOXL2(2325)
SEMA3A(2316) LOXL4(2271) SEMA3B(2250) MUC20(2130) HHIP(2103) PCSK9(2079) HPSE2(1779)
CR1L(1710) NPNT(1698) LGI1(1674) GLDN(1656) LGI2(1638) CDSN(1590) LIPG(1503)
NPTXR(1503) EDIL3(1443) HTRA1(1443) ERVV-1(1434) LPL(1428) ENPP7(1377) FGL2(1320)
FBLN7(1320) ASTL(1296) PLAU(1296) MFGE8(1164) ACP3(1161) HAPLN1(1065) SPESP1(1053)
APOH(1038) UPK3B(963) FSTL1(927) OSCAR(849) KLK2(786) BPIFA1(771) CRISP2(732)
SOD3(723) TIMP4(675) APCS(672) TIMP2(663) CLEC4D(648) IGLL1(642) PLAC1(639) TIMP3(636)
IL15(489) LY96(483) ECRG4(447) LYZL6(447) VPREB1(438) PATE1(381) DEFB118(372)
B2M(360) DEFB126(336) CXCL12(282)`

⚠️ **Three judgement calls in here:**
- **B2M** is β2-microglobulin — it is surface-*presented* as the obligate light chain of
  MHC-I, not simply secreted. Dropping it is defensible for a single-gene screen but note
  it isn't "just soluble."
- **PCSK9** is secreted, but it acts *on* a surface receptor (LDLR) and is itself a
  landed therapeutic target — may be worth keeping for other reasons.
- **CACNA2D2** is the α2δ auxiliary subunit of a Ca channel — GPI-anchored in some
  annotations; worth double-checking before dropping.

### 4b. The 229 requiring a co-receptor for surface expression — $20,514

Families: 87 receptor, 66 transporter, 59 misc, 17 enzyme. These are **obligate
multi-subunit complexes** — this is the highest-value filter conceptually because they
cannot reach the surface from single-gene overexpression:

- **Integrins** — the entire ITGA*/ITGB* set (ITGA1/3/4/9/10/11/2B/AE/AL/AM/AV/AX/AD, ITGB1/3/5/6/8)
- **GABA-A subunits** — GABRA3/A5, GABRB1/B2/B3, GABRG1/G2/G3, GABRD, GLRB, GABBR1/2
- **Nicotinic ACh subunits** — CHRNA1/2/3/5/10, CHRNB1/B2/B4, CHRND, CHRNE, CHRNG, CHRFAM7A
- **HLA class I & II** — HLA-B/C/G, HLA-DRA/DRB1/DRB3/DRB4/DRB5, HLA-DQA1/A2/B1/B2, HLA-DPA1/DPB1, CD74, B2M
- **CD3 / TCR complex** — CD3D/E/G, CD247, TRA, TRB, CD8B, CD79A/B
- **Na⁺/K⁺-ATPase** — ATP1A1/A2/A3, ATP1B1/B2/B3, ATP4B, ATP12A
- **CATSPER complex** — CATSPER1/2/3/4, CATSPERD/E/G/Q/T
- **Sarcoglycans** — SGCA/B/D/G, SSPN
- **Rh blood group** — RHD, RHCE, RHAG
- **Heteromeric amino-acid transporters** — SLC7A5/6/7/8/9/10/11/13, SLC3A1
- **LRRC8 volume channels**, **TREM1/2 + TYROBP/HCST adaptors**, **RAMP1/2/3**, **CD1A–D**,
  **NMDA/kainate subunits** (GRIN1/2A/2B/2C, GRIK5), **PKD1/PKD1L1/PKD1L3**, **CUBN/AMN**,
  **GP1BA/GP1BB/GP9**, **DUOX1/2**, **ABCC8/9**, **ABCG5/8**, **SCNN1B/G**, **OTOP1/2/3**

⚠️ **Do not drop these blindly.** Integrins and HLA are major delivery/adhesion targets.
The better move for high-value families is **co-express the obligate partner** (e.g. ITGB1
with the α subunits; B2M with HLA-I) rather than exclude them. Suggest splitting this
filter: *exclude* the ones you won't chase (CATSPER, sarcoglycan, taste, olfactory-adjacent)
and *rescue with a partner* the ones you care about (integrins, HLA, NMDA).

### 4c. The 16 with no accessible ECD — $1,256

`CATSPERT(5463) DIP2A(4716) PITPNM3(2925) ENOX1(1932) ENOX2(1833) ITPRIPL1(1668)
SPAG8(1458) AAMP(1305) NEU3(1287) PLSCR1(957) EBAG9(642) CCDC25(627) CAV1(537)
CAV3(456) S100A14(315) MYMX(255)`

Clean call — **CAV1/CAV3 are caveolins** (cytoplasmic-face hairpin, no ECD),
**PLSCR1** is a scramblase, **S100A14** is cytosolic-family. None present an
extracellular epitope to engage. Cheap to remove, and removing them avoids
guaranteed-negative wells.

### 4d. The 51 olfactory / taste receptors — $984

**All 51 are GPCRs**: 3 TAS1R (sweet/umami, 2.5 kb each), 27 TAS2R (bitter), 21 OR
(olfactory). Sizes 900–1,017 bp except the TAS1Rs.

`TAS1R3 TAS1R1 TAS1R2 | TAS2R1/4/5/7/8/9/10/13/14/19/20/30/31/36/38/39/40/43/45/46/50 |
OR11H4 OR56B1 OR5A2 OR14J1 OR2AT4 OR51E2 OR10H1 OR51E1 OR2B11 OR2AG2 OR2H1 OR3A1 OR2W3
OR4M1 OR1D2 OR2J2 OR2L13 OR2H2 OR7D4 OR14I1 OR4S2 OR2A7 OR4D1 OR4C13 OR2D2 OR5K1 OR8D1`

**Why they're a poor fit:** ORs are notoriously **retained in the ER** in heterologous
cells without accessory factors — RTP1/RTP2 chaperones are required for surface
trafficking (and RTP1 is itself in the canonical set). TAS2Rs have the same problem.
Without co-expressing the chaperones these will read as non-expressors regardless of
biology.

**But note they're only $984 total** — 0.4% of budget. The argument for dropping them is
**interpretability** (avoid 51 guaranteed-ambiguous wells), not cost. If you have plate
space, keeping them with RTP1/RTP2 co-expression would actually be an interesting
sub-experiment. **OR51E1/OR51E2 in particular are ectopically expressed in prostate
cancer** and have been pursued as targets — those two may be worth keeping regardless.

---

## 5. The giants — are they worth it?

**103 ORFs >5 kb = 886,773 bp = $62,074 (27.7% of budget).** But the framing "giants are
expensive" is only half right: at **~$600/gene** they're *cheaper per gene* than their
reputation. The real problem is **feasibility and error rate**, not price.

| group | n | bp | cost |
|---|---|---|---|
| **delivery / uptake receptors (high value)** | 14 | 152,388 | **$10,667** |
| **NOTCH family** | 4 | 28,062 | $1,964 |
| giant adhesion (FAT/DCHS/CELSR/TENM/CDH23/PCDH15/DSCAM/SDK/CSMD) | 19 | 179,931 | $12,595 |
| voltage-gated Ca channels (CACNA1×10) | 10 | 67,059 | $4,694 |
| voltage-gated Na channels (SCN×9) | 9 | 52,770 | $3,694 |
| TRP / PIEZO / NALCN | 6 | 37,908 | $2,654 |
| PTPR phosphatases | 7 | 44,199 | $3,094 |
| plexins | 6 | 34,809 | $2,437 |
| mucins | 5 | 96,564 | $6,759 |
| ECM / secreted-type | 7 | 61,794 | $4,326 |
| everything else | 16 | 131,289 | $9,190 |
| **TOTAL** | **103** | **886,773** | **$62,074** |

### The verdict: don't drop giants wholesale — drop by category

**Keep (high value, $12,632 for 18 genes):** **LRP1** (the canonical BBB/delivery
receptor), **LRP2/megalin**, **LRP1B**, **LRP4**, **SORL1**, **LY75/DEC-205** (classic
antigen-*delivery* receptor), **SIGLEC1/CD169** (macrophage uptake), **CSPG4/NG2**
(established ADC/CAR target), **STAB1/STAB2** (scavenger/clearance receptors), **CR1**,
**CUBN**, **MUC16/CA125** (clinical ADC target), **ITGB4**, **NOTCH1–4**.
**Omitting LRP1 or LY75 from a delivery-receptor library would be a serious gap** — these
are precisely the positive controls.

**Strong candidates to drop ($22,042 for 44 genes):**
- **Ion-channel giants — 25 genes, $11,042 (18%).** CACNA1×10, SCN×9, TRPM/PIEZO/NALCN×6.
  Multi-pass, small ECDs, poor heterologous surface expression, and not delivery targets.
  Best cost:benefit removal in the set.
- **Giant adhesion — 19 genes, $12,595 (20%).** FAT1–4, DCHS1/2, CELSR1–3, TENM2/3/4,
  CDH23, PCDH15, DSCAM/DSCAML1, SDK1/2, CSMD1. Genuine large-ECD surface proteins, but
  low delivery relevance per dollar and among the least synthesisable.

**Handle specially — the extreme tail:** **MUC16 (43.5 kb)** and **MUC3B (40.4 kb)**
**cannot be ordered as fragments at all** and duplication is futile (5.8% → 16.5% at 3×).
Their bulk is repetitive PTS/tandem-repeat domain — a **recombination hazard as well as a
synthesis one**. For a binding/internalisation screen the membrane-proximal region + TM is
what matters. **Truncate, or use a short isoform** (§5.7 found MUC4 isoforms at
1,125–1,176 aa vs 5,412 canonical — 79% shorter and *more* tractable).

### Scenario E changes the giants calculus

At AbstractBio-grade error via stitching, usable rises 79.3% → 94.7% overall. Giants are
where a single long fragment fails worst (52% clean at 10 kb), so they benefit most. **If
Genyro/Sidewinder quote under ~$44/piece, most giants become viable and the case for
dropping them weakens considerably** — another reason to get that quote before cutting.

---

## 5b. Twist clonal genes ($0.09/bp, sequence-perfect) — this changes the recommendation

*(Reading "$.09c/gene" as **$0.09/bp**, consistent with fragments at $0.07/bp. Flag if
you meant something else.)*

**Clonal genes are only 1.29× the fragment price but are sequence-verified perfect.** That
inverts the arithmetic, because the right metric is **cost per *usable* construct**, not
cost per bp:

- fragments: `0.07 × L / q` where `q` = clean fraction
- clonal: `0.09 × L`, `q = 1`
- Equal when **q = 0.07/0.09 = 0.778**

| criterion | fragments drop below 77.8% clean at | → clonal wins above |
|---|---|---|
| no protein-altering error | **2,273 bp** | **~2.3 kb** |
| no severe error | **3,852 bp** | **~3.85 kb** |

**Above ~2.3 kb, clonal is cheaper per usable construct than fragments — despite costing
29% more per base.** Given our median CDS is 1,506 bp and 870 ORFs are >1.8 kb, this
applies to a large slice of the library.

### Mixed strategies with clonal

| strategy | cost | usable | severe | $/usable |
|---|---|---|---|---|
| **D** — no clonal | $224,478 | 79.3% | 13.3% | $126 |
| **F1** clonal for >5 kb (103 giants) | $242,213 | 81.9% | 11.4% | $132 |
| **F2** clonal for >3.85 kb (187) — *severe crossover* | $249,460 | 83.3% | 10.5% | $133 |
| **F3** clonal for >2.27 kb (656) — *non-syn crossover* | **$276,416** | **89.0%** | **7.0%** | $138 |
| **F4** clonal for all >1.8 kb (870) | $285,072 | **90.9%** | **5.8%** | $140 |

**F3 is the sweet spot:** **+$51,938 over D (+23%) buys +9.7 points of usable library and
cuts the severe-error rate 13.3% → 7.0% (nearly halved).**

### The decisive comparison — fixing the long tail (870 ORFs, 3.03 Mbp)

| option | cost | usable | **$ per usable ORF** |
|---|---|---|---|
| fragments, 1 build | $212,080 | 70.0% | $348 |
| **fragments, 2 builds (duplication)** | **$424,161** | 89.5% | **$545** ← worst value |
| **CLONAL, sequence-perfect** | **$272,675** | **100.0%** | **$313** |
| stitched AbstractBio @ $44/piece | $159,305 | 91.0% | **$201** |

**⚠️ Clonal strictly dominates duplication: cheaper *and* perfect.** $272,675 for 100%
clean vs **$424,161 for 89.5%**. **Duplication should come off the table entirely** if
clonal is available at $0.09/bp — it was only ever a workaround for not having a
guaranteed-clean option.

**Clonal vs stitching:** stitching is still cheapest per usable ORF ($201 vs $313) *if*
assembly really lands near $44/piece — but clonal needs **no new vendor qualification, no
junction-fidelity risk, and no assembly step**, and its quality is contractual rather than
modelled. Sensible plan: **default to clonal for >2.3 kb now**, and treat stitching as the
optimisation to pursue once Genyro/Sidewinder are quoted and qualified.

⚠️ **Unknown: Twist's clonal-gene length ceiling.** Standard clonal gene products
typically top out in the low single-digit kb. **MUC16 (43.5 kb) and MUC3B (40.4 kb) remain
unorderable by any route** — they still need truncation or a short isoform regardless of
vendor. Confirm the clonal max length before assuming F3/F4 covers the whole tail.

## 5c. Elegen Infinia ($0.14–0.20/bp, 1:70,000) — ⚠️ UNVERIFIED, per user

At **1:70,000**, a mean 2,032 bp ORF is **97.6% clean** — excellent. But the comparison
that matters is **effective cost per usable base**:

| option | $/bp | clean @2 kb | **effective $/bp per usable construct** |
|---|---|---|---|
| Twist fragments | $0.07 | 80.2% | $0.0873 |
| **Twist clonal** | **$0.09** | **100%** | **$0.0900** |
| Elegen @ low end | $0.14 | 97.6% | $0.1434 |
| Elegen @ high end | $0.20 | 97.6% | $0.2049 |

**Clonal at $0.09/bp is cheaper per usable construct than Elegen at every length** —
$0.0900 vs $0.1434+. Elegen does beat *fragments* above **6,269 bp** ($0.14) or
**9,496 bp** ($0.20), but clonal beats both.

**So Elegen's role is specifically: where clonal cannot go.** Its value is entirely
contingent on **Twist's clonal length ceiling** — the open unknown from §5b. Elegen is
positioned as a long-linear-DNA vendor, so if clonal caps out at ~3.5–5 kb, Elegen is the
natural next tier for the 103 giants.

| strategy | cost | usable |
|---|---|---|
| **F3** clonal >2.27 kb (no Elegen) | **$276,416** | **89.0%** |
| G1 F3 + Elegen $0.14 above a **5 kb** clonal cap (103 → Elegen) | $320,755 | 88.6% |
| G2 F3 + Elegen $0.20 above a 5 kb cap | $373,961 | 88.6% |
| G3 F3 + Elegen $0.14 above a **3.5 kb** cap (238 → Elegen) | $348,191 | 88.3% |

Note G1–G3 all cost **more** than F3 for **marginally lower** usable output — because
Elegen is pricier than clonal and isn't perfect. **They are not upgrades; they are what you
pay if clonal has a ceiling.** So: **establish the clonal ceiling first.** If clonal
reaches 5 kb, F3 is the answer and Elegen isn't needed. If it caps at 3.5 kb, G3 costs
+$72k over F3.

Also worth asking Elegen about **turnaround time** — that can justify a premium
independently of cost-per-base.

## 5d. ⚠️ Is the AlphaMissense cutoff correct? — No, I'm not confident, and here's why

Fair challenge. Three distinct problems, in roughly increasing importance:

**(1) The threshold is DeepMind's, and it's arbitrary for us.** I used the pre-computed
`am_class = pathogenic`, whose score cutoff was chosen to hit **90% precision on ClinVar**.
That's a clinical-classification target, not ours. The **14% "ambiguous" bucket is simply
unassigned** — pushing it into "pathogenic" moves our figure 0.427 → **0.567**.

**(2) Possible upward bias for TM proteins specifically.** AlphaMissense leans heavily on
**evolutionary conservation (MSA depth)** plus AlphaFold structure. TM helices are among
the most conserved regions in the proteome, so random substitutions there score as
pathogenic. That means the **TM (0.463) > SP+TM (0.392)** ordering I reported may partly be
a *conservation artifact* rather than evidence that TM proteins are harder to traffic. It's
consistent with Taipale's mechanism, but it is not independent evidence for it.

**(3) The endpoint is wrong for our purpose — this is the real problem.** AlphaMissense
predicts **clinical pathogenicity**. We need **"does the protein still reach the surface and
present its epitope."** These diverge in both directions: a variant that abolishes catalytic
activity but traffics normally is *pathogenic yet fine for a binding screen*; a mildly
destabilising variant may be *clinically benign yet reduce surface levels*.

The better-matched endpoint is Taipale's mislocalization: **~1/6 of pathogenic missense
mislocalize**. Chaining that gives **0.427 × ~0.17 ≈ 0.07** — i.e. **~7% of random missense
may actually prevent surface display, not 43%.** That's a **6× difference** in the
parameter.

### The saving grace: it barely moves the answer

Severe-error rate under scenario D across the whole plausible range:

| missense-severity parameter | severe-error rate | source |
|---|---|---|
| **0.070** | **10.2%** | Taipale-chained (0.43 × ~1/6 mislocalise) — best endpoint match |
| 0.245 | 12.4% | VAMP-seq experimental (PTEN + TPMT abundance) |
| 0.320 | 13.3% | AlphaMissense genome-wide `pathogenic` |
| **0.427** | **14.6%** | AlphaMissense on **our** genes, `pathogenic` class ← used |
| *per-gene, not uniform* | *14.8%* | *same data, each gene weighted by its own fraction* |
| 0.567 | 16.2% | AlphaMissense counting `ambiguous` as pathogenic too |

**A 6× swing in the parameter (0.07 → 0.567) moves the severe-error rate only
10.2% → 16.2%** — because the **indel term (30% of errors, every one frameshifting) and the
nonsense term dominate**, and neither depends on missense severity.

**Conclusion:** the cutoff question is methodologically real but **practically minor —
worth ±3 percentage points**. The **Twist indel:substitution ratio is a far more important
unknown**, because it drives the term that actually dominates. Recommend: keep 0.427 as
the headline (it's measured on our actual genes), quote **10–16%** as the honest band, and
chase the indel ratio rather than refining the missense cutoff.

## 6. Consolidated recommendation

0. **Use Twist clonal genes ($0.09/bp, sequence-perfect) for everything >~2.3 kb — F3,
   $276,416, 89.0% usable.** Strongest single move available *today*, needs no new vendor,
   and **retires duplication entirely** (clonal is both cheaper and perfect vs 2×
   fragments). **First action: establish Twist's clonal length ceiling** — it determines
   whether Elegen is needed at all (§5c).
0b. **Elegen Infinia only if clonal has a ceiling.** Clonal is cheaper per usable base at
   every length ($0.090 vs $0.143+), so Elegen is the *fallback for lengths clonal can't
   reach*, not an upgrade. G1/G3 cost +$44k/+$72k over F3 for slightly *lower* usable
   output. Also ask about turnaround — speed may justify the premium independently.
   ⚠️ Elegen's $0.14-0.20/bp and 1:70,000 are **unverified**.
1. **Get Genyro / Sidewinder per-piece quotes.** Break-even is **$44/piece**; below it,
   scenario E beats D on cost *and* delivers 94.7% vs 79.3% usable with a 4× lower severe
   error rate. Best long-run option, but needs qualification — pursue in parallel with
   clonal as the interim default.
2. **Add the 330 likely-tier anchored genes — $19,931.** Best marginal value available:
   +15% targets for +9% cost.
3. **Manually review the 38 concerning rejects** (ANTXR2, EPHB1, TNFRSF11A, EFNB3/A2, LCT,
   SI, CPM, ART4 …). They fail on evidence density, not biology.
4. **Apply the cheap, unambiguous filters:** 16 no-ECD ($1,256) and 51 OR/TAS ($984) —
   for interpretability more than cost. Consider keeping OR51E1/E2.
5. **Split the 229 co-receptor set** — exclude the ones you won't chase; **rescue
   integrins and HLA by co-expressing the partner** rather than dropping them.
6. **Drop giants by category, not wholesale.** Cut ion-channel (25, $11k) and giant
   adhesion (19, $12.6k); **keep the 18 high-value delivery receptors ($12.6k)**.
   Truncate MUC16/MUC3B.
7. **If staying on scenario D:** target **85% no-severe with up to 2 builds — $394,235
   (1.76×)**. Do *not* chase 95%-no-non-synonymous; it costs $628k and still misses 247
   ORFs.

### Still open

- **Twist indel : substitution ratio** (assumed 30%) — **the biggest real unknown.** It
  drives the term that dominates the severe-error rate, and matters far more than
  refining the missense cutoff (§5d).
- **Twist clonal length ceiling** — gates the whole F3 vs G1/G3 decision.
- ✅ **AlphaMissense — DONE.** 28,366,506 variants over 2,235/2,243 genes. Our surfaceome
  is **more** missense-intolerant than the proteome average: **0.427 damaging vs 0.32
  genome-wide** (44.0% of our variants are `pathogenic`). By topology: TM 0.463 > GLOB
  0.458 > SP 0.412 > SP+TM 0.392. Raises the severe rate 13.3% → **14.8%** per-gene (14.6%
  if 0.427 is applied uniformly). But **the cutoff
  and endpoint are questionable — see §5d; honest band is 10-16%.**
  ⚠️ **CC BY-NC-SA 4.0, non-commercial** — needs a licence decision before any commercial
  or published use.
- **Genyro / Sidewinder specs** — length limits, junction fidelity, whether pieces are
  sequence-verified pre-assembly.
