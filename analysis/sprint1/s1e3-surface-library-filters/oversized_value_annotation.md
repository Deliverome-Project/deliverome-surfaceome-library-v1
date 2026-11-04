# >5 kb ORFs — value judgment for a truncated / minimal construct

The 107 ORFs dropped by the 5 kb packaging cap, annotated with my (Claude) judgment on
whether they're worth including **in a shortened form** for a surface-display / delivery
library, plus a candidate truncation strategy. Calls are quick expert triage, not final.

- **H (high)** — validated / delivery-relevant target with a plausible, often literature-backed truncation.
- **M (medium)** — plausibly useful; modular architecture makes a minimal construct feasible; less validated.
- **L (low)** — niche / structural / sensory / enzyme, or multi-pass with little targetable surface → skip for now.

Literature precedents for a shorter expressible variant are gathered per-gene by dispatched
research subagents (see the "Precedents" column / follow-up section).

## H — include via a truncated construct (dispatch precedent search)

| gene | aa | tier | tumor | why valuable | truncation strategy |
|---|---|---|---|---|---|
| **MUC16** | 14,507 | T3 | ✓ | CA125 — validated ovarian ADC/CAR target | membrane-proximal SEA domain(s) + TM (drop the ~60-repeat bulk); MUC16c114 / MUC-CD precedent |
| **LRP1** | 4,544 | T1 | ✓ | **transcytosis / delivery receptor** (angiopep BBB); scavenger | **minireceptor** — one ligand-binding cluster + TM + cytoplasmic NPxY |
| **LRP2 / megalin** | 4,655 | T1 | ✓ | renal endocytic / delivery receptor | LDLR-family minireceptor (one cluster + TM) |
| **NOTCH1** | 2,555 | T1 | ✓ | pan-cancer target | NRR + TM ± intracellular (ΔEGF "NEXT"/ΔE constructs are standard) |
| **NOTCH3** | 2,321 | T1 | ✓ | CADASIL + cancer | same NRR-based truncation |
| **ROS1** | 2,347 | T1 | ✓ | RTK oncogene (also fusion-driven) | membrane-proximal ectodomain + TM + kinase, or ECD-truncated |
| **ITGB4** | 1,822 | T1 | ✓ | carcinoma integrin, very high expression | ectodomain + TM (pairs with α6); or minimal β4 |
| **LY75 / DEC-205** | 1,722 | T1 | ✓ | dendritic-cell **endocytic / vaccine-delivery** receptor | subset of C-type lectin domains + TM + endocytic tail |
| **CSPG4 / NG2** | 2,322 | T3 | ✓ | melanoma / GBM **CAR** target, TPM 7,554 | minimal targetable domain + TM (CAR-antigen minigenes exist) |
| **STAB1** | 2,570 | T1 | ✓ | scavenger / **endocytic clearance** receptor (delivery) | FAS1/EGF-domain subset + TM |
| **SORL1** | 2,214 | T1 | ✓ | LDLR-family; Alzheimer's; sortilin-related | VPS10 + one LDLR cluster + TM |
| **LRP4** | 1,905 | T1 | ✓ | agrin/MuSK; myasthenia autoantigen | LDLa β-propeller subset + TM |
| **PLXNB1** (plexin family) | 2,135 | T1 | ✓ | semaphorin receptors — tumor / angiogenesis | Sema domain + TM (represents PLXNA1–4/B2/B3/D1) |
| **PTPRZ1** (RPTP family) | 2,315 | T3 | ✓ | GBM receptor-tyrosine-phosphatase | CAH + FN3 subset + TM (represents PTPRF/S/D) |
| **SIGLEC1 / CD169** | 1,709 | T3 | ✓ | macrophage **endocytic** delivery target | V-set + a few C2 domains + TM |
| **PCDH15 / CDH23** | 1,955 / 3,354 | T1 | ✓/– | Usher syndrome — **gene-therapy** relevant | EC-repeat subset / dual-vector (mini-PCDH15 precedents) |

## M — plausibly useful (modular, truncatable), lower priority

| gene | aa | tumor | note |
|---|---|---|---|
| NOTCH2 | 2,471 | ✓ | as NOTCH1/3 |
| STAB2 | 2,551 | – | scavenger/endocytic like STAB1 |
| FAT1 / FAT2 / FAT4 | 4,349–4,981 | ✓ | atypical cadherins / Hippo; tumor suppressors; repeat-deletable |
| MUC4 / MUC12 / MUC17 / MUC3A | 3,323–5,478 | ✓ | tumor-associated mucins; membrane-proximal + TM like MUC16 |
| PLXNA1–4 / PLXNB2 / PLXNB3 / PLXND1 | 1,871–1,925 | ✓ | semaphorin receptors (see PLXNB1) |
| PTPRF / PTPRS / PTPRD | 1,907–1,948 | ✓ | receptor phosphatases; neural/tumor |
| CR1 / CD35 | 2,039 | – | complement receptor; immune; CCP-module deletable |
| PIEZO1 / PIEZO2 | 2,521 / 2,752 | ✓ | mechanosensors; hard to truncate (multi-pass) but high interest |
| USH2A | 5,202 | – | Usher gene therapy (dual-vector) |
| DSCAM / DSCAML1 / SDK1 / SDK2 | 2,012–2,213 | ✓/– | Ig-superfamily neural adhesion |
| KIDINS220 | 1,771 | – | neurotrophin-signaling scaffold |
| PKD1 | 4,303 | ✓ | polycystin; PKD disease |
| ABCA1 / ABCA7 | 2,146–2,261 | – | lipid transporters; not classic antigens |
| CELSR1–3 | 2,923–3,312 | ✓/– | adhesion GPCRs (huge; repeat-deletable) |

## L — skip for now (niche / structural / sensory / multi-pass with little targetable surface)

- **Voltage-gated ion channels** (multi-pass, small extracellular targetable surface, poor truncation): `CACNA1A/B/C/D/E/F/G/H/I/S` (Cav), `SCN1A/2A/3A/4A/5A/7A/8A/9A/10A/11A` (Nav), `TRPM3/6/7`, `PKD1L2`, `PKDREJ`. *(Possible exception: SCN9A — a pain target — but not a surface-display-friendly antigen.)*
- **Brush-border / secreted enzymes**: `LCT` (lactase), `SI` (sucrase-isomaltase), `MGAM/MGAM2` (maltase), `PCSK5` (convertase).
- **Very large structural adhesion / matrix**: `MUC3B` (TPM 0), `MUC22`, `LRTM3`, `ZAN` (sperm), `PKHD1` (fibrocystin), `FRAS1`, `CSMD1/2/3`, `DCHS1/2`, `TENM1/2/4` (teneurins), `MEGF8`, `MALRD1`, `KIAA1549L`, `FER1L5` (LL), `PTPRZ2`.
- **ADGRV1** (6,306 aa adhesion GPCR) — largely structural; low priority unless a specific Usher-2 rationale.

## Summary
Of the 107, **~15 rate H** (worth a truncated construct now — headlined by **LRP1** for delivery and **MUC16/CA125**, **CSPG4/NG2**, **NOTCH1/3** for oncology), **~25 M**, and the rest **L** (mostly ion channels + brush-border enzymes + giant structural adhesion, where a minimal surface construct isn't worth it).
