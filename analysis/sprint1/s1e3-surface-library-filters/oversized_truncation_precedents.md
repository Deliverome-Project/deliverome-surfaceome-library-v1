# Truncation precedents for the >5 kb targets

Literature precedents for expressing a **shortened, surface-displayed (or functional)** version of
each high-value oversized target, gathered by per-gene research subagents (PubMed + web). Each
entry: the construct used and whether it retained surface expression / function. These de-risk a
"truncated-construct" sub-library for the giants dropped by the 5 kb ORF cap.

**General recipe that recurs across every family:** keep signal/anchor + the membrane-proximal
functional module + TM (± minimal cytoplasmic sorting/endocytosis motifs); delete the repetitive
bulk. Membrane anchoring is the make-or-break variable.

---

## LRP1 — minireceptor (single ligand-binding cluster + TM + tail) ✅ strong
- **Obermoeller-McCormick et al. 2001, J Cell Sci 114:899** (PMID 11181173) — mLRP2/3/4 (one cluster + TM + tail) fold, furin-cleave, reach the surface, bind ligand. **Cluster IV (mLRP4) is the workhorse; avoid cluster I** (misfolds).
- **Li et al. 2000, J Biol Chem 275:17187** (PMID 10747918) — the mLRP4 minireceptor endocytoses; the **YXXL + distal di-leucine** (not the NPxY motifs) is the dominant internalization signal → keep YXXL in the tail.
- **Jen et al. 2010, J Cell Sci 123:246** (PMID 20048341) — stable mLRP4 in neurons: ~84% surface+, binds RAP ~30×, endocytoses cargo to lysosomes.
- Caveat: correct cluster folding needed **RAP (LRPAP1) co-expression** → co-express or validate folding. Collapses ~13.6 kb → ~1.5–2 kb.

## LRP2 / megalin, LRP4, SORL1 — LDLR-family ✅ (megalin strong; SORL1 has a gap)
- **LRP2/megalin:** Takeda 2003 (PMID 12519751) + Nagai 2003 (PMID 14528014) + Marzolo 2003 (PMID 12694565) — a **one-cluster megalin minireceptor** (4th cluster + TM + tail) surface-expresses (biotinylation) and endocytoses ¹²⁵I-lactoferrin.
- **LRP4:** Zhang 2011 (PMID 21969364) + structure Zong 2012 (PMID 22302937) — N-terminal **LDLa + first β-propeller** fragment binds agrin and drives the surface LRP4–MuSK complex without TM/ICD.
- **SORL1/SorLA:** Jacobsen 2001 (PMID 11294867) — mature ectodomain surface-expresses + endocytoses (~10% surface); Andersen 2006 (PMID 16489755) — **CR-cluster** is the minimal APP-binding module. *Gap: no single published membrane-anchored VPS10+cluster minireceptor — combine the two.*

## NOTCH1/2/3 — NRR + TM cassette ✅ strong
- **Gordon 2009, PLoS One 4:e6613** (PMID 19701457) — minimal NOTCH1 & NOTCH2 = **NRR (LNR+HD) + TM + ICD**, surface-export measured directly.
- **Aste-Amézaga 2010, PLoS One 5:e9094** (PMID 20161710) — membrane-anchored minimal NOTCH1 (EGF-LBD or NRR) displayed for antibody FACS.
- **Li 2008, J Biol Chem 283:8046** (PMID 18182388) — truncated **NOTCH3** (leader + NRR + TM) surface-expressed in 293T.
- **Morsut/Lim 2016, Cell 164:780** (PMID 26830878) — **synNotch**: NOTCH1 NRR+TM core (~360 aa) is the standard compact surface-displayed chassis. NRR+TM ≈ 0.9–1.2 kb. Keep NRR intact (uncleaved ΔE/NEXT forms shed).

## ROS1 — ectodomain-deletion / ectodomain-swap ✅ (works; note it can activate kinase)
- **Jones et al. 2024, Sci Rep 14:22191** (PMID 39333184) — ECD-deleted ROS1 (drop ~1,800-aa ectodomain, keep membrane-proximal + TM + kinase) is membrane-expressed + active.
- **Xiong 1996, Mol Cell Biol 16:1509** (PMID 8657124) — EGFR-ectodomain↔c-Ros chimeras: swap the giant ectodomain for a small one → surface-expressed, ligand-activatable.
- Caveat: ectodomain truncation tends to make the kinase constitutively active → for a display antigen, use kinase-dead or drop the kinase (membrane-proximal ectodomain + TM only).

## ITGB4 (β4) — cytoplasmic-tail truncation ✅ strong
- **Spinardi 1995, J Cell Biol 129:473** (PMID 7721947) — **tail-less β4** (whole ~1,000-aa cytoplasmic domain deleted) reaches the surface, pairs with α6, binds laminins. Collapses ~5.5 kb → ~2.3 kb.
- **Niessen 1997, J Cell Sci 110:1705** (PMID 9264458) — residue-level map: truncate after ~aa 1355 and keep surface α6β4. Kim 1997 (PMID 8988037) + Murgia 1998 (PMID 9670011) confirm distal-tail-truncated β4 is the stably-tolerated, surface form. Natural short-tail isoform ITGB4E (PMC8895044).

## LY75 / DEC-205 — endocytic tail + TM (mini-ectodomain empirical) ⚠️ partial
- **Mahnke 2000, J Cell Biol 151:673** (PMID 11062267) — DEC-205 **TM + cytoplasmic tail** on a CD16 ectodomain: surface-expressed + endocytic; uptake driven by a coated-pit motif + distal **EDE triad**.
- **Shrimpton 2009** (PMID 19135256) — CTLD3+4 and CTLD9+10 carry ligand function (soluble Fc fusions). **Gully 2020** (PMID 33257321) cryo-EM: the 10 CTLDs form interlocked rings → don't excise internal CTLDs; keep contiguous sub-modules. *No published mini-DEC205 yet — design = intact CTLD block + TM + endocytic tail.*

## CSPG4 / NG2 — membrane-proximal D3 + TM ✅
- **Bluemel 2010, Cancer Immunol Immunother 59:1197** (PMID 20309546) — surface-displayed MCSP deletion mutants; **membrane-proximal D3 domain carries the most potent T-cell-engager epitope**, and smaller surface antigens lyse better → D3(+juxtamembrane)–TM cassette.
- **Nishiyama 1995, Mol Biol Cell 6:1819** (PMID 8590808) + **Stallcup 2001, J Cell Sci 114:2315** (PMID 11493670) — cytoplasmic-tail-deleted / ectodomain-truncated NG2 still surface-displays.

## STAB2 (stabilin-2) — the 190-kDa HARE ✅ strong (STAB1 gap)
- **Harris 2004, J Biol Chem 279:36201** (PMID 15208308) — recombinant **190-kDa HARE = C-terminal half of Stab2 (res ~1136–2551)** expressed alone: surface-expressed (~30–50%), fully endocytic for hyaluronan.
- **Harris 2006** (PMID 17145755) + **Pandey 2016, PLoS One 11:e0154124** (PMID 27100626) — maps the 4 cytoplasmic endocytic motifs (NPLY2519 key). Still ~4.25 kb → trim further. **STAB1: no published minimal construct** — graft the DDSLL tail motif onto the STAB2-style design.

## PCDH15 / CDH23 / USH2A — Usher gene-therapy mini-genes ✅ (PCDH15/USH2A) / ✗ (CDH23)
- **PCDH15:** Ivanchenko 2023, Nat Commun 14:2400 (PMID 37100771) — **mini-PCDH15** (delete 3–5 of 11 EC repeats) forms tip links + rescues hearing; dual-AAV full-length also works (PMID 38014037).
- **USH2A:** Pendse 2020 (bioRxiv 2020.02.04.934240) — **exon-13-deleted usherin** rescues cochlear+retinal phenotypes; ASO exon-skip QR-421a/ultevursen (Dulla 2021, PMID 33895329) in clinic.
- **CDH23:** triple-AAV of full-length **failed** to restore function (Yoshimura 2023, PMID 38132436); no working mini-CDH23 yet → least de-risked.
- Lesson: these large adhesion proteins tolerate **internal EC-repeat deletion** while keeping partner-binding/surface function.

## Mucins MUC1 / MUC4 / MUC17 (template for MUC16) — SEA + TM, drop the VNTR ✅ strong
- **MUC17:** Schneider 2019, Biochem J 476:2281 (PMID 31387973) — **MUC17-3TR** (3 of ~60 repeats + SEA + TM + tail): SEA autocleaves, surface-presented.
- **MUC4:** Moniaux 2007, Br J Cancer 97:345 (PMID 17595659) — **mini-MUC4** (~10% VNTR, keep NIDO/AMOP + vWD/EGF/TM/CT): membrane localization like WT.
- **MUC1:** Engelmann 2005, Glycobiology 15:1111 (PMID 15972891) — MUC1-M2 (truncated repeats + **SEA** + TM) recycles from the surface; natural MUC1/Y isoform lacks the VNTR entirely and is surface-expressed.
- **Caveats:** SEA autoproteolysis needs enough flanking context (don't over-trim near SEA); surface delivery can be N-glycosylation-dependent near the SEA/TM junction.

---

## PTPRZ1 / PTPRS / PTPRF (receptor phosphatases) — delete the Ig+FN3 array ✅
- **PTPRZ1:** Lorente 2005, Neuro-Oncology 7:154 (PMID 15831233) — the short RPTPβ isoform (**860 aa deleted from the ECD**, keep CAH+FN3+TM+2 PTP) surface-expresses in U87 GBM; Maeda 1998 (PMID 9660874) — short PTPζ-B at the neuronal surface.
- **PTPRS:** Han 2018, J Neurosci 38:6700 (PMID 29934346) — **ΔIg / ΔFN / ΔEcto** constructs all traffic to the HEK293T surface.
- **PTPRF/LAR:** Serra-Pagès 1994, J Biol Chem 269:23632 (PMID 8089133) — the P-subunit (short ECD stub + TM + 2 PTP) is surface-expressed; caveat: ectodomain **sheds/internalizes** (Aicher 1997, PMID 9245795).
- Strategy: the Ig+FN3 ectodomain array is heavily deletable; keep both PTP domains; ectodomain-swap chimera (CD45 precedent) is the fallback.

## SIGLEC1 / CD169 — V-set domain 1 + a few Ig + native TM/tail ⚠️ partial
- **Nath 1995, J Biol Chem 270:26184** (PMID 7592823) — **N-terminal V-set domain 1 alone** is necessary + sufficient for sialic-acid binding (soluble fusion).
- **Klaas 2023, Glycoconj J 40:213** (PMC10027830) — surface ectodomain-truncations: **≤4 Ig domains lose cell-ligand binding** (glycocalyx burial) → keep enough stalk (≈d1–d6).
- **Delputte 2011, PLoS One 6:e16827** (PMID 21347240) — full-length Sn is endocytic (delivers toxin/antigen into macrophages); endodomain deletion cut expression (Van Breedam 2013, PMID 23740482) → keep native TM + cytoplasmic tail.
- *No single short + surface + endocytosis-verified construct — design = d1–d~6 + native TM + tail.*

## MUC16 (CA125) — membrane-proximal C-terminal "cN" fragment (MUC-CD / MUC16ecto) ✅ strong, clinically de-risked
Full length ~14,507 aa; the field builds on the **retained C-terminal fragment** left on the cell after CA125 sheds. Residues counted from the C-terminus ("cN" = last N aa = ectodomain stub + TM + ~32-aa tail).
- **Rao 2015, PLoS One** (PMID 25965947) — the definitive truncation ladder: **MUC16c114** and **c344**-GFP stably surface-expressed (FACS). **c344 displays the CA125/OC125 epitope (and sheds it); c114 does not carry OC125 but is surface+ via the membrane-proximal 4H11 antibody.** 114 C-terminal aa suffice for surface trafficking + oncogenic signaling.
- **Thériault/Boivin 2011, Gynecol Oncol** (PMID 21421261) — N-terminal-deleted MUC16CTD reaches the surface (ER-trapped scFv blocks its delivery, confirming normal trafficking); **cytoplasmic-tail-deleted MUC16TMU still expresses** — the ~32-aa tail is dispensable for display, required only for intracellular signaling. (Das 2015, PMID 25691062 — 114 aa converges as the minimal functional unit.)
- **Cleavage caveat:** MUC16 is cut at a ~12-aa juxtamembrane motif (AWFPLDSNGTLP) in the acidifying Golgi (Das 2015, PMID 26044153) — include that motif + distal sequence and the distal part (incl. the CA125 epitope in c344) sheds. For a durably displayed surface antigen, stay membrane-proximal to the site or engineer it out.
- **Target validation (retained fragment is a real, internalizing surface antigen):** 4H11/9C9/4A5 mAbs bind + internalize (Rao 2010, PMID 20453816); **MUC-CD CAR-T** clears xenografts (Chekmasova 2010, PMID 20628030); **4H11-28z + IL-12 CAR** (Koneru 2015, PMID 25949921) advanced to a **Phase I trial of MUC-16ecto CAR-T** (Koneru 2015, PMID 25890361). 4H11 epitope crystallized to a **26-aa membrane-proximal segment** (Lee 2024, PMID 38374055).
- **Design:** a **cleavage-resistant c114-style construct** (retained stub + TM ± tail; ~114 aa) packages trivially and matches the most clinically validated MUC16 surface target; use the longer **c344** (one full SEA repeat) only if you must display the canonical CA125/OC125 epitope itself.

## Plexins (PLXNB1 headline; A1/A4/B2/B3/D1 family) — ΔEcto (TM + cytoplasmic) or small-ectodomain chimera ✅ strong
Ectodomain (~1,500 aa: Sema β-propeller + PSI + 3× IPT) is the entire size problem; TM + cytoplasmic GAP region is only ~600 aa (~1.9 kb). The ectodomain is an autoinhibitory clamp, so deleting it is well precedented.
- **PLXNB1 ΔEcto** (TM + full cytoplasmic C1–RBD–C2 GAP only) — **Oinuma/Negishi 2004, J Neurosci** (PMID 15601954): reaches the membrane, self-clusters ligand-independently, is a **constitutively active R-Ras GAP** (Rnd1-dependent). The exact PLXNB1 truncation — direct headline precedent. Confirmed for physiology (β1-integrin/migration) in Oinuma 2006 (PMID 16702230).
- **Founding "sema-minus":** Sema-deleted PlexinA1 (ΔSema) is surface-localized + constitutively active (Takahashi & Strittmatter 2001, PMID 11239433). **Modern, multi-family:** PLXNA1ΔECD + PLXNB2ΔECD both surface-express + signal (Pinskey/Giger 2022, eLife, PMID 36169302) — the trick transfers across A- and B-class.
- **Ligand-inducible chimera:** TrkA ecto+TM fused to Plexin-B2 cytoplasmic domain surface-expresses and signals on NGF (Perrot/Gutkind 2002, PMID 12183458) — swap the bulky ectodomain for a small dimerizable module for clean on-demand control.
- **Keep the C-terminal PDZ motif** — its deletion abolishes Rho output (Oinuma 2003, PMID 12730235; chimeric-B1 angiogenesis, Basile 2004, PMID 15289326).
- **Ligand-binding-only module:** PLXNB1 Sema+PSI (domains 1–2, ~first ~560 aa) folds and binds Sema4D as an isolated fragment (Janssen 2010, PMID 20877282) — anchor via TM if you want a binder/antigen instead of a signaler.
- **Caveats:** ΔEcto removes the autoinhibitory brace → **constitutively active**, not ligand-gated (structure: Kong 2016, PMID 27397516). Class-A plexins need **neuropilin-1** to bind class-3 semaphorins; **PLXNB1 binds membrane Sema4D directly (no neuropilin)** — an advantage for using B1 as the representative.
