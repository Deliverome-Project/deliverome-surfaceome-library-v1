# L1 design review: what is still open before the order

The [design brief](BRIEF.md) is the design. This note tracks only what remains open before
sending the order to synthesis, what was completed in the `2026-11-24-pre-freeze-c` rebuild, and the
known risks we are carrying knowingly.

Status is as of **2026-11-24**, against the order of record
(`l1_order_inserts_recoded.fasta`, cut **`2026-11-24-pre-freeze-c`**, SHA-256 `ef810d1d…`, built
from commit `7dc06001`): **2,417** genes, **2,762** active CDSs, **2,743** built CDSs (**19**
held), **3,375** ordered constructs across **36** 96-well plates, **3,375 / 3,375 (100%)**
ColabFold-refolded (`3,173` pass, `202` flagged with recorded dispositions, `0` pending) and
**3,375 / 3,375** passing plasmid QC. Every count is measured by
[`brief_numbers.py`](brief_numbers.py).

## Completed in the `2026-11-24-pre-freeze-c` rebuild

Every decision locked in on 2026-11-23 and every pre-order biological and sequence-motif audit
finding from 2026-11-24 is now **built, refolded, verified, and cut** into `2026-11-24-pre-freeze-c`
(no Modal compute or membership items remain queued):

| Area | Resolution built into `2026-11-24-pre-freeze-c` |
|---|---|
| **Domestication & synthesis motif hardening (`18` inserts cleaned)** | Eliminated all domestication-created `ATTAAA` poly(A) hexamers (`4` inserts: `KLRB1`, `SLC46A1` ×3), all 6-mer `GT[AG]AGT` consensus splice donors (`10` inserts: `EPHA1`, `FGFR1` ×2, `NRG1` ×3, `OR4F6`, `OR56B1`, `PLPPR1` ×2), and `4` fixable early (`≤ 60 nt`) ATGs (`EPHA10`, `HTR1B`, `KCNE3` ×2) via common synonymous human codons (`≥ 24%` usage); hardened `seqtools.domesticate()` and `synthesis_compliance.py` in CI (`0` poly(A) or 6-mer splice donors introduced across all 3,375 inserts). |
| **Pre-order biological audit (`14` internal CDSs fixed, `0` CDSs lost)** | Exhaustive audit across all built loop arms dropped 9 defective/redundant arms (`ANO4 @685` cytoplasmic `z = −37.9 Å`, 4 `SLC44A1` membrane-core arms, `LYPD6 @146` ω-proximal, 3 `<25 aa` close pairs on `HCN2, AQP4, SLC7A11`) and swapped 8 one-for-one to ColabFold-passing sites (`ATRN @1038`, `CDH3 @645`, `CDH5 @586` out of cleaved propeptides; `MDGA2 @151`, Modal-verified `BST2 @149` & `PRND @97` ≥12 aa upstream of GPI ω; `TRPV3 @615` & `XPR1 @503` onto distinct loops). All invariants enforced in CI. |
| **Gene list (`2,417` genes)** | Adopted `2026-11-23-abundance-v8` surface-proteomics tranche of 16 L1-fitting genes (`+ANO9, KLRG2, RNF150, TMEM117`; `-CBARP, FAM234B, LNPEP, PRRT3, SLC22A23, TMEM94`); dropped all 8 MHC class II chains (`HLA-DPB1, DQA1, DQB1, DRA, DRB1/3/4/5`), eliminating the 208-nt `HLA-DRB1/DRB3` swap pair. |
| **Segregating human LoF alleles** | Built high-confidence/0-diff `OR4Q2` (patch), `FCGR2C`, `OR10AC1`, `OR8J2`, and `SIGLEC16`; dropped lower-confidence `OR52P1` (5 unphased residue diffs). |
| **14 GTF-MANE vs. APPRIS genes** | Switched 5 to MANE (`APP, GRIN1, KCNMA1, NFASC, SIRPA`); kept 9 on APPRIS on PeptideAtlas support (`ADGRG6, ATP1A3, CSPG5, DPP10, DPP6, KCNQ2, ROBO2, SLC12A3, SLC4A7`). |
| **Isoform panel redraw (`150` forms)** | Relaxed reference N-term-side & TM1 matching in [`gate_isoform_cohort.py`](../gate_isoform_cohort.py) with short-fragment guard (`1,749/1,899` pass); redrawn across all 6 strata with fresh DeepTMHMM v1 & ColabFold structures. |
| **Literature ligand-contact screen (`556/556` internal genes)** | Finished the final 52 internal genes (100% coverage). Caught 4 hard vetoes (`SLC1A3 @156`, `STEAP4 @234`, `XKR8 @40` dropped; `SLC6A2 @197 → @195` swapped to 0-clash PASS) and flagged 7 soft/partner hits. |
| **25 steric-clash & propeptide upgrades** | Applied 21 zero-GPU internal upgrades (incl. `BACE1 @44 → @56`, `SCNN1A @168`, `PTCH1 @301`, `SLC4A7 @777`, `FCAMR @296`, `OTOP2 @215`, `SLCO2A1 @437`, `SLC46A2 @436`) + 4 Modal ColabFold swaps (`AREG @54 → @193` post-ADAM17 stalk, `FZD1 @243 → @588`, `OTOP3 @251 → @244`, `SLC6A20 @532 → @536`); clash validation reproduced offline by [`refold_check_clash_validation.py`](../refold_check_clash_validation.py). |
| **Topology & GPI resolutions** | Primary literature settled all 4 GPI conflicts (`CFC1, CRIPTO, RAET1G` dual TM+GPI, `GML`) and the 4 UniProt-discordant C-term genes (`LYSMD3, LYSMD4, SLC24A1` kept C-term; `TMEM222` reference moved to multipass residue 1). |
| **N-after-SP borderline rule retired** | Retired N-term-only `lddt_borderline_counted_as_pass` (`LDDT_MARGIN = 0.0`): `SNORC, SMO, SLITRK2` shifted to PASS sites, `CD300LB` rescued by Igκ PASS, `CIST1` & `FOLR1` built and flagged uniformly; bonus rescue `IL12RB1` cyno PASS. |
| **High & Not Accepted–only synthesis recoding** | Native codons preserved on Low/Medium inserts; `223` inserts (`798` codons across `181` CDSs: `780` High/Not-Accepted synthesis codons + `18` domestication motif-cleanup codons) are minimally recoded, leaving **`0` Not Accepted** and **`24` High** (protected native SP codons). |
| **Structure models & metadata** | Folded 24 missing natives + 63 tagged constructs on Modal (`2,762/2,762` active CDSs have native models; `3,375/3,375` constructs folded); added `design_id` to [`l1_rbp_motif_cassettes.tsv`](../l1_rbp_motif_cassettes.tsv); accepted `OR9G1/OR9G9` (102 nt) and `FCGR2A/FCGR2C` (174 nt) in CI. |
| **ISS barcode prefix standard & clean export** | Adopted **12-nt `d ≥ 3`** for future reservations in [`barcode_prefix_screen.py`](../barcode_prefix_screen.py) (+8,896 headroom; L1's 3,375 constructs satisfy 10-nt `d ≥ 3`). Clean-path export (`library/l1/*`) pushed to [Deliverome-Project/deliverome-surfaceome-library-v1](https://github.com/Deliverome-Project/deliverome-surfaceome-library-v1). |

## Open before order freeze

Only two scientific design questions and a small set of operational pre-order steps remain:

| Item | Type | What closes it |
|---|---|---|
| **Promoter: EFS (current) vs. EF1a** ([brief Q1](BRIEF.md#decisions-we-would-like-reviewers-to-weigh-in-on)) | Scientific decision | Matched-MOI test of a few dim-to-moderate library ORFs in direct `IRES2–PuroR` (EFS vs. EF1a), or reviewer consensus weighing EF1a's 1.5–2.7× display gain against +967 bp titer cost and intron A splicing risk. |
| **T7 0-nt spacer before `JT832`** ([brief Q2](BRIEF.md#decisions-we-would-like-reviewers-to-weigh-in-on)) | Scientific decision | Single-copy pilot or reviewer feedback on whether `STOP · CC · T7 · JT832` (3 nt `GGG` before `JT832`, 7 nt `GGGAGTC` before the 15-nt padlock arm) needs a short upstream spacer before freezing synthesis. |
| **Vendor re-quote on `2026-11-24-pre-freeze-c`** | External send | Submit `l1_order_inserts_recoded.fasta` (`3,375` inserts; `0` Not Accepted, `24` High in protected SP codons) across candidate synthesis vendors to confirm manufacturability, prep tier (`8–10 µg` target vs. `4 µg` floor), and turnaround. |
| **Concrete barcode assignment** | Freeze step | Assign 3,375 barcodes from the 3,644-barcode provisional block ([`l1_barcode_block_provisional.tsv`](../l1_barcode_block_provisional.tsv) + [`l1_barcode_block_extension_2026-11-15.tsv`](../l1_barcode_block_extension_2026-11-15.tsv)) to replace the `20-N` placeholder and release the 269 unused rows. |
| **Plate format & 6-plasmid titer ladder** ([brief section 10](BRIEF.md#10-how-the-library-is-ordered)) | Wet-lab / ordering | Choose low-binding V-bottom vs. skirted PCR plates; run the 6-plasmid (`p0028, p0043, p0027, p0037, p0024, p0029`, excluding `HEG1`) low-MOI qPCR/ddPCR length-titer ladder to set `FACTOR` in [`order_plan.py`](order_plan.py) before subplate pooling (`2`–`4` subpools/plate directly from `100 ng/µL`). |
| **3 native weak-SP cleavage calls (`TMEM81, MXRA8, IL1RAPL2`)** | Minor protein review | DeepTMHMM v1 and SignalP 6 disagree by > 5 aa; currently kept native and flagged. |
| **7 internal constructs splitting a predicted N-glycan sequon** | Minor protein review | Flagged in internal review; kept unless reviewer feedback prefers a shift. |
| **Licensing check on private S3 storage for DTU binaries** | Operations | Confirm whether licensed DeepTMHMM 1.0.24 / SignalP 6.0+h packages may be archived in private S3 alongside their recorded hashes. |

## Known risks carried knowingly

| Risk | Mitigation or empirical check |
|---|---|
| **Local synthesis rules approximate Twist's proprietary scorer** | Vendor re-quote of `l1_order_inserts_recoded.fasta` before ordering. |
| **216 chaperone-dependent CDSs + 31 K562-negative co-receptor genes** | Retained in the main pool and flagged outside the K562 display denominator until screened in lineage-matched or factor-complemented hosts. |
| **No within-gene Igκ vs. leaderless control on 567 multipass CDSs** | Evaluated against panel controls (`CCR5, OR51E2, CNR1`) and K562 display readout. |
| **9 of 21 multipass taste receptors & 14 class A GPCRs fail refold** (buried SpyTag003 or packed epitope; 10-aa spacer rescued 0/12) | Kept as built and flagged ([`multipass_refold_dispositions.tsv`](../multipass_refold_dispositions.tsv)); ALFA stays exposed (`RSA 0.56–0.74`) as display readout. |
| **Tag-site model v5 had 4 misnumbered training proteins (858 of 35,650 rows)** | All 11 affected L1 CDSs were redesigned on corrected features; full ranker retrain deferred to a future library release. |
| **Two native-ORF paralog pairs share > 100 nt before the barcode** (`OR9G1/OR9G9` 102 nt, `FCGR2A/FCGR2C` 174 nt) | Accepted as known covariates (shared sequence is native 3′ ORF, never recoded); actual recombination measured by long-read ORF–barcode sequencing. |
| **Non-assembly Type IIS enzymes** (`AarI/PaqCI` in 1,507, `BtgZI` in 1,419, `SapI` in 3,375 inserts) | Working as intended — cloning uses only `BbsI`, `Esp3I`, and `BsaI` (`0` internal sites across all 3,375 inserts). |
| **90 AFDB payload URLs unarchived** (`tag_refold_afdb_unarchived.tsv`, incl. 27 HTTP 404 in AFDB v6) | Every active L1 CDS (`2,762/2,762`) has a verified archived AFDB or ColabFold native model; every ordered construct (`3,375/3,375`) has its own ColabFold tagged fold. |
