# L1 design review: what is still open before the order

The [design brief](BRIEF.md) is the design. This note tracks only what remains open before
sending the order to synthesis, what was completed in the `2026-11-25-pre-freeze-e` rebuild, and the
known risks we are carrying knowingly.

Status is as of **2026-11-25**, against the order of record
(`l1_order_inserts_recoded.fasta`, cut **`2026-11-25-pre-freeze-e`**, SHA-256 `658682f8…`, built
from commit `4422a9c4`): **2,399** genes, **2,740** active CDSs, **2,730** built CDSs (**10**
held), **3,352** ordered constructs across **35** 96-well plates, **3,352 / 3,352 (100%)**
ColabFold-refolded (`3,153` pass, `199` flagged with recorded dispositions, `0` pending) and
**3,352 / 3,352** passing plasmid QC. Every count is measured by
[`brief_numbers.py`](brief_numbers.py).

## Completed in the `2026-11-25-pre-freeze-e` rebuild

Every membership tightening, new addition, biological audit fix, and sequence-motif hardening
finding is now **built, refolded, verified, and cut** into `2026-11-25-pre-freeze-e`
(no Modal compute or membership items remain queued):

| Area | Resolution built into `2026-11-25-pre-freeze-e` |
|---|---|
| **Uniform ≥ 30-aa ectodomain gate (`−20` genes / `−25` CDSs, `−9` held)** | Fixed `final_library.py` tier `LL` to enforce `ecd_ok(r)` uniformly across all non-GPI tiers and added a CI biological-invariant guardrail (`0` active non-GPI reference genes below 30 aa). Retired all 20 `<30 aa` `LL`-tier genes (`CATSPERH`, `CATSPERQ`, `CXorf66`, `MS4A5`, `RELL2`, `SMIM7`, `TEX38`, `TMEM134`, `TMEM139`, `TMEM169`, `TMEM174`, `TMEM210`, `TMEM222`, `TMEM238`, `TMEM266`, `TMEM268`, `TMEM273`, `TMEM86A`, `TMEM86B`, `TSPO2`) plus their 5 panel derivatives (`TMEM139` mouse, `TMEM210` cyno & mouse, `TMEM222` isoform, `TMEM86B` mouse), dropping held CDSs from `19 → 10` and shrinking the order from `36 → 35` plates. |
| **Surface-proteomics `likely`-tier additions (`+LNPEP`, `+PCDH11X`)** | Stopped excluding buildable surface-proteomics gap hits for `likely` tier / `mostly_intracellular` specificity (enforcing only the 5,000-nt ORF cap, which still excludes 7,476-nt `IGF2R`). Added `LNPEP` (`Q9UIQ6`, 3,078 nt, C-terminal, PASS refold `lDDT 0.966`) and `PCDH11X` (`Q9BZA7`, 4,044 nt, N after SP after Q24, SignalP 6 prob `0.9786`, Modal ColabFold PASS `lDDT 0.956`, `0` epitope contacts). |
| **`OR52P1` built from UniProt `Q8NH57`** | Restored `ENST00000438474` in [`allele_restorations.tsv`](../allele_restorations.tsv) (2 frameshift fixes + 5 coding substitutions) to encode UniProt `Q8NH57` (321 aa) **100% identically (`0` diffs)** while keeping ~98.5% native human codons; multipass Igκ Modal ColabFold refold PASSED (`lDDT 0.993`, `TM 0.997`, `RMSD 0.36 Å`, `0` epitope contacts). |
| **Domestication & synthesis motif hardening (`18` inserts cleaned)** | Eliminated all domestication-created `ATTAAA` poly(A) hexamers (`4` inserts: `KLRB1`, `SLC46A1` ×3), all 6-mer `GT[AG]AGT` consensus splice donors (`10` inserts: `EPHA1`, `FGFR1` ×2, `NRG1` ×3, `OR4F6`, `OR56B1`, `PLPPR1` ×2), and `4` fixable early (`≤ 60 nt`) ATGs (`EPHA10`, `HTR1B`, `KCNE3` ×2) via common synonymous human codons (`≥ 24%` usage); hardened `seqtools.domesticate()` and `synthesis_compliance.py` in CI (`0` poly(A) or 6-mer splice donors introduced across all 3,352 inserts). |
| **Pre-order biological audit (`14` internal CDSs fixed, `0` CDSs lost)** | Exhaustive audit across all built loop arms dropped 9 defective/redundant arms (`ANO4 @685` cytoplasmic `z = −37.9 Å`, 4 `SLC44A1` membrane-core arms, `LYPD6 @146` ω-proximal, 3 `<25 aa` close pairs on `HCN2, AQP4, SLC7A11`) and swapped 8 one-for-one to ColabFold-passing sites (`ATRN @1038`, `CDH3 @645`, `CDH5 @586` out of cleaved propeptides; `MDGA2 @151`, Modal-verified `BST2 @149` & `PRND @97` ≥12 aa upstream of GPI ω; `TRPV3 @615` & `XPR1 @503` onto distinct loops). All invariants enforced in CI. |
| **High & Not Accepted–only synthesis recoding** | Native codons preserved on Low/Medium inserts; `223` inserts (`798` codons across `181` CDSs: `780` High/Not-Accepted synthesis codons + `18` domestication motif-cleanup codons) are minimally recoded, leaving **`0` Not Accepted** and **`24` High** (protected native SP codons). |
| **Section 7 overhaul & topology schematic** | Rewrote Section 7 cleanly to eliminate GPI/override/Igκ repetition, added a five-panel membrane topology & placement schematic (`figures/placement_rules_schematic.png`), and verified 99.5% experimental topology concordance (`1,114/1,120`, `CATSPERH` gone). |
| **ISS barcode prefix standard & clean export** | Adopted **12-nt `d ≥ 3`** for future reservations in [`barcode_prefix_screen.py`](../barcode_prefix_screen.py) (+8,896 headroom; L1's 3,352 constructs satisfy 10-nt `d ≥ 3`). Clean-path export (`library/l1/*`) pushed to [Deliverome-Project/deliverome-surfaceome-library-v1](https://github.com/Deliverome-Project/deliverome-surfaceome-library-v1). |

## Open before order freeze

Only two scientific design questions and a small set of operational pre-order steps remain:

| Item | Type | What closes it |
|---|---|---|
| **Promoter: EFS (current) vs. EF1a** ([brief Q1](BRIEF.md#decisions-we-would-like-reviewers-to-weigh-in-on)) | Scientific decision | Matched-MOI test of a few dim-to-moderate library ORFs in direct `IRES2–PuroR` (EFS vs. EF1a), or reviewer consensus weighing EF1a's 1.5–2.7× display gain against +967 bp titer cost and intron A splicing risk. |
| **T7 0-nt spacer before `JT832`** ([brief Q2](BRIEF.md#decisions-we-would-like-reviewers-to-weigh-in-on)) | Scientific decision | Single-copy pilot or reviewer feedback on whether `STOP · CC · T7 · JT832` (3 nt `GGG` before `JT832`, 7 nt `GGGAGTC` before the 15-nt padlock arm) needs a short upstream spacer before freezing synthesis. |
| **Vendor re-quote on `2026-11-25-pre-freeze-e`** | External send | Submit `l1_order_inserts_recoded.fasta` (`3,352` inserts; `0` Not Accepted, `24` High in protected SP codons) across candidate synthesis vendors to confirm manufacturability, prep tier (`8–10 µg` target vs. `4 µg` floor), and turnaround. |
| **Concrete barcode assignment** | Freeze step | Assign 3,352 barcodes from the 3,644-barcode provisional block ([`l1_barcode_block_provisional.tsv`](../l1_barcode_block_provisional.tsv) + [`l1_barcode_block_extension_2026-11-15.tsv`](../l1_barcode_block_extension_2026-11-15.tsv)) to replace the `20-N` placeholder and release the 292 unused rows. |
| **Plate format & 6-plasmid titer ladder** ([brief section 10](BRIEF.md#10-how-the-library-is-ordered)) | Wet-lab / ordering | Choose low-binding V-bottom vs. skirted PCR plates; run the 6-plasmid (`p0028, p0043, p0027, p0037, p0024, p0029`, excluding `HEG1`) low-MOI qPCR/ddPCR length-titer ladder to set `FACTOR` in [`order_plan.py`](order_plan.py) before subplate pooling (`2`–`4` subpools/plate directly from `100 ng/µL`). |
| **3 native weak-SP cleavage calls (`TMEM81, MXRA8, IL1RAPL2`)** | Minor protein review | DeepTMHMM v1 and SignalP 6 disagree by > 5 aa; currently kept native and flagged. |
| **7 internal constructs splitting a predicted N-glycan sequon** | Minor protein review | Flagged in internal review; kept unless reviewer feedback prefers a shift. |
| **Licensing check on private S3 storage for DTU binaries** | Operations | Confirm whether licensed DeepTMHMM 1.0.24 / SignalP 6.0+h packages may be archived in private S3 alongside their recorded hashes. |

## Known risks carried knowingly

| Risk | Mitigation or empirical check |
|---|---|
| **Local synthesis rules approximate Twist's proprietary scorer** | Vendor re-quote of `l1_order_inserts_recoded.fasta` before ordering. |
| **217 chaperone-dependent CDSs + 31 K562-negative co-receptor genes** | Retained in the main pool and flagged outside the K562 display denominator until screened in lineage-matched or factor-complemented hosts. |
| **No within-gene Igκ vs. leaderless control on 565 multipass CDSs** | Evaluated against panel controls (`CCR5, OR51E2, CNR1`) and K562 display readout. |
| **9 of 21 multipass taste receptors & 14 class A GPCRs fail refold** (buried SpyTag003 or packed epitope; 10-aa spacer rescued 0/12) | Kept as built and flagged ([`multipass_refold_dispositions.tsv`](../multipass_refold_dispositions.tsv)); ALFA stays exposed (`RSA 0.56–0.74`) as display readout. |
| **Tag-site model v5 had 4 misnumbered training proteins (858 of 35,650 rows)** | All 11 affected L1 CDSs were redesigned on corrected features; full ranker retrain deferred to a future library release. |
| **Two native-ORF paralog pairs share > 100 nt before the barcode** (`OR9G1/OR9G9` 102 nt, `FCGR2A/FCGR2C` 174 nt) | Accepted as known covariates (shared sequence is native 3′ ORF, never recoded); actual recombination measured by long-read ORF–barcode sequencing. |
| **Non-assembly Type IIS enzymes** (`AarI/PaqCI` in 1,495, `BtgZI` in 1,413, `SapI` in 3,352 inserts) | Working as intended — cloning uses only `BbsI`, `Esp3I`, and `BsaI` (`0` internal sites across all 3,352 inserts). |
| **Unarchived historical AFDB payload URLs** (`tag_refold_afdb_unarchived.tsv`, incl. 27 HTTP 404 in AFDB v6) | Every active L1 CDS (`2,740/2,740`) has a verified archived AFDB or ColabFold native model; every ordered construct (`3,352/3,352`) has its own ColabFold tagged fold. |
