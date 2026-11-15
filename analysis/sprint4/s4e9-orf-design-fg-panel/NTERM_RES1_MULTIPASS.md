# N-terminal tag at residue 1 on multipass members: leader or no leader

**Scope.** The 563 L1 members that [`l1_tag_paths.tsv`](l1_tag_paths.tsv) places on the
rung "N at residue 1, multipass": the member's own DeepTMHMM v1 topology starts
extracellular, carries no signal peptide, and has more than one TM helix. For these the tag
cassette goes directly after Met1, so the cassette and the native N-terminal tail (the
residues before TM1) are translocated by TM1 itself, acting as a signal anchor. `tag_cassette`
is `deferred` for every one, so no DNA is bound to any choice made here.

**Leader plan.** The 495 GPCRs take the Igκ leader by default (Becca, 2026-11-10, on Aashish
Manglik's suggestion). For the 68 non-GPCRs, a literature review recommends Igκ for 63, moving
2 off this rung and holding 3 (see "Non-GPCRs: the literature review"); those verdicts await
Becca's confirmation. Each member's plan is in the `planned_leader` column of the per-member
table.

**Two snapshots.** The analysis sections ("Who is on the rung" through "Non-GPCRs: the
literature review") were written on 2026-11-10 against `main`: 563 members, with the
ortholog panel from before Compara 116. The per-member tables and the built arms are
regenerated from PR #283 at commit 13f3c1b8. That state has the Compara 116
ortholog panel, the gated isoform redraw, the reference fixes and the lane overrides, and it
puts 561 members on the rung. [`multipass_rung_changes.tsv`](multipass_rung_changes.tsv)
lists every member that left or entered, matched on `design_id`, with its reason and decided_by:
- 34 retired: 26 orthologs under the Compara 116 CDS policy, and 8 isoforms. Seven were removed
  by the gate or the redraw; the OPN3 isoform was dropped by the proteoform audit.
- 5 moved: ADGRF3 and GLP2R to "N after SP"; MFSD2A ×2 and AQP12B to internal (`lane_overrides.tsv`).
- 37 entered: 34 Compara 116 orthologs, and 3 isoforms (ABCC3, GPR161, HTR7).
- 2 routed in from "N after SP" on 2026-11-14: OR11H6 (`l1-p1381`) and TMEM63B (`l1-p2199`),
  which have no signal peptide (see "Decisions taken under Becca's delegation").
- 1 routed in from "C-terminal" on 2026-11-14: SLC41A2 (`l1-p1957`), whose C-terminus is
  cytoplasmic by experiment.
- 1 routed in from "N after SP" on 2026-11-15: human OR10G6 (`l1-p1362`), whose N-terminus
  Tedman displayed uncleaved. The rung now has 565 members.

**Terms used throughout.**

| Term | Meaning here |
|---|---|
| Cassette | The 51-aa L1 tag: GGGGS, SpyTag003, GGGGSGGGGS, P-ALFA-P, GGGGS. It carries 8 Lys/Arg and 5 Asp/Glu, net +3. It is used in full, leading GGGGS included: without that spacer, DeepTMHMM v1 called a cleavable signal peptide ending inside SpyTag003 on 12 of 62 residue-1 tag-alone arms, against 0 of 62 with it (the type III rung's spacer test). |
| Igκ leader | The murine Igκ signal peptide (`METDTLLLWVLLLWVPGSTG`) plus the +1 Asp of the mature κ chain, as built in [`build_gpcr_rescue_arms.py`](build_gpcr_rescue_arms.py). It is cleaved co-translationally, so the mature protein starts with Asp. |
| N-tail | Residues before the first v1 helix in the member's ordered protein. |
| Tagged tail | N-tail + 51: what TM1 has to translocate when there is no leader. |
| Charge difference (Δ) | Lys+Arg in the 15 residues after TM1 minus those in the 15 before it. The more positive flank faces the cytosol ([Hartmann et al., 1989](https://pubmed.ncbi.nlm.nih.gov/2762295/)), so Δ > 0 agrees with the N-out orientation this rung assumes. |
| Tedman screen | [Tedman et al. (2026)](https://pubmed.ncbi.nlm.nih.gov/42736269/): 766 canonical human GPCRs expressed one per cell in HEK293T, surface read by anti-HA immunostaining on unfixed cells. Class A and olfactory receptors carry HA immediately after Met1 with no leader; `_SP_HA` rows (68, class B/C/F) carry it after their native signal peptide. |

## Who is on the rung

GPCR membership is joined on Ensembl gene ID through `gene_resolver`, never on symbol, from two
sources: the Tedman Canonical sheet (Receptor Name with its `_ENST…` suffix stripped), and
UniProt keyword KW-0297 (G-protein coupled receptor, reviewed human; snapshot
[`uniprot_kw0297_human.tsv`](uniprot_kw0297_human.tsv)). Orthologs classify through the human
gene the manifest records for them.

| Members | library_reference | isoform | ortholog | Total |
|---|---:|---:|---:|---:|
| GPCR, in the Tedman screen | 462 | 5 | 22 | **489** |
| GPCR, not in the screen | 6 | 0 | 0 | **6** |
| Not a GPCR | 54 | 9 | 5 | **68** |
| **Total** | 522 | 14 | 27 | **563** |

A symbol join finds 486 in the screen. The Ensembl join adds three that the sheet names by an
older symbol: ACKR5 (`GPR182`, 6,317 a.u.), OR14L1 (`OR14L1P`, 769) and OR2I1 (`OR2I1P`,
132). The six GPCRs absent from the screen are ADGRF3, NPY4R, OPN1MW, OR51C1, OR8G3 and
TAS2R45. Of the 462 screened genes counted once through their library reference, 256 are class
A, 175 olfactory, 21 taste, 4 class B, 3 class C and 3 vomeronasal.

The 68 non-GPCRs are SLC5 (15 members), ABCC1/2/3/6/10, SLC10A1/2/6, SLC39A1/2/3/12,
SLC9A8/C1/C2/D1, SLC20A1/2, TTYH1/2/3, POPDC1/3, LMBR1L/LMBRD2, TLCD1/2, MFSD2A, SMIM7,
TMEM8B/170B/184A/198/205, and single members of other families (AQP12B, BEST1, GPM6B, KCNMA1,
KCNU1, MYMK, NIPA2, NIPAL4, PLPP2, PTCH1, SLC22A17, SLC31A1, SLC61A1, STRA6, TPRA1, TSPO2). The Tedman screen covers none of them,
so for these 68 the rung rests on topology alone.

N-tails are short for both groups: median 27 aa, interquartile range 24–40 for GPCRs and 9–39
for non-GPCRs. The tagged tail therefore has a median of 78 aa.

## What the evidence covers

### Tedman: the position, at 9–11 residues, without a leader

Tedman's class A and olfactory constructs put HA immediately after Met1, the junction this rung
uses, and read surface display without any exogenous leader. So the position is tested for the
489 screened GPCRs; the cassette size is not. HA with its start context is roughly 9–11
residues; the cassette is 51.

The screen also says where leaderless display stops working. Within class A, a longer native
N-tail goes with lower surface signal (Spearman ρ = −0.23, p = 1.5 × 10⁻⁴, n = 256 genes, one
library-reference member per gene); olfactory (ρ = −0.06, n = 175) and taste (n = 21) show no
trend over their narrow tail range. The trend is carried by the longest tails rather than spread
evenly:

| Native N-tail (aa) | Genes | Median Tedman intensity (a.u.) |
|---|---:|---:|
| 1–20 | 58 | 1,111 |
| 21–30 | 219 | 1,011 |
| 31–45 | 118 | 876 |
| 46–70 | 53 | 1,113 |
| 71–100 | 10 | 928 |
| > 100 | 4 | 164 |

Unit: gene (library-reference member); all screened classes; screen median across these 462 is
1,004 a.u. Of the five receptors with native tails of 98 aa or more, four read 121–199 a.u.:
ADRA1D (98 aa, 125), GPR135 (102, 198), CNR1 (116, 130) and GLP2R (173, 121). The fifth, CRHR2
(114, 2,385), starts with a pseudo signal peptide that is uncleaved and stays in the mature
receptor ([Rutz et al., 2006](https://pubmed.ncbi.nlm.nih.gov/16766521/)), so its tail is not
a plain leaderless tail.

That is the length bound this note uses: Tedman's leaderless HA constructs display at a typical
level with tagged tails up to about 110 residues (native ≤ 100 plus HA). With the 51-aa cassette,
any native tail of **60 aa or more** produces a tagged tail beyond that bound.

### The published mechanism

Leaderless GPCRs rely on TM1 to translocate their N-tail, and the literature ties failure of that
step to tail length and folding:

| Finding | Source |
|---|---|
| Across GPCRs, extracellular N-tails without a cleavable signal peptide are considerably shorter and carry fewer positive residues than tails preceded by one, suggesting different translocation mechanisms. | [Wallin & von Heijne, 1995](https://pubmed.ncbi.nlm.nih.gov/8577697/) |
| CB1's 116-aa N-tail is not efficiently translocated without a signal peptide, and the receptor is degraded by the proteasome. Adding a signal peptide or shortening the tail greatly enhances stability and surface expression without changing ligand binding. | [Andersson et al., 2003](https://pubmed.ncbi.nlm.nih.gov/12920192/) |
| Converting β2AR from a leaderless (type IIIb) to a signal-peptide (type IIIa) protein by adding an N-terminal cleavable signal sequence enhances its insertion into the ER membrane and expression of functional receptor. | [Guan, Kobilka & Kobilka, 1992](https://pubmed.ncbi.nlm.nih.gov/1331042/) |
| ETB without its native signal peptide is retained in the ER and nonfunctional; the signal peptide is required to translocate its N-tail. The authors propose signal peptides are needed where post-translational translocation of the N-terminus is impaired by stably folded domains. | [Köchl et al., 2002](https://pubmed.ncbi.nlm.nih.gov/11854280/) |
| A cleavable Lucy signal peptide raised olfactory receptor surface expression in HEK293T from 3 to 7 of 15 receptors without chaperones, and was cleaved from the mature protein without altering ligand binding. | [Shepard et al., 2013](https://pubmed.ncbi.nlm.nih.gov/23840901/) |

Every row is consistent with a leader helping and with harm being rare, but each tested a short
epitope or no tag. None tested a 51-aa cassette, and none compared leader against no leader on
a receptor that already displays well.

### The panel arms: what they decide

Eight GPCR arms are ordered and their maps are in `processed://constructs/ordered/`. Every one
carries the 51-aa cassette after the receptor's own Met1; the Igκ arms add the leader and Asp in
front of that Met1, so within each gene the leader is the only difference.

| Gene | Native N-tail (aa) | Tedman (a.u.) | Leader-free | Igκ |
|---|---:|---:|---|---|
| CNR1 | 116 | 130 | p0038 | p0041 |
| CCR5 | 30 | 302 | p0042 | p0043 |
| OR51E2 | 23 | 923 | p0044 | p0045 |
| ACKR3 | 45 | 2,200 | p0046 | — |
| AVPR2 | 38 | 6,367 | p0047 | — |

No display result for these arms is recorded in this repository as of 2026-11-10. When it
arrives, it settles three separate questions:

| Readout | What it decides |
|---|---|
| AVPR2 and ACKR3 leader-free | Whether the cassette displays at all after Met1 on receptors that express well with short tails. If both fail, the problem is the cassette, not the leader, and no leader choice fixes the rung. The cassette stays the full 51 aa, leading GGGGS included (Becca, 2026-11-11, "keep full"), so the next test would be its position, not its size. |
| CCR5 and OR51E2, leader-free vs Igκ | Whether the leader changes display inside the range Tedman already covers (tails 23–30, tagged 74–81). Igκ ≈ leader-free means no reason to add it to the 513 members without a flag below. Igκ well above leader-free means the 51-aa load matters even on short tails, and Igκ becomes the general default. |
| CNR1, leader-free vs Igκ | Whether a long tail needs the leader once the cassette is on it. This is the test for the 31 long-tail members below. |

Two things the panel cannot answer. No non-GPCR is in it, so the 68 non-GPCRs would inherit
whichever rule the GPCRs earn by analogy only. And no Igκ arm sits on a high expressor (AVPR2
and ACKR3 have none), so whether a leader ever lowers display on a receptor that is fine without
one stays untested.

### One arm per gene for the expression-rescue genes

24 of the L1's 25 expression-rescue genes belong to this rung (ADGRD2 is the exception). Each
rescue member was the Igκ leader fused to a GPCR that Tedman barely detects, meant to pair with
that gene's leader-free library reference as a within-gene test of the leader. The L1 is pooled
as one library, and Becca decided on 2026-11-13 that each gene keeps one member: "we are pooling
the whole L1 library together. why resuce pairs no i think we should make a single choice for
each L1". The per-gene choice is in "Expression-rescue genes: one arm each" below. The library
therefore carries no within-gene leader contrast, and the leader's effect stays untested in L1.

## Five subgroups, by what the leader decision turns on

Each member is assigned by [`nterm_res1_multipass.py`](nterm_res1_multipass.py), first rule
that matches. UniProt features (snapshot
[`uniprot_nterm_res1_features.tsv`](uniprot_nterm_res1_features.tsv)) describe the canonical's
first residues, so they count only when the member's first 20 residues are the canonical's: 51
members start elsewhere, among them the PTCH1 isoform `l1-p2491`, which begins at canonical
residue 152, after TM1 (100–117).

| Subgroup | Rule | GPCR | non-GPCR | Total |
|---|---|---:|---:|---:|
| 1 Orientation disputed | non-GPCR whose UniProt canonical starts cytoplasmic | 0 | 5 | 5 |
| 2 Possible signal peptide at residue 1 | UniProt or DeepTMHMM v2 calls one | 3 | 3 | 6 |
| 3 Short tail, cassette erodes TM1 charge bias | tail < 15 aa, Δ > 0 natively and ≤ 0 with the cassette | 4 | 4 | 8 |
| 4 Long tail | native tail ≥ 60 aa (tagged ≥ 111) | 25 | 6 | 31 |
| 5 No member-specific flag | none of the above | 463 | 50 | 513 |

### Can a leader flip TM1?

Not relative to the predicted topology. The leader is cleaved in the ER lumen and takes the
cassette and N-tail with it, so TM1 only has to stop transfer with its N-terminus out, which is
the orientation this rung already assumes. What a leader does is lock that orientation in. It
inverts a protein only where the prediction is wrong (subgroup 1), or where a native signal
sequence the prediction missed is pushed into the chain behind it and becomes an extra
hydrophobic segment (subgroup 2).

The leader-free cassette is the version that can flip TM1. It lengthens the tail TM1 must
translocate after it has inserted and adds net +3 to it. Where the native tail is shorter than
15 residues the cassette's charges also enter TM1's N-terminal flank (subgroup 3). Among the
non-GPCRs, that is TSPO2, MYMK, GPM6B and PLPP2; none of the 68 is a single signal-anchor
protein, since every one has at least two predicted helices.

### Subgroup 1: orientation disputed (5, all non-GPCR)

| Member | Gene | Tail (aa) | Tail net charge | Δ native | Evidence against an extracellular N-terminus |
|---|---|---:|---:|---:|---|
| `l1-p1213`, `l1-p2466` | MFSD2A | 46 | +5 | −6 | v1 calls 12 helices with an extracellular N-terminus, which puts the C-terminus outside too. The mouse cryo-EM structure has its post-TM12 helix on the intracellular surface ([Wood et al., 2021](https://pubmed.ncbi.nlm.nih.gov/34349262/)), and with an even helix count both termini then sit inside. The tail (`…QVKKEPKKKKQQ…`) and both flanks favor a cytoplasmic N-terminus. |
| `l1-p1988` | SLC5A4 | 28 | −3 | +1 | UniProt, by sequence analysis only. |
| `l1-p1948` | SLC39A2 | 8 | 0 | +2 | UniProt, by sequence analysis only. |
| `l1-p0128` | AQP12B | 7 | 0 | +3 | UniProt, by sequence analysis only. |

The literature review resolves four of the five. MFSD2A is mis-placed: a second structure, of
chicken MFSD2A, also has twelve helices ([Cater et al., 2021](https://pubmed.ncbi.nlm.nih.gov/34135507/)),
and the review places its glycosylated extracellular domain in the TM5–TM6 loop, which v1 calls
cytoplasmic. A tag after Met1 would be cytoplasmic, and a leader would invert the transporter.
SLC5A4 (SGLT3) shares SGLT1's 28-aa tail and 14-helix fold, and SGLT1's N-terminus is
extracellular by glycosylation scanning ([Turk et al., 1996](https://pubmed.ncbi.nlm.nih.gov/8567640/)).
For SLC39A2, the review reports that the live UniProt entry now calls residues 1–8
extracellular. AQP12B stays unresolved: v1 adds an N-terminal helix that the aquaporin fold
does not have.

### Subgroup 2: possible signal peptide at residue 1 (6)

| Member | Gene | Tail (aa) | Call | Tedman (a.u.), HA after Met1 |
|---|---|---:|---|---:|
| `l1-p0724` | GLP2R | 173 | v2: 66-aa signal peptide | 121 |
| `l1-p0475` | CRHR2 | 114 | UniProt: 1–19 | 2,385 |
| `l1-p1484` | OR52L1 | 44 | v2: 15-aa signal peptide | 1,018 |
| `l1-p2144` | TLCD1 | 7 | UniProt: 1–35 | — |
| `l1-p2054`, `l1-p2518` | SMIM7 | 4 | UniProt: 1–17 | — |

Tedman built CRHR2 and OR52L1 from the same transcripts we order, with HA after Met1, and both
displayed (2,385 and 1,018 a.u.). For CRHR2 that matches Rutz et al.'s finding that the segment
is an uncleaved pseudo signal peptide; the cyno ortholog `l1-p2558` (since dropped from the panel) carries the same N-terminus
and is assigned to subgroup 4 only because its UniProt entry is not the human one. So both stay
on the rung. Under the GPCR Igκ default, CRHR2's pseudo signal peptide would follow the leader and
the cassette as an internal hydrophobic segment. Rutz et al. show the segment cannot target the
ER on its own, but whether it inserts as an extra helix behind a leader is untested. GLP2R
displayed poorly in the same construct and has since left this rung for N after SP.

For TLCD1, the review found a leader-free N-terminal precedent. In the chicken ortholog
calfacilitin, a myc tag after residue 3 stained without detergent in HEK293T, so the first
hydrophobic segment is an uncleaved N-out TM1 rather than a signal peptide
([Papanayotou et al., 2013](https://pubmed.ncbi.nlm.nih.gov/23673622/), from the full text, as
reported by the review). For SMIM7 there is no literature. The KLKKK cluster right after its
hydrophobic segment favors a signal anchor over a cleaved signal peptide, and Igκ puts the
cassette outside under either reading.

### Subgroup 3: short tail, cassette erodes TM1 charge bias (8)

GPR119, GPR26, VN1R4, TAS2R4 (GPCRs), and TSPO2, MYMK, GPM6B (`l1-p2444`, isoform) and PLPP2
(`l1-p2486`, isoform). Native tails are 1–6 aa; with the cassette, Δ drops to between 0 and −2.
Δ alone is a weak predictor here: 177 of the 495 GPCRs already have Δ ≤ 0 natively and still
insert N-out, Tedman displayed GPR26 (1,099) and VN1R4 (1,158) with HA after Met1, and the panel's
shortest tail is OR51E2 at 23 aa, so no panel arm tests this case. A leader removes the question
for these eight, because TM1 no longer chooses its own orientation.

### Subgroup 4: long tail (31)

25 GPCRs and 6 non-GPCRs with native tails of 60 aa or more, listed in full in the per-member
table. Three are expression-rescue genes (CNR1, ADRA1D, GPR142), now carried by their Igκ + tag
reference alone. The other members of note:

| Member | Gene | Tail (aa) | Note |
|---|---|---:|---|
| `l1-p0067` | ADGRF3 | 707 | The ordered proteoform starts `MTTRKLSAHS…`, not at the canonical's 25-aa signal peptide, so it has no signal peptide and a 707-aa leaderless ectodomain. |
| `l1-p2491` | PTCH1 (isoform) | 285 | Starts at canonical residue 152, so it lacks TM1 and its N-tail is the first extracellular domain. |
| `l1-p2506` | SLC39A12 (isoform) | 233 | |
| `l1-p2205` | TMEM8B | 229 | |
| `l1-p0448` | CNR1 | 116 | Andersson's case; panel pair p0038/p0041. |

Leaderless translocation of the first four, with or without a cassette, is the case Köchl et al.
describe as needing a signal peptide. These are the members where a leader is most likely to
matter.

### Subgroup 5: no member-specific flag (513)

463 GPCRs, of which 458 are in the Tedman screen, and 50 non-GPCRs. For the GPCRs the position
is tested and the cassette size is not; 52 of them have tails of 46–59 aa, between the panel's
longest leader-free tail (ACKR3, 45) and the length bound. For the 50 non-GPCRs nothing beyond
topology bears on the question.

## Non-GPCRs: the literature review

Three reviewers covered the 68 non-GPCRs in parallel, one per group: the transporters (34), the
ABC transporters, channels and others (20), and the small, poorly studied proteins (14). Each
looked for experimental evidence of N-terminal orientation, a missed or mis-called signal
peptide, and any precedent for N-terminal tagging. Each then gave every member one verdict:
`igk`, `leader_free`, `off_rung` or `hold`. The table, with the evidence and PMIDs for every
member, is [`nterm_res1_nongpcr_review.tsv`](nterm_res1_nongpcr_review.tsv). The PMIDs cited
on this page were re-resolved in PubMed and checked against their abstracts.

| Verdict | Members | Confidence high / medium / low |
|---|---:|---|
| `igk` | 63 | 17 / 36 / 10 |
| `off_rung` | 2 | 2 / 0 / 0 |
| `hold` | 3 | 0 / 1 / 2 |
| `leader_free` | 0 | — |

No member came out `leader_free`, and the reason is the same throughout. Where the N-terminus
is reliably extracellular and residue 1 does not start a signal peptide, a leader has no
identified cost, and it removes the two leaderless risks: long tails that fail to translocate,
and short tails whose TM1 charge bias the cassette erodes. Leader-free tags have worked for some
of these families, though only small epitopes were tested:

| Protein | Leader-free N-terminal tag | Source |
|---|---|---|
| MRP1 (ABCC1) | HA at residue 4: functional, and stained on intact cells | [Kast & Gros, 1997](https://pubmed.ncbi.nlm.nih.gov/9334225/) |
| SGLT1 (SLC5A1) | His6 at residue 5: transport kept, and detected on non-permeabilized cells | [Lin et al., 1999](https://pubmed.ncbi.nlm.nih.gov/10441667/) |
| BK α (KCNMA1) | epitope-tagged constructs that place the N-terminus outside the cell | [Meera et al., 1997](https://pubmed.ncbi.nlm.nih.gov/9391153/) |
| Calfacilitin (TLCD1 ortholog) | myc after residue 3, stained without detergent | Papanayotou et al., 2013 (above) |

Structures and other experimental evidence put the N-terminus outside for the TTYH family
([Sukalskaia et al., 2021](https://pubmed.ncbi.nlm.nih.gov/34385445/)), NIPA2
([Karuppan et al., 2026](https://pubmed.ncbi.nlm.nih.gov/42481511/)) and TPRA1, whose N-terminal
glycans bind oncolytic virus M1 ([Hu et al., 2025](https://pubmed.ncbi.nlm.nih.gov/41053536/)).
Several members rest on family homology or prediction alone, and they carry `low` confidence
in the table: SLC9C1/C2, LMBR1L, LMBRD2, TSPO2, TMEM205, SMIM7, and the BEST1, PLPP2 and GPM6B
isoforms. The last three start at an internal methionine after canonical TM1, so they may not
fold without it.

| Member | Gene | Verdict | Why |
|---|---|---|---|
| `l1-p1213`, `l1-p2466` | MFSD2A | `off_rung` | Cytoplasmic N-terminus (subgroup 1). Its displayable site is the TM5–TM6 extracellular domain, which is an internal-loop option. |
| `l1-p0128` | AQP12B | `hold` | v1's N-terminal helix against the aquaporin fold; AQP12 stays intracellular when expressed. A glycosylation test of its NVS sequon at residues 5–7 would settle it. |
| `l1-p1994` | SLC61A1 | `hold` | Major facilitator superfamily (MFS) transporters are usually N-in, and v1's N-out call rests on a 13-residue TM1. The AlphaFold model or a surface-staining pilot would settle it. |
| `l1-p2516` | SLC9D1 (isoform) | `hold` | Starts at canonical residue 231 with no Met (`ITKREV…`), so "after Met1" is undefined. The canonical belongs on the "N after SP" rung; swap the member or confirm the transcript. |

The verdicts are recommendations, and none changes a placement in `l1_tag_paths.tsv`.

## Is the lane right for every member?

[NTERM_RES1_LANE_CHECK.md](NTERM_RES1_LANE_CHECK.md) tests each of the 561 members against the
lane's three conditions (extracellular residue 1, no cleavable signal peptide, at least two
helices), with DeepTMHMM v1 as the call of record and UniProt, DeepTMHMM v2, SignalP 6 and the
member's structure model as cross-checks. The lane holds for 560 and TSPO2 is unresolved. Its
earlier moves (MFSD2A ×2, AQP12B, GLP2R) are applied on #283 and off the rung, cyno NIPAL4 is
restarted at Met215 and holds, and SLC61A1 stays on the lane as the writer decided on 2026-11-13.

## Decisions taken under Becca's delegation

On 2026-11-12 Becca delegated the protein-level calls: "use your best judgment to complete a
new and improved evidence based l1 completely by the morning." This session made the calls
below on the evidence in the sections above. Each is recorded per member with its rejected
alternative, so any of them can be reversed.

| Call | Members | Evidence | Rejected alternative |
|---|---|---|---|
| Igκ + tag for the non-GPCRs | 58 | Literature review: N-out established where tested, and the leader removes the leaderless-translocation and TM1-flip risks at no identified cost. Five entrants without a review follow their family. | Tag alone |
| Igκ + tag for the rescue genes' references; rescue members retired | 23 | One arm per gene (Becca, 2026-11-13). The GPCR default makes the reference the same protein as its rescue twin, so the leader is kept and only the duplicate goes. | Keep the reference tag-alone beside its twin (no pairs in a pooled L1) |
| Tag alone for CRHR2 | 1 | Uncleaved pseudo signal peptide (Rutz et al., 2006); Tedman displayed it leader-free at 2,385 a.u. | Igκ + tag |
| Lane moves | 4 | MFSD2A: two structures. AQP12B: the aquaporin fold. GLP2R: a v2-called signal peptide. | Keep on this rung with Igκ |
| SLC61A1 kept on this rung | 2 | Decided by the manifest's writer (2026-11-13): DeepTMHMM v1 and v2, SignalP 6 and the AlphaFold model (13 crossings, pLDDT 89) all put the N-terminus outside, and the lane check rates the lane as holding. This session had argued for internal (the major facilitator N-in fold and a 13-residue TM1); that view is the rejected alternative. | Move to internal |
| Replace the member | 2 | Cyno NIPAL4 has a 214-residue Arg/Lys/Gly-rich extension before Met215, which v2 calls cytoplasmic. The OPN3 isoform keeps helices 1–4 only, and canonical OPN3 is already on this rung. | Build with Igκ + tag |
| Routed in from "N after SP": OR11H6, TMEM63B | 2 | Only DeepTMHMM v1 calls a signal peptide (residues 1–10) on either. SignalP 6 calls none, and DeepTMHMM v2 and UniProt put the N-terminus outside ahead of the first helix: OR11H6 1–43 extracellular with TM1 at 44–64 (a class A olfactory GPCR; Tedman displayed it tagged after Met1 at 1,075 a.u.); TMEM63B 1–40 extracellular with TM0 at 41–65 in the cryo-EM structure (PDB 8EHX; [Zheng et al., 2023](https://pubmed.ncbi.nlm.nih.gov/37543036/)). The Igκ swap of residues 1–10 deleted native extracellular residues, so both are built Igκ + tag at Met1 on the full native sequence (`lane_overrides.tsv`). Proposed by the N-terminal signal-peptide audit. | Keep on "N after SP" with the Igκ swap |
| Routed in from "C-terminal": SLC41A2 | 1 | Experimental topology ([Sahni et al., 2007](https://pubmed.ncbi.nlm.nih.gov/16984228/); UniProt Q96JW4, every domain ECO:0000269) puts residues 1–162 outside and the C-terminus (565–573) inside, the reverse of DeepTMHMM v1 and v2; the helix positions agree, only the orientation differs. A C-terminal tag would be cytoplasmic. Tag alone would leave TM1, whose flanks carry no charge bias (2 K/R each side), to translocate a 213-residue tail, so Igκ + tag at Met1 is the form that delivers the experimental orientation. Medium confidence: one study against both predictors; if they are right, the leader inverts the transporter. | Keep C-terminal; tag alone; an internal tag in the 428–469 loop |
| Routed in from "N after SP": human OR10G6 | 1 | Tedman displayed it with HA right after Met1 and no leader at 1,029 a.u., so its N-terminus is not cleaved; SignalP 6 and UniProt call no signal peptide (UniProt: 1–46 extracellular, TM1 47–67). DeepTMHMM v1 and v2 both call one on a 22-residue hydrophobic extension that mouse Or10g6 lacks. Built Igκ + tag as for OR11H6 (class A olfactory, GPCR default). Medium confidence. Becca's ortholog rule (2026-11-15: no ortholog built without a built human member) was the trigger: with the human built, mouse `l1-p3043` stays built. | Hold on "N after SP" (which would also hold mouse `l1-p3043`); tag alone, the form Tedman displayed |
| OR8G3 start recoded to ATG | 1 | A GTG start initiates poorly in a construct. | Keep in curation, unbuilt |

The single writer (l1-reviewer-package) applied the moves for MFSD2A ×2, AQP12B and GLP2R and the OR8G3 recode on PR #283, under its own delegated authority (`lane_overrides.tsv`). The proteoform audit dropped the OPN3 isoform. #283 also restarted cyno NIPAL4 at Met215 (dropping the 214-residue extension) and cyno CHRM3 and OR3A2 at their new Met1, so all three are built on the new sequences. SLC61A1 stays on this rung (the writer's decision, 2026-11-13). The GPCR default
remains Becca's own decision of 2026-11-10.

## Built arms

[`build_l1_multipass_arms.py`](build_l1_multipass_arms.py) builds every member of the rung. It
uses the layout of the type III and C-terminal builders: an L1 Esp3I wrapper, Kozak, an optional
Igκ leader plus Asp, the member's Met1, then the full 51-aa cassette on p003's codons, then the
native ORF. Outputs: [`constructs_l1_multipass.tsv`](constructs_l1_multipass.tsv), one row per
construct with its `status` and attribution; `constructs_l1_multipass.fasta` (inserts); and
`multipass_arms_mature.fasta` (mature proteins, leader removed, for refolding). Each construct
is keyed `<design_id>:<design>`.

| Design | Members | Who | Decided by |
|---|---:|---|---|
| Igκ + tag | 502 | GPCRs, including the 23 expression-rescue genes' references, OR11H6 (routed in 2026-11-14) and human OR10G6 (2026-11-15) | Becca, 2026-11-10 (GPCR default); one arm per gene, Becca 2026-11-13 |
| Igκ + tag | 62 | non-GPCRs: 55 with an `igk` review verdict (TMEM63B's and SLC41A2's added 2026-11-14), 5 entrants without a review (including cyno NIPAL4, restarted at Met215), and human and mouse SLC61A1 (kept on this rung by the writer) | this session, delegated; SLC61A1 by l1-reviewer-package |
| Tag alone | 1 | human CRHR2 | this session, delegated |
| Moved to another lane | 0 | MFSD2A ×2, AQP12B and GLP2R have already left the rung (`lane_overrides.tsv`) | this session, delegated; applied by l1-reviewer-package |

For each built construct, the build checks the encoded protein exactly, a single stop, and only
the two wrapper Esp3I sites, both alone and in the p0039 backbone. Rebuilding with `--check`
reproduces the record. ORFs come from the constructs record pinned in
`library_design.CONSTRUCTS_KEY`. OR8G3's GTG start is recoded to ATG in its arm, with no added
residue, matching the ATG fixes for SLC22A17 and TRPV6. Each row records its
`rejected_alternative`, so any call can be reversed.

## Expression-rescue genes: one arm each

Becca's rule (2026-11-13) is one member per gene in the pooled L1. The choice for each of the 25
genes is this session's, under her 2026-11-12 delegation, and is recorded gene by gene in
[`rescue_one_arm.tsv`](rescue_one_arm.tsv), written by [`rescue_one_arm.py`](rescue_one_arm.py).
The rescue members leave the order, so the L1 goes from 2,735 active members to 2,710; the
writer retires them in the manifest.

**Rule: each gene keeps its library reference, built with its lane's default design, and its
rescue member is retired.** For the 23 genes on this rung the default for a GPCR is Igκ + tag,
which is exactly the rescue arm's protein (Igκ + Asp + Met1 + cassette + receptor residues
2..end). So the leader the rescue set was chosen to test is kept, under the reference's
`design_id`, and only the duplicate goes. Refold cannot separate the two: once the leader is
cleaved they are the same protein, so the retired twin's verdict is the reference's. 19 pass;
CNR1 (tail only), GNRHR (epitope packed), TAS2R3 and TAS2R5 (fold disturbed) fail, as they did
leader-free.

| Gene | Kept | Its design | Refold | Retired | Why |
|---|---|---|---|---|---|
| 23 GPCRs on this rung | reference | Igκ + tag | 19 pass, 4 fail | rescue member (same protein) | the rule |
| ADGRD2 | `l1-p0061` | tag after its native 18-residue signal peptide (N-after-SP lane) | pass | `l1-p2361`, Igκ replacing the native signal peptide (pass) | Both pass. Nothing shows the native signal peptide underperforms Igκ, so the lane's design for a native signal peptide stands. |
| GLP2R | `l1-p0724` | Igκ replacing the signal peptide that DeepTMHMM v2 calls at 1–66 (N-after-SP lane) | pass | `l1-p2365`, Igκ + tag at Met1 (pass) | The rescue arm keeps residues 1–66 behind the tag, so a cleaved signal peptide would take the tag with it. Tedman's HA-after-Met1 construct read about 0.12× the median, as that would predict. |

Each row of `rescue_one_arm.tsv` names the rejected arm and what would change the call: for the 23,
evidence that the leader lowers display of a poorly displayed GPCR; for ADGRD2, evidence that its
own signal peptide targets poorly; for GLP2R, evidence that its N-terminal segment is not cleaved.

[`build_l1_rescue_arms.py`](build_l1_rescue_arms.py) and its tables (`constructs_l1_rescue.*`,
`rescue_arms_mature.fasta`) stay as the record of the retired members and are not ordered. The
within-pair barcode-swap risk and the linker and receptor recodes drawn against it no longer apply,
because no gene has two members.

## Refold results

l1-tag-refold folded every built arm against its native model: ColabFold model_1_ptm, 3
recycles, the native MSA reused, and its thresholds in `score_tagged_folds.py`. All 565 built arms
are scored: 493 pass and 72 fail. OR11H6 and TMEM63B, routed in on 2026-11-14, both pass
(native-region lDDT 0.983 and 0.967; SpyTag003/ALFA RSA 0.57/0.52 and 0.71/0.74). The refold session also folded the order of record, and its
verdicts for the 561 arms it covered are the same, 561 of 561. SLC41A2 and human OR10G6, routed in
later, both pass (native-region lDDT 0.993 and 0.974; SpyTag003/ALFA RSA 0.75/0.78 and
0.71/0.58), so OR10G6's tag-alone fallback is not needed.
[`multipass_refold_dispositions.py`](multipass_refold_dispositions.py) classifies each
failure and records a disposition per arm in
[`multipass_refold_dispositions.tsv`](multipass_refold_dispositions.tsv). It stops if the
verdicts do not cover exactly the built arms, or if any tagged sequence differs from the build.

**Failures cluster by receptor family, not by design.**

| Group | Fail / scored | Fail rate |
|---|---:|---:|
| Igκ + tag | 71 / 564 | 13% |
| Tag alone (CRHR2) | 1 / 1 | — |
| Olfactory GPCRs | 3 / 186 | 2% |
| Taste GPCRs | 9 / 22 | 41% |
| Other GPCRs | 57 / 295 | 19% |
| Non-GPCRs | 3 / 62 | 5% |

Unit: tagged arm. Failure is the refold session's verdict.

| Failure class | Arms | What it is | Disposition |
|---|---:|---|---|
| Buried SpyTag003 | 42 | SpyTag003 called buried (RSA < 0.25), mostly in a stretch the model is unsure of (SpyTag003 median pLDDT 40), with fewer than 5 confidently packed epitope residues. ALFA stays exposed. | Keep as built, flagged: read ALFA as the display control |
| Epitope packed | 13 | At least 5 epitope residues confidently packed on the receptor: GALR2, CX3CR1, HCRTR1 ×2, MRGPRX2, BDKRB1, MCHR2, NTSR2, PTGDR, P2RY4, GRPR, C5AR2, CCR4 (class A peptide and chemokine receptors) | Keep as built, flagged: read ALFA as the display control |
| Fold disturbed | 11 | The native region moved (lDDT < 0.85): seven taste receptors with 6–8-residue tails (TAS2R1, -5, -7, -8, -10, -14, -19), VN1R4, GPBAR1, MRGPRX4 and TMEM170B. The cassette sits next to TM1. | Keep as built, flagged: a 10-residue spacer rescued none of the 12 tested |
| Tail only | 6 | Only the member's own disordered N-tail next to the cassette moved; native region intact (lDDT ≥ 0.85): AGTR1, CNR1, OR5C1, OR8G3, OR9A2, TAS2R3 | Keep as built, flagged |

**Rule for all 72 failures: keep as built and flagged** (this session, delegated by Becca
2026-11-12; recorded 2026-11-14). The choice is between keeping an arm and holding it, because
none of the alternatives tested does better:

- The cassette size is fixed ("keep full").
- A 10-residue spacer between the cassette and native residue 2 rescued none of the 12 arms first
  classed fold-disturbed (`verdicts_multipass_spacer.tsv`; every one still fails on native-region
  lDDT). A 30-residue spacer made the snorkel arms worse.
- The ranker cannot vet an internal site for any of these members. None has a membrane-side
  mask or per-residue features, because none is on the internal rung, and GPCR loop sites would
  mostly fall back to snorkels.

No arm has a specific reason to hold: every verdict is a single-model prediction, and the pooled
screen measures display directly. The leader does not separate passes from failures: the four
rescue genes that fail (CNR1, GNRHR, TAS2R3, TAS2R5) fail both with and without it. Two of them
change class with the leader. GNRHR's Igκ arm has SpyTag003 buried without confident epitope
packing (epitope packed leader-free), and TAS2R3's native region is intact by a margin of 0.007
(tail only; fold disturbed leader-free).

Each row of `multipass_refold_dispositions.tsv` records the pattern, the metric it rests on, the
disposition, the rejected alternative, the bench result that would change it, decided_by,
implemented_by and the date. The bench result is a non-permeabilized stain of each arm:
- ALFA stains but SpyCatcher003 does not: a buried or packed SpyTag003 is real.
- Neither stains while total expression is present: a fold-disturbed or tail-only arm is
  misfolded or retained.

### Open design question for Becca: SpyTag003 buried at the mature N-terminus

In 53 of the 72 failing arms the refold calls SpyTag003 buried: the 42 "buried SpyTag003" arms,
10 of the 13 "epitope packed" and one fold-disturbed. They are almost all class A peptide and
chemokine GPCRs (50; the others are human and mouse SLC61A1, `l1-p1994` and `l1-p3105`, and 1
taste receptor). Their native N-tails run 6–114 aa (median 40). The call is weak for most: 48 of the 53
have SpyTag003 pLDDT below 50 (median 40.5), and only the 13 epitope-packed arms have five or
more epitope residues confidently packed. ALFA stays exposed in every one of the 53 (RSA
0.30–0.77, median 0.62; buried is < 0.25), so ALFA can serve as each arm's display control and
SpyCatcher003 capture is the readout at risk. They are kept as built.

The open question is whether L2 should move SpyTag003 for this receptor class: put ALFA first in
the cassette, or move SpyTag003 to the C-terminus or a loop. Or should we wait for the L1
ALFA-versus-SpyCatcher003 ratio on these 53 to show whether the burial is real?

## Per-member table

Attribution follows Becca's rule of 2026-11-10: note the agent behind every manual fix and its
reasoning. Nothing in this analysis overrides a value or uses a hand-picked member list. GPCR
status is computed from the screen and UniProt KW-0297, so there is no hand-added list of
unscreened GPCRs. Three of the nine GPCRs once counted by hand as missing from the screen are in
it under older symbols, as listed above. The subgroup rules and the cutoffs 15, 20 and 60 are
analysis choices made in this session and not yet reviewed by a person. Every row carries
`leader_subgroup_decided_by`, `leader_subgroup_implemented_by` and `leader_subgroup_reason`. The
GPCR leader default is Becca's decision. The non-GPCR verdicts come from a literature review
that cannot be regenerated by code, so the review table is committed with its PMIDs, and each of
its rows names who produced it and that it awaits Becca's confirmation.

[`nterm_res1_multipass.tsv`](nterm_res1_multipass.tsv), one row per member (565). Columns:
identity (`design_id`, the stable key across renumbering; `candidate_id`, `candidate_set`, `species`, `gene_symbol`, `ensembl_gene`); GPCR
evidence (`in_tedman`, `tedman_class`, `tedman_au`, `uniprot_kw0297`, `gpcr`, `group`); TM1
features from the v1 string (`n_tail_aa`, `tm1_start`, `tm1_end`, tail charges and N-glycosylation
sequons, `delta_charge_native`, `tagged_n_tail_aa`, `tagged_n_tail_net`, `delta_charge_tagged`);
cross-checks (`v2_signal_peptide_aa`, `uniprot_accession`, `uniprot_same_start`,
`uniprot_signal_end`, `uniprot_n_terminus`, `uniprot_n_terminus_evidence`); `leader_subgroup` with
its attribution; `planned_leader` with `planned_leader_decided_by`, `planned_leader_implemented_by`
and `planned_leader_reason`; `cassette` with `cassette_decided_by`, `cassette_implemented_by` and
`cassette_reason`; and the provenance pair `topology_version`, `sequence_sha256`.

Only v2's signal-peptide call is read. Its other state letters are not documented in this
repository, so v2's membrane side is not compared.

## Reproducing

Topology is DeepTMHMM v1 throughout (`topology.TOOL_V1`), matched by the sha256 of each member's
designed sequence against the snapshot pinned in `topology.SNAPSHOT_KEY`. DeepTMHMM v2 is read
only for its signal-peptide call and never replaces v1. Every external input is snapshotted with
its version:

| Input | Pinned as | Per-row fingerprint |
|---|---|---|
| Topology | `topology.SNAPSHOT_KEY` | `topology_version`, `sequence_sha256` |
| UniProt features and KW-0297 list | `uniprot_nterm_res1_features.tsv`, `uniprot_kw0297_human.tsv`; release in `uniprot_nterm_res1_snapshot.json` (2026_04, released 2026-11-04) | `uniprot_release`, `uniprot_entry_version`, `uniprot_sequence_version` |
| Tedman screen | `external://tedman-gpcr-surface-screen/2025-09-19/tedman2025_gpcr_screen_media-2.xlsx` (immutable) | — |
| Rung membership | `l1_tag_paths.tsv`, `candidate_manifest_tranches.tsv` | `candidate_id` |

    uv run --locked python nterm_res1_multipass.py                # offline, from the snapshots
    uv run --locked python nterm_res1_multipass.py --check-drift  # has any input moved?
    uv run --locked python nterm_res1_multipass.py --refresh      # re-fetch UniProt and rerun

Every output on this rung can be checked against a rebuild without writing anything:

    uv run --locked python nterm_res1_multipass.py --check         # the per-member table
    uv run --locked python build_l1_multipass_arms.py --check      # every arm, byte for byte
    uv run --locked python multipass_refold_dispositions.py --check  # the refold dispositions
    uv run --locked python rescue_one_arm.py                       # one arm per rescue gene

`--check-drift` re-fetches only the UniProt entry and sequence versions, and recomputes each
member's sequence hash and v1 topology version. It lists every member whose inputs changed, or
whose rung membership changed, and exits 1 if any did; that is the signal to rerun with
`--refresh`. On 2026-11-11 it reported no changes. The script stops if any member's designed
sequence no longer hashes to the value in `l1_tag_paths.tsv`.

## Decision log

| Date | Decision | Reason | Decided by / implemented by |
|---|---|---|---|
| 2026-11-10 | Rung characterized; leader decision set out as options A–D and left to Becca, with option A recommended | `tag_cassette` is deferred and the panel pairs have not read out | no decision taken: the leader choice is Becca's / Claude Code session l1-nterm-res1-multipass, PR #280 |
| 2026-11-10 | UniProt N-terminal features applied only to members that start where the canonical does | The PTCH1 isoform starts after TM1, so the canonical's experimentally cytoplasmic N-terminus does not describe it | no decision taken: analysis choice, not yet reviewed by Becca / Claude Code session l1-nterm-res1-multipass, PR #280 |
| 2026-11-10 | Long-tail bound set at native ≥ 60 aa | Tagged tail then exceeds the ~110 residues Tedman shows displaying leaderless; four of five screened receptors with tails ≥ 98 aa read 121–199 a.u. | no decision taken: analysis choice, not yet reviewed by Becca / Claude Code session l1-nterm-res1-multipass, PR #280 |
| 2026-11-10 | Member IDs given in the `l1-mNNNN` scheme | The renumbering on `main` replaced `l1-pNNNN`; a rerun gave the same 563 members, values and subgroups | no decision taken (ID update after the renumbering on `main`) / Claude Code session l1-nterm-res1-multipass, PR #280 |
| 2026-11-10 | Attribution columns added to the per-member table and to this log | Traceability for every analysis choice | Becca (rule for all sessions, relayed by l1-recode-dups) / Claude Code session l1-nterm-res1-multipass, PR #280 |
| 2026-11-10 | Igκ leader by default for the rung's 495 GPCRs | Becca, on Aashish Manglik's suggestion | Becca / Claude Code session l1-nterm-res1-multipass, PR #280 |
| 2026-11-10 | Non-GPCR verdicts from a literature review: 63 Igκ, 2 off the rung (MFSD2A), 3 on hold | Experimental orientation, signal-peptide and tagging evidence per member | recommendation, not yet confirmed by Becca / literature-review subagents of l1-nterm-res1-multipass, PR #280 |
| 2026-11-10 | Compara 116 ortholog panel checked and found not yet on `main` | 11 of `main`'s 200 ortholog keys are in `ortholog_panel_compara116.tsv`; PR #283 is open | no decision taken / Claude Code session l1-nterm-res1-multipass, PR #280 |
| 2026-11-11 | Topology fixed to DeepTMHMM v1 by default; UniProt release and per-entry versions pinned per row, with a `--check-drift` mode | Know when a UniProt change requires a rerun | Becca / Claude Code session l1-nterm-res1-multipass, PR #280 |
| 2026-11-11 | Full 51-aa cassette kept, leading GGGGS included; no trimmed cassette on this rung | Without the GGGGS, v1 called a cleavable signal peptide ending inside SpyTag003 on 12 of 62 residue-1 tag-alone arms, against 0 of 62 with it | Becca ("keep full") / Claude Code session l1-nterm-type3, PR #282; recorded here by l1-nterm-res1-multipass, PR #280 |
| 2026-11-12 | Isoform cohort to be gated before the draw (native ATG start; canonical N-terminal orientation and TM1 kept, by DeepTMHMM v1), then redrawn under the same seed, instead of swapping flagged isoforms one at a time | Per-member rejection reshuffled the strata and drew new forms with the same defects (no Met, fragments, no designed CDS); affects this rung's PTCH1, SLC9D1, OPN3, GPR161 and GPM6B isoforms | Becca ("Gate cohort, redraw") / Claude Code session l1-reviewer-package, PR #283; recorded here by l1-nterm-res1-multipass, PR #280 |
| 2026-11-12 | Arms built: Igκ + tag for 490 GPCRs (decided); both forms for 61 non-GPCRs and leader-free alternates for 25 GPCRs (proposed); 12 held | Becca's overnight request to build the rung's arms, building both forms where only a proposal exists | Becca (overnight plan, 2026-11-11; GPCR default 2026-11-10) / Claude Code session l1-nterm-res1-multipass, PR #280 |
| 2026-11-12 | Tables and arms regenerated on PR #283's state (9928f4b1): 566 members; 501 Igκ + tag decided, 142 proposed, 6 held; rung changes listed by design_id | #283's redraw, Compara 116 adoption and reference fixes landed; requested by l1-integration | Becca (redraw, fixes, Compara 116) / Claude Code session l1-nterm-res1-multipass, branch claude/l1-multipass-rekey-283 |
| 2026-11-12 | One design per member: 536 Igκ + tag, 24 tag alone, 6 moved to another lane, OR8G3 recoded to ATG; the proposed-alternative arms dropped | Becca's rule change: decide protein-level calls on the evidence instead of queueing them | l1-nterm-res1-multipass (delegated by Becca 2026-11-12, "use your best judgment") / Claude Code session l1-nterm-res1-multipass, PR #280 |
| 2026-11-12 | SLC61A1 (human, mouse) kept on this rung with Igκ; cyno NIPAL4 and the OPN3 isoform marked for replacement and not built | The lane check's v2, AlphaFold, charge and SignalP evidence outweighs the homology argument for moving SLC61A1; NIPAL4's cytoplasmic N-terminal extension and OPN3's 4-helix fragment | l1-nterm-res1-multipass (delegated by Becca 2026-11-12) on the evidence of l1-nterm-lane-check, PR #296 / Claude Code session l1-nterm-res1-multipass, PR #280 |
| 2026-11-12 | SLC61A1 (human, mouse) moved to internal after all, superseding the row above | With orientation disputed between v1/v2/AlphaFold and UniProt plus the major facilitator fold, a gated internal site is robust to either answer; reconciled with the single writer, which applied it | l1-nterm-res1-multipass and l1-reviewer-package (both delegated by Becca 2026-11-12) / applied by l1-reviewer-package, PR #283; recorded by l1-nterm-res1-multipass, PR #280 |
| 2026-11-12 | Rebuilt on PR #283's second pass (6bf341f2): 561 members; 534 Igκ + tag, 24 tag alone, SLC61A1 ×2 still to move, cyno NIPAL4 to replace | #283 applied the MFSD2A, AQP12B and GLP2R moves and the OR8G3 ATG recode, and dropped the OPN3 isoform | no new decision / Claude Code session l1-nterm-res1-multipass, PR #280 |
| 2026-11-12 | Rebuilt on PR #283 at 13f3c1b8: cyno NIPAL4 (now 404 aa from Met215), CHRM3 and OR3A2 built on their restarted sequences; 535 Igκ + tag, 24 tag alone, SLC61A1 ×2 still to move | #283 applied the NIPAL4 replacement by restarting it at Met215, which this session had decided; the replace decision is tied to the old sequence hash, so the new form builds | no new decision / Claude Code session l1-nterm-res1-multipass, PR #280 |
| 2026-11-12 | Refold failures dispositioned: all 91 non-passing arms kept as built; 14 flagged for SpyTag003 occlusion (read ALFA), 17 for local fold disturbance | Failures cluster by family (taste and non-olfactory class A), not by Igκ vs tag alone; most "buried" calls are low-confidence | l1-nterm-res1-multipass (delegated by Becca 2026-11-12) / Claude Code session l1-nterm-res1-multipass, PR #280 |
| 2026-11-12 | All 559 arms scored (489 pass, 70 fail); fold-disturbed split into native-region disturbance (12, under review with spacer refolds and the ranker) and tail-only (5, kept) | Tail-only failures compare two disordered N-tail conformations; the native fold is intact | l1-nterm-res1-multipass (delegated by Becca 2026-11-12) / Claude Code session l1-nterm-res1-multipass, PR #280 |
| 2026-11-11 | Lane checked for every member of the rung (rerun on `main` after #283 landed, 2026-11-13: 560 hold, 1 unresolved) | Becca asked for every member's lane to be checked, not only the flagged ones | recommendation, not yet confirmed by Becca / Claude Code session l1-nterm-lane-check |
| 2026-11-13 | SLC61A1 (human, mouse) kept on this rung and built with Igκ + tag, superseding the move to internal | The writer withdrew its earlier move: v1, v2, SignalP 6 and the AlphaFold model all call the N-terminus out, and the lane check holds; this session's MFS-fold argument is the rejected alternative | l1-reviewer-package (delegated by Becca 2026-11-12) / Claude Code session l1-nterm-res1-multipass, PR #302 |
| 2026-11-13 | No internal-site move for the 12 native-region-disturbed arms; decision waits on the spacer comparison refolds | The ranker cannot vet internal sites without a membrane-side mask and features (none of the 12 is on the internal rung), and expects GPCR loop sites to fall to snorkels | l1-nterm-res1-multipass (delegated by Becca 2026-11-12) / Claude Code session l1-nterm-res1-multipass, PR #302 |
| 2026-11-13 | 24 expression-rescue arms built as Igκ + tag; ADGRD2 held for its native signal peptide | Rescue members were Igκ + receptor with the cassette deferred; built as Igκ + tag they pair with their tag-alone references by the leader only | l1-reviewer-package (delegated by Becca 2026-11-12); ADGRD2 hold by l1-nterm-res1-multipass / Claude Code session l1-nterm-res1-multipass |
| 2026-11-13 | ADGRD2 rescue arm built with Igκ replacing its native 18-residue signal peptide | An internal native signal sequence behind the cassette could be cleaved with the tag | l1-reviewer-package (delegated by Becca 2026-11-12), on this session's recommendation / Claude Code session l1-nterm-res1-multipass |
| 2026-11-13 | Each rescue arm's last sense codon recoded synonymously (superseded by the next row) | Cut the tract ending at the barcode to 24–26 nt; flagged by l1-cterm-build's insert audit, PR #308 | l1-nterm-res1-multipass (delegated by Becca 2026-11-12) / Claude Code session l1-nterm-res1-multipass |
| 2026-11-13 | Rescue receptors recoded until no 80-nt identical run remains against the reference or another rescue arm (11–37 codons per arm, longest tract 75–79 nt) | The writer's rule measures the longest tract anywhere up to the barcode, not only the run ending at it; RT template switching needs only a long shared run | l1-reviewer-package (delegated by Becca 2026-11-12) / Claude Code session l1-nterm-res1-multipass |
| 2026-11-13 | TSPO2 kept on N at residue 1 with Igκ + tag | 5 helices put exactly one terminus outside; DeepTMHMM v1 and TSPO1 homology say N-out, only v2 says N-in, and the Igκ arm passes the refold. Low confidence; a glycosylation reporter or an N- vs C-tagged staining pair would settle it. Rejected: the C-terminal lane | l1-nterm-res1-multipass (delegated by Becca 2026-11-12) / Claude Code session l1-nterm-res1-multipass |
| 2026-11-13 | Rescue receptors restored to their designed DNA; each rescue arm instead gets its own linker-only cassette encoding (tags on p003's codons) | Becca: "sharing bases if that's the prot seq is fine I jsut want to diversify the linker since these proteins are at low abundance in the pool". The within-pair barcode-swap prediction is unchanged, because the linker is upstream of the shared receptor. Rejected: the receptor recode of the previous row | Becca, 2026-11-13 / Claude Code session l1-nterm-res1-multipass |
| 2026-11-13 | One arm per gene for the 25 expression-rescue genes: each keeps its library reference with its lane's default design (Igκ + tag for the 23 on this rung; ADGRD2 and GLP2R as on the N-after-SP lane), and the 25 rescue members are retired | Becca: "we are pooling the whole L1 library together. why resuce pairs no i think we should make a single choice for each L1". For the 23 the reference and its twin are the same protein; GLP2R's rescue arm would put the tag in front of a predicted cleaved signal peptide; ADGRD2's native signal peptide has nothing against it. Rejected: keeping both arms (supersedes the rows above on rescue pairs, their tag-alone references, recodes and linker encodings) | Becca, 2026-11-13 (the rule); per-gene choice by l1-nterm-res1-multipass, delegated by Becca 2026-11-12 / Claude Code session l1-nterm-res1-multipass |
| 2026-11-13 | Member ids re-keyed to the renumbering on `main` (#313); the prose cites design ids, which do not change | Ids only: every DNA and protein sequence is byte-identical. Mouse SLC61A1 (`l1-p3105`) now has a refold verdict (fail, low-confidence burial; kept as built) | no decision taken / Claude Code session l1-nterm-res1-multipass |
| 2026-11-14 | All 72 refold failures kept as built and flagged, each with its pattern, metric, rejected alternative and the bench result that would change it; buried SpyTag003 at the mature N-terminus (53 arms) left as an open design question for Becca | Spacer comparison closed (a 10-residue spacer rescued 0 of 12; 30 made snorkels worse); no internal site can be vetted; ALFA stays exposed in all 53 buried-SpyTag003 arms | l1-nterm-res1-multipass (delegated by Becca 2026-11-12) / Claude Code session l1-nterm-res1-multipass |
| 2026-11-14 | OR11H6 (`l1-p1381`) and TMEM63B (`l1-p2199`) routed in from "N after SP" and built Igκ + tag at Met1 on the full native sequence | Only DeepTMHMM v1 calls a signal peptide; SignalP 6, DeepTMHMM v2 and UniProt (for TMEM63B, the cryo-EM structure) put residues 1–43 and 1–40 outside ahead of the first helix, so the Igκ swap deleted native extracellular residues | l1-nterm-res1-multipass (delegated by Becca 2026-11-12), on the N-terminal signal-peptide audit's proposal / Claude Code session l1-nterm-res1-multipass |
| 2026-11-14 | SLC41A2 (`l1-p1957`) routed in from "C-terminal" and built Igκ + tag at Met1 | Experimental topology (Sahni 2007; UniProt ECO:0000269) puts the C-terminus inside and the 162-residue N-terminus outside, the reverse of both DeepTMHMM versions; the tail is too long for TM1 to translocate untargeted | l1-nterm-res1-multipass (delegated by Becca 2026-11-12), on the coordinator's rescue-round proposal / Claude Code session l1-nterm-res1-multipass |
| 2026-11-15 | Human OR10G6 (`l1-p1362`) routed in from "N after SP" and built Igκ + tag at Met1; mouse `l1-p3043` therefore stays built | Tedman's leader-free HA-after-Met1 construct displayed at 1,029 a.u., so the N-terminus is not cleaved; checked because Becca's ortholog rule would otherwise hold the mouse ortholog | l1-nterm-res1-multipass (delegated by Becca 2026-11-12), on the coordinator's check / Claude Code session l1-nterm-res1-multipass |
