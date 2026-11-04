# S1E3 partner audit — 2026-09-10

The reviewed amendment adds **ATP1B3/CD298 (T1), CDC50A/TMEM30A (T2), HFE (T1),
and P2RX6 (T2)**. The frozen 5,130-gene input reproduces the published 2,340-gene
membership exactly before the amendment. The revised membership is **2,344**, with
no removals or changes to existing tiers. These are computational inclusion decisions;
median RNA abundance is a screening proxy for endogenous partner availability, not
a demonstration of spare assembly capacity in every recipient cell line.

## Root cause and narrow rule

The old code has two independent problems. `coreceptor_ok` requires every
expression-known partner in a flat list to exceed median HPA cell-line TPM 5, even
when the list describes alternative subunits. Separately, `stable_complex_partner`
is always disqualifying, including for integral membrane subunits with their own
extracellular domain. Fixing only AND versus OR would leave CD298 excluded.

`coreceptor_rules.py` records six reviewed cases by the cohort's UniProt accession.
A requirement is an AND of roles, each satisfied by at least one explicitly supported
alternative. Only literature-supported sufficient routes are encoded: this does not
assert that every reported interaction partner can support surface trafficking.
Partner expression uses the saved resolver-backed Ensembl IDs and requires exactly
one matching HPA matrix gene; it never falls back to a display-symbol join. Missing,
non-finite, or ambiguous expression cannot satisfy a reviewed requirement. Changes
to the observed partner set fail closed until re-reviewed.

Only four reviewed integral membrane complexes may clear the stable-complex reason
veto. Their partner requirement must be satisfied and their canonical topology must
still contain a TM helix. Membership in canonical/likely, surface specificity,
accessibility, expression, ECD size, and ORF length remain necessary. No blanket
replacement of `all` with `any`, and no blanket acceptance of stable complexes.

## Evidence and affected proteins

| Protein | Sufficient endogenous route (median HPA cell-line TPM) | Old exclusion | Revised tier | Supporting evidence |
|---|---|---|---|---|
| ATP1B3 / CD298 | ATP1A1, 216.3 | Partner + stable-complex reason | T1 | Saved cohort: canonical, surface-dominant, single-pass type II; alpha/beta heterodimer. CD298 surface staining and alpha association: [Chiampanichayakul et al.](https://pubmed.ncbi.nlm.nih.gov/17176442/). |
| CDC50A / TMEM30A | ATP8B1, 5.55 **or** ATP8B2, 12.1 | Partner + stable-complex reason | T2 | P4-ATPase complexes traffic to PM; not every ATPase combination has the same destination. [van der Velden et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC3000991/) and [Bryde et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC3003355/). |
| P2RX6 | P2RX4, 9.7 | Partner + stable-complex reason | T2 | Functional P2X4/6 heteromers: [Lê et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC6793241/); human endothelial PM localization: [Glass et al.](https://pubmed.ncbi.nlm.nih.gov/12088286/). P2RX2 (0.1) is not a simultaneous requirement. |
| HFE | B2M, 1054.45 | Stable-complex reason only | T1 | Wild-type HFE associates with B2M and reaches the surface in transfected COS-7 cells: [Waheed et al.](https://pubmed.ncbi.nlm.nih.gov/9356458/). This does not rescue the trafficking-deficient C282Y variant. |

The CDS/isoform/ortholog selection policy has not changed. Additional topology-competent
isoforms are ATP1B3 1, CDC50A 2, HFE 8, P2RX6 1; orthologs are 1, 2, 2, 2 respectively.
ATP1B3's mouse ortholog remains below the existing 70% ECD-identity threshold.

| Count | Published | Revised | Change |
|---|---:|---:|---:|
| T1 genes | 964 | 966 | +2 |
| T2 genes | 878 | 880 | +2 |
| T3 genes | 286 | 286 | 0 |
| LL / SV / OE genes | 115 / 85 / 12 | 115 / 85 / 12 | 0 |
| Total genes | 2340 | 2344 | +4 |
| Surface-competent isoforms | 2333 | 2345 | +12 |
| Surface-competent orthologs | 3029 | 3036 | +7 |
| Oversized genes, outside library | 99 | 99 | 0 |

None of the four additions occurs in the published reference-coverage table. SurfaceBench
62/68, ViralZone 54/62, ADC 205/234, and TCE 25/30 therefore remain unchanged.

## Complete audit and deferred candidates

The audit enumerates **343 unique genes**: all **336 required-partner** records plus
seven additional stable-complex records. Older documentation's 333 required genes
describes an earlier cohort. `coreceptor_audit.tsv` gives every audited gene's accession,
partners/median TPM, all failed gates before and after, and tier outcome.

- **ATP1B1/B2, ATP1A1/A2/A3 and ATP12A:** already admitted by saved `INCLUDE_ALONE`
  exceptions, whose own rationale actually invokes common endogenous partners. They
  are affected by the same conceptual issue, but are not missing from the 2,340 list.
  This patch leaves their historical exceptions intact rather than recounting them as additions.
- **ATP8B1 and CDC50B/TMEM30B:** reviewed alternative routes now pass the partner gate,
  but both remain excluded by `mostly_intracellular`. ATP8B1's saved rationale conflates
  cytoplasmic protein-domain mass with subcellular location; the primary studies show
  a PM complex. That is a plausible upstream localization-annotation error for a separate
  review, not grounds to bypass this independent gate here.
- **NCSTN:** high-priority follow-up. PSEN1/2 are alternatives, but APH1 and PSENEN
  represent additional roles. The saved `APH1` token lacks an expression match; the
  old code ignores unknown partners. Establish the full assembly route before rescue.
- **SGCE/SGCZ:** the earlier exception passes them, but their documented partner lists
  contain restricted sarcoglycans. Reassess independent trafficking versus assembly
  requirements before removing their stable-complex exclusion.
- **KCNK15:** a common YWHAB chaperone does not establish that rare TASK channel
  partners are dispensable. Its historical exception is insufficient for a new rescue.
- **B2M:** still lacks its own conventional membrane anchor. Its association with
  HLA-E cannot satisfy this library's integral/GPI-anchored target design.
- **HLA-DQA1/DQB2/DRB4, SGCA/SGCD/SGCG/SSPN, TCR/CD3 and BCR components:**
  a common chaperone or one common subunit does not satisfy the full complex.
- **GRIN2A-D/GRIN3A-B:** remain excluded; their obligate GRIN1 partner has median 0.
- **RAMPs, neuronal Kv auxiliary subunits, ABCG5/8, SLC51B, SLC7A9:** the saved
  alternative lists offer no qualifying common partner. Alias/phrase ambiguities
  (for example GluA, GluK, APH1) remain review items, not positive evidence.

The exhaustive gate audit is not a claim that all 343 underlying literature annotations
were independently revalidated. Six routes were explicitly reviewed in this amendment.
Legacy unknown-partner skipping and older exceptions remain visible in the audit;
changing them globally could remove proteins and is outside this narrow correction.

## Reproduction and validation

Inputs and audit outputs are frozen in the approved team archive at
`processed://sprint1/s1e3-surface-library-filters/partner-audit-2026-09-10/`
and the accompanying `s1e3-partner-audit-reproducibility.zip` deliverable.
`inputs/` contains the cohort, resolver joins,
median HPA table, saved rescue lists, and published membership.
`coreceptor_audit_summary.json` records SHA-256 hashes. These are analysis snapshots;
source expression data remains at `external://expression/2026-08-04/`.
The canonical published membership is not overwritten by opening this PR.

Copy the bundle's `inputs/<relative-path>` files into this experiment directory, then run:

```sh
uv run pytest tests/test_s1e3_coreceptor_rules.py
uv run python analysis/sprint1/s1e3-surface-library-filters/audit_coreceptors.py
uv run python analysis/sprint1/s1e3-surface-library-filters/tier_plot.py
```

`audit_coreceptors.py --hydrate` downloads those same archived inputs and verifies
their hashes before running. Access requires the project's processed-bucket credentials.

The audit runs the pre-change pipeline from commit
`6b3d53ccab9b95246a0605b79427d69dca36aeb2` against the same inputs, asserts an exact
gene-and-tier match to the published baseline, then runs the amended pipeline. It checks
the exact four additions, zero removals/retiering, and 28 negative integration cases
covering other gates on all four added records. Unit tests cover alternatives, full
complexes, unknown/low/non-finite expression, changed partner lists, and missing anchors.

Notion experiment: https://www.notion.so/3b8ab86c768381bc8a73dc79ec72d832
