# Ortholog CDS eligibility and representative selection

Ortholog representative selection checks coding-sequence support before comparing
proteins with the intended human library form. The candidate pool includes all
high-confidence one-to-one mouse and cynomolgus Compara 116 relationships for the
library genes. Current L1/L2 records remain preserved; this revision produces
review candidates, not an activated manifest or ordered constructs.

## Eligibility and ranking

| Stage | Rule | Reason and retained evidence |
|---|---|---|
| Gene relationship | Use the frozen Compara 116 pool and stable human/target gene IDs | Preserve the full pool, including relationships without existing members and explicit ambiguous or unavailable mappings. |
| Transcript scope | Inspect the union of CDS, peptide and GTF models in each eligible target gene | Missing source records remain visible. Counts include noncoding GTF models and therefore are not counts of candidate proteins. |
| CDS identity | Require unambiguous, in-frame DNA with a terminal stop and exact translated-peptide agreement | Bind species, gene, transcript/protein versions, assembly, strand and exact sequence hashes. |
| Boundary evidence | Require a native ATG start, supported spliced start/stop coordinates, phase zero and no unresolved coding-boundary tags | An ATG alone does not establish a supported start. CDS NF means a coding boundary is unconfirmed, not proven truncation; mRNA NF remains advisory when independent coding checks pass. |
| Transcript class | Require protein-coding annotation for automatic selection | Decay-associated, loss-of-function and other unresolved models remain review cases; no biological inactivity is inferred. |
| Documented exceptions | Require an explicit checksum-bound exception registry | No implicit ATG/stop additions or blanket acceptance of uncertain models. The initial registry is empty. |
| Equivalent CDS encodings | Rank verified APPRIS evidence, Ensembl canonical status and transcript support level, then transcript ID | These candidates encode the same protein; annotation support replaces automatic reuse of the old CDS or lexical transcript choice. |
| Different eligible proteins | Maximize matching sequence blocks divided by the longer protein length against the actual human library protein | Retain the established correspondence metric. Within 0.001 of the maximum, prefer reviewed evidence, verified APPRIS principal status, remaining transcript-quality rank, current accession, then lexical accession. |
| Surface evidence | Resolve DeepTMHMM v1 and, for GPI genes, NetGPI 1.1 by exact protein sequence | A missing prediction is pending work, never a negative call or a reason to favor an already-predicted lower-ranking form. |

Boundary-supported eligibility does not prove biological full length, expression,
or identical domain architecture. Correspondence and topology evidence remain
available for scientific review; no new arbitrary protein-length cutoff is added.
A missing human reference, a gene with no eligible CDS, or a known ambiguous anchor
result remains an explicit sampling exclusion. Missing predictions for a provisional
winner block the full draw until the result permits rescoring.

## APPRIS availability and identity

The pinned Ensembl 116 GTF files contain no APPRIS tags for either species. Mouse
ranking instead uses the version-matched APPRIS `e116v50` assignments from release
`2026_08.v50`, verified against the corresponding CDS/GTF identity and peptide.
Stable transcript ID agreement alone is insufficient. The [official APPRIS
configuration](https://apprisws.bioinfo.cnio.es/config.json) lists mouse GRCm39 with
Ensembl 116; its macaque dataset is rhesus (*Macaca mulatta*, Mmul_10), not
cynomolgus (*Macaca fascicularis*). Rhesus annotations are not transferred to cyno.
The observed configuration is preserved as an immutable source snapshot.

| Evidence | Mouse | Cyno |
|---|---:|---:|
| Target genes | 1,876 | 1,862 |
| Candidate transcripts with any external APPRIS assignment | 6,657 | Unavailable |
| Verified principal transcripts | 4,956 | Unavailable |
| Verified principal transcripts passing CDS eligibility | 4,955 | Unavailable |
| Genes with an eligible verified principal transcript | 1,869 | Unavailable |

Six principal-ID matches fail sequence/source validation and contribute no APPRIS
preference. One verified principal transcript fails the CDS gate; its unconfirmed coding end
and missing stop-feature evidence remain recorded. Multiple transcripts can share a principal protein, so transcript
counts and gene/protein counts are separate. Missing APPRIS contributes no ranking
preference and is not an exclusion. Cyno therefore uses eligible CDS evidence,
protein correspondence, available reviewed evidence, canonical/support annotations,
and deterministic ties without an APPRIS contribution.

## Full-pool review and fixed panel design

The clean replay covers 20,730 mouse and 4,474 cyno transcript models. It retains
12,428 and 4,218 eligible transcripts, representing 8,054 and 4,181 gene/protein
groups. Five mouse and 54 cyno genes lack an eligible CDS.

Of 3,738 eligible gene relationships, 3,678 have a scored representative, 59 lack an
eligible CDS, and one lacks the human designed reference (OR4Q2). Relative to the
prior full-pool scoring, 3,491 representatives retain their protein and 187 change.
There are 3,643 sampleable representatives: 1,864 mouse and 1,779 cyno. Thirty-one
known non-surface results and four ambiguous anchor results remain explicit
exclusions. There are no pending winner predictions.

The 200-member review panel fills every quota without shortfalls. It retains 14
old gene/species entries with identical proteins and replaces 186 through the full
redraw. Of the removed entries, 183 remain CDS-eligible and remain the exact
full-pool representative; resampling explains their removal. The other three old
forms—cyno OPN3, ADCY2 and PRRG1—lack supported initiation. These panel changes
therefore do not count defective sequences or sequence repairs.

Independent checks confirm all 200 CDSs match native Ensembl 116 source bytes,
translate exactly to the selected protein hashes, and have supported coding
boundaries. All correspondence scores reproduce independently. No selected form
has an mRNA-boundary advisory. The largest length differences include cyno SLC8A1
(624 versus 973 human residues), IL31RA (598 versus 764), and MS4A10 (226 versus
267). These differences include internal unmatched sequence; domain equivalence
has not been established and remains part of pre-activation scientific review.

The draw retains the published correspondence cutoffs (mouse 0.8644; cyno 0.9595),
seed 20261011 for mouse and 20261012 for cyno, and the following quotas. Near means
at or above the cutoff using the recorded four-decimal score. Protein class comes
from the exact human reference topology, not the ortholog's predicted class.

| Species/band | GPI | Multipass | Type I | Type II | Type III | Total |
|---|---:|---:|---:|---:|---:|---:|
| Mouse near | 7 | 23 | 12 | 6 | 2 | 50 |
| Mouse diverged | 7 | 18 | 15 | 8 | 2 | 50 |
| Cyno near | 7 | 22 | 12 | 7 | 2 | 50 |
| Cyno diverged | 7 | 18 | 15 | 7 | 3 | 50 |

Candidates sort by stable human/target IDs and protein hash before each seeded
shuffle. This is a documented full redraw, not a claim to preserve the earlier
iteration order. Exact protein duplicates are excluded against fixed active L1
non-ortholog members, including its 25 constructed rescue arms, and across the
selected ortholog panel. Existing L2 ortholog reuse remains allowed. Shortfalls
block the panel without borrowing quotas or relaxing eligibility. A selected CDS
retains its transcript and DNA hash even when the protein matches an old member;
protein identity does not authorize silently substituting the old DNA.

## Reproduction and provenance

| Record | Location and purpose |
|---|---|
| Producing code | [Source revision 5f6f8870](https://github.com/Deliverome-Project/deliverome-analysis/tree/5f6f8870248538d24451b6230640e2becf06e62c); eligibility, scoring and draw write rule versions, commands and input/output hashes. |
| Mouse APPRIS input | [Official versioned principal assignments](https://apprisws.bioinfo.cnio.es/pub/releases/2026_08.v50/datafiles/mus_musculus/e116v50/appris_data.principal.txt); SHA-256 `1e9d500c61149d9a9fbe2af6b5657fa77bc72201d610357be435b59d2c0abe04`. |
| Pull request | [Ortholog CDS policy review](https://github.com/Deliverome-Project/deliverome-analysis/pull/262) |
| Reference source manifests | `ortholog_cds_eligibility_v2_inputs.json` and `ortholog_rebuild_inputs.json`; exact release-specific source keys and checksums. |
| Durable review | `processed://sprint4/s4e9-orf-design-fg-panel/library-design/ortholog-cds-eligibility/2026-11-05/final-v2/`; full transcript/candidate/pair ledgers, representative CDSs, APPRIS coverage, before/after decisions and clean run journal. |
| Prediction provenance | Exact submitted sequence, actual tool/model identity, run/job identifiers, raw output, NetGPI prediction ordinal and literal likelihood. Existing predictions and ad hoc curation remain separate records. |

The policy has 120 targeted tests and independent review covering exact
identity, APPRIS evidence binding, CDS integrity, pending predictions, duplicate
proteins and fixed quotas. Model runtime provenance records the observed academic
DeepTMHMM 1.0 identity and package hashes; a patch version is not inferred from a
previous export label. Historical more-specific labels remain preserved but are
not independently verified by the currently available package identity.

## Prediction persistence

The two Modal batches contain 58 successfully predicted proteins (39,859 residues),
within the approved limit of 75. Both use the same image definition and mounted
academic model volume. NetGPI 1.1 supplies 11 new exact-sequence predictions;
its job ID is `6AED3BD30003C19C5A3B159A`. Each record retains the exact submitted
protein, input/output ordinal, literal likelihood, raw output, and source run.

Both D1 stores contain the full append. Every imported topology and NetGPI field
matches the exported payload. All 2,359 prior NetGPI predictions, 2,187 separate
GPI curation records, and the previous 323-row topology release remain unchanged
in each store. This preserves the distinction between model output and ad hoc
curation. The new records use `topo_2026_11_05_ortholog_cds_policy` and the observed
`deeptmhmm-1.0-academic` label; the provenance index binds the runtime package and
weight hashes to each prediction.

| Durable record | Location |
|---|---|
| Prediction export and provenance index | `processed://sprint4/s4e9-orf-design-fg-panel/ortholog-cds-policy-predictions-2026-11-05/bundle/` |
| Exact prediction inputs and raw outputs | Source-specific URI and SHA-256 in `prediction_provenance_index.json` |
| Topology snapshot | `external://accessible-surfaceome/topology/2026-11-06-ortholog-cds-policy/topology.parquet`; 71,697 rows, preserving all 71,639 prior snapshot rows |
| Database verification and independent panel review | `processed://sprint4/s4e9-orf-design-fg-panel/library-design/ortholog-cds-eligibility/2026-11-05/final-v2-verification/` |

The post-import audit checks 65 affected served records and finds no stale features
or source-integrity gaps. The record-history sweep reports all 65 unchanged.

The default topology snapshot and active L1/L2 manifests remain unchanged. The
new snapshot is an explicit reusable input; adding its prediction cohorts to
species lookup permits exact-sequence reuse without substituting a different
species or protein.

## Open decisions

| Question | What resolves it |
|---|---|
| Which review revision becomes active? | Review the final before/after panel and explicit CDS choices before downstream construct regeneration and activation. |

## Decision log

| Date | Previous approach | Current approach | Reason |
|---|---|---|---|
| 2026-11-05 | Choose protein, then reuse or resolve a matching CDS | Gate all transcripts before protein ranking and retain an evidence-ranked encoding | Boundary support must affect the candidate pool before sampling. |
| 2026-11-05 | Mouse APPRIS evidence comes from a release-115 table | Use verified native release-116 assignments; cyno remains explicitly unavailable | Match the annotation release and avoid stable-ID-only or cross-species evidence transfer. |
| 2026-11-05 | All NF tags block the initial automatic gate | Coding NF remains blocking; mRNA NF is advisory under eligibility v2 | Unconfirmed untranslated transcript boundaries do not invalidate independently supported CDS boundaries; eight mouse transcripts are affected. |
| 2026-11-06 | Provisional pool awaits topology and anchor evidence | Full pool rescored with 58 topology and 11 NetGPI predictions; all 200 review slots filled | Missing predictions must not favor previously predicted candidates; exact native CDSs and fixed quotas pass independent review. |
