# L1/L2 ortholog work: agent handoff

Handoff state verified November 7, 2026, after PR 262. The completed code and review records are merged. No prediction jobs from this work are running. The remaining task is scientific review and downstream adoption of the proposed panel, not finishing the prediction batches.

## Start here

- Repository: https://github.com/Deliverome-Project/deliverome-analysis
- Experiment: `analysis/sprint4/s4e9-orf-design-fg-panel/`
- Read repository `AGENTS.md`, experiment `README.md`, and `ORTHOLOG_CDS_POLICY.md` before making changes.
- Final run index: `ortholog_cds_policy_final_run.json` in the experiment directory. It records results, database verification, snapshot identity, and durable archive locations.
- Start from current `origin/main` in your own suitable checkout. The completed analysis branch was `becca/ortholog-cds-eligibility`; its work is merged.
- This document is self-contained for a Claude agent: Codex chat history and local workspace access are not required. Retrieve large artifacts from the immutable S3 locations below.

## Merged work

| PR | Work | State |
|---|---|---|
| [256](https://github.com/Deliverome-Project/deliverome-analysis/pull/256) | Correct SpyTag003/ALFA tandem labels on loop constructs | Merged |
| [257](https://github.com/Deliverome-Project/deliverome-analysis/pull/257) | Account for all 451 isoform cohort exclusions | Merged |
| [258](https://github.com/Deliverome-Project/deliverome-analysis/pull/258) | Reproducible historical ortholog accession selection | Merged |
| [260](https://github.com/Deliverome-Project/deliverome-analysis/pull/260) | CDS boundary/completeness audit and evidence review | Merged; audit only |
| [262](https://github.com/Deliverome-Project/deliverome-analysis/pull/262) | Systematic CDS eligibility, APPRIS evidence, full redraw, prediction persistence | Merged; review panel only |

PR 262 merge commit: `713e8b6f2ee1e1767babe3234030a939adb8d635`.
Final analysis producing commit: `5f6f8870248538d24451b6230640e2becf06e62c`.
Final review-data commit: `0c2adc4cbc415e57d274c4c9e7a5586fd23fa5ad`.
120 targeted tests and all required CI checks passed. Independent agents reviewed policy semantics, exact final sequences, quotas, prediction export, and final documentation.

## Final policy and result

Use frozen Compara 116 relationships across the full library gene pool. Audit all transcript models before protein ranking. Automatic eligibility requires native ATG, supported coding-boundary coordinates, coherent source identities/versions, in-frame unambiguous CDS, terminal stop, no internal stop, and exact translated-peptide agreement. Coding NF tags block automatic selection; mRNA NF is advisory only when independent CDS checks pass. No implicit start/stop repairs are made. The exception registry is empty.

Mouse APPRIS uses verified release `2026_08.v50`, dataset `e116v50`. Ensembl 116 GTFs themselves contain no APPRIS tags for either species. Cyno has no applicable APPRIS dataset; rhesus assignments are not transferred. Missing APPRIS gives no ranking preference and is not an exclusion. Equivalent encodings rank annotation evidence; distinct eligible proteins retain the documented human-protein correspondence metric and tie rules.

Final counts: 25,204 transcript models; 3,678 scored gene-pair representatives; 3,643 sampleable representatives; 59 pairs without eligible CDS; one missing human reference (OR4Q2); 31 known non-surface and four ambiguous-anchor exclusions. No pending winner predictions.

The final panel fills all 200 slots (100 mouse, 100 cyno) and all 20 fixed quotas with the frozen correspondence cutoffs and seeds. All 200 CDSs match pinned native Ensembl source bytes and translate to the exact selected protein. There are no duplicate selected proteins or quota shortfalls.

The full redraw retains 14 old entries with identical proteins and replaces 186. Of those removed, 183 remain CDS-eligible and remain the exact full-pool representative: their removal is resampling, not a defect or sequence repair. The other three old forms are cyno OPN3, ADCY2, and PRRG1, which lack supported initiation.

## Predictions and preservation

- Completed 58 DeepTMHMM v1 predictions (39,859 residues), in two Modal batches of 57 and one. Same image definition and mounted model volume; package/runner/weight hashes recorded. The observed package reports 1.0 Academic, so new rows use `deeptmhmm-1.0-academic`. Do not claim the older 1.0.24 label was verified, or rewrite historical records to infer it.
- Completed 11 NetGPI 1.1 predictions; job `6AED3BD30003C19C5A3B159A`. Exact input sequences, run IDs, prediction/input ordinals, literal likelihoods, raw responses and hashes are retained.
- All results were appended to BOTH D1 stores and every imported field verified. Each store preserves the prior 2,359 NetGPI predictions, 2,187 independent GPI curation records, and prior 323-row topology release. NetGPI totals are now 2,370 per store.
- New topology release: `topo_2026_11_05_ortholog_cds_policy`; cohorts `mouse_ortholog_cds_policy` and `cyno_ortholog_cds_policy`.
- New immutable snapshot: `external://accessible-surfaceome/topology/2026-11-06-ortholog-cds-policy/topology.parquet`, SHA-256 `7683d75a3e04228d95fcfb8fe9b7e9b46b38c0ea254ca42550a666f5d8bbf216`, 71,697 rows. All 71,639 previous snapshot rows preserved. The default snapshot constant is unchanged.
- Post-import served-record audit: 65 affected records, zero stale features or integrity gaps. History sweep: all 65 unchanged. No record refresh was needed.
- Existing public D1 and the frozen snapshot differ by three historical rows; do not equate their totals. D1 topology totals after append: private 34,953, public 71,694.

## Durable records and local artifacts

| Record | Location |
|---|---|
| Complete clean replay, input references, commands, code archive, full ledgers | `processed://sprint4/s4e9-orf-design-fg-panel/library-design/ortholog-cds-eligibility/2026-11-05/final-v2/` |
| Predictions, raw responses, runtime identities and source-bound provenance | `processed://sprint4/s4e9-orf-design-fg-panel/ortholog-cds-policy-predictions-2026-11-05/`; exported D1 payload/index under `bundle/` |
| Before/after D1 verification, independent review code/report, audit/history logs, merge receipt and final CI | `processed://sprint4/s4e9-orf-design-fg-panel/library-design/ortholog-cds-eligibility/2026-11-05/final-v2-verification/` |

The final committed review files in this experiment directory are `ORTHOLOG_CDS_POLICY.md`, `ortholog_cds_policy_final_run.json`, `ortholog_panel_cds_eligible.tsv`, `ortholog_cds_panel_changes.tsv`, `ortholog_cds_final_before_after.tsv`, and `ortholog_cds_final_independent_review.json`. Exact selected CDSs are archived at `final-v2/results/panel/selected_cds.json.gz`; use the archive manifest to verify the stored and decompressed hashes. The archive includes producing source code and the replay command journal.

Earlier summaries and missing-prediction FASTA worklists remain historical inputs; they do NOT indicate outstanding jobs. Prefer files named final/v2 and the final run index.

D1 importer/audit tooling came from the clean accessible-surfaceome checkout at commit `484a19bc3da7f46a0d8e8f3e36ba991958729923`; exact script/lock hashes are in `final-v2-verification/external_toolchain.json`. Use a suitable checkout of that repository if further database work is needed. Never print credentials. Both D1 stores are shared production sources, not a main-versus-development split.

## Remaining work

1. Review biological correspondence before activating the proposed panel. CDS completeness does not prove preserved domains or functional equivalence. Explicit examples: cyno SLC8A1 624 aa versus the actual human library protein's 973 aa; IL31RA 598/764; MS4A10 226/267. Differences include internal unmatched sequence. Do not impose an arbitrary new cutoff or call these missing domains without evidence; document any new policy and rerun reproducibly. **Reviewed November 7, 2026:** [ORTHOLOG_CORRESPONDENCE_REVIEW.md](ORTHOLOG_CORRESPONDENCE_REVIEW.md) (code `review_ortholog_correspondence.py`, table `ortholog_correspondence_review.tsv`) finds 183 forms keep the human architecture, 11 should be replaced (all three named cases among them) and 6 need a decision; each flagged form has a reviewed next-in-draw replacement. Adopting those changes remains part of item 2.
2. Decide adoption scope with the user, then implement the exact selected CDSs through member-ID allocation, L1/L2 tranche updates, construct regeneration, and downstream audits. No canonical L1/L2 manifest or construct was activated by PR 262. Do not silently reuse an old CDS merely because it encodes the same protein.
3. Keep the earlier 61 diagnostic CDS initiation changes held. They were audit proposals, not approved production sequence repairs; this work does not authorize publishing engineered replacements or their predictions.
4. Recheck outstanding GitHub issues after adoption. Merging the audit/review PR is not proof that the original ortholog-adoption issue is closed.

## User authorizations and constraints

The user requested reproducible audits/filters/fixes, preservation of all NetGPI predictions and ad hoc annotations, D1 persistence, PRs, subagent work, and merging confident changes (including admin merge). Those completed actions need no repeat approval.

For this prediction revision, the user approved up to 75 paid Modal predictions: 58 completed, 17 unused. This is a bounded revision-specific approval, not blanket funding for another project or engineered constructs. Native public replacement sequences have approval for NetGPI checks through this rebuild. No jobs remain necessary for the completed review panel. Ask only for genuinely new scope, spending beyond the cap, or an unresolved scientific choice.

Preserve immutable S3 versions and existing database history. Keep exact sequence provenance separate from accession aliases and tool version separate from library tranche. Repository instructions require committed producing code, normal checks/hooks, and merge commits; do not bypass failing tests.
