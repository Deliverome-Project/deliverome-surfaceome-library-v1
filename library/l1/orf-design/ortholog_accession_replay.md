# Historical ortholog accession selection

Reapplying the original gene-symbol-to-UniProt accession rule against UniProt 2026_04 confirms 1,492 of 3,892 recorded successful picks. The other 2,400 successful picks remain unconfirmed by this release. This is an audit of the original mouse and cynomolgus cohorts; it makes no changes to the library or construct records.

| Cohort | Recorded successful picks | Same accession wins | Different accession wins | No exact gene-name candidate | Original unsuccessful records |
|---|---:|---:|---:|---:|---:|
| Mouse | 1,971 | 1,387 | 584 | 0 | 1 |
| Cynomolgus | 1,921 | 105 | 88 | 1,728 | 6 |
| Total | 3,892 | 1,492 | 672 | 1,728 | 7 |

The seven original unsuccessful records are retained in the review table but do not enter the confirmation denominator. Every confirmed accession also retains its recorded length; the audit does not compare amino-acid sequences.

## Selection and evidence

The pinned historical source ranks candidates by exact primary-symbol match, exact match to any gene-name token, reviewed status, length, then accession, in descending order. Reviewed status precedes length. A description of this rule as simply selecting the longest symbol-matched entry is inaccurate.

The offline search uses case-insensitive exact tokens from UniProt's complete gene-name field within the original taxon: mouse 10090 or cynomolgus 9541. The pinned taxonomy streams contain 261,799 and 62,196 entries respectively. Independent search response counts match these totals. Five representative exact searches reproduce the offline candidate accession sets: Ngfr and Nalcn in mouse, and PSEN1, CXCR6 and SLITRK1 in cynomolgus. These checks support the token interpretation; they are not a recovered historical search response for every gene.

The largest current candidate set contains 136 mouse entries or 11 cynomolgus entries. Recorded original counts reach 41 and 16 respectively. Neither cohort reaches the historical unpaginated 500-result search limit, so no row in this audit is unresolved because of that limit.

The large cynomolgus no-candidate group reflects the gene-name fields available in the pinned release. For 1,726 of those 1,728 successful historical picks, the recorded accession still exists but does not qualify for the queried symbol in this snapshot. For example, G7PVK2 exists at taxon 9541 with length 663 and empty gene-name fields; an independent exact SLITRK1 search returns zero results. The source evidence does not establish when or why the searchable annotation differed from the historical selection, and this audit does not infer a species-wide annotation change between particular releases.

## Permanent limits

The original live UniProt search responses and release headers were not saved. A later pinned release cannot recover their complete candidate sets, search-index state, or transient failures. A same-accession result confirms agreement with the pinned release, not the provenance of an unsaved response. A different-accession or no-candidate result does not demonstrate that the historical choice was wrong. Gene-symbol matching also does not independently establish the Ensembl gene-to-protein relationship; the original target and human gene identifiers are preserved without remapping.

## Materials and provenance

| Element | Record |
|---|---|
| Historical producing source | [Pinned original chooser](https://github.com/Deliverome-Project/accessible-surfaceome/blob/f9c73691889c1ddb50f63f482ed0c42c9cf7b40a/src/surface_proteome/candidates/download_deeptmhmm_sequence_sets.py) |
| Audit source and input manifest | [Committed replay script](https://github.com/Deliverome-Project/deliverome-analysis/blob/1dbbb68/library/l1/orf-design/replay_ortholog_accessions.py), `ortholog_accession_replay_inputs.json`; each raw object has an immutable key and SHA-256 digest |
| UniProt release | REST response headers report 2026_04, dated November 4, 2026; downloaded November 5. The FTP `current_release` release-date file still reports 2026_03, so it is not used as the release identity |
| Raw taxonomy and validation responses | `external://uniprot/2026_04/ortholog-accession-replay-2026-11-05/`; source URLs and headers are pinned with the input manifest |
| Original cohort metadata | Existing immutable `external://accessible-surfaceome/f9c73691/data/external/deeptmhmm_surfaceome_inputs/` objects, verified against their manifest digests |
| Derived audit | `processed://sprint4/s4e9-orf-design-fg-panel/ortholog_accession_replay.tsv`, `ortholog_accession_replay_summary.json`, and `ortholog_accession_replay_run.json`; duplicate immutable review bundle under `ortholog-accession-replay/2026-11-05/` |
| Scope request | [L1 library review issue 19](https://github.com/Deliverome-Project/deliverome-surfaceome-library-v1/issues/19) |
| Review record | [Accession audit pull request](https://github.com/Deliverome-Project/deliverome-analysis/pull/258), dependent on [ortholog rebuild review](https://github.com/Deliverome-Project/deliverome-analysis/pull/234) |

The run record gives the exact producing commit, input-manifest digest, output digests and rerun command. `--download-inputs` retrieves missing inputs through the existing S3 helper; all files are checksum-checked and every pinned REST response must declare the same release. Six offline tests cover the ranking, exact aliases, mismatches, result cap, taxonomy and duplicate checks, and length changes.

## Decision log

| Date | Decision | Reason |
|---|---|---|
| 2026-11-05 | Preserve the historical selection rule and all cohort records; report agreement against a newly pinned release separately | The original unsaved responses cannot be reconstructed, and a provenance audit must not silently replace library selections |
