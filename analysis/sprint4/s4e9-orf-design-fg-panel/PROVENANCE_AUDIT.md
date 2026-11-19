# Provenance of the s4e9 processed-bucket objects

Every object `deliverome_analysis.s3` writes carries a **provenance stamp**: S3 user metadata holding
the git commit, a `git-dirty` flag (true when the tracked tree had uncommitted changes at write time),
the creation time, the writing script, and the `uv.lock` / `renv.lock` hashes. A `git-dirty: true`
stamp records that the stamped commit alone does not reproduce the object; it does not say the object
is wrong.

This audit reads those stamps across all 14,270 objects under
`processed://sprint4/s4e9-orf-design-fg-panel/`, and verifies content where git can settle it. It was
prompted by a narrower worry — that a folder-wide upload had overwritten tables belonging to other
work, as one did in [PR #234](https://github.com/Deliverome-Project/deliverome-analysis/pull/234).
That did not happen. A different gap did turn up, in the models the public gene API serves.

## What the folder-wide pass actually did

Writes on 2026-11-15 fall into four groups. Times are UTC, as S3 reports them.

| Window | Objects | Writer | Commit | `git-dirty` | Keys touched |
|---|---|---|---|---|---|
| 08:33–08:42 | 1,455 | [publish_tagged_structures.py](https://github.com/Deliverome-Project/deliverome-analysis/blob/main/analysis/sprint4/s4e9-orf-design-fg-panel/tag_refold/publish_tagged_structures.py) on 1,444 | `a22130a2` on 1,444 | true on 1,445 | `predicted-structures/tagged/` |
| scattered | 328 | 5 named scripts, 43 unstamped | 20 commits | true on 280 | refold and verdict tables |
| 16:18–16:19 | 71 | no script stamp — the folder-wide shell script | `6e28f2ea` on all 71 | false on all 71 | experiment root, `reviewer-package/`, `reviewer-package/nterm/` |
| 23:18 | 1 | `processed.upload` per key | `36e2b609` | false | `library_structure_payload.tsv` |

Across all 1,855 writes that day the stamps name 28 distinct commits, every one of them reachable
and an ancestor of `origin/main`, so no object points at a commit that was lost or never merged.

The 16:18–16:19 group is the folder-wide pass: 71 keys spanning the whole experiment root, no `script`
value, one commit. Commit `6e28f2ea` is timestamped 16:18:08, so the upload ran seconds after it.

**It published exactly what that commit held.** All 71 objects are single-part uploads, so each ETag is
the MD5 of its content; every one of the 71 equals the MD5 of its blob at `6e28f2ea`, and every one
also equals its blob on `origin/main` today. With the 23:18 payload upload that is 72 of 72 objects
faithful to the commit they claim and current with the repository. No stale local artifact reached S3,
and nothing drifted from the record.

## The models the public API serves

[library_structure_payload.tsv](https://github.com/Deliverome-Project/deliverome-analysis/blob/main/analysis/sprint4/s4e9-orf-design-fg-panel/library_structure_payload.tsv)
is the table loaded into D1 `library_structure_public`, which the gene API reads. Its 2,899 rows draw
on two sources: 2,608 AlphaFold DB entries, which carry no key of ours, and 291 ColabFold models we
made and published ourselves. A row with `is_chosen = 1` is the **chosen model** — the one the API
returns first for that gene.

All 291 ColabFold keys resolve; none is dangling. Their stamps divide as follows.

| `git-dirty` | `is_chosen` | Rows | Consequence |
|---|---|---|---|
| false | 1 | 102 | Reproducible from the stamped commit |
| true | 1 | 187 | Served as a gene's chosen model, not reproducible from its commit |
| true | 0 | 2 | Published but not served |

So **187 genes' chosen structure comes from a dirty-tree publish** — one gene each, no duplicates.
The flaw is not diffuse: it tracks the publish run exactly. Four runs put these 291 models in S3, and
the dirty flag splits them cleanly in two.

| Commit | Published | Models | `git-dirty` | Chosen models affected |
|---|---|---|---|---|
| `cbeafee9` | 2026-10-14 | 41 | false | — |
| `71fd43eb` | 2026-10-15 | 184 | true | 182 |
| `c7b7792e` | 2026-11-10 | 61 | false | — |
| `9a1c2c91` | 2026-11-11 | 5 | true | 5 |

All four commits are on `main`. Remediation therefore scopes to two runs rather than 187 scattered
objects. This predates the index refresh that prompted the audit and was not introduced by it.

The 1,444 tagged structures from the 08:33 pass carry the same flaw at larger scale, and they are not
reached by the public API: no key under `predicted-structures/tagged/` appears anywhere in the payload.
Their exposure stays inside the tag-refold workstream.

## The predictor version is not on the objects

[publish_predicted_structures.py](https://github.com/Deliverome-Project/deliverome-analysis/blob/main/analysis/sprint4/s4e9-orf-design-fg-panel/publish_predicted_structures.py)
defines `COLABFOLD_VERSION = "1.6.2-cuda12"` and prints it to the terminal after uploading. Its
docstring says the published summary carries the version. It does not: none of the five published
`fold_summary.tsv` files contains the string, and their columns are prediction scores only
(`name`, `ok`, `mean_plddt`, `length`, and either `seconds` or the mismatched-fold columns).

The constant has not changed since it entered the file in `9547adce` on 2026-10-14, before the earliest
of these publishes, so the version behind all 291 models is recoverable from code history. It is not
pinned to the artifacts, which is what a reader six months out would need.

## Dirty, and correct

A `git-dirty` stamp is about the tree, not the molecule. For these models it was never going to be
the evidence either way: `fold_inputs_mismatched.fasta`, the amino-acid sequences ColabFold was
actually run on, is gitignored (`.gitignore:139`), so no commit of this repository — clean or dirty —
reproduces a fold from source. The question the stamp raises has to be answered on the coordinates.

[verify_colabfold_models.py](https://github.com/Deliverome-Project/deliverome-analysis/blob/main/analysis/sprint4/s4e9-orf-design-fg-panel/verify_colabfold_models.py)
answers it. For every chosen ColabFold model it reads the chain out of the published PDB, takes the
member's designed protein from `library_design.designed_proteins()` — the same helper the payload
builder uses, so this compares against the protein of record rather than a second opinion about what
that protein is — and requires exact equality. It also recomputes mean pLDDT from the C-alpha
B-factors, where ColabFold writes it, and compares that against the `fold_summary.tsv` published
beside the models.

| Check | Result |
|---|---|
| Chosen models compared | 289 of 289 |
| Sequence identical to the member's designed CDS | 289 of 289 |
| Model length agrees with the payload's `ordered_length_aa` | 289 of 289 |
| Largest gap between published mean pLDDT and the coordinates | 0.000 |

The models are the right protein, and every published score is the score its own coordinates carry.
This holds for the dirty and clean publishes alike, so **the dirty stamp did not produce a wrong
model**. It is worth more than a provenance note that these are the 200 members whose AlphaFold model
was a different protein — the set where an accession-keyed shortcut does the most damage, and the one
place a silent substitution would have been most expensive.

Equality is against the designed protein as the manifest defines it *today*, which is the stronger
test: it says the published coordinates describe the protein the library will order, not merely the
protein some earlier manifest named.

## Every ColabFold row carries its score

`library_structure_payload.tsv` records a mean pLDDT for all 291 ColabFold rows. pLDDT is how a
reader decides whether to trust a model's tag-site geometry, so a chosen structure without one offers
a structure that cannot be judged.

184 rows gained a score when `plddt_map()` began reading the `fold_summary.tsv` published beside the
models in S3: the 182 chosen models of the 2026-10-15 run, plus the two `stale_fold` alternates
(TSPAN10 `l1-p2284`, SLC22A17 `l1-p1882`). Each of the 289 chosen values agrees with the mean of its
own model's C-alpha B-factors to within 0.005.

The scores were always in S3. They were missing here because `plddt_map()` read the mismatched runs
from `~/.cache/deliverome-analysis/predicted_structures_mismatched/`, so whether a score reached the
table depended on whether the machine building it was the machine that had run the fold. That is the
same class of mistake as citing a construct map in `~/Downloads`: the durable copy exists in S3 and
the code reached past it to a local one.

## What is at stake for the 187

Re-folding was the obvious remedy before the models were checked, and the check rules it out. The
models come from ColabFold's MSA server, whose databases move, so a rerun would produce new
coordinates rather than confirm these; and there is now nothing to confirm, because the published
coordinates are already verified to encode the right protein with the scores they claim. Re-folding
would replace 187 verified models with unverified ones and re-derive every tag-site call resting on
them, to buy a cleaner stamp on an artifact whose inputs git does not hold anyway.

| Option | Cost | What it buys |
|---|---|---|
| Rebuild the payload to carry the 182 scores, and record the stamp gap | One table rebuild and a D1 reload; a provenance paragraph in the experiment README | The public record stops withholding confidence for 182 genes; coordinates and every downstream number stay as published |
| Re-fold the 187 and republish | ColabFold compute; the served structures change; every tag-site call resting on them is re-derived | A clean stamp, at the price of discarding verified models for unverified ones |
| Leave both as they are | None now | The 182 genes keep a chosen model with no confidence score, and the gap stays undiscoverable |

The first option is the one the evidence supports. The defect worth fixing is the missing scores, not
the models.

## What the reload did not verify

The D1 table is verified against the payload row by row: all 2,899 rows compared on
`mean_plddt` keyed by candidate, source, set, species and model URL, with zero mismatches, no
key present on one side only, no ColabFold row left without a score, and no `alphafold_db`
row given one. That result does not depend on the public API.

Three live checks are **not done**. They were attempted and abandoned, not passed:

| Check | State |
|---|---|
| TFRC returns only `l1-p2131` and `l1-p2617` | not done |
| `/v1/internalization/TFRC` returns 200 | not done |
| `/v1/catalog` carries `intern`, `intern_lit`, `intern_lit_grade` | not done |

`api.deliverome.org` stopped answering partway through the post-reload spot checks: one run
died on an empty response and a single gene fetch timed out at 90 s. No HTTP 429 was ever
returned, so the evidence in hand does not show rate limiting; a sibling session reports every
host on the zone refusing TCP on 443 and 80 while unrelated Cloudflare hosts answer normally,
which is the shape of an IP-level block rather than a throttle. Verification stopped there
rather than adding load to a zone under investigation.

Before that point the fix was confirmed reaching real traffic for three genes — NRG3 46.46,
ABCB4 80.41, TSPAN10 71.41 — each with the public answer matching the archive-bypass answer,
so the cached answer and the handler agreed.

## Open questions

| Question | Why it matters | What would resolve it |
|---|---|---|
| Did any 2026-11-15 write supersede a legitimately newer object from another session? | This is the original clobber worry. Content equality with `6e28f2ea` and `origin/main` shows the 71 objects are right *now*, but not what stood there at 16:18 | `s3:ListBucketVersions` on `deliverome-processed`. The bucket is versioned and scientist keys cannot delete versions, so the history exists; the per-scientist IAM policy just does not grant the list action |
| Should `publish_predicted_structures.py` stamp the predictor version? | Without it, a future reader cannot tie a model to the ColabFold build that made it, even with a clean commit | A decision to add it, then a re-stamp or a documented mapping for the 291 already published |
| Do the three live checks above pass? | The D1 table is right, but a cached or handler-level fault would still hide it from the page | Running them once the zone is reachable and Becca has cleared the traffic question |
| Should `fold_inputs_*.fasta` stay gitignored? | While it is, no commit reproduces a fold from source, so sequence verification is the only check available | A decision on whether the fold inputs belong in git or are adequately durable in S3 |
| Why did the 08:33 run publish from a tree 17 h behind its base? | `a22130a2` is dated 2026-11-14 15:51 while sibling sessions had already committed `d166b60c` and `e1dd5248`; a stale base plus a dirty tree is how two sessions overwrite each other | Reading that session's record |

## Next steps

1. Record the stamp gap for the two dirty publish runs in the experiment README, so a reader meets it
   without reading this audit.
2. Grant `s3:ListBucketVersions` to the scientist policy, or run the version comparison with an admin
   credential, to close the first open question.
3. Make `publish_predicted_structures.py` write `COLABFOLD_VERSION` into the uploaded summary, or
   correct its docstring to match what it does.
4. Decide whether `scripts/upload-experiment-outputs.sh` should refuse to run on a dirty tree. It
   behaved correctly here, and the one guard already in place — keeping fold outputs in
   `~/.cache/deliverome-analysis/predicted_structures`, outside the experiment folder, so the walk
   cannot republish them — is why the structure models were untouched by it.

## Decision log

| Date | Was | Now | Why |
|---|---|---|---|
| 2026-11-15 | The 16:18–16:19 folder-wide pass was suspected of overwriting other work's tables | The pass is verified faithful: 72 of 72 objects match both their stamped commit and `origin/main` | ETag-as-MD5 comparison against the git blob at `6e28f2ea` and at `origin/main` |
| 2026-11-15 | Provenance of the published ColabFold models was unexamined | 187 of the 291 are the chosen model for a gene and carry `git-dirty: true` | Stamp read on all 291 keys referenced by the payload |
| 2026-11-15 | A dirty stamp left it open whether the 187 models were the right protein | All 289 chosen ColabFold models are sequence-identical to their designed CDS, and every published mean pLDDT matches its coordinates to 0.000 | Coordinates compared against `library_design.designed_proteins()` by `verify_colabfold_models.py` |
| 2026-11-15 | Re-folding the 187 was one of three options | Ruled out: it would replace verified coordinates with unverified ones | The models passed sequence verification, so there is nothing for a rerun to confirm |
| 2026-11-15 | `plddt_map()` read the mismatched-run scores from `~/.cache/deliverome-analysis/` | It reads the `fold_summary.tsv` published beside the models in S3, and the committed summaries | A machine-local path made the score depend on which machine built the table; 182 chosen models had a score in S3 and none here |
| 2026-11-15 | 184 ColabFold rows carried no `mean_plddt` | All 291 carry one, each chosen value within 0.005 of its own coordinates | A chosen structure with no confidence score cannot be judged by the reader it is shown to |
| 2026-11-16 | `library_structure_public` served 184 rows with no score | Reloaded from the corrected payload: 2,899 rows, 2,710 members, 2,704 chosen, no ColabFold row without a score | The table is replace-in-place and the Worker reads it with no release filter, so the fix reaches the public page only on a reload |
