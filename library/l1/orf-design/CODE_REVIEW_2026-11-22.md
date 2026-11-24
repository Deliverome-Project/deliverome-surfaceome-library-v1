# L1 library code review, 2026-11-22

**Question.** Is the code exported to the private review repository
(`Deliverome-Project/deliverome-surfaceome-library-v1`) clear, modular and reproducible, and where can two
copies of one rule quietly disagree and produce a plausible wrong library?

**Scope.** Everything `reviewer-package/export_review_repo.py` exports, read at commit
`80e13804` (the tip of `claude/l1-review-export-scope` on 2026-11-22): the s4e9 rung builders
(`build_l1_*_arms.py`), the encoding issuers (`diversify_*.py`), the L1 checks (`l1_*.py`),
`synthesis_compliance.py`, `final_cut.py`, `library_pipeline.py`, `reviewer-package/*`,
`tag_refold/*`, the s5e13 internal-site scripts it pulls in, the scratch module
`library/l1/proteoform-selection/assemble_candidates.py` every builder imports, and the shared
modules in the export's import closure (`plasmid_qc`, `seqtools`, `library_design`,
`topology`, `constructs`, `tag_placement`, `membrane_side`, `genscript`, `primers` and others;
18 in all).

**Terms.**

| Term | Meaning here |
|---|---|
| Insert wrapper | The constant DNA around every insert's coding sequence: `ESP3I_5` + `KOZAK` + coding + `SPACER` + `FWD_PRIMER` + 20-N barcode slot + `REV_PRIMER` + `CS1` + `ESP3I_3` |
| Cassette | The 51-aa SpyTag003–ALFA tag, `GGGGS·SpyTag003·GGGGSGGGGS·P·ALFA·P·GGGGS` |
| p003 tag codons | SpyTag003 and P-ALFA-P keep p003's validated codons; only the four GGGGS linker blocks are recoded per arm |
| Drift | Two copies of one constant or rule that can change independently. "Differ today" means the copies already give different answers on some input |
| Pinned | A test or `--check` fails if the copies disagree |

## Findings ranked by drift risk and effort

Drift risk is high when the copy feeds the ordered DNA and nothing pins it; effort is the size
of the consolidation. Every proposed home is an existing module unless stated.

| # | Duplicated logic | Copies | Differ today? | Pinned? | Proposed home | Risk | Effort |
|---|---|---|---|---|---|---|---|
| 1 | Insert-wrapper assembly and its constants | 9 inline assemblies, constants in 3 modules | No | No | `plasmid_qc.wrap_insert()` + constants | High | Small |
| 2 | p003 tag-codon pinning and its assertion | 7 recomputations, 5 inline asserts | No | Partly | `plasmid_qc.p003_tag_codons()` | High | Small |
| 3 | Cassette peptide and tag-segment literals | 7 peptide, 4 segment copies | No | Partly | `plasmid_qc` (already holds them) | High | Small |
| 4 | Reading an insert's protein (Kozak → stop) | 3 functions, 2 hand-written codon tables | Yes (edge cases) | No | `plasmid_qc.translate_insert()` | Medium | Small |
| 5 | FASTA readers | 6 functions | Yes (header parsing, case, duplicates) | No | `seqtools.read_fasta()` | Medium | Small |
| 6 | Backbone placement check (`vector_check`) | 4 copies + `l1_pool_qc.assemble` | Return type only | No | `plasmid_qc.place_in_backbone()` | Medium | Small |
| 7 | Pinned S3 keys copied between files | order FASTA ×3, base CDS ×9, commit ×2 | No; 2 pin a superseded order | Partly | Each key owned by one module | Medium | Small |
| 8 | Igκ leader peptide literal | 2 literals + source | By design (± the +1 Asp) | No | Derive from `build_gpcr_rescue_arms.LEADER_AA` | Medium | Trivial |
| 9 | TSV, SHA-256 and suffix helpers | 16 table readers, 13 SHA helpers, 2 `common_suffix` | No | n/a | Leave, or `seqtools` for `common_suffix` | Low | Trivial |
| 10 | `issued()` / `issue()` in two diversifiers | 2 parallel implementations | Different domains | Each has `--check` | Leave | Low | Medium |

### 1. Insert-wrapper assembly

Every rung builder writes the wrapper out by hand:

| File | Lines |
|---|---|
| `build_l1_internal_arms.py` | 183–193 |
| `build_l1_nterm_arms.py` | 149–159, 171–181 (after a seam edit), 227–239 (`seam_fix.clean`) |
| `build_l1_cterm_arms.py` | 551–561 |
| `build_l1_multipass_arms.py` | 332–342 |
| `build_l1_type3_arms.py` | 307–317 |
| `build_l1_rescue_arms.py` | 219–229 |
| `l1_pool_qc.py` | 103–104 (`READOUT_BLOCK`, `OPENING`, the parser's view of the same layout) |

Two shorter forms are deliberate and stay: `assemble_candidates.py:47` and
`synthesis_compliance.py:749` assemble the pre-tag intermediate (no barcode block).

The constants live in three places with no test tying them together:
`build_constructs.py:53–75` (KOZAK, STOP, SPACER, FWD_PRIMER, REV_PRIMER, CS1),
`src/deliverome_analysis/plasmid_qc.py:34–40` (the same six plus `BARCODE_LEN`), and
`library/l1/proteoform-selection/assemble_candidates.py:25–36` (KOZAK, REV_PRIMER, CS1, and the
only copies of `ESP3I_5` / `ESP3I_3`). Builders take the core from `build_constructs`,
`l1_pool_qc` from `plasmid_qc`, `synthesis_compliance` from `assemble_candidates`. All three
agree today. A one-character change to any copy would build one set of inserts and check them
against another layout, and `l1_pool_qc` would report the mismatch as a construct defect rather
than a code defect.

The "exactly two Type IIS sites, both the wrapper's Esp3I" check that follows each assembly is
also repeated, in two forms: a returned error string (internal, nterm) and an assert (cterm,
multipass, type3, rescue). The rule is the same; only the failure path differs, and that
difference is meaningful (internal and nterm record a rejected site; the others never expect
one), so the consolidation keeps each caller's failure path and shares only the predicate.

**Fix.** `plasmid_qc` gains `ESP3I_5`, `ESP3I_3`, `BARCODE_SLOT`, `wrap_insert(coding)` and
`wrapper_sites_only(insert, enzymes)`; `build_constructs` and `assemble_candidates` re-export
from it; a test asserts the assembled layout byte for byte.

### 2. p003 tag-codon pinning

`seqtools.pin_codons(CASSETTE, CASSETTE_AA, TAG_SEGMENTS)` is recomputed at
`build_l1_internal_arms.py:253`, `build_l1_nterm_arms.py:281`, `build_l1_cterm_arms.py:734`,
`diversify_arm_cassettes.py:149` and `:172`, `l1_swap_risk.py:348` and `:389`, while
`diversify_cassettes.py:93` already wraps it as `tag_codons()` (used by `build_l1_rescue_arms.py:128`).
The check that an encoding keeps those codons is written inline at
`build_l1_internal_arms.py:174`, `build_l1_nterm_arms.py:134`, `diversify_arm_cassettes.py:151`
and `:196`, and `build_l1_cterm_arms.py:539`; `plasmid_qc` enforces the same rule a third way
(`require_p003_tag_codons`). This is the rule Becca set on 2026-10-29, and it is the one an
inconsistent copy would break silently: an arm whose tag codons drift still translates to the
right protein, so every protein-level check passes.

**Fix.** `plasmid_qc.p003_tag_codons()` and `plasmid_qc.keeps_p003_tag_codons(nt)`; all callers
use them; a test pins the codon map against `P003_CASSETTE_NT`.

### 3. Cassette peptide and tag-segment literals

The 51-aa peptide is typed out at `build_constructs.py:64`, `plasmid_qc.py:64`,
`tag_refold/build_refold_inputs.py:64`, `type3_nterm.py:51`,
`s5e13/build_internal_tagged.py:42` and `s5e13/l1_internal_additions.py:75`. The tag segments
are typed at `diversify_cassettes.py:82` (`TAG_SEGMENTS`), `plasmid_qc.py:83` (`TAG_PEPTIDES`),
`build_refold_inputs.py:65–66` (ALFA without its flanking prolines) and `annotate_insert.py:64–65`.
Pinned today: `tests/test_l1_pool_qc.py:45` (`P003_CASSETTE_NT == build_constructs.CASSETTE`),
`tests/test_l1_library_integrity.py:287` (that DNA is in the p003 map), and
`build_internal_tagged.py:68` (a source-text match against `build_constructs`). The refold and
type3 copies are pinned by nothing.

**Fix.** Every copy imports `plasmid_qc.P003_CASSETTE_PEPTIDE` / `TAG_PEPTIDES`;
`diversify_cassettes.TAG_SEGMENTS` becomes `tuple(qc.TAG_PEPTIDES.values())`. s5e13 files are
outside the main rebuild's set only if the main session confirms it.

### 4. Reading an insert's protein

| Function | Start | Unknown codon | Failure |
|---|---|---|---|
| `l1_consistency.translate_insert` (89) | after the full 22-nt KOZAK, asserts ATG | `seqtools.translate` rules | AssertionError |
| `tag_refold/build_refold_inputs.translate_insert` (176) | after `GCCGCCACC`; assumes the next codon is Met without checking | `X` | ValueError if no stop |
| `reviewer-package/fold_order_check.translate` (79) | after `GCCGCCACC`, requires ATG | KeyError (e.g. an `N`) | ValueError |

The two tag_refold/reviewer copies each carry their own codon table (`build_refold_inputs.py:166`,
`fold_order_check.py:64`) instead of `seqtools.translate`. They agree on every insert in the
current order; they differ on an insert with a second `GCCGCCACC` upstream, a non-ATG start, or
an ambiguous base, which are exactly the inputs a check exists to catch.

**Fix.** One `plasmid_qc.translate_insert(insert)` on the full Kozak with an ATG assertion,
built on `seqtools.translate`; the three callers use it.

### 5. FASTA readers

| Function | Key | Upper-cases | Duplicate header |
|---|---|---|---|
| `diversify_arm_cassettes.read_fasta` (70) | up to first `\|` or space | no | asserts |
| `l1_pool_qc.read_fasta` (124) | whole header, list of pairs | yes | kept |
| `build_constructs.read_fasta` (98) | up to first `\|` only | no | silently overwritten |
| `l1_consistency.fasta` (74) | up to first `\|` or space | no | silently overwritten |
| `tag_refold/build_refold_inputs.fasta` (199) | whole header, list of pairs | yes | kept |
| `tag_refold/verify_ordered_constructs.fasta` (78) | up to first `\|` or space | yes | silently overwritten |

`src/deliverome_analysis/search_fasta.parse_fasta` exists but is not in the export closure.
**Fix.** `seqtools.read_fasta(path)` returning ordered `(header, seq)` pairs and
`seqtools.fasta_by_id(path)` keyed on the first token that refuses duplicates; each caller
keeps its own case handling explicitly.

### 6. Backbone placement check

`vector_check` at `build_l1_multipass_arms.py:363`, `build_l1_cterm_arms.py:653`,
`build_l1_nterm_arms.py:260`, `build_l1_type3_arms.py:350` splices each insert into p0039 at
`assemble_candidates.SLOT_START/SLOT_LEN` and counts sites; `l1_pool_qc.assemble` (215) does
the same splice. The copies differ only in return type. All count sites rather than naming
them, so an extra BsaI that replaces a lost Esp3I would pass.
**Fix.** `plasmid_qc.place_in_backbone(backbone, insert)` and one `vector_check` that compares
site names, in a shared builder module (see Modularity).

### 7. Pinned keys copied between files

| Key | Copies | State |
|---|---|---|
| Order of record `synthesis/2026-11-21-snorkel-modules/l1_order_inserts_recoded.fasta` | `synthesis_compliance.py:104` (owner), `reviewer-package/fold_order_check.py:47`, `audit_domestication_barcodes.py:763` | agree; the next order cut must edit all three |
| Superseded order `2026-11-13-one-arm-renamed` | `tag_refold/build_refold_inputs.py:103`, `tag_refold/verify_ordered_constructs.py:51` | historical inputs of the refold run; `verify_ordered_constructs` reads as a live check but checks the 11-13 order |
| Base CDS record `library-candidates/candidate_cds.json` | `library_design.BASE_CDS_KEY` (owner), `l1_swap_risk.py:140–147` (×5), `build_refold_inputs.py:137` | literal copies of a constant that exists |
| `MAIN_COMMIT = "5b133d0b"` | `build_refold_inputs.py:132`, `native_models.py:52` | agree |

`build_refold_inputs.py:376–386` also reassigns the module global `ld.CDS_KEY` in a loop to read
older records. It restores it in `finally`, but any other caller of `library_design` running in
the same process during that loop reads the wrong record. **Fix.** Import the owner's constant;
pass the key as an argument to `ld.designed_proteins` instead of patching the global.

### 8. Igκ leader

`synthesis_compliance.py:124` (`METDTLLLWVLLLWVPGSTGD`, leader plus the +1 Asp, protected from
recoding), `tag_refold/verify_ordered_constructs.py:53` (`METDTLLLWVLLLWVPGSTG`, stripped to
get the mature chain) and `type3_nterm.py:59` (leader and Asp as two literals) all copy `build_gpcr_rescue_arms.LEADER_AA` (+ `PLUS1_AA`). The one-residue
difference is correct for each use, which is why it should be written as
`LEADER_AA + PLUS1_AA` and `LEADER_AA`, not as two literals.

### 9–10. Low-risk repetition

Five identical `read_table` in the builders, three `tsv`, three `_tsv`, three `read`; eight
string `sha` and five file `sha256`; `common_suffix` at `l1_swap_risk.py:75` and
`l1_swap_pool_model.py:204`; `FLANK = 150` in six scripts; `TAIL = STOP + SPACER + FWD_PRIMER`
at `build_l1_cterm_arms.py:98`, `l1_swap_risk.py:69` and inline at `diversify_cassettes.py:131,212`.
These are one-liners whose copies cannot disagree in a way that changes the library, so moving
them buys tidiness, not safety. `TAIL` and `common_suffix` are the exception worth taking: they
define the swap-tract measure, so they go with finding 1 and to `seqtools` respectively.
The per-script `SEED` values differ by design (each issue is its own draw) and stay.

## Modularity

- **Library constants live in CLI scripts.** Rung builders import the cassette and wrapper from
  `build_constructs.py` (the 2026-10 pilot builder, a CLI), the Igκ leader from
  `build_gpcr_rescue_arms.py`, `TAG_SEGMENTS` from `diversify_cassettes.py`, and `ENZYMES`,
  `ESP3I_5/3`, `SLOT_START/LEN` from `library/l1/proteoform-selection/assemble_candidates.py`.
  The order's DNA therefore depends on a scratch module reached by `sys.path.insert`.
  Findings 1–3 move these into `plasmid_qc`; the scripts keep re-exports so nothing else changes.
- **`sys.path` hacks.** 31 `sys.path.insert` calls across the scope (5 builders, `l1_pool_qc`,
  `synthesis_compliance`, 8 in `tag_refold/`, 9 in `reviewer-package/`, 2 in
  `s5e13/l1_internal_arms.py`, which also inserts `src/` ahead of the installed package).
  The scratch-path ones go away with findings 1–3; the folder-local ones are the repo's pattern
  for experiment scripts and are not worth changing in a review export.
- **Scripts that are also libraries.** `build_constructs`, `build_gpcr_rescue_arms`,
  `diversify_cassettes`, `diversify_arm_cassettes` (`arm_cassette`),
  `diversify_snorkel_modules` (`module_dna`), `l1_swap_risk` (`GATE`, `domesticate_junction`),
  `l1_consistency` and `l1_arm_ids` are imported by other scripts. Importing them is cheap
  (no work at import time was found), so this is a clarity cost, not a correctness one.
- **A shared builder module.** After findings 1, 2 and 6, what the six builders still repeat is
  `read_table`, the decision-table loading, and the TSV/FASTA renderers. A local
  `l1_build_common.py` in s4e9 would hold those; worth it only once the in-flight rebuild lands.

## Clarity

- **Functions too long to review in one sitting** (≥ 100 lines): `build_refold_inputs.build`
  (444), `s5e13/l1_internal_arms.main` (364), `reviewer-package/brief_numbers.measure` (300),
  `l1_completeness.build` (252), `signal_peptide_figures.draw_examples` (243),
  `score_tagged_folds.score` (197), `l1_swap_pool_model.main` (192), `l1_consistency.build`
  (147), `l1_tag_paths.main` (141), `l1_pool_qc.run` (133), every rung builder's `build`
  (98–126). The builders shrink by ~15 lines each with findings 1 and 6.
- **Stale comments.** Since per-arm linker encodings (2026-11-21) only the tag codons are p003's.
  `build_l1_multipass_arms.py:23` still says "cassette (51 aa, p003 codons)" and
  `build_l1_nterm_arms.py:165` says "the cassette's own codons (p003's) are never touched";
  both should say the tag codons are p003's and the linkers are the arm's own encoding. The
  other "p003" mentions (tag codons, the p003 tandem) are accurate. The one `l1-m` id
  (`l1_completeness.py:90`) quotes a historical record and is correct as written.
- **A committed table on retired ids.** `type3_arms.fasta` still names its 132 records by
  the retired `l1-m…` member ids (44 members); `type3_nterm.py` now writes the same sequences
  under `l1-p…` ids, so rerunning it rewrites the table. Its other output,
  `type3_nterm_features.tsv`, is current. Regenerating the table belongs with the rebuild.
- **`verify_ordered_constructs.py`** names itself a verification of the ordered constructs but
  pins the 2026-11-13 order (finding 7). Its docstring should say which order it verifies.

## Reproducibility

| Check | Finding |
|---|---|
| Network calls | Every fetching script in scope is listed in `external_inputs.tsv`; `tests/test_external_inputs.py` enforces it. No gap found. |
| RNG | Every draw is seeded (`SEED` per issuer, `random.Random(SEED)`). No gap found. |
| `--check` modes | All six builders, all four diversifiers, `l1_consistency`, `l1_completeness`, `synthesis_compliance` have one; `l1_pool_qc` has `--no-write`. **No `--check`:** `l1_proteoform_decisions`, `l1_insert_lengths`, `l1_swap_pool_model`, `l1_structure_gap` (writes outputs in place, not in `library_pipeline`); `l1_tag_paths` and `l1_cterm_decisions` are replayed by pipeline stages instead. |
| `--check` that writes | `l1_consistency.py --check` (main, 418–432) rewrites `l1_consistency.tsv` and its summary before deciding the exit code, so running the check changes the committed record; `l1_completeness` has a separate `--no-write`. Every other `--check` in scope compares without writing. The fix is to write only without `--check` and to compare against the committed file under it. |
| Broken `--check` at HEAD | `diversify_loop_cassettes.py --check` fails: "ACKR3_igk: encoding does not translate to the loop tag". Pre-existing; not investigated here. |
| Warm-state dependency | `l1_consistency --check` and `l1_completeness --check` stop with `NativesUnavailable` for the 24 co-receptor additions (l1-p3153…) until their native proteins are committed; expected mid-rebuild, but it means the consistency check is not runnable from the branch alone right now. |
| Ordering | Builders iterate `sorted(...)` over ids; no dict-order dependence found in the written tables. |
| Global mutation | `build_refold_inputs` patches `ld.CDS_KEY` (finding 7). |

## Tests: what guards what

| Invariant | Guard | Gap |
|---|---|---|
| Cassette DNA is p003's | `test_l1_pool_qc.py:45`, `test_l1_library_integrity.py:287` | peptide copies in refold, type3, s5e13 |
| Issued encodings keep p003 tag codons and translate | `test_l1_library_integrity.py:113–134` | the shared pinning function itself |
| No Type IIS beyond the wrapper | `test_l1_library_integrity.py:173, 380` | none |
| Wrapper constants equal across modules | none | finding 1 |
| One FASTA/insert-reading behavior | `test_l1_checks.py:34–44` (only `l1_consistency.translate_insert`) | findings 4–5 |
| Network inputs on record | `test_external_inputs.py` | none |
| Strict typing | mypy strict + ty on the 18 exported shared modules (`pyproject.toml:228–250`) | scripts are untyped by design |

## What phase 2 changes, and how it is proven

Each consolidation is one commit, moves logic into the module named above, leaves every script's
public names importable, and is proven output-identical by rerunning the affected `--check`
modes and comparing written outputs byte for byte, plus a test pinning the consolidated
invariant. The files touched depend on which ones the in-flight rebuild on
`claude/l1-review-export-scope` is not editing; the shared-module changes (`plasmid_qc`,
`seqtools`) and the tag_refold/reviewer-package readers are the likeliest to be free.

## Status of the fixes

| Finding | Done (branch `claude/l1-code-quality-review-8456e1`) | Waits for the rebuild to finish |
|---|---|---|
| 1 Wrapper | `plasmid_qc.wrap_insert`, `wrapper_sites_only`, `ESP3I_5/3`, `BARCODE_SLOT`, `CODING_START`, `STOP_TO_BARCODE`, slot and enzyme constants; every copy pinned by `tests/test_l1_shared_definitions.py` | builders, `l1_pool_qc`, `build_constructs`, `assemble_candidates` switch to them |
| 2 Tag codons | `plasmid_qc.p003_tag_codons`, `keeps_p003_tag_codons`; recomputations pinned equal | callers switch |
| 3 Cassette and tags | `type3_nterm`, `annotate_insert` import from `plasmid_qc` (outputs byte-identical); refold and s5e13 literals pinned | `build_refold_inputs`, s5e13 files |
| 4 Insert translation | `plasmid_qc.translate_insert`; the three copies pinned equal on all 3,304 committed inserts | callers switch |
| 5 FASTA readers | `seqtools.read_fasta`, `fasta_id`, `fasta_by_id`; the six copies pinned equal on the files each reads | callers switch |
| 6 Backbone check | `plasmid_qc.place_in_backbone`; `l1_pool_qc.assemble` pinned equal | `vector_check` copies |
| 7 Pinned keys | — | all |
| 8 Igk leader | `type3_nterm` imports it; the two remaining literals pinned to `build_gpcr_rescue_arms` | `synthesis_compliance`, `verify_ordered_constructs` |
| 9 Suffix | `seqtools.common_suffix`; both copies pinned equal | callers switch |
| `l1_consistency --check` writes | Fixed on `agy/l1-review-export-and-dedup` (`l1_consistency` also switched to `seqtools.fasta_by_id` + `plasmid_qc.translate_insert`) | — |
| `s5e13` export closure in `export_review_repo.py` | Fixed on `agy/l1-review-export-and-dedup`: `sibling_modules()` exports the 34 transitive `s5e13` sibling `.py` helpers; pinned by `tests/test_export_review_repo.py` | — |
| `synthesis_compliance` & `plot_rescue_selection` | Fixed on `agy/l1-review-export-and-dedup`: `synthesis_compliance` uses `plasmid_qc` (drops `scratch/` `sys.path` import), `seqtools.reverse_complement`, and `gpcr_arms.LEADER_AA + PLUS1_AA`; `plot_rescue_selection.swap_rate` delegates to `seqtools.barcode_swap_rate` | — |

## Follow-up notes (2026-11-22, second pass)

- **Merge note for the T7/Feldman readout update (`claude/l1-review-export-scope`):** PR #434 added
  `plasmid_qc.wrap_insert()` and `STOP_TO_BARCODE` using the pre-T7 readout
  (`SPACER + FWD_PRIMER + barcode + REV_PRIMER + CS1`). Once the uncommitted T7/Feldman readout
  (`T7_PROMOTER` + `FELDMAN_3PRIME`, `build_constructs.l1_readout`) lands, update
  `plasmid_qc.wrap_insert()`, `STOP_TO_BARCODE`, and `tests/test_l1_shared_definitions.py` in
  lockstep before switching the five `build_l1_*_arms.py` builders over.
- **`external_inputs.py` scan scope:** `external_inputs.py` scans only `s4e9/`. Now that the
  `s5e13` sibling closure is exported to `deliverome-surfaceome-library-v1`, note that
  `s5e13/hard_negative_benchmark.py:132` makes a live `urllib.request.urlopen` call outside that
  scan scope; track it in `external_inputs.tsv` if `s5e13` is replayed offline.
- **Pinned order FASTA & Compara-116 SHA-256s:** `synthesis_compliance.FASTA_KEY`,
  `fold_order_check.FASTA_KEY`, and `audit_domestication_barcodes.py` share the order-of-record
  key, and `adopt_ortholog_cds_policy.py` and `review_ortholog_correspondence.py` share six
  Compara-116 archive pins; `tests/test_l1_shared_definitions.py` pins these equal so an order
  re-cut cannot update one without the others.

## Decision log

- 2026-11-22: review written from `80e13804` at Becca's request ("review our code quality in the
  l1 repo … minimal repeated functions that can cause drift"); fixes wait on the main session's
  list of files the rebuild is not touching.
- 2026-11-22: the main session freed the shared modules, `type3_nterm.py`, `annotate_insert.py`
  and new s4e9 tests; everything else waits until the order is re-cut. The consolidation PR is
  stacked on `claude/l1-review-export-scope` and retargets `main` when that merges.
- 2026-11-22 (PR follow-up `agy/l1-review-export-and-dedup`): fixed the broken `s5e13` sibling
  import closure in `export_review_repo.py` (+ `tests/test_export_review_repo.py`), switched
  `l1_consistency.py` to `seqtools.fasta_by_id` / `plasmid_qc.translate_insert` and made its
  `--check` non-writing, switched `synthesis_compliance.py` off `scratch/` `assemble_candidates`
  to `plasmid_qc` + `seqtools.reverse_complement` + `build_gpcr_rescue_arms`, and switched
  `plot_rescue_selection.swap_rate` to `seqtools.barcode_swap_rate`.
