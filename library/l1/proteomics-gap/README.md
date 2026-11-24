# Surface-specific recurrence against the FG library

`surface_specific.py` (unchanged since sweep-v11; sha256 2a72ba80…) counts publications per
protein in three kinds of S3E8 input and applies rule A. `rule_variants.py` applies rules A and B
to its per-protein counts and diffs the not-in-library lists against a baseline folder.

| Rule | Definition |
|---|---|
| A | surface ≥ 20, surface − membrane ≥ 10, surface − global ≥ 10 |
| B | surface ≥ 20 and lod_membrane > 0.5; lod_global is reported, not filtered on |

| Snapshot | Publications (surface / membrane / global) | A: pass / not in L1/L2 | B: pass / not in L1/L2 |
|---|---|---|---|
| 2026-11-15-sweep-v11 | 75 / 76 / 18 | 308 / 37 | 355 / 44 |
| 2026-11-15-sweep-v12 | 76 / 76 / 18 | 310 / 37 | 347 / 44 |
| 2026-11-16-sweep-v13 | 77 / 76 / 18 | 341 / 40 | 364 / 42 |
| 2026-11-16-sweep-v15 | 71 / 76 / 18 | 269 / 35 | 347 / 44 |
| 2026-11-16-sweep-v16 | 71 / 77 / 18 | 244 / 33 | 334 / 42 |
| 2026-11-18-sweep-v20 | 73 / 77 / 18 | 278 / 34 | 349 / 44 |
| 2026-11-18-abundance-v4 | 84 / 77 / 18 | 422 / 60 | **394 / 51** |

The not-in-library lists are the same accessions in v11 and v12 under both rules
(`v12/rule_*_diff_vs_baseline.tsv`). sweep-v13 restores George 2026: rule A adds ITGAX, PTPRB and
PTPRD; rule B drops AGRN and F2 (`v13/rule_*_diff_vs_baseline.tsv`, both against v11). sweep-v15
applies the admission rules (annotation-selected tables never vote; imputation alone does not
hold): rule A adds KIAA1549 and drops AGRN, FAT1 and LNPEP against v11 (`v15/`). v14 was withdrawn. Becca has not chosen between A and B.
Library: origin/main `library/l1/orf-design/candidate_manifest_tranches.tsv`
(sha256 33beb879…), joined on Ensembl gene via `gene_resolver`.

    uv run --locked python surface_specific.py <D> <sweep data folder>   # D holds manifest.tsv
    uv run --locked python rule_variants.py <D> <sweep data folder> [<baseline D>]

sweep-v16 adds one membrane-fractionation publication (Zimmerman et al. 2025 xenografts), which raises
membrane counts: against v15, rule A drops VTN and ABCA7 and rule B drops F2 and PLAU
(`v16/rule_*_diff_vs_baseline.tsv` are against v15).

## An abundance snapshot needs its roles overlaid first (`abundance_roles.py`)

`abundance_screen.py` re-judges each screened input by the share of its total protein *amount*
in the surface reference, and writes the re-assigned `review_role` / `admitted` into
`manifest.json` **only** — it copies `detections_long.tsv` through unchanged, because no protein
list changes. `surface_specific.py` reads the roles from `detections_long.tsv`, so pointed at an
abundance folder directly it answers for the count screen alone.

**That failure is silent and looks like a result.** Run as-is against `2026-11-18-abundance-v4`,
rule B returns the same 44 proteins as the bare sweep with an all-`kept` diff, which reads as
"the abundance screen changed nothing" rather than "the roles never arrived". It changed 16
inputs, 11 of them from `below_adjusted_screen` or `provisional` into voting roles, which is what
takes surface publications from 73 to 84. So run the overlay first:

    uv run --locked python abundance_roles.py <sweep folder> <abundance folder> <out folder>
    uv run --locked python surface_specific.py <D> <out folder>
    uv run --locked python rule_variants.py <D> <out folder> [<baseline D>]

`v20/` is sweep-v20 as published (44 under rule B); **`v20ab/` is the current answer** —
abundance-v4's roles overlaid, 51 under rule B — and its diff is against `v20/`, so it isolates
the abundance screen's effect from the two new arms'. Against `v20`, rule B adds TGFB1, CACNA2D2,
TMEM51, SBSN, CBARP, PRRT3, FAM171B and SLC22A23 and drops VTN; against the 42 at v16, PLAU and
TENM2 are also new. sweep-v20 itself rebuilds Hoare 2016 from the authors' compiled sheet
(591 proteins, was 486) and adds Soh 2020 Table S6 as its own arm, and moves rule B from 42 to 44
on its own.

Two things to carry forward. The s4e9 reviewer brief cites the v16 numbers ("347 proteins qualify
and 42 are in neither L1 nor L2"), so it is now stale, and `s3e8_gap_gate_reasons.py` has not been
re-run — its per-gene table still covers the 42, leaving the 9 genes new since v16 with no gate
reason. Join on `uniprot_accession`, never the symbol: the sweep carries UniProt entry-name
mnemonics for several rows (`TMM94`, `F234B`, `K1549`, and now `TMM51`), which truncate the gene
name to five characters.
