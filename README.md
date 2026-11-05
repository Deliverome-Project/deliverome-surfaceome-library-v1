# L1 surface-protein ORF library: design review

This repository holds the design of Deliverome's L1 ORF library for external review, and
the code to rerun it. Start with the brief,
[`BRIEF.md`](analysis/sprint4/s4e9-orf-design-fg-panel/reviewer-package/BRIEF.md). The issues
still being worked on before the order are tracked in the
[cover note](analysis/sprint4/s4e9-orf-design-fg-panel/reviewer-package/COVER_NOTE.md).

## Rerun the design

You need [uv](https://docs.astral.sh/uv/) and the frozen input bundle: 85 files with a
SHA-256 each, attached to this repository's release
[`bundle-2026-11-05-r6`](https://github.com/Deliverome-Project/l1-library-review/releases/tag/bundle-2026-11-05-r6)
(1.07 GB tarball, SHA-256 `e98164e1ce1ef35aa9b81014edcb9657772a2232c1e1963b010a5482c2d840f1`).
After the download, no credentials or network access are used.

```bash
gh release download bundle-2026-11-05-r6 --repo Deliverome-Project/l1-library-review --dir ..
```

```bash
mkdir ../l1-reviewer-bundle && tar -xf ../l1-reviewer-bundle-2026-11-05-r6.tar -C ../l1-reviewer-bundle
```

```bash
uv sync --locked
```

```bash
cd analysis/sprint4/s4e9-orf-design-fg-panel/reviewer-package
uv run --locked python run_review.py --bundle ../../../../../l1-reviewer-bundle
```

The run checks the bundle, reruns the library's 15 decision stages and reports whether each
reproduces its published record, then recomputes every number in the brief. It takes
about three minutes. The package README beside the brief describes each stage and the
bundle's contents.

## Layout

The folders keep the paths they have in the lab's analysis repository, so the scripts and
the brief's links work unchanged.

| Path | What it is |
|---|---|
| `analysis/sprint4/s4e9-orf-design-fg-panel/` | the library pipeline, its committed tables, and `reviewer-package/` (brief, runner, bundle tools) |
| `analysis/sprint1/s1e3-surface-library-filters/` | the gene filter |
| `scratch/becca/appris-principal/` | proteoform selection, panels, CDS resolution and construct assembly |
| `src/deliverome_analysis/` | the shared utilities (S3 access, sequence tools, topology, plotting), complete |
| `historical/` | script versions two stages rerun as they were when their records were made |
| `EXPORT_MANIFEST.json` | which analysis-repo commit each managed file came from, with its SHA-256 |

## Adding to this repository

Add files freely. This repository is exported from Deliverome-Project/deliverome-analysis
(`export_review_repo.py`), and a re-export only updates the files listed in
`EXPORT_MANIFEST.json`. If you edit one of those, the next export leaves your version in
place and lists it, so the change can be carried back to the analysis repository by hand.
Files you add are never touched.

New utilities go in `src/deliverome_analysis/`, with dependencies declared in
`pyproject.toml` and locked in `uv.lock` (`uv add <package>`), never installed ad hoc.
