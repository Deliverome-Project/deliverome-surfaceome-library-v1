# N-terminal tag placement: reviewer package

A companion to the [L1 reviewer package](../README.md). It reruns, offline and without
Deliverome credentials, how the tag cassette is placed behind a signal peptide for the 1,300
L1 members on the "N after SP" rung. Read [`BRIEF.md`](BRIEF.md) first.

## What you need

| Item | Notes |
|---|---|
| This repository | the parent package's runner is imported (`../run_review.py`), so the folder layout must be kept |
| The frozen bundle | 9 files, about 65 MB, laid out as `s3/<external|processed>/<key>` and listed with SHA-256 in [`nterm_bundle_manifest.tsv`](nterm_bundle_manifest.tsv). Every key is distinct from the parent bundle's, except the topology snapshot and `candidate_cds.json`, which are byte-identical; it can be unpacked over the parent bundle. |
| [uv](https://docs.astral.sh/uv/) | installs the locked environment |
| No SignalP licence | SignalP 6.0 and the disorder predictors are not run; their outputs are pinned, keyed by the SHA-256 of the sequence predicted |

## Run it

```bash
uv sync --locked
```

```bash
cd analysis/sprint4/s4e9-orf-design-fg-panel/reviewer-package/nterm
uv run --locked python run_nterm_review.py --bundle /path/to/nterm-bundle
```

About 15 seconds. The runner:

1. checks every bundle file's SHA-256;
2. points the S3 helper at the bundle and makes real S3 unreachable;
3. reruns `l1_tag_paths.py` (each member's rung from DeepTMHMM v1 on its own designed
   protein) and compares it byte for byte with the committed table;
4. reruns [`nterm_placement.py`](nterm_placement.py) and compares its three outputs byte for
   byte, then prints the census agreement, the re-prediction verdict agreement and the frame
   checks.

## What `nterm_placement.py` checks

| Step | Recomputed | Checked against |
|---|---|---|
| (a) signal peptide | the SignalP 6.0 row describing each designed protein, chosen by sequence; the later of the SignalP 6.0 and DeepTMHMM v1 cleavage calls | the 2026-10-18 interval census, per member |
| (b) admissible interval | each predictor's fold onset, the median, the interval and the candidate sites ([`admissible_interval.py`](admissible_interval.py)) | the census: class, cleavage, onsets, interval and candidates, for all 1,291 shared members |
| (c) re-prediction | each candidate's tagged protein, its SHA-256, and the verdict recomputed from raw SignalP calls | the sweep's own verdicts, for all 7,929 tagged designs |
| frame | cleavage and site inside the protein; the residues at the junction; cassette DNA after nt `3n` translating to `native[:n] + cassette + native[n:]`; native and tagged sweep records matched by SHA-256 | must pass for every candidate site |

Outputs:

| File | One row per |
|---|---|
| [`nterm_members.tsv`](nterm_members.tsv) | member: route, calls, onsets, interval, class, candidates, re-prediction counts, frame result |
| [`nterm_positions.tsv`](nterm_positions.tsv) | candidate site: the residues either side, CDS nucleotide, tagged SHA-256, frame result, tagged SignalP call and verdict |
| [`nterm_numbers.json`](nterm_numbers.json) | every count in the brief |

## Optional: rerun SignalP 6.0 (licensed users)

```bash
uv run --locked python signalp_rerun.py fasta --out nterm_signalp_input.fasta
```

Run SignalP 6.0 on that file with `--organism eukarya --mode slow-sequential`, then:

```bash
uv run --locked python signalp_rerun.py compare --results sp6_out/prediction_results.txt
```

Checked against SignalP's own raw output from the sweep, the comparison reproduces all 9,046
calls and cleavage positions.

## Lab side: rebuilding the inputs

| Step | Command | Needs |
|---|---|---|
| freeze the predictor rows | `uv run --locked --group modal python snapshot_inputs.py --out DIR --publish DATE` | `.env` with read access to the surfaceome D1 database; Modal (read only) |
| collect the bundle | `uv run --locked python bundle_nterm.py collect --out DIR` | S3 credentials |
| hash and label it | `uv run --locked python bundle_nterm.py manifest --bundle DIR` | none |
| pack it | `uv run --locked python bundle_nterm.py pack --bundle DIR --out nterm-reviewer-bundle.tar` | none |
| regenerate outputs and figures | `uv run --locked python nterm_placement.py` then `uv run --locked python plot_nterm.py` | S3 credentials, or the bundle |

`snapshot_inputs.py` only reads from D1 and Modal; it publishes under this package's own
dated prefix and refuses to overwrite an existing key. If `nterm_placement.py`'s inputs
change, update the snapshot date in its constants and rerun `collect` and `manifest`.

## Sharing classes

| Class | Files | Meaning |
|---|---:|---|
| `licensed_tool_output` | 6 | SignalP 6.0 (DTU licence) outputs, and the disorder predictions frozen with them |
| `deliverome_record` | 2 | the interval census and `candidate_cds.json` |
| `deliverome_curation` | 1 | the DeepTMHMM v1 topology snapshot, also in the parent bundle |

Whether the `licensed_tool_output` files may go to outside reviewers is not yet decided. A
pack without them (`--exclude-class licensed_tool_output`) cannot reproduce steps (a)–(c).
