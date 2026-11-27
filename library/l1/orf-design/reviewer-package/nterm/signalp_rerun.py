"""Optional, for reviewers holding a SignalP 6.0 licence: rerun it and compare with the pins.

SignalP 6.0 is DTU-licensed, so the package ships its outputs rather than the tool. This
script lets a licensed reviewer check those outputs.

1. ``fasta`` writes every sequence the placement used from SignalP: each "N after SP"
   member's designed protein, plus each candidate tagged design that was re-predicted.
   Every record is named by its SHA-256, which is how the pins are keyed.
2. Run SignalP 6.0 yourself. The pins were made with ``--organism eukarya --mode
   slow-sequential`` (``signalp-6.0+h``), in batches of 64.
3. ``compare`` reads SignalP's ``prediction_results.txt`` and compares each call, and its
   cleavage position, with the pinned sweep.

    uv run --locked python signalp_rerun.py fasta --out nterm_signalp_input.fasta
    signalp6 --fastafile nterm_signalp_input.fasta --organism eukarya \\
        --mode slow-sequential --output_dir sp6_out
    uv run --locked python signalp_rerun.py compare --results sp6_out/prediction_results.txt

Expect small differences in the probabilities. The record shows SignalP's numbers depend on
batch composition: a two-shard rerun moved 13 probabilities by up to 1.0e-4. A cleavage
call at the margin can also move, because the position is an argmax. Compare the calls,
and treat probability differences below 1e-3 as equal.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import nterm_placement as npl  # noqa: E402


def pinned(processed) -> dict[str, dict]:
    """sha256 -> pinned SignalP call, for every native and tagged sequence the sweep scored."""
    return {r["sequence_sha256"]: r for r in npl.load_json(processed, npl.SWEEP_REPREDICTION_KEY)}


def write_fasta(out: Path) -> int:
    from deliverome_analysis.core.s3 import processed

    members, _ = npl.load_members(processed)
    pins = pinned(processed)
    with (HERE / "nterm_positions.tsv").open() as fh:
        tagged = [r for r in csv.DictReader(fh, delimiter="\t") if r["repredicted"] == "1"]
    seqs = {m["sha"]: m["sequence"] for m in members.values() if m["sha"] in pins}
    for r in tagged:
        m = members[r["candidate_id"]]
        n = int(r["after_residue"])
        seqs[r["tagged_sha256"]] = m["sequence"][:n] + npl.CASSETTE_AA + m["sequence"][n:]
    with out.open("w") as fh:
        for sha, s in sorted(seqs.items()):
            assert npl.sha256(s) == sha
            fh.write(f">{sha}\n{s}\n")
    print(f"wrote {out}: {len(seqs)} sequences")
    return 0


def compare(results: Path) -> int:
    from deliverome_analysis.core.s3 import processed

    pins = pinned(processed)
    seen = same_call = same_site = 0
    max_dp = 0.0
    differ = []
    for line in results.read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        cols = line.split("\t")
        sha, pred = cols[0].split()[0], cols[1]
        pin = pins.get(sha)
        if pin is None:
            continue
        seen += 1
        site = re.search(r"CS pos: (\d+)-", line)
        site = int(site.group(1)) if site else None
        sp_prob = float(cols[3])
        same_call += pred == pin["prediction"]
        same_site += site == pin["cleavage_site"]
        max_dp = max(max_dp, abs(sp_prob - pin["sp_probability"]))
        if pred != pin["prediction"] or site != pin["cleavage_site"]:
            differ.append(
                {
                    "sha256": sha,
                    "pinned": [pin["prediction"], pin["cleavage_site"]],
                    "rerun": [pred, site],
                }
            )
    print(f"{seen} sequences compared with the pins")
    print(f"  same SP/OTHER call: {same_call}; same cleavage position: {same_site}")
    print(f"  largest P(SP) difference: {max_dp:.2e}")
    if differ:
        print(json.dumps(differ[:20], indent=1))
    return 0 if seen and not differ else 1


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("fasta").add_argument("--out", type=Path, required=True)
    sub.add_parser("compare").add_argument("--results", type=Path, required=True)
    a = ap.parse_args(argv)
    return write_fasta(a.out) if a.cmd == "fasta" else compare(a.results)


if __name__ == "__main__":
    raise SystemExit(main())
