"""Domesticate every resolved candidate CDS and give it a library-candidate ID.

IDs are l1p0001 onwards -- a namespace of their own, so ~2,700 library
candidates do not swamp the P-series, which is for constructs that actually get
built. Four digits because the set runs past 999.

Only the 25 expression-rescue arms carry the Igk leader; every other candidate
is the plain domesticated ORF. Resumable: finished constructs are kept, so an
interrupted run picks up where it left off.
"""

from __future__ import annotations

import argparse
import csv
import json
import signal
import sys
from pathlib import Path

from deliverome_analysis.cloning import seqtools

ENZYMES = ["BbsI", "Esp3I", "BsaI"]
SEED = 20261011
_stop = False


def _on_sigint(_sig, _frm):
    global _stop
    _stop = True
    print("\n  interrupt received - finishing the current construct, then saving")


def domesticate_one(entry: dict, *, avoid_harmful_motifs: bool = False) -> dict | None:
    """The construct record for one resolved CDS, or None if domestication would change
    its protein.

    ``avoid_harmful_motifs`` defaults to False because ``candidate_constructs.json`` was
    domesticated before ``seqtools.domesticate`` added post-domestication motif scrubbing
    (2026-11-24); motif checks for the ordered library run downstream in
    ``synthesis_compliance.py``.
    """
    raw = entry["cds"]
    sites_before = len(seqtools.find_sites(raw, ENZYMES))
    seq, warnings = seqtools.domesticate(
        raw, ENZYMES, seed=SEED, avoid_harmful_motifs=avoid_harmful_motifs
    )
    if seqtools.translate(raw) != seqtools.translate(seq):
        return None
    return {
        "construct": seq,
        "sites_before": sites_before,
        "sites_after": len(seqtools.find_sites(seq, ENZYMES)),
        "edits": sum(1 for a, b in zip(raw, seq) if a != b),
        "warnings": "; ".join(warnings),
        "source": entry["origin"],
        "species": entry["species"],
        "transcript": entry["transcript"],
        "route": entry.get("route", ""),
        "variant_residues": entry.get("variant_residues", ""),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cds", default="candidate_cds.json")
    ap.add_argument("--rescue", default="rescue_constructs.tsv")
    ap.add_argument("--out", default="candidate_constructs.json")
    ap.add_argument("--limit", type=int, default=None)
    args = ap.parse_args()

    signal.signal(signal.SIGINT, _on_sigint)
    cds = json.loads(Path(args.cds).read_text())
    done = json.loads(Path(args.out).read_text()) if Path(args.out).exists() else {}

    # The rescue arms are already built as leader fusions; reuse them verbatim.
    rescue = {}
    if Path(args.rescue).exists():
        for r in csv.DictReader(open(args.rescue), delimiter="\t"):
            rescue[f"rescue:{r['symbol']}"] = r["construct"]

    giants = {}
    gp = Path("giant_constructs.json")
    if gp.exists():
        giants = json.loads(gp.read_text())
        done.update(giants)        # already domesticated by build_giant_orfs.py
    keys = sorted(cds) + sorted(rescue) + sorted(giants)
    todo = [k for k in keys if k not in done]
    if args.limit:
        todo = todo[: args.limit]
    print(f"  {len(keys)} candidates, {len(done)} already done, {len(todo)} to do")

    for n, key in enumerate(todo, 1):
        if _stop:
            break
        if key in rescue:
            seq = rescue[key]
            done[key] = {
                "construct": seq,
                "sites_before": None,
                "sites_after": len(seqtools.find_sites(seq, ENZYMES)),
                "edits": None,
                "warnings": "",
                "source": "rescue_arm",
                "species": "human",
                "transcript": "",
                "route": "",
                "variant_residues": "",
            }
        else:
            rec = domesticate_one(cds[key])
            if rec is None:
                print(f"  {key}: domestication changed the protein - skipped")
                continue
            done[key] = rec
        if n % 100 == 0:
            Path(args.out).write_text(json.dumps(done, indent=0, sort_keys=True))
            print(f"  {n}/{len(todo)} done ({len(done)} total)", flush=True)

    Path(args.out).write_text(json.dumps(done, indent=0, sort_keys=True))
    clean = sum(1 for v in done.values() if v["sites_after"] == 0)
    print(f"\n  {len(done)} constructs stored, {clean} clear all three enzymes")
    if len(done) == len(keys):
        print("  complete")
    else:
        print(f"  {len(keys) - len(done)} remaining - rerun to continue")
    return 0


if __name__ == "__main__":
    sys.exit(main())
