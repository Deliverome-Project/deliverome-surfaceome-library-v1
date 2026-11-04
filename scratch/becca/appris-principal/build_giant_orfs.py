"""Back-translate each truncated giant to a CDS, domesticate it, and assemble it.

The truncations are defined as protein, since that is how the literature
specifies them, so the ORF has to be built rather than taken from a transcript.
Codons come from the human usage table, then domestication clears the Golden
Gate enzymes exactly as it does for every other candidate.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from deliverome_analysis import seqtools

ENZYMES = ["BbsI", "Esp3I", "BsaI"]
SEED = 20261011


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--truncations", default="giant_truncations.tsv")
    ap.add_argument("--out", default="giant_constructs.json")
    args = ap.parse_args()

    usage = seqtools.codon_usage()
    best = {aa: max(cods, key=cods.get) for aa, cods in usage.items()}

    out = {}
    print(f"  {'gene':12} {'aa':>5} {'nt':>6} {'sites':>6} {'edits':>6}  protein preserved")
    for r in csv.DictReader(open(args.truncations), delimiter="\t"):
        prot = r["protein"]
        cds = "".join(best[a] for a in prot) + "TAA"
        before = len(seqtools.find_sites(cds, ENZYMES))
        tamed, warnings = seqtools.domesticate(cds, ENZYMES, seed=SEED)
        edits = sum(1 for a, b in zip(cds, tamed) if a != b)
        ok = seqtools.translate(tamed).rstrip("*") == prot
        after = len(seqtools.find_sites(tamed, ENZYMES))
        if not ok or after:
            print(f"  {r['gene']:12} FAILED protein={ok} sites_after={after}")
            continue
        out[f"giant:{r['gene']}"] = {
            "construct": tamed, "sites_before": before, "sites_after": after,
            "edits": edits, "warnings": "; ".join(warnings),
            "source": f"giant_truncation/{r['confidence']}", "species": "human",
            "transcript": r["accession"], "route": r["strategy"],
            "variant_residues": "",
        }
        print(f"  {r['gene']:12} {len(prot):5} {len(tamed):6} {before:6} {edits:6}  {ok}")

    Path(args.out).write_text(json.dumps(out, indent=0, sort_keys=True))
    print(f"\n  {len(out)} giant constructs -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
