"""Barcode/ORF swap risk for a set of constructs, checked before anything is ordered.

Lentiviral reverse transcription switches between the two co-packaged genomes. Where two
library members are identical the switch is seamless, so a barcode can end up attached to
a different member's ORF. Walton et al. fit that rate against the separation between the
two member-specific elements:

    Y(x) = 1/2 (1 - e^(-2kx)),  k = 5.3e-4 bp^-1

For this library the separation is the shared tract between where the ORFs stop differing
and where the barcode begins -- the constant linker/tag cassette, the stop, the pad and
the amplicon handle. Every base of it is a place a switch goes unnoticed.

The bar is **<=100 bp (5.0%)**, with **<=80 bp (4.1%)** preferred (Becca, 2026-10-07).
The script exits non-zero above the bar so an order cannot be staged past it.

    uv run --locked python check_swap_risk.py
    uv run --locked python check_swap_risk.py --max-distance 80
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

from deliverome_analysis.cloning import seqtools

HERE = Path(__file__).parent

#: The design files that describe what would actually be ordered. The N-terminal entry is
#: the disorder-anchored rebuild, NOT ``constructs_nterm.fasta``: that file places AXL's
#: cassette inside its own signal peptide (SignalP 6 p=0.9998 calls it 1-32 against the
#: UniProt cleavage site at 25), so the construct would lose the tag to cleavage and read
#: as a non-expressor. See rebuild_disorder_anchored.py.
FASTAS = (
    "constructs_cterm.fasta",
    "constructs_nterm_disorder_anchored.fasta",
    "constructs_loop.fasta",
    "constructs_gpcr_rescue.fasta",
)
PREFERRED, LIMIT = 80, 100


def load(path: Path) -> tuple[dict[str, tuple[str, str]], list[str]]:
    """``(name -> (insert, barcode), names skipped)``, barcode read from the FASTA header.

    A construct with no ``bc=`` in its header cannot be placed in this check at all, since
    the risk is a function of the distance from the ORF to the barcode. Those names are
    returned rather than dropped: a gate that quietly measures a subset of the order is
    worse than one that fails, because it still prints a verdict.
    """
    out: dict[str, tuple[str, str]] = {}
    skipped: list[str] = []
    name = bc = None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            header = line[1:]
            name = header.split("|")[0]
            m = re.search(r"bc=([ACGT]+)", header)
            bc = m.group(1) if m else None
            if bc:
                out[f"{name} [{path.stem}]"] = (bc, "")
            else:
                skipped.append(f"{name} [{path.stem}]")
        elif name and bc:
            k = f"{name} [{path.stem}]"
            out[k] = (out[k][0], out[k][1] + line.strip())
    return {k: (seq, bc) for k, (bc, seq) in out.items()}, skipped


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--max-distance", type=int, default=LIMIT)
    args = ap.parse_args(argv)

    members: dict[str, tuple[str, str]] = {}
    skipped: list[str] = []
    missing: list[str] = []
    for f in FASTAS:
        p = HERE / f
        if not p.exists():
            missing.append(f)
            continue
        found, no_bc = load(p)
        members.update(found)
        skipped.extend(no_bc)
    if len(members) < 2:
        raise SystemExit("need at least two barcoded constructs to check")

    risk = seqtools.barcode_swap_risk(members)
    scanned = len(FASTAS) - len(missing)
    print(f"{len(members)} barcoded construct(s) across {scanned} design file(s)")
    if missing:
        print(f"  design files absent: {', '.join(missing)}")
    if skipped:
        print(
            f"  NOT COVERED -- no barcode in the header, so the ORF-to-barcode distance is\n"
            f"  undefined for {len(skipped)}: {', '.join(skipped)}"
        )
    print(f"\nworst pair: {risk.describe()}")
    print(
        f"  bar {args.max_distance} bp = {seqtools.barcode_swap_rate(args.max_distance):.1%}"
        f"   preferred {PREFERRED} bp = {seqtools.barcode_swap_rate(PREFERRED):.1%}"
    )

    if risk.distance_bp > args.max_distance:
        shortfall = risk.distance_bp - args.max_distance
        raise SystemExit(
            f"\nFAIL: {shortfall} bp over the bar. The shared tract is the constant "
            f"linker/tag cassette; codon-diversifying it per member moves the last point "
            f"of difference toward the barcode and shortens the tract. See "
            f"seqtools.diversify_library_linkers."
        )
    if risk.distance_bp > PREFERRED:
        print(f"\nPASS, but {risk.distance_bp - PREFERRED} bp over the preferred bar")
    else:
        print("\nPASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
