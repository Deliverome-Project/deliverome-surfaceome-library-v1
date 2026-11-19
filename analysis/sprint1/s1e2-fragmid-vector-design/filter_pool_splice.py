"""s1e2 — apply the new splice-donor filter to the existing 30k pool, preserving indices.

`f_no_splice_donor` was added to the stack on 2026-09-26, after the pool was generated.
Re-running `generate_barcodes.py` with it would renumber the whole pool, and
`pool_index` is the stable key the P-series registry
(`barcode_assignments.tsv`) uses to record what is already spoken for — p0016..p0021 and
p0022..p0035 all reference it. So the pool is **filtered in place, not regenerated**:
survivors keep their original index and the failures are recorded as permanently retired.

Filtering can only ever remove members, and a subset of a distance-separated set is still
distance-separated, so the Levenshtein-3 / Hamming-4 guarantee survives untouched.

Files after this runs:

    barcodes_20nt_30k.tsv              the RESERVABLE pool — survivors only, 28,506,
                                       original indices preserved. Draw from this.
    barcodes_20nt_30k_retired_splice.tsv   the 1,494 removed, kept so the retirement is
                                       auditable rather than a silent deletion.
    barcodes_20nt_30k_asgenerated.tsv  the untouched 30,000 as the generator emitted them.

The reservable pool deliberately keeps the `barcodes_20nt_30k.tsv` name even though it no
longer holds 30,000: the P-series registry states that `pool_index` refers to that
filename, and every assignment from p0016 onward cites an index in it. Renaming would
break those references, and filtering in place keeps every surviving index meaning exactly
what it meant before.

No top-up is planned — the original generation was deliberately generous, and 28,506
covers the library as scoped.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from barcode_filter_stack import f_no_splice_donor

HERE = Path(__file__).resolve().parent
POOL = HERE / "barcodes_20nt_30k_asgenerated.tsv"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pool", type=Path, default=POOL)
    args = ap.parse_args(argv)

    rows = list(csv.DictReader(args.pool.open(), delimiter="\t"))
    keep = [r for r in rows if f_no_splice_donor(r["barcode"])]
    drop = [r for r in rows if not f_no_splice_donor(r["barcode"])]

    print(f"pool        {len(rows):>7,}")
    print(f"retired     {len(drop):>7,}  ({len(drop) / len(rows):.2%})  intrinsic GT[AG]AG")
    print(f"retained    {len(keep):>7,}  original indices preserved")
    print(f"as-generated{len(rows):>7,}  -> reservable pool is {len(drop):,} smaller")

    out = HERE / "barcodes_20nt_30k.tsv"
    with out.open("w", newline="") as h:
        w = csv.DictWriter(h, fieldnames=["index", "barcode"], delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows({"index": r["index"], "barcode": r["barcode"]} for r in keep)

    retired = HERE / "barcodes_20nt_30k_retired_splice.tsv"
    with retired.open("w", newline="") as h:
        w = csv.DictWriter(h, fieldnames=["index", "barcode"], delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows({"index": r["index"], "barcode": r["barcode"]} for r in drop)

    print(f"\nWrote {out.name} and {retired.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
