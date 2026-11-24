"""Correct the pool_index recorded for p0037 and p0038 in the barcode registry.

The pool file's `index` column is the s1e2 generator's original index, and the filter
stack removed some rows: index 23 does not exist, and 7, 33, 38 and 53 are likewise
absent. An allocator that walks the file row by row and records a running counter
therefore drifts from the first gap onward. Registry rows for pool indices 1-22 all
agree with the pool; p0037 and p0038, appended 2026-10-06, do not:

    p0037  recorded 23, its barcode ATAGTCAAGAGCGCACCAGA is at pool index 25
    p0038  recorded 24, its barcode ACAGTTAGTCGAGAAGCGCA is at pool index 26

The barcode is the authoritative value - it is what goes into the plasmid - so the
index is corrected to match it, not the other way round. Left uncorrected, the registry
claims 23 and 24 are spoken for (23 does not even exist) while 25 and 26 look free, so
the next reservation collides with a barcode already in a design.

Only the pool_index cell on those two rows changes; every other byte is preserved.

    uv run --locked python fix_pool_index_drift.py --dry-run
    uv run --locked python fix_pool_index_drift.py
"""

from __future__ import annotations

import argparse
import csv
import tempfile
from pathlib import Path

from deliverome_analysis.s3 import processed

POOL_KEY = "sprint1/s1e2-fragmid-vector-design/barcodes_20nt_30k.tsv"
REGISTRY_KEY = "sprint1/s1e2-fragmid-vector-design/barcode_assignments.tsv"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args(argv)

    rows = list(csv.DictReader(Path(processed.download(POOL_KEY)).open(), delimiter="\t"))
    by_barcode: dict[str, int] = {}
    for r in rows:
        if r.get("index"):
            by_barcode.setdefault(r["barcode"].strip().upper(), int(r["index"]))

    text = Path(processed.download(REGISTRY_KEY)).read_text()
    lines = text.splitlines()
    header = next(i for i, ln in enumerate(lines) if ln.startswith("barcode\t"))
    cols = lines[header].split("\t")
    bc_col, ix_col = cols.index("barcode"), cols.index("pool_index")

    out, fixed = list(lines), []
    for i in range(header + 1, len(lines)):
        if not lines[i].strip():
            continue
        cells = lines[i].split("\t")
        bc, claimed = cells[bc_col].strip().upper(), cells[ix_col].strip()
        if not claimed or not bc:
            continue
        actual = by_barcode.get(bc)
        if actual is None:
            raise SystemExit(f"{cells[2]}: barcode {bc} is not in the pool at all")
        if int(claimed) != actual:
            cells[ix_col] = str(actual)
            out[i] = "\t".join(cells)
            fixed.append((cells[2], claimed, actual, bc))

    if not fixed:
        print("every indexed row already agrees with the pool; nothing to do")
        return 0
    for construct, claimed, actual, bc in fixed:
        print(f"  {construct}: pool_index {claimed} -> {actual}   (barcode {bc} unchanged)")

    new = "\n".join(out) + ("\n" if text.endswith("\n") else "")
    changed = sum(1 for a, b in zip(text.splitlines(), out, strict=True) if a != b)
    if changed != len(fixed):
        raise SystemExit(f"expected {len(fixed)} changed line(s), got {changed}")
    for ln_old, ln_new in zip(text.splitlines(), out, strict=True):
        if ln_old != ln_new:
            o, n = ln_old.split("\t"), ln_new.split("\t")
            if [c for k, c in enumerate(o) if k != ix_col] != [
                c for k, c in enumerate(n) if k != ix_col
            ]:
                raise SystemExit("a cell other than pool_index changed; refusing to write")

    if args.dry_run:
        print("\ndry run — S3 not written")
        return 0

    with tempfile.TemporaryDirectory() as tmp:
        staged = Path(tmp) / "barcode_assignments.tsv"
        staged.write_text(new)
        processed.upload(staged, REGISTRY_KEY, content_type="text/tab-separated-values")
    back = Path(processed.download(REGISTRY_KEY)).read_text()
    if back != new:
        raise SystemExit("re-read does not match what was uploaded")
    body = [ln for ln in back.splitlines() if ln and not ln.startswith("#")]
    for r in csv.DictReader(body, delimiter="\t"):
        if r["pool_index"] and by_barcode[r["barcode"].strip().upper()] != int(r["pool_index"]):
            raise SystemExit(f"{r['construct']} still disagrees with the pool")
    print(f"\ncorrected {len(fixed)} row(s); every indexed row now agrees with the pool")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
