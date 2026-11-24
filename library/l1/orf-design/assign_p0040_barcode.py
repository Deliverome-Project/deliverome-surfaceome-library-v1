"""Give p0040 its own barcode, so it is distinguishable from p0016 by amplicon.

p0040 was built by copying p0016 and changing four bases in the target array, which
left it carrying p0016's barcode. Two constructs sharing a barcode cannot be told apart
by the amplicon readout, which is the whole point of the barcode, so they could only
ever be run in separate wells. One barcode fixes that.

Deliberately not `reserve_barcodes.py`. That script re-draws every entry in its
CONSTRUCTS list from a fixed starting P-number and silently reassigns barcodes that are
already registered while looking like a normal run. This asks for exactly one barcode,
for one named construct, and refuses if anything about the registry is not as expected.

Reuses that script's screens rather than restating them: a boundary-enzyme site created
against this experiment's junction context, and an intrinsic splice donor.

    uv run --locked python assign_p0040_barcode.py --dry-run
    uv run --locked python assign_p0040_barcode.py
"""

from __future__ import annotations

import argparse
import csv
import sys
import tempfile
import warnings
from pathlib import Path

from Bio import BiopythonWarning, SeqIO
from Bio.Seq import Seq

sys.path.insert(0, str(Path(__file__).parent))
from reserve_barcodes import POOL_KEY, REGISTRY_KEY, read_registry, screens_pass

from deliverome_analysis import constructs, seqtools
from deliverome_analysis.s3 import processed

CONSTRUCT = "p0040"
PARENT = "p0016"
ORF = "TFRC"
NOTE = "miR-L iFFL, PGK, intron in PuroR CDS, 4x17 sites, TACC spacer (no miR-374 8mer)"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args(argv)

    used_bc, used_ix, skipped, _done, _used_p = read_registry()
    print(f"registry: {len(used_bc)} barcode(s) in use, indices used {sorted(used_ix)}")
    print(f"          indices deliberately skipped: {skipped}")

    rows = list(csv.DictReader(Path(processed.download(POOL_KEY)).open(), delimiter="\t"))
    pool = {int(r["index"]): r["barcode"].strip().upper() for r in rows if r.get("index")}

    pick = None
    for ix in sorted(pool):
        if ix in used_ix or ix in skipped:
            continue
        ok, why = screens_pass(pool[ix])
        if not ok:
            print(f"  index {ix}: rejected, {why}")
            continue
        pick = (ix, pool[ix])
        break
    if pick is None:
        raise SystemExit("no pool index passes the screens")
    ix, bc = pick
    print(f"\nnext free index that passes both screens: {ix} -> {bc}")

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", BiopythonWarning)
        rec = SeqIO.read(processed.download(constructs.map_key(CONSTRUCT)), "genbank")
        parent = SeqIO.read(processed.download(constructs.map_key(PARENT)), "genbank")
    before, parent_seq = str(rec.seq).upper(), str(parent.seq).upper()

    feat = next(f for f in rec.features if f.qualifiers.get("label", [""])[0].startswith("20bp-BC"))
    old = str(feat.extract(rec.seq)).upper()
    if before.count(old) != 1:
        raise SystemExit(f"{CONSTRUCT}: barcode {old} is not unique in the plasmid")
    if old not in parent_seq:
        raise SystemExit(f"{CONSTRUCT}: its barcode is not the parent's, so this already ran")
    if bc in used_bc:
        raise SystemExit(f"barcode {bc} is already registered")

    after = before.replace(old, bc)
    problems = []
    if (
        len(after) != len(before)
        or sum(1 for a, b in zip(before, after, strict=True) if a != b) > 20
    ):
        problems.append("more than the barcode changed")
    if after.count(bc) != 1:
        problems.append("new barcode is not unique")
    if bc in parent_seq:
        problems.append("new barcode also occurs in the parent")
    for enz in ("BbsI", "Esp3I", "BsaI"):
        b, a = (len(seqtools.find_sites(s, [enz])) for s in (before, after))
        if a != b:
            problems.append(f"{enz} count changed {b} -> {a}")
    if problems:
        raise SystemExit(f"{CONSTRUCT}: refusing to write — {problems}")

    start = after.find(bc)
    print(f"swapped {old} -> {bc} at {start + 1:,}..{start + len(bc):,}")
    print(f"  p0016 keeps {old}; the two are now distinguishable by amplicon")
    feat.qualifiers["label"] = [f"20bp-BC-{ORF} [pool index {ix}]"]
    rec.seq = Seq(after)
    row = [
        bc,
        str(ix),
        CONSTRUCT,
        ORF,
        f"{start + 1}..{start + len(bc)}",
        "s1e2 pool",
        "designed",
        NOTE,
    ]

    if args.dry_run:
        print("\nregistry row that would be appended:\n  " + "\t".join(row))
        print("\ndry run — S3 not written")
        return 0

    # registry first: its own header says reserve before assigning
    text = Path(processed.download(REGISTRY_KEY)).read_text()
    if bc in text or f"\t{ix}\t" in text:
        raise SystemExit("registry changed under us; re-run")
    with tempfile.TemporaryDirectory() as tmp:
        staged = Path(tmp) / "barcode_assignments.tsv"
        staged.write_text(text.rstrip("\n") + "\n" + "\t".join(row) + "\n")
        processed.upload(staged, REGISTRY_KEY, content_type="text/tab-separated-values")
    back = Path(processed.download(REGISTRY_KEY)).read_text()
    if back.count(bc) != 1:
        raise SystemExit(f"registry re-read shows {back.count(bc)} rows for {bc}")
    print("registry appended and re-read: exactly one row")

    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / f"{CONSTRUCT}.gb"
        SeqIO.write(rec, path, "genbank")
        processed.upload(
            path, constructs.map_key(CONSTRUCT), content_type="chemical/seq-na-genbank"
        )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", BiopythonWarning)
        reread = str(
            SeqIO.read(processed.download(constructs.map_key(CONSTRUCT)), "genbank").seq
        ).upper()
    if reread != after:
        raise SystemExit("re-read does not match what was uploaded")
    print(f"{CONSTRUCT} written and re-read; barcode is pool index {ix}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
