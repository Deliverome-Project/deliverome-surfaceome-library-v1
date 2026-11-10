"""Provisionally reserve a block of pool barcodes for the L1 library.

L1 has 2,735 active members and no barcodes yet (``candidate_manifest_tranches.tsv`` reads
``not reserved``). Assigning barcodes to members waits until membership is frozen for the
order, but the barcodes themselves are reserved now, as one block in the shared P-series
registry, so that no panel construct drawn in the meantime can take one:

    processed://sprint1/s1e2-fragmid-vector-design/barcode_assignments.tsv

The block is larger than the library (default 3,200) so members added or swapped before
the order still draw from it. Unused rows are released at assignment, by appending a
release row; the registry is append-only.

Every barcode in the block passes, in index order from the s1e2 pool:

* ``reserve_barcodes.l1_screens_pass``: no BbsI / Esp3I / BsaI site, splice donor or
  junction miRNA seed in its real context, under **both** candidate 3' flanks (today's
  RT site and the T9A padlock draft), so no barcode is redrawn whichever flank is chosen;
* an RNA-binding-protein motif screen of the same context, under both flanks: no Pumilio
  element (``TGTA[ACT]ATA``), AU-rich element (``ATTTA``), CPEB element (``TTTTAT``,
  ``TTTTAAT``) or poly(A) signal (``AATAAA``, ``ATTAAA``) that overlaps the barcode, and no
  CFIm25 pair (two ``TGTA`` within 50 nt) that a barcode ``TGTA`` completes. The fixed
  flanks carry none of these, so any hit is the barcode's own or one it makes at a
  junction;
* the in situ prefix screen (distance 3 over the first 10 bases) against every barcode
  already in the registry and every other barcode in the block.

Writes ``l1_barcode_block_provisional.tsv`` here. ``--append`` then appends the block to the
registry with a conditional write (refused if the registry changed since it was read) and
re-reads it to verify every row landed once.

    uv run --locked python reserve_l1_block.py            # draw and write the table
    uv run --locked python reserve_l1_block.py --append   # also append to the registry
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from barcode_prefix_screen import PrefixScreen
from build_constructs import CS1, FWD_PRIMER
from reserve_barcodes import (
    L1_DRAW_FLANKS,
    POOL_KEY,
    REGISTRY_KEY,
    _locally_assigned,
    l1_screens_pass,
    read_registry,
)

from deliverome_analysis.s3 import processed

HERE = Path(__file__).parent
OUT = HERE / "l1_barcode_block_provisional.tsv"
DATE = "2026-11-10"
COLUMNS = [
    "barcode",
    "pool_index",
    "construct",
    "orf",
    "position",
    "assigned_by",
    "status",
    "notes",
]

#: RNA-binding-protein motifs, transcript sense (DNA alphabet).
RBP_MOTIFS = {
    "Pumilio PRE": re.compile(r"(?=(TGTA[ACT]ATA))"),
    "AU-rich element": re.compile(r"(?=(ATTTA))"),
    "CPEB element": re.compile(r"(?=(TTTTAT|TTTTAAT))"),
    "poly(A) signal": re.compile(r"(?=(AATAAA|ATTAAA))"),
}
TGTA = re.compile(r"(?=TGTA)")
CFIM_WINDOW = 50


def rbp_hits(bc: str, three_prime: str) -> list[str]:
    """RBP motifs that overlap the barcode in its context, plus CFIm25 pairs it completes."""
    context = FWD_PRIMER + bc + three_prime + CS1
    lo, hi = len(FWD_PRIMER), len(FWD_PRIMER) + len(bc)
    hits = []
    for name, rx in RBP_MOTIFS.items():
        for m in rx.finditer(context):
            start, end = m.start(), m.start() + len(m.group(1))
            if start < hi and end > lo:
                hits.append(f"{name} {m.group(1)}")
    tgta = [m.start() for m in TGTA.finditer(context)]
    for a in tgta:
        for b in tgta:
            if a < b <= a + CFIM_WINDOW and any(lo - 3 <= p < hi for p in (a, b)):
                hits.append(f"CFIm25 TGTA pair {a}/{b}")
    return hits


def rbp_pass(bc: str, flanks: tuple[str, ...] = L1_DRAW_FLANKS) -> tuple[bool, str]:
    for three_prime in flanks:
        hits = rbp_hits(bc, three_prime)
        if hits:
            return False, f"RBP motif {hits[0]}"
    return True, ""


def draw(n: int) -> tuple[list[tuple[int, str]], dict[str, int]]:
    used_bc, used_ix, skipped, _done, _p = read_registry()
    _local, local_bc = _locally_assigned()
    used_bc |= local_bc
    pool = {
        int(r["index"]): r["barcode"].strip().upper()
        for r in csv.DictReader(Path(processed.download(POOL_KEY)).open(), delimiter="\t")
        if r.get("index")
    }
    prefix = PrefixScreen(existing=used_bc)
    block, why = [], {}
    for idx in sorted(pool):
        if len(block) == n:
            break
        bc = pool[idx]
        if idx in used_ix or idx in skipped or bc in used_bc:
            why["already in use"] = why.get("already in use", 0) + 1
            continue
        for screen in (l1_screens_pass, rbp_pass):
            ok, reason = screen(bc)
            if not ok:
                key = reason.split(" ")[0] if reason.startswith("RBP") else reason.split(" [")[0]
                key = re.sub(r" [ACGT,]+$", "", key)
                why[key] = why.get(key, 0) + 1
                break
        else:
            if prefix.try_accept(bc):
                block.append((idx, bc))
            else:
                why["prefix distance"] = why.get("prefix distance", 0) + 1
    if len(block) < n:
        raise SystemExit(f"pool exhausted: only {len(block)} of {n} barcodes pass every screen")
    return block, why


def rows_for(block: list[tuple[int, str]]) -> list[dict[str, str]]:
    return [
        {
            "barcode": bc,
            "pool_index": str(idx),
            "construct": f"L1-res-{i:04d}",
            "orf": "",
            "position": "",
            "assigned_by": "s1e2 pool",
            "status": "provisional",
            "notes": (
                f"L1 library block reserved {DATE}; not yet assigned to a member; "
                "analysis/sprint4/s4e9-orf-design-fg-panel/reserve_l1_block.py; "
                "unused rows are released at assignment"
            ),
        }
        for i, (idx, bc) in enumerate(block, 1)
    ]


def append(rows: list[dict[str, str]]) -> None:
    etag = processed.etag(REGISTRY_KEY)
    text = Path(processed.download(REGISTRY_KEY)).read_text()
    body = [ln for ln in text.splitlines() if ln and not ln.startswith("#")]
    existing = list(csv.DictReader(body, delimiter="\t"))
    have_bc = {r["barcode"] for r in existing if r["barcode"]}
    have_ix = {r["pool_index"] for r in existing if r["pool_index"]}
    have_c = {r["construct"] for r in existing}
    for r in rows:
        for field, have in (("barcode", have_bc), ("pool_index", have_ix), ("construct", have_c)):
            if r[field] in have:
                raise SystemExit(f"{r['construct']}: {field} {r[field]} already in the registry")
    lines = ["\t".join(r[c] for c in COLUMNS) for r in rows]
    updated = (text if text.endswith("\n") else text + "\n") + "\n".join(lines) + "\n"
    with tempfile.TemporaryDirectory() as tmp:
        staged = Path(tmp) / "barcode_assignments.tsv"
        staged.write_text(updated)
        processed.upload_if(
            staged, REGISTRY_KEY, if_match=etag, content_type="text/tab-separated-values"
        )
    check = Path(processed.download(REGISTRY_KEY)).read_text()
    if not check.startswith(text.rstrip("\n")):
        raise SystemExit("verification failed: existing registry bytes changed")
    bad = [ln for ln in lines if check.count(ln) != 1]
    if bad:
        raise SystemExit(f"verification failed: {len(bad)} appended row(s) not present once")
    print(f"appended {len(lines)} rows to processed://{REGISTRY_KEY}; verified")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("-n", type=int, default=3200, help="block size (default 3,200)")
    ap.add_argument("--append", action="store_true", help="append the block to the registry")
    args = ap.parse_args(argv)

    if OUT.exists() and args.append:
        rows = list(csv.DictReader(OUT.open(), delimiter="\t"))
        print(f"appending the committed block in {OUT.name} ({len(rows)} rows)")
    else:
        block, why = draw(args.n)
        rows = rows_for(block)
        with OUT.open("w", newline="") as fh:
            w = csv.DictWriter(fh, COLUMNS, delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(rows)
        last = block[-1][0]
        print(f"drew {len(block):,} barcodes from pool indices up to {last:,}; wrote {OUT.name}")
        for k, v in sorted(why.items(), key=lambda kv: -kv[1]):
            print(f"  passed over, {k}: {v:,}")
    if args.append:
        append(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
