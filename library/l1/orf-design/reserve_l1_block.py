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

Attribution (Becca's convention of 2026-11-10): Becca Carlson decided on 2026-11-10 to
reserve 3,200 barcodes provisionally for L1 ("reserve 3200 provisionally now for the lib"),
and the RBP motif screen at her request of 2026-11-08; implemented by the Claude Code session
l1-reviewer-package, PR 283. The block regenerates from code only against the registry as it
stood before the block was appended (the appended rows now mark those barcodes as used), so
the committed table is the record; the registry rows are append-only and carry the date,
script and status in ``notes`` rather than separate attribution columns.

Writes ``l1_barcode_block_provisional.tsv`` here. ``--append`` then appends the block to the
registry with a conditional write (refused if the registry changed since it was read) and
re-reads it to verify every row landed once.

    uv run --locked python reserve_l1_block.py            # draw and write the table
    uv run --locked python reserve_l1_block.py --append   # also append to the registry

**Extension (2026-11-15).** With 2 or 3 arms per internal member L1 grew to about 3,265 inserts,
past the 3,200 block. Becca approved extending the reservation on 2026-11-15 ("Yes, extend
it"). The extension is a second appended block, not a regeneration: every pool barcode after
the block's last index (20,892) that passes the same screens, under both 3' flanks, with the
prefix screen run against the whole registry including the 3,200 block. That is every such
barcode the pool has left, so a library above the block plus extension would need a relaxed
prefix rule or a second pool. Its rows are named from ``L1-res-3201`` and recorded in
``l1_barcode_block_extension_2026-11-15.tsv``; like the block, it regenerates from code only
against the registry as it stood before it was appended. Decided by Becca Carlson
(2026-11-15); sized and implemented by the Claude Code session l1-reviewer-package at the
coordinator's request.

    uv run --locked python reserve_l1_block.py --extension            # draw, write the addendum
    uv run --locked python reserve_l1_block.py --extension --append   # also append it
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

from deliverome_analysis import rbp_motifs
from deliverome_analysis.s3 import processed

HERE = Path(__file__).parent
OUT = HERE / "l1_barcode_block_provisional.tsv"
DATE = "2026-11-10"
EXT_OUT = HERE / "l1_barcode_block_extension_2026-11-15.tsv"
EXT_DATE = "2026-11-15"
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

#: RNA-binding-protein motifs, transcript sense, from the shared screen
#: (``deliverome_analysis.rbp_motifs``): Pumilio PRE, AU-rich element (``ATTTA`` and the
#: ``TTATTTA[AT][AT]`` nonamer, which contains it), CPEB element, poly(A) signal, and CFIm25
#: ``TGTA``. The shared screen's splice-site motifs are left out: splice donors are already
#: screened by ``l1_screens_pass``, and its acceptor check is a heuristic that would reject
#: 98 barcodes of the block reserved on 2026-11-10 (all on that heuristic alone).
#: Here every hit that overlaps the barcode counts, flagged or not (a single ``ATTTA``
#: included), which is stricter than the construct gate.
#: Decided by: the reviewer-package session asked, 2026-11-10, to keep the motif set the block
#: was reserved with; the splice exclusion was the rbp-motif-screen session's choice the same
#: day, with no separate sign-off (shipped in PR #277, merged by Becca). Implemented by:
#: Claude Code session rbp-motif-screen, PR #277. The 3,200-barcode block was drawn before
#: this change and still passes it.
RBP_FAMILIES = {"Pumilio", "ARE", "CPEB", "PAS", "CFIm25"}
RBP_MOTIFS = tuple(m for m in rbp_motifs.MOTIFS if m.family in RBP_FAMILIES)
CFIM_WINDOW = rbp_motifs.CLUSTER_WINDOW
#: Prefix width L1's provisional block (2026-11-10) and extension (2026-11-15) were drawn under.
#: ``barcode_prefix_screen.PREFIX_BASES`` later moved to 12 for post-L1 reservations; pinning 10
#: here keeps L1's draw replayable (and 10-nt d >= 3 strictly implies 12-nt d >= 3).
L1_PREFIX_BASES = 10


def rbp_hits(bc: str, three_prime: str) -> list[str]:
    """RBP motifs that overlap the barcode in its context, plus CFIm25 pairs it completes.

    A pair is two ``TGTA`` within ``CFIM_WINDOW`` nt, one of them starting in the barcode or
    in the 3 nt before it (so it overlaps the barcode).
    """
    context = FWD_PRIMER + bc + three_prime + CS1
    lo, hi = len(FWD_PRIMER), len(FWD_PRIMER) + len(bc)
    hits = []
    for h in rbp_motifs.scan(context, motifs=RBP_MOTIFS, cluster_window=CFIM_WINDOW):
        if h.motif == "CFIm_TGTA":
            if h.cluster_size >= 2 and lo - 3 <= h.start < hi:
                hits.append(f"CFIm25 TGTA pair at {h.start}")
        elif h.start < hi and h.end > lo:
            hits.append(f"{h.motif} {h.match}")
    return hits


def rbp_pass(bc: str, flanks: tuple[str, ...] = L1_DRAW_FLANKS) -> tuple[bool, str]:
    for three_prime in flanks:
        hits = rbp_hits(bc, three_prime)
        if hits:
            return False, f"RBP motif {hits[0]}"
    return True, ""


def draw(
    n: int | None, after: int = -1, prefix_bases: int = L1_PREFIX_BASES
) -> tuple[list[tuple[int, str]], dict[str, int]]:
    """``n`` barcodes (``None``: every one that passes) from pool indices above ``after``."""
    used_bc, used_ix, skipped, _done, _p = read_registry()
    _local, local_bc = _locally_assigned()
    used_bc |= local_bc
    pool = {
        int(r["index"]): r["barcode"].strip().upper()
        for r in csv.DictReader(Path(processed.download(POOL_KEY)).open(), delimiter="\t")
        if r.get("index")
    }
    prefix = PrefixScreen(bases=prefix_bases, existing=used_bc)
    block, why = [], {}
    for idx in sorted(i for i in pool if i > after):
        if n is not None and len(block) == n:
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
    if n is not None and len(block) < n:
        raise SystemExit(f"pool exhausted: only {len(block)} of {n} barcodes pass every screen")
    return block, why


def rows_for(
    block: list[tuple[int, str]], first: int = 1, date: str = DATE, what: str = "block"
) -> list[dict[str, str]]:
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
                f"L1 library {what} reserved {date}; not yet assigned to a member; "
                "library/l1/orf-design/reserve_l1_block.py; "
                "unused rows are released at assignment"
            ),
        }
        for i, (idx, bc) in enumerate(block, first)
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
    ap.add_argument(
        "--extension",
        action="store_true",
        help="the 2026-11-15 extension: every passing barcode after the block, from L1-res-3201",
    )
    args = ap.parse_args(argv)
    out = EXT_OUT if args.extension else OUT

    if out.exists() and args.append:
        rows = list(csv.DictReader(out.open(), delimiter="\t"))
        print(f"appending the committed table {out.name} ({len(rows)} rows)")
    else:
        if args.extension:
            base = list(csv.DictReader(OUT.open(), delimiter="\t"))
            last = max(int(r["pool_index"]) for r in base)
            block, why = draw(None, after=last)
            rows = rows_for(block, first=len(base) + 1, date=EXT_DATE, what="block extension")
        else:
            block, why = draw(args.n)
            rows = rows_for(block)
        with out.open("w", newline="") as fh:
            w = csv.DictWriter(fh, COLUMNS, delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(rows)
        print(
            f"drew {len(block):,} barcodes from pool indices {block[0][0]:,} to "
            f"{block[-1][0]:,}; wrote {out.name}"
        )
        for k, v in sorted(why.items(), key=lambda kv: -kv[1]):
            print(f"  passed over, {k}: {v:,}")
    if args.append:
        append(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
