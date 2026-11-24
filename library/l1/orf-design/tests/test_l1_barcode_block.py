"""The provisional L1 barcode block passes every L1 screen, offline, over the committed table."""

from __future__ import annotations

import csv
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))

from barcode_prefix_screen import PrefixScreen  # noqa: E402
from reserve_barcodes import l1_screens_pass  # noqa: E402
from reserve_l1_block import OUT, rbp_hits, rbp_pass  # noqa: E402

PANEL = [
    "barcodes.tsv",
    "barcodes_loop.tsv",
    "barcodes_tfrc.tsv",
    "barcodes_gpcr.tsv",
    "barcodes_sp_swap.tsv",
]


def _block() -> list[dict[str, str]]:
    return list(csv.DictReader(OUT.open(), delimiter="\t"))


def test_block_is_unique_and_sized():
    rows = _block()
    assert len(rows) == 3200
    assert len({r["barcode"] for r in rows}) == len(rows)
    assert len({r["pool_index"] for r in rows}) == len(rows)
    assert all(r["status"] == "provisional" for r in rows)


def test_block_passes_l1_and_rbp_screens():
    for r in _block():
        assert l1_screens_pass(r["barcode"])[0], r
        assert rbp_pass(r["barcode"])[0], r


def test_block_is_prefix_distinct_from_itself_and_the_panel():
    panel = set()
    for name in PANEL:
        p = HERE / name
        if p.exists():
            panel |= {
                (r.get("barcode") or "").strip().upper()
                for r in csv.DictReader(p.open(), delimiter="\t")
                if len((r.get("barcode") or "").strip()) == 20
            }
    screen = PrefixScreen(existing=panel)
    for r in _block():
        assert r["barcode"] not in panel
        assert screen.try_accept(r["barcode"]), r


def test_rbp_screen_catches_planted_motifs():
    clean = _block()[0]["barcode"]
    for motif in ("TGTAAATA", "ATTTA", "TTTTAT", "AATAAA"):
        planted = motif + clean[len(motif) :]
        assert rbp_hits(planted, "")[:1], motif


# ------------------------------------------------- extension (Becca, 2026-11-15)


def _extension() -> list[dict[str, str]]:
    from reserve_l1_block import EXT_OUT

    return list(csv.DictReader(EXT_OUT.open(), delimiter="\t"))


def test_extension_continues_the_block():
    block, ext = _block(), _extension()
    assert len(ext) == 444
    assert [r["construct"] for r in ext] == [f"L1-res-{i:04d}" for i in range(3201, 3645)]
    last = max(int(r["pool_index"]) for r in block)
    assert min(int(r["pool_index"]) for r in ext) > last
    both = block + ext
    assert len({r["barcode"] for r in both}) == len(both)
    assert len({r["pool_index"] for r in both}) == len(both)
    assert all(r["status"] == "provisional" for r in ext)


def test_extension_passes_l1_and_rbp_screens():
    for r in _extension():
        assert l1_screens_pass(r["barcode"])[0], r
        assert rbp_pass(r["barcode"])[0], r


def test_extension_is_prefix_distinct_from_the_block_and_itself():
    screen = PrefixScreen(existing={r["barcode"] for r in _block()})
    for r in _extension():
        assert screen.try_accept(r["barcode"]), r
