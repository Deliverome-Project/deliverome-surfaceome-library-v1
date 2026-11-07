"""The prefix screen keeps barcodes distinct over the bases in situ sequencing reads."""

from __future__ import annotations

import csv
import random
import sys
from pathlib import Path

import pytest

S4E9 = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(S4E9))

from barcode_prefix_screen import PrefixScreen, min_prefix_distance  # noqa: E402

REPO = S4E9.parents[2]
POOL = REPO / "analysis/sprint1/s1e2-fragmid-vector-design/barcodes_20nt_30k.tsv"
L1_MEMBERS = 2735


def _pool() -> list[str]:
    with POOL.open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return [r["barcode"].strip().upper() for r in sorted(rows, key=lambda r: int(r["index"]))]


@pytest.mark.parametrize("bases,min_distance", [(6, 2), (6, 3), (8, 3)])
def test_screen_matches_brute_force(bases, min_distance):
    rng = random.Random(7)
    cands = ["".join(rng.choice("ACGT") for _ in range(20)) for _ in range(400)]
    screen = PrefixScreen(bases, min_distance)
    kept = [c for c in cands if screen.try_accept(c)]
    assert len(kept) > 1
    assert min_prefix_distance(kept, bases) >= min_distance
    # and nothing rejected could have been kept: each rejection is within reach of a kept one
    for c in cands:
        if c not in kept:
            near = min(sum(x != y for x, y in zip(c[:bases], k[:bases], strict=True)) for k in kept)
            assert near < min_distance


def test_existing_barcodes_are_avoided():
    existing = ["ACGTACGTACGTACGTACGT"]
    screen = PrefixScreen(10, 3, existing=existing)
    assert not screen.passes("ACGTACGTAAGTACGTACGT")  # 1 substitution in the prefix
    assert screen.passes("TTTTACGTACGTACGTACGT")  # 3 substitutions


def test_pool_supports_the_l1_draw():
    """The s1e2 pool yields 2,735 barcodes at distance >= 3 over the first 10 bases."""
    screen = PrefixScreen()
    kept = []
    for bc in _pool():
        if screen.try_accept(bc):
            kept.append(bc)
            if len(kept) == L1_MEMBERS:
                break
    assert len(kept) == L1_MEMBERS
    assert min_prefix_distance(kept, 10) >= 3


def test_unscreened_draw_collides():
    """Why the screen exists: the first 2,735 pool barcodes share 10-base prefixes."""
    assert min_prefix_distance(_pool()[:L1_MEMBERS], 10) == 0
