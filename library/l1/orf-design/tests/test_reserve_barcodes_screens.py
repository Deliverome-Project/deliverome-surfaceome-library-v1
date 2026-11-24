"""Barcode reservation screens: a miRNA seed made across a barcode junction is rejected.

The s1e2 generator screened each barcode for seeds on its own, so a seed can only be new
where the barcode meets its fixed flanks. These tests plant one there and check the screen
catches it, for today's 3' flank and for the T9A draft.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import pytest

S4E9 = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(S4E9))

from reserve_barcodes import (  # noqa: E402
    _SEEDS,
    FELDMAN_T9A_3PRIME,
    FWD_PRIMER,
    REV_PRIMER,
    junction_seed_hits,
    screens_pass,
)

POOL = S4E9.parents[1] / "l1/vector-design" / "barcodes_20nt_30k.tsv"
FLANKS = {"current": REV_PRIMER, "t9a-draft": FELDMAN_T9A_3PRIME}


def _clean_barcodes(three_prime: str, n: int = 200) -> list[str]:
    out = []
    with POOL.open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            bc = r["barcode"].strip().upper()
            if screens_pass(bc, three_prime)[0]:
                out.append(bc)
                if len(out) == n:
                    break
    return out


@pytest.mark.parametrize("name", sorted(FLANKS))
def test_seed_across_the_3prime_junction_is_rejected(name):
    flank = FLANKS[name]
    # a conserved seed whose last base is the flank's first, completed by the barcode's end
    seed = next(s for s in sorted(_SEEDS) if s[-1] == flank[0])
    for base in _clean_barcodes(flank):
        bc = base[:14] + seed[:6]
        if seed in bc or junction_seed_hits(bc, flank) != [seed]:
            continue
        ok, why = screens_pass(bc, flank)
        if why.startswith(("boundary", "splice")):
            continue  # the planted bases made a different defect; try another base
        assert not ok
        assert why == f"junction miRNA seed {seed}"
        return
    pytest.fail("no barcode could carry a seed that only the 3' junction completes")


def test_seed_across_the_5prime_junction_is_rejected():
    seed = next(s for s in sorted(_SEEDS) if s[0] == FWD_PRIMER[-1])
    for base in _clean_barcodes(REV_PRIMER):
        bc = seed[1:] + base[6:]
        if seed in bc or junction_seed_hits(bc) != [seed]:
            continue
        ok, why = screens_pass(bc)
        if why.startswith(("boundary", "splice")):
            continue
        assert not ok and why == f"junction miRNA seed {seed}"
        return
    pytest.fail("no barcode could carry a seed that only the 5' junction completes")


def test_windows_inside_the_barcode_are_not_rescreened():
    """Only junction-spanning windows are checked; the generator already cleared the rest."""
    seed = sorted(_SEEDS)[0]
    bc = "A" * 6 + seed + "A" * 7
    assert seed not in junction_seed_hits(bc)


def test_l1_draws_must_pass_both_flanks():
    """Until the 3' flank is chosen, an L1 barcode is rejected if either flank fails it."""
    from reserve_barcodes import L1_DRAW_FLANKS, l1_screens_pass

    assert set(L1_DRAW_FLANKS) == {REV_PRIMER, FELDMAN_T9A_3PRIME}
    # a seed completed only by the T9A flank's first base fails the L1 rule
    seed = next(
        s for s in sorted(_SEEDS) if s[-1] == FELDMAN_T9A_3PRIME[0] and s[-1] != REV_PRIMER[0]
    )
    for base in _clean_barcodes(REV_PRIMER):
        bc = base[:14] + seed[:6]
        if not screens_pass(bc, REV_PRIMER)[0] or junction_seed_hits(bc, FELDMAN_T9A_3PRIME) != [
            seed
        ]:
            continue
        ok, why = screens_pass(bc, FELDMAN_T9A_3PRIME)
        if why.startswith(("boundary", "splice")):
            continue
        assert not ok
        assert not l1_screens_pass(bc)[0]
        return
    pytest.fail("no barcode passes today's flank but fails only the T9A junction seed")
