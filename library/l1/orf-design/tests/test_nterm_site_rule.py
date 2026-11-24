"""The N-after-SP site rule: the admissible interval and the re-prediction-gated pick order."""

from __future__ import annotations

import sys
from pathlib import Path

NTERM = Path(__file__).resolve().parents[1] / "reviewer-package" / "nterm"
sys.path.insert(0, str(NTERM))

import admissible_interval as ai  # noqa: E402
import nterm_build_options as nbo  # noqa: E402


def _track(n: int, ordered_from: int) -> list[float]:
    """Disordered (0.9) up to residue ``ordered_from - 1``, ordered (0.1) from there on."""
    return [0.9 if i + 1 < ordered_from else 0.1 for i in range(n)]


def test_onset_is_first_residue_opening_five_ordered():
    assert ai.structure_onset(_track(80, 31), 20) == 31
    # a run of four ordered residues is not an onset
    t = _track(80, 40)
    t[24:28] = [0.1] * 4
    assert ai.structure_onset(t, 20) == 40


def test_interval_is_median_onset_minus_one_plus_tolerance_upper_middle():
    tracks = {"a": _track(80, 30), "b": _track(80, 34), "c": _track(80, 40)}
    iv = ai.tag_interval(20, tracks)
    assert iv.hi == 34 - 1 + ai.ONSET_TOLERANCE
    two = {"a": _track(80, 30), "b": _track(80, 34)}
    assert ai.tag_interval(20, two).hi == 34 - 1 + ai.ONSET_TOLERANCE  # upper middle


def test_fallback_interval_ends_at_c_plus_width_plus_tolerance():
    iv = ai.tag_interval(20, {})
    assert iv.hi == 20 + ai.FALLBACK_WIDTH + ai.ONSET_TOLERANCE


def test_no_room_gives_exactly_two_candidates():
    """Onset at c + 1 leaves c + 1 and c + 2: the 71 no-room Igκ swaps had nothing else to try."""
    iv = ai.tag_interval(20, {"a": _track(80, 21)})
    assert ai.candidate_positions(iv, 20) == [21, 22]


def test_candidates_exclude_the_cleavage_site_and_cap_at_twenty():
    iv = ai.tag_interval(20, {"a": _track(200, 120)})
    cands = ai.candidate_positions(iv, 20)
    assert cands[0] == 21 and cands[-1] == 20 + ai.CANDIDATE_SPAN


def test_cleavage_is_the_later_call():
    assert ai.cleavage_lo(22, 24) == 24
    assert ai.cleavage_lo(None, 24) == 24
    assert ai.cleavage_lo(None, None) is None


def _pos(n, strict, cp, verdict="preserved"):
    return {
        "after_residue": str(n),
        "within_strict_interval": "1" if strict else "0",
        "tagged_cleavage_probability": str(cp),
        "verdict": verdict,
    }


def test_pick_order_fold_then_softening_then_probability_then_offset():
    native = 0.95
    before_soft = _pos(25, True, 0.80)  # before the fold but softened by 0.15
    inside_firm = _pos(30, False, 0.95)
    assert nbo.gated_pick([before_soft, inside_firm], native) is before_soft
    firm = _pos(26, True, native - nbo.SOFTENING_BOUND)  # exactly at the bound is not softened
    assert nbo.gated_pick([before_soft, firm], native) is firm
    a, b = _pos(23, True, 0.93), _pos(22, True, 0.94)
    assert nbo.gated_pick([a, b], native) is b
    c, d = _pos(24, True, 0.93), _pos(22, True, 0.93)
    assert nbo.gated_pick([c, d], native) is d


def test_only_preserved_sites_are_picked():
    assert nbo.gated_pick([_pos(21, True, 0.9, verdict="cleavage_moved")], 0.95) is None


def test_one_softening_constant():
    """No module restates the bound as a number; every caller reads ``SOFTENING_BOUND``."""
    import ast

    for name in ("nterm_decide.py", "nterm_overnight.py", "nterm_site_picks.py"):
        tree = ast.parse((NTERM / name).read_text())
        literals = [
            n.lineno for n in ast.walk(tree) if isinstance(n, ast.Constant) and n.value == 0.07
        ]
        assert not literals, f"{name} restates the softening bound at lines {literals}"
