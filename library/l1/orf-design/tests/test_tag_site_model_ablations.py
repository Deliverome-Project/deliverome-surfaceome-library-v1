"""Unit tests for tag-site ablation bug fixes (contiguous smoothing and pool z-scoring)."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

S4E9 = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(S4E9))

import tag_site_model_ablations as A  # noqa: E402
import tag_site_negative_calibration as N  # noqa: E402


def test_smooth_contiguous_does_not_bleed_across_tm_gap() -> None:
    """Smoothing must not average residues from two separate extracellular loops."""
    df = pd.DataFrame(
        {
            "uniprot_acc": ["P1"] * 6,
            "residue": [10, 11, 12, 100, 101, 102],
        }
    )
    scores = np.array([0.0, 0.0, 0.0, 1.0, 1.0, 1.0])
    out = A.smooth_contiguous(df, scores, half=2)
    np.testing.assert_allclose(out[:3], [0.0, 0.0, 0.0])
    np.testing.assert_allclose(out[3:], [1.0, 1.0, 1.0])


def test_zscore_from_pool_matches_full_distribution_not_subsample() -> None:
    """Z-scoring uses the full extracellular pool mean/std rather than a biased subsample."""
    pool = pd.DataFrame({"uniprot_acc": ["P1"] * 4, "x": [0.0, 0.0, 2.0, 2.0]})
    sub = pd.DataFrame({"uniprot_acc": ["P1"] * 3, "x": [0.0, 0.0, 2.0]})
    out = A.zscore_from_pool(sub, pool, ["x"])
    expected_mu = pool["x"].mean()
    expected_sd = pool["x"].std()
    np.testing.assert_allclose(
        out["x"].to_numpy(), (sub["x"].to_numpy() - expected_mu) / expected_sd
    )


def test_add_region_geometry_bounds() -> None:
    """rel_pos_in_region is 0.0 at loop edges and 0.5 at the midpoint."""
    df = pd.DataFrame(
        {
            "uniprot_acc": ["P1"] * 5,
            "region_index": [1] * 5,
            "residue": [20, 21, 22, 23, 24],
        }
    )
    out = A.add_region_geometry(df)
    np.testing.assert_allclose(out["region_len"].to_numpy(), [5.0] * 5)
    np.testing.assert_allclose(out["rel_pos_in_region"].to_numpy(), [0.0, 0.25, 0.5, 0.25, 0.0])


def test_hedged_picks_forces_clean_pick2_when_pick1_is_soft_flagged() -> None:
    """When Pick 1 lands on a soft negative, Pick 2 must be clean even if a soft site scores higher."""
    res = np.array([10, 30, 50, 70])
    scores = np.array([0.9, 0.8, 0.7, 0.6])
    soft = np.array([True, True, False, False])
    picks = N.hedged_picks(res, scores, soft, k=3, sep=11)
    assert picks[0] == 0  # best overall (soft)
    assert picks[1] == 2  # forced clean Pick 2
    assert not soft[picks[1]]
