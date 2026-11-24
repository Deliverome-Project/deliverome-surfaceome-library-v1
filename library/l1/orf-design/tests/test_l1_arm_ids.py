"""Arm ids for every built insert: one per insert, permanent, from the per-arm table.

Synthetic: no S3, no committed tables beyond a temporary per-arm table. Implemented by the
Claude Code session l1-reviewer-package (2026-11-15).
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))

import l1_arm_ids as arms  # noqa: E402

from deliverome_analysis import library_design as ld  # noqa: E402

COLS = ["arm_id", "construct_key", "design_id", "arm_role", "arm_kind", "decision"]


@pytest.fixture
def decisions(tmp_path) -> dict:
    rows = [
        ["l1-p0325-a1", "l1-p0325:internal49", "l1-p0325", "loop_best", "loop", "build"],
        ["l1-p0325-a2", "l1-p0325:snorkel", "l1-p0325", "snorkel_companion", "snorkel", "build"],
        ["l1-p0001-a1", "l1-p0001:internal389", "l1-p0001", "loop_best", "loop", "build"],
        ["l1-p0001-a2", "l1-p0001:internal1012", "l1-p0001", "loop_2", "loop", "build"],
        ["l1-p1321-a1", "", "l1-p1321", "", "", "hold"],
        ["l1-p0009-a3", "l1-p0009:internal77", "l1-p0009", "loop_3", "loop", "hold"],
    ]
    path = tmp_path / "arms.tsv"
    with path.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(COLS)
        w.writerows(rows)
    return arms.load_decisions(path)


def test_single_arm_rungs_are_the_members_first_arm(decisions):
    for rung in ("nterm", "multipass", "type3", "cterm"):
        assert arms.arm_id_for(rung, {"placement": "C-terminal"}, "l1-p0042", decisions) == (
            "l1-p0042-a1"
        )


def test_internal_and_snorkel_arms_come_from_the_per_arm_table(decisions):
    row = {"construct_id": "l1-p0001:internal1012"}
    assert arms.arm_id_for("internal", row, "l1-p0001", decisions) == "l1-p0001-a2"
    snorkel = {"placement": "snorkel (TM module, then cassette)"}
    assert arms.arm_id_for("cterm", snorkel, "l1-p0325", decisions) == "l1-p0325-a2"


def test_an_arm_table_column_wins_but_must_belong_to_the_member(decisions):
    assert arms.arm_id_for("nterm", {"arm_id": "l1-p0042-a3"}, "l1-p0042", decisions) == (
        "l1-p0042-a3"
    )
    with pytest.raises(ValueError, match="not an arm of"):
        arms.arm_id_for("nterm", {"arm_id": "l1-p0043-a1"}, "l1-p0042", decisions)


def test_an_internal_or_snorkel_insert_the_table_does_not_build_fails(decisions):
    with pytest.raises(ValueError, match="no row"):
        arms.arm_id_for("internal", {"construct_id": "l1-p0002:internal5"}, "l1-p0002", decisions)
    snorkel = {"placement": "snorkel (TM module, then cassette)"}
    with pytest.raises(ValueError, match="no row"):  # held: a snorkel is never built alone
        arms.arm_id_for("cterm", snorkel, "l1-p1321", decisions)
    with pytest.raises(ValueError, match="decides 'hold'"):
        arms.arm_id_for("internal", {"construct_id": "l1-p0009:internal77"}, "l1-p0009", decisions)


def test_two_inserts_cannot_share_an_arm_id():
    arms.check_unique(["l1-p0001-a1", "l1-p0001-a2"])
    with pytest.raises(ValueError, match="more than one insert"):
        arms.check_unique(["l1-p0001-a1", "l1-p0002-a1", "l1-p0001-a1"])


def test_id_helpers_round_trip_and_refuse_member_ids():
    assert ld.arm_id("l1-p0325", 2) == "l1-p0325-a2"
    assert ld.design_of_arm("l1-p0325-a2") == "l1-p0325"
    assert ld.insert_id(7) == "l1-i0007"
    for bad in ("l1-m0325", "l1-p0325-a0", "l1-p0325:snorkel"):
        with pytest.raises(ValueError):
            ld.design_of_arm(bad)
    with pytest.raises(ValueError):
        ld.arm_id("l1-m0001", 1)
