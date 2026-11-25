"""Unit tests for swap-pool leader classification, order-plan ratios, and worked example maps."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest
from Bio import SeqIO

S4E9 = Path(__file__).resolve().parents[1]
REVIEWER = S4E9 / "reviewer-package"
for _p in (S4E9, REVIEWER):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

import diversify_arm_cassettes as dac  # noqa: E402
import export_review_repo as err  # noqa: E402
import l1_swap_pool_model as swap_model  # noqa: E402
import order_plan  # noqa: E402
import reserve_l1_block as rlb  # noqa: E402
from barcode_prefix_screen import PREFIX_BASES  # noqa: E402
from build_constructs import CASSETTE_AA  # noqa: E402


def _make_construct(prefix_len: int, body: str, *, has_cleavable_leader: bool) -> dict:
    prefix = "M" + "A" * (prefix_len - 1)
    protein = prefix + CASSETTE_AA + body + "*"
    return {
        "protein": protein,
        "cas_aa": prefix_len,
        "has_cleavable_leader": has_cleavable_leader,
    }


def _recombined_at_cassette(a: dict, b: dict) -> str:
    return a["protein"][: a["cas_aa"]] + b["protein"][b["cas_aa"] :]


class TestLeaderSwapClassification:
    def test_sp_to_sp_is_leader_swap(self) -> None:
        a = _make_construct(22, "LLLLLLLL", has_cleavable_leader=True)
        b = _make_construct(30, "VVVVVVVV", has_cleavable_leader=True)
        prot = _recombined_at_cassette(a, b)
        assert swap_model.leader_swap(prot, a, b) is True

    @pytest.mark.parametrize(
        ("a_leader", "b_leader", "a_prefix", "b_prefix"),
        [
            # Internal-loop construct (cas_aa <= 80) as barcode template b
            (True, False, 22, 40),
            # Internal-loop construct (cas_aa <= 80) donating uncleaved N-tail to SP construct b
            (False, True, 40, 22),
            # Leaderless type III / CRHR2 (cas_aa = 1) deleting SP of construct b
            (False, True, 1, 25),
            # SP construct adding leader onto leaderless type III / CRHR2
            (True, False, 25, 1),
            # HPN-like stalk insertion (cas_aa = 43, cterm rung, no cleavable leader)
            (True, False, 22, 43),
            # Both non-leader with short cas_aa <= 80
            (False, False, 26, 43),
            # Cleavable-leader construct with long propeptide (> LEADER_AA)
            (True, True, 120, 25),
            (True, True, 25, 120),
        ],
    )
    def test_non_cleavable_or_long_leader_is_not_leader_swap(
        self, a_leader: bool, b_leader: bool, a_prefix: int, b_prefix: int
    ) -> None:
        a = _make_construct(a_prefix, "LLLLLLLL", has_cleavable_leader=a_leader)
        b = _make_construct(b_prefix, "VVVVVVVV", has_cleavable_leader=b_leader)
        prot = _recombined_at_cassette(a, b)
        assert swap_model.leader_swap(prot, a, b) is False

    def test_cleavable_leader_keys_exclude_non_leader_constructs(self) -> None:
        cleavable = swap_model._cleavable_leader_keys()
        # 1,322 nterm + 567 multipass igk_tag + 9 type3 igk_tag = 1,898
        assert len(cleavable) == 1322 + 567 + 9
        # CRHR2 (leaderless multipass), HPN (cterm stalk), internal loops, type3 tag_alone excluded
        assert ("multipass", "l1-p0475") not in cleavable  # CRHR2
        assert ("cterm", "l1-p0864") not in cleavable  # HPN
        assert ("type3", "l1-p0181") not in cleavable  # BMPR1B (tag_alone)
        assert not any(rung in ("internal", "cterm") for rung, _ in cleavable)


class TestOrderPlanRatios:
    def test_factor_is_wired_into_sensitivity_and_summary(self) -> None:
        assert order_plan.FACTOR in order_plan.SENSITIVITY
        rows, summary = order_plan.plan(order_plan.inserts())
        assert rows
        key = f"pool_ug_total_per_pool_x{order_plan.FACTOR}"
        for scenario in summary["concentration_scenarios"].values():
            assert key in scenario

    def test_molar_vs_mass_volume_ratios_and_extra_dna(self) -> None:
        rows, summary = order_plan.plan(order_plan.inserts())
        shortest, median, longest = rows[0], rows[len(rows) // 2], rows[-1]
        ref = (
            summary["reference_insert"]["insert_nt"],
            summary["reference_insert"]["vector_bp"],
        )

        for f in order_plan.SENSITIVITY:
            w = summary["weighting"][f"x{f}"]
            pmol_l_vs_s = f ** ((longest["insert_nt"] - shortest["insert_nt"]) / 1000)
            pmol_l_vs_m = f ** ((longest["insert_nt"] - median["insert_nt"]) / 1000)
            assert w["longest_vs_shortest_pmol"] == pytest.approx(round(pmol_l_vs_s, 1))
            assert w["longest_vs_median_pmol"] == pytest.approx(round(pmol_l_vs_m, 1))

            # DNA mass/volume ratio multiplies molar ratio by plasmid-length (MW) ratio
            ng_s = order_plan.pool_ng(shortest["insert_nt"], shortest["vector_bp"], f, ref)
            ng_m = order_plan.pool_ng(median["insert_nt"], median["vector_bp"], f, ref)
            ng_l = order_plan.pool_ng(longest["insert_nt"], longest["vector_bp"], f, ref)
            assert ng_l / ng_s == pytest.approx(
                pmol_l_vs_s * (longest["vector_bp"] / shortest["vector_bp"])
            )
            assert ng_l / ng_m == pytest.approx(
                pmol_l_vs_m * (longest["vector_bp"] / median["vector_bp"])
            )
            assert ng_l / ng_s > pmol_l_vs_s
            assert ng_l / ng_m > pmol_l_vs_m

            # Ordered extra DNA scales assumed delivered extra DNA by ORDERED_UG / ASSUMED_UG (4/3)
            extra_assumed = sum(r[f"extra_ug_assumed_x{f}"] for r in rows)
            extra_ordered = sum(r[f"order_ug_x{f}"] - order_plan.ORDERED_UG for r in rows)
            assert w["extra_ug_assumed_total"] == pytest.approx(round(extra_assumed, 1))
            assert extra_ordered == pytest.approx(
                extra_assumed * order_plan.ORDERED_UG / order_plan.ASSUMED_UG, abs=0.2
            )

    def test_subpool_strategies_reduce_total_skew_monotonically(self) -> None:
        rows, _ = order_plan.plan(order_plan.inserts())
        sub = order_plan.subpool_strategies(rows)
        for f in order_plan.SENSITIVITY:
            s1 = sub["k1_per_plate"][f"x{f}"]["total_max_min_skew_fold"]
            s2 = sub["k2_per_plate"][f"x{f}"]["total_max_min_skew_fold"]
            s4 = sub["k4_per_plate"][f"x{f}"]["total_max_min_skew_fold"]
            assert s1 > s2 > s4 > 1.0


def test_axl_ires_puro_gb_puror_translation_starts_with_m() -> None:
    gb_path = S4E9 / "l1_example_AXL_ires_puro.gb"
    rec = SeqIO.read(gb_path, "genbank")
    puro = [f for f in rec.features if f.type == "CDS" and f.qualifiers.get("label") == ["PuroR"]]
    assert len(puro) == 1
    feat = puro[0]
    dna = feat.extract(rec.seq)
    assert str(dna).upper().startswith("ATG")
    expected_aa = str(dna.translate())
    assert expected_aa.startswith("MTEYKPT")
    assert feat.qualifiers["translation"] == [expected_aa]


def test_diversify_arm_cassettes_held_arms_and_clean_path_rewrite() -> None:
    assert dac.held_arms()
    rel = "library/l1/orf-design/diversify_arm_cassettes.py"
    rewritten = err.apply_clean_paths({rel: (S4E9 / "diversify_arm_cassettes.py").read_bytes()})
    text = rewritten["library/l1/orf-design/diversify_arm_cassettes.py"].decode("utf-8")
    assert ("sprint5/" + "s5e13-tag-site-benchmarking") not in text
    assert '"l1/tag-site-benchmarking" / "l1_internal_arm_decisions.tsv"' in text


def test_reserve_l1_block_pins_10nt_prefix_while_default_is_12nt() -> None:
    assert rlb.L1_PREFIX_BASES == 10
    assert PREFIX_BASES == 12
