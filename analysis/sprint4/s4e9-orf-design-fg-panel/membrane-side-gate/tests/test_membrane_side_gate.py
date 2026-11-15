"""Seeded cases for the L1 membrane-side gate, and consistency of its committed outputs.

Members are named by design id (``l1-pNNNN``), which survives renumbering. The two seeded picks are the cases that motivated the gate: ANO1 insert-after-698 (the
S5E13 v4 rank-1 pick for member l1-p0110; cytoplasmic in UniProt and in every mouse
TMEM16A structure) and XKR5 insert-after-62 (v4 rank-1 for l1-p2348). Both are placed
from the member's own AlphaFold DB model (C-alpha/C-beta copies in ``fixtures/``) with
the thresholds ``calibrate.py`` derived, so this runs offline.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from deliverome_analysis import membrane_side as ms

HERE = Path(__file__).resolve().parent
GATE = HERE.parent
FIX = HERE / "fixtures"
SEEDED = json.loads((FIX / "seeded_members.json").read_text())
CAL = json.loads((GATE / "calibration.json").read_text())
TH = ms.GateThresholds(**CAL["thresholds"])


def _junction(z: np.ndarray, n: int) -> float:
    return float(np.nanmean(z[n - 1 : n + 1]))


def _place(cid: str):
    m = SEEDED[cid]
    model = ms.read_model(FIX / f"AF-{m['uniprot_acc']}-F1-model_v6.ca_cb.pdb")
    rmap, match = ms.map_to_model(m["protein"], model.sequence)
    assert match == "identical"
    p = ms.place(model, ms.topology_on_model(m["topology"], rmap, len(model.sequence)))
    ok = ms.orientation_supported(p.topo_for, p.topo_against, p.kr_out, p.kr_in, TH)
    return m, model, p, ok


@pytest.mark.parametrize(("cid", "site"), [("l1-p0110", 698), ("l1-p2348", 62)])
def test_seeded_cytoplasmic_picks_fail(cid, site):
    m, model, p, ok = _place(cid)
    # DeepTMHMM v1 calls the junction extracellular: the error the gate exists to catch.
    assert m["topology"][site - 1] == "O"
    assert ok
    z = _junction(p.z, site)
    assert z < -TH.core
    assert ms.site_verdict(z, float(model.plddt[site - 1]), ok, TH) == "fail_cytoplasmic"


def test_ano1_replacement_candidates_pass():
    """ANO1's extracellular loops 541-565 and 806-891 are outside on the same placement."""
    _, model, p, ok = _place("l1-p0110")
    for site in (554, 823, 850, 874):
        z = _junction(p.z, site)
        assert ms.site_verdict(z, float(model.plddt[site - 1]), ok, TH) == "pass", site


def test_ano1_misassigned_segment_spans_both_sides():
    """v1's 628-731 'O' segment holds a real extracellular stretch and a cytoplasmic loop."""
    m, _, p, _ = _place("l1-p0110")
    seg = [i for i in range(627, 731) if m["topology"][i] == "O"]
    z = p.z[seg]
    assert (z >= TH.outside).any() and (z < -TH.core).any()


# --- committed outputs --------------------------------------------------------------

GATE_TSV = GATE / "l1_internal_pick_gate.tsv"
CHANGES_TSV = GATE / "l1_internal_pick_changes.tsv"


def test_committed_gate_table_has_seeded_failures():
    gate = pd.read_csv(GATE_TSV, sep="\t")
    for cid, site in (("l1-p0110", 698), ("l1-p2348", 62)):
        row = gate[
            (gate["design_id"] == cid)
            & (gate["source"] == "s5e13_v4")
            & (gate["insert_after_residue"] == site)
        ]
        assert len(row) == 1
        assert row["verdict"].iloc[0] == "fail_cytoplasmic"
    changes = pd.read_csv(CHANGES_TSV, sep="\t")
    assert {"l1-p0110", "l1-p2348"} <= set(changes["design_id"])


def test_every_internal_member_is_reported():
    summary = json.loads((GATE / "gate_summary.json").read_text())
    members = pd.read_csv(GATE / "l1_internal_member_gate.tsv", sep="\t")
    internal = members[members["placement"] == "internal"]
    assert len(internal) == summary["internal_members"]
    # No member without a model is passed: each is reported, with why.
    unscored = internal[internal["status"] != "placed"]
    assert set(unscored["candidate_id"]) == set(summary["internal_members_unscorable"])
    gate = pd.read_csv(GATE_TSV, sep="\t")
    assert not (
        gate["candidate_id"].isin(unscored["candidate_id"]) & gate["verdict"].eq("pass")
    ).any()


def test_thresholds_rest_on_their_evidence():
    ev = CAL["evidence"]
    assert TH.core == pytest.approx(ev["tmdet_half_thickness_median"], abs=0.05)
    assert TH.outside <= ev["published_displayed_multipass_z_percentiles"]["0"]


MASK_TSV = GATE / "l1_extracellular_mask_residues.tsv"


def test_residue_mask_covers_every_internal_member_and_matches_seeded_calls():
    mask = pd.read_csv(MASK_TSV, sep="\t", keep_default_na=False).set_index("design_id")
    members = pd.read_csv(GATE / "l1_internal_member_gate.tsv", sep="\t")
    assert set(mask.index) == set(members.loc[members["placement"] == "internal", "design_id"])
    usable = mask[mask["usable"] == 1]
    assert (usable["mask"].str.len() == usable["length_aa"]).all()
    # A member the gate cannot place, or whose face is unconfirmed, offers nothing to rank.
    assert (mask.loc[mask["status"] != "placed", "usable"] == 0).all()
    ano1 = mask.loc["l1-p0110", "mask"]
    assert ano1[697:699] != "EE"  # insert after 698: cytoplasmic
    assert ano1[822:824] == "EE"  # insert after 823: extracellular loop 806-891
    assert mask.loc["l1-p2348", "mask"][61:63] != "EE"  # XKR5 insert after 62


def test_junction_mask_reproduces_every_gate_verdict():
    """A ranker filtering on junction_mask == 'E' keeps exactly the sites this gate passes."""
    mask = pd.read_csv(MASK_TSV, sep="\t", keep_default_na=False).set_index("candidate_id")
    gate = pd.read_csv(GATE_TSV, sep="\t")
    scored = gate[gate["candidate_id"].isin(mask.index[mask["usable"] == 1])]
    scored = scored[
        scored["verdict"].isin(["pass", "review_interface", "fail_membrane", "fail_cytoplasmic"])
    ]
    for r in scored.itertuples():
        ch = mask.loc[r.candidate_id, "junction_mask"][int(r.insert_after_residue) - 1]
        if r.v1_state[0] == "O":
            assert (ch == "E") == (r.verdict == "pass"), (r.candidate_id, r.insert_after_residue)


def test_manual_partner_interface_exclusion_is_applied_and_attributed():
    """The beta-docking loop of ATP1A1 is excluded by hand, and the row says who decided it."""
    gate = pd.read_csv(GATE_TSV, sep="\t")
    row = gate[(gate["design_id"] == "l1-p0152") & (gate["insert_after_residue"] == 893)]
    assert len(row) == 1
    r = row.iloc[0]
    assert r["verdict"] == "fail_partner_interface"
    assert r["membrane_verdict"] == "pass"  # the membrane-side call is kept alongside
    assert "Becca, 2026-11-11" in r["decided_by"]
    assert "SYGQ" in r["reason"]
    mask = pd.read_csv(MASK_TSV, sep="\t", keep_default_na=False).set_index("design_id")
    assert set(mask.loc["l1-p0152", "junction_mask"][873:919]) == {"x"}
    # Sites outside the excluded loop are untouched: the M1-M2 site still passes.
    keep = gate[(gate["design_id"] == "l1-p0152") & (gate["insert_after_residue"] == 123)]
    assert keep["verdict"].iloc[0] == "pass"


def test_topology_decisions_for_mfsd2a_and_disp2():
    """MFSD2A is placed on its inverted v1 sides and usable; DISP2 keeps v1 and is released."""
    mask = pd.read_csv(MASK_TSV, sep="\t", keep_default_na=False).set_index("design_id")
    for d in ("l1-p1213", "l1-p2466"):
        r = mask.loc[d]
        assert r["topology_tool"] == "deeptmhmm_v1_inverted_by_override"
        assert r["usable"] == 1
    # The canonical's glycosylated loop (N230, N240) is now eligible.
    assert set(mask.loc["l1-p2466", "junction_mask"][229:240]) == {"E"}
    disp2 = mask.loc["l1-p0513"]
    assert disp2["topology_tool"] == "deeptmhmm_v1" and disp2["usable"] == 1
    over = pd.read_csv(GATE / "topology_overrides.tsv", sep="\t")
    assert over["decided_by"].str.contains("delegated by Becca 2026-11-12").all()


def test_uniprot_inversion_is_keyed_by_design_id():
    """Member ids are renumbered between panels; the calibration must not depend on them."""
    assert CAL["uniprot_inverted_design_ids"] == ["l1-p1995"]
    assert "uniprot_inverted" not in CAL


def test_snapshot_check_found_no_changed_calls():
    for run in ("a", "b"):
        d = json.loads(
            (GATE / "snapshot-check-2026-11-14" / f"diff_{run}_summary.json").read_text()
        )
        assert d["members_with_any_change"] == 0 and d["site_verdict_changes"] == 0


def test_vendored_s5e13_inputs_match_their_hashes():
    import hashlib

    v = GATE / "s5e13-c7be01a4"
    src = json.loads((v / "SOURCE.json").read_text())
    assert src["source_commit"].startswith("c7be01a4")
    for name, digest in src["sha256"].items():
        assert hashlib.sha256((v / name).read_bytes()).hexdigest() == digest


def test_no_ordered_site_fails_the_gate():
    t = pd.read_csv(GATE / "ordered-sites-2026-11-14" / "gate_on_ordered_sites.tsv", sep="\t")
    assert not t["verdict"].str.startswith("fail").any()
    assert not (t["verdict"] == "review_interface").any()
    assert not (t["verdict"] == "not_evaluated").any()
    atp = t[t["gene_symbol"].isin(["ATP1A1", "ATP1A2", "ATP1A3", "ATP1A4", "ATP12A"])]
    assert len(atp) == 5 and (atp["verdict"] == "pass").all()


def test_every_non_passing_ordered_site_has_a_disposition():
    d = GATE / "ordered-sites-2026-11-14"
    t = pd.read_csv(d / "gate_on_ordered_sites.tsv", sep="\t")
    disp = pd.read_csv(d / "dispositions.tsv", sep="\t")
    open_sites = set(t.loc[~t["verdict"].isin(["pass", "pass_snorkel"]), "design_id"])
    assert open_sites == set(disp["design_id"])
    assert disp["decided_by"].str.contains("delegated by Becca 2026-11-12").all()


def test_redesigned_ordered_sites_pass():
    """PR #331's seven redesigned members, gated on that PR's decisions file."""
    d = GATE / "ordered-sites-2026-11-14-pr331"
    t = pd.read_csv(d / "gate_on_ordered_sites.tsv", sep="\t").set_index("design_id")
    new_sites = {
        "l1-p0111": 868, "l1-p0178": 268, "l1-p0433": 53, "l1-p1773": 651,
        "l1-p1884": 59, "l1-p1972": 545, "l1-p2181": 55, "l1-p0153": 311, "l1-p1032": 212,
    }  # fmt: skip
    for design_id, site in new_sites.items():
        assert t.loc[design_id, "insert_after_residue"] == site
        assert t.loc[design_id, "verdict"] == "pass"
    disp = pd.read_csv(GATE / "ordered-sites-2026-11-14" / "dispositions.tsv", sep="\t")
    assert set(t.index[~t["verdict"].isin(["pass", "pass_snorkel"])]) == set(disp["design_id"])


def test_rescue_round_gate():
    """Rescue round (2026-11-14b): single-pass sites are never failed on model geometry."""
    t = pd.read_csv(GATE / "rescue-round-2026-11-14b" / "rescue_round_gate.tsv", sep="\t")
    assert len(t) == 74 and t["design_id"].nunique() == 39
    assert set(t["verdict"]) <= {"pass", "review_interface", "review_orientation",
                                 "review_single_pass_geometry", "not_placeable"}  # fmt: skip
    assert not t["verdict"].str.startswith("fail").any()
    single = t[t["n_tm"] < 2]
    assert set(single["basis"]) <= {"single_pass_tmdet", "not_placeable"}
    assert (single.loc[single["basis"] == "single_pass_tmdet", "placement_source"] == "tmdet").all()
    assert (t.loc[t["verdict"] == "not_placeable", "reason"].str.len() > 0).all()
    # VTCN1: no TmAlphaFold entry and one predicted TM, so not placeable either way.
    assert (t.loc[t["gene"] == "VTCN1", "verdict"] == "not_placeable").all()
