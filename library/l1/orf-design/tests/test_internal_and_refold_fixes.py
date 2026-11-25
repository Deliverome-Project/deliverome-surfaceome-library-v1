"""Regression tests for internal-site ranker, membrane-side mask, GPI buffer, and refold fixes."""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[1]
REPO = HERE.parents[2]
S5E13 = REPO / "library/l1/tag-site-benchmarking"
TAG_REFOLD = HERE / "tag_refold"
GATE_DIR = HERE / "membrane-side-gate"
NTERM_DIR = HERE / "reviewer-package/nterm"

for p in (HERE, S5E13, TAG_REFOLD, GATE_DIR, NTERM_DIR):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))


def test_masked_regions_does_not_split_on_review_band_dips():
    import membrane_mask as MM

    # Single extracellular loop with short 'r' review-band dips (e.g. TRPV3/AQP4/SLC7A11)
    mask = "....EEEEErrrrEEEE....EEEcEEE...EEEmEEE..."
    jmask = mask
    elig = [i + 1 for i, c in enumerate(jmask) if c == "E"]
    regs = MM.masked_regions(elig, (mask, jmask))
    reg_by_res = dict(zip(elig, regs, strict=True))

    # First loop bridged across 'rrrr' stays one region (region 0)
    assert reg_by_res[5] == reg_by_res[14] == 0
    # Transmembrane '.' splits to region 1; cytoplasmic 'c' splits to region 2
    assert reg_by_res[22] == 1
    assert reg_by_res[26] == 2
    # Membrane core 'm' splits to region 4
    assert reg_by_res[32] == 3
    assert reg_by_res[36] == 4


def test_membrane_mask_loads_r8_locally():
    import membrane_mask as MM

    assert MM.LOCAL.is_file()
    assert "2026-11-12-r8" in MM.MASK_KEY
    df = MM.load()
    assert len(df) > 500
    assert (df["design_id"] == "l1-p1961").any()


def test_v8_prodomino_fails_closed_on_unusable_mask_and_recomputes_regions():
    import design_l1_internal_v8_prodomino as V8

    mm = pd.DataFrame(
        [
            {
                "sequence_sha256": "sha_unusable",
                "usable": 0,
                "orientation_supported": False,
                "mask": "EEEEEEEE",
                "junction_mask": "EEEEEEEE",
            },
            {
                "sequence_sha256": "sha_usable",
                "usable": 1,
                "orientation_supported": True,
                "mask": "EErrEE..EE",
                "junction_mask": "EErrEE..EE",
            },
        ]
    )
    pool = pd.DataFrame(
        [{"sequence_sha256": "sha_unusable", "residue": r, "region_index": 0} for r in range(1, 9)]
        + [
            {"sequence_sha256": "sha_usable", "residue": r, "region_index": 99}
            for r in range(1, 11)
        ]
    )
    out = V8.apply_membrane_mask(pool, mm)

    # Unusable protein (usable=0 / orientation_supported=False) fails CLOSED
    unusable_rows = out[out["sequence_sha256"] == "sha_unusable"]
    assert not unusable_rows["eligible"].any()

    # Usable protein keeps only 'E' junctions and merges across 'rr' into region 0, splitting at '..'
    usable_elig = out[(out["sequence_sha256"] == "sha_usable") & out["eligible"]]
    reg_map = dict(zip(usable_elig["residue"], usable_elig["region_index"], strict=True))
    assert set(reg_map) == {1, 2, 5, 6, 9, 10}
    assert reg_map[1] == reg_map[5] == 0
    assert reg_map[9] == 1


def test_v8_prodomino_mutually_exclusive_eb_penalties_and_seam_windowing():
    import design_l1_internal_v8_prodomino as V8

    pool = pd.DataFrame(
        {
            "uniprot_acc": ["P00001"] * 12,
            "residue": list(range(1, 13)),
            # Residue 1: q90 (>=459), Residue 2: q75 only (400), others uncrowded
            "protrusion_14": [500.0, 400.0] + [100.0] * 10,
            # Residue 6 has d_other_domain=5.0 (<=6A), residue 12 has 7.0 (<=8A)
            "d_other_domain": [
                20.0,
                20.0,
                20.0,
                20.0,
                20.0,
                5.0,
                20.0,
                20.0,
                20.0,
                20.0,
                20.0,
                7.0,
            ],
        }
    )
    zeros = np.zeros(len(pool), dtype=bool)
    out, pen = V8.compute_soft_penalties(pool, zeros, zeros)

    # Crowding mutual exclusivity: residue 1 gets ONLY q90 (1.8071), NOT q75 + q90
    assert bool(out.loc[0, "crowded_q90"]) is True
    assert bool(out.loc[0, "crowded_q75"]) is False
    assert np.isclose(pen[0], 1.8071)

    # Residue 2 gets ONLY q75 (1.7569)
    assert bool(out.loc[1, "crowded_q75"]) is True
    assert bool(out.loc[1, "crowded_q90"]) is False
    assert np.isclose(pen[1], 1.7569)

    # Domain seam n-2..n+3 windowing around residue 6 (seam_6A):
    # A junction after n hits residue 6 iff n-2 <= 6 <= n+3  <=>  3 <= n <= 8.
    seam6_res = set(out.loc[out["seam_6A"], "residue"])
    assert seam6_res == {3, 4, 5, 6, 7, 8}
    # seam_8A excludes seam_6A (mutually exclusive): around residue 12 -> 9..12 (since n<=12)
    seam8_res = set(out.loc[out["seam_8A"], "residue"])
    assert seam8_res == {9, 10, 11, 12}
    assert not (out["seam_6A"] & out["seam_8A"]).any()
    assert np.isclose(pen[5], 0.5497)
    assert np.isclose(pen[11], 0.4774)


def test_v8_prodomino_pick_arms_snorkel_gating_and_deterministic_ties():
    import design_l1_internal_v8_prodomino as V8

    g = pd.DataFrame(
        {
            "residue": [10, 11, 12, 13, 50, 51, 52, 53, 90, 91, 92, 93],
            "region_index": [0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2],
            "eligible": [True] * 12,
            "hard_veto": [False] * 12,
            "s_rank": [0.8] * 12,  # exact ties -> must pick lowest residue per region
            "s_abs": [0.2] * 12,  # below pub_q25=0.5
            "plddt": [70.0] * 12,
            "prodomino": [0.5] * 12,
            "protrusion_14": [200.0] * 12,
            "d_other_domain": [20.0] * 12,
        }
    )
    # Cytoplasmic C-term -> 2 loop arms + snorkel companion
    res_cyto = V8.pick_arms(g, pub_med=0.7, pub_q25=0.5, cterm_cytoplasmic=True)
    assert res_cyto["snorkel_companion"] is True
    assert res_cyto["n_loop_arms"] == 2
    assert [p["residue"] for p in res_cyto["picks"]] == [10, 50]

    # Non-cytoplasmic C-term -> no snorkel companion, falls back to 3 loops
    res_ecto = V8.pick_arms(g, pub_med=0.7, pub_q25=0.5, cterm_cytoplasmic=False)
    assert res_ecto["snorkel_companion"] is False
    assert res_ecto["n_loop_arms"] == 3
    assert [p["residue"] for p in res_ecto["picks"]] == [10, 50, 90]


def test_rescue_round_enforces_gpi_omega_upstream_clearance(monkeypatch):
    import rescue_round as RR

    seq = "M" + "A" * 99  # 100 aa protein, omega at 80
    fake_uniprot = {
        "sequence": {"value": seq},
        "features": [
            {"type": "Signal", "location": {"start": {"value": 1}, "end": {"value": 10}}},
            {
                "type": "Lipidation",
                "description": "GPI-anchor amidated serine",
                "location": {"start": {"value": 80}, "end": {"value": 80}},
            },
            {"type": "Chain", "location": {"start": {"value": 11}, "end": {"value": 80}}},
        ],
    }
    monkeypatch.setattr(RR, "uniprot", lambda acc: fake_uniprot)
    m = pd.Series(
        {"identical_acc": "P99999", "gene_symbol": "FAKEGPI", "candidate_set": "library_reference"}
    )
    gpi = pd.DataFrame(
        columns=[
            "gene_symbol",
            "gpi_anchored",
            "omega_site_reference",
            "omega_form",
            "omega_site",
            "omega_aa",
        ]
    )
    res = RR.span(m, seq, {}, gpi)
    assert res["omega_member"] == 80
    assert RR.GPI_OMEGA_UPSTREAM_CLEARANCE == 12
    assert max(res["allowed"]) == 80 - 12
    # Positions omega-11 .. omega (69..80) must be excluded
    assert not (res["allowed"] & set(range(69, 81)))


def test_score_tagged_folds_discontiguous_helix_and_run_local_plddt():
    import score_tagged_folds as STF

    # Cassette indices 10..60 (len 51), ALFA at 25..40 (len 15)
    cas = np.arange(10, 61)
    a0 = 25
    alfa_len = 15

    # Case 1: 6 helical residues right before ALFA (19..24) + 5 right after ALFA (40..44)
    # Previously falsely concatenated to 11; now evaluated separately -> max_run=6, conf_run=0
    hel = np.zeros(80, dtype=bool)
    plddt = np.full(80, 90.0)
    hel[19:25] = True
    hel[40:45] = True
    run, conf_run = STF.cassette_non_alfa_helix(hel, plddt, cas, a0, alfa_len=alfa_len)
    assert run == 6
    assert conf_run == 0

    # Case 2: 10-aa low-pLDDT run (55.0) in pre-ALFA + 10-aa? No, 8-aa high-pLDDT (95.0) in post-ALFA
    # Global mean over all helical non_alfa would be >70, but the 10-aa run itself has mean 55 < 70
    hel2 = np.zeros(80, dtype=bool)
    plddt2 = np.full(80, 50.0)
    hel2[12:22] = True
    plddt2[12:22] = 55.0
    hel2[42:50] = True
    plddt2[42:50] = 95.0
    run2, conf_run2 = STF.cassette_non_alfa_helix(hel2, plddt2, cas, a0, alfa_len=alfa_len)
    assert run2 == 10
    assert conf_run2 == 0

    # Case 3: Genuine 10-aa run at pLDDT 80 in pre-ALFA + low-pLDDT 5-aa run in post-ALFA -> flagged
    hel3 = np.zeros(80, dtype=bool)
    plddt3 = np.full(80, 40.0)
    hel3[12:22] = True
    plddt3[12:22] = 80.0
    hel3[42:47] = True
    plddt3[42:47] = 30.0
    run3, conf_run3 = STF.cassette_non_alfa_helix(hel3, plddt3, cas, a0, alfa_len=alfa_len)
    assert run3 == 10
    assert conf_run3 == 10


def test_compare_internal_policies_window_and_tie_safe_ranking():
    import compare_internal_policies as CIP
    import ranker_hard_negative_layer as L

    # Contact at residue 23: insertion after n=20 has n-2..n+3 = 18..23 -> HIT!
    # Insertion after n=19 has 17..22 -> MISS! Insertion after n=26 has 24..29 -> MISS!
    ev = pd.DataFrame(
        {
            "uniprot_acc": ["P1"] * 10,
            "residue": list(range(18, 28)),
        }
    )
    contact = np.array([r == 23 for r in ev["residue"]], dtype=bool)
    veto = L.near(ev, pd.Series(contact), CIP.VETO_LO, hi=CIP.VETO_HI)
    vetoed_res = set(ev.loc[veto, "residue"])
    assert vetoed_res == {20, 21, 22, 23, 24, 25}

    # contact_metrics: top pick at 20 hits contact at 23 (n+3)
    scores = np.array([0.1, 0.1, 0.9, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1])
    cm = CIP.contact_metrics(ev, scores, contact)
    assert cm["contact_top1_within2"] == 1.0

    # Tie-safe percentile ranking: all-equal scores give identical percentile (55.0 for N=10)
    tied = np.full(len(ev), 0.5)
    cm_tied = CIP.contact_metrics(ev, tied, contact)
    assert cm_tied["contact_median_percentile"] == 55.0


def test_contiguous_smoothing_and_reference_zscore():
    import ranker_hard_negative_layer as L

    # Two loops separated by a TM helix gap (residues 1..4 and 50..53)
    df = pd.DataFrame(
        {
            "uniprot_acc": ["P1"] * 8,
            "residue": [1, 2, 3, 4, 50, 51, 52, 53],
            "plddt": [10.0, 20.0, 30.0, 40.0, 90.0, 90.0, 90.0, 90.0],
        }
    )
    scores = np.array([0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0])
    sm = L.smooth_contiguous(df, scores, half=2)
    assert np.allclose(sm[:4], 0.0)
    assert np.allclose(sm[4:], 1.0)

    # Reference z-score: subsampled train rows get the SAME z-scores as full ref rows
    sub = df.iloc[[0, 4]].copy()
    z_sub = L.zscore_from_reference(sub, df, ["plddt"])
    z_ref = L.zscore_from_reference(df, df, ["plddt"])
    assert np.allclose(z_sub["plddt"].to_numpy(), z_ref["plddt"].iloc[[0, 4]].to_numpy())


def test_topology_vs_experiment_deduplicates_pdb_case_across_sources():
    import topology_vs_experiment as TVE

    topo = "OOOOOOOOOO"
    e_uni = TVE.Entry(
        "uniprot",
        "P1",
        "P1_HUMAN",
        "G1",
        "A" * 10,
        [TVE.Region(1, 10, "out", "structure", "PDB:7EPT")],
    )
    e_tdb = TVE.Entry(
        "topdb",
        "P1_HUMAN",
        "P1_HUMAN",
        "G1",
        "A" * 10,
        [TVE.Region(1, 10, "out", "topdb:Xray/", "PDB:7ept")],
    )
    mapped = [
        TVE.Mapped(e_uni, "exact", 1.0, e_uni.regions, [], [], []),
        TVE.Mapped(e_tdb, "exact", 1.0, e_tdb.regions, [], [], []),
    ]
    calls, skipped = TVE._vote(topo, 5, (1, 10), mapped)
    assert not skipped
    assert len(calls) == 1
    c = TVE.call_site(topo, 5, mapped)
    assert c.side == "out"
    assert c.cites == ["PDB:7EPT"]


def test_nterm_decide_cleavage_site_and_flag_provenance():
    import nterm_decide as ND

    # igk_swap must never inherit native cleavage_site (e.g. 28)
    assert ND._opt_cleavage_site({"arm": "igk_swap", "cleavage_site": 20}, 28) == 20
    assert ND._opt_cleavage_site({"arm": "igk_swap"}, 28) == len(ND.LEADER_AA)
    assert ND._opt_cleavage_site({"arm": "native_site"}, 28) == 28

    # finalize must not overwrite Becca's 2026-11-22 flag_not_hold provenance with multipass rule
    row = {
        "candidate_id": "c1",
        "design_id": "l1-p9999",
        "arm": "native_site",
        "insert_after_residue": "25",
        "refold_verdict": "fail: native-region lDDT 0.80 < 0.85",
        "decided_by": ND.FLAG_RULE_BY,
        "evidence": "SignalP",
        "refold_disposition": "refold_flag=fail: native-region lDDT 0.80 < 0.85",
        "igk_replaces_residues": "",
    }
    res = ND.finalize([row], lambda *a, **k: {})[0]
    assert res["refold_rule"] == "built_and_flagged_refold_flag_never_hold"
    assert res["refold_rule_decided_by"] == ND.FLAG_RULE_BY
    assert res["refold_rule_date"] == "2026-11-22"


def test_membrane_side_gate_caches_calibration_and_imports_safely():
    import membrane_side_gate as MSG

    MSG._calibration.cache_clear()
    _ = MSG.thresholds()
    _ = MSG.thresholds()
    assert MSG._calibration.cache_info().hits >= 1

    # Calling _member_state with fresh empty globals still populates from members DataFrame
    members = pd.DataFrame(
        {
            "candidate_id": ["c1"],
            "design_id": ["l1-p0001"],
            "v2_termini_agree": [True],
        }
    )
    cid_map, v2_map = MSG._member_state(members, cid_of_design=None, v2_agrees=None)
    assert cid_map["l1-p0001"] == "c1"
    assert v2_map["c1"] is True


def test_l1_internal_arm_decisions_hygiene_and_audit_flags():
    path = S5E13 / "l1_internal_arm_decisions.tsv"
    with path.open() as fh:
        rows = {r["arm_id"]: r for r in csv.DictReader(fh, delimiter="\t")}

    # PLLP refold_disposition normalized
    assert rows["l1-p1652-a1"]["refold_disposition"] == "fail_flagged"
    assert "native-region lDDT" in rows["l1-p1652-a1"]["flag"]

    # 6 formerly-NaN v5 loop arms have veto_class == "none"
    for aid in (
        "l1-p0113-a5",
        "l1-p0398-a3",
        "l1-p1652-a1",
        "l1-p2081-a5",
        "l1-p3247-a2",
        "l1-p3285-a1",
    ):
        assert rows[aid]["veto_class"] == "none", aid

    # Defective/redundant pre-order audit arms were cleanly dropped or swapped out (2026-11-24)
    for dropped_or_replaced in (
        "l1-p0113-a6",
        "l1-p1961-a3",
        "l1-p1961-a4",
        "l1-p3281-a4",
        "l1-p3281-a5",
        "l1-p0163-a2",
        "l1-p0347-a5",
        "l1-p0349-a2",
        "l1-p1204-a3",
        "l1-p1675-a2",
        "l1-p1181-a1",
        "l1-p0188-a2",
        "l1-p2278-a3",
        "l1-p0846-a4",
        "l1-p2354-a5",
        "l1-p0131-a4",
        "l1-p2015-a3",
    ):
        assert dropped_or_replaced not in rows, dropped_or_replaced

    # Replacement arms are built with passing refolds
    for replacement in (
        "l1-p0163-a1",
        "l1-p0347-a4",
        "l1-p0349-a1",
        "l1-p1204-a2",
        "l1-p0188-a3",
        "l1-p1675-a3",
        "l1-p2278-a2",
        "l1-p2354-a3",
    ):
        assert rows[replacement]["decision"] == "build", replacement
        assert rows[replacement]["refold_verdict"] == "pass", replacement

    # Remaining usable=0 arms carry transparent audit_flag notes
    for aid in ("l1-p1961-a2", "l1-p0163-a1"):
        assert rows[aid]["decision"] == "build", aid
        assert "audit_flag:" in rows[aid]["flag"], aid


def test_cross_folder_paths_resolve_to_repo_files():
    import l1_internal_arms as LIA
    import v8_sources_verdicts as V8SV

    assert V8SV.ARMS.is_file()
    assert LIA.TAG_REFOLD_VERDICTS.is_file()
    assert (LIA.REPO / "src" / "deliverome_analysis").is_dir()
