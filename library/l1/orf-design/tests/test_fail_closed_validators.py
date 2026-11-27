"""Adversarial unit tests for fail-closed validators and builders in s4e9.

Covers silent join drops, duplicate arm/design IDs, NaN/missing values, hard vetoes,
malformed refold dispositions, GPI omega-proximal sites, fail-open fallbacks, and
cross-repo clean-path export compatibility across:
- l1_consistency.py
- l1_completeness.py
- l1_pool_qc.py
- build_l1_internal_arms.py
- build_l1_cterm_arms.py
- build_l1_nterm_arms.py
"""

from __future__ import annotations

import ast
import csv
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parents[1]
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
if str(HERE / "reviewer-package") not in sys.path:
    sys.path.insert(0, str(HERE / "reviewer-package"))

import build_l1_cterm_arms as b_cterm  # noqa: E402
import build_l1_internal_arms as b_int  # noqa: E402
import build_l1_nterm_arms as b_nterm  # noqa: E402
import l1_completeness as comp  # noqa: E402
import l1_consistency as cons  # noqa: E402
import l1_pool_qc as pool  # noqa: E402
import order_plan  # noqa: E402

from deliverome_analysis.cloning import seqtools  # noqa: E402

# --------------------------------------------------------------------------- helpers


def _write_tsv(path: Path, fieldnames: list[str], rows: list[dict]) -> Path:
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fieldnames})
    return path


def _valid_internal_pair() -> tuple[dict, dict]:
    built = {
        "arm_id": "l1-p0001-a1",
        "construct_id": "l1-p0001-a1",
        "design_id": "l1-p0001",
        "gene_symbol": "GENE1",
        "insert_nt": "1200",
        "native_aa": "300",
        "cassette_after_residue": "100",
    }
    decision = {
        "arm_id": "l1-p0001-a1",
        "design_id": "l1-p0001",
        "arm_kind": "loop",
        "decision": "build",
        "veto_class": "none",
        "refold_disposition": "pass",
        "flag": "",
    }
    return built, decision


# --------------------------------------------------------------------------- l1_consistency


def test_consistency_internal_decisions_accepts_valid_rows():
    built1, dec1 = _valid_internal_pair()
    built2 = {
        **built1,
        "arm_id": "l1-p0002-a1",
        "construct_id": "l1-p0002-a1",
        "design_id": "l1-p0002",
    }
    dec2 = {
        **dec1,
        "arm_id": "l1-p0002-a1",
        "design_id": "l1-p0002",
        "veto_class": "soft",
        "refold_disposition": "fail_flagged",
    }
    built3 = {
        **built1,
        "arm_id": "l1-p0003-a1",
        "construct_id": "l1-p0003-a1",
        "design_id": "l1-p0003",
    }
    dec3 = {
        **dec1,
        "arm_id": "l1-p0003-a1",
        "design_id": "l1-p0003",
        "veto_class": "none",
        "refold_disposition": "pending",
    }
    assert cons.check_internal_arm_decisions([built1, built2, built3], [dec1, dec2, dec3]) == {}


def test_consistency_internal_decisions_catches_silent_join_drop():
    built, _ = _valid_internal_pair()
    errs = cons.check_internal_arm_decisions([built], [])
    assert "l1-p0001-a1" in errs
    assert any("missing from l1_internal_arm_decisions.tsv" in e for e in errs["l1-p0001-a1"])


def test_consistency_internal_decisions_rejects_duplicate_arm_id():
    built, dec = _valid_internal_pair()
    # Duplicate in decision table raises loudly
    with pytest.raises(ValueError, match="duplicate arm_id"):
        cons.check_internal_arm_decisions([built], [dec, dict(dec)])
    # Duplicate in built table records mismatch
    errs = cons.check_internal_arm_decisions([built, dict(built)], [dec])
    assert any("duplicate arm_id" in e for e in errs["l1-p0001-a1"])


@pytest.mark.parametrize("bad_veto", ["hard", "HARD", "", "NaN", "nan", "n/a", "unknown"])
def test_consistency_internal_decisions_rejects_hard_or_nan_veto_class(bad_veto: str):
    built, dec = _valid_internal_pair()
    dec["veto_class"] = bad_veto
    errs = cons.check_internal_arm_decisions([built], [dec])
    assert "l1-p0001-a1" in errs
    assert any("veto_class" in e for e in errs["l1-p0001-a1"])


@pytest.mark.parametrize("bad_disp", ["", "NaN", "held", "fail", "fail: native-region lDDT 0.80"])
def test_consistency_internal_decisions_rejects_malformed_refold_disposition(bad_disp: str):
    built, dec = _valid_internal_pair()
    dec["refold_disposition"] = bad_disp
    errs = cons.check_internal_arm_decisions([built], [dec])
    assert "l1-p0001-a1" in errs
    assert any("refold_disposition" in e for e in errs["l1-p0001-a1"])


def test_consistency_internal_decisions_rejects_non_build_or_non_loop():
    built, dec = _valid_internal_pair()
    errs_hold = cons.check_internal_arm_decisions([built], [{**dec, "decision": "hold"}])
    assert any("decision=" in e for e in errs_hold["l1-p0001-a1"])
    errs_snork = cons.check_internal_arm_decisions([built], [{**dec, "arm_kind": "snorkel"}])
    assert any("arm_kind=" in e for e in errs_snork["l1-p0001-a1"])


def test_consistency_biology_audit_flags_gpi_omega_proximal_and_past_omega():
    gpi = [
        {
            "gene_symbol": "GPI1",
            "gpi_anchored": "1",
            "omega_site_reference": "150",
            "propeptide_aa": "20",
        },
        {
            "gene_symbol": "GPI_ALT",
            "gpi_anchored": "1",
            "omega_site_reference": "",
            "propeptide_aa": "25",
        },
        {
            "gene_symbol": "GPI_NO_OMEGA",
            "gpi_anchored": "1",
            "omega_site_reference": "",
            "propeptide_aa": "",
        },
    ]
    rows = [
        # Exactly 12 aa upstream -> allowed (150 - 138 = 12)
        {
            "arm_id": "l1-p0001-a1",
            "design_id": "l1-p0001",
            "gene_symbol": "GPI1",
            "insert_nt": "900",
            "native_aa": "170",
            "cassette_after_residue": "138",
        },
        # 11 aa upstream -> inside omega-11..omega-1 -> flagged
        {
            "arm_id": "l1-p0001-a2",
            "design_id": "l1-p0001",
            "gene_symbol": "GPI1",
            "insert_nt": "900",
            "native_aa": "170",
            "cassette_after_residue": "139",
        },
        # Past omega -> flagged
        {
            "arm_id": "l1-p0001-a3",
            "design_id": "l1-p0001",
            "gene_symbol": "GPI1",
            "insert_nt": "900",
            "native_aa": "170",
            "cassette_after_residue": "152",
        },
        # Fallback via native_aa - propeptide_aa = 200 - 25 = 175; site 170 -> clearance 5 -> flagged
        {
            "arm_id": "l1-p0002-a1",
            "design_id": "l1-p0002",
            "gene_symbol": "GPI_ALT",
            "insert_nt": "900",
            "native_aa": "200",
            "cassette_after_residue": "170",
        },
        # Missing omega on anchored GPI gene -> fail closed
        {
            "arm_id": "l1-p0003-a1",
            "design_id": "l1-p0003",
            "gene_symbol": "GPI_NO_OMEGA",
            "insert_nt": "900",
            "native_aa": "200",
            "cassette_after_residue": "50",
        },
    ]
    audit = cons.check_internal_biology(rows, gpi)
    assert audit["gpi_arms_checked"] == 5
    assert audit["gpi_omega_proximal_count"] == 4
    flagged = {h["arm_id"] for h in audit["gpi_omega_proximal_arms"]}
    assert flagged == {"l1-p0001-a2", "l1-p0001-a3", "l1-p0002-a1", "l1-p0003-a1"}


def test_consistency_biology_audit_rejects_nan_site_or_omega():
    gpi = [
        {
            "gene_symbol": "GPI1",
            "gpi_anchored": "1",
            "omega_site_reference": "NaN",
            "propeptide_aa": "",
        }
    ]
    row = {
        "arm_id": "l1-p0001-a1",
        "design_id": "l1-p0001",
        "gene_symbol": "GPI1",
        "insert_nt": "900",
        "native_aa": "170",
        "cassette_after_residue": "100",
    }
    with pytest.raises(ValueError, match="omega_site_reference"):
        cons.check_internal_biology([row], gpi)
    gpi_ok = [
        {
            "gene_symbol": "GPI1",
            "gpi_anchored": "1",
            "omega_site_reference": "150",
            "propeptide_aa": "20",
        }
    ]
    with pytest.raises(ValueError, match="cassette_after_residue"):
        cons.check_internal_biology([{**row, "cassette_after_residue": "NaN"}], gpi_ok)


def test_consistency_disposition_and_layout_and_fasta_fail_closed(tmp_path: Path):
    # Per-arm mode: missing arm_id never falls back to member-level disposition
    dispo = {
        "internal_arms": {"l1-p0001-a1": {"decision": "build", "flag": "ok"}},
        "internal": {"l1-p0001": {"refold_verdict": "fail"}},
    }
    assert (
        cons.disposition_of(
            "constructs_l1_internal", {"design_id": "l1-p0001", "arm_id": "l1-p0001-a99"}, dispo
        )
        == ""
    )

    # NaN and out-of-range cassette_after_residue in layout_check fail closed
    from build_constructs import CASSETTE_AA

    native = "MKTAYIAKQR"
    tagged = native[:5] + CASSETTE_AA + native[5:]
    assert cons.layout_check({"cassette_after_residue": "NaN"}, tagged, native).startswith(
        "mismatch"
    )
    assert cons.layout_check({"cassette_after_residue": "-1"}, tagged, native).startswith(
        "mismatch"
    )
    assert cons.layout_check({"cassette_after_residue": "999"}, tagged, native).startswith(
        "mismatch"
    )

    # Duplicate FASTA IDs raise ValueError
    fa = tmp_path / "dup.fasta"
    fa.write_text(">l1-p0001-a1|G1\nATGAAA\n>l1-p0001-a1|G1\nATGTTT\n")
    with pytest.raises(ValueError, match="duplicate FASTA"):
        cons.fasta(fa)


# --------------------------------------------------------------------------- l1_completeness


def test_completeness_fail_closed_checks(tmp_path: Path):
    # Strict by_design rejects missing design_id
    with pytest.raises(ValueError, match="missing design_id"):
        comp.by_design([{"_table": "constructs_l1_internal.tsv", "design_id": ""}], strict=True)

    # Duplicate arm_id in arm table raises
    with pytest.raises(ValueError, match="duplicate arm_id"):
        comp._check_arm_table_ids(
            "constructs_l1_internal.tsv",
            [
                {"arm_id": "l1-p0001-a1", "insert_nt": "900"},
                {"arm_id": "l1-p0001-a1", "insert_nt": "900"},
            ],
        )

    # Malformed decision or duplicate arm_id in arm_decisions raises
    bad_dec = _write_tsv(
        tmp_path / "bad_dec.tsv",
        ["arm_id", "design_id", "decision"],
        [{"arm_id": "l1-p0001-a1", "design_id": "l1-p0001", "decision": "NaN"}],
    )
    with pytest.raises(ValueError, match="unexpected decision"):
        comp._load_arm_built(bad_dec)

    dup_dec = _write_tsv(
        tmp_path / "dup_dec.tsv",
        ["arm_id", "design_id", "decision"],
        [
            {"arm_id": "l1-p0001-a1", "design_id": "l1-p0001", "decision": "build"},
            {"arm_id": "l1-p0001-a1", "design_id": "l1-p0001", "decision": "build"},
        ],
    )
    with pytest.raises(ValueError, match="duplicate arm_id"):
        comp._load_arm_built(dup_dec)

    # NaN / out-of-range cassette_after_residue returns False instead of clamping
    assert (
        comp.arm_matches_sequence(
            {"cassette_after_residue": "NaN", "protein_sha256": "abc"}, "MKTAY", "sha"
        )
        is False
    )
    assert (
        comp.arm_matches_sequence(
            {"cassette_after_residue": "999", "protein_sha256": "abc"}, "MKTAY", "sha"
        )
        is False
    )

    # Multipass arm missing from multipass_refold_dispositions does NOT silently default to pass
    # when tag_refold_verdicts shows a fail or is absent
    arm = {"construct_id": "l1-p0001:tag", "protein_sha256": "sha1"}
    assert "fail" in comp._multipass_verdict(arm, {}, {"sha1": {"verdict": "fail"}})
    assert comp._multipass_verdict(arm, {}, {"other_sha": {"verdict": "pass"}}) == "not recorded"
    assert comp._multipass_verdict(arm, {}, {"sha1": {"verdict": "pass"}}) == "pass"


# --------------------------------------------------------------------------- l1_pool_qc


def test_pool_qc_fail_closed_checks(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    # Missing protein hash and missing native are errors (not ok=True warnings)
    ins = {
        "arm_id": "l1-p0001-a1",
        "rung": "internal",
        "insert": pool.OPENING + "GCC" * 5 + "TAA" + pool.READOUT_BLOCK,
        "row": {"protein_sha256": "", "mature_sha256": ""},
    }
    findings = {f.check: f for f in pool.extra_findings(ins, native="")}
    assert not findings["protein_matches_table"].ok
    assert findings["protein_matches_table"].severity == "error"
    assert not findings["layout"].ok
    assert findings["layout"].severity == "error"

    # NaN or negative cassette_after_residue raises
    with pytest.raises(ValueError, match="cassette_after_residue"):
        pool.cassette_site(
            {"arm_id": "a1", "rung": "internal", "row": {"cassette_after_residue": "NaN"}}
        )
    with pytest.raises(ValueError, match="cassette_after_residue"):
        pool.cassette_site(
            {"arm_id": "a1", "rung": "internal", "row": {"cassette_after_residue": "-2"}}
        )

    # check_pool_completeness fails closed on NaN arm count and duplicate design_id
    inserts = [{"design_id": "l1-p0001", "arm_id": "l1-p0001-a1"}]
    res_nan = pool.check_pool_completeness(
        inserts,
        [],
        completeness_rows=[{"design_id": "l1-p0001", "status": "built", "arms": "NaN"}],
    )
    assert not res_nan.ok
    assert "l1-p0001" in res_nan.detail

    with pytest.raises(ValueError, match="duplicate design_id"):
        pool.check_pool_completeness(
            inserts,
            [],
            completeness_rows=[
                {"design_id": "l1-p0001", "status": "built", "arms": "1"},
                {"design_id": "l1-p0001", "status": "built", "arms": "1"},
            ],
        )

    # load_inserts fails loudly on silent TSV->FASTA row drop
    _write_tsv(
        tmp_path / "constructs_l1_type3.tsv",
        ["candidate_id", "design_id", "gene_symbol", "insert_nt"],
        [
            {
                "candidate_id": "l1-p0001",
                "design_id": "l1-p0001",
                "gene_symbol": "G1",
                "insert_nt": "100",
            },
            {
                "candidate_id": "l1-p0002",
                "design_id": "l1-p0002",
                "gene_symbol": "G2",
                "insert_nt": "100",
            },
        ],
    )
    (tmp_path / "constructs_l1_type3.fasta").write_text(">l1-p0001|G1\nATGAAATAA\n")
    monkeypatch.setattr(pool, "HERE", tmp_path)
    monkeypatch.setattr(pool, "RUNGS", {"type3": "constructs_l1_type3"})
    monkeypatch.setattr(pool.arm_ids, "load_decisions", lambda: {})
    with pytest.raises(ValueError, match="built row\\(s\\) missing from"):
        pool.load_inserts()


# --------------------------------------------------------------------------- builders


def test_build_internal_arms_fail_closed(tmp_path: Path):
    cols = [
        "arm_id",
        "design_id",
        "sequence_sha256",
        "arm_role",
        "arm_kind",
        "insert_after_residue",
        "decision",
        "veto_class",
        "refold_disposition",
        "flag",
    ]
    base = {
        "arm_id": "l1-p0001-a1",
        "design_id": "l1-p0001",
        "sequence_sha256": "abc",
        "arm_role": "loop_best",
        "arm_kind": "loop",
        "insert_after_residue": "50",
        "decision": "build",
        "veto_class": "none",
        "refold_disposition": "pass",
        "flag": "",
    }
    # Unknown decision (typo/NaN) must raise, never silently skip
    p_dec = _write_tsv(tmp_path / "bad_dec.tsv", cols, [{**base, "decision": "Build"}])
    with pytest.raises(AssertionError, match="unknown decision"):
        b_int.arm_calls(p_dec)

    # Hard veto_class must raise
    p_hard = _write_tsv(tmp_path / "hard.tsv", cols, [{**base, "veto_class": "hard"}])
    with pytest.raises(AssertionError, match="veto_class"):
        b_int.arm_calls(p_hard)

    # NaN veto_class must raise
    p_nan = _write_tsv(tmp_path / "nan_veto.tsv", cols, [{**base, "veto_class": "NaN"}])
    with pytest.raises(AssertionError, match="veto_class"):
        b_int.arm_calls(p_nan)

    # Malformed refold_disposition must raise
    p_rd = _write_tsv(tmp_path / "bad_rd.tsv", cols, [{**base, "refold_disposition": "fail: lddt"}])
    with pytest.raises(AssertionError, match="refold_disposition"):
        b_int.arm_calls(p_rd)

    # Fractional or NaN insert_after_residue in build_one raises ValueError
    with pytest.raises(ValueError, match="insert_after_residue"):
        b_int.build_one(
            {"design_id": "l1-p0001"}, {"insert_after_residue": "NaN"}, "ATGAAATAA", "MK", {}
        )
    with pytest.raises(ValueError, match="insert_after_residue"):
        b_int.build_one(
            {"design_id": "l1-p0001"}, {"insert_after_residue": "10.5"}, "ATGAAATAA", "MK", {}
        )


def test_build_cterm_arms_fail_closed():
    row = {
        "design_id": "l1-p9999",
        "gene_symbol": "G1",
        "species": "human",
        "uniprot_acc": "P00001",
        "cterm_tag": "kept",
        "flags": "",
        "uniprot_cterm_topology": "",
        "v2_cterm_side": "non-cytoplasmic",
        "uniprot_signal_peptide": "",
        "uniprot_tm_helices": "1",
        "stalk_site_after": "",
    }
    # Duplicate design_id in review raises
    with pytest.raises(ValueError, match="duplicate design_id"):
        b_cterm.dispositions([row, dict(row)], verdicts={}, stalk_verdicts={})

    # Unknown cterm_tag raises instead of silently treating as kept
    with pytest.raises(AssertionError, match="unexpected cterm_tag"):
        b_cterm.dispositions([{**row, "cterm_tag": "NaN"}], verdicts={}, stalk_verdicts={})


def test_build_nterm_arms_fail_closed(monkeypatch: pytest.MonkeyPatch):
    codon = {aa: c for c, aa in seqtools.GENETIC_CODE.items() if aa != "*"}
    protein = "MKTAYIAKQRQISFVK"
    orf = "".join(codon[a] for a in protein) + "TAA"

    # Unknown arm value raises
    with pytest.raises(AssertionError, match="unexpected arm"):
        b_nterm.build_one(
            {"candidate_id": "l1-p0001", "arm": "native", "native_after_residue": "5"},
            orf,
            protein,
            {},
        )
    # NaN or out-of-range native_after_residue raises
    with pytest.raises(ValueError, match="native_after_residue"):
        b_nterm.build_one(
            {"candidate_id": "l1-p0001", "arm": "native_site", "native_after_residue": "NaN"},
            orf,
            protein,
            {},
        )
    with pytest.raises(AssertionError, match="out of range"):
        b_nterm.build_one(
            {"candidate_id": "l1-p0001", "arm": "native_site", "native_after_residue": "999"},
            orf,
            protein,
            {},
        )

    # Unknown arm value in nterm_decisions.tsv fails build() loudly (no silent skip)
    monkeypatch.setattr(
        b_nterm,
        "read_table",
        lambda path: [{"candidate_id": "l1-p0001", "design_id": "l1-p0001", "arm": "typo_arm"}],
    )
    with pytest.raises(AssertionError, match="unknown arm values"):
        b_nterm.build()


# --------------------------------------------------------------------------- order_plan


def test_order_plan_fails_closed_on_stale_or_mismatched_numbering(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """order_plan must reject missing, extra, duplicate, or length-mismatched arms."""
    num_rows = [
        {
            "insert_id": "l1-i0001",
            "plate": "1",
            "well": "A01",
            "arm_id": "l1-p0001-a1",
            "design_id": "l1-p0001",
            "rung": "cterm",
            "ordered_insert_nt": "1200",
            "key": "lib:GENE1",
        },
        {
            "insert_id": "l1-i0002",
            "plate": "1",
            "well": "B01",
            "arm_id": "l1-p0002-a1",
            "design_id": "l1-p0002",
            "rung": "nterm",
            "ordered_insert_nt": "1500",
            "key": "lib:GENE2",
        },
    ]
    comp_rows = [
        {"arm_id": "l1-p0001-a1", "insert_id": "l1-i0001", "insert_nt": "1200"},
        {"arm_id": "l1-p0002-a1", "insert_id": "l1-i0002", "insert_nt": "1500"},
    ]

    # Exact match passes
    order_plan.verify_numbering_matches_compliance(num_rows, comp_rows)

    # Missing arm in l1_insert_numbering.tsv raises ValueError
    with pytest.raises(ValueError, match=r"missing from l1_insert_numbering\.tsv"):
        order_plan.verify_numbering_matches_compliance(num_rows[:1], comp_rows)

    # Extra arm in l1_insert_numbering.tsv raises ValueError
    with pytest.raises(ValueError, match=r"extra in l1_insert_numbering\.tsv"):
        order_plan.verify_numbering_matches_compliance(num_rows, comp_rows[:1])

    # Length mismatch between ordered_insert_nt and insert_nt raises ValueError
    stale_len = [{**num_rows[0], "ordered_insert_nt": "1203"}, num_rows[1]]
    with pytest.raises(ValueError, match="insert_nt mismatch"):
        order_plan.verify_numbering_matches_compliance(stale_len, comp_rows)

    # Duplicate arm_id raises ValueError
    with pytest.raises(ValueError, match="duplicate arm_id"):
        order_plan.verify_numbering_matches_compliance([num_rows[0], dict(num_rows[0])], comp_rows)

    # End-to-end inserts() fails closed on a stale l1_insert_numbering.tsv on disk
    num_path = _write_tsv(
        tmp_path / "l1_insert_numbering.tsv",
        ["insert_id", "plate", "well", "arm_id", "design_id", "rung", "ordered_insert_nt", "key"],
        stale_len,
    )
    comp_path = _write_tsv(
        tmp_path / "synthesis_compliance.tsv",
        ["arm_id", "insert_id", "insert_nt"],
        comp_rows,
    )
    manifest_path = _write_tsv(
        tmp_path / "candidate_manifest_tranches.tsv",
        ["design_id", "library", "retired_reason", "set", "insert_nt", "vector_bp"],
        [
            {
                "design_id": "l1-p0001",
                "library": "l1",
                "retired_reason": "",
                "set": "library_reference",
                "insert_nt": "1047",
                "vector_bp": "6000",
            },
            {
                "design_id": "l1-p0002",
                "library": "l1",
                "retired_reason": "",
                "set": "library_reference",
                "insert_nt": "1347",
                "vector_bp": "6300",
            },
        ],
    )
    monkeypatch.setattr(order_plan, "NUMBERING", num_path)
    monkeypatch.setattr(order_plan, "COMPLIANCE", comp_path)
    monkeypatch.setattr(order_plan, "MANIFEST", manifest_path)
    with pytest.raises(ValueError, match="insert_nt mismatch"):
        order_plan.inserts()


# --------------------------------------------------------------------------- cross-repo path portability


def test_hardened_scripts_export_cleanly_for_review_repo():
    import export_review_repo as err

    owned = [
        "library/l1/orf-design/l1_consistency.py",
        "library/l1/orf-design/l1_completeness.py",
        "library/l1/orf-design/l1_pool_qc.py",
        "library/l1/orf-design/build_l1_internal_arms.py",
        "library/l1/orf-design/build_l1_cterm_arms.py",
        "library/l1/orf-design/build_l1_nterm_arms.py",
    ]
    repo = HERE.parents[2]
    raw = {rel: (repo / rel).read_bytes() for rel in owned}
    remapped = err.apply_clean_paths(raw)
    assert len(remapped) == len(owned)
    for dest, data in remapped.items():
        assert dest.startswith("library/l1/orf-design/")
        text = data.decode("utf-8")
        ast.parse(text, filename=dest)
        assert "s5e13-tag-site-benchmarking" not in text
        assert "scratch/becca/" not in text
