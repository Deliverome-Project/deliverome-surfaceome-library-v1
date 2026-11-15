"""The whole-library L1 checks: per-arm consistency and per-member completeness.

Unit tests pin the pieces each check rests on (the insert translation, the layout rule); the
integration tests run both checks over the committed tables, so a rung table that drifts from
the manifest, the folds or the refold verdicts fails CI. Native proteins come from committed
tables; only a member missing from them needs S3, and then the test skips (never passes) when no
credentials are available, like the repo's other live tests.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))

import l1_completeness  # noqa: E402
import l1_consistency  # noqa: E402
from build_constructs import CASSETTE_AA, KOZAK  # noqa: E402

from deliverome_analysis import seqtools  # noqa: E402

CODON = {aa: c for c, aa in seqtools.GENETIC_CODE.items() if aa != "*"}


def _insert(protein: str) -> str:
    """A wrapped insert: Esp3I site, Kozak, the ORF, a stop, then a barcode-like tail."""
    return "CGTCTCA" + KOZAK + "".join(CODON[a] for a in protein) + "TAA" + "NNNNACGT"


def test_translate_insert_reads_kozak_atg_to_first_stop():
    protein = "MKTAYIAKQR"
    assert l1_consistency.translate_insert(_insert(protein)) == protein


def test_translate_insert_refuses_an_insert_without_kozak_atg():
    with pytest.raises(AssertionError):
        l1_consistency.translate_insert("CGTCTCA" + "ATGAAATAA")
    with pytest.raises(AssertionError):
        l1_consistency.translate_insert("CGTCTCA" + KOZAK + "GTGAAATAA")


def test_layout_check_accepts_the_cassette_at_the_recorded_site_only():
    native = "MKTAYIAKQRQISFVKSHFSRQ"
    row = {"cassette_after_residue": "5"}
    tagged = native[:5] + CASSETTE_AA + native[5:]
    assert l1_consistency.layout_check(row, tagged, native) == "ok"
    shifted = native[:6] + CASSETTE_AA + native[6:]
    assert l1_consistency.layout_check(row, shifted, native).startswith("mismatch")
    assert l1_consistency.layout_check(row, tagged + CASSETTE_AA, native).startswith("mismatch")


def test_layout_check_is_structural_for_snorkel_and_igk():
    native = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQ"
    snorkel = native + "LLVAAAIIGGS" + CASSETTE_AA
    assert l1_consistency.layout_check({"placement": "snorkel"}, snorkel, native).startswith("ok")
    igk = "METDTLLLWVLLLWVPGSTG" + CASSETTE_AA + native[25:]
    assert l1_consistency.layout_check({"design": "igk_tag"}, igk, native).startswith("ok")


def _offline_or_skip(fn):
    try:
        return fn()
    except l1_consistency.NativesUnavailable as exc:
        pytest.skip(f"needs S3 for a native protein not committed: {exc}")


@pytest.fixture(scope="module")
def consistency():
    return _offline_or_skip(l1_consistency.build)


def test_every_built_arm_is_consistent(consistency):
    rows, _ = consistency
    assert rows, "no built arms found"
    bad = [r for r in rows if r["mismatch"]]
    assert not bad, f"{len(bad)} mismatches, e.g. {bad[:3]}"


def test_every_built_arm_has_a_live_fold_and_a_recorded_verdict(consistency):
    rows, _ = consistency
    no_fold = [r["arm_id"] for r in rows if not r["fold"].startswith("folded")]
    no_verdict = [r["arm_id"] for r in rows if r["verdict"] not in ("pass", "fail")]
    assert not no_fold, f"{len(no_fold)} built arms without a live fold, e.g. {no_fold[:5]}"
    assert not no_verdict, f"{len(no_verdict)} built arms without a verdict: {no_verdict[:5]}"


def test_completeness_places_every_active_member_on_its_current_ids():
    rows, summary = _offline_or_skip(l1_completeness.build)
    assert summary["active_members"] == len(rows)
    assert summary["status"].get("unplaced", 0) == 0
    assert summary["stale_member_id"] == 0
    assert summary["stale_sequence"] == 0


def test_disposition_rules_per_rung():
    dispo = {
        "multipass": {"l1-p0001:tag_alone": "keep (SpyTag003 buried; bench decides)"},
        "nterm": {"l1-p0002": "lddt_borderline_counted_as_pass"},
        "cterm_keep": {"l1-p0003": ["keep the C-terminal tag"], "l1-p0004": ["keep after 235"]},
        "internal": {"l1-p0005": {"decision": "hold_no_passing_arm"}, "l1-p0006": {}},
    }
    d = l1_consistency.disposition_of
    assert d("constructs_l1_multipass", {"construct_id": "l1-p0001:tag_alone"}, dispo)
    assert not d("constructs_l1_multipass", {"construct_id": "l1-p0009:tag_alone"}, dispo)
    assert d("constructs_l1_type3", {"refold_disposition": "kept: helix flag"}, dispo)
    assert not d("constructs_l1_type3", {"refold_disposition": ""}, dispo)
    assert d("constructs_l1_nterm", {"design_id": "l1-p0002"}, dispo)
    c_term = {"design_id": "l1-p0003", "placement": "C-terminal"}
    assert d("constructs_l1_cterm", c_term, dispo)
    # A keep decision covers a re-sited arm only if it names that arm's site.
    stalk = {
        "design_id": "l1-p0003",
        "placement": "internal, stalk",
        "cassette_after_residue": "46",
    }
    assert not d("constructs_l1_cterm", stalk, dispo)
    furin = {"design_id": "l1-p0004", "placement": "internal", "cassette_after_residue": "235"}
    assert d("constructs_l1_cterm", furin, dispo)
    # A failing snorkel is disposed of only by the ranker's hold on its member.
    assert d("constructs_l1_cterm", {"design_id": "l1-p0005", "placement": "snorkel"}, dispo)
    assert not d("constructs_l1_cterm", {"design_id": "l1-p0006", "placement": "snorkel"}, dispo)


def test_every_failing_built_arm_has_a_recorded_disposition(consistency):
    rows, _ = consistency
    bare = [
        f"{r['gene_symbol']} {r['design_id']}"
        for r in rows
        if r["verdict"] == "fail" and r["disposition"].startswith("mismatch")
    ]
    assert not bare, f"{len(bare)} failing built arms without a disposition: {bare[:10]}"


ARM_COLS = ["arm_id", "design_id", "sequence_sha256", "arm_role", "arm_kind"]
ARM_COLS += ["insert_after_residue", "decision", "flag"]


def _arm_table(tmp_path, rows):
    import csv

    path = tmp_path / "arms.tsv"
    with path.open("w") as fh:
        w = csv.DictWriter(fh, fieldnames=ARM_COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, "") for c in ARM_COLS})
    return path


def test_per_arm_decisions_build_loops_and_reject_a_lone_snorkel(tmp_path):
    import build_l1_internal_arms as b

    ok = _arm_table(
        tmp_path,
        [
            {
                "arm_id": "l1-p0001-a1",
                "design_id": "l1-p0001",
                "arm_kind": "loop",
                "decision": "build",
            },
            {
                "arm_id": "l1-p0001-a2",
                "design_id": "l1-p0001",
                "arm_kind": "loop",
                "decision": "build",
            },
            {
                "arm_id": "l1-p0001-a3",
                "design_id": "l1-p0001",
                "arm_kind": "snorkel",
                "decision": "build",
            },
            {
                "arm_id": "l1-p0002-a1",
                "design_id": "l1-p0002",
                "arm_kind": "loop",
                "decision": "hold",
            },
        ],
    )
    calls = b.arm_calls(ok)
    assert [a["arm_id"] for a in calls["l1-p0001"]] == ["l1-p0001-a1", "l1-p0001-a2"]
    assert "l1-p0002" not in calls
    lone = _arm_table(
        tmp_path,
        [
            {
                "arm_id": "l1-p0003-a1",
                "design_id": "l1-p0003",
                "arm_kind": "snorkel",
                "decision": "build",
            }
        ],
    )
    with pytest.raises(AssertionError, match="only arm"):
        b.arm_calls(lone)
    four = _arm_table(
        tmp_path,
        [
            {
                "arm_id": f"l1-p0004-a{k}",
                "design_id": "l1-p0004",
                "arm_kind": "loop",
                "decision": "build",
            }
            for k in range(1, 5)
        ],
    )
    with pytest.raises(AssertionError, match="at most 3"):
        b.arm_calls(four)
    assert b.arm_calls(tmp_path / "absent.tsv") is None


def test_arm_shape_rules():
    shape = l1_completeness.arm_shape_problem
    loop, snorkel = {"placement": "internal"}, {"placement": "snorkel (TM module, then cassette)"}
    assert shape("internal", [loop, loop, snorkel]) == ""
    assert shape("internal", [snorkel]) == "snorkel_only"
    assert shape("N after SP", [loop, loop]) == "multi_arm_on_single_arm_rung"
    assert shape("N after SP", [loop]) == ""


def test_a_failing_internal_arm_is_disposed_only_by_its_own_flag():
    d = l1_consistency.disposition_of
    dispo = {
        "internal_arms": {
            "l1-p0001-a2": {"decision": "build", "flag": "built so the snorkel is not alone"},
            "l1-p0001-a3": {"decision": "build", "flag": ""},
        },
        "internal": {"l1-p0001": {"refold_verdict": "fail"}},
    }
    row = {"design_id": "l1-p0001", "arm_id": "l1-p0001-a2"}
    assert d("constructs_l1_internal", row, dispo).startswith("internal arm built with flag")
    assert d("constructs_l1_internal", {**row, "arm_id": "l1-p0001-a3"}, dispo) == ""


def test_orthologs_need_a_built_human_member_rule():
    rows = [
        {"design_id": "h1", "gene_symbol": "A", "species": "human", "status": "built"},
        {"design_id": "m1", "gene_symbol": "A", "species": "mouse", "status": "built"},
        {"design_id": "h2", "gene_symbol": "B", "species": "human", "status": "held"},
        {"design_id": "c2", "gene_symbol": "B", "species": "cyno", "status": "built"},
        {"design_id": "c3", "gene_symbol": "C", "species": "cyno", "status": "held"},
    ]
    # B's cyno is built while its only human member is held: not allowed.
    assert l1_completeness.orthologs_without_built_human(rows) == {"c2"}


def test_no_ortholog_is_built_without_a_built_human_member():
    # Becca, 2026-11-14: never substitute a mouse or cyno ortholog for a human member.
    _, summary = _offline_or_skip(l1_completeness.build)
    assert summary["orthologs_without_built_human"] == []


def test_every_member_with_a_built_arm_counts_as_built():
    # In per-arm mode the arm rows govern: a stale member-level hold (SLC35G2 and MS4A1 kept
    # hold_no_passing_arm in l1_internal_decisions.tsv after their arms were built under the
    # control-relative refold bar) must not hide built arms, 2026-11-15.
    import csv

    rows, _ = _offline_or_skip(l1_completeness.build)
    built = {r["design_id"] for r in rows if r["status"] == "built"}
    armed = set()
    for stem in l1_consistency.RUNGS:
        path = HERE / f"{stem}.tsv"
        if path.exists():
            with path.open() as fh:
                armed |= {
                    r["design_id"]
                    for r in csv.DictReader(fh, delimiter="\t")
                    if (r.get("insert_nt") or "").strip()
                }
    assert armed == built, f"armed but not built: {sorted(armed - built)[:5]}"
