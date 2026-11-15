"""Every Igκ swap on the N-after-SP rung names the evidence that residues 1..c are a cleaved
signal peptide (coordinator, 2026-11-14, after the fresh-eyes audit found 66 deletions resting
on DeepTMHMM v1 alone)."""

from __future__ import annotations

import csv
import sys
from pathlib import Path

NTERM = Path(__file__).resolve().parents[1] / "reviewer-package" / "nterm"
sys.path.insert(0, str(NTERM))


def _rows() -> list[dict]:
    with (NTERM / "nterm_decisions.tsv").open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _census() -> dict[str, dict]:
    # Keyed on design id: nterm_members.tsv is a dated census whose member ids are from the
    # numbering of its day, so a join on candidate_id names another protein after a renumber.
    with (NTERM / "nterm_members.tsv").open() as fh:
        return {r["design_id"]: r for r in csv.DictReader(fh, delimiter="\t")}


def test_every_igk_deletion_names_its_evidence():
    import nterm_decide as nd

    for r in _rows():
        if r["arm"] != "igk_swap":
            continue
        ev = r["igk_deletion_evidence"]
        assert ev, f"{r['candidate_id']} {r['gene_symbol']}: Igκ swap with no deletion evidence"
        c = r["igk_replaces_residues"].split("-")[1]
        assert f"deletes 1-{c}" in ev or ev.startswith("none:"), r["candidate_id"]
        if ev.startswith("none:"):
            assert r["design_id"] in nd.DEFERRED_ROUTES, (
                f"{r['candidate_id']}: Igκ deletion without evidence outside DEFERRED_ROUTES"
            )


def test_signalp_negative_members_need_a_second_source():
    """No Igκ swap rests on DeepTMHMM v1 where SignalP 6.0 calls no native signal peptide."""
    census = _census()
    for r in _rows():
        if r["arm"] != "igk_swap" or r["design_id"] in __import__("nterm_decide").DEFERRED_ROUTES:
            continue
        cls = census.get(r["design_id"], {}).get("class", "")
        ev = r["igk_deletion_evidence"]
        if cls == "refused_sp_disagree_deeptmhmm_v1_only" or not ev.startswith("SignalP 6.0"):
            assert ev.startswith(("second source:", "re-anchored:")), (
                f"{r['candidate_id']} {r['gene_symbol']}: SignalP-negative Igκ swap without a "
                f"second source ({ev!r})"
            )
            assert r["sp_review_verdict"] in ("cleaved_sp_confirmed", "unconfirmed")


def test_reanchors_rest_on_curated_or_published_evidence():
    for r in _rows():
        if r["arm"] == "igk_swap" and r["igk_deletion_evidence"].startswith("re-anchored:"):
            ev = r["igk_deletion_evidence"]
            assert "UniProt" in ev or "PMID" in ev, f"{r['candidate_id']}: re-anchor on {ev!r}"


def test_deferred_routes_are_few_and_named():
    import nterm_decide as nd

    assert len(nd.DEFERRED_ROUTES) <= 1
    for did, why in nd.DEFERRED_ROUTES.items():
        assert "refold" in why and "PR #" in why, did


def test_every_built_arm_that_fails_refold_carries_a_disposition():
    """A built arm with a failing refold verdict names the rule it is built under (2026-11-14)."""
    for r in _rows():
        if r["arm"] not in ("native_site", "igk_swap") or not r["refold_verdict"].startswith(
            "fail"
        ):
            continue
        assert r["refold_disposition"].startswith("built and flagged"), (
            f"{r['candidate_id']} {r['gene_symbol']}: built with a failing refold and no disposition"
        )
        assert "multipass rung" in r["refold_disposition"], r["candidate_id"]
        assert r["refold_rule"] == "built_and_flagged_multipass_disposition_rule", r["candidate_id"]


def test_no_ortholog_is_built_without_a_built_human_member():
    """Becca, 2026-11-14: never an ortholog in place of a human member."""
    import nterm_decide as nd

    rows = _rows()
    here = {
        r["gene_symbol"]
        for r in rows
        if r["arm"] in ("native_site", "igk_swap") and r["candidate_set"] != "ortholog"
    }
    built = here | nd.human_built_elsewhere()
    for r in rows:
        if r["candidate_set"] == "ortholog" and r["arm"] in ("native_site", "igk_swap"):
            assert r["gene_symbol"] in built, f"{r['design_id']} {r['gene_symbol']}: ortholog alone"
