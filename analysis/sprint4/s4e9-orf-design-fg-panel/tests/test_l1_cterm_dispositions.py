"""The C-terminal rung builds what its two rules say, and records every departure.

``build_l1_cterm_arms.dispositions()`` applies ``STALK_RULE`` and ``COUT_RULE`` to committed
tables. These tests check the rules on small cases, then check the built table against them:
no arm whose C-terminal refold fails is built without a recorded disposition, no arm is built
at a stalk site whose own refold fails, and no weak C-out call is built without a decision
row. Reads committed tables only. Members join on design_id.
"""

from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))

import build_l1_cterm_arms as arms  # noqa: E402


def _tsv(name: str) -> list[dict]:
    with (HERE / name).open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _by_design(name: str) -> dict[str, dict]:
    return {r["design_id"]: r for r in _tsv(name)}


def _decisions() -> dict[str, list[dict]]:
    out: dict[str, list[dict]] = {}
    for r in _tsv("l1_cterm_decisions.tsv"):
        # A scope names the member and its design id; since the per-insert cut they are equal.
        for d in dict.fromkeys(re.findall(r"\bl1-p[0-9]{4}\b", r["scope"])):
            out.setdefault(d, []).append(r)
    return out


#: Every table that holds built arms (one row per arm, keyed by design_id).
ARM_TABLES = (
    "constructs_l1_cterm.tsv",
    "constructs_l1_internal.tsv",
    "constructs_l1_multipass.tsv",
    "constructs_l1_nterm.tsv",
    "constructs_l1_rescue.tsv",
    "constructs_l1_type3.tsv",
)

REVIEW = {"design_id": "l1-p9999", "gene_symbol": "GENE", "stalk_site_after": "40"}
FAIL = {"verdict": "fail", "verdict_reasons": "6 epitope residues confidently packed"}
PASS = {"verdict": "pass", "verdict_reasons": ""}


def test_stalk_rule_clauses():
    rule = arms.refold_disposition
    assert rule(REVIEW, PASS, None, False).action == arms.CTERM
    # A failing C-terminal refold moves only to a stalk site that passes its own refold.
    moved = rule(REVIEW, FAIL, PASS, False)
    assert (moved.action, moved.site, moved.clause) == (arms.STALK, 40, "stalk_refold_pass")
    assert rule(REVIEW, FAIL, FAIL, False).action == arms.HOLD
    assert rule(REVIEW, FAIL, None, False).clause == "stalk_refold_pending"
    # A natively released ectodomain: the cassette stays at the C-terminus.
    assert rule(REVIEW, FAIL, PASS, True).clause == "released_ectodomain"
    # With a stalk site the trigger is the verdict, whatever its reason.
    fold = {"verdict": "fail", "verdict_reasons": "native-region lDDT 0.79 < 0.85"}
    assert rule(REVIEW, fold, PASS, False).action == rule(REVIEW, FAIL, PASS, False).action
    # With no stalk site the reason decides: packing alone keeps the tag, a failed fold holds.
    no_stalk = {**REVIEW, "stalk_site_after": ""}
    kept = rule(no_stalk, FAIL, None, False)
    assert (kept.action, kept.clause) == (arms.CTERM, "no_stalk_site")
    site = {"verdict": "fail", "verdict_reasons": "lDDT within 15 of the site 0.388 < 0.7"}
    for failed_fold in (fold, site):
        held = rule(no_stalk, failed_fold, None, False)
        assert (held.action, held.clause) == (arms.HOLD, "fold_fail_no_stalk")


def test_cout_rule_clauses():
    base = {
        "design_id": "l1-p9999",
        "flags": "",
        "uniprot_cterm_topology": "",
        "v2_cterm_side": "non-cytoplasmic",
        "uniprot_signal_peptide": "",
        "uniprot_tm_helices": "3",
    }
    rule = arms.cout_disposition
    assert rule(base).clause == "not_weak"
    assert rule({**base, "flags": "short_ecto_tail"}).action == arms.CTERM
    v2 = {**base, "flags": "v2_cterm_cytoplasmic", "v2_cterm_side": "cytoplasmic"}
    assert rule(v2).action == arms.HOLD
    assert rule({**v2, "uniprot_cterm_topology": "Extracellular 1-9"}).action == arms.CTERM
    assert rule({**v2, "flags": v2["flags"] + ";lane_override_topology"}).action == arms.CTERM
    # UniProt cytoplasmic: experimental holds; type I holds; weak holds; a bare prediction is
    # kept and flagged.
    experimental = {**base, "uniprot_cterm_topology": "Cytoplasmic 5-9 (ECO:0000269)"}
    assert rule(experimental).clause == "uniprot_experimental_cytoplasmic"
    predicted = {**base, "uniprot_cterm_topology": "Cytoplasmic 5-9 (ECO:0000255)"}
    assert rule(predicted).clause == "predicted_cytoplasmic_kept"
    assert rule(predicted).action == arms.CTERM
    assert rule({**predicted, "flags": "no_v2"}).clause == "uniprot_cytoplasmic"
    type_i = {**predicted, "uniprot_signal_peptide": "1-24", "uniprot_tm_helices": "1"}
    assert rule(type_i).clause == "uniprot_type_i"
    assert rule({**type_i, "uniprot_tm_helices": "7"}).action == arms.CTERM
    # An ortholog with no UniProt entry follows its human member only for those two holds.
    ortholog = {**base, "flags": "uniprot_not_checked"}
    assert rule(ortholog, type_i).clause == "ortholog_of_held_human"
    assert rule(ortholog, experimental).action == arms.HOLD
    assert rule(ortholog, predicted).action == arms.CTERM
    assert rule(ortholog, v2).action == arms.CTERM
    assert {rule(x).clause for x in (experimental, type_i, v2)} <= set(arms.COUT_HOLD_CLAUSES)


def test_built_arms_are_what_the_rules_give():
    review = _by_design("l1_cterm_review.tsv")
    built = {
        d: r
        for d, r in _by_design("constructs_l1_cterm.tsv").items()
        if not r["placement"].startswith("snorkel")
    }
    dispo = arms.dispositions()
    assert set(dispo) == set(review)
    placement = {
        arms.CTERM: "C-terminal",
        arms.STALK: "internal, stalk (refold)",
        arms.RELOCATE: "internal, before furin site",
    }
    for d, x in dispo.items():
        if x.action == arms.HOLD:
            assert d not in built, f"{d}: held by the rule but built"
            continue
        assert d in built, f"{d}: the rule builds it, but it has no arm"
        assert built[d]["placement"] == placement[x.action], d
        if x.action == arms.STALK:
            assert built[d]["cassette_after_residue"] == str(x.site), d
            assert x.clause == "stalk_refold_pass", f"{d}: built at a stalk site not yet refolded"


def test_no_failing_refold_is_built_without_a_recorded_disposition():
    built = _by_design("constructs_l1_cterm.tsv")
    verdicts = _by_design("l1_cterm_refold_verdicts.tsv")
    stalk = _by_design("l1_cterm_stalk_refold_verdicts.tsv")
    decisions = _decisions()
    dispo = arms.dispositions()
    # Members that have since left the C-terminal rung (the 2026-11-15 rescue round moved VTCN1,
    # PLVAP and TRPV1 to the internal rung) keep their old verdict rows but have no C-terminal
    # disposition; only members still on the rung are checked here.
    failing = [d for d, v in verdicts.items() if v["verdict"] == "fail" and d in dispo]
    assert failing, "no failing C-terminal refold: the verdict table was not read"
    for d in failing:
        rows = [r for r in decisions.get(d, []) if r["question"] == "tagged refold failed"]
        on_hold_for_cout = dispo[d].clause in arms.COUT_HOLD_CLAUSES
        assert rows or on_hold_for_cout, f"{d}: failing refold with no decision row"
        for r in rows:
            assert r["decided_by"] and r["reason"] and r["rejected_alternatives"], d
        if d in built and built[d]["placement"] == "C-terminal":
            # Still tagged at the failing C-terminus: only where the rule has no stalk site to
            # move to, or the ectodomain is released.
            assert dispo[d].clause in ("no_stalk_site", "released_ectodomain"), (d, dispo[d])
            assert rows[0]["decision"].startswith("keep"), d
    # A stalk-site arm is built only on a passing verdict for that exact protein.
    for d, r in built.items():
        if r["placement"] == "internal, stalk (refold)":
            assert stalk[d]["verdict"] == "pass", f"{d}: stalk-site refold {stalk[d]['verdict']}"
            assert stalk[d]["tagged_sha256"] == r["protein_sha256"], d
    for d, v in stalk.items():
        assert d in built or v["verdict"] == "fail", f"{d}: passing stalk-site design not built"
        assert d not in built or v["verdict"] == "pass", f"{d}: failing stalk-site design built"


def test_every_exception_to_the_stalk_rule_is_recorded_and_still_an_exception():
    review = _by_design("l1_cterm_review.tsv")
    verdicts = _by_design("l1_cterm_refold_verdicts.tsv")
    stalk = _by_design("l1_cterm_stalk_refold_verdicts.tsv")
    released = arms.released_genes(review.values())
    decisions = _decisions()
    dispo = arms.dispositions()
    couts = arms.cout_dispositions(review.values())
    differs = set()
    for d, r in review.items():
        if r["cterm_tag"] != "kept" or couts[d].action == arms.HOLD:
            continue
        rule = arms.refold_disposition(r, verdicts[d], stalk.get(d), r["gene_symbol"] in released)
        if rule.action != dispo[d].action:
            differs.add(d)
    assert differs == set(arms.RULE_EXCEPTIONS)
    for d in differs:
        rows = [r for r in decisions[d] if "exception to the stalk rule" in r["decision"]]
        assert len(rows) == 1 and rows[0]["decided_by"] == arms.COORDINATOR, d


def test_every_built_snorkel_design_is_one_snorkel_arm():
    # Multi-arm order rule (Becca, 2026-11-14): an internal member may carry loop arms in
    # constructs_l1_internal.tsv beside one snorkel in constructs_l1_cterm.tsv; no other member
    # has arms in two tables.
    import l1_snorkel

    where: dict[str, list[str]] = {}
    for table in ARM_TABLES:
        for r in _tsv(table):
            if (r.get("insert_nt") or "x").strip():
                where.setdefault(r["design_id"], []).append(f"{table}:{r.get('placement', '')}")
    designs = l1_snorkel.built_designs()
    for d in designs:
        snorkels = [w for w in where.get(d, []) if w.startswith("constructs_l1_cterm.tsv:snorkel")]
        assert len(snorkels) == 1, f"{d}: snorkel arms {snorkels}"
    for d, w in where.items():
        tables = {x.split(":")[0] for x in w}
        if len(tables) > 1:
            assert tables == {"constructs_l1_cterm.tsv", "constructs_l1_internal.tsv"}, (d, w)
            assert d in designs, f"{d}: arms in two tables but no built snorkel design"


def test_every_c_out_call_in_question_has_one_decision_row():
    review = _by_design("l1_cterm_review.tsv")
    built = _by_design("constructs_l1_cterm.tsv")
    decisions = _decisions()
    couts = arms.cout_dispositions(review.values())
    weak = {
        d
        for d, r in review.items()
        if set(r["flags"].split(";")) & (set(arms.WEAK_COUT_FLAGS) | {"lane_override_topology"})
    }
    cytoplasmic = {
        d for d, r in review.items() if r["uniprot_cterm_topology"].startswith("Cytoplasmic")
    }
    held = {d for d, x in couts.items() if x.action == arms.HOLD}
    # `held` may be empty: the C-out holds (VTCN1 x2, PLVAP, TRPV1) moved to the internal rung in
    # the 2026-11-15 rescue round. The other two sets prove the review was read.
    assert weak and cytoplasmic
    for d in weak | cytoplasmic | held:
        if review[d]["cterm_tag"] != "kept":
            continue
        rows = [
            r
            for r in decisions.get(d, [])
            if r["question"].startswith(("weak C-out call", "C-out call"))
        ]
        assert len(rows) == 1, f"{d}: {len(rows)} C-out decision rows"
        row = rows[0]
        assert row["decided_by"] == arms.COORDINATOR and row["implemented_by"], d
        assert row["evidence"] and row["rejected_alternatives"] and row["would_change_it"], d
        assert row["decision"].startswith("hold") == (d in held), d
        assert not (d in held and d in built), f"{d}: C-out unconfirmed, but an arm is built"
    # No member UniProt places in the cytoplasm on experimental evidence is built here.
    for d in cytoplasmic:
        if "ECO:0000269" in review[d]["uniprot_cterm_topology"]:
            assert d not in built, d
    # A lane override that makes a member C-terminal is a curated row for this exact protein.
    overrides = _by_design("lane_overrides.tsv")
    for d in weak:
        if "lane_override_topology" in review[d]["flags"].split(";"):
            assert overrides[d]["placement"] == "C-terminal", d
            assert overrides[d]["sequence_sha256"] == review[d]["sequence_sha256"], d


def test_committed_cterm_decisions_reproduce():
    # l1_cterm_decisions.tsv is generated from committed tables (no S3), so a change that moves a
    # member off the rung without regenerating it fails here, in CI, rather than only in
    # `library_pipeline.py verify --stage cterm_review`. #361 moved SLC41A2 to residue 1 and left
    # its C-out hold row behind, 2026-11-15.
    import l1_cterm_decisions

    assert l1_cterm_decisions.main(["--check"]) == 0
