"""The counts behind the brief's per-rung rationale (section 5), measured from the decision tables.

Each rung's rationale in ``BRIEF.md`` states a rule, the evidence behind it, the rejected
alternative, the measured counts and the known weak spots. This module measures the counts,
so none is typed by hand. It is called by ``brief_numbers.py`` (key ``rungs`` in
``brief_numbers.json``).

Every decision row is assigned to exactly one named class by the rule that produced it. A
row that matches no class, or more than one, stops the run: a new kind of decision must be
named here, and in the brief, before the brief can be regenerated.

Sources (all committed):

* N after SP: ``reviewer-package/nterm/nterm_decisions.tsv`` (``nterm_decide.py``)
* residue 1, type III: ``constructs_l1_type3.tsv`` (``build_l1_type3_arms.design_for``)
* residue 1, multipass: ``constructs_l1_multipass.tsv`` (``build_l1_multipass_arms.design_for``)
* internal and snorkel: S5E13 ``l1_internal_decisions.tsv`` (``l1_internal_decisions.py``),
  ``tag_site_model_v5.json``, ``compare_internal_policies.json``, ``snorkel_spacer_test.json``
  and the saved v5 training table on S3
* C-terminal: ``constructs_l1_cterm.tsv``, ``constructs_l1_cterm_summary.json``,
  ``l1_cterm_decisions.tsv`` (``l1_cterm_decisions.py``)
* every rung: ``l1_tag_paths.tsv``, ``lane_overrides.tsv``, ``fold_order_check.tsv`` and
  ``tag_refold_structures.tsv``

Requested by the coordinator (2026-11-14, Becca's review bar: "clearly explain your decision
making"); implemented by the Claude Code session l1-reviewer-package.
"""

from __future__ import annotations

import csv
import json
import sys
from collections import Counter
from pathlib import Path

from deliverome_analysis.cloning import plasmid_qc as qc

HERE = Path(__file__).resolve().parent
S4E9 = HERE.parent
S5E13 = S4E9.parents[1] / "l1/tag-site-benchmarking"
if str(S4E9) not in sys.path:
    sys.path.insert(0, str(S4E9))


#: The saved v5 training table (``tag_site_model_v5.py``), read for the misnumbered rows.
V5_TRAINING_KEY = "sprint5/s5e13-tag-site-benchmarking/v5/2026-11-14/v5_training_table.tsv"
#: Training proteins whose per-residue features were numbered on another isoform (S5E13
#: README, 2026-11-14, "v5 numbering gap"); retraining is deferred to L2.
V5_MISNUMBERED = ("NRXN1", "ASIC1", "SLC4A11", "SLC9A6")
IGK_MATURE_LEADER = qc.IGK_LEADER_AA


def _tsv(path: Path) -> list[dict]:
    with Path(path).open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def classify(rows: list[dict], classes: dict, what: str) -> dict[str, int]:
    """Count rows per class; each row must match exactly one class predicate."""
    counts = Counter()
    for r in rows:
        hits = [name for name, test in classes.items() if test(r)]
        if len(hits) != 1:
            raise SystemExit(f"{what}: row {r.get('design_id', r)} matches {hits or 'no class'}")
        counts[hits[0]] += 1
    return {name: counts[name] for name in classes}


def _starts(*prefixes: str):
    return lambda r: r["reason"].startswith(prefixes)


# --- N after SP -------------------------------------------------------------------------------
#: Reasons nterm_decide.py writes for members whose signal peptide only DeepTMHMM v1 calls
#: (second-source rule, PR #339; nterm/BRIEF.md section 7).
V1_ONLY = "DeepTMHMM v1 calls a signal peptide"
#: nterm_decide.py's rules (module docstring, "Decision rules"), by the reason each writes.
NTERM_RULES = {
    "native: re-prediction-gated site before the predicted fold (rule 1, Becca 2026-11-11)": lambda r: (
        r["arm"] == "native_site"
        and r["reason"] == "re-prediction-gated site before the predicted fold"
    ),
    "native: only preserved sites lie inside the predicted fold (rule 2)": lambda r: (
        r["arm"] == "native_site"
        and r["reason"].startswith(("only preserved sites", "no disorder track"))
    ),
    "native: cleavage shifts but stays upstream of the cassette (rule 3)": lambda r: (
        r["arm"] == "native_site" and r["reason"].startswith("cleavage shifts")
    ),
    "native: ortholog SignalP release admitted (rule 5)": lambda r: (
        r["arm"] == "native_site" and r["reason"].startswith("ortholog SignalP release")
    ),
    "native: round-2 SignalP on members with no prior call": lambda r: (
        r["arm"] == "native_site" and r["reason"].startswith("round-")
    ),
    "native: moved past a UniProt N-terminal pro-peptide": lambda r: (
        r["arm"] == "native_site" and r["reason"].startswith("cassette after")
    ),
    "native: next SignalP-passing site after a failed refold": lambda r: (
        r["arm"] == "native_site" and r["reason"].startswith("refold of")
    ),
    "igk: only DeepTMHMM v1 calls the signal peptide; a second source confirms it": lambda r: (
        r["arm"] == "igk_swap" and r["reason"].startswith(V1_ONLY) and "confirms" in r["reason"]
    ),
    "igk: only DeepTMHMM v1 calls it; a second source places the cleavage elsewhere": lambda r: (
        r["arm"] == "igk_swap" and r["reason"].startswith(V1_ONLY) and "elsewhere" in r["reason"]
    ),
    "igk: only DeepTMHMM v1 calls it; a signal anchor, its rung move deferred": lambda r: (
        r["arm"] == "igk_swap"
        and r["reason"].startswith(V1_ONLY)
        and "signal anchor" in r["reason"]
    ),
    "igk: native signal peptide fails at every site (rule 4)": lambda r: (
        r["arm"] == "igk_swap" and r["reason"].startswith("native signal peptide fails")
    ),
    "igk: round-2 SignalP, Igκ form keeps the leader's cleavage": lambda r: (
        r["arm"] == "igk_swap" and r["reason"].startswith("round-")
    ),
    "igk: moved past a UniProt N-terminal pro-peptide": lambda r: (
        r["arm"] == "igk_swap" and r["reason"].startswith("cassette after")
    ),
    "igk: next SignalP-passing form after a failed refold": lambda r: (
        r["arm"] == "igk_swap" and r["reason"].startswith("refold of")
    ),
    "no arm: refold fails and no other SignalP-passing form": lambda r: (
        r["arm"] == "no_nterm_arm" and r["reason"].startswith("refold failed")
    ),
    "no arm: no tagged form keeps a SignalP cleavage site": lambda r: (
        r["arm"] == "no_nterm_arm" and r["reason"].startswith("no tagged form")
    ),
    "no arm: the mature N-terminus after a pro-peptide has no passing form": lambda r: (
        r["arm"] == "no_nterm_arm" and r["reason"].startswith("cassette after")
    ),
    "no arm: only DeepTMHMM v1 calls a signal peptide and no second source confirms it": lambda r: (
        r["arm"] == "no_nterm_arm" and r["reason"].startswith(V1_ONLY)
    ),
    "no arm: replacement arm's refold pending": lambda r: (
        r["arm"] == "no_nterm_arm"
        and "refold pending" in r["reason"]
        and not r["reason"].startswith(V1_ONLY)
    ),
    "routed to another rung: a signal anchor or first TM helix, not a signal peptide": lambda r: (
        r["arm"] in ("routed_internal", "routed_multipass") and r["reason"].startswith(V1_ONLY)
    ),
    "igk: weak predicted signal peptide (cleavage-site probability < 0.5; Becca 2026-11-14)": lambda r: (
        r["arm"] == "igk_swap" and r["reason"].startswith("weak predicted SP")
    ),
    "native: refold fails at every SignalP-passing form; built and flagged (disposition rule)": lambda r: (
        r["arm"] == "native_site" and r["reason"].startswith("refold failed at every")
    ),
    "routed to the residue-1 multipass rung": lambda r: (
        r["arm"] == "routed_multipass" and r["reason"].startswith("routed to the residue-1")
    ),
    "no arm: ortholog rule, no human member of the gene is built": lambda r: (
        r["arm"] == "no_nterm_arm" and r["reason"].startswith("ortholog rule")
    ),
    "built and flagged: refold failed at every passing form (flag, never hold; Becca 2026-11-22)": lambda r: (
        r["arm"] in ("native_site", "igk_swap")
        and "built at the first pick, flagged" in r["reason"]
    ),
    "routed back to internal: no form keeps a cleavage site": lambda r: (
        r["arm"] == "routed_internal" and not r["reason"].startswith(V1_ONLY)
    ),
}

#: A flagged build keeps the reason its rule wrote and appends the flag, so it would also match
#: that rule: it is counted once, as flagged.
FLAGGED = (
    "built and flagged: refold failed at every passing form (flag, never hold; Becca 2026-11-22)"
)
NTERM_RULES = {
    name: (test if name == FLAGGED else (lambda r, t=test: t(r) and not NTERM_RULES[FLAGGED](r)))
    for name, test in NTERM_RULES.items()
}


def nterm() -> dict:
    d = _tsv(HERE / "nterm/nterm_decisions.tsv")
    built = [r for r in d if r["arm"] in ("native_site", "igk_swap")]
    return {
        "source": "reviewer-package/nterm/nterm_decisions.tsv",
        "decisions": len(d),
        "arm": dict(Counter(r["arm"] for r in d).most_common()),
        "by_rule": classify(d, NTERM_RULES, "nterm_decisions"),
        "built": len(built),
        "built_inside_predicted_fold": sum(r["inside_predicted_fold"] == "1" for r in built),
        "decided_by_becca_rule": sum(r["decided_by"].startswith("Becca") for r in d),
        # second-source rule (PR #339): swaps whose only confirmation is DeepTMHMM v2, and
        # re-anchored swaps whose new cleavage site rests on a sequence-analysis UniProt feature
        "igk_confirmed_by_deeptmhmm_v2_only": sorted(
            f"{r['gene_symbol']} ({r['candidate_set']})"
            for r in d
            if r["arm"] == "igk_swap"
            and r["reason"].startswith(V1_ONLY)
            and r["reason"].split("does not; ", 1)[-1].startswith("DeepTMHMM v2 SP")
            and "elsewhere" not in r["reason"]
        ),
        "igk_reanchored_on_sequence_analysis": sorted(
            r["gene_symbol"]
            for r in d
            if r["arm"] == "igk_swap"
            and "elsewhere" in r["reason"]
            and "ECO:0000255" in r["reason"]
        ),
        "refold_borderline_counted_as_pass": sorted(
            r["gene_symbol"] for r in d if r["refold_rule"] == "lddt_borderline_counted_as_pass"
        ),
    }


# --- Igκ leader: every rung that can carry one --------------------------------------------------
TYPE3_RULES = {
    "igk: GPI form with no signal peptide": lambda r: (
        r["design"] == "igk_tag" and r["reason"].startswith("zero-TM GPI")
    ),
    "igk: ordered form lacks the gene's signal peptide": lambda r: (
        r["design"] == "igk_tag" and r["reason"].startswith(("ordered form lacks", "IGSF5"))
    ),
    "igk: native N-tail of 20 residues or fewer": lambda r: (
        r["design"] == "igk_tag" and r["reason"].startswith("native N-tail")
    ),
    "igk: weak positive-inside drive (flank-charge difference <= 1)": lambda r: (
        r["design"] == "igk_tag" and r["reason"].startswith("flank-charge")
    ),
    "igk: internal signal-like segment (KLB)": lambda r: (
        r["design"] == "igk_tag" and r["reason"].startswith("internal signal-like")
    ),
    "tag alone: genuine type III with no identified risk": lambda r: r["design"] == "tag_alone",
}

MULTIPASS_RULES = {
    "igk: GPCR default (Becca 2026-11-10)": lambda r: (
        r["design"] == "igk_tag"
        and r["gpcr"] == "1"
        and r["reason"].startswith("GPCR on the multipass")
    ),
    "igk: GPCR reference of a retired expression-rescue member": lambda r: (
        r["design"] == "igk_tag" and r["gpcr"] == "1" and r["reason"].startswith("one arm per gene")
    ),
    "igk: non-GPCR, literature review verdict": lambda r: (
        r["design"] == "igk_tag"
        and r["gpcr"] == "0"
        and r["reason"].startswith("literature review verdict igk")
    ),
    "igk: non-GPCR, built as its reviewed family": lambda r: (
        r["design"] == "igk_tag"
        and r["gpcr"] == "0"
        and r["reason"].startswith("entered with the #283 redraw")
    ),
    "igk: non-GPCR, lane checked and kept (SLC61A1, TSPO2)": lambda r: (
        r["design"] == "igk_tag"
        and r["gpcr"] == "0"
        and not r["reason"].startswith(("literature review verdict", "entered with the #283"))
    ),
    "tag alone: CRHR2 pseudo signal peptide": lambda r: r["design"] == "tag_alone",
}


def igk() -> dict:
    nt = _tsv(HERE / "nterm/nterm_decisions.tsv")
    t3 = _tsv(S4E9 / "constructs_l1_type3.tsv")
    mp = _tsv(S4E9 / "constructs_l1_multipass.tsv")
    folds = _tsv(S4E9 / "fold_order_check.tsv")
    nt_built = [r for r in nt if r["arm"] in ("native_site", "igk_swap")]
    out = {
        "n_after_sp": {
            "igk_swap": sum(r["arm"] == "igk_swap" for r in nt_built),
            "native_site": sum(r["arm"] == "native_site" for r in nt_built),
        },
        "type3": {
            "igk_tag": sum(r["design"] == "igk_tag" for r in t3),
            "tag_alone": sum(r["design"] == "tag_alone" for r in t3),
            "by_rule": classify(t3, TYPE3_RULES, "constructs_l1_type3"),
        },
        "multipass": {
            "igk_tag": sum(r["design"] == "igk_tag" for r in mp),
            "tag_alone": sum(r["design"] == "tag_alone" for r in mp),
            "gpcr": sum(r["gpcr"] == "1" for r in mp),
            "non_gpcr": sum(r["gpcr"] == "0" for r in mp),
            "by_rule": classify(mp, MULTIPASS_RULES, "constructs_l1_multipass"),
        },
        "ordered_inserts_folded_as_igk_mature_chain": sum(
            r["n_terminal_peptide_not_folded"] == IGK_MATURE_LEADER for r in folds
        ),
    }
    # The three N-after-SP site rules compared, and native confidence vs re-prediction
    # (nterm/nterm_build_options.py over the re-predicted members).
    opts = json.loads((HERE / "nterm/nterm_build_options.json").read_text())
    rv = opts["rule_verdicts_over_repredicted_members"]
    conf = opts["native_cleavage_confidence_vs_outcome"]
    out["rule_comparison"] = {
        "repredicted_members": rv["n"],
        "preserved_cleavage_plus_1": rv["cleavage_plus_1"]["preserved"],
        "preserved_disorder_anchored": rv["disorder_anchored"]["preserved"],
        "preserved_reprediction_gated": opts["reprediction_gated"]["members_with_a_pick"],
    }
    out["native_confidence"] = {
        "no_candidate_preserved": sum(b["no_candidate_preserved"] for b in conf.values()),
        "no_candidate_preserved_native_ge_0_9": sum(
            b["no_candidate_preserved"] for k, b in conf.items() if float(k.split("-")[0]) >= 0.9
        ),
        "below_0_5": conf["0.00-0.50"]["members"],
        "below_0_5_with_a_preserved_candidate": conf["0.00-0.50"]["members"]
        - conf["0.00-0.50"]["no_candidate_preserved"],
    }
    out["members_with_igk_leader"] = (
        out["n_after_sp"]["igk_swap"] + out["type3"]["igk_tag"] + out["multipass"]["igk_tag"]
    )
    return out


# --- internal loops and snorkels -----------------------------------------------------------------
SNORKEL_TRIGGERS = {
    "every designed internal site fails refold": _starts(
        "every designed internal site fails refold"
    ),
    "every designed internal site sits on a ligand contact": _starts(
        "every designed internal site has a ligand contact"
    ),
    "no N-after-SP form and the gate rejects every junction": _starts("SignalP"),
    "v1 exposure too small for a loop, gate rejects all of it": lambda r: (
        not r["reason"].startswith(("every designed internal site", "SignalP"))
    ),
}


def internal() -> dict:
    from deliverome_analysis.core.s3 import processed

    d = _tsv(S5E13 / "l1_internal_decisions.tsv")
    model = json.loads((S5E13 / "tag_site_model_v5.json").read_text())
    policies = json.loads((S5E13 / "compare_internal_policies.json").read_text())["mean_over_draws"]
    vs_v4 = _tsv(S5E13 / "l1_internal_sites_vs_v4.tsv")
    training = _tsv(processed.download(V5_TRAINING_KEY))
    mis = [r for r in training if r["gene_symbol"] in V5_MISNUMBERED]
    loops = [r for r in d if r["decision"] == "internal_loop"]

    def metric(arm: str, m: str) -> float:
        return policies[f"all 176 | {arm}"][m]

    return {
        "source": "library/l1/tag-site-benchmarking/l1_internal_decisions.tsv",
        "members": len(d),
        "decision": dict(Counter(r["decision"] for r in d).most_common()),
        "loop_pick_rank": dict(Counter(r["pick_rank"] for r in loops).most_common()),
        "loop_confidence": dict(Counter(r["confidence"] for r in loops).most_common()),
        "loop_refold": dict(Counter(r["refold_verdict"] for r in loops).most_common()),
        "contact_vetoes": sum(bool(r["contact_veto"]) for r in d),
        "literature_pass": sum(
            "literature" in r["contact_sources_available"]
            and "not run" not in r["contact_sources_available"]
            for r in d
        ),
        "literature_not_run": sum("not run" in r["contact_sources_available"] for r in d),
        "vs_v4": dict(Counter(r["vs_v4"] for r in vs_v4).most_common()),
        "vs_v4_sites": len(vs_v4),
        "model": {
            "estimator": model["estimator"],
            "features": model["n_features"],
            "train_rows": model["n_train_rows"],
            "positives": model["n_positives"],
            "hard_negatives_w5": model["n_hard_negatives_in_frame"],
            "negatives_per_positive": model["negatives_per_positive"],
            "cv": "GroupKFold-5 by gene, three negative draws",
            "fitted_at_commit": model["git_commit"][:8],
        },
        "benchmark_mean_over_draws": {
            arm: {
                m: metric(arm, m)
                for m in (
                    "within10_top1",
                    "within10_top3",
                    "segment_first",
                    "contact_top3_any_within2",
                )
            }
            for arm in ("rf_positive_only", "rf_hn_w5 (v5)", "pu_rf", "rf_hn_w5_veto")
        },
        "training_misnumbered": {
            "proteins": sorted({r["gene_symbol"] for r in mis}),
            "proteins_total": len({r["gene_symbol"] for r in training}),
            "rows": len(mis),
            "rows_total": len(training),
            "positives": sum(r["label"] == "1" for r in mis),
            "positives_total": sum(r["label"] == "1" for r in training),
            "hard_negatives_w5": sum(r["weight"] == "5" for r in mis),
            "hard_negatives_w5_total": sum(r["weight"] == "5" for r in training),
        },
    }


def snorkel() -> dict:
    d = _tsv(S5E13 / "l1_internal_decisions.tsv")
    sn = [r for r in d if r["decision"] == "snorkel"]
    held = [r for r in d if r["decision"] == "hold_no_passing_arm"]
    spacer = json.loads((S5E13 / "snorkel_spacer_test.json").read_text())
    built = [
        r for r in _tsv(S4E9 / "constructs_l1_cterm.tsv") if r["placement"].startswith("snorkel")
    ]
    return {
        "decided": len(sn),
        "trigger": classify(sn, SNORKEL_TRIGGERS, "snorkel decisions"),
        # the refold of the ordered snorkel arm, from the refold table of record (the decision
        # table's snorkel_verdict column is not refreshed after later refolds; filter to rung ==
        # "cterm" because snorkel arms live in constructs_l1_cterm.tsv while multi-arm internal
        # members also carry separate internal-loop rows under the same design_id)
        "snorkel_refold": ordered_verdicts({r["design_id"] for r in built}, rung="cterm"),
        "snorkel_verdict_column": dict(
            Counter(r["snorkel_verdict"] or "none" for r in sn).most_common()
        ),
        "built_arms": len(built),
        "held_no_passing_arm": len(held),
        "held_snorkel_refold": dict(Counter(r["snorkel_verdict"] for r in held).most_common()),
        # every snorkel arm refolded in the decision round: the snorkel decisions (all pass) plus
        # the holds whose snorkel arm failed. Not every hold is a snorkel failure: three
        # 2026-11-22 holds (P2RX5, TMEM134, SLC22A9) never had a snorkel arm (blank
        # snorkel_verdict), so counting all holds overstated both numbers (fixed 2026-11-22).
        "snorkel_arms_refolded": len(sn) + sum(r["snorkel_verdict"] == "fail" for r in held),
        "snorkel_arms_failing": sum(r["snorkel_verdict"] == "fail" for r in held),
        "loops_flagged_snorkel_candidate": sum(
            r["snorkel_candidate"] == "True" for r in d if r["decision"] == "internal_loop"
        ),
        "snorkels_flagged_snorkel_candidate": sum(r["snorkel_candidate"] == "True" for r in sn),
        "spacer_test": {
            "change": spacer["change"],
            "failing_tested": len(spacer["failing"]),
            "controls": len(spacer["controls"]),
            "rescued": spacer["result"]["rescued"],
            "controls_still_pass": spacer["result"]["controls_still_pass"],
            "criterion": spacer["works"],
            "criterion_met": spacer["result"]["criterion_met"],
        },
    }


def ordered_verdicts(design_ids: set[str], *, rung: str | None = None) -> dict[str, int]:
    """Refold verdicts of the ordered inserts of these design ids (fold_order_check.tsv)."""
    verdict = {r["fold_id"]: r["verdict"] for r in _tsv(S4E9 / "tag_refold_structures.tsv")}
    return dict(
        Counter(
            verdict.get(r["fold_id"], "no verdict")
            for r in _tsv(S4E9 / "fold_order_check.tsv")
            if r["design_id"] in design_ids and (rung is None or r["rung"] == rung)
        ).most_common()
    )


# --- C-terminal -------------------------------------------------------------------------------
def cterm() -> dict:
    tp = _tsv(S4E9 / "l1_tag_paths.tsv")
    arms = _tsv(S4E9 / "constructs_l1_cterm.tsv")
    summary = json.loads((S4E9 / "constructs_l1_cterm_summary.json").read_text())
    decisions = _tsv(S4E9 / "l1_cterm_decisions.tsv")
    overrides = [r for r in _tsv(S4E9 / "lane_overrides.tsv")]
    on_rung = [r for r in tp if r["placement"] == "C-terminal"]
    placement = Counter(r["placement"] for r in arms)
    return {
        "members": len(on_rung),
        "built_arms_by_placement": dict(placement.most_common()),
        "inserts_ending_in_cassette": sum(
            p in ("C-terminal",) or p.startswith("snorkel") for p in (r["placement"] for r in arms)
        ),
        "stalk_or_furin_sites": sorted(
            f"{r['gene_symbol']} ({r['placement']})"
            for r in arms
            if r["placement"].startswith("internal")
        ),
        "held": sorted(h["gene_symbol"] for h in summary["held"]),
        "decision_questions": dict(Counter(r["question"] for r in decisions).most_common()),
        "gpi_rerouted_to_internal": sorted(
            r["gene_symbol"] for r in tp if r["gpi_rerouted"] == "1"
        ),
        "lane_overrides_to_cterm": sorted(
            r.get("gene_symbol", r.get("design_id", ""))
            for r in overrides
            if "C-terminal" in (r.get("placement", "") + r.get("to", ""))
        ),
    }


def refold_by_rung() -> dict:
    folds = _tsv(S4E9 / "fold_order_check.tsv")
    verdict = {r["fold_id"]: r["verdict"] for r in _tsv(S4E9 / "tag_refold_structures.tsv")}
    out: dict[str, Counter] = {}
    for r in folds:
        out.setdefault(r["rung"], Counter())[verdict.get(r["fold_id"], "no verdict")] += 1
    return {k: dict(v) for k, v in sorted(out.items())}


def measure() -> dict:
    return {
        "n_after_sp": nterm(),
        "igk": igk(),
        "internal": internal(),
        "snorkel": snorkel(),
        "cterm": cterm(),
        "refold_by_rung": refold_by_rung(),
    }


if __name__ == "__main__":
    print(json.dumps(measure(), indent=1))
