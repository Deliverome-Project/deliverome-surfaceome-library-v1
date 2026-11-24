"""Hard limit on function length across the L1 codebase.

Every function in the exported ``deliverome_analysis`` modules and the core L1 builders,
diversifiers, QC and synthesis scripts must be at most ``MAX_LINES`` (150 lines) with no
exceptions. Legacy functions elsewhere in ``s4e9`` and exported ``s5e13`` are capped by a
monotonic ratchet (``LEGACY_BUDGETS``): no new function may exceed ``MAX_LINES``, and no
existing long function may grow by even a single line.
"""

from __future__ import annotations

import ast
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
REPO = HERE.parents[2]
sys.path.insert(0, str(HERE / "reviewer-package"))

import export_review_repo as err  # noqa: E402

MAX_LINES = 150

#: Core L1 scripts where every function must be <= MAX_LINES with zero exceptions.
STRICT_CORE_SCRIPTS = {
    "build_l1_cterm_arms.py",
    "build_l1_nterm_arms.py",
    "build_l1_multipass_arms.py",
    "build_l1_type3_arms.py",
    "build_l1_internal_arms.py",
    "build_l1_rescue_arms.py",
    "l1_consistency.py",
    "l1_pool_qc.py",
    "synthesis_compliance.py",
    "diversify_cassettes.py",
    "diversify_arm_cassettes.py",
    "diversify_snorkel_modules.py",
}

#: Monotonic ratchet for existing long functions outside the strict core.
#: Never raise a number here; delete an entry once its function is refactored <= MAX_LINES.
LEGACY_BUDGETS: dict[str, int] = {
    "library/l1/orf-design/accessory_factor_dependence.py:main": 220,
    "library/l1/orf-design/apply_reference_fixes.py:reproduce": 159,
    "library/l1/orf-design/audit_compara116_non_atg.py:main": 164,
    "library/l1/orf-design/audit_l1_proteoforms.py:main": 416,
    "library/l1/orf-design/audit_library_completeness.py:audit_row": 174,
    "library/l1/orf-design/audit_library_completeness.py:main": 235,
    "library/l1/orf-design/audit_library_non_atg.py:inspect_members": 165,
    "library/l1/orf-design/audit_non_atg_alternatives.py:main": 331,
    "library/l1/orf-design/audit_ortholog_cds_eligibility.py:evaluate": 180,
    "library/l1/orf-design/audit_ortholog_cds_eligibility.py:main": 211,
    "library/l1/orf-design/build_library_structure_payload.py:main": 222,
    "library/l1/orf-design/coreceptor_host_line.py:main": 422,
    "library/l1/orf-design/draw_eligible_ortholog_panel.py:run": 215,
    "library/l1/orf-design/export_ortholog_eligibility_predictions.py:export": 219,
    "library/l1/orf-design/export_ortholog_predictions.py:export": 170,
    "library/l1/orf-design/l1_completeness.py:build": 252,
    "library/l1/orf-design/l1_swap_pool_model.py:main": 192,
    "library/l1/orf-design/membrane-side-gate/gate_rescue_round.py:main": 165,
    "library/l1/orf-design/nterm_res1_multipass.py:main": 199,
    "library/l1/orf-design/ortholog_architecture_gate.py:run": 284,
    "library/l1/orf-design/plot_order_batch_map.py:main": 163,
    "library/l1/orf-design/plot_rescue_selection.py:main": 179,
    "library/l1/orf-design/plot_vector_maps.py:main": 171,
    "library/l1/orf-design/rebuild_isoform_cohort.py:build": 217,
    "library/l1/orf-design/rebuild_order_maps.py:plans": 202,
    "library/l1/orf-design/redraw_panels.py:main": 212,
    "library/l1/orf-design/review_l1_cterm_members.py:main": 178,
    "library/l1/orf-design/review_ortholog_correspondence.py:run": 305,
    "library/l1/orf-design/review_priority_cds_forms.py:main": 301,
    "library/l1/orf-design/reviewer-package/brief_numbers.py:measure": 443,
    "library/l1/orf-design/reviewer-package/nterm/nterm_decide.py:apply_propeptides": 167,
    "library/l1/orf-design/reviewer-package/nterm/nterm_decide.py:apply_sp_disagree": 185,
    "library/l1/orf-design/reviewer-package/nterm/nterm_decide.py:decide": 359,
    "library/l1/orf-design/reviewer-package/nterm/nterm_placement.py:run": 179,
    "library/l1/orf-design/reviewer-package/signal_peptide_figures.py:draw_examples": 243,
    "library/l1/orf-design/reviewer-package/signal_peptide_figures.py:draw_tree": 179,
    "library/l1/orf-design/run_compara116_adoption.py:main": 165,
    "library/l1/orf-design/score_eligible_ortholog_pool.py:collect_jobs": 189,
    "library/l1/orf-design/score_eligible_ortholog_pool.py:run": 189,
    "library/l1/orf-design/score_ortholog_pool.py:run": 192,
    "library/l1/orf-design/tag_refold/build_refold_inputs.py:build": 444,
    "library/l1/orf-design/tag_refold/build_refold_inputs.py:main": 168,
    "library/l1/orf-design/tag_refold/publish_tagged_structures.py:main": 163,
    "library/l1/orf-design/tag_refold/score_tagged_folds.py:score": 197,
    "library/l1/orf-design/tag_site_model_ablations.py:main": 182,
    "library/l1/orf-design/tag_site_negative_calibration.py:main": 152,
    "library/l1/orf-design/tag_site_residue_scorer.py:main": 179,
    "library/l1/orf-design/tag_site_strict_oof_verify.py:main": 219,
    "library/l1/orf-design/topology_vs_experiment.py:main": 236,
    "library/l1/tag-site-benchmarking/annotation_rule.py:main": 229,
    "library/l1/tag-site-benchmarking/arm_revisions.py:main": 163,
    "library/l1/tag-site-benchmarking/corpus_tier_retrain.py:main": 226,
    "library/l1/tag-site-benchmarking/hard_negative_benchmark.py:main": 160,
    "library/l1/tag-site-benchmarking/l1_internal_arms.py:main": 364,
    "library/l1/tag-site-benchmarking/l1_internal_decisions.py:build": 173,
    "library/l1/tag-site-benchmarking/l1_loop_evidence.py:main": 166,
    "library/l1/tag-site-benchmarking/ranker_hard_negative_layer.py:main": 158,
    "library/l1/tag-site-benchmarking/ranker_pu_segment.py:main": 175,
    "library/l1/tag-site-benchmarking/score_l1_library.py:main": 377,
    "library/l1/tag-site-benchmarking/tag_site_model_v4.py:design": 303,
    "library/l1/tag-site-benchmarking/train_site_models.py:main": 194,
    "library/l1/tag-site-benchmarking/train_site_models.py:make_figures": 172,
}


def _exported_core_files() -> list[str]:
    files = err.collect("HEAD", include_historical=False)
    prefixes = (
        "src/deliverome_analysis/",
        "library/l1/orf-design/",
        "library/l1/tag-site-benchmarking/",
    )
    return [
        p
        for p in sorted(files)
        if p.endswith(".py")
        and p.startswith(prefixes)
        and "/tests/" not in p
        and not Path(p).name.startswith("test_")
        and (REPO / p).exists()
    ]


def test_function_length_hard_limit() -> None:
    strict_failures: list[str] = []
    ratchet_failures: list[str] = []
    observed_long: dict[str, int] = {}

    for rel in _exported_core_files():
        tree = ast.parse((REPO / rel).read_text())
        is_strict = (
            rel.startswith("src/deliverome_analysis/") or Path(rel).name in STRICT_CORE_SCRIPTS
        )
        for node in ast.walk(tree):
            if not isinstance(node, ast.FunctionDef | ast.AsyncFunctionDef):
                continue
            assert node.end_lineno is not None
            length = node.end_lineno - node.lineno + 1
            if length <= MAX_LINES:
                continue
            key = f"{rel}:{node.name}"
            if is_strict:
                strict_failures.append(f"{key} is {length}L (limit {MAX_LINES}L)")
            else:
                observed_long[key] = max(observed_long.get(key, 0), length)
                budget = LEGACY_BUDGETS.get(key)
                if budget is None:
                    ratchet_failures.append(
                        f"NEW long function {key}: {length}L (limit {MAX_LINES}L)"
                    )
                elif length > budget:
                    ratchet_failures.append(f"{key} grew from {budget}L to {length}L")

    assert not strict_failures, "Strict core functions exceed 150 lines:\n  " + "\n  ".join(
        strict_failures
    )
    assert not ratchet_failures, "Function-length limit exceeded:\n  " + "\n  ".join(
        ratchet_failures
    )
    stale = sorted(set(LEGACY_BUDGETS) - set(observed_long))
    assert not stale, (
        "Remove refactored functions from LEGACY_BUDGETS so the ratchet stays tight:\n  "
        + "\n  ".join(stale)
    )
