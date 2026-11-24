"""Python mirror of accessible-surfaceome PR #137's canonical/likely gates
(`feat: canonical state-dep escape gates on expression_level, not low_endo`,
branch claude/canonical-expression-gate). Keep byte-faithful to
viewer/lib/catalog-presets.ts @ PR137.

Changes vs PR130 (presets_pr130.py):
  1. Gate on effective_evidence_grade = evidence_grade_summary ?? evidence_grade
     (the synthesizer's holistic grade, which rescues genes the deterministic
     A1-only grade under-calls as `weak`). Requires the record to carry
     `evidence_grade_summary` (fetched by fetch_evidence_summary.py and injected
     into filters); falls back to the deterministic grade when absent.
  2. Canonical evidence FLOOR raised: `weak` and `conflicting` are now held OUT
     of canonical (PR130 admitted `weak` via a fail-closed guard).
  3. Canonical hard-excludes surface_call_reason in {tissue_restricted_surface,
     lysosomal_exocytosis} (non-classical surfacing -> Likely + facet).
  4. Canonical state-dependence escape gates on expression_level in
     {moderate, high} instead of low_endogenous_expression is False (breadth is
     a therapeutic-window property, not a surface-membership signal).
  5. Likely surface_accessibility floor raised to {high, moderate} ('low' dropped
     -> below-Likely `low` tier).
  6. Likely adds a `weak` structural carve-out: a weak-graded gene still passes
     if surface_call_reason is structurally unambiguous (STRUCTURAL_SURFACE_REASONS).
"""

from __future__ import annotations

from typing import Any

# surface_call_reason codes structurally unambiguous enough that a `weak`
# LITERATURE grade reflects thin coverage, not doubt about surface membership.
STRUCTURAL_SURFACE_REASONS: frozenset[str] = frozenset(
    {
        "classical_surface_receptor",
        "multipass_with_exposed_loops",
        "gpi_anchored",
        "tissue_restricted_surface",
    }
)

_FLOOR = ("direct_multi_method", "direct_single_method", "supportive_but_indirect")


def effective_evidence_grade(f: dict[str, Any]):
    """Synthesizer's holistic grade (evidence_grade_summary) if present, else the
    deterministic A1-only evidence_grade. Mirrors effectiveEvidenceGrade()."""
    return f.get("evidence_grade_summary") or f.get("evidence_grade")


def passes_canonical(f: dict[str, Any]) -> bool:
    return (
        f.get("surface_call_reason") != "tissue_restricted_surface"
        and f.get("surface_call_reason") != "lysosomal_exocytosis"
        and effective_evidence_grade(f) in _FLOOR
        and f.get("confidence") in ("high", "moderate")
        and f.get("surface_specificity") in ("surface_dominant", "mixed")
        and (
            f.get("state_dependence") in ("low", "moderate", "unclear")
            or f.get("expression_level") in ("moderate", "high")
        )
        and f.get("surface_accessibility") in ("high", "moderate")
        and f.get("evidence_density") in ("high", "moderate")
    )


def passes_likely(f: dict[str, Any]) -> bool:
    """Full Likely tier (Canonical ⊂ Likely). PR137's "Likely" chip filters to
    passes_likely_only; our analysis uses the full tier for membership labels."""
    g = effective_evidence_grade(f)
    grade_floor_ok = g in _FLOOR
    weak_structural_ok = g == "weak" and f.get("surface_call_reason") in STRUCTURAL_SURFACE_REASONS
    if not (grade_floor_ok or weak_structural_ok):
        return False
    if f.get("surface_specificity") not in ("surface_dominant", "mixed", "mostly_intracellular"):
        return False
    if f.get("surface_accessibility") not in ("high", "moderate"):  # PR137: 'low' dropped
        return False
    sd = f.get("state_dependence")
    if sd is not None and sd not in ("low", "moderate", "high", "unclear"):
        return False
    return True


def passes_likely_only(f: dict[str, Any]) -> bool:
    """Likely MINUS Canonical — the exclusive band the PR137 'Likely' chip shows."""
    return passes_likely(f) and not passes_canonical(f)
