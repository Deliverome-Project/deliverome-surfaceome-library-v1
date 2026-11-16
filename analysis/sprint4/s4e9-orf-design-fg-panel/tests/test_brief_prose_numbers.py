"""The numbers the reviewer brief states in prose equal the ones measured from the tables.

``reviewer-package/brief_numbers.py`` measures every count into ``brief_numbers.json``; the
brief's prose is written by hand. Each claim below is a fragment of the brief with ``{}``
where a measured number goes, and the path to that number in ``brief_numbers.json``. The test
fails when the brief no longer contains the fragment with the measured value, so a table
change that moves a number cannot leave the prose behind. Offline: it reads only committed
files.

Requested by the coordinator after the 2026-11-14 independent audit ("add a test that checks
prose numbers against tables"); implemented by the Claude Code session l1-reviewer-package.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

RP = Path(__file__).resolve().parents[1] / "reviewer-package"
BRIEF = (RP / "BRIEF.md").read_text()
NUMBERS = json.loads((RP / "brief_numbers.json").read_text())

#: (fragment with {} placeholders, paths into brief_numbers.json, one per placeholder)
CLAIMS = [
    # section 5 rungs
    (
        "| N after SP | the residue after the signal peptide is extracellular | predicted cleavage site | {} |",
        ["tag_paths.placement.N after SP"],
    ),
    (
        "| C-terminal | C-terminus extracellular | topology only | {} |",
        ["tag_paths.placement.C-terminal"],
    ),
    ("needs a loop or ectodomain site | {} |", ["tag_paths.placement.internal"]),
    (
        "{} are built ({} native, {} Igκ) and {} have no N-terminal arm",
        [
            "rungs.n_after_sp.built",
            "rungs.n_after_sp.arm.native_site",
            "rungs.n_after_sp.arm.igk_swap",
            "rungs.n_after_sp.arm.no_nterm_arm",
        ],
    ),
    ("otherwise. {} members carry it.", ["rungs.igk.members_with_igk_leader"]),
    (
        "Over the {} re-predicted members, {} have no",
        [
            "rungs.igk.rule_comparison.repredicted_members",
            "rungs.igk.native_confidence.no_candidate_preserved",
        ],
    ),
    (
        "candidate that keeps cleavage, and {} of those",
        ["rungs.igk.native_confidence.no_candidate_preserved_native_ge_0_9"],
    ),
    (
        "keeps cleavage for {} members, the most\ndisordered residue in the first 20 for {}, and the re-prediction-gated choice for {}.",
        [
            "rungs.igk.rule_comparison.preserved_cleavage_plus_1",
            "rungs.igk.rule_comparison.preserved_disorder_anchored",
            "rungs.igk.rule_comparison.preserved_reprediction_gated",
        ],
    ),
    # internal and snorkel
    ("**Counts.** {} members are tagged in a loop", ["rungs.internal.decision.internal_loop"]),
    (
        "It is set on {} of the {} loop members",
        ["rungs.snorkel.loops_flagged_snorkel_candidate", "rungs.internal.decision.internal_loop"],
    ),
    ("and on only {} of the 28 snorkels", ["rungs.snorkel.snorkels_flagged_snorkel_candidate"]),
    (
        # The decision round's denominator, which is not the order's: the order carries
        # rungs.snorkel.built_arms snorkel arms, pinned by the claim below.
        "{} of the {} arms refolded then failed",
        ["rungs.snorkel.snorkel_arms_failing", "rungs.snorkel.snorkel_arms_refolded"],
    ),
    (
        "The order carries {} snorkel arms, {} of which\npass",
        [
            "rungs.snorkel.built_arms",
            "refold.by_group.cterm: snorkel (TM module, then cassette).pass",
        ],
    ),
    (
        "Against v4, of the {} sites\ndesigned, {} are identical",
        ["rungs.internal.vs_v4_sites", "rungs.internal.vs_v4.identical"],
    ),
    (
        "are identical, {} within 10 residues, {} further away and {} on sequences v4",
        [
            "rungs.internal.vs_v4.within_10",
            "rungs.internal.vs_v4.further",
            "rungs.internal.vs_v4.no_v4_pick",
        ],
    ),
    (
        "5 of 176 positives and {} of {} rows",
        [
            "rungs.internal.training_misnumbered.rows",
            "rungs.internal.training_misnumbered.rows_total",
        ],
    ),
    # section 12 and the order
    ("| **All with a tag path** | **2,710** | **{}** |", ["order_plan.members"]),
    ("**{} pass, {} fail (", ["refold.pass", "refold.fail"]),
    ("it translates each of the {} inserts", ["refold.inserts"]),
    (
        "All {} match residue for residue, {} of them as the mature chain",
        ["refold.fold_check_ok", "refold.ordered_igk_leader_inserts"],
    ),
    (
        "The {} built members fill {} plates per set",
        ["order_plan.members", "order_plan.plates_per_set"],
    ),
    # swap risk
    (
        "| one shared cassette (p003 codons) | {} / 2,529 nt | {} |",
        [
            "swap_risk.schemes.shared.cterm_worst_partner_median",
            "swap_risk.schemes.shared.members_over_gate_100",
        ],
    ),
    ("For the {} members whose\ncassette ends the ORF", ["swap_risk.schemes.shared.cterm_members"]),
    # structure models and GPI
    (
        "| identical | {} | the model is of exactly this protein |",
        ["structure_match.verdict_active.identical"],
    ),
    (
        "| different | {} | the model is another proteoform",
        ["structure_match.verdict_active.different"],
    ),
    (
        "rule: {} genes, of which {} are NetGPI-anchored",
        ["gpi.gpi_set_genes", "gpi.basis.netgpi_anchored"],
    ),
]


def get(path: str):
    node = NUMBERS
    for part in path.split("."):
        node = node[part]
    return node


def fmt(v) -> str:
    if isinstance(v, float) and v.is_integer():
        v = int(v)
    return f"{v:,}" if isinstance(v, int) else str(v)


@pytest.mark.parametrize(("fragment", "paths"), CLAIMS, ids=[c[0][:40] for c in CLAIMS])
def test_brief_states_the_measured_number(fragment, paths):
    expected = fragment.format(*(fmt(get(p)) for p in paths))
    assert expected in BRIEF, f"brief should say: {expected!r}"
