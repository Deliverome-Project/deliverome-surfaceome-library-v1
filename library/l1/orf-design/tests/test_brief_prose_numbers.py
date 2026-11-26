# ruff: noqa: RUF001
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
    # section 7 placements
    (
        "| **1. N after SP** | Extracellular N-terminus with a cleavable signal peptide (includes 186 of 191 GPI CDSs) | After SP cleavage, verified by SignalP 6: native SP (`1,091`) or Igκ swap (`230`) | **{}** (1,321 / 1) |",
        ["tag_paths.placement.N after SP"],
    ),
    (
        "| **4. C-terminal** | Extracellular C-terminus, cytoplasmic N-terminus | Cassette at native C-terminus (`214`), or moved to stalk/pre-furin site (`3`) | **{}** (217 / 2) |",
        ["tag_paths.placement.C-terminal"],
    ),
    (
        "| **5. Internal loops (+ snorkels)** | Neither terminus freely extracellular (both ends cytoplasmic, C-GPI `BST2`, or terminal fallback) | 1–3 ranked extracellular loop sites (`1,066` arms) + C-terminal snorkel companions (`141`) | **{}** (585 / 7) |",
        ["tag_paths.placement.internal"],
    ),
    (
        # Reworded 2026-11-22: the brief now counts the decisions, not the built ones; built is
        # pinned by the tag-path table claim below (rungs.n_after_sp.built).
        "Of the {} CDSs decided here, {} keep their own signal\npeptide, {} have it replaced by "
        "the Igκ leader, {} gets no N-terminal construct",
        [
            "rungs.n_after_sp.decisions",
            "rungs.n_after_sp.arm.native_site",
            "rungs.n_after_sp.arm.igk_swap",
            "rungs.n_after_sp.arm.no_nterm_arm",
        ],
    ),
    ("| N after SP | {} | {} |", ["tag_paths.placement.N after SP", "rungs.n_after_sp.built"]),
    ("and never otherwise, and {} CDSs carry it:", ["rungs.igk.members_with_igk_leader"]),
    (
        "Over the {} re-predicted CDSs, {} have no",
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
        "keeps cleavage for {} CDSs, the most\ndisordered residue in the first 20 for {}, and the re-prediction-gated choice for {}.",
        [
            "rungs.igk.rule_comparison.preserved_cleavage_plus_1",
            "rungs.igk.rule_comparison.preserved_disorder_anchored",
            "rungs.igk.rule_comparison.preserved_reprediction_gated",
        ],
    ),
    # internal and snorkel
    (
        "The order carries {} snorkel constructs, {} of which\npass",
        [
            "rungs.snorkel.built_arms",
            "refold.by_group.cterm: snorkel (TM module, then cassette).pass",
        ],
    ),
    # section 8 and the order
    ("| **All with a tag path** | **{}** | **{}** |", ["tag_paths.members", "order_plan.members"]),
    ("**{} constructs: {} pass, {} fail, ", ["refold.inserts", "refold.pass", "refold.fail"]),
    ("translates each of the {} constructs", ["refold.inserts"]),
    (
        "All {} folded constructs match residue for residue, {} of them as the mature chain",
        ["refold.fold_check_ok", "refold.ordered_igk_leader_inserts"],
    ),
    (
        "The **{}** built CDSs fill **{}** plates",
        ["order_plan.members", "order_plan.plates_per_set"],
    ),
    # swap risk
    (
        "All **{}** CDSs exceeding our 100-nt shared-tract gate",
        ["swap_risk.schemes.linkers_only.members_over_gate_100"],
    ),
    # structure models and GPI
    (
        "| identical | {} | the model is of exactly this protein |",
        ["structure_match.verdict_active.identical"],
    ),
    (
        "| different | {} | the model is another form",
        ["structure_match.verdict_active.different"],
    ),
    ("| NetGPI calls a form anchored | {} | yes |", ["gpi.basis.netgpi_anchored"]),
    ("That gives {} GPI genes with {} L1 CDSs.", ["gpi.gpi_set_genes", "gpi.l1_gpi_members"]),
    # swap impurity by construct length (section 3)
    (
        "the median is {}% of a barcode's cells in the longest\nquarter of constructs and {}% in the shortest (Spearman {} against insert length). The {}\nbarcodes above 1%",
        [
            "swap_pool.by_length.longest_quarter_specific_median_pct",
            "swap_pool.by_length.shortest_quarter_specific_median_pct",
            "swap_pool.by_length.spearman_length_vs_specific",
            "swap_pool.by_length.barcodes_over_1pct",
        ],
    ),
    ("highest\nis {}%)", ["swap_pool.by_length.max_specific_pct"]),
    # proteoform audit and canonical proteins with no Ensembl transcript (section 5)
    (
        "(including {} genes where no Ensembl 115 transcript encodes the UniProt canonical\n  protein — {} same-length allele differences and {} length differences;",
        [
            "reference.canonical_without_transcript.genes",
            "reference.canonical_without_transcript.same_length_residues_differ",
            "reference.canonical_without_transcript.length_differs",
        ],
    ),
    (
        "{} of the\n  188 non-canonical picks tripped a structural/SP/TM/Pfam check and {} carried an unflagged > 10%",
        [
            "proteoform_audit.flagged_major",
            "proteoform_audit.unflagged_over_10pct",
        ],
    ),
    # The brief's PeptideAtlas evidence is now the 2026-11-22 review's full-accession query,
    # not the audit's pinned-pool query (which covered only some genes).
    (
        "querying PeptideAtlas build 607 across all {} reviewed UniProt accessions of {} differing genes\n"
        "  (`{}` peptides only for our form, `{}` only for canonical, `{}` both, `{}` neither)",
        [
            "proteoform_review.accessions_queried",
            "proteoform_review.reviewed",
            "proteoform_review.evidence_class.ours_supported",
            "proteoform_review.evidence_class.canonical_supported",
            "proteoform_review.evidence_class.both",
            "proteoform_review.evidence_class.neither",
        ],
    ),
    # isoform panel as drawn (section 5); refilled by the redraw
    (
        "| Isoforms | {} ({} genes) | 100 observed + 50 prediction-only | type I {}, multipass "
        "N-terminus cytoplasmic {}, multipass N-terminus extracellular {}, type II {}, GPI {}, "
        "type III {} |",
        [
            "isoform_panel.size",
            "isoform_panel.genes",
            "isoform_panel.by_class.type I",
            "isoform_panel.by_class.multipass N-in",
            "isoform_panel.by_class.multipass N-out",
            "isoform_panel.by_class.type II",
            "isoform_panel.by_class.GPI",
            "isoform_panel.by_class.type III",
        ],
    ),
    (
        "({} of the {} picks in the panel draw were replaced",
        ["isoform_panel.replaced", "isoform_panel.size"],
    ),
    # The class table since the 2026-11-22 redraw: multipass is two strata, the cohort and gate
    # are the rebuilt ones, and the panel column splits observed + prediction-only.
    *[
        (
            f"| {label} | {{}} | {{}} | {{}} | {{}} ({{}} + {{}}) | {{}} | {{}} |",
            [
                f"isoform_panel.classes.{c}.{k}"
                for k in (
                    "library_genes",
                    "cohort_forms",
                    "drawable_forms_gate_of_record",
                    "panel_forms",
                    "panel_forms_observed",
                    "panel_forms_prediction_only",
                    "uniform_expected_from_cohort",
                    "enrichment_vs_uniform_cohort",
                )
            ],
        )
        for label, c in (
            ("type I", "type I"),
            ("multipass, N-terminus extracellular", "multipass N-out"),
            ("multipass, N-terminus cytoplasmic", "multipass N-in"),
            ("type II", "type II"),
            ("GPI", "GPI"),
            ("type III", "type III"),
        )
    ],
    # native vs Igκ signal peptides (section 6)
    (
        "Of the CDSs tagged after a signal peptide, {} keep their own and\n{} have it replaced "
        "by the Igκ leader; at residue 1, {} multipass and {} type III CDSs",
        [
            "rungs.igk.n_after_sp.native_site",
            "rungs.igk.n_after_sp.igk_swap",
            "rungs.igk.multipass.igk_tag",
            "rungs.igk.type3.igk_tag",
        ],
    ),
    # section 3 swap risk
    (
        "Diversifying the four `GGGGS` linkers, as ordered, lowers\nthat to {}%",
        ["swap_pool.total_pct"],
    ),
    (
        "| the tag cassette (wrong-protein fusions) | every pair | 4.6% | {}% (the two tags only) |",
        ["swap_pool.shared_tag_chimera_pct"],
    ),
    (
        "the median is {}% of a barcode's cells and\n0.30% of the library ({}% chimeras, {}% mislabels)",
        [
            "swap_pool.specific_median_per_barcode_pct",
            "swap_pool.specific_chimera_pct",
            "swap_pool.specific_mislabel_pct",
        ],
    ),
    # section 9 cassette uniqueness and snorkel count
    (
        "all **{} cassettes are unique**",
        ["order_plan.inserts"],
    ),
    (
        "Each of the {} snorkel constructs",
        ["rungs.snorkel.built_arms"],
    ),
    # section 10 order plan & pooling
    (
        "(at 3.5×/kb: **{} µg** per pool at 100 ng/µL, **{} µg** at 75 ng/µL,\n**{} µg** at 50 ng/µL)",
        [
            ("order_plan", "concentration_scenarios", "100_ng_ul", "pool_ug_total_per_pool_x3.5"),
            ("order_plan", "concentration_scenarios", "75_ng_ul", "pool_ug_total_per_pool_x3.5"),
            ("order_plan", "concentration_scenarios", "50_ng_ul", "pool_ug_total_per_pool_x3.5"),
        ],
    ),
    (
        "Only the **{}** constructs over 4.3 kb exceed a\n  single 5-µL pooling well",
        ["order_plan.inserts_over_threshold"],
    ),
    (
        "| 2.0× | 13× / 18× | 31× / 50× | {} µg / 54.1 µg |",
        [("order_plan", "weighting", "x2.0", "extra_ug_assumed_total")],
    ),
    (
        "| 2.5× | 29× / 41× | 94× / 152× | {} µg / 74.9 µg |",
        [("order_plan", "weighting", "x2.5", "extra_ug_assumed_total")],
    ),
    (
        "| 3.5× | 101× / 141× | 495× / 801× | {} µg / 112.4 µg |",
        [("order_plan", "weighting", "x3.5", "extra_ug_assumed_total")],
    ),
    # section 11 reproducibility order count
    (
        "**Order of record (`2026-11-25-pre-freeze-e`).** **{}** constructs across **2,730** built CDSs",
        ["order_plan.inserts"],
    ),
]


def get(path: str | tuple[str, ...]):
    parts = path.split(".") if isinstance(path, str) else path
    node = NUMBERS
    for part in parts:
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
