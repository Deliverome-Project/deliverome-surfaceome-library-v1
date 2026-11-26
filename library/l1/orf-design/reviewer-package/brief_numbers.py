"""Every number in the L1 design brief, measured from the library's tables, plus its figures.

The brief (``BRIEF.md``) quotes counts; this script is where each one comes from, so a
count in the brief is never typed from memory. It reads the committed s4e9 tables and the
published S3 records only, so it runs the same against S3 or against the reviewer bundle
(``DELIVEROME_S3_MIRROR``, set by ``run_review.py``).

Counts are labeled with the tranche (L1), the cohort (active members, i.e. not retired as
a duplicate) and the topology tool (DeepTMHMM v1) they refer to.

    uv run --locked python brief_numbers.py            # numbers -> brief_numbers.json, figures
    uv run --locked python brief_numbers.py --no-figures
"""

# ruff: noqa: RUF001  (typographic multiplication sign, alpha, and en-dashes in figure labels are intentional)
from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter
from pathlib import Path

HERE = Path(__file__).parent
S4E9 = HERE.parent
PREFIX = "sprint4/s4e9-orf-design-fg-panel/"
OUT = HERE / "brief_numbers.json"
FIGURES = HERE / "figures"
#: The audit of the S5E13 literature positives (PR 236) and the corrected training corpus.
AUDIT_KEY = "scratch/becca/tag-site-corpus-audit/audit.tsv"
CORRECTED_KEY = "sprint5/s5e13-tag-site-benchmarking/corpus_audit_corpus.tsv"
GOLD_KEY = "scratch/becca/tag-site-corpus-audit/gold_standard.tsv"
AUDIT_DIR = "scratch/becca/tag-site-corpus-audit/"

PLACEMENTS = [
    "N after SP",
    "N at residue 1, multipass",
    "N at residue 1, single-pass (type III)",
    "C-terminal",
    "internal",
]
SETS = ["library_reference", "isoform", "ortholog"]


def _tsv(path: Path) -> list[dict]:
    with Path(path).open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def reference_selection() -> list[dict]:
    from deliverome_analysis.s3 import processed

    # The current record, read through the pipeline's own constant so the two cannot diverge.
    # (Until 2026-11-13 this named the round-1 record, "2026-11-12-reference-fixes/", and so
    # missed round 2; found by l1-proteoform-audit, PR #305.)
    sys.path.insert(0, str(S4E9))
    import library_pipeline as lp

    return _tsv(processed.download(PREFIX + lp.REFERENCE_RECORD + "reference_selection.tsv"))


def _ortholog_panel() -> dict:
    """The ortholog panel L1 orders, from the adopted CDS-policy draw."""
    # The ortholog panel L1 orders is the adopted CDS-policy draw, not the retired Compara 115
    # panel that ortholog_panel_summary.json describes (fixed 2026-11-14, l1-reviewer-package).
    final = json.loads((S4E9 / "ortholog_cds_final_panel_summary.json").read_text())
    adopted = _tsv(S4E9 / "ortholog_panel_adopted.tsv")
    return {
        "source": "ortholog_cds_final_panel_summary.json + ortholog_panel_adopted.tsv",
        "seed": final["policy"]["seed"],
        "cutoffs": final["policy"]["cutoffs"],
        "selected": len(adopted),
        "by_species": dict(Counter(r["species"] for r in adopted).most_common()),
        "by_topology_class": dict(Counter(r["topology_class"] for r in adopted).most_common()),
        "replaced_by_review": sum(bool(r["replaces"]) for r in adopted),
    }


def _gene_list() -> dict:
    """The frozen root gene list with its dated additions, withdrawals and exclusions."""
    root = json.loads(
        (S4E9.parents[2] / "library/l1/tag-site-architecture/library_members.json").read_text()
    )
    add_rows = _tsv(S4E9 / "gene_list_additions.tsv")
    added = [r["gene"] for r in add_rows if not (r.get("withdrawn") or "").strip()]
    withdrawn = [r["gene"] for r in add_rows if (r.get("withdrawn") or "").strip()]
    excluded = [r["gene"] for r in _tsv(S4E9 / "gene_list_exclusions.tsv")]
    excluded_set = set(excluded)
    excluded_from_root = sorted(g for g in root if g in excluded_set)
    active_root = [g for g in root if g not in excluded_set]
    return {
        "frozen_root_genes": len(root),
        "excluded_from_root": excluded_from_root,
        "active_root_genes": len(active_root),
        "added": added,
        "withdrawn_additions": withdrawn,
        "excluded_from_filter_output": excluded,
        "library_genes": len(active_root) + len(added),
    }


def _panel_pool() -> dict:
    """The pool the panels were drawn from (check_universe_coverage.py, verify stage 5b)."""
    pool = json.loads((S4E9 / "universe_coverage_summary.json").read_text())
    return {
        "library_genes": pool["library_genes"],
        "library_genes_in_april_set": pool["library_genes_in_april_set"],
        "mouse_eligible": pool["mouse"]["eligible_in_compara115"],
        "mouse_listed": pool["mouse"]["in_ortholog_list"],
        "mouse_never_queried": pool["mouse"]["eligible_not_listed"],
        "cyno_eligible": pool["cyno"]["eligible_in_compara115"],
        "cyno_listed": pool["cyno"]["in_ortholog_list"],
        "cyno_never_queried": pool["cyno"]["eligible_not_listed"],
        "isoform_genes_not_in_cohort": pool["isoform"]["isoform_gene_not_in_cohort"],
        "isoform_genes_outside_april_set": pool["isoform"]["of_which_outside_april_set"],
        "recorded_pairs_checked": pool["recorded_pairs_checked"],
    }


def measure() -> dict:
    n: dict = {}

    # --- gene list -----------------------------------------------------------------
    n["gene_list"] = _gene_list()

    # --- reference proteoform ladder (one row per library gene) ------------------------
    ref = reference_selection()
    decisions = Counter(r["decision"] for r in ref)
    disagree = {
        "resolved_mane",
        "weak_evidence_mane",
        "weak_evidence_appris",
        "unresolved_default_mane",
    }
    n["reference"] = {
        "genes": len(ref),
        "decision": dict(decisions.most_common()),
        "mane_appris_disagreements": sum(decisions[d] for d in disagree),
        "final_source": dict(Counter(r["selected_source"] for r in ref).most_common()),
        "decision_by_source": {
            d: dict(Counter(r["selected_source"] for r in ref if r["decision"] == d))
            for d in decisions
        },
        "needs_review": sum(r["needs_review"] == "True" for r in ref),
        "changes_uniprot_canonical": sum(r["changes_current_reference"] == "True" for r in ref),
        "named_overrides": sorted(r["symbol"] for r in ref if r["selected_source"] == "named"),
    }
    # Genes whose UniProt canonical protein no Ensembl 115 transcript translates to.
    no_tx = [r for r in ref if r["current_has_ensembl_transcript"] == "False"]
    n["reference"]["canonical_without_transcript"] = {
        "genes": len(no_tx),
        "same_length_residues_differ": sum(
            r["current_minus_selected"] in ("0", "0.0") for r in no_tx
        ),
        "length_differs": sum(r["current_minus_selected"] not in ("0", "0.0") for r in no_tx),
        "symbols": sorted(r["symbol"] for r in no_tx),
    }

    # --- tranches ------------------------------------------------------------------
    manifest = _tsv(S4E9 / "candidate_manifest_tranches.tsv")
    l1 = [r for r in manifest if r["library"] == "l1"]
    active = [r for r in l1 if not r["retired_reason"]]
    n["tranches"] = {
        "l1_rows": len(l1),
        "l1_retired_duplicates": len(l1) - len(active),
        "l1_active": len(active),
        "l1_active_by_set": dict(Counter(r["set"] for r in active).most_common()),
        "l2_active": sum(r["library"] == "l2" and not r["retired_reason"] for r in manifest),
        "all_rows_with_orf": sum(bool(r["orf_nt"]) for r in manifest),
        "l1_active_orf_nt_median": sorted(int(r["orf_nt"]) for r in active)[len(active) // 2],
    }
    # Library genes with no active reference member, and whether another member encodes them.
    built = {r["key"].split(":", 1)[1] for r in active if r["set"] == "library_reference"}
    also = {s for r in active for s in (r["also_encodes"] or "").split(";") if s}
    unbuilt = sorted({r["symbol"] for r in ref} - built)
    n["tranches"]["genes_without_reference_member"] = {
        "encoded_by_another_member": [g for g in unbuilt if g in also],
        "no_designed_cds": [g for g in unbuilt if g not in also],
    }

    # Domestication: synonymous edits that clear BbsI / Esp3I / BsaI from each ORF.
    edits = sorted(int(float(r["edits"])) for r in active if r["edits"])
    n["domestication"] = {
        "members_with_edit_count": len(edits),
        "edits_median": edits[len(edits) // 2],
        "edits_max": edits[-1],
        "members_needing_no_edit": sum(e == 0 for e in edits),
        "esp3i_sites_per_vector": dict(Counter(r["esp3i_sites"] for r in active)),
        "vector_bp_min_median_max": [
            min(int(r["vector_bp"]) for r in active),
            sorted(int(r["vector_bp"]) for r in active)[len(active) // 2],
            max(int(r["vector_bp"]) for r in active),
        ],
        # The manifest's has_igk_leader describes the untagged draft (only the retired rescue
        # arms had one there); the Igk leader of the ordered N-terminal and multipass arms is
        # counted from the order of record under "refold" ("ordered_igk_leader_inserts").
        "untagged_draft_igk_leader_members": sum(r["has_igk_leader"] == "True" for r in active),
    }
    # The pool the panels were drawn from (check_universe_coverage.py, verify stage 5b).
    n["panel_pool"] = _panel_pool()
    n["tranches"]["isoform_panel"] = json.loads((S4E9 / "isoform_panel_summary.json").read_text())
    n["tranches"]["ortholog_panel"] = _ortholog_panel()

    # --- tag paths (L1 active members with a designed CDS; DeepTMHMM v1) ----------------
    tp = _tsv(S4E9 / "l1_tag_paths.tsv")
    n["tag_paths"] = {
        "tool": "deeptmhmm_v1",
        "members": len(tp),
        "placement": {p: sum(r["placement"] == p for r in tp) for p in PLACEMENTS},
        "placement_by_set": {
            p: {s: sum(r["placement"] == p and r["candidate_set"] == s for r in tp) for s in SETS}
            for p in PLACEMENTS
        },
        "gpi_rerouted": sorted(r["gene_symbol"] for r in tp if r["gpi_rerouted"] == "1"),
        "topology_version": dict(Counter(r["topology_version"] for r in tp).most_common()),
        "topology_cohort": dict(Counter(r["topology_cohort"] for r in tp).most_common()),
        "internal_without_inter_tm_loop_ge_8": sorted(
            r["gene_symbol"]
            for r in tp
            if r["placement"] == "internal" and int(r["longest_inter_tm_ecto_run_aa"]) < 8
        ),
        "internal_inter_tm_run_median": sorted(
            int(r["longest_inter_tm_ecto_run_aa"]) for r in tp if r["placement"] == "internal"
        )[sum(r["placement"] == "internal" for r in tp) // 2],
    }

    # --- GPI ------------------------------------------------------------------------
    gpi = _tsv(S4E9 / "gpi_set.tsv")
    gpi_genes = {r["gene_symbol"] for r in gpi if r["gpi_anchored"] == "1"}
    gpi_members = [r for r in tp if r["gene_symbol"] in gpi_genes]
    n["gpi"] = {
        "gpi_set_genes": len(gpi_genes),
        "basis": dict(Counter(r["basis"] for r in gpi).most_common()),
        "conflict_genes": sorted(
            r["gene_symbol"] for r in gpi if r["basis"].startswith("conflict")
        ),
        "l1_gpi_members": len(gpi_members),
        "l1_gpi_members_by_rung_before_reroute": dict(
            Counter(r["rung_before_gpi"] for r in gpi_members).most_common()
        ),
        "l1_conflict_members": {
            r["gene_symbol"]: r["placement"]
            for r in gpi_members
            if r["gene_symbol"]
            in {g["gene_symbol"] for g in gpi if g["basis"].startswith("conflict")}
        },
    }

    # --- structure models vs the ordered protein ---------------------------------------
    pm = _tsv(S4E9 / "proteoform_match.tsv")
    active_ids = {r["candidate_id"] for r in active}
    pma = [r for r in pm if r["candidate_id"] in active_ids]
    n["structure_match"] = {
        "compared_rows": len(pm),
        "verdict_all_rows": dict(Counter(r["verdict"] for r in pm).most_common()),
        "verdict_active": dict(Counter(r["verdict"] for r in pma).most_common()),
        "different_active_by_set": dict(
            Counter(r["set"] for r in pma if r["verdict"] == "different").most_common()
        ),
        "different_same_length_active": sum(
            r["verdict"] == "different" and r["orf_aa"] == r["model_aa"] for r in pma
        ),
    }

    # --- refold verdicts of the ordered inserts (section 12) -----------------------------
    # One verdict per insert in the order of record: fold_order_check.tsv names the fold of the
    # exact ordered protein, tag_refold_structures.tsv carries its verdict. The C-terminal
    # build file holds C-terminal, snorkel and stalk arms, split by their recorded placement.
    folds = _tsv(S4E9 / "fold_order_check.tsv")
    verdict = {r["fold_id"]: r["verdict"] for r in _tsv(S4E9 / "tag_refold_structures.tsv")}
    cterm_place = {r["design_id"]: r["placement"] for r in _tsv(S4E9 / "constructs_l1_cterm.tsv")}
    by: dict[str, Counter] = {}
    for r in folds:
        group = r["rung"]
        if group == "cterm":
            group = f"cterm: {cterm_place.get(r['design_id'], '?')}"
        by.setdefault(group, Counter())[verdict.get(r["fold_id"], "no verdict")] += 1
    total = Counter()
    for c in by.values():
        total.update(c)
    n["refold"] = {
        "inserts": len(folds),
        "fold_check_ok": sum(r["status"] == "ok" for r in folds),
        "ordered_igk_leader_inserts": sum(
            r["n_terminal_peptide_not_folded"] == "METDTLLLWVLLLWVPGSTG" for r in folds
        ),
        "by_group": {g: dict(c) for g, c in sorted(by.items())},
        "pass": total["pass"],
        "fail": total["fail"],
        "pass_pct": round(100 * total["pass"] / len(folds), 1) if folds else None,
    }

    # --- per-rung rationale (section 5): rule classes, Igκ use, ranker, snorkel, C-terminal
    sys.path.insert(0, str(HERE))
    import rung_numbers

    n["rungs"] = rung_numbers.measure()

    # --- the ordering plan (order_plan.py)
    n["order_plan"] = json.loads((HERE / "l1_order_plan_summary.json").read_text())
    # Inserts per CDS in the order (an internally tagged CDS may be ordered at 2 or 3 sites).
    per_cds = Counter(r["design_id"] for r in _tsv(S4E9 / "l1_insert_numbering.tsv"))
    n["order_plan"]["cds_with_construct"] = len(per_cds)
    n["order_plan"]["cds_with_2_or_3_constructs"] = sum(v > 1 for v in per_cds.values())
    n["order_plan"]["constructs_of_cds_with_2_or_3"] = sum(v for v in per_cds.values() if v > 1)

    # --- the proteoform audit (audit_l1_proteoforms.py) on the reference CDSs that differ from the
    # UniProt canonical protein: what its flags caught, and what PeptideAtlas could test.
    audit_rows = [
        r
        for r in _tsv(S4E9 / "l1_proteoform_audit.tsv")
        if r["candidate_set"] == "library_reference"
        and r["canonical_aa"]
        and "UniProt_canonical" not in r["member_equals"]
    ]

    def _gap(r: dict) -> float:
        return abs(float(r["protein_aa"]) - float(r["canonical_aa"])) / float(r["canonical_aa"])

    unflagged_large = sorted(
        (r for r in audit_rows if r["flag_class"] != "proposed" and _gap(r) > 0.10),
        key=lambda r: r["gene_symbol"],
    )
    queried = [r for r in audit_rows if r["atlas_peptides_member_only"] != ""]
    canonical_only = [
        r
        for r in queried
        if r["atlas_peptides_member_only"] == "0"
        and r["atlas_peptides_canonical_only"].isdigit()
        and int(r["atlas_peptides_canonical_only"]) > 0
    ]
    n["proteoform_audit"] = {
        "flagged_major": sum(r["flag_class"] == "proposed" for r in audit_rows),
        "unflagged_over_10pct": len(unflagged_large),
        "unflagged_over_10pct_genes": [
            f"{r['gene_symbol']} {int(float(r['protein_aa']))}/{int(float(r['canonical_aa']))}"
            for r in unflagged_large
        ],
        "atlas_queried": len(queried),
        "atlas_canonical_only": len(canonical_only),
    }
    # --- the 2026-11-22 proteoform review (proteoform_review.py): PeptideAtlas build 607 queried
    # for every reviewed accession of each gene whose reference differed from the canonical.
    # The brief's section 4 states these; the audit's pinned-pool counts above are narrower.
    review = json.loads((S4E9 / "proteoform_review_2026-11-22_counts.json").read_text())
    n["proteoform_review"] = {
        "reviewed": review["reviewed"],
        "accessions_queried": review["accessions_queried"],
        "evidence_class": review["evidence_class"],
        "canonical_strong_ours_none": review["canonical_strong_ours_none"],
    }

    # --- the pooled barcode-swap model (l1_swap_pool_model.py), as percentages of integrants
    pool = json.loads((S4E9 / "l1_swap_pool_model_summary.json").read_text())
    frac = pool["library_fraction_wrong_protein"]
    pct = lambda x: round(100 * x, 1)  # noqa: E731
    n["swap_pool"] = {
        "pool_share_1_over": pool["pool_share_1_over"],
        "pool_weighting": pool["pool_weighting"],
        "sibling_tract_median_nt": int(pool["same_cds_pairs"]["suffix_tract_median_nt"]),
        "sibling_rate_median_pct": pct(pool["same_cds_pairs"]["per_pair_rate_median"]),
        "tail_mislabel_pct": pct(frac["background_mislabel"]),
        "shared_tag_chimera_pct": pct(frac["background_chimera"]),
        "specific_pct": pct(frac["specific_mislabel"] + frac["specific_chimera"]),
        "specific_chimera_pct": round(100 * frac["specific_chimera"], 2),
        "specific_mislabel_pct": round(100 * frac["specific_mislabel"], 2),
        "specific_median_per_barcode_pct": round(
            100 * pool["per_barcode_impurity"]["specific_part_median"], 2
        ),
        "total_pct": pct(frac["total"]),
    }
    # Specific-partner impurity against insert length: long families carry more of it.
    nt = {r["insert_id"]: int(r["insert_nt"]) for r in _tsv(HERE / "l1_order_plan.tsv")}
    per = sorted(
        (nt[r["insert_id"]], float(r["specific_mislabel"]) + float(r["specific_chimera"]))
        for r in _tsv(S4E9 / "l1_swap_pool_model.tsv")
    )
    lengths = sorted(x for x, _ in per)
    q1, q3 = lengths[len(lengths) // 4], lengths[3 * len(lengths) // 4]

    def _median(v: list[float]) -> float:
        v = sorted(v)
        return (v[(len(v) - 1) // 2] + v[len(v) // 2]) / 2

    def _ranks(v: list[float]) -> list[float]:
        order = sorted(range(len(v)), key=v.__getitem__)
        r = [0.0] * len(v)
        i = 0
        while i < len(v):
            j = i
            while j + 1 < len(v) and v[order[j + 1]] == v[order[i]]:
                j += 1
            for k in range(i, j + 1):
                r[order[k]] = (i + j) / 2
            i = j + 1
        return r

    def _spearman(a: list[float], b: list[float]) -> float:
        ra, rb = _ranks(a), _ranks(b)
        ma, mb = sum(ra) / len(ra), sum(rb) / len(rb)
        cov = sum((x - ma) * (y - mb) for x, y in zip(ra, rb, strict=True))
        var = (sum((x - ma) ** 2 for x in ra) * sum((y - mb) ** 2 for y in rb)) ** 0.5
        return cov / var

    n["swap_pool"]["by_length"] = {
        "longest_quarter_specific_median_pct": round(
            100 * _median([s for x, s in per if x >= q3]), 2
        ),
        "shortest_quarter_specific_median_pct": round(
            100 * _median([s for x, s in per if x <= q1]), 2
        ),
        "spearman_length_vs_specific": round(_spearman(*map(list, zip(*per, strict=True))), 2),
        "barcodes_over_1pct": sum(s > 0.01 for _, s in per),
        "max_specific_pct": pct(max(s for _, s in per)),
    }

    # --- literature tag-site positives behind the S5E13 ranker, read against their papers
    from deliverome_analysis.s3 import processed

    rows = _tsv(processed.download(AUDIT_KEY))
    corrected = _tsv(processed.download(CORRECTED_KEY))
    # The ranker's positives are the training rows; the audit also covers live API rows
    # outside training, reported separately so they never change the training counts.
    audit = [r for r in rows if r["in_training"] == "1"]
    outside = [r for r in rows if r["in_training"] != "1"]
    shifted = [r for r in audit if r["verdict"] == "confirmed_shifted"]
    n["tag_site_audit"] = {
        "positives": len(audit),
        "verdict": dict(Counter(r["verdict"] for r in audit).most_common()),
        "outside_training_rows": len(outside),
        "outside_training_verdict": dict(Counter(r["verdict"] for r in outside).most_common()),
        "outside_training_pinned_to_a_residue": sorted(
            f"{r['gene']} {int(float(r['corrected_residue']))}"
            for r in outside
            if not r["insert_after_residue"] and r["corrected_residue"]
        ),
        "shifted_by_numbering_basis": dict(
            Counter(r["paper_numbering_basis"] for r in shifted).most_common()
        ),
        "abcb1_shifts": sorted(
            int(float(r["corrected_residue"])) - int(r["insert_after_residue"])
            for r in shifted
            if r["gene"] == "ABCB1"
        ),
        "corrected_corpus_rows": len(corrected),
        "corrected_corpus_residue_positives": sum(
            r["verdict"] in ("confirmed_exact", "confirmed_shifted") for r in corrected
        ),
        "corrected_corpus_region_labels": sum(r["verdict"] == "region_only" for r in corrected),
    }
    gold = _tsv(processed.download(GOLD_KEY))
    shown = [r for r in gold if r["surface_shown"] == "yes"]
    n["tag_site_gold_standard"] = {
        "sites": len(gold),
        "proteins": len({r["uniprot_acc"] for r in gold}),
        "surface_shown_sites": len(shown),
        "surface_shown_proteins": len({r["uniprot_acc"] for r in shown}),
        "surface_negative_sites": sum(r["surface_shown"] == "no" for r in gold),
        "function_shown_sites": sum(r["function_shown"] == "yes" for r in gold),
        "training_proteins_covered": len(
            {r["uniprot_acc"] for r in audit} & {r["uniprot_acc"] for r in gold}
        ),
        "training_proteins": len({r["uniprot_acc"] for r in audit}),
        "sites_only_from_added": sorted(
            f"{r['gene']} {r['insert_after_residue']}"
            for r in gold
            if all("added" in sid for sid in r["site_ids"].split(";"))
        ),
    }
    failed = _tsv(processed.download(AUDIT_DIR + "failed_sites.tsv"))
    regions = _tsv(processed.download(AUDIT_DIR + "region_labels.tsv"))
    windowed = [r for r in regions if r["window_start"]]
    added = _tsv(processed.download(AUDIT_DIR + "added_sites.tsv"))
    n["tag_site_gold_standard"]["added_rows"] = len(added)
    n["tag_site_gold_standard"]["added_rows_confirmed"] = sum(
        r["verdict"].startswith("confirmed") for r in added
    )
    n["tag_site_gold_standard"]["failed_sites_excluded"] = len(failed)
    n["tag_site_gold_standard"]["region_labels"] = len(regions)
    n["tag_site_gold_standard"]["region_labels_with_window"] = len(windowed)
    n["tag_site_gold_standard"]["region_label_window_proteins"] = len(
        {r["uniprot_acc"] for r in windowed}
    )

    # --- cassette encoding and swap risk -------------------------------------------------
    swap = json.loads((S4E9 / "l1_swap_risk_summary.json").read_text())
    enc = _tsv(S4E9 / "l1_cassette_encodings.tsv")
    n["swap_risk"] = {
        "members": swap["members"],
        "gate_bp": swap["gate_bp"],
        "schemes": {
            k: {
                f: swap[k][f]
                for f in (
                    "worst_partner_tract_median",
                    "members_over_gate_100",
                    "cterm_members",
                    "cterm_worst_partner_median",
                    "cterm_worst_partner_max",
                )
            }
            for k in ("shared", "whole_cassette_widened", "linkers_only")
        },
        "whole_cassette_9nt_rule": swap["whole_cassette_draw"].split(":")[0],
        "issued_linker_only_encodings": sum(r["scheme"] == "linkers_only" for r in enc),
        "issue_batches": sorted({r["issued_batch"] for r in enc}),
        "cterm_junction_edits": swap["cterm_junction_edits"],
    }

    # --- the isoform panel as drawn and its class composition (isoform_panel_summary.json from
    # the sampler, isoform_architecture_classes.tsv from isoform_architecture_classes.py), so
    # the redraw in the rebuild refills section 5 by rerunning this script.
    panel = json.loads((S4E9 / "isoform_panel_summary.json").read_text())
    classes = {r["class"]: r for r in _tsv(S4E9 / "isoform_architecture_classes.tsv")}
    n["isoform_panel"] = {
        "size": panel["panel_size"],
        "genes": panel["genes"],
        "replaced": panel["replaced"],
        "refilled_after_draw": panel.get("refilled_after_draw", 0),
        "by_class": panel["by_class"],
        "classes": {
            c: {
                k: int(r[k]) if r[k].isdigit() else float(r[k])
                for k in (
                    "library_genes",
                    "cohort_forms",
                    "drawable_forms_gate_of_record",
                    "panel_forms",
                    "panel_forms_observed",
                    "panel_forms_prediction_only",
                    "uniform_expected_from_cohort",
                    "enrichment_vs_uniform_cohort",
                    "uniform_expected_from_drawable",
                    "enrichment_vs_uniform_drawable",
                )
            }
            for c, r in classes.items()
            if int(r["library_genes"]) > 0 and c != "unclassified"
        },
    }
    return n


# ------------------------------------------------------------------------------- figures


def figures(n: dict) -> None:
    import numpy as np

    from deliverome_analysis.plotting import (
        create_figure,
        get_color,
        save_figure,
        setup_plotting_style,
    )

    setup_plotting_style()

    # Figure 1 (Figure 6 in BRIEF.md): the reference ladder, one bar per decision, split by source.
    order = [
        "mane_appris_agree",
        "appris_only",
        "no_reference_uniprot_canonical",
        "both_observed",
        "resolved_mane",
        "weak_evidence_mane",
        "weak_evidence_appris",
        "unresolved_default_mane",
    ]
    labels = [
        "MANE = APPRIS",
        "APPRIS only",
        "neither (UniProt)",
        "disagree: both strongly observed",
        "disagree: strong MANE evidence",
        "disagree: weak MANE evidence",
        "disagree: weak APPRIS evidence",
        "disagree: no evidence",
    ]
    sources = ["mane", "appris", "uniprot_canonical", "named"]
    colors = {
        "mane": get_color("secondary"),
        "appris": get_color("accent"),
        "uniprot_canonical": get_color("quaternary"),
        "named": get_color("primary"),
    }
    by = n["reference"]["decision_by_source"]
    fig, ax = create_figure(figsize=(12, 7))
    y = np.arange(len(order))[::-1]
    left = np.zeros(len(order))
    for s in sources:
        vals = np.array([by.get(d, {}).get(s, 0) for d in order], dtype=float)
        ax.barh(y, vals, left=left, color=colors[s], label=f"taken: {s.replace('_', ' ')}")
        left += vals
    for yi, tot in zip(y, left, strict=True):
        ax.text(tot * 1.15, yi, f"{int(tot):,}", va="center")
    ax.set_xscale("log")
    ax.set_xlim(1, 2.0e4)
    ax.set_yticks(y, labels)
    ax.set_xlabel("library genes (log scale)")
    ax.legend(loc="lower right", frameon=False)
    save_figure(fig, "l1_reference_ladder", output_dir=FIGURES)

    # Figure 2: placement rung per member, by candidate set.
    set_colors = {
        "library_reference": get_color("secondary"),
        "isoform": get_color("accent"),
        "ortholog": get_color("quaternary"),
    }
    pbs = n["tag_paths"]["placement_by_set"]
    fig, ax = create_figure(figsize=(12, 6))
    y = np.arange(len(PLACEMENTS))[::-1]
    left = np.zeros(len(PLACEMENTS))
    for s in SETS:
        vals = np.array([pbs[p][s] for p in PLACEMENTS], dtype=float)
        ax.barh(y, vals, left=left, color=set_colors[s], label=s.replace("_", " "))
        left += vals
    for yi, tot in zip(y, left, strict=True):
        ax.text(tot + 15, yi, f"{int(tot):,}", va="center")
    ax.set_yticks(y, PLACEMENTS)
    ax.set_xlabel("L1 members (DeepTMHMM v1 topology of the ordered protein)")
    ax.set_xlim(0, left.max() * 1.15)
    ax.legend(loc="lower right", frameon=False)
    save_figure(fig, "l1_tag_rungs", output_dir=FIGURES)

    # Figure 3: worst-partner shared tract per member, shared cassette vs linkers-only.
    rows = _tsv(S4E9 / "l1_swap_risk.tsv")
    fig, ax = create_figure(figsize=(11, 6.5))
    for scheme, color, label in (
        ("shared", get_color("primary"), "one shared cassette (p003 codons)"),
        ("linkers_only", get_color("secondary"), "linkers recoded per member (issued)"),
    ):
        v = np.sort(np.array([int(r[f"tract_nt_{scheme}"]) for r in rows]))
        frac_over = 1.0 - np.arange(1, len(v) + 1) / len(v)
        ax.step(v, frac_over, where="post", color=color, label=label, lw=2.5)
    ax.axvline(n["swap_risk"]["gate_bp"], color=get_color("neutral"), ls="--", lw=1.5)
    ax.text(n["swap_risk"]["gate_bp"] * 1.08, 2e-3, "100-nt gate", color=get_color("neutral"))
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim(1 / len(rows), 1.0)
    ax.set_xlabel("longest tract shared with any other member, up to the barcode (nt)")
    ax.set_ylabel("fraction of members above")
    ax.legend(loc="upper right", frameon=False)
    save_figure(fig, "l1_swap_tracts", output_dir=FIGURES)

    _draw_candidate_sets_figure(FIGURES)
    _draw_promoter_comparison_figure(FIGURES)


def _draw_candidate_sets_figure(out_dir: Path) -> None:
    """Two-panel composition figure: counts and CDS length distributions by set and species."""
    import matplotlib.pyplot as plt
    import numpy as np

    from deliverome_analysis.plotting import get_color, save_figure

    rows = _tsv(S4E9 / "l1_tag_paths.tsv")
    groups = [
        ("library_reference", "human", "Human reference", get_color("secondary")),
        ("isoform", "human", "Human isoforms", get_color("accent")),
        ("ortholog", "mouse", "Mouse orthologs", get_color("quaternary")),
        ("ortholog", "cyno", "Cynomolgus orthologs", get_color("primary")),
    ]
    lengths_kb = [
        np.array(
            [
                int(r["length_aa"]) * 3 / 1000.0
                for r in rows
                if r["candidate_set"] == cs and r["species"] == sp
            ],
            dtype=float,
        )
        for cs, sp, _, _ in groups
    ]
    labels = [g[2] for g in groups]
    colors = [g[3] for g in groups]
    counts = [len(a) for a in lengths_kb]
    y = np.arange(len(groups))[::-1]

    fig, (ax_cnt, ax_len) = plt.subplots(
        1, 2, figsize=(13.5, 4.2), gridspec_kw={"width_ratios": [1.05, 1.35]}, sharey=True
    )
    ax_cnt.barh(y, counts, height=0.58, color=colors, edgecolor="#1F2933", linewidth=0.6)
    for yi, cnt in zip(y, counts, strict=True):
        ax_cnt.text(
            cnt + 45, yi, f"{cnt:,}", va="center", ha="left", fontsize=11, fontweight="semibold"
        )
    ax_cnt.set_yticks(y, labels, fontsize=11)
    ax_cnt.set_xlim(0, max(counts) * 1.2)
    ax_cnt.set_xlabel(f"Active CDSs (total {sum(counts):,})", fontsize=11)
    ax_cnt.text(
        0.0,
        1.04,
        "A  Library composition by candidate set",
        transform=ax_cnt.transAxes,
        ha="left",
        va="bottom",
        fontsize=12,
        fontweight="semibold",
    )
    ax_cnt.spines[["top", "right"]].set_visible(False)

    parts = ax_len.violinplot(lengths_kb, positions=y, vert=False, widths=0.65, showextrema=False)
    for body, col in zip(parts["bodies"], colors, strict=True):
        body.set_facecolor(col)
        body.set_edgecolor(col)
        body.set_alpha(0.32)
    for yi, arr, col in zip(y, lengths_kb, colors, strict=True):
        q25, med, q75 = np.percentile(arr, [25, 50, 75])
        ax_len.plot([q25, q75], [yi, yi], color=col, lw=4.5, solid_capstyle="round", zorder=3)
        ax_len.scatter([med], [yi], color="white", edgecolor=col, s=52, lw=1.8, zorder=4)
        ax_len.text(
            q75 + 0.14,
            yi + 0.24,
            f"median {med:.2f} kb",
            va="center",
            ha="left",
            fontsize=9.5,
            color="#334155",
        )
    ax_len.axvline(5.0, color="#64748B", ls="--", lw=1.2)
    ax_len.text(
        4.92, y[0] + 0.32, "5.0 kb cap", ha="right", va="center", fontsize=9.5, color="#64748B"
    )
    ax_len.set_xlim(0, 5.55)
    ax_len.set_xlabel("Native CDS length (kb)", fontsize=11)
    ax_len.text(
        0.0,
        1.04,
        "B  CDS length distribution (median and IQR)",
        transform=ax_len.transAxes,
        ha="left",
        va="bottom",
        fontsize=12,
        fontweight="semibold",
    )
    ax_len.spines[["top", "right"]].set_visible(False)

    fig.tight_layout(rect=(0, 0, 1, 0.94), w_pad=2.2)
    save_figure(fig, "l1_candidate_sets", output_dir=out_dir)


_PROMOTER_COLORS = {
    "ef1a": "#376B62",
    "efs": "#9E2A3B",
    "orig": "#ED9B96",
    "sv40": "#8B95A1",
    "pgk": "#E6A328",
}


def _draw_promoter_panel_a(ax) -> None:
    """Panel A: 72 h pre-selection mScarlet3 gMFI fold-change vs PGK (S5E11 gmfi.tsv)."""
    import numpy as np

    cells_a = ["HEK293T", "HeLa", "Nalm-6", "K562"]
    fold_a = {
        "EF1a": {
            "HEK293T": [10.83, 10.61, 11.59, 12.48],
            "HeLa": [5.21, 5.13, 5.52, 5.50],
            "Nalm-6": [3.43, 3.44, 3.13, 3.02],
            "K562": [2.51, 2.43, 2.25],
        },
        "EFS": {
            "HEK293T": [2.68, 2.91, 3.08, 3.51],
            "HeLa": [None, 2.12, 3.27, 3.07],
            "Nalm-6": [1.56, 1.63, 1.76, 1.79],
            "K562": [],
        },
        "EFS_orig": {
            "HEK293T": [2.02, 1.96, 1.81, 1.62],
            "HeLa": [None, 0.68, 0.44, 0.48],
            "Nalm-6": [1.11, 1.14, 1.08, 1.08],
            "K562": [1.09, 1.00, 0.90],
        },
    }
    markers = ["o", "s", "^", "D"]
    xa = np.arange(len(cells_a))
    w = 0.25
    series_a = [
        ("EF1a", -w, _PROMOTER_COLORS["ef1a"], 0.8),
        ("EFS", 0.0, _PROMOTER_COLORS["efs"], 1.8),
        ("EFS_orig", w, _PROMOTER_COLORS["orig"], 0.8),
    ]
    ax.axhline(1.0, color=_PROMOTER_COLORS["pgk"], ls="--", lw=1.8, zorder=1)
    for key, dx, col, lw in series_a:
        for i, cell in enumerate(cells_a):
            real = [(j, v) for j, v in enumerate(fold_a[key][cell]) if v is not None]
            if not real:
                ax.text(
                    xa[i] + dx,
                    1.4,
                    "n/a",
                    ha="center",
                    va="bottom",
                    rotation=90,
                    fontsize=9,
                    color="#64748B",
                )
                continue
            med = float(np.median([v for _, v in real]))
            ax.bar(
                xa[i] + dx,
                med - 1.0,
                bottom=1.0,
                width=w * 0.92,
                color=col,
                edgecolor="#1F2933",
                linewidth=lw,
                zorder=2,
            )
            offsets = np.linspace(-w * 0.26, w * 0.26, len(real))
            for (j, v), off in zip(real, offsets, strict=True):
                ax.scatter(
                    xa[i] + dx + off,
                    v,
                    marker=markers[j],
                    s=26,
                    facecolor="white",
                    edgecolor="#1F2933",
                    linewidth=0.9,
                    zorder=4,
                )
    ax.set_yscale("log")
    ax.set_ylim(0.35, 18)
    ax.set_yticks([0.5, 1, 2, 5, 10], ["0.5", "1", "2", "5", "10"])
    ax.set_xticks(xa, cells_a, fontsize=10.5)
    ax.set_ylabel("72 h mScarlet3 gMFI, fold vs. PGK (log)", fontsize=10.5)
    ax.text(
        0.0,
        1.04,
        "A  Pre-selection expression (72 h, 4 virus doses)",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=11,
        fontweight="semibold",
    )
    ax.spines[["top", "right"]].set_visible(False)


def _draw_promoter_panel_b(ax) -> None:
    """Panel B: Post-selection % SpyCatcher+ (37 C) + EF1a/EFS gMFI gain callout."""
    import numpy as np

    cells_bc = ["HEK293T", "K562", "Nalm-6"]
    spy_pct = {"EF1a": [95, 92, 94], "EFS": [83, 76, 87], "SV40": [60, 45, 27]}
    spy_gain = ["2.0×", "2.4×", "1.5×"]  # EF1a / EFS SpyCatcher gMFI at 37 C
    xb = np.arange(len(cells_bc))
    wb = 0.25
    for key, dx, col, lw in [
        ("EF1a", -wb, _PROMOTER_COLORS["ef1a"], 0.8),
        ("EFS", 0.0, _PROMOTER_COLORS["efs"], 1.8),
        ("SV40", wb, _PROMOTER_COLORS["sv40"], 0.8),
    ]:
        vals = spy_pct[key]
        ax.bar(
            xb + dx, vals, width=wb * 0.92, color=col, edgecolor="#1F2933", linewidth=lw, zorder=2
        )
        for i, v in enumerate(vals):
            ax.text(xb[i] + dx, v + 1.8, f"{v}%", ha="center", va="bottom", fontsize=9)
    for i, g in enumerate(spy_gain):
        ax.text(
            xb[i] - wb * 0.5,
            106,
            f"EF1α {g} gMFI",
            ha="center",
            va="bottom",
            fontsize=8.8,
            color=_PROMOTER_COLORS["ef1a"],
            fontweight="semibold",
        )
    ax.set_ylim(0, 118)
    ax.set_xticks(xb, cells_bc, fontsize=10.5)
    ax.set_ylabel("SpyCatcher+ live cells at 37 °C (%)", fontsize=10.5)
    ax.text(
        0.0,
        1.04,
        "B  Post-selection surface display (TfR1)",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=11,
        fontweight="semibold",
    )
    ax.spines[["top", "right"]].set_visible(False)


def _draw_promoter_panel_c(ax) -> None:
    """Panel C: Post-selection % internalized (37 C, repeat run; first run in diamonds)."""
    import numpy as np

    cells_bc = ["HEK293T", "K562", "Nalm-6"]
    pct_int = {"EF1a": [34, 45, 23], "EFS": [41, 53, 36], "SV40": [41, 41, 41]}
    pct_int_r1 = {"EF1a": [None, 44, 23], "EFS": [None, 57, 26], "SV40": [None, 42, 35]}
    xb = np.arange(len(cells_bc))
    wb = 0.25
    for key, dx, col, lw in [
        ("EF1a", -wb, _PROMOTER_COLORS["ef1a"], 0.8),
        ("EFS", 0.0, _PROMOTER_COLORS["efs"], 1.8),
        ("SV40", wb, _PROMOTER_COLORS["sv40"], 0.8),
    ]:
        vals = pct_int[key]
        ax.bar(
            xb + dx, vals, width=wb * 0.92, color=col, edgecolor="#1F2933", linewidth=lw, zorder=2
        )
        for i, v in enumerate(vals):
            r1 = pct_int_r1[key][i]
            top = max(v, r1 or 0)
            ax.text(xb[i] + dx, top + 2.2, f"{v}%", ha="center", va="bottom", fontsize=9)
            if r1 is not None:
                ax.scatter(
                    xb[i] + dx,
                    r1,
                    marker="D",
                    s=34,
                    facecolor="white",
                    edgecolor="#1F2933",
                    linewidth=1.1,
                    zorder=4,
                )
    ax.set_ylim(0, 72)
    ax.set_xticks(xb, cells_bc, fontsize=10.5)
    ax.set_ylabel("TfR1 internalized at 37 °C (%)", fontsize=10.5)
    ax.text(
        0.0,
        1.04,
        "C  Post-selection internalization (30 min)",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=11,
        fontweight="semibold",
    )
    ax.spines[["top", "right"]].set_visible(False)


def _promoter_legend_handles() -> list:
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    return [
        Patch(
            facecolor=_PROMOTER_COLORS["ef1a"],
            edgecolor="#1F2933",
            linewidth=0.8,
            label="EF1α + domesticated IRES2",
        ),
        Patch(
            facecolor=_PROMOTER_COLORS["efs"],
            edgecolor="#1F2933",
            linewidth=1.8,
            label="EFS + domesticated IRES2 (current baseline)",
        ),
        Patch(
            facecolor=_PROMOTER_COLORS["orig"],
            edgecolor="#1F2933",
            linewidth=0.8,
            label="EFS + undomesticated IRES2 (Panel A)",
        ),
        Patch(
            facecolor=_PROMOTER_COLORS["sv40"],
            edgecolor="#1F2933",
            linewidth=0.8,
            label="EFS + separate SV40–PuroR (Panels B–C)",
        ),
        Line2D(
            [0],
            [0],
            color=_PROMOTER_COLORS["pgk"],
            ls="--",
            lw=1.8,
            label="PGK baseline = 1 (Panel A)",
        ),
        Line2D(
            [0],
            [0],
            marker="D",
            color="none",
            markerfacecolor="white",
            markeredgecolor="#1F2933",
            markersize=6,
            label="First assay run (Panel C)",
        ),
    ]


def _draw_promoter_comparison_figure(out_dir: Path) -> None:
    """Three-panel promoter and selection comparison in the SAB-deck style."""
    import matplotlib.pyplot as plt

    from deliverome_analysis.plotting import save_figure

    fig, (ax_a, ax_b, ax_c) = plt.subplots(
        1, 3, figsize=(16.5, 5.2), gridspec_kw={"width_ratios": [1.25, 1.0, 1.0]}
    )
    _draw_promoter_panel_a(ax_a)
    _draw_promoter_panel_b(ax_b)
    _draw_promoter_panel_c(ax_c)
    fig.legend(
        handles=_promoter_legend_handles(),
        loc="lower center",
        bbox_to_anchor=(0.5, -0.02),
        ncol=3,
        frameon=False,
        fontsize=9.8,
    )
    fig.tight_layout(rect=(0, 0.10, 1, 0.94), w_pad=2.0)
    save_figure(fig, "l1_promoter_comparison", output_dir=out_dir)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--no-figures", action="store_true")
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args(argv)
    n = measure()
    import review_fingerprint as rf

    n["review_fingerprint"] = rf.compute(S4E9)
    args.out.write_text(json.dumps(n, indent=1) + "\n")
    print(f"wrote {args.out}")
    if not args.no_figures:
        figures(n)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
