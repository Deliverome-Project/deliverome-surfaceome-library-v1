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

from __future__ import annotations

import argparse
import csv
import json
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

    return _tsv(processed.download(PREFIX + "library-design/reference_selection.tsv"))


def measure() -> dict:
    n: dict = {}

    # --- gene list -----------------------------------------------------------------
    root = json.loads(
        (S4E9.parents[2] / "scratch/becca/tag-site-architecture/library_members.json").read_text()
    )
    added = [r["gene"] for r in _tsv(S4E9 / "gene_list_additions.tsv")]
    excluded = [r["gene"] for r in _tsv(S4E9 / "gene_list_exclusions.tsv")]
    n["gene_list"] = {
        "frozen_root_genes": len(root),
        "added": added,
        "excluded_from_filter_output": excluded,
        "library_genes": len(root) + len(added),
    }

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
        "igk_leader_members": sum(r["has_igk_leader"] == "True" for r in active),
    }
    for name in ("isoform", "ortholog"):
        s = json.loads((S4E9 / f"{name}_panel_summary.json").read_text())
        n["tranches"][f"{name}_panel"] = s

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

    # Figure 1: the reference ladder, one bar per decision, split by the source taken.
    order = [
        "mane_appris_agree",
        "appris_only",
        "no_reference_uniprot_canonical",
        "resolved_mane",
        "weak_evidence_mane",
        "weak_evidence_appris",
        "unresolved_default_mane",
    ]
    labels = [
        "MANE = APPRIS",
        "APPRIS only",
        "neither (UniProt)",
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


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--no-figures", action="store_true")
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args(argv)
    n = measure()
    args.out.write_text(json.dumps(n, indent=1) + "\n")
    print(f"wrote {args.out}")
    if not args.no_figures:
        figures(n)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
