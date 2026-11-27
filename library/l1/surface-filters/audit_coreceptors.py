"""Reproduce the published baseline and audit every required/stable-complex target.

Run after hydrating the pinned inputs described in CORECEPTOR_AUDIT.md:
uv run python library/l1/surface-filters/audit_coreceptors.py
"""

import argparse
import contextlib
import copy
import csv
import hashlib
import io
import json
import runpy
import shutil
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
BASELINE_REF = "6b3d53ccab9b95246a0605b79427d69dca36aeb2"
SNAPSHOT_PREFIX = "sprint1/s1e3-surface-library-filters/partner-audit-2026-09-10"
INPUTS = (
    "cohort.jsonl",
    "catalog_rows.json",
    "evidence_grade_summary.json",
    "coreceptor_partners.jsonl",
    "coreceptor_calls.tsv",
    "coreceptor_rationales.json",
    "oe_calls_136.tsv",
    "willow_screen_rescue_genes.json",
    "willow_screen_rescue_isoforms.json",
    "external/expression_tpm/resolved_expr.json",
    "external/expression_tpm/resolved_expr_full.json",
    "external/hpa_rna_celline/hpa_celline_per_gene.csv",
)


def run_pipeline(folder):
    with contextlib.redirect_stdout(io.StringIO()):
        return runpy.run_path(str(folder / "final_library.py"))


def failed_gates(state, record):
    filters = record["filters"]
    gene = record["gene"]
    membership = state["B"].passes_canonical(filters) or state["B"].passes_likely(filters)
    localization = (
        state["localization_ok"](record)
        if "localization_ok" in state
        else filters.get("surface_call_reason") not in state["B"].NONSURFACE_REASONS
    )
    gates = {
        "membership": membership,
        "accessibility": filters.get("surface_accessibility") != "no",
        "specificity": filters.get("surface_specificity") in ("surface_dominant", "mixed"),
        "anchor": state["anchored"](record),
        "partner": state["coreceptor_ok"](gene, filters),
        "expression": state["expr_ok"](gene, filters),
        "localization": localization,
        "orf": state["orf_ok"](record),
        "ecd": state["ecd_ok"](record),
    }
    return [name for name, passed in gates.items() if not passed]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hydrate", action="store_true", help="Download the frozen audit inputs")
    args = parser.parse_args()
    if args.hydrate:
        from deliverome_analysis.core.s3 import processed

        manifest = json.loads(
            processed.download(f"{SNAPSHOT_PREFIX}/coreceptor_audit_summary.json").read_text()
        )
        for name, expected_hash in manifest["input_sha256"].items():
            destination = HERE / name
            processed.download(f"{SNAPSHOT_PREFIX}/inputs/{name}", dest=destination)
            assert hashlib.sha256(destination.read_bytes()).hexdigest() == expected_hash, name
    sys.path.insert(0, str(HERE))
    from coreceptor_rules import REVIEWED_RULES

    from deliverome_analysis.core.provenance import read_pinned_text

    baseline_source = read_pinned_text(
        BASELINE_REF, (HERE / "final_library.py").relative_to(ROOT).as_posix(), repo=ROOT
    )
    with tempfile.TemporaryDirectory(prefix="s1e3-baseline-") as temporary:
        baseline_folder = Path(temporary)
        for name in INPUTS:
            if (HERE / name).exists():
                (baseline_folder / name).parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(HERE / name, baseline_folder / name)
        (baseline_folder / "final_library.py").write_text(baseline_source)
        baseline = run_pipeline(baseline_folder)
    current = run_pipeline(HERE)
    published_path = HERE / "published_library_membership.tsv"
    with published_path.open() as handle:
        published = {row["gene"]: row["tier"] for row in csv.DictReader(handle, delimiter="\t")}
    assert baseline["_lib_tier"] == published, "Baseline differs from the published membership"
    added = sorted(set(current["_lib_tier"]) - set(published))
    removed = sorted(set(published) - set(current["_lib_tier"]))
    assert added == ["ATP1B3", "CDC50A", "HFE", "P2RX6"], added
    assert not removed, removed
    assert all(current["_lib_tier"][gene] == tier for gene, tier in published.items())
    assert len(published) == 2340 and len(current["_lib_tier"]) == 2344
    negative_checks = 0
    for gene in added:
        original = current["_by_gene"][gene]
        for field, value in (
            ("surface_accessibility", "no"),
            ("surface_specificity", "mostly_intracellular"),
            ("deep_dive_tier", "excluded"),
            ("surface_call_reason", "secreted_only"),
        ):
            modified = copy.deepcopy(original)
            modified["filters"][field] = value
            assert not current["passes"](modified), (gene, field)
            negative_checks += 1
        for topology in (
            {"tm_helix_count": 0, "per_residue_topology": "O" * 300},
            {"tm_helix_count": 1, "per_residue_topology": "I" * 300},
            {"tm_helix_count": 1, "per_residue_topology": "O" * 2000},
        ):
            modified = copy.deepcopy(original)
            modified["canonical_topology"] = topology
            assert not current["passes"](modified), gene
            negative_checks += 1
    rows = []
    for record in current["recs"]:
        gene, filters = record["gene"], record["filters"]
        if (
            filters.get("co_receptor_dependency") != "required"
            and filters.get("surface_call_reason") != "stable_complex_partner"
        ):
            continue
        partners = current["part"].get(gene, {}).get("partners", [])
        rows.append(
            {
                "gene": gene,
                "uniprot": record.get("uniprot"),
                "reviewed_rule": record.get("uniprot") in REVIEWED_RULES,
                "partners_and_median_tpm": json.dumps(
                    {partner: current["partner_median"](partner) for partner in partners}
                ),
                "old_failed_gates": ";".join(failed_gates(baseline, record)),
                "new_failed_gates": ";".join(failed_gates(current, record)),
                "old_tier": published.get(gene, "excluded"),
                "new_tier": current["_lib_tier"].get(gene, "excluded"),
                "surface_isoforms": current["n_surface_isoforms"](record),
                "surface_orthologs": current["n_surface_orthologs"](record),
            }
        )
    with (HERE / "coreceptor_audit.tsv").open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(sorted(rows, key=lambda row: row["gene"]))
    summary = {
        "baseline_ref": BASELINE_REF,
        "cohort_genes": len(current["recs"]),
        "required_genes": sum(
            record["filters"].get("co_receptor_dependency") == "required"
            for record in current["recs"]
        ),
        "audited_genes": len(rows),
        "added": [row for row in rows if row["gene"] in added],
        "removed": removed,
        "retiered": [],
        "baseline_counts": baseline["_counts"],
        "current_counts": current["_counts"],
        "negative_integration_checks": negative_checks,
        "input_sha256": {
            name: hashlib.sha256((HERE / name).read_bytes()).hexdigest()
            for name in (*INPUTS, "published_library_membership.tsv")
            if (HERE / name).exists()
        },
    }
    (HERE / "coreceptor_audit_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(
        json.dumps(
            {
                key: value
                for key, value in summary.items()
                if key not in {"input_sha256", "baseline_counts", "current_counts"}
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
