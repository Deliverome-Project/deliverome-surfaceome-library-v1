"""Explain genes absent from the selected isoform cohort, using its frozen replay.

This audits historical selection; it never assigns new topology, changes evidence
thresholds, or redraws a panel. The historical gene/length/accession verdict join
is replayed solely to explain the original decision, not to transfer a prediction.

Run reproduce_isoform_universe.py --step isoforms first, then pass its WORK/iso
directory as --replay. The replay must reproduce the recorded cohort exactly.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

import library_pipeline as lp
import reproduce_isoform_universe as riu

REASONS = (
    "reference",
    "insufficient_evidence",
    "loses_all_tm_vs_reference",
    "exclude_lost_surface",
    "soluble_restriction",
)


def table(path: Path) -> list[dict[str, str]]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def decision(row: dict, verdict: dict, restricted: dict) -> str:
    """Original cohort filter precedence (6fb98222); no new eligibility rules."""
    if row["is_reference"] == "True":
        return "reference"
    observed = row["tier"] in {"expressed", "likely_expressed", "weak"} or (
        int(row["openprot_te_score"] or 0) > 0
    )
    if not (observed or float(row["trifid_score"] or 0) > 0.5):
        return "insufficient_evidence"
    if row["loses_all_tm_vs_reference"] == "True":
        return "loses_all_tm_vs_reference"
    if verdict.get("verdict") == "exclude_lost_surface":
        return "exclude_lost_surface"
    if restricted:
        return "soluble_restriction"
    raise ValueError(f"Unexplained cohort absence: {row}")


def explain(
    missing: set[str],
    scored: list[dict],
    verdicts: list[dict],
    restricted: list[dict],
    original_genes: set[str],
) -> tuple[list, list, dict]:
    verdict_by = {(r["gene"], r["length"], r["uniprot"]): r for r in verdicts}
    soluble_by = {(r["gene"], r["form"]): r for r in restricted}
    forms, by_gene = [], defaultdict(Counter)
    for row in scored:
        if row["gene"] not in missing:
            continue
        verdict = verdict_by.get((row["gene"], row["length"], row["uniprot"]), {})
        soluble = soluble_by.get((row["gene"], row["uniprot"] or row["ensembl_transcripts"]), {})
        reason = decision(row, verdict, soluble)
        by_gene[row["gene"]][reason] += 1
        forms.append(
            {
                **row,
                "exclusion_reason": reason,
                "historical_topology_verdict": verdict.get("verdict", ""),
                "historical_soluble_reason": soluble.get("reason", ""),
            }
        )
    if set(by_gene) != missing or not missing <= original_genes:
        raise ValueError("Missing gene was not offered to the original universe")
    genes = [
        {
            "gene": gene,
            "in_original_library": True,
            "in_scored_universe": True,
            **{reason + "_forms": by_gene[gene][reason] for reason in REASONS},
            "alternate_exclusion_reasons": ";".join(
                reason for reason in REASONS[1:] if by_gene[gene][reason]
            ),
        }
        for gene in sorted(missing)
    ]
    forms.sort(key=lambda r: tuple(r.values()))
    summary = {
        "missing_cohort_genes": len(missing),
        "missing_from_original_library": len(missing - original_genes),
        "missing_from_scored_universe": len(missing - set(by_gene)),
        "form_decisions": dict(sorted(sum(by_gene.values(), Counter()).items())),
        "gene_reason_combinations": dict(
            sorted(Counter(r["alternate_exclusion_reasons"] for r in genes).items())
        ),
        "unexplained": 0,
        "genes_without_reference_rows": sorted(g for g in missing if not by_gene[g]["reference"]),
    }
    return genes, forms, summary


def write_table(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def audit(replay: Path, out: Path) -> dict:
    out.mkdir(parents=True, exist_ok=True)
    run, inputs = replay / "run", replay / "in"
    # Verify immutable inputs even when the caller reuses an existing replay.
    spec = json.loads(riu.INPUTS.read_text())
    for name, entry in spec["files"].items():
        if riu.sha256(inputs / name) != entry["sha256"]:
            raise ValueError(f"Replay input checksum differs: {name}")
    pinned = lp.stage_inputs(["isoform_cohorts.tsv"], out / "inputs")["isoform_cohorts.tsv"]
    cohorts = table(pinned)
    if sorted(table(run / "isoform_cohorts.tsv"), key=lambda r: tuple(r.values())) != sorted(
        cohorts, key=lambda r: tuple(r.values())
    ):
        raise ValueError("Replay cohort differs from the immutable recorded cohort")
    current = set(json.loads(lp.gene_list(out).read_text()))
    with gzip.open(inputs / "uniprot_human_reviewed_isoforms.tsv.gz", "rt") as fh:
        has_isoform = {
            r["Gene Names (primary)"]
            for r in csv.DictReader(fh, delimiter="\t")
            if "-" in r["Entry"] and not r["Entry"].endswith("-1")
        }
    missing = (current & has_isoform) - {r["gene"] for r in cohorts}
    original = replay / "code/u/library/l1/tag-site-architecture/library_members.json"
    genes, forms, summary = explain(
        missing,
        table(run / "isoform_universe_scored.tsv"),
        table(run / "topology_verdicts.tsv"),
        table(run / "isoform_soluble_restricted.tsv"),
        set(json.loads(original.read_text())),
    )
    summary["historical_script_commits"] = riu.COMMITS
    summary["immutable_input_manifest_sha256"] = riu.sha256(riu.INPUTS)
    summary["historical_library_sha256"] = riu.sha256(original)
    summary["recorded_cohort_sha256"] = riu.sha256(pinned)
    summary["replay_limitation"] = (
        "Original Ensembl-only shortlist sequences can differ between same-length siblings; "
        "the replay retains the committed shortlist and historical verdict joins. "
        "These reason records do not establish exact-sequence topology for a new construct."
    )
    write_table(out / "isoform_coverage_genes.tsv", genes)
    write_table(out / "isoform_coverage_forms.tsv", forms)
    summary["output_sha256"] = {
        p.name: hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted(out.glob("isoform_coverage_*.tsv"))
    }
    (out / "isoform_coverage_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def publish(out: Path, prefix: str) -> None:
    """Create an immutable task-scoped review bundle; never update shared records."""
    from deliverome_analysis.core.s3 import processed

    root = "sprint4/s4e9-orf-design-fg-panel/library-design/isoform-coverage-audit/"
    if not prefix.startswith(root) or not prefix.endswith("/") or prefix == root:
        raise ValueError("Use a fresh dated task prefix under isoform-coverage-audit/")
    sources = [Path(__file__), lp.REPO / "tests/test_isoform_coverage.py", riu.INPUTS, lp.INPUTS]
    outputs = [
        out / name
        for name in (
            "isoform_coverage_genes.tsv",
            "isoform_coverage_forms.tsv",
            "isoform_coverage_summary.json",
        )
    ]
    for path in sources + outputs:
        processed.upload_if(path, prefix + path.name, if_match=None)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replay", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--publish-prefix", help="Fresh task-scoped S3 prefix (create only)")
    args = parser.parse_args()
    print(json.dumps(audit(args.replay, args.out), indent=2))
    if args.publish_prefix:
        publish(args.out, args.publish_prefix)


if __name__ == "__main__":
    main()
