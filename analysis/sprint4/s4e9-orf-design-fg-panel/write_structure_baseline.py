"""Regenerate the structure-coverage baseline the repo tests assert against.

`tests/test_library_structure_coverage.py` compares the committed measurements to this
file, so a new library member or a changed proteoform fails in CI rather than surfacing as
a tag site numbered against the wrong sequence. The baseline is derived, never hand-edited:
run this in the same commit that changes the library, and the diff shows exactly which
members moved.

    uv run --locked python write_structure_baseline.py
"""

from __future__ import annotations

import csv
import json
from datetime import UTC, datetime
from pathlib import Path

HERE = Path(__file__).parent
OUT = HERE / "structure_coverage_baseline.json"


def rows(name: str) -> list[dict]:
    path = HERE / name
    if not path.exists():
        raise SystemExit(f"{name} is missing — run the measurement scripts first")
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def truthy(value: str) -> bool:
    return value in ("True", "true", "1")


def main() -> int:
    l1 = rows("l1_structure_gap.tsv")
    l2 = rows("l2_structure_gap.tsv")
    match = rows("proteoform_match.tsv")
    folded = {r["name"]: truthy(r["ok"]) for r in rows("predicted_structures/fold_summary.tsv")}

    checked = {r["candidate_id"] for r in match}
    with_model = [r for r in l1 if truthy(r["has_alphafold"])]
    # A member with a model but no verdict could not be compared -- no designed CDS, or the
    # entry returned no sequence. Recorded so the completeness test stays exact.
    not_comparable = sorted(
        r["candidate_id"] for r in with_model if r["candidate_id"] not in checked
    )
    # A gap that turned out not to be one: the accession resolved on a second look, so there
    # is nothing to fold. Keeping these explicit stops them reading as unfinished work.
    no_fold_needed = sorted(
        {
            r["accession"]
            for r in l1
            if not truthy(r["has_alphafold"]) and r["accession"] and r["accession"] not in folded
        }
    )

    baseline = {
        "generated": datetime.now(UTC).strftime("%Y-%m-%d"),
        "note": (
            "Derived by write_structure_baseline.py; asserted by "
            "tests/test_library_structure_coverage.py. Do not hand-edit."
        ),
        "l1": {
            "candidates": len(l1),
            "with_model": len(with_model),
            "without_model": sum(1 for r in l1 if not truthy(r["has_alphafold"])),
            "identical": sum(1 for r in match if r["verdict"] == "identical"),
            "different": sorted(r["candidate_id"] for r in match if r["verdict"] == "different"),
            "contained": sorted(r["candidate_id"] for r in match if r["verdict"] == "contained"),
            "not_comparable": not_comparable,
            "no_fold_needed": no_fold_needed,
        },
        "l2": {
            "candidates": len(l2),
            "without_model": sum(1 for r in l2 if not truthy(r["has_alphafold"])),
        },
    }
    OUT.write_text(json.dumps(baseline, indent=2) + "\n")
    b = baseline["l1"]
    print(
        f"L1 {b['candidates']}: {b['identical']} identical, {len(b['different'])} different, "
        f"{len(b['contained'])} contained, {b['without_model']} without a model"
    )
    print(f"L2 {baseline['l2']['candidates']}: {baseline['l2']['without_model']} without a model")
    print(f"wrote {OUT.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
