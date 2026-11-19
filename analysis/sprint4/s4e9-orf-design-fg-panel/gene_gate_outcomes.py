"""Every cohort gene's outcome at the L1 gene-list gates, as a table a reviewer can read.

The design brief's section 2 states each gate as a rule and gives the counts that come out of
it. A reader who wants to know why one particular gene is absent, or which gate did the work,
cannot get that from a count. This writes it per gene: the gate that excluded it, or the tier
it entered the library by.

It reads the S1E3 filter's own ``GATES`` and ``first_failing_gate``, so the table cannot
disagree with the membership it documents; ``tests/test_gene_gate_outcomes.py`` asserts that
agreement on every gene.

The filter's inputs (``cohort.jsonl`` and its siblings) are fetched rather than committed, so
this runs where they are present: in a reviewer replay, which points the S3 helper at the
offline bundle, or after a deliberate refresh. It does not fetch anything itself — the cohort
endpoint is per-gene and fetching it in bulk has taken the public API down before.

    uv run --locked python gene_gate_outcomes.py
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

HERE = Path(__file__).parent
S1E3 = HERE.parents[1] / "sprint1/s1e3-surface-library-filters"
OUT = HERE / "gene_gate_outcomes.tsv"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--s1e3",
        type=Path,
        default=S1E3,
        help="the S1E3 folder with the filter's inputs staged. Defaults to the one in this "
        "repository; `library_pipeline.py verify --stage gene_list` stages a populated copy "
        "in its work directory, which is where this usually runs.",
    )
    a = ap.parse_args(argv)
    s1e3 = a.s1e3.resolve()
    if not (s1e3 / "cohort.jsonl").exists():
        raise SystemExit(
            f"{s1e3 / 'cohort.jsonl'} is missing. The filter's inputs are built, not committed. "
            "Run `library_pipeline.py verify --stage gene_list`, which builds them offline from "
            "the pinned release deposit, and pass its work directory's s1e3 folder as --s1e3."
        )
    sys.path.insert(0, str(s1e3))
    import final_library as F  # importing runs the filter, which is how it is meant to run

    tiers = {r["gene"]: F.tier(r) for r in F.recs if F.passes(r)}
    gate_names = [name for name, _pred in F.GATES]
    rows = []
    for r in F.recs:
        gene = r["gene"]
        gate = F.first_failing_gate(r)
        # One column per gate, so "the genes passing this filter" is a single-column filter.
        # The gates are applied in order and stop at the first failure, so a gene that failed
        # an earlier gate was never judged by the later ones: that is `not_reached`, which is
        # not the same claim as `fail`.
        per_gate, reached = {}, True
        for name in gate_names:
            if not reached:
                per_gate[name] = "not_reached"
            elif name == gate:
                per_gate[name] = "fail"
                reached = False
            else:
                per_gate[name] = "pass"
        rows.append(
            {
                "gene": gene,
                "uniprot_accession": r.get("uniprot") or "",
                "outcome": "excluded" if gate else "in_library",
                "tier": tiers.get(gene, ""),
                "gate_failed": gate or "",
                **per_gate,
            }
        )
    rows.sort(key=lambda d: (d["outcome"], d["gate_failed"], d["gene"]))

    with open(OUT, "w", newline="") as fh:
        # Explicit LF: csv defaults to \r\n, which git stores as \n, so the file would differ
        # from what this script writes and never reproduce byte for byte.
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    kept = sum(1 for d in rows if d["outcome"] == "in_library")
    print(f"wrote {OUT.name}: {len(rows)} cohort genes, {kept} in the library")
    by_gate: dict[str, int] = {}
    for d in rows:
        if d["gate_failed"]:
            by_gate[d["gate_failed"]] = by_gate.get(d["gate_failed"], 0) + 1
    for name, _pred in F.GATES:  # the brief's order, so the table reads against section 2
        if name in by_gate:
            print(f"  excluded by {name}: {by_gate[name]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
