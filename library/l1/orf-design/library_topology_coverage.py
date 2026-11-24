"""Which library members have a v1 DeepTMHMM topology for the exact sequence we order.

Answers one question per run, for one tranche, with one predictor, and says so in every
output column -- the three identities (library tranche / tool version / proteoform) have
been conflated before, each time producing a plausible wrong answer.

The protein is the designed CDS translated, which is an Ensembl-transcript product and
usually corresponds to a UniProt *isoform* rather than the canonical record. So a topology
row fetched by canonical accession is frequently about a different protein. Matching is on
the sequence throughout; see :mod:`deliverome_analysis.library_design`.

    uv run --locked python library_topology_coverage.py --library l1
    uv run --locked python library_topology_coverage.py --library l2 --tool deeptmhmm2
"""

from __future__ import annotations

import argparse
from pathlib import Path

from deliverome_analysis import library_design as ld
from deliverome_analysis import topology

HERE = Path(__file__).parent
MANIFEST = HERE / "candidate_manifest_tranches.tsv"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--library", choices=ld.LIBRARIES, required=True)
    ap.add_argument(
        "--tool",
        choices=(topology.TOOL_V1, topology.TOOL_V2),
        default=topology.PREFERRED_TOOL,
        help=f"DeepTMHMM tool version, NOT a library tranche (default {topology.PREFERRED_TOOL})",
    )
    ap.add_argument("--manifest", type=Path, default=MANIFEST)
    args = ap.parse_args(argv)

    members = ld.designed_proteins(args.manifest, library=args.library)
    print(f"{args.library}: {len(members)} members with a designed CDS")

    report = ld.topology_report(members, tool=args.tool)
    cov = ld.coverage(report)
    print(f"\ncoverage for tool={args.tool}:")
    print(cov.to_string(index=False))

    gap = report.loc[~report["found"]]
    out = HERE / f"{args.library}_{args.tool}_topology_gap.tsv"
    keep = [
        c
        for c in (
            "candidate_id",
            "library",
            "candidate_set",
            "species",
            "gene_symbol",
            "tool",
            "length",
            "sequence_sha256",
        )
        if c in gap.columns
    ]
    gap[keep].to_csv(out, sep="\t", index=False)
    print(
        f"\n{len(gap)} members need a {args.tool} run "
        f"({gap['sequence_sha256'].nunique()} distinct sequences) -> {out.name}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
