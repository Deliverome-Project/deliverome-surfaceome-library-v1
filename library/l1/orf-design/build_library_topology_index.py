"""Which topology row describes the protein each library member actually encodes.

The Worker cannot work this out at serve time. Its topology lookups are keyed by
``uniprot_acc`` within ``human_canonical`` / ``human_isoforms``, and a library member's
protein -- its designed CDS translated -- is an Ensembl-transcript product that usually
corresponds to a UniProt *isoform* record, not the canonical one. Resolving by accession
therefore returns a different protein's topology, with a plausible length and a plausible
helix count. Resolving by sequence is correct but needs a full scan of ~71k long TEXT rows,
which is fine offline and wrong in a request.

So the join is resolved here, once, by exact sequence, and published as an explicit
pointer: per candidate, the (cohort, topology_version, uniprot_acc_full) whose stored
sequence IS the ordered sequence. The Worker then joins on that pointer and never guesses.

Emits one row per library member that has a v1 topology for its exact designed sequence.
Members with none are omitted and counted -- they need a DeepTMHMM run, not a near-miss.

    uv run --locked python build_library_topology_index.py --release topo-2026-10-14
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from deliverome_analysis.cloning import library_design as ld
from deliverome_analysis.surfaceome import topology

HERE = Path(__file__).parent
MANIFEST = HERE / "candidate_manifest_tranches.tsv"
OUT = HERE / "library_topology_index.tsv"

COLUMNS = (
    "release",
    "candidate_id",
    "library",
    "candidate_set",
    "species",
    "gene_symbol",
    "topology_tool",
    "topology_cohort",
    "topology_version",
    "topology_acc",
    "ordered_length_aa",
    "sequence_sha256",
)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--release", required=True, help="e.g. topo-2026-10-14")
    ap.add_argument("--manifest", type=Path, default=MANIFEST)
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args(argv)

    rows: list[dict] = []
    skipped: dict[str, int] = {}
    for library in ld.LIBRARIES:
        members = ld.designed_proteins(args.manifest, library=library)
        report = ld.topology_report(members, tool=topology.TOOL_V1)
        for r in report.to_dict("records"):
            if not r["found"]:
                skipped[library] = skipped.get(library, 0) + 1
                continue
            rows.append(
                {
                    "release": args.release,
                    "candidate_id": r["candidate_id"],
                    "library": library,
                    "candidate_set": r["candidate_set"],
                    "species": r["species"],
                    "gene_symbol": r["gene_symbol"],
                    "topology_tool": topology.TOOL_V1,
                    "topology_cohort": r["cohort"],
                    "topology_version": r["topology_version"],
                    "topology_acc": r["uniprot_acc_full"],
                    "ordered_length_aa": r["length"],
                    "sequence_sha256": r["sequence_sha256"],
                }
            )
        print(f"{library}: {len(members)} members, {skipped.get(library, 0)} without v1")

    with args.out.open("w", newline="") as fh:
        # Explicit LF: csv defaults to CRLF, which rewrote every line of this table on
        # each run and left git normalising it back. A table that always looks changed is
        # a table nobody reads the diff of.
        w = csv.DictWriter(
            fh,
            fieldnames=COLUMNS,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        w.writeheader()
        w.writerows(rows)

    by_cohort: dict[str, int] = {}
    for r in rows:
        by_cohort[r["topology_cohort"]] = by_cohort.get(r["topology_cohort"], 0) + 1
    print(f"\nwrote {args.out.name}: {len(rows)} members with ordered-proteoform v1 topology")
    print("by source cohort:", dict(sorted(by_cohort.items(), key=lambda kv: -kv[1])))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
