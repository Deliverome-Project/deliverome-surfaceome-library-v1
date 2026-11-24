"""The internal-placement L1 members, as a vendorable TSV for literature enrichment.

One row per L1 member with no extracellular terminus, i.e. ``placement == "internal"``:
both termini cytoplasmic, so a terminal tag cannot be displayed and the site has to go in
a loop (or a snorkel). Built for the accessible-surfaceome tag-site enrichment sweep, which
needs to be reproducible from that repo alone, so the columns it joins on travel with it.

``topology_acc`` comes from the pointer index rather than the tag-path table, because the
whole point of that index is that a member's topology is NOT reachable from its gene's
canonical accession. A member added after the index release (SLC67A1, 2026-11-02) has
no pointer row; it is resolved here the same way the index is built -- by exact designed
sequence against v1 topology -- and the result must agree with the tag-path row's hash,
cohort and version, or the export stops. Rebuilding the index itself is a new release.

Two columns carry caveats that change how a row should be read, and both are data here
rather than a gene-name special case:

* ``gpi_rerouted = 1`` means the member reached "internal" because a GPI omega site removes
  its C-terminus, not because it is multipass. Its search space is the ectodomain, and
  ``longest_inter_tm_ecto_run_aa`` is 0 for it by construction -- that metric requires a TM
  on both sides, so it is meaningless for a single-pass protein and must not be read as
  "no site available".
* ``shed_risk`` is informational. Shedding never reroutes a placement, because no
  shedding-risk member has a usable inter-TM loop.

``longest_ecto_run_aa`` is the answer for a member whose site is not an inter-TM loop. It
is the longest extracellular run ANYWHERE, with ``longest_ecto_run_position`` saying
whether it sits before the first TM, between two, after the last, or in a protein with no
TM at all. For genuine multipass members it equals the inter-TM column; for the
``gpi_rerouted`` pair it is the ectodomain, which is where their site actually goes --
BST2 has 137 residues of it and RGMB 394, against the 0 the inter-TM column reports.

    uv run --locked python export_internal_members.py
"""

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path

HERE = Path(__file__).parent
TAG_PATHS = HERE / "l1_tag_paths.tsv"
POINTER_INDEX = HERE / "library_topology_index.tsv"
MANIFEST = HERE / "candidate_manifest_tranches.tsv"

COLUMNS = (
    "candidate_id",
    "library",
    "candidate_set",
    "species",
    "gene_symbol",
    "topology_acc",
    "topology_cohort",
    "topology_version",
    "sequence_sha256",
    "ordered_length_aa",
    "tm_count",
    "longest_inter_tm_ecto_run_aa",
    "longest_ecto_run_aa",
    "longest_ecto_run_start",
    "longest_ecto_run_position",
    "placement",
    "gpi_rerouted",
    "shed_risk",
)


def resolve_unindexed(rows: list[dict]) -> list[str]:
    """Fill ``topology_acc`` for rows the pointer index lacks, by exact designed sequence.

    Returns the candidate ids resolved this way. Each must match the tag-path row's
    ``sequence_sha256``, ``topology_cohort`` and ``topology_version``.
    """
    todo = {r["candidate_id"]: r for r in rows if not r["topology_acc"]}
    if not todo:
        return []
    from deliverome_analysis import library_design as ld
    from deliverome_analysis import topology

    libraries = {r["library"] for r in todo.values()}
    members: dict[str, dict[str, str]] = {}
    for library in sorted(libraries):
        designed = ld.designed_proteins(MANIFEST, library=library)
        members.update({cid: designed[cid] for cid in todo if cid in designed})
    if not members:
        return []
    report = ld.topology_report(members, tool=topology.TOOL_V1)
    for rec in report.to_dict("records"):
        row = todo[rec["candidate_id"]]
        if not rec["found"]:
            continue
        got = (rec["sequence_sha256"], rec["cohort"], rec["topology_version"])
        want = (row["sequence_sha256"], row["topology_cohort"], row["topology_version"])
        if got != want:
            raise SystemExit(f"{row['candidate_id']}: sequence resolve {got} != tag path {want}")
        row["topology_acc"] = rec["uniprot_acc_full"]
    return sorted(cid for cid, r in todo.items() if r["topology_acc"])


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=HERE / "l1_internal_members.tsv")
    args = ap.parse_args(argv)

    with POINTER_INDEX.open() as fh:
        acc = {r["candidate_id"]: r["topology_acc"] for r in csv.DictReader(fh, delimiter="\t")}

    rows = []
    with TAG_PATHS.open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["placement"] != "internal":
                continue
            rows.append(
                {
                    "candidate_id": r["candidate_id"],
                    "library": r["library"],
                    "candidate_set": r["candidate_set"],
                    "species": r["species"],
                    "gene_symbol": r["gene_symbol"],
                    "topology_acc": acc.get(r["candidate_id"], ""),
                    "topology_cohort": r["topology_cohort"],
                    "topology_version": r["topology_version"],
                    "sequence_sha256": r["sequence_sha256"],
                    "ordered_length_aa": r["length_aa"],
                    "tm_count": r["n_tm"],
                    "longest_inter_tm_ecto_run_aa": r["longest_inter_tm_ecto_run_aa"],
                    "longest_ecto_run_aa": r["longest_ecto_run_aa"],
                    "longest_ecto_run_start": r["longest_ecto_run_start"],
                    "longest_ecto_run_position": r["longest_ecto_run_position"],
                    "placement": r["placement"],
                    "gpi_rerouted": r["gpi_rerouted"],
                    "shed_risk": r["shed_risk"],
                }
            )

    resolved = resolve_unindexed(rows)
    for cid in resolved:
        r = next(x for x in rows if x["candidate_id"] == cid)
        print(
            f"not in {POINTER_INDEX.name}; resolved by sequence: {cid} {r['gene_symbol']} "
            f"-> {r['topology_acc']}"
        )

    missing = [r["candidate_id"] for r in rows if not r["topology_acc"]]
    if missing:
        raise SystemExit(f"no topology_acc for {len(missing)}: {missing[:5]}")

    with TAG_PATHS.open() as fh:
        internal = {
            r["sequence_sha256"]
            for r in csv.DictReader(fh, delimiter="\t")
            if r["placement"] == "internal"
        }
    exported = {r["sequence_sha256"] for r in rows}
    if exported != internal:
        raise SystemExit(f"internal sha256 not exported: {sorted(internal - exported)[:5]}")

    with args.out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    digest = hashlib.sha256(args.out.read_bytes()).hexdigest()
    print(f"wrote {args.out.name}: {len(rows)} rows")
    print(f"  sha256: {digest}")
    print(f"  distinct genes: {len({r['gene_symbol'] for r in rows})}")
    print(
        f"  gpi_rerouted (ectodomain search space, not a loop): "
        f"{sum(1 for r in rows if r['gpi_rerouted'] == '1')}"
    )
    print(f"  with a shed_risk: {sum(1 for r in rows if r['shed_risk'])}")
    thin = [
        r for r in rows if int(r["longest_inter_tm_ecto_run_aa"]) < 8 and r["gpi_rerouted"] != "1"
    ]
    print(
        f"  multipass with no inter-TM loop >= 8 aa: {len(thin)} "
        f"({', '.join(r['gene_symbol'] for r in thin)})"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
