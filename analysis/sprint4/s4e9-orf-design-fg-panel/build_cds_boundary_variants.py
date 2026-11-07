"""Build explicitly engineered CDS variants and exact-sequence prediction worklists.

Native files and current manifests are read-only. The output is a separate design
revision: an ATG edit does not certify completeness or transfer native predictions.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from deliverome_analysis import seqtools, topology

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE.parents[2] / "scratch/becca/appris-principal"))
from assemble_candidates import ENZYMES, assemble_one  # noqa: E402
from domesticate_candidates import domesticate_one  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    for name in ("native-cds", "native-constructs", "overrides", "backbone", "snapshot", "out"):
        ap.add_argument("--" + name, type=Path, required=True)
    ap.add_argument("--manifest", type=Path, action="append", required=True)
    ap.add_argument("--backbone-sha256", required=True)
    args = ap.parse_args()
    if hashlib.sha256(args.backbone.read_bytes()).hexdigest() != args.backbone_sha256:
        raise ValueError("Backbone checksum mismatch")
    native = json.loads(args.native_cds.read_text())
    cds = dict(native)
    constructs = json.loads(args.native_constructs.read_text())
    overrides = json.loads(args.overrides.read_text())
    members = {}
    for path in args.manifest:
        with path.open() as handle:
            for r in csv.DictReader(handle, delimiter="\t"):
                if not r.get("retired_reason"):
                    members.setdefault(r["key"], set()).add(r["candidate_id"])
    bb = str(SeqIO.read(args.backbone, "genbank").seq).upper()
    bb_sites = seqtools.find_sites(bb, ENZYMES)
    if len(bb) != 7672:
        raise ValueError("Unexpected backbone length")
    snapshot = pd.read_parquet(args.snapshot)
    records, sequences = [], {}
    for key, entry in sorted(overrides.items()):
        if key not in members:
            raise ValueError(f"No affected active member: {key}")
        edit = entry["boundary_design"]
        if native[key]["cds"] != edit["original_cds"]:
            raise ValueError(f"Native input changed: {key}")
        dna = entry["cds"]
        aa = seqtools.translate(dna)
        if not dna.startswith("ATG") or not aa.endswith("*") or "*" in aa[:-1]:
            raise ValueError(f"Invalid engineered ORF: {key}")
        built = domesticate_one(entry)
        if built is None or built["sites_after"] or seqtools.translate(built["construct"]) != aa:
            raise ValueError(f"Invalid domestication: {key}")
        cid = sorted(members[key])[0]
        _, row, clean = assemble_one(bb, bb_sites, cid, key, built)
        if not clean:
            raise ValueError(f"Invalid assembled restriction pattern: {key}")
        digest = topology.sequence_sha256(aa[:-1])
        name = "cdsfix_" + digest[:16]
        sequences[name] = aa[:-1]
        hit = topology.resolve(
            {name: aa[:-1]},
            tool=topology.TOOL_V1,
            cohorts=topology.cohorts_for_species(entry["species"]),
            snapshot=snapshot,
        ).iloc[0]
        records.append(
            {
                **row,
                "affected_member_ids": ";".join(sorted(members[key])),
                "prediction_id": name,
                "action": edit["action"],
                "native_protein_sha256": edit["original_protein_sha256"],
                "designed_protein_sha256": digest,
                "topology_found": bool(hit["found"]),
                "completeness": edit["completeness"],
                "expression_validated": False,
            }
        )
        cds[key] = entry
        constructs[key] = built
    args.out.mkdir(parents=True, exist_ok=False)
    for name, value in (
        ("engineered_candidate_cds.json", cds),
        ("engineered_candidate_constructs.json", constructs),
    ):
        (args.out / name).write_text(json.dumps(value, sort_keys=True) + "\n")
    with (args.out / "cds_boundary_variants.tsv").open("w") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(records[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(records)
    missing = {r["prediction_id"] for r in records if not r["topology_found"]}
    for name, ids in (
        ("engineered_proteins.fasta", set(sequences)),
        ("missing_deeptmhmm.fasta", missing),
    ):
        (args.out / name).write_text("".join(f">{n}\n{sequences[n]}\n" for n in sorted(ids)))
    summary = {
        "engineered_records": len(records),
        "distinct_proteins": len(sequences),
        "missing_topology": len(missing),
        "missing_topology_residues": sum(len(sequences[n]) for n in missing),
        "native_records_preserved_separately": True,
        "current_manifest_changed": False,
        "input_sha256": {
            str(p): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in [
                args.native_cds,
                args.native_constructs,
                args.overrides,
                args.backbone,
                args.snapshot,
                *args.manifest,
            ]
        },
    }
    duplicate_groups = {}
    for path in args.manifest:
        groups = defaultdict(list)
        with path.open() as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if row.get("retired_reason"):
                    continue
                protein = seqtools.translate(constructs[row["key"]]["construct"]).rstrip("*")
                groups[(row["library"], topology.sequence_sha256(protein))].append(row)
        duplicate_groups[path.name] = [
            [{"candidate_id": row["candidate_id"], "key": row["key"]} for row in group]
            for group in groups.values()
            if len(group) > 1 and any(row["key"] in overrides for row in group)
        ]
    summary["within_tranche_duplicates_involving_engineered_keys"] = duplicate_groups
    (args.out / "cds_boundary_build_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
