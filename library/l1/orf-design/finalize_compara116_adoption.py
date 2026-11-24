"""Combine a reviewed member map and exact constructs into an immutable candidate bundle.

Historical entries are preserved. New records may only append keys; an existing
key with different contents fails instead of silently changing an older member.
This is the pre-tag candidate stage, not a synthesis-ready order.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

from deliverome_analysis import seqtools


def read(path):
    with Path(path).open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def combine(old, selected, mapping):
    result = dict(old)
    for row in mapping:
        key, fid = row["member_key"], row["form_id"]
        value = selected[fid]
        if key in result and result[key] != value:
            raise ValueError(f"Historical record differs: {key}")
        result[key] = value
    return result


def audit(manifest, mapping, panel, cds, constructs):
    picks = {r["form_id"]: r for r in panel}
    if len(mapping) != 200 or {r["form_id"] for r in mapping} != set(picks):
        raise ValueError("Adoption must account for all 200 picks")
    active = [r for r in manifest if r["library"] == "l1" and not r.get("retired_reason")]
    orthologs = {r["candidate_id"]: r for r in active if r["set"] == "ortholog"}
    if set(orthologs) != {r["candidate_id"] for r in mapping}:
        raise ValueError("Active ortholog membership differs from selected member map")
    hashes = {}
    for r in active:
        key = r["key"]
        entry = constructs.get(key)
        if entry is None:
            raise ValueError(f"Active member lacks a constructed ORF: {key}")
        aa = seqtools.translate(entry["construct"]).rstrip("*")
        if not aa or "*" in aa:
            raise ValueError(f"Invalid constructed protein: {key}")
        digest = hashlib.sha256(aa.encode()).hexdigest()
        hashes.setdefault(digest, []).append(r["candidate_id"])
    for r in mapping:
        key = r["member_key"]
        expected = picks[r["form_id"]]["sequence_sha256"]
        aa = seqtools.translate(cds[key]["cds"])
        if not aa.endswith("*") or "*" in aa[:-1]:
            raise ValueError(f"Invalid selected CDS: {key}")
        if hashlib.sha256(aa[:-1].encode()).hexdigest() != expected:
            raise ValueError(f"Selected CDS identity mismatch: {key}")
        if hashes.get(expected) != [r["candidate_id"]]:
            raise ValueError(f"Selected protein missing or duplicated: {key}")
    blocked = [
        {
            "candidate_id": r["candidate_id"],
            "form_id": r["form_id"],
            "start_codon": cds[r["member_key"]]["cds"][:3],
        }
        for r in mapping
        if cds[r["member_key"]]["cds"][:3] != "ATG"
    ]
    return {
        "selected_members": len(mapping),
        "active_l1_members_checked": len(active),
        "selected_protein_duplicates": 0,
        "stage": "pre-tag candidate constructs",
        "non_atg_start_blockers": blocked,
        "expression_ready": False,
        "tags": "deferred",
        "barcodes": "not reserved",
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    for name in (
        "panel",
        "manifest",
        "member-map",
        "old-cds",
        "old-constructs",
        "selected-cds",
        "selected-constructs",
        "out",
    ):
        ap.add_argument("--" + name, type=Path, required=True)
    args = ap.parse_args()
    mapping = read(args.member_map)
    cds = combine(
        json.loads(args.old_cds.read_text()), json.loads(args.selected_cds.read_text()), mapping
    )
    constructs = combine(
        json.loads(args.old_constructs.read_text()),
        json.loads(args.selected_constructs.read_text()),
        mapping,
    )
    summary = audit(read(args.manifest), mapping, read(args.panel), cds, constructs)
    args.out.mkdir(parents=True, exist_ok=True)
    for name, value in (
        ("candidate_cds_compara116.json", cds),
        ("candidate_constructs_compara116.json", constructs),
    ):
        (args.out / name).write_text(json.dumps(value, sort_keys=True) + "\n")
    summary["input_sha256"] = {
        name: hashlib.sha256(path.read_bytes()).hexdigest()
        for name, path in vars(args).items()
        if name != "out"
    }
    summary["output_sha256"] = {
        p.name: hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted(args.out.glob("candidate_*_compara116.json"))
    }
    (args.out / "compara116_adoption_validation.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
