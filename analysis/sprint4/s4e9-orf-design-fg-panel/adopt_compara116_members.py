"""Plan and review Compara 116 panel adoption without overwriting a library record.

IDs belong to the tranche namespace, not the abandoned pre-tranche candidate manifest.
Every existing row remains in place; exact L1 proteins keep their IDs. L2 proteins
receive fresh L1 IDs and their L2 rows are retired as moved_to:<new ID>, following
redraw_panels.py. New proteins require genuine assembled construction rows.

The command writes a fresh review directory only. It never publishes or modifies
the input manifest, CDS, panels, reservations, or canonical S3 records.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

from deliverome_analysis import seqtools

HERE = Path(__file__).parent
ID = re.compile(r"^(l[12])-p([0-9]+)$")
CONSTRUCTION_COLUMNS = (
    "transcript",
    "route",
    "variant_residues",
    "orf_nt",
    "vector_bp",
    "sites_removed",
    "edits",
    "vector_sites",
    "esp3i_sites",
    "insert_nt",
    "has_igk_leader",
    "tag_cassette",
    "barcode",
)
RETIRED = "retired_reason"


def read(path: Path) -> list[dict]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def active(row: dict) -> bool:
    return not row.get(RETIRED, "").strip()


def protein_hash(entry: dict) -> str:
    protein = seqtools.translate(entry["cds"]).rstrip("*")
    if not protein or "*" in protein:
        raise ValueError("CDS has an empty protein or internal stop")
    return hashlib.sha256(protein.encode()).hexdigest()


def validate_ids(rows: list[dict]) -> set[str]:
    ids = [r["candidate_id"] for r in rows]
    if len(ids) != len(set(ids)):
        raise ValueError("Manifest candidate IDs are not unique")
    for row in rows:
        match = ID.fullmatch(row["candidate_id"])
        if not match or match[1] != row["library"]:
            raise ValueError(f"ID/library mismatch: {row['candidate_id']}")
    return set(ids)


def plan(
    panel: list[dict],
    manifest: list[dict],
    cds: dict,
    reserved_ids: set[str] | None = None,
    recorded_allocations: list[dict] | None = None,
) -> tuple[list[dict], dict]:
    """Allocate only above every live/retired/reserved L1 tranche ID."""
    all_ids = validate_ids(manifest)
    ledger = recorded_allocations or []
    if len({r["candidate_id"] for r in ledger}) != len(ledger) or len(
        {r["form_id"] for r in ledger}
    ) != len(ledger):
        raise ValueError("Allocation ledger IDs/form IDs are not unique")
    recorded = {r["form_id"]: r for r in ledger}
    reserved = all_ids | (reserved_ids or set()) | {r["candidate_id"] for r in ledger}
    # moved_to references are current tranche identities; prior_candidate_id is
    # deliberately excluded because it records the abandoned pre-tranche namespace.
    reserved |= {
        r[RETIRED].removeprefix("moved_to:")
        for r in manifest
        if r.get(RETIRED, "").startswith("moved_to:")
    }
    if any(ID.fullmatch(cid) is None for cid in reserved):
        raise ValueError("Invalid reserved tranche ID")
    next_n = 1 + max(
        (int(ID.fullmatch(cid)[2]) for cid in reserved if ID.fullmatch(cid)[1] == "l1"), default=0
    )
    max_before = next_n - 1
    by_protein = defaultdict(list)
    unavailable = []
    for row in manifest:
        if not active(row):
            continue
        if row["key"] not in cds:
            if row["set"] == "ortholog":
                raise ValueError(f"Active ortholog has no frozen CDS: {row['candidate_id']}")
            unavailable.append(row["candidate_id"])
            continue
        by_protein[protein_hash(cds[row["key"]])].append(row)
    if len({p["form_id"] for p in panel}) != len(panel):
        raise ValueError("Panel form IDs are not unique")
    if len({p["sequence_sha256"] for p in panel}) != len(panel):
        raise ValueError("Panel proteins are not unique")
    mapping = []
    for pick in panel:
        hits = by_protein.get(pick["sequence_sha256"], [])
        l1 = [r for r in hits if r["library"] == "l1"]
        if len(l1) > 1:
            raise ValueError(f"Multiple active L1 exact-protein members: {pick['form_id']}")
        if l1 and (l1[0]["set"] != "ortholog" or l1[0]["species"] != pick["species"]):
            raise ValueError(f"Selected protein duplicates a non-panel member: {pick['form_id']}")
        # Existing duplicate L2 rows are an unresolved identity decision, not an
        # invitation to choose a gene/construct by accession or table order.
        l2 = [r for r in hits if r["library"] == "l2"]
        if not l1 and len(l2) > 1:
            raise ValueError(f"Multiple active L2 exact-protein members: {pick['form_id']}")
        src = l1[0] if l1 else l2[0] if l2 else None
        if src and (src["species"] != pick["species"] or src["set"] != "ortholog"):
            raise ValueError(f"Existing member species/set differs: {pick['form_id']}")
        key = src["key"] if src else "ort:" + pick["form_id"]
        versioned = False
        if not src and any(r["key"] == key for r in manifest):
            species, gene, accession = pick["form_id"].split(":", 2)
            key = f"ort:{species}:{gene}:compara116:{accession}"
            versioned = True
            if any(r["key"] == key for r in manifest):
                raise ValueError(f"Versioned key already exists: {key}")
        if l1:
            cid, action = src["candidate_id"], "retained_l1"
        else:
            old = recorded.get(pick["form_id"])
            if old:
                if old["sequence_sha256"] != pick["sequence_sha256"] or old["member_key"] != key:
                    raise ValueError(f"Allocation ledger cannot be rebound: {pick['form_id']}")
                cid = old["candidate_id"]
                if not cid.startswith("l1-p") or cid in all_ids:
                    raise ValueError(f"Recorded allocation collides with manifest: {cid}")
            else:
                cid = f"l1-p{next_n:04d}"
                next_n += 1
                if cid in reserved:
                    raise ValueError(f"Allocated ID collides with reservation: {cid}")
                reserved.add(cid)
            action = "promoted_l2" if src else "new_form"
        mapping.append(
            {
                "form_id": pick["form_id"],
                "species": pick["species"],
                "human_symbol": pick["human_symbol"],
                "sequence_sha256": pick["sequence_sha256"],
                "candidate_id": cid,
                "member_key": key,
                "action": action,
                "versioned_key_collision": str(versioned),
                "source_candidate_id": src["candidate_id"] if src else "",
                "advisory_existing_candidate_ids": pick.get("existing_candidate_ids", ""),
            }
        )
    counts = dict(sorted(Counter(r["action"] for r in mapping).items()))
    return mapping, {
        "selected_forms": len(mapping),
        "actions": counts,
        "max_reserved_l1_before": max_before,
        "fresh_ids": [r["candidate_id"] for r in mapping if r["action"] != "retained_l1"],
        "legacy_pre_tranche_ids_are_separate_namespace": True,
        "active_nonortholog_members_without_cds": unavailable,
    }


def reconcile(manifest: list[dict], mapping: list[dict], built_rows: list[dict]) -> list[dict]:
    """Preserve existing rows, retire replaced/moved rows, append genuine builds."""
    cols = list(manifest[0])
    by_id = {r["candidate_id"]: r for r in manifest}
    built = {r["key"]: r for r in built_rows}
    if len(built) != len(built_rows):
        raise ValueError("Built row keys are not unique")
    keep = {r["candidate_id"] for r in mapping if r["action"] == "retained_l1"}
    moves = {
        r["source_candidate_id"]: r["candidate_id"] for r in mapping if r["action"] == "promoted_l2"
    }
    result = []
    for original in manifest:
        row = dict(original)
        if active(row) and row["library"] == "l1" and row["set"] == "ortholog":
            if row["candidate_id"] not in keep:
                row[RETIRED] = "dropped_from_panel:compara116"
        if row["candidate_id"] in moves:
            row[RETIRED] = "moved_to:" + moves[row["candidate_id"]]
        result.append(row)
    for entry in mapping:
        if entry["action"] == "retained_l1":
            continue
        if entry["action"] == "promoted_l2":
            row = dict(by_id[entry["source_candidate_id"]])
        else:
            if entry["member_key"] not in built:
                raise ValueError(f"No assembled row for {entry['member_key']}")
            construction = built[entry["member_key"]]
            if construction.get("sequence_sha256") != entry["sequence_sha256"]:
                raise ValueError(f"Assembled row protein hash differs: {entry['member_key']}")
            absent = [c for c in CONSTRUCTION_COLUMNS if not str(construction.get(c, "")).strip()]
            if absent:
                raise ValueError(f"Incomplete assembled row {entry['member_key']}: {absent}")
            if construction.get("candidate_id") not in (None, "", entry["candidate_id"]):
                raise ValueError(f"Assembled row candidate ID differs: {entry['member_key']}")
            if (
                construction.get("species") != entry["species"]
                or construction.get("set") != "ortholog"
            ):
                raise ValueError(f"Assembled row species/set differs: {entry['member_key']}")
            row = {c: str(construction.get(c, "")) for c in cols}
        row.update(
            candidate_id=entry["candidate_id"],
            key=entry["member_key"],
            library="l1",
            panel="ortholog_panel",
            prior_candidate_id=(entry["source_candidate_id"] or entry["candidate_id"]),
            retired_reason="",
        )
        result.append(row)
    validate_ids(result)
    selected = {r["candidate_id"] for r in mapping}
    actual = {
        r["candidate_id"]
        for r in result
        if active(r) and r["library"] == "l1" and r["set"] == "ortholog"
    }
    if selected != actual:
        raise ValueError("Adopted active L1 ortholog IDs differ from selected map")
    return result


def write(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panel", type=Path, default=HERE / "ortholog_panel_compara116.tsv")
    parser.add_argument("--manifest", type=Path, default=HERE / "candidate_manifest_tranches.tsv")
    parser.add_argument("--existing-cds", type=Path, required=True)
    parser.add_argument("--built-rows", type=Path)
    parser.add_argument("--reserved-ids", type=Path, action="append", default=[])
    parser.add_argument(
        "--allocation-ledger", type=Path, default=HERE / "compara116_id_allocations.tsv"
    )
    parser.add_argument("--plan-only", action="store_true")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    manifest = read(args.manifest)
    reservations = {r["candidate_id"] for path in args.reserved_ids for r in read(path)}
    ledger = read(args.allocation_ledger) if args.allocation_ledger.exists() else []
    mapping, summary = plan(
        read(args.panel), manifest, json.loads(args.existing_cds.read_text()), reservations, ledger
    )
    result = (
        None
        if args.plan_only
        else reconcile(manifest, mapping, read(args.built_rows) if args.built_rows else [])
    )
    args.out.mkdir(parents=True, exist_ok=False)
    write(args.out / "compara116_member_map.tsv", mapping)
    write(
        args.out / "compara116_id_allocations.tsv",
        [r for r in mapping if r["action"] != "retained_l1"],
    )
    if result is not None:
        write(args.out / "candidate_manifest_tranches_review.tsv", result)
        summary["retired_old_l1_orthologs"] = sum(
            r.get(RETIRED) == "dropped_from_panel:compara116" for r in result
        )
    summary["manifest_ready"] = result is not None
    paths = [args.panel, args.manifest, args.existing_cds, *args.reserved_ids]
    if args.built_rows:
        paths.append(args.built_rows)
    if args.allocation_ledger.exists():
        paths.append(args.allocation_ledger)
    summary["input_sha256"] = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    (args.out / "compara116_member_adoption_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
