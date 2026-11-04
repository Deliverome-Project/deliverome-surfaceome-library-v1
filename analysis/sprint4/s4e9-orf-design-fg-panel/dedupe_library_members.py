"""Retire library members whose ordered protein duplicates another member's.

A member's protein is its designed CDS translated (``library_design.designed_proteins``).
Two members that translate to the same protein give the screen one protein under two
barcodes: the readout cannot tell them apart, the pool spends a slot on nothing, and their
shared tract is the whole ORF, so their barcodes swap freely (``l1_swap_risk.py``). Every
identity upstream -- gene symbol, UniProt accession, transcript -- was different for each
pair found here, which is how each one got past the stage that drew it. This compares the
only identity that matters, the protein sequence.

**Retired, not removed.** A retired member keeps its ``candidate_id`` and its row, with
``retired_reason`` set; ids are never renumbered or reused, so every table already keyed
on them (tag paths, topology index, the SAB deck) still joins. Readers in
``library_design`` skip retired rows.

**Which member is kept** follows set precedence -- a library reference before an
expression-rescue arm, isoform, ortholog or giant truncation; L1 before L2 -- then fewer
residues off the source record, then the lower id. For each L1 group the result was checked
by hand and its reason is recorded in ``L1_DECISIONS``; an L1 group without one is an
error, so a new duplicate cannot be retired on the rule alone. L2 groups are retired by
the rule and categorised.

    uv run --locked python dedupe_library_members.py           # report, write nothing
    uv run --locked python dedupe_library_members.py --write    # update the manifest
    uv run --locked python dedupe_library_members.py --check    # fail on any active duplicate
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

from deliverome_analysis import library_design as ld
from deliverome_analysis import seqtools
from deliverome_analysis.s3 import processed

HERE = Path(__file__).parent
#: On a kept member: the other gene symbols whose retired member encoded the same protein,
#: so a readout for OPN1MW is read as OPN1MW / OPN1MW2 / OPN1MW3.
ALIAS_COLUMN = "also_encodes"
MANIFEST = HERE / "candidate_manifest_tranches.tsv"
OUT = HERE / "library_duplicate_members.tsv"
SUMMARY = HERE / "library_duplicate_summary.json"

SET_ORDER = {
    "library_reference": 0,
    "expression_rescue": 1,
    "isoform": 2,
    "ortholog": 3,
    "giant_truncation": 4,
}

#: Every L1 member retired here, the member it duplicates, and why. Checked by hand.
L1_DECISIONS: dict[str, tuple[str, str, str]] = {
    "l1-p2387": (
        "l1-p0017",
        "isoform_is_reference",
        "ACHE: the GPI rule moved the reference onto P22303-2 (the anchored H subunit); the "
        "isoform universe still took the MANE/UniProt-canonical T subunit as its reference, "
        "so it offered P22303-2 as an alternate. Both build ENST00000302913.",
    ),
    "l1-p2411": (
        "l1-p0263",
        "allele_variant_of_reference",
        "CD177: Q8N6Q3 differs from the GRCh38 reference at 3 residues (V184G, D204N, M205R), "
        "which UniProt annotates as sequence conflicts from one cDNA (CAC83724), not a splice "
        "form. The resolver's allele-variant route then built it from the reference's own "
        "transcript, so the member was the reference. Building the 3 residues would order a "
        "cDNA artefact, not an isoform.",
    ),
    "l1-p2585": (
        "l1-p1225",
        "ortholog_identical_to_human",
        "MMP16: the cyno ortholog G7PC63 is identical in protein to human MMP16; the panel "
        "exists to test sequence drift and this one has none. The near band admitted it "
        "because no step excluded 100% identity.",
    ),
    "l1-p1350": (
        "l1-p1349",
        "identical_paralog",
        "OPN1MW2 is a tandem copy of OPN1MW in the X-linked opsin array and encodes the same "
        "protein; keep OPN1MW, the array's first copy.",
    ),
    "l1-p1351": (
        "l1-p1349",
        "identical_paralog",
        "OPN1MW3, as OPN1MW2: the same protein as OPN1MW.",
    ),
    "l1-p1406": (
        "l1-p1404",
        "identical_paralog",
        "OR2A42 encodes the same protein as OR2A1 (a segmental duplication); keep OR2A1.",
    ),
    "l1-p2318": (
        "l1-p2317",
        "identical_paralog",
        "UPK3BL2 encodes the same protein as UPK3BL1 (a segmental duplication); keep UPK3BL1.",
    ),
}

CATEGORY_NOTE = {
    "isoform_is_reference": "isoform panel offered the selected reference as an alternate form",
    "allele_variant_of_reference": "resolved within the variant tolerance onto the same protein",
    "ortholog_identical_to_human": "ortholog protein identical to the human member",
    "identical_paralog": "a paralogous gene encoding the same protein",
    "identical_to_other_gene": "another gene's member (paralog or read-through) encodes the same protein",
    "duplicate_isoform": "two isoform keys resolving to the same protein",
}


def precedence(r: dict) -> tuple:
    vr = int(r["variant_residues"]) if r["variant_residues"].strip().isdigit() else 0
    return (r["library"] != ld.LIBRARY_L1, SET_ORDER.get(r["set"], 99), vr, r["candidate_id"])


def category(kept: dict, dup: dict) -> str:
    if dup["set"] == "ortholog" and dup["species"] != "human":
        return "ortholog_identical_to_human"
    if ld.gene_from_key(kept["key"]).upper() != ld.gene_from_key(dup["key"]).upper():
        return "identical_to_other_gene"
    if dup["set"] == "isoform" and kept["set"] == "isoform":
        return "duplicate_isoform"
    if dup["route"] == "allele variant" or dup["variant_residues"] not in ("", "0"):
        return "allele_variant_of_reference"
    return "isoform_is_reference"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--write", action="store_true")
    ap.add_argument("--check", action="store_true", help="exit 1 if any active member duplicates")
    args = ap.parse_args(argv)

    with open(processed.download(ld.CDS_KEY)) as fh:
        cds = json.load(fh)
    with MANIFEST.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        cols = list(reader.fieldnames or [])
        rows = list(reader)
    for col in (ld.RETIRED_COLUMN, ALIAS_COLUMN):
        if col not in cols:
            cols.append(col)
    by_id = {r["candidate_id"]: r for r in rows}

    groups: dict[str, list[dict]] = defaultdict(list)
    for r in rows:
        # A row moved to another id (repick_panel_members.py) is the same candidate, not a
        # duplicate of it, and keeps its own reason.
        if r["key"] in cds and not (r.get(ld.RETIRED_COLUMN) or "").startswith("moved_to:"):
            groups[seqtools.translate(cds[r["key"]]["cds"]).rstrip("*")].append(r)
    dups = sorted(
        (sorted(g, key=precedence) for g in groups.values() if len(g) > 1),
        key=lambda g: g[0]["candidate_id"],
    )

    if args.check:
        active = [
            [r["candidate_id"] for r in g if ld.is_active(r)]
            for g in dups
            if sum(ld.is_active(r) for r in g) > 1
        ]
        if active:
            raise SystemExit(f"{len(active)} protein(s) carried by >1 active member: {active}")
        print(f"no active duplicates ({len(dups)} identical groups, all resolved)")
        return 0

    table, retire = [], {}
    for g in dups:
        kept = g[0]
        for r in g[1:]:
            cid = r["candidate_id"]
            cat = category(kept, r)
            note = CATEGORY_NOTE[cat]
            if r["library"] == ld.LIBRARY_L1:
                if cid not in L1_DECISIONS:
                    raise SystemExit(
                        f"{cid} ({r['key']}) duplicates {kept['candidate_id']}: "
                        "add an L1_DECISIONS entry after checking it by hand"
                    )
                want, cat, note = L1_DECISIONS[cid]
                if want != kept["candidate_id"]:
                    raise SystemExit(
                        f"{cid}: L1_DECISIONS keeps {want}, the rule keeps "
                        f"{kept['candidate_id']}; reconcile them"
                    )
            retire[cid] = f"duplicate_of:{kept['candidate_id']} ({cat})"
            table.append(
                {
                    "candidate_id": cid,
                    "key": r["key"],
                    "library": r["library"],
                    "set": r["set"],
                    "kept_candidate_id": kept["candidate_id"],
                    "kept_key": kept["key"],
                    "kept_library": kept["library"],
                    "category": cat,
                    "protein_aa": len(seqtools.translate(cds[r["key"]]["cds"]).rstrip("*")),
                    "rationale": note,
                }
            )
    stale = set(L1_DECISIONS) - set(retire)
    if stale:
        raise SystemExit(f"L1_DECISIONS names members that no longer duplicate: {sorted(stale)}")

    def active_count(lib: str, retired: set[str]) -> int:
        return sum(1 for r in rows if r["library"] == lib and r["candidate_id"] not in retired)

    already = {r["candidate_id"] for r in rows if not ld.is_active(r)}
    summary = {
        "identical_groups": len(dups),
        "members_in_groups": sum(len(g) for g in dups),
        "groups_by_tranche": dict(
            Counter("+".join(sorted({r["library"] for r in g})) for g in dups)
        ),
        "retired": len(retire),
        "retired_by_tranche_and_category": {
            lib: dict(Counter(t["category"] for t in table if t["library"] == lib))
            for lib in ld.LIBRARIES
        },
        "members_before": {lib: active_count(lib, already - set(retire)) for lib in ld.LIBRARIES},
        "members_after": {lib: active_count(lib, already | set(retire)) for lib in ld.LIBRARIES},
        "l1_with_designed_cds_after": sum(
            1
            for r in rows
            if r["library"] == ld.LIBRARY_L1 and r["key"] in cds and r["candidate_id"] not in retire
        ),
    }
    print(json.dumps(summary, indent=1))
    for t in table:
        if t["library"] == ld.LIBRARY_L1:
            print(
                f"  retire {t['candidate_id']} {t['key']:<28} -> keep {t['kept_candidate_id']} "
                f"{t['kept_key']}  [{t['category']}]"
            )
    if not args.write:
        print("\ndry run -- re-run with --write")
        return 0

    aliases: dict[str, set[str]] = defaultdict(set)
    for t in table:
        kept_gene, gene = ld.gene_from_key(t["kept_key"]), ld.gene_from_key(t["key"])
        if gene.upper() != kept_gene.upper() and t["set"] != "ortholog":
            aliases[t["kept_candidate_id"]].add(gene)
    for r in rows:
        r[ld.RETIRED_COLUMN] = retire.get(r["candidate_id"], r.get(ld.RETIRED_COLUMN) or "")
        r[ALIAS_COLUMN] = ";".join(sorted(aliases.get(r["candidate_id"], ())))
    assert len({r["candidate_id"] for r in rows}) == len(rows) == len(by_id)
    with MANIFEST.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    with OUT.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(table[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(table)
    SUMMARY.write_text(json.dumps(summary, indent=1) + "\n")
    print(f"\nwrote {MANIFEST.name} ({len(retire)} retired), {OUT.name}, {SUMMARY.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
