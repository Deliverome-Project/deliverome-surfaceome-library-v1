"""Overlay an abundance snapshot's review roles onto its detections_long.tsv.

``abundance_screen.py`` re-judges each screened input by the share of its total protein
*amount* in the surface reference, and writes the re-assigned ``review_role`` / ``admitted``
into ``manifest.json`` only: it copies ``detections_long.tsv`` through **unchanged**, because
no protein list changes. ``surface_specific.py`` classifies a publication as surface /
membrane / global from ``detections_long.tsv``, so run against an abundance snapshot directly
it reads the *sweep's* roles and silently answers for the count screen alone.

That failure is quiet and looks like a result: on 2026-11-18, rule B against
``2026-11-18-abundance-v4`` returned the same 44 proteins as the bare sweep with an
all-``kept`` diff, which reads as "the abundance screen changed nothing" rather than "the
roles never arrived". It changed 16 inputs, 11 of them into voting roles.

This script writes a corrected input folder: the abundance snapshot's ``detections_long.tsv``
with each row's ``review_role`` / ``admitted`` taken from the abundance ``manifest.json`` by
``study_id``, plus symlinks to the symbol-resolution tables ``surface_specific.py`` also reads
(``protein_method_matrix.tsv``, ``novel_candidates.tsv``,
``clustering/candidate_review.tsv``), which only the sweep folder carries. It fails if a
``study_id`` in the detections is absent from the manifest, rather than leaving those rows on
the sweep's role.

    uv run --locked python abundance_roles.py <sweep folder> <abundance folder> <out folder>
    uv run --locked python surface_specific.py <D> <out folder>
    uv run --locked python rule_variants.py <D> <out folder> [<baseline D>]
"""

from __future__ import annotations

import csv
import json
import sys
from collections import Counter
from pathlib import Path

#: Tables surface_specific.py reads for resolved symbols; published by the sweep, not by
#: the abundance screen.
SYMBOL_TABLES = ("protein_method_matrix.tsv", "novel_candidates.tsv")
CLUSTERING_TABLE = "clustering/candidate_review.tsv"

csv.field_size_limit(10**9)


def main() -> None:
    sweep, abundance, out = (Path(a) for a in sys.argv[1:4])
    out.mkdir(parents=True, exist_ok=True)
    (out / "clustering").mkdir(exist_ok=True)

    manifest = json.loads((abundance / "manifest.json").read_text())
    roles = {d["study_id"]: d for d in manifest["datasets"]}

    with open(abundance / "detections_long.tsv", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fields = reader.fieldnames
        rows = list(reader)

    missing = {r["study_id"] for r in rows} - set(roles)
    if missing:
        raise SystemExit(
            f"{len(missing)} study_id(s) in detections_long.tsv are absent from "
            f"{abundance.name}/manifest.json, so their roles cannot be overlaid: "
            f"{sorted(missing)[:5]}"
        )

    changed: Counter[tuple[str, str]] = Counter()
    for row in rows:
        entry = roles[row["study_id"]]
        before, after = row["review_role"], entry["review_role"]
        if before != after or row["admitted"] != str(entry["admitted"]):
            changed[(before, after)] += 1
        row["review_role"] = after
        row["admitted"] = str(entry["admitted"])

    with open(out / "detections_long.tsv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    # Symlinked rather than copied: they are the sweep's bytes, unchanged, and a copy in a
    # second place is one more thing that can go stale.
    for name in SYMBOL_TABLES + (CLUSTERING_TABLE,):
        link = out / name
        link.unlink(missing_ok=True)
        link.symlink_to((sweep / name).resolve())

    print(f"{len(rows)} detection rows; roles from {abundance.name}/manifest.json")
    for (before, after), n in sorted(changed.items(), key=lambda kv: -kv[1]):
        print(f"  {before} -> {after}: {n} rows")
    print(f"wrote {out}/detections_long.tsv (+ symlinks to the sweep's symbol tables)")


if __name__ == "__main__":
    main()
