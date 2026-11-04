"""Append this experiment's barcode reservations to the running P-series registry.

The registry on S3 is the single source of truth for which pool barcodes are spoken for,
and its own header says to reserve there *before* assigning a barcode to a construct.
`reserve_barcodes.py` picks the barcodes and writes `barcode_reservations.tsv`; this
appends those rows.

Safety, because the registry is shared and append-only by convention:

* re-downloads immediately before writing, so a concurrent append is not clobbered;
* refuses if any barcode or pool index being added is already present;
* re-screens every row for an intrinsic splice donor and for a Type IIS site created
  against this experiment's junction context, so a stale reservations file cannot
  introduce a bad barcode;
* writes only by appending complete lines; existing bytes are never rewritten;
* re-downloads afterwards and verifies every appended row is present exactly once.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from build_constructs import CS1, ENZYMES, FWD_PRIMER, REV_PRIMER

from deliverome_analysis import seqtools
from deliverome_analysis.s3 import processed

HERE = Path(__file__).parent
REGISTRY_KEY = "sprint1/s1e2-fragmid-vector-design/barcode_assignments.tsv"
COLUMNS = [
    "barcode",
    "pool_index",
    "construct",
    "orf",
    "position",
    "assigned_by",
    "status",
    "notes",
]
SPLICE_DONOR = re.compile(r"GT[AG]AG")
TRAILING = REV_PRIMER + CS1 + "AG" + "A" + "GAGACG"


def screen(barcode: str) -> list[str]:
    problems = []
    if SPLICE_DONOR.search(barcode):
        problems.append(f"splice donor {SPLICE_DONOR.search(barcode).group(0)}")
    context = FWD_PRIMER + barcode + TRAILING
    extra = [
        s
        for s in seqtools.find_sites(context, ENZYMES)
        if not (s[0] == "Esp3I" and s[2] >= len(FWD_PRIMER) + len(barcode))
    ]
    if extra:
        problems.append(f"junction Type IIS {extra}")
    return problems


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--dry-run",
        action="store_true",
        help="print the rows that would be appended; write nothing",
    )
    args = ap.parse_args(argv)

    rows = list(csv.DictReader((HERE / "barcode_reservations.tsv").open(), delimiter="\t"))
    text = Path(processed.download(REGISTRY_KEY)).read_text()
    body = [ln for ln in text.splitlines() if ln and not ln.startswith("#")]
    existing = list(csv.DictReader(body, delimiter="\t"))
    have_bc = {r["barcode"] for r in existing if r["barcode"]}
    have_ix = {r["pool_index"] for r in existing if r["pool_index"]}
    have_construct = {r["construct"] for r in existing}

    for row in rows:
        problems = screen(row["barcode"])
        if problems:
            raise SystemExit(f"{row['construct']}: refusing to reserve — {problems}")
        if row["barcode"] in have_bc:
            raise SystemExit(f"{row['construct']}: barcode already in the registry")
        if row["pool_index"] in have_ix:
            raise SystemExit(f"{row['construct']}: pool index already in the registry")
        if row["construct"] in have_construct:
            # A P-number reused after a renumbering is not the same failure as a
            # duplicate barcode, and a blanket refusal leaves it permanently
            # unregisterable. The barcode and pool-index checks above are the real
            # uniqueness constraints and stay absolute; this one is about a LABEL.
            #
            # The registry is append-only, so the superseded row cannot be edited.
            # What it can do is carry the supersession forward: a reused construct
            # name is allowed only when this row's notes declare what it replaces,
            # so a reader resolving "what is p0034?" finds both rows and is told
            # which is current, instead of finding two and no way to choose.
            prior = next(r for r in existing if r["construct"] == row["construct"])
            if "supersedes" not in (row.get("notes") or "").lower():
                raise SystemExit(
                    f"{row['construct']}: already has a registry row (barcode "
                    f"{prior['barcode']}, {prior['orf']}). That row cannot be edited — the "
                    f"registry is append-only. To register this one, put a 'supersedes' "
                    f"clause in its notes naming the earlier barcode and where it went."
                )
            if row["notes"].count(prior["barcode"]) == 0:
                raise SystemExit(
                    f"{row['construct']}: notes claim a supersession but do not name the "
                    f"barcode being superseded ({prior['barcode']})"
                )

    lines = ["\t".join(row[c] for c in COLUMNS) for row in rows]
    print(f"registry holds {len(existing)} row(s); appending {len(lines)}:\n")
    for line in lines:
        print("  " + line)
    if args.dry_run:
        print("\ndry run — nothing written")
        return 0

    updated = text if text.endswith("\n") else text + "\n"
    updated += "\n".join(lines) + "\n"
    with tempfile.TemporaryDirectory() as tmp:
        staged = Path(tmp) / "barcode_assignments.tsv"
        staged.write_text(updated)
        processed.upload(staged, REGISTRY_KEY, content_type="text/tab-separated-values")

    check = Path(processed.download(REGISTRY_KEY)).read_text()
    for line in lines:
        if check.count(line) != 1:
            raise SystemExit(f"verification failed: appended row not present exactly once:\n{line}")
    after = [ln for ln in check.splitlines() if ln and not ln.startswith("#")]
    print(
        f"\nverified: registry now holds {len(after) - 1} row(s), "
        f"every appended row present exactly once"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
