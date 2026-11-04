"""Rebuild ``barcode_reservations.tsv`` as the COMPLETE set still owed to the registry.

``reserve_barcodes.py`` writes only the rows it drew on that run, and
``append_reservations.py`` reads nothing else -- so after a run that draws one arm, the
pending file names one arm and appending registers one arm. The other reservations are not
lost, they are simply invisible to the step whose job is to publish them. On 2026-10-18 the
file held a single row for ``AVPR2_plain_alfa``, a construct cut from the batch, while
**eight** built constructs had never reached the shared registry at all.

This derives the pending set instead of accumulating it: every construct in this folder's
barcode tables whose barcode the registry does not already hold. Deriving it means the file
cannot drift from what was built, and re-running after an append correctly yields nothing.

The registry is the shared s1e2 pool record, so a barcode missing from it is a barcode
another experiment could draw a second time. That is the whole risk being closed.

    uv run --locked python collect_reservations.py            # show what is owed
    uv run --locked python collect_reservations.py --write    # rewrite the pending file
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from reserve_barcodes import REGISTRY_KEY

from deliverome_analysis.s3 import processed

HERE = Path(__file__).parent
#: (table, gate column) -- the gate drops arms parked out of the order batch.
TABLES = (
    ("barcodes.tsv", None),
    ("barcodes_loop.tsv", None),
    ("barcodes_tfrc.tsv", "in_order_batch"),
    ("barcodes_gpcr.tsv", None),
    ("barcodes_sp_swap.tsv", None),
)


def registry_rows() -> list[dict]:
    text = Path(processed.download(REGISTRY_KEY)).read_text()
    body = [ln for ln in text.splitlines() if ln and not ln.startswith("#")]
    return list(csv.DictReader(body, delimiter="\t"))


def registry_barcodes() -> set[str]:
    return {r["barcode"].strip().upper() for r in registry_rows() if r["barcode"]}


def supersession_note(plasmid: str, prior: dict, now_at: str | None) -> str:
    """The clause ``append_reservations.py`` requires before a P-number may be reused.

    A renumbering can leave the registry holding a row that attaches a P-number to a
    construct that no longer carries it. The file is append-only, so that row cannot be
    corrected in place -- but a new row can say what it replaces, which is the difference
    between a reader finding two answers and being told which one is current.

    Derived, never typed: the prior barcode comes from the registry row being superseded
    and ``now_at`` from this folder's own tables, so the note cannot claim a supersession
    the data does not support.
    """
    where = f"that barcode is now {now_at}" if now_at else "that construct was never built"
    return (
        f"supersedes the earlier {plasmid} row (barcode {prior['barcode']}, "
        f"{prior['orf']}, pool index {prior['pool_index']}) -- {where}; "
        f"renumbered 2026-10-18, registry is append-only so the earlier row stands"
    )


def local_rows() -> list[dict]:
    out = []
    for table, gate in TABLES:
        path = HERE / table
        if not path.exists():
            continue
        for r in csv.DictReader(path.open(), delimiter="\t"):
            if gate and r.get(gate) != "1":
                continue
            bc = (r.get("barcode") or "").strip().upper()
            if not bc or set(bc) - set("ACGT"):
                continue  # a released/reassigned row carries prose, not a barcode
            arm = r.get("arm") or r.get("variant") or ""
            out.append(
                {
                    "barcode": bc,
                    "pool_index": r.get("pool_index", ""),
                    "construct": r["plasmid"],
                    "orf": r["gene"],
                    "position": "",
                    "assigned_by": "s1e2 pool",
                    "status": "designed",
                    "notes": f"{arm or 'nterm'}; analysis/sprint4/s4e9-orf-design-fg-panel; "
                    f"provisional P-number, not yet built",
                    "_table": table,
                }
            )
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--write", action="store_true")
    args = ap.parse_args(argv)

    reg = registry_rows()
    have = {r["barcode"].strip().upper() for r in reg if r["barcode"]}
    by_construct = {r["construct"].strip(): r for r in reg}
    rows = local_rows()
    owed = [r for r in rows if r["barcode"] not in have]

    # Where this folder reuses a P-number the registry already knows under a different
    # barcode, say so in the row rather than leaving it permanently unregisterable.
    bc_to_plasmid = {r["barcode"]: r["construct"] for r in rows}
    for r in owed:
        prior = by_construct.get(r["construct"])
        if prior and prior["barcode"].strip().upper() != r["barcode"]:
            r["notes"] += "; " + supersession_note(
                r["construct"], prior, bc_to_plasmid.get(prior["barcode"].strip().upper())
            )

    print(f"registry holds {len(have)} barcode(s); this folder has {len(rows)} in the order batch")
    print(f"owed to the registry: {len(owed)}")
    for r in sorted(owed, key=lambda x: x["construct"]):
        print(
            f"  {r['construct']:7} {r['orf']:8} idx {r['pool_index']:>3}  {r['barcode']}  ({r['_table']})"
        )
    if not owed:
        print("\nnothing owed — the registry already holds every barcode in the batch")

    if not args.write:
        print("\n(dry run) pass --write to rewrite barcode_reservations.tsv")
        return 0

    out = HERE / "barcode_reservations.tsv"
    cols = [k for k in rows[0] if not k.startswith("_")] if rows else []
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in sorted(owed, key=lambda x: x["construct"]):
            w.writerow({k: r[k] for k in cols})
    print(f"\nWrote {out.name} with {len(owed)} row(s) — append_reservations.py publishes them")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
