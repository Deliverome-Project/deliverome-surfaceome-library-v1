"""Allocate above both manifest history and published adoption reservations."""

from __future__ import annotations

import csv
import re
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path

HERE = Path(__file__).resolve().parent
#: Every published allocation ledger, oldest first. compara116's reserved l1-p2743..2931
#: for a panel that was never activated; ortholog_cds_policy's are live members (2026-11-10).
COMPARA116_RESERVATIONS = HERE / "compara116_id_allocations.tsv"
ORTHOLOG_CDS_POLICY_ALLOCATIONS = HERE / "ortholog_cds_policy_id_allocations.tsv"
DEFAULT_RESERVATIONS = (COMPARA116_RESERVATIONS, ORTHOLOG_CDS_POLICY_ALLOCATIONS)


def next_l1_member_number(
    rows: Iterable[Mapping[str, str]], reservation_paths: Sequence[Path] | None = None
) -> int:
    """Return the next unused L1 number, including retired and reserved IDs.

    Reservations are authoritative even before their adoption manifest becomes the
    default. A missing ledger is allowed for historical checkouts; a present but
    malformed ledger fails closed instead of making published IDs reusable.
    """
    numbers = []

    def consider(cid: str, *, reservation: bool = False):
        match = re.fullmatch(r"l1-p([0-9]+)", cid)
        if match:
            numbers.append(int(match[1]))
        elif reservation or cid.startswith("l1-"):
            raise ValueError(f"invalid reserved L1 candidate_id: {cid!r}")

    for row in rows:
        consider(row["candidate_id"])
    for path in DEFAULT_RESERVATIONS if reservation_paths is None else reservation_paths:
        path = Path(path)
        if not path.exists():
            continue
        with path.open() as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if "candidate_id" not in (reader.fieldnames or []):
                raise ValueError(f"reservation ledger has no candidate_id column: {path}")
            for row in reader:
                consider(row["candidate_id"], reservation=True)
    return max(numbers, default=0) + 1
