"""Allocate design ids above manifest history and published adoption reservations.

Design-stage scripts (adding genes, redrawing panels, adopting a panel, retiring duplicates)
allocate and refer to **design ids** (``l1-pNNNN``), never to member ids (``l1-mNNNN``),
which ``renumber_library_members.py`` reassigns whenever membership changes. They read the
manifest through :func:`to_design_space`, so ``candidate_id`` holds the design id while
they work; the manifest they write is unnumbered until the renumber stage runs again.
"""

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


RENUMBER_HINT = (
    "membership changed: rerun renumber_library_members.py --write to renumber the L1 "
    "members and re-key their records"
)


def to_design_space(rows: list[dict]) -> list[dict]:
    """Rows with ``candidate_id`` set to the design id, for design-stage code.

    A manifest written before the renumbering has no ``design_id`` column and is returned
    unchanged, so replays of earlier stages from their pinned manifests are byte-stable.
    """
    if not any("design_id" in r for r in rows):
        return rows
    out = []
    for r in rows:
        did = (r.get("design_id") or "").strip() or r["candidate_id"]
        if not re.fullmatch(r"l[12]-p[0-9]{4}", did):
            raise ValueError(f"row has no design id: {r['candidate_id']!r}")
        out.append({**r, "design_id": did, "candidate_id": did})
    return out


def next_l1_member_number(
    rows: Iterable[Mapping[str, str]], reservation_paths: Sequence[Path] | None = None
) -> int:
    """Return the next unused L1 design-id number, including retired and reserved IDs.

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
        # Design ids only: a member id (l1-mNNNN) is not an allocation.
        consider((row.get("design_id") or "").strip() or row["candidate_id"])
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
