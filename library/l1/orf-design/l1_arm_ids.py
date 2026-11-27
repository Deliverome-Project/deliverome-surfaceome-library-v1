"""Which arm of which member each built insert is: the arm id (``l1-pNNNN-aK``).

From the first per-insert cut an L1 member may be ordered as several inserts: an internal
member gets 1 to 3 loop arms, and a snorkel only beside a loop arm (Becca, 2026-11-14/15). The
order, the plate map, the synthesis FASTA, the fold check and swap risk therefore key on the
arm, not the member. The arm id is permanent: K is assigned once, in creation order, by the
session that designs the arm, and never reused or re-ranked (schema agreed 2026-11-14/15 by
the coordinator, l1-integration, ranker-pu-segment and l1-reviewer-package; see
``deliverome_analysis.cloning.library_design``).

Where each arm id comes from:

* an arm table row carries it: the ``arm_id`` column (snorkel rows of ``constructs_l1_cterm``)
  or a ``construct_id`` that is an arm id (``constructs_l1_internal``, one row per loop arm);
* an internal-rung arm (``constructs_l1_internal``, construct ``l1-pNNNN:internalN``) and a
  snorkel arm (``constructs_l1_cterm`` row whose placement is a snorkel, key
  ``l1-pNNNN:snorkel``) take it from the ranker's per-arm table
  (S5E13 ``l1_internal_arm_decisions.tsv``, column ``construct_key``);
* every other rung builds one arm per member, so its insert is the member's ``-a1``.

Nothing falls back: an internal or snorkel insert with no row in the per-arm table, a row
whose decision is not ``build``, or two inserts with one arm id all raise.

Implemented by the Claude Code session l1-reviewer-package (2026-11-15) at the coordinator's
request ("prepare the order side for one row per insert").
"""

from __future__ import annotations

import csv
from collections.abc import Iterable, Mapping
from pathlib import Path

from deliverome_analysis.cloning import library_design as ld

HERE = Path(__file__).resolve().parent
ARM_DECISIONS = HERE.parents[1] / "l1/tag-site-benchmarking" / "l1_internal_arm_decisions.tsv"
#: Rungs whose arm table holds exactly one arm per member.
SINGLE_ARM_RUNGS = ("nterm", "multipass", "type3", "cterm")


def is_snorkel(row: Mapping[str, str]) -> bool:
    return (row.get("placement") or "").startswith("snorkel")


def construct_key(rung: str, row: Mapping[str, str], design: str) -> str | None:
    """The per-arm table's ``construct_key`` for an internal or snorkel arm, else None."""
    if rung == "internal":
        return row["construct_id"]
    if rung == "cterm" and is_snorkel(row):
        return f"{design}:snorkel"
    return None


def load_decisions(path: Path = ARM_DECISIONS) -> dict[str, dict]:
    """``construct_key`` -> per-arm row, for every arm the ranker decided."""
    with path.open() as fh:
        rows = [r for r in csv.DictReader(fh, delimiter="\t") if r["construct_key"]]
    out: dict[str, dict] = {}
    for r in rows:
        if not ld.ARM_ID_RE.fullmatch(r["arm_id"]):
            raise ValueError(f"{path.name}: not an arm id: {r['arm_id']!r}")
        if ld.design_of_arm(r["arm_id"]) != r["design_id"]:
            raise ValueError(f"{path.name}: {r['arm_id']} is not an arm of {r['design_id']}")
        if r["construct_key"] in out:
            raise ValueError(f"{path.name}: construct_key {r['construct_key']} twice")
        out[r["construct_key"]] = r
    return out


def arm_id_for(rung: str, row: Mapping[str, str], design: str, decisions: Mapping) -> str:
    """The arm id of one built insert (see the module docstring for the rules)."""
    if row.get("arm_id"):
        arm = row["arm_id"]
    elif ld.ARM_ID_RE.fullmatch(row.get("construct_id") or ""):
        # The internal builder names each loop arm by its arm id (construct_id = arm_id).
        arm = row["construct_id"]
    else:
        key = construct_key(rung, row, design)
        if key is None:
            if rung not in SINGLE_ARM_RUNGS:
                raise ValueError(f"{rung}: unknown rung")
            return ld.arm_id(design, 1)
        dec = decisions.get(key)
        if dec is None:
            raise ValueError(f"{key}: built, but the per-arm table has no row for it")
        if dec["decision"] != "build":
            raise ValueError(f"{key}: built, but the per-arm table decides {dec['decision']!r}")
        arm = dec["arm_id"]
    if ld.design_of_arm(arm) != design:
        raise ValueError(f"{arm}: not an arm of {design}")
    return arm


def check_unique(arms: Iterable[str]) -> None:
    seen: set[str] = set()
    dup = sorted({a for a in arms if a in seen or seen.add(a)})
    if dup:
        raise ValueError(f"arm ids on more than one insert: {dup[:5]}")
