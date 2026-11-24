"""Replacing a panel pick the library cannot use, without disturbing any other pick.

Both panels are seeded, class-stratified draws. A pick can turn out unusable only after the
draw: no transcript encodes it, or its designed protein duplicates another member's. The
rule for that case is fixed here so a replacement is reproducible rather than chosen:

1. Draw exactly as the panel was first drawn. Every valid pick is therefore unchanged.
2. Then, in draw order, replace each rejected pick by continuing the **same seeded RNG**
   over its own stratum (tier/band and class): ``rng.choice`` over the stratum's
   candidates that were not picked, not rejected and not invalid, in input-file order.

Rejections come from a reviewed file (``panel_rejections.tsv`` in s4e9, one row per form
with its reason) plus a built-in validity test each sampler applies to its own candidates.
A replacement that later fails is added to the file and the panel redrawn; nothing else
moves, because replacements are drawn after every original pick.
"""

from __future__ import annotations

import csv
import random
from collections.abc import Callable
from pathlib import Path


def load_rejections(path: str | Path | None, panel: str) -> dict[str, str]:
    """``{form id: reason}`` for one panel (``isoform`` or ``ortholog``)."""
    if not path:
        return {}
    with open(path) as fh:
        return {r["form_id"]: r["reason"] for r in csv.DictReader(fh, delimiter="\t")
                if r["panel"] == panel}


def replace_rejected(
    rng: random.Random,
    picked: list[dict],
    stratum: list[dict],
    form_id: Callable[[dict], str],
    rejected: dict[str, str],
    invalid: Callable[[dict], str],
    used: set[str],
) -> tuple[list[dict], list[dict]]:
    """Picks with each rejected one replaced, and a log row per replacement.

    ``used`` holds every form id already in the panel, across strata, and is updated, so a
    form is never drawn twice. ``invalid`` returns a reason, or "" for a usable candidate.
    """
    out, log = [], []
    for p in picked:
        why = rejected.get(form_id(p)) or invalid(p)
        if not why:
            out.append(p)
            continue
        candidates = [
            c for c in stratum
            if form_id(c) not in used and form_id(c) not in rejected and not invalid(c)
        ]
        if not candidates:
            log.append({"rejected": form_id(p), "reason": why, "replacement": "",
                        "candidates": 0})
            continue
        new = rng.choice(candidates)
        used.add(form_id(new))
        out.append(new)
        log.append({"rejected": form_id(p), "reason": why, "replacement": form_id(new),
                    "candidates": len(candidates)})
    return out, log
