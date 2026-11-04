"""Thin, shared flow-cytometry conveniences for analysis code.

This is a small **lab convenience layer, not a cytometry engine**. Anything that is
genuinely engine/editor work — compensation, GatingML transforms, the gating strategy,
per-gate statistics, FlowJo workspace evaluation — lives in the pinned
``agentflow-cytometry`` package (imported as ``agentflow``) and must stay there
(see CLAUDE.md). What lives here is only the handful of primitives that experiment
*figure* scripts kept re-deriving by hand:

* ``read_fcs`` — read raw events out of an FCS file into a ``{channel: np.ndarray}`` dict,
  the exact ``flowio.FlowData`` → ``reshape`` → ``$PnS``/``$PnN`` index block that was
  copied into eight scripts. For full engine processing (compensation, transforms,
  gating) use ``agentflow`` instead; this is the plain raw-event reader.
* ``gmean`` — the geometric mean over positive events (the CytExpert statistic). Agentflow's
  ``engine.summarize`` emits median/mean but no geometric mean, so this stays here until it
  does; six identical copies existed across S3E9 and S4E1.
* ``asinh`` — a plain ``arcsinh(x / cofactor)`` **display** helper for plotting. It matches
  ``flowkit``'s ``AsinhTransform`` numerically but is not a GatingML transform; gating-time
  transforms belong to agentflow.
* ``cytoflex_live_singlets`` — the documented CytoFLEX scatter/singlet gate used as a
  reproducible, committed stand-in for the hand-drawn FlowJo polygons (so an analysis can
  satisfy the ``flowjo-workspace-not-required`` marker by citing this file).
* ``well_from_name`` — pull a plate well (e.g. ``A1``) out of an FCS sample name.

Heavy imports (``flowio``) are done lazily inside the functions that need them, so code
that only wants ``gmean``/``asinh``/``well_from_name`` does not require the optional
``flow`` dependency group, mirroring how ``deliverome_analysis.flow`` defers to agentflow.
"""

from __future__ import annotations

import re
from collections.abc import Callable
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from os import PathLike

# CytoFLEX exports FITC/APC on a ~150 cofactor by convention; callers pass their own.
DEFAULT_ASINH_COFACTOR = 150.0

# "01-Well-A1.fcs", "something_A01", etc. One capital row letter A-H then 1-2 digits.
WELL_REGEX = re.compile(r"(?<![A-Za-z])([A-H])(\d{1,2})(?![0-9])")


def gmean(values) -> float:
    """Geometric mean over the strictly-positive entries, as CytExpert computes it.

    Zero and negative events (below the fluorophore's background) are dropped rather than
    made to blow up the log; an empty positive set returns NaN.
    """
    values = np.asarray(values, dtype=float)
    positive = values[values > 0]
    if positive.size == 0:
        return float("nan")
    return float(np.exp(np.log(positive).mean()))


def asinh(x, cofactor: float = DEFAULT_ASINH_COFACTOR) -> np.ndarray:
    """``arcsinh(x / cofactor)`` — a display transform for plotting, not for gating.

    Numerically equal to ``flowkit.transforms.AsinhTransform`` at the same cofactor, but
    this is a plain helper: gating-time transforms belong to the agentflow engine.
    """
    return np.arcsinh(np.asarray(x, dtype=float) / cofactor)


def well_from_name(name: str, pattern: re.Pattern[str] = WELL_REGEX) -> str | None:
    """The plate well in an FCS sample name (``01-Well-A1.fcs`` → ``A1``), or None.

    The digits are not zero-padded, so ``A01`` and ``A1`` both return ``A1``.
    """
    match = pattern.search(name)
    if not match:
        return None
    return f"{match.group(1)}{int(match.group(2))}"


def read_fcs(
    source: str | PathLike,
    *,
    channels: list[str] | None = None,
    download: Callable[[str], Path] | None = None,
    prefer: str = "stain",
) -> dict[str, np.ndarray]:
    """Raw events from an FCS file as ``{channel_label: 1-D np.ndarray}``.

    This centralises the ``flowio.FlowData`` → ``np.reshape`` → ``$PnS``/``$PnN`` index block
    that was copied across the S3E9/S3E13 figure scripts. It returns **raw** events — no
    compensation or transform; for those use ``agentflow``.

    ``source`` is a local path, or an S3 key when ``download`` (e.g. ``s3.external.download``)
    is given to fetch it first. ``prefer`` picks the channel label: ``"stain"`` uses ``$PnS``
    (the marker name) and falls back to ``$PnN`` (the detector); ``"detector"`` uses ``$PnN``
    only — the MA900 export carries no ``$PnS``. ``channels`` restricts the returned dict.
    """
    import flowio

    path = Path(download(str(source))) if download is not None else Path(source)
    fcs = flowio.FlowData(str(path))
    matrix = np.reshape(fcs.events, (-1, fcs.channel_count))

    labels: list[str] = []
    for i in range(1, fcs.channel_count + 1):
        detector = fcs.text.get(f"p{i}n")
        stain = fcs.text.get(f"p{i}s")
        labels.append((stain or detector) if prefer == "stain" else detector)

    out = {label: matrix[:, i] for i, label in enumerate(labels) if label}
    if channels is not None:
        missing = [c for c in channels if c not in out]
        if missing:
            raise KeyError(f"{path.name} has no channel(s) {missing}; present: {sorted(out)}")
        out = {c: out[c] for c in channels}
    return out


def cytoflex_live_singlets(events: dict[str, np.ndarray]) -> np.ndarray:
    """Boolean mask reproducing the S3E9 ``all/Single Cells/live`` CytoFLEX hierarchy.

    The published hierarchy is hand-drawn in FlowJo and its polygons are not recoverable
    from the exported statistics, so this is a committed, reproducible stand-in — not a copy
    of the polygons. It is the gate recorded in code that lets an analysis cite this file in
    a ``flowjo-workspace-not-required`` marker. Expects raw ``FSC-A``/``FSC-H``/``SSC-A``.
    """
    fsc_a, fsc_h, ssc_a = events["FSC-A"], events["FSC-H"], events["SSC-A"]
    cells = (fsc_a > 2.2e5) & (fsc_a < 3e6) & (ssc_a > 3e4) & (ssc_a < 5e6)
    ratio = np.divide(fsc_h, fsc_a, out=np.zeros_like(fsc_h, dtype=float), where=fsc_a > 0)
    singlets = cells & (ratio > 0.40) & (ratio < 0.75)
    # Drop the low-FSC / high-SSC corner the exported "live" gate excludes.
    return singlets & (ssc_a < 10 ** (0.95 * np.log10(np.clip(fsc_a, 1, None)) + 1.05))
