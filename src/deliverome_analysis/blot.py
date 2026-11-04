"""Reusable helpers for western blot figures.

These are the general pieces of the blot-figure workflow described in
``docs/western-blot-figures.md``, factored out of the S4E12 figure script so the
standard is carried by code rather than only by a checklist. Everything here takes
arrays and numbers; nothing reads an experiment's files, panel definitions or S3 keys.

Four jobs, each of which has a specific way of going wrong:

``band_centroids``
    Locate bands along a lane axis. Segments contiguous runs above a detrended
    threshold instead of splitting at operator-chosen boundaries. Forcing a split
    point between two merged bands pushes both centroids outward, which in S4E12
    manufactured an apparent gel "smile" that was not there.

``fit_lane_centers`` / :class:`LaneComb`
    Turn band centroids into lane centers. A membrane cut into strips does not have
    one x origin -- the strips are physically separated and repositioned on the
    imager -- so a comb carries a per-strip ``offset`` and the same ``pitch``.

``find_rungs`` / ``assign_rungs``
    Recover a saturated ladder. A saturated band carries no intensity contrast but
    blooms laterally in proportion to its dye load, so rungs appear as maxima in the
    *width* of the saturated column. Which rung is which is then settled by searching
    assignments for maximum log-linearity, optionally checked against bands of known
    mass held out of the fit.

``normalise`` / ``style_panel``
    Render. The black point is the membrane background rather than the darkest noise,
    and the axes grid is turned off, because the shared plotting style enables it and
    a blot must not be crossed by rules.
"""

from __future__ import annotations

import itertools
from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
from scipy.ndimage import label, minimum_filter1d, uniform_filter1d

# LI-COR Chameleon Duo pre-stained ladder, kDa, high to low. Supplied as a default for
# `assign_rungs`; pass the series explicitly for any other ladder.
CHAMELEON_DUO_KDA = (260, 160, 125, 90, 70, 50, 40, 30, 25, 15, 10)


# ---------------------------------------------------------------------------
# Lane geometry
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class LaneComb:
    """Evenly spaced lane centers, with a per-strip x offset.

    ``offset`` exists because a membrane cut into strips and reimaged does not have a
    single x origin. Give each strip its own comb with the same ``pitch`` and its own
    ``offset`` rather than applying one map to both.
    """

    origin: float
    pitch: float
    n_lanes: int
    offset: float = 0.0

    def x(self, lane: int) -> float:
        """x center of ``lane`` (1-based)."""
        if not 1 <= lane <= self.n_lanes:
            raise ValueError(f"lane {lane} outside 1..{self.n_lanes}")
        return self.origin + self.pitch * (lane - 1) + self.offset

    def centers(self) -> list[float]:
        return [self.x(i) for i in range(1, self.n_lanes + 1)]

    def shifted(self, offset: float) -> LaneComb:
        """A copy of this comb with a different strip offset."""
        return LaneComb(self.origin, self.pitch, self.n_lanes, offset)


def band_centroids(
    profile: Sequence[float],
    *,
    x0: float = 0.0,
    rel_height: float = 0.30,
    min_width: int = 10,
    baseline_window: int = 121,
    smooth: int = 3,
) -> list[float]:
    """Intensity-weighted centroids of the bands in a 1-D lane profile.

    Each contiguous run above ``rel_height`` of the detrended maximum is one band, and
    its centroid is weighted by the detrended signal. No window boundaries are imposed,
    which is the point: splitting two merged bands at a chosen x biases both centroids
    away from the split and corrupts any comb fitted to them.
    """
    p = uniform_filter1d(np.asarray(profile, dtype=float), max(smooth, 1))
    detrended = np.clip(p - uniform_filter1d(minimum_filter1d(p, baseline_window), 61), 0, None)
    if not np.any(detrended > 0):
        return []
    runs, n = label(detrended > detrended.max() * rel_height)
    out: list[float] = []
    for k in range(1, n + 1):
        idx = np.flatnonzero(runs == k)
        if idx.size < min_width:
            continue
        w = detrended[idx]
        out.append(float(x0 + (w * idx).sum() / w.sum()))
    return out


def fit_lane_centers(
    centroids: Sequence[float], lanes: Sequence[int], n_lanes: int, *, offset: float = 0.0
) -> tuple[LaneComb, np.ndarray]:
    """Least-squares comb through ``centroids`` observed at 1-based ``lanes``.

    Returns the comb and the per-centroid residuals. Inspect the residuals: a comb that
    fits its own anchors to a few px can still be wrong if the lane assignment is wrong,
    so overlay the centers on the image before trusting any label.
    """
    x = np.asarray(centroids, dtype=float)
    i = np.asarray(lanes, dtype=float)
    if x.size != i.size:
        raise ValueError("centroids and lanes must be the same length")
    if x.size < 2:
        raise ValueError("need at least two anchors to fit a comb")
    pitch, intercept = np.polyfit(i, x, 1)
    comb = LaneComb(float(intercept + pitch), float(pitch), n_lanes, offset)
    residuals = x - np.array([comb.x(int(k)) - offset for k in i])
    return comb, residuals


# ---------------------------------------------------------------------------
# Saturated ladder
# ---------------------------------------------------------------------------


def find_rungs(
    image: np.ndarray,
    *,
    x_range: tuple[int, int],
    y_range: tuple[int, int],
    background: float,
    rel_threshold: float = 0.55,
    min_separation: int = 9,
    smooth: int = 3,
) -> list[int]:
    """Rung rows of a saturated ladder, from the width of its bloom.

    A saturated rung carries no usable intensity contrast, but it spreads sideways in
    proportion to how much dye is in it. Thresholding the ladder column and taking local
    maxima of the per-row width therefore recovers rung positions that an intensity
    profile cannot.

    ``background`` is the membrane level away from the ladder; the threshold sits
    ``rel_threshold`` of the way from it to the column maximum.
    """
    y0, y1 = y_range
    column = np.asarray(image, dtype=float)[y0:y1, x_range[0] : x_range[1]]
    threshold = background + rel_threshold * (column.max() - background)
    width = uniform_filter1d((column > threshold).sum(axis=1).astype(float), max(smooth, 1))
    peaks: list[int] = []
    for j in range(1, len(width) - 1):
        lo = max(0, j - min_separation)
        hi = min(len(width), j + min_separation + 1)
        if width[j] == width[lo:hi].max() and width[j] > width[lo:hi].min():
            if not peaks or j - (peaks[-1] - y0) >= min_separation:
                peaks.append(y0 + j)
    return peaks


@dataclass(frozen=True)
class RungStability:
    """How often each candidate rung row survives a sweep of windows and thresholds.

    A single run of :func:`find_rungs` is sensitive to the ladder column width and the
    threshold, so a rung found once is not necessarily a rung. Sweeping both and keeping
    the support fraction separates the rows that are really there from the ones that are
    an artifact of one setting. In S4E12 this is what showed that two hand-picked rungs
    never survived, while the rest held across 35 combinations.
    """

    support: dict[int, float]  # binned row -> fraction of sweeps that found it
    n_sweeps: int

    def support_for(self, y: float, tolerance: int = 6) -> float:
        """Best support among candidate rows within ``tolerance`` px of ``y``."""
        near = [f for row, f in self.support.items() if abs(row - y) <= tolerance]
        return max(near, default=0.0)


def rung_stability(
    image: np.ndarray,
    *,
    x_center: float,
    y_range: tuple[int, int],
    background: float,
    pads: Sequence[int] = (10, 12, 14, 16, 18),
    rel_thresholds: Sequence[float] = (0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70),
    bin_size: int = 4,
    **kwargs,
) -> RungStability:
    """Run :func:`find_rungs` over a sweep of column widths and thresholds.

    Returns how often each binned row was found. Use it to decide which rungs can be
    trusted without a human looking, and which need confirming.
    """
    votes: dict[int, int] = {}
    n = 0
    for pad in pads:
        for threshold in rel_thresholds:
            n += 1
            rows = find_rungs(
                image,
                x_range=(int(x_center - pad), int(x_center + pad + 2)),
                y_range=y_range,
                background=background,
                rel_threshold=threshold,
                **kwargs,
            )
            for row in rows:
                binned = row // bin_size * bin_size
                votes[binned] = votes.get(binned, 0) + 1
    return RungStability({row: v / n for row, v in votes.items()}, n)


def validate_fit(fit: LadderFit, known: Sequence[tuple[float, float, str]]) -> dict[str, float]:
    """Relative error of ``fit`` at bands of independently known mass.

    ``known`` is ``(y, expected_kDa, label)``. These bands must not be part of the fit --
    that is what makes them evidence. A high r2 alone is not: with three rungs a wrong
    assignment can still be near-perfectly log-linear.
    """
    return {label: (fit.mw(y) - expected) / expected for y, expected, label in known}


@dataclass(frozen=True)
class RungAssignment:
    """A mapping of ladder masses to image rows, with its log-linear fit quality."""

    rungs: dict[float, int]  # kDa -> y
    r2: float

    def fit(self) -> LadderFit:
        return fit_ladder(self.rungs)


def assign_rungs(
    positions: Sequence[int],
    series: Sequence[float],
    *,
    validators: Sequence[tuple[float, float]] = (),
    max_relative_error: float = 0.25,
) -> RungAssignment:
    """Decide which detected row is which ladder rung, by fit rather than by eye.

    Searches every way of assigning ``series`` masses (high to low) to an ordered subset
    of ``positions`` (top to bottom) and keeps the most log-linear. ``validators`` are
    ``(y, known_kDa)`` pairs that are NOT part of the fit -- bands whose mass is known
    independently, such as a loading control. An assignment whose calibration misses any
    validator by more than ``max_relative_error`` is rejected, which is what stops a
    merely self-consistent but wrong assignment from winning.
    """
    ys = sorted(int(p) for p in positions)
    masses = sorted((float(m) for m in series), reverse=True)
    if len(ys) < len(masses):
        raise ValueError(f"only {len(ys)} detected rows for {len(masses)} rungs")

    best: RungAssignment | None = None
    for combo in itertools.combinations(ys, len(masses)):
        rungs = dict(zip(masses, combo, strict=True))
        cal = fit_ladder(rungs)
        if any(abs(cal.mw(y) - known) / known > max_relative_error for y, known in validators):
            continue
        if best is None or cal.r2 > best.r2:
            best = RungAssignment(rungs, cal.r2)
    if best is None:
        raise ValueError("no rung assignment satisfied the validators")
    return best


@dataclass(frozen=True)
class LadderFit:
    """log10(mass) linear in image row, fitted to ladder rungs."""

    slope: float
    intercept: float
    r2: float

    def mw(self, y: float) -> float:
        """Mass in kDa at image row ``y``."""
        return float(10 ** (self.slope * y + self.intercept))

    def y(self, mw: float) -> float:
        """Image row for a mass in kDa."""
        return float((np.log10(mw) - self.intercept) / self.slope)


def fit_ladder(rungs: dict[float, int]) -> LadderFit:
    """Fit log10(kDa) against image row over ``{kDa: y}``."""
    if len(rungs) < 2:
        raise ValueError("need at least two rungs to calibrate")
    ys = np.array([float(y) for y in rungs.values()])
    logs = np.log10(np.array([float(m) for m in rungs]))
    slope, intercept = np.polyfit(ys, logs, 1)
    predicted = slope * ys + intercept
    ss_res = float(((logs - predicted) ** 2).sum())
    ss_tot = float(((logs - logs.mean()) ** 2).sum())
    r2 = 1.0 if ss_tot == 0 else 1.0 - ss_res / ss_tot
    return LadderFit(float(slope), float(intercept), r2)


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------


def normalise(
    panel: np.ndarray, *, black_percentile: float = 40.0, white_percentile: float = 99.8
) -> np.ndarray:
    """Fixed linear stretch of a blot panel to black-on-white in [0, 1].

    No gamma and no per-lane adjustment, so the mapping is reproducible and the panel
    stays quantitatively honest. The black point defaults to the 40th percentile because
    most of a blot is membrane: that percentile *is* the background, and mapping it to
    white is what makes a panel read as a blot. A 2nd-percentile black point sits in the
    noise and renders everything mid-grey.
    """
    arr = np.asarray(panel, dtype=float)
    lo, hi = np.percentile(arr, black_percentile), np.percentile(arr, white_percentile)
    if hi <= lo:
        raise ValueError("white percentile must exceed black percentile in intensity")
    return 1.0 - np.clip((arr - lo) / (hi - lo), 0.0, 1.0)


def style_panel(ax, *, line_width: float = 1.0, color: str = "black") -> None:
    """Box a blot panel: visible black spines, no x ticks, no grid.

    The shared plotting style sets ``axes.grid`` on, which draws rules across the
    membrane. This turns the grid off for one axes only and never touches rcParams.
    """
    ax.grid(False)
    ax.set_xticks([])
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(line_width)
        spine.set_color(color)
    ax.tick_params(axis="y", direction="out", length=4, width=line_width, color=color, pad=3)
