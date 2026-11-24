"""Tedman's admissible interval for a tag behind a signal peptide, reimplemented here.

The rule of record lives in the accessible-surfaceome repository
(``src/accessible_surfaceome/tag_sites/sp_interval.py``; read at ``dev`` 8547d572). This repo
must not import that private package, so the four pure functions the L1 interval census used
are restated below, unchanged in behaviour. The census itself ran the sibling module at
accessible-surfaceome commit ``1c89c8a19`` (branch ``claude/l1-post-sp-placement``, not pushed),
so equality with the census record, checked by ``nterm_placement.py`` for every member, is the
test that this copy and that one agree.

The rule (Tedman, Deep Receptor Scanning methods): call the cleavage site, find where
structure begins, and insert **between** the two.

* ``lo`` is the cleavage site: the **later** of the DeepTMHMM v1 and SignalP 6.0 calls,
  whichever answered (a cassette inside the signal peptide leaves with it, so err late).
* Each disorder track (higher = more disordered) gives an **onset**: the first residue after
  ``lo`` that begins five consecutive residues scored below 0.5.
* The interval is ``[lo, median_onset - 1 + 2]``: the two-residue tolerance covers six of the
  eight Tedman junctions that overshoot the strict interval. With no onset at all, the onset is
  taken as ``lo + fallback_width + 1``, so with the tolerance the interval ends at
  ``lo + fallback_width + tolerance``.
* Candidate insertion sites are the interval from ``lo + 1``, capped at ``lo + 20``. ``lo``
  itself is excluded: a cassette there forms the mature N-terminus, which none of Tedman's
  signal-peptide-bearing constructs does.

All positions are 1-based "insert after residue n" on the designed protein.
"""

from __future__ import annotations

from typing import NamedTuple

STRUCTURE_THRESHOLD = 0.5
STRUCTURE_RUN = 5
ONSET_TOLERANCE = 2
#: The census passed ``fallback_width=6`` (its docstring); the sibling module's default is 7.
FALLBACK_WIDTH = 6
CANDIDATE_SPAN = 20


class TagInterval(NamedTuple):
    lo: int
    hi: int
    onsets: dict[str, int | None]
    cleavage_site: int

    @property
    def width(self) -> int:
        return max(0, self.hi - self.lo)


def structure_onset(
    track: list[float],
    start: int,
    *,
    threshold: float = STRUCTURE_THRESHOLD,
    run: int = STRUCTURE_RUN,
) -> int | None:
    """First residue after ``start`` (1-based) opening ``run`` residues below ``threshold``."""
    if start < 0 or not track:
        return None
    for i in range(start, len(track) - run + 1):
        if all(track[j] < threshold for j in range(i, i + run)):
            return i + 1
    return None


def tag_interval(
    cleavage_site: int,
    tracks: dict[str, list[float]],
    *,
    tolerance: int = ONSET_TOLERANCE,
    fallback_width: int = FALLBACK_WIDTH,
) -> TagInterval:
    """``[cleavage_site, median onset - 1 + tolerance]``; the median is the upper middle."""
    onsets = {name: structure_onset(t, cleavage_site) for name, t in tracks.items()}
    found = sorted(v for v in onsets.values() if v is not None)
    onset = found[len(found) // 2] if found else cleavage_site + fallback_width + 1
    return TagInterval(
        lo=cleavage_site,
        hi=max(cleavage_site, onset - 1 + tolerance),
        onsets=onsets,
        cleavage_site=cleavage_site,
    )


def cleavage_lo(deeptmhmm_v1: int | None, signalp6: int | None) -> int | None:
    """The later of the two cleavage calls, using whichever answered."""
    answered = [c for c in (deeptmhmm_v1, signalp6) if c is not None]
    return max(answered) if answered else None


def candidate_positions(
    interval: TagInterval, cleavage: int, span: int = CANDIDATE_SPAN
) -> list[int]:
    """Insert-after positions inside the interval, past ``cleavage``, at most ``span`` past it."""
    return list(range(cleavage + 1, min(interval.hi, cleavage + span) + 1))
