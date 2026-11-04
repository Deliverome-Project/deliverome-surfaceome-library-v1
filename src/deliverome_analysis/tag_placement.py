"""Choosing where a tag cassette goes at a signal-peptide-bearing N-terminus.

The obvious rule -- splice at the signal-peptide cleavage site, so the mature protein
begins with the cassette -- is not what works. Of Tedman's 108 screen-validated HA
insertions whose accession carries a DeepTMHMM v1 signal peptide, **zero** sit at the
cleavage site: 107 are downstream, median +4 residues, range -1 to +19. Those junctions
land at roughly the 90th percentile of their own receptor's disorder, and our own 25
curated internal positive controls agree at the 73rd-83rd.

So the rule is: search a short window past the cleavage site and take the most flexible
residue in it.

Anchoring on flexibility rather than on the cleavage site also discharges a dependency we
cannot otherwise settle. DeepTMHMM v1 and SignalP 6 agree on the cleavage site for only 45
of 65 accessions in that set; the other 20 differ by 1-4 residues, which is the same scale
as the offset itself. "Cleavage site + 4" inherits that error -- a tag placed 4 past a site
called 4 early lands back at the cleavage site. "Most disordered residue in the window"
does not, because disorder is a property of the sequence: the same residue wins whichever
call sets the window's start.

Disorder tracks are passed in rather than fetched. They live in the surfaceome project's
``disorder_public`` (metapredict, NetSurfP-3.0, AlphaFold-disorder) and this package has no
business reaching into another project's database; keeping the scoring pure also makes it
testable without a network.
"""

from __future__ import annotations

from typing import NamedTuple

DEFAULT_WINDOW = 20
"""How far past the cleavage site to search, in residues.

Covers the full range Tedman's constructs occupy (max +19). Narrower would force some
genes onto a less flexible residue; wider leaves the region the signal peptide's own
geometry governs.
"""


class JunctionChoice(NamedTuple):
    """Where the cassette goes, and the evidence for it.

    ``junction`` follows the "after N" convention used throughout this project: the
    cassette is inserted between residues ``junction`` and ``junction + 1``, 1-indexed.
    """

    junction: int
    sp_end_used: int
    mean_percentile: float
    per_track_percentile: tuple[float, ...]
    window: tuple[int, int]

    def describe(self) -> str:
        lo, hi = self.window
        return (
            f"after residue {self.junction} (+{self.junction - self.sp_end_used} past the "
            f"cleavage site at {self.sp_end_used}, searched {lo}..{hi}), "
            f"{self.mean_percentile:.0f}th disorder percentile"
        )


def within_protein_percentile(track: list[float | None]) -> list[float]:
    """Each position's rank in its own protein, 0-100.

    Predictors disagree about absolute disorder far more than about rank -- metapredict
    puts Tedman's junctions at 0.61 and NetSurfP-3.0 at 0.93 for the same residues -- so
    scoring on rank is what lets three tracks be averaged without one dominating.
    ``None`` entries are treated as absent and score 0.
    """
    vals = [v for v in track if v is not None]
    if not vals:
        return [0.0] * len(track)
    n = len(vals)
    ordered = sorted(vals)
    out = []
    for v in track:
        if v is None:
            out.append(0.0)
            continue
        lo, hi = 0, n
        while lo < hi:  # bisect_left without the import
            mid = (lo + hi) // 2
            if ordered[mid] < v:
                lo = mid + 1
            else:
                hi = mid
        out.append(100.0 * lo / n)
    return out


def post_sp_junction(
    sequence: str,
    sp_end_calls: list[int],
    disorder_tracks: list[list[float | None]],
    *,
    window: int = DEFAULT_WINDOW,
) -> JunctionChoice:
    """Pick the cassette junction downstream of a signal peptide.

    ``sp_end_calls`` are each predictor's last signal-peptide residue; the **latest** is
    used. The failure modes are not symmetric: a cassette placed inside the signal peptide
    leaves with it at cleavage and the member reads as a non-expressor, indistinguishable
    from a real biological negative, whereas one placed a few residues late is still
    presented. Erring late is therefore the safe direction.

    ``disorder_tracks`` are per-residue, higher meaning more disordered. Pass
    ``1 - pLDDT`` for AlphaFold-disorder. A track may be **shorter** than the protein and
    is used as long as it reaches the end of the search window: NetSurfP-3.0 scores only
    the first 400 residues by design, and since the window sits at the N-terminus that
    truncation costs nothing here. A track that stops before the window ends is dropped
    rather than extrapolated. With no usable tracks the
    search cannot be made and the first mature residue is returned, which reproduces the
    old cleavage-site behaviour -- the caller should treat that as a degraded result, not a
    choice.
    """
    calls = [c for c in sp_end_calls if c and c > 0]
    if not calls:
        raise ValueError("no signal-peptide call: this is not a post-SP placement")
    sp_end = max(calls)
    if sp_end >= len(sequence):
        raise ValueError(f"signal peptide ends at {sp_end} but the protein is {len(sequence)} aa")

    lo = sp_end + 1
    hi = min(sp_end + window, len(sequence) - 1)
    if hi < lo:
        hi = lo
    usable = [t for t in disorder_tracks if len(t) >= hi]
    if not usable:
        return JunctionChoice(lo, sp_end, 0.0, (), (lo, hi))

    pct = [within_protein_percentile(t) for t in usable]
    best, best_score, best_each = lo, -1.0, ()
    for res in range(lo, hi + 1):
        each = tuple(p[res - 1] for p in pct)
        score = sum(each) / len(each)
        if score > best_score:
            best, best_score, best_each = res, score, each
    return JunctionChoice(best, sp_end, best_score, best_each, (lo, hi))
