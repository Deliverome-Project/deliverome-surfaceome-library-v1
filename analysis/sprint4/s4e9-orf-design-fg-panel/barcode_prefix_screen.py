"""Keep barcodes distinct over the bases in situ sequencing reads.

Optical pooled screening reads a barcode by in situ sequencing by synthesis, which reaches
only the first 9-12 bases (Feldman et al., Cell 2019). The s1e2 pool's barcodes are
distinct over all 20 nt but not over short prefixes: of its 28,506 barcodes, 1,244 pairs
share their first 10 bases, and any random draw of 2,735 contains such a pair. A draw
screened with ``PrefixScreen`` keeps every pair at least ``min_distance`` substitutions
apart over the first ``bases`` bases, so a member can still be told from every other
member by the bases in situ sequencing reads.

The default, minimum Hamming distance 3 over the first 10 bases, corrects one in situ read
error. The pool satisfies it for 2,735 barcodes (``tests/test_barcode_prefix_screen.py``).
The readout is assumed to start at the barcode's 5' end, next to the forward handle, as
Feldman's sequencing primer does.

How it works: two prefixes within distance ``min_distance - 1`` agree on at least
``bases - (min_distance - 1)`` positions, so they share a key made by masking some
``min_distance - 1`` positions. Each accepted barcode registers all its masked keys; a
candidate is rejected when any of its keys is already registered.
"""

from __future__ import annotations

import itertools
from collections.abc import Iterable

#: Bases read by in situ sequencing, and the minimum pairwise Hamming distance over them.
PREFIX_BASES = 10
PREFIX_MIN_DISTANCE = 3


class PrefixScreen:
    """Accept barcodes only if every pair stays ``min_distance`` apart over ``bases``."""

    def __init__(
        self,
        bases: int = PREFIX_BASES,
        min_distance: int = PREFIX_MIN_DISTANCE,
        existing: Iterable[str] = (),
    ) -> None:
        if not 1 <= min_distance <= bases:
            raise ValueError(f"min_distance {min_distance} must lie in 1..{bases}")
        self.bases = bases
        self.min_distance = min_distance
        self._masks = list(itertools.combinations(range(bases), min_distance - 1))
        self._keys: set[tuple[tuple[int, ...], str]] = set()
        # Barcodes already in use are registered without checking them against each
        # other: they are built, and a draw can only avoid them.
        for bc in existing:
            self._keys.update(self._keys_of(bc))

    def _keys_of(self, bc: str) -> list[tuple[tuple[int, ...], str]]:
        prefix = bc.upper()[: self.bases]
        if len(prefix) < self.bases:
            raise ValueError(f"barcode {bc!r} is shorter than {self.bases} bases")
        return [
            (mask, "".join("." if i in mask else b for i, b in enumerate(prefix)))
            for mask in self._masks
        ]

    def passes(self, bc: str) -> bool:
        """True if ``bc`` keeps the distance to every barcode accepted so far."""
        return not any(k in self._keys for k in self._keys_of(bc))

    def accept(self, bc: str) -> None:
        self._keys.update(self._keys_of(bc))

    def try_accept(self, bc: str) -> bool:
        if not self.passes(bc):
            return False
        self.accept(bc)
        return True


def min_prefix_distance(barcodes: Iterable[str], bases: int = PREFIX_BASES) -> int:
    """Smallest pairwise Hamming distance over the first ``bases`` bases (brute force)."""
    prefixes = [b.upper()[:bases] for b in barcodes]
    best = bases
    for a, b in itertools.combinations(prefixes, 2):
        d = sum(x != y for x, y in zip(a, b, strict=True))
        if d < best:
            best = d
            if best == 0:
                break
    return best
