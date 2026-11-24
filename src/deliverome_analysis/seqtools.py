"""Sequence-design helpers for Type IIS / Golden Gate construct assembly.

Gene symbol -> Ensembl stable ID -> codon-optimized or minimally-domesticated
CDS, plus a constraint-based DNA barcode generator.

A shared team utility with three independent features:

**1. Codon recoding** (:func:`recode_gene`, CLI ``recode``). Given a human
gene symbol it:

1. Resolves the symbol to an Ensembl **gene stable ID** (ENSG) on the latest
   GRCh38 assembly via the Ensembl REST API, and picks the **canonical**
   transcript (which matches MANE Select for the large majority of genes).
2. Downloads that transcript's coding sequence (CDS).
3. Runs a single `DNA Chisel <https://github.com/Edinburgh-Genome-Foundry/DnaChisel>`_
   optimization that jointly (a) codon-optimizes to a *Homo sapiens* usage
   profile, (b) removes every chosen restriction-enzyme recognition site
   (both strands), (c) keeps %GC within a windowed target, and (d) avoids
   homopolymer runs and repeated k-mers — all while an ``EnforceTranslation``
   constraint keeps the protein byte-for-byte identical.

**2. Minimal-edit domestication** (:func:`domesticate_gene` / :func:`domesticate`,
CLI ``domesticate``). For when you want the *native* codon usage left alone —
this only removes restriction sites, changing as few nucleotides as possible
(via DNA Chisel's ``AvoidChanges`` objective), with no codon optimization, no
GC windowing, and no repeat handling. ``domesticate()`` also works on
non-coding sequences (an IRES, a UTR) via ``preserve_translation=False``,
since there's no protein to keep constant there at all.

**3. Barcode generation** (:func:`generate_barcode`, CLI ``barcode``). Draws a
random DNA barcode of a chosen length that satisfies a %GC window, a
homopolymer-run limit, and freedom from chosen restriction-enzyme sites (plus
any other motif you want to avoid, e.g. the polyA signal ``AATAAA``) — for
tagging individual constructs (e.g. one barcode per receptor/ORF in a pooled
surface-protein screen) so each can be identified later by gDNA/amplicon
sequencing.

Design notes
------------
* Stable IDs, not symbols, are the durable identifier for recoding. The
  result records the ENSG (versioned + unversioned), the canonical ENST, and
  the ENSP so the provenance survives a symbol rename.
* Recoding and domestication both depend on ``dnachisel`` (+ its
  ``numpy``/``biopython`` transitive deps, already used elsewhere in this
  repo). An earlier version of this module did its own "most-frequent-codon +
  iterative site-scrub" pass for recoding; that greedy approach reliably
  clears sites but, for repeat-heavy stretches (e.g. Gly/Ser linkers), tends
  to pick the *same* synonymous codon for every occurrence of a residue — a
  monotonous, repeat-prone sequence that's a real synthesis/recombination
  liability. DNA Chisel solves codon usage, site avoidance, GC windowing, and
  repeat avoidance **jointly** instead of sequentially, which avoids that
  failure mode.
* ``domesticate()`` is the same joint-solving idea applied to a different
  goal: instead of optimizing codon usage, its only objective is minimizing
  total edits from the input. Useful whenever "codon-optimize this" isn't
  what you want — a native sequence you want to keep recognizable, or a
  structural RNA element (IRES, UTR) where codon optimization doesn't apply
  at all and blind editing risks disrupting secondary structure. It still
  cannot verify that RNA secondary structure survives an edit — only that
  translation is unchanged when ``preserve_translation=True``.
* Barcode generation is rejection sampling: draw a random sequence, check
  every constraint, retry on failure. All constraints here are independently
  satisfiable, so this converges in a handful of draws for typical
  parameters (e.g. a 20-nt barcode at 40-60% GC). This part remains
  stdlib-only (``random``) — no DNA Chisel dependency.

CLI
---
    uv run python -m deliverome_analysis.seqtools recode TF
    uv run python -m deliverome_analysis.seqtools recode TF --output tf.fasta --json tf.json
    uv run python -m deliverome_analysis.seqtools domesticate TF --enzymes BbsI,Esp3I,BsaI
    uv run python -m deliverome_analysis.seqtools barcode --length 20 --label TFRC
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
import os
import random
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from collections.abc import Iterable, Mapping
from dataclasses import asdict, dataclass, field
from typing import Any, NamedTuple, cast

import dnachisel as dc
import numpy as np

# --------------------------------------------------------------------------- #
# Genetic code + Homo sapiens codon usage
# --------------------------------------------------------------------------- #

GENETIC_CODE: dict[str, str] = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "TGT": "C",
    "TGC": "C",
    "TGA": "*",
    "TGG": "W",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}

# Type IIS enzymes whose recognition sites we scrub. Values are the 5'->3'
# recognition sequence; both strands are checked via the reverse complement.
ENZYME_SITES: dict[str, str] = {
    "BbsI": "GAAGAC",  # BpiI isoschizomer
    "Esp3I": "CGTCTC",  # BsmBI isoschizomer
    "BsaI": "GGTCTC",
}

_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def reverse_complement(seq: str) -> str:
    """Reverse complement of an uppercase DNA string."""
    return seq.translate(_COMPLEMENT)[::-1]


# Biopython's GenBank writer indents annotation continuation lines by 12 and wraps at
# 80, so a whitespace-free token longer than 68 cannot be wrapped. Biopython emits a
# "too long" BiopythonWarning and then writes it anyway, overflowing the line: the text
# survives a Biopython round trip, but the file is no longer within the 80-column GenBank
# convention, which stricter third-party parsers can reject. Repo-relative script paths
# in provenance notes are the usual offender.
GENBANK_ANNOTATION_WIDTH = 68


def overlong_genbank_words(text: str, width: int = GENBANK_ANNOTATION_WIDTH) -> list[str]:
    """Words in `text` too long for the GenBank writer to wrap within 80 columns.

    Call this on a provenance note before assigning it to ``record.annotations``.
    An empty list means every line of the written record stays within the GenBank
    column convention. Biopython itself still writes an over-long word, so this
    guards file portability, not loss of the text.
    """
    return [w for w in text.split() if len(w) > width]


# Template switching during lentiviral reverse transcription reassorts the two
# co-packaged genomes. Where two library members are identical, a switch is seamless, so
# any member-specific element separated from another by a shared tract can swap partners.
# Walton et al. measured that rate against separation and fit
#     Y(x) = 1/2 (1 - e^(-2kx)),  k = 5.3e-4 bp^-1
# across designs from a 19 bp spacer-iBAR pair up to multi-kb distal barcodes, where it
# approaches the 50% asymptote (half the genomes are the partner's).
# CROPseq-multi: a versatile solution for multiplexed perturbation and decoding in pooled
# CRISPR screens, bioRxiv 2024.03.17.585235.
WALTON_SWAP_K = 5.3e-4


def barcode_swap_rate(distance_bp: int, k: float = WALTON_SWAP_K) -> float:
    """Predicted fraction of barcodes that end up on the wrong member.

    ``distance_bp`` is the shared tract separating the two member-specific elements —
    for a barcoded ORF library, the sequence between where the ORFs stop differing and
    where the barcode begins. That is what a switch has to land in to go unnoticed.
    """
    if distance_bp < 0:
        raise ValueError("distance must not be negative")
    return 0.5 * (1 - math.exp(-2 * k * distance_bp))


class SwapRisk(NamedTuple):
    """Worst-case barcode/ORF reassortment across a library."""

    distance_bp: int
    rate: float
    members: tuple[str, str]

    def describe(self) -> str:
        a, b = self.members
        return (
            f"{self.distance_bp} bp shared between the ORF and the barcode "
            f"({a} vs {b}) -> {self.rate:.1%} predicted swapping"
        )


def longest_shared_tract(a: str, b: str) -> tuple[int, int, int]:
    """Longest run of identical bases anywhere in two sequences: (length, i, j).

    Continuous identity is what a template switch needs: a switch inside a shared run
    is seamless, and one outside it leaves a mismatch. So the longest continuous tract
    is the library's recombination substrate, wherever it falls -- a conserved region
    between two paralogues counts as much as a shared vector cassette.

    Binary search over the length with a set of substrings, not the O(n*m) dynamic
    program: at library scale the quadratic version is the difference between seconds
    and hours, and thousands of pairs have to be screened to find the few that matter.
    """
    lo, hi = 0, min(len(a), len(b))
    best = (0, 0, 0)
    while lo < hi:
        mid = (lo + hi + 1) // 2
        seen: dict[str, int] = {}
        for i in range(len(a) - mid + 1):
            seen.setdefault(a[i : i + mid], i)
        hit = None
        for j in range(len(b) - mid + 1):
            at = seen.get(b[j : j + mid])
            if at is not None:
                hit = (mid, at, j)
                break
        if hit:
            best, lo = hit, mid
        else:
            hi = mid - 1
    return best


def barcode_swap_risk(members: dict[str, tuple[str, str]], *, k: float = WALTON_SWAP_K) -> SwapRisk:
    """Worst-case swap risk over every pair in ``members``.

    ``members`` maps a name to ``(construct_sequence, barcode)``. For each pair the
    shared tract is the run of identical bases ending where each member's barcode
    begins; the pair with the longest such run sets the library's risk, because one
    recombinogenic pair is enough to mislabel cells.
    """
    if len(members) < 2:
        raise ValueError("need at least two members to have a pair to swap between")
    worst = (-1, ("", ""))
    names = sorted(members)
    for i, x in enumerate(names):
        sx, bx = members[x]
        ix = sx.find(bx)
        if ix < 0:
            raise ValueError(f"{x}: barcode {bx} is not in its construct")
        for y in names[i + 1 :]:
            sy, by = members[y]
            iy = sy.find(by)
            if iy < 0:
                raise ValueError(f"{y}: barcode {by} is not in its construct")
            n = 0
            while n < min(ix, iy) and sx[ix - 1 - n] == sy[iy - 1 - n]:
                n += 1
            if n > worst[0]:
                worst = (n, (x, y))
    return SwapRisk(worst[0], barcode_swap_rate(worst[0], k), worst[1])


_TRACT_BASE = np.uint64(0x9E3779B97F4A7C15)  # odd 64-bit multiplier; arithmetic wraps mod 2**64
_TRACT_CODE: dict[str, int] = {"A": 1, "C": 2, "G": 3, "T": 4}
_TRACT_LADDER_STEP = 16


class TractIndex:
    """Polynomial hashes of every substring of every DNA window, vectorized over the pool.

    With ``y[k] = x[k] * B**k`` and ``S`` its running sum, ``S[p+L] - S[p]`` is a substring's
    hash times ``B**p``; multiplying by ``B**(N-p)`` puts every substring on the same footing,
    so equal substrings get equal values wherever they sit. All arithmetic wraps mod 2**64.
    """

    def __init__(self, windows: list[str], owners: list[int]) -> None:
        self.windows = windows
        lens = np.array([len(w) for w in windows], dtype=np.int64)
        self.lens = lens
        self.starts = np.concatenate([[0], np.cumsum(lens + 1)[:-1]]).astype(np.int64)
        text = "\x00".join(windows) + "\x00"
        codes = np.frombuffer(text.encode(), dtype=np.uint8)
        x = np.zeros(len(codes), dtype=np.uint64)
        for ch, v in _TRACT_CODE.items():
            x[codes == ord(ch)] = v
        n = len(codes)
        with np.errstate(over="ignore"):
            self.pw = np.cumprod(np.full(n + 1, _TRACT_BASE, dtype=np.uint64))
            self.pw = np.concatenate([[np.uint64(1)], self.pw[:-1]])
            self.cum = np.concatenate([[np.uint64(0)], np.cumsum(x * self.pw[:n])])
        self.n = n
        self.window_at = np.repeat(np.arange(len(windows), dtype=np.int64), lens + 1)
        self.owner_at = np.repeat(np.asarray(owners, dtype=np.int64), lens + 1)
        self.pos = np.arange(n, dtype=np.int64) - np.repeat(self.starts, lens + 1)

    def table(
        self, length: int
    ) -> tuple[
        np.ndarray[Any, Any],
        np.ndarray[Any, Any],
        np.ndarray[Any, Any],
        np.ndarray[Any, Any],
    ]:
        """``(hash, window, owner, start)`` of every in-window L-mer, sorted by hash, owner."""
        p = np.arange(self.n - length + 1, dtype=np.int64)
        p = p[self.pos[p] + length <= self.lens[self.window_at[p]]]
        with np.errstate(over="ignore"):
            h = (self.cum[p + length] - self.cum[p]) * self.pw[self.n - p]
        own = self.owner_at[p]
        order = np.lexsort((own, h))
        return h[order], self.window_at[p][order], own[order], p[order]

    @staticmethod
    def groups(
        h: np.ndarray[Any, Any], own: np.ndarray[Any, Any]
    ) -> tuple[np.ndarray[Any, Any], np.ndarray[Any, Any]]:
        """Per element: its hash group's first and last index (groups sorted by owner)."""
        edge = np.concatenate([[True], h[1:] != h[:-1]])
        gid = np.cumsum(edge) - 1
        first = np.nonzero(edge)[0]
        last = np.concatenate([first[1:], [len(h)]]) - 1
        return first[gid], last[gid]

    def shared(self, length: int) -> np.ndarray[Any, Any]:
        """Per window: a window of a different owner sharing one of its L-mers, or -1."""
        h, win, own, _ = self.table(length)
        out = np.full(len(self.windows), -1, dtype=np.int64)
        if len(h):
            f, last = self.groups(h, own)
            multi = own[f] != own[last]
            other = np.where(own[f] != own, win[f], win[last])
            out[win[multi]] = other[multi]
        return out

    def confirm(self, i: int, j: int, length: int) -> bool:
        """Window ``j`` really holds a length-``length`` substring of window ``i``."""
        wi, wj = self.windows[i], self.windows[j]
        return any(wi[q : q + length] in wj for q in range(len(wi) - length + 1))


def longest_shared_across_pool(
    windows: list[str], owners: list[int] | None = None
) -> tuple[list[int], list[int | None]]:
    """Per window, the longest substring it shares with a window of another owner, and that partner.

    Pass ``owners`` to ignore matches between arms of the same member (same owner id); when
    omitted, every window is its own owner. Uses a coarse ladder pass followed by a parallel
    binary search over ``TractIndex``, and confirms every winning pair against the raw strings
    so a 64-bit hash collision raises rather than reporting a false match.
    """
    if not windows:
        return [], []
    own = list(range(len(windows))) if owners is None else owners
    if len(own) != len(windows):
        raise ValueError("owners must have the same length as windows")
    idx = TractIndex(windows, own)
    n = len(windows)
    lo = np.zeros(n, dtype=np.int64)
    hi = np.array([len(w) + 1 for w in windows], dtype=np.int64)
    partner = np.full(n, -1, dtype=np.int64)

    def trial(length: int, who: np.ndarray[Any, Any]) -> None:
        got = idx.shared(length)[who]
        ok = got >= 0
        lo[who[ok]] = length
        partner[who[ok]] = got[ok]
        hi[who[~ok]] = np.minimum(hi[who[~ok]], length)

    for length in range(_TRACT_LADDER_STEP, int(hi.max()), _TRACT_LADDER_STEP):
        who = np.nonzero((lo < length) & (length < hi))[0]
        if len(who):
            trial(length, who)
    while True:
        active = np.nonzero(hi - lo > 1)[0]
        if len(active) == 0:
            break
        mid = (lo[active] + hi[active]) // 2
        for m in np.unique(mid):
            trial(int(m), active[mid == m])
    out: list[int | None] = []
    for i in range(n):
        j = int(partner[i])
        if lo[i] > 0 and j >= 0 and not idx.confirm(i, j, int(lo[i])):
            raise AssertionError(f"hash collision between windows {i} and {j}")
        out.append(j if lo[i] > 0 and j >= 0 else None)
    return lo.tolist(), out


def forbidden_patterns(enzymes: list[str]) -> list[str]:
    """Both-strand recognition patterns for the requested enzymes."""
    pats: list[str] = []
    for name in enzymes:
        site = ENZYME_SITES[name]
        pats.append(site)
        rc = reverse_complement(site)
        if rc != site:
            pats.append(rc)
    return pats


# --------------------------------------------------------------------------- #
# Pure sequence operations
# --------------------------------------------------------------------------- #


def split_codons(cds: str) -> list[str]:
    """Split an in-frame CDS into codons, validating length."""
    cds = cds.upper()
    if len(cds) % 3 != 0:
        raise ValueError(f"CDS length {len(cds)} is not a multiple of 3")
    return [cds[i : i + 3] for i in range(0, len(cds), 3)]


def translate(cds: str) -> str:
    """Translate an in-frame CDS to a one-letter protein string."""
    return "".join(GENETIC_CODE[c] for c in split_codons(cds))


def read_fasta(path: str | os.PathLike[str]) -> list[tuple[str, str]]:
    """``[(header, sequence), ...]`` in file order: whole header without ``>``, case kept.

    Duplicate headers are kept, so a caller that needs unique ids uses :func:`fasta_by_id`.
    Sequence lines are stripped and joined; text before the first header is ignored.
    """
    out: list[tuple[str, str]] = []
    header: str | None = None
    buf: list[str] = []
    with open(path) as fh:
        for line in fh.read().splitlines():
            if line.startswith(">"):
                if header is not None:
                    out.append((header, "".join(buf)))
                header, buf = line[1:], []
            elif header is not None:
                buf.append(line.strip())
    if header is not None:
        out.append((header, "".join(buf)))
    return out


def fasta_id(header: str) -> str:
    """A FASTA record's id: the header up to the first ``|`` or whitespace."""
    return (header.split("|")[0].split() or [""])[0]


def fasta_by_id(path: str | os.PathLike[str]) -> dict[str, str]:
    """``{id: sequence}`` keyed by :func:`fasta_id`, case kept; a repeated id raises ``ValueError``."""
    out: dict[str, str] = {}
    for header, seq in read_fasta(path):
        key = fasta_id(header)
        if key in out:
            raise ValueError(f"{path}: duplicate FASTA id {key}")
        out[key] = seq
    return out


def common_suffix(a: str, b: str) -> int:
    """Length of the longest common suffix of two sequences.

    The swap-tract measure: how far a template switch just before the barcode can run back
    into two inserts without changing either.
    """
    n = 0
    while n < min(len(a), len(b)) and a[-1 - n] == b[-1 - n]:
        n += 1
    return n


def gc_content(seq: str) -> float:
    """Fraction of G/C bases (0-1)."""
    if not seq:
        return 0.0
    return (seq.count("G") + seq.count("C")) / len(seq)


def find_sites(seq: str, enzymes: list[str]) -> list[tuple[str, str, int]]:
    """Return ``(enzyme, strand, 0-based position)`` for every recognition site.

    ``strand`` is ``"+"`` when the forward recognition sequence matches and
    ``"-"`` when its reverse complement matches.
    """
    hits: list[tuple[str, str, int]] = []
    for name in enzymes:
        site = ENZYME_SITES[name]
        rc = reverse_complement(site)
        for pattern, strand in ((site, "+"), (rc, "-")):
            if strand == "-" and rc == site:
                continue
            start = 0
            while True:
                idx = seq.find(pattern, start)
                if idx == -1:
                    break
                hits.append((name, strand, idx))
                start = idx + 1
    return sorted(hits, key=lambda h: h[2])


CODON_TABLE_SPECIES = "h_sapiens"  # python_codon_tables key for Homo sapiens

# %GC window + homopolymer convention: matches the 35-65% / 50-nt-window
# rule used for the COMBINE domain-library synthesis (Moon, Herschl et al.,
# Nat Commun 2024) — see the follow-up note in PR #10 for the source.
DEFAULT_GC_MIN = 0.35
DEFAULT_GC_MAX = 0.65
DEFAULT_GC_WINDOW = 50
# Longest allowed run (forbids runs of 6+). Looser than the barcode
# generator's max_homopolymer=2: a barcode is synthetic and fully free to
# avoid any repeat, but a real protein's amino acid sequence can force short
# homopolymers regardless of codon choice (e.g. adjacent Lys/Asn codons
# share an "AA" prefix) — max_homopolymer=3 (forbid 4+) proved infeasible on
# real genes (e.g. TFRC) as a hard constraint; 5 solves cleanly in practice.
DEFAULT_MAX_HOMOPOLYMER = 5
DEFAULT_KMER_UNIQUENESS_K = 9


def dnachisel_optimize(
    cds: str,
    enzymes: list[str],
    gc_min: float = DEFAULT_GC_MIN,
    gc_max: float = DEFAULT_GC_MAX,
    gc_window: int = DEFAULT_GC_WINDOW,
    max_homopolymer: int = DEFAULT_MAX_HOMOPOLYMER,
    kmer_uniqueness_k: int = DEFAULT_KMER_UNIQUENESS_K,
    seed: int | None = None,
) -> tuple[str, list[str]]:
    """Codon-optimize + domesticate a CDS in one DNA Chisel solve.

    Unlike a sequential "optimize, then scrub sites" pass, this jointly
    satisfies codon usage, restriction-site avoidance, a windowed %GC target,
    and homopolymer avoidance, with k-mer uniqueness as a secondary objective
    (soft — see below). ``EnforceTranslation`` guarantees the protein is
    unchanged; the caller should still re-verify this (see ``recode_gene``),
    treating it as a correctness gate rather than trusting DNA Chisel alone.

    ``kmer_uniqueness_k`` is deliberately an *objective*, not a hard
    constraint: a hard uniqueness requirement is unsatisfiable for genuinely
    repetitive protein regions (e.g. a run of several Glycines has only 4
    synonymous codons to draw from), which would raise DNA Chisel's
    ``NoSolutionError`` for entirely legitimate sequences. As an objective it
    still meaningfully reduces repeats without ever blocking a solve.

    Returns the optimized CDS and a list of human-readable warnings (empty on
    full success; a failed hard constraint is reported here rather than
    raised, mirroring the shape callers expect from a domestication pass).
    """
    if seed is not None:
        np.random.seed(seed)

    unknown = set(enzymes) - set(ENZYME_SITES)
    if unknown:
        raise ValueError(f"Unknown enzyme(s): {sorted(unknown)}. Known: {sorted(ENZYME_SITES)}")

    homopolymer_run = max_homopolymer + 1  # AvoidPattern forbids this exact run length
    constraints = [
        dc.EnforceTranslation(),
        dc.EnforceGCContent(mini=gc_min, maxi=gc_max, window=gc_window),
        *(dc.AvoidPattern(dc.EnzymeSitePattern(name)) for name in enzymes),
        *(dc.AvoidPattern(dc.HomopolymerPattern(base, homopolymer_run)) for base in "ATGC"),
    ]
    objectives = [
        dc.CodonOptimize(species=CODON_TABLE_SPECIES, method="match_codon_usage"),
        dc.UniquifyAllKmers(k=kmer_uniqueness_k, boost=1.0),
    ]

    problem = dc.DnaOptimizationProblem(
        sequence=cds, constraints=constraints, objectives=objectives, logger=None
    )
    warnings: list[str] = []
    try:
        problem.resolve_constraints()
    except dc.NoSolutionError as err:
        warnings.append(f"DNA Chisel could not satisfy all constraints: {err}")

    try:
        problem.optimize()
    except dc.NoSolutionError as err:
        # Objective optimization can occasionally re-violate a constraint
        # that resolve_constraints() had already satisfied (e.g. a local
        # exhaustive search for CodonOptimize nudging a windowed %GC target
        # out of bounds). Re-resolve before giving up, since the sequence at
        # this point is otherwise almost fully optimized.
        warnings.append(f"DNA Chisel optimization step hit a conflict: {err}")
        try:
            problem.resolve_constraints()
        except dc.NoSolutionError as err2:
            warnings.append(f"DNA Chisel could not re-resolve after optimization: {err2}")

    return problem.sequence, warnings


def domesticate(
    sequence: str,
    enzymes: list[str],
    preserve_translation: bool = True,
    seed: int | None = None,
) -> tuple[str, list[str]]:
    """Remove restriction sites with the fewest possible nucleotide changes.

    Unlike :func:`dnachisel_optimize`, this does **not** codon-optimize —
    there is no ``CodonOptimize`` objective, no GC-content window, and no
    homopolymer/repeat handling. It only does two things: (1) forbid the
    given enzymes' recognition sites (both strands) as hard constraints, and
    (2) minimize the total number of nucleotide edits from the input via
    DNA Chisel's ``AvoidChanges`` objective. The result is the original
    sequence with the smallest possible set of point mutations needed to
    clear every site.

    Use this when you want to keep a native/original sequence as close to
    unchanged as possible — e.g. domesticating a real gene's native codon
    usage, or a non-coding element (an IRES, a UTR) where there's no protein
    to preserve and codon optimization doesn't apply at all.

    ``preserve_translation=True`` (the default) adds an ``EnforceTranslation``
    constraint, so this only works on an in-frame CDS in that mode; pass
    ``False`` for non-coding sequences (the input needn't be a multiple of 3).

    Returns the domesticated sequence and a list of human-readable warnings
    (empty on full success).
    """
    if seed is not None:
        np.random.seed(seed)

    unknown = set(enzymes) - set(ENZYME_SITES)
    if unknown:
        raise ValueError(f"Unknown enzyme(s): {sorted(unknown)}. Known: {sorted(ENZYME_SITES)}")

    constraints: list[Any] = [dc.AvoidPattern(dc.EnzymeSitePattern(name)) for name in enzymes]
    if preserve_translation:
        constraints.insert(0, dc.EnforceTranslation())

    problem = dc.DnaOptimizationProblem(
        sequence=sequence,
        constraints=constraints,
        objectives=[dc.AvoidChanges()],
        logger=None,
    )
    warnings: list[str] = []
    try:
        problem.resolve_constraints()
    except dc.NoSolutionError as err:
        warnings.append(f"DNA Chisel could not satisfy all constraints: {err}")

    try:
        problem.optimize()
    except dc.NoSolutionError as err:
        warnings.append(f"DNA Chisel optimization step hit a conflict: {err}")
        try:
            problem.resolve_constraints()
        except dc.NoSolutionError as err2:
            warnings.append(f"DNA Chisel could not re-resolve after optimization: {err2}")

    return problem.sequence, warnings


def domesticate_gene(
    symbol: str,
    enzymes: list[str] | None = None,
    species: str = "homo_sapiens",
    seed: int | None = None,
) -> RecodeResult:
    """End-to-end: symbol -> stable IDs -> native CDS -> minimal-edit domestication.

    Unlike :func:`recode_gene`, the native codon usage is left alone — this
    only removes restriction sites, via :func:`domesticate`, changing as few
    nucleotides as possible.
    """
    enzymes = enzymes or list(ENZYME_SITES)
    unknown = set(enzymes) - set(ENZYME_SITES)
    if unknown:
        raise ValueError(f"Unknown enzyme(s): {sorted(unknown)}. Known: {sorted(ENZYME_SITES)}")

    gene = resolve_gene(symbol, species=species)
    original_cds = fetch_cds(gene.ensembl_transcript_id)
    protein = translate(original_cds)

    domesticated, warnings = domesticate(
        original_cds, enzymes, preserve_translation=True, seed=seed
    )

    # Correctness gate: the protein must be byte-for-byte identical.
    if translate(domesticated) != protein:
        raise RuntimeError("BUG: domestication changed the protein sequence")

    return RecodeResult(
        gene=gene,
        protein=protein,
        original_cds=original_cds,
        optimized_cds=domesticated,
        enzymes=enzymes,
        sites_before=find_sites(original_cds, enzymes),
        sites_after=find_sites(domesticated, enzymes),
        gc_before=gc_content(original_cds),
        gc_after=gc_content(domesticated),
        warnings=warnings,
        method="domesticated",
    )


# --------------------------------------------------------------------------- #
# Barcode generation
# --------------------------------------------------------------------------- #

DEFAULT_BARCODE_AVOID_MOTIFS: tuple[str, ...] = (
    "AATAAA",  # canonical polyA signal - don't create one inside a barcode
)


def has_homopolymer_run(seq: str, max_run: int) -> bool:
    """True if ``seq`` contains a run of the same base longer than ``max_run``.

    ``max_run=2`` (the default used by :func:`generate_barcode`) disallows any
    run of 3+ identical bases, i.e. no ``AAA``/``TTT``/``GGG``/``CCC``.

    Every run is measured, including the runs of length 1 contributed by a
    base that differs from its neighbours. ``max_run=0`` is therefore
    unsatisfiable for any non-empty sequence, which is what makes it a valid
    way to express "these constraints cannot be met".
    """
    run = 0
    previous = ""
    for base in seq:
        run = run + 1 if base == previous else 1
        if run > max_run:
            return True
        previous = base
    return False


@dataclass
class BarcodeResult:
    """A generated barcode plus the constraints it was drawn to satisfy."""

    sequence: str
    length: int
    gc_content: float
    max_homopolymer: int
    enzymes: list[str]
    avoid_motifs: list[str]
    label: str | None = None
    attempts: int = 0

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    def fasta(self) -> str:
        header = f">barcode|{self.label or 'unlabeled'}|{self.length}nt|GC{self.gc_content:.2f}"
        return f"{header}\n{self.sequence}\n"


def generate_barcode(
    length: int = 20,
    gc_min: float = 0.40,
    gc_max: float = 0.60,
    max_homopolymer: int = 2,
    enzymes: list[str] | None = None,
    avoid_motifs: list[str] | None = None,
    avoid_start_atg: bool = True,
    existing: set[str] | None = None,
    label: str | None = None,
    max_attempts: int = 200_000,
    rng: random.Random | None = None,
) -> BarcodeResult:
    """Draw a random DNA barcode meeting length/GC/homopolymer/site constraints.

    Rejection sampling: draw ``length`` random bases, check every constraint,
    retry on failure. Fails loudly (``RuntimeError``) rather than silently
    returning a barcode that violates a constraint.

    Parameters
    ----------
    length:
        Barcode length in nt. 20 is a reasonable default for a pooled library
        of tens of thousands of constructs — long enough to keep collisions
        negligible even with heavy per-construct barcode redundancy, short
        enough to read cleanly in one amplicon.
    gc_min, gc_max:
        Inclusive %GC window (0-1). Extreme GC content hurts synthesis and
        sequencing quality; 0.40-0.60 is a safe general-purpose window.
    max_homopolymer:
        Longest allowed run of an identical base (``2`` disallows any run of
        3+). Homopolymers hurt synthesis fidelity and can cause sequencer
        slippage/miscalls.
    enzymes:
        Names from :data:`ENZYME_SITES` whose recognition sites (both
        strands) must not appear in the barcode. Defaults to ``["BbsI",
        "Esp3I"]`` to match this module's default Golden Gate enzymes.
    avoid_motifs:
        Additional literal motifs to exclude (checked as given; pass the
        reverse complement too if you need both strands covered). Defaults
        to the canonical polyA signal ``AATAAA``.
    avoid_start_atg:
        If true (default), reject barcodes starting with ``ATG`` so a
        barcode placed just after a stop codon can't read as a spurious
        start codon.
    existing:
        A set of previously used barcodes to avoid, for building up a
        non-colliding pool one call at a time.
    label:
        Free-text label recorded on the result (e.g. a gene/construct name).
    rng:
        Inject a seeded :class:`random.Random` for deterministic/testable
        output; a fresh one is created otherwise.
    """
    enzymes = enzymes if enzymes is not None else ["BbsI", "Esp3I"]
    avoid_motifs = (
        list(avoid_motifs) if avoid_motifs is not None else list(DEFAULT_BARCODE_AVOID_MOTIFS)
    )
    existing = existing or set()
    rng = rng or random.Random()

    unknown = set(enzymes) - set(ENZYME_SITES)
    if unknown:
        raise ValueError(f"Unknown enzyme(s): {sorted(unknown)}. Known: {sorted(ENZYME_SITES)}")
    if not 0.0 <= gc_min <= gc_max <= 1.0:
        raise ValueError(f"Invalid GC window: gc_min={gc_min}, gc_max={gc_max}")

    site_patterns = forbidden_patterns(enzymes)
    motif_patterns = [m.upper() for m in avoid_motifs]
    bases = "ACGT"

    for attempt in range(1, max_attempts + 1):
        seq = "".join(rng.choice(bases) for _ in range(length))
        gc = gc_content(seq)
        if not (gc_min <= gc <= gc_max):
            continue
        if has_homopolymer_run(seq, max_homopolymer):
            continue
        if avoid_start_atg and seq.startswith("ATG"):
            continue
        if any(p in seq for p in site_patterns):
            continue
        if any(m in seq for m in motif_patterns):
            continue
        if seq in existing:
            continue
        return BarcodeResult(
            sequence=seq,
            length=length,
            gc_content=gc,
            max_homopolymer=max_homopolymer,
            enzymes=enzymes,
            avoid_motifs=avoid_motifs,
            label=label,
            attempts=attempt,
        )

    raise RuntimeError(
        f"Could not generate a {length}-nt barcode meeting constraints after "
        f"{max_attempts} attempts; relax the GC window, homopolymer limit, "
        "or motif/enzyme list."
    )


# --------------------------------------------------------------------------- #
# Linker diversification (recombination avoidance)
# --------------------------------------------------------------------------- #

# Flexible linkers are the most repeat-prone sequence in a tagged construct.
# ``GGGGS`` back-translates to codons that all share a ``GG`` dinucleotide, so
# the naive encoding of ``GGGGS-tag-GGGGSGGGGS-tag-GGGGS`` repeats one 15-nt
# block four times inside 153 nt. Direct repeats of that size are a substrate
# for replication slippage during plasmid propagation, and for strand transfer
# during lentiviral reverse transcription.
#
# Two failure modes, which need different fixes.
#
# *Within* one construct, a repeated block can be deleted along with everything
# between two copies. :func:`diversify_linker` holds the longest internal
# direct repeat below ``max_repeat``.
#
# *Across* a pooled library, every member carrying a byte-identical cassette
# creates a long shared tract between each ORF and its barcode. Two genomes
# co-packaged in one lentiviral particle can template-switch inside that tract,
# which **uncouples the barcode from the ORF** — the failure that silently
# reassigns phenotypes to the wrong gene. :func:`diversify_library_linkers`
# gives every member its own encoding, under two constraints that target the
# two ways a tract can form:
#
# 1. **Core islands.** No k-mer lying wholly inside the cassette is shared with
#    any other member, so members whose ORFs differ cannot acquire a new shared
#    stretch from the cassette itself.
# 2. **Junction run-through.** For members whose flanking ORF sequence *is*
#    shared — paralogs, isoforms, orthologs — the flank identity cannot be
#    designed away; it is the gene. What can be designed is how far that tract
#    continues past the insertion point. Each member's cassette is required to
#    diverge from every same-flank partner within ``max_junction_identity`` nt
#    of each junction, which terminates the tract inside the cassette instead
#    of letting it run on toward the barcode.
#
# Constraint 2 is why flanking homology sets the priority order: members with
# the longest shared flanks have the least room and are served first, while the
# prefix space is still empty.
#
# Note the limit on what constraint 2 can achieve. Every glycine codon is
# ``GGN``, so a cassette beginning with Gly shares its first base with every
# other encoding no matter what is chosen; divergence can begin at base 2 but
# not base 1. Requiring zero shared junction-spanning k-mers is therefore
# unsatisfiable by construction, not merely hard — which is why the constraint
# is framed as bounded run-through rather than as k-mer exclusion.

DEFAULT_RECOMBINOGENIC_K = 20
"""Shared-tract length treated as recombinogenic.

RecA-dependent homologous recombination in *E. coli* needs a minimal efficient
processing segment in the low tens of bp and climbs steeply above it; retroviral
strand transfer is promoted by shorter tracts still. 20 nt is a deliberately
conservative working threshold for both, not a measured constant — lower it for
a more aggressive design, raise it if the constraint proves infeasible.
"""

DEFAULT_MAX_JUNCTION_IDENTITY = 8
"""How far a cassette may stay identical to a same-flank partner's, per junction.

The reachable floor is set by the codon structure at the cassette's own ends,
not by this number: with Gly at the boundary there are only ``4 ** (n // 3)``
distinct n-nt prefixes, so 8 affords 64 — enough to separate the paralog pairs
and small ortholog groups seen in practice, and not enough for a large family.
:meth:`LibraryLinkerPlan.audit` reports what was actually achieved.
"""

DEFAULT_MIN_CODON_USAGE = 0.10
"""Drop synonymous codons rarer than this fraction of their amino acid's usage.

Diversity is worth nothing if it is bought with codons the host translates
badly. The most frequent codon for an amino acid is always retained even if it
falls below the threshold, so the table can never be emptied.
"""

DEFAULT_LINKER_MAX_HOMOPOLYMER = 4
DEFAULT_MAX_INTERNAL_REPEAT = 9
DEFAULT_JUNCTION_SEED = 12
# Glycine codons are all GGN, so a Gly-rich linker has a hard %GC floor near
# 0.60 that no codon choice can go under (GGGGS bottoms out at 9/15). The
# window below brackets what is reachable for a whole cassette; the generic
# DEFAULT_GC_MIN/MAX used for entire CDSs is infeasible for a bare GS block.
# Call achievable_gc() before tightening either bound.
DEFAULT_LINKER_GC_MIN = 0.45
DEFAULT_LINKER_GC_MAX = 0.70


def codon_usage(species: str = CODON_TABLE_SPECIES) -> dict[str, dict[str, float]]:
    """Relative synonymous codon usage, as ``{amino_acid: {codon: fraction}}``."""
    import python_codon_tables

    return cast(dict[str, dict[str, float]], python_codon_tables.get_codons_table(species))


def synonymous_codons(
    min_usage: float = DEFAULT_MIN_CODON_USAGE,
    species: str = CODON_TABLE_SPECIES,
) -> dict[str, list[tuple[str, float]]]:
    """``{amino_acid: [(codon, usage_fraction), ...]}``, rare codons removed.

    The most frequent codon for each amino acid is kept unconditionally, so
    every amino acid retains at least one option however high ``min_usage`` is.
    """
    table = codon_usage(species)
    out: dict[str, list[tuple[str, float]]] = {}
    for aa, codons in table.items():
        ranked = sorted(codons.items(), key=lambda kv: kv[1], reverse=True)
        kept = [(c, f) for c, f in ranked if f >= min_usage] or ranked[:1]
        out[aa] = kept
    return out


def max_internal_repeat(seq: str) -> int:
    """Length of the longest substring occurring more than once in ``seq``.

    Direct repeats only — these are the homologous-recombination and slippage
    substrate. Inverted repeats (hairpins) are a different problem and are not
    measured here.
    """
    n = len(seq)
    if n < 2:
        return 0

    def repeats(k: int) -> bool:
        if k <= 0:
            return True
        seen: set[str] = set()
        for i in range(n - k + 1):
            sub = seq[i : i + k]
            if sub in seen:
                return True
            seen.add(sub)
        return False

    lo, hi = 0, n - 1
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if repeats(mid):
            lo = mid
        else:
            hi = mid - 1
    return lo


def core_kmers(core: str, k: int) -> set[str]:
    """Every k-mer lying wholly inside ``core``.

    Junction-spanning k-mers are deliberately excluded: they are partly fixed
    by the flanking ORF, so excluding them would make same-flank members
    unsatisfiable rather than merely constrained. Run-through across a junction
    is bounded by ``max_junction_identity`` instead.
    """
    return {core[i : i + k] for i in range(len(core) - k + 1)}


def pin_codons(reference: str, peptide: str, segments: Iterable[str]) -> dict[int, str]:
    """Codons to hold constant: every residue of each ``segments`` peptide, read off ``reference``.

    ``reference`` is an existing, validated encoding of ``peptide`` -- for the tag cassette,
    p003's -- and each segment (a tag such as SpyTag003 or ALFA) must occur exactly once in
    ``peptide``. The result is the ``fixed_codons`` argument of :func:`diversify_linker` and
    :func:`diversify_library_linkers`: those codons are copied verbatim into every
    encoding, so only the residues outside the segments are recoded.
    """
    reference = reference.upper()
    if len(reference) != 3 * len(peptide) or translate(reference) != peptide:
        raise ValueError("reference is not an encoding of peptide")
    fixed: dict[int, str] = {}
    for seg in segments:
        start = peptide.find(seg)
        if start < 0 or peptide.find(seg, start + 1) >= 0:
            raise ValueError(f"segment {seg!r} must occur exactly once in the peptide")
        for i in range(start, start + len(seg)):
            fixed[i] = reference[3 * i : 3 * i + 3]
    return fixed


def _check_fixed(peptide: str, fixed_codons: Mapping[int, str] | None) -> dict[int, str]:
    fixed = {int(i): c.upper() for i, c in (fixed_codons or {}).items()}
    for i, codon in fixed.items():
        if not 0 <= i < len(peptide):
            raise ValueError(f"fixed codon index {i} is outside the {len(peptide)}-aa peptide")
        if GENETIC_CODE.get(codon) != peptide[i]:
            raise ValueError(f"fixed codon {codon} at {i} does not encode {peptide[i]!r}")
    return fixed


def variable_kmers(dna: str, k: int, fixed_codons: Mapping[int, str] | None = None) -> set[str]:
    """:func:`core_kmers` less the windows lying wholly inside fixed codons.

    A window inside a pinned tag is identical in every encoding by design, so counting it
    would make every pair of members look equally overlapped and leave ``"minimize"``
    nothing to rank. Only the windows a codon choice can change are compared.
    """
    if not fixed_codons:
        return core_kmers(dna, k)
    pinned = bytearray(len(dna))
    for i in fixed_codons:
        pinned[3 * i : 3 * i + 3] = b"\x01\x01\x01"
    return {dna[i : i + k] for i in range(len(dna) - k + 1) if not all(pinned[i : i + k])}


def achievable_gc(
    peptide: str,
    min_usage: float = DEFAULT_MIN_CODON_USAGE,
    species: str = CODON_TABLE_SPECIES,
) -> tuple[float, float]:
    """Lowest and highest %GC any synonymous encoding of ``peptide`` can reach.

    Glycine is ``GGN``, so a Gly-rich linker cannot be brought down to a typical
    whole-CDS GC window however the codons are chosen. Call this before setting
    ``gc_min``/``gc_max`` rather than discovering it as a sampling failure.
    """
    table = synonymous_codons(min_usage=min_usage, species=species)
    lo = hi = 0
    for aa in peptide:
        counts = [c.count("G") + c.count("C") for c, _ in table[aa]]
        lo += min(counts)
        hi += max(counts)
    n = 3 * len(peptide)
    return (lo / n, hi / n) if n else (0.0, 0.0)


def _gc_tilted_weights(
    choices: list[list[str]],
    weights: list[list[float]],
    target_gc: float,
) -> list[list[float]]:
    """Reweight synonymous codons so the expected %GC lands on ``target_gc``.

    An exponential tilt ``w' = w * exp(-lam * gc)`` with one ``lam`` shared
    across every position: the minimum-relative-entropy adjustment to the usage
    distribution that meets the GC constraint in expectation, so it distorts
    codon usage as little as the target allows. Without it, rejection sampling
    on a Gly-rich peptide spends thousands of draws fighting the usage
    distribution's own GC bias (measured: ~18,000 draws and a 45% failure rate
    for 4x GGGGS at gc_max=0.65, versus ~13 draws with the tilt).
    """
    gcs = [[c.count("G") + c.count("C") for c in opts] for opts in choices]
    total = 3 * len(choices)

    def expected(lam: float) -> float:
        acc = 0.0
        for w, g in zip(weights, gcs, strict=True):
            tilted = [wi * math.exp(-lam * gi) for wi, gi in zip(w, g, strict=True)]
            z = sum(tilted)
            acc += sum(ti * gi for ti, gi in zip(tilted, g, strict=True)) / z
        return acc / total

    lo, hi = -20.0, 20.0
    if expected(hi) > target_gc or expected(lo) < target_gc:
        return weights  # target outside the reachable range; let rejection handle it
    for _ in range(60):
        mid = (lo + hi) / 2
        if expected(mid) > target_gc:
            lo = mid
        else:
            hi = mid
    lam = (lo + hi) / 2
    return [
        [wi * math.exp(-lam * gi) for wi, gi in zip(w, g, strict=True)]
        for w, g in zip(weights, gcs, strict=True)
    ]


def _end_codon_map(
    choices: list[list[str]], edge: int, from_start: bool
) -> tuple[dict[str, list[tuple[str, ...]]], list[int]]:
    """Map each reachable ``edge``-nt end of the encoding to the codons giving it."""
    ncod = -(-edge // 3)
    if ncod > len(choices):
        return {}, []
    idxs = list(range(ncod)) if from_start else list(range(len(choices) - ncod, len(choices)))
    opts = [choices[i] for i in idxs]
    total = 1
    for o in opts:
        total *= len(o)
    if total > 100_000:  # not worth enumerating; rejection sampling will cope
        return {}, []
    out: dict[str, list[tuple[str, ...]]] = {}
    for combo in itertools.product(*opts):
        joined = "".join(combo)
        key = joined[:edge] if from_start else joined[-edge:]
        out.setdefault(key, []).append(combo)
    return out, idxs


def end_capacity(
    peptide: str,
    max_junction_identity: int = DEFAULT_MAX_JUNCTION_IDENTITY,
    min_usage: float = DEFAULT_MIN_CODON_USAGE,
    species: str = CODON_TABLE_SPECIES,
) -> tuple[int, int]:
    """How many members sharing a flank this peptide can be given distinct ends for.

    Returns ``(distinct_5p_ends, distinct_3p_ends)``. This is the hard ceiling on
    the size of a paralog/isoform/ortholog group that
    :func:`diversify_library_linkers` can separate at a junction: a group larger
    than this cannot be given a unique run-through, whatever the sampling does.
    Raise ``max_junction_identity``, or put a non-Gly residue at the cassette
    boundary, to buy more.
    """
    table = synonymous_codons(min_usage=min_usage, species=species)
    choices = [[c for c, _ in table[aa]] for aa in peptide]
    edge = max_junction_identity + 1
    five, _ = _end_codon_map(choices, edge, True)
    three, _ = _end_codon_map(choices, edge, False)
    return len(five), len(three)


@dataclass
class LinkerVariant:
    """One synonymous encoding of a linker, plus the properties it was drawn for."""

    sequence: str
    peptide: str
    gc_content: float
    max_internal_repeat: int
    enzymes: list[str]
    recombinogenic_k: int
    shared_kmers: int = 0
    label: str | None = None
    attempts: int = 0

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _prepare_linker_codons(
    peptide: str,
    *,
    enzymes: list[str],
    gc_min: float,
    gc_max: float,
    kmer_policy: str,
    candidates: int,
    min_codon_usage: float,
    species: str,
    fixed_codons: Mapping[int, str] | None,
    pad5: str,
    pad3: str,
) -> tuple[
    dict[int, str],
    list[list[str]],
    list[list[float]],
    float,
    float,
    list[str],
]:
    """Validate ``diversify_linker`` arguments, compute codon choices/weights, and pre-flight ends."""
    unknown = set(enzymes) - set(ENZYME_SITES)
    if unknown:
        raise ValueError(f"Unknown enzyme(s): {sorted(unknown)}. Known: {sorted(ENZYME_SITES)}")
    if "*" in peptide:
        raise ValueError("peptide must not contain a stop codon")
    if not 0.0 <= gc_min <= gc_max <= 1.0:
        raise ValueError(f"Invalid GC window: gc_min={gc_min}, gc_max={gc_max}")
    if kmer_policy not in ("forbid", "minimize"):
        raise ValueError(f"kmer_policy must be 'forbid' or 'minimize', got {kmer_policy!r}")
    if candidates < 1:
        raise ValueError("candidates must be >= 1")

    table = synonymous_codons(min_usage=min_codon_usage, species=species)
    missing = sorted(set(peptide) - set(table))
    if missing:
        raise ValueError(f"No codons for residue(s): {missing}")
    fixed = _check_fixed(peptide, fixed_codons)

    choices = [
        [fixed[i]] if i in fixed else [c for c, _ in table[aa]] for i, aa in enumerate(peptide)
    ]
    usage = [[1.0] if i in fixed else [f for _, f in table[aa]] for i, aa in enumerate(peptide)]
    gcs = [[c.count("G") + c.count("C") for c in opts] for opts in choices]
    reach_lo = sum(min(g) for g in gcs) / (3 * len(peptide)) if peptide else 0.0
    reach_hi = sum(max(g) for g in gcs) / (3 * len(peptide)) if peptide else 0.0
    if reach_lo > gc_max or reach_hi < gc_min:
        raise ValueError(
            f"{peptide!r} can only reach %GC {reach_lo:.2f}-{reach_hi:.2f}, which does "
            f"not meet the requested window {gc_min:.2f}-{gc_max:.2f}. Glycine is GGN, "
            "so Gly-rich linkers have a high GC floor."
        )

    weights = _gc_tilted_weights(
        choices,
        usage,
        (max(gc_min, reach_lo) + min(gc_max, reach_hi)) / 2,
    )
    site_patterns = forbidden_patterns(enzymes)
    for side, pad, opts in (("5'", pad5, choices[0]), ("3'", pad3, choices[-1])):
        if not pad:
            continue
        blocked = [
            c
            for c in opts
            if any(s in (pad + c if side == "5'" else c + pad) for s in site_patterns)
        ]
        if len(blocked) == len(opts):
            raise RuntimeError(
                f"every codon available at the {side} end of {peptide!r} completes a Type IIS "
                f"site against the flanking sequence {pad!r}. No codon choice inside the "
                f"cassette can fix this — domesticate the flanking ORF with domesticate() so "
                f"the {len(pad)} bases at the junction change, then re-run."
            )
    return fixed, choices, weights, reach_lo, reach_hi, site_patterns


def _free_end_samplers(
    peptide: str,
    choices: list[list[str]],
    edge: int,
    max_junction_identity: int,
    avoid_prefixes: set[str],
    avoid_suffixes: set[str],
) -> tuple[
    list[str],
    dict[str, list[tuple[str, ...]]],
    list[int],
    list[str],
    dict[str, list[tuple[str, ...]]],
    list[int],
]:
    """Direct samplers for unclaimed 5'/3' boundary sequences when partners hold ends."""
    free5: list[str] = []
    free3: list[str] = []
    map5: dict[str, list[tuple[str, ...]]] = {}
    idx5: list[int] = []
    map3: dict[str, list[tuple[str, ...]]] = {}
    idx3: list[int] = []
    for side, is_prefix, avoid in (
        ("5'", True, avoid_prefixes),
        ("3'", False, avoid_suffixes),
    ):
        if not avoid:
            continue
        cmap, idxs = _end_codon_map(choices, edge, is_prefix)
        if not cmap:
            continue
        free = sorted(set(cmap) - avoid)
        if not free:
            raise RuntimeError(
                f"all {len(cmap)} reachable {edge}-nt {side} ends of {peptide!r} are "
                f"already taken by same-flank partners. That is the combinatorial "
                f"ceiling, not a sampling failure: raise max_junction_identity "
                f"(currently {max_junction_identity}) or change the residue at the "
                f"cassette's {side} boundary."
            )
        if is_prefix:
            free5, map5, idx5 = free, cmap, idxs
        else:
            free3, map3, idx3 = free, cmap, idxs
    if set(idx5) & set(idx3):  # the two forced ends overlap; fall back to rejection
        return [], {}, [], [], {}, []
    return free5, map5, idx5, free3, map3, idx3


def _linker_failure_message(
    peptide: str,
    max_attempts: int,
    rejected: dict[str, int],
    *,
    gc_min: float,
    gc_max: float,
    reach_lo: float,
    reach_hi: float,
    max_homopolymer: int,
    edge: int,
    max_junction_identity: int,
    max_repeat: int,
) -> str:
    """Diagnostic message naming the dominant rejection cause when sampling exhausts attempts."""
    culprit = max(rejected, key=lambda key: rejected[key])
    breakdown = ", ".join(f"{key}={count}" for key, count in rejected.items() if count)
    advice = {
        "gc": (
            f"every draw fell outside %GC {gc_min:.2f}-{gc_max:.2f}; this peptide can only "
            f"reach {reach_lo:.2f}-{reach_hi:.2f}"
        ),
        "homopolymer": f"max_homopolymer={max_homopolymer} is too tight for this peptide",
        "type_iis": (
            "every draw formed a Type IIS site across a junction. The cassette's first base "
            "is fixed by its first codon — every Gly codon is GGN — so a flank ending in the "
            "five bases that complete a site cannot be rescued by codon choice inside the "
            "cassette. Domesticate the flanking ORF with domesticate() and re-run"
        ),
        "end_taken": (
            f"the {edge}-nt ends still free were all rejected by the other constraints; raise "
            f"max_junction_identity (currently {max_junction_identity}) to widen the choice"
        ),
        "repeat": f"max_repeat={max_repeat} is too tight for a peptide this repetitive",
    }[culprit]
    return (
        f"Could not encode {peptide!r} ({len(peptide)} aa) in {max_attempts} attempts. "
        f"Rejections: {breakdown}. Dominant cause is {culprit!r}: {advice}."
    )


def diversify_linker(
    peptide: str,
    *,
    flank5: str = "",
    flank3: str = "",
    avoid_kmers: set[str] | None = None,
    avoid_prefixes: set[str] | None = None,
    avoid_suffixes: set[str] | None = None,
    fixed_codons: Mapping[int, str] | None = None,
    k: int = DEFAULT_RECOMBINOGENIC_K,
    kmer_policy: str = "forbid",
    candidates: int = 8,
    max_junction_identity: int = DEFAULT_MAX_JUNCTION_IDENTITY,
    enzymes: list[str] | None = None,
    gc_min: float = DEFAULT_LINKER_GC_MIN,
    gc_max: float = DEFAULT_LINKER_GC_MAX,
    max_homopolymer: int = DEFAULT_LINKER_MAX_HOMOPOLYMER,
    max_repeat: int = DEFAULT_MAX_INTERNAL_REPEAT,
    min_codon_usage: float = DEFAULT_MIN_CODON_USAGE,
    species: str = CODON_TABLE_SPECIES,
    label: str | None = None,
    max_attempts: int = 20_000,
    rng: random.Random | None = None,
) -> LinkerVariant:
    """Draw a synonymous DNA encoding of ``peptide`` that avoids shared tracts.

    Rejection sampling over GC-tilted, usage-weighted synonymous codons, in the same style as
    :func:`generate_barcode`: draw, check every constraint, retry; raises if none succeeds.
    Pass ``kmer_policy="minimize"`` to keep the best of ``candidates`` valid draws at library
    scale, ``avoid_prefixes``/``avoid_suffixes`` to prevent shared flanking tracts from running
    through the cassette, and ``fixed_codons`` (from :func:`pin_codons`) to hold tag codons
    constant while recoding only the linkers.
    """
    enzymes = enzymes if enzymes is not None else list(ENZYME_SITES)
    avoid_kmers = avoid_kmers or set()
    avoid_prefixes = avoid_prefixes or set()
    avoid_suffixes = avoid_suffixes or set()
    rng = rng or random.Random()
    pad5, pad3 = flank5.upper()[-5:], flank3.upper()[:5]
    edge = max_junction_identity + 1

    fixed, choices, weights, reach_lo, reach_hi, site_patterns = _prepare_linker_codons(
        peptide,
        enzymes=enzymes,
        gc_min=gc_min,
        gc_max=gc_max,
        kmer_policy=kmer_policy,
        candidates=candidates,
        min_codon_usage=min_codon_usage,
        species=species,
        fixed_codons=fixed_codons,
        pad5=pad5,
        pad3=pad3,
    )
    free5, map5, idx5, free3, map3, idx3 = _free_end_samplers(
        peptide, choices, edge, max_junction_identity, avoid_prefixes, avoid_suffixes
    )

    best: LinkerVariant | None = None
    pool = 0
    rejected = dict.fromkeys(("gc", "homopolymer", "type_iis", "end_taken", "repeat"), 0)

    for attempt in range(1, max_attempts + 1):
        codons = [
            rng.choices(opts, weights=w, k=1)[0] for opts, w in zip(choices, weights, strict=True)
        ]
        for free, cmap, idxs in ((free5, map5, idx5), (free3, map3, idx3)):
            if free:
                for pos, codon in zip(idxs, rng.choice(cmap[rng.choice(free)]), strict=True):
                    codons[pos] = codon
        dna = "".join(codons)
        gc = gc_content(dna)
        if not (gc_min <= gc <= gc_max):
            rejected["gc"] += 1
            continue
        if has_homopolymer_run(dna, max_homopolymer):
            rejected["homopolymer"] += 1
            continue
        if any(p in pad5 + dna + pad3 for p in site_patterns):
            rejected["type_iis"] += 1
            continue
        if (avoid_prefixes and dna[:edge] in avoid_prefixes) or (
            avoid_suffixes and dna[-edge:] in avoid_suffixes
        ):
            rejected["end_taken"] += 1
            continue
        rep = max_internal_repeat(dna)
        if rep > max_repeat:
            rejected["repeat"] += 1
            continue
        shared = len(variable_kmers(dna, k, fixed) & avoid_kmers) if avoid_kmers else 0
        if shared and kmer_policy == "forbid":
            continue

        found = LinkerVariant(
            sequence=dna,
            peptide=peptide,
            gc_content=gc,
            max_internal_repeat=rep,
            enzymes=enzymes,
            recombinogenic_k=k,
            shared_kmers=shared,
            label=label,
            attempts=attempt,
        )
        if shared == 0:
            return found
        if best is None or shared < best.shared_kmers:
            best = found
        pool += 1
        if pool >= candidates:
            return best

    if best is not None:
        return best
    raise RuntimeError(
        _linker_failure_message(
            peptide,
            max_attempts,
            rejected,
            gc_min=gc_min,
            gc_max=gc_max,
            reach_lo=reach_lo,
            reach_hi=reach_hi,
            max_homopolymer=max_homopolymer,
            edge=edge,
            max_junction_identity=max_junction_identity,
            max_repeat=max_repeat,
        )
    )


def _common_suffix(a: str, b: str) -> int:
    n = 0
    for x, y in zip(reversed(a), reversed(b), strict=False):
        if x != y:
            break
        n += 1
    return n


def _common_prefix(a: str, b: str) -> int:
    n = 0
    for x, y in zip(a, b, strict=False):
        if x != y:
            break
        n += 1
    return n


def _junction_pairs(
    members: dict[str, tuple[str, str]],
    seed: int,
) -> tuple[dict[tuple[str, str], int], dict[tuple[str, str], int]]:
    """Pairwise flank identity, for pairs that match at least ``seed`` nt.

    Members are bucketed on the ``seed``-mer touching the insertion point and
    compared only within buckets, which keeps this near-linear instead of
    all-pairs. Shorter matches are not recombinogenic on their own.
    """
    buckets: tuple[dict[str, list[str]], dict[str, list[str]]] = ({}, {})
    for name, (f5, f3) in members.items():
        if len(f5) >= seed:
            buckets[0].setdefault(f5[-seed:].upper(), []).append(name)
        if len(f3) >= seed:
            buckets[1].setdefault(f3[:seed].upper(), []).append(name)

    pairs: list[dict[tuple[str, str], int]] = [{}, {}]
    for side, (bucket, fn) in enumerate(
        ((buckets[0], _common_suffix), (buckets[1], _common_prefix))
    ):
        for group in bucket.values():
            for i, a in enumerate(group):
                for b in group[i + 1 :]:
                    m = fn(members[a][side].upper(), members[b][side].upper())
                    pairs[side][(a, b)] = pairs[side][(b, a)] = m
    return pairs[0], pairs[1]


@dataclass
class JunctionHomology:
    """How much sequence a library member shares with its nearest neighbour.

    ``spanning_match`` is the tract a *single* partner contributes on both sides
    of the insertion point at once — the identity that an identical cassette
    encoding would weld into one continuous stretch. That is the quantity that
    sets design priority. ``flank5_match`` and ``flank3_match`` are the best
    one-sided matches, which may be with different partners.
    """

    name: str
    flank5_match: int
    flank3_match: int
    spanning_match: int
    partner: str | None = None

    @property
    def risk(self) -> tuple[int, int]:
        """Sort key, descending: spanning tract first, then worst single side."""
        return (self.spanning_match, max(self.flank5_match, self.flank3_match))

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def junction_homology(
    members: dict[str, tuple[str, str]],
    seed: int = DEFAULT_JUNCTION_SEED,
) -> dict[str, JunctionHomology]:
    """Measure flanking sequence identity around each member's insertion point.

    ``members`` maps a name to ``(flank5, flank3)`` — the sequence immediately
    5' and 3' of where the cassette goes, in that member's own construct.
    Matches shorter than ``seed`` nt adjacent to the junction are reported as 0.
    """
    if seed < 1:
        raise ValueError("seed must be >= 1")
    pair5, pair3 = _junction_pairs(members, seed)

    partners: dict[str, set[str]] = {n: set() for n in members}
    for store in (pair5, pair3):
        for a, b in store:
            partners[a].add(b)

    out: dict[str, JunctionHomology] = {}
    for name in members:
        peers = partners[name]
        best_span, best_partner = 0, None
        for p in peers:
            span = pair5.get((name, p), 0) + pair3.get((name, p), 0)
            if span > best_span:
                best_span, best_partner = span, p
        out[name] = JunctionHomology(
            name=name,
            flank5_match=max((pair5.get((name, p), 0) for p in peers), default=0),
            flank3_match=max((pair3.get((name, p), 0) for p in peers), default=0),
            spanning_match=best_span,
            partner=best_partner,
        )
    return out


@dataclass
class LibraryLinkerPlan:
    """Per-member linker encodings plus the homology analysis that ordered them."""

    variants: dict[str, LinkerVariant]
    homology: dict[str, JunctionHomology]
    flanks: dict[str, tuple[str, str]]
    order: list[str]
    recombinogenic_k: int
    max_junction_identity: int
    seed: int
    peptide: str
    fixed_codons: dict[int, str] = field(default_factory=dict)

    def audit(self) -> dict[str, int]:
        """Recompute, from the issued sequences, what the plan actually achieved.

        Returned as measured worst cases rather than as the requested bounds, so
        a hand-edited or partially reused plan is checked rather than trusted:

        ``shared_core_kmers``
            largest number of wholly-internal k-mers any two cassettes share.
        ``longest_shared_substring``
            longest stretch any two cassettes have in common. This is the
            number to compare against ``recombinogenic_k``; the k-mer count
            above is only a proxy for how much they overlap in total.
        ``junction_runthrough_5`` / ``junction_runthrough_3``
            furthest a shared flanking tract continues into the cassette, over
            every pair of members whose flanks match on that side.
        ``spanning_tract``
            longest single continuous identity between any two members through
            the insertion point: flank + cassette + flank.
        ``max_internal_repeat``
            longest direct repeat inside any one cassette.

        With ``fixed_codons`` the pinned tags are shared by construction, so
        ``longest_shared_substring`` is at least the longest tag; the k-mer count
        excludes windows wholly inside them, as the draw did.
        """
        pair5, pair3 = _junction_pairs(self.flanks, self.seed)
        seq = {n: v.sequence for n, v in self.variants.items()}

        islands = 0
        seen: list[set[str]] = []
        for name in self.order:
            km = variable_kmers(seq[name], self.recombinogenic_k, self.fixed_codons)
            for other in seen:
                islands = max(islands, len(km & other))
            seen.append(km)

        # Longest substring any two cassettes share, found by bisecting on
        # length with one k-mer dictionary per probe — exact, and far cheaper
        # than comparing all pairs.
        def any_shared(width: int) -> bool:
            pool: set[str] = set()
            for s in seq.values():
                sub = {s[i : i + width] for i in range(len(s) - width + 1)}
                if sub & pool:
                    return True
                pool |= sub
            return False

        lo, hi = 0, min(len(s) for s in seq.values())
        while lo < hi:
            mid = (lo + hi + 1) // 2
            if any_shared(mid):
                lo = mid
            else:
                hi = mid - 1
        longest_shared = lo

        run5 = run3 = span = 0
        for (a, b), flank in pair5.items():
            core = _common_prefix(seq[a], seq[b])
            run5 = max(run5, core)
            if core < len(seq[a]):  # tract terminates inside the cassette
                span = max(span, flank + core)
            else:  # cassette identical end to end; tract continues into flank3
                span = max(span, flank + core + pair3.get((a, b), 0))
        for (a, b), _ in pair3.items():
            run3 = max(run3, _common_suffix(seq[a], seq[b]))

        return {
            "members": len(self.variants),
            "distinct_encodings": len({v.sequence for v in self.variants.values()}),
            "shared_core_kmers": islands,
            "longest_shared_substring": longest_shared,
            "junction_runthrough_5": run5,
            "junction_runthrough_3": run3,
            "spanning_tract": span,
            "max_internal_repeat": max(v.max_internal_repeat for v in self.variants.values()),
        }

    def to_rows(self) -> list[dict[str, Any]]:
        """One flat record per member, for writing to TSV."""
        rows = []
        for name in self.order:
            v, h = self.variants[name], self.homology[name]
            rows.append(
                {
                    "member": name,
                    "linker_dna": v.sequence,
                    "linker_aa": v.peptide,
                    "gc_content": round(v.gc_content, 4),
                    "max_internal_repeat": v.max_internal_repeat,
                    "draw_attempts": v.attempts,
                    "flank5_match": h.flank5_match,
                    "flank3_match": h.flank3_match,
                    "spanning_match": h.spanning_match,
                    "nearest_member": h.partner or "",
                }
            )
        return rows


def diversify_library_linkers(
    members: dict[str, tuple[str, str]],
    peptide: str,
    *,
    k: int = DEFAULT_RECOMBINOGENIC_K,
    seed: int = DEFAULT_JUNCTION_SEED,
    max_junction_identity: int = DEFAULT_MAX_JUNCTION_IDENTITY,
    kmer_policy: str = "minimize",
    rng: random.Random | None = None,
    existing: dict[str, str] | None = None,
    fixed_codons: Mapping[int, str] | None = None,
    **linker_kwargs: Any,
) -> LibraryLinkerPlan:
    """Issue one distinct cassette encoding per library member.

    Members are served in descending order of flanking homology, so the ones
    with the least sequence of their own to interrupt a shared tract choose from
    an empty space and the already-distinct ones take what is left. Each issued
    encoding must share no wholly-internal ``k``-mer with any other member, and
    must diverge from every same-flank partner within ``max_junction_identity``
    nt of each junction.

    ``members`` maps a name to ``(flank5, flank3)`` — the sequence either side of
    the insertion point in that member's own construct. Supply at least ``k`` nt
    of each so junction homology can be measured; 150 nt is ample.

    ``existing`` pins encodings that have already been issued: pass
    ``{name: dna}`` and those members keep exactly that sequence, while the rest
    are drawn around them. **Use it whenever a member is added to a library whose
    constructs are already built.** Without it every call re-derives the whole
    set, and because members are served in risk order, inserting one changes the
    order and the accumulated k-mer pool — so previously issued encodings come
    back different. Eight built cassettes plus two new members re-issued seven of
    the eight. A re-issued encoding silently diverges the design file from the
    DNA already in a freezer.

    ``fixed_codons`` (from :func:`pin_codons`) holds those residues on one encoding in
    every newly drawn member, so only the rest of the peptide varies. Pinned
    ``existing`` encodings are exempt: an already-built construct keeps whatever
    codons it was built with.

    Pinned sequences are not re-checked against the flank and enzyme constraints
    they were originally drawn under; they seed the k-mer pool and the junction
    exclusions so the *new* members avoid them. Call :meth:`LibraryLinkerPlan.audit`
    to measure what the combined set actually achieves.

    Raises ``RuntimeError`` naming the member that could not be satisfied, so a
    failure says which construct exhausted the space rather than only that
    something did. Check the result with :meth:`LibraryLinkerPlan.audit`.
    """
    pinned = dict(existing or {})
    stray = sorted(set(pinned) - set(members))
    if stray:
        raise ValueError(
            f"existing names are not members, so their flanks are unknown and the "
            f"new draws cannot avoid their junctions: {stray}"
        )
    for name, dna in pinned.items():
        if translate(dna) != peptide:
            raise ValueError(
                f"{name}: existing encoding translates to a different peptide than "
                f"the one being issued, so it is not a synonymous variant of it"
            )

    rng = rng or random.Random()
    fixed = _check_fixed(peptide, fixed_codons)
    homology = junction_homology(members, seed=seed)
    pair5, pair3 = _junction_pairs(members, seed)
    edge = max_junction_identity + 1

    issued_kmers: set[str] = set()
    variants: dict[str, LinkerVariant] = {}
    for name, dna in pinned.items():
        variants[name] = LinkerVariant(
            sequence=dna,
            peptide=peptide,
            gc_content=gc_content(dna),
            max_internal_repeat=max_internal_repeat(dna),
            enzymes=list(linker_kwargs.get("enzymes") or []),
            recombinogenic_k=k,
            label=name,
        )
        issued_kmers |= variable_kmers(dna, k, fixed)

    # Only undrawn members are ordered and served; a pinned one is never redrawn,
    # which is what makes the result stable as the library grows.
    order = sorted(
        (n for n in members if n not in pinned),
        key=lambda n: homology[n].risk,
        reverse=True,
    )
    for name in order:
        f5, f3 = members[name]
        taken5 = {variants[p].sequence[:edge] for p in variants if (name, p) in pair5}
        taken3 = {variants[p].sequence[-edge:] for p in variants if (name, p) in pair3}
        try:
            var = diversify_linker(
                peptide,
                flank5=f5,
                flank3=f3,
                avoid_kmers=issued_kmers,
                kmer_policy=kmer_policy,
                avoid_prefixes=taken5,
                avoid_suffixes=taken3,
                fixed_codons=fixed,
                k=k,
                max_junction_identity=max_junction_identity,
                label=name,
                rng=rng,
                **linker_kwargs,
            )
        except (RuntimeError, ValueError) as exc:
            raise RuntimeError(
                f"{name}: ran out of encoding space after {len(variants)} of "
                f"{len(members)} members were issued (spanning flank homology "
                f"{homology[name].spanning_match} nt with "
                f"{homology[name].partner or 'no partner'}; {len(taken5)} prefixes "
                f"and {len(taken3)} suffixes already taken by its own partners). {exc}"
            ) from exc
        variants[name] = var
        issued_kmers |= variable_kmers(var.sequence, k, fixed)

    return LibraryLinkerPlan(
        variants=variants,
        homology=homology,
        flanks=dict(members),
        order=list(pinned) + order,
        recombinogenic_k=k,
        max_junction_identity=max_junction_identity,
        seed=seed,
        peptide=peptide,
        fixed_codons=fixed,
    )


# --------------------------------------------------------------------------- #
# Ensembl REST access (latest GRCh38)
# --------------------------------------------------------------------------- #

ENSEMBL_REST = "https://rest.ensembl.org"
USER_AGENT = "deliverome-analysis seqtools (research; https://github.com/Deliverome-Project)"


def _ensembl_request(path: str, params: dict[str, str], accept: str) -> str:
    """GET an Ensembl REST endpoint with 429/5xx retry + Retry-After handling."""
    url = f"{ENSEMBL_REST}{path}"
    if params:
        url = f"{url}?{urllib.parse.urlencode(params)}"
    retryable = {429, 500, 502, 503, 504}
    last_err: Exception | None = None
    for attempt in range(5):
        req = urllib.request.Request(url, headers={"Accept": accept, "User-Agent": USER_AGENT})
        try:
            with urllib.request.urlopen(req, timeout=30) as resp:
                return cast(bytes, resp.read()).decode("utf-8")
        except urllib.error.HTTPError as err:
            last_err = err
            if err.code in retryable:
                wait = float(err.headers.get("Retry-After", 2**attempt))
                time.sleep(min(wait, 8.0))
                continue
            raise
        except urllib.error.URLError as err:
            last_err = err
            time.sleep(min(2**attempt, 8.0))
    raise RuntimeError(f"Ensembl REST request failed after retries: {url}") from last_err


@dataclass
class GeneRecord:
    """Stable-ID provenance for a resolved gene + its canonical transcript."""

    symbol: str
    ensembl_gene_id: str  # unversioned ENSG
    ensembl_gene_id_versioned: str
    ensembl_transcript_id: str  # canonical ENST (versioned)
    ensembl_protein_id: str | None
    assembly: str
    biotype: str


def resolve_gene(symbol: str, species: str = "homo_sapiens") -> GeneRecord:
    """Resolve a gene symbol to its ENSG + canonical transcript stable IDs."""
    payload = json.loads(
        _ensembl_request(
            f"/lookup/symbol/{species}/{symbol}",
            {"expand": "1"},
            accept="application/json",
        )
    )
    gene_id = payload["id"]
    version = payload.get("version")
    transcripts = payload.get("Transcript", [])
    canonical = next((t for t in transcripts if t.get("is_canonical")), None)
    if canonical is None:
        raise ValueError(f"No canonical transcript found for {symbol} ({gene_id})")
    enst = canonical["id"]
    if canonical.get("version") is not None:
        enst = f"{enst}.{canonical['version']}"
    translation = canonical.get("Translation") or {}
    return GeneRecord(
        symbol=symbol,
        ensembl_gene_id=gene_id,
        ensembl_gene_id_versioned=f"{gene_id}.{version}" if version is not None else gene_id,
        ensembl_transcript_id=enst,
        ensembl_protein_id=translation.get("id"),
        assembly=payload.get("assembly_name", "unknown"),
        biotype=canonical.get("biotype", "unknown"),
    )


def fetch_cds(transcript_id: str) -> str:
    """Fetch the CDS nucleotide sequence for an Ensembl transcript stable ID.

    The ``sequence/id`` endpoint rejects versioned IDs and always returns the
    current version, so any ``.N`` suffix is stripped before the request.
    """
    unversioned = transcript_id.split(".", 1)[0]
    fasta = _ensembl_request(
        f"/sequence/id/{unversioned}",
        {"type": "cds"},
        accept="text/x-fasta",
    )
    lines = [ln.strip() for ln in fasta.splitlines() if ln and not ln.startswith(">")]
    seq = "".join(lines).upper()
    if not seq or set(seq) - set("ACGTN"):
        raise ValueError(f"Unexpected CDS content for {transcript_id}")
    return seq


# --------------------------------------------------------------------------- #
# Orchestration
# --------------------------------------------------------------------------- #


@dataclass
class RecodeResult:
    """Everything a caller needs: stable IDs, sequences, and a change report."""

    gene: GeneRecord
    protein: str
    original_cds: str
    optimized_cds: str
    enzymes: list[str]
    sites_before: list[tuple[str, str, int]]
    sites_after: list[tuple[str, str, int]]
    gc_before: float
    gc_after: float
    warnings: list[str] = field(default_factory=list)
    method: str = "codon_optimized"  # or "domesticated" (see domesticate_gene)

    @property
    def ok(self) -> bool:
        """True when the protein is preserved and no forbidden sites remain."""
        return not self.sites_after and not self.warnings

    def to_dict(self) -> dict[str, Any]:
        d = asdict(self)
        d["gene"] = asdict(self.gene)
        d["ok"] = self.ok
        return d

    def fasta(self) -> str:
        header = (
            f">{self.gene.symbol}|{self.gene.ensembl_gene_id_versioned}|"
            f"{self.gene.ensembl_transcript_id}|{self.method}|"
            f"no_{'_'.join(self.enzymes)}"
        )
        return f"{header}\n{self.optimized_cds}\n"


def recode_gene(
    symbol: str,
    enzymes: list[str] | None = None,
    species: str = "homo_sapiens",
    seed: int | None = None,
) -> RecodeResult:
    """End-to-end: symbol -> stable IDs -> CDS -> DNA Chisel optimize + domesticate.

    ``seed`` seeds DNA Chisel's (numpy-backed) RNG for reproducible output —
    useful for tests; leave unset for normal use.
    """
    enzymes = enzymes or list(ENZYME_SITES)
    unknown = set(enzymes) - set(ENZYME_SITES)
    if unknown:
        raise ValueError(f"Unknown enzyme(s): {sorted(unknown)}. Known: {sorted(ENZYME_SITES)}")

    gene = resolve_gene(symbol, species=species)
    original_cds = fetch_cds(gene.ensembl_transcript_id)
    protein = translate(original_cds)

    optimized, warnings = dnachisel_optimize(original_cds, enzymes, seed=seed)

    # Correctness gate: the protein must be byte-for-byte identical. DNA
    # Chisel's EnforceTranslation constraint should already guarantee this;
    # treat any violation as a bug rather than trusting the library blindly.
    if translate(optimized) != protein:
        raise RuntimeError("BUG: codon optimization changed the protein sequence")

    return RecodeResult(
        gene=gene,
        protein=protein,
        original_cds=original_cds,
        optimized_cds=optimized,
        enzymes=enzymes,
        sites_before=find_sites(original_cds, enzymes),
        sites_after=find_sites(optimized, enzymes),
        gc_before=gc_content(original_cds),
        gc_after=gc_content(optimized),
        warnings=warnings,
    )


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="deliverome_analysis.seqtools",
        description="Gene-symbol codon recoding and constraint-based barcode generation.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    recode = subparsers.add_parser(
        "recode",
        help="Resolve a gene symbol to its Ensembl canonical CDS, codon-optimize it, "
        "and remove BbsI/Esp3I sites.",
    )
    recode.add_argument("symbol", help="Gene symbol, e.g. TF")
    recode.add_argument(
        "--enzymes",
        default="BbsI,Esp3I",
        help="Comma-separated enzymes to scrub (default: BbsI,Esp3I)",
    )
    recode.add_argument(
        "--species", default="homo_sapiens", help="Ensembl species (default: homo_sapiens)"
    )
    recode.add_argument("--output", help="Write optimized CDS FASTA to this path (default: stdout)")
    recode.add_argument(
        "--json", dest="json_path", help="Write the full report as JSON to this path"
    )

    domesticate_cmd = subparsers.add_parser(
        "domesticate",
        help="Resolve a gene symbol to its Ensembl canonical CDS and remove restriction "
        "sites with the fewest possible nucleotide changes (no codon optimization).",
    )
    domesticate_cmd.add_argument("symbol", help="Gene symbol, e.g. TF")
    domesticate_cmd.add_argument(
        "--enzymes",
        default="BbsI,Esp3I,BsaI",
        help="Comma-separated enzymes to scrub (default: BbsI,Esp3I,BsaI; "
        "Esp3I and BsmBI are the same recognition site)",
    )
    domesticate_cmd.add_argument(
        "--species", default="homo_sapiens", help="Ensembl species (default: homo_sapiens)"
    )
    domesticate_cmd.add_argument("--seed", type=int, help="Seed the RNG for reproducible output")
    domesticate_cmd.add_argument(
        "--output", help="Write domesticated CDS FASTA to this path (default: stdout)"
    )
    domesticate_cmd.add_argument(
        "--json", dest="json_path", help="Write the full report as JSON to this path"
    )

    barcode = subparsers.add_parser(
        "barcode",
        help="Generate a random DNA barcode meeting length/GC/homopolymer/enzyme constraints.",
    )
    barcode.add_argument(
        "--length", type=int, default=20, help="Barcode length in nt (default: 20)"
    )
    barcode.add_argument(
        "--gc-min", type=float, default=0.40, help="Minimum fraction GC, 0-1 (default: 0.40)"
    )
    barcode.add_argument(
        "--gc-max", type=float, default=0.60, help="Maximum fraction GC, 0-1 (default: 0.60)"
    )
    barcode.add_argument(
        "--max-homopolymer",
        type=int,
        default=2,
        help="Longest allowed run of an identical base (default: 2, i.e. no run of 3+)",
    )
    barcode.add_argument(
        "--enzymes",
        default="BbsI,Esp3I",
        help="Comma-separated enzyme sites to avoid (default: BbsI,Esp3I)",
    )
    barcode.add_argument(
        "--avoid-motifs",
        default="AATAAA",
        help="Comma-separated extra motifs to avoid (default: AATAAA, the polyA signal)",
    )
    barcode.add_argument(
        "--allow-atg-start",
        action="store_true",
        help="Allow the barcode to start with ATG (rejected by default)",
    )
    barcode.add_argument("--label", help="Free-text label recorded on the result, e.g. a gene name")
    barcode.add_argument("--seed", type=int, help="Seed the RNG for reproducible output")
    barcode.add_argument("--output", help="Write barcode FASTA to this path (default: stdout)")
    barcode.add_argument(
        "--json", dest="json_path", help="Write the full report as JSON to this path"
    )

    return parser


def _run_recode(args: argparse.Namespace) -> int:
    enzymes = [e.strip() for e in args.enzymes.split(",") if e.strip()]
    result = recode_gene(args.symbol, enzymes=enzymes, species=args.species)

    fasta = result.fasta()
    if args.output:
        with open(args.output, "w") as fh:
            fh.write(fasta)
    else:
        sys.stdout.write(fasta)

    if args.json_path:
        with open(args.json_path, "w") as fh:
            json.dump(result.to_dict(), fh, indent=2)

    g = result.gene
    print(
        f"\n{g.symbol}  {g.ensembl_gene_id_versioned}  {g.ensembl_transcript_id}"
        f"  ({g.assembly}, {g.biotype})",
        file=sys.stderr,
    )
    print(
        f"  protein: {len(result.protein)} aa | CDS: {len(result.optimized_cds)} nt | "
        f"GC {result.gc_before:.1%} -> {result.gc_after:.1%}",
        file=sys.stderr,
    )
    print(
        f"  sites removed: {len(result.sites_before)} -> {len(result.sites_after)} "
        f"({', '.join(enzymes)})",
        file=sys.stderr,
    )
    for w in result.warnings:
        print(f"  WARNING: {w}", file=sys.stderr)
    return 0 if result.ok else 1


def _run_domesticate(args: argparse.Namespace) -> int:
    enzymes = [e.strip() for e in args.enzymes.split(",") if e.strip()]
    result = domesticate_gene(args.symbol, enzymes=enzymes, species=args.species, seed=args.seed)

    fasta = result.fasta()
    if args.output:
        with open(args.output, "w") as fh:
            fh.write(fasta)
    else:
        sys.stdout.write(fasta)

    if args.json_path:
        with open(args.json_path, "w") as fh:
            json.dump(result.to_dict(), fh, indent=2)

    g = result.gene
    print(
        f"\n{g.symbol}  {g.ensembl_gene_id_versioned}  {g.ensembl_transcript_id}"
        f"  ({g.assembly}, {g.biotype})",
        file=sys.stderr,
    )
    n_changes = sum(
        1 for a, b in zip(result.original_cds, result.optimized_cds, strict=True) if a != b
    )
    print(
        f"  protein: {len(result.protein)} aa | CDS: {len(result.optimized_cds)} nt | "
        f"nucleotide changes: {n_changes} | GC {result.gc_before:.1%} -> {result.gc_after:.1%}",
        file=sys.stderr,
    )
    print(
        f"  sites removed: {len(result.sites_before)} -> {len(result.sites_after)} "
        f"({', '.join(enzymes)})",
        file=sys.stderr,
    )
    for w in result.warnings:
        print(f"  WARNING: {w}", file=sys.stderr)
    return 0 if result.ok else 1


def _run_barcode(args: argparse.Namespace) -> int:
    enzymes = [e.strip() for e in args.enzymes.split(",") if e.strip()]
    avoid_motifs = [m.strip() for m in args.avoid_motifs.split(",") if m.strip()]
    rng = random.Random(args.seed) if args.seed is not None else None

    result = generate_barcode(
        length=args.length,
        gc_min=args.gc_min,
        gc_max=args.gc_max,
        max_homopolymer=args.max_homopolymer,
        enzymes=enzymes,
        avoid_motifs=avoid_motifs,
        avoid_start_atg=not args.allow_atg_start,
        label=args.label,
        rng=rng,
    )

    fasta = result.fasta()
    if args.output:
        with open(args.output, "w") as fh:
            fh.write(fasta)
    else:
        sys.stdout.write(fasta)

    if args.json_path:
        with open(args.json_path, "w") as fh:
            json.dump(result.to_dict(), fh, indent=2)

    print(
        f"\nbarcode: {result.sequence}  ({result.length} nt, GC {result.gc_content:.1%}, "
        f"{result.attempts} attempt(s))",
        file=sys.stderr,
    )
    print(
        f"  constraints: GC [{args.gc_min:.0%}-{args.gc_max:.0%}], "
        f"max homopolymer {result.max_homopolymer}, avoiding {', '.join(result.enzymes)}"
        + (f" + {', '.join(result.avoid_motifs)}" if result.avoid_motifs else ""),
        file=sys.stderr,
    )
    return 0


def main(argv: list[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    if args.command == "recode":
        return _run_recode(args)
    if args.command == "domesticate":
        return _run_domesticate(args)
    if args.command == "barcode":
        return _run_barcode(args)
    raise AssertionError(f"unreachable: unknown command {args.command!r}")  # pragma: no cover


if __name__ == "__main__":
    raise SystemExit(main())
