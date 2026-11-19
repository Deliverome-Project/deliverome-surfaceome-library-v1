"""s1e2 — generate the production barcode set: 30,000 x 20 nt, full published filter stack.

This is the deliverable, not an analysis: it emits the actual barcodes to
`barcodes_20nt_30k.tsv`. Deterministic given SEED, so the file is reproducible from this
script alone.

FILTERS APPLIED (each traced in barcode_filter_stack.py, which documents the source and
whether we verified it):

  GC 40-60%                                     Hawkins/FREE, Shiroguchi
  no homopolymer triple (max run 2)             Hawkins/FREE
  no ATG anywhere                               ours -- avoids a spurious start codon
  no AATAAA + all 18 single-base point mutants  Slutskin
  no 3-bp self-complementarity, non-overlapping Hawkins/FREE (hairpin rule)
  no conserved+ miRNA 6-mer seed target site    TargetScanHuman 8.0, 207 sites. Slutskin
                                                defines the seed as bases 2-7, and every
                                                7mer-m8 target CONTAINS its 6-mer core, so
                                                the 6-mer rejects a strict superset. We were
                                                the less conservative of the two until this.
  no GT[AG]AG 5' splice-donor core            ours -- a cryptic donor in the 3'UTR would
                                                splice out the IRES and the marker. ADDED
                                                2026-09-26; see f_no_splice_donor
  no BbsI / Esp3I-BsmBI / BsaI, both strands    non-palindromic, so they need naming; every
                                                PALINDROMIC 6-cutter is already excluded as a
                                                side effect of the hairpin rule (see below)
  min Levenshtein distance 3                    Slutskin: differ "by 3 nt or more"
                                                (PMC5802814); corrects one edit incl. indels
  AND min Hamming distance 4                    matches the stricter substitution separation
                                                used by published UMI generators. Levenshtein
                                                3 alone permits Hamming-3 pairs -- the first
                                                build had 71 of them.

WHY BOTH DISTANCE METRICS. Hamming and Levenshtein are not nested in the useful direction.
For equal-length strings Levenshtein <= Hamming, so a Hamming-4 code gives NO indel guarantee
(two codewords related by an insertion plus a deletion can read as many substitutions but few
edits), while a Levenshtein-3 code permits Hamming-3 pairs. Published UMI generators typically
enforce Hamming only; FREE exists precisely because that leaves indels uncovered. Enforcing
both costs almost nothing here and closes both gaps.

PALINDROMIC RESTRICTION SITES COME FREE. Any palindromic 6-cutter P satisfies revcomp(P[0:3])
== P[3:6], i.e. its two halves are non-overlapping reverse complements -- exactly what the
3-bp self-complementarity rule rejects. A survey of 30 common 6-cutters found zero in the
first build, before any of them was named. Only non-palindromic sites (our Type IIS enzymes)
need an explicit blacklist.

NOT applied, with reasons in README 4.4: the GC-dinucleotide rule (costs 75% of survivors and
is a 2018 Illumina-chemistry concern) and all-human seed tiers (62%, mostly miRNAs not
expressed in any given line).

DISTANCE GUARANTEE. Enforced by deletion neighbourhoods: for equal-length strings, if
Levenshtein(a,b) <= k then a and b share a common subsequence with exactly k characters
deleted. Rejecting any candidate that shares a 2-deletion subsequence with an accepted one
therefore guarantees min distance >= 3. The condition is necessary but not sufficient, so the
filter over-rejects -- it never admits a too-close pair. A sampled pairwise check runs at the
end as an independent sanity test.

Run:
    uv run python analysis/sprint1/s1e2-fragmid-vector-design/generate_barcodes.py
"""

import csv
import itertools
import random
from pathlib import Path

from barcode_filter_stack import (
    f_gc,
    f_no_atg,
    f_no_enzyme,
    f_no_polya,
    f_no_selfcomp3,
    f_no_splice_donor,
    f_no_triples,
    load_seeds,
    make_seed_filter,
    polyadenylation_hexamers,
)

HERE = Path(__file__).resolve().parent

LENGTH = 20
TARGET = 30_000
MIN_DISTANCE = 3  # Levenshtein
MIN_HAMMING = 4
SEED = 20260812
OUT = HERE / f"barcodes_{LENGTH}nt_{TARGET // 1000}k.tsv"


def _pack(s):
    v = 0
    for ch in s:
        v = (v << 2) | "ACGT".index(ch)
    return v


def levenshtein(a, b):
    """Plain DP, used only for the sampled verification at the end."""
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i]
        for j, cb in enumerate(b, 1):
            cur.append(min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (ca != cb)))
        prev = cur
    return prev[-1]


def main():
    rng = random.Random(SEED)
    polya = polyadenylation_hexamers()
    seeds = load_seeds({"broadly_conserved", "conserved"}, 6)
    seed_ok = make_seed_filter(seeds, 6)

    def ok(bc):
        return (
            f_gc(bc)
            and f_no_triples(bc)
            and f_no_atg(bc)
            and f_no_polya(bc, polya)
            and f_no_enzyme(bc)
            and f_no_selfcomp3(bc)
            and f_no_splice_donor(bc)
            and seed_ok(bc)
        )

    print(f"generating {TARGET:,} x {LENGTH} nt, min Levenshtein {MIN_DISTANCE}")
    print(f"  miRNA 6-mer seed sites excluded: {len(seeds)} (conserved+, Slutskin definition)")
    print(f"  min Levenshtein {MIN_DISTANCE} AND min Hamming {MIN_HAMMING}")
    print(f"  polyA hexamers excluded:   {len(polya)}\n")

    k = MIN_DISTANCE - 1
    seen = set()
    barcodes = []
    # Hamming index. If Hamming(a,b) <= MIN_HAMMING-1 = 3, then splitting 20 positions into
    # 4 blocks of 5 guarantees by pigeonhole that at least one block is identical. Bucketing
    # on each block therefore finds every too-close pair while comparing almost nothing.
    blocks = [(0, 5), (5, 10), (10, 15), (15, 20)]
    hidx = [{} for _ in blocks]

    def hamming_ok(bc):
        for (bs, be), idx in zip(blocks, hidx, strict=True):
            for other in idx.get(bc[bs:be], ()):
                if sum(x != y for x, y in zip(bc, other, strict=True)) < MIN_HAMMING:
                    return False
        return True

    drawn = rejected_constraints = rejected_distance = rejected_hamming = 0
    while len(barcodes) < TARGET:
        drawn += 1
        bc = "".join(rng.choice("ACGT") for _ in range(LENGTH))
        if not ok(bc):
            rejected_constraints += 1
            continue
        nb = [
            _pack("".join(bc[i] for i in range(LENGTH) if i not in d))
            for d in itertools.combinations(range(LENGTH), k)
        ]
        if not seen.isdisjoint(nb):
            rejected_distance += 1
            continue
        if not hamming_ok(bc):
            rejected_hamming += 1
            continue
        seen.update(nb)
        for (bs, be), idx in zip(blocks, hidx, strict=True):
            idx.setdefault(bc[bs:be], []).append(bc)
        barcodes.append(bc)
        if len(barcodes) % 5000 == 0:
            print(f"  {len(barcodes):>6,} accepted   {drawn:>9,} drawn")

    print(f"\ndone: {len(barcodes):,} barcodes from {drawn:,} draws")
    print(
        f"  rejected by constraints: {rejected_constraints:>9,} "
        f"({100 * rejected_constraints / drawn:.2f}%)"
    )
    print(
        f"  rejected by Levenshtein: {rejected_distance:>9,} "
        f"({100 * rejected_distance / drawn:.2f}%)"
    )
    print(
        f"  rejected by Hamming:     {rejected_hamming:>9,} ({100 * rejected_hamming / drawn:.2f}%)"
    )

    gc = [(b.count("G") + b.count("C")) / LENGTH for b in barcodes]
    print(f"\n  GC: min {min(gc):.2f}  mean {sum(gc) / len(gc):.3f}  max {max(gc):.2f}")

    # every barcode must still pass every filter -- catches any bug in the accept path
    assert all(ok(b) for b in barcodes), "a written barcode fails the filter stack"
    assert len(set(barcodes)) == len(barcodes), "duplicate barcode"
    print("  all barcodes re-pass the full filter stack, and all are unique")

    # sampled independent distance check (full pairwise is 4.5e8 comparisons)
    n_pairs = 300_000
    worst = LENGTH
    for _ in range(n_pairs):
        a, b = rng.sample(barcodes, 2)
        worst = min(worst, levenshtein(a, b))
    print(
        f"  sampled {n_pairs:,} random pairs: minimum observed distance = {worst} "
        f"(guarantee is >= {MIN_DISTANCE})"
    )
    assert worst >= MIN_DISTANCE, "sampled pair violates the distance guarantee"

    with OUT.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["index", "barcode"])
        for i, b in enumerate(barcodes, 1):
            w.writerow([i, b])
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
