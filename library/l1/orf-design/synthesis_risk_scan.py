"""Manufacturability scan over the constructs staged for the vendor.

Answers one question before an order goes out: which regions will a gene synthesis
house flag as complex, and which of those are ours to fix. Vendors quote and schedule
off exactly these features, so finding them here is cheaper than finding them in a
"this region failed synthesis" email two weeks later.

A flagged region is only news if it is *new*. The lentiviral backbone is full of
extreme-GC stretches -- PuroR and PGK both sit near 80% -- but those exact bases are in
plasmids GenScript has already delivered to us and in Addgene deposits that exist as
physical DNA, so they are evidence that the region synthesises, not a warning that it
will not. Each hit is therefore checked against those proven sequences and reported as
inherited or new. Only the new ones need a decision.

Thresholds are the common vendor "complex sequence" bands, stated explicitly so they
can be argued with rather than trusted:

* **GC window** - GC below 25% or above 75% across a 50-nt window. Extreme windows
  drive oligo misannealing during assembly.
* **Direct repeats** - a substring of >=40 nt occurring more than once. Long direct
  repeats are both an assembly problem and a recombination substrate in E. coli and
  in the provirus.
* **Homopolymer** - a run of >=9 identical bases, which causes slippage.
* **Inverted repeats** - a stem of >=20 nt with a loop under 100 nt, i.e. a hairpin
  stable enough to stall polymerase. Perfect palindromes over ~30 bp are generally
  unstable in E. coli, so stem length is the number to read here.

Every hit is reported with the annotated feature it falls in, because "a 19% GC window
at 9,094" is only actionable once you know it is the miR-L target array.

    uv run --locked python synthesis_risk_scan.py
    uv run --locked python synthesis_risk_scan.py --window 100
    uv run --locked python synthesis_risk_scan.py --key constructs/ordered/p008.gb
"""

from __future__ import annotations

import argparse
import warnings
from pathlib import Path

from Bio import BiopythonWarning, SeqIO

from deliverome_analysis import constructs, seqtools

DESIGNED = (
    "p0016",
    "p0017",
    "p0018",
    "p0019",
    "p0020",
    "p0021",
    "p0039",
    "p0040",
)
# Sequences that exist as physical DNA: plasmids we ordered and received, and Addgene
# deposits. Anything found verbatim in one of these has demonstrably been synthesised.
PROVEN_PROCESSED = ("constructs/ordered/", "constructs/reference/addgene/")
# ordered/ holds every plasmid we have ordered, not only those delivered. These were ordered
# from GenScript but not yet received (Reagent Inventory Status "Ordered", read 2026-11-09),
# so their sequence is not yet evidence that it synthesises. Drop a name once it arrives.
NOT_YET_RECEIVED = frozenset(
    {"p0022", "p0023", "p0024", "p0025", "p0029", "p0038", "p0041"}
    | {f"p{n:04d}" for n in range(51, 76)}
)
PROVEN_EXTERNAL = ("du-2025-dimmer-mirna-circuits/2026-10-06/",)
GC_WINDOW = 50
GC_LOW, GC_HIGH = 0.25, 0.75
MIN_DIRECT_REPEAT = 40
MIN_HOMOPOLYMER = 9
MIN_STEM, MAX_LOOP = 20, 100
MAX_STEM_MISMATCHES = 4


def proven_sequences() -> dict[str, str]:
    """Both strands of every sequence we know has been built, keyed by name."""
    from deliverome_analysis.s3 import external, processed

    out: dict[str, str] = {}
    for prefix in PROVEN_PROCESSED:
        for key in processed.list_keys(prefix):
            if not key.endswith((".gb", ".gbk")) or Path(key).stem in NOT_YET_RECEIVED:
                continue
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", BiopythonWarning)
                rec = SeqIO.read(processed.download(key), "genbank")
            out[Path(key).stem] = str(rec.seq).upper()
    for prefix in PROVEN_EXTERNAL:
        for key in external.list_keys(prefix):
            if not key.endswith((".gb", ".gbk")):
                continue
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", BiopythonWarning)
                rec = SeqIO.read(external.download(key), "genbank")
            out[Path(key).stem] = str(rec.seq).upper()
    return out


def built_before(sub: str, proven: dict[str, str]) -> str | None:
    """Name of a proven sequence containing `sub` on either strand, if any.

    Both strands and a doubled sequence, so a circular plasmid's origin-spanning
    regions and minus-strand cassettes still match.
    """
    rc = seqtools.reverse_complement(sub)
    for name, seq in proven.items():
        doubled = seq + seq
        if sub in doubled or rc in doubled:
            return name
    return None


def novel_fraction(sub: str, proven: dict[str, str]) -> float:
    """Fraction of `sub` not covered by the longest stretch found in proven DNA.

    A window straddling a junction between two proven cassettes is almost entirely
    sequence that has been synthesised before, and saying so keeps a junction
    artefact from reading like a novel problem region.
    """
    best = 0
    for seq in proven.values():
        doubled = seq + seq
        for strand in (doubled, seqtools.reverse_complement(doubled)):
            lo, hi = best, len(sub)
            while lo < hi:
                mid = (lo + hi + 1) // 2
                if any(sub[i : i + mid] in strand for i in range(len(sub) - mid + 1)):
                    lo = mid
                else:
                    hi = mid - 1
            best = max(best, lo)
    return 1 - best / len(sub) if sub else 0.0


def merge(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Collapse overlapping or touching spans, so one GC valley is one row."""
    out: list[list[int]] = []
    for s, e in sorted(intervals):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [(s, e) for s, e in out]


def features_at(rec, start: int, end: int) -> str:
    """The most specific annotated feature overlapping [start, end), for the report."""
    hits = [
        f
        for f in rec.features
        if f.type != "source" and int(f.location.start) < end and int(f.location.end) > start
    ]
    if not hits:
        return "unannotated"
    best = min(hits, key=lambda f: int(f.location.end) - int(f.location.start))
    return best.qualifiers.get("label", [best.type])[0].split(" [")[0]


def gc_windows(seq: str, window: int) -> list[tuple[int, int, float]]:
    """Merged spans whose GC leaves the band, with the most extreme GC in each."""
    bad = [
        (i, i + window, seqtools.gc_content(seq[i : i + window]))
        for i in range(0, len(seq) - window + 1)
        if not GC_LOW <= seqtools.gc_content(seq[i : i + window]) <= GC_HIGH
    ]
    out = []
    for s, e in merge([(s, e) for s, e, _ in bad]):
        inside = [gc for bs, be, gc in bad if bs >= s and be <= e]
        out.append((s, e, max(inside, key=lambda g: abs(g - 0.5))))
    return out


def direct_repeats(seq: str, minimum: int) -> list[tuple[int, int, int]]:
    """Maximal repeated substrings as (first position, second position, length).

    Only the maximal one per repeated region: a 63-nt repeat also contains a 62-nt
    repeat one base along, and reporting every suffix buries the single real finding
    under twenty restatements of it.
    """
    seen: dict[str, int] = {}
    hits: list[tuple[int, int, int]] = []
    for i in range(len(seq) - minimum + 1):
        sub = seq[i : i + minimum]
        if sub not in seen:
            seen[sub] = i
            continue
        j = seen[sub]
        n = minimum
        while i + n < len(seq) and seq[j + n] == seq[i + n]:
            n += 1
        hits.append((j, i, n))
    # drop any hit that is contained in an earlier, longer one at the same offset
    hits.sort(key=lambda h: -h[2])
    kept: list[tuple[int, int, int]] = []
    for j, i, n in hits:
        if any(
            kj <= j and j + n <= kj + kn and ki <= i and i + n <= ki + kn and (i - j) == (ki - kj)
            for kj, ki, kn in kept
        ):
            continue
        kept.append((j, i, n))
    return sorted(kept)


def homopolymers(seq: str, minimum: int) -> list[tuple[int, str, int]]:
    out, i = [], 0
    while i < len(seq):
        j = i
        while j + 1 < len(seq) and seq[j + 1] == seq[i]:
            j += 1
        if j - i + 1 >= minimum:
            out.append((i, seq[i], j - i + 1))
        i = j + 1
    return out


def hairpins(seq: str, min_stem: int, max_loop: int) -> list[tuple[int, int, int]]:
    """Maximal inverted repeats as (stem start, stem end, stem length).

    Extension tolerates isolated mismatches, because stopping at the first one
    badly under-reports the thing that matters. The EF1a/EFS stem in p0018 is
    212 bp at 99% complementarity, but it carries two mismatches 150 bp in, so a
    perfect-match-only walk calls it 59 bp -- under a quarter of its real size,
    and on the wrong side of every published stability threshold.
    """
    rc = seqtools.reverse_complement
    index: dict[str, list[int]] = {}
    for i in range(len(seq) - min_stem + 1):
        index.setdefault(seq[i : i + min_stem], []).append(i)
    hits: list[tuple[int, int, int]] = []
    for i0 in range(len(seq) - min_stem + 1):
        for j0 in index.get(rc(seq[i0 : i0 + min_stem]), ()):
            if not 0 <= j0 - (i0 + min_stem) <= max_loop:
                continue
            i, j, n, budget = i0, j0, min_stem, 0
            while i - 1 >= 0 and j + n < len(seq):
                if seq[i - 1] == rc(seq[j + n]):
                    i, n = i - 1, n + 1
                elif budget < MAX_STEM_MISMATCHES:
                    budget, i, n = budget + 1, i - 1, n + 1
                else:
                    break
            hits.append((i, j + n, n))
    out: list[tuple[int, int, int]] = []
    for s, e, n in sorted(hits, key=lambda h: -h[2]):
        if any(os_ <= s and e <= oe for os_, oe, _ in out):
            continue
        out.append((s, e, n))
    return sorted(out)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--window", type=int, default=GC_WINDOW)
    ap.add_argument(
        "--key",
        action="append",
        default=None,
        help="processed:// GenBank key to scan; repeatable. Defaults to the s4e9 batch.",
    )
    args = ap.parse_args(argv)

    from deliverome_analysis.s3 import processed

    keys = args.key or [constructs.map_key(n) for n in DESIGNED]
    paths = [(Path(k).stem, Path(processed.download(k))) for k in keys]

    proven = proven_sequences()
    print(f"proven references: {len(proven)} sequence(s) that exist as physical DNA")
    print("  " + ", ".join(sorted(proven)))

    new_total = inherited_total = 0
    for name, path in paths:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", BiopythonWarning)
            rec = SeqIO.read(path, "genbank")
        seq = str(rec.seq).upper()
        rows: list[tuple[str, str, str, int, int]] = []

        for s, e, gc in gc_windows(seq, args.window):
            rows.append((f"GC {gc * 100:.0f}%", f"{s + 1:,}..{e:,}", features_at(rec, s, e), s, e))
        for j, i, n in direct_repeats(seq, MIN_DIRECT_REPEAT):
            rows.append(
                (
                    f"repeat {n} nt",
                    f"{j + 1:,} and {i + 1:,}",
                    features_at(rec, j, j + n) + " / " + features_at(rec, i, i + n),
                    i,
                    i + n,
                )
            )
        for i, base, n in homopolymers(seq, MIN_HOMOPOLYMER):
            rows.append((f"poly{base} x{n}", f"{i + 1:,}", features_at(rec, i, i + n), i, i + n))
        for s, e, n in hairpins(seq, MIN_STEM, MAX_LOOP):
            rows.append(
                (
                    f"hairpin stem {n} nt",
                    f"{s + 1:,}..{e:,} (loop {e - s - 2 * n})",
                    features_at(rec, s, e),
                    s,
                    e,
                )
            )

        fresh = [(k, w, f, built_before(seq[s:e], proven), seq[s:e]) for k, w, f, s, e in rows]
        news = [r for r in fresh if r[3] is None]
        old_ = [r for r in fresh if r[3] is not None]
        new_total += len(news)
        inherited_total += len(old_)

        gc = seqtools.gc_content(seq) * 100
        print(f"\n{name}  {len(seq):,} bp  GC {gc:.1f}%")
        print(f"    {len(old_)} flag(s) inherited from DNA that already exists -- not a new risk")
        if not news:
            print("    no new complex region")
        for kind, where, feat, _, sub in news:
            nov = novel_fraction(sub, proven)
            tag = "junction of proven parts" if nov < 0.15 else f"{nov:.0%} novel sequence"
            print(f"    NEW  {kind:<18} {where:<26} {feat}")
            print(f"         {'':<18} {tag}")

    print(
        f"\n{new_total} new flag(s), {inherited_total} inherited, across {len(paths)} construct(s)"
    )
    print(
        f"bands: GC outside {GC_LOW:.0%}-{GC_HIGH:.0%} over {args.window} nt; "
        f"direct repeat >={MIN_DIRECT_REPEAT} nt; homopolymer >={MIN_HOMOPOLYMER}; "
        f"hairpin stem >={MIN_STEM} nt with loop <={MAX_LOOP} nt"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
