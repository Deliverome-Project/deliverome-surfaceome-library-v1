"""Put PGK on the ORF transcript in p0018/p0019, breaking the EF1a/EFS palindrome.

`EFS_v1` is the EF1a core promoter with the intron removed: 209 of its 212 nt sit
inside `EF1a_v1`. In p0018/p0019 the regulator cassette is antisense, so the two face
each other and form a **212 bp inverted repeat, 99.1% complementary, 25 nt loop**, with
an unbroken 150 bp perfectly paired run. Perfect palindromes much over 30 bp are
generally unstable in E. coli, so the expected outcome is a plasmid that comes back
deleted, and the same hairpin would sit in the genomic RNA during reverse transcription.

Swapping the ORF-transcript promoter to `PGK_v1` removes the homology outright, since
PGK is unrelated to EF1a. It is the right move for *this pair only*: the PGK arms
(p0016/p0017/p0020/p0021) already run PGK antisense as the regulator promoter, so
putting PGK on their ORF transcript too would recreate the same palindrome, larger.

**This changes what p0016-vs-p0018 measures.** The pair was built to vary the regulator
promoter with the ORF promoter held at EFS. After this swap the two differ in both, so
the comparison reads as "PGK regulator / EFS output" against "EF1a regulator / PGK
output" rather than isolating the regulator. Becca's call, 2026-10-06, on the grounds
that a vector that will not clone measures nothing at all.

    uv run --locked python swap_ef1a_output_promoter.py --dry-run
    uv run --locked python swap_ef1a_output_promoter.py
"""

from __future__ import annotations

import argparse
import tempfile
import warnings
from pathlib import Path

from Bio import BiopythonWarning, SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation

from deliverome_analysis import seqtools
from deliverome_analysis.s3 import processed

KEY = "constructs/designed/{}.gb"
DESIGNS = ("p0018", "p0019")
DONOR = "p0016"  # carries PGK_v1, already in the series
OLD_LABEL, NEW_LABEL = "EFS_v1", "PGK_v1"
KOZAK_JUNCTION = "GACCGGTTCTAGAGCCGCCACC"  # promoter -> ORF, must survive untouched
NOTE = (
    "ORF-transcript promoter swapped {date} from EFS_v1 to PGK_v1 ({old} -> {new} nt). "
    "EFS is the intron-less EF1a core, so 209 of its 212 nt lie inside the antisense "
    "EF1a_v1 regulator promoter and the two formed a 212 bp inverted repeat, 99.1% "
    "complementary with a 25 nt loop - a palindrome well past the ~30 bp E. coli "
    "tolerates. PGK shares no homology with EF1a. Note this pair now differs from the "
    "PGK pair in both promoters, not just the regulator. See "
    "s4e9-orf-design-fg-panel/swap_ef1a_output_promoter.py"
)
DATE = "2026-10-06"


def feature(rec, prefix: str):
    hits = [f for f in rec.features if f.qualifiers.get("label", [""])[0].startswith(prefix)]
    if len(hits) != 1:
        raise SystemExit(f"expected exactly one {prefix!r} feature, found {len(hits)}")
    return hits[0]


def shift(rec, cut_start: int, cut_end: int, delta: int) -> None:
    """Move every feature lying wholly after the replaced span by `delta`."""
    for f in rec.features:
        a, b = int(f.location.start), int(f.location.end)
        if a >= cut_end:
            f.location = FeatureLocation(a + delta, b + delta, f.location.strand)
        elif b > cut_start and not (a >= cut_start and b <= cut_end):
            # a feature straddling the boundary would silently mis-span afterwards
            if f.type != "source":
                raise SystemExit(
                    f"feature {f.qualifiers.get('label', [f.type])[0]!r} straddles the "
                    f"promoter being replaced ({a + 1}..{b}); refusing to guess"
                )
            f.location = FeatureLocation(0, b + delta, f.location.strand)


def longest_hairpin(seq: str, min_stem: int = 20, max_loop: int = 100) -> int:
    """Longest inverted-repeat stem with a short loop, allowing a few mismatches."""
    rc = seqtools.reverse_complement
    index: dict[str, list[int]] = {}
    for i in range(len(seq) - min_stem + 1):
        index.setdefault(seq[i : i + min_stem], []).append(i)
    best = 0
    for i0 in range(len(seq) - min_stem + 1):
        for j0 in index.get(rc(seq[i0 : i0 + min_stem]), ()):
            if not 0 <= j0 - (i0 + min_stem) <= max_loop:
                continue
            i, n, budget = i0, min_stem, 0
            while i - 1 >= 0 and j0 + n < len(seq):
                if seq[i - 1] == rc(seq[j0 + n]):
                    i, n = i - 1, n + 1
                elif budget < 4:
                    budget, i, n = budget + 1, i - 1, n + 1
                else:
                    break
            best = max(best, n)
    return best


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args(argv)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", BiopythonWarning)
        donor = SeqIO.read(processed.download(KEY.format(DONOR)), "genbank")
    pgk_feat = feature(donor, NEW_LABEL)
    pgk = str(pgk_feat.extract(donor.seq)).upper()
    pgk_label = pgk_feat.qualifiers.get("label", [NEW_LABEL])[0]
    print(f"donor {DONOR}: {pgk_label}  {len(pgk)} nt, taken in its functional orientation")

    staged = []
    for name in DESIGNS:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", BiopythonWarning)
            rec = SeqIO.read(processed.download(KEY.format(name)), "genbank")
        before = str(rec.seq).upper()
        if NEW_LABEL in [f.qualifiers.get("label", [""])[0].split(" [")[0] for f in rec.features]:
            print(f"{name}: already on {NEW_LABEL}, nothing to do")
            continue
        efs = feature(rec, OLD_LABEL)
        a, b = int(efs.location.start), int(efs.location.end)
        if efs.location.strand != 1:
            raise SystemExit(f"{name}: expected the output promoter on the plus strand")
        old_len = b - a
        delta = len(pgk) - old_len

        after = before[:a] + pgk + before[b:]
        shift(rec, a, b, delta)
        efs.location = FeatureLocation(a, a + len(pgk), 1)
        efs.qualifiers["label"] = [pgk_label]
        efs.qualifiers["note"] = [
            "replaces EFS_v1, which shares 209 of its 212 nt with the antisense "
            "EF1a_v1 regulator promoter and formed a 212 bp inverted repeat with it"
        ]
        rec.seq = Seq(after)
        rec.features.sort(key=lambda f: (int(f.location.start), -int(f.location.end)))

        # --- checks, all of which must hold before anything is written ---
        problems = []
        if before[:a] != after[:a] or before[b:] != after[a + len(pgk) :]:
            problems.append("sequence outside the promoter changed")
        if after.count(KOZAK_JUNCTION) != before.count(KOZAK_JUNCTION):
            problems.append("the promoter->ORF junction was disturbed")
        for enz, want in (("BbsI", 0), ("Esp3I", 0)):
            if len(seqtools.find_sites(after, [enz])) != want:
                problems.append(f"{enz} count is not {want}")
        orf = feature(rec, "TFRC")
        cds = str(orf.extract(Seq(after))).upper()
        if before.count(cds) != 1 or after.count(cds) != 1 or len(cds) % 3:
            problems.append("the ORF is not intact and in frame at its new coordinates")
        hp_before, hp_after = longest_hairpin(before), longest_hairpin(after)
        if hp_after >= 60:
            problems.append(f"a {hp_after} nt hairpin stem remains")
        for lab in ("4x17nt", "scramble sites", "O-0013", "10_bc_CS1", "20bp-BC"):
            got = [f for f in rec.features if f.qualifiers.get("label", [""])[0].startswith(lab)]
            for f in got:
                if str(f.extract(Seq(after))).upper() not in before:
                    problems.append(f"feature {lab!r} no longer extracts its own bases")
        if problems:
            raise SystemExit(f"{name}: refusing to write — {problems}")

        print(
            f"{name}  {len(before):,} -> {len(after):,} bp ({delta:+,})   "
            f"longest hairpin stem {hp_before} -> {hp_after} nt   "
            f"AscI {before.count('GGCGCGCC')}->{after.count('GGCGCGCC')}  "
            f"FseI {before.count('GGCCGGCC')}->{after.count('GGCCGGCC')}"
        )
        comment = (
            str(rec.annotations.get("comment", ""))
            + "\n"
            + NOTE.format(date=DATE, old=old_len, new=len(pgk))
        ).strip()
        if long := seqtools.overlong_genbank_words(comment):
            raise SystemExit(f"{name}: words too wide for GenBank: {long}")
        rec.annotations["comment"] = comment
        staged.append((name, rec, after))

    if args.dry_run:
        print("\ndry run — S3 not written")
        return 0

    for name, rec, after in staged:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / f"{name}.gb"
            SeqIO.write(rec, path, "genbank")
            processed.upload(path, KEY.format(name), content_type="chemical/seq-na-genbank")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", BiopythonWarning)
            back = SeqIO.read(processed.download(KEY.format(name)), "genbank")
        if str(back.seq).upper() != after:
            raise SystemExit(f"{name}: re-read does not match what was uploaded")
        if longest_hairpin(str(back.seq).upper()) >= 60:
            raise SystemExit(f"{name}: re-read still carries a long hairpin")
    print(f"\nuploaded {len(staged)} record(s); re-read and confirmed palindrome-free")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
