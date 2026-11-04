"""Trim the Kozak out of the PuroR CDS features in p0016-p0021.

Every PuroR feature in the six miR-L designs starts 6 nt before its start codon, on the
``GCCACC`` Kozak. The sequence is right and the protein is right, but translating the
annotated feature from its first base yields ``ATM...`` instead of ``MTEY...`` — two
residues that do not exist. Any tool that trusts the CDS span reports the wrong protein.

This is annotation-only: **no base changes**. For each affected feature the biological
5' end moves 6 nt (these sit on the minus strand, so that is the plus-strand ``end``),
and the 6 nt removed are re-annotated as their own ``regulatory`` Kozak feature so the
information is kept rather than discarded.

It is the mirror of the p0039 defect: there the PuroR feature began 3 nt *after* its ATG
and excluded it; here it begins 6 nt *before* and swallows the Kozak.

    uv run --locked python fix_puror_annotation.py --dry-run
    uv run --locked python fix_puror_annotation.py
"""

from __future__ import annotations

import argparse
import tempfile
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature

from deliverome_analysis import seqtools
from deliverome_analysis.s3 import processed

DESIGNS = ("p0016", "p0017", "p0018", "p0019", "p0020", "p0021")
KEY = "constructs/designed/{}.gb"
KOZAK = "GCCACC"
NATIVE_START = "MTEYKPTVRL"
NOTE = (
    "PuroR CDS features re-cut {date} to begin at the start codon: each previously "
    "began 6 nt earlier, on the GCCACC Kozak, so translating the feature gave ATM... "
    "instead of MTEY.... Annotation only, no base changes; the Kozak is now its own "
    "feature. See s4e9-orf-design-fg-panel/fix_puror_annotation.py"
)
DATE = "2026-10-06"


def affected(rec) -> list[SeqFeature]:
    """PuroR features whose span opens on the Kozak rather than on the ATG."""
    out = []
    for f in rec.features:
        label = f.qualifiers.get("label", [""])[0]
        if not label.startswith("PuroR"):
            continue
        if str(f.extract(rec.seq)).upper().startswith(KOZAK + "ATG"):
            out.append(f)
    return out


def retrim(rec) -> tuple[list[SeqFeature], list[tuple[int, int]]]:
    """Move each affected feature's biological 5' end past the Kozak."""
    spans, changed = [], []
    for f in affected(rec):
        changed.append(f)
        s, e, strand = int(f.location.start), int(f.location.end), f.location.strand
        if strand == -1:
            kozak_span = (e - len(KOZAK), e)
            f.location = FeatureLocation(s, e - len(KOZAK), strand)
        else:
            kozak_span = (s, s + len(KOZAK))
            f.location = FeatureLocation(s + len(KOZAK), e, strand)
        spans.append(kozak_span)
    for s, e in dict.fromkeys(spans):  # one Kozak feature per distinct span
        rec.features.append(
            SeqFeature(
                FeatureLocation(s, e, -1 if rec.features else 1),
                type="regulatory",
                qualifiers={"label": ["Kozak (PuroR)"], "regulatory_class": ["other"]},
            )
        )
    rec.features.sort(key=lambda f: int(f.location.start))
    return changed, spans


def verify(rec, before_seq: str, changed: list[SeqFeature]) -> list[str]:
    """Only the features this script re-cut are checked. ``PuroR exon 2`` is a second
    exon and legitimately does not begin with a start codon."""
    bad = []
    if str(rec.seq).upper() != before_seq:
        bad.append("the sequence changed — this edit must be annotation-only")
    for f in changed:
        label = f.qualifiers.get("label", [""])[0]
        ex = str(f.extract(rec.seq)).upper()
        if not ex.startswith("ATG"):
            bad.append(f"{label!r} still does not begin with ATG: {ex[:9]}")
            continue
        aa = str(Seq(ex[: len(ex) // 3 * 3]).translate())
        if not aa.startswith(NATIVE_START[: len(aa)][:4]):
            bad.append(f"{label!r} translates {aa[:8]}, expected {NATIVE_START[:8]}")
    if not [f for f in rec.features if f.qualifiers.get("label", [""])[0] == "Kozak (PuroR)"]:
        bad.append("the Kozak was removed without being re-annotated")
    return bad


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args(argv)

    staged = []
    for name in DESIGNS:
        rec = SeqIO.read(processed.download(KEY.format(name)), "genbank")
        before = str(rec.seq).upper()
        hits = affected(rec)
        if not hits:
            print(f"{name}: PuroR features already start at the ATG, nothing to do")
            continue
        changed, spans = retrim(rec)
        problems = verify(rec, before, changed)
        if problems:
            raise SystemExit(f"{name}: refusing to write — {problems}")
        shown = [
            f"{f.qualifiers.get('label', [''])[0][:34]} -> "
            f"{str(f.extract(rec.seq))[:9].upper()} "
            f"({Seq(str(f.extract(rec.seq))[:9]).translate()!s})"
            for f in changed
        ]
        print(f"{name}  {len(changed)} feature(s) re-cut, Kozak re-annotated at {spans[0]}")
        for line in shown:
            print(f"     {line}")
        comment = (str(rec.annotations.get("comment", "")) + "\n" + NOTE.format(date=DATE)).strip()
        if long := seqtools.overlong_genbank_words(comment):
            raise SystemExit(f"{name}: GenBank would drop these words from the comment: {long}")
        rec.annotations["comment"] = comment
        staged.append((name, rec, before, changed))

    if args.dry_run:
        print("\ndry run — S3 not written")
        return 0

    for name, rec, before, changed in staged:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / f"{name}.gb"
            SeqIO.write(rec, path, "genbank")
            processed.upload(path, KEY.format(name), content_type="chemical/seq-na-genbank")
        back = SeqIO.read(processed.download(KEY.format(name)), "genbank")
        if str(back.seq).upper() != before:
            raise SystemExit(f"{name}: the uploaded sequence does not match the original")
        if verify(back, before, affected(back) or changed):
            raise SystemExit(f"{name}: re-read still fails verification")
    print(f"\nuploaded {len(staged)} record(s); sequences unchanged, CDS features re-read clean")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
