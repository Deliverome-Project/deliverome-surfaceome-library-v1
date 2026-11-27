"""Correct processed://constructs/ordered/p004.gb to the sequence GenScript actually built.

The record on S3 was a design map, not the as-built plasmid. GenScript's full-sequence
export for order item U431LEXFG0-2 differs from it at **11 bases**, and those 11 bases
change three cloning sites:

====================  ===================  ====================
region                map said             built is
====================  ===================  ====================
4,589..4,596          ``TTAATTAA`` PacI    ``CCTGCAGG`` SbfI
4,764 / 4,769         ``ATCGAT``  ClaI     ``CTCGAG``  XhoI
6,093                 ``CCCGGG``  XmaI     ``CCGCGG``  SacII
====================  ===================  ====================

So the plasmid is not ``paci_clai_site`` as its name says — those sites are SbfI and XhoI
in the tube. Type IIS content is unaffected (Esp3I/BbsI/BsaI are 0/0/3 either way), so
none of this touches Golden Gate; they are restriction-cloning edits Becca made to the
Addgene 236083 parent.

Same length and same rotation, and no annotated feature overlaps a changed base, so all
50 features transfer unchanged and only the sequence is replaced.

    uv run --locked python fix_p004_record.py --dry-run
    uv run --locked python fix_p004_record.py
"""

from __future__ import annotations

import argparse
import tempfile
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

from deliverome_analysis.core.s3 import processed

KEY = "constructs/ordered/p004.gb"
DEFAULT_AS_BUILT = (
    Path.home() / "Downloads" / "U431LEXFG0-2 full sequence.seq"
)  # allow-local-path: the vendor export lands in Downloads and has no S3 copy yet; --as-built overrides
SOURCE_NOTE = (
    "Sequence corrected 2026-10-06 to GenScript's as-built full-sequence export for "
    "order item U431LEXFG0-2. The previous record was a design map and differed at 11 "
    "bases: PacI->SbfI (4589..4596), ClaI->XhoI (4764/4769), XmaI->SacII (6093). "
    "Feature set carried over unchanged; no feature overlapped a changed base."
)
EXPECTED_DIFFS = 11


def load_as_built(path: Path) -> str:
    if not path.exists():
        raise SystemExit(f"as-built export not found: {path}")
    return "".join(c for c in path.read_text().upper() if c in "ACGT")


def build(as_built: Path = DEFAULT_AS_BUILT):
    rec = SeqIO.read(processed.download(KEY), "genbank")
    old = str(rec.seq).upper()
    built = load_as_built(as_built)

    if len(built) != len(old):
        raise SystemExit(f"length mismatch: record {len(old):,} vs as-built {len(built):,}")
    diffs = [i for i in range(len(old)) if old[i] != built[i]]
    if len(diffs) != EXPECTED_DIFFS:
        raise SystemExit(
            f"expected {EXPECTED_DIFFS} differences, found {len(diffs)} — the as-built file "
            "may be a different construct or a different rotation; stopping rather than "
            "overwriting the record"
        )
    touched = [
        f.qualifiers.get("label", ["?"])[0]
        for f in rec.features
        if f.type != "source"
        and any(int(f.location.start) <= d < int(f.location.end) for d in diffs)
    ]
    if touched:
        raise SystemExit(
            f"these features overlap a changed base and would need re-annotating: {touched}"
        )

    rec.seq = Seq(built)
    rec.annotations.setdefault("comment", "")
    rec.annotations["comment"] = (rec.annotations["comment"] + "\n" + SOURCE_NOTE).strip()
    return rec, old, built, diffs


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dry-run", action="store_true", help="show the change; write nothing")
    ap.add_argument(
        "--as-built",
        type=Path,
        default=DEFAULT_AS_BUILT,
        help="GenScript full-sequence export to correct the record against",
    )
    args = ap.parse_args(argv)

    rec, old, built, diffs = build(args.as_built)
    print(f"{KEY}  {len(built):,} bp, {len(rec.features)} features")
    print(f"  {len(diffs)} base(s) corrected:")
    for d in diffs:
        print(
            f"    {d + 1:>6,}  {old[d]} -> {built[d]}    …{old[max(0, d - 10) : d]}[{old[d]}/{built[d]}]{old[d + 1 : d + 11]}…"
        )

    if args.dry_run:
        print("\ndry run — S3 not written")
        return 0

    with tempfile.TemporaryDirectory() as tmp:
        staged = Path(tmp) / "p004.gb"
        SeqIO.write(rec, staged, "genbank")
        processed.upload(staged, KEY, content_type="chemical/seq-na-genbank")

    check = str(SeqIO.read(processed.download(KEY), "genbank").seq).upper()
    if check != built:
        raise SystemExit("re-read of the uploaded record does not match the as-built sequence")
    print(f"\nuploaded, and re-read confirms {len(check):,} bp matching the as-built export")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
