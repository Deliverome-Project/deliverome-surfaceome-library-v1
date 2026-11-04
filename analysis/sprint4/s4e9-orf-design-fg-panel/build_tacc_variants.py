"""Build the TACC-spacer variant of the p0016 miR-L target vector.

Du's 4x17nt array repeats `TACT` + the 17 nt complementary to the miR-L guide. Only
position 4 pairs with the guide, so it must not be A; T, C and G all work. T is the one
that also completes a conserved miR-374-5p 8mer (`TATTATAA`) across the spacer/site
junction, four times. Since we cannot size miR-374-5p's contribution in the target line
without measuring it, these variants let the experiment answer it instead of assuming:
same vectors, one base per repeat changed, four bases in the whole plasmid.

C over G: both drop every conserved+ seed match and both raise array GC 19.0 -> 23.8%,
but `TACG` introduces a CpG in each repeat (4 per array) and `TACC` introduces none.
CpG density is a silencing substrate in integrated lentiviral transgenes, which is the
same reason SFFV was rejected for this series, so C is the conservative pick.

This is an *additional* construct, not a replacement: p0016 stays as the
published-verbatim arm and the pair differs only in the spacer base.

Scoped to p0016 alone (Becca, 2026-10-06). Variants of p0018 and p0020 were built and
then dropped: the effect this tests is a threshold-region one that cannot move the
compensated setpoint, so one variant on the arm being run first answers it, and p0018
and p0020 are the two least-settled designs in the series.

    uv run --locked python build_tacc_variants.py --dry-run
    uv run --locked python build_tacc_variants.py
"""

from __future__ import annotations

import argparse
import tempfile
import warnings
from pathlib import Path

from Bio import BiopythonWarning, SeqIO
from Bio.Seq import Seq

from deliverome_analysis import seqtools
from deliverome_analysis.s3 import processed

KEY = "constructs/designed/{}.gb"
SITE = "ATTATAATGCTTATCTA"
OLD_SPACER, NEW_SPACER = "TACT", "TACC"
OLD_ARRAY, NEW_ARRAY = (OLD_SPACER + SITE) * 4, (NEW_SPACER + SITE) * 4
M374 = "TATTATA"
# parent -> new P-number. P-0039 is the highest registered, so these follow it.
VARIANTS = {"p0016": "p0040"}
NOTE = (
    "miR-L target-array spacer changed {date} from TACT to TACC, one base per 21-nt "
    "repeat, four in the plasmid. Derived from {parent}, which keeps Du's verbatim "
    "TACT array. TACT|ATTATA completes a conserved miR-374-5p 8mer (TATTATAA) at each "
    "spacer/site junction; the 17-mer alone carries no conserved seed match, so the "
    "site is made by the junction. Position 4 only has to avoid A to keep the guide "
    "from re-pairing, so C is as valid as T, drops all four miR-374 sites, and unlike "
    "G adds no CpG. See s4e9-orf-design-fg-panel/build_tacc_variants.py"
)
DATE = "2026-10-06"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args(argv)

    staged = []
    for parent, child in VARIANTS.items():
        # Checked here, not only in the write path: the child has since been given its
        # own barcode, so rebuilding it from the parent would silently put the parent's
        # barcode back. A dry run that claims it would proceed invites exactly that.
        if list(processed.list_keys(KEY.format(child))):
            print(f"{child}: already built from {parent}, nothing to do")
            continue
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", BiopythonWarning)
            rec = SeqIO.read(processed.download(KEY.format(parent)), "genbank")
        before = str(rec.seq).upper()
        if before.count(OLD_ARRAY) != 1:
            raise SystemExit(f"{parent}: expected one TACT array, found {before.count(OLD_ARRAY)}")
        after = before.replace(OLD_ARRAY, NEW_ARRAY)

        problems = []
        if len(after) != len(before):
            problems.append("length changed")
        if sum(1 for a, b in zip(before, after, strict=True) if a != b) != 4:
            problems.append("more than 4 bases changed")
        # the 17-mer also occurs in the miR-L hairpin arm, so the invariant is that
        # the count is unchanged, not that it equals the array's four copies
        if after.count(SITE) != before.count(SITE):
            problems.append("a miR-L site was destroyed")
        if after.count(NEW_ARRAY) != 1:
            problems.append("the re-spaced array is not present exactly once")
        if M374 in after:
            problems.append("a miR-374 site survives")
        for enz, want in (("BbsI", 0), ("Esp3I", 0)):
            if len(seqtools.find_sites(after, [enz])) != want:
                problems.append(f"{enz} is no longer {want}")
        for motif in ("AATAAA", "ATTAAA", "CCCGGG", "TTAATTAA"):
            if after.count(motif) != before.count(motif):
                problems.append(f"{motif} count changed")
        if after.count("GCGATCGC") != 1 or after.count("GTTTAAAC") != 1:
            problems.append("AsiSI/PmeI flanks disturbed")
        if problems:
            raise SystemExit(f"{parent}->{child}: refusing to write — {problems}")

        # retarget the array features onto the new bases
        for f in rec.features:
            lab = f.qualifiers.get("label", [""])[0]
            if lab.startswith("4x17nt"):
                f.qualifiers["label"] = [lab.replace("4x17nt", "4x17nt TACC-spacer", 1)]
            if lab.startswith("miR-L target site"):
                f.qualifiers["note"] = [
                    "17 nt complementary to the miR-L guide; preceded by a 4-nt TACC "
                    "spacer (Du uses TACT, which completes a miR-374-5p 8mer here)"
                ]
        rec.seq = Seq(after)
        rec.id = rec.name = child
        rec.description = f"{rec.description} | TACC spacer variant of {parent}"
        comment = (
            str(rec.annotations.get("comment", "")) + "\n" + NOTE.format(date=DATE, parent=parent)
        ).strip()
        if long := seqtools.overlong_genbank_words(comment):
            raise SystemExit(f"{child}: words too wide for GenBank: {long}")
        rec.annotations["comment"] = comment
        print(
            f"{parent} -> {child}  {len(after):,} bp  4 base(s) changed  "
            f"miR-374 sites {before.count(M374)} -> {after.count(M374)}  "
            f"array GC {seqtools.gc_content(OLD_ARRAY) * 100:.1f} -> "
            f"{seqtools.gc_content(NEW_ARRAY) * 100:.1f}%  CpG added 0"
        )
        staged.append((child, rec, after))

    if args.dry_run:
        print("\ndry run — S3 not written")
        return 0

    for child, rec, after in staged:
        if list(processed.list_keys(KEY.format(child))):
            raise SystemExit(f"{child}: a record appeared since the check; refusing to clobber")
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / f"{child}.gb"
            SeqIO.write(rec, path, "genbank")
            processed.upload(path, KEY.format(child), content_type="chemical/seq-na-genbank")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", BiopythonWarning)
            back = SeqIO.read(processed.download(KEY.format(child)), "genbank")
        if str(back.seq).upper() != after:
            raise SystemExit(f"{child}: re-read does not match what was uploaded")
    print(f"\nuploaded {len(staged)} new record(s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
