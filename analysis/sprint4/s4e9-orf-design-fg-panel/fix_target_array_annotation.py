"""Trim the flanking cloning sites out of the miR-L target-array feature.

The array feature in p0016, p0018 and p0020 spans 100 nt, but only the middle 84 nt is
the target array. The outer 16 nt are our own cloning sites: ``GCGATCGC`` (AsiSI) at the
5' end and ``GTTTAAAC`` (PmeI) at the 3'. The label asserts the span is *verbatim* from
Du et al., and for those 16 nt it is not.

The ``4x17nt`` part of the label is correct and is left alone, and the four sites it
names are now annotated individually so the map shows what the label claims. Each repeat is 21 nt and
the array is 4 of them, but only 17 nt of each is complementary to the miR-L guide — the
17-mer ``TAGATAAGCATTATAAT`` appears in the annotated hairpin — so the sites really are
4 x 17 nt with a 4-nt spacer after each.

Annotation only: **no base changes**. The two flanks are re-annotated as their own
features rather than discarded.

    uv run --locked python fix_target_array_annotation.py --dry-run
    uv run --locked python fix_target_array_annotation.py
"""

from __future__ import annotations

import argparse
import tempfile
from pathlib import Path
from typing import NamedTuple

from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation, SeqFeature

from deliverome_analysis import seqtools
from deliverome_analysis.s3 import processed

DESIGNS = ("p0016", "p0017", "p0018", "p0019", "p0020", "p0021")
KEY = "constructs/designed/{}.gb"
FLANK5, FLANK3 = ("GCGATCGC", "AsiSI"), ("GTTTAAAC", "PmeI")
SPACER = "TACT"  # the 4 nt that precede each site, giving the 21-nt period


class Arm(NamedTuple):
    """One of the two site blocks. Both share the 4x(4+17) nt architecture."""

    label_prefix: str  # how the array feature is labelled on the record
    site: str  # the 17 nt that repeat
    site_label: str  # what to call each copy
    site_note: str


ARMS = (
    Arm(
        "4x17nt",
        "ATTATAATGCTTATCTA",  # reverse complement of the guide-matching 17-mer
        "miR-L target site",
        "17 nt complementary to the miR-L guide",
    ),
    Arm(
        "scramble",
        "TGTAAATACTATTCTAT",  # same bases as the target site, shuffled
        "scramble site",
        "17 nt with no complementarity to the miR-L guide: the target site's own bases "
        "shuffled, so base composition and the repeat architecture are matched",
    ),
)

NOTE = (
    "{what} feature re-cut {date} to the 84 nt the array actually is: the span "
    "previously ran 100 nt and swallowed the flanking AsiSI and PmeI cloning sites, "
    "which are ours and not part of the array. Those 16 nt are now their own features, "
    "and the four 17-nt sites inside the array are annotated individually. Annotation "
    "only, no base changes. See s4e9-orf-design-fg-panel/fix_target_array_annotation.py"
)
DATE = "2026-10-06"


def site_label(arm: Arm, n: int) -> str:
    """The exact label for copy `n`.

    Matching per-site features on a *prefix* is a trap: "scramble site" is a prefix of
    the array's own label "scramble sites [...]", so a prefix test finds the array,
    concludes the sites already exist, and silently skips annotating them.
    """
    return f"{arm.site_label} {n} of 4"


def site_features(rec, arm: Arm) -> list:
    wanted = {site_label(arm, n) for n in (1, 2, 3, 4)}
    return [f for f in rec.features if f.qualifiers.get("label", [""])[0] in wanted]


def array_feature(rec, arm: Arm):
    """The array itself. Matched on the label prefix rather than on a substring like
    "miR-L target", which would also catch the per-site features added below."""
    return next(
        (
            f
            for f in rec.features
            if f.qualifiers.get("label", [""])[0].startswith(arm.label_prefix)
        ),
        None,
    )


def arm_of(rec) -> Arm | None:
    """Which arm this record is. The two are mutually exclusive by construction."""
    found = [a for a in ARMS if array_feature(rec, a) is not None]
    if len(found) > 1:
        raise SystemExit("record carries both a target and a scramble array")
    return found[0] if found else None


def retrim(rec, arm: Arm) -> list[str]:
    """Apply whichever of the two edits is still missing, so a re-run is safe."""
    f = array_feature(rec, arm)
    if f is None:
        return []
    array = (SPACER + arm.site) * 4
    seq = str(f.extract(rec.seq)).upper()
    did = []
    if seq != array:
        if not (seq.startswith(FLANK5[0]) and seq.endswith(FLANK3[0]) and array in seq):
            raise SystemExit(f"unexpected array content, refusing to guess: {seq[:20]}...")
        did.append("flanks")
    s, e, strand = int(f.location.start), int(f.location.end), f.location.strand
    n5, n3 = (len(FLANK5[0]), len(FLANK3[0])) if "flanks" in did else (0, 0)
    f.location = FeatureLocation(s + n5, e - n3, strand)
    flank_spans = ((FLANK5, (s, s + n5)), (FLANK3, (e - n3, e))) if "flanks" in did else ()
    for (_site, name), span in flank_spans:
        rec.features.append(
            SeqFeature(
                FeatureLocation(*span, strand),
                type="misc_feature",
                qualifiers={"label": [f"{name} site, flanking the array"]},
            )
        )
    if not site_features(rec, arm):
        annotate_sites(rec, arm, s + n5, strand)
        did.append("sites")
    rec.features.sort(key=lambda x: int(x.location.start))
    return did


def annotate_sites(rec, arm: Arm, array_start: int, strand: int) -> int:
    """Mark the four 17-nt sites inside the array.

    Each 21-nt repeat is a 4-nt ``TACT`` spacer followed by a 17-nt site; only the site
    carries the arm's identity. The spacers are left unannotated - four 4-nt features
    would be clutter - but the period is recorded on each site's note.
    """
    array = (SPACER + arm.site) * 4
    arr = str(rec.seq)[array_start : array_start + len(array)].upper()
    site = arm.site
    offsets = [i for i in range(len(arr) - len(site) + 1) if arr[i : i + len(site)] == site]
    if len(offsets) != 4:
        raise SystemExit(f"expected 4 sites inside the array, found {len(offsets)}")
    for n, off in enumerate(offsets, 1):
        rec.features.append(
            SeqFeature(
                FeatureLocation(array_start + off, array_start + off + len(site), strand),
                type="misc_feature",
                qualifiers={
                    "label": [site_label(arm, n)],
                    "note": [
                        f"{arm.site_note}; preceded by a {len(SPACER)}-nt {SPACER} spacer, "
                        f"giving a {len(SPACER) + len(site)}-nt repeat period"
                    ],
                },
            )
        )
    return len(offsets)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args(argv)

    staged = []
    for name in DESIGNS:
        rec = SeqIO.read(processed.download(KEY.format(name)), "genbank")
        before = str(rec.seq).upper()
        arm = arm_of(rec)
        if arm is None:
            print(f"{name}: no site array found, nothing to do")
            continue
        did = retrim(rec, arm)
        if not did:
            print(f"{name}: array span and sites already correct, nothing to do")
            continue
        if str(rec.seq).upper() != before:
            raise SystemExit(f"{name}: the sequence changed — this edit must be annotation-only")
        array = (SPACER + arm.site) * 4
        af = array_feature(rec, arm)
        got = str(af.extract(rec.seq)).upper()
        if got != array:
            raise SystemExit(f"{name}: re-cut span is not the 84-nt array: {len(got)} nt")
        print(
            f"{name}  {arm.site_label}s   changed: {'+'.join(did)}   array {len(got)} nt at "
            f"{int(af.location.start) + 1:,}..{int(af.location.end):,}"
        )
        sites = site_features(rec, arm)
        if len(sites) != 4:
            raise SystemExit(f"{name}: expected 4 per-site features, found {len(sites)}")
        print(f"     {got[:42]}")
        print(f"     {got[42:]}")
        print(f"     + AsiSI / PmeI annotated separately, + {len(sites)} individual 17-nt sites")
        comment = (
            str(rec.annotations.get("comment", ""))
            + "\n"
            + NOTE.format(date=DATE, what=af.qualifiers["label"][0].split(" [")[0])
        ).strip()
        if long := seqtools.overlong_genbank_words(comment):
            raise SystemExit(f"{name}: GenBank would drop these words from the comment: {long}")
        rec.annotations["comment"] = comment
        staged.append((name, rec, before, arm))

    if args.dry_run:
        print("\ndry run — S3 not written")
        return 0

    for name, rec, before, arm in staged:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / f"{name}.gb"
            SeqIO.write(rec, path, "genbank")
            processed.upload(path, KEY.format(name), content_type="chemical/seq-na-genbank")
        back = SeqIO.read(processed.download(KEY.format(name)), "genbank")
        if str(back.seq).upper() != before:
            raise SystemExit(f"{name}: uploaded sequence does not match the original")
        if str(array_feature(back, arm).extract(back.seq)).upper() != (SPACER + arm.site) * 4:
            raise SystemExit(f"{name}: re-read array feature is not the 84-nt array")
        if len(site_features(back, arm)) != 4:
            raise SystemExit(f"{name}: re-read record does not carry 4 per-site features")
    print(f"\nuploaded {len(staged)} record(s); sequences unchanged, array spans re-read clean")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
