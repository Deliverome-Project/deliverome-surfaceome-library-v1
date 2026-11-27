"""Annotate p0049 (GenScript's codon-optimised SLC9A3) with p0025's features, and publish it.

p0049 is p0025 with only the two SLC9A3 coding segments recoded by GenScript (quote
U7593ZDXG0 items 50/51, map received 2026-10-21 and unchanged in the 2026-10-24 revision).
GenScript's map annotates the whole insert as one ``Insert Sequence`` block, so it carries
no Kozak, ORF, cassette, barcode or 10x Capture Sequence features. Everything outside the two
recoded segments is byte-identical to p0025 and the two plasmids are the same length, so
p0025's features apply at the same coordinates; this copies them onto p0049's sequence.

Checked before anything is written:

* same length as p0025, and identical outside the two coding segments
  (bp 4838-4939 and 5093-7492);
* the Kozak-to-stop ORF translates to exactly p0025's protein;
* only the coding segments differ.

**Accepted as is (Becca, 2026-10-26).** The recoding introduced two BbsI sites and one BsaI
site inside the SLC9A3 ORF. p0049 exists only to compare codon-optimised against native
expression, so no Golden Gate step will use it, and the request (2026-10-22) for GenScript to
re-run the optimisation without them was dropped. The sites are annotated so nobody digests it by mistake. P-0075 (p0049 in the
p0050 backbone, ``build_p0050_set.py``) inherits them.

Input is GenScript's map as published on S3; the output replaces it at the same key, with
the sequence unchanged (the processed bucket is versioned, so the bare map stays readable).

    uv run --locked python annotate_p0049.py            # check and write locally
    uv run --locked python annotate_p0049.py --upload   # also publish to processed://
"""

from __future__ import annotations

import argparse
import copy
import re
import tempfile
from pathlib import Path

from Bio import SeqIO
from Bio.Restriction import BbsI, BsaI
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from deliverome_analysis.cloning import constructs

HERE = Path(__file__).parent
#: 1-based inclusive coding segments GenScript was asked to optimise (Becca to Mike, 2026-10-21).
RECODED = ((4838, 4939), (5093, 7492))
KOZAK_ATG = 4838


def _get(name: str, work: Path) -> SeqRecord:
    from deliverome_analysis.core import s3

    dest = work / f"{name}.gb"
    s3.processed.download(constructs.map_key(name), dest)
    return SeqIO.read(dest, "genbank")


def _protein(seq: str) -> str:
    return str(Seq(seq[KOZAK_ATG - 1 :]).translate(to_stop=True))


def build(p0049: SeqRecord, p0025: SeqRecord) -> tuple[SeqRecord, int]:
    a, b = str(p0049.seq).upper(), str(p0025.seq).upper()
    assert len(a) == len(b), f"length differs: {len(a)} vs {len(b)}"
    differing = [i + 1 for i in range(len(a)) if a[i] != b[i]]
    outside = [p for p in differing if not any(lo <= p <= hi for lo, hi in RECODED)]
    assert not outside, f"p0049 differs from p0025 outside the coding segments at {outside[:5]}"
    assert _protein(a) == _protein(b), "protein changed"

    rec = SeqRecord(p0049.seq, id="p0049", name="p0049")
    rec.annotations = {
        "molecule_type": "ds-DNA",
        "topology": "circular",
        "data_file_division": "SYN",
    }
    rec.description = (
        "p0049: SLC9A3, codon-optimised by GenScript, native SP 1-25, N-terminal SpyTag003-ALFA "
        "cassette after residue 34, in the p0039 backbone; otherwise identical to p0025"
    )
    for f in p0025.features:
        g = copy.deepcopy(f)
        label = g.qualifiers.get("label", [""])[0]
        if f.type == "CDS" and ("SLC9A3" in label or label.startswith("ORF")):
            g.qualifiers["note"] = [
                "codon-optimised by GenScript (quote U7593ZDXG0 item 50, map 2026-10-21); "
                "protein identical to p0025"
            ]
        rec.features.append(g)
    for enzyme in (BbsI, BsaI):
        for strand_site in {enzyme.site, str(Seq(enzyme.site).reverse_complement())}:
            for m in re.finditer(f"(?={strand_site})", a):
                if not any(lo <= m.start() + 1 <= hi for lo, hi in RECODED):
                    continue
                rec.features.append(
                    SeqFeature(
                        FeatureLocation(m.start(), m.start() + len(strand_site)),
                        type="misc_feature",
                        qualifiers={
                            "label": [f"{enzyme} site (introduced by codon optimisation)"],
                            "note": ["accepted 2026-10-26: p0049 is not used in Golden Gate"],
                        },
                    )
                )
    rec.features.sort(key=lambda f: int(f.location.start))
    return rec, len(differing)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument(
        "--out", type=Path, default=Path(tempfile.mkdtemp(prefix="p0049-")) / "p0049.gb"
    )
    ap.add_argument("--upload", action="store_true")
    args = ap.parse_args(argv)
    work = Path(tempfile.mkdtemp(prefix="p0049-in-"))
    p0049, p0025 = _get("p0049", work), _get("p0025", work)
    rec, n_diff = build(p0049, p0025)
    assert str(rec.seq).upper() == str(p0049.seq).upper(), "sequence must not change"
    SeqIO.write(rec, args.out, "genbank")
    added = [
        f.qualifiers["label"][0]
        for f in rec.features
        if "introduced" in f.qualifiers.get("label", [""])[0]
    ]
    print(
        f"p0049: {len(rec.seq)} bp, {n_diff} nt recoded, {len(rec.features)} features; flagged: {added}"
    )
    if args.upload:
        from deliverome_analysis.core import s3

        s3.processed.upload(args.out, constructs.map_key("p0049"), content_type="text/plain")
        print("uploaded processed://" + constructs.map_key("p0049"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
