"""Remove the inherited BbsI site from the six miR-L iFFL designs (p0016-p0021).

All six carry one BbsI site that came in with the bGH poly(A) from Addgene 235285 — it
is not a designed assembly junction, it was simply never domesticated. It is the only
thing stopping BbsI being usable as a second-tier Golden Gate enzyme across the whole
plasmid collection: p0039 is already BbsI-clean, and these six are not.

**Where the site sits.** The poly(A) is annotated on the **minus** strand, as part of the
antisense regulator cassette, so the functional reading runs the other way: the AATAAA
hexamer is 92 nt *upstream* of the BbsI site on the strand that matters, which puts the
site past both the cleavage position (~10-30 nt after the hexamer) and the GU/U-rich
downstream element. It is in the element's 3' tail, not outside it — low risk, but not
zero, and worth sequencing over if poly(A) readthrough ever shows up.

**Which base.** All 18 single-base changes that destroy the site are equally clean: none
creates a Golden Gate site, a spurious AATAAA, a GT[AG]AG splice donor or a homopolymer
run, and none gains any commercially available enzyme. The tie-break is therefore
conservatism, so the chosen edit is the transition ``G>A``, which also happens to be the
first base of the recognition sequence.

    uv run --locked python domesticate_mirna_bbsi.py --dry-run
    uv run --locked python domesticate_mirna_bbsi.py
"""

from __future__ import annotations

import argparse
import re
import tempfile
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

from deliverome_analysis import constructs, seqtools
from deliverome_analysis.s3 import processed

DESIGNS = ("p0016", "p0017", "p0018", "p0019", "p0020", "p0021")
NEW_BASE = "A"  # the transition; see the module docstring
SPLICE_DONOR = re.compile(r"GT[AG]AG")
NOTE = (
    "BbsI site in the bGH poly(A) domesticated {date}: single G>A transition at {pos}, "
    "92 nt downstream of the AATAAA hexamer on the functional (minus) strand. The site "
    "was inherited from Addgene 235285 and was not a designed junction. "
    "See s4e9-orf-design-fg-panel/domesticate_mirna_bbsi.py"
)
DATE = "2026-10-06"


def fix(name: str):
    rec = SeqIO.read(processed.download(constructs.map_key(name)), "genbank")
    seq = str(rec.seq).upper()
    sites = seqtools.find_sites(seq, ["BbsI"])
    if not sites:
        return rec, seq, seq, None
    if len(sites) != 1:
        raise SystemExit(f"{name}: expected one BbsI site, found {len(sites)}")
    pos = sites[0][2]
    mutated = seq[:pos] + NEW_BASE + seq[pos + 1 :]
    return rec, seq, mutated, pos


def check(name: str, before: str, after: str, pos: int) -> list[str]:
    bad = []
    if len(before) != len(after):
        bad.append("length changed")
    if sum(1 for a, b in zip(before, after, strict=True) if a != b) != 1:
        bad.append("more than one base changed")
    if seqtools.find_sites(after, ["BbsI"]):
        bad.append("BbsI still present")
    for enz in ("Esp3I", "BsaI"):
        if len(seqtools.find_sites(after, [enz])) != len(seqtools.find_sites(before, [enz])):
            bad.append(f"{enz} count changed")
    w0, w1 = before[pos - 40 : pos + 46], after[pos - 40 : pos + 46]
    if w1.count("AATAAA") != w0.count("AATAAA"):
        bad.append("created or destroyed an AATAAA")
    if len(SPLICE_DONOR.findall(w1)) != len(SPLICE_DONOR.findall(w0)):
        bad.append("created or destroyed a splice donor")
    if seqtools.has_homopolymer_run(w1[30:56], 4):
        bad.append("created a homopolymer run")
    return bad


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args(argv)

    staged = []
    for name in DESIGNS:
        rec, before, after, pos = fix(name)
        if pos is None:
            print(f"{name}: already BbsI-free, nothing to do")
            continue
        problems = check(name, before, after, pos)
        if problems:
            raise SystemExit(f"{name}: refusing to write — {problems}")
        print(
            f"{name}  {len(before):,} bp   BbsI at {pos + 1:,}  "
            f"{before[pos]}>{after[pos]}   …{before[pos - 8 : pos]}"
            f"[{before[pos : pos + 6]}->{after[pos : pos + 6]}]{before[pos + 6 : pos + 14]}…"
        )
        rec.seq = Seq(after)
        comment = (
            str(rec.annotations.get("comment", "")) + "\n" + NOTE.format(date=DATE, pos=pos + 1)
        ).strip()
        if long := seqtools.overlong_genbank_words(comment):
            raise SystemExit(f"{name}: GenBank would drop these words from the comment: {long}")
        rec.annotations["comment"] = comment
        staged.append((name, rec, after))

    if args.dry_run:
        print("\ndry run — S3 not written")
        return 0

    for name, rec, after in staged:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / f"{name}.gb"
            SeqIO.write(rec, path, "genbank")
            processed.upload(path, constructs.map_key(name), content_type="chemical/seq-na-genbank")
        back = str(SeqIO.read(processed.download(constructs.map_key(name)), "genbank").seq).upper()
        if back != after:
            raise SystemExit(f"{name}: re-read does not match what was uploaded")
    print(f"\nuploaded {len(staged)} record(s); each re-read and confirmed BbsI-free")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
