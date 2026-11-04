"""Move the six miR-L designs onto the new amplicon handle, so O-0013 binds them.

p0016-p0021 still carry the original ``JT_fw_XmaI_mod`` forward handle. That handle is a
poor primer in its own right — a 10-bp palindrome ``GACCCGGGTC``, self-dimer
-13.7 kcal/mol, 74% GC — and it is the same XmaI palindrome that blocked O-0009/O-0010
in s4e11. The s4e9 constructs already moved to ``AGTCGTTAGCTGCCGTTCC`` (O-0013), which
has no palindrome and a Tm of 67.5 C.

Both handles are 19 nt, so this is a straight substitution: no coordinate shifts, and
every downstream annotation keeps its position.

**The 2-nt pad moves to ``CC`` too**, so the whole P-series shares one amplicon layout
(Becca, 2026-10-06). That is a deliberate trade: it takes the stop codon's +4 base from a
purine to a pyrimidine, ``TAAG`` -> ``TAAC``, which is the weaker termination context. TAA
is still the strongest stop codon and readthrough stays low, but these are C-terminally
tagged constructs, so any readthrough product carries an extension through the barcode.

    uv run --locked python retrofit_mirna_handle.py --dry-run
    uv run --locked python retrofit_mirna_handle.py
"""

from __future__ import annotations

import argparse
import re
import tempfile
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

from deliverome_analysis import seqtools
from deliverome_analysis.s3 import processed

DESIGNS = ("p0016", "p0017", "p0018", "p0019", "p0020", "p0021")
KEY = "constructs/designed/{}.gb"
OLD_HANDLE = "TGGACCCGGGTCAGGCTCC"
NEW_HANDLE = "AGTCGTTAGCTGCCGTTCC"  # O-0013, TFRC_bc_fwd_JT832
OLD_PAD, NEW_PAD = "GG", "CC"  # the 2 nt between the stop codon and the handle
WPRE_REV = "AAAGCCATACGGGAAGCAATAGC"  # the matched reverse primer
TYPE_IIS = ["BbsI", "Esp3I", "BsaI"]
SPLICE_DONOR = re.compile(r"GT[AG]AG")
NOTE = (
    "Forward amplicon handle swapped {date} from JT_fw_XmaI_mod "
    "(TGGACCCGGGTCAGGCTCC) to AGTCGTTAGCTGCCGTTCC (O-0013) and the 2-nt pad from "
    "GG to CC, matching the s4e9 constructs. Stop context becomes TAAC. "
    "See s4e9-orf-design-fg-panel/retrofit_mirna_handle.py"
)
DATE = "2026-10-06"


def check(before: str, after: str) -> list[str]:
    bad = []
    if len(before) != len(after):
        bad.append("length changed")
    if after.count(NEW_HANDLE) != 1:
        bad.append(f"new handle appears {after.count(NEW_HANDLE)} times, expected 1")
    if OLD_HANDLE in after:
        bad.append("old handle still present")
    for enz in TYPE_IIS:
        b, a = len(seqtools.find_sites(before, [enz])), len(seqtools.find_sites(after, [enz]))
        if a != b:
            bad.append(f"{enz} count changed {b} -> {a}")
    if after.count("CCCGGG") > before.count("CCCGGG"):
        bad.append("introduced an XmaI site")
    i = after.find(NEW_HANDLE)
    w_b, w_a = before[i - 30 : i + 50], after[i - 30 : i + 50]
    if w_a.count("AATAAA") != w_b.count("AATAAA"):
        bad.append("created or destroyed an AATAAA near the junction")
    if len(SPLICE_DONOR.findall(w_a)) != len(SPLICE_DONOR.findall(w_b)):
        bad.append("created or destroyed a splice donor near the junction")
    if seqtools.reverse_complement(WPRE_REV) not in after:
        bad.append("the matched reverse primer site is missing")
    return bad


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args(argv)

    staged = []
    for name in DESIGNS:
        rec = SeqIO.read(processed.download(KEY.format(name)), "genbank")
        before = str(rec.seq).upper()
        after, did = before, []
        if OLD_HANDLE in after:
            if after.count(OLD_HANDLE) != 1:
                raise SystemExit(
                    f"{name}: expected one old handle, found {after.count(OLD_HANDLE)}"
                )
            after = after.replace(OLD_HANDLE, NEW_HANDLE)
            did.append("handle")
        h = after.find(NEW_HANDLE)
        if h < 0:
            raise SystemExit(f"{name}: new handle not found after the swap")
        if after[h - len(OLD_PAD) : h] == OLD_PAD:
            after = after[: h - len(OLD_PAD)] + NEW_PAD + after[h:]
            did.append("pad")
        if not did:
            print(f"{name}: already on the new handle and pad, nothing to do")
            continue
        problems = check(before, after)
        if problems:
            raise SystemExit(f"{name}: refusing to write — {problems}")
        pos = after.find(NEW_HANDLE)
        print(
            f"{name}  {len(after):,} bp   changed: {'+'.join(did):<11} "
            f"XmaI {before.count('CCCGGG')}->{after.count('CCCGGG')}   "
            f"stop+4 {before[pos - len(NEW_PAD)]}->{after[pos - len(NEW_PAD)]}   "
            f"…{after[pos - 7 : pos]}[{NEW_HANDLE[:8]}…]…"
        )
        # the annotation labels the old oligo by name; retitle it
        for feat in rec.features:
            if feat.qualifiers.get("label", [""])[0].startswith("JT_fw"):
                feat.qualifiers["label"] = ["O-0013 TFRC_bc_fwd_JT832"]
        rec.seq = Seq(after)
        comment = (str(rec.annotations.get("comment", "")) + "\n" + NOTE.format(date=DATE)).strip()
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
            processed.upload(path, KEY.format(name), content_type="chemical/seq-na-genbank")
        back = str(SeqIO.read(processed.download(KEY.format(name)), "genbank").seq).upper()
        if back != after:
            raise SystemExit(f"{name}: re-read does not match what was uploaded")
    print(f"\nuploaded {len(staged)} record(s); each re-read and confirmed on the new handle")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
