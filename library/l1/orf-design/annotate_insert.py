"""Emit an annotated GenBank for any designed insert, so the layout is readable in a map.

The construct FASTAs carry sequence only. This reconstructs the block structure from the
same constants the builders use — not from a hand-typed coordinate table — so an
annotation can never drift from what was actually built. Every block is located by
searching for its own sequence and the spans are asserted to tile the insert exactly.

    uv run --locked python annotate_insert.py AXL
    uv run --locked python annotate_insert.py TFRC_linker_4x --out-dir ~/Downloads
"""

from __future__ import annotations

import argparse
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord
from build_constructs import (
    CASSETTE_AA,
    CS1,
    ENZYMES,
    FELDMAN_3PRIME,
    FWD_PRIMER,
    KOZAK,
    L1_AFTER_BARCODE,
    L1_BEFORE_BARCODE,
    REV_PRIMER,
    SPACER,
    STOP,
    T7_PROMOTER,
    read_fasta,
)

from deliverome_analysis import plasmid_qc as qc
from deliverome_analysis import seqtools

HERE = Path(__file__).parent
#: The N-terminal entry is the disorder-anchored rebuild, NOT constructs_nterm.fasta.
#: That file is superseded: it places AXL's cassette inside AXL's own signal peptide, so a
#: map assembled from it loses the tag to cleavage. Every GenBank map in the p0022-p0029
#: batch was built from it, and all seven whose junction moved carry the old insert.
SOURCES = (
    "constructs_nterm_disorder_anchored.fasta",
    "constructs_loop.fasta",
    "constructs_tag_arms.fasta",
    # The GPCR rescue arms. Omitting this file is why the three recycled slots
    # (p0034 AVPR2, p0037 CNR1, p0038 CCR5) had no GenBank map and had never been
    # through plasmid_qc: the assembler could not see their inserts, so
    # rebuild_order_maps.py reported "unknown construct" and the arms went barcoded
    # but unchecked. An ordered construct no checker has seen is the gap these
    # modules exist to close.
    "constructs_gpcr_rescue.fasta",
    # The signal-peptide swap arms, for the same reason: an insert this list cannot
    # see is an insert plasmid_qc never checks, and p0048 went barcoded and
    # unpublished on its first rebuild until this line existed.
    "constructs_sp_swap.fasta",
)
BUILD_DATE = "06-OCT-2026"


# The 51-aa terminal cassette is not the only one this panel builds. The loop arms carry
# a 19-aa ALFA-only tag, because a 51-aa insert does not belong in a short extracellular
# loop, and the TFRC linker arms lengthen the inter-tag spacer to 3x or 4x GGGGS. Looking
# only for the canonical string made this function raise StopIteration on six of the
# fifteen inserts, so every variant the panel builds is enumerated here.
SPYTAG003 = qc.TAG_PEPTIDES["SpyTag003"]
ALFA_P = qc.TAG_PEPTIDES["P-ALFA-P"]
G4S = "GGGGS"
LOOP_TAG_AA = qc.LOOP_CASSETTE_PEPTIDE  # build_loop_inserts.LOOP_TAG_AA

CASSETTE_VARIANTS = [
    ("ALFA/SpyTag cassette", CASSETTE_AA),
    *(
        (f"ALFA/SpyTag cassette ({n}x GGGGS spacer)", G4S + SPYTAG003 + G4S * n + ALFA_P + G4S)
        for n in (3, 4)
    ),
    ("ALFA tag (loop)", LOOP_TAG_AA),
]
# The canonical cassette must decompose into the same parts the variants are built from,
# or the 3x/4x spacers would be assembled out of pieces that do not match what was built.
assert G4S + SPYTAG003 + G4S * 2 + ALFA_P + G4S == CASSETTE_AA


def blocks(insert: str) -> list[tuple[str, int, int, str]]:
    """(label, start, end, type) for every block, located by sequence, not by arithmetic."""
    # An L1 insert ends in the L1 readout region (T7 promoter and Feldman's 3' flank, Becca
    # 2026-11-22, build_constructs.l1_readout); every earlier insert in the legacy layout.
    l1 = insert.endswith(L1_AFTER_BARCODE)
    tail = (
        len(L1_BEFORE_BARCODE) + 20 + len(L1_AFTER_BARCODE)
        if l1
        else len(STOP) + len(SPACER) + len(FWD_PRIMER) + 20 + len(REV_PRIMER) + len(CS1)
    )
    cds_end = len(insert) - tail
    cds = insert[len(KOZAK) : cds_end]
    found = None
    for label, cass_aa in CASSETTE_VARIANTS:
        width = len(cass_aa) * 3
        at = next(
            (
                i
                for i in range(0, len(cds) - width + 1, 3)
                if seqtools.translate(cds[i : i + width]) == cass_aa
            ),
            None,
        )
        if at is not None:
            found = (label, at, width)
            break
    if found is None:
        raise SystemExit("no known tag cassette found — add its variant to CASSETTE_VARIANTS")
    cass_label, cass, cassette_nt = found
    o = len(KOZAK)
    out = [("Kozak (AgeI-XbaI)", 0, len(KOZAK), "regulatory")]
    if cass > 0:  # a cassette at offset 0 leaves no N-terminal ORF block to annotate
        out.append(("ORF, N-terminal", o, o + cass, "CDS"))
    out.append((cass_label, o + cass, o + cass + cassette_nt, "CDS"))
    rest = len(cds) - cass - cassette_nt
    if rest:
        out.append(("ORF, remainder", o + cass + cassette_nt, cds_end, "CDS"))
    p = cds_end
    barcode_at = p + len(L1_BEFORE_BARCODE if l1 else STOP + SPACER + FWD_PRIMER)
    barcode = ("barcode", insert[barcode_at:][:20], "misc_feature")
    layout = (
        (
            ("stop", STOP, "regulatory"),
            ("pad", SPACER, "misc_feature"),
            ("T7 promoter", T7_PROMOTER, "promoter"),
            ("5' handle", FWD_PRIMER, "primer_bind"),
            barcode,
            ("Feldman padlock + RT site", FELDMAN_3PRIME, "primer_bind"),
            ("10x Capture Sequence 1", CS1, "misc_feature"),
        )
        if l1
        else (
            ("stop", STOP, "regulatory"),
            ("pad", SPACER, "misc_feature"),
            ("5' handle", FWD_PRIMER, "primer_bind"),
            barcode,
            ("gap + RT site", REV_PRIMER, "primer_bind"),
            ("10x Capture Sequence 1", CS1, "misc_feature"),
        )
    )
    for label, seq, kind in layout:
        out.append((label, p, p + len(seq), kind))
        p += len(seq)
    if p != len(insert):
        raise RuntimeError(f"blocks cover {p} nt but the insert is {len(insert)}")
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("construct")
    ap.add_argument("--out-dir", type=Path, default=HERE)
    args = ap.parse_args(argv)

    inserts: dict[str, str] = {}
    for name in SOURCES:
        if (HERE / name).exists():
            inserts.update(read_fasta(HERE / name))
    if args.construct not in inserts:
        raise SystemExit(f"unknown construct {args.construct!r}; have: {sorted(inserts)}")
    insert = inserts[args.construct]

    rec = SeqRecord(
        Seq(insert),
        id=args.construct,
        name=args.construct[:16],
        description=f"s4e9 synthesis insert, {len(insert)} bp linear",
        annotations={"molecule_type": "ds-DNA", "topology": "linear", "date": BUILD_DATE},
    )
    for label, s, e, kind in blocks(insert):
        rec.features.append(
            SeqFeature(FeatureLocation(s, e), type=kind, qualifiers={"label": [label]})
        )
    sites = seqtools.find_sites(insert, ENZYMES)
    if sites:
        raise RuntimeError(f"{args.construct}: Type IIS site(s) present: {sites}")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    gbk = args.out_dir.expanduser() / f"{args.construct}_insert.gb"
    fa = args.out_dir.expanduser() / f"{args.construct}_insert.fasta"
    SeqIO.write(rec, gbk, "genbank")
    fa.write_text(
        f">{args.construct} s4e9 insert | {len(insert)} bp\n"
        + "\n".join(insert[i : i + 60] for i in range(0, len(insert), 60))
        + "\n"
    )
    print(f"{args.construct}  {len(insert):,} bp, {len(rec.features)} features, 0 Type IIS")
    for label, s, e, _ in blocks(insert):
        print(f"  {s + 1:>6,}..{e:<6,} {e - s:>5} nt  {label}")
    print(f"\nWrote {gbk} and {fa}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
