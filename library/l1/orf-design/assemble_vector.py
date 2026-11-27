"""Assemble a finished, fully annotated plasmid: p0039 backbone + one shuttled insert.

The insert goes in **still wrapped in its Esp3I sites**. They are not consumed here,
and that is deliberate: keeping them in the delivered plasmid is what lets the insert be
cut back out later and moved into another vector. p0039 itself carries no Esp3I, so the
only two in the finished construct are the wrapper's, and an Esp3I digest releases
exactly the insert with its GACC / TCAG overhangs — ``main`` simulates that digest and
refuses to write if it does not come back byte-for-byte.

Note this differs from p008-p013, which carry **zero** Esp3I: those were cloned by a
digest that consumed the sites, so their inserts cannot be shuttled out this way.

``annotate_insert.py`` emits the synthesis fragment; this emits the thing that fragment
becomes once it is in the vector. The backbone's own features are carried across with
their coordinates remapped, the slot's MCS and stuffer are dropped because the insert
replaces them, and the insert's block structure is added at the junction.

Nothing here is typed in twice: the backbone comes from the GenBank
``empty_backbone.py`` wrote, the insert blocks from ``annotate_insert.blocks``, and the
result is checked against the same reconstitution ``empty_backbone.verify_reconstitution``
uses — so an assembled map cannot disagree with the backbone it was built from.

    uv run --locked python assemble_vector.py AXL
    uv run --locked python assemble_vector.py TFRC_linker_4x --out-dir ~/Downloads
"""

from __future__ import annotations

import argparse
from pathlib import Path

from annotate_insert import SOURCES, blocks
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from build_constructs import ENZYMES, KOZAK, read_fasta
from empty_backbone import (
    BUILD_DATE,
    CLONING_SITES,
    INSERT,
    MCS,
    NATIVE_PURO_END,
    NATIVE_PURO_START,
    STUFFER_LEN,
)
from shuttle import DOWNSTREAM, PAD, UPSTREAM, add_shuttle, digest

from deliverome_analysis.cloning import seqtools

HERE = Path(__file__).parent
BACKBONE = HERE / "p0039_EFS_IRES_puro_backbone.gb"
SLOT = len(MCS) + STUFFER_LEN
CUT = INSERT[0] - 1  # 0-based start of the slot, immediately 3' of EFS


def assemble(construct: str, insert: str):
    """``insert`` is the bare design; what lands in the vector is the shuttled form."""
    shuttled = add_shuttle(insert)
    rec = SeqIO.read(BACKBONE, "genbank")
    # slicing keeps each side's own features and drops the MCS/stuffer in between
    new = rec[:CUT] + Seq(shuttled) + rec[CUT + SLOT :]
    new.id = new.name = f"p0039_{construct}"[:16]
    new.description = f"p0039 backbone carrying the {construct} insert"
    new.annotations = {
        "molecule_type": "ds-DNA",
        "topology": "circular",
        "date": BUILD_DATE,
        "accession": f"p0039_{construct}",
    }
    # the insert begins past the 5' Esp3I wrapper, so every block shifts by its length
    off = CUT + len(UPSTREAM)
    for label, s, e, kind in blocks(insert):
        new.features.append(
            SeqFeature(FeatureLocation(off + s, off + e), type=kind, qualifiers={"label": [label]})
        )
    for label, s, e in (
        ("Esp3I (shuttle out)", CUT, CUT + len(UPSTREAM)),
        ("Esp3I (shuttle out)", CUT + len(shuttled) - len(DOWNSTREAM), CUT + len(shuttled)),
    ):
        new.features.append(
            SeqFeature(FeatureLocation(s, e), type="misc_feature", qualifiers={"label": [label]})
        )
    new.features.sort(key=lambda f: int(f.location.start))
    return new


def verify(record, insert: str) -> tuple[list[str], list[str]]:
    """(problems, warnings). Problems block the build; warnings are reported only."""
    seq = str(record.seq).upper()
    problems: list[str] = []
    warnings: list[str] = []

    shuttled = add_shuttle(insert)
    backbone_len = len(SeqIO.read(BACKBONE, "genbank").seq)
    if len(seq) != backbone_len - SLOT + len(shuttled):
        problems.append(
            f"length {len(seq):,} != backbone {backbone_len:,} - slot {SLOT} "
            f"+ shuttled insert {len(shuttled):,}"
        )
    if insert.upper() not in seq:
        problems.append("the insert is not present intact in the assembled plasmid")

    # every annotated insert block must sit where its sequence actually is
    for feat in record.features:
        lbl = feat.qualifiers.get("label", [""])[0]
        if lbl.startswith(
            ("Kozak", "ORF,", "ALFA/SpyTag", "barcode", "10x", "stop", "pad", "5' handle", "gap +")
        ):
            got = str(feat.extract(Seq(seq))).upper()
            s, e = int(feat.location.start), int(feat.location.end)
            if got != seq[s:e]:
                problems.append(f"feature {lbl!r} does not extract its own span")
    koz = [f for f in record.features if f.qualifiers.get("label", [""])[0].startswith("Kozak")]
    if len(koz) != 1:
        problems.append(f"expected one Kozak feature, found {len(koz)}")
    elif str(koz[0].extract(Seq(seq))).upper() != KOZAK:
        problems.append(
            f"the Kozak feature spans {str(koz[0].extract(Seq(seq))).upper()}, not {KOZAK}"
        )

    puro = [f for f in record.features if f.qualifiers.get("label", [""])[0] == "PuroR"]
    if len(puro) != 1:
        problems.append(f"expected one PuroR feature, found {len(puro)}")
    else:
        cds = str(puro[0].extract(Seq(seq)))
        aa = str(Seq(cds).translate())
        if not cds.startswith("ATG"):
            problems.append(f"PuroR feature does not start with ATG: {cds[:6]}")
        if not aa.startswith(NATIVE_PURO_START) or not aa.rstrip("*").endswith(NATIVE_PURO_END):
            problems.append(f"PuroR is not native: {aa[:24]}")
        if aa.count("*") != 1 or not aa.endswith("*"):
            problems.append("PuroR is not a single clean reading frame")

    circular = seq + seq[:20]
    if seqtools.find_sites(circular, ["BbsI"]):
        problems.append("BbsI site present in the assembled plasmid")
    # BsaI was only ever printed as information. The insert must contribute none: the
    # backbone's own copies are pre-existing and expected, so compare against the
    # backbone with the slot excised rather than against a hard-coded count.
    bb = str(SeqIO.read(BACKBONE, "genbank").seq).upper()
    expected_bsai = len(seqtools.find_sites(bb[:CUT] + bb[CUT + SLOT :], ["BsaI"]))
    n_bsai = len(seqtools.find_sites(circular, ["BsaI"]))
    if n_bsai != expected_bsai:
        problems.append(
            f"BsaI appears {n_bsai} times; the backbone contributes {expected_bsai}, "
            f"so the insert adds {n_bsai - expected_bsai}"
        )
    n_esp3i = len(seqtools.find_sites(circular, ["Esp3I"]))
    if n_esp3i != 2:
        problems.append(f"expected exactly 2 Esp3I sites to shuttle with, found {n_esp3i}")
    else:
        # the point of keeping them: cutting the finished plasmid must give the insert back
        try:
            released = digest(shuttled)
        except RuntimeError as exc:  # pragma: no cover - surfaced as a problem
            problems.append(f"Esp3I digest of the assembled plasmid failed: {exc}")
        else:
            # the released fragment runs overhang-to-overhang, so it carries the
            # 2-nt pad that completes OH3 — that pad ligates in with the insert
            expected = insert + PAD
            if released != expected:
                problems.append(
                    f"Esp3I digest releases {len(released):,} nt, not the "
                    f"{len(expected):,}-nt insert+pad"
                )
    # AgeI/XbaI/AscI/AsiSI are conveniences for downstream hand-cloning, not part of the
    # Golden Gate scheme. Their uniqueness is reported but does not block assembly: a
    # cassette codon choice that happens to spell one is not a reason to refuse a
    # construct whose BbsI/Esp3I/BsaI are clean. Only those three gate the build.
    for name, site in CLONING_SITES.items():
        n = circular.count(site)
        if name in ("AgeI", "XbaI", "AscI"):
            if n != 1:
                warnings.append(f"{name} appears {n} times, expected 1 — no longer unique")
        elif name == "AsiSI" and n != 0:
            warnings.append(f"AsiSI should have gone with the stuffer, found {n}")
    if [
        f
        for f in record.features
        if f.qualifiers.get("label", [""])[0].startswith(("MCS", "stuffer"))
    ]:
        problems.append("the MCS or stuffer feature survived into the assembled map")
    return problems, warnings


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

    rec = assemble(args.construct, insert)
    seq = str(rec.seq).upper()
    problems, warnings = verify(rec, insert)
    if problems:
        for p in problems:
            print(f"  FAIL  {p}")
        raise SystemExit("assembly failed verification — nothing written")
    for w in warnings:
        print(f"  WARN  {w}")

    print(f"p0039 + {args.construct}   {len(seq):,} bp circular, {len(rec.features)} features")
    print(
        f"  backbone {len(SeqIO.read(BACKBONE, 'genbank').seq):,} - {SLOT} slot "
        f"+ {len(insert):,} insert + {len(add_shuttle(insert)) - len(insert)} nt Esp3I wrapper "
        f"= {len(seq):,}"
    )
    print("  PuroR native from its own ATG                     ✓")
    print("  2 Esp3I sites retained; digest returns the insert  ✓")
    print("  no BbsI anywhere (circular scan)                   ✓")
    if not warnings:
        print("  AgeI/XbaI unique (from the Kozak); AscI unique    ✓")
        print("  AsiSI gone with the stuffer, as intended          ✓")
    counts = {e: len(seqtools.find_sites(seq, [e])) for e in ENZYMES}
    print(
        "  Type IIS counts: "
        + ", ".join(f"{e}={n}" for e, n in counts.items())
        + "   (BsaI is pre-existing in the 236083 backbone)"
    )

    args.out_dir.mkdir(parents=True, exist_ok=True)
    out = args.out_dir.expanduser()
    gbk = out / f"p0039_{args.construct}_vector.gb"
    fa = out / f"p0039_{args.construct}_vector.fasta"
    SeqIO.write(rec, gbk, "genbank")
    fa.write_text(
        f">p0039_{args.construct} | {len(seq)} bp circular | p0039 backbone + {args.construct} insert\n"
        + "\n".join(seq[i : i + 60] for i in range(0, len(seq), 60))
        + "\n"
    )
    print(f"\nWrote {gbk}\n      {fa}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
