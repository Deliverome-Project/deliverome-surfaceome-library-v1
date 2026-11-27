"""Full GenBank maps of one L1 member, AXL, in the two vector layouts Figure 2a compares.

Worked examples for reviewers, not registered constructs: they carry no P-number and the
barcode is the 20-N placeholder every L1 insert has until barcodes are assigned.

- ``l1_example_AXL_ires_puro.gb`` -- the library layout: the AXL L1 insert in the p0039
  backbone, ``EFS - insert - IRES2 - PuroR - WPRE``, one transcript.
- ``l1_example_AXL_sv40_puro.gb`` -- the alternative not used: the same insert with PuroR from
  its own SV40 promoter and no IRES, the P-0050 layout (P-0008's selection unit).

Nothing is re-implemented here:

- the insert is the AXL member's record (``l1-p0168`` in ``constructs_l1_nterm.fasta``),
  already in the new L1 readout layout (``build_constructs.l1_readout``: stop, CC, T7 promoter,
  5' handle, barcode, Feldman padlock construct + RT site with T->A at 9, 10x CS1);
- the IRES-puro map is ``assemble_vector.assemble`` and must pass ``assemble_vector.verify``;
- the SV40-puro map is ``build_p0050_set.build`` (P-0039 construct 5' of AscI, p0050 from AscI
  onward, as GenScript builds the P-0050 set) and must pass ``build_p0050_set.verify``.

This script only refines the insert's labels: the Feldman flank is split into its padlock
construct and RT site, the handle is named for its primer (JT832), the barcode is marked as a
placeholder, and the tagged ORF gets one translated CDS feature. Every added feature is checked
to extract the constant it names, and the translated ORF must hash to the member's recorded
``protein_sha256``.

    uv run --locked python l1_example_vector_maps.py            # build, verify, write
    uv run --locked python l1_example_vector_maps.py --upload   # also publish to constructs/scratch/
"""

from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import itertools
import tempfile
from pathlib import Path

import assemble_vector
import build_p0050_set
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from build_constructs import (
    CS1,
    FELDMAN_3PRIME,
    FWD_PRIMER,
    KOZAK,
    T7_PROMOTER,
    read_fasta,
)
from shuttle import DOWNSTREAM, UPSTREAM

from deliverome_analysis.cloning import constructs, seqtools

HERE = Path(__file__).parent
MEMBER, GENE = "l1-p0168", "AXL"
SOURCE_FASTA = HERE / "constructs_l1_nterm.fasta"
SOURCE_TSV = HERE / "constructs_l1_nterm.tsv"
#: Feldman's 17-nt padlock construct, then his RT site with T->A at position 9.
PADLOCK, RT_SITE = FELDMAN_3PRIME[:17], FELDMAN_3PRIME[17:]
BUILD_DATE = "22-NOV-2026"  # L1 readout adopted (build_constructs.l1_readout)
OUT = {
    "ires": HERE / f"l1_example_{GENE}_ires_puro.gb",
    "sv40": HERE / f"l1_example_{GENE}_sv40_puro.gb",
}
#: Descriptive names, never a P-number: these are worked examples, not registered constructs.
KEYS = {
    "ires": "constructs/scratch/l1-example-axl-ires-puro.gb",
    "sv40": "constructs/scratch/l1-example-axl-sv40-puro.gb",
}
DESCRIPTION = {
    "ires": (
        "Worked example, not a registered construct: L1 AXL insert (l1-p0168, barcode 20 N "
        "placeholder) in the p0039 library layout, EFS - insert - IRES2 - PuroR - WPRE"
    ),
    "sv40": (
        "Worked example, not a registered construct: L1 AXL insert (l1-p0168, barcode 20 N "
        "placeholder) in the P-0050 layout, PuroR from its own SV40 promoter, no IRES"
    ),
}


def label(f) -> str:
    return f.qualifiers.get("label", [""])[0]


def member_insert() -> tuple[str, str]:
    """The AXL member's bare insert (wrapper removed) and its recorded protein hash."""
    wrapped = read_fasta(SOURCE_FASTA)[MEMBER]
    assert wrapped.startswith(UPSTREAM + KOZAK) and wrapped.endswith(CS1 + DOWNSTREAM), MEMBER
    rows = [
        r for r in csv.DictReader(SOURCE_TSV.open(), delimiter="\t") if r["candidate_id"] == MEMBER
    ]
    assert len(rows) == 1 and rows[0]["gene_symbol"] == GENE, rows
    return wrapped[len(UPSTREAM) : -len(DOWNSTREAM)], rows[0]["protein_sha256"]


def refine(rec, protein_sha: str):
    """Split and rename the insert's readout labels; add one translated CDS for the ORF."""
    seq = str(rec.seq).upper()
    feats = []
    for f in rec.features:
        lab = label(f)
        s, e = int(f.location.start), int(f.location.end)
        if lab == "Feldman padlock + RT site":
            assert seq[s:e] == FELDMAN_3PRIME
            m = s + len(PADLOCK)
            feats += [
                SeqFeature(FeatureLocation(s, m, 1), type="misc_feature",
                           qualifiers={"label": ["padlock construct (Feldman, 17 nt)"]}),
                SeqFeature(FeatureLocation(m, e, 1), type="primer_bind",
                           qualifiers={"label": ["RT site (Feldman, T9A)"],
                                       "note": ["Feldman's RT site with T->A at position 9"]}),
            ]  # fmt: skip
            continue
        if lab == "5' handle":
            assert seq[s:e] == FWD_PRIMER
            f.qualifiers["label"] = ["5' handle (JT832)"]
        elif lab == "barcode":
            assert seq[s:e] == "N" * 20
            f.qualifiers["label"] = ["barcode (20 N placeholder)"]
        elif lab == "T7 promoter" and seq[s:e] == T7_PROMOTER:  # the insert's, not the backbone's
            f.qualifiers["label"] = ["T7 promoter (readout)"]
        feats.append(f)
    rec.features = feats

    koz = next(f for f in rec.features if label(f).startswith("Kozak"))
    stop = next(f for f in rec.features if label(f) == "stop")
    a, b = int(koz.location.end), int(stop.location.end)
    protein = str(Seq(seq[a:b]).translate())
    assert protein.endswith("*") and protein.count("*") == 1, "ORF is not one clean frame"
    got = hashlib.sha256(protein.rstrip("*").encode()).hexdigest()
    if got != protein_sha:
        raise SystemExit(f"{MEMBER}: translated ORF hashes to {got}, record says {protein_sha}")
    rec.features.append(
        SeqFeature(
            FeatureLocation(a, b, 1),
            type="CDS",
            qualifiers={
                "label": [f"{GENE} (tagged, N after signal peptide)"],
                "translation": [protein.rstrip("*")],
                "note": [f"L1 member {MEMBER}; protein_sha256 {protein_sha}"],
            },
        )
    )
    rec.features.sort(key=lambda f: int(f.location.start))
    return protein.rstrip("*")


def merge_puro(rec) -> None:
    """p0050 annotates PuroR's ATG apart from its body; give it one CDS so viewers translate it."""
    parts = sorted((f for f in rec.features if label(f) == "PuroR" and f.type == "CDS"),
                   key=lambda f: int(f.location.start))  # fmt: skip
    if len(parts) < 2:
        return
    s, e = int(parts[0].location.start), int(parts[-1].location.end)
    rec.features = [f for f in rec.features if f not in parts]
    rec.features.append(SeqFeature(FeatureLocation(s, e, 1), type="CDS",
                                   qualifiers={"label": ["PuroR"]}))  # fmt: skip
    rec.features.sort(key=lambda f: int(f.location.start))


def check(rec, protein: str, n_esp3i_expected: int = 2) -> list[str]:
    """Checks shared by both maps; returns human-readable lines, raises on failure."""
    seq = str(rec.seq).upper()
    lines = []
    want = {
        "T7 promoter (readout)": T7_PROMOTER,
        "5' handle (JT832)": FWD_PRIMER,
        "barcode (20 N placeholder)": "N" * 20,
        "padlock construct (Feldman, 17 nt)": PADLOCK,
        "RT site (Feldman, T9A)": RT_SITE,
        "10x Capture Sequence 1": CS1,
        "Kozak (AgeI-XbaI)": KOZAK,
    }
    for lab, s in want.items():
        hits = [f for f in rec.features if label(f) == lab]
        assert len(hits) == 1, f"{lab}: {len(hits)} features"
        assert str(hits[0].extract(rec.seq)).upper() == s, f"{lab} does not extract its constant"
    # the readout runs in this order, abutting, from the stop codon to CS1
    order = ["stop", "pad", *list(want)[:6]]
    spans = [next(f for f in rec.features if label(f) == lab).location for lab in order]
    assert all(int(x.end) == int(y.start) for x, y in itertools.pairwise(spans)), (
        "readout not contiguous"
    )
    lines.append(
        "readout: stop | CC | T7 | 5' handle (JT832) | barcode | padlock | RT site (T9A) | CS1, abutting  ok"
    )

    cds = next(f for f in rec.features if label(f).startswith(f"{GENE} (tagged"))
    assert str(cds.extract(rec.seq).translate()).rstrip("*") == protein, "protein changed"
    lines.append(f"{GENE} protein {len(protein)} aa, identical to the member record  ok")

    esp = [f for f in rec.features if label(f) == "Esp3I (shuttle out)"]
    circ = seq + seq[:20]
    n = len(seqtools.find_sites(circ, ["Esp3I"]))
    assert n == n_esp3i_expected and len(esp) == 2, f"Esp3I sites {n}, annotated {len(esp)}"
    assert str(esp[0].extract(rec.seq)).upper().startswith("CGTCTC")
    assert str(esp[1].extract(rec.seq)).upper().endswith("GAGACG")
    lines.append(
        "Esp3I: exactly 2 sites (circular scan), both the insert's inward-facing wrapper  ok"
    )
    assert not seqtools.find_sites(circ, ["BbsI"]), "BbsI present"
    lines.append("BbsI: none  ok")
    for name in ("PuroR",):
        p = [f for f in rec.features if label(f) == name and f.type == "CDS"]
        assert len(p) == 1, f"{name}: {len(p)} CDS features"
        aa = str(p[0].extract(rec.seq).translate())
        assert aa.startswith(build_p0050_set.PURO_START) and aa.rstrip("*").endswith(
            build_p0050_set.PURO_END
        ), "PuroR not native"
    lines.append("PuroR native, one CDS  ok")
    return lines


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument(
        "--upload", action="store_true", help="publish each map to processed://constructs/scratch/"
    )
    args = ap.parse_args(argv)

    insert, protein_sha = member_insert()

    # the committed backbone assemble_vector uses must be the registered p0039 map
    bb_local = str(SeqIO.read(assemble_vector.BACKBONE, "genbank").seq).upper()
    bb_s3 = str(SeqIO.read(constructs.fetch("p0039"), "genbank").seq).upper()
    assert bb_local == bb_s3, (
        "committed p0039 backbone differs from processed://" + constructs.map_key("p0039")
    )

    # (a) library layout: assemble_vector, verified by its own checks
    ires = assemble_vector.assemble(f"{GENE}_L1", insert)
    problems, warnings = assemble_vector.verify(ires, insert)
    if problems:
        raise SystemExit("assemble_vector.verify failed:\n  " + "\n  ".join(problems))
    protein = refine(ires, protein_sha)
    ires.id = ires.name = "l1_AXL_ires"
    ires.description = DESCRIPTION["ires"]
    ires.annotations["accession"] = "l1_example_AXL_ires_puro"

    # (b) the alternative: the P-0050 builder on the finished library-layout construct
    with tempfile.TemporaryDirectory() as tmp:
        p0050 = SeqIO.read(build_p0050_set._fetch("p0050", Path(tmp)), "genbank")
    sv40 = build_p0050_set.build(copy.deepcopy(ires), p0050, "l1_AXL_sv40")
    build_p0050_set.verify(sv40, ires, p0050)
    merge_puro(sv40)
    sv40.description = DESCRIPTION["sv40"]
    sv40.annotations["accession"] = "l1_example_AXL_sv40_puro"

    for rec in (ires, sv40):
        rec.annotations["date"] = BUILD_DATE

    report = {}
    for kind, rec in (("ires", ires), ("sv40", sv40)):
        lines = check(rec, protein)
        labels = {label(f) for f in rec.features}
        need = ["EFS_v1", "WPRE", "5' LTR (truncated)", "3' LTR (Delta-U3)", "cPPT/CTS"]
        need += ["IRES2_v1"] if kind == "ires" else ["SV40 enhancer", "SV40 promoter"]
        missing = [n for n in need if n not in labels]
        assert not missing, f"{kind}: unannotated {missing}"
        if kind == "sv40":
            assert "IRES2_v1" not in labels
        lines.append("annotated: " + ", ".join(need))
        report[kind] = lines
        SeqIO.write(rec, OUT[kind], "genbank")

    n_ires, n_sv40 = len(ires.seq), len(sv40.seq)
    assert n_sv40 == n_ires - 320, (n_ires, n_sv40)  # 620-bp IRES2 out, 300-bp SV40 in
    for kind, rec in (("ires", ires), ("sv40", sv40)):
        print(f"{OUT[kind].name}: {len(rec.seq):,} bp circular, {len(rec.features)} features")
        for line in report[kind]:
            print(f"  {line}")
    print(
        f"  insert {len(insert):,} nt bare, {len(UPSTREAM) + len(insert) + len(DOWNSTREAM):,} nt with the Esp3I wrapper"
    )
    print(f"  assemble_vector.verify and build_p0050_set.verify passed; {len(warnings)} warning(s)")
    for w in warnings:
        print(f"  WARN {w}")

    if args.upload:
        from deliverome_analysis.core import s3

        for kind, key in KEYS.items():
            if s3.processed.exists(key):
                raise SystemExit(f"processed://{key} already exists; not overwriting")
            s3.processed.upload(OUT[kind], key, content_type="text/plain")
            print(f"  uploaded processed://{key}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
