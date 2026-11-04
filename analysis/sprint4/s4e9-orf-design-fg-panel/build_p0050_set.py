"""Build the P-0050-backbone copies of the order batch: P-0051 .. P-0075.

GenScript makes these by cloning the 363-bp "50 insert" (the P-0008 SV40 enhancer/promoter
plus the first 50 bp of PuroR) into each finished P-0039 construct, replacing IRES2 and the
PuroR 5' end -- five premium synthesis items, each cloned into up to five vectors -- 25 slots, which is
exactly the 24 constructs plus the codon-optimised SLC9A3 (p0049 -> P-0075) (Mike Cohen,
2026-10-21, quote U7593ZDXG0). So the product of each one is **the P-0039 construct with
everything from AscI onward taken from p0050**::

    p00xx (P-0039 bb)   …EFS | Esp3I·Kozak·ORF·cassette·stop·pad·handle·barcode·Esp3I | AscI · IRES2 · PuroR · WPRE…
    p0050               …EFS | AgeI·XbaI·AsiSI·stuffer                                 | AscI · SV40 · PuroR · WPRE…
    p00yy (this)        …EFS | Esp3I·Kozak·ORF·cassette·stop·pad·handle·barcode·Esp3I | AscI · SV40 · PuroR · WPRE…

That is exactly the ``p0024_in_p0050_EXAMPLE`` map sent to GenScript on 2026-10-21 and
returned unchanged as quote item 55: the 5' flank, insert and both Esp3I shuttle sites come
from the P-0039 construct, and the SV40-PuroR unit, WPRE, 3' LTR and plasmid backbone from
p0050. It is 320 bp shorter than the P-0039 version (620-bp IRES2 out, 300-bp SV40 in).

Every map is checked before it is written: AgeI and AscI unique; the sequence 5' of AscI
identical to the source construct; the sequence from AscI onward identical to p0050; the
Kozak-to-stop ORF translating to the same protein as the source; PuroR present as the native
``MTEYKP…RKPGA``; and 0 Esp3I sites outside the insert's own pair.

Inputs are the designed maps already on S3 (``processed://constructs/designed/p00xx.gb``),
falling back to a committed copy only when S3 is unreachable. Numbers are assigned in
plasmid order from P-0051, matching ``p0050_set.tsv`` (the table sent to GenScript) and the
Reagent Inventory sheet; the script refuses to run if the table and the computed assignment
disagree, so the numbers cannot drift from what was registered.

    uv run --locked python build_p0050_set.py --out-dir ~/Downloads/p0050-set
    uv run --locked python build_p0050_set.py --upload      # publish to processed://constructs/designed/
"""

from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import re
import sys
import tempfile
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

HERE = Path(__file__).parent
TABLE = HERE / "p0050_set.tsv"
BACKBONE = "p0050"
FIRST_NUMBER = 51
AGEI, ASCI, ESP3I = "ACCGGT", "GGCGCGCC", "CGTCTC"
PURO_START, PURO_END = "MTEYKPTVRLATRDDVPRAVR", "TWCMTRKPGA"
#: The P-0039 constructs on quote U7593ZDXG0, in plasmid order.
SOURCES = [
    "p0022",
    "p0023",
    "p0024",
    "p0025",
    "p0026",
    "p0027",
    "p0028",
    "p0029",
    "p0030",
    "p0031",
    "p0032",
    "p0033",
    "p0034",
    "p0035",
    "p0037",
    "p0038",
    "p0041",
    "p0042",
    "p0043",
    "p0044",
    "p0045",
    "p0046",
    "p0047",
    "p0048",
    "p0049",  # codon-optimised SLC9A3 as GenScript built it; annotated by annotate_p0049.py
]


def _count(seq: str, site: str) -> int:
    rc = str(Seq(site).reverse_complement())
    n = seq.count(site)
    return n if rc == site else n + seq.count(rc)


def _fetch(name: str, workdir: Path) -> Path:
    """Download ``constructs/designed/<name>.gb`` from S3; fall back to a committed copy."""
    key = f"constructs/designed/{name}.gb"
    dest = workdir / f"{name}.gb"
    try:
        from deliverome_analysis import s3

        s3.processed.download(key, dest)
        return dest
    except Exception as exc:  # any S3 failure falls back, loudly
        local = next(HERE.glob(f"{name}_*.gb"), None) or (HERE / f"{name}.gb")
        if not local.exists():
            raise FileNotFoundError(f"{key} unreachable ({exc}) and no committed copy") from exc
        print(f"  [warn] {key}: {exc}; using committed copy {local.name}", file=sys.stderr)
        return local


def _orf_protein(seq: str) -> str:
    a = seq.index(AGEI)
    start = seq.index("GCCGCCACCATG", a) + 9
    prot = str(Seq(seq[start:]).translate(to_stop=True))
    return prot


def build(source: SeqRecord, backbone: SeqRecord, new_name: str) -> SeqRecord:
    s, b = str(source.seq).upper(), str(backbone.seq).upper()
    cut_s, cut_b = s.index(ASCI), b.index(ASCI)
    out_seq = s[:cut_s] + b[cut_b:]
    rec = SeqRecord(Seq(out_seq), id=new_name, name=new_name)
    rec.annotations = {
        "molecule_type": "ds-DNA",
        "topology": "circular",
        "data_file_division": "SYN",
    }
    what = re.sub(r"\s+in the p0039 backbone", "", source.description.split(":", 1)[-1]).strip()
    rec.description = (
        f"{new_name}: {what} -- the {source.name} insert in the {BACKBONE} backbone "
        f"(EFS ORF; IRES2-PuroR replaced by the P-0008 SV40 enhancer/promoter -> PuroR)"
    )
    shift = cut_s - cut_b
    for f in source.features:
        if int(f.location.end) <= cut_s:
            rec.features.append(copy.deepcopy(f))
    for f in backbone.features:
        if int(f.location.start) >= cut_b:
            g = copy.deepcopy(f)
            g.location = FeatureLocation(
                int(f.location.start) + shift, int(f.location.end) + shift, f.location.strand
            )
            rec.features.append(g)
    rec.features.append(
        SeqFeature(
            FeatureLocation(cut_s, cut_s + len(ASCI)),
            type="misc_feature",
            qualifiers={"label": "AscI (P-0039 / p0050 junction)"},
        )
    )
    rec.features.sort(key=lambda f: int(f.location.start))
    return rec


def verify(rec: SeqRecord, source: SeqRecord, backbone: SeqRecord) -> None:
    s, b, o = str(source.seq).upper(), str(backbone.seq).upper(), str(rec.seq).upper()
    assert _count(o, AGEI) == 1 and _count(o, ASCI) == 1, "AgeI/AscI must be unique"
    cut = o.index(ASCI)
    assert o[:cut] == s[: s.index(ASCI)], "5' of AscI must equal the source construct"
    assert o[cut:] == b[b.index(ASCI) :], "AscI onward must equal p0050"
    assert len(o) == len(s) - 320, f"expected {len(s) - 320} bp, got {len(o)}"
    assert _orf_protein(o) == _orf_protein(s), "ORF protein changed"
    puro = o.index("ATG" + "ACCGAGTACAAGCCCACGGTG")
    prot = str(Seq(o[puro:]).translate(to_stop=True))
    assert prot.startswith(PURO_START) and prot.endswith(PURO_END), "PuroR not native"
    assert _count(o, ESP3I) == 2, "exactly the insert's own Esp3I pair"


def expected_table() -> list[dict[str, str]]:
    rows = []
    for i, src in enumerate(SOURCES):
        rows.append({"source": src, "plasmid": f"p{FIRST_NUMBER + i:04d}"})
    return rows


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--out-dir", type=Path, default=HERE / "p0050_set")
    ap.add_argument(
        "--upload", action="store_true", help="publish to processed://constructs/designed/"
    )
    args = ap.parse_args(argv)

    with TABLE.open() as fh:
        table = list(csv.DictReader(fh, delimiter="\t"))
    exp = expected_table()
    got = [{"source": r["source"], "plasmid": r["plasmid"]} for r in table]
    if got != exp:
        sys.exit(f"{TABLE.name} disagrees with the computed assignment; fix one before building")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix="p0050set-"))
    backbone = SeqIO.read(_fetch(BACKBONE, work), "genbank")
    manifest = []
    for row in table:
        source = SeqIO.read(_fetch(row["source"], work), "genbank")
        rec = build(source, backbone, row["plasmid"])
        verify(rec, source, backbone)
        out = args.out_dir / f"{row['plasmid']}.gb"
        SeqIO.write(rec, out, "genbank")
        sha = hashlib.sha256(out.read_bytes()).hexdigest()
        manifest.append((row["plasmid"], row["source"], len(rec.seq), sha))
        print(
            f"  {row['source']} -> {row['plasmid']}  {len(rec.seq):>6} bp  {row['genscript_item_name']}"
        )
        if args.upload:
            from deliverome_analysis import s3

            s3.processed.upload(
                out, f"constructs/designed/{row['plasmid']}.gb", content_type="text/plain"
            )
    with (args.out_dir / "MANIFEST.tsv").open("w") as fh:
        fh.write("plasmid\tsource\tlength_bp\tsha256\n")
        for m in manifest:
            fh.write("\t".join(map(str, m)) + "\n")
    print(
        f"{len(manifest)} maps written to {args.out_dir}" + (" and uploaded" if args.upload else "")
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
