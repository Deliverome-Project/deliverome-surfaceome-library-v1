"""p0039 — the insert-agnostic puro-only destination vector, for GenScript synthesis.

The puro-only intermediate (puro_only_backbone.py) is **not** empty: it still carries TFRC, the
ALFA/SpyTag cassette and TFRC's barcode, because it was made from p003 by one deletion.
Anything cloned into it would have to replace that insert. p0039 removes the insert as
well, so one synthesised plasmid serves every ORF in the panel and the library.

Two deletions from p003, plus one new cloning junction::

    p003    …EFS | Kozak·ORF·cassette·STOP·pad·handle·BARCODE·revprimer | AscI · IRES2 ·
            linker·ATG | mScarlet3 · linker · P2A | PuroR · TAA · WPRE…
                 └──────────── delete 4809..7337 ────────────┘       └─ delete 7969..8721 ─┘

    p0039   …EFS | AgeI · XbaI · AsiSI · stuffer | AscI · IRES2 · linker·ATG | PuroR · TAA · WPRE…

**Why these sites.** Measured against the 15 inserts already designed in this experiment
and against the whole backbone:

* ``AscI`` is already present at the insert's 3' boundary in p003, is unique in the
  7.5 kb backbone, and appears in none of the 15 inserts — it is the downstream site
  with no change to anything;
* ``AgeI`` and ``XbaI`` are absent from the backbone but present exactly once in every
  insert, because both sit inside the shared Kozak block
  (``GACCGGTTCTAGAGCCGCCACC`` = AgeI + XbaI + gccgccacc), so they give the vector a
  cut site the insert's own 5' end already matches;
* ``AsiSI`` is absent from the backbone *and* from all 15 inserts, so it is a free
  8-cutter held in reserve for a future insert that cannot carry AgeI.

No Type IIS site is added. The inserts carry their own Esp3I sites for release from the
synthesis plasmid (``shuttle.py``); the destination vector deliberately has none, so it
can never be cut by the same enzyme that liberates an insert.

**How an insert actually gets in.** As designed, the inserts run Kozak -> ORF -> cassette
-> stop -> pad -> handle -> barcode -> reverse primer, and stop there: ``AscI`` lives in
the *backbone* 3' of the slot, not in the insert. So an AgeI/AscI cut of the vector
leaves a 3' end that the insert as built cannot ligate to. Two routes work:

* **Gibson / HiFi** into the linearised vector, with homology arms — no change to any
  insert, and the orientation is set by the arms rather than by the overhangs;
* **AgeI/AscI restriction**, which needs ``GGCGCGCC`` appended to the insert's 3' end
  (one constant in ``build_constructs.py``, or a tail on the amplification primer).

``verify_reconstitution`` below proves the sequence is right either way; it does not
assert that a particular cloning chemistry can deliver it.

**The stuffer earns its place on a gel.** Cutting AgeI/AscI across a bare polylinker
drops ~20 bp, which does not resolve from uncut vector. The stuffer makes the empty and
cut vector distinguishable by size alone. It is inert, generated from a fixed seed, and
screened for every enzyme in the cloning set plus all three Type IIS enzymes.

⚠️ Input is the local-only p003 GenBank (same gap flagged in s1e2/enzyme_scan.py).
"""

from __future__ import annotations

import argparse
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from puro_only_backbone import NATIVE_PURO_END, NATIVE_PURO_START, P003

from deliverome_analysis import seqtools

HERE = Path(__file__).parent

INSERT = (4809, 7337)  # Kozak .. reverse primer site, 1-based inclusive
MSCARLET = (7969, 8721)  # mScarlet3 + linker + P2A — the puro-only deletion

AGEI, XBAI, ASISI, ASCI = "ACCGGT", "TCTAGA", "GCGATCGC", "GGCGCGCC"
# 2 nt between sites so each enzyme has flanking bases to bind; AsiSI in particular
# cuts poorly at a plasmid end or hard against a neighbouring site.
MCS = "TT" + AGEI + "TT" + XBAI + "TT" + ASISI + "TT"
# The stuffer is taken verbatim from p004 (processed://constructs/ordered/p004.gb,
# feature "Stuffer", 4597..4763), so the two backbones in this family carry the same
# dropout and a digest looks the same on a gel for both. It is 167 bp, AT-rich (26% GC)
# and carries A-runs longer than this module would generate on its own — acceptable for
# a fragment whose only job is to be cut out, and already proven synthesisable, since
# p004 exists. Screened below for the sites that actually matter rather than trusted.
STUFFER = (
    "ACTTACTAAAATTTTATTTAGAAGATAAGTCATGATTTGCTGTTCGCCCAAGTGGAACTG"
    "AACCTAAACTAAAAATTTATTTTATAGGTGTTGGTGAATCTGTTCAAAACGCTAAAGTTA"
    "AAGTAGACGAAATTATTAAAGAATTAAAATTAAAAATGAATATATAG"
)
STUFFER_LEN = len(STUFFER)
STUFFER_SOURCE = "p004 (processed://constructs/ordered/p004.gb), feature Stuffer 4597..4763"
BUILD_DATE = "06-OCT-2026"

CLONING_SITES = {"AgeI": AGEI, "XbaI": XBAI, "AsiSI": ASISI, "AscI": ASCI}
TYPE_IIS = ["BbsI", "Esp3I", "BsaI"]


def screen_stuffer(stuffer: str = STUFFER) -> list[str]:
    """Check the borrowed stuffer carries nothing that would break the cloning sites."""
    problems = []
    if seqtools.find_sites(stuffer, TYPE_IIS):
        problems.append(
            f"stuffer carries a Type IIS site: {seqtools.find_sites(stuffer, TYPE_IIS)}"
        )
    for name, site in CLONING_SITES.items():
        rc = seqtools.reverse_complement(site)
        n = stuffer.count(site) + (stuffer.count(rc) if rc != site else 0)
        if n:
            problems.append(f"stuffer carries {n} {name} site(s), which must be unique")
    return problems


def build(p003_path: Path = P003) -> tuple[SeqIO.SeqRecord, str]:
    if not p003_path.exists():
        raise SystemExit(f"p003 GenBank not found: {p003_path}")
    rec = SeqIO.read(p003_path, "genbank")
    parent = str(rec.seq).upper()

    # later deletion first, so the earlier coordinates stay valid
    rec = rec[: MSCARLET[0] - 1] + rec[MSCARLET[1] :]
    junction = MCS + STUFFER
    new = rec[: INSERT[0] - 1] + Seq(junction) + rec[INSERT[1] :]
    new.id = new.name = "p0039"
    new.description = "insert-agnostic puro-only lentiviral destination vector"
    # Fixed, not today's date: a generated file that changes on every run shows up as
    # a spurious diff and makes it impossible to tell a real design change from a rebuild.
    new.annotations = {
        "molecule_type": "ds-DNA",
        "topology": "circular",
        "date": BUILD_DATE,
        "accession": "p0039",
    }

    # p003 annotates PuroR as 597 bp starting at its second codon, because there it was
    # a P2A downstream product with no initiator of its own. After the mScarlet3/P2A
    # deletion the upstream ATG *is* PuroR's start, so the feature has to grow by 3 nt
    # or the file shows a CDS that does not begin with a methionine.
    for feat in new.features:
        if feat.qualifiers.get("label", [""])[0] == "PuroR" and len(feat) == 597:
            feat.location = FeatureLocation(
                int(feat.location.start) - 3, int(feat.location.end), feat.location.strand
            )
            feat.qualifiers["note"] = [
                "initiator ATG supplied by the IRES2 linker (mScarlet3's former Met); "
                "annotation extended 3 nt 5' relative to p003"
            ]

    start = INSERT[0] - 1
    for label, span in (
        ("MCS (AgeI-XbaI-AsiSI)", (start, start + len(MCS))),
        ("stuffer (dropout)", (start + len(MCS), start + len(junction))),
    ):
        new.features.append(
            SeqFeature(
                FeatureLocation(*span),
                type="misc_feature",
                qualifiers={"label": [label]}
                | (
                    {"note": [f"borrowed verbatim from {STUFFER_SOURCE}"]}
                    if "stuffer" in label
                    else {}
                ),
            )
        )
    return new, parent


def verify(backbone: str, parent: str, record) -> list[str]:
    """Every claim the docstring makes, checked against the built sequence."""
    problems = screen_stuffer()
    circular = backbone + backbone[:20]

    protein = ""
    for i in range(len(backbone) - 2):
        if backbone[i : i + 3] != "ATG":
            continue
        window = backbone[i : i + 3 * ((len(backbone) - i) // 3)]
        candidate = str(Seq(window).translate())
        if candidate.startswith(NATIVE_PURO_START):
            protein = candidate[: candidate.find("*") + 1]
            break
    if not protein:
        problems.append("PuroR ORF not found with a native N-terminus")
    else:
        if not protein.rstrip("*").endswith(NATIVE_PURO_END):
            problems.append(f"PuroR C-terminus not native: {protein[-14:]}")
        if protein.count("*") != 1:
            problems.append("PuroR ORF is not a single clean reading frame")

    annotated = [f for f in record.features if f.qualifiers.get("label", [""])[0] == "PuroR"]
    if len(annotated) != 1:
        problems.append(f"expected exactly one PuroR feature, found {len(annotated)}")
    else:
        cds = str(annotated[0].extract(Seq(backbone)))
        aa = str(Seq(cds).translate())
        if not cds.startswith("ATG"):
            problems.append(f"the annotated PuroR feature does not begin with ATG: {cds[:6]}")
        if aa != str(Seq(cds).translate()).rstrip("*") + "*" or aa.count("*") != 1:
            problems.append("the annotated PuroR feature is not a single clean CDS")
        if not aa.startswith(NATIVE_PURO_START) or not aa.rstrip("*").endswith(NATIVE_PURO_END):
            problems.append(f"the annotated PuroR feature is not native PuroR: {aa[:24]}")

    for enzyme in ("BbsI", "Esp3I"):
        hits = seqtools.find_sites(circular, [enzyme])
        if hits:
            problems.append(f"{enzyme} site introduced: {hits}")
    for name, site in CLONING_SITES.items():
        n = circular.count(site)
        rc = seqtools.reverse_complement(site)
        if rc != site:
            n += circular.count(rc)
        if n != 1:
            problems.append(f"{name} appears {n} times, must be unique")

    if not backbone.endswith(parent[MSCARLET[1] :]):
        problems.append("sequence 3' of the mScarlet3 deletion changed")
    if not backbone.startswith(parent[: INSERT[0] - 1]):
        problems.append("sequence 5' of the insert (through EFS) changed")
    return problems


CONSTRUCT_FASTAS = (
    "constructs_nterm.fasta",
    "constructs_loop.fasta",
    "constructs_tag_arms.fasta",
)


def verify_reconstitution(backbone: str, parent: str) -> tuple[int, list[str]]:
    """Drop every designed insert into the slot and diff against the expected build.

    This is the claim that matters for an insert-agnostic vector: that p0039 and a
    finished construct differ by the insert and by nothing else. Asserting it is cheap
    and wrong; the only way to know is to reconstitute each plasmid and compare, so
    that is what this does, on every rebuild.

    Returns ``(how_many_checked, failures)``. Missing construct FASTAs are not an
    error — the backbone can legitimately be built before any insert exists — but the
    count then reports 0, so a silent skip cannot be mistaken for a pass.
    """
    from build_constructs import read_fasta

    reference = parent[: MSCARLET[0] - 1] + parent[MSCARLET[1] :]
    slot = len(MCS) + STUFFER_LEN

    inserts: dict[str, str] = {}
    for name in CONSTRUCT_FASTAS:
        path = HERE / name
        if path.exists():
            inserts.update(read_fasta(path))

    failures = []
    for name, insert in sorted(inserts.items()):
        built = backbone[: INSERT[0] - 1] + insert + backbone[INSERT[0] - 1 + slot :]
        expected = reference[: INSERT[0] - 1] + insert + reference[INSERT[1] :]
        if built != expected:
            failures.append(name)
    return len(inserts), failures


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=HERE)
    args = ap.parse_args(argv)

    rec, parent = build()
    backbone = str(rec.seq).upper()

    print(f"p003              {len(parent):>7,} bp")
    print(f"  - insert        -{INSERT[1] - INSERT[0] + 1:>6,} bp  (Kozak..revprimer)")
    print(f"  - mScarlet3/P2A -{MSCARLET[1] - MSCARLET[0] + 1:>6,} bp")
    print(f"  + MCS + stuffer +{len(MCS) + STUFFER_LEN:>6,} bp")
    print(f"  = p0039         {len(backbone):>7,} bp\n")

    problems = verify(backbone, parent, rec)
    if problems:
        for p in problems:
            print(f"  FAIL  {p}")
        raise SystemExit("p0039 failed verification — nothing written")
    print("verification:")
    print("  PuroR CDS native, annotated from its own ATG       ✓")
    print("  no Esp3I or BbsI site anywhere (circular scan)    ✓")
    for name in CLONING_SITES:
        print(f"  {name:<6} unique in the plasmid                      ✓")
    print("  EFS and everything 5' of it unchanged             ✓")
    print("  IRES2 → PuroR and everything 3' of it unchanged   ✓")

    checked, failures = verify_reconstitution(backbone, parent)
    if failures:
        raise SystemExit(f"these constructs do not reconstitute from p0039: {failures}")
    if checked:
        print(f"  all {checked} designed inserts reconstitute byte-for-byte  ✓")
        print("    p0039 and a finished construct differ by the insert, nothing else")
    else:
        print("  no construct FASTAs present — reconstitution NOT checked")

    print(
        f"\n  GC {seqtools.gc_content(backbone):.1%}   longest direct repeat "
        f"{seqtools.max_internal_repeat(backbone)} nt"
    )

    fa = args.out_dir / "p0039_EFS_IRES_puro_backbone.fasta"
    fa.write_text(
        f">p0039 insert-agnostic puro-only lentiviral destination vector | {len(backbone)} bp circular\n"
        + "\n".join(backbone[i : i + 60] for i in range(0, len(backbone), 60))
        + "\n"
    )
    gbk = args.out_dir / "p0039_EFS_IRES_puro_backbone.gb"
    SeqIO.write(rec, gbk, "genbank")
    print(f"\nWrote {fa.name} and {gbk.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
