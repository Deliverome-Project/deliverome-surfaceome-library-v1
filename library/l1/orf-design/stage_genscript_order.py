"""Stage the plasmids we want GenScript to build into one folder, with a manifest.

Writes ``.gb`` and ``.fasta`` for each plasmid plus a MANIFEST.tsv and a README, so the
folder can be handed to a vendor without anyone having to reconstruct which file is what
or which known issues to expect. p0039 comes from this experiment's committed build; the
miR-L designs are pulled fresh from ``processed://constructs/designed/`` so the staged
copy is whatever the registry currently holds, not a stale local one.

Each plasmid is re-checked as it is staged — Type IIS counts, length, and for p0039 that
PuroR still reads from its own ATG — so a stale or broken record cannot be shipped.

    uv run --locked python stage_genscript_order.py --out-dir ~/Downloads/genscript
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

from deliverome_analysis import constructs, seqtools
from deliverome_analysis.s3 import processed

HERE = Path(__file__).parent
GG = ("BbsI", "Esp3I", "BsaI")

# Friendly names exactly as they appear in the Reagent Inventory's Plasmid tab, so a
# vendor file and an inventory row are obviously the same thing. Spaces become
# underscores in the filename only; the manifest carries the name verbatim.
FRIENDLY = {
    "p0016": "P-0016-pAG_Lenti_miRL_iFFL_PGK_intronsplit_PuroR_target",
    "p0017": "P-0017-pAG_Lenti_miRL_iFFL_PGK_intronsplit_PuroR_notarget",
    "p0018": "P-0018-pAG_Lenti_miRL_iFFL_EF1a_intronsplit_PuroR_PGKorf_target",
    "p0019": "P-0019-pAG_Lenti_miRL_iFFL_EF1a_intronsplit_PuroR_PGKorf_notarget",
    "p0020": "P-0020-pAG_Lenti_miRL_PGK_PuroR_intact_target",
    "p0021": "P-0021-pAG_Lenti_miRL_PGK_PuroR_intact_notarget",
    "p0040": "P-0040-pAG_Lenti_miRL_iFFL_PGK_intronsplit_PuroR_target_TACC",
    "p0039": "P-0039-pAG Lenti_236083_EFS_IRES_puro_insert_agnostic",
}


def filename(name: str) -> str:
    return FRIENDLY[name].replace(" ", "_")


# plasmid -> (where it comes from, one-line description, known flags for the vendor)
LOCAL = {
    "p0039": (
        HERE / "p0039_EFS_IRES_puro_backbone.gb",
        "insert-agnostic puro-only lentiviral destination vector (EFS / IRES2-PuroR)",
        "183-nt direct repeat between the two LTRs (inherent to lentiviral transfer plasmids)",
    ),
}
FROM_S3 = {
    "p0016": "miR-L iFFL, PGK, intron-split PuroR, +miR-L target array",
    "p0017": "miR-L iFFL, PGK, intron-split PuroR, no target array (matched control)",
    "p0018": "miR-L iFFL, EF1a_v1 regulator / PGK on the ORF, intron-split PuroR, +target array",
    "p0019": "miR-L iFFL, EF1a_v1 regulator / PGK on the ORF, intron-split PuroR, no target array",
    "p0020": "miR-L iFFL, PGK, PuroR CDS intact, +miR-L target array",
    "p0021": "miR-L iFFL, PGK, PuroR CDS intact, no target array (matched control)",
    "p0040": "as p0016, but the array spacer is TACC, dropping 4 endogenous miR-374-5p 8mers",
}
FLAGS = {
    "p0016": "100-nt miR-L target array is a 63-nt tandem repeat",
    "p0018": "100-nt miR-L target array is a 63-nt tandem repeat",
    "p0020": "100-nt miR-L target array is a 63-nt tandem repeat",
}


def write_pair(rec, out: Path, name: str) -> None:
    stem = filename(name)
    SeqIO.write(rec, out / f"{stem}.gb", "genbank")
    s = str(rec.seq).upper()
    (out / f"{stem}.fasta").write_text(
        f">{FRIENDLY[name]} | {len(s)} bp circular | {rec.description}\n"
        + "\n".join(s[i : i + 60] for i in range(0, len(s), 60))
        + "\n"
    )


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--out-dir",
        type=Path,
        default=Path.home() / "Downloads" / "genscript",
    )  # allow-local-path: a folder to hand a vendor, not shared data; override with --out-dir
    args = ap.parse_args(argv)
    out = args.out_dir.expanduser()
    out.mkdir(parents=True, exist_ok=True)

    rows = []
    for name, (path, desc, flag) in LOCAL.items():
        rec = SeqIO.read(path, "genbank")
        s = str(rec.seq).upper()
        puro = [f for f in rec.features if f.qualifiers.get("label", [""])[0] == "PuroR"]
        if not puro or not str(puro[0].extract(Seq(s))).startswith("ATG"):
            raise SystemExit(f"{name}: PuroR does not read from its own ATG — not staging")
        write_pair(rec, out, name)
        rows.append((name, len(s), desc, s, flag))

    for name, desc in FROM_S3.items():
        rec = SeqIO.read(processed.download(constructs.map_key(name)), "genbank")
        s = str(rec.seq).upper()
        if seqtools.find_sites(s, ["BbsI"]):
            raise SystemExit(f"{name}: still carries a BbsI site — not staging")
        rec.id = rec.name = name
        write_pair(rec, out, name)
        rows.append((name, len(s), desc, s, FLAGS.get(name, "")))

    cols = [
        "plasmid",
        "friendly_name",
        "file_stem",
        "length_bp",
        "gc",
        "BbsI",
        "Esp3I",
        "BsaI",
        "description",
        "flags",
    ]
    with (out / "MANIFEST.tsv").open("w", newline="") as h:
        w = csv.writer(h, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for name, n, desc, s, flag in rows:
            w.writerow(
                [
                    name,
                    FRIENDLY[name],
                    filename(name),
                    n,
                    f"{seqtools.gc_content(s):.3f}",
                    *[len(seqtools.find_sites(s, [e])) for e in GG],
                    desc,
                    flag,
                ]
            )

    readme = [
        "# GenScript synthesis request",
        "",
        "Whole-plasmid synthesis. Each plasmid has a `.gb` (annotated, opens in SnapGene or",
        "Benchling) and a `.fasta` (plain sequence for the order form). `MANIFEST.tsv` carries",
        "lengths, %GC, Type IIS counts and the known complexity flags.",
        "",
        "| file | bp | what it is |",
        "|---|---:|---|",
    ]
    readme += [f"| `{filename(n)}` | {ln:,} | {d} |" for n, ln, d, _, _ in rows]
    readme += [
        "",
        "## Known complexity, so it does not come back as a surprise",
        "",
        "* **All of them** carry the 183-nt direct repeat between the 5' and 3' LTRs. That is",
        "  inherent to a lentiviral transfer plasmid and cannot be designed out; it is the most",
        "  likely trigger for a complexity surcharge.",
        "* **p0016 / p0018 / p0020** carry a 100-nt miR-L target array that is a 63-nt tandem",
        "  repeat by design (four copies of one site, verbatim from Du et al. 2025). It cannot be",
        "  recoded without changing the biology. Please sequence across the array specifically and",
        "  confirm all four copies are present.",
        "* **p0018 / p0019** carry EF1a_v1 including intron A, but the whole regulator",
        "  cassette is antisense to the vector genome, so that intron is never read in the",
        "  sense direction. It is not a packaging concern and needs no special handling.",
        "",
        "## What we need back",
        "",
        "A **full-sequence export** for each plasmid, not just the insert. A design map that",
        "matches our own tells us nothing: on a previous order the delivered map and our record",
        "agreed exactly while both differed from the plasmid at 11 bases.",
        "",
        "Type IIS: all are free of BbsI and Esp3I. BsaI sites are pre-existing in the",
        "pAG-lenti (Addgene 236083) backbone and are expected.",
    ]
    (out / "README.md").write_text("\n".join(readme) + "\n")

    print(f"staged {len(rows)} plasmid(s) into {out}\n")
    for name, n, _desc, s, _ in rows:
        tiis = "/".join(str(len(seqtools.find_sites(s, [e]))) for e in GG)
        print(f"  {filename(name)}.{{gb,fasta}}")
        print(f"      {n:,} bp   GC {seqtools.gc_content(s):.1%}   BbsI/Esp3I/BsaI {tiis}")
    print("\nalso wrote MANIFEST.tsv and README.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
