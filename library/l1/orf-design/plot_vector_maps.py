"""Linear comparative map of the miR-L vector series, drawn over the integrating region.

A circular map per plasmid hides the only thing that matters here: these ten vectors
share a backbone and differ in a handful of cassettes. Drawing the **LTR-to-LTR region**
-- the part that actually integrates -- as one track per construct, aligned at the 5'
LTR, puts those differences side by side.

Everything is read back out of the GenBank records on S3, so the figure cannot drift
from what is staged for the vendor.

Two things the drawing has to get right or it misleads. The regulator cassette is
**antisense**, so each feature is drawn as an arrow pointing the way it is transcribed;
a flat box would imply the whole provirus reads left to right. And the intron carrying
the miR-L hairpin sits *inside* the PuroR CDS in six of the vectors, so it is drawn
nested on top of PuroR rather than suppressed as an overlap -- that nesting is exactly
what distinguishes those arms from the 3'UTR ones.

A second panel zooms the 3' cassette, because the site array that the whole experiment
turns on is ~100 bp in an 8 kb genome and is invisible at full scale.

    uv run --locked python plot_vector_maps.py
"""

from __future__ import annotations

import argparse
import re
import warnings
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from Bio import BiopythonWarning, SeqIO

from deliverome_analysis.cloning import constructs
from deliverome_analysis.core.plotting import (
    get_categorical_colors,
    save_figure,
    setup_plotting_style,
)
from deliverome_analysis.core.s3 import processed

HERE = Path(__file__).parent

# construct -> (row label, one-line descriptor)
ROWS = [
    ("p0039", "P-0039", "destination vector, insert-agnostic"),
    ("p0016", "P-0016", "PGK reg · intron in PuroR CDS · target (TACT)"),
    ("p0040", "P-0040", "PGK reg · intron in PuroR CDS · target (TACC, own barcode)"),
    ("p0017", "P-0017", "PGK reg · intron in PuroR CDS · scramble"),
    ("p0018", "P-0018", "EF1a reg / PGK ORF · intron in CDS · target (TACT)"),
    ("p0019", "P-0019", "EF1a reg / PGK ORF · intron in CDS · scramble"),
    ("p0020", "P-0020", "PGK reg · intron in PuroR 3'UTR · target (TACT)"),
    ("p0021", "P-0021", "PGK reg · intron in PuroR 3'UTR · scramble"),
]

# label prefix -> (class, short name on the figure). Order matters: first match wins.
CLASSES = [
    ("5' LTR", "ltr", "LTR"),
    ("3' LTR", "ltr", "LTR"),
    ("HIV-1 Psi", "vector", None),
    ("RRE", "vector", None),
    ("cPPT", "vector", None),
    ("WPRE", "vector", "WPRE"),
    ("bGH poly(A)", "term", "pA"),
    ("PGK_v1", "promoter", "PGK"),
    ("EFS_v1", "promoter", "EFS"),
    ("EF1a_v1", "promoter", "EF1a"),
    ("PuroR", "marker", "PuroR"),
    ("intron + miR-L", "mirna", "intron+miR-L"),
    ("IRES2", "vector", "IRES"),
    ("TFRC", "orf", "TFRC"),
    ("Linker_STOP", "tag", "tags"),
    ("4x17nt", "sites", "4x17 target"),
    ("scramble sites", "sites", "scramble"),
    ("O-0013", "amplicon", None),
    ("20bp-BC", "amplicon", None),
    ("10_bc_CS1", "amplicon", "BC / CS1"),
    ("stuffer", "cloning", "stuffer"),
    ("MCS", "cloning", "MCS"),
]
CLASS_ORDER = [
    ("ltr", "LTR"),
    ("vector", "lentiviral elements"),
    ("promoter", "promoter"),
    ("orf", "ORF (TFRC)"),
    ("tag", "SpyTag/ALFA cassette"),
    ("amplicon", "handle · barcode · 10x CS1"),
    ("sites", "miR-L site array"),
    ("marker", "PuroR"),
    ("mirna", "intron + miR-L hairpin"),
    ("term", "poly(A)"),
    ("cloning", "cloning stuffer / MCS"),
]


NESTED = ("intron + miR-L",)  # lies inside another feature and must not be hidden


def arrow(ax, a, b, y, color, strand, h=0.52, z=2):
    """A feature as a directional arrow, so antisense cassettes read as antisense."""
    tip = min((b - a) * 0.45, 70)
    if strand == -1:
        pts = [(b, y + h / 2), (a + tip, y + h / 2), (a, y), (a + tip, y - h / 2), (b, y - h / 2)]
    else:
        pts = [(a, y - h / 2), (b - tip, y - h / 2), (b, y), (b - tip, y + h / 2), (a, y + h / 2)]
    ax.add_patch(
        mpatches.Polygon(pts, closed=True, facecolor=color, edgecolor="white", lw=0.5, zorder=z)
    )


def classify(label: str):
    for prefix, cls, short in CLASSES:
        if label.startswith(prefix):
            return cls, short
    return None, None


def genome(rec):
    """(start, end) of the LTR-to-LTR region: what actually integrates."""
    five = next(f for f in rec.features if f.qualifiers.get("label", [""])[0].startswith("5' LTR"))
    three = next(f for f in rec.features if f.qualifiers.get("label", [""])[0].startswith("3' LTR"))
    return int(five.location.start), int(three.location.end)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=HERE / "figures")
    args = ap.parse_args(argv)

    setup_plotting_style()
    colors = dict(
        zip([c for c, _ in CLASS_ORDER], get_categorical_colors(len(CLASS_ORDER)), strict=True)
    )

    records = {}
    for name, _, _ in ROWS:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", BiopythonWarning)
            records[name] = SeqIO.read(processed.download(constructs.map_key(name)), "genbank")

    span = max(genome(r)[1] - genome(r)[0] for r in records.values())
    fig, (ax, az) = plt.subplots(
        2,
        1,
        figsize=(13.5, 0.62 * len(ROWS) + 5.2),
        gridspec_kw={"height_ratios": [len(ROWS), len(ROWS) * 0.72], "hspace": 0.38},
    )

    for row, (name, pnum, desc) in enumerate(ROWS):
        rec = records[name]
        g0, g1 = genome(rec)
        y = len(ROWS) - row - 1
        ax.plot([0, g1 - g0], [y, y], color="0.75", lw=1.0, zorder=1, solid_capstyle="butt")
        drawn, nested = [], []
        for f in sorted(rec.features, key=lambda f: -(int(f.location.end) - int(f.location.start))):
            lab = f.qualifiers.get("label", [""])[0]
            cls, short = classify(lab)
            if cls is None:
                continue
            a, b = int(f.location.start) - g0, int(f.location.end) - g0
            if b <= 0 or a >= g1 - g0:
                continue
            if lab.startswith(NESTED):
                nested.append((a, b, cls, short, f.location.strand))
                continue
            if any(not (b <= pa or a >= pb) for pa, pb in drawn):
                continue
            drawn.append((a, b))
            arrow(ax, a, b, y, colors[cls], f.location.strand)
            if short and (b - a) > span * 0.045:
                ax.text(
                    (a + b) / 2,
                    y,
                    short,
                    ha="center",
                    va="center",
                    fontsize=6.4,
                    color="white",
                    zorder=4,
                )
        for a, b, cls, _short, strand in nested:
            arrow(ax, a, b, y, colors[cls], strand, h=0.26, z=5)
        ax.text(-span * 0.012, y, pnum, ha="right", va="center", fontsize=8, fontweight="bold")
        ax.text(-span * 0.012, y - 0.30, desc, ha="right", va="center", fontsize=6.2, color="0.35")

    ax.set_xlim(-span * 0.30, span * 1.02)
    ax.set_ylim(-0.9, len(ROWS) - 0.3)
    ax.set_yticks([])
    ax.set_xlabel("position in the integrating region, 5' LTR to 3' LTR (bp)")
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.legend(
        handles=[mpatches.Patch(facecolor=colors[c], label=lab) for c, lab in CLASS_ORDER],
        loc="upper center",
        bbox_to_anchor=(0.5, -0.13),
        ncol=4,
        frameon=False,
        fontsize=7.5,
    )
    # --- zoom: the 3' cassette, where the ten vectors actually differ ---
    ZOOM = (
        "O-0013",
        "20bp-BC",
        "Potential_revprimer",
        "RT ",
        "10_bc_CS1",
        "AsiSI",
        "4x17nt",
        "scramble sites",
        "PmeI",
        "miR-L target site",
        "scramble site",
    )
    zspan = 0
    for name, _pnum, _desc in ROWS:
        rec = records[name]
        hits = [f for f in rec.features if f.qualifiers.get("label", [""])[0].startswith(ZOOM)]
        if not hits:
            continue
        z0 = min(int(f.location.start) for f in hits)
        zspan = max(zspan, max(int(f.location.end) for f in hits) - z0)
    for row, (name, pnum, _desc) in enumerate(ROWS):
        rec = records[name]
        hits = [f for f in rec.features if f.qualifiers.get("label", [""])[0].startswith(ZOOM)]
        y = len(ROWS) - row - 1
        if not hits:
            az.text(
                zspan / 2,
                y,
                "no barcode/site cassette in the destination vector",
                ha="center",
                va="center",
                fontsize=6.8,
                color="0.45",
                style="italic",
            )
            az.text(-zspan * 0.012, y, pnum, ha="right", va="center", fontsize=8, fontweight="bold")
            continue
        z0 = min(int(f.location.start) for f in hits)
        az.plot([0, zspan], [y, y], color="0.85", lw=0.8, zorder=1)
        for f in sorted(hits, key=lambda f: int(f.location.start)):
            lab = f.qualifiers.get("label", [""])[0]
            # Match the numbered copies EXACTLY: "scramble site" is a prefix of the
            # array's own label "scramble sites [...]", so a prefix test would skip
            # the array itself and leave the control rows blank.
            if re.fullmatch(r"(miR-L target|scramble) site \d of 4", lab):
                continue  # drawn as part of the parent array, marked pale below
            cls, _ = classify(lab)
            a, b = int(f.location.start) - z0, int(f.location.end) - z0
            col = colors[cls] if cls else "0.6"
            az.add_patch(
                mpatches.Rectangle(
                    (a, y - 0.26), b - a, 0.52, facecolor=col, edgecolor="white", lw=0.5, zorder=2
                )
            )
        arr = [
            f
            for f in hits
            if f.qualifiers.get("label", [""])[0].startswith(("4x17nt", "scramble sites"))
        ]
        for f in arr:
            a = int(f.location.start) - z0
            for k in range(4):
                s = a + 4 + k * 21
                az.add_patch(
                    mpatches.Rectangle(
                        (s, y - 0.26),
                        17,
                        0.52,
                        facecolor="white",
                        alpha=0.38,
                        edgecolor="none",
                        zorder=3,
                    )
                )
        az.text(-zspan * 0.012, y, pnum, ha="right", va="center", fontsize=8, fontweight="bold")
        spacer = "TACC" if name == "p0040" else "TACT"
        kind = "scramble" if name in ("p0017", "p0019", "p0021") else f"target, {spacer} spacer"
        az.text(zspan * 1.015, y, kind, ha="left", va="center", fontsize=6.8, color="0.35")
    az.set_xlim(-zspan * 0.30, zspan * 1.30)
    az.set_ylim(-0.9, len(ROWS) - 0.3)
    az.set_yticks([])
    az.set_xlabel(
        "3' cassette detail, bp from the amplicon handle  (pale bands = the four 17-nt "
        "sites; grey = AsiSI/PmeI flanks and the RT site)"
    )
    for side in ("top", "right", "left"):
        az.spines[side].set_visible(False)

    save_figure(fig, "s4e9_vector_maps", output_dir=args.out_dir)
    print(f"wrote {args.out_dir / 's4e9_vector_maps.png'}")
    for name, pnum, _ in ROWS:
        g0, g1 = genome(records[name])
        print(f"  {pnum}  plasmid {len(records[name].seq):,} bp   LTR-to-LTR {g1 - g0:,} bp")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
