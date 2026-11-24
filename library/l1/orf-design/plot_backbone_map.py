"""Circular map of p0039, the insert-agnostic puro-only destination vector.

Reads the GenBank written by empty_backbone.py rather than recomputing the sequence,
so the map can only ever show what was actually built. Unique cloning sites are read
back out of the sequence, not typed in, for the same reason.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO

from deliverome_analysis import seqtools
from deliverome_analysis.plotting import (
    get_categorical_colors,
    save_figure,
    setup_plotting_style,
)

HERE = Path(__file__).parent
GBK = HERE / "p0039_EFS_IRES_puro_backbone.gb"

# Features worth showing, in plasmid order. Anything smaller is noise on a 7.8 kb map.
SHOW = [
    ("5' LTR (truncated)", "LTR"),
    ("HIV-1 Psi", "element"),
    ("RRE", "element"),
    ("cPPT/CTS", "element"),
    ("EFS_v1", "promoter"),
    ("MCS (AgeI-XbaI-AsiSI)", "cloning"),
    ("stuffer (dropout)", "cloning"),
    ("IRES2_v1", "element"),
    ("PuroR", "cds"),
    ("WPRE", "element"),
    ("3' LTR (Delta-U3)", "LTR"),
    ("SV40 poly(A) signal", "element"),
    ("AmpR", "cds"),
    ("ori", "element"),
]
SITES = {"AgeI": "ACCGGT", "XbaI": "TCTAGA", "AsiSI": "GCGATCGC", "AscI": "GGCGCGCC"}


def collect(rec) -> list[tuple[str, str, int, int]]:
    seen, out = set(), []
    for label, kind in SHOW:
        for f in rec.features:
            name = f.qualifiers.get("label", [""])[0]
            if name != label or name in seen:
                continue
            seen.add(name)
            out.append((name, kind, int(f.location.start), int(f.location.end)))
            break
    return out


def arc(ax, start, end, total, radius, width, color, label, zorder=2):
    t0, t1 = 90 - 360 * start / total, 90 - 360 * end / total
    ax.add_patch(
        mpatches.Wedge(
            (0, 0),
            radius + width / 2,
            t1,
            t0,
            width=width,
            facecolor=color,
            edgecolor="white",
            linewidth=0.8,
            zorder=zorder,
        )
    )
    mid = np.deg2rad((t0 + t1) / 2)
    return mid, label


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=HERE / "figures")
    args = ap.parse_args(argv)

    rec = SeqIO.read(GBK, "genbank")
    seq = str(rec.seq).upper()
    total = len(seq)
    feats = collect(rec)
    missing = [lbl for lbl, _ in SHOW if lbl not in {f[0] for f in feats}]
    if missing:
        raise RuntimeError(f"features not found in the GenBank: {missing}")

    setup_plotting_style()
    fig = plt.figure(figsize=(17, 9.5))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.05, 1.0], wspace=0.05)
    ax = fig.add_subplot(gs[0, 0])
    axl = fig.add_subplot(gs[0, 1])

    palette = get_categorical_colors(6)
    colour = {
        "LTR": "#9aa0a6",
        "element": palette[0],
        "promoter": palette[2],
        "cds": palette[1],
        "cloning": palette[3],
    }

    R, W = 1.0, 0.11
    ax.add_patch(
        mpatches.Wedge(
            (0, 0), R + W / 2, 0, 360, width=W, facecolor="#ececec", edgecolor="none", zorder=1
        )
    )
    suppress = {"MCS (AgeI-XbaI-AsiSI)", "stuffer (dropout)"}
    rung = 0
    for name, kind, s, e in feats:
        mid, label = arc(ax, s, e, total, R, W, colour[kind], name, zorder=2)
        if name in suppress:  # too small to label here; shown in the blow-up instead
            continue
        lr = 1.34 + 0.21 * (rung % 2)  # stagger: adjacent small features collide at one radius
        rung += 1
        x, y = lr * np.cos(mid), lr * np.sin(mid)
        ax.annotate(
            f"{label}\n{e - s:,} bp",
            xy=((R + W / 2) * np.cos(mid), (R + W / 2) * np.sin(mid)),
            xytext=(x, y),
            ha="left" if x > 0 else "right",
            va="center",
            fontsize=10.5,
            arrowprops={
                "arrowstyle": "-",
                "lw": 0.7,
                "color": "#999999",
                "shrinkA": 0,
                "shrinkB": 2,
            },
        )

    # the four unique cutters sit inside 330 bp, so one callout, not four radial lines
    sites = {}
    for enzyme, site in SITES.items():
        if seq.count(site) != 1:
            raise RuntimeError(f"{enzyme} is not unique in the built plasmid")
        sites[enzyme] = seq.find(site) + 1
    lo, hi = min(sites.values()), max(sites.values())
    ang = np.deg2rad(90 - 360 * ((lo + hi) / 2) / total)
    ax.annotate(
        "cloning region\n"
        + "   ".join(f"{e} {p:,}" for e, p in sorted(sites.items(), key=lambda kv: kv[1])),
        xy=((R - W / 2) * np.cos(ang), (R - W / 2) * np.sin(ang)),
        xytext=(0.0, -1.82),
        ha="center",
        va="center",
        fontsize=10.5,
        zorder=6,
        bbox={"boxstyle": "round,pad=0.4", "fc": "white", "ec": "#333333", "lw": 0.9},
        arrowprops={"arrowstyle": "-", "lw": 0.9, "color": "#333333"},
    )

    ax.text(0, 0.16, "p0039", ha="center", va="center", fontsize=30, weight="bold")
    ax.text(0, 0.02, f"{total:,} bp", ha="center", va="center", fontsize=15)
    ax.text(
        0,
        -0.12,
        f"GC {seqtools.gc_content(seq):.0%}",
        ha="center",
        va="center",
        fontsize=12,
        color="#666666",
    )
    ax.set_xlim(-1.95, 1.95)
    ax.set_ylim(-2.05, 1.62)
    ax.set_aspect("equal")
    ax.axis("off")

    # ---- linear blow-up of the EFS -> IRES2 junction -------------------------
    efs = next(f for f in feats if f[0] == "EFS_v1")
    ires = next(f for f in feats if f[0] == "IRES2_v1")
    x0, x1 = efs[3] - 230, ires[2] + 230
    mcs_end = next(f for f in feats if f[0].startswith("MCS"))[3]
    stuf = next(f for f in feats if f[0].startswith("stuffer"))
    blocks = [
        ("EFS_v1\npromoter", efs[2], efs[3], colour["promoter"]),
        ("MCS", next(f for f in feats if f[0].startswith("MCS"))[2], mcs_end, colour["cloning"]),
        ("stuffer\n(dropout)", stuf[2], stuf[3], "#c9c2e8"),
        ("IRES2_v1", ires[2], ires[3], colour["element"]),
    ]
    for label, s, e, c in blocks:
        axl.add_patch(mpatches.Rectangle((s, 0), e - s, 1, facecolor=c, edgecolor="white", lw=1.2))
        # EFS and IRES2 run off both ends of this window, so centre each label on the
        # part that is actually visible or it lands outside the axes entirely.
        axl.text(
            (max(s, x0) + min(e, x1)) / 2,
            0.5,
            label,
            ha="center",
            va="center",
            fontsize=10.5,
            color="white",
            weight="bold",
        )
    for i, (enzyme, pos) in enumerate(sorted(sites.items(), key=lambda kv: kv[1])):
        y = 1.30 + 0.40 * i
        axl.plot([pos, pos], [1.0, y], color="#222222", lw=1.1)
        axl.text(pos, y + 0.06, f"{enzyme}\n{pos:,}", ha="center", va="bottom", fontsize=10.5)
    axl.annotate(
        "",
        xy=(mcs_end, -0.28),
        xytext=(sites["AscI"] + 8, -0.28),
        arrowprops={"arrowstyle": "<->", "lw": 1.1, "color": "#333333"},
    )
    axl.text(
        (mcs_end + sites["AscI"]) / 2,
        -0.52,
        f"{sites['AscI'] + 8 - mcs_end} bp drops out on an AgeI/AscI cut",
        ha="center",
        va="top",
        fontsize=10.5,
    )
    axl.text(x0 + 4, 3.25, "cloning junction, to scale", fontsize=12.5, weight="bold")
    axl.set_xlim(x0, x1)
    axl.set_ylim(-1.0, 3.5)
    axl.axis("off")

    save_figure(fig, "p0039_plasmid_map", output_dir=args.out_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
