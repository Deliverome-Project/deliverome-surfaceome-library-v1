"""Schematic of a loop construct and a snorkel construct (reviewer brief, "Snorkels").

Two panels, not to scale and not data:

- (a) a four-helix multipass protein with both termini cytoplasmic, tagged by inserting the
  SpyTag003-ALFA cassette into an extracellular loop;
- (b) the same protein with the snorkel module appended to its cytoplasmic C-terminus: Gly/Ser
  linker, mouse PDGFRB juxtamembrane residues and TM helix in reversed orientation, a second
  linker, then the cassette, which the added helix carries outside
  ([Brown et al., 2013](https://pubmed.ncbi.nlm.nih.gov/24023844/)).

Every residue count in the labels is read from the code that builds the constructs, never
typed in: the four module parts from ``snorkel_design.py`` (s5e13; parsed with ``ast`` so that
importing its model stack is not needed), checked against ``l1_snorkel.PREFIX_AA`` (the module
the order is built with), and the cassette and its two tags from
``deliverome_analysis.plasmid_qc``.

    uv run --locked python plot_snorkel_schematic.py
"""

from __future__ import annotations

import argparse
import ast
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch, Rectangle

from deliverome_analysis.cloning.plasmid_qc import P003_CASSETTE_PEPTIDE, TAG_PEPTIDES
from deliverome_analysis.core.plotting import (
    COLORS,
    SEQUENTIAL_PALETTES,
    save_figure,
    setup_plotting_style,
)

HERE = Path(__file__).resolve().parent
SNORKEL_DESIGN = HERE.parents[1] / "l1" / "tag-site-benchmarking" / "snorkel_design.py"
OUT_DIR = HERE / "reviewer-package" / "figures"
BASENAME = "snorkel_schematic"

INK = COLORS["dark"]
GREY = COLORS["neutral"]
PROTEIN = COLORS["primary"]
CASSETTE = COLORS["secondary"]  # the internal-tag color of the L1 construct figure
SNORKEL_TM = COLORS["accent"]
JUXTA = SEQUENTIAL_PALETTES["amber"][2]
MEMBRANE = SEQUENTIAL_PALETTES["amber"][5]
FS = 10.5
CAS_W = 3.1  # cassette box width, in axis units


def module_parts() -> dict[str, str]:
    """The snorkel module's four parts, read from snorkel_design.py and checked against the build."""
    tree = ast.parse(SNORKEL_DESIGN.read_text())
    names = ("L1_GS", "L1_JUXTA", "TMD", "L2")
    found = {
        t.id: node.value.value
        for node in tree.body
        if isinstance(node, ast.Assign) and isinstance(node.value, ast.Constant)
        for t in node.targets
        if isinstance(t, ast.Name) and t.id in names
    }
    missing = set(names) - set(found)
    if missing:
        raise SystemExit(f"snorkel_design.py no longer defines {sorted(missing)}")
    from l1_snorkel import PREFIX_AA

    if "".join(found[n] for n in names) != PREFIX_AA:
        raise SystemExit("snorkel_design.py parts do not spell l1_snorkel.PREFIX_AA")
    for tag in TAG_PEPTIDES.values():
        assert tag in P003_CASSETTE_PEPTIDE
    return found


def wiggle(ax, p0, p1, n_waves=4, amp=0.12, **kw):
    """A flexible linker: a sine wave from p0 to p1."""
    p0, p1 = np.asarray(p0, float), np.asarray(p1, float)
    t = np.linspace(0, 1, 200)
    d = p1 - p0
    normal = np.array([-d[1], d[0]]) / np.hypot(*d)
    pts = p0 + np.outer(t, d) + np.outer(amp * np.sin(2 * np.pi * n_waves * t), normal)
    ax.plot(pts[:, 0], pts[:, 1], **kw)


def helix(ax, x, y0, y1, color, w=0.62, label=None, text_color="white"):
    ax.add_patch(
        FancyBboxPatch(
            (x - w / 2, y0),
            w,
            y1 - y0,
            boxstyle="round,pad=0,rounding_size=0.25",
            facecolor=color,
            edgecolor="none",
            zorder=3,
        )
    )
    if label:
        ax.text(
            x,
            (y0 + y1) / 2,
            label,
            ha="center",
            va="center",
            color=text_color,
            fontsize=FS - 1,
            fontweight="bold",
            zorder=4,
        )


def loop(ax, x0, x1, y, height, **kw):
    """A loop between two helix ends, bulging up (height > 0) or down."""
    t = np.linspace(0, np.pi, 100)
    xs = x0 + (x1 - x0) * (1 - np.cos(t)) / 2
    ax.plot(xs, y + height * np.sin(t), **kw)


def cassette(ax, x, y, w=CAS_W, h=0.62, text="SpyTag003\u2013ALFA"):
    ax.add_patch(
        FancyBboxPatch(
            (x - w / 2, y),
            w,
            h,
            boxstyle="round,pad=0,rounding_size=0.15",
            facecolor=CASSETTE,
            edgecolor="none",
            zorder=5,
        )
    )
    ax.text(
        x,
        y + h / 2,
        text,
        ha="center",
        va="center",
        color="white",
        fontsize=FS - 1.5,
        fontweight="bold",
        zorder=6,
    )


def draw_protein(ax, cassette_in_loop: bool):
    """Four helices; ECL1 (H1-H2) and ECL2 (H3-H4) outside, N- and C-termini inside."""
    xs = [1.4, 2.5, 3.6, 4.7]
    lo, hi = 3.75, 6.25
    line = dict(color=PROTEIN, lw=2.4, solid_capstyle="round", zorder=2)
    for i, x in enumerate(xs, 1):
        helix(ax, x, lo, hi, PROTEIN, label=f"H{i}")
    loop(ax, xs[0], xs[1], hi, 0.55, **line)  # ECL1, short
    loop(ax, xs[1], xs[2], lo, -0.55, **line)  # ICL
    if cassette_in_loop:
        loop(ax, xs[2], xs[3], hi, 1.55, **line)
        cassette(ax, (xs[2] + xs[3]) / 2, hi + 1.25)
    else:
        loop(ax, xs[2], xs[3], hi, 0.75, **line)
    # N-terminus, inside
    ax.plot([xs[0], xs[0] - 0.35], [lo, lo - 0.85], **line)
    ax.text(xs[0] - 0.42, lo - 1.0, "N", ha="center", va="top", fontsize=FS, color=INK)
    return xs, lo, hi, line


def membrane(ax, xmax, label=True):
    ax.add_patch(Rectangle((0, 4.0), xmax, 2.0, facecolor=MEMBRANE, edgecolor="none", zorder=0))
    ax.text(0.1, 7.6, "out", ha="left", va="center", fontsize=FS, color=GREY, style="italic")
    ax.text(0.1, 1.3, "in", ha="left", va="center", fontsize=FS, color=GREY, style="italic")
    if not label:
        return
    ax.text(
        xmax - 0.1,
        5.0,
        "plasma\nmembrane",
        ha="right",
        va="center",
        fontsize=FS - 1.5,
        color=GREY,
        linespacing=1.0,
    )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = ap.parse_args()

    parts = module_parts()
    n_link1, n_juxta = len(parts["L1_GS"]), len(parts["L1_JUXTA"])
    n_tm, n_link2 = len(parts["TMD"]), len(parts["L2"])
    n_cas = len(P003_CASSETTE_PEPTIDE)
    n_module = n_link1 + n_juxta + n_tm + n_link2 + n_cas

    setup_plotting_style()
    fig, (ax_a, ax_b) = plt.subplots(
        1, 2, figsize=(11.5, 4.6), gridspec_kw={"width_ratios": [7.0, 11.0], "wspace": 0.04}
    )
    for ax, xmax in ((ax_a, 7.0), (ax_b, 11.0)):
        ax.set_xlim(0, xmax)
        ax.set_ylim(0.6, 9.6)
        ax.set_aspect("equal")
        ax.axis("off")
        membrane(ax, xmax, label=ax is ax_a)

    # (a) cassette inserted into an extracellular loop
    xs, lo, hi, line = draw_protein(ax_a, cassette_in_loop=True)
    ax_a.plot([xs[3], xs[3] + 0.35], [lo, lo - 0.85], **line)
    ax_a.text(xs[3] + 0.42, lo - 1.0, "C", ha="center", va="top", fontsize=FS, color=INK)
    ax_a.text(
        (xs[2] + xs[3]) / 2,
        hi + 2.05,
        f"cassette ({n_cas} aa) in a loop",
        ha="center",
        va="bottom",
        fontsize=FS,
        color=INK,
    )
    ax_a.text(0.0, 9.6, "a", fontsize=FS + 3, fontweight="bold", va="top", color=INK)

    # (b) the snorkel module on the cytoplasmic C-terminus
    xs, lo, hi, line = draw_protein(ax_b, cassette_in_loop=False)
    c_end = (xs[3] + 0.35, lo - 0.85)
    ax_b.plot([xs[3], c_end[0]], [lo, c_end[1]], **line)
    ax_b.text(c_end[0] - 0.05, c_end[1] - 0.15, "C", ha="right", va="top", fontsize=FS, color=INK)
    x_tm = 6.9
    juxta_start = (x_tm - 0.45, 2.55)
    juxta_end = (x_tm, lo - 0.05)
    link = dict(color=GREY, lw=1.8, zorder=2)
    wiggle(ax_b, c_end, juxta_start, n_waves=3, **link)
    ax_b.plot(
        *zip(juxta_start, juxta_end, strict=True),
        color=JUXTA,
        lw=5,
        solid_capstyle="round",
        zorder=3,
    )
    helix(ax_b, x_tm, lo, hi, SNORKEL_TM, label="TM", text_color=INK)
    cas_x, cas_y = x_tm + 0.6, hi + 1.3
    wiggle(ax_b, (x_tm, hi), (cas_x, cas_y), n_waves=2, **link)
    cassette(ax_b, cas_x, cas_y)

    lab = dict(fontsize=FS - 0.5, color=INK, va="center", ha="left")
    arrow = dict(arrowstyle="-", color=GREY, lw=0.8)
    ax_b.annotate(
        f"Gly/Ser linker ({n_link1} aa)",
        xy=(5.75, 2.55),
        xytext=(4.0, 1.05),
        arrowprops=arrow,
        **lab,
    )
    ax_b.annotate(
        f"PDGFRB juxtamembrane ({n_juxta} aa)",
        xy=(6.68, 2.95),
        xytext=(7.55, 2.2),
        arrowprops=arrow,
        **lab,
    )
    ax_b.annotate(
        f"mouse PDGFRB TM helix\n({n_tm} aa), reversed",
        xy=(x_tm + 0.32, 5.0),
        xytext=(7.55, 5.0),
        arrowprops=arrow,
        linespacing=1.1,
        **lab,
    )
    ax_b.annotate(
        f"linker ({n_link2} aa)", xy=(x_tm + 0.12, 6.9), xytext=(7.55, 7.0), arrowprops=arrow, **lab
    )
    ax_b.annotate(
        f"cassette ({n_cas} aa)",
        xy=(cas_x + CAS_W / 2, cas_y + 0.31),
        xytext=(cas_x + CAS_W / 2 + 0.25, cas_y + 0.31),
        arrowprops=arrow,
        **lab,
    )
    ax_b.text(
        (xs[2] + xs[3]) / 2,
        hi + 1.0,
        "loops untouched",
        ha="center",
        va="bottom",
        fontsize=FS,
        color=GREY,
        style="italic",
    )
    ax_b.text(
        x_tm + 0.6,
        9.25,
        f"snorkel module: {n_module} aa after the C-terminus",
        ha="center",
        va="center",
        fontsize=FS,
        color=INK,
    )
    ax_b.text(0.0, 9.6, "b", fontsize=FS + 3, fontweight="bold", va="top", color=INK)

    save_figure(fig, BASENAME, output_dir=args.out_dir)
    plt.close(fig)
    print(f"wrote {args.out_dir / (BASENAME + '.png')} (module {n_module} aa)")


if __name__ == "__main__":
    main()
