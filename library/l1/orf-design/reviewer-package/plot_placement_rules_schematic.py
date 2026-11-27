# ruff: noqa: RUF001
"""Five-panel schematic of L1 tag placement rules by membrane topology.

Reads live CDS counts from ``l1_tag_paths.tsv`` and ``constructs_l1_*.tsv`` so the
figure always agrees with the committed order of record.

    uv run --locked python plot_placement_rules_schematic.py
"""

from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch, Rectangle

from deliverome_analysis.core.plotting import (
    COLORS,
    SEQUENTIAL_PALETTES,
    save_figure,
    setup_plotting_style,
)

HERE = Path(__file__).resolve().parent
S4E9 = HERE.parent
OUT_DIR = HERE / "figures"
BASENAME = "placement_rules_schematic"

INK = COLORS["dark"]
GREY = COLORS["neutral"]
PROTEIN = COLORS["primary"]
CASSETTE = COLORS["secondary"]
LEADER = COLORS["accent"]
MEMBRANE = SEQUENTIAL_PALETTES["amber"][5]
MEM_EDGE = SEQUENTIAL_PALETTES["amber"][3]
GPI_COLOR = COLORS["highlight"] if "highlight" in COLORS else "#8E44AD"


def load_counts() -> dict[str, int]:
    paths = list(csv.DictReader((S4E9 / "l1_tag_paths.tsv").open(), delimiter="\t"))
    by_p = Counter(r["placement"] for r in paths)
    nterm = list(csv.DictReader((S4E9 / "constructs_l1_nterm.tsv").open(), delimiter="\t"))
    mp = list(csv.DictReader((S4E9 / "constructs_l1_multipass.tsv").open(), delimiter="\t"))
    t3 = list(csv.DictReader((S4E9 / "constructs_l1_type3.tsv").open(), delimiter="\t"))
    ct = list(csv.DictReader((S4E9 / "constructs_l1_cterm.tsv").open(), delimiter="\t"))
    intl = list(csv.DictReader((S4E9 / "constructs_l1_internal.tsv").open(), delimiter="\t"))
    return {
        "nterm_cds": by_p["N after SP"],
        "nterm_native": sum(1 for r in nterm if r.get("arm") == "native_site"),
        "nterm_igk": sum(1 for r in nterm if r.get("arm") == "igk_swap"),
        "mp_cds": by_p["N at residue 1, multipass"],
        "mp_igk": sum(1 for r in mp if r.get("design") == "igk_tag"),
        "t3_cds": by_p["N at residue 1, single-pass (type III)"],
        "t3_alone": sum(1 for r in t3 if r.get("design") == "tag_alone"),
        "t3_igk": sum(1 for r in t3 if r.get("design") == "igk_tag"),
        "ct_cds": by_p["C-terminal"],
        "ct_built": sum(1 for r in ct if "snorkel" not in r.get("placement", "").lower()),
        "int_cds": by_p["internal"],
        "int_built_cds": len({r["design_id"] for r in intl}),
        "int_loops": len(intl),
        "snorkels": sum(1 for r in ct if "snorkel" in r.get("placement", "").lower()),
    }


def draw_helix(
    ax, x: float, y0: float = -0.78, y1: float = 0.78, w: float = 0.34, color: str = PROTEIN
) -> None:
    ax.add_patch(
        FancyBboxPatch(
            (x - w / 2, y0),
            w,
            y1 - y0,
            boxstyle="round,pad=0.03,rounding_size=0.14",
            facecolor=color,
            edgecolor="white",
            linewidth=0.9,
            zorder=4,
        )
    )


def draw_cassette(ax, x: float, y: float, w: float = 1.15, h: float = 0.40) -> None:
    ax.add_patch(
        FancyBboxPatch(
            (x - w / 2, y - h / 2),
            w,
            h,
            boxstyle="round,pad=0.02,rounding_size=0.10",
            facecolor=CASSETTE,
            edgecolor="white",
            linewidth=1.0,
            zorder=6,
        )
    )
    ax.text(
        x,
        y,
        "Tag",
        ha="center",
        va="center",
        color="white",
        fontsize=8.5,
        fontweight="bold",
        zorder=7,
    )


def draw_cleaved_leader(
    ax, x: float, y: float, label: str, w: float = 0.92, h: float = 0.34, color: str = GREY
) -> None:
    ax.add_patch(
        FancyBboxPatch(
            (x - w / 2, y - h / 2),
            w,
            h,
            boxstyle="round,pad=0.02,rounding_size=0.08",
            facecolor="white",
            edgecolor=color,
            linewidth=1.1,
            linestyle="--",
            zorder=5,
        )
    )
    ax.text(
        x,
        y,
        label,
        ha="center",
        va="center",
        color=color,
        fontsize=7.5,
        fontweight="bold",
        zorder=6,
    )


def draw_arc(
    ax, x0: float, x1: float, y_base: float, height: float, color: str = PROTEIN, lw: float = 2.0
) -> None:
    t = np.linspace(0, np.pi, 60)
    xm = 0.5 * (x0 + x1)
    xr = 0.5 * (x1 - x0)
    xs = xm - xr * np.cos(t)
    ys = y_base + height * np.sin(t)
    ax.plot(xs, ys, color=color, lw=lw, solid_capstyle="round", zorder=3)


def draw_panel_frame(ax, title: str, count_line: str, rule_line: str) -> None:
    ax.set_xlim(-2.1, 2.1)
    ax.set_ylim(-2.55, 2.75)
    ax.axis("off")
    # Lipid bilayer band
    ax.add_patch(
        Rectangle(
            (-2.05, -0.68), 4.1, 1.36, facecolor=MEMBRANE, edgecolor="none", alpha=0.42, zorder=1
        )
    )
    ax.plot([-2.05, 2.05], [0.68, 0.68], color=MEM_EDGE, lw=1.0, alpha=0.8, zorder=2)
    ax.plot([-2.05, 2.05], [-0.68, -0.68], color=MEM_EDGE, lw=1.0, alpha=0.8, zorder=2)
    # Header above diagram
    ax.text(0, 2.58, title, ha="center", va="top", fontsize=9.8, fontweight="bold", color=INK)
    ax.text(
        0, 2.18, count_line, ha="center", va="top", fontsize=8.4, color=CASSETTE, fontweight="bold"
    )
    # Concise rule caption below diagram
    ax.text(0, -1.55, rule_line, ha="center", va="top", fontsize=7.8, color=INK, linespacing=1.32)


def draw_cut_arrow(ax, x: float, color: str) -> None:
    ax.annotate(
        "",
        xy=(x, 1.44),
        xytext=(x, 1.72),
        arrowprops={"arrowstyle": "-|>", "color": color, "lw": 1.1},
    )
    ax.text(x, 1.78, "cut", ha="center", va="bottom", fontsize=7.0, color=color)


def _draw_n_after_sp(ax, c: dict[str, int]) -> None:
    draw_panel_frame(
        ax,
        "1. N after SP",
        f"{c['nterm_cds']:,} CDSs",
        f"Native SP ({c['nterm_native']:,}) or Igκ ({c['nterm_igk']})\nIncludes 186 GPI CDSs\n(1 held)",
    )
    draw_helix(ax, 0.35)
    ax.plot([0.35, 0.35, -0.35], [0.78, 1.25, 1.25], color=PROTEIN, lw=2.2, zorder=3)
    draw_cassette(ax, -0.60, 1.25, w=0.95)
    draw_cleaved_leader(ax, -1.52, 1.25, "SP", w=0.62, color=GREY)
    ax.plot([-1.18, -1.07], [1.25, 1.25], color=GREY, lw=1.2, ls=":", zorder=4)
    draw_cut_arrow(ax, -1.12, GREY)
    ax.plot([0.35, 0.35, 0.95], [-0.78, -1.18, -1.18], color=PROTEIN, lw=2.2, zorder=3)
    ax.text(
        1.12, -1.18, "C", ha="left", va="center", fontsize=8.5, fontweight="bold", color=PROTEIN
    )
    ax.text(
        -2.00,
        0.84,
        "Extracellular",
        ha="left",
        va="bottom",
        fontsize=7.5,
        color=GREY,
        style="italic",
    )
    ax.text(
        -2.00, -0.86, "Cytoplasmic", ha="left", va="top", fontsize=7.5, color=GREY, style="italic"
    )


def _draw_multipass(ax, c: dict[str, int]) -> None:
    draw_panel_frame(
        ax,
        "2. Multipass (N-out)",
        f"{c['mp_cds']:,} CDSs",
        f"Cleaved Igκ + tag ({c['mp_igk']}/{c['mp_cds']})\nAssists TM1 translocation\n& boosts GPCR expression",
    )
    xs = [-0.55, 0.15, 0.85]
    for x in xs:
        draw_helix(ax, x)
    draw_arc(ax, xs[0], xs[1], -0.78, -0.40)
    draw_arc(ax, xs[1], xs[2], 0.78, 0.42)
    ax.plot([xs[0], xs[0]], [0.78, 1.25], color=PROTEIN, lw=2.2, zorder=3)
    draw_cassette(ax, -0.62, 1.25, w=0.95)
    draw_cleaved_leader(ax, -1.52, 1.25, "Igκ", w=0.64, color=LEADER)
    ax.plot([-1.19, -1.09], [1.25, 1.25], color=LEADER, lw=1.2, ls=":", zorder=4)
    draw_cut_arrow(ax, -1.14, LEADER)
    ax.plot([xs[2], xs[2], 1.35], [-0.78, -1.18, -1.18], color=PROTEIN, lw=2.2, zorder=3)
    ax.text(
        1.45, -1.18, "C", ha="left", va="center", fontsize=8.5, fontweight="bold", color=PROTEIN
    )


def main() -> int:
    setup_plotting_style()
    c = load_counts()

    fig, axes = plt.subplots(1, 5, figsize=(13.6, 4.35), gridspec_kw={"wspace": 0.08})

    _draw_n_after_sp(axes[0], c)
    _draw_multipass(axes[1], c)

    # 3. N at residue 1, single-pass (type III)
    ax = axes[2]
    draw_panel_frame(
        ax,
        "3. Type III (1-TM N-out)",
        f"{c['t3_cds']} CDSs",
        f"Tag at Met1 ({c['t3_alone']} alone;\n{c['t3_igk']} with Igκ for short\ntail / charge-bias risk)",
    )
    draw_helix(ax, 0.25)
    ax.plot([0.25, 0.25, -0.35], [0.78, 1.25, 1.25], color=PROTEIN, lw=2.2, zorder=3)
    draw_cassette(ax, -0.55, 1.25, w=1.0)
    ax.text(
        -1.18, 1.25, "M1", ha="right", va="center", fontsize=8.0, fontweight="bold", color=PROTEIN
    )
    ax.plot([0.25, 0.25, 0.90], [-0.78, -1.18, -1.18], color=PROTEIN, lw=2.2, zorder=3)
    ax.text(
        1.02, -1.18, "C", ha="left", va="center", fontsize=8.5, fontweight="bold", color=PROTEIN
    )

    # 4. C-terminal (Type II / C-out)
    ax = axes[3]
    draw_panel_frame(
        ax,
        "4. C-terminal",
        f"{c['ct_cds']} CDSs",
        f"Tag at native C-terminus\n({c['ct_built']} built, {c['ct_cds'] - c['ct_built']} held;\n3 stalk/furin moves)",
    )
    draw_helix(ax, -0.35)
    ax.plot([-0.35, -0.35, -0.95], [-0.78, -1.18, -1.18], color=PROTEIN, lw=2.2, zorder=3)
    ax.text(
        -1.08, -1.18, "N", ha="right", va="center", fontsize=8.5, fontweight="bold", color=PROTEIN
    )
    ax.plot([-0.35, -0.35, 0.45], [0.78, 1.25, 1.25], color=PROTEIN, lw=2.2, zorder=3)
    draw_cassette(ax, 0.72, 1.25, w=1.0)

    # 5. Internal loops (+ optional snorkel)
    ax = axes[4]
    draw_panel_frame(
        ax,
        "5. Internal loops",
        f"{c['int_cds']} CDSs",
        f"1–3 ranked loops ({c['int_loops']:,} arms)\n+ {c['snorkels']} C-term snorkels\n({c['int_cds'] - c['int_built_cds']} held; incl. BST2)",
    )
    ixs = [-1.05, -0.35, 0.35, 1.05]
    for x in ixs:
        draw_helix(ax, x, w=0.30)
    draw_arc(ax, ixs[0], ixs[1], 0.78, 0.62)
    draw_arc(ax, ixs[1], ixs[2], -0.78, -0.38)
    draw_arc(ax, ixs[2], ixs[3], 0.78, 0.42)
    draw_cassette(ax, -0.70, 1.44, w=0.92, h=0.38)
    ax.plot([ixs[0], ixs[0], -1.45], [-0.78, -1.15, -1.15], color=PROTEIN, lw=2.0, zorder=3)
    ax.text(
        -1.55, -1.15, "N", ha="right", va="center", fontsize=8.2, fontweight="bold", color=PROTEIN
    )
    ax.plot([ixs[3], ixs[3], 1.45], [-0.78, -1.15, -1.15], color=PROTEIN, lw=2.0, zorder=3)
    ax.text(
        1.55, -1.15, "C", ha="left", va="center", fontsize=8.2, fontweight="bold", color=PROTEIN
    )

    # Legend strip at bottom
    fig.text(
        0.5,
        0.025,
        "First matching rule assigns each CDS (DeepTMHMM v1 on exact built sequence, plus 33 structural/experimental overrides).  "
        "Teal box = 51-aa SpyTag003–ALFA cassette;  dashed box = cleaved signal peptide / Igκ leader.",
        ha="center",
        va="bottom",
        fontsize=8.2,
        color=GREY,
    )
    fig.subplots_adjust(left=0.015, right=0.985, top=0.95, bottom=0.10)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    save_figure(fig, OUT_DIR / BASENAME)
    plt.close(fig)
    print(f"wrote {OUT_DIR / (BASENAME + '.png')} ({c})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
