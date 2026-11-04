"""Does swapping the AlphaFold model for a deposited one change which sites pass?

Plots each ordered site on the two axes the `surface_loop` lane actually gates on --
contiguous coil-run length and window-max RSA -- once from the AlphaFold model and once
from the deposited structure, with an arrow between. Signals for both come from
`experimental_signals.py`, which runs the same pydssp and the same freesasa over the two
inputs so only the coordinates differ.

The shaded corner is the region that passes (`MIN_LOOP_LEN` 8, `RSA_MIN` 0.40). A site
whose arrow leaves the corner is one the predicted structure let through and the
experimental structure would not.

    uv run --locked python plot_experimental_vs_alphafold.py
"""

from __future__ import annotations

import csv
from pathlib import Path

from matplotlib.patches import Rectangle

from deliverome_analysis.plotting import (
    COLORS,
    create_figure,
    get_color,
    save_figure,
    setup_plotting_style,
)

HERE = Path(__file__).parent
TABLE = HERE / "experimental_vs_alphafold.tsv"
MIN_LOOP_LEN, RSA_MIN = 8, 0.40

#: label offsets, hand-tuned so nothing collides
NUDGE = {
    "P79 (p0032)": (8, -20),
    "Y161 (p0033)": (10, 14),
    "N171 (retired)": (12, -6),
    "P66 (p0030)": (10, 8),
    "E193 (retired)": (6, 10),
    "G320 (p0031)": (8, -20),
}


def main() -> int:
    rows = [r for r in csv.DictReader(TABLE.open(), delimiter="\t") if r["residue"]]
    setup_plotting_style()
    fig, ax = create_figure(figsize=(9.6, 6.0))

    # Only the upper-right corner passes; shading the whole column would imply a
    # site with a long coil run but no exposure also passes, which it does not.
    ax.add_patch(Rectangle((MIN_LOOP_LEN, RSA_MIN), 30, 1.2, facecolor=COLORS["light"], zorder=0))
    ax.axvline(MIN_LOOP_LEN, color=COLORS["neutral"], lw=1.1, ls="--", zorder=1)
    ax.axhline(RSA_MIN, color=COLORS["neutral"], lw=1.1, ls="--", zorder=1)
    ax.annotate(
        f"passes:  coil run ≥ {MIN_LOOP_LEN}   and   window RSA ≥ {RSA_MIN}",
        (29.4, 1.02),
        ha="right",
        fontsize=8.6,
        color=COLORS["neutral"],
    )

    for r in rows:
        x0, y0 = int(r["af_loop_len"]), float(r["af_window_rsa"])
        x1, y1 = int(r["exp_loop_len"]), float(r["exp_window_rsa"])
        moved = (x0, y0) != (x1, y1)
        if moved:
            ax.annotate(
                "",
                xy=(x1, y1),
                xytext=(x0, y0),
                arrowprops={
                    "arrowstyle": "-|>",
                    "color": COLORS["line"],
                    "lw": 2.0,
                    "shrinkA": 7,
                    "shrinkB": 8,
                },
            )
        ax.scatter(
            [x0],
            [y0],
            s=96,
            facecolors="none",
            linewidths=1.6,
            edgecolors=COLORS["neutral"],
            zorder=4,
        )
        flips = r["af_passes"] == "True" and r["exp_passes"] != "True"
        ax.scatter(
            [x1],
            [y1],
            s=132,
            zorder=5,
            linewidths=0.9,
            edgecolors=COLORS["dark"],
            color=get_color("primary") if flips else COLORS["secondary"],
        )
        dx, dy = NUDGE.get(r["site"], (8, 8))
        ax.annotate(
            f"{r['gene']} {r['site']}",
            (x1, y1),
            textcoords="offset points",
            xytext=(dx, dy),
            fontsize=8.4,
            color=COLORS["dark"],
        )

    ax.scatter(
        [],
        [],
        s=96,
        facecolors="none",
        edgecolors=COLORS["neutral"],
        linewidths=1.6,
        label="from the AlphaFold model",
    )
    ax.scatter(
        [],
        [],
        s=132,
        color=COLORS["secondary"],
        edgecolors=COLORS["dark"],
        linewidths=0.9,
        label="from the deposited structure",
    )
    ax.scatter(
        [],
        [],
        s=132,
        color=get_color("primary"),
        edgecolors=COLORS["dark"],
        linewidths=0.9,
        label="passed on AlphaFold, fails on the structure",
    )
    ax.set_xlim(-1.6, 30)
    ax.set_ylim(0, 1.1)
    ax.set_xlabel("contiguous coil run containing the site (residues)")
    ax.set_ylabel("window-max relative solvent accessibility")
    ax.legend(loc="lower right", fontsize=8.2, frameon=False)
    fig.tight_layout()
    save_figure(fig, "s4e9_experimental_vs_alphafold", output_dir=HERE / "figures")
    print("wrote figures/s4e9_experimental_vs_alphafold.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
