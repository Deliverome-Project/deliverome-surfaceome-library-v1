"""Why these three GPCRs, and what moving the two loop sites bought.

Panel A places the five panel GPCRs in the distribution of Tedman's whole screen, so the
selection argument is visible rather than asserted: CNR1 and CCR5 sit in the weak tail,
and OR51E2 does not. Panel B shows the barcode-swap tract for the two same-gene loop
pairs before and after the site move.

Every point in A is one receptor construct, so the screen is drawn as the individual
observations it is rather than as a summary.

    uv run --locked python plot_rescue_selection.py
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import pandas as pd

from deliverome_analysis.plotting import (
    COLORS,
    create_figure,
    get_color,
    save_figure,
    setup_plotting_style,
)
from deliverome_analysis.s3 import external

HERE = Path(__file__).parent
SCREEN_KEY = "tedman-gpcr-surface-screen/2025-09-19/tedman2025_gpcr_screen_media-2.xlsx"
ARMS = ("CNR1", "CCR5", "OR51E2")
PANEL = ("CNR1", "CCR5", "OR51E2", "ACKR3", "AVPR2")

#: Same-gene loop pair, shared tract before and after the 2026-10-13 site move.
#: Measured with seqtools.longest_shared_tract on the assembled inserts.
SWAP = [
    ("SLC2A4\nP66 vs E193 → G320", 972, 762),
    ("MS4A1\nY161 vs N171 → P79", 505, 432),
]
BAR_BP = 100  # check_swap_risk.py's bar


def ordinal(n: int) -> str:
    """1st/2nd/3rd/4th, with the 11-13 exception that naive suffixing gets wrong."""
    if 11 <= n % 100 <= 13:
        return f"{n}th"
    return f"{n}{ {1: 'st', 2: 'nd', 3: 'rd'}.get(n % 10, 'th') }".replace(" ", "")


def swap_rate(bp: int) -> float:
    """Walton's fit: Y(x) = 1/2 (1 - e^(-2kx)), k = 5.3e-4 bp^-1."""
    return 0.5 * (1.0 - np.exp(-2 * 5.3e-4 * bp))


def main() -> int:
    setup_plotting_style()
    screen = pd.read_excel(external.download(SCREEN_KEY), sheet_name="Canonical")
    au = pd.to_numeric(screen["Immunostaining Intensity"], errors="coerce")
    screen = screen.assign(au=au).dropna(subset=["au"])
    ours = {r.gene_id: float(r.au) for r in screen[screen["gene_id"].isin(PANEL)].itertuples()}

    fig, (ax, bx) = create_figure(1, 2, figsize=(13.0, 5.2))

    # --- A: the screen as an ECDF, so a percentile is read off directly ------
    vals = np.sort(screen["au"].to_numpy())
    frac = np.arange(1, len(vals) + 1) / len(vals) * 100.0
    ax.step(
        vals,
        frac,
        where="post",
        lw=1.8,
        color=COLORS["neutral"],
        label=f"Tedman screen (n={len(vals)} receptors)",
    )
    med = float(np.median(vals))
    ax.axhline(50, color=COLORS["line"], lw=1.0, zorder=0)
    ax.axvline(med, color=COLORS["dark"], lw=1.0, ls="--", label=f"screen median {med:.0f} a.u.")
    offsets = {
        "CNR1": (16, -34),
        "CCR5": (16, 14),
        "OR51E2": (14, -36),
        "ACKR3": (-10, 22),
        "AVPR2": (-46, -40),
    }
    for gene in PANEL:
        v = ours[gene]
        pct = float((vals < v).mean() * 100.0)
        picked = gene in ARMS
        ax.scatter(
            [v],
            [pct],
            s=132 if picked else 74,
            color=get_color("primary") if picked else COLORS["secondary"],
            edgecolors=COLORS["dark"],
            linewidths=0.9,
            zorder=6,
            marker="o" if picked else "s",
        )
        ax.annotate(
            f"{gene}  {v:.0f} a.u.\n{ordinal(round(pct))} pct",
            (v, pct),
            textcoords="offset points",
            xytext=offsets[gene],
            fontsize=8.4,
            zorder=7,
            color=COLORS["dark"] if picked else COLORS["neutral"],
        )
    ax.scatter(
        [],
        [],
        s=132,
        color=get_color("primary"),
        edgecolors=COLORS["dark"],
        linewidths=0.9,
        label="rescue arm built",
    )
    ax.scatter(
        [],
        [],
        s=74,
        marker="s",
        color=COLORS["secondary"],
        edgecolors=COLORS["dark"],
        linewidths=0.9,
        label="panel, no arm",
    )
    ax.set_xscale("log")
    ax.set_xlabel("surface immunostaining (a.u.)")
    ax.set_ylabel("percentile of the screen")
    ax.set_ylim(-4, 108)
    ax.legend(loc="upper left", fontsize=7.8, frameon=False)

    # --- B: swap tract before vs after --------------------------------------
    for i, (_label, before, after) in enumerate(SWAP):
        bx.annotate(
            "",
            xy=(after, i),
            xytext=(before, i),
            arrowprops={
                "arrowstyle": "-|>",
                "color": COLORS["line"],
                "lw": 2.6,
                "shrinkA": 6,
                "shrinkB": 7,
            },
        )
        bx.scatter(
            [before],
            [i],
            s=96,
            color=COLORS["neutral"],
            edgecolors=COLORS["dark"],
            linewidths=0.9,
            zorder=4,
        )
        bx.scatter(
            [after],
            [i],
            s=124,
            color=get_color("primary"),
            edgecolors=COLORS["dark"],
            linewidths=0.9,
            zorder=5,
        )
        bx.annotate(
            f"{before} bp / {swap_rate(before):.1%}",
            (before, i),
            textcoords="offset points",
            xytext=(6, 16),
            ha="left",
            fontsize=8.2,
            color=COLORS["neutral"],
        )
        bx.annotate(
            f"{after} bp / {swap_rate(after):.1%}",
            (after, i),
            textcoords="offset points",
            xytext=(-6, -24),
            ha="right",
            fontsize=8.2,
            color=COLORS["dark"],
        )
    bx.axvline(
        BAR_BP,
        color=COLORS["danger"],
        lw=1.2,
        ls=":",
        label=f"bar {BAR_BP} bp = {swap_rate(BAR_BP):.1%}",
    )
    bx.scatter(
        [],
        [],
        s=96,
        color=COLORS["neutral"],
        edgecolors=COLORS["dark"],
        linewidths=0.9,
        label="before the site move",
    )
    bx.scatter(
        [],
        [],
        s=124,
        color=get_color("primary"),
        edgecolors=COLORS["dark"],
        linewidths=0.9,
        label="after",
    )
    bx.set_yticks(np.arange(len(SWAP)))
    bx.set_yticklabels([s.replace("\n", "  ") for s, _, _ in SWAP], fontsize=8.4)
    bx.set_ylim(-0.75, len(SWAP) - 0.25)
    bx.set_xlim(0, 1180)
    bx.set_xlabel("longest tract shared by the pair (bp)")
    bx.legend(loc="center right", fontsize=7.8, frameon=False)

    fig.tight_layout()
    save_figure(fig, "s4e9_rescue_selection", output_dir=HERE / "figures")

    # Tables live at the experiment root, not under figures/ — that is the documented
    # layout, and .gitignore only un-ignores *.tsv at the root, so a table written
    # beside the PNG is silently untracked.
    with (HERE / "rescue_selection_values.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["panel", "label", "metric", "value"])
        for gene in PANEL:
            w.writerow(["A", gene, "immunostaining_au", f"{ours[gene]:.1f}"])
        w.writerow(["A", "screen", "median_au", f"{med:.1f}"])
        w.writerow(["A", "screen", "n_receptors", len(screen)])
        for label, before, after in SWAP:
            flat = label.replace("\n", " ")
            w.writerow(["B", flat, "shared_tract_bp_before", before])
            w.writerow(["B", flat, "shared_tract_bp_after", after])
    print("wrote figures/s4e9_rescue_selection.png and rescue_selection_values.tsv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
