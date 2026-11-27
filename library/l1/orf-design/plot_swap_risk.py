"""Every construct pair in the order, against Walton's barcode-swap curve.

The point of the figure is the separation: after per-construct cassette diversification,
no two *different* members are at risk, and every pair over the bar is two constructs of
one member — where the shared tract is the member's own ORF and cannot be recoded apart.

One point per pair, so all 210 are shown rather than a summary.

    uv run --locked python plot_swap_risk.py
"""

from __future__ import annotations

import csv
import itertools
from pathlib import Path

import numpy as np

from deliverome_analysis.cloning import seqtools
from deliverome_analysis.core.plotting import (
    COLORS,
    create_figure,
    get_color,
    save_figure,
    setup_plotting_style,
)

HERE = Path(__file__).parent
FASTAS = (
    "constructs_nterm_disorder_anchored.fasta",
    "constructs_loop.fasta",
    "constructs_tag_arms.fasta",
    "constructs_cterm.fasta",
    "constructs_gpcr_rescue.fasta",
)
BAR_BP = 100
PLACEHOLDER = "N" * 20


def load() -> dict[str, str]:
    out: dict[str, str] = {}
    for fname in FASTAS:
        path = HERE / fname
        if not path.exists():
            continue
        name = None
        for line in path.read_text().splitlines():
            if line.startswith(">"):
                name = line[1:].split("|")[0]
                out[name] = ""
            elif name:
                out[name] += line.strip()
    # The barcode placeholder is not sequence; leaving it in would read as identity.
    return {k: v.replace(PLACEHOLDER, "") for k, v in out.items()}


def main() -> int:
    setup_plotting_style()
    cons = load()
    pairs = []
    for a, b in itertools.combinations(sorted(cons), 2):
        bp = seqtools.longest_shared_tract(cons[a], cons[b])[0]
        pairs.append((bp, a, b, a.split("_")[0] == b.split("_")[0]))

    same = [p for p in pairs if p[3]]
    cross = [p for p in pairs if not p[3]]

    fig, ax = create_figure(figsize=(9.6, 5.6))
    grid = np.linspace(1, 2600, 600)
    ax.plot(
        grid,
        [seqtools.barcode_swap_rate(x) * 100 for x in grid],
        lw=1.6,
        color=COLORS["neutral"],
        zorder=1,
        label="Walton fit  Y = 0.5(1 - exp(-2kx)), k = 5.3e-4",
    )
    ax.axvline(
        BAR_BP,
        color=COLORS["danger"],
        lw=1.3,
        ls=":",
        zorder=2,
        label=f"bar {BAR_BP} bp = {seqtools.barcode_swap_rate(BAR_BP):.1%}",
    )

    rng = np.random.default_rng(0)
    cx = np.array([p[0] for p in cross], dtype=float)
    ax.scatter(
        cx + rng.uniform(-6, 6, len(cx)),
        [seqtools.barcode_swap_rate(v) * 100 for v in cx],
        s=26,
        alpha=0.55,
        linewidths=0,
        color=COLORS["secondary"],
        zorder=3,
        label=f"different members (n={len(cross)}, all under the bar)",
    )
    for bp, _a, _b, _ in same:
        ax.scatter(
            [bp],
            [seqtools.barcode_swap_rate(bp) * 100],
            s=118,
            color=get_color("primary"),
            edgecolors=COLORS["dark"],
            linewidths=0.9,
            zorder=5,
        )
    ax.scatter(
        [],
        [],
        s=118,
        color=get_color("primary"),
        edgecolors=COLORS["dark"],
        linewidths=0.9,
        label=f"two constructs of one member (n={len(same)})",
    )

    labels = {
        2373: ("TFRC linker arms\n(3 pairs)", -12, 14, "right"),
        1490: ("CNR1 plain vs igk", -12, 2, "right"),
        1130: ("CCR5", 10, -4, "left"),
        1034: ("OR51E2", 4, -24, "left"),
        762: ("SLC2A4 loop pair", 12, -6, "left"),
        432: ("MS4A1 loop pair", 12, -8, "left"),
    }
    for bp, (text, dx, dy, ha) in labels.items():
        ax.annotate(
            text,
            (bp, seqtools.barcode_swap_rate(bp) * 100),
            textcoords="offset points",
            xytext=(dx, dy),
            ha=ha,
            fontsize=8.4,
            color=COLORS["dark"],
        )

    ax.set_xlabel("longest tract the pair shares (bp)")
    ax.set_ylabel("predicted barcode/ORF swapping (%)")
    ax.set_xlim(-40, 2700)
    ax.set_ylim(-2, 52)
    ax.legend(loc="lower right", fontsize=8.0, frameon=False)
    fig.tight_layout()
    save_figure(fig, "s4e9_swap_risk_pairs", output_dir=HERE / "figures")

    with (HERE / "swap_risk_pairs.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["construct_a", "construct_b", "same_member", "shared_bp", "swap_rate"])
        for bp, a, b, is_same in sorted(pairs, reverse=True):
            w.writerow([a, b, is_same, bp, f"{seqtools.barcode_swap_rate(bp):.4f}"])
    over = sum(1 for p in cross if p[0] > BAR_BP)
    print(
        f"{len(pairs)} pairs: {len(same)} same-member, {len(cross)} cross-gene "
        f"({over} over the {BAR_BP} bp bar; max {max(p[0] for p in cross)} bp)"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
