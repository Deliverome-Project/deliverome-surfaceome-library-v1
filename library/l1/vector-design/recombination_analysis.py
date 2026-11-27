"""s1e2 — lentiviral template-switching (ORF↔barcode decoupling) risk vs span length.

Applies the Walton, Qin & Blainey (CROPseq-multi, bioRxiv 2024.03.17.585235) rate
model to our construct:

    P(>=1 switch over L bp) = 1 - (1 - p)^L

with p = 5.3e-4 per bp (point estimate), reported range 4.9e-4 - 13.5e-4. Their model
is fit as a function of the length of *homologous* intervening sequence, which is our
case: the linker/tag scaffold between the variable region and the barcode is constant
across library members, hence 100% homologous between two co-packaged genomes.

The figure shows the point estimate with the reported rate range as a shaded band
(i.e. the "error bars" on the published rate), marks our candidate spans, and marks
the 5% threshold.

Also reports barcode-collision risk at the real target library size (30,000).

Run:
    uv run python library/l1/vector-design/recombination_analysis.py
"""

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

from deliverome_analysis.core.plotting import COLORS, save_figure, setup_plotting_style

HERE = Path(__file__).resolve().parent

# Walton, Qin & Blainey per-bp template-switching probability.
P_POINT = 5.3e-4
P_LOW = 4.9e-4
P_HIGH = 13.5e-4

ACCEPTABLE = 0.05  # 5% threshold used throughout the doc

# Candidate spans in our construct.
SPANS = {
    "stop → barcode (naive)": 21,
    "linker start → barcode (CORRECT)": 177,
    "+ snorkel +60 bp": 237,
    "+ snorkel +120 bp": 297,
}

TARGET_LIBRARY = 30_000  # confirmed target library size


def p_switch(length, p):
    return 1.0 - (1.0 - p) ** length


def threshold_length(p, target=ACCEPTABLE):
    return np.log(1 - target) / np.log(1 - p)


def main() -> None:
    setup_plotting_style()

    lengths = np.arange(1, 401)
    point = p_switch(lengths, P_POINT) * 100
    low = p_switch(lengths, P_LOW) * 100
    high = p_switch(lengths, P_HIGH) * 100

    fig, ax = plt.subplots(figsize=(8, 5.5))

    # mathtext for the exponents — the bundled Manrope lacks U+207B (superscript minus).
    ax.fill_between(
        lengths,
        low,
        high,
        color=COLORS["primary"],
        alpha=0.15,
        linewidth=0,
        label=r"reported rate range (4.9-13.5 $\times$ 10$^{-4}$/bp)",
    )
    ax.plot(
        lengths,
        point,
        color=COLORS["primary"],
        lw=2.5,
        label=r"point estimate (5.3 $\times$ 10$^{-4}$/bp)",
    )

    ax.axhline(ACCEPTABLE * 100, ls="--", color="0.35", lw=1.5, zorder=0)
    ax.annotate(
        "5% threshold",
        xy=(398, ACCEPTABLE * 100),
        xytext=(0, 5),
        textcoords="offset points",
        ha="right",
        fontsize="small",
        color="0.35",
    )

    # Mark our candidate spans. Numbered markers + a compact side table keeps the
    # 177/237/297 labels from colliding.
    legend_lines = []
    for i, (label, L) in enumerate(SPANS.items(), start=1):
        y = p_switch(L, P_POINT) * 100
        ylo = p_switch(L, P_LOW) * 100
        yhi = p_switch(L, P_HIGH) * 100
        ax.errorbar(
            L,
            y,
            yerr=[[y - ylo], [yhi - y]],
            fmt="o",
            color=COLORS["secondary"],
            markersize=10,
            capsize=5,
            capthick=2,
            lw=2,
            zorder=5,
        )
        ax.annotate(
            str(i),
            xy=(L, y),
            xytext=(0, 0),
            textcoords="offset points",
            ha="center",
            va="center",
            fontsize="xx-small",
            color="white",
            weight="bold",
            zorder=6,
        )
        legend_lines.append(f"{L} bp → {y:.1f}%  ·  {label}")

    ax.set_xlabel("homologous span between variable region and barcode (bp)")
    ax.set_ylabel("P(≥1 template switch) (%)")
    ax.set_xlim(0, 400)
    ax.set_ylim(0, 42)

    # Fold the span annotations into the legend so nothing can overlap the error bars.
    handles, labels = ax.get_legend_handles_labels()
    for line in legend_lines:
        handles.append(
            Line2D(
                [],
                [],
                marker="o",
                color=COLORS["secondary"],
                linestyle="none",
                markersize=7,
            )
        )
        labels.append(line)
    ax.legend(
        handles,
        labels,
        frameon=False,
        fontsize="x-small",
        loc="upper left",
        handletextpad=0.7,
        labelspacing=0.75,
    )

    save_figure(fig, "s1e2_recombination_risk", output_dir=str(HERE / "figures"))

    # ---- table out ----
    rows = []
    for label, L in SPANS.items():
        rows.append(
            {
                "span": label.replace("\n", " "),
                "bp": L,
                "p_switch_low_4.9e-4": f"{p_switch(L, P_LOW) * 100:.2f}%",
                "p_switch_point_5.3e-4": f"{p_switch(L, P_POINT) * 100:.2f}%",
                "p_switch_high_13.5e-4": f"{p_switch(L, P_HIGH) * 100:.2f}%",
            }
        )
    out = HERE / "recombination_risk_by_span.tsv"
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {out}")

    print("\n--- 5% threshold length ---")
    for name, p in [("low 4.9e-4", P_LOW), ("point 5.3e-4", P_POINT), ("high 13.5e-4", P_HIGH)]:
        print(f"  {name}: {threshold_length(p):.0f} bp")

    print(f"\n--- expected mispaired clones at target library size N={TARGET_LIBRARY:,} ---")
    for label, L in SPANS.items():
        f = p_switch(L, P_POINT)
        print(
            f"  {label.replace(chr(10), ' '):<38} {L:>4} bp: "
            f"{f * 100:>5.2f}%  ≈ {f * TARGET_LIBRARY:>8,.0f} of {TARGET_LIBRARY:,} clones affected"
        )

    # Barcode collision at the real library size (20-nt barcode, effective space from
    # barcode_length_analysis.py).
    eff_space_20nt = 6.541e11
    frac = 1 - (1 - 1 / eff_space_20nt) ** (TARGET_LIBRARY - 1)
    print(f"\n--- 20-nt barcode collision at N={TARGET_LIBRARY:,} ---")
    print(f"  effective space ≈ {eff_space_20nt:.3e}")
    print(f"  expected collision fraction ≈ {frac:.2e}  (~{frac * TARGET_LIBRARY:.4f} clones)")
    print("  => negligible; 20 nt is far more than sufficient for a 30k library.")


if __name__ == "__main__":
    main()
