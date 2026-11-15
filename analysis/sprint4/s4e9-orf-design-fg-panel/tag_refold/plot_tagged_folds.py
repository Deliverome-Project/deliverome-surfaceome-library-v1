"""Figures for the tagged-construct refolds, from ``tag_refold_metrics.tsv``.

1. ``s4e9_tag_refold_lddt``: native-region lDDT (left) and lDDT within 15 residues of the
   insertion site (right), one point per construct, by rung, beside the untagged controls,
   with the pass thresholds.
2. ``s4e9_tag_refold_epitopes``: SpyTag003 vs ALFA relative solvent accessibility, one
   point per construct, with the buried threshold.

Implemented by Claude Code session l1-tag-refold.

    uv run --locked python analysis/sprint4/s4e9-orf-design-fg-panel/tag_refold/plot_tagged_folds.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from score_tagged_folds import EPITOPE_RSA_PASS, LDDT_PASS, LOCAL_LDDT_PASS  # noqa: E402

from deliverome_analysis.plotting import (  # noqa: E402
    CATEGORICAL_PALETTE,
    COLORS,
    save_figure,
    setup_plotting_style,
)

ORDER = [
    ("control", "untagged\ncontrols"),
    ("sp_control", "untagged\nSP controls"),
    ("nterm", "N after SP"),
    ("cterm", "C-terminal"),
    ("type3", "type III"),
    ("multipass", "res-1\nmultipass"),
    ("loop", "internal\nloop"),
    ("rescue", "expression\nrescue"),
    ("multipass_spacer", "multipass\n+spacer"),
    ("cterm_spacer", "snorkel\n+spacer"),
]


def strip(ax, df: pd.DataFrame, col: str, threshold: float, ylabel: str) -> None:
    rng = np.random.default_rng(0)
    present = [
        (k, lab)
        for k, lab in ORDER
        if (df["rung"] == k).any() and df.loc[df["rung"] == k, col].notna().any()
    ]
    for x, (k, _) in enumerate(present):
        v = df.loc[df["rung"] == k, col].dropna().to_numpy(float)
        color = (
            COLORS["neutral"]
            if "control" in k
            else CATEGORICAL_PALETTE[x % len(CATEGORICAL_PALETTE)]
        )
        ax.scatter(
            x + rng.uniform(-0.28, 0.28, len(v)), v, s=8, alpha=0.6, color=color, linewidths=0
        )
        ax.hlines(np.median(v), x - 0.35, x + 0.35, color=COLORS["dark"], lw=1.5)
        ax.text(
            x,
            1.02,
            f"n={len(v)}",
            ha="center",
            va="bottom",
            fontsize=9,
            transform=ax.get_xaxis_transform(),
        )
    ax.axhline(threshold, color=COLORS["primary"], ls="--", lw=1)
    ax.set_xticks(range(len(present)), [lab for _, lab in present], fontsize=9)
    ax.set_ylabel(ylabel)


def main() -> None:
    import matplotlib.pyplot as plt

    df = pd.read_csv(HERE.parent / "tag_refold_metrics.tsv", sep="\t")
    df = df[df["verdict"].isin(["pass", "fail", "control"]) & (df["kind"] != "superseded")]
    out = HERE.parent / "figures"
    setup_plotting_style(font_scale=1.0)

    fig, axes = plt.subplots(1, 2, figsize=(15, 5.2), gridspec_kw={"width_ratios": [1.15, 1]})
    strip(axes[0], df, "lddt_confident", LDDT_PASS, r"native-region C$\alpha$ lDDT vs native model")
    strip(axes[1], df[df["rung"].str.contains("control") == False], "lddt_site", LOCAL_LDDT_PASS,  # noqa: E712
          r"C$\alpha$ lDDT within 15 residues of the site")  # fmt: skip
    for ax in axes:
        ax.set_ylim(0.4, 1.01)
    fig.tight_layout()
    save_figure(fig, "s4e9_tag_refold_lddt", output_dir=out)

    t = df[df["verdict"].isin(["pass", "fail"])]
    fig, ax = plt.subplots(figsize=(6.5, 6))
    present = [(k, lab) for k, lab in ORDER if (t["rung"] == k).any()]
    for x, (k, lab) in enumerate(present):
        s = t[t["rung"] == k]
        ax.scatter(s["spytag003_rsa"], s["alfa_rsa"], s=10, alpha=0.6, linewidths=0,
                   color=CATEGORICAL_PALETTE[(x + 2) % len(CATEGORICAL_PALETTE)],
                   label=f"{lab.replace(chr(10), ' ')} (n={len(s)})")  # fmt: skip
    ax.axvline(EPITOPE_RSA_PASS, color=COLORS["primary"], ls="--", lw=1)
    ax.axhline(EPITOPE_RSA_PASS, color=COLORS["primary"], ls="--", lw=1)
    ax.set_xlabel("SpyTag003 mean relative SASA")
    ax.set_ylabel("ALFA mean relative SASA")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.legend(fontsize=8, frameon=False, loc="upper left")
    fig.tight_layout()
    save_figure(fig, "s4e9_tag_refold_epitopes", output_dir=out)


if __name__ == "__main__":
    main()
