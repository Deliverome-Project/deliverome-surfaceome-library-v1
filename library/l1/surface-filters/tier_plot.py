#!/usr/bin/env python3
"""Grouped bar chart: T1/T2/T3 with hues for genes / surface isoforms / surface
orthologs. Run: uv run python tier_plot.py
"""

import os

import final_library as F  # reuses the final cut + per-tier aggregation (F.agg)
import matplotlib.pyplot as plt
import numpy as np

from deliverome_analysis.core.plotting import CATEGORICAL_PALETTE, save_figure, setup_plotting_style

HERE = os.path.dirname(os.path.abspath(__file__))
setup_plotting_style()

tiers = ["T1", "T2", "T3", "LL", "SV", "OE"]
tier_labels = [
    "T1\nCore",
    "T2\nExtended",
    "T3\nFlagged",
    "LL\nLow-lit",
    "SV\nScreen-val",
    "OE\nOverexpr",
]
# olfactory/taste footnote: how T3's chaperone-dependent GPCRs are identified
OLF_NOTE = (
    "T3 includes olfactory / taste / vomeronasal GPCRs (UniProt family or HGNC gene-group "
    "name contains 'olfactory'/'taste'/'vomeronasal') — RTP-chaperone-dependent, so flagged."
)
series = [("genes", "Genes"), ("iso", "Surface isoforms"), ("orth", "Surface orthologs")]

x = np.arange(len(tiers))
w = 0.26
# ---- figure 1: grouped bars per tier (hue = category) + overarching "All tiers" total ----
groups = [*tiers, "OV", "All"]
group_labels = [*tier_labels, "OV\nOversized", "All\ntiers"]


def _gval(g, key):
    return sum(F.agg[t][key] for t in tiers) if g == "All" else F.agg[g][key]


xg = np.arange(len(groups))
w = 0.27
fig, ax = plt.subplots(figsize=(15, 5.4))
for i, (key, label) in enumerate(series):
    vals = [_gval(g, key) for g in groups]
    xs = xg + (i - 1) * w
    ax.bar(xs, vals, w, label=label, color=CATEGORICAL_PALETTE[i])
    for xx, v in zip(xs, vals, strict=False):
        ax.text(xx, v, str(v), ha="center", va="bottom", fontsize=7, rotation=90)
ax.axvline(
    len(tiers) - 0.5, color="0.6", linestyle="--", linewidth=1
)  # separate the overarching total
ax.set_xticks(xg)
ax.set_xticklabels(group_labels)
ax.set_ylabel("count")
ax.set_xlabel("tier (confidence of autonomous surface presentation)")
ax.legend(title="")
fig.text(0.5, -0.02, OLF_NOTE, ha="center", va="top", fontsize=7.5, color="0.4", wrap=True)
save_figure(fig, "tier_gene_isoform_ortholog", output_dir=os.path.join(HERE, "figures"))

# second view: genes with >=1 surface isoform / ortholog (coverage, not totals)
fig2, ax2 = plt.subplots(figsize=(8.5, 4.8))
series2 = [
    ("genes", "Genes"),
    ("iso_genes", "Genes w/ ≥1 surf isoform"),
    ("orth_genes", "Genes w/ ≥1 surf ortholog"),
]
for i, (key, label) in enumerate(series2):
    vals = [F.agg[t][key] for t in tiers]
    xs = x + (i - 1) * w
    ax2.bar(xs, vals, w, label=label, color=CATEGORICAL_PALETTE[i])
    for xx, v in zip(xs, vals, strict=False):
        ax2.text(xx, v, str(v), ha="center", va="bottom")
ax2.set_xticks(x)
ax2.set_xticklabels(tier_labels)
ax2.set_ylabel("genes")
ax2.set_xlabel("tier")
ax2.legend(title="")
save_figure(fig2, "tier_variant_coverage", output_dir=os.path.join(HERE, "figures"))

# third view: x = category, bars hued by tier (T1/T2/T3) + an overarching "All tiers" total bar
cats = [("genes", "Genes"), ("iso", "Surface isoforms"), ("orth", "Surface orthologs")]
bars = [
    ("T1", "T1 — Core", CATEGORICAL_PALETTE[0]),
    ("T2", "T2 — Extended", CATEGORICAL_PALETTE[1]),
    ("T3", "T3 — Flagged", CATEGORICAL_PALETTE[2]),
    ("LL", "LL — Low-lit", CATEGORICAL_PALETTE[3]),
    ("SV", "SV — Screen-val", CATEGORICAL_PALETTE[4]),
    ("OE", "OE — Overexpr-rescue", CATEGORICAL_PALETTE[5]),
    ("OV", "OV — Oversized (>5kb, sep.)", "0.6"),
    ("All", "All tiers (total, excl. OV)", "0.35"),
]
xc = np.arange(len(cats))
wc = 0.1
fig3, ax3 = plt.subplots(figsize=(15, 5.6))
for i, (t, label, color) in enumerate(bars):
    vals = [sum(F.agg[tt][k] for tt in tiers) if t == "All" else F.agg[t][k] for k, _ in cats]
    xs = xc + (i - 3.5) * wc
    ax3.bar(xs, vals, wc, label=label, color=color)
    for xx, v in zip(xs, vals, strict=False):
        ax3.text(xx, v, str(v), ha="center", va="bottom", fontsize=6.5, rotation=90)
ax3.set_xticks(xc)
ax3.set_xticklabels([lbl for _, lbl in cats])
ax3.set_ylabel("count")
# legend above the axes so it never clips the tall "All tiers" bars
ax3.legend(title="", ncol=4, loc="lower center", bbox_to_anchor=(0.5, 1.01), frameon=False)
fig3.text(0.5, -0.02, OLF_NOTE, ha="center", va="top", fontsize=7.5, color="0.4", wrap=True)
save_figure(fig3, "category_by_tier_with_total", output_dir=os.path.join(HERE, "figures"))
print("wrote tier_gene_isoform_ortholog + tier_variant_coverage + category_totals_by_tier")
