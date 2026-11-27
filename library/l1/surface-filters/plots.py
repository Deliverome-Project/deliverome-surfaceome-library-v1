#!/usr/bin/env python3
"""Plot isoform + ortholog outcome counts from the candidate TSVs.

Run:  uv run python plots.py
"""

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from deliverome_analysis.core.plotting import CATEGORICAL_PALETTE, save_figure, setup_plotting_style

HERE = os.path.dirname(os.path.abspath(__file__))
iso = pd.read_csv(os.path.join(HERE, "isoform_candidates.tsv"), sep="\t")
orth = pd.read_csv(os.path.join(HERE, "ortholog_candidates.tsv"), sep="\t")

setup_plotting_style()


def label_bars(ax, xs, ys):
    for x, y in zip(xs, ys, strict=False):
        ax.text(x, y, f"{int(y)}", ha="center", va="bottom")


# ---- isoforms: topology status + recommendation ----
topo_order = ["surface", "secreted", "intracellular"]
rec_order = ["include", "flag", "drop"]
topo = iso["topology_status"].value_counts().reindex(topo_order).fillna(0)
rec = iso["recommendation"].value_counts().reindex(rec_order).fillna(0)

fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.2))
c = CATEGORICAL_PALETTE
a1.bar(topo_order, topo.values, color=[c[1], c[2], c[0]])
a1.set_ylabel("isoforms")
a1.set_xlabel("topology status  (inclusion gate)")
label_bars(a1, range(len(topo_order)), topo.values)

a2.bar(rec_order, rec.values, color=[c[1], c[3], c[0]])
a2.set_ylabel("isoforms")
a2.set_xlabel(
    "recommendation  (include = surface & same topology · "
    "flag = surface, TM changed · drop = secreted/intracellular)"
)
label_bars(a2, range(len(rec_order)), rec.values)
save_figure(fig, "isoform_outcomes", output_dir=os.path.join(HERE, "figures"))

# ---- orthologs: recommendation grouped by species ----
oorder = ["surrogate", "evaluate", "drop"]
piv = (
    orth.pivot_table(
        index="recommendation", columns="species", values="gene", aggfunc="count", fill_value=0
    )
    .reindex(oorder)
    .fillna(0)
)
mouse = piv["mouse"].values if "mouse" in piv else np.zeros(len(oorder))
cyno = piv["cyno"].values if "cyno" in piv else np.zeros(len(oorder))

fig2, ax = plt.subplots(figsize=(7.5, 4.2))
x = np.arange(len(oorder))
w = 0.38
ax.bar(x - w / 2, mouse, w, label="mouse", color=c[0])
ax.bar(x + w / 2, cyno, w, label="cyno", color=c[1])
ax.set_xticks(x)
ax.set_xticklabels(oorder)
ax.set_ylabel("orthologs")
ax.set_xlabel(
    "recommendation  (surrogate = faithful & topology-matched · "
    "evaluate = diverges · drop = not surface)"
)
ax.legend()
label_bars(ax, x - w / 2, mouse)
label_bars(ax, x + w / 2, cyno)
save_figure(fig2, "ortholog_outcomes", output_dir=os.path.join(HERE, "figures"))

print("isoform topology:", topo.to_dict())
print("isoform recommendation:", rec.to_dict())
print("ortholog by species:\n", piv)
print("wrote isoform_outcomes + ortholog_outcomes figures")
