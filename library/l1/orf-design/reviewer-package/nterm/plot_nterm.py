"""The N-terminal brief's two figures, from the committed outputs of ``nterm_placement.py``.

    uv run --locked python plot_nterm.py

1. ``nterm_room_before_fold``: per L1 "N after SP" member with an interval, how many
   insertion sites lie strictly between the cleavage site and the predicted fold onset.
2. ``nterm_reprediction_by_offset``: per tagged design, what SignalP 6.0 says about the
   native signal peptide, by how far past the cleavage site the cassette sits.
"""

from __future__ import annotations

import collections
import csv
from pathlib import Path

import numpy as np
from matplotlib.patches import Patch

from deliverome_analysis.plotting import create_figure, get_color, save_figure, setup_plotting_style

HERE = Path(__file__).resolve().parent
FIGURES = HERE / "figures"
VERDICTS = [
    ("preserved", "preserved", "secondary"),
    ("cleavage_moved", "cleavage moved, still before the cassette", "accent"),
    ("cleavage_crosses_junction", "cleavage at or past the junction", "quaternary"),
    ("signal_peptide_lost", "signal peptide lost", "primary"),
]


def _tsv(path: Path) -> list[dict]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main() -> int:
    setup_plotting_style()
    members = _tsv(HERE / "nterm_members.tsv")
    positions = _tsv(HERE / "nterm_positions.tsv")

    # Figure 1: strict room, 0..20 (20 is also the candidate cap, so it collects the tail).
    placed = [m for m in members if m["class"] in ("placeable", "degenerate")]
    room = collections.Counter(int(m["n_candidates_strict"]) for m in placed)
    x = np.arange(0, 21)
    vals = np.array([room.get(i, 0) for i in x])
    fig, ax = create_figure(figsize=(12, 6))
    colors = [get_color("primary") if i == 0 else get_color("secondary") for i in x]
    ax.bar(x, vals, color=colors, width=0.8)
    for xi, v in zip(x, vals, strict=True):
        if v:
            ax.text(xi, v + 6, f"{v}", ha="center", va="bottom", fontsize=12)
    ax.set_xticks([0, 5, 10, 15, 20], ["0", "5", "10", "15", "20 (cap)"])
    ax.set_xlabel("insertion sites between cleavage and predicted fold onset")
    ax.set_ylabel("L1 members")
    ax.set_ylim(0, vals.max() * 1.12)
    ax.legend(
        handles=[
            Patch(
                color=get_color("primary"), label="no room: fold starts at the mature N-terminus"
            ),
            Patch(color=get_color("secondary"), label="at least one site before the fold"),
        ],
        loc="upper right",
        frameon=False,
    )
    save_figure(fig, "nterm_room_before_fold", output_dir=FIGURES)

    # Figure 2: verdict share by offset past the cleavage site, tagged designs only.
    rep = [p for p in positions if p["repredicted"] == "1"]
    by = collections.defaultdict(collections.Counter)
    for p in rep:
        by[int(p["offset_past_cleavage"])][p["verdict"]] += 1
    offs = np.arange(1, 21)
    totals = np.array([sum(by[o].values()) for o in offs], dtype=float)
    fig, ax = create_figure(figsize=(12, 6.5))
    bottom = np.zeros(len(offs))
    for key, label, color in VERDICTS:
        share = np.array([by[o][key] for o in offs]) / totals
        ax.bar(offs, share, bottom=bottom, color=get_color(color), label=label, width=0.8)
        bottom += share
    for o, t in zip(offs, totals, strict=True):
        ax.text(o, 1.01, f"{int(t)}", ha="center", va="bottom", fontsize=10, rotation=90)
    ax.set_xticks([1, 5, 10, 15, 20])
    ax.set_xlabel("cassette inserted this many residues past the cleavage site")
    ax.set_ylabel("share of tagged designs")
    ax.set_ylim(0, 1.12)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2, frameon=False)
    save_figure(fig, "nterm_reprediction_by_offset", output_dir=FIGURES)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
