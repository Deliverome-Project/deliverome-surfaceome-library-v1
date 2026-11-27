"""Are the four loop insertion sites actually in loops? Checked against deposited structures.

The `surface_loop` lane gates on secondary structure taken from an AlphaFold model. For
both of these genes experimental structures exist, and they are the better evidence. This
draws the extracellular regions with each deposited structure's own DSSP assignment
(served by PDBe) so the question "is my site in a helix" is answered by looking.

Experimental structures used: CD20 8vgn (2.5 A), 8vgo (2.6 A), 6vja (3.3 A); GLUT4 7wsm
(3.25 A), 7wsn (3.31 A) -- all cryo-EM, so extramembrane helical extensions may in part be
stabilised by the solubilisation and bound Fabs. Treat a helix call as a strong risk
signal, not proof the position cannot take an insertion.

    uv run --locked python plot_loop_sites_vs_experiment.py
"""

from __future__ import annotations

import csv
import json
import urllib.request
from pathlib import Path

from deliverome_analysis.core.plotting import (
    COLORS,
    create_figure,
    get_color,
    save_figure,
    setup_plotting_style,
)

HERE = Path(__file__).parent
PDBE = "https://www.ebi.ac.uk/pdbe/api"

TARGETS = {
    "MS4A1 (CD20)": {
        "acc": "P11836",
        "pdbs": [("8vgn", "2.5 A"), ("8vgo", "2.6 A"), ("6vja", "3.3 A")],
        "extracellular": [(79, 84), (142, 188)],
        "sites": {79: "P79\n(p0032)", 161: "Y161\n(p0033)", 171: "N171\nretired"},
        "span": (70, 200),
    },
    "SLC2A4 (GLUT4)": {
        "acc": "P14672",
        "pdbs": [("7wsm", "3.25 A"), ("7wsn", "3.31 A")],
        "extracellular": [(46, 81), (193, 201), (309, 323)],
        "sites": {66: "P66\n(p0030)", 193: "E193\nretired", 320: "G320\n(p0031)"},
        "span": (40, 330),
    },
}


def fetch(url: str):
    return json.load(urllib.request.urlopen(url))


def helices(pdb: str, acc: str) -> list[tuple[int, int]]:
    best = fetch(f"{PDBE}/mappings/best_structures/{acc}")[acc]
    chain = next(r["chain_id"] for r in best if r["pdb_id"] == pdb)
    ss = fetch(f"{PDBE}/pdb/entry/secondary_structure/{pdb}")
    out = []
    for mol in ss[pdb]["molecules"]:
        for c in mol["chains"]:
            if c["chain_id"] != chain:
                continue
            for kind in ("helices", "strands"):
                for seg in c["secondary_structure"].get(kind, []):
                    out.append((seg["start"]["residue_number"], seg["end"]["residue_number"]))
    return sorted(out)


def main() -> int:
    setup_plotting_style()
    fig, axes = create_figure(2, 1, figsize=(12.4, 7.2))
    rows_out = []

    for ax, (title, spec) in zip(axes, TARGETS.items(), strict=True):
        lo, hi = spec["span"]
        tracks = [(pdb, res, helices(pdb, spec["acc"])) for pdb, res in spec["pdbs"]]

        for a, b in spec["extracellular"]:
            ax.axvspan(a, b, color=COLORS["light"], zorder=0)
            ax.annotate(
                f"extracellular {a}-{b}",
                ((a + b) / 2, len(tracks) - 0.18),
                ha="center",
                fontsize=7.8,
                color=COLORS["neutral"],
            )

        for i, (pdb, res, segs) in enumerate(tracks):
            y = len(tracks) - 1 - i
            ax.hlines(y, lo, hi, color=COLORS["line"], lw=2.2, zorder=1)
            for s, e in segs:
                if e < lo or s > hi:
                    continue
                ax.plot(
                    [max(s, lo), min(e, hi)],
                    [y, y],
                    color=COLORS["secondary"],
                    lw=9.5,
                    zorder=2,
                    solid_capstyle="butt",
                )
            ax.annotate(
                f"{pdb}  {res}",
                (lo - 2, y),
                ha="right",
                va="center",
                fontsize=8.2,
                color=COLORS["dark"],
            )

        for pos, label in spec["sites"].items():
            n_hel = sum(1 for _, _, segs in tracks if any(s <= pos <= e for s, e in segs))
            retired = "retired" in label
            colour = COLORS["neutral"] if retired else get_color("primary")
            ax.vlines(
                pos,
                -0.55,
                len(tracks) - 0.55,
                color=colour,
                lw=1.7,
                ls="--" if retired else "-",
                zorder=4,
            )
            ax.annotate(
                f"{label}\n{n_hel}/{len(tracks)} helix",
                (pos, -0.62),
                ha="center",
                va="top",
                fontsize=8.0,
                color=colour,
            )
            rows_out.append([title, pos, label.replace("\n", " "), n_hel, len(tracks)])

        ax.set_xlim(lo - 26, hi + 4)
        ax.set_ylim(-2.25, len(tracks) + 0.35)
        ax.set_yticks([])
        # The house style strips axes titles, so the panel is named inside it.
        ax.annotate(
            title,
            (0.004, 0.97),
            xycoords="axes fraction",
            ha="left",
            va="top",
            fontsize=10.5,
            fontweight="bold",
            color=COLORS["dark"],
        )
        ax.set_xlabel("residue")
        for side in ("left", "right", "top"):
            ax.spines[side].set_visible(False)

    axes[0].hlines(
        [],
        [],
        [],
        color=COLORS["secondary"],
        lw=9.5,
        label="helix or strand in the deposited structure",
    )
    axes[0].legend(loc="upper right", fontsize=8.0, frameon=False)
    fig.tight_layout()
    save_figure(fig, "s4e9_loop_sites_vs_experiment", output_dir=HERE / "figures")

    with (HERE / "loop_sites_vs_experiment.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["protein", "residue", "site", "structures_helix", "structures_total"])
        w.writerows(rows_out)
    print("wrote figures/s4e9_loop_sites_vs_experiment.png and loop_sites_vs_experiment.tsv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
