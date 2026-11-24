"""The bottom line of barcode swapping in L1, with and without codon diversification.

Becca, 2026-11-22, on the two-panel swap figure: "we need to know the bottom line for the lib with
or without codon diversification." One stacked bar per scenario, all on the order of record
(``audit_domestication_barcodes.RECODED_KEY``, T7 promoter included) and all from the same pooled
model (``l1_swap_pool_model.py``, Walton et al. 2024 fit):

* no codon diversification: ``swap_nothing_diversified_scenario.json``;
* linkers diversified, the order of record: ``l1_swap_pool_model_summary.json``;
* linkers and tags diversified (not chosen): ``swap_tags_recoded_scenario.json``.

Each bar is the share of integrants whose barcode names the wrong protein, split by the source of
the shared sequence, using the model's ``library_fraction_wrong_protein`` terms:

* shared tail before the barcode, T7 included = ``background_mislabel`` (a switch in the constant
  tail leaves the other construct's intact protein under this barcode);
* tag cassette fusions = ``background_chimera`` (a switch inside a shared tag fuses two ORFs);
* specific partners = ``specific_mislabel + specific_chimera`` (siblings, isoforms, orthologs).

The script checks that each background term sits mostly in the region its label names
(``wrong_protein_by_block_region``), so the labels cannot drift from the numbers. Every number
drawn is read from the JSONs; none is typed in. Leader swaps (same mature protein behind another
leader) are not counted, as in the model.

    uv run --locked python plot_swap_bottom_line.py
"""

from __future__ import annotations

import json
from pathlib import Path

from matplotlib.patches import Patch

from deliverome_analysis.plotting import (
    COLORS,
    SEQUENTIAL_PALETTES,
    create_figure,
    save_figure,
    setup_plotting_style,
)

HERE = Path(__file__).resolve().parent
OUT_DIR = HERE / "reviewer-package" / "figures"
SCENARIOS = (  # top to bottom
    ("No codon diversification", "swap_nothing_diversified_scenario.json"),
    ("Linkers diversified\n(order of record)", "l1_swap_pool_model_summary.json"),
    (
        "Linkers and tags diversified\n(not chosen: tags not shown to display)",
        "swap_tags_recoded_scenario.json",
    ),
)
SEGMENTS = (  # (label, colour, terms summed)
    (
        "Shared tail before the barcode, T7 included",
        SEQUENTIAL_PALETTES["teal"][3],
        ("background_mislabel",),
    ),
    ("Tag cassette fusions", COLORS["primary"], ("background_chimera",)),
    (
        "Specific partners (siblings, isoforms, orthologs)",
        COLORS["accent"],
        ("specific_mislabel", "specific_chimera"),
    ),
)
#: the block region each background segment is named for, and the least share of it that must lie
#: there for the label to hold
REGION_CHECK = {
    "background_mislabel": (("background|mislabel|3' tail",), 0.95),
    "background_chimera": (("background|chimera|tag cassette",), 0.85),
}
INSIDE_MIN = 0.6  # percentage points: segments narrower than this are labelled beside the bar


def load(fname: str) -> dict[str, float]:
    d = json.loads((HERE / fname).read_text())
    frac, regions = d["library_fraction_wrong_protein"], d["wrong_protein_by_block_region"]
    for term, (keys, floor) in REGION_CHECK.items():
        share = sum(regions[k] for k in keys) / frac[term]
        assert share >= floor, f"{fname}: only {share:.0%} of {term} lies in {keys}"
    vals = {label: 100 * sum(frac[t] for t in terms) for label, _c, terms in SEGMENTS}
    total = 100 * frac["total"]
    assert abs(sum(vals.values()) - total) < 1e-9, f"{fname}: segments do not sum to the total"
    return {**vals, "total": total}


def main() -> int:
    setup_plotting_style()
    data = [(name, load(f)) for name, f in SCENARIOS]
    fig, ax = create_figure(figsize=(9, 3.6))
    ys = list(range(len(data)))[::-1]
    xmax = max(v["total"] for _n, v in data)
    for y, (_name, v) in zip(ys, data, strict=True):
        left, beside = 0.0, []
        for label, colour, _terms in SEGMENTS:
            w = v[label]
            ax.barh(y, w, left=left, height=0.55, color=colour, edgecolor="white", linewidth=0.8)
            if w >= INSIDE_MIN:
                ax.text(
                    left + w / 2,
                    y,
                    f"{w:.2f}%",
                    ha="center",
                    va="center",
                    color="white",
                    fontsize=10,
                )
            else:
                beside.append(f"{label.split(' (')[0].split(',')[0].lower()} {w:.2f}%")
            left += w
        ax.text(
            left + 0.08,
            y + 0.06,
            f"{v['total']:.2f}%",
            ha="left",
            va="bottom",
            fontsize=12,
            fontweight="bold",
            color=COLORS["dark"],
        )
        if beside:
            ax.text(
                left + 0.08,
                y - 0.04,
                "incl. " + "; ".join(beside),
                ha="left",
                va="top",
                fontsize=8.5,
                color=COLORS["neutral"],
            )
    ax.set_yticks(ys, [n for n, _v in data])
    ax.set_xlim(0, xmax * 1.45)
    ax.set_ylim(-0.6, len(data) - 0.4)
    ax.set_xlabel("Integrants whose barcode names the wrong protein (%)")
    ax.grid(axis="y", visible=False)
    ax.spines[["top", "right"]].set_visible(False)
    handles = [Patch(facecolor=c, label=label) for label, c, _t in SEGMENTS]
    ax.legend(
        handles=handles,
        loc="lower left",
        bbox_to_anchor=(0.0, 1.0),
        ncol=3,
        frameon=False,
        fontsize=9,
        handlelength=1.2,
        columnspacing=1.2,
    )
    ax.text(
        0.0,
        -0.3,
        "Upper bounds: every virion is assumed to carry two different constructs. "
        "All bars: order of record 2026-11-22-rebuild, same model and tail.",
        ha="left",
        va="top",
        fontsize=8.5,
        color=COLORS["neutral"],
        transform=ax.transAxes,
    )
    save_figure(fig, "l1_swap_bottom_line", output_dir=OUT_DIR)
    for name, v in data:
        parts = ", ".join(f"{label.split(' (')[0]} {v[label]:.2f}%" for label, _c, _t in SEGMENTS)
        print(f"{name.replace(chr(10), ' ')}: total {v['total']:.2f}% ({parts})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
