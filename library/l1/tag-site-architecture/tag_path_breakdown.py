"""Where the tag cassette can go on each canonical library member, and what that rests on.

The placement ladder is ecto-N -> ecto-C -> internal loop -> snorkel, but the four rungs
are not equally safe, and two of them split on a distinction the rung name hides:

* **N after signal-peptide cleavage** is the bulk of the library. The tag enters an
  already-translocated mature protein, so the only question is where the peptidase cuts.
  Placed too early the tag leaves with the signal and the member reads as a non-expressor,
  which is indistinguishable in the data from a real biological negative.
* **N at residue 1** has no cleaved signal, so the tag sits in the N-tail that the signal
  anchor has to translocate. For the 7TM-dominated multipass members this is Tedman's own
  construct design, validated across hundreds of GPCRs. For single-pass type III members
  it is not: one TM carries the whole topogenic decision, and the charge balance across it
  sets the orientation, so a 51-aa cassette on the N-tail can invert the protein.
* **C-terminal** is only reached when there is no accessible N-terminus at all.
* **internal** needs an AlphaFold model and per-residue structural signals per gene.

Counts come from the committed DeepTMHMM predictions, not from UniProt subcellular text.
Reads the per-member table ``derive_tag_paths.py`` writes, so the figure reproduces from
this repo alone without a checkout of the surfaceome repo beside it.

    uv run --locked python library/l1/tag-site-architecture/tag_path_breakdown.py
"""

from __future__ import annotations

import csv
from pathlib import Path

from deliverome_analysis.plotting import get_categorical_colors, save_figure, setup_plotting_style

HERE = Path(__file__).parent
ROWS = HERE / "tag_path_assignments.tsv"
TABLE = HERE / "tag_path_breakdown.tsv"

# label -> (predicate, what the placement rests on)
GROUPS = (
    ("N-term, after signal peptide", lambda r: r["rung"] == "N after SP",
     "predicted cleavage site"),
    ("N-term at residue 1, multipass",
     lambda r: r["rung"] == "N at residue 1" and int(r["n_tm"]) >= 2, "Tedman GPCR screen"),
    ("N-term at residue 1, single-pass (type III)",
     lambda r: r["rung"] == "N at residue 1" and int(r["n_tm"]) <= 1, "topogenesis, untested"),
    ("C-terminal", lambda r: r["rung"] == "C-terminal", "topology only"),
    ("internal loop or snorkel", lambda r: r["rung"] == "internal", "AlphaFold model per gene"),
)


def main() -> int:
    rows = list(csv.DictReader(ROWS.open(), delimiter="\t"))
    counts = [(label, sum(1 for r in rows if pred(r)), rests) for label, pred, rests in GROUPS]
    total = len(rows)
    with TABLE.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["placement", "members", "pct_of_library", "rests_on"])
        for label, n, rests in counts:
            w.writerow([label, n, f"{100 * n / total:.1f}", rests])

    setup_plotting_style()
    import matplotlib.pyplot as plt

    colors = get_categorical_colors(len(counts))
    fig, ax = plt.subplots(figsize=(9.6, 4.6))
    ys = range(len(counts))[::-1]
    for y, (label, n, rests), c in zip(ys, counts, colors, strict=True):
        ax.barh(y, n, color=c, height=0.62)
        ax.text(n + 18, y, f"{n:,}  ({100 * n / total:.1f}%)", va="center", fontsize=9)
        ax.text(18, y - 0.34, rests, va="top", fontsize=7.5, color="0.4", style="italic")
    ax.set_yticks(list(ys))
    ax.set_yticklabels([label for label, _n, _r in counts], fontsize=9)
    ax.set_xlabel(f"canonical library members (n={total:,})")
    ax.set_xlim(0, max(n for _l, n, _r in counts) * 1.26)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0)
    fig.tight_layout()
    save_figure(fig, "tag_path_breakdown", output_dir=HERE / "figures")

    for label, n, rests in counts:
        print(f"  {label:<44} {n:>5,}  {100 * n / total:>5.1f}%   {rests}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
