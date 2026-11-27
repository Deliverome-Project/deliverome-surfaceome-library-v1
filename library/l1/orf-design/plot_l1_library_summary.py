"""L1 library in numbers, for the 2026-11-06 SAB deck.

Three figures, each counted from the tranche manifest and the tag-path table rather than
copied from a write-up:

* ``s4e9_l1_composition`` -- L1 members by candidate set, orthologs split by species.
* ``s4e9_l1_tag_placement`` -- where the tag goes, per member, by candidate set.
* ``s4e9_library_tranches`` -- L1 against L2, by candidate set.

The unit throughout is the library **member** (one ordered construct), not the gene. The
25 expression-rescue arms are L1 members but are not in ``l1_tag_paths.tsv`` (their tag
site is fixed by the arm design), so the placement figure counts 2,700, not 2,725.

    uv run --locked python plot_l1_library_summary.py
"""

from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt

from deliverome_analysis.cloning import library_design as ld
from deliverome_analysis.core.plotting import COLORS, save_figure, setup_plotting_style
from deliverome_analysis.core.s3 import processed

HERE = Path(__file__).parent
FIGURES = HERE / "figures"
TRANCHE_KEY = "sprint4/s4e9-orf-design-fg-panel/library-candidates/candidate_manifest_tranches.tsv"
TAG_PATHS = HERE / "l1_tag_paths.tsv"

SET_LABELS = {
    "library_reference": "Library reference\n(one proteoform per gene)",
    "isoform": "Alternative isoforms",
    "ortholog": "Orthologs",
    "expression_rescue": "GPCR expression-rescue arms",
    "giant_truncation": "Giant-protein truncations",
}
SET_COLORS = {
    "library_reference": COLORS["secondary"],
    "isoform": COLORS["quaternary"],
    "ortholog": COLORS["accent"],
    "expression_rescue": COLORS["primary"],
    "giant_truncation": COLORS["neutral"],
}
PLACEMENTS = [
    ("N after SP", "N-term, after signal peptide"),
    ("internal", "Internal loop or snorkel"),
    ("N at residue 1, multipass", "N-term at residue 1, multipass"),
    ("C-terminal", "C-terminal"),
    ("N at residue 1, single-pass (type III)", "N-term at residue 1, type III"),
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def load_tranches() -> list[dict[str, str]]:
    """S3 is the record; the committed copy is the offline fallback."""
    try:
        return read_tsv(processed.download(TRANCHE_KEY))
    except Exception as exc:  # any S3 failure falls back the same way
        print(f"S3 unavailable ({exc.__class__.__name__}); using committed copy")
        return read_tsv(HERE / "candidate_manifest_tranches.tsv")


def plot_composition(rows: list[dict[str, str]]) -> None:
    l1 = [r for r in rows if r["library"] == ld.LIBRARY_L1]
    order = ["library_reference", "ortholog", "isoform", "expression_rescue"]
    by_set = Counter(r["set"] for r in l1)
    species = Counter(r["species"] for r in l1 if r["set"] == "ortholog")

    fig, ax = plt.subplots(figsize=(9, 4.6))
    y = list(range(len(order)))[::-1]
    for yi, s in zip(y, order, strict=True):
        ax.barh(yi, by_set[s], color=SET_COLORS[s], height=0.65)
        label = f"{by_set[s]:,}"
        if s == "ortholog":
            # the bar is too short to split legibly, so the species split goes in the value
            label += f"  ({species['cyno']} cynomolgus, {species['mouse']} mouse)"
        ax.text(by_set[s] + 25, yi, label, va="center", fontsize=15)
    ax.set_yticks(y)
    ax.set_yticklabels([SET_LABELS[s] for s in order], fontsize=14)
    ax.set_xlabel(f"L1 members (total {len(l1):,})", fontsize=15)
    ax.set_xlim(0, max(by_set.values()) * 1.15)
    ax.tick_params(axis="x", labelsize=13)
    ax.grid(axis="y", visible=False)
    save_figure(fig, "s4e9_l1_composition", output_dir=FIGURES)
    plt.close(fig)
    print("L1 by set:", dict(by_set), "orthologs by species:", dict(species))


def plot_tag_placement(paths: list[dict[str, str]]) -> None:
    sets = ["library_reference", "isoform", "ortholog"]
    counts = Counter((r["placement"], r["candidate_set"]) for r in paths)
    totals = Counter(r["placement"] for r in paths)
    unknown = set(totals) - {p for p, _ in PLACEMENTS}
    assert not unknown, f"unmapped placements: {unknown}"

    fig, ax = plt.subplots(figsize=(9, 4.6))
    y = list(range(len(PLACEMENTS)))[::-1]
    for yi, (key, _) in zip(y, PLACEMENTS, strict=True):
        left = 0
        for s in sets:
            n = counts[(key, s)]
            ax.barh(yi, n, left=left, color=SET_COLORS[s], height=0.65)
            left += n
        share = totals[key] / len(paths)
        ax.text(left + 15, yi, f"{totals[key]:,} ({share:.0%})", va="center", fontsize=14)
    ax.set_yticks(y)
    ax.set_yticklabels([label for _, label in PLACEMENTS], fontsize=14)
    ax.set_xlabel(f"L1 members with a tag path (n = {len(paths):,})", fontsize=15)
    ax.set_xlim(0, max(totals.values()) * 1.3)
    ax.tick_params(axis="x", labelsize=13)
    ax.grid(axis="y", visible=False)
    handles = [plt.Rectangle((0, 0), 1, 1, color=SET_COLORS[s]) for s in sets]
    ax.legend(
        handles,
        ["Library reference", "Isoform", "Ortholog"],
        fontsize=12,
        loc="lower right",
        frameon=False,
    )
    save_figure(fig, "s4e9_l1_tag_placement", output_dir=FIGURES)
    plt.close(fig)
    print("placements:", dict(totals))


def plot_tranches(rows: list[dict[str, str]]) -> None:
    order = ["library_reference", "isoform", "ortholog", "expression_rescue", "giant_truncation"]
    counts = Counter((r["library"], r["set"]) for r in rows)
    libs = [ld.LIBRARY_L1, ld.LIBRARY_L2]

    fig, ax = plt.subplots(figsize=(9, 3.6))
    y = [1, 0]
    for yi, lib in zip(y, libs, strict=True):
        left = 0
        for s in order:
            n = counts[(lib, s)]
            if n:
                ax.barh(yi, n, left=left, color=SET_COLORS[s], height=0.6)
            left += n
        ax.text(left + 40, yi, f"{left:,}", va="center", fontsize=15)
    ax.set_yticks(y)
    ax.set_yticklabels(["L1 (designed)", "L2 (next tranche)"], fontsize=14)
    ax.set_xlabel("Library members", fontsize=15)
    ax.set_xlim(0, sum(counts[(ld.LIBRARY_L2, s)] for s in order) * 1.15)
    ax.tick_params(axis="x", labelsize=13)
    ax.grid(axis="y", visible=False)
    present = [s for s in order if any(counts[(lib, s)] for lib in libs)]
    handles = [plt.Rectangle((0, 0), 1, 1, color=SET_COLORS[s]) for s in present]
    labels = [SET_LABELS[s].split("\n")[0] for s in present]
    ax.legend(
        handles,
        labels,
        fontsize=11,
        ncol=3,
        loc="upper center",
        bbox_to_anchor=(0.45, -0.32),
        frameon=False,
    )
    save_figure(fig, "s4e9_library_tranches", output_dir=FIGURES)
    plt.close(fig)
    print("tranches:", {lib: sum(counts[(lib, s)] for s in order) for lib in libs})


def main() -> int:
    setup_plotting_style()
    rows = load_tranches()
    paths = read_tsv(TAG_PATHS)
    assert {r["library"] for r in paths} == {ld.LIBRARY_L1}
    plot_composition(rows)
    plot_tag_placement(paths)
    plot_tranches(rows)
    l1 = [r for r in rows if r["library"] == ld.LIBRARY_L1]
    print(
        "L1 status:",
        "barcode reserved",
        sum(r["barcode"] != "not reserved" for r in l1),
        "| tag cassette",
        dict(Counter(r["tag_cassette"] for r in l1)),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
