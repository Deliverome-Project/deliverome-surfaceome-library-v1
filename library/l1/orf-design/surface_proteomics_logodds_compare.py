"""Does the surface-proteomics challenge set change if rule B uses the global log-odds?

The challenge set behind section 2's "the gates have been tested" paragraph is the list of
proteins that pass the cross-method survey's rule B and are in neither L1 nor L2. Rule B is

    n_surface >= 20  and  lod_membrane > 0.5

where, per protein, ``k_x`` is the number of publications of kind ``x`` that detect it and
``N_x`` is the number of publications of that kind in the snapshot, and

    lod_x = logit((k_surface + 0.5) / (N_surface + 1)) - logit((k_x + 0.5) / (N_x + 1)).

Three kinds of publication are counted (``surface_specific.py`` in
``library/l1/proteomics-gap/``): **surface** = admitted intact-cell surface-labeling
inputs, **membrane** = membrane-fractionation inputs, **global** = post-lysis whole-lysate
glycoproteomes. ``lod_global`` is already a column of the published tables; rule B reports it
and never filters on it. This script answers "what if it did":

========== ===========================================================
variant    filter (always with n_surface >= 20)
========== ===========================================================
membrane   lod_membrane > 0.5 (rule B as published)
global     lod_global > 0.5
both       lod_membrane > 0.5 and lod_global > 0.5
pooled     log-odds against membrane + global publications pooled > 0.5
========== ===========================================================

``pooled`` is the other reasonable reading of "global": one background rate from every
non-surface publication rather than a separate rate per kind.

Input, committed so the comparison regenerates offline: ``surface_proteomics_counts.tsv``,
written by ``surface_proteomics_snapshot.py`` from the S3 keys pinned in
``surface_proteomics_inputs.json`` (per-protein publication counts of each pinned snapshot, for
proteins with at least 20 surface publications, plus ``N`` per kind). The figure shows the
snapshot the record names as ``current``. Library membership is that table's ``library``
column, which joins on Ensembl gene via ``gene_resolver`` against
``candidate_manifest_tranches.tsv``.

Run from the repo root:

    uv run --locked python library/l1/orf-design/surface_proteomics_logodds_compare.py
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from deliverome_analysis.core.plotting import save_figure, setup_plotting_style

HERE = Path(__file__).resolve().parent
COUNTS = HERE / "surface_proteomics_counts.tsv"
OUT = HERE / "surface_proteomics_logodds_compare.tsv"
SUMMARY = HERE / "surface_proteomics_logodds_summary.tsv"
FIGURE = HERE / "figures" / "surface_proteomics_logodds_compare"

INPUTS = HERE / "surface_proteomics_inputs.json"
MIN_SURFACE = 20
LOD_MIN = 0.5
VARIANTS = ("membrane", "global", "both", "pooled")


def logit(p: float) -> float:
    return math.log(p / (1 - p))


def lod(k_s: int, n_s: int, k_x: int, n_x: int) -> float:
    return logit((k_s + 0.5) / (n_s + 1)) - logit((k_x + 0.5) / (n_x + 1))


def score(df: pd.DataFrame) -> pd.DataFrame:
    out = []
    for _, sub in df.groupby("snapshot", sort=False):
        n = {k: int(sub[f"N_{k}"].iloc[0]) for k in ("surface", "membrane", "global")}
        sub = sub.copy()
        sub["lod_membrane"] = [
            lod(s, n["surface"], m, n["membrane"])
            for s, m in zip(sub.n_surface_pubs, sub.n_membrane_pubs, strict=True)
        ]
        sub["lod_global"] = [
            lod(s, n["surface"], g, n["global"])
            for s, g in zip(sub.n_surface_pubs, sub.n_global_pubs, strict=True)
        ]
        sub["lod_pooled"] = [
            lod(s, n["surface"], m + g, n["membrane"] + n["global"])
            for s, m, g in zip(
                sub.n_surface_pubs, sub.n_membrane_pubs, sub.n_global_pubs, strict=True
            )
        ]
        out.append(sub)
    df = pd.concat(out)
    df["in_library"] = df.library.astype(str).str.match(r"l[12]")
    df["pass_membrane"] = df.lod_membrane > LOD_MIN
    df["pass_global"] = df.lod_global > LOD_MIN
    df["pass_both"] = df.pass_membrane & df.pass_global
    df["pass_pooled"] = df.lod_pooled > LOD_MIN
    for c in ("lod_membrane", "lod_global", "lod_pooled"):
        df[c] = df[c].round(3)
    return df


def figure(df: pd.DataFrame, current: str) -> None:
    import numpy as np

    setup_plotting_style()
    sub = df[df.snapshot.eq(current)].copy()
    in_v1 = sub.library.astype(str).str.contains(r"\bl1\b")
    pass_mem = sub.pass_membrane
    outside_pass = (~in_v1) & pass_mem
    outside_fail = (~in_v1) & (~pass_mem)

    fig = plt.figure(figsize=(10.6, 4.6))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.38, 1.0], hspace=0.42, wspace=0.28)
    ax = fig.add_subplot(gs[:, 0])
    ax_m = fig.add_subplot(gs[0, 1])
    ax_g = fig.add_subplot(gs[1, 1])

    c_v1 = "#0F766E"
    c_out_bg = "#C8C2BE"
    c_out_hi = "#B3243B"

    # Left panel: joint scatter with shaded membrane-passing region
    ax.axvspan(LOD_MIN, 4.0, color=c_v1, alpha=0.055, zorder=0)
    ax.scatter(
        sub.lod_membrane[outside_fail],
        sub.lod_global[outside_fail],
        s=9,
        c=c_out_bg,
        alpha=0.55,
        linewidths=0,
        label=f"Outside v1, membrane ≤ 0.5 (n = {int(outside_fail.sum()):,})",
        zorder=1,
    )
    ax.scatter(
        sub.lod_membrane[in_v1],
        sub.lod_global[in_v1],
        s=12,
        c=c_v1,
        alpha=0.58,
        linewidths=0,
        label=f"In v1 library (n = {int(in_v1.sum()):,})",
        zorder=2,
    )
    ax.scatter(
        sub.lod_membrane[outside_pass],
        sub.lod_global[outside_pass],
        s=26,
        c=c_out_hi,
        edgecolors="white",
        linewidths=0.5,
        alpha=0.92,
        label=f"Outside v1, membrane > 0.5 (n = {int(outside_pass.sum()):,})",
        zorder=3,
    )
    ax.axvline(LOD_MIN, color="#1f1718", lw=1.0, ls="--", zorder=4)
    ax.axhline(LOD_MIN, color="#6f5d5a", lw=0.9, ls=":", zorder=4)

    n_ul = int((outside_fail & sub.pass_global).sum())
    n_r_v1 = int((in_v1 & pass_mem).sum())
    ax.text(
        -2.55,
        2.95,
        f"{n_ul:,} non-v1 proteins pass lysate\nbackground but fail membrane",
        fontsize=8.0,
        color="#524644",
        va="top",
        ha="left",
        bbox=dict(boxstyle="round,pad=0.28", fc="white", ec="#e0d9d4", alpha=0.94),
    )
    ax.text(
        3.48,
        2.95,
        f"Membrane > 0.5\n{n_r_v1:,} in v1 · {int(outside_pass.sum())} outside",
        fontsize=8.0,
        color=c_v1,
        fontweight="bold",
        va="top",
        ha="right",
        bbox=dict(boxstyle="round,pad=0.28", fc="white", ec="#c7dedb", alpha=0.94),
    )
    ax.text(
        0.0,
        1.03,
        "A. Background comparison across recurring surface-proteomics proteins",
        transform=ax.transAxes,
        fontsize=9.5,
        fontweight="bold",
        color="#1f1718",
        ha="left",
        va="bottom",
    )
    ax.set_xlim(-2.7, 3.6)
    ax.set_ylim(-1.85, 3.35)
    ax.set_xlabel("Surface vs. membrane fractionation (log-odds)", fontsize=10.5)
    ax.set_ylabel("Surface vs. whole-lysate glycoproteome (log-odds)", fontsize=10.5)
    ax.tick_params(labelsize=9)
    ax.legend(frameon=True, facecolor="white", edgecolor="#e0d9d4", fontsize=7.8, loc="lower left")

    # Right panels: 1D distributions comparing separation of v1 vs outside-v1
    bins = np.linspace(-2.5, 3.5, 42)
    for panel, col, label_txt, ls in [
        (ax_m, "lod_membrane", "B. Membrane-fractionation background (clean separation)", "--"),
        (ax_g, "lod_global", "C. Whole-lysate glycoproteome background (heavy overlap)", ":"),
    ]:
        panel.hist(
            sub.loc[~in_v1, col],
            bins=bins,
            density=True,
            color="#B3243B",
            alpha=0.28,
            label="Outside v1",
            edgecolor="none",
        )
        panel.hist(
            sub.loc[in_v1, col],
            bins=bins,
            density=True,
            color=c_v1,
            alpha=0.48,
            label="In v1 library",
            edgecolor="none",
        )
        panel.axvline(LOD_MIN, color="#1f1718", lw=0.95, ls=ls)
        panel.text(
            0.0,
            1.05,
            label_txt,
            transform=panel.transAxes,
            fontsize=9.0,
            fontweight="bold",
            color="#1f1718",
            ha="left",
            va="bottom",
        )
        panel.set_xlim(-2.6, 3.5)
        panel.set_ylabel("Density", fontsize=9.5)
        panel.tick_params(labelsize=8.5)
    ax_m.legend(frameon=False, fontsize=8.0, loc="upper right")
    ax_g.set_xlabel("Log-odds (intact-cell surface vs. background)", fontsize=10.5)

    fig.subplots_adjust(left=0.075, right=0.98, top=0.90, bottom=0.13)
    save_figure(fig, FIGURE)


def main() -> None:
    current = json.loads(INPUTS.read_text())["current"]
    df = score(pd.read_csv(COUNTS, sep="\t", keep_default_na=False))

    summary = []
    for snap in df.snapshot.unique():
        sub = df[df.snapshot.eq(snap)]
        base = set(sub.uniprot_accession[sub.pass_membrane & ~sub.in_library])
        for v in VARIANTS:
            hit = sub[sub[f"pass_{v}"]]
            gap = set(hit.uniprot_accession[~hit.in_library])
            sym = dict(zip(sub.uniprot_accession, sub.gene_symbol, strict=True))
            summary.append(
                {
                    "snapshot": snap,
                    "variant": v,
                    "pass": len(hit),
                    "not_in_library": len(gap),
                    "gained_vs_membrane": ", ".join(sorted(sym[a] for a in gap - base)),
                    "lost_vs_membrane": ", ".join(sorted(sym[a] for a in base - gap)),
                }
            )
    summary = pd.DataFrame(summary)
    summary.to_csv(SUMMARY, sep="\t", index=False)
    cols = [
        "snapshot",
        "uniprot_accession",
        "gene_symbol",
        "catalog_tier",
        "n_surface_pubs",
        "n_membrane_pubs",
        "n_global_pubs",
        "lod_membrane",
        "lod_global",
        "lod_pooled",
        "pass_membrane",
        "pass_global",
        "pass_both",
        "pass_pooled",
        "in_library",
        "library",
    ]
    keep = df[[f"pass_{v}" for v in VARIANTS]].any(axis=1)
    df[keep][cols].sort_values(["snapshot", "gene_symbol"]).to_csv(OUT, sep="\t", index=False)
    pd.set_option("display.width", 250, "display.max_colwidth", 120)
    print(summary.to_string(index=False))
    figure(df, current)


if __name__ == "__main__":
    main()
