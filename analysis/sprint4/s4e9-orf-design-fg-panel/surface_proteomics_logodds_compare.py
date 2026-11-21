"""Does the surface-proteomics challenge set change if rule B uses the global log-odds?

The challenge set behind section 2's "the gates have been tested" paragraph is the list of
proteins that pass the cross-method survey's rule B and are in neither L1 nor L2. Rule B is

    n_surface >= 20  and  lod_membrane > 0.5

where, per protein, ``k_x`` is the number of publications of kind ``x`` that detect it and
``N_x`` is the number of publications of that kind in the snapshot, and

    lod_x = logit((k_surface + 0.5) / (N_surface + 1)) - logit((k_x + 0.5) / (N_x + 1)).

Three kinds of publication are counted (``surface_specific.py`` in
``scratch/becca/s3e8-library-gap/``): **surface** = admitted intact-cell surface-labeling
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
the per-protein publication counts of each snapshot reduced to proteins with at least 20
surface publications (no protein below that can pass any variant), plus ``N`` per kind.
``--refresh <snapshot>=<surface_counts_all_proteins.tsv>`` rebuilds it from the full tables
``surface_specific.py`` writes. Library membership is that table's ``library`` column, which
joins on Ensembl gene via ``gene_resolver`` against ``candidate_manifest_tranches.tsv``.

Run from the repo root:

    uv run --locked python analysis/sprint4/s4e9-orf-design-fg-panel/surface_proteomics_logodds_compare.py
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from deliverome_analysis.plotting import get_color, save_figure, setup_plotting_style

HERE = Path(__file__).resolve().parent
COUNTS = HERE / "surface_proteomics_counts.tsv"
OUT = HERE / "surface_proteomics_logodds_compare.tsv"
SUMMARY = HERE / "surface_proteomics_logodds_summary.tsv"
FIGURE = HERE / "figures" / "surface_proteomics_logodds_compare"

#: Publications per kind in each snapshot, as rule_variants.py counts them.
N_PUBS = {
    "sweep-v16": {"surface": 71, "membrane": 77, "global": 18},
    "sweep-v20": {"surface": 73, "membrane": 77, "global": 18},
    "abundance-v4": {"surface": 84, "membrane": 77, "global": 18},
}
SNAPSHOTS = list(N_PUBS)
MIN_SURFACE = 20
LOD_MIN = 0.5
VARIANTS = ("membrane", "global", "both", "pooled")


def logit(p: float) -> float:
    return math.log(p / (1 - p))


def lod(k_s: int, n_s: int, k_x: int, n_x: int) -> float:
    return logit((k_s + 0.5) / (n_s + 1)) - logit((k_x + 0.5) / (n_x + 1))


def refresh(pairs: list[str]) -> None:
    frames = []
    for pair in pairs:
        snap, path = pair.split("=", 1)
        df = pd.read_csv(path, sep="\t", keep_default_na=False)
        df = df[df.n_surface_pubs >= MIN_SURFACE]
        df.insert(0, "snapshot", snap)
        frames.append(
            df[
                [
                    "snapshot",
                    "uniprot_accession",
                    "gene_symbol",
                    "catalog_tier",
                    "n_surface_pubs",
                    "n_membrane_pubs",
                    "n_global_pubs",
                    "library",
                ]
            ]
        )
    pd.concat(frames).to_csv(COUNTS, sep="\t", index=False)
    print(f"wrote {COUNTS.name}")


def score(df: pd.DataFrame) -> pd.DataFrame:
    out = []
    for snap, sub in df.groupby("snapshot", sort=False):
        n = N_PUBS[snap]
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


def figure(df: pd.DataFrame) -> None:
    setup_plotting_style()
    sub = df[df.snapshot.eq("abundance-v4")]
    fig, ax = plt.subplots(figsize=(8, 7))
    groups = [
        ("in L1 or L2", sub.in_library, get_color("neutral"), 6, 0.35),
        ("not in L1 or L2", ~sub.in_library, get_color("primary"), 10, 0.6),
    ]
    for label, mask, color, size, alpha in groups:
        ax.scatter(
            sub.lod_membrane[mask],
            sub.lod_global[mask],
            s=size,
            c=color,
            alpha=alpha,
            label=label,
            linewidths=0,
        )
    ax.axvline(LOD_MIN, color="0.3", lw=0.8, ls="--")
    ax.axhline(LOD_MIN, color="0.3", lw=0.8, ls=":")
    ax.set_xlabel("log-odds, surface vs membrane fractionation")
    ax.set_ylabel("log-odds, surface vs whole-lysate glycoproteome")
    handles, labels = ax.get_legend_handles_labels()
    handles += [
        plt.Line2D([], [], color="0.3", lw=0.8, ls="--"),
        plt.Line2D([], [], color="0.3", lw=0.8, ls=":"),
    ]
    labels += ["membrane threshold 0.5 (rule B)", "global threshold 0.5"]
    ax.legend(handles, labels, frameon=False, fontsize=7, loc="lower right")
    save_figure(fig, FIGURE)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--refresh", nargs="+", metavar="SNAPSHOT=PATH")
    args = ap.parse_args()
    if args.refresh:
        refresh(args.refresh)
    df = score(pd.read_csv(COUNTS, sep="\t", keep_default_na=False))

    summary = []
    for snap in SNAPSHOTS:
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
    figure(df)


if __name__ == "__main__":
    main()
