"""s1e2 — head-to-head: (Twist fragments + Twist clonal) vs (AbstractBio + stitching).

Two families of synthesis strategy, each with one free parameter, plotted as cost-vs-quality
frontiers so the trade can be read off directly:

  F-series  short end per scenario D, then **Twist clonal above a length threshold T**
            (sequence-perfect, $0.09/bp). Free parameter: T.
  E-series  every ORF >1.2 kb split into <=1.2 kb **AbstractBio pieces** (1:30,000) and
            stitched. Free parameter: the per-piece assembly fee, which is UNQUOTED.
  S-series  ** Synplexity clonal genes: FLAT $50/gene, hard cap 3 kb ** — flat-fee pricing
            covers 1,897 of 2,243 ORFs (84.6%); the 346 ORFs >3 kb must route elsewhere,
            and that tail is where all remaining cost and error now live.

⚠️ ASSUMPTIONS on Synplexity, flagged for confirmation:
   * "clonal gene" is taken to mean **sequence-verified, i.e. error rate 0**, as for Twist
     clonal. If it is only *mostly* clean, every S-series quality number moves.
   * $50 is assumed **per gene, flat, at any length up to the 3 kb cap** — no per-bp
     component and no minimum order / per-plate structure.
   * The 3 kb cap is treated as hard.
   An earlier version of this analysis modelled the flat fee with NO length cap; those
   numbers (H-series) were unrealistic and are superseded by the S-series below.

Everything else mirrors synthesis_cost_analysis.py exactly (same gene set, same error model,
ORF-only costs, no cassette).

Run:
    uv run python library/l1/vector-design/strategy_comparison.py
"""

import csv
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from deliverome_analysis.core.plotting import COLORS, save_figure, setup_plotting_style

HERE = Path(__file__).resolve().parent

# ---- vendor parameters (identical to synthesis_cost_analysis.py) ----
P_POOL, P_ABIO, P_FRAG, P_CLONAL = 0.00647, 0.01000, 0.07, 0.09
POOL_CAP, ABIO_CAP = 1800, 1200
E_POOL, E_ABIO, E_FRAG = 1 / 3000, 1 / 30000, 1 / 7500

# ---- error model ----
F_MIS, F_NON = 392 / 549, 23 / 549
INDEL_F = 0.30
MIS_SEV = 0.32  # keep the headline model comparable; §5.5 shows 0.427 shifts all curves
# together by ~1.5 points, so it does not change any ranking here.

# ---- Synplexity clonal genes ----
SYN_FLAT = 50.0  # $ per gene, flat, regardless of length
SYN_CAP = 3000  # hard length cap -- covers 1,897 / 2,243 ORFs (84.6%)


def p_alter(p):
    return p * ((1 - INDEL_F) * (F_MIS + F_NON) + INDEL_F)


def p_severe(p):
    return p * ((1 - INDEL_F) * (F_NON + F_MIS * MIS_SEV) + INDEL_F)


def load_lengths():
    rows = json.loads((HERE / "canonical_topology.json").read_text())
    return [r["len"] * 3 + 3 for r in rows]


def quality(lengths_and_errors):
    """(frac with no protein-altering error, frac with a severe error)."""
    ok = sev = 0.0
    for L, err in lengths_and_errors:
        ok += (1 - p_alter(err)) ** L
        sev += 1 - (1 - p_severe(err)) ** L
    n = len(lengths_and_errors)
    return ok / n, sev / n


def short_end(L):
    """Scenario-D routing below 1.8 kb: AbstractBio <=1.2k, pools 1.2-1.8k."""
    if L <= ABIO_CAP:
        return L * P_ABIO, E_ABIO
    return L * P_POOL, E_POOL


def f_series(cds, threshold):
    """Short end per D; Twist clonal at/above `threshold`, Twist fragments below it."""
    cost = 0.0
    qe = []
    for L in cds:
        if L <= POOL_CAP and L < threshold:
            c, e = short_end(L)
        elif L >= threshold:
            c, e = L * P_CLONAL, 0.0  # sequence-verified
        else:
            c, e = L * P_FRAG, E_FRAG
        cost += c
        qe.append((L, e))
    ok, sev = quality(qe)
    return dict(cost=cost, ok=ok, sev=sev)


def s_series(cds, syn_from, tail="clonal", fee=44.0, tail_rate=None):
    """Synplexity (flat, perfect) for `syn_from` <= L <= SYN_CAP; cheap per-bp below;
    `tail` handles L > SYN_CAP.

    syn_from -- length at/above which we pay Synplexity's $50 instead of a cheap per-bp
                vendor. Below it, scenario-D routing (AbstractBio <=1.2k, pools 1.2-1.8k,
                fragments in between if syn_from > 1.8k).
    tail     -- 'clonal' | 'fragments' | 'elegen' | 'stitch'
    """
    cost = 0.0
    qe = []
    for L in cds:
        if L > SYN_CAP:
            if tail == "clonal":
                c, e = L * P_CLONAL, 0.0
            elif tail == "fragments":
                c, e = L * P_FRAG, E_FRAG
            elif tail == "elegen":
                c, e = L * (tail_rate or 0.14), 1 / 70000
            elif tail == "stitch":
                npieces = math.ceil(L / ABIO_CAP)
                c, e = L * P_ABIO + npieces * fee, E_ABIO
            else:
                raise ValueError(tail)
        elif L >= syn_from:
            c, e = SYN_FLAT, 0.0  # sequence-verified clonal
        elif L <= ABIO_CAP:
            c, e = L * P_ABIO, E_ABIO
        elif L <= POOL_CAP:
            c, e = L * P_POOL, E_POOL
        else:
            c, e = L * P_FRAG, E_FRAG
        cost += c
        qe.append((L, e))
    ok, sev = quality(qe)
    return dict(cost=cost, ok=ok, sev=sev)


def e_series(cds, fee):
    """Every ORF >1.2 kb split into <=1.2 kb AbstractBio pieces, then stitched at `fee`."""
    cost = 0.0
    pieces = 0
    qe = []
    for L in cds:
        if L <= ABIO_CAP:
            cost += L * P_ABIO
        else:
            n = math.ceil(L / ABIO_CAP)
            pieces += n
            cost += L * P_ABIO + n * fee
        qe.append((L, E_ABIO))
    ok, sev = quality(qe)
    return dict(cost=cost, ok=ok, sev=sev, pieces=pieces)


def main() -> None:
    setup_plotting_style()
    cds = load_lengths()
    print(f"n = {len(cds)} ORFs, {sum(cds):,} bp total\n")

    # ---------- named scenarios ----------
    named = {}
    named["D  baseline (no clonal)"] = f_series(cds, threshold=10**9)
    named["F3 clonal >2.27 kb"] = f_series(cds, threshold=2273)
    named["F4 clonal >1.8 kb"] = f_series(cds, threshold=POOL_CAP + 1)
    named["E  stitch @ $44/piece"] = e_series(cds, 44)
    named["E  stitch @ $100/piece"] = e_series(cds, 100)
    named["S1 Synplexity all <=3 kb + clonal tail"] = s_series(cds, 1, tail="clonal")
    named["S2 Synplexity all <=3 kb + fragments tail"] = s_series(cds, 1, tail="fragments")
    named["S3 Synplexity 714 bp-3 kb + clonal tail"] = s_series(cds, 714, tail="clonal")
    named["S4 Synplexity 1.8-3 kb + clonal tail"] = s_series(cds, POOL_CAP + 1, tail="clonal")
    named["S7 AbstractBio <=1.2 kb + Synplexity 1.2-3 kb + clonal tail"] = s_series(
        cds, ABIO_CAP + 1, tail="clonal"
    )
    named["S5 Synplexity all <=3 kb + Elegen tail"] = s_series(cds, 1, tail="elegen")
    named["S6 Synplexity all <=3 kb + stitched tail @$44"] = s_series(
        cds, 1, tail="stitch", fee=44.0
    )

    rows = []
    print(f"{'scenario':46s} {'cost':>11s} {'usable':>8s} {'severe':>8s} {'$/usable':>10s}")
    for k, v in named.items():
        per = v["cost"] / (v["ok"] * len(cds))
        print(
            f"{k:46s} ${v['cost']:>10,.0f} {100 * v['ok']:>7.1f}% "
            f"{100 * v['sev']:>7.1f}% ${per:>9,.0f}"
        )
        rows.append(
            dict(
                scenario=k,
                cost=round(v["cost"]),
                pct_usable=round(100 * v["ok"], 1),
                pct_severe=round(100 * v["sev"], 1),
                dollars_per_usable_orf=round(per),
            )
        )

    tsv = HERE / "strategy_comparison.tsv"
    with tsv.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {tsv}")

    # ---------- what the 3 kb cap does ----------
    n_in = sum(1 for L in cds if L <= SYN_CAP)
    tail = [L for L in cds if L > SYN_CAP]
    print(f"\n=== the {SYN_CAP:,} bp Synplexity cap ===")
    print(
        f"  covered  : {n_in:,} ORFs ({100 * n_in / len(cds):.1f}%) at ${SYN_FLAT:.0f} each "
        f"= ${n_in * SYN_FLAT:,.0f}  -- sequence-perfect"
    )
    print(
        f"  tail     : {len(tail):,} ORFs ({100 * len(tail) / len(cds):.1f}%), "
        f"{sum(tail):,} bp, mean {sum(tail) / len(tail):,.0f} bp, max {max(tail):,} bp"
    )
    print("  tail cost by route:")
    for nm, rate in (
        ("Twist fragments $0.07", P_FRAG),
        ("Twist clonal $0.09", P_CLONAL),
        ("Elegen $0.14", 0.14),
        ("Elegen $0.20", 0.20),
    ):
        print(f"    {nm:24s} ${rate * sum(tail):>10,.0f}")

    print("\n  break-even: a $50 flat fee equals a per-bp vendor at...")
    for nm, rate in (
        ("AbstractBio $0.01", P_ABIO),
        ("pools $0.00647", P_POOL),
        ("Twist fragments $0.07", P_FRAG),
        ("Twist clonal $0.09", P_CLONAL),
    ):
        print(f"    {nm:24s} L = {SYN_FLAT / rate:>7,.0f} bp")
    print("  -> Synplexity beats fragments above 714 bp and clonal above 556 bp, but the")
    print("     cheap short-end vendors are still far cheaper in dollars below ~1.8 kb.")

    print(f"\n=== where to start paying ${SYN_FLAT:.0f} instead of a cheap per-bp vendor ===")
    print(f"{'syn_from':>9s} {'cost':>11s} {'usable':>8s} {'severe':>8s}  (clonal tail)")
    for t in (1, 400, 714, 1000, 1201, 1801, 2273, 3001):
        v = s_series(cds, t, tail="clonal")
        lab = "never" if t > SYN_CAP else f"{t:,}"
        print(f"{lab:>9s} ${v['cost']:>10,.0f} {100 * v['ok']:>7.1f}% {100 * v['sev']:>7.1f}%")

    # ---------- frontiers ----------
    thresholds = np.unique(np.concatenate([np.arange(600, 6000, 50), [10**9]]))
    f_curve = [f_series(cds, threshold=int(t)) for t in thresholds]
    syn_from = np.unique(np.concatenate([np.arange(1, SYN_CAP + 1, 25), [SYN_CAP + 1]]))
    s_curve = [s_series(cds, int(t), tail="clonal") for t in syn_from]
    fees = np.arange(0, 131, 5)
    e_curve = [e_series(cds, float(f)) for f in fees]

    # ---------- figure ----------
    fig, axes = plt.subplots(1, 3, figsize=(17, 5.6))
    C_F, C_E, C_H = COLORS["primary"], COLORS["secondary"], COLORS["quaternary"]

    for ax, key, ylab in (
        (axes[0], "ok", "% of ORFs with NO protein-altering error"),
        (axes[1], "sev", "% of ORFs carrying a severe error"),
    ):
        ax.plot(
            [v["cost"] / 1000 for v in f_curve],
            [100 * v[key] for v in f_curve],
            color=C_F,
            lw=2.6,
            label="Twist fragments + clonal $0.09/bp\n(sweeping the clonal threshold)",
        )
        ax.plot(
            [v["cost"] / 1000 for v in s_curve],
            [100 * v[key] for v in s_curve],
            color=C_H,
            lw=2.8,
            ls="--",
            label=(
                f"Synplexity ${SYN_FLAT:.0f}/gene flat to {SYN_CAP // 1000} kb + clonal tail\n"
                "(sweeping where the flat fee starts)"
            ),
        )
        ax.plot(
            [v["cost"] / 1000 for v in e_curve],
            [100 * v[key] for v in e_curve],
            color=C_E,
            lw=2.6,
            label="AbstractBio pieces + stitching\n(sweeping the per-piece fee)",
        )
        for label, mark in (
            ("D  baseline (no clonal)", "o"),
            ("F3 clonal >2.27 kb", "s"),
            ("E  stitch @ $44/piece", "^"),
            ("E  stitch @ $100/piece", "v"),
            ("S1 Synplexity all <=3 kb + clonal tail", "D"),
            ("S4 Synplexity 1.8-3 kb + clonal tail", "P"),
        ):
            v = named[label]
            col = C_E if label.startswith("E") else (C_H if label.startswith("S") else C_F)
            ax.plot(
                [v["cost"] / 1000],
                [100 * v[key]],
                mark,
                color=col,
                ms=9,
                mec="white",
                mew=1.2,
                zorder=5,
            )
            ax.annotate(
                label.split()[0],
                xy=(v["cost"] / 1000, 100 * v[key]),
                xytext=(7, 5),
                textcoords="offset points",
                fontsize=9,
                fontweight="bold",
                color=col,
            )
        ax.set_xlabel("cost ($k, ORF only)")
        ax.set_ylabel(ylab)
        if key == "ok":
            ax.legend(frameon=False, fontsize=7, loc="lower right")
            ax.annotate(
                "the stitching curve is horizontal by construction:\n"
                "the assembly fee moves cost, not error rate",
                xy=(0.03, 0.62),
                xycoords="axes fraction",
                fontsize=7.5,
                color=C_E,
                style="italic",
            )
        ax.annotate(
            "better ↑" if key == "ok" else "better ↓",
            xy=(0.02, 0.95 if key == "ok" else 0.05),
            xycoords="axes fraction",
            fontsize=8,
            color="0.4",
        )

    # panel 3: $ per usable ORF
    ax = axes[2]
    order = [
        "D  baseline (no clonal)",
        "F3 clonal >2.27 kb",
        "E  stitch @ $44/piece",
        "S1 Synplexity all <=3 kb + clonal tail",
        "S7 AbstractBio <=1.2 kb + Synplexity 1.2-3 kb + clonal tail",
        "S6 Synplexity all <=3 kb + stitched tail @$44",
    ]
    short = ["D", "F3", "E@44", "S1", "S7", "S6"]
    cols = [C_F, C_F, C_E, C_H, C_H, C_H]
    vals = [named[k]["cost"] / (named[k]["ok"] * len(cds)) for k in order]
    ax.bar(np.arange(len(vals)), vals, 0.62, color=cols)
    for i, v in enumerate(vals):
        ax.text(i, v + 5, f"${v:,.0f}", ha="center", fontsize=8.5)
    ax.set_xticks(np.arange(len(vals)))
    ax.set_xticklabels(short, fontsize=10)
    ax.set_ylabel("$ per USABLE ORF")
    ax.set_ylim(0, max(vals) * 1.18)

    save_figure(fig, "s1e2_strategy_comparison", output_dir=str(HERE / "figures"))
    print("wrote figures/s1e2_strategy_comparison.{svg,png}")


if __name__ == "__main__":
    main()
