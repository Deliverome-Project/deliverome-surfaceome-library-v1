"""s1e2 — synthesis cost + error burden for the canonical surfaceome ORF library.

Gene set: the **2,243 canonical** genes from accessible-surfaceome PR #130
(`passes_canonical` over the 5,130-gene deep-dive cohort). Lengths are real
`topology_public.protein_length` values pulled from the read-only public D1
mirror; all 2,243 matched (0 unmatched). CDS = 3*aa + 3 (stop).

Vendor parameters (as supplied):
  * Twist gene/oligo pools   $0.00647/bp, <= 1.8 kb, error ~1:3000
  * Twist gene fragments     $0.07/bp,               error ~1:7500
  * Abstract Bio             $0.01/bp,    <= ~1.2 kb, error 1:20k - 1:40k

IMPORTANT: costs here are **ORF only** — no tags, barcode, linkers or adapters.
The cassette scenario adds those to show the threshold-crossing penalty.

Run:
    uv run python analysis/sprint1/s1e2-fragmid-vector-design/synthesis_cost_analysis.py
"""

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from deliverome_analysis.plotting import COLORS, save_figure, setup_plotting_style

HERE = Path(__file__).resolve().parent

# ---- vendor parameters ----
P_POOL, P_ABIO, P_FRAG = 0.00647, 0.01000, 0.07
P_CLONAL = 0.09  # Twist clonal genes -- sequence-perfect, so error rate is 0
POOL_CAP, ABIO_CAP = 1800, 1200
E_POOL, E_ABIO, E_FRAG = 1 / 3000, 1 / 30000, 1 / 7500  # AbstractBio = 1:20k-1:40k midpoint

# ---- error model ----
F_MIS, F_NON = 392 / 549, 23 / 549  # exact single-nt spectrum, 61 sense codons
INDEL_F = 0.30  # UNVERIFIED for Twist -- biggest lever on the severe column
MIS_SEV = 0.32  # AlphaMissense: 32% of all possible human missense likely pathogenic

CASSETTE_BP = 275  # GS linkers + SpyTag003 + ALFA + barcode + primer sites

TOPOLOGY_JSON = HERE / "canonical_topology.json"


def load_lengths():
    """CDS lengths for the 2,243 canonical genes (falls back to the cached copy)."""
    rows = json.loads(TOPOLOGY_JSON.read_text())
    for r in rows:
        r["cds"] = r["len"] * 3 + 3
    return rows


def p_alter(p, indel_f=INDEL_F):
    return p * ((1 - indel_f) * (F_MIS + F_NON) + indel_f)


def p_severe(p, indel_f=INDEL_F):
    return p * ((1 - indel_f) * (F_NON + F_MIS * MIS_SEV) + indel_f)


def route(cds, abio_cap=ABIO_CAP, pool_cap=POOL_CAP, use_abio=True, use_pool=True):
    """Return (price, error_rate) for one ORF under the hybrid routing."""
    if use_abio and cds <= abio_cap:
        return P_ABIO, E_ABIO
    if use_pool and cds <= pool_cap:
        return P_POOL, E_POOL
    return P_FRAG, E_FRAG


SCENARIOS = {
    "A  all fragments": dict(use_abio=False, use_pool=False),
    "B  Twist pools ≤1.8k\n   + fragments": dict(use_abio=False, use_pool=True),
    "C  AbstractBio ≤1.2k\n   + fragments": dict(use_abio=True, use_pool=False),
    "D  AbstractBio ≤1.2k\n   + pools 1.2-1.8k\n   + fragments": dict(use_abio=True, use_pool=True),
}


def evaluate(rows, cassette=0, copies=1, **kw):
    cost = 0.0
    ok = 0.0
    sev = 0.0
    for r in rows:
        L = r["cds"] + cassette
        price, err = route(L, **kw)
        cost += L * price * copies
        ok += (1 - p_alter(err)) ** L
        sev += 1 - (1 - p_severe(err)) ** L
    n = len(rows)
    return dict(cost=cost, frac_ok=ok / n, frac_severe=sev / n, n=n)


def main() -> None:
    setup_plotting_style()
    rows = load_lengths()
    cds = np.array([r["cds"] for r in rows])
    print(
        f"n = {len(rows)} ORFs   total = {cds.sum():,} bp   mean = {cds.mean():.0f}  median = {np.median(cds):.0f}"
    )

    # ---------- table ----------
    out = []
    for name, kw in SCENARIOS.items():
        base = evaluate(rows, **kw)
        cass = evaluate(rows, cassette=CASSETTE_BP, **kw)
        x2 = evaluate(rows, copies=2, **kw)
        out.append(
            dict(
                scenario=name.replace("\n", " ").strip(),
                cost_orf_only=round(base["cost"]),
                cost_with_cassette=round(cass["cost"]),
                cost_2x_variants=round(x2["cost"]),
                pct_no_protein_altering=round(100 * base["frac_ok"], 1),
                usable_orfs=round(base["frac_ok"] * base["n"]),
                pct_with_severe_error=round(100 * base["frac_severe"], 1),
            )
        )
        print(
            f"  {name.replace(chr(10), ' '):<52} ${base['cost']:>10,.0f} | "
            f"+cassette ${cass['cost']:>10,.0f} | 2x ${x2['cost']:>10,.0f} | "
            f"ok {100 * base['frac_ok']:.1f}%"
        )

    tsv = HERE / "synthesis_cost_by_scenario.tsv"
    with tsv.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(out)
    print(f"wrote {tsv}")

    # ---------- figure: 3 panels ----------
    fig, axes = plt.subplots(1, 3, figsize=(16.5, 5.6))
    labels = [n.strip().split()[0] for n in SCENARIOS]  # A / B / C / D
    xs = np.arange(len(labels))

    # panel 1: cost, three bar groups
    ax = axes[0]
    base_costs = [o["cost_orf_only"] / 1000 for o in out]
    cass_costs = [o["cost_with_cassette"] / 1000 for o in out]
    x2_costs = [o["cost_2x_variants"] / 1000 for o in out]
    w = 0.27
    ax.bar(xs - w, base_costs, w, label="ORF only", color=COLORS["primary"])
    ax.bar(xs, cass_costs, w, label=f"+ cassette (+{CASSETTE_BP} bp)", color=COLORS["secondary"])
    ax.bar(xs + w, x2_costs, w, label="2 variants per ORF", color=COLORS["quaternary"])
    for i, v in enumerate(base_costs):
        ax.text(i - w, v + 8, f"${v:.0f}k", ha="center", fontsize=8)
    for i, v in enumerate(x2_costs):
        ax.text(i + w, v + 8, f"${v:.0f}k", ha="center", fontsize=8)
    ax.set_xticks(xs)
    ax.set_xticklabels(labels, fontsize=11)
    ax.set_ylabel("synthesis cost ($ thousands)")
    ax.legend(frameon=False, fontsize=8)
    ax.text(
        0.0,
        -0.20,
        "A all fragments   ·   B pools ≤1.8k + frag   ·   C AbstractBio ≤1.2k + frag\n"
        "D AbstractBio ≤1.2k + pools 1.2-1.8k + frag",
        transform=ax.transAxes,
        fontsize=8,
        color="0.35",
        va="top",
    )

    # panel 2: usable ORFs
    ax = axes[1]
    okp = [o["pct_no_protein_altering"] for o in out]
    sevp = [o["pct_with_severe_error"] for o in out]
    ax.bar(xs - 0.19, okp, 0.38, label="no protein-altering error", color=COLORS["primary"])
    ax.bar(xs + 0.19, sevp, 0.38, label="has highly-detrimental error", color=COLORS["accent"])
    for i, v in enumerate(okp):
        ax.text(i - 0.19, v + 1.2, f"{v:.0f}%", ha="center", fontsize=8)
    for i, v in enumerate(sevp):
        ax.text(i + 0.19, v + 1.2, f"{v:.0f}%", ha="center", fontsize=8)
    ax.set_xticks(xs)
    ax.set_xticklabels(labels, fontsize=11)
    ax.set_ylabel("% of the 2,243 ORFs")
    ax.set_ylim(0, 100)
    ax.legend(frameon=False, fontsize=8)

    # panel 3: cost concentration by length
    ax = axes[2]
    order = np.sort(cds)[::-1]
    cum_bp = np.cumsum(order) / order.sum() * 100
    rank_pct = np.arange(1, len(order) + 1) / len(order) * 100
    ax.plot(rank_pct, cum_bp, color=COLORS["primary"], lw=2.6)
    ax.plot([0, 100], [0, 100], ls=":", color="0.6", lw=1.4)
    for k, lab in ((50, "top 50"), (100, "top 100"), (500, "top 500")):
        xp = k / len(order) * 100
        yp = cum_bp[k - 1]
        ax.plot([xp], [yp], "o", color=COLORS["secondary"], ms=8)
        ax.annotate(
            f"{lab}\n{yp:.0f}% of bp",
            xy=(xp, yp),
            xytext=(9, -14),
            textcoords="offset points",
            fontsize=8,
            color=COLORS["secondary"],
        )
    ax.set_xlabel("ORFs, ranked longest → shortest (%)")
    ax.set_ylabel("cumulative share of total bp (%)")
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)

    save_figure(fig, "s1e2_synthesis_cost", output_dir=str(HERE / "figures"))
    print("wrote figures/s1e2_synthesis_cost.{svg,png}")


if __name__ == "__main__":
    main()
