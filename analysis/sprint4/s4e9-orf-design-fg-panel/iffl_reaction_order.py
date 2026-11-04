"""Why a fully complementary site represses hard but does not dosage-compensate.

Du et al. Fig 1F: a single 21-nt fully complementary miR-L site cuts target expression
by 1-2 orders of magnitude yet leaves it proportional to dosage. Multimerising a weaker
17-nt site to four copies flattens it. That looks paradoxical -- the stronger site
regulates less well -- until you look at what the removal rate is proportional to.

In an IFFL the miRNA and the target are made from the same DNA, so **both** scale with
copy number N. Everything then turns on the reaction order of miRNA-mediated removal:

* **Second order**, rate = k x [miRNA] x [mRNA]. Doubling N doubles both, so removal
  goes up four-fold while production only doubles. Steady state aN/(d + k b N)
  saturates: dosage-invariant.
* **Zero order in mRNA**, rate = kcat x [miRNA]. The RISC pool is the bottleneck, so
  removal scales with N exactly as production does. Their difference also scales with
  N: a parallel line on a log-log plot. Strong repression, no compensation.

Full complementarity drives Ago2 slicing. Each engagement is a committed bind-cut-release
cycle, targets are never limiting, and the enzyme runs saturated -- zero order in mRNA.
Partial complementarity cannot slice; repression needs the RISC to find and occupy the
message, binding is the rate-limiting step, and the rate tracks the product of the two
concentrations -- second order. Hence a sweet spot: too much pairing saturates the
catalyst, too little does not repress, and multimerising recovers magnitude without
leaving the unsaturated regime. Consistent with their "most sensitive at 16-20 nt".

`--endogenous` applies the same reasoning to the miR-374-5p sites the TACT spacer
creates, and the answer is reassuring. miR-374 is non-slicing, so its removal is second
order too -- but its concentration is **fixed**, not scaling with N. Algebraically that
makes it indistinguishable from basal decay: it adds to d and nothing else. The
compensated plateau is A/(k*b) regardless of d, so **miR-374 cannot move the setpoint**.
It only raises the dosage at which compensation engages, which suppresses output below
that threshold and is invisible above it.

The caveat the algebra does not carry: miR-374 occupying a site also blocks miR-L from
it, which would reduce the effective k rather than raise d, and that *would* move the
plateau. It stays small because miR-L pairs over 17 nt and slices while miR-374 holds
8 nt, so miR-L wins the shared site on affinity and turnover, not on numbers.

This is a reading of the published behaviour, not a fit to their model; Du's own
treatment is in their Data S2.

    uv run --locked python iffl_reaction_order.py
    uv run --locked python iffl_reaction_order.py --endogenous
"""

from __future__ import annotations

import argparse
import math

A, D, B = 1.0, 1.0, 1.0  # production, basal decay, miRNA production, all per copy


def second_order(n: float, k: float) -> float:
    """Removal = k x miRNA x mRNA. Saturates -> dosage compensation."""
    return A * n / (D + k * B * n)


def zero_order(n: float, kcat: float) -> float:
    """Removal = kcat x miRNA, capped by the RISC pool. Stays proportional to dosage."""
    return max(0.0, (A - kcat * B) * n / D)


def endogenous_report() -> None:
    """What a fixed endogenous miRNA pool does to the dosage response."""
    k = 30.0
    print("A fixed endogenous pool adds to the effective decay constant and nothing else.\n")
    cols = (1.0, 2.0, 5.0, 10.0)
    print(f"{'copies':>7}" + "".join(f"{f'd_eff={d}':>13}" for d in cols) + "   output lost at 1x")
    for n in (1, 3, 10, 100):
        vals = [second_order_d(n, k, d) for d in cols]
        lost = f"{100 * (1 - vals[-1] / vals[0]):.0f}% at d_eff=10" if n == 1 else ""
        print(f"{n:>7}" + "".join(f"{v:>13.4f}" for v in vals) + f"   {lost}")
    print(f"\n  plateau as copies -> infinity: {A / (k * B):.4f} for every d_eff")
    print("  so the setpoint is untouched; only the threshold moves, and only")
    print("  below that threshold is any output lost. Even a 10-fold inflated")
    print("  decay costs under a quarter of the signal at one copy.")


def second_order_d(n: float, k: float, d_eff: float) -> float:
    return A * n / (d_eff + k * B * n)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--endogenous", action="store_true", help=endogenous_report.__doc__)
    ap.add_argument("--k", type=float, default=30.0, help="second-order rate constant")
    ap.add_argument("--kcat", type=float, default=0.97, help="per-miRNA turnover")
    args = ap.parse_args(argv)
    if args.endogenous:
        endogenous_report()
        return 0

    copies = [1, 3, 10, 30, 100, 300]
    print(
        f"{'copies':>7} {'unregulated':>12} {'4x partial':>12} {'1x perfect':>12}   fold vs unregulated"
    )
    for n in copies:
        u, p, s = A * n / D, second_order(n, args.k), zero_order(n, args.kcat)
        print(
            f"{n:>7} {u:>12.1f} {p:>12.3f} {s:>12.3f}   partial {u / p:>8.1f}x   perfect {u / s:>5.1f}x"
        )

    print()
    for label, f in (
        ("4x partial, second order", lambda n: second_order(n, args.k)),
        ("1x perfect, zero order  ", lambda n: zero_order(n, args.kcat)),
    ):
        lo, hi = 10, 300
        slope = (math.log10(f(hi)) - math.log10(f(lo))) / (math.log10(hi) - math.log10(lo))
        print(f"  log-log slope, {lo}->{hi} copies, {label}: {slope:+.2f}")
    print("  0 = dosage-invariant, 1 = fully dosage-proportional")
    print()
    print("  Note both give the same ~30x repression at one copy. The mechanisms are")
    print("  indistinguishable from a single-dosage measurement; only the slope separates")
    print("  them, which is why a dosage axis is the whole experiment.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
