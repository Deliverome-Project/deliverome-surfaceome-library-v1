"""s1e2 — barcode length justification.

Why 20 nt? The collision math below shows 18 nt is already more than sufficient;
we standardize on 20 nt to match prevailing ORF-barcode practice (many ORF-BC
constructs use 20+ nt). This quantifies the trade-off behind the choice
(see README.md).

For each candidate length L we compute:

  * raw sequence space            S_raw = 4**L
  * effective space               S_eff = p(L) * 4**L
      where p(L) is the fraction of random L-mers passing our barcode design
      constraints (GC window, homopolymer-run cap, no start ATG, no AATAAA polyA
      signal, no BbsI/Esp3I recognition site on either strand) — estimated by
      Monte Carlo. These are the SAME constraints enforced by the production
      generator in `deliverome_analysis.codon_opt` (PR #10); we re-implement the
      check locally so this analysis stands alone on a main-based branch.
  * P(any collision)              1 - exp(-N*(N-1) / (2*S_eff))     (birthday)
  * expected collision fraction   ~ 1 - (1 - 1/S_eff)**(N-1)  ~= N/S_eff
      (the fraction of barcodes that share their code with >=1 other clone — the
       metric that actually matters, since colliding barcodes are just dropped)

for realistic library sizes N. Run:

    uv run python library/l1/vector-design/barcode_length_analysis.py
"""

import csv
import math
import random
from pathlib import Path

import matplotlib.pyplot as plt

from deliverome_analysis.core.plotting import CATEGORICAL_PALETTE, save_figure, setup_plotting_style

HERE = Path(__file__).resolve().parent

# Design constraints (mirror deliverome_analysis.codon_opt barcode generator).
GC_MIN, GC_MAX = 0.40, 0.60
MAX_HOMOPOLYMER = 3
# Recognition sites to avoid on both strands.
_SITES = ["GAAGAC", "GTCTTC", "CGTCTC", "GAGACG"]  # BbsI + Esp3I/BsmBI, both strands
_POLYA = "AATAAA"

LENGTHS = list(range(10, 25))
LIBRARY_SIZES = [10_000, 100_000, 1_000_000, 10_000_000]
CHOSEN_LENGTH = 20
ACCEPTABLE_COLLISION_FRACTION = 0.01  # 1%
MC_SAMPLES = 40_000
SEED = 20260722


def _passes_constraints(bc: str) -> bool:
    gc = (bc.count("G") + bc.count("C")) / len(bc)
    if not (GC_MIN <= gc <= GC_MAX):
        return False
    run = 1
    for i in range(1, len(bc)):
        run = run + 1 if bc[i] == bc[i - 1] else 1
        if run > MAX_HOMOPOLYMER:
            return False
    if bc.startswith("ATG"):
        return False
    if _POLYA in bc:
        return False
    return all(site not in bc for site in _SITES)


def _effective_fraction(length: int, rng: random.Random) -> float:
    hits = sum(
        _passes_constraints("".join(rng.choice("ACGT") for _ in range(length)))
        for _ in range(MC_SAMPLES)
    )
    return hits / MC_SAMPLES


def _p_any_collision(n: int, space: float) -> float:
    if space <= 0:
        return 1.0
    return 1.0 - math.exp(-n * (n - 1) / (2.0 * space))


def _expected_collision_fraction(n: int, space: float) -> float:
    if space <= 0:
        return 1.0
    return 1.0 - (1.0 - 1.0 / space) ** (n - 1)


def main() -> None:
    setup_plotting_style()
    rng = random.Random(SEED)

    rows = []
    eff_by_length = {}
    for length in LENGTHS:
        p = _effective_fraction(length, rng)
        s_raw = 4.0**length
        s_eff = p * s_raw
        eff_by_length[length] = s_eff
        for n in LIBRARY_SIZES:
            rows.append(
                {
                    "barcode_length": length,
                    "raw_space": f"{s_raw:.3e}",
                    "constraint_pass_fraction": f"{p:.3f}",
                    "effective_space": f"{s_eff:.3e}",
                    "library_size": n,
                    "p_any_collision": f"{_p_any_collision(n, s_eff):.3e}",
                    "expected_collision_fraction": f"{_expected_collision_fraction(n, s_eff):.3e}",
                }
            )

    out_tsv = HERE / "barcode_length_collision.tsv"
    with out_tsv.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {out_tsv}")

    # Report the chosen length.
    print(f"\n--- barcode length = {CHOSEN_LENGTH} nt ---")
    s_eff_18 = eff_by_length[CHOSEN_LENGTH]
    print(f"effective space ~ {s_eff_18:.3e}")
    for n in LIBRARY_SIZES:
        frac = _expected_collision_fraction(n, s_eff_18)
        print(f"  N={n:>10,}: expected collision fraction = {frac:.2e}")

    # Figure: expected collision fraction vs length, one line per library size.
    fig, ax = plt.subplots(figsize=(7, 5))
    for i, n in enumerate(LIBRARY_SIZES):
        y = [_expected_collision_fraction(n, eff_by_length[length]) for length in LENGTHS]
        color = CATEGORICAL_PALETTE[i % len(CATEGORICAL_PALETTE)]
        ax.plot(LENGTHS, y, marker="o", color=color, label=f"N = {n:,}")
    ax.axhline(ACCEPTABLE_COLLISION_FRACTION, ls="--", color="0.4", lw=1.5, zorder=0)
    ax.axvline(CHOSEN_LENGTH, ls=":", color="0.4", lw=1.5, zorder=0)
    ax.set_yscale("log")
    ax.set_xlabel("barcode length (nt)")
    ax.set_ylabel("expected fraction of barcodes lost to collision")
    ax.legend(title="library size", frameon=False, fontsize="small")
    save_figure(fig, "s1e2_barcode_length", output_dir=str(HERE / "figures"))
    print(f"wrote {HERE / 'figures' / 's1e2_barcode_length.svg'}")


if __name__ == "__main__":
    main()
