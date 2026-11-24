"""s1e2 — barcode design: miRNA-seed exclusion, edit-distance capacity, and the UMI budget.

Answers three questions that `barcode_length_analysis.py` (length/collision only) does not:

1. **Are we excluding miRNA seed sites?**  We were not. Our barcode is *transcribed* — it
   sits between the ORF stop and IRES2, an intercistronic region that behaves like a 3'UTR
   for cistron 1 — so a barcode that happens to present a miRNA 7mer-m8 site can get its own
   transcript repressed. The failure mode is a **barcode-specific expression penalty that
   mimics a weak surface expressor**, which is very hard to separate from real biology.
   Seed sites come from TargetScanHuman 8.0 (`mirna_seed_sites.tsv`, tier-selectable).

2. **How many barcodes fit at a given length with a real error-correcting distance?**  The
   existing generator imposes **no distance constraint at all** — barcodes are just random
   passing sequences, so there is no error-correction guarantee. Here we build sets with a
   guaranteed **minimum Levenshtein distance**, which is what lets a sequencing read with an
   indel or substitution still be assigned unambiguously.

3. **Can we shorten the barcode to free bases for a UMI?**  Swept below.

4. **How wide should an added molecular UMI be?**  Sized on the number of items sharing one
   barcode, reported as the *undercount* -- the fraction of distinct items lost to UMI
   collision. Feeds the S1E6 ELN page.

METHOD NOTE — how the distance guarantee is enforced cheaply. Pairwise Levenshtein over
20,000 codewords is ~2e8 DP computations, far too slow. Instead we use deletion
neighbourhoods: for equal-length strings, if `Levenshtein(a, b) <= k` then `a` and `b`
share at least one common subsequence with exactly `k` characters deleted. So we keep the
set of all (d-1)-deletion subsequences of accepted codewords and reject any candidate that
collides. Because the condition is necessary but not sufficient, this **over-rejects**: the
accepted set is guaranteed to have min distance >= d, and the achieved counts below are a
**lower bound** on true capacity. Conservative in the safe direction for a design decision.

Run:
    uv run python library/l1/vector-design/barcode_design.py
"""

import csv
import itertools
import math
import random
from pathlib import Path

import matplotlib.pyplot as plt

from deliverome_analysis.plotting import (
    CATEGORICAL_PALETTE,
    COLORS,
    save_figure,
    setup_plotting_style,
)

HERE = Path(__file__).resolve().parent
SEED_TSV = HERE / "mirna_seed_sites.tsv"

# ---- existing design constraints (mirrors barcode_length_analysis.py) ----
GC_MIN, GC_MAX = 0.40, 0.60
MAX_HOMOPOLYMER = 3
_SITES = ["GAAGAC", "GTCTTC", "CGTCTC", "GAGACG"]  # BbsI + Esp3I/BsmBI, both strands
_POLYA = "AATAAA"

# ---- the design point under test ----
CURRENT_LENGTH = 20  # what the README currently specifies
TARGET_N = 20_000  # library size the user asked about
LENGTHS = [14, 15, 16, 17, 18, 19, 20]
DISTANCES = [2, 3, 4]  # d=3 corrects 1 edit; d=5 corrects 2 but is not tractable here
SEED = 20260811

# Which conservation tiers to exclude. 'broadly_conserved' + 'conserved' is the
# recommended default -- see the README for why 'all' over-filters.
TIER_SETS = {
    "none": set(),
    "broadly_conserved": {"broadly_conserved"},
    "conserved+": {"broadly_conserved", "conserved"},
    "all_human": {"broadly_conserved", "conserved", "poorly_conserved", "unassigned"},
}
DEFAULT_TIERS = "conserved+"


def load_seed_sites(tiers: set[str]) -> set[str]:
    """DNA 7-mer target sites to exclude, for the requested conservation tiers."""
    if not tiers:
        return set()
    if not SEED_TSV.exists():
        raise SystemExit(
            f"{SEED_TSV} not found. Regenerate from TargetScanHuman 8.0:\n"
            "  curl -sLO https://www.targetscan.org/vert_80/vert_80_data_download/"
            "miR_Family_Info.txt.zip && unzip -o miR_Family_Info.txt.zip\n"
            "then keep Species ID 9606 and reverse-complement the Seed+m8 column."
        )
    out = set()
    with SEED_TSV.open() as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("target_site"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3 and parts[2] in tiers:
                out.add(parts[0])
    return out


def make_checker(seed_sites: set[str], forbid_internal_atg: bool = True):
    """Constraint predicate. `forbid_internal_atg` fixes a real gap: the shipped check is
    `bc.startswith("ATG")`, i.e. position 0 only, so ~26% of passing barcodes still carry an
    internal ATG despite the README calling this 'start-ATG avoidance'."""

    def ok(bc: str) -> bool:
        gc = (bc.count("G") + bc.count("C")) / len(bc)
        if not (GC_MIN <= gc <= GC_MAX):
            return False
        run = 1
        for i in range(1, len(bc)):
            run = run + 1 if bc[i] == bc[i - 1] else 1
            if run > MAX_HOMOPOLYMER:
                return False
        if forbid_internal_atg:
            if "ATG" in bc:
                return False
        elif bc.startswith("ATG"):
            return False
        if _POLYA in bc:
            return False
        if any(s in bc for s in _SITES):
            return False
        return not any(s in bc for s in seed_sites)

    return ok


def _pack(s: str) -> int:
    """2-bit encode, so deletion neighbourhoods live in a set of ints not strings."""
    v = 0
    for ch in s:
        v = (v << 2) | (0 if ch == "A" else 1 if ch == "C" else 2 if ch == "G" else 3)
    return v


def _deletion_neighbourhood(bc: str, k: int) -> list[int]:
    n = len(bc)
    return [
        _pack("".join(bc[i] for i in range(n) if i not in drop))
        for drop in itertools.combinations(range(n), k)
    ]


def build_code(length: int, min_distance: int, ok, rng, target: int, max_tries: int):
    """Greedy construction with a guaranteed min Levenshtein distance. Returns
    (n_accepted, n_candidates_drawn)."""
    k = min_distance - 1
    seen: set[int] = set()
    accepted = 0
    tries = 0
    while accepted < target and tries < max_tries:
        tries += 1
        bc = "".join(rng.choice("ACGT") for _ in range(length))
        if not ok(bc):
            continue
        if k == 0:
            nb = [_pack(bc)]
        else:
            nb = _deletion_neighbourhood(bc, k)
        if seen.isdisjoint(nb):
            seen.update(nb)
            accepted += 1
    return accepted, tries


def main() -> None:
    setup_plotting_style()
    rng = random.Random(SEED)

    # ---------- 1. what the seed filter costs ----------
    print("=== 1. miRNA seed-site exclusion (TargetScanHuman 8.0, human) ===")
    print(f"current recommended barcode length: {CURRENT_LENGTH} nt\n")
    print(
        f"{'tier set':20s} {'sites':>7s} {'pass frac':>11s} {'effective space':>17s} {'vs 20k':>10s}"
    )
    pass_frac = {}
    for name, tiers in TIER_SETS.items():
        sites = load_seed_sites(tiers)
        ok = make_checker(sites)
        n_mc = 60_000
        hits = sum(
            ok("".join(rng.choice("ACGT") for _ in range(CURRENT_LENGTH))) for _ in range(n_mc)
        )
        p = hits / n_mc
        pass_frac[name] = p
        eff = p * 4.0**CURRENT_LENGTH
        print(f"{name:20s} {len(sites):>7,} {p:>11.4f} {eff:>17.3e} {eff / TARGET_N:>9.2e}x")
    print(
        "\n  Every option leaves >1e7x headroom over a 20,000-clone library, so the seed\n"
        "  filter is effectively free on collision grounds. 'all_human' over-filters on the\n"
        "  strength of ~1,600 unassigned-tier miRNAs that mostly are not expressed."
    )

    # ---------- 2. internal-ATG gap in the shipped checker ----------
    sites_def = load_seed_sites(TIER_SETS[DEFAULT_TIERS])
    ok_loose = make_checker(sites_def, forbid_internal_atg=False)
    n_mc = 60_000
    loose = [
        bc
        for bc in ("".join(rng.choice("ACGT") for _ in range(CURRENT_LENGTH)) for _ in range(n_mc))
        if ok_loose(bc)
    ]
    with_atg = sum("ATG" in bc for bc in loose)
    print("\n=== 2. the shipped ATG check is position-0 only ===")
    print(
        f"  of {len(loose):,} barcodes passing the shipped constraints, "
        f"{with_atg:,} ({100 * with_atg / len(loose):.1f}%) contain an INTERNAL ATG"
    )

    # ---------- 3. can we shorten to free bases for a UMI? ----------
    print(f"\n=== 3. edit-distance capacity: can {TARGET_N:,} barcodes fit in fewer bases? ===")
    print("  (achieved counts are LOWER bounds -- the filter over-rejects; see module docstring)")
    ok = make_checker(sites_def)
    rows = []
    print(
        f"\n{'len':>4s} {'d':>3s} {'achieved':>10s} {'target met':>11s} {'UMI bases free':>15s} {'UMI codes':>10s}"
    )
    results = {}
    for length in LENGTHS:
        for d in DISTANCES:
            target = TARGET_N
            max_tries = 60 * TARGET_N
            got, tries = build_code(length, d, ok, rng, target, max_tries)
            free = CURRENT_LENGTH - length
            met = got >= TARGET_N
            results[(length, d)] = got
            print(
                f"{length:>4d} {d:>3d} {got:>10,} {'YES' if met else 'no':>11s} "
                f"{free:>15d} {4**free if free else 0:>10,}"
            )
            rows.append(
                dict(
                    barcode_length=length,
                    min_levenshtein=d,
                    barcodes_achieved=got,
                    target=TARGET_N,
                    target_met=int(met),
                    candidates_drawn=tries,
                    umi_bases_freed=free,
                    umi_codes=4**free if free else 0,
                    seed_tier_excluded=DEFAULT_TIERS,
                )
            )

    out = HERE / "barcode_design_capacity.tsv"
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {out}")

    # ---------- 4. the UMI itself ----------
    print("\n=== 4. what the freed bases can actually do ===")
    print(f"  ORF identity needs >= {math.ceil(math.log(2243, 4))} nt RAW to index 2,243 ORFs")
    print("  -- so 3-4 nt CANNOT encode ORF identity; it can only serve as a molecular UMI.")
    for n in (3, 4, 5, 6):
        print(f"    {n} nt -> {4**n:>6,} molecular UMIs per barcode")

    # ---------- 5. UMI width, for the S1E6 optional addition ----------
    # A vector-encoded degenerate field counts independent integration events; it cannot
    # collapse PCR/RT duplicates, because every molecule amplified from one provirus carries
    # the same vector sequence. Either way the width is set by how many distinct items share
    # a single barcode (M) and how much undercounting is tolerable.
    print("\n=== 5. UMI width vs items sharing one barcode (undercount %) ===")
    Ms = [13, 50, 100, 300, 1000, 3000]
    print(f"{'UMI nt':>7} {'space':>11} " + " ".join(f"{'M=' + str(m):>9}" for m in Ms))
    umi_rows = []
    for k in range(5, 13):
        S = 4**k
        cells = []
        for M in Ms:
            under = 1 - (S * (1 - math.exp(-M / S))) / M
            cells.append(under)
            umi_rows.append(
                dict(
                    umi_nt=k, umi_space=S, items_per_barcode=M, undercount_pct=round(100 * under, 3)
                )
            )
        print(f"{k:>7} {S:>11,} " + " ".join(f"{100 * c:8.2f}%" for c in cells))

    print("\n  largest M holding undercount under 1%:")
    for k in (6, 7, 8, 9, 10, 12):
        S = 4**k
        lo, hi = 1, 10**7
        while lo < hi:
            mid = (lo + hi + 1) // 2
            lo, hi = (mid, hi) if 1 - (S * (1 - math.exp(-mid / S))) / mid < 0.01 else (lo, mid - 1)
        print(f"    {k:>2} nt -> M = {lo:,}")

    # Placement: bases inserted BETWEEN the ORF and the barcode lengthen the
    # template-switching span; bases placed 3' of the barcode do not.
    P_SWITCH, SPAN = 5.3e-4, 177
    print(f"\n  placement, template-switching mispairing (span {SPAN} bp today):")
    for lab, extra in (
        ("3' of the barcode (distal)", 0),
        ("between ORF and barcode, 6 nt", 6),
        ("between ORF and barcode, 8 nt", 8),
        ("between ORF and barcode, 10 nt", 10),
    ):
        sp = SPAN + extra
        print(f"    {lab:34s} span {sp:>4} bp -> {100 * (1 - (1 - P_SWITCH) ** sp):5.2f}%")

    # Local vs global uniqueness. An ORF-specific barcode partitions the library, so the
    # appended field only needs to be unique WITHIN one ORF (M = representation depth).
    # A single tag carrying the whole job would have to stay unique across every cell.
    N_ORF_CEILING, REPRESENTATION = 30_000, 500
    total_cells = N_ORF_CEILING * REPRESENTATION

    def min_width(M, bar=0.01):
        for k in range(1, 25):
            S = 4**k
            if 1 - (S * (1 - math.exp(-M / S))) / M < bar:
                return k
        return None

    print("\n  local vs global uniqueness (width for <1% undercount):")
    print(
        f"    one global tag over all {total_cells:,} cells        -> {min_width(total_cells):>2} nt"
    )
    print(
        f"    ORF-specific barcode + local field, M={REPRESENTATION:<4}      -> "
        f"{min_width(REPRESENTATION):>2} nt"
    )
    print("    -> partitioning by ORF is what makes the appended field cheap.")

    print(f"\n  plasmid-pool complexity caps resolvable clones (cells/pair = {REPRESENTATION}/B):")
    for B in (1, 10, 50, 100, 500):
        print(
            f"    B={B:>4} field values per ORF -> {N_ORF_CEILING * B:>12,} plasmid clones, "
            f"{REPRESENTATION / B:>6.1f} cells share a pair"
        )

    # Two statistics for the same collisions. Ambiguous fraction = share of clones sharing a
    # tag with >=1 other; undercount = share of distinct clones lost. undercount ~ ambig/2.
    print("\n  ambiguous fraction vs undercount at M=500 (same collisions, different question):")
    for k in (5, 6, 7, 8, 9, 10):
        S = 4**k
        M = REPRESENTATION
        ambig = 1 - math.exp(-(M - 1) / S)
        under = 1 - (S * (1 - math.exp(-M / S))) / M
        print(
            f"    {k:>2} nt  ambiguous {100 * ambig:6.2f}%   undercount {100 * under:6.2f}%   "
            f"colliding pairs {M * (M - 1) / 2 / S:6.1f}"
        )

    # A degenerate window cannot be filtered the way the designed barcode is: every member of
    # N-k is in the pool. Shortening is the only lever, and it acts on all three liabilities.
    print("\n  liabilities a DEGENERATE field cannot be filtered for:")
    print(f"    {'nt':>4} {'homopolymer>3':>14} {'BbsI/Esp3I':>12} {'conserved seed':>15}")
    seed7 = load_seed_sites({"broadly_conserved", "conserved"})
    for k in range(4, 11):
        n_mc, hp, st, sd = 200_000, 0, 0, 0
        for _ in range(n_mc):
            f = "".join(rng.choice("ACGT") for _ in range(k))
            run = mx = 1
            for i in range(1, k):
                run = run + 1 if f[i] == f[i - 1] else 1
                mx = max(mx, run)
            hp += mx > 3
            st += any(x in f for x in _SITES)
            sd += any(f[i : i + 7] in seed7 for i in range(max(0, k - 6)))
        print(
            f"    {k:>4} {100 * hp / n_mc:13.2f}% {100 * st / n_mc:11.2f}% {100 * sd / n_mc:14.2f}%"
        )
    print("    a 6-nt window cannot contain a 7-mer seed internally; junction-spanning seeds")
    print("    survive at any width and are set by the flanks, so check those separately.")

    umi_out = HERE / "umi_width.tsv"
    with umi_out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(umi_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(umi_rows)
    print(f"\nwrote {umi_out}")

    # ---------- figure ----------
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.4))
    ax = axes[0]
    for i, d in enumerate(DISTANCES):
        ys = [results[(L, d)] for L in LENGTHS]
        ax.plot(
            LENGTHS,
            ys,
            marker="o",
            color=CATEGORICAL_PALETTE[i % len(CATEGORICAL_PALETTE)],
            label=f"min Levenshtein {d}" + (" (corrects 1 edit)" if d == 3 else ""),
        )
    ax.axhline(TARGET_N, ls="--", color="0.35", lw=1.6)
    ax.annotate(
        f"target {TARGET_N:,}",
        xy=(LENGTHS[0], TARGET_N),
        xytext=(2, 6),
        textcoords="offset points",
        fontsize=9,
        color="0.35",
    )
    ax.axvline(CURRENT_LENGTH, ls=":", color="0.5", lw=1.4)
    ax.set_yscale("log")
    ax.set_xlabel("barcode length (nt)")
    ax.set_ylabel("barcodes built (lower bound)")
    ax.legend(frameon=False, fontsize=8, loc="lower right")

    ax = axes[1]
    names = list(TIER_SETS)
    vals = [100 * pass_frac[n] for n in names]
    ax.bar(range(len(names)), vals, 0.6, color=COLORS["primary"])
    for i, v in enumerate(vals):
        ax.text(i, v + 1.0, f"{v:.0f}%", ha="center", fontsize=9)
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels([n.replace("_", "\n") for n in names], fontsize=9)
    ax.set_ylabel(f"% of random {CURRENT_LENGTH}-mers passing all constraints")
    ax.set_ylim(0, max(vals) * 1.2)
    ax.annotate(
        "miRNA seed tiers excluded",
        xy=(0.5, -0.20),
        xycoords="axes fraction",
        ha="center",
        fontsize=8,
        color="0.4",
    )

    save_figure(fig, "s1e2_barcode_design", output_dir=str(HERE / "figures"))
    print("wrote figures/s1e2_barcode_design.{svg,png}")


if __name__ == "__main__":
    main()
