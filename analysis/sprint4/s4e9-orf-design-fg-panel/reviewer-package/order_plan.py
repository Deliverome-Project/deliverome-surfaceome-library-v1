"""How the L1 library is ordered: plate layout, volumes, and length-weighted pooling amounts.

Every active L1 member is ordered as purified plasmid in 96-well V-bottom plates, in order of
synthesized-insert length (shortest first), so a plate holds members of similar length. The
order number (``L1-O0001`` ...) follows that order; the member's ``candidate_id`` is kept, since
every library record is keyed on it.

Amounts (Becca, 2026-11-10): 4 ug ordered per member, 3 ug assumed delivered, normalized to
100 ng/uL, so 30 uL per member, split into three plate sets with the same layout:

| Plate set | Volume per well | DNA at 100 ng/uL |
|---|---:|---:|
| Pooling A, B | 5 uL each | 0.5 ug each |
| Stock | 20 uL | 2.0 ug |

The vendor may deliver below 100 ng/uL. The volumes are fixed, so a lower concentration lowers
the DNA in every well in proportion; the summary reports each plate set's DNA, and each pool's
total, at 100, 75 and 50 ng/uL. Pooling weights are relative, so the members needing extra DNA
do not change with concentration.

Length-weighted pooling follows Moon et al. 2026 (Nat Commun, PMC13319239), who pooled
``0.2 * 3.5^(L/1000)`` pmol of each plasmid, L the effector length, to offset lower lentiviral
titer for longer genomes (Sweeney & Vink 2021, PMC8141603: 3-4x lower functional titer per
extra kb of packaged genome). Our backbone is constant, so differences in insert length are
differences in packaged length, and the weight per member is ``factor^((L - L_ref)/1000)`` in
pmol, converted to mass with the member's own vector length. The reference is set so that a
member with a 4.3-kb insert needs exactly one pooling well (500 ng) per pool; members above
4.3 kb need more than a well holds, which is what the extra order covers. Shorter members need
less than a well, down to nanoliters, which is a pipetting question the plan reports.

The per-kb factor is uncertain for our length range (Sweeney & Vink's data are shallower,
~2x/kb, below ~7 kb packaged), so the plan is computed at 3.5 (Moon et al.) and at 2.0 and 2.5
for comparison.

    uv run --locked python order_plan.py      # -> l1_order_plan.tsv, l1_order_plan_summary.json, figure
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

HERE = Path(__file__).parent
S4E9 = HERE.parent
MANIFEST = S4E9 / "candidate_manifest_tranches.tsv"
OUT = HERE / "l1_order_plan.tsv"
SUMMARY = HERE / "l1_order_plan_summary.json"
PLATES = HERE / "l1_order_pooling_plates.tsv"
FIGURES = HERE / "figures"

ORDERED_UG, ASSUMED_UG, CONC_NG_UL = 4.0, 3.0, 100.0
STOCK_UL, POOL_UL, N_POOLS = 20.0, 5.0, 2
CONC_SCENARIOS = (100.0, 75.0, 50.0)  # ng/uL: assumed, and two lower deliveries
EXTRA_THRESHOLD_NT = 4300  # members whose insert is longer get extra DNA for pooling
FACTOR = 3.5  # fold lower functional titer per extra kb (Moon et al. 2026)
SENSITIVITY = (2.0, 2.5, 3.5)
WELLS = [f"{r}{c:02d}" for c in range(1, 13) for r in "ABCDEFGH"]  # column-major: A01, B01, ...
NG_PER_PMOL_PER_BP = 0.65  # 650 g/mol per bp of double-stranded DNA
MIN_HAND_UL = 0.5  # smallest volume worth pipetting by hand
#: Plates whose members are pooled well by well at their own volume rather than at one
#: per-plate dilution: the long tail, where the weight spread within a plate is widest and
#: the members needing extra DNA sit.
PER_WELL_PLATES = 2


def members() -> list[dict]:
    with MANIFEST.open() as fh:
        rows = [
            r
            for r in csv.DictReader(fh, delimiter="\t")
            if r["library"] == "l1" and not r["retired_reason"]
        ]
    return sorted(rows, key=lambda r: (int(float(r["insert_nt"])), r["candidate_id"]))


def pool_ng(insert_nt: int, vector_bp: int, factor: float, ref: tuple[int, int]) -> float:
    """ng of a member per pool so that pmol follows factor^(L/1000), normalized to ``ref``.

    ``ref`` is (insert_nt, vector_bp) of a member that gets exactly one pooling well.
    """
    ref_insert, ref_vector = ref
    ref_pmol = POOL_UL * CONC_NG_UL / (ref_vector * NG_PER_PMOL_PER_BP)
    pmol = ref_pmol * factor ** ((insert_nt - ref_insert) / 1000)
    return pmol * vector_bp * NG_PER_PMOL_PER_BP


def plan(rows: list[dict]) -> tuple[list[dict], dict]:
    available_ul = ASSUMED_UG * 1000 / CONC_NG_UL
    spare_ul = available_ul - STOCK_UL - N_POOLS * POOL_UL
    assert spare_ul >= 0, "plate sets exceed the assumed delivery"
    # Reference: a 4.3-kb insert in a vector of the median length for inserts near 4.3 kb.
    near = sorted(
        int(float(r["vector_bp"]))
        for r in rows
        if abs(int(float(r["insert_nt"])) - EXTRA_THRESHOLD_NT) <= 200
    )
    ref = (EXTRA_THRESHOLD_NT, near[len(near) // 2])
    out = []
    for i, r in enumerate(rows):
        ins, vec = int(float(r["insert_nt"])), int(float(r["vector_bp"]))
        row = {
            "order_id": f"L1-O{i + 1:04d}",
            "plate": i // 96 + 1,
            "well": WELLS[i % 96],
            "candidate_id": r["candidate_id"],
            "key": r["key"],
            "set": r["set"],
            "insert_nt": ins,
            "vector_bp": vec,
        }
        for f in SENSITIVITY:
            ng = pool_ng(ins, vec, f, ref)
            row[f"pool_ng_per_pool_x{f}"] = round(ng, 2)
            row[f"pool_ul_per_pool_x{f}"] = round(ng / CONC_NG_UL, 4)
            # DNA beyond the standard pooling wells, ordered on top of the 4 ug.
            extra_ug = max(0.0, N_POOLS * (ng - POOL_UL * CONC_NG_UL)) / 1000
            row[f"extra_ug_assumed_x{f}"] = round(extra_ug, 3)
            row[f"order_ug_x{f}"] = round(ORDERED_UG + extra_ug * ORDERED_UG / ASSUMED_UG, 2)
        out.append(row)

    def stats(f: float) -> dict:
        ul = [r[f"pool_ul_per_pool_x{f}"] for r in out]
        extra = [r for r in out if r[f"extra_ug_assumed_x{f}"] > 0]
        longest, shortest = out[-1], out[0]
        return {
            "members_needing_extra": len(extra),
            "extra_ug_assumed_total": round(sum(r[f"extra_ug_assumed_x{f}"] for r in extra), 1),
            "order_ug_max": max(r[f"order_ug_x{f}"] for r in out),
            "pool_ul_shortest": shortest[f"pool_ul_per_pool_x{f}"],
            "pool_ul_median": sorted(ul)[len(ul) // 2],
            "pool_ul_longest": longest[f"pool_ul_per_pool_x{f}"],
            "members_below_hand_pipetting": sum(v < MIN_HAND_UL for v in ul),
            "pool_ug_total_per_pool": round(
                sum(r[f"pool_ng_per_pool_x{f}"] for r in out) / 1000, 1
            ),
            "longest_vs_shortest_pmol": round(
                f ** ((longest["insert_nt"] - shortest["insert_nt"]) / 1000), 1
            ),
            "longest_vs_median_pmol": round(
                f ** ((longest["insert_nt"] - out[len(out) // 2]["insert_nt"]) / 1000), 1
            ),
        }

    summary = {
        "members": len(out),
        "plates_per_set": math.ceil(len(out) / 96),
        "plate_sets": {
            "stock_ul": STOCK_UL,
            "pooling_ul_each": POOL_UL,
            "pooling_sets": N_POOLS,
            "spare_ul": spare_ul,
        },
        "concentration_scenarios": {
            f"{c:g}_ng_ul": {
                "pooling_well_ug": round(POOL_UL * c / 1000, 3),
                "stock_well_ug": round(STOCK_UL * c / 1000, 3),
                "delivered_ug": round(available_ul * c / 1000, 3),
                "pool_ug_total_per_pool_x3.5": round(
                    sum(r["pool_ng_per_pool_x3.5"] for r in out) * c / CONC_NG_UL / 1000, 1
                ),
            }
            for c in CONC_SCENARIOS
        },
        "ordered_ug": ORDERED_UG,
        "assumed_ug": ASSUMED_UG,
        "conc_ng_ul": CONC_NG_UL,
        "insert_nt_min_median_max": [
            out[0]["insert_nt"],
            out[len(out) // 2]["insert_nt"],
            out[-1]["insert_nt"],
        ],
        "vector_bp_min_max": [min(r["vector_bp"] for r in out), max(r["vector_bp"] for r in out)],
        "members_insert_over_threshold": sum(r["insert_nt"] > EXTRA_THRESHOLD_NT for r in out),
        "reference_member": {"insert_nt": ref[0], "vector_bp": ref[1]},
        "weighting": {f"x{f}": stats(f) for f in SENSITIVITY},
    }
    return out, summary


def pooling_plates(rows: list[dict]) -> list[dict]:
    """Per pooling plate: one dilution so a fixed 5 uL from every well delivers the weighting.

    Plates are length-ordered, so each holds a narrow length range and a single concentration
    (the plate median's weighted mass per 5 uL) approximates every member's target; the
    within-plate error is the spread of the members' targets around it. The last
    ``PER_WELL_PLATES`` plates are pooled per well at each member's own volume instead.
    """
    n_plates = max(r["plate"] for r in rows)
    out = []
    for p in range(1, n_plates + 1):
        pr = [r for r in rows if r["plate"] == p]
        row = {
            "plate": p,
            "members": len(pr),
            "insert_nt_min": pr[0]["insert_nt"],
            "insert_nt_max": pr[-1]["insert_nt"],
            "mode": "per_well" if p > n_plates - PER_WELL_PLATES else "plate_dilution",
        }
        for f in SENSITIVITY:
            ng = sorted(r[f"pool_ng_per_pool_x{f}"] for r in pr)
            med = ng[len(ng) // 2]
            row[f"conc_ng_ul_x{f}"] = round(min(med / POOL_UL, CONC_NG_UL), 3)
            row[f"max_error_fold_x{f}"] = round(max(ng[-1] / med, med / ng[0]), 2)
        out.append(row)
    return out


def figure(rows: list[dict]) -> None:
    from deliverome_analysis.plotting import (
        create_figure,
        get_color,
        save_figure,
        setup_plotting_style,
    )

    setup_plotting_style()
    fig, ax = create_figure(figsize=(11, 6.5))
    x = [r["insert_nt"] for r in rows]
    for f, color in (
        (3.5, get_color("primary")),
        (2.5, get_color("accent")),
        (2.0, get_color("secondary")),
    ):
        ax.scatter(
            x,
            [r[f"pool_ul_per_pool_x{f}"] for r in rows],
            s=6,
            color=color,
            label=f"{f}\u00d7 lower titer per kb",
            alpha=0.8,
        )
    ax.axhline(POOL_UL, color=get_color("neutral"), ls="--", lw=1.5)
    ax.text(300, POOL_UL * 1.15, "one pooling well (5 \u00b5L)", color=get_color("neutral"))
    ax.axhline(MIN_HAND_UL, color=get_color("neutral"), ls=":", lw=1.5)
    ax.text(300, MIN_HAND_UL * 1.15, "0.5 \u00b5L", color=get_color("neutral"))
    ax.axvline(EXTRA_THRESHOLD_NT, color=get_color("neutral"), ls="--", lw=1)
    ax.set_yscale("log")
    ax.set_xlabel("synthesized insert length (nt)")
    ax.set_ylabel("\u00b5L per pool at 100 ng/\u00b5L")
    ax.legend(loc="lower right", bbox_to_anchor=(0.76, 0.0), frameon=False, markerscale=4)
    save_figure(fig, "l1_pooling_volume_by_length", output_dir=FIGURES)


def main() -> int:
    rows, summary = plan(members())
    with OUT.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    plates = pooling_plates(rows)
    with PLATES.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(plates[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(plates)
    dil = [p for p in plates if p["mode"] == "plate_dilution"]
    summary["pooling_plates"] = {
        f"x{f}": {
            "plate_dilution_plates": len(dil),
            "conc_ng_ul_range": [
                min(p[f"conc_ng_ul_x{f}"] for p in dil),
                max(p[f"conc_ng_ul_x{f}"] for p in dil),
            ],
            "max_error_fold_plate_dilution": max(p[f"max_error_fold_x{f}"] for p in dil),
            "max_error_fold_excluding_plate_1": max(p[f"max_error_fold_x{f}"] for p in dil[1:]),
            "per_well_members": sum(p["members"] for p in plates if p["mode"] == "per_well"),
        }
        for f in SENSITIVITY
    }
    SUMMARY.write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps(summary, indent=1))
    figure(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
