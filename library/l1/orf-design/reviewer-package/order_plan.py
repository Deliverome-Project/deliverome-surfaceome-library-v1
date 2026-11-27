"""How the L1 library is ordered: plate layout, volumes, and length-weighted pooling amounts.

Every built L1 insert is ordered as purified plasmid in 96-well V-bottom plates, in order of
ordered (tagged, synthesis-recoded) insert length (shortest first), so a plate holds inserts of
similar length; held members have no insert and are not ordered. A member may be ordered as
several inserts (an internal member's 2 or 3 arms; Becca, 2026-11-14/15), so the plate
position belongs to the insert: its insert id (``l1-i0001`` ...) is assigned by
``renumber_library_members.py`` (``l1_insert_numbering.tsv``) with the same sort, and this plan
asserts the two agree. Each row carries the insert's permanent arm id (``l1-pNNNN-aK``) and the
member's design id (``l1-pNNNN``). Amounts below are per insert (per well).

Amounts (decided by Becca Carlson, 2026-11-07, plate sets changed 2026-11-10: "5ul x 2
plates + 20ul. assumed at 100ng/ul but could be lower"; implemented by the Claude Code
session l1-reviewer-package, PRs 243 and 283): 4 ug ordered per member, 3 ug assumed delivered, normalized to
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
NUMBERING = S4E9 / "l1_insert_numbering.tsv"
COMPLIANCE = S4E9 / "synthesis_compliance.tsv"
OUT = HERE / "l1_order_plan.tsv"
SUMMARY = HERE / "l1_order_plan_summary.json"
PLATES = HERE / "l1_order_pooling_plates.tsv"
FIGURES = HERE / "figures"

ORDERED_UG, ASSUMED_UG, CONC_NG_UL = 4.0, 3.0, 100.0
STOCK_UL, POOL_UL, N_POOLS = 20.0, 5.0, 2
CONC_SCENARIOS = (100.0, 75.0, 50.0)  # ng/uL: assumed, and two lower deliveries
EXTRA_THRESHOLD_NT = 4300  # members whose insert is longer get extra DNA for pooling
FACTOR = 3.5  # fold lower functional titer per extra kb (Moon et al. 2026)
SENSITIVITY = tuple(sorted({2.0, 2.5, FACTOR}))
WELLS = [f"{r}{c:02d}" for c in range(1, 13) for r in "ABCDEFGH"]  # column-major: A01, B01, ...
NG_PER_PMOL_PER_BP = 0.65  # 650 g/mol per bp of double-stranded DNA
MIN_HAND_UL = 0.5  # smallest volume worth pipetting by hand
#: Plates whose members are pooled well by well at their own volume rather than at one
#: per-plate dilution: the long tail, where the weight spread within a plate is widest and
#: the members needing extra DNA sit.
PER_WELL_PLATES = 2


def verify_numbering_matches_compliance(
    numbered: list[dict[str, str]], compliance: list[dict[str, str]]
) -> None:
    """Fail closed if ``l1_insert_numbering.tsv`` and ``synthesis_compliance.tsv`` disagree
    on the set of ``(arm_id, insert_nt)`` pairs."""

    def _index(rows: list[dict[str, str]], length_col: str, label: str) -> dict[str, int]:
        out: dict[str, int] = {}
        for row in rows:
            arm_id = row["arm_id"].strip()
            if not arm_id:
                raise ValueError(f"{label} contains a row with an empty arm_id")
            if arm_id in out:
                raise ValueError(f"duplicate arm_id {arm_id!r} in {label}")
            try:
                nt = int(row[length_col])
            except (KeyError, TypeError, ValueError) as exc:
                raise ValueError(
                    f"invalid {length_col}={row.get(length_col)!r} for {arm_id!r} in {label}"
                ) from exc
            out[arm_id] = nt
        return out

    numbered_by_arm = _index(numbered, "ordered_insert_nt", "l1_insert_numbering.tsv")
    compliance_by_arm = _index(compliance, "insert_nt", "synthesis_compliance.tsv")

    missing = sorted(set(compliance_by_arm) - set(numbered_by_arm))
    extra = sorted(set(numbered_by_arm) - set(compliance_by_arm))
    length_mismatches = sorted(
        (arm, numbered_by_arm[arm], compliance_by_arm[arm])
        for arm in set(numbered_by_arm) & set(compliance_by_arm)
        if numbered_by_arm[arm] != compliance_by_arm[arm]
    )
    if missing or extra or length_mismatches:
        parts: list[str] = []
        if missing:
            parts.append(f"missing from l1_insert_numbering.tsv ({len(missing)}): {missing[:5]}")
        if extra:
            parts.append(f"extra in l1_insert_numbering.tsv ({len(extra)}): {extra[:5]}")
        if length_mismatches:
            parts.append(
                f"insert_nt mismatch vs synthesis_compliance.tsv ({len(length_mismatches)}): "
                f"{length_mismatches[:5]}"
            )
        raise ValueError(
            "l1_insert_numbering.tsv does not match synthesis_compliance.tsv; "
            "rerun renumber_library_members.py --write and synthesis_compliance.py ("
            + "; ".join(parts)
            + ")"
        )


def inserts() -> list[dict]:
    """Ordered L1 inserts, in plate order, with their ordered insert and vector lengths.

    Plates and pooling weights use the ordered insert: tagged and synthesis-recoded, as in the
    order of record (``synthesis_compliance.tsv``), not the untagged draft length in the
    manifest (decided 2026-11-14 by the coordinator under Becca's delegation). The vector grows
    by the same number of bases. Held members have no insert and are not ordered.
    """
    with NUMBERING.open() as fh:
        numbered = list(csv.DictReader(fh, delimiter="\t"))
    with COMPLIANCE.open() as fh:
        compliance = list(csv.DictReader(fh, delimiter="\t"))
    verify_numbering_matches_compliance(numbered, compliance)
    with MANIFEST.open() as fh:
        member = {
            r["design_id"]: r
            for r in csv.DictReader(fh, delimiter="\t")
            if r["library"] == "l1" and not r["retired_reason"]
        }
    rows = []
    for n in numbered:
        m = member[n["design_id"]]
        ordered = int(n["ordered_insert_nt"])
        rows.append(
            {
                **n,
                "set": m["set"],
                "insert_nt": ordered,
                "vector_bp": int(float(m["vector_bp"])) + ordered - int(float(m["insert_nt"])),
            }
        )
    order = sorted(rows, key=lambda r: (r["insert_nt"], r["arm_id"]))
    # The insert id is the plate position; a mismatch means the renumber stage is stale.
    assert [r["insert_id"] for r in order] == [f"l1-i{i:04d}" for i in range(1, len(rows) + 1)], (
        "insert ids are not in plate order; rerun renumber_library_members.py --write"
    )
    return order


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
        assert (int(r["plate"]), r["well"]) == (i // 96 + 1, WELLS[i % 96]), r["insert_id"]
        row = {
            "plate": i // 96 + 1,
            "well": WELLS[i % 96],
            "insert_id": r["insert_id"],
            "arm_id": r["arm_id"],
            "design_id": r["design_id"],
            "rung": r["rung"],
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
            "inserts_needing_extra": len(extra),
            "extra_ug_assumed_total": round(sum(r[f"extra_ug_assumed_x{f}"] for r in extra), 1),
            "order_ug_max": max(r[f"order_ug_x{f}"] for r in out),
            "pool_ul_shortest": shortest[f"pool_ul_per_pool_x{f}"],
            "pool_ul_median": sorted(ul)[len(ul) // 2],
            "pool_ul_longest": longest[f"pool_ul_per_pool_x{f}"],
            "inserts_below_hand_pipetting": sum(v < MIN_HAND_UL for v in ul),
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
        "inserts": len(out),
        "members": len({r["design_id"] for r in out}),
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
                f"pool_ug_total_per_pool_x{FACTOR}": round(
                    sum(r[f"pool_ng_per_pool_x{FACTOR}"] for r in out) * c / CONC_NG_UL / 1000, 1
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
        "inserts_over_threshold": sum(r["insert_nt"] > EXTRA_THRESHOLD_NT for r in out),
        "reference_insert": {"insert_nt": ref[0], "vector_bp": ref[1]},
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
            "inserts": len(pr),
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


def subpool_strategies(rows: list[dict], *, ks: tuple[int, ...] = (1, 2, 4)) -> dict:
    """Direct 100 ng/uL subplate pooling (2 uL/well into K contiguous subpools per plate).

    Each 96-well plate is split into ``K`` contiguous length-ordered column blocks (12, 6, or 3
    columns). Every well contributes a uniform 2 uL directly from the 100 ng/uL stock (no per-well
    dilution); each subpool tube is then pipetted into the master pool at its geometric-mean target
    weight (minimizing max fold error vs. ideal). Subpool tubes needing < 2.0 uL into the master
    pool are diluted 1:10 first so every transfer stays >= 2 uL.
    """
    n_plates = max(r["plate"] for r in rows)
    result: dict[str, dict] = {}
    for k in ks:
        by_factor: dict[str, dict] = {}
        for f in SENSITIVITY:
            col_ng = f"pool_ng_per_pool_x{f}"
            ratios: list[float] = []
            vols: list[float] = []
            worst_plate, worst_err = 0, 0.0
            for p in range(1, n_plates + 1):
                pr = [r for r in rows if r["plate"] == p]
                chunks = [pr[i * len(pr) // k : (i + 1) * len(pr) // k] for i in range(k)]
                for ch in chunks:
                    if not ch:
                        continue
                    ngs = [float(r[col_ng]) for r in ch]
                    gmean = math.sqrt(min(ngs) * max(ngs))
                    vols.append(gmean * len(ch) / CONC_NG_UL)
                    for v in ngs:
                        ratio = gmean / v  # realized / ideal
                        ratios.append(ratio)
                        err = max(ratio, 1.0 / ratio)
                        if err > worst_err:
                            worst_err, worst_plate = err, p
            abs_err = sorted(max(x, 1.0 / x) for x in ratios)
            by_factor[f"x{f}"] = {
                "subpool_tubes": len(vols),
                "worst_from_target_fold": round(worst_err, 2),
                "worst_plate": worst_plate,
                "total_max_min_skew_fold": round(max(ratios) / min(ratios), 2),
                "p95_from_target_fold": round(abs_err[int(0.95 * len(abs_err))], 2),
                "master_pool_ul_min": round(min(vols), 2),
                "master_pool_ul_max": round(max(vols), 1),
                "tubes_below_2ul": sum(v < 2.0 for v in vols),
            }
        result[f"k{k}_per_plate"] = by_factor
    return result


def figure(rows: list[dict]) -> None:
    from deliverome_analysis.core.plotting import (
        create_figure,
        save_figure,
        setup_plotting_style,
    )

    setup_plotting_style()
    fig, (ax1, ax2) = create_figure(figsize=(13.2, 5.4), ncols=2)

    lengths_kb = [r["insert_nt"] / 1000.0 for r in rows]
    med_kb = lengths_kb[len(lengths_kb) // 2]
    n_plates = max(r["plate"] for r in rows)

    # Helper: realized/ideal representation ratio per construct for K subpools/plate at factor f
    def _subpool_ratios(k: int, f: float) -> list[float]:
        col_ng = f"pool_ng_per_pool_x{f}"
        out_r: list[float] = []
        for p in range(1, n_plates + 1):
            pr = [r for r in rows if r["plate"] == p]
            chunks = [pr[i * len(pr) // k : (i + 1) * len(pr) // k] for i in range(k)]
            for ch in chunks:
                if not ch:
                    continue
                ngs = [float(r[col_ng]) for r in ch]
                gmean = math.sqrt(min(ngs) * max(ngs))
                out_r.extend(gmean / v for v in ngs)
        return out_r

    # Panel A: Why calibrate slope first — expected library skew vs. ideal (true slope = 3.5x/kb)
    equimolar_rel = [FACTOR ** (-(x - med_kb)) for x in lengths_kb]
    overcorr_rel = [(FACTOR / 2.0) ** (x - med_kb) for x in lengths_kb]
    k2_ratios = _subpool_ratios(2, FACTOR)

    ax1.axhline(1.0, color="#6f5d5a", ls="--", lw=1.1, zorder=1)
    ax1.axhspan(1 / 1.35, 1.35, color="#e8f3ef", alpha=0.8, zorder=0)
    ax1.plot(
        lengths_kb,
        equimolar_rel,
        color="#b83b4a",
        lw=2.3,
        label="No length weighting (equimolar; 495\u00d7 span)",
        zorder=2,
    )
    ax1.plot(
        lengths_kb,
        overcorr_rel,
        color="#d48a29",
        lw=2.1,
        ls="--",
        label="Uncalibrated 3.5\u00d7/kb if true slope is 2.0\u00d7/kb (16\u00d7 overcorrection)",
        zorder=3,
    )
    ax1.scatter(
        lengths_kb,
        k2_ratios,
        s=6,
        color="#2c7a72",
        alpha=0.8,
        label="Calibrated 2 subpools/plate (1.68\u00d7 total max/min skew)",
        zorder=4,
    )
    ax1.set_yscale("log")
    ax1.set_ylim(0.006, 22)
    ax1.set_xlim(0.35, 5.65)
    ax1.set_xlabel("Synthesized insert length (kb)", fontsize=12)
    ax1.set_ylabel("Expected representation vs. ideal (fold)", fontsize=12)
    ax1.text(
        0.0,
        1.03,
        "A.  Expected library skew: uncalibrated slope vs. calibrated subplate pooling",
        transform=ax1.transAxes,
        fontsize=10.5,
        fontweight="bold",
        color="#2a2422",
    )
    ax1.legend(
        loc="lower left",
        frameon=True,
        facecolor="white",
        framealpha=0.92,
        fontsize=8.6,
        markerscale=2.5,
    )

    # Panel B: Pooling granularity comparison (1 vs 2 vs 4 subpools/plate at 3.5x/kb, 100 ng/uL)
    plate_pos = [(r["plate"] - 1) + (WELLS.index(r["well"]) / 96.0) + 1.0 for r in rows]
    k1_ratios = _subpool_ratios(1, FACTOR)
    k4_ratios = _subpool_ratios(4, FACTOR)

    ax2.axhline(1.0, color="#6f5d5a", ls="--", lw=1.1, zorder=1)
    ax2.scatter(
        plate_pos,
        k1_ratios,
        s=8,
        color="#c8b9b3",
        alpha=0.7,
        label=f"1 pool/plate ({n_plates} tubes): worst {max(max(k1_ratios), 1 / min(k1_ratios)):.2f}\u00d7, total skew {max(k1_ratios) / min(k1_ratios):.2f}\u00d7",
        zorder=2,
    )
    ax2.scatter(
        plate_pos,
        k2_ratios,
        s=8,
        color="#d48a29",
        alpha=0.75,
        label=f"2 subpools/plate ({2 * n_plates} tubes): worst {max(max(k2_ratios), 1 / min(k2_ratios)):.2f}\u00d7, total skew {max(k2_ratios) / min(k2_ratios):.2f}\u00d7",
        zorder=3,
    )
    ax2.scatter(
        plate_pos,
        k4_ratios,
        s=7,
        color="#2c7a72",
        alpha=0.85,
        label=f"4 subpools/plate ({4 * n_plates} tubes): worst {max(max(k4_ratios), 1 / min(k4_ratios)):.2f}\u00d7, total skew {max(k4_ratios) / min(k4_ratios):.2f}\u00d7",
        zorder=4,
    )
    ax2.set_xlim(0.5, n_plates + 1.2)
    ax2.set_ylim(0.55, 1.78)
    ax2.set_xlabel("96-well stock plate (shortest \u2192 longest)", fontsize=12)
    ax2.set_ylabel("Realized / ideal DNA (fold, 3.5\u00d7/kb)", fontsize=12)
    ax2.text(
        0.0,
        1.03,
        "B.  Discretization error at 2 \u00b5L/well (100 ng/\u00b5L): 1 vs. 2 vs. 4 subpools/plate",
        transform=ax2.transAxes,
        fontsize=10.5,
        fontweight="bold",
        color="#2a2422",
    )
    ax2.legend(
        loc="upper left",
        frameon=True,
        facecolor="white",
        framealpha=0.92,
        fontsize=8.6,
        markerscale=2.5,
    )

    fig.tight_layout(w_pad=2.8)
    save_figure(fig, "l1_pooling_volume_by_length", output_dir=FIGURES)


def main() -> int:
    rows, summary = plan(inserts())
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
            "per_well_inserts": sum(p["inserts"] for p in plates if p["mode"] == "per_well"),
        }
        for f in SENSITIVITY
    }
    summary["subpool_strategies"] = subpool_strategies(rows)
    SUMMARY.write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps(summary, indent=1))
    figure(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
