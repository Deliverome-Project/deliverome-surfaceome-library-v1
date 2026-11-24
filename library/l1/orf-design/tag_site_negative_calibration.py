#!/usr/bin/env python3
"""Principled PU weighting for candidate negative classes (interfaces, seams, protrusion).

In positive-unlabeled (PU) data we only observe ~149 confirmed working sites, so a positive-only
benchmark can measure the COST of avoiding a negative class (confirmed sites pushed out of top-3)
but never the BENEFIT (silent failures avoided on unlabeled residues).

This script calibrates picking-time penalties for candidate negative classes in two ways:
  1. Depletion ratios (RR_c) on the 149 confirmed sites vs extracellular background of the same
     covered proteins, yielding parameter-free log-odds penalties lambda_c = max(0, -log(RR_c))
     and Fisher exact p-values.
  2. Soft-penalty Pareto sweeps on out-of-fold log-odds scores of the fixed ranker
     (allx_z + contiguous smooth), measuring risk-exposure reduction vs confirmed-site recall
     (gene- and family-held-out, unweighted and L1-class-mix weighted), plus a combined
     Empirical-Bayes calibrated policy and multi-construct hedging rule.

Outputs:
    tag-site-models/negative_calibration.json  (committed summary record)
    tag-site-models/negative_depletion.tsv     (local gitignored table)
    tag-site-models/negative_pareto.tsv        (local gitignored table)

    uv run --locked python library/l1/orf-design/tag_site_negative_calibration.py
"""

from __future__ import annotations

import json
import math
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import tag_site_model_ablations as A  # noqa: E402
import tag_site_model_benchmark as B  # noqa: E402

sys.path.insert(0, str(B.MODELS))
import ranker_hard_negative_layer as L  # noqa: E402
import ranker_pu_segment as R  # noqa: E402

OUT = B.OUT
VETO_LO, VETO_HI = 2, 3
LAMBDAS = (0.0, 0.25, 0.5, 1.0, 2.0, 4.0, 1e6)


def window_mask(df: pd.DataFrame, flag: np.ndarray) -> np.ndarray:
    """True where a flagged residue of the same protein lies in n-VETO_LO..n+VETO_HI."""
    out = np.zeros(len(df), bool)
    res = df["residue"].astype(int).to_numpy()
    for _, idx in df.groupby("uniprot_acc").indices.items():
        bad = res[idx][flag[idx]]
        if len(bad):
            d = bad[None, :] - res[idx][:, None]
            out[idx] = ((d >= -VETO_LO) & (d <= VETO_HI)).any(1)
    return out


def build_negative_classes(
    ev: pd.DataFrame,
) -> tuple[pd.DataFrame, dict[str, tuple[np.ndarray, np.ndarray]]]:
    """Return annotated ev and dict[class_name -> (flag_window, covered_mask)].

    `covered_mask` marks rows belonging to proteins where the annotation source exists
    (e.g., proteins with an oligomer model for d_oligomer_partner, or multi-domain proteins
    for d_other_domain), so depletion ratios are computed on the covered denominator.
    """
    L.INTERFACES = "all"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        h = L.constraints(ev, R.load_screens())

    struct = pd.read_csv(B.MODELS / "data" / "structure_global_features.tsv", sep="\t")
    h = h.merge(
        struct[
            [
                "uniprot_acc",
                "residue",
                "d_other_domain",
                "d_oligomer_partner",
                "protrusion_14",
                "protrusion_14_oligomer",
            ]
        ],
        on=["uniprot_acc", "residue"],
        how="left",
    )

    all_cov = np.ones(len(h), bool)
    dom_cov = h["d_other_domain"].notna().to_numpy()
    olig_cov = h["d_oligomer_partner"].notna().to_numpy()
    prot_cov = h["protrusion_14"].notna().to_numpy()
    prot_olig_cov = h["protrusion_14_oligomer"].notna().to_numpy()

    q75_prot = float(h["protrusion_14"].dropna().quantile(0.75))
    q90_prot = float(h["protrusion_14"].dropna().quantile(0.90))
    q75_olig = float(h["protrusion_14_oligomer"].dropna().quantile(0.75))

    raw_defs: dict[str, tuple[np.ndarray, np.ndarray, bool]] = {
        # name: (residue_flag, covered_mask, apply_window)
        "ligand_contact": (h["hn_ligand_contact"].fillna(False).to_numpy(bool), all_cov, True),
        "partner_interface": (
            h["hn_partner_interface"].fillna(False).to_numpy(bool),
            all_cov,
            True,
        ),
        "homomer_interface_5A": (
            h["hn_homomer_interface"].fillna(False).to_numpy(bool),
            olig_cov,
            True,
        ),
        "homomer_near_8A": (
            (h["d_oligomer_partner"] <= 8.0).fillna(False).to_numpy(bool),
            olig_cov,
            True,
        ),
        "domain_seam_6A": (
            (h["d_other_domain"] <= 6.0).fillna(False).to_numpy(bool),
            dom_cov,
            True,
        ),
        "domain_seam_8A": (
            (h["d_other_domain"] <= 8.0).fillna(False).to_numpy(bool),
            dom_cov,
            True,
        ),
        "crowded_protrusion_q75": (
            (h["protrusion_14"] >= q75_prot).fillna(False).to_numpy(bool),
            prot_cov,
            False,
        ),
        "crowded_protrusion_q90": (
            (h["protrusion_14"] >= q90_prot).fillna(False).to_numpy(bool),
            prot_cov,
            False,
        ),
        "crowded_oligomer_q75": (
            (h["protrusion_14_oligomer"] >= q75_olig).fillna(False).to_numpy(bool),
            prot_olig_cov,
            False,
        ),
        "buried_core": (h["hn_buried"].fillna(False).to_numpy(bool), all_cov, True),
        "disulfide_cys": (h["hn_disulfide_cys"].fillna(False).to_numpy(bool), all_cov, True),
        "glycan_n": (h["hn_glycan_n"].fillna(False).to_numpy(bool), all_cov, True),
        "sequon_n": (h["hn_sequon_n"].fillna(False).to_numpy(bool), all_cov, True),
        "screen_intolerant": (
            h["hn_screen_intolerant"].fillna(False).to_numpy(bool),
            all_cov,
            True,
        ),
    }

    classes = {
        name: (window_mask(h, flag) if use_win else flag, cov)
        for name, (flag, cov, use_win) in raw_defs.items()
    }
    return h, classes


def compute_depletion(
    ev: pd.DataFrame,
    on_gold: np.ndarray,
    classes: dict[str, tuple[np.ndarray, np.ndarray]],
) -> pd.DataFrame:
    """Relative risk RR_c = P(in class | positive) / P(in class | background) on covered rows."""
    rows = []
    for name, (flag, cov) in classes.items():
        pos = on_gold & cov
        bg = (~on_gold) & cov
        n_pos = int(pos.sum())
        n_bg = int(bg.sum())
        k_pos = int((flag & pos).sum())
        k_bg = int((flag & bg).sum())
        pos_rate = k_pos / n_pos if n_pos else np.nan
        bg_rate = k_bg / n_bg if n_bg else np.nan
        rr = (pos_rate / bg_rate) if (n_pos and k_bg > 0) else np.nan
        # Laplace-smoothed RR (+0.5 / +1 pseudocount on positives) for finite log-odds penalty
        pos_smooth = (k_pos + 0.5) / (n_pos + 1.0) if n_pos else np.nan
        bg_smooth = (k_bg + 0.5) / (n_bg + 1.0) if n_bg else np.nan
        rr_smooth = pos_smooth / bg_smooth if (n_pos and n_bg) else np.nan
        _, pval = fisher_exact([[k_pos, n_pos - k_pos], [k_bg, n_bg - k_bg]])
        log_odds = max(0.0, -math.log(rr_smooth)) if rr_smooth and rr_smooth > 0 else 0.0
        if rr_smooth < 0.15:
            tier = "hard_veto"
        elif rr_smooth < 0.70:
            tier = "soft_penalty"
        else:
            tier = "do_not_penalize"
        rows.append(
            {
                "negative_class": name,
                "n_pos_covered": n_pos,
                "k_pos_in_class": k_pos,
                "pos_rate": round(pos_rate, 4),
                "bg_rate": round(bg_rate, 4),
                "rr": round(rr, 4) if not np.isnan(rr) else np.nan,
                "rr_smooth": round(rr_smooth, 4),
                "fisher_p": round(float(pval), 5),
                "log_odds_penalty": round(log_odds, 4),
                "recommended_tier": tier,
            }
        )
    return pd.DataFrame(rows).sort_values("rr_smooth").reset_index(drop=True)


def to_logit(p: np.ndarray, eps: float = 1e-4) -> np.ndarray:
    clipped = np.clip(p, eps, 1.0 - eps)
    return np.log(clipped / (1.0 - clipped))


def hedged_picks(
    res: np.ndarray,
    scores: np.ndarray,
    soft_flag: np.ndarray,
    k: int = 3,
    sep: int = 11,
) -> list[int]:
    """Pick 1 is best overall; if Pick 1 is soft-flagged, force Pick 2 to be clean."""
    order = np.argsort(-scores, kind="stable")
    chosen: list[int] = []
    for i in order:
        if np.isnan(scores[i]):
            continue
        if len(chosen) == 1 and soft_flag[chosen[0]] and soft_flag[i]:
            continue
        if all(abs(int(res[i]) - int(res[j])) >= sep for j in chosen):
            chosen.append(int(i))
            if len(chosen) == k:
                break
    # Fallback if no clean site remained for Pick 2
    if len(chosen) < k:
        for i in order:
            if np.isnan(scores[i]) or int(i) in chosen:
                continue
            if all(abs(int(res[i]) - int(res[j])) >= sep for j in chosen):
                chosen.append(int(i))
                if len(chosen) == k:
                    break
    return chosen


def eval_hedged(
    ev: pd.DataFrame,
    scores: np.ndarray,
    soft_flag: np.ndarray,
    gold_sites: dict[str, np.ndarray],
    cls: dict[str, str],
    comp: dict[str, int],
) -> dict:
    """Evaluate hedged picking: top-1/top-3 recall + fraction with ANY or ALL picks exposed."""
    rows = []
    any_exp, all_exp = [], []
    tmp = ev[["uniprot_acc", "residue"]].assign(s=scores)
    for acc, idx in tmp.groupby("uniprot_acc").indices.items():
        s = tmp["s"].to_numpy()[idx]
        ok = ~np.isnan(s)
        if not ok.any():
            continue
        idx, s = idx[ok], s[ok]
        res = tmp["residue"].to_numpy()[idx]
        sf = soft_flag[idx]
        picks = hedged_picks(res, s, sf)
        picked_res = res[picks]
        g = gold_sites[acc]
        t1 = bool(np.abs(g - picked_res[0]).min() <= B.FLANK)
        t3 = bool(any(np.abs(g - r).min() <= B.FLANK for r in picked_res))
        rows.append(
            {
                "uniprot_acc": acc,
                "top1_within10": t1,
                "top3_within10": t3,
                "auroc": np.nan,
            }
        )
        any_exp.append(bool(sf[picks].any()))
        all_exp.append(bool(sf[picks].all()))
    pp = pd.DataFrame(rows)
    wtd = A.l1_weighted_metrics(pp, cls, comp)
    return {
        "top1_within10": round(float(pp["top1_within10"].mean()), 4),
        "top3_within10": round(float(pp["top3_within10"].mean()), 4),
        "l1w_top1_within10": wtd["l1w_top1_within10"],
        "l1w_top3_within10": wtd["l1w_top3_within10"],
        "any_pick_soft_exposed": round(float(np.mean(any_exp)), 4),
        "all_picks_soft_exposed": round(float(np.mean(all_exp)), 4),
    }


def main() -> int:
    tr, allx = B.inputs()
    gold, regions, audit = B.truth()
    gold_sites = gold.groupby("uniprot_acc")["residue"].apply(lambda v: np.array(v)).to_dict()

    ev = allx[allx["uniprot_acc"].isin(gold_sites)].reset_index(drop=True)
    have = set(zip(ev["uniprot_acc"], ev["residue"], strict=True))
    gold_in = gold[
        [(a, r) in have for a, r in zip(gold["uniprot_acc"], gold["residue"], strict=True)]
    ]
    gold_sites = gold_in.groupby("uniprot_acc")["residue"].apply(lambda v: np.array(v)).to_dict()
    ev = ev[ev["uniprot_acc"].isin(gold_sites)].reset_index(drop=True)

    key = ev.set_index(["uniprot_acc", "residue"]).index
    on_gold = key.isin(list(zip(gold_in["uniprot_acc"], gold_in["residue"], strict=True)))

    ev_h, classes = build_negative_classes(ev)
    dep = compute_depletion(ev_h, on_gold, classes)
    print("\n=== PU DEPLETION RATIOS ON 149 CONFIRMED SITES ===")
    print(dep.to_string(index=False))

    # Fit out-of-fold fixed ranker (04_fix_allx_zscore_contig) on gene and family splits
    conf, w_base = B.confirmed_frame(tr, allx, gold, regions, audit)
    w5 = np.where(w_base > 1.0, 5.0, 1.0)
    cls = A.protein_class(allx)
    comp = A.l1_composition()
    genes = sorted(set(conf["gene_symbol"]) | set(ev["gene_symbol"]))
    splits = {sc: B.fold_map(genes, sc) for sc in ("gene", "family")}

    oof_logits = {}
    for sc, folds in splits.items():
        p = A.oof_rf(conf, w5, ev, allx, folds, B.FEATS, B.CATEGORICAL, "allx_z", "contiguous")
        oof_logits[sc] = to_logit(p)

    # Individual Pareto sweeps for candidate soft/uncertain negative classes
    sweep_classes = [
        "partner_interface",
        "homomer_interface_5A",
        "domain_seam_6A",
        "domain_seam_8A",
        "crowded_protrusion_q75",
        "crowded_protrusion_q90",
    ]
    pareto_rows = []
    for sc in ("gene", "family"):
        base = oof_logits[sc]
        for cname in sweep_classes:
            flag, _cov = classes[cname]
            for lam in LAMBDAS:
                pen = base - lam * flag.astype(float)
                pp = B.per_protein(ev, pen, gold_sites)
                summ = B.summarise(pp)
                wtd = A.l1_weighted_metrics(pp, cls, comp)
                exp = A.exposure_rates(ev, pen, {cname: flag})[f"exp_{cname}"]
                pareto_rows.append(
                    {
                        "split": sc,
                        "negative_class": cname,
                        "lambda": lam,
                        "top1_within10": summ["top1_within10"],
                        "top3_within10": summ["top3_within10"],
                        "l1w_top1_within10": wtd["l1w_top1_within10"],
                        "l1w_top3_within10": wtd["l1w_top3_within10"],
                        "exposure": exp,
                    }
                )

    pareto_df = pd.DataFrame(pareto_rows)

    # Combined Empirical-Bayes calibrated policy:
    #   hard veto on tier == "hard_veto"
    #   additive log-odds penalty (-log(RR_smooth)) on tier == "soft_penalty"
    hard_names = dep.loc[dep["recommended_tier"] == "hard_veto", "negative_class"].tolist()
    soft_rows = dep.loc[dep["recommended_tier"] == "soft_penalty"]
    hard_mask = np.zeros(len(ev), bool)
    for n in hard_names:
        hard_mask |= classes[n][0]
    soft_penalty = np.zeros(len(ev), float)
    soft_any = np.zeros(len(ev), bool)
    for r in soft_rows.itertuples(index=False):
        f = classes[r.negative_class][0]
        soft_penalty += float(r.log_odds_penalty) * f.astype(float)
        soft_any |= f

    # Compare policies: baseline, hard-only, EB calibrated (hard + soft log-odds), and hedged
    eval_flags = {
        "hard_any": hard_mask,
        "soft_any": soft_any,
        "binding": classes["ligand_contact"][0]
        | classes["partner_interface"][0]
        | classes["homomer_interface_5A"][0],
        "seam_8A": classes["domain_seam_8A"][0],
        "crowded_q75": classes["crowded_protrusion_q75"][0],
    }
    combined_summary = []
    hedging_summary = {}
    for sc in ("gene", "family"):
        base = oof_logits[sc]
        policies = {
            "unpenalized_allxz_rf": base,
            "hard_veto_only": np.where(hard_mask, -1e9, base),
            "eb_calibrated_hard_plus_soft": np.where(hard_mask, -1e9, base - soft_penalty),
            "eb_calibrated_soft_x0.5": np.where(hard_mask, -1e9, base - 0.5 * soft_penalty),
        }
        for pname, score in policies.items():
            pp = B.per_protein(ev, score, gold_sites)
            summ = B.summarise(pp)
            wtd = A.l1_weighted_metrics(pp, cls, comp)
            exp = A.exposure_rates(ev, score, eval_flags)
            combined_summary.append({"split": sc, "policy": pname, **summ, **wtd, **exp})
        hedging_summary[sc] = eval_hedged(
            ev, np.where(hard_mask, -1e9, base), soft_any, gold_sites, cls, comp
        )

    comb_df = pd.DataFrame(combined_summary)
    print("\n=== COMBINED EMPIRICAL-BAYES POLICIES ===")
    print(
        comb_df[
            [
                "split",
                "policy",
                "top1_within10",
                "top3_within10",
                "l1w_top1_within10",
                "l1w_top3_within10",
                "exp_hard_any",
                "exp_binding",
                "exp_seam_8A",
                "exp_crowded_q75",
            ]
        ].to_string(index=False)
    )
    print("\n=== MULTI-CONSTRUCT HEDGING (HARD VETO + FORCE CLEAN PICK 2 IF PICK 1 IS SOFT) ===")
    print(json.dumps(hedging_summary, indent=2))

    OUT.mkdir(exist_ok=True)
    dep.to_csv(OUT / "negative_depletion.tsv", sep="\t", index=False)
    pareto_df.to_csv(OUT / "negative_pareto.tsv", sep="\t", index=False)
    payload = {
        "n_confirmed_sites": len(gold_in),
        "n_proteins": int(ev["uniprot_acc"].nunique()),
        "hard_veto_classes": hard_names,
        "soft_penalty_classes": {
            r.negative_class: float(r.log_odds_penalty) for r in soft_rows.itertuples(index=False)
        },
        "depletion": dep.to_dict(orient="records"),
        "combined_policies": comb_df.to_dict(orient="records"),
        "hedging": hedging_summary,
        "pareto": pareto_df.to_dict(orient="records"),
    }
    (OUT / "negative_calibration.json").write_text(json.dumps(payload, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
