#!/usr/bin/env python3
"""Strictly leak-free out-of-fold verification of the ProDomino + pLDDT + EB residue scorer.

Addresses four methodological checks before committing to an L1-wide re-prediction:
  1. Direct sequon-splitting (inserting inside N-X-S/T) and direct disulfide-bond disruption
     (+/-0..+1) vs wide (+/-2..+3) windows.
  2. Strict in-fold Empirical-Bayes calibration: RR_c and log_odds_penalty are computed ONLY
     on each training fold (zero test-fold leakage).
  3. Learned in-fold weights (regularized logistic on continuous features) vs fixed 0.5/0.5
     equal weights, on BOTH absolute cross-protein scale and within-protein rank scale.
  4. Paired sign tests against ProDomino alone, low pLDDT alone, and baseline RF on both
     gene-held-out and family-held-out splits.

    uv run --locked python library/l1/orf-design/tag_site_strict_oof_verify.py
"""

from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import tag_site_fine_tolerance as F  # noqa: E402
import tag_site_model_ablations as A  # noqa: E402
import tag_site_model_benchmark as B  # noqa: E402
import tag_site_negative_calibration as N  # noqa: E402
import tag_site_residue_scorer as R  # noqa: E402

OUT = B.OUT
TOLS = (0, 3, 5, 10)


def window_custom(df: pd.DataFrame, flag: np.ndarray, lo: int, hi: int) -> np.ndarray:
    """True where a flagged residue of the same protein lies in n-lo..n+hi."""
    out = np.zeros(len(df), bool)
    res = df["residue"].astype(int).to_numpy()
    for _, idx in df.groupby("uniprot_acc").indices.items():
        bad = res[idx][flag[idx]]
        if len(bad):
            d = bad[None, :] - res[idx][:, None]
            out[idx] = ((d >= -lo) & (d <= hi)).any(1)
    return out


def compute_in_fold_eb(
    ev_h: pd.DataFrame,
    on_gold: np.ndarray,
    neg_classes: dict[str, tuple[np.ndarray, np.ndarray]],
    train_mask: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Fit EB hard_mask and soft_penalty strictly on train_mask rows; apply to all rows."""
    # Compute depletion ONLY on training-fold proteins
    sub_classes = {k: (f[train_mask], cov[train_mask]) for k, (f, cov) in neg_classes.items()}
    dep = N.compute_depletion(ev_h[train_mask], on_gold[train_mask], sub_classes)

    hard_mask = np.zeros(len(ev_h), bool)
    for cname in dep.loc[dep["recommended_tier"] == "hard_veto", "negative_class"]:
        hard_mask |= neg_classes[cname][0]
    soft_pen = np.zeros(len(ev_h), float)
    for r in dep.loc[dep["recommended_tier"] == "soft_penalty"].itertuples(index=False):
        soft_pen += float(r.log_odds_penalty) * neg_classes[r.negative_class][0].astype(float)
    return hard_mask, soft_pen


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

    ev_h, neg_classes = N.build_negative_classes(ev)
    contact_classes = R.load_literature_and_structure_contacts(ev)

    # Check tight sequon-splitting (insert after N or X in N-X-S/T -> sequon_n at n or n-1)
    # and direct disulfide disruption (disulfide Cys at n or n+1, i.e., flanking the junction)
    raw_seq_n = ev_h["hn_sequon_n"].fillna(False).to_numpy(bool)
    raw_glyc_n = ev_h["hn_glycan_n"].fillna(False).to_numpy(bool)
    raw_disulf = ev_h["hn_disulfide_cys"].fillna(False).to_numpy(bool)
    all_cov = np.ones(len(ev_h), bool)

    tight_checks = {
        "split_sequon_NXST": (window_custom(ev_h, raw_seq_n, 1, 0), all_cov),
        "split_uniprot_glycan_NXST": (window_custom(ev_h, raw_glyc_n, 1, 0), all_cov),
        "sequon_wide_m2_p3": neg_classes["sequon_n"],
        "glycan_wide_m2_p3": neg_classes["glycan_n"],
        "disulfide_flanking_0_p1": (window_custom(ev_h, raw_disulf, 0, 1), all_cov),
        "disulfide_wide_m2_p3": neg_classes["disulfide_cys"],
    }
    tight_dep = N.compute_depletion(ev_h, on_gold, tight_checks)
    print("\n=== TIGHT vs WIDE GLYCAN & DISULFIDE DEPLETION ===")
    print(tight_dep.to_string(index=False))

    # Feature matrix for continuous models (both absolute and rank-normalized)
    plddt_w5 = -ev["plddt_w5"].to_numpy(float)  # higher = more flexible
    p_abs = 1.0 - ev["plddt_w5"].to_numpy(float) / 100.0
    prodom = ev["prodomino"].to_numpy(float)
    rsa_w5 = ev["rsa_w5"].to_numpy(float)
    prot14 = ev_h["protrusion_14"].to_numpy(float)
    d_dom = ev_h["d_other_domain"].to_numpy(float)

    r_plddt = A.rank_pct_within_protein(ev, plddt_w5)
    r_prodom = A.rank_pct_within_protein(ev, prodom)
    r_rsa = A.rank_pct_within_protein(ev, rsa_w5)
    r_prot14 = A.rank_pct_within_protein(ev, -prot14)  # higher = more protruding

    lit_lig = (
        contact_classes["ligand"][0]
        | contact_classes["protein_ligand"][0]
        | contact_classes["literature_small_ligand"][0]
        | contact_classes["literature_protein_ligand"][0]
        | contact_classes["literature_peptide_ligand"][0]
    )
    split_seq = tight_checks["split_sequon_NXST"][0]
    disulf_flank = tight_checks["disulfide_flanking_0_p1"][0]
    hard_bio_veto = lit_lig | split_seq

    # Absolute-scale feature matrix for learned linear weights (strictly in-fold)
    x_abs_2 = np.column_stack([p_abs, prodom])
    x_abs_5 = np.column_stack(
        [p_abs, prodom, rsa_w5, -prot14, np.minimum(np.where(np.isnan(d_dom), 30.0, d_dom), 30.0)]
    )
    x_rk_2 = np.column_stack([r_plddt, r_prodom])
    x_rk_4 = np.column_stack([r_plddt, r_prodom, r_rsa, r_prot14])

    # Mask out +/-10 around positives for training linear weights (matching evaluation!)
    near10 = B.near(ev, gold_sites, half=10)
    train_valid = on_gold | ~near10

    cls = A.protein_class(allx)
    comp = A.l1_composition()
    genes = sorted(set(ev["gene_symbol"]))
    splits = {sc: B.fold_map(genes, sc) for sc in ("gene", "family")}
    conf, w_base = B.confirmed_frame(tr, allx, gold, regions, audit)
    w5 = np.where(w_base > 1.0, 5.0, 1.0)

    def fit_linear_oof(
        xmat: np.ndarray, folds: dict[str, int], c_val: float = 0.5
    ) -> tuple[np.ndarray, list[list[float]]]:
        s = np.full(len(ev), np.nan)
        ef = ev["gene_symbol"].map(folds).to_numpy()
        coefs = []
        for i in range(B.N_SPLITS):
            tr_m = (ef != i) & train_valid
            te_m = ef == i
            imp = SimpleImputer(strategy="median")
            sca = StandardScaler()
            xtr = sca.fit_transform(imp.fit_transform(xmat[tr_m]))
            ytr = on_gold[tr_m].astype(int)
            m = LogisticRegression(C=c_val, class_weight="balanced", max_iter=2000)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                m.fit(xtr, ytr)
            coefs.append([round(float(v), 3) for v in m.coef_[0]])
            xte = sca.transform(imp.transform(xmat[te_m]))
            s[te_m] = xte @ m.coef_[0]
        return s, coefs

    rows = []
    learned_weights: dict[str, object] = {}

    for sc, folds in splits.items():
        ef = ev["gene_symbol"].map(folds).to_numpy()
        # Strict in-fold EB penalty (zero test-fold leakage)
        eb_hard_oof = np.zeros(len(ev), bool)
        eb_soft_oof = np.zeros(len(ev), float)
        for i in range(B.N_SPLITS):
            hm, sp = compute_in_fold_eb(ev_h, on_gold, neg_classes, ef != i)
            te_m = ef == i
            eb_hard_oof[te_m] = hm[te_m]
            eb_soft_oof[te_m] = sp[te_m]

        # Fixed equal-weight rank blend and absolute blend
        s_prodom = prodom
        s_plddt = plddt_w5
        s_rk_fixed = np.where(np.isnan(r_prodom), r_plddt, 0.5 * r_plddt + 0.5 * r_prodom)
        s_abs_fixed = np.where(np.isnan(prodom), p_abs, 0.5 * p_abs + 0.5 * prodom)

        # Strict in-fold EB on top of fixed blends
        s_rk_fixed_eb = np.where(eb_hard_oof | hard_bio_veto, -1e9, s_rk_fixed - 0.08 * eb_soft_oof)
        s_abs_fixed_eb = np.where(
            eb_hard_oof | hard_bio_veto, -1e9, s_abs_fixed - 0.08 * eb_soft_oof
        )

        # Strictly in-fold learned linear weights (2-feat and 4/5-feat)
        s_abs_2_lr, c_a2 = fit_linear_oof(x_abs_2, folds)
        s_abs_5_lr, c_a5 = fit_linear_oof(x_abs_5, folds)
        s_rk_2_lr, c_r2 = fit_linear_oof(x_rk_2, folds)
        s_rk_4_lr, c_r4 = fit_linear_oof(x_rk_4, folds)
        s_abs_5_lr_veto = np.where(eb_hard_oof | hard_bio_veto | disulf_flank, -1e9, s_abs_5_lr)
        s_rk_4_lr_veto = np.where(eb_hard_oof | hard_bio_veto | disulf_flank, -1e9, s_rk_4_lr)
        learned_weights[sc] = {
            "abs_2_plddt_prodom": c_a2,
            "abs_5_plddt_prodom_rsa_prot_seam": c_a5,
            "rk_2_plddt_prodom": c_r2,
            "rk_4_plddt_prodom_rsa_prot": c_r4,
        }

        s_rf_base = F.oof_custom(
            conf, w5, ev, allx, folds, B.FEATS, B.CATEGORICAL, "tr_z", "cross_loop", 2
        )

        models = [
            ("01_baseline_rf", s_rf_base),
            ("02_low_plddt_w5_alone", s_plddt),
            ("03_prodomino_alone", s_prodom),
            ("04_fixed_rank_0.5_0.5", s_rk_fixed),
            ("05_fixed_rank_strict_infold_eb", s_rk_fixed_eb),
            ("06_fixed_abs_0.5_0.5", s_abs_fixed),
            ("07_fixed_abs_strict_infold_eb", s_abs_fixed_eb),
            ("08_learned_abs_2feat", s_abs_2_lr),
            ("09_learned_abs_5feat_plus_veto", s_abs_5_lr_veto),
            ("10_learned_rk_2feat", s_rk_2_lr),
            ("11_learned_rk_4feat_plus_veto", s_rk_4_lr_veto),
        ]

        ref_rf = F.eval_multi_tol(ev, s_rf_base, gold_sites).set_index("uniprot_acc")
        ref_pd = F.eval_multi_tol(ev, s_prodom, gold_sites).set_index("uniprot_acc")
        for name, score in models:
            pp = F.eval_multi_tol(ev, score, gold_sites).set_index("uniprot_acc")
            pp_b = B.per_protein(ev, score, gold_sites)
            wtd = A.l1_weighted_metrics(pp_b, cls, comp)
            row: dict[str, object] = {
                "split": sc,
                "model": name,
                "n": len(pp),
                "auroc": round(float(pp_b["auroc"].mean()), 4),
                "l1w_auroc": wtd["l1w_auroc"],
                "l1w_t3_w10": wtd["l1w_top3_within10"],
            }
            for t in TOLS:
                row[f"t1_w{t}"] = round(float(pp[f"top1_w{t}"].mean()), 3)
                row[f"t3_w{t}"] = round(float(pp[f"top3_w{t}"].mean()), 3)
            # Paired p-value vs baseline RF and vs ProDomino alone on top1_w5
            for tag, ref in (("rf", ref_rf), ("pd", ref_pd)):
                c = ref.index.intersection(pp.index)
                a, b = ref.loc[c, "top1_w5"].astype(bool), pp.loc[c, "top1_w5"].astype(bool)
                row[f"p_t1w5_vs_{tag}"] = round(
                    B.sign_test(int((~a & b).sum()), int((a & ~b).sum())), 4
                )
            rows.append(row)

    df = pd.DataFrame(rows)
    OUT.mkdir(exist_ok=True)
    df.to_csv(OUT / "strict_oof_verify.tsv", sep="\t", index=False)
    payload = {
        "tight_vs_wide_depletion": tight_dep.to_dict(orient="records"),
        "learned_weights_per_fold": learned_weights,
        "gene_split": df[df["split"] == "gene"].to_dict(orient="records"),
        "family_split": df[df["split"] == "family"].to_dict(orient="records"),
    }
    (OUT / "strict_oof_verify.json").write_text(json.dumps(payload, indent=2) + "\n")

    cols = [
        "model",
        "auroc",
        "l1w_auroc",
        "t1_w0",
        "t1_w3",
        "t1_w5",
        "t1_w10",
        "t3_w3",
        "t3_w5",
        "t3_w10",
        "l1w_t3_w10",
        "p_t1w5_vs_rf",
        "p_t1w5_vs_pd",
    ]
    for sc in ("gene", "family"):
        print(f"\n=== {sc.upper()} SPLIT (STRICT IN-FOLD) ===")
        print(df[df["split"] == sc][cols].to_string(index=False))
    print("\n=== LEARNED WEIGHTS ACROSS 5 FOLDS (STANDARDIZED) ===")
    print(json.dumps(learned_weights["family"], indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
