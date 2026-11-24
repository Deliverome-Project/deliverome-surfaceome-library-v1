#!/usr/bin/env python3
"""Production candidate L1 internal-site designer: ProDomino + pLDDT + EB penalties + vetoes.

Combines:
  - Within-protein residue ranker: 0.50 * rank_pct(-plddt_w5) + 0.50 * rank_pct(prodomino)
    minus Empirical-Bayes soft penalties (protrusion_14 crowding, domain_seam_6A/8A,
    wide disulfide +/-2..+3, buried core, literature partner).
  - Hard vetoes: membrane-side gate (junction_mask == 'E'), structure + transferred +
    literature ligands (n-2..n+3), direct flanking disulfide Cys (0..+1), split N-glycan
    sequon (N|X|S/T).
  - Cross-protein absolute scorer (for the 1 / 2 / 3 construct + snorkel order rule):
    5-feature standardized logistic regression on [1 - plddt_w5/100, prodomino, rsa_w5,
    -protrusion_14, min(d_other_domain, 30)], calibrated out-of-fold on the 143 confirmed
    hosting regions.

Outputs:
    tag-site-models/l1_v8_prodomino_design_summary.json  (committed summary record)
    tag-site-models/l1_v8_prodomino_internal_designs.tsv  (per-CDS design & comparison vs v5; published to S3)

    uv run --locked python library/l1/orf-design/design_l1_internal_v8_prodomino.py
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
REPO = HERE.parents[2]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(REPO / "library/l1/tag-site-benchmarking"))

import membrane_mask as MM  # noqa: E402
import tag_site_model_ablations as A  # noqa: E402
import tag_site_model_benchmark as B  # noqa: E402
import tag_site_negative_calibration as N  # noqa: E402
import tag_site_strict_oof_verify as S  # noqa: E402

DATA = HERE / "tag-site-models" / "data"
OUT = HERE / "tag-site-models"
MIN_REGION_LEN = 4
MIN_SEPARATION = 25
CLEAR_MARGIN = 0.90
Q75_PROT = 380.0
Q90_PROT = 459.0


def fit_absolute_model(
    allx: pd.DataFrame, gold: pd.DataFrame
) -> tuple[SimpleImputer, StandardScaler, LogisticRegression, float, float]:
    """Fit 5-feature absolute scorer on confirmed sites and compute OOF region-ceiling median/q25."""
    gold_sites = gold.groupby("uniprot_acc")["residue"].apply(lambda v: np.array(v)).to_dict()
    ev = allx[allx["uniprot_acc"].isin(gold_sites)].reset_index(drop=True)
    have = set(zip(ev["uniprot_acc"], ev["residue"], strict=True))
    gold_in = gold[
        [(a, r) in have for a, r in zip(gold["uniprot_acc"], gold["residue"], strict=True)]
    ]
    gold_sites = gold_in.groupby("uniprot_acc")["residue"].apply(lambda v: np.array(v)).to_dict()
    ev = ev[ev["uniprot_acc"].isin(gold_sites)].reset_index(drop=True)

    struct = pd.read_csv(B.MODELS / "data" / "structure_global_features.tsv", sep="\t")
    ev = ev.merge(
        struct[["uniprot_acc", "residue", "protrusion_14", "d_other_domain"]],
        on=["uniprot_acc", "residue"],
        how="left",
    )
    key = ev.set_index(["uniprot_acc", "residue"]).index
    on_gold = key.isin(list(zip(gold_in["uniprot_acc"], gold_in["residue"], strict=True)))
    near10 = B.near(ev, gold_sites, half=10)
    train_valid = on_gold | ~near10

    def feat_mat(df: pd.DataFrame) -> np.ndarray:
        p_abs = 1.0 - df["plddt_w5"].to_numpy(float) / 100.0
        prodom = df["prodomino"].to_numpy(float)
        rsa = df["rsa_w5"].to_numpy(float)
        prot = -df["protrusion_14"].to_numpy(float)
        ddom = df["d_other_domain"].to_numpy(float)
        ddom = np.minimum(np.where(np.isnan(ddom), 30.0, ddom), 30.0)
        return np.column_stack([p_abs, prodom, rsa, prot, ddom])

    x_ev = feat_mat(ev)
    genes = sorted(set(ev["gene_symbol"]))
    folds = B.fold_map(genes, "family")
    ef = ev["gene_symbol"].map(folds).to_numpy()
    oof_p = np.full(len(ev), np.nan)
    for i in range(B.N_SPLITS):
        tr_m = (ef != i) & train_valid
        te_m = ef == i
        imp = SimpleImputer(strategy="median")
        sca = StandardScaler()
        xtr = sca.fit_transform(imp.fit_transform(x_ev[tr_m]))
        m = LogisticRegression(C=0.5, class_weight="balanced", max_iter=2000)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            m.fit(xtr, on_gold[tr_m].astype(int))
        oof_p[te_m] = m.predict_proba(sca.transform(imp.transform(x_ev[te_m])))[:, 1]

    ev_oof = ev.assign(oof=oof_p)
    rkey = ["uniprot_acc", "region_index"]
    reg_max = ev_oof.groupby(rkey, dropna=False)["oof"].max().rename("region_ceiling").reset_index()
    pos_regs = ev_oof.loc[on_gold, rkey].drop_duplicates()
    hosting = pos_regs.merge(reg_max, on=rkey, how="inner")["region_ceiling"].dropna().to_numpy()
    med, q25 = float(np.median(hosting)), float(np.quantile(hosting, 0.25))

    # Full fit on all confirmed sites
    imp_f = SimpleImputer(strategy="median")
    sca_f = StandardScaler()
    x_tr = sca_f.fit_transform(imp_f.fit_transform(x_ev[train_valid]))
    m_f = LogisticRegression(C=0.5, class_weight="balanced", max_iter=2000)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m_f.fit(x_tr, on_gold[train_valid].astype(int))
    return imp_f, sca_f, m_f, med, q25


def pick_arms(
    g: pd.DataFrame,
    pub_med: float,
    pub_q25: float,
) -> dict:
    """Apply L1's 1/2/3 construct + snorkel order rule on eligible, non-vetoed junctions."""
    elig = g[g["eligible"] & ~g["hard_veto"]].copy()
    if elig.empty:
        return {"status": "no_eligible_site", "n_arms": 0, "picks": []}

    # Prefer regions of length >= MIN_REGION_LEN if any exist
    r_lens = elig.groupby("region_index")["residue"].transform("count")
    if (r_lens >= MIN_REGION_LEN).any():
        elig = elig[r_lens >= MIN_REGION_LEN].copy()

    # Region ceilings from absolute scorer
    reg_ceil = elig.groupby("region_index")["s_abs"].max().sort_values(ascending=False)
    best_ceil = float(reg_ceil.iloc[0])
    runner_ceil = float(reg_ceil.iloc[1]) if len(reg_ceil) >= 2 else 0.0

    snorkel_companion = best_ceil < pub_q25
    if snorkel_companion:
        n_target = 2  # + 1 snorkel companion = 3
    elif best_ceil < pub_med:
        n_target = 3
    elif runner_ceil >= CLEAR_MARGIN * best_ceil:
        n_target = 2
    else:
        n_target = 1

    # Greedy picks separated by >= MIN_SEPARATION (or in different regions)
    order = elig.sort_values("s_rank", ascending=False)
    picks = []
    for r in order.itertuples(index=False):
        if all(
            (r.region_index != p["region_index"])
            or (abs(int(r.residue) - p["residue"]) >= MIN_SEPARATION)
            for p in picks
        ):
            picks.append(
                {
                    "residue": int(r.residue),
                    "region_index": int(r.region_index) if pd.notna(r.region_index) else -1,
                    "s_rank": round(float(r.s_rank), 4),
                    "s_abs": round(float(r.s_abs), 4),
                    "plddt": float(r.plddt) if pd.notna(r.plddt) else np.nan,
                    "prodomino": float(r.prodomino) if pd.notna(r.prodomino) else np.nan,
                    "protrusion_14": float(r.protrusion_14)
                    if pd.notna(r.protrusion_14)
                    else np.nan,
                    "d_other_domain": float(r.d_other_domain)
                    if pd.notna(r.d_other_domain)
                    else np.nan,
                }
            )
            if len(picks) == n_target:
                break

    return {
        "status": "designed",
        "best_ceiling": round(best_ceil, 4),
        "n_loop_arms": len(picks),
        "snorkel_companion": snorkel_companion,
        "picks": picks,
    }


def _load_backfilled_pool(allx: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    pd_path = (
        OUT / "l1_prodomino_scores.tsv"
        if (OUT / "l1_prodomino_scores.tsv").exists()
        else DATA / "l1_prodomino_scores.tsv"
    )
    sg_path = (
        OUT / "l1_structure_global_features.tsv"
        if (OUT / "l1_structure_global_features.tsv").exists()
        else DATA / "l1_structure_global_features.tsv"
    )
    pd_new = pd.read_csv(pd_path, sep="\t")
    sg_new = pd.read_csv(sg_path, sep="\t")

    p56_new = pd_new[pd_new["uniprot_acc"] == "P56750"][["residue", "prodomino"]].rename(
        columns={"prodomino": "prodom_new"}
    )
    p56_old = allx[allx["uniprot_acc"] == "P56750"][["residue", "prodomino"]].rename(
        columns={"prodomino": "prodom_old"}
    )
    par = p56_old.merge(p56_new, on="residue", how="inner").dropna()
    parity_max_diff = (
        float(np.abs(par["prodom_new"] - par["prodom_old"]).max()) if len(par) else None
    )
    parity_corr = float(par["prodom_new"].corr(par["prodom_old"])) if len(par) else None
    print(
        f"ProDomino parity check (CLDN17 P56750, n={len(par)}): "
        f"max_abs_diff={parity_max_diff:.5f}, corr={parity_corr:.6f}"
    )

    dec = pd.read_csv(B.MODELS / "l1_internal_decisions.tsv", sep="\t")
    dec = dec[dec["decision"].str.startswith("internal_loop")].copy()
    tgt = pd.read_csv(B.MODELS / "l1_adopted_targets.tsv", sep="\t")
    sha_to_acc = dict(zip(tgt["sequence_sha256"], tgt["uniprot_acc"], strict=False))
    acc_to_sha = dict(zip(tgt["uniprot_acc"], tgt["sequence_sha256"], strict=False))
    dec["uniprot_acc"] = dec["sequence_sha256"].map(sha_to_acc)
    dec = dec.dropna(subset=["uniprot_acc"]).reset_index(drop=True)

    mem = pd.read_csv(DATA / "l1_member_extracellular.tsv", sep="\t", low_memory=False)
    pool = pd.concat(
        [mem, allx[~allx["uniprot_acc"].isin(set(mem["uniprot_acc"]))]], ignore_index=True
    )
    pool = pool[pool["uniprot_acc"].isin(set(dec["uniprot_acc"]))].copy().reset_index(drop=True)
    pool["residue"] = pool["residue"].astype(int)
    pool["sequence_sha256"] = pool["uniprot_acc"].map(acc_to_sha)

    pd_map = (
        pd_new[["sequence_sha256", "residue", "prodomino"]]
        .drop_duplicates(["sequence_sha256", "residue"])
        .rename(columns={"prodomino": "prodomino_backfill"})
    )
    pool = pool.merge(pd_map, on=["sequence_sha256", "residue"], how="left")
    pool["prodomino"] = pool["prodomino"].fillna(pool["prodomino_backfill"])

    sg_map = sg_new[
        ["sequence_sha256", "residue", "protrusion_14", "d_other_domain"]
    ].drop_duplicates(["sequence_sha256", "residue"])
    pool = pool.merge(sg_map, on=["sequence_sha256", "residue"], how="left")

    mm = MM.load()
    jmask_by_sha = dict(
        zip(
            mm.loc[mm["usable"] == 1, "sequence_sha256"],
            mm.loc[mm["usable"] == 1, "junction_mask"],
            strict=False,
        )
    )
    pool["eligible"] = [
        True
        if (jm := jmask_by_sha.get(sha)) is None
        else (1 <= res <= len(jm) and jm[res - 1] == MM.ELIGIBLE)
        for sha, res in zip(pool["sequence_sha256"], pool["residue"], strict=True)
    ]
    parity = {
        "n_residues": len(par),
        "max_abs_diff": round(parity_max_diff, 5) if parity_max_diff is not None else None,
        "pearson_r": round(parity_corr, 6) if parity_corr is not None else None,
    }
    return dec, pool, parity


def _annotate_and_score_pool(pool: pd.DataFrame, imp_f, sca_f, m_f) -> pd.DataFrame:
    c_df = pd.read_csv(
        B.MODELS / "data" / "contacts_by_protein.tsv",
        sep="\t",
        usecols=["sequence_sha256", "residue", "class"],
        low_memory=False,
    ).dropna(subset=["sequence_sha256", "residue"])
    c_df["residue"] = c_df["residue"].astype(int)
    lig_classes = {
        "ligand",
        "protein_ligand",
        "literature_small_ligand",
        "literature_protein_ligand",
        "literature_peptide_ligand",
    }
    lit_part_classes = {"literature_partner", "literature_homomer"}
    lig_keys = set(
        zip(
            c_df.loc[c_df["class"].isin(lig_classes), "sequence_sha256"],
            c_df.loc[c_df["class"].isin(lig_classes), "residue"],
            strict=True,
        )
    )
    lp_keys = set(
        zip(
            c_df.loc[c_df["class"].isin(lit_part_classes), "sequence_sha256"],
            c_df.loc[c_df["class"].isin(lit_part_classes), "residue"],
            strict=True,
        )
    )
    pkeys = list(zip(pool["sequence_sha256"], pool["residue"], strict=True))
    lig_win = N.window_mask(pool, np.array([k in lig_keys for k in pkeys], bool))
    lp_win = N.window_mask(pool, np.array([k in lp_keys for k in pkeys], bool))

    u = pd.read_csv(B.MODELS / "uniprot_site_features.tsv", sep="\t")
    dis_keys = set()
    for r in u[u["type"] == "DISULFID"].itertuples(index=False):
        for col in ("begin", "end"):
            v = getattr(r, col, None)
            if pd.notna(v):
                try:
                    dis_keys.add((r.uniprot_acc, int(v)))
                except ValueError, TypeError:
                    pass
    ukeys = list(zip(pool["uniprot_acc"], pool["residue"], strict=True))
    raw_dis = np.array([k in dis_keys for k in ukeys], bool)
    dis_flank = S.window_custom(pool, raw_dis, 0, 1)
    dis_wide = N.window_mask(pool, raw_dis)

    raw_seq = pool["glyc_sequon_here"].fillna(0).to_numpy(float) > 0.5
    split_seq = S.window_custom(pool, raw_seq, 1, 0)
    buried = (pool["rsa"].fillna(1.0) < 0.1) & (pool["plddt"].fillna(0.0) >= 80.0)
    buried_win = N.window_mask(pool, buried.to_numpy(bool))

    pool["hard_veto"] = lig_win | dis_flank | split_seq
    pool["lig_win"] = lig_win
    pool["dis_flank"] = dis_flank
    pool["split_seq"] = split_seq
    pool["crowded_q75"] = (pool["protrusion_14"] >= Q75_PROT).fillna(False)
    pool["crowded_q90"] = (pool["protrusion_14"] >= Q90_PROT).fillna(False)
    pool["seam_6A"] = (pool["d_other_domain"] <= 6.0).fillna(False)
    pool["seam_8A"] = (pool["d_other_domain"] <= 8.0).fillna(False)

    r_plddt = A.rank_pct_within_protein(pool, -pool["plddt_w5"].to_numpy(float))
    r_prodom = A.rank_pct_within_protein(pool, pool["prodomino"].to_numpy(float))
    base_rank = np.where(
        np.isnan(r_prodom),
        r_plddt,
        np.where(np.isnan(r_plddt), r_prodom, 0.5 * r_plddt + 0.5 * r_prodom),
    )
    soft_pen = (
        1.7569 * pool["crowded_q75"].to_numpy(float)
        + 1.8071 * pool["crowded_q90"].to_numpy(float)
        + 0.5497 * pool["seam_6A"].to_numpy(float)
        + 0.4774 * pool["seam_8A"].to_numpy(float)
        + 0.7190 * dis_wide.astype(float)
        + 1.3774 * buried_win.astype(float)
    )
    pool["s_rank"] = base_rank - 0.08 * soft_pen - 0.05 * lp_win.astype(float)

    x_abs = np.column_stack(
        [
            1.0 - pool["plddt_w5"].to_numpy(float) / 100.0,
            pool["prodomino"].to_numpy(float),
            pool["rsa_w5"].to_numpy(float),
            -pool["protrusion_14"].to_numpy(float),
            np.minimum(np.where(pool["d_other_domain"].isna(), 30.0, pool["d_other_domain"]), 30.0),
        ]
    )
    pool["s_abs"] = m_f.predict_proba(sca_f.transform(imp_f.transform(x_abs)))[:, 1]
    return pool


def _design_members(
    dec: pd.DataFrame, pool: pd.DataFrame, pub_med: float, pub_q25: float
) -> pd.DataFrame:
    design_rows = []
    for r in dec.itertuples(index=False):
        g = pool[pool["uniprot_acc"] == r.uniprot_acc]
        if g.empty:
            continue
        res_v8 = pick_arms(g, pub_med, pub_q25)
        old_site = int(r.insert_after_residue)
        old_row = g[g["residue"] == old_site]
        old_r = old_row.iloc[0] if not old_row.empty else None

        if res_v8["status"] != "designed" or not res_v8["picks"]:
            design_rows.append(
                {
                    "gene": r.gene,
                    "uniprot_acc": r.uniprot_acc,
                    "sequence_sha256": r.sequence_sha256,
                    "status": res_v8["status"],
                    "old_v5_site": old_site,
                    "v8_top1": np.nan,
                    "n_loop_arms": 0,
                    "snorkel_companion": False,
                }
            )
            continue

        p1 = res_v8["picks"][0]
        p_all = [p["residue"] for p in res_v8["picks"]]
        new_r = g[g["residue"] == p1["residue"]].iloc[0]
        design_rows.append(
            {
                "gene": r.gene,
                "uniprot_acc": r.uniprot_acc,
                "sequence_sha256": r.sequence_sha256,
                "status": "designed",
                "old_v5_site": old_site,
                "v8_top1": p1["residue"],
                "v8_picks": ";".join(str(x) for x in p_all),
                "n_loop_arms": res_v8["n_loop_arms"],
                "snorkel_companion": res_v8["snorkel_companion"],
                "best_ceiling": res_v8["best_ceiling"],
                "dist_to_v8_top1": abs(p1["residue"] - old_site),
                "dist_to_v8_top3": min(abs(x - old_site) for x in p_all),
                "delta_plddt": float(new_r["plddt"] - old_r["plddt"])
                if old_r is not None
                else np.nan,
                "delta_prodomino": float(new_r["prodomino"] - old_r["prodomino"])
                if old_r is not None
                else np.nan,
                "delta_protrusion_14": float(new_r["protrusion_14"] - old_r["protrusion_14"])
                if old_r is not None
                else np.nan,
                "old_crowded_q75": bool(old_r["crowded_q75"]) if old_r is not None else False,
                "new_crowded_q75": bool(new_r["crowded_q75"]),
                "old_seam_8A": bool(old_r["seam_8A"]) if old_r is not None else False,
                "new_seam_8A": bool(new_r["seam_8A"]),
                "old_hard_veto": bool(old_r["hard_veto"]) if old_r is not None else False,
                "new_hard_veto": bool(new_r["hard_veto"]),
            }
        )
    return pd.DataFrame(design_rows)


def _build_summary(pool: pd.DataFrame, ddf: pd.DataFrame, parity: dict) -> dict:
    ok = ddf[ddf["status"] == "designed"]
    mv = ok[ok["dist_to_v8_top1"] > 0]

    def pct(col: str) -> float:
        return round(
            float(pool.groupby("uniprot_acc")[col].apply(lambda s: s.notna().any()).mean() * 100), 1
        )

    return {
        "prodomino_parity_cldn17": parity,
        "l1_feature_coverage_after_backfill": {
            "n_proteins": int(pool["uniprot_acc"].nunique()),
            "pct_with_plddt": pct("plddt_w5"),
            "pct_with_prodomino": pct("prodomino"),
            "pct_with_protrusion_14": pct("protrusion_14"),
            "pct_with_d_other_domain": pct("d_other_domain"),
        },
        "v8_prodomino_designs": {
            "n_scored": len(ddf),
            "n_designed": len(ok),
            "n_no_eligible_site": int((ddf["status"] != "designed").sum()),
            "loop_arm_counts": ok["n_loop_arms"].value_counts().sort_index().to_dict(),
            "n_snorkel_companions": int(ok["snorkel_companion"].sum()),
            "total_loop_constructs": int(ok["n_loop_arms"].sum()),
        },
        "vs_v5_comparison": {
            "exact_same_top1": int((ok["dist_to_v8_top1"] == 0).sum()),
            "shift_1_to_3": int(
                ((ok["dist_to_v8_top1"] >= 1) & (ok["dist_to_v8_top1"] <= 3)).sum()
            ),
            "shift_4_to_10": int(
                ((ok["dist_to_v8_top1"] >= 4) & (ok["dist_to_v8_top1"] <= 10)).sum()
            ),
            "shift_gt_10": int((ok["dist_to_v8_top1"] > 10).sum()),
            "old_outside_v8_picks_w10": int((ok["dist_to_v8_top3"] > 10).sum()),
            "moved_median_delta_plddt": round(float(mv["delta_plddt"].median()), 2),
            "moved_frac_lower_plddt": round(float((mv["delta_plddt"] < 0).mean()), 3),
            "moved_median_delta_prodomino": round(float(mv["delta_prodomino"].median()), 4),
            "moved_frac_higher_prodomino": round(float((mv["delta_prodomino"] > 0).mean()), 3),
            "moved_median_delta_protrusion_14": round(float(mv["delta_protrusion_14"].median()), 1),
            "moved_frac_less_crowded": round(float((mv["delta_protrusion_14"] < 0).mean()), 3),
            "old_top1_crowded_q75": int(ok["old_crowded_q75"].sum()),
            "new_top1_crowded_q75": int(ok["new_crowded_q75"].sum()),
            "old_top1_seam_8A": int(ok["old_seam_8A"].sum()),
            "new_top1_seam_8A": int(ok["new_seam_8A"].sum()),
            "old_top1_hard_veto": int(ok["old_hard_veto"].sum()),
            "new_top1_hard_veto": int(ok["new_hard_veto"].sum()),
        },
    }


def main() -> int:
    _tr, allx = B.inputs()
    gold, _regions, _audit = B.truth()
    imp_f, sca_f, m_f, pub_med, pub_q25 = fit_absolute_model(allx, gold)
    print(f"Absolute model OOF hosting-region ceilings: median={pub_med:.4f}, q25={pub_q25:.4f}")

    dec, pool, parity = _load_backfilled_pool(allx)
    pool = _annotate_and_score_pool(pool, imp_f, sca_f, m_f)
    ddf = _design_members(dec, pool, pub_med, pub_q25)
    ddf.to_csv(OUT / "l1_v8_prodomino_internal_designs.tsv", sep="\t", index=False)

    summary = _build_summary(pool, ddf, parity)
    (OUT / "l1_v8_prodomino_design_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
