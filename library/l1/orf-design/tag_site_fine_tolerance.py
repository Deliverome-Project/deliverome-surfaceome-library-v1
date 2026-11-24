#!/usr/bin/env python3
"""Fine-grained residue-level tolerance evaluation (+/-0, +/-3, +/-5, +/-10) and L1 diff.

Tests whether modeling fixes and candidate improvements before re-predicting L1 actually
improve residue-level precision (exact, +/-3, +/-5) in addition to loop-level (+/-10),
with paired sign tests against the baseline on both gene- and family-held-out splits.

Also diffs the fixed ranker against L1's currently decided loop sites to measure how many
sites move (0, 1-3, 4-10, >10 residues) and whether small shifts improve local biophysics
(pLDDT, RSA, dist_to_tm).

Outputs:
    tag-site-models/fine_tolerance.json  (committed summary record)
    tag-site-models/fine_tolerance.tsv   (local gitignored table)

    uv run --locked python library/l1/orf-design/tag_site_fine_tolerance.py
"""

from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import tag_site_model_ablations as A  # noqa: E402
import tag_site_model_benchmark as B  # noqa: E402
import tag_site_negative_calibration as N  # noqa: E402

OUT = B.OUT
TOLS = (0, 3, 5, 10)
MEMBER_FEATURES = "l1_member_extracellular.tsv"


def make_rf(feats: list[str], cats: list[str], max_depth: int = 4, min_leaf: int = 12) -> Pipeline:
    parts = [("num", SimpleImputer(strategy="median"), feats)]
    if cats:
        parts.append(("cat", OneHotEncoder(handle_unknown="ignore", min_frequency=5), cats))
    clf = RandomForestClassifier(
        n_estimators=800,
        max_depth=max_depth,
        min_samples_leaf=min_leaf,
        class_weight="balanced_subsample",
        random_state=B.RF_SEED,
        n_jobs=-1,
    )
    return Pipeline([("pre", ColumnTransformer(parts)), ("clf", clf)])


def oof_custom(
    frame: pd.DataFrame,
    w: np.ndarray,
    ev: pd.DataFrame,
    allx: pd.DataFrame,
    folds: dict[str, int],
    feats: list[str],
    cats: list[str],
    norm: str = "allx_z",
    smooth_mode: str = "contiguous",
    smooth_half: int = 2,
    max_depth: int = 4,
    min_leaf: int = 12,
) -> np.ndarray:
    s = np.full(len(ev), np.nan)
    tf = frame["gene_symbol"].map(folds).to_numpy()
    ef = ev["gene_symbol"].map(folds).to_numpy()
    if norm == "allx_z":
        f_norm = A.zscore_from_pool(frame, allx, feats)
        e_norm = A.zscore_from_pool(ev, allx, feats)
    elif norm == "tr_z":
        f_norm = B.zscore(frame, feats)
        e_norm = B.zscore(ev, feats)
    else:
        f_norm, e_norm = frame, ev
    for i in range(B.N_SPLITS):
        m = make_rf(feats, cats, max_depth=max_depth, min_leaf=min_leaf)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            m.fit(
                f_norm[tf != i],
                frame.loc[tf != i, "label"].to_numpy(),
                clf__sample_weight=w[tf != i],
            )
        held = ef == i
        if held.any():
            s[held] = m.predict_proba(e_norm[held])[:, 1]
    if smooth_mode == "cross_loop":
        return B.smooth(ev, s, smooth_half)
    if smooth_mode == "contiguous":
        return A.smooth_contiguous(ev, s, smooth_half)
    return s


def confirmed_frame_mask(
    tr: pd.DataFrame,
    allx: pd.DataFrame,
    gold: pd.DataFrame,
    regions: pd.DataFrame,
    audit: pd.DataFrame,
    mask_half: int,
) -> tuple[pd.DataFrame, np.ndarray]:
    """Same as B.confirmed_frame, with configurable negative-exclusion half-width."""
    rng = np.random.default_rng(B.SEED)
    gold_sites = gold.groupby("uniprot_acc")["residue"].apply(lambda v: np.array(v)).to_dict()
    cited = audit.dropna(subset=["insert_after_residue"])
    cited_sites = (
        cited.groupby("uniprot_acc")["insert_after_residue"]
        .apply(lambda v: np.array(v, int))
        .to_dict()
    )
    neg = tr[tr["label"] == 0]
    drop = (
        B.near(neg, gold_sites, half=mask_half)
        | B.near(neg, cited_sites, half=mask_half)
        | B.in_windows(neg, regions)
    )
    neg = neg[~drop]
    key = allx.set_index(["uniprot_acc", "residue"])
    pos_idx = [k for k in zip(gold["uniprot_acc"], gold["residue"], strict=True) if k in key.index]
    pos = key.loc[pos_idx].reset_index().assign(label=1, weight=1.0)
    new = sorted(set(pos["uniprot_acc"]) - set(tr["uniprot_acc"]))
    extra = []
    for a in new:
        pool = allx[allx["uniprot_acc"] == a]
        pool = pool[~B.near(pool, gold_sites, half=mask_half) & ~B.in_windows(pool, regions)]
        k = min(len(pool), 200 * int((pos["uniprot_acc"] == a).sum()))
        extra.append(pool.iloc[rng.choice(len(pool), size=k, replace=False)].assign(label=0))
    cols = ["uniprot_acc", "gene_symbol", "residue", *B.FEATS, *B.CATEGORICAL, "label", "weight"]
    extra_df = pd.concat(extra).assign(weight=1.0) if extra else pd.DataFrame(columns=cols)
    frame = pd.concat([pos[cols], neg[cols], extra_df[cols]], ignore_index=True)
    frame = frame.astype({"label": int, "weight": float, "residue": int})
    return frame, frame["weight"].to_numpy(float)


def eval_multi_tol(
    ev: pd.DataFrame,
    scores: np.ndarray,
    gold_sites: dict[str, np.ndarray],
) -> pd.DataFrame:
    """Per-protein hits at exact (+/-0), +/-3, +/-5, +/-10 for top-1 and top-3 (sep=11)."""
    rows = []
    tmp = ev[["uniprot_acc", "residue"]].assign(s=scores)
    for acc, idx in tmp.groupby("uniprot_acc").indices.items():
        s = tmp["s"].to_numpy()[idx]
        ok = ~np.isnan(s)
        if not ok.any():
            continue
        idx, s = idx[ok], s[ok]
        res = tmp["residue"].to_numpy()[idx]
        picks = B.picks(res, s)
        picked = res[picks]
        g = gold_sites[acc]
        d1 = int(np.abs(g - picked[0]).min())
        d3 = int(min(np.abs(g - r).min() for r in picked))
        row: dict[str, object] = {"uniprot_acc": acc, "d1": d1, "d3": d3}
        for t in TOLS:
            row[f"top1_w{t}"] = d1 <= t
            row[f"top3_w{t}"] = d3 <= t
        rows.append(row)
    return pd.DataFrame(rows)


def l1_site_diff(
    conf: pd.DataFrame,
    w5: np.ndarray,
    allx: pd.DataFrame,
) -> dict:
    """Compare full-fit allx_z contiguous RF vs L1's currently decided loop sites."""
    mem_path = B.DATA / MEMBER_FEATURES
    if not mem_path.is_file():
        return {"skipped": f"{MEMBER_FEATURES} not cached locally"}
    dec = pd.read_csv(B.MODELS / "l1_internal_decisions.tsv", sep="\t")
    dec = dec[dec["decision"].str.startswith("internal_loop")]
    tgt = pd.read_csv(B.MODELS / "l1_adopted_targets.tsv", sep="\t")
    acc_map = dict(zip(tgt["sequence_sha256"], tgt["uniprot_acc"], strict=False))
    dec = dec.assign(uniprot_acc=dec["sequence_sha256"].map(acc_map)).dropna(subset=["uniprot_acc"])

    mem = pd.read_csv(mem_path, sep="\t", low_memory=False)
    for c in B.CATEGORICAL:
        mem[c] = mem[c].fillna("").astype(str)
    mem["residue"] = mem["residue"].astype(int)
    pool = pd.concat(
        [mem, allx[~allx["uniprot_acc"].isin(set(mem["uniprot_acc"]))]], ignore_index=True
    )
    pool = pool[pool["uniprot_acc"].isin(set(dec["uniprot_acc"]))].reset_index(drop=True)

    m_fixed = make_rf(B.FEATS, B.CATEGORICAL)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m_fixed.fit(
            A.zscore_from_pool(conf, allx, B.FEATS),
            conf["label"].to_numpy(),
            clf__sample_weight=w5,
        )
    s_raw = m_fixed.predict_proba(A.zscore_from_pool(pool, pool, B.FEATS))[:, 1]
    s_sm = A.smooth_contiguous(pool, s_raw, 2)
    pool = pool.assign(s_fixed=s_sm)

    rows = []
    for r in dec.itertuples(index=False):
        g = pool[pool["uniprot_acc"] == r.uniprot_acc]
        if g.empty:
            continue
        site_old = int(r.insert_after_residue)
        old_row = g[g["residue"] == site_old]
        if old_row.empty:
            continue
        res = g["residue"].to_numpy()
        sc = g["s_fixed"].to_numpy()
        picks = B.picks(res, sc)
        top1_new = int(res[picks[0]])
        top3_new = [int(res[i]) for i in picks]
        new_row = g[g["residue"] == top1_new].iloc[0]
        old_r = old_row.iloc[0]
        dist_top1 = abs(top1_new - site_old)
        dist_top3 = min(abs(x - site_old) for x in top3_new)
        rows.append(
            {
                "gene": r.gene,
                "uniprot_acc": r.uniprot_acc,
                "old_site": site_old,
                "new_top1": top1_new,
                "dist_to_new_top1": dist_top1,
                "dist_to_new_top3": dist_top3,
                "delta_plddt": float(new_row["plddt"] - old_r["plddt"]),
                "delta_rsa": float(new_row["rsa"] - old_r["rsa"]),
                "delta_dist_to_tm": float(new_row["dist_to_tm"] - old_r["dist_to_tm"])
                if pd.notna(new_row["dist_to_tm"]) and pd.notna(old_r["dist_to_tm"])
                else np.nan,
            }
        )
    df = pd.DataFrame(rows)
    small = df[(df["dist_to_new_top1"] >= 1) & (df["dist_to_new_top1"] <= 10)]
    large = df[df["dist_to_new_top1"] > 10]
    return {
        "n_scored": len(df),
        "exact_same_top1": int((df["dist_to_new_top1"] == 0).sum()),
        "shift_1_to_3": int(((df["dist_to_new_top1"] >= 1) & (df["dist_to_new_top1"] <= 3)).sum()),
        "shift_4_to_10": int(
            ((df["dist_to_new_top1"] >= 4) & (df["dist_to_new_top1"] <= 10)).sum()
        ),
        "shift_gt_10": len(large),
        "old_outside_new_top3_w10": int((df["dist_to_new_top3"] > 10).sum()),
        "small_shift_1_10_medians": {
            "n": len(small),
            "median_delta_plddt": round(float(small["delta_plddt"].median()), 3)
            if len(small)
            else None,
            "median_delta_rsa": round(float(small["delta_rsa"].median()), 4)
            if len(small)
            else None,
            "median_delta_dist_to_tm": round(float(small["delta_dist_to_tm"].median()), 2)
            if len(small)
            else None,
            "frac_lower_plddt": round(float((small["delta_plddt"] < 0).mean()), 3)
            if len(small)
            else None,
            "frac_higher_rsa": round(float((small["delta_rsa"] > 0).mean()), 3)
            if len(small)
            else None,
        },
        "large_shift_gt10_medians": {
            "n": len(large),
            "median_delta_plddt": round(float(large["delta_plddt"].median()), 3)
            if len(large)
            else None,
            "median_delta_rsa": round(float(large["delta_rsa"].median()), 4)
            if len(large)
            else None,
            "frac_lower_plddt": round(float((large["delta_plddt"] < 0).mean()), 3)
            if len(large)
            else None,
            "frac_higher_rsa": round(float((large["delta_rsa"] > 0).mean()), 3)
            if len(large)
            else None,
        },
    }


def main() -> int:
    tr, allx = B.inputs()
    struct = pd.read_csv(B.MODELS / "data" / "structure_global_features.tsv", sep="\t")
    allx = allx.merge(
        struct[["uniprot_acc", "residue", "protrusion_14", "d_other_domain"]],
        on=["uniprot_acc", "residue"],
        how="left",
    )
    gold, regions, audit = B.truth()
    gold_sites = gold.groupby("uniprot_acc")["residue"].apply(lambda v: np.array(v)).to_dict()

    ev = allx[allx["uniprot_acc"].isin(gold_sites)].reset_index(drop=True)
    have = set(zip(ev["uniprot_acc"], ev["residue"], strict=True))
    gold_in = gold[
        [(a, r) in have for a, r in zip(gold["uniprot_acc"], gold["residue"], strict=True)]
    ]
    gold_sites = gold_in.groupby("uniprot_acc")["residue"].apply(lambda v: np.array(v)).to_dict()
    ev = ev[ev["uniprot_acc"].isin(gold_sites)].reset_index(drop=True)

    conf, w_base = B.confirmed_frame(tr, allx, gold, regions, audit)
    conf = conf.merge(
        struct[["uniprot_acc", "residue", "protrusion_14", "d_other_domain"]],
        on=["uniprot_acc", "residue"],
        how="left",
    )
    w5 = np.where(w_base > 1.0, 5.0, 1.0)
    conf_m3, w_m3 = confirmed_frame_mask(tr, allx, gold, regions, audit, mask_half=3)
    w5_m3 = np.where(w_m3 > 1.0, 5.0, 1.0)

    # Precompute EB penalties & full veto
    key = ev.set_index(["uniprot_acc", "residue"]).index
    on_gold = key.isin(list(zip(gold_in["uniprot_acc"], gold_in["residue"], strict=True)))
    ev_h, classes = N.build_negative_classes(ev.drop(columns=["protrusion_14", "d_other_domain"]))
    dep = N.compute_depletion(ev_h, on_gold, classes)
    hard_mask = np.zeros(len(ev), bool)
    for cname in dep.loc[dep["recommended_tier"] == "hard_veto", "negative_class"]:
        hard_mask |= classes[cname][0]
    soft_pen = np.zeros(len(ev), float)
    for r in dep.loc[dep["recommended_tier"] == "soft_penalty"].itertuples(index=False):
        soft_pen += float(r.log_odds_penalty) * classes[r.negative_class][0].astype(float)
    _bind, full_veto, _flags = A.annotate_veto_flags(ev)

    genes = sorted(set(conf["gene_symbol"]) | set(ev["gene_symbol"]))
    splits = {sc: B.fold_map(genes, sc) for sc in ("gene", "family")}

    plddt_raw = -ev["plddt"].to_numpy(float)
    plddt_w5 = -ev["plddt_w5"].to_numpy(float)
    prodomino = ev["prodomino"].to_numpy(float)
    r_plddt = A.rank_pct_within_protein(ev, A.smooth_contiguous(ev, plddt_w5, 2))
    r_prodom = A.rank_pct_within_protein(ev, prodomino)
    blend_plddt_prodom = np.where(np.isnan(r_prodom), r_plddt, 0.5 * r_plddt + 0.5 * r_prodom)

    feats_ps = [*B.FEATS, "protrusion_14", "d_other_domain"]
    summary_rows = []

    for sc, folds in splits.items():
        s_base = oof_custom(
            conf, w5, ev, allx, folds, B.FEATS, B.CATEGORICAL, "tr_z", "cross_loop", 2
        )
        s_sm2 = oof_custom(
            conf, w5, ev, allx, folds, B.FEATS, B.CATEGORICAL, "allx_z", "contiguous", 2
        )
        s_sm1 = oof_custom(
            conf, w5, ev, allx, folds, B.FEATS, B.CATEGORICAL, "allx_z", "contiguous", 1
        )
        s_sm0 = oof_custom(conf, w5, ev, allx, folds, B.FEATS, B.CATEGORICAL, "allx_z", "none", 0)
        s_d6 = oof_custom(
            conf, w5, ev, allx, folds, B.FEATS, B.CATEGORICAL, "allx_z", "contiguous", 1, 6, 6
        )
        s_m3 = oof_custom(
            conf_m3, w5_m3, ev, allx, folds, B.FEATS, B.CATEGORICAL, "allx_z", "contiguous", 1
        )
        s_ps = oof_custom(
            conf, w5, ev, allx, folds, feats_ps, B.CATEGORICAL, "allx_z", "contiguous", 1
        )
        s_eb_sm2 = np.where(hard_mask, -1e9, N.to_logit(s_sm2) - soft_pen)
        s_eb_sm1 = np.where(hard_mask, -1e9, N.to_logit(s_sm1) - soft_pen)
        s_eb_sm0 = np.where(hard_mask, -1e9, N.to_logit(s_sm0) - soft_pen)
        s_veto = np.where(full_veto & ~np.isnan(s_sm2), -1e9, s_sm2)
        r_eb = A.rank_pct_within_protein(ev, s_eb_sm1)
        s_eb_prodom = np.where(np.isnan(r_prodom), r_eb, 0.5 * r_eb + 0.5 * r_prodom)

        models = [
            ("01_baseline_rf_sm2", s_base),
            ("02_allxz_contig_sm2", s_sm2),
            ("03_allxz_contig_sm1", s_sm1),
            ("04_allxz_nosmooth", s_sm0),
            ("05_allxz_rf_d6_l6_sm1", s_d6),
            ("06_allxz_mask3_sm1", s_m3),
            ("07_allxz_protrusion_seam_feat", s_ps),
            ("08_allxz_sm2_eb_calibrated", s_eb_sm2),
            ("09_allxz_sm1_eb_calibrated", s_eb_sm1),
            ("10_allxz_sm0_eb_calibrated", s_eb_sm0),
            ("11_allxz_sm2_full_veto", s_veto),
            ("12_low_plddt_raw_residue", plddt_raw),
            ("13_low_plddt_w5", plddt_w5),
            ("14_plddt_prodomino_blend", blend_plddt_prodom),
            ("15_allxz_eb_prodomino_blend", s_eb_prodom),
        ]

        ref_pp = eval_multi_tol(ev, s_base, gold_sites).set_index("uniprot_acc")
        for name, score in models:
            pp = eval_multi_tol(ev, score, gold_sites).set_index("uniprot_acc")
            row: dict[str, object] = {"split": sc, "model": name, "n": len(pp)}
            for t in TOLS:
                row[f"t1_w{t}"] = round(float(pp[f"top1_w{t}"].mean()), 3)
                row[f"t3_w{t}"] = round(float(pp[f"top3_w{t}"].mean()), 3)
            # Paired sign test vs baseline on t1_w3, t1_w5, t3_w5
            common = ref_pp.index.intersection(pp.index)
            for mcol in ("top1_w3", "top1_w5", "top3_w5"):
                a = ref_pp.loc[common, mcol].astype(bool)
                b = pp.loc[common, mcol].astype(bool)
                wins, losses = int((~a & b).sum()), int((a & ~b).sum())
                row[f"p_{mcol}"] = round(B.sign_test(wins, losses), 4)
            summary_rows.append(row)

    df = pd.DataFrame(summary_rows)
    l1_diff = l1_site_diff(conf, w5, allx)

    OUT.mkdir(exist_ok=True)
    df.to_csv(OUT / "fine_tolerance.tsv", sep="\t", index=False)
    payload = {
        "tolerances": list(TOLS),
        "gene_split": df[df["split"] == "gene"].to_dict(orient="records"),
        "family_split": df[df["split"] == "family"].to_dict(orient="records"),
        "l1_site_diff": l1_diff,
    }
    (OUT / "fine_tolerance.json").write_text(json.dumps(payload, indent=2) + "\n")

    cols = [
        "model",
        "t1_w0",
        "t1_w3",
        "t1_w5",
        "t1_w10",
        "t3_w0",
        "t3_w3",
        "t3_w5",
        "t3_w10",
        "p_top1_w5",
        "p_top3_w5",
    ]
    for sc in ("gene", "family"):
        print(f"\n=== {sc.upper()} SPLIT (t1/t3 at +/-0, +/-3, +/-5, +/-10) ===")
        print(df[df["split"] == sc][cols].to_string(index=False))
    print("\n=== L1 BUILT LOOP SITE DIFF (OLD DECIDED vs FIXED ALLX_Z RF) ===")
    print(json.dumps(l1_diff, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
