#!/usr/bin/env python3
"""Diagnostic ablations and model-improvement tests for the internal tag-site ranker.

Tests specific modeling flaws and improvement hypotheses on the corpus-audit confirmed sites
(149 sites on 108 proteins), evaluated out-of-fold under both gene-grouped and HGNC
family-grouped 5-fold CV, reported unweighted and reweighted to L1's internal class mix:

  Bug fixes tested:
    1. Contiguous-only smoothing (never averaging across TM-helix boundaries).
    2. Pre-subsampling within-protein z-scoring (computed on all extracellular residues
       rather than on the hard-negative-enriched training frame).
    3. Dropping the collinear AlphaFold disorder block (disorder, disorder_w5, disorder_w15).

  Modeling & feature improvements tested:
    4. Adding zero-missingness multipass loop/TM geometry features (loop_len, dist_to_tm,
       region_len, rel_pos_in_region).
    5. Reweighting training proteins to L1's internal class mix (95.7% multipass vs 69.4%
       in the training corpus).
    6. Regularized L2 logistic regression (avoiding balanced_subsample x sample_weight=5
       leaf distortion).
    7. Masked within-protein pairwise logistic ranking (excluding +/-10 residues around
       positives from negative pairs, fixing learning_to_rank.py's unmasked pairs).
    8. Decoupling intrinsic tolerance from hard-negative exclusions (applying deterministic
       binding/structural vetoes on top of the tolerance scorer).
    9. Blending low pLDDT with ProDomino (ESM-2 based) to estimate the PLM ceiling.

Outputs:
    tag-site-models/ablations_summary.json  (committed summary record)
    tag-site-models/ablations_summary.tsv   (local gitignored table)
    tag-site-models/ablations_by_class.tsv  (local gitignored table)

    uv run --locked python library/l1/orf-design/tag_site_model_ablations.py
"""

from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import tag_site_model_benchmark as B  # noqa: E402

sys.path.insert(0, str(B.MODELS))
import ranker_hard_negative_layer as L  # noqa: E402
import ranker_pu_segment as R  # noqa: E402

OUT = B.OUT
DISORDER = {"disorder", "disorder_w5", "disorder_w15"}
FEATS_NO_DISORDER = [f for f in B.FEATS if f not in DISORDER]
GEOM_FEATS = ["loop_len", "dist_to_tm", "region_len", "rel_pos_in_region"]
CLASSES = ("single-pass", "multipass, 2-6 TM", "multipass, 7+ TM", "no TM helix")
BINDING = ["hn_ligand_contact", "hn_partner_interface", "hn_homomer_interface"]
STRUCTURAL = ["hn_disulfide_cys", "hn_glycan_n", "hn_sequon_n", "hn_buried", "hn_screen_intolerant"]
VETO_LO, VETO_HI = 2, 3


def add_region_geometry(df: pd.DataFrame) -> pd.DataFrame:
    """Derive region_len and relative distance from loop ends (0 at edge, 0.5 at midpoint)."""
    out = df.copy()
    g = out.groupby(["uniprot_acc", "region_index"], dropna=False)["residue"]
    rmin = g.transform("min")
    rmax = g.transform("max")
    rlen = (rmax - rmin + 1).astype(float)
    edge_dist = np.minimum(out["residue"] - rmin, rmax - out["residue"]).astype(float)
    out["region_len"] = rlen
    out["rel_pos_in_region"] = np.where(rlen > 1, edge_dist / (rlen - 1.0), 0.0)
    return out


def smooth_contiguous(df: pd.DataFrame, s: np.ndarray, half: int = B.SMOOTH_HALF) -> np.ndarray:
    """Average scores over +/-half residues ONLY within contiguous residue runs (no cross-TM)."""
    if half <= 0:
        return s.copy()
    out = s.copy()
    tmp = df[["uniprot_acc", "residue"]].assign(s=s)
    for _, g in tmp.groupby("uniprot_acc"):
        g = g.sort_values("residue")
        res = g["residue"].to_numpy()
        vals = g["s"].to_numpy(dtype=float)
        idx = g.index.to_numpy()
        # Split into contiguous runs where residue increments by 1
        breaks = np.where(np.diff(res) != 1)[0] + 1
        starts = np.concatenate([[0], breaks])
        ends = np.concatenate([breaks, [len(res)]])
        for lo, hi in zip(starts, ends, strict=True):
            seg = vals[lo:hi]
            if len(seg) <= 1 or np.isnan(seg).all():
                out[idx[lo:hi]] = seg
                continue
            # NaN-aware moving average within contiguous segment
            sm = np.empty_like(seg)
            for k in range(len(seg)):
                w = seg[max(0, k - half) : min(len(seg), k + half + 1)]
                sm[k] = np.nanmean(w) if not np.isnan(w).all() else np.nan
            out[idx[lo:hi]] = sm
    return out


def zscore_from_pool(df: pd.DataFrame, pool: pd.DataFrame, feats: list[str]) -> pd.DataFrame:
    """Z-score df using per-protein mean/std computed on the full extracellular pool."""
    stats = pool.groupby("uniprot_acc")[feats].agg(["mean", "std"])
    out = df.copy()
    acc = out["uniprot_acc"]
    for f in feats:
        mu = acc.map(stats[(f, "mean")])
        sd = acc.map(stats[(f, "std")]).replace(0, np.nan)
        out[f] = (out[f] - mu) / sd
    return out


def protein_class(allx: pd.DataFrame) -> dict[str, str]:
    tm = allx.groupby("uniprot_acc")["tm_helix_count"].max()
    out = {}
    for a, n in tm.items():
        n = 0 if n != n else int(n)
        out[a] = (
            CLASSES[3] if n == 0 else CLASSES[0] if n == 1 else CLASSES[1] if n <= 6 else CLASSES[2]
        )
    return out


def l1_composition() -> dict[str, int]:
    tp = pd.read_csv(HERE / "l1_tag_paths.tsv", sep="\t", low_memory=False)
    n = tp.loc[tp["placement"] == "internal", "n_tm"].fillna(0).astype(int)
    lab = np.select([n == 0, n == 1, n <= 6], [CLASSES[3], CLASSES[0], CLASSES[1]], CLASSES[2])
    return pd.Series(lab).value_counts().reindex(CLASSES, fill_value=0).astype(int).to_dict()


def l1_weighted_metrics(pp: pd.DataFrame, cls: dict[str, str], comp: dict[str, int]) -> dict:
    d = pp.assign(protein_class=pp["uniprot_acc"].map(cls))
    have = [c for c in CLASSES if (d["protein_class"] == c).any()]
    w = np.array([comp[c] for c in have], dtype=float)
    w /= w.sum()
    res = {}
    for metric in ("top1_within10", "top3_within10", "auroc"):
        vals = [d.loc[d["protein_class"] == c, metric].astype(float).mean() for c in have]
        res[f"l1w_{metric}"] = round(float(np.dot(w, vals)), 4)
    return res


def annotate_veto_flags(rows: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, dict[str, np.ndarray]]:
    L.INTERFACES = "all"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        h = L.constraints(rows, R.load_screens())
    for c in BINDING + STRUCTURAL:
        h[c] = h[c].fillna(False).astype(bool)

    def win(flag: np.ndarray) -> np.ndarray:
        out = np.zeros(len(h), bool)
        res = h["residue"].astype(int).to_numpy()
        for _, idx in h.groupby("uniprot_acc").indices.items():
            bad = res[idx][flag[idx]]
            if len(bad):
                d = bad[None, :] - res[idx][:, None]
                out[idx] = ((d >= -VETO_LO) & (d <= VETO_HI)).any(1)
        return out

    binding = win(h[BINDING].any(axis=1).to_numpy())
    structural = win(h[STRUCTURAL].any(axis=1).to_numpy())
    full = binding | structural
    flags = {
        "binding": binding,
        "ligand": win(h["hn_ligand_contact"].to_numpy()),
        "partner": win(h["hn_partner_interface"].to_numpy()),
        "homomer": win(h["hn_homomer_interface"].to_numpy()),
        "structural": structural,
    }
    return binding, full, flags


def exposure_rates(rows: pd.DataFrame, scores: np.ndarray, flags: dict[str, np.ndarray]) -> dict:
    tmp = rows[["uniprot_acc", "residue"]].assign(s=scores)
    hits = {k: [] for k in flags}
    for _, idx in tmp.groupby("uniprot_acc").indices.items():
        s = tmp["s"].to_numpy()[idx]
        ok = ~np.isnan(s)
        if not ok.any():
            continue
        idx, s = idx[ok], s[ok]
        top = idx[B.picks(tmp["residue"].to_numpy()[idx], s)]
        for k, f in flags.items():
            hits[k].append(bool(f[top].any()))
    return {f"exp_{k}": round(float(np.mean(v)), 4) for k, v in hits.items()}


def make_logistic(feats: list[str], cats: list[str], c_val: float = 0.2) -> Pipeline:
    num_pipe = Pipeline([("imp", SimpleImputer(strategy="median")), ("sca", StandardScaler())])
    parts = [("num", num_pipe, feats)]
    if cats:
        parts.append(("cat", OneHotEncoder(handle_unknown="ignore", min_frequency=5), cats))
    clf = LogisticRegression(
        C=c_val, class_weight="balanced", max_iter=3000, random_state=B.RF_SEED
    )
    return Pipeline([("pre", ColumnTransformer(parts)), ("clf", clf)])


def oof_rf(
    frame: pd.DataFrame,
    w: np.ndarray,
    ev: pd.DataFrame,
    allx: pd.DataFrame,
    folds: dict[str, int],
    feats: list[str],
    cats: list[str],
    norm: str,
    smooth_mode: str,
) -> np.ndarray:
    s = np.full(len(ev), np.nan)
    tf = frame["gene_symbol"].map(folds).to_numpy()
    ef = ev["gene_symbol"].map(folds).to_numpy()
    if norm == "allx_z":
        f_norm = zscore_from_pool(frame, allx, feats)
        e_norm = zscore_from_pool(ev, allx, feats)
    elif norm == "tr_z":
        f_norm = B.zscore(frame, feats)
        e_norm = B.zscore(ev, feats)
    else:
        f_norm, e_norm = frame, ev
    for i in range(B.N_SPLITS):
        m = B.model(feats, cats)
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
        return B.smooth(ev, s, B.SMOOTH_HALF)
    if smooth_mode == "contiguous":
        return smooth_contiguous(ev, s, B.SMOOTH_HALF)
    return s


def oof_logistic(
    frame: pd.DataFrame,
    w: np.ndarray,
    ev: pd.DataFrame,
    allx: pd.DataFrame,
    folds: dict[str, int],
    feats: list[str],
    cats: list[str],
    norm: str,
    smooth_mode: str,
) -> np.ndarray:
    s = np.full(len(ev), np.nan)
    tf = frame["gene_symbol"].map(folds).to_numpy()
    ef = ev["gene_symbol"].map(folds).to_numpy()
    if norm == "allx_z":
        f_norm = zscore_from_pool(frame, allx, feats)
        e_norm = zscore_from_pool(ev, allx, feats)
    else:
        f_norm, e_norm = frame, ev
    for i in range(B.N_SPLITS):
        m = make_logistic(feats, cats)
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
    if smooth_mode == "contiguous":
        return smooth_contiguous(ev, s, B.SMOOTH_HALF)
    return s


def oof_pairwise_masked(
    frame: pd.DataFrame,
    ev: pd.DataFrame,
    allx: pd.DataFrame,
    folds: dict[str, int],
    feats: list[str],
    pairs_per_pos: int = 200,
) -> np.ndarray:
    """Within-protein pairwise logistic ranking on the masked confirmed frame (+/-10 excluded)."""
    s = np.full(len(ev), np.nan)
    tf = frame["gene_symbol"].map(folds).to_numpy()
    ef = ev["gene_symbol"].map(folds).to_numpy()
    f_norm = zscore_from_pool(frame, allx, feats)
    e_norm = zscore_from_pool(ev, allx, feats)
    rng = np.random.default_rng(B.SEED)

    for i in range(B.N_SPLITS):
        tr_sub = f_norm[tf != i].reset_index(drop=True)
        imp, sca = SimpleImputer(strategy="median"), StandardScaler()
        xtr = sca.fit_transform(imp.fit_transform(tr_sub[feats].to_numpy(dtype=float)))
        ytr = tr_sub["label"].to_numpy()
        acc = tr_sub["uniprot_acc"].to_numpy()

        diffs, labels = [], []
        for a in np.unique(acc[ytr == 1]):
            pi = np.where((acc == a) & (ytr == 1))[0]
            ni = np.where((acc == a) & (ytr == 0))[0]
            if not len(ni):
                continue
            for p in pi:
                take = rng.choice(ni, size=min(pairs_per_pos, len(ni)), replace=False)
                d = xtr[p] - xtr[take]
                diffs.append(d)
                labels.append(np.ones(len(take)))
                diffs.append(-d)
                labels.append(np.zeros(len(take)))
        m = LogisticRegression(C=0.2, fit_intercept=False, max_iter=3000, random_state=B.RF_SEED)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            m.fit(np.vstack(diffs), np.concatenate(labels))

        held = ef == i
        if held.any():
            xev = sca.transform(imp.transform(e_norm.loc[held, feats].to_numpy(dtype=float)))
            s[held] = xev @ m.coef_[0]
    return smooth_contiguous(ev, s, B.SMOOTH_HALF)


def rank_pct_within_protein(df: pd.DataFrame, scores: np.ndarray) -> np.ndarray:
    """Convert raw scores to within-protein percentile rank in [0, 1] (higher = better)."""
    out = np.full(len(df), np.nan)
    tmp = df[["uniprot_acc"]].assign(s=scores)
    for _, idx in tmp.groupby("uniprot_acc").indices.items():
        s = tmp["s"].to_numpy()[idx]
        ok = ~np.isnan(s)
        if ok.sum() < 2:
            continue
        r = pd.Series(s[ok]).rank(pct=True).to_numpy()
        out[idx[ok]] = r
    return out


def class_reweights(frame: pd.DataFrame, cls: dict[str, str], comp: dict[str, int]) -> np.ndarray:
    """Per-row protein-class multiplier to match L1's internal class mix."""
    c_series = frame["uniprot_acc"].map(cls)
    pos_counts = frame.loc[frame["label"] == 1, "uniprot_acc"].map(cls).value_counts()
    tot_pos = pos_counts.sum()
    tot_l1 = sum(comp.values())
    mult = {}
    for c in CLASSES:
        p_frac = pos_counts.get(c, 0) / tot_pos if tot_pos else 0.0
        l1_frac = comp.get(c, 0) / tot_l1 if tot_l1 else 0.0
        mult[c] = (l1_frac / p_frac) if p_frac > 0 else 1.0
    return c_series.map(mult).fillna(1.0).to_numpy(dtype=float)


def main() -> int:
    tr, allx = B.inputs()
    allx = add_region_geometry(allx)
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
    geom_cols = [
        "uniprot_acc",
        "residue",
        "loop_len",
        "dist_to_tm",
        "region_len",
        "rel_pos_in_region",
    ]
    conf = conf.merge(allx[geom_cols], on=["uniprot_acc", "residue"], how="left")
    w5 = np.where(w_base > 1.0, 5.0, 1.0)
    w1 = np.ones_like(w5)

    cls = protein_class(allx)
    comp = l1_composition()
    cw = class_reweights(conf, cls, comp)
    _binding, full_veto, flags = annotate_veto_flags(ev)

    genes = sorted(set(conf["gene_symbol"]) | set(ev["gene_symbol"]))
    splits = {sc: B.fold_map(genes, sc) for sc in ("gene", "family")}

    # Precompute simple baselines
    plddt_cross = B.smooth(ev, -ev["plddt_w5"].to_numpy(float), B.SMOOTH_HALF)
    plddt_contig = smooth_contiguous(ev, -ev["plddt_w5"].to_numpy(float), B.SMOOTH_HALF)
    plddt_raw = -ev["plddt_w5"].to_numpy(float)
    prodomino = ev["prodomino"].to_numpy(float)

    # Blend: rank-normalized low pLDDT + ProDomino (falling back to pLDDT where ProDomino is NaN)
    r_plddt = rank_pct_within_protein(ev, plddt_contig)
    r_prodom = rank_pct_within_protein(ev, prodomino)
    blend_plddt_prodom = np.where(np.isnan(r_prodom), r_plddt, 0.5 * r_plddt + 0.5 * r_prodom)

    summary_rows = []
    class_rows = []

    for split_name, folds in splits.items():
        # Evaluate each model variant
        s_base = oof_rf(conf, w5, ev, allx, folds, B.FEATS, B.CATEGORICAL, "tr_z", "cross_loop")
        s_contig = oof_rf(conf, w5, ev, allx, folds, B.FEATS, B.CATEGORICAL, "tr_z", "contiguous")
        s_nosm = oof_rf(conf, w5, ev, allx, folds, B.FEATS, B.CATEGORICAL, "tr_z", "none")
        s_allxz = oof_rf(conf, w5, ev, allx, folds, B.FEATS, B.CATEGORICAL, "allx_z", "contiguous")
        s_nodis = oof_rf(
            conf, w5, ev, allx, folds, FEATS_NO_DISORDER, B.CATEGORICAL, "allx_z", "contiguous"
        )
        s_glob_nodis = oof_rf(
            conf, w5, ev, allx, folds, FEATS_NO_DISORDER, B.CATEGORICAL, "raw", "none"
        )
        s_geom = oof_rf(
            conf,
            w5,
            ev,
            allx,
            folds,
            FEATS_NO_DISORDER + GEOM_FEATS,
            B.CATEGORICAL,
            "raw",
            "contiguous",
        )
        s_geom_z = oof_rf(
            conf,
            w5,
            ev,
            allx,
            folds,
            FEATS_NO_DISORDER + GEOM_FEATS,
            B.CATEGORICAL,
            "allx_z",
            "contiguous",
        )
        s_geom_l1w = oof_rf(
            conf,
            w5 * cw,
            ev,
            allx,
            folds,
            FEATS_NO_DISORDER + GEOM_FEATS,
            B.CATEGORICAL,
            "raw",
            "contiguous",
        )
        s_logit = oof_logistic(
            conf,
            w1,
            ev,
            allx,
            folds,
            FEATS_NO_DISORDER + GEOM_FEATS,
            B.CATEGORICAL,
            "allx_z",
            "contiguous",
        )
        s_pair = oof_pairwise_masked(conf, ev, allx, folds, FEATS_NO_DISORDER + GEOM_FEATS)

        candidates = [
            ("01_baseline_confirmed_rf", s_base),
            ("02_fix_contiguous_smooth", s_contig),
            ("03_fix_no_smooth", s_nosm),
            ("04_fix_allx_zscore_contig", s_allxz),
            ("05_drop_disorder_allxz_contig", s_nodis),
            ("06_global_nosmooth_drop_disorder", s_glob_nodis),
            ("07_global_nodis_plus_geom", s_geom),
            ("07b_allxz_nodis_plus_geom", s_geom_z),
            ("08_global_nodis_geom_l1weights", s_geom_l1w),
            ("09_l2_logistic_nodis_geom", s_logit),
            ("10_pairwise_masked_nodis_geom", s_pair),
            ("11_low_plddt_cross_smooth", plddt_cross),
            ("12_low_plddt_contig_smooth", plddt_contig),
            ("13_low_plddt_raw_w5", plddt_raw),
            (
                "14_low_plddt_full_veto",
                np.where(full_veto & ~np.isnan(plddt_contig), -1e9, plddt_contig),
            ),
            (
                "14b_allxz_rf_full_veto",
                np.where(full_veto & ~np.isnan(s_allxz), -1e9, s_allxz),
            ),
            ("15_nodis_rf_full_veto", np.where(full_veto & ~np.isnan(s_nodis), -1e9, s_nodis)),
            ("16_geom_rf_full_veto", np.where(full_veto & ~np.isnan(s_geom), -1e9, s_geom)),
            ("17_logit_full_veto", np.where(full_veto & ~np.isnan(s_logit), -1e9, s_logit)),
            ("18_plddt_prodomino_blend", blend_plddt_prodom),
            (
                "19_plddt_prodomino_full_veto",
                np.where(full_veto & ~np.isnan(blend_plddt_prodom), -1e9, blend_plddt_prodom),
            ),
        ]

        for name, sc in candidates:
            pp = B.per_protein(ev, sc, gold_sites)
            summ = B.summarise(pp)
            wtd = l1_weighted_metrics(pp, cls, comp)
            exp = exposure_rates(ev, sc, flags)
            summary_rows.append({"model": name, "split": split_name, **summ, **wtd, **exp})
            if split_name == "gene":
                d = pp.assign(protein_class=pp["uniprot_acc"].map(cls))
                for c in CLASSES:
                    sub = d[d["protein_class"] == c]
                    if len(sub):
                        class_rows.append({"model": name, "protein_class": c, **B.summarise(sub)})

    summ_df = pd.DataFrame(summary_rows)
    cls_df = pd.DataFrame(class_rows)
    OUT.mkdir(exist_ok=True)
    summ_df.to_csv(OUT / "ablations_summary.tsv", sep="\t", index=False)
    cls_df.to_csv(OUT / "ablations_by_class.tsv", sep="\t", index=False)

    payload = {
        "n_proteins": int(ev["uniprot_acc"].nunique()),
        "n_confirmed_sites": len(gold_in),
        "l1_class_composition": comp,
        "gene_split": summ_df[summ_df["split"] == "gene"].to_dict(orient="records"),
        "family_split": summ_df[summ_df["split"] == "family"].to_dict(orient="records"),
    }
    (OUT / "ablations_summary.json").write_text(json.dumps(payload, indent=2) + "\n")

    for sc in ("gene", "family"):
        print(f"\n=== {sc.upper()} SPLIT ===")
        cols = [
            "model",
            "top1_within10",
            "top3_within10",
            "mean_auroc",
            "l1w_top1_within10",
            "l1w_top3_within10",
            "exp_binding",
            "exp_structural",
        ]
        print(summ_df[summ_df["split"] == sc][cols].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
