#!/usr/bin/env python3
"""Can a fused model rank published internal tag sites above matched controls?

Trains and cross-validates several classifiers on the feature table from
``export_site_features.py``: published internal junctions (label 1) against
region-matched extracellular residues on the same protein (label 0).

**Read every number here as ranking against matched controls, not as accuracy.**
The negatives are untried residues, not failures — this is positive-unlabeled
data, so precision and calibration are uninterpretable and are not reported.
AUROC answers the one question the matching makes answerable: given a published
site and a control drawn from the same protein and the same region class, how
often does the model put the published one first? Chance is 0.5 by construction.

**Folds are grouped by gene.** A protein contributes its positives and all of
their matched negatives to exactly one fold. Splitting by site would put
near-identical rows on both sides of the split and inflate every model alike.

**One pooled model, topology class as a feature.** With 185 positives — 122
multipass, 48 type I, and single digits elsewhere (``type_ii`` does not occur at
all) — a model per class would be fitting noise in most cells. Class enters as a
feature and comes back out as a reported stratum.

    uv run --locked python train_site_models.py
    uv run --locked python train_site_models.py --figures
"""

from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import HistGradientBoostingClassifier, RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[3] / "src"))

N_SPLITS = 5
SEED = 20261013

# Below this many positives a cell is reported, not fitted. Declared up front so
# it is not chosen after seeing which cells would have looked good.
MIN_POSITIVES = 40

NUMERIC = [
    "prodomino",
    "prodomino_w5",
    "indeli_e",
    "indeli_e_w5",
    "epictope",
    "epictope_w5",
    "disorder",
    "disorder_w5",
    "disorder_w15",
    "rsa",
    "rsa_w5",
    "plddt",
    "plddt_w5",
    "loop_len",
    "tm_helix_count",
    "signal_peptide_length",
    "ecd_length_residues",
    "icd_length_residues",
    "dist_to_tm",
    "dist_to_terminus",
    "rel_position",
    "protein_length",
    # N-glycosylation sequons, computed from sequence so coverage is complete.
    "glyc_sequon_here",
    "dist_to_glyc",
    "n_glyc_w10",
    # Structure-derived ligand / binder contacts, and predicted pockets.
    "contact_here",
    "dist_to_contact",
    "n_contact_w10",
    "contact_ecd_here",
    "dist_to_contact_ecd",
    "dist_to_pocket",
    # Homomer state is protein-level only; per-residue interface burial would
    # need the Schweke oligomer models and is not represented.
    "homomer_stoichiometry",
    "homomer_higher_order",
    # Local sequence character of the host stretch.
    "kd_hydropathy_w7",
    "pro_gly_frac_w7",
    "cys_in_w10",
    # Sparse: conservation exists only where a deterministic lane claims the
    # residue (~5% of rows), so this cannot stand in for EpicTope's conservation.
    "lane_conservation",
    # EpicTope's column entropy, per residue. Low entropy = conserved. Depth
    # is a feature in its own right: entropy over 7 aligned sequences is a
    # weaker statement than entropy over 7 when some are missing.
    "cons_entropy",
    "cons_entropy_w5",
    "cons_entropy_w15",
    "msa_depth",
    "lane_disorder",
    "lane_surface_loop",
    "lane_terminal",
    "lane_snorkel",
]
CATEGORICAL = ["ss3", "topology_class", "region_class"]

# The features that separate published sites from matched controls on their own
# (univariate |AUROC - 0.5| >= 0.09, measured on this table). Everything else in
# NUMERIC is kept for the record but carries no univariate signal here, and at
# 176 positives a column that carries none is a column that costs.
LEAN_NUMERIC = [
    "plddt",
    "plddt_w5",
    "prodomino",
    "prodomino_w5",
    "epictope",
    "epictope_w5",
    "rsa",
    "rsa_w5",
    "disorder",
    "disorder_w5",
    "disorder_w15",
    "loop_len",
    "cons_entropy",
    "cons_entropy_w5",
    "cons_entropy_w15",
    "kd_hydropathy_w7",
]

# ---------------------------------------------------------------------------
# Ablation blocks. "Does adding these features help?" is only answerable by
# adding them to a fixed baseline and holding everything else constant —
# comparing `full` against `lean` does not answer it, because `lean` also drops
# baseline columns.
# ---------------------------------------------------------------------------
BASELINE_NUMERIC = [
    "prodomino",
    "prodomino_w5",
    "indeli_e",
    "indeli_e_w5",
    "epictope",
    "epictope_w5",
    "disorder",
    "disorder_w5",
    "disorder_w15",
    "rsa",
    "rsa_w5",
    "plddt",
    "plddt_w5",
    "loop_len",
    "tm_helix_count",
    "signal_peptide_length",
    "ecd_length_residues",
    "icd_length_residues",
    "dist_to_tm",
    "dist_to_terminus",
    "rel_position",
    "protein_length",
    "lane_disorder",
    "lane_surface_loop",
    "lane_terminal",
    "lane_snorkel",
]

# INDELi-E is excluded from the model that SCORES anything. Dropping it
# improves out-of-fold AUROC by 0.0006 — it earns nothing — and the reason is
# structural: its stored arrays are saturated (89.8% of residues are byte 255),
# so it carries 104 distinct values across the 64% of training rows it covers,
# against ProDomino's 777 over 90%. The logistic imputes a missing numeric with
# the column median, so for everything it does not cover it contributes a
# constant with a fitted coefficient. See published_predictor_ablation.py.
#
# BASELINE_NUMERIC keeps it, because the benchmark and ablation scripts have to
# be able to measure what it was worth. Nothing that produces a score should
# use that list.
EXCLUDED_FROM_SCORING = ["indeli_e", "indeli_e_w5"]
FUSED_NUMERIC = [f for f in BASELINE_NUMERIC if f not in EXCLUDED_FROM_SCORING]

CONSERVATION_BLOCK = ["cons_entropy", "cons_entropy_w5", "cons_entropy_w15", "msa_depth"]
ANNOTATION_BLOCK = [
    "glyc_sequon_here",
    "dist_to_glyc",
    "n_glyc_w10",
    "contact_here",
    "dist_to_contact",
    "n_contact_w10",
    "contact_ecd_here",
    "dist_to_contact_ecd",
    "dist_to_pocket",
    "homomer_stoichiometry",
    "homomer_higher_order",
    "lane_conservation",
]
SEQUENCE_BLOCK = ["kd_hydropathy_w7", "pro_gly_frac_w7", "cys_in_w10"]

# The main scoreboard runs ONE feature set: the baseline, which the ablation
# below shows is the best of them. Comparing feature sets there as well would
# answer the same question twice, in the weaker unpaired form.
FEATURE_SETS = {"baseline": BASELINE_NUMERIC}

ABLATION = {
    "baseline": BASELINE_NUMERIC,
    "+ conservation": BASELINE_NUMERIC + CONSERVATION_BLOCK,
    "+ glycan/contact/homomer": BASELINE_NUMERIC + ANNOTATION_BLOCK,
    "+ sequence character": BASELINE_NUMERIC + SEQUENCE_BLOCK,
    "+ everything": BASELINE_NUMERIC + CONSERVATION_BLOCK + ANNOTATION_BLOCK + SEQUENCE_BLOCK,
    "lean": LEAN_NUMERIC,
}

# Single-feature rankers. No fitting: the raw score IS the prediction, which is
# exactly how these predictors were compared in the s5e13 scoreboard.
BASELINES = {
    "ProDomino (alone)": "prodomino",
    "INDELi-E (alone)": "indeli_e",
    "EpicTope (alone)": "epictope",
    "ours — disorder (alone)": "disorder",
    "ours — RSA (alone)": "rsa",
}

SUBSETS = {
    # name: (predicate on the positives, what it is)
    "all internal": (None, "every published internal junction"),
    "surface_and_function": (
        lambda d: d["validation_level"] == "surface_and_function",
        "tag shown on the surface AND protein function measured intact",
    ),
    "surface-confirmed": (
        lambda d: d["validation_level"].isin(["surface_and_function", "surface_only"]),
        "tag shown on the extracellular face, function measured or not",
    ),
    "position validated": (
        lambda d: d["position_evidence"] == "validated",
        "the paper pins this residue, rather than the region",
    ),
    "exact-site evidence": (
        lambda d: d["evidence_type"] == "published tag insertion at this exact site",
        "insertion reported at this very residue",
    ),
}


def _xgboost_or_none():
    """xgboost, or None with one explanatory line.

    The macOS wheel dynamically links Homebrew's libomp, which is not a Python
    dependency and cannot be pinned in uv.lock. HistGradientBoosting is the same
    algorithm family and vendors its own OpenMP, so the run proceeds without
    xgboost rather than stopping — but it says so, because "xgboost" silently
    absent from a scoreboard reads as "xgboost lost".
    """
    try:
        import xgboost as xgb

        return xgb
    except Exception as e:
        msg = str(e).splitlines()[0] if str(e) else type(e).__name__
        if not getattr(_xgboost_or_none, "_warned", False):
            print(
                f"  NOTE: xgboost unavailable ({msg.strip() or type(e).__name__}); "
                "its slot is left out of the scoreboard.\n"
                "        Enable it with: brew install libomp",
                file=sys.stderr,
            )
            _xgboost_or_none._warned = True  # type: ignore[attr-defined]
        return None


def models(n_pos: int, numeric: list[str] | None = None) -> dict:
    """The fitted model zoo.

    Tree depths are kept small on purpose: with ~100-185 positives a deep
    ensemble memorises the proteins rather than the signal.
    """
    xgb = _xgboost_or_none()
    num = numeric or NUMERIC

    pre_dense = ColumnTransformer(
        [
            (
                "num",
                Pipeline([("imp", SimpleImputer(strategy="median")), ("sc", StandardScaler())]),
                num,
            ),
            ("cat", OneHotEncoder(handle_unknown="ignore", min_frequency=5), CATEGORICAL),
        ]
    )
    # Trees take NaN natively, so a protein that was never run through INDELi-E
    # keeps a missing value instead of being handed the median as if measured.
    pre_sparse = ColumnTransformer(
        [
            ("num", "passthrough", num),
            ("cat", OneHotEncoder(handle_unknown="ignore", min_frequency=5), CATEGORICAL),
        ]
    )
    zoo = {
        "logistic": Pipeline(
            [
                ("pre", pre_dense),
                ("clf", LogisticRegression(max_iter=5000, C=0.5, random_state=SEED)),
            ]
        ),
        "random forest": Pipeline(
            [
                (
                    "pre",
                    ColumnTransformer(
                        [
                            ("num", SimpleImputer(strategy="median"), num),
                            (
                                "cat",
                                OneHotEncoder(handle_unknown="ignore", min_frequency=5),
                                CATEGORICAL,
                            ),
                        ]
                    ),
                ),
                (
                    "clf",
                    RandomForestClassifier(
                        n_estimators=600,
                        min_samples_leaf=5,
                        max_features="sqrt",
                        random_state=SEED,
                        n_jobs=-1,
                    ),
                ),
            ]
        ),
        "hist gradient boosting": Pipeline(
            [
                ("pre", pre_sparse),
                (
                    "clf",
                    HistGradientBoostingClassifier(
                        max_depth=3,
                        max_iter=300,
                        learning_rate=0.05,
                        min_samples_leaf=10,
                        random_state=SEED,
                    ),
                ),
            ]
        ),
    }
    if xgb is not None:
        zoo["xgboost"] = Pipeline(
            [
                ("pre", pre_sparse),
                (
                    "clf",
                    xgb.XGBClassifier(
                        n_estimators=400,
                        max_depth=3,
                        learning_rate=0.05,
                        subsample=0.8,
                        colsample_bytree=0.8,
                        reg_lambda=2.0,
                        min_child_weight=5,
                        eval_metric="logloss",
                        random_state=SEED,
                        n_jobs=-1,
                    ),
                ),
            ]
        )
    return zoo


def oof_predictions(df: pd.DataFrame, model) -> tuple[np.ndarray, list[float], list[float]]:
    """Out-of-fold scores plus per-fold AUROC / average precision.

    Returns OOF scores aligned to ``df``'s order so strata can be scored from
    one honest set of held-out predictions.
    """
    y = df["label"].to_numpy()
    groups = df["gene_symbol"].to_numpy()
    oof = np.full(len(df), np.nan)
    aurocs: list[float] = []
    aps: list[float] = []
    n_splits = min(N_SPLITS, len(np.unique(groups)))
    for tr, te in GroupKFold(n_splits=n_splits).split(df, y, groups):
        if len(np.unique(y[tr])) < 2:
            continue
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            model.fit(df.iloc[tr], y[tr])
            p = model.predict_proba(df.iloc[te])[:, 1]
        oof[te] = p
        if len(np.unique(y[te])) == 2:
            aurocs.append(roc_auc_score(y[te], p))
            aps.append(average_precision_score(y[te], p))
    return oof, aurocs, aps


def baseline_scores(df: pd.DataFrame, col: str) -> np.ndarray:
    """Raw feature as a ranker; missing values get the median.

    A predictor that never scored a protein genuinely cannot rank it — handing
    it the median says exactly that, rather than quietly dropping those rows and
    scoring it on an easier subset than everyone else.
    """
    s = pd.to_numeric(df[col], errors="coerce")
    return s.fillna(s.median()).to_numpy()


def fold_metrics_from_scores(
    df: pd.DataFrame, scores: np.ndarray
) -> tuple[list[float], list[float]]:
    y = df["label"].to_numpy()
    groups = df["gene_symbol"].to_numpy()
    aurocs, aps = [], []
    n_splits = min(N_SPLITS, len(np.unique(groups)))
    for _, te in GroupKFold(n_splits=n_splits).split(df, y, groups):
        if len(np.unique(y[te])) == 2:
            aurocs.append(roc_auc_score(y[te], scores[te]))
            aps.append(average_precision_score(y[te], scores[te]))
    return aurocs, aps


def mean_sd(v: list[float]) -> str:
    if not v:
        return "n/a"
    return f"{np.mean(v):.3f} ± {np.std(v, ddof=1) if len(v) > 1 else 0:.3f}"


def subset_frame(df: pd.DataFrame, pred) -> pd.DataFrame:
    """Keep matching positives plus only the negatives matched to them."""
    if pred is None:
        return df
    pos = df[(df["label"] == 1) & pred(df)]
    keep_ids = set(pos["site_id"])
    neg = df[(df["label"] == 0) & df["matched_to"].isin(keep_ids)]
    return pd.concat([pos, neg]).sort_values(["gene_symbol", "residue"]).reset_index(drop=True)


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--features", default="site_features.tsv")
    ap.add_argument("--figures", action="store_true", help="render the comparison figures")
    args = ap.parse_args()

    path = HERE / args.features
    if not path.is_file():
        raise SystemExit(
            f"{path.name} not found. Generate it first:\n"
            "  uv run --locked python export_site_features.py --out site_features.tsv "
            "--env-from <checkout>/.env"
        )
    df = pd.read_csv(path, sep="\t", dtype={"ss3": str, "matched_to": str, "site_id": str})
    for c in CATEGORICAL:
        df[c] = df[c].fillna("").astype(str)
    for c in ("validation_level", "position_evidence", "evidence_type"):
        df[c] = df[c].fillna("").astype(str)

    print(
        f"{len(df)} rows: {int(df['label'].sum())} positives, "
        f"{int((df['label'] == 0).sum())} matched negatives, "
        f"{df['gene_symbol'].nunique()} proteins\n"
    )

    records: list[dict] = []
    oof_store: dict[str, np.ndarray] = {}

    # ---------------- per subset ----------------
    for sub_name, (pred, blurb) in SUBSETS.items():
        d = subset_frame(df, pred)
        n_pos = int(d["label"].sum())
        n_prot = d["gene_symbol"].nunique()
        print(f"=== {sub_name} — {blurb}")
        print(
            f"    {n_pos} positives / {int((d['label'] == 0).sum())} negatives, {n_prot} proteins"
        )
        if n_pos < MIN_POSITIVES:
            print(f"    SKIPPED: fewer than {MIN_POSITIVES} positives — reported, not fitted.\n")
            records.append(
                {
                    "subset": sub_name,
                    "model": "(not fitted)",
                    "n_pos": n_pos,
                    "n_proteins": n_prot,
                    "auroc": np.nan,
                    "auroc_sd": np.nan,
                    "ap": np.nan,
                    "n_folds": 0,
                    "auroc_folds": "",
                }
            )
            continue

        rows = []
        for name, col in BASELINES.items():
            cov = d[col].notna().mean()
            a, p = fold_metrics_from_scores(d, baseline_scores(d, col))
            rows.append((name, a, p, cov))
        # Both feature sets, every time: whether more features help is the
        # question, so it is answered in the same table rather than asserted.
        for fs_name, fs_cols in FEATURE_SETS.items():
            for name, model in models(n_pos, fs_cols).items():
                label = name if len(FEATURE_SETS) == 1 else f"{name} [{fs_name}]"
                oof, a, p = oof_predictions(d, model)
                rows.append((label, a, p, 1.0))
                if sub_name == "all internal":
                    oof_store[label] = oof

        for name, a, p, cov in sorted(rows, key=lambda r: -(np.mean(r[1]) if r[1] else 0)):
            cov_s = f"  (feature on {cov:.0%} of rows)" if cov < 1 else ""
            print(f"      {name:<30} AUROC {mean_sd(a)}   AP {mean_sd(p)}{cov_s}")
            records.append(
                {
                    "subset": sub_name,
                    "model": name,
                    "n_pos": n_pos,
                    "n_proteins": n_prot,
                    "auroc": np.mean(a) if a else np.nan,
                    "auroc_sd": (np.std(a, ddof=1) if len(a) > 1 else 0.0) if a else np.nan,
                    "ap": np.mean(p) if p else np.nan,
                    "n_folds": len(a),
                    # Individual folds, so the figure can show the
                    # observations rather than only their mean.
                    "auroc_folds": ";".join(f"{v:.4f}" for v in a),
                }
            )
        print()

    pd.DataFrame(records).to_csv(HERE / "model_scoreboard.tsv", sep="\t", index=False)
    print(f"wrote model_scoreboard.tsv ({len(records)} rows)")

    # ---------------- strata, from one set of held-out predictions ----------------
    best = (
        max(
            oof_store,
            key=lambda k: roc_auc_score(
                df["label"].to_numpy()[~np.isnan(oof_store[k])],
                oof_store[k][~np.isnan(oof_store[k])],
            ),
        )
        if oof_store
        else None
    )
    strat_records: list[dict] = []
    if best:
        oof = oof_store[best]
        ok = ~np.isnan(oof)
        print(f"\n=== strata for the best pooled model ({best}), out-of-fold ===")
        for col in ("topology_class", "region_class", "validation_level", "position_evidence"):
            print(f"  -- by {col}")
            for val, g in df[ok].groupby(df.loc[ok, col]):
                # A stratum defined on a positive-only column has no negatives of
                # its own; score its positives against all matched negatives.
                if col in ("validation_level", "position_evidence"):
                    pos_idx = g.index[g["label"] == 1]
                    if not len(pos_idx):
                        continue
                    neg_idx = df.index[
                        (df["label"] == 0) & df["matched_to"].isin(df.loc[pos_idx, "site_id"])
                    ]
                    idx = pos_idx.union(neg_idx)
                    idx = idx[~np.isnan(oof[idx])]
                else:
                    idx = g.index
                y = df.loc[idx, "label"].to_numpy()
                if len(np.unique(y)) < 2 or not len(idx):
                    continue
                a = roc_auc_score(y, oof[idx])
                n_p = int(y.sum())
                flag = "" if n_p >= MIN_POSITIVES else f"   (n={n_p}: indicative only)"
                print(f"       {str(val)[:34]:<36} AUROC {a:.3f}  n_pos={n_p}{flag}")
                strat_records.append(
                    {
                        "model": best,
                        "stratum_kind": col,
                        "stratum": val,
                        "auroc": a,
                        "n_pos": n_p,
                        "n_rows": len(idx),
                    }
                )
        pd.DataFrame(strat_records).to_csv(HERE / "model_strata.tsv", sep="\t", index=False)
        print(f"\nwrote model_strata.tsv ({len(strat_records)} rows)")

    # ---------------- what the model leans on ----------------
    d_all = df
    zoo = models(int(d_all["label"].sum()))
    which = "xgboost" if "xgboost" in zoo else "hist gradient boosting"
    try:
        m = zoo[which]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            m.fit(d_all, d_all["label"].to_numpy())
        names = [n.split("__", 1)[-1] for n in m.named_steps["pre"].get_feature_names_out()]
        if hasattr(m.named_steps["clf"], "feature_importances_"):
            score = m.named_steps["clf"].feature_importances_
            kind = "gain"
        else:
            # HistGradientBoosting exposes no gain; permutation on the fitted
            # model is the closest readable substitute.
            from sklearn.inspection import permutation_importance

            Xt = m.named_steps["pre"].transform(d_all)
            r = permutation_importance(
                m.named_steps["clf"],
                Xt,
                d_all["label"].to_numpy(),
                n_repeats=10,
                random_state=SEED,
                scoring="roc_auc",
            )
            score, kind = r.importances_mean, "permutation ΔAUROC"
        imp = (
            pd.DataFrame({"feature": names, "importance": score, "kind": kind})
            .sort_values("importance", ascending=False)
            .reset_index(drop=True)
        )
        imp.to_csv(HERE / "model_feature_importance.tsv", sep="\t", index=False)
        print(
            f"\n=== {which} {kind}, top 15 (fitted on everything; for reading, not for scoring) ==="
        )
        for _, r in imp.head(15).iterrows():
            print(f"   {r['feature']:<28} {r['importance']:.4f}")
    except Exception as e:
        print(f"\nfeature importance skipped: {type(e).__name__}: {e}")

    run_ablation(df)

    if args.figures:
        make_figures(pd.DataFrame(records), pd.DataFrame(strat_records), df, oof_store)
    return 0


def run_ablation(df: pd.DataFrame) -> None:
    """Does each feature block help, added to a fixed baseline?

    Folds are **paired**: GroupKFold is deterministic given the same rows and
    groups, so every feature set is scored on identical splits and the per-fold
    difference is a like-for-like comparison. That matters at five folds — the
    unpaired sd (~0.05) swamps differences the paired delta resolves cleanly,
    because the folds themselves differ far more than the feature sets do.
    """
    print("\n" + "=" * 72)
    print("ABLATION — does each block help, added to the same baseline?")
    print("=" * 72)

    records: list[dict] = []
    for sub_name, (pred, _) in SUBSETS.items():
        if sub_name not in ("all internal", "surface_and_function", "exact-site evidence"):
            continue
        d = subset_frame(df, pred)
        n_pos = int(d["label"].sum())
        print(f"\n--- {sub_name} ({n_pos} positives)")
        for model_name in ("random forest", "xgboost", "logistic"):
            base_folds: list[float] | None = None
            print(f"    {model_name}")
            for fs_name, cols in ABLATION.items():
                zoo = models(n_pos, cols)
                if model_name not in zoo:
                    continue
                _, a, _ = oof_predictions(d, zoo[model_name])
                if fs_name == "baseline":
                    base_folds = a
                    print(f"      {fs_name:<26} AUROC {mean_sd(a)}")
                    records.append(
                        {
                            "subset": sub_name,
                            "model": model_name,
                            "feature_set": fs_name,
                            "auroc": np.mean(a),
                            "delta_vs_baseline": 0.0,
                            "delta_sd": 0.0,
                            "n_folds": len(a),
                            "n_pos": n_pos,
                            # The baseline is its own reference, so every fold
                            # ties. Written explicitly: an empty trailing field
                            # makes a ragged TSV once whitespace is stripped.
                            "folds_improved": 0,
                        }
                    )
                    continue
                if base_folds and len(base_folds) == len(a):
                    deltas = [x - b for x, b in zip(a, base_folds, strict=False)]
                    md = float(np.mean(deltas))
                    sd = float(np.std(deltas, ddof=1)) if len(deltas) > 1 else 0.0
                    wins = sum(1 for x in deltas if x > 0)
                    mark = "  <-- helps" if md > 0 and wins >= len(deltas) - 1 else ""
                    print(
                        f"      {fs_name:<26} AUROC {mean_sd(a)}   "
                        f"paired Δ {md:+.4f} ± {sd:.4f}  ({wins}/{len(deltas)} folds up){mark}"
                    )
                    records.append(
                        {
                            "subset": sub_name,
                            "model": model_name,
                            "feature_set": fs_name,
                            "auroc": np.mean(a),
                            "delta_vs_baseline": md,
                            "delta_sd": sd,
                            "n_folds": len(a),
                            "n_pos": n_pos,
                            "folds_improved": wins,
                        }
                    )
    pd.DataFrame(records).to_csv(HERE / "model_ablation.tsv", sep="\t", index=False)
    print(f"\nwrote model_ablation.tsv ({len(records)} rows)")


def make_figures(
    scoreboard: pd.DataFrame,
    strata: pd.DataFrame,
    df: pd.DataFrame,
    oof_store: dict[str, np.ndarray],
) -> None:
    """House style only: setup_plotting_style() then save_figure()."""
    import matplotlib.pyplot as plt

    from deliverome_analysis.core.plotting import save_figure, setup_plotting_style

    setup_plotting_style()
    figdir = HERE / "figures"

    # --- model comparison on the full internal set, one point per CV fold ---
    # Dot plot, not bars: AUROC is referenced to 0.5, and a bar drawn from a
    # truncated axis makes a 0.05 difference look like a doubling.
    sub = scoreboard[(scoreboard["subset"] == "all internal") & scoreboard["auroc"].notna()]
    sub = sub.sort_values("auroc")
    fig, ax = plt.subplots(figsize=(7.4, 4.6))
    ypos = np.arange(len(sub))
    for y, folds in zip(ypos, sub["auroc_folds"], strict=False):
        vals = [float(v) for v in str(folds).split(";") if v]
        if vals:
            ax.scatter(
                vals,
                np.full(len(vals), y),
                s=16,
                zorder=2,
                alpha=0.55,
                facecolor="none",
                edgecolor="0.35",
                linewidth=0.8,
            )
    ax.errorbar(
        sub["auroc"], ypos, xerr=sub["auroc_sd"], fmt="o", ms=6, lw=1.2, capsize=3, zorder=3
    )
    ax.axvline(0.5, ls="--", lw=1.0, color="0.35")
    ax.text(
        0.505,
        len(sub) - 0.25,
        "chance (matched controls)",
        fontsize=7,
        color="0.35",
        va="center",
        ha="left",
    )
    ax.set_yticks(ypos)
    ax.set_yticklabels(sub["model"])
    ax.set_ylim(-0.6, len(sub) - 0.1)
    ax.set_xlabel("AUROC, published internal site vs matched control")
    ax.set_xlim(0.44, 0.85)
    save_figure(fig, figdir / "s5e13_model_auroc")
    plt.close(fig)

    # --- AUROC by subset, for the fitted models only ---
    fitted = scoreboard[
        scoreboard["auroc"].notna() & ~scoreboard["model"].str.contains(r"\(alone\)")
    ]
    if not fitted.empty:
        fig, ax = plt.subplots(figsize=(7.8, 4.2))
        subsets = list(dict.fromkeys(fitted["subset"]))
        mods = list(dict.fromkeys(fitted["model"]))
        w = 0.78 / max(1, len(mods))
        x = np.arange(len(subsets))
        for i, mod in enumerate(mods):
            g = fitted[fitted["model"] == mod].set_index("subset").reindex(subsets)
            xs = x + i * w - 0.39 + w / 2
            for xi, folds in zip(xs, g["auroc_folds"].fillna(""), strict=False):
                vals = [float(v) for v in str(folds).split(";") if v]
                if vals:
                    ax.scatter(
                        np.full(len(vals), xi),
                        vals,
                        s=7,
                        zorder=2,
                        alpha=0.5,
                        facecolor="none",
                        edgecolor="0.35",
                        linewidth=0.5,
                    )
            ax.errorbar(
                xs,
                g["auroc"],
                yerr=g["auroc_sd"],
                fmt="o",
                ms=4.5,
                lw=1.0,
                capsize=2,
                label=mod,
                zorder=3,
            )
        ax.axhline(0.5, ls="--", lw=1.0, color="0.35")
        ax.set_xticks(x)
        npos = {s: int(fitted[fitted["subset"] == s]["n_pos"].iloc[0]) for s in subsets}
        ax.set_xticklabels([f"{s}\n(n={npos[s]})" for s in subsets], fontsize=7)
        ax.set_ylabel("AUROC")
        ax.set_ylim(0.44, 0.92)
        ax.legend(fontsize=7, frameon=False, ncol=2)
        save_figure(fig, figdir / "s5e13_model_auroc_by_subset")
        plt.close(fig)

    # --- strata for the best model ---
    if not strata.empty:
        s = strata[strata["stratum_kind"].isin(["topology_class", "region_class"])].copy()
        if not s.empty:
            s["lab"] = s["stratum"].astype(str) + " (n=" + s["n_pos"].astype(str) + ")"
            s = s.sort_values("auroc")
            fig, ax = plt.subplots(figsize=(6.8, 3.8))
            # Open markers for cells under the declared floor: they are shown so
            # the reader sees every stratum, not because they can carry a claim.
            thin = s["n_pos"] < MIN_POSITIVES
            ax.scatter(s.loc[~thin, "auroc"], np.arange(len(s))[~thin.to_numpy()], s=48, zorder=3)
            ax.scatter(
                s.loc[thin, "auroc"],
                np.arange(len(s))[thin.to_numpy()],
                s=48,
                zorder=3,
                facecolor="none",
                edgecolor="0.4",
                linewidth=1.1,
            )
            ax.axvline(0.5, ls="--", lw=1.0, color="0.35")
            ax.set_yticks(np.arange(len(s)))
            ax.set_yticklabels(s["lab"], fontsize=8)
            ax.set_xlabel("AUROC, out-of-fold")
            ax.set_xlim(0.44, 0.85)
            save_figure(fig, figdir / "s5e13_model_auroc_by_stratum")
            plt.close(fig)
    # --- feature ablation: every model x subset cell, against the baseline ---
    apath = HERE / "model_ablation.tsv"
    if apath.is_file():
        abl = pd.read_csv(apath, sep="\t")
        abl = abl[abl["feature_set"] != "baseline"]
        if not abl.empty:
            blocks = list(dict.fromkeys(abl["feature_set"]))
            fig, ax = plt.subplots(figsize=(7.4, 4.0))
            mods = list(dict.fromkeys(abl["model"]))
            jit = np.linspace(-0.17, 0.17, max(1, len(mods)))
            for mi, mod in enumerate(mods):
                g = abl[abl["model"] == mod]
                ax.scatter(
                    g["delta_vs_baseline"],
                    [blocks.index(b) + jit[mi] for b in g["feature_set"]],
                    s=26,
                    label=mod,
                    zorder=3,
                )
            # The mean over all cells, so the eye is not drawn to one outlier.
            m = abl.groupby("feature_set")["delta_vs_baseline"].mean()
            ax.scatter(
                [m[b] for b in blocks],
                range(len(blocks)),
                marker="D",
                s=46,
                facecolor="none",
                edgecolor="0.15",
                linewidth=1.3,
                zorder=4,
                label="mean of cells",
            )
            ax.axvline(0, ls="-", lw=1.3, color="0.45")
            ax.text(0.001, -0.62, "baseline", fontsize=7, color="0.4", ha="left")
            ax.set_yticks(range(len(blocks)))
            ax.set_yticklabels(blocks, fontsize=8)
            ax.invert_yaxis()
            ax.set_xlabel("change in AUROC vs baseline (paired folds)")
            ax.legend(fontsize=7, frameon=False, loc="lower left")
            save_figure(fig, figdir / "s5e13_feature_ablation")
            plt.close(fig)

    print("\nwrote figures to figures/")


if __name__ == "__main__":
    raise SystemExit(main())
