#!/usr/bin/env python3
"""A tag-site model that owes nothing to ProDomino, EpicTope or INDELi-E.

Depending on third-party predictors costs us three ways. **Coverage**: they ran
on the literature cohort only, so 335 of the 512 panel proteins have none and
the model scores those on an imputed constant. **Licence**: ProDomino is
CC-BY-NC-4.0. **Control**: we cannot fix, retrain or extend something we did not
build.

So this searches for a model over our own features only. The published
predictors are not excluded by principle — they are excluded, and then we look
at what it costs.

The feature pool is wider than the fused set ever used. Blocks:

* **structure** — pLDDT, RSA (freesasa on the AlphaFold model)
* **disorder** — metapredict, three windows
* **conservation** — per-residue entropy and MSA depth. This is the signal
  EpicTope is largely built from, computed by us instead of imported.
* **sequence character** — hydropathy, Pro/Gly fraction, cysteines in a
  10-residue window. The last is the only thing in the model that can see a
  disulfide at all.
* **contacts** — ligand/binding-site proximity. CD20 is the case for this:
  nothing else in the model knows a residue is an antibody epitope.
* **glycosylation** — sequon presence and distance
* **geometry** — loop length, distance to TM and to the terminus, position
* **topology** — TM count, ECD/ICD lengths, signal-peptide length

Three axes are searched together, because they interact: the feature set, the
estimator, and whether features are re-expressed **within each protein**. That
last one matters most in principle — the task is to rank residues inside one
protein, so between-protein variation is noise the model cannot use.

    uv run --locked python independent_model_search.py
    uv run --locked python independent_model_search.py --stage full
"""

from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import (
    ExtraTreesClassifier,
    HistGradientBoostingClassifier,
    RandomForestClassifier,
    StackingClassifier,
)
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[3] / "src"))

from train_site_models import CATEGORICAL, N_SPLITS, SEED  # noqa: E402

STRUCTURE = ["plddt", "plddt_w5", "rsa", "rsa_w5"]
DISORDER = ["disorder", "disorder_w5", "disorder_w15"]
CONSERVATION = ["cons_entropy", "cons_entropy_w5", "cons_entropy_w15", "msa_depth"]
SEQUENCE = ["kd_hydropathy_w7", "pro_gly_frac_w7", "cys_in_w10"]
CONTACTS = [
    "contact_here",
    "dist_to_contact",
    "n_contact_w10",
    "contact_ecd_here",
    "dist_to_contact_ecd",
    "dist_to_pocket",
]
GLYC = ["glyc_sequon_here", "dist_to_glyc", "n_glyc_w10"]
GEOMETRY = ["loop_len", "dist_to_tm", "dist_to_terminus", "rel_position"]
TOPOLOGY = [
    "tm_helix_count",
    "signal_peptide_length",
    "ecd_length_residues",
    "icd_length_residues",
    "protein_length",
    "n_extracellular_regions",
]
HOMOMER = ["homomer_stoichiometry", "homomer_higher_order"]

SETS: dict[str, list[str]] = {
    "structure+disorder": STRUCTURE + DISORDER,
    "+ sequence": STRUCTURE + DISORDER + SEQUENCE,
    "+ conservation": STRUCTURE + DISORDER + SEQUENCE + CONSERVATION,
    "+ contacts": STRUCTURE + DISORDER + SEQUENCE + CONTACTS,
    "+ geometry": STRUCTURE + DISORDER + SEQUENCE + GEOMETRY,
    "everything of ours": (
        STRUCTURE
        + DISORDER
        + SEQUENCE
        + CONSERVATION
        + CONTACTS
        + GLYC
        + GEOMETRY
        + TOPOLOGY
        + HOMOMER
    ),
    "ours, no conservation": (
        STRUCTURE + DISORDER + SEQUENCE + CONTACTS + GLYC + GEOMETRY + TOPOLOGY + HOMOMER
    ),
}


def within_protein(df: pd.DataFrame, feats: list[str], how: str) -> pd.DataFrame:
    if how == "raw":
        return df
    out = df.copy()
    g = out.groupby("uniprot_acc")
    for f in feats:
        if f not in out.columns:
            continue
        if how == "rank":
            out[f] = g[f].rank(pct=True)
        else:
            mu = g[f].transform("mean")
            sd = g[f].transform("std").replace(0, np.nan)
            out[f] = (out[f] - mu) / sd
    return out


def make(kind: str, feats: list[str]) -> Pipeline:
    dense = kind in ("logistic", "stack")
    num = (
        Pipeline([("imp", SimpleImputer(strategy="median")), ("sc", StandardScaler())])
        if dense
        else SimpleImputer(strategy="median")
    )
    pre = ColumnTransformer(
        [
            ("num", num, feats),
            ("cat", OneHotEncoder(handle_unknown="ignore", min_frequency=5), CATEGORICAL),
        ]
    )
    clf = {
        "logistic": LogisticRegression(max_iter=5000, C=0.5, random_state=SEED),
        "rf": RandomForestClassifier(
            n_estimators=600,
            max_depth=6,
            min_samples_leaf=8,
            class_weight="balanced_subsample",
            random_state=SEED,
            n_jobs=-1,
        ),
        "rf_deep": RandomForestClassifier(
            n_estimators=800,
            max_depth=12,
            min_samples_leaf=4,
            class_weight="balanced_subsample",
            random_state=SEED,
            n_jobs=-1,
        ),
        "rf_shallow": RandomForestClassifier(
            n_estimators=800,
            max_depth=4,
            min_samples_leaf=12,
            class_weight="balanced_subsample",
            random_state=SEED,
            n_jobs=-1,
        ),
        "extra_trees": ExtraTreesClassifier(
            n_estimators=800,
            max_depth=8,
            min_samples_leaf=6,
            class_weight="balanced_subsample",
            random_state=SEED,
            n_jobs=-1,
        ),
        "hgb": HistGradientBoostingClassifier(
            max_depth=3, learning_rate=0.05, max_iter=400, random_state=SEED
        ),
        "stack": StackingClassifier(
            estimators=[
                (
                    "rf",
                    RandomForestClassifier(
                        n_estimators=400,
                        max_depth=6,
                        min_samples_leaf=8,
                        class_weight="balanced_subsample",
                        random_state=SEED,
                        n_jobs=-1,
                    ),
                ),
                (
                    "hgb",
                    HistGradientBoostingClassifier(
                        max_depth=3, learning_rate=0.05, max_iter=250, random_state=SEED
                    ),
                ),
                ("lr", LogisticRegression(max_iter=5000, C=0.5, random_state=SEED)),
            ],
            final_estimator=LogisticRegression(max_iter=2000, random_state=SEED),
            cv=3,
            n_jobs=1,
        ),
    }[kind]
    return Pipeline([("pre", pre), ("clf", clf)])


def smooth(frame: pd.DataFrame, scores: np.ndarray, half: int) -> np.ndarray:
    if half <= 0:
        return scores
    out = scores.copy()
    tmp = frame[["uniprot_acc", "residue"]].copy()
    tmp["s"] = scores
    for _, g in tmp.groupby("uniprot_acc"):
        g = g.sort_values("residue")
        vals = g["s"].to_numpy()
        pad = np.pad(vals, half, mode="edge")
        ker = np.ones(2 * half + 1) / (2 * half + 1)
        out[g.index.to_numpy()] = np.convolve(pad, ker, mode="valid")
    return out


def load():
    df = pd.read_csv(HERE / "site_features.tsv", sep="\t", low_memory=False)
    allx = pd.read_csv(HERE / "all_extracellular_targets.tsv", sep="\t", low_memory=False)
    for d in (df, allx):
        for c in ("ss3", "topology_class", "region_class"):
            d[c] = d[c].fillna("").astype(str)
        d["gene_symbol"] = d["gene_symbol"].fillna(d["uniprot_acc"]).astype(str)
    truth = (
        df[df["label"] == 1]
        .groupby("uniprot_acc")["residue"]
        .apply(lambda v: {int(x) for x in v})
        .to_dict()
    )
    ev = allx[allx["uniprot_acc"].isin(truth)].copy().reset_index(drop=True)
    ev_y = np.array(
        [
            1 if int(r) in truth.get(a, ()) else 0
            for a, r in zip(ev["uniprot_acc"], ev["residue"], strict=False)
        ]
    )
    return df, ev, ev_y


def evaluate(df, ev, ev_y, feats, kind, norm, smoothing=(0, 2, 5)):
    y = df["label"].to_numpy()
    groups = df["gene_symbol"].to_numpy()
    tr_df = within_protein(df, feats, norm)
    ev_df = within_protein(ev, feats, norm)
    gg = ev["gene_symbol"].to_numpy()
    m = make(kind, feats)
    oof = np.full(len(df), np.nan)
    sc = np.full(len(ev), np.nan)
    for tr, te in GroupKFold(n_splits=min(N_SPLITS, len(np.unique(groups)))).split(
        tr_df, y, groups
    ):
        if len(np.unique(y[tr])) < 2:
            continue
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            m.fit(tr_df.iloc[tr], y[tr])
        oof[te] = m.predict_proba(tr_df.iloc[te])[:, 1]
        mask = np.isin(gg, list(set(tr_df.iloc[te]["gene_symbol"])))
        if mask.any():
            sc[mask] = m.predict_proba(ev_df[mask])[:, 1]
    ok, ok2 = ~np.isnan(oof), ~np.isnan(sc)
    a = float(roc_auc_score(y[ok], oof[ok]))
    out = []
    for half in smoothing:
        s2 = smooth(ev, sc, half)
        out.append((half, a, float(roc_auc_score(ev_y[ok2], s2[ok2])), sc))
    return out


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--stage", choices=["sets", "full"], default="full")
    ap.add_argument("--out", default="independent_model_search.tsv")
    args = ap.parse_args()

    df, ev, ev_y = load()
    print(
        f"{len(df)} matched rows, {int(df['label'].sum())} positives; "
        f"{len(ev)} extracellular residues on {ev['uniprot_acc'].nunique()} proteins\n"
    )

    rows = []
    print("stage 1 — feature sets, random forest, raw / within-protein rank / z")
    for sname, feats in SETS.items():
        for norm in ("raw", "rank", "z"):
            for half, a, t, _ in evaluate(df, ev, ev_y, feats, "rf", norm, (0,)):
                rows.append(
                    {
                        "stage": "sets",
                        "feature_set": sname,
                        "n_features": len(feats),
                        "model": "rf",
                        "normalisation": norm,
                        "smoothing_half": half,
                        "matched_auroc": round(a, 4),
                        "all_extracellular_auroc": round(t, 4),
                    }
                )
                print(
                    f"  {sname:<24} {norm:<5} {len(feats):>2}f   matched {a:.4f}   all-EC {t:.4f}",
                    flush=True,
                )

    if args.stage == "full":
        best = sorted(rows, key=lambda r: -r["all_extracellular_auroc"])[:2]
        print("\nstage 2 — the two best sets, every estimator, every normalisation")
        for b in best:
            feats = SETS[b["feature_set"]]
            for kind in (
                "logistic",
                "rf",
                "rf_shallow",
                "rf_deep",
                "extra_trees",
                "hgb",
                "stack",
            ):
                for norm in ("raw", "rank", "z"):
                    for half, a, t, _ in evaluate(df, ev, ev_y, feats, kind, norm):
                        rows.append(
                            {
                                "stage": "full",
                                "feature_set": b["feature_set"],
                                "n_features": len(feats),
                                "model": kind,
                                "normalisation": norm,
                                "smoothing_half": half,
                                "matched_auroc": round(a, 4),
                                "all_extracellular_auroc": round(t, 4),
                            }
                        )
                    best_half = max(
                        [
                            r
                            for r in rows
                            if r["stage"] == "full"
                            and r["model"] == kind
                            and r["normalisation"] == norm
                            and r["feature_set"] == b["feature_set"]
                        ],
                        key=lambda r: r["all_extracellular_auroc"],
                    )
                    print(
                        f"  {b['feature_set']:<24} {kind:<12} {norm:<5} "
                        f"best smooth ±{best_half['smoothing_half']}   "
                        f"matched {best_half['matched_auroc']:.4f}   "
                        f"all-EC {best_half['all_extracellular_auroc']:.4f}",
                        flush=True,
                    )

    out = pd.DataFrame(rows).sort_values("all_extracellular_auroc", ascending=False)
    out.to_csv(HERE / args.out, sep="\t", index=False)
    print("\ntop 12 by vs-all-extracellular:")
    print(out.head(12).to_string(index=False))
    print(f"\nwrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
