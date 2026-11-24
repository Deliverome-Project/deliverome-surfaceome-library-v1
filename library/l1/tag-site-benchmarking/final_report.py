#!/usr/bin/env python3
"""The morning report: every candidate model, fold by fold, against the baselines.

Each variant is scored three ways, because no single number decides it:

1. **AUROC, 5-fold GroupKFold by gene**, vs every extracellular residue of the
   protein, with negatives within +/-10 of a published site masked out. Fold
   values are kept, not just the mean -- with 176 positives the spread is wide
   enough that a difference under ~0.05 is not a difference.
2. **Top-k recall**: how often the published site lands in the protein's top 1,
   3 or 10 residues. This is closer to how the model is used than AUROC is --
   you order a handful of sites, not a ranking.
3. **SLC6A8**, which is held out in the strongest sense: it contributes zero
   rows to any training table, and Jon's HA substitution at 193-201 works in
   the lab. Where the model puts 193 is the only check here against a real
   result rather than against publication.

ProDomino and EpicTope are scored identically, on the residues they cover.

    uv run --locked python final_report.py
"""

from __future__ import annotations

import json
import sys
import warnings
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[3] / "src"))

from independent_model_search import SETS, load, make, smooth, within_protein  # noqa: E402
from train_site_models import FUSED_NUMERIC, N_SPLITS  # noqa: E402
from two_tier_model import (  # noqa: E402
    CONTACTS,
    GEOMETRY,
    GLYC,
    HOMOMER,
    MASK_HALF,
    add_disulfide,
    build_global,
)

SEED = 20261014
RATIO = 200
OURS = SETS["+ conservation"]
SLC6A8 = "P48029"
JON = range(193, 202)


def variants(cys_new):
    v = {
        "ours 14": (OURS, "hgb"),
        "ours + glyc (17)": (OURS + GLYC, "hgb"),
        "ours + glyc + disulfide (22)": (OURS + GLYC + cys_new, "hgb"),
        "everything ours (34)": (OURS + GLYC + CONTACTS + HOMOMER + GEOMETRY + cys_new, "hgb"),
        "ours + glyc, forest": (OURS + GLYC, "rf_shallow"),
        "fused (with third-party)": (FUSED_NUMERIC, "hgb"),
    }
    return v


def fold_eval(train, ev, ev_y, feats, kind):
    y = train["label"].to_numpy()
    groups = train["gene_symbol"].to_numpy()
    tr = within_protein(train, feats, "z")
    ed = within_protein(ev, feats, "z")
    gg = ev["gene_symbol"].to_numpy()
    sc = np.full(len(ev), np.nan)
    fold = np.full(len(ev), -1)
    for i, (a, b) in enumerate(GroupKFold(n_splits=N_SPLITS).split(tr, y, groups)):
        if len(np.unique(y[a])) < 2:
            continue
        m = make(kind, feats)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            m.fit(tr.iloc[a], y[a])
        mask = np.isin(gg, list(set(tr.iloc[b]["gene_symbol"])))
        if mask.any():
            sc[mask] = m.predict_proba(ed[mask])[:, 1]
            fold[mask] = i
    return smooth(ev, sc, 2), fold


def mask_keep(ev, ev_y, half=MASK_HALF):
    keep = np.ones(len(ev), bool)
    byacc: dict[str, list[int]] = {}
    for a, r, lab in zip(ev["uniprot_acc"], ev["residue"], ev_y, strict=False):
        if lab == 1:
            byacc.setdefault(a, []).append(int(r))
    for i, (a, r, lab) in enumerate(zip(ev["uniprot_acc"], ev["residue"], ev_y, strict=False)):
        if lab != 1 and byacc.get(a) and min(abs(int(r) - p) for p in byacc[a]) <= half:
            keep[i] = False
    return keep


def metrics(ev, ev_y, sc, fold, keep):
    ok = keep & ~np.isnan(sc)
    folds = []
    for f in sorted(set(fold)):
        if f < 0:
            continue
        k = ok & (fold == f)
        if k.sum() > 50 and len(np.unique(ev_y[k])) == 2:
            folds.append(roc_auc_score(ev_y[k], sc[k]))
    folds = np.array(folds)
    # top-k recall, per protein
    hits = {1: 0, 3: 0, 10: 0}
    n_prot = 0
    d = pd.DataFrame({"acc": ev["uniprot_acc"], "res": ev["residue"], "y": ev_y, "s": sc})[ok]
    for _, g in d.groupby("acc"):
        if g["y"].sum() == 0:
            continue
        n_prot += 1
        order = g.sort_values("s", ascending=False)
        pos_ranks = np.where(order["y"].to_numpy() == 1)[0] + 1
        for k in hits:
            if pos_ranks.min() <= k:
                hits[k] += 1
    return folds, {k: v / max(n_prot, 1) for k, v in hits.items()}, n_prot


def slc6a8_rank(ev, sc):
    m = (ev["uniprot_acc"] == SLC6A8).to_numpy() & ~np.isnan(sc)
    if m.sum() == 0:
        return None
    g = pd.DataFrame({"res": ev["residue"][m], "s": sc[m]}).sort_values("s", ascending=False)
    g["rank"] = np.arange(1, len(g) + 1)
    r193 = g.loc[g["res"] == 193, "rank"]
    best_jon = g[g["res"].isin(JON)]["rank"].min() if len(g[g["res"].isin(JON)]) else np.nan
    return {
        "n": len(g),
        "top_residue": int(g.iloc[0]["res"]),
        "rank_193": int(r193.iloc[0]) if len(r193) else None,
        "best_rank_in_193_201": int(best_jon) if best_jon == best_jon else None,
    }


def main() -> int:
    df, ev, ev_y = load()
    allx = pd.read_csv(HERE / "all_extracellular_targets.tsv", sep="\t", low_memory=False)
    for c in ("ss3", "topology_class", "region_class"):
        allx[c] = allx[c].fillna("").astype(str)
    allx["gene_symbol"] = allx["gene_symbol"].fillna(allx["uniprot_acc"]).astype(str)
    seqs = json.loads((HERE / "sequences_cache.json").read_text())
    df, cys_new = add_disulfide(df, seqs)
    ev, _ = add_disulfide(ev, seqs)
    allx, _ = add_disulfide(allx, seqs)

    # SLC6A8 is not in the evaluation cohort (it has no published site), so it
    # is scored separately by a model fitted on everything.
    slc = allx[allx["uniprot_acc"] == SLC6A8].copy().reset_index(drop=True)
    print(f"SLC6A8 rows available for scoring: {len(slc)}")

    pos = df[df["label"] == 1].copy()
    train = build_global(pos, allx, np.random.default_rng(SEED), RATIO)
    keep = mask_keep(ev, ev_y)
    print(
        f"evaluation: {len(ev)} residues, {int(ev_y.sum())} positives; "
        f"masking +/-{MASK_HALF} drops {int((~keep).sum())} negatives\n"
    )

    rows = []
    for name, (feats, kind) in variants(cys_new).items():
        feats = [f for f in feats if f in train.columns and f in ev.columns]
        sc, fold = fold_eval(train, ev, ev_y, feats, kind)
        folds, topk, _n = metrics(ev, ev_y, sc, fold, keep)
        m = make(kind, feats)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            m.fit(within_protein(train, feats, "z"), train["label"].to_numpy())
            s_slc = m.predict_proba(within_protein(slc, feats, "z"))[:, 1]
        s_slc = smooth(slc, s_slc, 2)
        jon = slc6a8_rank(slc.assign(uniprot_acc=SLC6A8), s_slc)
        rows.append(
            {
                "model": name,
                "n_features": len(feats),
                "estimator": kind,
                "auroc_mean": round(float(folds.mean()), 4),
                "auroc_sd": round(float(folds.std(ddof=1)), 4),
                "folds": ";".join(f"{v:.3f}" for v in folds),
                "top1": round(topk[1], 3),
                "top3": round(topk[3], 3),
                "top10": round(topk[10], 3),
                "slc6a8_top_residue": jon["top_residue"],
                "slc6a8_rank_193": jon["rank_193"],
                "slc6a8_best_rank_in_jon_site": jon["best_rank_in_193_201"],
            }
        )
        print(
            f"  {name:<30} {folds.mean():.4f} +/- {folds.std(ddof=1):.4f}   "
            f"top1 {topk[1]:.2f} top3 {topk[3]:.2f}   "
            f"SLC6A8 best={jon['top_residue']} rank(193)={jon['rank_193']}",
            flush=True,
        )

    for name, col in (("ProDomino", "prodomino"), ("EpicTope", "epictope")):
        sc = ev[col].to_numpy(dtype=float)
        fold = np.zeros(len(ev), int)
        folds, topk, _n = metrics(ev, ev_y, sc, fold, keep & ev[col].notna().to_numpy())
        s_slc = slc[col].to_numpy(dtype=float) if col in slc.columns else np.full(len(slc), np.nan)
        jon = (
            slc6a8_rank(slc.assign(uniprot_acc=SLC6A8), s_slc)
            if not np.all(np.isnan(s_slc))
            else None
        )
        rows.append(
            {
                "model": name,
                "n_features": 1,
                "estimator": "published",
                "auroc_mean": round(float(folds.mean()), 4) if len(folds) else None,
                "auroc_sd": None,
                "folds": ";".join(f"{v:.3f}" for v in folds),
                "top1": round(topk[1], 3),
                "top3": round(topk[3], 3),
                "top10": round(topk[10], 3),
                "slc6a8_top_residue": jon["top_residue"] if jon else None,
                "slc6a8_rank_193": jon["rank_193"] if jon else None,
                "slc6a8_best_rank_in_jon_site": jon["best_rank_in_193_201"] if jon else None,
            }
        )
        print(
            f"  {name:<30} {rows[-1]['auroc_mean']}   top1 {topk[1]:.2f} top3 {topk[3]:.2f}"
            f"   SLC6A8 {'n/a' if not jon else jon['rank_193']}"
        )

    out = pd.DataFrame(rows)
    out.to_csv(HERE / "final_report.tsv", sep="\t", index=False)
    (HERE / "final_report.json").write_text(
        json.dumps(
            {
                "generated_at": datetime.now(UTC).isoformat(),
                "cv": f"GroupKFold by gene, {N_SPLITS} folds",
                "evaluation": f"vs all extracellular residues, negatives within +/-{MASK_HALF} masked",
                "negatives_per_positive": RATIO,
            },
            indent=2,
        )
        + "\n"
    )
    print("\nwrote final_report.tsv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
