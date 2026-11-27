#!/usr/bin/env python3
"""Tag-site model v5: v4 plus upweighted hard negatives, designed on the adopted L1 panel.

What changes from v4, and nothing else:

    training      the same 176 published sites and the same 200-per-site negatives
                  (seed 20261014), plus every hard negative on a published protein
                  that v4's draw left out, with all hard negatives weighted
                  HN_WEIGHT: buried core, disulfide Cys, N-glycan N, sequon N,
                  insertion-intolerant junctions, and D1 ligand contacts
                  (`ranker_hard_negative_layer --ligands`)
    calibration   the published-site references (sites, protein ceilings, region
                  ceilings) are recomputed out of fold with the same weights

Estimator, 13 features, smoothing and the design rules are v4's. Designs go through
`tag_site_model_v4.design`, so the membrane-side mask applies (`membrane_mask.py`).

Why: on one gene -> fold map and three negative draws, this model raises within-10
top-1 from 0.467 to 0.491 on the shipped labels and from 0.462 to 0.494 on the
audit-corrected ones, gaining 3-6 proteins against 1-2 lost on each draw, with
top-3 never down (`refit_corrected_hard_negative.py`).

Targets are the adopted L1 panel (Compara 116 ortholog adoption, PR #283, at
ADOPTED_COMMIT): its `placement == internal` members, mapped by `sequence_sha256`
onto the accession their per-residue features are filed under. A member with no
per-residue features cannot be scored and is listed, not designed.

Decided by Becca: hard negatives (2026-11-03/04), ligand contacts among them
(2026-11-05), shipping the refit model (2026-11-11). Implemented by the
ranker-pu-segment session (Claude Code), PR #245.

    uv run --locked python tag_site_model_v5.py            # fit, then design the adopted panel
"""

from __future__ import annotations

import argparse
import json
import sys
import warnings
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.model_selection import GroupKFold

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[3] / "src"))

import ranker_hard_negative_layer as L  # noqa: E402
import ranker_pu_segment as R  # noqa: E402
import tag_site_model_v4 as V4  # noqa: E402
from independent_model_search import make, smooth, within_protein  # noqa: E402
from two_tier_model import build_global  # noqa: E402

HN_WEIGHT = 5
# Loop-length floor of 4: a loop of 4 or more eligible residues competes on its score;
# 1-3 residue stretches are picked only when nothing longer exists. Decided by Becca
# 2026-11-11 ("set to 4, some min makes sense"), after dropping v4's floor of 10 let
# single-residue stretches win.
MIN_REGION_LEN = 4
ADOPTED_COMMIT = "9a1c2c91"  # on main; read with git show, no branch fetch
ADOPTED_PATH = "library/l1/orf-design/l1_tag_paths.tsv"
MODEL_Z = HERE / "tag_site_model_v5_rank.joblib"
MODEL_RAW = HERE / "tag_site_model_v5_absolute.joblib"
SPEC = HERE / "tag_site_model_v5.json"
TARGETS = HERE / "l1_adopted_targets.tsv"
NO_FEATURES = HERE / "l1_adopted_no_features.tsv"
OUT = "tag_sites_v5_adopted.tsv"
EXTRA = ("slc67a1_extracellular.tsv",)


def training():
    """v4's frame plus upweighted hard negatives, and the per-row weights."""
    L.INTERFACES = "ligands"
    pos, allx_raw, _ev, _y = R.frame()
    allx_hn = L.constraints(allx_raw, R.load_screens())
    base = build_global(pos, allx_raw, np.random.default_rng(V4.SEED), V4.RATIO)
    tr, is_hn, n_added = L.upweight_frame(base, pos, allx_hn)
    w = np.where(is_hn, float(HN_WEIGHT), 1.0)
    return pos, allx_raw, tr, w, int(is_hn.sum()), n_added


def published_oof_weighted(train, w, allx, feats):
    """`tag_site_model_v4.published_oof` with the training weights passed to each fold."""
    y = train["label"].to_numpy()
    groups = train["gene_symbol"].to_numpy()
    oof = np.full(len(train), np.nan)
    pool = (
        allx[allx["uniprot_acc"].isin(set(train.loc[y == 1, "uniprot_acc"]))]
        .copy()
        .reset_index(drop=True)
    )
    pool_oof = np.full(len(pool), np.nan)
    pool_gene = pool["gene_symbol"].to_numpy()
    for a, b in GroupKFold(n_splits=V4.N_SPLITS).split(train, y, groups):
        m = make(V4.KIND, feats)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            m.fit(train.iloc[a], y[a], clf__sample_weight=w[a])
            oof[b] = m.predict_proba(train.iloc[b])[:, 1]
            held = set(train.iloc[b]["gene_symbol"]) - set(train.iloc[a]["gene_symbol"])
            mask = np.isin(pool_gene, list(held))
            if mask.any():
                pool_oof[mask] = m.predict_proba(pool[mask])[:, 1]
    pos = (y == 1) & ~np.isnan(oof)
    pool = pool.assign(oof=smooth(pool, pool_oof, 2)).dropna(subset=["oof"])
    ceilings = pool.groupby("uniprot_acc")["oof"].max()
    key = ["uniprot_acc", "region_index"]
    region_max = pool.groupby(key, dropna=False)["oof"].max().rename("region_ceiling")
    sites = train.loc[y == 1, key].drop_duplicates()
    hosting = sites.merge(region_max.reset_index(), on=key, how="inner")
    return (
        np.sort(oof[pos]),
        np.sort(ceilings.to_numpy()),
        np.sort(hosting["region_ceiling"].to_numpy()),
    )


def fit():
    _pos, allx, tr, w, n_hn, n_added = training()
    feats = [f for f in V4.FEATS if f in tr.columns]
    y = tr["label"].to_numpy()
    mz, mraw = make(V4.KIND, feats), make(V4.KIND, feats)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        mz.fit(within_protein(tr, feats, "z"), y, clf__sample_weight=w)
        mraw.fit(tr, y, clf__sample_weight=w)
    pub_oof, pub_ceilings, pub_regions = published_oof_weighted(tr, w, allx, feats)
    bench = pd.read_csv(HERE / "refit_corrected_hard_negative.tsv", sep="\t")
    bench = bench[
        (bench["training_labels"] == "shipped 176") & (bench["label_set"] == "shipped 176")
    ]
    cv = {
        f"{hn}_{m}": round(float(bench.loc[bench["hard_negatives"] == hn, m].mean()), 4)
        for hn in bench["hard_negatives"].unique()
        for m in ("within10_top1", "within10_top3", "exact_top1", "segment_first")
    }
    import sklearn

    q = ("0.25", "0.5", "0.75")
    spec = {
        "fitted_at": datetime.now(UTC).isoformat(),
        "git_commit": V4.git_sha(),
        "supersedes": "tag_site_model_v4 (same features and estimator; no hard negatives)",
        "estimator": "RandomForestClassifier(n_estimators=800, max_depth=4, "
        "min_samples_leaf=12, class_weight='balanced_subsample'), fitted with sample weights",
        "features": feats,
        "n_features": len(feats),
        "hard_negatives": "buried core (RSA<0.1, pLDDT>=80), UniProt disulfide Cys, UniProt "
        "N-glycan N, sequon N, domain-insertion-intolerant junctions, D1 ligand contacts",
        "hard_negative_weight": HN_WEIGHT,
        "n_hard_negatives_in_frame": n_hn,
        "n_hard_negatives_added_to_v4_frame": n_added,
        "rank_model": "features z-scored within each protein; ranks residues within a protein",
        "absolute_model": "raw features; comparable across proteins",
        "calibration": "published sites scored out-of-fold (GroupKFold-5 by gene), weighted",
        "smoothing_half": 2,
        "negatives_per_positive": V4.RATIO,
        "seed": V4.SEED,
        "n_train_rows": len(tr),
        "n_positives": int(y.sum()),
        "published_absolute_quantiles_oof": {
            k: round(float(np.quantile(pub_oof, float(k))), 6) for k in q
        },
        "published_protein_ceiling_quantiles_oof": {
            k: round(float(np.quantile(pub_ceilings, float(k))), 6) for k in q
        },
        "published_region_ceiling_quantiles_oof": {
            k: round(float(np.quantile(pub_regions, float(k))), 6) for k in q
        },
        "benchmark_mean_over_three_negative_draws": cv,
        "benchmark": "refit_corrected_hard_negative.py: shipped labels, fixed gene folds",
        "decided_by": "Becca: hard negatives 2026-11-03/04, ligand contacts 2026-11-05, "
        "shipping the refit 2026-11-11",
        "implemented_by": "ranker-pu-segment session (Claude Code), PR #245",
        "sklearn_version": sklearn.__version__,
        "numpy_version": np.__version__,
        "published_absolute_scores_oof": [round(float(v), 6) for v in pub_oof],
        "published_protein_ceilings_oof": [round(float(v), 6) for v in pub_ceilings],
        "published_region_ceilings_oof": [round(float(v), 6) for v in pub_regions],
    }
    return mz, mraw, spec, pub_ceilings, pub_regions


def adopted_tag_paths() -> pd.DataFrame:
    from deliverome_analysis.core.provenance import read_pinned_text

    local = HERE / "data" / f"l1_tag_paths_adopted_{ADOPTED_COMMIT}.tsv"
    if not local.is_file():
        repo = HERE.parents[2]
        txt = read_pinned_text(ADOPTED_COMMIT, ADOPTED_PATH, repo=repo)
        local.parent.mkdir(parents=True, exist_ok=True)
        local.write_text(txt)
    return pd.read_csv(local, sep="\t", low_memory=False)


def adopted_targets() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Internal members of the adopted panel, on the accession their features are filed under."""
    tp = adopted_tag_paths()
    internal = tp[tp["placement"] == "internal"]
    rec = pd.read_csv(HERE / "l1_member_reconciliation.tsv", sep="\t")
    tgt = pd.read_csv(HERE / "l1_member_targets.tsv", sep="\t")
    frame = rec[["candidate_id", "sequence_sha256"]].merge(
        tgt[["candidate_id", "uniprot_acc"]], on="candidate_id"
    )
    sha_to_acc = dict(zip(frame["sequence_sha256"], frame["uniprot_acc"], strict=False))
    for extra in EXTRA:  # members exported after the 566 (SLC67A1, PR #281)
        e = pd.read_csv(HERE / extra, sep="\t", usecols=["uniprot_acc"], low_memory=False)
        accs = set(e["uniprot_acc"])
        for _, r in internal[internal["gene_symbol"].isin({"SLC67A1"})].iterrows():
            for a in accs:
                sha_to_acc.setdefault(r["sequence_sha256"], a)
    t = internal.assign(uniprot_acc=internal["sequence_sha256"].map(sha_to_acc))
    t = t.rename(columns={"candidate_id": "member_id", "gene_symbol": "gene"})
    cols = [
        "gene",
        "uniprot_acc",
        "member_id",
        "candidate_set",
        "topology_cohort",
        "sequence_sha256",
    ]
    have, missing = t[t["uniprot_acc"].notna()][cols], t[t["uniprot_acc"].isna()][cols]
    return have, missing


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--refit", action="store_true", help="refit even if the model files exist")
    ap.add_argument(
        "--snorkel-below-pct",
        type=float,
        default=25.0,
        help="flag snorkel candidates below this percentile of published loops (default 25)",
    )
    args = ap.parse_args()
    import joblib

    if args.refit or not MODEL_Z.is_file():
        mz, mraw, spec, pub_ceilings, pub_regions = fit()
        joblib.dump(mz, MODEL_Z)
        joblib.dump(mraw, MODEL_RAW)
        SPEC.write_text(json.dumps(spec, indent=2) + "\n")
        print(
            f"fitted v5 on {spec['n_train_rows']} rows ({spec['n_hard_negatives_in_frame']} hard "
            f"negatives at weight {HN_WEIGHT})"
        )
        print(f"  published REGION ceilings, OOF: {spec['published_region_ceiling_quantiles_oof']}")
    else:
        mz, mraw = joblib.load(MODEL_Z), joblib.load(MODEL_RAW)
        spec = json.loads(SPEC.read_text())
        pub_ceilings = np.array(spec["published_protein_ceilings_oof"])
        pub_regions = np.array(spec["published_region_ceilings_oof"])

    have, missing = adopted_targets()
    have.to_csv(TARGETS, sep="\t", index=False)
    missing.to_csv(NO_FEATURES, sep="\t", index=False)
    print(
        f"adopted panel: {len(have) + len(missing)} internal members; {len(have)} with features, "
        f"{len(missing)} without (listed in {NO_FEATURES.name})"
    )
    V4.design(
        mz,
        mraw,
        spec,
        pub_ceilings,
        pub_regions,
        TARGETS.name,
        OUT,
        EXTRA,
        use_mask=True,
        skip_refused=True,
        snorkel_below_pct=args.snorkel_below_pct,
        min_region_len=MIN_REGION_LEN,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
