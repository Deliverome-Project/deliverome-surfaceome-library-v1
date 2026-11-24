#!/usr/bin/env python3
"""The production tag-site model: 13 features, no dependency on other predictors.

Supersedes `tag_site_model_v3.py` in two ways, both of which change numbers a
reader would otherwise take at face value.

**It drops the four conservation features.** They were never available for the
library: `compute_conservation.py` only ever aligned the 133 proteins carrying a
published site, so 55,453 of the 63,017 extracellular residues of the 515 -- 88%
-- reached the v3 model with conservation imputed, against 14% of its training
rows. `conservation_gap.py` measures what that costs: nothing, 0.8043 -> 0.8036
AUROC. So the features go, because a model that asks for a column its targets do
not have is one imputation away from a silent shift, and the shift lands hardest
on exactly the cross-protein comparison the absolute score exists to make.

**It calibrates against OUT-OF-FOLD published scores.** v3 scored the 176
published sites with a model fitted on those same sites and used their quantiles
as the bar the library had to clear. That bar was inflated: the model recognises
its own training positives. Here each positive is scored by a fold that never
saw its gene, which is the same footing the library is on.

Two models are fitted on the same rows and the same features:

    rank model      features z-scored WITHIN each protein. Ranking residues
                    against their own protein's spread is what lifted AUROC
                    from 0.66 to 0.75; it is the right model for "where in THIS
                    protein", and its score means nothing across proteins.
    absolute model  the same features, raw. Every input is on a scale that
                    means the same thing in any protein -- RSA, pLDDT,
                    disorder, hydropathy, distance to a glycosylation sequon --
                    so its probability IS comparable across proteins, and can
                    be read against the published sites' own distribution.

    uv run --locked python tag_site_model_v4.py --fit
    uv run --locked python tag_site_model_v4.py --design \
        --targets internal_tag_targets_515.tsv --out tag_sites_v4_515.tsv
"""

from __future__ import annotations

import argparse
import json
import subprocess
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

from final_report import MASK_HALF, fold_eval, mask_keep, metrics  # noqa: E402
from independent_model_search import (  # noqa: E402
    CONSERVATION,
    SETS,
    load,
    make,
    smooth,
    within_protein,
)
from independent_model_search import load as load_eval  # noqa: E402
from two_tier_model import GLYC, build_global  # noqa: E402

SEED = 20261014
RATIO = 200
N_SPLITS = 5
FEATS = [f for f in SETS["+ conservation"] + GLYC if f not in CONSERVATION]
KIND = "rf_shallow"
MIN_SEPARATION = 25  # residues; two picks closer than this are one decision
CLEAR_MARGIN = 0.90  # a runner-up region below this fraction of the ceiling is not in contention
MIN_REGION_LEN = 10  # an extracellular region shorter than this is not a separate place to tag
MODEL_Z = HERE / "tag_site_model_v4_rank.joblib"
MODEL_RAW = HERE / "tag_site_model_v4_absolute.joblib"
SPEC = HERE / "tag_site_model_v4.json"


def git_sha() -> str:
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=HERE, text=True).strip()
    except Exception:
        return "unknown"


def training_frame():
    df, _, _ = load()
    allx = pd.read_csv(HERE / "all_extracellular_targets.tsv", sep="\t", low_memory=False)
    for c in ("ss3", "topology_class", "region_class"):
        allx[c] = allx[c].fillna("").astype(str)
    allx["gene_symbol"] = allx["gene_symbol"].fillna(allx["uniprot_acc"]).astype(str)
    pos = df[df["label"] == 1].copy()
    return pos, allx, build_global(pos, allx, np.random.default_rng(SEED), RATIO)


def published_oof(train, allx, feats):
    """Out-of-fold scores for the published cohort, per site AND per protein.

    Two references come out of this, and only one of them is a fair bar.

    `sites` is each published residue scored by a fold that never saw its gene.
    Reading the library's BEST residue against that is not like-for-like: the
    library number is a maximum over several hundred extracellular residues and
    a maximum is pulled upwards by sheer count, while a published site is one
    nominated residue. Comparing the two makes every target look good.

    `ceilings` fixes that by taking the same maximum on both sides -- for each
    published protein, the best out-of-fold score anywhere in its extracellular
    residues. A library protein clearing the median published CEILING is then
    saying something real: this protein offers at least as good a best option as
    half the proteins where somebody did successfully put a tag in.
    """
    y = train["label"].to_numpy()
    groups = train["gene_symbol"].to_numpy()
    oof = np.full(len(train), np.nan)
    # smooth() writes back by label, so the pool needs a clean 0..n-1 index.
    pool = (
        allx[allx["uniprot_acc"].isin(set(train.loc[y == 1, "uniprot_acc"]))]
        .copy()
        .reset_index(drop=True)
    )
    pool_oof = np.full(len(pool), np.nan)
    pool_gene = pool["gene_symbol"].to_numpy()
    for a, b in GroupKFold(n_splits=N_SPLITS).split(train, y, groups):
        if len(np.unique(y[a])) < 2:
            continue
        m = make(KIND, feats)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            m.fit(train.iloc[a], y[a])
            oof[b] = m.predict_proba(train.iloc[b])[:, 1]
            held = set(train.iloc[b]["gene_symbol"]) - set(train.iloc[a]["gene_symbol"])
            mask = np.isin(pool_gene, list(held))
            if mask.any():
                pool_oof[mask] = m.predict_proba(pool[mask])[:, 1]
    pos = (y == 1) & ~np.isnan(oof)
    pool = pool.assign(oof=smooth(pool, pool_oof, 2)).dropna(subset=["oof"])
    ceilings = pool.groupby("uniprot_acc")["oof"].max()
    # Third reference, and the one the ordering rule uses. A REGION's ceiling is
    # what our per-region numbers are, so the bar they are held to has to be the
    # same quantity: the best score inside the region that actually hosts a
    # published site. Held against whole-protein ceilings instead, a single
    # region can essentially never clear the bar and the rule degenerates.
    key = ["uniprot_acc", "region_index"]
    region_max = pool.groupby(key, dropna=False)["oof"].max().rename("region_ceiling")
    sites = train.loc[y == 1, key].drop_duplicates()
    hosting = sites.merge(region_max.reset_index(), on=key, how="inner")
    return (
        np.sort(oof[pos]),
        np.sort(ceilings.to_numpy()),
        np.sort(hosting["region_ceiling"].to_numpy()),
    )


def cv_metrics(train, feats) -> dict:
    """The cross-validated metric, computed rather than quoted.

    These numbers were literals in the spec below, copied from
    `conservation_gap.py`. Nothing recomputed them, so the model and its
    reported performance could drift apart silently -- growing the corpus and
    refitting would leave the recorded AUROC sitting still. Worse, the string
    describing the evaluation named +/-10 masking, which the frame this file
    builds does not do, so the two did not even measure the same thing.

    This reuses `final_report`'s own evaluation so the number stays comparable
    to the one it replaces: score every extracellular residue of every
    published protein out-of-fold, drop negatives within MASK_HALF of a true
    site, then AUROC per fold. The corpus size travels with it, so a stale
    snapshot is visible in the spec rather than inferred.
    """
    _, ev, ev_y = load_eval()
    sc, fold = fold_eval(train, ev, ev_y, feats, KIND)
    folds, topk, n_prot = metrics(ev, ev_y, sc, fold, mask_keep(ev, ev_y))
    return {
        "cv_auroc_mean": round(float(folds.mean()), 4),
        "cv_auroc_sd": round(float(folds.std(ddof=1)), 4),
        "cv_auroc_folds": [round(float(f), 4) for f in folds],
        "cv_top1": round(topk[1], 4),
        "cv_top3": round(topk[3], 4),
        "cv_top10": round(topk[10], 4),
        "cv_n_proteins_scored": n_prot,
        "corpus_n_positives": int((train["label"] == 1).sum()),
        "corpus_n_proteins": int(train.loc[train["label"] == 1, "uniprot_acc"].nunique()),
        "corpus_snapshot": site_features_mtime(),
    }


def site_features_mtime() -> str:
    """When the corpus snapshot was exported, so staleness is on the record."""
    f = HERE / "site_features.tsv"
    return datetime.fromtimestamp(f.stat().st_mtime, UTC).isoformat() if f.exists() else "unknown"


def fit_both():
    pos, allx, train = training_frame()
    feats = [f for f in FEATS if f in train.columns]
    y = train["label"].to_numpy()

    mz = make(KIND, feats)
    mraw = make(KIND, feats)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        mz.fit(within_protein(train, feats, "z"), y)
        mraw.fit(train, y)
        pub_in = np.sort(mraw.predict_proba(pos)[:, 1])
    pub_oof, pub_ceilings, pub_regions = published_oof(train, allx, feats)
    cv = cv_metrics(train, feats)
    import sklearn

    spec = {
        "fitted_at": datetime.now(UTC).isoformat(),
        "git_commit": git_sha(),
        "supersedes": "tag_site_model_v3 (17 features, in-sample calibration)",
        "estimator": "RandomForestClassifier(n_estimators=800, max_depth=4, "
        "min_samples_leaf=12, class_weight='balanced_subsample')",
        "features": feats,
        "n_features": len(feats),
        "dropped_vs_v3": CONSERVATION,
        "why_dropped": "unavailable for 88% of library residues; costs 0.0007 AUROC",
        "rank_model": "features z-scored within each protein; ranks residues within a protein",
        "absolute_model": "raw features; comparable across proteins",
        "calibration": "published sites scored out-of-fold (GroupKFold-5 by gene)",
        "smoothing_half": 2,
        "negatives_per_positive": RATIO,
        "seed": SEED,
        "n_train_rows": len(train),
        "n_positives": int(y.sum()),
        "n_published_scored_oof": len(pub_oof),
        "published_absolute_quantiles_oof": {
            q: round(float(np.quantile(pub_oof, float(q))), 6) for q in ("0.25", "0.5", "0.75")
        },
        "published_absolute_quantiles_in_sample": {
            q: round(float(np.quantile(pub_in, float(q))), 6) for q in ("0.25", "0.5", "0.75")
        },
        "n_published_proteins_with_ceiling": len(pub_ceilings),
        "published_protein_ceiling_quantiles_oof": {
            q: round(float(np.quantile(pub_ceilings, float(q))), 6) for q in ("0.25", "0.5", "0.75")
        },
        "n_published_regions": len(pub_regions),
        "published_region_ceiling_quantiles_oof": {
            q: round(float(np.quantile(pub_regions, float(q))), 6) for q in ("0.25", "0.5", "0.75")
        },
        **cv,
        "cv": f"5-fold GroupKFold by gene, vs all extracellular residues, "
        f"+/-{MASK_HALF} masking, computed here",
        "sklearn_version": sklearn.__version__,
        "numpy_version": np.__version__,
    }
    return mz, mraw, spec, pub_oof, pub_ceilings, pub_regions


def _fill(picks: list, rows, limit: int) -> None:
    """Take residues from `rows` until `limit`, keeping picks a separation apart."""
    for _, r in rows.iterrows():
        if len(picks) >= limit:
            return
        if any(abs(r["residue"] - q["residue"]) < MIN_SEPARATION for q in picks):
            continue
        picks.append(r)


def design(
    mz,
    mraw,
    spec,
    pub_ceilings,
    pub_regions,
    targets_file,
    out,
    extra_tables=(),
    use_mask=True,
    skip_refused=False,
    snorkel_below_pct=25.0,
    min_region_len=MIN_REGION_LEN,
    override_tables=(),
):
    allx = pd.read_csv(HERE / "all_extracellular_targets.tsv", sep="\t", low_memory=False)
    # `override_tables` REPLACES every row of a protein it lists (keep="first" would keep
    # the old rows). Used for the 2026-11-14 proteoform re-export of 15 proteins whose
    # cohort rows were numbered on another isoform; default empty, so nothing else moves.
    for ov in map(Path, override_tables):
        o = pd.read_csv(ov if ov.is_absolute() else HERE / ov, sep="\t", low_memory=False)
        allx = allx[~allx["uniprot_acc"].isin(set(o["uniprot_acc"]))]
        allx = pd.concat([allx, o[allx.columns]], ignore_index=True)
    for c in ("ss3", "topology_class", "region_class"):
        allx[c] = allx[c].fillna("").astype(str)
    allx["gene_symbol"] = allx["gene_symbol"].fillna(allx["uniprot_acc"]).astype(str)
    # Proteoforms the library orders that are not canonical live in their own
    # table: export_site_features.py reaches them only through the additive
    # full-accession pass, and folding them into the cohort table would change
    # a file other results already depend on.
    # `extra_tables` (--extra-extracellular) adds a member exported after these
    # tables were frozen -- SLC67A1, 2026-11-02 -- without rewriting either of them.
    # A protein already present keeps its existing rows (keep="first").
    for extra in (HERE / "l1_member_extracellular.tsv", *map(Path, extra_tables)):
        if not extra.is_absolute():
            extra = HERE / extra
        if not extra.is_file():
            continue
        e = pd.read_csv(extra, sep="\t", low_memory=False)
        for c in ("ss3", "topology_class", "region_class"):
            e[c] = e[c].fillna("").astype(str)
        e["gene_symbol"] = e["gene_symbol"].fillna(e["uniprot_acc"]).astype(str)
        allx = pd.concat([allx, e[allx.columns]], ignore_index=True)
        allx = allx.drop_duplicates(["uniprot_acc", "residue"], keep="first")
    targets = pd.read_csv(HERE / targets_file, sep="\t")
    sym = dict(zip(targets["uniprot_acc"], targets["gene"], strict=False))
    scored = (
        allx[allx["uniprot_acc"].isin(set(targets["uniprot_acc"]))].copy().reset_index(drop=True)
    )
    missing = sorted(set(targets["uniprot_acc"]) - set(scored["uniprot_acc"]))
    if missing:
        print(f"no extracellular rows for {len(missing)}: {missing[:8]}")
    feats = spec["features"]
    trained_genes = set(
        pd.read_csv(HERE / "site_features.tsv", sep="\t", usecols=["uniprot_acc", "label"]).query(
            "label == 1"
        )["uniprot_acc"]
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        scored["rank_score"] = smooth(
            scored, mz.predict_proba(within_protein(scored, feats, "z"))[:, 1], 2
        )
        scored["absolute_score"] = smooth(scored, mraw.predict_proba(scored)[:, 1], 2)
    # The membrane-side mask (membrane_mask.py; decided by Becca 2026-11-10/11,
    # implemented by the ranker-pu-segment session, PR #245). Scores are computed on
    # every v1 O residue as before; only junctions the mask calls extracellular stay
    # eligible, and regions are rebuilt on those, so a v1 segment that merges two
    # helices and a cytoplasmic loop (the anoctamins) splits where it should. A
    # target with no usable mask is refused, never designed on v1 O alone.
    if use_mask:
        import membrane_mask

        scored, refused = membrane_mask.apply(scored)
        if len(refused):
            refused_path = HERE / (Path(out).stem + "_refused.tsv")
            refused.to_csv(refused_path, sep="\t", index=False)
            print(f"\nREFUSED {len(refused)} target(s) with no usable membrane-side mask:")
            print(refused.to_string(index=False))
            print(f"written to {refused_path.name}")
            if not skip_refused:
                raise membrane_mask.MaskMissing(
                    f"{len(refused)} target(s) have no usable membrane-side mask; "
                    "rerun with --skip-refused to design the rest and leave these undesigned"
                )
    else:
        print("\nWARNING: designing WITHOUT the membrane-side mask (--no-mask); not for L1 members")
    # Regions are judged against regions, proteins against proteins. Both
    # references are out-of-fold and both are the same quantity as the number
    # they are compared with.
    q25, med = (float(np.quantile(pub_regions, q)) for q in (0.25, 0.5))
    print(
        f"published REGION ceilings, out-of-fold: 25% {q25:.4f}  50% {med:.4f}  "
        f"(n={len(pub_regions)} regions that host a published site)"
    )

    # pLDDT and RSA ARE the model (permutation importance: 0.044 against 0.017
    # for everything else). Where they are absent the imputer fills the training
    # median, and the forest returns a confident-looking number built on nothing.
    # That is the conservation mistake, so it is marked rather than hidden.
    cover = scored.groupby("uniprot_acc")[["plddt", "rsa"]].apply(
        lambda g: float(min(g["plddt"].notna().mean(), g["rsa"].notna().mean()))
    )
    rows = []
    for acc, g in scored.groupby("uniprot_acc"):
        structure_cover = float(cover.get(acc, 0.0))
        g = g.sort_values("rank_score", ascending=False)
        per_region = (
            g.groupby("region_index", dropna=False)
            .head(1)
            .sort_values("rank_score", ascending=False)
        )
        # Two different questions, two different representatives. WHICH residue
        # to pick inside a region is the rank model's job; whether the region is
        # worth anything at all is the absolute model's.
        region_len = g.groupby("region_index", dropna=False)["residue"].size()
        region_ceiling = g.groupby("region_index", dropna=False)["absolute_score"].max()
        # A six-residue extracellular stub is not somewhere a second construct
        # can go. CD20 has exactly two regions -- a 6-residue loop at 79-84 and
        # the 47-residue loop at 142-188 -- and counting the stub made it look
        # like a protein with an alternative to fall back on. Only regions long
        # enough to host an insertion count as separate options; the floor is
        # empirical, and 5 of the 143 regions that host a published site fall
        # below it, so a short region is still eligible to be PICKED when it is
        # all a protein has.
        # `min_region_len` 1 lets every loop compete on its score (v5; Becca 2026-11-11:
        # a short exposed stretch is fine when the model rates it). The v4 default of 10
        # ignored a better-scoring short loop whenever a longer one existed.
        viable = region_len[region_len >= min_region_len].index
        n_viable = len(viable)
        scoreable = region_ceiling.loc[region_ceiling.index.isin(viable)]
        if scoreable.empty:
            # Nothing is long enough. Score against every region so the protein
            # still gets a ceiling and a pick, but `n_viable` stays 0 -- reusing
            # len(scoreable) here reported the total region count and made a
            # protein with four 6-residue stubs look like it had four places to
            # tag.
            scoreable = region_ceiling
        ordered = scoreable.sort_values(ascending=False)
        ceiling = float(ordered.iloc[0])
        runner_up = float(ordered.iloc[1]) if len(ordered) > 1 else None
        good_regions = int((scoreable >= med).sum())
        ceiling_pct = float((pub_ceilings < ceiling).mean())
        # A second construct buys something only when the model cannot choose.
        # A region that clears the bar with nothing else within 10% of it is not
        # an ambiguous call, and ordering a backup there spends a slot on a site
        # the model already ranks well below the first.
        uncontested = runner_up is None or runner_up < CLEAR_MARGIN * ceiling
        if ceiling < q25:
            n_order, why = (
                3,
                "best region below the 25th percentile of regions hosting a published site",
            )
        elif good_regions == 0:
            n_order, why = 3, "no region clears the median published region"
        elif uncontested:
            n_order, why = 1, "one region clears the median and nothing else comes within 10% of it"
        else:
            # A protein with several good regions is still a protein the model
            # cannot choose for, so the count of good regions does not override
            # a close runner-up. Only the margin decides.
            n_order, why = (
                2,
                "the runner-up region is within 10% of the best; the model cannot choose",
            )
        # Asking for three sites on a protein with one place to put them is not
        # a recommendation. Cap at the regions that can actually take a tag,
        # plus at most one extra residue inside a region long enough to hold two
        # a separation apart.
        # How many well-separated residues a region can hold, capped at two: a
        # second position in the same loop is a real hedge, a third is not a
        # different decision. A `>= 2 * MIN_SEPARATION` test was too blunt --
        # CD20's 47-residue loop holds 162 and 187 comfortably.
        room = sum(
            min(1 + (int(region_len.get(ri, 0)) - 1) // MIN_SEPARATION, 2) for ri in scoreable.index
        )
        if n_order > room:
            n_order = max(room, 1)
            why = f"{why}; only {n_order} place(s) to put a tag on this protein"

        # Fill ONLY from regions long enough to be a separate place, then cap
        # the recommendation at what they delivered.
        #
        # `room` is computed before placement and assumes a region of length L
        # can hold 1 + (L-1)//MIN_SEPARATION sites, which is true only if they
        # are packed to the ends. The first pick is placed by RANK, so it often
        # lands mid-region and leaves nothing a separation away: SLC2A4's region
        # 0 spans 46-81 and room says 2, but 67 was chosen and no residue in the
        # region is 25 away from it. The shortfall then fell through to regions
        # the floor exists to exclude, and nothing re-capped n_order because the
        # pick count had been satisfied.
        #
        # That contradicted this rule's own intent -- a short region should be
        # picked only when it is all a protein has -- in every case it fired:
        # all 19 below-floor picks were on proteins that DID have a viable
        # region. A protein with no viable region at all still falls back,
        # because then the stub genuinely is the only option.
        picks: list = []
        viable_set = set(scoreable.index)
        in_viable = per_region["region_index"].isin(viable_set)
        sources = (
            [per_region[in_viable], g[g["region_index"].isin(viable_set)]]
            if len(viable_set)
            else [per_region, g]
        )
        for source in sources:
            _fill(picks, source, n_order)
        if len(picks) < n_order:
            n_order = len(picks)
            why = f"{why}; only {n_order} site(s) far enough apart exist on this protein"
        picks = picks[:n_order] if n_order >= 1 else picks
        # Which residue inside a region is the rank model's call, but the
        # region's ceiling is the absolute model's, and for about one region in
        # sixteen they disagree by more than 0.05 -- the rank pick is then
        # credited with a ceiling it is not the one to deliver. MS4A1 is the
        # clearest case: the rank model takes 162 at 0.552 while the region's
        # absolute peak is 145 at 0.633. Carry both so the alternative is on
        # the table rather than buried.
        abs_peak = {}
        for ri, sub in g.groupby("region_index", dropna=False):
            top = sub.loc[sub["absolute_score"].idxmax()]
            abs_peak[ri] = (int(top["residue"]), float(top["absolute_score"]))
        for rank, r in enumerate(picks, 1):
            peak_res, peak_abs = abs_peak.get(r["region_index"], (None, float("nan")))
            rows.append(
                {
                    "gene": sym.get(acc, acc),
                    "uniprot_acc": acc,
                    "rank": rank,
                    "insert_after_residue": int(r["residue"]),
                    "region_index": int(r["region_index"]) if pd.notna(r["region_index"]) else None,
                    "v1_region_index": (
                        int(r["v1_region_index"])
                        if use_mask and pd.notna(r.get("v1_region_index"))
                        else None
                    ),
                    "membrane_side_mask": "applied" if use_mask else "NOT APPLIED",
                    "rank_score": round(float(r["rank_score"]), 4),
                    "absolute_score": round(float(r["absolute_score"]), 4),
                    "region_ceiling": round(ceiling, 4),
                    "region_absolute_peak_residue": peak_res,
                    "region_absolute_peak_score": round(peak_abs, 4),
                    "absolute_peak_gap": round(peak_abs - float(r["absolute_score"]), 4),
                    "beats_median_published": bool(r["absolute_score"] >= med),
                    "protein_ceiling_pct_of_published": round(ceiling_pct, 4),
                    "best_region_pct_of_published_regions": round(
                        float((pub_regions < ceiling).mean()), 4
                    ),
                    # Poor everywhere: even the best eligible loop sits below the cutoff
                    # percentile of loops that host a published site, so an internal site
                    # is a weak bet and a terminus or snorkel construct is the alternative
                    # to consider. Decided by Becca 2026-11-11; cutoff set per run.
                    "snorkel_candidate": bool(
                        100 * float((pub_regions < ceiling).mean()) < snorkel_below_pct
                    ),
                    "snorkel_cutoff_pct": snorkel_below_pct,
                    "min_region_len": min_region_len,
                    "n_regions_clearing_median": good_regions,
                    "n_extracellular_regions": int(g["region_index"].nunique()),
                    "n_regions_long_enough_to_tag": n_viable,
                    "this_region_length": int(region_len.get(r["region_index"], 0)),
                    "shares_a_region_with_another_pick": bool(
                        sum(q["region_index"] == r["region_index"] for q in picks) > 1
                    ),
                    "recommended_sites_to_order": n_order,
                    "runner_up_region_ceiling": None if runner_up is None else round(runner_up, 4),
                    "top_region_uncontested": bool(uncontested),
                    "why": why,
                    "in_training_corpus": acc in trained_genes,
                    "structure_feature_coverage": round(structure_cover, 3),
                    "structure_features_imputed": bool(structure_cover < 0.8),
                    # Reported, never fitted. `dist_to_tm` is NOT one of the 13
                    # features, so the model cannot prefer a mid-loop residue to
                    # one pressed against the membrane and will not warn you. The
                    # separation rule makes that worse: forcing a second pick 25
                    # residues from the first pushes it to the loop edge, which
                    # is how MS4A1 ends up recommending 187 at a clearance of 2.
                    # 29 of 167 published sites sit within 5 of a TM, so this is
                    # a caution to weigh rather than a veto.
                    "dist_to_tm": None if pd.isna(r.get("dist_to_tm")) else float(r["dist_to_tm"]),
                    "loop_len": r["loop_len"],
                    "rsa": r["rsa"],
                    "plddt": r["plddt"],
                    "disorder": r["disorder"],
                    "topology_class": r["topology_class"],
                    "tm_helix_count": r["tm_helix_count"],
                }
            )
    out_df = pd.DataFrame(rows)
    out_df.to_csv(HERE / out, sep="\t", index=False)
    first = out_df[out_df["rank"] == 1]
    print(f"\n{first['gene'].nunique()} proteins designed, {len(out_df)} sites")
    print("\nhow many sites to order:")
    print(first.groupby(["recommended_sites_to_order", "why"]).size().to_string())
    print(
        f"\nsnorkel candidates (best eligible loop below the {snorkel_below_pct:g}th percentile "
        f"of published loops): {int(first['snorkel_candidate'].sum())} of {len(first)}"
    )
    print(
        f"\nproteins with at least one region at or above the median published region: "
        f"{int((first['region_ceiling'] >= med).sum())} of {len(first)}"
    )
    return out_df


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--fit", action="store_true")
    ap.add_argument("--design", action="store_true")
    ap.add_argument("--targets", default="internal_tag_targets_515.tsv")
    ap.add_argument("--out", default="tag_sites_v4_515.tsv")
    ap.add_argument(
        "--no-mask",
        action="store_true",
        help="design without the membrane-side mask; only for non-L1 targets the mask does not cover",
    )
    ap.add_argument(
        "--skip-refused",
        action="store_true",
        help="design the targets that have a usable mask and list the rest as refused",
    )
    ap.add_argument(
        "--min-region-len",
        type=int,
        default=MIN_REGION_LEN,
        help=f"loops shorter than this compete only when nothing longer exists "
        f"(default {MIN_REGION_LEN}, the shipped v4 rule; 1 = every loop competes)",
    )
    ap.add_argument(
        "--snorkel-below-pct",
        type=float,
        default=25.0,
        help="flag a protein as a snorkel candidate when its best eligible loop scores below "
        "this percentile of loops hosting a published site (default 25)",
    )
    ap.add_argument(
        "--extra-extracellular",
        action="append",
        default=[],
        help="additional per-residue extracellular table(s) to score from (repeatable)",
    )
    args = ap.parse_args()

    import joblib

    if args.fit or not MODEL_Z.is_file():
        mz, mraw, spec, pub_oof, pub_ceilings, pub_regions = fit_both()
        joblib.dump(mz, MODEL_Z)
        joblib.dump(mraw, MODEL_RAW)
        spec["published_absolute_scores_oof"] = [round(float(v), 6) for v in pub_oof]
        spec["published_protein_ceilings_oof"] = [round(float(v), 6) for v in pub_ceilings]
        spec["published_region_ceilings_oof"] = [round(float(v), 6) for v in pub_regions]
        SPEC.write_text(json.dumps(spec, indent=2) + "\n")
        print(f"fitted on {spec['n_train_rows']} rows, {spec['n_positives']} positives")
        print(f"  features ({spec['n_features']}): {', '.join(spec['features'])}")
        print(f"  published SITES, out-of-fold:    {spec['published_absolute_quantiles_oof']}")
        print(
            f"  published SITES, in-sample (v3): {spec['published_absolute_quantiles_in_sample']}"
        )
        print(
            "  published PROTEIN CEILINGS, OOF: "
            f"{spec['published_protein_ceiling_quantiles_oof']} over "
            f"{spec['n_published_proteins_with_ceiling']} proteins -- the bar the library is held to"
        )
    else:
        mz, mraw = joblib.load(MODEL_Z), joblib.load(MODEL_RAW)
        spec = json.loads(SPEC.read_text())
        pub_ceilings = np.array(spec["published_protein_ceilings_oof"])
        pub_regions = np.array(spec["published_region_ceilings_oof"])

    if args.design:
        design(
            mz,
            mraw,
            spec,
            pub_ceilings,
            pub_regions,
            args.targets,
            args.out,
            args.extra_extracellular,
            use_mask=not args.no_mask,
            skip_refused=args.skip_refused,
            snorkel_below_pct=args.snorkel_below_pct,
            min_region_len=args.min_region_len,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
