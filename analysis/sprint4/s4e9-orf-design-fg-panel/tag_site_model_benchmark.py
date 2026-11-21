#!/usr/bin/env python3
"""Internal tag-site models scored on the confirmed literature sites, held out by protein and by family.

Question: how well does each internal tag-site model find a site that a paper confirms,
on proteins it was not trained on?

Truth is the corpus audit's gold standard (`scratch/becca/tag-site-corpus-audit/gold_standard.tsv`):
156 internal tag sites on 113 proteins that a reader confirmed against the cited paper, at the
paper's residue in human canonical numbering. Every model below is scored on the same proteins.

Models (13 numeric per-residue features from the AlphaFold model and sequence, plus three
one-hot categories, the shipped ranker's estimator and features throughout):

    shipped ranker            the L1 internal-site ranker as fitted: the 176 literature labels
                              as collected, 200 unlabelled extracellular residues per label,
                              hard negatives (buried core, disulfide Cys, N-glycan and sequon N,
                              insertion-intolerant junctions, ligand contacts) at weight 5
    no hard-negative weight   the same frame with every row at weight 1
    retrained on confirmed    the same estimator and negatives, with the confirmed sites as
                              positives; residues within 10 of any site the literature agent
                              cited, or inside a region-only label, are neither positive nor
                              negative

Positive-unlabelled (PU) policies, bagged as in Mordelet & Vert (2014) and as `pu_rf` in
s5e13 `compare_internal_policies.py`: 50 forests (the same estimator), each fitted on the
positives against a fresh draw of one unlabelled extracellular residue per positive from the
positives' proteins, scores averaged over the 50. Features are z-scored within each protein over
all its extracellular residues before bagging, since a bag holds too few residues per protein
to z-score in.

    PU, no hard negatives           the 176 literature labels; hard-negative residues are
                                    left in the unlabelled pool (s5e13's pu_rf)
    PU + hard negatives x5 / x1     every bag also holds all 2,141 hard negatives as fixed
                                    negatives at weight 5 (or 1); they leave the unlabelled pool
    PU + hard negatives x5, confirmed
                                    positives are the confirmed sites; hard negatives and the
                                    unlabelled pool withhold the same residues as the confirmed
                                    retrain above (within 10 of a cited site, region labels)

Published predictors, scored as stored (no training): ProDomino, EpicTope and INDELi-E, plus
three single features (disorder, solvent accessibility, low pLDDT) as naive baselines.

Splits: 5-fold GroupKFold by gene, and 5-fold GroupKFold by HGNC gene group so that paralogs
(KCNK, KCNJ, GABA-A subunits, ...) are never on both sides of a fold. A trained model scores a
protein only from a fit that did not see its gene (or its family).

Metrics, per protein over all its extracellular residues (the residues the ranker chooses
among), scores smoothed over +/-2 residues as in production:

    top-1 within 10   the highest-scoring residue lies within 10 of a confirmed site
    top-3 within 10   one of the three highest-scoring residues at least 11 apart does
    rank percentile   where the best-scoring confirmed residue falls (0 = top)
    AUROC             confirmed residues against residues more than 10 from any confirmed site

Feature importance is a block ablation of the shipped ranker under the gene split.

    uv run --locked python analysis/sprint4/s4e9-orf-design-fg-panel/tag_site_model_benchmark.py
    uv run --locked python analysis/sprint4/s4e9-orf-design-fg-panel/tag_site_model_benchmark.py --figures-only
"""

from __future__ import annotations

import hashlib
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
MODELS = REPO / "analysis/sprint5/s5e13-tag-site-benchmarking"
AUDIT = REPO / "scratch/becca/tag-site-corpus-audit"
OUT = HERE / "tag-site-models"
DATA = OUT / "data"
FIG = HERE / "reviewer-package" / "figures"
S3_INPUTS = "sprint5/s5e13-tag-site-benchmarking/v5/2026-11-14/"
S3_OUT = "sprint4/s4e9-orf-design-fg-panel/tag-site-models/2026-11-21/"
INPUT_SHA256 = {
    "v5_training_table.tsv": "e7d948688cf47f1cbf1ad3a6ff39a4e256f1822a6091c9b3174122568e10d23c",
    "all_extracellular_targets.tsv": "d3dcdac545b1bfe322df2965f5510c45c68eee780fb4b8c08628ea21462ffefe",
}

SEED = 20261014  # negative draw and bootstrap
RF_SEED = 20261013  # the shipped forest's random_state (train_site_models.SEED)
N_SPLITS = 5
FLANK = 10
SMOOTH_HALF = 2
FEATS = [
    "plddt",
    "plddt_w5",
    "rsa",
    "rsa_w5",
    "disorder",
    "disorder_w5",
    "disorder_w15",
    "kd_hydropathy_w7",
    "pro_gly_frac_w7",
    "cys_in_w10",
    "glyc_sequon_here",
    "dist_to_glyc",
    "n_glyc_w10",
]
CATEGORICAL = ["ss3", "topology_class", "region_class"]
BLOCKS = {
    "structure model (pLDDT, RSA)": ["plddt", "plddt_w5", "rsa", "rsa_w5"],
    "disorder (3 windows)": ["disorder", "disorder_w5", "disorder_w15"],
    "sequence (hydropathy, Pro/Gly, Cys)": ["kd_hydropathy_w7", "pro_gly_frac_w7", "cys_in_w10"],
    "N-glycosylation (3)": ["glyc_sequon_here", "dist_to_glyc", "n_glyc_w10"],
    "categories (SS, topology, region)": CATEGORICAL,
}
TRAINED = (
    "L1 site ranker",
    "ranker, hard negatives not up-weighted",
    "ranker retrained on confirmed sites",
    "PU, no hard negatives",
    "PU + hard negatives at weight 5",
    "PU + hard negatives at weight 1",
    "PU + hard negatives at weight 5, confirmed sites",
)
PU_ARMS = TRAINED[3:]
PU_BAGS = 50  # s5e13 ranker_pu_segment.PU_BAGS
PU_BAG_RATIO = 1  # unlabelled residues per positive in each bag (ranker_pu_segment.PU_BAG_RATIO)
HN_WEIGHT = 5.0  # the shipped ranker's hard-negative weight
#: Each PU arm and published predictor is also compared, paired, against these.
REFERENCES = (TRAINED[0], "PU, no hard negatives", "ProDomino", "EpicTope")
PREDICTORS = {
    "ProDomino": ("prodomino", 1),
    "EpicTope": ("epictope", 1),
    "INDELi-E": ("indeli_e", 1),
    "disorder alone": ("disorder_w5", 1),
    "solvent accessibility alone": ("rsa_w5", 1),
    "low pLDDT alone": ("plddt_w5", -1),
}
GABA_FAMILY = "Gamma-aminobutyric acid type A receptor subunits"


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def inputs() -> tuple[pd.DataFrame, pd.DataFrame]:
    """The shipped ranker's pinned training table and extracellular feature table."""
    DATA.mkdir(parents=True, exist_ok=True)
    for name, want in INPUT_SHA256.items():
        p = DATA / name
        if not p.is_file():
            from deliverome_analysis.s3 import processed

            processed.download(S3_INPUTS + name, p)
        got = sha256(p)
        assert got == want, f"{name}: SHA-256 {got} != pinned {want}"
    allx = pd.read_csv(DATA / "all_extracellular_targets.tsv", sep="\t", low_memory=False)
    for c in CATEGORICAL:
        allx[c] = allx[c].fillna("").astype(str)
    allx["gene_symbol"] = allx["gene_symbol"].fillna(allx["uniprot_acc"]).astype(str)
    allx["residue"] = allx["residue"].astype(int)
    allx = allx.drop_duplicates(["uniprot_acc", "residue"]).reset_index(drop=True)
    tr = pd.read_csv(DATA / "v5_training_table.tsv", sep="\t")
    n = len(tr)
    tr = tr.merge(
        allx[["uniprot_acc", "residue", *CATEGORICAL]], on=["uniprot_acc", "residue"], how="left"
    )
    assert len(tr) == n and tr["ss3"].notna().all(), "training rows missing from feature table"
    return tr, allx


def model(feats: list[str], cats: list[str]) -> Pipeline:
    """The shipped ranker's estimator (`independent_model_search.make('rf_shallow')`)."""
    parts = [("num", SimpleImputer(strategy="median"), feats)]
    if cats:
        parts.append(("cat", OneHotEncoder(handle_unknown="ignore", min_frequency=5), cats))
    clf = RandomForestClassifier(
        n_estimators=800,
        max_depth=4,
        min_samples_leaf=12,
        class_weight="balanced_subsample",
        random_state=RF_SEED,
        n_jobs=-1,
    )
    return Pipeline([("pre", ColumnTransformer(parts)), ("clf", clf)])


def zscore(df: pd.DataFrame, feats: list[str]) -> pd.DataFrame:
    """Features z-scored within each protein: the rank model's input."""
    out = df.copy()
    g = out.groupby("uniprot_acc")
    for f in feats:
        mu, sd = g[f].transform("mean"), g[f].transform("std").replace(0, np.nan)
        out[f] = (out[f] - mu) / sd
    return out


def smooth(df: pd.DataFrame, s: np.ndarray, half: int = SMOOTH_HALF) -> np.ndarray:
    out = s.copy()
    tmp = df[["uniprot_acc", "residue"]].assign(s=s)
    for _, g in tmp.groupby("uniprot_acc"):
        g = g.sort_values("residue")
        pad = np.pad(g["s"].to_numpy(), half, mode="edge")
        out[g.index.to_numpy()] = np.convolve(pad, np.ones(2 * half + 1) / (2 * half + 1), "valid")
    return out


def truth() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    gold = pd.read_csv(AUDIT / "gold_standard.tsv", sep="\t")
    gold = gold.rename(columns={"insert_after_residue": "residue"})
    regions = pd.read_csv(AUDIT / "region_labels.tsv", sep="\t")
    audit = pd.read_csv(AUDIT / "audit.tsv", sep="\t")
    return gold, regions, audit


def near(df: pd.DataFrame, sites: dict[str, np.ndarray], half: int = FLANK) -> np.ndarray:
    out = np.zeros(len(df), bool)
    for i, (a, r) in enumerate(zip(df["uniprot_acc"], df["residue"], strict=True)):
        s = sites.get(a)
        if s is not None and len(s) and np.abs(s - r).min() <= half:
            out[i] = True
    return out


def in_windows(df: pd.DataFrame, regions: pd.DataFrame) -> np.ndarray:
    w = regions.dropna(subset=["window_start", "window_end"])
    out = np.zeros(len(df), bool)
    for a, lo, hi in zip(w["uniprot_acc"], w["window_start"], w["window_end"], strict=True):
        out |= (df["uniprot_acc"].to_numpy() == a) & df["residue"].between(lo, hi).to_numpy()
    return out


def confirmed_frame(tr, allx, gold, regions, audit) -> tuple[pd.DataFrame, np.ndarray]:
    """Shipped frame with confirmed sites as positives; cited and region residues withheld."""
    rng = np.random.default_rng(SEED)
    gold_sites = gold.groupby("uniprot_acc")["residue"].apply(lambda v: np.array(v)).to_dict()
    cited = audit.dropna(subset=["insert_after_residue"])
    cited_sites = (
        cited.groupby("uniprot_acc")["insert_after_residue"]
        .apply(lambda v: np.array(v, int))
        .to_dict()
    )
    neg = tr[tr["label"] == 0]
    drop = near(neg, gold_sites) | near(neg, cited_sites) | in_windows(neg, regions)
    neg = neg[~drop]
    key = allx.set_index(["uniprot_acc", "residue"])
    pos_idx = [k for k in zip(gold["uniprot_acc"], gold["residue"], strict=True) if k in key.index]
    pos = key.loc[pos_idx].reset_index().assign(label=1, weight=1.0)
    # Proteins with a confirmed site but no rows in the shipped frame (the six sites the agent
    # had filed under another gene): unlabelled negatives drawn as the shipped frame drew them.
    new = sorted(set(pos["uniprot_acc"]) - set(tr["uniprot_acc"]))
    extra = []
    for a in new:
        pool = allx[allx["uniprot_acc"] == a]
        pool = pool[~near(pool, gold_sites) & ~in_windows(pool, regions)]
        k = min(len(pool), 200 * int((pos["uniprot_acc"] == a).sum()))
        extra.append(pool.iloc[rng.choice(len(pool), size=k, replace=False)].assign(label=0))
    cols = ["uniprot_acc", "gene_symbol", "residue", *FEATS, *CATEGORICAL, "label", "weight"]
    extra_df = pd.concat(extra).assign(weight=1.0) if extra else pd.DataFrame(columns=cols)
    frame = pd.concat([pos[cols], neg[cols], extra_df[cols]], ignore_index=True)
    frame = frame.astype({"label": int, "weight": float, "residue": int})
    return frame, frame["weight"].to_numpy(float)


def families(genes: list[str]) -> dict[str, str]:
    fam = (
        pd.read_csv(MODELS / "training_gene_families.tsv", sep="\t")
        .set_index("gene_symbol")["family"]
        .to_dict()
    )
    out = {}
    for g in genes:
        if g in fam:
            out[g] = fam[g]
        elif g.startswith("GABR"):
            out[g] = GABA_FAMILY
        else:
            out[g] = g
    return out


def fold_map(genes: list[str], scheme: str) -> dict[str, int]:
    genes = sorted(set(genes))
    groups = [families(genes)[g] for g in genes] if scheme == "family" else genes
    out = {}
    for i, (_a, b) in enumerate(GroupKFold(n_splits=N_SPLITS).split(genes, groups=groups)):
        for j in b:
            out[genes[j]] = i
    return out


def fit(frame, w, feats, cats):
    m = model(feats, cats)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m.fit(zscore(frame, feats), frame["label"].to_numpy(), clf__sample_weight=w)
    return m


def oof_scores(frame, w, ev, folds, feats=FEATS, cats=CATEGORICAL) -> np.ndarray:
    """Rank-model score of every evaluation residue from a fit that never saw its group."""
    s = np.full(len(ev), np.nan)
    tf = frame["gene_symbol"].map(folds).to_numpy()
    ef = ev["gene_symbol"].map(folds).to_numpy()
    for i in range(N_SPLITS):
        m = fit(frame[tf != i], w[tf != i], feats, cats)
        held = ef == i
        if held.any():
            s[held] = m.predict_proba(zscore(ev[held], feats))[:, 1]
    return smooth(ev, s)


def pu_sets(
    allx: pd.DataFrame,
    pos_keys: list[tuple[str, int]],
    hn_keys: set[tuple[str, int]],
    withheld: np.ndarray | None = None,
) -> pd.DataFrame:
    """Rows for a PU fit: role p (positive), h (hard negative) or u (unlabelled pool).

    Features are z-scored within each protein over all its extracellular residues. The
    unlabelled pool is drawn from the positives' proteins only; positive rows repeat as often
    as the label list repeats them. `withheld` (aligned with `allx`) removes residues from the
    unlabelled pool. Pass no hard negatives to leave them in the unlabelled pool.
    """
    pos_prots = {a for a, _ in pos_keys}
    keep = allx["uniprot_acc"].isin(pos_prots | {a for a, _ in hn_keys}).to_numpy()
    wh = np.zeros(int(keep.sum()), bool) if withheld is None else np.asarray(withheld)[keep]
    full = zscore(allx[keep], FEATS).reset_index(drop=True)
    key = list(zip(full["uniprot_acc"], full["residue"], strict=True))
    pset = set(pos_keys)
    is_p = np.array([k in pset for k in key])
    is_h = np.array([k in hn_keys for k in key]) & ~is_p
    idx = pd.Series(np.arange(len(full)), index=pd.MultiIndex.from_tuples(key))
    pos = full.iloc[idx.loc[pos_keys].to_numpy()].assign(role="p")
    hn = full[is_h].assign(role="h")
    in_pos = full["uniprot_acc"].isin(pos_prots).to_numpy()
    unl = full[~is_p & ~is_h & ~wh & in_pos].assign(role="u")
    assert len(hn) == len(hn_keys - pset), "hard negatives missing from the feature table"
    return pd.concat([pos, hn, unl], ignore_index=True)


def oof_pu(sets, ev, folds, hn_weight, seed=SEED) -> np.ndarray:
    """Bagged PU, out of fold: each bag is the training positives, a fresh unlabelled draw of
    PU_BAG_RATIO per positive and (hn_weight > 0) every training hard negative at hn_weight."""
    rng = np.random.default_rng(seed)
    s = np.full(len(ev), np.nan)
    tf = sets["gene_symbol"].map(folds).to_numpy()
    ef = ev["gene_symbol"].map(folds).to_numpy()
    assert not (np.isnan(tf).any() or np.isnan(ef).any()), "gene without a fold"
    role = sets["role"].to_numpy()
    for i in range(N_SPLITS):
        held = ef == i
        if not held.any():
            continue
        pi = np.where((tf != i) & (role == "p"))[0]
        hi = np.where((tf != i) & (role == "h"))[0] if hn_weight else np.array([], int)
        ui = np.where((tf != i) & (role == "u"))[0]
        n = min(len(ui), PU_BAG_RATIO * len(pi))
        target = zscore(ev[held], FEATS)
        acc = np.zeros(int(held.sum()))
        for _ in range(PU_BAGS):
            rows = np.concatenate([pi, hi, rng.choice(ui, size=n, replace=False)])
            y = np.r_[np.ones(len(pi)), np.zeros(len(hi) + n)].astype(int)
            w = np.r_[np.ones(len(pi)), np.full(len(hi), float(hn_weight or 1)), np.ones(n)]
            m = model(FEATS, CATEGORICAL)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                m.fit(sets.iloc[rows], y, clf__sample_weight=w)
            acc += m.predict_proba(target)[:, 1]
        s[held] = acc / PU_BAGS
    return smooth(ev, s)


def pu_frames(tr, allx, conf, gold, regions, audit) -> dict[str, tuple[pd.DataFrame, float]]:
    """The PU arms' row sets and hard-negative weights (0 = no hard negatives)."""
    lit_pos = list(
        zip(*(tr.loc[tr["label"] == 1, c] for c in ("uniprot_acc", "residue")), strict=True)
    )
    lit_hn = tr[(tr["label"] == 0) & (tr["weight"] == HN_WEIGHT)]
    lit_hn = set(zip(lit_hn["uniprot_acc"], lit_hn["residue"], strict=True))
    conf_pos = conf[conf["label"] == 1]
    conf_pos = list(zip(conf_pos["uniprot_acc"], conf_pos["residue"], strict=True))
    conf_hn = conf[(conf["label"] == 0) & (conf["weight"] == HN_WEIGHT)]
    conf_hn = set(zip(conf_hn["uniprot_acc"], conf_hn["residue"], strict=True))
    gold_sites = gold.groupby("uniprot_acc")["residue"].apply(lambda v: np.array(v)).to_dict()
    cited = audit.dropna(subset=["insert_after_residue"])
    cited_sites = (
        cited.groupby("uniprot_acc")["insert_after_residue"]
        .apply(lambda v: np.array(v, int))
        .to_dict()
    )
    withheld = near(allx, gold_sites) | near(allx, cited_sites) | in_windows(allx, regions)
    lit_with_hn = pu_sets(allx, lit_pos, lit_hn)
    out = {
        "PU, no hard negatives": (pu_sets(allx, lit_pos, set()), 0.0),
        "PU + hard negatives at weight 5": (lit_with_hn, HN_WEIGHT),
        "PU + hard negatives at weight 1": (lit_with_hn, 1.0),
        "PU + hard negatives at weight 5, confirmed sites": (
            pu_sets(allx, conf_pos, conf_hn, withheld),
            HN_WEIGHT,
        ),
    }
    assert tuple(out) == PU_ARMS
    for name, (sets, hw) in out.items():
        n = sets["role"].value_counts().to_dict()
        print(f"{name}: {n} (hard-negative weight {hw:g})")
    return out


def picks(res: np.ndarray, s: np.ndarray, k: int = 3, sep: int = FLANK + 1) -> list[int]:
    out: list[int] = []
    for i in np.argsort(-s, kind="stable"):
        if all(abs(res[i] - res[j]) >= sep for j in out):
            out.append(i)
        if len(out) == k:
            break
    return out


def per_protein(ev, scores, gold_sites) -> pd.DataFrame:
    rows = []
    tmp = ev[["uniprot_acc", "gene_symbol", "residue"]].assign(s=scores)
    for a, g in tmp.groupby("uniprot_acc"):
        g = g.dropna(subset=["s"])
        if g.empty:
            continue
        res, s = g["residue"].to_numpy(), g["s"].to_numpy()
        sites = gold_sites[a]
        d = np.abs(res[:, None] - sites[None, :]).min(1)
        exact, far = d == 0, d > FLANK
        top = picks(res, s)
        order = (-s).argsort(kind="stable").argsort(kind="stable")
        rank_pct = order[exact].min() / max(len(s) - 1, 1) if exact.any() else np.nan
        auc = (
            roc_auc_score(exact[exact | far], s[exact | far])
            if exact.any() and far.any()
            else np.nan
        )
        rows.append(
            {
                "uniprot_acc": a,
                "gene_symbol": g["gene_symbol"].iloc[0],
                "n_residues": len(s),
                "top1_within10": bool(d[top[0]] <= FLANK),
                "top3_within10": bool((d[top] <= FLANK).any()),
                "rank_pct": rank_pct,
                "auroc": auc,
                "chance_top1_within10": float((d <= FLANK).mean()),
            }
        )
    return pd.DataFrame(rows)


def summarise(pp: pd.DataFrame) -> dict:
    return {
        "proteins": len(pp),
        "top1_within10": round(float(pp["top1_within10"].mean()), 3),
        "top3_within10": round(float(pp["top3_within10"].mean()), 3),
        "median_rank_pct": round(float(pp["rank_pct"].median()), 3),
        "mean_auroc": round(float(pp["auroc"].mean()), 3),
        "chance_top1_within10": round(float(pp["chance_top1_within10"].mean()), 3),
    }


def check_shipped(tr: pd.DataFrame) -> float:
    """Refit on the whole pinned frame and compare with the committed shipped rank model."""
    import joblib

    m = fit(tr, tr["weight"].to_numpy(float), FEATS, CATEGORICAL)
    saved = joblib.load(MODELS / "tag_site_model_v5_rank.joblib")
    z = zscore(tr, FEATS)
    d = float(np.abs(m.predict_proba(z)[:, 1] - saved.predict_proba(z)[:, 1]).max())
    assert d < 1e-9, f"refit differs from the committed shipped model by {d}"
    return d


def figures(pp_all: pd.DataFrame, summary: pd.DataFrame, ablation: pd.DataFrame) -> list[Path]:
    import matplotlib.pyplot as plt

    from deliverome_analysis.plotting import save_figure, setup_plotting_style

    setup_plotting_style()
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    FIG.mkdir(parents=True, exist_ok=True)
    out = []

    # Figure A: top-1 and top-3 within 10, per model, gene vs family split, one protein = one unit.
    order = [*TRAINED, *PREDICTORS]
    sub = summary[summary["model"].isin(order)]
    fig, axes = plt.subplots(1, 2, figsize=(11, 5.0), sharey=True)
    y = np.arange(len(order))
    for ax, metric, label in (
        (
            axes[0],
            "top1_within10",
            "Top pick within 10 residues\n(fraction of proteins)",
        ),
        (
            axes[1],
            "top3_within10",
            "Any of top 3 picks within 10\n(fraction of proteins)",
        ),
    ):
        for k, (scheme, mk) in enumerate((("gene", "o"), ("family", "s"))):
            v = sub[sub["split"] == scheme].set_index("model").reindex(order)
            ax.errorbar(
                v[metric],
                y + (k - 0.5) * 0.25,
                xerr=[v[metric] - v[f"{metric}_lo"], v[f"{metric}_hi"] - v[metric]],
                fmt=mk,
                color=colors[k],
                label=f"held out by {scheme}",
                capsize=2,
            )
        chance = sub[(sub["split"] == "gene") & (sub["model"] == TRAINED[0])]
        if metric == "top1_within10":
            ax.axvline(
                chance["chance_top1_within10"].iloc[0],
                color="grey",
                ls="--",
                lw=1,
                label="random residue",
            )
        ax.set_xlabel(label)
        ax.set_xlim(0, 1)
    axes[0].set_yticks(y)
    axes[0].set_yticklabels(
        [
            f"{m}  (n={int(sub[(sub.model == m) & (sub.split == 'gene')]['proteins'].iloc[0])})"
            for m in order
        ]
    )
    axes[0].invert_yaxis()
    axes[0].axhline(len(TRAINED) - 0.5, color="grey", lw=0.6)
    fig.legend(
        *axes[0].get_legend_handles_labels(),
        loc="upper center",
        ncol=3,
        frameon=False,
        bbox_to_anchor=(0.5, 1.04),
    )
    fig.tight_layout()
    save_figure(fig, "tag_site_models_heldout", output_dir=FIG)
    out.append(FIG / "tag_site_models_heldout.png")
    plt.close(fig)

    # Figure B: cumulative rank percentile of the best confirmed residue, family split.
    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    shown = [
        "L1 site ranker",
        "ranker retrained on confirmed sites",
        "ProDomino",
        "EpicTope",
        "low pLDDT alone",
    ]
    for k, m in enumerate(shown):
        v = pp_all[(pp_all["model"] == m) & (pp_all["split"] == "family")]["rank_pct"].dropna()
        x = np.sort(v.to_numpy())
        ax.step(
            np.r_[0, x],
            np.r_[0, np.arange(1, len(x) + 1) / len(x)],
            where="post",
            color=colors[k % len(colors)],
            ls="-" if m in TRAINED else "--",
            label=f"{m} (n={len(x)})",
        )
    ax.plot([0, 1], [0, 1], color="grey", lw=0.8, ls="--", label="random")
    ax.set_xlabel(
        "Rank of the best-scoring confirmed residue\n(fraction of the protein's extracellular residues)"
    )
    ax.set_ylabel("Fraction of proteins")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.legend(frameon=False, fontsize=7, loc="lower right")
    fig.tight_layout()
    save_figure(fig, "tag_site_models_rank_cdf", output_dir=FIG)
    out.append(FIG / "tag_site_models_rank_cdf.png")
    plt.close(fig)

    # Figure C: block ablation of the shipped ranker, gene split, change against the full model.
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.6), sharey=True)
    ab = ablation[ablation["block"] != "none"].set_index("block")
    full = ablation[ablation["block"] == "none"].iloc[0]
    yy = np.arange(len(ab))
    for ax, metric, label in (
        (axes[0], "mean_auroc", "Change in mean\nper-protein AUROC"),
        (
            axes[1],
            "top3_within10",
            "Change in fraction with one of\nthe top 3 picks within 10 residues",
        ),
    ):
        ax.barh(yy, ab[metric] - full[metric], color=colors[0])
        ax.axvline(0, color="grey", lw=0.8)
        ax.set_xlabel(label)
    axes[0].set_yticks(yy)
    axes[0].set_yticklabels([f"without {b}" for b in ab.index])
    axes[0].invert_yaxis()
    fig.tight_layout()
    save_figure(fig, "tag_site_models_feature_blocks", output_dir=FIG)
    out.append(FIG / "tag_site_models_feature_blocks.png")
    plt.close(fig)
    return out


def bootstrap(pp: pd.DataFrame, metric: str, n: int = 2000) -> tuple[float, float]:
    rng = np.random.default_rng(SEED)
    v = pp[metric].astype(float).dropna().to_numpy()
    b = rng.choice(v, size=(n, len(v)), replace=True).mean(1)
    return float(np.quantile(b, 0.025)), float(np.quantile(b, 0.975))


def sign_test(a: int, b: int) -> float:
    """Exact two-sided sign test on discordant pairs."""
    from math import comb

    n, k = a + b, min(a, b)
    return min(1.0, 2 * sum(comb(n, i) for i in range(k + 1)) / 2**n) if n else 1.0


def paired(pp_all: pd.DataFrame, split: str = "gene") -> pd.DataFrame:
    """Each model against each reference on the proteins both score.

    Against the shipped ranker every model is compared; against the other REFERENCES the
    shipped ranker and the PU arms are. Top-1 and top-3 count the proteins only one of the two
    gets within 10; AUROC counts the proteins where one has the higher per-protein AUROC (ties
    and proteins without an AUROC dropped). Each is an exact two-sided sign test.
    """
    g = pp_all[pp_all["split"] == split]
    rows = []
    for ref in REFERENCES:
        r = g[g["model"] == ref].set_index("uniprot_acc")
        others = g["model"].unique() if ref == TRAINED[0] else [TRAINED[0], *PU_ARMS]
        for m in [x for x in others if x != ref]:
            o = g[g["model"] == m].set_index("uniprot_acc")
            c = r.index.intersection(o.index)
            for metric in ("top1_within10", "top3_within10", "auroc"):
                if metric == "auroc":
                    a, b = r.loc[c, metric].astype(float), o.loc[c, metric].astype(float)
                    ok = a.notna() & b.notna()
                    a, b = a[ok], b[ok]
                    only_ref, only_other = int((a > b).sum()), int((b > a).sum())
                else:
                    a, b = r.loc[c, metric].astype(bool), o.loc[c, metric].astype(bool)
                    only_ref, only_other = int((a & ~b).sum()), int((~a & b).sum())
                rows.append(
                    {
                        "reference": ref,
                        "model": m,
                        "split": split,
                        "metric": metric,
                        "proteins": len(a),
                        "reference_value": round(float(a.mean()), 3),
                        "other": round(float(b.mean()), 3),
                        "only_reference_better": only_ref,
                        "only_other_better": only_other,
                        "sign_test_p": round(sign_test(only_ref, only_other), 4),
                    }
                )
    return pd.DataFrame(rows)


def main() -> int:
    import sys

    if "--figures-only" in sys.argv:  # redraw from the tables a full run wrote
        rd = {"sep": "\t"}
        figures(
            pd.read_csv(OUT / "tag_site_models_heldout_per_protein.tsv", **rd),
            pd.read_csv(OUT / "tag_site_models_heldout.tsv", **rd),
            pd.read_csv(OUT / "tag_site_models_feature_blocks.tsv", **rd),
        )
        return 0
    tr, allx = inputs()
    d = check_shipped(tr)
    print(
        f"refit on the pinned frame matches the committed shipped rank model (max |diff| {d:.1e})"
    )
    gold, regions, audit = truth()
    gold_sites = gold.groupby("uniprot_acc")["residue"].apply(lambda v: np.array(v)).to_dict()
    ev = allx[allx["uniprot_acc"].isin(gold_sites)].reset_index(drop=True)
    have = set(zip(ev["uniprot_acc"], ev["residue"], strict=True))
    gold_in = gold[
        [(a, r) in have for a, r in zip(gold["uniprot_acc"], gold["residue"], strict=True)]
    ]
    gold_sites = gold_in.groupby("uniprot_acc")["residue"].apply(lambda v: np.array(v)).to_dict()
    ev = ev[ev["uniprot_acc"].isin(gold_sites)].reset_index(drop=True)
    print(
        f"gold standard: {len(gold)} sites on {gold['uniprot_acc'].nunique()} proteins; "
        f"{len(gold_in)} sites on {len(gold_sites)} proteins are extracellular residues in the "
        "feature table and are scored"
    )

    conf, conf_w = confirmed_frame(tr, allx, gold, regions, audit)
    frames = {
        "L1 site ranker": (tr, tr["weight"].to_numpy(float)),
        "ranker, hard negatives not up-weighted": (tr, np.ones(len(tr))),
        "ranker retrained on confirmed sites": (conf, conf_w),
    }
    pu = pu_frames(tr, allx, conf, gold, regions, audit)
    genes = sorted(
        set(tr["gene_symbol"])
        | set(conf["gene_symbol"])
        | set(ev["gene_symbol"])
        | {g for sets, _ in pu.values() for g in sets["gene_symbol"]}
    )
    pp_rows, summary_rows = [], []
    for scheme in ("gene", "family"):
        folds = fold_map(genes, scheme)
        for name, (frame, w) in frames.items():
            pp = per_protein(ev, oof_scores(frame, w, ev, folds), gold_sites)
            pp_rows.append(pp.assign(model=name, split=scheme))
            print(scheme, name, summarise(pp))
        for name, (sets, hw) in pu.items():
            pp = per_protein(ev, oof_pu(sets, ev, folds, hw), gold_sites)
            pp_rows.append(pp.assign(model=name, split=scheme))
            print(scheme, name, summarise(pp), flush=True)
        if scheme == "gene":  # stored predictor scores do not depend on the split
            for name, (col, sign) in PREDICTORS.items():
                s = ev[col].to_numpy(float) * sign
                pp = per_protein(ev, smooth(ev, s), gold_sites)
                for sc in ("gene", "family"):
                    pp_rows.append(pp.assign(model=name, split=sc))
                print("stored", name, summarise(pp))
    pp_all = pd.concat(pp_rows, ignore_index=True)
    for (m, sc), pp in pp_all.groupby(["model", "split"], sort=False):
        row = {"model": m, "split": sc, **summarise(pp)}
        for metric in ("top1_within10", "top3_within10", "auroc"):
            row[f"{metric}_lo"], row[f"{metric}_hi"] = bootstrap(pp, metric)
        summary_rows.append(row)
    summary = pd.DataFrame(summary_rows)

    folds = fold_map(genes, "gene")
    frame, w = frames["L1 site ranker"]
    ab_rows = [
        {
            "block": "none",
            **summarise(pp_all[(pp_all.model == TRAINED[0]) & (pp_all.split == "gene")]),
        }
    ]
    for block, cols in BLOCKS.items():
        feats = [f for f in FEATS if f not in cols]
        cats = [c for c in CATEGORICAL if c not in cols]
        pp = per_protein(ev, oof_scores(frame, w, ev, folds, feats, cats), gold_sites)
        ab_rows.append({"block": block, **summarise(pp)})
        print("without", block, summarise(pp))
    ablation = pd.DataFrame(ab_rows)

    OUT.mkdir(exist_ok=True)
    paths = {
        "summary": OUT / "tag_site_models_heldout.tsv",
        "per_protein": OUT / "tag_site_models_heldout_per_protein.tsv",
        "ablation": OUT / "tag_site_models_feature_blocks.tsv",
        "paired": OUT / "tag_site_models_paired.tsv",
    }
    pd.concat([paired(pp_all, "gene"), paired(pp_all, "family")]).to_csv(
        paths["paired"], sep="\t", index=False
    )
    summary.to_csv(paths["summary"], sep="\t", index=False, float_format="%.4g")
    pp_all.to_csv(paths["per_protein"], sep="\t", index=False, float_format="%.4g")
    ablation.to_csv(paths["ablation"], sep="\t", index=False, float_format="%.4g")
    figs = figures(pp_all, summary, ablation)
    meta = {
        "inputs": {
            **{
                k: {"sha256": v, "s3": f"processed://{S3_INPUTS}{k}"}
                for k, v in INPUT_SHA256.items()
            },
            "gold_standard.tsv": sha256(AUDIT / "gold_standard.tsv"),
            "region_labels.tsv": sha256(AUDIT / "region_labels.tsv"),
            "audit.tsv": sha256(AUDIT / "audit.tsv"),
        },
        "refit_vs_committed_shipped_rank_model_max_abs_diff": d,
        "gold_sites_scored": len(gold_in),
        "gold_proteins_scored": len(gold_sites),
        "confirmed_frame": {
            "rows": len(conf),
            "positives": int(conf["label"].sum()),
            "proteins_with_positive": int(conf.loc[conf.label == 1, "uniprot_acc"].nunique()),
        },
        "pu": {
            "bags": PU_BAGS,
            "unlabelled_per_positive_per_bag": PU_BAG_RATIO,
            "seed": SEED,
            "arms": {
                name: {
                    "hard_negative_weight": hw,
                    "rows_by_role": {
                        {"p": "positive", "h": "hard_negative", "u": "unlabelled_pool"}[k]: int(v)
                        for k, v in sets["role"].value_counts().items()
                    },
                }
                for name, (sets, hw) in pu.items()
            },
        },
        "splits": {
            "gene": "GroupKFold-5 by gene symbol",
            "family": "GroupKFold-5 by HGNC gene group (training_gene_families.tsv; GABA-A "
            "subunits added to their group)",
        },
        "outputs": {k: str(v.relative_to(REPO)) for k, v in paths.items()},
        "figures": [str(p.relative_to(REPO)) for p in figs if p is not None],
    }
    (OUT / "tag_site_models_heldout.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(summary.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
