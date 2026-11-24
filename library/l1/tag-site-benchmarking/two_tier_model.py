#!/usr/bin/env python3
"""A two-tier model: choose the loop, then choose the residue in it.

The design decision is two-stage and the model was not. Picking a tag site means
first deciding **which extracellular region** takes the cassette, then **where in
that region** to put it. Those are different questions with different features:
whether a loop can host a 30-40 aa insert is about the loop (its length, how
exposed and flexible it is overall, how conserved); where in the loop is about
the residue.

Tier 1 scores regions. Tier 2 scores residues within a region, with every
feature expressed relative to that region. A residue's final score is the
product, so a good residue in a hopeless loop cannot win.

Also tested here, in the final configuration rather than the early one:
**disulfide geometry** (is the residue inside a cysteine pair, how far to the
nearest cysteine), **glycosylation** sequons, **ligand and PPI contacts**, and
**homo-oligomer state**. Those blocks were tried once, early, before
within-protein normalisation and before training on the evaluation
distribution, and dismissed on that basis — which was premature.

Scored vs all extracellular residues with negatives within +/-10 of a published
site masked out: a residue five away shares nearly every windowed feature with
the site and would very likely take the tag itself, so counting it as a
negative charges the model for being right.

    uv run --locked python two_tier_model.py
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import Pipeline

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[3] / "src"))

from independent_model_search import SETS, load, make, smooth, within_protein  # noqa: E402
from train_site_models import N_SPLITS  # noqa: E402

SEED = 20261014
OURS = SETS["+ conservation"]
GLYC = ["glyc_sequon_here", "dist_to_glyc", "n_glyc_w10"]
CONTACTS = [
    "contact_here",
    "dist_to_contact",
    "n_contact_w10",
    "contact_ecd_here",
    "dist_to_contact_ecd",
    "dist_to_pocket",
]
HOMOMER = ["homomer_stoichiometry", "homomer_higher_order"]
GEOMETRY = ["loop_len", "dist_to_tm", "dist_to_terminus", "rel_position"]
MASK_HALF = 10
RATIO = 200


def add_disulfide(d: pd.DataFrame, seqs: dict[str, str]) -> tuple[pd.DataFrame, list[str]]:
    """Where the residue sits relative to the cysteines of its own region.

    `cys_in_w10` already counts cysteines nearby. What it cannot say is whether
    an insertion would lengthen a loop that two cysteines hold shut, which is
    the thing that actually constrains a tag.
    """
    out = d.copy()
    n_cys, inside, d_left, d_right = [], [], [], []
    for acc, res in zip(out["uniprot_acc"], out["residue"], strict=False):
        s = seqs.get(acc)
        if not s:
            n_cys.append(np.nan)
            inside.append(np.nan)
            d_left.append(np.nan)
            d_right.append(np.nan)
            continue
        cys = [i + 1 for i, ch in enumerate(s) if ch == "C"]
        before = [c for c in cys if c < res]
        after = [c for c in cys if c > res]
        n_cys.append(len(cys))
        inside.append(1.0 if before and after else 0.0)
        d_left.append(res - before[-1] if before else np.nan)
        d_right.append(after[0] - res if after else np.nan)
    out["n_cys_protein"] = n_cys
    out["between_cys"] = inside
    out["dist_cys_left"] = d_left
    out["dist_cys_right"] = d_right
    out["cys_pair_span"] = np.array(d_left, float) + np.array(d_right, float)
    return out, ["n_cys_protein", "between_cys", "dist_cys_left", "dist_cys_right", "cys_pair_span"]


def region_table(d: pd.DataFrame, feats: list[str]) -> pd.DataFrame:
    """One row per extracellular region, aggregated from its residues."""
    g = d.groupby(["uniprot_acc", "region_index"], dropna=False)
    agg = g.agg(
        n_residues=("residue", "size"),
        start=("residue", "min"),
        end=("residue", "max"),
        label=("label", "max") if "label" in d.columns else ("residue", "size"),
        gene_symbol=("gene_symbol", "first"),
        topology_class=("topology_class", "first"),
    ).reset_index()
    for f in feats:
        if f not in d.columns:
            continue
        agg[f"{f}_max"] = g[f].max().to_numpy()
        agg[f"{f}_mean"] = g[f].mean().to_numpy()
    agg["region_len"] = agg["end"] - agg["start"] + 1
    return agg


def masked_auc(ev, ev_y, sc, half=MASK_HALF):
    keep = np.ones(len(ev), bool)
    pos_by_acc: dict[str, list[int]] = {}
    for a, r, lab in zip(ev["uniprot_acc"], ev["residue"], ev_y, strict=False):
        if lab == 1:
            pos_by_acc.setdefault(a, []).append(int(r))
    for i, (a, r, lab) in enumerate(zip(ev["uniprot_acc"], ev["residue"], ev_y, strict=False)):
        if lab == 1:
            continue
        near = pos_by_acc.get(a)
        if near and min(abs(int(r) - p) for p in near) <= half:
            keep[i] = False
    ok = keep & ~np.isnan(sc)
    return float(roc_auc_score(ev_y[ok], sc[ok]))


def build_global(pos, allx, rng, ratio):
    keys = set(zip(pos["uniprot_acc"], pos["residue"], strict=False))
    pool = allx[allx["uniprot_acc"].isin(set(pos["uniprot_acc"]))].copy()
    pool = pool[~pool.apply(lambda r: (r["uniprot_acc"], r["residue"]) in keys, axis=1)]
    n = min(len(pool), ratio * len(pos))
    neg = pool.iloc[rng.choice(len(pool), size=n, replace=False)].copy()
    neg["label"] = 0
    cols = [c for c in pos.columns if c in neg.columns]
    return pd.concat([pos[cols], neg[cols]], ignore_index=True)


def residue_oof(train, ev, feats, kind="hgb"):
    y = train["label"].to_numpy()
    groups = train["gene_symbol"].to_numpy()
    tr = within_protein(train, feats, "z")
    ed = within_protein(ev, feats, "z")
    gg = ev["gene_symbol"].to_numpy()
    sc = np.full(len(ev), np.nan)
    for a, b in GroupKFold(n_splits=min(N_SPLITS, len(np.unique(groups)))).split(tr, y, groups):
        if len(np.unique(y[a])) < 2:
            continue
        m = make(kind, feats)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            m.fit(tr.iloc[a], y[a])
        mask = np.isin(gg, list(set(tr.iloc[b]["gene_symbol"])))
        if mask.any():
            sc[mask] = m.predict_proba(ed[mask])[:, 1]
    return smooth(ev, sc, 2)


def region_oof(train_regions, ev_regions, feats, fold_of_gene=None):
    """Tier 1: probability that a region hosts a site at all.

    `fold_of_gene` (gene -> fold) pins the split to one fixed elsewhere, so a
    region model can be compared with a residue model on identical folds. Left
    as None, the regions are split by their own GroupKFold, as before.
    """
    y = train_regions["label"].to_numpy().astype(int)
    groups = train_regions["gene_symbol"].to_numpy()
    cols = [c for c in feats if c in train_regions.columns and c in ev_regions.columns]
    out = np.full(len(ev_regions), np.nan)
    key_ev = list(zip(ev_regions["uniprot_acc"], ev_regions["region_index"], strict=False))
    if fold_of_gene is None:
        splits = GroupKFold(n_splits=min(N_SPLITS, len(np.unique(groups)))).split(
            train_regions, y, groups
        )
    else:
        # a gene absent from the map (-1) is training-only: in every fold, scored in none
        f = np.array([fold_of_gene.get(g, -1) for g in groups])
        splits = [(np.where(f != k)[0], np.where(f == k)[0]) for k in sorted(set(f) - {-1})]
    for a, b in splits:
        if len(np.unique(y[a])) < 2:
            continue
        pipe = Pipeline(
            [
                ("imp", SimpleImputer(strategy="median")),
                (
                    "clf",
                    RandomForestClassifier(
                        n_estimators=600,
                        max_depth=4,
                        min_samples_leaf=6,
                        class_weight="balanced_subsample",
                        random_state=SEED,
                        n_jobs=-1,
                    ),
                ),
            ]
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            pipe.fit(train_regions.iloc[a][cols], y[a])
        held = set(train_regions.iloc[b]["gene_symbol"])
        mask = np.array([g in held for g in ev_regions["gene_symbol"]])
        if mask.any():
            out[mask] = pipe.predict_proba(ev_regions[mask][cols])[:, 1]
    return pd.Series(out, index=pd.MultiIndex.from_tuples(key_ev)), cols


def main() -> int:
    df, ev, ev_y = load()
    allx = pd.read_csv(HERE / "all_extracellular_targets.tsv", sep="\t", low_memory=False)
    for c in ("ss3", "topology_class", "region_class"):
        allx[c] = allx[c].fillna("").astype(str)
    allx["gene_symbol"] = allx["gene_symbol"].fillna(allx["uniprot_acc"]).astype(str)

    # sequences for the disulfide geometry
    seqs: dict[str, str] = {}
    cache = HERE / "sequences_cache.json"
    if cache.exists():
        import json

        seqs = json.loads(cache.read_text())
    df, cys_new = add_disulfide(df, seqs)
    ev, _ = add_disulfide(ev, seqs)
    allx, _ = add_disulfide(allx, seqs)
    have_cys = bool(seqs)
    print(
        f"disulfide features: {'built from cached sequences' if have_cys else 'NO SEQUENCES — skipped'}"
    )

    pos = df[df["label"] == 1].copy()
    rng = np.random.default_rng(SEED)
    train = build_global(pos, allx, rng, RATIO)

    blocks = {
        "ours (14)": OURS,
        "+ glycosylation": OURS + GLYC,
        "+ contacts / PPI": OURS + CONTACTS,
        "+ homomer": OURS + HOMOMER,
        "+ geometry": OURS + GEOMETRY,
        "+ everything": OURS + GLYC + CONTACTS + HOMOMER + GEOMETRY,
    }
    if have_cys:
        blocks["+ disulfide"] = OURS + cys_new
        blocks["+ everything"] = OURS + GLYC + CONTACTS + HOMOMER + GEOMETRY + cys_new

    print(f"\nsingle-tier, masked +/-{MASK_HALF}:")
    results = {}
    for name, feats in blocks.items():
        feats = [f for f in feats if f in train.columns and f in ev.columns]
        sc = residue_oof(train, ev, feats)
        a = masked_auc(ev, ev_y, sc)
        results[name] = (a, feats, sc)
        print(f"  {name:<20} {len(feats):>2}f   AUROC {a:.4f}", flush=True)

    best_name = max(results, key=lambda k: results[k][0])
    best_feats = results[best_name][1]
    print(f"\nbest single-tier: {best_name} ({results[best_name][0]:.4f})")

    print("\ntwo-tier:")
    # Region features have to be computed over ALL residues of a region, not
    # over the sampled training subset. At 1:200 a training region holds a
    # handful of sampled residues while an evaluation region holds ~75, so a
    # "max" means something different on each side and the model inverts --
    # which is exactly what the first version did, scoring 0.30 on a broadcast
    # of a region classifier that is 0.85 on its own.
    pos_keys = set(zip(pos["uniprot_acc"], pos["residue"], strict=False))
    full = allx[allx["uniprot_acc"].isin(set(pos["uniprot_acc"]))].copy()
    full["label"] = [
        1 if (a, r) in pos_keys else 0
        for a, r in zip(full["uniprot_acc"], full["residue"], strict=False)
    ]
    tr_reg = region_table(full, best_feats)
    ev_reg = region_table(ev.assign(label=ev_y), best_feats)
    reg_feats = [c for c in tr_reg.columns if c.endswith(("_max", "_mean"))] + [
        "n_residues",
        "region_len",
    ]
    reg_scores, used = region_oof(tr_reg, ev_reg, reg_feats)
    print(
        f"  tier 1: {len(used)} region features, "
        f"{int(tr_reg['label'].sum())} positive regions of {len(tr_reg)}"
    )
    ok = reg_scores.notna()
    if ok.any():
        ra = roc_auc_score(ev_reg["label"].to_numpy()[ok.to_numpy()], reg_scores[ok].to_numpy())
        print(f"  tier 1 region AUROC {ra:.4f}")
    key = list(zip(ev["uniprot_acc"], ev["region_index"], strict=False))
    rmap = reg_scores.to_dict()
    rvec = np.array([rmap.get(k, np.nan) for k in key])
    res_sc = results[best_name][2]

    # A product is the wrong combination, and the reason is worth stating.
    # Tier 1 is good at regions (0.85) but broadcasting it to residues inverts:
    # positive residues average 0.70 while negative residues average 0.86,
    # because the big exposed regions tier 1 likes are exactly the ones holding
    # most of the negatives. Region-level AUROC weights each region once;
    # residue-level weights each residue once, and the size weighting flips the
    # sign. So tier 1 belongs as a GATE -- keep the top regions, rank residues
    # inside them -- not as a multiplier.
    def gate(k: int) -> np.ndarray:
        out = res_sc.copy()
        keep = np.zeros(len(ev), bool)
        top = {
            acc: set(
                ev_reg[ev_reg["uniprot_acc"] == acc]
                .assign(
                    s=[
                        rmap.get((acc, ri), np.nan)
                        for ri in ev_reg[ev_reg["uniprot_acc"] == acc]["region_index"]
                    ]
                )
                .nlargest(k, "s")["region_index"]
            )
            for acc in ev["uniprot_acc"].unique()
        }
        for i, (acc, ri) in enumerate(zip(ev["uniprot_acc"], ev["region_index"], strict=False)):
            keep[i] = ri in top.get(acc, set())
        out[~keep] = out[~keep] * 0.1
        return out

    for label, combo in (
        ("product", rvec * res_sc),
        ("region only", rvec),
        ("residue only", res_sc),
        ("gate: top-1 region", gate(1)),
        ("gate: top-2 regions", gate(2)),
        ("gate: top-3 regions", gate(3)),
    ):
        a = masked_auc(ev, ev_y, combo)
        print(f"  {label:<26} AUROC {a:.4f}")
        results[f"two-tier: {label}"] = (a, best_feats, combo)

    rows = [
        {"variant": k, "n_features": len(v[1]), "auroc_masked10": round(v[0], 4)}
        for k, v in results.items()
    ]
    pd.DataFrame(rows).sort_values("auroc_masked10", ascending=False).to_csv(
        HERE / "two_tier_model.tsv", sep="\t", index=False
    )
    best = max(results, key=lambda k: results[k][0])
    print(f"\nBEST OVERALL: {best}  {results[best][0]:.4f}")
    print("wrote two_tier_model.tsv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
