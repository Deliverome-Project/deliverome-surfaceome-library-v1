#!/usr/bin/env python3
"""The internal tag-site ranker trained on the confirmed literature sites, scored held out.

Question: trained only on internal tag sites a reader confirmed against the paper (the corpus
audit's gold standard, `library/l1/tag-site-corpus-audit/gold_standard.tsv`), how often does
the ranker put one of its three picks within 10 residues of a confirmed site on a protein it never
saw, overall and by protein class, and which design choices matter?

Every trained model here uses the confirmed sites as positives (`tag_site_model_benchmark.
confirmed_frame`): 200 unlabelled extracellular residues per positive from the same proteins, the
shipped frame's hard negatives, and residues within 10 of any site the literature agent cited, or
inside a region-only label, withheld from both classes. The variants change one choice each:

    confirmed-site ranker     features z-scored within each protein, scores averaged over +/-2
                              residues, hard negatives at weight 5
    global features           raw features, not z-scored
    no smoothing              each residue's own score
    hard negatives x1         hard negatives at weight 1
    no hard negatives         hard-negative rows removed
    + UniProt N-glycans       three more features from UniProt CARBOHYD N-linked annotations

Published predictors (ProDomino, EpicTope, INDELi-E) and single features are scored as stored.

Outputs (tag-site-models/, tables also on S3 under S3_OUT):
    confirmed_ranker_heldout.tsv            per model x split: top-3 / top-1 within 10, median
                                            rank of the best confirmed residue, mean AUROC
    confirmed_ranker_per_protein.tsv        per protein x model x split
    confirmed_ranker_by_class.tsv           per protein class x model (gene split)
    confirmed_ranker_training_by_class.tsv  the training frame by protein class
    confirmed_ranker_feature_blocks.tsv     block ablation of the confirmed-site ranker
    confirmed_ranker_paired.tsv             each model against the confirmed-site ranker
    confirmed_ranker_vs_176.tsv             proteins where the 176-label ranker and the
                                            confirmed-site ranker disagree, with site verdicts
    confirmed_ranker_l1_agreement.tsv       per L1 internal CDS: the decided site against the
                                            confirmed-site ranker's top 3
    confirmed_ranker.json                   run record
Figures (reviewer-package/figures/): tag_site_ranker_heldout.png, tag_site_ranker_by_class.png,
tag_site_ranker_feature_blocks.png

    uv run --locked python library/l1/orf-design/tag_site_ranker_confirmed.py
    uv run --locked python library/l1/orf-design/tag_site_ranker_confirmed.py --figures-only
    uv run --locked python library/l1/orf-design/tag_site_ranker_confirmed.py --counts
"""

from __future__ import annotations

import hashlib
import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import tag_site_model_benchmark as B  # noqa: E402

OUT = B.OUT
FIG = B.FIG
S3_OUT = "sprint4/s4e9-orf-design-fg-panel/tag-site-models/confirmed-ranker/2026-11-22/"
MEMBER_FEATURES = "l1_member_extracellular.tsv"
MEMBER_SHA256 = "376461da576005b72d6bf672547282a109fbced17ca4ce3ac55903e24dea84e4"
PRIMARY = "confirmed-site ranker"
VARIANTS = {
    # name: (norm, smoothing half-width, hard-negative weight or None to drop, extra features)
    PRIMARY: ("z", 2, 5.0, False),
    "global features": ("raw", 2, 5.0, False),
    "no smoothing": ("z", 0, 5.0, False),
    "global features, no smoothing": ("raw", 0, 5.0, False),
    "hard negatives x1": ("z", 2, 1.0, False),
    "no hard negatives": ("z", 2, None, False),
    "+ UniProt N-glycans": ("z", 2, 5.0, True),
}
PUBLISHED = ("ProDomino", "EpicTope", "INDELi-E")
GLYCO = ["uniprot_glycan_here", "dist_to_uniprot_glycan", "n_uniprot_glycan_w10"]
CLASSES = ("single-pass", "multipass, 2-6 TM", "multipass, 7+ TM", "no TM helix")


def protein_class(allx: pd.DataFrame) -> dict[str, str]:
    """Protein class from the DeepTMHMM v1 helix count carried in the feature table."""
    tm = allx.groupby("uniprot_acc")["tm_helix_count"].max()
    out = {}
    for a, n in tm.items():
        n = 0 if n != n else int(n)
        out[a] = (
            CLASSES[3] if n == 0 else CLASSES[0] if n == 1 else CLASSES[1] if n <= 6 else CLASSES[2]
        )
    return out


def add_uniprot_glycans(df: pd.DataFrame) -> pd.DataFrame:
    """UniProt CARBOHYD N-linked positions as three features, on the canonical numbering."""
    u = pd.read_csv(B.MODELS / "uniprot_site_features.tsv", sep="\t")
    covered = set(u["uniprot_acc"])
    g = u[(u["type"] == "CARBOHYD") & u["description"].fillna("").str.contains("N-linked")]
    sites = (
        g.assign(begin=pd.to_numeric(g["begin"], errors="coerce"))
        .dropna(subset=["begin"])
        .groupby("uniprot_acc")["begin"]
        .apply(lambda v: np.array(sorted(v), int))
        .to_dict()
    )
    here, dist, n10 = [], [], []
    for a, r in zip(df["uniprot_acc"], df["residue"], strict=True):
        if a not in covered:
            here.append(np.nan), dist.append(np.nan), n10.append(np.nan)
            continue
        s = sites.get(a, np.array([], int))
        d = np.abs(s - r)
        here.append(float((d == 0).any()))
        dist.append(float(min(d.min(), 200)) if len(d) else 200.0)
        n10.append(float((d <= 10).sum()))
    return df.assign(
        uniprot_glycan_here=here, dist_to_uniprot_glycan=dist, n_uniprot_glycan_w10=n10
    )


def norm(df: pd.DataFrame, feats: list[str], how: str) -> pd.DataFrame:
    return B.zscore(df, feats) if how == "z" else df


def fit(frame, w, feats, cats, how):
    m = B.model(feats, cats)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m.fit(norm(frame, feats, how), frame["label"].to_numpy(), clf__sample_weight=w)
    return m


def oof(frame, w, ev, folds, feats, cats, how, half) -> np.ndarray:
    s = np.full(len(ev), np.nan)
    tf = frame["gene_symbol"].map(folds).to_numpy()
    ef = ev["gene_symbol"].map(folds).to_numpy()
    for i in range(B.N_SPLITS):
        m = fit(frame[tf != i], w[tf != i], feats, cats, how)
        held = ef == i
        if held.any():
            s[held] = m.predict_proba(norm(ev[held], feats, how))[:, 1]
    return B.smooth(ev, s, half) if half else s


def variant_frame(conf, w, hn):
    if hn is None:
        keep = w <= 1.0
        return conf[keep].reset_index(drop=True), np.ones(int(keep.sum()))
    return conf, np.where(w > 1.0, hn, 1.0)


def summary_rows(pp_all: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (m, sc), pp in pp_all.groupby(["model", "split"], sort=False):
        row = {"model": m, "split": sc, **B.summarise(pp)}
        for metric in ("top3_within10", "top1_within10"):
            row[f"{metric}_lo"], row[f"{metric}_hi"] = B.bootstrap(pp, metric)
        rows.append(row)
    return pd.DataFrame(rows)


def by_class(pp_all: pd.DataFrame, cls: dict[str, str]) -> pd.DataFrame:
    g = pp_all[pp_all["split"] == "gene"].assign(protein_class=lambda d: d["uniprot_acc"].map(cls))
    rows = []
    for (m, c), pp in g.groupby(["model", "protein_class"], sort=False):
        row = {"model": m, "protein_class": c, **B.summarise(pp)}
        for metric in ("top3_within10", "top1_within10"):
            row[f"{metric}_lo"], row[f"{metric}_hi"] = B.bootstrap(pp, metric)
        rows.append(row)
    return pd.DataFrame(rows)


def training_by_class(conf: pd.DataFrame, w: np.ndarray, cls: dict[str, str]) -> pd.DataFrame:
    d = conf.assign(protein_class=conf["uniprot_acc"].map(cls), hard=w > 1.0)
    rows = []
    for c in CLASSES:
        x = d[d["protein_class"] == c]
        rows.append(
            {
                "protein_class": c,
                "proteins_with_a_confirmed_site": int(x.loc[x.label == 1, "uniprot_acc"].nunique()),
                "confirmed_sites": int((x.label == 1).sum()),
                "unlabelled_residues": int(((x.label == 0) & ~x.hard).sum()),
                "hard_negatives": int(x.hard.sum()),
            }
        )
    return pd.DataFrame(rows)


def paired_rows(pp_all: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for split in ("gene", "family"):
        g = pp_all[pp_all["split"] == split]
        ref = g[g["model"] == PRIMARY].set_index("uniprot_acc")
        for m in [x for x in g["model"].unique() if x != PRIMARY]:
            o = g[g["model"] == m].set_index("uniprot_acc")
            c = ref.index.intersection(o.index)
            for metric in ("top3_within10", "top1_within10"):
                a, b = ref.loc[c, metric].astype(bool), o.loc[c, metric].astype(bool)
                x, y = int((a & ~b).sum()), int((~a & b).sum())
                rows.append(
                    {
                        "model": m,
                        "split": split,
                        "metric": metric,
                        "proteins": len(c),
                        "confirmed_site_ranker": round(float(a.mean()), 3),
                        "other": round(float(b.mean()), 3),
                        "only_confirmed_site_ranker": x,
                        "only_other": y,
                        "sign_test_p": round(B.sign_test(x, y), 4),
                    }
                )
    return pd.DataFrame(rows)


def l1_agreement(conf, w, allx, feats) -> pd.DataFrame:
    """The decided L1 loop site against the full-fit confirmed-site ranker's top 3 picks.

    The decided sites were placed by the 176-label ranker inside the membrane-side mask; this
    ranks every extracellular residue of the same protein with the confirmed-site ranker (no
    mask), so it measures how far the two rankers agree, not a redesign.
    """
    m = fit(conf, w, feats, B.CATEGORICAL, "z")
    dec = pd.read_csv(B.MODELS / "l1_internal_decisions.tsv", sep="\t")
    dec = dec[dec["decision"].str.startswith("internal_loop")]
    tgt = pd.read_csv(B.MODELS / "l1_adopted_targets.tsv", sep="\t")
    acc = dict(zip(tgt["sequence_sha256"], tgt["uniprot_acc"], strict=False))
    dec = dec.assign(uniprot_acc=dec["sequence_sha256"].map(acc)).dropna(subset=["uniprot_acc"])
    mem = pd.read_csv(B.DATA / MEMBER_FEATURES, sep="\t", low_memory=False)
    for c in B.CATEGORICAL:
        mem[c] = mem[c].fillna("").astype(str)
    mem["residue"] = mem["residue"].astype(int)
    pool = pd.concat(
        [mem, allx[~allx["uniprot_acc"].isin(set(mem["uniprot_acc"]))]], ignore_index=True
    )
    pool = pool[pool["uniprot_acc"].isin(set(dec["uniprot_acc"]))].reset_index(drop=True)
    s = B.smooth(pool, m.predict_proba(B.zscore(pool, feats))[:, 1], B.SMOOTH_HALF)
    v5 = B.smooth(pool, _v5_rank().predict_proba(B.zscore(pool, B.FEATS))[:, 1], B.SMOOTH_HALF)
    pool = pool.assign(s=s, v5=v5)
    rows = []
    for r in dec.itertuples(index=False):
        g = pool[pool["uniprot_acc"] == r.uniprot_acc]
        if g.empty:
            continue
        res, sc = g["residue"].to_numpy(), g["s"].to_numpy()
        top = B.picks(res, sc)
        site = int(r.insert_after_residue)
        order = (-sc).argsort(kind="stable").argsort(kind="stable")
        here = np.flatnonzero(res == site)
        rows.append(
            {
                "member_id": r.member_id,
                "gene": r.gene,
                "uniprot_acc": r.uniprot_acc,
                "decided_site": site,
                "confirmed_ranker_top3": ";".join(str(int(res[i])) for i in top),
                "decided_site_within10_of_top3": bool(np.abs(res[top] - site).min() <= B.FLANK),
                "decided_site_is_top3": bool((res[top] == site).any()),
                "decided_site_rank_pct": round(order[here[0]] / max(len(sc) - 1, 1), 4)
                if len(here)
                else np.nan,
                "spearman_vs_176_ranker": round(
                    float(pd.Series(sc).corr(pd.Series(g["v5"].to_numpy()), method="spearman")), 4
                ),
                "n_residues": len(g),
            }
        )
    return pd.DataFrame(rows)


def _v5_rank():
    import joblib

    return joblib.load(B.MODELS / "tag_site_model_v5_rank.joblib")


def vs_176(tr, ev, folds, gold_sites, gold, pp_conf) -> pd.DataFrame:
    """Where the 176-label ranker and the confirmed-site ranker disagree, and on which sites."""
    s = B.smooth(
        ev, oof(tr, tr["weight"].to_numpy(float), ev, folds, B.FEATS, B.CATEGORICAL, "z", 0)
    )
    pp176 = B.per_protein(ev, s, gold_sites).set_index("uniprot_acc")
    pc = pp_conf.set_index("uniprot_acc")
    verdict = gold.groupby("uniprot_acc")["verdict"].apply(lambda v: ";".join(sorted(set(v))))
    in_tr = gold.groupby("uniprot_acc")["in_training"].max()
    rows = []
    for a in pc.index.intersection(pp176.index):
        for metric in ("top3_within10", "top1_within10"):
            x, y = bool(pp176.loc[a, metric]), bool(pc.loc[a, metric])
            if x != y:
                rows.append(
                    {
                        "uniprot_acc": a,
                        "gene_symbol": pc.loc[a, "gene_symbol"],
                        "metric": metric,
                        "hit_by": "176-label ranker" if x else "confirmed-site ranker",
                        "site_verdicts": verdict.get(a, ""),
                        "in_176": int(in_tr.get(a, 0)),
                    }
                )
    return pd.DataFrame(rows)


def l1_composition() -> dict[str, int]:
    """L1 CDSs on the internal placement by protein class (DeepTMHMM v1 helix count, n_tm)."""
    tp = pd.read_csv(HERE / "l1_tag_paths.tsv", sep="\t", low_memory=False)
    n = tp.loc[tp["placement"] == "internal", "n_tm"].fillna(0).astype(int)
    lab = np.select([n == 0, n == 1, n <= 6], [CLASSES[3], CLASSES[0], CLASSES[1]], CLASSES[2])
    return pd.Series(lab).value_counts().reindex(CLASSES, fill_value=0).astype(int).to_dict()


def composition(pp_all: pd.DataFrame, cls: dict[str, str], comp: dict[str, int]) -> pd.DataFrame:
    """Test proteins (gold standard) against L1 internal CDSs, by protein class."""
    test = pp_all[(pp_all.model == PRIMARY) & (pp_all.split == "gene")]["uniprot_acc"].map(cls)
    t = test.value_counts().reindex(CLASSES, fill_value=0)
    tot = sum(comp.values())
    return pd.DataFrame(
        {
            "protein_class": CLASSES,
            "test_proteins": [int(t[c]) for c in CLASSES],
            "test_fraction": [round(float(t[c] / t.sum()), 4) for c in CLASSES],
            "l1_internal_cds": [comp[c] for c in CLASSES],
            "l1_fraction": [round(comp[c] / tot, 4) for c in CLASSES],
        }
    )


def weighted(
    pp_all: pd.DataFrame, cls: dict[str, str], comp: dict[str, int], n: int = 2000
) -> pd.DataFrame:
    """Hit rates reweighted to the L1 internal class mix; stratified bootstrap over proteins.

    Each class's hit rate is weighted by its share of L1 internal CDSs. A model that scores no
    protein of a class has that class's weight spread over the classes it does score.
    """
    rng = np.random.default_rng(B.SEED)
    rows = []
    d = pp_all.assign(protein_class=pp_all["uniprot_acc"].map(cls))
    for (m, sc), g in d.groupby(["model", "split"], sort=False):
        row = {"model": m, "split": sc, "proteins": len(g)}
        for metric in ("top3_within10", "top1_within10"):
            per = {c: g.loc[g.protein_class == c, metric].astype(float).to_numpy() for c in CLASSES}
            have = [c for c in CLASSES if len(per[c])]
            w = np.array([comp[c] for c in have], float)
            w /= w.sum()
            est = float(sum(wi * per[c].mean() for wi, c in zip(w, have, strict=True)))
            boots = np.zeros(n)
            for c, wi in zip(have, w, strict=True):
                v = per[c]
                boots += wi * rng.choice(v, size=(n, len(v)), replace=True).mean(1)
            row[metric] = round(est, 4)
            row[f"{metric}_lo"], row[f"{metric}_hi"] = (
                round(float(np.quantile(boots, 0.025)), 4),
                round(float(np.quantile(boots, 0.975)), 4),
            )
            row[f"{metric}_unweighted"] = round(float(g[metric].astype(float).mean()), 4)
        ch = {
            c: g.loc[g.protein_class == c, "chance_top1_within10"].astype(float).mean()
            for c in have
        }
        row["chance_top1_within10"] = round(
            float(sum(wi * ch[c] for wi, c in zip(w, have, strict=True))), 4
        )
        rows.append(row)
    return pd.DataFrame(rows)


def _ablation_figure(ablation: pd.DataFrame, colors: list) -> None:
    """Figure 11: the change in each measure when one feature block is removed."""
    import matplotlib.pyplot as plt

    from deliverome_analysis.plotting import save_figure

    fig, axes = plt.subplots(1, 2, figsize=(11, 3.6), sharey=True)
    ab = ablation[ablation["block"] != "none"].set_index("block")
    full = ablation[ablation["block"] == "none"].iloc[0]
    yy = np.arange(len(ab))
    for ax, metric, label in (
        (
            axes[0],
            "top3_within10",
            "Change in fraction with one of\nthe top 3 picks within 10 residues",
        ),
        (axes[1], "mean_auroc", "Change in mean\nper-protein AUROC"),
    ):
        ax.barh(yy, ab[metric] - full[metric], color=colors[0])
        ax.axvline(0, color="grey", lw=0.8)
        ax.set_xlabel(label)
    axes[0].set_yticks(yy)
    axes[0].set_yticklabels([f"without {b}" for b in ab.index])
    axes[0].invert_yaxis()
    fig.tight_layout()
    save_figure(fig, "tag_site_ranker_feature_blocks", output_dir=FIG)
    plt.close(fig)


def figures(summary, cls_tab, ablation, wtd, comp_tab) -> list[Path]:
    import matplotlib.pyplot as plt

    from deliverome_analysis.plotting import save_figure, setup_plotting_style

    setup_plotting_style()
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    FIG.mkdir(parents=True, exist_ok=True)
    out = []

    order = [*VARIANTS, *B.PREDICTORS]
    fig, axes = plt.subplots(1, 2, figsize=(12, 7.0), sharey=True)
    y = np.arange(len(order))
    n = summary[summary["split"] == "gene"].set_index("model")["proteins"]
    for ax, m, label in (
        (
            axes[0],
            "top3_within10",
            "Any of the top 3 picks within 10 residues\nof a confirmed site (fraction of proteins)",
        ),
        (
            axes[1],
            "top1_within10",
            "Top pick within 10 residues\nof a confirmed site (fraction of proteins)",
        ),
    ):
        for k, (scheme, mk) in enumerate((("gene", "o"), ("family", "s"))):
            v = summary[summary["split"] == scheme].set_index("model").reindex(order)
            ax.errorbar(
                v[m],
                y + (k - 0.5) * 0.28,
                xerr=[v[m] - v[f"{m}_lo"], v[f"{m}_hi"] - v[m]],
                fmt=mk,
                color=colors[k],
                capsize=2,
                label=f"held out by {scheme}",
            )
        vw = wtd[wtd["split"] == "gene"].set_index("model").reindex(order)
        ax.errorbar(
            vw[m],
            y + 0.42,
            xerr=[vw[m] - vw[f"{m}_lo"], vw[f"{m}_hi"] - vw[m]],
            fmt="D",
            color=colors[2],
            capsize=2,
            markersize=5,
            label="held out by gene, weighted to L1's class mix",
        )
        ax.axhline(len(VARIANTS) - 0.5, color="grey", lw=0.6)
        ax.set_xlim(0, 1)
        ax.set_xlabel(label)
    chance = summary[(summary["split"] == "gene") & (summary["model"] == PRIMARY)]
    axes[1].axvline(
        chance["chance_top1_within10"].iloc[0], color="grey", ls="--", lw=1, label="random residue"
    )
    axes[0].set_yticks(y)
    axes[0].set_yticklabels([f"{m}  (n={int(n[m])})" for m in order])
    for lab, m in zip(axes[0].get_yticklabels(), order, strict=True):
        if m in PUBLISHED:
            lab.set_fontweight("bold")
    axes[0].invert_yaxis()
    fig.tight_layout()
    h, lab = axes[1].get_legend_handles_labels()
    fig.legend(h, lab, loc="lower center", bbox_to_anchor=(0.6, 1.0), ncol=2, frameon=False)
    save_figure(fig, "tag_site_ranker_heldout", output_dir=FIG)
    out.append(FIG / "tag_site_ranker_heldout.png")
    plt.close(fig)

    shown = [PRIMARY, "ProDomino", "EpicTope", "low pLDDT alone"]
    fig, axes = plt.subplots(1, 3, figsize=(17, 4.6), sharey=True)
    width = 0.2
    x = np.arange(len(CLASSES))
    ct = comp_tab.set_index("protein_class").reindex(CLASSES)
    for k, (col, ncol, lab_) in enumerate(
        (
            ("test_fraction", "test_proteins", "gold-standard test proteins"),
            ("l1_fraction", "l1_internal_cds", "L1 internal CDSs"),
        )
    ):
        axes[0].bar(x + (k - 0.5) * 0.35, ct[col], 0.35, color=["0.6", "0.25"][k], label=lab_)
        for xi, (fr, nn) in enumerate(zip(ct[col], ct[ncol], strict=True)):
            axes[0].text(x[xi] + (k - 0.5) * 0.35, fr + 0.01, f"{int(nn)}", ha="center", fontsize=7)
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(CLASSES, fontsize=8)
    axes[0].set_ylabel("Fraction")
    axes[0].legend(frameon=False, fontsize=8, loc="upper left")
    axes[0].set_title("")
    for ax, m_, label in (
        (axes[1], "top3_within10", "Any of the top 3 picks within 10"),
        (axes[2], "top1_within10", "Top pick within 10"),
    ):
        for k, m in enumerate(shown):
            v = cls_tab[cls_tab["model"] == m].set_index("protein_class").reindex(CLASSES)
            ax.bar(
                x + (k - 1.5) * width,
                v[m_],
                width,
                yerr=[v[m_] - v[f"{m_}_lo"], v[f"{m_}_hi"] - v[m_]],
                color=colors[k],
                capsize=2,
                label=m,
            )
            for xi, nn in enumerate(v["proteins"]):
                if nn == nn:
                    ax.text(
                        x[xi] + (k - 1.5) * width,
                        0.02,
                        f"{int(nn)}",
                        ha="center",
                        fontsize=6,
                        color="white",
                    )
        ax.set_xticks(x)
        ax.set_xticklabels(CLASSES, fontsize=8)
        ax.set_ylim(0, 1)
        ax.set_xlabel(label)
    rand = cls_tab[cls_tab["model"] == PRIMARY].set_index("protein_class").reindex(CLASSES)
    axes[2].scatter(
        x,
        rand["chance_top1_within10"],
        marker="_",
        s=400,
        color="black",
        zorder=3,
        label="random residue",
    )
    fig.tight_layout()
    h, lab = axes[2].get_legend_handles_labels()
    fig.legend(
        h, lab, frameon=False, fontsize=8, ncol=5, loc="lower center", bbox_to_anchor=(0.5, 1.0)
    )
    save_figure(fig, "tag_site_ranker_by_class", output_dir=FIG)
    out.append(FIG / "tag_site_ranker_by_class.png")
    plt.close(fig)

    _ablation_figure(ablation, colors)
    out.append(FIG / "tag_site_ranker_feature_blocks.png")
    return out


def construct_counts() -> pd.DataFrame:
    """CDSs on the internal placement (l1_tag_paths.tsv) by number of built constructs."""
    tp = pd.read_csv(HERE / "l1_tag_paths.tsv", sep="\t", low_memory=False)
    internal = set(tp.loc[tp["placement"] == "internal", "sequence_sha256"])
    a = pd.read_csv(B.MODELS / "l1_internal_arm_decisions.tsv", sep="\t")
    a = a[a["sequence_sha256"].isin(internal)]
    built = a[a["decision"] == "build"]
    per = built.groupby("sequence_sha256").agg(
        constructs=("arm_id", "size"),
        loops=("arm_kind", lambda v: int((v != "snorkel").sum())),
        snorkel=("arm_kind", lambda v: int((v == "snorkel").sum())),
    )
    verdict = (
        a.sort_values("arm_role").drop_duplicates("sequence_sha256").set_index("sequence_sha256")
    )["order_verdict"]
    rows = [{"what": "CDSs on the internal placement", "value": len(internal)}]
    rows.append({"what": "with at least one built construct", "value": len(per)})
    rows.append({"what": "held (nothing built)", "value": len(internal) - len(per)})
    rows.append({"what": "built loop constructs", "value": int(per["loops"].sum())})
    rows.append({"what": "built snorkel constructs", "value": int(per["snorkel"].sum())})
    for k, v in per["constructs"].value_counts().sort_index().items():
        rows.append({"what": f"CDSs with {k} built construct(s)", "value": int(v)})
    for k, v in per["loops"].value_counts().sort_index().items():
        rows.append({"what": f"CDSs with {k} built loop construct(s)", "value": int(v)})
    for k, v in verdict.reindex(per.index).fillna("none recorded").value_counts().items():
        rows.append({"what": f"built CDSs with order verdict {k}", "value": int(v)})
    return pd.DataFrame(rows)


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def class_mix_tables(pp_all, cls):
    """Class composition and L1-weighted hit rates, written beside the other tables."""
    comp = l1_composition()
    wtd, comp_tab = weighted(pp_all, cls, comp), composition(pp_all, cls, comp)
    wtd.to_csv(OUT / "confirmed_ranker_l1_weighted.tsv", sep="\t", index=False)
    comp_tab.to_csv(OUT / "confirmed_ranker_class_composition.tsv", sep="\t", index=False)
    return wtd, comp_tab


def main() -> int:
    paths = {
        k: OUT / f"confirmed_ranker_{k}.tsv"
        for k in (
            "heldout",
            "per_protein",
            "by_class",
            "training_by_class",
            "feature_blocks",
            "paired",
            "vs_176",
            "l1_agreement",
        )
    }
    if "--counts" in sys.argv:
        c = construct_counts()
        c.to_csv(OUT / "confirmed_ranker_construct_counts.tsv", sep="\t", index=False)
        print(c.to_string(index=False))
        return 0
    if "--figures-only" in sys.argv:
        rd = {"sep": "\t"}
        _tr, allx = B.inputs()
        cls_tab = by_class(pd.read_csv(paths["per_protein"], **rd), protein_class(allx))
        cls_tab.to_csv(paths["by_class"], sep="\t", index=False, float_format="%.4g")
        pp_all = pd.read_csv(paths["per_protein"], **rd)
        wtd, comp_tab = class_mix_tables(pp_all, protein_class(allx))
        figures(
            pd.read_csv(paths["heldout"], **rd),
            cls_tab,
            pd.read_csv(paths["feature_blocks"], **rd),
            wtd,
            comp_tab,
        )
        return 0

    tr, allx = B.inputs()
    p = B.DATA / MEMBER_FEATURES
    if not p.is_file():
        from deliverome_analysis.s3 import processed

        processed.download(B.S3_INPUTS + MEMBER_FEATURES, p)
    assert sha256(p) == MEMBER_SHA256, f"{MEMBER_FEATURES} does not match its pinned SHA-256"
    gold, regions, audit = B.truth()
    gold_sites = gold.groupby("uniprot_acc")["residue"].apply(lambda v: np.array(v)).to_dict()
    ev = allx[allx["uniprot_acc"].isin(gold_sites)].reset_index(drop=True)
    have = set(zip(ev["uniprot_acc"], ev["residue"], strict=True))
    gold_in = gold[
        [(a, r) in have for a, r in zip(gold["uniprot_acc"], gold["residue"], strict=True)]
    ]
    gold_sites = gold_in.groupby("uniprot_acc")["residue"].apply(lambda v: np.array(v)).to_dict()
    ev = ev[ev["uniprot_acc"].isin(gold_sites)].reset_index(drop=True)
    conf, conf_w = B.confirmed_frame(tr, allx, gold, regions, audit)
    cls = protein_class(allx)
    print(
        f"confirmed frame: {len(conf)} rows, {int(conf.label.sum())} confirmed sites, "
        f"{int((conf_w > 1).sum())} hard negatives; scored: {len(gold_in)} sites on {len(gold_sites)} proteins"
    )

    conf_g, ev_g = add_uniprot_glycans(conf), add_uniprot_glycans(ev)
    genes = sorted(set(conf["gene_symbol"]) | set(ev["gene_symbol"]))
    pp_rows = []
    for scheme in ("gene", "family"):
        folds = B.fold_map(genes, scheme)
        for name, (how, half, hn, glyco) in VARIANTS.items():
            fr, w = variant_frame(conf_g if glyco else conf, conf_w, hn)
            feats = B.FEATS + GLYCO if glyco else B.FEATS
            e = ev_g if glyco else ev
            pp = B.per_protein(
                ev, oof(fr, w, e, folds, feats, B.CATEGORICAL, how, half), gold_sites
            )
            pp_rows.append(pp.assign(model=name, split=scheme))
            print(scheme, name, B.summarise(pp))
        if scheme == "gene":
            gene_folds = folds
    for name, (col, sign) in B.PREDICTORS.items():
        pp = B.per_protein(ev, B.smooth(ev, ev[col].to_numpy(float) * sign), gold_sites)
        for sc in ("gene", "family"):
            pp_rows.append(pp.assign(model=name, split=sc))
        print("stored", name, B.summarise(pp))
    pp_all = pd.concat(pp_rows, ignore_index=True)
    summary = summary_rows(pp_all)
    cls_tab = by_class(pp_all, cls)
    train_cls = training_by_class(conf, conf_w, cls)

    ab_rows = [
        {
            "block": "none",
            **B.summarise(pp_all[(pp_all.model == PRIMARY) & (pp_all.split == "gene")]),
        }
    ]
    w5 = np.where(conf_w > 1.0, 5.0, 1.0)
    for block, cols in B.BLOCKS.items():
        feats = [f for f in B.FEATS if f not in cols]
        cats = [c for c in B.CATEGORICAL if c not in cols]
        pp = B.per_protein(ev, oof(conf, w5, ev, gene_folds, feats, cats, "z", 2), gold_sites)
        ab_rows.append({"block": block, **B.summarise(pp)})
        print("without", block, B.summarise(pp))
    ablation = pd.DataFrame(ab_rows)

    pp_conf = pp_all[(pp_all.model == PRIMARY) & (pp_all.split == "gene")]
    disagree = vs_176(tr, ev, gene_folds, gold_sites, gold, pp_conf)
    agree = l1_agreement(conf, w5, allx, B.FEATS)

    OUT.mkdir(exist_ok=True)
    for key, df in (
        ("heldout", summary),
        ("per_protein", pp_all),
        ("by_class", cls_tab),
        ("training_by_class", train_cls),
        ("feature_blocks", ablation),
        ("paired", paired_rows(pp_all)),
        ("vs_176", disagree),
        ("l1_agreement", agree),
    ):
        df.to_csv(paths[key], sep="\t", index=False, float_format="%.4g")
    wtd, comp_tab = class_mix_tables(pp_all, cls)
    figs = figures(summary, cls_tab, ablation, wtd, comp_tab)
    meta = {
        "inputs": {
            **{k: v for k, v in B.INPUT_SHA256.items()},
            MEMBER_FEATURES: MEMBER_SHA256,
            "gold_standard.tsv": sha256(B.AUDIT / "gold_standard.tsv"),
            "region_labels.tsv": sha256(B.AUDIT / "region_labels.tsv"),
            "audit.tsv": sha256(B.AUDIT / "audit.tsv"),
            "l1_internal_decisions.tsv": sha256(B.MODELS / "l1_internal_decisions.tsv"),
        },
        "confirmed_frame": {
            "rows": len(conf),
            "confirmed_sites": int(conf["label"].sum()),
            "proteins": int(conf.loc[conf.label == 1, "uniprot_acc"].nunique()),
            "hard_negatives": int((conf_w > 1).sum()),
        },
        "scored": {"sites": len(gold_in), "proteins": len(gold_sites)},
        "l1_agreement": {
            "cds": len(agree),
            "decided_site_within10_of_top3": int(agree["decided_site_within10_of_top3"].sum()),
            "decided_site_is_top3": int(agree["decided_site_is_top3"].sum()),
            "median_spearman_vs_176_ranker": float(agree["spearman_vs_176_ranker"].median()),
        },
        "outputs": {k: str(v.relative_to(B.REPO)) for k, v in paths.items()},
        "figures": [str(f.relative_to(B.REPO)) for f in figs],
        "s3": f"processed://{S3_OUT}",
    }
    (OUT / "confirmed_ranker.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(summary.to_string(index=False))
    print(cls_tab.to_string(index=False))
    print(train_cls.to_string(index=False))
    print(json.dumps(meta["l1_agreement"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
