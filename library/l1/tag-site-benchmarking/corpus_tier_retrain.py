#!/usr/bin/env python3
"""Does training on only the better-evidenced positives improve the v4 ranker?

The 176 published sites carry evidence tiers that the production model ignores:
`evidence_type` (exact site, same loop, or a different kind of insertion),
`position_evidence` (residue pinned by the paper or inferred), `validation_level`
(surface shown, function measured) and the agent's `confidence`. An earlier
stage of this experiment (`train_site_models.py`, matched negatives, AUROC
only) found that models fitted on the better-validated subsets scored higher.
That was never tested under the production protocol, where it matters, nor on
top-1, which is the number that decides an order.

Here each subset of the 176 -- defined exactly as `train_site_models.SUBSETS`
defines them, plus `confidence == "high"` -- trains the rank model under the
production protocol unchanged: the 13 features of `tag_site_model_v4.py`,
within-protein z-scoring, 200 negatives per positive drawn from the proteins
that carry a kept positive (a nested subsample of the all-176 negatives, same
seed), +/-10 masking, +/-2 smoothing, and the SAME five GroupKFold-by-gene
splits for every variant. A dropped positive is neither a positive nor a
negative: it leaves the training frame, and if its protein has no kept
positive the protein's negatives leave with it. The fold assignment is fixed
from the all-176 frame and asserted identical across variants, and the
all-176 variant is asserted to reproduce `final_report.fold_eval` exactly, so
every per-fold difference is paired.

Two evaluation designs, because they answer different questions:

    fixed   every variant is scored on all 176 positives of the held-out genes,
            so the comparison measures the quality of the training labels only.
    own     every variant is scored on its own subset's positives, as the
            earlier benchmark was. A lower-tier positive on a kept protein is
            then neither a positive nor a negative: it and its +/-10 flank are
            withheld from AUROC and masked top-1, and it earns no credit under
            within-10. The all-176 model is scored on the same relabelled rows,
            so the paired delta still compares two models on one evaluation.

Three recalls follow `structure_global_ablation.py`: AUROC and `top1` are the
spec's (negatives within +/-10 of a published site removed before ranking);
`top1_within10` keeps them and credits a first pick within 10 residues of a
published site, the same decision under MIN_SEPARATION 25. Strata follow
`region_class_strata.py` (`region_class` x TM-helix count, cells under 8
proteins dropped). The five EndoNB sites are scored with each gene held out,
as in `endonb_scored.py`.

    uv run --locked python corpus_tier_retrain.py
"""

from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[3] / "src"))

from endonb_scored import leave_gene_out, rank_of  # noqa: E402
from final_report import fold_eval, mask_keep, metrics  # noqa: E402
from independent_model_search import load, make, smooth, within_protein  # noqa: E402
from region_class_strata import topo  # noqa: E402
from tag_site_model_v4 import FEATS, KIND, N_SPLITS, RATIO, SEED, training_frame  # noqa: E402
from train_site_models import SUBSETS as TRAIN_SUBSETS  # noqa: E402

WITHIN = 10  # residues; a pick this close to a published site is the same decision
MIN_PROTEINS = 8  # a stratum smaller than this is not reported, as in region_class_strata
BASE = "all internal"

# The tiers, in the order the README tables use. Predicates come from
# train_site_models so the subsets are the ones the earlier stage scored.
SUBSETS: dict[str, tuple] = {
    BASE: TRAIN_SUBSETS["all internal"],
    "exact-site evidence": TRAIN_SUBSETS["exact-site evidence"],
    "position validated": TRAIN_SUBSETS["position validated"],
    "surface-confirmed": TRAIN_SUBSETS["surface-confirmed"],
    "surface_and_function": TRAIN_SUBSETS["surface_and_function"],
    "confidence high": (
        lambda d: d["confidence"] == "high",
        "the literature agent's own confidence",
    ),
}
STRATA = [
    ("all", None, None),
    ("disordered", "disordered", None),
    ("ordered_surface", "ordered_surface", None),
    *[
        (f"{rc} / {tp}", rc, tp)
        for rc in ("disordered", "ordered_surface")
        for tp in ("single-pass", "multipass 2-6", "multipass 7+")
    ],
]


def gene_folds(train: pd.DataFrame) -> dict[str, int]:
    """Fold of every training gene, from the all-176 frame -- the production split."""
    y = train["label"].to_numpy()
    groups = train["gene_symbol"].to_numpy()
    out: dict[str, int] = {}
    for i, (_a, b) in enumerate(GroupKFold(n_splits=N_SPLITS).split(train, y, groups)):
        for g in set(groups[b]):
            out[g] = i
    return out


def fold_eval_fixed(train, ev, feats, kind, fold_of_gene):
    """`final_report.fold_eval` with the fold assignment taken as given.

    The split is by gene, decided once on the all-176 frame; a subset frame
    trains fold i on its own rows outside fold i and scores the same held-out
    genes. With the full frame this reproduces `fold_eval` exactly (asserted
    in main), so the base numbers are the production numbers.
    """
    y = train["label"].to_numpy()
    tr = within_protein(train, feats, "z")
    ed = within_protein(ev, feats, "z")
    trf = train["gene_symbol"].map(fold_of_gene).to_numpy()
    evf = ev["gene_symbol"].map(fold_of_gene).fillna(-1).astype(int).to_numpy()
    assert not np.isnan(trf.astype(float)).any(), "a training gene has no fold"
    sc = np.full(len(ev), np.nan)
    fold = np.full(len(ev), -1)
    for i in range(N_SPLITS):
        a = trf != i
        if len(np.unique(y[a])) < 2:
            continue
        m = make(kind, feats)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            m.fit(tr[a], y[a])
        mask = evf == i
        if mask.any():
            sc[mask] = m.predict_proba(ed[mask])[:, 1]
            fold[mask] = i
    return smooth(ev, sc, 2), fold


def subset_frame(train: pd.DataFrame, pos: pd.DataFrame, pred) -> pd.DataFrame:
    """The training frame for one tier: its positives, 200 negatives each.

    Negatives are a nested subsample of the all-176 negatives (same seed),
    restricted to proteins that keep at least one positive, so a protein whose
    every site is dropped leaves training entirely rather than staying as
    all-negative rows. Dropped positives are not in the negative pool either:
    `build_global` never drew them.
    """
    if pred is None:
        return train
    keep = pos[pred(pos)]
    keys = set(zip(keep["uniprot_acc"], keep["residue"], strict=False))
    is_pos = train["label"] == 1
    tp = train[
        is_pos
        & pd.Series(
            list(zip(train["uniprot_acc"], train["residue"], strict=False)), index=train.index
        ).isin(keys)
    ]
    assert len(tp) == len(keep), (len(tp), len(keep))
    neg = train[~is_pos & train["uniprot_acc"].isin(set(tp["uniprot_acc"]))]
    n = min(len(neg), RATIO * len(tp))
    rng = np.random.default_rng(SEED)
    neg = neg.iloc[np.sort(rng.choice(len(neg), size=n, replace=False))]
    return pd.concat([tp, neg], ignore_index=True)


def per_protein(d: pd.DataFrame) -> pd.DataFrame:
    """One row per protein with a positive: masked top-1, within-10 top-1, fold.

    `d` has columns acc, res, y, s, keep, fold. Masked top-1 ranks only kept
    residues and the positives; within-10 ranks everything and credits a
    first pick within WITHIN residues of any positive.
    """
    rows = []
    for acc, g in d.groupby("acc"):
        truth = g.loc[g["y"] == 1, "res"].to_numpy()
        if len(truth) == 0:
            continue
        masked = g[g["keep"] | (g["y"] == 1)].sort_values("s", ascending=False)
        full = g.sort_values("s", ascending=False)
        first = int(full["res"].iloc[0])
        rows.append(
            {
                "acc": acc,
                "fold": int(g["fold"].iloc[0]),
                "n_residues": len(g),
                "top1": bool(masked["y"].iloc[0] == 1),
                "top1_within10": bool(np.abs(truth - first).min() <= WITHIN),
                "rank_within10": int(
                    np.argmax(
                        np.abs(full["res"].to_numpy()[:, None] - truth[None, :]).min(1) <= WITHIN
                    )
                    + 1
                ),
            }
        )
    return pd.DataFrame(rows).set_index("acc")


def evaluate(ev, ev_y, sc, fold, keep):
    """Production metrics plus the per-protein table, on one score vector."""
    folds, topk, n_prot = metrics(ev, ev_y, sc, fold, keep)
    d = pd.DataFrame(
        {
            "acc": ev["uniprot_acc"],
            "res": ev["residue"],
            "y": ev_y,
            "s": sc,
            "keep": keep,
            "fold": fold,
        }
    )
    d = d[~np.isnan(sc)]
    pp = per_protein(d)
    assert len(pp) == n_prot, (len(pp), n_prot)
    return folds, topk, pp, d


def stratum_metrics(d: pd.DataFrame) -> tuple[dict, pd.DataFrame] | None:
    ok = d[d["keep"] | (d["y"] == 1)]
    if len(np.unique(ok["y"])) < 2:
        return None
    pp = per_protein(d)
    if len(pp) < MIN_PROTEINS:
        return None
    return {
        "proteins": len(pp),
        "sites": int(d["y"].sum()),
        "residues_per_protein": round(len(d) / len(pp)),
        "auroc": round(float(roc_auc_score(ok["y"], ok["s"])), 4),
        "top1": round(float(pp["top1"].mean()), 3),
        "top1_within10": round(float(pp["top1_within10"].mean()), 3),
    }, pp


def paired(folds, base_folds, pp, base_pp) -> dict:
    """Paired deltas against the base: per fold for AUROC and both recalls, per protein for within-10."""
    delta = folds - base_folds
    out = {
        "delta_auroc_mean": round(float(delta.mean()), 4),
        "delta_auroc_sd": round(float(delta.std(ddof=1)), 4),
        "folds_auroc_up": int((delta > 0).sum()),
        "folds_auroc_down": int((delta < 0).sum()),
    }
    both = pp.join(base_pp, lsuffix="", rsuffix="_base", how="inner")
    for col in ("top1", "top1_within10"):
        by_fold = both.groupby("fold")[[col, f"{col}_base"]].mean()
        df = by_fold[col] - by_fold[f"{col}_base"]
        out[f"delta_{col}"] = round(float(both[col].mean() - both[f"{col}_base"].mean()), 3)
        out[f"folds_{col}_up"] = int((df > 0).sum())
        out[f"folds_{col}_down"] = int((df < 0).sum())
        out[f"proteins_gained_{col}"] = int((both[col] & ~both[f"{col}_base"]).sum())
        out[f"proteins_lost_{col}"] = int((~both[col] & both[f"{col}_base"]).sum())
    return out


def own_labels(ev, ev_y, keep, pos_sub: pd.DataFrame):
    """Rows and labels for the own-subset design.

    Evaluation proteins are those with a kept positive; a lower-tier positive
    on such a protein is withheld (keep False, y 0) rather than counted as a
    negative. The +/-10 flank around it is already withheld by `mask_keep`,
    which was computed on all 176.
    """
    keys = set(zip(pos_sub["uniprot_acc"], pos_sub["residue"], strict=False))
    y_sub = np.array(
        [
            1 if (a, int(r)) in keys else 0
            for a, r in zip(ev["uniprot_acc"], ev["residue"], strict=False)
        ]
    )
    rows = ev["uniprot_acc"].isin(set(pos_sub["uniprot_acc"])).to_numpy()
    keep_sub = keep & ~((ev_y == 1) & (y_sub == 0))
    return rows, y_sub, keep_sub


def endonb_ranks(train_sub, allx, pos, feats) -> list[dict]:
    """Rank of each EndoNB site with its gene held out, on this training frame."""
    sites = json.loads((HERE / "endonb_sites.json").read_text())
    seen: set[str] = set()
    wanted = []
    for s in sites:
        if s["proteoform"] in set(allx["uniprot_acc"]) and s["gene"] not in seen:
            seen.add(s["gene"])
            wanted.append(s)
    out = []
    for s in wanted:
        acc, res = s["proteoform"], s["insert_after_residue"]
        sub = leave_gene_out(train_sub, allx, acc, feats)
        r_rank, _ = rank_of(sub, res, "rank_score")
        out.append(
            {
                "gene": s["gene"],
                "uniprot_acc": acc,
                "endonb_residue_ours": res,
                "extracellular_residues": len(sub),
                "rank": r_rank,
                "model_top_residue": int(sub["residue"].iloc[0]),
                "gene_in_this_training_frame": bool(
                    acc in set(train_sub.loc[train_sub["label"] == 1, "uniprot_acc"])
                ),
            }
        )
    return out


def main() -> int:
    pos, allx, train = training_frame()
    _, ev, ev_y = load()
    feats = [f for f in FEATS if f in train.columns]
    assert len(feats) == len(FEATS), "a production feature is missing from the training frame"
    pos = pos.reset_index(drop=True)
    n_pos_all = int((train["label"] == 1).sum())
    print(
        f"corpus: {n_pos_all} positives on {pos['uniprot_acc'].nunique()} proteins; "
        f"{len(train)} training rows; {len(ev)} evaluation residues"
    )
    for c in ("evidence_type", "position_evidence", "validation_level", "confidence"):
        print(f"  {c}: {pos[c].value_counts().to_dict()}")

    # ITGB5 is in neither the library nor the corpus; endonb_scored.py keeps its
    # rows beside the script.
    allx_endo = allx
    extra = HERE / "itgb5_extracellular.tsv"
    if extra.is_file():
        e = pd.read_csv(extra, sep="\t", low_memory=False)
        for c in ("ss3", "topology_class", "region_class"):
            e[c] = e[c].fillna("").astype(str)
        e["gene_symbol"] = e["gene_symbol"].fillna(e["uniprot_acc"]).astype(str)
        allx_endo = pd.concat([allx, e], ignore_index=True)

    fold_of_gene = gene_folds(train)
    keep = mask_keep(ev, ev_y)
    strata = pd.DataFrame(
        {
            "region_class": ev["region_class"].fillna("").astype(str),
            "topology": ev["tm_helix_count"].map(topo),
        }
    )

    # The base variant must reproduce the production evaluation exactly.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc_ref, fold_ref = fold_eval(train, ev, ev_y, feats, KIND)
    sc_base, fold_base = fold_eval_fixed(train, ev, feats, KIND, fold_of_gene)
    assert np.array_equal(fold_ref, fold_base), "fixed folds differ from final_report.fold_eval"
    assert np.allclose(sc_ref, sc_base, equal_nan=True), (
        "base scores differ from the production run"
    )
    print("base variant reproduces final_report.fold_eval: folds identical, scores identical")

    results, strata_rows, endo_rows = [], [], []
    hits: dict[str, pd.DataFrame] = {}
    scores: dict[str, np.ndarray] = {}
    frames: dict[str, pd.DataFrame] = {}
    for name, (pred, what) in SUBSETS.items():
        tr = subset_frame(train, pos, pred)
        frames[name] = tr
        pos_sub = pos if pred is None else pos[pred(pos)]
        if name == BASE:
            sc, fold = sc_base, fold_base
        else:
            sc, fold = fold_eval_fixed(tr, ev, feats, KIND, fold_of_gene)
        assert np.array_equal(fold, fold_base), "fold assignment differs between variants"
        scores[name] = sc
        n_pos = int((tr["label"] == 1).sum())
        n_prot_tr = tr.loc[tr["label"] == 1, "uniprot_acc"].nunique()
        print(
            f"\n{name:<22} {n_pos:>3} positives on {n_prot_tr:>3} proteins, "
            f"{len(tr) - n_pos} negatives ({what})"
        )

        for design in ("fixed", "own"):
            if design == "fixed":
                rows = np.ones(len(ev), bool)
                y_d, keep_d = ev_y, keep
            else:
                rows, y_d, keep_d = own_labels(ev, ev_y, keep, pos_sub)
            ev_d = ev[rows].reset_index(drop=True)
            folds, topk, pp, d = evaluate(ev_d, y_d[rows], sc[rows], fold[rows], keep_d[rows])
            # the base model on the same evaluation rows and labels
            bfolds, btopk, bpp, bd = evaluate(
                ev_d, y_d[rows], sc_base[rows], fold[rows], keep_d[rows]
            )
            key = f"{name}::{design}"
            hits[key] = pp
            row = {
                "variant": name,
                "design": design,
                "n_train_positives": n_pos,
                "n_train_proteins": n_prot_tr,
                "n_train_negatives": len(tr) - n_pos,
                "n_eval_sites": int(y_d[rows].sum()),
                "n_eval_proteins": len(pp),
                "auroc_mean": round(float(folds.mean()), 4),
                "auroc_sd": round(float(folds.std(ddof=1)), 4),
                "top1": round(topk[1], 3),
                "top3": round(topk[3], 3),
                "top10": round(topk[10], 3),
                "top1_within10": round(float(pp["top1_within10"].mean()), 3),
                "base_auroc_mean": round(float(bfolds.mean()), 4),
                "base_top1": round(btopk[1], 3),
                "base_top1_within10": round(float(bpp["top1_within10"].mean()), 3),
                **paired(folds, bfolds, pp, bpp),
                "folds": [round(float(f), 4) for f in folds],
                "base_folds": [round(float(f), 4) for f in bfolds],
            }
            results.append(row)
            print(
                f"  {design:<5} AUROC {row['auroc_mean']:.4f} ± {row['auroc_sd']:.4f} "
                f"({row['delta_auroc_mean']:+.4f}, {row['folds_auroc_up']}/{N_SPLITS} folds up)  "
                f"top-1 {row['top1']:.3f} ({row['delta_top1']:+.3f})  "
                f"within-10 {row['top1_within10']:.3f} ({row['delta_top1_within10']:+.3f}; "
                f"+{row['proteins_gained_top1_within10']}/-{row['proteins_lost_top1_within10']} "
                f"proteins)  on {row['n_eval_sites']} sites / {row['n_eval_proteins']} proteins"
            )
            # strata on the same scores
            s = strata[rows].reset_index(drop=True)
            okm = ~np.isnan(sc[rows])
            s = s[okm].reset_index(drop=True)
            d = d.reset_index(drop=True)
            bd = bd.reset_index(drop=True)
            for label, rc, tp in STRATA:
                m = np.ones(len(d), bool)
                if rc is not None:
                    m &= (s["region_class"] == rc).to_numpy()
                if tp is not None:
                    m &= (s["topology"] == tp).to_numpy()
                r = stratum_metrics(d[m])
                br = stratum_metrics(bd[m])
                if r is None or br is None:
                    continue
                sm, spp = r
                bm, bspp = br
                both = spp.join(bspp, lsuffix="", rsuffix="_base", how="inner")
                strata_rows.append(
                    {
                        "variant": name,
                        "design": design,
                        "stratum": label,
                        **sm,
                        "base_auroc": bm["auroc"],
                        "base_top1": bm["top1"],
                        "base_top1_within10": bm["top1_within10"],
                        "delta_auroc": round(sm["auroc"] - bm["auroc"], 4),
                        "delta_top1": round(sm["top1"] - bm["top1"], 3),
                        "delta_top1_within10": round(sm["top1_within10"] - bm["top1_within10"], 3),
                        "proteins_gained_within10": int(
                            (both["top1_within10"] & ~both["top1_within10_base"]).sum()
                        ),
                        "proteins_lost_within10": int(
                            (~both["top1_within10"] & both["top1_within10_base"]).sum()
                        ),
                    }
                )

        for e in endonb_ranks(tr, allx_endo, pos, feats):
            endo_rows.append({"variant": name, **e})

    out = pd.DataFrame(results)
    out.drop(columns=["folds", "base_folds"]).to_csv(
        HERE / "corpus_tier_retrain.tsv", sep="\t", index=False
    )
    st = pd.DataFrame(strata_rows)
    st.to_csv(HERE / "corpus_tier_retrain_strata.tsv", sep="\t", index=False)
    for design in ("fixed", "own"):
        print(f"\nby stratum, {design} evaluation (top-1 within 10):")
        print(
            st[st["design"] == design]
            .pivot(index="stratum", columns="variant", values="top1_within10")[list(SUBSETS)]
            .to_string()
        )
    endo = pd.DataFrame(endo_rows)
    endo.to_csv(HERE / "corpus_tier_retrain_endonb.tsv", sep="\t", index=False)
    print("\nEndoNB, gene held out (rank of the knock-in site on the rank score):")
    print(endo.pivot(index="gene", columns="variant", values="rank")[list(SUBSETS)].to_string())

    # per-protein hits on the fixed design, for the README
    gene = ev.drop_duplicates("uniprot_acc").set_index("uniprot_acc")["gene_symbol"]
    first_stratum = (
        pd.DataFrame(
            {
                "acc": ev["uniprot_acc"],
                "topology": strata["topology"],
                "region_class": strata["region_class"],
            }
        )
        .drop_duplicates("acc")
        .set_index("acc")
    )
    ph = pd.DataFrame(
        {
            "gene": gene,
            "topology": first_stratum["topology"],
            "region_class_first_residue": first_stratum["region_class"],
        }
    )
    for name in SUBSETS:
        pp = hits[f"{name}::fixed"]
        ph[f"within10::{name}"] = pp["top1_within10"]
        ph[f"rank_within10::{name}"] = pp["rank_within10"]
    ph = ph.dropna(subset=[f"within10::{BASE}"])
    ph.index.name = "uniprot_acc"
    ph.to_csv(HERE / "corpus_tier_retrain_protein_hits.tsv", sep="\t")

    (HERE / "corpus_tier_retrain.json").write_text(
        json.dumps(
            {
                "corpus_n_positives": n_pos_all,
                "corpus_n_proteins": int(pos["uniprot_acc"].nunique()),
                "protocol": "5-fold GroupKFold by gene fixed on the all-176 frame and asserted "
                "identical across variants; 13 features z-scored within protein; 200 "
                "negatives per positive, nested in the all-176 negatives, from proteins with "
                "a kept positive; +/-10 masking for AUROC and top1, within-10 credit for "
                "top1_within10, scores smoothed +/-2; rank model only",
                "subsets": {k: v[1] for k, v in SUBSETS.items()},
                "designs": {
                    "fixed": "scored on all 176 positives of the held-out genes",
                    "own": "scored on the subset's own positives; lower-tier positives on kept "
                    "proteins withheld from AUROC and masked top-1, no credit under within-10",
                },
                "results": results,
            },
            indent=2,
        )
        + "\n"
    )
    print(
        "\nwrote corpus_tier_retrain.{tsv,json}, corpus_tier_retrain_strata.tsv, "
        "corpus_tier_retrain_endonb.tsv, corpus_tier_retrain_protein_hits.tsv"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
