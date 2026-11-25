#!/usr/bin/env python3
"""Retrain the v4 ranker with labels that do not assume the published residue is the only good one.

A published site is one place a tag worked, not the only place it would. v4's
training frame (`two_tier_model.build_global`) excludes only the exact published
residue from the negative pool, so the site's neighbours and the rest of its
loop are sampled as negatives, while the evaluation masks negatives within
+/-10. This script fits five arms with the v4 estimator and features on one
fixed gene -> fold map and scores each the same way.

    v4               build_global as shipped (the reproduction check)
    mask10           negatives exclude residues within +/-10 of any published site
    segment_excluded negatives exclude every residue of a segment that hosts a site
    pu_bagged        bagged PU learning (Mordelet & Vert 2014): K forests, each on
                     the positives against a fresh random draw of unlabelled
                     residues, probabilities averaged
    two_stage        stage 1 ranks segments (region_table + region_oof from
                     two_tier_model, label = segment hosts a site); stage 2 orders
                     residues inside each segment by the v4 rank score
    hard_negative    published sites against near-certain failures only (buried
                     core, disulfide cysteines, N-glycan and sequon N, junctions
                     a domain-insertion screen calls intolerant); every other
                     residue is unlabelled and dropped
    hard_negative_two_stage
                     two_stage with stage 1 trained on host segments against
                     segments whose residues are ALL hard negatives, and stage 2
                     the hard_negative residue model

Why bagged PU rather than Elkan-Noto: Elkan-Noto divides the classifier's
probability by a constant c = P(labelled | positive), a monotone rescaling that
cannot change any within-protein ranking, and every metric here is a ranking.
It also assumes labelled positives are a random sample of all positives, which
published sites are not. Bagging changes the fitted model itself: an untested
good residue sits in only some of the unlabelled draws, so its pull on the
boundary is averaged down rather than fixed in.

The gene -> fold map is taken from GroupKFold on the v4 frame, the split
`final_report.fold_eval` makes, and reused for every arm; GroupKFold balances
groups by row count, so letting each arm split its own frame would change the
folds whenever the negatives change.

Metrics, out of fold, on the full cohort (176 sites, 127 proteins) and on the
111 exact-site positives, overall and by TM class:

    auroc_masked      per-fold residue AUROC, negatives within +/-10 masked
                      (final_report.metrics); by class, pooled within the class
    masked_top1       top residue after masking is a published site
    within10_top1/3   separated picks (>=25 apart) within 10 of a site
    segment_first     segment of the top residue hosts a site; multipass
                      proteins with >=2 segments of >=4 residues
    site_percentile   published residue's within-protein percentile against the
                      masked negatives (median over sites)

and, gene held out as `endonb_scored.py` does it, the five EndoNB knock-ins
and every screened protein with features (`dms_screen_labels.tsv`, from
`collect_screen_labels.py`): AUROC tolerant vs intolerant on its extracellular
junctions, and how many of the arm's three picks land on an intolerant one.

    uv run --locked python ranker_pu_segment.py               # all arms
    uv run --locked python ranker_pu_segment.py --arm mask10  # one arm
"""

from __future__ import annotations

import argparse
import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(REPO / "src"))

from endonb_scored import TOL, leave_gene_out, rank_of  # noqa: E402
from final_report import fold_eval, mask_keep  # noqa: E402
from independent_model_search import load, make, smooth, within_protein  # noqa: E402
from tag_site_model_v4 import FEATS, KIND, MIN_SEPARATION, N_SPLITS, RATIO, SEED  # noqa: E402
from two_tier_model import MASK_HALF, build_global, region_oof, region_table  # noqa: E402

EXACT = "published tag insertion at this exact site"
ARMS = (
    "v4",
    "mask10",
    "segment_excluded",
    "pu_bagged",
    "two_stage",
    "hard_negative",
    "hard_negative_two_stage",
)
TWO_STAGE_BASE = {"two_stage": "v4", "hard_negative_two_stage": "hard_negative"}
CLASSES = ("single-pass", "multipass 2-6", "multipass 7+")
HARD_RSA, HARD_PLDDT = 0.1, 80  # buried core: exposed to under a tenth, in a confident model
# domain-insertion screens (dms_screen_labels.tsv): ASIC1, KCNA3, P2RX3, Kir2.1
SCREENED_INSERTION = ("P78348", "P22001", "P56373", "P63252")
MIN_SEGMENT = 4  # residues; a shorter extracellular segment is not a separate choice
PU_BAGS = 50
PU_BAG_RATIO = 1  # unlabelled residues per positive in each bag (Mordelet & Vert use |P|)
OUT = HERE / "ranker_pu_segment.tsv"
OUT_PROT = HERE / "ranker_pu_segment_per_protein.tsv"
OUT_ENDONB = HERE / "ranker_pu_segment_endonb.tsv"
OUT_SCREENS = HERE / "ranker_pu_segment_screens.tsv"
OUT_JSON = HERE / "ranker_pu_segment.json"


def tm_class(n) -> str:
    n = int(n or 0)
    return "single-pass" if n <= 1 else ("multipass 2-6" if n <= 6 else "multipass 7+")


# --------------------------------------------------------------------------- frames


def frame():
    df, ev, ev_y = load()
    allx = pd.read_csv(HERE / "all_extracellular_targets.tsv", sep="\t", low_memory=False)
    for c in ("ss3", "topology_class", "region_class"):
        allx[c] = allx[c].fillna("").astype(str)
    allx["gene_symbol"] = allx["gene_symbol"].fillna(allx["uniprot_acc"]).astype(str)
    pos = df[df["label"] == 1].copy()
    return pos, allx, ev, ev_y


def sample_negatives(pos, allx, rng, ratio, drop):
    """build_global with an extra exclusion: `drop` marks pool rows that may not be negatives."""
    keys = set(zip(pos["uniprot_acc"], pos["residue"], strict=False))
    inpos = allx["uniprot_acc"].isin(set(pos["uniprot_acc"])).to_numpy()
    pool = allx[inpos & ~drop].copy()
    pool = pool[~pool.apply(lambda r: (r["uniprot_acc"], r["residue"]) in keys, axis=1)]
    n = min(len(pool), ratio * len(pos))
    neg = pool.iloc[rng.choice(len(pool), size=n, replace=False)].copy()
    neg["label"] = 0
    cols = [c for c in pos.columns if c in neg.columns]
    return pd.concat([pos[cols], neg[cols]], ignore_index=True)


def near_site(allx, pos, half=MASK_HALF):
    sites = pos.groupby("uniprot_acc")["residue"].apply(lambda v: np.array(sorted(v))).to_dict()
    out = np.zeros(len(allx), bool)
    for i, (a, r) in enumerate(zip(allx["uniprot_acc"], allx["residue"], strict=False)):
        s = sites.get(a)
        if s is not None and np.abs(s - int(r)).min() <= half:
            out[i] = True
    return out


def in_host_segment(allx, pos):
    host = set(zip(pos["uniprot_acc"], pos["region_index"], strict=False))
    return np.array(
        [(a, ri) in host for a, ri in zip(allx["uniprot_acc"], allx["region_index"], strict=False)]
    )


def training_set(arm, pos, allx):
    rng = np.random.default_rng(SEED)
    if arm in ("v4", "two_stage"):
        return build_global(pos, allx, rng, RATIO)
    if arm == "mask10":
        return sample_negatives(pos, allx, rng, RATIO, near_site(allx, pos))
    if arm == "segment_excluded":
        return sample_negatives(pos, allx, rng, RATIO, in_host_segment(allx, pos))
    if arm in ("hard_negative", "hard_negative_two_stage"):
        return hard_negative_set(pos, allx)
    raise ValueError(arm)


# --------------------------------------------------------------------------- constraints


def sequences() -> dict[str, str]:
    seqs = json.loads((HERE / "sequences_cache.json").read_text())
    uni = HERE / "uniprot_site_features.json"
    if uni.is_file():
        seqs = {**json.loads(uni.read_text())["sequences"], **seqs}
    return seqs


def load_screens() -> pd.DataFrame:
    return pd.read_csv(HERE / "dms_screen_labels.tsv", sep="\t")


def uniprot_positions(uni: pd.DataFrame):
    """acc -> {residue: partner} for disulfides, and acc -> sets for the other types."""
    u = uni.copy()
    u["begin"] = pd.to_numeric(u["begin"], errors="coerce")
    u["end"] = pd.to_numeric(u["end"], errors="coerce")
    partner: dict[str, dict[int, str]] = {}
    for a, b, e in u.loc[u["type"] == "DISULFID", ["uniprot_acc", "begin", "end"]].itertuples(
        index=False
    ):
        if b != b or e != e:
            continue
        b, e = int(b), int(e)
        d = partner.setdefault(a, {})
        if b == e:
            d[b] = "interchain"
        else:
            d[b], d[e] = str(e), str(b)
    glyc = u[(u["type"] == "CARBOHYD") & u["description"].fillna("").str.contains("N-linked")]
    glycan = glyc.groupby("uniprot_acc")["begin"].apply(lambda v: {int(x) for x in v if x == x})
    func = u[u["type"].isin(["BINDING", "ACT_SITE", "METAL", "SITE"])]
    sites: dict[str, dict[int, str]] = {}
    for a, t, b, e, desc in func[["uniprot_acc", "type", "begin", "end", "description"]].itertuples(
        index=False
    ):
        if b != b:
            continue
        for r in range(int(b), int(e if e == e else b) + 1):
            sites.setdefault(a, {})[r] = f"{t}{':' + desc if desc == desc and desc else ''}"
    return partner, glycan.to_dict(), sites


def annotate_constraints(rows: pd.DataFrame, screens: pd.DataFrame) -> pd.DataFrame:
    """One boolean column per near-certain reason a residue cannot take a tag.

    Residues inside the membrane are not among them: every row scored here is
    extracellular by construction (DeepTMHMM v1), and the per-residue features
    were only exported for extracellular residues, so a TM negative could only
    teach "is extracellular", which the candidate set already encodes.
    """
    uni = pd.read_csv(HERE / "uniprot_site_features.tsv", sep="\t")
    partner, glycan, func = uniprot_positions(uni)
    seqs = sequences()
    intol = screens[(screens["kind"] == "insertion") & (screens["tolerance_class"] == "intolerant")]
    intol_keys = set(zip(intol["uniprot_acc"], intol["residue"].astype(int), strict=False))
    out = rows.copy()
    acc = out["uniprot_acc"].to_numpy()
    res = out["residue"].astype(int).to_numpy()
    out["hn_buried"] = (
        (out["rsa"] < HARD_RSA) & (out["plddt"] >= HARD_PLDDT)
    ).to_numpy()  # NaN compares False
    out["hn_disulfide_cys"] = [r in partner.get(a, {}) for a, r in zip(acc, res, strict=False)]
    out["disulfide_partner"] = [
        partner.get(a, {}).get(r, "") for a, r in zip(acc, res, strict=False)
    ]
    out["hn_glycan_n"] = [r in glycan.get(a, set()) for a, r in zip(acc, res, strict=False)]

    def sequon(a, r):
        s = seqs.get(a, "")
        return r + 2 <= len(s) and s[r - 1] == "N" and s[r] != "P" and s[r + 1] in "ST"

    out["hn_sequon_n"] = [sequon(a, r) for a, r in zip(acc, res, strict=False)]
    out["hn_screen_intolerant"] = [(a, r) in intol_keys for a, r in zip(acc, res, strict=False)]
    out["functional_site"] = [func.get(a, {}).get(r, "") for a, r in zip(acc, res, strict=False)]
    out["hard_negative"] = out[[c for c in out.columns if c.startswith("hn_")]].any(axis=1)
    return out


def hard_negative_set(pos, allx):
    """Published sites against near-certain failures only; everything else is dropped.

    The pool is the published proteins plus the insertion-screened ones, whose
    intolerant junctions are the only measured failures available.
    """
    keys = set(zip(pos["uniprot_acc"], pos["residue"], strict=False))
    prots = set(pos["uniprot_acc"]) | set(SCREENED_INSERTION)
    pool = allx[allx["uniprot_acc"].isin(prots) & allx["hard_negative"]]
    pool = pool[
        [
            (a, int(r)) not in keys
            for a, r in zip(pool["uniprot_acc"], pool["residue"], strict=False)
        ]
    ]
    neg = pool.assign(label=0)
    cols = [c for c in pos.columns if c in neg.columns]
    return pd.concat([pos[cols], neg[cols]], ignore_index=True)


# --------------------------------------------------------------------------- fitting


def fixed_splits(groups, fold_of_gene):
    # genes outside the published cohort (screened proteins in the hard-negative
    # arm) carry no positive and are never evaluated here: training in every fold
    f = np.array([fold_of_gene.get(g, -1) for g in groups])
    return [(np.where(f != k)[0], np.where(f == k)[0]) for k in range(N_SPLITS)]


def oof_residue(train, ev, feats, fold_of_gene):
    """final_report.fold_eval with the split pinned to `fold_of_gene`."""
    y = train["label"].to_numpy()
    tr = within_protein(train, feats, "z")
    ed = within_protein(ev, feats, "z")
    gg = ev["gene_symbol"].to_numpy()
    sc = np.full(len(ev), np.nan)
    fold = np.full(len(ev), -1)
    for i, (a, b) in enumerate(fixed_splits(train["gene_symbol"].to_numpy(), fold_of_gene)):
        m = make(KIND, feats)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            m.fit(tr.iloc[a], y[a])
        mask = np.isin(gg, list(set(tr.iloc[b]["gene_symbol"])))
        sc[mask] = m.predict_proba(ed[mask])[:, 1]
        fold[mask] = i
    return smooth(ev, sc, 2), fold


def pu_frame(pos, allx, genes_out=()):
    """Positives plus EVERY unlabelled residue of their proteins, z-scored on that full set.

    z-scoring has to happen before bagging: a bag holds about one unlabelled
    residue per protein, so per-protein statistics computed inside it are noise.
    """
    keys = set(zip(pos["uniprot_acc"], pos["residue"], strict=False))
    p = pos[~pos["gene_symbol"].isin(genes_out)]
    u = allx[allx["uniprot_acc"].isin(set(p["uniprot_acc"]))]
    u = u[[(a, int(r)) not in keys for a, r in zip(u["uniprot_acc"], u["residue"], strict=False)]]
    u = u.assign(label=0)
    cols = [c for c in p.columns if c in u.columns]
    return pd.concat([p[cols], u[cols]], ignore_index=True)


def pu_fit_predict(tr_z, y, targets_z, feats, rng, bags=PU_BAGS, ratio=PU_BAG_RATIO):
    pi = np.where(y == 1)[0]
    ui = np.where(y == 0)[0]
    n = min(len(ui), ratio * len(pi))
    acc = [np.zeros(len(t)) for t in targets_z]
    for _ in range(bags):
        rows = np.concatenate([pi, rng.choice(ui, size=n, replace=False)])
        m = make(KIND, feats)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            m.fit(tr_z.iloc[rows], y[rows])
        for k, t in enumerate(targets_z):
            acc[k] += m.predict_proba(t)[:, 1]
    return [a / bags for a in acc]


def oof_pu(pos, allx, ev, feats, fold_of_gene):
    full = pu_frame(pos, allx)
    tr_z = within_protein(full, feats, "z")
    y = full["label"].to_numpy()
    ed = within_protein(ev, feats, "z")
    gg = ev["gene_symbol"].to_numpy()
    sc = np.full(len(ev), np.nan)
    fold = np.full(len(ev), -1)
    rng = np.random.default_rng(SEED)
    for i, (a, b) in enumerate(fixed_splits(full["gene_symbol"].to_numpy(), fold_of_gene)):
        mask = np.isin(gg, list(set(full.iloc[b]["gene_symbol"])))
        (p,) = pu_fit_predict(tr_z.iloc[a].reset_index(drop=True), y[a], [ed[mask]], feats, rng)
        sc[mask] = p
        fold[mask] = i
        print(f"    pu fold {i}: {len(a)} training rows, {mask.sum()} scored", flush=True)
    return smooth(ev, sc, 2), fold


def segment_tables(pos, rows, hard_labels=False):
    """region_table over every residue of each segment; label = segment hosts a site.

    With `hard_labels`, a segment is a negative only if EVERY residue in it is a
    hard negative; segments that are neither are unlabelled and left out of
    training (label NaN). Evaluation still scores every segment.
    """
    full = rows.copy()
    keys = set(zip(pos["uniprot_acc"], pos["residue"], strict=False))
    full["label"] = [
        1 if (a, int(r)) in keys else 0
        for a, r in zip(full["uniprot_acc"], full["residue"], strict=False)
    ]
    reg = region_table(full, FEATS)
    feats = [c for c in reg.columns if c.endswith(("_max", "_mean"))] + ["n_residues", "region_len"]
    if hard_labels:
        allhard = full.groupby(["uniprot_acc", "region_index"], dropna=False)["hard_negative"].all()
        h = [
            bool(allhard.get(k, False))
            for k in zip(reg["uniprot_acc"], reg["region_index"], strict=False)
        ]
        reg["label"] = np.where(reg["label"] == 1, 1.0, np.where(h, 0.0, np.nan))
    return reg, feats


def lexicographic(ev, seg_score, res_score):
    """Residues of the best segment first, by residue score; then the next segment.

    Returned as a within-protein percentile in (0, 1], so it can be pooled across
    proteins for AUROC; the order inside each protein is the whole content.
    """
    d = pd.DataFrame(
        {
            "acc": ev["uniprot_acc"].to_numpy(),
            "seg": seg_score,
            "res": res_score,
        }
    )
    d["seg"] = d["seg"].fillna(-1.0)
    out = np.full(len(d), np.nan)
    for _, g in d.groupby("acc"):
        ok = g["res"].notna()
        if not ok.any():
            continue
        gg = g[ok].sort_values(["seg", "res"], ascending=False)
        n = len(gg)
        out[gg.index.to_numpy()] = 1 - np.arange(n) / n
    return out


def oof_two_stage(pos, ev, v4_scores, v4_fold, fold_of_gene, train_rows=None):
    """Stage 1 out of fold on segments, then residues inside segments by `v4_scores`.

    `train_rows` (hard-negative variant) are the residues whose segments train
    stage 1, labelled with `hard_labels`; by default stage 1 trains on the
    evaluation proteins' own segments with label = hosts a site.
    """
    reg, rfeats = segment_tables(pos, ev)
    if train_rows is None:
        tr_reg = reg
    else:
        tr_reg, _ = segment_tables(pos, train_rows, hard_labels=True)
        tr_reg = tr_reg[tr_reg["label"].notna()].reset_index(drop=True)
    seg, _ = region_oof(tr_reg, reg, rfeats, fold_of_gene=fold_of_gene)
    seg_map = seg.to_dict()
    seg_vec = np.array(
        [
            seg_map.get((a, ri), np.nan)
            for a, ri in zip(ev["uniprot_acc"], ev["region_index"], strict=False)
        ]
    )
    reg = reg.assign(
        stage1=[seg_map.get(k) for k in zip(reg["uniprot_acc"], reg["region_index"], strict=False)]
    )
    ok = reg["stage1"].notna()
    seg_auc = roc_auc_score(reg.loc[ok, "label"], reg.loc[ok, "stage1"])
    return lexicographic(ev, seg_vec, v4_scores), v4_fold, float(seg_auc)


# --------------------------------------------------------------------------- metrics


def separated_top(res, sc, k=3, sep=MIN_SEPARATION):
    out: list[int] = []
    for i in np.argsort(-sc, kind="stable"):
        if len(out) >= k:
            break
        r = int(res[i])
        if any(abs(r - q) < sep for q in out):
            continue
        out.append(r)
    return out


def per_protein(ev, sc, truth, all_sites, prot_class):
    """One row per protein in `truth`: hits, segment choice, site percentiles."""
    seglen = ev.groupby(["uniprot_acc", "region_index"]).size()
    rows, site_rows = [], []
    d = ev[["uniprot_acc", "residue", "region_index"]].assign(s=sc)
    for acc, g in d.groupby("uniprot_acc"):
        if acc not in truth:
            continue
        g = g.dropna(subset=["s"])
        sites = np.array(sorted(truth[acc]))
        every = np.array(sorted(all_sites[acc]))
        res = g["residue"].to_numpy().astype(int)
        s = g["s"].to_numpy()
        dist_all = np.abs(res[:, None] - every[None, :]).min(1)
        neg = dist_all > MASK_HALF  # masked negatives: far from every published site
        is_site = np.isin(res, sites)
        keep = neg | is_site
        top_masked = res[keep][np.argmax(s[keep])]
        picks = separated_top(res, s)
        host = set(g.loc[g["residue"].isin(sites), "region_index"])
        top_seg = g["region_index"].to_numpy()[np.argmax(s)]
        nseg = int((seglen.loc[acc] >= MIN_SEGMENT).sum())
        cls = prot_class[acc]
        eligible = cls != "single-pass" and nseg >= 2
        pct = []
        for site in sites:
            j = np.where(res == site)[0]
            if not len(j):
                continue
            v = s[j[0]]
            n = int(neg.sum()) + 1
            r = int((s[neg] > v).sum()) + 1
            pct.append(100 * (1 - (r - 1) / n))
            site_rows.append({"uniprot_acc": acc, "residue": int(site), "percentile": pct[-1]})
        rows.append(
            {
                "uniprot_acc": acc,
                "gene_symbol": ev.loc[ev["uniprot_acc"] == acc, "gene_symbol"].iloc[0],
                "tm_class": cls,
                "n_sites": len(sites),
                "n_residues": len(g),
                "n_segments_ge4": nseg,
                "masked_top1": bool(top_masked in sites),
                "within10_top1": bool(np.abs(picks[0] - sites).min() <= MASK_HALF),
                "within10_top3": bool(min(np.abs(p - sites).min() for p in picks) <= MASK_HALF),
                "segment_eligible": eligible,
                "segment_first": bool(top_seg in host) if eligible else np.nan,
                "site_percentile_median": float(np.median(pct)) if pct else np.nan,
                "top_residue": int(res[np.argmax(s)]),
                "picks": ",".join(map(str, picks)),
                "sites": ",".join(map(str, sites)),
            }
        )
    return pd.DataFrame(rows), pd.DataFrame(site_rows)


def fold_aucs(ev_y, sc, fold, keep):
    ok = keep & ~np.isnan(sc)
    out = []
    for f in range(N_SPLITS):
        k = ok & (fold == f)
        if k.sum() > 50 and len(np.unique(ev_y[k])) == 2:
            out.append(roc_auc_score(ev_y[k], sc[k]))
    return np.array(out)


def summarise(arm, cohort, ev, ev_y_c, keep_c, sc, fold, prot, sites, prot_class):
    rows = []
    cls_vec = ev["uniprot_acc"].map(prot_class).to_numpy()
    for cls in ("all", *CLASSES):
        p = prot if cls == "all" else prot[prot["tm_class"] == cls]
        st = sites[sites["uniprot_acc"].isin(set(p["uniprot_acc"]))]
        cm = np.ones(len(ev), bool) if cls == "all" else (cls_vec == cls)
        ok = keep_c & cm & ~np.isnan(sc) & np.isin(ev["uniprot_acc"], list(p["uniprot_acc"]))
        if cls == "all":
            fa = fold_aucs(ev_y_c, np.where(ok, sc, np.nan), fold, ok)
            auc, auc_sd = float(fa.mean()), float(fa.std(ddof=1))
        else:
            auc = float(roc_auc_score(ev_y_c[ok], sc[ok])) if ev_y_c[ok].sum() else np.nan
            auc_sd = np.nan
        seg = p[p["segment_eligible"]]
        rows.append(
            {
                "arm": arm,
                "cohort": cohort,
                "tm_class": cls,
                "n_proteins": len(p),
                "n_sites": int(p["n_sites"].sum()),
                "auroc_masked": round(auc, 4),
                "auroc_masked_fold_sd": round(auc_sd, 4) if auc_sd == auc_sd else np.nan,
                "masked_top1": round(p["masked_top1"].mean(), 3),
                "within10_top1": round(p["within10_top1"].mean(), 3),
                "within10_top3": round(p["within10_top3"].mean(), 3),
                "n_segment_eligible": len(seg),
                "segment_first": round(seg["segment_first"].astype(float).mean(), 3)
                if len(seg)
                else np.nan,
                "site_percentile_median": round(float(st["percentile"].median()), 1),
                "site_percentile_top10pct": round(float((st["percentile"] >= 90).mean()), 3),
            }
        )
    return rows


# --------------------------------------------------------------------------- held-out proteins


def endonb_sites(allx):
    sites = json.loads((HERE / "endonb_sites.json").read_text())
    seen, wanted = set(), []
    for s in sites:
        if s["proteoform"] in set(allx["uniprot_acc"]) and s["gene"] not in seen:
            seen.add(s["gene"])
            wanted.append(s)
    return wanted


def with_itgb5(allx):
    extra = HERE / "itgb5_extracellular.tsv"
    if not extra.is_file():
        return allx
    e = pd.read_csv(extra, sep="\t", low_memory=False)
    for c in ("ss3", "topology_class", "region_class"):
        e[c] = e[c].fillna("").astype(str)
    e["gene_symbol"] = e["gene_symbol"].fillna(e["uniprot_acc"]).astype(str)
    return pd.concat([allx, e], ignore_index=True)


def region_fit_predict(tr_reg, ev_reg, cols):
    """region_oof's estimator, fitted once on every training region."""
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.impute import SimpleImputer
    from sklearn.pipeline import Pipeline

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
        pipe.fit(tr_reg[cols], tr_reg["label"].astype(int))
    return dict(zip(ev_reg["region_index"], pipe.predict_proba(ev_reg[cols])[:, 1], strict=False))


def held_out(arm, acc, ctx):
    """Every extracellular residue of `acc`, scored by `arm` refitted without its gene.

    `rank_score` is the arm's within-protein score (for the two-stage arms, the
    lexicographic percentile). Same leave-gene-out as `endonb_scored.py`.
    """
    key = (arm, acc)
    if key in ctx["cache"]:
        return ctx["cache"][key]
    allx_x, pos, feats = ctx["allx_x"], ctx["pos"], ctx["feats"]
    held = set(allx_x.loc[allx_x["uniprot_acc"] == acc, "gene_symbol"]) | {acc}
    p = pos[~pos["gene_symbol"].isin(held) & (pos["uniprot_acc"] != acc)]
    base = TWO_STAGE_BASE.get(arm, arm)
    if base == "pu_bagged":
        full = pu_frame(p, allx_x)
        sub = allx_x[allx_x["uniprot_acc"] == acc].copy().reset_index(drop=True)
        (pr,) = pu_fit_predict(
            within_protein(full, feats, "z"),
            full["label"].to_numpy(),
            [within_protein(sub, feats, "z")],
            feats,
            np.random.default_rng(SEED),
        )
        sub["rank_score"] = smooth(sub, pr, 2)
    else:
        sub = leave_gene_out(ctx["train"][base], allx_x, acc, feats)
    if arm in TWO_STAGE_BASE:
        sub = sub.sort_values("residue").reset_index(drop=True)
        if arm == "two_stage":
            rows = allx_x[allx_x["uniprot_acc"].isin(set(p["uniprot_acc"]))]
            tr_reg, rfeats = segment_tables(p, rows)
        else:
            prots = (set(p["uniprot_acc"]) | set(SCREENED_INSERTION)) - {acc}
            rows = allx_x[allx_x["uniprot_acc"].isin(prots) & ~allx_x["gene_symbol"].isin(held)]
            tr_reg, rfeats = segment_tables(p, rows, hard_labels=True)
            tr_reg = tr_reg[tr_reg["label"].notna()]
        ev_reg, _ = segment_tables(p, sub)
        seg = region_fit_predict(tr_reg, ev_reg, rfeats)
        seg_vec = np.array([seg.get(ri, np.nan) for ri in sub["region_index"]])
        sub["rank_score"] = lexicographic(sub, seg_vec, sub["rank_score"].to_numpy())
    ctx["cache"][key] = sub
    return sub


def endonb_arm(arm, ctx):
    rows = []
    for s in endonb_sites(ctx["allx_x"]):
        acc, res = s["proteoform"], s["insert_after_residue"]
        sub = held_out(arm, acc, ctx)
        n = len(sub)
        r, _ = rank_of(sub, res, "rank_score")
        top = sub.sort_values("rank_score", ascending=False).iloc[0]
        rows.append(
            {
                "arm": arm,
                "gene": s["gene"],
                "uniprot_acc": acc,
                "endonb_residue_ours": res,
                "extracellular_residues": n,
                "rank": r,
                "percentile": round(100 * (1 - (r - 1) / n), 1) if r else None,
                "model_top_residue": int(top["residue"]),
                "n_segments": int(sub["region_index"].nunique()),
            }
        )
    return rows


def screen_arm(arm, ctx, screens):
    """Measured tolerance on extracellular junctions of screened proteins, gene held out.

    AUROC is tolerant (1) against intolerant (0) on the arm's within-protein
    score; `top3_on_intolerant` counts the arm's three separated picks that land
    on a junction the screen calls intolerant, out of the picks the screen labels.
    """
    rows = []
    have = set(ctx["allx_x"]["uniprot_acc"])
    for (sid, kind, gene, acc), lab in screens.groupby(
        ["screen_id", "kind", "gene", "uniprot_acc"]
    ):
        if acc not in have:
            continue
        sub = held_out(arm, acc, ctx)
        m = sub[["residue", "rank_score"]].merge(
            lab[["residue", "tolerance_class"]].drop_duplicates("residue"), on="residue"
        )
        y = (m["tolerance_class"] == "tolerant").astype(int).to_numpy()
        auc = roc_auc_score(y, m["rank_score"]) if 0 < y.sum() < len(y) else np.nan
        picks = separated_top(sub["residue"].to_numpy(), sub["rank_score"].to_numpy())
        cls = dict(zip(m["residue"], m["tolerance_class"], strict=False))
        labelled = [p for p in picks if p in cls]
        rows.append(
            {
                "arm": arm,
                "screen_id": sid,
                "kind": kind,
                "gene": gene,
                "uniprot_acc": acc,
                "in_published_corpus": acc in set(ctx["pos"]["uniprot_acc"]),
                "n_extracellular_labelled": len(m),
                "n_tolerant": int(y.sum()),
                "n_intolerant": int(len(y) - y.sum()),
                "auroc_tolerant_vs_intolerant": round(float(auc), 4) if auc == auc else np.nan,
                "top3_picks": ",".join(map(str, picks)),
                "top3_labelled": len(labelled),
                "top3_on_intolerant": sum(cls[p] == "intolerant" for p in labelled),
            }
        )
    return rows


# --------------------------------------------------------------------------- main


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--arm", choices=ARMS, action="append")
    ap.add_argument("--skip-heldout", action="store_true", help="skip EndoNB and screens")
    args = ap.parse_args()
    arms = args.arm or list(ARMS)

    pos, allx_raw, ev, ev_y = frame()
    screens = load_screens()
    allx = annotate_constraints(allx_raw, screens)
    feats = [f for f in FEATS if f in pos.columns]
    assert len(feats) == 13, feats
    prot_class = pos.groupby("uniprot_acc")["tm_helix_count"].first().map(tm_class).to_dict()

    base = build_global(pos, allx, np.random.default_rng(SEED), RATIO)
    yb = base["label"].to_numpy()
    fold_of_gene: dict[str, int] = {}
    for i, (_a, b) in enumerate(
        GroupKFold(n_splits=N_SPLITS).split(base, yb, base["gene_symbol"].to_numpy())
    ):
        for g in set(base.iloc[b]["gene_symbol"]):
            fold_of_gene[g] = i
    print(f"fold map: {len(fold_of_gene)} genes over {N_SPLITS} folds", flush=True)

    # cohorts: truth sets and masks
    all_sites = pos.groupby("uniprot_acc")["residue"].apply(lambda v: {int(x) for x in v}).to_dict()
    ex = pos[pos["evidence_type"] == EXACT]
    exact_sites = (
        ex.groupby("uniprot_acc")["residue"].apply(lambda v: {int(x) for x in v}).to_dict()
    )
    keep_all = mask_keep(ev, ev_y)
    ev_y_ex = np.array(
        [
            1 if int(r) in exact_sites.get(a, ()) else 0
            for a, r in zip(ev["uniprot_acc"], ev["residue"], strict=False)
        ]
    )
    # exact cohort: positives = exact sites; every published site's +/-10 (and any
    # non-exact site itself) is neither class
    keep_ex = keep_all & ((ev_y == 0) | (ev_y_ex == 1))
    cohorts = {
        "all 176": (all_sites, ev_y, keep_all),
        "exact-site 111": (exact_sites, ev_y_ex, keep_ex),
    }
    ctx = {
        "pos": pos,
        "feats": feats,
        "allx_x": annotate_constraints(with_itgb5(allx_raw), screens),
        "train": {},
        "cache": {},
    }

    summary, prot_rows, endonb_rows, screen_rows, meta = [], [], [], [], {}
    scores: dict[str, tuple[np.ndarray, np.ndarray]] = {}

    def residue_arm(name):
        if name not in scores:
            tr = training_set(name, pos, allx)
            ctx["train"][name] = tr
            scores[name] = oof_residue(tr, ev, feats, fold_of_gene)
        return scores[name]

    for arm in arms:
        print(f"\n== {arm}", flush=True)
        if arm == "pu_bagged":
            sc, fold = oof_pu(pos, allx, ev, feats, fold_of_gene)
            meta[arm] = {"bags": PU_BAGS, "unlabelled_per_positive_per_bag": PU_BAG_RATIO}
        elif arm in TWO_STAGE_BASE:
            b = TWO_STAGE_BASE[arm]
            r_sc, r_fold = residue_arm(b)
            train_rows = None
            if arm == "hard_negative_two_stage":
                prots = set(pos["uniprot_acc"]) | set(SCREENED_INSERTION)
                train_rows = allx[allx["uniprot_acc"].isin(prots)]
                treg, _ = segment_tables(pos, train_rows, hard_labels=True)
                meta[arm] = {
                    "stage1_positive_segments": int((treg["label"] == 1).sum()),
                    "stage1_all_hard_negative_segments": int((treg["label"] == 0).sum()),
                    "stage1_unlabelled_segments": int(treg["label"].isna().sum()),
                }
            sc, fold, seg_auc = oof_two_stage(pos, ev, r_sc, r_fold, fold_of_gene, train_rows)
            meta.setdefault(arm, {})["stage1_segment_auroc"] = round(seg_auc, 4)
            print(f"    stage 1 segment AUROC {seg_auc:.4f}  {meta[arm]}", flush=True)
        else:
            sc, fold = residue_arm(arm)
            train = ctx["train"][arm]
            npos = int((train["label"] == 1).sum())
            meta[arm] = {"n_positive": npos, "n_negative": int(len(train) - npos)}
            if arm == "hard_negative":
                hn = allx[
                    allx["uniprot_acc"].isin(set(train["uniprot_acc"])) & allx["hard_negative"]
                ]
                meta[arm]["negatives_by_reason"] = {
                    c: int(hn[c].sum()) for c in hn.columns if c.startswith("hn_")
                }
                meta[arm]["unlabelled"] = "dropped"
            if arm == "v4":
                ref, _ = fold_eval(train, ev, ev_y, feats, KIND)
                same = bool(np.allclose(ref, sc, equal_nan=True))
                meta[arm]["matches_final_report_fold_eval"] = same
                assert same, "fixed folds do not reproduce final_report.fold_eval"
            print(f"    {meta[arm]}", flush=True)
        scores[arm] = (sc, fold)
        for cohort, (truth, yc, keep) in cohorts.items():
            prot, sites = per_protein(ev, sc, truth, all_sites, prot_class)
            summary += summarise(arm, cohort, ev, yc, keep, sc, fold, prot, sites, prot_class)
            prot_rows.append(prot.assign(arm=arm, cohort=cohort))
            r = summary[-4]
            print(
                f"    {cohort:<15} AUROC {r['auroc_masked']:.4f}  masked top1 {r['masked_top1']:.3f}  "
                f"w10 top1 {r['within10_top1']:.3f} top3 {r['within10_top3']:.3f}  "
                f"segment {r['segment_first']:.3f} (n={r['n_segment_eligible']})  "
                f"pct {r['site_percentile_median']}",
                flush=True,
            )
        if not args.skip_heldout:
            if arm not in ctx["train"] and TWO_STAGE_BASE.get(arm, arm) != "pu_bagged":
                residue_arm(TWO_STAGE_BASE.get(arm, arm))
            rows = endonb_arm(arm, ctx)
            endonb_rows += rows
            print(
                "    EndoNB "
                + "  ".join(f"{r['gene']} {r['rank']}/{r['extracellular_residues']}" for r in rows),
                flush=True,
            )
            srows = screen_arm(arm, ctx, screens)
            screen_rows += srows
            ins = [s for s in srows if s["kind"] == "insertion"]
            print(
                "    insertion screens "
                + "  ".join(
                    f"{s['screen_id'].split('_')[0]}/{s['screen_id'].split('_')[1][:6]} "
                    f"AUROC {s['auroc_tolerant_vs_intolerant']} "
                    f"intol {s['top3_on_intolerant']}/{s['top3_labelled']}"
                    for s in ins
                ),
                flush=True,
            )

    pd.DataFrame(summary).to_csv(OUT, sep="\t", index=False)
    pd.concat(prot_rows).to_csv(OUT_PROT, sep="\t", index=False)
    if endonb_rows:
        pd.DataFrame(endonb_rows).to_csv(OUT_ENDONB, sep="\t", index=False)
    if screen_rows:
        pd.DataFrame(screen_rows).to_csv(OUT_SCREENS, sep="\t", index=False)
    OUT_JSON.write_text(
        json.dumps(
            {
                "estimator": KIND,
                "features": feats,
                "seed": SEED,
                "negatives_per_positive": RATIO,
                "cv": f"GroupKFold by gene, {N_SPLITS} folds, gene->fold fixed from the v4 frame",
                "mask_half": MASK_HALF,
                "min_separation": MIN_SEPARATION,
                "min_segment": MIN_SEGMENT,
                "endonb_tolerance": TOL,
                "hard_negative_rule": {
                    "buried": f"rsa < {HARD_RSA} and plddt >= {HARD_PLDDT}",
                    "disulfide_cys": "UniProt DISULFID, both partners",
                    "glycan_n": "UniProt CARBOHYD N-linked",
                    "sequon_n": "N of N-X-S/T, X != P, from the sequence",
                    "screen_intolerant": "intolerant junction in a domain-insertion screen",
                    "membrane": "not used: every scored residue is extracellular by construction",
                },
                "arms": meta,
            },
            indent=2,
        )
        + "\n"
    )
    print(
        f"\nwrote {OUT.name}, {OUT_PROT.name}, {OUT_ENDONB.name}, {OUT_SCREENS.name}, {OUT_JSON.name}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
