#!/usr/bin/env python3
"""Hard negatives as a separate layer on top of v4, not as a replacement for its negatives.

A published site shows a tag worked there; it says nothing about where one would
fail. `ranker_pu_segment.py`'s hard-negative arm replaced v4's negatives with the
near-certain failures and dropped everything else. Here v4's training is kept and
the hard negatives are added on top, three ways:

    upweight_w{1,5,20}  v4's frame plus every hard negative on the published
                        proteins as an extra negative, weighted w (weight 1 for
                        the rest); all of v4's unlabelled negatives are kept
    veto_pm2            v4's out-of-fold scores; any residue within 2 of a hard
                        negative is removed before picking
    soft_x0.5, soft_x0.2  the same residues have their score multiplied instead
    pbad                a separate forest, hard negatives against every other
                        residue (P(bad)), on the absolute features; final score
                        v4 x (1 - P(bad))
    pu_bagged           ranker_pu_segment's bagged PU arm, for reference: 50 forests,
                        each on the positives against a fresh draw of one
                        unlabelled residue per positive
    pu_hard             bagged PU with the hard negatives kept apart: each bag is
                        the positives, one ordinary unlabelled residue per positive
                        and one hard negative per positive, so half of every bag's
                        negatives are known failures and half are the usual draw

Picks still sit at least 25 apart. Every arm also reports exact top-1: the
protein's top residue IS a published site, nothing masked. Hard negatives are defined as in
`ranker_pu_segment.annotate_constraints`: buried core (RSA < 0.1, pLDDT >= 80),
UniProt disulfide cysteines and partners, UniProt N-linked glycan N and sequon N,
and junctions a domain-insertion screen calls intolerant. Same gene -> fold map,
metrics and cohorts as `ranker_pu_segment.py`.

Two counts come first, because they bound what any of this can do:

    how many published sites (and their +/-10 windows) sit on a hard negative,
    by class -- if published sites sit on them often, a veto removes answers;
    how many of v4's top-1 / top-3 picks across the 566 L1 members land on one
    -- the size of the change a veto makes where the benchmark cannot score it.

    uv run --locked python ranker_hard_negative_layer.py
"""

from __future__ import annotations

import argparse
import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.model_selection import GroupKFold

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(REPO / "src"))

import ranker_pu_segment as R  # noqa: E402
from final_report import mask_keep  # noqa: E402
from independent_model_search import make  # noqa: E402
from l1_loop_evidence import member_rows  # noqa: E402
from l1_loop_evidence import score as full_fit_score  # noqa: E402
from two_tier_model import build_global  # noqa: E402

VETO_HALF = 2  # the scores are smoothed over +/-2, so a hard negative taints its neighbours
WEIGHTS = (1, 5, 20)
SOFT = (0.5, 0.2)
REASONS = {
    "buried core": ["hn_buried"],
    "disulfide Cys": ["hn_disulfide_cys"],
    "glycan N / sequon N": ["hn_glycan_n", "hn_sequon_n"],
    "insertion-intolerant": ["hn_screen_intolerant"],
    "any": ["hard_negative"],
}
INTERFACE_REASONS = {
    "ligand contact": ["hn_ligand_contact"],
    "protein-partner interface": ["hn_partner_interface"],
    "homomer interface": ["hn_homomer_interface"],
}
HOMOMER_CONTACT = 5.0  # Angstrom to the partner chain in the Schweke model: in contact
STRUCTURE_KEY = "sprint5/s5e13-tag-site-benchmarking/structure_global_features.tsv"
INTERFACES: str | bool = False  # "ligands" or "all", set by --ligands / --interfaces
OUT = HERE / "ranker_hard_negative_layer.tsv"
OUT_PAIRED = HERE / "ranker_hard_negative_layer_paired.tsv"
OUT_OVERLAP = HERE / "ranker_hard_negative_layer_overlap.tsv"
OUT_L1 = HERE / "ranker_hard_negative_layer_l1_picks.tsv"
OUT_JSON = HERE / "ranker_hard_negative_layer.json"


def structure_features() -> pd.DataFrame:
    """PR #223's per-residue structure table (d_oligomer_partner), cached under data/."""
    local = HERE / "data" / "structure_global_features.tsv"
    if not local.is_file():
        from deliverome_analysis import s3

        local.parent.mkdir(parents=True, exist_ok=True)
        s3.processed.download(STRUCTURE_KEY, str(local))
    return pd.read_csv(local, sep="\t", usecols=["uniprot_acc", "residue", "d_oligomer_partner"])


def constraints(rows: pd.DataFrame, screens: pd.DataFrame) -> pd.DataFrame:
    """annotate_constraints, plus interface classes when --ligands or --interfaces is set.

    Ligand contact and protein-partner interface come from D1's
    `contact_observation`, classified by `fetch_interface_contacts.py` and moved
    onto each L1 member's ordered protein by `map_contacts_to_members.py`
    (sequence alignment, identical residues only); antibody epitopes are never
    used. Homomer interface is a residue within HOMOMER_CONTACT of the partner
    chain in the Schweke 2024 model (PR #223), computed for the published corpus.
    """
    out = R.annotate_constraints(rows, screens)
    if not INTERFACES:
        return out
    from map_contacts_to_members import contact_keys

    key = list(zip(out["uniprot_acc"], out["residue"].astype(int), strict=False))
    lig = contact_keys({"ligand"})
    out["hn_ligand_contact"] = [k in lig for k in key]
    if INTERFACES == "all":
        partner = contact_keys({"protein_partner"})
        out["hn_partner_interface"] = [k in partner for k in key]
        sf = structure_features()
        d = out[["uniprot_acc", "residue"]].merge(sf, on=["uniprot_acc", "residue"], how="left")
        out["hn_homomer_interface"] = (d["d_oligomer_partner"] <= HOMOMER_CONTACT).to_numpy()
    out["hard_negative"] = out[[c for c in out.columns if c.startswith("hn_")]].any(axis=1)
    return out


def flags(rows: pd.DataFrame, cols: list[str]) -> pd.Series:
    return rows[cols].any(axis=1)


def near(rows: pd.DataFrame, flag: pd.Series, half: int, hi: int | None = None) -> np.ndarray:
    """True where a flagged residue of the same protein lies in [n - half, n + hi].

    When ``hi`` is None, uses a symmetric +/-``half`` window around residue ``n``.
    For an insertion after residue ``n`` (junction between ``n`` and ``n+1``), pass
    ``hi=half + 1`` so the window covers ``n - half .. n + half + 1`` (e.g. ``n-2..n+3``).
    """
    upper = half if hi is None else hi
    out = np.zeros(len(rows), bool)
    res = rows["residue"].astype(int).to_numpy()
    f = flag.to_numpy(dtype=bool)
    for _, idx in rows.groupby("uniprot_acc").indices.items():
        bad = res[idx][f[idx]]
        if len(bad):
            d = bad[None, :] - res[idx][:, None]
            out[idx] = ((d >= -half) & (d <= upper)).any(1)
    return out


def smooth_contiguous(frame: pd.DataFrame, scores: np.ndarray, half: int = 2) -> np.ndarray:
    """Smooth scores within each contiguous extracellular residue run (never across TM gaps)."""
    out = np.asarray(scores, dtype=float).copy()
    if half <= 0:
        return out
    ker = np.ones(2 * half + 1) / (2 * half + 1)
    tmp = frame[["uniprot_acc", "residue"]].copy()
    tmp["s"] = out
    for _, g in tmp.groupby("uniprot_acc"):
        g = g.sort_values("residue", kind="mergesort")
        res = g["residue"].astype(int).to_numpy()
        vals = g["s"].to_numpy(dtype=float)
        idx = g.index.to_numpy()
        splits = np.where(np.diff(res) != 1)[0] + 1
        start = 0
        for end in [*splits.tolist(), len(res)]:
            seg = vals[start:end]
            pad = np.pad(seg, half, mode="edge")
            out[idx[start:end]] = np.convolve(pad, ker, mode="valid")
            start = end
    return out


def zscore_from_reference(df: pd.DataFrame, ref: pd.DataFrame, feats: list[str]) -> pd.DataFrame:
    """Standardize ``feats`` per protein using mean/std from the full extracellular ``ref``.

    Prevents train/eval z-score shift when ``df`` is a subsampled training frame and
    ``ref`` holds all extracellular residues. Proteins absent from ``ref`` fall back to
    their own per-protein mean/std in ``df``.
    """
    out = df.copy()
    cols = [f for f in feats if f in out.columns and f in ref.columns]
    if not cols:
        return out
    ref_g = ref.groupby("uniprot_acc")[cols]
    ref_mu = ref_g.mean()
    ref_sd = ref_g.std().replace(0, np.nan)
    own_g = out.groupby("uniprot_acc")
    for f in cols:
        mu = out["uniprot_acc"].map(ref_mu[f])
        sd = out["uniprot_acc"].map(ref_sd[f])
        missing = ~out["uniprot_acc"].isin(ref_mu.index)
        if missing.any():
            mu = mu.where(~missing, own_g[f].transform("mean"))
            sd = sd.where(~missing, own_g[f].transform("std").replace(0, np.nan))
        out[f] = (out[f] - mu) / sd
    return out


# --------------------------------------------------------------------------- counts


def overlap_table(ev, pos) -> pd.DataFrame:
    """Published sites on, near, or with a hard negative inside their +/-10 window."""
    rows = []
    key = ev.set_index(["uniprot_acc", "residue"])
    for cohort, p in (
        ("exact-site 111", pos[pos["evidence_type"] == R.EXACT]),
        ("all 176", pos),
    ):
        for reason, cols in REASONS.items():
            f = flags(ev, cols)
            on = key.index.isin(list(zip(p["uniprot_acc"], p["residue"], strict=False)))
            sel = ev[on].copy()
            sel_flag = f[on].to_numpy()
            pm2 = near(ev, f, VETO_HALF)[on]
            pm10 = near(ev, f, R.MASK_HALF)[on]
            rows.append(
                {
                    "cohort": cohort,
                    "reason": reason,
                    "n_sites": len(sel),
                    "site_is_hard_negative": int(sel_flag.sum()),
                    "hard_negative_within_2": int(pm2.sum()),
                    "hard_negative_within_10": int(pm10.sum()),
                    "frac_within_2": round(float(pm2.mean()), 3),
                    "frac_within_10": round(float(pm10.mean()), 3),
                }
            )
    return pd.DataFrame(rows)


def l1_pick_table(pos, allx_raw, screens) -> tuple[pd.DataFrame, dict]:
    """Where v4's three picks land on the 566 members, refitted on all 176 sites."""
    allx = constraints(allx_raw, screens)
    targets = pd.read_csv(HERE / "l1_member_targets.tsv", sep="\t")
    rows = constraints(member_rows(allx_raw), screens)
    rows = rows[rows["uniprot_acc"].isin(set(targets["uniprot_acc"]))].reset_index(drop=True)
    rows["v4"] = full_fit_score("v4", pos, allx, rows)
    pm2 = near(rows, rows["hard_negative"], VETO_HALF)
    rows["near_hn"] = pm2
    out = []
    uni = set(pd.read_csv(HERE / "uniprot_site_features.tsv", sep="\t")["uniprot_acc"])
    for t in targets.itertuples(index=False):
        g = rows[rows["uniprot_acc"] == t.uniprot_acc]
        if g.empty:
            continue
        picks = R.separated_top(g["residue"].to_numpy(), g["v4"].to_numpy())
        gv = g.set_index("residue")
        veto = g[~g["near_hn"]]
        vpicks = R.separated_top(veto["residue"].to_numpy(), veto["v4"].to_numpy())
        rec = {
            "candidate_id": t.candidate_id,
            "gene": t.gene,
            "uniprot_acc": t.uniprot_acc,
            "has_uniprot_annotation": t.uniprot_acc in uni,
            "picks": ",".join(map(str, picks)),
            "picks_after_veto": ",".join(map(str, vpicks)),
            "top1_changes_under_veto": bool(vpicks[:1] != picks[:1]),
        }
        for reason, cols in REASONS.items():
            tag = reason.split(" ")[0].replace("-", "_").replace("/", "_").lower()
            f = flags(gv, cols)
            fn = pd.Series(near(g, flags(g, cols), VETO_HALF), index=g["residue"].to_numpy())
            rec[f"top1_on_{tag}"] = bool(f.loc[picks[0]])
            rec[f"top1_within2_{tag}"] = bool(fn.loc[picks[0]])
            rec[f"top3_on_{tag}"] = int(sum(bool(f.loc[p]) for p in picks))
            rec[f"top3_within2_{tag}"] = int(sum(bool(fn.loc[p]) for p in picks))
        out.append(rec)
    df = pd.DataFrame(out)
    summary = {
        "members": len(df),
        "members_without_uniprot_annotation": int((~df["has_uniprot_annotation"]).sum()),
        "picks": int(df["picks"].str.count(",").add(1).sum()),
        "top1_changes_under_veto": int(df["top1_changes_under_veto"].sum()),
    }
    for reason in REASONS:
        tag = reason.split(" ")[0].replace("-", "_").replace("/", "_").lower()
        summary[reason] = {
            "top1_on": int(df[f"top1_on_{tag}"].sum()),
            "top1_within_2": int(df[f"top1_within2_{tag}"].sum()),
            "top3_picks_on": int(df[f"top3_on_{tag}"].sum()),
            "top3_picks_within_2": int(df[f"top3_within2_{tag}"].sum()),
        }
    return df, summary


# --------------------------------------------------------------------------- arms


def oof_weighted(train, w, ev, feats, fold_of_gene):
    """oof_residue with per-row sample weights passed to the forest."""
    y = train["label"].to_numpy()
    tr = zscore_from_reference(train, ev, feats)
    ed = zscore_from_reference(ev, ev, feats)
    gg = ev["gene_symbol"].to_numpy()
    sc = np.full(len(ev), np.nan)
    fold = np.full(len(ev), -1)
    for i, (a, _b) in enumerate(R.fixed_splits(train["gene_symbol"].to_numpy(), fold_of_gene)):
        m = make(R.KIND, feats)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            m.fit(tr.iloc[a], y[a], clf__sample_weight=w[a])
        # score every evaluation gene of this fold, including a gene the training frame
        # lacks (the corrected corpus drops 21 proteins' positives entirely)
        mask = np.array([fold_of_gene.get(g, -1) == i for g in gg])
        sc[mask] = m.predict_proba(ed[mask])[:, 1]
        fold[mask] = i
    return smooth_contiguous(ev, sc, 2), fold


def upweight_frame(base, pos, allx):
    """v4's frame, plus every hard negative on a published protein that it did not draw."""
    hn = allx[allx["uniprot_acc"].isin(set(pos["uniprot_acc"])) & allx["hard_negative"]]
    pkeys = set(zip(pos["uniprot_acc"], pos["residue"], strict=False))
    hn = hn[
        [(a, int(r)) not in pkeys for a, r in zip(hn["uniprot_acc"], hn["residue"], strict=False)]
    ]
    hkeys = set(zip(hn["uniprot_acc"], hn["residue"].astype(int), strict=False))
    in_base = set(zip(base["uniprot_acc"], base["residue"].astype(int), strict=False))
    extra = hn[
        [(a, int(r)) not in in_base for a, r in zip(hn["uniprot_acc"], hn["residue"], strict=False)]
    ]
    extra = extra.assign(label=0)[[c for c in base.columns if c in extra.columns]]
    tr = pd.concat([base, extra], ignore_index=True)
    is_hn = np.array(
        [
            (a, int(r)) in hkeys and lab == 0
            for a, r, lab in zip(tr["uniprot_acc"], tr["residue"], tr["label"], strict=False)
        ]
    )
    return tr, is_hn, len(extra)


def oof_pbad(ev, feats, fold_of_gene):
    """P(bad): hard negatives against every other residue, absolute features, out of fold."""
    y = ev["hard_negative"].astype(int).to_numpy()
    gg = ev["gene_symbol"].to_numpy()
    out = np.full(len(ev), np.nan)
    for a, b in R.fixed_splits(gg, fold_of_gene):
        m = make(R.KIND, feats)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            m.fit(ev.iloc[a], y[a])
        out[b] = m.predict_proba(ev.iloc[b])[:, 1]
    return out


def oof_pu_hard(pos, allx, ev, feats, fold_of_gene, bags=R.PU_BAGS):
    """Bagged PU in which every bag also holds hard negatives, |P| of each kind."""
    full = R.pu_frame(pos, allx)
    hk = set(
        zip(
            allx.loc[allx["hard_negative"], "uniprot_acc"],
            allx.loc[allx["hard_negative"], "residue"].astype(int),
            strict=False,
        )
    )
    is_hn = np.array(
        [
            lab == 0 and (a, int(r)) in hk
            for a, r, lab in zip(full["uniprot_acc"], full["residue"], full["label"], strict=False)
        ]
    )
    tr_z = zscore_from_reference(full, ev, feats)
    y = full["label"].to_numpy()
    ed = zscore_from_reference(ev, ev, feats)
    gg = ev["gene_symbol"].to_numpy()
    sc = np.full(len(ev), np.nan)
    fold = np.full(len(ev), -1)
    rng = np.random.default_rng(R.SEED)
    for i, (a, b) in enumerate(R.fixed_splits(full["gene_symbol"].to_numpy(), fold_of_gene)):
        pi = a[y[a] == 1]
        hi = a[(y[a] == 0) & is_hn[a]]
        ui = a[(y[a] == 0) & ~is_hn[a]]
        mask = np.isin(gg, list(set(full.iloc[b]["gene_symbol"])))
        acc = np.zeros(int(mask.sum()))
        for _ in range(bags):
            rows = np.concatenate(
                [
                    pi,
                    rng.choice(ui, size=min(len(ui), len(pi)), replace=False),
                    rng.choice(hi, size=min(len(hi), len(pi)), replace=False),
                ]
            )
            m = make(R.KIND, feats)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                m.fit(tr_z.iloc[rows], y[rows])
            acc += m.predict_proba(ed[mask])[:, 1]
        sc[mask] = acc / bags
        fold[mask] = i
    return smooth_contiguous(ev, sc, 2), fold, int(is_hn.sum())


def exact_top1(prot: pd.DataFrame) -> pd.Series:
    """The protein's single top residue is a published site, nothing masked."""
    return pd.Series(
        [str(t) in s.split(",") for t, s in zip(prot["top_residue"], prot["sites"], strict=False)],
        index=prot.index,
    )


# --------------------------------------------------------------------------- main


def main() -> int:
    global INTERFACES, OUT, OUT_PAIRED, OUT_OVERLAP, OUT_L1, OUT_JSON
    ap = argparse.ArgumentParser()
    group = ap.add_mutually_exclusive_group()
    group.add_argument(
        "--ligands", action="store_true", help="add D1 ligand contacts to the hard negatives"
    )
    group.add_argument(
        "--interfaces",
        action="store_true",
        help="add ligand contacts, protein-partner and homomer interfaces to the hard negatives",
    )
    args = ap.parse_args()
    if args.ligands or args.interfaces:
        INTERFACES = "all" if args.interfaces else "ligands"
        tag = "_interfaces" if args.interfaces else "_ligands"
        REASONS.update(
            INTERFACE_REASONS if args.interfaces else {"ligand contact": ["hn_ligand_contact"]}
        )
        REASONS["any"] = REASONS.pop("any")
        OUT, OUT_PAIRED, OUT_OVERLAP, OUT_L1, OUT_JSON = (
            p.with_name(p.stem + tag + p.suffix)
            for p in (OUT, OUT_PAIRED, OUT_OVERLAP, OUT_L1, OUT_JSON)
        )
    pos, allx_raw, ev, ev_y = R.frame()
    screens = R.load_screens()
    allx = constraints(allx_raw, screens)
    ev = constraints(ev, screens)
    feats = [f for f in R.FEATS if f in pos.columns]
    prot_class = pos.groupby("uniprot_acc")["tm_helix_count"].first().map(R.tm_class).to_dict()

    base = build_global(pos, allx_raw, np.random.default_rng(R.SEED), R.RATIO)
    fold_of_gene: dict[str, int] = {}
    for i, (_a, b) in enumerate(
        GroupKFold(n_splits=R.N_SPLITS).split(base, base["label"], base["gene_symbol"])
    ):
        for g in set(base.iloc[b]["gene_symbol"]):
            fold_of_gene[g] = i

    overlap = overlap_table(ev, pos)
    overlap.to_csv(OUT_OVERLAP, sep="\t", index=False)
    print(overlap.to_string(index=False), flush=True)

    l1, l1_summary = l1_pick_table(pos, allx_raw, screens)
    l1.to_csv(OUT_L1, sep="\t", index=False)
    print(json.dumps(l1_summary, indent=1), flush=True)

    all_sites = pos.groupby("uniprot_acc")["residue"].apply(lambda v: {int(x) for x in v}).to_dict()
    ex = pos[pos["evidence_type"] == R.EXACT]
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
    keep_ex = keep_all & ((ev_y == 0) | (ev_y_ex == 1))
    cohorts = {
        "all 176": (all_sites, ev_y, keep_all),
        "exact-site 111": (exact_sites, ev_y_ex, keep_ex),
    }

    v4, fold = R.oof_residue(base, ev, feats, fold_of_gene)
    pm2 = near(ev, ev["hard_negative"], VETO_HALF)
    arms: dict[str, np.ndarray] = {"v4": v4}
    meta: dict = {"veto_half": VETO_HALF, "residues_within_2_of_hard_negative": int(pm2.sum())}

    tr, is_hn, n_extra = upweight_frame(base, pos, allx)
    meta["upweight"] = {"hard_negatives_in_frame": int(is_hn.sum()), "added_to_v4_frame": n_extra}
    for w in WEIGHTS:
        weights = np.where(is_hn, float(w), 1.0)
        arms[f"upweight_w{w}"], _ = oof_weighted(tr, weights, ev, feats, fold_of_gene)
        print(f"upweight w={w} done", flush=True)

    arms["veto_pm2"] = np.where(pm2, -1.0, v4)
    for m in SOFT:
        arms[f"soft_x{m}"] = np.where(pm2, v4 * m, v4)

    pbad = oof_pbad(ev, feats, fold_of_gene)
    from sklearn.metrics import roc_auc_score

    meta["pbad_auroc_oof"] = round(float(roc_auc_score(ev["hard_negative"], pbad)), 4)
    arms["pbad"] = v4 * (1 - pbad)
    arms["pu_bagged"], _ = R.oof_pu(pos, allx_raw, ev, feats, fold_of_gene)
    arms["pu_hard"], _, n_hn = oof_pu_hard(pos, allx, ev, feats, fold_of_gene)
    meta["pu_hard"] = {"hard_negatives_in_unlabelled_pool": n_hn, "bags": R.PU_BAGS}

    summary, prot_rows = [], []
    for arm, sc in arms.items():
        for cohort, (truth, yc, keep) in cohorts.items():
            prot, sites = R.per_protein(ev, sc, truth, all_sites, prot_class)
            rows = R.summarise(arm, cohort, ev, yc, keep, sc, fold, prot, sites, prot_class)
            ex1 = exact_top1(prot)
            for r in rows:
                sel = prot["tm_class"] == r["tm_class"] if r["tm_class"] != "all" else ex1 == ex1
                r["exact_top1"] = round(float(ex1[sel].mean()), 3)
            summary += rows
            prot = prot.assign(exact_top1=ex1)
            prot_rows.append(prot.assign(arm=arm, cohort=cohort))
        r = [s for s in summary if s["arm"] == arm and s["tm_class"] == "all"]
        print(
            f"{arm:<14} "
            + "  ".join(
                f"{x['cohort']}: exact {x['exact_top1']:.3f} w10 {x['within10_top1']:.3f} "
                f"top3 {x['within10_top3']:.3f} "
                f"masked {x['masked_top1']:.3f} seg {x['segment_first']:.3f}"
                for x in r
            ),
            flush=True,
        )
    pd.DataFrame(summary).to_csv(OUT, sep="\t", index=False)

    prot = pd.concat(prot_rows)
    paired = []
    for cohort in cohorts:
        q = prot[prot["cohort"] == cohort]
        ref = q[q["arm"] == "v4"].set_index("uniprot_acc")
        for arm in arms:
            if arm == "v4":
                continue
            a = q[q["arm"] == arm].set_index("uniprot_acc")
            row = {"cohort": cohort, "arm": arm}
            for m in (
                "exact_top1",
                "masked_top1",
                "within10_top1",
                "within10_top3",
                "segment_first",
            ):
                x, y = a[m].astype(float), ref[m].astype(float)
                ok = x.notna() & y.notna()
                row[f"{m}_gained"] = int(((x > y) & ok).sum())
                row[f"{m}_lost"] = int(((x < y) & ok).sum())
            paired.append(row)
    pd.DataFrame(paired).to_csv(OUT_PAIRED, sep="\t", index=False)
    OUT_JSON.write_text(
        json.dumps(
            {
                "estimator": R.KIND,
                "seed": R.SEED,
                "cv": "GroupKFold by gene, 5 folds, gene->fold fixed from the v4 frame",
                "weights": list(WEIGHTS),
                "soft_multipliers": list(SOFT),
                "l1_picks": l1_summary,
                **meta,
            },
            indent=2,
        )
        + "\n"
    )
    print(json.dumps(meta, indent=1))
    print(
        f"wrote {OUT.name}, {OUT_PAIRED.name}, {OUT_OVERLAP.name}, {OUT_L1.name}, {OUT_JSON.name}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
