#!/usr/bin/env python3
"""Which hard negatives the model of record trains on: a benchmark that cannot pick a winner alone.

Becca, 2026-11-15: "even if v5 doesn't do better with lit it is still good to consider adding
hard negatives given the sparse data we are training on" and "we want to consider ligand and co
receptors / homomer contacts as well". The coordinator, under her delegation: the retrained model
is adopted on hypothesis grounds unless it clearly hurts, meaning the bootstrap interval of its
top-3 or region-hit difference from b0 lies wholly below zero.

Training data: v5's frame (176 published internal sites on 127 proteins, 200 unlabelled
extracellular residues per site, three draws), with the 5 positives on 4 proteins (NRXN1, ASIC1,
SLC4A11, SLC9A6) that were numbered on another isoform re-featured from the proteoform-fixed
export (c7be01a4 exporter, data/allx_isoform_fix15_2026-11-14.tsv).

Hard-negative sets, each added to the training frame at its own weight:

    b0          v5's classes (buried core, disulfide Cys, glycan/sequon N, screen-intolerant,
                D1 ligand contacts) at w_d1
    b1          b0 + literature ligand contacts with structure evidence, at w_x
    b2          b0 + every verified literature ligand contact, at w_x
    c_partner   b0 + D1 protein-partner (co-receptor) interface residues, at w_x
    d_homomer   c_partner + homomer interface residues (Schweke assembly models, the corpus
                proteins only), at w_x
    e           b1 + partner + homomer: every class Becca asked to consider, at w_x

w_d1 is 5 (v5) or 10; w_x is half of w_d1 (lower than D1, as proposed). Every set is also scored
veto-only: b0's scores with the extra class's residues removed before picking (window n-2..n+3,
as in the decisions). Floors: a random extracellular residue, and the midpoint of the longest
extracellular region.

Metrics, out of fold, gene-grouped and HGNC-family-grouped 5-fold: within-10 top-1 and top-3,
region hit top-1 and top-3 (a published site lies in a picked region). Mean and SD over the
three negative draws; 95% bootstrap interval over proteins (2,000 resamples) for each variant's
difference from b0 at the same weight. The literature reaches a published protein only where its
L1 member's ordered protein is that protein's sequence (identical), so its effect is measured on
few proteins; the number is reported.

    uv run --locked python hard_negative_benchmark.py            # writes hard_negative_benchmark.*
    uv run --locked python hard_negative_benchmark.py --quick    # one draw, gene folds only
"""

from __future__ import annotations

import argparse
import json
import sys
import time
import urllib.request
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[3] / "src"))

import ranker_hard_negative_layer as L  # noqa: E402
import ranker_pu_segment as R  # noqa: E402
from corpus_tier_retrain import gene_folds  # noqa: E402
from sklearn.model_selection import GroupKFold  # noqa: E402
from two_tier_model import build_global  # noqa: E402

DRAWS = (R.SEED, 1, 2)
OUT = HERE / "hard_negative_benchmark.tsv"
OUT_PROT = HERE / "hard_negative_benchmark_per_protein.tsv"
OUT_JSON = HERE / "hard_negative_benchmark.json"
FAMILIES = HERE / "training_gene_families.tsv"
FIXED = HERE / "data" / "allx_isoform_fix15_2026-11-14.tsv"
LIT = HERE / "data" / "lit_ligand_contacts_2026-11-14.tsv"
VETO = (2, 3)  # residues n-2 .. n+3 for an insertion after n (contact_veto.VETO_HALF)
BOOT = 2000


# --------------------------------------------------------------------------- data
def fix_isoform_positives(pos, allx, ev):
    """Replace the misnumbered proteins' rows with the proteoform-fixed export."""
    fx = pd.read_csv(FIXED, sep="\t", low_memory=False)
    for c in ("ss3", "topology_class", "region_class"):
        fx[c] = fx[c].fillna("").astype(str)
    fx["gene_symbol"] = fx["gene_symbol"].fillna(fx["uniprot_acc"]).astype(str)
    accs = set(fx["uniprot_acc"]) & set(pos["uniprot_acc"])
    allx = pd.concat([allx[~allx["uniprot_acc"].isin(accs)], fx[allx.columns]], ignore_index=True)
    ev_cols = [c for c in ev.columns if c in fx.columns]
    ev = pd.concat(
        [ev[~ev["uniprot_acc"].isin(accs)], fx[fx["uniprot_acc"].isin(accs)][ev_cols]],
        ignore_index=True,
    )
    # only the model's inputs and the region index move; labels and ids stay the positives' own
    feat_cols = [c for c in (*R.FEATS, "region_index") if c in pos.columns and c in fx.columns]
    key = fx.set_index(["uniprot_acc", "residue"])
    n_fixed = 0
    for i in pos.index[pos["uniprot_acc"].isin(accs)]:
        k = (pos.at[i, "uniprot_acc"], int(pos.at[i, "residue"]))
        if k in key.index:
            for c in feat_cols:
                pos.at[i, c] = key.at[k, c]
            n_fixed += 1
    return pos, allx.reset_index(drop=True), ev.reset_index(drop=True), sorted(accs), n_fixed


def literature_on_published(ev, structure_only: bool) -> set[tuple[str, int]]:
    """Verified literature ligand contacts mapped onto published proteins (identical sequence)."""
    lit = pd.read_csv(LIT, sep="\t")
    lit = lit[lit["citation_verified"].astype(str).str.lower() == "true"]
    lit = lit[lit["contact_class"].isin({"small_ligand", "peptide_ligand", "protein_ligand"})]
    if structure_only:
        lit = lit[lit["evidence"] == "structure"]
    t = pd.read_csv(HERE / "l1_adopted_targets.tsv", sep="\t")
    acc_of = dict(zip(t["sequence_sha256"], t["uniprot_acc"], strict=False))
    keys = {
        (acc_of[s], int(r))
        for s, r in zip(lit["sequence_sha256"], lit["residue"], strict=False)
        if s in acc_of
    }
    return keys & set(zip(ev["uniprot_acc"], ev["residue"].astype(int), strict=False))


def families(genes) -> dict[str, str]:
    """HGNC gene group per gene symbol (first listed), cached in training_gene_families.tsv."""
    have = (
        pd.read_csv(FAMILIES, sep="\t").set_index("gene_symbol")["family"].to_dict()
        if FAMILIES.is_file()
        else {}
    )
    for g in sorted(set(genes) - set(have)):
        req = urllib.request.Request(
            f"https://rest.genenames.org/fetch/symbol/{g}",
            headers={"Accept": "application/json", "User-Agent": "deliverome-analysis s5e13"},
        )
        docs = json.loads(urllib.request.urlopen(req, timeout=60).read())["response"]["docs"]
        grp = (docs[0].get("gene_group") or [None])[0] if docs else None
        have[g] = grp or f"(no HGNC group) {g}"
        time.sleep(0.15)
    pd.DataFrame(sorted(have.items()), columns=["gene_symbol", "family"]).to_csv(
        FAMILIES, sep="\t", index=False
    )
    return have


def family_folds(base, fam) -> dict[str, int]:
    """gene -> fold, with every gene of an HGNC family in the same fold."""
    y = base["label"].to_numpy()
    groups = base["gene_symbol"].map(fam).fillna(base["gene_symbol"]).to_numpy()
    out: dict[str, int] = {}
    for i, (_a, b) in enumerate(GroupKFold(n_splits=R.N_SPLITS).split(base, y, groups)):
        for g in set(base["gene_symbol"].to_numpy()[b]):
            out[g] = i
    return out


# --------------------------------------------------------------------------- scoring
def weighted_frame(base, pos, allx, extra: np.ndarray, w_d1: float, w_x: float):
    """v5's upweighting with a second class of hard negatives at its own weight."""
    hn_any = allx["hard_negative"].to_numpy() | extra
    a2 = allx.assign(hard_negative=hn_any)
    tr, is_hn, _ = L.upweight_frame(base, pos, a2)
    keyx = set(
        zip(allx.loc[extra, "uniprot_acc"], allx.loc[extra, "residue"].astype(int), strict=False)
    )
    keyd = set(
        zip(
            allx.loc[allx["hard_negative"], "uniprot_acc"],
            allx.loc[allx["hard_negative"], "residue"].astype(int),
            strict=False,
        )
    )
    k = list(zip(tr["uniprot_acc"], tr["residue"].astype(int), strict=False))
    w = np.ones(len(tr))
    for i, kk in enumerate(k):
        if is_hn[i]:
            w[i] = w_d1 if kk in keyd else (w_x if kk in keyx else w_d1)
    return tr, w


def vetoed(ev, flag: np.ndarray) -> np.ndarray:
    """True where a flagged residue lies in n-2 .. n+3 of the junction after residue n."""
    out = np.zeros(len(ev), bool)
    res = ev["residue"].astype(int).to_numpy()
    for _acc, idx in ev.groupby("uniprot_acc").indices.items():
        bad = res[idx][flag[idx]]
        if len(bad):
            d = bad[None, :] - res[idx][:, None]
            out[idx] = ((d >= -VETO[0]) & (d <= VETO[1])).any(1)
    return out


def floors(ev, seed: int) -> dict[str, np.ndarray]:
    rng = np.random.default_rng(seed)
    rand = rng.random(len(ev))
    mid = np.zeros(len(ev))
    for _acc, g in ev.groupby("uniprot_acc"):
        sizes = g.groupby("region_index")["residue"].agg(["size", "min", "max"])
        best = sizes["size"].idxmax()
        m = (sizes.loc[best, "min"] + sizes.loc[best, "max"]) / 2
        mid[g.index] = -np.abs(g["residue"].to_numpy() - m) + 1000 * (
            g["region_index"].to_numpy() == best
        )
    return {"floor_random": rand, "floor_longest_loop_mid": mid}


def per_protein(ev, sc, all_sites, prot_class):
    prot, _ = R.per_protein(ev, sc, all_sites, all_sites, prot_class)
    reg = ev.set_index(["uniprot_acc", "residue"])["region_index"]
    hits = []
    for r in prot.itertuples(index=False):
        picks = [int(p) for p in str(r.picks).split(",") if p]
        sites = [int(s) for s in str(r.sites).split(",") if s]
        host = {reg.get((r.uniprot_acc, s)) for s in sites} - {None}
        pr = [reg.get((r.uniprot_acc, p)) for p in picks]
        hits.append((bool(pr and pr[0] in host), bool(any(p in host for p in pr[:3]))))
    prot["region_top1"] = [h[0] for h in hits]
    prot["region_top3"] = [h[1] for h in hits]
    return prot


METRICS = ("within10_top1", "within10_top3", "region_top1", "region_top3")


def bootstrap(a: pd.DataFrame, b: pd.DataFrame, m: str, seed=0) -> tuple[float, float, float]:
    """Mean difference a-b over proteins with a 95% bootstrap interval (paired by protein)."""
    j = a.set_index("uniprot_acc")[m].astype(float) - b.set_index("uniprot_acc")[m].astype(float)
    j = j.dropna().to_numpy()
    rng = np.random.default_rng(seed)
    bs = [rng.choice(j, len(j)).mean() for _ in range(BOOT)]
    return float(j.mean()), float(np.quantile(bs, 0.025)), float(np.quantile(bs, 0.975))


# --------------------------------------------------------------------------- main
def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true", help="one draw, gene folds, w_d1 5 only")
    args = ap.parse_args(argv)
    draws = DRAWS[:1] if args.quick else DRAWS
    weights = (5.0,) if args.quick else (5.0, 10.0)
    schemes = ("gene",) if args.quick else ("gene", "family")

    L.INTERFACES = "all"  # flags ligand, partner and homomer classes
    pos, allx_raw, ev, _ = R.frame()
    pos, allx_raw, ev, fixed_accs, n_fixed = fix_isoform_positives(pos, allx_raw, ev)
    screens = R.load_screens()
    allx = L.constraints(allx_raw, screens)
    ev = L.constraints(ev, screens)
    # b0 = v5's classes only: the interface flags must not count as hard negatives there
    v5_cols = [c for c in allx.columns if c.startswith("hn_")]
    v5_cols = [c for c in v5_cols if c not in ("hn_partner_interface", "hn_homomer_interface")]
    for df in (allx, ev):
        df["hard_negative"] = df[v5_cols].any(axis=1)
    feats = [f for f in R.FEATS if f in pos.columns]
    prot_class = pos.groupby("uniprot_acc")["tm_helix_count"].first().map(R.tm_class).to_dict()
    all_sites = pos.groupby("uniprot_acc")["residue"].apply(lambda v: {int(x) for x in v}).to_dict()

    def keyset(df):
        return list(zip(df["uniprot_acc"], df["residue"].astype(int), strict=False))

    lit_s, lit_a = literature_on_published(ev, True), literature_on_published(ev, False)
    extra_ev = {
        "b0": np.zeros(len(ev), bool),
        "b1_lit_structure": np.array([k in lit_s for k in keyset(ev)]),
        "b2_lit_all": np.array([k in lit_a for k in keyset(ev)]),
        "c_partner": ev["hn_partner_interface"].to_numpy(),
        "d_partner_homomer": (ev["hn_partner_interface"] | ev["hn_homomer_interface"]).to_numpy(),
    }
    extra_ev["e_lit_structure_partner_homomer"] = (
        extra_ev["b1_lit_structure"] | extra_ev["d_partner_homomer"]
    )
    extra_allx = {
        "b0": np.zeros(len(allx), bool),
        "b1_lit_structure": np.array([k in lit_s for k in keyset(allx)]),
        "b2_lit_all": np.array([k in lit_a for k in keyset(allx)]),
        "c_partner": allx["hn_partner_interface"].to_numpy(),
        "d_partner_homomer": (
            allx["hn_partner_interface"] | allx["hn_homomer_interface"]
        ).to_numpy(),
    }
    extra_allx["e_lit_structure_partner_homomer"] = (
        extra_allx["b1_lit_structure"] | extra_allx["d_partner_homomer"]
    )
    # how many published sites each class would veto (n-2 .. n+3), as the coordinator asked
    pub_keys = {(a, r) for a, rs in all_sites.items() for r in rs}
    is_pub = np.array([k in pub_keys for k in keyset(ev)])
    classes_ev = {
        "d1_ligand (in b0)": ev["hn_ligand_contact"].to_numpy(),
        **{k: v for k, v in extra_ev.items() if k != "b0"},
    }
    published_vetoed = {k: int((vetoed(ev, v) & is_pub).sum()) for k, v in classes_ev.items()}
    reach = {
        k: len({a for a, f in zip(ev["uniprot_acc"], v, strict=False) if f})
        for k, v in extra_ev.items()
    }
    fam = families(pos["gene_symbol"].unique())

    rows, prot_rows = [], []
    for scheme in schemes:
        for draw in draws:
            base = build_global(pos, allx_raw, np.random.default_rng(draw), R.RATIO)
            fog = (
                gene_folds(build_global(pos, allx_raw, np.random.default_rng(R.SEED), R.RATIO))
                if scheme == "gene"
                else family_folds(base, fam)
            )
            arms = {}
            for wd in weights:
                for name in extra_allx:
                    t0 = time.time()
                    tr, w = weighted_frame(base, pos, allx, extra_allx[name], wd, wd / 2)
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        sc, _ = L.oof_weighted(tr, w, ev, feats, fog)
                    arms[f"{name} train w{int(wd)}"] = sc
                    if name != "b0":
                        b0 = arms[f"b0 train w{int(wd)}"]
                        arms[f"{name} veto-only w{int(wd)}"] = np.where(
                            vetoed(ev, extra_ev[name]), -1.0, b0
                        )
                    print(scheme, draw, name, wd, f"{time.time() - t0:.0f}s", flush=True)
            arms.update(floors(ev, draw))
            for arm, sc in arms.items():
                prot = per_protein(ev, sc, all_sites, prot_class)
                prot_rows.append(prot.assign(scheme=scheme, draw=draw, arm=arm))
                rows.append(
                    {"scheme": scheme, "draw": draw, "arm": arm, "proteins": len(prot)}
                    | {m: round(float(prot[m].mean()), 4) for m in METRICS}
                )
    t = pd.DataFrame(rows)
    t.to_csv(OUT, sep="\t", index=False)
    prot = pd.concat(prot_rows, ignore_index=True)
    prot.to_csv(OUT_PROT, sep="\t", index=False)

    summary = t.groupby(["scheme", "arm"], sort=False)[list(METRICS)].agg(["mean", "std"]).round(4)
    diffs = []
    for scheme in schemes:
        for arm in t["arm"].unique():
            wd = arm.rsplit(" w", 1)[-1] if " w" in arm else "5"
            ref = f"b0 train w{wd}"
            if arm == ref:
                continue
            for m in METRICS:
                ds = []
                for draw in draws:
                    q = prot[(prot["scheme"] == scheme) & (prot["draw"] == draw)]
                    ds.append(bootstrap(q[q["arm"] == arm], q[q["arm"] == ref], m, seed=draw))
                d = np.array(ds)
                diffs.append(
                    {
                        "scheme": scheme,
                        "arm": arm,
                        "vs": ref,
                        "metric": m,
                        "mean_diff": round(float(d[:, 0].mean()), 4),
                        "ci_low": round(float(d[:, 1].mean()), 4),
                        "ci_high": round(float(d[:, 2].mean()), 4),
                        "clearly_hurts": bool(d[:, 2].mean() < 0),
                        "clearly_helps": bool(d[:, 1].mean() > 0),
                    }
                )
    dd = pd.DataFrame(diffs)
    dd.to_csv(HERE / "hard_negative_benchmark_diffs.tsv", sep="\t", index=False)
    OUT_JSON.write_text(
        json.dumps(
            {
                "decision_rule": "adopted unless the 95% bootstrap interval of its top-3 or "
                "region-hit difference from b0 lies wholly below zero (coordinator under Becca's "
                "delegation, 2026-11-15); weights chosen out of fold, reported",
                "positives": len(pos),
                "proteins": int(pos["uniprot_acc"].nunique()),
                "fixed_isoform_proteins": fixed_accs,
                "fixed_positive_rows": n_fixed,
                "draws": list(draws),
                "weights_d1": list(weights),
                "weight_extra": "half of the D1 weight",
                "veto_window": "n-2 .. n+3",
                "folds": {
                    "gene": "GroupKFold-5 by gene (v4's map)",
                    "family": "GroupKFold-5 by HGNC gene group (training_gene_families.tsv)",
                },
                "published_proteins_reached": reach,
                "published_sites_vetoed_n-2..n+3": published_vetoed,
                "published_sites": len(pub_keys),
                "bootstrap_resamples": BOOT,
                "n_families": len(set(fam.values())),
            },
            indent=2,
        )
        + "\n"
    )
    print(summary.to_string())
    print(dd[dd["metric"].isin(["within10_top3", "region_top3"])].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
