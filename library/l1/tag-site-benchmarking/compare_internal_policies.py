#!/usr/bin/env python3
"""Which internal-site policy: RF on positives only, PU-RF, or PU-RF with a contact veto?

Becca, 2026-11-12: "For internal sites we will want to consider both random forest positive
prediction (maybe with PU) plus hard negatives like the literature agent ligand sites."
This scores five policies on the same gene -> fold map, out of fold, on three draws:

    rf_positive_only   v4's frame (published sites against a random draw of unlabelled
                       extracellular residues), through the weighted fit with every weight 1
    rf_hn_w5 (v5)      the same frame plus the ligand-inclusive hard negatives at weight 5,
                       the shipped model
    pu_rf              bagged PU (Mordelet & Vert 2014): 50 forests, each on the positives
                       and a fresh draw of one unlabelled residue per positive
    pu_rf_veto         pu_rf, then any residue within 2 of a ligand contact is removed
                       before picking (the scores are smoothed over +/-2)
    rf_hn_w5_veto      v5 with the same veto

Ligand contacts on the published proteins are D1's structure-derived extracellular ligand
contacts (`fetch_interface_contacts.py`). With `--literature PATH`, contacts the literature
agents (PR #255) verified are added to the veto; they only reach published proteins that
are also L1 members.

Two questions, two metric sets:

    tolerated sites   the 176 published sites (127 proteins) and the 111 exact-site subset:
                      within-10 top-1, within-10 top-3, exact top-1, right segment first
    contact sites     on published proteins with an extracellular ligand contact, held
                      out by gene: how often the first pick or any of the top 3 lands
                      within 2 of a contact, and the contacts' median within-protein
                      percentile (lower is better; a policy that avoids them ranks them low)

The veto arms remove the vetoed contacts by construction, so their contact numbers are not
evidence that the scorer avoids contacts. Compare the scorers on the contact metrics
without a veto, and the veto arms on what the veto costs at tolerated sites.

    uv run --locked python compare_internal_policies.py   # writes compare_internal_policies.*
"""

from __future__ import annotations

import argparse
import json
import sys
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
from independent_model_search import smooth, within_protein  # noqa: E402
from two_tier_model import build_global  # noqa: E402

DRAWS = (R.SEED, 1, 2)  # v4's negative draw, then two more (as refit_corrected_hard_negative.py)
HN_WEIGHT = 5
VETO_HALF = L.VETO_HALF
OUT = HERE / "compare_internal_policies.tsv"
OUT_PROT = HERE / "compare_internal_policies_per_protein.tsv"
OUT_JSON = HERE / "compare_internal_policies.json"


def oof_pu_seeded(pos, allx, ev, feats, fog, seed):
    """R.oof_pu with an explicit seed, so the PU arm can be repeated on three draws."""
    full = R.pu_frame(pos, allx)
    tr_z = within_protein(full, feats, "z")
    y = full["label"].to_numpy()
    ed = within_protein(ev, feats, "z")
    gg = ev["gene_symbol"].to_numpy()
    sc = np.full(len(ev), np.nan)
    rng = np.random.default_rng(seed)
    for a, b in R.fixed_splits(full["gene_symbol"].to_numpy(), fog):
        mask = np.isin(gg, list(set(full.iloc[b]["gene_symbol"])))
        (p,) = R.pu_fit_predict(tr_z.iloc[a].reset_index(drop=True), y[a], [ed[mask]], feats, rng)
        sc[mask] = p
    return smooth(ev, sc, 2)


def literature_keys(path: str | None, ev: pd.DataFrame) -> set[tuple[str, int]]:
    """Verified literature contacts on published proteins, as (accession, residue)."""
    if not path:
        return set()
    lit = pd.read_csv(path, sep="\t")
    lit = lit[lit["citation_verified"].astype(str).str.lower() == "true"]
    # the literature table is numbered on L1 ordered proteins; it reaches a published
    # protein only where that protein's ordered sequence is its canonical (same sha)
    t = pd.read_csv(HERE / "l1_adopted_targets.tsv", sep="\t")
    acc_of = dict(zip(t["sequence_sha256"], t["uniprot_acc"], strict=False))
    keys = {
        (acc_of[s], int(r))
        for s, r in zip(lit["sequence_sha256"], lit["residue"], strict=False)
        if s in acc_of
    }
    present = set(zip(ev["uniprot_acc"], ev["residue"].astype(int), strict=False))
    return keys & present


def contact_metrics(ev, sc, contact: np.ndarray) -> dict:
    d = ev[["uniprot_acc", "residue"]].assign(s=sc, c=contact)
    top1, top3, pct, n = 0, 0, [], 0
    for _acc, g in d.dropna(subset=["s"]).groupby("uniprot_acc"):
        if not g["c"].any():
            continue
        n += 1
        res = g["residue"].to_numpy().astype(int)
        s = g["s"].to_numpy()
        bad = res[g["c"].to_numpy()]
        picks = R.separated_top(res, s)
        top1 += bool(np.abs(picks[0] - bad).min() <= VETO_HALF)
        top3 += any(np.abs(p - bad).min() <= VETO_HALF for p in picks)
        order = s.argsort().argsort() / max(len(s) - 1, 1) * 100
        pct.append(float(np.median(order[g["c"].to_numpy()])))
    return {
        "contact_proteins": n,
        "contact_top1_within2": round(top1 / n, 3) if n else np.nan,
        "contact_top3_any_within2": round(top3 / n, 3) if n else np.nan,
        "contact_median_percentile": round(float(np.median(pct)), 1) if pct else np.nan,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--literature", help="verified literature contacts TSV (PR #255)")
    args = ap.parse_args()

    L.INTERFACES = "ligands"
    pos, allx_raw, ev, _ev_y = R.frame()
    screens = R.load_screens()
    allx = L.constraints(allx_raw, screens)
    ev = L.constraints(ev, screens)
    feats = [f for f in R.FEATS if f in pos.columns]
    prot_class = pos.groupby("uniprot_acc")["tm_helix_count"].first().map(R.tm_class).to_dict()
    all_sites = pos.groupby("uniprot_acc")["residue"].apply(lambda v: {int(x) for x in v}).to_dict()
    ex = pos[pos["evidence_type"] == R.EXACT]
    exact_sites = (
        ex.groupby("uniprot_acc")["residue"].apply(lambda v: {int(x) for x in v}).to_dict()
    )

    lit = literature_keys(args.literature, ev)
    key = list(zip(ev["uniprot_acc"], ev["residue"].astype(int), strict=False))
    contact = ev["hn_ligand_contact"].to_numpy() | np.array([k in lit for k in key])
    veto = L.near(ev, pd.Series(contact, index=ev.index), VETO_HALF)

    rows, prot_rows = [], []
    fog = None
    for draw in DRAWS:
        base = build_global(pos, allx_raw, np.random.default_rng(draw), R.RATIO)
        if fog is None:
            fog = gene_folds(base)  # fixed from v4's shipped draw, reused for every draw
        tr, is_hn, _ = L.upweight_frame(base, pos, allx)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            rf, fold = L.oof_weighted(base, np.ones(len(base)), ev, feats, fog)
            v5, _ = L.oof_weighted(tr, np.where(is_hn, float(HN_WEIGHT), 1.0), ev, feats, fog)
        pu = oof_pu_seeded(pos, allx_raw, ev, feats, fog, draw)
        arms = {
            "rf_positive_only": rf,
            "rf_hn_w5 (v5)": v5,
            "pu_rf": pu,
            "pu_rf_veto": np.where(veto, -1.0, pu),
            "rf_hn_w5_veto": np.where(veto, -1.0, v5),
        }
        for arm, sc in arms.items():
            cm = contact_metrics(ev, sc, contact)
            for cohort, truth in (("all 176", all_sites), ("exact-site 111", exact_sites)):
                prot, _sites = R.per_protein(ev, sc, truth, all_sites, prot_class)
                ex1 = L.exact_top1(prot)
                rows.append(
                    {
                        "draw": draw,
                        "arm": arm,
                        "cohort": cohort,
                        "proteins": len(prot),
                        "within10_top1": round(float(prot["within10_top1"].mean()), 3),
                        "within10_top3": round(float(prot["within10_top3"].mean()), 3),
                        "exact_top1": round(float(ex1.mean()), 3),
                        "segment_first": round(float(prot["segment_first"].dropna().mean()), 3),
                        **cm,
                    }
                )
                prot_rows.append(prot.assign(draw=draw, arm=arm, cohort=cohort, exact_top1=ex1))
            print(draw, arm, rows[-2], flush=True)
        del fold
    t = pd.DataFrame(rows)
    t.to_csv(OUT, sep="\t", index=False)
    prot = pd.concat(prot_rows)
    prot.to_csv(OUT_PROT, sep="\t", index=False)

    mean = (
        t.drop(columns="draw")
        .groupby(["cohort", "arm"], sort=False)
        .mean(numeric_only=True)
        .round(3)
    )
    paired = []
    for cohort in ("all 176", "exact-site 111"):
        for draw in DRAWS:
            q = prot[(prot["cohort"] == cohort) & (prot["draw"] == draw)]
            ref = q[q["arm"] == "rf_hn_w5 (v5)"].set_index("uniprot_acc")
            for arm in t["arm"].unique():
                if arm == "rf_hn_w5 (v5)":
                    continue
                a = q[q["arm"] == arm].set_index("uniprot_acc")
                r = {"cohort": cohort, "draw": draw, "arm": arm}
                for m in ("within10_top1", "within10_top3"):
                    x, y = a[m].astype(float), ref[m].astype(float)
                    r[f"{m}_gained_vs_v5"] = int((x > y).sum())
                    r[f"{m}_lost_vs_v5"] = int((x < y).sum())
                paired.append(r)
    OUT_JSON.write_text(
        json.dumps(
            {
                "estimator": R.KIND,
                "features": feats,
                "draws": list(DRAWS),
                "cv": "GroupKFold by gene, 5 folds, gene->fold fixed from the v4 draw",
                "pu_bags": R.PU_BAGS,
                "hn_weight": HN_WEIGHT,
                "veto_half": VETO_HALF,
                "contact_residues_on_published_proteins": int(contact.sum()),
                "literature_contacts_on_published_proteins": len(lit),
                "residues_vetoed": int(veto.sum()),
                "mean_over_draws": {
                    f"{c} | {a}": v for (c, a), v in mean.to_dict(orient="index").items()
                },
                "paired_vs_v5": paired,
                "decided_by": "ranker-pu-segment (delegated by Becca 2026-11-12)",
                "implemented_by": "ranker-pu-segment session (Claude Code), PR #245",
            },
            indent=2,
        )
        + "\n"
    )
    print(mean.to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
