#!/usr/bin/env python3
"""Is the model better at inter-TM loops than at disordered ectodomain stems?

The single-pass result -- top-1 0.233 against 0.463 for 7+ TM proteins -- has two
candidate explanations that the topology split cannot separate. Either long
ectodomains are simply a bigger haystack, or the model is genuinely better at a
constrained inter-TM loop than at a long disordered stem where many residues look
alike.

`region_class` separates them: every extracellular residue is already labelled
`disordered` or `ordered_surface` by the metapredict cut the feature export
applies. Scoring the published cohort out-of-fold once and splitting the SAME
scores by that label answers it directly, and splitting again by topology shows
whether the two explanations are the same explanation.

    uv run --locked python region_class_strata.py
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[3] / "src"))

from final_report import fold_eval, mask_keep  # noqa: E402
from independent_model_search import load  # noqa: E402
from tag_site_model_v4 import FEATS, KIND, training_frame  # noqa: E402


def topo(tm) -> str:
    if pd.isna(tm):
        return "unknown"
    tm = int(tm)
    return "single-pass" if tm <= 1 else ("multipass 2-6" if tm <= 6 else "multipass 7+")


def score_block(g: pd.DataFrame) -> dict | None:
    if len(np.unique(g["y"])) < 2:
        return None
    hits, n_prot = 0, 0
    for _, gg in g.groupby("acc"):
        if gg["y"].sum() == 0:
            continue
        n_prot += 1
        order = gg.sort_values("s", ascending=False)
        hits += int(np.where(order["y"].to_numpy() == 1)[0].min()) + 1 <= 1
    p = g[g["y"] == 1]
    return {
        "proteins": n_prot,
        "sites": int(g["y"].sum()),
        "residues_per_protein": round(len(g) / max(n_prot, 1)),
        "auroc": round(roc_auc_score(g["y"], g["s"]), 4),
        "top1": round(hits / max(n_prot, 1), 3),
        "median_disorder": round(float(p["disorder"].median()), 3),
        "median_rsa": round(float(p["rsa"].median()), 3),
        "median_loop_len": round(float(p["loop_len"].median()), 1),
    }


def main() -> int:
    _df, ev, ev_y = load()
    _p, _a, train = training_frame()
    keep = mask_keep(ev, ev_y)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc, _fold = fold_eval(train, ev, ev_y, FEATS, KIND)
    ok = keep & ~np.isnan(sc)
    d = pd.DataFrame(
        {
            "acc": ev["uniprot_acc"],
            "y": ev_y,
            "s": sc,
            "region_class": ev["region_class"].fillna("").astype(str),
            "topology": ev["tm_helix_count"].map(topo),
            "disorder": ev["disorder"],
            "rsa": ev["rsa"],
            "loop_len": ev["loop_len"],
        }
    )[ok]
    d = d[d["region_class"] != ""]

    rows = []
    for name, g in [*d.groupby("region_class"), ("all", d)]:
        r = score_block(g)
        if r:
            rows.append({"split": "region class", "stratum": name, **r})
    for (rc, tp), g in d.groupby(["region_class", "topology"]):
        r = score_block(g)
        if r and r["proteins"] >= 8:
            rows.append({"split": "both", "stratum": f"{rc} / {tp}", **r})
    out = pd.DataFrame(rows)
    out.to_csv(HERE / "region_class_strata.tsv", sep="\t", index=False)
    print(out.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
