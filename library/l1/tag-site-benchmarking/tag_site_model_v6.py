#!/usr/bin/env python3
"""Tag-site model v6: v5 with the hard negatives Becca asked to consider, and the positives fixed.

Becca, 2026-11-15: "even if v5 doesn't do better with lit it is still good to consider adding
hard negatives given the sparse data we are training on" and "we want to consider ligand and co
receptors / homomer contacts as well". The coordinator, under her delegation: adopt the
retrained model on hypothesis grounds unless it clearly hurts. hard_negative_benchmark.py found
that no variant clearly hurts or helps (bootstrap over proteins, gene- and family-grouped folds),
so v6 trains on all of them (variant e):

    D1 hard negatives (v5's: buried core, disulfide Cys, glycan/sequon N, screen-intolerant,
    D1 ligand contacts)                                           weight 5 (chosen out of fold)
    + literature ligand contacts with structure evidence (PR #255)  weight 2.5
    + D1 protein-partner interface residues                         weight 2.5
    + homomer interface residues (Schweke assembly models)          weight 2.5

and the 5 positives on 4 proteins that v5 had numbered on another isoform are re-featured from
the proteoform-fixed export. Same features, estimator, seed and negative draw as v5.

    uv run --locked python tag_site_model_v6.py   # writes tag_site_model_v6_{rank,absolute}.joblib
"""

from __future__ import annotations

import json
import sys
import warnings
from datetime import UTC, datetime
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[3] / "src"))

import hard_negative_benchmark as H  # noqa: E402
import ranker_hard_negative_layer as L  # noqa: E402
import ranker_pu_segment as R  # noqa: E402
import tag_site_model_v4 as V4  # noqa: E402
import tag_site_model_v5 as V5  # noqa: E402
from independent_model_search import make, within_protein  # noqa: E402
from two_tier_model import build_global  # noqa: E402

W_D1, W_EXTRA = 5.0, 2.5
MODEL_Z = HERE / "tag_site_model_v6_rank.joblib"
MODEL_RAW = HERE / "tag_site_model_v6_absolute.joblib"
SPEC = HERE / "tag_site_model_v6.json"
MIN_REGION_LEN = V5.MIN_REGION_LEN


def training():
    """v5's frame with the fixed positives, plus the extra hard-negative classes at W_EXTRA."""
    L.INTERFACES = "all"
    pos, allx_raw, ev, _ = R.frame()
    pos, allx_raw, ev, fixed_accs, n_fixed = H.fix_isoform_positives(pos, allx_raw, ev)
    allx = L.constraints(allx_raw, R.load_screens())
    v5_cols = [
        c
        for c in allx.columns
        if c.startswith("hn_") and c not in ("hn_partner_interface", "hn_homomer_interface")
    ]
    allx["hard_negative"] = allx[v5_cols].any(axis=1)
    keys = list(zip(allx["uniprot_acc"], allx["residue"].astype(int), strict=False))
    lit = H.literature_on_published(ev, structure_only=True)
    extra = (
        np.array([k in lit for k in keys])
        | allx["hn_partner_interface"].to_numpy()
        | allx["hn_homomer_interface"].to_numpy()
    )
    base = build_global(pos, allx_raw, np.random.default_rng(V4.SEED), V4.RATIO)
    tr, w = H.weighted_frame(base, pos, allx, extra, W_D1, W_EXTRA)
    counts = {
        **{c: int(allx.loc[allx["hard_negative"], c].sum()) for c in v5_cols},
        "literature_structure_ligand": int(extra.sum() and sum(k in lit for k in keys)),
        "hn_partner_interface": int(allx["hn_partner_interface"].sum()),
        "hn_homomer_interface": int(allx["hn_homomer_interface"].sum()),
    }
    meta = {
        "fixed_isoform_proteins": fixed_accs,
        "fixed_positive_rows": n_fixed,
        "n_weighted_d1": int((w == W_D1).sum()),
        "n_weighted_extra": int((w == W_EXTRA).sum()),
        "hard_negative_residues_on_training_proteins_by_class": counts,
    }
    return pos, allx_raw, tr, w, meta


def fit():
    _pos, allx, tr, w, meta = training()
    feats = [f for f in V4.FEATS if f in tr.columns]
    y = tr["label"].to_numpy()
    mz, mraw = make(V4.KIND, feats), make(V4.KIND, feats)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        mz.fit(within_protein(tr, feats, "z"), y, clf__sample_weight=w)
        mraw.fit(tr, y, clf__sample_weight=w)
    pub_oof, pub_ceilings, pub_regions = V5.published_oof_weighted(tr, w, allx, feats)
    import sklearn

    q = ("0.25", "0.5", "0.75")
    spec = {
        "fitted_at": datetime.now(UTC).isoformat(),
        "git_commit": V4.git_sha(),
        "supersedes": "tag_site_model_v5 (same features, estimator, seed and draw)",
        "estimator": "RandomForestClassifier(n_estimators=800, max_depth=4, "
        "min_samples_leaf=12, class_weight='balanced_subsample'), fitted with sample weights",
        "features": feats,
        "n_features": len(feats),
        "n_train_rows": len(tr),
        "hard_negatives": "v5's classes at weight 5; literature structure-evidence ligand "
        "contacts, D1 protein-partner interfaces and Schweke homomer interfaces at weight 2.5",
        "weights": {"d1": W_D1, "extra": W_EXTRA},
        "chosen_by": "hard_negative_benchmark.py: no variant clearly hurts (95% bootstrap over "
        "proteins); adopted on Becca's hypothesis-driven standard (2026-11-15) via the "
        "coordinator; D1 weight 5 over 10 chosen out of fold on within-10 top-1",
        **meta,
        "seed": int(V4.SEED),
        "ratio": int(V4.RATIO),
        "min_region_len": MIN_REGION_LEN,
        "published_absolute_scores_oof": [round(float(v), 6) for v in pub_oof],
        "published_protein_ceilings_oof": [round(float(v), 6) for v in pub_ceilings],
        "published_region_ceilings_oof": [round(float(v), 6) for v in pub_regions],
        "published_region_ceiling_quantiles_oof": {
            k: round(float(np.quantile(pub_regions, float(k))), 4) for k in q
        },
        "sklearn_version": sklearn.__version__,
        "numpy_version": np.__version__,
        "decided_by": "coordinator under Becca's 2026-11-12 delegation, on her 2026-11-15 words",
        "implemented_by": "ranker-pu-segment session (Claude Code)",
    }
    return mz, mraw, spec


def main() -> int:
    import joblib

    mz, mraw, spec = fit()
    joblib.dump(mz, MODEL_Z)
    joblib.dump(mraw, MODEL_RAW)
    SPEC.write_text(json.dumps(spec, indent=2) + "\n")
    print(
        f"fitted v6 on {spec['n_train_rows']} rows; D1 hard negatives {spec['n_weighted_d1']}, "
        f"extra {spec['n_weighted_extra']}; region ceilings OOF "
        f"{spec['published_region_ceiling_quantiles_oof']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
