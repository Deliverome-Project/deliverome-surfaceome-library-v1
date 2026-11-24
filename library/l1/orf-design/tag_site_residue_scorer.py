#!/usr/bin/env python3
"""Residue-level tag-site scorer: continuous tie-breaking, ProDomino/PLM, and literature contacts.

Tests two questions before re-predicting on full L1:
  1. Depletion ratios (RR_c) for every contact class in `contacts_by_protein.tsv`,
     explicitly separating literature-derived contacts (`literature_small_ligand`,
     `literature_homomer`, `literature_partner`, `literature_protein_ligand`,
     `literature_peptide_ligand`) from structure-derived contacts (`ligand`,
     `protein_ligand`, `protein_partner`, and `antibody` as an exposed-loop control).
  2. Whether blending the z-score-fixed RF with continuous residue-level signals
     (available on ~96% of L1 today: `-plddt_w5`, `rsa_w5`, plus `ProDomino` where
     available) and Empirical-Bayes penalties improves BOTH fine residue-level
     precision (exact, +/-3, +/-5) AND loop-level recall (+/-10), on gene and family
     splits (unweighted + L1-weighted + exposure).

Outputs:
    tag-site-models/residue_scorer.json  (committed summary record)
    tag-site-models/residue_scorer.tsv   (local gitignored table)

    uv run --locked python library/l1/orf-design/tag_site_residue_scorer.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import tag_site_fine_tolerance as F  # noqa: E402
import tag_site_model_ablations as A  # noqa: E402
import tag_site_model_benchmark as B  # noqa: E402
import tag_site_negative_calibration as N  # noqa: E402

OUT = B.OUT
TOLS = (0, 3, 5, 10)
CONTACT_CLASSES = (
    "ligand",
    "protein_ligand",
    "protein_partner",
    "literature_small_ligand",
    "literature_protein_ligand",
    "literature_peptide_ligand",
    "literature_partner",
    "literature_homomer",
    "antibody",
)


def load_literature_and_structure_contacts(
    ev: pd.DataFrame,
) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """Map contacts_by_protein.tsv onto ev (via canonical_acc or sequence_sha256->uniprot_acc)."""
    cpath = B.MODELS / "data" / "contacts_by_protein.tsv"
    c = pd.read_csv(
        cpath,
        sep="\t",
        usecols=["sequence_sha256", "canonical_acc", "residue", "class", "source"],
        low_memory=False,
    )
    tgt = pd.read_csv(B.MODELS / "l1_adopted_targets.tsv", sep="\t")
    sha_to_acc = dict(zip(tgt["sequence_sha256"], tgt["uniprot_acc"], strict=False))
    c["acc"] = c["canonical_acc"].fillna(c["sequence_sha256"].map(sha_to_acc))
    c = c.dropna(subset=["acc", "residue"]).copy()
    c["residue"] = c["residue"].astype(int)
    covered_accs = set(c["acc"].unique())
    cov_mask = ev["uniprot_acc"].isin(covered_accs).to_numpy()
    ev_keys = list(zip(ev["uniprot_acc"], ev["residue"].astype(int), strict=True))

    out: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for cls_name in CONTACT_CLASSES:
        sub = c[c["class"] == cls_name]
        kset = set(zip(sub["acc"], sub["residue"], strict=True))
        raw = np.array([k in kset for k in ev_keys], dtype=bool)
        win = N.window_mask(ev, raw)
        out[cls_name] = (win, cov_mask)
    return out


def main() -> int:
    tr, allx = B.inputs()
    gold, regions, audit = B.truth()
    gold_sites = gold.groupby("uniprot_acc")["residue"].apply(lambda v: np.array(v)).to_dict()

    ev = allx[allx["uniprot_acc"].isin(gold_sites)].reset_index(drop=True)
    have = set(zip(ev["uniprot_acc"], ev["residue"], strict=True))
    gold_in = gold[
        [(a, r) in have for a, r in zip(gold["uniprot_acc"], gold["residue"], strict=True)]
    ]
    gold_sites = gold_in.groupby("uniprot_acc")["residue"].apply(lambda v: np.array(v)).to_dict()
    ev = ev[ev["uniprot_acc"].isin(gold_sites)].reset_index(drop=True)

    key = ev.set_index(["uniprot_acc", "residue"]).index
    on_gold = key.isin(list(zip(gold_in["uniprot_acc"], gold_in["residue"], strict=True)))

    # 1. Depletion ratios for literature-derived and structure-derived protein contacts
    contact_classes = load_literature_and_structure_contacts(ev)
    dep_contacts = N.compute_depletion(ev, on_gold, contact_classes)
    print("\n=== CONTACT DEPLETION ON CONFIRMED SITES (contacts_by_protein.tsv, +/-2..+3) ===")
    print(dep_contacts.to_string(index=False))

    # 2. Structural/protrusion/seam EB penalties from negative_calibration
    ev_h, neg_classes = N.build_negative_classes(ev)
    dep_neg = N.compute_depletion(ev_h, on_gold, neg_classes)
    hard_mask = np.zeros(len(ev), bool)
    for cname in dep_neg.loc[dep_neg["recommended_tier"] == "hard_veto", "negative_class"]:
        hard_mask |= neg_classes[cname][0]
    soft_pen = np.zeros(len(ev), float)
    for r in dep_neg.loc[dep_neg["recommended_tier"] == "soft_penalty"].itertuples(index=False):
        soft_pen += float(r.log_odds_penalty) * neg_classes[r.negative_class][0].astype(float)

    # Hard ligand + protein-ligand + literature_small_ligand/protein_ligand/peptide_ligand veto
    lit_lig_mask = (
        contact_classes["ligand"][0]
        | contact_classes["protein_ligand"][0]
        | contact_classes["literature_small_ligand"][0]
        | contact_classes["literature_protein_ligand"][0]
        | contact_classes["literature_peptide_ligand"][0]
    )
    lit_partner_mask = (
        contact_classes["literature_partner"][0] | contact_classes["literature_homomer"][0]
    )
    all_bind_mask = (
        lit_lig_mask
        | contact_classes["protein_partner"][0]
        | lit_partner_mask
        | neg_classes["homomer_interface_5A"][0]
    )
    eval_flags = {
        "lit_lig": lit_lig_mask,
        "lit_partner": lit_partner_mask,
        "all_binding": all_bind_mask,
        "seam_8A": neg_classes["domain_seam_8A"][0],
        "crowded_q75": neg_classes["crowded_protrusion_q75"][0],
    }

    # 3. Continuous residue-level blends on gene & family splits
    conf, w_base = B.confirmed_frame(tr, allx, gold, regions, audit)
    w5 = np.where(w_base > 1.0, 5.0, 1.0)
    cls = A.protein_class(allx)
    comp = A.l1_composition()
    genes = sorted(set(conf["gene_symbol"]) | set(ev["gene_symbol"]))
    splits = {sc: B.fold_map(genes, sc) for sc in ("gene", "family")}

    # Continuous signals (within-protein percentile ranks in [0, 1])
    plddt_w5 = -ev["plddt_w5"].to_numpy(float)
    rsa_w5 = ev["rsa_w5"].to_numpy(float)
    prodom = ev["prodomino"].to_numpy(float)

    r_plddt = A.rank_pct_within_protein(ev, plddt_w5)
    r_rsa = A.rank_pct_within_protein(ev, rsa_w5)
    r_prodom = A.rank_pct_within_protein(ev, prodom)

    # Convert EB soft penalty to rank-compatible scale (normalize logit penalty by ~4.0)
    eb_rank_pen = 0.08 * soft_pen
    lit_soft_pen = 0.05 * lit_partner_mask.astype(float)

    rows = []
    for sc, folds in splits.items():
        s_base = F.oof_custom(
            conf, w5, ev, allx, folds, B.FEATS, B.CATEGORICAL, "tr_z", "cross_loop", 2
        )
        s_allxz = F.oof_custom(
            conf, w5, ev, allx, folds, B.FEATS, B.CATEGORICAL, "allx_z", "contiguous", 2
        )
        r_rf = A.rank_pct_within_protein(ev, s_allxz)

        # L1-ready blends (ZERO ProDomino needed — uses only features on 96% of L1 today)
        b1_rf_plddt = np.where(np.isnan(r_plddt), r_rf, 0.50 * r_rf + 0.50 * r_plddt)
        b2_rf_plddt_rsa = np.where(
            np.isnan(r_plddt) | np.isnan(r_rsa),
            b1_rf_plddt,
            0.50 * r_rf + 0.35 * r_plddt + 0.15 * r_rsa,
        )
        b3_rf_plddt_eb = np.where(
            hard_mask | lit_lig_mask, -1e9, b1_rf_plddt - eb_rank_pen - lit_soft_pen
        )
        b4_rf_plddt_rsa_eb = np.where(
            hard_mask | lit_lig_mask, -1e9, b2_rf_plddt_rsa - eb_rank_pen - lit_soft_pen
        )

        # ProDomino blends (requires running ProDomino on remaining ~360 L1 internal CDSs)
        c1_plddt_prodom = np.where(np.isnan(r_prodom), r_plddt, 0.50 * r_plddt + 0.50 * r_prodom)
        c2_plddt_prodom_eb = np.where(
            hard_mask | lit_lig_mask, -1e9, c1_plddt_prodom - eb_rank_pen - lit_soft_pen
        )
        c3_3way = np.where(
            np.isnan(r_prodom),
            b1_rf_plddt,
            0.34 * r_rf + 0.33 * r_plddt + 0.33 * r_prodom,
        )
        c4_3way_eb = np.where(hard_mask | lit_lig_mask, -1e9, c3_3way - eb_rank_pen - lit_soft_pen)

        candidates = [
            ("A0_baseline_rf", s_base),
            ("A1_allxz_rf", s_allxz),
            ("B1_L1ready_rf_plddt", b1_rf_plddt),
            ("B2_L1ready_rf_plddt_rsa", b2_rf_plddt_rsa),
            ("B3_L1ready_rf_plddt_eb_litveto", b3_rf_plddt_eb),
            ("B4_L1ready_rf_plddt_rsa_eb_litveto", b4_rf_plddt_rsa_eb),
            ("C1_prodom_plddt", c1_plddt_prodom),
            ("C2_prodom_plddt_eb_litveto", c2_plddt_prodom_eb),
            ("C3_prodom_3way_rf_plddt", c3_3way),
            ("C4_prodom_3way_eb_litveto", c4_3way_eb),
        ]

        ref_pp = F.eval_multi_tol(ev, s_base, gold_sites).set_index("uniprot_acc")
        for name, score in candidates:
            pp = F.eval_multi_tol(ev, score, gold_sites).set_index("uniprot_acc")
            pp_b = B.per_protein(ev, score, gold_sites)
            wtd = A.l1_weighted_metrics(pp_b, cls, comp)
            exp = A.exposure_rates(ev, score, eval_flags)
            row: dict[str, object] = {"split": sc, "model": name, "n": len(pp)}
            for t in TOLS:
                row[f"t1_w{t}"] = round(float(pp[f"top1_w{t}"].mean()), 3)
                row[f"t3_w{t}"] = round(float(pp[f"top3_w{t}"].mean()), 3)
            row["l1w_t1_w10"] = wtd["l1w_top1_within10"]
            row["l1w_t3_w10"] = wtd["l1w_top3_within10"]
            # Paired sign test vs baseline at top1_w5 and top1_w10
            common = ref_pp.index.intersection(pp.index)
            for mcol in ("top1_w3", "top1_w5", "top1_w10", "top3_w5"):
                a = ref_pp.loc[common, mcol].astype(bool)
                b = pp.loc[common, mcol].astype(bool)
                wins, losses = int((~a & b).sum()), int((a & ~b).sum())
                row[f"wl_{mcol}"] = f"{wins}/{losses}"
                row[f"p_{mcol}"] = round(B.sign_test(wins, losses), 4)
            row.update(exp)
            rows.append(row)

    df = pd.DataFrame(rows)
    OUT.mkdir(exist_ok=True)
    df.to_csv(OUT / "residue_scorer.tsv", sep="\t", index=False)
    payload = {
        "contact_depletion": dep_contacts.to_dict(orient="records"),
        "gene_split": df[df["split"] == "gene"].to_dict(orient="records"),
        "family_split": df[df["split"] == "family"].to_dict(orient="records"),
    }
    (OUT / "residue_scorer.json").write_text(json.dumps(payload, indent=2) + "\n")

    cols = [
        "model",
        "t1_w0",
        "t1_w3",
        "t1_w5",
        "t1_w10",
        "t3_w3",
        "t3_w5",
        "t3_w10",
        "l1w_t3_w10",
        "p_top1_w5",
        "exp_lit_lig",
        "exp_seam_8A",
        "exp_crowded_q75",
    ]
    for sc in ("gene", "family"):
        print(f"\n=== {sc.upper()} SPLIT ===")
        print(df[df["split"] == sc][cols].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
