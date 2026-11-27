#!/usr/bin/env python3
"""Reproduce steric-clash (<2.5 A), +/-8 own_ref lDDT, and cassette_contact_res validation AUROCs.

Reads the ColabFold validation PDBs from ``refold_check_validation.py`` (``refold-validation-2026-11-20``)
and ``refold_check_validation_metrics_own_ref.tsv``, computes per-site heavy-atom cassette-native
steric clashes (<2.5 A, excluding the covalent junction residues +/-1) and +/-8 site lDDT against
each protein's own untagged ColabFold fold, and writes:

* ``refold_check_validation_clashes.tsv`` (per-site metrics)
* ``refold_check_validation_clashes_auroc.tsv`` (AUROC, 95% protein-bootstrap CI, Mann-Whitney p)

With ``--report-only`` (offline mode), skips PDB parsing and recomputes summary AUROCs and clash
enrichment rates directly from the committed ``refold_check_validation_clashes.tsv``.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE / "tag_refold"))

import refold_check_validation as rcv  # noqa: E402
import structure_metrics as sm  # noqa: E402

from deliverome_analysis.surfaceome import membrane_side as ms  # noqa: E402

CLASH_CUTOFF_A = 2.5
SITE_WINDOW_SHORT = 8
CASSETTE_LEN = 51
OUT_SCORES = HERE / "refold_check_validation_clashes.tsv"
OUT_AUROC = HERE / "refold_check_validation_clashes_auroc.tsv"
FOLDS_DIR = rcv.FOLDS


def _heavy_atoms_by_residue(pdb_path: Path) -> list[np.ndarray]:
    """Return list of (N_heavy, 3) float arrays indexed by 0-based residue order."""
    by_res: list[list[list[float]]] = []
    seen_res: dict[tuple[str, int, str], int] = {}
    for line in pdb_path.read_text().splitlines():
        if not line.startswith("ATOM  "):
            continue
        atom_name = line[12:16].strip()
        elem = line[76:78].strip() if len(line) >= 78 else ""
        if elem == "H" or (not elem and atom_name.startswith("H")):
            continue
        res_key = (line[21], int(line[22:26]), line[26])
        idx = seen_res.get(res_key)
        if idx is None:
            idx = len(by_res)
            seen_res[res_key] = idx
            by_res.append([])
        by_res[idx].append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return [np.asarray(coords, dtype=float) for coords in by_res]


def count_steric_clashes(
    pdb_path: Path, insert_after: int, cutoff: float = CLASH_CUTOFF_A
) -> tuple[int, int]:
    """Return (clash_residue_pairs, clashing_cassette_res) with any heavy-atom pair < cutoff."""
    res_atoms = _heavy_atoms_by_residue(pdb_path)
    c0 = insert_after
    c1 = insert_after + CASSETTE_LEN
    cas_indices = range(c0, min(c1, len(res_atoms)))
    nat_indices = [
        i for i in range(len(res_atoms)) if (i < c0 - 1 or i > c1) and len(res_atoms[i]) > 0
    ]
    if not nat_indices:
        return 0, 0
    nat_coords = np.vstack([res_atoms[i] for i in nat_indices])
    nat_owner = np.concatenate([np.full(len(res_atoms[i]), i, dtype=int) for i in nat_indices])
    cutoff2 = cutoff * cutoff
    total_res_pairs = 0
    clashing_cas_res = 0
    for ci in cas_indices:
        ca = res_atoms[ci]
        if len(ca) == 0:
            continue
        diff = ca[:, None, :] - nat_coords[None, :, :]
        d2 = np.einsum("ijk,ijk->ij", diff, diff)
        hit_atoms = np.any(d2 < cutoff2, axis=0)
        hit_nat_res = np.unique(nat_owner[hit_atoms])
        if len(hit_nat_res):
            total_res_pairs += len(hit_nat_res)
            clashing_cas_res += 1
    return total_res_pairs, clashing_cas_res


def lddt_site_window(
    tagged_pdb: Path, native_pdb: Path, insert_after: int, window: int = SITE_WINDOW_SHORT
) -> float:
    """C-alpha lDDT of native residues within +/-window of insert_after (scored vs native_pdb)."""
    t = ms.read_model(tagged_pdb)
    n = ms.read_model(native_pdb)
    c0 = insert_after
    c1 = insert_after + CASSETTE_LEN
    nat_idx = np.arange(len(n.ca))
    tag_idx = np.where(nat_idx < c0, nat_idx, nat_idx + (c1 - c0))
    if len(tag_idx) != len(n.ca) or int(tag_idx[-1]) >= len(t.ca):
        return float("nan")
    per_res = sm.lddt_ca(n.ca, t.ca[tag_idx])
    lo = max(0, insert_after - window)
    hi = min(len(n.ca), insert_after + window)
    mask = np.zeros(len(n.ca), dtype=bool)
    mask[lo:hi] = True
    mask &= n.plddt >= 70.0
    vals = per_res[mask]
    vals = vals[np.isfinite(vals)]
    return float(np.mean(vals)) if len(vals) else float("nan")


def compute_scores(folds_dir: Path = FOLDS_DIR) -> pd.DataFrame:
    df = rcv.load("own_ref")
    rows = []
    for r in df.to_dict("records"):
        tag_pdb = folds_dir / f"{r['fold_id']}.pdb"
        nat_pdb = folds_dir / f"n.{r['native_sha256'][:16]}.pdb"
        if not tag_pdb.is_file() or not nat_pdb.is_file():
            continue
        pos = int(r["insert_after_native_residue"])
        clashes, clash_res = count_steric_clashes(tag_pdb, pos)
        lddt8 = lddt_site_window(tag_pdb, nat_pdb, pos, window=SITE_WINDOW_SHORT)
        rows.append(
            {
                "fold_id": r["fold_id"],
                "uniprot_acc": r["uniprot_acc"],
                "gene_symbol": r["gene_symbol"],
                "class": r["class"],
                "insert_after_native_residue": pos,
                "native_aa": int(r["native_aa"]),
                "lddt_confident_own_ref": round(float(r["lddt_confident"]), 4),
                "lddt_site_8_own_ref": round(lddt8, 4) if np.isfinite(lddt8) else np.nan,
                "cassette_contact_res": float(r["cassette_contact_res"]),
                "clash_pairs_2p5A": clashes,
                "clash_cassette_res_2p5A": clash_res,
            }
        )
    out = pd.DataFrame(rows)
    out.to_csv(OUT_SCORES, sep="\t", index=False)
    return out


def summarize(df: pd.DataFrame) -> pd.DataFrame:
    work = df.copy()
    work["neg_cassette_contact_res"] = -work["cassette_contact_res"]
    work["neg_clash_pairs_2p5A"] = -work["clash_pairs_2p5A"]
    metrics = [
        ("lddt_confident_own_ref", "native-region lDDT (own_ref)"),
        ("lddt_site_8_own_ref", "site lDDT (+/-8, own_ref)"),
        ("neg_cassette_contact_res", "cassette_contact_res (negated)"),
        ("neg_clash_pairs_2p5A", "steric clashes <2.5 A (negated)"),
    ]
    out_rows = []
    for neg in ("known_failure", "decoy_core", "decoy_hardneg"):
        sub = work[work["class"].isin(["positive", neg])]
        for col, label in metrics:
            pos_v = sub.loc[sub["class"] == "positive", col].dropna().to_numpy(float)
            neg_v = sub.loc[sub["class"] == neg, col].dropna().to_numpy(float)
            auc = rcv.auroc(pos_v, neg_v)
            lo, hi = rcv.boot_auroc(sub, col, neg) if np.isfinite(auc) else (np.nan, np.nan)
            pval = float(mannwhitneyu(pos_v, neg_v, alternative="two-sided").pvalue)
            out_rows.append(
                {
                    "negative_class": neg,
                    "score": label,
                    "n_positive": len(pos_v),
                    "n_negative": len(neg_v),
                    "auroc": round(auc, 3),
                    "boot95_lo": round(lo, 3),
                    "boot95_hi": round(hi, 3),
                    "mwu_p": round(pval, 5),
                }
            )
    aucs = pd.DataFrame(out_rows)
    aucs.to_csv(OUT_AUROC, sep="\t", index=False)
    print(aucs.to_string(index=False))
    for thresh in (5, 8):
        rates = {}
        for cls in ("positive", "known_failure", "decoy_core"):
            s = work[work["class"] == cls]
            k = int((s["clash_pairs_2p5A"] >= thresh).sum())
            rates[cls] = (k, len(s), k / len(s) if len(s) else 0.0)
        enr = (
            rates["known_failure"][2] / rates["positive"][2]
            if rates["positive"][2]
            else float("nan")
        )
        print(
            f"clashes >= {thresh}: pos={rates['positive'][0]}/{rates['positive'][1]} "
            f"({rates['positive'][2]:.1%}), fail={rates['known_failure'][0]}/{rates['known_failure'][1]} "
            f"({rates['known_failure'][2]:.1%}), core={rates['decoy_core'][0]}/{rates['decoy_core'][1]} "
            f"({rates['decoy_core'][2]:.1%}), fail/pos enrichment={enr:.1f}x"
        )
    return aucs


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--report-only", action="store_true", help="summarize committed scores TSV offline"
    )
    ap.add_argument("--folds-dir", type=Path, default=FOLDS_DIR)
    args = ap.parse_args(argv)
    df = pd.read_csv(OUT_SCORES, sep="\t") if args.report_only else compute_scores(args.folds_dir)
    summarize(df)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
