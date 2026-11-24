#!/usr/bin/env python3
"""Per-region v5 scores for chosen L1 internal members (input to the agent review).

v5's design step keeps only its picks. This recomputes, for each member, every residue's v5
scores exactly as tag_site_model_v4.design does (both forests on the member's features,
z-scored within the protein, smoothed over +/-2), applies the membrane-side mask r8 (the gate
the internal decisions use) and summarises each eligible region:

    region_index (r8), start, end, length, passes the loop floor (V5.MIN_REGION_LEN)
    region_ceiling (best absolute score) and its residue, with the published-region percentile
    the top three junctions by rank_score, at least 5 residues apart
    the eligible junction ranges
    the v5 design picks that fall in the region (rank and site), so the agent can join on them

The design tables' own `region_index` was computed under the mask each design ran with (r4 for
the adopted panel); this table's regions are r8's, and picks are matched by residue. A member
with no usable r8 row keeps DeepTMHMM v1's regions, and its `mask` column says why.

    uv run --locked python v5_region_table.py --pilot   # the agent pilot's members
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[3] / "src"))

import membrane_mask  # noqa: E402
import tag_site_model_v5 as V5  # noqa: E402
from rescue_round import score  # noqa: E402
from site_count_proposal import designs  # noqa: E402

OUT = HERE / "v5_regions_pilot.tsv"
MASK_R8 = HERE / "data" / "mask_r8" / "l1_extracellular_mask_residues.tsv"
#: every feature table a v5 design read; later ones replace a protein's rows (override)
BASE = ["all_extracellular_targets.tsv", "l1_member_extracellular.tsv", "slc67a1_extracellular.tsv"]
EXTRA = [
    "l1_ortholog_member_extracellular.tsv",
    "l1_start_trim_member_extracellular.tsv",
    "l1_routed_in_member_extracellular.tsv",
    "entpd2_extracellular.tsv",
]
OVERRIDE = ["data/allx_isoform_fix15_2026-11-14.tsv"]
PROVENANCE = [
    "l1_ortholog_member_features.tsv",
    "l1_start_trim_member_features.tsv",
    "l1_routed_in_member_features.tsv",
]
SEPARATION = 5


def features() -> pd.DataFrame:
    parts = []
    for f in BASE + EXTRA:
        p = HERE / f
        if p.is_file():
            parts.append(pd.read_csv(p, sep="\t", low_memory=False))
    allx = pd.concat(parts, ignore_index=True).drop_duplicates(["uniprot_acc", "residue"])
    for f in OVERRIDE:
        o = pd.read_csv(HERE / f, sep="\t", low_memory=False)
        allx = pd.concat([allx[~allx["uniprot_acc"].isin(set(o["uniprot_acc"]))], o])
    return allx.reset_index(drop=True)


def ranges(res) -> str:
    res = sorted(res)
    out, s, p = [], None, None
    for r in res:
        if s is None:
            s = p = r
        elif r == p + 1:
            p = r
        else:
            out.append(f"{s}-{p}")
            s = p = r
    if s is not None:
        out.append(f"{s}-{p}")
    return ";".join(out)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pilot", action="store_true", help="all_bad, snorkel and held members")
    ap.add_argument("--genes", help="comma-separated gene symbols (internal members)")
    ap.add_argument("--out", default=str(OUT))
    args = ap.parse_args(argv)
    m = pd.read_csv(HERE / "site_count_proposal.tsv", sep="\t")
    if args.genes:
        m = m[m["gene"].isin(args.genes.split(","))]
    if args.pilot:
        m = m[(m["all_bad"] == True) | m["decision"].isin(["snorkel", "hold_no_passing_arm"])]  # noqa: E712
    dec = pd.read_csv(HERE / "l1_internal_decisions.tsv", sep="\t").set_index("design_id")
    d = designs()
    frame = (
        d.dropna(subset=["uniprot_acc"])
        .drop_duplicates("sequence_sha256")
        .set_index("sequence_sha256")["uniprot_acc"]
        .to_dict()
    )
    # orthologs, start-trimmed and routed-in members are keyed on member_id@commit
    for f in PROVENANCE:
        if (HERE / f).is_file():
            pv = pd.read_csv(HERE / f, sep="\t")
            for did, key in zip(pv["design_id"], pv["key"], strict=False):
                if did in dec.index:
                    frame.setdefault(dec.loc[did, "sequence_sha256"], key)
    spec = json.loads(V5.SPEC.read_text())
    pub = np.sort(np.array(spec["published_region_ceilings_oof"]))
    want = {}
    for did in m["design_id"]:
        sha = dec.loc[did, "sequence_sha256"]
        if sha in frame:
            want[frame[sha]] = (did, sha)
    allx = features()
    s = score(allx[allx["uniprot_acc"].isin(want)])
    r8 = pd.read_csv(MASK_R8, sep="\t").drop_duplicates("sequence_sha256")
    r8 = r8.set_index("sequence_sha256")
    rows = []
    for acc, g0 in s.groupby("uniprot_acc"):
        did, sha = want[acc]
        if sha in r8.index and int(r8.loc[sha, "usable"]) == 1:
            mk = (r8.loc[sha, "mask"], r8.loc[sha, "junction_mask"])
            g = g0[membrane_mask.eligible_after(mk, g0["residue"])].copy()
            g["region_index"] = membrane_mask.masked_regions(g["residue"], mk)
            mask_note = "r8"
        else:
            g = g0.copy()  # regions stay DeepTMHMM v1's
            mask_note = (
                "not applied: no r8 row"
                if sha not in r8.index
                else f"not applied: r8 unusable ({r8.loc[sha, 'status']})"
            )
        picks = d[d["sequence_sha256"] == sha][["rank", "insert_after_residue"]]
        for ri, rg in g.groupby("region_index"):
            rg = rg.sort_values("rank_score", ascending=False)
            top = []
            for x in rg.itertuples(index=False):
                if all(abs(x.residue - t.residue) >= SEPARATION for t in top):
                    top.append(x)
                if len(top) == 3:
                    break
            peak = rg.loc[rg["absolute_score"].idxmax()]
            inreg = picks[picks["insert_after_residue"].isin(set(rg["residue"]))]
            row = {
                "design_id": did,
                "frame_accession": acc,
                "sequence_sha256": sha,
                "region_index": int(ri),
                "mask": mask_note,
                "start": int(rg["residue"].min()),
                "end": int(rg["residue"].max()),
                "length": len(rg),
                "passes_loop_floor": len(rg) >= V5.MIN_REGION_LEN,
                "region_ceiling": round(float(peak["absolute_score"]), 4),
                "region_absolute_peak_residue": int(peak["residue"]),
                "region_ceiling_published_pct": round(
                    float(np.searchsorted(pub, peak["absolute_score"]) / len(pub)) * 100, 1
                ),
                "eligible_junctions": ranges(rg["residue"]),
                "v5_picks_in_region": ";".join(
                    f"rank{int(p.rank)}@{int(p.insert_after_residue)}"
                    for p in inreg.itertuples(index=False)
                ),
            }
            for k, t in enumerate(top, 1):
                row[f"top{k}_residue"] = int(t.residue)
                row[f"top{k}_rank_score"] = round(float(t.rank_score), 4)
                row[f"top{k}_absolute_score"] = round(float(t.absolute_score), 4)
            rows.append(row)
    out = pd.DataFrame(rows).sort_values(["design_id", "start"])
    out.to_csv(args.out, sep="\t", index=False)
    missing = sorted(set(m["design_id"]) - set(out["design_id"]))
    print(
        f"{out['design_id'].nunique()} members, {len(out)} regions -> {args.out}; "
        f"mask: {out.drop_duplicates('design_id')['mask'].value_counts().to_dict()}; "
        f"no extracellular features for {len(missing)}: {missing}"
    )
    print(f"published-region q25 {np.quantile(pub, 0.25):.4f} (the all_bad cutoff)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
