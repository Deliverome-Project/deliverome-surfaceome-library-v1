#!/usr/bin/env python3
"""Per-region and per-segment tables for every L1 internal member (input to the agent review).

Exports what the internal decisions already used; no new modelling. Keyed on design_id.

`l1_internal_regions.tsv`, one row per extracellular region of each member:
    region bounds and length, eligible junctions, whether it passes the loop floor
    (V5.MIN_REGION_LEN), v5 and v6 region ceilings (best absolute score) and best rank
    residues, the v5 ceiling's percentile among published regions, the member's
    order-rule verdict, and which arms sit in the region.
    Regions are mask r8's (membrane-side gate) where the member has a usable r8 row.
    Otherwise (`region_basis`) they are the extracellular runs of the member's feature table:
    the rescue round's span rules for single-pass and GPI members, UniProt topology for the
    two ENTPD2 members.

`l1_internal_segments.tsv`, one row per DeepTMHMM v1 segment (topology of record) of each
member's ordered protein: the segment's v1 state (O outside, I inside, M helix, S signal
peptide, B beta), bounds, and how many of its residues and junctions the r8 mask makes
eligible, with the arms inside it. A member whose ordered protein has no v1 row has none.

    uv run --locked python l1_region_tables.py
"""

from __future__ import annotations

import json
import sys
from itertools import groupby
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[3] / "src"))

import l1_internal_decisions as L  # noqa: E402
import membrane_mask  # noqa: E402
import rescue_round as RR  # noqa: E402
import tag_site_model_v5 as V5  # noqa: E402
import tag_site_model_v6 as V6  # noqa: E402
from redesign_all import frames  # noqa: E402
from rescue_round import score  # noqa: E402
from v5_region_table import SEPARATION, features, ranges  # noqa: E402

from deliverome_analysis import topology  # noqa: E402

ARMS = HERE / "l1_internal_arm_decisions.tsv"
MEMBERS = HERE / "l1_internal_arm_members.tsv"
MASK_R8 = HERE / "data" / "mask_r8" / "l1_extracellular_mask_residues.tsv"
OUT_REGIONS = HERE / "l1_internal_regions.tsv"
OUT_SEGMENTS = HERE / "l1_internal_segments.tsv"
OUT_JSON = HERE / "l1_region_tables.json"


def arm_label(x) -> str:
    site = "snorkel" if x.arm_kind == "snorkel" else f"after{int(x.insert_after_residue)}"
    return f"{x.arm_id}@{site}:{x.decision}"


def main() -> int:
    arms = pd.read_csv(ARMS, sep="\t")
    mem = pd.read_csv(MEMBERS, sep="\t")
    dec = pd.read_csv(L.OUT, sep="\t")
    _i, _d, prot = L.head_tables()
    sha_of = arms.drop_duplicates("design_id").set_index("design_id")["sequence_sha256"]
    mem = mem[mem["design_id"].isin(sha_of.index)].copy()
    mem["sequence_sha256"] = mem["design_id"].map(sha_of)

    # feature rows for every member: the frames the designs read, rescue members on theirs
    fr = frames(dec)
    want = {fr[s]: s for s in mem["sequence_sha256"] if s in fr}
    allx = features()
    pv = pd.read_csv(RR.PROVENANCE, sep="\t")
    rk = dict(zip(pv["key"], pv["sequence_sha256"], strict=True))
    rf = pd.read_csv(RR.FEATURES, sep="\t", low_memory=False)
    want = {k: s for k, s in want.items() if s not in set(rk.values())} | rk
    allx = pd.concat([allx[~allx["uniprot_acc"].isin(rk)], rf[rf["uniprot_acc"].isin(rk)]])
    sub = allx[allx["uniprot_acc"].isin(want)]
    s5, s6 = score(sub, V5), score(sub, V6)
    s5["sha"] = s5["uniprot_acc"].map(want)
    s5["v6_rank"] = s6["rank_score"].to_numpy()
    s5["v6_abs"] = s6["absolute_score"].to_numpy()
    spec = json.loads(V5.SPEC.read_text())
    pub = np.sort(np.array(spec["published_region_ceilings_oof"]))
    r8 = pd.read_csv(MASK_R8, sep="\t").drop_duplicates("sequence_sha256")
    r8 = r8.set_index("sequence_sha256")
    rescue_shas = set(rk.values())

    built = arms[arms["decision"] == "build"]
    rows, seg_rows = [], []
    by_sha = dict(tuple(s5.groupby("sha")))
    for m in mem.itertuples(index=False):
        sha = m.sequence_sha256
        g0 = by_sha.get(sha)
        my = arms[arms["design_id"] == m.design_id]
        loops = my[my["arm_kind"] == "loop"]
        member = {
            "design_id": m.design_id, "member_id": m.member_id, "gene": m.gene,
            "sequence_sha256": sha, "member_decision": m.decision,
            "member_status": m.status, "order_rule_sites": m.recommended_sites_to_order,
            "order_rule_why": m.why, "all_bad": m.all_bad,
            "arms_built": int(m.arms_built),
        }  # fmt: skip
        has_r8 = sha in r8.index and int(r8.loc[sha, "usable"]) == 1
        if g0 is not None:
            if has_r8:
                mk = (r8.loc[sha, "mask"], r8.loc[sha, "junction_mask"])
                g = g0[membrane_mask.eligible_after(mk, g0["residue"])].copy()
                g["region"] = membrane_mask.masked_regions(g["residue"], mk)
                basis = "mask r8"
            else:
                g = g0.assign(region=g0["region_index"])
                basis = (
                    "rescue-round span rules (no r8 row)" if sha in rescue_shas
                    else "feature table's extracellular runs (no usable r8 row)"
                )  # fmt: skip
            for ri, rg in g.groupby("region"):
                res = set(rg["residue"].astype(int))
                rg = rg.sort_values("rank_score", ascending=False)
                top = []
                for x in rg.itertuples(index=False):
                    if all(abs(x.residue - t.residue) >= SEPARATION for t in top):
                        top.append(x)
                    if len(top) == 3:
                        break
                p5 = rg.loc[rg["absolute_score"].idxmax()]
                p6 = rg.loc[rg["v6_abs"].idxmax()]
                inreg = loops[loops["insert_after_residue"].isin(res)]
                rows.append({
                    **member, "region_basis": basis, "region_index": int(ri),
                    "start": int(rg["residue"].min()), "end": int(rg["residue"].max()),
                    "n_eligible_junctions": len(rg),
                    "passes_loop_floor": len(rg) >= V5.MIN_REGION_LEN,
                    "eligible_junctions": ranges(rg["residue"]),
                    "v5_region_ceiling": round(float(p5["absolute_score"]), 4),
                    "v5_ceiling_residue": int(p5["residue"]),
                    "v5_ceiling_published_pct": round(
                        float(np.searchsorted(pub, p5["absolute_score"]) / len(pub)) * 100, 1),
                    "v5_top_rank_residues": ";".join(
                        f"{int(t.residue)}:{t.rank_score:.3f}" for t in top),
                    "v6_region_ceiling": round(float(p6["v6_abs"]), 4),
                    "v6_ceiling_residue": int(p6["residue"]),
                    "v6_best_rank_residue": int(rg.loc[rg["v6_rank"].idxmax(), "residue"]),
                    "arms_in_region": ";".join(arm_label(x) for x in inreg.itertuples()),
                    "built_arms_in_region": int((inreg["decision"] == "build").sum()),
                })  # fmt: skip
        else:
            rows.append({**member, "region_basis": "no feature rows", "region_index": None})

    # segments on DeepTMHMM v1, the topology of record
    seqs = {m.design_id: prot[m.sequence_sha256] for m in mem.itertuples(index=False)}
    topo = topology.resolve(seqs, tool=topology.TOOL_V1).set_index("key")
    for m in mem.itertuples(index=False):
        t = topo.loc[m.design_id]
        my = built[built["design_id"] == m.design_id]
        sha = m.sequence_sha256
        has_r8 = sha in r8.index and int(r8.loc[sha, "usable"]) == 1
        mask = r8.loc[sha, "mask"] if has_r8 else None
        junc = r8.loc[sha, "junction_mask"] if has_r8 else None
        if not bool(t["found"]):
            seg_rows.append({"design_id": m.design_id, "gene": m.gene, "sequence_sha256": sha,
                             "topology_tool": topology.TOOL_V1, "state": None,
                             "note": "no DeepTMHMM v1 row for the ordered protein"})  # fmt: skip
            continue
        s = str(t["per_residue_topology"])
        pos = 1
        for k, (state, run) in enumerate(groupby(s)):
            n = len(list(run))
            a, b = pos, pos + n - 1
            pos += n
            inside = my[my["insert_after_residue"].between(a, b)]
            seg_rows.append({
                "design_id": m.design_id, "gene": m.gene, "sequence_sha256": sha,
                "topology_tool": topology.TOOL_V1,
                "topology_version": t.get("topology_version"),
                "segment_index": k, "state": state, "start": a, "end": b, "length": n,
                "mask_r8_eligible_residues": None if mask is None else mask[a - 1:b].count("E"),
                "mask_r8_eligible_junctions": None if junc is None else junc[a - 1:b].count("E"),
                "mask_r8_string": None if mask is None else mask[a - 1:b],
                "built_arms_in_segment": ";".join(arm_label(x) for x in inside.itertuples()),
            })  # fmt: skip

    reg = pd.DataFrame(rows)
    seg = pd.DataFrame(seg_rows)
    reg.to_csv(OUT_REGIONS, sep="\t", index=False)
    seg.to_csv(OUT_SEGMENTS, sep="\t", index=False)
    summary = {
        "members": int(mem["design_id"].nunique()),
        "members_with_regions": int(reg.dropna(subset=["region_index"])["design_id"].nunique()),
        "regions": int(reg["region_index"].notna().sum()),
        "region_basis": reg.drop_duplicates("design_id")["region_basis"].value_counts().to_dict(),
        "members_with_v1_segments": int(seg.dropna(subset=["state"])["design_id"].nunique()),
        "segments": int(seg["state"].notna().sum()),
        "built_loop_arms_placed_in_a_region": int(reg["built_arms_in_region"].fillna(0).sum()),
        "built_loop_arms": int((built["arm_kind"] == "loop").sum()),
    }
    OUT_JSON.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
