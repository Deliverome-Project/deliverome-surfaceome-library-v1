#!/usr/bin/env python3
"""How many internal arms each L1 member would get under the multi-site rule (a proposal, not a build).

Becca, 2026-11-15 00:25 (relayed by the coordinator): "where tags don't look good we had planned
2 or 3 tag sites based on delta in score between best and worst and versus library median as
well as if they are all bad perhaps a tag and snorkel. never use the snorkel alone."

The rule already exists in v5's design step, inherited from v4 (2026-10-15, README "How many
sites to order"; `recommended_sites_to_order` in every tag_sites_v5_*.tsv), and the L1 build
used only rank 1. This counts, on the current internal decisions, how many of each member's
recommended sites survive what every arm must pass:

    contact veto   no ligand contact within the +/-2 window (contact_veto.vetoed)
    gate           the junction is E in the membrane-side mask r8 (members without an r8 row
                   are reported as gate-unknown)
    refold         the cached verdict for exactly that tagged protein (refold_verdicts_seen.tsv);
                   unfolded sites are reported as refold-pending, not as passing

and lists, for the snorkel-only and held members, their best loop arm and snorkel verdicts.

    uv run --locked python site_count_proposal.py    # writes site_count_proposal*.tsv/.json
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[3] / "src"))

import contact_veto as CV  # noqa: E402
import tag_site_model_v5 as V5  # noqa: E402
from build_internal_tagged import CASSETTE_AA  # noqa: E402
from l1_internal_decisions import head_tables  # noqa: E402
from snorkel_design import MODULE  # noqa: E402

DESIGNS = [  # later tables replace earlier ones for the same protein, as in the decisions
    "tag_sites_v5_adopted.tsv",
    "tag_sites_v5_orthologs.tsv",
    "data/tag_sites_v5_start_trims.tsv",
    "tag_sites_v5_single_pass.tsv",
    "tag_sites_v5_orientation_unconfirmed.tsv",
    "tag_sites_v5_routed_in.tsv",
    "tag_sites_v5_isoform_fix.tsv",
]
MASK = HERE / "data" / "mask_r8" / "l1_extracellular_mask_residues.tsv"
OUT = HERE / "site_count_proposal.tsv"
OUT_ARMS = HERE / "site_count_proposal_arms.tsv"
OUT_SNORKEL = HERE / "site_count_proposal_snorkel_and_held.tsv"
OUT_JSON = HERE / "site_count_proposal.json"


def designs() -> pd.DataFrame:
    t = pd.read_csv(HERE / "l1_adopted_targets.tsv", sep="\t")
    a2s = dict(zip(t["uniprot_acc"], t["sequence_sha256"], strict=False))
    tabs = []
    for i, f in enumerate(DESIGNS):
        x = pd.read_csv(HERE / f, sep="\t")
        if "sequence_sha256" not in x:
            x["sequence_sha256"] = x["uniprot_acc"].map(a2s)
        tabs.append(x.assign(_src=f, _p=i))
    d = pd.concat(tabs, ignore_index=True).dropna(subset=["sequence_sha256"])
    keep = d.groupby("sequence_sha256")["_p"].transform("max")
    return d[d["_p"] == keep]


def main() -> int:
    dec = pd.read_csv(HERE / "l1_internal_decisions.tsv", sep="\t")
    d = designs()
    _i, _ids, prot = head_tables()
    mask = pd.read_csv(MASK, sep="\t").set_index("sequence_sha256")
    seen = pd.read_csv(HERE / "refold_verdicts_seen.tsv", sep="\t").set_index("tagged_sha256")
    contacts = CV.member_contacts(dec[["design_id", "sequence_sha256"]], prot, download=False)
    spec = json.loads(V5.SPEC.read_text())
    pub = np.array(spec["published_region_ceilings_oof"])
    q25, med = (float(np.quantile(pub, q)) for q in (0.25, 0.5))

    arms, members = [], []
    for r in dec.itertuples(index=False):
        g = d[d["sequence_sha256"] == r.sequence_sha256].sort_values("rank")
        if g.empty or r.sequence_sha256 not in prot:
            members.append({"design_id": r.design_id, "gene": r.gene, "decision": r.decision,
                            "note": "no v5 design row"})  # fmt: skip
            continue
        seq = prot[r.sequence_sha256]
        top = g.iloc[0]
        n_rec = int(top["recommended_sites_to_order"])
        m = mask.loc[r.sequence_sha256] if r.sequence_sha256 in mask.index else None
        ok = []
        for x in g.itertuples(index=False):
            n = int(x.insert_after_residue)
            veto = bool(CV.vetoed(contacts, r.sequence_sha256, n))
            if m is None or int(m["usable"]) != 1:
                gate = "unknown"
            else:
                gate = "pass" if m["junction_mask"][n - 1] == "E" else "fail"
            tsha = hashlib.sha256((seq[:n] + CASSETTE_AA + seq[n:]).encode()).hexdigest()
            refold = seen.loc[tsha, "verdict"] if tsha in seen.index else "not_folded"
            arms.append({"design_id": r.design_id, "gene": r.gene, "decision": r.decision,
                         "rank": int(x.rank), "insert_after_residue": n,
                         "region_index": x.region_index,
                         "rank_score": x.rank_score, "absolute_score": x.absolute_score,
                         "region_ceiling": x.region_ceiling, "contact_veto": veto,
                         "gate_r8": gate, "refold": refold, "tagged_sha256": tsha})  # fmt: skip
            if not veto and gate != "fail" and refold != "fail":
                ok.append((n, refold))
        ceiling, runner = float(top["region_ceiling"]), top["runner_up_region_ceiling"]
        members.append(
            {
                "design_id": r.design_id,
                "gene": r.gene,
                "decision": r.decision,
                "best_region_ceiling": round(ceiling, 4),
                "runner_up_region_ceiling": runner,
                "best_minus_runner_up": None if pd.isna(runner) else round(ceiling - runner, 4),
                "recommended_sites_to_order": n_rec,
                "why": top["why"],
                "all_bad": ceiling < q25,
                "designed_sites": len(g),
                "sites_clear_of_veto_gate_refold": len(ok),
                "sites_refold_pass": sum(v == "pass" for _n, v in ok),
            }
        )
    a, mm = pd.DataFrame(arms), pd.DataFrame(members)
    a.to_csv(OUT_ARMS, sep="\t", index=False)
    mm.to_csv(OUT, sep="\t", index=False)

    # snorkel-only and held members: best loop arm and the snorkel arm
    rows = []
    for r in dec[dec["decision"].isin(["snorkel", "hold_no_passing_arm"])].itertuples():
        x = a[(a["design_id"] == r.design_id)]
        gp = x[(x["gate_r8"] != "fail") & ~x["contact_veto"]].sort_values("rank")
        best = gp.iloc[0] if len(gp) else (x.sort_values("rank").iloc[0] if len(x) else None)
        snork = hashlib.sha256((prot.get(r.sequence_sha256, "") + MODULE).encode()).hexdigest()
        rows.append(
            {
                "design_id": r.design_id,
                "gene": r.gene,
                "decision": r.decision,
                "snorkel_verdict": r.snorkel_verdict,
                "snorkel_tagged_sha256": snork,
                "best_loop_site": None if best is None else int(best["insert_after_residue"]),
                "best_loop_rank": None if best is None else int(best["rank"]),
                "best_loop_gate_r8": None if best is None else best["gate_r8"],
                "best_loop_contact_veto": None if best is None else bool(best["contact_veto"]),
                "best_loop_refold": None if best is None else best["refold"],
                "n_designed_loop_sites": len(x),
                "n_loop_sites_gate_ok_and_clean": len(gp),
                "n_loop_sites_refold_pass": int((x["refold"] == "pass").sum()),
            }
        )
    pd.DataFrame(rows).to_csv(OUT_SNORKEL, sep="\t", index=False)
    OUT_JSON.write_text(
        json.dumps(
            {
                "published_region_ceiling_q25": round(q25, 4),
                "published_region_ceiling_median": round(med, 4),
                "library_median_best_ceiling": round(float(mm["best_region_ceiling"].median()), 4),
                "members": len(mm),
            },
            indent=2,
        )
        + "\n"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
