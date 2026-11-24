#!/usr/bin/env python3
"""Revise the built v5 arms under the 2026-11-15 contact rulings, with v6 as the scorer of record.

The coordinator, under Becca's delegation (2026-11-15): the library is built on v5's picks,
which all have refold verdicts; v6 does not measurably improve on them, and site choice within a
region is weakly determined (198 of 575 rank-1 sites move between two statistically
indistinguishable models). v6 scores every site the rulings below force, and adds an arm where it
disagrees with v5 about the best region.

    hard veto       structure or transferred ligand contacts, and protein_ligand partners (SHH
                    and the like), at residues n-2 .. n+3 of an insert after n: the arm is
                    replaced by v6's best clean junction in the same region, else in another
                    region apart from the member's other arms; with none, it is dropped
    soft veto       structure protein-partner contacts, literature homomer contacts and
                    literature ligand contacts without structure evidence: replaced when a clean
                    junction in the same region scores at least 90% of the arm's own v6 rank score;
                    otherwise dropped if the junction itself is a contact residue (distance 0) and
                    the member keeps another clean arm; otherwise kept, flagged with the contact
    flag only       literature partner rows
    DNA             a site whose cassette seam creates a restriction site the builder cannot
                    recode (KCND1 l1-p0999 after 286, Esp3I) is replaced as under the hard veto
    v6 region       where v6's best region differs from v5's (both under r8), v6's best clean site
                    in that region is added when the member has fewer than three loop arms;
                    otherwise the disagreement is flagged on the member

Why partner and homomer contacts are soft: a veto on them would remove 10 (partner) or 17
(partner and homomer) of the 176 published, working internal tag sites; the ligand veto removes 2
(hard_negative_benchmark.json).

Writes arm_revisions.tsv (one row per action) and arm_revisions_new_constructs.tsv (the tagged
proteins that still need a refold). Nothing is built here: l1_internal_arms.py applies the
revisions once their verdicts exist.

    uv run --locked python arm_revisions.py
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

import membrane_mask  # noqa: E402
import tag_site_model_v6 as V6  # noqa: E402
from build_internal_tagged import CASSETTE_AA  # noqa: E402
from l1_internal_decisions import head_tables  # noqa: E402
from redesign_all import frames  # noqa: E402
from rescue_round import score  # noqa: E402
from v5_region_table import features  # noqa: E402

OUT = HERE / "arm_revisions.tsv"
OUT_NEW = HERE / "arm_revisions_new_constructs.tsv"
OUT_JSON = HERE / "arm_revisions.json"
MASK_R8 = HERE / "data" / "mask_r8" / "l1_extracellular_mask_residues.tsv"
LIG = {"literature_small_ligand", "literature_peptide_ligand", "literature_protein_ligand"}
SEPARATION = 25
DNA_CONFLICT = {("l1-p0999", 286): "cassette seam creates an Esp3I site (l1-integration, #360)"}
DECIDED_BY = "coordinator under Becca's 2026-11-12 delegation, 2026-11-15 contact and v6 rulings"


def contact_tiers(c: pd.DataFrame) -> pd.DataFrame:
    lit = c["class"].str.startswith("literature_")
    tier = np.select(
        [
            c["class"].isin(["ligand", "protein_ligand"]),
            c["class"].isin(LIG) & (c["evidence"] == "structure"),
            (c["class"] == "protein_partner") & ~lit,
            c["class"].isin(LIG) | (c["class"] == "literature_homomer"),
            c["class"] == "literature_partner",
        ],
        ["hard", "hard", "soft", "soft", "flag"],
        default="none",
    )
    return c.assign(tier=tier)[lambda d: d["tier"] != "none"]


def hits(ct: pd.DataFrame, sha: str, n: int) -> pd.DataFrame:
    x = ct[(ct["sequence_sha256"] == sha) & ct["residue"].between(n - 2, n + 3)]
    d = (x["residue"] - n).where(x["residue"] <= n, x["residue"] - (n + 1)).abs()
    return x.assign(distance=d)


def main() -> int:
    arms = pd.read_csv(HERE / "l1_internal_arm_decisions.tsv", sep="\t")
    dec = pd.read_csv(HERE / "l1_internal_decisions.tsv", sep="\t")
    loops = arms[(arms["decision"] == "build") & (arms["arm_kind"] == "loop")].copy()
    ct = contact_tiers(pd.read_csv(HERE / "contacts_by_protein.tsv", sep="\t", low_memory=False))
    _i, _d, prot = head_tables()
    seen = set(pd.read_csv(HERE / "refold_verdicts_seen.tsv", sep="\t")["tagged_sha256"])
    seen |= set(pd.read_csv(HERE / "data" / "verdicts_loop.tsv", sep="\t")["tagged_sha256"])
    # v6 scores on every residue of every internal member, regions under r8
    fr = frames(dec)
    allx = features()
    want = {fr[s]: s for s in loops["sequence_sha256"].unique() if s in fr}
    sc = score(allx[allx["uniprot_acc"].isin(want)], V6)
    r8 = pd.read_csv(MASK_R8, sep="\t").drop_duplicates("sequence_sha256")
    r8 = r8.set_index("sequence_sha256")
    cand = {}
    for acc, g in sc.groupby("uniprot_acc"):
        sha = want[acc]
        if sha in r8.index and int(r8.loc[sha, "usable"]) == 1:
            mk = (r8.loc[sha, "mask"], r8.loc[sha, "junction_mask"])
            g = g[membrane_mask.eligible_after(mk, g["residue"])].copy()
            g["region"] = membrane_mask.masked_regions(g["residue"], mk)
        else:
            g = g.assign(region=g["region_index"])
        cand[sha] = g[["residue", "region", "rank_score", "absolute_score"]].reset_index(drop=True)

    from annotation_rule import Checker

    checker = Checker()

    def clean(sha, n):
        return hits(ct[ct["tier"].isin(["hard", "soft"])], sha, n).empty

    def region_of(sha, n):
        g = cand.get(sha)
        if g is None:
            return None
        r = g.loc[g["residue"] == n, "region"]
        return int(r.iloc[0]) if len(r) else None

    def v6_score(sha, n):
        g = cand.get(sha)
        r = None if g is None else g.loc[g["residue"] == n, "rank_score"]
        return float(r.iloc[0]) if r is not None and len(r) else np.nan

    def best_clean(sha, region, others, same_region_only, exclude=()):
        g = cand.get(sha)
        if g is None:
            return None
        g = g.sort_values("rank_score", ascending=False)
        for x in g.itertuples(index=False):
            if same_region_only and x.region != region:
                continue
            n = int(x.residue)
            if n in exclude:
                continue
            if not clean(sha, n):
                continue
            if checker.blocked(sha, n):  # annotation rule (annotation_rule.py), 2026-11-15
                continue
            if any((region_of(sha, m) == x.region) and abs(m - n) < SEPARATION for m in others):
                continue
            return n, int(x.region), float(x.rank_score)
        return None

    rows = []
    for did, g in loops.groupby("design_id"):
        sha = g["sequence_sha256"].iloc[0]
        sites = [int(n) for n in g["insert_after_residue"]]
        keep_sites = list(sites)
        for r in g.itertuples(index=False):
            n = int(r.insert_after_residue)
            others = [m for m in keep_sites if m != n]
            h = hits(ct, sha, n)
            hard = h[h["tier"] == "hard"]
            soft = h[h["tier"] == "soft"]
            flag = h[h["tier"] == "flag"]
            reg = region_of(sha, n)
            base = {"design_id": did, "gene": r.gene, "sequence_sha256": sha,
                    "arm_id": r.arm_id, "old_site": n, "old_region": reg,
                    "old_v6_rank_score": round(v6_score(sha, n), 4)}  # fmt: skip

            def describe(x):
                return "; ".join(
                    f"{t.source}:{t.partner} at {t.residue} (distance {t.distance})"
                    for t in x.drop_duplicates(["residue", "source"]).itertuples(index=False)
                )

            if (did, n) in DNA_CONFLICT or len(hard):
                why = DNA_CONFLICT.get((did, n)) or f"hard contact: {describe(hard)}"
                ex = (n,)  # never the site being replaced
                rep = best_clean(sha, reg, others, True, ex) or best_clean(
                    sha, reg, others, False, ex
                )
                if rep:
                    keep_sites = [*others, rep[0]]
                    rows.append({**base, "action": "replace", "new_site": rep[0],
                                 "new_region": rep[1], "new_v6_rank_score": round(rep[2], 4),
                                 "reason": why})  # fmt: skip
                else:
                    keep_sites = others
                    rows.append({**base, "action": "drop", "reason": f"{why}; no clean site"})
                continue
            if len(soft):
                own = v6_score(sha, n)
                rep = best_clean(sha, reg, others, True)
                if rep and np.isfinite(own) and rep[2] >= 0.9 * own:
                    keep_sites = [*others, rep[0]]
                    rows.append({**base, "action": "replace", "new_site": rep[0],
                                 "new_region": rep[1], "new_v6_rank_score": round(rep[2], 4),
                                 "reason": f"soft contact: {describe(soft)}; a clean site in "
                                 "the same region scores within 10%"})  # fmt: skip
                    continue
                on_site = (soft["distance"] == 0).any()
                other_clean = [m for m in others if clean(sha, m)]
                if on_site and other_clean:
                    keep_sites = others
                    rows.append({**base, "action": "drop",
                                 "reason": f"soft contact at the junction itself: "
                                 f"{describe(soft)}; the member keeps another clean arm"})  # fmt: skip
                    continue
                rows.append({**base, "action": "keep_flagged",
                             "reason": f"soft contact: {describe(soft)}; no clean site within "
                             "10% in the region"})  # fmt: skip
                continue
            if len(flag):
                rows.append({**base, "action": "keep_flagged",
                             "reason": f"literature partner (flag only): {describe(flag)}"})  # fmt: skip
        # v6 region-level disagreement
        rows_member = [x for x in rows if x["design_id"] == did]
        n_loops = len(sites) - sum(x["action"] == "drop" for x in rows_member)
        rows.extend(region_disagreement(did, sha, r.gene, n_loops, keep_sites, best_clean,
                                        region_of))  # fmt: skip

    rev = pd.DataFrame(rows)
    new = []
    for x in rev[rev["action"].isin(["replace", "add"])].itertuples(index=False):
        n = int(x.new_site)
        seq = prot[x.sequence_sha256]
        t = hashlib.sha256((seq[:n] + CASSETTE_AA + seq[n:]).encode()).hexdigest()
        if t not in seen:
            new.append({"design_id": x.design_id, "gene": x.gene,
                        "sequence_sha256": x.sequence_sha256, "insert_after_residue": n,
                        "action": x.action, "tagged_sha256": t, "length_aa": len(seq) + 51,
                        "tagged_protein": seq[:n] + CASSETTE_AA + seq[n:]})  # fmt: skip
    nc = pd.DataFrame(new).drop_duplicates("tagged_sha256")
    rev = rev.assign(decided_by=DECIDED_BY, implemented_by="ranker-pu-segment session",
                     decided_on="2026-11-15")  # fmt: skip
    rev.to_csv(OUT, sep="\t", index=False)
    nc.to_csv(OUT_NEW, sep="\t", index=False)
    aa = int(nc["length_aa"].sum()) if len(nc) else 0
    summary = {
        "actions": rev["action"].value_counts().to_dict(),
        "members_with_v6_region_disagreement": int(
            rev.loc[rev["action"].isin(["add", "flag_member"]), "design_id"].nunique()
        ),
        "new_constructs_to_fold": len(nc),
        "residues_to_fold": aa,
        "cost_usd_range": [round(aa * 0.000062, 2), round(aa * 0.00008, 2)],
    }
    OUT_JSON.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return 0


def region_disagreement(did, sha, gene, n_loops, keep_sites, best_clean, region_of):
    """v6's best region vs v5's (both designed under r8 by redesign_all.py)."""
    a = pd.read_csv(HERE / "tag_sites_v5_all_r8.tsv", sep="\t")
    b = pd.read_csv(HERE / "tag_sites_v6_all_r8.tsv", sep="\t")
    r5 = a[(a["design_id"] == did) & (a["rank"] == 1)]
    r6 = b[(b["design_id"] == did) & (b["rank"] == 1)]
    if r5.empty or r6.empty:
        return []
    s5, s6 = int(r5["insert_after_residue"].iloc[0]), int(r6["insert_after_residue"].iloc[0])
    g5, g6 = region_of(sha, s5), region_of(sha, s6)
    if g5 is None or g6 is None or g5 == g6:
        return []
    if any(region_of(sha, m) == g6 for m in keep_sites):
        return []  # an arm already tests v6's region
    base = {"design_id": did, "gene": gene, "sequence_sha256": sha, "old_region": g5}
    if n_loops >= 3:
        return [{**base, "action": "flag_member", "new_region": g6,
                 "reason": f"v6's best region ({g6}, site {s6}) differs from v5's ({g5}); "
                 "member already has three arms"}]  # fmt: skip
    rep = best_clean(sha, g6, keep_sites, True)
    if rep is None:
        return [{**base, "action": "flag_member", "new_region": g6,
                 "reason": f"v6's best region ({g6}) differs from v5's ({g5}); no clean site "
                 "in it"}]  # fmt: skip
    return [{**base, "action": "add", "new_site": rep[0], "new_region": rep[1],
             "new_v6_rank_score": round(rep[2], 4),
             "reason": f"v6's best region ({g6}) differs from v5's ({g5}): a second model's "
             "hypothesis, added as an arm"}]  # fmt: skip


if __name__ == "__main__":
    raise SystemExit(main())
