#!/usr/bin/env python3
"""One row per extracellular segment of every L1 internal-tag member, with its evidence.

Meant to sit in front of whoever makes the final loop and residue choice (a
scientist or an agent): the ranker's opinion of each segment, and next to it
every piece of negative evidence present in that segment, as explicit columns
rather than folded into a score.

The scores come from one arm of `ranker_pu_segment.py`, refitted on all 176
published sites (default: the recommended arm). Members that carry a published
site are scored in-sample; `in_published_corpus` marks them.

    l1_loop_evidence.tsv          one row per (member, extracellular segment)
    l1_loop_evidence.json         column dictionary, rules and provenance

    uv run --locked python l1_loop_evidence.py [--arm pu_bagged]
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

from independent_model_search import make, smooth, within_protein  # noqa: E402
from ranker_pu_segment import (  # noqa: E402
    ARMS,
    HARD_PLDDT,
    HARD_RSA,
    KIND,
    SEED,
    TWO_STAGE_BASE,
    annotate_constraints,
    frame,
    lexicographic,
    load_screens,
    pu_fit_predict,
    pu_frame,
    region_fit_predict,
    segment_tables,
    training_set,
)
from tag_site_model_v4 import FEATS  # noqa: E402

OUT = HERE / "l1_loop_evidence.tsv"
META = HERE / "l1_loop_evidence.json"
TM_ADJACENT = 3  # residues from a TM helix end; an insert here pulls on the helix
LOW_RSA = 0.25  # relative solvent accessibility below which a residue counts as poorly exposed

COLUMNS = {
    "member_id": "L1 member id after the 2026-11-07 renumbering (l1-m...), joined by the ordered "
    "protein's sequence_sha256; use this one",
    "candidate_id": "the pre-renumbering id (l1-p...) the S5E13 tables carry",
    "gene": "gene symbol",
    "uniprot_acc": "accession whose numbering every position in this row uses",
    "candidate_set": "library_reference / ortholog / isoform",
    "in_published_corpus": "this accession carries a published site, so its scores are in-sample",
    "segment_index": "extracellular segment (maximal DeepTMHMM v1 'O' run), 0-based along the chain",
    "start": "first residue of the segment",
    "end": "last residue of the segment",
    "length": "residues in the segment",
    "n_segments": "extracellular segments in this protein",
    "segment_score": "arm's score for the segment (two-stage arms: stage-1 probability; "
    "otherwise the best residue score in the segment)",
    "segment_rank": "rank of this segment within the protein by segment_score (1 = ranker's choice)",
    "best_residue": "highest-scoring residue in the segment; the tag goes after it",
    "best_residue_score": "its within-protein score",
    "best_residue_percentile": "its percentile among all extracellular residues of the protein",
    "published_sites": "published insertion sites in this segment (corpus members only)",
    "tm_adjacent": f"residues within {TM_ADJACENT} of a TM helix end, as ranges",
    "disulfide_cys": "UniProt DISULFID cysteines in the segment, as Cys>partner "
    "('interchain' when the partner is on another chain)",
    "glycan_n": "UniProt N-linked glycosylated asparagines in the segment",
    "sequon_n": "N of every N-X-S/T sequon (X != P) in the segment, from the sequence",
    "functional_sites": "UniProt BINDING / ACT_SITE / METAL / SITE positions in the segment, with type",
    "ligand_contacts": "residues contacting a small molecule, peptide or endogenous ligand in a "
    "deposited structure (D1 contact_observation), moved onto the ordered protein by sequence "
    "alignment (identical residues only), as ranges; counted as a hard negative",
    "partner_interface": "residues contacting another protein (co-receptor, hetero-complex) in a "
    "deposited structure, from the same D1 table, as ranges",
    "antibody_epitopes": "residues an antibody, Fab, VHH or scFv contacts: exposed surface the "
    "protein does not need, shown for reference and not counted as negative evidence",
    "ligand_contacts_transferred": "PREDICTED ligand contacts: a relative's structure-derived ligand "
    "contact moved across by local alignment at >= 50% identity (transfer_ligand_contacts.py; "
    "precision 0.68, recall 0.37 on 72 proteins with their own structures hidden), as "
    "residue<donor accession; not counted in hard_negative_fraction",
    "screen_intolerant": "junctions a domain-insertion screen calls intolerant, as ranges "
    "(screened proteins only; blank elsewhere means not screened, not tolerant)",
    "median_plddt": "median AlphaFold pLDDT over the segment",
    "max_disorder": "maximum metapredict disorder over the segment",
    "membrane_side_eligible": "junctions the membrane-side mask calls extracellular (insert-after "
    "residue n with mask[n-1] and mask[n] both 'E'; membrane_mask.py, s4e9 gate PR #291), as "
    "ranges; 'NO USABLE MASK' when the member has none, in which case nothing should be designed",
    "buried_fraction": f"fraction of residues with RSA < {HARD_RSA} and pLDDT >= {HARD_PLDDT}",
    "low_rsa_fraction": f"fraction of residues with RSA < {LOW_RSA}",
    "hard_negative_fraction": "fraction of residues carrying any hard-negative reason "
    "(buried, disulfide Cys, glycan N, sequon N, screen-intolerant, ligand contact)",
}


def ranges(vals) -> str:
    v = sorted({int(x) for x in vals})
    if not v:
        return ""
    out, a, b = [], v[0], v[0]
    for x in v[1:]:
        if x == b + 1:
            b = x
            continue
        out.append(f"{a}" if a == b else f"{a}-{b}")
        a = b = x
    out.append(f"{a}" if a == b else f"{a}-{b}")
    return ",".join(out)


def member_rows(allx: pd.DataFrame) -> pd.DataFrame:
    """allx plus the non-canonical proteoforms, as tag_site_model_v4.design reads them."""
    e = pd.read_csv(HERE / "l1_member_extracellular.tsv", sep="\t", low_memory=False)
    for c in ("ss3", "topology_class", "region_class"):
        e[c] = e[c].fillna("").astype(str)
    e["gene_symbol"] = e["gene_symbol"].fillna(e["uniprot_acc"]).astype(str)
    cols = [c for c in allx.columns if c in e.columns]
    rows = pd.concat([allx[cols], e[cols]], ignore_index=True)
    return rows.drop_duplicates(["uniprot_acc", "residue"], keep="first").reset_index(drop=True)


def score(arm: str, pos, allx, rows) -> np.ndarray:
    """The arm refitted on every published site, applied to `rows`."""
    feats = [f for f in FEATS if f in pos.columns]
    base = TWO_STAGE_BASE.get(arm, arm)
    rz = within_protein(rows, feats, "z")
    if base == "pu_bagged":
        full = pu_frame(pos, allx)
        (p,) = pu_fit_predict(
            within_protein(full, feats, "z"),
            full["label"].to_numpy(),
            [rz],
            feats,
            np.random.default_rng(SEED),
        )
    else:
        tr = training_set(base, pos, allx)
        m = make(KIND, feats)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            m.fit(within_protein(tr, feats, "z"), tr["label"].to_numpy())
            p = m.predict_proba(rz)[:, 1]
    return smooth(rows, p, 2)


def stage1(arm: str, pos, allx, rows) -> dict:
    if arm == "two_stage":
        tr = allx[allx["uniprot_acc"].isin(set(pos["uniprot_acc"]))]
        tr_reg, rf = segment_tables(pos, tr)
    else:
        from ranker_pu_segment import SCREENED_INSERTION

        prots = set(pos["uniprot_acc"]) | set(SCREENED_INSERTION)
        tr_reg, rf = segment_tables(pos, allx[allx["uniprot_acc"].isin(prots)], hard_labels=True)
        tr_reg = tr_reg[tr_reg["label"].notna()]
    out = {}
    for acc, g in rows.groupby("uniprot_acc"):
        ev_reg, _ = segment_tables(pos, g)
        for ri, v in region_fit_predict(tr_reg, ev_reg, rf).items():
            out[(acc, ri)] = v
    return out


TAG_PATHS_KEY = (
    "sprint4/s4e9-orf-design-fg-panel/library-design/member-renumbering/2026-11-07/l1_tag_paths.tsv"
)


def tag_paths() -> pd.DataFrame:
    """The 2026-11-07 tag-path table: current member ids, and every internal-placement member."""
    local = HERE / "data" / "l1_tag_paths_2026-11-07.tsv"
    if not local.is_file():
        from deliverome_analysis import s3

        local.parent.mkdir(parents=True, exist_ok=True)
        s3.processed.download(TAG_PATHS_KEY, str(local))
    return pd.read_csv(local, sep="\t", low_memory=False)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--arm", choices=ARMS, default="pu_bagged")
    args = ap.parse_args()

    pos, allx_raw, _, _ = frame()
    screens = load_screens()
    allx = annotate_constraints(allx_raw, screens)
    targets = pd.read_csv(HERE / "l1_member_targets.tsv", sep="\t")
    rows = annotate_constraints(member_rows(allx_raw), screens)
    rows = rows[rows["uniprot_acc"].isin(set(targets["uniprot_acc"]))].reset_index(drop=True)
    import membrane_mask

    member_mask, _refused = membrane_mask.masks_for(sorted(set(rows["uniprot_acc"])))
    # contacts already moved onto each member's ordered protein (map_contacts_to_members.py)
    from map_contacts_to_members import contact_keys

    contact_sets = {c: contact_keys({c}) for c in ("ligand", "protein_partner", "antibody")}
    tr = pd.read_csv(HERE / "member_transferred_ligand_contacts.tsv", sep="\t")
    tr = tr[tr["passes_threshold"]]
    transferred: dict[tuple[str, int], set[str]] = {}
    for a, r, d in zip(tr["member_acc"], tr["residue"], tr["donor_acc"], strict=False):
        transferred.setdefault((a, int(r)), set()).add(d)
    key = zip(rows["uniprot_acc"], rows["residue"].astype(int), strict=False)
    rows["hn_ligand_contact"] = [k in contact_sets["ligand"] for k in key]
    rows["hard_negative"] = rows["hard_negative"] | rows["hn_ligand_contact"]
    missing = sorted(set(targets["uniprot_acc"]) - set(rows["uniprot_acc"]))
    print(f"{rows['uniprot_acc'].nunique()} members with extracellular rows; missing {missing}")

    rows["score"] = score(args.arm, pos, allx, rows)
    if args.arm in TWO_STAGE_BASE:
        seg = stage1(args.arm, pos, allx, rows)
        segv = np.array(
            [
                seg.get(k, np.nan)
                for k in zip(rows["uniprot_acc"], rows["region_index"], strict=False)
            ]
        )
        rows["order"] = lexicographic(rows, segv, rows["score"].to_numpy())
        rows["seg_stage1"] = segv
    else:
        rows["order"] = rows["score"]
    rows["pct"] = rows.groupby("uniprot_acc")["order"].rank(pct=True) * 100

    pub = pos.groupby("uniprot_acc")["residue"].apply(lambda v: sorted(int(x) for x in v)).to_dict()
    out = []
    for t in targets.itertuples(index=False):
        g = rows[rows["uniprot_acc"] == t.uniprot_acc]
        if g.empty:
            continue
        segs = []
        for ri, s in g.groupby("region_index"):
            best = s.loc[s["order"].idxmax()]
            dis = s[s["hn_disulfide_cys"]]
            segs.append(
                {
                    "candidate_id": t.candidate_id,
                    "gene": t.gene,
                    "uniprot_acc": t.uniprot_acc,
                    "candidate_set": t.candidate_set,
                    "in_published_corpus": t.uniprot_acc in pub,
                    "segment_index": int(ri),
                    "start": int(s["residue"].min()),
                    "end": int(s["residue"].max()),
                    "length": len(s),
                    "n_segments": int(g["region_index"].nunique()),
                    "segment_score": round(
                        float(s["seg_stage1"].iloc[0] if "seg_stage1" in s else s["score"].max()),
                        4,
                    ),
                    "best_residue": int(best["residue"]),
                    "best_residue_score": round(float(best["score"]), 4),
                    "best_residue_percentile": round(float(best["pct"]), 1),
                    "published_sites": ",".join(
                        str(p)
                        for p in pub.get(t.uniprot_acc, [])
                        if s["residue"].min() <= p <= s["residue"].max()
                    ),
                    "tm_adjacent": ranges(s.loc[s["dist_to_tm"] <= TM_ADJACENT, "residue"]),
                    "disulfide_cys": ",".join(
                        f"C{int(r)}>{p if p == 'interchain' else 'C' + p}"
                        for r, p in zip(dis["residue"], dis["disulfide_partner"], strict=False)
                    ),
                    "glycan_n": ",".join(str(int(r)) for r in s.loc[s["hn_glycan_n"], "residue"]),
                    "sequon_n": ",".join(str(int(r)) for r in s.loc[s["hn_sequon_n"], "residue"]),
                    "functional_sites": ";".join(
                        f"{int(r)} {f}"
                        for r, f in zip(s["residue"], s["functional_site"], strict=False)
                        if f
                    ),
                    **{
                        col: ranges(
                            r
                            for r in s["residue"]
                            if (t.uniprot_acc, int(r)) in contact_sets.get(cls, set())
                        )
                        for col, cls in (
                            ("ligand_contacts", "ligand"),
                            ("partner_interface", "protein_partner"),
                            ("antibody_epitopes", "antibody"),
                        )
                    },
                    "ligand_contacts_transferred": ",".join(
                        f"{int(r)}<{'/'.join(sorted(transferred[(t.uniprot_acc, int(r))]))}"
                        for r in s["residue"]
                        if (t.uniprot_acc, int(r)) in transferred
                    ),
                    "screen_intolerant": ranges(s.loc[s["hn_screen_intolerant"], "residue"]),
                    "buried_fraction": round(float(s["hn_buried"].mean()), 3),
                    "low_rsa_fraction": round(float((s["rsa"] < LOW_RSA).mean()), 3),
                    "hard_negative_fraction": round(float(s["hard_negative"].mean()), 3),
                    "membrane_side_eligible": (
                        ranges(
                            r
                            for r, e in zip(
                                s["residue"],
                                membrane_mask.eligible_after(
                                    member_mask[t.uniprot_acc], s["residue"]
                                ),
                                strict=False,
                            )
                            if e
                        )
                        if t.uniprot_acc in member_mask
                        else "NO USABLE MASK"
                    ),
                    "median_plddt": round(float(s["plddt"].median()), 1),
                    "max_disorder": round(float(s["disorder"].max()), 3),
                }
            )
        order = sorted(range(len(segs)), key=lambda i: -segs[i]["segment_score"])
        for rank, i in enumerate(order, 1):
            segs[i]["segment_rank"] = rank
        out += segs
    tp = tag_paths()
    recon = pd.read_csv(HERE / "l1_member_reconciliation.tsv", sep="\t")
    new_id = recon.set_index("candidate_id")["sequence_sha256"].map(
        tp.drop_duplicates("sequence_sha256").set_index("sequence_sha256")["candidate_id"]
    )
    for r in out:
        r["member_id"] = new_id.get(r["candidate_id"])
    internal = tp[tp["placement"] == "internal"]
    unscored = internal[~internal["sequence_sha256"].isin(set(recon["sequence_sha256"]))]
    df = pd.DataFrame(out)[list(COLUMNS)]
    df.to_csv(OUT, sep="\t", index=False)
    META.write_text(
        json.dumps(
            {
                "arm": args.arm,
                "fitted_on": f"{len(pos)} published sites, {pos['uniprot_acc'].nunique()} proteins",
                "members": int(df["candidate_id"].nunique()),
                "internal_members_in_tag_paths_2026_11_07": len(internal),
                "internal_members_not_scored": unscored[
                    ["candidate_id", "gene_symbol", "topology_cohort"]
                ].to_dict(orient="records"),
                "not_scored_reason": "absent from l1_internal_members.tsv, so no S5E13 per-residue "
                "feature rows were ever exported for it",
                "segments": len(df),
                "columns": COLUMNS,
            },
            indent=2,
        )
        + "\n"
    )
    print(f"wrote {OUT.name}: {len(df)} segments on {df['candidate_id'].nunique()} members")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
