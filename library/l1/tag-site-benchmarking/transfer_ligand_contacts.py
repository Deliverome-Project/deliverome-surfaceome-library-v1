#!/usr/bin/env python3
"""Ligand contacts transferred from structured relatives, for members that have none.

430 of the 566 L1 members have no D1 contact of their own: nobody has solved
them with a ligand bound. Many have a relative that has been solved (another
GPCR, transporter or channel of the same family), and a ligand pocket is
usually conserved across a family. This transfers those contacts by alignment,
and measures how often a transferred contact is real before anything uses it.

Donors are every ligand contact in D1's active contact release (all 1,680
accessions, not only ours), classified as in `fetch_interface_contacts.py`
(ligand class only: no antibodies, no partner proteins), with contexts that
cannot be extracellular dropped.

For each recipient:

1. Shortlist donors by shared 3-mers (Jaccard), keep the best TOP_DONORS.
2. Align recipient and donor locally (BLOSUM62, gap open -10, extend -0.5).
3. For every donor contact the alignment reaches, record the recipient residue,
   the alignment's overall identity, the identity of the +/-5 aligned window
   around the contact, and whether the residue itself is identical.

Validation (`ligand_transfer_validation.tsv`): proteins with their own ligand
contacts on DeepTMHMM v1 extracellular residues are scored with their own
accession, and any donor at >= SELF_IDENTITY overall identity, removed, so the
answer is never copied from the protein itself. Precision and recall are
measured on extracellular residues, over a sweep of identity thresholds.

Application (`member_transferred_ligand_contacts.tsv`): every L1 member's
ordered protein, with the donor that is the member's own canonical accession
excluded (its contacts are already placed exactly by
`map_contacts_to_members.py`). Each row keeps its alignment metrics, so a reader
can apply a stricter threshold; `passes_threshold` marks the rule chosen from
the validation.

    uv run --locked python transfer_ligand_contacts.py
"""

from __future__ import annotations

import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Align
from Bio.Align import substitution_matrices

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[3] / "src"))

from fetch_uniprot_site_features import ACCESSION  # noqa: E402
from fetch_uniprot_site_features import fetch as fetch_features  # noqa: E402
from map_contacts_to_members import ordered_proteins  # noqa: E402

DONORS = HERE / "data" / "d1_ligand_contacts_all.tsv"  # gitignored cache; S3 holds the copy
OUT_VAL = HERE / "ligand_transfer_validation.tsv"
OUT = HERE / "member_transferred_ligand_contacts.tsv"
OUT_JSON = HERE / "ligand_transfer.json"
TOP_DONORS = 8
WINDOW = 5
SELF_IDENTITY = 0.95  # in validation, a donor this close is the same protein
GLOBAL_GRID = (0.25, 0.35, 0.5, 0.7)
WINDOW_GRID = (0.0, 0.5, 0.7)
K = 3


def donor_contacts() -> pd.DataFrame:
    if DONORS.is_file():
        return pd.read_csv(DONORS, sep="\t")
    from deliverome_analysis.env import load_env

    load_env()
    from fetch_interface_contacts import classify
    from score_l1_library import d1

    rel = d1("SELECT release_id FROM contact_active_release")[0]["release_id"]
    rows, off = [], 0
    while True:
        chunk = d1(
            "SELECT uniprot_acc, source, context, observation_json FROM contact_observation "
            "WHERE release_id = ? AND context NOT LIKE 'non_extracellular%' "
            "ORDER BY observation_id LIMIT 2000 OFFSET ?",
            [rel, off],
        )
        for r in chunk:
            o = json.loads(r["observation_json"] or "{}")
            if classify(r["uniprot_acc"], r["source"], o) != "ligand":
                continue
            for p in o.get("positions", []):
                if isinstance(p, (int, float)):
                    rows.append(
                        {
                            "donor_acc": r["uniprot_acc"],
                            "donor_residue": int(p),
                            "source": r["source"],
                            "partner": str(o.get("partner") or "").replace("\t", " "),
                            "pdb": o.get("pdb") or "",
                        }
                    )
        off += len(chunk)
        if len(chunk) < 2000:
            break
    df = pd.DataFrame(rows).drop_duplicates()
    DONORS.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(DONORS, sep="\t", index=False)
    return df


def kmers(s: str) -> set[str]:
    return {s[i : i + K] for i in range(len(s) - K + 1)}


def aligner() -> Align.PairwiseAligner:
    al = Align.PairwiseAligner(mode="local")
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    al.open_gap_score, al.extend_gap_score = -10, -0.5
    return al


def transfer(rec_seq: str, donors: list[tuple[str, str, list[int]]], al) -> list[dict]:
    """Every donor contact the local alignment reaches, with its alignment metrics."""
    out = []
    for dacc, dseq, positions in donors:
        a = al.align(rec_seq, dseq)[0]
        pairs = []  # (rec_pos, donor_pos), 1-based, aligned columns only
        for (r0, r1), (d0, _d1) in zip(*a.aligned, strict=False):
            pairs += [(r0 + k + 1, d0 + k + 1) for k in range(r1 - r0)]
        if not pairs:
            continue
        same = np.array([rec_seq[r - 1] == dseq[d - 1] for r, d in pairs])
        ident = float(same.mean())
        cover = len(pairs) / len(rec_seq)
        d2i = {d: i for i, (_r, d) in enumerate(pairs)}
        for p in set(positions):
            i = d2i.get(p)
            if i is None:
                continue
            lo, hi = max(0, i - WINDOW), min(len(pairs), i + WINDOW + 1)
            out.append(
                {
                    "residue": pairs[i][0],
                    "donor_acc": dacc,
                    "donor_residue": p,
                    "alignment_identity": round(ident, 3),
                    "alignment_coverage": round(cover, 3),
                    "window_identity": round(float(same[lo:hi].mean()), 3),
                    "residue_identical": bool(same[i]),
                }
            )
    return out


def shortlist(seq: str, donor_k: dict[str, set[str]], exclude: set[str]) -> list[str]:
    q = kmers(seq)
    scored = [(len(q & k) / max(len(q | k), 1), a) for a, k in donor_k.items() if a not in exclude]
    scored.sort(reverse=True)
    return [a for s, a in scored[:TOP_DONORS] if s > 0.02]


def passes(df: pd.DataFrame, g: float, w: float) -> pd.Series:
    return (df["alignment_identity"] >= g) & (df["window_identity"] >= w)


def main() -> int:
    dc = donor_contacts()
    dc = dc[dc["donor_acc"].map(lambda a: bool(ACCESSION.match(str(a))))]
    print(f"donor ligand contacts: {len(dc)} residues on {dc['donor_acc'].nunique()} accessions")
    dpos = dc.groupby("donor_acc")["donor_residue"].apply(lambda v: sorted(set(int(x) for x in v)))
    dseq = {}
    for a in dpos.index:
        d = fetch_features(a) or {}
        if d.get("sequence"):
            dseq[a] = d["sequence"]
    donor_k = {a: kmers(s) for a, s in dseq.items()}
    al = aligner()

    # ---------------------------------------------------------------- validation
    allx = pd.read_csv(
        HERE / "all_extracellular_targets.tsv",
        sep="\t",
        usecols=["uniprot_acc", "residue"],
        low_memory=False,
    )
    ext = allx.groupby("uniprot_acc")["residue"].apply(lambda v: set(int(x) for x in v)).to_dict()
    own = {
        a: set(p) & ext.get(a, set())
        for a, p in dpos.items()
        if a in ext and set(p) & ext.get(a, set())
    }
    val = []
    for acc, truth in sorted(own.items()):
        seq = dseq.get(acc)
        if not seq:
            continue
        cands = shortlist(seq, donor_k, {acc})
        rows = transfer(seq, [(a, dseq[a], dpos[a]) for a in cands], al)
        rows = [r for r in rows if r["alignment_identity"] < SELF_IDENTITY]
        for r in rows:
            r["recipient"] = acc
            r["extracellular"] = r["residue"] in ext[acc]
            r["true_contact"] = r["residue"] in truth
        val += rows
    v = pd.DataFrame(val)
    sweep = []
    n_truth = sum(len(t) for a, t in own.items() if a in dseq)
    n_ext = sum(len(ext[a]) for a in own if a in dseq)
    base_rate = n_truth / n_ext
    for g in GLOBAL_GRID:
        for w in WINDOW_GRID:
            sel = v[passes(v, g, w) & v["extracellular"]].drop_duplicates(["recipient", "residue"])
            tp = int(sel["true_contact"].sum())
            sweep.append(
                {
                    "min_alignment_identity": g,
                    "min_window_identity": w,
                    "predicted_residues": len(sel),
                    "true_positive": tp,
                    "precision": round(tp / max(len(sel), 1), 3),
                    "recall": round(tp / max(n_truth, 1), 3),
                    "proteins_with_prediction": int(sel["recipient"].nunique()),
                    "base_rate": round(base_rate, 4),
                }
            )
    sw = pd.DataFrame(sweep)
    sw.to_csv(OUT_VAL, sep="\t", index=False)
    print(f"validation: {len(own)} proteins, {n_truth} true extracellular ligand contacts")
    print(sw.to_string(index=False))
    # the loosest rule whose precision is at least half
    ok = sw[sw["precision"] >= 0.5].sort_values("recall", ascending=False)
    chosen = ok.iloc[0] if len(ok) else sw.sort_values("precision", ascending=False).iloc[0]
    g, w = float(chosen["min_alignment_identity"]), float(chosen["min_window_identity"])
    print(f"chosen: identity >= {g}, window identity >= {w}")

    # ---------------------------------------------------------------- application
    members = pd.read_csv(HERE / "member_interface_mapping.tsv", sep="\t")
    recon = pd.read_csv(HERE / "l1_member_reconciliation.tsv", sep="\t")
    prot = ordered_proteins()
    seq_of = dict(zip(recon["candidate_id"], recon["sequence_sha256"].map(prot), strict=False))
    out = []
    for m in members.itertuples(index=False):
        seq = seq_of[m.candidate_id]
        excl = {m.canonical_acc} if isinstance(m.canonical_acc, str) else set()
        cands = shortlist(seq, donor_k, excl)
        for r in transfer(seq, [(a, dseq[a], dpos[a]) for a in cands], al):
            out.append({"candidate_id": m.candidate_id, "member_acc": m.member_acc, **r})
    t = pd.DataFrame(out)
    t["passes_threshold"] = passes(t, g, w)
    meta = dc.drop_duplicates(["donor_acc", "donor_residue"]).set_index(
        ["donor_acc", "donor_residue"]
    )[["partner", "pdb"]]
    t = t.join(meta, on=["donor_acc", "donor_residue"])
    t.to_csv(OUT, sep="\t", index=False)
    best = t[t["passes_threshold"]]
    by_member = defaultdict(set)
    for c, r in zip(best["candidate_id"], best["residue"], strict=False):
        by_member[c].add(r)
    OUT_JSON.write_text(
        json.dumps(
            {
                "donors": {"residues": len(dc), "accessions": int(dc["donor_acc"].nunique())},
                "top_donors": TOP_DONORS,
                "window": WINDOW,
                "validation_proteins": len(own),
                "validation_true_contacts": n_truth,
                "base_rate": round(base_rate, 4),
                "chosen_rule": {"min_alignment_identity": g, "min_window_identity": w},
                "members_with_transferred_contacts": len(by_member),
                "transferred_residues": int(sum(len(v) for v in by_member.values())),
            },
            indent=2,
        )
        + "\n"
    )
    print(json.loads(OUT_JSON.read_text()))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
