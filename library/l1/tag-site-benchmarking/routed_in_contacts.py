#!/usr/bin/env python3
"""Ligand contacts for members routed into the internal rung after the contacts were built.

`member_interface_contacts.tsv` and `member_transferred_ligand_contacts.tsv` cover the internal
rung as it stood on 2026-11-07. A member routed in later (`l1_internal_decisions.ROUTED_IN`)
gets the same two sources here, built the same way, so the ±2 contact veto treats it like any
other member:

    d1_ligand            its own D1 ligand contacts (the donor table of transfer_ligand_contacts,
                         every ligand-class contact in D1's active release), used only where the
                         ordered protein is the accession's own sequence (identical), so the
                         numbering is the same
    transferred_ligand   contacts transferred from relatives by transfer_ligand_contacts's own
                         shortlist, alignment and threshold (alignment identity >= 0.5, the rule
                         that file chose: precision 0.68 on held-out proteins), its own
                         accession excluded as a donor

    uv run --locked python routed_in_contacts.py      # writes routed_in_ligand_contacts.tsv
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[3] / "src"))

import transfer_ligand_contacts as T  # noqa: E402
from l1_internal_decisions import ROUTED_IN, full_head_paths, head_tables  # noqa: E402

#: transfer_ligand_contacts's chosen rule (ligand_transfer_validation.tsv row 6)
MIN_IDENTITY, MIN_WINDOW = 0.5, 0.0
OUT = HERE / "routed_in_ligand_contacts.tsv"


def contact_rows(members) -> list[dict]:
    """Own and transferred ligand contacts for (design_id, gene, sha256, accession, sequence)."""
    dc = T.donor_contacts()
    dc = dc[dc["donor_acc"].map(lambda a: bool(T.ACCESSION.match(str(a))))]
    dpos = dc.groupby("donor_acc")["donor_residue"].apply(lambda v: sorted({int(x) for x in v}))
    dseq = {}
    for a in dpos.index:
        d = T.fetch_features(a) or {}
        if d.get("sequence"):
            dseq[a] = d["sequence"]
    donor_k = {a: T.kmers(s) for a, s in dseq.items()}
    al = T.aligner()
    rows = []
    for did, gene, sha, acc, seq in members:
        # a designed ORF (l1-p...) has no accession, so no own D1 contacts
        own = (T.fetch_features(acc) or {}).get("sequence") if T.ACCESSION.match(acc) else None
        if own == seq:
            for res in dpos.get(acc, []):
                rows.append(
                    {"design_id": did, "gene": gene, "sequence_sha256": sha, "residue": res,
                     "source": "d1_ligand", "donor_acc": acc, "alignment_identity": 1.0}
                )  # fmt: skip
        cands = T.shortlist(seq, donor_k, {acc})
        t = pd.DataFrame(T.transfer(seq, [(a, dseq[a], dpos[a]) for a in cands], al))
        if len(t):
            t = t[T.passes(t, MIN_IDENTITY, MIN_WINDOW)]
            for x in t.itertuples(index=False):
                rows.append(
                    {"design_id": did, "gene": gene, "sequence_sha256": sha,
                     "residue": int(x.residue), "source": "transferred_ligand",
                     "donor_acc": x.donor_acc, "alignment_identity": x.alignment_identity}
                )  # fmt: skip
        print(f"{gene} ({did}): {len(dpos.get(acc, []))} own D1 ligand contacts "
              f"(sequence identical: {own == seq}); donors shortlisted {len(cands)}; "
              f"transferred passing {len(t) if len(t) else 0}")  # fmt: skip
    return rows


COLS = ["design_id", "gene", "sequence_sha256", "residue", "source", "donor_acc",
        "alignment_identity"]  # fmt: skip


def main() -> int:
    _internal, _ids, prot = head_tables()
    tp = full_head_paths()
    members = []
    for did, r in ROUTED_IN.items():
        sha = tp[tp["design_id"] == did].iloc[0]["sequence_sha256"]
        members.append((did, r["gene"], sha, r["sequence_accession"], prot[sha]))
    pd.DataFrame(contact_rows(members), columns=COLS).drop_duplicates(
        ["sequence_sha256", "residue", "source"]
    ).to_csv(OUT, sep="\t", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
