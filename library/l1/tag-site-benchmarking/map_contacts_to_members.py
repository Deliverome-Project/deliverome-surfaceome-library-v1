#!/usr/bin/env python3
"""Put D1 interface contacts onto the protein each L1 member actually orders.

D1's contacts (`interface_contacts.tsv`) are positions in the UniProt numbering
of a human canonical accession. An L1 member's per-residue rows, scores and
designs are numbered on its ORDERED protein: the designed CDS translated, read
from `candidate_cds.json` and identified by the member's `sequence_sha256`. For
452 of 566 members that protein is the canonical sequence itself. For the rest
(isoforms, cyno and mouse orthologs, designed ORFs) it is a different protein,
and a position copied across by number lands on the wrong residue.

So every member is mapped the same way, by sequence and never by number:

1. Find the gene's human canonical accession: the member accession for
   `human_canonical`; the accession without its isoform suffix for
   `human_isoforms`; otherwise the cohort's canonical accession for the
   resolver-approved symbol (`gene_resolver.to_canonical`), and failing that the
   reviewed human UniProt entry for the gene (EBI Proteins API).
2. Align the canonical sequence to the ordered protein (global, BLOSUM62, gap
   open -10, extend -0.5).
3. Transfer a contact only where the two aligned residues are IDENTICAL. A
   contact on a residue the ordered protein lacks or has changed is dropped and
   counted, not moved to a neighbor.

Writes:

    member_interface_contacts.tsv   candidate_id, member_acc (the per-residue frame),
                                    residue (member numbering), canonical_acc,
                                    canonical_residue, class, source, context,
                                    partner, pdb
    member_interface_mapping.tsv    per member: frame, canonical accession and how it
                                    was found, identity, contacts transferred / dropped

    uv run --locked python map_contacts_to_members.py
"""

from __future__ import annotations

import hashlib
import json
import sys
import urllib.request
from pathlib import Path

import pandas as pd
from Bio import Align
from Bio.Align import substitution_matrices
from Bio.Seq import Seq

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[3] / "src"))

from fetch_uniprot_site_features import ACCESSION  # noqa: E402
from fetch_uniprot_site_features import fetch as fetch_features  # noqa: E402

from deliverome_analysis.gene_resolver import to_canonical  # noqa: E402

CDS_KEY = "sprint4/s4e9-orf-design-fg-panel/library-candidates/candidate_cds.json"
OUT = HERE / "member_interface_contacts.tsv"
OUT_MAP = HERE / "member_interface_mapping.tsv"
SEARCH = "https://www.ebi.ac.uk/proteins/api/proteins?gene={g}&taxid=9606&reviewed=true"


def ordered_proteins() -> dict[str, str]:
    local = HERE / "data" / "candidate_cds.json"
    if not local.is_file():
        from deliverome_analysis import s3

        local.parent.mkdir(parents=True, exist_ok=True)
        s3.processed.download(CDS_KEY, str(local))
    out = {}
    for v in json.loads(local.read_text()).values():
        if v.get("cds"):
            p = str(Seq(v["cds"]).translate()).rstrip("*")
            out[hashlib.sha256(p.encode()).hexdigest()] = p
    return out


def reviewed_human(gene: str) -> str | None:
    req = urllib.request.Request(
        SEARCH.format(g=gene),
        headers={"Accept": "application/json", "User-Agent": "deliverome-analysis s5e13"},
    )
    with urllib.request.urlopen(req, timeout=60) as r:
        hits = json.loads(r.read())
    exact = [
        h["accession"]
        for h in hits
        if "-" not in h["accession"]  # the search also returns isoform entries
        and any(
            n.get("value") == gene
            for g in h.get("gene", [])
            for n in [g.get("name", {}), *g.get("synonyms", [])]
        )
    ]
    return exact[0] if len(exact) == 1 else None


def canonical_accessions(members: pd.DataFrame) -> pd.DataFrame:
    allx = pd.read_csv(
        HERE / "all_extracellular_targets.tsv",
        sep="\t",
        usecols=["gene_symbol", "uniprot_acc"],
        low_memory=False,
    ).drop_duplicates()
    hc = members.loc[
        members["topology_cohort"] == "human_canonical", ["gene_symbol", "topology_acc"]
    ]
    pool = pd.concat([hc.rename(columns={"topology_acc": "uniprot_acc"}), allx]).dropna()
    # human only: the cohort table also carries the cyno and mouse member accessions
    nonhuman = set(members.loc[members["species"] != "human", "topology_acc"])
    pool = pool[
        pool["uniprot_acc"].map(lambda a: bool(ACCESSION.match(a)))
        & ~pool["uniprot_acc"].isin(nonhuman)
    ]
    sym = to_canonical(sorted(set(members["gene_symbol"]) | set(pool["gene_symbol"])))
    pool["sym"] = pool["gene_symbol"].map(sym)
    by_sym = pool.groupby("sym")["uniprot_acc"].agg(lambda v: sorted(set(v)))
    acc, how = [], []
    for m in members.itertuples(index=False):
        if m.topology_cohort == "human_canonical":
            acc.append(m.topology_acc)
            how.append("member accession")
        elif m.topology_cohort == "human_isoforms":
            acc.append(m.topology_acc.split("-")[0])
            how.append("isoform base accession")
        elif (
            isinstance(by_sym.get(sym.get(m.gene_symbol)), list)
            and len(by_sym[sym[m.gene_symbol]]) == 1
        ):
            acc.append(by_sym[sym[m.gene_symbol]][0])
            how.append("cohort canonical for approved symbol")
        else:
            a = reviewed_human(sym.get(m.gene_symbol) or m.gene_symbol)
            acc.append(a)
            how.append("reviewed human UniProt entry" if a else "not found")
    return members.assign(canonical_acc=acc, canonical_source=how)


def aligner() -> Align.PairwiseAligner:
    al = Align.PairwiseAligner(mode="global")
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    al.open_gap_score, al.extend_gap_score = -10, -0.5
    return al


def position_map(src: str, dst: str, al) -> tuple[dict[int, int], float]:
    """1-based src -> dst positions where the aligned residues are identical."""
    if src == dst:
        return {i: i for i in range(1, len(src) + 1)}, 1.0
    a = al.align(src, dst)[0]
    out = {}
    for (s0, s1), (d0, _d1) in zip(*a.aligned, strict=False):
        for k in range(s1 - s0):
            if src[s0 + k] == dst[d0 + k]:
                out[s0 + k + 1] = d0 + k + 1
    return out, len(out) / max(len(dst), 1)


def contact_keys(classes: set[str]) -> set[tuple[str, int]]:
    """(accession, residue) for contacts of `classes`, in every frame this folder uses.

    Canonical keys come from `interface_contacts.tsv` (the published corpus and the
    canonical L1 members are numbered on them); member keys come from
    `member_interface_contacts.tsv`, already moved onto each member's ordered
    protein. A canonical accession never names a non-canonical frame, so the union
    cannot put a contact on the wrong protein.
    """
    keys: set[tuple[str, int]] = set()
    for f, acc in ((HERE / "interface_contacts.tsv", "uniprot_acc"), (OUT, "member_acc")):
        d = pd.read_csv(f, sep="\t")
        d = d[d["class"].isin(classes)]
        keys |= set(zip(d[acc], d["residue"].astype(int), strict=False))
    return keys


def main() -> int:
    members = pd.read_csv(HERE / "l1_member_reconciliation.tsv", sep="\t")
    prot = ordered_proteins()
    members["ordered"] = members["sequence_sha256"].map(prot)
    assert members["ordered"].notna().all(), "an ordered protein is missing from candidate_cds.json"
    members = canonical_accessions(members)
    # the per-residue frame: member rows are keyed by the l1 target accession
    targets = pd.read_csv(HERE / "l1_member_targets.tsv", sep="\t")
    frame = dict(zip(targets["candidate_id"], targets["uniprot_acc"], strict=False))
    ic = pd.read_csv(HERE / "interface_contacts.tsv", sep="\t")
    by_acc = {a: g for a, g in ic.groupby("uniprot_acc")}
    al = aligner()
    rows, summary = [], []
    for m in members.itertuples(index=False):
        rec = {
            "candidate_id": m.candidate_id,
            "gene": m.gene_symbol,
            "topology_cohort": m.topology_cohort,
            "member_acc": frame[m.candidate_id],
            "canonical_acc": m.canonical_acc,
            "canonical_source": m.canonical_source,
        }
        contacts = by_acc.get(m.canonical_acc)
        if not isinstance(m.canonical_acc, str):
            summary.append({**rec, "status": "no canonical accession"})
            continue
        d = fetch_features(m.canonical_acc) or {}
        src = d.get("sequence", "")
        if not src:
            summary.append({**rec, "status": "no canonical sequence"})
            continue
        pmap, ident = position_map(src, m.ordered, al)
        n_in = 0 if contacts is None else len(contacts.drop_duplicates(["residue", "class"]))
        moved = 0
        if contacts is not None:
            for c in contacts.itertuples(index=False):
                t = pmap.get(int(c.residue))
                if t is None:
                    continue
                moved += 1
                rows.append(
                    {
                        "candidate_id": m.candidate_id,
                        "member_acc": frame[m.candidate_id],
                        "residue": t,
                        "canonical_acc": m.canonical_acc,
                        "canonical_residue": int(c.residue),
                        "class": c[ic.columns.get_loc("class")],
                        "source": c.source,
                        "context": c.context,
                        "partner": c.partner,
                        "pdb": c.pdb,
                    }
                )
        summary.append(
            {
                **rec,
                "status": "identical sequence" if src == m.ordered else "aligned",
                "canonical_length": len(src),
                "ordered_length": len(m.ordered),
                "identical_residue_fraction": round(ident, 4),
                "contact_observations": 0 if contacts is None else len(contacts),
                "contact_residue_classes": n_in,
                "transferred": moved,
                "dropped": (0 if contacts is None else len(contacts)) - moved,
            }
        )
    out = pd.DataFrame(rows).drop_duplicates()
    out.to_csv(OUT, sep="\t", index=False)
    sm = pd.DataFrame(summary)
    sm.to_csv(OUT_MAP, sep="\t", index=False)
    print(sm.groupby(["topology_cohort", "status"]).size().to_string())
    print(sm.groupby("canonical_source").size().to_string())
    print(
        f"contacts transferred {int(sm['transferred'].sum())}, dropped {int(sm['dropped'].sum())}; "
        f"members with any: {out['candidate_id'].nunique()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
