#!/usr/bin/env python3
"""The internal-loop refold queue for l1-tag-refold, capped at its budget.

l1-tag-refold asked for about one site per member, and two for the top 100. This takes:

    rank 1          every designed internal-rung member: the 532 in
                    `l1_internal_proposed_sites.tsv`, plus the 27 newly adopted orthologs
                    designed on PR #298 (`tag_sites_v5_orthologs.tsv` at ORTHOLOG_COMMIT)
    rank 2          the 100 members where the model is least sure which loop to use, by the
                    runner-up loop's ceiling as a fraction of the best loop's; a second refold
                    there decides the most

Each row carries the member id, the insertion residue in the member's designed protein, the
rank, and the SHA-256 of the tagged protein (51-aa cassette inserted after that residue).

    uv run --locked python refold_queue.py        # writes l1_internal_refold_queue.tsv
"""

from __future__ import annotations

import hashlib
import io
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[3] / "src"))

from build_internal_tagged import CASSETTE_AA  # noqa: E402
from map_contacts_to_members import ordered_proteins  # noqa: E402

#: PR #298's ortholog designs at the commit the committed queue records ("PR #298 @ 31b3ed6f"),
#: which is on main; the session branch it was first read from no longer exists
ORTHOLOG_COMMIT = "31b3ed6f"
ORTHOLOG_PATH = "library/l1/tag-site-benchmarking/tag_sites_v5_orthologs.tsv"
TOP_N_SECOND = 100
OUT = HERE / "l1_internal_refold_queue.tsv"


def ortholog_designs() -> tuple[pd.DataFrame, str]:
    from deliverome_analysis.provenance import read_pinned_text

    repo = HERE.parents[2]
    txt = read_pinned_text(ORTHOLOG_COMMIT, ORTHOLOG_PATH, repo=repo)
    return pd.read_csv(io.StringIO(txt), sep="\t"), ORTHOLOG_COMMIT


def current_ordered_proteins() -> dict[str, str]:
    """The adopted panel's designed-CDS record (2026-11-10 ortholog CDS policy, the one
    PR #283 and the membrane-side gate pin), translated, keyed by sha256."""
    import json

    from Bio.Seq import Seq

    local = HERE / "data" / "candidate_cds_current.json"
    if not local.is_file():
        from deliverome_analysis import s3

        s3.processed.download(
            "sprint4/s4e9-orf-design-fg-panel/library-candidates/"
            "2026-11-10-ortholog-cds-policy/candidate_cds.json",
            str(local),
        )
    out = {}
    for v in json.loads(local.read_text()).values():
        if v.get("cds"):
            q = str(Seq(v["cds"]).translate()).rstrip("*")
            out[hashlib.sha256(q.encode()).hexdigest()] = q
    return out


def main() -> int:
    mine = pd.read_csv(HERE / "l1_internal_proposed_sites.tsv", sep="\t")
    orth, orth_commit = ortholog_designs()
    cols = [
        "member_id",
        "sequence_sha256",
        "gene",
        "rank",
        "insert_after_residue",
        "region_ceiling",
        "runner_up_region_ceiling",
        "snorkel_candidate",
    ]
    designs = pd.read_csv(HERE / "tag_sites_v5_adopted.tsv", sep="\t")
    mine = mine.merge(
        designs[["uniprot_acc", "rank", "region_ceiling", "runner_up_region_ceiling"]],
        left_on=["frame_accession", "rank"],
        right_on=["uniprot_acc", "rank"],
    )
    allsites = pd.concat(
        [mine[cols].assign(source="PR #245"), orth[cols].assign(source=f"PR #298 @ {orth_commit}")],
        ignore_index=True,
    )
    r1 = allsites[allsites["rank"] == 1].copy()
    r1["ambiguity"] = r1["runner_up_region_ceiling"] / r1["region_ceiling"]
    has2 = set(allsites.loc[allsites["rank"] == 2, "member_id"])
    second = (
        r1[r1["member_id"].isin(has2)]
        .sort_values("ambiguity", ascending=False)
        .head(TOP_N_SECOND)["member_id"]
    )
    q = pd.concat(
        [
            allsites[allsites["rank"] == 1],
            allsites[(allsites["rank"] == 2) & allsites["member_id"].isin(set(second))],
        ]
    )
    prot = ordered_proteins()
    prot.update(current_ordered_proteins())  # the adopted orthologs postdate the cached copy
    tagged = []
    for r in q.itertuples(index=False):
        p = prot[r.sequence_sha256]
        n = int(r.insert_after_residue)
        tagged.append(hashlib.sha256((p[:n] + CASSETTE_AA + p[n:]).encode()).hexdigest())
    q = q.assign(tagged_protein_sha256=tagged).rename(columns={"member_id": "candidate_id"})
    q = q[
        [
            "candidate_id",
            "gene",
            "rank",
            "insert_after_residue",
            "sequence_sha256",
            "tagged_protein_sha256",
            "snorkel_candidate",
            "source",
        ]
    ].sort_values(["rank", "candidate_id"])
    q.to_csv(OUT, sep="\t", index=False)
    print(
        f"{q['candidate_id'].nunique()} members, {len(q)} constructs "
        f"({int((q['rank'] == 1).sum())} rank 1, {int((q['rank'] == 2).sum())} rank 2); "
        f"orthologs from {orth_commit}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
