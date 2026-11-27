#!/usr/bin/env python3
"""Proposed internal tag sites for the L1 internal rung, as tagged protein sequences for refolding.

Takes the v5 masked design of the adopted L1 panel (`tag_sites_v5_adopted.tsv`), and for each
member builds the ordered protein (designed CDS translated, matched by `sequence_sha256`) with
the full 51-aa cassette inserted after the proposed residue:

    GGGGS | SpyTag003 | GGGGSGGGGS | P-ALFA v1.5-P | GGGGS
    = s4e9 build_constructs.CASSETTE_AA, checked against the adopted panel's copy at runtime

Picks sent for refolding: rank 1 for every designed member, and rank 2 where the design
recommends two or more sites. Every pick is a PROPOSAL for Becca, not a final choice.

    l1_internal_proposed_sites.tsv     one row per proposed site, with member_id, design_id,
                                       scores, snorkel flag, refold priority and attribution
    l1_internal_tagged_constructs.fasta  the tagged proteins, header `member_id|design_id|rank|after`

Snorkel: `snorkel_candidate` marks a member whose best membrane-side-eligible loop scores
below the 25th percentile of loops hosting a published site. For those, a snorkel construct
(a transmembrane domain plus the tag appended to the cytoplasmic C-terminus; Brown et al.
2013, PMID 24023844) is the alternative to weigh; the internal pick is still proposed.

    uv run --locked python build_internal_tagged.py
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[3] / "src"))

import membrane_mask  # noqa: E402
from map_contacts_to_members import ordered_proteins  # noqa: E402
from tag_site_model_v5 import ADOPTED_COMMIT  # noqa: E402

from deliverome_analysis.cloning import plasmid_qc as qc  # noqa: E402

CASSETTE_AA = qc.P003_CASSETTE_PEPTIDE
DESIGNS = HERE / "tag_sites_v5_adopted.tsv"
OUT = HERE / "l1_internal_proposed_sites.tsv"
FASTA = HERE / "l1_internal_tagged_constructs.fasta"
DECIDED = (
    "Becca: v5 model and mask (2026-11-11), loop floor 4 (2026-11-11), snorkel flag "
    "(2026-11-11), SYGQ exclusion (2026-11-11); picks are proposals"
)
IMPLEMENTED = "ranker-pu-segment session (Claude Code), PR #245"


def check_cassette() -> None:
    """The cassette must be the s4e9 constant; fail loudly if the adopted panel disagrees."""
    try:
        src = subprocess.check_output(
            [
                "git",
                "show",
                f"{ADOPTED_COMMIT}:library/l1/orf-design/build_constructs.py",
            ],
            cwd=HERE.parents[2],
            text=True,
        )
    except subprocess.CalledProcessError:
        print("note: could not read s4e9 build_constructs.py to cross-check the cassette")
        return
    assert f'CASSETTE_AA = "{CASSETTE_AA}"' in src, "cassette differs from s4e9 CASSETTE_AA"


def main() -> int:
    assert len(CASSETTE_AA) == 51
    check_cassette()
    d = pd.read_csv(DESIGNS, sep="\t")
    t = pd.read_csv(HERE / "l1_adopted_targets.tsv", sep="\t")
    sha = dict(zip(t["uniprot_acc"], t["sequence_sha256"], strict=False))
    mask = membrane_mask.load()
    ids = mask.drop_duplicates("sequence_sha256").set_index("sequence_sha256")[
        ["candidate_id", "design_id"]
    ]
    prot = ordered_proteins()
    d["sequence_sha256"] = d["uniprot_acc"].map(sha)
    d["member_id"] = d["sequence_sha256"].map(ids["candidate_id"])
    d["design_id"] = d["sequence_sha256"].map(ids["design_id"])
    assert d["sequence_sha256"].notna().all() and d["member_id"].notna().all()
    send = (d["rank"] == 1) | ((d["rank"] == 2) & (d["recommended_sites_to_order"] >= 2))
    d["refold_priority"] = d["rank"].where(send)
    rows, fa = [], []
    for r in d.itertuples(index=False):
        p = prot[r.sequence_sha256]
        n = int(r.insert_after_residue)
        assert 1 <= n < len(p)
        tagged = p[:n] + CASSETTE_AA + p[n:]
        assert len(tagged) == len(p) + 51 and tagged[n : n + 51] == CASSETTE_AA
        rows.append(
            {
                "member_id": r.member_id,
                "design_id": r.design_id,
                "gene": r.gene,
                "frame_accession": r.uniprot_acc,
                "sequence_sha256": r.sequence_sha256,
                "rank": r.rank,
                "insert_after_residue": n,
                "junction": f"{p[n - 1]}{n}|{p[n]}{n + 1}",
                "rank_score": r.rank_score,
                "absolute_score": r.absolute_score,
                "best_region_pct_of_published_regions": r.best_region_pct_of_published_regions,
                "this_region_length": r.this_region_length,
                "recommended_sites_to_order": r.recommended_sites_to_order,
                "why": r.why,
                "snorkel_candidate": r.snorkel_candidate,
                "refold_priority": r.refold_priority,
                "native_length": len(p),
                "tagged_length": len(tagged),
                "status": "proposed",
                "decided_by": DECIDED,
                "implemented_by": IMPLEMENTED,
            }
        )
        if r.refold_priority == r.refold_priority:  # not NaN
            fa.append(f">{r.member_id}|{r.design_id}|rank{r.rank}|after{n}\n{tagged}\n")
    out = pd.DataFrame(rows)
    out.to_csv(OUT, sep="\t", index=False)
    FASTA.write_text("".join(fa))
    print(
        f"{out['member_id'].nunique()} members, {len(out)} proposed sites; "
        f"{len(fa)} tagged constructs for refolding "
        f"({int((out['refold_priority'] == 1).sum())} rank 1, "
        f"{int((out['refold_priority'] == 2).sum())} rank 2); "
        f"snorkel candidates {int(out.loc[out['rank'] == 1, 'snorkel_candidate'].sum())}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
