"""Build the expression-rescue ORFs: Igk leader fused to each receptor.

Each construct is the murine Ig kappa-chain signal peptide, the native +1 Asp
that every published Igk vector keeps, and then the receptor's own CDS with its
start codon intact, so the mature protein is the native receptor preceded by a
single Asp.

Domestication runs on the **assembled** fusion rather than on the receptor
alone, because a site can be created at the seam where two clean sequences meet
-- the same reason barcodes are screened against their flanking primer sites
rather than in isolation.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from deliverome_analysis import seqtools

ENZYMES = ["BbsI", "Esp3I", "BsaI"]
LEADER = "ATGGAGACAGACACACTCCTGCTATGGGTACTGCTGCTCTGGGTTCCAGGTTCCACTGGT"
PLUS1 = "GAC"  # native first residue of the mature kappa chain; see igk_leader.tsv
SEED = 20261011


def seam_sites(seq: str, a: int, b: int) -> list[tuple[str, str, int]]:
    """Sites that straddle the cut between a and b, so neither part owns them.

    ``find_sites`` reports 0-based starts, so a six-base site straddling the
    boundary starts somewhere in the six bases before it.
    """
    return [s for s in seqtools.find_sites(seq, ENZYMES) if a - 7 < s[2] < b + 1]


def check_seam_detector() -> None:
    """A seam check that reports nothing is only worth believing if it can fire.

    Builds one fusion with a site deliberately straddling the junction and one
    with a site far from it, and asserts the detector separates them.
    """
    boundary = len(LEADER + PLUS1)
    tail = "GCCGCCGCC" * 40
    straddling = LEADER + "GGT" + "CTC" + tail      # GGTCTC spans the boundary
    distant = LEADER + PLUS1 + tail[:300] + "GGTCTC" + tail[300:]
    if len(seam_sites(straddling, boundary, boundary)) != 1:
        raise AssertionError("seam detector missed a site straddling the junction")
    if seam_sites(distant, boundary, boundary):
        raise AssertionError("seam detector counted a site far from the junction")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--arms", default="expression_rescue_arms.tsv")
    ap.add_argument("--cds", default="rescue_cds.json")
    ap.add_argument("--out", default="rescue_constructs.tsv")
    args = ap.parse_args()

    check_seam_detector()
    store = json.loads(Path(args.cds).read_text())
    with open(args.arms) as fh:
        arms = list(csv.DictReader(fh, delimiter="\t"))

    rows, missing = [], []
    for arm in arms:
        sym = arm["symbol"]
        rec = store.get(sym)
        if rec is None:
            missing.append(sym)
            continue

        fusion = LEADER + PLUS1 + rec["cds"]
        before = seqtools.find_sites(fusion, ENZYMES)
        seam_before = seam_sites(fusion, len(LEADER + PLUS1), len(LEADER + PLUS1))

        tamed, warnings = seqtools.domesticate(fusion, ENZYMES, seed=SEED)
        after = seqtools.find_sites(tamed, ENZYMES)
        edits = sum(1 for x, y in zip(fusion, tamed) if x != y)
        leader_touched = tamed[: len(LEADER)] != LEADER

        # the protein must be untouched, leader + Asp + native receptor
        want = seqtools.translate(LEADER + PLUS1 + rec["cds"])
        got = seqtools.translate(tamed)
        assert want == got, f"{sym}: domestication changed the protein"

        rows.append(
            {
                "symbol": sym,
                "gpcr_class": arm["gpcr_class"],
                "arm": arm["arm"],
                "transcript": rec["transcript_id"],
                "receptor_nt": len(rec["cds"]),
                "construct_nt": len(tamed),
                "sites_before": len(before),
                "sites_at_seam": len(seam_before),
                "sites_after": len(after),
                "edits": edits,
                "leader_modified": leader_touched,
                "warnings": "; ".join(warnings) or "",
                "construct": tamed,
            }
        )
        seam_note = f"  SEAM x{len(seam_before)}" if seam_before else ""
        lead_note = "  LEADER EDITED" if leader_touched else ""
        print(
            f"  {sym:9} {len(rec['cds']):5}nt  sites {len(before)}->{len(after)}"
            f"  edits {edits:3}{seam_note}{lead_note}"
        )

    if missing:
        print(f"\nno CDS for: {', '.join(missing)}")

    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    clean = sum(1 for r in rows if r["sites_after"] == 0)
    print(f"\nwrote {args.out} ({len(rows)} constructs, {clean} clear all three enzymes)")
    print(f"seam-created sites found in {sum(1 for r in rows if r['sites_at_seam'])} constructs")
    print(f"leader modified in {sum(1 for r in rows if r['leader_modified'])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
