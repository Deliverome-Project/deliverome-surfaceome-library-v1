"""Pick the library receptors that get an exogenous Igk leader, and which also get recoded.

Tedman's screen measures native surface immunostaining for 763 canonical GPCRs; 526 are
library members and 276 sit below his own "expressed" bar, the 67th percentile at 1,025 a.u.
Two interventions are on the table for the worst of them, and this picks who receives each.

**Every class is represented, and the draw is balanced across them rather than ranked
straight.** Olfactory receptors dominate the bottom of the screen -- 17 of the 30 lowest --
so a straight ranking returns an olfactory experiment and says nothing about the rest. An
equal-as-possible split across classes, capped by what each has available and filled from
the lowest expressor upward within each, keeps the olfactory arm informative without letting
it crowd out the others.

Olfactory receptors are the class where a leader is least likely to be the fix: in the Lv
2016 insect-cell screen they were rescued by an N-terminal BRIL fusion rather than by a
signal peptide, so the fusion partner, not the leader, is what moved them. That makes them
worth testing -- a leader-only arm is a real test of whether the leader alone does anything
for the hardest class -- but not worth 17 of 25 slots.

**A receptor that already has a signal peptide does not get a second one**, and one whose
N-terminus is cytoplasmic must never get one -- the leader would invert the topology and
throw the cytoplasmic sorting determinants outside. Both are checked against DeepTMHMM
rather than assumed from the class.

The published support for the leader is CB1: adding a signal peptide greatly enhances its
stability and surface expression, because its 116-residue N-tail is poorly translocated
without one (Andersson et al. 2003, PMID 12920192).

Recoding is a second, independent intervention, so the subset that receives both is spread
across the intensity range rather than taken from the extreme tail: confining it to the
worst receptors would confound "does recoding add anything" with a floor effect, and leave
no receptor where the leader alone was nearly enough.

    uv run --locked python select_expression_rescue_arms.py
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from collections import defaultdict

import pandas as pd

from deliverome_analysis.core.s3 import external
from library import library_members
from sample_isoform_panel import _topology_sources

SCREEN = "tedman-gpcr-surface-screen/2025-09-19/tedman2025_gpcr_screen_media-2.xlsx"
EXPRESSED_PERCENTILE = 0.67      # Tedman's own "expressed" bar
EXCLUDE_CLASSES: set[str] = set()
N_LEADER = 25
N_RECODE = 10


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--leader", type=int, default=N_LEADER)
    ap.add_argument("--recode", type=int, default=N_RECODE)
    ap.add_argument("--out", default="expression_rescue_arms.tsv")
    ap.add_argument("--surfaceome-repo", default="~/Git/accessible-surfaceome",
                    help="checkout (or staged snapshot) holding data/processed/deeptmhmm")
    ap.add_argument("--library", default=None, help="gene-list JSON (default: the frozen root)")
    ap.add_argument("--reference", default="reference_selection.tsv")
    args = ap.parse_args()

    ted = pd.read_excel(external.download(SCREEN), sheet_name="Canonical")[
        ["gene_id", "GPCR Class", "Immunostaining Intensity"]
    ].dropna(subset=["gene_id"])
    ted.columns = ["symbol", "gpcr_class", "intensity"]
    bar = float(ted.intensity.quantile(EXPRESSED_PERCENTILE))

    lib = library_members(args.library) if args.library else library_members()
    by_acc, by_ref_gene = _topology_sources(repo=args.surfaceome_repo)
    selection = {
        r["symbol"]: r for r in csv.DictReader(open(args.reference), delimiter="\t")
    }

    def topology(gene: str) -> dict | None:
        row = selection.get(gene)
        if row is None:
            return None
        acc = row.get("selected_uniprot_equivalent") or row.get("current_accession") or ""
        return by_acc.get(acc) or by_acc.get(acc.split("-")[0]) or by_ref_gene.get(gene)

    rows, skipped = [], {"not_in_library": 0, "above_bar": 0, "excluded_class": 0,
                         "has_signal_peptide": 0, "n_terminus_cytoplasmic": 0,
                         "no_topology": 0}
    for r in ted.sort_values("intensity").itertuples():
        if r.symbol not in lib:
            skipped["not_in_library"] += 1
            continue
        if r.intensity >= bar:
            skipped["above_bar"] += 1
            continue
        if r.gpcr_class in EXCLUDE_CLASSES:
            skipped["excluded_class"] += 1
            continue
        topo = topology(r.symbol)
        if topo is None:
            skipped["no_topology"] += 1
            continue
        if topo["sp"] > 0:
            skipped["has_signal_peptide"] += 1
            continue
        if topo["n_term"] != "extracellular":
            skipped["n_terminus_cytoplasmic"] += 1
            continue
        rows.append({"symbol": r.symbol, "gpcr_class": r.gpcr_class,
                     "intensity": round(r.intensity, 1), "tm_helices": topo["tm"]})

    eligible = rows                        # already ascending by intensity
    # Equal-as-possible across classes, capped by availability; the shortfall from thin
    # classes is handed back to the ones that still have candidates, largest pool first.
    by_class: dict[str, list[dict]] = defaultdict(list)
    for row in eligible:
        by_class[row["gpcr_class"]].append(row)
    quota = {c: 0 for c in by_class}
    left = args.leader
    while left > 0:
        open_classes = [c for c in by_class if quota[c] < len(by_class[c])]
        if not open_classes:
            break
        share = max(left // len(open_classes), 1)
        for c in sorted(open_classes, key=lambda c: -len(by_class[c])):
            if left <= 0:
                break
            take = min(share, len(by_class[c]) - quota[c], left)
            quota[c] += take
            left -= take
    leader = sorted(
        (row for c, n in quota.items() for row in by_class[c][:n]),
        key=lambda r: r["intensity"],
    )
    # Spread the recode subset across the leader set's intensity range.
    step = max(len(leader) // args.recode, 1)
    recode = {leader[i]["symbol"] for i in list(range(0, len(leader), step))[: args.recode]}

    for row in leader:
        row["arm"] = "igk_leader + recoded" if row["symbol"] in recode else "igk_leader"
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(leader[0]), delimiter="\t")
        w.writeheader()
        w.writerows(leader)

    summary = {
        "tedman_expressed_bar_au": round(bar, 1),
        "in_library_and_below_bar": sum(1 for r in ted.itertuples()
                                        if r.symbol in lib and r.intensity < bar),
        "eligible_after_filters": len(eligible),
        "leader_arms": len(leader),
        "also_recoded": len(recode),
        "skipped": skipped,
        "classes": pd.Series([r["gpcr_class"] for r in leader]).value_counts().to_dict(),
    }
    Path("expression_rescue_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"wrote {args.out} ({len(leader)} receptors, {len(recode)} of them recoded)")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
