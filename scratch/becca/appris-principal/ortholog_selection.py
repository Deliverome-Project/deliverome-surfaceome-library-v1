"""Re-pick the mouse ortholog proteoform, and compare three ways of doing it.

The one-to-one high-confidence Compara mapping is **gene level** and is not in question.
What is in question is which proteoform was taken inside each confirmed ortholog gene, and
the answer turns out to be the weakest available rule: the current cohort took the
**longest symbol-matched entry, reviewed or not**. 654 of 1,971 mouse picks are unreviewed
TrEMBL, and for 543 of those a Swiss-Prot entry existed and was passed over; against the
longest reviewed proteoform the selection ran longer in 265 cases and shorter in 13.

So this is the same length tiebreak that made UniProt's human canonical unreliable, minus
the reviewed-entry floor. Three alternatives are scored here:

* **reviewed** — prefer the Swiss-Prot canonical for the gene.
* **appris** — the mouse APPRIS principal, selected per species.
* **correspondence** — the candidate whose sequence best matches the **human reference we
  selected**, by alignment. An ortholog construct exists to be the counterpart of a human
  one, so resembling that human form is the property that actually matters, and selecting
  on it makes "human P1 and mouse P1 are different splice forms" impossible by construction.

Scored on identity to the human reference and on length ratio, since no ground truth exists.

    uv run --locked python ortholog_selection.py
"""

from __future__ import annotations

import argparse
import csv
import json
import statistics
from collections import Counter, defaultdict
from difflib import SequenceMatcher
from pathlib import Path

from select_principal import (
    APPRIS_RELEASE,
    load_uniprot_proteoforms,
    parse_appris,
    parse_ensembl_pep,
    parse_uniprot_fasta,
    select_principal,
)


def identity(a: str, b: str) -> float:
    """Fraction of the longer sequence covered by matching blocks."""
    sm = SequenceMatcher(None, a, b, autojunk=False)
    matched = sum(bl.size for bl in sm.get_matching_blocks())
    return matched / max(len(a), len(b))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--surfaceome-repo", default="~/Git/accessible-surfaceome")
    ap.add_argument("--cache-dir", default="~/.cache/deliverome/appris")
    ap.add_argument("--universe", default="isoform_universe_scored.tsv")
    args = ap.parse_args()

    repo = Path(args.surfaceome_repo).expanduser()
    cache = Path(args.cache_dir).expanduser()
    seqdir = repo / "data/external/deeptmhmm_surfaceome_inputs"

    human_ref = {
        r["gene"]: r for r in csv.DictReader(open(args.universe), delimiter="\t")
        if r["is_reference"] == "True"
    }
    # The selected human reference sequence, from the cohort or the UniProt stream.
    human_up = load_uniprot_proteoforms(cache / "uniprot_human_reviewed_isoforms.tsv.gz")
    human_seq: dict[str, str] = {}
    for gene, r in human_ref.items():
        for p in human_up.get(gene, []):
            if p["accession"] == r["uniprot"]:
                human_seq[gene] = p["seq"]
                break

    mouse_cohort = parse_uniprot_fasta(
        seqdir / "mouse_ortholog_one2one_highconf_non_hla.fasta"
    )
    mouse_up = load_uniprot_proteoforms(cache / "uniprot_mouse_reviewed_isoforms.tsv.gz")
    meta = [
        m for m in csv.DictReader(
            (seqdir / "mouse_ortholog_one2one_highconf_non_hla_metadata.csv").open()
        )
        if m.get("status") == "ok"
    ]
    appris_m = parse_appris(cache / f"appris_mus_musculus_{APPRIS_RELEASE}.txt")
    pep_m = parse_ensembl_pep(cache / "Mus_musculus.GRCm39.pep.all.fa.gz")
    appris_by_symbol: dict[str, list[dict]] = defaultdict(list)
    for rows in appris_m.values():
        for r in rows:
            appris_by_symbol[r["symbol"]].append(r)

    rows, skipped = [], Counter()
    for m in meta:
        human_symbol = (m["query_input_gene_symbols"] or "").split(";")[0].strip()
        mouse_symbol = m["target_gene_symbol"]
        ref = human_seq.get(human_symbol)
        if not ref:
            skipped["no human reference sequence"] += 1
            continue
        current_acc = m["selected_uniprot_accession"]
        current = mouse_cohort.get(current_acc)
        if not current:
            skipped["current pick not in the cohort FASTA"] += 1
            continue

        candidates: dict[str, str] = {current_acc: current["seq"]}
        for p in mouse_up.get(mouse_symbol, []):
            candidates[p["accession"]] = p["seq"]
        sel = select_principal(appris_by_symbol.get(mouse_symbol, []), pep_m)
        appris_seq = sel["seq"] if sel and sel["seq"] else None

        scored = {a: identity(ref, s) for a, s in candidates.items()}
        best_acc = max(scored, key=lambda a: scored[a])
        reviewed = [p for p in mouse_up.get(mouse_symbol, []) if not p["is_isoform"]]
        reviewed_acc = reviewed[0]["accession"] if reviewed else ""

        def ratio(s: str) -> float:
            return min(len(s), len(ref)) / max(len(s), len(ref))

        rows.append(
            {
                "human_symbol": human_symbol,
                "mouse_symbol": mouse_symbol,
                "human_ref_length": len(ref),
                "current_acc": current_acc,
                "current_reviewed": m["selected_reviewed"],
                "current_identity": round(scored[current_acc], 4),
                "current_length_ratio": round(ratio(current["seq"]), 4),
                "reviewed_acc": reviewed_acc,
                "reviewed_identity": round(scored.get(reviewed_acc, 0.0), 4)
                if reviewed_acc else "",
                "appris_identity": round(identity(ref, appris_seq), 4) if appris_seq else "",
                "correspondence_acc": best_acc,
                "correspondence_identity": round(scored[best_acc], 4),
                "n_candidates": len(candidates),
                "correspondence_changes_pick": best_acc != current_acc,
            }
        )

    with open("ortholog_selection_mouse.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote ortholog_selection_mouse.tsv ({len(rows)} orthologs)")
    for k, n in skipped.most_common():
        print(f"  skipped — {k}: {n}")

    def med(key):
        vals = [r[key] for r in rows if r[key] != ""]
        return round(statistics.median(vals), 4) if vals else None

    summary = {
        "orthologs_scored": len(rows),
        "median_identity_to_human_reference": {
            "current (longest symbol match)": med("current_identity"),
            "reviewed (Swiss-Prot canonical)": med("reviewed_identity"),
            "appris (mouse principal)": med("appris_identity"),
            "correspondence (best match)": med("correspondence_identity"),
        },
        "correspondence_changes_the_pick": sum(
            1 for r in rows if r["correspondence_changes_pick"]
        ),
        "current_pick_unreviewed": sum(
            1 for r in rows if r["current_reviewed"] == "unreviewed"
        ),
        "correspondence_beats_current_by_over_5pct": sum(
            1 for r in rows
            if r["correspondence_identity"] - r["current_identity"] > 0.05
        ),
    }
    Path("ortholog_selection_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
