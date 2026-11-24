"""Deduplicated gold standard of literature tag sites that survive the audit.

Reads ``audit.tsv`` and keeps rows whose verdict is ``confirmed_exact`` or
``confirmed_shifted``; the residue is ``corrected_residue`` where the audit set one,
else the row's ``insert_after_residue`` (human canonical ``insert_after`` numbering).
Five confirmed sites are excluded because the paper shows the tagged construct failed
(``FAILED_SITES``): the position is right, but the tag was not tolerated there, so they are
written to ``failed_sites.tsv`` as candidate negatives rather than positives.
Rows that land on the same (accession, residue) from different papers are merged:
site ids and PMIDs are ``;``-joined, the verdict is ``confirmed_exact`` if any row is,
and ``surface_shown`` / ``function_shown`` take the best value across rows
(``yes`` > ``not_reported`` > ``no``).

Run::

    uv run --locked python library/l1/tag-site-corpus-audit/gold_standard.py \
        --audit library/l1/tag-site-corpus-audit/audit.tsv \
        --out library/l1/tag-site-corpus-audit/gold_standard.tsv
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

RANK = {"yes": 0, "not_reported": 1, "no": 2, "": 3}
# Confirmed positions where the paper reports the construct did not reach the surface or
# could not be labeled; all five are substitutions, not insertions.
FAILED_SITES = {
    "SLC6A3-internal-288-lit": "tetraCys EL3 construct failed to mature and traffic to the cell surface",
    "SLC6A3-internal-387-lit": "tetraCys EL4 construct failed to mature and traffic to the cell surface",
    "NFASC-internal-533-lit": "amber mutant excluded by the authors: poor expression, weak/absent labeling",
    "NFASC-internal-570-lit": "amber mutant excluded by the authors: poor expression, weak/absent labeling",
    "TRPV1-internal-464-lit": "ncAA mutant was not efficiently labelled",
}
COLS = ["uniprot_acc", "gene", "insert_after_residue", "verdict", "surface_shown",
        "function_shown", "pmid", "site_ids", "in_training", "n_rows"]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--audit", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--added", type=Path, default=None,
                    help="added_sites.tsv: real sites recovered from papers the agent had filed "
                         "under the wrong gene (not live API rows); same verdict/surface/function columns")
    args = ap.parse_args()

    groups: dict[tuple[str, int], list[dict]] = {}
    failed = []
    with args.audit.open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["verdict"] not in ("confirmed_exact", "confirmed_shifted"):
                continue
            residue = int(r["corrected_residue"] or r["insert_after_residue"])
            if r["site_id"] in FAILED_SITES:
                failed.append({"uniprot_acc": r["uniprot_acc"], "gene": r["gene"],
                               "insert_after_residue": residue, "site_id": r["site_id"],
                               "pmid": r["pmid"], "why_failed": FAILED_SITES[r["site_id"]]})
                continue
            groups.setdefault((r["uniprot_acc"], residue), []).append(r)
    assert len(failed) == len(FAILED_SITES), "a FAILED_SITES id is missing or no longer confirmed"
    if args.added is not None:
        with args.added.open() as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r["verdict"] not in ("confirmed_exact", "confirmed_shifted"):
                    continue
                r = {**r, "site_id": r.get("site_id") or f"{r['gene']}-added-{r['insert_after_residue']}",
                     "in_training": "0", "corrected_residue": ""}
                groups.setdefault((r["uniprot_acc"], int(r["insert_after_residue"])), []).append(r)
    with args.out.with_name("failed_sites.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(failed[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(failed)

    out = []
    for (acc, residue), rs in sorted(groups.items(), key=lambda kv: (kv[1][0]["gene"], kv[0][1])):
        out.append({
            "uniprot_acc": acc, "gene": rs[0]["gene"], "insert_after_residue": residue,
            "verdict": "confirmed_exact" if any(r["verdict"] == "confirmed_exact" for r in rs)
            else "confirmed_shifted",
            "surface_shown": min((r["surface_shown"] for r in rs), key=lambda v: RANK.get(v, 3)),
            "function_shown": min((r["function_shown"] for r in rs), key=lambda v: RANK.get(v, 3)),
            "pmid": ";".join(dict.fromkeys(r["pmid"] for r in rs if r["pmid"])),
            "site_ids": ";".join(r["site_id"] for r in rs),
            "in_training": int(any(r["in_training"] == "1" for r in rs)),
            "n_rows": len(rs),
        })
    with args.out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out)
    n_prot = len({r["uniprot_acc"] for r in out})
    n_surf = sum(1 for r in out if r["surface_shown"] == "yes")
    n_surf_prot = len({r["uniprot_acc"] for r in out if r["surface_shown"] == "yes"})
    print(f"{sum(len(v) for v in groups.values())} audit rows -> {len(out)} sites, {n_prot} proteins; "
          f"surface shown {n_surf} sites / {n_surf_prot} proteins; "
          f"surface negative {sum(1 for r in out if r['surface_shown'] == 'no')}; "
          f"function shown {sum(1 for r in out if r['function_shown'] == 'yes')}; "
          f"merged sites {sum(1 for r in out if r['n_rows'] > 1)}; "
          f"failed sites excluded {len(failed)}")


if __name__ == "__main__":
    main()
