"""Merge per-batch audit verdicts into ``audit.tsv`` and print the tallies for the README.

Inputs: ``sites.tsv`` from ``fetch_sites.py`` (the 191 live rows, 176 flagged
``in_training``) and one or more verdict TSVs produced by hand/agent reading
(``--verdicts``); ``--overrides`` are later readings (a cited construct paper, a
supplement) whose non-empty fields replace the earlier verdict for the same site. Each verdict row must name a ``site_id``; the merge keeps the
API fields from ``sites.tsv`` and the audit fields from the verdict files, and
fails loudly if a training site has no verdict or has two.

Run::

    uv run --locked python library/l1/tag-site-corpus-audit/merge_audit.py \
        --sites library/l1/tag-site-corpus-audit/data/sites.tsv \
        --verdicts /path/batch*.tsv --out library/l1/tag-site-corpus-audit/audit.tsv
"""

from __future__ import annotations

import argparse
import collections
import csv
from pathlib import Path

AUDIT_COLS = ["pmid_resolves_to_cited_paper", "paper_quote", "paper_numbering_basis",
              "species_of_construct", "verdict", "corrected_residue", "surface_shown",
              "function_shown", "notes"]
VERDICTS = ["confirmed_exact", "confirmed_shifted", "region_only", "wrong_species_numbering",
            "not_supported", "terminal_misfiled", "paper_inaccessible"]
OUT_COLS = ["site_id", "gene", "uniprot_acc", "in_training", "insert_after_residue", "tag_type",
            "evidence_type", "position_evidence", "validation_level", "pmid", *AUDIT_COLS]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sites", required=True, type=Path)
    ap.add_argument("--verdicts", required=True, nargs="+", type=Path)
    ap.add_argument("--overrides", nargs="*", default=[], type=Path,
                    help="verdict TSVs from a later reading that replace earlier verdicts for the same site_id")
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    with args.sites.open() as fh:
        sites = {r["site_id"]: r for r in csv.DictReader(fh, delimiter="\t")}
    verdicts: dict[str, dict] = {}
    for p in args.verdicts:
        with p.open() as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                sid = r["site_id"].strip()
                if sid in verdicts:
                    raise SystemExit(f"duplicate verdict for {sid} in {p}")
                if r["verdict"] not in VERDICTS:
                    raise SystemExit(f"bad verdict {r['verdict']!r} for {sid} in {p}")
                verdicts[sid] = r

    for p in args.overrides:
        with p.open() as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                sid = r["site_id"].strip()
                if sid not in verdicts:
                    raise SystemExit(f"override for {sid} in {p} but no earlier verdict")
                if r["verdict"] not in VERDICTS:
                    raise SystemExit(f"bad verdict {r['verdict']!r} for {sid} in {p}")
                verdicts[sid] = {**verdicts[sid], **{k: v for k, v in r.items() if v != ""}}

    training = [s for s in sites.values() if s["in_training"] == "1"]
    missing = [s["site_id"] for s in training if s["site_id"] not in verdicts]
    if missing:
        raise SystemExit(f"{len(missing)} training sites without a verdict: {missing}")
    extra = sorted(set(verdicts) - set(sites))
    if extra:
        raise SystemExit(f"verdicts for unknown sites: {extra}")

    out_rows = []
    for s in sorted(sites.values(), key=lambda x: (x["gene"], int(x["insert_after_residue"] or -1))):
        v = verdicts.get(s["site_id"])
        if v is None:
            continue  # live rows outside the training set were not audited
        row = {c: s.get(c, "") for c in OUT_COLS if c not in AUDIT_COLS}
        row.update({c: (v.get(c) or "").replace("\t", " ").replace("\n", " ").strip() for c in AUDIT_COLS})
        out_rows.append(row)
    with args.out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=OUT_COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out_rows)

    tr = [r for r in out_rows if r["in_training"] == "1"]
    print(f"{len(tr)} training sites audited\n")
    print("verdict\tn")
    for k, n in collections.Counter(r["verdict"] for r in tr).most_common():
        print(f"{k}\t{n}")
    print("\nevidence_type x verdict")
    tab = collections.Counter((r["evidence_type"], r["verdict"]) for r in tr)
    for et in sorted({r["evidence_type"] for r in tr}):
        print(et, {v: tab[(et, v)] for v in VERDICTS if tab[(et, v)]})
    print("\nposition_evidence x verdict")
    tab = collections.Counter((r["position_evidence"], r["verdict"]) for r in tr)
    for et in sorted({r["position_evidence"] for r in tr}):
        print(et, {v: tab[(et, v)] for v in VERDICTS if tab[(et, v)]})
    print("\nvalidation_level x verdict")
    tab = collections.Counter((r["validation_level"], r["verdict"]) for r in tr)
    for et in sorted({r["validation_level"] for r in tr}):
        print(et, {v: tab[(et, v)] for v in VERDICTS if tab[(et, v)]})
    print("\nsurface_shown", collections.Counter(r["surface_shown"] for r in tr))
    print("function_shown", collections.Counter(r["function_shown"] for r in tr))
    print("pmid_resolves", collections.Counter(r["pmid_resolves_to_cited_paper"] for r in tr))
    print("species", collections.Counter(r["species_of_construct"] for r in tr))
    claimed_fn = [r for r in tr if r["validation_level"] in ("surface_and_function", "function_only", "function_perturbed")]
    print("function claimed but not shown:", sum(1 for r in claimed_fn if r["function_shown"] != "yes"), "of", len(claimed_fn))
    claimed_sf = [r for r in tr if r["validation_level"] in ("surface_and_function", "surface_only")]
    print("surface claimed but not shown:", sum(1 for r in claimed_sf if r["surface_shown"] != "yes"), "of", len(claimed_sf))


if __name__ == "__main__":
    main()
