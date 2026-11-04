#!/usr/bin/env python3
"""Join deterministic TPM (GTEx/HPA-tissue/HPA-celline/TCGA) to cohort + partner
genes, canonicalizing symbols with the repo gene resolver (MyGene) and joining on
Ensembl gene id (alias-proof). See external/expression_tpm/PROVENANCE.md.

Outputs:
  (default) external/expression_tpm/resolved_expr.json -- cohort + co-receptor
      partners, {symbol: {ensembl, max_expression_any, max_tpm_gtex,
      max_ntpm_hpa_tissue, max_tpm_hpa_celline, max_tpm_tcga}}, joined against
      genomewide_expression_tpm.csv (PRE-FILTERED to >=10 TPM).
  (--full)  external/expression_tpm/resolved_expr_full.json -- cohort only,
      flat {symbol: max_expression_any}, joined against the UNFILTERED
      genomewide_expression_tpm_rederived.csv.

Why both: the >=10 table has no row for a gene below 10 TPM, so the default
join records null there. That reads correctly as 'fails' for a >=10 gate, but
is unusable once the gate drops to >=5 (Willow validation, 2026-08-14) -- null
then means 'below 10', not 'no data'. --full supplies the real values.

Run: uv run python build_expression_join.py [--full]
"""
# ruff: noqa: E741  (exploratory scratch script)

import argparse
import csv
import json
import os
import re

from deliverome_analysis.gene_resolver import resolve

HERE = os.path.dirname(os.path.abspath(__file__))
TPM_CSV = os.path.join(HERE, "external/expression_tpm/genomewide_expression_tpm.csv")
CACHE = os.path.join(HERE, "external/mygene_symbol_resolution")
# The two passes resolve different symbol sets, and the resolver keeps one response per
# cache directory, so they each get their own: a pinned pair then reruns offline.
CACHE_FULL = os.path.join(CACHE, "full")
OUT = os.path.join(HERE, "external/expression_tpm/resolved_expr.json")
# Unfiltered re-derivation (65,098 rows, no >=10 cut) + its cohort-only flat map.
FULL_TPM_CSV = os.path.join(HERE, "external/expression_tpm/genomewide_expression_tpm_rederived.csv")
FULL_OUT = os.path.join(HERE, "external/expression_tpm/resolved_expr_full.json")
COHORT = os.path.join(HERE, "cohort.jsonl")
NUM = (
    "max_tpm_gtex",
    "max_ntpm_hpa_tissue",
    "max_tpm_hpa_celline",
    "max_tpm_tcga",
    "max_expression_any",
)


def _index_tpm(path):
    """symbol -> vals and ensembl-base -> vals, for one TPM table.

    The two tables disagree on the ensembl column name: the >=10 export uses
    ``ensembl_bases`` (plural, possibly multi-valued), the full re-derivation
    uses ``ensembl_base``. Accept either.
    """
    tpm_by_ens, tpm_by_sym = {}, {}
    for row in csv.DictReader(open(path)):
        vals = {}
        for c in NUM:
            try:
                vals[c] = float(row[c])
            except TypeError, ValueError:
                vals[c] = None
        sym = (row.get("gene_symbol_canonical") or "").strip()
        if sym:
            tpm_by_sym[sym] = vals
        raw_ens = row.get("ensembl_bases") or row.get("ensembl_base") or ""
        for e in str(raw_ens).replace(";", ",").replace("|", ",").split(","):
            e = e.strip().split(".")[0]
            if e:
                tpm_by_ens[e] = vals
    return tpm_by_ens, tpm_by_sym


def _lookup(sym, resolved, tpm_by_ens, tpm_by_sym):
    r = resolved.get(sym, {})
    for e in r.get("ensembl", []):
        if e in tpm_by_ens:
            return r, tpm_by_ens[e]
    return r, (tpm_by_sym.get(r.get("symbol")) or tpm_by_sym.get(sym))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--full",
        action="store_true",
        help=(
            "Write resolved_expr_full.json instead: cohort-only, flat "
            "{symbol: max_expression_any}, joined against the UNFILTERED "
            "re-derived table. The default table is pre-filtered to >=10 TPM, so "
            "a gene below 10 has no row there and lands as null -- fine for a >=10 "
            "gate, useless for the >=5 gate the Willow validation moved to. This "
            "mode is what gives those genes real values."
        ),
    )
    args = ap.parse_args()

    cohort_records = [json.loads(l) for l in open(COHORT)]
    cohort = [r["gene"] for r in cohort_records]

    if args.full:
        # Cohort only: the full map exists to let the expression gate drop below
        # 10 TPM for cohort genes; partner ubiquity is a separate question and is
        # answered from the HPA cell-line median table, not from here.
        symbols = sorted(set(cohort))
        tpm_by_ens, tpm_by_sym = _index_tpm(FULL_TPM_CSV)
        resolved = resolve(symbols, cache_dir=CACHE_FULL)
        out = {}
        for sym in symbols:
            _r, vals = _lookup(sym, resolved, tpm_by_ens, tpm_by_sym)
            out[sym] = (vals or {}).get("max_expression_any")
        json.dump(out, open(FULL_OUT, "w"))
        cov = sum(1 for v in out.values() if v is not None)
        below10 = sum(1 for v in out.values() if v is not None and v < 10)
        print(
            f"[--full] {len(out)} cohort symbols -> {os.path.basename(FULL_OUT)}; "
            f"TPM coverage {cov}/{len(out)}; {below10} below 10 TPM "
            f"(these are the ones the >=10 table cannot represent)"
        )
        return

    # Co-receptor partner symbols, whose TPM the ubiquity rescue needs.
    #
    # coreceptor_partners.jsonl was the output of a curated pass that read the
    # upstream free-text partner field and resolved it to real symbols. That
    # side-file is gone (not in git, S3, or on disk), but the same pass's outputs
    # ARE in processed/: coreceptor_rationales.json (143 genes, clean `partners`
    # lists) and coreceptor_calls.tsv (146 genes; final_library.py needs it for
    # the gate itself).
    #
    # Do NOT resolve the cohort record's raw `co_receptor_partners` wholesale:
    # it is the un-curated text the pass exists to clean up (1,219 strings, 235
    # of which are not symbols — "14-3-3 proteins", "AMPA receptors (GluA
    # subunits)", "Kir6.2"). Use it only for genes the curated pass never
    # covered, and only for symbol-shaped values, which recovers real partners
    # like GRIN1 and CD147 without importing prose.
    #
    # Verified: this reproduces the committed library to within 6 genes
    # (2340 vs 2346), and those 6 are GRIN2A-D / GRIN3A/B, explained by upstream
    # change rather than by the partner set — see the experiment README.
    SYMBOL_SHAPED = re.compile(r"^[A-Z][A-Z0-9]{1,9}(-[A-Z0-9]{1,4})?$")
    partners: set[str] = set()
    curated_genes: set[str] = set()

    legacy = os.path.join(HERE, "coreceptor_partners.jsonl")
    rationales = os.path.join(HERE, "coreceptor_rationales.json")
    if os.path.exists(legacy):
        for line in open(legacy):
            rec = json.loads(line)
            curated_genes.add(rec.get("gene"))
            partners.update(rec.get("partners") or [])
    elif os.path.exists(rationales):
        for gene, entry in json.load(open(rationales)).items():
            curated_genes.add(gene)
            partners.update(entry.get("partners") or [])
        calls = os.path.join(HERE, "coreceptor_calls.tsv")
        if os.path.exists(calls):
            with open(calls) as fh:
                curated_genes.update(row["gene"] for row in csv.DictReader(fh, delimiter="\t"))
    else:
        raise SystemExit(
            "missing the curated co-receptor partner set. Fetch from S3:\n"
            "  processed://sprint1/s1e3-surface-library-filters/"
            "coreceptor_rationales.json\n"
            "  processed://sprint1/s1e3-surface-library-filters/"
            "coreceptor_calls.tsv\n"
            "Do NOT substitute the cohort record's raw co_receptor_partners "
            "wholesale — it is un-curated free text."
        )

    for rec in cohort_records:
        if rec["gene"] in curated_genes:
            continue
        partners.update(
            p for p in (rec.get("co_receptor_partners") or []) if p and SYMBOL_SHAPED.match(p)
        )

    partners = {p for p in partners if p}
    symbols = sorted(set(cohort) | partners)

    resolved = resolve(symbols, cache_dir=CACHE)
    tpm_by_ens, tpm_by_sym = _index_tpm(TPM_CSV)

    out = {}
    for sym in symbols:
        r, vals = _lookup(sym, resolved, tpm_by_ens, tpm_by_sym)
        out[sym] = {"ensembl": r.get("ensembl", []), **(vals or {c: None for c in NUM})}
    json.dump(out, open(OUT, "w"))
    cov = sum(1 for s in cohort if out[s]["max_expression_any"] is not None)
    print(f"resolved+joined {len(out)} symbols; cohort TPM coverage {cov}/{len(cohort)}")


if __name__ == "__main__":
    main()
