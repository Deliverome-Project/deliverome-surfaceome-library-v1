#!/usr/bin/env python3
"""Fetch the genome-wide surfaceome catalog (19,324 rows) from the public API.

Writes catalog_rows.json — the row list final_library.py (LL tier: UniProt
db-bit) and reference_coverage.py (universe + uniprot->symbol map) read. Each
row carries: symbol, n_sources, db (bitmask over db_keys), uniprot, tr, deep_dive,
sb, intern, ddf, deep_dive_tier, deep_dive_facet, low_lit_uniprot.

Read-only, public, edge-cached endpoint (same as fetch_cohort.py). Stdlib-only.
Run: uv run python fetch_catalog.py   (or plain `python fetch_catalog.py`)
"""

import json
import os
import urllib.request

BASE = "https://api.deliverome.org/surfaceome/v1"
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "catalog_rows.json")


def main():
    req = urllib.request.Request(
        f"{BASE}/catalog", headers={"User-Agent": "deliverome-analysis/surface-library"}
    )
    with urllib.request.urlopen(req, timeout=60) as resp:
        d = json.load(resp)
    rows = d["rows"]  # the row LIST, not the {db_keys,...,rows} envelope
    with open(OUT, "w") as fh:
        json.dump(rows, fh)
    dbk = d.get("db_keys")
    print(f"wrote {len(rows)} catalog rows -> {OUT}")
    print(f"db_keys (bit order): {dbk}")
    print(f"row0 keys: {list(rows[0].keys())}")


if __name__ == "__main__":
    main()
