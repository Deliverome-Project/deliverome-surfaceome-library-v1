"""Snapshot the curated surface architecture attribute from D1.

``surface_annotation.annotation_json`` carries ``filters.subcategory`` — the architecture
call (``GPI_anchored``, ``single_pass_T1``, ``single_pass_T2``, ``multi_pass``, ``GPCR``,
``tetraspanin``, ``other``) and ``filters.surface_call_reason``, both produced by the
annotation run with cited evidence behind them. It covers all 5,332 annotated genes, which
includes every one of the library's 2,369.

That attribute is the answer to "is this protein GPI-anchored", so nothing here predicts it
from sequence. The snapshot exists only so the selection runs offline and so a reviewer can
see exactly which call each gene had when the library was picked.

``json_extract`` over the whole table exceeds D1's limits, hence the paging.

    uv run --locked python fetch_surface_subcategory.py
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

DATABASE = "surfaceome_public"
WRANGLER_CWD = Path("~/Git/accessible-surfaceome").expanduser()
PAGE = 400
COLUMNS = ("gene_symbol", "uniprot_acc", "subcategory", "surface_call_reason")

SQL = """SELECT gene_symbol, uniprot_acc,
       json_extract(annotation_json,'$.filters.subcategory') AS subcategory,
       json_extract(annotation_json,'$.filters.surface_call_reason') AS surface_call_reason
FROM (SELECT * FROM surface_annotation ORDER BY gene_symbol LIMIT {page} OFFSET {offset});"""


def query(sql: str) -> list[dict]:
    proc = subprocess.run(
        ["npx", "--yes", "wrangler", "d1", "execute", DATABASE, "--remote", "--json",
         f"--command={sql}"],
        capture_output=True, text=True, cwd=WRANGLER_CWD,
    )
    if proc.returncode != 0:
        sys.stderr.write((proc.stderr or proc.stdout)[-1500:])
        raise SystemExit(proc.returncode)
    return json.loads(proc.stdout)[0]["results"]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="surface_subcategory.tsv")
    args = ap.parse_args()

    rows: list[dict] = []
    offset = 0
    while True:
        page = query(SQL.format(page=PAGE, offset=offset))
        if not page:
            break
        rows.extend(page)
        offset += PAGE
    rows.sort(key=lambda r: r["gene_symbol"] or "")

    with open(args.out, "w") as fh:
        fh.write("\t".join(COLUMNS) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in COLUMNS) + "\n")
    n_gpi = sum(r["subcategory"] == "GPI_anchored" for r in rows)
    print(f"wrote {args.out} ({len(rows)} genes, {n_gpi} GPI_anchored)")


if __name__ == "__main__":
    main()
