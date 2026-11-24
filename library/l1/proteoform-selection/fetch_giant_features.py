"""Fetch UniProt sequence features for the oversized genes we intend to truncate.

A truncation is only safe if it is anchored on real coordinates -- the signal
peptide, the transmembrane span, the domain boundaries -- rather than on a
residue number copied out of a paper and hoped to still apply to the canonical
isoform. This pulls those features so every cut can be checked against them.
"""

from __future__ import annotations

import argparse
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

API = "https://rest.uniprot.org/uniprotkb/search"
WANTED = {"Signal", "Transmembrane", "Domain", "Repeat", "Topological domain", "Chain",
          "Region", "Site", "Propeptide"}


def fetch(gene: str) -> dict | None:
    q = f'gene_exact:{gene} AND organism_id:9606 AND reviewed:true'
    url = f"{API}?query={urllib.parse.quote(q)}&format=json&size=1"
    req = urllib.request.Request(url, headers={"User-Agent": "deliverome-analysis (research)"})
    with urllib.request.urlopen(req, timeout=60) as fh:
        data = json.load(fh)
    hits = data.get("results") or []
    if not hits:
        return None
    e = hits[0]
    feats = []
    for f in e.get("features", []):
        if f.get("type") not in WANTED:
            continue
        loc = f.get("location", {})
        s, t = loc.get("start", {}).get("value"), loc.get("end", {}).get("value")
        if s is None or t is None:
            continue
        feats.append({"type": f["type"], "start": s, "end": t,
                      "description": f.get("description", "")})
    return {
        "accession": e["primaryAccession"],
        "length": e["sequence"]["length"],
        "sequence": e["sequence"]["value"],
        "features": feats,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--genes", nargs="+", required=True)
    ap.add_argument("--out", default="giant_features.json")
    args = ap.parse_args()
    cache = Path(args.out)
    store = json.loads(cache.read_text()) if cache.exists() else {}
    for g in args.genes:
        if g in store:
            continue
        try:
            rec = fetch(g)
        except Exception as exc:                      # noqa: BLE001
            print(f"  {g:10} FAILED {type(exc).__name__}: {exc}")
            continue
        if rec is None:
            print(f"  {g:10} no reviewed human entry")
            continue
        store[g] = rec
        tm = [f for f in rec["features"] if f["type"] == "Transmembrane"]
        print(f"  {g:10} {rec['accession']:10} {rec['length']:6} aa  "
              f"{len(rec['features']):3} features, {len(tm)} TM")
        cache.write_text(json.dumps(store, indent=1, sort_keys=True))
        time.sleep(0.3)
    cache.write_text(json.dumps(store, indent=1, sort_keys=True))
    print(f"\n  {len(store)} genes cached in {cache}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
