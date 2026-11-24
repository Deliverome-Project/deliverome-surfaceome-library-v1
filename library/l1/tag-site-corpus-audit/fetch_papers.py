"""Resolve every cited PMID/DOI behind the literature tag sites and pull open-access text.

Reads ``api_rows.json`` (from ``fetch_sites.py``) and writes, under ``--out``:

* ``papers.tsv``  one row per cited PMID or DOI: resolved title, first author, year,
  journal, PMCID, Europe PMC open-access flag, and whether full text was fetched.
* ``fulltext/<id>.xml`` Europe PMC full-text XML for open-access papers (untrusted
  downloaded data; read it with ``python -I`` from outside this folder).
* ``abstracts/<id>.txt`` the abstract for every resolved paper.

Run::

    uv run --locked python library/l1/tag-site-corpus-audit/fetch_papers.py \
        --api-rows library/l1/tag-site-corpus-audit/data/api_rows.json \
        --out library/l1/tag-site-corpus-audit/data
"""

from __future__ import annotations

import argparse
import csv
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/"
UA = {"User-Agent": "deliverome-analysis tag-site audit"}


def get(url: str, retries: int = 4) -> bytes | None:
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=90) as r:
                return r.read()
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None
            if attempt == retries - 1:
                raise
        except Exception:  # noqa: BLE001
            if attempt == retries - 1:
                raise
        time.sleep(2 * (attempt + 1))
    return None


def search(query: str) -> dict | None:
    url = EPMC + "search?" + urllib.parse.urlencode(
        {"query": query, "format": "json", "resultType": "core", "pageSize": 1})
    data = json.loads(get(url) or b"{}")
    hits = data.get("resultList", {}).get("result", [])
    return hits[0] if hits else None


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--api-rows", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()
    (args.out / "fulltext").mkdir(parents=True, exist_ok=True)
    (args.out / "abstracts").mkdir(parents=True, exist_ok=True)

    rows = json.loads(args.api_rows.read_text())
    ids: dict[str, dict] = {}
    for r in rows:
        for s in r.get("sources") or []:
            key = f"PMID:{s['pmid']}" if s.get("pmid") else s.get("citation", "")
            ids.setdefault(key, {"sites": [], "cited_title": s.get("title"), "cited_authors": s.get("authors")})
            ids[key]["sites"].append(r["site_id"])

    cols = ["cited_id", "n_sites", "cited_title", "cited_authors", "resolved_pmid", "resolved_doi",
            "resolved_title", "first_author", "year", "journal", "pmcid", "is_open_access",
            "fulltext_fetched", "sites"]
    out = []
    for key, info in sorted(ids.items()):
        if key.startswith("PMID:"):
            hit = search(f"EXT_ID:{key[5:]} AND SRC:MED")
        elif key.startswith("DOI "):
            hit = search(f'DOI:"{key[4:]}"')
        else:
            hit = None
        rec = {"cited_id": key, "n_sites": len(info["sites"]), "cited_title": info["cited_title"],
               "cited_authors": info["cited_authors"], "sites": ";".join(info["sites"]),
               "fulltext_fetched": 0}
        if hit:
            rec.update({"resolved_pmid": hit.get("pmid"), "resolved_doi": hit.get("doi"),
                        "resolved_title": hit.get("title"), "first_author": hit.get("authorString", "")[:80],
                        "year": hit.get("pubYear"), "journal": hit.get("journalTitle"),
                        "pmcid": hit.get("pmcid"), "is_open_access": hit.get("isOpenAccess")})
            safe = key.replace(":", "_").replace("/", "_").replace(" ", "_")
            if hit.get("abstractText"):
                (args.out / "abstracts" / f"{safe}.txt").write_text(hit["abstractText"])
            target = args.out / "fulltext" / f"{safe}.xml"
            if hit.get("pmcid") and not target.exists():
                try:
                    xml = get(EPMC + f"{hit['pmcid']}/fullTextXML", retries=2)
                except Exception as e:  # noqa: BLE001 - embargoed/non-OA PMC entries return 5xx
                    print(f"{key}: full text not retrievable ({e})")
                    xml = None
                if xml:
                    target.write_bytes(xml)
            rec["fulltext_fetched"] = int(target.exists())
        out.append(rec)
        time.sleep(0.3)

    with (args.out / "papers.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out)
    n_ft = sum(r["fulltext_fetched"] for r in out)
    print(f"{len(out)} cited ids, {sum(1 for r in out if r.get('resolved_pmid'))} resolved, {n_ft} full texts")


if __name__ == "__main__":
    main()
