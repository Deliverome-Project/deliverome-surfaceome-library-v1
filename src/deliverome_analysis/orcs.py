"""Fetch and parse BioGRID ORCS genome-wide CRISPR screen hits.

A thin, reusable reader for the BioGRID ORCS webservice — pull the hit table of any
published screen by its ORCS screen id, as a tidy ``[gene, screen_id, score, hit]`` frame.

    from deliverome_analysis import orcs
    hits = orcs.parse_hits(orcs.fetch_screen("1234"), "1234")

REST: ``https://orcsws.thebiogrid.org/screen/{id}/?accesskey=...&format=json&hit=yes``.
The response is a JSON list of hit records; each carries OFFICIAL_SYMBOL, SCORE.N and a HIT flag.
Needs an ORCS access key in ``.env`` as ``BIOGRID_ORCS_ACCESS_KEY`` (or ``ORCS_ACCESS_KEY``). NB:
ORCS has its OWN key registry at https://orcsws.thebiogrid.org — a main BioGRID webservice key is
rejected. Resolve/join the returned symbols through ``deliverome_analysis.gene_resolver`` before
keying them against another dataset — ORCS reports display symbols, which drift.
"""

import os

import pandas as pd
import requests

WS = "https://orcsws.thebiogrid.org"


def _access_key() -> str:
    from deliverome_analysis.env import load_env

    load_env()  # idempotent; ensures .env is in os.environ regardless of caller
    key = os.environ.get("BIOGRID_ORCS_ACCESS_KEY") or os.environ.get("ORCS_ACCESS_KEY")
    if not key:
        raise RuntimeError(
            "No ORCS key. Set BIOGRID_ORCS_ACCESS_KEY in .env (register at "
            "https://orcsws.thebiogrid.org — separate from the main BioGRID key)."
        )
    return key.strip()


def fetch_screen(screen_id: str, hits_only: bool = True) -> dict:
    key = _access_key()
    params = {"accesskey": key, "format": "json"}
    if hits_only:
        params["hit"] = "yes"
    r = requests.get(f"{WS}/screen/{screen_id}/", params=params, timeout=90)
    r.raise_for_status()
    data = r.json()
    if isinstance(data, dict) and data.get("STATUS") == "ERROR":
        raise RuntimeError(f"ORCS error: {data.get('MESSAGE')}")
    return data


def parse_hits(payload, screen_id: str) -> pd.DataFrame:
    records = payload.values() if isinstance(payload, dict) else payload
    rows = []
    for r in records:
        if not isinstance(r, dict):
            continue
        gene = r.get("OFFICIAL_SYMBOL") or r.get("official_symbol")
        rows.append(
            {
                "gene": gene,
                "screen_id": str(screen_id),
                "score": r.get("SCORE.1") or r.get("score.1"),
                "hit": r.get("HIT") or r.get("hit"),
            }
        )
    # Explicit columns so an empty screen returns a well-formed frame, not a bare (column-less)
    # DataFrame that would make the dropna below raise KeyError.
    df = pd.DataFrame(rows, columns=["gene", "screen_id", "score", "hit"])
    return df.dropna(subset=["gene"]).reset_index(drop=True)
