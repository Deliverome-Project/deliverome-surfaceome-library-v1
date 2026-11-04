"""Read full therapeutic-antibody heavy chains from KEGG DRUG.

Thera-SAbDab carries only Fv sequences, so the CH3 allotype (DEL vs EEM at EU
356/358) and the hinge are not readable from it. KEGG DRUG (https://rest.kegg.jp)
stores the **full** heavy and light chains for many approved biologics, with no
login, from which both are directly readable. Coverage is partial — older mAbs
often have no SEQUENCE block — so callers must handle ``None``.

This module is the thin, reusable layer: resolve a name to a KEGG drug id, fetch
and cache the record, and extract the heavy chain, the CH3 allotype, and the
hinge. Trend analysis lives in experiments, not here.

Allotype is read from the invariant ``TLPPSR`` motif that immediately precedes
EU 356: the next three residues are 356/357/358, ``DEL`` (G1m1) or ``EEM``
(nG1m1). The hinge is the segment between the CH1 tail and the CH2 start
(``APE..GG``).
"""

from __future__ import annotations

import re
import urllib.parse
import urllib.request
from pathlib import Path

REST = "https://rest.kegg.jp"
_ALLOTYPE_RE = re.compile(r"TLPPSR(...)")
_HEAVY_SPLIT = "(Heavy chain)"
_LIGHT_SPLIT = "(Light"


def find_drug_id(name: str, timeout: int = 20) -> str | None:
    """Resolve a therapeutic name to its first KEGG DRUG id (``D#####``)."""
    text = _get(f"{REST}/find/drug/{urllib.parse.quote(name.lower())}", timeout)
    if not text.strip():
        return None
    return text.split("\t", 1)[0].replace("dr:", "").strip() or None


def fetch_record(drug_id: str, cache_dir: str | Path | None = None, timeout: int = 20) -> str:
    """Fetch a KEGG DRUG record, caching to ``cache_dir`` when given."""
    if cache_dir is not None:
        cached = Path(cache_dir) / f"{drug_id}.txt"
        if cached.exists():
            return cached.read_text()
    text = _get(f"{REST}/get/{drug_id}", timeout)
    if cache_dir is not None and text:
        cached = Path(cache_dir) / f"{drug_id}.txt"
        cached.parent.mkdir(parents=True, exist_ok=True)
        cached.write_text(text)
    return text


def heavy_chain(record: str) -> str | None:
    """Extract the full heavy-chain amino-acid sequence from a KEGG record."""
    if _HEAVY_SPLIT in record:
        block = record.split(_HEAVY_SPLIT, 1)[1].split(_LIGHT_SPLIT, 1)[0]
        seq = re.sub(r"[^A-Z]", "", block)
        return seq or None
    # Single-chain SEQUENCE block fallback: accept only if it reaches the Fc.
    m = re.search(r"\nSEQUENCE\s+(.+?)\n[A-Z]", record, re.S)
    if m:
        seq = re.sub(r"[^A-Z]", "", m.group(1))
        if "TLPPSR" in seq:
            return seq
    return None


def allotype(heavy: str | None) -> str | None:
    """Classify the CH3 allotype at EU 356/358: 'DEL', 'EEM', or the raw triple."""
    if not heavy:
        return None
    m = _ALLOTYPE_RE.search(heavy)
    if not m:
        return None
    triple = m.group(1)
    return {"DEL": "DEL", "EEM": "EEM"}.get(triple, triple)


def hinge(heavy: str | None) -> str | None:
    """Return the hinge segment (CH1 tail → CH2 start ``APExxGG``), if present."""
    if not heavy:
        return None
    m = re.search(r"(EPK[A-Z]{0,30}?DKTHTCPPCP)APE[LA]{2}GG", heavy)
    return m.group(1) if m else None


def _get(url: str, timeout: int) -> str:
    req = urllib.request.Request(url, headers={"User-Agent": "deliverome-analysis"})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            return resp.read().decode()
    except Exception:
        return ""
