"""Resolve cell-line names to Cellosaurus CVCL ids + parse culture metadata.

A reusable reader for the Cellosaurus REST API — turn a cell-line *name* into its canonical
record and pull the fields we reuse across experiments (CVCL accession, disease, derived-from
site, doubling time, medium).

    from deliverome_analysis import cellosaurus
    rec = cellosaurus.resolve("K562")
    meta = cellosaurus.parse_metadata(rec)   # {'cvcl_id': 'CVCL_0004', 'doubling_time_h': 19.5, ...}

Never hand-type CVCL accessions — always resolve by name via the API. ``resolve`` aggregates
across hyphenation variants (Cellosaurus stores ``K562`` as ``K-562``) and prefers the record whose
RECOMMENDED name matches, so it picks the canonical line over derivatives that merely contain the
query. Medium is not a standard Cellosaurus field (usually absent); fill it from ATCC/literature
downstream.
"""

import re

import requests

API = "https://api.cellosaurus.org"

_FIELDS = "id,ac,ca,di,derived-from-site,doubling-time,cc,sy"


def _norm(s: str) -> str:
    return re.sub(r"[^a-z0-9]", "", str(s).lower())


def _query_variants(name: str) -> list[str]:
    """Cellosaurus canonical names are hyphenated (K562 -> K-562); q=K562 misses them.
    Try the raw name plus a letter/digit-hyphenated variant."""
    variants = [name]
    hyph = re.sub(r"([A-Za-z])(\d)", r"\1-\2", name)
    if hyph != name:
        variants.append(hyph)
    return variants


def search(name: str, rows: int = 20) -> dict:
    r = requests.get(
        f"{API}/search/cell-line",
        params={"q": name, "format": "json", "fields": _FIELDS, "rows": rows},
        timeout=30,
    )
    r.raise_for_status()
    return r.json()


def records(payload: dict) -> list[dict]:
    return payload.get("Cellosaurus", {}).get("cell-line-list", [])


def _recommended_name(rec: dict) -> str:
    for nm in rec.get("name-list", []):
        if isinstance(nm, dict) and nm.get("type") == "identifier":
            return nm.get("value", "")
    nl = rec.get("name-list", [])
    return (nl[0].get("value") if nl and isinstance(nl[0], dict) else "") or ""


def resolve(name: str, rows: int = 30) -> dict | None:
    """Return the canonical Cellosaurus record for `name`.

    Aggregates records across hyphenation variants and picks the one whose RECOMMENDED
    name matches (normalized) — which selects the canonical line over derivatives that
    merely contain the query in their name.
    """
    target = _norm(name)
    seen: dict[str, dict] = {}
    for q in _query_variants(name):
        for rec in records(search(q, rows=rows)):
            acc = _accession(rec)
            seen.setdefault(acc, rec)
    exact = [r for r in seen.values() if _norm(_recommended_name(r)) == target]
    if exact:
        # shortest recommended name = the base line, not a longer-named variant
        return min(exact, key=lambda r: len(_recommended_name(r)))
    return next(iter(seen.values()), None)


def _accession(rec: dict) -> str:
    accs = rec.get("accession-list", [])
    for a in accs:
        if a.get("type") == "primary":
            return a["value"]
    return accs[0]["value"] if accs else ""


def _doubling_time(rec: dict) -> tuple[str | None, float | None]:
    lst = rec.get("doubling-time-list") or []
    if not lst:
        return None, None
    raw = lst[0].get("doubling-time-value")
    nums = [float(x) for x in re.findall(r"\d+(?:\.\d+)?", raw or "")]
    hours = (sum(nums) / len(nums)) if nums else None  # midpoint of a range
    return raw, hours


def _derived_site(rec: dict) -> str | None:
    lst = rec.get("derived-from-site-list") or []
    return (lst[0].get("site") or {}).get("value") if lst else None


def _disease(rec: dict) -> str | None:
    lst = rec.get("disease-list") or []
    if not lst:
        return None
    d = lst[0]
    return d.get("value") or (d.get("cv-term") or {}).get("value") or d.get("label")


def _medium(rec: dict) -> str | None:
    for cc in rec.get("comment-list", []):
        if cc.get("category") in ("Medium", "Culture conditions"):
            return cc.get("value")
    return None


def parse_metadata(rec: dict) -> dict:
    raw_dt, hours = _doubling_time(rec)
    return {
        "cvcl_id": _accession(rec),
        "category": rec.get("category"),
        "disease": _disease(rec),
        "derived_from_site": _derived_site(rec),
        "doubling_time_raw": raw_dt,
        "doubling_time_h": hours,
        "medium": _medium(rec),
    }
