"""Gene-symbol resolver (MyGene.info) — repo-wide canonicalization.

Maps gene symbols / aliases / previous symbols to the HGNC-approved symbol +
Ensembl gene id(s), so analyses in this repo join to external gene-keyed data
(expression, orthology, …) consistently and alias-drift (e.g. ``CDC50A`` ->
approved ``TMEM30A``) doesn't cause silent join misses.

**Why this and not ``import accessible_surfaceome``:** the accessible-surfaceome
catalog resolves identity with MyGene.info via
``accessible_surfaceome.merge.gene_symbols.resolve_gene_symbols_with_mygene``.
That function is DataFrame-shaped and tied to their merge pipeline, and the
package is a heavy private dependency. Rather than couple to it, this helper
reuses the **same MyGene config** (scopes / fields / species below) so results
match the catalog, with only ``mygene`` as a dependency.

Results cache to a JSON keyed by the input-symbol-set hash (re-runs are free),
mirroring the catalog's caching convention.
"""

from __future__ import annotations

import hashlib
import json
import os
from datetime import UTC, datetime
from pathlib import Path

# Mirror accessible_surfaceome.merge.gene_symbols (keep in sync if they change).
_SCOPES = ["symbol", "alias", "prev_symbol"]
_FIELDS = ["ensembl.gene", "symbol"]
_SPECIES = "human"


def _default_cache_dir() -> Path:
    root = os.environ.get("DELIVEROME_CACHE_DIR")
    base = Path(root) if root else Path.cwd() / "data" / "external"
    return base / "mygene_symbol_resolution"


def resolve(symbols, cache_dir: str | os.PathLike | None = None) -> dict[str, dict]:
    """Resolve gene symbols/aliases to canonical symbol + Ensembl gene ids.

    Returns ``{query_symbol: {"symbol": approved|None, "ensembl": [ids],
    "ambiguous": bool, "notfound": bool}}``. Deterministic in the symbol set;
    the MyGene response is cached on disk keyed by the input hash.
    """
    symbols = sorted({s for s in symbols if s})
    cache_dir = Path(cache_dir) if cache_dir else _default_cache_dir()
    cache_path = cache_dir / "mygene_response.json"
    input_hash = hashlib.sha256(json.dumps(symbols, separators=(",", ":")).encode()).hexdigest()

    hits = None
    if cache_path.exists():
        try:
            cached = json.loads(cache_path.read_text())
            if cached.get("input_sha256") == input_hash:
                hits = cached.get("raw_hits")
        except json.JSONDecodeError, OSError:
            hits = None
    if hits is None:
        import mygene  # imported lazily so the dep is only needed when resolving

        hits = mygene.MyGeneInfo().querymany(
            symbols, scopes=_SCOPES, fields=_FIELDS, species=_SPECIES
        )
        cache_dir.mkdir(parents=True, exist_ok=True)
        cache_path.write_text(
            json.dumps(
                {
                    "input_sha256": input_hash,
                    "fetched_at_utc": datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ"),
                    "raw_hits": hits,
                },
                indent=2,
            )
        )

    # Group every hit by its query before choosing, because the choice depends
    # on the whole candidate set — see _select_hits.
    grouped: dict[str, list[dict]] = {}
    for hit in hits:
        grouped.setdefault(hit.get("query"), []).append(hit)

    out: dict[str, dict] = {}
    for q, q_hits in grouped.items():
        found = [h for h in q_hits if not h.get("notfound")]
        if not found:
            out[q] = {"symbol": None, "ensembl": [], "ambiguous": False, "notfound": True}
            continue

        chosen, ambiguous = _select_hits(q, found)
        symbols = [h.get("symbol") for h in chosen if h.get("symbol")]
        ensembl: set[str] = set()
        for hit in chosen:
            ens = hit.get("ensembl")
            if isinstance(ens, dict):
                ensembl.add(ens.get("gene"))
            elif isinstance(ens, list):
                for e in ens:
                    if isinstance(e, dict):
                        ensembl.add(e.get("gene"))

        out[q] = {
            "symbol": symbols[0] if symbols else None,
            "ensembl": sorted(x for x in ensembl if x),
            "ambiguous": ambiguous,
            "notfound": False,
        }
    return out


def _ensembl_ids(hit: dict) -> set[str]:
    """Ensembl gene ids on one hit. MyGene returns a dict or a list of dicts."""
    ens = hit.get("ensembl")
    if isinstance(ens, dict):
        return {ens.get("gene")} - {None}
    if isinstance(ens, list):
        return {e.get("gene") for e in ens if isinstance(e, dict)} - {None}
    return set()


def _select_hits(query: str, found: list[dict]) -> tuple[list[dict], bool]:
    """Pick which MyGene hits define a query's record.

    **An exact match on the approved symbol wins outright.** MyGene ranks
    `querymany` results by text relevance, not by which scope matched, so a
    query that IS an approved symbol can be outranked by an unrelated gene
    that merely retired the same string as an alias. Real cases, all with the
    wrong gene scoring higher than the right one:

    ==========  ==================================  =====================
    Query       Outranking hit                      Why it collides
    ==========  ==================================  =====================
    ``CNR1``    ``PCDHA4`` / ``PCDHA@``             protocadherin alphas
                                                    were once "cadherin-
                                                    related neuronal
                                                    receptors"
    ``CNR2``    ``PCDHA6``                          same
    ``HRH1``    ``DHX8``                            once "human RNA
                                                    helicase 1"
    ``ADRA1A``  ``ADRA1D``                          shared alias history
    ==========  ==================================  =====================

    Taking the top hit made the resolver return an unrelated gene, and
    unioning every hit's Ensembl ids put unrelated genes into the join key —
    silently, since ``to_ensembl`` / ``to_canonical`` drop the ambiguity flag.
    An approved symbol must therefore beat somebody else's retired alias.

    Returns the hits to use and whether the query is genuinely ambiguous. An
    exact match is *not* ambiguous: the query is an approved symbol, and the
    other candidates are only historical aliases of other genes.
    """
    exact = [h for h in found if h.get("symbol") == query]

    # Prefer exact hits that actually carry a join key. MyGene sometimes keeps a
    # stale record whose `symbol` is a retired string with no Ensembl mapping,
    # alongside the live gene: `LGTN` (renamed EIF2D) and `TRB` (the TCR-beta
    # locus, distinct from THRB) both do this. Letting such a hit win returns a
    # record with no Ensembl id, which joins to nothing — the same silent drop
    # this function exists to prevent, just in the other direction.
    exact_with_key = [h for h in exact if _ensembl_ids(h)]
    if exact_with_key:
        return exact_with_key, len({h.get("symbol") for h in exact_with_key}) > 1
    if exact and any(_ensembl_ids(h) for h in found):
        # The approved-symbol hit has no join key but another candidate does.
        # Fall back so the row can still join, and flag it: we are no longer
        # certain the key belongs to the queried symbol.
        return found, True
    if exact:
        return exact, len({h.get("symbol") for h in exact}) > 1

    # No hit is an approved symbol, so the query is an alias or previous
    # symbol. Keep MyGene's ranking, but flag disagreement between hits.
    return found, len({h.get("symbol") for h in found if h.get("symbol")}) > 1


def to_ensembl(symbols, cache_dir=None) -> dict[str, list[str]]:
    """Convenience: ``{symbol: [ensembl_gene_ids]}``."""
    return {q: r["ensembl"] for q, r in resolve(symbols, cache_dir).items()}


def to_canonical(symbols, cache_dir=None) -> dict[str, str | None]:
    """Convenience: ``{symbol: approved_symbol_or_None}``."""
    return {q: r["symbol"] for q, r in resolve(symbols, cache_dir).items()}
