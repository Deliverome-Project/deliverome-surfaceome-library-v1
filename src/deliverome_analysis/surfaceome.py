"""Accessible-surfaceome membership — the repo-wide way to restrict to surface proteins.

**Read this before filtering any dataset to "surface proteins."**

There are two very different lists in ``deliverome-external`` under
``accessible-surfaceome/`` and it is easy to grab the wrong one:

* ``2026-07-10/candidate_universe.tsv`` with ``in_db_union == 1`` -> **5,680**
  accessions. This is the **pre-curation candidate net**: everything flagged by
  *any one* of the 7 merge sources (UniProt/GO/SURFY/CSPA/HPA + DeepTMHMM/
  COMPARTMENTS). It is the *input* to curation, not a surface-protein list. It
  contains ~1,078 proteins later curated **"no" (not surface)**, ~973 low-
  confidence, and hundreds never deep-dived at all — including secreted ECM
  (COL1A2, VCAN, PRELP, …). Do **not** use it as "the surfaceome."
* The **curated tiers** — ``canonical`` (**1,782**) and ``likely`` (1,243) —
  are the deep-dive-certified "this *is* a cell-surface protein" calls. Their
  union (``canonical + likely`` = **3,025**) is the membership base every
  surface-restricted analysis in this repo should use.

This helper serves the curated tiers, keyed on the stable ``uniprot_accession``,
so nobody has to re-derive them or reach for ``in_db_union`` again.

Usage
-----
    from deliverome_analysis.surfaceome import surface_accessions, canonical_accessions

    surf = surface_accessions()             # canonical + likely — the default (3,025)
    canon = canonical_accessions()          # canonical only (1,782)

    # restrict your own frame (join on the stable uniprot accession, not symbol):
    df = df[df["uniprot_accession"].isin(surf)]
    # or, isoform-suffix-tolerant:
    from deliverome_analysis.surfaceome import filter_to_surface
    df = filter_to_surface(df, "uniprot_accession")

Source of truth / versioning
----------------------------
Backed by a pinned snapshot of the live catalog
(``api.deliverome.org/surfaceome/v1/catalog``) exported to
``deliverome-external`` — see :data:`CATALOG_KEY`. **To move to a newer catalog
release, bump that one constant**; every caller follows. Reads go through the
cached S3 helper (``external.download``), so repeat calls are free.
"""

from __future__ import annotations

import json
import re
from functools import lru_cache

import pandas as pd

# --- pinned catalog version. Bump THIS to adopt a newer accessible-surfaceome release. ---
# Snapshot of the live catalog (universe cu_2026_05_12_post_silentdrop_fix).
# Tier counts at this pin: canonical=1782, likely=1243, low=973, no=1078, uncertain=54.
CATALOG_KEY = "accessible-surfaceome/api-catalog/2026-08-26/catalog.json"

# ``deep_dive_tier`` values are MUTUALLY EXCLUSIVE labels in the catalog (unlike
# the API's ``passes_likely`` preset, where canonical ⊂ likely). So the surface
# membership base is the *union* of the two positive tiers.
SURFACE_TIERS: tuple[str, ...] = ("canonical", "likely")

_ISOFORM_SUFFIX = re.compile(r"-\d+$")


def base_accession(acc: str) -> str:
    """Strip a UniProt isoform suffix (``P12345-2`` -> ``P12345``).

    The catalog and the ``candidate_universe`` join key are both base (isoform-
    stripped) primary accessions, so match on this form.
    """
    return _ISOFORM_SUFFIX.sub("", str(acc))


@lru_cache(maxsize=1)
def load_catalog() -> pd.DataFrame:
    """Load the pinned catalog snapshot as a tidy per-gene DataFrame.

    Columns: ``symbol``, ``uniprot_accession`` (base), ``tier``
    (``deep_dive_tier``: canonical/likely/low/no/uncertain), ``n_sources``.
    Only rows that were deep-dived (have a tier) are returned; the ~14k genome-
    wide rows with no deep dive are dropped. Cached in-process.
    """
    from deliverome_analysis.s3 import external

    path = external.download(CATALOG_KEY)
    with open(path) as fh:
        payload = json.load(fh)

    records = []
    for row in payload.get("rows", []):
        tier = row.get("deep_dive_tier")
        acc = row.get("uniprot")
        if tier is None or not acc:
            continue
        records.append(
            {
                "symbol": row.get("symbol"),
                "uniprot_accession": base_accession(acc),
                "tier": tier,
                "n_sources": row.get("n_sources"),
            }
        )
    # Kept faithful to the catalog (5,130 deep-dived gene rows; tier counts
    # canonical=1,782 / likely=1,243). A handful of paralog pairs share one
    # UniProt accession (PTPRZ1/PTPRZ2 -> P23471, olfactory/defensin families),
    # so ``accessions()`` — a set of accessions — is one or two smaller than the
    # gene-row count (canonical distinct-accession = 1,781). That is expected:
    # joins key on the accession, which cannot separate such paralogs anyway.
    return pd.DataFrame.from_records(records).reset_index(drop=True)


@lru_cache(maxsize=1)
def catalog_symbols() -> dict[str, str]:
    """Base UniProt accession -> gene symbol for every catalog row (~19k genes).

    For display (labels, hover text) only: join on the accession, never on the symbol.
    Unlike :func:`load_catalog` this keeps rows that were never deep-dived.
    """
    from deliverome_analysis.s3 import external

    with open(external.download(CATALOG_KEY)) as fh:
        rows = json.load(fh).get("rows", [])
    return {
        base_accession(r["uniprot"]): r["symbol"]
        for r in rows
        if r.get("uniprot") and r.get("symbol")
    }


def membership(tiers: tuple[str, ...] = SURFACE_TIERS) -> pd.DataFrame:
    """Catalog rows whose ``tier`` is in ``tiers`` (default: canonical + likely)."""
    df = load_catalog()
    return df[df["tier"].isin(tiers)].reset_index(drop=True)


def accessions(tiers: tuple[str, ...] = SURFACE_TIERS) -> frozenset[str]:
    """Set of base UniProt accessions in ``tiers`` (default: canonical + likely)."""
    return frozenset(membership(tiers)["uniprot_accession"])


def surface_accessions() -> frozenset[str]:
    """The recommended default surface set: canonical + likely.

    3,025 gene rows -> 3,024 distinct accessions (a paralog pair shares one).
    """
    return accessions(SURFACE_TIERS)


def canonical_accessions() -> frozenset[str]:
    """Canonical tier only — the high-confidence certified surface set.

    1,782 gene rows -> 1,781 distinct accessions (PTPRZ1/PTPRZ2 share P23471).
    """
    return accessions(("canonical",))


def filter_to_surface(
    df: pd.DataFrame,
    accession_col: str,
    tiers: tuple[str, ...] = SURFACE_TIERS,
) -> pd.DataFrame:
    """Return ``df`` restricted to rows whose accession is a surface protein.

    Isoform-suffix tolerant: ``P12345-2`` matches on ``P12345``. Joins on the
    stable UniProt accession — never on a bare gene symbol (symbols drift).
    """
    allowed = accessions(tiers)
    mask = df[accession_col].map(lambda a: base_accession(a) in allowed)
    return df[mask].reset_index(drop=True)
