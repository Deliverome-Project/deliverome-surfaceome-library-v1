"""Classify protein groups as human, contaminant, iRT standard or watchlist.

Searches against the published database (``external://search-fasta/...``) report contaminant
entries as ``Cont_<accession>`` and the iRT kit as ``Biognosys_iRT``. Every precursor's
protein group gets one class:

=============  ===============================================================================
``irt``        the iRT standard entry only: a spike-in, never a sample protein
``contaminant`` every accession is a contaminant entry
``mixed``      human and contaminant accessions share the precursor's peptides; the human
               protein is real, but some of its signal may come from the contaminant
``watchlist``  a human protein the contaminant libraries list (keratins, KRTAPs, filaggrin,
               ...); it may be endogenous or dust, so it is flagged, not removed
``human``      everything else
=============  ===============================================================================

The trimmer keeps ``contaminant`` and ``irt`` precursors out of the matrix limpa quantifies,
keeps ``mixed`` and ``watchlist`` in it, and records every class in the QC table so QC can
report how much of each run's signal they carry.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

CONTAMINANT_PREFIXES = ("Cont_", "CON__")  # our database / Frankenfield; MaxQuant convention
IRT_ACCESSIONS = frozenset({"Biognosys_iRT"})
# Human proteins the contaminant libraries list, from the published search database.
WATCHLIST_KEY = (
    "search-fasta/human-sp-iso_cont-universal-and-cell-culture_irt/2026_03/"
    "human_contaminant_watchlist.tsv"
)
CLASSES = ["human", "watchlist", "mixed", "contaminant", "irt"]
EXCLUDE_FROM_QUANT = ("contaminant", "irt")


def _base(acc: str) -> str:
    return acc.split("-")[0]


def is_contaminant(accession: str) -> bool:
    return accession.startswith(CONTAMINANT_PREFIXES)


def classify_group(accessions: str, watchlist: frozenset[str] = frozenset()) -> str:
    """One protein group (``;``-joined accessions) -> its class (see module doc)."""
    accs = [a.strip() for a in str(accessions).split(";") if a.strip()]
    if accs and all(a in IRT_ACCESSIONS for a in accs):
        return "irt"
    cont = [is_contaminant(a) for a in accs]
    if accs and all(cont):
        return "contaminant"
    if any(cont):
        return "mixed"
    if watchlist and any(_base(a) in watchlist for a in accs):
        return "watchlist"
    return "human"


def classify(accessions: pd.Series, watchlist: frozenset[str] = frozenset()) -> pd.Series:
    """Vectorised :func:`classify_group` over unique groups."""
    uniq = {g: classify_group(g, watchlist) for g in pd.unique(accessions.astype(str))}
    return accessions.astype(str).map(uniq)


def load_watchlist(key: str = WATCHLIST_KEY) -> frozenset[str]:
    """Human accessions on the published contaminant watchlist (isoform suffix removed)."""
    from deliverome_analysis.s3 import external

    w = pd.read_csv(Path(external.download(key)), sep="\t")
    return frozenset(_base(str(a)) for a in w["human_accession"].dropna())
