"""The GPI-anchored gene set.

The primary source is the curated architecture attribute already on D1:
``surface_annotation.annotation_json -> filters.subcategory == "GPI_anchored"``. It is an
evidence-backed call made per gene during annotation, it covers all 2,369 library genes,
and it is the same field the rest of the surfaceome database bins architecture by. Nothing
here predicts anchoring from sequence — the attribute exists, so it is read, not re-derived.
``fetch_surface_subcategory.py`` refreshes the local snapshot.

Two smaller sources are unioned in because the D1 set is not quite a superset:

* UniProt keyword KW-0336. Seven of its genes are not ``GPI_anchored`` on D1, three of them
  in the library (BST2, NCAM1, TNFRSF10C). Each has both an anchored and a transmembrane
  form, so neither call is wrong — the union just keeps the anchored reading available.
``gpi_exclusions.tsv`` subtracts the other way, for genes a source calls anchored that the
evidence contradicts. It is deliberately narrow: a gene only leaves when the sources
disagree with each other, never on one predictor's say-so. Where NetGPI is alone in
dissenting it is NetGPI that is wrong -- ``netgpi_blind_spots.tsv`` records those, and they
stay in the set.

* ``gpi_supplement.tsv``, one cited row per gene, for anchors both sources missed. The
  alpha2delta subunits are why it exists: Davies et al. 2010 (PNAS, PMID 20080692) show GPI
  anchoring for alpha2delta-1, -2 and -3 directly, while UniProt still annotates them as
  single-pass type I. D1 carries alpha2delta-1 and -3 but files alpha2delta-2 as ``other``.

This matters beyond bookkeeping. A GPI protein has a signal peptide and no transmembrane
helix, so DeepTMHMM calls it ``SP`` — correctly. Without the anchor attribute those genes
read as secreted, which is how CACNA2D3 came to be described as a case of the predictor
being wrong when the predictor was right and the annotation was stale.

The supplement is deliberately a file with a citation per row, not a hardcoded list: adding
a gene is a claim about the literature and should be reviewable as one. Membership in the
library is still curation's decision, not this module's — ``gpi_genes()`` only answers
whether a protein is anchored, which is what every topology guard in this analysis needs.
"""

from __future__ import annotations

import csv
import difflib
from pathlib import Path

KEYWORD_FILE = "uniprot_human_gpi.tsv"
SUPPLEMENT_FILE = Path(__file__).with_name("gpi_supplement.tsv")
SUBCATEGORY_FILE = Path(__file__).with_name("surface_subcategory.tsv")
EXCLUSION_FILE = Path(__file__).with_name("gpi_exclusions.tsv")
GPI_SUBCATEGORY = "GPI_anchored"


def _column(path: Path, key: str, where=None) -> set[str]:
    if not path.exists():
        return set()
    with path.open() as fh:
        return {
            (row.get(key) or "").strip()
            for row in csv.DictReader(fh, delimiter="\t")
            if where is None or where(row)
        } - {""}


def surface_subcategory(path: Path | str = SUBCATEGORY_FILE) -> dict[str, str]:
    """gene -> architecture subcategory, as the annotation run called it."""
    path = Path(path)
    if not path.exists():
        return {}
    with path.open() as fh:
        return {
            (r["gene_symbol"] or "").strip(): (r["subcategory"] or "").strip()
            for r in csv.DictReader(fh, delimiter="\t")
        }


#: The GPI set of record: NetGPI 1.1 per form unioned with curation under the CONFIDENT
#: rule, built by s4e9 ``build_gpi_set.py`` (Becca, 2026-11-01: NetGPI per form is
#: authoritative for every stage). Within the 2,368-gene library it is the curated union
#: below plus GFRAL and RAET1G.
GPI_SET = Path(__file__).parents[3] / "analysis/sprint4/s4e9-orf-design-fg-panel/gpi_set.tsv"


def gpi_genes(gpi_set: str | Path = GPI_SET) -> set[str]:
    """Gene symbols GPI-anchored in any form, from the NetGPI set of record.

    Raises if the file is missing: an empty set would silently reclassify every GPI gene.
    """
    path = Path(gpi_set)
    if not path.exists():
        raise FileNotFoundError(f"GPI set of record not found: {path}")
    out = _column(path, "gene_symbol", lambda r: (r.get("gpi_anchored") or "") == "1")
    if not out:
        raise ValueError(f"{path} lists no GPI-anchored gene")
    return out


def curated_gpi_genes(
    cache_dir: str | Path = "~/.cache/deliverome/appris",
    supplement: str | Path = SUPPLEMENT_FILE,
    subcategory: str | Path = SUBCATEGORY_FILE,
    exclusions: str | Path = EXCLUSION_FILE,
) -> set[str]:
    """The curated union: the D1 attribute, plus KW-0336 and the supplement.

    No longer the set of record. It still defines which genes' forms are sent to NetGPI
    (``build_*_gpi_fasta.py``), and it is what every stage used until 2026-11-01.
    """
    return (
        _column(Path(subcategory), "gene_symbol",
                lambda r: (r.get("subcategory") or "").strip() == GPI_SUBCATEGORY)
        | _column(Path(cache_dir).expanduser() / KEYWORD_FILE, "Gene Names (primary)")
        | _column(Path(supplement), "gene")
    ) - _column(Path(exclusions), "gene")


# The anchor attribute is annotated per gene, from the canonical human form, so an isoform
# inherits it only if it keeps the C-terminus that carries the signal. 30 residues is the
# measured span of the signal (UniProt's omega-site annotations put it at 18-39, median 25).
CTERM_WINDOW = 30
CTERM_MIN_IDENTITY = 0.5


def anchor_retained(
    form_seq: str | None,
    reference_seq: str | None,
    window: int = CTERM_WINDOW,
    min_identity: float = CTERM_MIN_IDENTITY,
) -> bool | None:
    """Does this form keep the reference's C-terminal anchor region? None if unknown.

    **Within one species only.** An isoform either keeps the anchor-bearing exon or swaps it
    for another, so the comparison is bimodal: across the 62 GPI-gene isoforms with a
    UniProt sequence, 46 score >= 0.8 and 14 score < 0.3, with 2 in between. The low group
    is the known soluble isoforms -- secreted DAF (CD55-2), soluble uPAR (PLAUR-2), CD59-2,
    the ACHE readthrough form -- so the test recovers its own positive controls.

    It does **not** transfer across species, and the cutoff above must not be applied to an
    ortholog. The GPI signal is conserved as architecture, not as sequence: human NYX versus
    mouse Nyx, which differ in anchoring (PMID 15905181), scores 0.37, while human versus
    mouse CD55, which are both anchored, scores 0.43. The two populations overlap, so no
    identity threshold separates them.
    """
    if not form_seq or not reference_seq:
        return None
    x, y = reference_seq[-window:], form_seq[-window:]
    sm = difflib.SequenceMatcher(None, x, y, autojunk=False)
    identity = sum(b.size for b in sm.get_matching_blocks()) / max(len(x), len(y))
    return identity >= min_identity
