"""Peptide- and protein-level identification metrics for DIA searches, defined once.

Every QC page, QC record and experiment script computes these numbers through this module, so a
metric means the same thing in every run and can be compared across runs in
``processed://db/dia_qc.sqlite``. Inputs and outputs are pydantic models: the digest and
detectability rules a search used (:class:`DigestRules`), and one validated record of metrics
per run (:class:`RunPeptideMetrics`).

Definitions
-----------
- **Peptide**: a distinct stripped sequence. A modified and an unmodified form of one sequence,
  or one sequence at two charges, count once. I and L are not merged.
- **Precursor**: sequence + modifications + charge.
- **Protein group**: the search engine's group (alphaDIA ``pg.name``), named by its lead
  accession; the lead accession (isoform suffix kept) stands for the group's sequence.
- **Missed cleavage**: an internal K or R not followed by P (:func:`missed_cleavages`).
- **Points per peak**: MS2 data points across a precursor's elution peak, estimated as
  ``GAUSSIAN_BASE_PER_FWHM x FWHM / cycle time`` (a Gaussian's full width at 4 sigma is 1.699 x
  its FWHM).
- **Theoretical tryptic peptides** of a protein: its full trypsin digest (cleave after K/R unless
  followed by P) with ``missed_cleavages`` = 0, length within ``[min_length, max_length]``.
- **Detectable tryptic peptides**: the theoretical ones whose m/z, with carbamidomethyl C, falls
  inside ``mz_range`` at one of ``charges``. ``exclude_glyco_sequons`` (off by default) also drops
  peptides with an N-glycosylation sequon (N-x-S/T, x not P) from numerator and denominator.
  It is off because most sequons are unoccupied: in S3E12 HEK runs, sequon peptides of
  non-surface proteins were identified at 0.88x the rate of other peptides (surface: 0.33x).
- **Tryptic coverage** of a protein group: identified peptides that equal a theoretical
  (detectable) peptide of its lead accession, divided by the number of theoretical
  (detectable) peptides. Peptides with a missed cleavage therefore never count.
- **Gene-unique peptide**: occurs in exactly one gene's sequences in the search FASTA (isoforms of
  one gene count as one gene), using the search's own digest rules (missed cleavages included).
- **Surface**: accessible-surfaceome tier canonical or likely (the repo default), by the lead
  accession's base (isoform suffix stripped).
"""

from __future__ import annotations

import functools
import re
from collections import defaultdict
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from pydantic import BaseModel, ConfigDict, Field, model_validator

GAUSSIAN_BASE_PER_FWHM = 1.699  # 4 sigma / (2 sqrt(2 ln 2) sigma)
PROTON = 1.007276
WATER = 18.010565
# Monoisotopic residue masses; C carries the fixed carbamidomethyl (+57.021464).
RESIDUE_MASS: dict[str, float] = {
    "G": 57.021464, "A": 71.037114, "S": 87.032028, "P": 97.052764, "V": 99.068414,
    "T": 101.047679, "C": 160.030649, "L": 113.084064, "I": 113.084064, "N": 114.042927,
    "D": 115.026943, "Q": 128.058578, "K": 128.094963, "E": 129.042593, "M": 131.040485,
    "H": 137.058912, "F": 147.068414, "R": 156.101111, "Y": 163.063329, "W": 186.079313,
    "U": 150.953636, "O": 237.147727,
}  # fmt: skip
SURFACE_TIERS = ("canonical", "likely")
TIERS = ("canonical", "likely", "low / uncertain", "not surface")
_CLEAVE = re.compile(r"(?<=[KR])(?!P)")
N_GLYCO_SEQUON = re.compile(r"N[^P][ST]")


class DigestRules(BaseModel):
    """How peptides were generated and which a search could see (the library's settings)."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    enzyme: str = "trypsin"
    missed_cleavages: int = Field(1, ge=0, le=3, description="for the gene-unique map")
    min_length: int = Field(7, ge=1)
    max_length: int = Field(35, ge=1)
    charges: tuple[int, ...] = (2, 3, 4)
    mz_range: tuple[float, float] = (380.0, 980.0)
    exclude_glyco_sequons: bool = Field(
        False, description="detectable peptides exclude N-x-S/T sequons (glycans not searched)"
    )

    @model_validator(mode="after")
    def _check(self) -> DigestRules:
        if self.enzyme != "trypsin":
            raise ValueError("only trypsin is implemented")
        if self.min_length > self.max_length:
            raise ValueError("min_length > max_length")
        if not self.charges or min(self.charges) < 1:
            raise ValueError("charges must be positive")
        if self.mz_range[0] >= self.mz_range[1]:
            raise ValueError("mz_range must be (low, high)")
        return self


class RunPeptideMetrics(BaseModel):
    """One run's identification metrics (see the module docstring for every definition)."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    run: str
    protein_groups: int = Field(ge=0)
    precursors: int = Field(ge=0)
    peptides: int = Field(ge=0)
    gene_unique_peptides: int | None = Field(None, ge=0)
    peptides_per_protein_median: float | None = Field(None, ge=0)
    proteins_2plus_peptides: int = Field(ge=0)
    surface_protein_groups: int | None = Field(None, ge=0)
    peptides_per_surface_protein_median: float | None = Field(None, ge=0)
    surface_proteins_2plus_peptides: int | None = Field(None, ge=0)
    missed_cleavage_pct: float = Field(ge=0, le=100)
    missed_cleavage_surface_pct: float | None = Field(None, ge=0, le=100)
    fwhm_median_s: float | None = Field(None, ge=0)
    points_per_peak_median: float | None = Field(None, ge=0)
    points_per_peak_by_tier: dict[str, float] = Field(default_factory=dict)
    tryptic_coverage_surface_pct: float | None = Field(None, ge=0, le=100)
    tryptic_coverage_nonsurface_pct: float | None = Field(None, ge=0, le=100)
    detectable_coverage_surface_pct: float | None = Field(None, ge=0, le=100)
    detectable_coverage_nonsurface_pct: float | None = Field(None, ge=0, le=100)

    @model_validator(mode="after")
    def _consistent(self) -> RunPeptideMetrics:
        if self.peptides > self.precursors:
            raise ValueError("more peptides than precursors")
        if self.proteins_2plus_peptides > self.protein_groups:
            raise ValueError("more 2-peptide proteins than protein groups")
        unknown = set(self.points_per_peak_by_tier) - set(TIERS)
        if unknown:
            raise ValueError(f"unknown tiers {unknown}")
        return self

    def flat(self) -> dict[str, float]:
        """Numeric metrics as ``{name: value}`` for a QC record (None values left out)."""
        out: dict[str, float] = {}
        for k, v in self.model_dump(exclude={"run", "points_per_peak_by_tier"}).items():
            if v is not None:
                out[k] = float(v)
        for t, v in self.points_per_peak_by_tier.items():
            out[f"points_per_peak_{tier_key(t)}"] = float(v)
        return out


def tier_key(tier: str) -> str:
    """``low / uncertain`` -> ``low_uncertain`` (the QC record's metric-name spelling)."""
    return tier.replace(" / ", "_").replace(" ", "_")


# ----------------------------------------------------------------------------------- sequence
def cleave(seq: str) -> list[str]:
    """Full trypsin digest: cut after K or R unless followed by P."""
    return [p for p in _CLEAVE.split(seq) if p]


def digest(seq: str, rules: DigestRules, missed: int | None = None) -> list[str]:
    """Peptides of ``seq`` with up to ``missed`` (default ``rules.missed_cleavages``) missed
    cleavages, within the length limits, in sequence order (duplicates kept)."""
    m = rules.missed_cleavages if missed is None else missed
    pieces = cleave(seq)
    out = []
    for i in range(len(pieces)):
        for j in range(i, min(i + m + 1, len(pieces))):
            pep = "".join(pieces[i : j + 1])
            if rules.min_length <= len(pep) <= rules.max_length:
                out.append(pep)
    return out


def missed_cleavages(seq: str) -> int:
    """Internal K/R not followed by P (the C-terminal residue does not count)."""
    return sum(1 for i, a in enumerate(seq[:-1]) if a in "KR" and seq[i + 1] != "P")


def peptide_mass(seq: str) -> float:
    """Monoisotopic neutral mass with carbamidomethyl C; unknown residues raise KeyError."""
    return sum(RESIDUE_MASS[a] for a in seq) + WATER


def mz(seq: str, charge: int) -> float:
    return (peptide_mass(seq) + charge * PROTON) / charge


def has_glyco_sequon(seq: str) -> bool:
    """N-x-S/T with x not P (the N-glycosylation motif)."""
    return N_GLYCO_SEQUON.search(seq) is not None


def detectable(seq: str, rules: DigestRules) -> bool:
    """Inside ``rules.mz_range`` at one of ``rules.charges``; with
    ``rules.exclude_glyco_sequons``, also free of an N-glycosylation sequon."""
    if rules.exclude_glyco_sequons and has_glyco_sequon(seq):
        return False
    try:
        mass = peptide_mass(seq)
    except KeyError:  # X, B, Z ...: not predictable
        return False
    lo, hi = rules.mz_range
    return any(lo <= (mass + z * PROTON) / z <= hi for z in rules.charges)


def points_per_peak(fwhm_s: float, cycle_time_s: float) -> float:
    if cycle_time_s <= 0:
        raise ValueError("cycle time must be positive")
    return GAUSSIAN_BASE_PER_FWHM * fwhm_s / cycle_time_s


# ---------------------------------------------------------------------------------- reference
@dataclass(frozen=True)
class ProteinReference:
    """What the metrics need from the search FASTA: per accession its theoretical and detectable
    tryptic peptides, and per peptide the number of genes it occurs in."""

    rules: DigestRules
    theoretical: Mapping[str, frozenset[str]]
    detectable: Mapping[str, frozenset[str]]
    genes_per_peptide: Mapping[str, int] = field(repr=False)

    @classmethod
    def from_records(
        cls, records: Iterable[tuple[str, str, str]], rules: DigestRules | None = None
    ) -> ProteinReference:
        """``records``: (accession, gene, sequence) per FASTA entry."""
        rules = rules or DigestRules()
        theo: dict[str, frozenset[str]] = {}
        det: dict[str, frozenset[str]] = {}
        genes: dict[str, set[str]] = defaultdict(set)
        for acc, gene, seq in records:
            full = frozenset(digest(seq, rules, missed=0))
            theo[acc] = full
            det[acc] = frozenset(p for p in full if detectable(p, rules))
            for p in digest(seq, rules):
                genes[p].add(gene)
        return cls(rules, theo, det, {p: len(g) for p, g in genes.items()})

    @classmethod
    def from_fasta_text(cls, text: str, rules: DigestRules | None = None) -> ProteinReference:
        from deliverome_analysis.search_fasta import accession, parse_fasta

        def records() -> Iterable[tuple[str, str, str]]:
            for header, seq in parse_fasta(text):
                acc = accession(header)
                m = re.search(r"\bGN=(\S+)", header)
                yield acc, (m[1] if m else acc.split("-")[0]), seq

        return cls.from_records(records(), rules)


DEFAULT_SEARCH_FASTA = "human-sp-iso_cont-universal-and-cell-culture_irt/2026_03"


@functools.lru_cache(maxsize=2)
def default_reference(name: str = DEFAULT_SEARCH_FASTA) -> ProteinReference | None:
    """The published search FASTA's reference (cached), or None when S3 is unreachable."""
    from deliverome_analysis.search_fasta import download

    try:
        return ProteinReference.from_fasta_text(download(name).read_text())
    except Exception:  # no credentials / offline: FASTA-based metrics are left out
        return None


# ------------------------------------------------------------------------------------ metrics
PRECURSOR_COLUMNS = ("sequence", "mods", "charge", "protein_group")


def run_metrics(
    prec: pd.DataFrame,
    *,
    run: str,
    tiers: Mapping[str, str] | None = None,
    cycle_time_s: float | None = None,
    reference: ProteinReference | None = None,
) -> RunPeptideMetrics:
    """Metrics of one run's identified target precursors (already FDR-filtered).

    ``prec`` columns: ``sequence`` (stripped), ``mods``, ``charge``, ``protein_group`` (lead
    accession first, ``;``-joined), optionally ``fwhm_s``. ``tiers`` maps base accession ->
    surfaceome tier; without it the surface metrics are None. ``reference`` adds gene-unique
    peptides and coverage; ``cycle_time_s`` adds points per peak.
    """
    missing = [c for c in PRECURSOR_COLUMNS if c not in prec]
    if missing:
        raise ValueError(f"missing columns {missing}")
    p = prec.copy()
    p["lead"] = p["protein_group"].astype(str).str.split(";").str[0]
    p["tier"] = (
        p["lead"].str.split("-").str[0].map(dict(tiers)).fillna("not surface")
        if tiers is not None
        else "not surface"
    )
    peps = p.drop_duplicates("sequence")
    mc = peps["sequence"].map(missed_cleavages) > 0
    groups = p.groupby("protein_group").agg(
        seqs=("sequence", lambda s: frozenset(s)), lead=("lead", "first"), tier=("tier", "first")
    )
    groups["n"] = groups["seqs"].map(len)
    gsurf = groups["tier"].isin(SURFACE_TIERS)

    def med(s: pd.Series) -> float | None:
        return float(s.median()) if len(s.dropna()) else None

    values: dict[str, object] = {
        "run": run,
        "protein_groups": len(groups),
        "precursors": len(p.drop_duplicates(["sequence", "mods", "charge"])),
        "peptides": len(peps),
        "peptides_per_protein_median": med(groups["n"]),
        "proteins_2plus_peptides": int((groups["n"] >= 2).sum()),
        "missed_cleavage_pct": float(100 * mc.mean()) if len(peps) else 0.0,
    }
    if tiers is not None:
        psurf = peps["tier"].isin(SURFACE_TIERS)
        values |= {
            "surface_protein_groups": int(gsurf.sum()),
            "peptides_per_surface_protein_median": med(groups.loc[gsurf, "n"]),
            "surface_proteins_2plus_peptides": int((gsurf & (groups["n"] >= 2)).sum()),
            "missed_cleavage_surface_pct": float(100 * mc[psurf].mean()) if psurf.any() else None,
        }
    if "fwhm_s" in p:
        values["fwhm_median_s"] = med(p["fwhm_s"])
        if cycle_time_s:
            ppp = GAUSSIAN_BASE_PER_FWHM * p["fwhm_s"] / cycle_time_s
            values["points_per_peak_median"] = med(ppp)
            values["points_per_peak_by_tier"] = {
                t: float(v) for t, v in ppp.groupby(p["tier"]).median().items() if pd.notna(v)
            }
    if reference is not None:
        values["gene_unique_peptides"] = int(
            (peps["sequence"].map(reference.genes_per_peptide).fillna(0) == 1).sum()
        )
        for kind, table in (
            ("tryptic", reference.theoretical),
            ("detectable", reference.detectable),
        ):
            cov = np.array(
                [
                    len(s & t) / len(t) if (t := table.get(lead)) else np.nan
                    for s, lead in zip(groups["seqs"], groups["lead"], strict=True)
                ],
                dtype=float,
            )
            if tiers is not None:
                for name, mask in (
                    ("surface", gsurf.to_numpy()),
                    ("nonsurface", ~gsurf.to_numpy()),
                ):
                    c = cov[mask]
                    c = c[~np.isnan(c)]
                    values[f"{kind}_coverage_{name}_pct"] = (
                        float(100 * np.median(c)) if len(c) else None
                    )
    return RunPeptideMetrics.model_validate(values)
