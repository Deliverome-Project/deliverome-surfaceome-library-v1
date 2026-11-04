"""Loader for Thera-SAbDab, the Therapeutic Structural Antibody Database (OPIG).

Thera-SAbDab (https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/therasabdab/)
is a living catalogue of ~1,100 therapeutic antibodies with, per entry, the
WHO-INN name, format, CH1 isotype, highest clinical phase, Fv sequences, PDB
hits, INN proposal/recommendation years, target, and sponsor.

Scope note (read before using for Fc work): the sequences are **Fv only** (the
variable domains, ~110-140 aa). There is **no constant region, hinge, or G1m
allotype** in the catalogue — `ch1_isotype` is an isotype label (G1/G2/G4), not
an allotype. Isotype and format *trends over time* are answerable from this data;
per-antibody allotype (DEL vs EEM) and hinge are not.

The catalogue updates in place upstream; pin an analysis by saving the CSV
(`download()`) and recording the fetch date, or mirror it to
`external://therasabdab/<fetch-date>/`.
"""

from __future__ import annotations

import csv
import re
import urllib.request
from dataclasses import dataclass
from pathlib import Path

CATALOG_URL = (
    "https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/static/downloads/"
    "TheraSAbDab_SeqStruc_OnlineDownload.csv"
)

_YEAR_RE = re.compile(r"(?:19|20)\d\d")


@dataclass(frozen=True)
class Therapeutic:
    """One Thera-SAbDab entry. `raw` keeps the full original row."""

    name: str
    fmt: str
    ch1_isotype: str
    highest_trial: str
    status: str
    heavy: str
    light: str
    year_proposed: int | None
    year_recommended: int | None
    target: str
    companies: str
    raw: dict

    @property
    def isotype(self) -> str:
        """Normalise CH1 isotype to IgG1/IgG2/IgG3/IgG4/mixed/other/na."""
        v = (self.ch1_isotype or "").strip()
        if ";" in v:  # bispecific with two chains
            return "mixed"
        v = v.upper()
        if v in ("", "NA", "NONE"):
            return "na"
        for n in ("1", "2", "3", "4"):
            if v.startswith(f"G{n}"):
                return f"IgG{n}"
        return "other"

    @property
    def year(self) -> int | None:
        """Best available INN year: proposed, else recommended."""
        return self.year_proposed if self.year_proposed is not None else self.year_recommended

    @property
    def is_whole_mab(self) -> bool:
        return "whole mab" in (self.fmt or "").lower()

    @property
    def is_bispecific(self) -> bool:
        return "bispecific" in (self.fmt or "").lower()

    @property
    def is_approved(self) -> bool:
        return "approved" in (self.highest_trial or "").lower()


def _year(value: str) -> int | None:
    m = _YEAR_RE.search(value or "")
    return int(m.group(0)) if m else None


def download(dest: str | Path, url: str = CATALOG_URL, timeout: int = 120) -> Path:
    """Fetch the Thera-SAbDab CSV to ``dest`` and return its path.

    The catalogue is a living resource; record the fetch date for reproducibility.
    """
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    req = urllib.request.Request(url, headers={"User-Agent": "deliverome-analysis"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        dest.write_bytes(resp.read())
    return dest


def _col(fieldnames: list[str], *needles: str) -> str:
    """Find the column whose name contains all needles (case-insensitive)."""
    for name in fieldnames:
        low = name.lower()
        if all(n.lower() in low for n in needles):
            return name
    raise KeyError(f"no column matching {needles} in {fieldnames}")


def load(path: str | Path) -> list[Therapeutic]:
    """Parse a Thera-SAbDab CSV into :class:`Therapeutic` records."""
    with Path(path).open(newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh)
        fns = reader.fieldnames or []
        c_name = _col(fns, "therapeutic")
        c_fmt = _col(fns, "format")
        c_iso = _col(fns, "isotype")
        c_trial = _col(fns, "highest")
        c_status = _col(fns, "est", "status")
        c_heavy = _col(fns, "heavysequence")  # the non-bispecific heavy column sorts first
        c_light = _col(fns, "lightsequence")
        c_prop = _col(fns, "year", "proposed")
        c_rec = _col(fns, "year", "recommended")
        c_target = _col(fns, "target")
        c_comp = _col(fns, "companies")
        out = []
        for row in reader:
            if not (row.get(c_name) or "").strip():
                continue
            out.append(
                Therapeutic(
                    name=row[c_name].strip(),
                    fmt=row[c_fmt].strip(),
                    ch1_isotype=row[c_iso].strip(),
                    highest_trial=row[c_trial].strip(),
                    status=row[c_status].strip(),
                    heavy=row[c_heavy].strip(),
                    light=row[c_light].strip(),
                    year_proposed=_year(row[c_prop]),
                    year_recommended=_year(row[c_rec]),
                    target=row[c_target].strip(),
                    companies=row[c_comp].strip(),
                    raw=dict(row),
                )
            )
    return out
