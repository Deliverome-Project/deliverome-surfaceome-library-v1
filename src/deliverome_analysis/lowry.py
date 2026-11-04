"""Reusable analysis helpers for Lowry plate-reader standard curves.

The module parses SpectraMax UTF-16 text exports, extracts serial-dilution
standards and blanks, fits quadratic calibrations, tests lack of fit
against replicate error, and converts unknown absorbance values to concentration.

The lab-default assay is the Bio-Rad DC Protein Assay. Bio-Rad specifies a
0.2-1.5 mg/mL working range for its microplate procedure. Future standard
curves and diluted unknowns are restricted to that interval;
the plate reader's wider optical range does not expand the assay's validated
concentration range.
"""

from __future__ import annotations

import csv
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import numpy as np
from scipy import stats

BIO_RAD_DC_WORKING_RANGE_MG_ML = (0.2, 1.5)
"""Manufacturer-specified microplate working range for the Bio-Rad DC assay."""

BIO_RAD_DC_PROTOCOL_URL = (
    "https://www.bio-rad.com/sites/default/files/webroot/web/pdf/lsr/literature/Bulletin_6837.pdf"
)
"""Bio-Rad DC Protein Assay quick guide supporting the working-range constant."""


@dataclass(frozen=True)
class LowryCurve:
    """Concentrations, replicate absorbances, and pooled blank measurements."""

    concentrations_mg_ml: np.ndarray
    absorbances: np.ndarray
    blanks: np.ndarray

    @property
    def blank_mean(self) -> float:
        return float(self.blanks.mean())

    @property
    def blank_sd(self) -> float:
        return float(self.blanks.std(ddof=1))

    @property
    def corrected_absorbances(self) -> np.ndarray:
        return self.absorbances - self.blank_mean


@dataclass(frozen=True)
class CalibrationFit:
    """A fitted calibration and the statistics needed to evaluate it."""

    model: Literal["quadratic"]
    degree: int
    minimum_mg_ml: float
    maximum_mg_ml: float
    coefficients: np.ndarray
    coefficient_standard_errors: np.ndarray
    coefficient_ci95_low: np.ndarray
    coefficient_ci95_high: np.ndarray
    coefficient_p_values: np.ndarray
    blank_mean: float
    blank_sd: float
    n_levels: int
    n_measurements: int
    r_squared_replicates: float
    adjusted_r_squared_replicates: float
    r_squared_level_means: float
    rmse_absorbance: float
    lack_of_fit_f: float
    lack_of_fit_p_value: float

    def predict_corrected_absorbance(self, concentration_mg_ml: float | np.ndarray):
        """Predict blank-corrected absorbance from concentration."""
        return np.polyval(self.coefficients, concentration_mg_ml)

    def predict_raw_absorbance(self, concentration_mg_ml: float | np.ndarray):
        """Predict raw absorbance from concentration."""
        return self.predict_corrected_absorbance(concentration_mg_ml) + self.blank_mean

    def concentration_from_absorbance(
        self, raw_absorbance: float, dilution_factor: float = 1.0
    ) -> float:
        """Convert raw absorbance to original-sample concentration.

        For a quadratic calibration, the physically relevant real root inside
        the fitted concentration interval is selected. A value outside the
        calibrated response range raises ``ValueError`` instead of extrapolating.
        """
        if dilution_factor <= 0:
            raise ValueError("dilution_factor must be positive")
        target = raw_absorbance - self.blank_mean
        roots = np.roots([*self.coefficients[:-1], self.coefficients[-1] - target])
        valid = [
            float(root.real)
            for root in roots
            if abs(root.imag) < 1e-9
            and self.minimum_mg_ml - 1e-12 <= root.real <= self.maximum_mg_ml + 1e-12
        ]
        if len(valid) != 1:
            raise ValueError(
                "Absorbance does not map uniquely inside the fitted concentration range; "
                "dilute or concentrate the sample and read it again."
            )
        return valid[0] * dilution_factor

    def as_record(self) -> dict[str, float | int | str]:
        """Return a flat record suitable for a TSV output."""
        record: dict[str, float | int | str] = {
            "model": self.model,
            "range_mg_ml": f"{self.minimum_mg_ml:g}-{self.maximum_mg_ml:g}",
            "n_levels": self.n_levels,
            "n_measurements": self.n_measurements,
            "r_squared_replicates": self.r_squared_replicates,
            "adjusted_r_squared_replicates": self.adjusted_r_squared_replicates,
            "r_squared_level_means": self.r_squared_level_means,
            "rmse_absorbance": self.rmse_absorbance,
            "lack_of_fit_f": self.lack_of_fit_f,
            "lack_of_fit_p_value": self.lack_of_fit_p_value,
            "pooled_blank_mean": self.blank_mean,
            "pooled_blank_sd": self.blank_sd,
        }
        powers = range(self.degree, -1, -1)
        for power, value, standard_error, low, high, p_value in zip(
            powers,
            self.coefficients,
            self.coefficient_standard_errors,
            self.coefficient_ci95_low,
            self.coefficient_ci95_high,
            self.coefficient_p_values,
            strict=True,
        ):
            label = "intercept" if power == 0 else f"coefficient_x{power}"
            record[label] = value
            record[f"{label}_standard_error"] = standard_error
            record[f"{label}_ci95_low"] = low
            record[f"{label}_ci95_high"] = high
            record[f"{label}_p_value"] = p_value
        return record


def read_spectramax_txt(path: str | Path) -> dict[str, float]:
    """Read well values from a SpectraMax UTF-16 tab-delimited text export."""
    path = Path(path)
    with path.open(encoding="utf-16", newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    header_index = next(
        i for i, row in enumerate(rows) if len(row) > 2 and row[1].startswith("Temperature")
    )
    header = rows[header_index]
    values = next(row for row in rows[header_index + 1 :] if len(row) >= len(header) - 1)
    return {
        well: float(value)
        for well, value in zip(header[2:], values[2:], strict=False)
        if well and value
    }


def serial_dilution_concentrations(
    starting_mg_ml: float, number_of_levels: int, dilution_factor: float = 2.0
) -> np.ndarray:
    """Construct concentrations for a serial dilution, highest first."""
    if starting_mg_ml <= 0 or number_of_levels < 1 or dilution_factor <= 1:
        raise ValueError("Use a positive start, at least one level, and dilution_factor > 1")
    return starting_mg_ml / dilution_factor ** np.arange(number_of_levels)


def _read_wells(
    wells: dict[str, float], names: Iterable[str], description: str, hint: str = ""
) -> np.ndarray:
    """Look up well names, reporting every missing one at once.

    The bare ``KeyError`` that a direct lookup raises names only the first
    absent well and gives no clue that the plate simply uses a different
    layout, which is the usual cause.
    """
    names = list(names)
    missing = [name for name in names if name not in wells]
    if missing:
        # Show up to a full plate row, so a 9- or 12-well blank block or a
        # whole standard row is listed in full rather than clipped.
        shown = ", ".join(missing[:12])
        if len(missing) > 12:
            shown += f", ... ({len(missing)} total)"
        raise ValueError(
            f"{description} absent from the plate export: {shown}. "
            f"Wells present: {', '.join(sorted(wells)[:8])}, ...{hint}"
        )
    return np.array([wells[name] for name in names])


def extract_serial_dilution_curve(
    wells: dict[str, float],
    *,
    starting_mg_ml: float,
    number_of_levels: int,
    dilution_factor: float = 2.0,
    replicate_rows: str = "ABC",
    first_standard_column: int = 1,
    blank_columns: tuple[int, ...] | None = None,
    blank_wells: Sequence[str] | None = None,
) -> LowryCurve:
    """Extract a row-replicated serial dilution and pooled blanks from a plate.

    Blank position is **required** and has no default -- plates vary, and a
    wrong blank silently biases every concentration the curve reports rather
    than failing. Give exactly one of:

    ``blank_wells``
        An explicit list of well names -- any layout at all: a single row
        (``[f"D{i}" for i in range(1, 10)]``), one column, a scattered set, or
        a different number of blanks from the standards.
    ``blank_columns``
        Shorthand for a full block, read as ``replicate_rows`` x these columns;
        ``(10, 11, 12)`` over rows ``ABC`` means A10:C12.
    """
    if blank_columns is not None and blank_wells is not None:
        raise ValueError("Pass either blank_columns or blank_wells, not both")
    if blank_columns is None and blank_wells is None:
        raise ValueError(
            "Blank position is required: pass blank_wells=[...] for any well layout, "
            "or blank_columns=(...) for a replicate_rows x columns block. "
            "There is no default -- check the plate export's header."
        )
    if blank_wells is not None and len(blank_wells) == 0:
        raise ValueError("blank_wells must name at least one well")

    concentrations = serial_dilution_concentrations(
        starting_mg_ml, number_of_levels, dilution_factor
    )
    standard_columns = range(first_standard_column, first_standard_column + number_of_levels)
    absorbances = np.stack(
        [
            _read_wells(
                wells,
                (f"{row}{column}" for row in replicate_rows),
                f"Standard wells for column {column}",
                hint=" Check replicate_rows and first_standard_column.",
            )
            for column in standard_columns
        ]
    )

    if blank_wells is not None:
        names = list(blank_wells)
        hint = ""
    else:
        names = [f"{row}{column}" for column in blank_columns for row in replicate_rows]
        hint = " Pass blank_wells=[...] if the blanks are not a rows x columns block."
    blanks = _read_wells(wells, names, "Blank wells", hint=hint)

    return LowryCurve(concentrations, absorbances, blanks)


def load_serial_dilution_curve(path: str | Path, **kwargs) -> LowryCurve:
    """Convenience wrapper combining SpectraMax parsing and curve extraction."""
    return extract_serial_dilution_curve(read_spectramax_txt(path), **kwargs)


def fit_lowry_curve(
    curve: LowryCurve,
    *,
    minimum_mg_ml: float = BIO_RAD_DC_WORKING_RANGE_MG_ML[0],
    maximum_mg_ml: float = BIO_RAD_DC_WORKING_RANGE_MG_ML[1],
) -> CalibrationFit:
    """Fit the lab-default Bio-Rad DC quadratic calibration.

    The selection envelope must remain within Bio-Rad's validated 0.2-1.5 mg/mL
    microplate working range. The returned calibration bounds are tightened to
    the lowest and highest standards actually present, preventing extrapolation
    into an unmeasured part of that envelope. At least three distinct in-range
    standard levels are required. With exactly three levels, the quadratic is
    saturated at the level means, so lack-of-fit statistics are returned as
    ``nan``.
    """
    validated_minimum, validated_maximum = BIO_RAD_DC_WORKING_RANGE_MG_ML
    if minimum_mg_ml < validated_minimum or maximum_mg_ml > validated_maximum:
        raise ValueError(
            "Bio-Rad DC calibrations must stay within the validated "
            f"{validated_minimum:g}-{validated_maximum:g} mg/mL range"
        )
    if minimum_mg_ml >= maximum_mg_ml:
        raise ValueError("minimum_mg_ml must be less than maximum_mg_ml")

    model: Literal["quadratic"] = "quadratic"
    degree = 2
    included = (curve.concentrations_mg_ml >= minimum_mg_ml) & (
        curve.concentrations_mg_ml <= maximum_mg_ml
    )
    levels = curve.concentrations_mg_ml[included]
    values = curve.corrected_absorbances[included]
    parameter_count = degree + 1
    if len(np.unique(levels)) < parameter_count:
        raise ValueError(
            "A quadratic calibration needs at least three distinct standard "
            "levels inside 0.2-1.5 mg/mL"
        )

    x = np.repeat(levels, values.shape[1])
    y = values.reshape(-1)
    design = np.vander(x, parameter_count)
    coefficients, _, _, _ = np.linalg.lstsq(design, y, rcond=None)
    predictions = design @ coefficients
    residuals = y - predictions
    n = len(y)
    residual_df = n - parameter_count
    residual_ss = float(np.square(residuals).sum())
    residual_variance = residual_ss / residual_df
    covariance = residual_variance * np.linalg.inv(design.T @ design)
    standard_errors = np.sqrt(np.diag(covariance))
    t_values = coefficients / standard_errors
    coefficient_p_values = 2 * stats.t.sf(np.abs(t_values), residual_df)
    t_critical = stats.t.ppf(0.975, residual_df)

    total_ss = float(np.square(y - y.mean()).sum())
    r_squared = 1 - residual_ss / total_ss
    level_means = values.mean(axis=1)
    level_predictions = np.polyval(coefficients, levels)
    level_total_ss = float(np.square(level_means - level_means.mean()).sum())
    level_r_squared = 1 - float(np.square(level_means - level_predictions).sum()) / level_total_ss

    replicate_count = values.shape[1]
    pure_error_ss = float(np.square(values - level_means[:, None]).sum())
    lack_of_fit_ss = replicate_count * float(np.square(level_means - level_predictions).sum())
    lack_df = len(levels) - parameter_count
    pure_error_df = len(levels) * (replicate_count - 1)
    if lack_df > 0:
        lack_f = (lack_of_fit_ss / lack_df) / (pure_error_ss / pure_error_df)
        lack_p_value = stats.f.sf(lack_f, lack_df, pure_error_df)
    else:
        lack_f = float("nan")
        lack_p_value = float("nan")

    return CalibrationFit(
        model=model,
        degree=degree,
        minimum_mg_ml=float(levels.min()),
        maximum_mg_ml=float(levels.max()),
        coefficients=coefficients,
        coefficient_standard_errors=standard_errors,
        coefficient_ci95_low=coefficients - t_critical * standard_errors,
        coefficient_ci95_high=coefficients + t_critical * standard_errors,
        coefficient_p_values=coefficient_p_values,
        blank_mean=curve.blank_mean,
        blank_sd=curve.blank_sd,
        n_levels=len(levels),
        n_measurements=n,
        r_squared_replicates=r_squared,
        adjusted_r_squared_replicates=1 - (1 - r_squared) * (n - 1) / residual_df,
        r_squared_level_means=level_r_squared,
        rmse_absorbance=np.sqrt(residual_variance),
        lack_of_fit_f=lack_f,
        lack_of_fit_p_value=lack_p_value,
    )
