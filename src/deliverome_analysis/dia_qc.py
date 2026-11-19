"""QC metrics and figures for Spectronaut DIA data, before and after limpa.

Two stages, matching where problems show up:

* **Pre-limpa** (from the Spectronaut Normal Report, filtered exactly as limpa filters it:
  ``EG.Qvalue`` / ``PG.Qvalue`` <= cutoff, ``EG.IsImputed`` dropped, zero intensities missing):
  IDs per run, intensity distributions, data completeness, missingness vs intensity,
  replicate CVs, retention-time error and gradient coverage; and, when the export carries
  them, peak width, points per peak, signal-to-noise and the MS1 mass correction Spectronaut
  applied (from ``FG.PrecMz`` / ``FG.PrecMzCalibrated``, already in the factory report).
* **Post-limpa** (from a :class:`deliverome_analysis.limpa.LimpaResult`): the fitted detection
  probability curve against observed missingness, protein standard error vs abundance,
  run-run correlation (PCA from 4 runs), and agreement with Spectronaut's own PG.Quantity.

Usage
-----
    from deliverome_analysis.dia_qc import dia_qc_report
    from deliverome_analysis.limpa import run_limpa_spectronaut

    trim_spectronaut_report(report, "trimmed")          # once; QC then reads only this
    res = run_limpa_spectronaut(matrix_dir="trimmed")
    qc = dia_qc_report(matrix_dir="trimmed", output_dir="figures", prefix="s9e1_qc",
                       limpa_result=res)
    qc["run_summary"]        # per-run table: IDs, median intensity, RT error, ...
    qc["figures"]            # name -> saved PNG path

Individual ``plot_*`` functions return a Figure and assume ``setup_plotting_style()`` has
already been called; ``dia_qc_report`` calls it for you.

CVs are computed only **within replicate groups** (``R.Condition`` with >= 2 runs, or a
``group`` column you pass) -- a CV across different conditions measures biology, not
reproducibility. Contaminant tracking is deliberately absent: it needs a curated contaminant
accession list, which should come from a source file, not be typed into code.
"""

from __future__ import annotations

import os
import re
from pathlib import Path

import numpy as np
import pandas as pd

from deliverome_analysis.spectronaut_trim import PRECURSOR_QC_COLUMNS

_REPORT_COLUMNS = {
    "R.FileName": "run",
    "R.Condition": "condition",
    "PG.ProteinGroups": "protein",
    "PEP.StrippedSequence": "peptide",
    "EG.ModifiedSequence": "modseq",
    "FG.Charge": "charge",
    "EG.TotalQuantity (Settings)": "intensity",
    "PG.Quantity": "pg_quantity",
    "EG.ApexRT": "apex_rt",
    "EG.RTPredicted": "rt_predicted",
    "EG.Qvalue": "eg_q",
    "PG.Qvalue": "pg_q",
    "EG.IsImputed": "imputed",
    "PEP.NrOfMissedCleavages": "missed_cleavages",
    "PG.Genes": "gene",
    **PRECURSOR_QC_COLUMNS,
}
# Peak-shape metrics plotted per run when present: column -> axis label, log scale.
PEAK_METRICS = {
    "fwhm": ("Peak FWHM (min)", False),
    "peak_width": ("Peak width, start to end (min)", False),
    "points_per_peak": ("Data points per peak (MS2)", False),
    "points_per_peak_ms1": ("Data points per peak (MS1)", False),
    "signal_to_noise": ("Signal-to-noise", True),
    "ms1_calibration_ppm": ("MS1 mass correction applied (ppm)", False),
    "irt_residual": ("Empirical - predicted iRT", False),
    "tailing_factor": ("Mean tailing factor at FWHM", False),
    "cscore": ("Identification score (Cscore)", False),
    "ms1_mass_error_ppm": ("Observed MS1 mass error (ppm)", False),
    "peak_gauss_score": ("Peak Gaussian-shape score", False),
    "fragment_coelution": ("Fragment co-elution", False),
    "ion_ladder_score": ("Ion ladder (b/y) score", False),
}
_CENTRED = {"ms1_calibration_ppm", "irt_residual", "ms1_mass_error_ppm"}
_REQUIRED = [
    "R.FileName",
    "PG.ProteinGroups",
    "EG.ModifiedSequence",
    "FG.Charge",
    "EG.TotalQuantity (Settings)",
    "EG.Qvalue",
    "PG.Qvalue",
]


# ----------------------------------------------------------------------------------------------
# loading
# ----------------------------------------------------------------------------------------------


def load_spectronaut_precursors(
    report: str | os.PathLike, q_cutoff: float = 0.01, protein_id: str = "PG.ProteinGroups"
) -> pd.DataFrame:
    """Precursor x run rows from a Spectronaut Normal Report, filtered like limpa's reader.

    Returns columns: run, condition, protein, peptide, precursor, intensity, log2_intensity,
    pg_quantity, apex_rt, rt_predicted, rt_error (those the report carries). ``precursor``
    is ``<EG.ModifiedSequence>.<FG.Charge>``, the same key limpa uses.
    """
    header = pd.read_csv(report, sep="\t", nrows=0).columns
    cols = dict(_REPORT_COLUMNS)
    if protein_id != "PG.ProteinGroups":
        cols.pop("PG.ProteinGroups")
        cols[protein_id] = "protein"
    required = [protein_id if c == "PG.ProteinGroups" else c for c in _REQUIRED]
    missing = [c for c in required if c not in header]
    if missing:
        raise ValueError(f"Report is missing required column(s): {missing}")
    use = [c for c in cols if c in header]
    df = pd.read_csv(report, sep="\t", usecols=use).rename(columns=cols)

    keep = (df["eg_q"] <= q_cutoff) & (df["pg_q"] <= q_cutoff) & (df["intensity"] > 0)
    if "imputed" in df:
        keep &= ~df["imputed"].astype(bool)
    df = df.loc[keep].drop(columns=[c for c in ("eg_q", "pg_q", "imputed") if c in df])
    df["precursor"] = df["modseq"] + "." + df["charge"].astype(str)
    df["log2_intensity"] = np.log2(df["intensity"])
    if {"prec_mz", "prec_mz_calibrated"} <= set(df.columns):
        df["ms1_calibration_ppm"] = (df["prec_mz_calibrated"] - df["prec_mz"]) / df["prec_mz"] * 1e6
    df = df.drop(columns=[c for c in ("prec_mz", "prec_mz_calibrated") if c in df])
    if "condition" not in df:
        df["condition"] = "all"
    from deliverome_analysis.contaminants import classify

    df["protein_class"] = classify(df["protein"]).to_numpy()
    return _derive(df.drop(columns=["modseq", "charge"]).reset_index(drop=True))


def load_trimmed_precursors(matrix_dir: str | os.PathLike) -> pd.DataFrame:
    """The same table as :func:`load_spectronaut_precursors`, from a trimmed directory.

    Reads ``precursor_qc.parquet`` written by
    :func:`deliverome_analysis.spectronaut_trim.trim_spectronaut_report`, which already applied
    the q-value / imputation filters, so the full report is never reopened.
    """
    d = Path(matrix_dir)
    path = d / "precursor_qc.parquet"
    if not path.exists():
        raise FileNotFoundError(f"{path} not found; re-run trim_spectronaut_report on the report")
    df = pd.read_parquet(path)
    ann = pd.read_parquet(d / "precursors.parquet")
    if (d / "excluded_precursors.parquet").exists():  # contaminant/iRT rows kept out of limpa
        ann = pd.concat([ann, pd.read_parquet(d / "excluded_precursors.parquet")])
    ann = ann.set_index("precursor")
    df["protein"] = df["precursor"].map(ann["protein_id"])
    if "PEP.StrippedSequence" in ann:
        df["peptide"] = df["precursor"].map(ann["PEP.StrippedSequence"])
    if "PG.Genes" in ann:
        df["gene"] = df["precursor"].map(ann["PG.Genes"])
    if "PG.ProteinAccessions" in ann:  # stable IDs for tier and name lookups
        df["accessions"] = df["precursor"].map(ann["PG.ProteinAccessions"])
    if "protein_class" in ann:
        df["protein_class"] = df["precursor"].map(ann["protein_class"])
    samples = pd.read_csv(d / "samples.tsv", sep="\t", dtype=str).set_index("run")
    cond = samples["R.Condition"] if "R.Condition" in samples else pd.Series("all", samples.index)
    df["condition"] = df["run"].map(cond).fillna("all")
    df["intensity"] = np.exp2(df["log2_intensity"])
    if "PEP.NrOfMissedCleavages" in ann:
        mc = pd.to_numeric(ann["PEP.NrOfMissedCleavages"], errors="coerce")
        df["missed_cleavages"] = df["precursor"].map(mc)
    return _derive(df)


def _derive(df: pd.DataFrame) -> pd.DataFrame:
    if {"apex_rt", "rt_predicted"} <= set(df.columns):
        df["rt_error"] = df["apex_rt"] - df["rt_predicted"]
    if {"irt_empirical", "irt_predicted"} <= set(df.columns):
        # Per-run iRT recalibration residual. Without iRT spike-ins the scale is anchored on
        # endogenous library peptides, so this measures RT reproducibility, not drift vs. a standard.
        df["irt_residual"] = df["irt_empirical"] - df["irt_predicted"]
    if "peak_width" not in df and {"start_rt", "end_rt"} <= set(df.columns):
        df["peak_width"] = df["end_rt"] - df["start_rt"]
    return df


def run_order(prec: pd.DataFrame) -> list[str]:
    """Runs sorted by condition, then name (stable, readable axis order)."""
    runs = prec[["run", "condition"]].drop_duplicates()
    return runs.sort_values(["condition", "run"])["run"].tolist()


def _common_prefix(strings: list[str]) -> str:
    """The leading characters every string shares (run names, not paths)."""
    first, last = min(strings), max(strings)
    n = next((i for i, (a, b) in enumerate(zip(first, last, strict=False)) if a != b), len(first))
    return first[:n]


def short_run_names(runs: list[str]) -> dict[str, str]:
    """Strip the prefix/suffix every run name shares (cut at _ - . boundaries) for axis labels."""
    if len(runs) < 2:
        return {r: r for r in runs}
    pre = _common_prefix(runs)
    suf = _common_prefix([r[::-1] for r in runs])[::-1]
    pre = pre[: max((m.end() for m in re.finditer(r"[_\-. ]", pre)), default=0)]
    cut = [m.start() for m in re.finditer(r"[_\-. ]", suf)]
    suf = suf[cut[0] :] if cut else ""
    out = {r: r[len(pre) : len(r) - len(suf)] or r for r in runs}
    return out if len(set(out.values())) == len(runs) else {r: r for r in runs}


def run_labels(runs: list[str], condition: dict | pd.Series | None = None) -> dict[str, str]:
    """Axis labels: :func:`short_run_names`, plus ``(condition)`` when runs span conditions."""
    names = short_run_names(runs)
    if condition is None:
        return names
    cond = {r: str(condition[r]) for r in runs}
    if len(set(cond.values())) < 2:
        return names
    return {r: f"{names[r]} ({cond[r]})" for r in runs}


# ----------------------------------------------------------------------------------------------
# pre-limpa metrics
# ----------------------------------------------------------------------------------------------


def run_summary(prec: pd.DataFrame) -> pd.DataFrame:
    """Per-run IDs, intensity and RT-error summary, in :func:`run_order`."""
    g = prec.groupby("run")
    out = pd.DataFrame(
        {
            "condition": g["condition"].first(),
            "precursors": g["precursor"].nunique(),
            "peptides": g["peptide"].nunique() if "peptide" in prec else np.nan,
            "protein_groups": g["protein"].nunique(),
            "median_log2_intensity": g["log2_intensity"].median(),
        }
    )
    if "rt_error" in prec:
        out["median_abs_rt_error_min"] = g["rt_error"].apply(lambda s: s.abs().median())
    for col in PEAK_METRICS:
        if col in prec:
            out[f"median_{col}"] = g[col].median()
    if "protein_class" in prec:
        tot = g["intensity"].sum()
        for c in ("contaminant", "watchlist", "irt"):
            out[f"{c}_signal_pct"] = (
                (100 * prec[prec["protein_class"] == c].groupby("run")["intensity"].sum() / tot)
                .reindex(out.index)
                .fillna(0.0)
            )
    if "missed_cleavages" in prec:
        out["frac_missed_cleavage"] = g["missed_cleavages"].apply(lambda s: (s > 0).mean())
    n_prec = prec["precursor"].nunique()
    out["precursor_completeness"] = out["precursors"] / n_prec
    return out.loc[run_order(prec)]


def precursor_completeness(prec: pd.DataFrame) -> pd.DataFrame:
    """Per precursor: runs observed, fraction of runs observed, mean observed log2 intensity."""
    n_runs = prec["run"].nunique()
    g = prec.groupby("precursor")
    out = pd.DataFrame({"n_runs": g["run"].nunique(), "mean_log2": g["log2_intensity"].mean()})
    out["frac_observed"] = out["n_runs"] / n_runs
    return out


def replicate_groups(prec: pd.DataFrame, groups: pd.Series | None = None) -> dict[str, list[str]]:
    """Replicate groups with >= 2 runs: ``groups`` (run -> label) or Spectronaut R.Condition."""
    labels = groups if groups is not None else prec.groupby("run")["condition"].first()
    out = {str(k): sorted(v.index) for k, v in labels.groupby(labels)}
    return {k: v for k, v in out.items() if len(v) >= 2}


def precursor_cv(prec: pd.DataFrame, groups: pd.Series | None = None) -> pd.DataFrame:
    """Linear-scale CV of each precursor seen in every run of its replicate group.

    Returns columns: group, precursor, cv (fraction), mean_log2. Empty if no group has >= 2 runs.
    """
    rows = []
    for name, runs in replicate_groups(prec, groups).items():
        wide = prec[prec["run"].isin(runs)].pivot_table(
            index="precursor", columns="run", values="intensity", aggfunc="first"
        )
        wide = wide.dropna()
        if wide.empty:
            continue
        cv = wide.std(axis=1, ddof=1) / wide.mean(axis=1)
        rows.append(
            pd.DataFrame(
                {
                    "group": name,
                    "precursor": wide.index,
                    "cv": cv.values,
                    "mean_log2": np.log2(wide).mean(axis=1).values,
                }
            )
        )
    cols = ["group", "precursor", "cv", "mean_log2"]
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(columns=cols)


def rt_consistency(
    prec: pd.DataFrame, *, bins: int = 50, k: float = 5.0, min_flag_min: float = 0.1
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Run-to-run retention-time differences per precursor, with outliers flagged.

    For each precursor and run: ``delta_min`` = apex RT minus the median apex RT of the same
    precursor across runs (with 2 runs: minus the other run). Each run's smooth shift (median delta within ``bins``
    quantile bins of RT, interpolated) is the ``trend_min``: LC drift or a gradient offset,
    which affects every peptide alike. ``resid_min`` = delta - trend is what is left for the
    individual peptide; it is flagged when ``|resid| > max(k * robust SD of the run, min_flag_min)``.
    A lone outlier usually means a different peak was picked (interference or a wrong
    match), not an LC problem.

    Returns (per precursor x run table, per run summary). Needs >= 2 runs and ``apex_rt``.
    """
    wide = prec.pivot_table(index="precursor", columns="run", values="apex_rt", aggfunc="first")
    runs = list(wide.columns)
    rows, summ = [], []
    # >= 3 runs: the median over all runs is robust to one run picking a different peak, so
    # only that run is flagged. 2 runs: there is no majority, both runs see the same |delta|.
    ref_all = wide.median(axis=1) if len(runs) >= 3 else None
    for r in runs:
        ref = ref_all if ref_all is not None else wide.drop(columns=r).median(axis=1)
        d = pd.DataFrame({"apex_rt": wide[r], "ref_rt": ref}).dropna()
        if len(d) < 30:
            continue
        nb = max(3, min(bins, len(d) // 20))
        d["delta_min"] = d["apex_rt"] - d["ref_rt"]
        q = np.unique(np.quantile(d["ref_rt"], np.linspace(0, 1, nb + 1)))
        idx = np.clip(np.searchsorted(q, d["ref_rt"], side="right") - 1, 0, len(q) - 2)
        centers = d.groupby(idx)["ref_rt"].median().to_numpy()
        medians = d.groupby(idx)["delta_min"].median().to_numpy()
        d["trend_min"] = np.interp(d["ref_rt"], centers, medians)
        d["resid_min"] = d["delta_min"] - d["trend_min"]
        mad = 1.4826 * np.median(np.abs(d["resid_min"] - np.median(d["resid_min"])))
        cut = max(k * mad, min_flag_min)
        d["flag"] = d["resid_min"].abs() > cut
        d["run"] = r
        rows.append(d.reset_index())
        summ.append(
            {
                "run": r,
                "precursors": len(d),
                "median_shift_s": 60 * d["delta_min"].median(),
                "drift_range_s": 60 * (medians.max() - medians.min()),
                "residual_sd_s": 60 * mad,
                "flag_cutoff_s": 60 * cut,
                "flagged": int(d["flag"].sum()),
                "flagged_pct": 100 * d["flag"].mean(),
            }
        )
    out = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    if not out.empty:
        info = prec.drop_duplicates("precursor").set_index("precursor")
        for c in ("protein", "peptide"):
            if c in info:
                out[c] = out["precursor"].map(info[c])
    return out, pd.DataFrame(summ)


def class_signal(prec: pd.DataFrame) -> pd.DataFrame:
    """Per run: share (%) of summed precursor intensity and precursor counts by protein class.

    Classes come from ``deliverome_analysis.contaminants`` (human, watchlist, mixed,
    contaminant, irt). Needs the ``protein_class`` column the trimmer writes.
    """
    from deliverome_analysis.contaminants import CLASSES

    g = prec.groupby(["run", "protein_class"])
    inten = g["intensity"].sum().unstack(fill_value=0).reindex(columns=CLASSES, fill_value=0)
    n = g["precursor"].nunique().unstack(fill_value=0).reindex(columns=CLASSES, fill_value=0)
    pct = 100 * inten.div(inten.sum(axis=1), axis=0)
    out = pct.add_suffix("_signal_pct").join(n.add_suffix("_precursors"))
    return out.loc[[r for r in run_order(prec) if r in out.index]]


def top_contaminants(prec: pd.DataFrame, n: int = 30) -> pd.DataFrame:
    """Contaminant, mixed and watchlist protein groups ranked by median share of run signal."""
    p = prec[prec["protein_class"].isin(["contaminant", "mixed", "watchlist"])]
    if p.empty:
        return pd.DataFrame(columns=["protein", "protein_class", "median_signal_pct", "runs"])
    total = prec.groupby("run")["intensity"].sum()
    s = p.groupby(["protein", "protein_class", "run"])["intensity"].sum().reset_index()
    s["pct"] = 100 * s["intensity"] / s["run"].map(total)
    out = s.groupby(["protein", "protein_class"]).agg(
        median_signal_pct=("pct", "median"), runs=("run", "nunique")
    )
    return out.reset_index().sort_values("median_signal_pct", ascending=False).head(n)


# ----------------------------------------------------------------------------------------------
# surface (membrane) proteins by accessible-surfaceome tier
# ----------------------------------------------------------------------------------------------

# Display classes, most to least confident. A protein group takes the best tier of any member.
TIER_CLASSES = ["canonical", "likely", "low / uncertain", "not surface"]
SURFACE_CLASSES = TIER_CLASSES[:3]
_TIER_RANK = {"canonical": 0, "likely": 1, "low": 2, "uncertain": 2}


def surface_annotation(
    prec: pd.DataFrame,
    tiers: dict[str, str] | None = None,
    symbols: dict[str, str] | None = None,
) -> pd.DataFrame:
    """Per protein group: ``tier`` (one of :data:`TIER_CLASSES`), ``symbol``, ``accessions``.

    Joins on base UniProt accession (``PG.ProteinAccessions`` when the export has it, else the
    protein id). ``tiers`` maps base accession -> catalog ``deep_dive_tier`` and ``symbols``
    base accession -> gene symbol (display only); both default to the pinned
    accessible-surfaceome catalog (:mod:`deliverome_analysis.surfaceome`, needs S3).
    """
    from deliverome_analysis import surfaceome

    if tiers is None:
        cat = surfaceome.load_catalog()
        tiers = dict(zip(cat["uniprot_accession"], cat["tier"], strict=True))
    if symbols is None:
        symbols = surfaceome.catalog_symbols()
    groups = prec.drop_duplicates("protein")
    accs_of = groups["accessions"] if "accessions" in groups else groups["protein"]
    rows = []
    for pid, accs in zip(groups["protein"], accs_of, strict=True):
        bases = list(
            dict.fromkeys(
                surfaceome.base_accession(a.strip()) for a in str(accs).split(";") if a.strip()
            )
        )
        ranks = [_TIER_RANK[tiers[b]] for b in bases if tiers.get(b) in _TIER_RANK]
        syms = dict.fromkeys(symbols[b] for b in bases if symbols.get(b))
        rows.append(
            (
                pid,
                TIER_CLASSES[min(ranks)] if ranks else "not surface",
                ";".join(syms),
                ";".join(bases),
            )
        )
    return pd.DataFrame(rows, columns=["protein", "tier", "symbol", "accessions"]).set_index(
        "protein"
    )


def _with_tier(prec: pd.DataFrame, ann: pd.DataFrame) -> pd.DataFrame:
    t = prec["protein"].map(ann["tier"]).fillna("not surface")
    return prec.assign(tier=pd.Categorical(t, TIER_CLASSES))


def surface_by_tier(prec: pd.DataFrame, ann: pd.DataFrame) -> dict[str, pd.DataFrame]:
    """Per run x tier: identifications, signal shares and digestion.

    Keys: ``peptides`` and ``protein_groups`` (counts), ``signal_ms2`` and, when the export
    has ``FG.MS1Quantity``, ``signal_ms1`` (% of the run's summed precursor signal),
    ``missed`` (% of precursors with at least one missed cleavage), and the medians
    ``points_per_peak`` and ``fwhm`` (peak width, minutes) when the export has them. Rows follow
    :func:`run_order`, columns :data:`TIER_CLASSES`.
    """
    p = _with_tier(prec, ann)
    order = run_order(prec)

    def by(col, fn):
        return p.groupby(["run", "tier"], observed=False)[col].agg(fn).unstack().reindex(order)

    out = {"peptides": by("peptide" if "peptide" in p else "precursor", "nunique")}
    out["protein_groups"] = by("protein", "nunique")
    ms2 = "ms2_quantity" if "ms2_quantity" in p and p["ms2_quantity"].notna().any() else "intensity"
    for key, col in (("signal_ms2", ms2), ("signal_ms1", "ms1_quantity")):
        if col in p and p[col].notna().any():
            tot = p.groupby("run")[col].sum().reindex(order)
            out[key] = 100 * by(col, "sum").div(tot, axis=0)
    if "missed_cleavages" in p and p["missed_cleavages"].notna().any():
        p["mc_any"] = (p["missed_cleavages"] > 0).astype(float).where(p["missed_cleavages"].notna())
        out["missed"] = 100 * by("mc_any", "mean")
        out["missed_all"] = 100 * p.groupby("run")["mc_any"].mean().reindex(order)
    for key in ("points_per_peak", "fwhm"):
        if key in p and p[key].notna().any():
            out[key] = by(key, "median")
    return out


def top_surface_proteins(prec: pd.DataFrame, ann: pd.DataFrame, n: int = 25) -> pd.DataFrame:
    """Per run and surface tier, the ``n`` groups with the most signal: run, tier, protein, share."""
    p = _with_tier(prec, ann)
    p = p[p["tier"].isin(SURFACE_CLASSES)]
    tot = prec.groupby("run")["intensity"].sum()
    g = p.groupby(["run", "tier", "protein"], observed=True)["intensity"].sum().reset_index()
    # reindex, not .map: with no surface group g is empty, and pandas 3 keeps the str dtype of
    # an empty `run` column through .map, so the division would fail on str / float.
    g["share_pct"] = 100 * g["intensity"] / tot.reindex(g["run"]).to_numpy()
    g = g.sort_values("intensity", ascending=False).groupby(["run", "tier"], observed=True).head(n)
    return g.drop(columns="intensity").reset_index(drop=True)


def plot_class_signal(cs: pd.DataFrame, prec: pd.DataFrame):
    """Stacked bars: % of each run's precursor signal from contaminant, mixed, watchlist, iRT."""
    plt = _plt()
    from deliverome_analysis.plotting import COLORS

    names = run_labels(list(cs.index), prec.groupby("run")["condition"].first())
    cols = [
        ("contaminant", COLORS["primary"]),
        ("mixed", COLORS["accent"]),
        ("watchlist", COLORS["quaternary"]),
        ("irt", COLORS["secondary"]),
    ]
    fig, ax = plt.subplots(figsize=(_run_axis_width(len(cs)), 6))
    bottom = np.zeros(len(cs))
    for c, color in cols:
        v = cs[f"{c}_signal_pct"].to_numpy()
        ax.bar(range(len(cs)), v, bottom=bottom, color=color, width=0.7, label=c)
        bottom += v
    ax.set_xticks(range(len(cs)), [names[r] for r in cs.index], rotation=30, ha="right")
    ax.set_ylabel("Share of precursor signal (%)")
    ax.legend(frameon=False, fontsize="x-small")
    return fig


def charge_state_rt_spread(prec: pd.DataFrame, *, min_flag_min: float = 0.1) -> pd.DataFrame:
    """Modified peptides seen at >= 2 charge states in a run: apex RT spread across charges.

    The charge states of one peptide co-elute; a spread above ``min_flag_min`` (6 s by
    default) suggests one of them was matched to the wrong peak.
    """
    # group by *modified* sequence: differently modified forms (oxidised Met, N-terminal
    # acetylation) of one stripped peptide genuinely elute apart
    d = prec.dropna(subset=["apex_rt"]).assign(
        modified_peptide=lambda x: x["precursor"].str.rsplit(".", n=1).str[0]
    )
    g = d.groupby(["run", "modified_peptide"])["apex_rt"]
    out = pd.DataFrame({"n_charges": g.size(), "spread_min": g.max() - g.min()})
    out = out[out["n_charges"] >= 2].reset_index()
    out["flag"] = out["spread_min"] > min_flag_min
    return out


def plot_rt_consistency(rtc: pd.DataFrame, prec: pd.DataFrame):
    """Per run: RT difference vs the other runs across the gradient (trend + 5-95% band),
    flagged precursors marked."""
    plt = _plt()
    order = [r for r in run_order(prec) if r in set(rtc["run"])]
    names = run_labels(order, prec.groupby("run")["condition"].first())
    colors = dict(zip(order, _palette(len(order)), strict=True))
    fig, ax = plt.subplots(figsize=(10, 6))
    for r in order:
        d = rtc[rtc["run"] == r].sort_values("ref_rt")
        b = pd.qcut(d["ref_rt"], max(3, min(40, len(d) // 20)), duplicates="drop")
        g = d.groupby(b, observed=True)
        x = g["ref_rt"].median()
        ax.fill_between(
            x,
            60 * g["delta_min"].quantile(0.05),
            60 * g["delta_min"].quantile(0.95),
            color=colors[r],
            alpha=0.15,
            lw=0,
        )
        ax.plot(x, 60 * g["delta_min"].median(), color=colors[r], lw=2.5, label=names[r])
        f = d[d["flag"]]
        ax.scatter(f["ref_rt"], 60 * f["delta_min"], s=6, color=colors[r], alpha=0.5)
    ax.axhline(0, color="black", lw=0.8)
    lim = 60 * np.nanpercentile(np.abs(rtc["delta_min"]), 99.5)
    ax.set_ylim(-lim, lim)
    ax.set_xlabel("Retention time in the other runs (min)")
    ax.set_ylabel("This run minus other runs (s)")
    ax.legend(frameon=False, fontsize="x-small")
    return fig


def conditional_detection(prec: pd.DataFrame, bins: int = 30, min_count: int = 50) -> pd.DataFrame:
    """Empirical detection curve: P(detected in another run | observed at log2 x in this run).

    For every precursor observed in run i at log2 intensity x, count how many of the other
    runs also detected it, then average within bins of x. Unlike "fraction of runs detected",
    this has no 1/n floor, so it is meaningful from 2 runs up. It estimates the quantity
    limpa's DPC models, with one bias: x is a noisy single-run value, so at the low end it
    reads slightly above the fitted curve. Returns columns log2_intensity, p_detected, n.
    """
    wide = prec.pivot_table(
        index="precursor", columns="run", values="log2_intensity", aggfunc="first"
    )
    obs = wide.notna().to_numpy()
    vals = wide.to_numpy()
    n_runs = obs.shape[1]
    if n_runs < 2:
        return pd.DataFrame(columns=["log2_intensity", "p_detected", "n"])
    others_detected = obs.sum(axis=1, keepdims=True) - obs  # per (precursor, run i)
    x = vals[obs]
    y = others_detected[obs] / (n_runs - 1)
    edges = np.linspace(np.nanmin(x), np.nanmax(x), bins + 1)
    idx = np.clip(np.digitize(x, edges) - 1, 0, bins - 1)
    d = pd.DataFrame({"bin": idx, "x": x, "y": y})
    b = d.groupby("bin").agg(
        log2_intensity=("x", "mean"), p_detected=("y", "mean"), n=("y", "size")
    )
    return b[b["n"] >= min_count].reset_index(drop=True)


# ----------------------------------------------------------------------------------------------
# plotting helpers
# ----------------------------------------------------------------------------------------------


def _plt():
    import matplotlib.pyplot as plt

    return plt


def _palette(n: int) -> list[str]:
    from deliverome_analysis.plotting import CATEGORICAL_PALETTE

    return [CATEGORICAL_PALETTE[i % len(CATEGORICAL_PALETTE)] for i in range(n)]


def _condition_colors(prec_or_summary: pd.DataFrame) -> dict[str, str]:
    conds = sorted(prec_or_summary["condition"].astype(str).unique())
    return dict(zip(conds, _palette(len(conds)), strict=True))


def _run_axis_width(n_runs: int) -> float:
    return max(6.0, 1.1 * n_runs + 3)


# ----------------------------------------------------------------------------------------------
# pre-limpa figures
# ----------------------------------------------------------------------------------------------


def plot_ids_per_run(summary: pd.DataFrame):
    """One panel per level (precursors, peptides, protein groups): IDs per run at 1% FDR."""
    plt = _plt()

    levels = [c for c in ("precursors", "peptides", "protein_groups") if summary[c].notna().all()]
    names = run_labels(list(summary.index), summary["condition"])
    cmap = _condition_colors(summary)
    colors = [cmap[str(c)] for c in summary["condition"]]
    x = np.arange(len(summary))
    fig, axes = plt.subplots(1, len(levels), figsize=(len(levels) * (0.7 * len(summary) + 3), 6))
    for ax, lvl in zip(np.atleast_1d(axes), levels, strict=True):
        ax.bar(x, summary[lvl], color=colors, width=0.7)
        ax.set_xticks(x, [names[r] for r in summary.index], rotation=30, ha="right")
        ax.set_ylabel(f"{lvl.replace('_', ' ').capitalize()} (1% FDR)")
        ax.ticklabel_format(axis="y", style="sci", scilimits=(-3, 4))
    fig.tight_layout()
    return fig


def plot_intensity_distributions(prec: pd.DataFrame):
    """Per-run boxplots of precursor log2 intensity (normalization / loading check)."""
    plt = _plt()
    order = run_order(prec)
    cmap = _condition_colors(prec)
    cond = prec.groupby("run")["condition"].first()
    names = run_labels(order, cond)
    data = [prec.loc[prec["run"] == r, "log2_intensity"].values for r in order]
    fig, ax = plt.subplots(figsize=(_run_axis_width(len(order)), 6))
    bp = ax.boxplot(data, patch_artist=True, showfliers=False, widths=0.6)
    for patch, r in zip(bp["boxes"], order, strict=True):
        patch.set_facecolor(cmap[str(cond[r])])
        patch.set_alpha(0.8)
    for med in bp["medians"]:
        med.set_color("black")
    ax.set_xticks(range(1, len(order) + 1), [names[r] for r in order], rotation=30, ha="right")
    ax.set_ylabel("Precursor log2 intensity")
    _condition_legend(ax, cmap)
    return fig


def _condition_legend(ax, cmap: dict[str, str]) -> None:
    if len(cmap) < 2:
        return
    from matplotlib.patches import Patch

    ax.legend(
        handles=[Patch(color=c, label=k) for k, c in cmap.items()], frameon=False, fontsize="small"
    )


def plot_completeness(compl: pd.DataFrame, n_runs: int):
    """Bars: number of precursors observed in exactly k of the n runs."""
    plt = _plt()
    from deliverome_analysis.plotting import COLORS

    counts = compl["n_runs"].value_counts().reindex(range(1, n_runs + 1), fill_value=0)
    fig, ax = plt.subplots(figsize=(max(6.0, 0.6 * n_runs + 4), 6))
    ax.bar(counts.index, counts.values, color=COLORS["secondary"])
    ax.set_xlabel("Runs in which precursor is observed")
    ax.set_ylabel("Precursors")
    ax.set_xticks(counts.index if n_runs <= 20 else ax.get_xticks())
    full = counts.get(n_runs, 0) / max(counts.sum(), 1)
    ax.text(
        0.02,
        0.97,
        f"{full:.0%} observed in all {n_runs} runs",
        transform=ax.transAxes,
        va="top",
        fontsize="small",
    )
    return fig


def plot_detection_reproducibility(det: pd.DataFrame):
    """P(precursor detected in another run | seen at log2 x in this run), binned.

    A pre-limpa QC view: "if I see this precursor at intensity x, will I see it again?".
    Not an estimate of limpa's DPC -- a single run's value includes noise, so low bins are
    dominated by precursors that dipped in one run but are strong elsewhere (see
    :func:`plot_dpc_fit` for limpa's own diagnostic).
    """
    plt = _plt()
    from deliverome_analysis.plotting import COLORS

    fig, ax = plt.subplots(figsize=(7.5, 6))
    ax.scatter(
        det["log2_intensity"],
        det["p_detected"],
        s=40 + 160 * det["n"] / det["n"].max(),
        color=COLORS["secondary"],
        alpha=0.85,
    )
    ax.set_xlabel("log2 intensity in one run")
    ax.set_ylabel("P(detected in another run)")
    ax.set_ylim(-0.02, 1.05)
    return fig


def plot_dpc_fit(result, bins: int = 20, min_count: int = 10):
    """limpa's DPC diagnostic: detected fraction vs limpa's estimated mean log2 intensity.

    Uses the precursors dpc() was fitted on, with x = limpa's complete-data mean estimate
    (which accounts for the runs where the precursor was missing), as in ``limpa::plotDPC``.
    Light points are individual precursors (jittered); large points are bin means.
    """
    plt = _plt()
    from deliverome_analysis.limpa import detection_probability
    from deliverome_analysis.plotting import COLORS

    pts = result.dpc_points.dropna(subset=["mu"])
    frac = pts["n_detected"] / pts["nsamples"]
    rng = np.random.default_rng(0)
    jitter = rng.uniform(-0.5, 0.5, len(pts)) / pts["nsamples"].iloc[0] * 0.6
    edges = np.linspace(pts["mu"].min(), pts["mu"].max(), bins + 1)
    idx = np.clip(np.digitize(pts["mu"], edges) - 1, 0, bins - 1)
    b = pd.DataFrame({"bin": idx, "mu": pts["mu"], "frac": frac}).groupby("bin")
    b = b.agg(mu=("mu", "mean"), frac=("frac", "mean"), n=("frac", "size"))
    b = b[b["n"] >= min_count]
    fig, ax = plt.subplots(figsize=(7.5, 6))
    ax.scatter(
        pts["mu"],
        frac + jitter,
        s=5,
        alpha=0.25,
        color=COLORS["neutral"],
        rasterized=True,
        label="precursors (jittered)",
    )
    ax.scatter(
        b["mu"],
        b["frac"],
        s=40 + 160 * b["n"] / b["n"].max(),
        color=COLORS["primary"],
        alpha=0.9,
        label="bin mean",
    )
    x = np.linspace(pts["mu"].min(), pts["mu"].max(), 200)
    b0, b1 = result.dpc
    ax.plot(
        x,
        detection_probability(x, b0, b1),
        color=COLORS["dark"],
        lw=2.5,
        label=f"limpa DPC ($\\beta_0$={b0:.2f}, $\\beta_1$={b1:.2f})",
    )
    ax.set_xlabel("limpa estimated mean log2 intensity")
    ax.set_ylabel("Fraction of runs detected")
    ax.set_ylim(-0.12, 1.12)
    ax.legend(frameon=False, fontsize="x-small", loc="lower right", markerscale=1)
    return fig


def plot_precursor_cv(cv: pd.DataFrame):
    """Histogram of replicate CVs (%), one line per replicate group, with medians."""
    plt = _plt()
    groups = sorted(cv["group"].unique())
    colors = dict(zip(groups, _palette(len(groups)), strict=True))
    fig, ax = plt.subplots(figsize=(7.5, 6))
    bins = np.linspace(0, 100, 51)
    lines = []
    for gname in groups:
        pct = 100 * cv.loc[cv["group"] == gname, "cv"].clip(upper=1.0)
        ax.hist(pct, bins=bins, histtype="step", lw=2.5, color=colors[gname])
        med = pct.median()
        ax.axvline(med, color=colors[gname], ls="--", lw=1.5)
        lines.append(f"{gname}: median {med:.1f}%, {(pct < 20).mean():.0%} < 20% (n={len(pct):,})")
    lines.append("CVs >= 100% are piled into the last bin")
    ax.set_xlabel("Precursor CV across replicate runs (%)")
    ax.set_ylabel("Precursors")
    ax.text(
        0.98,
        0.97,
        "\n".join(lines),
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize="x-small",
    )
    return fig


def plot_run_violins(prec: pd.DataFrame, column: str, ylabel: str, *, centred=False, log=False):
    """Per-run violins of one precursor-level metric, median marked, coloured by condition.

    ``centred`` puts zero mid-axis (for errors); ``log`` plots log10 of a positive metric. The
    axis is clipped to the central 98% so a few outliers do not flatten every violin.
    """
    plt = _plt()
    order = run_order(prec)
    cmap = _condition_colors(prec)
    cond = prec.groupby("run")["condition"].first()
    names = run_labels(order, cond)
    vals = prec[column].astype(float)
    if log:
        vals = np.log10(vals.where(vals > 0))
    if centred:
        hi = np.nanpercentile(np.abs(vals), 99)
        lo, ylim = -hi, (-hi, hi)
    else:
        lo, hi = np.nanpercentile(vals, [1, 99])
        pad = 0.05 * (hi - lo or 1)
        ylim = (lo - pad, hi + pad)
    # Medians come from all values; the density is estimated only inside the plotted range,
    # because heavy tails inflate the SD-based bandwidth and smear the core into a diamond.
    groups = [vals[(prec["run"] == r).to_numpy()].dropna() for r in order]
    data = [g[(g >= lo) & (g <= hi)].values for g in groups]
    fig, ax = plt.subplots(figsize=(_run_axis_width(len(order)), 6))
    vp = ax.violinplot(data, showextrema=False)
    for body, r in zip(vp["bodies"], order, strict=True):
        body.set_facecolor(cmap[str(cond[r])])
        body.set_alpha(0.7)
    meds = [g.median() for g in groups]
    ax.hlines(
        meds,
        np.arange(1, len(order) + 1) - 0.25,
        np.arange(1, len(order) + 1) + 0.25,
        color="black",
        lw=2,
    )
    ax.set_ylim(*ylim)
    if centred:
        ax.axhline(0, color="black", lw=0.8)
    ax.set_xticks(range(1, len(order) + 1), [names[r] for r in order], rotation=30, ha="right")
    ax.set_ylabel(f"log10 {ylabel.lower()}" if log else ylabel)
    return fig


def plot_rt_error(prec: pd.DataFrame):
    """Per-run distribution of apex RT minus Spectronaut's calibrated predicted RT (min)."""
    return plot_run_violins(prec, "rt_error", "Apex RT - predicted RT (min)", centred=True)


def plot_peak_metric(prec: pd.DataFrame, column: str):
    """Per-run violins of one :data:`PEAK_METRICS` column (peak width, points, S/N, ppm)."""
    label, log = PEAK_METRICS[column]
    return plot_run_violins(prec, column, label, centred=column in _CENTRED, log=log)


def plot_missed_cleavages(prec: pd.DataFrame):
    """Per-run share of identified precursors with 0, 1 and 2+ missed cleavages."""
    plt = _plt()
    order = run_order(prec)
    names = run_labels(order, prec.groupby("run")["condition"].first())
    mc = prec["missed_cleavages"].clip(upper=2).map({0: "0", 1: "1", 2: "2+"})
    frac = pd.crosstab(prec["run"], mc, normalize="index").reindex(order).fillna(0)
    frac = frac.reindex(columns=[c for c in ("0", "1", "2+") if c in frac.columns])
    fig, ax = plt.subplots(figsize=(_run_axis_width(len(order)), 6))
    bottom = np.zeros(len(order))
    for col, color in zip(frac.columns, _palette(len(frac.columns)), strict=True):
        ax.bar(range(len(order)), 100 * frac[col], bottom=bottom, color=color, width=0.7, label=col)
        bottom += 100 * frac[col].to_numpy()
    ax.set_xticks(range(len(order)), [names[r] for r in order], rotation=30, ha="right")
    ax.set_ylabel("Identified precursors (%)")
    ax.legend(title="Missed cleavages", frameon=False, fontsize="x-small", loc="lower right")
    return fig


def plot_gradient_coverage(prec: pd.DataFrame):
    """Precursor identifications across the gradient (apex RT), one line per run."""
    plt = _plt()
    order = run_order(prec)
    names = run_labels(order, prec.groupby("run")["condition"].first())
    colors = dict(zip(order, _palette(len(order)), strict=True))
    bins = np.linspace(prec["apex_rt"].min(), prec["apex_rt"].max(), 61)
    fig, ax = plt.subplots(figsize=(9, 6))
    for r in order:
        ax.hist(
            prec.loc[prec["run"] == r, "apex_rt"],
            bins=bins,
            histtype="step",
            lw=2,
            color=colors[r],
            label=names[r],
        )
    ax.set_xlabel("Apex retention time (min)")
    ax.set_ylabel("Precursors per bin")
    if len(order) <= 12:
        ax.legend(frameon=False, fontsize="x-small")
    return fig


# ----------------------------------------------------------------------------------------------
# post-limpa figures
# ----------------------------------------------------------------------------------------------


def _nprec_bucket(n: pd.Series) -> pd.Series:
    return pd.cut(n, [0, 1, 3, 9, np.inf], labels=["1", "2-3", "4-9", "10+"])


def plot_protein_se(result):
    """limpa protein standard error vs mean log2 expression, coloured by precursors per protein."""
    plt = _plt()
    mean = result.protein_log2.mean(axis=1)
    se = result.protein_se.mean(axis=1)
    bucket = _nprec_bucket(result.proteins["NPrec"].reindex(mean.index))
    cats = list(bucket.cat.categories)
    colors = dict(zip(cats, _palette(len(cats)), strict=True))
    fig, ax = plt.subplots(figsize=(7.5, 6))
    for c in cats:
        m = bucket == c
        ax.scatter(
            mean[m],
            se[m],
            s=6,
            alpha=0.5,
            color=colors[c],
            rasterized=True,
            label=f"{c} precursors (n={int(m.sum()):,})",
        )
    ax.set_xlabel("Protein log2 expression (limpa, mean over runs)")
    ax.set_ylabel("Standard error (mean over runs)")
    ax.legend(frameon=False, fontsize="x-small", markerscale=3)
    return fig


def plot_run_correlation(result, prec: pd.DataFrame | None = None):
    """Pearson correlation of limpa protein log2 expression between runs.

    With ``prec``, only proteins detected in every run are used: limpa fills a protein's
    missing runs from the DPC, and those model values (which can sit many log2 units from the
    observed run) would otherwise dominate the correlation. Values are on limpa's scale,
    whose between-protein spread is narrower than Spectronaut's, so expect lower r than
    Spectronaut's own run correlation for the same data; compare runs within a dataset.
    """
    plt = _plt()
    from matplotlib.colors import LinearSegmentedColormap

    from deliverome_analysis.plotting import get_sequential_palette

    samples = result.samples.assign(_cond=result.samples.get("R.Condition", "all"))
    order = samples.sort_values(["_cond"]).index.tolist()
    order = sorted(order, key=lambda r: (str(samples.loc[r, "_cond"]), r))
    mat = result.protein_log2[order]
    if prec is not None:
        seen = prec.groupby("protein")["run"].nunique()
        mat = mat.loc[mat.index.intersection(seen.index[seen == len(order)])]
    corr = mat.corr()
    names = run_labels(order, samples["_cond"])
    n = len(order)
    cmap = LinearSegmentedColormap.from_list("teal", get_sequential_palette("teal")[::-1])
    fig, ax = plt.subplots(figsize=(max(6.0, 0.8 * n + 3), max(5.0, 0.8 * n + 2)))
    lo = min(0.9, float(np.nanmin(corr.values)))
    im = ax.imshow(corr.values, cmap=cmap, vmin=lo, vmax=1)
    ax.set_xticks(range(n), [names[r] for r in order], rotation=45, ha="right")
    ax.set_yticks(range(n), [names[r] for r in order])
    if n <= 12:
        for i in range(n):
            for j in range(n):
                v = corr.values[i, j]
                ax.text(
                    j,
                    i,
                    f"{v:.3f}",
                    ha="center",
                    va="center",
                    fontsize="x-small",
                    color="white" if v > lo + 0.6 * (1 - lo) else "black",
                )
    label = "Pearson r (limpa scale)"
    if prec is not None:
        label += f"\n{len(mat):,} proteins detected in all runs"
    fig.colorbar(im, ax=ax, label=label, fraction=0.046, pad=0.04)
    ax.grid(False)
    return fig


def plot_pca(result):
    """PCA of limpa protein log2 expression (runs >= 4), coloured by Spectronaut condition."""
    plt = _plt()
    x = result.protein_log2.T
    x = x - x.mean(axis=0)
    u, s, _ = np.linalg.svd(x.values, full_matrices=False)
    pcs = u[:, :2] * s[:2]
    var = s**2 / np.sum(s**2)
    cond = result.samples["R.Condition"].astype(str) if "R.Condition" in result.samples else None
    labels = cond.reindex(x.index) if cond is not None else pd.Series("all", index=x.index)
    cmap = dict(zip(sorted(labels.unique()), _palette(labels.nunique()), strict=True))
    fig, ax = plt.subplots(figsize=(7.5, 6))
    for k, c in cmap.items():
        m = (labels == k).values
        ax.scatter(pcs[m, 0], pcs[m, 1], s=120, color=c, label=k, edgecolor="white")
    ax.set_xlabel(f"PC1 ({var[0]:.0%})")
    ax.set_ylabel(f"PC2 ({var[1]:.0%})")
    if len(cmap) > 1:
        ax.legend(frameon=False, fontsize="small")
    return fig


def protein_agreement(prec: pd.DataFrame, result) -> pd.DataFrame:
    """Long table of limpa vs Spectronaut PG.Quantity (log2) per protein x run.

    ``*_centred`` columns subtract each protein's mean over the runs where *both* tools have a
    value -- i.e. between-run changes, which is what differential expression uses.
    """
    lim = (
        result.protein_log2.rename_axis("protein")
        .reset_index()
        .melt(id_vars="protein", var_name="run", value_name="limpa_log2")
    )
    sn = prec.groupby(["run", "protein"])["pg_quantity"].first().reset_index()
    sn["spectronaut_log2"] = np.log2(sn["pg_quantity"].where(sn["pg_quantity"] > 0))
    out = lim.merge(sn[["run", "protein", "spectronaut_log2"]], on=["run", "protein"], how="left")
    ann = result.proteins[["NPrec", "PropObs"]].rename_axis("protein").reset_index()
    out = out.merge(ann, on="protein", how="left")
    both = out["spectronaut_log2"].notna()
    for col in ("limpa_log2", "spectronaut_log2"):
        mean = out[both].groupby("protein")[col].transform("mean")
        out.loc[both, col.replace("_log2", "_centred")] = out.loc[both, col] - mean
    return out


def plot_limpa_vs_spectronaut(agree: pd.DataFrame):
    """Left: absolute levels, coloured by precursors per protein. Right: between-run changes."""
    plt = _plt()
    both = agree.dropna(subset=["spectronaut_log2"])
    bucket = _nprec_bucket(both["NPrec"])
    cats = list(bucket.cat.categories)
    colors = dict(zip(cats, _palette(len(cats)), strict=True))
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6.5))
    for c in cats[::-1]:
        m = bucket == c
        for ax, suffix in ((ax1, "_log2"), (ax2, "_centred")):
            ax.scatter(
                both.loc[m, "spectronaut" + suffix],
                both.loc[m, "limpa" + suffix],
                s=4,
                alpha=0.35,
                color=colors[c],
                rasterized=True,
                label=f"{c} precursors",
            )

    def r(a, b):
        ok = both[[a, b]].dropna()
        return np.corrcoef(ok[a], ok[b])[0, 1] if len(ok) > 2 else np.nan

    only_limpa = agree["spectronaut_log2"].isna().sum()
    ax1.text(
        0.02,
        0.97,
        f"Pearson r = {r('spectronaut_log2', 'limpa_log2'):.3f}\n"
        "offset grows with precursors per protein:\nthe tools use different absolute scales",
        transform=ax1.transAxes,
        va="top",
        fontsize="x-small",
    )
    ax1.set_xlabel("Spectronaut log2 PG.Quantity")
    ax1.set_ylabel("limpa log2 protein expression")
    ax1.legend(frameon=False, fontsize="x-small", loc="lower right", markerscale=4)
    lim = np.nanpercentile(np.abs(both[["spectronaut_centred", "limpa_centred"]].values), 99.5)
    ax2.plot([-lim, lim], [-lim, lim], color="black", lw=0.8)
    ax2.set_xlim(-lim, lim)
    ax2.set_ylim(-lim, lim)
    ax2.text(
        0.02,
        0.97,
        f"Pearson r = {r('spectronaut_centred', 'limpa_centred'):.3f}\n"
        f"{only_limpa:,} protein x run values lack a Spectronaut\nquantity (limpa estimates them)",
        transform=ax2.transAxes,
        va="top",
        fontsize="x-small",
    )
    ax2.set_xlabel("Spectronaut log2, centred per protein")
    ax2.set_ylabel("limpa log2, centred per protein")
    fig.tight_layout()
    return fig


# ----------------------------------------------------------------------------------------------
# orchestration
# ----------------------------------------------------------------------------------------------


def dia_qc_report(
    report: str | os.PathLike | None = None,
    output_dir: str | os.PathLike = "figures",
    prefix: str = "dia_qc",
    limpa_result=None,
    q_cutoff: float = 0.01,
    groups: pd.Series | None = None,
    *,
    matrix_dir: str | os.PathLike | None = None,
) -> dict:
    """Compute all QC metrics and save one PNG per check (via ``save_figure``).

    Args:
        report: Spectronaut Normal Report path. Give this or ``matrix_dir``, not both.
        output_dir: where figures are written.
        prefix: filename prefix, e.g. ``"s9e1_qc"``.
        limpa_result: a LimpaResult from :func:`deliverome_analysis.limpa.run_limpa_spectronaut`
            (same report, same q_cutoff). Without it only the pre-limpa figures are made.
        q_cutoff: q-value threshold, matching the limpa run.
        groups: optional run -> replicate-group labels for CVs (default: R.Condition).
        matrix_dir: a ``trim_spectronaut_report`` output directory. Reads its small
            ``precursor_qc.parquet`` instead of the report; ``q_cutoff`` is then the one the
            trim used.

    Returns:
        dict with ``run_summary``, ``completeness``, ``cv``, ``agreement`` (or None) tables
        and ``figures`` (name -> PNG path). Checks that do not apply (e.g. CV with no
        replicates, PCA with < 4 runs) are listed under ``skipped`` with the reason.
    """
    from deliverome_analysis.plotting import save_figure, setup_plotting_style

    setup_plotting_style()
    plt = _plt()
    out = Path(output_dir)
    if (report is None) == (matrix_dir is None):
        raise ValueError("Give exactly one of report= or matrix_dir=")
    if matrix_dir is not None:
        prec = load_trimmed_precursors(matrix_dir)
    else:
        prec = load_spectronaut_precursors(report, q_cutoff=q_cutoff)
    n_runs = prec["run"].nunique()
    summary = run_summary(prec)
    compl = precursor_completeness(prec)
    cv = precursor_cv(prec, groups)
    figures: dict[str, Path] = {}
    skipped: dict[str, str] = {}

    def save(name, fig):
        save_figure(fig, f"{prefix}_{name}", output_dir=out)
        plt.close(fig)
        figures[name] = out / f"{prefix}_{name}.png"

    save("ids_per_run", plot_ids_per_run(summary))
    save("intensity_distributions", plot_intensity_distributions(prec))
    save("completeness", plot_completeness(compl, n_runs))
    det = conditional_detection(prec)
    if det.empty:
        skipped["detection_reproducibility"] = "needs >= 2 runs"
    else:
        save("detection_reproducibility", plot_detection_reproducibility(det))
    if cv.empty:
        skipped["precursor_cv"] = "no replicate group with >= 2 runs"
    else:
        save("precursor_cv", plot_precursor_cv(cv))
    if "rt_error" in prec:
        save("rt_error", plot_rt_error(prec))
        save("gradient_coverage", plot_gradient_coverage(prec))
    else:
        skipped["rt_error"] = "report lacks EG.ApexRT / EG.RTPredicted"
    contam = top_contam = None
    if "protein_class" in prec:
        contam, top_contam = class_signal(prec), top_contaminants(prec)
        if (contam.filter(like="_precursors").drop(columns="human_precursors").sum().sum()) > 0:
            save("protein_class_signal", plot_class_signal(contam, prec))
        else:
            skipped["protein_class_signal"] = "no contaminant, watchlist or iRT precursors"
    else:
        skipped["protein_class_signal"] = "trimmed without protein classes"
    rtc = rt_summary = charge_spread = None
    if "apex_rt" in prec and n_runs >= 2:
        rtc, rt_summary = rt_consistency(prec)
        charge_spread = charge_state_rt_spread(prec)
        if not rtc.empty:
            save("rt_consistency", plot_rt_consistency(rtc, prec))
    else:
        skipped["rt_consistency"] = "needs EG.ApexRT and >= 2 runs"
    if "missed_cleavages" in prec:
        save("missed_cleavages", plot_missed_cleavages(prec))
    else:
        skipped["missed_cleavages"] = "export lacks PEP.NrOfMissedCleavages"
    for col in PEAK_METRICS:
        if col == "peak_width" and "fwhm" in prec:
            continue  # FWHM is the better width measure; show one of the two
        if col in prec and prec[col].notna().any():
            save(col, plot_peak_metric(prec, col))
        else:
            skipped[col] = "export lacks this column (see docs/spectronaut-lean-report.md)"

    agreement = None
    if limpa_result is not None:
        if limpa_result.dpc_points is not None:
            save("dpc_fit", plot_dpc_fit(limpa_result))
        else:
            skipped["dpc_fit"] = "LimpaResult has no dpc_points (older limpa driver output)"
        save("protein_se", plot_protein_se(limpa_result))
        save("run_correlation", plot_run_correlation(limpa_result, prec))
        if n_runs >= 4:
            save("pca", plot_pca(limpa_result))
        else:
            skipped["pca"] = f"only {n_runs} runs (needs >= 4)"
        if "pg_quantity" in prec:
            agreement = protein_agreement(prec, limpa_result)
            save("limpa_vs_spectronaut", plot_limpa_vs_spectronaut(agreement))
        else:
            skipped["limpa_vs_spectronaut"] = "report lacks PG.Quantity"

    return {
        "run_summary": summary,
        "completeness": compl,
        "detection": det,
        "cv": cv,
        "agreement": agreement,
        "protein_class_signal": contam,
        "top_contaminants": top_contam,
        "rt_consistency": rtc,
        "rt_summary": rt_summary,
        "charge_state_rt": charge_spread,
        "figures": figures,
        "skipped": skipped,
    }
