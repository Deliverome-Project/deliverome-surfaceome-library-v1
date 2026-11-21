"""Turn an alphaDIA output folder into the trimmed directory the DIA QC suite reads.

:func:`deliverome_analysis.spectronaut_trim.trim_spectronaut_report` defines a small on-disk
contract (``precursor_log2.parquet``, ``precursors.parquet``, ``excluded_precursors.parquet``,
``precursor_qc.parquet``, ``samples.tsv``, ``manifest.json``). Everything downstream --
``dia_qc.dia_qc_report(matrix_dir=...)``, the dashboard, ``dia_qc_pipeline``, ``qcdb`` and
``limpa.run_limpa_spectronaut(matrix_dir=...)`` -- reads only that contract. This module writes
the same files from alphaDIA 2.x's ``precursors.parquet`` (or ``.tsv``), so the whole suite
runs on alphaDIA searches unchanged.

Filtering mirrors the Spectronaut path: precursor q <= ``q_cutoff`` and protein-group
q <= ``q_cutoff``, targets only, intensity > 0. alphaDIA does not impute.

Column mapping (alphaDIA, checked against ``alphadia/constants/keys.py`` at v2.4.0):

=============================  ===============================  ==============================
alphaDIA                       trimmed contract                 note
=============================  ===============================  ==============================
``raw.name``                   ``run``
``precursor.sequence``         ``PEP.StrippedSequence``
``sequence + mods + charge``   ``precursor``                    ``SEQ[mods@sites].z``
``pg.name``                    ``protein_id``                   alphaDIA's protein group
``pg.proteins``                ``PG.ProteinAccessions``         every accession the precursor
                                                                maps to (``;``-joined); used to
                                                                classify contaminants
``pg.genes``                   ``PG.Genes``                     UniProt entry names
                                                                (e.g. ``HNRPD_HUMAN``), not
                                                                gene symbols
``precursor.intensity``        ``log2_intensity``               log2
``pg.intensity``               ``pg_quantity``
``precursor.rt.observed``      ``apex_rt``                      seconds -> minutes
``precursor.rt.calibrated``    ``rt_predicted``                 seconds -> minutes
``precursor.rt.fwhm``          ``fwhm``                         seconds -> minutes
``precursor.mz.observed/lib``  ``ms1_mass_error_ppm``           (observed - library) / library;
                                                                omitted when all zero (no MS1
                                                                calibration in the search)
(computed from sequence)       ``PEP.NrOfMissedCleavages``      trypsin rule, K/R not before P
``fwhm`` + cycle time          ``points_per_peak``              MS2 points across the whole
                                                                peak, 1.699 x FWHM / cycle time
                                                                (Gaussian: 4 sigma = 1.699 FWHM);
                                                                needs the run's cycle time from
                                                                the runner's ``benchmark.json``
=============================  ===============================  ==============================

``rt_error`` downstream is therefore alphaDIA's RT calibration residual, the analogue of
Spectronaut's ApexRT - RTPredicted. ``points_per_peak`` is an estimate from alphaDIA's fitted
FWHM, comparable to Spectronaut's ``EG.DatapointsPerPeak`` (data points across the peak) but not
counted point by point. Spectronaut-only metrics (Cscore, iRT, S/N, tailing, Gaussian score)
have no alphaDIA equivalent and are left out; QC reports them as not available.

    uv run python -m deliverome_analysis.alphadia_trim <alphadia-output-dir> trimmed/
"""

from __future__ import annotations

import argparse
import json
import os
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from deliverome_analysis.peptide_metrics import GAUSSIAN_BASE_PER_FWHM  # one definition

StrPath = str | os.PathLike[str]
RunInfo = dict[str, dict[str, Any]]  # run -> fields from the search record (search_run_info)

COLUMNS = {
    "raw.name": "run",
    "precursor.sequence": "sequence",
    "precursor.mods": "mods",
    "precursor.mod_sites": "mod_sites",
    "precursor.charge": "charge",
    "precursor.intensity": "intensity",
    "precursor.qval": "precursor_q",
    "precursor.decoy": "decoy",
    "pg.name": "pg_name",
    "pg.proteins": "proteins",
    "pg.genes": "genes",
    "pg.qval": "pg_q",
    "pg.intensity": "pg_quantity",
    "precursor.rt.observed": "apex_rt",
    "precursor.rt.calibrated": "rt_predicted",
    "precursor.rt.fwhm": "fwhm",
    "precursor.mz.library": "mz_library",
    "precursor.mz.observed": "mz_observed",
}
REQUIRED = [
    "raw.name",
    "precursor.sequence",
    "precursor.charge",
    "precursor.intensity",
    "precursor.qval",
    "pg.proteins",
]
SECONDS_TO_MINUTES = ("apex_rt", "rt_predicted", "fwhm")
QC_COLUMNS = (
    "pg_quantity",
    "apex_rt",
    "rt_predicted",
    "fwhm",
    "points_per_peak",
    "ms1_mass_error_ppm",
)


@dataclass
class AlphaDIATrimResult:
    """Summary of one :func:`trim_alphadia_output` call."""

    out_dir: Path
    n_runs: int
    n_precursors: int
    n_values: int
    dropped_q: int
    seconds: float


def find_precursor_table(alphadia_dir: StrPath) -> Path:
    """alphaDIA's precursor table in ``alphadia_dir`` (parquet preferred over tsv)."""
    d = Path(alphadia_dir)
    for name in ("precursors.parquet", "precursors.tsv"):
        if (d / name).exists():
            return d / name
    raise FileNotFoundError(f"No precursors.parquet or precursors.tsv in {d}")


def is_alphadia_output(path: StrPath) -> bool:
    """True when ``path`` is a folder holding an alphaDIA 2.x precursor table."""
    try:
        table = find_precursor_table(path)
    except FileNotFoundError, NotADirectoryError:
        return False
    return set(REQUIRED) <= set(_columns(table))


def _columns(path: Path) -> list[str]:
    if path.suffix == ".parquet":
        return list(pq.read_schema(path).names)
    return list(pd.read_csv(path, sep="\t", nrows=0).columns)


def tryptic_missed_cleavages(sequences: pd.Series) -> pd.Series:
    """Internal K/R not followed by P; the C-terminal residue is never counted."""
    return sequences.str[:-1].str.count(r"[KR](?!P)").astype("int64")


def precursor_ids(
    seq: pd.Series, mods: pd.Series, sites: pd.Series, charge: pd.Series
) -> pd.Series:
    """``SEQ.z`` for unmodified precursors, ``SEQ[mods@sites].z`` otherwise."""
    mods = mods.fillna("").astype(str)
    sites = sites.fillna("").astype(str)
    tag = ("[" + mods + "@" + sites + "]").where(mods != "", "")
    return seq + tag + "." + charge.astype(int).astype(str)


def search_run_info(alphadia_dir: StrPath) -> RunInfo:
    """Per run: acquisition (raw header) and search library, from the runner's benchmark.json.

    deliverome-internal's ``run_benchmark.py`` writes ``benchmark.json`` next to the alphaDIA
    output folder. Returns ``{run: {acquired, instrument_*, file_bytes, search_raw_sha256,
    search_raw_bytes, run_uid, search_library, search_library_sha256}}`` (run = raw file
    stem, as in ``raw.name``; ``search_raw_*`` = the raw file the search read), or
    ``{}`` when there is no record. Lets a QC page carry the same acquisition key
    (instrument serial | acquisition time) as the search index without reading raw files.
    Also carries the run's cycle time (raw method) and, from ``alphadia/stat.tsv``, the MS2/MS1
    mass calibration (bias in ppm, variance; MS1 is None when alphaDIA did not calibrate on
    MS1) and the search tolerances.
    """
    record_path = Path(alphadia_dir).parent / "benchmark.json"
    if not record_path.exists():
        return {}
    stats: dict[str, dict[str, float | None]] = {}
    stat_path = Path(alphadia_dir) / "stat.tsv"
    if stat_path.exists():  # alphaDIA's per-run search statistics: mass calibration, tolerances
        for _, row in pd.read_csv(stat_path, sep="\t").iterrows():

            def num(col: str, row: Any = row) -> float | None:
                v = row.get(col)
                return float(v) if v is not None and pd.notna(v) else None

            stats[str(row.get("raw.name"))] = {
                "ms2_mass_bias_ppm": num("calibration.ms2_bias"),
                "ms2_mass_variance": num("calibration.ms2_variance"),
                "ms1_mass_bias_ppm": num("calibration.ms1_bias"),
                "ms1_mass_variance": num("calibration.ms1_variance"),
                "ms2_tolerance_ppm": num("optimization.ms2_error"),
                "ms1_tolerance_ppm": num("optimization.ms1_error"),
            }
    rec = json.loads(record_path.read_text())
    lib = rec.get("library") or {}
    out: RunInfo = {}
    for item in rec.get("inputs", []):
        acq = item.get("acquisition") or {}
        out[Path(item["path"]).stem] = {
            **acq,
            "file_bytes": item.get("bytes"),
            # the file the search actually read, checked against raw:// by the QC index
            "search_raw_sha256": item.get("sha256"),
            "search_raw_bytes": item.get("bytes"),
            "run_uid": item.get("run_uid"),
            "search_library": Path(lib["path"]).name if lib.get("path") else "predicted",
            "search_library_sha256": lib.get("sha256"),
            "cycle_time_s": (item.get("method") or {}).get("cycle_time_s"),
            **stats.get(Path(item["path"]).stem, {}),
        }
    return out


def trim_alphadia_output(
    alphadia_dir: StrPath | list[StrPath] | tuple[StrPath, ...],
    out_dir: StrPath,
    q_cutoff: float = 0.01,
    *,
    conditions: dict[str, str] | None = None,
    watchlist: frozenset[str] | None = None,
    exclude_classes: tuple[str, ...] = ("contaminant", "irt"),
) -> AlphaDIATrimResult:
    """Write the trimmed QC/limpa directory for alphaDIA output folder(s) (see module doc).

    ``alphadia_dir`` may be a list of folders, e.g. one alphaDIA search per raw file; their
    precursor tables are concatenated so cross-run QC (CVs, RT consistency, completeness)
    works. Each folder keeps its own FDR and protein grouping: a precursor's protein group
    in ``precursors.parquet`` is the one from the first folder that reports it.

    ``conditions`` maps run name -> condition (``R.Condition`` in ``samples.tsv``); runs not
    listed get ``"all"``. ``watchlist`` and ``exclude_classes`` behave exactly as in
    :func:`~deliverome_analysis.spectronaut_trim.trim_spectronaut_report`.
    """
    t0 = time.time()
    dirs = list(alphadia_dir) if isinstance(alphadia_dir, (list, tuple)) else [alphadia_dir]
    sources = [find_precursor_table(d) for d in dirs]
    tables = []
    for src in sources:
        header = _columns(src)
        missing = [c for c in REQUIRED if c not in header]
        if missing:
            raise ValueError(f"{src} is missing required alphaDIA column(s): {missing}")
        use = [c for c in COLUMNS if c in header]
        if src.suffix == ".parquet":
            tables.append(pd.read_parquet(src, columns=use))
        else:
            tables.append(pd.read_csv(src, sep="\t", usecols=use))
    df = pd.concat(tables, ignore_index=True).rename(columns=COLUMNS)
    n_in = len(df)

    keep = (df["precursor_q"] <= q_cutoff) & (df["intensity"] > 0)
    if "pg_q" in df:
        keep &= df["pg_q"] <= q_cutoff
    if "decoy" in df:
        keep &= df["decoy"].astype(int) == 0
    df = df.loc[keep.fillna(False)].copy()
    if df.empty:
        raise ValueError("No alphaDIA precursors pass the q-value filters")
    dropped_q = n_in - len(df)

    blank = pd.Series("", index=df.index)
    df["precursor"] = precursor_ids(
        df["sequence"], df.get("mods", blank), df.get("mod_sites", blank), df["charge"]
    )
    if df.duplicated(["precursor", "run"]).any():
        raise ValueError("alphaDIA output has duplicate (precursor, run) rows after filtering")
    notes: list[str] = []
    run_info: RunInfo = {}
    for d in dirs:
        run_info |= search_run_info(d)
    cycle = {r: v.get("cycle_time_s") for r, v in run_info.items() if v.get("cycle_time_s")}
    if "fwhm" in df and cycle:
        # FWHM is still in seconds here
        df["points_per_peak"] = (
            GAUSSIAN_BASE_PER_FWHM * df["fwhm"] / df["run"].map(cycle).astype(np.float64)
        )
        missing_cycle = sorted(set(df["run"]) - set(cycle))
        if missing_cycle:
            notes.append(f"points_per_peak missing for runs without a cycle time: {missing_cycle}")
    elif "fwhm" in df:
        notes.append("points_per_peak omitted: no cycle time in the runs' benchmark.json")
    for col in SECONDS_TO_MINUTES:
        if col in df:
            df[col] = df[col] / 60.0
    if {"mz_library", "mz_observed"} <= set(df.columns):
        mz = df["mz_library"].to_numpy(np.float64)
        ppm = (df["mz_observed"].to_numpy(np.float64) - mz) / mz * 1e6
        if np.any(ppm[np.isfinite(ppm)] != 0):
            df["ms1_mass_error_ppm"] = ppm
        else:
            # alphaDIA 2.4.0 copies the library m/z into precursor.mz.observed when the search
            # does not calibrate on MS1 (stat.tsv: calibration.ms1_bias = NaN), so a 0 ppm
            # "error" would be an artefact, not a measurement.
            notes.append("ms1_mass_error_ppm omitted: precursor.mz.observed equals the library m/z")
    df["log2_intensity"] = np.log2(df["intensity"].to_numpy(np.float64))

    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    # Per-precursor annotation, first row seen (as the Spectronaut trimmer does).
    first = df.drop_duplicates("precursor")
    ann = pd.DataFrame(
        {
            "precursor": first["precursor"].to_numpy(),
            # pg.name is the protein group alphaDIA reports (pg.matrix rows); older outputs
            # without it fall back to the precursor's accession list.
            "protein_id": first.get("pg_name", first["proteins"]).astype(str).to_numpy(),
            "PG.ProteinAccessions": first["proteins"].astype(str).to_numpy(),
            "PEP.StrippedSequence": first["sequence"].to_numpy(),
            "PEP.NrOfMissedCleavages": tryptic_missed_cleavages(first["sequence"]).to_numpy(),
        }
    )
    if "genes" in first:
        ann["PG.Genes"] = first["genes"].to_numpy()
    from deliverome_analysis.contaminants import classify

    ann["protein_class"] = classify(
        ann["PG.ProteinAccessions"], watchlist or frozenset()
    ).to_numpy()
    in_matrix = ~ann["protein_class"].isin(exclude_classes).to_numpy()

    runs = list(pd.unique(df["run"]))
    wide = df.pivot(index="precursor", columns="run", values="log2_intensity").reindex(
        index=ann["precursor"], columns=runs
    )
    matrix_ann = ann[in_matrix].reset_index(drop=True)
    wide_kept = wide.loc[matrix_ann["precursor"]]
    pq.write_table(
        pa.table(
            {
                "precursor": matrix_ann["precursor"].tolist(),
                **{r: wide_kept[r].to_numpy(np.float64) for r in runs},
            }
        ),
        out / "precursor_log2.parquet",
        compression="zstd",
    )
    pq.write_table(
        pa.Table.from_pandas(matrix_ann, preserve_index=False), out / "precursors.parquet"
    )
    pq.write_table(
        pa.Table.from_pandas(ann[~in_matrix].reset_index(drop=True), preserve_index=False),
        out / "excluded_precursors.parquet",
    )

    qc = df[["precursor", "run", "log2_intensity"]].copy()
    for col in QC_COLUMNS:
        if col in df:
            qc[col] = df[col].to_numpy(np.float32)
    pq.write_table(
        pa.Table.from_pandas(qc.reset_index(drop=True), preserve_index=False),
        out / "precursor_qc.parquet",
        compression="zstd",
    )

    conditions = conditions or {}
    pd.DataFrame({"run": runs, "R.Condition": [conditions.get(r, "all") for r in runs]}).to_csv(
        out / "samples.tsv", sep="\t", index=False
    )

    n_values = int(np.isfinite(wide_kept.to_numpy()).sum())
    result = AlphaDIATrimResult(
        out, len(runs), int(in_matrix.sum()), n_values, int(dropped_q), round(time.time() - t0, 1)
    )
    stat = sources[0].stat()
    run_info = {r: v for r, v in run_info.items() if r in set(runs)}
    (out / "manifest.json").write_text(
        json.dumps(
            {
                "source_engine": "alphadia",
                "source_report": str(sources[0].resolve()),
                "source_reports": [str(p.resolve()) for p in sources],
                "source_bytes": sum(p.stat().st_size for p in sources),
                "source_mtime": time.strftime("%Y-%m-%dT%H:%M:%S", time.gmtime(stat.st_mtime)),
                "protein_id_column": "pg.proteins",
                "q_cutoff": q_cutoff,
                "imputed_flag_present": False,
                "qc_columns": {c: c for c in QC_COLUMNS if c in df},
                "protein_classes": ann["protein_class"].value_counts().to_dict(),
                "excluded_from_matrix": {
                    "classes": list(exclude_classes),
                    "precursors": int((~in_matrix).sum()),
                },
                "notes": notes,
                "watchlist_size": len(watchlist or ()),
                "n_runs": result.n_runs,
                "n_precursors": result.n_precursors,
                "n_values": result.n_values,
                "dropped_failed_q": result.dropped_q,
                "dropped_imputed": 0,
                "run_info": run_info,
                "pyarrow_version": pa.__version__,
                "seconds": result.seconds,
            },
            indent=2,
        )
        + "\n"
    )
    return result


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description="alphaDIA output folder -> trimmed QC directory")
    ap.add_argument(
        "alphadia_dir", nargs="+", help="alphaDIA output folder(s) (each holds precursors.parquet)"
    )
    ap.add_argument("--out", dest="out_dir", required=True, help="trimmed directory to write")
    ap.add_argument("--q-cutoff", type=float, default=0.01)
    ap.add_argument(
        "--no-watchlist", action="store_true", help="skip the S3 contaminant watchlist (offline)"
    )
    a = ap.parse_args(argv)
    watchlist = None
    if not a.no_watchlist:
        from deliverome_analysis.contaminants import load_watchlist

        watchlist = load_watchlist()
    r = trim_alphadia_output(a.alphadia_dir, a.out_dir, a.q_cutoff, watchlist=watchlist)
    print(
        f"{r.n_runs} run(s), {r.n_precursors} precursors in matrix, {r.n_values} values, "
        f"{r.dropped_q} rows failed q <= {a.q_cutoff} -> {r.out_dir} ({r.seconds} s)"
    )


if __name__ == "__main__":
    main()
