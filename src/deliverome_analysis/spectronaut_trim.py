"""Trim a Spectronaut Normal Report to the compact matrix limpa needs.

A full "BGS Factory Report (Normal)" carries ~80 columns and one row per precursor per run:
about 310 MB per Astral run, so ~120 GB at 384 runs -- too big to load as a long table.
limpa only needs a precursor x run log2 intensity matrix plus the protein each precursor
belongs to. This module **streams** the report in blocks (memory bounded by the output
matrix, not the report) and writes:

    <out_dir>/precursor_log2.parquet   precursor x run log2 intensity, NaN = not detected
    <out_dir>/precursors.parquet       precursor, protein id, accessions/genes (if exported)
    <out_dir>/samples.tsv              run, condition, replicate (from R.* columns)
    <out_dir>/precursor_qc.parquet     long precursor x run table of the QC columns
                                       (``PRECURSOR_QC_COLUMNS``) the report carries, for
                                       ``deliverome_analysis.dia_qc`` -- so QC never reopens
                                       the full report
    <out_dir>/manifest.json            input file, size, mtime, filters, counts, versions

Filtering is identical to ``limpa::readSpectronaut``: EG.Qvalue and PG.Qvalue <= q_cutoff,
values flagged EG.IsImputed dropped, zero intensities treated as missing, intensities from
``EG.TotalQuantity (Settings)``. Precursor ids are ``<EG.ModifiedSequence>.<FG.Charge>`` (limpa's
key), and values are float64 so results match reading the report directly.

Pass the output directory to ``run_limpa_spectronaut(matrix_dir=...)``.

Exporting less in the first place: Spectronaut's ``.sne`` is a proprietary binary that only
Spectronaut reads, so trimming at the source means a report schema containing just
``EXPORT_COLUMNS`` (save one from Spectronaut's report editor and use it for CLI exports).
docs/spectronaut-lean-report.md is the how-to; print the list with
``python -m deliverome_analysis.spectronaut_trim --columns``.
"""

from __future__ import annotations

import json
import os
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.csv as pacsv
import pyarrow.parquet as pq

REQUIRED = [
    "R.FileName",
    "EG.ModifiedSequence",
    "FG.Charge",
    "EG.TotalQuantity (Settings)",
    "EG.Qvalue",
    "PG.Qvalue",
]
OPTIONAL = ["EG.IsImputed", "R.Condition", "R.Replicate", "PG.ProteinAccessions", "PG.Genes"]
# Per-precursor annotation carried into precursors.parquet (first value seen).
# FG.XICDBID is Spectronaut's key into its XIC database export (one SQLite per raw file).
_ANNOTATION = ["PEP.StrippedSequence", "PEP.NrOfMissedCleavages", "FG.XICDBID"]
# Per precursor-run columns carried into precursor_qc.parquet when the report has them, mapped
# to the names deliverome_analysis.dia_qc uses. The first block is in Spectronaut's factory
# Normal report; the peak-shape block needs the lean template (Spectronaut 21 names -- confirm
# the spelling in the report editor; a column that is absent is simply skipped).
PRECURSOR_QC_COLUMNS = {
    "PG.Quantity": "pg_quantity",
    "EG.ApexRT": "apex_rt",
    "EG.RTPredicted": "rt_predicted",
    "FG.PrecMz": "prec_mz",
    "FG.PrecMzCalibrated": "prec_mz_calibrated",
    "EG.iRTEmpirical": "irt_empirical",
    "EG.iRTPredicted": "irt_predicted",
    "EG.Cscore": "cscore",
    # lean template: peak shape, observed mass error and identification subscores. Names
    # checked against the strings in the Spectronaut 21.1 libraries (2026-10-07).
    "EG.StartRT": "start_rt",
    "EG.EndRT": "end_rt",
    "EG.FWHM": "fwhm",
    "EG.DatapointsPerPeak": "points_per_peak",
    "EG.DatapointsPerPeak (MS1)": "points_per_peak_ms1",
    "EG.SignalToNoise": "signal_to_noise",
    "EG.MeanTailingFactor": "tailing_factor",
    "FG.CalibratedMassAccuracy (PPM)": "ms1_mass_error_ppm",
    "EG.PeakGaussScore": "peak_gauss_score",
    "FG.FragmentCoelution": "fragment_coelution",
    "EG.IonLadderScore": "ion_ladder_score",
    # MS1 and MS2 precursor areas, for signal shares by level (e.g. the surface-protein
    # share of MS1 vs MS2 signal). Names checked the same way.
    "FG.MS1Quantity": "ms1_quantity",
    "FG.MS2Quantity": "ms2_quantity",
}
# Columns worth exporting from Spectronaut for limpa + deliverome_analysis.dia_qc.
EXPORT_COLUMNS = list(
    dict.fromkeys(
        [
            *REQUIRED,
            "PG.ProteinGroups",
            *OPTIONAL,
            *_ANNOTATION,
            *PRECURSOR_QC_COLUMNS,
        ]
    )
)
# Fragment-level template (one row per fragment ion): b/y coverage, MS2 mass error, and the
# ion labels for Spectronaut's XIC database. Names checked like EXPORT_COLUMNS.
FRAGMENT_EXPORT_COLUMNS = [
    "R.FileName",
    "EG.ModifiedSequence",
    "FG.Charge",
    "FG.XICDBID",
    "F.FrgIon",
    "F.FrgType",
    "F.FrgNum",
    "F.Charge",
    "F.FrgLossType",
    "F.FrgMz",
    "F.CalibratedMassAccuracy (PPM)",
    "F.PeakArea",
    "F.NormalizedPeakArea",
    "F.ExcludedFromQuantification",
    "F.PossibleInterference",
    "F.InterferenceScore",
]


def report_header(path: str | os.PathLike) -> list[str]:
    """Column names of a Spectronaut report (TSV or Parquet) without reading the rows."""
    path = Path(path)
    if path.suffix == ".parquet":
        return pq.read_schema(path).names
    return pacsv.open_csv(
        path,
        read_options=pacsv.ReadOptions(block_size=1 << 20),
        parse_options=pacsv.ParseOptions(delimiter="\t"),
    ).schema.names


def check_export(path: str | os.PathLike, template: str = "lean") -> dict[str, list[str]]:
    """Compare an export's columns with a template's list: ``missing`` and ``extra``.

    Run this on the first export made with a hand-built ``.rs`` template: a column that was
    not ticked, or ticked under a slightly different name, shows up as missing.
    """
    want = EXPORT_COLUMNS if template == "lean" else FRAGMENT_EXPORT_COLUMNS
    have = report_header(path)
    return {
        "missing": [c for c in want if c not in have],
        "extra": [c for c in have if c not in want],
    }


_NUMERIC = {"FG.Charge", "EG.TotalQuantity (Settings)", "EG.Qvalue", "PG.Qvalue"}


@dataclass
class TrimResult:
    out_dir: Path
    n_runs: int
    n_precursors: int
    n_values: int
    dropped_q: int
    dropped_imputed: int
    seconds: float


def trim_spectronaut_report(
    report: str | os.PathLike,
    out_dir: str | os.PathLike,
    q_cutoff: float = 0.01,
    protein_id: str = "PG.ProteinGroups",
    block_size: int = 256 << 20,
    *,
    watchlist: frozenset[str] | None = None,
    exclude_classes: tuple[str, ...] = ("contaminant", "irt"),
) -> TrimResult:
    """Stream ``report`` into a wide precursor x run log2 matrix for limpa (see module doc).

    Each precursor's protein group is classified (``deliverome_analysis.contaminants``):
    groups whose class is in ``exclude_classes`` (contaminant-only groups and the iRT standard,
    by default) are left out of the matrix limpa quantifies and listed in
    ``excluded_precursors.parquet``; the QC table keeps every precursor with its class.
    ``watchlist``: human accessions to flag (e.g. ``contaminants.load_watchlist()``).
    """
    t0 = time.time()
    report = Path(report)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    header = pacsv.open_csv(
        report,
        read_options=pacsv.ReadOptions(block_size=1 << 20),
        parse_options=pacsv.ParseOptions(delimiter="\t"),
    ).schema.names
    required = [*REQUIRED, protein_id]
    missing = [c for c in required if c not in header]
    if missing:
        raise ValueError(f"Report is missing required column(s): {missing}")
    qc_cols = [c for c in PRECURSOR_QC_COLUMNS if c in header]
    extra = [c for c in [*OPTIONAL, *_ANNOTATION] if c in header]
    use = list(dict.fromkeys(required + extra + qc_cols))
    text = {c: pa.string() for c in use if c not in _NUMERIC and c not in qc_cols}
    reader = pacsv.open_csv(
        report,
        read_options=pacsv.ReadOptions(block_size=block_size),
        parse_options=pacsv.ParseOptions(delimiter="\t"),
        convert_options=pacsv.ConvertOptions(
            include_columns=use,
            column_types={
                **text,
                "FG.Charge": pa.int32(),
                "EG.TotalQuantity (Settings)": pa.float64(),
                "EG.Qvalue": pa.float64(),
                "PG.Qvalue": pa.float64(),
                **{c: pa.float64() for c in qc_cols},
            },
            strings_can_be_null=True,
        ),
    )

    prec_index: dict[str, int] = {}
    run_index: dict[str, int] = {}
    ann_rows: list[pd.DataFrame] = []
    run_rows: list[pd.DataFrame] = []
    rows, cols, vals = [], [], []
    dropped_q = dropped_imp = 0
    # QC rows stream straight to disk so memory stays bounded by the matrix at any run count.
    qc_writer = None
    qc_path = out / "precursor_qc.parquet"

    for batch in reader:
        tb = pa.Table.from_batches([batch])
        ok_q = pc.and_(
            pc.less_equal(tb["EG.Qvalue"], q_cutoff), pc.less_equal(tb["PG.Qvalue"], q_cutoff)
        )
        ok_q = pc.fill_null(ok_q, False)
        dropped_q += len(tb) - pc.sum(ok_q).as_py() if len(tb) else 0
        keep = ok_q
        if "EG.IsImputed" in tb.column_names:
            imp = pc.equal(pc.utf8_lower(pc.cast(tb["EG.IsImputed"], pa.string())), "true")
            imp = pc.fill_null(imp, False)
            dropped_imp += pc.sum(pc.and_(keep, imp)).as_py() or 0
            keep = pc.and_(keep, pc.invert(imp))
        keep = pc.and_(keep, pc.fill_null(pc.greater(tb["EG.TotalQuantity (Settings)"], 0), False))
        df = tb.filter(keep).to_pandas()
        if df.empty:
            continue
        prec = df["EG.ModifiedSequence"] + "." + df["FG.Charge"].astype(str)

        new_p = pd.unique(prec[~prec.isin(prec_index.keys())])
        if len(new_p):
            first = df.loc[~prec.duplicated() & prec.isin(new_p)]
            ann = pd.DataFrame({"precursor": prec[first.index]})
            for c in [protein_id, "PG.ProteinAccessions", "PG.Genes", *_ANNOTATION]:
                if c in df:
                    ann[c] = first[c].values
            ann_rows.append(ann)
            base = len(prec_index)
            prec_index.update({p: base + i for i, p in enumerate(ann["precursor"])})
        run = df["R.FileName"]
        new_r = pd.unique(run[~run.isin(run_index.keys())])
        if len(new_r):
            first = df.loc[~run.duplicated() & run.isin(new_r)]
            info = pd.DataFrame({"run": first["R.FileName"].values})
            for c in ("R.Condition", "R.Replicate"):
                if c in df:
                    info[c] = first[c].values
            run_rows.append(info)
            base = len(run_index)
            run_index.update({r: base + i for i, r in enumerate(info["run"])})

        rows.append(prec.map(prec_index).to_numpy(np.int32))
        cols.append(run.map(run_index).to_numpy(np.int32))
        vals.append(np.log2(df["EG.TotalQuantity (Settings)"].to_numpy(np.float64)))

        qc = pd.DataFrame({"precursor": prec.values, "run": run.values, "log2_intensity": vals[-1]})
        for c in qc_cols:
            if c not in ("FG.PrecMz", "FG.PrecMzCalibrated"):
                qc[PRECURSOR_QC_COLUMNS[c]] = df[c].to_numpy(np.float32)
        if {"FG.PrecMz", "FG.PrecMzCalibrated"} <= set(qc_cols):
            # MS1 mass correction Spectronaut applied, in ppm; computed in float64 because
            # float32 m/z is only good to ~0.06 ppm, the same order as the offsets themselves.
            mz = df["FG.PrecMz"].to_numpy(np.float64)
            qc["ms1_calibration_ppm"] = (
                (df["FG.PrecMzCalibrated"].to_numpy(np.float64) - mz) / mz * 1e6
            ).astype(np.float32)
        qc = qc.drop_duplicates(["precursor", "run"])  # first wins, as in the matrix
        qc_tb = pa.Table.from_pandas(qc, preserve_index=False)
        if qc_writer is None:
            qc_writer = pq.ParquetWriter(qc_path, qc_tb.schema, compression="zstd")
        qc_writer.write_table(qc_tb)

    if qc_writer is not None:
        qc_writer.close()
    if not rows:
        raise ValueError("No values pass the q-value / imputation filters")
    r = np.concatenate(rows)
    c = np.concatenate(cols)
    v = np.concatenate(vals)
    mat = np.full((len(prec_index), len(run_index)), np.nan)
    # Spectronaut writes each precursor once per run; if a schema repeats rows (e.g. fragment
    # columns), keep the first value, matching how a long-format reader would see one value.
    order = np.lexsort((np.arange(len(r)), c, r))
    first = np.ones(len(r), bool)
    first[1:] = (r[order][1:] != r[order][:-1]) | (c[order][1:] != c[order][:-1])
    sel = order[first]
    mat[r[sel], c[sel]] = v[sel]

    runs = list(run_index)
    precs = list(prec_index)
    ann = pd.concat(ann_rows, ignore_index=True).rename(columns={protein_id: "protein_id"})
    from deliverome_analysis.contaminants import classify

    acc_col = "PG.ProteinAccessions" if "PG.ProteinAccessions" in ann else "protein_id"
    ann["protein_class"] = classify(ann[acc_col], watchlist or frozenset()).to_numpy()
    keep = ~ann["protein_class"].isin(exclude_classes).to_numpy()
    wide = pa.table(
        {
            "precursor": [p for p, k in zip(precs, keep, strict=True) if k],
            **{run: mat[keep, j] for j, run in enumerate(runs)},
        }
    )
    pq.write_table(wide, out / "precursor_log2.parquet", compression="zstd")
    # limpa needs precursors.parquet to match the matrix row for row
    pq.write_table(
        pa.Table.from_pandas(ann[keep], preserve_index=False), out / "precursors.parquet"
    )
    pq.write_table(
        pa.Table.from_pandas(ann[~keep], preserve_index=False),
        out / "excluded_precursors.parquet",
    )
    pd.concat(run_rows, ignore_index=True).to_csv(out / "samples.tsv", sep="\t", index=False)

    stat = report.stat()
    result = TrimResult(
        out,
        len(runs),
        int(keep.sum()),
        int(np.isfinite(mat[keep]).sum()),
        int(dropped_q),
        int(dropped_imp),
        round(time.time() - t0, 1),
    )
    (out / "manifest.json").write_text(
        json.dumps(
            {
                "source_report": str(report.resolve()),
                "source_bytes": stat.st_size,
                "source_mtime": time.strftime("%Y-%m-%dT%H:%M:%S", time.gmtime(stat.st_mtime)),
                "protein_id_column": protein_id,
                "q_cutoff": q_cutoff,
                "imputed_flag_present": "EG.IsImputed" in use,
                "qc_columns": {c: PRECURSOR_QC_COLUMNS[c] for c in qc_cols},
                "protein_classes": ann["protein_class"].value_counts().to_dict(),
                "excluded_from_matrix": {
                    "classes": list(exclude_classes),
                    "precursors": int((~keep).sum()),
                },
                "watchlist_size": len(watchlist or ()),
                "n_runs": result.n_runs,
                "n_precursors": result.n_precursors,
                "n_values": result.n_values,
                "dropped_failed_q": result.dropped_q,
                "dropped_imputed": result.dropped_imputed,
                "pyarrow_version": pa.__version__,
                "seconds": result.seconds,
            },
            indent=2,
        )
        + "\n"
    )
    return result


if __name__ == "__main__":
    import sys

    args = sys.argv[1:]
    if args == ["--columns"]:
        print("\n".join(EXPORT_COLUMNS))
    elif args == ["--fragment-columns"]:
        print("\n".join(FRAGMENT_EXPORT_COLUMNS))
    elif len(args) in (2, 3) and args[0] == "--check":
        tmpl = args[2] if len(args) == 3 else "lean"
        res = check_export(args[1], tmpl)
        print(f"{tmpl} template: {len(res['missing'])} missing, {len(res['extra'])} extra")
        for c in res["missing"]:
            print(f"  MISSING  {c}")
        for c in res["extra"]:
            print(f"  extra    {c}")
        sys.exit(1 if res["missing"] else 0)
    else:
        sys.exit(
            "usage: python -m deliverome_analysis.spectronaut_trim "
            "--columns | --fragment-columns | --check <report> [lean|fragment]"
        )
