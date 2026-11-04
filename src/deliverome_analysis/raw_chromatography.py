"""Run-level chromatography from Thermo .raw files: TIC, spray and LC traces.

The Spectronaut report says nothing about the chromatography itself. These read it from the
raw file so QC can show the TIC up front and tell instrument problems apart:

* :func:`scan_traces` -- one row per scan: RT, MS order, TIC, base-peak intensity and ion
  injection (fill) time. Reads scan headers only, never the spectra, so it is fast.
* :func:`status_log` -- the instrument status log (spray voltage/current, source
  temperatures, vacuum ...) sampled on an RT grid, numeric entries only.
* :func:`device_traces` -- every non-MS device the acquisition recorded (LC pump pressure
  and flow on Evosep / Vanquish, UV ...), one long table.

Which status keys and devices exist depends on the instrument and LC method, so every
function returns whatever is there and :func:`chromatography` reports what it found.
Run on Modal next to the cached raw files: ``scripts/raw_chromatography_modal.py``.
Requires ``alpharaw`` (Thermo RawFileReader through pythonnet/mono).
"""

from __future__ import annotations

import os
import re

import numpy as np
import pandas as pd

_INJECT = re.compile(r"ion injection time", re.IGNORECASE)
_NUM = re.compile(r"^[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?$")


def _open(raw_path: str | os.PathLike):
    from alpharaw.raw_access.pythermorawfilereader import RawFileReader

    return RawFileReader(str(raw_path))


def _injection_label(r, scan: int) -> str | None:
    labels = list(r.source.GetTrailerExtraInformation(scan).Labels)
    return next((lab for lab in labels if _INJECT.search(lab)), None)


def scan_traces(raw_path: str | os.PathLike, *, injection_every: int = 1) -> pd.DataFrame:
    """Per scan: ``scan, rt, ms_order, tic, base_peak, inject_ms``.

    ``inject_ms`` is the ion injection (fill) time from the scan trailer, read for MS1 scans
    and every ``injection_every``-th MS2 scan (NaN otherwise); a rising fill time at constant
    TIC target means less ion current is arriving.
    """
    r = _open(raw_path)
    try:
        first, last = r.FirstSpectrumNumber, r.LastSpectrumNumber
        label = _injection_label(r, first)
        rows = []
        for s in range(first, last + 1):
            stats = r.source.GetScanStatsForScanNumber(s)
            order = r.GetMSOrderForScanNum(s)
            inj = np.nan
            if label and (order == 1 or (s - first) % injection_every == 0):
                trailer = r.source.GetTrailerExtraInformation(s)
                vals = dict(zip(trailer.Labels, trailer.Values, strict=False))
                v = str(vals.get(label, "")).strip()
                inj = float(v) if _NUM.match(v) else np.nan
            rows.append((s, r.RTFromScanNum(s), order, stats.TIC, stats.BasePeakIntensity, inj))
    finally:
        r.Close()
    return pd.DataFrame(rows, columns=["scan", "rt", "ms_order", "tic", "base_peak", "inject_ms"])


def status_log(raw_path: str | os.PathLike, *, step_min: float = 0.2) -> pd.DataFrame:
    """Numeric status-log entries every ``step_min`` minutes: long table ``rt, key, value``."""
    r = _open(raw_path)
    try:
        start, end = r.RTFromScanNum(r.FirstSpectrumNumber), r.RTFromScanNum(r.LastSpectrumNumber)
        rows = []
        for rt in np.arange(start, end, step_min):
            try:
                entry = r.GetStatusLogForRetentionTime(float(rt))
            except Exception:
                continue
            for key, val in entry.items():
                v = str(val).strip().split(" ")[0]
                if _NUM.match(v):
                    rows.append((float(rt), str(key).strip().rstrip(":"), float(v)))
    finally:
        r.Close()
    return pd.DataFrame(rows, columns=["rt", "key", "value"])


def device_traces(raw_path: str | os.PathLike) -> pd.DataFrame:
    """Every non-MS device trace: long table ``device, name, index, rt, value``.

    Thermo stores LC (pump pressure, flow), UV and other analog channels as extra
    "instruments" in the raw file; each point's value is the scan-stats TIC field.
    """
    from ThermoFisher.CommonCore.Data.Business import Device

    r = _open(raw_path)
    rows = []
    try:
        for kind in ("Analog", "MSAnalog", "Pda", "UV", "Other"):
            dev = getattr(Device, kind, None)
            if dev is None:
                continue
            n = r.source.GetInstrumentCountOfType(dev)
            for i in range(1, n + 1):
                try:
                    r.source.SelectInstrument(dev, i)
                    info = r.source.GetInstrumentData()
                    name = str(info.Name or kind)
                    labels = [str(c) for c in (info.ChannelLabels or []) if str(c)]
                    label = labels[0] if len(labels) == 1 else name
                    h = r.source.RunHeaderEx
                    for j in range(h.FirstSpectrum, h.LastSpectrum + 1):
                        rows.append(
                            (
                                kind,
                                label,
                                i,
                                r.source.RetentionTimeFromScanNumber(j),
                                r.source.GetScanStatsForScanNumber(j).TIC,
                            )
                        )
                except Exception:
                    continue
    finally:
        try:
            r.source.SelectInstrument(Device.MS, 1)
        finally:
            r.Close()
    return pd.DataFrame(rows, columns=["device", "name", "index", "rt", "value"])


def run_info(raw_path: str | os.PathLike) -> dict:
    """Acquisition metadata from the raw-file header: when, on which instrument, how big."""
    r = _open(raw_path)
    info = {"file_bytes": os.path.getsize(raw_path)}
    try:
        info["acquired"] = r.GetCreationDate()
        try:
            d = r.source.GetInstrumentData()
            info |= {
                "instrument_name": str(d.Name or ""),
                "instrument_model": str(d.Model or ""),
                "instrument_serial": str(d.SerialNumber or ""),
                "instrument_software": str(d.SoftwareVersion or ""),
            }
        except Exception as exc:  # older files: keep the date even without instrument data
            info["instrument_error"] = type(exc).__name__
    finally:
        r.Close()
    return info


def chromatography(raw_path: str | os.PathLike) -> dict:
    """All three tables plus a summary of what the file contained."""
    scans = scan_traces(raw_path)
    status = status_log(raw_path)
    devices = device_traces(raw_path)
    found = {
        "info": run_info(raw_path),
        "n_scans": len(scans),
        "ms_orders": sorted(int(o) for o in scans["ms_order"].unique()),
        "has_injection_time": bool(scans["inject_ms"].notna().any()),
        "status_keys": sorted(status["key"].unique().tolist()),
        "devices": sorted(
            {f"{d}:{n}" for d, n in zip(devices["device"], devices["name"], strict=True)}
        ),
    }
    return {"scans": scans, "status": status, "devices": devices, "found": found}


# ----------------------------------------------------------------------------------------------
# reading the extracted traces, and per-run diagnostics (pure pandas; no raw files needed)
# ----------------------------------------------------------------------------------------------


def load_chromatography(directory: str | os.PathLike) -> dict:
    """The tables written by ``scripts/raw_chromatography_modal.py`` (missing ones are None)."""
    import json
    from pathlib import Path

    d = Path(directory)
    out = {}
    for name in ("scans", "status", "devices"):
        f = d / f"{name}.parquet"
        out[name] = pd.read_parquet(f) if f.exists() else None
    found = d / "found.json"
    out["found"] = json.loads(found.read_text()) if found.exists() else {}
    return out


def binned_trace(df: pd.DataFrame, value: str, *, bins: int = 1500) -> pd.DataFrame:
    """Mean of ``value`` in ``bins`` equal RT bins: a light curve for plotting (rt, value)."""
    d = df.dropna(subset=[value])
    if d.empty:
        return pd.DataFrame(columns=["rt", value])
    edges = np.linspace(d["rt"].min(), d["rt"].max(), bins + 1)
    idx = np.clip(np.searchsorted(edges, d["rt"], side="right") - 1, 0, bins - 1)
    g = d.groupby(idx)[value].mean()
    return pd.DataFrame({"rt": (edges[g.index] + edges[g.index + 1]) / 2, value: g.to_numpy()})


def _rolling_median(t: np.ndarray, y: np.ndarray, width_min: float) -> np.ndarray:
    """Centred rolling median over a time window (scans are not evenly spaced in RT)."""
    s = pd.Series(y, index=pd.to_timedelta(t, unit="min"))
    fwd = s.rolling(pd.Timedelta(minutes=width_min / 2), min_periods=1).median().to_numpy()
    bwd = s[::-1].rolling(pd.Timedelta(minutes=width_min / 2), min_periods=1).median()
    return np.nanmean([fwd, bwd.to_numpy()[::-1]], axis=0)


def tic_metrics(scans: pd.DataFrame) -> pd.DataFrame:
    """Per run, from MS1 scans: TIC level, elution window and stability.

    ``elution_start``/``elution_end``: where the smoothed TIC is above 20 % of its maximum.
    Within that window: ``tic_median``, ``tic_area`` (TIC x min), ``tic_centroid_rt`` (signal-
    weighted mean RT), ``roughness`` (median |log10 change| between consecutive MS1 scans:
    spray flicker shows here first), ``dropouts`` (MS1 scans below 30 % of the local
    1-minute median) and ``dropout_min``, and ``inject_ms`` (median MS1 fill time).
    """
    rows = []
    for run, d in scans[scans["ms_order"] == 1].groupby("run", sort=False):
        d = d.sort_values("rt")
        t, y = d["rt"].to_numpy(), d["tic"].to_numpy(dtype=float)
        smooth = _rolling_median(t, y, 0.5)
        on = smooth >= 0.2 * np.nanmax(smooth)
        if not on.any():
            continue
        lo, hi = t[on].min(), t[on].max()
        w = (t >= lo) & (t <= hi)
        tw, yw = t[w], y[w]
        dt = np.diff(tw, prepend=tw[0])
        local = _rolling_median(tw, yw, 1.0)
        drop = yw < 0.3 * local
        logy = np.log10(np.clip(yw, 1, None))
        rows.append(
            {
                "run": run,
                "elution_start": lo,
                "elution_end": hi,
                "tic_median": float(np.median(yw)),
                "tic_area": float(np.sum(yw * dt)),
                "tic_centroid_rt": float(np.sum(tw * yw) / np.sum(yw)),
                "roughness": float(np.median(np.abs(np.diff(logy)))) if len(logy) > 2 else np.nan,
                "dropouts": int(drop.sum()),
                "dropout_min": float(np.sum(dt[drop])),
                "inject_ms": float(d["inject_ms"].median()) if "inject_ms" in d else np.nan,
            }
        )
    return pd.DataFrame(rows).set_index("run") if rows else pd.DataFrame()


def scan_metrics(scans: pd.DataFrame) -> pd.DataFrame:
    """Per run: scan counts, MS1 cycle time, MS2 fill time and median TIC per MS level."""
    rows = {}
    for run, d in scans.groupby("run", sort=False):
        ms1, ms2 = d[d["ms_order"] == 1], d[d["ms_order"] == 2]
        cycle = np.diff(np.sort(ms1["rt"].to_numpy())) * 60
        rows[run] = {
            "ms1_scans": len(ms1),
            "ms2_scans": len(ms2),
            "cycle_time_s": float(np.median(cycle)) if cycle.size else np.nan,
            "ms1_tic_median_all": float(ms1["tic"].median()) if len(ms1) else np.nan,
            "ms2_tic_median": float(ms2["tic"].median()) if len(ms2) else np.nan,
            "ms2_inject_ms": float(ms2["inject_ms"].median()) if len(ms2) else np.nan,
        }
    return pd.DataFrame.from_dict(rows, orient="index")


_SPRAY = re.compile(r"^spray current \(", re.IGNORECASE)  # "Spray Current (µA)"
# Status-log entries worth plotting for diagnosis (the log has ~200, mostly electronics).
DIAGNOSTIC_STATUS = re.compile(
    r"^(spray current|spray voltage|ion transfer tube temp|vaporizer temp|sheath gas pressure|"
    r"aux gas pressure|sweep gas pressure)",
    re.IGNORECASE,
)
_PRESSURE = re.compile(r"pressure", re.IGNORECASE)


def trace_metrics(
    status: pd.DataFrame | None,
    devices: pd.DataFrame | None,
    windows: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Per run: spray noise (status log) and LC pressure deviation from the batch median.

    ``spray_noise_pct``: short-term flicker of the spray current within the elution window
    (``windows``: :func:`tic_metrics` output) -- median absolute deviation from a 2-minute
    rolling median, as % of the current. The current falls steadily as the gradient's organic
    share rises, so its plain CV would measure the gradient, not the spray; and the logged
    "Spray Current std. dev." is 0 throughout on the Astral, so it is not used. ``pressure_dev_pct``: the largest gap between a run's
    pump-pressure trace and the median trace of all runs, as % of the median pressure (a clog
    rises, a leak falls). Columns are absent when the raw files carry neither.
    """
    out = pd.DataFrame()
    if status is not None and not status.empty:
        cur = status[status["key"].str.contains(_SPRAY)]
        if not cur.empty:
            noise = {}
            for run, d in cur.sort_values("rt").groupby("run"):
                if windows is not None and run in windows.index:
                    w = windows.loc[run]
                    d = d[(d["rt"] >= w["elution_start"]) & (d["rt"] <= w["elution_end"])]
                d = d[d["value"] > 0]
                if len(d) < 10:
                    continue
                t, y = d["rt"].to_numpy(), d["value"].to_numpy(dtype=float)
                trend = _rolling_median(t, y, 2.0)
                noise[run] = 100 * float(np.median(np.abs(y - trend))) / float(np.median(y))
            if noise:
                out["spray_noise_pct"] = pd.Series(noise)
    if devices is not None and not devices.empty:
        pr = devices[devices["name"].str.contains(_PRESSURE)]
        if not pr.empty:
            name = pr.groupby("name").size().idxmax()  # the main pump channel
            pr = pr[pr["name"] == name]
            grid = np.linspace(pr["rt"].min(), pr["rt"].max(), 400)
            curves = {
                run: np.interp(grid, d["rt"].to_numpy(), d["value"].to_numpy())
                for run, d in pr.sort_values("rt").groupby("run")
            }
            med = np.median(np.vstack(list(curves.values())), axis=0)
            scale = np.median(np.abs(med)) or 1.0
            out = out.join(
                pd.DataFrame(
                    {
                        "pressure_dev_pct": {
                            r: 100 * float(np.max(np.abs(c - med))) / scale
                            for r, c in curves.items()
                        },
                        "pump_pressure_max": pr.groupby("run")["value"].max(),
                    }
                ),
                how="outer",
            )
    return out


# Heuristic flags. Each compares a run with the median of the other runs, so it needs at least
# three runs, and says what pattern it saw and what usually causes it. They are prompts to look
# at the traces, not diagnoses: record the confirmed cause in the run table to build history.
_RULES = [
    (
        "spray",
        lambda m, ref: (
            m["roughness"] > 2 * ref["roughness"] or m["dropouts"] > max(3 * ref["dropouts"], 20)
        ),
        "Spiky TIC or signal dropouts: unstable electrospray (emitter wear or clog, bubbles, "
        "spray voltage/position). Check the TIC trace for flicker and drops.",
    ),
    (
        "spray",
        lambda m, ref: m.get("spray_noise_pct", np.nan) > 2 * ref.get("spray_noise_pct", np.inf),
        "Spray current flickers more than in the other runs: emitter wear, a partial clog or "
        "bubbles.",
    ),
    (
        "column / LC",
        lambda m, ref: (
            m.get("fwhm_s", np.nan) > 1.25 * ref.get("fwhm_s", np.inf)
            or m.get("tailing", np.nan) > 1.2 * ref.get("tailing", np.inf)
        ),
        "Broader or tailing peaks with a normal TIC: column ageing or void, dead volume in a "
        "fitting, or column temperature.",
    ),
    (
        "column / LC",
        lambda m, ref: (
            abs(m.get("rt_shift_min", 0)) > 0.5
            or abs(m["tic_centroid_rt"] - ref["tic_centroid_rt"]) > 1.0
        ),
        "Elution shifted against the other runs: gradient delay, pump flow, trap/column "
        "equilibration or a leak. Compare the pressure trace if the file has one.",
    ),
    (
        "column / LC",
        lambda m, ref: m.get("pressure_dev_pct", np.nan) > 10,
        "LC pressure departs from the other runs by more than 10 %: rising suggests a clog, "
        "falling a leak.",
    ),
    (
        "sample / load",
        lambda m, ref: (
            m["tic_area"] < 0.5 * ref["tic_area"] and m["roughness"] <= 1.5 * ref["roughness"]
        ),
        "Much less signal with a smooth TIC: less material reached the source (injection volume, "
        "sample loss, autosampler pickup) rather than a spray problem.",
    ),
]


def run_diagnostics(metrics: pd.DataFrame) -> pd.DataFrame:
    """Add ``flags`` (categories) and ``hints`` (plain-language reasons) to per-run metrics.

    ``metrics``: one row per run, from :func:`tic_metrics` joined with :func:`trace_metrics`
    and, optionally, peak metrics (``fwhm_s``, ``tailing``, ``rt_shift_min``). Each run is
    compared with the median of the other runs, so fewer than three runs get no flags.
    """
    out = metrics.copy()
    flags, hints = [], []
    for run in out.index:
        others = out.drop(index=run)
        if len(others) < 2:
            flags.append("")
            hints.append("")
            continue
        ref = others.median(numeric_only=True)
        m = out.loc[run]
        hit = [(cat, text) for cat, rule, text in _RULES if bool(rule(m, ref))]
        flags.append(", ".join(dict.fromkeys(c for c, _ in hit)))
        hints.append(" ".join(t for _, t in hit))
    out["flags"] = flags
    out["hints"] = hints
    return out
