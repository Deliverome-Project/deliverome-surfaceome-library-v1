"""Extract peak traces (XICs) for chosen peptides straight from Thermo .raw files.

The backup route to peak-level data when Spectronaut's own XIC database is not available:
for each target precursor, read only the DIA MS2 scans whose isolation window contains it,
inside a retention-time window around Spectronaut's apex, and pull the intensity of its b/y
fragment ions (plus the MS1 isotope trace). Output is one long table that
:func:`plot_precursor_xics` draws.

Reading .raw needs Thermo's RawFileReader through ``alpharaw`` + pythonnet + mono/.NET, which
the repo does not install; ``scripts/raw_xic_modal.py`` runs this on Modal next to S3, since
an Astral raw file is ~18 GB. The mass and target functions here are pure Python and tested
offline.

Masses: monoisotopic residue masses and the three modifications Spectronaut reports in our
searches. ``precursor_mz`` reproduces Spectronaut's ``FG.PrecMz`` (checked on the Astral export),
which is the guard against a wrong constant.
"""

from __future__ import annotations

import os
import re
from pathlib import Path

import numpy as np
import pandas as pd

PROTON = 1.007276467
WATER = 18.010564684
# Monoisotopic residue masses (Unimod / NIST).
RESIDUE = {
    "G": 57.021464,
    "A": 71.037114,
    "S": 87.032028,
    "P": 97.052764,
    "V": 99.068414,
    "T": 101.047679,
    "C": 103.009185,
    "L": 113.084064,
    "I": 113.084064,
    "N": 114.042927,
    "D": 115.026943,
    "Q": 128.058578,
    "K": 128.094963,
    "E": 129.042593,
    "M": 131.040485,
    "H": 137.058912,
    "F": 147.068414,
    "R": 156.101111,
    "Y": 163.063329,
    "W": 186.079313,
    "U": 150.953636,
}
# Spectronaut modification names -> monoisotopic delta (Unimod 4, 35, 1).
MODS = {
    "Carbamidomethyl (C)": 57.021464,
    "Oxidation (M)": 15.994915,
    "Acetyl (Protein N-term)": 42.010565,
}
C13 = 1.003355  # isotope spacing


def parse_modseq(modseq: str) -> tuple[str, list[float]]:
    """``_[Acetyl (Protein N-term)]MC[Carbamidomethyl (C)]K_`` -> ("MCK", per-residue deltas).

    An N-terminal modification is added to the first residue. Unknown modification names
    raise, so a new one is added to ``MODS`` deliberately rather than silently ignored.
    """
    s = modseq.strip("_")
    seq, deltas, pending = [], [], 0.0
    for tok in re.finditer(r"\[([^\]]+)\]|([A-Z])", s):
        name, aa = tok.groups()
        if name is not None:
            if name not in MODS:
                raise ValueError(f"Unknown modification {name!r} in {modseq}; add it to MODS")
            if seq:
                deltas[-1] += MODS[name]
            else:
                pending += MODS[name]
        else:
            if aa not in RESIDUE:
                raise ValueError(f"Unknown residue {aa!r} in {modseq}")
            seq.append(aa)
            deltas.append(pending)
            pending = 0.0
    return "".join(seq), deltas


def _residue_masses(modseq: str) -> np.ndarray:
    seq, deltas = parse_modseq(modseq)
    return np.array([RESIDUE[a] for a in seq]) + np.array(deltas)


def precursor_mz(modseq: str, charge: int) -> float:
    return (_residue_masses(modseq).sum() + WATER + charge * PROTON) / charge


def fragment_ions(
    modseq: str, precursor_charge: int, *, min_len: int = 2, max_frag_charge: int = 2
) -> pd.DataFrame:
    """b and y ions (charge 1..min(max_frag_charge, precursor_charge - 1 or 1)) with their m/z."""
    m = _residue_masses(modseq)
    n = len(m)
    rows = []
    for z in range(1, max(1, min(max_frag_charge, precursor_charge - 1)) + 1):
        for i in range(min_len, n):
            b = (m[:i].sum() + z * PROTON) / z
            y = (m[n - i :].sum() + WATER + z * PROTON) / z
            sfx = "" if z == 1 else f"^{z}"
            rows += [(f"b{i}{sfx}", "b", i, z, b), (f"y{i}{sfx}", "y", i, z, y)]
    return pd.DataFrame(rows, columns=["ion", "type", "number", "charge", "mz"])


def split_precursor(precursor: str) -> tuple[str, int]:
    """limpa / trimmer key ``<EG.ModifiedSequence>.<FG.Charge>`` -> (modseq, charge)."""
    modseq, charge = precursor.rsplit(".", 1)
    return modseq, int(charge)


def select_targets(
    matrix_dir: str | os.PathLike,
    protein_ids: list[str],
    *,
    per_protein: int = 4,
) -> pd.DataFrame:
    """Precursors to extract for each protein: the most intense ones seen in every run.

    Reads the trimmer's small tables. Returns precursor, protein, modseq, charge, mz and
    one apex RT per run (columns ``apex_rt:<run>``).
    """
    d = Path(matrix_dir)
    qc = pd.read_parquet(
        d / "precursor_qc.parquet", columns=["precursor", "run", "log2_intensity", "apex_rt"]
    )
    ann = pd.read_parquet(d / "precursors.parquet").set_index("precursor")
    qc["protein"] = qc["precursor"].map(ann["protein_id"])
    qc = qc[qc["protein"].isin(protein_ids)]
    n_runs = qc["run"].nunique()
    g = qc.groupby(["protein", "precursor"])
    stats = pd.DataFrame({"n": g["run"].nunique(), "mean": g["log2_intensity"].mean()})
    stats = stats[stats["n"] == n_runs].reset_index()
    stats = stats.sort_values(["protein", "mean"], ascending=[True, False])
    stats = stats.groupby("protein").head(per_protein)
    apex = qc.pivot_table(index="precursor", columns="run", values="apex_rt", aggfunc="first")
    apex.columns = [f"apex_rt:{c}" for c in apex.columns]
    out = stats.merge(apex, left_on="precursor", right_index=True)
    parts = out["precursor"].map(split_precursor)
    out["modseq"] = [p[0] for p in parts]
    out["charge"] = [p[1] for p in parts]
    out["mz"] = [precursor_mz(m, z) for m, z in parts]
    missing = sorted(set(protein_ids) - set(out["protein"]))
    if missing:
        raise ValueError(f"No precursor seen in every run for: {missing}")
    return out.reset_index(drop=True)


# ----------------------------------------------------------------------------------------------
# raw-file access (needs alpharaw + pythonnet + mono/.NET; runs on Modal)
# ----------------------------------------------------------------------------------------------


class _Rec(dict):
    __getattr__ = dict.__getitem__  # t.mz and t["apex_rt:<run>"] both work


def _records(df: pd.DataFrame) -> list[_Rec]:
    return [_Rec(r) for r in df.to_dict("records")]


def scan_index(raw_path: str | os.PathLike) -> pd.DataFrame:
    """One row per scan: scan, rt (min), ms_order, isolation window low/high (MS2)."""
    from alpharaw.raw_access.pythermorawfilereader import RawFileReader

    r = RawFileReader(str(raw_path))
    try:
        rows = []
        for s in range(r.FirstSpectrumNumber, r.LastSpectrumNumber + 1):
            order = r.GetMSOrderForScanNum(s)
            if order == 2:
                c = r.GetPrecursorMassForScanNum(s, 0)
                w = r.GetIsolationWidthForScanNum(s, 0)
                lo, hi = c - w / 2, c + w / 2
            else:
                lo = hi = np.nan
            rows.append((s, r.RTFromScanNum(s), order, lo, hi))
    finally:
        r.Close()
    return pd.DataFrame(rows, columns=["scan", "rt", "ms_order", "iso_low", "iso_high"])


def _intensity_at(mz: np.ndarray, inten: np.ndarray, targets: np.ndarray, ppm: float):
    """Max centroid intensity within +/- ppm of each target m/z (0 if none); mz sorted."""
    tol = targets * ppm * 1e-6
    lo = np.searchsorted(mz, targets - tol, "left")
    hi = np.searchsorted(mz, targets + tol, "right")
    return np.array([inten[a:b].max() if b > a else 0.0 for a, b in zip(lo, hi, strict=True)])


def extract_xics(
    raw_path: str | os.PathLike,
    targets: pd.DataFrame,
    run: str,
    *,
    index: pd.DataFrame | None = None,
    rt_window: float = 0.5,
    ppm: float = 10.0,
    n_isotopes: int = 3,
) -> pd.DataFrame:
    """Fragment (MS2) and isotope (MS1) traces for each target precursor in one raw file.

    ``targets`` comes from :func:`select_targets`; its ``apex_rt:<run>`` column centres the
    window. Returns precursor, ms_level, ion, ion_mz, rt, intensity.
    """
    from alpharaw.raw_access.pythermorawfilereader import RawFileReader

    idx = index if index is not None else scan_index(raw_path)
    ms1 = idx[idx["ms_order"] == 1]
    ms2 = idx[idx["ms_order"] == 2]
    r = RawFileReader(str(raw_path))
    out = []
    try:
        for t in _records(targets):
            apex = t[f"apex_rt:{run}"]
            win = (idx["rt"] >= apex - rt_window) & (idx["rt"] <= apex + rt_window)
            frags = fragment_ions(t.modseq, t.charge)
            sel2 = ms2[win.loc[ms2.index] & (ms2["iso_low"] <= t.mz) & (ms2["iso_high"] >= t.mz)]
            for s, rt in zip(sel2["scan"], sel2["rt"], strict=True):
                mz, inten = (np.asarray(a) for a in r.GetCentroidMassListFromScanNum(int(s)))
                o = np.argsort(mz)
                vals = _intensity_at(mz[o], inten[o], frags["mz"].to_numpy(), ppm)
                out.append(
                    pd.DataFrame(
                        {
                            "precursor": t.precursor,
                            "ms_level": 2,
                            "ion": frags["ion"],
                            "ion_mz": frags["mz"],
                            "rt": rt,
                            "intensity": vals,
                        }
                    )
                )
            iso = t.mz + np.arange(n_isotopes) * C13 / t.charge
            sel1 = ms1[win.loc[ms1.index]]
            for s, rt in zip(sel1["scan"], sel1["rt"], strict=True):
                mz, inten = (np.asarray(a) for a in r.GetCentroidMassListFromScanNum(int(s)))
                o = np.argsort(mz)
                vals = _intensity_at(mz[o], inten[o], iso, ppm)
                out.append(
                    pd.DataFrame(
                        {
                            "precursor": t.precursor,
                            "ms_level": 1,
                            "ion": [f"M+{k}" for k in range(n_isotopes)],
                            "ion_mz": iso,
                            "rt": rt,
                            "intensity": vals,
                        }
                    )
                )
    finally:
        r.Close()
    res = (
        pd.concat(out, ignore_index=True)
        if out
        else pd.DataFrame(columns=["precursor", "ms_level", "ion", "ion_mz", "rt", "intensity"])
    )
    res.insert(0, "run", run)
    return res


# ----------------------------------------------------------------------------------------------
# plotting
# ----------------------------------------------------------------------------------------------


def top_fragments(xic: pd.DataFrame, n: int = 6) -> pd.DataFrame:
    """Keep each precursor's n fragments with the most summed MS2 signal across runs."""
    ms2 = xic[xic["ms_level"] == 2]
    tot = ms2.groupby(["precursor", "ion"])["intensity"].sum()
    keep = tot[tot > 0].groupby(level=0, group_keys=False).nlargest(n).index
    return ms2.set_index(["precursor", "ion"]).loc[keep].reset_index()


def plot_precursor_xics(
    xic: pd.DataFrame,
    targets: pd.DataFrame,
    protein: str,
    *,
    title: str | None = None,
    n_fragments: int = 6,
):
    """Grid for one protein: rows = its precursors, columns = runs; top fragment traces with
    Spectronaut's apex marked. Call ``setup_plotting_style()`` first."""
    import matplotlib.pyplot as plt

    from deliverome_analysis.dia_qc import run_labels
    from deliverome_analysis.plotting import CATEGORICAL_PALETTE

    tg = targets[targets["protein"] == protein]
    runs = sorted(xic["run"].unique())
    labels = run_labels(runs)
    frag = top_fragments(xic[xic["precursor"].isin(tg["precursor"])], n_fragments)
    fig, axes = plt.subplots(
        len(tg),
        len(runs),
        figsize=(4.2 * len(runs), 3.0 * len(tg)),
        squeeze=False,
        sharex="row",
        sharey="row",
    )
    for i, t in enumerate(_records(tg)):
        ions = list(dict.fromkeys(frag.loc[frag["precursor"] == t.precursor, "ion"]))
        colors = {
            ion: CATEGORICAL_PALETTE[k % len(CATEGORICAL_PALETTE)] for k, ion in enumerate(ions)
        }
        for j, run in enumerate(runs):
            ax = axes[i, j]
            d = frag[(frag["precursor"] == t.precursor) & (frag["run"] == run)]
            for ion in ions:
                e = d[d["ion"] == ion].sort_values("rt")
                ax.plot(e["rt"], e["intensity"], color=colors[ion], lw=1.6, label=ion)
            apex = t.get(f"apex_rt:{run}")
            if apex is not None:
                ax.axvline(apex, color="black", ls="--", lw=0.8)
            if i == 0:
                ax.text(
                    0.5,
                    1.02,
                    labels[run],
                    transform=ax.transAxes,
                    ha="center",
                    va="bottom",
                    fontsize="small",
                )
            if j == 0:
                ax.set_ylabel(f"{t.modseq.strip('_')}\n{t.charge}+", fontsize="x-small")
                ax.ticklabel_format(axis="y", style="sci", scilimits=(-3, 4))
            if i == len(tg) - 1:
                ax.set_xlabel("RT (min)")
        axes[i, -1].legend(frameon=False, fontsize="xx-small", loc="upper right")
    if title:
        # setup_plotting_style() disables set_title/suptitle by project convention; a header is
        # drawn as figure text only when the caller asks for one (e.g. figures shared outside
        # a captioned document).
        fig.tight_layout(rect=(0, 0, 1, 1 - 0.35 / fig.get_figheight()))
        fig.text(
            0.5,
            1 - 0.12 / fig.get_figheight(),
            title,
            ha="center",
            va="top",
            fontsize="large",
            fontweight="bold",
        )
    else:
        fig.tight_layout()
    return fig
