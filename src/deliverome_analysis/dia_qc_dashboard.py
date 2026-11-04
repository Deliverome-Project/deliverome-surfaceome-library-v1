"""Interactive DIA QC dashboard: one self-contained HTML page per YAML layout.

A layout (``dia_qc_layouts/*.yaml``, or your own file) names the panels and their order; this
module turns it into an HTML page of interactive Plotly charts (hover, zoom, toggle runs) and
opens it in the browser. Layouts are fixed: a panel whose columns the export lacks stays in
place and says what is missing, so the page doubles as a coverage map of what this export
can and cannot answer.

    from deliverome_analysis.dia_qc_dashboard import dia_qc_dashboard
    path, coverage = dia_qc_dashboard(matrix_dir="trimmed", layout="default", output="qc.html")

or from the shell:

    uv run python -m deliverome_analysis.dia_qc_dashboard trimmed --layout instrument

Layouts shipped: ``default`` (everything), ``instrument`` (chromatography, mass accuracy, peak
sampling) and ``sample`` (depth, reproducibility, digestion, identification evidence).

Only per-run summaries are embedded (quantiles, histograms, counts), never the precursor
table, so a page stays small at any run count. Plotly.js loads from jsDelivr, so viewing
needs a network connection. Static per-check PNGs for papers and Notion still come from
:func:`deliverome_analysis.dia_qc.dia_qc_report`.
"""

from __future__ import annotations

import html
import json
import math
import os
import webbrowser
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from deliverome_analysis import dia_qc as qc

LAYOUT_DIR = Path(__file__).parent / "dia_qc_layouts"
PLOTLY_JS = "https://cdn.jsdelivr.net/npm/plotly.js-dist-min@2.35.2/plotly.min.js"
GRID = "rgba(31,23,24,0.07)"
# Dark, high-contrast hover labels (ink background, white text) on every chart.
HOVER = {
    "bgcolor": "#1f1718",
    "bordercolor": "#1f1718",
    "font": {"family": "Manrope, Inter, sans-serif", "size": 12.5, "color": "#ffffff"},
    "align": "left",
}


@dataclass(frozen=True)
class Metric:
    kind: str  # ids | dist | completeness | detection | cv | gradient | missed | planned
    title: str
    column: str | None = None
    label: str = ""
    centred: bool = False
    log: bool = False
    note: str = ""  # one-line reading guide shown under the chart
    needs: str = ""  # what the export must carry for this panel


_LEAN = "the lean report template (docs/spectronaut-lean-report.md)"

METRICS: dict[str, Metric] = {
    "protein_explorer": Metric(
        "explorer",
        "Protein explorer",
        note="Pick a run: points are proteins, placed by abundance and by how far this run sits "
        "from the other runs (its replicate group when there is one). Click a point or a row to "
        "see that protein's precursors in every run.",
    ),
    "ids_precursors": Metric("ids", "Precursors per run", "precursors", "Precursors (1% FDR)"),
    "ids_peptides": Metric("ids", "Peptides per run", "peptides", "Peptides (1% FDR)"),
    "ids_protein_groups": Metric(
        "ids", "Protein groups per run", "protein_groups", "Protein groups (1% FDR)"
    ),
    "completeness_runs": Metric(
        "compl_runs",
        "Data completeness per run",
        note="Each point is a run: the share of everything identified in the batch that this "
        "run detected. Boxes summarise the runs of a sample group; colour is the group, hover "
        "for the run, click it for details.",
    ),
    "completeness": Metric(
        "completeness",
        "Detection overlap",
        note="Share of precursors detected in k of n runs.",
    ),
    "charge_states": Metric(
        "charge",
        "Precursor charge states",
        note="Share of each run's identified precursors at each charge. A shift towards 1+ or "
        "towards high charges between runs points at spray conditions or sample cleanup.",
    ),
    "surface_protein_groups": Metric(
        "surf_ids",
        "Surface protein groups",
        "protein_groups",
        "Protein groups",
        note="Protein groups by accessible-surfaceome tier (deep green: canonical, light green: "
        "likely, light yellow: low / uncertain). A group takes "
        "the best tier of any of its proteins. Click a bar to list its proteins.",
        needs="the accessible-surfaceome catalog (S3) and PG.ProteinAccessions",
    ),
    "surface_peptides": Metric(
        "surf_ids",
        "Surface peptides",
        "peptides",
        "Peptides",
        note="Peptides of surface protein groups, by tier.",
        needs="the accessible-surfaceome catalog (S3) and PG.ProteinAccessions",
    ),
    "surface_signal": Metric(
        "surf_signal",
        "Surface share of signal",
        note="Share of each run's summed precursor signal from surface protein groups, by tier, "
        "for MS2 and (when exported) MS1. A few abundant proteins dominate: click a bar.",
        needs="the accessible-surfaceome catalog; FG.MS1Quantity for the MS1 share",
    ),
    "surface_missed": Metric(
        "surf_missed",
        "Missed cleavages by tier",
        note="Precursors with at least one missed cleavage, per tier; the line is all precursors. "
        "Surface proteins digest differently (glycosylation, membrane-embedded segments).",
        needs="PEP.NrOfMissedCleavages and the accessible-surfaceome catalog",
    ),
    "surface_points_per_peak": Metric(
        "surf_points",
        "Points per peak by tier (MS2)",
        "points_per_peak",
        "Data points",
        note="Sampling of each elution peak, by tier. Below ~6 points quantification degrades.",
        needs=_LEAN,
    ),
    "tic_ms1": Metric(
        "tic",
        "TIC (MS1)",
        "1",
        "Total ion current",
        note="MS1 total ion current against retention time, one line per run. Look for "
        "dropouts and spikes (spray), a shifted or compressed profile (LC, column) and the "
        "overall level (load). Drag to zoom, double-click to reset; click a line for that run.",
        needs="the raw files (scripts/raw_chromatography_modal.py)",
    ),
    "tic_ms2": Metric(
        "tic",
        "TIC (MS2)",
        "2",
        "Total ion current",
        note="Fragment-scan total ion current; follows MS1 when the instrument is healthy.",
        needs="the raw files (scripts/raw_chromatography_modal.py)",
    ),
    "injection_time": Metric(
        "inject",
        "MS1 injection time",
        note="Ion fill time per MS1 scan. It rises when less ion current arrives (weaker spray, "
        "less sample) and caps at the method maximum.",
        needs="the raw files (scripts/raw_chromatography_modal.py)",
    ),
    "instrument_traces": Metric(
        "traces",
        "Instrument and LC traces",
        note="Whatever the raw files logged: spray voltage/current, source temperatures, LC pump "
        "pressure. Pick a trace from the menu.",
        needs="raw files whose status log or LC device records them",
    ),
    "qc_history": Metric(
        "history",
        "Runs over time",
        note="This page's runs (large points) against every earlier run in the shared QC index, "
        "by acquisition date. Pick a metric from the menu. Slow drifts (column ageing, falling "
        "sensitivity) show here before they show within one batch.",
        needs="the shared QC index (python -m deliverome_analysis.qcdb build dia_qc)",
    ),
    "chromatography_diagnosis": Metric(
        "diagnosis",
        "Run diagnosis",
        note="Per-run chromatography metrics with heuristic flags. Each run is compared with the "
        "median of the others (needs 3+ runs). Flags are prompts to look at the traces, not "
        "verdicts: confirm the cause and record it in run_diagnostics.tsv to build history.",
        needs="the raw files (scripts/raw_chromatography_modal.py)",
    ),
    "detection_reproducibility": Metric(
        "detection",
        "Detection reproducibility",
        note="Chance a precursor seen at this intensity in one run is seen in the others.",
    ),
    "intensity": Metric("dist", "Precursor intensity", "log2_intensity", "log2 intensity"),
    "precursor_cv": Metric(
        "cv",
        "Replicate CV",
        note="Within replicate groups only; needs at least 2 runs per group.",
    ),
    "rt_error": Metric(
        "dist",
        "Retention-time error",
        "rt_error",
        "Apex RT - predicted RT (min)",
        centred=True,
        needs="EG.ApexRT and EG.RTPredicted",
    ),
    "irt_residual": Metric(
        "dist",
        "iRT recalibration residual",
        "irt_residual",
        "Empirical - predicted iRT",
        centred=True,
        note="Anchored on endogenous library peptides unless iRT standards were searched.",
        needs="EG.iRTEmpirical and EG.iRTPredicted",
    ),
    "rt_consistency": Metric(
        "rtcons",
        "Retention time vs other runs",
        note="Each run's shift across the gradient (line: median, band: 5th-95th percentile) "
        "relative to the median of all runs. A flat line near zero is a stable LC; a lone "
        "peptide far from its run's line is flagged as a likely wrong peak.",
        needs="EG.ApexRT and >= 2 runs",
    ),
    "charge_state_rt": Metric(
        "chargert",
        "Charge-state co-elution",
        note="Spread in apex RT between charge states of the same modified peptide in one run. "
        "They should co-elute; spreads above 6 s are flagged.",
        needs="EG.ApexRT",
    ),
    "rt_flagged": Metric(
        "rtflag",
        "RT outliers",
        note="Precursors whose RT disagrees with the other runs after removing the run's own "
        "shift (|residual| above 5 robust SD and at least 6 s). Click a row for the protein's "
        "precursors.",
        needs="EG.ApexRT and >= 2 runs",
    ),
    "dpc_fit": Metric(
        "limpa_dpc",
        "Detection probability curve (limpa)",
        note="Fraction of runs a precursor is detected in, against limpa's estimated mean "
        "intensity; the line is the fitted curve limpa uses instead of imputing.",
        needs="a limpa run (limpa_result)",
    ),
    "protein_se": Metric(
        "limpa_se",
        "Protein standard error (limpa)",
        note="limpa's per-protein uncertainty against abundance; proteins with few precursors "
        "and low intensity carry the most.",
        needs="a limpa run (limpa_result)",
    ),
    "run_correlation": Metric(
        "limpa_corr",
        "Run-to-run correlation (limpa)",
        note="Pearson r of limpa protein values between runs, on proteins detected in every "
        "run (model-filled values would dominate otherwise).",
        needs="a limpa run (limpa_result)",
    ),
    "pca": Metric(
        "limpa_pca",
        "PCA of runs (limpa)",
        note="Runs on the first two principal components of limpa protein values.",
        needs="a limpa run with at least 4 runs",
    ),
    "limpa_vs_spectronaut": Metric(
        "limpa_agree",
        "limpa vs Spectronaut",
        note="Between-run change per protein (each protein centred on its mean), limpa against "
        "Spectronaut PG.Quantity. Shading: all values; dots: the 300 protein-run values where the "
        "two disagree most, named on hover. Click a dot for the protein.",
        needs="a limpa run and PG.Quantity",
    ),
    "contaminant_signal": Metric(
        "contam",
        "Contaminant and standard signal",
        note="Share of each run's precursor signal from contaminant-only groups (left out of "
        "limpa), groups shared with a contaminant (kept, flagged), watchlist human proteins "
        "(kept, flagged) and the iRT standard (left out of limpa).",
        needs="a search against the published database (Cont_ / Biognosys_iRT entries)",
    ),
    "top_contaminants": Metric(
        "contamtable",
        "Top contaminants",
        note="Contaminant, shared and watchlist groups by median share of run signal. Click a "
        "row for details.",
        needs="a search against the published database (Cont_ / Biognosys_iRT entries)",
    ),
    "gradient_coverage": Metric("gradient", "Gradient coverage", "apex_rt", needs="EG.ApexRT"),
    "irt_standards": Metric(
        "planned",
        "iRT standards",
        note="Retention and intensity of spiked iRT/PRTC peptides across runs.",
        needs="a spiked iRT or PRTC standard, with its sequences in the search FASTA",
    ),
    "ms1_calibration_ppm": Metric(
        "dist",
        "MS1 mass correction",
        "ms1_calibration_ppm",
        "Correction applied (ppm)",
        centred=True,
        note="The calibration Spectronaut applied; drift between runs shows here.",
        needs="FG.PrecMz and FG.PrecMzCalibrated",
    ),
    "ms1_mass_error": Metric(
        "dist",
        "MS1 mass error",
        "ms1_mass_error_ppm",
        "Observed error (ppm)",
        centred=True,
        needs=_LEAN,
    ),
    "peak_gauss_score": Metric(
        "dist", "Peak shape (Gaussian fit)", "peak_gauss_score", "Score", needs=_LEAN
    ),
    "fragment_coelution": Metric(
        "dist",
        "Fragment co-elution",
        "fragment_coelution",
        "Score",
        needs=_LEAN,
        note="How tightly a precursor's fragment traces rise and fall together.",
    ),
    "ion_ladder_score": Metric(
        "dist",
        "Ion ladder (b/y)",
        "ion_ladder_score",
        "Score",
        needs=_LEAN,
        note="Spectronaut's per-precursor b/y ladder score; exact ion counts need the fragment report.",
    ),
    "ms2_mass_error": Metric(
        "planned",
        "MS2 mass error",
        note="Observed fragment mass error per run.",
        needs="the fragment-level report (F.CalibratedMassAccuracy (PPM))",
    ),
    "missed_cleavages": Metric("missed", "Missed cleavages", needs="PEP.NrOfMissedCleavages"),
    "fwhm": Metric("dist", "Peak width (FWHM)", "fwhm", "FWHM (min)", needs=_LEAN),
    "points_per_peak": Metric(
        "dist", "Points per peak (MS2)", "points_per_peak", "Data points", needs=_LEAN
    ),
    "points_per_peak_ms1": Metric(
        "dist", "Points per peak (MS1)", "points_per_peak_ms1", "Data points", needs=_LEAN
    ),
    "tailing_factor": Metric(
        "dist", "Peak tailing", "tailing_factor", "Tailing factor at FWHM", needs=_LEAN
    ),
    "signal_to_noise": Metric(
        "dist", "Signal-to-noise", "signal_to_noise", "S/N", log=True, needs=_LEAN
    ),
    "cscore": Metric("dist", "Identification score", "cscore", "Cscore", needs="EG.Cscore"),
    "fragment_coverage": Metric(
        "planned",
        "b/y ion coverage",
        note="Fragment ions observed and used per precursor, by ion type.",
        needs="the fragment-level report (F.FrgType, F.FrgNum)",
    ),
}


# Scientist-facing explanations behind each panel's (i) icon: what is measured, what healthy
# looks like, and what a problem looks like. Ranges are rules of thumb, not acceptance criteria.
INFO: dict[str, str] = {
    "protein_explorer": "Each point is a protein group. Across: how abundant it is on average. Up "
    "or down: how much higher or lower it is in the selected run than in the other runs. Most "
    "proteins should sit near zero; a run with many proteins far from zero has a global problem "
    "(load, digestion, a different sample). Colour is the surfaceome tier.",
    "ids_precursors": "Precursors (a peptide at one charge state) identified at 1% FDR in each run. "
    "The most sensitive depth measure. Replicates should agree within a few percent; a run "
    "clearly lower than its replicates lost signal somewhere (load, spray, column).",
    "ids_peptides": "Distinct peptide sequences identified at 1% FDR per run.",
    "ids_protein_groups": "Protein groups identified at 1% FDR per run. A group collects proteins "
    "that the identified peptides cannot tell apart (isoforms, close paralogs). Less sensitive "
    "than precursors: a run can lose many precursors before it loses protein groups.",
    "completeness_runs": "For each run, the share of everything identified anywhere in the batch "
    "that this run also identified. In a good replicate set every run is above ~90% at the "
    "protein level and close to its replicates; one low point is a weak run.",
    "completeness": "How many precursors were found in all runs, in all but one, and so on. A tall "
    "right-hand bar means consistent detection; a lot in the left bars means many identifications "
    "near the detection limit or one run behaving differently.",
    "charge_states": "Share of identified precursors at each charge. Tryptic peptides are mostly 2+ "
    "and 3+. More 1+ can mean poor desolvation or contaminants (polymers); a shift to high charges "
    "can follow changes in spray conditions or solvent. Runs of the same sample type should match.",
    "surface_protein_groups": "Protein groups whose protein is a cell-surface protein in the "
    "accessible-surfaceome catalog, split by confidence: canonical (curated, certain), likely, and "
    "low/uncertain. For whole-cell lysate, surface proteins are a small, low-abundance fraction, "
    "so these counts are a direct readout of how deep the run reaches into the membrane proteome.",
    "surface_peptides": "Peptides from surface protein groups, by tier. Membrane proteins give "
    "fewer tryptic peptides (transmembrane segments, glycosylated ectodomains), so expect fewer "
    "peptides per protein than for soluble proteins.",
    "surface_signal": "What fraction of each run's total precursor signal comes from surface "
    "proteins, by tier. In lysate a few very abundant proteins can dominate a tier (in HeLa, "
    "enolase and keratin 8 in 'likely'), so click a bar to see who carries the signal. Surface enrichment "
    "protocols should raise the canonical share markedly.",
    "surface_missed": "Fraction of precursors with at least one missed tryptic cleavage, per tier, "
    "against the whole run (dotted line). Membrane proteins often digest less completely. A rise "
    "across all tiers points at digestion (enzyme ratio, time, denaturant), not the instrument.",
    "surface_points_per_peak": "How many measurements describe each elution peak, by tier. "
    "Below about 6 points the peak shape, and so the quantity, becomes unreliable; low-abundance "
    "surface proteins are the first to suffer.",
    "tic_ms1": "Total ion current: everything the mass spectrometer sees in each MS1 scan, along "
    "the gradient, one line per run. Overlaid replicates should have the same shape and height. "
    "Short spikes or sudden drops point at the electrospray; a stretched, shifted or flattened "
    "profile points at the LC or column; a uniformly lower curve means less material.",
    "tic_ms2": "The same as the MS1 TIC for fragment scans. It should track MS1; if MS1 looks normal "
    "but MS2 does not, look at fragmentation or the DIA window settings.",
    "injection_time": "How long the instrument fills the ion trap for each MS1 scan to reach its "
    "target. It goes up when fewer ions arrive (weaker spray, less sample) and flattens at the "
    "method maximum. Long fill times at the start and end of the run are normal (little elutes).",
    "instrument_traces": "Signals the instrument and LC logged during the run. Spray current falls "
    "steadily as the gradient becomes more organic; flicker or drops on top of that trend mean an "
    "unstable spray. LC pump pressure should follow the same curve every run: a gradual rise "
    "suggests a clogging column or frit, a drop a leak.",
    "chromatography_diagnosis": "Summary numbers for each run's chromatography, compared with the "
    "other runs. Roughness and dropouts measure TIC flicker; spray noise measures spray-current "
    "flicker; pressure deviation compares LC pressure traces. Flags are rule-of-thumb prompts to "
    "look at the traces, not diagnoses. Record what a flagged run really was in "
    "run_diagnostics.tsv so the rules can be checked against experience.",
    "detection_reproducibility": "If a precursor is seen in one run at a given intensity, how "
    "likely is it to be seen in the other runs? Should approach 1 at high intensity and fall at "
    "low intensity; a curve that stays low points at inconsistent detection.",
    "intensity": "Distribution of precursor intensities in each run. Replicates should overlap; a "
    "whole-distribution shift is a loading or sensitivity difference.",
    "precursor_cv": "Coefficient of variation of each precursor's intensity across the runs of a "
    "replicate group. Lower is better: for technical replicates, a median below ~10-15% is good, "
    "and most precursors should be below 20%.",
    "rt_error": "Measured apex retention time minus the time Spectronaut predicted from its "
    "library. Centred on zero and narrow is good; a wide or shifted distribution means the "
    "chromatography differs from the library's.",
    "irt_residual": "How well each precursor fits the run's retention-time calibration (empirical "
    "minus predicted iRT). Narrow is good; wide means irregular chromatography.",
    "rt_consistency": "How each run's retention times differ from the other runs along the "
    "gradient. A flat line near zero is a stable LC; an offset is a shifted gradient; a slope or "
    "curve means the elution was stretched or compressed in parts of the gradient.",
    "charge_state_rt": "The same peptide at different charge states should elute at the same "
    "time. Spreads above a few seconds suggest a wrong peak pick for one of them.",
    "rt_flagged": "Precursors whose retention time disagrees with the other runs once the run's "
    "overall shift is removed: likely a different (wrong) peak was picked. Click one to inspect.",
    "dpc_fit": "limpa's detection probability curve: how likely a precursor is to be detected, "
    "given how abundant it is. limpa uses it to account for missing values instead of imputing "
    "them. A steep curve means missingness is mostly about low abundance, as expected.",
    "protein_se": "limpa's standard error for each protein's value: the uncertainty it carries "
    "forward into differential testing. Higher for low-abundance proteins with few precursors.",
    "run_correlation": "How similar the runs' protein values are (Pearson r). Technical replicates "
    "should correlate very highly (above ~0.95); a run that correlates less with all others "
    "stands apart.",
    "pca": "Runs placed by the main patterns in their protein values. Replicates should cluster; "
    "the separation should follow the experimental groups, not run order or batch.",
    "limpa_vs_spectronaut": "Whether limpa's protein fold changes between runs agree with "
    "Spectronaut's own protein quantities. Points near the diagonal mean the two agree.",
    "contaminant_signal": "Share of each run's signal from contaminant proteins (keratins from "
    "dust, trypsin, BSA from serum), from proteins on the contaminant watchlist, and from the "
    "spiked iRT standard. A run with much more contaminant signal had a handling problem.",
    "top_contaminants": "The contaminants carrying the most signal. Useful for finding the source "
    "(skin keratins, serum albumin, enzyme).",
    "gradient_coverage": "Where along the gradient peptides elute. A good method spreads them "
    "across most of the gradient; a pile-up at one end wastes instrument time.",
    "irt_standards": "Retention and intensity of spiked iRT standard peptides: the most direct "
    "check of LC stability across runs, when the standards are spiked and searched.",
    "ms1_calibration_ppm": "The mass correction Spectronaut applied to each precursor. Stable "
    "within and between runs is good; drift means the mass calibration moved during the run.",
    "ms1_mass_error": "Remaining precursor mass error after calibration. Should be centred on zero "
    "and within a few ppm.",
    "peak_gauss_score": "How closely each elution peak resembles an ideal Gaussian. Lower scores "
    "mean distorted peaks (tailing, splitting, interference).",
    "fragment_coelution": "How tightly a precursor's fragment ions rise and fall together. High "
    "means a clean identification; low suggests interference from another peptide.",
    "ion_ladder_score": "Spectronaut's score for how continuous the b/y fragment-ion series is. "
    "Higher means better-supported identifications.",
    "ms2_mass_error": "Fragment-ion mass error; needs the fragment-level report.",
    "missed_cleavages": "Share of identified precursors with 0, 1 or 2+ missed trypsin cleavages. "
    "Mostly 0 is a complete digest; a rising share of missed cleavages means incomplete digestion.",
    "fwhm": "Peak width at half height. Narrow, consistent peaks mean good chromatography; widening "
    "across runs is an early sign of column ageing or dead volume.",
    "points_per_peak": "Number of MS2 measurements across each elution peak. Below about 6 the "
    "peak, and so the quantity, is poorly defined.",
    "points_per_peak_ms1": "Number of MS1 measurements across each elution peak.",
    "tailing_factor": "Peak asymmetry: 1 is symmetric, above ~1.5 tails noticeably. Tailing that "
    "grows over runs points at the column or a fitting.",
    "signal_to_noise": "Peak height relative to background noise. Higher is better; a drop across "
    "the board means lower sensitivity.",
    "cscore": "Spectronaut's identification score. Higher means more confident identifications; "
    "compare distributions between runs rather than absolute values.",
    "fragment_coverage": "Which fragment ions (b and y) were observed per precursor; needs the "
    "fragment-level report.",
}


# ----------------------------------------------------------------------------------------------
# layouts
# ----------------------------------------------------------------------------------------------


def available_layouts() -> list[str]:
    return sorted(p.stem for p in LAYOUT_DIR.glob("*.yaml"))


def load_layout(layout: str | os.PathLike) -> dict:
    """A shipped layout by name (``default``) or a YAML file path; panel names are validated."""
    path = Path(layout)
    if not path.suffix:
        path = LAYOUT_DIR / f"{layout}.yaml"
    if not path.exists():
        raise FileNotFoundError(f"No layout {layout!r}; shipped: {available_layouts()}")
    spec = yaml.safe_load(path.read_text())
    unknown = [p for s in spec.get("sections", []) for p in s.get("panels", []) if p not in METRICS]
    if unknown:
        raise ValueError(f"Unknown panel(s) {unknown} in {path}; known: {sorted(METRICS)}")
    spec.setdefault("title", "DIA QC")
    spec.setdefault("columns", 2)
    return spec


# ----------------------------------------------------------------------------------------------
# per-panel figures (Plotly figure dicts built from per-run summaries)
# ----------------------------------------------------------------------------------------------


def _clean(obj):
    """JSON-safe: NaN/inf -> None, numpy scalars -> Python."""
    if isinstance(obj, dict):
        return {k: _clean(v) for k, v in obj.items()}
    if isinstance(obj, list | tuple):
        return [_clean(v) for v in obj]
    if isinstance(obj, np.generic):
        obj = obj.item()
    if isinstance(obj, float) and not math.isfinite(obj):
        return None
    return obj


class _Ctx:
    def __init__(
        self,
        prec: pd.DataFrame,
        groups: pd.Series | None,
        limpa=None,
        ann: pd.DataFrame | None = None,
        chrom: dict | None = None,
        history: pd.DataFrame | None = None,
    ):
        self.prec = prec
        self.limpa = limpa
        self.groups = groups
        self.ann = ann  # surface tier per protein group (dia_qc.surface_annotation)
        self.chrom = chrom  # raw_chromatography.load_chromatography output
        self.history = history  # qcdb latest_run_metrics (earlier runs), or None
        self._surface = None
        self._diag = None
        self.order = qc.run_order(prec)
        self.cond = prec.groupby("run")["condition"].first()
        self.labels = qc.run_labels(self.order, self.cond)
        conds = sorted(self.cond.astype(str).unique())
        pal = qc._palette(max(len(conds), len(self.order)))
        self.cond_color = dict(zip(conds, pal, strict=False))
        self.run_color = dict(zip(self.order, pal, strict=False))
        self._summary = None
        self._rt = None

    @property
    def rt(self) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        """(per precursor x run RT table, per run summary, charge-state spreads), cached."""
        if self._rt is None:
            rtc, summ = qc.rt_consistency(self.prec)
            self._rt = (rtc, summ, qc.charge_state_rt_spread(self.prec))
        return self._rt

    @property
    def summary(self) -> pd.DataFrame:
        if self._summary is None:
            self._summary = qc.run_summary(self.prec)
        return self._summary

    def color(self, run: str) -> str:
        return self.cond_color[str(self.cond[run])]

    @property
    def surface(self) -> dict:
        if self._surface is None:
            self._surface = qc.surface_by_tier(self.prec, self.ann)
        return self._surface

    @property
    def diagnosis(self) -> pd.DataFrame:
        """Per-run chromatography metrics + flags (raw traces joined with peak metrics)."""
        if self._diag is None:
            from deliverome_analysis import raw_chromatography as rc

            m = rc.tic_metrics(self.chrom["scans"])
            m = m.join(rc.scan_metrics(self.chrom["scans"]), how="left")
            tm = rc.trace_metrics(self.chrom.get("status"), self.chrom.get("devices"), m)
            if not tm.empty:
                m = m.join(tm, how="left")
            p = self.prec
            for col, name, scale in (("fwhm", "fwhm_s", 60.0), ("tailing_factor", "tailing", 1.0)):
                if col in p and p[col].notna().any():
                    m[name] = p.groupby("run")[col].median().reindex(m.index) * scale
            if "apex_rt" in p and len(self.order) > 1:
                summ = self.rt[1]
                if "median_shift_s" in summ:
                    m["rt_shift_min"] = (
                        summ.set_index("run")["median_shift_s"].reindex(m.index) / 60
                    )
            ids = self.summary["precursors"].reindex(m.index)
            m.insert(0, "precursors", ids)
            self._diag = rc.run_diagnostics(m.reindex([r for r in self.order if r in m.index]))
        return self._diag

    def run_index(self, run: str) -> int:
        return self.order.index(run)

    def protein_name(self, pid) -> str:
        """Display name: gene symbol(s) from the tier annotation, else the protein id."""
        if self.ann is not None and pid in self.ann.index and self.ann.at[pid, "symbol"]:
            return f"{self.ann.at[pid, 'symbol']} ({pid})"
        return str(pid)


def _layout(xlabel: str = "", ylabel: str = "", **kw) -> dict:
    lay = {
        "margin": {"l": 60, "r": 16, "t": 10, "b": 70},
        "height": 340,
        "font": {"family": "Manrope, Inter, sans-serif", "size": 12, "color": "#1f1718"},
        "paper_bgcolor": "rgba(0,0,0,0)",
        "plot_bgcolor": "rgba(0,0,0,0)",
        "colorway": ["#bc3c4c", "#3d6b60", "#c07830", "#5848a8", "#6e1428", "#7aab9f", "#f4aa28"],
        "hoverlabel": HOVER,
        # run labels like "01" must stay categories, not become numbers
        "xaxis": {
            "title": {"text": xlabel},
            "gridcolor": GRID,
            "automargin": True,
            "linecolor": "#e6dad4",
        },
        "yaxis": {"title": {"text": ylabel}, "gridcolor": GRID, "zerolinecolor": "#cfc4bf"},
        "legend": {"orientation": "h", "y": -0.25},
        "hovermode": "closest",
    }
    lay.update(kw)
    return lay


def _cat(lay: dict) -> dict:
    lay["xaxis"]["type"] = "category"
    return lay


def _fig_ids(ctx: _Ctx, m: Metric) -> dict:
    s = ctx.summary
    x = [ctx.labels[r] for r in s.index]
    trace = {
        "type": "bar",
        "x": x,
        "y": s[m.column].tolist(),
        "marker": {"color": [ctx.color(r) for r in s.index]},
        "hovertemplate": "%{x}<br>%{y:,}<extra></extra>",
    }
    return {"data": [trace], "layout": _cat(_layout("", m.label, showlegend=False))}


def _fig_dist(ctx: _Ctx, m: Metric) -> dict:
    p = ctx.prec
    vals = p[m.column].astype(float)
    if m.log:
        vals = vals.where(vals > 0)
    data = []
    for r in ctx.order:
        v = vals[(p["run"] == r).to_numpy()].dropna()
        q = np.percentile(v, [5, 25, 50, 75, 95]) if len(v) else [np.nan] * 5
        data.append(
            {
                "type": "box",
                "name": ctx.labels[r],
                "x": [ctx.labels[r]],
                "lowerfence": [q[0]],
                "q1": [q[1]],
                "median": [q[2]],
                "q3": [q[3]],
                "upperfence": [q[4]],
                "mean": [float(v.mean()) if len(v) else np.nan],
                "boxmean": True,
                "marker": {"color": ctx.color(r)},
                "line": {"color": ctx.color(r)},
                "fillcolor": ctx.color(r) + "99",
                "hoverinfo": "y+name",
            }
        )
    yaxis = {"title": {"text": m.label}, "gridcolor": GRID, "zerolinecolor": "#80706a"}
    if m.log:
        yaxis["type"] = "log"
    return {"data": data, "layout": _cat(_layout("", m.label, showlegend=False, yaxis=yaxis))}


def _fig_completeness(ctx: _Ctx, m: Metric) -> dict:
    compl = qc.precursor_completeness(ctx.prec)
    n = len(ctx.order)
    counts = compl["n_runs"].value_counts().reindex(range(1, n + 1), fill_value=0)
    trace = {
        "type": "bar",
        "x": [f"{k}/{n}" for k in counts.index],
        "y": (100 * counts / counts.sum()).round(2).tolist(),
        "customdata": counts.tolist(),
        "marker": {"color": "#3D6B60"},
        "hovertemplate": "detected in %{x} runs<br>%{y:.1f}% (%{customdata:,})<extra></extra>",
    }
    return {
        "data": [trace],
        "layout": _layout("Runs detected", "Precursors (%)", showlegend=False),
    }


def _fig_detection(ctx: _Ctx, m: Metric) -> dict:
    det = qc.conditional_detection(ctx.prec)
    trace = {
        "type": "scatter",
        "mode": "lines+markers",
        "x": det["log2_intensity"].round(3).tolist(),
        "y": det["p_detected"].round(4).tolist(),
        "customdata": det["n"].tolist(),
        "line": {"color": "#BC3C4C"},
        "hovertemplate": "log2 %{x:.1f}<br>P = %{y:.3f} (n=%{customdata:,})<extra></extra>",
    }
    lay = _layout("log2 intensity in one run", "P(detected in other runs)", showlegend=False)
    lay["yaxis"]["range"] = [0, 1.02]
    return {"data": [trace], "layout": lay}


def _fig_cv(ctx: _Ctx, m: Metric) -> dict | None:
    cv = qc.precursor_cv(ctx.prec, ctx.groups)
    if cv.empty:
        return None
    edges = np.linspace(0, 100, 51)
    data = []
    for i, (g, d) in enumerate(sorted(cv.groupby("group"))):
        pct = 100 * d["cv"].clip(upper=1.0)
        h, _ = np.histogram(pct, bins=edges)
        data.append(
            {
                "type": "scatter",
                "mode": "lines",
                "line": {"shape": "hvh", "width": 2, "color": qc._palette(i + 1)[i]},
                "x": ((edges[:-1] + edges[1:]) / 2).tolist(),
                "y": h.tolist(),
                "name": f"{g}: median {pct.median():.1f}%, {(pct < 20).mean():.0%} < 20%",
                "hovertemplate": "CV %{x:.0f}%<br>%{y:,} precursors<extra></extra>",
            }
        )
    return {"data": data, "layout": _layout("Precursor CV (%)", "Precursors")}


def _fig_gradient(ctx: _Ctx, m: Metric) -> dict:
    p = ctx.prec
    edges = np.linspace(p["apex_rt"].min(), p["apex_rt"].max(), 61)
    data = []
    for r in ctx.order:
        h, _ = np.histogram(p.loc[p["run"] == r, "apex_rt"], bins=edges)
        data.append(
            {
                "type": "scatter",
                "mode": "lines",
                "line": {"shape": "hvh", "width": 2, "color": ctx.run_color[r]},
                "x": ((edges[:-1] + edges[1:]) / 2).round(3).tolist(),
                "y": h.tolist(),
                "name": ctx.labels[r],
                "hovertemplate": "%{x:.1f} min<br>%{y:,}<extra>%{fullData.name}</extra>",
            }
        )
    return {"data": data, "layout": _layout("Apex retention time (min)", "Precursors per bin")}


def _fig_missed(ctx: _Ctx, m: Metric) -> dict:
    p = ctx.prec
    mc = p["missed_cleavages"].clip(upper=2).map({0: "0", 1: "1", 2: "2+"})
    frac = pd.crosstab(p["run"], mc, normalize="index").reindex(ctx.order).fillna(0)
    data = []
    for i, col in enumerate(c for c in ("0", "1", "2+") if c in frac.columns):
        data.append(
            {
                "type": "bar",
                "x": [ctx.labels[r] for r in ctx.order],
                "y": (100 * frac[col]).round(2).tolist(),
                "name": f"{col} missed",
                "marker": {"color": qc._palette(3)[i]},
                "hovertemplate": "%{x}<br>%{y:.1f}%<extra>%{fullData.name}</extra>",
            }
        )
    return {"data": data, "layout": _cat(_layout("", "Identified precursors (%)", barmode="stack"))}


# ----------------------------------------------------------------------------------------------
# completeness per run, charge states, surface tiers, chromatography
# ----------------------------------------------------------------------------------------------

# Deep green canonical, light green likely, light yellow low / uncertain, grey not surface.
TIER_COLORS = {
    "canonical": "#1F6B45",
    "likely": "#8FD0A0",
    "low / uncertain": "#F2D96B",
    "not surface": "#CFC4BF",
}


def _fig_compl_runs(ctx: _Ctx, m: Metric) -> dict | None:
    p = ctx.prec
    if len(ctx.order) < 2:
        return None
    levels = {"Precursors": "precursor", "Peptides": "peptide", "Protein groups": "protein"}
    rows = []
    for name, col in levels.items():
        if col not in p:
            continue
        total = p[col].nunique()
        per = p.groupby("run")[col].nunique()
        for r in ctx.order:
            rows.append((name, r, 100 * per.get(r, 0) / total))
    df = pd.DataFrame(rows, columns=["level", "run", "pct"])
    data = []
    for g in sorted(ctx.cond.astype(str).unique()):
        runs = [r for r in ctx.order if str(ctx.cond[r]) == g]
        d = df[df["run"].isin(runs)]
        data.append(
            {
                "type": "box",
                "name": g,
                "x": d["level"].tolist(),
                "y": d["pct"].round(3).tolist(),
                "text": [ctx.labels[r] for r in d["run"]],
                "customdata": [["run", ctx.run_index(r)] for r in d["run"]],
                "boxpoints": "all",
                "jitter": 0.4,
                "pointpos": 0,
                "marker": {
                    "color": ctx.cond_color[g],
                    "size": 9,
                    "line": {"width": 1, "color": "#fff"},
                },
                "line": {"color": ctx.cond_color[g]},
                "fillcolor": ctx.cond_color[g] + "33",
                "hovertemplate": "%{text}<br>%{x}: %{y:.1f}%<extra>%{fullData.name}</extra>",
            }
        )
    lay = _cat(_layout("", "Detected in this run (% of batch)", boxmode="group"))
    lay["showlegend"] = len(data) > 1
    return {"data": data, "layout": lay}


def _fig_charge(ctx: _Ctx, m: Metric) -> dict | None:
    z = pd.to_numeric(ctx.prec["precursor"].str.rsplit(".", n=1).str[-1], errors="coerce")
    if z.isna().all():
        return None
    lab = z.clip(upper=5).map(lambda v: "5+" if v >= 5 else f"{int(v)}+" if v == v else None)
    frac = pd.crosstab(ctx.prec["run"], lab, normalize="index").reindex(ctx.order).fillna(0)
    cols = sorted(frac.columns, key=lambda c: int(c.rstrip("+")))
    shades = ["#F4C070", "#3D6B60", "#BC3C4C", "#5848A8", "#1F1718"]
    data = [
        {
            "type": "bar",
            "x": [ctx.labels[r] for r in ctx.order],
            "y": (100 * frac[c]).round(2).tolist(),
            "name": c,
            "customdata": [["run", i] for i in range(len(ctx.order))],
            "marker": {"color": shades[min(int(c.rstrip("+")), 5) - 1]},
            "hovertemplate": "%{x}<br>%{y:.1f}% at " + c + "<extra></extra>",
        }
        for c in cols
    ]
    return {"data": data, "layout": _cat(_layout("", "Precursors (%)", barmode="stack"))}


def _tier_bars(ctx: _Ctx, table: pd.DataFrame, total=None, fmt=",.0f", x=None, offset=None) -> list:
    out = []
    xs = x or [ctx.labels[r] for r in ctx.order]
    for t in qc.SURFACE_CLASSES:
        vals = table[t].reindex(ctx.order)
        pct = (
            (100 * vals / total.reindex(ctx.order)).round(1).tolist() if total is not None else None
        )
        out.append(
            {
                "type": "bar",
                "x": xs,
                "y": vals.round(3).tolist(),
                "name": t,
                "legendgroup": t,
                "showlegend": offset in (None, 0),
                "offsetgroup": offset,
                "marker": {"color": TIER_COLORS[t]},
                "customdata": [["tier", i, t] for i in range(len(ctx.order))],
                "text": pct,
                "textposition": "none",
                "hovertemplate": "%{x}<br>"
                + t
                + ": %{y:"
                + fmt
                + "}"
                + (" (%{text}% of all)" if total is not None else "")
                + "<extra></extra>",
            }
        )
    return out


def _fig_surf_ids(ctx: _Ctx, m: Metric) -> dict:
    tab = ctx.surface[m.column]
    data = _tier_bars(ctx, tab, total=tab.sum(axis=1))
    return {"data": data, "layout": _cat(_layout("", m.label, barmode="stack"))}


def _fig_surf_signal(ctx: _Ctx, m: Metric) -> dict:
    s = ctx.surface
    levels = [(k, lab) for k, lab in (("signal_ms2", "MS2"), ("signal_ms1", "MS1")) if k in s]
    data = []
    for li, (key, lab) in enumerate(levels):
        x = [ctx.labels[r] for r in ctx.order]
        if len(levels) > 1:
            x = [x, [lab] * len(x)]
        data += _tier_bars(ctx, s[key], fmt=".2f", x=x, offset=li if len(levels) > 1 else None)
    for t in data:
        t["hovertemplate"] = t["hovertemplate"].replace("}<extra>", "}% of run signal<extra>")
    lay = _layout("", "% of run signal", barmode="stack")
    if len(levels) == 1:
        lay = _cat(lay)
        lay["annotations"] = [
            {
                "text": "MS2 signal (EG.TotalQuantity); MS1 needs FG.MS1Quantity in the export",
                "xref": "paper",
                "yref": "paper",
                "x": 0,
                "y": 1.06,
                "showarrow": False,
                "font": {"size": 11, "color": "#6F5D5A"},
                "xanchor": "left",
            }
        ]
    return {"data": data, "layout": lay}


def _fig_surf_missed(ctx: _Ctx, m: Metric) -> dict | None:
    s = ctx.surface
    if "missed" not in s:
        return None
    x = [ctx.labels[r] for r in ctx.order]
    data = [
        {
            "type": "bar",
            "x": x,
            "y": s["missed"][t].reindex(ctx.order).round(2).tolist(),
            "name": t,
            "marker": {"color": TIER_COLORS[t]},
            "customdata": [
                ["tier", i, t] if t != "not surface" else ["run", i] for i in range(len(x))
            ],
            "hovertemplate": "%{x}<br>" + t + ": %{y:.1f}%<extra></extra>",
        }
        for t in qc.TIER_CLASSES
    ]
    data.append(
        {
            "type": "scatter",
            "mode": "lines+markers",
            "x": x,
            "y": s["missed_all"].round(2).tolist(),
            "name": "all precursors",
            "line": {"color": "#1F1718", "dash": "dot"},
            "customdata": [["run", i] for i in range(len(x))],
            "hovertemplate": "%{x}<br>all: %{y:.1f}%<extra></extra>",
        }
    )
    return {
        "data": data,
        "layout": _cat(_layout("", "Precursors with a missed cleavage (%)", barmode="group")),
    }


def _fig_surf_points(ctx: _Ctx, m: Metric) -> dict:
    p = qc._with_tier(ctx.prec, ctx.ann)
    data = []
    for t in qc.TIER_CLASSES:
        xs, q = [], {k: [] for k in ("lowerfence", "q1", "median", "q3", "upperfence")}
        for r in ctx.order:
            v = p.loc[(p["run"] == r) & (p["tier"] == t), m.column].dropna()
            if v.empty:
                continue
            xs.append(ctx.labels[r])
            for k, val in zip(q, np.percentile(v, [5, 25, 50, 75, 95]), strict=True):
                q[k].append(float(val))
        if xs:
            data.append(
                {
                    "type": "box",
                    "name": t,
                    "x": xs,
                    **q,
                    "marker": {"color": TIER_COLORS[t]},
                    "line": {"color": TIER_COLORS[t]},
                    "fillcolor": TIER_COLORS[t] + "88",
                }
            )
    return {"data": data, "layout": _cat(_layout("", m.label, boxmode="group"))}


def _fig_tic(ctx: _Ctx, m: Metric) -> dict | None:
    from deliverome_analysis import raw_chromatography as rc

    sc = ctx.chrom["scans"]
    sc = sc[sc["ms_order"] == int(m.column)]
    if sc.empty:
        return None
    data = []
    for r in ctx.order:
        d = sc[sc["run"] == r]
        if d.empty:
            continue
        b = rc.binned_trace(d, "tic")
        data.append(
            {
                "type": "scattergl",
                "mode": "lines",
                "x": b["rt"].round(3).tolist(),
                "y": b["tic"].round(0).tolist(),
                "name": ctx.labels[r],
                "line": {"color": ctx.run_color[r], "width": 1.3},
                "customdata": [["run", ctx.run_index(r)]] * len(b),
                "hovertemplate": "%{x:.2f} min<br>TIC %{y:.3s}<extra>%{fullData.name}</extra>",
            }
        )
    lay = _layout("Retention time (min)", m.label, height=320)
    lay["yaxis"]["exponentformat"] = "SI"
    return {"data": data, "layout": lay}


def _fig_inject(ctx: _Ctx, m: Metric) -> dict | None:
    from deliverome_analysis import raw_chromatography as rc

    sc = ctx.chrom["scans"]
    sc = sc[(sc["ms_order"] == 1) & sc["inject_ms"].notna()]
    if sc.empty:
        return None
    data = []
    for r in ctx.order:
        d = sc[sc["run"] == r]
        if d.empty:
            continue
        b = rc.binned_trace(d, "inject_ms", bins=600)
        data.append(
            {
                "type": "scattergl",
                "mode": "lines",
                "x": b["rt"].round(3).tolist(),
                "y": b["inject_ms"].round(3).tolist(),
                "name": ctx.labels[r],
                "line": {"color": ctx.run_color[r], "width": 1.3},
                "customdata": [["run", ctx.run_index(r)]] * len(b),
                "hovertemplate": "%{x:.2f} min<br>%{y:.2f} ms<extra>%{fullData.name}</extra>",
            }
        )
    return {
        "data": data,
        "layout": _layout("Retention time (min)", "Injection time (ms)", height=320),
    }


def _fig_traces(ctx: _Ctx, m: Metric) -> dict | None:
    """One trace per (signal, run); a dropdown picks the signal."""
    from deliverome_analysis import raw_chromatography as rc

    series = []
    st = ctx.chrom.get("status")
    if st is not None and not st.empty:
        for key, d in st.groupby("key"):
            if rc.DIAGNOSTIC_STATUS.search(key) and d["value"].nunique() > 1:
                series.append((f"status: {key}", d.rename(columns={"value": "v"})))
    dv = ctx.chrom.get("devices")
    if dv is not None and not dv.empty:
        for name, d in dv.groupby("name"):
            if d["value"].nunique() > 1:
                series.append((f"LC/device: {name}", d.rename(columns={"value": "v"})))
    if not series:
        return None
    pri = ["Pump_Pressure", "Spray Current (", "Spray Voltage"]
    series.sort(
        key=lambda kv: (next((i for i, k in enumerate(pri) if k in kv[0]), len(pri)), kv[0])
    )
    data, owner = [], []
    for si, (_name, d) in enumerate(series):
        for r in ctx.order:
            dr = d[d["run"] == r]
            if dr.empty:
                continue
            b = rc.binned_trace(dr, "v", bins=600)
            data.append(
                {
                    "type": "scattergl",
                    "mode": "lines",
                    "x": b["rt"].round(3).tolist(),
                    "y": b["v"].round(4).tolist(),
                    "name": ctx.labels[r],
                    "visible": si == 0,
                    "line": {"color": ctx.run_color[r], "width": 1.3},
                    "customdata": [["run", ctx.run_index(r)]] * len(b),
                    "hovertemplate": "%{x:.2f} min<br>%{y:.4g}<extra>%{fullData.name}</extra>",
                }
            )
            owner.append(si)
    buttons = [
        {
            "label": name,
            "method": "update",
            "args": [
                {"visible": [o == si for o in owner]},
                {"yaxis.title.text": name.split(": ", 1)[1]},
            ],
        }
        for si, (name, _) in enumerate(series)
    ]
    lay = _layout("Retention time (min)", series[0][0].split(": ", 1)[1], height=340)
    lay["updatemenus"] = [
        {"buttons": buttons, "x": 0, "y": 1.18, "xanchor": "left", "showactive": True}
    ]
    lay["margin"]["t"] = 40
    return {"data": data, "layout": lay}


HISTORY_METRICS = [
    ("precursors", "Precursors"),
    ("protein_groups", "Protein groups"),
    ("protein_groups_canonical", "Canonical surface protein groups"),
    ("tic_area", "TIC area"),
    ("roughness", "TIC roughness"),
    ("spray_noise_pct", "Spray noise (%)"),
    ("pump_pressure_max", "Max pump pressure"),
    ("fwhm_s", "Median FWHM (s)"),
    ("rt_shift_min", "RT shift (min)"),
    ("missed_cleavage_pct", "Missed cleavage (%)"),
    ("cycle_time_s", "Cycle time (s)"),
    ("inject_ms", "MS1 fill time (ms)"),
]


def _history_frame(ctx: _Ctx) -> pd.DataFrame:
    """Earlier runs from the index plus this page's runs (run, metric, value, when, current)."""
    found = (ctx.chrom or {}).get("found") or {}
    cur = run_metrics(ctx)
    cur["when"] = cur["run"].map(lambda r: ((found.get(r) or {}).get("info") or {}).get("acquired"))
    cur["current"] = True
    frames = [cur]
    h = ctx.history
    if h is not None and not h.empty:
        from deliverome_analysis.qcdb import run_uid

        # exclude this page's own acquisitions (by instrument serial + acquisition time when
        # known, so two instruments' "Sample_01" are not confused), then add the rest
        mine = {run_uid({"run": r, **((found.get(r) or {}).get("info") or {})}) for r in ctx.order}
        key = h["run_uid"] if "run_uid" in h else "name:" + h["run"]
        h = h[~key.isin(mine)].assign(current=False)
        frames.append(
            h.rename(columns={"acquired": "when"})[["run", "metric", "value", "when", "current"]]
        )
    out = pd.concat(frames, ignore_index=True).dropna(subset=["when"])
    out["when"] = pd.to_datetime(out["when"], errors="coerce", utc=True, format="ISO8601")
    return out.dropna(subset=["when"])


def _fig_history(ctx: _Ctx, m: Metric) -> dict | None:
    df = _history_frame(ctx)
    metrics = [(k, lab) for k, lab in HISTORY_METRICS if (df["metric"] == k).any()]
    if not metrics:
        return None
    data, owner = [], []
    for mi, (key, _lab) in enumerate(metrics):
        d = df[df["metric"] == key].sort_values("when")
        for current in (False, True):
            dd = d[d["current"] == current]
            if dd.empty:
                continue
            data.append(
                {
                    "type": "scatter",
                    "mode": "markers",
                    "x": dd["when"].dt.strftime("%Y-%m-%d %H:%M").tolist(),
                    "y": dd["value"].round(4).tolist(),
                    "text": dd["run"].tolist(),
                    "name": "this page" if current else "earlier runs",
                    "visible": mi == 0,
                    "showlegend": mi == 0,
                    "marker": {
                        "size": 11 if current else 7,
                        "color": "#BC3C4C" if current else "#80706A",
                        "opacity": 1 if current else 0.6,
                        "line": {"width": 1, "color": "#fff"},
                    },
                    "customdata": [
                        ["run", ctx.run_index(r)] if r in ctx.order else None for r in dd["run"]
                    ],
                    "hovertemplate": "%{text}<br>%{x}<br>%{y:.4g}<extra></extra>",
                }
            )
            owner.append(mi)
    buttons = [
        {
            "label": lab,
            "method": "update",
            "args": [{"visible": [o == mi for o in owner]}, {"yaxis.title.text": lab}],
        }
        for mi, (_k, lab) in enumerate(metrics)
    ]
    lay = _layout("Acquired", metrics[0][1], height=340)
    lay["xaxis"]["type"] = "date"
    lay["updatemenus"] = [{"buttons": buttons, "x": 0, "y": 1.18, "xanchor": "left"}]
    lay["margin"]["t"] = 40
    return {"data": data, "layout": lay}


_EMPTY_REASON = {
    "contam": "No contaminant, watchlist or iRT precursors: was this searched without the "
    "published search database (contaminants and iRT)?",
    "limpa_pca": "PCA needs at least 4 runs.",
    "limpa_agree": "Needs PG.Quantity in the export.",
    "cv": "Needs at least 2 runs in a replicate group.",
    "chargert": "No modified peptide was seen at two charge states in one run.",
    "rtcons": "Needs at least 2 runs with enough shared precursors.",
    "history": "No acquisition dates: pass chromatography_dir (raw-file headers) to place runs in time.",
}
_BOX_NOTE = "Box: quartiles; whiskers: 5th-95th percentile; dashed line: mean."


def _fig_rtcons(ctx: _Ctx, m: Metric) -> dict | None:
    rtc, summ, _ = ctx.rt
    if rtc.empty:
        return None
    s_by = summ.set_index("run")
    data = []
    for r in [r for r in ctx.order if r in set(rtc["run"])]:
        d = rtc[rtc["run"] == r]
        b = pd.qcut(d["ref_rt"], max(3, min(40, len(d) // 20)), duplicates="drop")
        g = d.groupby(b, observed=True)
        x = g["ref_rt"].median().round(3).tolist()
        lo = (60 * g["delta_min"].quantile(0.05)).round(2).tolist()
        hi = (60 * g["delta_min"].quantile(0.95)).round(2).tolist()
        c = ctx.run_color[r]
        data.append(
            {
                "type": "scatter",
                "mode": "lines",
                "x": x + x[::-1],
                "y": hi + lo[::-1],
                "fill": "toself",
                "fillcolor": c + "26",
                "line": {"width": 0},
                "hoverinfo": "skip",
                "showlegend": False,
            }
        )
        sm = s_by.loc[r]
        data.append(
            {
                "type": "scatter",
                "mode": "lines",
                "x": x,
                "y": (60 * g["delta_min"].median()).round(2).tolist(),
                "line": {"color": c, "width": 2.5},
                "name": f"{ctx.labels[r]}: shift {sm['median_shift_s']:+.1f} s, "
                f"{int(sm['flagged']):,} flagged ({sm['flagged_pct']:.1f}%)",
                "hovertemplate": "%{x:.1f} min: %{y:+.1f} s<extra>%{fullData.name}</extra>",
            }
        )
    lim = float(60 * np.nanpercentile(np.abs(rtc["delta_min"]), 99))
    lay = _layout("Retention time, median of runs (min)", "This run minus median (s)")
    lay["yaxis"]["range"] = [-max(lim, 3), max(lim, 3)]
    return {"data": data, "layout": lay}


def _fig_chargert(ctx: _Ctx, m: Metric) -> dict | None:
    _, _, cs = ctx.rt
    if cs.empty:
        return None
    sec = (60 * cs["spread_min"]).clip(upper=60)
    edges = np.arange(0, 61, 1.0)
    h, _ = np.histogram(sec, bins=edges)
    trace = {
        "type": "bar",
        "x": (edges[:-1] + 0.5).tolist(),
        "y": h.tolist(),
        "marker": {"color": ["#3D6B60" if e < 6 else "#BC3C4C" for e in edges[:-1]]},
        "hovertemplate": "%{x:.0f} s: %{y:,} peptides<extra></extra>",
    }
    lay = _layout(
        "Apex RT spread across charge states (s; 60 = 60 or more)",
        "Peptide-run pairs",
        showlegend=False,
    )
    lay["yaxis"]["type"] = "log"
    lay["annotations"] = [
        {
            "x": 0.98,
            "y": 0.95,
            "xref": "paper",
            "yref": "paper",
            "xanchor": "right",
            "showarrow": False,
            "font": {"size": 11},
            "text": f"{int(cs['flag'].sum()):,} of {len(cs):,} ({cs['flag'].mean():.1%}) above 6 s",
        }
    ]
    return {"data": [trace], "layout": lay}


def _fig_limpa_dpc(ctx: _Ctx, m: Metric) -> dict:
    from deliverome_analysis.limpa import detection_probability

    pts = ctx.limpa.dpc_points.dropna(subset=["mu"])
    frac = pts["n_detected"] / pts["nsamples"]
    edges = np.linspace(pts["mu"].min(), pts["mu"].max(), 31)
    idx = np.clip(np.digitize(pts["mu"], edges) - 1, 0, 29)
    b = (
        pd.DataFrame({"bin": idx, "mu": pts["mu"], "frac": frac})
        .groupby("bin")
        .agg(mu=("mu", "mean"), frac=("frac", "mean"), n=("frac", "size"))
    )
    b = b[b["n"] >= 10]
    b0, b1 = ctx.limpa.dpc
    x = np.linspace(pts["mu"].min(), pts["mu"].max(), 200)
    data = [
        {
            "type": "scatter",
            "mode": "markers",
            "x": b["mu"].round(3).tolist(),
            "y": b["frac"].round(4).tolist(),
            "customdata": b["n"].tolist(),
            "marker": {"color": "#BC3C4C", "size": 9},
            "name": "observed (bin mean)",
            "hovertemplate": "log2 %{x:.1f}: %{y:.3f} (n=%{customdata:,})<extra></extra>",
        },
        {
            "type": "scatter",
            "mode": "lines",
            "x": x.round(3).tolist(),
            "y": detection_probability(x, b0, b1).round(4).tolist(),
            "line": {"color": "#1F1718", "width": 2.5},
            "name": f"limpa DPC (b0={b0:.2f}, b1={b1:.2f})",
        },
    ]
    lay = _layout("limpa estimated mean log2 intensity", "Fraction of runs detected")
    lay["yaxis"]["range"] = [0, 1.05]
    return {"data": data, "layout": lay}


def _fig_limpa_se(ctx: _Ctx, m: Metric) -> dict:
    r = ctx.limpa
    mean, se = r.protein_log2.mean(axis=1), r.protein_se.mean(axis=1)
    nprec = r.proteins["NPrec"].reindex(mean.index)
    bucket = pd.cut(nprec, [0, 1, 3, 9, np.inf], labels=["1", "2-3", "4-9", "10+"])
    data = []
    for i, c in enumerate(bucket.cat.categories):
        sel = (bucket == c).to_numpy()
        if not sel.any():
            continue
        data.append(
            {
                "type": "scattergl",
                "mode": "markers",
                "x": mean[sel].round(2).tolist(),
                "y": se[sel].round(3).tolist(),
                "text": list(mean.index[sel]),
                "name": f"{c} precursors (n={int(sel.sum()):,})",
                "marker": {"size": 4, "opacity": 0.5, "color": qc._palette(4)[i]},
                "hovertemplate": "%{text}<br>log2 %{x:.2f}, SE %{y:.3f}<extra></extra>",
            }
        )
    return {
        "data": data,
        "layout": _layout(
            "Protein log2 (limpa, mean over runs)", "Standard error (mean over runs)"
        ),
    }


def _fig_limpa_corr(ctx: _Ctx, m: Metric) -> dict:
    mat = ctx.limpa.protein_log2.reindex(columns=ctx.order)
    seen = ctx.prec.groupby("protein")["run"].nunique()
    mat = mat.loc[mat.index.intersection(seen.index[seen == len(ctx.order)])]
    corr = mat.corr()
    names = [ctx.labels[r] for r in ctx.order]
    trace = {
        "type": "heatmap",
        "z": corr.round(4).values.tolist(),
        "x": names,
        "y": names,
        "colorscale": [[0, "#F3ECE5"], [1, "#3D6B60"]],
        "zmin": float(min(0.9, np.nanmin(corr.values))),
        "zmax": 1,
        "text": corr.round(3).values.tolist(),
        "texttemplate": "%{text}",
        "hovertemplate": "%{y} vs %{x}: r = %{z:.4f}<extra></extra>",
        "colorbar": {"title": {"text": "r"}},
    }
    lay = _cat(_layout("", "", showlegend=False))
    lay["yaxis"]["type"] = "category"
    lay["annotations"] = [
        {
            "x": 0,
            "y": 1.08,
            "xref": "paper",
            "yref": "paper",
            "showarrow": False,
            "xanchor": "left",
            "font": {"size": 11},
            "text": f"{len(mat):,} proteins detected in all runs",
        }
    ]
    return {"data": [trace], "layout": lay}


def _fig_limpa_pca(ctx: _Ctx, m: Metric) -> dict | None:
    mat = ctx.limpa.protein_log2.reindex(columns=ctx.order).dropna()
    if mat.shape[1] < 4:
        return None
    x = (mat.T - mat.T.mean()).to_numpy()
    u, sv, _ = np.linalg.svd(x, full_matrices=False)
    var = sv**2 / (sv**2).sum()
    pcs = u[:, :2] * sv[:2]
    data = [
        {
            "type": "scatter",
            "mode": "markers+text",
            "x": [float(pcs[i, 0])],
            "y": [float(pcs[i, 1])],
            "text": [ctx.labels[r]],
            "textposition": "top center",
            "marker": {"size": 12, "color": ctx.color(r)},
            "name": ctx.labels[r],
        }
        for i, r in enumerate(ctx.order)
    ]
    return {
        "data": data,
        "layout": _layout(f"PC1 ({var[0]:.0%})", f"PC2 ({var[1]:.0%})", showlegend=False),
    }


def _fig_limpa_agree(ctx: _Ctx, m: Metric) -> dict | None:
    if "pg_quantity" not in ctx.prec:
        return None
    a = qc.protein_agreement(ctx.prec, ctx.limpa).dropna(
        subset=["limpa_centred", "spectronaut_centred"]
    )
    if a.empty:
        return None
    lim = float(np.nanpercentile(np.abs(a[["limpa_centred", "spectronaut_centred"]]), 99.5))
    edges = np.linspace(-lim, lim, 81)
    h, _, _ = np.histogram2d(a["spectronaut_centred"], a["limpa_centred"], bins=[edges, edges])
    c = (edges[:-1] + edges[1:]) / 2
    r = float(np.corrcoef(a["spectronaut_centred"], a["limpa_centred"])[0, 1])
    trace = {
        "type": "heatmap",
        "x": c.round(3).tolist(),
        "y": c.round(3).tolist(),
        "z": np.where(h.T > 0, np.log10(h.T), np.nan).round(3).tolist(),
        "colorscale": [[0, "#F3ECE5"], [1, "#BC3C4C"]],
        "showscale": False,
        "hoverinfo": "skip",
        "name": "all protein-run values",
    }
    # Name the proteins where the two disagree most: those are what is worth a look.
    a = a.assign(gap=(a["limpa_centred"] - a["spectronaut_centred"]).abs())
    top = a.nlargest(300, "gap")
    names = [ctx.protein_name(p) for p in top["protein"]]
    pts = {
        "type": "scattergl",
        "mode": "markers",
        "x": top["spectronaut_centred"].round(3).tolist(),
        "y": top["limpa_centred"].round(3).tolist(),
        "text": names,
        "customdata": [["pid", str(p)] for p in top["protein"]],
        "meta": [ctx.labels.get(r, r) for r in top["run"]],
        "marker": {
            "size": 6,
            "color": "#3D6B60",
            "opacity": 0.85,
            "line": {"width": 0.5, "color": "#fff"},
        },
        "name": "300 most discordant",
        "hovertemplate": "<b>%{text}</b><br>run %{meta}<br>Spectronaut %{x:+.2f}, limpa %{y:+.2f}"
        "<extra>click for details</extra>",
    }
    lay = _layout("Spectronaut, centred log2", "limpa, centred log2", showlegend=False)
    lay["annotations"] = [
        {
            "x": 0.02,
            "y": 0.98,
            "xref": "paper",
            "yref": "paper",
            "xanchor": "left",
            "showarrow": False,
            "font": {"size": 11},
            "text": f"Pearson r = {r:.3f} on {len(a):,} protein-run values",
        }
    ]
    lay["shapes"] = [
        {
            "type": "line",
            "x0": -lim,
            "y0": -lim,
            "x1": lim,
            "y1": lim,
            "line": {"color": "#6F5D5A", "dash": "dash", "width": 1},
        }
    ]
    return {"data": [trace, pts], "layout": lay}


def _has_classes(ctx: _Ctx) -> bool:
    p = ctx.prec
    return "protein_class" in p and p["protein_class"].ne("human").any()


def _fig_contam(ctx: _Ctx, m: Metric) -> dict | None:
    if not _has_classes(ctx):
        return None
    cs = qc.class_signal(ctx.prec)
    x = [ctx.labels[r] for r in cs.index]
    colors = {
        "contaminant": "#BC3C4C",
        "mixed": "#F4AA28",
        "watchlist": "#8878C8",
        "irt": "#3D6B60",
    }
    data = [
        {
            "type": "bar",
            "x": x,
            "y": cs[f"{c}_signal_pct"].round(4).tolist(),
            "name": c,
            "customdata": cs[f"{c}_precursors"].tolist(),
            "marker": {"color": col},
            "hovertemplate": "%{x}: %{y:.3f}% (%{customdata:,} precursors)"
            "<extra>%{fullData.name}</extra>",
        }
        for c, col in colors.items()
    ]
    return {
        "data": data,
        "layout": _cat(_layout("", "Share of precursor signal (%)", barmode="stack")),
    }


_BUILDERS = {
    "ids": _fig_ids,
    "dist": _fig_dist,
    "completeness": _fig_completeness,
    "detection": _fig_detection,
    "cv": _fig_cv,
    "gradient": _fig_gradient,
    "missed": _fig_missed,
    "rtcons": _fig_rtcons,
    "contam": _fig_contam,
    "chargert": _fig_chargert,
    "limpa_dpc": _fig_limpa_dpc,
    "limpa_se": _fig_limpa_se,
    "limpa_corr": _fig_limpa_corr,
    "limpa_pca": _fig_limpa_pca,
    "limpa_agree": _fig_limpa_agree,
    "compl_runs": _fig_compl_runs,
    "charge": _fig_charge,
    "surf_ids": _fig_surf_ids,
    "surf_signal": _fig_surf_signal,
    "surf_missed": _fig_surf_missed,
    "surf_points": _fig_surf_points,
    "tic": _fig_tic,
    "inject": _fig_inject,
    "traces": _fig_traces,
    "history": _fig_history,
}
_REQUIRED_COLUMN = {
    "gradient": "apex_rt",
    "missed": "missed_cleavages",
    "rtcons": "apex_rt",
    "chargert": "apex_rt",
    "rtflag": "apex_rt",
    "contam": "protein_class",
    "contamtable": "protein_class",
}


_CHROM_KINDS = ("tic", "inject", "traces", "diagnosis")


def panel_status(
    prec: pd.DataFrame, name: str, limpa=None, *, ann=None, chrom=None
) -> tuple[str, str]:
    """('shown' | 'missing' | 'planned', reason) for one panel on this precursor table."""
    m = METRICS[name]
    if m.kind.startswith("surf"):
        if ann is None:
            return "missing", (
                "Needs the accessible-surfaceome tiers; the catalog could not be loaded "
                "(S3 credentials?) or surface=False was passed."
            )
        col = (
            m.column
            if m.kind == "surf_points"
            else "missed_cleavages"
            if m.kind == "surf_missed"
            else None
        )
        if col and (col not in prec or prec[col].notna().sum() == 0):
            return "missing", f"Not in this export; needs {m.needs}."
        return "shown", ""
    if m.kind in _CHROM_KINDS:
        if not chrom or chrom.get("scans") is None:
            return "missing", f"Needs {m.needs}, passed as chromatography_dir."
        return "shown", ""
    if m.kind.startswith("limpa"):
        if limpa is None:
            return "missing", "Needs a limpa run: pass limpa_result (see scratch/becca/astral-xic)."
        if m.kind == "limpa_dpc" and limpa.dpc_points is None:
            return "missing", "This limpa output has no DPC points."
        return "shown", ""
    if m.kind == "planned":
        return "planned", f"Needs {m.needs}."
    col = m.column if m.kind == "dist" else _REQUIRED_COLUMN.get(m.kind)
    if col and (col not in prec or prec[col].notna().sum() == 0):
        return "missing", f"Not in this export; needs {m.needs or col}."
    return "shown", ""


# ----------------------------------------------------------------------------------------------
# protein explorer (client-side: data embedded once, all interaction in the page)
# ----------------------------------------------------------------------------------------------


# Precursor x run values embedded for the click-through; above this the page would pass
# ~50 MB, so the explorer keeps protein-level values and says so.
MAX_EMBEDDED_PRECURSOR_VALUES = 3_000_000


def _explorer_data(ctx: _Ctx, limpa_result=None, names: dict | None = None) -> dict:
    """Protein x run values, display names and per-protein precursor traces, compactly.

    Protein values: limpa's DPC-Quant when a LimpaResult is given, else Spectronaut's
    PG.Quantity (log2), else the median precursor log2 intensity per protein and run.
    """
    p, runs = ctx.prec, ctx.order
    if limpa_result is not None:
        prot = limpa_result.protein_log2.reindex(columns=runs)
        source = "limpa DPC-Quant (log2)"
    elif "pg_quantity" in p and p["pg_quantity"].notna().any():
        pg = p.dropna(subset=["pg_quantity"]).query("pg_quantity > 0")
        prot = np.log2(
            pg.pivot_table(index="protein", columns="run", values="pg_quantity", aggfunc="first")
        ).reindex(columns=runs)
        source = "Spectronaut PG.Quantity (log2)"
    else:
        prot = p.pivot_table(
            index="protein", columns="run", values="log2_intensity", aggfunc="median"
        ).reindex(columns=runs)
        source = "median precursor log2 intensity"
    prot = prot.dropna(how="all")
    pidx = {pid: i for i, pid in enumerate(prot.index)}
    genes = p.groupby("protein")["gene"].first() if "gene" in p else pd.Series(dtype=str)
    ann = ctx.ann
    label, tier, sym, acc = [], [], [], []
    for pid in prot.index:
        a = ann.loc[pid] if ann is not None and pid in ann.index else None
        nm = (
            (names or {}).get(pid)
            or (a["symbol"] if a is not None and a["symbol"] else None)
            or (genes.get(pid) if len(genes) else None)
        )
        label.append(f"{nm} ({pid})" if nm and nm != pid else str(pid))
        tier.append(a["tier"] if a is not None else "")
        sym.append(a["symbol"] if a is not None else (nm or ""))
        acc.append(a["accessions"] if a is not None else str(pid))
    wide = (
        p[p["protein"].isin(pidx)]
        .pivot_table(
            index=["protein", "precursor"], columns="run", values="log2_intensity", aggfunc="first"
        )
        .reindex(columns=runs)
    )
    nprec = wide.groupby(level=0).size().reindex(prot.index).fillna(0).astype(int)
    detail = wide.size <= MAX_EMBEDDED_PRECURSOR_VALUES
    if not detail:
        wide = wide.iloc[:0]
    r2 = lambda a: np.round(a.astype(float), 2)  # noqa: E731 - 0.01 log2 is ample for display
    return _clean(
        {
            "runs": [ctx.labels[r] for r in runs],
            "group": [str(ctx.cond[r]) for r in runs],
            "source": source,
            "detail": bool(detail),
            "pid": [str(x) for x in prot.index],
            "label": label,
            "tier": tier,
            "sym": sym,
            "acc": acc,
            "nprec": nprec.tolist(),
            "val": r2(prot.to_numpy()).tolist(),
            "prec_p": [pidx[i] for i in wide.index.get_level_values(0)],
            "prec_id": [str(x) for x in wide.index.get_level_values(1)],
            "prec_val": r2(wide.to_numpy()).tolist(),
            "prec_flag": _prec_flags(ctx, wide),
        }
    )


def _prec_flags(ctx: _Ctx, wide: pd.DataFrame) -> list[str]:
    """Per embedded precursor: the run labels where its RT was flagged ('' if none)."""
    if "apex_rt" not in ctx.prec or ctx.prec["run"].nunique() < 2:
        return [""] * len(wide)
    rtc = ctx.rt[0]
    if rtc.empty:
        return [""] * len(wide)
    f = (
        rtc[rtc["flag"]]
        .groupby("precursor")["run"]
        .agg(lambda r: ", ".join(ctx.labels[x] for x in r))
    )
    return [f.get(p, "") for p in wide.index.get_level_values(1)]


def _inspect_data(ctx: _Ctx) -> dict:
    """What the click-through inspector shows for a run, and the proteins behind a tier bar."""
    s = ctx.summary
    surf = ctx.surface if ctx.ann is not None else {}
    diag = (
        _diagnosis_frame(ctx)
        if ctx.chrom and ctx.chrom.get("scans") is not None
        else pd.DataFrame()
    )
    runs = []
    for r in ctx.order:
        rows = [
            ["Sample group", str(ctx.cond[r])],
            *[
                [lab, f"{int(s.loc[r, c]):,}"]
                for c, lab in (
                    ("precursors", "Precursors"),
                    ("peptides", "Peptides"),
                    ("protein_groups", "Protein groups"),
                )
                if c in s
            ],
        ]
        if surf:
            for t in qc.SURFACE_CLASSES:
                rows.append([f"{t} protein groups", f"{int(surf['protein_groups'].loc[r, t]):,}"])
            rows.append(
                [
                    "Surface share of MS2 signal",
                    f"{surf['signal_ms2'].loc[r, qc.SURFACE_CLASSES].sum():.1f}%",
                ]
            )
            if "signal_ms1" in surf:
                rows.append(
                    [
                        "Surface share of MS1 signal",
                        f"{surf['signal_ms1'].loc[r, qc.SURFACE_CLASSES].sum():.1f}%",
                    ]
                )
        if "missed_cleavages" in ctx.prec:
            rows.append(
                [
                    "Missed cleavage (any)",
                    f"{100 * (ctx.prec.loc[ctx.prec['run'] == r, 'missed_cleavages'] > 0).mean():.1f}%",
                ]
            )
        flags = hints = ""
        if not diag.empty and r in diag.index:
            d = diag.loc[r]
            for c, lab, fmt in _DIAG_COLS[1:]:
                if c in d and pd.notna(d[c]):
                    rows.append([lab, fmt.format(d[c])])
            flags, hints = d["flags"], d["hints"]
        runs.append(
            {"label": ctx.labels[r], "run": r, "rows": rows, "flags": flags, "hints": hints}
        )
    top = {}
    if surf:
        t = qc.top_surface_proteins(ctx.prec, ctx.ann)
        for (run, tier), g in t.groupby(["run", "tier"], observed=True):
            top.setdefault(str(ctx.run_index(run)), {})[str(tier)] = [
                [str(pid), round(float(sh), 3)]
                for pid, sh in zip(g["protein"], g["share_pct"], strict=True)
            ]
    return _clean({"runs": runs, "top": top, "tierColors": TIER_COLORS})


_EXPLORER_HTML = """
<div class="card wide"><h3>Protein explorer {info}</h3>
<div class="controls"><label>Run <select id="ex-run"></select></label>
<label>Find <input id="ex-q" placeholder="gene or accession"></label>
<label class="chk"><input type="checkbox" id="ex-surf"> Surface proteins only</label>
<span class="note" id="ex-src"></span></div>
<div class="ex-grid"><div id="ex-scatter"></div>
<div><div class="note" id="ex-cap"></div><div class="scroll ex-table"><table id="ex-tbl"></table></div></div></div>
<div id="ex-detail" class="ex-detail"><p class="note">Click a protein to see its precursors.</p></div>
{note}</div>"""

_EXPLORER_JS = r"""
(function(){const E=EXPLORER,I=INSPECT,$=id=>document.getElementById(id);
const TC=I.tierColors,esc=t=>String(t).replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const pill=t=>t?'<span class="tier" style="background:'+TC[t]+(t=='canonical'?'':';color:#1F1718')+'">'+esc(t)+'</span>':'';
// ---------------- inspector (right-hand drawer; every chart's points open it) ----------------
const ins=$('inspector');
function open(h){$('ins-body').innerHTML=h;ins.classList.add('open');ins.scrollTop=0;openedAt=Date.now()}
let openedAt=0;const close=()=>ins.classList.remove('open');$('ins-close').onclick=close;
// a click anywhere outside the drawer closes it (ignoring the click that just opened it)
document.addEventListener('click',e=>{if(ins.classList.contains('open')&&!ins.contains(e.target)&&Date.now()-openedAt>350)close()});
document.addEventListener('keydown',e=>{if(e.key=='Escape')ins.classList.remove('open')});
function link(href,txt){return '<a href="'+href+'" target="_blank" rel="noopener">'+esc(txt)+'</a>'}
function protView(i){if(!E)return;const runs=E.runs,v=E.val[i];
 const syms=(E.sym[i]||'').split(';').filter(Boolean),accs=(E.acc[i]||'').split(';').filter(Boolean);
 const m=v.filter(x=>x!=null),mean=m.length?m.reduce((a,b)=>a+b,0)/m.length:null;
 let h='<p class="eyebrow">Protein group</p><h2>'+esc(syms.join(' / ')||E.pid[i])+'</h2>'+pill(E.tier[i])+
  '<dl><dt>Group</dt><dd>'+esc(E.pid[i])+'</dd><dt>UniProt</dt><dd>'+accs.map(a=>link('https://www.uniprot.org/uniprotkb/'+a,a)).join(', ')+'</dd>'+
  (syms.length?'<dt>Surfaceome</dt><dd>'+syms.map(s=>link('https://surfaceome.deliverome.org/'+encodeURIComponent(s)+'/',s)).join(', ')+'</dd>':'')+
  '<dt>Precursors</dt><dd>'+E.nprec[i]+'</dd><dt>Seen in</dt><dd>'+m.length+' of '+runs.length+' runs</dd></dl>'+
  '<table><tr><th>Run</th><th>log2</th><th>vs mean</th></tr>'+runs.map((r,j)=>'<tr><td>'+esc(r)+'</td><td class="num">'+
  (v[j]==null?'&ndash;':v[j].toFixed(2))+'</td><td class="num">'+(v[j]==null||mean==null?'':(v[j]-mean>0?'+':'')+(v[j]-mean).toFixed(2))+'</td></tr>').join('')+'</table>'+
  '<p class="note">Values: '+esc(E.source)+'.</p>';
 const rows=[];E.prec_p.forEach((p,k)=>{if(p==i)rows.push(k)});
 if(E.detail&&rows.length)h+='<h3 class="sub">Precursors</h3><div id="ins-lines"></div><p class="note">Grey: each precursor '+
  '(gaps = not detected); red: the protein value. Hover a line for its precursor.</p>';
 open(h);
 if(E.detail&&rows.length){const lines=rows.map(k=>({type:'scatter',mode:'lines+markers',x:runs,y:E.prec_val[k],
  name:E.prec_id[k],line:{color:'#CFC4BF',width:1},marker:{size:4},connectgaps:false,
  hovertemplate:esc(E.prec_id[k])+'<br>%{x}: %{y:.2f}<extra></extra>'}));
  lines.push({type:'scatter',mode:'lines+markers',x:runs,y:v,name:'protein',line:{color:'#BC3C4C',width:3},marker:{size:7},
   hovertemplate:'protein<br>%{x}: %{y:.2f}<extra></extra>'});
  Plotly.newPlot('ins-lines',lines,Object.assign(baseLayout(),{height:260,margin:{l:44,r:8,t:6,b:40},showlegend:false,
   xaxis:{type:'category'},yaxis:{title:{text:'log2'},gridcolor:GRID}}),{responsive:true,displaylogo:false,displayModeBar:false})}
}
function runView(r){const R=I.runs[r];if(!R)return;
 let h='<p class="eyebrow">Run</p><h2>'+esc(R.label)+'</h2><p class="note mono">'+esc(R.run)+'</p>'+
  (R.flags?'<div class="flags">'+R.flags.split(',').map(f=>'<span class="flag">'+esc(f.trim())+'</span>').join('')+'</div><p>'+esc(R.hints)+'</p>':'')+
  '<dl>'+R.rows.map(([k,v])=>'<dt>'+esc(k)+'</dt><dd>'+esc(v)+'</dd>').join('')+'</dl>';
 open(h)}
function tierView(r,t){const rows=((I.top[r]||{})[t])||[];const idx=E?Object.fromEntries(E.pid.map((p,i)=>[p,i])):{};
 open('<p class="eyebrow">'+esc(I.runs[r].label)+'</p><h2>'+pill(t)+' proteins</h2><p class="note">Top '+rows.length+
  ' by share of the run\'s precursor signal. Click one for details.</p><table class="click"><tr><th>Protein</th><th>% of signal</th></tr>'+
  rows.map(([p,sh])=>'<tr data-pid="'+esc(p)+'"><td>'+esc(E&&idx[p]!=null?E.label[idx[p]]:p)+'</td><td class="num">'+sh.toFixed(3)+'</td></tr>').join('')+'</table>');
 $('ins-body').querySelector('table').onclick=e=>{const tr=e.target.closest('tr[data-pid]');if(tr)window.inspect(['pid',tr.dataset.pid])}}
window.inspect=c=>{if(!Array.isArray(c))return;
 if(c[0]=='run')runView(c[1]);else if(c[0]=='tier')tierView(c[1],c[2]);else if(c[0]=='prot')protView(c[1]);
 else if(c[0]=='pid'&&E){const i=E.pid.indexOf(c[1]);if(i>=0)protView(i)}};
for(const id of Object.keys(FIGS)){const el=$(id);if(el&&el.on)el.on('plotly_click',ev=>window.inspect(ev.points[0].customdata))}
// ---------------- protein explorer ----------------
if(!E||!$('ex-run'))return;const R=E.runs.length,N=E.pid.length,sel=$('ex-run');
E.runs.forEach((r,i)=>sel.add(new Option(r+(new Set(E.group).size>1?' ('+E.group[i]+')':''),i)));
$('ex-src').textContent='Protein values: '+E.source+'.';
const mean=a=>{let s=0,n=0;for(const v of a)if(v!=null){s+=v;n++}return n?s/n:null};
function refIdx(r){const same=[...Array(R).keys()].filter(j=>j!=r&&E.group[j]==E.group[r]);
 return same.length?same:[...Array(R).keys()].filter(j=>j!=r)}
let cur=null;const tiers=[...new Set(E.tier)].filter(Boolean);const order=Object.keys(TC).filter(t=>tiers.includes(t));
function dev(r){const ref=refIdx(r);return E.val.map(v=>{const o=mean(ref.map(j=>v[j]));
 return v[r]==null||o==null?null:v[r]-o})}
function draw(){const r=+sel.value,d=dev(r),q=$('ex-q').value.trim().toLowerCase(),only=$('ex-surf')&&$('ex-surf').checked;
 const x=E.val.map(mean);let idx=[...Array(N).keys()].filter(i=>d[i]!=null);
 if(only)idx=idx.filter(i=>E.tier[i]&&E.tier[i]!='not surface');
 const hit=q?idx.filter(i=>E.label[i].toLowerCase().includes(q)):[];
 const groups=order.length?order.map(t=>[t,idx.filter(i=>E.tier[i]==t)]):[['proteins',idx]];
 const trs=groups.map(([t,ii])=>({type:'scattergl',mode:'markers',x:ii.map(i=>x[i]),y:ii.map(i=>d[i]),
  customdata:ii.map(i=>['prot',i]),text:ii.map(i=>E.label[i]),name:t,
  marker:{size:t=='not surface'?4:6,color:TC[t]||'#3D6B60',opacity:t=='not surface'?.35:.85},
  hovertemplate:'%{text}<br>'+t+'<br>mean %{x:.2f}, this run %{y:+.2f}<extra></extra>'}));
 trs.push({type:'scattergl',mode:'markers',x:hit.map(i=>x[i]),y:hit.map(i=>d[i]),customdata:hit.map(i=>['prot',i]),
  text:hit.map(i=>E.label[i]),marker:{size:10,color:'#F4AA28',line:{color:'#1F1718',width:1}},
  hovertemplate:'%{text}<extra></extra>',name:'search',showlegend:false});
 const pk=cur==null||d[cur]==null?[]:[cur];
 trs.push({type:'scattergl',mode:'markers',x:pk.map(i=>x[i]),y:pk.map(i=>d[i]),customdata:pk.map(i=>['prot',i]),
  text:pk.map(i=>E.label[i]),marker:{size:14,color:'#BC3C4C',line:{color:'#1F1718',width:1.5}},
  hovertemplate:'%{text}<extra></extra>',name:'selected',showlegend:false});
 Plotly.react('ex-scatter',trs,Object.assign(baseLayout(),{height:440,margin:{l:60,r:10,t:10,b:50},showlegend:order.length>0,
  legend:{orientation:'h',y:1.08},
  xaxis:{title:{text:'Mean protein log2 across runs'},gridcolor:GRID},
  yaxis:{title:{text:'This run minus other runs (log2)'},gridcolor:GRID,zerolinecolor:'#80706A'}}),{responsive:true,displaylogo:false});
 const pool=(q?hit:idx).slice().sort((a,b)=>Math.abs(d[b])-Math.abs(d[a])).slice(0,30);
 $('ex-cap').textContent=(q?pool.length+' match(es) for "'+q+'"':'30 '+(only?'surface ':'')+'proteins furthest from the other runs')+
  ' in run '+E.runs[r]+'. Click a row or a point.';
 $('ex-tbl').innerHTML='<tr><th>Protein</th><th>Tier</th><th>Prec.</th><th>This run &minus; others</th></tr>'+
  pool.map(i=>'<tr data-i="'+i+'"'+(i==cur?' class="on"':'')+'><td>'+esc(E.label[i])+'</td><td>'+pill(E.tier[i])+'</td><td class="num">'+E.nprec[i]+
  '</td><td class="num">'+(d[i]>0?'+':'')+d[i].toFixed(2)+'</td></tr>').join('');
}
function detail(i){cur=i;draw();const rows=[];E.prec_p.forEach((p,k)=>{if(p==i)rows.push(k)});
 if(!E.detail){$('ex-detail').innerHTML='<h3>'+esc(E.label[i])+'</h3><p class="note">Precursor detail '+
  'is not embedded for a dataset this large.</p>';return}
 const runs=E.runs;const lines=rows.map(k=>({type:'scatter',mode:'lines+markers',x:runs,y:E.prec_val[k],
  name:E.prec_id[k],line:{color:'#CFC4BF',width:1},marker:{size:4},connectgaps:false,
  hovertemplate:esc(E.prec_id[k])+'<br>%{x}: %{y:.2f}<extra></extra>'}));
 lines.push({type:'scatter',mode:'lines+markers',x:runs,y:E.val[i],name:'protein',
  line:{color:'#BC3C4C',width:3},marker:{size:8},hovertemplate:'protein<br>%{x}: %{y:.2f}<extra></extra>'});
 $('ex-detail').innerHTML='<h3>'+esc(E.label[i])+' '+pill(E.tier[i])+'</h3><p class="note">'+rows.length+
  ' precursors (grey; gaps = not detected) and the protein value (red). Hover a line for its precursor.</p>'+
  '<div id="ex-lines"></div><div class="scroll ex-table"><table id="ex-ptbl"></table></div>';
 Plotly.newPlot('ex-lines',lines,Object.assign(baseLayout(),{height:380,margin:{l:60,r:10,t:10,b:60},showlegend:false,
  xaxis:{type:'category'},yaxis:{title:{text:'log2 intensity'},gridcolor:GRID}}),{responsive:true,displaylogo:false});
 $('ex-ptbl').innerHTML='<tr><th>Precursor</th>'+runs.map(r=>'<th>'+esc(r)+'</th>').join('')+'<th>RT flagged in</th></tr>'+
  rows.map(k=>'<tr'+(E.prec_flag[k]?' class="rtflag"':'')+'><td>'+esc(E.prec_id[k])+'</td>'+E.prec_val[k].map(v=>'<td class="num">'+(v==null?'&ndash;':v.toFixed(2))+'</td>').join('')+'<td>'+esc(E.prec_flag[k]||'')+'</td></tr>').join('');
}
sel.onchange=draw;$('ex-q').oninput=draw;if($('ex-surf'))$('ex-surf').onchange=draw;
$('ex-tbl').onclick=e=>{const tr=e.target.closest('tr[data-i]');if(tr){detail(+tr.dataset.i);protView(+tr.dataset.i)}};
draw();$('ex-scatter').on('plotly_click',ev=>{const c=ev.points[0].customdata;if(c&&c[0]=='prot'){detail(c[1]);protView(c[1])}});

})();
"""


# ----------------------------------------------------------------------------------------------
# page
# ----------------------------------------------------------------------------------------------

_CSS = """
/* Design tokens from the accessible-surfaceome viewer (viewer/app/design-tokens.css). */
:root{--ink:#1f1718;--muted:#6f5d5a;--muted-soft:#80706a;--line:rgba(31,23,24,.1);--line-soft:rgba(31,23,24,.06);
--bg:#fefefe;--bg-warm:#f3ece5;--card:#fff;--accent:#bc3c4c;--accent-deep:#6e1428;--accent-blush:#fde8e6;
--teal:#3d6b60;--teal-tint:#cce8e4;--amber:#c07830;--amber-tint:#faecd4;--lav:#5848a8;--lav-tint:#e4e0f8;
--ok:#2e7a55;--warn:#8c4210;--link:#3a2888;--radius:1rem;--radius-sm:.5rem;
--shadow:0 4px 12px rgba(31,23,24,.05);--shadow-lg:0 24px 60px rgba(31,23,24,.1);
--sans:"Manrope","Inter",system-ui,sans-serif;--display:"Playfair Display",Georgia,serif}
*{box-sizing:border-box}html{scroll-behavior:smooth}
body{margin:0;background:var(--bg);color:var(--ink);font:400 14px/1.55 var(--sans);-webkit-font-smoothing:antialiased}
a{color:var(--link)}a:hover{color:var(--accent-deep)}
header{background:linear-gradient(180deg,var(--bg-warm),var(--bg));padding:36px 32px 18px;border-bottom:1px solid var(--line-soft)}
.eyebrow{font-size:.72rem;letter-spacing:.14em;text-transform:uppercase;color:var(--accent-deep);font-weight:600;margin:0 0 6px}
h1{font:500 2.3rem/1.1 var(--display);letter-spacing:-.02em;margin:0 0 8px}
.meta{color:var(--muted);font-size:13px}
nav.toc{position:sticky;top:0;z-index:5;background:rgba(254,254,254,.92);backdrop-filter:blur(8px);
border-bottom:1px solid var(--line-soft);padding:8px 32px;display:flex;gap:6px;flex-wrap:wrap}
nav.toc a{font-size:12.5px;text-decoration:none;color:var(--muted);padding:4px 10px;border-radius:999px}
nav.toc a:hover{background:var(--accent-blush);color:var(--accent-deep)}
main{padding:8px 32px 60px;max-width:1500px}
section{scroll-margin-top:52px}
h2{font:600 1.3rem/1.25 var(--sans);letter-spacing:-.01em;margin:34px 0 12px}
.grid{display:grid;grid-template-columns:repeat(var(--cols),minmax(0,1fr));gap:16px}
.card{background:var(--card);border:1px solid var(--line);border-radius:var(--radius);padding:16px 18px;box-shadow:var(--shadow)}
.card h3{font-size:14.5px;font-weight:600;margin:0 0 6px}
.note{color:var(--muted);font-size:12.5px;margin:6px 0 0}
.gap{display:flex;flex-direction:column;justify-content:center;min-height:150px;box-shadow:none;
background:repeating-linear-gradient(135deg,#fff,#fff 10px,#faf6f3 10px,#faf6f3 20px)}
.pill,.tier,.flag{display:inline-block;font-size:11px;font-weight:600;padding:2px 9px;border-radius:999px;margin-right:4px}
.tier{color:#fff}.flag{background:var(--amber-tint);color:var(--warn)}.ok{color:var(--ok);font-size:12px}
.shown{background:#e4f1ea;color:var(--ok)}.missing{background:var(--amber-tint);color:var(--warn)}
.planned{background:var(--lav-tint);color:var(--lav)}
table{border-collapse:collapse;width:100%;font-size:12.5px;background:var(--card)}
th,td{border-bottom:1px solid var(--line-soft);padding:6px 9px;text-align:left;vertical-align:top}
th{color:var(--muted);font-weight:600;background:#fcfaf8;position:sticky;top:0}
td.num{text-align:right;font-variant-numeric:tabular-nums}td.hint{color:var(--muted);min-width:260px}
.scroll{overflow:auto;border:1px solid var(--line);border-radius:var(--radius-sm)}
details summary{cursor:pointer;color:var(--muted);margin:10px 0}
.wide{grid-column:1/-1}.controls{display:flex;gap:16px;align-items:center;flex-wrap:wrap;margin:8px 0}
.controls select,.controls input[type=text],.controls input:not([type]){font:inherit;padding:5px 9px;border:1px solid var(--line);border-radius:var(--radius-sm)}
.chk{display:flex;gap:6px;align-items:center;color:var(--muted)}
.ex-grid{display:grid;grid-template-columns:3fr 2fr;gap:16px}.ex-table{max-height:420px;overflow:auto}
tr[data-i],tr[data-pid],tr[data-run]{cursor:pointer}
tr[data-i]:hover,tr[data-pid]:hover,tr[data-run]:hover,tr.on{background:var(--accent-blush)}
tr.rtflag td{color:#9c0006}
.ex-detail{margin-top:14px;border-top:1px solid var(--line-soft);padding-top:12px}
#inspector{position:fixed;top:0;right:0;height:100vh;width:min(420px,92vw);background:var(--card);z-index:20;
border-left:1px solid var(--line);box-shadow:var(--shadow-lg);padding:22px 24px 40px;overflow:auto;
transform:translateX(105%);transition:transform .22s ease}
#inspector.open{transform:none}
#inspector h2{margin:4px 0 10px}#ins-close{position:absolute;top:14px;right:16px;border:0;background:none;
font-size:22px;cursor:pointer;color:var(--muted)}
#inspector dl{display:grid;grid-template-columns:auto 1fr;gap:4px 14px;font-size:13px;margin:14px 0}
#inspector dt{color:var(--muted)}#inspector dd{margin:0;font-variant-numeric:tabular-nums}
.mono{font-family:ui-monospace,Menlo,monospace;font-size:11.5px;word-break:break-all}
.hint-bar{font-size:12.5px;color:var(--muted);margin:10px 0 0}
.info{display:inline-flex;align-items:center;justify-content:center;width:16px;height:16px;border-radius:50%;
border:1px solid var(--muted-soft);color:var(--muted);font:600 10.5px var(--sans);cursor:help;position:relative;
vertical-align:1px;margin-left:4px}
.info:hover,.info:focus{background:var(--accent-blush);color:var(--accent-deep);border-color:var(--accent);outline:none}
.info:hover::after,.info:focus::after{content:attr(data-tip);position:absolute;left:0;top:22px;z-index:30;
width:min(340px,80vw);background:#fff;color:var(--ink);border:1px solid var(--line);border-radius:var(--radius-sm);
box-shadow:var(--shadow-lg);padding:10px 12px;font:400 12.5px/1.5 var(--sans);text-align:left;white-space:normal}
.card:has(.info:hover),.card:has(.info:focus){z-index:10;position:relative}
#inspector h3.sub{font-size:13px;margin:16px 0 4px}
@media (max-width:860px){.grid,.ex-grid{grid-template-columns:1fr}header,main,nav.toc{padding-left:16px;padding-right:16px}}
"""


def _info(name: str) -> str:
    text = INFO.get(name)
    if not text:
        return ""
    return (
        f"<span class='info' tabindex='0' role='note' aria-label='{html.escape(text)}'"
        f" data-tip='{html.escape(text)}'>i</span>"
    )


def _fmt(v) -> str:
    if isinstance(v, float | np.floating):
        if not math.isfinite(v):
            return ""
        return f"{v:,.0f}" if abs(v) >= 1000 else f"{v:.3g}"
    if isinstance(v, int | np.integer):
        return f"{v:,}"
    return html.escape(str(v))


def _table(df: pd.DataFrame, index_name: str = "") -> str:
    head = "".join(f"<th>{html.escape(str(c))}</th>" for c in [index_name, *df.columns])
    rows = []
    for idx, row in df.iterrows():
        cells = [f"<td>{html.escape(str(idx))}</td>"]
        for v in row:
            num = isinstance(v, int | float | np.integer | np.floating)
            cells.append(f"<td class='num'>{_fmt(v)}</td>" if num else f"<td>{_fmt(v)}</td>")
        rows.append("<tr>" + "".join(cells) + "</tr>")
    return f"<div class='scroll'><table><tr>{head}</tr>{''.join(rows)}</table></div>"


def _contam_table(ctx: _Ctx) -> str:
    top = qc.top_contaminants(ctx.prec)
    if top.empty:
        return ""
    rows = "".join(
        f"<tr data-pid='{html.escape(str(r.protein))}'><td>{html.escape(ctx.protein_name(r.protein))}</td>"
        f"<td>{html.escape(r.protein_class)}</td><td class='num'>{r.median_signal_pct:.3f}</td>"
        f"<td class='num'>{r.runs}</td></tr>"
        for r in top.itertuples()
    )
    return (
        "<div class='scroll ex-table'><table id='contam-tbl'><tr><th>Protein group</th>"
        "<th>Class</th><th>Median % of run signal</th><th>Runs</th></tr>"
        f"{rows}</table></div><script>document.getElementById('contam-tbl').onclick=e=>{{"
        "const t=e.target.closest('tr[data-pid]');if(t&&window.inspect)"
        "window.inspect(['pid',t.dataset.pid])};</script>"
    )


_DIAG_COLS = [
    ("precursors", "Precursors", "{:,.0f}"),
    ("tic_area_rel", "TIC area vs others", "{:.0f}%"),
    ("roughness", "TIC roughness", "{:.3f}"),
    ("dropouts", "Dropout scans", "{:,.0f}"),
    ("spray_noise_pct", "Spray noise", "{:.1f}%"),
    ("pressure_dev_pct", "Pressure deviation", "{:.1f}%"),
    ("fwhm_s", "Median FWHM", "{:.1f} s"),
    ("tailing", "Median tailing", "{:.2f}"),
    ("rt_shift_min", "RT shift", "{:+.2f} min"),
    ("tic_centroid_rt", "TIC centroid", "{:.1f} min"),
    ("inject_ms", "MS1 fill time", "{:.2f} ms"),
    ("ms2_inject_ms", "MS2 fill time", "{:.2f} ms"),
    ("cycle_time_s", "Cycle time", "{:.2f} s"),
    ("pump_pressure_max", "Max pump pressure", "{:.0f}"),
]


def _diagnosis_frame(ctx: _Ctx) -> pd.DataFrame:
    d = ctx.diagnosis.copy()
    if d.empty:
        return d
    med = d["tic_area"].median()
    d["tic_area_rel"] = 100 * d["tic_area"] / med if med else np.nan
    return d


def _diagnosis_table(ctx: _Ctx) -> str:
    d = _diagnosis_frame(ctx)
    if d.empty:
        return ""
    cols = [c for c in _DIAG_COLS if c[0] in d and d[c[0]].notna().any()]
    head = "".join(f"<th>{h}</th>" for _, h, _ in cols)
    rows = []
    for run, r in d.iterrows():
        cells = "".join(
            f"<td class='num'>{f.format(r[c]) if pd.notna(r[c]) else '&ndash;'}</td>"
            for c, _, f in cols
        )
        flags = (
            "".join(
                f"<span class='flag'>{html.escape(f.strip())}</span>"
                for f in r["flags"].split(",")
                if f.strip()
            )
            or "<span class='ok'>no flags</span>"
        )
        rows.append(
            f"<tr data-run='{ctx.run_index(run)}'><td>{html.escape(ctx.labels[run])}</td><td>{flags}</td>"
            f"{cells}<td class='hint'>{html.escape(r['hints'])}</td></tr>"
        )
    return (
        "<div class='scroll'><table id='diag-tbl'><tr><th>Run</th><th>Flags</th>"
        f"{head}<th>What it suggests</th></tr>{''.join(rows)}</table></div>"
        "<script>document.getElementById('diag-tbl').onclick=e=>{const t=e.target.closest("
        "'tr[data-run]');if(t&&window.inspect)window.inspect(['run',+t.dataset.run])};</script>"
    )


def _rt_flag_table(ctx: _Ctx, n: int = 50) -> str:
    rtc = ctx.rt[0]
    if rtc.empty:
        return ""
    f = rtc[rtc["flag"]].copy()
    f["resid_s"] = 60 * f["resid_min"]
    f = f.reindex(f["resid_s"].abs().sort_values(ascending=False).index).head(n)
    head = "".join(
        f"<th>{c}</th>"
        for c in ("Precursor", "Protein", "Run", "RT (min)", "Residual (s)", "Flagged in")
    )
    n_runs = rtc[rtc["flag"]].groupby("precursor")["run"].nunique()
    rows = []
    for r in f.itertuples():
        raw = str(getattr(r, "protein", ""))
        pid = html.escape(raw)
        rows.append(
            f"<tr data-pid='{pid}'><td>{html.escape(r.precursor)}</td>"
            f"<td>{html.escape(ctx.protein_name(raw))}</td>"
            f"<td>{html.escape(ctx.labels[r.run])}</td><td class='num'>{r.ref_rt:.2f}</td>"
            f"<td class='num'>{r.resid_s:+.1f}</td><td class='num'>{n_runs[r.precursor]} run(s)</td></tr>"
        )
    total = int(rtc["flag"].sum())
    return (
        f"<p class='note'>{total:,} precursor-run values flagged; the {len(f)} largest shown."
        "</p><div class='scroll ex-table'><table id='rt-tbl'><tr>"
        f"{head}</tr>{''.join(rows)}</table></div>"
        "<script>document.getElementById('rt-tbl').onclick=e=>{const t=e.target.closest("
        "'tr[data-pid]');if(t&&window.inspect)window.inspect(['pid',t.dataset.pid])};</script>"
    )


def build_dashboard_html(
    prec: pd.DataFrame,
    layout: str | os.PathLike | dict = "default",
    *,
    title: str | None = None,
    source: str = "",
    groups: pd.Series | None = None,
    limpa_result=None,
    names: dict | None = None,
    ann: pd.DataFrame | None = None,
    chrom: dict | None = None,
    history: pd.DataFrame | None = None,
) -> tuple[str, pd.DataFrame]:
    """HTML for one layout plus the coverage table (section, panel, status, reason).

    ``names`` optionally maps protein id -> display name for the explorer (display only;
    nothing joins on it).
    """
    spec = layout if isinstance(layout, dict) else load_layout(layout)
    ctx = _Ctx(prec, groups, limpa_result, ann=ann, chrom=chrom, history=history)
    figs: dict[str, dict] = {}
    explorer = None
    cov, sections_html = [], []
    for si, sec in enumerate(spec["sections"]):
        cards = []
        for pi, name in enumerate(sec["panels"]):
            m = METRICS[name]
            status, reason = panel_status(prec, name, limpa_result, ann=ann, chrom=chrom)
            if m.kind == "explorer":
                if len(ctx.order) < 2:
                    status, reason = "missing", "Needs at least 2 runs."
                else:
                    explorer = _explorer_data(ctx, limpa_result, names)
                    cov.append(
                        {
                            "section": sec["title"],
                            "panel": m.title,
                            "status": status,
                            "reason": reason,
                        }
                    )
                    cards.append(
                        _EXPLORER_HTML.format(
                            note=f"<p class='note'>{html.escape(m.note)}</p>", info=_info(name)
                        )
                    )
                    continue
            if m.kind == "contamtable" and status == "shown":
                table = _contam_table(ctx) if _has_classes(ctx) else ""
                if table:
                    cov.append(
                        {
                            "section": sec["title"],
                            "panel": m.title,
                            "status": status,
                            "reason": reason,
                        }
                    )
                    cards.append(
                        f"<div class='card wide'><h3>{html.escape(m.title)} {_info(name)}</h3>"
                        f"<p class='note'>{html.escape(m.note)}</p>{table}</div>"
                    )
                    continue
                status, reason = "missing", _EMPTY_REASON["contam"]
            if m.kind == "diagnosis" and status == "shown":
                table = _diagnosis_table(ctx)
                if table:
                    cov.append(
                        {
                            "section": sec["title"],
                            "panel": m.title,
                            "status": status,
                            "reason": reason,
                        }
                    )
                    cards.append(
                        f"<div class='card wide'><h3>{html.escape(m.title)} {_info(name)}</h3>"
                        f"<p class='note'>{html.escape(m.note)}</p>{table}</div>"
                    )
                    continue
                status, reason = "missing", "No MS1 scans in the chromatography tables."
            if m.kind == "rtflag" and status == "shown":
                table = _rt_flag_table(ctx)
                if table:
                    cov.append(
                        {
                            "section": sec["title"],
                            "panel": m.title,
                            "status": status,
                            "reason": reason,
                        }
                    )
                    cards.append(
                        f"<div class='card wide'><h3>{html.escape(m.title)} {_info(name)}</h3>"
                        f"<p class='note'>{html.escape(m.note)}</p>{table}</div>"
                    )
                    continue
                status, reason = "missing", "Needs at least 2 runs."
            fig = (
                _BUILDERS[m.kind](ctx, m)
                if status == "shown"
                and m.kind not in ("explorer", "rtflag", "contamtable", "diagnosis")
                else None
            )
            if status == "shown" and fig is None:
                status, reason = "missing", _EMPTY_REASON.get(m.kind, "No data for this panel.")
            cov.append(
                {"section": sec["title"], "panel": m.title, "status": status, "reason": reason}
            )
            text = " ".join(t for t in (m.note, _BOX_NOTE if m.kind == "dist" else "") if t)
            note = f"<p class='note'>{html.escape(text)}</p>" if text else ""
            if fig is None:
                cards.append(
                    f"<div class='card gap'><h3>{html.escape(m.title)} {_info(name)} "
                    f"<span class='pill {status}'>{status}</span></h3>"
                    f"<p class='note'>{html.escape(reason)}</p>{note}</div>"
                )
                continue
            div = f"p{si}_{pi}"
            figs[div] = _clean(fig)
            cards.append(
                f"<div class='card'><h3>{html.escape(m.title)} {_info(name)}</h3><div id='{div}'></div>{note}</div>"
            )
        sections_html.append(
            f"<section id='s{si}'><h2>{html.escape(sec['title'])}</h2><div class='grid'>{''.join(cards)}</div>"
            "</section>"
        )
    coverage = pd.DataFrame(cov)
    n_shown = int((coverage["status"] == "shown").sum())

    summary = ctx.summary.copy()
    summary.index = [ctx.labels[r] for r in summary.index]
    summary = summary.drop(columns=[c for c in ("condition",) if summary[c].nunique() < 2])
    cov_tbl = coverage.assign(
        status=coverage["status"].map(lambda s: f"<span class='pill {s}'>{s}</span>")
    )
    cov_rows = "".join(
        f"<tr><td>{html.escape(r.section)}</td><td>{html.escape(r.panel)}</td>"
        f"<td>{r.status}</td><td>{html.escape(r.reason)}</td></tr>"
        for r in cov_tbl.itertuples()
    )
    if explorer is None and len(ctx.order) >= 2:
        explorer = _explorer_data(ctx, limpa_result, names)  # the inspector uses it too
    inspect = _inspect_data(ctx)
    run_rows = "".join(
        f"<tr data-run='{i}'><td>{html.escape(str(lab))}</td>"
        + "".join(
            f"<td class='num'>{_fmt(v)}</td>"
            if isinstance(v, int | float | np.integer | np.floating)
            else f"<td>{_fmt(v)}</td>"
            for v in summary.iloc[i]
        )
        + "</tr>"
        for i, lab in enumerate(summary.index)
    )
    run_head = "".join(f"<th>{html.escape(str(c))}</th>" for c in ["run", *summary.columns])
    nav = "".join(
        f"<a href='#s{si}'>{html.escape(sec['title'])}</a>"
        for si, sec in enumerate(spec["sections"])
    )
    engine = ""
    if limpa_result is not None:
        e = str(limpa_result.summary.get("engine") or "R limpa")
        e = "limpa-rs (Rust)" if "rust" in e.lower() else e
        engine = f" · protein quantification: {html.escape(e)} {html.escape(str(limpa_result.summary.get('limpa_version', '')))}"
    page_title = title or spec["title"]
    body = f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>{html.escape(page_title)}</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link href="https://fonts.googleapis.com/css2?family=Manrope:wght@400;500;600;700&family=Playfair+Display:wght@500&display=swap" rel="stylesheet">
<script src="{PLOTLY_JS}"></script><style>{_CSS}</style></head>
<body style="--cols:{int(spec["columns"])}"><header><p class="eyebrow">Deliverome · DIA QC</p>
<h1>{html.escape(page_title)}</h1>
<div class="meta">{html.escape(source)} · {len(ctx.order)} runs · {prec["precursor"].nunique():,}
precursors · {n_shown} of {len(coverage)} panels shown{engine}</div>
<p class="hint-bar">Click any point, bar or table row to open its details on the right.</p></header>
<nav class="toc"><a href="#runs">Runs</a>{nav}</nav><main>
<section id="runs"><h2>Runs</h2><div class='scroll'><table id='run-tbl'><tr>{run_head}</tr>{run_rows}</table></div>
<details><summary>Coverage: what this export can and cannot show</summary>
<div class="scroll"><table><tr><th>Section</th><th>Panel</th><th>Status</th><th>Needs</th></tr>
{cov_rows}</table></div></details></section>
{"".join(sections_html)}</main>
<aside id="inspector" aria-label="Details"><button id="ins-close" aria-label="Close">&times;</button>
<div id="ins-body"></div></aside>
<script>const GRID={json.dumps(GRID)};
const baseLayout=()=>({{font:{{family:'Manrope, Inter, sans-serif',size:12,color:'#1f1718'}},paper_bgcolor:'rgba(0,0,0,0)',
plot_bgcolor:'rgba(0,0,0,0)',hovermode:'closest',hoverlabel:{json.dumps(HOVER)}}});
const FIGS={json.dumps(figs)};
for (const [id,f] of Object.entries(FIGS)) Plotly.newPlot(id,f.data,f.layout,
{{responsive:true,displaylogo:false}});
const EXPLORER={json.dumps(explorer)};const INSPECT={json.dumps(inspect)};{_EXPLORER_JS}
document.getElementById('run-tbl').onclick=e=>{{const t=e.target.closest('tr[data-run]');
if(t&&window.inspect)window.inspect(['run',+t.dataset.run])}};</script></body></html>"""
    return body, coverage


# ----------------------------------------------------------------------------------------------
# run metrics and the QC record (one JSON per page, indexed by deliverome_analysis.qcdb)
# ----------------------------------------------------------------------------------------------


def run_metrics(ctx: _Ctx) -> pd.DataFrame:
    """Long table (run, metric, value) of every per-run number this page knows."""
    p = ctx.prec
    rows, seen = [], set()

    def add(run, metric, value):  # first source wins: e.g. precursors is in two tables
        if value is not None and pd.notna(value) and (run, metric) not in seen:
            seen.add((run, metric))
            rows.append((run, metric, float(value)))

    s = ctx.summary
    for col in s.columns:
        if pd.api.types.is_numeric_dtype(s[col]):
            for run, v in s[col].items():
                add(run, col, v)
    g = p.groupby("run")
    per_run = {"summed_intensity": g["intensity"].sum()}
    if "apex_rt" in p:
        per_run["median_precursor_rt"] = g["apex_rt"].median()
    if "peptide" in p:
        plen = p.drop_duplicates(["run", "peptide"]).assign(n=lambda d: d["peptide"].str.len())
        per_run["peptide_length_mean"] = plen.groupby("run")["n"].mean()
        per_run["peptide_length_median"] = plen.groupby("run")["n"].median()
    z = pd.to_numeric(p["precursor"].str.rsplit(".", n=1).str[-1], errors="coerce")
    per_run["charge_mean"] = z.groupby(p["run"]).mean()
    if "missed_cleavages" in p:
        per_run["missed_cleavage_pct"] = 100 * (p["missed_cleavages"] > 0).groupby(p["run"]).mean()
    per_run["completeness_pct"] = 100 * g["precursor"].nunique() / p["precursor"].nunique()
    for metric, series in per_run.items():
        for run, v in series.items():
            add(run, metric, v)
    if ctx.ann is not None:
        surf = ctx.surface
        for key in ("protein_groups", "peptides", "signal_ms2", "signal_ms1", "missed"):
            if key in surf:
                for t in qc.TIER_CLASSES:
                    name = f"{key}_{t.replace(' / ', '_').replace(' ', '_')}"
                    for run, v in surf[key][t].items():
                        add(run, name, v)
    if ctx.chrom and ctx.chrom.get("scans") is not None:
        d = ctx.diagnosis
        for col in d.columns:
            if pd.api.types.is_numeric_dtype(d[col]):
                for run, v in d[col].items():
                    add(run, col, v)
    return pd.DataFrame(rows, columns=["run", "metric", "value"])


def qc_record(
    ctx: _Ctx,
    *,
    source: str,
    title: str | None,
    layout: str | os.PathLike,
    output: Path,
    limpa_result=None,
) -> dict:
    """Everything the shared QC index needs about this page: runs, metrics, YAML, versions."""
    import hashlib
    import importlib.metadata as md
    import time

    from deliverome_analysis.provenance import provenance

    spec_path = Path(layout) if Path(layout).suffix else LAYOUT_DIR / f"{layout}.yaml"
    yaml_text = spec_path.read_text()
    created = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    found = (ctx.chrom or {}).get("found") or {}
    info_keys = (
        "acquired",
        "instrument_model",
        "instrument_name",
        "instrument_serial",
        "file_bytes",
    )
    runs = []
    for r in ctx.order:
        info = (found.get(r) or {}).get("info") or {}
        runs.append(
            {"run": r, "label": ctx.labels[r], "condition": str(ctx.cond[r]),
             **{k: info.get(k) for k in info_keys}}
        )  # fmt: skip
    diagnosis = []
    if ctx.chrom and ctx.chrom.get("scans") is not None:
        d = ctx.diagnosis
        diagnosis = [
            {"run": r, "flags": d.at[r, "flags"], "hints": d.at[r, "hints"]} for r in d.index
        ]
    lim = None
    if limpa_result is not None:
        lim = {
            k: limpa_result.summary.get(k)
            for k in ("engine", "limpa_version", "n_proteins", "dpc_beta0", "dpc_beta1", "n_runs")
        }
    m = run_metrics(ctx)
    return _clean(
        {
            "schema_version": 1,  # qcdb.RECORD_SCHEMA_VERSION
            "kind": "dia_qc",
            "record_id": hashlib.sha256(f"{source}|{output.name}|{created}".encode()).hexdigest()[
                :16
            ],
            "created": created,
            "source": source,
            "title": title,
            "page": output.name,
            "layout": str(layout),
            "layout_yaml": yaml_text,
            "layout_sha256": hashlib.sha256(yaml_text.encode()).hexdigest(),
            "provenance": provenance(),
            "versions": {
                pkg: md.version(pkg) for pkg in ("deliverome-analysis", "pandas", "numpy", "duckdb")
            },
            "limpa": lim,
            "runs": runs,
            "metrics": [list(t) for t in m.itertuples(index=False)],
            "diagnosis": diagnosis,
        }
    )


def _report_name(path) -> str:
    """File name without the S3 cache's ``.<etag>`` suffix (``x.tsv.3fe6...-112`` -> ``x.tsv``)."""
    name = Path(path).name
    for ext in (".tsv", ".csv", ".txt"):
        if ext + "." in name:
            return name[: name.index(ext + ".") + len(ext)]
    return name


def run_record(
    *,
    source: str,
    layout: str | os.PathLike,
    coverage: pd.DataFrame,
    output: Path,
    matrix_dir=None,
    limpa_result=None,
) -> dict:
    """What produced a dashboard page: input, layout (with a hash of its YAML), code version."""
    import hashlib
    import importlib.metadata as md
    import subprocess
    import time

    spec_path = Path(layout) if Path(layout).suffix else LAYOUT_DIR / f"{layout}.yaml"
    try:
        commit = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            cwd=Path(__file__).parent,
            check=True,
        ).stdout.strip()
        dirty = bool(
            subprocess.run(
                ["git", "status", "--porcelain", "--", "."],  # this package dir (cwd)
                capture_output=True,
                text=True,
                cwd=Path(__file__).parent,
                check=True,
            ).stdout.strip()
        )
    except OSError, subprocess.CalledProcessError:
        commit, dirty = None, None
    rec = {
        "created": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "output": output.name,
        "source": source,
        "matrix_dir": str(matrix_dir) if matrix_dir else None,
        "layout": str(layout),
        "layout_sha256": hashlib.sha256(spec_path.read_bytes()).hexdigest(),
        "git_commit": commit,
        "src_uncommitted_changes": dirty,
        "versions": {p: md.version(p) for p in ("deliverome-analysis", "pandas", "numpy")},
        "plotly_js": PLOTLY_JS,
        "limpa": (limpa_result.summary if limpa_result is not None else None),
        "panels": coverage.to_dict(orient="records"),
    }
    return rec


def dia_qc_dashboard(
    matrix_dir: str | os.PathLike | None = None,
    report: str | os.PathLike | None = None,
    *,
    layout: str | os.PathLike = "default",
    output: str | os.PathLike | None = None,
    title: str | None = None,
    groups: pd.Series | None = None,
    q_cutoff: float = 0.01,
    open_browser: bool = True,
    limpa_result=None,
    names: dict | None = None,
    chromatography_dir: str | os.PathLike | None = None,
    surface: bool | pd.DataFrame = True,
    history: bool | pd.DataFrame = True,
) -> tuple[Path, pd.DataFrame]:
    """Write the dashboard for one layout and (by default) open it in the browser.

    Give a ``trim_spectronaut_report`` directory (preferred: reads only its small QC table)
    or a report. Returns the HTML path and the coverage table.

    ``chromatography_dir``: output of ``scripts/raw_chromatography_modal.py`` (TIC, spray and
    LC traces from the raw files); also writes ``<page>.run_diagnostics.tsv`` with blank
    ``confirmed_cause`` / ``notes`` columns for recording what a flagged run turned out to be.
    ``surface``: annotate protein groups with accessible-surfaceome tiers (needs S3), or pass
    a ready :func:`deliverome_analysis.dia_qc.surface_annotation` table, or False.
    ``history``: earlier runs for the "Runs over time" panel, from the shared QC index
    (:mod:`deliverome_analysis.qcdb`, needs S3), or a ``latest_run_metrics`` frame, or False.
    Also writes ``<page>.qcrecord.json``, the record the index is built from.
    """
    if (report is None) == (matrix_dir is None):
        raise ValueError("Give exactly one of matrix_dir= or report=")
    if matrix_dir is not None:
        prec, source = qc.load_trimmed_precursors(matrix_dir), str(matrix_dir)
        man = Path(matrix_dir) / "manifest.json"
        if man.exists():
            source = _report_name(json.loads(man.read_text())["source_report"])
    else:
        prec = qc.load_spectronaut_precursors(report, q_cutoff=q_cutoff)
        source = _report_name(report)
    ann = None
    if isinstance(surface, pd.DataFrame):
        ann = surface
    elif surface:
        try:
            ann = qc.surface_annotation(prec)
        except Exception as exc:  # no S3 credentials, offline: surface panels say what is missing
            print(
                f"surface tiers unavailable ({type(exc).__name__}: {exc}); skipping surface panels"
            )
    hist = history if isinstance(history, pd.DataFrame) else None
    if history is True:
        try:
            from deliverome_analysis import qcdb

            con = qcdb.connect("dia_qc")
            hist = pd.read_sql("select * from latest_run_metrics", con)
            con.close()
        except Exception as exc:  # no index yet, or offline: the panel shows this page only
            print(f"QC history unavailable ({type(exc).__name__}); showing this page's runs only")
    chrom = None
    if chromatography_dir is not None:
        from deliverome_analysis import raw_chromatography as rc

        chrom = rc.load_chromatography(chromatography_dir)
    layout_name = Path(layout).stem
    out = Path(output) if output else Path(f"dia_qc_{layout_name}.html")
    page, coverage = build_dashboard_html(
        prec,
        layout,
        title=title,
        source=source,
        groups=groups,
        limpa_result=limpa_result,
        names=names,
        ann=ann,
        chrom=chrom,
        history=hist,
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(page)
    ctx = _Ctx(prec, groups, limpa_result, ann=ann, chrom=chrom)
    if chrom and chrom.get("scans") is not None:
        diag = _diagnosis_frame(ctx)
        if not diag.empty:
            diag.assign(confirmed_cause="", notes="").to_csv(
                out.with_suffix(".run_diagnostics.tsv"),
                sep="\t",
                index_label="run",
                float_format="%.4g",
            )
    # The record the shared QC index (deliverome_analysis.qcdb) is built from.
    record = qc_record(
        ctx, source=source, title=title, layout=layout, output=out, limpa_result=limpa_result
    )
    out.with_suffix(".qcrecord.json").write_text(json.dumps(record, indent=1, default=str) + "\n")
    rec = run_record(
        source=source,
        layout=layout,
        coverage=coverage,
        output=out,
        matrix_dir=matrix_dir,
        limpa_result=limpa_result,
    )
    out.with_suffix(".run.json").write_text(json.dumps(rec, indent=2, default=str) + "\n")
    if open_browser:
        webbrowser.open(out.resolve().as_uri())
    return out, coverage


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("matrix_dir", help="trim_spectronaut_report output directory")
    ap.add_argument("--layout", default="default", help=f"one of {available_layouts()} or a path")
    ap.add_argument("--out", default=None)
    ap.add_argument("--no-open", action="store_true")
    ap.add_argument("--limpa-dir", default=None, help="outdir of an earlier limpa run")
    ap.add_argument(
        "--chromatography-dir", default=None, help="output of scripts/raw_chromatography_modal.py"
    )
    ap.add_argument("--no-surface", action="store_true", help="skip accessible-surfaceome tiers")
    a = ap.parse_args()
    lim = None
    if a.limpa_dir:
        from deliverome_analysis.limpa import load_limpa_result

        lim = load_limpa_result(a.limpa_dir)
    path, cov = dia_qc_dashboard(
        a.matrix_dir,
        layout=a.layout,
        output=a.out,
        open_browser=not a.no_open,
        limpa_result=lim,
        chromatography_dir=a.chromatography_dir,
        surface=not a.no_surface,
    )
    print(path)
    print(cov["status"].value_counts().to_string())
