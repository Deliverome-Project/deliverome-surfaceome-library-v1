"""limpa protein quantification + differential expression for Spectronaut DIA reports.

`limpa <https://bioconductor.org/packages/limpa>`_ (Li, Cobbold & Smyth, 2025) is an R /
Bioconductor package from the limma authors. Instead of imputing missing precursor values it
fits a **detection probability curve (DPC)** -- the lower a precursor's true intensity, the
more likely it is to go undetected -- and uses it to estimate a complete protein x run matrix
*with a standard error for every value*. Those standard errors become precision weights in a
limma differential-expression fit.

This module shells out to ``Rscript`` (driver: ``r/limpa_spectronaut.R`` next to this file)
and returns pandas DataFrames, so callers stay in Python.

Setup (once per machine)
------------------------
Needs R (the version pinned in ``renv.lock``) plus the R packages it pins (limpa, limma, ...),
installed into the repo's renv project library::

    scripts/install-r-deps.sh

Check with :func:`limpa_available`. Set ``DELIVEROME_RSCRIPT`` to use a non-PATH ``Rscript``.

Input: what to export from Spectronaut
--------------------------------------
A **Normal Report** (long format: one row per precursor per run) -- the "BGS Factory Report
(Normal)" works as-is. Required columns: ``R.FileName``, ``EG.ModifiedSequence``,
``FG.Charge``, ``EG.TotalQuantity (Settings)``, ``EG.Qvalue``, ``PG.Qvalue`` and the protein
ID column (default ``PG.ProteinGroups``). Strongly recommended: ``EG.IsImputed``.

* Values with ``EG.Qvalue`` or ``PG.Qvalue`` above ``q_cutoff`` become missing.
* Values flagged ``EG.IsImputed == True`` are **dropped** -- limpa must see real missingness,
  not Spectronaut's imputed fill-ins. Without that column they can't be removed, and you get
  a warning; set Spectronaut's imputing strategy to "No Imputing" to be safe either way.
* Do not export a protein-level (PG) report: limpa needs precursor-level values.

Usage
-----
    from deliverome_analysis.s3 import processed
    from deliverome_analysis.limpa import run_limpa_spectronaut

    report = processed.download("sprintN/sNeM-slug/Report.tsv")

    # Quantification only (any number of runs >= 2):
    res = run_limpa_spectronaut(report)
    res.protein_log2      # protein x run, log2, no NAs
    res.protein_se        # matching standard errors
    res.proteins          # protein_id, PG.ProteinAccessions, NPrec, PropObs

    # Differential expression -- needs replicate runs per group:
    samples = pd.DataFrame({"run": [...], "group": ["ctrl", "ctrl", "ctrl", "trt", "trt", "trt"]})
    res = run_limpa_spectronaut(report, samples=samples, formula="~ group")
    res.de["grouptrt"]    # limma topTable: logFC, AveExpr, t, P.Value, adj.P.Val, B

Protein IDs are Spectronaut protein groups (``;``-joined UniProt accessions). Join to other
data on accessions, never gene symbols -- see ``deliverome_analysis.surfaceome``.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

R_SCRIPT = Path(__file__).parent / "r" / "limpa_spectronaut.R"
COMPARATORS = ("maxlfq_proda",)
# R packages are pinned by renv (renv.lock at the repo root). Running Rscript from the repo
# root lets the root .Rprofile activate that project library instead of the user library.
REPO_ROOT = Path(__file__).resolve().parents[2]
INSTALL_HINT = "Run scripts/install-r-deps.sh to restore the pinned R packages (renv.lock)."

# limpa-rs: Rust port of limpa 1.4.2's per-protein quantification (Deliverome-Project/limpa-rs),
# identical to R limpa on our Astral export to 7e-8 log2 and ~11x faster (65 s vs ~5 h at 384
# runs on an expanded workload). It pins exact pandas/numpy versions, so it never shares our
# environment: it runs from its checksum-verified prebuilt wheel in an isolated, throwaway uv
# environment, with its own pinned R library (``python -m deliverome_analysis.limpa setup-rs``).
LIMPA_RS_VERSION = "0.1.0rc1"
LIMPA_RS_RELEASE = "https://github.com/Deliverome-Project/limpa-rs/releases/download/v0.1.0-rc.1/"
LIMPA_RS_WHEELS = {  # (system, machine) -> (wheel, sha256 from the release's SHA256SUMS)
    ("Darwin", "arm64"): (
        "limpa_rs-0.1.0rc1-py3-none-macosx_11_0_arm64.whl",
        "c7ea3b4b14a60ff3ac197394ea42f87b0cedc4b7966fe6aafc339d1a04b66562",
    ),
    ("Linux", "x86_64"): (
        "limpa_rs-0.1.0rc1-py3-none-linux_x86_64.whl",
        "0b9b36100aaef01f16e8b0ef5e9daa9156955948babc6e9b0af6898ea4c6d909",
    ),
}
LIMPA_RS_SEED = 141  # limpa-rs's default DPC subsampling seed; recorded in the summary
ENGINES = ("rust", "r")


def _limpa_rs_wheel() -> Path:
    """Download (once) and verify the limpa-rs wheel for this platform; return its path."""
    import hashlib
    import platform
    import urllib.request

    key = (platform.system(), platform.machine())
    if key not in LIMPA_RS_WHEELS:
        raise RuntimeError(
            f"limpa-rs has no prebuilt wheel for {key}; use engine='r' (identical results)."
        )
    name, sha = LIMPA_RS_WHEELS[key]
    cache = Path(os.environ.get("XDG_CACHE_HOME", Path.home() / ".cache")) / "deliverome-analysis"
    path = cache / "limpa-rs" / name
    if not path.exists() or hashlib.sha256(path.read_bytes()).hexdigest() != sha:
        path.parent.mkdir(parents=True, exist_ok=True)
        tmp = path.with_suffix(".part")
        urllib.request.urlretrieve(LIMPA_RS_RELEASE + name, tmp)
        got = hashlib.sha256(tmp.read_bytes()).hexdigest()
        if got != sha:
            tmp.unlink()
            raise RuntimeError(f"limpa-rs wheel checksum mismatch: {got} != {sha}")
        tmp.replace(path)
    return path


def _limpa_rs_cmd(*python_args: str) -> list[str]:
    uv = shutil.which("uv")
    if not uv:
        raise RuntimeError("engine='rust' needs uv on PATH (it runs limpa-rs in its own env)")
    return [
        uv,
        "run",
        "--isolated",
        "--no-project",
        "--python",
        "3.14",
        "--with",
        str(_limpa_rs_wheel()),
        "python",
        *python_args,
    ]


def limpa_rs_setup() -> str:
    """One-time: restore limpa-rs's own pinned R environment into its user cache."""
    proc = subprocess.run(
        _limpa_rs_cmd("-m", "limpa_rs", "setup-r"), capture_output=True, text=True
    )
    if proc.returncode != 0:
        raise RuntimeError(f"limpa-rs setup failed:\n{proc.stderr.strip()[-2000:]}")
    return proc.stdout.strip()


def limpa_rs_available() -> bool:
    """True when the limpa-rs wheel fits this platform and its R environment is set up."""
    try:
        proc = subprocess.run(
            _limpa_rs_cmd("-m", "limpa_rs", "doctor"), capture_output=True, text=True, timeout=300
        )
    except RuntimeError, OSError, subprocess.TimeoutExpired:
        return False
    return proc.returncode == 0 and "not configured" not in (proc.stdout + proc.stderr)


_RS_DRIVER = """
import json, sys
import pandas as pd
from limpa_rs import run_limpa_spectronaut
a = json.loads(sys.argv[1])
if a.get("samples"):
    a["samples"] = pd.read_csv(a["samples"], sep="\\t")
run_limpa_spectronaut(**a)
"""


def _run_limpa_rs(src: Path, is_report: bool, out: Path, **kw) -> None:
    args = {
        ("report" if is_report else "matrix_dir"): str(src),
        "outdir": str(out),
        "engine": "rust",
        "seed": LIMPA_RS_SEED,
        **kw,
    }
    proc = subprocess.run(
        _limpa_rs_cmd("-c", _RS_DRIVER, json.dumps(args)), capture_output=True, text=True
    )
    if proc.returncode != 0:
        msg = proc.stderr.strip()
        if "not configured" in msg or "setup-r" in msg:
            msg += "\nRun once: uv run python -m deliverome_analysis.limpa setup-rs"
        raise RuntimeError(f"limpa-rs failed:\n{msg[-3000:]}")


@dataclass
class LimpaResult:
    """Outputs of one limpa run. All matrices are indexed by ``protein_id``."""

    protein_log2: pd.DataFrame
    protein_se: pd.DataFrame
    proteins: pd.DataFrame
    samples: pd.DataFrame
    summary: dict
    de: dict[str, pd.DataFrame] = field(default_factory=dict)
    dpc_points: pd.DataFrame | None = None  # precursor, mu, n_detected, nsamples (DPC fit set)
    maxlfq_log2: pd.DataFrame | None = None  # MaxLFQ protein x run with NAs (with ``compare``)
    # comparator DE tables: method -> {coefficient/contrast -> topTable-shaped DataFrame}
    compare: dict[str, dict[str, pd.DataFrame]] = field(default_factory=dict)
    outdir: Path | None = None

    @property
    def dpc(self) -> tuple[float, float]:
        """(beta0, beta1) of the fitted detection probability curve."""
        return self.summary["dpc_beta0"], self.summary["dpc_beta1"]


def detection_probability(log2_intensity, beta0: float, beta1: float):
    """P(precursor detected | true log2 intensity) under the fitted DPC (a logistic curve)."""
    x = np.asarray(log2_intensity, dtype=float)
    return 1.0 / (1.0 + np.exp(-(beta0 + beta1 * x)))


def _rscript() -> str | None:
    """Rscript to use: DELIVEROME_RSCRIPT if it is an executable, else Rscript on PATH."""
    env = os.environ.get("DELIVEROME_RSCRIPT")
    if env and shutil.which(env):
        return shutil.which(env)
    return shutil.which("Rscript")


def _r_cwd() -> Path | None:
    """Repo root if it carries the renv lockfile (editable install), else inherit the cwd."""
    return REPO_ROOT if (REPO_ROOT / "renv.lock").exists() else None


def limpa_available() -> bool:
    """True if Rscript exists and can load limpa."""
    rs = _rscript()
    if not rs:
        return False
    proc = subprocess.run(
        [rs, "-e", "quit(status = !requireNamespace('limpa', quietly = TRUE))"],
        capture_output=True,
        text=True,
        cwd=_r_cwd(),
    )
    return proc.returncode == 0


def run_limpa_spectronaut(
    report: str | os.PathLike | None = None,
    *,
    matrix_dir: str | os.PathLike | None = None,
    cores: int = 1,
    outdir: str | os.PathLike | None = None,
    protein_id: str = "PG.ProteinGroups",
    q_cutoff: float = 0.01,
    samples: pd.DataFrame | None = None,
    formula: str | None = None,
    contrasts: list[str] | None = None,
    compare: list[str] | None = None,
    verbose: bool = False,
    engine: str = "rust",
) -> LimpaResult:
    """Quantify proteins with limpa from a Spectronaut Normal Report, optionally test DE.

    Args:
        report: local path to the Spectronaut Normal Report TSV (download S3 keys first with
            ``processed.download(key)``). Give this or ``matrix_dir``.
        matrix_dir: output of :func:`deliverome_analysis.spectronaut_trim.trim_spectronaut_report`
            -- the compact precursor matrix; much faster and lighter for large cohorts.
        cores: worker processes for the per-protein fits (forked, so macOS/Linux only). Agrees
            with ``cores=1`` to optimizer tolerance (max 4e-4 log2 on the liver export), not
            bit-for-bit: limpa's starting values come from a whole-input imputation.
        outdir: where the R driver writes its TSVs. Default: a temp dir (outputs are still
            returned as DataFrames); pass a path to keep them.
        protein_id: report column used to group precursors into proteins.
        q_cutoff: EG.Qvalue / PG.Qvalue threshold; values above it are treated as missing.
        samples: optional run-level table with a ``run`` column matching ``R.FileName`` plus
            any covariates for ``formula``. Spectronaut's own ``R.Condition`` /
            ``R.Replicate`` are always available without it.
        formula: R model formula over the sample columns, e.g. ``"~ group"`` or
            ``"~ 0 + group"``. Omit for quantification only.
        contrasts: limma contrasts over the design's column names, e.g. ``["grouptrt-groupctrl"]``
            (use with a ``~ 0 + ...`` formula). Without it, one table per non-intercept
            coefficient.
        compare: comparator pipelines run on the same filtered precursors. Currently
            ``["maxlfq_proda"]`` (MaxLFQ via iq, then proDA), the runner-up kept for ongoing
            evaluation; the matrix is returned as ``maxlfq_log2`` and tables under
            ``compare["maxlfq_proda"]``. Needs ``formula`` for DE; without it only MaxLFQ is
            returned. (The full six-method comparison is S5E7: analysis/sprint5/s5e7-dia-limpa-benchmark/.)
        verbose: echo R's progress messages.
        engine: ``"rust"`` (default): limpa-rs, the Rust port of limpa's quantification,
            identical to R limpa to ~1e-7 log2 and an order of magnitude faster; run in its own
            isolated environment (one-time ``setup-rs``). ``"r"``: limpa itself via this repo's
            renv. Detection curve, hyperparameters and DE are R limpa/limma in both.

    Returns:
        LimpaResult with protein_log2, protein_se, proteins, samples, summary and de tables.
    """
    # Validate arguments before touching the filesystem or R, so bad calls fail the same way
    # on machines without R (e.g. CI).
    if (report is None) == (matrix_dir is None):
        raise ValueError("give exactly one of report or matrix_dir")
    if contrasts and not formula:
        raise ValueError("contrasts need a formula")
    if engine not in ENGINES:
        raise ValueError(f"engine must be one of {ENGINES}")
    unknown = set(compare or ()) - set(COMPARATORS)
    if unknown:
        raise ValueError(f"unknown compare method(s) {sorted(unknown)}; choose from {COMPARATORS}")
    # absolute paths: R runs from the repo root, not the caller's cwd
    src = Path(report if report is not None else matrix_dir).resolve()
    if not src.exists():
        raise FileNotFoundError(src)
    keep = outdir is not None
    if engine == "rust":
        # limpa-rs publishes its output directory atomically and refuses to overwrite one
        base = Path(outdir).resolve() if keep else Path(tempfile.mkdtemp(prefix="limpa-rs-"))
        out = base if keep else base / "out"
        if out.exists() and any(out.iterdir()):
            raise FileExistsError(f"{out} is not empty; limpa-rs will not overwrite results")
        if out.exists():
            out.rmdir()
        samples_path = None
        if samples is not None:
            samples_path = (base if not keep else out.parent) / f".{out.name}_samples.tsv"
            samples.to_csv(samples_path, sep="\t", index=False)
        try:
            _run_limpa_rs(
                src,
                report is not None,
                out,
                cores=int(cores),
                protein_id=protein_id,
                q_cutoff=q_cutoff,
                samples=str(samples_path) if samples_path else None,
                formula=formula,
                contrasts=contrasts,
                compare=compare,
                verbose=verbose,
            )
            res = _read_outputs(out, keep)
            res.summary["engine"] = f"limpa-rs {LIMPA_RS_VERSION} (rust)"
            return res
        finally:
            if samples_path is not None:
                samples_path.unlink(missing_ok=True)
            if not keep:
                shutil.rmtree(base, ignore_errors=True)

    rs = _rscript()
    if not rs:
        raise RuntimeError(f"Rscript not found on PATH. {INSTALL_HINT}")

    out = Path(outdir).resolve() if keep else Path(tempfile.mkdtemp(prefix="limpa-"))
    out.mkdir(parents=True, exist_ok=True)

    cmd = [
        rs,
        str(R_SCRIPT),
        f"--report={src}" if report is not None else f"--matrix-dir={src}",
        f"--cores={int(cores)}",
        f"--outdir={out}",
        f"--protein-id={protein_id}",
        f"--q-cutoff={q_cutoff}",
    ]
    if samples is not None:
        if "run" not in samples.columns:
            raise ValueError("samples needs a 'run' column matching Spectronaut R.FileName")
        samples_path = out / "samples_input.tsv"
        samples.to_csv(samples_path, sep="\t", index=False)
        cmd.append(f"--samples={samples_path}")
    if formula:
        cmd.append(f"--formula={formula}")
    if contrasts:
        cmd.append(f"--contrasts={','.join(contrasts)}")
    if compare:
        cmd.append(f"--compare={','.join(compare)}")

    proc = subprocess.run(cmd, capture_output=True, text=True, cwd=_r_cwd())
    if verbose and proc.stderr:
        print(proc.stderr)
    if proc.returncode != 0:
        msg = proc.stderr.strip()
        # Keep R's error itself, not the package-loading chatter before it.
        if "Error" in msg:
            msg = msg[msg.index("Error") :].replace("Execution halted", "").strip()
        if "there is no package called" in msg:
            msg += f"\n{INSTALL_HINT}"
        raise RuntimeError(f"limpa failed:\n{msg}")
    for line in proc.stderr.splitlines():
        if line.startswith("WARNING:"):
            print(line)

    try:
        res = _read_outputs(out, keep)
        res.summary["engine"] = "R limpa"
        return res
    finally:
        if not keep:
            shutil.rmtree(out, ignore_errors=True)


def load_limpa_result(outdir: str | os.PathLike) -> LimpaResult:
    """Re-load a LimpaResult from the ``outdir`` of an earlier ``run_limpa_spectronaut`` call."""
    return _read_outputs(Path(outdir).resolve(), True)


def _read_outputs(out: Path, keep: bool) -> LimpaResult:
    def tsv(name: str, index: str | None = "protein_id") -> pd.DataFrame:
        df = pd.read_csv(out / name, sep="\t", keep_default_na=True, na_values=["NA"])
        return df.set_index(index) if index else df

    summary = json.loads((out / "summary.json").read_text())
    de = {name: tsv(fname) for name, fname in (summary.get("de_tables") or {}).items()}
    compare: dict[str, dict[str, pd.DataFrame]] = {}
    for key, fname in (summary.get("compare_tables") or {}).items():
        method, coef = key.split("::", 1)
        compare.setdefault(method, {})[coef] = tsv(fname)
    maxlfq = out / "maxlfq_log2.tsv"
    return LimpaResult(
        protein_log2=tsv("protein_log2.tsv"),
        protein_se=tsv("protein_se.tsv"),
        proteins=tsv("protein_annotation.tsv"),
        samples=tsv("samples.tsv", index="run"),
        summary=summary,
        de=de,
        dpc_points=tsv("dpc_points.tsv", index=None) if (out / "dpc_points.tsv").exists() else None,
        maxlfq_log2=tsv("maxlfq_log2.tsv") if maxlfq.exists() else None,
        compare=compare,
        outdir=out if keep else None,
    )


if __name__ == "__main__":
    import sys

    if sys.argv[1:] == ["setup-rs"]:
        print(limpa_rs_setup())
    elif sys.argv[1:] == ["doctor-rs"]:
        print("limpa-rs ready" if limpa_rs_available() else "limpa-rs not ready: run setup-rs")
    else:
        sys.exit("usage: python -m deliverome_analysis.limpa setup-rs | doctor-rs")
