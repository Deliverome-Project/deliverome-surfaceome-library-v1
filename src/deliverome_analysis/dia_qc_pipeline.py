"""One command from a Spectronaut export to a published DIA QC page.

Chains the reusable steps so a new dataset needs no glue code:

1. **trim** a Spectronaut Normal Report, or an alphaDIA output folder, to the small tables QC
   and limpa read (:func:`deliverome_analysis.spectronaut_trim.trim_spectronaut_report` /
   :func:`deliverome_analysis.alphadia_trim.trim_alphadia_output`, with the contaminant
   watchlist when S3 is reachable). Skipped when given an already-trimmed directory.
2. **limpa** protein quantification (limpa-rs by default; :mod:`deliverome_analysis.limpa`).
3. **chromatography** from the Thermo .raw files on Modal (TIC, spray, LC pressure;
   ``scripts/raw_chromatography_modal.py``), when ``raw_prefix`` is given.
4. **dashboard** (interactive HTML, run record, run_diagnostics.tsv) and, optionally, the
   static per-check PNGs.
5. **publish** everything under a ``processed://`` experiment prefix, then rebuild the shared
   QC index (:mod:`deliverome_analysis.qcdb`) so the runs appear in history queries.

Each step reuses its output when it already exists in ``out_dir``, so a rerun only redoes what
is missing (delete a subfolder to force it).

    from deliverome_analysis.dia_qc_pipeline import run_dia_qc
    run_dia_qc("report.tsv", "qc/", raw_prefix="external://thermo/astral-zoom/20260916/",
               title="Astral HeLa", publish_to="sprint5/s5e7-dia-limpa-benchmark/qc/2026-10-08/")

or from the shell:

    uv run python -m deliverome_analysis.dia_qc_pipeline report.tsv --out qc/ \\
        --raw-prefix external://thermo/astral-zoom/20260916/ --title "Astral HeLa"
"""

from __future__ import annotations

import os
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
MODAL_SCRIPT = REPO_ROOT / "scripts" / "raw_chromatography_modal.py"


@dataclass
class QCRun:
    """Where each output of :func:`run_dia_qc` landed."""

    out_dir: Path
    trimmed: Path
    limpa: Path | None = None
    chromatography: Path | None = None
    dashboard: Path | None = None
    figures: Path | None = None
    published: list[str] = field(default_factory=list)


def _is_trimmed(path: Path) -> bool:
    return path.is_dir() and (path / "precursor_qc.parquet").exists()


def trim_step(source: str | os.PathLike, out_dir: Path) -> Path:
    """Trimmed directory for ``source`` (a report file, or a directory that is already trimmed)."""
    src = Path(source)
    if _is_trimmed(src):
        return src
    dest = out_dir / "trimmed"
    if _is_trimmed(dest):
        return dest
    from deliverome_analysis.alphadia_trim import is_alphadia_output, trim_alphadia_output
    from deliverome_analysis.spectronaut_trim import trim_spectronaut_report

    watchlist = None
    try:
        from deliverome_analysis.contaminants import load_watchlist

        watchlist = load_watchlist()
    except Exception as exc:  # no S3: trim without flagging watchlist proteins
        print(f"contaminant watchlist unavailable ({type(exc).__name__}); trimming without it")
    if is_alphadia_output(src):  # an alphaDIA output folder instead of a Spectronaut report
        trim_alphadia_output(src, dest, watchlist=watchlist)
    else:
        trim_spectronaut_report(src, dest, watchlist=watchlist)
    return dest


def limpa_step(trimmed: Path, out_dir: Path, *, engine: str = "rust", cores: int = 4) -> Path:
    """limpa protein quantification into ``out_dir/limpa`` (reused when complete)."""
    dest = out_dir / "limpa"
    if (dest / "protein_log2.tsv").exists():
        return dest
    from deliverome_analysis.limpa import run_limpa_spectronaut

    run_limpa_spectronaut(matrix_dir=trimmed, cores=cores, outdir=dest, engine=engine)
    return dest


def chromatography_step(trimmed: Path, out_dir: Path, raw_prefix: str) -> Path:
    """TIC, spray and LC traces from the raw files on Modal into ``out_dir/chromatography``.

    Raw files are found as ``<raw_prefix><R.FileName>.raw``; credentials stay local (Modal
    receives presigned URLs). Needs ``uv sync --group modal`` and a Modal login.
    """
    dest = out_dir / "chromatography"
    if (dest / "scans.parquet").exists():
        return dest
    if not MODAL_SCRIPT.exists():
        raise FileNotFoundError(f"{MODAL_SCRIPT} not found; run from a repository checkout")
    cmd = [
        "uv", "run", "--group", "modal", "modal", "run", str(MODAL_SCRIPT),
        "--matrix-dir", str(trimmed), "--raw-prefix", raw_prefix, "--out-dir", str(dest),
    ]  # fmt: skip
    subprocess.run(cmd, check=True, cwd=REPO_ROOT)
    return dest


def publish_step(run: QCRun, prefix: str) -> list[str]:
    """Upload every output under ``processed://<prefix>`` (an experiment-folder key)."""
    from deliverome_analysis.s3 import processed

    prefix = prefix.rstrip("/") + "/"
    keys = []
    for sub in (run.limpa, run.chromatography, run.figures):
        if sub is None:
            continue
        for f in sorted(p for p in sub.iterdir() if p.is_file()):
            keys.append(processed.upload(f, f"{prefix}{sub.name}/{f.name}"))
    if run.dashboard is not None:
        for f in sorted(run.dashboard.parent.glob(run.dashboard.stem + ".*")):
            keys.append(processed.upload(f, prefix + f.name))
    return keys


def run_dia_qc(
    source: str | os.PathLike,
    out_dir: str | os.PathLike,
    *,
    raw_prefix: str | None = None,
    chromatography_dir: str | os.PathLike | None = None,
    limpa: bool = True,
    engine: str = "rust",
    cores: int = 4,
    layout: str = "default",
    title: str | None = None,
    figures: bool = False,
    publish_to: str | None = None,
    open_browser: bool = False,
    history: bool = True,
) -> QCRun:
    """Run the whole QC (see the module docstring) and return where everything is.

    ``source``: a Spectronaut Normal Report or a ``trim_spectronaut_report`` directory.
    ``raw_prefix``: S3 prefix of the raw files (``external://...`` / ``raw://...``) to read the
    chromatography on Modal; or pass ``chromatography_dir`` if that was already extracted.
    ``publish_to``: a ``processed://`` key prefix (``sprint{N}/s{N}e{M}-slug/...``) to upload to.
    ``history``: plot earlier runs from the shared QC index (needs S3).
    """
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    run = QCRun(out_dir=out, trimmed=trim_step(source, out))
    if limpa:
        run.limpa = limpa_step(run.trimmed, out, engine=engine, cores=cores)
    if chromatography_dir is not None:
        run.chromatography = Path(chromatography_dir)
    elif raw_prefix:
        run.chromatography = chromatography_step(run.trimmed, out, raw_prefix)

    from deliverome_analysis.dia_qc_dashboard import dia_qc_dashboard

    lim = None
    if run.limpa is not None:
        from deliverome_analysis.limpa import load_limpa_result

        lim = load_limpa_result(run.limpa)
    run.dashboard, _ = dia_qc_dashboard(
        run.trimmed,
        layout=layout,
        output=out / "dia_qc.html",
        title=title,
        open_browser=open_browser,
        limpa_result=lim,
        chromatography_dir=run.chromatography,
        history=history,
    )
    if figures:
        from deliverome_analysis.dia_qc import dia_qc_report

        run.figures = out / "figures"
        dia_qc_report(matrix_dir=run.trimmed, output_dir=run.figures, limpa_result=lim)
    if publish_to:
        run.published = publish_step(run, publish_to)
        try:  # refresh the shared QC index so this run shows up in history queries
            from deliverome_analysis import qcdb

            qcdb.build("dia_qc")
        except Exception as exc:
            print(f"QC index not rebuilt ({type(exc).__name__}: {exc}); run: qcdb build dia_qc")
    return run


def main(argv: list[str] | None = None) -> None:
    import argparse

    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("source", help="Spectronaut Normal Report (.tsv) or a trimmed directory")
    ap.add_argument("--out", required=True, help="output directory")
    ap.add_argument("--raw-prefix", default=None, help="S3 prefix of the .raw files (Modal step)")
    ap.add_argument("--chromatography-dir", default=None, help="already-extracted traces")
    ap.add_argument("--no-limpa", action="store_true")
    ap.add_argument("--engine", default="rust", choices=["rust", "r"])
    ap.add_argument("--cores", type=int, default=4)
    ap.add_argument("--layout", default="default")
    ap.add_argument("--title", default=None)
    ap.add_argument("--figures", action="store_true", help="also write the static PNGs")
    ap.add_argument("--publish-to", default=None, help="processed:// key prefix to upload to")
    ap.add_argument("--open", action="store_true", help="open the page in the browser")
    a = ap.parse_args(argv)
    run = run_dia_qc(
        a.source,
        a.out,
        raw_prefix=a.raw_prefix,
        chromatography_dir=a.chromatography_dir,
        limpa=not a.no_limpa,
        engine=a.engine,
        cores=a.cores,
        layout=a.layout,
        title=a.title,
        figures=a.figures,
        publish_to=a.publish_to,
        open_browser=a.open,
    )
    print(f"dashboard: {run.dashboard}")
    for k in run.published:
        print(f"processed://{k}")


if __name__ == "__main__":
    main(sys.argv[1:])
