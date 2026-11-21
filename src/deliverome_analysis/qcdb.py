"""Shared, queryable SQLite indexes of QC and flow runs, rebuilt from the experiment folders.

The experiment folders on ``processed://`` stay the source of truth; nothing is edited in a
database. Each run already leaves its files there:

* **DIA QC** -- the dashboard writes ``<page>.qcrecord.json`` next to the page: per-run
  metrics (identifications, surface tiers, chromatography diagnosis, acquisition metadata),
  the layout YAML text and hash, and the code/package versions. A scientist's verdict on a
  flagged run goes in the sibling ``<page>.run_diagnostics.tsv`` (``confirmed_cause``,
  ``notes``) and is picked up on the next build.
* **Flow** -- every run published by :func:`deliverome_analysis.flow.publish_run` has a
  ``publication.json`` listing its ``gates.yaml``, ``run.json`` and ``summary.csv``.
* **Search** -- every DIA search run by deliverome-internal's ``workflows/alphadia_qc``
  runner leaves a ``benchmark.json`` (inputs and checksums, the computer it ran on, runtime,
  CPU time, peak and mean memory, software versions) and alphaDIA's ``alphadia/stat.tsv``
  (identifications), uploaded by ``scripts/publish_alphadia_runs.py``. One row per search, so
  runtime and resource use can be compared across machines, libraries and file sizes. Scans
  ``scratch/`` as well as the experiment folders, since QC batches often have no experiment.
  Each search is also checked against the archived raw file (``raw_provenance``: was the file
  the search read the same size and SHA-256 as the copy on ``raw://``? see
  :mod:`deliverome_analysis.raw_provenance`); the DIA QC index carries the same status per run.

``build`` writes one SQLite file per database to ``processed://db/<name>.sqlite``:

* **incremental** -- each run ("unit") is fingerprinted by the S3 ETags of the files it is read
  from, all taken from one listing. Units whose fingerprint is unchanged are kept from the
  previous snapshot; only new or changed ones are read, and deleted ones are dropped.
  ``full=True`` rebuilds from scratch (also used when the snapshot's schema is older).
* **resilient** -- each unit is ingested inside its own savepoint. One that fails (a malformed
  record, a summary that breaks the flow format contract) is rolled back, listed in
  ``build_errors`` with the reason, and retried on the next build; the rest still land.
* **safe under concurrency** -- the snapshot is replaced only if it is still the version the
  build started from (S3 ``If-Match``). If another build uploaded in between, this one starts
  over on top of the newer snapshot (a few retries), so no published run is dropped.

Read with the shared helpers (each returns a DataFrame from the current snapshot, ETag-cached):

    from deliverome_analysis import qcdb
    qcdb.qc_metrics(["precursors", "tic_area"], instrument="OA10120", wide=True)
    qcdb.qc_runs(since="2026-01-01")
    qcdb.qc_diagnoses()                               # flagged runs and recorded verdicts
    qcdb.flow_runs(experiment="s3e2-test-titer-flow")
    qcdb.flow_populations(experiment="s3e2-test-titer-flow", gate="singlets")  # + metadata
    qcdb.flow_stats(gate="singlets", stat="median", wide=True)
    qcdb.search_runs(host="DESKTOP-6IHF2SA")          # runtime, CPU, memory per search
    qcdb.search_with_qc(["precursors"], library=...)  # each search + its QC page's metrics
    qcdb.read("flow", "select ...", params)           # anything else, read-only
    con = qcdb.connect("dia_qc")                      # or a raw read-only sqlite3 connection

    uv run python -m deliverome_analysis.qcdb build all          # exits 1 if any unit failed
    uv run python -m deliverome_analysis.qcdb query flow "select * from build_errors"

Event-level flow reanalysis (re-gating, new statistics) runs on the events themselves through
DuckDB: :mod:`deliverome_analysis.flowdb`.
"""

from __future__ import annotations

import csv
import hashlib
import io
import json
import math
import os
import re
import shutil
import sqlite3
import sys
import tempfile
import time
from pathlib import Path

import pandas as pd

SCHEMA_VERSION = 9  # of the SQLite files; a different value forces a full rebuild
RECORD_SCHEMA_VERSION = 1  # of <page>.qcrecord.json, written by dia_qc_dashboard.qc_record
DATABASES = ("dia_qc", "flow", "search")
# ``build all`` order: search first, so the DIA QC index reads this build's raw provenance.
BUILD_ORDER = ("search", "flow", "dia_qc")
QC_RECORD_SUFFIX = ".qcrecord.json"
DIAGNOSIS_SUFFIX = ".run_diagnostics.tsv"
FLOW_MANIFEST = "publication.json"
SCAN_PREFIX = "sprint"  # experiment folders: sprint{N}/s{N}e{M}-slug/...
SCRATCH_PREFIX = "scratch/"
# Where each database looks. DIA QC pages and searches also live in scratch/ (QC batches with
# no experiment folder yet); for the QC history an experiment-folder record of an acquisition
# always outranks a scratch one (see latest_runs). Flow stays limited to experiment folders.
SCAN_PREFIXES = {
    "dia_qc": (SCAN_PREFIX, SCRATCH_PREFIX),
    "flow": (SCAN_PREFIX,),
    "search": (SCAN_PREFIX, SCRATCH_PREFIX),
}
SEARCH_RECORD = "benchmark.json"
SEARCH_STATS = "alphadia/stat.tsv"
SEARCH_CONFIG = "config.yaml"
UPLOAD_ATTEMPTS = 3
# Replaced searches moved aside by scripts/publish_alphadia_runs.py
# (<method>/_superseded/<run>/<old start>/): kept for the record, never indexed.
SUPERSEDED_DIR = "_superseded"

# The agentflow output format this module reads. agentflow's tests/test_output_contract.py pins
# the same shape on the producing side; change both together.
FLOW_SUMMARY_COLUMNS = frozenset(
    {"sample_id", "gate", "parent", "count", "parent_count", "percent_parent", "percent_total"}
)
FLOW_SUMMARY_OPTIONAL = frozenset({"reviewed", "experiment_label", "is_example"})
FLOW_STAT_COLUMN = re.compile(
    r"^(finite_signal_count|(?:median|mean|min|max|p05|p25|p75|p95|sd)_signal):(?P<det>.+)$"
)
FLOW_RUN_KEYS = frozenset({"inputs", "agentflow"})


def db_key(name: str) -> str:
    if name not in DATABASES:
        raise ValueError(f"unknown database {name!r}; one of {DATABASES}")
    return f"db/{name}.sqlite"


def _bucket():
    from deliverome_analysis.s3 import processed

    return processed


def _raw_bucket():
    from deliverome_analysis.s3 import raw

    return raw


def _read_text(bucket, key: str) -> str:
    return Path(bucket.download(key)).read_text()


def _sha256(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


def _now() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


# ----------------------------------------------------------------------------------------------
# DIA QC
# ----------------------------------------------------------------------------------------------

_QC_SCHEMA = """
create table records (
  record_id text primary key, record_key text, folder text, created text, source text,
  title text, layout text, layout_sha256 text, layout_yaml text, git_commit text,
  git_dirty text, uv_lock_sha256 text, renv_lock_sha256 text, limpa_engine text,
  limpa_version text, versions_json text, schema_version integer,
  scope text,  -- 'experiment' (sprint*/ folder) or 'scratch'
  -- peptide_metrics.DEFINITIONS_VERSION and DigestRules the record's numbers were computed
  -- under (NULL: a record from before versioning, definitions unknown)
  definitions_version text, digest_rules_json text);

-- run_uid identifies an acquisition across datasets: instrument serial + acquisition time from
-- the raw-file header, or the run name alone when the record has no raw-file metadata.
create table runs (
  record_id text, run text, run_uid text, label text, condition text, acquired text,
  instrument_model text, instrument_name text, instrument_serial text, file_bytes integer,
  search_library text, search_library_sha256 text,
  -- the raw file the run's search read (alphaDIA trims), and how it compares with the copy
  -- on raw:// (deliverome_analysis.raw_provenance). The raw_* columns are re-derived from
  -- search.sqlite on every build, so they follow the search index.
  search_raw_sha256 text, search_raw_bytes integer,
  raw_provenance text, raw_s3_key text, raw_s3_bytes integer, raw_s3_sha256 text,
  primary key (record_id, run));
-- definitions_version: the peptide_metrics version for metrics computed by that module, ''
-- for a metric it does not define, NULL when the record predates versioning (unknown).
create table run_metrics (
  record_id text, run text, metric text, value real, definitions_version text);
create table diagnoses (
  record_id text, run text, flags text, hints text, confirmed_cause text, notes text);
create index run_metrics_metric on run_metrics(metric);
create index runs_uid on runs(run_uid);
-- The newest record per acquisition, experiment folders before scratch: the views to use for
-- history plots and queries.
create view latest_runs as
  select * from (
    select r.*, rec.created, rec.scope, rec.definitions_version,
           row_number() over (partition by r.run_uid
                              order by rec.scope = 'experiment' desc, rec.created desc) as rank
    from runs r join records rec using (record_id))
  where rank = 1;
create view latest_run_metrics as
  select m.*, lr.run_uid, lr.acquired, lr.instrument_serial, lr.condition,
         lr.definitions_version as record_definitions_version
  from run_metrics m join latest_runs lr using (record_id, run);
"""


def run_uid(run: dict) -> str:
    serial, acquired = run.get("instrument_serial"), run.get("acquired")
    return f"{serial}|{acquired}" if serial and acquired else f"name:{run['run']}"


def _diagnosis_key(record_key: str) -> str:
    return record_key[: -len(QC_RECORD_SUFFIX)] + DIAGNOSIS_SUFFIX


def _qc_units(listing: dict[str, str]) -> dict[str, str]:
    """unit key -> fingerprint (record ETag + the editable diagnosis TSV's ETag)."""
    return {
        k: f"{etag}|{listing.get(_diagnosis_key(k), '')}"
        for k, etag in listing.items()
        if k.endswith(QC_RECORD_SUFFIX)
    }


def _diagnosis_overrides(bucket, record_key: str, listing: dict[str, str]) -> dict[str, dict]:
    """confirmed_cause / notes from the editable TSV next to a record, by run."""
    key = _diagnosis_key(record_key)
    if key not in listing:
        return {}
    out = {}
    for row in csv.DictReader(io.StringIO(_read_text(bucket, key)), delimiter="\t"):
        out[row.get("run", "")] = {
            "confirmed_cause": (row.get("confirmed_cause") or "").strip(),
            "notes": (row.get("notes") or "").strip(),
        }
    return out


def _qc_delete(con: sqlite3.Connection, unit: str) -> None:
    ids = [r[0] for r in con.execute("select record_id from records where record_key = ?", [unit])]
    for t in ("runs", "run_metrics", "diagnoses", "records"):
        con.executemany(f"delete from {t} where record_id = ?", [(i,) for i in ids])


def _qc_ingest(con: sqlite3.Connection, bucket, key: str, listing: dict[str, str]) -> bool:
    rec = json.loads(_read_text(bucket, key))
    if rec.get("kind") != "dia_qc":
        raise ValueError(f"not a dia_qc record (kind={rec.get('kind')!r})")
    if rec.get("schema_version", 0) > RECORD_SCHEMA_VERSION:
        raise ValueError(f"record schema {rec['schema_version']} is newer than this reader")
    rid = rec["record_id"]
    v, prov, lim = rec.get("versions", {}), rec.get("provenance", {}), rec.get("limpa") or {}
    defs = rec.get("peptide_metrics")  # absent in records from before versioning
    con.execute(
        "insert into records values (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
        (
            rid, key, key.rsplit("/", 1)[0], rec.get("created"), rec.get("source"),
            rec.get("title"), rec.get("layout"), rec.get("layout_sha256"), rec.get("layout_yaml"),
            prov.get("git-commit"), prov.get("git-dirty"), prov.get("uv-lock-sha256"),
            prov.get("renv-lock-sha256"), lim.get("engine"), lim.get("limpa_version"),
            json.dumps(v), rec.get("schema_version"),
            "scratch" if key.startswith(SCRATCH_PREFIX) else "experiment",
            defs["definitions_version"] if defs else None,
            json.dumps(defs["digest_rules"], sort_keys=True) if defs else None,
        ),
    )  # fmt: skip
    for r in rec.get("runs", []):
        con.execute(
            "insert into runs values (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
            (
                rid, r["run"], run_uid(r), r.get("label"), r.get("condition"), r.get("acquired"),
                r.get("instrument_model"), r.get("instrument_name"), r.get("instrument_serial"),
                r.get("file_bytes"), r.get("search_library"), r.get("search_library_sha256"),
                r.get("search_raw_sha256"), r.get("search_raw_bytes"), None, None, None, None,
            ),
        )  # fmt: skip
    metrics = rec.get("metrics", [])
    if any(len(m) != 3 or not isinstance(m[2], int | float) for m in metrics):
        raise ValueError("metrics must be [run, metric, number] triples")
    versioned = set(defs.get("metrics", [])) if defs else set()

    def version_of(metric: str) -> str | None:
        if not defs:
            return None
        return defs["definitions_version"] if metric in versioned else ""

    con.executemany(
        "insert into run_metrics values (?,?,?,?,?)",
        [(rid, r, m, v, version_of(m)) for r, m, v in metrics],
    )
    verdicts = _diagnosis_overrides(bucket, key, listing)
    for d in rec.get("diagnosis", []):
        o = verdicts.get(d["run"], {})
        con.execute(
            "insert into diagnoses values (?,?,?,?,?,?)",
            (rid, d["run"], d.get("flags"), d.get("hints"), o.get("confirmed_cause", ""),
             o.get("notes", "")),
        )  # fmt: skip
    return True


def _search_raw_index(bucket) -> list:
    """Raw-file facts of every indexed search (from the current search.sqlite), or []."""
    from deliverome_analysis.raw_provenance import SearchRaw

    if not bucket.etag(db_key("search")):
        return []
    con = sqlite3.connect(f"file:{bucket.download(db_key('search'))}?mode=ro", uri=True)
    try:
        cols = {r[1] for r in con.execute("pragma table_info(search_runs)")}
        if "raw_provenance" not in cols:  # a snapshot from before raw provenance existed
            return []
        rows = con.execute(
            "select run_uid, raw_sha256, raw_s3_key, raw_s3_bytes, raw_s3_sha256, raw_provenance "
            "from search_runs"
        ).fetchall()
    finally:
        con.close()
    return [SearchRaw(*r) for r in rows]


def qc_raw_provenance(searched_sha256, searched_bytes, uid, searches) -> tuple:
    """(status, raw key, raw bytes, raw SHA-256) of one QC run against its archived raw file.

    ``searched_sha256`` / ``searched_bytes`` describe the file the run's search read (an
    alphaDIA trim's ``run_info``); the archived copy comes from the search index entries of
    the same file or acquisition (``searches``, see :func:`_search_raw_index`). Without the
    search's own checksum and size the status is ``unknown`` (unless the raw file is
    missing): a Spectronaut page or an older QC record cannot say which file it saw.
    """
    from deliverome_analysis.raw_provenance import archived_object, classify

    obj = archived_object(searched_sha256, uid, searches)
    status = classify(searched_sha256, searched_bytes, obj)
    if obj is None:
        return status, None, None, None
    return status, obj.key, obj.size, obj.sha256


def raw_provenance_of_runs(runs: dict[str, dict], *, bucket=None) -> dict[str, str]:
    """{run: raw_provenance status} for a QC page's runs, from the current search index.

    ``runs`` maps run -> an alphaDIA trim's ``run_info`` entry (``search_raw_sha256``,
    ``search_raw_bytes``, ``run_uid``). Used by the dashboard to flag a run whose search read
    a different file than the one archived on raw://. Needs S3 (reads search.sqlite).
    """
    searches = _search_raw_index(bucket or _bucket())
    return {
        run: qc_raw_provenance(
            info.get("search_raw_sha256"), info.get("search_raw_bytes"), info.get("run_uid"),
            searches,
        )[0]
        for run, info in runs.items()
    }  # fmt: skip


def _qc_finalize(con: sqlite3.Connection, bucket, raw_bucket) -> dict:
    """Re-derive every QC run's raw provenance from the current search index (cheap, local)."""
    searches = _search_raw_index(bucket)
    rows = con.execute(
        "select record_id, run, run_uid, search_raw_sha256, search_raw_bytes from runs"
    ).fetchall()
    updates = [
        (*qc_raw_provenance(sha, size, uid, searches), rid, run)
        for rid, run, uid, sha, size in rows
    ]
    con.execute("begin")
    con.executemany(
        "update runs set raw_provenance = ?, raw_s3_key = ?, raw_s3_bytes = ?, raw_s3_sha256 = ? "
        "where record_id = ? and run = ?",
        updates,
    )
    con.execute("commit")
    return _provenance_counts(con, "runs")


def _provenance_counts(con: sqlite3.Connection, table: str) -> dict:
    return dict(
        con.execute(
            f"select coalesce(raw_provenance, 'unknown'), count(*) from {table} group by 1"
        ).fetchall()
    )


# ----------------------------------------------------------------------------------------------
# flow (agentflow)
# ----------------------------------------------------------------------------------------------

_FLOW_SCHEMA = """
create table flow_runs (
  run_id text primary key, unit_key text, experiment text, prefix text, agentflow_version text,
  agentflow_commit text, recipe_sha256 text, gates_sha256 text, gates_yaml text,
  input_count integer, signal_space text, layout_version integer, events_key text,
  versions_json text, producer text);
create table flow_inputs (run_id text, sample_id text, sha256 text, signal_space text);
create table populations (
  run_id text, sample_id text, gate text, parent text, count integer, parent_count integer,
  percent_parent real, percent_total real, reviewed integer, experiment_label text);
create table sample_metadata (run_id text, sample_id text, key text, value text);
create table signal_stats (
  run_id text, sample_id text, gate text, detector text, stat text, value real);
create index signal_stats_lookup on signal_stats(run_id, gate, detector, stat);
"""

_POP_COLUMNS = (
    "sample_id", "gate", "parent", "count", "parent_count", "percent_parent", "percent_total",
    "reviewed", "experiment_label",
)  # fmt: skip


def _manifest_file(manifest: dict, name: str) -> str | None:
    """Key of the published file whose path ends with ``name`` (layouts differ by version)."""
    hits = [f["key"] for f in manifest.get("files", []) if f["key"].endswith("/" + name)]
    return sorted(hits, key=len)[0] if hits else None


def check_flow_summary(summary: pd.DataFrame) -> None:
    """Raise ValueError if a summary.csv breaks the agentflow format contract this reader uses."""
    missing = FLOW_SUMMARY_COLUMNS - set(summary.columns)
    if missing:
        raise ValueError(f"summary.csv lacks required columns {sorted(missing)}")
    known = FLOW_SUMMARY_COLUMNS | FLOW_SUMMARY_OPTIONAL
    odd = [
        c
        for c in summary.columns
        if c not in known and not c.startswith("metadata:") and not FLOW_STAT_COLUMN.match(c)
    ]
    if odd:
        raise ValueError(f"summary.csv has columns outside the format contract: {odd[:5]}")
    bad = [
        c
        for c in summary.columns
        if FLOW_STAT_COLUMN.match(c) and not pd.api.types.is_numeric_dtype(summary[c])
    ]
    if bad:
        raise ValueError(f"summary.csv statistic columns are not numeric: {bad[:5]}")


def _summary_long(run_id: str, summary: pd.DataFrame):
    """populations, long signal_stats and sample_metadata rows from an agentflow summary.csv."""
    pops = []
    for _, r in summary.iterrows():
        row = [run_id] + [r.get(c) for c in _POP_COLUMNS]
        row[8] = int(bool(row[8])) if pd.notna(row[8]) else None  # reviewed
        pops.append(tuple(None if (isinstance(x, float) and pd.isna(x)) else x for x in row))
    meta_cols = [c for c in summary.columns if c.startswith("metadata:")]
    meta = (
        summary.drop_duplicates("sample_id")
        .melt(id_vars=["sample_id"], value_vars=meta_cols, var_name="key")
        .dropna(subset=["value"])
    )
    metadata = [
        (run_id, s, k.split(":", 1)[1], str(v))
        for s, k, v in zip(meta["sample_id"], meta["key"], meta["value"], strict=True)
    ]
    stat_cols = [c for c in summary.columns if FLOW_STAT_COLUMN.match(c)]
    long = summary.melt(id_vars=["sample_id", "gate"], value_vars=stat_cols, var_name="col")
    long = long.dropna(subset=["value"])
    stat_det = long["col"].str.split(":", n=1, expand=True)
    long["stat"] = stat_det[0].str.removesuffix("_signal")
    long["detector"] = stat_det[1]
    stats = list(
        zip(
            [run_id] * len(long), long["sample_id"], long["gate"], long["detector"], long["stat"],
            long["value"].astype(float), strict=True,
        )
    )  # fmt: skip
    return pops, stats, metadata


def _flow_units(listing: dict[str, str]) -> dict[str, str]:
    """unit key (publication.json) -> fingerprint (ETags of every file in its folder)."""
    units = {}
    for k in listing:
        if k.endswith("/" + FLOW_MANIFEST):
            folder = k.rsplit("/", 1)[0] + "/"
            tags = sorted(f"{f}={e}" for f, e in listing.items() if f.startswith(folder))
            units[k] = _sha256("|".join(tags))
    return units


def _flow_delete(con: sqlite3.Connection, unit: str) -> None:
    ids = [r[0] for r in con.execute("select run_id from flow_runs where unit_key = ?", [unit])]
    for t in ("flow_inputs", "populations", "sample_metadata", "signal_stats", "flow_runs"):
        con.executemany(f"delete from {t} where run_id = ?", [(i,) for i in ids])


def _flow_ingest(con: sqlite3.Connection, bucket, key: str, listing: dict[str, str]) -> bool:
    manifest = json.loads(_read_text(bucket, key))
    if manifest.get("status") != "complete":
        return False  # an unfinished publication is not a run yet
    prefix = manifest.get("prefix") or key.rsplit("/", 1)[0]
    run_json = _manifest_file(manifest, "run.json")
    summary_key = _manifest_file(manifest, "summary.csv")
    if not (run_json and summary_key):
        raise ValueError("publication lists no run.json or summary.csv")
    prov = json.loads(_read_text(bucket, run_json))
    missing = FLOW_RUN_KEYS - set(prov)
    if missing:
        raise ValueError(f"run.json lacks {sorted(missing)}")
    summary = pd.read_csv(bucket.download(summary_key))
    check_flow_summary(summary)
    gates_key = _manifest_file(manifest, "gates.yaml") or _manifest_file(
        manifest, "reproducibility.yaml"
    )
    gates = _read_text(bucket, gates_key) if gates_key else ""
    af = prov.get("agentflow") or manifest.get("agentflow") or {}
    parts = prefix.split("/")
    con.execute(
        "insert into flow_runs values (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
        (
            prefix, key, parts[1] if len(parts) > 1 else "", prefix, af.get("version"),
            af.get("git_commit"), prov.get("recipe_sha256"), _sha256(gates) if gates else None,
            gates, len(prov.get("inputs", [])), prov.get("signal_space"),
            manifest.get("layout_version"), _manifest_file(manifest, "events.parquet"),
            json.dumps(prov.get("versions", {})),
            # who drew the gates: agentflow, flowjo (applied by flow_import), or script
            prov.get("producer", "agentflow"),
        ),
    )  # fmt: skip
    con.executemany(
        "insert into flow_inputs values (?,?,?,?)",
        [
            (prefix, i.get("sample_id"), i.get("sha256"), i.get("signal_space"))
            for i in prov["inputs"]
        ],
    )
    pops, stats, metadata = _summary_long(prefix, summary)
    con.executemany(f"insert into populations values ({','.join('?' * 10)})", pops)
    con.executemany("insert into sample_metadata values (?,?,?,?)", metadata)
    con.executemany("insert into signal_stats values (?,?,?,?,?,?)", stats)
    return True


# ----------------------------------------------------------------------------------------------
# build / connect
# ----------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------
# search runs (alphaDIA runtime and resources)
# ----------------------------------------------------------------------------------------------

_SEARCH_SCHEMA = """
create table search_runs (
  run_key text primary key, folder text, run text, run_uid text, acquired text,
  instrument_serial text, instrument_model text, status text, started_utc text,
  engine text, engine_version text, raw_file text, raw_bytes integer, raw_sha256 text,
  raw_storage text, n_raw integer, library text, library_sha256 text, fasta_sha256 text,
  config text, config_sha256 text, threads integer,
  host text, cpu_model text, physical_cpus integer, logical_cpus integer, ram_gib real,
  os text, wsl integer,
  wall_seconds real, cpu_user_seconds real, cpu_system_seconds real, mean_cores_busy real,
  peak_pss_gib real, mean_pss_gib real, peak_system_used_gib real, mean_system_used_gib real,
  max_process_rss_gib real, stage_seconds real,
  precursors integer, protein_groups integer, gradient_seconds real,
  ms2_mass_bias_ppm real, ms2_mass_variance real, ms1_mass_bias_ppm real, ms1_mass_variance real,
  ms2_tolerance_ppm real, ms1_tolerance_ppm real,
  repo_commit text, repo_dirty integer, uv_lock_sha256 text, versions_json text,
  method_file text, scans integer, run_minutes real, cycle_time_s real,
  faims_attached integer, faims_cvs text,
  ms1_analyzer text, ms1_resolution real, ms1_mz_low real, ms1_mz_high real, ms1_max_it_ms real,
  ms2_analyzer text, dia_windows integer, dia_mz_low real, dia_mz_high real,
  dia_width_min real, dia_width_median real, dia_width_max real, dia_scheme_sha256 text,
  ms2_nce text, ms2_max_it_ms text, spray_voltage_v real, method_json text,
  -- the archived copy of the searched raw file on raw:// and how it compares
  -- (deliverome_analysis.raw_provenance); filled after ingest, see _search_finalize
  raw_s3_key text, raw_s3_bytes integer, raw_s3_sha256 text, raw_provenance text);
-- raw:// objects already looked up; immutable, so a known object is never HEADed again
create table raw_objects (
  raw_key text primary key, bytes integer, sha256 text, etag text, checked text);
create index search_runs_host on search_runs(host);
create index search_runs_uid on search_runs(run_uid);
create view acquisitions as  -- one row per acquisition with its method (any search of it)
select run_uid, min(run) as run, acquired, instrument_serial, instrument_model, method_file,
       scans, run_minutes, cycle_time_s, faims_attached, faims_cvs, ms1_analyzer,
       ms1_resolution, ms1_mz_low, ms1_mz_high, ms1_max_it_ms, ms2_analyzer, dia_windows,
       dia_mz_low, dia_mz_high, dia_width_min, dia_width_median, dia_width_max,
       dia_scheme_sha256, ms2_nce, ms2_max_it_ms, spray_voltage_v, method_json
from search_runs where method_json is not null and run_uid is not null group by run_uid;
"""
_SEARCH_COLUMNS = 79


def _method_columns(m: dict | None) -> tuple:
    """search_runs method columns from a benchmark input's ``method`` (raw_method.py)."""
    if not m or "ms2" not in m:
        return (
            (m or {}).get("sample", {}).get("method_file"),
            *([None] * 21),
            (json.dumps(m) if m else None),
        )
    ms1 = (m.get("ms1") or [{}])[0]
    ms2, faims = m["ms2"], m.get("faims") or {}
    iso = ms2.get("isolation_mz") or [None, None]
    rng = ms1.get("range") or [None, None]
    method_file = (m.get("sample") or {}).get("method_file")
    join = lambda xs: ";".join(f"{x:g}" for x in xs) if xs else None  # noqa: E731
    return (
        method_file.replace("\\", "/").split("/")[-1] if method_file else None,
        m.get("scans"), m.get("rt_end_min"), m.get("cycle_time_s"),
        None if faims.get("attached") is None else int(faims["attached"]),
        join(faims.get("cvs")),
        ms1.get("analyzer"), ms1.get("resolution"), rng[0], rng[1], ms1.get("max_it_ms"),
        ";".join(ms2.get("analyzers") or []) or None, ms2.get("windows"), iso[0], iso[1],
        ms2.get("width_min"), ms2.get("width_median"), ms2.get("width_max"),
        ms2.get("scheme_sha256"), join(ms2.get("nce")), join(ms2.get("max_it_ms")),
        (m.get("status_mid_run") or {}).get("spray_voltage_v"), json.dumps(m),
    )  # fmt: skip


def _search_units(listing: dict[str, str]) -> dict[str, str]:
    """unit key (the benchmark.json) -> fingerprint (its ETag + stat.tsv and config.yaml ETags)."""
    out = {}
    for k, etag in listing.items():
        if k.endswith("/" + SEARCH_RECORD):
            folder = k[: -len(SEARCH_RECORD)]
            out[k] = (
                f"{etag}|{listing.get(folder + SEARCH_STATS, '')}"
                f"|{listing.get(folder + SEARCH_CONFIG, '')}"
            )
    return out


def _search_delete(con: sqlite3.Connection, unit: str) -> None:
    con.execute("delete from search_runs where run_key = ?", [unit])


def _search_finalize(con: sqlite3.Connection, bucket, raw_bucket) -> dict:
    """Check each search's raw input against raw:// (new and still-undecided rows only).

    A row whose status is final (match / unverified / mismatch / size_mismatch) is never
    looked up again, and a raw object seen before comes from ``raw_objects``, so a rebuild
    only HEADs raw files of new searches. ``unknown`` and ``raw_missing`` rows are retried
    each build (a late raw upload, a credential that gained access).
    """
    from deliverome_analysis.raw_provenance import FINAL_STATUSES, RawLookup, RawObject, classify

    cache = [
        RawObject(key=k, exists=True, size=b, sha256=s, etag=e)
        for k, b, s, e in con.execute("select raw_key, bytes, sha256, etag from raw_objects")
    ]
    lookup = RawLookup(raw_bucket, cache)
    final = ",".join(f"'{s}'" for s in sorted(FINAL_STATUSES))
    todo = con.execute(
        "select run_key, folder, raw_file, raw_sha256, raw_bytes, n_raw from search_runs "
        f"where raw_provenance is null or raw_provenance not in ({final}) order by run_key"
    ).fetchall()
    con.execute("begin")
    for run_key, folder, raw_file, sha, size, n_raw in todo:
        key, obj = lookup.lookup(folder, raw_file) if raw_file else (None, None)
        # raw_bytes sums all inputs: it describes this raw file only for one-file searches
        status = classify(sha, size if n_raw == 1 else None, obj)
        con.execute(
            "update search_runs set raw_s3_key = ?, raw_s3_bytes = ?, raw_s3_sha256 = ?, "
            "raw_provenance = ? where run_key = ?",
            (key, obj.size if obj else None, obj.sha256 if obj else None, status, run_key),
        )
    now = _now()
    con.executemany(
        "insert or ignore into raw_objects values (?,?,?,?,?)",
        [(o.key, o.size, o.sha256, o.etag, now) for o in lookup.objects.values()],
    )
    con.execute("commit")
    if lookup.denied:
        print(f"search: raw:// not readable, provenance left unknown ({lookup.denied})")
    return _provenance_counts(con, "search_runs") | {"raw_lookups": lookup.calls}


def _search_ingest(con: sqlite3.Connection, bucket, key: str, listing: dict[str, str]) -> bool:
    rec = json.loads(_read_text(bucket, key))
    if "timing" not in rec or "inputs" not in rec:
        raise ValueError("not a search benchmark record (no timing/inputs)")
    folder = key[: -len(SEARCH_RECORD)]
    stats = {}
    if folder + SEARCH_STATS in listing:
        row = next(
            csv.DictReader(io.StringIO(_read_text(bucket, folder + SEARCH_STATS)), delimiter="\t")
        )

        def num(col):
            v = row.get(col)
            try:
                f = float(v)
            except TypeError, ValueError:
                return None
            return None if math.isnan(f) else f

        stats = {
            "precursors": int(float(row["search.precursors"])),
            "protein_groups": int(float(row["search.proteins"])),
            "gradient_seconds": float(row["raw.gradient_length"]),
            # mass accuracy: alphaDIA's calibration bias (ppm) and variance; MS1 stays None when
            # the search did not calibrate on MS1. Tolerances are the optimised search windows.
            "mass": (
                num("calibration.ms2_bias"), num("calibration.ms2_variance"),
                num("calibration.ms1_bias"), num("calibration.ms1_variance"),
                num("optimization.ms2_error"), num("optimization.ms1_error"),
            ),
        }  # fmt: skip
    config_sha = None
    if folder + SEARCH_CONFIG in listing:  # the exact config the search ran with
        config_sha = hashlib.sha256(
            Path(bucket.download(folder + SEARCH_CONFIG)).read_bytes()
        ).hexdigest()
    dirty = rec.get("software", {}).get("repo_dirty")
    dirty = None if dirty is None else int(bool(dirty))
    t, h, sw = rec["timing"], rec.get("host", {}), rec.get("software", {})
    pk = sw.get("packages", {})
    inp = rec["inputs"]
    lib, fasta = rec.get("library") or {}, rec.get("fasta") or {}
    acq = inp[0].get("acquisition") or {}
    con.execute(
        f"insert into search_runs values ({','.join('?' * _SEARCH_COLUMNS)})",
        (
            key, folder.rstrip("/"), rec.get("run_id"),
            inp[0].get("run_uid") or run_uid({"run": rec.get("run_id"), **acq}),
            acq.get("acquired"), acq.get("instrument_serial"), acq.get("instrument_model"),
            rec.get("status"), rec.get("started_utc"),
            "alphadia" if "alphadia" in pk else None, pk.get("alphadia"),
            Path(inp[0]["path"]).name, sum(i.get("bytes") or 0 for i in inp),
            inp[0].get("sha256"), rec.get("raw_storage"), len(inp),
            Path(lib["path"]).name if lib.get("path") else "predicted", lib.get("sha256"),
            fasta.get("sha256"), Path(rec.get("config_source") or "").name or None,
            config_sha, rec.get("threads"),
            h.get("hostname"), h.get("cpu_model"), h.get("physical_cpus"), h.get("logical_cpus"),
            h.get("ram_gib"), h.get("os"), int(bool(h.get("wsl"))),
            t.get("wall_seconds"), t.get("cpu_user_seconds"), t.get("cpu_system_seconds"),
            t.get("mean_cores_busy"), t.get("peak_pss_gib"), t.get("mean_pss_gib"),
            t.get("peak_system_used_gib"), t.get("mean_system_used_gib"),
            t.get("max_single_process_rss_gib"), t.get("stage_seconds"),
            stats.get("precursors"), stats.get("protein_groups"), stats.get("gradient_seconds"),
            *stats.get("mass", (None,) * 6),
            sw.get("repo_commit"), dirty, sw.get("uv_lock_sha256"),
            json.dumps(pk),
            *_method_columns(inp[0].get("method")),
            None, None, None, None,  # raw provenance: _search_finalize
        ),
    )  # fmt: skip
    return True


def _no_finalize(con: sqlite3.Connection, bucket, raw_bucket) -> dict:
    return {}


# name -> (schema, units, ingest, delete, finalize). ``finalize`` runs after every build on the
# whole database (kept units included) and returns counts for the build statistics.
_SPECS = {
    "dia_qc": (_QC_SCHEMA, _qc_units, _qc_ingest, _qc_delete, _qc_finalize),
    "flow": (_FLOW_SCHEMA, _flow_units, _flow_ingest, _flow_delete, _no_finalize),
    "search": (_SEARCH_SCHEMA, _search_units, _search_ingest, _search_delete, _search_finalize),
}
# Lookup caches carried over when a snapshot is rebuilt from scratch (schema change, --full).
_CARRIED_CACHES = {"search": ("raw_objects",)}
_MAIN_TABLE = {"dia_qc": "records", "flow": "flow_runs", "search": "search_runs"}
_BOOKKEEPING = """
create table if not exists sources (unit_key text primary key, fingerprint text, ingested text);
create table if not exists build_errors (unit_key text, error text, at text);
create table if not exists meta (key text primary key, value text);
"""


def _schema_version(con: sqlite3.Connection) -> int:
    try:
        row = con.execute("select value from meta where key = 'schema_version'").fetchone()
    except sqlite3.DatabaseError:
        return 0
    return int(row[0]) if row else 0


def _start(name: str, path: Path, bucket, etag: str | None, full: bool) -> sqlite3.Connection:
    """The previous snapshot to update, or a fresh database when there is none / it is stale.

    A fresh database keeps the previous snapshot's lookup caches (``_CARRIED_CACHES``), so a
    schema change does not mean looking every raw file up again.
    """
    previous = None
    if etag:
        previous = path.with_suffix(".previous")
        shutil.copy(bucket.download(db_key(name)), previous)
    if previous and not full:
        shutil.copy(previous, path)
        con = sqlite3.connect(path)
        if _schema_version(con) == SCHEMA_VERSION:
            return con
        con.close()  # older schema: rebuild everything
    path.unlink(missing_ok=True)
    con = sqlite3.connect(path)
    con.executescript(_SPECS[name][0])
    if previous:
        _carry_caches(con, name, previous)
    return con


def _carry_caches(con: sqlite3.Connection, name: str, previous: Path) -> None:
    if not _CARRIED_CACHES.get(name):
        return
    con.execute("attach database ? as old", [str(previous)])  # a local copy, never uploaded
    try:
        for table in _CARRIED_CACHES.get(name, ()):
            old = {
                r[0] for r in con.execute("select name from pragma_table_info(?, 'old')", [table])
            }
            new = [r[1] for r in con.execute(f"pragma table_info({table})")]
            if old and set(new) <= old:
                cols = ",".join(new)
                con.execute(
                    f"insert or ignore into {table} ({cols}) select {cols} from old.{table}"
                )
        con.commit()
    finally:
        con.execute("detach database old")


def _update(
    con: sqlite3.Connection, name: str, bucket, listing: dict[str, str], raw_bucket=None
) -> dict:
    _, units_of, ingest, delete, finalize = _SPECS[name]
    con.isolation_level = None  # explicit transactions: one savepoint per unit
    con.executescript(_BOOKKEEPING)
    con.execute("delete from build_errors")
    con.commit()
    units = units_of(listing)
    have = dict(con.execute("select unit_key, fingerprint from sources"))
    stats = {"kept": 0, "ingested": 0, "skipped": 0, "removed": 0, "failed": 0}
    for gone in set(have) - set(units):
        delete(con, gone)
        con.execute("delete from sources where unit_key = ?", [gone])
        stats["removed"] += 1
    con.commit()
    for unit, fp in sorted(units.items()):
        if have.get(unit) == fp:
            stats["kept"] += 1
            continue
        con.execute("savepoint unit")
        try:
            delete(con, unit)
            produced = ingest(con, bucket, unit, listing)
            con.execute("insert or replace into sources values (?,?,?)", [unit, fp, _now()])
            con.execute("release unit")
            stats["ingested" if produced else "skipped"] += 1
        except Exception as exc:  # one bad unit must not stop the index for everyone
            con.execute("rollback to unit")
            con.execute("release unit")
            con.execute("delete from sources where unit_key = ?", [unit])  # retry next build
            con.execute(
                "insert into build_errors values (?,?,?)",
                [unit, f"{type(exc).__name__}: {exc}", _now()],
            )
            stats["failed"] += 1
        con.commit()
    provenance = finalize(con, bucket, raw_bucket)
    n = con.execute(f"select count(*) from {_MAIN_TABLE[name]}").fetchone()[0]
    con.executemany(
        "insert or replace into meta values (?,?)",
        [
            ("database", name),
            ("schema_version", str(SCHEMA_VERSION)),
            ("built", _now()),
            ("records", str(n)),
            ("scanned_prefix", ",".join(SCAN_PREFIXES[name])),
            *((f"last_build_{k}", str(v)) for k, v in stats.items()),
            *((f"raw_provenance_{k}", str(v)) for k, v in provenance.items()),
        ],
    )
    con.commit()
    return stats | {"records": n, "raw_provenance": provenance}


def _is_conflict(exc: Exception) -> bool:
    code = str(getattr(exc, "response", {}).get("Error", {}).get("Code", ""))
    return code in {"PreconditionFailed", "ConditionalRequestConflict", "412", "409"}


def build(
    name: str,
    *,
    bucket=None,
    out: str | os.PathLike | None = None,
    upload: bool = True,
    full: bool = False,
    raw_bucket=None,
) -> dict:
    """Update one database from the experiment folders (see the module docstring).

    Returns build statistics (kept / ingested / removed / failed units, records, path, and
    ``raw_provenance`` counts by status). With ``upload`` the snapshot is replaced only if
    nobody else replaced it meanwhile; on a conflict the build restarts from the newer
    snapshot, up to ``UPLOAD_ATTEMPTS`` times. ``raw_bucket`` is where searched raw files are
    archived (default ``deliverome_analysis.s3.raw``; only the search index uses it).
    """
    if bucket is None:
        bucket = _bucket()
        if raw_bucket is None:
            raw_bucket = _raw_bucket()
    workdir = Path(tempfile.mkdtemp(prefix="qcdb-"))
    for attempt in range(1, UPLOAD_ATTEMPTS + 1):
        etag = bucket.etag(db_key(name))
        listing = {}
        for prefix in SCAN_PREFIXES[name]:
            listing |= bucket.list_etags(prefix=prefix)
        listing = {k: v for k, v in listing.items() if not _superseded(k)}
        path = workdir / f"{name}.sqlite"
        con = _start(name, path, bucket, etag, full)
        try:
            stats = _update(con, name, bucket, listing, raw_bucket)
        finally:
            con.close()
        if out:
            Path(out).parent.mkdir(parents=True, exist_ok=True)
            shutil.copy(path, out)
        stats["path"] = str(out or path)
        if not upload:
            break
        try:
            bucket.upload_if(
                path, db_key(name), if_match=etag, content_type="application/vnd.sqlite3"
            )
            stats["uploaded"] = f"processed://{db_key(name)}"
            break
        except Exception as exc:
            if not _is_conflict(exc) or attempt == UPLOAD_ATTEMPTS:
                raise
            print(f"{name}: another build uploaded first; rebuilding on top of it ({attempt})")
    print(
        f"{name}: {stats['records']} records ({stats['ingested']} read, {stats['kept']} kept, "
        f"{stats['removed']} removed, {stats['failed']} failed)"
        + (f" -> {stats['uploaded']}" if "uploaded" in stats else "")
    )
    _warn_provenance(name, stats.get("raw_provenance") or {})
    return stats


def _warn_provenance(name: str, counts: dict) -> None:
    """Say loudly when indexed numbers come from a raw file that is not the archived one."""
    from deliverome_analysis.raw_provenance import PROBLEM_STATUSES

    bad = {k: v for k, v in counts.items() if k in PROBLEM_STATUSES and v}
    if bad:
        print(
            f"::warning::{name}: runs whose searched raw file differs from raw:// or is "
            f"missing there: {bad}; see the raw_provenance column"
        )


def _superseded(key: str) -> bool:
    return f"/{SUPERSEDED_DIR}/" in f"/{key}"


def connect(name: str, *, bucket=None, path: str | os.PathLike | None = None) -> sqlite3.Connection:
    """Read-only connection to the current snapshot (or a local file given as ``path``)."""
    local = Path(path) if path else Path((bucket or _bucket()).download(db_key(name)))
    return sqlite3.connect(f"file:{local}?mode=ro", uri=True)


# ----------------------------------------------------------------------------------------------
# readers: the common questions as DataFrames, so nobody re-derives the SQL
# ----------------------------------------------------------------------------------------------


def read(name: str, sql: str, params: list | tuple = (), *, bucket=None, path=None) -> pd.DataFrame:
    """Run one read-only query against a database snapshot."""
    con = connect(name, bucket=bucket, path=path)
    try:
        return pd.read_sql(sql, con, params=list(params))
    finally:
        con.close()


def _where(clauses: dict[str, object]) -> tuple[str, list]:
    """(" where a = ? and b >= ?", params) from {"a = ?": v, ...}, skipping None values."""
    used = {k: v for k, v in clauses.items() if v is not None}
    if not used:
        return "", []
    return " where " + " and ".join(used), list(used.values())


def qc_runs(*, instrument: str | None = None, since: str | None = None, **kw) -> pd.DataFrame:
    """One row per acquisition (newest record): run, run_uid, acquired, instrument, record."""
    where, params = _where({"instrument_serial = ?": instrument, "acquired >= ?": since})
    return read("dia_qc", f"select * from latest_runs{where} order by acquired", params, **kw)


UNVERSIONED = "unversioned"  # how qc_metrics names records from before definitions versioning


def qc_metrics(
    metrics: list[str] | None = None,
    *,
    instrument: str | None = None,
    since: str | None = None,
    wide: bool = False,
    definitions_version: str | None = None,
    allow_mixed_definitions: bool = False,
    **kw,
) -> pd.DataFrame:
    """Per-acquisition QC metrics (long: run, metric, value, acquired, definitions_version,
    ...; ``wide``: one row per acquisition, one column per metric). ``metrics`` limits to
    those names.

    Metrics defined by :mod:`deliverome_analysis.peptide_metrics` carry the
    ``DEFINITIONS_VERSION`` they were computed under (``''`` for other metrics, NULL for
    records from before versioning). ``definitions_version`` keeps only numbers known to be
    computed under that version (plus metrics the module does not define). A ``wide`` table
    puts one metric's values side by side, so it refuses (ValueError) to mix definition
    versions of a metric unless ``allow_mixed_definitions``; it has a ``definitions_version``
    column per acquisition either way.
    """
    from deliverome_analysis.peptide_metrics import metric_names

    where, params = _where({"instrument_serial = ?": instrument, "acquired >= ?": since})
    if metrics:
        where += (" and " if where else " where ") + f"metric in ({','.join('?' * len(metrics))})"
        params += list(metrics)
    df = read("dia_qc", f"select * from latest_run_metrics{where} order by acquired", params, **kw)
    if definitions_version is not None and not df.empty:
        v = df["definitions_version"]
        df = df[
            (v == definitions_version) | (v == "") | (v.isna() & ~df["metric"].isin(metric_names()))
        ].reset_index(drop=True)
    if not wide or df.empty:
        return df
    versioned = df[df["definitions_version"] != ""]  # NULL counts as its own version
    mixed = (
        versioned.assign(v=versioned["definitions_version"].fillna(UNVERSIONED))
        .groupby("metric")["v"]
        .unique()
    )
    mixed = mixed[mixed.map(len) > 1]
    if len(mixed) and not allow_mixed_definitions:
        raise ValueError(
            "metrics computed under different peptide_metrics definitions would share a "
            f"column: {dict(mixed.map(sorted))}. Pass definitions_version=... to keep one, "
            "or allow_mixed_definitions=True (the definitions_version column then says which)."
        )
    df["definitions_version"] = df["record_definitions_version"].fillna(UNVERSIONED)
    index = ["run_uid", "run", "acquired", "instrument_serial", "condition", "definitions_version"]
    df[index] = df[index].fillna("")  # pivot_table silently drops rows with an empty key
    return df.pivot_table(
        index=index, columns="metric", values="value", aggfunc="first"
    ).reset_index()


def qc_diagnoses(*, flagged_only: bool = True, **kw) -> pd.DataFrame:
    """Diagnosis flags, hints and recorded verdicts for the newest record of each acquisition."""
    sql = (
        "select lr.run_uid, lr.run, lr.acquired, lr.instrument_serial, d.flags, d.hints, "
        "d.confirmed_cause, d.notes, r.folder from diagnoses d "
        "join latest_runs lr using (record_id, run) join records r using (record_id)"
    )
    if flagged_only:
        sql += " where d.flags != '' or d.confirmed_cause != ''"
    return read("dia_qc", sql + " order by lr.acquired", **kw)


def flow_runs(*, experiment: str | None = None, **kw) -> pd.DataFrame:
    """Published flow runs (without the gate YAML text; select it from flow_runs if needed)."""
    where, params = _where({"experiment = ?": experiment})
    cols = (
        "run_id, experiment, producer, agentflow_version, agentflow_commit, recipe_sha256, gates_sha256, "
        "input_count, signal_space, events_key"
    )
    return read("flow", f"select {cols} from flow_runs{where} order by run_id", params, **kw)


def search_runs(*, host: str | None = None, status: str | None = "ok", **kw) -> pd.DataFrame:
    """One row per DIA search: input, library, computer, runtime, CPU, peak/mean memory, IDs."""
    where, params = _where({"host = ?": host, "status = ?": status})
    return read("search", f"select * from search_runs{where} order by started_utc", params, **kw)


def search_with_qc(
    metrics: list[str] | None = None,
    *,
    host: str | None = None,
    library: str | None = None,
    wide: bool = True,
    bucket=None,
    search_path=None,
    qc_path=None,
) -> pd.DataFrame:
    """Each search next to the QC metrics of the same acquisition and the same library.

    Joins ``search.sqlite`` (runtime, CPU, memory, computer, IDs) with ``dia_qc.sqlite`` on
    ``run_uid`` (instrument serial | acquisition time from the raw header). A QC run that
    names its search library only matches searches with that library (by SHA-256), so a raw
    file searched with two libraries pairs each search with its own QC page; QC runs without
    a library (Spectronaut pages) match every search of the acquisition. When several QC
    records match, the experiment-folder one wins, then the newest. ``wide``: one row per
    search with a ``qc:<metric>`` column per metric; otherwise long (metric, value).
    """
    if search_path is None:
        search_path = (bucket or _bucket()).download(db_key("search"))
    if qc_path is None:
        qc_path = (bucket or _bucket()).download(db_key("dia_qc"))
    con = sqlite3.connect(f"file:{search_path}?mode=ro", uri=True)
    try:
        con.execute("attach database ? as qc", [f"file:{qc_path}?mode=ro"])
        where, params = _where({"s.host = ?": host, "s.library = ?": library})
        if metrics:
            where += (" and " if where else " where ") + (
                f"m.metric in ({','.join('?' * len(metrics))})"
            )
            params += list(metrics)
        sql = f"""
            with matched as (
              select s.run_key, r.record_id, r.run,
                     row_number() over (
                       partition by s.run_key
                       order by rec.scope = 'experiment' desc, rec.created desc) as rank
              from search_runs s
              join qc.runs r on r.run_uid = s.run_uid
              join qc.records rec using (record_id)
              where r.search_library_sha256 is null
                 or r.search_library_sha256 = s.library_sha256)
            select s.*, rec.folder as qc_folder,
                   rec.definitions_version as qc_definitions_version, m.metric, m.value
            from search_runs s
            join matched x on x.run_key = s.run_key and x.rank = 1
            join qc.records rec on rec.record_id = x.record_id
            join qc.run_metrics m on m.record_id = x.record_id and m.run = x.run
            {where} order by s.started_utc"""
        df = pd.read_sql(sql, con, params=params)
    finally:
        con.close()
    if not wide or df.empty:
        return df
    base = [c for c in df.columns if c not in ("metric", "value")]
    qc_wide = df.pivot_table(index="run_key", columns="metric", values="value", aggfunc="first")
    qc_wide.columns = [f"qc:{c}" for c in qc_wide.columns]
    return df[base].drop_duplicates("run_key").merge(qc_wide.reset_index(), on="run_key")


def _with_metadata(df: pd.DataFrame, metadata: bool, kw: dict) -> pd.DataFrame:
    if not metadata or df.empty:
        return df
    md = read("flow", "select run_id, sample_id, key, value from sample_metadata", **kw)
    if md.empty:
        return df
    wide = md.pivot_table(
        index=["run_id", "sample_id"], columns="key", values="value", aggfunc="first"
    ).reset_index()
    return df.merge(wide, on=["run_id", "sample_id"], how="left")


def flow_populations(
    *,
    experiment: str | None = None,
    run_id: str | None = None,
    gate: str | None = None,
    metadata: bool = True,
    **kw,
) -> pd.DataFrame:
    """Population counts and percentages per run x sample x gate, with sample metadata columns."""
    where, params = _where(
        {"f.experiment = ?": experiment, "p.run_id = ?": run_id, "p.gate = ?": gate}
    )
    df = read(
        "flow",
        f"select f.experiment, p.* from populations p join flow_runs f using (run_id){where}",
        params,
        **kw,
    )
    return _with_metadata(df, metadata, kw)


def flow_stats(
    *,
    experiment: str | None = None,
    run_id: str | None = None,
    gate: str | None = None,
    detector: str | None = None,
    stat: str | None = "median",
    wide: bool = False,
    metadata: bool = False,
    **kw,
) -> pd.DataFrame:
    """Signal statistics (long: run_id, sample_id, gate, detector, stat, value). ``wide``: one
    column per detector (pick a single ``stat``)."""
    where, params = _where(
        {
            "f.experiment = ?": experiment,
            "s.run_id = ?": run_id,
            "s.gate = ?": gate,
            "s.detector = ?": detector,
            "s.stat = ?": stat,
        }
    )
    df = read(
        "flow",
        f"select f.experiment, s.* from signal_stats s join flow_runs f using (run_id){where}",
        params,
        **kw,
    )
    if wide and not df.empty:
        if df["stat"].nunique() > 1:
            raise ValueError("wide=True needs a single stat (e.g. stat='median')")
        df = df.pivot_table(
            index=["experiment", "run_id", "sample_id", "gate"],
            columns="detector",
            values="value",
            aggfunc="first",
        ).reset_index()
    return _with_metadata(df, metadata, kw)


def main(argv: list[str] | None = None) -> int:
    import argparse

    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("build", help="update and upload; exits 1 if any unit failed")
    b.add_argument("name", choices=[*DATABASES, "all"])
    b.add_argument("--full", action="store_true", help="rebuild from scratch")
    b.add_argument("--no-upload", action="store_true")
    b.add_argument("--out", default=None, help="also keep a copy in this directory")
    q = sub.add_parser("query", help="run SQL against the current snapshot")
    q.add_argument("name", choices=DATABASES)
    q.add_argument("sql")
    a = ap.parse_args(argv)
    if a.cmd == "query":
        con = connect(a.name)
        print(pd.read_sql(a.sql, con).to_string(index=False))
        con.close()
        return 0
    failed = 0
    for name in BUILD_ORDER if a.name == "all" else (a.name,):
        out = Path(a.out) / f"{name}.sqlite" if a.out else None
        stats = build(name, out=out, upload=not a.no_upload, full=a.full)
        if stats["failed"]:
            con = connect(name, path=stats["path"])
            for unit, err in con.execute("select unit_key, error from build_errors"):
                print(f"::error::{name}: {unit}: {err}")
            con.close()
            failed += stats["failed"]
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
