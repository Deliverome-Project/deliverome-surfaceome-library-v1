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
import os
import re
import shutil
import sqlite3
import sys
import tempfile
import time
from pathlib import Path

import pandas as pd

SCHEMA_VERSION = 3  # of the SQLite files; a different value forces a full rebuild
RECORD_SCHEMA_VERSION = 1  # of <page>.qcrecord.json, written by dia_qc_dashboard.qc_record
DATABASES = ("dia_qc", "flow")
QC_RECORD_SUFFIX = ".qcrecord.json"
DIAGNOSIS_SUFFIX = ".run_diagnostics.tsv"
FLOW_MANIFEST = "publication.json"
SCAN_PREFIX = "sprint"  # experiment folders: sprint{N}/s{N}e{M}-slug/...
UPLOAD_ATTEMPTS = 3

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
  limpa_version text, versions_json text, schema_version integer);
-- run_uid identifies an acquisition across datasets: instrument serial + acquisition time from
-- the raw-file header, or the run name alone when the record has no raw-file metadata.
create table runs (
  record_id text, run text, run_uid text, label text, condition text, acquired text,
  instrument_model text, instrument_name text, instrument_serial text, file_bytes integer,
  primary key (record_id, run));
create table run_metrics (record_id text, run text, metric text, value real);
create table diagnoses (
  record_id text, run text, flags text, hints text, confirmed_cause text, notes text);
create index run_metrics_metric on run_metrics(metric);
create index runs_uid on runs(run_uid);
-- The newest record per acquisition: the views to use for history plots and queries.
create view latest_runs as
  select * from (
    select r.*, rec.created,
           row_number() over (partition by r.run_uid order by rec.created desc) as rank
    from runs r join records rec using (record_id))
  where rank = 1;
create view latest_run_metrics as
  select m.*, lr.run_uid, lr.acquired, lr.instrument_serial, lr.condition
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
    con.execute(
        "insert into records values (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
        (
            rid, key, key.rsplit("/", 1)[0], rec.get("created"), rec.get("source"),
            rec.get("title"), rec.get("layout"), rec.get("layout_sha256"), rec.get("layout_yaml"),
            prov.get("git-commit"), prov.get("git-dirty"), prov.get("uv-lock-sha256"),
            prov.get("renv-lock-sha256"), lim.get("engine"), lim.get("limpa_version"),
            json.dumps(v), rec.get("schema_version"),
        ),
    )  # fmt: skip
    for r in rec.get("runs", []):
        con.execute(
            "insert into runs values (?,?,?,?,?,?,?,?,?,?)",
            (
                rid, r["run"], run_uid(r), r.get("label"), r.get("condition"), r.get("acquired"),
                r.get("instrument_model"), r.get("instrument_name"), r.get("instrument_serial"),
                r.get("file_bytes"),
            ),
        )  # fmt: skip
    metrics = rec.get("metrics", [])
    if any(len(m) != 3 or not isinstance(m[2], int | float) for m in metrics):
        raise ValueError("metrics must be [run, metric, number] triples")
    con.executemany(
        "insert into run_metrics values (?,?,?,?)", [(rid, r, m, v) for r, m, v in metrics]
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

_SPECS = {
    "dia_qc": (_QC_SCHEMA, _qc_units, _qc_ingest, _qc_delete),
    "flow": (_FLOW_SCHEMA, _flow_units, _flow_ingest, _flow_delete),
}
_MAIN_TABLE = {"dia_qc": "records", "flow": "flow_runs"}
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
    """The previous snapshot to update, or a fresh database when there is none / it is stale."""
    if etag and not full:
        shutil.copy(bucket.download(db_key(name)), path)
        con = sqlite3.connect(path)
        if _schema_version(con) == SCHEMA_VERSION:
            return con
        con.close()  # older schema: rebuild everything
    path.unlink(missing_ok=True)
    con = sqlite3.connect(path)
    con.executescript(_SPECS[name][0])
    return con


def _update(con: sqlite3.Connection, name: str, bucket, listing: dict[str, str]) -> dict:
    _, units_of, ingest, delete = _SPECS[name]
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
    n = con.execute(f"select count(*) from {_MAIN_TABLE[name]}").fetchone()[0]
    con.executemany(
        "insert or replace into meta values (?,?)",
        [
            ("database", name),
            ("schema_version", str(SCHEMA_VERSION)),
            ("built", _now()),
            ("records", str(n)),
            ("scanned_prefix", SCAN_PREFIX),
            *((f"last_build_{k}", str(v)) for k, v in stats.items()),
        ],
    )
    con.commit()
    return stats | {"records": n}


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
) -> dict:
    """Update one database from the experiment folders (see the module docstring).

    Returns build statistics (kept / ingested / removed / failed units, records, path). With
    ``upload`` the snapshot is replaced only if nobody else replaced it meanwhile; on a
    conflict the build restarts from the newer snapshot, up to ``UPLOAD_ATTEMPTS`` times.
    """
    bucket = bucket or _bucket()
    workdir = Path(tempfile.mkdtemp(prefix="qcdb-"))
    for attempt in range(1, UPLOAD_ATTEMPTS + 1):
        etag = bucket.etag(db_key(name))
        listing = bucket.list_etags(prefix=SCAN_PREFIX)
        path = workdir / f"{name}.sqlite"
        con = _start(name, path, bucket, etag, full)
        try:
            stats = _update(con, name, bucket, listing)
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
    return stats


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


def qc_metrics(
    metrics: list[str] | None = None,
    *,
    instrument: str | None = None,
    since: str | None = None,
    wide: bool = False,
    **kw,
) -> pd.DataFrame:
    """Per-acquisition QC metrics (long: run, metric, value, acquired, ...; ``wide``: one row
    per acquisition, one column per metric). ``metrics`` limits to those names."""
    where, params = _where({"instrument_serial = ?": instrument, "acquired >= ?": since})
    if metrics:
        where += (" and " if where else " where ") + f"metric in ({','.join('?' * len(metrics))})"
        params += list(metrics)
    df = read("dia_qc", f"select * from latest_run_metrics{where} order by acquired", params, **kw)
    if not wide or df.empty:
        return df
    index = ["run_uid", "run", "acquired", "instrument_serial", "condition"]
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
    for name in DATABASES if a.name == "all" else (a.name,):
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
