"""SQL over every published flow run's events, for reanalysis (re-gating, new statistics).

The flow index (:mod:`deliverome_analysis.qcdb`, ``processed://db/flow.sqlite``) holds runs,
gate recipes, versions, populations and summary statistics. The events themselves stay in each
run's ``events.parquet`` on S3 -- hundreds of millions of rows across runs, so they are not
copied into a database. DuckDB queries them in place, reading only the columns and row groups a
query needs:

    from deliverome_analysis import flowdb
    con = flowdb.connect()                      # index tables + an `events` view
    con.sql("select run_id, count(*) from flow_runs join populations using (run_id) group by 1")
    con.sql('''
        select run_id, sample_id, count(*) filter (where "gate:singlets") as singlets,
               median("signal:FL5-A") filter (where "gate:singlets") as fl5_median
        from events where run_id like '%s3e2%' group by 1, 2
    ''').df()
    df = flowdb.events(con, run_id, columns=["signal:FL5-A"], gate="singlets")  # -> pandas
    flowdb.query('select sample_id, count(*) from events group by 1', runs=[run_id])  # one-shot

Summary-level reads (populations, statistics, metadata, QC history) need no DuckDB: use the
readers in :mod:`deliverome_analysis.qcdb` (``flow_populations``, ``flow_stats``, ...).

Event columns: ``run_id``, ``sample_id``, ``event_index``, ``metadata:*`` (well, plate, group
...), ``raw:<detector>`` (as acquired), ``signal:<detector>`` (in the run's signal space,
compensated when the recipe compensates) and one boolean ``gate:<name>`` per gate. Detectors and
gate names differ between runs; columns are matched by name and are NULL where a run lacks them.

S3 access uses DuckDB's ``httpfs`` extension (downloaded by DuckDB on first use) with the same
credentials as :mod:`deliverome_analysis.s3`; ``local=True`` downloads the files instead.
"""

from __future__ import annotations

import os
import sqlite3
from pathlib import Path

import pandas as pd

from deliverome_analysis import qcdb


def _s3_secret(con) -> None:
    """S3 credentials for DuckDB: the personal S3_* keys, else the standard AWS chain (CI role)."""
    from deliverome_analysis.env import load_env

    load_env()
    con.execute("INSTALL httpfs; LOAD httpfs;")
    region = os.environ.get("S3_REGION") or os.environ.get("AWS_REGION") or "us-east-1"
    if os.environ.get("S3_ACCESS_KEY_ID") and os.environ.get("S3_SECRET_ACCESS_KEY"):
        con.execute(
            "CREATE OR REPLACE SECRET deliverome_s3 (TYPE s3, KEY_ID ?, SECRET ?, REGION ?)",
            [os.environ["S3_ACCESS_KEY_ID"], os.environ["S3_SECRET_ACCESS_KEY"], region],
        )
    else:
        con.execute("INSTALL aws; LOAD aws;")
        con.execute(
            "CREATE OR REPLACE SECRET deliverome_s3 (TYPE s3, PROVIDER credential_chain, REGION ?)",
            [region],
        )


def _load_index(con, index: sqlite3.Connection) -> None:
    tables = [r[0] for r in index.execute("select name from sqlite_master where type='table'")]
    for t in tables:
        df = pd.read_sql(f"select * from {t}", index)
        con.register(f"_{t}_df", df)
        con.execute(f"create or replace table {t} as select * from _{t}_df")
        con.unregister(f"_{t}_df")


def connect(
    *,
    index_path: str | os.PathLike | None = None,
    local: bool = False,
    runs: list[str] | None = None,
    bucket=None,
):
    """DuckDB connection with the flow index tables and an ``events`` view over all runs.

    ``runs``: limit the ``events`` view to these run_ids (faster when you know which you need).
    ``local``: download each events.parquet (cached) instead of reading it from S3.
    """
    import duckdb

    con = duckdb.connect()
    with qcdb.connect("flow", path=index_path, bucket=bucket) as index:
        _load_index(con, index)
    files = con.execute(
        "select run_id, events_key from flow_runs where events_key is not null"
    ).fetchall()
    if runs is not None:
        files = [(r, k) for r, k in files if r in set(runs)]
    if not files:
        return con
    if local:
        b = bucket or qcdb._bucket()
        uris = [(r, str(Path(b.download(k)).resolve())) for r, k in files]
    else:
        _s3_secret(con)
        name = (bucket or qcdb._bucket())._name
        uris = [(r, f"s3://{name}/{k}") for r, k in files]
    con.register("_event_files_df", pd.DataFrame(uris, columns=["run_id", "uri"]))
    con.execute("create or replace table event_files as select * from _event_files_df")
    con.unregister("_event_files_df")
    paths = ", ".join("'" + u.replace("'", "''") + "'" for _, u in uris)
    con.execute(
        f"""create or replace view events as
        select f.run_id, e.* exclude (filename)
        from read_parquet([{paths}], union_by_name = true, filename = true) e
        join event_files f on e.filename = f.uri"""
    )
    return con


def events(
    con,
    run_id: str,
    *,
    sample_id: str | None = None,
    columns: list[str] | None = None,
    gate: str | None = None,
) -> pd.DataFrame:
    """One run's events as pandas: optional sample, columns and gate (members only)."""
    cols = ["sample_id", "event_index", *(columns or [])]
    select = ", ".join(f'"{c}"' for c in cols) if columns else "*"
    where = ["run_id = ?"]
    params: list = [run_id]
    if sample_id is not None:
        where.append("sample_id = ?")
        params.append(sample_id)
    if gate is not None:
        where.append(f'"gate:{gate}"')
    return con.execute(f"select {select} from events where {' and '.join(where)}", params).df()


def gate_counts(con, run_id: str) -> pd.DataFrame:
    """Events per sample in every gate, recomputed from the events (cross-checks populations)."""
    gates = [c[0] for c in con.execute("describe events").fetchall() if c[0].startswith("gate:")]
    sums = ", ".join(f'count(*) filter (where "{g}") as "{g[5:]}"' for g in gates)
    wide = con.execute(
        f"select sample_id, {sums} from events where run_id = ? group by sample_id", [run_id]
    ).df()
    long = wide.melt(id_vars="sample_id", var_name="gate", value_name="count")
    return long[long["count"] > 0].reset_index(drop=True)


def query(sql: str, params: list | None = None, **connect_kw) -> pd.DataFrame:
    """One-shot SQL over the flow index tables and the ``events`` view, as a DataFrame.

    ``connect_kw`` go to :func:`connect` -- pass ``runs=[...]`` to read only those runs' events.
    For several queries, keep a connection from :func:`connect` instead.
    """
    con = connect(**connect_kw)
    try:
        return con.execute(sql, params or []).df()
    finally:
        con.close()
