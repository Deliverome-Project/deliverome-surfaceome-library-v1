#!/usr/bin/env python3
"""Saved copy of the Cloudflare D1 ``surfaceome_public`` tables S5E13 reads, and the query path.

Ten-plus S5E13 scripts read per-residue features (topology, disorder, RSA, pLDDT), the
published tag sites and the homomer list from D1 with credentials. D1 is live: it gains rows,
and nothing outside it records what a script saw. This module freezes the five tables those
scripts query into one SQLite file, published once to an append-only key with its SHA-256,
and makes every caller read that file by default.

    TABLES      the ten in TABLES below. Whole tables, not a row subset: several scripts scan
                by cohort or version, so only the whole table replays every query exactly.
                D1 is SQLite, so the scripts' own SQL runs on the file unchanged; the table and
                index definitions are copied from D1's sqlite_master and rows go in rowid order.

Live D1 is used only behind an explicit flag, ``S5E13_D1_LIVE=1`` in the environment (with
the CLOUDFLARE_* credentials). Without the flag a script never touches the network for D1.

    uv run --locked python d1_snapshot.py restore     # fetch + verify the saved copy into data/
    uv run --locked python d1_snapshot.py verify      # check the local copy's SHA-256
    S5E13_D1_LIVE=1 uv run --locked python d1_snapshot.py build   # (once) dump live D1

Built 2026-11-14 by the ranker-pu-segment session for the L1 reproducibility audit (B6.2).
"""

from __future__ import annotations

import gzip
import hashlib
import json
import os
import shutil
import sqlite3
import sys
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
#: every D1 table an S5E13 script queries (found by matching every table in D1's sqlite_master
#: against the scripts' source). The first dump (2026-11-14) held only the first five and was
#: superseded the same day by 2026-11-14b, which adds the other five.
TABLES = (
    "topology_public",
    "disorder_public",
    "tag_site_public",
    "schweke_homomer_public",
    "library_structure_public",
    "tag_site_prediction",
    "contact_active_release",
    "contact_observation",
    "surface_bind_site",
    "candidate_universe_public",  # surfaceome-union-ok: whole-table copy for the saved D1 dump
)
SNAPSHOT_DATE = "2026-11-14b"
#: external:// is append-only (no delete); a fresh dated prefix is never rewritten
KEY = f"accessible-surfaceome/d1/{SNAPSHOT_DATE}/surfaceome_public_s5e13.sqlite.gz"
LOCAL = HERE / "data" / f"d1_surfaceome_public_{SNAPSHOT_DATE}.sqlite"
MANIFEST = HERE / f"d1_snapshot_{SNAPSHOT_DATE}.json"
#: written into the manifest: the earlier key is incomplete and nothing reads it
SUPERSEDES = {
    "key": "external://accessible-surfaceome/d1/2026-11-14/surfaceome_public_s5e13.sqlite.gz",
    "status": "incomplete and unused",
    "why": "holds only five of the ten queried tables (topology_public, disorder_public, "
    "tag_site_public, schweke_homomer_public, library_structure_public); missing "
    "tag_site_prediction, contact_active_release, contact_observation, surface_bind_site and "
    "candidate_universe_public, so export_site_features failed on it. The key is append-only "  # surfaceome-union-ok: prose naming a table, not a filter
    "and stays in the bucket; no script reads it",
}
LIVE_ENV = "S5E13_D1_LIVE"
PAGE = 250


def live() -> bool:
    """True only when the caller explicitly asked for live D1."""
    return os.environ.get(LIVE_ENV, "") == "1"


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def _manifest() -> dict:
    return json.loads(MANIFEST.read_text())


def _src() -> None:
    p = str(HERE.parents[3] / "src")
    if p not in sys.path:
        sys.path.insert(0, p)


# ---------------------------------------------------------------------------
# live D1 (flagged)
# ---------------------------------------------------------------------------
def live_query(sql: str, params: list | None = None, tries: int = 4) -> list[dict]:
    if not live():
        raise RuntimeError(f"live D1 is off; set {LIVE_ENV}=1 to query it")
    _src()
    from deliverome_analysis.env import load_env

    load_env()
    acct = os.environ["CLOUDFLARE_ACCOUNT_ID"]
    token = os.environ["CLOUDFLARE_API_TOKEN"]
    db = os.environ["CLOUDFLARE_D1_SURFACEOME_PUBLIC_ID"]
    url = f"https://api.cloudflare.com/client/v4/accounts/{acct}/d1/database/{db}/query"
    body = json.dumps({"sql": sql, "params": params or []}).encode()
    last: Exception | None = None
    for _ in range(tries):
        try:
            req = urllib.request.Request(
                url,
                data=body,
                headers={"Authorization": f"Bearer {token}", "Content-Type": "application/json"},
            )
            with urllib.request.urlopen(req, timeout=180) as r:
                payload = json.loads(r.read().decode())
            if not payload.get("success"):
                raise RuntimeError(f"D1 error: {payload.get('errors')}")
            return payload["result"][0]["results"]
        except Exception as e:  # retried, then re-raised
            last = e
    raise RuntimeError(f"D1 query failed after {tries} tries: {last}")


# ---------------------------------------------------------------------------
# the saved copy (default)
# ---------------------------------------------------------------------------
_CONN: sqlite3.Connection | None = None


def connect() -> sqlite3.Connection:
    """Read-only connection to the verified local copy; restores it if absent."""
    global _CONN
    if _CONN is None:
        if not LOCAL.is_file():
            restore()
        _CONN = sqlite3.connect(f"file:{LOCAL}?mode=ro", uri=True)
        _CONN.row_factory = sqlite3.Row
    return _CONN


def _as_d1(v):
    # D1's HTTP API returns JSON, which writes a whole-number REAL as an integer (80, not
    # 80.0). Match it, so pandas dtypes and printed values equal a live run's.
    return int(v) if isinstance(v, float) and v.is_integer() else v


def snapshot_query(sql: str, params: list | None = None) -> list[dict]:
    cur = connect().execute(sql, params or [])
    cols = [d[0] for d in cur.description or []]
    return [{c: _as_d1(v) for c, v in zip(cols, r, strict=True)} for r in cur.fetchall()]


def query(sql: str, params: list | None = None) -> list[dict]:
    """The query every S5E13 script uses: the saved copy, or live D1 when flagged."""
    return live_query(sql, params) if live() else snapshot_query(sql, params)


def restore() -> Path:
    """Download the saved copy, check both SHA-256s, unpack into data/."""
    _src()
    from deliverome_analysis import s3

    m = _manifest()
    gz = LOCAL.with_suffix(".sqlite.gz")
    LOCAL.parent.mkdir(parents=True, exist_ok=True)
    s3.external.download(KEY, str(gz))
    if _sha256(gz) != m["gz_sha256"]:
        raise SystemExit(f"{KEY}: SHA-256 mismatch on the download")
    tmp = LOCAL.with_suffix(".sqlite.part")
    with gzip.open(gz, "rb") as fi, open(tmp, "wb") as fo:
        shutil.copyfileobj(fi, fo)
    if _sha256(tmp) != m["sqlite_sha256"]:
        tmp.unlink()
        raise SystemExit(f"{KEY}: unpacked SQLite does not match its SHA-256")
    tmp.replace(LOCAL)
    gz.unlink()
    return LOCAL


def verify() -> None:
    m = _manifest()
    got = _sha256(LOCAL)
    if got != m["sqlite_sha256"]:
        raise SystemExit(f"{LOCAL.name}: {got} != {m['sqlite_sha256']}")
    con = connect()
    for t, n in m["row_counts"].items():
        k = con.execute(f"SELECT COUNT(*) FROM {t}").fetchone()[0]
        if k != n:
            raise SystemExit(f"{t}: {k} rows, manifest says {n}")
    print(f"{LOCAL.name}: sha256 and row counts match the manifest")


def build() -> None:
    """Dump the five tables from live D1 (flag required) and write the manifest."""
    if not live():
        raise SystemExit(f"build reads live D1: set {LIVE_ENV}=1")
    LOCAL.parent.mkdir(parents=True, exist_ok=True)
    if LOCAL.exists():
        LOCAL.unlink()
    out = sqlite3.connect(LOCAL)
    counts = {}
    for t in TABLES:
        ddl = live_query("SELECT type, name, sql FROM sqlite_master WHERE tbl_name = ?", [t])
        for d in sorted(ddl, key=lambda d: d["type"] != "table"):
            if d["sql"]:
                out.execute(d["sql"])
        n = live_query(f"SELECT COUNT(*) AS n FROM {t}")[0]["n"]
        cols = [c["name"] for c in live_query(f"PRAGMA table_info({t})")]
        ins = f"INSERT INTO {t} ({','.join(cols)}) VALUES ({','.join('?' * len(cols))})"
        got = 0
        while got < n:
            rows = live_query(f"SELECT * FROM {t} ORDER BY rowid LIMIT ? OFFSET ?", [PAGE, got])
            if not rows:
                break
            out.executemany(ins, [[r[c] for c in cols] for r in rows])
            got += len(rows)
            print(f"  {t}: {got}/{n}", end="\r", file=sys.stderr)
        out.commit()
        k = out.execute(f"SELECT COUNT(*) FROM {t}").fetchone()[0]
        again = live_query(f"SELECT COUNT(*) AS n FROM {t}")[0]["n"]
        if not (k == n == again):
            raise SystemExit(f"{t}: dumped {k}, D1 had {n} then {again}; rerun")
        counts[t] = k
        print(f"  {t}: {k} rows", file=sys.stderr)
    out.execute("VACUUM")
    out.close()
    gz = LOCAL.with_suffix(".sqlite.gz")
    with open(LOCAL, "rb") as fi, gzip.GzipFile(gz, "wb", mtime=0) as fo:
        shutil.copyfileobj(fi, fo)
    MANIFEST.write_text(
        json.dumps(
            {
                "source": "Cloudflare D1 surfaceome_public (read-only)",
                "dumped_on": SNAPSHOT_DATE,
                "key": f"external://{KEY}",
                "tables": list(TABLES),
                "row_counts": counts,
                "sqlite_sha256": _sha256(LOCAL),
                "gz_sha256": _sha256(gz),
                "built_by": "d1_snapshot.py build (ranker-pu-segment session)",
                "supersedes": SUPERSEDES,
            },
            indent=2,
        )
        + "\n"
    )
    print(f"wrote {LOCAL.name}, {gz.name} and {MANIFEST.name}")


def publish() -> None:
    """Upload the gzip to its append-only key; refuses if the key already exists."""
    _src()
    from deliverome_analysis import s3

    gz = LOCAL.with_suffix(".sqlite.gz")
    if _sha256(gz) != _manifest()["gz_sha256"]:
        raise SystemExit("local gzip does not match the manifest")
    if s3.external.exists(KEY):
        raise SystemExit(f"external://{KEY} already exists; dated keys are never rewritten")
    s3.external.upload(str(gz), KEY)
    print(f"external://{KEY}")


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "verify"
    {"build": build, "publish": publish, "restore": restore, "verify": verify}[cmd]()
