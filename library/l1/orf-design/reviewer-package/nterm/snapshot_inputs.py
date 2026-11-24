"""Freeze the predictor rows that N-terminal placement reads, keyed by the sequence predicted.

Lab-side, run once; reviewers never run it. It needs the repo ``.env`` (Cloudflare D1, read
only) and Modal (read only), and it writes nothing to D1 or Modal.

What it freezes, for every current L1 member on the "N after SP" rung:

* **SignalP 6.0 and disorder rows from the accessible-surfaceome public D1 database**, at the
  versions the 2026-10-18 interval census pinned (``signalp_public`` sp6_2026_09_28,
  ``disorder_public`` dis_2026_09_28, and the ``*_designed_public`` tables for sequences
  re-predicted under cohort ``l1_needs_prediction``).
* **The sequence each accession-keyed row was predicted on.** ``signalp_public`` and
  ``disorder_public`` carry an accession but no sequence. Their sweeps read the sequence
  from ``topology_public`` by ``uniprot_acc_full`` (accessible-surfaceome
  ``modal/signalp6_app.py`` and ``modal/disorder_app.py``), so every ``topology_public``
  sequence stored under that accession is frozen too. The runner then decides by
  sequence, never by accession, whether a row describes the designed protein.
* **The panel's tagged re-prediction** (``canary_rp_1792371023``: the eight disorder-anchored
  panel arms, tagged and native). On S3 only its submission was published, so its
  ``reprediction.json`` is copied from the ``reprediction-runs`` Modal volume.

**Attribution** (convention of 2026-11-10). Copying the canary's ``reprediction.json`` off the
Modal volume was a one-off, implemented by the Claude Code session ``nterm-reviewer-package``
(PR #254, 2026-11-04). The decision was the agent's, so that FOLR1's 0.982 to 0.981 could be
re-read; Becca did not review it. It cannot be regenerated from S3 alone: the source is Modal
working storage. The copy is kept under this package's dated prefix, and the original stays at
``reprediction-runs/canary_rp_1792371023/`` on the volume. The version pins are the census's
(PR #191), not a new choice.

Which accessions are fetched: the census's own lookup key for members it covered, plus
every accession whose ``topology_public`` sequence hashes to the member's designed protein
(this is how the members added since the census are reached). Being fetched does not make
a row usable; that is the runner's decision, by sequence.

    uv run --locked python snapshot_inputs.py --out ~/nterm-snapshot
    uv run --locked python snapshot_inputs.py --out ~/nterm-snapshot --publish 2026-11-04

``--publish`` uploads the two files under this package's own prefix
(``processed://sprint4/s4e9-orf-design-fg-panel/reviewer-package/nterm/<date>/``) and refuses
to overwrite a key that exists.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import subprocess
import sys
import tempfile
import time
import urllib.error
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
S4E9 = HERE.parents[1]
REPO = S4E9.parents[2]
sys.path.insert(0, str(S4E9))

PREFIX = "sprint4/s4e9-orf-design-fg-panel/"
CENSUS_KEY = PREFIX + "triage/l1_interval_census.tsv"
PUBLISH_PREFIX = PREFIX + "reviewer-package/nterm/"
ROWS_NAME = "predictor_rows.jsonl.gz"
CANARY_RUN = "canary_rp_1792371023"
CANARY_NAME = f"{CANARY_RUN}.reprediction.json"

#: The census's pins. The ``*_public`` tables also hold ortholog runs; an unpinned read
#: would collapse onto an arbitrary one.
SIGNALP_VERSION = "sp6_2026_09_28"
DISORDER_VERSION = "dis_2026_09_28"
SIGNALP_DESIGNED_VERSION = "sp6_2026_10_18"
DISORDER_DESIGNED_VERSION = "dis_2026_10_18"
DESIGNED_COHORT = "l1_needs_prediction"
HUMAN_COHORTS = ("human_canonical", "human_isoforms")

SIGNALP_FIELDS = (
    "prediction sp_probability other_probability cleavage_site cleavage_probability "
    "n_region_start n_region_end h_region_start h_region_end c_region_start c_region_end "
    "organism mode tool_version protein_length"
).split()
DISORDER_FIELDS = (
    "predictor protein_length window_residues sequence_substituted disorder_hex plddt_hex "
    "tool_version"
).split()


def d1_query():
    """Read-only SELECT against the surfaceome public D1 database, from the repo ``.env``."""
    from deliverome_analysis.env import load_env

    load_env()
    account = os.environ["CLOUDFLARE_ACCOUNT_ID"]
    token = os.environ["CLOUDFLARE_API_TOKEN"]
    db = os.environ["CLOUDFLARE_D1_SURFACEOME_PUBLIC_ID"]
    url = f"https://api.cloudflare.com/client/v4/accounts/{account}/d1/database/{db}/query"

    def query(sql: str, params: list) -> list[dict]:
        assert sql.lstrip().lower().startswith("select"), "D1 is read-only here"
        req = urllib.request.Request(
            url,
            data=json.dumps({"sql": sql, "params": params}).encode(),
            headers={"Authorization": f"Bearer {token}", "Content-Type": "application/json"},
        )
        for attempt in range(8):
            try:
                with urllib.request.urlopen(req, timeout=120) as r:
                    payload = json.load(r)
                break
            except urllib.error.HTTPError as exc:  # 429 under a burst of reads; back off
                if exc.code not in (429, 500, 502, 503) or attempt == 7:
                    raise
                time.sleep(2**attempt)
        time.sleep(0.2)
        if not payload.get("success"):
            raise RuntimeError(f"D1 error: {payload.get('errors')}")
        return payload["result"][0]["results"]

    return query


def _chunks(items, n):
    items = sorted(items)
    for i in range(0, len(items), n):
        yield items[i : i + n]


def members() -> dict[str, dict]:
    """Current L1 "N after SP" members with their designed protein and its SHA-256."""
    from deliverome_analysis import library_design as ld
    from deliverome_analysis import topology

    proteins = ld.designed_proteins(S4E9 / "candidate_manifest_tranches.tsv", library="l1")
    with (S4E9 / "l1_tag_paths.tsv").open() as fh:
        lane = [r for r in csv.DictReader(fh, delimiter="\t") if r["placement"] == "N after SP"]
    out = {}
    for r in lane:
        seq = proteins[r["candidate_id"]]["sequence"]
        sha = topology.sequence_sha256(seq)
        assert sha == r["sequence_sha256"], f"{r['candidate_id']}: l1_tag_paths sha is stale"
        out[r["candidate_id"]] = {"sha": sha, "sequence": seq}
    return out


def accessions_to_fetch(lane: dict[str, dict]) -> set[str]:
    """Census lookup keys, plus every human accession whose sequence is a member's protein."""
    from deliverome_analysis import topology
    from deliverome_analysis.s3 import processed

    with open(processed.download(CENSUS_KEY)) as fh:
        census = list(csv.DictReader(fh, delimiter="\t"))
    accs = {r["lookup_key"] for r in census if r["lookup_kind"] == "accession" and r["lookup_key"]}
    snap = topology.load_snapshot()
    shas = {m["sha"] for m in lane.values()}
    hit = snap[snap["sequence_sha256"].isin(shas) & snap["cohort"].isin(HUMAN_COHORTS)]
    return accs | set(hit["uniprot_acc_full"])


def fetch(lane: dict[str, dict], accs: set[str], query) -> list[dict]:
    records: list[dict] = []
    for chunk in _chunks(accs, 20):
        marks = ",".join("?" * len(chunk))
        for r in query(
            "select uniprot_acc_full, topology_version, cohort, sequence from topology_public "
            f"where uniprot_acc_full in ({marks})",
            list(chunk),
        ):
            seq = r.pop("sequence")
            records.append(
                {
                    "kind": "input_sequence",
                    "source": "topology_public",
                    "key": r["uniprot_acc_full"],
                    "version": r["topology_version"],
                    "cohort": r["cohort"],
                    "sequence_sha256": hashlib.sha256(seq.encode()).hexdigest(),
                    "sequence": seq,
                }
            )
    for chunk in _chunks(accs, 40):
        marks = ",".join("?" * len(chunk))
        for r in query(
            f"select uniprot_acc_full, {', '.join(SIGNALP_FIELDS)} from signalp_public "
            f"where signalp_version = ? and uniprot_acc_full in ({marks})",
            [SIGNALP_VERSION, *chunk],
        ):
            key = r.pop("uniprot_acc_full")
            records.append(
                {
                    "kind": "signalp",
                    "source": "signalp_public",
                    "key": key,
                    "version": SIGNALP_VERSION,
                }
                | r
            )
    for chunk in _chunks(accs, 12):
        marks = ",".join("?" * len(chunk))
        for r in query(
            f"select uniprot_acc_full, {', '.join(DISORDER_FIELDS)} from disorder_public "
            f"where disorder_version = ? and uniprot_acc_full in ({marks})",
            [DISORDER_VERSION, *chunk],
        ):
            key = r.pop("uniprot_acc_full")
            records.append(
                {
                    "kind": "disorder",
                    "source": "disorder_public",
                    "key": key,
                    "version": DISORDER_VERSION,
                }
                | r
            )
    shas = {m["sha"] for m in lane.values()}
    for chunk in _chunks(shas, 40):
        marks = ",".join("?" * len(chunk))
        for r in query(
            f"select sequence_sha256, {', '.join(SIGNALP_FIELDS)} from signalp_designed_public "
            f"where signalp_version = ? and cohort = ? and sequence_sha256 in ({marks})",
            [SIGNALP_DESIGNED_VERSION, DESIGNED_COHORT, *chunk],
        ):
            key = r.pop("sequence_sha256")
            records.append(
                {
                    "kind": "signalp",
                    "source": "signalp_designed_public",
                    "key": key,
                    "version": SIGNALP_DESIGNED_VERSION,
                }
                | r
            )
    for chunk in _chunks(shas, 12):
        marks = ",".join("?" * len(chunk))
        for r in query(
            f"select sequence_sha256, {', '.join(DISORDER_FIELDS)} from disorder_designed_public "
            f"where disorder_version = ? and cohort = ? and sequence_sha256 in ({marks})",
            [DISORDER_DESIGNED_VERSION, DESIGNED_COHORT, *chunk],
        ):
            key = r.pop("sequence_sha256")
            records.append(
                {
                    "kind": "disorder",
                    "source": "disorder_designed_public",
                    "key": key,
                    "version": DISORDER_DESIGNED_VERSION,
                }
                | r
            )
    return records


def write_rows(records: list[dict], path: Path) -> None:
    """Sorted JSON lines, gzipped with a zero mtime, so a rerun on the same data is identical."""

    def order(r: dict) -> tuple:
        return tuple(
            r.get(k, "") for k in ("kind", "source", "key", "version", "predictor", "cohort")
        )

    body = "".join(json.dumps(r, sort_keys=True) + "\n" for r in sorted(records, key=order))
    with open(path, "wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as gz:
        gz.write(body.encode())


def copy_canary(dest: Path) -> None:
    """The panel canary's reprediction.json, from the reprediction-runs Modal volume."""
    with tempfile.TemporaryDirectory() as tmp:
        subprocess.run(
            ["modal", "volume", "get", "reprediction-runs", f"{CANARY_RUN}/reprediction.json", tmp],
            check=True,
            capture_output=True,
        )
        dest.write_bytes((Path(tmp) / "reprediction.json").read_bytes())


def publish(out: Path, date: str) -> None:
    from deliverome_analysis.s3 import processed

    for name in (ROWS_NAME, CANARY_NAME):
        key = f"{PUBLISH_PREFIX}{date}/{name}"
        if processed.exists(key):
            raise SystemExit(f"processed://{key} exists; publish under a new date")
        processed.upload(out / name, key)
        print(f"uploaded processed://{key}")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--publish", metavar="DATE", help="upload under nterm/<DATE>/")
    args = ap.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)

    lane = members()
    accs = accessions_to_fetch(lane)
    print(f"{len(lane)} N-after-SP members; {len(accs)} accessions to fetch")
    records = fetch(lane, accs, d1_query())
    write_rows(records, args.out / ROWS_NAME)
    kinds: dict[str, int] = {}
    for r in records:
        kinds[f"{r['kind']}:{r['source']}"] = kinds.get(f"{r['kind']}:{r['source']}", 0) + 1
    print(json.dumps(kinds, indent=1))
    copy_canary(args.out / CANARY_NAME)
    if args.publish:
        publish(args.out, args.publish)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
