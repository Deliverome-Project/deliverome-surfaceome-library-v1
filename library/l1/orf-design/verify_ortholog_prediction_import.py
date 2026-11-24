"""SELECT-only before/after verification of the ortholog prediction D1 append.

No database clients are opened until main() is explicitly invoked. Reports contain
counts, keys and content hashes, never credentials or complete prediction payloads.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import sys
import time
from pathlib import Path

OLD_VERSION = "topo_2026_11_05_l1_compara116"
NEW_VERSION = "topo_2026_11_05_ortholog_cds_policy"


def file_sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def content_sha(value):
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode()
    ).hexdigest()


def keyed_hashes(rows, keys):
    result = {}
    for row in rows:
        identity = json.dumps([row[k] for k in keys], separators=(",", ":"))
        if identity in result:
            raise ValueError("Database returned duplicate natural keys")
        result[identity] = content_sha(row)
    return result


class Reader:
    def __init__(self, client, pause):
        self.client = client
        self.pause = pause

    def select(self, sql, params=None):
        if not sql.lstrip().upper().startswith("SELECT ") or ";" in sql.rstrip(";"):
            raise ValueError("Verification permits one SELECT statement only")
        time.sleep(self.pause)
        return self.client.query(sql, params or [])


def load_importers(root):
    sys.path.insert(0, str(root / "src"))
    from accessible_surfaceome.cloud.d1_client import D1Client, D1Config
    from accessible_surfaceome.cloud.gpi_predictions import PREDICTION_COLS, validate_prediction
    from accessible_surfaceome.env import load_env

    path = root / "scripts/cloud/upload_topology_to_d1.py"
    spec = importlib.util.spec_from_file_location("verified_topology_uploader", path)
    uploader = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = uploader
    spec.loader.exec_module(uploader)
    return D1Client, D1Config, load_env, PREDICTION_COLS, validate_prediction, uploader


def expected_export(directory, columns, validate_gpi, uploader):
    topology = [
        json.loads(line)
        for line in (directory / "topology_records.jsonl").read_text().splitlines()
        if line.strip()
    ]
    payload = json.loads((directory / "gpi_records.json").read_text())
    if payload["curation"]:
        raise ValueError("This append must not contain curation changes")
    expected_topology = {}
    for row in topology:
        if (
            row["topology_version"] != NEW_VERSION
            or row["tool_version"] != "deeptmhmm-1.0-academic"
        ):
            raise ValueError("Export version/model differs from the approved native append")
        expected = dict(zip(uploader.COLS, uploader._row_to_params(row), strict=True))
        key = (expected["cohort"], expected["uniprot_acc_full"])
        if key in expected_topology:
            raise ValueError("Export contains duplicate D1 topology keys")
        expected_topology[key] = expected
    expected_gpi, number_keys = {}, set()
    for row in payload["predictions"]:
        validate_gpi(row)
        key = (row["run_id"], row["record_id"])
        number = (row["run_id"], row["prediction_number"])
        if key in expected_gpi or number in number_keys:
            raise ValueError("Export contains duplicate D1 NetGPI keys")
        number_keys.add(number)
        expected_gpi[key] = {c: row[c] for c in columns}
    if not expected_topology or not expected_gpi:
        raise ValueError("Expected a nonempty topology and NetGPI append")
    return expected_topology, expected_gpi


def capture(reader, old_version):
    topology_groups = reader.select(
        "SELECT topology_version, cohort, tool_version, COUNT(*) AS n "
        "FROM topology_public GROUP BY topology_version, cohort, tool_version"
    )
    old_rows = reader.select(
        "SELECT * FROM topology_public WHERE topology_version = ?", [old_version]
    )
    old_release = reader.select(
        "SELECT * FROM topology_release WHERE topology_version = ?", [old_version]
    )
    new_rows = reader.select(
        "SELECT * FROM topology_public WHERE topology_version = ?", [NEW_VERSION]
    )
    new_release = reader.select(
        "SELECT * FROM topology_release WHERE topology_version = ?", [NEW_VERSION]
    )
    gpi_groups = reader.select(
        "SELECT run_id, COUNT(*) AS n FROM netgpi_prediction GROUP BY run_id"
    )
    netgpi = {}
    for group in sorted(gpi_groups, key=lambda r: r["run_id"]):
        rows = reader.select("SELECT * FROM netgpi_prediction WHERE run_id = ?", [group["run_id"]])
        netgpi[group["run_id"]] = {
            "n": len(rows),
            "exact": sum(r["sequence_resolution"] == "exact_input" for r in rows),
            "rows": keyed_hashes(rows, ("run_id", "record_id")),
        }
    curation = {}
    groups = reader.select(
        "SELECT source_sha256, source_kind, COUNT(*) AS n "
        "FROM gpi_curation_record GROUP BY source_sha256, source_kind"
    )
    for group in sorted(groups, key=lambda r: (r["source_sha256"], r["source_kind"])):
        rows = reader.select(
            "SELECT * FROM gpi_curation_record WHERE source_sha256 = ? AND source_kind = ?",
            [group["source_sha256"], group["source_kind"]],
        )
        key = json.dumps([group["source_sha256"], group["source_kind"]], separators=(",", ":"))
        curation[key] = {
            "n": len(rows),
            "rows": keyed_hashes(rows, ("source_sha256", "record_number")),
        }
    return {
        "topology_total": sum(r["n"] for r in topology_groups),
        "topology_groups": {
            json.dumps(
                [r["topology_version"], r["cohort"], r["tool_version"]], separators=(",", ":")
            ): r["n"]
            for r in topology_groups
        },
        "old_topology_rows": keyed_hashes(
            old_rows, ("topology_version", "cohort", "uniprot_acc_full")
        ),
        "old_release_sha256": content_sha(old_release),
        "new_topology_rows": new_rows,
        "new_release": new_release,
        "netgpi_total": sum(g["n"] for g in netgpi.values()),
        "netgpi_runs": netgpi,
        "curation_total": sum(g["n"] for g in curation.values()),
        "curation_groups": curation,
    }


def verify_after(reader, current, before, topology, netgpi, columns):
    for field in ("old_topology_rows", "old_release_sha256", "curation_groups", "curation_total"):
        if current[field] != before[field]:
            raise ValueError(f"Historical records changed: {field}")
    old_groups = {
        k: n for k, n in current["topology_groups"].items() if json.loads(k)[0] != NEW_VERSION
    }
    if old_groups != before["topology_groups"]:
        raise ValueError("Historical topology group counts changed")
    stored = {(r["cohort"], r["uniprot_acc_full"]): r for r in current["new_topology_rows"]}
    if len(current["new_topology_rows"]) != len(topology) or set(stored) != set(topology):
        raise ValueError("Stored new topology key set differs from the complete export")
    for key, expected in topology.items():
        row = stored[key]
        if any(row.get(c) != value for c, value in expected.items()):
            raise ValueError(f"Stored topology columns differ from export: {key}")
        if any(value is not None for c, value in row.items() if c.startswith("dtm2_")):
            raise ValueError("Unexpected DeepTMHMM2 values in a v1 row")
    if current["topology_total"] != before["topology_total"] + len(topology):
        raise ValueError("Topology total changed by an unexpected amount")
    release = current["new_release"]
    cohorts = {r["cohort"] for r in topology.values()}
    if (
        len(release) != 1
        or release[0]["n_rows"] != len(topology)
        or set(release[0]["cohorts_present"].split(",")) != cohorts
        or release[0]["deeptmhmm_version"] != "deeptmhmm-1.0-academic"
    ):
        raise ValueError("New release metadata does not describe the complete export")
    new_runs = {key[0] for key in netgpi}
    existing_runs = {k: v for k, v in current["netgpi_runs"].items() if k not in new_runs}
    if existing_runs != before["netgpi_runs"]:
        raise ValueError("Historical NetGPI run rows changed")
    for run in sorted(new_runs):
        rows = reader.select("SELECT * FROM netgpi_prediction WHERE run_id = ?", [run])
        found = {(r["run_id"], r["record_id"]): r for r in rows}
        expected = {k: r for k, r in netgpi.items() if k[0] == run}
        if len(rows) != len(expected) or set(found) != set(expected):
            raise ValueError("Stored NetGPI keys differ from export")
        for key, row in expected.items():
            if any(found[key].get(c) != row[c] for c in columns):
                raise ValueError(f"Stored NetGPI columns differ from export: {key}")
    if current["netgpi_total"] != before["netgpi_total"] + len(netgpi):
        raise ValueError("NetGPI total changed by an unexpected amount")


def affected_genes(directory):
    index = json.loads((directory / "prediction_provenance_index.json").read_text())
    members = [m for key in ("topology", "netgpi") for row in index[key] for m in row["members"]]
    symbols = sorted({r["gene"] for r in members if r.get("gene")})
    stable_ids = sorted({r["human_ensembl_gene"] for r in members if r.get("human_ensembl_gene")})
    return symbols, stable_ids


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--surfaceome-root", type=Path, required=True)
    parser.add_argument("--export-dir", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument(
        "--before", type=Path, help="compare after-import state against this before report"
    )
    parser.add_argument("--old-topology-version", default=OLD_VERSION)
    parser.add_argument("--pause", type=float, default=0.4, help="seconds before each SELECT")
    args = parser.parse_args(argv)
    if args.out.exists():
        raise FileExistsError("Refusing to overwrite an existing verification report")
    if args.pause < 0.4:
        raise ValueError("Keep D1 requests paced at no more than 2.5 per second")
    args.export_dir = args.export_dir.resolve()
    D1Client, D1Config, load_env, columns, validate_gpi, uploader = load_importers(
        args.surfaceome_root.resolve()
    )
    topology, netgpi = expected_export(args.export_dir, columns, validate_gpi, uploader)
    export_hashes = {
        name: file_sha(args.export_dir / name)
        for name in (
            "topology_records.jsonl",
            "gpi_records.json",
            "prediction_provenance_index.json",
        )
    }
    before = json.loads(args.before.read_text()) if args.before else None
    if before and (
        before["status"] != "passed"
        or before["mode"] != "before"
        or before["export_sha256"] != export_hashes
    ):
        raise ValueError("Before report is not successful or belongs to different export bytes")
    if before and before["old_topology_version"] != args.old_topology_version:
        raise ValueError("Before report uses a different historical topology version")
    report = {
        "mode": "after" if before else "before",
        "status": "running",
        "databases": {},
        "export_sha256": export_hashes,
        "expected_topology": len(topology),
        "expected_netgpi": len(netgpi),
        "old_topology_version": args.old_topology_version,
        "verification_script_sha256": file_sha(__file__),
        "importer_sha256": {
            "topology_uploader": file_sha(
                args.surfaceome_root / "scripts/cloud/upload_topology_to_d1.py"
            ),
            "gpi_validator": file_sha(
                args.surfaceome_root / "src/accessible_surfaceome/cloud/gpi_predictions.py"
            ),
        },
    }
    symbols, stable_ids = affected_genes(args.export_dir)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    (args.out.parent / "affected_human_genes.json").write_text(json.dumps(symbols, indent=2) + "\n")
    (args.out.parent / "affected_human_ensembl_genes.json").write_text(
        json.dumps(stable_ids, indent=2) + "\n"
    )
    try:
        load_env()
        targets = (
            ("surfaceome_agents", D1Config.from_env()),
            ("surfaceome_public", D1Config.from_env_public()),
        )
        if targets[0][1].database_id == targets[1][1].database_id:
            raise ValueError("Private/public D1 configurations must identify distinct databases")
        for name, config in targets:
            with D1Client(config) as client:
                reader = Reader(client, args.pause)
                state = capture(reader, args.old_topology_version)
                if not state["old_topology_rows"]:
                    raise ValueError(f"{name}: historical topology version is absent")
                if before:
                    verify_after(
                        reader, state, before["databases"][name], topology, netgpi, columns
                    )
                elif (
                    state["new_topology_rows"]
                    or state["new_release"]
                    or {k[0] for k in netgpi} & state["netgpi_runs"].keys()
                ):
                    raise ValueError(
                        f"{name}: new prediction version/run is not empty before import"
                    )
                # Compact report: preserve full-row hashes rather than complete sequences.
                state["new_topology_rows"] = keyed_hashes(
                    state["new_topology_rows"], ("cohort", "uniprot_acc_full")
                )
                report["databases"][name] = state
                print(f"{name}: {report['mode']} verification passed", flush=True)
        report["status"] = "passed"
    except Exception as error:
        report["status"] = "failed"
        report["error"] = f"{type(error).__name__}: {error}"
        raise
    finally:
        args.out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
