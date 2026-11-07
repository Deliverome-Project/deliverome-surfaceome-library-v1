"""Export new native ortholog predictions into the existing D1 import schemas.

All validation and output are offline. Runtime-bound v1 supplements are mandatory;
tool versions come from those records. NetGPI exports bind the completed native
wrapper receipt to actual submitted sequences and literal DTU output. No D1 writes,
submissions, snapshot activation or legacy/curation changes occur here.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import defaultdict
from datetime import date
from pathlib import Path

import run_ortholog_netgpi as netgpi
import validated_v1_supplements as v1

TOPOLOGY_VERSION = "topo_2026_11_05_ortholog_cds_policy"
MEMBER_FIELDS = (
    "gene",
    "species",
    "human_ensembl_gene",
    "ortholog_ensembl_gene",
    "ortholog_symbols",
    "picked_acc",
    "picked_source",
    "picked_cds_transcript",
    "picked_cds_sha256",
    "picked_sequence_sha256",
    "candidate_id",
    "member_key",
    "existing_candidate_ids",
)


def derive(topology: str) -> dict:
    """Match the canonical D1 parser: terminal S is skipped, M/B is indeterminate."""
    orientation = {"O": "extracellular", "I": "cytoplasmic"}
    first = topology.lstrip("S")[:1]
    last = topology.rstrip("S")[-1:]
    return {
        "tm_helix_count": len(re.findall("M+", topology)),
        "beta_strand_count": len(re.findall("B+", topology)),
        "signal_peptide_length": len(topology) - len(topology.lstrip("S")),
        "ecd_length_residues": topology.count("O"),
        "icd_length_residues": topology.count("I"),
        "n_terminal_orientation": orientation.get(first, "indeterminate"),
        "c_terminal_orientation": orientation.get(last, "indeterminate"),
    }


def claim_netgpi_keys(prediction: dict, record_keys: set, ordinal_keys: set):
    """Both independent D1 uniqueness constraints must hold before any import."""
    record_key = (prediction["run_id"], prediction["record_id"])
    ordinal_key = (prediction["run_id"], prediction["prediction_number"])
    if record_key in record_keys or ordinal_key in ordinal_keys:
        raise ValueError("Duplicate NetGPI record or prediction ordinal within a run")
    record_keys.add(record_key)
    ordinal_keys.add(ordinal_key)


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def local_path(root: Path, name: str) -> Path:
    path = Path(name)
    if path.is_absolute():
        raise ValueError("Export manifest paths must be relative to root")
    path = (root / path).resolve()
    if not path.is_relative_to(root):
        raise ValueError("Export source escapes its archive root")
    return path


def pick_members(root: Path, tables: list[str]) -> dict:
    result = defaultdict(dict)
    for name in tables:
        path = local_path(root, name)
        with path.open() as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                digest = row.get("picked_sequence_sha256", "")
                if not digest:
                    continue
                if not re.fullmatch(r"[0-9a-f]{64}", digest) or row.get("species") not in {
                    "mouse",
                    "cyno",
                }:
                    raise ValueError("Pick identity lacks exact sequence/species provenance")
                if not row.get("picked_acc") or not row.get("ortholog_ensembl_gene"):
                    raise ValueError("Pick lacks its accession or stable target gene")
                member = {key: row.get(key, "") for key in MEMBER_FIELDS}
                result[digest][json.dumps(member, sort_keys=True)] = member
    return {digest: list(records.values()) for digest, records in result.items()}


def completed_netgpi(root: Path, run: dict, source) -> tuple[list[dict], dict]:
    receipt_path = local_path(root, run["submission"])
    receipt = json.loads(receipt_path.read_text())
    if (
        receipt.get("status") != "completed"
        or not receipt.get("job_id")
        or receipt.get("tool_version") != "NetGPI-1.1"
        or receipt.get("execute_approved") is not True
        or receipt.get("service") != "https://services.healthtech.dtu.dk/cgi-bin/webface2.cgi"
        or receipt.get("config") != "/var/www/services/services/NetGPI-1.1/webface.cf"
    ):
        raise ValueError("NetGPI receipt is not a completed approved native1.1 run")
    if receipt.get("reported_tool_versions") not in ([], ["1.1"]):
        raise ValueError("NetGPI reported version conflicts with the recorded service")
    batch_root = receipt_path.parent
    if not {"submitted.fasta", "raw_output_protein_type.txt", "predictions.json"} <= set(
        receipt["output_sha256"]
    ):
        raise ValueError("NetGPI receipt must checksum the input, raw output and predictions")
    for relative, expected in receipt["output_sha256"].items():
        path = local_path(batch_root, relative)
        if sha(path) != expected:
            raise ValueError(f"NetGPI recorded output checksum differs: {relative}")
        source(path)
    fasta = batch_root / "submitted.fasta"
    raw = batch_root / "raw_output_protein_type.txt"
    prediction_path = local_path(root, run["predictions"])
    if prediction_path != (batch_root / "predictions.json").resolve():
        raise ValueError("NetGPI receipt and prediction file belong to different runs")
    submitted = v1.submitted_fasta(fasta)
    if sha(fasta) != receipt["input_sha256"] or len(submitted) != receipt["sequence_count"]:
        raise ValueError("NetGPI actual submitted FASTA/count differs from its receipt")
    records = receipt["records"]
    if len(records) != len(submitted) or {r["record_id"] for r in records} != set(submitted):
        raise ValueError("NetGPI receipt does not exactly cover submitted records")
    for ordinal, (identity, protein) in enumerate(submitted.items(), 1):
        matches = [r for r in records if r["record_id"] == identity]
        if (
            len(matches) != 1
            or matches[0]["input_number"] != ordinal
            or matches[0]["sequence_sha256"] != identity
            or matches[0]["sequence_length"] != len(protein)
        ):
            raise ValueError("NetGPI input ordinal/sequence/length differs from submission")
    raw_text = raw.read_text()
    raw_rows = [line for line in raw_text.splitlines() if line and not line.startswith("#")]
    predictions = json.loads(prediction_path.read_text())
    if predictions != netgpi.parse_predictions(raw_rows, receipt):
        raise ValueError("NetGPI predictions differ from literal numbered DTU output")
    if receipt.get("prediction_count") != len(predictions):
        raise ValueError("NetGPI completed count differs from its actual predictions")
    source(receipt_path)
    raw_uri = source(raw)
    header = "\n".join(line for line in raw_text.splitlines() if line.startswith("#"))
    expanded = [
        {
            **prediction,
            "sequence": submitted[prediction["record_id"]],
            "raw_header": header,
            "source_uri": raw_uri,
            "source_sha256": sha(raw),
        }
        for prediction in predictions
    ]
    provenance = {
        "run_id": receipt["run_id"],
        "job_id": receipt["job_id"],
        "tool_version": receipt["tool_version"],
        "submission_uri": source(receipt_path),
        "submission_sha256": sha(receipt_path),
        "input_sha256": receipt["input_sha256"],
        "prediction_file_uri": source(prediction_path),
        "prediction_file_sha256": sha(prediction_path),
        "receipt": receipt,
    }
    return expanded, provenance


def export(root: Path, manifest_path: Path, out: Path, import_code_sha: str) -> dict:
    root = root.resolve()
    manifest = json.loads(manifest_path.read_text())
    if manifest.get("topology_version") != TOPOLOGY_VERSION:
        raise ValueError("Use the dedicated immutable ortholog CDS-policy topology version")
    date.fromisoformat(manifest["date"])
    if not re.fullmatch(r"[0-9a-f]{40}", import_code_sha):
        raise ValueError("Import provenance requires a full Git commit")
    prefix = manifest["s3_prefix"]
    if not prefix.startswith("sprint4/s4e9-orf-design-fg-panel/") or not prefix.endswith("/"):
        raise ValueError("Use a fresh experiment-scoped S3 prefix")
    sources = {}

    def source(path: Path) -> str:
        path = path.resolve()
        if not path.is_relative_to(root):
            raise ValueError("Bound source is outside the archive root")
        name = path.relative_to(root).as_posix()
        actual = sha(path)
        previous = sources.get(name)
        if previous and previous["sha256"] != actual:
            raise ValueError("Source changed while exporting")
        sources[name] = {"sha256": actual, "source": "processed://" + prefix + name}
        return sources[name]["source"]

    members = pick_members(root, manifest["pick_tables"])
    for name in manifest["pick_tables"]:
        source(local_path(root, name))
    topology_runs = manifest.get("topology", [])
    outputs = [local_path(root, run["output"]) for run in topology_runs]
    records = [local_path(root, run["run_record"]) for run in topology_runs]
    calls, validation = v1.load(outputs, records)
    if any(r["state"] != "duplicate_consistent" for r in validation["duplicate_reports"]):
        raise ValueError(
            "Conflicting per-residue calls require separate version/curation; refuses silent selection"
        )
    run_sources = {}
    for record_path in records:
        run = json.loads(record_path.read_text())
        if run["run_id"] in run_sources:
            raise ValueError("Duplicate topology run ID across listed rounds")
        bound = {
            key: v1.bound_file(record_path, run[key])
            for key in ("submitted_fasta", "output", "runtime_identity")
        }
        run_sources[run["run_id"]] = {
            "run_record_uri": source(record_path),
            "run_record_sha256": sha(record_path),
            **{
                key: {"source_uri": source(path), "sha256": sha(path)}
                for key, path in bound.items()
            },
        }
        for name in run.get("extra_source_files", []):
            source(local_path(record_path.parent, name))
    topology_rows, topology_index = {}, {}
    for species in ("mouse", "cyno"):
        for digest, call in sorted(calls[species].items()):
            matching = [m for m in members.get(digest, []) if m["species"] == species]
            if not matching:
                raise ValueError(
                    f"Predicted protein lacks matching source species/pick identity: {species}/{digest}"
                )
            versions = {e["tool_version"] for e in call["provenance"]}
            runtime_versions = {e["model_identity"]["tool_version"] for e in call["provenance"]}
            if len(versions) != 1 or versions != runtime_versions:
                raise ValueError("Validated run/runtime model versions disagree")
            tool_version = next(iter(versions))
            for member in matching:
                accession = member["picked_acc"]
                row = {
                    "topology_version": TOPOLOGY_VERSION,
                    "cohort": species + "_ortholog_cds_policy",
                    "hgnc_id": None,
                    "uniprot_accession": accession.split("-")[0],
                    "uniprot_accession_full": accession,
                    "isoform_id": accession,
                    "gene_symbol": member["ortholog_symbols"],
                    "species": "cynomolgus" if species == "cyno" else "mouse",
                    "is_canonical": 0,
                    "sequence": call["sequence"],
                    "protein_length": len(call["sequence"]),
                    "deeptmhmm_label": call["label"],
                    "per_residue_topology": call["topology"],
                    "predicted_surface_membrane": int(call["label"] in {"TM", "SP+TM"}),
                    "predicted_secreted": int(call["label"] == "SP"),
                    "tool_version": tool_version,
                    "retrieved_at": manifest["date"],
                    **derive(call["topology"]),
                }
                key = (TOPOLOGY_VERSION, row["cohort"], accession)
                if key in topology_rows and topology_rows[key] != row:
                    raise ValueError(
                        f"One D1 topology key has inconsistent sequence/identity: {key}"
                    )
                topology_rows[key] = row
                index = topology_index.setdefault(
                    key,
                    {
                        "row_key": list(key),
                        "sequence_sha256": digest,
                        "members": [],
                        "prediction_evidence": [],
                    },
                )
                if member not in index["members"]:
                    index["members"].append(member)
                for evidence in call["provenance"]:
                    item = {k: v for k, v in evidence.items() if k != "model_identity"}
                    item["source_files"] = run_sources[evidence["run_id"]]
                    if item not in index["prediction_evidence"]:
                        index["prediction_evidence"].append(item)
    gpi_predictions, gpi_index = [], []
    gpi_run_ids, gpi_record_keys, gpi_ordinal_keys = set(), set(), set()
    for run in manifest.get("netgpi", []):
        predictions, provenance = completed_netgpi(root, run, source)
        if provenance["run_id"] in gpi_run_ids:
            raise ValueError("Duplicate NetGPI run ID across listed receipts")
        gpi_run_ids.add(provenance["run_id"])
        for p in predictions:
            digest = p["sequence_sha256"]
            matching = members.get(digest, [])
            if not matching:
                raise ValueError("Exact NetGPI protein lacks source ortholog identity")
            receipt_record = next(
                r for r in provenance["receipt"]["records"] if r["record_id"] == p["record_id"]
            )
            submitted_species = {
                e["pick"]["species"] for e in receipt_record["native_ortholog_provenance"]
            }
            if not submitted_species <= {m["species"] for m in matching}:
                raise ValueError("NetGPI submitted species provenance differs from pick tables")
            for entry in receipt_record["native_ortholog_provenance"]:
                submitted_member = entry["pick"]
                identity_fields = (
                    "species",
                    "gene",
                    "ortholog_ensembl_gene",
                    "picked_acc",
                    "picked_cds_transcript",
                    "picked_cds_sha256",
                )
                if not any(
                    all(member.get(k) == submitted_member.get(k, "") for k in identity_fields)
                    for member in matching
                ):
                    raise ValueError(
                        "NetGPI submitted gene/CDS/member provenance differs from pick tables"
                    )
            key = (p["run_id"], p["record_id"], p["prediction_number"])
            claim_netgpi_keys(p, gpi_record_keys, gpi_ordinal_keys)
            gpi_predictions.append(
                {
                    "run_id": p["run_id"],
                    "record_id": p["record_id"],
                    "prediction_number": p["prediction_number"],
                    "tool_version": p["tool_version"],
                    "job_id": p["job_id"],
                    "sequence_sha256": digest,
                    "sequence": p["sequence"],
                    "sequence_resolution": "exact_input",
                    "sequence_length": p["sequence_length"],
                    "prediction": p["prediction"],
                    "likelihood_text": p["likelihood_text"],
                    "omega_position": int(p["omega_position_text"])
                    if p["omega_position_text"].isdigit()
                    else None,
                    "omega_residue": p["omega_residue"],
                    "raw_line": p["raw_line"],
                    "raw_header": p["raw_header"],
                    "source_uri": p["source_uri"],
                    "source_sha256": p["source_sha256"],
                    "import_code_sha": import_code_sha,
                    "members_json": json.dumps(matching, sort_keys=True),
                }
            )
            gpi_index.append(
                {
                    "prediction_key": list(key),
                    "input_number": p["input_number"],
                    "omega_position_text": p["omega_position_text"],
                    "members": matching,
                    "submission_source": {k: v for k, v in provenance.items() if k != "receipt"},
                    "submitted_record_provenance": receipt_record,
                }
            )
    for relative in manifest.get("extra_source_files", []):
        source(local_path(root, relative))
    out.mkdir(parents=True, exist_ok=False)
    (out / "topology_records.jsonl").write_text(
        "".join(json.dumps(row, sort_keys=True) + "\n" for _, row in sorted(topology_rows.items()))
    )
    (out / "gpi_records.json").write_text(
        json.dumps({"predictions": gpi_predictions, "curation": []}, indent=2) + "\n"
    )
    index = {
        "topology": [v for _, v in sorted(topology_index.items())],
        "netgpi": gpi_index,
        "validated_topology_runs": validation,
        "sources": sources,
        "export_manifest_sha256": sha(manifest_path),
        "export_manifest": manifest,
        "import_code_sha": import_code_sha,
        "legacy_curation": "not read or modified",
    }
    (out / "prediction_provenance_index.json").write_text(json.dumps(index, indent=2) + "\n")
    summary = {
        "topology_version": TOPOLOGY_VERSION,
        "topology_rows": len(topology_rows),
        "netgpi_predictions": len(gpi_predictions),
        "curation_records": 0,
        "topology_tool_versions": sorted({r["tool_version"] for r in topology_rows.values()}),
        "netgpi_tool_versions": sorted({r["tool_version"] for r in gpi_predictions}),
        "source_files": sources,
        "activation": False,
        "output_sha256": {p.name: sha(p) for p in sorted(out.iterdir())},
    }
    (out / "prediction_export_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("root", "manifest", "out"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--import-code-sha", required=True)
    args = parser.parse_args()
    summary = export(args.root, args.manifest, args.out, args.import_code_sha)
    print(json.dumps({k: v for k, v in summary.items() if k != "source_files"}, indent=2))


if __name__ == "__main__":
    main()
