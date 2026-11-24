"""Export exact prediction records and independent GPI curation for D1 import.

The run manifest lists each submitted FASTA, raw output and external job ID.
Historical NetGPI outputs and curation are archived separately. Missing legacy
input sequences and job IDs stay unknown: an accession and a length are not an
exact-sequence provenance record.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import defaultdict
from pathlib import Path

from publish_topology_supplement import derive


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fasta(path):
    result = {}
    name, sequence = None, []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                result[name] = "".join(sequence)
            name, sequence = line[1:].strip(), []
        elif line.strip():
            sequence.append(line.strip())
    if name is not None:
        result[name] = "".join(sequence)
    return result


def export(root, manifest_path, out, import_code_sha):
    manifest = json.loads(manifest_path.read_text())
    prefix = manifest["s3_prefix"]
    members = defaultdict(dict)
    for name in manifest["pick_tables"]:
        for row in csv.DictReader((root / name).open(), delimiter="\t"):
            member = {
                k: row[k]
                for k in (
                    "gene",
                    "species",
                    "human_ensembl_gene",
                    "ortholog_ensembl_gene",
                    "ortholog_symbols",
                    "picked_acc",
                )
            }
            members[row["picked_sequence_sha256"]][json.dumps(member, sort_keys=True)] = member
    source_files = {}

    def source(path):
        relative = str(path.relative_to(root))
        source_files[relative] = {"sha256": sha(path), "source": "processed:" + prefix + relative}
        return source_files[relative]["source"]

    predictions = []
    for run in manifest["netgpi"]:
        path = root / run["output"]
        inputs = fasta(root / run["fasta"]) if run.get("fasta") else {}
        if run.get("fasta"):
            source(root / run["fasta"])
        text = path.read_text()
        header = "\n".join(line for line in text.splitlines() if line.startswith("#"))
        version = re.search(r"NetGPI\s+([\d.]+)", header)
        if version is None:
            raise ValueError(f"{path.name}: no recorded NetGPI version")
        seen = set()
        number = 0
        for line in text.splitlines():
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) != 6 or fields[0] in seen:
                raise ValueError(f"{path.name}: malformed or duplicate NetGPI output ID")
            seen.add(fields[0])
            number += 1
            sequence = inputs.get(fields[0])
            if inputs and sequence is None:
                raise ValueError("prediction ID absent from the actual submitted FASTA")
            if sequence is not None and len(sequence) != int(fields[1]):
                raise ValueError("predicted length differs from submitted sequence")
            digest = hashlib.sha256(sequence.encode()).hexdigest() if sequence is not None else None
            if run.get("digest_ids") and digest != fields[0]:
                raise ValueError("submitted digest ID does not match its sequence")
            predictions.append(
                {
                    "run_id": run["run_id"],
                    "record_id": fields[0],
                    "prediction_number": number,
                    "tool_version": "NetGPI-" + version.group(1),
                    "job_id": run.get("job_id"),
                    "sequence_sha256": digest,
                    "sequence": sequence,
                    "sequence_resolution": "exact_input"
                    if sequence is not None
                    else "unresolved_legacy",
                    "sequence_length": int(fields[1]),
                    "prediction": fields[2],
                    "likelihood_text": fields[4],
                    "omega_position": int(fields[3]) if fields[3].isdigit() else None,
                    "omega_residue": fields[5],
                    "raw_line": line,
                    "raw_header": header,
                    "source_uri": source(path),
                    "source_sha256": sha(path),
                    "import_code_sha": import_code_sha,
                    "members_json": json.dumps(list(members[digest].values()), sort_keys=True),
                }
            )
        if inputs and seen != set(inputs):
            raise ValueError("NetGPI output does not cover every submitted sequence")
    curation = []
    for entry in manifest["curation"]:
        path = root / entry["path"]
        for i, row in enumerate(csv.DictReader(path.open(), delimiter="\t"), 1):
            if any(row.get(key) != value for key, value in entry.get("where", {}).items()):
                continue
            curation.append(
                {
                    "source_sha256": sha(path),
                    "record_number": i,
                    "source_kind": entry["kind"],
                    "gene_symbol": row.get("gene_symbol")
                    or row.get("gene")
                    or row.get("human_symbol")
                    or row.get("Gene Names (primary)"),
                    "source_uri": source(path),
                    "import_code_sha": import_code_sha,
                    "record_json": json.dumps(row, sort_keys=True),
                }
            )
    topology = {}
    for run in manifest["topology"]:
        path = root / run["output"]
        submitted = fasta(root / run["fasta"])
        source(path)
        source(root / run["fasta"])
        lines = [line.strip() for line in path.read_text().splitlines() if line.strip()]
        if len(lines) % 3:
            raise ValueError("incomplete DeepTMHMM record")
        seen = set()
        for i in range(0, len(lines), 3):
            identity, label = lines[i][1:].rsplit("|", 1)
            digest, label = identity.strip(), label.strip()
            sequence, topo = lines[i + 1 : i + 3]
            if (
                submitted.get(digest) != sequence
                or hashlib.sha256(sequence.encode()).hexdigest() != digest
            ):
                raise ValueError("DeepTMHMM result differs from its submitted protein")
            if len(topo) != len(sequence):
                raise ValueError("DeepTMHMM topology length differs")
            seen.add(digest)
            if not members[digest]:
                raise ValueError("predicted protein lacks its source ortholog identity")
            for member in members[digest].values():
                row = {
                    "topology_version": manifest["topology_version"],
                    "cohort": member["species"] + "_ortholog_repick",
                    "hgnc_id": None,
                    "uniprot_accession": member["picked_acc"].split("-")[0],
                    "uniprot_accession_full": member["picked_acc"],
                    "isoform_id": member["picked_acc"],
                    "gene_symbol": member["ortholog_symbols"],
                    "species": {"cyno": "cynomolgus", "mouse": "mouse"}[member["species"]],
                    "is_canonical": 0,
                    "sequence": sequence,
                    "protein_length": len(sequence),
                    "deeptmhmm_label": label,
                    "per_residue_topology": topo,
                    "predicted_surface_membrane": int(label in {"TM", "SP+TM", "BETA"}),
                    "predicted_secreted": int(label == "SP"),
                    "tool_version": "deeptmhmm-1.0.24",
                    "retrieved_at": manifest["date"],
                    **derive(topo),
                }
                key = (row["cohort"], row["uniprot_accession_full"])
                if key in topology and topology[key]["sequence"] != sequence:
                    raise ValueError(f"one prediction key refers to different sequences: {key}")
                topology[key] = row
        if seen != set(submitted):
            raise ValueError("DeepTMHMM output does not cover its entire submitted FASTA")
    out.mkdir(parents=True, exist_ok=True)
    (out / "gpi_records.json").write_text(
        json.dumps({"predictions": predictions, "curation": curation}, indent=1) + "\n"
    )
    (out / "topology_records.jsonl").write_text(
        "".join(json.dumps(row, sort_keys=True) + "\n" for row in topology.values())
    )
    summary = {
        "netgpi_predictions": len(predictions),
        "curation_records": len(curation),
        "netgpi_exact_input_sequences": sum(r["sequence"] is not None for r in predictions),
        "netgpi_unresolved_legacy_sequences": sum(r["sequence"] is None for r in predictions),
        "topology_rows": len(topology),
        "topology_version": manifest["topology_version"],
        "sources": source_files,
    }
    (out / "prediction_export_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return {k: v for k, v in summary.items() if k != "sources"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--import-code-sha", required=True)
    args = parser.parse_args()
    print(json.dumps(export(args.root, args.manifest, args.out, args.import_code_sha), indent=2))


if __name__ == "__main__":
    main()
