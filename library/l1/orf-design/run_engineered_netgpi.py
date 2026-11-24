"""Prepare, or explicitly execute an approved, sequence-pinned NetGPI 1.1 batch.

Preparation is offline. Execution delegates to the existing DTU helper while
retaining the exact upload, HTTP responses, job ID and numbered predictions.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
HELPER = HERE.parents[2] / "library/l1/proteoform-selection/netgpi_submit.py"
EXPECTED_FASTA_SHA256 = "b5961bbb973bdada02d299d2c73cac6925eed4f56d4d04078508ef0f270d1c6d"
EXPECTED_COUNT = 61


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def helper_module(path: Path):
    spec = importlib.util.spec_from_file_location("netgpi_helper", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    if module.CGI != "https://services.healthtech.dtu.dk/cgi-bin/webface2.cgi":
        raise ValueError("Unexpected submission service")
    if module.CONFIG != "/var/www/services/services/NetGPI-1.1/webface.cf":
        raise ValueError("Unexpected NetGPI service version")
    return module


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--variants", type=Path, required=True)
    parser.add_argument("--helper", type=Path, default=HELPER)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument(
        "--execute-approved", action="store_true", help="submit only after explicit approval"
    )
    args = parser.parse_args(argv)
    data = args.fasta.read_bytes()
    if sha(data) != EXPECTED_FASTA_SHA256:
        raise ValueError("FASTA differs from the prepared 61-sequence payload")
    helper = helper_module(args.helper)
    recs = helper.records(args.fasta)
    seqs = dict(recs)
    if len(recs) != EXPECTED_COUNT or len(seqs) != EXPECTED_COUNT:
        raise ValueError("Expected 61 uniquely named sequences")
    with args.variants.open() as handle:
        variants = {r["prediction_id"]: r for r in csv.DictReader(handle, delimiter="\t")}
    if set(variants) != set(seqs):
        raise ValueError("Prediction IDs differ between FASTA and design records")
    records = []
    for ordinal, (name, sequence) in enumerate(recs, 1):
        digest = sha(sequence.encode())
        if digest != variants[name]["designed_protein_sha256"]:
            raise ValueError(f"Engineered protein hash differs: {name}")
        records.append(
            {
                "input_number": ordinal,
                "record_id": name,
                "sequence_sha256": digest,
                "sequence_length": len(sequence),
                "members": variants[name]["affected_member_ids"],
                "native_protein_sha256": variants[name]["native_protein_sha256"],
            }
        )
    metadata_path = args.out / "submission.json"
    if not args.execute_approved:
        args.out.mkdir(parents=True, exist_ok=False)
        (args.out / "submitted.fasta").write_bytes(data)
        metadata = {
            "status": "prepared_awaiting_approval",
            "job_id": None,
            "run_id": "engineered-cds-netgpi11-" + EXPECTED_FASTA_SHA256[:12],
            "tool_version": "NetGPI-1.1",
            "format": "short",
            "sequence_count": len(recs),
            "input_sha256": EXPECTED_FASTA_SHA256,
            "variants_sha256": sha(args.variants.read_bytes()),
            "helper_sha256": sha(args.helper.read_bytes()),
            "wrapper_sha256": sha(Path(__file__).read_bytes()),
            "service": helper.CGI,
            "config": helper.CONFIG,
            "result_url_template": helper.RESULT,
            "records": records,
        }
        metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
        print(f"Prepared {len(recs)} sequences; no submission made: {metadata_path}")
        return
    metadata = json.loads(metadata_path.read_text())
    if metadata["status"] != "prepared_awaiting_approval" or metadata.get("job_id"):
        raise ValueError("This prepared batch was already submitted or is not ready")
    for key, actual in (
        ("input_sha256", sha((args.out / "submitted.fasta").read_bytes())),
        ("variants_sha256", sha(args.variants.read_bytes())),
        ("helper_sha256", sha(args.helper.read_bytes())),
        ("wrapper_sha256", sha(Path(__file__).read_bytes())),
    ):
        if metadata[key] != actual:
            raise ValueError(f"Prepared provenance changed: {key}")
    http = args.out / "http"
    http.mkdir(exist_ok=False)
    exchanges = []

    def recorded_curl(arguments):
        number = len(exchanges) + 1
        stem = http / f"{number:04d}"
        command = ["curl", "-sS", *arguments]
        # Keep the submission response body that the original helper discards.
        if "-o" in command:
            command[command.index("-o") + 1] = str(stem.with_suffix(".body.bin"))
        result = subprocess.run(command, capture_output=True, check=False)
        stem.with_suffix(".stdout.bin").write_bytes(result.stdout)
        stem.with_suffix(".stderr.bin").write_bytes(result.stderr)
        exchanges.append({"number": number, "argv": command, "returncode": result.returncode})
        if "-w" in arguments and result.stdout[-3:].isdigit():
            stem.with_suffix(".body.bin").write_bytes(result.stdout[:-3])
            if result.stdout[-3:] == b"200":
                (args.out / "raw_output_protein_type.txt").write_bytes(result.stdout[:-3])
        (args.out / "http_exchanges.json").write_text(json.dumps(exchanges, indent=2) + "\n")
        return result.returncode, result.stdout.decode(), result.stderr.decode()

    helper.curl = recorded_curl
    metadata["status"] = "submitting_approved_batch"
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
    job = helper.submit(args.out / "submitted.fasta")
    metadata.update(
        {"job_id": job, "result_url": helper.RESULT.format(job=job), "status": "submitted"}
    )
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
    print(f"Submitted {len(recs)} sequences, NetGPI 1.1 job {job}", flush=True)
    rows = helper.collect(job, len(recs))
    predictions, seen = [], set()
    for number, line in enumerate(rows, 1):
        fields = line.split("\t")
        if len(fields) != 6 or fields[0] not in seqs or fields[0] in seen:
            raise ValueError("Malformed, unknown or duplicate NetGPI output record")
        seen.add(fields[0])
        if int(fields[1]) != len(seqs[fields[0]]):
            raise ValueError("NetGPI length differs from the submitted engineered sequence")
        predictions.append(
            {
                "run_id": metadata["run_id"],
                "job_id": job,
                "tool_version": "NetGPI-1.1",
                "prediction_number": number,
                "record_id": fields[0],
                "sequence_sha256": sha(seqs[fields[0]].encode()),
                "sequence_length": int(fields[1]),
                "prediction": fields[2],
                "omega_position_text": fields[3],
                "likelihood_text": fields[4],
                "omega_residue": fields[5],
                "raw_line": line,
            }
        )
    if seen != set(seqs):
        raise ValueError("NetGPI predictions do not cover exactly the 61 inputs")
    (args.out / "netgpi_predictions.txt").write_text(
        "# NetGPI 1.1 short output\n" + "\n".join(rows) + "\n"
    )
    (args.out / "predictions.json").write_text(json.dumps(predictions, indent=2) + "\n")
    metadata.update({"status": "completed", "prediction_count": len(predictions)})
    metadata["output_sha256"] = {
        str(p.relative_to(args.out)): sha(p.read_bytes())
        for p in sorted(args.out.rglob("*"))
        if p.is_file() and p != metadata_path
    }
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
    print(f"Recorded {len(predictions)} predictions with exact inputs, job ID and raw responses")


if __name__ == "__main__":
    main()
