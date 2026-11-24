"""Prepare or explicitly execute an approved native-ortholog NetGPI1.1 review batch.

The existing engineered wrapper supplies the pinned DTU helper/endpoint checks;
its 61-record constants and outputs are untouched. Preparation is offline and
binds actual FASTA bytes/count to the native CDS-eligible scoring picks and ledger.
Future replacement batches use fresh directories and the same checks.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
import subprocess
from datetime import UTC, datetime
from decimal import Decimal, InvalidOperation
from pathlib import Path

import run_engineered_netgpi as common

APPROVAL = "Approve NetGPI checks through this rebuild"
AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")
PICK_FIELDS = (
    "gene",
    "species",
    "compara_release",
    "human_ensembl_gene",
    "ortholog_ensembl_gene",
    "picked_acc",
    "picked_source",
    "picked_cds_transcript",
    "picked_cds_sha256",
    "picked_topology_label",
    "topology_tool",
    "selection_state",
    "human_reference_sha256",
)


def read_table(path: Path) -> list[dict]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def source_info() -> dict:
    root = Path(__file__).resolve().parent
    commit = subprocess.run(
        ["git", "-C", str(root), "rev-parse", "HEAD"], check=True, capture_output=True, text=True
    ).stdout.strip()
    dirty = subprocess.run(
        ["git", "-C", str(root), "status", "--porcelain"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    return {"commit": commit, "dirty": bool(dirty)}


def batch_records(fasta: Path, picks: Path, ledger: Path, helper) -> list[dict]:
    sequences = helper.records(fasta)
    if not sequences or len({name for name, _ in sequences}) != len(sequences):
        raise ValueError("Expected nonempty uniquely named FASTA records")
    for name, sequence in sequences:
        if name != common.sha(sequence.encode()) or set(sequence) - AMINO_ACIDS:
            raise ValueError(f"FASTA ID/alphabet differs from exact native protein: {name}")
    hashes = {name for name, _ in sequences}
    selected = [r for r in read_table(picks) if r.get("picked_sequence_sha256") in hashes]
    evidence = [r for r in read_table(ledger) if r.get("protein_sha256") in hashes]
    records = []
    for ordinal, (name, sequence) in enumerate(sequences, 1):
        matches = [r for r in selected if r["picked_sequence_sha256"] == name]
        if not matches:
            raise ValueError(f"Submitted protein is not a provisional native ortholog pick: {name}")
        provenance = []
        for pick in sorted(matches, key=lambda r: (r["species"], r["gene"], r["picked_acc"])):
            if (
                pick["species"] not in {"mouse", "cyno"}
                or pick.get("compara_release") != "116"
                or pick.get("human_is_gpi") != "True"
                or pick.get("selection_state") not in {"needs_netgpi", "needs_netgpi_and_v1"}
                or int(pick["picked_length"]) != len(sequence)
            ):
                raise ValueError(f"Pick is not a pending native ortholog anchor check: {name}")
            if not pick.get("picked_cds_transcript") or not re.fullmatch(
                r"[0-9a-f]{64}", pick.get("picked_cds_sha256", "")
            ):
                raise ValueError(f"Pick lacks exact CDS provenance: {name}")
            aliases = [
                r
                for r in evidence
                if r["species"] == pick["species"]
                and r["gene"] == pick["gene"]
                and r["ortholog_ensembl_gene"] == pick["ortholog_ensembl_gene"]
            ]
            winners = [
                r
                for r in aliases
                if r["candidate_acc"] == pick["picked_acc"]
                and r["cds_eligibility"] == "eligible"
                and r["winner_alias"] == "True"
                and r["winner_protein"] == "True"
            ]
            if (
                len(winners) != 1
                or winners[0]["representative_transcript"] != pick["picked_cds_transcript"]
            ):
                raise ValueError(f"Scoring ledger does not prove the exact eligible winner: {name}")
            provenance.append(
                {
                    "pick": {key: pick.get(key, "") for key in PICK_FIELDS},
                    "winner_ledger": winners[0],
                    "same_protein_aliases": aliases,
                    "member_identity_status": "provisional ortholog; no new member ID allocated by this batch",
                    "recorded_member_references": {
                        key: value
                        for key, value in pick.items()
                        if "candidate_id" in key or "member_id" in key
                    },
                }
            )
        records.append(
            {
                "input_number": ordinal,
                "record_id": name,
                "sequence_sha256": name,
                "sequence_length": len(sequence),
                "native_ortholog_provenance": provenance,
            }
        )
    return records


def fingerprints(fasta: Path, picks: Path, ledger: Path, helper_path: Path) -> dict:
    return {
        "input_sha256": common.sha(fasta.read_bytes()),
        "picks_sha256": common.sha(picks.read_bytes()),
        "ledger_sha256": common.sha(ledger.read_bytes()),
        "helper_sha256": common.sha(helper_path.read_bytes()),
        "recording_helper_sha256": common.sha(Path(common.__file__).read_bytes()),
        "wrapper_sha256": common.sha(Path(__file__).read_bytes()),
    }


def prepare(
    fasta: Path, picks: Path, ledger: Path, out: Path, helper_path: Path = common.HELPER
) -> dict:
    helper = common.helper_module(helper_path)
    records = batch_records(fasta, picks, ledger, helper)
    pins = fingerprints(fasta, picks, ledger, helper_path)
    metadata = {
        "status": "prepared_for_explicit_execution",
        "job_id": None,
        "run_id": "native-ortholog-netgpi11-" + pins["input_sha256"][:12],
        "tool_version": "NetGPI-1.1",
        "version_provenance": "pinned DTU service config and result URL",
        "format": "short",
        "sequence_count": len(records),
        "records": records,
        **pins,
        "service": helper.CGI,
        "config": helper.CONFIG,
        "result_url_template": helper.RESULT,
        "approval_reference": APPROVAL,
        "approval_scope": "public native ortholog replacement candidates through this rebuild",
        "execute_approved": False,
        "preparation_source": source_info(),
        "source_files": {"fasta": fasta.name, "picks": picks.name, "ledger": ledger.name},
        "prediction_number_definition": "one-based data-row ordinal in this DTU job's short output; input_number is recorded separately",
    }
    out.mkdir(parents=True, exist_ok=False)
    (out / "submitted.fasta").write_bytes(fasta.read_bytes())
    (out / "submission.json").write_text(json.dumps(metadata, indent=2) + "\n")
    return metadata


def validate_prepared(
    fasta: Path, picks: Path, ledger: Path, out: Path, helper_path: Path
) -> tuple[dict, object]:
    metadata = json.loads((out / "submission.json").read_text())
    if metadata["status"] != "prepared_for_explicit_execution" or metadata.get("job_id"):
        raise ValueError(
            "Prepared batch is already submitted or not ready; refuses duplicate submission"
        )
    helper = common.helper_module(helper_path)
    for key, actual in fingerprints(fasta, picks, ledger, helper_path).items():
        if metadata.get(key) != actual:
            raise ValueError(f"Prepared batch provenance changed: {key}")
    if common.sha((out / "submitted.fasta").read_bytes()) != metadata["input_sha256"]:
        raise ValueError("Prepared FASTA copy changed")
    if metadata["records"] != batch_records(fasta, picks, ledger, helper):
        raise ValueError("Prepared scoring/sequence records changed")
    return metadata, helper


def recorded_curl(out: Path, runner):
    """Same response-capture pattern as run_engineered_netgpi, with GET headers too."""
    http = out / "http"
    http.mkdir(exist_ok=False)
    exchanges = []

    def curl(arguments):
        number = len(exchanges) + 1
        stem = http / f"{number:04d}"
        command = ["curl", "-sS", *arguments]
        if "-o" in command:
            command[command.index("-o") + 1] = str(stem.with_suffix(".body.bin"))
        if "-D" not in command:
            command.extend(["-D", str(stem.with_suffix(".headers.bin"))])
        result = runner(command, capture_output=True, check=False)
        stem.with_suffix(".stdout.bin").write_bytes(result.stdout)
        stem.with_suffix(".stderr.bin").write_bytes(result.stderr)
        if "-D" in arguments and arguments[arguments.index("-D") + 1] == "-":
            stem.with_suffix(".headers.bin").write_bytes(result.stdout)
        exchange = {"number": number, "argv": command, "returncode": result.returncode}
        if "-w" in arguments and result.stdout[-3:].isdigit():
            exchange["http_status_text"] = result.stdout[-3:].decode()
            stem.with_suffix(".body.bin").write_bytes(result.stdout[:-3])
            if result.stdout[-3:] == b"200":
                (out / "raw_output_protein_type.txt").write_bytes(result.stdout[:-3])
        exchanges.append(exchange)
        (out / "http_exchanges.json").write_text(json.dumps(exchanges, indent=2) + "\n")
        return (
            result.returncode,
            result.stdout.decode(errors="replace"),
            result.stderr.decode(errors="replace"),
        )

    return curl


def parse_predictions(rows: list[str], metadata: dict) -> list[dict]:
    inputs = {r["record_id"]: r for r in metadata["records"]}
    if len(rows) != len(inputs):
        raise ValueError("NetGPI output count differs from submitted sequence count")
    predictions, seen = [], set()
    for number, line in enumerate(rows, 1):
        fields = line.split("\t")
        if len(fields) != 6 or fields[0] not in inputs or fields[0] in seen:
            raise ValueError("Malformed, unknown or duplicate NetGPI output ID")
        name = fields[0]
        if int(fields[1]) != inputs[name]["sequence_length"]:
            raise ValueError("NetGPI output length differs from exact submitted native protein")
        try:
            likelihood = Decimal(fields[4])
        except InvalidOperation as exc:
            raise ValueError("Malformed literal NetGPI likelihood") from exc
        if not likelihood.is_finite() or not 0 <= likelihood <= 1:
            raise ValueError("NetGPI likelihood outside [0,1]")
        seen.add(name)
        predictions.append(
            {
                "run_id": metadata["run_id"],
                "job_id": metadata["job_id"],
                "tool_version": "NetGPI-1.1",
                "prediction_number": number,
                "input_number": inputs[name]["input_number"],
                "record_id": name,
                "sequence_sha256": name,
                "sequence_length": int(fields[1]),
                "prediction": fields[2],
                "omega_position_text": fields[3],
                "likelihood_text": fields[4],
                "omega_residue": fields[5],
                "raw_line": line,
            }
        )
    return predictions


def execute(
    metadata: dict, helper, out: Path, *, approved: bool = False, runner=subprocess.run
) -> dict:
    if not approved:
        raise ValueError("Explicit approved execution flag is required")
    metadata.update(
        status="submitting_approved_batch",
        execute_approved=True,
        execution_source=source_info(),
        execution_started_utc=datetime.now(UTC).isoformat(),
    )
    path = out / "submission.json"

    def save():
        path.write_text(json.dumps(metadata, indent=2) + "\n")

    save()
    helper.curl = recorded_curl(out, runner)
    try:
        job = helper.submit(out / "submitted.fasta")
        metadata.update(job_id=job, result_url=helper.RESULT.format(job=job), status="submitted")
        save()
        print(
            f"Submitted {metadata['sequence_count']} native proteins; NetGPI1.1 job {job}",
            flush=True,
        )
        rows = helper.collect(job, metadata["sequence_count"])
        raw_path = out / "raw_output_protein_type.txt"
        raw = raw_path.read_text()
        reported = sorted(
            set(re.findall(r"NetGPI[\s-]+(?:version\s+)?(\d+\.\d+)", raw, flags=re.I))
        )
        if any(version != "1.1" for version in reported):
            raise ValueError("Raw output reports a different NetGPI version")
        metadata["reported_tool_versions"] = reported
        predictions = parse_predictions(rows, metadata)
        (out / "netgpi_predictions.txt").write_text(
            "# NetGPI 1.1 short output\n" + "\n".join(rows) + "\n"
        )
        (out / "predictions.json").write_text(json.dumps(predictions, indent=2) + "\n")
        metadata.update(status="completed", prediction_count=len(predictions))
    except BaseException as exc:
        metadata.update(
            status="failed_after_submission_attempt", error_type=type(exc).__name__, error=str(exc)
        )
        save()
        raise
    metadata["output_sha256"] = {
        str(p.relative_to(out)): common.sha(p.read_bytes())
        for p in sorted(out.rglob("*"))
        if p.is_file() and p != path
    }
    save()
    return metadata


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("fasta", "picks", "ledger", "out"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--helper", type=Path, default=common.HELPER)
    parser.add_argument("--execute-approved", action="store_true")
    args = parser.parse_args(argv)
    if args.execute_approved:
        metadata, helper = validate_prepared(
            args.fasta, args.picks, args.ledger, args.out, args.helper
        )
        execute(metadata, helper, args.out, approved=True)
    else:
        metadata = prepare(args.fasta, args.picks, args.ledger, args.out, args.helper)
        print(
            f"Prepared {metadata['sequence_count']} native proteins offline: {args.out / 'submission.json'}"
        )


if __name__ == "__main__":
    main()
