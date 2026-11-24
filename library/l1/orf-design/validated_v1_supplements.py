"""Validate explicit DeepTMHMMv1 run provenance without altering a snapshot.

Topology alphabets alone cannot identify the predictor. Every supplement is
bound to submitted sequences, runtime model identity and a committed v1 runner.
"""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path

import score_ortholog_pool as legacy

SCHEMA = "ortholog-v1-run-1"
TOOL_VERSION = "deeptmhmm-1.0-academic"
RUNNER_PATH = "library/l1/proteoform-selection/modal/deeptmhmm_v1_app.py"
LABELS = {"TM", "SP+TM", "SP", "GLOB", "BETA"}
ROOT = Path(__file__).resolve().parents[3]


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def bound_file(record_path, spec):
    path = Path(spec["path"])
    if path.is_absolute():
        raise ValueError("Run-record file paths must be relative to its directory")
    path = (record_path.parent / path).resolve()
    if (
        not re.fullmatch(r"[0-9a-f]{64}", spec.get("sha256", ""))
        or sha256(path.read_bytes()) != spec["sha256"]
    ):
        raise ValueError(f"Run-record file checksum mismatch: {path.name}")
    return path


def committed_runner_blob(commit):
    from deliverome_analysis.provenance import read_pinned_bytes

    if not re.fullmatch(r"[0-9a-f]{40}", commit):
        raise ValueError("Runner provenance requires its full committed Git revision")
    try:
        return read_pinned_bytes(commit, RUNNER_PATH, repo=ROOT)
    except FileNotFoundError as error:
        raise ValueError("Committed v1 runner blob is unavailable") from error


def submitted_fasta(path):
    records, name, parts = {}, None, []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                if name in records:
                    raise ValueError("Duplicate submitted FASTA digest")
                records[name] = "".join(parts)
            name, parts = line[1:].strip(), []
        elif line.strip():
            if name is None:
                raise ValueError("Submitted FASTA sequence has no digest header")
            parts.append(line.strip())
    if name is not None:
        if name in records:
            raise ValueError("Duplicate submitted FASTA digest")
        records[name] = "".join(parts)
    if not records or any(sha256(seq.encode()) != name for name, seq in records.items()):
        raise ValueError("Submitted FASTA headers must be exact full sequence digests")
    return records


def validate_runtime(identity, runner, record):
    fingerprints = identity.get("package_files_sha256", {})
    models = [f"deeptmhmm_cv_{i}.model" for i in range(5)]
    checkpoints = identity.get("esm_checkpoint_files", [])
    if (
        identity.get("schema_version") != "deeptmhmm-academic-sweep-1"
        or identity.get("version_source")
        != 'predict.py ArgumentParser description="DeepTMHMM v1.0"'
        or identity.get("package_predict_sha256") != fingerprints.get("predict.py")
        or identity.get("model_files") != models
        or not checkpoints
        or any(not path.endswith(".pt") for path in checkpoints)
        or any(path not in fingerprints for path in models + checkpoints)
        or any(not re.fullmatch(r"[0-9a-f]{64}", digest) for digest in fingerprints.values())
        or identity.get("submitted_fasta_sha256") != record["submitted_fasta"]["sha256"]
        or identity.get("output_sha256") != record["output"]["sha256"]
        or identity.get("failed_shards") != 0
        or not identity.get("shards")
    ):
        raise ValueError("Runtime identity lacks exact package weights and sweep provenance")
    common = {
        key: identity[key]
        for key in (
            "tool_version",
            "model_version",
            "version_source",
            "runner_sha256",
            "package_predict_sha256",
            "package_files_sha256",
            "model_files",
            "esm_checkpoint_files",
        )
    }
    for shard in identity["shards"]:
        if (
            shard.get("ok") is not True
            or shard.get("runtime_identity") != common
            or not shard.get("tag")
            or shard.get("n") != len(shard.get("sequence_sha256s", []))
            or not re.fullmatch(r"[0-9a-f]{64}", shard.get("submitted_fasta_sha256", ""))
            or not re.fullmatch(r"[0-9a-f]{64}", shard.get("output_sha256", ""))
        ):
            raise ValueError("Runtime shard lacks successful, consistent model provenance")


def load(supplements, run_records):
    if len(supplements) != len(run_records):
        raise ValueError("Each v1 supplement requires one corresponding v1 run record")
    by_species = {"mouse": {}, "cyno": {}}
    runs, duplicate_reports = [], []
    for output_path, record_path in zip(supplements, run_records, strict=True):
        record = json.loads(record_path.read_text())
        if record.get("schema_version") != SCHEMA or record.get("tool_version") != TOOL_VERSION:
            raise ValueError("Supplement lacks supported DeepTMHMMv1 run provenance")
        if (
            not isinstance(record.get("run_id"), str)
            or not record["run_id"]
            or not isinstance(record.get("job_id"), str)
            or not record["job_id"]
        ):
            raise ValueError("Supplement provenance requires run and external job IDs")
        runner = record["runner"]
        if runner.get("path") != RUNNER_PATH or sha256(
            committed_runner_blob(runner["git_commit"])
        ) != runner.get("sha256"):
            raise ValueError("Supplement did not originate from its committed v1 runner")
        fasta_path = bound_file(record_path, record["submitted_fasta"])
        bound_output = bound_file(record_path, record["output"])
        if bound_output != output_path.resolve():
            raise ValueError("v1 supplement/run-record output pairing mismatch")
        identity_path = bound_file(record_path, record["runtime_identity"])
        identity = json.loads(identity_path.read_text())
        if (
            identity.get("tool_version") != TOOL_VERSION
            or identity.get("model_version") != "1.0"
            or not isinstance(identity.get("version_source"), str)
            or not identity["version_source"]
            or identity.get("runner_sha256") != runner["sha256"]
            or not re.fullmatch(r"[0-9a-f]{64}", identity.get("package_predict_sha256", ""))
        ):
            raise ValueError(
                "Runtime model identity does not prove the expected v1 package/version"
            )
        validate_runtime(identity, runner, record)
        submitted = submitted_fasta(fasta_path)
        shard_sequences = [sha for shard in identity["shards"] for sha in shard["sequence_sha256s"]]
        if len(shard_sequences) != len(set(shard_sequences)) or set(shard_sequences) != set(
            submitted
        ):
            raise ValueError("Runtime shards do not exactly cover submitted sequences")
        labels = legacy.read_v1_supplement(output_path)
        if set(labels) != set(submitted) or set(labels.values()) - LABELS:
            raise ValueError("v1 output does not exactly cover its submitted FASTA/classes")
        raw_lines = output_path.read_text().splitlines(keepends=True)
        if len(raw_lines) != 3 * len(submitted) or any(not line.strip() for line in raw_lines):
            raise ValueError("Supplement must preserve the runner's exact three-line records")
        raw_records = {}
        for index in range(0, len(raw_lines), 3):
            digest = raw_lines[index][1:].rsplit("|", 1)[0].strip()
            raw_records[digest] = "".join(raw_lines[index : index + 3])
        for shard in identity["shards"]:
            digests = shard["sequence_sha256s"]
            fasta_text = "".join(f">{digest}\n{submitted[digest]}\n" for digest in digests)
            output_text = "".join(raw_records[digest] for digest in digests)
            if (
                sha256(fasta_text.encode()) != shard["submitted_fasta_sha256"]
                or sha256(output_text.encode()) != shard["output_sha256"]
            ):
                raise ValueError("Runtime shard hashes disagree with submitted/output records")
        members = {sha: set() for sha in submitted}
        for member in record["members"]:
            sha, species = member["sequence_sha256"], member["species"]
            if sha not in submitted or species not in by_species:
                raise ValueError(
                    "v1 run membership has an unsubmitted sequence or unsupported species"
                )
            members[sha].add(species)
        if any(not species for species in members.values()):
            raise ValueError("Every v1 submitted sequence requires explicit species membership")
        lines = [line.strip() for line in output_path.read_text().splitlines() if line.strip()]
        for i in range(0, len(lines), 3):
            sha = lines[i][1:].rsplit("|", 1)[0].strip()
            sequence, topology = lines[i + 1 : i + 3]
            if sequence != submitted[sha]:
                raise ValueError("v1 output protein differs from submitted input")
            evidence = {
                "run_id": record["run_id"],
                "job_id": record["job_id"],
                "tool_version": TOOL_VERSION,
                "output_record_number": i // 3 + 1,
                "runner": runner,
                "run_record_sha256": sha256(record_path.read_bytes()),
                "submitted_fasta_sha256": record["submitted_fasta"]["sha256"],
                "output_sha256": record["output"]["sha256"],
                "runtime_identity_sha256": record["runtime_identity"]["sha256"],
                "model_identity": identity,
            }
            for species in sorted(members[sha]):
                previous = by_species[species].get(sha)
                if previous and previous["label"] != labels[sha]:
                    raise ValueError(f"Conflicting v1 supplement labels for {species}/{sha}")
                if previous:
                    duplicate_reports.append(
                        {
                            "species": species,
                            "sequence_sha256": sha,
                            "state": "duplicate_consistent"
                            if previous["topology"] == topology
                            else "same_label_different_topology",
                            "run_id": record["run_id"],
                        }
                    )
                    previous["provenance"].append(evidence)
                else:
                    by_species[species][sha] = {
                        "label": labels[sha],
                        "topology": topology,
                        "sequence": sequence,
                        "provenance": [evidence],
                    }
        runs.append(
            {
                "run_id": record["run_id"],
                "job_id": record["job_id"],
                "tool_version": TOOL_VERSION,
                "output": record["output"],
                "submitted_fasta": record["submitted_fasta"],
                "runtime_identity": record["runtime_identity"],
                "runner": runner,
                "run_record_sha256": sha256(record_path.read_bytes()),
                "proteins": len(submitted),
            }
        )
    return by_species, {"runs": runs, "duplicate_reports": duplicate_reports}


def overlay(species, report_rows, supplements):
    """Use validated new labels only when the frozen same-species v1 call is absent."""
    labels, evidence, reconciliation = {}, {}, []
    for row in report_rows:
        sha = row["key"]
        label = row.get("deeptmhmm_label", "") if row["found"] else ""
        new = supplements.get(species, {}).get(sha)
        if new and label and label != new["label"]:
            raise ValueError(
                f"Conflicting frozen/new v1 labels for {species}/{sha}: {label} versus {new['label']}"
            )
        if new:
            old_topology = row.get("per_residue_topology") if row["found"] else None
            reconciliation.append(
                {
                    "species": species,
                    "sequence_sha256": sha,
                    "state": "new_validated_v1" if not label else "frozen_label_retained",
                    "same_topology": old_topology == new["topology"]
                    if isinstance(old_topology, str)
                    else None,
                    "frozen_label": label,
                    "supplement_label": new["label"],
                }
            )
        evidence[sha] = {
            "label_source": "frozen_v1_snapshot"
            if label
            else "validated_v1_supplement"
            if new
            else "pending",
            "supplement_runs": new["provenance"] if new else [],
        }
        labels[sha] = label or (new["label"] if new else "")
    return labels, evidence, reconciliation
