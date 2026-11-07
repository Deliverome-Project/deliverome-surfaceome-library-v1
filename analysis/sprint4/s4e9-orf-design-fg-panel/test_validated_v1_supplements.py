import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
spec = importlib.util.spec_from_file_location(
    "validated_supplements", HERE / "validated_v1_supplements.py"
)
supplements = importlib.util.module_from_spec(spec)
spec.loader.exec_module(supplements)


def hash_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fixture(folder, monkeypatch, *, label="TM", topology="MM", run_id="round1"):
    folder.mkdir(parents=True)
    runner = b"committed v1 runner"
    runner_sha = hashlib.sha256(runner).hexdigest()
    monkeypatch.setattr(supplements, "committed_runner_blob", lambda _commit: runner)
    protein = "MA"
    sha = hashlib.sha256(protein.encode()).hexdigest()
    fasta = folder / "submitted.fasta"
    fasta.write_text(f">{sha}\n{protein}\n")
    output = folder / "predictions.3line"
    output.write_text(f">{sha}|{label}\n{protein}\n{topology}\n")
    identity = folder / "runtime.json"
    common = {
        "tool_version": "deeptmhmm-1.0-academic",
        "model_version": "1.0",
        "version_source": 'predict.py ArgumentParser description="DeepTMHMM v1.0"',
        "runner_sha256": runner_sha,
        "package_predict_sha256": hashlib.sha256(b"package predict.py").hexdigest(),
        "package_files_sha256": {
            name: hashlib.sha256(name.encode()).hexdigest()
            for name in ["predict.py", "esm.pt"] + [f"deeptmhmm_cv_{i}.model" for i in range(5)]
        },
        "model_files": [f"deeptmhmm_cv_{i}.model" for i in range(5)],
        "esm_checkpoint_files": ["esm.pt"],
    }
    common["package_files_sha256"]["predict.py"] = common["package_predict_sha256"]
    identity.write_text(
        json.dumps(
            {
                **common,
                "schema_version": "deeptmhmm-academic-sweep-1",
                "submitted_fasta_sha256": hash_file(fasta),
                "output_sha256": hash_file(output),
                "failed_shards": 0,
                "shards": [
                    {
                        "tag": "s0",
                        "n": 1,
                        "ok": True,
                        "runtime_identity": common,
                        "submitted_fasta_sha256": hash_file(fasta),
                        "output_sha256": hash_file(output),
                        "sequence_sha256s": [sha],
                    }
                ],
            }
        )
    )
    record = {
        "schema_version": "ortholog-v1-run-1",
        "run_id": run_id,
        "tool_version": "deeptmhmm-1.0-academic",
        "job_id": "ap-test",
        "runner": {"path": supplements.RUNNER_PATH, "git_commit": "a" * 40, "sha256": runner_sha},
        "submitted_fasta": {"path": fasta.name, "sha256": hash_file(fasta)},
        "output": {"path": output.name, "sha256": hash_file(output)},
        "runtime_identity": {"path": identity.name, "sha256": hash_file(identity)},
        "members": [{"sequence_sha256": sha, "species": "mouse"}],
    }
    run = folder / "run.json"
    run.write_text(json.dumps(record))
    return output, run, sha, record


def test_valid_supplement_is_exact_and_species_scoped(tmp_path, monkeypatch):
    output, run, sha, _ = fixture(tmp_path / "r1", monkeypatch)
    labels, report = supplements.load([output], [run])
    assert labels["mouse"][sha]["label"] == "TM"
    assert sha not in labels["cyno"]
    assert report["runs"][0]["proteins"] == 1
    assert labels["mouse"][sha]["provenance"][0]["output_record_number"] == 1


def test_supplement_without_paired_provenance_fails(tmp_path):
    with pytest.raises(ValueError, match="corresponding"):
        supplements.load([tmp_path / "out"], [])


@pytest.mark.parametrize("field", ["tool", "runner", "output_hash", "runtime", "membership"])
def test_incorrect_predictor_or_provenance_is_rejected(tmp_path, monkeypatch, field):
    output, run, sha, record = fixture(tmp_path / "r1", monkeypatch)
    if field == "tool":
        record["tool_version"] = "deeptmhmm-2"
    elif field == "runner":
        record["runner"]["path"] = "scripts/deeptmhmm_v2_app.py"
    elif field == "output_hash":
        record["output"]["sha256"] = "0" * 64
    elif field == "runtime":
        identity = run.parent / "runtime.json"
        data = json.loads(identity.read_text())
        data["model_version"] = "2.0"
        identity.write_text(json.dumps(data))
        record["runtime_identity"]["sha256"] = hash_file(identity)
    else:
        record["members"] = [{"sequence_sha256": sha, "species": "rhesus"}]
    run.write_text(json.dumps(record))
    with pytest.raises(ValueError):
        supplements.load([output], [run])


def test_complete_submitted_coverage_and_residue_format_are_required(tmp_path, monkeypatch):
    output, run, _, record = fixture(tmp_path / "r1", monkeypatch, topology="MX")
    with pytest.raises(ValueError, match="topology"):
        supplements.load([output], [run])
    output.write_text("")
    record["output"]["sha256"] = hash_file(output)
    identity = run.parent / "runtime.json"
    data = json.loads(identity.read_text())
    data["output_sha256"] = hash_file(output)
    identity.write_text(json.dumps(data))
    record["runtime_identity"]["sha256"] = hash_file(identity)
    run.write_text(json.dumps(record))
    with pytest.raises(ValueError, match="cover"):
        supplements.load([output], [run])


def test_new_call_does_not_override_conflicting_frozen_label(tmp_path, monkeypatch):
    output, run, sha, _ = fixture(tmp_path / "r1", monkeypatch)
    calls, _ = supplements.load([output], [run])
    with pytest.raises(ValueError, match="Conflicting frozen/new"):
        supplements.overlay("mouse", [{"key": sha, "found": True, "deeptmhmm_label": "SP"}], calls)
    labels, evidence, report = supplements.overlay(
        "mouse",
        [{"key": sha, "found": True, "deeptmhmm_label": "TM", "per_residue_topology": "MI"}],
        calls,
    )
    assert labels[sha] == "TM"
    assert evidence[sha]["label_source"] == "frozen_v1_snapshot"
    assert report[0]["same_topology"] is False
    assert supplements.overlay("cyno", [{"key": sha, "found": False}], calls)[0][sha] == ""


def test_repeated_rounds_preserve_or_reject_actual_conflicts(tmp_path, monkeypatch):
    first, first_run, _, _ = fixture(tmp_path / "r1", monkeypatch)
    second, second_run, _, _ = fixture(tmp_path / "r2", monkeypatch, topology="MI", run_id="round2")
    _, report = supplements.load([first, second], [first_run, second_run])
    assert report["duplicate_reports"][0]["state"] == "same_label_different_topology"
    third, third_run, _, _ = fixture(
        tmp_path / "r3", monkeypatch, label="SP", topology="SS", run_id="round3"
    )
    with pytest.raises(ValueError, match="Conflicting v1 supplement"):
        supplements.load([first, third], [first_run, third_run])


@pytest.mark.parametrize("field", ["weights", "checkpoint", "shard", "coverage"])
def test_incomplete_runtime_evidence_fails_closed(tmp_path, monkeypatch, field):
    output, run, _, record = fixture(tmp_path / "r1", monkeypatch)
    path = run.parent / "runtime.json"
    data = json.loads(path.read_text())
    if field == "weights":
        del data["package_files_sha256"]["deeptmhmm_cv_0.model"]
    elif field == "checkpoint":
        data["esm_checkpoint_files"] = []
    elif field == "shard":
        data["shards"][0]["ok"] = False
    else:
        data["shards"][0]["sequence_sha256s"] = ["f" * 64]
    path.write_text(json.dumps(data))
    record["runtime_identity"]["sha256"] = hash_file(path)
    run.write_text(json.dumps(record))
    with pytest.raises(ValueError):
        supplements.load([output], [run])


def test_runner_observes_version_and_hashes_every_weight_without_importing_modal(tmp_path):
    import ast
    import re
    from functools import lru_cache

    runner = HERE.parents[2] / "scratch/becca/appris-principal/modal/deeptmhmm_v1_app.py"
    tree = ast.parse(runner.read_text())
    functions = [
        node
        for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name in {"_sha_file", "_runtime_identity"}
    ]
    scope = {
        "Path": Path,
        "hashlib": hashlib,
        "lru_cache": lru_cache,
        "re": re,
        "PACKAGE": str(tmp_path),
        "__file__": str(runner),
    }
    for name in ["esm1b.pt"] + [f"deeptmhmm_cv_{i}.model" for i in range(5)]:
        (tmp_path / name).write_bytes(name.encode())
    predict = tmp_path / "predict.py"
    predict.write_text('parser = argparse.ArgumentParser(description="DeepTMHMM v1.0")')
    exec(compile(ast.Module(body=functions, type_ignores=[]), str(runner), "exec"), scope)
    identity = scope["_runtime_identity"]()
    assert identity["tool_version"] == "deeptmhmm-1.0-academic"
    assert identity["package_files_sha256"]["esm1b.pt"] == hash_file(tmp_path / "esm1b.pt")
    scope["_runtime_identity"].cache_clear()
    predict.write_text('parser = argparse.ArgumentParser(description="DeepTMHMM v2.0")')
    with pytest.raises(ValueError, match=r"v1\.0"):
        scope["_runtime_identity"]()
