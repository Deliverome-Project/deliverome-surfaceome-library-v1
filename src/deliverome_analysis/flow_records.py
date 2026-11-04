"""Small Git-reviewable records derived from completed flow runs, not new gate state."""

import csv
import hashlib
import json
import os
import re
import shlex
import subprocess
import tempfile
from pathlib import Path

import yaml

#: Shape of `run-records/<id>.yaml`. Anything that writes one reads this rather than repeating
#: the number, so a bump cannot land in one writer and not another. `deliverome_analysis.ma900`
#: writes the same shape for sorts the instrument ran, and tests/test_ma900.py pins this value --
#: changing it there is the prompt to check every writer still emits what the new version needs.
RUN_RECORD_SCHEMA_VERSION = 2


def digest(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def run_directory(output):
    root = Path(output).resolve()
    return root if (root / "run.json").is_file() else root / "run"


def repository(path):
    path = Path(path).resolve()
    directory = path if path.is_dir() else path.parent
    result = subprocess.run(
        ["git", "-C", str(directory), "rev-parse", "--show-toplevel"],
        capture_output=True,
        text=True,
        check=False,
    )
    return Path(result.stdout.strip()) if result.returncode == 0 else None


def relative(path, root):
    path = Path(path).resolve()
    return path.relative_to(root).as_posix() if root and path.is_relative_to(root) else None


def source_file(path, root, role):
    return {"role": role, "path": relative(path, root), "sha256": digest(path)}


def capture_execution(
    samples, recipe, output, experiment, run_id, *, local_only=False, archive=False, sources=()
):
    """Capture before analysis; never infer execution identity when publishing old runs."""
    root = repository(recipe)
    identity = {"git_commit": None, "working_tree_dirty": None}
    if root:
        identity["git_commit"] = subprocess.check_output(
            ["git", "-C", str(root), "rev-parse", "HEAD"],
            text=True,
        ).strip()
    paths = [(Path(recipe), "recipe"), (Path(samples), "samples")]
    if root:
        paths.extend(
            (root / name, name) for name in ("uv.lock", "pyproject.toml") if (root / name).is_file()
        )
    paths.extend(
        (p, "input_manifest")
        for p in sorted(Path(samples).resolve().parent.glob("input_manifest.*"))
        if p.suffix in {".txt", ".csv", ".tsv", ".yaml", ".yml", ".json"}
    )
    paths.extend(
        (p, "analysis_script") for p in sorted(Path(samples).resolve().parent.glob("*.py"))
    )
    with Path(samples).open(newline="") as stream:
        for row in csv.DictReader(stream):
            if row.get("compensation_path"):
                paths.append(
                    (Path(samples).resolve().parent / row["compensation_path"], "compensation")
                )
    paths.extend((Path(p), "additional_source") for p in sources)
    files = []
    seen = set()
    for path, role in paths:
        if path.resolve() not in seen:
            files.append(source_file(path, root, role))
            seen.add(path.resolve())
    if root:
        # Generated output directories must not be mistaken for changed source.
        scoped_paths = sorted(
            {"src", "pyproject.toml", "uv.lock"} | {f["path"] for f in files if f["path"]}
        )
        identity["working_tree_dirty"] = bool(
            subprocess.check_output(
                ["git", "-C", str(root), "status", "--porcelain", "--", *scoped_paths],
                text=True,
            ).strip()
        )
    args = [
        "uv",
        "run",
        "--locked",
        "--group",
        "flow",
        "python",
        "-m",
        "deliverome_analysis.flow",
        "run",
        "--experiment",
        experiment,
        "--samples",
        relative(samples, root),
        "--recipe",
        relative(recipe, root),
        "--out",
        relative(output, root),
        "--run-id",
        run_id,
    ]
    if local_only:
        args.append("--local-only")
    if archive:
        args.append("--zip")
    for path in sources:
        args.extend(["--source", relative(path, root)])
    return {
        "schema_version": 1,
        "analysis_source": identity,
        "source_files": files,
        "rerun": {
            "working_directory": "repository root",
            "command": shlex.join(args) if all(a is not None for a in args) else None,
            "note": "Equivalent wrapper command; use a NEW output directory and run ID.",
        },
    }


def save_execution(output, context):
    run = run_directory(output)
    provenance = json.loads((run / "run.json").read_text())
    for role, field in [("recipe", "recipe_sha256"), ("samples", "samples_sha256")]:
        saved = next(f for f in context["source_files"] if f["role"] == role)
        if saved["sha256"] != provenance.get(field):
            raise ValueError(f"{role} differs from the source captured before execution")
    (run / "execution-context.json").write_text(json.dumps(context, indent=2) + "\n")


def default_record_directory(output, experiment):
    root = repository(output)
    experiment = experiment.removeprefix("analysis/").strip("/")
    if root and (root / "analysis" / experiment).is_dir():
        return root / "analysis" / experiment / "run-records"
    parent = Path(output).resolve().parent
    # Keep the record outside the ignored generated runs/ folder, even without Git.
    return (parent.parent if parent.name == "runs" else parent) / "run-records"


def write_record(
    output,
    experiment,
    run_id,
    *,
    published_uri=None,
    publication=None,
    publication_sha256=None,
    record_dir=None,
):
    """Use immutable run provenance; current checkout is never treated as historical proof."""
    if published_uri and (publication is None or not publication_sha256):
        raise ValueError("A published record requires its verified completion manifest")
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", run_id) or run_id in {".", ".."}:
        raise ValueError("Run ID must be one path component")
    run = run_directory(output)
    provenance = json.loads((run / "run.json").read_text())
    if not provenance.get("inputs"):
        raise ValueError("Run has no recorded inputs")
    context_path = run / "execution-context.json"
    context = json.loads(context_path.read_text()) if context_path.exists() else None
    files = context["source_files"] if context else []
    if context:
        for role, key in [("recipe", "recipe_sha256"), ("samples", "samples_sha256")]:
            matches = [f for f in files if f["role"] == role]
            if len(matches) != 1 or matches[0]["sha256"] != provenance.get(key):
                raise ValueError("Execution context does not match the completed run")
    warnings = []
    if not context:
        warnings.append(
            "Legacy run: analysis commit, lockfile and original command were not captured."
        )
    if context and not context["analysis_source"].get("git_commit"):
        warnings.append("Analysis Git identity was unavailable at execution.")
    if context and context["analysis_source"].get("working_tree_dirty"):
        warnings.append(
            "Analysis checkout had uncommitted files; commit alone does not reproduce it."
        )
    if any(f["path"] is None for f in files):
        warnings.append(
            "Some sources were outside the repository; move them into the experiment for future runs."
        )
    snapshots = {
        name: {"sha256": digest(run / name)}
        for name in (
            "run.json",
            "recipe.json",
            "samples.csv",
            "reproducibility.yaml",
            "execution-context.json",
        )
        if (run / name).is_file()
    }
    if publication is not None:
        mapping = publication.get("source_path_mapping", {})
        prefix = publication["prefix"]
        for name, entry in snapshots.items():
            original = (run / name).relative_to(Path(output).resolve()).as_posix()
            target = mapping.get(original)
            match = next(
                (f for f in publication["files"] if f["key"] == f"{prefix}/{target}"), None
            )
            if not match or match["sha256"] != entry["sha256"]:
                raise ValueError(f"Local {name} does not match the published run")
    record = {
        "schema_version": RUN_RECORD_SCHEMA_VERSION,
        "experiment": experiment.removeprefix("analysis/").strip("/"),
        "run_id": run_id,
        "agentflow": provenance.get("agentflow"),
        "analysis_source": context["analysis_source"] if context else None,
        "source_files": files,
        "run_sha256": snapshots["run.json"]["sha256"],
        "input_count": len(provenance["inputs"]),
        "publication": {"uri": published_uri, "manifest_sha256": publication_sha256}
        if published_uri
        else None,
        "rerun": context["rerun"] if context else None,
        "warnings": warnings,
    }
    if published_uri:
        destination = (
            Path(record_dir) if record_dir else default_record_directory(output, experiment)
        )
        path = destination / f"{run_id}.yaml"
    else:
        # Local exploration must not create extra Git artifacts. Publication creates
        # the single compact Git record; this preview stays with generated outputs.
        destination = run
        path = run / "run-record.yaml"
    destination.mkdir(parents=True, exist_ok=True)
    text = yaml.safe_dump(record, sort_keys=False, allow_unicode=True)
    if path.exists():
        if path.read_text() != text:
            raise ValueError(
                "Run record already exists with different content; keep historical records immutable"
            )
        return path
    # Publish a complete file without overwriting a concurrent/historical record.
    with tempfile.NamedTemporaryFile(mode="w", dir=destination, delete=False) as stream:
        temporary = Path(stream.name)
        try:
            stream.write(text)
            stream.flush()
            os.fsync(stream.fileno())
        except BaseException:
            temporary.unlink(missing_ok=True)
            raise
    try:
        os.link(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)
    return path
