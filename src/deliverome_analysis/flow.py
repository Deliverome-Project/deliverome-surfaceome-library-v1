"""Run and publish Agentflow analyses under an experiment's processed/agentflow prefix."""

import argparse
import base64
import hashlib
import json
import posixpath
import re
import shutil
import tempfile
import zipfile
from datetime import UTC, datetime
from pathlib import Path
from urllib.parse import urlsplit
from uuid import uuid4

from .flow_records import capture_execution, run_directory, save_execution, write_record
from .s3 import processed


def run_prefix(experiment, run_id=None):
    experiment = experiment.removeprefix("analysis/").strip("/")
    match = re.fullmatch(r"sprint(\d+)/s(\d+)e\d+-[a-z0-9]+(?:-[a-z0-9]+)*", experiment)
    if not match or match[1] != match[2]:
        raise ValueError("Use the existing Notion-backed experiment path: sprintN/sNeM-slug")
    run_id = run_id or datetime.now(UTC).strftime("%Y%m%dT%H%M%SZ-") + uuid4().hex[:8]
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", run_id) or run_id in {".", ".."}:
        raise ValueError("Run ID must be one path component")
    return f"{experiment}/agentflow/{run_id}"


def checksum(path):
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def organize_workspace(output, destination):
    """Create a browsable copy; preserve original inputs and rewrite report links."""
    root, destination = Path(output).resolve(), Path(destination)
    run = root if (root / "run.json").is_file() else root / "run"
    destination.mkdir(parents=True, exist_ok=False)
    mapping = {}
    for path in sorted(root.rglob("*")):
        if not path.is_file():
            continue
        if path.is_symlink() or path.suffix.lower() == ".fcs":
            raise ValueError(
                "Publish outputs only; keep FCS inputs and symlinks outside the output folder"
            )
        if path == run / "run-record.yaml":
            continue  # Local-only preview; the published Git record is generated separately.
        if path.suffix.lower() == ".zip":
            continue  # Optional archive is regenerated from the organized workspace.
        original = path.relative_to(root)
        relative = path.relative_to(run) if path.is_relative_to(run) else original
        name = path.name
        if name == "reproducibility.yaml" and path.parent == run:
            target = Path("gates.yaml")
        elif name.endswith(".gatingml.xml"):
            target = Path("gatingml") / name
        elif path.suffix.lower() in {".png", ".pdf", ".svg"}:
            if original.parts[0] == "figures":
                target = original
            else:
                folder = "acquisition" if name.startswith("time-") else "gates"
                target = Path("figures") / folder / name
        elif path.suffix.lower() == ".json" or path.suffix.lower() in {".yaml", ".yml"}:
            target = Path("provenance") / relative
        elif name == "README.md":
            target = Path("provenance/source-README.md")
        else:
            target = relative
        target_path = destination / target
        if target_path.exists() and checksum(target_path) != checksum(path):
            raise ValueError(f"Output paths collide at {target}")
        target_path.parent.mkdir(parents=True, exist_ok=True)
        digest = checksum(path)
        shutil.copy2(path, target_path)
        if checksum(target_path) != digest or checksum(path) != digest:
            raise ValueError("Output changed while organizing workspace")
        mapping[str(original)] = str(target)
    for source, target in mapping.items():
        if not source.endswith(".html"):
            continue

        def rewrite(match, source=source, target=target):
            value = match[3]
            parsed = urlsplit(value)
            if parsed.scheme or parsed.netloc or not parsed.path:
                return match[0]
            original = posixpath.normpath(posixpath.join(posixpath.dirname(source), parsed.path))
            if original not in mapping:
                if match[1].lower().startswith("src"):
                    raise ValueError(f"Missing report image: {original}")
                return match[0]
            image = destination / mapping[original]
            if match[1].lower().startswith("src") and image.suffix.lower() == ".png":
                encoded = base64.b64encode(image.read_bytes()).decode("ascii")
                return match[1] + match[2] + "data:image/png;base64," + encoded + match[2]
            url = posixpath.relpath(mapping[original], posixpath.dirname(target) or ".")
            if parsed.query:
                url += "?" + parsed.query
            if parsed.fragment:
                url += "#" + parsed.fragment
            return match[1] + match[2] + url + match[2]

        path = destination / target
        path.write_text(re.sub(r"((?:href|src)=)([\"'])(.*?)(?:\2)", rewrite, path.read_text()))
    (destination / "README.md").write_text(
        "# Agentflow analysis workspace\n\n"
        "Start with report.html. Gates are in gates.yaml; samples.csv records sample identities and source locations.\n\n"
        "- summary.csv: population counts and signal statistics\n"
        "- events.parquet: all events and gate membership (already column-compressed)\n"
        "- quality.csv: acquisition checks\n"
        "- figures/: summary figures, gates/, and acquisition/ plots\n"
        "- gatingml/: per-sample gate interchange files\n"
        "- provenance/: original run metadata and input recipe snapshots\n"
        "- publication.json: completion record, checksums and original-to-published path mapping\n\n"
        "Report images are embedded: download report.html alone to view its plots offline. Download the whole folder for linked tables, gates and provenance. "
        "Private S3 is storage, not a publicly hosted website. No ZIP is required.\n\n"
        "For reruns, retrieve raw FCS inputs from their sample-sheet S3 locations, update local fcs_path entries, "
        "and use the saved gates.yaml with the recorded Agentflow version. Raw FCS inputs are not bundled. "
        "Publication is not scientific approval; original metadata interpretations remain as recorded.\n"
    )
    return mapping


def publish_run(output, experiment, run_id=None, *, bucket=None, archive=False, record_dir=None):
    """Publish a browsable completed run; optional ZIP and completion manifest last."""
    bucket = processed if bucket is None else bucket
    root = Path(output).resolve()
    run = root if (root / "run.json").is_file() else root / "run"
    if not all(
        (run / name).is_file()
        for name in ["run.json", "events.parquet", "summary.csv", "reproducibility.yaml"]
    ):
        raise ValueError(
            "Publication requires a completed Agentflow run with events, summary and provenance"
        )
    provenance = json.loads((run / "run.json").read_text())
    if not provenance.get("inputs"):
        raise ValueError("Run has no recorded inputs")
    prefix = run_prefix(experiment, run_id)
    if bucket.list_keys(prefix=prefix + "/"):
        raise ValueError("Run prefix already exists; choose a new run ID")
    manifest = {
        "status": "complete",
        "layout_version": 2,
        "prefix": prefix,
        "agentflow": provenance.get("agentflow"),
        "input_count": len(provenance["inputs"]),
        "files": [],
        "note": "Saved run, not scientific approval. Browse or download without extracting an archive.",
    }
    with tempfile.TemporaryDirectory(prefix="agentflow-publish-") as temporary:
        organized = Path(temporary) / "workspace"
        manifest["source_path_mapping"] = organize_workspace(root, organized)
        if (organized / "publication.json").exists():
            raise ValueError("Output paths collide with publication destinations")
        if archive:
            with zipfile.ZipFile(
                Path(temporary) / "agentflow-workspace.zip", "w", compression=zipfile.ZIP_DEFLATED
            ) as handle:
                for path in sorted(organized.rglob("*")):
                    if path.is_file():
                        handle.write(path, path.relative_to(organized))
            shutil.move(
                Path(temporary) / "agentflow-workspace.zip", organized / "agentflow-workspace.zip"
            )
        for path in sorted(organized.rglob("*")):
            if not path.is_file():
                continue
            relative = str(path.relative_to(organized))
            digest = checksum(path)
            bucket.upload(path, prefix + "/" + relative)
            manifest["files"].append(
                {"key": prefix + "/" + relative, "sha256": digest, "bytes": path.stat().st_size}
            )
        completion = Path(temporary) / "publication.json"
        completion.write_text(json.dumps(manifest, indent=2) + "\n")
        bucket.upload(completion, prefix + "/publication.json", content_type="application/json")
        record = write_record(
            output,
            experiment,
            prefix.rsplit("/", 1)[1],
            published_uri="s3://deliverome-processed/" + prefix + "/",
            publication=manifest,
            publication_sha256=checksum(completion),
            record_dir=record_dir,
        )
        print(
            json.dumps(
                {
                    "run_record": str(record),
                    "action": "Review and commit this record; keep its referenced source files in Git",
                }
            )
        )
    # Refresh the shared flow index in the bucket just published to (never another one: tests
    # publish to fakes). A failure must not fail a finished publication; the nightly rebuild
    # catches up.
    if hasattr(bucket, "list_etags"):
        try:
            from deliverome_analysis import qcdb

            qcdb.build("flow", bucket=bucket)
        except Exception as exc:
            print(f"flow index not rebuilt ({type(exc).__name__}: {exc}); run: qcdb build flow")
    return "s3://deliverome-processed/" + prefix + "/"


def run_experiment(
    samples,
    recipe,
    output,
    experiment,
    *,
    run_id=None,
    local_only=False,
    archive=False,
    sources=(),
    record_dir=None,
):
    """Run the pinned Agentflow engine and publish by default; no Qt dependency."""
    from agentflow.batch import run_batch

    prefix = run_prefix(experiment, run_id)
    if not local_only and processed.list_keys(prefix=prefix + "/"):
        raise ValueError("Run prefix already exists; choose a new run ID")
    sources = tuple(sources)
    options = {"local_only": local_only, "archive": archive, "sources": sources}
    identifier = prefix.rsplit("/", 1)[1]
    context = capture_execution(samples, recipe, output, experiment, identifier, **options)
    statistics = run_batch(samples, recipe, output)
    if context != capture_execution(samples, recipe, output, experiment, identifier, **options):
        raise ValueError(
            "Analysis sources or Git identity changed during execution; do not publish this run"
        )
    save_execution(output, context)
    destination = (
        None
        if local_only
        else publish_run(output, experiment, identifier, archive=archive, record_dir=record_dir)
    )
    if local_only:
        record = write_record(output, experiment, identifier, record_dir=record_dir)
        print(json.dumps({"run_record": str(record), "status": "local_only"}))
    return statistics, destination


def record_existing_run(output, experiment, published_uri, *, bucket=None, record_dir=None):
    """Backfill a Git record after verifying the existing S3 publication, without uploading."""
    bucket = processed if bucket is None else bucket
    parsed = urlsplit(published_uri)
    if (
        parsed.scheme != "s3"
        or parsed.netloc != "deliverome-processed"
        or parsed.query
        or parsed.fragment
    ):
        raise ValueError("Use the published s3://deliverome-processed/.../ run URI")
    prefix = parsed.path.strip("/")
    identifier = prefix.rsplit("/", 1)[-1]
    if run_prefix(experiment, identifier) != prefix:
        raise ValueError("Published URI does not belong to this experiment")
    marker_path = bucket.download(prefix + "/publication.json")
    manifest = json.loads(Path(marker_path).read_text())
    if manifest.get("status") != "complete" or manifest.get("prefix") != prefix:
        raise ValueError("Publication is not complete or does not match the requested run")
    run = run_directory(output)
    if not (run / "run.json").is_file():
        raise ValueError("Provide the original completed local run directory")
    return write_record(
        output,
        experiment,
        identifier,
        published_uri="s3://deliverome-processed/" + prefix + "/",
        publication=manifest,
        publication_sha256=checksum(marker_path),
        record_dir=record_dir,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    for name in ["run", "publish", "record"]:
        command = sub.add_parser(name)
        command.add_argument("--experiment", required=True)
        if name != "record":
            command.add_argument("--run-id")
        command.add_argument(
            "--record-dir", help="Override the Git record folder for published runs"
        )
        if name != "record":
            command.add_argument(
                "--zip", action="store_true", help="Also create an optional workspace ZIP"
            )
        command.add_argument("--out", required=True, help="Local run directory")
        if name == "record":
            command.add_argument("--published-uri", required=True)
        if name == "run":
            command.add_argument("--samples", required=True)
            command.add_argument("--recipe", required=True)
            command.add_argument("--local-only", action="store_true")
            command.add_argument(
                "--source",
                action="append",
                default=[],
                help="Additional source/control configuration to hash and list for Git",
            )
    args = parser.parse_args()
    if args.command == "record":
        record = record_existing_run(
            args.out, args.experiment, args.published_uri, record_dir=args.record_dir
        )
        print(json.dumps({"run_record": str(record), "action": "Review and commit"}))
        return
    if args.command == "publish":
        destination = publish_run(
            args.out, args.experiment, args.run_id, archive=args.zip, record_dir=args.record_dir
        )
    else:
        _, destination = run_experiment(
            args.samples,
            args.recipe,
            args.out,
            args.experiment,
            run_id=args.run_id,
            local_only=args.local_only,
            archive=args.zip,
            sources=args.source,
            record_dir=args.record_dir,
        )
    print(json.dumps({"output": args.out, "published": destination}))


if __name__ == "__main__":
    main()
