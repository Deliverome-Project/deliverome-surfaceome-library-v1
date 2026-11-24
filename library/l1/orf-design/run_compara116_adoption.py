"""Replay the pinned Compara116 pre-tag candidate adoption into a fresh directory.

No downloads, canonical writes, barcode reservations or prediction submissions occur.
Inputs must already be staged, including the historical records and pinned backbone.
The run journal records each command, its status, checksums and scientific counts.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
OLD_CONSTRUCTS_SHA256 = "cdc5fa5e119e7590d910b44db03138a829f7a59e6d8ec430d8ebbe9f21008e3f"
BACKBONE_SHA256 = "57915c1a5db2f0856a6bf6ae1e1d1d9f71661fb37cf68a4938e74cae289485b4"
SCRIPTS = (
    "resolve_compara116_cds.py",
    "adopt_compara116_members.py",
    "build_compara116_constructs.py",
    "finalize_compara116_adoption.py",
    "member_id_allocations.py",
)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def git(*args):
    result = subprocess.run(["git", *args], cwd=REPO, capture_output=True, text=True, check=True)
    return result.stdout.strip()


def verify_pin(path: Path, expected: str):
    actual = digest(path)
    if actual != expected:
        raise ValueError(f"Input digest mismatch: {path.name}: expected {expected}, got {actual}")
    return actual


def run_stage(name: str, argv: list[str], journal: dict, out: Path):
    stage = {"name": name, "argv": argv, "status": "running"}
    journal["stages"].append(stage)
    (out / "run.json").write_text(json.dumps(journal, indent=2) + "\n")
    result = subprocess.run(argv, cwd=REPO, capture_output=True, text=True, check=False)
    (out / f"{name}.stdout.log").write_text(result.stdout)
    (out / f"{name}.stderr.log").write_text(result.stderr)
    stage["returncode"] = result.returncode
    stage["status"] = "passed" if result.returncode == 0 else "failed"
    if result.returncode:
        stage["error"] = result.stderr.strip() or result.stdout.strip()
        raise RuntimeError(f"{name} failed with exit code {result.returncode}")


def counts(out: Path, name: str):
    files = {
        "resolve": ("resolved/compara116_cds_provenance.json", "summary"),
        "plan": ("plan/compara116_member_adoption_summary.json", None),
        "construct": ("constructs/compara116_construct_provenance.json", None),
        "reconcile": ("review/compara116_member_adoption_summary.json", None),
        "finalize": ("bundle/compara116_adoption_validation.json", None),
    }
    path, child = files[name]
    value = json.loads((out / path).read_text())
    return value[child] if child else value


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--input-manifest", type=Path, default=HERE / "compara116_cds_inputs.json")
    for name in ("old-manifest", "old-cds", "old-constructs", "backbone", "out"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--allocation-ledger", type=Path)
    args = parser.parse_args(argv)
    for name, value in vars(args).items():
        if isinstance(value, Path):
            setattr(args, name, value.resolve())
    args.out.mkdir(parents=True, exist_ok=False)
    journal = {
        "scope": "pre-tag Compara116 candidate adoption; no translation-end repairs",
        "status": "validating",
        "stages": [],
        "input_sha256": {},
        "source_commit": git("rev-parse", "HEAD"),
        "source_dirty": bool(git("status", "--porcelain")),
        "script_sha256": {},
    }
    try:
        script_paths = [Path(__file__), *(HERE / name for name in SCRIPTS)]
        script_paths += [
            REPO / "library/l1/proteoform-selection" / name
            for name in ("assemble_candidates.py", "domesticate_candidates.py")
        ]
        journal["script_sha256"] = {str(p.relative_to(REPO)): digest(p) for p in script_paths}
        manifest = json.loads(args.input_manifest.read_text())
        journal["input_sha256"]["input_manifest"] = digest(args.input_manifest)
        for name, item in manifest["files"].items():
            path = args.input_dir / name
            if path.resolve().parent != args.input_dir:
                raise ValueError("Input manifest file names must be direct children of input-dir")
            journal["input_sha256"][f"staged/{name}"] = verify_pin(path, item["sha256"])
        explicit = {
            "old_manifest": (
                args.old_manifest,
                manifest["files"]["candidate_manifest_tranches.tsv"]["sha256"],
            ),
            "old_cds": (args.old_cds, manifest["files"]["candidate_cds.json"]["sha256"]),
            "old_constructs": (args.old_constructs, OLD_CONSTRUCTS_SHA256),
            "backbone": (args.backbone, BACKBONE_SHA256),
        }
        for label, (path, expected) in explicit.items():
            journal["input_sha256"][label] = verify_pin(path, expected)
        if args.allocation_ledger:
            journal["input_sha256"]["allocation_ledger"] = digest(args.allocation_ledger)
        panel = args.input_dir / "ortholog_panel_compara116.tsv"
        resolved, plan, built, review, bundle = (
            args.out / name for name in ("resolved", "plan", "constructs", "review", "bundle")
        )
        # An absent explicit ledger is deliberate. Never pick up an ambient repository
        # default whose appearance would silently change this run's allocation inputs.
        ledger = args.allocation_ledger or args.out / "no-prior-allocation-ledger.tsv"
        adoption = [
            "--panel",
            panel,
            "--manifest",
            args.old_manifest,
            "--existing-cds",
            args.old_cds,
        ]
        commands = [
            (
                "resolve",
                "resolve_compara116_cds.py",
                [
                    "--manifest",
                    args.input_manifest,
                    "--input-dir",
                    args.input_dir,
                    "--out-dir",
                    resolved,
                ],
            ),
            (
                "plan",
                "adopt_compara116_members.py",
                [
                    *adoption,
                    "--allocation-ledger",
                    ledger,
                    "--plan-only",
                    "--out",
                    plan,
                ],
            ),
            (
                "construct",
                "build_compara116_constructs.py",
                [
                    "--panel",
                    panel,
                    "--cds",
                    resolved / "compara116_candidate_cds.json",
                    "--backbone",
                    args.backbone,
                    "--backbone-sha256",
                    BACKBONE_SHA256,
                    "--member-map",
                    plan / "compara116_member_map.tsv",
                    "--out",
                    built,
                ],
            ),
            (
                "reconcile",
                "adopt_compara116_members.py",
                [
                    *adoption,
                    "--allocation-ledger",
                    plan / "compara116_id_allocations.tsv",
                    "--built-rows",
                    built / "compara116_assembly_rows.tsv",
                    "--out",
                    review,
                ],
            ),
            (
                "finalize",
                "finalize_compara116_adoption.py",
                [
                    "--panel",
                    panel,
                    "--manifest",
                    review / "candidate_manifest_tranches_review.tsv",
                    "--member-map",
                    review / "compara116_member_map.tsv",
                    "--old-cds",
                    args.old_cds,
                    "--old-constructs",
                    args.old_constructs,
                    "--selected-cds",
                    resolved / "compara116_candidate_cds.json",
                    "--selected-constructs",
                    built / "compara116_candidate_constructs.json",
                    "--out",
                    bundle,
                ],
            ),
        ]
        journal["status"] = "running"
        for name, script, arguments in commands:
            command = [sys.executable, str(HERE / script), *(str(x) for x in arguments)]
            run_stage(name, command, journal, args.out)
            journal["stages"][-1]["counts"] = counts(args.out, name)
        if (plan / "compara116_member_map.tsv").read_bytes() != (
            review / "compara116_member_map.tsv"
        ).read_bytes():
            raise ValueError("Reconciliation changed the allocated member map")
        journal["status"] = "passed"
    except Exception as exc:
        journal["status"] = "failed"
        journal["error"] = f"{type(exc).__name__}: {exc}"
        raise
    finally:
        journal["output_sha256"] = {
            str(p.relative_to(args.out)): digest(p)
            for p in sorted(args.out.rglob("*"))
            if p.is_file() and p != args.out / "run.json"
        }
        (args.out / "run.json").write_text(json.dumps(journal, indent=2) + "\n")
    print(f"All five stages passed; journal: {args.out / 'run.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
