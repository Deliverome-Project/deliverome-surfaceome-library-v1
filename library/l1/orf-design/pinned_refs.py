"""Every git commit an s4e9 script pins must be reachable from the checkout's history.

Stages replay committed records with ``git show <ref>:<path>`` (synthesis compliance reads the
rung arm builds at ``INPUT_REF``; the fold check reads the refold table at ``REFOLD_COMMIT``).
A ref that lives only on an unmerged branch makes the order unbuildable from ``main``: the
independent audit of 2026-11-14 found ``INPUT_REF`` pointing at the draft integration branch,
so ``main`` could not rebuild 1,817 of the 2,678 ordered inserts.

``check()`` fails if any pinned ref is not an ancestor of ``HEAD``. It is the pipeline stage
``pinned_refs`` (which needs full history, as ``verify`` does anyway) and is also run by
``tests/test_pinned_refs_on_main.py``, which skips on CI's one-commit checkout.

Requested by the coordinator (2026-11-14); implemented by the Claude Code session
l1-reviewer-package.

    uv run --locked python pinned_refs.py
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
#: The two experiment folders whose scripts build the order: S4E9 and the S5E13 internal-site
#: decisions it reads.
ROOTS = (HERE, HERE.parents[1] / "l1" / "tag-site-benchmarking")
#: ``NAME_REF = "<sha>"`` (also ``_COMMIT``, ``_SHA``, ``_AT``), and ``commit="<sha>"`` keywords.
PIN = re.compile(
    r'^\s*([A-Z][A-Z0-9_]*(?:REF|COMMIT|SHA|_AT))\s*[:=][^"\n]*"([0-9a-f]{7,40})"', re.M
)
KEYWORD_PIN = re.compile(r'\b(commit|ref|rev)\s*=\s*"([0-9a-f]{7,40})"')
#: A script that names a session branch reads whatever that branch holds today, or nothing
#: once it is deleted. A string literal naming one, outside a comment, is refused.
BRANCH = re.compile(r"^[^#\n]*[\"\'](?:origin/)?(claude/[A-Za-z0-9._/-]+)[\"\']", re.M)
#: Branch references known on 2026-11-14 and owned by another session, with the owner. Listed
#: so the check stays strict for everything else; remove each entry when its script is fixed.
BRANCH_PENDING = {
    "../../sprint5/s5e13-tag-site-benchmarking/compute_conservation.py": "ranker-pu-segment",
    "../../sprint5/s5e13-tag-site-benchmarking/run_predictors_for_target.py": "ranker-pu-segment",
}

#: Pins that record what an earlier check read, kept as history, with the reason. Anything
#: that feeds the order must not be listed here.
HISTORICAL = {
    ("tag_refold/verify_ordered_constructs.py", "INTEGRATION_COMMIT"): (
        "l1-tag-refold's 2026-11-12 proof read the draft integration build; superseded by "
        "reviewer-package/fold_order_check.py, which reads the order of record. The ten rung "
        "files at that commit are saved on S3 with SHA-256 (save_rung_builds.py), so the proof "
        "survives the branch"
    ),
    ("membrane-side-gate/gate_inputs.py", "S5E13_COMMIT"): (
        "names where the gate's three vendored S5E13 files came from (PR #180, closed); the "
        "files are committed in membrane-side-gate/ranker-c7be01a4/ and refused unless their "
        "SHA-256 matches SOURCE.json, so nothing is replayed from that commit"
    ),
    ("publish_viewer_release.py", "AGENT_PICKS_COMMIT"): (
        "post-cut surfaceome viewer D1 publication tool reading an exploratory agent-picks table "
        "from PR #233 (claude/agentic-tag-sites-c26525); does not feed the L1 order"
    ),
}


def git(*args: str) -> subprocess.CompletedProcess:
    return subprocess.run(["git", *args], cwd=HERE, capture_output=True, text=True)


def is_shallow() -> bool:
    return git("rev-parse", "--is-shallow-repository").stdout.strip() == "true"


def scripts() -> list[tuple[str, Path]]:
    out = []
    for root in ROOTS:
        for path in sorted(root.rglob("*.py")):
            if {"tests", "review-repo", "data"} & set(path.relative_to(root).parts):
                continue
            out.append((Path(os.path.relpath(path, HERE)).as_posix(), path))
    return out


def pins() -> list[tuple[str, str, str]]:
    out = []
    for rel, path in scripts():
        text = path.read_text(errors="ignore")
        out += [(rel, name, ref) for name, ref in PIN.findall(text)]
        out += [(rel, f"{kw}=", ref) for kw, ref in KEYWORD_PIN.findall(text)]
    return out


def branch_refs() -> dict[str, list[str]]:
    """Scripts that name a session branch as a git ref, keyed by path."""
    out = {}
    for rel, path in scripts():
        found = BRANCH.findall(path.read_text(errors="ignore"))
        if found:
            out[rel] = sorted(set(found))
    return out


def _is_exported_review_repo() -> bool:
    return (REPO / "EXPORT_MANIFEST.json").is_file()


def _check_exported_snapshots() -> str:
    """Verify materialized pinned snapshots in an exported review repo (no git history needed)."""
    branches = branch_refs()
    new = sorted(set(branches) - set(BRANCH_PENDING))
    if new:
        raise SystemExit(
            f"scripts read a session branch by name: {[(b, branches[b]) for b in new]}"
        )
    snap_root = REPO / "historical" / "git-snapshots"
    manifest_path = snap_root / "SHA256SUMS.json"
    if not manifest_path.is_file():
        raise SystemExit(
            f"exported review repo is missing {manifest_path.relative_to(REPO)}; "
            "re-export with historical git snapshots enabled"
        )
    manifest: dict[str, str] = json.loads(manifest_path.read_text())
    if not manifest:
        raise SystemExit(f"{manifest_path.relative_to(REPO)} is empty")

    missing_files, bad_hashes = [], []
    for key, expected in sorted(manifest.items()):
        target = snap_root / key
        if not target.is_file():
            missing_files.append(key)
        elif hashlib.sha256(target.read_bytes()).hexdigest() != expected:
            bad_hashes.append(key)

    # Check that every (commit, path) declared in export_review_repo.PINNED_GIT_SNAPSHOTS is
    # materialized and recorded in SHA256SUMS.json.
    missing_declared: list[str] = []
    exporter_dir = HERE / "reviewer-package"
    if (exporter_dir / "export_review_repo.py").is_file():
        if str(exporter_dir) not in sys.path:
            sys.path.insert(0, str(exporter_dir))
        import export_review_repo as err

        for commit, rel in err.PINNED_GIT_SNAPSHOTS:
            candidates = (f"{commit}/{rel}", f"{commit}/{err._remap_path(rel)}")
            if not any(c in manifest and (snap_root / c).is_file() for c in candidates):
                missing_declared.append(f"{commit}:{rel}")

    # Also check that every active pinned ref used by a script calling read_pinned_* has at
    # least one materialized snapshot directory under historical/git-snapshots/<ref>/.
    script_text = {rel: path.read_text(errors="ignore") for rel, path in scripts()}
    missing_commits: list[str] = []
    for rel, name, ref in pins():
        if (rel, name) in HISTORICAL:
            continue
        text = script_text.get(rel, "")
        if "read_pinned_" not in text:
            continue
        if name != "commit=" and not re.search(
            rf"read_pinned_\w+\([^)]*\b{re.escape(name)}\b", text
        ):
            continue
        if not (snap_root / ref).is_dir():
            missing_commits.append(f"{rel}:{name}={ref}")

    problems = []
    if missing_files:
        problems.append(f"missing snapshot files: {missing_files[:5]}")
    if bad_hashes:
        problems.append(f"SHA-256 mismatch: {bad_hashes[:5]}")
    if missing_declared:
        problems.append(f"declared snapshots not materialized: {missing_declared[:5]}")
    if missing_commits:
        problems.append(f"pinned commits without materialized snapshots: {missing_commits[:5]}")
    if problems:
        raise SystemExit("exported pinned snapshots failed: " + "; ".join(problems))
    return (
        f"exported review repo: all {len(manifest)} materialized git snapshots in "
        f"historical/git-snapshots/ exist and match SHA256SUMS.json"
    )


def check() -> str:
    if _is_exported_review_repo():
        return _check_exported_snapshots()
    if is_shallow():
        raise SystemExit(
            "shallow clone: pinned-ref ancestry needs full history (git fetch --unshallow)"
        )
    branches = branch_refs()
    new = sorted(set(branches) - set(BRANCH_PENDING))
    if new:
        raise SystemExit(
            f"scripts read a session branch by name: {[(b, branches[b]) for b in new]}"
        )
    found = pins()
    assert found, "no pinned refs found; the pattern no longer matches the scripts"
    missing, not_ancestor = [], []
    historical_count = 0
    for rel, name, ref in found:
        if (rel, name) in HISTORICAL:
            historical_count += 1
            continue
        if git("cat-file", "-e", f"{ref}^{{commit}}").returncode:
            missing.append(f"{rel}:{name}={ref}")
        elif git("merge-base", "--is-ancestor", ref, "HEAD").returncode:
            not_ancestor.append(f"{rel}:{name}={ref}")
    if missing or not_ancestor:
        raise SystemExit(
            f"pinned refs not reachable from HEAD: missing {missing}; not ancestors {not_ancestor}"
        )
    pending = sorted(set(BRANCH_PENDING) & set(branches))
    return (
        f"all {len(found) - historical_count} pinned git refs in S4E9 and S5E13 are ancestors of "
        f"HEAD ({historical_count} historical pin(s) listed with a reason); "
        f"{len(pending)} script(s) still name a session branch, owned elsewhere: {pending}"
    )


if __name__ == "__main__":
    print(check())
