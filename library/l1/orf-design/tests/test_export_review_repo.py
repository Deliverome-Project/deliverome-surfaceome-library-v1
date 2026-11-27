"""The L1 review-repo export is import-closed: no exported script misses a local or package import."""

from __future__ import annotations

import ast
import hashlib
import json
import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE / "reviewer-package"))

import export_review_repo as err  # noqa: E402


def _is_exported_review_repo() -> bool:
    return (err.REPO / err.MANIFEST).is_file()


@pytest.mark.skipif(
    _is_exported_review_repo(), reason="exporter self-test requires the source repository"
)
def test_exported_scripts_are_import_closed_at_head():
    files = err.collect("HEAD", include_historical=False)
    tracked = set(err.git("ls-tree", "-r", "--name-only", "HEAD").decode().splitlines())
    exported = set(files)

    # 1. Every same-directory .py module imported by any exported script is also exported.
    missing_siblings: dict[str, list[str]] = {}
    for path, data in sorted(files.items()):
        if not path.endswith(".py") or path.startswith(err.PKG + "/"):
            continue
        sibs = err._imported_siblings(path, data, tracked)
        miss = sorted(sibs - exported)
        if miss:
            missing_siblings[path] = miss
    assert not missing_siblings, f"exported scripts miss sibling imports: {missing_siblings}"

    # 2. The s5e13 internal-arm script's sibling helpers (which used to be missing) are included.
    for helper in (
        "library/l1/tag-site-benchmarking/contact_veto.py",
        "library/l1/tag-site-benchmarking/tag_site_model_v4.py",
        "library/l1/tag-site-benchmarking/tag_site_model_v5.py",
    ):
        assert helper in exported, f"expected {helper} in review export"

    # 3. Every deliverome_analysis module imported by any exported script is exported.
    pkg_mods = {
        p[len(err.PKG) + 1 :].split("/")[0].removesuffix(".py")
        for p in exported
        if p.startswith(err.PKG + "/")
    }
    missing_pkg: dict[str, list[str]] = {}
    for path, data in sorted(files.items()):
        if not path.endswith(".py"):
            continue
        in_pkg = path.startswith(err.PKG + "/")
        needed = err._imported_modules(data, in_package=in_pkg)
        miss = sorted(m for m in needed if m not in pkg_mods and f"{err.PKG}/{m}.py" in tracked)
        if miss:
            missing_pkg[path] = miss
    assert not missing_pkg, f"exported scripts miss package modules: {missing_pkg}"


def _is_shallow() -> bool:
    return err.git("rev-parse", "--is-shallow-repository").decode().strip() == "true"


@pytest.mark.skipif(
    _is_exported_review_repo(), reason="exporter self-test requires the source repository"
)
def test_clean_paths_export_strips_internal_markers_and_stays_valid():
    files = err.collect("HEAD", clean_paths=True, include_historical=not _is_shallow())
    internal_dir = re.compile(r"(?:^|/)sprint\d+/|(?:^|/)s\d+e\d+-|(?:^|/)scratch/becca/")
    bad_paths = [p for p in files if internal_dir.search(p)]
    assert not bad_paths, f"exported paths still carry internal markers: {bad_paths[:10]}"

    # Every clean-path exported Python file parses cleanly with ast.parse.
    for path, data in sorted(files.items()):
        if path.endswith(".py"):
            ast.parse(data, filename=path)

    # S3 keys (without leading analysis/) remain intact so the offline bundle still replays.
    syn = files["library/l1/orf-design/synthesis_compliance.py"].decode()
    assert '"sprint4/s4e9-orf-design-fg-panel/synthesis/' in syn
    # Repo-relative references to s5e13 were rewritten to library/l1/tag-site-benchmarking.
    internal_arms = files["library/l1/orf-design/build_l1_internal_arms.py"].decode()
    assert '"l1" / "tag-site-benchmarking"' in internal_arms
    assert "s5e13-tag-site-benchmarking" not in internal_arms

    # Historical commit tarballs also remap member paths so run_stage.py resolves dest / riu.AP.
    if not _is_shallow():
        import io
        import tarfile

        with tarfile.open(fileobj=io.BytesIO(files["historical/88786350.tar"])) as tar:
            names = set(tar.getnames())
        assert "library/l1/proteoform-selection/isoform_universe.py" in names
        assert "library/l1/tag-site-architecture/library_members.json" in names
        assert not [n for n in names if "scratch/becca/" in n]


@pytest.mark.skipif(
    _is_exported_review_repo(), reason="exporter self-test requires the source repository"
)
def test_exported_clean_paths_repo_runs_offline_tests(tmp_path: Path):
    """End-to-end replay smoke test: export HEAD with clean paths into a fresh
    directory (no .git history, no warm cache, dead S3 endpoint) and run its
    offline test suite."""
    repo = tmp_path / "l1-review-export"
    full_history = not _is_shallow()
    err.export(
        repo, "HEAD", overwrite_edited=True, clean_paths=True, include_historical=full_history
    )

    # Generated container, CI, docs, and manifest files are present and digest-pinned.
    for rel in ("Dockerfile", ".github/workflows/ci.yml", "README.md", "EXPORT_MANIFEST.json"):
        assert (repo / rel).is_file(), f"missing generated file in export: {rel}"

    dockerfile = (repo / "Dockerfile").read_text()
    assert "@sha256:" in dockerfile
    assert "ghcr.io/astral-sh/uv:" in dockerfile
    assert "--no-group flow" in dockerfile

    # Run exported offline tests inside the exported tree with no git history,
    # no warm cache, and an unreachable S3 endpoint.
    fake_home = tmp_path / "empty-home"
    fake_home.mkdir()
    env = {
        "PATH": os.environ.get("PATH", ""),
        "HOME": str(fake_home),
        "PYTHONPATH": str(repo / "src"),
        "S3_ENDPOINT_URL": "http://127.0.0.1:9",
        "PYTHONHASHSEED": "0",
    }
    test_targets = [
        "tests/test_provenance_snapshots.py",
        "tests/test_s3_offline_mirror.py",
        "library/l1/orf-design/tests/test_brief_no_experiment_ids.py",
        "library/l1/orf-design/tests/test_l1_shared_definitions.py",
        "library/l1/orf-design/tests/test_mint_zenodo_archive.py",
    ]
    res = subprocess.run(
        [sys.executable, "-m", "pytest", "-q", *test_targets],
        cwd=repo,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )
    assert res.returncode == 0, f"exported offline tests failed:\n{res.stdout}\n{res.stderr}"

    # Every relative markdown link in the exported repo must resolve (mirrors review-repo CI).
    # Static root files (LICENSE, LICENSING.md, CITATION.cff) live in the destination review repo
    # rather than being emitted by export_review_repo.py into an empty tmp_path.
    static_review_repo_root_files = {"LICENSE", "LICENSING.md", "CITATION.cff"}
    link_re = re.compile(r"!?\[[^\]]*\]\(([^)]+)\)")
    broken_links = []
    for md in sorted(repo.rglob("*.md")):
        if any(p in {".git", ".venv", "historical", "review-repo"} for p in md.parts):
            continue
        for raw in link_re.findall(md.read_text(errors="ignore")):
            target = raw.split("#", 1)[0].strip().strip("<>")
            if not target or target.startswith(("http://", "https://", "mailto:", "#")):
                continue
            if md == repo / "README.md" and target in static_review_repo_root_files:
                continue
            if not (md.parent / target).exists():
                broken_links.append(f"{md.relative_to(repo)}: {target}")
    assert not broken_links, "Broken relative markdown links in export:\n  " + "\n  ".join(
        broken_links[:20]
    )

    # Every materialized pinned git snapshot is present and matches SHA256SUMS.json, and
    # provenance.read_pinned_text resolves them inside the exported repo with zero .git history.
    if full_history:
        snap_root = repo / "historical" / "git-snapshots"
        manifest = json.loads((snap_root / "SHA256SUMS.json").read_text())
        assert len(manifest) == len(err.PINNED_GIT_SNAPSHOTS)
        for key, expected in manifest.items():
            actual = hashlib.sha256((snap_root / key).read_bytes()).hexdigest()
            assert actual == expected, f"snapshot checksum mismatch: {key}"

        probe = subprocess.run(
            [
                sys.executable,
                "-c",
                (
                    "from pathlib import Path; "
                    "from deliverome_analysis.provenance import read_pinned_text; "
                    "t1 = read_pinned_text('8432c010', 'library/l1/orf-design/tag_refold_structures.tsv', repo=Path('.')); "
                    "t2 = read_pinned_text('8432c0104c8fdaa57e78173d18fc5bb7064d91bd', 'library/l1/orf-design/constructs_l1_cterm.tsv', repo=Path('.')); "
                    "assert len(t1) > 100 and len(t2) > 100"
                ),
            ],
            cwd=repo,
            env=env,
            capture_output=True,
            text=True,
            check=False,
        )
        assert probe.returncode == 0, f"offline pinned-snapshot read failed:\n{probe.stderr}"


def test_no_exported_script_calls_raw_git_show():
    """Every exported script must read historical commits via provenance.read_pinned_*,
    never raw ``git show`` subprocess calls that break without .git history."""
    if _is_exported_review_repo():
        manifest = json.loads((err.REPO / err.MANIFEST).read_text())["files"]
        files = {
            p: (err.REPO / p).read_bytes()
            for p in manifest
            if p.endswith(".py") and (err.REPO / p).is_file()
        }
    else:
        files = err.collect("HEAD", include_historical=False)
    raw_git_show = re.compile(r'["\']git["\'][^\n]*["\']show["\']')
    allowed = {
        "library/l1/orf-design/reviewer-package/export_review_repo.py",
        "src/deliverome_analysis/provenance.py",
    }
    offenders: list[str] = []
    for path, data in sorted(files.items()):
        if not path.endswith(".py") or path in allowed or path.startswith("historical/"):
            continue
        if raw_git_show.search(data.decode("utf-8", errors="replace")):
            offenders.append(path)
    assert not offenders, f"exported scripts still call raw git show: {offenders}"
