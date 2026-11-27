"""Export the L1 review repository from a commit of this one.

The review repo holds everything a reviewer needs to rerun the L1 design, and is also a
working repo: files can be added to it, and the shared utilities edited, as review needs.
So the export manages a known set of files and leaves everything else alone:

* **Managed files** are listed in ``review_repo_files.txt`` (every repository file the
  offline replay opens, found by tracing it; see ``--trace``), the folders in ``ALWAYS``,
  the files that define the environment, and the modules of the shared
  ``deliverome_analysis`` package that the exported code imports (directly or through
  another package module; see ``package_modules``). The rest of the package (DIA
  proteomics, search indexing, flow cytometry, blots) has nothing to do with L1 and stays
  out. They are
  taken from a commit (``git archive``), never from the working tree, so an export always
  corresponds to a commit of this repository.
* **Historical scripts**: the old script versions that two stages rerun are exported as
  ``historical/<commit>.tar``, made by the same ``git archive`` command those stages use
  (``run_stage.py`` reads them).
* ``EXPORT_MANIFEST.json`` in the review repo records the source commit and every managed
  path with its SHA-256. A re-export updates managed files, removes ones no longer
  exported, and **skips any managed file edited in the review repo since the last export**
  (listed, so it can be reconciled by hand; ``--overwrite-edited`` replaces them). Files the
  export did not write are never touched.

    uv run --locked python export_review_repo.py --dest ~/git/deliverome-surfaceome-library-v1
    uv run --locked python export_review_repo.py --dest ... --ref claude/l1-reviewer-package-aa447f
    uv run --locked python export_review_repo.py --trace trace.log   # refresh the file list
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import io
import json
import shutil
import subprocess
import sys
import tarfile
from pathlib import Path

HERE = Path(__file__).parent
S4E9 = HERE.parent
REPO = S4E9.parents[2]
FILES = HERE / "review_repo_files.txt"
TEMPLATES = HERE / "review-repo"
MANIFEST = "EXPORT_MANIFEST.json"

#: Exported whole, on top of the traced list: every script and table of the library's own
#: folders (not only those the replay runs), and the environment. The shared package is not
#: listed: only the modules this code imports are exported (``package_modules``).
ALWAYS = [
    "library/l1/orf-design",
    "library/l1/surface-filters",
    # The vector and the barcode set the brief describes in sections 8 and 9 were designed
    # here, and a reviewer cannot check the barcode filter stack without it.
    "library/l1/vector-design",
    "library/l1/proteoform-selection",
    "library/l1/tag-site-architecture",
    # The internalization-assay schematic in the brief's opening (committed drawing script).
    "library/l1/figures",
    # The surface-proteomics gap analysis the gene section links (lists and the code that builds them).
    "library/l1/proteomics-gap",
    # The internal-site model, its corpus audit truth tables, and its construct rule.
    "library/l1/tag-site-corpus-audit",
    "library/l1/tag-site-benchmarking/compare_internal_policies.py",
    "library/l1/tag-site-benchmarking/hard_negative_benchmark.py",
    "library/l1/tag-site-benchmarking/isoform_fix_redesign.tsv",
    "library/l1/tag-site-benchmarking/compare_internal_policies.json",
    "library/l1/tag-site-benchmarking/l1_internal_additions.py",
    "library/l1/tag-site-benchmarking/l1_internal_arm_decisions.tsv",
    "library/l1/tag-site-benchmarking/l1_internal_arm_ids.tsv",
    "library/l1/tag-site-benchmarking/l1_internal_arms.py",
    "library/l1/tag-site-benchmarking/l1_internal_decisions.py",
    "library/l1/tag-site-benchmarking/l1_internal_sites_vs_v4.tsv",
    "library/l1/tag-site-benchmarking/l1_snorkel_designs.tsv",
    "library/l1/tag-site-benchmarking/ortholog_member_disorder_modal.py",
    "library/l1/tag-site-benchmarking/snorkel_design.py",
    "library/l1/tag-site-benchmarking/snorkel_spacer_test.json",
    "library/l1/tag-site-benchmarking/tag_site_model_v5.json",
    "library/l1/tag-site-benchmarking/tag_site_model_v6.py",
    "library/l1/tag-site-benchmarking/training_gene_families.tsv",
    "library/l1/tag-site-benchmarking/uniprot_site_features.tsv",
    "pyproject.toml",
    "uv.lock",
    ".python-version",
    "tests/test_l1_core_modules.py",
    "tests/test_provenance_snapshots.py",
    "tests/test_s3_offline_mirror.py",
]
PKG = "src/deliverome_analysis"
S5E13 = "library/l1/tag-site-benchmarking"
#: Files this package writes itself, from ``review-repo/`` templates (path in dest -> template).
GENERATED = {
    "README.md": "README.md",
    ".gitignore": "gitignore",
    ".dockerignore": "dockerignore",
    "Dockerfile": "Dockerfile",
    ".github/workflows/ci.yml": "ci.yml",
}
#: Never exported, even if a trace opened them.
EXCLUDE_PREFIXES = (".env", ".venv/", ".git/", "historical/")
#: Clean 3-level directory mapping for the exported review repository so reviewers never see
#: internal sprint/experiment IDs (``s4e9``, ``s5e13``, ``s1e2``, ``s1e3``, ``s3e8``) or
#: personal ``scratch/<name>/`` folders. Every target preserves the 3-level depth so
#: ``Path(__file__).parents[2]`` still resolves to the repository root; S3 keys (without the
#: leading ``analysis/`` or ``scratch/``) are left untouched so the offline S3 bundle replays.
CLEAN_PATH_MAP: dict[str, str] = {
    "library/l1/orf-design": "library/l1/orf-design",
    "library/l1/tag-site-benchmarking": "library/l1/tag-site-benchmarking",
    "library/l1/vector-design": "library/l1/vector-design",
    "library/l1/surface-filters": "library/l1/surface-filters",
    "library/l1/proteomics-gap": "library/l1/proteomics-gap",
    "library/l1/proteoform-selection": "library/l1/proteoform-selection",
    "library/l1/tag-site-architecture": "library/l1/tag-site-architecture",
    "library/l1/tag-site-corpus-audit": "library/l1/tag-site-corpus-audit",
    "library/l1/figures": "library/l1/figures",
}
_CLEAN_REL_SUBS: tuple[tuple[str, str], ...] = (
    ('"library" / "l1" / "orf-design"', '"library" / "l1" / "orf-design"'),
    ("'library' / 'l1' / 'orf-design'", "'library' / 'l1' / 'orf-design'"),
    ('"l1" / "orf-design"', '"l1" / "orf-design"'),
    ("'l1' / 'orf-design'", "'l1' / 'orf-design'"),
    ('/ "l1/orf-design"', '/ "l1/orf-design"'),
    ("/ 'l1/orf-design'", "/ 'l1/orf-design'"),
    ('"l1" / "tag-site-benchmarking"', '"l1" / "tag-site-benchmarking"'),
    ("'l1' / 'tag-site-benchmarking'", "'l1' / 'tag-site-benchmarking'"),
    ('"l1/tag-site-benchmarking"', '"l1/tag-site-benchmarking"'),
    ("'l1/tag-site-benchmarking'", "'l1/tag-site-benchmarking'"),
    ('"l1/vector-design"', '"l1/vector-design"'),
    ("'l1/vector-design'", "'l1/vector-design'"),
    ('"l1/surface-filters"', '"l1/surface-filters"'),
    ("'l1/surface-filters'", "'l1/surface-filters'"),
    ("ranker-c7be01a4", "ranker-c7be01a4"),
)
_CLEAN_MD_REL_SUBS: tuple[tuple[str, str], ...] = (
    ("../../../sprint1/s1e2-fragmid-vector-design/", "../../vector-design/"),
    ("../../../sprint1/s1e3-surface-library-filters/", "../../surface-filters/"),
    ("../../../sprint5/s5e13-tag-site-benchmarking/", "../../tag-site-benchmarking/"),
    ("../../sprint1/s1e2-fragmid-vector-design/", "../vector-design/"),
    ("../../sprint1/s1e3-surface-library-filters/", "../surface-filters/"),
    ("../../sprint5/s5e13-tag-site-benchmarking/", "../tag-site-benchmarking/"),
)
_CLEAN_SUBDIR_MAP: dict[str, str] = {
    "ranker-c7be01a4": "ranker-c7be01a4",
}
#: S3 keys under ``processed://scratch/becca/tag-site-corpus-audit/`` must keep their
#: original S3 prefix so the offline bundle replays, while repo-relative paths like
#: ``REPO / "library/l1/tag-site-corpus-audit"`` remain remapped.
_RESTORE_S3_SUBS: tuple[tuple[str, str], ...] = (
    (
        "processed/scratch/becca/tag-site-corpus-audit",
        "processed/scratch/becca/tag-site-corpus-audit",
    ),
    (
        "processed://scratch/becca/tag-site-corpus-audit",
        "processed://scratch/becca/tag-site-corpus-audit",
    ),
    ('= "scratch/becca/tag-site-corpus-audit/', '= "scratch/becca/tag-site-corpus-audit/'),
    (
        '(processed, "scratch/becca/tag-site-corpus-audit/',
        '(processed, "scratch/becca/tag-site-corpus-audit/',
    ),
)
_REWRITE_SUFFIXES = (".py", ".md", ".txt", ".json")

_S4E9 = "library/l1/orf-design"
_S5E13 = "library/l1/tag-site-benchmarking"
_S1E3 = "library/l1/surface-filters"
#: Every historical ``(commit, repo_rel_path)`` file read by exported scripts via
#: ``deliverome_analysis.provenance.read_pinned_{bytes,text}``. Materialized into
#: ``historical/git-snapshots/<commit>/<path>`` with ``SHA256SUMS.json`` so exported
#: repos and Zenodo archives replay offline with zero ``.git`` history.
PINNED_GIT_SNAPSHOTS: tuple[tuple[str, str], ...] = (
    # synthesis_compliance.py & audit_domestication_barcodes.py (2026-11-25-pre-freeze-d)
    ("4422a9c4c6e2a5569e8a38a7190ef12046e5ddf2", f"{_S4E9}/constructs_l1_nterm.fasta"),
    ("4422a9c4c6e2a5569e8a38a7190ef12046e5ddf2", f"{_S4E9}/constructs_l1_nterm.tsv"),
    ("4422a9c4c6e2a5569e8a38a7190ef12046e5ddf2", f"{_S4E9}/constructs_l1_multipass.fasta"),
    ("4422a9c4c6e2a5569e8a38a7190ef12046e5ddf2", f"{_S4E9}/constructs_l1_multipass.tsv"),
    ("4422a9c4c6e2a5569e8a38a7190ef12046e5ddf2", f"{_S4E9}/constructs_l1_type3.fasta"),
    ("4422a9c4c6e2a5569e8a38a7190ef12046e5ddf2", f"{_S4E9}/constructs_l1_type3.tsv"),
    ("4422a9c4c6e2a5569e8a38a7190ef12046e5ddf2", f"{_S4E9}/constructs_l1_internal.fasta"),
    ("4422a9c4c6e2a5569e8a38a7190ef12046e5ddf2", f"{_S4E9}/constructs_l1_internal.tsv"),
    ("4422a9c4c6e2a5569e8a38a7190ef12046e5ddf2", f"{_S4E9}/constructs_l1_cterm.fasta"),
    ("4422a9c4c6e2a5569e8a38a7190ef12046e5ddf2", f"{_S4E9}/constructs_l1_cterm.tsv"),
    ("4422a9c4c6e2a5569e8a38a7190ef12046e5ddf2", f"{_S4E9}/candidate_manifest_tranches.tsv"),
    ("4422a9c4c6e2a5569e8a38a7190ef12046e5ddf2", f"{_S4E9}/l1_cassette_encodings.tsv"),
    ("599ba2483347d9c80f4e3dc8f801c337f727dfd7", f"{_S4E9}/synthesis_recodes.tsv"),
    # reviewer-package/fold_order_check.py
    ("4422a9c4", f"{_S4E9}/tag_refold_structures.tsv"),
    # Prior cut snapshots retained for historical comparison & exporter smoke tests
    ("8432c0104c8fdaa57e78173d18fc5bb7064d91bd", f"{_S4E9}/constructs_l1_cterm.tsv"),
    ("8432c010", f"{_S4E9}/tag_refold_structures.tsv"),
    # s1e3/audit_coreceptors.py
    ("6b3d53ccab9b95246a0605b79427d69dca36aeb2", f"{_S1E3}/final_library.py"),
    # membrane-side-gate/gate_inputs.py & gate_on_ordered_sites.py
    ("6bf341f231e0841ffea13e589c4119a909fde93f", f"{_S4E9}/l1_tag_paths.tsv"),
    ("6bf341f231e0841ffea13e589c4119a909fde93f", f"{_S4E9}/candidate_manifest_tranches.tsv"),
    ("6bf341f231e0841ffea13e589c4119a909fde93f", f"{_S4E9}/library_structure_payload.tsv"),
    ("34910ec68217eff67e10a8299524ae93bb210b8a", f"{_S4E9}/l1_tag_paths.tsv"),
    ("34910ec68217eff67e10a8299524ae93bb210b8a", f"{_S4E9}/candidate_manifest_tranches.tsv"),
    ("34910ec68217eff67e10a8299524ae93bb210b8a", f"{_S4E9}/library_structure_payload.tsv"),
    ("debbc418", f"{_S5E13}/l1_internal_decisions.tsv"),
    # multipass_rung_changes.py
    ("fbd83ab2", f"{_S4E9}/nterm_res1_multipass.tsv"),
    ("fbd83ab2", f"{_S4E9}/constructs_l1_multipass.tsv"),
    # tag_refold/build_refold_inputs.py, native_models.py & verify_ordered_constructs.py
    ("5b97ddfc557de0a1c78737b100daf16a1f734202", f"{_S4E9}/candidate_manifest_tranches.tsv"),
    ("5b133d0b", f"{_S4E9}/candidate_manifest_tranches.tsv"),
    ("5b133d0b", f"{_S4E9}/l1_tag_paths.tsv"),
    ("5b133d0b", f"{_S4E9}/library_structure_payload.tsv"),
    ("3ef7f910aff5510c452632b378a2985992436cb1", f"{_S4E9}/constructs_l1_cterm.fasta"),
    ("3ef7f910aff5510c452632b378a2985992436cb1", f"{_S4E9}/constructs_l1_cterm.tsv"),
    ("59050daa", f"{_S4E9}/l1_cassette_encodings.tsv"),
    ("8d4b935e22b10357b6cb45cacd615fc4c8a3c71e", f"{_S4E9}/synthesis_recodes.tsv"),
    # s5e13 internal-site scripts
    ("798d8717", f"{_S5E13}/l1_internal_arm_decisions.tsv"),
    ("89226cb2", f"{_S5E13}/l1_internal_arm_decisions.tsv"),
    ("74492351", f"{_S4E9}/l1_tag_paths.tsv"),
    ("74492351", f"{_S4E9}/l1_member_renumbering.tsv"),
    ("dc0ab0fc", f"{_S5E13}/tag_sites_v5_orthologs.tsv"),
    ("31b3ed6f", f"{_S5E13}/tag_sites_v5_orthologs.tsv"),
    ("9a1c2c91", f"{_S4E9}/l1_tag_paths.tsv"),
)


def git(*args: str) -> bytes:
    return subprocess.run(["git", "-C", str(REPO), *args], check=True, capture_output=True).stdout


def git_cat_batch(specs: list[str]) -> list[bytes]:
    """Fetch many ``<rev>:<path>`` blobs in a single ``git cat-file --batch`` call."""
    if not specs:
        return []
    proc = subprocess.run(
        ["git", "-C", str(REPO), "cat-file", "--batch"],
        input=("\n".join(specs) + "\n").encode("utf-8"),
        check=True,
        capture_output=True,
    )
    buf = memoryview(proc.stdout)
    pos = 0
    out: list[bytes] = []
    for spec in specs:
        nl = proc.stdout.index(b"\n", pos)
        header = bytes(buf[pos:nl]).decode("utf-8")
        pos = nl + 1
        parts = header.split()
        if len(parts) < 3 or parts[1] == "missing":
            raise RuntimeError(f"git cat-file --batch missing object: {spec!r} ({header})")
        size = int(parts[2])
        out.append(bytes(buf[pos : pos + size]))
        pos += size + 1  # trailing LF after each object
    return out


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def trace_to_list(log: Path) -> int:
    """Turn an audit-hook trace (one absolute path per line) into review_repo_files.txt."""
    tracked = set(git("ls-files").decode().splitlines())
    seen = set()
    for line in log.read_text().splitlines():
        try:
            rel = Path(line).resolve().relative_to(REPO.resolve()).as_posix()
        except ValueError:
            continue
        if rel in tracked and not rel.startswith(EXCLUDE_PREFIXES):
            seen.add(rel)
    # Package modules follow from the imports (``package_modules``), so are not listed here.
    listed = sorted(p for p in seen if not p.startswith(PKG + "/"))
    FILES.write_text("\n".join(listed) + "\n")
    print(f"wrote {FILES.name}: {len(listed)} traced files outside the package")
    return 0


def historical_commits() -> list[str]:
    sys.path.insert(0, str(S4E9))
    import reproduce_isoform_universe as riu

    return sorted(set(riu.COMMITS.values()))


def archive(ref: str, paths: list[str]) -> dict[str, bytes]:
    """``paths`` (files or folders) at ``ref`` -> {path: bytes}, via ``git archive``."""
    blob = git("archive", "--format=tar", ref, "--", *sorted(set(paths)))
    out: dict[str, bytes] = {}
    with tarfile.open(fileobj=io.BytesIO(blob)) as tar:
        for m in tar.getmembers():
            if m.isfile() and not m.name.startswith(EXCLUDE_PREFIXES):
                out[m.name] = tar.extractfile(m).read()
    return out


def _imported_modules(source: bytes, in_package: bool) -> set[str]:
    """Top-level ``deliverome_analysis`` modules a Python file imports, at any depth."""
    try:
        tree = ast.parse(source)
    except SyntaxError:
        return set()
    found: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for a in node.names:
                if a.name.startswith("deliverome_analysis."):
                    found.add(a.name.split(".")[1])
        elif isinstance(node, ast.ImportFrom):
            if node.level == 0 and node.module == "deliverome_analysis":
                found |= {a.name for a in node.names}
            elif node.level == 0 and (node.module or "").startswith("deliverome_analysis."):
                found.add(node.module.split(".")[1])
            elif node.level == 1 and in_package:
                found |= (
                    {node.module.split(".")[0]} if node.module else {a.name for a in node.names}
                )
    return found


def package_modules(ref: str, files: dict[str, bytes]) -> dict[str, bytes]:
    """The package files the exported code needs: ``__init__.py`` and every module reachable
    by import from ``files`` (the non-package exports), followed through the package itself.

    Found by parsing imports rather than listed by hand, so a new import in an L1 script
    brings its module along and a module L1 stops using drops out. A module named by a
    string (``importlib``) would be missed; nothing in the exported code does that.
    """
    pkg = archive(ref, [PKG])
    by_name: dict[str, list[str]] = {}
    for path in pkg:
        name = path[len(PKG) + 1 :].split("/")[0].removesuffix(".py")
        by_name.setdefault(name, []).append(path)
    todo = set()
    for path, data in files.items():
        if path.endswith(".py") and b"deliverome_analysis" in data:
            todo |= _imported_modules(data, in_package=False)
    needed: set[str] = set()
    while todo:
        name = todo.pop()
        if name in needed or name not in by_name:
            continue
        needed.add(name)
        for path in by_name[name]:
            if path.endswith(".py"):
                todo |= _imported_modules(pkg[path], in_package=True)
    keep = {"__init__"} | needed
    return {p: d for name in sorted(keep) for p, d in ((p, pkg[p]) for p in by_name.get(name, []))}


def _imported_siblings(path: str, source: bytes, tracked: set[str]) -> set[str]:
    """Same-directory ``.py`` modules imported by ``path`` that exist in ``tracked``."""
    try:
        tree = ast.parse(source)
    except SyntaxError:
        return set()
    parent = Path(path).parent
    found: set[str] = set()
    for node in ast.walk(tree):
        names: list[str] = []
        if isinstance(node, ast.Import):
            names = [a.name.split(".")[0] for a in node.names]
        elif isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
            names = [node.module.split(".")[0]]
        for name in names:
            cand = (parent / f"{name}.py").as_posix()
            if cand != path and cand in tracked and not cand.startswith(EXCLUDE_PREFIXES):
                found.add(cand)
    return found


def sibling_modules(ref: str, files: dict[str, bytes]) -> dict[str, bytes]:
    """Same-directory ``.py`` helpers transitively imported by any exported script at ``ref``.

    When ``ALWAYS`` names individual scripts from another experiment (e.g. ``s5e13``'s
    internal-site ranker and arm builder) rather than the whole folder, those scripts
    import local sibling modules via ``sys.path.insert(0, str(HERE))``. Following those
    imports keeps the exported scripts runnable in the review repository without exporting
    unrelated data tables from that folder.
    """
    tracked = set(git("ls-tree", "-r", "--name-only", ref).decode().splitlines())
    out: dict[str, bytes] = {}
    seen = set(files)
    # Only directories that still have unexported tracked .py files can yield new siblings.
    unexported_dirs = {
        Path(p).parent
        for p in tracked
        if p.endswith(".py")
        and p not in seen
        and not p.startswith(PKG + "/")
        and not p.startswith(EXCLUDE_PREFIXES)
    }
    todo = {
        p
        for p in files
        if p.endswith(".py") and not p.startswith(PKG + "/") and Path(p).parent in unexported_dirs
    }
    while todo:
        batch: set[str] = set()
        for path in sorted(todo):
            data = files.get(path) if path in files else out[path]
            for sib in _imported_siblings(path, data, tracked):
                if sib not in seen:
                    batch.add(sib)
        if not batch:
            break
        ordered = sorted(batch)
        blobs = git_cat_batch([f"{ref}:{p}" for p in ordered])
        fetched = dict(zip(ordered, blobs, strict=True))
        out.update(fetched)
        seen |= set(fetched)
        todo = set(fetched)
    return out


def package_paths(ref: str = "HEAD") -> set[str]:
    """The package paths an export of ``ref`` publishes (for ``render_brief``'s link check)."""
    wanted = [p for p in FILES.read_text().split() if p] + ALWAYS
    out = archive(ref, wanted)
    out.update(sibling_modules(ref, out))
    return set(package_modules(ref, out))


def _remap_path(rel: str) -> str:
    for old, new in CLEAN_PATH_MAP.items():
        if rel == old or rel.startswith(old + "/"):
            rel = new + rel[len(old) :]
            break
    for old, new in _CLEAN_SUBDIR_MAP.items():
        rel = rel.replace(f"/{old}/", f"/{new}/")
    return rel


_CLEAN_SCRATCH_EXCLUDES: tuple[str, ...] = tuple(
    sorted(new for old, new in CLEAN_PATH_MAP.items() if old.startswith("scratch/"))
)


_DEAD_MYPY_OVERRIDES: tuple[str, ...] = (
    "[[tool.mypy.overrides]]\n"
    'module = ["deliverome_analysis.peptide_metrics", "deliverome_analysis.alphadia_trim"]\n'
    "strict = true\n\n",
    '[[tool.mypy.overrides]]\nmodule = "deliverome_analysis.raw_provenance"\nstrict = true\n\n',
    "[[tool.mypy.overrides]]\n"
    "# QC metric modules: chromatography (TIC, spray, LC), surface topology, and the pipeline\n"
    "module = [\n"
    '    "deliverome_analysis.raw_chromatography",\n'
    '    "deliverome_analysis.surface_topology",\n'
    '    "deliverome_analysis.dia_qc_pipeline",\n'
    "]\n"
    "strict = true\n\n",
)


def _rewrite_pyproject_clean_paths(text: str) -> str:
    """Update ruff excludes, pytest testpaths, and mypy overrides in ``pyproject.toml``
    for clean ``library/l1/*`` exports without invalidating ``uv.lock``."""
    excludes = [
        '"scratch"',
        '"analysis/*/data"',
        '"library/l1/*/data"',
        *(f'"{p}"' for p in _CLEAN_SCRATCH_EXCLUDES),
    ]
    replacement = "extend-exclude = [\n" + "".join(f"    {e},\n" for e in excludes) + "]"
    text = text.replace('extend-exclude = ["scratch", "analysis/*/data"]', replacement)
    text = text.replace(
        'testpaths = ["tests", "analysis"]',
        'testpaths = ["tests", "library/l1/orf-design/tests"]',
    )
    # Strip mypy overrides for non-exported lab modules (safe: uv.lock ignores [tool.mypy];
    # [tool.uv.sources] and dependency-groups stay intact so uv sync --locked / --frozen passes).
    for block in _DEAD_MYPY_OVERRIDES:
        text = text.replace(block, "")
    return text


def _ruff_format_bytes(path: str, data: bytes) -> bytes:
    """Re-format a clean-path-rewritten Python file with ruff if available so CI format checks pass."""
    ruff = shutil.which("ruff")
    if not ruff:
        return data
    proc = subprocess.run(
        [ruff, "format", "--config", str(REPO / "pyproject.toml"), "--stdin-filename", path, "-"],
        input=data,
        capture_output=True,
        check=False,
    )
    return proc.stdout if proc.returncode == 0 else data


def apply_clean_paths(files: dict[str, bytes]) -> dict[str, bytes]:
    """Remap internal ``analysis/sprintN/sNeM-*`` and ``scratch/*`` paths to clean ``library/l1/*``
    paths, rewriting repo-relative references inside text files while leaving S3 keys intact."""
    out: dict[str, bytes] = {}
    scratch_prefixes = tuple(f"{p}/" for p in _CLEAN_SCRATCH_EXCLUDES)
    for rel, data in files.items():
        dest = _remap_path(rel)
        if dest == "pyproject.toml":
            text = data.decode("utf-8", errors="surrogateescape")
            data = _rewrite_pyproject_clean_paths(text).encode("utf-8", errors="surrogateescape")
        elif dest.endswith(_REWRITE_SUFFIXES):
            orig_text = data.decode("utf-8", errors="surrogateescape")
            text = orig_text
            for old, new in CLEAN_PATH_MAP.items():
                text = text.replace(old, new)
            if dest.endswith(".py"):
                for old, new in _CLEAN_REL_SUBS:
                    text = text.replace(old, new)
                for old, new in _RESTORE_S3_SUBS:
                    text = text.replace(old, new)
            elif dest.endswith(".md"):
                for old, new in _CLEAN_MD_REL_SUBS:
                    text = text.replace(old, new)
            data = text.encode("utf-8", errors="surrogateescape")
            if dest.endswith(".py") and text != orig_text and not dest.startswith(scratch_prefixes):
                data = _ruff_format_bytes(dest, data)
        out[dest] = data
    return out


def pinned_git_snapshots(*, clean_paths: bool = False) -> dict[str, bytes]:
    """Materialize every pinned ``(commit, rel_path)`` file under ``historical/git-snapshots/``
    with a ``SHA256SUMS.json`` manifest so exported repos never need git history."""
    out: dict[str, bytes] = {}
    manifest: dict[str, str] = {}
    specs = [f"{commit}:{rel}" for commit, rel in PINNED_GIT_SNAPSHOTS]
    blobs = git_cat_batch(specs)
    for (commit, rel), data in zip(PINNED_GIT_SNAPSHOTS, blobs, strict=True):
        dest_rel = _remap_path(rel) if clean_paths else rel
        key = f"{commit}/{dest_rel}"
        out[f"historical/git-snapshots/{key}"] = data
        manifest[key] = sha256(data)
    out["historical/git-snapshots/SHA256SUMS.json"] = (
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    ).encode()
    return out


def _remap_tar_paths(blob: bytes) -> bytes:
    """Remap file paths inside ``historical/<commit>.tar`` using :func:`_remap_path` so
    ``run_stage.py`` resolves ``dest / riu.AP`` in a clean-paths exported review repo."""
    buf = io.BytesIO()
    with (
        tarfile.open(fileobj=io.BytesIO(blob)) as src,
        tarfile.open(fileobj=buf, mode="w", format=tarfile.PAX_FORMAT) as dst,
    ):
        for m in src.getmembers():
            if not m.isfile():
                continue
            fh = src.extractfile(m)
            assert fh is not None
            data = fh.read()
            info = tarfile.TarInfo(_remap_path(m.name))
            info.size = len(data)
            info.mtime = 0
            info.mode = m.mode
            dst.addfile(info, io.BytesIO(data))
    return buf.getvalue()


def collect(
    ref: str, *, clean_paths: bool = False, include_historical: bool = True
) -> dict[str, bytes]:
    """Every managed path -> its bytes at ``ref``."""
    wanted = [p for p in FILES.read_text().split() if p] + ALWAYS
    out = archive(ref, wanted)
    out.update(sibling_modules(ref, out))
    out.update(package_modules(ref, out))
    if include_historical:
        commits = historical_commits()
        import reproduce_isoform_universe as riu

        for commit in commits:
            tar_blob = git(
                "archive",
                commit,
                riu.AP,
                "library/l1/tag-site-architecture/library_members.json",
            )
            if clean_paths:
                tar_blob = _remap_tar_paths(tar_blob)
            out[f"historical/{commit}.tar"] = tar_blob
    gen_items = list(GENERATED.items())
    gen_specs = [f"{ref}:{(TEMPLATES / tpl).relative_to(REPO).as_posix()}" for _, tpl in gen_items]
    for (dest, _), data in zip(gen_items, git_cat_batch(gen_specs), strict=True):
        out[dest] = data
    if clean_paths:
        out = apply_clean_paths(out)
    if include_historical:
        out.update(pinned_git_snapshots(clean_paths=clean_paths))
    return out


def bundle_tables(bundle: Path, commit: str) -> dict[str, bytes]:
    """The small tables ``data_locations.tsv`` (at ``commit``) assigns to the review repository.

    They come from the frozen bundle, not from git: they are data, kept out of this
    repository's history, and each is checked against its recorded SHA-256 before it is copied.
    """
    import csv
    import io as _io

    rel = (HERE / "data_locations.tsv").relative_to(REPO).as_posix()
    rows = csv.DictReader(_io.StringIO(git("show", f"{commit}:{rel}").decode()), delimiter="\t")
    out: dict[str, bytes] = {}
    for r in rows:
        if r["home"] != "review_repo":
            continue
        data = (bundle / "s3" / r["bucket"] / r["key"]).read_bytes()
        if sha256(data) != r["sha256"]:
            raise SystemExit(f"{r['bucket']}/{r['key']}: SHA-256 differs from data_locations.tsv")
        out[r["location"]] = data
    return out


def export(
    dest: Path,
    ref: str,
    overwrite_edited: bool,
    bundle: Path | None = None,
    *,
    clean_paths: bool = True,
    include_historical: bool = True,
) -> int:
    commit = git("rev-parse", ref).decode().strip()
    files = collect(commit, clean_paths=clean_paths, include_historical=include_historical)
    if bundle is not None:
        files.update(bundle_tables(bundle, commit))
    old = {}
    if (dest / MANIFEST).exists():
        old = json.loads((dest / MANIFEST).read_text())["files"]

    edited = []
    for rel, digest in old.items():
        path = dest / rel
        if path.exists() and sha256(path.read_bytes()) != digest:
            edited.append(rel)
    written = removed = 0
    for rel, data in sorted(files.items()):
        path = dest / rel
        if rel in edited and not overwrite_edited:
            continue
        if path.exists() and rel not in old and path.read_bytes() != data:
            edited.append(rel)  # a file the reviewer added under a now-managed path
            if not overwrite_edited:
                continue
        path.parent.mkdir(parents=True, exist_ok=True)
        if not path.exists() or path.read_bytes() != data:
            path.write_bytes(data)
            written += 1
    for rel in sorted(set(old) - set(files)):
        if rel not in edited and (dest / rel).exists():
            (dest / rel).unlink()
            removed += 1

    kept = set(edited) if not overwrite_edited else set()
    manifest = {
        "source_repository": "Deliverome-Project/deliverome-analysis",
        "source_commit": commit,
        "source_ref": ref,
        "files": {
            rel: (old[rel] if rel in kept and rel in old else sha256(data))
            for rel, data in sorted(files.items())
            if not (rel in kept and rel not in old)
        },
    }
    (dest / MANIFEST).write_text(json.dumps(manifest, indent=1) + "\n")
    mods = sorted(Path(r).stem for r in files if r.startswith(PKG + "/") and r.endswith(".py"))
    print(f"package modules exported ({len(mods)}): {', '.join(mods)}")
    print(
        f"exported {len(files)} managed files from {commit[:10]} to {dest}: "
        f"{written} written, {removed} removed"
    )
    if edited:
        verb = "overwritten" if overwrite_edited else "left as edited (not updated)"
        print(f"{len(edited)} file(s) edited in the review repo, {verb}:")
        for rel in edited:
            print(f"  {rel}")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--dest", type=Path)
    ap.add_argument("--ref", default="HEAD", help="commit or branch to export (default HEAD)")
    ap.add_argument(
        "--bundle",
        type=Path,
        help="the collected bundle folder; its small tables are copied to data/ (public_data.py)",
    )
    ap.add_argument("--overwrite-edited", action="store_true")
    ap.add_argument(
        "--keep-internal-paths",
        action="store_true",
        help="keep analysis/sprintN/sNeM-* and scratch/* paths instead of remapping to library/l1/*",
    )
    ap.add_argument("--trace", type=Path, help="rebuild review_repo_files.txt from a trace log")
    a = ap.parse_args(argv)
    if a.trace:
        return trace_to_list(a.trace)
    if not a.dest:
        ap.error("--dest is required")
    bundle = a.bundle.expanduser() if a.bundle else None
    return export(
        a.dest.expanduser(),
        a.ref,
        a.overwrite_edited,
        bundle,
        clean_paths=not a.keep_internal_paths,
    )


if __name__ == "__main__":
    raise SystemExit(main())
