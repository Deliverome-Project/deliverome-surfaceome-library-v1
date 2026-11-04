"""Export the L1 review repository from a commit of this one.

The review repo holds everything a reviewer needs to rerun the L1 design, and is also a
working repo: files can be added to it, and the shared utilities edited, as review needs.
So the export manages a known set of files and leaves everything else alone:

* **Managed files** are listed in ``review_repo_files.txt`` (every repository file the
  offline replay opens, found by tracing it; see ``--trace``) plus the whole
  ``deliverome_analysis`` package and the files that define the environment. They are
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

    uv run --locked python export_review_repo.py --dest ~/git/l1-library-review
    uv run --locked python export_review_repo.py --dest ... --ref claude/l1-reviewer-package-aa447f
    uv run --locked python export_review_repo.py --trace trace.log   # refresh the file list
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
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
#: folders (not only those the replay runs), the shared utilities, and the environment.
ALWAYS = [
    "analysis/sprint4/s4e9-orf-design-fg-panel",
    "analysis/sprint1/s1e3-surface-library-filters",
    "scratch/becca/appris-principal",
    "scratch/becca/tag-site-architecture",
    "src/deliverome_analysis",
    "pyproject.toml",
    "uv.lock",
    ".python-version",
    "tests/test_s3_offline_mirror.py",
]
#: Files this package writes itself, from ``review-repo/`` templates (path in dest -> template).
GENERATED = {"README.md": "README.md", ".gitignore": "gitignore"}
#: Never exported, even if a trace opened them.
EXCLUDE_PREFIXES = (".env", ".venv/", ".git/", "historical/")


def git(*args: str) -> bytes:
    return subprocess.run(["git", "-C", str(REPO), *args], check=True, capture_output=True).stdout


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
    # The package itself is exported whole, so it need not be listed file by file.
    listed = sorted(p for p in seen if not p.startswith("src/deliverome_analysis/"))
    FILES.write_text("\n".join(listed) + "\n")
    print(f"wrote {FILES.name}: {len(listed)} traced files outside the package")
    return 0


def historical_commits() -> list[str]:
    sys.path.insert(0, str(S4E9))
    import reproduce_isoform_universe as riu

    return sorted(set(riu.COMMITS.values()))


def collect(ref: str) -> dict[str, bytes]:
    """Every managed path -> its bytes at ``ref``."""
    wanted = [p for p in FILES.read_text().split() if p] + ALWAYS
    blob = git("archive", "--format=tar", ref, "--", *sorted(set(wanted)))
    out: dict[str, bytes] = {}
    with tarfile.open(fileobj=io.BytesIO(blob)) as tar:
        for m in tar.getmembers():
            if m.isfile() and not m.name.startswith(EXCLUDE_PREFIXES):
                out[m.name] = tar.extractfile(m).read()
    commits = historical_commits()
    import reproduce_isoform_universe as riu

    for commit in commits:
        out[f"historical/{commit}.tar"] = git(
            "archive",
            commit,
            riu.AP,
            "scratch/becca/tag-site-architecture/library_members.json",
        )
    for dest, template in GENERATED.items():
        out[dest] = git("show", f"{ref}:{(TEMPLATES / template).relative_to(REPO).as_posix()}")
    return out


def export(dest: Path, ref: str, overwrite_edited: bool) -> int:
    commit = git("rev-parse", ref).decode().strip()
    files = collect(commit)
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
    ap.add_argument("--overwrite-edited", action="store_true")
    ap.add_argument("--trace", type=Path, help="rebuild review_repo_files.txt from a trace log")
    a = ap.parse_args(argv)
    if a.trace:
        return trace_to_list(a.trace)
    if not a.dest:
        ap.error("--dest is required")
    return export(a.dest.expanduser(), a.ref, a.overwrite_edited)


if __name__ == "__main__":
    raise SystemExit(main())
