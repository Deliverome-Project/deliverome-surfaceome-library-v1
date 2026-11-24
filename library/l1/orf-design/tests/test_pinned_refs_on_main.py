"""Every git commit an s4e9 script pins must be reachable from this checkout's history.

The check lives in ``pinned_refs.py`` (also the pipeline stage ``pinned_refs``). CI checks out
one commit, where ancestry cannot be decided; this test then skips and says so. It runs
wherever the full history is present.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pinned_refs

_IS_EXPORTED_REVIEW_REPO = (Path(__file__).resolve().parents[4] / "EXPORT_MANIFEST.json").is_file()


@pytest.mark.skipif(
    _IS_EXPORTED_REVIEW_REPO, reason="source-repo git history check (not in exported review repo)"
)
def test_pinned_refs_are_ancestors_of_head():
    if pinned_refs.is_shallow():
        pytest.skip("shallow clone: ancestry needs full history (git fetch --unshallow)")
    try:
        pinned_refs.check()
    except SystemExit as exc:
        pytest.fail(str(exc))


@pytest.mark.skipif(
    _IS_EXPORTED_REVIEW_REPO, reason="source-repo branch-ref check (not in exported review repo)"
)
def test_no_script_reads_a_session_branch_by_name():
    """Runs without history: a session branch named as a git ref is refused, except the
    entries in BRANCH_PENDING (owned by another session), which must still need fixing."""
    found = pinned_refs.branch_refs()
    new = sorted(set(found) - set(pinned_refs.BRANCH_PENDING))
    assert not new, f"scripts read a session branch by name: {new}"
    fixed = sorted(set(pinned_refs.BRANCH_PENDING) - set(found))
    assert not fixed, f"fixed: remove from pinned_refs.BRANCH_PENDING: {fixed}"
