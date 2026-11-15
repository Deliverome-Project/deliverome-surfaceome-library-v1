"""No committed s4e9 JSON record names a path on one person's machine.

A provenance record that points at ``/Users/...`` or a scratch folder cannot be followed by
anyone else, and leaks a home-directory layout. Records name files relative to the repository.
Added after the 2026-11-14 independent audit (item Q7) by the Claude Code session
l1-reviewer-package.
"""

from __future__ import annotations

import re
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
MACHINE = re.compile(r'"(?:/Users/|/home/|/private/tmp/|/var/folders/|~/)')


def test_json_records_carry_no_machine_paths():
    files = subprocess.run(
        ["git", "ls-files", "*.json"], cwd=HERE, capture_output=True, text=True, check=True
    ).stdout.split()
    bad = [f for f in files if MACHINE.search((HERE / f).read_text(errors="ignore"))]
    assert not bad, f"machine-local paths in committed records: {bad}"
