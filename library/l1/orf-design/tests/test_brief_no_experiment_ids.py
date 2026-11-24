"""The reviewer documents never show an internal experiment ID (``s4e9``, ``S5E13``, ...).

Becca, 2026-11-22: "i told you not to use internal snen ids in this document". Reviewers do not
know our sprint/experiment numbering, so an ID in the text is noise at best. Link targets may still
point into the repository's experiment folders (a link has to resolve), so only the text a reader
sees is checked: prose, link text and code spans, with every Markdown link target removed. Fenced
code blocks are commands a reviewer runs, so their paths must stay exact and are not checked.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

RP = Path(__file__).resolve().parents[1] / "reviewer-package"
DOCS = ["BRIEF.md", "COVER_NOTE.md", "DESIGN.md", "nterm/BRIEF.md"]
EXPERIMENT_ID = re.compile(r"\bs\d+e\d+\b", re.IGNORECASE)
LINK_TARGET = re.compile(r"\]\([^)]*\)")
FENCED = re.compile(r"^```.*?^```", re.MULTILINE | re.DOTALL)
S3_URI = re.compile(r"`(?:processed|external|raw)://[^`]+`")


def visible_text(md: str) -> str:
    return S3_URI.sub("", LINK_TARGET.sub("]", FENCED.sub("", md)))


@pytest.mark.parametrize("name", DOCS)
def test_no_experiment_id_in_visible_text(name: str) -> None:
    text = visible_text((RP / name).read_text())
    hits = [
        f"line {i}: {line.strip()[:100]}"
        for i, line in enumerate(text.splitlines(), 1)
        if EXPERIMENT_ID.search(line)
    ]
    assert not hits, f"{name} shows internal experiment IDs:\n" + "\n".join(hits[:20])


def test_link_targets_are_not_checked() -> None:
    assert not EXPERIMENT_ID.search(visible_text("[the filter](../../sprint1/s1e3-x/f.py)"))
