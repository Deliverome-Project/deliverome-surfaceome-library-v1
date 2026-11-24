"""Tests that changing any L1 protein sequence or tag location marks the review out of date."""

from __future__ import annotations

import json
import sys
from pathlib import Path

S4E9 = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(S4E9 / "reviewer-package"))

import review_fingerprint as rf  # noqa: E402


def _write_tsv(path: Path, header: list[str], rows: list[list[str]]) -> None:
    lines = ["\t".join(header), *("\t".join(r) for r in rows)]
    path.write_text("\n".join(lines) + "\n")


def _seed_mini_library(root: Path) -> Path:
    _write_tsv(
        root / "l1_tag_paths.tsv",
        ["candidate_id", "placement", "lane_override", "length_aa", "signal_peptide_aa"],
        [
            ["l1-p0001", "internal", "", "400", "20"],
            ["l1-p0002", "C-terminal", "", "300", "0"],
        ],
    )
    _write_tsv(
        root / "constructs_l1_internal.tsv",
        [
            "construct_id",
            "design_id",
            "placement",
            "cassette_after_residue",
            "native_aa",
            "protein_aa",
        ],
        [["l1-p0001-a1", "l1-p0001", "internal", "142", "MAAA", "MAATAGAA"]],
    )
    _write_tsv(
        root / "constructs_l1_cterm.tsv",
        [
            "candidate_id",
            "design_id",
            "placement",
            "cassette_after_residue",
            "native_aa",
            "protein_aa",
        ],
        [["l1-p0002", "l1-p0002", "C-terminal", "300", "MCCC", "MCCCTAG"]],
    )
    brief = root / "brief_numbers.json"
    fp = rf.compute(root)
    brief.write_text(
        json.dumps({"tag_paths": {"placement": fp["placements"]}, "review_fingerprint": fp})
    )
    return brief


def test_review_fingerprint_passes_when_unchanged(tmp_path: Path) -> None:
    brief = _seed_mini_library(tmp_path)
    assert rf.check_brief_freshness(tmp_path, brief) == []


def test_review_fingerprint_flags_protein_sequence_change(tmp_path: Path) -> None:
    brief = _seed_mini_library(tmp_path)
    _write_tsv(
        tmp_path / "constructs_l1_internal.tsv",
        [
            "construct_id",
            "design_id",
            "placement",
            "cassette_after_residue",
            "native_aa",
            "protein_aa",
        ],
        [["l1-p0001-a1", "l1-p0001", "internal", "142", "MAAX", "MAATAGAX"]],
    )
    problems = rf.check_brief_freshness(tmp_path, brief)
    assert len(problems) == 1
    assert "l1-p0001" in problems[0]


def test_review_fingerprint_flags_tag_location_change(tmp_path: Path) -> None:
    brief = _seed_mini_library(tmp_path)
    _write_tsv(
        tmp_path / "constructs_l1_internal.tsv",
        [
            "construct_id",
            "design_id",
            "placement",
            "cassette_after_residue",
            "native_aa",
            "protein_aa",
        ],
        [["l1-p0001-a1", "l1-p0001", "internal", "210", "MAAA", "MAATAGAA"]],
    )
    problems = rf.check_brief_freshness(tmp_path, brief)
    assert len(problems) == 1
    assert "l1-p0001" in problems[0]


def test_review_fingerprint_flags_placement_count_drift(tmp_path: Path) -> None:
    brief = _seed_mini_library(tmp_path)
    _write_tsv(
        tmp_path / "l1_tag_paths.tsv",
        ["candidate_id", "placement", "lane_override", "length_aa", "signal_peptide_aa"],
        [
            ["l1-p0001", "internal", "", "400", "20"],
            ["l1-p0002", "C-terminal", "", "300", "0"],
            ["l1-p0003", "N after SP", "", "250", "18"],
        ],
    )
    problems = rf.check_brief_freshness(tmp_path, brief)
    assert any("tag_paths.placement out of date" in p for p in problems)
    assert any("added=1 (l1-p0003)" in p for p in problems)


def test_brief_numbers_matches_live_tables() -> None:
    """Committed brief_numbers.json must match the live L1 tag-path and construct tables."""
    problems = rf.check_brief_freshness()
    assert not problems, (
        "Review brief is out of date — rerun brief_numbers.py:\n  - " + "\n  - ".join(problems)
    )
