"""Fingerprint every L1 member's protein sequence and tag location for review freshness.

Computes a deterministic per-``design_id`` SHA-256 from the committed placement table
(``l1_tag_paths.tsv``) and all five built-arm tables (``constructs_l1_*.tsv``), plus an
overall library digest. Comparing a stamped fingerprint (stored in ``brief_numbers.json``)
against the live tables immediately flags any protein or tag-location change as out of date
and names the exact ``design_id``s affected.

    uv run --locked python library/l1/orf-design/reviewer-package/review_fingerprint.py --check
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
S4E9 = HERE.parent
BRIEF_NUMBERS = HERE / "brief_numbers.json"

ARM_TABLES: tuple[tuple[str, str], ...] = (
    ("nterm", "constructs_l1_nterm.tsv"),
    ("multipass", "constructs_l1_multipass.tsv"),
    ("type3", "constructs_l1_type3.tsv"),
    ("internal", "constructs_l1_internal.tsv"),
    ("cterm", "constructs_l1_cterm.tsv"),
)
TAG_PATHS = "l1_tag_paths.tsv"


def _read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _sha256_json(obj: Any) -> str:
    payload = json.dumps(obj, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def compute(s4e9_dir: Path = S4E9) -> dict[str, Any]:
    """Return placement counts, per-member fingerprints, and an overall library digest."""
    tag_rows = _read_tsv(s4e9_dir / TAG_PATHS)
    placements = dict(sorted(Counter(r.get("placement", "") for r in tag_rows).items()))

    by_member: dict[str, dict[str, Any]] = {}
    for r in tag_rows:
        cid = r.get("candidate_id", "")
        if not cid:
            continue
        by_member.setdefault(cid, {"tag_path": None, "arms": []})["tag_path"] = {
            "placement": r.get("placement", ""),
            "lane_override": r.get("lane_override", ""),
            "length_aa": r.get("length_aa", ""),
            "signal_peptide_aa": r.get("signal_peptide_aa", ""),
            "sequence_sha256": r.get("sequence_sha256", ""),
        }

    for rung, filename in ARM_TABLES:
        for r in _read_tsv(s4e9_dir / filename):
            did = r.get("design_id") or r.get("candidate_id") or ""
            if not did:
                continue
            arm_rec = {
                "rung": rung,
                "construct_id": r.get("construct_id") or r.get("arm") or "",
                "placement": r.get("placement") or r.get("design") or "",
                "cassette_after_residue": r.get("cassette_after_residue")
                or r.get("insert_after_residue")
                or "",
                "native_aa": r.get("native_aa", ""),
                "protein_aa": r.get("protein_aa", ""),
                "native_sha256": r.get("native_sha256", ""),
                "protein_sha256": r.get("protein_sha256", ""),
                "mature_sha256": r.get("mature_sha256", ""),
            }
            by_member.setdefault(did, {"tag_path": None, "arms": []})["arms"].append(arm_rec)

    member_hashes: dict[str, str] = {}
    for did in sorted(by_member):
        entry = by_member[did]
        entry["arms"].sort(
            key=lambda a: (
                a["rung"],
                a["construct_id"],
                a["placement"],
                a["cassette_after_residue"],
                a["protein_sha256"],
                a["native_sha256"],
                a["mature_sha256"],
                a["protein_aa"],
            )
        )
        member_hashes[did] = _sha256_json(entry)[:16]

    return {
        "n_tag_paths": len(tag_rows),
        "placements": placements,
        "n_members_fingerprinted": len(member_hashes),
        "library_sha256": _sha256_json(member_hashes),
        "members": member_hashes,
    }


def diff_fingerprints(recorded: dict[str, Any], live: dict[str, Any]) -> dict[str, list[str]]:
    """Return design_ids added, removed, or changed between two fingerprints."""
    old_m: dict[str, str] = recorded.get("members") or {}
    new_m: dict[str, str] = live.get("members") or {}
    old_keys, new_keys = set(old_m), set(new_m)
    return {
        "added": sorted(new_keys - old_keys),
        "removed": sorted(old_keys - new_keys),
        "changed": sorted(k for k in (old_keys & new_keys) if old_m[k] != new_m[k]),
    }


def check_brief_freshness(
    s4e9_dir: Path = S4E9, brief_numbers_path: Path = BRIEF_NUMBERS
) -> list[str]:
    """Return human-readable staleness problems between ``brief_numbers.json`` and live TSVs."""
    brief = json.loads(brief_numbers_path.read_text())
    live = compute(s4e9_dir)
    problems: list[str] = []

    recorded_placements = (brief.get("tag_paths") or {}).get("placement") or {}
    if recorded_placements != live["placements"]:
        problems.append(
            f"tag_paths.placement out of date: brief_numbers.json has {recorded_placements}, "
            f"live l1_tag_paths.tsv has {live['placements']}"
        )

    recorded_fp = brief.get("review_fingerprint")
    if isinstance(recorded_fp, dict) and recorded_fp.get("library_sha256"):
        if recorded_fp["library_sha256"] != live["library_sha256"]:
            delta = diff_fingerprints(recorded_fp, live)
            parts = [
                f"{k}={len(v)} ({', '.join(v[:5])}{'…' if len(v) > 5 else ''})"
                for k, v in delta.items()
                if v
            ]
            problems.append("protein/tag-site fingerprint out of date: " + "; ".join(parts))
    return problems


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--check", action="store_true", help="exit non-zero if brief_numbers.json is stale"
    )
    args = ap.parse_args(argv)
    if args.check:
        problems = check_brief_freshness()
        if problems:
            print("REVIEW OUT OF DATE:")
            for p in problems:
                print(f"  - {p}")
            return 1
        print("review is up to date with live L1 tables")
        return 0
    fp = compute()
    print(json.dumps({k: v for k, v in fp.items() if k != "members"}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
