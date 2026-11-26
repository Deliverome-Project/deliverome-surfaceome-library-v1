"""Every file of the reviewer bundle has exactly one home a reviewer can reach.

``reviewer-package/data_locations.tsv`` assigns each row of ``bundle_manifest.tsv`` to the
public R2 bucket, the review repository's ``data/`` folder, or the private release
(``public_data.py``). These tests fail if a file is missing, duplicated, or recorded with a
different size or SHA-256 than the manifest, so the complete dataset cannot silently shrink.

The review-repository check runs where ``data/`` exists (the exported repository); the online
check of the public bucket runs only with ``L1_CHECK_PUBLIC_DATA=1``, since it needs network.
"""

from __future__ import annotations

import csv
import hashlib
import os
import sys
from pathlib import Path

import pytest

RP = Path(__file__).resolve().parents[1] / "reviewer-package"
REPO = Path(__file__).resolve().parents[4]
HOMES = {"public_r2", "public_source", "review_repo", "private_release"}


def _rows(name: str) -> list[dict[str, str]]:
    with (RP / name).open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


MANIFEST = _rows("bundle_manifest.tsv")
LOCATIONS = _rows("data_locations.tsv")


def _key(r: dict[str, str]) -> tuple[str, str]:
    return r["bucket"], r["key"]


def test_every_manifest_file_has_exactly_one_home() -> None:
    placed = [_key(r) for r in LOCATIONS]
    dupes = {k for k in placed if placed.count(k) > 1}
    assert not dupes, f"files with more than one home: {sorted(dupes)[:5]}"
    missing = {_key(r) for r in MANIFEST} - set(placed)
    extra = set(placed) - {_key(r) for r in MANIFEST}
    assert not missing, f"{len(missing)} bundle files have no home, e.g. {sorted(missing)[:5]}"
    assert not extra, f"{len(extra)} homes for files not in the bundle, e.g. {sorted(extra)[:5]}"


def test_sizes_and_hashes_match_the_manifest() -> None:
    by_key = {_key(r): r for r in MANIFEST}
    bad = [
        _key(r)
        for r in LOCATIONS
        if (r["bytes"], r["sha256"]) != (by_key[_key(r)]["bytes"], by_key[_key(r)]["sha256"])
    ]
    assert not bad, f"{len(bad)} rows disagree with bundle_manifest.tsv, e.g. {bad[:5]}"


def test_homes_follow_the_rules() -> None:
    sys.path.insert(0, str(RP))
    import public_data

    for r in LOCATIONS:
        assert r["home"] in HOMES, r
        home, license_ = public_data.home_of(r)
        assert (r["home"], r["license"]) == (home, license_), (
            f"{_key(r)}: rerun public_data.py plan"
        )
        if r["home"] in ("public_r2", "public_source"):
            assert r["license"], f"{_key(r)}: a public file needs a stated license"
            assert r["location"].startswith("https://"), f"{_key(r)}: a public file needs a URL"
        if r["home"] == "review_repo":
            assert int(r["bytes"]) <= public_data.REPO_MAX_BYTES


@pytest.mark.skipif(not (REPO / "data").is_dir(), reason="runs in the review repository")
def test_review_repo_tables_are_present_and_intact() -> None:
    bad = []
    for r in LOCATIONS:
        if r["home"] != "review_repo":
            continue
        path = REPO / r["location"]
        if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != r["sha256"]:
            bad.append(r["location"])
    assert not bad, f"{len(bad)} committed tables missing or changed, e.g. {bad[:5]}"


def test_every_active_afdb_model_is_archived_and_bundled() -> None:
    """Every active AlphaFold DB model in library_manifest_l1.tsv must be archived and bundled."""
    s4e9 = RP.parent
    arch_rows = list(csv.DictReader((s4e9 / "tag_refold_afdb_archive.tsv").open(), delimiter="\t"))
    archived = {r["archive_key"].removeprefix("external://") for r in arch_rows}
    bundled = {_key(r) for r in MANIFEST}
    placed = {_key(r) for r in LOCATIONS}

    mf_rows = list(csv.DictReader((s4e9 / "library_manifest_l1.tsv").open(), delimiter="\t"))
    active_accs = sorted(
        {
            r["structure_acc"]
            for r in mf_rows
            if r["structure_source"] == "alphafold_db" and r["structure_acc"]
        }
    )
    missing_arch, missing_bundle, missing_loc = [], [], []
    for acc in active_accs:
        key = f"alphafold-db/v6/AF-{acc}-F1-model_v6.pdb"
        if key not in archived:
            missing_arch.append(acc)
        if ("external", key) not in bundled:
            missing_bundle.append(acc)
        if ("external", key) not in placed:
            missing_loc.append(acc)
    assert not missing_arch, (
        f"{len(missing_arch)} active AFDB models not archived: {missing_arch[:5]}"
    )
    assert not missing_bundle, (
        f"{len(missing_bundle)} active AFDB models missing from bundle_manifest.tsv: {missing_bundle[:5]}"
    )
    assert not missing_loc, (
        f"{len(missing_loc)} active AFDB models missing from data_locations.tsv: {missing_loc[:5]}"
    )


@pytest.mark.skipif(
    os.environ.get("L1_CHECK_PUBLIC_DATA") != "1", reason="set L1_CHECK_PUBLIC_DATA=1 (network)"
)
def test_every_public_object_is_online() -> None:
    sys.path.insert(0, str(RP))
    import public_data

    assert public_data.check(workers=4) == 0
