"""Tests for offline pinned git-snapshot reading in deliverome_analysis.provenance."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from deliverome_analysis import provenance as prov


def test_read_pinned_bytes_from_materialized_snapshot_without_git(tmp_path: Path) -> None:
    commit = "8432c010"
    rel = "library/l1/orf-design/example.tsv"
    payload = b"id\tval\nl1-p0001\t42\n"
    digest = hashlib.sha256(payload).hexdigest()

    snap_root = tmp_path / prov.PINNED_SNAPSHOT_DIR
    snap_file = snap_root / commit / rel
    snap_file.parent.mkdir(parents=True)
    snap_file.write_bytes(payload)
    (snap_root / prov.PINNED_MANIFEST_NAME).write_text(json.dumps({f"{commit}/{rel}": digest}))

    assert prov.read_pinned_bytes(commit, rel, repo=tmp_path) == payload
    assert prov.read_pinned_text(commit, rel, repo=tmp_path, expected_sha256=digest) == (
        payload.decode("utf-8")
    )


def test_read_pinned_bytes_rejects_tampered_snapshot(tmp_path: Path) -> None:
    commit = "8432c010"
    rel = "library/l1/orf-design/example.tsv"
    payload = b"id\tval\nl1-p0001\t42\n"
    digest = hashlib.sha256(payload).hexdigest()

    snap_root = tmp_path / prov.PINNED_SNAPSHOT_DIR
    snap_file = snap_root / commit / rel
    snap_file.parent.mkdir(parents=True)
    snap_file.write_bytes(b"tampered")
    (snap_root / prov.PINNED_MANIFEST_NAME).write_text(json.dumps({f"{commit}/{rel}": digest}))

    with pytest.raises(ValueError, match="checksum mismatch"):
        prov.read_pinned_bytes(commit, rel, repo=tmp_path)


def test_read_pinned_bytes_rejects_expected_sha256_mismatch(tmp_path: Path) -> None:
    commit = "8432c010"
    rel = "library/l1/orf-design/example.tsv"
    payload = b"hello\n"
    snap_file = tmp_path / prov.PINNED_SNAPSHOT_DIR / commit / rel
    snap_file.parent.mkdir(parents=True)
    snap_file.write_bytes(payload)

    with pytest.raises(ValueError, match="checksum mismatch"):
        prov.read_pinned_bytes(commit, rel, repo=tmp_path, expected_sha256="0" * 64)


def test_read_pinned_bytes_missing_without_git_raises_file_not_found(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError, match="not found"):
        prov.read_pinned_bytes("deadbeef", "missing/file.tsv", repo=tmp_path)


def test_sha256_file_and_verify_s3_objects(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    good = tmp_path / "good.txt"
    good.write_bytes(b"hello world\n")
    good_sha = hashlib.sha256(b"hello world\n").hexdigest()
    assert prov.sha256_file(good) == good_sha

    bad = tmp_path / "bad.txt"
    bad.write_bytes(b"different\n")

    class _FakeBucket:
        def __init__(self, files: dict[str, Path]) -> None:
            self._files = files

        def download(self, key: str) -> Path:
            if key not in self._files:
                raise FileNotFoundError(key)
            return self._files[key]

    import deliverome_analysis.s3 as s3_mod

    monkeypatch.setattr(s3_mod, "external", _FakeBucket({"ok": good, "tampered": bad}))
    monkeypatch.setattr(s3_mod, "processed", _FakeBucket({}))

    missing, changed = prov.verify_s3_objects(
        [
            ("external", "ok", good_sha),
            ("external", "tampered", good_sha),
            ("processed", "absent", good_sha),
        ],
        max_workers=2,
    )
    assert len(missing) == 1 and "processed://absent" in missing[0]
    assert changed == ["external://tampered"]
