"""Offline unit tests for Zenodo archive staging, verification, and mocked upload."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path
from typing import Any

import pytest

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE / "reviewer-package"))

import mint_zenodo_archive as mza  # noqa: E402


def _make_dummy_repo(path: Path) -> None:
    path.mkdir(parents=True)
    subprocess.run(["git", "init", "-q", str(path)], check=True)
    subprocess.run(["git", "-C", str(path), "config", "user.email", "test@example.com"], check=True)
    subprocess.run(["git", "-C", str(path), "config", "user.name", "Test"], check=True)
    (path / "README.md").write_text("# Dummy review repo\n")
    subprocess.run(["git", "-C", str(path), "add", "README.md"], check=True)
    subprocess.run(["git", "-C", str(path), "commit", "-q", "-m", "init"], check=True)


def test_stage_and_verify_archive_and_detect_corruption(tmp_path: Path) -> None:
    repo = tmp_path / "review-repo"
    _make_dummy_repo(repo)
    bundle = tmp_path / "l1-reviewer-bundle.tar.gz"
    bundle.write_bytes(b"fake-tar-gz-bytes")
    image = tmp_path / "l1-image.tar.gz"
    image.write_bytes(b"fake-oci-image")

    out = tmp_path / "staged"
    sums = mza.stage_archive(
        review_repo=repo,
        bundle_tar=bundle,
        image_tar=image,
        out_dir=out,
        version="1.0.0",
    )
    assert set(sums) == {
        "deliverome-surfaceome-library-v1.bundle",
        "l1-reviewer-bundle.tar.gz",
        "deliverome-surfaceome-library-v1-image.tar.gz",
        "RERUN.md",
        "zenodo_metadata.json",
    }
    assert mza.verify_archive(out) == sums

    # Flip one byte in the staged bundle and assert verify_archive fails loudly.
    (out / "l1-reviewer-bundle.tar.gz").write_bytes(b"corrupted-bytes!!")
    with pytest.raises(ValueError, match=r"checksum mismatch: l1-reviewer-bundle\.tar\.gz"):
        mza.verify_archive(out)


class _FakeResp:
    def __init__(self, payload: dict[str, Any]) -> None:
        self._payload = payload

    def raise_for_status(self) -> None:
        return None

    def json(self) -> dict[str, Any]:
        return self._payload


class _FakeSession:
    def __init__(self) -> None:
        self.calls: list[tuple[str, str]] = []
        self.uploaded: list[str] = []

    def post(self, url: str, **_kw: Any) -> _FakeResp:
        self.calls.append(("POST", url))
        return _FakeResp(
            {
                "id": 999,
                "links": {
                    "bucket": "https://sandbox.zenodo.org/api/files/bucket-1",
                    "html": "https://sandbox.zenodo.org/deposit/999",
                },
            }
        )

    def put(self, url: str, **_kw: Any) -> _FakeResp:
        self.calls.append(("PUT", url))
        if "/files/bucket-1/" in url:
            self.uploaded.append(url.rsplit("/", 1)[-1])
        return _FakeResp({"id": 999, "links": {"html": "https://sandbox.zenodo.org/deposit/999"}})


def test_upload_archive_stays_draft_by_default_and_uploads_sums(tmp_path: Path) -> None:
    repo = tmp_path / "review-repo"
    _make_dummy_repo(repo)
    bundle = tmp_path / "l1-reviewer-bundle.tar.gz"
    bundle.write_bytes(b"bundle")
    out = tmp_path / "staged"
    mza.stage_archive(review_repo=repo, bundle_tar=bundle, out_dir=out, version="1.0.0")

    fake = _FakeSession()
    dep = mza.upload_archive(out, token="tok", sandbox=True, publish=False, session=fake)
    assert dep["id"] == 999
    assert "SHA256SUMS" in fake.uploaded
    assert "deliverome-surfaceome-library-v1.bundle" in fake.uploaded
    assert not any(url.endswith("/actions/publish") for _, url in fake.calls)
