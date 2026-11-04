"""The offline mirror: replay reads from a local folder, never S3, and refuse writes.

Runs without credentials; the record-mode test mocks the S3 client.
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from deliverome_analysis import s3


@pytest.fixture
def mirror(tmp_path: Path, monkeypatch) -> Path:
    root = tmp_path / "mirror"
    (root / "external/uniprot/2026_03").mkdir(parents=True)
    (root / "external/uniprot/2026_03/human.tsv").write_text("a\tb\n1\t2\n")
    (root / "processed/sprint4/s4e9-x").mkdir(parents=True)
    (root / "processed/sprint4/s4e9-x/t.tsv").write_text("x\n")
    monkeypatch.setenv("DELIVEROME_S3_MIRROR", str(root))
    monkeypatch.delenv("DELIVEROME_S3_MIRROR_MODE", raising=False)
    return root


def _no_client():
    raise AssertionError("the mirror must not build an S3 client")


def test_reads_come_from_the_mirror(mirror: Path, tmp_path: Path) -> None:
    with patch.object(s3.external, "_client_factory", _no_client):
        assert s3.external.download("uniprot/2026_03/human.tsv").read_text() == "a\tb\n1\t2\n"
        dest = s3.external.download("uniprot/2026_03/human.tsv", tmp_path / "copy.tsv")
        assert dest.read_text().startswith("a\tb")
        assert s3.external.read_tsv("uniprot/2026_03/human.tsv")["b"].tolist() == [2]
        assert s3.external.exists("uniprot/2026_03/human.tsv")
        assert not s3.external.exists("uniprot/2026_03/mouse.tsv")
        assert s3.external.list_keys("uniprot/") == ["uniprot/2026_03/human.tsv"]


def test_buckets_are_kept_apart(mirror: Path) -> None:
    with patch.object(s3.processed, "_client_factory", _no_client):
        assert not s3.processed.exists("uniprot/2026_03/human.tsv")
        assert s3.processed.list_keys() == ["sprint4/s4e9-x/t.tsv"]


def test_missing_key_names_the_mirror(mirror: Path) -> None:
    with pytest.raises(FileNotFoundError, match="offline mirror"):
        s3.external.download("uniprot/2026_03/absent.tsv")


def test_uploads_are_refused(mirror: Path, tmp_path: Path) -> None:
    f = tmp_path / "f.tsv"
    f.write_text("x\n")
    with pytest.raises(RuntimeError, match="offline replay"):
        s3.processed.upload(f, "sprint4/s4e9-x/new.tsv")
    with pytest.raises(RuntimeError, match="offline replay"):
        s3.external.upload(f, "uniprot/2026_03/new.tsv")


def test_record_mode_copies_what_is_read(tmp_path: Path, monkeypatch) -> None:
    root = tmp_path / "mirror"
    monkeypatch.setenv("DELIVEROME_S3_MIRROR", str(root))
    monkeypatch.setenv("DELIVEROME_S3_MIRROR_MODE", "record")
    monkeypatch.setattr(s3, "_CACHE_ROOT", tmp_path / "cache")
    fake = MagicMock()
    fake.head_object.return_value = {"ETag": '"abc"'}
    fake.download_file.side_effect = lambda b, k, path: Path(path).write_text("payload\n")
    with patch.object(s3.external, "_client_factory", lambda: fake):
        s3.external.download("ensembl/release-115/x.fa")
    assert (root / "external/ensembl/release-115/x.fa").read_text() == "payload\n"


def test_unknown_mode_is_an_error(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setenv("DELIVEROME_S3_MIRROR", str(tmp_path))
    monkeypatch.setenv("DELIVEROME_S3_MIRROR_MODE", "write")
    with pytest.raises(ValueError, match="read' or 'record"):
        s3.external.exists("a/b/c")
