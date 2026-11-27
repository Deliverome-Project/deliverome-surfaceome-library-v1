"""Focused offline unit tests for L1 core-module bug fixes."""

from __future__ import annotations

import email.message
import hashlib
import importlib.util
import io
import json
import threading
import urllib.error
from pathlib import Path
from types import ModuleType
from unittest.mock import MagicMock, patch

import numpy as np
import pytest

from deliverome_analysis.cloning import plasmid_qc as qc
from deliverome_analysis.cloning import seqtools
from deliverome_analysis.core import s3
from deliverome_analysis.surfaceome import membrane_side as ms

REPO_ROOT = Path(__file__).resolve().parents[1]
REVIEWER_PKG = REPO_ROOT / "library" / "l1" / "orf-design" / "reviewer-package"


def _load_reviewer_module(name: str) -> ModuleType:
    path = REVIEWER_PKG / f"{name}.py"
    spec = importlib.util.spec_from_file_location(f"_test_{name}", path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _write_tsv(path: Path, header: list[str], rows: list[list[str]]) -> None:
    lines = ["\t".join(header), *("\t".join(r) for r in rows)]
    path.write_text("\n".join(lines) + "\n")


# --------------------------------------------------------------------------- #
# 1. review_fingerprint: detect same-length sequence edits via content digests
# --------------------------------------------------------------------------- #


def test_review_fingerprint_detects_same_length_protein_edit(tmp_path: Path) -> None:
    rf = _load_reviewer_module("review_fingerprint")
    _write_tsv(
        tmp_path / "l1_tag_paths.tsv",
        [
            "candidate_id",
            "placement",
            "lane_override",
            "length_aa",
            "signal_peptide_aa",
            "sequence_sha256",
        ],
        [["l1-p0001", "internal", "", "400", "20", "aaa111"]],
    )
    _write_tsv(
        tmp_path / "constructs_l1_internal.tsv",
        [
            "construct_id",
            "design_id",
            "placement",
            "cassette_after_residue",
            "native_aa",
            "protein_aa",
            "native_sha256",
            "protein_sha256",
        ],
        [["l1-p0001-a1", "l1-p0001", "internal", "142", "400", "451", "aaa111", "bbb222"]],
    )
    brief = tmp_path / "brief_numbers.json"
    fp = rf.compute(tmp_path)
    brief.write_text(
        json.dumps({"tag_paths": {"placement": fp["placements"]}, "review_fingerprint": fp})
    )
    assert rf.check_brief_freshness(tmp_path, brief) == []

    # Mutate protein_sha256 while keeping integer lengths native_aa=400 and protein_aa=451 equal.
    _write_tsv(
        tmp_path / "constructs_l1_internal.tsv",
        [
            "construct_id",
            "design_id",
            "placement",
            "cassette_after_residue",
            "native_aa",
            "protein_aa",
            "native_sha256",
            "protein_sha256",
        ],
        [["l1-p0001-a1", "l1-p0001", "internal", "142", "400", "451", "aaa111", "ccc333"]],
    )
    problems = rf.check_brief_freshness(tmp_path, brief)
    assert len(problems) == 1
    assert "l1-p0001" in problems[0]


def test_review_fingerprint_detects_same_length_tag_path_sequence_edit(tmp_path: Path) -> None:
    rf = _load_reviewer_module("review_fingerprint")
    _write_tsv(
        tmp_path / "l1_tag_paths.tsv",
        [
            "candidate_id",
            "placement",
            "lane_override",
            "length_aa",
            "signal_peptide_aa",
            "sequence_sha256",
        ],
        [["l1-p0002", "C-terminal", "", "300", "0", "sha_v1"]],
    )
    brief = tmp_path / "brief_numbers.json"
    fp = rf.compute(tmp_path)
    brief.write_text(
        json.dumps({"tag_paths": {"placement": fp["placements"]}, "review_fingerprint": fp})
    )

    _write_tsv(
        tmp_path / "l1_tag_paths.tsv",
        [
            "candidate_id",
            "placement",
            "lane_override",
            "length_aa",
            "signal_peptide_aa",
            "sequence_sha256",
        ],
        [["l1-p0002", "C-terminal", "", "300", "0", "sha_v2"]],
    )
    problems = rf.check_brief_freshness(tmp_path, brief)
    assert len(problems) == 1
    assert "l1-p0002" in problems[0]


# --------------------------------------------------------------------------- #
# 2. plasmid_qc.check_shuttle_sites: circular origin-wrapping deduplication
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("rot", [0, 1, 2, 3, 4, 5, 6, -1, -2, -3, -4, -5])
def test_check_shuttle_sites_counts_origin_and_near_origin_sites_once(rot: int) -> None:
    # Two Esp3I sites: forward CGTCTC at index 0, reverse GAGACG in the middle.
    base = "CGTCTC" + "A" * 40 + "GAGACG" + "T" * 40
    rotated = base[rot:] + base[:rot]
    finding = qc.check_shuttle_sites(rotated)
    assert finding.ok, f"rot={rot}: {finding.detail}"
    assert "2 Esp3I" in finding.detail


def test_check_shuttle_sites_handles_short_and_wrong_count_sequences() -> None:
    # A single site starting at index 0 must report 1, not 2 from wrap duplication.
    one_at_zero = "CGTCTC" + "A" * 50
    f_one = qc.check_shuttle_sites(one_at_zero)
    assert not f_one.ok
    assert f_one.detail.startswith("1 Esp3I")

    # Empty or sub-site-length inputs return 0 cleanly.
    f_empty = qc.check_shuttle_sites("")
    assert not f_empty.ok
    assert f_empty.detail.startswith("0 Esp3I")
    f_short = qc.check_shuttle_sites("CGTCT")
    assert not f_short.ok
    assert f_short.detail.startswith("0 Esp3I")


# --------------------------------------------------------------------------- #
# 3. seqtools.LibraryLinkerPlan.audit: 3'-flank spanning_tract coverage
# --------------------------------------------------------------------------- #


def _make_plan(
    variants: dict[str, str],
    flanks: dict[str, tuple[str, str]],
    *,
    seed: int = 8,
) -> seqtools.LibraryLinkerPlan:
    lv = {
        name: seqtools.LinkerVariant(
            sequence=dna,
            peptide=seqtools.translate(dna),
            gc_content=seqtools.gc_content(dna),
            max_internal_repeat=0,
            enzymes=(),
            recombinogenic_k=15,
            attempts=1,
        )
        for name, dna in variants.items()
    }
    return seqtools.LibraryLinkerPlan(
        variants=lv,
        homology=seqtools.junction_homology(flanks, seed=seed),
        flanks=flanks,
        order=list(variants),
        recombinogenic_k=15,
        max_junction_identity=9,
        seed=seed,
        peptide=seqtools.translate(next(iter(variants.values()))),
    )


def test_library_linker_plan_audit_counts_3prime_flank_spanning_tract() -> None:
    shared3 = "TAACCAGTCGTTAGCTGCCGTTCC"  # 24 nt shared on 3' side only
    flanks = {
        "a": ("AAAAAAAAAAAAAAAAAAAAAAAA", shared3),
        "b": ("CCCCCCCCCCCCCCCCCCCCCCCC", shared3),
    }
    # Cassettes share a 6-nt suffix ("GGCAGC") but differ at their 5' end.
    cas_a = "GGTGGTGGTGGCAGC"
    cas_b = "GGAGGAGGCGGCAGC"
    plan = _make_plan({"a": cas_a, "b": cas_b}, flanks)
    report = plan.audit()
    assert report["junction_runthrough_5"] == 0
    assert report["junction_runthrough_3"] == 6
    assert report["spanning_tract"] == len(shared3) + 6

    # Identical cassettes with 3'-flank-only homology span flank3 + full cassette.
    plan_same = _make_plan({"a": cas_a, "b": cas_a}, flanks)
    report_same = plan_same.audit()
    assert report_same["spanning_tract"] == len(shared3) + len(cas_a)


def test_library_linker_plan_audit_picks_max_of_5prime_and_3prime_tracts() -> None:
    shared5 = "CGATCGATCGAT"  # 12 nt
    shared3 = "TTGACCAGTCAGTTGACCAGTCAGTTGACC"  # 30 nt
    flanks = {
        "a": ("AAAA" + shared5, shared3 + "AAAA"),
        "b": ("TTTT" + shared5, shared3 + "TTTT"),
    }
    # 5' common prefix = 4 ("GGAG"), 3' common suffix = 12 ("GGCGGCGGCAGC") -> 3' span wins.
    cas_a = "GGAGGGGGCGGCGGCAGC"
    cas_b = "GGAGTCGGCGGCGGCAGC"
    plan = _make_plan({"a": cas_a, "b": cas_b}, flanks)
    report = plan.audit()
    assert report["junction_runthrough_5"] == 4
    assert report["junction_runthrough_3"] == 12
    assert report["spanning_tract"] == len(shared3) + 12


# --------------------------------------------------------------------------- #
# 4. membrane_side: NaN/non-finite pLDDT fail-closed & Kabsch superpose
# --------------------------------------------------------------------------- #


def test_site_verdict_rejects_non_finite_plddt_and_z() -> None:
    th = ms.GateThresholds(
        outside=15.0,
        margin=3.0,
        core=10.0,
        min_plddt=50.0,
        min_topology_ratio=2.0,
        max_kr_excess_outside=3,
    )
    for bad in (float("nan"), float("inf"), float("-inf")):
        assert ms.site_verdict(-23.0, bad, True, th) == "unscorable"
        assert ms.site_verdict(20.0, bad, True, th) == "unscorable"
        assert ms.site_verdict(bad, 90.0, True, th) == "unscorable"
    # Finite low pLDDT keeps documented hold behavior on failing heights.
    assert ms.site_verdict(-23.0, 30.0, True, th) == "hold_low_plddt"


def test_superpose_handles_planar_and_degenerate_inputs() -> None:
    # Planar points under a 3D rotation
    planar = np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 3.0, 0.0], [2.0, 3.0, 0.0]])
    rot_true = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    shifted = planar @ rot_true.T + np.array([1.0, -2.0, 4.0])
    r, t, rmsd = ms.superpose(planar, shifted)
    assert np.isclose(np.linalg.det(r), 1.0)
    assert np.allclose(r @ r.T, np.eye(3), atol=1e-12)
    assert rmsd == pytest.approx(0.0, abs=1e-12)
    assert np.allclose(planar @ r.T + t, shifted, atol=1e-12)

    # Planar reflection is realizable as a proper 3D rotation (det == +1)
    reflected = planar * np.array([-1.0, 1.0, 1.0])
    r_ref, _, rmsd_ref = ms.superpose(planar, reflected)
    assert np.isclose(np.linalg.det(r_ref), 1.0)
    assert rmsd_ref == pytest.approx(0.0, abs=1e-12)

    # Coincident points still yield a proper rotation matrix
    coincident = np.zeros((4, 3))
    r_zero, t_zero, rmsd_zero = ms.superpose(coincident, coincident + 5.0)
    assert np.isclose(np.linalg.det(r_zero), 1.0)
    assert np.allclose(t_zero, [5.0, 5.0, 5.0])
    assert rmsd_zero == pytest.approx(0.0, abs=1e-12)

    # Empty or wrong-shape inputs raise ValueError rather than returning NaN
    with pytest.raises(ValueError, match="non-empty"):
        ms.superpose(np.empty((0, 3)), np.empty((0, 3)))
    with pytest.raises(ValueError, match="non-empty"):
        ms.superpose(np.zeros((3, 2)), np.zeros((3, 2)))


# --------------------------------------------------------------------------- #
# 5. s3.exists (404 vs 403/network) & atomic s3.download unique temp files
# --------------------------------------------------------------------------- #


class _FakeClientError(Exception):
    def __init__(self, code: str, status: int) -> None:
        super().__init__(f"{code} ({status})")
        self.response = {
            "Error": {"Code": code},
            "ResponseMetadata": {"HTTPStatusCode": status},
        }


def test_s3_exists_returns_false_only_on_404_and_reraises_403_or_network(monkeypatch) -> None:
    monkeypatch.delenv("DELIVEROME_S3_MIRROR", raising=False)

    fake = MagicMock()
    with patch.object(s3.external, "_client_factory", lambda: fake):
        fake.head_object.return_value = {"ETag": '"abc"'}
        assert s3.external.exists("a/b.tsv") is True

        for err in (
            _FakeClientError("404", 404),
            _FakeClientError("NoSuchKey", 404),
            _FakeClientError("NotFound", 404),
            _FakeClientError("Unknown", 404),
        ):
            fake.head_object.side_effect = err
            assert s3.external.exists("a/missing.tsv") is False

        fake.head_object.side_effect = _FakeClientError("403", 403)
        with pytest.raises(_FakeClientError):
            s3.external.exists("a/forbidden.tsv")

        fake.head_object.side_effect = ConnectionError("network unreachable")
        with pytest.raises(ConnectionError, match="network unreachable"):
            s3.external.exists("a/net.tsv")


def test_s3_download_uses_unique_temp_files_and_cleans_up_on_error(
    tmp_path: Path, monkeypatch
) -> None:
    monkeypatch.delenv("DELIVEROME_S3_MIRROR", raising=False)
    monkeypatch.setattr(s3, "_CACHE_ROOT", tmp_path / "cache")

    barrier = threading.Barrier(2)
    seen_targets: list[str] = []
    lock = threading.Lock()

    def concurrent_download(_bucket: str, _key: str, dest: str) -> None:
        with lock:
            seen_targets.append(dest)
        barrier.wait(timeout=5)
        Path(dest).write_bytes(b"clean-payload")

    fake = MagicMock()
    fake.head_object.return_value = {"ETag": '"etag1"'}
    fake.download_file.side_effect = concurrent_download

    results: list[Path] = []
    errors: list[BaseException] = []

    def worker() -> None:
        try:
            results.append(s3.external.download("ref/v1/data.bin"))
        except BaseException as exc:
            errors.append(exc)

    with patch.object(s3.external, "_client_factory", lambda: fake):
        t1, t2 = threading.Thread(target=worker), threading.Thread(target=worker)
        t1.start()
        t2.start()
        t1.join(timeout=5)
        t2.join(timeout=5)

    assert not errors
    assert len(seen_targets) == 2
    assert seen_targets[0] != seen_targets[1]
    assert len(results) == 2 and results[0] == results[1]
    assert results[0].read_bytes() == b"clean-payload"
    assert list(results[0].parent.glob("*.part")) == []

    # Failed download removes its .part file and never leaves a cached entry.
    def failing_download(_bucket: str, _key: str, dest: str) -> None:
        Path(dest).write_bytes(b"truncated")
        raise RuntimeError("socket dropped")

    fake.head_object.return_value = {"ETag": '"etag2"'}
    fake.download_file.side_effect = failing_download
    with (
        patch.object(s3.external, "_client_factory", lambda: fake),
        pytest.raises(RuntimeError, match="socket dropped"),
    ):
        s3.external.download("ref/v1/broken.bin")
    cache_dir = tmp_path / "cache" / s3.external.name / "ref" / "v1"
    assert not any(p.name.startswith("broken.bin") for p in cache_dir.iterdir())
    assert list(cache_dir.glob("*.part")) == []


# --------------------------------------------------------------------------- #
# 6. public_data: HTTP 429/5xx retry with Retry-After and atomic writes
# --------------------------------------------------------------------------- #


class _FakeHTTPResponse(io.BytesIO):
    def __enter__(self) -> _FakeHTTPResponse:
        return self

    def __exit__(self, *_args: object) -> None:
        self.close()


class _FailingStreamResponse:
    def __init__(self, first_chunk: bytes) -> None:
        self._sent = False
        self._first = first_chunk

    def read(self, _n: int = -1) -> bytes:
        if not self._sent:
            self._sent = True
            return self._first
        raise ConnectionResetError("connection reset mid-stream")

    def __enter__(self) -> _FailingStreamResponse:
        return self

    def __exit__(self, *_args: object) -> None:
        return None


def _http_error(code: int, retry_after: str | None = None) -> urllib.error.HTTPError:
    hdrs = email.message.Message()
    if retry_after is not None:
        hdrs["Retry-After"] = retry_after
    return urllib.error.HTTPError("https://example.test/obj", code, "err", hdrs, io.BytesIO(b""))


def test_public_data_download_retries_429_and_5xx_with_retry_after(
    tmp_path: Path, monkeypatch
) -> None:
    pd_mod = _load_reviewer_module("public_data")
    payload = b"hello public bundle"
    calls = 0
    sleeps: list[float] = []

    def fake_urlopen(_req: object, timeout: int = 900) -> _FakeHTTPResponse:
        nonlocal calls
        calls += 1
        if calls == 1:
            raise _http_error(429, retry_after="3")
        if calls == 2:
            raise _http_error(503)
        return _FakeHTTPResponse(payload)

    monkeypatch.setattr(pd_mod.urllib.request, "urlopen", fake_urlopen)
    monkeypatch.setattr(pd_mod.time, "sleep", sleeps.append)

    dest = tmp_path / "out" / "file.bin"
    pd_mod._download_file("https://example.test/obj", dest)

    assert calls == 3
    assert sleeps == [3.0, 2.0]
    assert dest.read_bytes() == payload
    assert list(dest.parent.glob("*.part")) == []


def test_public_data_download_is_atomic_and_does_not_retry_404(tmp_path: Path, monkeypatch) -> None:
    pd_mod = _load_reviewer_module("public_data")
    sleeps: list[float] = []
    monkeypatch.setattr(pd_mod.time, "sleep", sleeps.append)

    # 404 fails immediately with zero retries.
    def fake_404(_req: object, timeout: int = 900) -> _FakeHTTPResponse:
        raise _http_error(404)

    monkeypatch.setattr(pd_mod.urllib.request, "urlopen", fake_404)
    dest = tmp_path / "missing.bin"
    with pytest.raises(urllib.error.HTTPError) as exc_info:
        pd_mod._download_file("https://example.test/missing", dest)
    assert exc_info.value.code == 404
    assert sleeps == []
    assert not dest.exists()

    # Mid-stream failure leaves no partial file at dest and no leftover .part temp files.
    monkeypatch.setattr(
        pd_mod.urllib.request,
        "urlopen",
        lambda _req, timeout=900: _FailingStreamResponse(b"partial-bytes"),
    )
    partial_dest = tmp_path / "partial.bin"
    with pytest.raises(ConnectionResetError, match="mid-stream"):
        pd_mod._download_file("https://example.test/partial", partial_dest)
    assert not partial_dest.exists()
    assert list(tmp_path.glob("*.part")) == []

    # fetch() end-to-end with verified SHA-256 succeeds cleanly.
    good = b"verified-content"
    digest = hashlib.sha256(good).hexdigest()
    monkeypatch.setattr(
        pd_mod,
        "locations",
        lambda *_homes: [
            {
                "bucket": "external",
                "key": "demo/item.txt",
                "bytes": str(len(good)),
                "sha256": digest,
                "home": "public_r2",
                "license": "CC0",
                "location": "https://example.test/item.txt",
            }
        ],
    )
    monkeypatch.setattr(
        pd_mod.urllib.request,
        "urlopen",
        lambda _req, timeout=900: _FakeHTTPResponse(good),
    )
    out = tmp_path / "bundle"
    assert pd_mod.fetch(out, workers=2) == 0
    assert (out / "s3" / "external" / "demo" / "item.txt").read_bytes() == good


def test_domesticate_avoids_creating_polya_and_splice_donor_motifs() -> None:
    # Exact contexts where naive domestication created ATTAAA (KLRB1/SLC46A1) and GTGAGT (EPHA1)
    cases = [
        ("ATGAGACCGGGTCTCTTAAACTAA", ("BsaI",)),  # GGTCTC + TTAAAC (ATTAAA trap)
        ("ATGTTTGGTGAGACCTTAAAGTAA", ("BsaI",)),  # GAGACC + TTAAAG (ATTAAA trap)
        ("ATGACCGTCTCTGAGTGGTAA", ("BsaI", "Esp3I")),  # GTCTCT + GAGTGG (GTGAGT trap)
    ]
    for cds, enzymes in cases:
        dom, warnings = seqtools.domesticate(cds, list(enzymes), preserve_translation=True, seed=0)
        assert not warnings
        assert seqtools.translate(dom) == seqtools.translate(cds)
        assert not seqtools.find_sites(dom, list(enzymes))
        assert "AATAAA" not in dom and "ATTAAA" not in dom
        assert "GTAAGT" not in dom and "GTGAGT" not in dom
