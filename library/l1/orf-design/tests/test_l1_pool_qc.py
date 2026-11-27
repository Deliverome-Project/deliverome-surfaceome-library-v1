"""Plasmid QC over the whole L1 pool (``l1_pool_qc.py``).

The unit tests prove each pool-level measure on synthetic inserts, including that it fails on
the defect it exists for. The integration test runs the whole QC over the committed arm tables
and requires zero error-level failures; it is an expected failure until every upstream arm has
its own linker encoding (``l1_arm_cassette_encodings.tsv``), because until then 2,942 inserts
carry p003's one shared cassette DNA (Becca, 2026-11-21: "we still don't want shared linkers
in the middle").
"""

from __future__ import annotations

import random
import sys
from collections import Counter
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))

import l1_pool_qc as pool  # noqa: E402

from deliverome_analysis.cloning import plasmid_qc as qc  # noqa: E402
from deliverome_analysis.cloning import seqtools  # noqa: E402


def _rand(n: int, rng: random.Random) -> str:
    return "".join(rng.choice("ACGT") for _ in range(n))


def _brute(windows: list[str], owners: list[int], i: int) -> int:
    best = 0
    for j, w in enumerate(windows):
        if owners[j] == owners[i]:
            continue
        best = max(best, seqtools.longest_shared_tract(windows[i], w)[0])
    return best


def test_the_plasmid_qc_constant_is_build_constructs_cassette():
    from build_constructs import CASSETTE, CASSETTE_AA

    assert qc.P003_CASSETTE_NT == CASSETTE
    assert qc.P003_CASSETTE_PEPTIDE == CASSETTE_AA


def test_the_synthesis_gc_ceiling_is_the_order_s():
    import synthesis_compliance

    assert qc.SYNTHESIS_GC_MAX == synthesis_compliance.GC50_MAX


def test_longest_shared_matches_brute_force_and_ignores_the_same_owner():
    rng = random.Random(7)
    core = _rand(70, rng)
    windows = [
        _rand(90, rng) + core + _rand(90, rng),
        _rand(60, rng) + core[:45] + _rand(100, rng),
        _rand(200, rng),
        _rand(40, rng) + core + _rand(40, rng),  # same owner as window 0
        _rand(120, rng) + core[10:62] + _rand(30, rng),
    ]
    owners = [0, 1, 2, 0, 3]
    tract, partner = pool.longest_shared(windows, owners)
    assert tract == [_brute(windows, owners, i) for i in range(len(windows))]
    assert tract[0] == 52 and partner[0] == 4  # not window 3, which is the same member's
    for i, j in enumerate(partner):
        if j is not None:
            assert owners[j] != owners[i]


def test_a_whole_shared_cassette_shows_as_a_tract_at_least_its_length():
    rng = random.Random(11)
    windows = [_rand(150, rng) + qc.P003_CASSETTE_NT + _rand(150, rng) for _ in range(3)]
    tract, _ = pool.longest_shared(windows, [0, 1, 2])
    assert min(tract) >= len(qc.P003_CASSETTE_NT)


def test_cassette_sharing_fails_shared_dna_and_names_p003s():
    holders = Counter([qc.P003_CASSETTE_NT, qc.P003_CASSETTE_NT, "OWN"])
    shared = pool.cassette_sharing(qc.P003_CASSETTE_NT, holders)
    assert not shared.ok and "p003's shared encoding" in shared.detail
    assert pool.cassette_sharing("OWN", holders).ok
    assert not pool.cassette_sharing("", holders).ok


def test_cassette_window_finds_the_cassette_and_stops_at_the_barcode():
    rng = random.Random(3)
    head = "ATG" + "GCC" * 20
    insert = (
        pool.OPENING[:-3] + head + qc.P003_CASSETTE_NT + "TAA" + pool.READOUT_BLOCK + _rand(10, rng)
    )
    cas, window = pool.cassette_window(insert)
    assert cas == qc.P003_CASSETTE_NT
    # the cassette sits 92 nt in, so the window runs from the start to the barcode slot
    assert window == insert[: insert.index("N")]


def test_assembly_puts_the_insert_in_the_slot_and_keeps_the_backbone():
    bb, feats = pool.backbone()
    assert {"EFS_v1", "PuroR", "WPRE"} <= set(feats)
    ins, _ = pool.load_inserts()
    one = ins[0]
    vec = pool.assemble(bb, one["insert"])
    assert (
        vec[qc.BACKBONE_SLOT_START : qc.BACKBONE_SLOT_START + len(one["insert"])] == one["insert"]
    )
    spec = qc.Spec(plasmid="x", designed_insert=one["insert"], backbone=bb, backbone_features=feats)
    assert qc.check_backbone(vec, spec).ok
    assert qc.check_backbone_features(vec, spec).ok


def test_an_arm_without_an_issued_encoding_fails_both_linker_checks():
    """The verdict the pool QC gives the old order: p003's shared DNA in an upstream arm."""
    bb, feats = pool.backbone()
    ins, _ = pool.load_inserts()
    up = next(i for i in ins if i["rung"] == "nterm")
    cas, _ = pool.cassette_window(up["insert"])
    shared = up["insert"].replace(cas, qc.P003_CASSETTE_NT)
    up = {**up, "insert": shared}
    spec = pool.spec_for(up, bb, feats, {}, {}, {})
    assert spec.cassette_nt == qc.P003_CASSETTE_NT
    found = {f.check: f for f in qc.check_plasmid(pool.assemble(bb, shared), spec)}
    assert not found["diversified_linker"].ok
    assert not found["linkers_not_shared"].ok


@pytest.mark.xfail(
    not pool.ARM_ENCODINGS.exists(),
    reason="linkers not yet recoded in N-terminal/internal arms",
    strict=True,
)
def test_the_whole_pool_passes_plasmid_qc():
    rows, summary = pool.run()
    assert summary["failing_checks"] == [], {
        c: summary["checks"][c] for c in summary["failing_checks"] if c in summary["checks"]
    } or summary["pool_checks"]
    assert summary["inserts"] == len(rows) == summary["cassette_dna_variants"]
    assert summary["largest_shared_cassette_group"] == 1
    # every per-insert check ran on every insert (none silently skipped), except the
    # membrane anchor, which needs topology the arm tables do not carry
    for check, counts in summary["checks"].items():
        if check != "membrane_anchor":
            assert counts.get(pool.SKIP, 0) == 0, (check, counts)
