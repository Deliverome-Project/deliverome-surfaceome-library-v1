"""Invariants of domestication and barcodes that no other test checked.

The audit (``audit_domestication_barcodes.py``) recomputes them from the final inserts and the
recoded order of record. These tests read its committed outputs and the committed barcode
block, so they run offline in CI:

* every domestication edit is synonymous, and the inserts carry no unlogged native change;
* every final insert has exactly the two wrapper Esp3I sites, the designed Kozak, one stop, the
  protein its arm table records, and one barcode slot;
* the reserved barcode block meets the s1e2 pool's own claims (GC 40-60%, homopolymer at most 2,
  no ATG, full-length Hamming distance at least 4), and no barcode occurs in the backbone;
* the Twist recodes change no protein, cassette, readout block or barcode slot, add no gated
  Type IIS site, and lengthen no swap tract past the bar.
"""

from __future__ import annotations

import csv
import json
import re
from pathlib import Path

import numpy as np

from deliverome_analysis import seqtools

S4E9 = Path(__file__).resolve().parents[1]
GATE = 100


def _tsv(name: str) -> list[dict]:
    with (S4E9 / name).open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _audit() -> dict:
    return json.loads((S4E9 / "l1_domestication_barcode_audit.json").read_text())


def test_every_logged_domestication_edit_is_synonymous():
    edits = _tsv("l1_domestication_edits.tsv")
    assert edits, "no domestication edits logged"
    bad = [
        (e["design_id"], e["codon"])
        for e in edits
        if seqtools.translate(e["old_codon"]) != seqtools.translate(e["new_codon"])
        or e["synonymous"] != "True"
    ]
    assert not bad, f"nonsynonymous domestication edits: {bad[:5]}"
    assert _audit()["domestication"].get("translation_changed", 0) == 0
    assert _audit()["domestication"].get("edit_count_mismatch", 0) == 0


def test_every_final_insert_passes_the_gated_checks():
    rows = _tsv("l1_insert_audit.tsv")
    assert len(rows) == _audit()["inserts"]
    for col in ("wrapper_only", "kozak_and_start", "single_stop_at_end", "protein_matches_table"):
        bad = [r["id"] for r in rows if r[col] != "True"]
        assert not bad, f"{col} fails for {bad[:5]}"
    assert all(r["barcode_slots"] == "1" for r in rows)


def test_inserts_carry_no_unlogged_native_change():
    # The only native change an arm makes outside its cassette or leader is a logged junction
    # edit (constructs_l1_cterm.tsv junction_edit, e.g. SLC41A2 GAC->GAT).
    # Joined on the permanent design id: member ids follow plate order and move on a renumber
    # (l1-reviewer-package, 2026-11-13, after the 25 rescue members were retired).
    cterm = {r["design_id"]: r for r in _tsv("constructs_l1_cterm.tsv")}
    for r in _tsv("l1_insert_audit.tsv"):
        off = r["native_nt_not_verbatim"]
        if off in ("", "0"):
            continue
        assert off != "-1", f"{r['id']}: leader arm is not leader + domesticated ORF"
        d = r["design_id"]
        assert d in cterm and cterm[d]["junction_edit"], f"{r['id']} ({d}): unlogged change"


def test_barcode_block_meets_the_pool_claims():
    codes = [r["barcode"] for r in _tsv("l1_barcode_block_provisional.tsv")]
    assert all(0.40 <= seqtools.gc_content(c) <= 0.60 for c in codes)
    assert all(max(len(m.group(0)) for m in re.finditer(r"(.)\1*", c)) <= 2 for c in codes)
    assert not any("ATG" in c for c in codes)
    lut = {b: i for i, b in enumerate("ACGT")}
    arr = np.array([[lut[b] for b in c] for c in codes], dtype=np.int8)
    for i in range(len(arr) - 1):
        assert (arr[i + 1 :] != arr[i]).sum(axis=1).min() >= 4, f"{codes[i]}: Hamming < 4"


def test_no_barcode_occurs_in_the_backbone():
    text = (S4E9 / "p0039_EFS_IRES_puro_backbone.fasta").read_text().splitlines()
    bb = "".join(line.strip() for line in text if not line.startswith(">")).upper()
    hay = bb + bb + "|" + seqtools.reverse_complement(bb + bb)
    hits = [r["barcode"] for r in _tsv("l1_barcode_block_provisional.tsv") if r["barcode"] in hay]
    assert not hits, f"barcodes found in the backbone: {hits[:5]}"
    assert _audit()["barcodes"]["found_in_an_insert_or_backbone"] == []


def test_recoded_order_keeps_every_invariant():
    rec = _audit()["recoded"]
    allowed = {"new unscreened Type IIS site (AarI/PaqCI/BtgZI/SapI)"}
    assert set(rec["problems"]) <= allowed, rec["problems"]
    assert rec["cterm_terminal_linker_blocks_shared"] == 0
    # The pinned order is the current build: every insert matched, none rebuilt since.
    assert rec["recoded_inserts"] == _audit()["inserts"]
    assert rec["not_in_original_set"] == 0
    assert rec.get("rebuilt_since_order_not_compared", []) == []
    rows = _tsv("l1_insert_audit.tsv")
    grew = [
        r["id"]
        for r in rows
        if int(r["recoded_worst_partner_tract_nt"]) > GATE >= int(r["worst_partner_tract_nt"])
    ]
    assert not grew, f"recodes pushed a swap tract over {GATE} nt: {grew[:5]}"


def test_retired_rescue_members_are_not_audited_as_ordered():
    """The 25 expression-rescue members (l1-p2361..l1-p2385) are retired (Becca 2026-11-13,
    one arm per gene; PR #312). None may come back as an audited insert, and the stale order of
    record must say how many it still holds rather than pass silently."""
    retired = {f"l1-p{n:04d}" for n in range(2361, 2386)}
    rows = _tsv("l1_insert_audit.tsv")
    assert not {r["design_id"] for r in rows} & retired
    assert "expression rescue" not in _audit()["by_rung"]
    rec = _audit()["recoded"]
    assert rec["retired_members_still_in_order"] in (0, len(retired))
