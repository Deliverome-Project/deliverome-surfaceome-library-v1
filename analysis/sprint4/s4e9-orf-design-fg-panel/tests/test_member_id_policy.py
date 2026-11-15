"""Design ids are the join key; member ids survive only in dated records.

From the first per-insert cut (2026-11-15) the plate position belongs to an insert
(``l1-iNNNN``, ``l1_insert_numbering.tsv``) and a member is named by its design id everywhere
current, so any ``l1-m`` id outside the dated records is stale. The rule below is the one that
held while member ids were plate positions; the checks still apply, with the current "member
id" of every member now equal to its design id.

Decision (coordinator under Becca's 2026-11-12 delegation, 2026-11-13): a member id
(``l1-mNNNN``) is the plate position and moves at every membership change until the order
freeze. It is authoritative only in ``l1_member_renumbering.tsv``, the order plan and the
records the renumber stage re-keys. Dated audit and screen tables keep the member id that was
current when they were written. They are not re-keyed, so they must never be joined on it.
The stable key is the design id (``l1-pNNNN``).

These tests check the data side of that rule, cheaply:

* every s4e9 table that carries a member id also carries a design id, unless it is listed in
  :data:`NO_DESIGN_ID` with the stable key to use instead;
* every design id in those tables names a real design;
* the tables the renumber stage keeps current agree with ``l1_member_renumbering.tsv``.

Implemented by the Claude Code session l1-reviewer-package.
"""

from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))

import renumber_library_members as rn  # noqa: E402

MEMBER = re.compile(r"^l1-m\d{4}$")
DESIGN = re.compile(r"^l[12]-p\d{4}$")

#: Tables that carry a member id but no design id, with the stable key a join must use.
NO_DESIGN_ID = {
    "l1_rbp_motif_cassettes.tsv": "none recorded; regenerate keyed on design id before use",
    "membrane-side-gate/calibration_tmdet.tsv": "gene_symbol (calibration set, not a join input)",
}
#: Records the renumber stage re-keys on every membership change, so their member ids are
#: always current, but which carry no design id column. Join them to design ids through
#: ``l1_member_renumbering.tsv``. Listed by name (independent audit, 2026-11-14) so that a new
#: table cannot pass by being a renumber record.
CURRENT_WITHOUT_DESIGN_ID: set[str] = set()  # none since the per-insert cut: they hold design ids
#: Dated records whose member ids are from the numbering of their day (24 on 2026-11-14,
#: after the ordered-length renumber). Each carries a design id, which is the join key; their
#: member ids are never re-keyed. A table not listed here that names a member by an id that is
#: not current fails :func:`test_stale_member_ids_only_in_dated_records`.
DATED = {
    "l1_domestication_edits.tsv",
    "l1_proteoform_audit.tsv",
    "l1_proteoform_decisions.tsv",
    "l1_rbp_motif_members.tsv",
    # The ranker's dated snorkel list, committed byte for byte under its pinned SHA-256
    # (l1_snorkel.COMMITTED_COPIES), so its member ids cannot be re-keyed; read on design_id.
    "l1_snorkel_designs_2026-11-14b.tsv",
    "l1_reference_signal_peptide_audit.tsv",
    "membrane-side-gate/calibration_uniprot_side.tsv",
    "membrane-side-gate/l1_extracellular_mask.tsv",
    "membrane-side-gate/l1_extracellular_mask_residues.tsv",
    "membrane-side-gate/l1_internal_member_gate.tsv",
    "membrane-side-gate/l1_internal_pick_changes.tsv",
    "membrane-side-gate/l1_internal_pick_gate.tsv",
    "membrane-side-gate/l1_snorkel_candidates.tsv",
    "membrane-side-gate/ordered-sites-2026-11-14-pr331/gate_on_ordered_sites.tsv",
    "membrane-side-gate/ordered-sites-2026-11-14/gate_on_ordered_sites.tsv",
    "nterm_res1_lane_check.tsv",
    "nterm_res1_lane_review.tsv",
    "nterm_res1_lane_screen.tsv",
    "nterm_res1_nongpcr_review.tsv",
    "nterm_res1_signalp6.tsv",
    "proteoform_decisions_applied.tsv",
    "reference_form_screen.tsv",
    "stop_codon_fix.tsv",
    "tag_refold_construct_check.tsv",
    "twist_complexity_map.tsv",
}
#: Tables kept on current member ids: the renumber stage's records plus the map and plan.
CURRENT = set(rn.RECORDS) | {
    "l1_member_renumbering.tsv",
    "reviewer-package/l1_order_plan.tsv",
    "candidate_manifest_tranches.tsv",
}


def _tables() -> dict[str, list[dict[str, str]]]:
    out = {}
    for path in sorted(HERE.rglob("*.tsv")):
        rel = path.relative_to(HERE).as_posix()
        if rel.startswith(("tests/", "review-repo/")) or "/review-repo/" in rel:
            continue
        with path.open() as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        if rows and any(MEMBER.match(v or "") for r in rows[:500] for v in r.values()):
            out[rel] = rows
    return out


def _design_cols(rows: list[dict[str, str]]) -> list[str]:
    return [c for c in rows[0] if c.endswith("design_id")]


def test_member_id_tables_carry_a_design_id():
    tables = _tables()
    current = sorted(r for r, rows in tables.items() if not _design_cols(rows) and r in CURRENT)
    assert current == sorted(CURRENT_WITHOUT_DESIGN_ID), (
        "renumber records without a design id changed; list them in CURRENT_WITHOUT_DESIGN_ID "
        f"or add a design_id column: {current}"
    )
    assert CURRENT_WITHOUT_DESIGN_ID <= set(rn.RECORDS), "only renumber records may be listed"
    missing = sorted(r for r, rows in tables.items() if not _design_cols(rows) and r not in CURRENT)
    assert missing == sorted(NO_DESIGN_ID), (
        "tables with member ids but no design id (add a design_id column, or list the stable "
        f"key in NO_DESIGN_ID): {sorted(set(missing) - set(NO_DESIGN_ID))}; "
        f"listed but now fixed: {sorted(set(NO_DESIGN_ID) - set(missing))}"
    )


def test_design_ids_name_real_designs():
    with (HERE / "candidate_manifest_tranches.tsv").open() as fh:
        known = {r["design_id"] for r in csv.DictReader(fh, delimiter="\t")}
    bad = {}
    for rel, rows in _tables().items():
        for col in _design_cols(rows):
            unknown = {r[col] for r in rows if DESIGN.match(r[col] or "")} - known
            if unknown:
                bad[f"{rel}:{col}"] = sorted(unknown)[:3]
    assert not bad, f"design ids that name no design: {bad}"


def test_current_tables_agree_with_the_renumbering():
    """The manifest keys every member by its design id, and the order plan's plate positions
    are the insert numbering's (per-insert cut: ``l1_insert_numbering.tsv``)."""
    with (HERE / "candidate_manifest_tranches.tsv").open() as fh:
        manifest = list(csv.DictReader(fh, delimiter="\t"))
    assert all(r["candidate_id"] == r["design_id"] for r in manifest), "manifest names a member id"
    with (HERE / "l1_insert_numbering.tsv").open() as fh:
        numbered = [(r["insert_id"], r["arm_id"], r["plate"], r["well"])
                    for r in csv.DictReader(fh, delimiter="\t")]  # fmt: skip
    with (HERE / "reviewer-package/l1_order_plan.tsv").open() as fh:
        plan = [(r["insert_id"], r["arm_id"], r["plate"], r["well"])
                for r in csv.DictReader(fh, delimiter="\t")]  # fmt: skip
    assert plan == numbered, "order plan disagrees with l1_insert_numbering.tsv"


def test_stale_member_ids_only_in_dated_records():
    """A table whose member ids disagree with the current numbering must be a listed dated
    record (joined on its design id); every other table is current."""
    with (HERE / "l1_member_renumbering.tsv").open() as fh:
        member = {r["design_id"]: r["candidate_id"] for r in csv.DictReader(fh, delimiter="\t")}
    stale = set()
    for rel, rows in _tables().items():
        dcols = _design_cols(rows)
        mcols = [c for c in rows[0] if c in ("candidate_id", "member_id")]
        if not dcols or not mcols:
            continue
        if any(
            r[dcols[0]] in member
            and MEMBER.match(r[mcols[0]] or "")
            and r[mcols[0]] != member[r[dcols[0]]]
            for r in rows
        ):
            stale.add(rel)
    assert stale <= DATED, f"stale member ids outside the dated records: {sorted(stale - DATED)}"
    assert not stale & CURRENT, f"current tables with stale ids: {sorted(stale & CURRENT)}"


#: Folders whose scripts join their tables to each other and to current tables on
#: ``candidate_id``; every table in them must be on the current numbering
#: (``rekey_rung_tables.DESIGN_KEYED`` keeps them there). Added 2026-11-14 after such a join named
#: the wrong protein (PSCA in ``tests/test_nterm_igk_evidence.py``).
JOINED_ON_MEMBER_ID = ("reviewer-package/nterm/",)


def test_tables_joined_on_member_id_are_current():
    with (HERE / "l1_member_renumbering.tsv").open() as fh:
        member = {r["design_id"]: r["candidate_id"] for r in csv.DictReader(fh, delimiter="\t")}
    stale = []
    for rel, rows in _tables().items():
        if not rel.startswith(JOINED_ON_MEMBER_ID) or not _design_cols(rows):
            continue
        d = _design_cols(rows)[0]
        for c in ("candidate_id", "member_id"):
            if c in rows[0] and any(
                r[d] in member and MEMBER.match(r[c] or "") and r[c] != member[r[d]] for r in rows
            ):
                stale.append(f"{rel}:{c}")
    assert not stale, f"tables joined on member id carry stale ids (run rekey_rung_tables): {stale}"
