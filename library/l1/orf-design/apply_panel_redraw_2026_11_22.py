"""Apply the 2026-11-22 isoform redraw and ortholog feature rule to the L1 manifest.

**Decided by Becca Carlson, 2026-11-22** ("yes we should always rebuild with the latest"; the
multipass split; the narrower feature rule). Implemented by Claude Code session
l1-review-export-scope (worker fork), 2026-11-22.

Inputs are the two committed proposals:

* ``isoform_panel_redraw_2026-11-22.tsv`` (``redraw_isoform_panel_2026_11_22.py``, the shared
  seeded stream, the method of record). Each dropped L1 isoform member is retired
  ``dropped_from_panel:isoform_redraw_2026-11-22 (Becca 2026-11-22)``. Each added form is
  moved from its L2 row (retired ``moved_to:<new id>``) to a new L1 row, as
  ``redraw_panels.py`` does: ``design_id`` = ``candidate_id`` = the new id,
  ``prior_candidate_id`` = the L2 id, ``panel`` = ``isoform_panel``.
* ``ortholog_rule_2026-11-22.tsv`` (``ortholog_rule_2026_11_22.py``). A drawn form the narrower
  rule keeps is restored when an active L2 member already has its exact protein and CDS (as
  ``adopt_ortholog_cds_policy.py`` promotes one); its replacement is retired
  ``replaced_by_original:ortholog_rule_2026-11-22``. ``HOLD`` lists the forms left on their
  replacement, with the reason.

New ids come from ``member_id_allocations.next_l1_member_number`` on the manifest as it stands,
isoforms first in table order, then orthologs; no existing row is renumbered, and a retired id
is never reused. The script is idempotent: a change already in the manifest is verified, not
repeated. ``panel_redraw_2026-11-22_changes.tsv`` lists every change. Not touched: the panel
tables of record (``isoform_panel.tsv``, ``ortholog_panel_adopted.tsv``, ``panel_inputs.json``
pins), the designed-CDS and construct records (every promoted key is already in both), and the
rung tables. Member ids (``l1-mNNNN``) follow from ``renumber_library_members.py``.

    uv run --locked python apply_panel_redraw_2026_11_22.py           # dry run: print the changes
    uv run --locked python apply_panel_redraw_2026_11_22.py --write   # write the manifest
    uv run --locked python apply_panel_redraw_2026_11_22.py --check   # fail unless applied
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from member_id_allocations import RENUMBER_HINT, next_l1_member_number  # noqa: E402
from redraw_panels import iso_keys  # noqa: E402

MANIFEST = HERE / "candidate_manifest_tranches.tsv"
ISOFORM = HERE / "isoform_panel_redraw_2026-11-22.tsv"
ORTHOLOG = HERE / "ortholog_rule_2026-11-22.tsv"
CDS_ELIGIBLE = HERE / "ortholog_panel_cds_eligible.tsv"
CHANGES = HERE / "panel_redraw_2026-11-22_changes.tsv"
ISO_RETIRED = "dropped_from_panel:isoform_redraw_2026-11-22 (Becca 2026-11-22)"
ORT_RETIRED = "replaced_by_original:ortholog_rule_2026-11-22"
DECIDED_BY = "Becca Carlson, 2026-11-22"
IMPLEMENTED_BY = "Claude Code session l1-review-export-scope (worker fork), 2026-11-22"
#: Forms the narrower rule would restore that stay on their replacement, and why.
HOLD = {
    "cyno:SLC8A1:ENSMFAP00000038334": (
        "held pending Becca's decision: keeps SP, N-side, TM1 and anchor, but lacks TM6-TM10 "
        "and the exchanger domain (coordinator, 2026-11-22)"
    ),
    "cyno:ADGRB2:A0A2K5USV4": (
        "held: no designed-CDS record or manifest member encodes this form, so restoring it "
        "needs a new CDS and construct record (not published here)"
    ),
}
CHANGE_COLUMNS = [
    "panel",
    "change",
    "candidate_id",
    "key",
    "prior_candidate_id",
    "form_id",
    "detail",
    "decided_by",
    "implemented_by",
]


def read(path: Path) -> list[dict]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def active(r: dict) -> bool:
    return not (r.get("retired_reason") or "").strip()


#: Retirements recorded after this redraw by later, separate decisions; a row so retired still
#: counts as this redraw's change having been made (e.g. the MANE-GTF reference switch,
#: ``apply_mane_gtf.py``, retired six panel members that became identical to a reference).
LATER_RETIREMENTS = ("duplicates_reference:",)


def made(r: dict) -> bool:
    return active(r) or r["retired_reason"].startswith(LATER_RETIREMENTS)


def retire(row: dict, reason: str) -> bool:
    """Set ``reason``; True if it changed. Refuses to overwrite a different retirement."""
    if row["retired_reason"] == reason:
        return False
    if row["retired_reason"]:
        raise SystemExit(f"{row['candidate_id']}: already retired ({row['retired_reason']})")
    row["retired_reason"] = reason
    return True


def promote(rows: list[dict], src: dict, panel: str, next_n: list[int]) -> tuple[dict, bool]:
    """The L1 row moved from L2 row ``src``: the existing one if already moved, else a new one."""
    done = [
        r
        for r in rows
        if r["library"] == "l1" and made(r) and r["prior_candidate_id"] == src["candidate_id"]
    ]
    if done:
        assert len(done) == 1 and done[0]["key"] == src["key"], src["candidate_id"]
        return done[0], False
    assert src["library"] == "l2", f"{src['candidate_id']}: not an L2 row"
    new_id = f"l1-p{next_n[0]:04d}"
    next_n[0] += 1
    if active(src):
        retire(src, f"moved_to:{new_id}")
    # else: the L2 row already moved once, to an L1 id since retired; that record stays, and
    # the new L1 row names the L2 row as its prior.
    row = {
        **src,
        "candidate_id": new_id,
        "design_id": new_id,
        "library": "l1",
        "panel": panel,
        "prior_candidate_id": src["candidate_id"],
        "retired_reason": "",
    }
    rows.append(row)
    return row, True


def apply(rows: list[dict]) -> tuple[list[dict], int]:
    """Apply both proposals to ``rows`` in place. Returns (change log, number of edits)."""
    next_n = [next_l1_member_number(rows)]
    by_id = {r["candidate_id"]: r for r in rows}
    log: list[dict] = []
    edits = 0

    def note(panel, change, row, form_id, detail):
        log.append(
            {
                "panel": panel,
                "change": change,
                "candidate_id": row["candidate_id"] if row else "",
                "key": row["key"] if row else "",
                "prior_candidate_id": row.get("prior_candidate_id", "") if row else "",
                "form_id": form_id,
                "detail": detail,
                "decided_by": DECIDED_BY,
                "implemented_by": IMPLEMENTED_BY,
            }
        )

    # ---- isoform panel ------------------------------------------------------------------
    plan = read(ISOFORM)
    for p in plan:
        if p["status"] != "dropped":
            continue
        row = by_id[p["design_id"]]
        assert row["set"] == "isoform" and row["library"] == "l1", p["design_id"]
        edits += retire(row, ISO_RETIRED)
        note("isoform", "retired", row, p["form_id"], ISO_RETIRED)
    for p in plan:
        if p["status"] != "added":
            continue
        keys = iso_keys(p)
        srcs = [
            r
            for r in rows
            if r["library"] == "l2"
            and r["set"] == "isoform"
            and r["key"] in keys
            and (active(r) or r["retired_reason"].startswith("moved_to:"))
        ]
        assert len(srcs) == 1, f"{p['form_id']}: {len(srcs)} L2 rows"
        src = srcs[0]
        earlier = "" if active(src) else src["retired_reason"].split(":", 1)[1]
        row, changed = promote(rows, src, "isoform_panel", next_n)
        by_id[row["candidate_id"]] = row
        edits += changed
        detail = f"moved from L2 {src['candidate_id']}"
        if earlier and earlier != row["candidate_id"]:
            detail += f"; its earlier L1 id {earlier} stays retired"
        note("isoform", "added", row, p["form_id"], detail)
    for p in plan:
        if p["status"] == "kept":
            row = by_id[p["design_id"]]
            assert made(row) and row["library"] == "l1", p["design_id"]

    # ---- ortholog rule ------------------------------------------------------------------
    eligible = {r["form_id"]: r for r in read(CDS_ELIGIBLE)}
    for o in read(ORTHOLOG):
        if not o["panel_change"].startswith(("restore ", "held on ")):
            continue
        f = o["form_id"]
        repl = by_id[o["replacement_design_id"]]
        if f in HOLD:
            assert active(repl), f"{repl['candidate_id']}: held replacement is retired"
            note("ortholog", "held", repl, f, HOLD[f])
            continue
        l2 = by_id[eligible[f]["existing_candidate_ids"].strip()]
        row, changed = promote(rows, l2, "ortholog_panel", next_n)
        by_id[row["candidate_id"]] = row
        edits += changed
        note(
            "ortholog", "added", row, f, f"restored drawn form; moved from L2 {l2['candidate_id']}"
        )
        edits += retire(repl, ORT_RETIRED)
        note("ortholog", "retired", repl, o["adopted_form_id"], ORT_RETIRED)

    ids = [r["candidate_id"] for r in rows]
    assert len(ids) == len(set(ids)), "candidate ids are not unique"
    return log, edits


def verify_records(rows: list[dict], log: list[dict]) -> None:
    """Every added member's key is in the designed-CDS and construct records."""
    from deliverome_analysis import library_design as ld
    from deliverome_analysis.s3 import processed

    with open(processed.download(ld.CDS_KEY)) as fh:
        cds = json.load(fh)
    with open(processed.download(ld.CONSTRUCTS_KEY)) as fh:
        constructs = json.load(fh)
    missing = [
        c["key"]
        for c in log
        if c["change"] == "added" and (c["key"] not in cds or c["key"] not in constructs)
    ]
    if missing:
        raise SystemExit(f"added members with no CDS or construct record: {missing}")


def text(rows: list[dict], cols: list[str]) -> str:
    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=cols, delimiter="\t", lineterminator="\n")
    w.writeheader()
    for r in rows:
        w.writerow({c: r.get(c, "") for c in cols})
    return buf.getvalue()


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    mode = ap.add_mutually_exclusive_group()
    mode.add_argument("--write", action="store_true")
    mode.add_argument("--check", action="store_true", help="fail unless already applied")
    a = ap.parse_args(argv)
    with MANIFEST.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        cols, rows = list(reader.fieldnames or []), list(reader)
    log, edits = apply(rows)
    verify_records(rows, log)
    counts = {
        (c["panel"], c["change"]): sum(
            1 for x in log if (x["panel"], x["change"]) == (c["panel"], c["change"])
        )
        for c in log
    }
    for (panel, change), n in sorted(counts.items()):
        print(f"{panel:9s} {change:8s} {n}")
    print(f"{edits} manifest edits needed")
    if a.check:
        same = edits == 0 and CHANGES.exists() and CHANGES.read_text() == text(log, CHANGE_COLUMNS)
        print("applied and matches the record" if same else "NOT applied, or the log differs")
        return 0 if same else 1
    if not a.write:
        print("dry run -- re-run with --write")
        return 0
    MANIFEST.write_text(text(rows, cols))
    CHANGES.write_text(text(log, CHANGE_COLUMNS))
    print(f"wrote {MANIFEST.name}, {CHANGES.name}")
    print(RENUMBER_HINT)
    return 0


if __name__ == "__main__":
    sys.exit(main())
