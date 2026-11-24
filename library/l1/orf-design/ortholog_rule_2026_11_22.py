"""The ortholog feature rule narrowed to display-deciding features (proposed, 2026-11-22).

**Decided by Becca Carlson, 2026-11-22**: "shouldn't some features be ok to use as long as they
are compatible with display? we are fine with diversity". Implemented by Claude Code session
l1-review-export-scope (worker fork), 2026-11-22.

The correspondence review (``ORTHOLOG_CORRESPONDENCE_REVIEW.md``) replaced a drawn ortholog
form when any human feature (signal peptide, any TM helix, N-terminal side, GPI anchor or a
Pfam domain across an insertion or deletion) had no counterpart; eleven of the 200 drawn forms
were replaced on adoption (``ortholog_panel_adoption_judgements.tsv``). Under the narrower rule
a form is replaced only when it loses a feature that decides display:

* the signal peptide (``signal_peptide_lost``);
* the N-terminal side (``n_terminal_orientation_changed``);
* the first TM helix (``TM1`` among ``human_tm_helices_unmatched``);
* the membrane anchor (``gpi_anchor_lost``, or every human TM helix gone).

A later helix or a Pfam domain is no longer a reason. This script applies that rule to every
one of the 200 drawn forms in ``ortholog_correspondence_review.tsv`` (the evidence is reused
unchanged; nothing is re-predicted), and writes, for each form the adoption replaced, whether
the drawn form would now be kept, and the resulting diff of the adopted ortholog panel.

Design ids: a replacement that leaves carries its L1 id (which would be retired
``dropped_from_panel``); a restored form carries its L2 member id when it has one and a
*proposed* next-free L1 id, numbered after the isoform redraw's proposals
(``isoform_panel_redraw_2026-11-22.tsv``) so the two proposals never claim the same number.
Nothing is allocated and the manifest, CDS records and rung tables are not touched.

    uv run --locked python ortholog_rule_2026_11_22.py           # write ortholog_rule_2026-11-22.tsv
    uv run --locked python ortholog_rule_2026_11_22.py --check   # replay, write nothing
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import sys
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
REVIEW = HERE / "ortholog_correspondence_review.tsv"
JUDGEMENTS = HERE / "ortholog_panel_adoption_judgements.tsv"
#: The 2026-11-10 adoption this rule is evaluated against (``adopt_ortholog_cds_policy.py``);
#: ``ortholog_panel_adopted.tsv`` itself now carries the rule's restorations.
ADOPTED = HERE / "ortholog_panel_adopted_2026-11-10.tsv"
CDS_ELIGIBLE = HERE / "ortholog_panel_cds_eligible.tsv"
MANIFEST = HERE / "candidate_manifest_tranches.tsv"
ISOFORM_REDRAW = HERE / "isoform_panel_redraw_2026-11-22.tsv"
OUT = HERE / "ortholog_rule_2026-11-22.tsv"
OUT_SUMMARY = HERE / "ortholog_rule_2026-11-22_summary.json"
DECIDED_BY = "Becca Carlson, 2026-11-22 (narrower feature rule: display-deciding features only)"
IMPLEMENTED_BY = "Claude Code session l1-review-export-scope (worker fork), 2026-11-22"
COLUMNS = [
    "form_id",
    "species",
    "band",
    "topology_class",
    "earlier_rule",
    "adoption_decision",
    "adopted_form_id",
    "human_features_lost",
    "display_features_lost",
    "narrow_rule",
    "panel_change",
    "restored_form_design_id",
    "restored_form_id_status",
    "replacement_design_id",
    "replacement_id_status",
    "restored_form_cds",
    "decided_by",
    "implemented_by",
]


def read(path: Path) -> list[dict]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def true(v: str) -> bool:
    return v.strip() == "True"


def display_losses(r: dict) -> list[str]:
    """The display-deciding human features the form lacks (empty: keep under the narrow rule)."""
    out = []
    if true(r["signal_peptide_lost"]):
        out.append("signal peptide")
    if true(r["n_terminal_orientation_changed"]):
        out.append(f"N-terminal side ({r['human_n_orientation']} -> {r['ortholog_n_orientation']})")
    unmatched = [s.split(":")[0] for s in r["human_tm_helices_unmatched"].split(";") if s]
    if "TM1" in unmatched:
        out.append("first TM helix")
    human_tm = int(r["human_tm_helices"] or 0)
    if true(r["gpi_anchor_lost"]) or (human_tm and len(unmatched) == human_tm):
        out.append("membrane anchor")
    return out


def all_losses(r: dict) -> list[str]:
    flags = (
        "signal_peptide_lost",
        "tm_helix_lost",
        "n_terminal_orientation_changed",
        "gpi_anchor_lost",
        "pfam_domain_lost_with_indel",
    )
    out = [f for f in flags if true(r[f])]
    if r["human_tm_helices_unmatched"]:
        out.append(f"unmatched {r['human_tm_helices_unmatched']}")
    return out


def build() -> tuple[list[dict], dict]:
    review = {r["form_id"]: r for r in read(REVIEW)}
    judged = {j["form_id"]: j for j in read(JUDGEMENTS)}
    adopted = {r["form_id"] for r in read(ADOPTED)}
    eligible = {r["form_id"]: r for r in read(CDS_ELIGIBLE)}
    manifest = read(MANIFEST)
    by_key: dict[str, list[dict]] = {}
    for m in manifest:
        by_key.setdefault(m["key"], []).append(m)

    # Next free L1 number after both the manifest and the isoform redraw's proposals.
    from member_id_allocations import next_l1_member_number

    next_free = next_l1_member_number(manifest)
    if ISOFORM_REDRAW.exists():
        iso = [
            int(r["design_id"].split("-p")[1])
            for r in read(ISOFORM_REDRAW)
            if r["status"] == "added" and r["design_id"]
        ]
        next_free = max([next_free - 1, *iso]) + 1

    rows: list[dict] = []
    for f, r in sorted(review.items()):
        j = judged.get(f)
        lost = display_losses(r)
        row = {
            "form_id": f,
            "species": r["species"],
            "band": r["band"],
            "topology_class": r["topology_class"],
            "earlier_rule": r["recommendation"],
            "adoption_decision": j["decision"] if j else "keep (not flagged)",
            "adopted_form_id": j["adopted_form_id"] if j else f,
            "human_features_lost": ";".join(all_losses(r)),
            "display_features_lost": ";".join(lost),
            "narrow_rule": "replace" if lost else "keep",
            "decided_by": DECIDED_BY,
            "implemented_by": IMPLEMENTED_BY,
        }
        assert (f in adopted) == (not j or j["decision"] != "replace"), f
        if j and j["decision"] == "replace":
            if lost:
                row["panel_change"] = "none: still replaced"
            else:
                row["panel_change"] = f"restore {f}; drop {j['adopted_form_id']}"
                repl_now = [
                    m
                    for m in by_key.get(f"ort:{j['adopted_form_id']}", [])
                    if m["library"] == "l1"
                    and m["retired_reason"] in ("", "replaced_by_original:ortholog_rule_2026-11-22")
                ]
                assert len(repl_now) == 1, f"{j['adopted_form_id']}: no single L1 row"
                repl = repl_now
                row["replacement_design_id"] = repl[0]["design_id"]
                row["replacement_id_status"] = "existing; would be retired dropped_from_panel"
                e = eligible[f]
                row["restored_form_cds"] = (
                    f"{e['cds_eligibility_status']} ({e['encoding_source']}, "
                    f"{e['ensembl_transcripts']})"
                )
                l2 = e["existing_candidate_ids"].strip()
                now = [
                    m
                    for m in manifest
                    if m["library"] == "l1"
                    and not m["retired_reason"]
                    and l2
                    and m["prior_candidate_id"] == l2
                ]
                if now:
                    row["restored_form_design_id"] = now[0]["design_id"]
                    row["restored_form_id_status"] = (
                        f"allocated by apply_panel_redraw_2026_11_22.py; moved from L2 member {l2}"
                    )
                    row["replacement_id_status"] = (
                        f"retired ({repl_now[0]['retired_reason']})"
                        if repl_now and repl_now[0]["retired_reason"]
                        else "existing; would be retired"
                    )
                else:
                    from apply_panel_redraw_2026_11_22 import HOLD

                    if f in HOLD:
                        row["panel_change"] = f"held on {j['adopted_form_id']}: {HOLD[f]}"
                        row["replacement_id_status"] = "existing; stays active (held)"
                        row["restored_form_id_status"] = "none: held, so no id is allocated"
                    else:
                        row["restored_form_design_id"] = f"l1-p{next_free:04d}"
                        row["restored_form_id_status"] = (
                            "proposed next free L1 id (not allocated)"
                            + (
                                f"; moved from L2 member {l2}"
                                if l2
                                else "; no L2 member holds this form, so its CDS member is new"
                            )
                        )
                        next_free += 1
        elif lost and j:
            # An explicit keep stands: the narrow rule never overrides a recorded decision.
            row["panel_change"] = (
                f"none: kept by the recorded decision ({j['decided_by']}); "
                "the narrow rule alone would replace it"
            )
        elif lost:
            row["panel_change"] = "would now be replaced (not flagged by the earlier review)"
        else:
            row["panel_change"] = "none"
        rows.append(row)

    replaced = [r for r in rows if r["adoption_decision"] == "replace"]
    restored = [r for r in replaced if r["narrow_rule"] == "keep"]
    summary = {
        "decided_by": DECIDED_BY,
        "implemented_by": IMPLEMENTED_BY,
        "status": "proposed; manifest, CDS records and rung tables unchanged",
        "drawn_forms": len(rows),
        "earlier_rule": dict(Counter(r["earlier_rule"] for r in rows)),
        "narrow_rule": dict(Counter(r["narrow_rule"] for r in rows)),
        "replaced_on_adoption": len(replaced),
        "replaced_and_still_replaced": sorted(
            r["form_id"] for r in replaced if r["narrow_rule"] == "replace"
        ),
        "replaced_and_now_kept": sorted(r["form_id"] for r in restored),
        "kept_on_adoption_but_narrow_rule_replaces": sorted(
            r["form_id"]
            for r in rows
            if r["adoption_decision"] != "replace" and r["narrow_rule"] == "replace"
        ),
        "panel_diff": {
            "kept": len(rows) - sum(r["panel_change"].startswith("restore") for r in restored),
            "dropped": sum(r["panel_change"].startswith("restore") for r in restored),
            "added": sum(r["panel_change"].startswith("restore") for r in restored),
        },
        "held": sorted(r["form_id"] for r in restored if r["panel_change"].startswith("held")),
        "retired_replacements": sorted(
            r["replacement_design_id"] for r in restored if r["panel_change"].startswith("restore")
        ),
        "restored_l1_ids": sorted(
            r["restored_form_design_id"]
            for r in restored
            if r["restored_form_id_status"].startswith("allocated")
        ),
    }
    return rows, summary


def tsv_text(rows: list[dict]) -> str:
    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
    w.writeheader()
    for r in rows:
        w.writerow({c: r.get(c, "") for c in COLUMNS})
    return buf.getvalue()


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--check", action="store_true", help="replay; fail if the outputs differ")
    a = ap.parse_args(argv)
    rows, summary = build()
    text, js = tsv_text(rows), json.dumps(summary, indent=2) + "\n"
    print(js)
    if a.check:
        same = OUT.exists() and OUT.read_text() == text and OUT_SUMMARY.read_text() == js
        print("matches the record" if same else "DIFFERS from the record")
        return 0 if same else 1
    OUT.write_text(text)
    OUT_SUMMARY.write_text(js)
    print(f"wrote {OUT.name}, {OUT_SUMMARY.name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
