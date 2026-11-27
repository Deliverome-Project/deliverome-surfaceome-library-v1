"""Write the L1 panel tables of record from the panels as applied on 2026-11-22.

The 2026-11-22 changes were applied to the manifest by three scripts, none of which writes the
panel tables of record (each says so):

* ``apply_panel_redraw_2026_11_22.py`` applied the isoform redraw
  (``isoform_panel_redraw_2026-11-22.tsv``, multipass split by N-terminal side) and the
  ortholog narrower feature rule (``ortholog_rule_2026-11-22.tsv``: four drawn forms restored,
  two held);
* ``refill_isoform_panel_2026_11_22.py`` refilled the six isoform slots the MANE-GTF reference
  switch emptied, continuing the draw's seeded RNG (``panel_refill_2026-11-22_changes.tsv``).

This script writes the tables that describe the panels L1 now orders:

* ``isoform_panel.tsv``, ``isoform_panel_replacements.tsv``, ``isoform_panel_summary.json``:
  the applied draw, replayed with the committed sampler exactly as the refill script replays it
  (the replay must give the 150 applied picks), with each emptied slot's row replaced by its
  refill. The replacements table is the draw's, plus one row per refill; the summary is the
  sampler's, with ``genes``, ``by_class`` and ``already_in_cohort`` recounted on the final
  panel and the refills counted separately (``refilled_after_draw``).
* ``ortholog_panel_adopted.tsv``: the 2026-11-10 adoption
  (``ortholog_panel_adopted_2026-11-10.tsv``, ``adopt_ortholog_cds_policy.py``) with each form
  the rule restored back in its replacement's place, its row taken from the drawn panel
  (``ortholog_panel_cds_eligible.tsv``). Held forms stay on their replacement.

The previous records keep dated names (``isoform_panel_2026-11-12.tsv`` etc.,
``ortholog_panel_adopted_2026-11-10.tsv``): the scripts that reproduce those draws compare
against them. Every output is checked against the manifest: the active L1 isoform and ortholog
members must be exactly the panels written here.

**Decided by** Becca Carlson, 2026-11-22 (the redraw, the refill, the narrower feature rule).
Implemented by Claude Code session l1-review-export-scope (brief-reconciliation fork),
2026-11-22. Changes no manifest row, CDS record or S3 key.

    uv run --locked python write_panels_of_record.py           # write the four tables
    uv run --locked python write_panels_of_record.py --check   # fail if they differ
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import sys
import tempfile
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import redraw_isoform_panel_2026_11_22 as rd  # noqa: E402
import refill_isoform_panel_2026_11_22 as rf  # noqa: E402
from redraw_panels import iso_keys  # noqa: E402

MANIFEST = HERE / "candidate_manifest_tranches.tsv"
ISO_PANEL = HERE / "isoform_panel.tsv"
ISO_REPLACEMENTS = HERE / "isoform_panel_replacements.tsv"
ISO_SUMMARY = HERE / "isoform_panel_summary.json"
ORT_RULE = HERE / "ortholog_rule_2026-11-22.tsv"
ORT_ADOPTED_2026_11_10 = HERE / "ortholog_panel_adopted_2026-11-10.tsv"
ORT_DRAWN = HERE / "ortholog_panel_cds_eligible.tsv"
ORT_ADOPTED = HERE / "ortholog_panel_adopted.tsv"
REFILL_REASON = (
    "its protein became its gene's reference (MANE GTF rerun 2026-11-22); "
    "refilled after the draw by the same seeded RNG"
)


def read(path: Path) -> list[dict]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def tsv_text(rows: list[dict], cols: list[str]) -> str:
    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=cols, delimiter="\t", lineterminator="\n")
    w.writeheader()
    for r in rows:
        w.writerow({c: r.get(c, "") for c in cols})
    return buf.getvalue()


def isoform_panel(work: Path, manifest: list[dict]) -> tuple[list[dict], list[dict], dict]:
    """(panel rows, replacement rows, summary) of the applied draw plus refills."""
    extra = json.loads(rd.OUT_SUMMARY.read_text())["rejected_as_unbuildable_in_redraw"]
    panel, repl, summary = rd.draw(work / "draw", split=True, extra=extra)
    applied = {r["form_id"] for r in rd.read(rf.REDRAW) if r["status"] != "dropped"}
    assert {rf.form_id(r) for r in panel} == applied, "the replay is not the applied draw"
    _, calls, _ = rf.replay(work / "replay", extra)
    candidates = {rf.form_id(c): c for call in calls for c in call["stratum"]}

    by_id = {m["design_id"]: m for m in manifest}
    refills = [c for c in read(rf.CHANGES) if c["change"] == "added"]
    panel = list(panel)
    extra_repl = []
    for c in refills:
        emptied = by_id[c["refills"]]
        slot = [i for i, r in enumerate(panel) if emptied["key"] in iso_keys(r)]
        assert len(slot) == 1, f"{c['refills']}: {len(slot)} panel rows"
        old = panel[slot[0]]
        new = {**candidates[c["form_id"]], "panel_tier": old["panel_tier"]}
        assert f"{new['panel_tier']} / {new['topology_class']}" == c["stratum"], c["form_id"]
        panel[slot[0]] = new
        extra_repl.append(
            {
                "panel_tier": old["panel_tier"],
                "topology_class": old["topology_class"],
                "rejected": rf.form_id(old),
                "reason": REFILL_REASON,
                "replacement": c["form_id"],
                "candidates": c["detail"].split("draw (")[1].split(" candidates")[0],
            }
        )
    panel.sort(key=lambda r: (r["panel_tier"], r["topology_class"], r["gene"]))
    summary = {
        **summary,
        "by_class": dict(Counter(r["topology_class"] for r in panel).most_common()),
        "genes": len({r["gene"] for r in panel}),
        "already_in_cohort": sum(1 for r in panel if r["in_current_cohort"] == "True"),
        "refilled_after_draw": len(refills),
        "record": (
            "isoform_panel_redraw_2026-11-22.tsv as applied (apply_panel_redraw_2026_11_22.py) "
            "plus panel_refill_2026-11-22_changes.tsv; written by write_panels_of_record.py"
        ),
    }
    return panel, list(repl) + extra_repl, summary


def ortholog_panel() -> list[dict]:
    """The 2026-11-10 adoption with the rule's restored forms back in place."""
    adopted = read(ORT_ADOPTED_2026_11_10)
    drawn = {r["form_id"]: r for r in read(ORT_DRAWN)}
    pos = {r["form_id"]: i for i, r in enumerate(adopted)}
    restored = 0
    for r in read(ORT_RULE):
        if not r["panel_change"].startswith("restore "):
            continue
        d, i = drawn[r["form_id"]], pos[r["adopted_form_id"]]
        old = adopted[i]
        assert (d["species"], d["band"], d["topology_class"]) == (
            old["species"],
            old["band"],
            old["topology_class"],
        ), f"{r['form_id']}: restored form is not in its replacement's stratum"
        adopted[i] = {
            "species": d["species"],
            "band": d["band"],
            "topology_class": d["topology_class"],
            "human_symbol": d["human_symbol"],
            "form_id": d["form_id"],
            "picked_acc": d["picked_acc"],
            "picked_score": d["picked_score"],
            "sequence_sha256": d["sequence_sha256"],
            "human_reference_sha256": d["human_reference_sha256"],
            "human_ensembl_gene": d["human_ensembl_gene"],
            "ortholog_ensembl_gene": d["ortholog_ensembl_gene"],
            "cds_transcript": d["ensembl_transcripts"],
            "selected_cds_sha256": d["selected_cds_sha256"],
            "cds_source": "final-v2 selected_cds",
            "replaces": "",
            "decided_by": r["decided_by"],
            "implemented_by": r["implemented_by"],
            "decision_reason": (
                f"Restored by the narrower feature rule: the drawn form loses no "
                f"display-deciding feature ({r['human_features_lost']} only); "
                f"{r['adopted_form_id']} (L1 {r['replacement_design_id']}) is retired"
            ),
        }
        restored += 1
    assert restored == 4, f"{restored} forms restored, expected 4"
    return adopted


def check_manifest(manifest: list[dict], iso: list[dict], ort: list[dict]) -> None:
    active = [m for m in manifest if m["library"] == "l1" and not m["retired_reason"]]
    iso_members = {m["key"] for m in active if m["set"] == "isoform"}
    iso_hit = [sorted(iso_keys(r) & iso_members) for r in iso]
    assert all(len(h) == 1 for h in iso_hit), "an isoform panel row matches no single member"
    assert {h[0] for h in iso_hit} == iso_members, "active L1 isoform members != the panel"
    # An ortholog member is keyed ``ort:<form id>``, except where its L2 row was keyed by
    # another accession of the same protein (two of the restored forms); those match by
    # protein. Protein alone is not enough: the four start-trimmed cyno forms
    # (``ortholog_start_trims.tsv``) are built shorter than the panel's drawn protein.
    import hashlib

    from deliverome_analysis.cloning import library_design as ld

    designed = ld.designed_proteins(MANIFEST, library=ld.LIBRARY_L1)
    ort_active = {m["key"]: m for m in active if m["set"] == "ortholog"}
    by_protein = {
        hashlib.sha256(designed[m["candidate_id"]]["sequence"].encode()).hexdigest(): k
        for k, m in ort_active.items()
        if m["candidate_id"] in designed
    }
    matched = [
        f"ort:{r['form_id']}"
        if f"ort:{r['form_id']}" in ort_active
        else by_protein.get(r["sequence_sha256"], "")
        for r in ort
    ]
    assert "" not in matched, "an ortholog panel row matches no active member"
    assert sorted(matched) == sorted(ort_active), "active L1 ortholog members != the panel"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--check", action="store_true", help="fail if the tables differ")
    ap.add_argument("--work", type=Path, default=None)
    a = ap.parse_args(argv)
    work = a.work or Path(tempfile.mkdtemp(prefix="panels-of-record-"))
    manifest = read(MANIFEST)
    iso, repl, summary = isoform_panel(work, manifest)
    ort = ortholog_panel()
    check_manifest(manifest, iso, ort)
    iso_cols = list(read(HERE / "isoform_panel_2026-11-12.tsv")[0])
    ort_cols = list(read(ORT_ADOPTED_2026_11_10)[0])
    repl_cols = ["panel_tier", "topology_class", "rejected", "reason", "replacement", "candidates"]
    out = {
        ISO_PANEL: tsv_text(iso, iso_cols),
        ISO_REPLACEMENTS: tsv_text(repl, repl_cols),
        ISO_SUMMARY: json.dumps(summary, indent=2) + "\n",
        ORT_ADOPTED: tsv_text(ort, ort_cols),
    }
    differs = [p.name for p, t in out.items() if not p.exists() or p.read_text() != t]
    print(
        f"isoform panel: {len(iso)} forms, {summary['genes']} genes, {summary['by_class']}; "
        f"{summary['replaced']} replaced in the draw, {summary['refilled_after_draw']} refilled"
    )
    print(f"ortholog panel: {len(ort)} forms, {sum(bool(r['replaces']) for r in ort)} replaced")
    if a.check:
        print("differs: " + ", ".join(differs) if differs else "all four tables match")
        return 1 if differs else 0
    for p, t in out.items():
        p.write_text(t)
    print("wrote " + ", ".join(p.name for p in out))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
