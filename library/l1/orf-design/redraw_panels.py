"""Redraw the L1 isoform panel by rule, and bring the manifest into line.

**Isoform-only since 2026-11-12.** The L1 ortholog panel is the Compara 116 adoption of
2026-11-10 (``adopt_ortholog_cds_policy.py``), the ortholog record of truth; rerunning the
original ortholog sampler here would undo it. This script no longer runs or reconciles the
ortholog panel, keeps the ortholog rows of ``panel_rejections.tsv`` as they are, and asserts
that every active L1 ortholog row of the manifest is byte-identical before and after.

**The cohort gate** (Becca Carlson, 2026-11-12, "Gate cohort, redraw"; implemented by the
Claude Code session l1-reviewer-package (subagent), PR #283): the isoform sampler reads
``isoform_gate.tsv`` (``gate_isoform_cohort.py``) through its validity hook, so failing
picks are replaced and failing forms are never drawn. ``panel_changes.tsv`` carries
``date`` / ``decided_by`` / ``implemented_by`` / ``reason`` per change (the manifest schema
has no attribution columns).

The panels are seeded, class-stratified draws (``library/l1/proteoform-selection/
sample_isoform_panel.py``, ``sample_ortholog_panel.py``). Both now read a pinned input
snapshot from S3 and, run with ``--original``, reproduce the published 150 and 200 byte for
byte. A pick the library cannot use is replaced by continuing the same seeded draw within
its stratum (``panel_draw.py``), so every valid pick keeps its place.

A pick is rejected when:

| Reason | Detected by |
|---|---|
| it is its gene's selected reference form | the isoform sampler (reference_selection.tsv) |
| its protein is identical to human (correspondence 1.0) | the ortholog sampler |
| no transcript encodes it, so it has no designed CDS | this script: no manifest key |
| its designed protein duplicates another active member's | this script, on the translated CDS |

The last two are only knowable after the draw, so this script writes them to
``panel_rejections.tsv`` and reruns both samplers until a pass adds no rejection. Each
pass is deterministic, so the whole loop is.

The manifest is then reconciled to the redrawn panels without renumbering: a pick that is
already an L1 member keeps its id; a pick that was an L2 member gets the next free L1 id
(in panel order) and its L2 row is retired ``moved_to:<id>``; an L1 panel member no longer
drawn is retired ``dropped_from_panel``. ``panel_changes.tsv`` lists every change.

    uv run --locked python redraw_panels.py            # dry run
    uv run --locked python redraw_panels.py --write    # write tables + manifest
    uv run --locked python redraw_panels.py --write --publish   # and the panel tables on S3
"""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import proteoform_decisions as pd_
from member_id_allocations import RENUMBER_HINT, next_l1_member_number, to_design_space

from deliverome_analysis.cloning import library_design as ld
from deliverome_analysis.cloning import seqtools
from deliverome_analysis.core.s3 import processed

HERE = Path(__file__).parent
REPO = HERE.parents[2]
SAMPLERS = REPO / "library/l1/proteoform-selection"
MANIFEST = HERE / "candidate_manifest_tranches.tsv"
REJECTIONS = HERE / "panel_rejections.tsv"
CHANGES = HERE / "panel_changes.tsv"
GATE = HERE / "isoform_gate.tsv"
#: sampler output -> committed file. This is the 2026-11-12 method; its draw is kept under dated
#: names because ``isoform_panel*.tsv/json`` now hold the applied 2026-11-22 redraw plus refills
#: (``write_panels_of_record.py``). A --write here no longer overwrites that record.
PANEL_FILES = {
    "isoform_panel.tsv": "isoform_panel_2026-11-12.tsv",
    "isoform_panel_replacements.tsv": "isoform_panel_replacements_2026-11-12.tsv",
    "isoform_panel_summary.json": "isoform_panel_summary_2026-11-12.json",
}
DATE = "2026-11-12"
DECIDED_BY = "Becca Carlson, 2026-11-12 (Gate cohort, redraw)"
IMPLEMENTED_BY = "Claude Code session l1-reviewer-package (subagent), PR #283"
PRIOR_ATTRIBUTION = {
    "date": "2026-11-01",
    "decided_by": "Becca Carlson; implementer not recorded (before the 2026-11-10 "
    "attribution convention)",
    "implemented_by": "not recorded",
    "reason": "replacement drawn by rule",
}
CHANGE_COLUMNS = [
    "change",
    "candidate_id",
    "key",
    "detail",
    "date",
    "decided_by",
    "implemented_by",
    "reason",
]
DESIGN_PREFIX = "sprint4/s4e9-orf-design-fg-panel/library-design/"
MAX_PASSES = 10
REJECTION_COLUMNS = [
    "panel",
    "form_id",
    "reason",
    "decided_by",
    "implemented_by",
    "alternative",
    "sources",
]
RULE_ROW = {
    "decided_by": "rule (redraw_panels.py; Becca Carlson, 2026-11-01)",
    "implemented_by": "redraw_panels.py",
    "alternative": "",
    "sources": "",
}
#: Domain checks made here on the two fragments named first (UniProt REST features read
#: 2026-11-11; each sequence identical to the pinned UniProt 2026-10-11 table).
VERIFIED = {
    "ABCC3:O15438-5": "verified by l1-reviewer-package: O15438-5 is canonical O15438 1-546 "
    "plus a 26-aa novel tail (VSP_039041); it lacks UniProt's ABC transporter 1 (629-851), ABC "
    "transmembrane type-1 2 (971-1252) and ABC transporter 2 (1291-1523) domains and TM "
    "helices 11-17, so it has neither nucleotide-binding domain nor the second TMD",
    "OPN3:ENST00000469376": "verified by l1-reviewer-package: the 189-aa form is canonical "
    "Q9H1Y3 1-181 plus 8 divergent residues; it keeps TM1-4 (41-177) and lacks TM5-7 "
    "(202-312), including the retinal-binding Lys299 (N6-retinylidene lysine)",
}


def explicit_rejections() -> dict[tuple[str, str], dict]:
    """Isoform forms rejected by decision (``proteoform_decisions_applied.tsv``)."""
    out = {}
    for key, r in pd_.by_key().items():
        if r["decision"] != "redraw" or not key.startswith("iso:"):
            continue
        fid = key.split(":", 1)[1]
        a = pd_.attribution(key)
        src = a["sources"] + (f"; {VERIFIED[fid]}" if fid in VERIFIED else "")
        out[("isoform", fid)] = {
            "reason": f"explicit rejection: {a['reason']}",
            "decided_by": a["decided_by"],
            "implemented_by": a["implemented_by"],
            "alternative": a["alternative"],
            "sources": src,
        }
    return out


def read(path) -> list[dict]:
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def iso_id(r: dict) -> str:
    return f"{r['gene']}:{r['uniprot'].strip() or r['ensembl_transcripts']}"


def ort_id(r: dict) -> str:
    return f"{r['species']}:{r['human_symbol']}:{r['picked_acc']}"


def iso_keys(r: dict) -> set[str]:
    """Manifest keys an isoform panel row may resolve to (split_library_tranches rule)."""
    keys = {f"iso:{r['gene']}:{r['uniprot']}"} if r["uniprot"].strip() else set()
    keys |= {f"iso:{r['gene']}:{t.strip()}" for t in r["ensembl_transcripts"].split(";") if t}
    return keys


def run_samplers(workdir: Path, rejections: Path) -> None:
    subprocess.run(
        [
            sys.executable,
            str(SAMPLERS / "sample_isoform_panel.py"),
            "--snapshot",
            "--reject",
            str(rejections),
            "--gate",
            str(GATE),
        ],
        cwd=workdir,
        check=True,
        capture_output=True,
        env={"PYTHONPATH": str(SAMPLERS), **__import__("os").environ},
    )


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--write", action="store_true")
    ap.add_argument("--publish", action="store_true")
    args = ap.parse_args(argv)

    with open(processed.download(ld.CDS_KEY)) as fh:
        cds = json.load(fh)
    with MANIFEST.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        cols = list(reader.fieldnames or [])
        rows = to_design_space(list(reader))
    by_key: dict[str, list[dict]] = {}
    for r in rows:
        by_key.setdefault(r["key"], []).append(r)

    def protein(key: str) -> str:
        return seqtools.translate(cds[key]["cds"]).rstrip("*")

    gate = {g["form_id"]: g for g in read(GATE)}

    def attribution(key: str) -> dict:
        fid = key.split(":", 1)[1]
        exp = rejections.get(("isoform", fid))
        if exp and exp["reason"].startswith("explicit rejection"):
            return {
                "date": DATE,
                "decided_by": exp["decided_by"],
                "implemented_by": exp["implemented_by"],
                "reason": exp["reason"],
            }
        why = next(
            (
                f"gate: {g['reason']}"
                for fid, g in gate.items()
                if g["verdict"] == "fail" and key in iso_keys(g)
            ),
            "",
        )
        if why:
            return {
                "date": DATE,
                "decided_by": DECIDED_BY,
                "implemented_by": IMPLEMENTED_BY,
                "reason": why,
            }
        # Neither gated nor rejected: a pick that moved because the seeded replacement draw
        # continues in draw order after the explicit rejections of 2026-11-12.
        return {
            "date": DATE,
            "decided_by": "rule (panel_draw.py replacement within the stratum; Becca Carlson, "
            "2026-11-01), after the 2026-11-12 explicit rejections",
            "implemented_by": IMPLEMENTED_BY,
            "reason": "drawn by rule as the replacement in its stratum, or displaced by it",
        }

    def ortholog_l1(rs: list[dict]) -> list[str]:
        return [
            "\t".join(r.get(c, "") for c in cols)
            for r in rs
            if r["library"] == ld.LIBRARY_L1 and r["set"] == "ortholog" and ld.is_active(r)
        ]

    orthologs_before = ortholog_l1(rows)
    # Proteins of every active member outside the isoform panel (the adopted ortholog panel
    # included), which no pick may repeat.
    fixed_proteins = {
        protein(r["key"]): r["candidate_id"]
        for r in rows
        if r["library"] == ld.LIBRARY_L1
        and ld.is_active(r)
        and r["set"] != "isoform"
        and r["key"] in cds
    }

    rejections: dict[tuple[str, str], dict] = {}
    if REJECTIONS.exists():
        rejections = {
            (r["panel"], r["form_id"]): {
                c: r.get(c) or RULE_ROW.get(c, "") for c in REJECTION_COLUMNS[2:]
            }
            for r in read(REJECTIONS)
        }
    for k, v in explicit_rejections().items():
        rejections.setdefault(k, v)
    # The original ortholog panel's rejections are kept as recorded, never re-derived.

    work = Path(tempfile.mkdtemp(prefix="panel-redraw-"))
    rej_path = work / "panel_rejections.tsv"
    for n in range(1, MAX_PASSES + 1):
        with rej_path.open("w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(REJECTION_COLUMNS)
            w.writerows(
                [p, f, *(v[c] for c in REJECTION_COLUMNS[2:])]
                for (p, f), v in sorted(rejections.items())
            )
        run_samplers(work, rej_path)

        picks: list[tuple[str, str, str]] = []  # (panel, form id, manifest key or "")
        for r in read(work / "isoform_panel.tsv"):
            hit = sorted(iso_keys(r) & set(by_key))
            picks.append(("isoform", iso_id(r), hit[0] if hit else ""))

        new: dict[tuple[str, str], str] = {}
        seen = dict(fixed_proteins)
        # Panel order is (tier/species, band, class, gene): earlier picks keep a shared protein.
        for panel, fid, key in picks:
            if not key or key not in cds:
                new[(panel, fid)] = "no designed CDS (no transcript encodes the drawn form)"
                continue
            p = protein(key)
            if p in seen:
                new[(panel, fid)] = f"designed protein duplicates {seen[p]}"
                continue
            seen[p] = by_key[key][0]["candidate_id"] if by_key[key] else key
        new = {k: v for k, v in new.items() if k not in rejections}
        print(f"pass {n}: {len(picks)} picks, {len(new)} new rejection(s)")
        for (panel, fid), why in sorted(new.items()):
            print(f"  reject {panel} {fid}: {why}")
        if not new:
            break
        rejections.update({k: {**RULE_ROW, "reason": why} for k, why in new.items()})
    else:
        raise SystemExit(f"no fixed point after {MAX_PASSES} passes")

    # ---- reconcile the manifest -------------------------------------------------------
    target = [key for _, _, key in picks]
    target_set = set(target)
    l1_panel = [r for r in rows if r["library"] == ld.LIBRARY_L1 and r["set"] == "isoform"]
    next_n = next_l1_member_number(rows)
    changes, added = [], []
    for r in l1_panel:
        if ld.is_active(r) and r["key"] not in target_set:
            r[ld.RETIRED_COLUMN] = "dropped_from_panel"
            changes.append(
                {
                    "change": "dropped",
                    "candidate_id": r["candidate_id"],
                    "key": r["key"],
                    "detail": "no longer drawn",
                    **attribution(r["key"]),
                }
            )
    active_l1 = {r["key"] for r in l1_panel if ld.is_active(r)}
    for panel, _fid, key in picks:
        if key in active_l1:
            continue
        src = next((r for r in by_key[key] if r["library"] == ld.LIBRARY_L2), None)
        if src is None:
            raise SystemExit(f"{key}: drawn but neither an active L1 nor an L2 member")
        new_id = f"l1-p{next_n:04d}"
        next_n += 1
        src[ld.RETIRED_COLUMN] = f"moved_to:{new_id}"
        added.append(
            {
                **src,
                "candidate_id": new_id,
                **({ld.DESIGN_ID_COLUMN: new_id} if ld.DESIGN_ID_COLUMN in cols else {}),
                "library": ld.LIBRARY_L1,
                "panel": f"{panel}_panel",
                "prior_candidate_id": src["candidate_id"],
                ld.RETIRED_COLUMN: "",
            }
        )
        changes.append(
            {
                "change": "added",
                "candidate_id": new_id,
                "key": key,
                "detail": f"from {src['candidate_id']}; replacement drawn by rule",
                **attribution(key),
            }
        )
    assert ortholog_l1(rows + added) == orthologs_before, "an L1 ortholog row changed"
    print(
        f"\nmanifest: {len(added)} added to L1, "
        f"{sum(c['change'] == 'dropped' for c in changes)} dropped"
    )
    for c in changes:
        print(f"  {c['change']:<8} {c['candidate_id']} {c['key']}  ({c['detail']})")

    if not args.write:
        print("\ndry run -- re-run with --write")
        return 0
    with MANIFEST.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows + added)
    print(RENUMBER_HINT)
    shutil.copy(rej_path, REJECTIONS)
    for name, dest in PANEL_FILES.items():
        # The pinned shortlist code writes CRLF; git stores LF. Copying the CRLF bytes and
        # publishing them put S3 out of step with every fresh checkout (2026-11-03).
        (HERE / dest).write_bytes((work / name).read_bytes().replace(b"\r\n", b"\n"))
    # Earlier changes are kept; a rerun that recomputes a change replaces its row.
    redone = {(c["change"], c["candidate_id"]) for c in changes}
    prior = [
        {**PRIOR_ATTRIBUTION, **c}
        for c in read(CHANGES)
        if (c["change"], c["candidate_id"]) not in redone
    ]
    with CHANGES.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=CHANGE_COLUMNS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(prior + changes)
    print(f"wrote {MANIFEST.name}, {REJECTIONS.name}, {CHANGES.name} and the panel tables")
    if args.publish:
        for name in ("isoform_panel.tsv",):
            processed.upload(str(HERE / PANEL_FILES[name]), DESIGN_PREFIX + name)
            print(f"published processed://{DESIGN_PREFIX}{name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
