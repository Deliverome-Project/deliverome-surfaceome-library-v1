"""Redraw the L1 isoform and ortholog panels by rule, and bring the manifest into line.

The panels are seeded, class-stratified draws (``scratch/becca/appris-principal/
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

from member_id_allocations import next_l1_member_number

from deliverome_analysis import library_design as ld
from deliverome_analysis import seqtools
from deliverome_analysis.s3 import processed

HERE = Path(__file__).parent
REPO = HERE.parents[2]
SAMPLERS = REPO / "scratch/becca/appris-principal"
MANIFEST = HERE / "candidate_manifest_tranches.tsv"
REJECTIONS = HERE / "panel_rejections.tsv"
CHANGES = HERE / "panel_changes.tsv"
PANEL_FILES = (
    "isoform_panel.tsv",
    "ortholog_panel.tsv",
    "isoform_panel_replacements.tsv",
    "ortholog_panel_replacements.tsv",
    "isoform_panel_summary.json",
    "ortholog_panel_summary.json",
)
DESIGN_PREFIX = "sprint4/s4e9-orf-design-fg-panel/library-design/"
MAX_PASSES = 10


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
    for script in ("sample_isoform_panel.py", "sample_ortholog_panel.py"):
        subprocess.run(
            [sys.executable, str(SAMPLERS / script), "--snapshot", "--reject", str(rejections)],
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
        rows = list(reader)
    by_key: dict[str, list[dict]] = {}
    for r in rows:
        by_key.setdefault(r["key"], []).append(r)

    def protein(key: str) -> str:
        return seqtools.translate(cds[key]["cds"]).rstrip("*")

    # Proteins of every active member outside the two panels, which no pick may repeat.
    fixed_proteins = {
        protein(r["key"]): r["candidate_id"]
        for r in rows
        if r["library"] == ld.LIBRARY_L1
        and ld.is_active(r)
        and r["set"] not in ("isoform", "ortholog")
        and r["key"] in cds
    }

    rejections: dict[tuple[str, str], str] = {}
    if REJECTIONS.exists():
        rejections = {(r["panel"], r["form_id"]): r["reason"] for r in read(REJECTIONS)}

    work = Path(tempfile.mkdtemp(prefix="panel-redraw-"))
    rej_path = work / "panel_rejections.tsv"
    for n in range(1, MAX_PASSES + 1):
        with rej_path.open("w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(["panel", "form_id", "reason"])
            w.writerows([p, f, why] for (p, f), why in sorted(rejections.items()))
        run_samplers(work, rej_path)

        picks: list[tuple[str, str, str]] = []  # (panel, form id, manifest key or "")
        for r in read(work / "isoform_panel.tsv"):
            hit = sorted(iso_keys(r) & set(by_key))
            picks.append(("isoform", iso_id(r), hit[0] if hit else ""))
        for r in read(work / "ortholog_panel.tsv"):
            key = f"ort:{r['species']}:{r['human_symbol']}:{r['picked_acc']}"
            picks.append(("ortholog", ort_id(r), key if key in by_key else ""))

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
        rejections.update(new)
    else:
        raise SystemExit(f"no fixed point after {MAX_PASSES} passes")

    # ---- reconcile the manifest -------------------------------------------------------
    target = [key for _, _, key in picks]
    target_set = set(target)
    l1_panel = [
        r for r in rows if r["library"] == ld.LIBRARY_L1 and r["set"] in ("isoform", "ortholog")
    ]
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
            }
        )
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
    shutil.copy(rej_path, REJECTIONS)
    for name in PANEL_FILES:
        # The pinned shortlist code writes CRLF; git stores LF. Copying the CRLF bytes and
        # publishing them put S3 out of step with every fresh checkout (2026-11-03).
        (HERE / name).write_bytes((work / name).read_bytes().replace(b"\r\n", b"\n"))
    with CHANGES.open("w", newline="") as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=["change", "candidate_id", "key", "detail"],
            delimiter="\t",
            lineterminator="\n",
        )
        w.writeheader()
        w.writerows(changes)
    print(f"wrote {MANIFEST.name}, {REJECTIONS.name}, {CHANGES.name} and the panel tables")
    if args.publish:
        for name in ("isoform_panel.tsv", "ortholog_panel.tsv"):
            processed.upload(str(HERE / name), DESIGN_PREFIX + name)
            print(f"published processed://{DESIGN_PREFIX}{name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
