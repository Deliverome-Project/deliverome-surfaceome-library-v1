"""Refill the isoform-panel slots emptied by the MANE-GTF reference switch (2026-11-22).

``apply_mane_gtf.py`` moved six library references onto proteins that were already isoform-panel
members, and retired those members ``duplicates_reference:...`` (ADGRB2, LSR, MUC13, PCDH9,
STEAP3, TMEM134). Each slot is refilled by the panel's own replacement rule
(``panel_draw.py``): continue the **same seeded RNG** over the slot's stratum, ``rng.choice``
over the candidates not picked, not rejected and not invalid.

**Why the refill continues the stream after the whole draw** instead of re-running the draw with
six more rejections. Rerunning makes the RNG call ``choice`` once per extra rejection inside
the strata where the six sit. Every replacement drawn later in the stream then shifts. Run on
2026-11-22, that moved 7 members who were never retired (ATP1A1, CD55, GPM6B, MOSMO, NTM,
SLC43A2, TMPRSS13). The draw is replayed exactly as applied (its RNG state verified by
reproducing the 150 picks of ``isoform_panel_redraw_2026-11-22.tsv``). The six refills are
then drawn from the RNG's state at the end of that draw, in draw order. No other member moves.

A candidate is skipped (and the next one drawn) when no designed CDS in the current record
encodes it, or its designed protein repeats an active L1 member's protein: any reference
(now the MANE-GTF ones) or another panel member. Each refill moves from its L2 row
(retired ``moved_to:<new id>``) to a new L1 row, id from ``next_l1_member_number`` on the
manifest as read. ``panel_refill_2026-11-22_changes.tsv`` records each refill and each skip.
The manifest is read once, re-read and compared just before the single write (row count and
bytes), and written once.

**Decided by** Becca Carlson, 2026-11-22 (the redraw); the refill was requested by the
coordinating session. Implemented by Claude Code session l1-review-export-scope (worker fork).

    uv run --locked python refill_isoform_panel_2026_11_22.py           # dry run
    uv run --locked python refill_isoform_panel_2026_11_22.py --write   # write the manifest
    uv run --locked python refill_isoform_panel_2026_11_22.py --check   # fail unless applied
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import random
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import apply_panel_redraw_2026_11_22 as apply_mod  # noqa: E402
import redraw_isoform_panel_2026_11_22 as rd  # noqa: E402

MANIFEST = HERE / "candidate_manifest_tranches.tsv"
REDRAW = HERE / "isoform_panel_redraw_2026-11-22.tsv"
CHANGES = HERE / "panel_refill_2026-11-22_changes.tsv"
EMPTIED = "duplicates_reference:"
DECIDED_BY = "Becca Carlson, 2026-11-22 (panel redraw; refill of slots emptied by MANE-GTF)"
IMPLEMENTED_BY = "Claude Code session l1-review-export-scope (worker fork), 2026-11-22"
COLUMNS = [
    "change",
    "candidate_id",
    "key",
    "prior_candidate_id",
    "form_id",
    "stratum",
    "refills",
    "detail",
    "v1_topology",
    "decided_by",
    "implemented_by",
]


def form_id(r: dict) -> str:
    return f"{r['gene']}:{r['uniprot'].strip() or r['ensembl_transcripts']}"


def replay(work: Path, extra: dict[str, str]) -> tuple[random.Random, list[dict], list[dict]]:
    """The applied draw, replayed. Returns (its RNG at the end, each stratum call, the panel)."""
    import sample_isoform_panel as s

    rngs: list[random.Random] = []
    calls: list[dict] = []
    base_random, base_replace = s.random, s.replace_rejected

    class Recording(random.Random):
        def __init__(self, seed):
            super().__init__(seed)
            rngs.append(self)

    def recording_replace(rng, picked, stratum, fid, rejected, invalid, used):
        out, log = base_replace(rng, picked, stratum, fid, rejected, invalid, used)
        calls.append(
            {"out": out, "stratum": stratum, "rejected": rejected, "invalid": invalid, "used": used}
        )
        return out, log

    s.random = type("R", (), {"Random": Recording})
    s.replace_rejected = recording_replace
    try:
        panel, _, _ = rd.draw(work, split=True, extra=extra)
    finally:
        s.random, s.replace_rejected = base_random, base_replace
    assert len(rngs) == 1, "the sampler made more than one RNG"
    return rngs[0], calls, panel


def plan(work: Path, rows: list[dict]) -> tuple[list[dict], list[dict]]:
    """(refill decisions, change-log rows) for the slots emptied on ``rows``."""
    from redraw_panels import iso_keys

    from deliverome_analysis import library_design as ld
    from deliverome_analysis import seqtools, topology
    from deliverome_analysis.s3 import processed

    summary = json.loads(rd.OUT_SUMMARY.read_text())
    extra = summary["rejected_as_unbuildable_in_redraw"]
    rng, calls, panel = replay(work, extra)
    applied = {r["form_id"] for r in rd.read(REDRAW) if r["status"] != "dropped"}
    assert {form_id(r) for r in panel} == applied, "the replay is not the applied draw"

    with open(processed.download(ld.CDS_KEY)) as fh:
        cds = json.load(fh)
    protein = lambda k: seqtools.translate(cds[k]["cds"]).rstrip("*")  # noqa: E731
    by_key: dict[str, list[dict]] = {}
    for m in rows:
        by_key.setdefault(m["key"], []).append(m)
    iso_l1 = [m for m in rows if m["set"] == "isoform" and m["library"] == "l1"]
    emptied = {
        m["key"]: m
        for m in iso_l1
        if m["retired_reason"].startswith(EMPTIED) and m["panel"] == "isoform_panel"
    }
    # Slots already refilled by an earlier run of this script, keyed by the emptied member.
    done = {}
    if CHANGES.exists():
        done = {c["refills"]: c for c in rd.read(CHANGES) if c["change"] == "added"}
    taken = {
        protein(m["key"]): m["design_id"]
        for m in rows
        if m["library"] == "l1"
        and not m["retired_reason"]
        and m["key"] in cds
        and m["design_id"] not in {c["candidate_id"] for c in done.values()}
    }

    decisions, log = [], []
    for call in calls:
        for p in call["out"]:
            hit = sorted(iso_keys(p) & set(emptied))
            if not hit:
                continue
            slot = emptied[hit[0]]
            stratum = f"{rd.TIER_NAME[p['cohort_tier']]} / {p['topology_class']}"
            while True:
                cands = [
                    c
                    for c in call["stratum"]
                    if form_id(c) not in call["used"]
                    and form_id(c) not in call["rejected"]
                    and not call["invalid"](c)
                ]
                if not cands:
                    raise SystemExit(f"{form_id(p)}: stratum {stratum} has no candidate left")
                new = rng.choice(cands)
                call["used"].add(form_id(new))
                keys = sorted(iso_keys(new) & set(by_key))
                why = ""
                if not keys or keys[0] not in cds:
                    why = "no designed CDS (no transcript encodes the drawn form)"
                elif protein(keys[0]) in taken:
                    why = f"designed protein duplicates {taken[protein(keys[0])]}"
                if why:
                    log.append(
                        {
                            "change": "skipped",
                            "form_id": form_id(new),
                            "key": keys[0] if keys else "",
                            "stratum": stratum,
                            "refills": slot["design_id"],
                            "detail": why,
                        }
                    )
                    continue
                taken[protein(keys[0])] = keys[0]
                decisions.append(
                    {
                        "slot": slot,
                        "form": new,
                        "key": keys[0],
                        "stratum": stratum,
                        "sequence": protein(keys[0]),
                        "candidates": len(cands),
                    }
                )
                break
    assert len(decisions) == len(emptied), (
        f"{len(emptied)} emptied slots, {len(decisions)} refilled"
    )
    topo = topology.resolve({d["key"]: d["sequence"] for d in decisions}, tool=topology.TOOL_V1)
    found = {r["key"]: r for r in topo.to_dict("records")}
    for d in decisions:
        t = found[d["key"]]
        d["v1"] = f"yes ({t['topology_version']}, {t['cohort']})" if t["found"] else "no"
    return decisions, log


def apply(rows: list[dict], decisions: list[dict], log: list[dict]) -> tuple[list[dict], int]:
    next_n = [apply_mod.next_l1_member_number(rows)]
    out, edits = [], 0
    for d in decisions:
        srcs = [r for r in rows if r["library"] == "l2" and r["key"] == d["key"]]
        assert len(srcs) == 1, f"{d['key']}: {len(srcs)} L2 rows"
        row, changed = apply_mod.promote(rows, srcs[0], "isoform_panel", next_n)
        edits += changed
        out.append(
            {
                "change": "added",
                "candidate_id": row["candidate_id"],
                "key": row["key"],
                "prior_candidate_id": row["prior_candidate_id"],
                "form_id": form_id(d["form"]),
                "stratum": d["stratum"],
                "refills": d["slot"]["design_id"],
                "detail": (
                    f"refills {d['slot']['design_id']} ({d['slot']['key']}, retired "
                    f"{d['slot']['retired_reason'].split(' (')[0]}): next form in the same "
                    f"stratum by the panel's seeded RNG, continued after the applied draw "
                    f"({d['candidates']} candidates); moved from L2 {srcs[0]['candidate_id']}"
                ),
                "v1_topology": d["v1"],
            }
        )
    full = [{**c, "decided_by": DECIDED_BY, "implemented_by": IMPLEMENTED_BY} for c in out + log]
    return full, edits


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
    ap.add_argument("--work", type=Path, default=None)
    a = ap.parse_args(argv)
    raw = MANIFEST.read_bytes()
    reader = csv.DictReader(io.StringIO(raw.decode()), delimiter="\t")
    cols, rows = list(reader.fieldnames or []), list(reader)
    n_before = len(rows)
    decisions, log = plan(a.work or Path(tempfile.mkdtemp(prefix="isoform-refill-")), rows)
    changes, edits = apply(rows, decisions, log)
    for c in changes:
        print(f"{c['change']:8s} {c['candidate_id'] or '-':9s} {c['form_id']:40s} {c['detail']}")
        if c.get("v1_topology"):
            print(f"         v1 topology: {c['v1_topology']}")
    print(f"{edits} manifest edits needed; rows {n_before} -> {len(rows)}")
    if a.check:
        same = edits == 0 and CHANGES.exists() and CHANGES.read_text() == text(changes, COLUMNS)
        print("applied and matches the record" if same else "NOT applied, or the log differs")
        return 0 if same else 1
    if not a.write:
        print("dry run -- re-run with --write")
        return 0
    now = MANIFEST.read_bytes()
    n_now = sum(1 for _ in csv.DictReader(io.StringIO(now.decode()), delimiter="\t"))
    if now != raw or n_now != n_before:
        raise SystemExit(f"manifest changed while planning ({n_before} -> {n_now} rows); rerun")
    assert len(rows) == n_before + edits, "row count does not match the planned additions"
    MANIFEST.write_text(text(rows, cols))
    CHANGES.write_text(text(changes, COLUMNS))
    print(f"wrote {MANIFEST.name} ({len(rows)} rows), {CHANGES.name}")
    print(apply_mod.RENUMBER_HINT)
    return 0


if __name__ == "__main__":
    sys.exit(main())
