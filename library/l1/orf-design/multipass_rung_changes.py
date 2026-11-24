"""What #283's redraw and reference fixes changed on the "N at residue 1, multipass" rung.

Compares the rung and its built arms as recorded on PR #280 (commit `BEFORE`, main's
numbering) with the current tables (after merging #283). Members are matched on `design_id`,
which both numberings share. Writes `multipass_rung_changes.tsv`, one row per design_id that
left, entered, or stayed with a changed construct. Nothing is dropped silently: a member that
left carries the reason from #283's own change records (`panel_changes.tsv`,
`reference_changes.tsv`, `lane_overrides.tsv`) or its new lane in `l1_tag_paths.tsv`.

Attribution: implemented by the Claude Code session l1-nterm-res1-multipass, at the request of
l1-integration (rekey after #283, 2026-11-12). The changes themselves were decided by Becca
("Gate cohort, redraw"; "Send all fixes", 2026-11-12) and applied by l1-reviewer-package on
PR #283.

    uv run --locked python multipass_rung_changes.py
"""

from __future__ import annotations

import csv
import io
from pathlib import Path

HERE = Path(__file__).parent
REL = "library/l1/orf-design/"
#: PR #280's head before #283 was merged: the last state in main's numbering.
BEFORE = "fbd83ab2"
OUT = HERE / "multipass_rung_changes.tsv"
RUNG = "N at residue 1, multipass"
IMPLEMENTED_BY = "Claude Code session l1-nterm-res1-multipass (rekey after PR #283)"
DELEGATED = 'l1-nterm-res1-multipass (delegated by Becca 2026-11-12, "use your best judgment")'


def at(commit: str, name: str) -> list[dict]:
    from deliverome_analysis.provenance import read_pinned_text

    text = read_pinned_text(commit, f"{REL}{name}", repo=HERE.parents[2])
    return list(csv.DictReader(io.StringIO(text), delimiter="\t"))


def now(name: str) -> list[dict]:
    with (HERE / name).open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main() -> int:
    old = {r["design_id"]: r for r in at(BEFORE, "nterm_res1_multipass.tsv")}
    new = {r["design_id"]: r for r in now("nterm_res1_multipass.tsv")}
    old_c = {r["construct_id"]: r for r in at(BEFORE, "constructs_l1_multipass.tsv")}
    new_c = {r["construct_id"]: r for r in now("constructs_l1_multipass.tsv")}
    manifest = {r["design_id"]: r for r in now("candidate_manifest_tranches.tsv")}
    paths = {r["design_id"]: r for r in now("l1_tag_paths.tsv") if r.get("design_id")}
    if not paths:  # l1_tag_paths is keyed on candidate_id; map through the manifest
        by_cid = {r["candidate_id"]: r for r in now("l1_tag_paths.tsv")}
        paths = {
            d: by_cid[m["candidate_id"]] for d, m in manifest.items() if m["candidate_id"] in by_cid
        }
    why: dict[str, dict] = {}
    for r in now("panel_changes.tsv"):
        why.setdefault(r["candidate_id"], r)
    for r in now("reference_changes.tsv"):
        why[r["candidate_id"]] = r
    for r in now("lane_overrides.tsv"):  # keyed on design_id; applied by the writer on #283
        why[r["design_id"]] = r

    rows = []
    for d in sorted(set(old) | set(new)):
        o, n = old.get(d), new.get(d)
        if o and not n:
            src = why.get(d, {})
            lane = paths.get(d, {}).get("placement", "")
            if lane:
                change, reason = (
                    f"moved to {lane}",
                    src.get("reason", "") or "new lane in l1_tag_paths",
                )
            elif manifest.get(d, {}).get("retired_reason"):
                change, reason = "retired", manifest[d]["retired_reason"]
                if "ortholog_cds_policy" in reason:
                    src = {"decided_by": "Becca (Compara 116 adoption)", **src}
            else:
                change, reason = "retired", src.get("reason", "not on the L1 roster")
            decided_by = src.get("decided_by", "")
        elif n and not o:
            src = why.get(d, {})
            change, reason = "entered", src.get("reason", "") or src.get("detail", "")
            decided_by = src.get("decided_by", "")
            if not reason and n["candidate_set"] == "ortholog":
                reason = "Compara 116 ortholog panel adopted on PR #283"
                decided_by = decided_by or "Becca (Compara 116 adoption)"
        else:

            def built(cs: dict, d: str = d) -> str:
                return ",".join(
                    sorted(f"{c.split(':')[1]}/{r['status']}" for c, r in cs.items()
                           if c.startswith(d + ":"))
                )  # fmt: skip

            before, after = built(old_c), built(new_c)
            if before == after and o["candidate_id"] == n["candidate_id"]:
                continue
            reason = f"{o['candidate_id']} -> {n['candidate_id']}"
            decided_by = ""
            if before == after:
                change = "renumbered"
            else:
                change = f"design changed: {before or 'none'} -> {after}"
                if d == "l1-p2477":
                    reason += "; passes the isoform gate, so it is built"
                    decided_by = 'Becca, 2026-11-12 ("Gate cohort, redraw")'
                else:
                    reason += "; one design decided per member under Becca's 2026-11-12 delegation"
                    decided_by = DELEGATED
        rows.append(
            {
                "design_id": d,
                "gene_symbol": (n or o)["gene_symbol"],
                "candidate_set": (n or o)["candidate_set"],
                "species": (n or o)["species"],
                "candidate_id_before": o["candidate_id"] if o else "",
                "candidate_id_now": n["candidate_id"] if n else "",
                "change": change,
                "reason": reason,
                "decided_by": decided_by,
                "implemented_by": IMPLEMENTED_BY,
            }
        )
    cols = list(rows[0])
    with OUT.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    counts: dict[str, int] = {}
    for r in rows:
        k = r["change"].split(":")[0]
        counts[k] = counts.get(k, 0) + 1
    print(f"rung before {len(old)}, now {len(new)}; wrote {OUT.name}: {counts}")
    for r in rows:
        if r["change"] != "renumbered":
            print(
                f"  {r['design_id']} {r['gene_symbol']:10s} {r['candidate_set']:17s} {r['change']} | {r['reason'][:90]}"
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
