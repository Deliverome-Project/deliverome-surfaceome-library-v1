"""Rerun the isoform gate on the cohort of record once every form has a DeepTMHMM v1 row.

The gate of record (``isoform_gate.tsv``, 2026-11-12) passed 802 of the 1,894 cohort forms and
failed 931 only because the form had no v1 prediction. ``isoform_v1_topology.py`` predicted
those 955 proteins on 2026-11-20 and published them as ``(topo_2026_11_20_isoform_cohort,
human_isoforms)`` in a new snapshot. This reruns ``gate_isoform_cohort.gate`` unchanged against
that snapshot and writes ``isoform_gate_full_v1.tsv``; the gate of record is not touched.

It also checks the new v1 calls against the October DeepTMHMM v1 run of the cohort shortlist
(``scratch/becca/appris-principal/topology_shortlist.3line``), which covered many of the same
proteins but was never loaded into the topology snapshot: an identical sequence should get an
identical per-residue string.

    uv run --locked python regate_isoform_forms.py

Decision: Becca Carlson, 2026-11-20 (review comment 27: "well then we need to address this
first"). Implemented by the Claude Code session l1-isoform-universe (subagent of the brief
review).
"""

from __future__ import annotations

import csv
import json
import tempfile
from collections import Counter
from pathlib import Path

import gate_isoform_cohort as g

HERE = Path(__file__).resolve().parent
RUN = HERE / "topology_runs" / "2026-11-20_isoform_cohort"
OUT = HERE / "isoform_gate_full_v1.tsv"
SUMMARY = HERE / "isoform_gate_full_v1_summary.json"
OCTOBER = HERE.parents[2] / "scratch/becca/appris-principal/topology_shortlist.3line"
DECIDED_BY = "Becca Carlson, 2026-11-20 (review comment 27: 'address this first')"
IMPLEMENTED_BY = "Claude Code session l1-isoform-universe (subagent of the brief review)"


def three_line(path: Path) -> dict[str, str]:
    lines, out = path.read_text().splitlines(), {}
    for i in range(len(lines) - 2):
        if lines[i].startswith(">"):
            out[lines[i + 1].strip()] = lines[i + 2].strip()
    return out


def main() -> int:
    pub = json.loads(RUN.with_suffix(".published.json").read_text())
    inp = g.load_inputs(
        Path(tempfile.mkdtemp(prefix="isoform-regate-")),
        snapshot_key=pub["snapshot_key"],
        snapshot_sha256=pub["snapshot_sha256"],
    )
    rows = g.gate(inp)
    for r in rows:
        r |= {
            "decided_by": DECIDED_BY,
            "implemented_by": IMPLEMENTED_BY,
            "decision_reason": g.DECISION_REASON,
        }
    OUT.write_text(g.tsv_text(rows))
    with (HERE / "isoform_gate.tsv").open() as fh:
        before = {r["form_id"]: r for r in csv.DictReader(fh, delimiter="\t")}
    moves = Counter((before[r["form_id"]]["verdict"], r["verdict"]) for r in rows)
    newly = [
        r
        for r in rows
        if before[r["form_id"]]["verdict"] == "fail"
        and "unverifiable" in before[r["form_id"]]["reason"]
    ]
    new_fail_first = Counter(r["reason"].split("; ")[0] for r in newly if r["verdict"] == "fail")
    # October run vs this run, sequence by sequence.
    now, octo = three_line(RUN.with_suffix(".3line")), three_line(OCTOBER)
    shared = set(now) & set(octo)
    summary = g.summarise(rows) | {
        "decision": DECIDED_BY,
        "implemented_by": IMPLEMENTED_BY,
        "topology_snapshot": {"key": pub["snapshot_key"], "sha256": pub["snapshot_sha256"]},
        "gate_of_record": {
            "pass": sum(r["verdict"] == "pass" for r in before.values()),
            "fail": sum(r["verdict"] == "fail" for r in before.values()),
        },
        "verdict_moves": {f"{a} -> {b}": n for (a, b), n in sorted(moves.items())},
        "previously_unverifiable": {
            "forms": len(newly),
            "now_pass": sum(r["verdict"] == "pass" for r in newly),
            "now_fail": sum(r["verdict"] == "fail" for r in newly),
            "now_fail_first_reason": dict(new_fail_first.most_common()),
        },
        "pass_by_tier": {
            t: sum(r["verdict"] == "pass" for r in rows if r["cohort_tier"] == t)
            for t in ("observed", "full_only")
        },
        "october_shortlist_run_check": {
            "sequences_in_both_runs": len(shared),
            "identical_per_residue_string": sum(now[s] == octo[s] for s in shared),
        },
    }
    SUMMARY.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
