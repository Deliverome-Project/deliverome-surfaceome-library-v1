"""Verdicts from the 2026-11-10 SignalP 6.0 run: the five new members and the Igκ swaps.

Reads ``signalp_inputs.tsv`` (what each sequence is) and the run's ``calls.json`` (pinned on
S3 by ``signalp6_modal.py``; run ID below) and applies the sweep's verdict rule to every
tagged design against its own untagged parent: native for the new members, the Igκ-swapped
protein for the Igκ arms.

    uv run --locked python analyse_signalp_run.py

**Attribution** (convention of 2026-11-10): analysis only, no decision taken. Written by the
Claude Code session ``nterm-reviewer-package`` (PR #254).
"""

from __future__ import annotations

import collections
import csv
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import nterm_placement as npl  # noqa: E402

RUN_ID = "nterm_sp6_2026_11_10"
CALLS_KEY = (
    f"sprint4/s4e9-orf-design-fg-panel/reviewer-package/nterm/signalp6-runs/{RUN_ID}/calls.json"
)
PARENT = {"tagged": "native", "igk_tagged": "igk_native"}


def main() -> int:
    from deliverome_analysis.s3 import processed

    with (HERE / "signalp_inputs.tsv").open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    calls = npl.load_json(processed, CALLS_KEY)
    parents = {
        (r["candidate_id"], r["arm"]): calls[r["sequence_sha256"]]
        for r in rows
        if r["arm"] in PARENT.values()
    }
    out = []
    for r in rows:
        call = calls[r["sequence_sha256"]]
        row = {
            k: r[k]
            for k in (
                "candidate_id",
                "design_id",
                "gene_symbol",
                "set",
                "arm",
                "after_residue",
                "native_after_residue",
            )
        }
        row |= {
            "prediction": call["prediction"],
            "cleavage_site": call["cleavage_site"] or "",
            "cleavage_probability": call["cleavage_probability"] or "",
            "sp_probability": call["sp_probability"],
        }
        if r["arm"] in PARENT:
            parent = parents[(r["candidate_id"], PARENT[r["arm"]])]
            row["verdict"] = (
                "parent_has_no_signal_peptide"
                if parent["prediction"] != "SP"
                else npl.verdict(parent, call, int(r["after_residue"]))
            )
        else:
            row["verdict"] = ""
        out.append(row)
    with (HERE / "signalp_run_verdicts.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out)

    summary: dict = {"run_id": RUN_ID, "calls_key": CALLS_KEY}
    for s, parent_arm, tag_arm in (
        ("new_member", "native", "tagged"),
        ("igk_swap", "igk_native", "igk_tagged"),
    ):
        par = [r for r in out if r["set"] == s and r["arm"] == parent_arm]
        tag = [r for r in out if r["set"] == s and r["arm"] == tag_arm]
        per = collections.defaultdict(collections.Counter)
        for r in tag:
            per[r["candidate_id"]][r["verdict"]] += 1
        summary[s] = {
            "members": len(par),
            "parent_calls": dict(collections.Counter(r["prediction"] for r in par)),
            "parent_cleavage_sites": dict(
                sorted(collections.Counter(str(r["cleavage_site"]) for r in par).items())
            ),
            "tagged_designs": len(tag),
            "tagged_verdicts": dict(sorted(collections.Counter(r["verdict"] for r in tag).items())),
            "members_with_a_preserved_site": sum(1 for c in per.values() if c["preserved"]),
            "members_with_none_preserved": sorted(
                f"{cid} {next(r['gene_symbol'] for r in par if r['candidate_id'] == cid)}"
                for cid, c in per.items()
                if not c["preserved"]
            ),
        }
    (HERE / "signalp_run_summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(
        json.dumps(
            {
                k: {kk: vv for kk, vv in v.items() if kk != "members_with_none_preserved"}
                if isinstance(v, dict)
                else v
                for k, v in summary.items()
            },
            indent=1,
        )
    )
    print({s: len(summary[s]["members_with_none_preserved"]) for s in ("new_member", "igk_swap")})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
