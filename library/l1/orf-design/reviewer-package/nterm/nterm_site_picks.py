"""One insertion site per L1 "N after SP" member, by the re-prediction-gated rule, with status.

Becca chose the rule on 2026-11-11 ("Re-prediction-gated"): among a member's candidate sites,
take one whose tagged design keeps the native cleavage site (SignalP 6.0); prefer sites strictly
before the predicted fold onset, then a cleavage-site probability that falls by no more than
0.07, then the highest tagged cleavage-site probability, then the smallest offset
(``nterm_build_options.gated_pick``). The cassette is the 51-aa p003 tandem (Becca, 2026-11-10).

Each member gets one row and a ``status``:

| status | meaning |
|---|---|
| ``ready`` | a gated site before the fold, from the sweep or the 2026-11-10 run |
| ``held_inside_fold`` | the only preserved sites are inside the predicted fold (decision 3 in BUILD_OPTIONS.md) |
| ``held_moved_upstream`` | no preserved site; cleavage moves but stays upstream of the cassette (not yet accepted) |
| ``held_igk`` | native signal fails at every site; the Igκ swap is recommended, not yet decided (decision 4) |
| ``refused`` | no interval: no SignalP row, or the two signal-peptide calls disagree (decision 5) |

For a picked site it records the residues either side, the CDS nucleotide the cassette follows
(``3n``), and the SHA-256 of the tagged protein that SignalP scored, so a builder can check it
inserts what was predicted. No insert is assembled here: N-after-SP members have no issued
cassette encoding yet (see BUILD_OPTIONS.md).

**Attribution** (convention of 2026-11-10): the rule was decided by Becca on 2026-11-11; this
table is implemented by the Claude Code session ``nterm-reviewer-package`` (PR #254). Every row
carries ``decided_by``, ``implemented_by`` and ``reason``.

    uv run --locked python nterm_site_picks.py
"""

from __future__ import annotations

import collections
import csv
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import nterm_build_options as nbo  # noqa: E402
import nterm_placement as npl  # noqa: E402

OUT = HERE / "nterm_site_picks.tsv"
SUMMARY = HERE / "nterm_site_picks_summary.json"
DECIDED_BY = "Becca, 2026-11-11: re-prediction-gated site rule; 2026-11-10: p003 tandem"
IMPLEMENTED_BY = "Claude Code session nterm-reviewer-package (PR #254)"
STATUS_OF = {
    "preserved_before_fold": "ready",
    "preserved_inside_fold_only": "held_inside_fold",
    "cleavage_moved_upstream_only": "held_moved_upstream",
    "signal_peptide_lost_everywhere": "held_igk",
    "cleavage_crosses_junction_everywhere": "held_igk",
}
REASON = {
    "ready": "gated site before the predicted fold; SignalP keeps the native cleavage site",
    "held_inside_fold": "every preserved site lies inside the predicted fold; awaiting decision 3",
    "held_moved_upstream": "no site keeps the cleavage site exactly; awaiting acceptance of an upstream shift",
    "held_igk": "native signal peptide fails at every site; Igκ swap recommended, awaiting decision 4",
    "refused": "no admissible interval; awaiting decision 5",
}


def moved_pick(positions: list[dict], cleavage: int) -> dict | None:
    """A proposed site for a member whose cleavage only ever moves upstream of the cassette.

    Proposed by the Claude Code session nterm-reviewer-package (PR #254, 2026-11-12), not
    decided: among candidates whose tagged cleavage moves but stays before the insertion,
    prefer sites before the predicted fold, then the smallest shift of the cleavage site, then
    the highest tagged cleavage-site probability, then the smallest offset.
    """
    moved = [p for p in positions if p["verdict"] == "cleavage_moved" and p["tagged_cleavage"]]
    if not moved:
        return None
    return min(
        moved,
        key=lambda p: (
            p["within_strict_interval"] != "1",
            abs(int(p["tagged_cleavage"]) - cleavage),
            -float(p["tagged_cleavage_probability"]),
            int(p["after_residue"]),
        ),
    )


def _tsv(path: Path) -> list[dict]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main() -> int:
    from deliverome_analysis.s3 import processed

    members, _ = npl.load_members(processed)
    opts = {r["candidate_id"]: r for r in _tsv(HERE / "nterm_build_options.tsv")}
    positions = collections.defaultdict(list)
    for p in _tsv(HERE / "nterm_positions.tsv"):
        positions[p["candidate_id"]].append(p)
    # Members added after the sweep are re-predicted by nterm_overnight.py (group new_member),
    # in the pinned runs it names; earlier ones (2026-11-10) were retired by PR #283.
    late = collections.defaultdict(dict)
    late_native = {}
    for r in _tsv(HERE / "overnight_verdicts.tsv"):
        if r["group"] == "new_member":
            if r["arm"] == "native":
                late_native[r["candidate_id"]] = float(r["cp"] or 0)
            else:
                late[r["candidate_id"]][r["after_residue"]] = r | {"cleavage_probability": r["cp"]}

    rows = []
    for cid in sorted(members):
        m, o = members[cid], opts[cid]
        sub = o["subgroup"]
        pos = positions.get(cid, [])
        if sub == "not_repredicted":
            for p in pos:
                v = late[cid][p["after_residue"]]
                # tagged_cleavage too, so moved_pick can place a member whose cleavage only
                # moves upstream (CD300LB, 2026-11-21; no earlier new member needed it)
                p.update(
                    verdict=v["verdict"],
                    tagged_cleavage_probability=v["cleavage_probability"],
                    tagged_cleavage=v["cleavage"],
                )
            native_cp = late_native[cid]
            verdicts = {p["verdict"] for p in pos}
            if "preserved" in verdicts:
                strict = any(
                    p["verdict"] == "preserved" and p["within_strict_interval"] == "1" for p in pos
                )
                sub = "preserved_before_fold" if strict else "preserved_inside_fold_only"
            elif verdicts == {"signal_peptide_lost"}:
                sub = "signal_peptide_lost_everywhere"
            elif "cleavage_moved" in verdicts:
                sub = "cleavage_moved_upstream_only"
            else:
                sub = "cleavage_crosses_junction_everywhere"
        else:
            native_cp = float(o["native_cleavage_probability"] or 0)
        status = "refused" if sub.startswith("refused") else STATUS_OF[sub]
        pick = nbo.gated_pick([p for p in pos if p["verdict"]], native_cp) if pos else None
        if pick is None and sub == "cleavage_moved_upstream_only":
            pick = moved_pick(pos, int(o["cleavage_site"]))
        row = {
            "candidate_id": cid,
            "design_id": m["design_id"],
            "gene_symbol": m["gene_symbol"],
            "candidate_set": m["candidate_set"],
            "subgroup": sub,
            "status": status,
            "cleavage_site": o["cleavage_site"],
            "insert_after_residue": pick["after_residue"] if pick else "",
            "offset_past_cleavage": pick["offset_past_cleavage"] if pick else "",
            "residue_before": pick["residue_before"] if pick else "",
            "residue_after": pick["residue_after"] if pick else "",
            "cds_insert_after_nt": pick["cds_insert_nt"] if pick else "",
            "inside_predicted_fold": "" if not pick else int(pick["within_strict_interval"] != "1"),
            "native_cleavage_probability": native_cp or "",
            "tagged_cleavage_probability": pick["tagged_cleavage_probability"] if pick else "",
            "softened_over_0_07": ""
            if not pick
            else int(native_cp - float(pick["tagged_cleavage_probability"]) > nbo.SOFTENING_BOUND),
            "tagged_protein_sha256": pick["tagged_sha256"] if pick else "",
            "tagged_cleavage_site": pick.get("tagged_cleavage", "") if pick else "",
            "pick_rule": (
                ""
                if not pick
                else "moved_upstream (proposed)"
                if sub == "cleavage_moved_upstream_only"
                else "reprediction_gated"
            ),
            "decided_by": (
                "proposed, not decided: nterm-reviewer-package (PR #254), 2026-11-12"
                if sub == "cleavage_moved_upstream_only"
                else DECIDED_BY
            ),
            "implemented_by": IMPLEMENTED_BY,
            "reason": REASON[status],
        }
        rows.append(row)

    with OUT.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    ready = [r for r in rows if r["status"] == "ready"]
    summary = {
        "tranche": "l1",
        "members": len(rows),
        "status": dict(sorted(collections.Counter(r["status"] for r in rows).items())),
        "ready_offset_histogram": dict(
            sorted(collections.Counter(int(r["offset_past_cleavage"]) for r in ready).items())
        ),
        "ready_softened_over_0_07": sum(r["softened_over_0_07"] == 1 for r in ready),
        "ready_by_set": dict(
            sorted(collections.Counter(r["candidate_set"] for r in ready).items())
        ),
    }
    SUMMARY.write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps(summary, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
