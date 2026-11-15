"""Where the N-terminal cassette would go for each L1 "N after SP" member, under three rules.

Reads the committed outputs of ``nterm_placement.py`` (``nterm_members.tsv``,
``nterm_positions.tsv``) and the pinned re-prediction sweep for the native calls. Runs no
predictor and changes no design record. It answers two questions for the build decision:

1. **Which rule, per subgroup.** For every member with re-predicted candidates it applies

   * ``cleavage_plus_1``: insert after the first mature residue (``c + 1``), the default
     ``build_constructs.py`` used before the panel was reworked;
   * ``disorder_anchored``: the most disordered residue in ``c + 1 .. c + 20``
     (``tag_placement.post_sp_junction``), the rule ``rebuild_disorder_anchored.py`` applied
     to the panel;
   * ``reprediction_gated``: among candidates whose tagged SignalP 6.0 call keeps the native
     cleavage site, prefer those strictly before the fold onset, then those whose
     cleavage-site probability falls by no more than 0.07 (the record's gradeable bound),
     then the highest tagged cleavage-site probability, then the smallest offset;

   and records each rule's verdict, so the rules can be compared on what SignalP says.
2. **Whether weak or uncertain native signal peptides predict failure**, the evidence for
   replacing them with the Igκ leader (``build_sp_swap_arms.py``).

**Attribution** (convention of 2026-11-10). No rule here has been adopted: the three rules,
the 0.07 softening bound (the record's calibration limit, from the PR #191 handoff) and the
split that treats an upstream-only cleavage shift as usable (following the record's judgement
of HEG1) are proposals from the Claude Code session ``nterm-reviewer-package`` (PR #254,
2026-11-10), awaiting Becca. Each row's ``rule_status`` says so. The overrides that
``nterm_placement.py`` applied to a member (``override``, ``decided_by``, ``implemented_by``,
``override_reason``) are carried through unchanged.

Writes ``nterm_build_options.tsv`` (one row per member) and ``nterm_build_options.json``.

    uv run --locked python nterm_build_options.py
"""

from __future__ import annotations

import collections
import csv
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
S4E9 = HERE.parents[1]
sys.path.insert(0, str(HERE))

import nterm_placement as npl  # noqa: E402

SOFTENING_BOUND = 0.07
PRESERVED = "preserved"
CONF_BINS = [(0.0, 0.5), (0.5, 0.8), (0.8, 0.9), (0.9, 0.95), (0.95, 1.0001)]


def _tsv(path: Path) -> list[dict]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _bin(x: float) -> str:
    for lo, hi in CONF_BINS:
        if lo <= x < hi:
            return f"{lo:.2f}-{min(hi, 1):.2f}"
    raise ValueError(x)


def pick_key(before_fold: bool, native_cp: float, tagged_cp: float, offset: int) -> tuple:
    """Sort key of the re-prediction-gated rule (Becca, 2026-11-11); the smallest key wins.

    Among sites that keep the native cleavage site: a site before the predicted fold, then one
    whose cleavage-site probability falls by no more than ``SOFTENING_BOUND``, then the
    highest tagged probability, then the smallest offset. The one definition every caller
    uses (``gated_pick``, ``nterm_overnight._pick``, ``nterm_decide`` round 2).
    """
    return (not before_fold, native_cp - tagged_cp > SOFTENING_BOUND, -tagged_cp, offset)


def gated_pick(positions: list[dict], native_cp: float) -> dict | None:
    """The re-prediction-gated choice among preserved candidates, or None."""
    ok = [p for p in positions if p["verdict"] == PRESERVED]
    if not ok:
        return None
    return min(
        ok,
        key=lambda p: pick_key(
            p["within_strict_interval"] == "1",
            native_cp,
            float(p["tagged_cleavage_probability"]),
            int(p["after_residue"]),
        ),
    )


def subgroup(m: dict, positions: list[dict]) -> str:
    if m["class"].startswith("refused"):
        return m["class"]
    if not positions or not any(p["repredicted"] == "1" for p in positions):
        return "not_repredicted"
    verdicts = {p["verdict"] for p in positions}
    if PRESERVED in verdicts:
        strict_ok = any(
            p["verdict"] == PRESERVED and p["within_strict_interval"] == "1" for p in positions
        )
        return "preserved_before_fold" if strict_ok else "preserved_inside_fold_only"
    if verdicts == {"signal_peptide_lost"}:
        return "signal_peptide_lost_everywhere"
    if "cleavage_moved" in verdicts:
        # The site shifts but stays upstream of the cassette, so the cassette is still
        # displayed (the record judges HEG1's +1 shift benign on this ground).
        return "cleavage_moved_upstream_only"
    return "cleavage_crosses_junction_everywhere"


def main() -> int:
    from deliverome_analysis.s3 import processed

    members = _tsv(HERE / "nterm_members.tsv")
    by_member = collections.defaultdict(list)
    for p in _tsv(HERE / "nterm_positions.tsv"):
        by_member[p["candidate_id"]].append(p)
    sweep = npl.load_json(processed, npl.SWEEP_REPREDICTION_KEY)
    # The sweep's member_id is the design id (l1-pNNNN), from before the renumbering.
    native = {r["member_id"]: r for r in sweep if r["role"] == "native"}

    rows = []
    for m in members:
        cid = m["candidate_id"]
        pos = sorted(by_member.get(cid, []), key=lambda p: int(p["after_residue"]))
        nat = native.get(m["design_id"]) or native.get(m["swept_under_duplicate"].split(",")[0])
        sg = subgroup(m, pos)
        row = {
            "candidate_id": cid,
            "design_id": m["design_id"],
            "gene_symbol": m["gene_symbol"],
            "candidate_set": m["candidate_set"],
            "class": m["class"],
            "subgroup": sg,
            "cleavage_site": m["lo"],
            "native_sp_probability": nat["sp_probability"] if nat else "",
            "native_cleavage_probability": nat["cleavage_probability"] if nat else "",
            "n_candidates": m["n_candidates"],
            "n_preserved": m["n_preserved"],
            "n_moved": m["n_moved"],
            "n_crosses_junction": m["n_crosses_junction"],
            "n_sp_lost": m["n_sp_lost"],
            "rule_status": "proposal, no decision taken (nterm-reviewer-package, PR #254)",
            "override": m["override"],
            "decided_by": m["decided_by"],
            "implemented_by": m["implemented_by"],
            "override_reason": m["override_reason"],
        }
        rep = [p for p in pos if p["repredicted"] == "1"]
        if rep and nat:
            at = {int(p["after_residue"]): p for p in rep}
            c1 = rep[0]
            dpick = int(m["disorder_pick"]) if m["disorder_pick"] else None
            g = gated_pick(rep, nat["cleavage_probability"])
            row |= {
                "cleavage_plus_1_site": c1["after_residue"],
                "cleavage_plus_1_verdict": c1["verdict"],
                "disorder_anchored_site": dpick or "",
                "disorder_anchored_verdict": (
                    at[dpick]["verdict"] if dpick in at else "outside_interval"
                ),
                "gated_site": g["after_residue"] if g else "",
                "gated_offset": (int(g["after_residue"]) - int(m["lo"])) if g else "",
                "gated_inside_fold": "" if not g else int(g["within_strict_interval"] != "1"),
                "gated_cleavage_probability": g["tagged_cleavage_probability"] if g else "",
                "gated_softened": (
                    ""
                    if not g
                    else int(
                        nat["cleavage_probability"] - float(g["tagged_cleavage_probability"])
                        > SOFTENING_BOUND
                    )
                ),
            }
        rows.append(row)

    cols = sorted(
        {k for r in rows for k in r}, key=lambda k: list(rows[0]).index(k) if k in rows[0] else 99
    )
    cols = list(rows[0]) + [c for c in cols if c not in rows[0]]
    with (HERE / "nterm_build_options.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n", restval="")
        w.writeheader()
        w.writerows(rows)

    rep_rows = [r for r in rows if r.get("cleavage_plus_1_verdict")]
    rule_verdicts = {
        rule: dict(sorted(collections.Counter(r[f"{rule}_verdict"] for r in rep_rows).items()))
        for rule in ("cleavage_plus_1", "disorder_anchored")
    }
    gated = [r for r in rep_rows if r["gated_site"] != ""]
    conf = collections.defaultdict(collections.Counter)
    for r in rep_rows:
        b = _bin(r["native_cleavage_probability"])
        conf[b]["members"] += 1
        conf[b]["cleavage_plus_1_not_preserved"] += r["cleavage_plus_1_verdict"] != PRESERVED
        conf[b]["no_candidate_preserved"] += r["n_preserved"] == "0"
        conf[b]["signal_peptide_lost_everywhere"] += (
            r["subgroup"] == "signal_peptide_lost_everywhere"
        )
    with (S4E9 / "constructs_sp_swap.tsv").open() as fh:
        swap = list(csv.DictReader(fh, delimiter="\t"))
    out = {
        "tranche": "l1",
        "members": len(rows),
        "subgroups": dict(sorted(collections.Counter(r["subgroup"] for r in rows).items())),
        "rule_verdicts_over_repredicted_members": rule_verdicts | {"n": len(rep_rows)},
        "reprediction_gated": {
            "members_with_a_pick": len(gated),
            "pick_before_fold": sum(r["gated_inside_fold"] == 0 for r in gated),
            "pick_inside_fold_tolerance": sum(r["gated_inside_fold"] == 1 for r in gated),
            "pick_softened_over_bound": sum(r["gated_softened"] == 1 for r in gated),
            "pick_is_cleavage_plus_1": sum(
                r["gated_site"] == r["cleavage_plus_1_site"] for r in gated
            ),
            "pick_is_disorder_anchored": sum(
                str(r["gated_site"]) == str(r["disorder_anchored_site"]) for r in gated
            ),
            "offset_histogram": dict(
                sorted(collections.Counter(r["gated_offset"] for r in gated).items())
            ),
        },
        "native_cleavage_confidence_vs_outcome": {b: dict(c) for b, c in sorted(conf.items())},
        "exceptions": {
            sg: sorted(
                f"{r['candidate_id']} {r['gene_symbol']}" for r in rows if r["subgroup"] == sg
            )
            for sg in (
                "signal_peptide_lost_everywhere",
                "cleavage_moved_upstream_only",
                "cleavage_crosses_junction_everywhere",
                "not_repredicted",
            )
        },
        "igk_swap_arms_built": [
            {
                k: s[k]
                for k in (
                    "construct",
                    "gene",
                    "native_cleavage_p",
                    "igk_cleavage_p",
                    "reprediction_run",
                )
            }
            for s in swap
        ],
    }
    (HERE / "nterm_build_options.json").write_text(json.dumps(out, indent=1, sort_keys=True) + "\n")
    print(json.dumps({k: v for k, v in out.items() if k != "exceptions"}, indent=1))
    print({k: len(v) for k, v in out["exceptions"].items()})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
