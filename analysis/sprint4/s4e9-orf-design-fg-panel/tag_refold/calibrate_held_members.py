"""Control-relative refold verdicts for held internal members (rescue round, 2026-11-14).

Small 4-TM folds (MS4A family, EMPs, PMP22, tetraspanin-likes) failed the absolute
native-region lDDT bar (0.85) with every tagged arm tried. This asks whether the bar or the tag
is at fault: each member's untagged designed protein is folded with exactly the tagged arms'
pipeline (``tag_refold_calibration_controls.tsv`` -> ``build_refold_inputs.py`` -> rung
``control``) and scored by ``score_tagged_folds.py`` against the same native reference.

Rule (coordinator under Becca's delegation, 2026-11-14; rejected alternative: lowering 0.85
for every member): where the untagged control itself scores below 0.85, the absolute bar is
uninformative for that member, and its native-region bar becomes
``max(control lDDT - 0.05, 0.70)``. The site bar (0.70) and the epitope rules are unchanged.
Where the control reaches 0.85, nothing changes and the member stays held.

An arm passes control-relative when its only failure reasons are the native-region lDDT one
and its lDDT reaches the member's bar. The original verdicts are not touched; this writes

- ``tag_refold_control_relative.tsv``: one row per tagged arm of a calibrated member, with
  its original verdict beside the control-relative one;
- ``tag_refold_held_calibration.tsv``: one row per member, with the control, the bar, the
  best arm and which arms, if any, pass under the rule.

Implemented by Claude Code session l1-tag-refold.

    uv run --locked python analysis/sprint4/s4e9-orf-design-fg-panel/tag_refold/calibrate_held_members.py
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
EXP = HERE.parent
sys.path.insert(0, str(HERE))

from score_tagged_folds import LDDT_PASS, LOCAL_LDDT_PASS  # noqa: E402

DECIDED = "coordinator under Becca's delegation, 2026-11-14"
FLOOR = 0.70
MARGIN = 0.05


def tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def main() -> None:
    members = tsv(EXP / "tag_refold_calibration_controls.tsv")
    queue = tsv(EXP / "tag_refold_queue.tsv")
    metrics = {r["fold_id"]: r for r in tsv(EXP / "tag_refold_metrics.tsv")}
    arms_out, member_out = [], []
    for mem in members:
        d = mem["design_id"]
        ctrl = [q for q in queue if q["design_id"] == d and q["kind"] == "control"]
        if len(ctrl) != 1:
            raise SystemExit(f"{d}: no control row in the queue")
        c = metrics.get(ctrl[0]["fold_id"], {})
        if not c.get("lddt_confident"):
            # reported, not guessed: the member keeps its original verdicts until scored
            member_out.append({"design_id": d, "gene_symbol": mem["gene_symbol"],
                               "species": mem["species"], "control_fold_id": ctrl[0]["fold_id"],
                               "rule_applies": "control not scored yet"})  # fmt: skip
            continue
        cl = float(c["lddt_confident"])
        relative = cl < LDDT_PASS
        bar = round(max(cl - MARGIN, FLOOR), 3) if relative else LDDT_PASS
        arms = [q for q in queue if q["design_id"] == d and q["kind"] in ("tagged", "superseded")
                and q["fold_id"] in metrics and metrics[q["fold_id"]].get("lddt_confident")]  # fmt: skip
        rows = []
        for q in arms:
            m = metrics[q["fold_id"]]
            reasons = [x for x in m["verdict_reasons"].split("; ") if x]
            other = [x for x in reasons if not x.startswith("native-region lDDT")]
            lddt = float(m["lddt_confident"])
            if not relative:
                cr = m["verdict"]
            else:
                cr = "pass" if not other and lddt >= bar else "fail"
            why = other + (
                [f"native-region lDDT {lddt} < control-relative bar {bar}"] if lddt < bar else []
            )
            rows.append({
                "design_id": d, "gene_symbol": mem["gene_symbol"], "species": mem["species"],
                "fold_id": q["fold_id"], "kind": q["kind"], "rung": q["rung"],
                "source_id": q["source_id"], "insert_after_native_residue": q["insert_after_native_residue"],
                "tagged_sha256": q["tagged_sha256"], "lddt_confident": m["lddt_confident"],
                "lddt_site": m["lddt_site"], "spytag003_rsa": m["spytag003_rsa"],
                "alfa_rsa": m["alfa_rsa"], "epitope_confident_contact_res": m["epitope_confident_contact_res"],
                "original_verdict": m["verdict"], "original_reasons": m["verdict_reasons"],
                "control_fold_id": ctrl[0]["fold_id"], "control_lddt": c["lddt_confident"],
                "native_bar": bar, "rule_applies": "yes" if relative else "no (control reaches 0.85)",
                "control_relative_verdict": cr,
                "control_relative_reasons": "; ".join(why) if cr == "fail" else "",
                "decided_by": DECIDED, "implemented_by": "l1-tag-refold",
            })  # fmt: skip
        rows.sort(key=lambda r: -float(r["lddt_confident"]))
        arms_out += rows
        best = rows[0] if rows else {}
        live = [r for r in rows if r["kind"] == "tagged"]
        passing = [r for r in live if r["control_relative_verdict"] == "pass"]
        member_out.append({
            "design_id": d, "gene_symbol": mem["gene_symbol"], "species": mem["species"],
            "control_fold_id": ctrl[0]["fold_id"], "control_lddt": c["lddt_confident"],
            "control_reference": c.get("native_model", ""), "native_bar": bar,
            "rule_applies": "yes" if relative else "no (control reaches 0.85; member stays held)",
            "arms_scored": len(rows), "live_arms": len(live),
            "best_arm": best.get("source_id", ""), "best_arm_rung": best.get("rung", ""),
            "best_arm_lddt": best.get("lddt_confident", ""), "best_arm_site_lddt": best.get("lddt_site", ""),
            "best_arm_spytag003_rsa": best.get("spytag003_rsa", ""), "best_arm_alfa_rsa": best.get("alfa_rsa", ""),
            "best_arm_epitope_contacts": best.get("epitope_confident_contact_res", ""),
            "best_arm_original_reasons": best.get("original_reasons", ""),
            "passing_arms_control_relative": ", ".join(
                f"{r['source_id']} ({r['rung']}, lDDT {r['lddt_confident']}, {r['tagged_sha256'][:12]})" for r in passing
            ),
            "decided_by": DECIDED, "implemented_by": "l1-tag-refold",
        })  # fmt: skip
    write(EXP / "tag_refold_control_relative.tsv", arms_out)
    keys = max(member_out, key=len).keys()
    write(
        EXP / "tag_refold_held_calibration.tsv",
        [{k: m.get(k, "") for k in keys} for m in member_out],
    )
    n = sum(1 for m in member_out if m.get("passing_arms_control_relative"))
    print(
        f"{len(member_out)} members, {len(arms_out)} arms; {n} members with a passing live arm under the rule"
    )
    print(f"site bar {LOCAL_LDDT_PASS} and epitope rules unchanged")


if __name__ == "__main__":
    main()
