"""Why each held N-after-SP member is held, its best candidate, and what the rule still allows.

For the rescue round Becca asked for on 2026-11-14 ("let's find replacements where things were
held if you can and it makes sense"), relayed by the coordinator. Per held member it lists every
tagged form that keeps a SignalP cleavage site upstream of the cassette (native candidates from
the sweep and the pinned runs, Igκ forms from the round-2 and overnight inputs), each form's
refold verdict and metrics from ``tag_refold_verdicts.tsv``, and:

- ``unfolded_passing_forms``: SignalP-passing forms never folded (rescue route a);
- ``igk_route``: whether an Igκ form is open under the evidence rule (rescue route b);
- ``near_miss``: one failing criterion only, an lDDT within 0.03 of its bar or 5-6 packed
  epitope residues. Reported only: the bars are not lowered.

Joined on design id (member ids are renumbered between cuts). Written by the Claude Code
session nterm-reviewer-package, 2026-11-14.

    uv run --locked python hold_review.py
"""

from __future__ import annotations

import collections
import csv
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
S4E9 = HERE.parents[1]
sys.path.insert(0, str(HERE))

import nterm_decide as nd  # noqa: E402
import nterm_overnight as no  # noqa: E402
import nterm_placement as npl  # noqa: E402

OUT = HERE / "hold_review.tsv"
BARS = {"native": 0.85, "site": 0.70}


def _tsv(path: Path) -> list[dict]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _f(x: str) -> float | None:
    try:
        return float(x)
    except TypeError, ValueError:
        return None


def shortfall(v: dict) -> tuple[int, float, str]:
    """(failing criteria, worst normalised shortfall, summary) for one verdict row."""
    reasons = [r.strip() for r in (v.get("verdict_reasons") or "").split(";") if r.strip()]
    worst, notes = 0.0, []
    for r in reasons:
        m = re.search(r"([0-9.]+) < ([0-9.]+)", r)
        if m:
            got, bar = float(m.group(1)), float(m.group(2))
            worst = max(worst, (bar - got) / bar)
            notes.append(f"{r}")
            continue
        m = re.search(r"(\d+) epitope residues confidently packed", r)
        if m:
            worst = max(worst, (int(m.group(1)) - 4) / 5)
            notes.append(r)
            continue
        worst = max(worst, 1.0)
        notes.append(r)
    return len(reasons), worst, "; ".join(notes)


def near_miss(v: dict) -> bool:
    reasons = [r.strip() for r in (v.get("verdict_reasons") or "").split(";") if r.strip()]
    if len(reasons) != 1:
        return False
    r = reasons[0]
    m = re.search(r"([0-9.]+) < ([0-9.]+)", r)
    if m and "lDDT" in r:
        return float(m.group(2)) - float(m.group(1)) <= 0.03
    m = re.search(r"(\d+) epitope residues confidently packed", r)
    return bool(m) and int(m.group(1)) <= 6


def main() -> int:
    from deliverome_analysis.s3 import processed

    calls = no.prior_calls(processed)
    members, _ = npl.load_members(processed)
    by_design = {m["design_id"]: m for m in members.values()}
    verdicts = {}
    for r in _tsv(S4E9 / "tag_refold_verdicts.tsv"):
        if r["verdict"] in ("pass", "fail"):
            verdicts.setdefault(r["tagged_sha256"], r)
    census = {r["design_id"]: r for r in _tsv(HERE / "nterm_members.tsv")}
    review = {r["design_id"]: r for r in _tsv(HERE / "sp_disagree_review.tsv")}
    forms = collections.defaultdict(dict)  # design -> sha -> form
    for p in _tsv(HERE / "nterm_positions.tsv"):
        if p["verdict"] == "preserved":
            forms[p["design_id"]][p["tagged_sha256"]] = {
                "arm": "native",
                "n": int(p["after_residue"]),
                "cp": p["tagged_cleavage_probability"],
            }
    for path in (no.INPUTS_TSV, nd.INPUTS_TSV):
        rows = _tsv(path)
        parent = {
            (r["design_id"], r["arm"]): r for r in rows if r["arm"] in ("native", "igk_native")
        }
        for t in rows:
            if t["arm"] not in ("tagged", "igk_tagged"):
                continue
            par = parent.get((t["design_id"], "native" if t["arm"] == "tagged" else "igk_native"))
            pc, tc = (
                calls.get(par["sequence_sha256"]) if par else None,
                calls.get(t["sequence_sha256"]),
            )
            if not (pc and tc and pc["prediction"] == "SP"):
                continue
            n = int(t["after_residue"])
            if npl.verdict(pc, tc, n) == npl.PRESERVED and n > pc["cleavage_site"]:
                forms[t["design_id"]][t["sequence_sha256"]] = {
                    "arm": "igk" if t["arm"] == "igk_tagged" else "native",
                    "n": int(t["native_after_residue"] or n),
                    "cp": tc["cleavage_probability"],
                }
    out = []
    for d in _tsv(HERE / "nterm_decisions.tsv"):
        if d["arm"] != "no_nterm_arm":
            continue
        did, reason = d["design_id"], d["reason"]
        fs = forms.get(did, {})
        folded = {s: verdicts[s] for s in fs if s in verdicts}
        unfolded = [s for s in fs if s not in verdicts]
        fails = {s: v for s, v in folded.items() if v["verdict"] == "fail"}
        passes = [s for s, v in folded.items() if v["verdict"] == "pass"]
        best = min(fails, key=lambda s: shortfall(fails[s])[:2]) if fails else None
        bv = fails.get(best, {})
        c = census.get(did, {})
        native_sp = c.get("signalp_prediction") == "SP" or (
            did in by_design
            and (calls.get(nd.sha(by_design[did]["sequence"]), {}).get("prediction") == "SP")
        )
        if "pro-peptide" in reason:
            why = "pro-peptide relocation " + (
                "loses SignalP cleavage" if "SignalP" in reason.split(";")[-1] else "fails refold"
            )
        elif "refold failed" in reason:
            why = "refold failed at every SignalP-passing form"
        elif "signal anchor" in reason:
            why = "evidence rule: signal anchor / first TM"
        elif "only DeepTMHMM v2" in reason or "no independent source" in reason:
            why = "evidence rule: no curated cleaved signal peptide"
        elif "no tagged form" in reason:
            why = "no tagged form keeps a SignalP cleavage site"
        else:
            why = "other: " + reason[:80]
        igk = (
            "open: native SignalP SP and no Igκ form folded"
            if native_sp and not any(f["arm"] == "igk" for f in fs.values()) and "refold" not in why
            else "not applicable"
            if "refold" in why
            else "closed: "
            + ("native SignalP calls no signal peptide" if not native_sp else "Igκ forms tried")
        )
        out.append(
            {
                "design_id": did,
                "candidate_id": d["candidate_id"],
                "gene_symbol": d["gene_symbol"],
                "candidate_set": d["candidate_set"],
                "hold_reason": why,
                "native_signalp_sp": int(bool(native_sp)),
                "sp_review_verdict": review.get(did, {}).get("verdict", ""),
                "signalp_passing_forms": len(fs),
                "folded_pass": len(passes),
                "folded_fail": len(fails),
                "unfolded_passing_forms": len(unfolded),
                "best_form": (f"{fs[best]['arm']} after {fs[best]['n']}" if best else ""),
                "best_lddt_native": bv.get("lddt_confident", ""),
                "best_lddt_site": bv.get("lddt_site", ""),
                "best_epitope_packed": bv.get("epitope_confident_contact_res", ""),
                "best_fail_reasons": bv.get("verdict_reasons", ""),
                "near_miss": int(bool(bv) and near_miss(bv)),
                "igk_route": igk,
                "decision_reason": reason,
            }
        )
    with OUT.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out)
    print(len(out), collections.Counter(r["hold_reason"] for r in out))
    print(
        "near misses",
        sum(r["near_miss"] for r in out),
        "with unfolded forms",
        sum(r["unfolded_passing_forms"] > 0 for r in out),
    )
    print("passing forms present", [r["gene_symbol"] for r in out if r["folded_pass"]])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
