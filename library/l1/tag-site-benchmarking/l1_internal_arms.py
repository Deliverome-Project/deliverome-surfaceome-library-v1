#!/usr/bin/env python3
"""Internal-rung arms: how many tag sites each member gets, and which (the S5E13 order rule).

`l1_internal_decisions.py` decides each member's first site. This applies the order rule from the
v5 design (`recommended_sites_to_order`, written by tag_site_model_v4.design since 2026-10-15)
and Becca's 2026-11-15 additions, as the coordinator adopted them under her delegation:

    order rule     best region ceiling (absolute forest, within the r8-gated regions) against the
                   out-of-fold ceilings of regions hosting a published internal tag (v5: median
                   0.6951, q25 0.5574): below the median -> 3 sites; at or above it with the
                   runner-up region >= 0.90 of the best -> 2; otherwise 1. Capped by room (a
                   second site in a different region, or >= 25 residues away in the same one).
    all bad        best ceiling < q25: the two best loop arms plus a snorkel where the C-terminus
                   is cytoplasmic and the snorkel passes refold; otherwise a third loop. Never
                   more than three arms.
    never alone    a snorkel is never a member's only arm. A snorkel member gets its loop arms by
                   the same rule; if none passes refold, its best contact-clean, gate-passing
                   loop arm is built beside the snorkel and flagged; with no clean loop at all
                   the member is held and listed for Becca.
    every arm      passes the contact veto, the membrane-side gate and refold (the flagged loop
                   beside a snorkel is the one exception, and it is labelled)
    held (22)      unchanged until the control-relative refold bar arrives (coordinator,
                   2026-11-15); no arm is built for them here

Thresholds: the v4 session's choice, 2026-10-15, adopted by the coordinator 2026-11-15, not yet
confirmed by Becca. Schema agreed with l1-integration and l1-reviewer-package (2026-11-15): one
row per arm; `arm_id` `<design_id>-aK`, K assigned in creation order and never reused (every
pre-existing arm is -a1); `arm_role` loop_best / loop_2 / loop_3 / snorkel_companion; `decision`
build / hold; `agent_review` and `agent_reason` are left empty for the agent session's advisory
verdict. `construct_key` keeps the builder's `<design_id>:internal<N>` / `<design_id>:snorkel`.

    uv run --locked python l1_internal_arms.py     # writes l1_internal_arm_decisions.tsv/.json
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(REPO / "src"))

import contact_veto as CV  # noqa: E402
import l1_internal_decisions as L  # noqa: E402
import tag_site_model_v4 as V4  # noqa: E402
import tag_site_model_v5 as V5  # noqa: E402
from build_internal_tagged import CASSETTE_AA  # noqa: E402
from site_count_proposal import designs  # noqa: E402
from snorkel_design import MODULE  # noqa: E402

#: While a revision waits on refold verdicts, the table is written as a pending copy: the s4e9
#: build pins the canonical l1_internal_arm_decisions.tsv by SHA-256, so the canonical files stay
#: as the build last pinned them until the revision is final and the build re-pins (coordinator,
#: 2026-11-15). Set PENDING = False to promote.
PENDING = False
_SUFFIX = "_pending" if PENDING else ""
OUT = HERE / f"l1_internal_arm_decisions{_SUFFIX}.tsv"
OUT_JSON = HERE / f"l1_internal_arm_decisions{_SUFFIX}.json"
OUT_MEMBERS = HERE / f"l1_internal_arm_members{_SUFFIX}.tsv"
MODEL_VERSION = "v5"
MAX_ARMS = 3
CONTROL_RELATIVE = (
    HERE / "data" / "refold_2026-11-15" / "l1-tag-refold-control-relative-arms-2026-11-15.tsv"
)
DECIDED_BY = (
    "order rule: v4 session 2026-10-15; adopted with Becca's 2026-11-15 additions (all bad -> "
    "tag + snorkel, never a snorkel alone, cap 3) by the coordinator under her delegation"
)
FLAG_LOOP = (
    "loop arm fails refold; built so the snorkel is not the member's only arm (Becca 2026-11-15)"
)
#: gate verdicts from the membrane-side gate session for sites with no r8 row (routed-in members)
GATE_SESSION = {
    ("l1-p0567", 34): "pass",
    ("l1-p0567", 458): "pass",
    ("l1-p2987", 33): "pass",
    ("l1-p2987", 458): "pass",
}


#: members built on another rung, so their internal rows are holds (coordinator, 2026-11-15)
BUILT_ELSEWHERE = {"l1-p1362": "built on the residue-1 multipass rung, #361"}
#: refold verdicts of record (s4e9 tag_refold tables on main), keyed by tagged_sha256
TAG_REFOLD_VERDICTS = REPO / "library/l1/orf-design" / "tag_refold_verdicts.tsv"


#: the annotation rule's per-arm verdicts (annotation_rule.py, on the #370 table)
ANNOTATION_RULE = HERE / "annotation_rule_arms.tsv"
ANNOTATION_DECIDED = (
    "coordinator under Becca's delegation, 2026-11-15: disulfide Cys, insertion-intolerant "
    "position or a split annotated sequon in n-2 ... n+3 is replaced; a glycan or predicted "
    "sequon beside the junction is flagged"
)


def apply_annotation_rule(arms: list, prot: dict) -> list:
    """Replace, drop or flag built loop arms under the annotation rule (annotation_rule.py).

    A replacement is taken only when its refold verdict is pass (s4e9 tables or the s5e13
    verdicts seen); until then the arm stays, flagged as pending. With no replacement, or a
    failed one, the arm is held when the member keeps another built arm, and otherwise kept,
    flagged for Becca.
    """
    if not ANNOTATION_RULE.is_file():
        return arms
    rule = pd.read_csv(ANNOTATION_RULE, sep="\t")
    rule = rule[rule["case"].isin(["replace", "flag"])]
    by_site = {(r.design_id, int(r.insert_after_residue)): r for r in rule.itertuples()}
    v = pd.read_csv(TAG_REFOLD_VERDICTS, sep="\t", low_memory=False)
    verdict = dict(zip(v["tagged_sha256"], v["verdict"], strict=False))
    seen = pd.read_csv(L.VERDICTS_SEEN, sep="\t")
    verdict |= {k: x for k, x in zip(seen["tagged_sha256"], seen["verdict"], strict=False)
                if k not in verdict}  # fmt: skip
    builds = {}
    for x in arms:
        # count loops only: dropping a member's last loop would leave a snorkel alone
        if x.get("decision") == "build" and x.get("arm_kind") == "loop":
            builds[x["design_id"]] = builds.get(x["design_id"], 0) + 1
    out = []
    for x in arms:
        n = x.get("insert_after_residue")
        key = (x["design_id"], int(n)) if x.get("arm_kind") == "loop" and n == n and n else None
        r = by_site.get(key) if x.get("decision") == "build" else None
        if r is None:
            out.append(x)
            continue
        flag = x.get("flag") if isinstance(x.get("flag"), str) else ""
        if r.case == "flag":
            out.append({**x, "flag": "; ".join(filter(None, [flag,
                        f"annotation (flag only): {r.annotations}"]))})  # fmt: skip
            continue
        why = f"annotation rule: {r.annotations}"
        rep = r.replacement if r.replacement == r.replacement else None
        vv = verdict.get(r.replacement_tagged_sha256, "not_folded") if rep else None
        if rep and vv == "pass":
            m = int(rep)
            seq = prot[x["sequence_sha256"]]
            t = hashlib.sha256((seq[:m] + CASSETTE_AA + seq[m:]).encode()).hexdigest()
            out.append({**x, "insert_after_residue": m, "construct_key": f"{x['design_id']}"
                        f":internal{m}", "refold_tagged_sha256": t, "model_version": "v6",
                        "flag": "; ".join(filter(None, [flag, f"replaced {int(n)} -> {m}: "
                                                 f"{why} ({r.replacement_basis})",
                                                 getattr(r, "replacement_flag", "")
                                                 if isinstance(getattr(r, "replacement_flag", ""),
                                                               str) else ""])),
                        "reason": why, "decided_by": ANNOTATION_DECIDED})  # fmt: skip
        elif rep and vv == "not_folded":
            out.append({**x, "flag": "; ".join(filter(None, [flag, f"{why}; replacement after "
                        f"{int(rep)} pending refold"]))})  # fmt: skip
        elif builds[x["design_id"]] > 1:
            res = "no clean eligible site" if not rep else f"replacement after {int(rep)} {vv}"
            out.append(
                {
                    **x,
                    "decision": "hold",
                    "reason": f"{why}; {res}; the member keeps another loop arm, so this one is dropped",
                    "decided_by": ANNOTATION_DECIDED,
                }
            )
            builds[x["design_id"]] -= 1  # fmt: skip
        else:
            res = "no clean eligible site" if not rep else f"replacement after {int(rep)} {vv}"
            out.append({**x, "flag": "; ".join(filter(None, [flag, f"{why}; {res}; the "
                        "member's only loop arm, kept for Becca (never a snorkel alone)"]))})  # fmt: skip
    return out


#: arms with a hard-tier contact in the final literature table (lit_final_check.py), resolved
#: without folding (coordinator, 2026-11-15): dropped where the member keeps another loop arm,
#: otherwise kept, flagged and listed for Becca
LIT_HARD = {
    ("l1-p0177", 264): "small ligand at 267 (structure)",
    ("l1-p0512", 263): "protein ligand at 261-265 (structure)",
    ("l1-p0995", 408): "small ligand at 409 and 411 (structure)",
}
LIT_HARD_DECIDED = (
    "coordinator under Becca's delegation, 2026-11-15: hard literature contacts from the final "
    "table are resolved without folding"
)


def apply_lit_hard(arms: list) -> list:
    loops = {}
    for x in arms:
        if x.get("decision") == "build" and x.get("arm_kind") == "loop":
            loops[x["design_id"]] = loops.get(x["design_id"], 0) + 1
    out = []
    for x in arms:
        n = x.get("insert_after_residue")
        k = (x["design_id"], int(n)) if x.get("decision") == "build" and n == n and n else None
        if k not in LIT_HARD:
            out.append(x)
            continue
        why = f"hard literature contact (final table 2026-11-15): {LIT_HARD[k]}"
        if loops[x["design_id"]] > 1:
            loops[x["design_id"]] -= 1
            out.append({**x, "decision": "hold", "decided_by": LIT_HARD_DECIDED,
                        "reason": f"{why}; the member keeps another loop arm, so this one is "
                        "dropped"})  # fmt: skip
        else:
            f = x.get("flag") if isinstance(x.get("flag"), str) else ""
            out.append({**x, "flag": "; ".join(filter(None, [f, f"{why}; the member's only "
                        "loop arm, kept for Becca"]))})  # fmt: skip
    return out


#: agent-reviewer sites adopted under the coordinator's rulings (agent_site_adopt.py)
AGENT_ADOPTIONS = HERE / "agent_site_adoptions.tsv"


def apply_agent_adoptions(arms: list, members: list) -> list:
    """Replace, substitute or add loop arms from the agent reviewer's passing sites.

    replace     the built loop at old_site moves to new_site (it loses the old site's flag)
    substitute  the member's refold-failing loop beside its snorkel moves to new_site
    add         a new loop arm; a held member's placeholder row is replaced and it is built
    """
    if not AGENT_ADOPTIONS.is_file():
        return arms
    ad = pd.read_csv(AGENT_ADOPTIONS, sep="\t")
    ad = ad[ad["action"].isin(["replace", "substitute", "add"])]
    out = list(arms)

    def moved(x, r, old):
        n = int(r.new_site)
        return {**x, "insert_after_residue": n, "construct_key": f"{r.design_id}:internal{n}",
                "refold_tagged_sha256": r.tagged_sha256, "model_version": "v6",
                "flag": r.flag if isinstance(r.flag, str) else "",
                "reason": f"agent site: {old} -> {n}; {r.reason}", "decided_by": r.decided_by,
                "agent_review": f"agent reviewer pilot 2026-11-15: proposed after {n}",
                "agent_reason": r.agent_reason}  # fmt: skip

    for r in ad.itertuples(index=False):
        mine = [i for i, x in enumerate(out) if x["design_id"] == r.design_id]
        loops = [i for i in mine if out[i].get("decision") == "build"
                 and out[i].get("arm_kind") == "loop"]  # fmt: skip
        if r.action == "replace":
            i = [k for k in loops if int(out[k]["insert_after_residue"]) == int(r.old_site)]
            assert len(i) == 1, (r.design_id, r.old_site)
            out[i[0]] = moved(out[i[0]], r, int(r.old_site))
        elif r.action == "substitute":
            failing = [k for k in loops if "fails refold" in str(out[k].get("flag"))
                       or out[k].get("refold") == "fail"]  # fmt: skip
            assert len(failing) == 1, (r.design_id, [out[k].get("flag") for k in loops])
            k = failing[0]
            out[k] = moved(out[k], r, int(out[k]["insert_after_residue"]))
        else:
            held = [k for k in mine if out[k].get("decision") == "hold"]
            base = out[(loops or held or mine)[0]]
            new = moved({k: v for k, v in base.items() if k != "arm_id"}, r, "new arm")
            new.update(decision="build", arm_kind="loop", snorkel_design=None,
                       order_verdict=base.get("order_verdict"), design_rank=None,
                       region_index=None, contact_veto=False, gate="pass",
                       refold_verdict="pass")  # fmt: skip
            if not loops and held:  # a held member's placeholder gives way to its arm
                for k in sorted(held, reverse=True):
                    out.pop(k)
                for mm in members:
                    if mm["design_id"] == r.design_id:
                        mm["status"], mm["note"] = "built", r.reason
            out.append(new)
    # loop roles in arm order (again, after the agent's changes)
    roles = {}
    for x in out:
        if x.get("decision") == "build" and x.get("arm_kind") == "loop":
            k = roles.get(x["design_id"], 0)
            x["arm_role"] = ["loop_best", "loop_2", "loop_3"][min(k, 2)]
            roles[x["design_id"]] = k + 1
    return out


def refold_of_record(a: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Each arm's verdict from the refold tables, joined on the folded protein's sha.

    `refold_verdict` is the absolute verdict (native-region lDDT 0.85, site lDDT 0.70) of the
    construct folded with exactly this sha; `refold_disposition` says why a built arm stands:
    pass, pass_control_relative (the untagged control folds below 0.85 and the arm clears
    max(control - 0.05, 0.70); coordinator, 2026-11-15) or fail_flagged (a coordinator
    exception or the loop that keeps a snorkel company; the flag says which). Verdicts that
    l1-tag-refold has published to verdicts_loop but not yet to the s4e9 tables fall back to
    `refold_verdicts_seen.tsv`; a built arm with no verdict anywhere is returned as unfolded.
    """
    v = pd.read_csv(TAG_REFOLD_VERDICTS, sep="\t", low_memory=False)
    v = v[v["verdict"].isin(["pass", "fail"])].drop_duplicates("tagged_sha256")
    v = v.set_index("tagged_sha256")
    seen = pd.read_csv(L.VERDICTS_SEEN, sep="\t").drop_duplicates("tagged_sha256")
    seen = seen.set_index("tagged_sha256")
    rel = pd.read_csv(CONTROL_RELATIVE, sep="\t")
    rel = rel[rel["rule_applies"] == "yes"].drop_duplicates("tagged_sha256")
    rel = rel.set_index("tagged_sha256")
    verdict, source, reasons, disp, flags = [], [], [], [], []
    for x in a.itertuples(index=False):
        t, flag = x.refold_tagged_sha256, x.flag if isinstance(x.flag, str) else ""
        if isinstance(t, str) and t in v.index:
            vv, src, why = v.at[t, "verdict"], f"s4e9 tag_refold_verdicts:{v.at[t, 'source_id']}", (
                v.at[t, "verdict_reasons"])  # fmt: skip
        elif isinstance(t, str) and t in seen.index:
            vv, src, why = seen.at[t, "verdict"], "s5e13 refold_verdicts_seen (verdicts_loop)", (
                seen.at[t, "verdict_reasons"])  # fmt: skip
        else:
            vv, src, why = "not_folded", "", ""
        d = ""
        if x.decision == "build":
            if vv == "pass":
                d = "pass"
            elif t in rel.index and rel.at[t, "control_relative_verdict"] == "pass":
                d = "pass_control_relative"
                note = (f"passes the control-relative bar (untagged control lDDT "
                        f"{rel.at[t, 'control_lddt']:.3f}, bar {rel.at[t, 'native_bar']:.3f}); "
                        f"fails the absolute bar: {why}")  # fmt: skip
                flag = "; ".join(filter(None, [note, flag]))
            elif vv == "fail":
                d = "fail_flagged"
            else:
                d = "not_folded"
        verdict.append(vv)
        source.append(src)
        reasons.append(why if isinstance(why, str) else "")
        disp.append(d)
        flags.append(flag)
    a = a.assign(refold_verdict=verdict, refold_disposition=disp, refold_source=source,
                 refold_reasons=reasons, flag=flags)  # fmt: skip
    cols = list(a.columns)
    for c in ("refold_disposition", "refold_source", "refold_reasons"):
        cols.remove(c)
        cols.insert(cols.index("refold_verdict") + 1, c)
    a = a[cols]
    unfolded = a[(a["decision"] == "build") & (a["refold_disposition"] == "not_folded")]
    return a, unfolded


def tagged(seq: str, n: int) -> str:
    return hashlib.sha256((seq[:n] + CASSETTE_AA + seq[n:]).encode()).hexdigest()


REVISIONS = HERE / "arm_revisions.tsv"
ARM_IDS = HERE / "l1_internal_arm_ids.tsv"  # registry: an arm_id is never reused
PROVISIONAL_COMMIT = "89226cb2"  # the provisional table l1-integration test-built (#360)
#: coordinator exceptions, 2026-11-15
SPECIAL = {
    # CD20: two loop hypotheses beside the snorkel although both fail their bar
    "l1-p1249": {"loops": [(186, "fold near-miss: native-region lDDT 0.770 against its "
                                 "control-relative bar 0.795"),
                           (161, "fails site lDDT 0.49")],
                 "why": "coordinator exception to 'do not build arms that fail their bar': the "
                 "alternative is holding CD20 or a lone snorkel"},
    # snorkel members whose only loop sits on a literature-tier (soft) ligand contact
    "l1-p1543": {"loops": [(247, "literature, inferred evidence: MFaN (allosteric modulator) and "
                                 "carvacrol (agonist) at distance 0-2")],
                 "why": "soft contact under the 2026-11-15 split; loop passes refold"},
    "l1-p1894": {"loops": [(201, "literature, mutagenesis evidence: L-ascorbate (apparent Km) "
                                 "at distance 1")],
                 "why": "soft contact under the 2026-11-15 split; loop passes refold"},
}  # fmt: skip


def apply_revisions(arms: list, members: list, prot: dict, seen: pd.Series) -> list:
    """arm_revisions.tsv applied to the loop arms; SPECIAL members rebuilt."""
    import hashlib as _h

    rev = pd.read_csv(REVISIONS, sep="\t") if REVISIONS.is_file() else pd.DataFrame()
    out = []
    by_member = {}
    for x in arms:
        by_member.setdefault(x["design_id"], []).append(x)
    for did, lst in by_member.items():
        r = rev[rev["design_id"] == did] if len(rev) else rev
        if did in SPECIAL:
            sp = SPECIAL[did]
            base = next(x for x in lst)
            sha = base["sequence_sha256"]
            seq = prot[sha]
            new = []
            for n, why in sp["loops"]:
                t = _h.sha256((seq[:n] + CASSETTE_AA + seq[n:]).encode()).hexdigest()
                new.append({**{k: v for k, v in base.items() if k not in ("arm_id",)},
                            "arm_kind": "loop", "arm_role": "", "insert_after_residue": n,
                            "construct_key": f"{did}:internal{n}", "decision": "build",
                            "refold_verdict": seen.get(t, "not_folded"),
                            "refold_tagged_sha256": t, "flag": why, "reason": sp["why"],
                            "snorkel_design": ""})  # fmt: skip
            snork = _h.sha256((seq + MODULE).encode()).hexdigest()
            new.append({**{k: v for k, v in base.items() if k not in ("arm_id",)},
                        "arm_kind": "snorkel", "arm_role": "snorkel_companion",
                        "insert_after_residue": None, "construct_key": f"{did}:snorkel",
                        "decision": "build", "refold_verdict": "pass",
                        "refold_tagged_sha256": snork,
                        "snorkel_design": "l1_snorkel_designs/2026-11-14b", "flag": "",
                        "reason": sp["why"]})  # fmt: skip
            out.extend(new)
            continue
        keep = []
        for x in lst:
            if x.get("decision") != "build" or x.get("arm_kind") != "loop":
                keep.append(x)
                continue
            n = int(x["insert_after_residue"])
            hit = r[r["old_site"] == n] if len(r) else r
            act = hit["action"].iloc[0] if len(hit) else None
            why = hit["reason"].iloc[0] if len(hit) else ""
            if act == "drop":
                continue
            if act == "keep_flagged":
                x = {**x, "flag": "; ".join(filter(None, [x.get("flag") or "", why]))}
            if act == "replace":
                m = int(hit["new_site"].iloc[0])
                seq = prot[x["sequence_sha256"]]
                t = _h.sha256((seq[:m] + CASSETTE_AA + seq[m:]).encode()).hexdigest()
                v = seen.get(t, "not_folded")
                hard = why.startswith("hard") or "Esp3I" in why
                if v == "pass":
                    x = {**x, "insert_after_residue": m, "construct_key": f"{did}:internal{m}",
                         "refold_verdict": v, "refold_tagged_sha256": t,
                         "flag": f"replaced {n} -> {m}: {why}", "reason": why}  # fmt: skip
                elif hard:
                    continue
                else:
                    x = {**x, "flag": f"{why}; replacement {m} {v}, so the site stays"}
            keep.append(x)
        member_flags = []
        if len(r) and (r["action"] == "flag_member").any():
            member_flags.append(r.loc[r["action"] == "flag_member", "reason"].iloc[0])
        for a in r[r["action"] == "add"].itertuples(index=False) if len(r) else []:
            m = int(a.new_site)
            base = next((x for x in keep if x.get("arm_kind") == "loop"), None)
            if base is None:
                continue
            # the cap counts every built arm, the snorkel companion included (arm_revisions
            # counted loops only, so AQP3 and SLC15A2 were proposed a fourth arm)
            if sum(x.get("decision") == "build" for x in keep) >= MAX_ARMS:
                member_flags.append(f"{a.reason}; not added at {m}: member already has "
                                    f"{MAX_ARMS} arms")  # fmt: skip
                continue
            seq = prot[base["sequence_sha256"]]
            t = _h.sha256((seq[:m] + CASSETTE_AA + seq[m:]).encode()).hexdigest()
            v = seen.get(t, "not_folded")
            if v == "pass":
                keep.append({**{k: v for k, v in base.items() if k != "arm_id"},
                             "insert_after_residue": m, "construct_key": f"{did}:internal{m}",
                             "region_index": a.new_region, "design_rank": None,
                             "refold_verdict": "pass", "refold_tagged_sha256": t,
                             "flag": "", "reason": a.reason, "model_version": "v6"})  # fmt: skip
            else:
                member_flags.append(f"{a.reason}; not added at {m}: refold {v}")
        if member_flags:
            msg = "; ".join(member_flags)
            for x in keep:
                if x.get("arm_kind") == "loop":
                    x["flag"] = "; ".join(filter(None, [x.get("flag") or "", msg]))
        out.extend(keep)
    # loop roles in arm order
    roles = {}
    for x in out:
        if x.get("decision") == "build" and x.get("arm_kind") == "loop":
            k = roles.get(x["design_id"], 0)
            x["arm_role"] = ["loop_best", "loop_2", "loop_3"][min(k, 2)]
            roles[x["design_id"]] = k + 1
    for mm in members:
        if mm["design_id"] in SPECIAL:
            mm["status"], mm["note"] = "built", SPECIAL[mm["design_id"]]["why"]
    return out


def assign_arm_ids(arms: list) -> list:
    """Stable ids from the registry; a new construct gets the member's next unused K."""
    reg = (
        pd.read_csv(ARM_IDS, sep="\t")
        if ARM_IDS.is_file()
        else pd.DataFrame(columns=["arm_id", "design_id", "construct_key"])
    )
    if reg.empty:
        # seed from the provisional per-arm table (the version l1-integration test-built), so
        # every construct already there keeps its id and a dropped arm's K is retired
        import io

        from deliverome_analysis.provenance import read_pinned_text

        rel = "library/l1/tag-site-benchmarking/l1_internal_arm_decisions.tsv"
        txt = read_pinned_text(PROVISIONAL_COMMIT, rel, repo=HERE.parents[2])
        old = pd.read_csv(io.StringIO(txt), sep="\t")
        old["construct_key"] = old["construct_key"].fillna("hold")
        reg = old[["arm_id", "design_id", "construct_key"]].drop_duplicates("arm_id")
        reg.to_csv(ARM_IDS, sep="\t", index=False)
    pairs = zip(reg["design_id"], reg["construct_key"], strict=True)
    known = dict(zip(pairs, reg["arm_id"], strict=True))
    next_k = {}
    for a in reg["arm_id"]:
        d, k = a.rsplit("-a", 1)
        next_k[d] = max(next_k.get(d, 0), int(k))
    new_rows = []
    for x in arms:
        key = (x["design_id"], x.get("construct_key") or "hold")
        if key in known:
            x["arm_id"] = known[key]
            continue
        k = next_k.get(x["design_id"], 0) + 1
        next_k[x["design_id"]] = k
        x["arm_id"] = f"{x['design_id']}-a{k}"
        known[key] = x["arm_id"]
        new_rows.append({"arm_id": x["arm_id"], "design_id": key[0], "construct_key": key[1]})
    if new_rows:
        pd.concat([reg, pd.DataFrame(new_rows)], ignore_index=True).to_csv(
            ARM_IDS, sep="\t", index=False
        )
    return arms


RESCUE_ATTRIBUTION = (
    "rescue round, coordinator under Becca's 2026-11-12 delegation and her 2026-11-14 request"
)
RESCUE_GATE = HERE / "data" / "rescue_round_gate.tsv"  # membrane-side gate, PR #354
RESCUE_REFOLD = (
    HERE / "data" / "refold_2026-11-15" / "l1-tag-refold-rescue-74-verdicts-2026-11-15.tsv"
)
SITE_CAVEAT = {
    ("l1-p2342", 40): "rests on the N-terminus being resolved: 11 residues from a possible "
    "signal anchor (7-29); site 187 is preferred",
    ("l1-p2628", 40): "rests on the N-terminus being resolved: 11 residues from a possible "
    "signal anchor (7-29); site 187 is preferred",
}


def rescue_arms(arms: list, members: list) -> None:
    """Held non-internal members routed internal by the rescue round (rescue_round.py)."""
    sites = pd.read_csv(HERE / "rescue_round_sites.tsv", sep="\t")
    gate = pd.read_csv(RESCUE_GATE, sep="\t").set_index(["design_id", "insert_after_residue"])
    ref = pd.read_csv(RESCUE_REFOLD, sep="\t").drop_duplicates("tagged_sha256")
    ref = ref.set_index("tagged_sha256")
    built_genes_human = {
        x["gene"] for x in arms if x.get("decision") == "build"
    }  # human internal members already built
    per = {}
    for r in sites.sort_values(["design_id", "rescue_rank"]).itertuples(index=False):
        n = int(r.insert_after_residue)
        gv = gate.loc[(r.design_id, n), "verdict"] if (r.design_id, n) in gate.index else "missing"
        rv = ref.loc[r.tagged_sha256, "verdict"] if r.tagged_sha256 in ref.index else "not_folded"
        flags = []
        if gv in ("not_placeable", "review_single_pass_geometry"):
            flags.append(f"gate {gv}: rests on refold and distance from the membrane")
        if (r.design_id, n) in SITE_CAVEAT:
            flags.append(SITE_CAVEAT[(r.design_id, n)])
        ok = rv == "pass" and gv in ("pass", "not_placeable", "review_single_pass_geometry")
        per.setdefault(r.design_id, []).append((r, n, gv, rv, ok, "; ".join(flags)))
    # a gene's human member counts as built unless it is held: a held internal member, or a
    # rescue member with no passing site. Members on other rungs (N-terminal, C-terminal,
    # multipass) are built there, which the ortholog rule accepts (PCDHB16 human l1-p1602).
    tp = L.full_head_paths()
    held = {x["design_id"] for x in arms if x.get("decision") == "hold"}
    held |= {d for d, lst in per.items() if not any(x[4] for x in lst)}
    human = tp[tp["species"] == "human"]
    human_built = set(built_genes_human) | set(
        human.loc[~human["design_id"].isin(held), "gene_symbol"]
    )
    for did, lst in per.items():
        r0 = lst[0][0]
        passing = [x for x in lst if x[4]]
        ortholog_blocked = r0.species != "human" and r0.gene not in human_built
        status = "built" if passing and not ortholog_blocked else "held"
        note = (
            "ortholog: its human member is not built (Becca 2026-11-15)"
            if ortholog_blocked
            else ("" if passing else "no rescue site passes refold and gate")
        )
        rows = passing if status == "built" else []
        for k, (r, n, gv, rv, _ok, flag) in enumerate(rows, 1):
            arms.append(
                {
                    "arm_id": f"{did}-a{k}",
                    "construct_key": f"{did}:internal{n}",
                    "design_id": did,
                    "member_id": r.member_id,
                    "gene": r.gene,
                    "sequence_sha256": r.sequence_sha256,
                    "arm_role": ["loop_best", "loop_2", "loop_3"][k - 1],
                    "arm_kind": "loop",
                    "insert_after_residue": n,
                    "decision": "build",
                    "order_verdict": "rescue_round",
                    "model_version": MODEL_VERSION,
                    "design_rank": int(r.rescue_rank),
                    "contact_veto": False,
                    "gate": gv,
                    "refold_verdict": rv,
                    "refold_tagged_sha256": r.tagged_sha256,
                    "flag": flag,
                    "reason": f"routed internal by the rescue round; {r.model_scope}",
                    "agent_review": "",
                    "agent_reason": "",
                    "decided_by": RESCUE_ATTRIBUTION,
                    "implemented_by": "ranker-pu-segment session (Claude Code)",
                    "decided_on": "2026-11-15",
                }
            )
        if status == "held":
            arms.append(
                {
                    "arm_id": f"{did}-a1",
                    "design_id": did,
                    "member_id": r0.member_id,
                    "gene": r0.gene,
                    "sequence_sha256": r0.sequence_sha256,
                    "decision": "hold",
                    "order_verdict": "rescue_round",
                    "model_version": MODEL_VERSION,
                    "reason": note,
                    "decided_by": RESCUE_ATTRIBUTION,
                    "implemented_by": "ranker-pu-segment session (Claude Code)",
                    "decided_on": "2026-11-15",
                }
            )
        members.append(
            {
                "design_id": did,
                "member_id": r0.member_id,
                "gene": r0.gene,
                "decision": "rescue_round",
                "arms_built": len(rows),
                "has_snorkel": False,
                "has_flagged_loop": any(x[5] for x in rows),
                "status": status,
                "note": note,
            }
        )


def scores_and_veto_class(a: pd.DataFrame, dec: pd.DataFrame) -> pd.DataFrame:
    """Both models' scores at each loop junction, and the contact tier the junction falls in.

    v5 placed the built sites; v6 is the scorer of record for replacement and added sites
    (`model_version`). Scores are tag_site_model_v4.design's (rank forest z-scored within the
    protein, absolute forest raw, both smoothed over +/-2). `veto_class` is the strongest tier
    among contacts at n-2 ... n+3 (arm_revisions.contact_tiers): hard, soft, flag or none.
    """
    import tag_site_model_v6 as V6
    from arm_revisions import contact_tiers, hits
    from redesign_all import frames
    from rescue_round import score
    from v5_region_table import features

    fr = frames(dec)
    loops = a[(a["arm_kind"] == "loop") & a["insert_after_residue"].notna()]
    want = {fr[s]: s for s in loops["sequence_sha256"].unique() if s in fr}
    allx = features()
    # rescue-round members (single-pass and GPI) have their own feature table, keyed as in
    # rescue_round.py; their junctions lie outside v5's training distribution
    import rescue_round as RR

    if RR.FEATURES.is_file() and RR.PROVENANCE.is_file():
        pv = pd.read_csv(RR.PROVENANCE, sep="\t")
        # rescue members are not in the internal decisions table; the feature provenance
        # carries their sequence
        rk = dict(zip(pv["key"], pv["sequence_sha256"], strict=True))
        rf = pd.read_csv(RR.FEATURES, sep="\t", low_memory=False)
        rf = rf[rf["uniprot_acc"].isin(rk)]
        # a rescue member is scored on its rescue features only (its frame key, if any, is
        # the canonical row the rescue round replaced)
        want = {k: s for k, s in want.items() if s not in set(rk.values())} | rk
        allx = pd.concat([allx[~allx["uniprot_acc"].isin(rk)], rf], ignore_index=True)
    sub = allx[allx["uniprot_acc"].isin(want)]
    look = {}
    for tag, model in (("v5", V5), ("v6", V6)):
        s = score(sub, model)
        s["sha"] = s["uniprot_acc"].map(want)
        look[tag] = s.set_index(["sha", "residue"])[["rank_score", "absolute_score"]]
    ct = contact_tiers(pd.read_csv(HERE / "contacts_by_protein.tsv", sep="\t", low_memory=False))
    order = {"hard": 0, "soft": 1, "flag": 2}
    cols = {c: [] for c in ("v5_rank_score", "v5_absolute_score", "v6_rank_score",
                            "v6_absolute_score", "veto_class", "veto_contacts")}  # fmt: skip
    for x in a.itertuples(index=False):
        n = x.insert_after_residue
        if x.arm_kind != "loop" or pd.isna(n):
            for c in cols:
                cols[c].append("n/a" if c.startswith("veto") else None)
            continue
        n = int(n)
        for tag in ("v5", "v6"):
            k = (x.sequence_sha256, n)
            hit = k in look[tag].index
            for f in ("rank_score", "absolute_score"):
                v = look[tag].loc[k, f] if hit else np.nan
                v = float(v.iloc[0]) if isinstance(v, pd.Series) else float(v)
                cols[f"{tag}_{f}"].append(round(v, 4) if v == v else None)
        h = hits(ct, x.sequence_sha256, n)
        if h.empty:
            cols["veto_class"].append("none")
            cols["veto_contacts"].append("")
        else:
            cols["veto_class"].append(min(h["tier"], key=order.get))
            cols["veto_contacts"].append(
                "; ".join(sorted({f"{t}:{c}@{int(r)}" for t, c, r in
                                  zip(h["tier"], h["class"], h["residue"], strict=True)}))
            )  # fmt: skip
    a = a.assign(**cols)
    # put the new columns beside the model and veto columns they explain
    lead = list(a.columns[: a.columns.get_loc("region_ceiling") + 1])
    rest = [c for c in a.columns if c not in lead and c not in cols]
    return a[lead + list(cols) + rest]


def main() -> int:
    dec = pd.read_csv(L.OUT, sep="\t")
    d = designs()
    _i, _ids, prot = L.head_tables()
    mask = pd.read_csv(HERE / "data" / "mask_r8" / "l1_extracellular_mask_residues.tsv", sep="\t")
    mask = mask.set_index("sequence_sha256")
    seen = pd.read_csv(L.VERDICTS_SEEN, sep="\t").set_index("tagged_sha256")["verdict"]
    # snorkel arms are scored in the C-terminal verdict table, downloaded by the decisions build
    ct = HERE / "data" / "verdicts_cterm.tsv"
    if ct.is_file():
        c = pd.read_csv(ct, sep="\t")
        c = c[c["verdict"].isin(["pass", "fail"])].drop_duplicates("tagged_sha256")
        snork_seen = c.set_index("tagged_sha256")["verdict"]
        seen = pd.concat([seen, snork_seen[~snork_seen.index.isin(seen.index)]])
    # held members whose untagged control folds below 0.85 are judged against max(control - 0.05,
    # 0.70) (coordinator, 2026-11-15); their constructs carry that verdict instead
    rel = pd.read_csv(CONTROL_RELATIVE, sep="\t")
    rel = rel[rel["rule_applies"] == "yes"]
    control_of = rel.drop_duplicates("design_id").set_index("design_id")["control_lddt"].to_dict()
    rel_verdict = rel.set_index("tagged_sha256")["control_relative_verdict"]
    snorkel_rel = rel[rel["rung"] == "cterm"].set_index("design_id")["control_relative_verdict"]
    seen = pd.concat([rel_verdict, seen[~seen.index.isin(rel_verdict.index)]])
    contacts = CV.member_contacts(dec[["design_id", "sequence_sha256"]], prot, download=False)
    # evidence per arm (not all of it vetoes): nearest ligand, partner and homomer contact
    groups = CV.classified_contacts(dec[["design_id", "sequence_sha256"]], prot)
    import d1_snapshot

    schweke = {
        r["gene_symbol"]
        for r in d1_snapshot.query("SELECT gene_symbol FROM schweke_homomer_public")
    }

    def evidence(sha, n, gene):
        out = {}
        for grp in ("ligand", "partner", "homomer"):
            out.update(CV.nearest(groups, sha, n, grp) if n else {
                f"nearest_{grp}_distance": None, f"nearest_{grp}_residue": None,
                f"nearest_{grp}_source": "", f"nearest_{grp}_evidence": ""})  # fmt: skip
        out["homomer_assembly_structure"] = (
            "Schweke assembly model exists; interface residues not computed for members"
            if gene in schweke
            else "no assembly structure looked at: absence of homomer contacts is not evidence"
        )
        return out

    spec = json.loads(V5.SPEC.read_text())
    q25 = float(np.quantile(np.array(spec["published_region_ceilings_oof"]), 0.25))

    def gate(did, sha, n):
        if (did, n) in GATE_SESSION:
            return GATE_SESSION[(did, n)]
        if sha not in mask.index or int(mask.loc[sha, "usable"]) != 1:
            return "unknown"
        return "pass" if mask.loc[sha, "junction_mask"][n - 1] == "E" else "fail"

    def apart(n, reg, chosen):
        return all(reg != r or abs(n - m) >= V4.MIN_SEPARATION for m, r in chosen)

    arms, members = [], []
    for r in dec.itertuples(index=False):
        sha, did = r.sequence_sha256, r.design_id
        seq = prot.get(sha)
        g = d[d["sequence_sha256"] == sha].sort_values("rank")
        top = g.iloc[0] if len(g) else None
        n_rec = int(top["recommended_sites_to_order"]) if top is not None else 1
        ceiling = float(top["region_ceiling"]) if top is not None else np.nan
        all_bad = bool(ceiling < q25)
        region_of = dict(
            zip(g["insert_after_residue"].astype(int), g["region_index"], strict=False)
        )
        status, note, rows = "built", "", []

        def loop_row(
            n, rank, verdict_override=None, flag="", seq=seq, sha=sha, did=did, region_of=region_of
        ):
            t = tagged(seq, n)
            return {
                "arm_type": "loop",
                "insert_after_residue": n,
                "design_rank": rank,
                "region_index": region_of.get(n),
                "contact_veto": bool(CV.vetoed(contacts, sha, n)),
                "gate": gate(did, sha, n),
                "refold": verdict_override or seen.get(t, "not_folded"),
                "tagged_sha256": t,
                "flag": flag,
            }

        def clean_loops(g=g, loop_row=loop_row):
            """Design sites passing veto, gate and refold, in rank order, mutually apart."""
            out = []
            for x in g.itertuples(index=False):
                n = int(x.insert_after_residue)
                row = loop_row(n, int(x.rank))
                if row["contact_veto"] or row["gate"] == "fail" or row["refold"] != "pass":
                    continue
                if apart(n, row["region_index"], [(a["insert_after_residue"], a["region_index"])
                                                  for a in out]):  # fmt: skip
                    out.append(row)
            return out

        if r.decision == "hold_no_passing_arm":
            ctrl = control_of.get(did)
            loops = clean_loops() if ctrl is not None else []
            if loops:
                rows = loops[: min(n_rec, MAX_ARMS)]
                note = f"passes control-relative bar (control {ctrl:.3f})"
            elif ctrl is not None and snorkel_rel.get(did) == "pass":
                cands = [
                    loop_row(int(x.insert_after_residue), int(x.rank))
                    for x in g.itertuples(index=False)
                ]
                cands = [c for c in cands if not c["contact_veto"] and c["gate"] != "fail"]
                if cands:
                    best = cands[0]
                    best["flag"] = FLAG_LOOP
                    snork = {
                        "arm_type": "snorkel", "insert_after_residue": None, "design_rank": None,
                        "region_index": None, "contact_veto": False, "gate": "n/a",
                        "refold": "pass",
                        "tagged_sha256": hashlib.sha256((seq + MODULE).encode()).hexdigest(),
                        "flag": "",
                    }  # fmt: skip
                    rows = [best, snork]
                    note = f"snorkel passes control-relative bar (control {ctrl:.3f}); no loop does"
                else:
                    status, note = (
                        "held",
                        "no contact-clean loop arm; a snorkel is never built alone",
                    )
            else:
                status = "held"
                note = (
                    "no arm passes the refold bar: its untagged control folds at or above 0.85, "
                    "so the absolute bar holds"
                    if ctrl is None
                    else f"no arm passes the control-relative bar (control {ctrl:.3f})"
                )
        elif r.decision == "route_to_cterm_experimental_topology" or seq is None:
            status, note = "not_internal", r.decision
        elif r.decision == "snorkel":
            snork = {
                "arm_type": "snorkel",
                "insert_after_residue": None,
                "design_rank": None,
                "region_index": None,
                "contact_veto": False,
                "gate": "n/a",
                "refold": r.snorkel_verdict,
                "tagged_sha256": hashlib.sha256((seq + MODULE).encode()).hexdigest(),
                "flag": "",
            }
            loops = clean_loops()[: MAX_ARMS - 1]
            if not loops:
                cands = [
                    loop_row(int(x.insert_after_residue), int(x.rank))
                    for x in g.itertuples(index=False)
                ]
                cands = [c for c in cands if not c["contact_veto"] and c["gate"] != "fail"]
                if cands:
                    best = cands[0]
                    best["flag"] = FLAG_LOOP
                    loops = [best]
            if not loops:
                status, note = "held", "no contact-clean loop arm; a snorkel is never built alone"
            else:
                rows = loops + ([snork] if snork["refold"] == "pass" else [])
        else:  # internal loops (incl. routed-in members and the orientation-unconfirmed one)
            first = loop_row(int(r.insert_after_residue), int(r.pick_rank))
            first["refold"] = (
                r.refold_verdict if r.refold_verdict in ("pass", "fail") else (first["refold"])
            )
            chosen = [first]
            for c in clean_loops():
                if len(chosen) >= MAX_ARMS:
                    break
                if apart(c["insert_after_residue"], c["region_index"],
                         [(a["insert_after_residue"], a["region_index"]) for a in chosen]):  # fmt: skip
                    chosen.append(c)
            target = min(n_rec, MAX_ARMS)
            if all_bad:
                ok, _why = L._cyto_c_terminus(seq, L._gate_confirmed(sha))
                st = hashlib.sha256((seq + MODULE).encode()).hexdigest()
                sv = seen.get(st, "not_folded")
                if ok and sv == "pass":
                    rows = [
                        *chosen[:2],
                        {"arm_type": "snorkel", "insert_after_residue": None,
                         "design_rank": None, "region_index": None, "contact_veto": False,
                         "gate": "n/a", "refold": sv, "tagged_sha256": st, "flag": ""}
                    ]  # fmt: skip
                else:
                    rows = chosen[:MAX_ARMS]
                    note = (
                        f"all bad; snorkel {'not possible' if not ok else 'refold ' + sv}, so a "
                        "third loop"
                    )
            else:
                rows = chosen[:target]
        verdict = (
            "all_bad" if all_bad else {1: "order_1", 2: "order_2", 3: "order_3"}.get(n_rec, "")
        )
        loop_roles = iter(["loop_best", "loop_2", "loop_3"])
        for k, a in enumerate(rows, 1):
            snorkel = a["arm_type"] == "snorkel"
            key = f"{did}:snorkel" if snorkel else f"{did}:internal{a['insert_after_residue']}"
            arms.append(
                {
                    "arm_id": f"{did}-a{k}",
                    "construct_key": key,
                    "design_id": did,
                    "member_id": r.member_id,
                    "gene": r.gene,
                    "sequence_sha256": sha,
                    "arm_role": "snorkel_companion" if snorkel else next(loop_roles),
                    "arm_kind": a["arm_type"],
                    "insert_after_residue": a["insert_after_residue"],
                    "snorkel_design": "l1_snorkel_designs/2026-11-14b" if snorkel else "",
                    "decision": "build",
                    "order_verdict": verdict,
                    "order_why": None if top is None else top["why"],
                    "model_version": MODEL_VERSION,
                    "design_rank": a["design_rank"],
                    "region_index": a["region_index"],
                    "region_ceiling": None if top is None else round(ceiling, 4),
                    "contact_veto": a["contact_veto"],
                    "gate": a["gate"],
                    "refold_verdict": a["refold"],
                    "refold_tagged_sha256": a["tagged_sha256"],
                    "flag": a["flag"],
                    **evidence(sha, a["insert_after_residue"], r.gene),
                    "reason": note,
                    "agent_review": "",
                    "agent_reason": "",
                    "decided_by": DECIDED_BY,
                    "implemented_by": "ranker-pu-segment session (Claude Code)",
                    "decided_on": "2026-11-15",
                }
            )
        if not rows and r.decision in ("hold_no_passing_arm", "snorkel"):
            arms.append(
                {
                    "arm_id": f"{did}-a1",
                    "construct_key": "",
                    "design_id": did,
                    "member_id": r.member_id,
                    "gene": r.gene,
                    "sequence_sha256": sha,
                    "arm_role": "",
                    "arm_kind": "",
                    "decision": "hold",
                    "order_verdict": "",
                    "model_version": MODEL_VERSION,
                    "reason": note,
                    "agent_review": "",
                    "agent_reason": "",
                    "decided_by": DECIDED_BY,
                    "implemented_by": "ranker-pu-segment session (Claude Code)",
                    "decided_on": "2026-11-15",
                }
            )
        members.append(
            {
                "design_id": did,
                "member_id": r.member_id,
                "gene": r.gene,
                "decision": r.decision,
                "best_region_ceiling": round(ceiling, 4) if ceiling == ceiling else None,
                "recommended_sites_to_order": n_rec,
                "why": None if top is None else top["why"],
                "all_bad": all_bad,
                "arms_built": len(rows),
                "has_snorkel": any(a["arm_type"] == "snorkel" for a in rows),
                "has_flagged_loop": any(a["flag"] for a in rows),
                "status": status if rows or status != "built" else "held",
                "note": note,
            }
        )
    rescue_arms(arms, members)
    arms = apply_revisions(arms, members, prot, seen)
    arms = apply_annotation_rule(arms, prot)
    arms = apply_agent_adoptions(arms, members)
    arms = apply_lit_hard(arms)
    arms = assign_arm_ids(arms)
    a, m = pd.DataFrame(arms), pd.DataFrame(members)
    for did, why in BUILT_ELSEWHERE.items():
        a.loc[a["design_id"] == did, ["decision", "reason"]] = ["hold", why]
        m.loc[m["design_id"] == did, ["status", "note"]] = ["not_internal", why]
    a, unfolded = refold_of_record(a)
    stale = HERE / "l1_internal_arms_unfolded.tsv"
    stale.unlink(missing_ok=True)  # written only while some built arm has no verdict
    if len(unfolded):
        unfolded.to_csv(stale, sep="\t", index=False)
        print(f"{len(unfolded)} built arms have no refold verdict -> l1_internal_arms_unfolded.tsv")
    # the member summary follows the final arms (after revisions and the coordinator exceptions)
    b = a[a["decision"] == "build"]
    m["arms_built"] = m["design_id"].map(b.groupby("design_id").size()).fillna(0).astype(int)
    m["has_snorkel"] = m["design_id"].isin(b.loc[b["arm_kind"] == "snorkel", "design_id"])
    m["has_flagged_loop"] = m["design_id"].isin(
        b.loc[b["flag"].fillna("").astype(str) != "", "design_id"]
    )
    m.loc[(m["arms_built"] == 0) & (m["status"] == "built"), "status"] = "held"
    a = scores_and_veto_class(a, dec)
    # contact evidence at each arm's final site: arms moved by a revision, the annotation rule or
    # an agent adoption carried their old site's values until here
    ev = [evidence(x.sequence_sha256, int(x.insert_after_residue)
                   if x.arm_kind == "loop" and x.insert_after_residue == x.insert_after_residue
                   else None, x.gene) for x in a.itertuples(index=False)]  # fmt: skip
    ev = pd.DataFrame(ev, index=a.index)
    # final literature table (2026-11-15, lit_final_check.py): a contact that newly lands in an
    # arm's window is a flag only; no replacement in this pass (coordinator, 2026-11-15)
    lit = HERE / "lit_final_arms.tsv"
    if lit.is_file():
        lf = pd.read_csv(lit, sep="\t")
        lf = lf[lf["new_tier"].fillna("none") != "none"]
        note = {(r.design_id, int(r.insert_after_residue)):
                f"literature 2026-11-15 (final table, new): {r.new_contacts}; flag only, to "
                "resolve before the order" for r in lf.itertuples()}  # fmt: skip
        for i, x in a.iterrows():
            n = x["insert_after_residue"]
            k = (x["design_id"], int(n)) if x["decision"] == "build" and n == n else None
            if k in note:
                f = x["flag"] if isinstance(x["flag"], str) else ""
                a.at[i, "flag"] = "; ".join(filter(None, [f, note[k]]))
    a = a.drop(columns=[c for c in a.columns if c.startswith("nearest_")], errors="ignore")
    at = a.columns.get_loc("flag") + 1
    for c in ev.columns:
        if c.startswith("nearest_"):
            a.insert(at, c, ev[c])
            at += 1
    a.to_csv(OUT, sep="\t", index=False)
    m.to_csv(OUT_MEMBERS, sep="\t", index=False)  # one row per member
    built = m[m["status"] == "built"]
    summary = {
        "members": len(m),
        "arms_by_member": built["arms_built"].value_counts().sort_index().to_dict(),
        "members_with_snorkel": int(built["has_snorkel"].sum()),
        "members_with_flagged_loop": int(built["has_flagged_loop"].sum()),
        "held": int((m["status"] == "held").sum()),
        "internal_inserts": int((a["decision"] == "build").sum()),
        "snorkel_members": {
            k: int(v)
            for k, v in m[m["decision"] == "snorkel"]
            .apply(
                lambda x: (
                    "held"
                    if x["status"] == "held"
                    else (
                        "snorkel + flagged loop"
                        if x["has_flagged_loop"]
                        else "snorkel + passing loop"
                    )
                ),
                axis=1,
            )
            .value_counts()
            .items()
        },
        "q25_bar": round(q25, 4),
    }
    summary["arms_by_member"] = {str(k): int(v) for k, v in summary["arms_by_member"].items()}
    OUT_JSON.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
