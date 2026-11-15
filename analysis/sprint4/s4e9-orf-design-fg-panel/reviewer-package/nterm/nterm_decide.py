"""Decide every L1 N-after-SP member's N-terminal arm, on the evidence, under Becca's delegation.

Becca, 2026-11-12: "try to not queue protein specific things for management decisions or
anything that would block you from fully completing things. use your best judgment to complete
a new and improved evidence based l1 completely by the morning." Coordinator plan
``OVERNIGHT-L1-2026-11-12.md``, "RULE CHANGE". So nothing is held: every member gets one of

* ``native_site``: the p003 tandem after residue ``n`` of the native protein;
* ``igk_swap``: the Igκ leader plus Asp replace residues ``1..c``, and the tandem goes after
  the same mature residue;
* ``no_nterm_arm``: no tagged form keeps a secretory signal. No N-terminal arm is built, and
  the member is flagged to the manifest writer.

Two subcommands:

``inputs``
    Round-2 SignalP sequences, for every member with no passing form yet. Only sequences not
    already in a pinned run are written.

    * Members with no SignalP or disorder row anywhere (53 after PR #283): the native protein,
      and the tandem at the interval's own fallback span (``admissible_interval``,
      ``FALLBACK_WIDTH``; no track means no onset). Both native and Igκ forms are written.
    * Ortholog-release members called OTHER, predictor-disagreement members with no tracks,
      and orthologs or new members whose native forms all fail: Igκ forms on the same span or
      candidates.
    * ATRN, EPHA6, NRG2: routed here from the internal rung by ranker-pu-segment (PR #245,
      delegated decision, 2026-11-12). Native calls and Igκ forms at residues 15 to 45, because
      their cleavage is unknown until SignalP runs. The span used is chosen after the run.

``decide``
    One row per member (``nterm_decisions.tsv``), each with the arm, site, evidence, the
    rejected alternative, ``decided_by``, ``implemented_by`` and ``reason``.

Decision rules, in order of the evidence available:

1. **Preserved before the fold:** the re-prediction-gated site. Becca's rule, 2026-11-11.
2. **Preserved only inside the fold:** native site. Among the preserved candidates, the
   choice runs: cleavage-site probability softened by no more than 0.07, then the highest
   mean within-protein disorder percentile at the junction (metapredict, NetSurfP-3.0,
   AlphaFold-disorder = 1 - pLDDT), then the highest tagged cleavage probability. No refold
   verdict existed when this ran; ``refold_verdict`` is ``pending`` and a refold can
   overturn the call. Rejected: the Igκ swap, which does not move the cassette out of the fold.
3. **Cleavage moves but stays upstream of the cassette:** accept the shift
   (``nterm_site_picks.moved_pick``). The cassette is still at the mature N-terminus.
4. **Native signal fails at every candidate:** Igκ swap where an Igκ-tagged form keeps the
   leader's cleavage, at the preserved site chosen as in rule 2.
5. **Refused, ortholog release admitted:** gated native site; failing that, the Igκ swap.
6. **No passing form at all:** ``no_nterm_arm``.

    uv run --locked python nterm_decide.py inputs
    uv run --locked python nterm_decide.py decide
"""

from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
S4E9 = HERE.parents[1]
sys.path.insert(0, str(S4E9))
sys.path.insert(0, str(HERE))

import admissible_interval as ai  # noqa: E402
import nterm_build_options as nbo  # noqa: E402
import nterm_overnight as no  # noqa: E402
import nterm_placement as npl  # noqa: E402
import propeptide_check as ppc  # noqa: E402
from build_gpcr_rescue_arms import LEADER_AA, PLUS1_AA  # noqa: E402

from deliverome_analysis import tag_placement  # noqa: E402

INPUTS_TSV = HERE / "round2_inputs.tsv"
NEW_FASTA = HERE / "round2_inputs_new.fasta"
DECISIONS = HERE / "nterm_decisions.tsv"
SUMMARY = HERE / "nterm_decisions_summary.json"
#: Routed here from the internal rung by ranker-pu-segment (PR #245, delegated 2026-11-12);
#: a member already on the N-after-SP rung by the current table is skipped (EPHA6 after #283).
ROUTED = {"ATRN": "l1-p0163", "EPHA6": "l1-p0578", "NRG2": "l1-p1321"}  # design ids
ROUTED_SPAN = range(15, 46)
DECIDED_BY = 'nterm-reviewer-package (delegated by Becca 2026-11-12, "use your best judgment")'
BECCA_RULE = "Becca, 2026-11-11 (re-prediction-gated site rule)"
IMPLEMENTED_BY = "Claude Code session nterm-reviewer-package (PR #254)"
HEAD = len(LEADER_AA) + len(PLUS1_AA)


def sha(s: str) -> str:
    return hashlib.sha256(s.encode()).hexdigest()


def _tsv(path: Path) -> list[dict]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def fallback_candidates(c: int) -> list[tuple[int, bool]]:
    """The interval's fallback span when no disorder track exists: c+1 .. c+FALLBACK+2."""
    interval = ai.tag_interval(c, {})
    strict = ai.tag_interval(c, {}, tolerance=0)
    return [(n, n <= strict.hi) for n in ai.candidate_positions(interval, c)]


def all_members(processed) -> dict[str, dict]:
    """The N-after-SP lane plus the three members the ranker routed here."""
    members, _ = npl.load_members(processed)
    import nterm_placement as _np  # noqa: F401  (already imported; keeps load_members cached)

    from deliverome_analysis import library_design as ld

    proteins = ld.designed_proteins(S4E9 / "candidate_manifest_tranches.tsv", library=ld.LIBRARY_L1)
    cds = ld.designed_cds(S4E9 / "candidate_manifest_tranches.tsv", library=ld.LIBRARY_L1)
    design = ld.design_ids(S4E9 / "candidate_manifest_tranches.tsv")
    by_design = {d: c for c, d in design.items() if c in proteins}
    for gene, did in ROUTED.items():
        cid = by_design[did]
        if cid in members:  # already on the rung by the current table (EPHA6 after #283)
            continue
        p = proteins[cid]
        members[cid] = {
            "candidate_id": cid,
            "design_id": did,
            "gene_symbol": gene,
            "species": p["species"],
            "candidate_set": p["candidate_set"],
            "sequence": p["sequence"],
            "sha": sha(p["sequence"]),
            "cds": cds[cid],
            "v1_cleavage": None,
            "routed": True,
        }
    return members


def inputs() -> int:
    from deliverome_analysis.s3 import processed

    members = all_members(processed)
    opts = {r["candidate_id"]: r for r in _tsv(HERE / "nterm_build_options.tsv")}
    props = {r["candidate_id"]: r for r in _tsv(HERE / "overnight_proposals.tsv")}
    notes = json.loads((HERE / "overnight_input_notes.json").read_text())["notes"]
    note = {n["candidate_id"]: n for n in notes}
    over = collections.defaultdict(list)
    for r in _tsv(no.INPUTS_TSV):
        if r["arm"] in ("tagged", "igk_tagged"):
            over[r["candidate_id"]].append((int(r["native_after_residue"]), r["strict"] == "1"))
    rows = []
    for cid, m in sorted(members.items()):
        if m.get("routed"):
            rows += no.native_rows(cid, m, 0, [(n, False) for n in ROUTED_SPAN], "routed")
            continue
        sub = opts[cid]["subgroup"]
        p = props.get(cid)
        n_ = note.get(cid, {}).get("note", "")
        c = m["v1_cleavage"]
        if p and p["proposal"] == "hold" and p["group"] in ("new_member", "ortholog"):
            c = int(p["parent_cleavage"] or opts[cid]["cleavage_site"] or c)
            rows += no.igk_rows(cid, m, c, sorted(over[cid]), "igk_after_native_fail")
        elif sub == "refused_no_signalp_row" and n_ == "no ortholog row":
            if c is None:
                continue
            cands = fallback_candidates(c)
            rows += no.native_rows(cid, m, c, cands, "no_row_native")
            rows += no.igk_rows(cid, m, c, cands, "no_row_igk")
        elif (sub == "refused_no_signalp_row" and n_.startswith("refused_sp_disagree")) or (
            sub == "refused_sp_disagree_deeptmhmm_v1_only" and n_ == "no tracks"
        ):
            rows += no.igk_rows(cid, m, c, fallback_candidates(c), "igk_fallback")
    for r in rows:
        r["design_id"] = members[r["candidate_id"]]["design_id"]
        r["gene_symbol"] = members[r["candidate_id"]]["gene_symbol"]
        r["sequence_sha256"] = sha(r["sequence"])
    # the routed members also get Igκ forms over the same span, cleaving after their SP end,
    # which is unknown: write them after the native run decides c (second pass of `inputs`).
    known = no.prior_calls(processed)
    rows += routed_igk(members, known)
    for r in rows:
        r.setdefault("sequence_sha256", sha(r["sequence"]))
        r.setdefault("design_id", members[r["candidate_id"]]["design_id"])
        r.setdefault("gene_symbol", members[r["candidate_id"]]["gene_symbol"])
    cols = [
        "group",
        "candidate_id",
        "design_id",
        "gene_symbol",
        "arm",
        "after_residue",
        "native_after_residue",
        "cleavage_site",
        "strict",
        "sequence_sha256",
    ]
    with INPUTS_TSV.open("w", newline="") as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=cols,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
            restval="",
        )
        w.writeheader()
        w.writerows(rows)
    seen = set()
    with NEW_FASTA.open("w") as fh:
        for r in rows:
            if r["sequence_sha256"] in seen or r["sequence_sha256"] in known:
                continue
            seen.add(r["sequence_sha256"])
            fh.write(f">{r['sequence_sha256']}\n{r['sequence']}\n")
    print(collections.Counter((r["group"], r["arm"]) for r in rows))
    print(f"{len(rows)} rows; {len(seen)} not yet predicted -> {NEW_FASTA.name}")
    return 0


def routed_igk(members: dict, known: dict) -> list[dict]:
    """Igκ forms for a routed member, once its native SignalP call gives the cleavage site."""
    out = []
    for cid, m in members.items():
        if not m.get("routed"):
            continue
        call = known.get(m["sha"])
        if not call or call["prediction"] != "SP":
            continue
        c = call["cleavage_site"]
        out += no.igk_rows(cid, m, c, fallback_candidates(c), "routed_igk")
    return out


# ---- decide ----


def junction_disorder(m: dict, base: dict, n: int) -> float | None:
    """Mean within-protein disorder percentile at residue n over the member's tracks."""
    _, tracks = npl.evidence(npl.route(m, None, base), base)
    vals = [
        tag_placement.within_protein_percentile(t)[n - 1] for t in tracks.values() if len(t) >= n
    ]
    return sum(vals) / len(vals) if vals else None


def decide() -> int:
    from deliverome_analysis.s3 import processed

    members = all_members(processed)
    base = npl.load_snapshot(processed)
    calls = no.prior_calls(processed)
    picks = {r["candidate_id"]: r for r in _tsv(HERE / "nterm_site_picks.tsv")}
    props = {r["candidate_id"]: r for r in _tsv(HERE / "overnight_proposals.tsv")}
    positions = collections.defaultdict(list)
    for p in _tsv(HERE / "nterm_positions.tsv"):
        positions[p["candidate_id"]].append(p)
    r2 = collections.defaultdict(list)
    for r in _tsv(INPUTS_TSV):
        r2[r["candidate_id"]].append(r)
    # The native residue c the Igκ leader replaces (1..c), from the rows that built the forms.
    replaced = {}
    for path in (no.INPUTS_TSV, INPUTS_TSV):
        for r in _tsv(path):
            if r["arm"] == "igk_native":
                replaced[r["candidate_id"]] = r["cleavage_site"]

    def row(
        cid,
        arm,
        n,
        native_n,
        c,
        inside,
        cp,
        tagged_sha,
        decided_by,
        reason,
        alternative,
        evidence,
        refold="pending",
    ):
        m = members[cid]
        return {
            "candidate_id": cid,
            "design_id": m["design_id"],
            "gene_symbol": m["gene_symbol"],
            "candidate_set": m["candidate_set"],
            "arm": arm,
            "insert_after_residue": n,
            "igk_replaces_residues": f"1-{replaced[cid]}"
            if arm == "igk_swap" and cid in replaced
            else "",
            "native_after_residue": native_n,
            "cleavage_site": c,
            "inside_predicted_fold": inside,
            "tagged_cleavage_probability": cp,
            "tagged_sha256": tagged_sha,
            "evidence": evidence,
            "refold_verdict": refold,
            "decided_by": decided_by,
            "implemented_by": IMPLEMENTED_BY,
            "reason": reason,
            "rejected_alternative": alternative,
        }

    def best_round2(cid: str, parent_arm: str, tagged_arm: str):
        rows = [r for r in r2[cid] if r["arm"] in (parent_arm, tagged_arm)]
        par = next((r for r in rows if r["arm"] == parent_arm), None)
        if not par or calls[par["sequence_sha256"]]["prediction"] != "SP":
            return None, None
        pc = calls[par["sequence_sha256"]]
        ok = []
        for t in rows:
            if t["arm"] != tagged_arm:
                continue
            call = calls[t["sequence_sha256"]]
            v = npl.verdict(pc, call, int(t["after_residue"]))
            # The cassette must follow the cleavage site of the parent's own SignalP call: these
            # candidates were laid out from DeepTMHMM v1's cleavage alone, and where SignalP
            # cuts later (SMO ortholog 34 vs 32, ADAM33 isoform 29 vs 27) the first candidates
            # sit at or before the cut. The rule of record takes the later call and excludes it.
            if v == "preserved" and int(t["after_residue"]) > pc["cleavage_site"]:
                key = nbo.pick_key(
                    t["strict"] == "1",
                    pc["cleavage_probability"] or 0,
                    call["cleavage_probability"] or 0,
                    int(t["native_after_residue"] or t["after_residue"]),
                )
                ok.append((key, t, call))
        if not ok:
            return None, pc
        best = min(ok, key=lambda x: x[0])
        return (None, None, None, best[1], best[2]), pc

    out = []
    for cid, m in sorted(members.items()):
        pk, pr = picks.get(cid), props.get(cid)
        status = pk["status"] if pk else "routed"
        if status == "ready":
            out.append(
                row(
                    cid,
                    "native_site",
                    pk["insert_after_residue"],
                    pk["insert_after_residue"],
                    pk["cleavage_site"],
                    0,
                    pk["tagged_cleavage_probability"],
                    pk["tagged_protein_sha256"],
                    BECCA_RULE,
                    "re-prediction-gated site before the predicted fold",
                    "none needed",
                    "SignalP keeps the native cleavage; site before fold",
                    refold="pending",
                )
            )
            continue
        if status == "held_inside_fold" and not any(
            p["verdict"] == "preserved" for p in positions[cid]
        ):
            # re-predicted after the sweep (nterm_overnight new_member): its gated pick stands
            out.append(
                row(
                    cid,
                    "native_site",
                    pk["insert_after_residue"],
                    pk["insert_after_residue"],
                    pk["cleavage_site"],
                    1,
                    pk["tagged_cleavage_probability"],
                    pk["tagged_protein_sha256"],
                    DECIDED_BY,
                    "only preserved sites are inside the predicted fold "
                    "(re-predicted after the sweep); gated pick",
                    "Igκ swap (same fold problem) or hold for refold",
                    "overnight re-prediction",
                )
            )
            continue
        if status == "held_inside_fold":
            native_cp = float(pk["native_cleavage_probability"] or 0)
            ok = [p for p in positions[cid] if p["verdict"] == "preserved"]
            scored = []
            for p in ok:
                n = int(p["after_residue"])
                d = junction_disorder(m, base, n)
                soft = native_cp - float(p["tagged_cleavage_probability"]) > nbo.SOFTENING_BOUND
                scored.append((soft, -(d or 0), -float(p["tagged_cleavage_probability"]), n, p, d))
            soft, _nd, _, n, p, d = min(scored, key=lambda x: x[:4])
            out.append(
                row(
                    cid,
                    "native_site",
                    n,
                    n,
                    pk["cleavage_site"],
                    1,
                    p["tagged_cleavage_probability"],
                    p["tagged_sha256"],
                    DECIDED_BY,
                    "only preserved sites are inside the predicted fold; chosen by "
                    "unsoftened cleavage, then junction disorder percentile "
                    f"({d:.0f})"
                    if d is not None
                    else "no disorder track",
                    "Igκ swap (same fold problem) or hold for refold",
                    f"{len(ok)} preserved sites; softened={int(soft)}",
                )
            )
            continue
        if status == "held_moved_upstream":
            out.append(
                row(
                    cid,
                    "native_site",
                    pk["insert_after_residue"],
                    pk["insert_after_residue"],
                    pk["cleavage_site"],
                    pk["inside_predicted_fold"],
                    pk["tagged_cleavage_probability"],
                    pk["tagged_protein_sha256"],
                    DECIDED_BY,
                    f"cleavage shifts to {pk['tagged_cleavage_site']} but stays upstream "
                    "of the cassette, which stays at the mature N-terminus",
                    "Igκ swap, which changes the mature N-terminus too",
                    "SignalP: shifted, not lost",
                )
            )
            continue
        if pr and pr["proposal"] != "hold":
            arm = pr["proposal"]
            rule = (
                "native signal peptide fails at every site; Igκ-tagged form keeps the "
                "leader's cleavage"
                if arm == "igk_swap"
                else "ortholog SignalP release admitted; gated site"
                if pr["group"] == "ortholog"
                else "member new to the rung, re-predicted; gated site"
            )
            out.append(
                row(
                    cid,
                    arm,
                    pr["insert_after_residue"],
                    pr["native_after_residue"],
                    pr["parent_cleavage"],
                    pr["inside_predicted_fold"],
                    pr["tagged_cleavage_probability"],
                    pr["tagged_sha256"],
                    DECIDED_BY,
                    rule,
                    "no N-terminal arm" if arm == "igk_swap" else "Igκ swap",
                    f"{pr['n_preserved']}/{pr['n_tagged']} tagged forms preserved",
                )
            )
            continue
        # round 2
        decided = False
        for parent_arm, tagged_arm, arm in (
            ("native", "tagged", "native_site"),
            ("igk_native", "igk_tagged", "igk_swap"),
        ):
            best, pc = best_round2(cid, parent_arm, tagged_arm)
            if best:
                t, call = best[3], best[4]
                out.append(
                    row(
                        cid,
                        arm,
                        t["after_residue"],
                        t["native_after_residue"],
                        pc["cleavage_site"],
                        int(t["strict"] != "1"),
                        call["cleavage_probability"],
                        t["sequence_sha256"],
                        DECIDED_BY,
                        f"round-2 SignalP ({t['group']}): tagged form keeps the cleavage",
                        "Igκ swap" if arm == "native_site" else "no N-terminal arm",
                        "interval from the fallback span (no disorder track)"
                        if t["group"].startswith(("no_row", "igk_fallback", "routed"))
                        else "candidates as before",
                    )
                )
                decided = True
                break
        if not decided:
            out.append(
                row(
                    cid,
                    "no_nterm_arm",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                    DECIDED_BY,
                    "no tagged form, native or Igκ, keeps a SignalP cleavage site",
                    "build anyway at the gated or v1 site",
                    (pr or {}).get("reason", "") or "no passing form",
                    refold="",
                )
            )
    out, requests = apply_refold(out, members, base, calls, positions, r2, replaced, row)
    out, more, need = apply_propeptides(out, members, calls, row)
    requests += more
    out = apply_weak_sp(out, row)
    out, more_sp = apply_sp_disagree(out, members, calls, row)
    requests += more_sp
    asked = {r["tagged_sha256"] for r in requests}
    for d in out:  # every built form still without a verdict goes to refold
        if (
            d["arm"] in ("native_site", "igk_swap")
            and d["refold_verdict"].startswith("pending")
            and d["tagged_sha256"] not in asked
        ):
            asked.add(d["tagged_sha256"])
            requests.append(
                {
                    "candidate_id": d["candidate_id"],
                    "design_id": d["design_id"],
                    "arm": d["arm"],
                    "insert_after_residue": d["insert_after_residue"],
                    "tagged_sha256": d["tagged_sha256"],
                    "replaces": "",
                }
            )
    out = finalize(out, row)
    out = apply_ortholog_rule(out, row)
    out = [
        d | {"arm": "routed_multipass", "reason": RETIRED[d["design_id"]]}
        if d["design_id"] in RETIRED and d["arm"] == "no_nterm_arm"
        else d
        for d in out
    ]
    if not requests:  # nothing left to fold: leave an empty table, never a stale one
        REFOLD_REQUESTS.write_text(
            "candidate_id\tdesign_id\tarm\tinsert_after_residue\ttagged_sha256\treplaces\n"
        )
    if requests:
        with REFOLD_REQUESTS.open("w", newline="") as fh:
            w = csv.DictWriter(
                fh, fieldnames=list(requests[0]), delimiter="\t", lineterminator="\n"
            )
            w.writeheader()
            w.writerows(requests)
    if not ALL_ALTERNATIVES:
        REFOLD_ALTERNATIVES.write_text(
            "candidate_id\tdesign_id\tarm\tinsert_after_residue\ttagged_sha256\treplaces\n"
        )
    if ALL_ALTERNATIVES:
        with REFOLD_ALTERNATIVES.open("w", newline="") as fh:
            w = csv.DictWriter(
                fh, fieldnames=list(ALL_ALTERNATIVES[0]), delimiter="\t", lineterminator="\n"
            )
            w.writeheader()
            w.writerows(ALL_ALTERNATIVES)
    print(f"{len(ALL_ALTERNATIVES)} further untested alternatives -> {REFOLD_ALTERNATIVES.name}")
    PROPEP_FASTA.write_text("".join(f">{k}\n{v}\n" for k, v in need.items()))
    print(f"{len(need)} relocated forms not yet called by SignalP -> {PROPEP_FASTA.name}")
    # Keep the existing table's row order (by design id); new rows go last. The order carries
    # no meaning, and keeping it makes each change's diff only its own rows.
    if DECISIONS.exists():
        prior = {r["design_id"]: i for i, r in enumerate(_tsv(DECISIONS))}
        out.sort(key=lambda r: (prior.get(r["design_id"], len(prior)), r["design_id"]))
    with DECISIONS.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out)
    summary = {
        "members": len(out),
        "members_on_rung": sum(not r["arm"].startswith("routed_") for r in out),
        "arms": dict(sorted(collections.Counter(r["arm"] for r in out).items())),
        "inside_predicted_fold": dict(
            sorted(
                collections.Counter(f"{r['arm']}:{r['inside_predicted_fold']}" for r in out).items()
            )
        ),
        "decided_by": dict(collections.Counter(r["decided_by"] for r in out)),
        "routed_internal": sorted(
            f"{r['candidate_id']} {r['gene_symbol']} {r['arm']}"
            for r in out
            if r["arm"].startswith("routed_")
        ),
        "no_nterm_arm": sorted(
            f"{r['candidate_id']} {r['gene_symbol']}" for r in out if r["arm"] == "no_nterm_arm"
        ),
        "calls_from_runs": list(no.RUN_IDS),
        "refold_verdicts": dict(
            sorted(collections.Counter(r["refold_verdict"].split(":")[0] for r in out).items())
        ),
        "refold_requests": len(requests),
        "relocated_past_propeptide": sum("pro-peptide" in r["evidence"] for r in out),
        "refold_verdicts_sha256": VERDICTS_USED.get("sha256", ""),
        "refold_verdicts_source": VERDICTS_KEY,
        "hold_pending": HOLD_PENDING["on"],
        "uniprot_snapshot": json.loads((HERE / "uniprot_features_snapshot.json").read_text())[
            "release"
        ],
    }
    SUMMARY.write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps({k: v for k, v in summary.items() if k != "no_nterm_arm"}, indent=1))
    print(f"no_nterm_arm: {len(summary['no_nterm_arm'])}")
    return 0


# ---- refold verdicts ----

VERDICTS_KEY = (
    "sprint4/s4e9-orf-design-fg-panel/predicted-structures/tagged/2026-11-12/verdicts_nterm.tsv"
)
REFOLD_REQUESTS = HERE / "refold_requests.tsv"
#: SHA-256 of the verdicts file this run read; ``decision_inputs.py`` freezes it by that hash.
VERDICTS_USED: dict[str, str] = {}


def load_verdicts(processed) -> dict[str, dict]:
    """l1-tag-refold's verdicts keyed by tagged-protein SHA-256 (empty until it publishes)."""
    if Path(VERDICTS_KEY).exists():  # a committed table, e.g. main's tag_refold_verdicts.tsv
        path = VERDICTS_KEY
    elif not processed.exists(VERDICTS_KEY):
        return {}
    else:
        path = processed.download(VERDICTS_KEY)
    VERDICTS_USED["sha256"] = hashlib.sha256(Path(path).read_bytes()).hexdigest()
    with open(path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    key = next(k for k in ("protein_sha256", "tagged_sha256") if k in rows[0])
    col = next(k for k in ("verdict", "pass_fail", "result", "pass") if k in rows[0])
    out = {}
    # A form can be scored under several rungs; the N-terminal row wins.
    rows.sort(key=lambda r: str(r.get("rung", "")).startswith("nterm"))
    for r in rows:
        v = str(r[col]).strip().lower()
        if v not in ("pass", "fail"):  # not_folded, no_native_model: no verdict yet
            continue
        why = r.get("verdict_reasons") or r.get("fail_reasons") or ""
        out[r[key]] = {
            "row": r,
            "pass": v == "pass" or borderline(why),
            "why": f"borderline, counted as pass ({why})"
            if v == "fail" and borderline(why)
            else why,
        }
    return out


#: A fail whose only failing criteria are lDDT values this close to their bar is a pass.
#: Decided by nterm-reviewer-package (delegated by Becca 2026-11-12): the verdict rests on one
#: ColabFold model (model_1_ptm, 3 recycles), whose lDDT moves by more than 0.01 between runs,
#: so FOLR1 at 0.698 against 0.70 is not evidence against the site. Epitope packing, a buried
#: tag and a helix inside the cassette stay decisive at any margin.
LDDT_MARGIN = 0.01
_LDDT = re.compile(r"lDDT[^;]*?([0-9.]+) < ([0-9.]+)")


def borderline(why: str) -> bool:
    parts = [p.strip() for p in why.split(";") if p.strip()]
    if not parts:
        return False
    for p in parts:
        m = _LDDT.search(p)
        if not m or float(m.group(2)) - float(m.group(1)) > LDDT_MARGIN + 1e-9:
            return False
    return True


def alternatives(cid, m, base, calls, positions, r2, replaced) -> list[dict]:
    """Every SignalP-passing tagged form of a member, best first: native, then Igκ."""
    opts = []
    for p in positions.get(cid, []):
        if p["verdict"] == "preserved":
            n = int(p["after_residue"])
            opts.append(
                {
                    "arm": "native_site",
                    "n": n,
                    "native_n": n,
                    "sha": p["tagged_sha256"],
                    "strict": p["within_strict_interval"] == "1",
                    "cp": float(p["tagged_cleavage_probability"] or 0),
                }
            )
    for path in (no.INPUTS_TSV, INPUTS_TSV):
        rows = [r for r in _tsv(path) if r["candidate_id"] == cid]
        parents = {
            r["arm"]: calls.get(r["sequence_sha256"])
            for r in rows
            if r["arm"] in ("native", "igk_native")
        }
        for t in rows:
            par = parents.get("native" if t["arm"] == "tagged" else "igk_native")
            if t["arm"] not in ("tagged", "igk_tagged") or not par or par["prediction"] != "SP":
                continue
            call = calls.get(t["sequence_sha256"])
            if (
                call
                and npl.verdict(par, call, int(t["after_residue"])) == "preserved"
                and int(t["after_residue"]) > par["cleavage_site"]
            ):
                opts.append(
                    {
                        "arm": "native_site" if t["arm"] == "tagged" else "igk_swap",
                        "n": int(t["after_residue"]),
                        "native_n": int(t["native_after_residue"] or t["after_residue"]),
                        "sha": t["sequence_sha256"],
                        "strict": t["strict"] == "1",
                        "cp": call["cleavage_probability"] or 0,
                    }
                )
    seen, ranked = set(), []
    for o in sorted(
        opts,
        key=lambda o: (
            o["arm"] != "native_site",
            not o["strict"],
            -(junction_disorder(m, base, o["native_n"]) or 0),
            -o["cp"],
            o["native_n"],
        ),
    ):
        if o["sha"] not in seen:
            seen.add(o["sha"])
            ranked.append(o)
    return ranked


def apply_refold(out, members, base, calls, positions, r2, replaced, row):
    """Record refold verdicts; replace a failed arm with its next untested alternative.

    Decided by nterm-reviewer-package (delegated by Becca 2026-11-12): a failed refold (native
    fold lost or cassette buried, l1-tag-refold's criteria) overturns the arm. The next
    SignalP-passing alternative not already failed takes its place, pending its own refold
    (``refold_requests.tsv``). A member with nothing left gets ``no_nterm_arm``.
    """
    from deliverome_analysis.s3 import processed

    verdicts = load_verdicts(processed)
    requests, result = [], []
    for d in out:
        v = verdicts.get(d["tagged_sha256"])
        if not v:
            result.append(d)
            continue
        if v["pass"]:
            result.append(
                d
                | {"refold_verdict": "pass" if not v["why"].startswith("borderline") else v["why"]}
            )
            continue
        cid = d["candidate_id"]
        failed = {s for s, x in verdicts.items() if not x["pass"]}
        untested = [
            o
            for o in alternatives(cid, members[cid], base, calls, positions, r2, replaced)
            if o["sha"] not in failed
        ]
        # an alternative already folded and passing comes first; the order is otherwise kept
        untested.sort(key=lambda o: not (verdicts.get(o["sha"]) or {}).get("pass"))
        alt = untested[0] if untested else None
        # Every remaining alternative goes to refold at once, so one round settles the member
        # instead of one residue per round (2026-11-12: 45 of 59 next-site arms failed again).
        for o in untested[1:]:
            if o["sha"] not in verdicts and (o["arm"] != "igk_swap" or cid in replaced):
                ALL_ALTERNATIVES.append(
                    {
                        "candidate_id": cid,
                        "design_id": members[cid]["design_id"],
                        "arm": o["arm"],
                        "insert_after_residue": o["n"],
                        "tagged_sha256": o["sha"],
                        "replaces": d["tagged_sha256"],
                    }
                )
        if alt is None:
            pick = disposition_pick(
                d,
                alternatives(cid, members[cid], base, calls, positions, r2, replaced),
                verdicts,
                replaced,
            )
            if pick:
                o, vrow, klass = pick
                new = row(
                    cid,
                    o["arm"],
                    o["n"],
                    o["native_n"],
                    d["cleavage_site"],
                    int(not o["strict"]),
                    o["cp"],
                    o["sha"],
                    DISPOSITION_BY,
                    f"refold failed at every SignalP-passing form; built and flagged under the "
                    f"multipass rung's disposition rule ({klass}): {DISPOSITION_TEXT[klass]}",
                    "hold the member (every form fails refold)",
                    "l1-tag-refold verdict + disposition rule",
                    refold=f"fail: {verdicts[o['sha']]['why']}",
                )
                result.append(new | {"refold_disposition": disposition_text(vrow, klass)})
                continue
            result.append(
                row(
                    cid,
                    "no_nterm_arm",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                    DECIDED_BY,
                    f"refold failed ({v['why']}) and no other SignalP-passing form",
                    f"keep the failed {d['arm']} after {d['insert_after_residue']}",
                    "l1-tag-refold verdict",
                    refold=f"fail: {v['why']}",
                )
            )
            continue
        if alt["arm"] == "igk_swap" and cid not in replaced:
            result.append(d | {"refold_verdict": f"fail: {v['why']}"})
            continue
        new = row(
            cid,
            alt["arm"],
            alt["n"],
            alt["native_n"],
            d["cleavage_site"],
            int(not alt["strict"]),
            alt["cp"],
            alt["sha"],
            DECIDED_BY,
            f"refold of {d['arm']} after {d['insert_after_residue']} failed "
            f"({v['why']}); next SignalP-passing form",
            f"the failed {d['arm']} after {d['insert_after_residue']}",
            "l1-tag-refold verdict + SignalP",
            refold=(
                "pending (replaces a failed arm)"
                if alt["sha"] not in verdicts
                else "pass"
                if not verdicts[alt["sha"]]["why"].startswith("borderline")
                else verdicts[alt["sha"]]["why"]
            ),
        )
        result.append(new)
        if alt["sha"] in verdicts:
            continue
        requests.append(
            {
                "candidate_id": cid,
                "design_id": members[cid]["design_id"],
                "arm": alt["arm"],
                "insert_after_residue": alt["n"],
                "tagged_sha256": alt["sha"],
                "replaces": d["tagged_sha256"],
            }
        )
    return result, requests


SP_REVIEW = HERE / "sp_disagree_review.tsv"
#: Re-anchor an Igκ swap at an independently confirmed cleavage site (UniProt SIGNAL, else
#: DeepTMHMM v2) when the DeepTMHMM v1 site is not confirmed. False holds those members instead.
REANCHOR = {"on": True}


def apply_sp_disagree(out, members, calls, row):
    """Apply ``sp_disagree_review.tsv`` to class ``refused_sp_disagree_deeptmhmm_v1_only``.

    Rule decided by the coordinator under Becca's 2026-11-12 delegation (2026-11-14): the Igκ
    swap deletes residues 1..c on DeepTMHMM v1's call alone, so it stays only where an
    independent source (UniProt SIGNAL on the same N-terminal sequence, or DeepTMHMM v2)
    confirms a cleaved signal peptide within 3 residues of c. A signal anchor or first TM helix
    is not deleted: the member is held here with the rung it fits. A cleaved signal peptide at
    another site re-anchors the swap there (implemented by nterm-reviewer-package, proposed to
    the coordinator 2026-11-14); a member with no independent evidence is held.
    """
    import build_constructs as bc

    review = {r["candidate_id"]: r for r in _tsv(SP_REVIEW)}
    for r in review.values():
        ev = r["cleaved_sp_evidence"] or (
            f"{r['reanchor_source']} [{r['reanchor_evidence']}]" if r["reanchor_cleavage"] else ""
        )
        if r.get("human_ortholog_signal") and "human ortholog" not in ev:
            ev = f"{ev}; {r['human_ortholog_signal']}" if ev else r["human_ortholog_signal"]
        SP_EVIDENCE[r["candidate_id"]] = ev
    census = {r["candidate_id"]: r for r in _tsv(HERE / "nterm_members.tsv")}
    for cid, m in members.items():
        c = census.get(cid, {})
        if c.get("signalp_prediction") == "SP":
            NATIVE_SP[cid] = (c["signalp_cleavage"], c["signalp_cleavage_probability"])
            continue
        call = calls.get(sha(m["sequence"]))
        if call and call["prediction"] == "SP":
            NATIVE_SP[cid] = (call["cleavage_site"], round(call["cleavage_probability"] or 0, 4))
    result, requests = [], []
    for d in out:
        r = review.get(d["candidate_id"])
        d = d | {"sp_review_verdict": "", "sp_sources_checked": ""}
        if not r:
            result.append(d)
            continue
        d |= {"sp_review_verdict": r["verdict"], "sp_sources_checked": r["sources_checked"]}
        c1 = r["deeptmhmm_v1_cleavage"]
        lead = f"DeepTMHMM v1 calls a signal peptide 1-{c1} that SignalP 6.0 does not"
        cid = d["candidate_id"]
        if r["verdict"] == "cleaved_sp_confirmed":
            if d["arm"] == "igk_swap":
                d["reason"] = (
                    f"{lead}; {r['cleaved_sp_evidence']} confirms a cleaved signal peptide, so "
                    "the Igκ leader replaces it"
                )
                d["decided_by"] = SP_REVIEW_BY
            result.append(d)
            continue
        if r["verdict"] == "signal_anchor_or_tm" and d["design_id"] in DEFERRED_ROUTES:
            d["reason"] = (
                f"{lead}; it is a signal anchor ({r['anchor_or_tm_evidence']}): "
                + DEFERRED_ROUTES[d["design_id"]]
            )
            d["decided_by"] = SP_REVIEW_BY
            result.append(d)
            continue
        if r["verdict"] == "signal_anchor_or_tm" or not (REANCHOR["on"] and r["reanchor_cleavage"]):
            why = (
                f"{lead}; it is a signal anchor or first TM helix ({r['anchor_or_tm_evidence']}), "
                f"so residues 1-{c1} stay: route to the "
                f"{r['route_to'] or ROUTE_OVERRIDE.get(d['design_id'], 'appropriate')} rung"
                if r["verdict"] == "signal_anchor_or_tm"
                else f"{lead}; "
                + (
                    r.get("hold_reason")
                    or "no independent source confirms a cleaved signal peptide"
                )
            )
            # Retired here and built on the rung that fits: OR11H6 and TMEM63B by the multipass
            # builder (PR #337), ENTPD2 by the internal-rung ranker. Those sessions own the lane
            # override; CATSPERH stays held (coordinator, 2026-11-14).
            arm = (
                ROUTED_ARM.get(
                    r["route_to"] or ROUTE_OVERRIDE.get(d["design_id"], ""), "no_nterm_arm"
                )
                if r["verdict"] == "signal_anchor_or_tm"
                else "no_nterm_arm"
            )
            result.append(
                row(
                    cid,
                    arm,
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                    SP_REVIEW_BY,
                    why,
                    f"Igκ swap deleting residues 1-{c1}",
                    "second-source review",
                    refold="",
                )
                | {k: d[k] for k in ("sp_review_verdict", "sp_sources_checked")}
            )
            continue
        c2 = int(r["reanchor_cleavage"])
        seq = members[cid]["sequence"]
        n = int(d["native_after_residue"] or 0)
        n = n if n > c2 else c2 + 1
        parent = LEADER_AA + PLUS1_AA + seq[c2:]
        q = HEAD + (n - c2)
        tagged = parent[:q] + bc.CASSETTE_AA + parent[q:]
        pc, tc = calls.get(sha(parent)), calls.get(sha(tagged))
        if not (pc and tc and npl.verdict(pc, tc, q) == npl.PRESERVED and q > pc["cleavage_site"]):
            result.append(
                row(
                    cid,
                    "no_nterm_arm",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                    SP_REVIEW_BY,
                    f"{lead}; re-anchored at {r['reanchor_source']} the Igκ form is not "
                    "confirmed by SignalP",
                    f"Igκ swap deleting residues 1-{c1}",
                    "second-source review + SignalP",
                    refold="",
                )
                | {k: d[k] for k in ("sp_review_verdict", "sp_sources_checked")}
            )
            continue
        replaced_by = {cid: c2}
        new = row(
            cid,
            "igk_swap",
            q,
            n,
            pc["cleavage_site"],
            d["inside_predicted_fold"],
            tc["cleavage_probability"],
            sha(tagged),
            SP_REVIEW_BY,
            f"{lead}; {r['reanchor_source']} [{r['reanchor_evidence']}] puts the cleaved "
            f"signal peptide elsewhere, so the Igκ leader replaces residues 1-{c2} instead",
            f"Igκ swap deleting residues 1-{c1} (DeepTMHMM v1 site)",
            "second-source review + SignalP",
            refold="pending",
        )
        new["igk_replaces_residues"] = f"1-{replaced_by[cid]}"
        new |= {k: d[k] for k in ("sp_review_verdict", "sp_sources_checked")}
        verdicts = load_verdicts_cached()
        v = verdicts.get(sha(tagged))
        if v:
            new["refold_verdict"] = (
                ("pass" if not v["why"].startswith("borderline") else v["why"])
                if v["pass"]
                else f"fail: {v['why']}"
            )
            if not v["pass"]:
                new = row(
                    cid,
                    "no_nterm_arm",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                    SP_REVIEW_BY,
                    f"{lead}; re-anchored at {r['reanchor_source']} the Igκ form fails "
                    f"refold ({v['why']})",
                    f"Igκ swap deleting residues 1-{c1}",
                    "second-source review + l1-tag-refold verdict",
                    refold=f"fail: {v['why']}",
                ) | {k: d[k] for k in ("sp_review_verdict", "sp_sources_checked")}
        else:
            requests.append(
                {
                    "candidate_id": cid,
                    "design_id": members[cid]["design_id"],
                    "arm": "igk_swap",
                    "insert_after_residue": q,
                    "tagged_sha256": sha(tagged),
                    "replaces": d["tagged_sha256"],
                }
            )
        result.append(new)
    return result, requests


#: Signal anchors whose Igκ arm stays until the receiving rung's arm is built and passes refold
#: (agreed with the ranker session, 2026-11-14). The only rows allowed an Igκ deletion with no
#: cleaved-signal-peptide evidence; ``tests/test_nterm_igk_evidence.py`` holds them to this list.
DEFERRED_ROUTES: dict[str, str] = {}
#: Signal anchors whose receiving rung is known from that rung's own build, where the topology
#: string gives no side: cyno ENTPD2 is built internal after 458 by the ranker (PR #341,
#: gate +21.2 Å, refold pass, 2026-11-14), as human ENTPD2 l1-p0567 is.
ROUTE_OVERRIDE = {"l1-p2987": "internal"}
ROUTED_ARM = {"residue-1 multipass": "routed_multipass", "internal": "routed_internal"}
SP_REVIEW_BY = "coordinator under Becca's 2026-11-12 delegation, 2026-11-14 (second-source rule)"
_VERDICTS_CACHE: dict = {}


def load_verdicts_cached() -> dict:
    if "v" not in _VERDICTS_CACHE:
        from deliverome_analysis.s3 import processed

        _VERDICTS_CACHE["v"] = load_verdicts(processed)
    return _VERDICTS_CACHE["v"]


HOLD_PENDING = {"on": False}
HOLD_PENDING_BY = (
    "coordinator decision under Becca's delegation, 2026-11-14 (finalize #254 for merge)"
)


def finalize(out, row):
    """Hold arms with no refold verdict (``--hold-pending``) and spell out the refold rule.

    Every row gets the rule its ``refold_verdict`` rests on, so a reader (or l1_consistency)
    sees where this rung's call differs from l1-tag-refold's by rule rather than by omission:
    ``refold_rule``, ``refold_rule_metric``, ``refold_rule_threshold``,
    ``refold_rule_decided_by``, ``refold_rule_date``.
    """
    result = []
    routed = set(ROUTED.values())
    for d in out:
        if d["design_id"] in routed and d["arm"] == "no_nterm_arm":
            # #245 proposed moving these onto this rung; no N-after-SP form keeps a SignalP
            # cleavage site, so they stay on the internal rung, where l1_tag_paths.tsv places
            # and builds them. They are not N-after-SP members and not holds.
            d = d | {
                "arm": "routed_internal",
                "reason": "routed back to internal (#245): no N-after-SP form, native or Igκ, "
                "keeps a SignalP cleavage site; built on the internal rung (l1_tag_paths.tsv)",
            }
        v = d["refold_verdict"]
        if HOLD_PENDING["on"] and d["arm"] != "no_nterm_arm" and v.startswith("pending"):
            replacing = "replaces" in v
            d = row(
                d["candidate_id"],
                "no_nterm_arm",
                "",
                "",
                "",
                "",
                "",
                "",
                HOLD_PENDING_BY,
                "replacement arm refold pending" if replacing else "arm refold pending",
                f"build the unfolded {d['arm']} after {d['insert_after_residue']}",
                d["evidence"],
                refold=v,
            )
        rule = {
            "refold_rule": "",
            "refold_rule_metric": "",
            "refold_rule_threshold": "",
            "refold_rule_decided_by": "",
            "refold_rule_date": "",
        }
        v = d["refold_verdict"]
        if v.startswith("borderline"):
            m = _LDDT.search(v)
            rule = {
                "refold_rule": "lddt_borderline_counted_as_pass",
                "refold_rule_metric": v[v.index("(") + 1 : v.rindex("<")].strip(),
                "refold_rule_threshold": f"{m.group(2)} with margin {LDDT_MARGIN}",
                "refold_rule_decided_by": DECIDED_BY,
                "refold_rule_date": "2026-11-13",
            }
        elif d["decided_by"] == HOLD_PENDING_BY:
            rule = {
                "refold_rule": "hold_without_refold_verdict",
                "refold_rule_metric": "",
                "refold_rule_threshold": "",
                "refold_rule_decided_by": HOLD_PENDING_BY,
                "refold_rule_date": "2026-11-14",
            }
        elif v.startswith(("pass", "fail")):
            rule = {
                "refold_rule": "l1-tag-refold verdict",
                "refold_rule_metric": "",
                "refold_rule_threshold": "",
                "refold_rule_decided_by": "l1-tag-refold (PR #301)",
                "refold_rule_date": "2026-11-12",
            }
        if d.get("refold_disposition"):
            rule = {
                "refold_rule": "built_and_flagged_multipass_disposition_rule",
                "refold_rule_metric": d["refold_verdict"].removeprefix("fail: "),
                "refold_rule_threshold": "fold intact (native-region lDDT >= 0.85, site lDDT >= "
                "0.70 where scored), ALFA RSA >= 0.25; class epitope_packed or "
                "low_confidence_burial",
                "refold_rule_decided_by": DISPOSITION_BY,
                "refold_rule_date": "2026-11-14",
            }
        result.append(
            d
            | rule
            | {
                "igk_deletion_evidence": deletion_evidence(d),
                "refold_disposition": d.get("refold_disposition", ""),
                "weak_sp_flag": d.get("weak_sp_flag", ""),
            }
        )
    return result


def deletion_evidence(d: dict) -> str:
    """What shows that residues 1..c an Igκ swap deletes are a cleaved signal peptide."""
    if d["arm"] != "igk_swap":
        return ""
    if d["design_id"] in DEFERRED_ROUTES:
        return "none: signal anchor, kept until the receiving rung's arm passes refold"
    c = d["igk_replaces_residues"].split("-")[1]
    if d.get("sp_review_verdict") == "cleaved_sp_confirmed":
        return f"second source: {SP_EVIDENCE.get(d['candidate_id'], '')} (deletes 1-{c})"
    if d.get("sp_review_verdict") == "unconfirmed":  # re-anchored on curated evidence
        return f"re-anchored: {SP_EVIDENCE.get(d['candidate_id'], '')} (deletes 1-{c})"
    call = NATIVE_SP.get(d["candidate_id"])
    if call:
        return f"SignalP 6.0 native SP 1-{call[0]} (p={call[1]}) (deletes 1-{c})"
    return ""


SP_EVIDENCE: dict[str, str] = {}
NATIVE_SP: dict[str, tuple] = {}


ALL_ALTERNATIVES: list[dict] = []
REFOLD_ALTERNATIVES = HERE / "refold_requests_alternatives.tsv"


# ---- weak predicted signal peptide -> Igκ ----

WEAK_SP = HERE / "weak_sp.tsv"
WEAK_SP_BY = (
    "coordinator under Becca's delegation, 2026-11-15 (weak predicted SP: native or tagged "
    "cleavage-site probability < 0.5, proposed by nterm-reviewer-package, not yet confirmed "
    "by Becca)"
)


def apply_weak_sp(out, row):
    """Swap a weak-SP native arm for its Igκ form when that form keeps the leader's cleavage and
    passes refold (or the disposition rule); otherwise keep the native arm, flagged.

    The member list and Igκ forms are ``weak_sp.tsv`` (``weak_sp.py``, 2026-11-15), applied only
    while the member's native arm is the one recorded there (``native_tagged_sha256``).
    """
    if not WEAK_SP.exists():
        return out
    weak = {r["design_id"]: r for r in _tsv(WEAK_SP)}
    verdicts = load_verdicts_cached()
    result = []
    for d in out:
        w = weak.get(d["design_id"])
        d = d | {"weak_sp_flag": ""}
        if not w or d["arm"] != "native_site" or d["tagged_sha256"] != w["native_tagged_sha256"]:
            result.append(d)
            continue
        cuts = (
            f"cut calls: SignalP {w['signalp_cut']} (p={w['native_cleavage_probability']}), "
            f"DeepTMHMM v1 {w['deeptmhmm_v1_cut'] or 'none'}; c is the later call"
        )
        if w["status"]:
            result.append(d | {"weak_sp_flag": f"weak predicted SP; {cuts}; {w['status']}"})
            continue
        if w["signalp_verdict"] != "preserved":
            result.append(
                d
                | {
                    "weak_sp_flag": "weak predicted SP; Igκ form failed SignalP "
                    f"({w['signalp_igk_tagged']}, {w['signalp_verdict']}); native arm kept"
                }
            )
            continue
        v = verdicts.get(w["tagged_sha256"])
        k = qualifies(v["row"]) if v and not v["pass"] else ""
        if not v or (not v["pass"] and not k):
            why = "refold pending" if not v else f"refold failed ({v['why']})"
            result.append(
                d | {"weak_sp_flag": f"weak predicted SP; Igκ form {why}; native arm kept"}
            )
            continue
        cp = w["signalp_igk_tagged"].split("p=")[-1]
        new = row(
            d["candidate_id"],
            "igk_swap",
            int(w["insert_after_residue"]),
            int(w["native_after_residue"]),
            20,
            d["inside_predicted_fold"],
            cp,
            w["tagged_sha256"],
            WEAK_SP_BY,
            f"weak predicted SP (native cleavage-site probability "
            f"{w['native_cleavage_probability']}, tagged {w['tagged_cleavage_probability']}); "
            f"{cuts}; the Igκ leader replaces {w['igk_replaces']} and keeps its cleavage",
            f"the native arm after {d['native_after_residue']} (tagged cleavage probability "
            f"{d['tagged_cleavage_probability']})",
            "weak-SP rule + SignalP + l1-tag-refold verdict",
            refold="pass"
            if v["pass"] and not v["why"].startswith("borderline")
            else v["why"]
            if v["pass"]
            else f"fail: {v['why']}",
        )
        new["igk_replaces_residues"] = w["igk_replaces"]
        if k:
            new["refold_disposition"] = disposition_text(v["row"], k)
        result.append(new | {"weak_sp_flag": f"weak predicted SP; {cuts}; swapped to Igκ"})
    return result


#: Held members another rung has taken, retired here before that rung's lane override reaches
#: main (after which they leave this rung's membership altogether). That rung's override is the
#: record; this rung adds none.
RETIRED = {
    "l1-p1362": "routed to the residue-1 multipass rung (Igκ + tag at Met1; lane override on "
    "branch claude/l1-res1-slc41a2, 33f26feb): SignalP 6.0 and UniProt call no signal peptide, "
    "and the leader-free N-terminally tagged form was displayed in Tedman's screen",
}


# ---- ortholog rule ----

#: Becca, 2026-11-14 (relayed by the coordinator): "Never substitute a mouse or cyno ortholog
#: for a human and don't build orthologs for genes with no good human proteoform." An ortholog
#: arm is built only if a human member of its gene is built, on this rung or another
#: (``constructs_l1_*.tsv`` other than this rung's own).
ORTHOLOG_RULE_BY = "Becca, 2026-11-14 (ortholog rule)"


def human_built_elsewhere() -> set[str]:
    genes = set()
    for f in sorted(S4E9.glob("constructs_l1_*.tsv")):
        if f.name in ("constructs_l1_nterm.tsv", "constructs_l1_type3_proteins.tsv"):
            continue
        for r in _tsv(f):
            if r.get("candidate_set") != "ortholog" and r.get("species", "human") == "human":
                genes.add(r["gene_symbol"])
    return genes


def apply_ortholog_rule(out, row):
    here = {
        r["gene_symbol"]
        for r in out
        if r["arm"] in ("native_site", "igk_swap") and r["candidate_set"] != "ortholog"
    }
    built = here | human_built_elsewhere()
    result = []
    for d in out:
        if (
            d["candidate_set"] == "ortholog"
            and d["arm"] in ("native_site", "igk_swap")
            and d["gene_symbol"] not in built
        ):
            blank = dict.fromkeys(d, "")
            d = row(
                d["candidate_id"],
                "no_nterm_arm",
                "",
                "",
                "",
                "",
                "",
                "",
                ORTHOLOG_RULE_BY,
                "ortholog rule: no human member of the gene is built, so the ortholog is not "
                "built in its place",
                f"build the ortholog's {d['arm']} after {d['insert_after_residue']}",
                d["evidence"],
                refold=d["refold_verdict"],
            )
            d = blank | d
        result.append(d)
    return result


# ---- built-and-flagged dispositions (the multipass rung's rule) ----

#: Decided by the coordinator under Becca's 2026-11-12 delegation and her 2026-11-14 rescue
#: request: this rung takes the residue-1 multipass rung's disposition rule
#: (``multipass_refold_dispositions.py``) verbatim, restricted to failures with the fold intact.
DISPOSITION_BY = (
    "coordinator under Becca's 2026-11-12 delegation and 2026-11-14 rescue request "
    "(multipass disposition rule)"
)
CONFIDENT_EPITOPE = 5  # multipass_refold_dispositions.CONFIDENT_EPITOPE
LOW_PLDDT = 50
RSA_EXPOSED = 0.25  # the refold session's burial threshold
DISPOSITION_TEXT = {
    "epitope_packed": "the model packs epitope residues on the native domain; fold intact; "
    "same rule as the multipass rung; read ALFA as the display control",
    "low_confidence_burial": "SpyTag003 called buried with no confident epitope packing; fold "
    "intact; same rule as the multipass rung; read ALFA as the display control",
}
BENCH = (
    "bench SpyCatcher003:ALFA capture ratio on non-permeabilized cells: ALFA stains but "
    "SpyCatcher003 does not means SpyTag003 is occluded; move it for L2"
)


def failure_class(r: dict) -> str:
    """``multipass_refold_dispositions.klass``, restated (same order, same thresholds)."""
    reasons = str(r.get("verdict_reasons", ""))
    if float(r.get("epitope_confident_contact_res") or 0) >= CONFIDENT_EPITOPE:
        return "epitope_packed"
    if "native-region lDDT" in reasons:
        return "fold_disturbed"
    if "within 15 of the site" in reasons:
        return "tail_local"
    return "low_confidence_burial"


def qualifies(r: dict) -> str:
    """The class if the multipass rule builds it here with the fold intact, else ''."""
    reasons = str(r.get("verdict_reasons", ""))
    if "helix inside the cassette" in reasons:  # no multipass arm builds this reason class
        return ""
    k = failure_class(r)
    if k not in ("epitope_packed", "low_confidence_burial"):
        return ""
    if k == "low_confidence_burial" and "SpyTag003 buried" not in reasons:
        return ""
    lc, ls, alfa = (r.get(c) for c in ("lddt_confident", "lddt_site", "alfa_rsa"))
    if not lc or float(lc) < 0.85 or (ls not in ("", None) and float(ls) < 0.70):
        return ""
    if not alfa or float(alfa) < RSA_EXPOSED or "ALFA buried" in reasons:
        return ""
    return k


def disposition_pick(d, alts, verdicts, replaced):
    """Best failed SignalP-passing form the rule builds: fewest packed epitope residues, then
    the highest site lDDT."""
    forms = [
        *alts,
        {
            "arm": d["arm"],
            "n": d["insert_after_residue"],
            "native_n": d["native_after_residue"],
            "sha": d["tagged_sha256"],
            "strict": d["inside_predicted_fold"] != "1",
            "cp": d["tagged_cleavage_probability"],
        },
    ]
    ok = []
    for o in forms:
        v = verdicts.get(o["sha"])
        if not v or v["pass"] or (o["arm"] == "igk_swap" and d["candidate_id"] not in replaced):
            continue
        k = qualifies(v["row"])
        if k:
            r = v["row"]
            ok.append(
                (
                    float(r.get("epitope_confident_contact_res") or 0),
                    -float(r.get("lddt_site") or 1),
                    o,
                    r,
                    k,
                )
            )
    if not ok:
        return None
    _, _, o, r, k = min(ok, key=lambda x: x[:2])
    return o, r, k


def disposition_text(r: dict, k: str) -> str:
    return (
        f"built and flagged ({k}): model packs {int(float(r.get('epitope_confident_contact_res') or 0))} "
        f"epitope residues on the native domain; SpyTag003 RSA {r.get('spytag003_rsa')}, ALFA RSA "
        f"{r.get('alfa_rsa')}; native-region lDDT {r.get('lddt_confident')}, site lDDT "
        f"{r.get('lddt_site') or 'not scored'}; same rule as the multipass rung. Would change: "
        f"{BENCH}"
    )


# ---- N-terminal pro-peptides ----

PROPEP_FASTA = HERE / "propeptide_inputs_new.fasta"


def apply_propeptides(out, members, calls, row):
    """Move a cassette that sits inside a constitutively removed N-terminal pro-peptide.

    Decided by nterm-reviewer-package (delegated by Becca 2026-11-12), during the overnight
    review of every design step: the site rule knew the signal peptide and disorder only, so
    on 63 arms the cassette sat inside a pro-domain that furin or a related convertase removes
    on the way to the surface, taking the cassette with it. The cassette moves to just after
    the pro-peptide (UniProt PROPEP, ``propeptide_check.py``), i.e. to the N-terminus of the
    mature protein, which keeps the R-X-[KR]-R site whole. MMP14 and BACE1 junctions at that
    position are in the literature tag-site corpus. The arm (native or Igκ) is kept; the new
    form must still keep the SignalP cleavage site and pass refold. Activation-only
    pro-peptides (``ppc.ACTIVATION_ONLY``) and RGMA (``ppc.STAYS_ATTACHED``) keep their site.
    """
    import build_constructs as bc

    from deliverome_analysis.s3 import processed

    entries = ppc.entries_by_gene()
    verdicts = load_verdicts(processed)
    mtsv = {r["candidate_id"]: r for r in _tsv(HERE / "nterm_members.tsv")}
    hi = {k: r["interval_hi"] for k, r in mtsv.items()}
    result, requests, need = [], [], {}
    for d in out:
        gene = d["gene_symbol"]
        if (
            d["arm"] == "no_nterm_arm"
            or gene in ppc.ACTIVATION_ONLY
            or gene in ppc.STAYS_ATTACHED
            or d["candidate_id"] not in members
        ):
            result.append(d)
            continue
        cid, seq = d["candidate_id"], members[d["candidate_id"]]["sequence"]
        c, n = ppc.native_cleavage(d), int(d["native_after_residue"])
        _, props = ppc.nterm_propeptides(seq, entries.get(gene, []), c)
        hit = [p for p in props if p[2] != "unmapped" and p[0] - 1 <= n < p[1]]
        if not hit:
            result.append(d)
            continue
        a, b = min(p[0] for p in hit), max(p[1] for p in hit)
        if d["arm"] == "igk_swap":
            parent, q = LEADER_AA + PLUS1_AA + seq[c:], HEAD + (b - c)
        else:
            parent, q = seq, b
        tagged = parent[:q] + bc.CASSETTE_AA + parent[q:]
        t_sha = sha(tagged)
        pcall, tcall = calls.get(sha(parent)), calls.get(t_sha)
        if (
            pcall is None
            and d["arm"] == "native_site"
            and mtsv.get(cid, {}).get("signalp_cleavage")
        ):
            # the native SignalP call comes from the frozen D1 snapshot, not a pinned run
            pcall = {"prediction": "SP", "cleavage_site": int(mtsv[cid]["signalp_cleavage"])}
        why = f"cassette after {n} lies in the N-terminal pro-peptide {a}-{b} (UniProt, removed in the mature protein)"
        inside = int(bool(hi.get(cid)) and b > int(hi[cid] or 0))
        if tcall is None or pcall is None:
            need[t_sha] = tagged
            result.append(
                row(
                    cid,
                    d["arm"],
                    q,
                    b,
                    d["cleavage_site"],
                    inside,
                    "",
                    t_sha,
                    DECIDED_BY,
                    why + "; moved to the mature N-terminus, awaiting SignalP",
                    f"keep the cassette after {n}, inside the pro-peptide",
                    "UniProt pro-peptide",
                    refold="pending",
                )
            )
            continue
        if npl.verdict(pcall, tcall, q) != npl.PRESERVED:
            result.append(
                row(
                    cid,
                    "no_nterm_arm",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                    DECIDED_BY,
                    why + "; at the mature N-terminus SignalP no longer keeps the cleavage site",
                    f"keep the cassette after {n}, inside the pro-peptide (lost on maturation)",
                    "UniProt pro-peptide + SignalP",
                    refold="",
                )
            )
            continue
        v = verdicts.get(t_sha)
        k = qualifies(v["row"]) if v and not v["pass"] else ""
        if k:  # the multipass rung's disposition rule, as for any other failed form
            result.append(
                row(
                    cid,
                    d["arm"],
                    q,
                    b,
                    d["cleavage_site"],
                    inside,
                    tcall["cleavage_probability"],
                    t_sha,
                    DISPOSITION_BY,
                    why + "; moved to the mature N-terminus, where refold fails; built and "
                    f"flagged under the multipass rung's disposition rule ({k}): "
                    + DISPOSITION_TEXT[k],
                    "hold the member (the relocated form fails refold)",
                    "UniProt pro-peptide + SignalP + l1-tag-refold verdict + disposition rule",
                    refold=f"fail: {v['why']}",
                )
                | {"refold_disposition": disposition_text(v["row"], k)}
            )
            continue
        if v and not v["pass"]:
            result.append(
                row(
                    cid,
                    "no_nterm_arm",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                    DECIDED_BY,
                    why + f"; at the mature N-terminus the refold failed ({v['why']})",
                    f"keep the cassette after {n}, inside the pro-peptide (lost on maturation)",
                    "UniProt pro-peptide + l1-tag-refold verdict",
                    refold=f"fail: {v['why']}",
                )
            )
            continue
        result.append(
            row(
                cid,
                d["arm"],
                q,
                b,
                d["cleavage_site"],
                inside,
                tcall["cleavage_probability"],
                t_sha,
                DECIDED_BY,
                why + "; moved to the mature N-terminus, SignalP keeps the cleavage site",
                f"keep the cassette after {n}, inside the pro-peptide (lost on maturation)",
                "UniProt pro-peptide + SignalP",
                refold="pass" if v else "pending",
            )
        )
        if not v:
            requests.append(
                {
                    "candidate_id": cid,
                    "design_id": members[cid]["design_id"],
                    "arm": d["arm"],
                    "insert_after_residue": q,
                    "tagged_sha256": t_sha,
                    "replaces": d["tagged_sha256"],
                }
            )
    return result, requests, need


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("cmd", choices=["inputs", "decide"])
    ap.add_argument(
        "--verdicts",
        help="refold-verdict key to read instead of l1-tag-refold's live file, e.g. a frozen "
        "copy named in decision_inputs.json, to reproduce an earlier decision set",
    )
    ap.add_argument(
        "--hold-pending",
        action="store_true",
        help="hold every built arm that has no refold verdict (finalizing for merge)",
    )
    a = ap.parse_args(argv)
    HOLD_PENDING["on"] = a.hold_pending
    if a.verdicts:
        global VERDICTS_KEY
        VERDICTS_KEY = a.verdicts
    return inputs() if a.cmd == "inputs" else decide()


if __name__ == "__main__":
    raise SystemExit(main())
