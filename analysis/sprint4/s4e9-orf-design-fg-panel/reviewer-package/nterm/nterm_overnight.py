"""Overnight 2026-11-12: re-prediction inputs for the held and refused N-after-SP groups.

Becca, 2026-11-11: "for n term proteins run signal reprediction and disorder etc as we had
planned to select the site then also refold those". Coordinator plan
``OVERNIGHT-L1-2026-11-12.md``; this session's Modal cap is $25, for SignalP 6.0 only here.

Three steps, each a subcommand:

``snapshot``
    Read-only D1: the ortholog releases ``signalp_public`` sp6_2026_09_28_ortho and
    ``disorder_public`` dis_ortho_2026_09_28, for ortholog accessions whose ``topology_public``
    sequence is a refused member's designed protein (matched by SHA-256), with those input
    sequences. The census read only the human releases, so 121 members were refused for want
    of a SignalP row. Published under ``nterm/2026-11-12/`` (this package's prefix).
``inputs``
    Groups every held or refused member from the current ``nterm_build_options.tsv``, so it
    reruns unchanged when the membership changes (it did on 2026-11-12, PR #283):

    * ``new_member``: members with candidates but no sweep record. The native protein and
      the tandem at every candidate.
    * ``igk_failing_native``: native signal lost, or cleaved past the junction, at every
      candidate. The Igκ leader plus Asp replaces residues ``1..c``, and the tandem goes after
      the same mature residue as each native candidate.
    * ``ortholog``: refused for want of a human SignalP row but covered by the ortholog
      release. The interval is recomputed exactly as ``nterm_placement.place`` does.
    * ``igk_disagree``: DeepTMHMM v1 calls a cleaved signal peptide (its ``S`` run is followed
      by an extracellular residue, so not a signal anchor) and SignalP 6.0 does not. The
      interval uses DeepTMHMM v1's cleavage and the member's own disorder tracks; then the
      same Igκ construction.

    Every SignalP call already made is reused by sequence SHA-256 (runs in ``RUN_IDS``); only
    unseen sequences go to ``overnight_inputs_new.fasta`` for a new pinned run.
``analyse``
    Verdicts against each sequence's own untagged parent in the same run (the sweep's rule),
    and one proposed site per member.

**Attribution** (convention of 2026-11-10). Groups and rules are proposals by the Claude Code
session ``nterm-reviewer-package`` (PR #254, 2026-11-12), for Becca to decide. Admitting the
ortholog release and using DeepTMHMM v1's cleavage for the 57 are both choices she has not yet
made.

    uv run --locked python nterm_overnight.py snapshot --publish
    uv run --locked python nterm_overnight.py inputs
    uv run --locked python nterm_overnight.py analyse
"""

from __future__ import annotations

import argparse
import collections
import csv
import gzip
import hashlib
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
S4E9 = HERE.parents[1]
sys.path.insert(0, str(S4E9))
sys.path.insert(0, str(HERE))

import admissible_interval as ai  # noqa: E402
import build_constructs as bc  # noqa: E402
import nterm_build_options as nbo  # noqa: E402
import nterm_placement as npl  # noqa: E402
import snapshot_inputs as si  # noqa: E402
from build_gpcr_rescue_arms import LEADER_AA, PLUS1_AA  # noqa: E402

DATE = "2026-11-12"
#: The second ortholog snapshot (after #283 changed the membership); the first,
#: ``ortholog_rows_2026-11-12.jsonl.gz``, stays on S3 as the pre-#283 record.
ORTHO_ROWS = f"ortholog_rows_{DATE}c.jsonl.gz"  # after #283 6bf341f2
ORTHO_KEY = f"{npl.PREFIX}reviewer-package/nterm/{DATE}/{ORTHO_ROWS}"
SIGNALP_ORTHO = "sp6_2026_09_28_ortho"
DISORDER_ORTHO = "dis_ortho_2026_09_28"
ORTHO_COHORTS = ("cyno_ortholog", "cyno_ortholog_repick", "mouse_ortholog", "mouse_ortholog_repick")
INPUTS_TSV = HERE / "overnight_inputs.tsv"
INPUTS_FASTA = HERE / "overnight_inputs.fasta"
#: Every pinned SignalP run whose calls are reused, matched by sequence SHA-256. A run that
#: has not been published yet is skipped.
RUN_IDS = (
    "nterm_sp6_2026_11_10",
    "nterm_sp6_2026_11_12",
    "nterm_sp6_2026_11_12b",
    "nterm_sp6_2026_11_12c",
    "nterm_sp6_2026_11_12d",
    "nterm_sp6_2026_11_12e",
    "nterm_sp6_2026_11_12f",
    "nterm_sp6_2026_11_14g",
    "nterm_sp6_2026_11_14h",
    "nterm_sp6_2026_11_14i",
    "nterm_sp6_2026_11_14j",
    "nterm_sp6_2026_11_15a",
)
NEW_FASTA = HERE / "overnight_inputs_new.fasta"


def calls_key(run_id: str) -> str:
    return f"{npl.PREFIX}reviewer-package/nterm/signalp6-runs/{run_id}/calls.json"


PROPOSED_BY = "proposed, not decided: nterm-reviewer-package (PR #254), 2026-11-12"
IMPLEMENTED_BY = "Claude Code session nterm-reviewer-package (PR #254)"


def sha(s: str) -> str:
    return hashlib.sha256(s.encode()).hexdigest()


def _tsv(path: Path) -> list[dict]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def refused(cls: str) -> dict[str, dict]:
    return {r["candidate_id"]: r for r in _tsv(HERE / "nterm_members.tsv") if r["class"] == cls}


# ---- snapshot ----


def snapshot(publish: bool) -> int:
    from deliverome_analysis import topology
    from deliverome_analysis.s3 import processed

    lane = si.members()
    shas = {lane[c]["sha"] for c in refused("refused_no_signalp_row")}
    snap = topology.load_snapshot()
    hit = snap[snap["sequence_sha256"].isin(shas) & snap["cohort"].isin(ORTHO_COHORTS)]
    accs = sorted(set(hit["uniprot_acc_full"]))
    q = si.d1_query()
    records = []
    for chunk in si._chunks(accs, 20):
        marks = ",".join("?" * len(chunk))
        for r in q(
            "select uniprot_acc_full, topology_version, cohort, sequence from "
            f"topology_public where uniprot_acc_full in ({marks})",
            list(chunk),
        ):
            s = r.pop("sequence")
            records.append(
                {
                    "kind": "input_sequence",
                    "source": "topology_public",
                    "key": r["uniprot_acc_full"],
                    "version": r["topology_version"],
                    "cohort": r["cohort"],
                    "sequence_sha256": sha(s),
                    "sequence": s,
                }
            )
    for chunk in si._chunks(accs, 40):
        marks = ",".join("?" * len(chunk))
        for r in q(
            f"select uniprot_acc_full, {', '.join(si.SIGNALP_FIELDS)} from signalp_public "
            f"where signalp_version = ? and uniprot_acc_full in ({marks})",
            [SIGNALP_ORTHO, *chunk],
        ):
            key = r.pop("uniprot_acc_full")
            records.append(
                {
                    "kind": "signalp",
                    "source": "signalp_public",
                    "key": key,
                    "version": SIGNALP_ORTHO,
                }
                | r
            )
    for chunk in si._chunks(accs, 12):
        marks = ",".join("?" * len(chunk))
        for r in q(
            f"select uniprot_acc_full, {', '.join(si.DISORDER_FIELDS)} from "
            f"disorder_public where disorder_version = ? and uniprot_acc_full in ({marks})",
            [DISORDER_ORTHO, *chunk],
        ):
            key = r.pop("uniprot_acc_full")
            records.append(
                {
                    "kind": "disorder",
                    "source": "disorder_public",
                    "key": key,
                    "version": DISORDER_ORTHO,
                }
                | r
            )
    out = HERE / "refold_batches" / ORTHO_ROWS
    out.parent.mkdir(exist_ok=True)
    si.write_rows(records, out)
    print(f"{len(accs)} ortholog accessions; {collections.Counter(r['kind'] for r in records)}")
    if publish:
        if processed.exists(ORTHO_KEY):
            raise SystemExit(f"processed://{ORTHO_KEY} exists")
        processed.upload(out, ORTHO_KEY)
        print(f"uploaded processed://{ORTHO_KEY}")
    return 0


# ---- inputs ----


def ortho_evidence(member: dict, rows: list[dict]):
    """SignalP row and tracks from an ortholog accession whose every input sequence is the member's."""
    seqs, sp, dis = collections.defaultdict(set), {}, collections.defaultdict(dict)
    for r in rows:
        if r["kind"] == "input_sequence":
            seqs[r["key"]].add(r["sequence_sha256"])
        elif r["kind"] == "signalp":
            sp[r["key"]] = r
        else:
            dis[r["key"]][r["predictor"]] = r
    exact = sorted(a for a, s in seqs.items() if s == {member["sha"]} and a in sp)
    if not exact:
        return None, {}, ""
    acc = exact[0]
    tracks = {}
    for pred, row in sorted(dis[acc].items()):
        t = npl.decode_track(row)
        if t:
            tracks[pred] = t
    return sp[acc], tracks, acc


def native_rows(cid, m, c, cands: list[tuple[int, bool]], group: str) -> list[dict]:
    out = [
        {
            "group": group,
            "candidate_id": cid,
            "arm": "native",
            "after_residue": "",
            "native_after_residue": "",
            "cleavage_site": c,
            "strict": "",
            "sequence": m["sequence"],
        }
    ]
    for n, strict in cands:
        out.append(
            {
                "group": group,
                "candidate_id": cid,
                "arm": "tagged",
                "after_residue": n,
                "native_after_residue": n,
                "cleavage_site": c,
                "strict": int(strict),
                "sequence": m["sequence"][:n] + bc.CASSETTE_AA + m["sequence"][n:],
            }
        )
    return out


def igk_rows(cid, m, c, cands: list[tuple[int, bool]], group: str) -> list[dict]:
    native = m["sequence"]
    igk = LEADER_AA + PLUS1_AA + native[c:]
    head = len(LEADER_AA) + len(PLUS1_AA)
    out = [
        {
            "group": group,
            "candidate_id": cid,
            "arm": "igk_native",
            "after_residue": "",
            "native_after_residue": "",
            "cleavage_site": c,
            "strict": "",
            "sequence": igk,
        }
    ]
    for n, strict in cands:
        q = head + (n - c)
        assert igk[head:q] == native[c:n]
        out.append(
            {
                "group": group,
                "candidate_id": cid,
                "arm": "igk_tagged",
                "after_residue": q,
                "native_after_residue": n,
                "cleavage_site": c,
                "strict": int(strict),
                "sequence": igk[:q] + bc.CASSETTE_AA + igk[q:],
            }
        )
    return out


def prior_calls(processed, include_new: bool = True) -> dict:
    """Calls from every published run in RUN_IDS, keyed by sequence SHA-256."""
    calls = {}
    for run in RUN_IDS:
        if processed.exists(calls_key(run)):
            calls |= npl.load_json(processed, calls_key(run))
    return calls


def inputs() -> int:
    """Group every held or refused member from the current tables; write only unseen sequences."""
    from deliverome_analysis.s3 import processed

    members, _ = npl.load_members(processed)
    base = npl.load_snapshot(processed)
    with gzip.open(processed.download(ORTHO_KEY), "rt") as fh:
        ortho = [json.loads(line) for line in fh]
    opts = {r["candidate_id"]: r for r in _tsv(HERE / "nterm_build_options.tsv")}
    pos = collections.defaultdict(list)
    for p in _tsv(HERE / "nterm_positions.tsv"):
        pos[p["candidate_id"]].append((int(p["after_residue"]), p["within_strict_interval"] == "1"))
    rows, notes = [], []
    for cid in sorted(members):
        m, sub = members[cid], opts[cid]["subgroup"]
        cands = sorted(pos.get(cid, []))
        if sub == "not_repredicted":
            rows += native_rows(cid, m, int(opts[cid]["cleavage_site"]), cands, "new_member")
        elif sub in ("signal_peptide_lost_everywhere", "cleavage_crosses_junction_everywhere"):
            rows += igk_rows(cid, m, int(opts[cid]["cleavage_site"]), cands, "igk_failing_native")
        elif sub == "refused_no_signalp_row":
            signalp, tracks, acc = ortho_evidence(m, ortho)
            if signalp is None:
                notes.append({"candidate_id": cid, "group": "ortholog", "note": "no ortholog row"})
                continue
            pl = npl.place(m, signalp, tracks)
            notes.append(
                {"candidate_id": cid, "group": "ortholog", "note": pl["class"], "accession": acc}
            )
            if pl["candidates"]:
                rows += native_rows(
                    cid,
                    m,
                    pl["lo"],
                    [(n, n <= pl["median_onset"] - 1) for n in pl["candidates"]],
                    "ortholog",
                )
        elif sub == "refused_sp_disagree_deeptmhmm_v1_only":
            c = m["v1_cleavage"]
            _, tracks = npl.evidence(npl.route(m, None, base), base)
            if not tracks:
                notes.append({"candidate_id": cid, "group": "igk_disagree", "note": "no tracks"})
                continue
            strict = ai.tag_interval(c, tracks, tolerance=0)
            cands = ai.candidate_positions(ai.tag_interval(c, tracks), c)
            notes.append(
                {
                    "candidate_id": cid,
                    "group": "igk_disagree",
                    "note": f"interval {c}-{strict.hi + 2}",
                }
            )
            rows += igk_rows(cid, m, c, [(n, n <= strict.hi) for n in cands], "igk_disagree")
    for r in rows:
        r["design_id"] = members[r["candidate_id"]]["design_id"]
        r["gene_symbol"] = members[r["candidate_id"]]["gene_symbol"]
        r["sequence_sha256"] = sha(r["sequence"])
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
    known = prior_calls(processed)
    seen = set()
    with NEW_FASTA.open("w") as fh:
        for r in rows:
            if r["sequence_sha256"] in seen or r["sequence_sha256"] in known:
                continue
            seen.add(r["sequence_sha256"])
            fh.write(f">{r['sequence_sha256']}\n{r['sequence']}\n")
    (HERE / "overnight_input_notes.json").write_text(
        json.dumps(
            {
                "by_group_note": dict(
                    collections.Counter(f"{n['group']}:{n['note'].split(' ')[0]}" for n in notes)
                ),
                "notes": notes,
            },
            indent=1,
        )
        + "\n"
    )
    print(collections.Counter((r["group"], r["arm"]) for r in rows))
    print(collections.Counter(f"{n['group']}:{n['note'].split(' ')[0]}" for n in notes))
    print(f"{len(rows)} rows; {len(seen)} sequences not yet predicted -> {NEW_FASTA.name}")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("snapshot").add_argument("--publish", action="store_true")
    sub.add_parser("inputs")
    sub.add_parser("analyse")
    a = ap.parse_args(argv)
    if a.cmd == "snapshot":
        return snapshot(a.publish)
    if a.cmd == "inputs":
        return inputs()
    return analyse()


def _pick(tagged: list[dict], parent_cp: float) -> dict | None:
    """Preserved sites only, ordered by ``nterm_build_options.pick_key``."""
    ok = [t for t in tagged if t["verdict"] == "preserved"]
    if not ok:
        return None
    return min(
        ok,
        key=lambda t: nbo.pick_key(
            t["strict"] == "1", parent_cp, t["cp"] or 0, int(t["native_after_residue"])
        ),
    )


def analyse() -> int:
    """Verdicts against each parent in the same calls store, and one proposed site per member."""
    from deliverome_analysis.s3 import processed

    members, _ = npl.load_members(processed)
    calls = prior_calls(processed, include_new=True)
    by_member = collections.defaultdict(lambda: {"parent": None, "tagged": []})
    verdict_rows = []
    for r in _tsv(INPUTS_TSV):
        c = calls[r["sequence_sha256"]]
        rec = {
            **r,
            "prediction": c["prediction"],
            "cleavage": c["cleavage_site"],
            "cp": c["cleavage_probability"],
            "sp_probability": c["sp_probability"],
        }
        key = "parent" if r["arm"] in ("native", "igk_native") else "tagged"
        if key == "parent":
            by_member[r["candidate_id"]]["parent"] = rec
        else:
            by_member[r["candidate_id"]]["tagged"].append(rec)
        verdict_rows.append(rec)
    for d in by_member.values():
        par = d["parent"]
        for t in d["tagged"]:
            t["verdict"] = (
                "parent_has_no_signal_peptide"
                if par["prediction"] != "SP"
                else npl.verdict(
                    # the parent's own SignalP call, not the design cleavage in the input table
                    {"prediction": par["prediction"], "cleavage_site": par["cleavage"]},
                    {"prediction": t["prediction"], "cleavage_site": t["cleavage"]},
                    int(t["after_residue"]),
                )
            )
    proposals = []
    for cid, d in sorted(by_member.items()):
        par, m = d["parent"], members[cid]
        pick = _pick(d["tagged"], par["cp"] or 0) if par["prediction"] == "SP" else None
        group = par["group"]
        igk = group.startswith("igk")
        proposals.append(
            {
                "candidate_id": cid,
                "design_id": m["design_id"],
                "gene_symbol": m["gene_symbol"],
                "group": group,
                "proposal": ("igk_swap" if igk else "native_site") if pick else "hold",
                "insert_after_residue": pick["after_residue"] if pick else "",
                "native_after_residue": pick["native_after_residue"] if pick else "",
                "inside_predicted_fold": "" if not pick else int(pick["strict"] != "1"),
                "parent_cleavage": par["cleavage"] or "",
                "parent_cleavage_probability": par["cp"] or "",
                "tagged_cleavage_probability": pick["cp"] if pick else "",
                "n_tagged": len(d["tagged"]),
                "n_preserved": sum(t["verdict"] == "preserved" for t in d["tagged"]),
                "tagged_sha256": pick["sequence_sha256"] if pick else "",
                "decided_by": PROPOSED_BY,
                "implemented_by": IMPLEMENTED_BY,
                "reason": (
                    (
                        "Igκ leader replaces 1..c; SignalP keeps the leader's cleavage with the "
                        "tandem in place"
                    )
                    if igk
                    else (
                        "re-predicted native site (new member)"
                        if group == "new_member"
                        else "ortholog SignalP release admitted; gated site"
                    )
                )
                if pick
                else (
                    "parent has no SignalP signal peptide"
                    if par["prediction"] != "SP"
                    else "no tagged form keeps the cleavage site"
                ),
            }
        )
    with (HERE / "overnight_proposals.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(proposals[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(proposals)
    with (HERE / "overnight_verdicts.tsv").open("w", newline="") as fh:
        cols = [
            "group",
            "candidate_id",
            "design_id",
            "gene_symbol",
            "arm",
            "after_residue",
            "native_after_residue",
            "strict",
            "prediction",
            "cleavage",
            "cp",
            "sp_probability",
            "verdict",
            "sequence_sha256",
        ]
        w = csv.DictWriter(
            fh,
            fieldnames=cols,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
            restval="",
        )
        w.writeheader()
        w.writerows(verdict_rows)
    summary = {
        "runs": list(RUN_IDS),
        "proposals": dict(
            sorted(collections.Counter(f"{p['group']}:{p['proposal']}" for p in proposals).items())
        ),
        "proposal_inside_fold": dict(
            sorted(
                collections.Counter(
                    f"{p['group']}:{p['inside_predicted_fold']}"
                    for p in proposals
                    if p["proposal"] != "hold"
                ).items()
            )
        ),
        "tagged_verdicts": dict(
            sorted(
                collections.Counter(
                    f"{t['group']}:{t['verdict']}" for d in by_member.values() for t in d["tagged"]
                ).items()
            )
        ),
        "holds": sorted(
            f"{p['candidate_id']} {p['gene_symbol']} ({p['group']}: {p['reason']})"
            for p in proposals
            if p["proposal"] == "hold"
        ),
    }
    (HERE / "overnight_summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps({k: v for k, v in summary.items() if k != "holds"}, indent=1))
    print(f"holds: {len(summary['holds'])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
