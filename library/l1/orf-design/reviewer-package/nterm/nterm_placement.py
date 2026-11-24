"""Reproduce N-terminal tag placement for every L1 "N after SP" member, and check its frame.

Offline, from pinned inputs only (run it through ``run_nterm_review.py`` to make that
enforced). For each current L1 member whose DeepTMHMM v1 topology puts a signal peptide in
front of an extracellular N-terminus, it recomputes the three steps of placement and checks
each against the record:

(a) **Signal peptide.** The SignalP 6.0 call that describes the member's own designed
    protein, chosen by sequence (below), and the DeepTMHMM v1 signal-peptide length.
(b) **Admissible interval.** The later of the two cleavage calls, the three disorder
    predictors' structure onsets, Tedman's interval and its candidate insertion sites
    (``admissible_interval.py``). Compared, member by member, with the 2026-10-18 interval
    census (``processed://.../triage/l1_interval_census.tsv``).
(c) **Re-prediction on the tagged protein.** Every candidate is built as
    ``native[:n] + cassette + native[n:]``, hashed, and looked up in the 9,046-sequence
    SignalP sweep ``sweep_rp_1792381414``. The verdict is recomputed from the raw native and
    tagged calls and compared with the sweep's own.

**Which predictor rows describe a member: by sequence, never by accession or length.** The
accession-keyed rows (``signalp_public``, ``disorder_public``) were predicted on the
``topology_public`` sequence stored under that accession; the snapshot carries every such
sequence. A row is used as

* ``exact``: every stored input sequence for that accession equals the designed protein; or
  a ``*_designed_public`` row whose ``sequence_sha256`` is the designed protein's;
* ``prefix``: the census's own accession, whose input sequence first differs from the
  designed protein at residue ``d``; SignalP is used only if ``d > 70`` (SignalP reads the
  first 70 residues), and disorder tracks are cut before ``d``, as the census did;
* otherwise not at all.

AlphaFold-disorder is ``1 - pLDDT`` from the AlphaFold DB model of the accession; the model's
own sequence is not in the snapshot, so for that track the match is to the accession, not
proven on the sequence. It is reported as such.

**Frame checks**, per candidate site (no pick is changed; failures are reported):

* the cleavage site and the insertion residue lie inside the designed protein;
* the residues either side of the junction are recorded, read from the designed protein;
* inserting the p003 cassette DNA at nt ``3n`` of the designed CDS translates to exactly
  ``native[:n] + cassette + native[n:]`` with no internal stop;
* the sweep's native record for the member is the designed protein (SHA-256), and a tagged
  record is found only by the SHA-256 of that exact tagged protein.

**Attribution of the judgements this file applies** (convention of 2026-11-10: Becca asked
that every manual fix name who decided it, which agent implemented it, and why). The rows
each one touches carry ``override``, ``decided_by``, ``implemented_by`` and
``override_reason`` columns, so the attribution survives when a table is read alone.

* ``duplicate_sweep_record``: a tagged design that a retired duplicate of the same protein
  submitted under its own design id is accepted as this member's re-prediction, because the
  sweep keeps one record per sequence. Decided and implemented by the Claude Code session
  ``nterm-reviewer-package`` (PR #254, 2026-11-04); not reviewed by Becca. It touches 3
  members (ACHE, CD177, UPK3BL1).
* ``prefix_accession_rows``: an accession whose sequence first differs from the designed
  protein after residue 70 lends its SignalP call, and its disorder tracks cut before the
  difference. This is the interval census's rule (PR #191, 2026-10-18), restated here by
  ``nterm-reviewer-package`` (PR #254). It touches 9 members.

Everything else is computed from the pinned inputs with no hand-made choice. The rule itself
(Tedman's interval, the later cleavage call, the +20 cap) is the record's, as
``admissible_interval.py`` cites.

Writes ``nterm_members.tsv``, ``nterm_positions.tsv`` and ``nterm_numbers.json`` to ``--out``
(default: this folder).

    uv run --locked python nterm_placement.py
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

from deliverome_analysis import library_design as ld  # noqa: E402
from deliverome_analysis import seqtools, tag_placement  # noqa: E402

PREFIX = "sprint4/s4e9-orf-design-fg-panel/"
#: 2026-11-21: the 2026-11-12b rows unchanged, plus the D1 rows (same pinned versions) for the 20
#: human accessions of the co-receptor additions on this rung (snapshot_inputs.fetch, merged by
#: Claude Code session l1-review-export-scope (subagent); Becca, 2026-11-21: "rebuild now").
#: 2026-11-22: plus the D1 rows (same pinned versions) for the round-2 additions on this rung
#: and then for the 2026-11-22 isoform redraw (2026-11-22b) (add_members.py snapshot; Claude Code
#: session l1-review-export-scope (subagent)).
#: 2026-11-22c: plus the D1 rows for the 25 members whose reference moved to the GTF MANE Select
#: protein on this rung and the ITGA6, PAM and KIT refills (add_members.py snapshot; Claude Code
#: session l1-review-export-scope (subagent)).
#: 2026-11-23: plus the D1 rows for the 23 new/switched members on this rung in the
#: 2026-11-23-pre-freeze membership cut (add_members.py snapshot).
SNAPSHOT = PREFIX + "reviewer-package/nterm/2026-11-23-pre-freeze/"
ROWS_KEY = SNAPSHOT + "predictor_rows.jsonl.gz"
CANARY_KEY = SNAPSHOT + "canary_rp_1792371023.reprediction.json"
CENSUS_KEY = PREFIX + "triage/l1_interval_census.tsv"
SWEEP = PREFIX + "reprediction/sweep_rp_1792381414/"
SWEEP_SUBMISSION_KEY = SWEEP + "submission.json"
SWEEP_REPREDICTION_KEY = SWEEP + "reprediction.json"
CANARY_SUBMISSION_KEY = PREFIX + "reprediction/canary_rp_1792371023/submission.json"
#: SLC9A3's junction sweep (26/34/40/46), which holds the as-built 34.
SLC9A3_CANARY_KEY = PREFIX + "reprediction/canary_rp_1792374944/results.json"

LANE = "N after SP"
IMPLEMENTED_BY = "Claude Code session nterm-reviewer-package (PR #254)"
OVERRIDES = {
    "duplicate_sweep_record": (
        "nterm-reviewer-package (agent), 2026-11-04; not reviewed by Becca",
        "the sweep keeps one record per sequence, and a retired duplicate of this member's "
        "identical protein submitted it under its own design id",
    ),
    "prefix_accession_rows": (
        "interval census rule, PR #191 (2026-10-18)",
        "the accession's sequence first differs from the designed protein after residue 70, "
        "so SignalP's window is identical and disorder tracks are cut before the difference",
    ),
}


def attribution(kinds: list[str]) -> dict:
    """``override`` / ``decided_by`` / ``implemented_by`` / ``override_reason`` for a row."""
    return {
        "override": ";".join(kinds),
        "decided_by": " | ".join(OVERRIDES[k][0] for k in kinds),
        "implemented_by": IMPLEMENTED_BY if kinds else "",
        "override_reason": " | ".join(OVERRIDES[k][1] for k in kinds),
    }


SIGNALP_WINDOW = 70
PREDICTORS = ("alphafold-disorder", "metapredict", "netsurfp-3.0")
CASSETTE_NT = bc.CASSETTE
CASSETTE_AA = bc.CASSETTE_AA

PLACEABLE, DEGENERATE = "placeable", "degenerate"
R_NO_SIGNALP = "refused_no_signalp_row"
R_DISAGREE_DTM = "refused_sp_disagree_deeptmhmm_v1_only"
R_DISAGREE_SP = "refused_sp_disagree_signalp_only"
R_NO_TRACK = "refused_no_disorder_track"
R_NO_CLEAVAGE = "refused_no_cleavage_call"
R_EMPTY = "refused_empty_candidate_set"

LOST, PRESERVED = "signal_peptide_lost", "preserved"
CROSSES, MOVED = "cleavage_crosses_junction", "cleavage_moved"


def sha256(s: str) -> str:
    return hashlib.sha256(s.encode()).hexdigest()


def first_difference(a: str, b: str) -> int | None:
    """1-based first residue where two sequences differ (or the shorter one ends); None if equal."""
    if a == b:
        return None
    for i, (x, y) in enumerate(zip(a, b, strict=False)):
        if x != y:
            return i + 1
    return min(len(a), len(b)) + 1


# ---- inputs ----


def load_json(bucket, key):
    with open(bucket.download(key)) as fh:
        return json.load(fh)


def load_snapshot(processed) -> dict:
    """The frozen D1 rows, indexed for routing by sequence."""
    with gzip.open(processed.download(ROWS_KEY), "rt") as fh:
        rows = [json.loads(line) for line in fh]
    out = {
        "input_sequences": collections.defaultdict(dict),
        "signalp": {},
        "disorder": collections.defaultdict(dict),
    }
    for r in rows:
        if r["kind"] == "input_sequence":
            out["input_sequences"][r["key"]][r["sequence_sha256"]] = r["sequence"]
        elif r["kind"] == "signalp":
            out["signalp"][(r["source"], r["key"])] = r
        elif r["kind"] == "disorder":
            out["disorder"][(r["source"], r["key"])][r["predictor"]] = r
    out["n_rows"] = len(rows)
    return out


def decode_track(row: dict) -> list[float] | None:
    """A disorder row's per-residue track; AlphaFold-disorder's is 1 - pLDDT from plddt_hex."""
    if row.get("sequence_substituted"):
        return None
    if row["predictor"] == "alphafold-disorder":
        h = row.get("plddt_hex")
        return [1.0 - b / 255.0 for b in bytes.fromhex(h)] if h else None
    h = row.get("disorder_hex")
    return [b / 255.0 for b in bytes.fromhex(h)] if h else None


# ---- (a) and (b): route by sequence, then the interval ----


def route(member: dict, census_row: dict | None, snap: dict) -> dict:
    """Which frozen rows describe this member's protein, decided on the sequence."""
    seq, sha = member["sequence"], member["sha"]
    exact = sorted(
        acc for acc, seqs in snap["input_sequences"].items() if seqs and set(seqs) == {sha}
    )
    census_key = census_row["lookup_key"] if census_row else ""
    census_kind = census_row["lookup_kind"] if census_row else ""
    if census_kind == "accession" and census_key in exact:
        return {"route": "exact_accession", "source_key": census_key, "first_difference": None}
    if census_kind != "designed" and exact:
        return {"route": "exact_accession", "source_key": exact[0], "first_difference": None}
    if ("signalp_designed_public", sha) in snap["signalp"] or (
        "disorder_designed_public",
        sha,
    ) in snap["disorder"]:
        return {"route": "exact_designed", "source_key": sha, "first_difference": None}
    if census_kind == "accession" and census_key in snap["input_sequences"]:
        seqs = snap["input_sequences"][census_key]
        if len(seqs) == 1:
            d = first_difference(next(iter(seqs.values())), seq)
            return {"route": "prefix_accession", "source_key": census_key, "first_difference": d}
        return {"route": "ambiguous_accession", "source_key": census_key, "first_difference": None}
    return {"route": "none", "source_key": "", "first_difference": None}


def evidence(rt: dict, snap: dict) -> tuple[dict | None, dict[str, list[float]]]:
    """(SignalP row, disorder tracks) admitted by the route, tracks cut at a first difference."""
    if rt["route"] in ("none", "ambiguous_accession"):
        return None, {}
    designed = rt["route"] == "exact_designed"
    sp_src = "signalp_designed_public" if designed else "signalp_public"
    dis_src = "disorder_designed_public" if designed else "disorder_public"
    signalp = snap["signalp"].get((sp_src, rt["source_key"]))
    tracks = {}
    for pred, row in sorted(snap["disorder"].get((dis_src, rt["source_key"]), {}).items()):
        t = decode_track(row)
        if t:
            tracks[pred] = t
    d = rt["first_difference"]
    if d is not None:
        if d <= SIGNALP_WINDOW:
            signalp = None
        tracks = {k: v[: d - 1] for k, v in tracks.items()}
    return signalp, tracks


def input_lengths(member: dict, rt: dict, snap: dict) -> tuple[bool | None, int | None]:
    """For an exact route: do the rows' recorded lengths equal the designed protein's?

    SignalP's ``protein_length``, each disorder row's ``protein_length``, and its scored
    window (the whole protein, or the first 400 residues for NetSurfP-3.0). Returns
    (ok, AlphaFold-disorder track length); the AlphaFold track's length is the model's, the
    one length check available for a track whose input sequence is not in the snapshot.
    """
    if not rt["route"].startswith("exact"):
        return None, None
    n = len(member["sequence"])
    designed = rt["route"] == "exact_designed"
    sp = snap["signalp"].get(
        ("signalp_designed_public" if designed else "signalp_public", rt["source_key"])
    )
    ok = sp is None or sp["protein_length"] == n
    af_len = None
    dis = snap["disorder"].get(
        ("disorder_designed_public" if designed else "disorder_public", rt["source_key"]), {}
    )
    for pred, row in dis.items():
        t = decode_track(row)
        if pred == "alphafold-disorder":
            af_len = len(t) if t else None
            continue
        want = min(n, 400) if pred == "netsurfp-3.0" else n
        ok &= row["protein_length"] == n and (t is None or len(t) == want)
    return ok, af_len


def place(member: dict, signalp: dict | None, tracks: dict[str, list[float]]) -> dict:
    """The census rule for one member: refusals first, then interval and candidates."""
    dtm = member["v1_cleavage"]
    sp_called = bool(signalp) and signalp["prediction"] == "SP"
    sp_cut = signalp["cleavage_site"] if sp_called else None
    reasons = []
    if signalp is None:
        reasons.append(R_NO_SIGNALP)
    elif dtm is not None and not sp_called:
        reasons.append(R_DISAGREE_DTM)
    elif sp_called and dtm is None:
        reasons.append(R_DISAGREE_SP)
    if not tracks:
        reasons.append(R_NO_TRACK)
    lo = ai.cleavage_lo(dtm, sp_cut)
    if lo is None and not reasons:
        reasons.append(R_NO_CLEAVAGE)
    out = {
        "signalp_prediction": signalp["prediction"] if signalp else "",
        "signalp_cleavage": sp_cut,
        "signalp_cleavage_probability": signalp.get("cleavage_probability") if signalp else None,
        "signalp_sp_probability": signalp.get("sp_probability") if signalp else None,
        "lo": lo,
        "lo_source": (
            ""
            if lo is None
            else "both"
            if dtm == sp_cut
            else "deeptmhmm_v1"
            if sp_cut is None or (dtm or 0) > sp_cut
            else "signalp6"
        ),
        "tracks_answered": ";".join(sorted(tracks)),
        "onsets": "",
        "median_onset": None,
        "interval_lo": None,
        "interval_hi": None,
        "candidates": [],
        "n_candidates_strict": 0,
        "all_reasons": ";".join(reasons),
    }
    if lo is None or not tracks:
        out["class"] = reasons[0]
        return out
    interval = ai.tag_interval(lo, tracks)
    strict = ai.tag_interval(lo, tracks, tolerance=0)
    observed = any(v is not None for v in interval.onsets.values())
    cands = ai.candidate_positions(interval, lo)
    # The census records each predictor's onset even on a refused row; nothing else.
    out["onsets"] = ";".join(f"{k}={v}" for k, v in sorted(interval.onsets.items()))
    if reasons:
        out["class"] = reasons[0]
        return out
    out |= {
        "median_onset": strict.hi + 1 if observed else None,
        "interval_lo": interval.lo,
        "interval_hi": interval.hi,
        "candidates": cands,
        "n_candidates_strict": len(ai.candidate_positions(strict, lo)),
    }
    out["class"] = (
        R_EMPTY if not cands else DEGENERATE if not out["n_candidates_strict"] else PLACEABLE
    )
    return out


def disorder_pick(member: dict, signalp: dict | None, tracks: dict[str, list[float]]):
    """The older POST_SP_PLACEMENT rule (most disordered residue in lo+1..lo+20), for comparison."""
    calls = [c for c in (member["v1_cleavage"], signalp and signalp.get("cleavage_site")) if c]
    if not calls or not tracks:
        return None
    try:
        return tag_placement.post_sp_junction(
            member["sequence"], calls, list(tracks.values())
        ).junction
    except ValueError:
        return None


# ---- (c) re-prediction ----


def verdict(native: dict, tagged: dict, after: int) -> str:
    """The sweep's rule, restated: lost > preserved > crosses the junction > moved."""
    if tagged["prediction"] != "SP":
        return LOST
    if tagged["cleavage_site"] == native["cleavage_site"]:
        return PRESERVED
    if tagged["cleavage_site"] >= after:
        return CROSSES
    return MOVED


# ---- frame checks ----


def frame_check(member: dict, after: int, lo: int | None) -> tuple[list[str], str]:
    """Failures for one candidate site (empty = pass) and the tagged protein."""
    seq, cds = member["sequence"], member["cds"]
    fails = []
    if lo is not None and not 0 < lo < len(seq):
        fails.append(f"cleavage {lo} outside 1..{len(seq) - 1}")
    if not 0 < after < len(seq):
        fails.append(f"insertion after {after} outside 1..{len(seq) - 1}")
    tagged = seq[:after] + CASSETTE_AA + seq[after:]
    cut = 3 * after
    tagged_cds = cds[:cut] + CASSETTE_NT + cds[cut:]
    protein = seqtools.translate(tagged_cds)
    # The terminal stop is checked once per member (cds_ends_in_stop), not here: three L1
    # CDSs have none, and that is a fact about the ORF, not about where the cassette goes.
    if protein.removesuffix("*") != tagged:
        fails.append(
            f"cassette DNA at nt {cut} does not translate to native[:n]+cassette+native[n:]"
        )
    if "*" in protein.removesuffix("*"):
        fails.append("internal stop")
    return fails, tagged


# ---- the whole run ----


def load_members(processed) -> dict[str, dict]:
    with (S4E9 / "l1_tag_paths.tsv").open() as fh:
        paths = list(csv.DictReader(fh, delimiter="\t"))
    lane = {r["candidate_id"]: r for r in paths if r["placement"] == LANE}
    proteins = ld.designed_proteins(S4E9 / "candidate_manifest_tranches.tsv", library=ld.LIBRARY_L1)
    cds = ld.designed_cds(S4E9 / "candidate_manifest_tranches.tsv", library=ld.LIBRARY_L1)
    # Members are numbered l1-mNNNN by plate position since 2026-11-10; the census and the
    # re-prediction sweep predate that and key on the design id (l1-pNNNN) instead.
    design = ld.design_ids(S4E9 / "candidate_manifest_tranches.tsv")
    out = {}
    for cid, r in sorted(lane.items()):
        seq = proteins[cid]["sequence"]
        if sha256(seq) != r["sequence_sha256"]:
            raise SystemExit(f"{cid}: l1_tag_paths.tsv sha256 is not the designed protein's")
        out[cid] = {
            "candidate_id": cid,
            "design_id": design[cid],
            "gene_symbol": r["gene_symbol"],
            "species": r["species"],
            "candidate_set": r["candidate_set"],
            "sequence": seq,
            "sha": r["sequence_sha256"],
            "cds": cds[cid],
            "v1_cleavage": int(r["signal_peptide_aa"]) or None,
        }
    return out, collections.Counter(r["placement"] for r in paths)


def run(out_dir: Path) -> dict:
    from deliverome_analysis.s3 import processed

    members, rungs = load_members(processed)
    snap = load_snapshot(processed)
    with open(processed.download(CENSUS_KEY)) as fh:
        census = {r["candidate_id"]: r for r in csv.DictReader(fh, delimiter="\t")}
    submission = load_json(processed, SWEEP_SUBMISSION_KEY)
    repred = load_json(processed, SWEEP_REPREDICTION_KEY)

    # The sweep's records, checked to be what they claim: each sequence hashes to its key.
    seq_of = {}
    for s in submission:
        if sha256(s["sequence"]) != s["sequence_sha256"]:
            raise SystemExit(
                f"sweep submission: {s['member_id']} sequence does not hash to its key"
            )
        seq_of[s["sequence_sha256"]] = s["sequence"]
    by_sha = {r["sequence_sha256"]: r for r in repred}
    sweep_native = {r["member_id"]: r for r in repred if r["role"] == "native"}
    sweep_positions = collections.defaultdict(set)
    for r in repred:
        if r["role"] == "tagged":
            sweep_positions[r["member_id"]].add(r["after_residue"])

    m_rows, p_rows, mismatches = [], [], collections.Counter()
    for cid, m in members.items():
        did = m["design_id"]
        c = census.get(did)
        rt = route(m, c, snap)
        signalp, tracks = evidence(rt, snap)
        pl = place(m, signalp, tracks)
        lengths_ok, af_len = input_lengths(m, rt, snap)
        pick = disorder_pick(m, signalp, tracks) if pl["class"] in (PLACEABLE, DEGENERATE) else None

        # Census agreement, field by field, so a disagreement says what disagrees.
        agree = {}
        if c:
            agree = {
                "class": pl["class"] == c["class"],
                "lo": str(pl["lo"] or "") == c["lo"],
                "onsets": pl["onsets"] == c["onsets"],
                "interval": (str(pl["interval_lo"] or ""), str(pl["interval_hi"] or ""))
                == (c["interval_lo"], c["interval_hi"]),
                "candidates": ",".join(map(str, pl["candidates"])) == c["candidates"],
                "sha": m["sha"] == c["sequence_sha256"],
            }
            for k, ok in agree.items():
                mismatches[k] += not ok

        nat = sweep_native.get(did)
        native_is_designed = bool(nat) and nat["sequence_sha256"] == m["sha"]
        counts = collections.Counter()
        member_fails = []
        swept_as: set[str] = set()
        for p in pl["candidates"]:
            fails, tagged = frame_check(m, p, pl["lo"])
            t = by_sha.get(sha256(tagged)) if native_is_designed else None
            # A retired duplicate of the same protein submitted the same tagged sequence, and
            # the sweep keeps one record per sequence; it is the same prediction.
            twin = (
                t is not None
                and t["member_id"] != did
                and (sweep_native.get(t["member_id"], {}).get("sequence_sha256") == m["sha"])
            )
            if t is not None and t["after_residue"] != p:
                fails.append(
                    f"sweep record {t['member_id']}@{t['after_residue']} for this sequence"
                )
            if t is not None and t["member_id"] != did and not twin:
                fails.append(f"sweep record belongs to {t['member_id']}, a different protein")
            if twin:
                swept_as.add(t["member_id"])
            if t is not None and seq_of.get(t["sequence_sha256"]) != tagged:
                fails.append("sweep tagged sequence differs from native[:n]+cassette+native[n:]")
            v = verdict(nat, t, p) if t is not None else ""
            counts[v or "not_repredicted"] += 1
            member_fails += [f"after {p}: {f}" for f in fails]
            p_rows.append(
                {
                    "candidate_id": cid,
                    "design_id": did,
                    "gene_symbol": m["gene_symbol"],
                    "after_residue": p,
                    "offset_past_cleavage": p - pl["lo"],
                    "within_strict_interval": int(p <= pl["median_onset"] - 1)
                    if pl["median_onset"]
                    else "",
                    "residue_before": m["sequence"][p - 1],
                    "residue_after": m["sequence"][p] if p < len(m["sequence"]) else "",
                    "cds_insert_nt": 3 * p,
                    "tagged_sha256": sha256(tagged),
                    "frame_ok": int(not fails),
                    "frame_failures": "; ".join(fails),
                    **attribution(
                        (["duplicate_sweep_record"] if twin else [])
                        + (["prefix_accession_rows"] if rt["route"] == "prefix_accession" else [])
                    ),
                    "repredicted": int(t is not None),
                    "tagged_prediction": t["prediction"] if t else "",
                    "tagged_cleavage": t["cleavage_site"]
                    if t and t["cleavage_site"] is not None
                    else "",
                    "tagged_cleavage_probability": t["cleavage_probability"]
                    if t and t["cleavage_probability"] is not None
                    else "",
                    "tagged_sp_probability": t["sp_probability"] if t else "",
                    "verdict": v,
                    "sweep_verdict": t["comparison"]["verdict"] if t else "",
                    "verdict_agrees": "" if not t else int(v == t["comparison"]["verdict"]),
                }
            )
        swept = sweep_positions.get(did, set()) | {
            pos for other in swept_as for pos in sweep_positions.get(other, set())
        }
        m_rows.append(
            {
                "candidate_id": cid,
                "design_id": did,
                "gene_symbol": m["gene_symbol"],
                "species": m["species"],
                "candidate_set": m["candidate_set"],
                "length_aa": len(m["sequence"]),
                "sequence_sha256": m["sha"],
                "route": rt["route"],
                "source_key": rt["source_key"],
                "first_difference": rt["first_difference"] or "",
                "alphafold_disorder_used": int("alphafold-disorder" in tracks),
                "alphafold_track_residues": af_len or "",
                "input_lengths_ok": "" if lengths_ok is None else int(lengths_ok),
                "deeptmhmm_v1_cleavage": m["v1_cleavage"] or "",
                "signalp_prediction": pl["signalp_prediction"],
                "signalp_cleavage": pl["signalp_cleavage"] or "",
                "signalp_cleavage_probability": pl["signalp_cleavage_probability"] or "",
                "lo": pl["lo"] or "",
                "lo_source": pl["lo_source"],
                "tracks_answered": pl["tracks_answered"],
                "onsets": pl["onsets"],
                "median_onset": pl["median_onset"] or "",
                "interval_lo": pl["interval_lo"] or "",
                "interval_hi": pl["interval_hi"] or "",
                "class": pl["class"],
                "all_reasons": pl["all_reasons"],
                "candidates": ",".join(map(str, pl["candidates"])),
                "n_candidates": len(pl["candidates"]),
                "n_candidates_strict": pl["n_candidates_strict"],
                "disorder_pick": pick or "",
                "disorder_pick_is_candidate": "" if pick is None else int(pick in pl["candidates"]),
                "in_census": int(bool(c)),
                "census_class": c["class"] if c else "",
                "census_agrees": "" if not c else int(all(agree.values())),
                "sweep_native_is_designed": "" if not nat else int(native_is_designed),
                "swept_under_duplicate": ",".join(sorted(swept_as)),
                "cds_ends_in_stop": int(seqtools.translate(m["cds"][-3:]) == "*"),
                "sweep_positions_match_candidates": ""
                if not nat
                else int(swept == set(pl["candidates"])),
                "n_repredicted": sum(n for k, n in counts.items() if k != "not_repredicted"),
                "n_preserved": counts[PRESERVED],
                "n_moved": counts[MOVED],
                "n_crosses_junction": counts[CROSSES],
                "n_sp_lost": counts[LOST],
                "frame_ok": int(not member_fails),
                "frame_failures": " | ".join(member_fails),
                **attribution(
                    (["duplicate_sweep_record"] if swept_as else [])
                    + (["prefix_accession_rows"] if rt["route"] == "prefix_accession" else [])
                ),
            }
        )

    numbers = summarize(m_rows, p_rows, rungs, mismatches, snap, census, members)
    numbers["panel"] = panel_check(processed, members)
    write_tsv(out_dir / "nterm_members.tsv", m_rows)
    write_tsv(out_dir / "nterm_positions.tsv", p_rows)
    (out_dir / "nterm_numbers.json").write_text(
        json.dumps(numbers, indent=1, sort_keys=True) + "\n"
    )
    return numbers


def summarize(m_rows, p_rows, rungs, mismatches, snap, census, members) -> dict:
    cls = collections.Counter(r["class"] for r in m_rows)
    in_census = [r for r in m_rows if r["in_census"]]
    has_cands = [r for r in m_rows if r["n_candidates"]]
    rep = [r for r in p_rows if r["repredicted"]]
    full = [r for r in has_cands if r["n_repredicted"] == r["n_candidates"]]
    part = [r for r in has_cands if 0 < r["n_repredicted"] < r["n_candidates"]]
    none = [r for r in has_cands if r["n_repredicted"] == 0]
    survivors = [r for r in has_cands if r["n_repredicted"] and r["n_sp_lost"] < r["n_repredicted"]]
    all_lost = [r for r in full if r["n_sp_lost"] == r["n_candidates"]]
    by_offset = collections.defaultdict(collections.Counter)
    for r in rep:
        by_offset[r["offset_past_cleavage"]][r["verdict"]] += 1
    picks = [r for r in m_rows if r["disorder_pick"] != ""]
    return {
        "tranche": "l1",
        "topology_tool": "deeptmhmm_v1",
        "rungs_l1": dict(sorted(rungs.items())),
        "nafter_sp_members": len(m_rows),
        "routes": dict(sorted(collections.Counter(r["route"] for r in m_rows).items())),
        "prefix_routes_first_difference": sorted(
            (r["candidate_id"], r["first_difference"])
            for r in m_rows
            if r["route"] == "prefix_accession"
        ),
        "alphafold_disorder_used": sum(r["alphafold_disorder_used"] for r in m_rows),
        "alphafold_track_length_differs": sorted(
            (r["candidate_id"], r["length_aa"], r["alphafold_track_residues"])
            for r in m_rows
            if r["alphafold_track_residues"] != ""
            and r["alphafold_track_residues"] != r["length_aa"]
        ),
        "input_lengths_differ": sorted(
            r["candidate_id"] for r in m_rows if r["input_lengths_ok"] == 0
        ),
        "classes": dict(sorted(cls.items())),
        "census": {
            "members_in_census": len(in_census),
            # The census is keyed by design id (l1-pNNNN); members are listed by member id.
            "census_design_ids_not_in_lane_now": sorted(
                set(census) - {m["design_id"] for m in members.values()}
            ),
            "lane_members_not_in_census": sorted(
                f"{c} ({m['design_id']})"
                for c, m in members.items()
                if m["design_id"] not in census
            ),
            "agree_all_fields": sum(r["census_agrees"] == 1 for r in in_census),
            "field_mismatches": dict(sorted(mismatches.items())),
            "disagreeing_members": sorted(
                r["candidate_id"] for r in in_census if r["census_agrees"] == 0
            ),
        },
        "lo_source": dict(sorted(collections.Counter(r["lo_source"] for r in has_cands).items())),
        "candidates": {
            "members_with_candidates": len(has_cands),
            "positions": len(p_rows),
            "per_member": dict(
                sorted(collections.Counter(r["n_candidates"] for r in has_cands).items())
            ),
        },
        "disorder_pick": {
            "members": len(picks),
            "is_a_candidate": sum(r["disorder_pick_is_candidate"] == 1 for r in picks),
        },
        "reprediction": {
            "positions_repredicted": len(rep),
            "positions_not_repredicted": len(p_rows) - len(rep),
            "members_fully_covered": len(full),
            "members_partly_covered": len(part),
            "members_not_covered": len(none),
            "members_not_covered_ids": sorted(r["candidate_id"] for r in none),
            "verdicts": dict(sorted(collections.Counter(r["verdict"] for r in rep).items())),
            "verdict_agrees_with_sweep": sum(r["verdict_agrees"] == 1 for r in rep),
            "members_with_a_surviving_position": len(survivors),
            "members_losing_every_position": len(all_lost),
            "by_offset": {str(k): dict(sorted(v.items())) for k, v in sorted(by_offset.items())},
            "sweep_native_not_designed": sorted(
                r["candidate_id"] for r in m_rows if r["sweep_native_is_designed"] == 0
            ),
            "sweep_positions_differ": sorted(
                r["candidate_id"] for r in m_rows if r["sweep_positions_match_candidates"] == 0
            ),
        },
        "cds_without_stop": sorted(
            (r["candidate_id"], r["gene_symbol"]) for r in m_rows if not r["cds_ends_in_stop"]
        ),
        "swept_under_duplicate": sorted(
            (r["candidate_id"], r["swept_under_duplicate"])
            for r in m_rows
            if r["swept_under_duplicate"]
        ),
        "frame": {
            "positions_checked": len(p_rows),
            "positions_failing": sum(1 for r in p_rows if not r["frame_ok"]),
            "members_failing": sorted(r["candidate_id"] for r in m_rows if not r["frame_ok"]),
        },
        "snapshot_rows": snap["n_rows"],
    }


def panel_check(processed, members) -> list[dict]:
    """The eight panel arms: tagged sequence, native vs tagged SignalP, and L1 identity."""
    sub = load_json(processed, CANARY_SUBMISSION_KEY)
    rp = {
        r["sequence_sha256"]: r | {"run": "canary_rp_1792371023"}
        for r in load_json(processed, CANARY_KEY)
    }
    for r in load_json(processed, SLC9A3_CANARY_KEY):
        rp.setdefault(r["sequence_sha256"], r | {"run": "canary_rp_1792374944"})
    orfs = bc.read_fasta(S4E9 / "orf_domesticated.fasta")
    with (S4E9 / "constructs_nterm_disorder_anchored.tsv").open() as fh:
        built = {
            r["gene"]: int(r["insert_after_residue"]) for r in csv.DictReader(fh, delimiter="\t")
        }
    l1_by_gene = collections.defaultdict(list)
    for m in members.values():
        l1_by_gene[m["gene_symbol"]].append(m)
    out = []
    for gene, n in sorted(built.items()):
        native = seqtools.translate(orfs[gene]).rstrip("*")
        tagged = native[:n] + CASSETTE_AA + native[n:]
        nat = next(
            (
                rp[s["sequence_sha256"]]
                for s in sub
                if s["member_id"] == gene and s["role"] == "native"
            ),
            None,
        )
        tag = rp.get(sha256(tagged))
        same = [m["candidate_id"] for m in l1_by_gene.get(gene, []) if m["sha"] == sha256(native)]
        out.append(
            {
                "gene": gene,
                "insert_after_residue": n,
                "junction": f"{native[n - 1]}{n}|{native[n]}",
                "canary_native_is_panel_orf": bool(nat)
                and nat["sequence_sha256"] == sha256(native),
                "tagged_found_by_sha": tag is not None,
                "tagged_run": tag and tag["run"],
                "native_cleavage": nat and nat["cleavage_site"],
                "native_cleavage_probability": nat and nat["cleavage_probability"],
                "tagged_prediction": tag and tag["prediction"],
                "tagged_cleavage": tag and tag["cleavage_site"],
                "tagged_cleavage_probability": tag and tag["cleavage_probability"],
                "l1_members_with_this_protein": same,
            }
        )
    return out


def write_tsv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", type=Path, default=HERE)
    args = ap.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)
    n = run(args.out)
    print(json.dumps({k: n[k] for k in ("nafter_sp_members", "routes", "classes")}, indent=1))
    print(
        json.dumps(
            n["census"] | {"disagreeing_members": n["census"]["disagreeing_members"][:10]}, indent=1
        )
    )
    print(
        json.dumps({k: v for k, v in n["reprediction"].items() if k != "by_offset"}, indent=1)[
            :2000
        ]
    )
    print(json.dumps(n["frame"], indent=1)[:1500])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
