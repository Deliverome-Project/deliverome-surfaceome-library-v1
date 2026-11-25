"""Per-member signal-peptide status and UniProt differences for the L1 review viewer.

One row per active L1 member of ``candidate_manifest_tranches.tsv``, keyed on ``design_id``
and on ``sequence_sha256`` -- the sha256 of the ordered native protein (the member's
designed CDS translated, :func:`ld.designed_proteins`), hashed with
:func:`topology.sequence_sha256`. That is the same protein ``build_library_structure_payload``
sizes its rows from, and the same hash ``l1_tag_paths.tsv`` and ``l1_proteoform_audit.tsv``
carry; the script checks both.

Columns
-------
``sp_status``
    How the built construct treats the N-terminus:

    * ``native``   -- the member keeps its own signal peptide (N-after-SP native-site arm,
      or an internal / C-terminal arm on a member with a cleaved signal peptide);
    * ``igk_swap`` -- an Igκ leader replaces the native signal peptide;
    * ``igk_res1`` -- an Igκ leader is added at residue 1 of a member with no native signal
      peptide (multipass and type III ``igk_tag`` arms);
    * ``none``     -- no signal peptide and no leader.

    Read from the construct tables (``constructs_l1_{nterm,internal,multipass,type3,cterm,
    rescue}.tsv``), ``reviewer-package/nterm/nterm_decisions.tsv`` and the manifest's
    ``has_igk_leader``. Blank when they disagree or when the member has no built or decided
    arm yet; ``sp_flags`` says which, and ``--check`` lists them. Nothing is guessed.
``native_sp_last_residue``
    Last residue of the cleaved native signal peptide, member numbering, filled even when an
    Igκ leader replaces it. In order of preference: the residues an Igκ swap deletes
    (``igk_replaces_residues``); the native cleavage site the N-terminal review decided
    (``nterm_site_picks.tsv``); otherwise the end of the leading ``S`` run of the member's
    own DeepTMHMM v1 per-residue topology (``topology.resolve``, ``deeptmhmm_v1``). Blank when
    a review overruled the v1 call (signal anchor / first TM helix, or "no signal peptide"
    lane overrides) or when the member is built leader-free at residue 1.
``sp_evidence``
    Short label of what decided the status and the SP end, taken from the decision tables'
    reason / evidence columns.
``uniprot_canonical_acc``, ``uniprot_canonical_length``, ``matches_uniprot_canonical``
    The gene's UniProt canonical from the pinned 2026-10-11 snapshot the reference ladder used
    (``library/l1/proteoform-selection/candidate_inputs.uniprot_tables``). Human members take
    the accession from ``l1_proteoform_audit.tsv`` (``canonical_acc``, joined on design id);
    orthologs take their own species' record named in the member key (base accession), or
    the record whose sequence is identical to the member; otherwise blank. No network fetch.
``uniprot_diff_json``
    Differences of the member against that canonical from a global alignment (Biopython
    ``PairwiseAligner``, BLOSUM62, open -10, extend -0.5): a JSON list of
    ``{kind, member_start, member_end, uniprot_start, uniprot_end, member_aa, uniprot_aa}``,
    1-based inclusive; the side with no residues has ``null`` coordinates and ``""``.
    ``kind`` is ``substitution`` (a run of mismatched aligned columns), ``insertion`` /
    ``deletion`` (internal residues only in the member / only in UniProt), ``n_extension`` /
    ``c_extension`` (member residues beyond UniProt's ends) or ``n_truncation`` /
    ``c_truncation`` (UniProt residues beyond the member's ends). ``[]`` when identical,
    blank when there is no canonical.
``sp_flags``
    ``;``-joined notes for review: source disagreements, undecided members, and v1 signal
    peptides that a decision overruled.

The manifest is read once (it is copied to a temporary file that every helper reads) and
re-counted at the end; the run fails if its row count changed meanwhile.

    uv run --locked python l1_member_sp_uniprot.py           # write the TSV
    uv run --locked python l1_member_sp_uniprot.py --check   # rebuild, compare, spot-check
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import re
import sys
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import library_pipeline as lp  # noqa: E402

sys.path.insert(0, str(lp.AP))
import candidate_inputs  # noqa: E402
from Bio.Align import PairwiseAligner, substitution_matrices  # noqa: E402

from deliverome_analysis import library_design as ld  # noqa: E402
from deliverome_analysis import topology  # noqa: E402

csv.field_size_limit(sys.maxsize)

MANIFEST = HERE / "candidate_manifest_tranches.tsv"
OUT = HERE / "l1_member_sp_uniprot.tsv"
NTERM = HERE / "reviewer-package" / "nterm"
SIGNALP_CACHE = NTERM / "nterm_signalp6_built_calls.tsv"
RUNGS = ("nterm", "internal", "multipass", "type3", "cterm", "rescue")

COLS = [
    "design_id",
    "gene_symbol",
    "species",
    "candidate_set",
    "sequence_sha256",
    "length_aa",
    "sp_status",
    "native_sp_last_residue",
    "sp_evidence",
    "sp_concern",
    "sp_metrics_json",
    "uniprot_canonical_acc",
    "uniprot_canonical_length",
    "matches_uniprot_canonical",
    "uniprot_diff_json",
    "sp_flags",
]

#: Phrases in a member's lane override saying its v1 signal peptide is not a cleaved one.
NOT_CLEAVED_SP = re.compile(r"signal anchor|no signal peptide|N-terminus is not cleaved", re.I)
#: An ortholog's signal peptide named in a lane override as supporting evidence.
ORTHOLOG_SP = re.compile(r"\b(mouse|cyno|chicken|rat)\b[^.;]*signal peptide", re.I)
UNIPROT_ACC = re.compile(
    r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})(-\d+)?$"
)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def by_id(path: Path, key: str = "candidate_id") -> dict[str, dict[str, str]]:
    out: dict[str, dict[str, str]] = {}
    for r in read_tsv(path):
        if r[key] in out:
            raise SystemExit(f"{path.name}: duplicate {key} {r[key]}")
        out[r[key]] = r
    return out


# ---------------------------------------------------------------- signal peptide


def v1_sp_end(per_residue: str) -> int | None:
    """End of the leading ``S`` run of a v1 per-residue string, or None."""
    m = re.match(r"S+", per_residue or "")
    return m.end() if m else None


def construct_arms() -> dict[str, list[tuple[str, str]]]:
    """``{design_id: [(rung, design), ...]}`` from every L1 construct table."""
    arms: dict[str, list[tuple[str, str]]] = defaultdict(list)
    for rung in RUNGS:
        for r in read_tsv(HERE / f"constructs_l1_{rung}.tsv"):
            design = r.get("arm") or r.get("design") or r.get("decision") or "c_terminal"
            arms[r["design_id"]].append((rung, design))
    return arms


def status_of_arm(rung: str, design: str, has_sp: bool) -> str:
    """The sp_status one built arm implies."""
    if rung == "nterm":
        return {"native_site": "native", "igk_swap": "igk_swap"}[design]
    if rung in ("multipass", "type3", "rescue"):
        if design in ("igk_tag",):
            return "igk_res1"
        if design == "igk_sp_swap":
            return "igk_swap"
        if design == "tag_alone":
            return "none"
        raise SystemExit(f"unknown {rung} design {design!r}")
    # internal and C-terminal arms leave the native N-terminus as it is
    return "native" if has_sp else "none"


def swap_evidence(dec: dict[str, str], lane: str) -> str:
    """Short label for an Igκ swap: why the native SP was replaced, and where it ends."""
    reason = dec["reason"]
    if (
        reason.startswith("native signal peptide fails at every site")
        or "igk_after_native_fail" in reason
    ):
        why = "SignalP 6 predicts tagged native SP uncleaved"
    elif reason.startswith("weak predicted SP"):
        why = "weak-SP rule (SignalP 6 native/tagged cleavage)"
    elif "that SignalP 6.0 does not" in reason:
        why = "SignalP 6 calls no native SP; second source confirms a cleaved SP"
    elif "no_row_igk" in reason:
        why = "SignalP 6: no native row; Igκ form keeps cleavage"
    elif reason.startswith("refold of igk_swap"):
        why = "refold verdict"
    else:
        why = dec["evidence"]
    ev = dec["igk_deletion_evidence"]
    sources = []
    if ev.startswith("SignalP 6.0 native SP"):
        p = re.search(r"\(p=([0-9.]+)\)", ev)
        sources.append(f"SignalP 6.0 (p={p.group(1)})" if p else "SignalP 6.0")
    else:
        if "DeepTMHMM v2" in ev:
            sources.append("DeepTMHMM v2")
        if "human ortholog UniProt" in ev:
            sources.append("human ortholog UniProt SIGNAL")
        elif "UniProt" in ev and "SIGNAL" in ev:
            sources.append("UniProt SIGNAL")
        if "PMID" in ev:
            sources.append("literature")
    m = ORTHOLOG_SP.search(lane)
    if m:
        sources.append(f"{m.group(1).lower()} ortholog")
    end = " + ".join(sources) or ev
    if ev.startswith("re-anchored"):
        end += " (re-anchored)"
    return f"{why}; SP end: {end}"


def native_evidence(dec: dict[str, str]) -> str:
    ev = dec["evidence"]
    if "pro-peptide" in ev:
        return "SignalP 6 keeps native cleavage; UniProt pro-peptide"
    if "refold" in ev:
        return f"SignalP 6 keeps native cleavage; {ev}"
    return "SignalP 6 keeps native cleavage"


def _opt_float(v: object) -> float | None:
    if v in (None, ""):
        return None
    return round(float(v), 4)


def _opt_int(v: object) -> int | None:
    if v in (None, ""):
        return None
    return int(float(v))


def load_signalp_cache() -> dict[str, dict[str, str]]:
    """Local cache of SignalP 6 calls from pinned S3 runs, keyed by sequence SHA-256."""
    if not SIGNALP_CACHE.exists():
        return {}
    return by_id(SIGNALP_CACHE, "sequence_sha256")


def refresh_signalp_cache(manifest_text: str) -> int:
    """Pull SignalP 6 calls for active N-after-SP native & built-tagged SHAs from pinned runs."""
    sys.path.insert(0, str(NTERM))
    import nterm_overnight as no

    from deliverome_analysis.s3 import processed

    with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as fh:
        fh.write(manifest_text)
        snap = Path(fh.name)
    try:
        rows = list(csv.DictReader(io.StringIO(manifest_text), delimiter="\t"))
        members = [r for r in rows if r["library"] == ld.LIBRARY_L1 and ld.is_active(r)]
        prot = ld.designed_proteins(snap, library=ld.LIBRARY_L1)
    finally:
        snap.unlink()
    decisions = by_id(NTERM / "nterm_decisions.tsv", "design_id")
    wanted: set[str] = set()
    for r in members:
        did = ld.design_id(r)
        dec = decisions.get(did)
        if not dec or dec["arm"] not in ("native_site", "igk_swap"):
            continue
        if r["candidate_id"] in prot:
            wanted.add(topology.sequence_sha256(prot[r["candidate_id"]]["sequence"]))
        if dec.get("tagged_sha256"):
            wanted.add(dec["tagged_sha256"])

    calls = no.prior_calls(processed)
    cols = [
        "sequence_sha256",
        "prediction",
        "cleavage_site",
        "sp_probability",
        "cleavage_probability",
    ]
    out_rows = []
    for sha in sorted(wanted):
        c = calls.get(sha)
        if not c:
            continue
        out_rows.append(
            {
                "sequence_sha256": sha,
                "prediction": c.get("prediction") or "",
                "cleavage_site": "" if c.get("cleavage_site") is None else str(c["cleavage_site"]),
                "sp_probability": ""
                if c.get("sp_probability") is None
                else str(c["sp_probability"]),
                "cleavage_probability": ""
                if c.get("cleavage_probability") is None
                else str(c["cleavage_probability"]),
            }
        )
    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=cols, delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(out_rows)
    SIGNALP_CACHE.write_text(buf.getvalue())
    print(f"wrote {SIGNALP_CACHE.name}: {len(out_rows)} SignalP 6 calls")
    return 0


def sp_metrics_for(
    dec: dict[str, str] | None,
    nsha: str,
    v1_end: int | None,
    mem: dict[str, str],
    opt: dict[str, str],
    dis: dict[str, str],
    pos_by_sha: dict[str, dict[str, str]],
    sp_calls: dict[str, dict[str, str]],
) -> tuple[str, str]:
    """Return ``(sp_concern, sp_metrics_json)`` for built N-after-SP constructs."""
    if not dec or dec["arm"] not in ("native_site", "igk_swap"):
        return "", ""
    arm = dec["arm"]
    tsha = dec.get("tagged_sha256", "")
    cn = sp_calls.get(nsha) or {}
    ct = sp_calls.get(tsha) or {}
    pt = pos_by_sha.get(tsha) or {}

    dtm_cut = _opt_int(mem.get("deeptmhmm_v1_cleavage") or dis.get("deeptmhmm_v1_cleavage"))
    if dtm_cut is None:
        dtm_cut = v1_end
    pre_pred = cn.get("prediction") or mem.get("signalp_prediction") or ""
    pre_cut = _opt_int(cn.get("cleavage_site") or mem.get("signalp_cleavage"))
    pre_sp_p = _opt_float(cn.get("sp_probability") or opt.get("native_sp_probability"))
    pre_cs_p = _opt_float(
        cn.get("cleavage_probability")
        or mem.get("signalp_cleavage_probability")
        or opt.get("native_cleavage_probability")
    )
    post_pred = ct.get("prediction") or pt.get("tagged_prediction") or ""
    post_cut = _opt_int(ct.get("cleavage_site") or pt.get("tagged_cleavage"))
    post_sp_p = _opt_float(ct.get("sp_probability") or pt.get("tagged_sp_probability"))
    post_cs_p = _opt_float(
        ct.get("cleavage_probability")
        or pt.get("tagged_cleavage_probability")
        or dec.get("tagged_cleavage_probability")
    )
    spread = abs(pre_cut - dtm_cut) if (pre_cut is not None and dtm_cut is not None) else None
    delta = (
        round(post_cs_p - pre_cs_p, 4)
        if (arm == "native_site" and pre_cs_p is not None and post_cs_p is not None)
        else None
    )

    reasons: list[str] = []
    if arm == "native_site":
        if pre_pred != "SP":
            reasons.append("SignalP 6 does not call a native SP")
        if spread is not None and spread > 3:
            reasons.append(f"DeepTMHMM v1 vs SignalP 6 cut spread {spread} aa")
        if pre_cs_p is not None and pre_cs_p < 0.5:
            reasons.append(f"weak pre-tag cleavage (p={pre_cs_p:.2f})")
        if post_cs_p is not None and post_cs_p < 0.5:
            reasons.append(f"weak post-tag cleavage (p={post_cs_p:.2f})")
        if post_sp_p is not None and post_sp_p < 0.8:
            reasons.append(f"low post-tag SP prob ({post_sp_p:.2f})")
        if delta is not None and delta < -0.07:
            reasons.append(f"post-tag cleavage softened ({delta:+.2f})")
        if pre_cut is not None and post_cut is not None and post_cut != pre_cut:
            reasons.append(f"post-tag cleavage shifted ({pre_cut}->{post_cut})")
        if dec.get("weak_sp_flag"):
            reasons.append("weak-SP review flag")
    else:
        if pre_pred != "SP":
            reasons.append("SignalP 6 calls no native SP (second-source Igκ swap)")
        if spread is not None and spread > 3:
            reasons.append(f"native DeepTMHMM v1 vs SignalP 6 cut spread {spread} aa")
        if "re-anchored" in (dec.get("igk_deletion_evidence") or ""):
            reasons.append("Igκ deletion boundary re-anchored")
        if post_cs_p is not None and post_cs_p < 0.8:
            reasons.append(f"weak post-tag Igκ cleavage (p={post_cs_p:.2f})")
        if post_cut is not None and post_cut != 20:
            reasons.append(f"Igκ cleavage shifted ({post_cut})")

    concern = "concerning" if reasons else "clean"
    payload = {
        "deeptmhmm_v1_cut": dtm_cut,
        "signalp_pre_pred": pre_pred or None,
        "signalp_pre_cut": pre_cut,
        "signalp_pre_sp_prob": pre_sp_p,
        "signalp_pre_cs_prob": pre_cs_p,
        "signalp_post_pred": post_pred or None,
        "signalp_post_cut": post_cut,
        "signalp_post_sp_prob": post_sp_p,
        "signalp_post_cs_prob": post_cs_p,
        "cut_spread_aa": spread,
        "cleavage_prob_delta": delta,
        "post_leader": "igk" if arm == "igk_swap" else "native",
        "insert_after_residue": _opt_int(dec.get("native_after_residue")),
        "igk_replaces_residues": dec.get("igk_replaces_residues") or None,
        "concern": concern,
        "concern_reasons": reasons,
    }
    return concern, json.dumps(payload, separators=(",", ":"))


# ---------------------------------------------------------------- UniProt


def load_uniprot(species: str) -> dict[str, str]:
    """``{accession: sequence}`` from the pinned snapshot tables for one species."""
    out: dict[str, str] = {}
    for path in candidate_inputs.uniprot_tables(species):
        with gzip.open(path, "rt") if path.suffix == ".gz" else path.open() as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                acc, seq = row.get("Entry"), row.get("Sequence")
                if acc and seq:
                    out[acc] = seq.strip().upper()
    return out


def canonical_for(
    member: dict[str, str],
    seq: str,
    audit: dict[str, dict[str, str]],
    uniprot: dict[str, dict[str, str]],
    by_seq: dict[str, dict[str, list[str]]],
) -> tuple[str, str]:
    """(canonical accession, how it was found) or ("", reason)."""
    species = topology.normalize_species(member["species"])
    table = uniprot.get(species, {})
    key_acc = member["key"].split(":")[-1]
    if species == "human":
        acc = (audit.get(member["design_id"]) or {}).get("canonical_acc", "")
        if acc:
            acc = acc.split("-")[0]
            return (acc, "audit") if acc in table else ("", f"audit acc {acc} not in snapshot")
        if UNIPROT_ACC.match(key_acc):
            base = key_acc.split("-")[0]
            if base in table:
                return base, "member key"
    else:
        if UNIPROT_ACC.match(key_acc):
            base = key_acc.split("-")[0]
            if base in table:
                return base, "member key"
    hits = sorted({a.split("-")[0] for a in by_seq.get(species, {}).get(seq, [])})
    hits = [a for a in hits if a in table]
    if len(hits) == 1:
        return hits[0], "sequence identity"
    if len(hits) > 1:
        return "", f"sequence matches {len(hits)} records"
    return "", "no record"


def aligner() -> PairwiseAligner:
    a = PairwiseAligner()
    a.mode = "global"
    a.substitution_matrix = substitution_matrices.load("BLOSUM62")
    a.open_gap_score = -10
    a.extend_gap_score = -0.5
    return a


def diff(member: str, ref: str, al: PairwiseAligner) -> list[dict]:
    """Differences of ``member`` against ``ref`` from one optimal global alignment."""
    if member == ref:
        return []
    aln = al.align(member, ref)[0]
    m_row, r_row = aln[0], aln[1]
    cols = []  # (kind, member_pos|None, ref_pos|None)
    mi = ri = 0
    for a, b in zip(m_row, r_row, strict=True):
        mp = rp = None
        if a != "-":
            mi += 1
            mp = mi
        if b != "-":
            ri += 1
            rp = ri
        if mp and rp:
            kind = "match" if a == b else "substitution"
        elif mp:
            kind = "insertion"
        else:
            kind = "deletion"
        cols.append((kind, mp, rp, a, b))
    # leading / trailing one-sided gaps are extensions or truncations
    n = len(cols)
    first = next(i for i, c in enumerate(cols) if c[1] and c[2])
    last = n - 1 - next(i for i, c in enumerate(reversed(cols)) if c[1] and c[2])
    for i in range(first):
        cols[i] = ("n_extension" if cols[i][1] else "n_truncation", *cols[i][1:])
    for i in range(last + 1, n):
        cols[i] = ("c_extension" if cols[i][1] else "c_truncation", *cols[i][1:])

    out: list[dict] = []
    run: list = []

    def flush() -> None:
        if not run:
            return
        mps = [c[1] for c in run if c[1]]
        rps = [c[2] for c in run if c[2]]
        out.append(
            {
                "kind": run[0][0],
                "member_start": mps[0] if mps else None,
                "member_end": mps[-1] if mps else None,
                "uniprot_start": rps[0] if rps else None,
                "uniprot_end": rps[-1] if rps else None,
                "member_aa": "".join(c[3] for c in run if c[1]),
                "uniprot_aa": "".join(c[4] for c in run if c[2]),
            }
        )
        run.clear()

    for c in cols:
        if c[0] == "match" or (run and run[0][0] != c[0]):
            flush()
        if c[0] != "match":
            run.append(c)
    flush()
    return out


# ---------------------------------------------------------------- build


def _sp_call(
    r: dict,
    did: str,
    seq: str,
    v1,
    decisions: dict,
    picks: dict,
    tag_paths: dict,
    arms: dict,
    flags: list[str],
    stats: Counter,
) -> tuple[str, int | None, str]:
    """The member's signal-peptide status, last SP residue and evidence; appends to ``flags``."""
    # --- signal peptide
    v1_end = None
    v1_found = bool(seq) and did in v1.index and bool(v1.loc[did, "found"])
    if v1_found:
        v1_end = v1_sp_end(v1.loc[did, "per_residue_topology"])
    elif seq:
        flags.append("no DeepTMHMM v1 row for this sequence")
    dec = decisions.get(did)
    pick = picks.get(did)
    lane = (tag_paths.get(r["candidate_id"]) or {}).get("lane_override", "")
    built = arms.get(did, [])

    sp_end: int | None = None
    evidence = ""
    overruled = ""
    if dec and dec["arm"] == "igk_swap":
        sp_end = int(dec["igk_replaces_residues"].split("-")[-1])
        evidence = swap_evidence(dec, lane)
        if pick and pick["cleavage_site"] and int(pick["cleavage_site"]) != sp_end:
            stats["swap_end_differs_from_site_pick"] += 1
    elif dec and dec["arm"] == "native_site":
        if pick and pick["cleavage_site"]:
            sp_end = int(pick["cleavage_site"])
        else:
            sp_end = int(dec["cleavage_site"])
        evidence = native_evidence(dec)
    elif dec and dec.get("sp_review_verdict") == "signal_anchor_or_tm":
        overruled = "second-source review: signal anchor / TM"
    elif NOT_CLEAVED_SP.search(lane):
        overruled = "lane override: " + NOT_CLEAVED_SP.search(lane).group(0)
    elif any(rung in ("multipass", "type3") for rung, _ in built):
        overruled = "built at residue 1 (no cleaved SP)"
    if sp_end is None and not overruled and v1_end:
        sp_end = v1_end
        evidence = "DeepTMHMM v1"
        if dec and dec.get("sp_review_verdict") == "cleaved_sp_confirmed":
            evidence = "DeepTMHMM v1; second source confirms cleaved SP"
    if overruled and v1_end:
        flags.append(f"v1 SP 1-{v1_end} overruled ({overruled})")
    has_sp = sp_end is not None

    statuses = {status_of_arm(rung, d, has_sp) for rung, d in built}
    if dec and dec["arm"] in ("native_site", "igk_swap"):
        statuses.add(status_of_arm("nterm", dec["arm"], has_sp))
    elif dec and dec["arm"] == "routed_multipass" and not built:
        flags.append("routed to residue-1 multipass but no multipass arm built")
    if r.get("has_igk_leader", "").strip().lower() in ("true", "1"):
        statuses.add("igk_in_base_cds")
    status = ""
    if len(statuses) == 1:
        status = statuses.pop()
        if status == "igk_in_base_cds":
            status = ""
            flags.append("manifest has_igk_leader but no Igκ arm")
    elif len(statuses) > 1:
        flags.append("sources disagree: " + ",".join(sorted(statuses)))
    else:
        placement = (tag_paths.get(r["candidate_id"]) or {}).get("placement", "")
        if placement in ("internal", "C-terminal"):
            status = "native" if has_sp else "none"
            flags.append(f"no arm built yet; {placement} rung keeps the native N-terminus")
        else:
            flags.append(f"undecided: no arm built (placement {placement or 'none'})")
    if status == "igk_res1" and has_sp:
        flags.append("igk_res1 but a native SP is called")
    if status == "native" and not has_sp:
        flags.append("native but no SP called")
    if not evidence:
        if overruled:
            evidence = overruled
        elif not v1_found:
            evidence = "no DeepTMHMM v1 row"
        else:
            evidence = "DeepTMHMM v1: no SP"
    return status, sp_end, evidence


def build(manifest_text: str) -> tuple[list[dict[str, str]], dict]:
    with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as fh:
        fh.write(manifest_text)
        snap = Path(fh.name)
    try:
        rows = list(csv.DictReader(io.StringIO(manifest_text), delimiter="\t"))
        members = [r for r in rows if r["library"] == ld.LIBRARY_L1 and ld.is_active(r)]
        for r in members:
            r["design_id"] = ld.design_id(r)
        prot = ld.designed_proteins(snap, library=ld.LIBRARY_L1)
    finally:
        snap.unlink()

    seqs = {
        r["design_id"]: prot[r["candidate_id"]]["sequence"]
        for r in members
        if r["candidate_id"] in prot
    }
    no_cds = [r["design_id"] for r in members if r["candidate_id"] not in prot]

    v1 = topology.resolve(seqs, tool=topology.TOOL_V1).set_index("key")
    audit = by_id(HERE / "l1_proteoform_audit.tsv", "design_id")
    tag_paths = by_id(HERE / "l1_tag_paths.tsv")
    decisions = by_id(NTERM / "nterm_decisions.tsv", "design_id")
    picks = by_id(NTERM / "nterm_site_picks.tsv", "design_id")
    nterm_members = by_id(NTERM / "nterm_members.tsv", "design_id")
    nterm_opts = by_id(NTERM / "nterm_build_options.tsv", "design_id")
    sp_disagree = by_id(NTERM / "sp_disagree_review.tsv", "design_id")
    pos_by_sha = {
        r["tagged_sha256"]: r
        for r in read_tsv(NTERM / "nterm_positions.tsv")
        if r.get("tagged_sha256")
    }
    sp_calls = load_signalp_cache()
    arms = construct_arms()

    uniprot = {sp: load_uniprot(sp) for sp in ("human", "mouse", "cyno")}
    by_seq: dict[str, dict[str, list[str]]] = {}
    for sp, table in uniprot.items():
        idx: dict[str, list[str]] = defaultdict(list)
        for acc, s in table.items():
            idx[s].append(acc)
        by_seq[sp] = idx
    al = aligner()

    stats: Counter = Counter()
    sha_mismatch: list[str] = []
    out: list[dict[str, str]] = []
    for r in members:
        did = r["design_id"]
        seq = seqs.get(did, "")
        flags: list[str] = []
        sha = topology.sequence_sha256(seq) if seq else ""
        for src, row in (
            ("l1_tag_paths", tag_paths.get(r["candidate_id"])),
            ("audit", audit.get(did)),
        ):
            if row and row.get("sequence_sha256") and sha and row["sequence_sha256"] != sha:
                sha_mismatch.append(f"{did}:{src}")
                flags.append(f"sha differs from {src}")

        status, sp_end, evidence = _sp_call(
            r, did, seq, v1, decisions, picks, tag_paths, arms, flags, stats
        )
        v1_end = (
            v1_sp_end(v1.loc[did, "per_residue_topology"])
            if (seq and did in v1.index and bool(v1.loc[did, "found"]))
            else None
        )
        sp_concern, sp_metrics = sp_metrics_for(
            decisions.get(did),
            sha,
            v1_end,
            nterm_members.get(did) or {},
            nterm_opts.get(did) or {},
            sp_disagree.get(did) or {},
            pos_by_sha,
            sp_calls,
        )
        if sp_concern:
            stats[f"sp_concern:{sp_concern}"] += 1

        # --- UniProt
        acc, how = canonical_for(r, seq, audit, uniprot, by_seq) if seq else ("", "no CDS")
        ref = uniprot[topology.normalize_species(r["species"])].get(acc, "") if acc else ""
        stats[f"canonical_via:{how if acc else 'none'}"] += 1
        if acc:
            d = diff(seq, ref, al)
            diff_json = json.dumps(d, separators=(",", ":"))
            match = "1" if seq == ref else "0"
        else:
            diff_json = match = ""
            if how not in ("no record",):
                flags.append(f"no canonical: {how}")

        out.append(
            {
                "design_id": did,
                "gene_symbol": ld.gene_from_key(r["key"]),
                "species": r["species"],
                "candidate_set": r["set"],
                "sequence_sha256": sha,
                "length_aa": str(len(seq)) if seq else "",
                "sp_status": status,
                "native_sp_last_residue": str(sp_end) if sp_end else "",
                "sp_evidence": evidence,
                "sp_concern": sp_concern,
                "sp_metrics_json": sp_metrics,
                "uniprot_canonical_acc": acc,
                "uniprot_canonical_length": str(len(ref)) if ref else "",
                "matches_uniprot_canonical": match,
                "uniprot_diff_json": diff_json,
                "sp_flags": "; ".join(flags),
            }
        )
    summary = {
        "manifest_rows": len(rows),
        "members": len(members),
        "no_designed_cds": no_cds,
        "sha_mismatch": sha_mismatch,
        **stats,
    }
    return out, summary


def render(rows: list[dict[str, str]]) -> str:
    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=COLS, delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(rows)
    return buf.getvalue()


def spot_checks(rows: list[dict[str, str]]) -> list[str]:
    """The cases the viewer relies on; returns failures."""
    fails = []
    by_gene = defaultdict(list)
    for r in rows:
        by_gene[(r["gene_symbol"], r["species"], r["candidate_set"])].append(r)
    glp2r = by_gene[("GLP2R", "human", "library_reference")]
    if (
        not glp2r
        or glp2r[0]["sp_status"] != "igk_swap"
        or glp2r[0]["native_sp_last_residue"] != "66"
    ):
        fails.append(f"GLP2R: expected igk_swap / 66, got {glp2r[:1]}")
    elif (
        "DeepTMHMM v2" not in glp2r[0]["sp_evidence"]
        or "mouse ortholog" not in glp2r[0]["sp_evidence"]
    ):
        fails.append(f"GLP2R evidence: {glp2r[0]['sp_evidence']}")
    for gene, want in (("ADRB2", "igk_res1"), ("CD52", "native")):
        hit = by_gene.get((gene, "human", "library_reference"))
        if hit and hit[0]["sp_status"] != want:
            fails.append(f"{gene}: expected {want}, got {hit[0]['sp_status']}")
    for gene, want_concern in (
        ("CD52", "clean"),
        ("TMEM81", "concerning"),
        ("GLP2R", "concerning"),
    ):
        hit = by_gene.get((gene, "human", "library_reference"))
        if hit and hit[0]["sp_concern"] != want_concern:
            fails.append(
                f"{gene}: expected sp_concern={want_concern}, got {hit[0]['sp_concern']!r}"
            )
    for r in rows:
        if r["sp_status"] not in ("", "native", "igk_swap", "igk_res1", "none"):
            fails.append(f"{r['design_id']}: bad sp_status {r['sp_status']}")
        if r["sp_concern"] not in ("", "clean", "concerning"):
            fails.append(f"{r['design_id']}: bad sp_concern {r['sp_concern']}")
        if r["sp_metrics_json"]:
            try:
                m = json.loads(r["sp_metrics_json"])
                if not isinstance(m, dict) or m.get("concern") != r["sp_concern"]:
                    fails.append(f"{r['design_id']}: sp_metrics_json concern mismatch")
            except json.JSONDecodeError:
                fails.append(f"{r['design_id']}: bad sp_metrics_json")
        if (
            r["uniprot_canonical_acc"]
            and r["matches_uniprot_canonical"] == "1"
            and r["uniprot_diff_json"] != "[]"
        ):
            fails.append(f"{r['design_id']}: identical but diff not []")
    return fails


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument(
        "--check", action="store_true", help="rebuild and compare with the committed TSV"
    )
    ap.add_argument(
        "--refresh-signalp-cache",
        action="store_true",
        help="refresh reviewer-package/nterm/nterm_signalp6_built_calls.tsv from pinned S3 runs",
    )
    args = ap.parse_args()

    text = MANIFEST.read_text()
    if args.refresh_signalp_cache:
        return refresh_signalp_cache(text)
    digest = hashlib.sha256(text.encode()).hexdigest()[:12]
    rows, summary = build(text)
    now = MANIFEST.read_text()
    n_then = text.count("\n")
    n_now = now.count("\n")
    if n_now != n_then:
        raise SystemExit(
            f"manifest row count changed during the run ({n_then} -> {n_now} lines); rerun"
        )
    if hashlib.sha256(now.encode()).hexdigest()[:12] != digest:
        print(
            "note: manifest content changed during the run (row count unchanged)", file=sys.stderr
        )

    fails = spot_checks(rows)
    for f in fails:
        print("FAIL", f, file=sys.stderr)
    status = Counter(r["sp_status"] or "(blank)" for r in rows)
    differ = sum(r["matches_uniprot_canonical"] == "0" for r in rows)
    no_can = sum(not r["uniprot_canonical_acc"] for r in rows)
    flagged = [r for r in rows if r["sp_flags"]]
    print(f"manifest {digest}: {summary['members']} active L1 members")
    print("sp_status:", dict(status))
    print(f"differ from UniProt canonical: {differ}; no canonical: {no_can}")
    print(f"flagged: {len(flagged)}; sha mismatches: {len(summary['sha_mismatch'])}")
    print({k: v for k, v in summary.items() if k not in ("no_designed_cds", "sha_mismatch")})

    text_out = render(rows)
    if args.check:
        if not OUT.exists() or OUT.read_text() != text_out:
            print(f"FAIL {OUT.name} is stale; rerun without --check", file=sys.stderr)
            return 1
        return 1 if fails else 0
    OUT.write_text(text_out)
    print(f"wrote {OUT.name}")
    return 1 if fails else 0


if __name__ == "__main__":
    raise SystemExit(main())
