"""Membrane-side gate for internal (loop) tag sites in the L1 library.

Every internal tag site proposed for an L1 member -- the S5E13 v4 ranker's picks, the
built panel loop arms, and the tag-site API's deterministic ``disorder`` and
``surface_loop`` lanes -- is placed on the member's own structure model, in a membrane
fitted to that model (``deliverome_analysis.membrane_side``), and given a verdict from
the junction's height above the mid-plane. Thresholds come from ``calibration.json``,
which ``calibrate.py`` derives from TMDET placements, UniProt experimental topology and
the depths of published working tag sites; none is set here.

DeepTMHMM v1 stays the topology of record. The gate is a separate layer on top of it: it
never relabels a residue, it only refuses a site the member's model puts on the wrong
side of, or inside, the membrane. Its second output, the model-corrected extracellular
mask, is the same evidence expressed per residue, for re-ranking a member whose picks all
fail.

Writes, in this folder:

* ``l1_internal_pick_gate.tsv``   -- one row per (member, pick): depth, pLDDT, verdict.
* ``l1_internal_pick_changes.tsv`` -- the shipped picks the gate would change, with the
  proposed action. Proposed only: the manifests and cassette tables are not touched.
* ``l1_internal_member_gate.tsv``  -- one row per internal member: model, placement,
  side signals, and its gate status (scorable or why not).
* ``l1_extracellular_mask.tsv``    -- per member, DeepTMHMM v1 ``O`` segments and the
  parts of them the model places outside the membrane.
* ``gate_summary.json``            -- the counts the README quotes.

    uv run --locked python membrane_side_gate.py
"""

from __future__ import annotations

import functools
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import gate_inputs as gi  # noqa: E402
from placements import member_placements  # noqa: E402

from deliverome_analysis import membrane_side as ms  # noqa: E402

#: Who implemented this gate, written into every output row (Becca's 2026-11-10 rule:
#: "for any manual fixes please note the agent who made the fix / adjustment and reasoning").
IMPLEMENTED_BY = "l1-membrane-side-gate session (Claude Code), PR #291"
#: Who decided the verdict rule. Becca set the requirement on 2026-11-10 (no shipped pick on
#: the cytoplasmic side or in the membrane core; every cutoff justified from evidence); the
#: cutoffs themselves are derived by calibrate.py, not chosen by hand.
DECIDED_BY = "Becca, 2026-11-10 (requirement); cutoffs derived by calibrate.py"

#: Two picks closer than this are one decision (S5E13 v4's MIN_SEPARATION).
MIN_SEPARATION = 25
#: An extracellular run shorter than this is not a separate place to tag (v4's MIN_REGION_LEN).
MIN_REGION_LEN = 10


@functools.lru_cache(maxsize=1)
def _calibration() -> dict:
    return json.loads((HERE / "calibration.json").read_text())


def thresholds() -> ms.GateThresholds:
    return ms.GateThresholds(**_calibration()["thresholds"])


#: design_id -> candidate_id for the members of this run; filled by main().
CID_OF_DESIGN: dict[str, str] = {}


def uniprot_inverted(cid_of_design: dict[str, str] | None = None) -> set[str]:
    """Members (candidate ids) whose face UniProt experimental topology (ECO:0000269)
    contradicts. Recorded by design id in calibration.json, which survives renumbering."""
    mapping = CID_OF_DESIGN if cid_of_design is None else cid_of_design
    return {mapping[d] for d in _calibration()["uniprot_inverted_design_ids"] if d in mapping}


#: Manual exclusions of partner-assembly interfaces, one row per gene, with who decided,
#: who implemented and why (`manual_exclusions.tsv`). Decided by Becca on 2026-11-11 for
#: the beta-subunit docking loop of the P-type ATPase alpha subunits; the extent (the whole
#: H7-H8 loop) was chosen by the l1-membrane-side-gate session (PR #291) from solved
#: beta-contact structures. The region is recomputed per member from the rule, so it
#: follows the member's own numbering.
EXCLUSIONS_TSV = HERE / "manual_exclusions.tsv"


def manual_exclusions(members: pd.DataFrame) -> dict[str, dict]:
    """candidate_id -> {start, end, ...attribution} for each member a manual exclusion hits.

    The excluded region is the DeepTMHMM v1 O segment that contains the gene's motif in
    the member's own designed protein. A member whose protein lacks the motif, or whose
    motif is not in a v1 O segment, gets no exclusion and is reported by the caller.
    """
    ex = pd.read_csv(EXCLUSIONS_TSV, sep="\t").set_index("gene_symbol")
    out = {}
    for m in members.itertuples():
        if m.gene_symbol not in ex.index:
            continue
        e = ex.loc[m.gene_symbol]
        at = m.protein.find(e["motif"])
        seg = [(a, b) for a, b in ms.segments(m.topology, "O") if at >= 0 and a <= at + 1 <= b]
        if not seg:
            out[m.candidate_id] = dict(start=None, end=None, **e.to_dict())
            continue
        a, b = seg[0]
        out[m.candidate_id] = dict(start=a, end=b, motif_at=at + 1, **e.to_dict())
    return out


#: candidate_id -> whether DeepTMHMM v2 independently calls both termini on the sides the
#: (possibly overridden) v1 string does. Filled by main() from the members table.
V2_AGREES: dict[str, bool | None] = {}


def _member_state(
    members: pd.DataFrame,
    cid_of_design: dict[str, str] | None = None,
    v2_agrees: dict[str, bool | None] | None = None,
) -> tuple[dict[str, str], dict[str, bool | None]]:
    """Resolve explicit or module-level member orientation state, deriving from ``members`` if unset."""
    if cid_of_design is None:
        cid_of_design = CID_OF_DESIGN or (
            dict(zip(members["design_id"], members["candidate_id"], strict=True))
            if {"design_id", "candidate_id"}.issubset(members.columns)
            else {}
        )
    if v2_agrees is None:
        v2_agrees = V2_AGREES or (
            dict(zip(members["candidate_id"], members["v2_termini_agree"], strict=True))
            if {"candidate_id", "v2_termini_agree"}.issubset(members.columns)
            else {}
        )
    return cid_of_design, v2_agrees


def orientation_ok(
    cid: str,
    p: dict,
    th,
    inverted: set[str],
    v2_agrees: dict[str, bool | None] | None = None,
) -> bool:
    """The chosen face is backed by the topology vote and not contradicted.

    A face fails when experimental UniProt topology inverts it, or when the topology vote is
    weak. An excess of flank Lys/Arg on the extracellular face contradicts it only when
    DeepTMHMM v2 does not independently confirm v1's termini: the positive-inside rule is
    there to catch a v1 inversion, and some proteins carry genuinely K/R-rich extracellular
    loops (STEAP3, whose v1, v2 and UniProt topologies agree). Decided 2026-11-12 by the
    membrane-side-gate session (delegated by Becca 2026-11-12, "use your best judgment");
    implemented by the l1-membrane-side-gate session, PR #291. MFSD2A, which this rule was
    built to catch, is still caught: v2 disagrees with v1 there.
    """
    if cid in inverted:
        return False
    if p["topo_for"] < th.min_topology_ratio * max(p["topo_against"], 1):
        return False
    kr_contradicts = (p["kr_out"] - p["kr_in"]) > th.max_kr_excess_outside
    agrees = V2_AGREES if v2_agrees is None else v2_agrees
    return not (kr_contradicts and agrees.get(cid) is not True)


def junction(z: np.ndarray, n: int) -> float:
    """Height of an insert-after-``n`` junction (1-based): mean of residues n and n+1."""
    vals = [z[n - 1]] + ([z[n]] if n < len(z) else [])
    vals = [v for v in vals if not np.isnan(v)]
    return float(np.mean(vals)) if vals else float("nan")


# --- picks --------------------------------------------------------------------------


def collect_picks(members: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """Every proposed internal site, carried onto the member's own residue numbering."""
    by_cid = members.set_index("candidate_id")
    internal = set(members.loc[members["placement"] == "internal", "candidate_id"])
    rows, notes = [], {}

    v4 = gi.v4_picks()
    notes["v4_rows"] = len(v4)
    notes["v4_rows_member_no_longer_internal"] = int((~v4["candidate_id"].isin(internal)).sum())
    for r in v4[v4["candidate_id"].isin(internal)].itertuples():
        if r.s5e13_sequence_sha256 != by_cid.loc[r.candidate_id, "sequence_sha256"]:
            # Numbered on a protein this member no longer has: not scored on the new one.
            rows.append(
                dict(
                    candidate_id=r.candidate_id,
                    source="s5e13_v4",
                    site_id=f"{r.design_id}-v4-r{r.rank}",
                    rank=int(r.rank),
                    shipped=1,
                    api_residue=int(r.insert_after_residue),
                    numbering_check="PROTEIN_CHANGED_SINCE_S5E13",
                )
            )
            continue
        rows.append(
            dict(
                candidate_id=r.candidate_id,
                source="s5e13_v4",
                site_id=f"{r.design_id}-v4-r{r.rank}",
                rank=int(r.rank),
                shipped=1,
                insert_after_residue=int(r.insert_after_residue),
                numbering_check="v1_O_check",
            )
        )

    # Built loop arms: residue identity on the designed protein is checked, not assumed.
    tp = gi.tag_paths()
    ref = tp[tp["candidate_set"] == "library_reference"].set_index("gene_symbol")["candidate_id"]
    for r in gi.loop_arms().itertuples():
        cid = ref.get(r.gene)
        if cid is None or cid not in by_cid.index:
            notes.setdefault("loop_arms_without_member", []).append(r.gene)
            continue
        aa = by_cid.loc[cid, "protein"]
        n = int(r.insert_after_residue)
        ok = aa[n - 1] == r.residue_before and (n >= len(aa) or aa[n] == r.residue_after)
        rows.append(
            dict(
                candidate_id=cid,
                source="built_loop_arm",
                site_id=f"{r.gene}-{r.variant}",
                rank=np.nan,
                shipped=1,
                insert_after_residue=n,
                numbering_check="residue_identity_ok" if ok else "RESIDUE_MISMATCH",
            )
        )

    # API deterministic lanes, numbered on the human canonical; carried by sequence.
    for cid in sorted(internal):
        m = by_cid.loc[cid]
        for s in gi.api_sites(m["gene_symbol"]):
            acc = s["api_uniprot_acc"]
            canon = gi.uniprot_sequence(acc) if acc else None
            n_api = int(s["insert_after_residue"])
            base = dict(
                candidate_id=cid,
                source=f"api_{s['det_path']}",
                site_id=f"{s['site_id']}@{cid}",
                rank=np.nan,
                shipped=0,
                api_residue=n_api,
            )
            if canon is None or not 1 <= n_api <= len(canon):
                rows.append(base | dict(numbering_check="no_canonical_sequence"))
                continue
            if canon[n_api - 1] != s.get("residue_before"):
                rows.append(base | dict(numbering_check="API_RESIDUE_MISMATCH"))
                continue
            cmap, _ = ms.map_to_model(canon, m["protein"])
            j = cmap[n_api - 1]
            if j is None:
                rows.append(base | dict(numbering_check="not_in_member_protein"))
                continue
            rows.append(
                base | dict(insert_after_residue=j + 1, numbering_check="residue_identity_ok")
            )
    picks = pd.DataFrame(rows)
    notes["api_fetched"] = gi.api_fetch_date()
    return picks, notes


# --- verdicts -----------------------------------------------------------------------


def score_picks(
    picks: pd.DataFrame,
    members: pd.DataFrame,
    placed: dict,
    th,
    *,
    cid_of_design: dict[str, str] | None = None,
    v2_agrees: dict[str, bool | None] | None = None,
) -> pd.DataFrame:
    exclusions = manual_exclusions(members)
    cid_map, v2_map = _member_state(members, cid_of_design, v2_agrees)
    inverted = uniprot_inverted(cid_map)
    by_cid = members.set_index("candidate_id")
    out = []
    for r in picks.to_dict("records"):
        m = by_cid.loc[r["candidate_id"]]
        p = placed[r["candidate_id"]]
        n = r.get("insert_after_residue")
        r |= dict(
            design_id=m["design_id"],
            gene_symbol=m["gene_symbol"],
            candidate_set=m["candidate_set"],
            species=m["species"],
            library="l1",
            topology_tool=m["topology_layer"],
            topology_version=m["topology_version"],
            topology_cohort=m["topology_cohort"],
        )
        if n is None or (isinstance(n, float) and np.isnan(n)):
            out.append(r | dict(verdict="unscorable", reason=r["numbering_check"]))
            continue
        n = int(r["insert_after_residue"])
        r["insert_after_residue"] = n
        r["residue_before"] = m["protein"][n - 1]
        r["residue_after"] = m["protein"][n] if n < len(m["protein"]) else ""
        r["v1_state"] = m["topology"][n - 1] + (m["topology"][n] if n < len(m["topology"]) else "")
        if r["numbering_check"] == "v1_O_check":
            r["numbering_check"] = "v1_O" if m["topology"][n - 1] == "O" else "NOT_V1_O"
        if p["status"] != "placed":
            out.append(r | dict(verdict="unscorable", reason=p["status"]))
            continue
        z = junction(p["z"], n)
        pl = float(np.nanmean(p["plddt"][n - 1 : n + 1]))
        orient_ok = orientation_ok(r["candidate_id"], p, th, inverted, v2_map)
        v = ms.site_verdict(z, pl, orient_ok, th)
        row = r | dict(
            model_source=p["model_source"],
            model_match=p["model_match"],
            height_A=round(z, 1),
            plddt=round(pl, 1),
            membrane_verdict=v,
            verdict=v,
            reason=f"junction at {z:+.1f} A, pLDDT {pl:.0f}; thresholds in calibration.json",
        )
        e = exclusions.get(r["candidate_id"])
        if e and e["start"] is not None and e["start"] <= n <= e["end"]:
            row |= dict(
                verdict="fail_partner_interface",
                reason=f"manual exclusion: {e['motif']} loop {e['start']}-{e['end']}; {e['basis']}",
                decided_by=e["decided_by"],
                implemented_by=e["implemented_by"],
            )
        out.append(row)
    cols = [
        "candidate_id", "design_id", "gene_symbol", "candidate_set", "species", "library",
        "source",
        "site_id", "rank", "shipped", "insert_after_residue", "residue_before",
        "residue_after", "api_residue", "numbering_check", "v1_state", "topology_tool",
        "topology_version", "topology_cohort", "model_source", "model_match", "height_A",
        "plddt", "membrane_verdict", "verdict", "reason", "decided_by", "implemented_by",
    ]  # fmt: skip
    df = pd.DataFrame(out)
    for c, default in (("decided_by", DECIDED_BY), ("implemented_by", IMPLEMENTED_BY)):
        df[c] = df[c].fillna(default) if c in df.columns else default
    return df[[c for c in cols if c in df.columns]].sort_values(
        ["candidate_id", "shipped", "source", "rank", "insert_after_residue"],
        ascending=[True, False, True, True, True],
    )


# --- member table and corrected mask ------------------------------------------------


def member_table(
    members: pd.DataFrame,
    placed: dict,
    th,
    *,
    cid_of_design: dict[str, str] | None = None,
    v2_agrees: dict[str, bool | None] | None = None,
) -> pd.DataFrame:
    cid_map, v2_map = _member_state(members, cid_of_design, v2_agrees)
    inverted = uniprot_inverted(cid_map)
    rows = []
    for m in members.itertuples():
        p = placed[m.candidate_id]
        base = dict(
            candidate_id=m.candidate_id,
            design_id=m.design_id,
            gene_symbol=m.gene_symbol,
            candidate_set=m.candidate_set,
            species=m.species,
            placement=m.placement,
            topology_tool=m.topology_layer,
            topology_version=m.topology_version,
            topology_cohort=m.topology_cohort,
            status=p["status"],
            model_source=p.get("model_source", ""),
            model_match=p.get("model_match", ""),
            placement_source=p.get("placement_source", ""),
            tmdet_rmsd=p.get("tmdet_rmsd", np.nan),
        )
        if p["status"] == "placed":
            base |= dict(
                topology_votes_for=p["topo_for"],
                topology_votes_against=p["topo_against"],
                kr_flank_extracellular=p["kr_out"],
                kr_flank_cytoplasmic=p["kr_in"],
                orientation_supported=orientation_ok(m.candidate_id, p, th, inverted, v2_map),
                uniprot_experimental_inverted=m.candidate_id in inverted,
            )
        rows.append(base)
    return pd.DataFrame(rows)


def corrected_mask(members: pd.DataFrame, placed: dict, th) -> pd.DataFrame:
    """v1 ``O`` segments, and the runs inside them the model places outside the membrane."""
    rows = []
    for m in members.itertuples():
        p = placed[m.candidate_id]
        if m.placement != "internal" or p["status"] != "placed":
            continue
        for a, b in ms.segments(m.topology, "O"):
            z = p["z"][a - 1 : b]
            keep = "".join("o" if (not np.isnan(v) and v >= th.outside) else "." for v in z)
            runs = [(a + s - 1, a + e - 1) for s, e in ms.segments(keep, "o")]
            rows.append(
                dict(
                    candidate_id=m.candidate_id,
                    design_id=m.design_id,
                    gene_symbol=m.gene_symbol,
                    v1_segment=f"{a}-{b}",
                    v1_length=b - a + 1,
                    min_height_A=round(float(np.nanmin(z)), 1) if np.isfinite(z).any() else np.nan,
                    max_height_A=round(float(np.nanmax(z)), 1) if np.isfinite(z).any() else np.nan,
                    residues_outside=sum(e - s + 1 for s, e in runs),
                    model_outside_runs=";".join(f"{s}-{e}" for s, e in runs),
                    longest_outside_run=max((e - s + 1 for s, e in runs), default=0),
                )
            )
    return pd.DataFrame(rows)


#: Per-residue mask alphabet, one character per residue of the member's designed protein.
MASK_ALPHABET = {
    "E": "eligible: DeepTMHMM v1 O, and the model puts it at or above `outside`",
    "r": "v1 O in the review band (outside - margin <= z < outside)",
    "m": "v1 O inside the membrane core (-core <= z < outside - margin)",
    "c": "v1 O on the cytoplasmic side (z < -core)",
    "u": "v1 O with no model coordinate (residue not in the model)",
    "i": "v1 I or M, but the model puts it at or above `outside` at pLDDT >= 70 (review, not eligible)",
    ".": "not v1 O, and not outside on the model",
    "x": "v1 O, excluded by a manual partner-interface exclusion (manual_exclusions.tsv)",
}


def residue_mask(
    members: pd.DataFrame,
    placed: dict,
    th,
    *,
    cid_of_design: dict[str, str] | None = None,
    v2_agrees: dict[str, bool | None] | None = None,
) -> pd.DataFrame:
    exclusions = manual_exclusions(members)
    cid_map, v2_map = _member_state(members, cid_of_design, v2_agrees)
    """One row per internal member: a per-residue extracellular mask for site rankers.

    Keyed by ``sequence_sha256`` of the designed protein, the join a ranker must use. A
    member with ``usable`` = 0 (not placed, or its extracellular face is unconfirmed) has
    no mask a ranker may rank on; its sites must be held, not ranked on v1 alone. An
    insertion after residue n is eligible when residues n and n+1 are both ``E``; both at
    or above ``outside`` means the junction passes this gate. ``junction_mask`` is the
    rule rankers should use: character n says whether "insert after n" passes this gate
    (``E``), sits in the review band (``r``), the membrane (``m``) or the cytoplasm (``c``),
    or lacks a coordinate (``u``); ``.`` where v1 does not call residue n ``O``.
    """
    inverted = uniprot_inverted(cid_map)
    rows = []
    for m in members.itertuples():
        if m.placement != "internal":
            continue
        p = placed[m.candidate_id]
        base = dict(
            candidate_id=m.candidate_id,
            design_id=m.design_id,
            gene_symbol=m.gene_symbol,
            candidate_set=m.candidate_set,
            species=m.species,
            sequence_sha256=m.sequence_sha256,
            length_aa=len(m.protein),
            topology_tool=m.topology_layer,
            topology_version=m.topology_version,
            topology_cohort=m.topology_cohort,
            status=p["status"],
            placement_source=p.get("placement_source", ""),
        )
        if p["status"] != "placed":
            rows.append(
                base | dict(orientation_supported=False, usable=0, mask="", junction_mask="")
            )
            continue
        ok = orientation_ok(m.candidate_id, p, th, inverted, v2_map)
        z, pl, t = p["z"], p["plddt"], m.topology
        chars = []
        for i, state in enumerate(t):
            zi = z[i]
            if state == "O":
                if np.isnan(zi):
                    chars.append("u")
                elif zi >= th.outside:
                    chars.append("E")
                elif zi >= th.outside - th.margin:
                    chars.append("r")
                elif zi >= -th.core:
                    chars.append("m")
                else:
                    chars.append("c")
            elif not np.isnan(zi) and zi >= th.outside and pl[i] >= ms.FIT_MIN_PLDDT:
                chars.append("i")
            else:
                chars.append(".")
        e = exclusions.get(m.candidate_id)
        if e and e["start"] is not None:
            for i in range(e["start"] - 1, e["end"]):
                chars[i] = "x"
        mask = "".join(chars)
        assert len(mask) == len(m.protein)
        # Junction n (insert after residue n) judged exactly as the gate judges a site:
        # the mean height of residues n and n+1 against the calibrated cut-offs, for a
        # residue n that v1 calls O. This is the string a ranker should filter on, so a
        # ranker's eligible site and this gate's pass are the same thing by construction.
        jchars = []
        for n in range(1, len(t) + 1):
            if t[n - 1] != "O":
                jchars.append(".")
                continue
            zj = junction(z, n)
            if np.isnan(zj):
                jchars.append("u")
            elif zj >= th.outside:
                jchars.append("E")
            elif zj >= th.outside - th.margin:
                jchars.append("r")
            elif zj >= -th.core:
                jchars.append("m")
            else:
                jchars.append("c")
        if e and e["start"] is not None:
            for i in range(e["start"] - 1, e["end"]):
                if jchars[i] != ".":
                    jchars[i] = "x"
        junction_mask = "".join(jchars)
        rows.append(
            base
            | dict(
                orientation_supported=ok,
                usable=int(ok),
                n_v1_O=t.count("O"),
                n_eligible=mask.count("E"),
                n_removed=sum(mask.count(c) for c in "mcu"),
                n_review=mask.count("r"),
                n_i_outside=mask.count("i"),
                n_eligible_junctions=junction_mask.count("E"),
                mask=mask,
                junction_mask=junction_mask,
                decided_by=DECIDED_BY,
                implemented_by=IMPLEMENTED_BY,
            )
        )
    return pd.DataFrame(rows)


# --- proposed changes ---------------------------------------------------------------


#: A v1 O segment with at most this fraction of its residues at or above `outside` "barely
#: leaves the membrane". Set by the l1-membrane-side-gate session (PR #291), 2026-11-11,
#: from the observed gap among the membrane-core failures: none has a segment between
#: 0.25 (MIP) and 0.39 (KCNJ1, RNFT2). Not derived by code.
BARELY_OUTSIDE_FRACTION = 0.25


def failure_cause(topology: str, z: np.ndarray, n: int, th) -> dict:
    """Why a failing site fails, from the v1 O segment holding residue n, on the model."""
    seg = [(a, b) for a, b in ms.segments(topology, "O") if a <= n <= b]
    if not seg:
        return dict(cause="site not in a v1 O segment")
    a, b = seg[0]
    zz = z[a - 1 : b]
    out = [a + i for i, v in enumerate(zz) if not np.isnan(v) and v >= th.outside]
    cyto = int(np.sum(zz < -th.core))
    frac = len(out) / (b - a + 1)
    if cyto >= 3 and out:
        cause = "v1 segment spans both faces (helices missed)"
    elif cyto >= 3:
        cause = "v1 segment on the wrong side"
    elif frac <= BARELY_OUTSIDE_FRACTION:
        cause = "v1 segment barely leaves the membrane"
    else:
        cause = "site at the base of a loop that reaches outside"
    return dict(
        cause=cause,
        v1_segment=f"{a}-{b}",
        segment_fraction_outside=round(frac, 2),
        segment_residues_cytoplasmic=cyto,
        residues_to_nearest_outside=min((abs(n - i) for i in out), default=np.nan),
    )


def proposed_changes(
    gate: pd.DataFrame, mask: pd.DataFrame, members: pd.DataFrame, placed: dict, th
) -> pd.DataFrame:
    """For every shipped pick that does not pass, what would replace it.

    A failing pick gets a proposed action; a pick for review or on hold gets "review" or
    "hold", with ``if_refused`` saying what would follow if it were refused. The other
    shipped picks of the member are kept if they pass. A replacement is taken,
    in order, from the member's API deterministic sites that pass and sit at least
    MIN_SEPARATION from every kept pick (deepest first); otherwise the member needs a
    re-rank on its corrected extracellular mask, whose runs are listed.
    """
    by_cid = members.set_index("candidate_id")
    rows = []
    for cid, g in gate.groupby("candidate_id"):
        shipped = g[g["shipped"] == 1]
        bad = shipped[shipped["verdict"] != "pass"]
        if bad.empty:
            continue
        kept = shipped[shipped["verdict"] == "pass"]
        taken = list(kept["insert_after_residue"].astype(int))
        alts = g[(g["shipped"] == 0) & (g["verdict"] == "pass")].sort_values(
            "height_A", ascending=False
        )
        mrow = mask[mask["candidate_id"] == cid]
        runs = ";".join(x for x in mrow["model_outside_runs"] if x)
        for b in bad.itertuples():
            alt = None
            for a in alts.itertuples():
                if all(abs(int(a.insert_after_residue) - t) >= MIN_SEPARATION for t in taken):
                    alt = a
                    break
            if b.verdict == "unscorable":
                fallback = f"hold: member not placed ({b.reason})"
            elif not kept.empty and alt is None:
                fallback = "drop; the member keeps its other passing pick(s)"
            elif alt is not None:
                fallback = (
                    f"replace with API site {alt.site_id.split('@')[0]}, "
                    f"insert after member residue {int(alt.insert_after_residue)}"
                )
            else:
                fallback = "re-rank on the corrected extracellular mask"
            # Only a fail is a change; a review or hold is a question, and the fallback
            # is what to do if the answer is "refuse the pick".
            is_fail = b.verdict.startswith("fail_")
            action = fallback if is_fail or b.verdict == "unscorable" else b.verdict.split("_")[0]
            if is_fail and alt is not None:
                taken.append(int(alt.insert_after_residue))
            rows.append(
                dict(
                    candidate_id=cid,
                    design_id=b.design_id,
                    gene_symbol=b.gene_symbol,
                    candidate_set=b.candidate_set,
                    source=b.source,
                    site_id=b.site_id,
                    rank=b.rank,
                    insert_after_residue=b.insert_after_residue,
                    height_A=b.height_A,
                    plddt=b.plddt,
                    verdict=b.verdict,
                    other_shipped_picks_passing=len(kept),
                    proposed_action=action,
                    if_refused=fallback,
                    replacement_insert_after=(
                        int(alt.insert_after_residue) if alt is not None else np.nan
                    ),
                    replacement_height_A=(alt.height_A if alt is not None else np.nan),
                    corrected_mask_runs=runs,
                    **(
                        dict(cause="manual partner-interface exclusion")
                        if b.verdict == "fail_partner_interface"
                        else failure_cause(
                            by_cid.loc[cid, "topology"],
                            placed[cid]["z"],
                            int(b.insert_after_residue),
                            th,
                        )
                        if b.verdict.startswith("fail_")
                        else {}
                    ),
                    # Proposed, not applied: Becca approves changes to shipped picks.
                    decided_by="pending Becca's approval",
                    implemented_by=f"proposed by {IMPLEMENTED_BY}; not applied",
                    reason=f"{b.verdict}: {b.reason}" if b.reason else b.verdict,
                )
            )
    return pd.DataFrame(rows)


def main(argv: list[str] | None = None) -> int:
    import argparse

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=HERE, help="where to write the tables")
    ap.add_argument("--panel-commit", help="read members at this commit instead of the pin")
    ap.add_argument("--cds-key", help="designed-CDS record to read instead of the pin")
    ap.add_argument("--snapshot-key", help="DeepTMHMM topology snapshot (external://) to read")
    ap.add_argument("--snapshot-sha256", help="SHA-256 the snapshot must have")
    args = ap.parse_args(argv)
    out = args.out_dir
    out.mkdir(parents=True, exist_ok=True)
    pins = gi.configure(
        panel_commit=args.panel_commit,
        cds_key=args.cds_key,
        snapshot_key=args.snapshot_key,
        snapshot_sha256=args.snapshot_sha256,
    )
    th = thresholds()
    tp = gi.tag_paths()
    ref = tp[tp["candidate_set"] == "library_reference"].set_index("gene_symbol")["candidate_id"]
    extra = {ref[g] for g in gi.loop_arms()["gene"] if g in ref.index}
    members = gi.internal_members(extra)
    placed = member_placements(members)
    V2_AGREES.update(dict(zip(members["candidate_id"], members["v2_termini_agree"], strict=True)))
    CID_OF_DESIGN.update(dict(zip(members["design_id"], members["candidate_id"], strict=True)))

    picks, notes = collect_picks(members)
    gate = score_picks(picks, members, placed, th)
    mtab = member_table(members, placed, th)
    mask = corrected_mask(members, placed, th)
    changes = proposed_changes(gate, mask, members, placed, th)

    gate.to_csv(out / "l1_internal_pick_gate.tsv", sep="\t", index=False)
    mtab.to_csv(out / "l1_internal_member_gate.tsv", sep="\t", index=False)
    mask.to_csv(out / "l1_extracellular_mask.tsv", sep="\t", index=False)
    rmask = residue_mask(members, placed, th)
    rmask.to_csv(out / "l1_extracellular_mask_residues.tsv", sep="\t", index=False)
    (out / "l1_extracellular_mask_residues.alphabet.json").write_text(
        json.dumps(MASK_ALPHABET, indent=1) + "\n"
    )
    changes.to_csv(out / "l1_internal_pick_changes.tsv", sep="\t", index=False)

    internal = mtab[mtab["placement"] == "internal"]
    shipped = gate[gate["shipped"] == 1]
    summary = dict(
        library="l1",
        topology_tool="deeptmhmm_v1",
        **pins,
        panel_commit_resolved=gi.PANEL_COMMIT,
        s5e13_commit=gi.S5E13_COMMIT,
        topology_snapshot=gi.TOPOLOGY_SNAPSHOT_KEY,
        thresholds=th.__dict__,
        internal_members=len(internal),
        internal_members_by_status=internal["status"].value_counts().to_dict(),
        internal_members_by_placement_source=internal["placement_source"]
        .fillna("none")
        .replace("", "none")
        .value_counts()
        .to_dict(),
        internal_members_by_set=internal["candidate_set"].value_counts().to_dict(),
        internal_members_unscorable=internal.loc[
            internal["status"] != "placed", "candidate_id"
        ].tolist(),
        internal_members_with_v4_pick=int(
            gate.loc[gate["source"] == "s5e13_v4", "candidate_id"].nunique()
        ),
        orientation_unsupported=internal.loc[
            internal["orientation_supported"].eq(False), "candidate_id"
        ].tolist(),
        picks_by_source_verdict={
            f"{s}|{v}": int(n) for (s, v), n in gate.groupby(["source", "verdict"]).size().items()
        },
        shipped_picks=len(shipped),
        shipped_by_verdict=shipped["verdict"].value_counts().to_dict(),
        shipped_rank1_not_pass=shipped[(shipped["rank"] == 1) & (shipped["verdict"] != "pass")][
            ["candidate_id", "gene_symbol", "insert_after_residue", "height_A", "verdict"]
        ].to_dict("records"),
        changes_by_action=changes["proposed_action"].str.split(" ").str[0].value_counts().to_dict()
        if len(changes)
        else {},
        notes=notes,
    )
    (out / "gate_summary.json").write_text(json.dumps(summary, indent=1, default=str) + "\n")
    print(
        json.dumps(
            {
                k: summary[k]
                for k in ("internal_members_by_status", "shipped_by_verdict", "changes_by_action")
            },
            indent=1,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
