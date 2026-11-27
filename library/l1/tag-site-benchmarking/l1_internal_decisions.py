#!/usr/bin/env python3
"""One decided internal-rung tag placement per L1 member, with reason and the rejected alternative.

Becca delegated protein-level calls on 2026-11-12 ("use your best judgment to complete a new
and improved evidence based l1 completely by the morning"). Every row here is such a call:
`decided_by` names this session, never Becca, and each row keeps the alternative so it can be
reversed. One row per internal-rung member of the adopted panel (571).

Rules, in order:

    no usable mask, fewer than two confident helices  (ATRN, EPHA6, NRG2, RGMB)
        -> route_to_nterm_after_sp: no internal site. Single-pass or GPI-anchored proteins
           have no two-sided loop to tag; the N-terminus after the signal peptide is
           extracellular, and rung moves belong to the manifest writer.
    no usable mask, orientation unconfirmed, no eligible junction  (TMEM169, TMEM174)
        -> snorkel: the v1 extracellular exposure is 2 and 12 residues and the gate rejects
           all of it. A snorkel tag is built for multipass proteins with small loops
           (Brown et al. 2013, PMID 24023844).
    no usable mask, orientation unconfirmed, eligible junctions present
        -> internal_loop_orientation_unconfirmed: designed on the mask as computed, since v1
           is the topology of record and the gate only could not confirm the face. SLC44A1
           and TMC2 lose most v1 O residues to the gate, a sign v1 may be inverted, so they
           are flagged for a rung check.
    designed, including snorkel candidates
        -> internal_loop at the v5 rank-1 pick. Snorkel candidates stay on the loop: by
           construction a quarter of loops that hosted a working published tag score below
           the 25th percentile, the library has no snorkel builder, and a refold failure
           would change the call.

A refold verdict, when l1-tag-refold reports, overrides the pick: if rank 1 fails and rank 2
passes, the site moves to rank 2.

    uv run --locked python l1_internal_decisions.py      # writes l1_internal_decisions.tsv
"""

from __future__ import annotations

import io
import json
import sys
from pathlib import Path

import joblib
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[3] / "src"))

import membrane_mask  # noqa: E402
import tag_site_model_v4 as V4  # noqa: E402
import tag_site_model_v5 as V5  # noqa: E402
from refold_queue import ORTHOLOG_PATH  # noqa: E402

OUT = HERE / "l1_internal_decisions.tsv"
#: --check: no network, verdicts from the committed cache only, outputs to CHECK_DIR
OFFLINE = False
CHECK_DIR = HERE / "data" / "check"
OUT_JSON = HERE / "l1_internal_decisions.json"
#: the adopted panel's final head (PR #283): membership and member ids come from here
MEMBERSHIP_COMMIT = (
    "74492351"  # main after #254; same 581 members as 6bf341f2, member ids renumbered
)
CDS_FINAL_KEY = (
    "sprint4/s4e9-orf-design-fg-panel/library-candidates/2026-11-12-reference-fixes/"
    "candidate_cds.json"
)
WINDOW = 10  # residues either side of a junction that must align identically to transfer it
UNCONFIRMED_DESIGN = "tag_sites_v5_orientation_unconfirmed.tsv"
ISOFORM_FIX_DESIGN = "tag_sites_v5_isoform_fix.tsv"
DECIDED_BY = 'ranker-pu-segment (delegated by Becca 2026-11-12, "use your best judgment")'
IMPLEMENTED_BY = "ranker-pu-segment session (Claude Code), PR #245"
SINGLE_PASS_OR_GPI = {"ATRN", "EPHA6", "NRG2", "RGMB"}
INVERSION_SUSPECT = {"SLC44A1", "TMC2"}


def design_unconfirmed(mask_df: pd.DataFrame, accs: set[str]) -> pd.DataFrame:
    """v5 design of the orientation-unconfirmed members on their masks as computed."""
    m = mask_df.copy()
    sha = membrane_mask.accession_to_sha()
    targets = pd.read_csv(HERE / "l1_adopted_targets.tsv", sep="\t")
    sub = targets[targets["uniprot_acc"].isin(accs)]
    m.loc[m["sequence_sha256"].isin({sha[a] for a in accs}), "usable"] = 1
    sub.to_csv(HERE / "data" / "unconfirmed_targets.tsv", sep="\t", index=False)
    real_load = membrane_mask.load
    membrane_mask.load = lambda: m  # these rows only; recorded in the decision reason
    try:
        spec = json.loads(V5.SPEC.read_text())
        V4.design(
            joblib.load(V5.MODEL_Z),
            joblib.load(V5.MODEL_RAW),
            spec,
            np.array(spec["published_protein_ceilings_oof"]),
            np.array(spec["published_region_ceilings_oof"]),
            "data/unconfirmed_targets.tsv",
            UNCONFIRMED_DESIGN,
            V5.EXTRA,
            use_mask=True,
            skip_refused=True,
            min_region_len=V5.MIN_REGION_LEN,
        )
    finally:
        membrane_mask.load = real_load
    return pd.read_csv(HERE / UNCONFIRMED_DESIGN, sep="\t")


LEFT: list = []

#: Build-driven switches: design_id -> (rank to use, reason). Each is a decision under the
#: 2026-11-12 delegation, raised by the arm builder (l1-integration, PR #299).
BUILD_OVERRIDES = {
    "l1-p1253": (
        2,
        "rank 1 (91) completes an Esp3I site at the cassette seam, and the builder never recodes "
        "a seam silently (l1-integration, PR #299). Rank 2 (167) builds cleanly, sits in the "
        "protein's best-scoring loop, and has the higher absolute score (0.668 vs 0.520)",
    ),
}


def apply_build_overrides(out, designs):
    for did, (rank, why) in BUILD_OVERRIDES.items():
        i = out.index[out["design_id"] == did]
        if not len(i):
            continue
        sha = out.loc[i[0], "sequence_sha256"]
        alt = designs[(designs["sequence_sha256"] == sha) & (designs["rank"] == rank)].iloc[0]
        old = int(out.loc[i[0], "insert_after_residue"])
        for col in ("rank_score", "absolute_score", "this_region_length"):
            out.loc[i, col] = alt[col]
        out.loc[i, "insert_after_residue"] = int(alt["insert_after_residue"])
        out.loc[i, "pick_rank"] = rank
        out.loc[i, "reason"] = why
        out.loc[i, "rejected_alternative"] = f"rank 1 at residue {old} (Esp3I site at the seam)"
    return out


def rekey_to_head(out, designs, orth, unconf):
    """Final-head membership and ids; transfer a site to sequences decided under no other rule."""
    internal, design_of, prot = head_tables()
    by_sha = out.set_index("sequence_sha256")
    picks = pd.concat(
        [
            designs[
                [
                    "sequence_sha256",
                    "gene",
                    "rank",
                    "insert_after_residue",
                    "rank_score",
                    "absolute_score",
                    "best_region_pct_of_published_regions",
                    "this_region_length",
                    "snorkel_candidate",
                ]
            ],
            orth[
                [
                    "sequence_sha256",
                    "gene",
                    "rank",
                    "insert_after_residue",
                    "rank_score",
                    "absolute_score",
                    "best_region_pct_of_published_regions",
                    "this_region_length",
                    "snorkel_candidate",
                ]
            ],
        ]
    )
    rows = []
    for m in internal.itertuples(index=False):
        s = m.sequence_sha256
        ids = {"member_id": m.candidate_id, "design_id": design_of.get(m.candidate_id)}
        if s in by_sha.index:
            r = by_sha.loc[s].to_dict()
            rows.append({**r, **ids, "sequence_sha256": s})
            continue
        dst = prot[s]
        done = False
        src_rows = out[
            (out["gene"] == m.gene_symbol) & out["decision"].str.startswith("internal_loop")
        ]
        for src in src_rows.itertuples(index=False):
            cand = picks[picks["sequence_sha256"] == src.sequence_sha256].sort_values("rank")
            for c in cand.itertuples(index=False):
                n2 = transfer(prot[src.sequence_sha256], dst, int(c.insert_after_residue))
                if n2 is None:
                    continue
                rows.append(
                    {
                        **ids,
                        "gene": m.gene_symbol,
                        "sequence_sha256": s,
                        "candidate_set": m.candidate_set,
                        "decision": "internal_loop_transferred",
                        "insert_after_residue": n2,
                        "pick_rank": int(c.rank),
                        "rank_score": c.rank_score,
                        "absolute_score": c.absolute_score,
                        "best_region_pct_of_published_regions": c.best_region_pct_of_published_regions,
                        "this_region_length": c.this_region_length,
                        "snorkel_candidate": bool(c.snorkel_candidate),
                        "confidence": "low",
                        "reason": f"new sequence at {MEMBERSHIP_COMMIT} (no features or mask row yet); "
                        f"site transferred from {src.member_id} (same gene) rank {int(c.rank)} "
                        f"residue {int(c.insert_after_residue)}: the junction and +/-{WINDOW} "
                        "residues align identically, so it is the same loop. Rescore with v5 once "
                        "per-residue features and a membrane-side mask row exist",
                        "rejected_alternative": "leave undecided until features exist; or snorkel",
                    }
                )
                done = True
                break
            if done:
                break
        if not done:
            c_side = v1_c_terminal_side(dst)
            base = {
                **ids,
                "gene": m.gene_symbol,
                "sequence_sha256": s,
                "candidate_set": m.candidate_set,
                "confidence": "low",
            }
            if c_side == "I":
                rows.append(
                    {
                        **base,
                        "decision": "snorkel",
                        "reason": f"new sequence at {MEMBERSHIP_COMMIT} with no features or mask "
                        "row, and no decided site of the same gene aligns identically over the "
                        "junction window; DeepTMHMM v1 puts the C-terminus in the cytoplasm, so "
                        "a snorkel (Brown et al. 2013, PMID 24023844) displays the tag outside. "
                        "Rescore with v5 once features and a mask row exist",
                        "rejected_alternative": "an untransferred guess at an internal site",
                    }
                )
            else:
                rows.append(
                    {
                        **base,
                        "decision": "hold_pending_gate",
                        "reason": f"new sequence at {MEMBERSHIP_COMMIT} with no features or mask "
                        "row, and no decided site of the same gene aligns identically; DeepTMHMM "
                        f"v1 puts the C-terminus on side {c_side!r}, so a snorkel would put the "
                        "tag inside the cell. Held until the membrane-side gate writes a mask row "
                        "and v5 can score it",
                        "rejected_alternative": "snorkel (wrong side under v1); an internal site "
                        "ranked without the mask (the ranker refuses unmasked targets)",
                    }
                )
    head = set(internal["sequence_sha256"])
    LEFT.extend(
        out.loc[~out["sequence_sha256"].isin(head), ["member_id", "gene", "decision"]].to_dict(
            orient="records"
        )
    )
    return pd.DataFrame(rows)


def v1_c_terminal_side(seq: str) -> str | None:
    """DeepTMHMM v1 label of the last residue of exactly this sequence (None if no v1 row)."""
    from deliverome_analysis.surfaceome import topology

    t = topology.for_sequence(seq)
    if t is None:
        return None
    s = t["per_residue_topology"] if isinstance(t, dict) else t.per_residue_topology
    return s[-1] if s else None


#: The gate's mask at MEMBERSHIP_COMMIT (PR #291, 8325d6d0). The v5 designs were made on r4;
#: r4 -> r5 changes no designed member's eligibility (only ATRN and NRG2, both routed to the
#: N-terminus), so r5 is used here only to check rows decided without the mask.
GATE_R5_KEY = (
    "sprint4/s4e9-orf-design-fg-panel/membrane-side-gate/2026-11-12-r8/"
    "l1_extracellular_mask_residues.tsv"
)


#: v5 designs of the 11 start-trimmed sequences on mask r5 (PR #298 at 5e0e747b, by the
#: l1-ortholog-features session). They replace the transfer / snorkel / hold fallbacks.
START_TRIMS_KEY = "sprint5/s5e13-tag-site-benchmarking/tag_sites_v5_start_trims.tsv"
#: b7a54e89 adds DISP2 and MFSD2A x2 on mask r8; earlier 6615b4a8, a2ef41ca, 5e0e747b
START_TRIMS_COMMIT = "b7a54e89"


def apply_late_designs(out: pd.DataFrame) -> pd.DataFrame:
    """Members designed after the main run take their v5 rank-1 pick, as every other does."""
    local = HERE / "data" / "tag_sites_v5_start_trims.tsv"
    if not local.is_file():
        from deliverome_analysis.core import s3

        s3.processed.download(START_TRIMS_KEY, str(local))
    late = pd.read_csv(local, sep="\t")
    r1 = late[late["rank"] == 1].set_index("sequence_sha256")
    for i in out.index[out["sequence_sha256"].isin(r1.index)]:
        r = r1.loc[out.loc[i, "sequence_sha256"]]
        was = out.loc[i, "decision"]
        if bool(r["snorkel_candidate"]):
            reason = (
                "v5 rank-1 eligible junction on the gate mask (r5-r7); kept on the loop although it is a "
                "snorkel candidate, as for every designed member. A failed refold moves it to "
                "rank 2 or to a snorkel"
            )
            alt, conf = "snorkel construct", "low"
        else:
            reason = "v5 rank-1 eligible junction under the gate mask (r5-r7)"
            alt, conf = "v5 rank-2 pick", "medium"
        if float(r.get("plddt", 100)) < 50:
            reason += f"; the pick sits at pLDDT {float(r['plddt']):.0f}, so the structure evidence is thin"
            conf = "low"
        reason += f" (designed late, PR #298 at {START_TRIMS_COMMIT}; replaces the {was} fallback)"
        for k, v in site_row(
            r.to_dict() | {"rank": r["rank"]}, "internal_loop", reason, alt, conf
        ).items():
            out.loc[i, k] = v
    return out


#: Members routed to N-after-SP that the N-terminal rung cannot build: SignalP 6 calls no
#: signal peptide on any native or Igk form (#254, nterm-reviewer-package; relayed by
#: l1-integration 2026-11-12).
NO_SP_ON_NTERM_REVIEW = {"ATRN", "NRG2"}


#: Single-pass members whose extracellular face UniProt fixes without the gate (which needs two
#: confident helices to fit a membrane): gene -> (accession, UniProt extracellular span, why).
SINGLE_PASS_ECTO = {
    "ATRN": (
        "O75882",
        (84, 1279),
        "UniProt O75882 (identical sequence): signal peptide 1-28 and propeptide 29-83 are "
        "cleaved, so the hydrophobic 66-78 that DeepTMHMM v1 counts as a second TM is not in "
        "the mature protein; extracellular 84-1279, TM 1280-1300, cytoplasmic 1301-1429, and "
        "the 15 residues after the TM carry 5 K/R",
    ),
}
SINGLE_PASS_DESIGN = HERE / "tag_sites_v5_single_pass.tsv"

#: Members routed INTO the internal rung after the membership commit, keyed by design_id.
#: Each is designed with v5 unchanged, without the mask (it has no row in r8), kept to the
#: UniProt ectodomain; then the contact veto and the refold rule apply as for every member.
ROUTED_IN = {
    "l1-p0567": {
        "gene": "ENTPD2",
        "uniprot_acc": "Q9Y5L3",
        "ectodomain": (29, 462),
        "why": "SignalP 6 calls no signal peptide; UniProt Q9Y5L3 (identical sequence) has "
        "TMs 8-28 and 463-483 with both termini cytoplasmic and the ectodomain 29-462 between "
        "them, and DeepTMHMM v2 places the first TM at 8-25. DeepTMHMM v1 reads 1-28 as a "
        "signal peptide, so the earlier N-after-SP route (Igk swap) put the tag on a "
        "cytoplasmic N-terminus",
        "decided_by": "coordinator, under Becca's 2026-11-12 delegation, on the N-terminal "
        "session's sp_disagree_review",
        "decided_on": "2026-11-14",
        "features": "entpd2_extracellular.tsv",
        "sequence_accession": "Q9Y5L3",
        "gate": "the membrane-side gate session (2026-11-14, on UniProt topology and the AFDB v6 "
        "model, TMDET frame) passes both sites: after 34 at +31.8 A, after 458 at +18.4 A "
        "(cutoff 14.0 A); every junction 29-462 passes, the lowest +14.1 A",
    },
    "l1-p2987": {
        "gene": "ENTPD2",
        "uniprot_acc": "l1-m1445@74492351",  # features keyed on the member (ortholog recipe)
        "sequence_accession": "A0A2K5VH27",
        "cohort": "cyno_ortholog",
        "ectodomain": (29, 462),
        "why": "cyno ortholog of ENTPD2, on its own evidence: no native SignalP 6 call (class "
        "refused_no_signalp_row in every nterm run; the 0.98 was the Igk leader of the swapped "
        "construct, confirmed by nterm-reviewer-package 2026-11-14); DeepTMHMM v1 SP 1-28 only; "
        "DeepTMHMM v2 on the designed cyno "
        "sequence calls TMs 9-26 and 466-484, both termini cytoplasmic and no signal peptide; the "
        "protein is 96.8% identical to human ENTPD2, aligns gap-free, residues 1-30 identical, so "
        "human TM 8-28 and ectodomain 29-462 map to the same positions. DeepTMHMM v1 reads 1-28 "
        "as a signal peptide, as for human",
        "decided_by": "coordinator, under Becca's 2026-11-12 delegation, on the N-terminal "
        "session's sp_disagree_review (cyno applied on its own sequence evidence)",
        "decided_on": "2026-11-14",
        "features": "l1_routed_in_member_extracellular.tsv",
        "gate": "the membrane-side gate session (2026-11-14, on the DeepTMHMM v2 topology and the "
        "member's ColabFold model, helix-axis frame; AFDB v6 A0A2K5VH27 as cross-check) passes both "
        "sites: after 33 at +31.0 A, after 458 at +21.2 A (cutoff 14.0 A); every junction 27-464 "
        "passes",
    },
}
ROUTED_IN_DESIGN = HERE / "tag_sites_v5_routed_in.tsv"


def design_routed_in() -> pd.DataFrame:
    """v5 (unchanged) for ROUTED_IN members, without the mask, kept to the UniProt ectodomain."""
    _internal, _ids, prot = head_tables()
    tp = full_head_paths()
    rows = []
    for did, r in ROUTED_IN.items():
        m = tp[tp["design_id"] == did].iloc[0]
        rows.append(
            {
                "gene": r["gene"],
                "uniprot_acc": r["uniprot_acc"],
                "member_id": m["candidate_id"],
                "candidate_set": m["candidate_set"],
                "topology_cohort": r.get("cohort", "human_canonical"),
                "sequence_sha256": m["sequence_sha256"],
            }
        )
    pd.DataFrame(rows).to_csv(HERE / "data" / "routed_in_targets.tsv", sep="\t", index=False)
    spec = json.loads(V5.SPEC.read_text())
    d = V4.design(
        joblib.load(V5.MODEL_Z),
        joblib.load(V5.MODEL_RAW),
        spec,
        np.array(spec["published_protein_ceilings_oof"]),
        np.array(spec["published_region_ceilings_oof"]),
        "data/routed_in_targets.tsv",
        str(ROUTED_IN_DESIGN.relative_to(HERE)),
        V5.EXTRA,
        use_mask=False,  # no r8 row exists for these sequences; the face comes from UniProt
        skip_refused=True,
        min_region_len=V5.MIN_REGION_LEN,
        override_tables=tuple(r["features"] for r in ROUTED_IN.values()),
    )
    span = {r["uniprot_acc"]: r["ectodomain"] for r in ROUTED_IN.values()}
    inside = d.apply(
        lambda r: (
            span[r["uniprot_acc"]][0] <= r["insert_after_residue"] < span[r["uniprot_acc"]][1]
        ),
        axis=1,
    )
    d = d[inside].copy()
    d["sequence_sha256"] = d["uniprot_acc"].map(
        {t["uniprot_acc"]: t["sequence_sha256"] for t in rows}
    )
    for t in rows:
        assert prot[t["sequence_sha256"]], t
    d.to_csv(ROUTED_IN_DESIGN, sep="\t", index=False)
    return d


def full_head_paths() -> pd.DataFrame:
    """Every member at MEMBERSHIP_COMMIT (any placement), with main's member and design ids."""
    repo = HERE.parents[2]
    base = "library/l1/orf-design/"

    def show(path):
        from deliverome_analysis.core.provenance import read_pinned_text

        return pd.read_csv(
            io.StringIO(read_pinned_text(MEMBERSHIP_COMMIT, f"{base}{path}", repo=repo)),
            sep="\t",
            low_memory=False,
        )

    tp = show("l1_tag_paths.tsv")
    ren = show("l1_member_renumbering.tsv")
    return tp.merge(ren[["candidate_id", "design_id"]], on="candidate_id", how="left")


def apply_routed_in(out: pd.DataFrame) -> pd.DataFrame:
    """Add each ROUTED_IN member as an internal-loop row on its v5 rank-1 ectodomain site."""
    d = design_routed_in()
    tp = full_head_paths()
    add = []
    for did, r in ROUTED_IN.items():
        if (out["design_id"] == did).any():
            continue
        m = tp[tp["design_id"] == did].iloc[0]
        pick = d[d["uniprot_acc"] == r["uniprot_acc"]].sort_values("rank").iloc[0].to_dict()
        add.append(
            {
                "member_id": m["candidate_id"],
                "design_id": did,
                "gene": r["gene"],
                "sequence_sha256": m["sequence_sha256"],
                "candidate_set": m["candidate_set"],
                **site_row(
                    pick,
                    "internal_loop",
                    f"routed in from N after SP: {r['why']}. v5 rank-{int(pick['rank'])} site in "
                    f"the ectodomain {r['ectodomain'][0]}-{r['ectodomain'][1]}, scored without the "
                    f"mask (no r8 row for this sequence); {r['gate']}",
                    "N after SP with an Igk swap (the tag would sit on a cytoplasmic N-terminus)",
                    "low",
                ),
            }
        )
    return pd.concat([out, pd.DataFrame(add)], ignore_index=True) if add else out


def design_single_pass() -> pd.DataFrame:
    """v5 without the mask, for SINGLE_PASS_ECTO members only, kept to the UniProt ectodomain."""
    targets = pd.read_csv(HERE / "l1_adopted_targets.tsv", sep="\t")
    accs = {a for a, _span, _why in SINGLE_PASS_ECTO.values()}
    targets[targets["uniprot_acc"].isin(accs)].to_csv(
        HERE / "data" / "single_pass_targets.tsv", sep="\t", index=False
    )
    spec = json.loads(V5.SPEC.read_text())
    d = V4.design(
        joblib.load(V5.MODEL_Z),
        joblib.load(V5.MODEL_RAW),
        spec,
        np.array(spec["published_protein_ceilings_oof"]),
        np.array(spec["published_region_ceilings_oof"]),
        "data/single_pass_targets.tsv",
        str(SINGLE_PASS_DESIGN.relative_to(HERE)),
        V5.EXTRA,
        use_mask=False,  # no mask can exist for one helix; the face comes from UniProt
        skip_refused=True,
        min_region_len=V5.MIN_REGION_LEN,
    )
    span = {g: sp for g, (_a, sp, _w) in SINGLE_PASS_ECTO.items()}
    inside = d.apply(
        lambda r: span[r["gene"]][0] <= r["insert_after_residue"] < span[r["gene"]][1], axis=1
    )
    assert inside.all(), d[~inside]
    return d


def nterm_fallback(out: pd.DataFrame) -> pd.DataFrame:
    """Snorkel the N-after-SP members that have no signal peptide, if the C-terminus is inside."""
    import re

    from deliverome_analysis.surfaceome import topology

    _internal, _ids, prot = head_tables()
    for i in out.index[
        (out["decision"] == "route_to_nterm_after_sp") & out["gene"].isin(NO_SP_ON_NTERM_REVIEW)
    ]:
        gene = out.loc[i, "gene"]
        if gene in SINGLE_PASS_ECTO:
            sp = design_single_pass()
            r = sp[(sp["gene"] == gene) & (sp["rank"] == 1)].iloc[0].to_dict()
            row = site_row(
                r,
                "internal_loop",
                SINGLE_PASS_ECTO[gene][2]
                + ". v5 rank-1 in that ectodomain, scored without the mask because a single helix "
                "cannot be placed by the gate; the site is "
                f"{int(r['dist_to_tm'])} residues from the TM, far from the loop-base and short-loop "
                f"cases the mask exists to catch. Pick pLDDT {float(r['plddt']):.0f}. Raised by "
                "the manifest writer (#283) via l1-integration",
                "snorkel (adds a second TM; the ectodomain face is not in doubt)",
                "low",
            )
            for k, v in row.items():
                out.loc[i, k] = v
            continue
        seq = prot[out.loc[i, "sequence_sha256"]]
        t = topology.for_sequence(seq)
        topo = t["per_residue_topology"]
        last_tm_end = max(m.end() for m in re.finditer("M+", topo))
        inside = sum(c in "KR" for c in seq[last_tm_end : last_tm_end + 15])
        if topo[-1] != "I" or inside < 3:
            out.loc[i, "decision"] = "no_arm"
            out.loc[i, "reason"] = (
                "no signal peptide (SignalP 6, #254) and the C-terminal side is not supported as "
                "cytoplasmic, so neither an N-terminal nor a snorkel arm displays the tag outside"
            )
            continue
        out.loc[i, "decision"] = "snorkel"
        out.loc[i, "confidence"] = "medium"
        out.loc[i, "reason"] = (
            "SignalP 6 calls no signal peptide on any native or Igk form (#254), so there is no "
            "N-after-SP arm; the membrane-side gate finds fewer than two confident helices, so "
            "the ranker refuses an internal site. The C-terminus is cytoplasmic on two "
            f"independent signals: DeepTMHMM v1 ({t['tool_version']}) ends it 'I', and the 15 "
            f"residues after the last TM carry {inside} K/R (positive-inside). A snorkel (Brown "
            "et al. 2013, PMID 24023844) needs only that"
            + (
                "; the ectodomain is shed (shed_risk high, strong evidence; UniProt O14511 chain "
                "112-404 is the released Neuregulin-2), and v5's best ectodomain site (153, "
                "pLDDT 33) lies in that chain, so its tag would leave with the shed ectodomain; "
                "the snorkel keeps the tag membrane-anchored"
                if out.loc[i, "gene"] == "NRG2"
                else ""
            )
        )
        out.loc[i, "rejected_alternative"] = (
            "an internal site in the ectodomain ranked without the mask (refused by rule); no arm"
        )
    return out


#: l1-tag-refold's verdicts on cassette-inserted folds (PR #254 family), joined on the tagged
#: protein's sha256. A snapshot: the table refreshes as folds land, so rerun to pick them up.
VERDICTS_KEY = (
    "sprint4/s4e9-orf-design-fg-panel/predicted-structures/tagged/2026-11-12/verdicts_loop.tsv"
)
SNORKEL_VERDICTS_KEY = (
    "sprint4/s4e9-orf-design-fg-panel/predicted-structures/tagged/2026-11-12/verdicts_cterm.tsv"
)
REFOLD_PENDING = HERE / "l1_refold_alternatives_requested.tsv"
VERDICTS_SEEN = HERE / "refold_verdicts_seen.tsv"


def _tagged_sha(prot: dict, sha: str, n: int) -> str:
    import hashlib

    from build_internal_tagged import CASSETTE_AA

    p = prot[sha]
    return hashlib.sha256((p[:n] + CASSETTE_AA + p[n:]).encode()).hexdigest()


def _gate_confirmed(sha: str) -> bool:
    m = pd.read_csv(HERE / "data" / "mask_r8" / "l1_extracellular_mask_residues.tsv", sep="\t")
    row = m[m["sequence_sha256"] == sha]
    return bool(len(row)) and int(row["usable"].iloc[0]) == 1


def _cyto_c_terminus(seq: str, gate_ok: bool = False) -> tuple[bool, str]:
    """v1 says the C-terminus is inside, and either the gate confirmed the member's face from
    its structure or the positive-inside flank (>= 3 K/R in 15 residues) agrees."""
    import re

    from deliverome_analysis.surfaceome import topology

    t = topology.for_sequence(seq)
    if t is None:
        return False, "no exact-sequence DeepTMHMM v1 row"
    topo = t["per_residue_topology"]
    end = max(m.end() for m in re.finditer("M+", topo))
    kr = sum(c in "KR" for c in seq[end : end + 15])
    ev = f"v1 C-terminus {topo[-1]!r}, {kr} K/R after the last TM" + (
        ", face confirmed by the membrane-side gate (r8 usable)" if gate_ok else ""
    )
    return topo[-1] == "I" and (gate_ok or kr >= 3), ev


#: The internal-site policy, chosen by compare_internal_policies.py (2026-11-12): v5 plus a
#: pick-time veto, so no site sits within contact_veto.VETO_HALF of a ligand contact.
POLICY = (
    "v5 (RF + ligand-inclusive hard negatives, w=5) + ligand-contact veto, residues n-2..n+3 "
    "for an insertion after n"
)
POLICY_DATE = "2026-11-12"


def _history(out: pd.DataFrame, i, old, new, why: str, date: str = POLICY_DATE) -> None:
    # a member routed in later is decided under its own date, not the 2026-11-12 policy date
    routed = ROUTED_IN.get(out.loc[i, "design_id"]) if "design_id" in out else None
    if routed and date == POLICY_DATE:
        date = routed["decided_on"]
    entry = f"{date}: {old} -> {new}: {why}"
    prev = out.loc[i, "site_history"] if "site_history" in out else None
    out.loc[i, "site_history"] = f"{prev} | {entry}" if isinstance(prev, str) and prev else entry


def human_members_of(out: pd.DataFrame) -> dict[str, list[str]]:
    """design_id of each non-human member in `out` -> the design_ids of its gene's human members."""
    tp = full_head_paths()
    human = tp[tp["species"] == "human"].groupby("gene_symbol")["design_id"].apply(sorted)
    nh = tp[(tp["species"] != "human") & tp["design_id"].isin(set(out["design_id"]))]
    return {d: human.get(g, []) for d, g in zip(nh["design_id"], nh["gene_symbol"], strict=True)}


def apply_contact_veto(
    out: pd.DataFrame, alternatives: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Move any decided site within VETO_HALF of a ligand contact to the next clean site.

    Also returns the alternatives with every vetoed site removed, so a later refold move
    can only land on a clean site."""
    import contact_veto as C

    _internal, _ids, prot = head_tables()
    # a member not in the 2026-11-07 contact tables gets its human member's (an ortholog) or any
    # same-gene member's D1 and transferred contacts by alignment; the label says what was
    # actually looked up per member (contact_veto.contact_lookup)
    contacts, lookup = C.contact_lookup(
        out,
        prot,
        download=not OFFLINE,
        human_of=human_members_of(out),
        routed_ids=set(ROUTED_IN),
        gene_of=dict(zip(out["design_id"], out["gene"], strict=True)),
    )
    alts = alternatives.dropna(subset=["insert_after_residue"])
    out["contact_veto"] = ""
    label = {r["design_id"]: C.sources_available(r) for _, r in lookup.iterrows()}
    out["contact_sources_available"] = out["design_id"].map(label)
    for i in out.index[out["insert_after_residue"].notna()]:
        sh, n = out.loc[i, "sequence_sha256"], int(out.loc[i, "insert_after_residue"])
        src = C.vetoed(contacts, sh, n)
        if not src:
            continue
        mine = alts[(alts["sequence_sha256"] == sh) & (alts["insert_after_residue"] != n)]
        mine = mine.drop_duplicates("insert_after_residue").sort_values("rank")
        clean = [
            r
            for _, r in mine.iterrows()
            if not C.vetoed(contacts, sh, int(r["insert_after_residue"]))
        ]
        why = f"ligand contact within {C.VETO_HALF} of the junction ({', '.join(src)})"
        out.loc[i, "contact_veto"] = f"after {n}: {', '.join(src)}"
        if clean:
            a = clean[0]
            new = int(a["insert_after_residue"])
            row = site_row(
                a,
                "internal_loop",
                f"contact veto moved it: the site after {n} has a {why}; v5 rank "
                f"{int(a['rank'])} after {new} has none",
                f"the vetoed site after {n}",
                "medium",
            )
            for k, val in row.items():
                out.loc[i, k] = val
            _history(out, i, n, new, why)
        else:
            ok, ev = _cyto_c_terminus(prot[sh], _gate_confirmed(sh))
            if ok:
                for c in ("insert_after_residue", "pick_rank", "rank_score", "absolute_score"):
                    out.loc[i, c] = np.nan
                out.loc[i, "decision"] = "snorkel"
                out.loc[i, "confidence"] = "medium"
                out.loc[i, "reason"] = (
                    f"every designed internal site has a {why.split(' (')[0]}, so a snorkel "
                    f"(Brown et al. 2013, PMID 24023844) avoids the binding site; {ev}"
                )
                out.loc[i, "rejected_alternative"] = f"the vetoed internal site after {n}"
                _history(out, i, n, "snorkel", why)
            else:
                out.loc[i, "reason"] = (
                    f"{out.loc[i, 'reason']}. Kept despite a {why}: no clean designed site and a "
                    f"snorkel is not supported ({ev})"
                )
    keep = [
        not C.vetoed(contacts, sh, int(n))
        for sh, n in zip(alts["sequence_sha256"], alts["insert_after_residue"], strict=True)
    ]
    return out, alts[keep]


def apply_refold(out: pd.DataFrame, alternatives: pd.DataFrame) -> pd.DataFrame:
    """Fold verdicts in: a failed pick moves to a passing alternative, waits on an unfolded one,
    or, with none left, becomes a snorkel when the C-terminus is cytoplasmic."""
    from deliverome_analysis.core import s3

    local = HERE / "data" / "verdicts_loop.tsv"
    if OFFLINE:
        latest = pd.DataFrame(columns=["tagged_sha256", "verdict"])
    else:
        s3.processed.download(VERDICTS_KEY, str(local))  # always the latest snapshot
        latest = pd.read_csv(local, sep="\t")
    latest = latest[latest["verdict"].isin(["pass", "fail"])]
    # The live table is rewritten whole on every publish and has dropped rows before (52
    # alternative verdicts vanished between etags 8501b370 and d6662204). Keep every pass/fail
    # verdict ever seen, newest first, in a committed cache, so a partial republish cannot
    # silently undo a decision.
    keep = ["tagged_sha256", "source_id", "design_id", "verdict", "verdict_reasons"]
    if OFFLINE:
        merged = pd.read_csv(VERDICTS_SEEN, sep="\t")
    elif VERDICTS_SEEN.is_file():
        seen = pd.read_csv(VERDICTS_SEEN, sep="\t")
        merged = pd.concat([latest[keep], seen[keep]]).drop_duplicates("tagged_sha256")
    else:
        merged = latest[keep].drop_duplicates("tagged_sha256")
    if not OFFLINE:
        merged.sort_values("tagged_sha256").to_csv(VERDICTS_SEEN, sep="\t", index=False)
    v = merged.set_index("tagged_sha256")
    _internal, _ids, prot = head_tables()
    alts = alternatives.dropna(subset=["insert_after_residue"]).copy()
    alts = alts[alts["sequence_sha256"].isin(prot)]
    alts["tsha"] = [
        _tagged_sha(prot, sh, int(n))
        for sh, n in zip(alts["sequence_sha256"], alts["insert_after_residue"], strict=True)
    ]
    alts["verdict"] = alts["tsha"].map(v["verdict"])
    pending: list[dict] = []
    out["refold_verdict"] = "not_sent"
    for i in out.index[out["insert_after_residue"].notna()]:
        sh, n = out.loc[i, "sequence_sha256"], int(out.loc[i, "insert_after_residue"])
        verdict = v["verdict"].get(_tagged_sha(prot, sh, n), "not_in_table")
        out.loc[i, "refold_verdict"] = verdict
        if verdict == "not_in_table":
            pending.append(
                {
                    "design_id": out.loc[i, "design_id"],
                    "gene": out.loc[i, "gene"],
                    "sequence_sha256": sh,
                    "rank": int(out.loc[i, "pick_rank"]),
                    "insert_after_residue": n,
                    "tagged_protein_sha256": _tagged_sha(prot, sh, n),
                    "why": "decided site not yet folded",
                }
            )
        if verdict != "fail":
            continue
        why = v.loc[_tagged_sha(prot, sh, n), "verdict_reasons"]
        mine = alts[(alts["sequence_sha256"] == sh) & (alts["insert_after_residue"] != n)]
        mine = mine.drop_duplicates("insert_after_residue").sort_values("rank")
        passing = mine[mine["verdict"] == "pass"]
        untested = mine[~mine["verdict"].isin(["pass", "fail"])]
        if len(passing):
            a = passing.iloc[0]
            row = site_row(
                a,
                "internal_loop",
                f"refold moved it: the site after {n} fails ({why}); v5 rank {int(a['rank'])} "
                f"after {int(a['insert_after_residue'])} passes",
                f"the failed site after {n}",
                "medium",
            )
            for k, val in row.items():
                out.loc[i, k] = val
            out.loc[i, "refold_verdict"] = "pass"
            _history(out, i, n, int(a["insert_after_residue"]), f"refold fails ({why})")
        elif len(untested):
            out.loc[i, "refold_verdict"] = "fail_alternative_pending"
            out.loc[i, "reason"] = (
                f"{out.loc[i, 'reason']}. Refold fails this site ({why}); the alternative(s) "
                f"after {', '.join(str(int(x)) for x in untested['insert_after_residue'])} are "
                "queued for refolding, and the pick moves when one passes"
            )
            for r in untested.itertuples(index=False):
                pending.append(
                    {
                        "design_id": out.loc[i, "design_id"],
                        "gene": out.loc[i, "gene"],
                        "sequence_sha256": sh,
                        "rank": int(r.rank),
                        "insert_after_residue": int(r.insert_after_residue),
                        "tagged_protein_sha256": r.tsha,
                        "why": f"the decided site after {n} fails refold",
                    }
                )
        else:
            ok, ev = _cyto_c_terminus(prot[sh], _gate_confirmed(sh))
            tried = ", ".join(str(int(x)) for x in mine["insert_after_residue"]) or "none designed"
            if ok:
                for c in ("insert_after_residue", "pick_rank", "rank_score", "absolute_score"):
                    out.loc[i, c] = np.nan
                out.loc[i, "decision"] = "snorkel"
                out.loc[i, "confidence"] = "medium"
                out.loc[i, "reason"] = (
                    f"every designed internal site fails refold (after {n}: {why}; "
                    f"alternatives: {tried}), so a snorkel (Brown et al. 2013, PMID 24023844) "
                    f"displays the tag without opening a loop; {ev}"
                )
                out.loc[i, "rejected_alternative"] = f"the failed internal site after {n}"
                _history(out, i, n, "snorkel", f"refold fails ({why}), no passing alternative")
            else:
                out.loc[i, "refold_verdict"] = "fail_no_alternative"
                out.loc[i, "reason"] = (
                    f"{out.loc[i, 'reason']}. Refold fails this site ({why}) and every "
                    f"alternative ({tried}); kept because a snorkel is not supported ({ev})"
                )
    # Cumulative: l1-tag-refold builds its queue from this file, and an emptied copy once
    # dropped 52 served requests from its publish. Never shrink it; mark what is still open.
    new = pd.DataFrame(pending)
    old = (
        pd.read_csv(REFOLD_PENDING, sep="\t")
        if REFOLD_PENDING.is_file() and REFOLD_PENDING.stat().st_size > 1
        else pd.DataFrame()
    )
    allreq = pd.concat([old, new]).drop_duplicates("tagged_protein_sha256", keep="last")
    if len(allreq):
        allreq["still_pending"] = allreq["tagged_protein_sha256"].isin(
            set(new["tagged_protein_sha256"]) if len(new) else set()
        )
        allreq.to_csv(REFOLD_PENDING, sep="\t", index=False)
    return out


#: Off until the pre-registered spacer test (snorkel_spacer_test.json) reads out: the
#: coordinator's 2026-11-13 decision keeps the 22 failing snorkels as built until then.
HOLD_FAILING_SNORKELS = True  # spacer test failed 2026-11-13: 0/5 rescued, 0/2 controls


def apply_snorkel_verdicts(out: pd.DataFrame) -> pd.DataFrame:
    """A snorkel whose own refold fails has no arm left: hold it, with the evidence.

    Rule (coordinator, under Becca's 2026-11-12 delegation, 2026-11-13): a failing arm moves
    to the next v5 site that passes the contact veto; with none, the member is held. Every
    snorkel here became one only after each designed internal site failed refold or sat on a
    ligand contact, so a failing snorkel has no site left to move to.
    """
    if not HOLD_FAILING_SNORKELS:
        return out
    import hashlib

    from snorkel_design import MODULE

    local = HERE / "data" / "verdicts_cterm.tsv"
    if not OFFLINE:
        from deliverome_analysis.core import s3

        s3.processed.download(SNORKEL_VERDICTS_KEY, str(local))
        latest = pd.read_csv(local, sep="\t")
        latest = latest[latest["verdict"].isin(["pass", "fail"])]
        keep = ["tagged_sha256", "source_id", "design_id", "verdict", "verdict_reasons"]
        seen = pd.read_csv(VERDICTS_SEEN, sep="\t")
        merged = pd.concat([latest[keep], seen[keep]]).drop_duplicates("tagged_sha256")
        merged.sort_values("tagged_sha256").to_csv(VERDICTS_SEEN, sep="\t", index=False)
    v = pd.read_csv(VERDICTS_SEEN, sep="\t").set_index("tagged_sha256")
    _internal, _ids, prot = head_tables()
    out["snorkel_verdict"] = ""
    for i in out.index[out["decision"] == "snorkel"]:
        sha = hashlib.sha256((prot[out.loc[i, "sequence_sha256"]] + MODULE).encode()).hexdigest()
        if sha not in v.index:
            out.loc[i, "snorkel_verdict"] = "not_folded"
            continue
        verdict, why = v.loc[sha, "verdict"], v.loc[sha, "verdict_reasons"]
        out.loc[i, "snorkel_verdict"] = verdict
        if verdict != "fail":
            continue
        out.loc[i, "decision"] = "hold_no_passing_arm"
        out.loc[i, "confidence"] = "high"
        out.loc[i, "rejected_alternative"] = (
            f"snorkel arm (refold fails: {why}); {out.loc[i, 'rejected_alternative']}"
        )
        out.loc[i, "reason"] = (
            f"the snorkel arm fails refold ({why}), and no internal site is left: "
            f"{out.loc[i, 'reason']}. Held with no arm"
        )
        _history(
            out,
            i,
            "snorkel",
            "hold",
            f"snorkel refold fails ({why}); a 30-residue L1 spacer rescued 0 of 5 "
            "(snorkel_spacer_test.json)",
            date="2026-11-13",
        )
    return out


#: Members whose experimental topology contradicts v1 at the decided site (2026-11-13).
EXPERIMENTAL_TOPOLOGY = {
    "l1-p1995": (
        "route_to_cterm_experimental_topology",
        "UniProt Q14CX5 (identical sequence) has experimental topology (ECO:0000269; Anzai & "
        "Matsumura 2019, PMID 31615651) with both termini extracellular, the reverse of "
        "DeepTMHMM v1. The v5 site after 446 lies in its cytoplasmic segment 430-466, so the "
        "tag would sit inside the cell. Its C-terminus (486-517) is extracellular, so a plain "
        "C-terminal tag displays outside; the rung move is the manifest writer's",
        "the internal site after 446 (cytoplasmic under the experimental topology); an internal "
        "site in the experimental extracellular loop 224-264 (v1 calls it inside, so v5 never "
        "scored it)",
    ),
}


def apply_experimental_topology(out: pd.DataFrame) -> pd.DataFrame:
    for did, (decision, reason, alt) in EXPERIMENTAL_TOPOLOGY.items():
        for i in out.index[out["design_id"] == did]:
            old = out.loc[i, "insert_after_residue"]
            for c in ("insert_after_residue", "pick_rank", "rank_score", "absolute_score"):
                out.loc[i, c] = np.nan
            out.loc[i, "decision"] = decision
            out.loc[i, "confidence"] = "high"
            out.loc[i, "reason"] = reason
            out.loc[i, "rejected_alternative"] = alt
            _history(out, i, int(old) if old == old else "none", "C-terminal rung", "experimental "
                     "topology puts the site inside the cell", date="2026-11-13")  # fmt: skip
    return out


def gate_check(out: pd.DataFrame) -> pd.DataFrame:
    """Hold fallback rows (transferred sites, fallback snorkels) the r5 gate contradicts."""
    local = HERE / "data" / "mask_r8" / "l1_extracellular_mask_residues.tsv"
    if not local.is_file():
        from deliverome_analysis.core import s3

        local.parent.mkdir(parents=True, exist_ok=True)
        s3.processed.download(GATE_R5_KEY, str(local))
    m = pd.read_csv(local, sep="\t").set_index("sequence_sha256")
    fallback = out["decision"].isin(["internal_loop_transferred", "hold_pending_gate"]) | (
        (out["decision"] == "snorkel") & out["reason"].str.startswith("new sequence")
    )
    # r8 (PR #291 at 419409f1) confirms the face of members the gate had left unconfirmed:
    # their rows become ordinary internal loops when the decided junction passes.
    for i in out.index[out["decision"] == "internal_loop_orientation_unconfirmed"]:
        sh = out.loc[i, "sequence_sha256"]
        if sh in m.index and m.loc[sh, "usable"] == 1:
            n = int(out.loc[i, "insert_after_residue"])
            if m.loc[sh, "junction_mask"][n - 1] == "E":
                out.loc[i, "decision"] = "internal_loop"
                out.loc[i, "confidence"] = "medium"
                out.loc[i, "reason"] = (
                    f"{out.loc[i, 'reason']}. Mask r8 confirms the face (positive-inside now "
                    "read on confident residues, and v2 agrees with v1's termini) and the "
                    f"junction after {n} passes"
                )
    for i in out.index[fallback]:
        r = out.loc[i]
        if r["sequence_sha256"] not in m.index:
            continue
        g = m.loc[r["sequence_sha256"]]
        if g["usable"] != 1 and g["status"] == "placed":
            out.loc[i, "decision"] = "hold_pending_topology"
            out.loc[i, "reason"] = (
                "mask r5 (PR #291) places the protein but cannot confirm v1's face: the "
                "positive-inside check contradicts it, a possible whole-protein inversion. No "
                "internal site and no snorkel until the topology is checked against a structure"
            )
            out.loc[i, "rejected_alternative"] = f"the fallback ({r['decision']}) it replaces"
        elif r["decision"] == "internal_loop_transferred" and g["usable"] == 1:
            n = int(r["insert_after_residue"])
            if g["junction_mask"][n - 1] != "E":
                out.loc[i, "decision"] = "hold_pending_design"
                out.loc[i, "reason"] = (
                    f"the transferred junction after {n} fails mask r5 "
                    f"(junction_mask {g['junction_mask'][n - 1]!r}); held until v5 designs this "
                    "sequence on r5"
                )
                out.loc[i, "rejected_alternative"] = f"the transferred site after {n}"
        if out.loc[i, "decision"].startswith("hold"):
            for c in ("insert_after_residue", "pick_rank", "rank_score", "absolute_score"):
                out.loc[i, c] = np.nan
    return out


#: PR #298's 27-ortholog designs, pinned (was the moving branch head)
ORTHOLOG_COMMIT = "dc0ab0fc"


def ortholog_designs() -> pd.DataFrame:
    from deliverome_analysis.core.provenance import read_pinned_text

    repo = HERE.parents[2]
    txt = read_pinned_text(ORTHOLOG_COMMIT, ORTHOLOG_PATH, repo=repo)
    return pd.read_csv(io.StringIO(txt), sep="\t")


def site_row(r, kind, reason, alternative, confidence):
    return {
        "decision": kind,
        "insert_after_residue": int(r["insert_after_residue"]),
        "pick_rank": int(r["rank"]),
        "rank_score": r["rank_score"],
        "absolute_score": r["absolute_score"],
        "best_region_pct_of_published_regions": r["best_region_pct_of_published_regions"],
        "this_region_length": r["this_region_length"],
        "snorkel_candidate": bool(r["snorkel_candidate"]),
        "confidence": confidence,
        "reason": reason,
        "rejected_alternative": alternative,
    }


def head_tables():
    """Final-head internal members, their ids, and every designed protein needed to align."""
    import hashlib

    from Bio.Seq import Seq

    from deliverome_analysis.core.provenance import read_pinned_text

    repo = HERE.parents[2]
    base = "library/l1/orf-design/"
    tp = pd.read_csv(
        io.StringIO(read_pinned_text(MEMBERSHIP_COMMIT, f"{base}l1_tag_paths.tsv", repo=repo)),
        sep="\t",
        low_memory=False,
    )
    ren = pd.read_csv(
        io.StringIO(
            read_pinned_text(MEMBERSHIP_COMMIT, f"{base}l1_member_renumbering.tsv", repo=repo)
        ),
        sep="\t",
    )
    local = HERE / "data" / "candidate_cds_2026-11-12.json"
    if not local.is_file():
        from deliverome_analysis.core import s3

        s3.processed.download(CDS_FINAL_KEY, str(local))
    prot = {}
    for path in (
        HERE / "data" / "candidate_cds_2026-11-12-ortholog-start-trims.json",
        HERE / "data" / "candidate_cds_2026-11-12-reference-fixes-2.json",
        local,
        HERE / "data" / "candidate_cds_current.json",
        HERE / "data" / "candidate_cds.json",
    ):
        if not path.is_file():  # a missing table would silently drop proteins
            raise SystemExit(f"missing input {path.relative_to(HERE)}; run --restore")
        for v in json.loads(path.read_text()).values():
            if v.get("cds"):
                q = str(Seq(v["cds"]).translate()).rstrip("*")
                prot.setdefault(hashlib.sha256(q.encode()).hexdigest(), q)
    # The CDS record of record (library_design.CDS_KEY): proteins of members whose reference
    # moved since the saved inputs (the 2026-11-22 GTF MANE switches, e.g. P2RX5), so their own
    # 2026-11-07 contacts can still be carried onto the protein now ordered.
    from deliverome_analysis.cloning import library_design as ld

    for v in ld.load_designed_cds().values():
        if v.get("cds"):
            q = str(Seq(v["cds"]).translate()).rstrip("*")
            prot.setdefault(hashlib.sha256(q.encode()).hexdigest(), q)
    internal = tp[tp["placement"] == "internal"]
    design_of = dict(zip(ren["candidate_id"], ren["design_id"], strict=False))
    return internal, design_of, prot


def transfer(src: str, dst: str, n: int) -> int | None:
    """Junction n of `src` in `dst` numbering, if residues n-WINDOW..n+1+WINDOW align identically."""
    from Bio import Align
    from Bio.Align import substitution_matrices

    al = Align.PairwiseAligner(mode="global")
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    al.open_gap_score, al.extend_gap_score = -10, -0.5
    a = al.align(src, dst)[0]
    pos = {}
    for (s0, s1), (d0, _d1) in zip(*a.aligned, strict=False):
        for k in range(s1 - s0):
            pos[s0 + k + 1] = d0 + k + 1
    lo, hi = max(1, n - WINDOW), min(len(src), n + 1 + WINDOW)
    for i in range(lo, hi + 1):
        j = pos.get(i)
        if j is None or src[i - 1] != dst[j - 1] or (i > lo and pos.get(i - 1) != j - 1):
            return None
    return pos[n]


#: every input the decisions read, hashed into the JSON so --check can prove they are unchanged
INPUT_FILES = (
    "tag_sites_v5_isoform_fix.tsv",
    "data/allx_isoform_fix15_2026-11-14.tsv",
    "member_interface_contacts.tsv",
    "member_transferred_ligand_contacts.tsv",
    "l1_internal_members.tsv",
    "data/lit_ligand_contacts_2026-11-14.tsv",
    "data/lit_ligand_coverage_2026-11-14.tsv",
    "l1_adopted_targets.tsv",
    "tag_sites_v5_adopted.tsv",
    "tag_site_model_v5_rank.joblib",
    "tag_site_model_v5_absolute.joblib",
    "tag_site_model_v5.json",
    "refold_verdicts_seen.tsv",
    "data/tag_sites_v5_start_trims.tsv",
    "data/mask_r4/l1_extracellular_mask_residues.tsv",
    "data/mask_r8/l1_extracellular_mask_residues.tsv",
    "data/candidate_cds_2026-11-12-ortholog-start-trims.json",
    "data/candidate_cds_2026-11-12-reference-fixes-2.json",
    "data/candidate_cds_2026-11-12.json",
    "data/candidate_cds_current.json",
    "data/candidate_cds.json",
    "data/l1_tag_paths_adopted_9a1c2c91.tsv",
    "l1_member_extracellular.tsv",
    "all_extracellular_targets.tsv",
    "slc67a1_extracellular.tsv",
    "entpd2_extracellular.tsv",
    "l1_routed_in_member_extracellular.tsv",
    "routed_in_ligand_contacts.tsv",
)


def _sha256(path: Path) -> str:
    import hashlib

    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _topology_snapshot() -> str:
    from deliverome_analysis.surfaceome import topology

    return str(getattr(topology, "SNAPSHOT_KEY", "unknown"))


#: Saved copies of the gitignored data/ inputs, one dated key each, with the SHA-256 recorded
#: in INPUTS_MANIFEST. `--restore` fetches and verifies them; nothing else downloads them.
INPUTS_MANIFEST = HERE / "l1_internal_decisions_inputs.tsv"
INPUTS_KEY_PREFIX = "sprint5/s5e13-tag-site-benchmarking/decision_inputs/2026-11-14/"


def missing_inputs() -> list[str]:
    return [f for f in INPUT_FILES if not (HERE / f).is_file()]


def require_inputs() -> None:
    """A missing input is an error, never a silently smaller run."""
    gone = missing_inputs()
    if gone:
        raise SystemExit(
            "missing inputs (run `l1_internal_decisions.py --restore` to fetch the saved "
            "copies):\n  " + "\n  ".join(gone)
        )


def restore_inputs() -> int:
    """Fetch every data/ input from its saved key and check it against the manifest."""
    from deliverome_analysis.core import s3

    man = pd.read_csv(INPUTS_MANIFEST, sep="\t", dtype=str).fillna("")
    bad = []
    for r in man.itertuples(index=False):
        local = HERE / r.path
        if local.is_file() and _sha256(local) == r.sha256:
            continue
        if not r.s3_key:
            bad.append(f"{r.path}: missing or changed, and has no saved key ({r.origin})")
            continue
        local.parent.mkdir(parents=True, exist_ok=True)
        s3.processed.download(r.s3_key, str(local))
        if _sha256(local) != r.sha256:
            bad.append(f"{r.path}: processed://{r.s3_key} does not match its SHA-256")
        else:
            print(f"restored {r.path}")
    if bad:
        raise SystemExit("\n".join(bad))
    require_inputs()
    print(f"all {len(INPUT_FILES)} inputs present and verified")
    return 0


def input_hashes() -> dict:
    require_inputs()
    return {f: _sha256(HERE / f) for f in INPUT_FILES}


def check() -> int:
    """Rebuild the decisions offline into CHECK_DIR and compare byte for byte."""
    global OFFLINE, OUT, OUT_JSON, REFOLD_PENDING, UNCONFIRMED_DESIGN, SINGLE_PASS_DESIGN
    global ROUTED_IN_DESIGN
    import shutil

    require_inputs()
    committed_json = json.loads(OUT_JSON.read_text())
    recorded = committed_json.get("inputs", {})
    now = input_hashes()
    changed = sorted(f for f in recorded if now.get(f) != recorded[f])
    snap = committed_json.get("pinned", {}).get("topology_snapshot")
    if snap and snap != _topology_snapshot():
        changed.append(f"topology snapshot {snap} -> {_topology_snapshot()}")
    CHECK_DIR.mkdir(parents=True, exist_ok=True)
    committed = {
        "decisions": OUT,
        "decisions_json": OUT_JSON,
        "requests": REFOLD_PENDING,
        "unconfirmed": HERE / UNCONFIRMED_DESIGN,
        "single_pass": SINGLE_PASS_DESIGN,
        "routed_in": ROUTED_IN_DESIGN,
    }
    OFFLINE = True
    OUT, OUT_JSON = CHECK_DIR / OUT.name, CHECK_DIR / OUT_JSON.name
    shutil.copy(REFOLD_PENDING, CHECK_DIR / REFOLD_PENDING.name)  # cumulative: start from it
    REFOLD_PENDING = CHECK_DIR / REFOLD_PENDING.name
    UNCONFIRMED_DESIGN = str((CHECK_DIR / Path(UNCONFIRMED_DESIGN).name).relative_to(HERE))
    SINGLE_PASS_DESIGN = CHECK_DIR / SINGLE_PASS_DESIGN.name
    ROUTED_IN_DESIGN = CHECK_DIR / ROUTED_IN_DESIGN.name
    build()
    rebuilt = {
        "decisions": OUT,
        "decisions_json": OUT_JSON,
        "requests": REFOLD_PENDING,
        "unconfirmed": HERE / UNCONFIRMED_DESIGN,
        "single_pass": SINGLE_PASS_DESIGN,
        "routed_in": ROUTED_IN_DESIGN,
    }
    bad = [k for k in committed if committed[k].read_bytes() != rebuilt[k].read_bytes()]
    print(f"inputs changed since the committed run: {changed or 'none'}")
    print(f"outputs differing from the committed files: {bad or 'none'}")
    return 1 if (bad or changed) else 0


def main(argv=None) -> int:
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--check",
        action="store_true",
        help="rebuild offline (cached refold verdicts, pinned inputs) and compare byte for byte",
    )
    ap.add_argument(
        "--restore",
        action="store_true",
        help="fetch the saved copies of the gitignored data/ inputs and verify their SHA-256",
    )
    args = ap.parse_args(argv)
    if args.restore:
        return restore_inputs()
    return check() if args.check else build()


def build() -> int:
    require_inputs()
    mask = membrane_mask.load()
    ids = mask.drop_duplicates("sequence_sha256").set_index("sequence_sha256")
    targets = pd.read_csv(HERE / "l1_adopted_targets.tsv", sep="\t")
    acc_sha = dict(zip(targets["uniprot_acc"], targets["sequence_sha256"], strict=False))
    designs = pd.read_csv(HERE / "tag_sites_v5_adopted.tsv", sep="\t")
    # 2026-11-14: 11 members redesigned on proteoform-fixed features (redesign_isoform_fix.py);
    # their designs replace the adopted-panel rows, then every rule below applies as usual
    fix = HERE / ISOFORM_FIX_DESIGN
    if fix.is_file():
        f = pd.read_csv(fix, sep="\t")
        designs = pd.concat(
            [designs[~designs["uniprot_acc"].isin(set(f["uniprot_acc"]))], f[designs.columns]],
            ignore_index=True,
        )
    designs["sequence_sha256"] = designs["uniprot_acc"].map(acc_sha)
    orth = ortholog_designs()
    unusable = mask[mask["usable"] != 1]
    unconfirmed = unusable[(unusable["status"] == "placed") & (unusable["n_eligible"] > 0)]
    unconf_accs = {a for a, s in acc_sha.items() if s in set(unconfirmed["sequence_sha256"])}
    unconf = design_unconfirmed(mask, unconf_accs)
    unconf["sequence_sha256"] = unconf["uniprot_acc"].map(acc_sha)

    tp = pd.read_csv(
        HERE / "data" / f"l1_tag_paths_adopted_{V5.ADOPTED_COMMIT}.tsv", sep="\t", low_memory=False
    )
    internal = tp[tp["placement"] == "internal"]
    rows = []
    for m in internal.itertuples(index=False):
        s = m.sequence_sha256
        base = {
            "member_id": m.candidate_id,
            "design_id": ids["design_id"].get(s),
            "gene": m.gene_symbol,
            "sequence_sha256": s,
            "candidate_set": m.candidate_set,
        }
        mrow = ids.loc[s] if s in ids.index else None
        d = designs[(designs["sequence_sha256"] == s) & (designs["rank"] == 1)]
        if d.empty:
            d = orth[(orth["sequence_sha256"] == s) & (orth["rank"] == 1)]
        u = unconf[(unconf["sequence_sha256"] == s) & (unconf["rank"] == 1)]
        if m.gene_symbol in SINGLE_PASS_OR_GPI and (mrow is None or mrow["usable"] != 1):
            rows.append(
                {
                    **base,
                    "decision": "route_to_nterm_after_sp",
                    "confidence": "high",
                    "reason": "the membrane-side gate finds fewer than two confident helices: a "
                    "single-pass or GPI-anchored protein has no two-sided loop to tag, and its "
                    "N-terminus after the signal peptide is extracellular; rung move to be applied "
                    "by the manifest writer",
                    "rejected_alternative": "internal loop on DeepTMHMM v1 O alone (the ranker refuses "
                    "targets without a usable mask)",
                }
            )
        elif not u.empty:
            r = u.iloc[0]
            flag = m.gene_symbol in INVERSION_SUSPECT
            rows.append(
                {
                    **base,
                    **site_row(
                        r,
                        "internal_loop_orientation_unconfirmed",
                        "the gate could not confirm the face; designed on its mask as computed, "
                        "because v1 is the topology of record"
                        + (
                            "; most v1 O residues fail the gate, a sign v1 may be inverted: "
                            "rung check requested"
                            if flag
                            else ""
                        ),
                        "snorkel construct (no builder in the library)",
                        "low",
                    ),
                }
            )
        elif mrow is not None and mrow["usable"] != 1:
            rows.append(
                {
                    **base,
                    "decision": "snorkel",
                    "confidence": "medium",
                    "reason": f"v1 extracellular exposure {int(mrow['n_v1_O'])} residues and the gate "
                    "rejects all of it; a snorkel tag (TM + tag on the cytoplasmic C-terminus, Brown "
                    "et al. 2013, PMID 24023844) is built for multipass proteins with small loops",
                    "rejected_alternative": "no tag; or an internal site the gate places in or below "
                    "the membrane",
                }
            )
        elif not d.empty:
            r = d.iloc[0]
            if bool(r["snorkel_candidate"]):
                reason = (
                    "v5 rank-1 eligible junction; kept on the loop although it is a snorkel "
                    "candidate: a quarter of loops hosting working published tags score below the "
                    "25th percentile, and the library has no snorkel builder. A failed refold "
                    "moves it to rank 2 or to a snorkel"
                )
                alt, conf = "snorkel construct", "low"
            else:
                reason = "v5 rank-1 eligible junction under membrane-side mask r4"
                alt, conf = "v5 rank-2 pick", "medium"
            rows.append({**base, **site_row(r, "internal_loop", reason, alt, conf)})
        else:
            rows.append(
                {
                    **base,
                    "decision": "unresolved",
                    "reason": "no design and no mask rule applies",
                }
            )
    out = pd.DataFrame(rows)
    out = rekey_to_head(out, designs, orth, unconf)
    out = apply_build_overrides(out, designs)
    out = apply_late_designs(out)
    out = apply_routed_in(out)
    out = gate_check(out)
    out = nterm_fallback(out)
    late = pd.read_csv(HERE / "data" / "tag_sites_v5_start_trims.tsv", sep="\t")
    single = pd.read_csv(SINGLE_PASS_DESIGN, sep="\t")
    single["sequence_sha256"] = single["uniprot_acc"].map(acc_sha)
    cols = ["sequence_sha256", "rank", "insert_after_residue", "rank_score", "absolute_score",
            "best_region_pct_of_published_regions", "this_region_length", "snorkel_candidate"]  # fmt: skip
    routed = pd.read_csv(ROUTED_IN_DESIGN, sep="\t")
    alternatives = pd.concat([x[cols] for x in (designs, orth, late, unconf, single, routed)])
    out["site_history"] = ""
    out, clean_alternatives = apply_contact_veto(out, alternatives)
    out = apply_refold(out, clean_alternatives)
    out = apply_snorkel_verdicts(out)
    out = apply_experimental_topology(out)
    out = out.assign(
        decided_by=DECIDED_BY,
        implemented_by=IMPLEMENTED_BY,
        decided_on="2026-11-12",
    )
    for did, r in ROUTED_IN.items():
        out.loc[out["design_id"] == did, "decided_by"] = r["decided_by"]
        out.loc[out["design_id"] == did, "decided_on"] = r["decided_on"]
    out.to_csv(OUT, sep="\t", index=False)
    OUT_JSON.write_text(
        json.dumps(
            {
                "membership_commit": MEMBERSHIP_COMMIT,
                "members": len(out),
                "decisions": out["decision"].value_counts().to_dict(),
                "left_internal_rung": LEFT,
                "decided_by": DECIDED_BY,
                "implemented_by": IMPLEMENTED_BY,
                "policy": POLICY,
                "policy_chosen_by": "compare_internal_policies.py (RF positive-only vs PU-RF vs "
                "PU-RF + veto vs v5 vs v5 + veto, 3 negative draws)",
                "contact_veto": out["contact_veto"].ne("").sum().item(),
                "pinned": {
                    "membership_commit": MEMBERSHIP_COMMIT,
                    "ortholog_designs_commit": ORTHOLOG_COMMIT,
                    "start_trim_designs_commit": START_TRIMS_COMMIT,
                    "gate_mask": GATE_R5_KEY,
                    "refold_verdicts": "refold_verdicts_seen.tsv (cumulative cache)",
                    "topology_snapshot": _topology_snapshot(),
                },
                "inputs": input_hashes(),
            },
            indent=2,
        )
        + "\n"
    )
    print(out["decision"].value_counts().to_string())
    print(out.groupby("decision")["confidence"].value_counts().to_string())
    assert (~out["decision"].isin(["unresolved"])).all(), out[out["decision"] == "unresolved"]
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
