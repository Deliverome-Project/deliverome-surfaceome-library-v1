#!/usr/bin/env python3
"""Internal-rung sites for members added to L1 after the internal decisions were built.

The co-receptor rule Becca adopted on 2026-11-21 admits 25 buildable genes ("Yes, rebuild
now"); one of them, SCNN1B (ENaC beta, l1-p3174), lands on the internal rung (DeepTMHMM v1: two
TM helices, both termini cytoplasmic, one 452-residue extracellular loop). It is designed here
the way every internal member was, and its rows are APPENDED to ``l1_internal_decisions.tsv``
and ``l1_internal_arm_decisions.tsv``; every existing row is asserted byte-identical, because
rerunning ``l1_internal_decisions.py`` / ``l1_internal_arms.py`` from scratch re-reads moving
inputs (the latest refold-verdict snapshots) that this addition must not drag in.

    features   ``ortholog_member_features.py --set coreceptor_additions``: the designed
               protein, its pinned DeepTMHMM v1 row, the AlphaFold DB v6 model of P51168
               (identical to the designed protein, residue for residue; AFDB served it on
               2026-11-21), AlphaFold-disorder and metapredict recomputed on Modal and equal
               byte for byte to D1's stored rows, NetSurfP-3.0 from D1
    design     v5 unchanged (``tag_site_model_v4.design`` with the v5 models), as for the
               ROUTED_IN members: without the r8 membrane mask, which has no row for this
               sequence, and kept to the v1 extracellular span (the exported residues are
               exactly the v1 ``O`` residues)
    veto       the +/-2 ligand-contact veto on D1 ligand contacts and contacts transferred by
               ``transfer_ligand_contacts`` (``routed_in_contacts.contact_rows``, the routed-in
               members' recipe); protein-partner contacts (here the ENaC alpha/gamma interfaces)
               are reported, never vetoed, as for every member
    arms       the order rule of ``l1_internal_arms.py``: the design's
               ``recommended_sites_to_order`` (best region ceiling against published regions),
               capped at three, sites at least ``MIN_SEPARATION`` apart; a second or third site
               must pass refold. No refold has run yet, so every candidate is listed for the one
               refold batch of the additions and every arm row carries ``refold_verdict``
               pending; arms after the first are confirmed or dropped when verdicts arrive.

Decided by: Claude Code session l1-review-export-scope (subagent), under Becca's 2026-11-21
"Yes, rebuild now"; the rules are the ranker sessions' (PR #245) and the coordinator's adoption
of the order rule (2026-11-15).

**Member sets** (``SETS``, ``--set``). ``coreceptor_additions`` is the set above.
``additions_2026_11_22`` is the internal-rung members of the 2026-11-22 additions (round-2 genes,
the isoform redraw and its refills, the ortholog restorations) and the members whose reference
moved to the GTF MANE Select protein (``mane_gtf_changes_2026-11-22.tsv``), with features from
``ortholog_member_features.py --set additions_2026_11_22``. A member already decided on ANOTHER
protein (a MANE switch) has its rows replaced: rows are keyed by ``sequence_sha256``, the old
rows leave the tables and the new go at the end; every other member's row is kept byte for byte.
A member already decided on THIS protein is skipped.

**Deltas** (``--apply-delta FILE``). A decision delta from another session (rows in the two
tables' columns, ``table`` = ``decisions`` or ``arms``) replaces the named members' rows the
same way, so a later redesign can land on top of these appends.

**Refold verdicts** (``--apply-verdicts``) follow Becca's 2026-11-22 rule ("Yes, build them":
refold failures flag, never hold): a first arm that fails on fold or site lDDT stays built with
``flag`` = ``refold fail: <reasons>``; a mechanistic failure (an epitope packed on the protein,
a confident helix inside the cassette, a buried site) keeps the rung's rule and is held for the
next site; a second or third arm must pass (order rule) and is dropped otherwise.

    uv run --locked python l1_internal_additions.py --set S          # append rows, write requests
    uv run --locked python l1_internal_additions.py --set S --check  # recompute, no writes
    uv run --locked python l1_internal_additions.py --set S --apply-verdicts  # after the refold
    uv run --locked python l1_internal_additions.py --apply-delta delta.tsv
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import joblib
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[3] / "src"))

import tag_site_model_v4 as V4  # noqa: E402
import tag_site_model_v5 as V5  # noqa: E402

from deliverome_analysis import plasmid_qc as qc  # noqa: E402

PANEL = HERE.parents[1] / "l1" / "orf-design"
DECISIONS = HERE / "l1_internal_decisions.tsv"
ARMS = HERE / "l1_internal_arm_decisions.tsv"
MAX_ARMS = 3
FLAG_BY = "Becca Carlson, 2026-11-22 ('Yes, build them'): refold failures flag, never hold"

#: Per member set: its feature tables, outputs and attribution.
SETS = {
    "coreceptor_additions": dict(
        features="l1_coreceptor_addition_member_extracellular.tsv",
        provenance=HERE / "l1_coreceptor_addition_member_features.tsv",
        targets=HERE / "data" / "coreceptor_addition_targets.tsv",
        design=HERE / "tag_sites_v5_coreceptor_additions.tsv",
        contacts=HERE / "coreceptor_addition_ligand_contacts.tsv",
        requests=HERE / "refold_requests_coreceptor_additions.tsv",
        decided_on="2026-11-21",
        decided_by=(
            "Claude Code session l1-review-export-scope (subagent), under Becca's 2026-11-21 "
            '"Yes, rebuild now"; site rules of PR #245, order rule adopted by the coordinator '
            "2026-11-15"
        ),
        implemented_by="Claude Code session l1-review-export-scope (subagent), 2026-11-21",
        why="co-receptor addition (Becca, 2026-11-21)",
        batch="tagged/2026-11-21-coreceptor",
    ),
    "additions_2026_11_22": dict(
        features="l1_additions_2026_11_22_member_extracellular.tsv",
        provenance=HERE / "l1_additions_2026_11_22_member_features.tsv",
        targets=HERE / "data" / "additions_2026_11_22_targets.tsv",
        design=HERE / "tag_sites_v5_additions_2026_11_22.tsv",
        contacts=HERE / "additions_2026_11_22_ligand_contacts.tsv",
        requests=HERE / "refold_requests_additions_2026_11_22.tsv",
        decided_on="2026-11-22",
        decided_by=(
            "Claude Code session l1-review-export-scope (subagent), under Becca's 2026-11-21 "
            '"Yes, rebuild now" and 2026-11-22 "Yes, in this build"; site rules of PR #245, '
            "order rule adopted by the coordinator 2026-11-15"
        ),
        implemented_by="Claude Code session l1-review-export-scope (subagent), 2026-11-22",
        why="2026-11-22 addition or GTF MANE Select switch (Becca, 2026-11-21/22)",
        batch="tagged/2026-11-22-all-rungs",
    ),
}
#: Set by configure(); the module-level names the rest of the script reads.
FEATURES = PROVENANCE = TARGETS = DESIGN = CONTACTS = REQUESTS = None
DECIDED_ON = DECIDED_BY = IMPLEMENTED_BY = WHY = BATCH = None


def configure(name: str) -> None:
    global FEATURES, PROVENANCE, TARGETS, DESIGN, CONTACTS, REQUESTS
    global DECIDED_ON, DECIDED_BY, IMPLEMENTED_BY, WHY, BATCH
    s = SETS[name]
    FEATURES, PROVENANCE, TARGETS = s["features"], s["provenance"], s["targets"]
    DESIGN, CONTACTS, REQUESTS = s["design"], s["contacts"], s["requests"]
    DECIDED_ON, DECIDED_BY, IMPLEMENTED_BY = s["decided_on"], s["decided_by"], s["implemented_by"]
    WHY, BATCH = s["why"], s["batch"]


def model_note(m) -> str:
    """Which model the member's features were computed on, from the provenance row."""
    if m.model_source == "alphafold_db":
        acc = str(m.model_ref).rsplit("AF-", 1)[-1].split("-F1-")[0]
        return f"features on the AlphaFold DB v6 model of {acc}, identical to the designed protein"
    if str(m.model_source).startswith("colabfold"):
        return f"features on our ColabFold model of the designed protein ({m.model_ref})"
    return "features without a structure model (sequence predictors only)"


def replace_rows(path: Path, new: pd.DataFrame, cols, ids: set | None = None) -> int:
    """Drop the rows of the members in ``new`` (or ``ids``) and append ``new``; every other row
    verbatim."""
    old = path.read_text().splitlines()
    col = old[0].split("\t").index("design_id")
    ids = set(new["design_id"]) if ids is None else set(ids)
    kept = [ln for ln in old[1:] if ln.split("\t")[col] not in ids]
    table = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    out = pd.concat([table[~table["design_id"].isin(ids)], new[cols]], ignore_index=True)
    out.to_csv(path, sep="\t", index=False, lineterminator="\n")
    now = path.read_text().splitlines()
    assert now[0] == old[0] and now[1 : 1 + len(kept)] == kept, f"{path.name}: another row changed"
    return len(old) - 1 - len(kept)


CASSETTE_AA = qc.P003_CASSETTE_PEPTIDE


def members() -> pd.DataFrame:
    p = pd.read_csv(PROVENANCE, sep="\t")
    # a member with no UniProt accession of its own sequence has none (not NaN)
    p["sequence_accession"] = p["sequence_accession"].fillna("").astype(str)
    tp = pd.read_csv(PANEL / "l1_tag_paths.tsv", sep="\t")
    p = p.merge(
        tp[["candidate_id", "placement", "candidate_set"]],
        left_on="design_id",
        right_on="candidate_id",
        validate="one_to_one",
    )
    assert (p["placement"] == "internal").all(), p[["design_id", "placement"]]
    act = pd.read_csv(PANEL / "candidate_manifest_tranches.tsv", sep="\t", dtype=str).fillna("")
    retired = set(act.loc[act["retired_reason"] != "", "candidate_id"])
    return p[~p["design_id"].isin(retired)].reset_index(drop=True)


def protein_of(m) -> str:
    from deliverome_analysis import library_design as ld

    prot = ld.designed_proteins(PANEL / "candidate_manifest_tranches.tsv", library=ld.LIBRARY_L1)
    seq = prot[m["design_id"]]["sequence"]
    assert hashlib.sha256(seq.encode()).hexdigest() == m["sequence_sha256"], m["design_id"]
    return seq


#: UniProt 2026_04 SIGNAL/PROPEP/CHAIN features, pinned by the N-terminal package (its copy of
#: external://uniprot/2026_04/nterm-propeptides/uniprot_features.tsv).
UNIPROT_FEATURES = PANEL / "reviewer-package/nterm/uniprot_features.tsv"


def mature_only(mem: pd.DataFrame) -> str:
    """The set's feature table without residues inside a UniProt propeptide (cleaved in the
    secretory pathway, so a tag there leaves the mature protein), for members whose sequence IS
    that entry's. Written beside the data and returned as the table the design reads.

    Added 2026-11-22: CDH3's GTF MANE Select form (P22223, PROPEP 25..107) drew its rank-1 site
    after residue 86, inside the propeptide."""
    import re

    u = pd.read_csv(UNIPROT_FEATURES, sep="\t", dtype=str).fillna("")
    seq = {r["design_id"]: protein_of(r) for _, r in mem.iterrows()}
    feat = pd.read_csv(HERE / FEATURES, sep="\t", low_memory=False)
    drop = pd.Series(False, index=feat.index)
    for m in mem.itertuples():
        e = u[(u["Entry"] == m.sequence_accession) & (u["Sequence"] == seq[m.design_id])]
        for prop in e["Propeptide"]:
            for lo, hi in re.findall(r"PROPEP (\d+)\.\.(\d+)", prop):
                hit = (feat["uniprot_acc"] == m.key) & feat["residue"].between(int(lo), int(hi))
                print(
                    f"{m.design_id} {m.gene}: propeptide {lo}-{hi} excluded ({hit.sum()} residues)"
                )
                drop |= hit
    out = HERE / "data" / f"{Path(FEATURES).stem}_mature.tsv"
    feat[~drop].to_csv(out, sep="\t", index=False)
    return str(out.relative_to(HERE))


def design(mem: pd.DataFrame) -> pd.DataFrame:
    global FEATURES
    FEATURES = mature_only(mem)
    TARGETS.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(
        {
            "gene": mem["gene"],
            "uniprot_acc": mem["key"],
            "member_id": mem["design_id"],
            "candidate_set": mem["candidate_set"],
            "topology_cohort": mem["topology_cohort"],
            "sequence_sha256": mem["sequence_sha256"],
        }
    ).to_csv(TARGETS, sep="\t", index=False)
    spec = json.loads(V5.SPEC.read_text())
    d = V4.design(
        joblib.load(V5.MODEL_Z),
        joblib.load(V5.MODEL_RAW),
        spec,
        np.array(spec["published_protein_ceilings_oof"]),
        np.array(spec["published_region_ceilings_oof"]),
        str(TARGETS.relative_to(HERE)),
        str(DESIGN.relative_to(HERE)),
        V5.EXTRA,
        use_mask=False,  # no r8 row for these sequences; the face is v1's extracellular span
        skip_refused=True,
        min_region_len=V5.MIN_REGION_LEN,
        override_tables=(FEATURES,),
    )
    feat = pd.read_csv(HERE / FEATURES, sep="\t", usecols=["uniprot_acc", "residue"])
    ecto = feat.groupby("uniprot_acc")["residue"].agg(["min", "max"])
    inside = d.apply(
        lambda r: (
            ecto.loc[r["uniprot_acc"], "min"]
            <= r["insert_after_residue"]
            < ecto.loc[r["uniprot_acc"], "max"]
        ),
        axis=1,
    )
    assert inside.all(), d[~inside]
    d["sequence_sha256"] = d["uniprot_acc"].map(
        dict(zip(mem["key"], mem["sequence_sha256"], strict=True))
    )
    d.to_csv(DESIGN, sep="\t", index=False)
    return d


def ligand_contacts(mem: pd.DataFrame, seqs: dict[str, str]) -> pd.DataFrame:
    import routed_in_contacts as RC

    rows = RC.contact_rows(
        [
            (m.design_id, m.gene, m.sequence_sha256, m.sequence_accession, seqs[m.design_id])
            for m in mem.itertuples()
        ]
    )
    c = pd.DataFrame(rows)
    c.to_csv(CONTACTS, sep="\t", index=False)
    return c


def vetoed(c: pd.DataFrame, did: str, n: int) -> list[str]:
    """Ligand contacts within +/-2 of the junction (n | n+1), as contact_veto.vetoed reads them."""
    if c.empty:
        return []
    mine = c[c["design_id"] == did]
    hit = mine[(mine["residue"] >= n - 1) & (mine["residue"] <= n + 2)]
    return sorted({f"{r.source}:{int(r.residue)}" for r in hit.itertuples()})


def partner_contacts(acc: str) -> pd.DataFrame:
    """PDB interface residues of this accession with other proteins (evidence, never vetoed)."""
    ic = pd.read_csv(HERE / "interface_contacts.tsv", sep="\t")
    cols = list(ic.columns)
    acc_col = next(c for c in ("uniprot_acc", "member_acc", "acc") if c in cols)
    return ic[ic[acc_col] == acc]


def _member_arms(
    m,
    seq: str,
    chosen: list[dict],
    n_rec: int,
    ceiling: float,
    why_count: str,
    reason: str,
    arms: pd.DataFrame,
    issued_arms: set,
) -> tuple[int, list[dict]]:
    """The member's loop arms by the order rule, numbered after every arm id it has had."""
    did, sha = m.design_id, m.sequence_sha256
    new_arms = []
    # arm ids are permanent and never reused: a member decided again on a new protein
    # numbers its arms after every arm id it has had (per-arm table, issued encodings)
    k0 = max(
        [int(a.rsplit("-a", 1)[1]) for a in arms["arm_id"] if a.startswith(f"{did}-a")]
        + [int(a.rsplit("-a", 1)[1]) for a in issued_arms if a.startswith(f"{did}-a")]
        + [0]
    )
    for k, r in enumerate(chosen[:n_rec], start=1):
        xr = r["x"]
        role = "loop_best" if k == 1 else f"loop_{k}"
        arm = {col: "" for col in arms.columns}
        arm.update(
            {
                "arm_id": f"{did}-a{k0 + k}",
                "construct_key": f"{did}:internal{r['n']}",
                "design_id": did,
                "member_id": did,
                "gene": m.gene,
                "sequence_sha256": sha,
                "arm_role": role,
                "arm_kind": "loop",
                "insert_after_residue": str(r["n"]),
                "decision": "build",
                "order_verdict": f"{n_rec} site(s)",
                "order_why": why_count,
                "model_version": "v5",
                "design_rank": str(r["rank"]),
                "region_index": str(r["region"]),
                "region_ceiling": f"{ceiling:.6g}",
                "v5_rank_score": f"{float(xr.rank_score):.6g}",
                "v5_absolute_score": f"{float(xr.absolute_score):.6g}",
                "contact_veto": "False",
                "gate": "unknown (no r8 row; v1 extracellular span)",
                "refold_verdict": "pending",
                "refold_tagged_sha256": hashlib.sha256(
                    (seq[: r["n"]] + CASSETTE_AA + seq[r["n"] :]).encode()
                ).hexdigest(),
                "flag": ""
                if k == 1
                else "kept only if it passes refold (order rule: a second or third site "
                "must pass refold)",
                "reason": reason if k == 1 else f"order rule site {k} of {n_rec}; {why_count}",
                "decided_by": DECIDED_BY,
                "implemented_by": IMPLEMENTED_BY,
                "decided_on": DECIDED_ON,
            }
        )
        new_arms.append(arm)
    return k0, new_arms


def _refold_requests(m, seq: str, chosen: list[dict], n_rec: int, k0: int) -> list[dict]:
    """One fold batch for the member: the order rule's sites plus two spares."""
    did, sha = m.design_id, m.sequence_sha256
    requests = []
    for r in chosen[: n_rec + 2]:  # the order rule's sites plus two spares, one fold batch
        t = seq[: r["n"]] + CASSETTE_AA + seq[r["n"] :]
        requests.append(
            {  # l1_refold_alternatives_requested.tsv's columns, then the extras
                "design_id": did,
                "gene": m.gene,
                "sequence_sha256": sha,
                "rank": r["rank"],
                "insert_after_residue": r["n"],
                "tagged_protein_sha256": hashlib.sha256(t.encode()).hexdigest(),
                "why": "order-rule site" if chosen.index(r) < n_rec else "spare",
                "still_pending": True,
                "arm_id": f"{did}-a{k0 + chosen.index(r) + 1}" if chosen.index(r) < n_rec else "",
                "native_model": f"{m.model_source} {m.model_ref} (identical sequence)",
                "tagged_protein": t,
            }
        )
    return requests


def build(check: bool) -> int:
    mem = members()
    seqs = {m["design_id"]: protein_of(m) for _, m in mem.iterrows()}
    d = design(mem)
    c = ligand_contacts(mem, seqs)
    spec = json.loads(V5.SPEC.read_text())
    pub = np.array(spec["published_region_ceilings_oof"])
    q25, med = float(np.quantile(pub, 0.25)), float(np.median(pub))

    dec = pd.read_csv(DECISIONS, sep="\t", dtype=str, keep_default_na=False)
    arms = pd.read_csv(ARMS, sep="\t", dtype=str, keep_default_na=False)
    issued_arms = set(
        pd.read_csv(PANEL / "l1_arm_cassette_encodings.tsv", sep="\t", usecols=["arm_id"])["arm_id"]
    )
    new_dec, new_arms, requests = [], [], []
    for m in mem.itertuples():
        did, sha, seq = m.design_id, m.sequence_sha256, seqs[m.design_id]
        same = dec[(dec["design_id"] == did) & (dec["sequence_sha256"] == sha)]
        if len(same):  # decided on this protein already: nothing to do
            print(f"{did} {m.gene}: already decided on this protein, skipped")
            continue
        replaced = dec.loc[dec["design_id"] == did, "sequence_sha256"].tolist()
        g = d[d["sequence_sha256"] == sha].sort_values("rank")
        top = g.iloc[0]
        feat = pd.read_csv(HERE / FEATURES, sep="\t", usecols=["uniprot_acc", "residue"])
        span = feat.groupby("uniprot_acc")["residue"].agg(lambda r: (int(r.min()), int(r.max())))
        n_rec = min(int(top["recommended_sites_to_order"]), MAX_ARMS)
        ceiling = float(top["region_ceiling"])
        clean, chosen = [], []
        for x in g.itertuples(index=False):
            n = int(x.insert_after_residue)
            v = vetoed(c, did, n)
            row = dict(n=n, rank=int(x.rank), region=x.region_index, veto=v, x=x)
            clean.append(row)
            if v:
                continue
            if all(
                r["region"] != row["region"] or abs(n - r["n"]) >= V4.MIN_SEPARATION for r in chosen
            ):
                chosen.append(row)
        if not chosen:
            # Every designed site has a ligand contact at the junction. No loop arm; a snorkel is
            # never an arm alone (Becca, 2026-11-15), so the member is held for review.
            new_dec.append(
                {c: "" for c in dec.columns}
                | {
                    "member_id": did, "design_id": did, "gene": m.gene,
                    "candidate_set": m.candidate_set, "decision": "hold_no_passing_arm",
                    "confidence": "low", "sequence_sha256": sha,
                    "reason": f"{WHY}: every v5-designed site has a ligand contact within 2 "
                    "residues of the junction; no loop arm, and a snorkel is never an arm alone",
                    "contact_veto": "; ".join(
                        f"rank {r['rank']} after {r['n']} vetoed: {', '.join(r['veto'])}"
                        for r in clean
                    ),
                    "contact_sources_available": "d1_ligand, transferred_ligand (routed-in recipe)",
                    "decided_by": DECIDED_BY, "implemented_by": IMPLEMENTED_BY,
                    "decided_on": DECIDED_ON,
                }
            )  # fmt: skip
            print(f"{did} {m.gene}: every designed site contact-vetoed; held")
            continue
        pick = chosen[0]
        x = pick["x"]
        pct = float(getattr(x, "best_region_pct_of_published_regions", np.nan))
        why_count = (
            f"best region ceiling {ceiling:.3f} against published regions (median {med:.4f}, "
            f"q25 {q25:.4f}) -> {int(top['recommended_sites_to_order'])} site(s) by the order rule"
        )
        reason = (
            f"{WHY}: v5 rank-{pick['rank']} site in the "
            f"extracellular span {span[m.key][0]}-{span[m.key][1]}, scored without the r8 mask "
            f"(no row for this sequence; the span is DeepTMHMM v1's extracellular residues); "
            f"{model_note(m)}; {why_count}"
        )
        if replaced:
            reason += (
                f"; replaces the decision made on its earlier protein "
                f"(sequence_sha256 {replaced[0][:12]}), whose reference moved"
            )
        alt = chosen[1] if len(chosen) > 1 else None
        new_dec.append(
            {
                "member_id": did,
                "design_id": did,
                "gene": m.gene,
                "candidate_set": m.candidate_set,
                "decision": "internal_loop",
                "insert_after_residue": str(pick["n"]),
                "pick_rank": str(pick["rank"]),
                "rank_score": f"{float(x.rank_score):.6g}" if hasattr(x, "rank_score") else "",
                "absolute_score": f"{float(x.absolute_score):.6g}"
                if hasattr(x, "absolute_score")
                else "",
                "best_region_pct_of_published_regions": "" if np.isnan(pct) else f"{pct:g}",
                "this_region_length": str(getattr(x, "region_length", "")),
                "snorkel_candidate": "False",
                "confidence": "low",
                "reason": reason,
                "rejected_alternative": (
                    f"rank-{alt['rank']} site after {alt['n']}" if alt else "no second clean site"
                ),
                "sequence_sha256": sha,
                "site_history": "",
                "contact_veto": "; ".join(
                    f"rank {r['rank']} after {r['n']} vetoed: {', '.join(r['veto'])}"
                    for r in clean
                    if r["veto"]
                )
                or "none",
                "contact_sources_available": "d1_ligand, transferred_ligand (routed-in recipe)",
                "refold_verdict": "pending",
                "snorkel_verdict": "",
                "decided_by": DECIDED_BY,
                "implemented_by": IMPLEMENTED_BY,
                "decided_on": DECIDED_ON,
            }
        )
        k0, member_arms = _member_arms(
            m, seq, chosen, n_rec, ceiling, why_count, reason, arms, issued_arms
        )
        new_arms += member_arms
        requests += _refold_requests(m, seq, chosen, n_rec, k0)
        partners = partner_contacts(m.sequence_accession)  # member-side rows; often none
        print(
            f"{did} {m.gene}: ceiling {ceiling:.3f} (median {med:.4f}, q25 {q25:.4f}) -> "
            f"{n_rec} site(s); chosen {[(r['n'], r['rank']) for r in chosen[:n_rec]]}; "
            f"vetoed {[(r['n'], r['veto']) for r in clean if r['veto']]}; "
            f"{len(c)} ligand contacts, {len(partners)} PDB interface rows (not vetoed)"
        )

    req = pd.DataFrame(requests)
    if check or not new_dec:
        if new_dec:
            print(pd.DataFrame(new_dec).T.to_string())
        return 0
    n_dec = replace_rows(DECISIONS, pd.DataFrame(new_dec), list(dec.columns))
    n_arm = replace_rows(
        ARMS,
        pd.DataFrame(new_arms, columns=list(arms.columns)),
        list(arms.columns),
        ids={r["design_id"] for r in new_dec},
    )
    req.to_csv(REQUESTS, sep="\t", index=False, lineterminator="\n")
    print(
        f"appended {len(new_dec)} decision(s), {len(new_arms)} arm(s) (replacing {n_dec} "
        f"decision and {n_arm} arm rows of earlier proteins); {len(req)} refold requests"
    )
    return 0


def apply_delta(path: Path) -> int:
    """Replace members' rows with a decision delta (column ``table``: decisions or arms)."""
    d = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    for table, dest in (("decisions", DECISIONS), ("arms", ARMS)):
        rows = d[d["table"] == table].drop(columns="table")
        if rows.empty:
            continue
        cols = list(pd.read_csv(dest, sep="\t", nrows=0).columns)
        assert set(rows.columns) <= set(cols), f"{path.name}: unknown {table} columns"
        rows = rows.reindex(columns=cols, fill_value="")
        print(f"{dest.name}: {replace_rows(dest, rows, cols)} row(s) replaced, {len(rows)} added")
    return 0


def mechanistic(reasons: str) -> bool:
    """A failure that says the tag itself is wrong there, not that the fold scored low."""
    return "epitope" in reasons or "helix inside the cassette" in reasons or "buried" in reasons


def apply_verdicts() -> int:
    """Fill the set's refold cells from ``refold_verdicts_seen.tsv`` (the cumulative cache).

    Only rows still ``pending`` and appended by this set change. A pass is recorded as ``pass``
    with its source. Under Becca's 2026-11-22 rule a first arm failing on fold or site lDDT is
    built with a flag; a mechanistic failure, or a failed second or third arm, is held (the
    member's next site is then the rung's call). Applied 2026-11-21 after the co-receptor refold
    batch (processed://.../predicted-structures/tagged/2026-11-21-coreceptor/).
    """
    seen = pd.read_csv(HERE / "refold_verdicts_seen.tsv", sep="\t", dtype=str).fillna("")
    by_sha = {r.tagged_sha256: r for r in seen.itertuples()}
    dec = pd.read_csv(DECISIONS, sep="\t", dtype=str, keep_default_na=False)
    arms = pd.read_csv(ARMS, sep="\t", dtype=str, keep_default_na=False)
    mine = arms["decided_by"].eq(DECIDED_BY) & arms["refold_verdict"].eq("pending")
    n = 0
    for i in arms.index[mine]:
        v = by_sha.get(arms.at[i, "refold_tagged_sha256"])
        if v is None:
            continue
        src = f"s4e9 tag_refold_verdicts:{v.source_id} ({BATCH})"
        first = arms.at[i, "arm_role"] == "loop_best"
        if v.verdict == "pass":
            disp, flag = "pass", arms.at[i, "flag"]
        elif first and not mechanistic(v.verdict_reasons):
            disp, flag = "flag_built", f"refold fail: {v.verdict_reasons} (built; {FLAG_BY})"
        else:
            disp, flag = "held", f"refold fail: {v.verdict_reasons}"
            arms.at[i, "decision"] = "hold"
        arms.loc[i, ["refold_verdict", "refold_reasons", "refold_source", "refold_disposition"]] = [
            v.verdict, v.verdict_reasons if v.verdict != "pass" else "", src, disp,
        ]  # fmt: skip
        arms.at[i, "flag"] = flag
        if first:
            dec.loc[dec["design_id"].eq(arms.at[i, "design_id"]), "refold_verdict"] = v.verdict
        n += 1
    # A member of this set left with no arm to build (its only site failed mechanistically) is
    # held, as the rung holds any member without a passing or flaggable arm.
    for i in dec.index[dec["decided_by"].eq(DECIDED_BY) & dec["decision"].eq("internal_loop")]:
        did = dec.at[i, "design_id"]
        if not (arms["design_id"].eq(did) & arms["decision"].eq("build")).any():
            dec.at[i, "decision"] = "hold_no_passing_arm"
            dec.at[i, "reason"] += "; its only designed site failed refold mechanistically: held"
            print(f"{did} {dec.at[i, 'gene']}: no arm left to build, held")
    dec.to_csv(DECISIONS, sep="\t", index=False, lineterminator="\n")
    arms.to_csv(ARMS, sep="\t", index=False, lineterminator="\n")
    print(f"applied {n} refold verdict(s)")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--set", default="coreceptor_additions", choices=sorted(SETS))
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--apply-verdicts", action="store_true")
    ap.add_argument("--apply-delta", type=Path)
    a = ap.parse_args()
    configure(a.set)
    if a.apply_delta:
        return apply_delta(a.apply_delta)
    return apply_verdicts() if a.apply_verdicts else build(a.check)


if __name__ == "__main__":
    raise SystemExit(main())
