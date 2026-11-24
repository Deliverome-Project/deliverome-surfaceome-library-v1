#!/usr/bin/env python3
"""Rescue round 2026-11-14: internal ectodomain sites for held L1 members not on the internal rung.

Becca, 2026-11-14: "let's find replacements where things were held if you can and it makes
sense." The coordinator, under Becca's 2026-11-12 delegation, asked this session for the held
members that are NOT on the internal rung (`rescue_round_members.tsv`, from the coordinator's
l1-held-members-rescue-round-2026-11-14 list; VTCN1 human and cyno and SLC41A2 are handled
elsewhere and excluded there). For each member:

    features   its own ordered protein, built by the ortholog recipe
               (`ortholog_member_features.py --set rescue`): identical AlphaFold DB model or our
               ColabFold fold, stored D1 predictor rows where the sequence is identical
    span       the residues a site may sit in, from UniProt topology of the identical entry,
               or of the gene's human canonical mapped by global alignment:
                 extracellular topological domains, minus the signal peptide, any cleaved
                 propeptide, and the first SP_MARGIN residues of the mature chain (a site
                 there is the held N-after-SP placement);
                 GPI-anchored: only N-terminal of the omega site (gpi_set.tsv at branch
                 claude/l1-cterm-audit-followups-2, ad4bdaec);
                 type II stalk members (CD70, BRICD5, CD209, SCARA5): the folded extracellular
                 domain only (UniProt Domain features, collagen-like and coiled-coil excluded),
                 since their stalk sites already failed refold
    score      v5 unchanged: both forests on the member's features, z-scored within protein and
               smoothed over +/-2 exactly as tag_site_model_v4.design does. Rank 1 is the best
               rank-model junction in the region with the highest absolute-model ceiling; rank 2
               the best junction in the next region long enough to hold one, else a second
               junction at least MIN_SEPARATION from rank 1 in the same region
    veto       no site within +/-2 of a ligand contact (own D1 contacts on an identical
               accession, plus contacts transferred at identity >= 0.5), as for every member
    skip       a member with no junction inside its span, or whose span holds no region of at
               least MIN_REGION_LEN residues (v5's loop floor)

v5 was trained on internal loops of multipass proteins. A site on a single-pass or GPI
ectodomain is outside its training distribution: the score orders residues by the same
features, but the evidence that a site works rests on the membrane-side gate and the refold.
Every output row says so.

    uv run --locked python rescue_round.py          # writes rescue_round_sites.tsv and .json
"""

from __future__ import annotations

import hashlib
import json
import sys
import warnings
from pathlib import Path

import joblib
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[3] / "src"))

import tag_site_model_v4 as V4  # noqa: E402
import tag_site_model_v5 as V5  # noqa: E402
from build_internal_tagged import CASSETTE_AA  # noqa: E402
from independent_model_search import within_protein  # noqa: E402

MEMBERS = HERE / "rescue_round_members.tsv"
FEATURES = HERE / "l1_rescue_member_extracellular.tsv"
PROVENANCE = HERE / "l1_rescue_member_features.tsv"
UNIPROT = HERE / "data" / "uniprot_rescue"  # cached UniProt JSON (fetched 2026-11-14)
GPI_SET = HERE / "data" / "gpi_set_cterm_audit_2.tsv"  # gpi_set.tsv at ad4bdaec
GPI_SET_COMMIT = "ad4bdaec"
OUT = HERE / "rescue_round_sites.tsv"
CONTACTS = HERE / "rescue_round_ligand_contacts.tsv"
OUT_JSON = HERE / "rescue_round_sites.json"
STALK = {"CD70", "BRICD5", "CD209", "SCARA5"}
#: not a folded domain: a type II stalk member keeps only its globular domain(s)
UNFOLDED_DOMAIN = ("collagen", "coiled")
#: residues of the mature N-terminus excluded: a site this close to the signal-peptide
#: cleavage is the held N-after-SP placement under another name, and can alter cleavage
SP_MARGIN = 10
#: VTCN1 (coordinator, 2026-11-14): the N-terminal segment's topology is unresolved. SignalP 6
#: calls no signal peptide, both DeepTMHMM versions call ~7-29 a signal anchor, and UniProt's
#: signal peptide 1-24 and TM 260-280 are predictions; one report (PMID 12818166, mouse,
#: unverified here) describes it as GPI-anchored. The site must be outside under either reading,
#: so it is kept inside the Ig domains, 10 residues past the anchor's end (29) and at least 30
#: residues N-terminal of 259, ahead of any plausible omega site. Q7Z7D3 numbering.
VTCN1_WINDOW = (40, 229)
TOPOLOGY_NOTE = {
    "VTCN1": "membrane topology of the N-terminal segment is unresolved (cleaved signal peptide "
    "or signal anchor; C-terminal TM or possibly GPI); the site lies in the Ig domains within "
    "Q7Z7D3 40-229, outside under either reading"
}
OOD = (
    "v5 was trained on internal loops of multipass proteins; a single-pass or GPI ectodomain "
    "site is outside its training distribution, so this site rests on the membrane-side gate "
    "and the refold"
)
ATTRIBUTION = (
    "rescue round, coordinator under Becca's 2026-11-12 delegation and her 2026-11-14 request"
)


def _aligner():
    from Bio import Align
    from Bio.Align import substitution_matrices

    al = Align.PairwiseAligner(mode="global")
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    al.open_gap_score, al.extend_gap_score = -10, -0.5
    return al


def position_map(ref: str, member: str) -> dict[int, int]:
    """1-based ref residue -> member residue, for every aligned (not gapped) pair."""
    if ref == member:
        return {i: i for i in range(1, len(ref) + 1)}
    a = _aligner().align(ref, member)[0]
    out = {}
    for (r0, r1), (m0, _m1) in zip(*a.aligned, strict=False):
        for k in range(r1 - r0):
            out[r0 + k + 1] = m0 + k + 1
    return out


def uniprot(acc: str) -> dict | None:
    p = UNIPROT / f"{acc}.json"
    return json.loads(p.read_text()) if p.is_file() else None


def features_of(rec: dict) -> dict:
    """Extracellular segments, signal end, GPI omega and domains from a UniProt record."""
    f = {"ext": [], "tm": [], "signal_end": 0, "omega": None, "domains": [], "chain": None,
         "propeptides": []}  # fmt: skip
    for x in rec.get("features", []):
        t, loc = x["type"], x["location"]
        s, e = loc["start"].get("value"), loc["end"].get("value")
        d = (x.get("description") or "").lower()
        if s is None or e is None:
            continue
        if t == "Topological domain" and d.startswith("extracellular"):
            f["ext"].append((s, e))
        elif t == "Transmembrane":
            f["tm"].append((s, e))
        elif t == "Signal":
            f["signal_end"] = max(f["signal_end"], e)
        elif t == "Lipidation" and "gpi" in d:
            f["omega"] = s
        elif t == "Propeptide":
            f["propeptides"].append((s, e))
        elif t == "Domain":
            f["domains"].append((s, e, x.get("description") or ""))
        elif t == "Chain" and f["chain"] is None:
            f["chain"] = (s, e)
    return f


def span(m: pd.Series, seq: str, hum: dict, gpi: pd.DataFrame) -> dict:
    """Allowed residues on the member, with the evidence used to get them."""
    acc = str(m["identical_acc"])
    own = None if acc.startswith("l1-") else uniprot(acc.split(";")[0])
    use_own = own is not None and own.get("sequence", {}).get("value") == seq
    if use_own and not (features_of(own)["ext"] or features_of(own)["omega"]):
        use_own = False  # an entry without topology (TrEMBL) cannot place the ectodomain
    ref_acc = acc if use_own else hum.get(m["gene_symbol"])
    rec = own if use_own else (uniprot(ref_acc) if ref_acc else None)
    if rec is None:
        return {"ref": None, "allowed": set(), "why": "no UniProt record with topology"}
    f = features_of(rec)
    ref_seq = rec["sequence"]["value"]
    pm = position_map(ref_seq, seq)
    ident = sum(ref_seq[r - 1] == seq[q - 1] for r, q in pm.items()) / len(seq)
    g = gpi[(gpi["gene_symbol"] == m["gene_symbol"]) & (gpi["gpi_anchored"] == 1)]
    omega_ref, omega_member, omega_basis = f["omega"], None, ""
    if len(g):
        r = g.iloc[0]
        ref_w = r["omega_site_reference"]
        if str(r["omega_form"]).endswith(acc) and pd.notna(r["omega_site"]):
            omega_member, omega_basis = int(r["omega_site"]), "gpi_set omega_form is the member"
        elif (
            m["candidate_set"] != "ortholog"
            and pd.notna(ref_w)
            and int(ref_w) <= len(seq)
            and seq[int(ref_w) - 1] == r["omega_aa"]
        ):
            # gpi_set's reference coordinate is on the library's reference form, the member
            omega_member, omega_basis = int(ref_w), "gpi_set omega_site_reference on the member"
        elif pd.notna(ref_w):
            omega_ref = int(ref_w)
    if omega_member is None and omega_ref is not None:
        omega_member = pm.get(omega_ref)
        omega_basis = "omega mapped from the reference by alignment" if omega_member else ""
    ext = f["ext"]
    if not ext:  # GPI or secreted entries carry no topological domain: the mature chain
        start = (f["chain"] or (f["signal_end"] + 1, len(ref_seq)))[0]
        ext = [(max(start, f["signal_end"] + 1), len(ref_seq))]
    if m["gene_symbol"] == "VTCN1":
        lo, hi = VTCN1_WINDOW
        ext = [(max(s, lo), min(e, hi)) for s, e, _d in f["domains"] if min(e, hi) >= max(s, lo)]
    if m["gene_symbol"] in STALK:
        ext = [
            (max(s, a), min(e, b))
            for s, e, d in f["domains"]
            for a, b in ext
            if min(e, b) >= max(s, a) and not any(u in d.lower() for u in UNFOLDED_DOMAIN)
        ]
    allowed = set()
    for s, e in ext:
        for r in range(max(s, f["signal_end"] + 1), e + 1):
            if r in pm:
                allowed.add(pm[r])
    # a cleaved propeptide takes the tag with it; an N-terminal one also moves the mature start
    for s, e in f["propeptides"]:
        allowed -= {pm[r] for r in range(s, e + 1) if r in pm}
    mature = f["signal_end"]
    for s, e in sorted(f["propeptides"]):
        if s <= mature + 1:
            mature = max(mature, e)
    if mature:
        start = pm.get(mature, mature)
        allowed = {r for r in allowed if r > start + SP_MARGIN}
    if omega_member is not None:
        allowed = {r for r in allowed if r < omega_member}
    return {
        "ref": ref_acc,
        "ref_is_identical": bool(use_own),
        "ref_identity": round(ident, 3),
        "omega_member": omega_member,
        "omega_basis": omega_basis,
        "signal_end_ref": f["signal_end"],
        "ext_ref": ";".join(f"{s}-{e}" for s, e in ext),
        "allowed": allowed,
        "why": "",
    }


def regions(res: list[int]) -> list[list[int]]:
    out: list[list[int]] = []
    for r in sorted(res):
        if out and r == out[-1][-1] + 1:
            out[-1].append(r)
        else:
            out.append([r])
    return out


def score(features: pd.DataFrame, model=V5) -> pd.DataFrame:
    """A model's scores on every row (v5 by default), exactly as tag_site_model_v4.design."""
    spec = json.loads(model.SPEC.read_text())
    feats = spec["features"]
    mz, mraw = joblib.load(model.MODEL_Z), joblib.load(model.MODEL_RAW)
    s = features.copy().reset_index(drop=True)
    for c in ("ss3", "topology_class", "region_class"):
        s[c] = s[c].fillna("").astype(str)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        s["rank_score"] = V4.smooth(s, mz.predict_proba(within_protein(s, feats, "z"))[:, 1], 2)
        s["absolute_score"] = V4.smooth(s, mraw.predict_proba(s)[:, 1], 2)
    return s


def pick(g: pd.DataFrame, allowed: set[int], vetoed) -> tuple[list[dict], str]:
    """Up to two sites inside `allowed`, v5's region-then-residue rule, contact veto applied."""
    g = g[g["residue"].isin(allowed) & g["residue"].add(1).isin(allowed)]
    if g.empty:
        if len(allowed) <= V5.MIN_REGION_LEN:
            return [], f"span too short: {len(allowed)} residue(s) after the exclusions"
        return [], (
            f"no scored junction in the {len(allowed)}-residue span: the residues DeepTMHMM v1 "
            "calls extracellular, which carry v5's features, lie outside it"
        )
    g = g[~g["residue"].map(vetoed)]
    if g.empty:
        return [], "every junction inside the span is within the ligand-contact veto"
    regs = regions(g["residue"].tolist())
    viable = [r for r in regs if len(r) >= V5.MIN_REGION_LEN]
    if not viable:
        return [], f"no region of at least {V5.MIN_REGION_LEN} residues inside the span"
    ceil = sorted(viable, key=lambda r: -g[g["residue"].isin(r)]["absolute_score"].max())
    first = g[g["residue"].isin(ceil[0])].sort_values("rank_score", ascending=False).iloc[0]
    picks = [first]
    if len(ceil) > 1:
        picks.append(
            g[g["residue"].isin(ceil[1])].sort_values("rank_score", ascending=False).iloc[0]
        )
    else:
        far = g[
            g["residue"].isin(ceil[0])
            & ((g["residue"] - first["residue"]).abs() >= V4.MIN_SEPARATION)
        ]
        if len(far):
            picks.append(far.sort_values("rank_score", ascending=False).iloc[0])
    return [p.to_dict() for p in picks], ""


def main() -> int:
    import routed_in_contacts as RC
    from l1_internal_decisions import head_tables

    m = pd.read_csv(MEMBERS, sep="\t").query("in_rescue")
    prov = pd.read_csv(PROVENANCE, sep="\t").set_index("design_id")
    feat = pd.read_csv(FEATURES, sep="\t", low_memory=False)
    hum = json.loads((UNIPROT / "_human_canonical.json").read_text())
    gpi = pd.read_csv(GPI_SET, sep="\t")
    _i, _d, prot = head_tables()
    scored = score(feat)
    tp = pd.read_csv(
        HERE.parents[1] / "sprint4" / "s4e9-orf-design-fg-panel" / "l1_tag_paths.tsv",
        sep="\t",
        usecols=["candidate_id", "shed_risk"],
    )
    shed_of = dict(zip(tp["candidate_id"], tp["shed_risk"], strict=True))
    shed = {d: shed_of.get(mid, "") for d, mid in zip(m["design_id"], m["member_id"], strict=True)}
    crow = pd.DataFrame(
        RC.contact_rows(
            [
                (r.design_id, r.gene_symbol, r.sequence_sha256, str(r.identical_acc),
                 prot[r.sequence_sha256])
                for r in m.itertuples()
            ]
        ),
        columns=RC.COLS,
    ).drop_duplicates(["sequence_sha256", "residue", "source"])  # fmt: skip
    crow.to_csv(CONTACTS, sep="\t", index=False)
    contacts = crow.groupby("sequence_sha256")["residue"].apply(set).to_dict()
    rows, summary = [], []
    for r in m.itertuples(index=False):
        key = prov.loc[r.design_id, "key"]
        seq = prot[r.sequence_sha256]
        sp = span(pd.Series(r._asdict()), seq, hum, gpi)
        g = scored[scored["uniprot_acc"] == key]
        c = contacts.get(r.sequence_sha256, set())

        def vetoed(n, c=c):
            return any(n - 2 <= x <= n + 3 for x in c)

        picks, why = pick(g, sp["allowed"], vetoed) if sp["allowed"] else ([], sp["why"])
        base = {
            "design_id": r.design_id,
            "member_id": r.member_id,
            "gene": r.gene_symbol,
            "species": r.species,
            "sequence_sha256": r.sequence_sha256,
            "placement_of_record": r.placement_of_record,
            "n_tm": r.n_tm,
            "span_reference": sp.get("ref"),
            "span_reference_identical": sp.get("ref_is_identical"),
            "span_reference_identity": sp.get("ref_identity"),
            "span_reference_extracellular": sp.get("ext_ref"),
            "omega_site_member": sp.get("omega_member"),
            "omega_basis": sp.get("omega_basis"),
            "n_allowed_residues": len(sp["allowed"]),
            "n_ligand_contacts": len(c),
            "shed_risk": shed.get(r.design_id, ""),
            "topology_note": TOPOLOGY_NOTE.get(r.gene_symbol, ""),
            "model_scope": OOD if r.n_tm <= 1 else "multipass; v5's training domain",
            "proposed_by": ATTRIBUTION,
            "implemented_by": "ranker-pu-segment session (Claude Code)",
            "proposed_on": "2026-11-14",
        }
        if not picks:
            summary.append({**base, "status": "skipped", "skip_reason": why})
            continue
        summary.append({**base, "status": "scored", "skip_reason": ""})
        for k, p in enumerate(picks, 1):
            n = int(p["residue"])
            tagged = seq[:n] + CASSETTE_AA + seq[n:]
            rows.append(
                {
                    **base,
                    "rescue_rank": k,
                    "insert_after_residue": n,
                    "junction": f"{seq[n - 1]}|{seq[n]}",
                    "rank_score": round(float(p["rank_score"]), 4),
                    "absolute_score": round(float(p["absolute_score"]), 4),
                    "plddt": p.get("plddt"),
                    "rsa": p.get("rsa"),
                    "tagged_sha256": hashlib.sha256(tagged.encode()).hexdigest(),
                }
            )
    sites = pd.DataFrame(rows)
    summ = pd.DataFrame(summary)
    sites.to_csv(OUT, sep="\t", index=False)
    summ.to_csv(HERE / "rescue_round_members_scored.tsv", sep="\t", index=False)
    counts = {
        "members_in_rescue": len(m),
        "scored": int((summ["status"] == "scored").sum()),
        "skipped": int((summ["status"] == "skipped").sum()),
        "sites_proposed": len(sites),
    }
    OUT_JSON.write_text(
        json.dumps({**counts, "gpi_set_commit": GPI_SET_COMMIT, "model_scope": OOD}, indent=2)
        + "\n"
    )
    print(json.dumps(counts))
    print(summ[summ["status"] == "skipped"][["gene", "design_id", "skip_reason"]].to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
