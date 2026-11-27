#!/usr/bin/env python3
"""Built loop arms in re-entrant pore loops, functional sites and channel signature motifs.

The agent review found KCNC3 `l1-p0997-a1` inside the GYGD selectivity filter: a tag in a
potassium channel's re-entrant P-loop that passed the membrane-side gate, the contact veto and
refold. Re-entrant loops (K+ channel P-loops, glutamate-receptor and P2X pore loops, aquaporin
NPA half-helices, CLC and transporter hairpins) are not extracellular loops. Coordinator,
2026-11-15: count every built loop arm, plus the pending revisions, against

    (a) intramem       inside or within 2 residues of a UniProt INTRAMEM span
    (b) functional     on or within 2 residues of an ACT_SITE, BINDING or SITE feature
    (c) pore_motif     inside or within 2 residues of a MOTIF or REGION whose description names a
                       selectivity filter, pore, gate or signature motif
    (d) k_filter       potassium channels without such annotation: inside or within 2 residues of
                       a T-x-G-[YF]-G filter run or the 12-residue pore helix N-terminal to it,
                       read from the sequence

For each hit, what the r8 mask (`mask`, `junction_mask` at n) and the arm's gate column recorded.
Features come from the EBI Proteins API (one JSON per accession, cached in
data/uniprot_features_pore/) and are mapped onto each ordered protein by global alignment,
keeping identical residues only. "Within 2" uses the contact veto's window, n-2 ... n+3 for an
insert after n: a span [b, e] hits when it overlaps that window.

    uv run --locked python pore_rule.py   # writes pore_rule_arms.tsv and pore_rule.json
"""

from __future__ import annotations

import json
import re
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[3] / "src"))

TYPES = ["INTRAMEM", "MOTIF", "REGION", "ACT_SITE", "BINDING", "SITE"]
API = "https://www.ebi.ac.uk/proteins/api/features/{acc}?types=" + ",".join(TYPES)
CACHE = HERE / "data" / "uniprot_features_pore"
PORE_WORDS = re.compile(
    r"selectivity filter|pore|gate|gating|signature|sequence motif", re.IGNORECASE
)
K_FILTER = re.compile(r"T.G[YF]G")
PORE_HELIX = 12
WINDOW = 2
ARMS = HERE / "l1_internal_arm_decisions.tsv"
PENDING = HERE / "l1_internal_arm_decisions_pending.tsv"
MASK_R8 = HERE / "data" / "mask_r8" / "l1_extracellular_mask_residues.tsv"
OUT = HERE / "pore_rule_arms.tsv"
OUT_JSON = HERE / "pore_rule.json"
CLASSES = ("intramem", "functional", "pore_motif", "k_filter")


def fetch(acc: str) -> dict | None:
    hit = CACHE / f"{acc}.json"
    if hit.is_file():
        return json.loads(hit.read_text()) or None
    req = urllib.request.Request(
        API.format(acc=acc),
        headers={"Accept": "application/json", "User-Agent": "deliverome-analysis s5e13"},
    )
    d = None
    for attempt in range(4):
        try:
            with urllib.request.urlopen(req, timeout=60) as r:
                d = json.loads(r.read())
            break
        except urllib.error.HTTPError as e:
            if e.code == 404:
                break
            time.sleep(2 + 3 * attempt)
        except urllib.error.URLError, TimeoutError:
            time.sleep(2 + 3 * attempt)
    else:
        raise RuntimeError(f"EBI Proteins API failed for {acc}")
    CACHE.mkdir(parents=True, exist_ok=True)
    hit.write_text(json.dumps(d or {}))
    time.sleep(0.3)
    return d


def overlaps(b: int, e: int, n: int) -> int | None:
    """Distance from span [b, e] to the junction after n, or None beyond the window."""
    lo, hi = n - WINDOW, n + 1 + WINDOW
    if e < lo or b > hi:
        return None
    if b <= n + 1 and e >= n:
        return 0
    return (n - e) if e < n else (b - (n + 1))


def human_accessions() -> dict[str, str]:
    """gene -> human canonical accession (D1 human_canonical cohort, as contacts_by_protein)."""
    import d1_snapshot

    out = {}
    for r in d1_snapshot.query(
        "SELECT gene_symbol, uniprot_acc FROM topology_public WHERE cohort = 'human_canonical' "
        "ORDER BY topology_version"
    ):
        out.setdefault(r["gene_symbol"], r["uniprot_acc"])
    return out


def spans(acc: str, member: str) -> dict[str, list[tuple[int, int, str]]]:
    """class -> [(begin, end, label)] on the member's ordered protein.

    `acc` is the member's own accession, or for a designed ORF or a non-human ortholog the
    gene's human canonical accession; either way features move only onto identical residues.
    """
    from fetch_uniprot_site_features import ACCESSION
    from rescue_round import position_map

    out = {c: [] for c in CLASSES}
    d = fetch(acc) if acc and ACCESSION.match(acc) else None
    annotated_pore = False
    if d and d.get("sequence"):
        ref = d["sequence"]
        pmap = position_map(ref, member)

        def mp(p):
            m = pmap.get(p)
            return m if m is not None and ref[p - 1] == member[m - 1] else None

        for f in d.get("features", []):
            try:
                b, e = int(f["begin"]), int(f["end"])
            except KeyError, ValueError:
                continue
            desc = str(f.get("description", ""))
            mb, me = mp(b), mp(e)
            if mb is None or me is None:
                continue
            lab = f"{f['type']} {b}-{e} {desc}".strip()
            if f["type"] == "INTRAMEM":
                out["intramem"].append((mb, me, lab))
            elif f["type"] in ("ACT_SITE", "BINDING", "SITE"):
                out["functional"].append((mb, me, lab))
            elif f["type"] in ("MOTIF", "REGION") and PORE_WORDS.search(desc):
                out["pore_motif"].append((mb, me, lab))
                annotated_pore = True
    return out, annotated_pore


def main() -> int:
    import l1_internal_decisions as L

    from deliverome_analysis.surfaceome import topology

    tabs = [pd.read_csv(ARMS, sep="\t").assign(table="canonical")]
    if PENDING.is_file():
        tabs.append(pd.read_csv(PENDING, sep="\t").assign(table="pending"))
    arms = pd.concat(tabs)
    built = arms[(arms["decision"] == "build") & (arms["arm_kind"] == "loop")]
    built = built.drop_duplicates(["design_id", "insert_after_residue"])
    _i, _d, prot = L.head_tables()
    m = built.drop_duplicates("design_id")
    seqs = dict(zip(m["design_id"], m["sequence_sha256"].map(prot), strict=True))
    t = topology.resolve(seqs).set_index("key")
    r8 = pd.read_csv(MASK_R8, sep="\t").drop_duplicates("sequence_sha256")
    r8 = r8.set_index("sequence_sha256")
    from fetch_uniprot_site_features import ACCESSION

    human = human_accessions()
    gene_of = dict(zip(built["design_id"], built["gene"], strict=False))
    ann, basis = {}, {}
    for did in seqs:
        acc = t.loc[did, "uniprot_acc"] if bool(t.loc[did, "found"]) else None
        basis[did] = "own accession"
        if not (acc and ACCESSION.match(str(acc))):
            acc = human.get(gene_of[did])
            basis[did] = f"via human canonical {acc}" if acc else "no accession"
        sp, has_pore = spans(acc, seqs[did])
        # a non-human member's own entry is often unannotated: add the human canonical's
        # features, moved onto the member's identical residues
        hum = human.get(gene_of[did])
        species = str(t.loc[did, "species"]) if bool(t.loc[did, "found"]) else ""
        if hum and hum != acc and species and species != "human":
            extra, extra_pore = spans(hum, seqs[did])
            for c in CLASSES:
                sp[c] += extra[c]
            has_pore = has_pore or extra_pore
            basis[did] += f" + human canonical {hum}"
        gene = str(t.loc[did, "gene_symbol"]) if bool(t.loc[did, "found"]) else ""
        if gene.startswith("KCN") and not has_pore:
            for mo in K_FILTER.finditer(seqs[did]):
                b, e = mo.start() + 1, mo.end()
                sp["k_filter"].append((max(1, b - PORE_HELIX), e,
                                       f"sequence {mo.group()} at {b}-{e} + pore helix"))  # fmt: skip
        ann[did] = sp
    rows = []
    for x in built.itertuples(index=False):
        n = int(x.insert_after_residue)
        sha = x.sequence_sha256
        has = sha in r8.index and int(r8.loc[sha, "usable"]) == 1
        row = {"arm_id": x.arm_id, "table": x.table, "design_id": x.design_id, "gene": x.gene,
               "feature_basis": basis[x.design_id],
               "insert_after_residue": n,
               "member_built_arms": int(((arms["table"] == x.table) & (arms["design_id"]
                                         == x.design_id) & (arms["decision"] == "build")).sum()),
               "mask_at_n": r8.loc[sha, "mask"][n - 1] if has else "",
               "mask_at_n1": r8.loc[sha, "mask"][n] if has and n < len(r8.loc[sha, "mask"]) else "",
               "junction_mask_at_n": r8.loc[sha, "junction_mask"][n - 1] if has else "",
               "r8_placement_source": r8.loc[sha, "placement_source"] if has else "no r8 row",
               "gate": x.gate, "refold_verdict": x.refold_verdict}  # fmt: skip
        hits = []
        for c in CLASSES:
            hs = [(overlaps(b, e, n), lab) for b, e, lab in ann[x.design_id][c]]
            hs = sorted((d, lab) for d, lab in hs if d is not None)
            row[c] = "; ".join(f"{lab} (distance {d})" for d, lab in hs)
            if hs:
                hits.append(c)
        row["classes"] = ";".join(hits)
        row["case"] = (
            "replace" if set(hits) & {"intramem", "functional", "pore_motif"}
            else "flag" if hits else "none"
        )  # fmt: skip
        rows.append(row)
    out = pd.DataFrame(rows)
    out.to_csv(OUT, sep="\t", index=False)
    canon = out[out["table"] == "canonical"]
    summary = {
        "built_loop_arms_canonical": len(canon),
        "pending_only_arms": int((out["table"] == "pending").sum()),
        "hits_by_class_canonical": {c: int((canon[c] != "").sum()) for c in CLASSES},
        "hits_by_class_pending_only": {
            c: int((out.loc[out["table"] == "pending", c] != "").sum()) for c in CLASSES
        },
        "feature_basis": out.drop_duplicates("design_id")["feature_basis"]
        .str.split(" ")
        .str[0]
        .value_counts()
        .to_dict(),
        "replace_case_arms": int((out["case"] == "replace").sum()),
        "replace_case_only_arm": int(
            ((out["case"] == "replace") & (out["member_built_arms"] == 1)).sum()
        ),
        "flag_case_arms": int((out["case"] == "flag").sum()),
        "hit_arms_junction_mask_E": int(
            ((out["case"] != "none") & (out["junction_mask_at_n"] == "E")).sum()
        ),
        "hit_arms_gate": out.loc[out["case"] != "none", "gate"].value_counts().to_dict(),
    }
    OUT_JSON.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
