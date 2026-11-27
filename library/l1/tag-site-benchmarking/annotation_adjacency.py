#!/usr/bin/env python3
"""How many built internal loop arms sit next to a disulfide cysteine, glycan or intolerant site.

The agent reviewer rejects a junction on or within 2 residues of a named disulfide cysteine,
glycosylation site or insertion-intolerant position. Those are v5/v6 training hard negatives
(`ranker_pu_segment.annotate_constraints`) but not a pick-time veto, so a built arm can sit
beside one. This counts them across every built internal loop arm, to decide whether it
should become a pick-time rule.

Annotations, each moved onto the member's ordered protein by global alignment (identical
residues only, `rescue_round.position_map`):

    disulfide_cys     UniProt DISULFID (both cysteines; an interchain bond is one residue)
    glycan_n          UniProt CARBOHYD, N-linked
    sequon_n          N-X-S/T (X not P) read on the ordered protein itself
    screen_intolerant insertion-screen intolerant junctions (`dms_screen_labels.tsv`), which
                      exist only for the four screened proteins

Distance follows the contact veto: for an insert after n the flanking residues are n and n+1,
distance 0 is either of them, and "within 2" is residues n-2 ... n+3.

UniProt features come from the EBI Proteins API cache of `fetch_uniprot_site_features.py`
(data/uniprot_features/, one JSON per accession), fetched for any accession not yet cached.

    uv run --locked python annotation_adjacency.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[3] / "src"))

import l1_internal_decisions as L  # noqa: E402
from annotation_rule import distance  # noqa: E402
from fetch_uniprot_site_features import ACCESSION, fetch  # noqa: E402
from ranker_pu_segment import load_screens  # noqa: E402
from rescue_round import position_map  # noqa: E402

from deliverome_analysis.surfaceome import topology  # noqa: E402

ARMS = HERE / "l1_internal_arm_decisions.tsv"
OUT = HERE / "annotation_adjacency.tsv"
OUT_JSON = HERE / "annotation_adjacency.json"
CLASSES = ("disulfide_cys", "glycan_n", "sequon_n", "screen_intolerant")
WITHIN = 2


def annotations(acc: str, member: str, intol: dict) -> dict[str, dict[int, str]]:
    """class -> {member residue: label}."""
    out = {c: {} for c in CLASSES}
    for i in range(len(member) - 2):
        if member[i] == "N" and member[i + 1] != "P" and member[i + 2] in "ST":
            out["sequon_n"][i + 1] = f"N{i + 1} sequon"
    d = fetch(acc) if acc and ACCESSION.match(acc) else None  # designed ORFs have none
    if not d or not d.get("sequence"):
        return out
    pmap = position_map(d["sequence"], member)
    ref = d["sequence"]

    def put(cls, pos, label):
        m = pmap.get(pos)
        if m is not None and ref[pos - 1] == member[m - 1]:
            out[cls][m] = label

    for f in d.get("features", []):
        try:
            b, e = int(f["begin"]), int(f["end"])
        except KeyError, ValueError:
            continue
        if f["type"] == "DISULFID":
            kind = "interchain" if b == e else f"to C{e}" if b else ""
            put("disulfide_cys", b, f"C{b} disulfide {kind}".strip())
            if e != b:
                put("disulfide_cys", e, f"C{e} disulfide to C{b}")
        elif f["type"] == "CARBOHYD" and "N-linked" in str(f.get("description", "")):
            put("glycan_n", b, f"N{b} {f.get('description', '')}".strip())
    for r in intol.get(acc, ()):
        put("screen_intolerant", r, f"after {r} insertion-intolerant")
    return out


def main() -> int:
    arms = pd.read_csv(ARMS, sep="\t")
    built = arms[(arms["decision"] == "build") & (arms["arm_kind"] == "loop")].copy()
    _i, _d, prot = L.head_tables()
    m = built.drop_duplicates("design_id")
    t = topology.resolve({d: prot[s] for d, s in zip(m["design_id"], m["sequence_sha256"],
                                                      strict=True)}).set_index("key")  # fmt: skip
    sc = load_screens()
    sc = sc[(sc["kind"] == "insertion") & (sc["tolerance_class"] == "intolerant")]
    intol = sc.groupby("uniprot_acc")["residue"].apply(lambda s: {int(x) for x in s}).to_dict()
    ann = {}
    for d, s in zip(m["design_id"], m["sequence_sha256"], strict=True):
        acc = t.loc[d, "uniprot_acc"] if bool(t.loc[d, "found"]) else None
        ann[d] = (acc, annotations(acc, prot[s], intol))
    rows = []
    for x in built.itertuples(index=False):
        n = int(x.insert_after_residue)
        acc, a = ann[x.design_id]
        row = {"arm_id": x.arm_id, "design_id": x.design_id, "gene": x.gene,
               "uniprot_acc": acc, "insert_after_residue": n,
               "member_built_loop_arms": int((built["design_id"] == x.design_id).sum()),
               "member_built_arms": int(((arms["design_id"] == x.design_id)
                                         & (arms["decision"] == "build")).sum())}  # fmt: skip
        hit = []
        for c in CLASSES:
            near = sorted(((distance(r, n), r, lab) for r, lab in a[c].items()))
            dmin = near[0][0] if near else None
            row[f"{c}_nearest_distance"] = dmin
            row[f"{c}_nearest"] = near[0][2] if near else ""
            row[f"{c}_within_{WITHIN}"] = dmin is not None and dmin <= WITHIN
            if dmin is not None and dmin <= WITHIN:
                hit.append(f"{c}:{near[0][2]} (distance {dmin})")
        row["any_within_2"] = bool(hit)
        row["hits"] = "; ".join(hit)
        rows.append(row)
    out = pd.DataFrame(rows)
    out.to_csv(OUT, sep="\t", index=False)
    w = out[out["any_within_2"]]
    uniprot_cls = ["disulfide_cys", "glycan_n", "screen_intolerant"]
    annotated = out[[f"{c}_within_{WITHIN}" for c in uniprot_cls]].any(axis=1)
    summary = {
        "built_loop_arms": len(out),
        "members": int(out["design_id"].nunique()),
        "members_without_uniprot_entry": int(
            out.drop_duplicates("design_id")["uniprot_acc"].isna().sum()
        ),
        "window": "residues n-2 ... n+3 for an insert after n (as the contact veto)",
        "arms_within_2_by_class": {c: int(out[f"{c}_within_{WITHIN}"].sum()) for c in CLASSES},
        "arms_within_2_of_annotated_disulfide_glycan_or_intolerant": int(annotated.sum()),
        "arms_within_2_of_any_including_unannotated_sequons": len(w),
        "of_which_the_members_only_loop_arm": int(
            (annotated & (out["member_built_loop_arms"] == 1)).sum()
        ),
        "of_which_the_members_only_arm": int((annotated & (out["member_built_arms"] == 1)).sum()),
        "members_with_an_annotated_hit": int(out.loc[annotated, "design_id"].nunique()),
        "on_the_flanking_residue_itself": {
            c: int((out[f"{c}_nearest_distance"] == 0).sum()) for c in CLASSES
        },
    }
    OUT_JSON.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
