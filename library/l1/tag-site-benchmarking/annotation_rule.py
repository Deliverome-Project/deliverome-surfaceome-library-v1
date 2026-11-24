#!/usr/bin/env python3
"""Pick-time rule for junctions beside annotated disulfide cysteines, glycans and intolerant sites.

Coordinator under Becca's delegation, 2026-11-15, on the count in `annotation_adjacency.py`
(80 of 1,065 built loop arms within 2 residues of such a residue; 31 of them a member's only
arm). Evidence beside the ruling: 12 of the 111 published working sites lie within 2 residues of
a v5 training hard negative, nearly all beside a sequon asparagine, so proximity to a glycan is
tolerated in practice; a disulfide cysteine or an insertion-intolerant position is not.

For a cassette inserted after residue n the window is residues n-2 ... n+3, as the contact veto.

| Case | Annotation in the window | Rule |
|---|---|---|
| replace | UniProt disulfide cysteine | replace (below) |
| replace | insertion-intolerant junction from an insertion screen | replace |
| replace | the junction splits an annotated N-glycan sequon (insert after N or after X of N-X-S/T) | replace: the glycan is removed by construction |
| flag | an annotated N-glycan beside the junction, sequon intact | flag only |
| flag | the junction splits a predicted (unannotated) sequon whose N is extracellular | flag "splits a predicted sequon"; no replacement on prediction alone |
| flag | a predicted sequon beside the junction | flag only |

Replace: the best clean site at least MIN_GAP residues from the offending residue, in the same
region within 10% of the arm's v6 rank score if one exists, otherwise the best clean site in any
eligible region. Clean means: eligible under mask r8, no hard or soft contact in the window, no
replace-case annotation, and not within MIN_SEPARATION of another built arm in the same region.
The replacement is refolded; if it fails and the member has another arm the arm is dropped; if it
is the member's only arm it is kept, flagged, and listed for Becca.

    uv run --locked python annotation_rule.py   # writes annotation_rule_arms.tsv and the fold batch
"""

from __future__ import annotations

import hashlib
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[3] / "src"))


def distance(r: int, n: int) -> int:
    """Residues from residue r to the junction after n (0 on either flanking residue)."""
    return 0 if r in (n, n + 1) else (n - r if r < n else r - (n + 1))


WINDOW = 2
MIN_GAP = 5
MIN_SEPARATION = 25  # the order rule's distance between two arms in one region
WITHIN_SCORE = 0.9
REPLACE_CLASSES = ("disulfide_cys", "screen_intolerant")
#: feature spans from pore_rule.member_spans (coordinator, 2026-11-15, after KCNC3's GYGD arm):
#: INTRAMEM, ACT_SITE/BINDING/SITE and pore/filter/gate/signature MOTIF or REGION are replaced;
#: a sequence-only K+ filter run (no annotation) is flagged for review
SPAN_REPLACE = ("intramem", "functional", "pore_motif")
SPAN_FLAG = ("k_filter",)
NEAR = 5
#: the per-arm table the rule is applied to: #370 as merged. Pinned, because l1_internal_arms.py
#: applies this script's output to that same table, and a rerun on the revised table would
#: no longer see the arms it moved.
ARMS_COMMIT = "798d8717"
ARMS_REL = "library/l1/tag-site-benchmarking/l1_internal_arm_decisions.tsv"
MASK_R8 = HERE / "data" / "mask_r8" / "l1_extracellular_mask_residues.tsv"
OUT = HERE / "annotation_rule_arms.tsv"
OUT_FOLD = HERE / "annotation_rule_to_fold.tsv"
OUT_FASTA = HERE / "annotation_rule_to_fold.fasta"
OUT_JSON = HERE / "annotation_rule.json"
USD_PER_RESIDUE = (0.000062, 0.00008)
#: name of the next replacement batch (the first, annotation_rule_to_fold.tsv, was submitted
#: on 2026-11-15)
BATCH = "2026-11-15b"


def sequon_at(seq: str, p: int) -> bool:
    """N-X-S/T (X not P) starting at residue p (1-based)."""
    return p + 2 <= len(seq) and seq[p - 1] == "N" and seq[p] != "P" and seq[p + 1] in "ST"


def classify(ann: dict, seq: str, n: int, mask: str | None = None) -> list[dict]:
    """Every annotation the rule acts on for an insert after n, with its case."""
    out = []
    for cls in REPLACE_CLASSES:
        for r, lab in ann[cls].items():
            d = distance(r, n)
            if d <= WINDOW:
                out.append({"case": "replace", "class": cls, "residue": r, "distance": d,
                            "label": lab})  # fmt: skip
    for r, lab in ann["glycan_n"].items():
        d = distance(r, n)
        if n in (r, r + 1) and sequon_at(seq, r):
            out.append({"case": "replace", "class": "splits_annotated_sequon", "residue": r,
                        "distance": d, "label": lab})  # fmt: skip
        elif d <= WINDOW:
            out.append({"case": "flag", "class": "glycan_n", "residue": r, "distance": d,
                        "label": lab})  # fmt: skip
    for r, lab in ann["sequon_n"].items():
        if r in ann["glycan_n"]:
            continue
        d = distance(r, n)
        outside = mask is None or (1 <= r <= len(mask) and mask[r - 1] == "E")
        if n in (r, r + 1) and outside:
            out.append({"case": "flag", "class": "splits_predicted_sequon", "residue": r,
                        "distance": d, "label": lab})  # fmt: skip
        elif d <= WINDOW:
            out.append({"case": "flag", "class": "sequon_n", "residue": r, "distance": d,
                        "label": lab})  # fmt: skip
    return out + classify_spans(ann, n)


def member_annotations(acc, gene: str, species: str, seq: str, intol: dict, human: dict) -> dict:
    """Residue annotations and feature spans for one ordered protein.

    Features come from the member's own UniProt entry; a designed ORF (no entry) or a
    non-human member also takes the gene's human canonical entry. Every feature moves onto the
    ordered protein by global alignment, identical residues only.
    """
    from annotation_adjacency import annotations
    from fetch_uniprot_site_features import ACCESSION
    from pore_rule import K_FILTER, PORE_HELIX, spans

    accs = [acc] if acc and ACCESSION.match(str(acc)) else []
    hum = human.get(gene)
    if hum and hum not in accs and (not accs or species != "human"):
        accs.append(hum)
    ann = annotations(None, seq, intol)  # sequons only
    sp_all = {c: [] for c in (*SPAN_REPLACE, *SPAN_FLAG)}
    has_pore = False
    for a in accs:
        one = annotations(a, seq, intol)
        for k, v in one.items():
            ann[k] = {**ann.get(k, {}), **v}
        sp, pore = spans(a, seq)
        for c in sp_all:
            sp_all[c] += sp.get(c, [])
        has_pore = has_pore or pore
    if str(gene).startswith("KCN") and not has_pore:
        for mo in K_FILTER.finditer(seq):
            b, e = mo.start() + 1, mo.end()
            sp_all["k_filter"].append((max(1, b - PORE_HELIX), e,
                                       f"sequence {mo.group()} at {b}-{e} + pore helix"))  # fmt: skip
    ann["spans"] = sp_all
    ann["sources"] = accs
    return ann


def span_gap(b: int, e: int, k: int) -> int:
    """Residues between the span [b, e] and the junction after k (0 when they touch)."""
    if b <= k + 1 and e >= k:
        return 0
    return (k - e) if e < k else (b - (k + 1))


def residue_gap(b: int, e: int, k: int) -> int:
    """|k - residue| to the nearest residue of [b, e] (0 inside): the MIN_GAP measure, as the
    first rule used it for single residues (ABCG2 after 598 is 5 from C603)."""
    return 0 if b <= k <= e else min(abs(k - b), abs(k - e))


PORE_ADJACENT = 10
PORE_SPANS = ("intramem", "pore_motif")


def pore_adjacent(ann: dict, k: int) -> str:
    """'pore-adjacent: N residues from <feature>' within PORE_ADJACENT of an INTRAMEM or
    named pore/filter span, else ''. The residues between a pore helix and the next helix are
    covered where UniProt annotates the pore helix as INTRAMEM."""
    best = None
    for cls in PORE_SPANS:
        for b, e, lab in ann.get("spans", {}).get(cls, []):
            if cls != "intramem" and not re.search(r"pore|filter", lab, re.IGNORECASE):
                continue  # e.g. SCNN1A's GRIP "gating ... by proteolysis" region is not a pore
            g = residue_gap(b, e, k)
            if g <= PORE_ADJACENT and (best is None or g < best[0]):
                best = (g, lab)
    return f"pore-adjacent: {best[0]} residues from {best[1]}" if best else ""


def classify_spans(ann: dict, n: int) -> list[dict]:
    out = []
    for case, classes in (("replace", SPAN_REPLACE), ("flag", SPAN_FLAG)):
        for cls in classes:
            for b, e, lab in ann.get("spans", {}).get(cls, []):
                d = span_gap(b, e, n)
                if d <= WINDOW:
                    out.append({"case": case, "class": cls, "residue": b, "end": e,
                                "distance": d, "label": lab})  # fmt: skip
                elif case == "replace" and d <= NEAR:
                    # 3 to 5 residues from a functional or pore feature: flag only
                    out.append({"case": "flag", "class": f"{cls}_near", "residue": b, "end": e,
                                "distance": d, "label": lab})  # fmt: skip
    return out


def blocked(ann: dict, seq: str, n: int, mask: str | None = None) -> bool:
    """True where the rule forbids a pick at n (any replace-case annotation)."""
    return any(v["case"] == "replace" for v in classify(ann, seq, n, mask))


def pinned_arms() -> pd.DataFrame:
    import io

    from deliverome_analysis.provenance import read_pinned_text

    txt = read_pinned_text(ARMS_COMMIT, ARMS_REL, repo=HERE.parents[2])
    return pd.read_csv(io.StringIO(txt), sep="\t")


class Checker:
    """Pick-time test for any junction of an internal member: `checker.blocked(sha, n)`.

    Annotations are built on first use per member (UniProt features from the EBI cache,
    mapped onto the ordered protein; insertion-screen labels; sequons read on the protein).
    """

    def __init__(self):
        import l1_internal_decisions as L
        from ranker_pu_segment import load_screens

        _i, _d, self.prot = L.head_tables()
        r8 = pd.read_csv(MASK_R8, sep="\t").drop_duplicates("sequence_sha256")
        self.r8 = r8.set_index("sequence_sha256")
        sc = load_screens()
        sc = sc[(sc["kind"] == "insertion") & (sc["tolerance_class"] == "intolerant")]
        self.intol = (
            sc.groupby("uniprot_acc")["residue"].apply(lambda s: {int(x) for x in s}).to_dict()
        )
        self.ann: dict[str, dict] = {}
        from pore_rule import human_accessions

        self.human = human_accessions()
        # the library's gene for each ordered protein (an ortholog's v1 row carries its entry
        # name, e.g. A0A2K5W8Z7_MACFA, where the gene symbol should be)
        dec = pd.read_csv(L.OUT, sep="\t")
        self.gene = dict(zip(dec["sequence_sha256"], dec["gene"], strict=False))

    def _ann(self, sha: str) -> dict:
        if sha not in self.ann:
            from deliverome_analysis import topology

            t = topology.resolve({"m": self.prot[sha]}).iloc[0]
            found = bool(t["found"])
            self.ann[sha] = member_annotations(
                t["uniprot_acc"] if found else None,
                self.gene.get(sha) or (str(t["gene_symbol"]) if found else ""),
                str(t["species"]) if found else "",
                self.prot[sha], self.intol, self.human,
            )  # fmt: skip
        return self.ann[sha]

    def _mask(self, sha: str) -> str | None:
        ok = sha in self.r8.index and int(self.r8.loc[sha, "usable"]) == 1
        return self.r8.loc[sha, "mask"] if ok else None

    def classify(self, sha: str, n: int) -> list[dict]:
        return classify(self._ann(sha), self.prot[sha], int(n), self._mask(sha))

    def blocked(self, sha: str, n: int) -> bool:
        return any(v["case"] == "replace" for v in self.classify(sha, n))


def describe(vs: list[dict]) -> str:
    return "; ".join(f"{v['class']}: {v['label']} (distance {v['distance']})" for v in vs)


def main() -> int:
    import l1_internal_decisions as L
    import rescue_round as RR
    import tag_site_model_v6 as V6
    from arm_revisions import CASSETTE_AA, contact_tiers, hits
    from l1_region_tables import frames, membrane_mask
    from ranker_pu_segment import load_screens
    from rescue_round import score
    from v5_region_table import features

    from deliverome_analysis import topology

    arms = pinned_arms()
    built = arms[(arms["decision"] == "build") & (arms["arm_kind"] == "loop")].copy()
    dec = pd.read_csv(L.OUT, sep="\t")
    _i, _d, prot = L.head_tables()
    r8 = pd.read_csv(MASK_R8, sep="\t").drop_duplicates("sequence_sha256")
    r8 = r8.set_index("sequence_sha256")
    m = built.drop_duplicates("design_id")
    seqs = dict(zip(m["design_id"], m["sequence_sha256"].map(prot), strict=True))
    t = topology.resolve(seqs).set_index("key")
    sc = load_screens()
    sc = sc[(sc["kind"] == "insertion") & (sc["tolerance_class"] == "intolerant")]
    intol = sc.groupby("uniprot_acc")["residue"].apply(lambda s: {int(x) for x in s}).to_dict()

    def mask_of(sha):
        ok = sha in r8.index and int(r8.loc[sha, "usable"]) == 1
        return (r8.loc[sha, "mask"], r8.loc[sha, "junction_mask"]) if ok else None

    from pore_rule import human_accessions

    human = human_accessions()
    ann, cls_rows = {}, []
    for x in built.itertuples(index=False):
        if x.design_id not in ann:
            tr = t.loc[x.design_id]
            found = bool(tr["found"])
            ann[x.design_id] = member_annotations(
                tr["uniprot_acc"] if found else None,
                x.gene,  # the library's gene (an ortholog's v1 row names its entry, not a gene)
                str(tr["species"]) if found else "",
                seqs[x.design_id], intol, human,
            )  # fmt: skip
        mk = mask_of(x.sequence_sha256)
        vs = classify(ann[x.design_id], seqs[x.design_id], int(x.insert_after_residue),
                      mk[0] if mk else None)  # fmt: skip
        cls_rows.append((x, vs))

    to_fix = [(x, vs) for x, vs in cls_rows if any(v["case"] == "replace" for v in vs)]
    # v6 scores on every eligible junction of the members that need a replacement
    fr = frames(dec)
    shas = {x.sequence_sha256 for x, _ in to_fix}
    want = {fr[s]: s for s in shas if s in fr}
    allx = features()
    pv = pd.read_csv(RR.PROVENANCE, sep="\t")
    rk = {k: s for k, s in zip(pv["key"], pv["sequence_sha256"], strict=True) if s in shas}
    rf = pd.read_csv(RR.FEATURES, sep="\t", low_memory=False)
    want = {k: s for k, s in want.items() if s not in set(rk.values())} | rk
    allx = pd.concat([allx[~allx["uniprot_acc"].isin(rk)], rf[rf["uniprot_acc"].isin(rk)]])
    s6 = score(allx[allx["uniprot_acc"].isin(want)], V6)
    s6["sha"] = s6["uniprot_acc"].map(want)
    cand = {}
    for sha, g in s6.groupby("sha"):
        mk = mask_of(sha)
        if mk:
            g = g[membrane_mask.eligible_after(mk, g["residue"])].copy()
            g["region"] = membrane_mask.masked_regions(g["residue"], mk)
        else:
            g = g.assign(region=g["region_index"])
        cand[sha] = g.sort_values("rank_score", ascending=False)
    ct = contact_tiers(pd.read_csv(HERE / "contacts_by_protein.tsv", sep="\t", low_memory=False))
    ct = ct[ct["tier"].isin(["hard", "soft"])]
    # every verdict already on record (s4e9 tables, then the s5e13 verdicts seen)
    vt = pd.read_csv(HERE.parents[1] / "sprint4" / "s4e9-orf-design-fg-panel"
                     / "tag_refold_verdicts.tsv", sep="\t", low_memory=False)  # fmt: skip
    vs_seen = pd.read_csv(L.VERDICTS_SEEN, sep="\t")
    verdict = dict(zip(vs_seen["tagged_sha256"], vs_seen["verdict"], strict=False))
    verdict |= dict(zip(vt["tagged_sha256"], vt["verdict"], strict=False))
    seen = {k for k, x in verdict.items() if x in ("pass", "fail")}

    rows = []
    for x, vs in cls_rows:
        n = int(x.insert_after_residue)
        rep = [v for v in vs if v["case"] == "replace"]
        row = {"arm_id": x.arm_id, "design_id": x.design_id, "gene": x.gene,
               "sequence_sha256": x.sequence_sha256, "insert_after_residue": n,
               "member_built_arms": int(((arms["design_id"] == x.design_id)
                                         & (arms["decision"] == "build")).sum()),
               "case": "replace" if rep else ("flag" if vs else "none"),
               "classes": ";".join(sorted({v["class"] for v in vs})),
               "annotations": describe(vs), "v6_rank_score": x.v6_rank_score}  # fmt: skip
        if rep:
            g = cand.get(x.sequence_sha256)
            seq, mk = seqs[x.design_id], mask_of(x.sequence_sha256)
            bad = [(v["residue"], v.get("end", v["residue"])) for v in rep]
            others = built[(built["design_id"] == x.design_id) & (built["arm_id"] != x.arm_id)]
            region = None
            if g is not None and n in set(g["residue"]):
                region = int(g.loc[g["residue"] == n, "region"].iloc[0])
            pick, basis = None, "no clean eligible site"
            if g is not None:
                ok = []
                for c in g.itertuples(index=False):
                    k = int(c.residue)
                    if k == n or min(residue_gap(b, e, k) for b, e in bad) < MIN_GAP:
                        continue
                    if blocked(ann[x.design_id], seq, k, mk[0] if mk else None):
                        continue
                    if not hits(ct, x.sequence_sha256, k).empty:
                        continue
                    t_k = hashlib.sha256((seq[:k] + CASSETTE_AA + seq[k:]).encode()).hexdigest()
                    if verdict.get(t_k) == "fail":  # already folded and failed
                        continue
                    near = [
                        int(o) for o in others["insert_after_residue"]
                        if abs(int(o) - k) < MIN_SEPARATION
                        and int(g.loc[g["residue"] == int(o), "region"].iloc[0]
                                if (g["residue"] == int(o)).any() else -1) == int(c.region)
                    ]  # fmt: skip
                    if near:
                        continue
                    ok.append(c)
                same = [c for c in ok if region is not None and int(c.region) == region]
                own = float(x.v6_rank_score) if x.v6_rank_score == x.v6_rank_score else np.nan
                within = [c for c in same if own == own and c.rank_score >= WITHIN_SCORE * own]

                def unflagged(cs, mk=mk, seq=seq, did=x.design_id):
                    # within a tier, a site with no flag-case annotation (a glycan or predicted
                    # sequon beside it) beats a higher-scoring one that has one
                    free = [
                        c
                        for c in cs
                        if not classify(ann[did], seq, int(c.residue), mk[0] if mk else None)
                    ]
                    return (free or cs)[0]  # fmt: skip

                if within:
                    pick, basis = unflagged(within), "same region, within 10% of the arm (v6)"
                elif ok:
                    pick, basis = unflagged(ok), "best clean site in any eligible region (v6)"
                # a pick within PORE_ADJACENT of a pore span is still pore-adjacent: prefer a
                # clean site away from the pore scoring within 10% of it (coordinator, 2026-11-15)
                if pick is not None and pore_adjacent(ann[x.design_id], int(pick.residue)):
                    away = [c for c in ok if not pore_adjacent(ann[x.design_id], int(c.residue))
                            and int(c.region) != int(pick.region)
                            and c.rank_score >= WITHIN_SCORE * pick.rank_score]  # fmt: skip
                    if away:
                        pick = unflagged(away)
                        basis = ("away from the pore: another region, within 10% of the best "
                                 "pore-adjacent site (v6)")  # fmt: skip
            if pick is not None:
                k = int(pick.residue)
                tagged = hashlib.sha256((seq[:k] + CASSETTE_AA + seq[k:]).encode()).hexdigest()
                row |= {"replacement": k, "replacement_region": int(pick.region),
                        "arm_region": region,
                        "replacement_v6_rank_score": round(float(pick.rank_score), 4),
                        "replacement_basis": basis, "replacement_tagged_sha256": tagged,
                        "replacement_flag": pore_adjacent(ann[x.design_id], k),
                        "replacement_already_folded": tagged in seen,
                        "replacement_verdict_on_record": verdict.get(tagged, ""),
                        "length_aa": len(seq) + len(CASSETTE_AA),
                        "tagged_protein": seq[:k] + CASSETTE_AA + seq[k:]}  # fmt: skip
            else:
                row |= {"replacement": None, "arm_region": region, "replacement_basis": basis}
        rows.append(row)
    out = pd.DataFrame(rows)
    out.drop(columns=["tagged_protein"], errors="ignore").to_csv(OUT, sep="\t", index=False)
    rp = out[out["case"] == "replace"]
    fold = rp[rp["replacement"].notna() & ~rp["replacement_already_folded"].fillna(False)
              .astype(bool)].drop_duplicates("replacement_tagged_sha256")  # fmt: skip
    # a submitted batch is a record: never rewritten; constructs not yet submitted go to the
    # next dated batch
    batches = sorted(HERE.glob("annotation_rule_to_fold*.tsv"))
    sent = set()
    for f in batches:
        sent |= set(pd.read_csv(f, sep="\t")["replacement_tagged_sha256"])
    fold = fold[~fold["replacement_tagged_sha256"].isin(sent)]
    out_fold, out_fasta = OUT_FOLD, OUT_FASTA
    if batches:
        out_fold = HERE / f"annotation_rule_to_fold_{BATCH}.tsv"
        out_fasta = HERE / f"annotation_rule_to_fold_{BATCH}.fasta"
        if out_fold.is_file():
            fold = fold.iloc[0:0]
    if len(fold):
        fold[["design_id", "gene", "arm_id", "sequence_sha256", "insert_after_residue", "replacement",
              "replacement_tagged_sha256", "length_aa", "tagged_protein"]].to_csv(
            out_fold, sep="\t", index=False)  # fmt: skip
    if len(fold):
        with out_fasta.open("w") as fh:
            for f in fold.itertuples(index=False):
                fh.write(f">{f.design_id}:after{int(f.replacement)} gene={f.gene} "
                       f"request=annotation_rule_replace replaces={f.arm_id} "
                       f"tagged_sha256={f.replacement_tagged_sha256} "
                       f"native_sha256={f.sequence_sha256}\n{f.tagged_protein}\n")  # fmt: skip
    flags = out[out["case"] == "flag"]
    aa = int(fold["length_aa"].sum()) if len(fold) else 0

    def has(df, c):
        return int(df["classes"].str.split(";").apply(lambda s: c in s).sum())

    glyc = out["classes"].str.contains("glycan_n|splits_annotated_sequon", regex=True)
    summary = {
        "built_loop_arms": len(out),
        "replace_case_arms": len(rp),
        "replace_by_class": {c: has(rp, c) for c in
                             ("disulfide_cys", "screen_intolerant", "splits_annotated_sequon",
                              *SPAN_REPLACE)},
        "flag_by_class": {c: has(out[out["case"] == "flag"], c) for c in
                          ("glycan_n", "splits_predicted_sequon", "sequon_n", *SPAN_FLAG,
                           *(f"{c}_near" for c in SPAN_REPLACE))},
        "replace_members_only_arm": int((rp["member_built_arms"] == 1).sum()),
        "annotated_glycan_arms": int(glyc.sum()),
        "annotated_glycan_split_sequon": has(out, "splits_annotated_sequon"),
        "annotated_glycan_beside_only": int((glyc & ~out["classes"].str.contains(
            "splits_annotated_sequon")).sum()),
        "predicted_sequon_split_extracellular": has(out, "splits_predicted_sequon"),
        "predicted_sequon_beside_only": has(out, "sequon_n"),
        "flag_only_arms": len(flags),
        "replacement_found": int(rp["replacement"].notna().sum()),
        "replacement_basis": rp["replacement_basis"].value_counts().to_dict(),
        "replacement_already_folded": int(rp["replacement_already_folded"].fillna(False)
                                          .astype(bool).sum()),
        "constructs_to_fold": len(fold),
        "residues_to_fold": aa,
        "cost_usd_range": [round(aa * USD_PER_RESIDUE[0], 2), round(aa * USD_PER_RESIDUE[1], 2)],
    }  # fmt: skip
    OUT_JSON.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
