"""s1e2 — audit of deep-dived genes with a membrane anchor that were REJECTED.

Runs PR #130's passes_canonical / passes_likely over the 5,130-gene deep-dive cohort,
takes the ~2,262 rejects, and keeps those with >=1 TM helix or a GPI anchor (1,129).
Reports why each failed and surfaces the concerning subset: high/moderate confidence +
surface_dominant|mixed + accessible, which fail only on evidence_grade=weak /
evidence_density=low -- i.e. corpus-coverage failures rather than biology failures.

See REPORT_library_scoping.md section 1. Read-only: D1Client.public().
"""

import importlib.util
import json
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, "src")
try:
    from accessible_surfaceome.env import load_env

    load_env()
except Exception:
    pass
from accessible_surfaceome.cloud.d1_client import D1Client

# The committed copy in S1E3, not the /tmp scratch copy this was first run against.
PRESETS = Path(__file__).resolve().parents[1] / "s1e3-surface-library-filters" / "presets_pr130.py"
spec = importlib.util.spec_from_file_location("p", PRESETS)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

F = [
    "evidence_grade",
    "confidence",
    "surface_specificity",
    "state_dependence",
    "surface_accessibility",
    "evidence_density",
    "low_endogenous_expression",
    "subcategory",
    "surface_call_reason",
    "requires_coreceptor_for_expression",
    "ecd_accessibility_class",
    "llm_family",
]
sel = ", ".join(f"json_extract(annotation_json,'$.filters.{f}') AS {f}" for f in F)
with D1Client.public() as c:
    rows = c.query(f"SELECT gene_symbol, uniprot_acc, {sel} FROM surface_annotation")
    topo = c.query("""SELECT gene_symbol, uniprot_acc, protein_length, tm_helix_count, deeptmhmm_label,
                           ecd_length_residues, topology_version
                    FROM topology_public WHERE cohort='human_canonical' AND is_canonical=1""")


def ded(rs, k):
    b = {}
    for r in rs:
        kk = r[k]
        if kk not in b or (r.get("topology_version") or "") > (b[kk].get("topology_version") or ""):
            b[kk] = r
    return list(b.values())


topo = ded(topo, "uniprot_acc")
ta = {r["uniprot_acc"]: r for r in topo}
ts = {r["gene_symbol"]: r for r in topo}


def norm(v):
    return bool(v) if v in (0, 1) else v


recs = {}
for r in rows:
    f = {
        k: (
            norm(r[k])
            if k in ("low_endogenous_expression", "requires_coreceptor_for_expression")
            else r[k]
        )
        for k in F
    }
    recs[r["gene_symbol"]] = dict(f, acc=r["uniprot_acc"])


def gpi(f):
    return f.get("subcategory") == "GPI_anchored" or f.get("surface_call_reason") == "gpi_anchored"


canon = set()
likely = set()
for g, f in recs.items():
    if m.passes_canonical(f):
        canon.add(g)
    elif m.passes_likely(f):
        likely.add(g)
rejected = set(recs) - canon - likely
print(
    f"deep-dive cohort={len(recs)}  canonical={len(canon)}  likely-only={len(likely)}  REJECTED={len(rejected)}"
)

# rejected WITH a membrane anchor
anch = []
for g in rejected:
    f = recs[g]
    t = ta.get(f["acc"]) or ts.get(g)
    tm = (t or {}).get("tm_helix_count") or 0
    if tm >= 1 or gpi(f):
        anch.append(
            dict(
                gene=g,
                tm=tm,
                gpi=gpi(f),
                acc=f["acc"],
                cds=((t or {}).get("protein_length") or 0) * 3 + 3,
                label=(t or {}).get("deeptmhmm_label"),
                ecd=(t or {}).get("ecd_length_residues"),
                sub=f.get("subcategory"),
                conf=f.get("confidence"),
                eg=f.get("evidence_grade"),
                spec=f.get("surface_specificity"),
                acc_cls=f.get("surface_accessibility"),
                dens=f.get("evidence_density"),
                sd=f.get("state_dependence"),
                fam=f.get("llm_family"),
            )
        )
print(f"\nREJECTED with >=1 TM helix or GPI anchor: {len(anch)}")
print(f"  TM>=1: {sum(1 for a in anch if a['tm'] >= 1)}   GPI: {sum(1 for a in anch if a['gpi'])}")
print("\n  why were they rejected? (the failing gate)")


def fails(f):
    out = []
    if f.get("evidence_grade") not in (
        "direct_multi_method",
        "direct_single_method",
        "supportive_but_indirect",
        "weak",
    ):
        out.append("evidence_grade=conflicting")
    if f.get("confidence") not in ("high", "moderate"):
        out.append(f"confidence={f.get('confidence')}")
    if f.get("surface_specificity") not in ("surface_dominant", "mixed", "mostly_intracellular"):
        out.append(f"specificity={f.get('surface_specificity')}")
    if f.get("surface_accessibility") not in ("high", "moderate", "low"):
        out.append(f"accessibility={f.get('surface_accessibility')}")
    return out or ["(passes likely gates - check evidence_grade/density)"]


ctr = Counter()
for a in anch:
    for r_ in fails(recs[a["gene"]]):
        ctr[r_] += 1
for k, v in ctr.most_common(12):
    print(f"    {k:<46} {v}")
print("\n  by subcategory:", dict(Counter(a["sub"] for a in anch).most_common(10)))
print("  by confidence :", dict(Counter(a["conf"] for a in anch).most_common()))
print("  by specificity:", dict(Counter(a["spec"] for a in anch).most_common()))
print("  by accessibility:", dict(Counter(a["acc_cls"] for a in anch).most_common()))
print("\n  20 largest rejected-with-anchor (potentially concerning):")
for a in sorted(anch, key=lambda x: -(x["ecd"] or 0))[:20]:
    print(
        f"    {a['gene']:<12} TM={a['tm']:<2} gpi={str(a['gpi'])[:1]} ecd={a['ecd']:<5} {a['label'] or '?':<6} "
        f"conf={a['conf']:<9} spec={str(a['spec'])[:18]:<18} acc={str(a['acc_cls'])[:9]}"
    )
json.dump(anch, open("/tmp/rejected_anchored.json", "w"))
