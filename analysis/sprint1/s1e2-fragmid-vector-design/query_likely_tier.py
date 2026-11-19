"""s1e2 — 'likely'-tier surfaceome genes that carry a membrane anchor.

Runs PR #130's `passes_likely` over the deep-dive cohort, takes the genes that are
likely-but-not-canonical, and keeps those with >=1 TM helix OR GPI anchoring -- i.e.
things that could actually be screened. Then applies the same screen-ability filters
used on the canonical set (no co-receptor requirement, ECD != none, <=5 kb, no OR/TAS).

Result (README 5.9): 625 likely-only -> 444 with an anchor -> 330 after filters,
499,710 bp, ~$19,931 under scenario-D routing. Best marginal value in the analysis.

Read-only: uses D1Client.public(). Needs accessible-surfaceome on sys.path + its .env.
"""

import importlib.util
import json
import re
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
]
sel = ", ".join(f"json_extract(annotation_json,'$.filters.{f}') AS {f}" for f in F)
with D1Client.public() as c:
    rows = c.query(f"SELECT gene_symbol, uniprot_acc, {sel} FROM surface_annotation")
    topo = c.query("""SELECT gene_symbol, uniprot_acc, protein_length, tm_helix_count,
                           deeptmhmm_label, ecd_length_residues, topology_version
                    FROM topology_public WHERE cohort='human_canonical' AND is_canonical=1""")


def ded(rs, k):
    b = {}
    for r in rs:
        kk = r[k]
        if kk not in b or (r.get("topology_version") or "") > (b[kk].get("topology_version") or ""):
            b[kk] = r
    return list(b.values())


topo = ded(topo, "uniprot_acc")
t_by_acc = {r["uniprot_acc"]: r for r in topo}
t_by_sym = {r["gene_symbol"]: r for r in topo}


def norm(v):
    return bool(v) if v in (0, 1) else v


canon = set()
likely = set()
facet = {}
for r in rows:
    f = {
        k: (
            norm(r[k])
            if k in ("low_endogenous_expression", "requires_coreceptor_for_expression")
            else r[k]
        )
        for k in F
    }
    facet[r["gene_symbol"]] = dict(f, acc=r["uniprot_acc"])
    if m.passes_canonical(f):
        canon.add(r["gene_symbol"])
    elif m.passes_likely(f):
        likely.add(r["gene_symbol"])
print(f"canonical={len(canon)}  likely-only (not canonical)={len(likely)}")

P_POOL, P_ABIO, P_FRAG = 0.00647, 0.01, 0.07


def cost(Ls):
    return sum(L * (P_ABIO if L <= 1200 else (P_POOL if L <= 1800 else P_FRAG)) for L in Ls)


def gpi(g):
    f = facet.get(g, {})
    return f.get("subcategory") == "GPI_anchored" or f.get("surface_call_reason") == "gpi_anchored"


cands = []
for g in likely:
    f = facet[g]
    t = t_by_acc.get(f["acc"]) or t_by_sym.get(g)
    if not t or not t["protein_length"]:
        continue
    tm = t["tm_helix_count"] or 0
    if tm >= 1 or gpi(g):
        cands.append(
            dict(
                gene=g,
                cds=t["protein_length"] * 3 + 3,
                tm=tm,
                gpi=gpi(g),
                label=t["deeptmhmm_label"],
                ecd=t["ecd_length_residues"],
                coreq=f.get("requires_coreceptor_for_expression"),
                ecdclass=f.get("ecd_accessibility_class"),
                sub=f.get("subcategory"),
            )
        )
print(f"\nlikely-only with >=1 TM helix OR GPI-anchored: {len(cands)}")
print(f"  total bp {sum(c['cds'] for c in cands):,}   cost ${cost([c['cds'] for c in cands]):,.0f}")
print("  subcategory:", dict(Counter(c["sub"] for c in cands).most_common(8)))
print("  TM>=1:", sum(1 for c in cands if c["tm"] >= 1), " GPI:", sum(1 for c in cands if c["gpi"]))

# apply the same screen-ability filters we used on canonical
clean = [c for c in cands if not c["coreq"] and c["ecdclass"] != "none" and c["cds"] <= 5000]
ORTAS = re.compile(r"^(OR\d|TAS1R|TAS2R|VN1R)")
clean = [c for c in clean if not ORTAS.match(c["gene"])]
print("\nafter the same screen-ability filters (no co-receptor req, ECD!=none, <=5kb, no OR/TAS):")
print(
    f"  {len(clean)} genes  {sum(c['cds'] for c in clean):,} bp  ${cost([c['cds'] for c in clean]):,.0f}"
)
print("\n  20 shortest (cheapest to add):")
for c in sorted(clean, key=lambda x: x["cds"])[:20]:
    print(
        f"    {c['gene']:<12} {c['cds']:>6} bp  TM={c['tm']:<2} gpi={str(c['gpi'])[:1]}  {c['label']:<6} ecd={c['ecd']}"
    )
json.dump(clean, open("/tmp/likely_add.json", "w"))
