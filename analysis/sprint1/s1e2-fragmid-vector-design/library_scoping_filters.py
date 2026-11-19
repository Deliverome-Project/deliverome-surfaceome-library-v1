"""s1e2 — library scoping filters for the 2,243 canonical surfaceome ORFs.

Computes the cost of removing categories that cannot be scored in a functional-genomics
surface screen (README 5.6). Filters:
  1. secreted-only and NOT GPI-anchored  (GPI proteins are DeepTMHMM 'SP' but ARE
     surface-displayed, so they must be spared -- detected via the deep-dive
     subcategory=GPI_anchored / surface_call_reason=gpi_anchored facets)
  2. GLOB (no TM helix, no extracellular residues)
  3. requires_coreceptor_for_expression -- cannot reach the surface from single-gene
     overexpression, so would score as a false negative
  4. ecd_accessibility_class = none -- nothing extracellular to engage
  5. olfactory / taste receptors -- ER-retained in heterologous cells without accessory
     factors (RTP1/RTP2 for ORs)
  6. CDS > 5 kb -- beyond the practical single-fragment ceiling

Also applies the isoform rule: drop an isoform when the canonical has >=1 TM helix and
the isoform has 0 (a truncation that lost membrane anchoring is not a testable surface
variant).

Reads /tmp/facets.json + /tmp/canon_topology.json produced by query_isoforms_orthologs.py
and the facet query; needs accessible-surfaceome on sys.path for the D1 client.
"""

import json
import re
import sys

sys.path.insert(0, "src")
try:
    from accessible_surfaceome.env import load_env

    load_env()
except Exception:
    pass
from accessible_surfaceome.cloud.d1_client import D1Client

facets = json.load(open("/tmp/facets.json"))
topo = json.load(open("/tmp/canon_topology.json"))
for r in topo:
    r["cds"] = r["len"] * 3 + 3
byg = {r["gene"]: r for r in topo}
P_POOL, P_ABIO, P_FRAG = 0.00647, 0.01, 0.07


def cost(rows):
    return sum(
        r["cds"] * (P_ABIO if r["cds"] <= 1200 else (P_POOL if r["cds"] <= 1800 else P_FRAG))
        for r in rows
    )


base = cost(topo)
print(f"BASELINE: {len(topo)} ORFs  ${base:,.0f}\n")
print("SRC in canonical set?", "SRC" in byg)


def gpi(g):
    f = facets.get(g, {})
    return f.get("subcategory") == "GPI_anchored" or f.get("surface_call_reason") == "gpi_anchored"


FILTERS = {}
# 1. secreted-only, NOT GPI-anchored
FILTERS["1 secreted-only & NOT GPI-anchored"] = [
    r for r in topo if r["label"] == "SP" and not gpi(r["gene"])
]
# 2. GLOB (no TM, no ECD)
FILTERS["2 GLOB (no TM, no ECD)"] = [r for r in topo if r["label"] == "GLOB"]
# 3. requires a co-receptor to reach the surface
FILTERS["3 requires_coreceptor_for_expression"] = [
    r
    for r in topo
    if facets.get(r["gene"], {}).get("requires_coreceptor_for_expression") in (1, "1", True)
]
# 4. no extracellular domain to engage
FILTERS["4 ecd_accessibility_class = none"] = [
    r for r in topo if facets.get(r["gene"], {}).get("ecd_accessibility_class") == "none"
]
FILTERS["4b ecd_accessibility_class none+minimal"] = [
    r
    for r in topo
    if facets.get(r["gene"], {}).get("ecd_accessibility_class") in ("none", "minimal")
]
# 5. olfactory / taste receptors (poor heterologous surface expression)
ORTAS = re.compile(r"^(OR\d|OR\d+[A-Z]|TAS1R|TAS2R|VN1R)")
FILTERS["5 olfactory + taste receptors"] = [r for r in topo if ORTAS.match(r["gene"])]
# 6. beyond practical fragment ceiling
FILTERS["6 CDS > 5 kb"] = [r for r in topo if r["cds"] > 5000]

print(f"\n{'filter':<44}{'ORFs':>6}{'bp':>12}{'cost saved':>13}")
print("-" * 76)
for k, v in FILTERS.items():
    print(f"{k:<44}{len(v):>6}{sum(r['cds'] for r in v):>12,}{cost(v):>13,.0f}")

# union of the recommended set (1,2,3,4,5,6)
rec_keys = [
    "1 secreted-only & NOT GPI-anchored",
    "2 GLOB (no TM, no ECD)",
    "3 requires_coreceptor_for_expression",
    "4 ecd_accessibility_class = none",
    "5 olfactory + taste receptors",
    "6 CDS > 5 kb",
]
drop = {r["gene"] for k in rec_keys for r in FILTERS[k]}
keep = [r for r in topo if r["gene"] not in drop]
print(f"\nUNION of recommended filters: drop {len(drop)} ORFs -> keep {len(keep)}")
print(
    f"  cost ${cost(keep):,.0f}   (saves ${base - cost(keep):,.0f}, {100 * (base - cost(keep)) / base:.1f}%)"
)
print(
    f"\n  GPI-anchored SPARED by filter 1: {sum(1 for r in topo if r['label'] == 'SP' and gpi(r['gene']))}"
)
json.dump(sorted(drop), open("/tmp/drop_genes.json", "w"))
json.dump([r["gene"] for r in keep], open("/tmp/keep_genes.json", "w"))

# ---- ISOFORM filter: canonical has TM, isoform lost it ----
print("\n" + "=" * 76)
print("ISOFORMS — drop those that lose the TM anchor relative to canonical")
with D1Client.public() as c:
    iso = c.query("""SELECT gene_symbol, uniprot_acc, uniprot_acc_full, protein_length,
                          ecd_length_residues, deeptmhmm_label, tm_helix_count, topology_version
                   FROM topology_public WHERE cohort='human_isoforms'""")
    cn = c.query("""SELECT uniprot_acc, protein_length, tm_helix_count, deeptmhmm_label, topology_version
                  FROM topology_public WHERE cohort='human_canonical' AND is_canonical=1""")


def ded(rows, k):
    b = {}
    for r in rows:
        kk = r[k]
        if kk not in b or (r.get("topology_version") or "") > (b[kk].get("topology_version") or ""):
            b[kk] = r
    return list(b.values())


iso = ded(iso, "uniprot_acc_full")
cn = ded(cn, "uniprot_acc")
g2a = json.load(open("/tmp/canon.json"))
g2a = {g: a for g, a in g2a if a}
accs = set(g2a.values())
cbа = {r["uniprot_acc"]: r for r in cn if r["uniprot_acc"] in accs}
cands = []
for r in iso:
    b = r["uniprot_acc"]
    if b not in cbа or r["deeptmhmm_label"] == "SP":
        continue
    cl = cbа[b]["protein_length"]
    il = r["protein_length"]
    if not cl or not il:
        continue
    d = abs(il - cl)
    pct = 100 * d / cl
    if pct < 20 or d < 100:
        continue
    cands.append(
        dict(
            gene=r["gene_symbol"],
            iso=r["uniprot_acc_full"],
            cl=cl,
            il=il,
            d=d,
            pct=pct,
            c_tm=cbа[b]["tm_helix_count"],
            i_tm=r["tm_helix_count"],
            c_lab=cbа[b]["deeptmhmm_label"],
            i_lab=r["deeptmhmm_label"],
            ecd=r["ecd_length_residues"],
        )
    )
print(f"  starting tier (>=20% & >=100 aa, non-SP): {len(cands)} isoforms")
lost = [c for c in cands if (c["c_tm"] or 0) >= 1 and (c["i_tm"] or 0) == 0]
kept = [c for c in cands if not ((c["c_tm"] or 0) >= 1 and (c["i_tm"] or 0) == 0)]
Lk = [c["il"] * 3 + 3 for c in kept]
Ll = [c["il"] * 3 + 3 for c in lost]


def costL(L):
    return sum(x * (P_ABIO if x <= 1200 else (P_POOL if x <= 1800 else P_FRAG)) for x in L)


print(f"  DROP (canonical has TM, isoform has 0 TM): {len(lost)}  ${costL(Ll):,.0f}")
print(
    f"  KEEP: {len(kept)} isoforms across {len({c['gene'] for c in kept})} genes  {sum(Lk):,} bp  ${costL(Lk):,.0f}"
)
print("\n  examples dropped (TM lost):")
for c in sorted(lost, key=lambda x: -x["d"])[:8]:
    print(
        f"    {c['gene']:<10} {c['iso']:<14} {c['cl']}->{c['il']} aa  TM {c['c_tm']}->{c['i_tm']}  {c['c_lab']}->{c['i_lab']}"
    )
