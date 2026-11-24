"""s1e2 — pull the canonical set, isoforms and orthologs from the read-only public D1.

Regenerates the inputs behind READMEs 5.4 / 5.6 / 5.7:
  * the 2,243 canonical genes (accessible-surfaceome PR #130 `passes_canonical`)
  * protein lengths + DeepTMHMM topology (`topology_public`)
  * non-secreted isoforms differing substantially from canonical
  * mouse / cyno orthologs with ECD percent identity (`compara_ortholog_ecd`)

Requires the accessible-surfaceome repo on sys.path and its .env (Cloudflare D1
credentials). Read-only: uses `D1Client.public()`, the column-whitelisted mirror.

TWO GOTCHAS this script handles, both of which silently corrupt naive queries:
  1. `topology_public` has MULTIPLE `topology_version` rows per sequence -> dedupe
     to the latest per `uniprot_acc_full`, or every count doubles.
  2. Orthologs must join via `compara_ortholog_ecd.ortholog_uniprot_acc`. Matching on
     upper-cased ortholog gene symbols finds 302 mouse orthologs; the key finds 1,777.
"""

import json
import sys

sys.path.insert(0, "src")
try:
    from accessible_surfaceome.env import load_env

    load_env()
except Exception:
    pass
from accessible_surfaceome.cloud.d1_client import D1Client

canon = json.load(open("/tmp/canon.json"))
gene2acc = {g: a for g, a in canon if a}
accs = set(gene2acc.values())
P_POOL, P_ABIO, P_FRAG = 0.00647, 0.01, 0.07


def cost_of(L):
    return sum(x * (P_ABIO if x <= 1200 else (P_POOL if x <= 1800 else P_FRAG)) for x in L)


with D1Client.public() as c:
    print("=== topology_version values ===")
    for r in c.query(
        "SELECT cohort, topology_version, COUNT(*) n FROM topology_public GROUP BY cohort, topology_version"
    ):
        print("  ", r)
    canon_rows = c.query("""SELECT uniprot_acc, protein_length, ecd_length_residues, topology_version
                          FROM topology_public WHERE cohort='human_canonical' AND is_canonical=1""")
    iso = c.query("""SELECT gene_symbol, uniprot_acc, uniprot_acc_full, protein_length,
                          ecd_length_residues, deeptmhmm_label, topology_version
                   FROM topology_public WHERE cohort='human_isoforms'""")
    orth = c.query("""SELECT cohort, uniprot_acc, protein_length, deeptmhmm_label, topology_version
                    FROM topology_public WHERE cohort IN ('mouse_ortholog','cyno_ortholog')""")
    ecd = c.query("""SELECT human_uniprot_acc, species, ortholog_uniprot_acc, ecd_pct_identity
                   FROM compara_ortholog_ecd WHERE ecd_pct_identity IS NOT NULL""")


def dedupe(rows, key):
    best = {}
    for r in rows:
        k = r[key]
        if k not in best or (r.get("topology_version") or "") > (
            best[k].get("topology_version") or ""
        ):
            best[k] = r
    return list(best.values())


canon_rows = dedupe(canon_rows, "uniprot_acc")
iso = dedupe(iso, "uniprot_acc_full")
orth = dedupe(orth, "uniprot_acc")
cb = {r["uniprot_acc"]: r for r in canon_rows if r["uniprot_acc"] in accs}
print(f"\ndeduped: canonical={len(cb)}  isoforms={len(iso)}  orthologs={len(orth)}")

print("\n" + "=" * 72)
print("1) ISOFORMS — our genes, non-secreted-only, substantially different from canonical")
cands = []
for r in iso:
    b = r["uniprot_acc"]
    if b not in cb or r["deeptmhmm_label"] == "SP":
        continue
    cl, il = cb[b]["protein_length"], r["protein_length"]
    if not cl or not il:
        continue
    d = abs(il - cl)
    pct = 100 * d / cl
    cands.append(
        dict(
            gene=r["gene_symbol"],
            iso=r["uniprot_acc_full"],
            canon_len=cl,
            iso_len=il,
            d_abs=d,
            d_pct=pct,
            label=r["deeptmhmm_label"],
            ecd=r["ecd_length_residues"],
            canon_ecd=cb[b]["ecd_length_residues"],
        )
    )
print(
    f"  non-SP isoforms of our genes (deduped): {len(cands)} across {len({c['gene'] for c in cands})} genes"
)
for pc, ab in ((10, 50), (20, 100), (30, 150)):
    s = [c for c in cands if c["d_pct"] >= pc and c["d_abs"] >= ab]
    L = [c["iso_len"] * 3 + 3 for c in s]
    print(
        f"   >={pc}% AND >={ab} aa : {len(s):5d} isoforms ({len({c['gene'] for c in s})} genes)  "
        f"{sum(L):>9,} bp  ${cost_of(L):>10,.0f}"
    )
sel = [c for c in cands if c["d_pct"] >= 20 and c["d_abs"] >= 100]
json.dump(sel, open("/tmp/iso_substantial.json", "w"))

print("\n" + "=" * 72)
print("2) ORTHOLOGS — joined via compara_ortholog_ecd.ortholog_uniprot_acc (proper key)")
olen = {r["uniprot_acc"]: r for r in orth}
for sp in ("mouse", "cynomolgus"):
    recs = [r for r in ecd if r["species"] == sp and r["human_uniprot_acc"] in accs]
    print(f"\n  {sp}: {len(recs)} ECD-identity records for our genes")
    for cut in (0, 70, 80, 90, 95):
        keep = []
        for r in recs:
            if cut and (r["ecd_pct_identity"] or 0) < cut:
                continue
            t = olen.get(r["ortholog_uniprot_acc"])
            if t and t["protein_length"] and t["deeptmhmm_label"] != "SP":
                keep.append(t["protein_length"])
        L = [x * 3 + 3 for x in keep]
        lbl = "all (non-SP)" if cut == 0 else f">={cut}% ECD id"
        print(f"    {lbl:<16} {len(keep):5d} orthologs  {sum(L):>9,} bp  ${cost_of(L):>10,.0f}")
