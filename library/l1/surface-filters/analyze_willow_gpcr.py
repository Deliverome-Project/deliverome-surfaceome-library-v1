#!/usr/bin/env python3
"""Validate the surface-library filter against Willow's GPCR surface-expression
screen (Tedman et al. 2025, bioRxiv 2025.09.19.677468).

For canonical GPCRs (clean gene-symbol join) and their isoforms (approximate
gene+length join, since the screen uses Ensembl transcript IDs and our catalog
uses UniProt isoform IDs), report:
  - how many screened GPCRs are in our cohort and pass our filters (by class/tier),
  - how many screened isoforms "make it through" (parent in library + our topology
    calls the matched isoform surface),
  - the discordances: screen shows well-expressed at the surface but we excluded.

Run: uv run --locked python analyze_willow_gpcr.py
"""
# ruff: noqa: E402, E741

import json
import os

import build_library as B
import final_library as F
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
XLSX = os.path.join(HERE, "external/willow_gpcr_screen/tedman2025_gpcr_screen_media-2.xlsx")

recs = {json.loads(l)["gene"]: json.loads(l) for l in open(os.path.join(HERE, "cohort.jsonl"))}


# ---- filter logic: reuse final_library's LIVE gate stack (single source of truth,
# so this validation cannot drift from the real filter). gate_results composes the
# per-gate booleans from final_library's atomic gate functions; membership + ecd are
# now included (the old local copy omitted them and used TPM>=10 / no allowlist). ----
def gate_results(g, r):
    f = r["filters"]
    return {
        "membership": B.passes_canonical(f) or B.passes_likely(f),
        "acc": f.get("surface_accessibility") != "no",
        "spec": f.get("surface_specificity") in ("surface_dominant", "mixed"),
        "anchored": F.anchored(r),
        "coreceptor": F.coreceptor_ok(g, f),
        "expr": F.expr_ok(g, f),
        "reason": f.get("surface_call_reason") not in B.NONSURFACE_REASONS,
        "orf": F.orf_ok(r),
        "ecd": F.ecd_ok(r),
    }


def passes(g, r):
    return all(gate_results(g, r).values())


def tier(r):
    return F.tier(r)


# ======================= CANONICAL =======================
can = pd.read_excel(XLSX, sheet_name="Canonical")
can = can[can["Receptor Name"].notna()].copy()
# one row per gene: max surface immunostaining across its constructs
g_int = can.groupby("gene_id")["Immunostaining Intensity"].max()
g_cls = can.groupby("gene_id")["GPCR Class"].first()
genes = sorted(g_int.index)

WELL = 1000.0  # surface-positive threshold: immunostaining >= 1000 = positive, else negative
STRONG = 1155.4  # screen 75th pctile (secondary "strongly expressed" marker)

rows = []
for g in genes:
    r = recs.get(g)
    intens = float(g_int[g])
    cls = g_cls[g]
    olf = ("Olfactory" in str(cls)) or (cls in ("Taste", "Vomeronasal"))
    if r is None:
        rows.append(
            dict(
                gene=g,
                cls=cls,
                olf=olf,
                intens=intens,
                in_cohort=False,
                passes=None,
                tier=None,
                fails="",
            )
        )
        continue
    gr = gate_results(g, r)
    rows.append(
        dict(
            gene=g,
            cls=cls,
            olf=olf,
            intens=intens,
            in_cohort=True,
            passes=all(gr.values()),
            tier=tier(r) if all(gr.values()) else None,
            fails=",".join(k for k, v in gr.items() if not v),
        )
    )
cdf = pd.DataFrame(rows)

n = len(cdf)
in_coh = cdf["in_cohort"].sum()
passed = (cdf["passes"] == True).sum()  # noqa: E712
print(f"===== CANONICAL GPCRs: {n} screened genes =====")
print(f"in our cohort (accessible-surfaceome 5,130): {in_coh}  |  not in cohort: {n - in_coh}")
print(
    f"PASS our filters: {passed}  ({passed / n * 100:.0f}% of all screened, {passed / in_coh * 100:.0f}% of in-cohort)"
)
print(f"in cohort but FILTERED OUT: {in_coh - passed}")

print("\n-- by GPCR class --")
for olf_flag, label in [(False, "non-olfactory (A/B/C/F)"), (True, "olfactory/taste/vomeronasal")]:
    sub = cdf[cdf["olf"] == olf_flag]
    p = (sub["passes"] == True).sum()  # noqa: E712
    ic = sub["in_cohort"].sum()
    print(f"  {label:32}: screened {len(sub):3}  in-cohort {ic:3}  pass {p:3}")

print("\n-- tier of the passers --")
print(cdf[cdf["passes"] == True]["tier"].value_counts().to_dict())  # noqa: E712

print("\n-- why the in-cohort GPCRs were filtered out (gate failure counts; a gene can fail >1) --")
from collections import Counter

fc = Counter()
for s in cdf[(cdf["in_cohort"]) & (cdf["passes"] == False)]["fails"]:  # noqa: E712
    for k in s.split(","):
        fc[k] += 1
print(dict(fc.most_common()))

# ---- FALSE NEGATIVES: well surface-expressed in the screen but excluded ----
excluded = cdf[(cdf["passes"] != True)].copy()  # noqa: E712  (fail OR not-in-cohort)
fn = excluded[excluded["intens"] >= WELL].sort_values("intens", ascending=False)
print(
    f"\n===== WELL-EXPRESSED (immunostaining >= {WELL:.0f}, screen median) but NOT in our library: {len(fn)} ====="
)
print(
    f"   of which strongly expressed (>= {STRONG:.0f}, 75th pctile): {(fn['intens'] >= STRONG).sum()}"
)
print(
    f"   split: not-in-cohort {(~fn['in_cohort']).sum()}  |  in-cohort-but-filtered {fn['in_cohort'].sum()}"
)
print("\n  top 30 by surface intensity:")
print(f"  {'gene':10}{'class':16}{'intensity':>10}  {'status'}")
for _, x in fn.head(30).iterrows():
    status = "NOT in cohort" if not x["in_cohort"] else f"filtered: {x['fails']}"
    print(f"  {x['gene']:10}{x['cls']!s:16}{x['intens']:>10.0f}  {status}")

cdf.to_csv(os.path.join(HERE, "willow_canonical_vs_filter.csv"), index=False)
fn.to_csv(os.path.join(HERE, "willow_wellexpressed_excluded.csv"), index=False)

# ======================= ISOFORMS =======================
iso = pd.read_excel(XLSX, sheet_name="Isoforms")
iso["gene"] = iso["Isoform Transcript"].str.split("_").str[0]


def our_isoform_topo(gene, length_aa):
    """Approx-match a screen isoform (gene + AA length) to our catalog isoform by
    closest length (<=5 aa), return our topology_status ('surface'/...), or None."""
    r = recs.get(gene)
    if not r:
        return None, None
    sub = r["filters"].get("subcategory")
    best, bestd = None, 999
    for it in r.get("isoforms", []):
        L = len(it.get("per_residue_topology") or "")
        d = abs(L - length_aa)
        if d < bestd:
            bestd, best = d, it
    if best is None or bestd > 5:
        return None, bestd if best else None
    st = B.topology_status(
        best.get("tm_helix_count"),
        best.get("signal_peptide_length"),
        best.get("per_residue_topology"),
        sub,
    )
    return st, bestd


irows = []
for _, x in iso.iterrows():
    g = x["gene"]
    r = recs.get(g)
    parent_pass = passes(g, r) if r else False
    st, d = our_isoform_topo(g, int(x["length (AA)"]))
    made_through = bool(parent_pass and st == "surface")
    irows.append(
        dict(
            transcript=x["Isoform Transcript"],
            gene=g,
            length=int(x["length (AA)"]),
            surf=float(x["Surface Expression"]),
            pct_change=float(x["% change in PME"]),
            parent_in_cohort=r is not None,
            parent_pass=parent_pass,
            our_iso_match=st is not None,
            our_topo=st,
            made_through=made_through,
        )
    )
idf = pd.DataFrame(irows)

print(f"\n\n===== ISOFORMS: {len(idf)} screened transcripts, {idf['gene'].nunique()} genes =====")
print(
    f"parent gene in our cohort: {idf['parent_in_cohort'].sum()}  |  parent passes filters: {idf['parent_pass'].sum()}"
)
print(f"matched to one of our catalog isoforms (gene+len<=5aa): {idf['our_iso_match'].sum()}")
print(
    f"'make it through' (parent passes AND our topology calls it surface): {idf['made_through'].sum()}"
)
print(
    f"screen says isoform loses surface vs canonical (% change < 0): {(idf['pct_change'] < 0).sum()} / {len(idf)}"
)

iwell = idf[idf["surf"] >= WELL].sort_values("surf", ascending=False)
imiss = iwell[~iwell["made_through"]]
print(f"\n-- well surface-expressed isoforms (Surface Expression >= {WELL:.0f}): {len(iwell)} --")
print(f"   of those, did NOT make it through our library: {len(imiss)}")
print("   top 20 well-expressed isoforms we don't capture:")
print(f"   {'transcript':30}{'gene':9}{'surf':>8}{'%chg':>8}  {'reason we miss it'}")
for _, x in imiss.head(20).iterrows():
    if not x["parent_in_cohort"]:
        why = "parent not in cohort"
    elif not x["parent_pass"]:
        why = "parent gene filtered out"
    elif not x["our_iso_match"]:
        why = "isoform absent from our catalog"
    elif x["our_topo"] != "surface":
        why = f"our topology calls it '{x['our_topo']}'"
    else:
        why = "?"
    print(f"   {x['transcript']:30}{x['gene']:9}{x['surf']:>8.0f}{x['pct_change']:>8.0f}  {why}")

idf.to_csv(os.path.join(HERE, "willow_isoforms_vs_filter.csv"), index=False)
print(
    "\nSaved: willow_canonical_vs_filter.csv, willow_wellexpressed_excluded.csv, willow_isoforms_vs_filter.csv"
)
