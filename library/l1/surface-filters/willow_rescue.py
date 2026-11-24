#!/usr/bin/env python3
"""Screen-based inclusion + isoform concordance from Willow's GPCR screen.

Produces:
  - willow_screen_rescue_genes.json : canonical GPCRs surface-expressed in the
    screen (immunostaining >= SCREEN_MIN) that our filters excluded -> add to the
    library as a flagged "SV" (screen-validated) set.  final_library.py reads this.
  - willow_screen_rescue_isoforms.json : native (non-SP-HA) isoforms surface-
    expressed in the screen whose parent gene is in the library.
  - willow_isoform_concordance.csv : per screen-isoform, ENST->UniProt-isoform
    (via the map), our catalog topology call vs the screen surface readout.

Run: uv run --locked python willow_rescue.py
"""
# ruff: noqa: E741, RUF003  (∪ used in comments)

import csv
import json
import os

import build_library as B
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
XLSX = os.path.join(HERE, "external/willow_gpcr_screen/tedman2025_gpcr_screen_media-2.xlsx")

recs = {json.loads(l)["gene"]: json.loads(l) for l in open(os.path.join(HERE, "cohort.jsonl"))}
# inject the synthesizer's holistic grade so PR137 membership evaluates faithfully
_egs_path = os.path.join(HERE, "evidence_grade_summary.json")
if os.path.exists(_egs_path):
    _egs = json.load(open(_egs_path))
    for _r in recs.values():
        if _r["filters"].get("evidence_grade_summary"):
            continue  # refreshed cohort already carries it
        _s = _egs.get(_r["gene"])
        if _s and _s != "__ERROR__":
            _r["filters"]["evidence_grade_summary"] = _s
expr = json.load(open(os.path.join(HERE, "external/expression_tpm/resolved_expr.json")))
mxfull = json.load(open(os.path.join(HERE, "external/expression_tpm/resolved_expr_full.json")))
cohort_genes = set(recs)
emap = json.load(open(os.path.join(HERE, "willow_enst_uniprot_isoform_map.json")))
# Co-receptor partners from the cohort record (API field captured by fetch_cohort's
# compact()); legacy coreceptor_partners.jsonl side-file is the fallback.
part = {
    g: {"partners": r["co_receptor_partners"]}
    for g, r in recs.items()
    if r.get("co_receptor_partners")
}
if not part and os.path.exists(os.path.join(HERE, "coreceptor_partners.jsonl")):
    part = {
        o["gene"]: o
        for o in (json.loads(l) for l in open(os.path.join(HERE, "coreceptor_partners.jsonl")))
        if "partners" in o
    }
med_ens, med_sym = {}, {}
for r in csv.DictReader(
    open(os.path.join(HERE, "external/hpa_rna_celline/hpa_celline_per_gene.csv"))
):
    try:
        med_ens[r["ensembl_base"]] = float(r["median_tpm_celline"])
        med_sym[r["gene_name"]] = float(r["median_tpm_celline"])
    except ValueError:
        pass

SCREEN_MIN = (
    1000.0  # surface-positive threshold: anti-HA immunostaining >= 1000 = positive, else negative
)


def mx(g):
    return mxfull.get(g)  # full re-derived reference


def _psym(p):
    return str(p).replace("(", " ").split()[0] if str(p).split() else str(p)


def partner_median(p):
    for e in (expr.get(p) or {}).get("ensembl", []):
        if e.split(".")[0] in med_ens:
            return med_ens[e.split(".")[0]]
    return med_sym.get(p)


def gate_results(g, r):
    f = r["filters"]
    anchored = (
        f.get("subcategory") == "GPI_anchored"
        or ((r.get("canonical_topology") or {}).get("tm_helix_count") or 0) >= 1
    )
    if f.get("co_receptor_dependency") != "required":
        cok = True
    else:
        o = part.get(g)
        ps = o["partners"] if o else []
        non_self = [
            p
            for p in ps
            if not (_psym(p) == g or "homo" in str(p).lower() or "self" in str(p).lower())
        ]
        if o and ps and not non_self:
            cok = True  # all partners self/homo-oligomer -> satisfied
        else:  # co-deliverability rejected -> partner must be endogenously ubiquitous
            meds = [m for m in (partner_median(p) for p in non_self) if m is not None]
            cok = bool(meds) and all(m >= 5 for m in meds)
    eok = f.get("low_endogenous_expression") is False or (mx(g) is not None and mx(g) >= 5)
    t = (r.get("canonical_topology") or {}).get("per_residue_topology")
    return {
        "membership": B.passes_canonical(f) or B.passes_likely(f),  # PR137 canonical ∪ likely
        "acc": f.get("surface_accessibility") != "no",
        "spec": f.get("surface_specificity") in ("surface_dominant", "mixed"),
        "anchored": anchored,
        "coreceptor": cok,
        "expr": eok,
        "reason": f.get("surface_call_reason") not in B.NONSURFACE_REASONS,
        "orf": t is None or (len(t) * 3 + 3) <= 5000,
    }


def passes(g, r):
    return r is not None and all(gate_results(g, r).values())


def tier_arch(r):
    f = r["filters"]
    arch = B.arch_risk(
        f, r.get("canonical_topology") or {}, r.get("uniprot_family"), r.get("hgnc_gene_groups")
    )
    constitutive = f.get("surface_call_reason") not in (
        "cell_state_induced",
        "lysosomal_exocytosis",
    )
    return (
        "T1"
        if (arch == "low" and constitutive)
        else ("T2" if (arch == "med" and constitutive) else "T3")
    )


# ---------------- canonical rescue (surface immunostaining) ----------------
# Screen GPCRs our library excluded that reach the surface in the screen
# (anti-HA immunostaining >= SCREEN_MIN) -> add as a flagged screen-validated tier.
can = pd.read_excel(XLSX, sheet_name="Canonical")
can = can[can["Receptor Name"].notna()]
g_int = can.groupby("gene_id")["Immunostaining Intensity"].max()
g_cls = can.groupby("gene_id")["GPCR Class"].first()

rescue = {}
not_in_cohort = []
for g, intens in g_int.items():
    if intens < SCREEN_MIN:
        continue
    r = recs.get(g)
    if r is None:
        not_in_cohort.append((g, float(intens), g_cls[g]))
        continue
    if passes(g, r):
        continue  # already in library
    gr = gate_results(g, r)
    rescue[g] = {
        "immunostaining": round(float(intens), 1),
        "gpcr_class": g_cls[g],
        "olfactory": ("Olfactory" in str(g_cls[g])) or (g_cls[g] in ("Taste", "Vomeronasal")),
        "prior_fail_gates": [k for k, v in gr.items() if not v],
        "arch_tier": tier_arch(r),
    }

json.dump(rescue, open(os.path.join(HERE, "willow_screen_rescue_genes.json"), "w"), indent=1)
olf = sum(1 for v in rescue.values() if v["olfactory"])
print(f"===== CANONICAL screen-rescue (immunostaining >= {SCREEN_MIN}) =====")
print(
    f"genes added to library (in cohort, previously excluded): {len(rescue)}  [{olf} olfactory/taste, {len(rescue) - olf} non-olf]"
)
print(
    f"screen-positive but NOT in our cohort (catalog gap, not added): {len(not_in_cohort)} -> {[g for g, _, _ in not_in_cohort]}"
)

# ---------------- isoform concordance via map ----------------
iso = pd.read_excel(XLSX, sheet_name="Isoforms")
iso["gene"] = iso["Isoform Transcript"].str.split("_").str[0]
iso["enst"] = iso["Isoform Transcript"].str.extract(r"(ENST\d+)")
iso["sp_ha"] = iso["Isoform Transcript"].str.contains("_SP_HA")

lib_genes = {g for g, r in recs.items() if passes(g, r)} | set(rescue)

conc_rows, rescue_iso = [], {}
matched = both_surface = topo_disagree = 0
for _, x in iso.iterrows():
    g, enst = x["gene"], x["enst"]
    m = emap.get(enst) or {}
    iso_id = m.get("isoform_id")
    our = None
    if iso_id and g in recs:
        sub = recs[g]["filters"].get("subcategory")
        for it in recs[g].get("isoforms", []):
            if it.get("isoform_id") == iso_id:
                our = B.topology_status(
                    it.get("tm_helix_count"),
                    it.get("signal_peptide_length"),
                    it.get("per_residue_topology"),
                    sub,
                )
                break
    screen_surface = x["Surface Expression"] >= SCREEN_MIN
    if our is not None:
        matched += 1
        if screen_surface and our == "surface":
            both_surface += 1
        if screen_surface and our != "surface":
            topo_disagree += 1
    conc_rows.append(
        dict(
            transcript=x["Isoform Transcript"],
            gene=g,
            enst=enst,
            uniprot_isoform=iso_id,
            sp_ha=bool(x["sp_ha"]),
            screen_surf=round(float(x["Surface Expression"]), 1),
            pct_change=round(float(x["% change in PME"]), 1),
            in_our_catalog=our is not None,
            our_topology=our,
        )
    )
    # native isoform, surface-positive, parent in library, our catalog doesn't already count it surface
    if screen_surface and not x["sp_ha"] and g in lib_genes and our != "surface":
        rescue_iso[x["Isoform Transcript"]] = dict(
            gene=g,
            enst=enst,
            uniprot_isoform=iso_id,
            screen_surf=round(float(x["Surface Expression"]), 1),
            in_our_catalog=our is not None,
            our_topology=our,
        )

pd.DataFrame(conc_rows).to_csv(os.path.join(HERE, "willow_isoform_concordance.csv"), index=False)
json.dump(rescue_iso, open(os.path.join(HERE, "willow_screen_rescue_isoforms.json"), "w"), indent=1)

print("\n===== ISOFORM concordance (via ENST->UniProt-isoform map) =====")
print(
    f"screen isoforms: {len(iso)} | mapped to a UniProt isoform id: {iso['enst'].map(lambda e: bool((emap.get(e) or {}).get('isoform_id'))).sum()}"
)
print(f"present in OUR catalog (evaluable): {matched}")
print(f"  of the evaluable, screen-surface & our-topology-surface agree: {both_surface}")
print(f"  screen-surface but our topology says NON-surface (disagreement): {topo_disagree}")
print(f"native (non-SP-HA) surface-positive isoforms to ADD for library genes: {len(rescue_iso)}")
print(
    f"screen isoforms that are Ensembl-only (no UniProt isoform, can't be in our catalog): {(~iso['enst'].isin([e for e, v in emap.items() if v.get('isoform_id')])).sum()}"
)
