#!/usr/bin/env python3
"""Build the surface-protein library shortlist from the deep-dive cohort.

Reads cohort.jsonl (from fetch_cohort.py), applies the PR130-pinned membership
presets + the displayable red-lines/flags/tiers (Component 5A), the isoform
topology filter (5B) and the ortholog topology filter (5C), and writes the
ranked TSVs + a funnel figure.

Run:  uv run python build_library.py
"""
# ruff: noqa: RUF001  (unicode in docstrings; exploratory scratch script)

import json
import os

import matplotlib.pyplot as plt
import pandas as pd

# Membership now comes from the API-served `deep_dive_tier` tag (deepDiveTier() on
# accessible-surfaceome main, where PR137 is merged). We read the tag directly rather
# than recompute the preset predicate; the preset mirror is kept only as a fallback for
# older cohort snapshots fetched before the API served the tag.
from presets_pr137 import passes_canonical as _preset_canonical
from presets_pr137 import passes_likely as _preset_likely

from deliverome_analysis.plotting import (
    CATEGORICAL_PALETTE,
    save_figure,
    setup_plotting_style,
)


def passes_canonical(f):
    """Canonical membership. Prefer the API's `deep_dive_tier`; fall back to the
    local PR137 preset mirror when the tag is absent (pre-tag cohort snapshots)."""
    t = f.get("deep_dive_tier")
    return t == "canonical" if t is not None else _preset_canonical(f)


def passes_likely(f):
    """Full Likely tier (canonical ⊂ likely) — API tag "canonical" or "likely",
    else the preset mirror. deep_dive_tier=="likely" is the exclusive Likely band."""
    t = f.get("deep_dive_tier")
    return t in ("canonical", "likely") if t is not None else _preset_likely(f)


HERE = os.path.dirname(os.path.abspath(__file__))
COHORT = os.path.join(HERE, "cohort.jsonl")

# ---- vocab ----------------------------------------------------------------
NONSURFACE_REASONS = {
    "stable_complex_partner",
    "secreted_only",
    "inner_leaflet_anchored",
    "cytoplasmic",
    "nuclear",
    "mitochondrial_internal",
    "endomembrane_resident",
    "nuclear_envelope",
    "pmhc_only_intracellular",
}
INDUCED_REASONS = {
    "cell_state_induced",
    "lysosomal_exocytosis",
    "dual_localization",
    "tissue_restricted_surface",
}
LOW_ARCH = {"single_pass_T1", "single_pass_T2", "GPI_anchored", "tetraspanin"}


# ---- flags ----------------------------------------------------------------
def is_olf_taste(uniprot_family, groups):
    txt = (uniprot_family or "").lower() + " " + " ".join(groups or []).lower()
    return any(k in txt for k in ("olfactory", "taste", "vomeronasal"))


def arch_risk(f, ct, uniprot_family, groups):
    """Topology-first (uses real canonical_topology.tm_helix_count).

    Acceptance floor = GPI-anchored OR >=1 TM helix (Becca). Buckets:
      low  = GPI or single-pass (1 TMD)          -> T1
      med  = >=2 TMD (multi-pass, GPCR, tetraspanin) -> T2
      olf_taste = olfactory/taste GPCR (RTP-chaperone-dependent) -> T3
      non_anchored = neither GPI nor any TMD (won't attach alone) -> T3
    """
    sub = f.get("subcategory")
    tm = (ct or {}).get("tm_helix_count") or 0
    if sub == "GPCR" and is_olf_taste(uniprot_family, groups):
        return "olf_taste"
    if sub != "GPI_anchored" and tm < 1:
        return "non_anchored"
    if sub == "GPI_anchored" or tm == 1:
        return "low"
    return "med"


def conditionality(f):
    # Constitutive = has a surface baseline (not induced off a low/absent one),
    # and the presentation mechanism itself isn't induction/regulated exocytosis.
    # Deliberately NOT gated on induction_trigger — a gene can be *further*
    # inducible (tumor/immune) yet still present constitutively (ICAM1-class),
    # which PR130 keeps via low_endogenous_expression=False.
    ok = f.get("low_endogenous_expression") is False and f.get("surface_call_reason") not in {
        "cell_state_induced",
        "lysosomal_exocytosis",
    }
    return "constitutive" if ok else "conditional"


def retention(f):
    # Surface-dominant = predominantly at the surface. A shed/secreted form can
    # still coexist (it's a fraction), so those are carried as separate flags,
    # not folded into retention.
    return "retained" if f.get("surface_specificity") == "surface_dominant" else "leaky"


def expr_prior(f):
    lvl = {"high": 2, "moderate": 1, "low": 0, "absent": 0}.get(f.get("expression_level"), 0)
    brd = 1 if f.get("expression_breadth") in {"pan_tissue", "broad"} else 0
    return lvl + brd


def membership(f):
    if passes_canonical(f):
        return "canonical"
    if passes_likely(f):
        return "likely_only"
    return None


def red_lines_pass(f):
    if f.get("co_receptor_dependency") == "required":
        return False
    if f.get("surface_call_reason") in NONSURFACE_REASONS:
        return False
    return True


def assign_tier(mem, f, arch, cond):
    # Core acceptance (shared by T1/T2): canonical membership, constitutive,
    # surface_specificity in {surface_dominant, mixed} (already guaranteed by
    # passes_canonical), co_receptor_dependency in {none, modulatory}.
    # T1 vs T2 differ only on architecture (low vs med).
    core = (
        mem == "canonical"
        and cond == "constitutive"
        and f.get("surface_specificity") in {"surface_dominant", "mixed"}
        and f.get("co_receptor_dependency") in {"none", "modulatory"}
    )
    if core and arch == "low":
        return "T1"
    if core and arch == "med":
        return "T2"
    return "T3"


# ---- topology -------------------------------------------------------------
def topology_status(tm_count, sp_len, topo, parent_subcat=None, deeptmhmm_label=None):
    if deeptmhmm_label:
        lab = deeptmhmm_label.upper()
        if lab in {"TM", "SP+TM", "BETA"}:
            return "surface"
        if lab == "SP":
            return "secreted"
        if lab == "GLOB":
            return "intracellular"
    if parent_subcat == "GPI_anchored":
        return "surface"
    topo = topo or ""
    has_mem = ("M" in topo) or ("B" in topo) or ((tm_count or 0) >= 1)
    if has_mem:
        return "surface"
    if ("O" in topo) or ((sp_len or 0) > 0):
        return "secreted"
    return "intracellular"


# ---- main -----------------------------------------------------------------
def main():
    records = [json.loads(line) for line in open(COHORT) if line.strip()]
    print(f"loaded {len(records)} cohort records")

    # --- validation: reproduce PR130's recast counts from live data ---
    raw_canonical = sum(1 for r in records if passes_canonical(r["filters"]))
    raw_likely = sum(1 for r in records if passes_likely(r["filters"]))
    print(
        f"VALIDATION raw passes_canonical={raw_canonical} (PR130 was 2243; PR137 differs) | "
        f"raw passes_likely={raw_likely}"
    )

    disp_rows, iso_rows, orth_rows = [], [], []
    for r in records:
        f = r["filters"]
        mem = membership(f)
        if mem is None or not red_lines_pass(f):
            continue
        uf, grp = r.get("uniprot_family"), r.get("hgnc_gene_groups")
        arch = arch_risk(f, r.get("canonical_topology") or {}, uf, grp)
        cond = conditionality(f)
        ret = retention(f)
        cof = f.get("co_receptor_dependency")
        tier = assign_tier(mem, f, arch, cond)
        disp_rows.append(
            {
                "gene": r["gene"],
                "uniprot": r["uniprot"],
                "tier": tier,
                "membership": mem,
                "arch_risk": arch,
                "conditionality": cond,
                "retention": ret,
                "co_receptor_dependency": cof,
                "subcategory": f.get("subcategory"),
                "surface_call_reason": f.get("surface_call_reason"),
                "state_dependence": f.get("state_dependence"),
                "expr_prior": expr_prior(f),
                "has_shed_form": f.get("has_shed_form"),
                "has_secreted_form": f.get("has_secreted_form"),
                "expression_level": f.get("expression_level"),
                "expression_breadth": f.get("expression_breadth"),
                "evidence_density": f.get("evidence_density"),
                "ecd_accessibility_class": f.get("ecd_accessibility_class"),
                "confidence": f.get("confidence"),
                "evidence_grade": f.get("evidence_grade"),
            }
        )

        ctm = (r.get("canonical_topology") or {}).get("tm_helix_count")
        parent_sub = f.get("subcategory")
        # isoforms (5B)
        for iso in r.get("isoforms", []):
            st = topology_status(
                iso.get("tm_helix_count"),
                iso.get("signal_peptide_length"),
                iso.get("per_residue_topology"),
                parent_sub,
            )
            tmd = (iso.get("tm_helix_count") or 0) - (ctm or 0)
            if st != "surface":
                rec = "drop"
            elif tmd != 0:
                rec = "flag"
            else:
                rec = "include"
            iso_rows.append(
                {
                    "gene": r["gene"],
                    "isoform_id": iso.get("isoform_id"),
                    "topology_status": st,
                    "tm_delta": tmd,
                    "ecd_pct_id": iso.get("ecd_pct_identity_to_canonical"),
                    "full_length_pct_id": iso.get("full_length_pct_identity_to_canonical"),
                    "recommendation": rec,
                    "tier": tier,
                }
            )
        # orthologs (5C)
        for species in ("mouse", "cyno"):
            for o in r.get("orthologs", {}).get(species, []):
                st = topology_status(
                    o.get("tm_helix_count"),
                    None,
                    o.get("per_residue_topology"),
                    parent_sub,
                    o.get("deeptmhmm_label"),
                )
                ecd_id = o.get("ecd_pct_identity_to_human_canonical")
                tmd = (o.get("tm_helix_count") or 0) - (ctm or 0)
                diverges = (
                    (ecd_id is not None and ecd_id < 80)
                    or o.get("tm_absent_from_model")
                    or tmd != 0
                )
                if st != "surface":
                    rec = "drop"
                elif diverges:
                    rec = "evaluate"
                else:
                    rec = "surrogate"
                orth_rows.append(
                    {
                        "gene": r["gene"],
                        "species": species,
                        "ortholog_symbol": o.get("ortholog_symbol"),
                        "topology_status": st,
                        "tm_delta": tmd,
                        "ecd_pct_id": ecd_id,
                        "full_length_pct_id": o.get("full_length_pct_identity_to_human_canonical"),
                        "tm_absent_from_model": o.get("tm_absent_from_model"),
                        "recommendation": rec,
                        "tier": tier,
                    }
                )

    disp = pd.DataFrame(disp_rows).sort_values(
        ["tier", "expr_prior", "gene"], ascending=[True, False, True]
    )
    iso = pd.DataFrame(iso_rows)
    orth = pd.DataFrame(orth_rows)
    disp.to_csv(os.path.join(HERE, "displayable_tiers.tsv"), sep="\t", index=False)
    iso.to_csv(os.path.join(HERE, "isoform_candidates.tsv"), sep="\t", index=False)
    orth.to_csv(os.path.join(HERE, "ortholog_candidates.tsv"), sep="\t", index=False)

    # ---- counts ----
    n_cohort = len(records)
    n_union = len(disp)
    tier_counts = disp["tier"].value_counts().to_dict()
    mem_counts = disp["membership"].value_counts().to_dict()
    lines = [
        "# Surface library — first-pass counts\n",
        f"- cohort (deep-dive): **{n_cohort}**",
        f"- displayable union (canonical ∪ likely, red-lines passed): **{n_union}**",
        f"  - by membership: canonical={mem_counts.get('canonical', 0)}, "
        f"likely_only={mem_counts.get('likely_only', 0)}",
        f"- **T1** (core will-present): **{tier_counts.get('T1', 0)}**",
        f"- **T2** (extended): **{tier_counts.get('T2', 0)}**",
        f"- **T3** (full flagged): **{tier_counts.get('T3', 0)}**",
        "",
        "## arch_risk × tier",
        disp.pivot_table(
            index="arch_risk", columns="tier", values="gene", aggfunc="count", fill_value=0
        ).to_string(),
        "",
        "## Isoform candidates (topology)",
        (iso["topology_status"].value_counts().to_string() if len(iso) else "none"),
        f"isoform recommendations: {iso['recommendation'].value_counts().to_dict() if len(iso) else {}}",
        "",
        "## Ortholog candidates (topology)",
        f"ortholog recommendations: {orth['recommendation'].value_counts().to_dict() if len(orth) else {}}",
    ]
    open(os.path.join(HERE, "counts.md"), "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))

    # ---- funnel figure ----
    setup_plotting_style()
    stages = ["cohort", "union", "T1+T2+T3", "T2+T1", "T1"]
    vals = [
        n_cohort,
        n_union,
        tier_counts.get("T1", 0) + tier_counts.get("T2", 0) + tier_counts.get("T3", 0),
        tier_counts.get("T1", 0) + tier_counts.get("T2", 0),
        tier_counts.get("T1", 0),
    ]
    fig, ax = plt.subplots(figsize=(7, 4))
    colors = (CATEGORICAL_PALETTE * 3)[: len(stages)]
    ax.barh(range(len(stages)), vals, color=colors)
    ax.set_yticks(range(len(stages)))
    ax.set_yticklabels(stages)
    ax.invert_yaxis()
    for i, v in enumerate(vals):
        ax.text(v, i, f" {v}", va="center")
    ax.set_xlabel("genes")
    save_figure(fig, "surface_library_funnel", output_dir=os.path.join(HERE, "figures"))
    print(
        "wrote displayable_tiers.tsv, isoform_candidates.tsv, ortholog_candidates.tsv, counts.md, funnel"
    )


if __name__ == "__main__":
    main()
