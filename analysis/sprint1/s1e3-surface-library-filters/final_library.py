#!/usr/bin/env python3
"""Final surface-library cut: apply the full filter set, tier by architecture +
presentation, and count genes / surface-competent isoforms / surface-competent
orthologs per tier. Run: uv run python final_library.py
"""
# ruff: noqa: E402, E741, RUF003  (exploratory scratch script; ∪ used in comments)

import csv
import json
import os

import build_library as B  # reuse arch_risk, topology_status, is_olf_taste, NONSURFACE_REASONS, membership
from coreceptor_rules import integral_complex_ok, reviewed_partner_ok

HERE = os.path.dirname(os.path.abspath(__file__))
recs = [json.loads(l) for l in open(os.path.join(HERE, "cohort.jsonl"))]

# PR137 gates on effective_evidence_grade = evidence_grade_summary ?? evidence_grade.
# The refreshed cohort carries evidence_grade_summary in filters directly; for older
# snapshots that lack it, fall back to the per-gene fetch (fetch_evidence_summary.py).
_egs_path = os.path.join(HERE, "evidence_grade_summary.json")
if os.path.exists(_egs_path):
    _egs = json.load(open(_egs_path))
    for r in recs:
        if r["filters"].get("evidence_grade_summary"):
            continue  # already present from the refreshed cohort
        s = _egs.get(r["gene"])
        if s and s != "__ERROR__":
            r["filters"]["evidence_grade_summary"] = s
expr = json.load(open(os.path.join(HERE, "external/expression_tpm/resolved_expr.json")))
# Full re-derived expression (NO >=10 pre-filter) — real max TPM incl. low expressers,
# so the expression gate can use a lower cutoff without lacking data for <10-TPM genes.
mxfull = json.load(open(os.path.join(HERE, "external/expression_tpm/resolved_expr_full.json")))
cohort_genes = {r["gene"] for r in recs}
cohort_accessions = {record["gene"]: record.get("uniprot") for record in recs}
# Co-receptor partners come straight from the cohort record — the API field
# accessibility_risks.co_receptor_requirements.partners, captured by fetch_cohort's
# compact(). Falls back to the legacy coreceptor_partners.jsonl side-file only for
# older cohort snapshots that predate the captured field.
part = {
    r["gene"]: {"partners": r["co_receptor_partners"]}
    for r in recs
    if r.get("co_receptor_partners")
}
if not part and os.path.exists(os.path.join(HERE, "coreceptor_partners.jsonl")):
    part = {
        o["gene"]: o
        for o in (json.loads(l) for l in open(os.path.join(HERE, "coreceptor_partners.jsonl")))
        if "partners" in o
    }

# Co-receptor INCLUDE_ALONE allowlist: subunits whose SAVED deep-dive rationale shows
# they reach the surface without a heterologous partner (chaperone/scaffold/homodimer).
# Source: coreceptor_calls.tsv (per-gene calls, uploaded to processed s1e3). Gate is
# otherwise unchanged — genuine obligate pairs stay excluded (no co-delivery).
CORECEPTOR_ALLOWLIST = set()
_cc = os.path.join(HERE, "coreceptor_calls.tsv")
if os.path.exists(_cc):
    for _row in csv.DictReader(open(_cc), delimiter="\t"):
        if _row.get("call") == "INCLUDE_ALONE":
            CORECEPTOR_ALLOWLIST.add(_row["gene"])

# mostly_intracellular compartment rescue: admit a mostly_intracellular gene if its
# PRIMARY subcellular compartment is endosome/lysosome/PM (cycles to the surface —
# FcRn, GLUT4, RANKL) but NOT ER/Golgi (CANX, GRP94, FURIN — genuine residents).
# Prefer the refreshed cohort's primary_compartment; fall back to the older side file.
_comp = {r["gene"]: r.get("primary_compartment") for r in recs if r.get("primary_compartment")}
if not _comp and os.path.exists(os.path.join(HERE, "mostly_intracellular_compartments.json")):
    _comp = json.load(open(os.path.join(HERE, "mostly_intracellular_compartments.json")))
_ER_GOLGI = ("er", "endoplasmic", "golgi", "ergic", "reticulum")
_SURF_CYCLE = ("endosome", "lysosome", "plasma", "cell membrane", "recycl", "vesicle")


def compartment_surface(g):
    c = (_comp.get(g) or "").lower()
    if any(k in c for k in _ER_GOLGI):
        return False
    return any(k in c for k in _SURF_CYCLE)


# ECD size gate: drop genes whose extracellular domain is < MIN_ECD_RES residues
# ('minimal' category — no displayable/targetable ectodomain: FXYD, IFITM, CMTM,
# TRAT1, tail-anchored). GPI-anchored proteins are exempt (mature protein is the ECD).
MIN_ECD_RES = 30


def ecd_residues(r):
    return ((r.get("canonical_topology") or {}).get("per_residue_topology") or "").count("O")


def ecd_ok(r):
    if r["filters"].get("subcategory") == "GPI_anchored":
        return True
    return ecd_residues(r) >= MIN_ECD_RES


# partner median cell-line TPM from the re-derived per-gene cell-line summary
med_ens, med_sym = {}, {}
for r in csv.DictReader(
    open(os.path.join(HERE, "external/hpa_rna_celline/hpa_celline_per_gene.csv"))
):
    try:
        m = float(r["median_tpm_celline"])
    except ValueError:
        continue
    med_ens[r["ensembl_base"]] = m
    med_sym[r["gene_name"]] = m


def partner_median(p):
    for e in (expr.get(p) or {}).get("ensembl", []):
        if e.split(".")[0] in med_ens:
            return med_ens[e.split(".")[0]]
    return med_sym.get(p)


def reviewed_partner_median(partner):
    identifiers = {
        identifier.split(".")[0]
        for identifier in (expr.get(partner) or {}).get("ensembl", [])
        if identifier.split(".")[0] in med_ens
    }
    return med_ens.get(next(iter(identifiers))) if len(identifiers) == 1 else None


def mx(g):
    return mxfull.get(g)  # full re-derived reference (real value, incl. <10 TPM)


def _partner_sym(p):
    p = str(p)
    return p.replace("(", " ").split()[0] if p.split() else p


def _is_self_partner(g, p):
    # A homo-oligomer / self partner: the gene makes its own partner, so a
    # required "co-receptor" that is itself is always satisfied (fixes BDKRB1,
    # OTOP1/2/3, PCDHGB7 — never exclude a homodimer on co-receptor grounds).
    return _partner_sym(p) == g or "homo" in str(p).lower() or "self" in str(p).lower()


CODELIVER = False  # co-deliverability rejected: a partner that isn't endogenously expressed in
# common cell lines still won't chaperone the gene to the surface, so keep the ubiquity rule.


def coreceptor_ok(g, f):
    reviewed = reviewed_partner_ok(
        cohort_accessions.get(g), part.get(g, {}).get("partners", []), reviewed_partner_median
    )
    if reviewed is not None:
        return reviewed
    if g in CORECEPTOR_ALLOWLIST:  # rationale shows it traffics without the partner
        return True
    if f.get("co_receptor_dependency") != "required":
        return True
    o = part.get(g)
    if not o or not o["partners"]:
        return False
    non_self = [p for p in o["partners"] if not _is_self_partner(g, p)]
    if not non_self:  # all partners are self/homo-oligomer -> satisfied
        return True
    if CODELIVER and all(_partner_sym(p) in cohort_genes for p in non_self):
        return True  # co-deliverable: every required partner is itself a library gene
    meds = [m for m in (partner_median(p) for p in non_self) if m is not None]
    return bool(meds) and all(m >= 5 for m in meds)  # ubiquitous = median cell-line TPM >= 5


def anchored(r):
    f = r["filters"]
    return (
        f.get("subcategory") == "GPI_anchored"
        or ((r.get("canonical_topology") or {}).get("tm_helix_count") or 0) >= 1
    )


def localization_ok(record):
    reason = record["filters"].get("surface_call_reason")
    if reason == "stable_complex_partner":
        return integral_complex_ok(
            record, part.get(record["gene"], {}).get("partners", []), reviewed_partner_median
        )
    return reason not in B.NONSURFACE_REASONS


EXPR_MIN_TPM = 5  # relaxed from 10; full reference has real low-TPM values now


def expr_ok(g, f):
    return f.get("low_endogenous_expression") is False or (
        mx(g) is not None and mx(g) >= EXPR_MIN_TPM
    )


# Lentiviral packaging: LTR-LTR budget ~8.5 kb. Current p003 construct -> ORFs up
# to ~4 kb; with mScarlet removed ~4.7 kb. Provisional cap 5 kb. ORF nt = aa*3 + stop.
MAX_ORF_NT = 5000


def orf_nt(r):
    t = (r.get("canonical_topology") or {}).get("per_residue_topology")
    return len(t) * 3 + 3 if t else None


def orf_ok(r):
    n = orf_nt(r)
    return n is None or n <= MAX_ORF_NT


def passes(r):
    g, f = r["gene"], r["filters"]
    return (
        # Membership base: library = (canonical ∪ likely) with the surface filters
        # on top (PR137 gates via presets_pr137, evaluated on the injected
        # evidence_grade_summary). Genes that clear the surface filters but are
        # neither canonical nor likely (mostly `surface_accessibility=low`) are
        # NOT included — they fall to the below-Likely tier.
        (B.passes_canonical(f) or B.passes_likely(f))
        and f.get("surface_accessibility") != "no"
        # compartment rescue dropped: mostly_intracellular genes stay out (OXTR/ACKR2
        # fall back to the SV screen tier; the rest were mostly non-targetable).
        and f.get("surface_specificity") in ("surface_dominant", "mixed")
        and anchored(r)
        and coreceptor_ok(g, f)
        and expr_ok(g, f)
        and localization_ok(r)
        and orf_ok(r)
        and ecd_ok(r)  # drop minimal-ECD (<30 extracellular residues) — no displayable domain
    )


def tier(r):
    f = r["filters"]
    arch = B.arch_risk(
        f, r.get("canonical_topology") or {}, r.get("uniprot_family"), r.get("hgnc_gene_groups")
    )
    constitutive = f.get("surface_call_reason") not in (
        "cell_state_induced",
        "lysosomal_exocytosis",
    )
    if arch == "low" and constitutive:
        return "T1"
    if arch == "med" and constitutive:
        return "T2"
    return "T3"  # high-risk arch (olfactory/taste) or conditional presentation


def n_surface_isoforms(r):
    sub = r["filters"].get("subcategory")
    return sum(
        1
        for iso in r.get("isoforms", [])
        if B.topology_status(
            iso.get("tm_helix_count"),
            iso.get("signal_peptide_length"),
            iso.get("per_residue_topology"),
            sub,
        )
        == "surface"
    )


# Orthologs are cross-species surrogates for preclinical binder testing; a binder against
# the human ECD only cross-reacts if the ortholog ECD is conserved. Drop orthologs whose
# ECD % identity to human canonical is below this floor (therapeutically irrelevant surrogates).
ORTHO_ECD_MIN = 70


def n_surface_orthologs(r):
    sub = r["filters"].get("subcategory")
    n = 0
    for sp in ("mouse", "cyno"):
        for o in r.get("orthologs", {}).get(sp, []):
            if (o.get("ecd_pct_identity_to_human_canonical") or 0) < ORTHO_ECD_MIN:
                continue
            if (
                B.topology_status(
                    o.get("tm_helix_count"),
                    None,
                    o.get("per_residue_topology"),
                    sub,
                    o.get("deeptmhmm_label"),
                )
                == "surface"
            ):
                n += 1
    return n


passing = [r for r in recs if passes(r)]
from collections import defaultdict

agg = defaultdict(lambda: {"genes": 0, "iso": 0, "iso_genes": 0, "orth": 0, "orth_genes": 0})
_lib_tier = {}  # gene -> tier, for the exported library membership list
for r in passing:
    t = tier(r)
    _lib_tier[r["gene"]] = t
    a = agg[t]
    a["genes"] += 1
    ni, no = n_surface_isoforms(r), n_surface_orthologs(r)
    a["iso"] += ni
    a["orth"] += no
    a["iso_genes"] += ni > 0
    a["orth_genes"] += no > 0

# Screen-validated GPCR set, loaded EARLY so the LL tier can defer to it: SV's direct
# screen evidence supersedes the low-lit UniProt proxy when a gene qualifies for both
# (previously such genes were counted in both LL and SV, inflating the total).
_by_gene = {r["gene"]: r for r in recs}
_sv = json.load(open(os.path.join(HERE, "willow_screen_rescue_genes.json")))
_sv_iso = json.load(open(os.path.join(HERE, "willow_screen_rescue_isoforms.json")))

# --- Low-lit / UniProt-positive deterministic additions (tier LL; TPM kept, no family exclusion) ---
# Non-canonical, UniProt-surface (db bit 0), n_papers_found < 100, that our LLM gates
# excluded but pass the deterministic filters (membrane-anchored + max TPM >= 10).
_rows = json.load(open(os.path.join(HERE, "catalog_rows.json")))
_uni = {r["symbol"]: bool(r["db"] & 1) for r in _rows}
_current = {r["gene"] for r in passing}
for r in recs:
    g, f = r["gene"], r["filters"]
    if g in _current or g in _sv or not _uni.get(g):
        continue  # skip core passers and screen-validated genes (SV supersedes LL)
    n = f.get("n_papers_found")
    if n is None or n >= 100 or B.passes_canonical(f):
        continue
    if not anchored(r) or (mx(g) or 0) < EXPR_MIN_TPM or not orf_ok(r):
        continue
    _lib_tier[g] = "LL"
    a = agg["LL"]
    a["genes"] += 1
    ni, no = n_surface_isoforms(r), n_surface_orthologs(r)
    a["iso"] += ni
    a["orth"] += no
    a["iso_genes"] += ni > 0
    a["orth_genes"] += no > 0

# --- Screen-validated surface GPCRs (tier SV) from Willow's screen (Tedman 2025) ---
# GPCRs our filters excluded but that show clear cell-surface expression in the screen
# (immunostaining >= median). Direct experimental evidence of surface competence when
# expressed, so included; flagged in a distinct tier, not merged into T1-T3. Also carries
# the count of native screen-validated isoform variants (Ensembl-only, absent from our
# UniProt-isoform catalog). Both lists are built by willow_rescue.py.
_passing_genes = {r["gene"] for r in passing}
for g in _sv:
    r = _by_gene.get(g)
    if r is None or g in _passing_genes:
        continue  # already in core (compartment rescue / co-receptor allowlist) — don't double-count
    _lib_tier[g] = "SV"
    a = agg["SV"]
    a["genes"] += 1
    ni, no = n_surface_isoforms(r), n_surface_orthologs(r)
    a["iso"] += ni
    a["orth"] += no
    a["iso_genes"] += ni > 0
    a["orth_genes"] += no > 0

# --- Overexpression-surface rescued (tier OE): mostly_intracellular genes with a GENUINE
# functional PM pool on overexpression — the DLL3 / GLUT4 / FcRn pattern. Subagent-classified
# RESCUE from the saved overexpression_surface_localization_observed rationale + surface-
# specificity rationale (oe_calls_136.tsv). Blocked by the surface_specificity gate alone;
# included as a flagged tier, not merged into T1-T3. ec_sites in oe_calls_136.tsv marks which
# carry an extracellular targetable patch (display) vs expression-only (GLUT4/AQP2-type). ---
_oe = [
    row["gene"]
    for row in csv.DictReader(open(os.path.join(HERE, "oe_calls_136.tsv")), delimiter="\t")
    if row["call"] == "RESCUE"
]
for g in _oe:
    r = _by_gene.get(g)
    if r is None or g in _lib_tier:
        continue  # already counted in core / LL / SV
    _lib_tier[g] = "OE"
    a = agg["OE"]
    a["genes"] += 1
    ni, no = n_surface_isoforms(r), n_surface_orthologs(r)
    a["iso"] += ni
    a["orth"] += no
    a["iso_genes"] += ni > 0
    a["orth_genes"] += no > 0


# --- Oversized (>5 kb ORF): clears every gate EXCEPT the ORF cap. NOT in the library
# (handled separately as truncated constructs) but tracked as a distinct bar for reference. ---
def _passes_except_orf(r):
    g, f = r["gene"], r["filters"]
    return (
        (B.passes_canonical(f) or B.passes_likely(f))
        and f.get("surface_accessibility") != "no"
        and f.get("surface_specificity") in ("surface_dominant", "mixed")
        and anchored(r)
        and coreceptor_ok(g, f)
        and expr_ok(g, f)
        and localization_ok(r)
        and ecd_ok(r)
    )


for r in recs:
    if _passes_except_orf(r) and not orf_ok(r):
        a = agg["OV"]
        a["genes"] += 1
        ni, no = n_surface_isoforms(r), n_surface_orthologs(r)
        a["iso"] += ni
        a["orth"] += no
        a["iso_genes"] += ni > 0
        a["orth_genes"] += no > 0

_total = len(passing) + agg["LL"]["genes"] + agg["SV"]["genes"] + agg["OE"]["genes"]
print(
    f"TOTAL: {len(passing)} core + {agg['LL']['genes']} low-lit/UniProt "
    f"+ {agg['SV']['genes']} screen-validated GPCR + {agg['OE']['genes']} overexpr-rescue "
    f"= {_total} genes (+{len(_sv_iso)} screen-validated isoform variants)\n"
)
print(
    f"{'tier':6}{'genes':>8}{'surf_isoforms':>15}{'(genes w/iso)':>15}{'surf_orthologs':>16}{'(genes w/orth)':>16}"
)
tot = {"genes": 0, "iso": 0, "iso_genes": 0, "orth": 0, "orth_genes": 0}
for t in ("T1", "T2", "T3", "LL", "SV", "OE"):
    a = agg[t]
    for k in tot:
        tot[k] += a[k]
    print(
        f"{t:6}{a['genes']:>8}{a['iso']:>15}{a['iso_genes']:>15}{a['orth']:>16}{a['orth_genes']:>16}"
    )
print(
    f"{'ALL':6}{tot['genes']:>8}{tot['iso']:>15}{tot['iso_genes']:>15}{tot['orth']:>16}{tot['orth_genes']:>16}"
)
print(
    f"{'OV':6}{agg['OV']['genes']:>8}{agg['OV']['iso']:>15}{agg['OV']['iso_genes']:>15}{agg['OV']['orth']:>16}{agg['OV']['orth_genes']:>16}  (oversized >5kb — NOT in library, truncation list)"
)
_counts = {t: dict(agg[t]) for t in ("T1", "T2", "T3", "LL", "SV", "OE")}
_counts["ALL"] = dict(tot)
_counts["OV"] = dict(agg["OV"])  # oversized, tracked separately (not summed into ALL)
json.dump(_counts, open(os.path.join(HERE, "tier_counts.json"), "w"), indent=1)
# export the full library membership (gene -> tier) for cross-referencing / downstream use
with open(os.path.join(HERE, "library_membership.tsv"), "w") as _fh:
    _fh.write("gene\ttier\n")
    for _g, _t in sorted(_lib_tier.items()):
        _fh.write(f"{_g}\t{_t}\n")

# membership / rescue provenance among passers
from collections import Counter

mem = Counter()
tpm_rescued = corec_rescued = 0
for r in passing:
    f = r["filters"]
    m = "canonical" if B.passes_canonical(f) else ("likely" if B.passes_likely(f) else "neither")
    mem[m] += 1
    if f.get("low_endogenous_expression") is True and (mx(r["gene"]) or 0) >= 10:
        tpm_rescued += 1
    if f.get("co_receptor_dependency") == "required":
        corec_rescued += 1
print(f"\nmembership of passers: {dict(mem)}")
print(
    f"expression-rescued (low_endog but TPM>=10): {tpm_rescued} | co-receptor-rescued (required+ubiquitous): {corec_rescued}"
)
