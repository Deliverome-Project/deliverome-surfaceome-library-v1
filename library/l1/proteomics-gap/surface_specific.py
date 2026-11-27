"""Surface-specific recurrence, S3E8 sweep-v11, against the FG library (L1/L2).

Counts publications per protein in three non-overlapping kinds of S3E8 input:
  surface  = admitted intact-cell surface-labeling inputs (core + provisional; they vote)
  membrane = membrane-fractionation inputs (complementary; never vote)
  global   = post-lysis whole-lysate glycoproteomes (complementary; the only whole-lysate data in S3E8)
Rule (Becca, 2026-11-15): surface >= 20, surface - membrane >= 10, surface - global >= 10.
Also reports per-kind detection fractions, because the three kinds have different numbers of
publications. Library join: Ensembl gene id via gene_resolver; accession -> symbol from S3E8.
"""
import sys, csv
from collections import defaultdict
import pandas as pd
from deliverome_analysis.cloning import library_design as ld
from deliverome_analysis.surfaceome import gene_resolver, surfaceome
D, S = sys.argv[1], sys.argv[2]
man = pd.read_csv(f"{D}/manifest.tsv", sep="\t", dtype=str)
man = man[man.species.eq("human") & man.retired_reason.fillna("").eq("")]
man["sym"] = man.key.map(ld.gene_from_key)
kinds = {"surface": defaultdict(set), "membrane": defaultdict(set), "global": defaultdict(set)}
pubs = {k: set() for k in kinds}
symbol, raw_symbol = {}, {}
for r in csv.DictReader(open(f"{S}/detections_long.tsv"), delimiter="\t"):
    if r["method_group"] == "membrane-fractionation":
        k = "membrane"
    elif r["review_role"] == "complementary":
        k = "global"
    elif r["admitted"] == "True":
        k = "surface"
    else:
        continue  # held / below-screen surface inputs: neither counted nor penalized
    kinds[k][r["uniprot_accession"]].add(r["publication_id"]); pubs[k].add(r["publication_id"])
    raw_symbol.setdefault(r["uniprot_accession"], r["gene_symbol"])
# Resolved symbols first: detections_long carries the source table's own label, often a
# protein description ("CD151 antigen"), which would miss the library join.
for f in ("protein_method_matrix.tsv", "novel_candidates.tsv", "clustering/candidate_review.tsv"):
    for r in csv.DictReader(open(f"{S}/{f}"), delimiter="\t"):
        if r["gene_symbol"]: symbol.setdefault(r["uniprot_accession"], r["gene_symbol"])
print("publications per kind:", {k: len(v) for k, v in pubs.items()})
cr = {r["uniprot_accession"]: r for r in csv.DictReader(open(f"{S}/clustering/candidate_review.tsv"), delimiter="\t")}
canon, surf = surfaceome.canonical_accessions(), surfaceome.surface_accessions()
ens = gene_resolver.to_ensembl(sorted(set(man.sym) | set(symbol.values())), cache_dir=f"{D}/mygene")
lib = {}
for lb, s in zip(man.library, man.sym):
    for g in ens.get(s) or []: lib.setdefault(g, set()).add(lb)
rows = []
for a, sp in kinds["surface"].items():
    s = symbol.get(a, "")
    resolved = bool(s)
    n_s, n_m, n_g = len(sp), len(kinds["membrane"][a]), len(kinds["global"][a])
    c = cr.get(a, {})
    rows.append({"uniprot_accession": a, "gene_symbol": s or raw_symbol.get(a, ""), "symbol_resolved": resolved,
        "catalog_tier": "canonical" if a in canon else "likely" if a in surf else (c.get("tier") or "not in catalog"),
        "n_surface_pubs": n_s, "n_membrane_pubs": n_m, "n_global_pubs": n_g,
        "surface_minus_membrane": n_s - n_m, "surface_minus_global": n_s - n_g,
        "frac_surface": round(n_s / len(pubs["surface"]), 3), "frac_membrane": round(n_m / len(pubs["membrane"]), 3),
        "frac_global": round(n_g / len(pubs["global"]), 3),
        "n_tm": c.get("n_tm", ""), "has_gpi": c.get("has_gpi", ""), "contaminant_flag": c.get("contaminant_flag", ""),
        "protein_name": c.get("protein_name", ""),
        "library": ",".join(sorted(set().union(*[lib.get(g, set()) for g in (ens.get(s) or [])]))) if s else "?symbol"})
df = pd.DataFrame(rows)
df.to_csv(f"{D}/surface_counts_all_proteins.tsv", sep="\t", index=False)
hit = df[(df.n_surface_pubs >= 20) & (df.surface_minus_membrane >= 10) & (df.surface_minus_global >= 10)]
hit = hit.sort_values("n_surface_pubs", ascending=False)
print("pass rule:", len(hit), "| in library:", hit.library.str.match(r"l[12]").sum(),
      "| not in library:", hit.library.eq("").sum(), "| unknown symbol:", hit.library.eq("?symbol").sum())
print("not-in-library by tier:", hit[~hit.library.str.match(r"l[12]")].catalog_tier.value_counts().to_dict())
hit.to_csv(f"{D}/surface_specific_all.tsv", sep="\t", index=False)
out = hit[~hit.library.str.match(r"l[12]")]
out.to_csv(f"{D}/surface_specific_not_in_library.tsv", sep="\t", index=False)
pd.set_option("display.width", 250)
print(out[["gene_symbol","catalog_tier","n_surface_pubs","n_membrane_pubs","n_global_pubs","n_tm","has_gpi","contaminant_flag","protein_name"]].to_string(index=False, max_colwidth=42))
for g in ("FAM234A",):
    print(df[df.gene_symbol.eq(g)][["gene_symbol","catalog_tier","n_surface_pubs","n_membrane_pubs","n_global_pubs","library"]].to_string(index=False))
