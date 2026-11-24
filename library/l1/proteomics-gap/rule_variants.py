"""Rules A and B side by side on surface_specific.py's per-protein counts, with diffs to a baseline.

Run surface_specific.py first (unchanged); it writes <D>/surface_counts_all_proteins.tsv.
  A (Becca, 2026-11-15): surface >= 20, surface - membrane >= 10, surface - global >= 10.
  B (log-odds): surface >= 20 and lod_membrane > 0.5, where
    lod_x = logit((k_s + 0.5) / (N_s + 1)) - logit((k_x + 0.5) / (N_x + 1)),
    k = a protein's publications of a kind, N = all publications of that kind.
  lod_global (against the whole-lysate glycoproteomes) is reported, never filtered on: a
  glycoproteome filter drops surface glycoproteins such as KIT, VCAM1 and DPP4.
N per kind is counted from detections_long.tsv with surface_specific.py's own classification.

    uv run --locked python rule_variants.py <D> <sweep data folder> [<baseline D>]
"""

import csv
import math
import sys

import pandas as pd

D, S = sys.argv[1], sys.argv[2]
BASE = sys.argv[3] if len(sys.argv) > 3 else None
pubs = {"surface": set(), "membrane": set(), "global": set()}
for r in csv.DictReader(open(f"{S}/detections_long.tsv"), delimiter="\t"):
    if r["method_group"] == "membrane-fractionation":
        pubs["membrane"].add(r["publication_id"])
    elif r["review_role"] == "complementary":
        pubs["global"].add(r["publication_id"])
    elif r["admitted"] == "True":
        pubs["surface"].add(r["publication_id"])
N = {k: len(v) for k, v in pubs.items()}


def logit(p: float) -> float:
    return math.log(p / (1 - p))


def lod(k_s: int, k_x: int, n_x: int) -> float:
    return logit((k_s + 0.5) / (N["surface"] + 1)) - logit((k_x + 0.5) / (n_x + 1))


df = pd.read_csv(f"{D}/surface_counts_all_proteins.tsv", sep="\t", keep_default_na=False)
df["lod_membrane"] = [
    round(lod(s, m, N["membrane"]), 3) for s, m in zip(df.n_surface_pubs, df.n_membrane_pubs)
]
df["lod_global"] = [
    round(lod(s, g, N["global"]), 3) for s, g in zip(df.n_surface_pubs, df.n_global_pubs)
]
df["in_library"] = df.library.astype(str).str.match(r"l[12]")
rule = {
    "A": (df.n_surface_pubs >= 20) & (df.surface_minus_membrane >= 10) & (df.surface_minus_global >= 10),
    "B": (df.n_surface_pubs >= 20) & (df.lod_membrane > 0.5),
}
print("publications per kind:", N)
cols = ["uniprot_accession", "gene_symbol", "catalog_tier", "n_surface_pubs", "n_membrane_pubs",
        "n_global_pubs", "lod_membrane", "lod_global", "library"]  # fmt: skip
for name, mask in rule.items():
    hit = df[mask].sort_values("n_surface_pubs", ascending=False)
    hit.to_csv(f"{D}/rule_{name}_pass.tsv", sep="\t", index=False)
    hit[~hit.in_library].to_csv(f"{D}/rule_{name}_not_in_library.tsv", sep="\t", index=False)
    print(f"rule {name}: {len(hit)} pass, {(~hit.in_library).sum()} not in L1/L2")
    if not BASE:
        continue
    old = pd.read_csv(f"{BASE}/rule_{name}_not_in_library.tsv", sep="\t", keep_default_na=False)
    new = hit[~hit.in_library]
    o, n = old.set_index("uniprot_accession"), new.set_index("uniprot_accession")
    rows = []
    for a in sorted(set(o.index) | set(n.index)):
        r = {"uniprot_accession": a, "gene_symbol": (n if a in n.index else o).loc[a, "gene_symbol"]}
        r["change"] = "added" if a not in o.index else "dropped" if a not in n.index else "kept"
        for c in ("n_surface_pubs", "n_membrane_pubs", "n_global_pubs", "lod_membrane"):
            r[f"{c}_old"] = o.loc[a, c] if a in o.index else ""
            r[f"{c}_new"] = n.loc[a, c] if a in n.index else ""
        if a not in n.index:  # why it left: rule now fails, or now in the library
            row = df[df.uniprot_accession.eq(a)]
            r["reason"] = (
                "not detected in this snapshot" if row.empty
                else "now in L1/L2" if row.in_library.iloc[0] else "now fails the rule"
            )
        else:
            r["reason"] = ""
        r["counts_changed"] = r["change"] == "kept" and any(
            str(r[f"{c}_old"]) != str(r[f"{c}_new"])
            for c in ("n_surface_pubs", "n_membrane_pubs", "n_global_pubs")
        )
        rows.append(r)
    diff = pd.DataFrame(rows)
    diff.to_csv(f"{D}/rule_{name}_diff_vs_baseline.tsv", sep="\t", index=False)
    print(
        f"  vs baseline: {len(o)} -> {len(n)} not in library; added {diff.change.eq('added').sum()}, "
        f"dropped {diff.change.eq('dropped').sum()}, counts changed {diff.counts_changed.sum()}"
    )
pd.set_option("display.width", 200)
for g in ("TFRC", "LAMP1", "CD44", "CANX", "HSPA5", "FAM234A", "KIT", "VCAM1", "DPP4"):
    print(df[df.gene_symbol.eq(g)][cols].to_string(index=False, header=False))
