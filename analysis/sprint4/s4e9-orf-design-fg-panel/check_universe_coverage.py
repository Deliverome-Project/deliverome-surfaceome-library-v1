"""Which library genes the isoform and ortholog universes could draw from, and why not the rest.

Both universes were built from accessible-surfaceome files frozen at commit f9c73691
(2026-04-27), whose human query set is deliverome-internal's expression-filtered
"surfaceome_expressed" list (2,359 genes). The library gene list came later, from the S1E3
filter over surfaceome release 1.3.0. Later Compara pulls (2026-06-01, 2026-09-20/21, and
the 2026-10-12 consolidation onto Compara 115) went into the surfaceome D1 tables, not into
the f9c73691 files the library reads. So a library gene outside the April set was never
queried for orthologs, and its isoforms were never offered to the isoform universe.

This measures that, from pinned inputs only:

* the April human query set: `human_canonical_non_hla.fasta` at f9c73691;
* the recorded ortholog lists: `{mouse,cyno}_ortholog_one2one_highconf_non_hla_metadata.csv`
  at f9c73691;
* Ensembl Compara release 115 protein homology dumps (human file for cyno, mouse file for
  mouse), filtered to `ortholog_one2one` with `is_high_confidence == 1`, the rule the April
  BioMart pull applied;
* gene symbol -> Ensembl gene from the Ensembl 115 human peptide FASTA, primary assembly;
* the isoform universe as recorded (`isoform_cohorts.tsv`) and UniProt 2026_03 reviewed
  isoforms.

It asserts every recorded ortholog pair is one-to-one high-confidence in Compara 115, and
writes one row per library gene. Nothing here changes the universes or the panels: the gap
is recorded, not repaired (Becca, 2026-11-05).

    uv run --locked python check_universe_coverage.py
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
import tempfile
from collections import defaultdict
from pathlib import Path

import library_pipeline as lp

HERE = Path(__file__).parent
OUT = HERE / "universe_coverage.tsv"
SUMMARY = HERE / "universe_coverage_summary.json"
INPUTS = [
    "human_canonical_non_hla.fasta",
    "mouse_ortholog_one2one_highconf_non_hla_metadata.csv",
    "cyno_ortholog_one2one_highconf_non_hla_metadata.csv",
    "Compara.115.homo_sapiens.protein_default.homologies.tsv.gz",
    "Compara.115.mus_musculus.protein_default.homologies.tsv.gz",
    "Homo_sapiens.GRCh38.pep.all.fa.gz",
    "uniprot_human_reviewed_isoforms.tsv.gz",
    "isoform_cohorts.tsv",
]


def one2one_hc(path: Path, other: str, human_side: str) -> dict[str, set[str]]:
    """Human gene -> its one-to-one high-confidence ortholog gene(s) in ``other``."""
    pairs: dict[str, set[str]] = defaultdict(set)
    with gzip.open(path, "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["homology_type"] != "ortholog_one2one" or r["is_high_confidence"] != "1":
                continue
            if human_side == "query" and r["homology_species"] == other:
                pairs[r["gene_stable_id"]].add(r["homology_gene_stable_id"])
            elif human_side == "homology" and r["homology_species"] == "homo_sapiens":
                pairs[r["homology_gene_stable_id"]].add(r["gene_stable_id"])
    return pairs


def symbol_to_genes(pep: Path) -> dict[str, set[str]]:
    out: dict[str, set[str]] = defaultdict(set)
    with gzip.open(pep, "rt") as fh:
        for line in fh:
            if not line.startswith(">") or not line.split()[2].startswith("chromosome:"):
                continue
            tags = dict(t.split(":", 1) for t in line[1:].split()[1:] if ":" in t)
            if tags.get("gene_symbol"):
                out[tags["gene_symbol"]].add(tags["gene"].split(".")[0])
    return out


def measure(paths: dict[str, Path], genes: list[str]) -> tuple[list[dict], dict]:
    april = {
        m.group(1)
        for line in paths["human_canonical_non_hla.fasta"].open()
        if line.startswith(">") and (m := re.search(r"GN=(\S+)", line))
    }
    recorded = {}
    for sp in ("mouse", "cyno"):
        with paths[f"{sp}_ortholog_one2one_highconf_non_hla_metadata.csv"].open() as fh:
            recorded[sp] = list(csv.DictReader(fh))
    compara = {
        "cyno": one2one_hc(
            paths["Compara.115.homo_sapiens.protein_default.homologies.tsv.gz"],
            "macaca_fascicularis",
            "query",
        ),
        "mouse": one2one_hc(
            paths["Compara.115.mus_musculus.protein_default.homologies.tsv.gz"],
            "mus_musculus",
            "homology",
        ),
    }
    not_hc = [
        (sp, r["query_input_gene_symbols"], r["query_ensembl_gene_id"], r["target_ensembl_gene_id"])
        for sp in ("mouse", "cyno")
        for r in recorded[sp]
        if r["target_ensembl_gene_id"] not in compara[sp].get(r["query_ensembl_gene_id"], set())
    ]
    assert not not_hc, f"recorded pairs not one-to-one high-confidence in Compara 115: {not_hc[:5]}"
    listed = {sp: {r["query_ensembl_gene_id"] for r in recorded[sp]} for sp in recorded}
    sym = symbol_to_genes(paths["Homo_sapiens.GRCh38.pep.all.fa.gz"])
    with gzip.open(paths["uniprot_human_reviewed_isoforms.tsv.gz"], "rt") as fh:
        has_isoform = {
            r["Gene Names (primary)"]
            for r in csv.DictReader(fh, delimiter="\t")
            if "-" in r["Entry"] and not r["Entry"].endswith("-1") and r["Gene Names (primary)"]
        }
    with paths["isoform_cohorts.tsv"].open() as fh:
        in_iso_universe = {r["gene"] for r in csv.DictReader(fh, delimiter="\t")}

    rows = []
    for g in sorted(genes):
        ensg = sym.get(g, set())
        row = {
            "gene": g,
            "ensembl_gene": ";".join(sorted(ensg)),
            "in_april_query_set": g in april,
        }
        for sp in ("mouse", "cyno"):
            eligible = any(e in compara[sp] for e in ensg)
            row[f"{sp}_one2one_hc_compara115"] = eligible
            row[f"{sp}_in_ortholog_list"] = bool(ensg & listed[sp])
        row["has_uniprot_isoform"] = g in has_isoform
        row["in_isoform_universe"] = g in in_iso_universe
        rows.append(row)

    def count(pred):
        return sum(1 for r in rows if pred(r))

    summary = {
        "library_genes": len(rows),
        "april_query_set": len(april),
        "library_genes_in_april_set": count(lambda r: r["in_april_query_set"]),
        "library_genes_without_ensembl_gene": count(lambda r: not r["ensembl_gene"]),
        "recorded_pairs_checked": {sp: len(recorded[sp]) for sp in recorded},
    }
    for sp in ("mouse", "cyno"):
        missing = [
            r for r in rows if r[f"{sp}_one2one_hc_compara115"] and not r[f"{sp}_in_ortholog_list"]
        ]
        summary[sp] = {
            "eligible_in_compara115": count(lambda r, s=sp: r[f"{s}_one2one_hc_compara115"]),
            "in_ortholog_list": count(lambda r, s=sp: r[f"{s}_in_ortholog_list"]),
            "eligible_not_listed": len(missing),
            "eligible_not_listed_outside_april_set": sum(
                1 for r in missing if not r["in_april_query_set"]
            ),
        }
    iso_missing = [r for r in rows if r["has_uniprot_isoform"] and not r["in_isoform_universe"]]
    summary["isoform"] = {
        "genes_with_uniprot_isoform": count(lambda r: r["has_uniprot_isoform"]),
        "in_isoform_universe": count(lambda r: r["in_isoform_universe"]),
        "isoform_gene_not_in_universe": len(iso_missing),
        "of_which_outside_april_set": sum(1 for r in iso_missing if not r["in_april_query_set"]),
    }
    return rows, summary


def write(rows: list[dict], summary: dict, out: Path, summary_out: Path) -> None:
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    summary_out.write_text(json.dumps(summary, indent=1) + "\n")


def run(out: Path = OUT, summary_out: Path = SUMMARY) -> dict:
    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        paths = lp.stage_inputs(INPUTS, work / "in")
        genes = json.loads(lp.gene_list(work).read_text())
        rows, summary = measure(paths, genes)
    write(rows, summary, out, summary_out)
    return summary


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=OUT)
    ap.add_argument("--summary", type=Path, default=SUMMARY)
    args = ap.parse_args(argv)
    print(json.dumps(run(args.out, args.summary), indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
