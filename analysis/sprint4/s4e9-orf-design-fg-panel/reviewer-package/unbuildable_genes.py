"""Why eight library genes have no L1 member: what Ensembl 115 offers against the protein chosen.

CDS resolution (``scratch/becca/appris-principal/resolve_candidate_cds.py``) accepts a
transcript only if its translation equals the selected reference protein, allowing up to
three residue differences between an Ensembl translation and a UniProt sequence of the
same length (reference-allele differences). For each gene it marked unresolved, this
lists the selected UniProt protein and every Ensembl 115 translation of the gene, with its
transcript biotype, length, and how it differs.

    uv run --locked python unbuildable_genes.py      # -> unbuildable_genes.tsv
"""

from __future__ import annotations

import csv
import difflib
import gzip
import sys
from pathlib import Path

HERE = Path(__file__).parent
REPO = HERE.parents[3]
AP = REPO / "scratch/becca/appris-principal"
sys.path.insert(0, str(AP))

import candidate_inputs as inputs  # noqa: E402

OUT = HERE / "unbuildable_genes.tsv"


def unresolved_genes() -> list[str]:
    with (AP / "candidate_cds_coverage.tsv").open() as fh:
        return [
            r["key"].split(":", 1)[1]
            for r in csv.DictReader(fh, delimiter="\t")
            if r["group"] == "library" and r["status"] == "unresolved"
        ]


def uniprot_by_gene(genes: set[str]) -> dict[str, tuple[str, str]]:
    """gene -> (accession, sequence) from the pinned reviewed UniProt table."""
    out = {}
    with gzip.open(inputs.uniprot_tables("human")[0], "rt") as fh:
        head = next(fh).rstrip("\n").split("\t")
        acc, names, seq = (head.index(c) for c in ("Entry", "Gene Names", "Sequence"))
        for line in fh:
            f = line.rstrip("\n").split("\t")
            for g in f[names].split():
                if g in genes and g not in out:
                    out[g] = (f[acc], f[seq])
    return out


def ensembl_translations(genes: set[str]) -> dict[str, list[dict]]:
    """gene -> every Ensembl 115 peptide record for that gene symbol."""
    out: dict[str, list[dict]] = {g: [] for g in genes}
    with gzip.open(inputs.pep_fasta("human"), "rt") as fh:
        rec, seq = None, []
        for line in [*fh, ">"]:
            if line.startswith(">"):
                if rec and rec["gene"] in genes:
                    rec["seq"] = "".join(seq)
                    out[rec["gene"]].append(rec)
                if line.strip() == ">":
                    break
                fields = dict(
                    kv.split(":", 1) for kv in line[1:].split() if ":" in kv and "=" not in kv
                )
                rec = {
                    "transcript": fields.get("transcript", "").split(".")[0],
                    "biotype": fields.get("transcript_biotype", ""),
                    "gene": fields.get("gene_symbol", ""),
                }
                seq = []
            else:
                seq.append(line.strip())
    return out


def describe(ref: str, seq: str) -> str:
    if seq == ref:
        return "identical"
    s = difflib.SequenceMatcher(None, ref, seq, autojunk=False)
    same = sum(b.size for b in s.get_matching_blocks())
    diffs = sum(1 for a, b in zip(ref, seq, strict=False) if a != b)
    if len(seq) == len(ref):
        return f"same length, {diffs} residue(s) differ"
    return f"{len(seq) - len(ref):+d} aa; {same}/{len(ref)} reference residues aligned"


def main() -> int:
    genes = unresolved_genes()
    uni = uniprot_by_gene(set(genes))
    ens = ensembl_translations(set(genes))
    rows = []
    for g in genes:
        acc, ref = uni.get(g, ("", ""))
        if not ens[g]:
            rows.append(
                {
                    "gene": g,
                    "uniprot": acc,
                    "uniprot_aa": len(ref),
                    "transcript": "",
                    "biotype": "",
                    "ensembl_aa": "",
                    "vs_uniprot": "no Ensembl 115 translation",
                }
            )
        for t in sorted(ens[g], key=lambda r: r["transcript"]):
            rows.append(
                {
                    "gene": g,
                    "uniprot": acc,
                    "uniprot_aa": len(ref),
                    "transcript": t["transcript"],
                    "biotype": t["biotype"],
                    "ensembl_aa": len(t["seq"]),
                    "vs_uniprot": describe(ref, t["seq"]) if ref else "no reviewed UniProt entry",
                }
            )
    with OUT.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    for r in rows:
        print("\t".join(str(v) for v in r.values()))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
