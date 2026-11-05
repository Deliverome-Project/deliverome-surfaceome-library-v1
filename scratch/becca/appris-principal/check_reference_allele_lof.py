"""Evidence for every `reference_allele_lof` override, recomputed from the pinned releases.

A `reference_allele_lof` override says two things, and both are checked here rather than
taken on trust from the override's free-text basis:

1. The GRCh38 primary-assembly copy of the gene cannot be built: every primary-assembly
   transcript is `protein_coding_LoF` in Ensembl 115 or does not translate without an
   internal stop.
2. The transcript the override builds from (on an alternative-haplotype or patch scaffold)
   translates, residue for residue, to the UniProt 2026_03 sequence the override names.

Inputs are only the pinned Ensembl 115 human cds FASTA (`candidate_inputs.cds_fasta`), the
pinned UniProt 2026_03 table (`resolve_candidate_cds.load_uniprot`) and
`reference_overrides.tsv`, so the result is the same years from now. It writes one row per
transcript of each such gene and exits non-zero if either claim fails.

    uv run --locked python check_reference_allele_lof.py
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

from Bio.Seq import Seq

import candidate_inputs as inputs
from resolve_candidate_cds import load_uniprot, read_fasta

HERE = Path(__file__).parent
OVERRIDES = HERE / "reference_overrides.tsv"
OUT = HERE / "reference_allele_lof_evidence.tsv"
#: The transcript each override builds from, as resolved by `resolve_candidate_cds.py`
#: (route "protein match") and recorded in the s4e9 manifest.
BUILT_FROM = {
    "KIR2DS4": "ENST00000617825",
    "OR12D1": "ENST00000703726",
    "OR5AL1": "ENST00000710330",
}


def field(head: str, key: str) -> str:
    for tok in head.split():
        if tok.startswith(key + ":"):
            return tok.split(":", 1)[1]
    return ""


def main() -> int:
    with OVERRIDES.open() as fh:
        overrides = [
            r for r in csv.DictReader(fh, delimiter="\t") if r["basis"] == "reference_allele_lof"
        ]
    cds = read_fasta(inputs.cds_fasta("human"))
    uniprot = load_uniprot("human")
    rows, failures = [], []
    for ov in overrides:
        gene, acc = ov["gene"], ov["form"]
        target = uniprot.get(acc)
        primary_buildable = []
        for enst, (head, seq) in sorted(cds.items()):
            if field(head, "gene_symbol") != gene:
                continue
            region = head.split()[2]
            prot = str(Seq(seq[: len(seq) // 3 * 3]).translate())
            body = prot[:-1] if prot.endswith("*") else prot
            row = {
                "gene": gene,
                "transcript": enst,
                "region": region,
                "primary_assembly": region.startswith("chromosome:"),
                "biotype": field(head, "transcript_biotype"),
                "translated_aa": len(body),
                "internal_stops": body.count("*"),
                "equals_uniprot": body == target,
                "built_from": enst == BUILT_FROM.get(gene),
            }
            rows.append(row)
            if row["primary_assembly"] and row["biotype"] == "protein_coding" and not row["internal_stops"]:
                primary_buildable.append(enst)
        if primary_buildable:
            failures.append(f"{gene}: GRCh38 primary copy is buildable ({primary_buildable})")
        built = [r for r in rows if r["gene"] == gene and r["built_from"]]
        if not built or not built[0]["equals_uniprot"]:
            failures.append(f"{gene}: {BUILT_FROM.get(gene)} does not translate to {acc}")
    with OUT.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    for ov in overrides:
        gene = ov["gene"]
        prim = [r for r in rows if r["gene"] == gene and r["primary_assembly"]]
        print(
            f"{gene}: primary copies "
            + ", ".join(f"{r['transcript']} {r['biotype']} {r['internal_stops']} stops" for r in prim)
            + f"; builds {BUILT_FROM.get(gene)} = {ov['form']}"
        )
    print(f"wrote {OUT.name}: {len(rows)} transcripts")
    for f in failures:
        print("FAIL", f)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
