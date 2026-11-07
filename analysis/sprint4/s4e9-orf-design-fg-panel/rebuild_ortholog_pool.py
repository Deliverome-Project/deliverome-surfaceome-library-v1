"""Rebuild gene-level ortholog eligibility over every current library gene.

Read the union of all three species' Compara protein homology dumps: Ensembl's
species dumps are arbitrary subsets, and the first species in a row is arbitrary.
Keep one-to-one, high-confidence human/mouse and human/cynomolgus pairs. Missing
gene mappings, ambiguous mappings and missing protein translations remain explicit
rows, rather than silently shrinking the pool. The April query set is not a filter.

This is the upstream pool, before proteoform selection, topology and panel sampling.
It does not treat an existing manifest entry as a prerequisite for eligibility.

    uv run --locked python rebuild_ortholog_pool.py --inputs DIR --out DIR

DIR contains library_members.json (the frozen root plus additions), the current
candidate_manifest_tranches.tsv, and the release-specific peptide and homology files
named by input_names(). No live API call is made during rebuilding.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

SPECIES = {"mouse": "mus_musculus", "cyno": "macaca_fascicularis"}
PEPTIDES = {
    "human": "Homo_sapiens.GRCh38.pep.all.fa.gz",
    "mouse": "Mus_musculus.GRCm39.pep.all.fa.gz",
    "cyno": "Macaca_fascicularis.Macaca_fascicularis_6.0.pep.all.fa.gz",
}
MANIFEST = Path(__file__).with_name("ortholog_rebuild_inputs.json")


def verify_inputs(inputs: Path, manifest: Path = MANIFEST) -> dict:
    spec = json.loads(manifest.read_text())
    for name, entry in spec["files"].items():
        path = inputs / name
        if not path.is_file():
            raise ValueError(f"missing pinned input: {name}; use --stage to fetch it")
        if hashlib.sha256(path.read_bytes()).hexdigest() != entry["sha256"]:
            raise ValueError(f"{name}: checksum differs from {manifest.name}")
    return spec


def stage_inputs(inputs: Path, manifest: Path = MANIFEST) -> None:
    from deliverome_analysis import s3

    spec = json.loads(manifest.read_text())
    inputs.mkdir(parents=True, exist_ok=True)
    for name, entry in spec["files"].items():
        path = inputs / name
        if path.exists() and hashlib.sha256(path.read_bytes()).hexdigest() == entry["sha256"]:
            continue
        bucket, key = entry["source"].split(":", 1)
        getattr(s3, bucket).download(key, path)
    verify_inputs(inputs, manifest)


def input_names(release: int) -> list[str]:
    return [
        "library_members.json",
        "candidate_manifest_tranches.tsv",
        *PEPTIDES.values(),
        *[
            f"Compara.{release}.{sp}.protein_default.homologies.tsv.gz"
            for sp in ("homo_sapiens", *SPECIES.values())
        ],
    ]


def read_peptides(path: Path) -> list[dict]:
    """Keep stable gene/transcript/protein IDs and exact sequences from Ensembl."""
    rows = []
    current = None
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if current is not None:
                    current["sequence"] = "".join(current["sequence"]).rstrip("*")
                    rows.append(current)
                fields = line[1:].split()
                tags = dict(f.split(":", 1) for f in fields[1:] if ":" in f)
                current = {
                    "protein": fields[0].split(".")[0],
                    "gene": tags["gene"].split(".")[0],
                    "transcript": tags["transcript"].split(".")[0],
                    "symbol": tags.get("gene_symbol", ""),
                    "primary_assembly": any(
                        f.startswith(("chromosome:", "primary_assembly:")) for f in fields
                    ),
                    "sequence": [],
                }
            elif current is not None:
                current["sequence"].append(line.strip())
    if current is not None:
        current["sequence"] = "".join(current["sequence"]).rstrip("*")
        rows.append(current)
    return rows


def read_pairs(paths: list[Path]) -> dict[str, dict[str, set[str]]]:
    """Normalize both orientations; union dumps and deduplicate repeated pairs."""
    pairs = {s: defaultdict(set) for s in SPECIES}
    species_key = {v: k for k, v in SPECIES.items()}
    for path in paths:
        with gzip.open(path, "rt") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                if row["homology_type"] != "ortholog_one2one":
                    continue
                if row["is_high_confidence"] != "1":
                    continue
                left, right = row["species"], row["homology_species"]
                if left == "homo_sapiens" and right in species_key:
                    pairs[species_key[right]][row["gene_stable_id"]].add(
                        row["homology_gene_stable_id"]
                    )
                elif right == "homo_sapiens" and left in species_key:
                    pairs[species_key[left]][row["homology_gene_stable_id"]].add(
                        row["gene_stable_id"]
                    )
    return pairs


def human_gene_map(genes: list[str], peptides: list[dict], manifest: list[dict]) -> dict:
    """Prefer the library reference's transcript; resolve other genes by unique symbol.

    Transcript identity survives symbol changes (including SLC22A18 -> SLC67A1).
    A conflict is reported and is never resolved by picking the first gene.
    """
    by_transcript, by_symbol = defaultdict(set), defaultdict(set)
    for row in peptides:
        if not row["primary_assembly"]:
            continue
        by_transcript[row["transcript"]].add(row["gene"])
        by_symbol[row["symbol"]].add(row["gene"])
    transcript_for_gene = defaultdict(set)
    for row in manifest:
        if row["set"] == "library_reference" and row["species"] == "human":
            if row.get("retired_reason") or not row.get("transcript"):
                continue
            if row["key"].startswith("lib:"):
                transcript_for_gene[row["key"].split(":", 1)[1]].add(
                    row["transcript"].split(".")[0]
                )
    result = {}
    for gene in genes:
        tx_genes = set().union(*(by_transcript[t] for t in transcript_for_gene[gene]))
        sym_genes = by_symbol[gene]
        candidates = tx_genes or sym_genes
        # Pseudoautosomal symbols have both an X and a Y copy. A unique selected
        # transcript disambiguates those; a disjoint symbol mapping is a conflict.
        conflict = bool(tx_genes and sym_genes and not tx_genes <= sym_genes)
        result[gene] = {
            "human_ensembl_gene": ";".join(sorted(candidates)),
            "gene_mapping": (
                "conflict"
                if conflict
                else "ambiguous"
                if len(candidates) > 1
                else "missing"
                if not candidates
                else "selected_transcript"
                if tx_genes
                else "unique_primary_assembly_symbol"
            ),
        }
    return result


def build(genes, mapping, pairs, peptides, old_pairs=None, release=116):
    """One accounting row per (library gene, species); no silent exclusions."""
    old_pairs = old_pairs or {}
    by_gene = {s: defaultdict(list) for s in SPECIES}
    for species in SPECIES:
        for row in peptides[species]:
            if row["primary_assembly"]:
                by_gene[species][row["gene"]].append(row)
    rows, candidates = [], {}
    for gene in sorted(genes):
        info = mapping[gene]
        ensg = info["human_ensembl_gene"]
        for species in SPECIES:
            targets = pairs[species].get(ensg, set())
            target = next(iter(targets)) if len(targets) == 1 else ""
            forms = by_gene[species].get(target, [])
            status = "eligible"
            if info["gene_mapping"] in {"missing", "ambiguous", "conflict"}:
                status = "human_gene_mapping_" + info["gene_mapping"]
            elif not targets:
                status = "no_one2one_high_confidence_pair"
            elif len(targets) != 1:
                status = "ambiguous_ortholog_gene"
            elif not forms:
                status = "ortholog_translation_missing"
            if status == "eligible":
                for form in forms:
                    candidates[(species, form["protein"])] = form
            rows.append(
                {
                    "gene": gene,
                    "species": species,
                    "compara_release": release,
                    **info,
                    "ortholog_ensembl_gene": ";".join(sorted(targets)),
                    "ortholog_symbols": ";".join(sorted({r["symbol"] for r in forms} - {""})),
                    "in_april_ortholog_list": ensg in old_pairs.get(species, set()),
                    "n_ensembl_translations": len(forms),
                    "status": status,
                }
            )
    return rows, candidates


def write_tsv(path, rows):
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def run(inputs: Path, out: Path, release: int, manifest: Path = MANIFEST):
    spec = verify_inputs(inputs, manifest)
    if spec["compara_release"] != release:
        raise ValueError("requested Compara release differs from the input manifest")
    files = {n: inputs / n for n in input_names(release)}
    missing = [n for n, p in files.items() if not p.is_file()]
    if missing:
        raise ValueError(f"missing inputs: {missing}")
    genes = json.loads(files["library_members.json"].read_text())
    if len(genes) != len(set(genes)):
        raise ValueError("library gene list contains duplicates")
    manifest = list(csv.DictReader(files["candidate_manifest_tranches.tsv"].open(), delimiter="\t"))
    peptides = {s: read_peptides(files[n]) for s, n in PEPTIDES.items()}
    mapping = human_gene_map(genes, peptides["human"], manifest)
    pairs = read_pairs([p for n, p in files.items() if n.startswith("Compara.")])
    old = {}
    for species in SPECIES:
        name = f"{species}_ortholog_one2one_highconf_non_hla_metadata.csv"
        path = inputs / name
        if path.exists():
            files[name] = path
            old[species] = {r["query_ensembl_gene_id"] for r in csv.DictReader(path.open())}
    rows, candidates = build(genes, mapping, pairs, peptides, old, release)
    assert len(rows) == len(genes) * len(SPECIES)
    out.mkdir(parents=True, exist_ok=True)
    write_tsv(out / "ortholog_pool.tsv", rows)
    with (out / "ortholog_candidates.fasta").open("w") as fh:
        for (species, protein), row in sorted(candidates.items()):
            fh.write(f">{species}|{row['gene']}|{protein}|{row['transcript']}\n{row['sequence']}\n")
    summary = {
        "compara_release": release,
        "library_genes": len(genes),
        "species": {
            s: {
                "status": dict(Counter(r["status"] for r in rows if r["species"] == s)),
                "eligible_not_in_april_list": sum(
                    r["species"] == s
                    and r["status"] == "eligible"
                    and not r["in_april_ortholog_list"]
                    for r in rows
                )
                if s in old
                else None,
            }
            for s in SPECIES
        },
        "ensembl_candidate_proteins": len(candidates),
        "stage": "gene-level eligibility; proteoform selection and panel redraw pending",
        "input_sha256": {n: hashlib.sha256(p.read_bytes()).hexdigest() for n, p in files.items()},
    }
    (out / "ortholog_pool_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--release", type=int, default=116)
    parser.add_argument("--manifest", type=Path, default=MANIFEST)
    parser.add_argument("--stage", action="store_true", help="fetch checksummed inputs from S3")
    args = parser.parse_args()
    if args.stage:
        stage_inputs(args.inputs, args.manifest)
    print(json.dumps(run(args.inputs, args.out, args.release, args.manifest), indent=2))


if __name__ == "__main__":
    main()
