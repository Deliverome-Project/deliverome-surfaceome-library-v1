"""Record annotation and N-terminal evidence for selected non-ATG cyno forms.

Read-only: no candidate CDS or selected-protein record is modified. Alignment
results describe homology; they do not establish an unobserved native start.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
from collections import defaultdict
from pathlib import Path

import resolve_compara116_cds as resolver
from Bio import Align
from Bio.Align import substitution_matrices
from Bio.Seq import Seq

GTF = "Macaca_fascicularis.Macaca_fascicularis_6.0.116.gtf.gz"
GENES = {"OPN3", "ADCY2", "PRRG1"}


def attributes(text):
    result = defaultdict(list)
    for key, value in re.findall(r'(\S+) "([^"]*)"', text):
        result[key].append(value)
    return dict(result)


def transcript_evidence(features):
    transcript = next(f for f in features if f["feature"] == "transcript")
    plus = transcript["strand"] == "+"
    coding = sorted(
        (f for f in features if f["feature"] == "CDS"), key=lambda f: f["start"], reverse=not plus
    )
    tags = sorted({t for f in features for t in f["attributes"].get("tag", [])})
    first = coding[0] if coding else None
    start = transcript["start"] if plus else transcript["end"]
    cds_start = (first["start"] if plus else first["end"]) if first else None
    first_bases = []
    for feature in coding:
        positions = (
            range(feature["start"], feature["end"] + 1)
            if plus
            else range(feature["end"], feature["start"] - 1, -1)
        )
        first_bases.extend(list(positions)[: 3 - len(first_bases)])
        if len(first_bases) == 3:
            break
    annotated_start = sorted(
        {
            position
            for f in features
            if f["feature"] == "start_codon"
            for position in range(f["start"], f["end"] + 1)
        }
    )
    return {
        "transcript": transcript,
        "tags": tags,
        "cds_start_NF_tag_present": "cds_start_NF" in tags,
        "cds_end_NF_tag_present": "cds_end_NF" in tags,
        "five_prime_coordinate": start,
        "cds_five_prime_coordinate": cds_start,
        "cds_begins_at_transcript_boundary": start == cds_start,
        "first_cds_feature": first,
        "first_three_CDS_genomic_positions": first_bases,
        "annotated_start_codon_positions": annotated_start,
        "start_codon_at_CDS_start": bool(
            annotated_start and sorted(first_bases) == annotated_start
        ),
        "start_codon_features": [f for f in features if f["feature"] == "start_codon"],
        "stop_codon_features": [f for f in features if f["feature"] == "stop_codon"],
        "five_prime_utr_features": [f for f in features if f["feature"] == "five_prime_utr"],
        "interpretation_limit": "Absent NF tags are not evidence of a complete CDS start.",
    }


def read_gtf(path, wanted_genes):
    output = defaultdict(list)
    wanted = re.compile(
        r'gene_id "(?:' + "|".join(re.escape(g) for g in sorted(wanted_genes)) + r')"'
    )
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith("#") or not wanted.search(line):
                continue
            fields = line.rstrip().split("\t")
            attrs = attributes(fields[8])
            if not set(attrs.get("gene_id", [])) & wanted_genes or "transcript_id" not in attrs:
                continue
            output[attrs["transcript_id"][0]].append(
                {
                    "chromosome": fields[0],
                    "source": fields[1],
                    "feature": fields[2],
                    "start": int(fields[3]),
                    "end": int(fields[4]),
                    "strand": fields[6],
                    "phase": fields[7],
                    "attributes": attrs,
                }
            )
    return output


def alignment_evidence(reference, query):
    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score, aligner.extend_gap_score = -10, -0.5
    alignment = aligner.align(reference, query)[0]
    coordinates = alignment.coordinates.tolist()
    return {
        "reference_length": len(reference),
        "query_length": len(query),
        "reference_first_80": reference[:80],
        "query_first_80": query[:80],
        "score": alignment.score,
        "coordinates": coordinates,
        "first_aligned_reference_offset": int(alignment.aligned[0][0][0]),
        "first_aligned_query_offset": int(alignment.aligned[1][0][0]),
        "n_terminal_alignment": "\n".join(str(alignment).splitlines()[:7]),
        "methionine_added_query_first_80": "M" + query[:79],
        "methionine_added_length": len(query) + 1,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--download-inputs", action="store_true")
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text())
    paths = {}
    for name, item in manifest["files"].items():
        path = args.input_dir / name
        if not path.exists() and args.download_inputs:
            from deliverome_analysis.core import s3

            bucket, key = item["source"].split(":", 1)
            path.parent.mkdir(parents=True, exist_ok=True)
            getattr(s3, bucket).download(key.removeprefix("//"), path)
        if hashlib.sha256(path.read_bytes()).hexdigest() != item["sha256"]:
            raise ValueError(f"Input digest mismatch: {name}")
        paths[name] = path
    resolver.verify_published_checksum(paths[GTF], paths["cyno_gtf_CHECKSUMS"])
    panel = [
        r
        for r in csv.DictReader(paths["ortholog_panel_compara116.tsv"].open(), delimiter="\t")
        if r["species"] == "cyno" and r["human_symbol"] in GENES
    ]
    if len(panel) != 3:
        raise ValueError("Expected exactly the three selected cyno cases")
    wanted_genes = {r["ortholog_ensembl_gene"] for r in panel}
    gtf = read_gtf(paths[GTF], wanted_genes)
    cyno_peptides = defaultdict(list)
    for header, protein in resolver.fasta(paths[f"{resolver.ASSEMBLIES['cyno']}.pep.all.fa.gz"]):
        tags = resolver.tags(header)
        gene = tags["gene"].split(".")[0]
        if gene in wanted_genes:
            transcript = tags["transcript"].split(".")[0]
            cyno_peptides[gene].append(
                {
                    "protein_id": header.split()[0],
                    "transcript": tags["transcript"],
                    "sequence": protein.removesuffix("*"),
                    "annotation": transcript_evidence(gtf[transcript]),
                }
            )
    mouse_peptides = defaultdict(list)
    for header, protein in resolver.fasta(paths[f"{resolver.ASSEMBLIES['mouse']}.pep.all.fa.gz"]):
        tags = resolver.tags(header)
        symbol = tags.get("gene_symbol", "").upper()
        if symbol in GENES:
            mouse_peptides[symbol].append(
                {
                    "protein_id": header.split()[0],
                    "transcript": tags["transcript"],
                    "sequence": protein.removesuffix("*"),
                }
            )
    existing = json.loads(paths["candidate_cds.json"].read_text())
    records = []
    for row in panel:
        gene = row["human_symbol"]
        candidates = cyno_peptides[row["ortholog_ensembl_gene"]]
        selected = next(
            p
            for p in candidates
            if resolver.sha256(p["sequence"]) == row["sequence_sha256"]
            and p["transcript"].split(".")[0] in row["ensembl_transcripts"].split(";")
        )
        sequence = selected["sequence"]
        accession = {"ADCY2": "A0A2K5TV06", "OPN3": "A0A2K5VXZ0", "PRRG1": "A0A2K5UYI3"}[gene]
        uniprot = json.loads(paths[f"{accession}.json"].read_text())
        if resolver.sha256(uniprot["sequence"]["value"]) != row["sequence_sha256"]:
            raise ValueError(f"UniProt sequence does not match selected model: {gene}")
        headers = dict(
            line.split(":", 1)
            for line in paths[f"{accession}.headers"].read_text().lower().splitlines()
            if ":" in line
        )
        if headers.get("x-uniprot-release", "").strip() != manifest["uniprot_release"].lower():
            raise ValueError("Unexpected UniProt annotation release")
        references = []
        for key, value in sorted(existing.items()):
            if key != "lib:" + gene and not key.startswith("ort:mouse:" + gene + ":"):
                continue
            protein = str(Seq(value["cds"]).translate()).removesuffix("*")
            references.append(
                {
                    "candidate_key": key,
                    "transcript": value["transcript"],
                    "species": value["species"],
                    "cds_starts_atg": value["cds"].startswith("ATG"),
                    "protein_sha256": resolver.sha256(protein),
                    "comparison": alignment_evidence(protein, sequence),
                }
            )
        records.append(
            {
                "form_id": row["form_id"],
                "selected_protein_sha256": row["sequence_sha256"],
                "selected": selected,
                "uniprot_annotation": {
                    "accession": accession,
                    "release": manifest["uniprot_release"],
                    "entry_type": uniprot.get("entryType"),
                    "fragment_flag": uniprot.get("proteinDescription", {}).get("flag"),
                    "sequence_cautions": [
                        c
                        for c in uniprot.get("comments", [])
                        if c.get("commentType") == "SEQUENCE CAUTION"
                    ],
                    "nonterminal_or_conflict_features": [
                        f
                        for f in uniprot.get("features", [])
                        if f.get("type")
                        in {"Non-terminal residue", "Sequence conflict", "Non-adjacent residues"}
                    ],
                    "limit": "Missing fragment/caution labels do not prove a complete native N terminus.",
                },
                "designed_reference_comparisons": references,
                "mouse_release116_models": [
                    {**p, "comparison": alignment_evidence(p["sequence"], sequence)}
                    for p in mouse_peptides[gene]
                ],
                "other_cyno_release116_models": [p for p in candidates if p is not selected],
                "unmodified": True,
            }
        )
    output = {
        "inputs": manifest,
        "method": {
            "alignment": "Global pairwise alignment, BLOSUM62, gap open -10, extension -0.5; first optimum returned by Biopython.",
            "coordinates": "Zero-based protein offsets; all exact prefixes, lengths, and complete alignment paths are retained.",
            "limits": "Homology suggests missing or discordant residues; it does not prove the native cyno initiation site. No CDS or protein is modified.",
        },
        "records": records,
    }
    args.out_dir.mkdir(parents=True, exist_ok=True)
    (args.out_dir / "compara116_non_atg_evidence.json").write_text(
        json.dumps(output, indent=2) + "\n"
    )
    # Preserve the literal selected-transcript annotations for independent inspection.
    selected_transcripts = {r["selected"]["transcript"].split(".")[0] for r in records}
    (args.out_dir / "compara116_non_atg_gtf_features.json").write_text(
        json.dumps({t: gtf[t] for t in sorted(selected_transcripts)}, indent=2) + "\n"
    )
    print(
        json.dumps(
            [
                {
                    "form_id": r["form_id"],
                    "transcript": r["selected"]["transcript"],
                    "length": len(r["selected"]["sequence"]),
                    "start_codon_features": len(
                        r["selected"]["annotation"]["start_codon_features"]
                    ),
                    "cds_start_NF": r["selected"]["annotation"]["cds_start_NF_tag_present"],
                    "cds_at_transcript_boundary": r["selected"]["annotation"][
                        "cds_begins_at_transcript_boundary"
                    ],
                }
                for r in records
            ],
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
