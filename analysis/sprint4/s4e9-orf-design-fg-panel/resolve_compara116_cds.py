"""Resolve the selected Compara116 panel to exact protein-verified coding DNA.

Reuse a selected form's existing CDS unchanged when its translation is exact.
Otherwise require an explicitly selected release116 transcript, the selected
stable gene, and exact peptide identity. No length or variant-tolerance route.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from collections import Counter
from pathlib import Path

from Bio.Seq import Seq

STOPS = {"TAA", "TAG", "TGA"}
ASSEMBLIES = {"mouse": "Mus_musculus.GRCm39", "cyno": "Macaca_fascicularis.Macaca_fascicularis_6.0"}


def sha256(value: str) -> str:
    return hashlib.sha256(value.encode()).hexdigest()


def verify_published_checksum(path: Path, checksum_path: Path):
    """Check Ensembl's published BSD sum and 1024-byte block count."""
    rows = [line.split() for line in checksum_path.read_text().splitlines()]
    expected = next((r for r in rows if len(r) == 3 and r[2] == path.name), None)
    if expected is None:
        raise ValueError(f"Missing published checksum for {path.name}")
    data = path.read_bytes()
    checksum = 0
    for byte in data:
        checksum = (((checksum >> 1) | ((checksum & 1) << 15)) + byte) & 0xFFFF
    if (checksum, (len(data) + 1023) // 1024) != (int(expected[0]), int(expected[1])):
        raise ValueError(f"Published Ensembl checksum mismatch: {path.name}")


def fasta(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    header, sequence = None, []
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(sequence).upper()
                header, sequence = line[1:].strip(), []
            else:
                sequence.append(line.strip())
    if header is not None:
        yield header, "".join(sequence).upper()


def tags(header: str) -> dict[str, str]:
    return dict(token.split(":", 1) for token in header.split()[1:] if ":" in token)


def verify_cds(cds: str, protein: str, *, allow_added_stop: bool = True):
    """Return usable DNA and explicit stop provenance, or reject the candidate."""
    cds = cds.upper()
    if not cds or len(cds) % 3 or set(cds) - set("ACGT"):
        return None
    translated = str(Seq(cds).translate())
    native_stop = cds[-3:] if cds[-3:] in STOPS else ""
    if native_stop:
        translated = translated[:-1]
    if "*" in translated or translated != protein:
        return None
    if not native_stop and not allow_added_stop:
        return None
    # Same explicit stop policy as the existing CDS resolver: supply TAA when
    # the release CDS omits its terminal stop, never call it a native codon.
    return {
        "cds": cds if native_stop else cds + "TAA",
        "native_stop": native_stop,
        "stop_added": "" if native_stop else "TAA",
        "protein_sha256": sha256(translated),
    }


def selected_peptides(panel: list[dict], peptide_paths: dict[str, Path]):
    wanted = {s: set() for s in ASSEMBLIES}
    for row in panel:
        wanted[row["species"]].update(row["ensembl_transcripts"].split(";"))
    records = {s: {} for s in ASSEMBLIES}
    for species, path in peptide_paths.items():
        for header, protein in fasta(path):
            fields = tags(header)
            transcript = fields["transcript"].split(".")[0]
            if transcript not in wanted[species]:
                continue
            value = {
                "gene": fields["gene"].split(".")[0],
                "protein_id": header.split()[0],
                "transcript_version": fields["transcript"],
                "protein": protein.removesuffix("*"),
            }
            previous = records[species].get(transcript)
            if previous and previous != value:
                raise ValueError(f"Conflicting peptide records: {transcript}")
            records[species][transcript] = value
    return records


def resolve(panel, existing, manifest_rows, peptides, cds_records):
    id_to_key = {r["candidate_id"]: r["key"] for r in manifest_rows}
    output, provenance = {}, []
    for row in panel:
        form = row["form_id"]
        if form in output:
            raise ValueError(f"Duplicate panel form: {form}")
        species = row["species"]
        selected = []
        for transcript in sorted(set(row["ensembl_transcripts"].split(";")) - {""}):
            pep = peptides[species].get(transcript)
            if (
                pep
                and pep["gene"] == row["ortholog_ensembl_gene"]
                and sha256(pep["protein"]) == row["sequence_sha256"]
            ):
                selected.append((transcript, pep))
        if not selected:
            raise ValueError(
                f"No selected release116 transcript has the exact target protein/gene: {form}"
            )
        proteins = {pep["protein"] for _, pep in selected}
        if len(proteins) != 1:
            raise ValueError(f"Conflicting selected proteins: {form}")
        protein = next(iter(proteins))
        hit, source_key, route, selected_transcript = None, "", "", ""
        keys = sorted(
            {
                id_to_key[cid]
                for cid in row["existing_candidate_ids"].split(";")
                if cid and cid in id_to_key
            }
        )
        for key in keys:
            prior = existing.get(key)
            if not prior or prior["species"] != species:
                continue
            checked = verify_cds(prior["cds"], protein, allow_added_stop=False)
            if checked:
                if prior.get("stop_added"):
                    if prior["stop_added"] != prior["cds"][-3:]:
                        raise ValueError(f"Existing stop provenance disagrees with DNA: {key}")
                    checked["native_stop"], checked["stop_added"] = "", prior["stop_added"]
                # Preserve prior schema and values exactly; selected-sequence
                # identity is recorded separately from legacy variant metadata.
                output[form] = dict(prior)
                hit, source_key, route = checked, key, "exact existing CDS reuse"
                selected_transcript = prior["transcript"]
                break
        if hit is None:
            for transcript, pep in selected:
                record = cds_records[species].get(transcript)
                if (
                    not record
                    or record["gene"] != pep["gene"]
                    or record["transcript_version"] != pep["transcript_version"]
                ):
                    continue
                checked = verify_cds(record["cds"], protein)
                if not checked:
                    continue
                value = {
                    "species": species,
                    "transcript": transcript,
                    "cds": checked["cds"],
                    "origin": "ortholog",
                    "route": "Ensembl116 exact protein",
                    "variant_residues": 0,
                }
                if checked["stop_added"]:
                    value["stop_added"] = checked["stop_added"]
                output[form] = value
                hit, route, selected_transcript = checked, "release116 exact CDS", transcript
                break
        if hit is None:
            raise ValueError(f"No exact unambiguous CDS for selected protein: {form}")
        provenance.append(
            {
                "form_id": form,
                "species": species,
                "human_symbol": row["human_symbol"],
                "selected_gene": row["ortholog_ensembl_gene"],
                "cds_transcript": selected_transcript,
                "selected_release116_transcripts": ";".join(t for t, _ in selected),
                "selected_release116_transcript_versions": ";".join(
                    p["transcript_version"] for _, p in selected
                ),
                "selected_release116_proteins": ";".join(p["protein_id"] for _, p in selected),
                "route": route,
                "source_candidate_key": source_key,
                "selected_protein_sha256": row["sequence_sha256"],
                "translated_protein_sha256": hit["protein_sha256"],
                "protein_identity": "exact",
                "cds_sha256": sha256(output[form]["cds"]),
                "cds_nt": len(output[form]["cds"]),
                "native_terminal_stop": hit["native_stop"],
                "stop_added": hit["stop_added"],
                "stop_added_inherited": bool(source_key and hit["stop_added"]),
                "starts_atg": output[form]["cds"].startswith("ATG"),
                "record_preserved_unchanged": bool(source_key),
            }
        )
    return output, provenance


def cds_by_transcript(path: Path, wanted: set[str]):
    output = {}
    for header, dna in fasta(path):
        transcript = header.split()[0].split(".")[0]
        if transcript not in wanted:
            continue
        value = {
            "gene": tags(header)["gene"].split(".")[0],
            "cds": dna,
            "header": header,
            "transcript_version": header.split()[0],
        }
        if transcript in output and output[transcript] != value:
            raise ValueError(f"Conflicting CDS records: {transcript}")
        output[transcript] = value
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--download-inputs", action="store_true")
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text())
    if manifest["ensembl_release"] != 116:
        raise ValueError("This resolver requires the selected release116 geneset")
    paths = {}
    for name, item in manifest["files"].items():
        path = args.input_dir / name
        if not path.exists() and args.download_inputs:
            from deliverome_analysis import s3

            bucket, key = item["source"].split(":", 1)
            path.parent.mkdir(parents=True, exist_ok=True)
            getattr(s3, bucket).download(key.removeprefix("//"), path)
        if hashlib.sha256(path.read_bytes()).hexdigest() != item["sha256"]:
            raise ValueError(f"Input digest mismatch: {name}")
        paths[name] = path
    panel = list(csv.DictReader(paths["ortholog_panel_compara116.tsv"].open(), delimiter="\t"))
    existing = json.loads(paths["candidate_cds.json"].read_text())
    members = list(csv.DictReader(paths["candidate_manifest_tranches.tsv"].open(), delimiter="\t"))
    peptides = selected_peptides(
        panel, {s: paths[f"{a}.pep.all.fa.gz"] for s, a in ASSEMBLIES.items()}
    )
    cds = {
        s: cds_by_transcript(paths[f"{a}.cds.all.fa.gz"], set(peptides[s]))
        for s, a in ASSEMBLIES.items()
    }
    for species, assembly in ASSEMBLIES.items():
        verify_published_checksum(
            paths[f"{assembly}.cds.all.fa.gz"], paths[f"{species}_cds_CHECKSUMS"]
        )
    output, rows = resolve(panel, existing, members, peptides, cds)
    summary = {
        "selected": len(panel),
        "resolved": len(output),
        "routes": dict(Counter(r["route"] for r in rows)),
        "stop_added": sum(bool(r["stop_added"]) for r in rows),
        "stop_added_new": sum(
            bool(r["stop_added"]) and not r["stop_added_inherited"] for r in rows
        ),
        "non_atg": sum(not r["starts_atg"] for r in rows),
        "input_manifest_sha256": hashlib.sha256(args.manifest.read_bytes()).hexdigest(),
    }
    args.out_dir.mkdir(parents=True, exist_ok=True)
    (args.out_dir / "compara116_candidate_cds.json").write_text(json.dumps(output, indent=2) + "\n")
    (args.out_dir / "compara116_cds_provenance.json").write_text(
        json.dumps({"inputs": manifest, "records": rows, "summary": summary}, indent=2) + "\n"
    )
    with (args.out_dir / "compara116_cds_coverage.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
