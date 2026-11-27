"""Read-only, evidence-bound recommendations for seven priority L1 forms.

Reference proteins, ortholog counterparts and an intentional alternative-form
screen require different decisions. This report changes no sequence or member.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from collections import defaultdict
from pathlib import Path

import audit_compara116_non_atg as features
import audit_library_non_atg as baseline
import audit_non_atg_alternatives as alternatives
import resolve_compara116_cds as resolver
from Bio.Seq import Seq

from deliverome_analysis.core.provenance import provenance

HUMAN = {"lib:OR8G3": "P0DMU2", "lib:SLC22A17": "Q8WUG5-2", "lib:TRPV6": "Q9H1D0"}
KEYS = [
    *HUMAN,
    "ort:cyno:OPN3:A0A2K5VXZ0",
    "ort:cyno:ADCY2:A0A2K5TV06",
    "ort:cyno:PRRG1:ENSMFAP00000048244",
    "iso:SLC9D1:ENST00000622371",
]


def intended_comparison(protein, intended):
    return {
        "ordinary_translation_matches": protein == intended,
        "first_codon_ATG_replacement_matches": "M" + protein[1:] == intended,
        "extra_M_prepend_matches": "M" + protein == intended,
        "native_translation_sha256": resolver.sha256(protein),
        "intended_protein_sha256": resolver.sha256(intended),
        "native_length": len(protein),
        "intended_length": len(intended),
        "aligned_comparison": alternatives.compare(intended, protein),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--download-inputs", action="store_true")
    args = parser.parse_args()
    if args.out_dir.exists():
        raise FileExistsError(f"Output directory already exists: {args.out_dir}")
    pins = json.loads(args.manifest.read_text())
    paths = baseline.load_paths(pins, args.input_dir, args.download_inputs)
    native = json.loads(paths["candidate_cds.json"].read_text())
    old = json.loads(paths["baseline_evidence.json"].read_text())
    old = {r["member_key"]: r for r in old["records"]}
    adjudication = json.loads(paths["alternative_adjudication.json"].read_text())
    adjudication = {r["member_key"]: r for r in adjudication["records"]}
    cases = json.loads(paths["compara_evidence.json"].read_text())
    cases = {r["form_id"]: r for r in cases["records"]}
    references = {
        r["symbol"]: r
        for r in csv.DictReader(paths["reference_selection.tsv"].open(), delimiter="\t")
    }
    human_genes = {old[k]["gene"] for k in HUMAN}
    gtf = features.read_gtf(paths["Homo_sapiens.GRCh38.115.gtf.gz"], human_genes)
    models = defaultdict(list)
    for header, dna in resolver.fasta(paths["Homo_sapiens.GRCh38.cds.all.fa.gz"]):
        gene = resolver.tags(header)["gene"].split(".")[0]
        if gene not in human_genes:
            continue
        transcript = header.split()[0].split(".")[0]
        if len(dna) % 3 or set(dna) - set("ACGT"):
            continue
        protein = str(Seq(dna).translate()).removesuffix("*")
        if "*" in protein:
            continue
        annotation = features.transcript_evidence(gtf[transcript]) if transcript in gtf else None
        models[gene].append(
            {
                "transcript": transcript,
                "transcript_version": header.split()[0],
                "start_codon": dna[:3],
                "protein": protein,
                "annotation": annotation,
                "cds_sha256": resolver.sha256(dna),
                "transcript_biotype": resolver.tags(header).get("transcript_biotype"),
            }
        )
    output = []
    for key in KEYS:
        symbol = key.split(":")[2 if key.startswith("ort:") else 1]
        record = {
            "member_key": key,
            "symbol": symbol,
            "approved_for_canonical_activation": False,
            "sequence_changed": False,
            "review_status": "pending_scientific_review",
        }
        if key in HUMAN:
            accession = HUMAN[key]
            entry_acc = accession.split("-")[0]
            entry = json.loads(paths[f"{entry_acc}.json"].read_text())
            headers = dict(
                line.split(":", 1)
                for line in paths[f"{entry_acc}.headers"].read_text().lower().splitlines()
                if ":" in line
            )
            if headers.get("x-uniprot-release", "").strip() != pins["uniprot_release"].lower():
                raise ValueError("Unexpected UniProt release")
            intended = (
                next(resolver.fasta(paths[f"{accession}.fasta"]))[1]
                if "-" in accession
                else entry["sequence"]["value"]
            )
            protein = str(Seq(native[key]["cds"]).translate()).removesuffix("*")
            comparison = intended_comparison(protein, intended)
            matching = []
            for model in models[old[key]["gene"]]:
                if model["protein"] == intended or "M" + model["protein"][1:] == intended:
                    matching.append(
                        {k: v for k, v in model.items() if k != "protein"}
                        | {
                            "ordinary_translation_matches_intended": model["protein"] == intended,
                            "first_M_normalized_matches_intended": "M" + model["protein"][1:]
                            == intended,
                            "annotated_ATG_exact_intended": model["start_codon"] == "ATG"
                            and model["protein"] == intended
                            and alternatives.annotated_start_stop(model["annotation"]),
                        }
                    )
            exact = [m["transcript"] for m in matching if m["annotated_ATG_exact_intended"]]
            record.update(
                {
                    "form_role": "human_reference",
                    "selected_transcript": native[key]["transcript"],
                    "selected_start_codon": native[key]["cds"][:3],
                    "native_cds_sha256": resolver.sha256(native[key]["cds"]),
                    "selected_annotation": old[key]["annotation"],
                    "reference_selection": references[symbol],
                    "intended_database_accession": accession,
                    "uniprot_entry_type": entry.get("entryType"),
                    "intended_sequence_comparison": comparison,
                    "intended_sequence_matching_transcripts": matching,
                    "better_supported_exact_intended_ATG_transcripts": exact,
                    "curated_caution_evidence": [
                        {
                            "comment_type": c["commentType"],
                            "evidence_codes": sorted(
                                {
                                    e["evidenceCode"]
                                    for t in c.get("texts", [])
                                    for e in t.get("evidences", [])
                                }
                            ),
                            "source_comment_index": i,
                        }
                        for i, c in enumerate(entry.get("comments", []))
                        if c["commentType"] in {"CAUTION", "SEQUENCE CAUTION"}
                    ],
                }
            )
            if not comparison["first_codon_ATG_replacement_matches"]:
                raise ValueError(f"Priority recommendation premise fails: {key}")
            if symbol == "TRPV6":
                recommendation = "supported_initiator_codon_replacement"
                confidence = "Experimental non-AUG initiator-M evidence plus exact intended-protein identity after ACG-to-ATG substitution."
                reason = f"Keep the selected {len(protein)}-aa form; replace the first codon, not prepend a residue or switch to the downstream shorter AUG form."
                record["experimental_support"] = baseline.TRPV6_SUPPORT
            elif symbol == "SLC22A17":
                recommendation = (
                    "annotation_and_curated_sequence_supported_replacement_pending_review"
                )
                confidence = "Annotated initiation position and reviewed isoform identity; the UniProt non-AUG statement is curator inference, not direct initiator-M validation here."
                reason = f"CTG-to-ATG produces the exact selected {accession} {len(intended)}-aa protein. Do not substitute the different {entry['sequence']['length']}-aa canonical isoform."
            else:
                recommendation = "curated_intended_protein_start_variant_requires_conflict_review"
                if not matching or any(
                    m["transcript_biotype"] != "protein_coding_LoF" for m in matching
                ):
                    raise ValueError("OR8G3 LoF-versus-initiation review premise changed")
                confidence = f"Reviewed sequence supports an M-first {len(intended)}-aa protein; UniProt suggests possible GUG initiation by curator inference, while the matching Ensembl transcripts are protein_coding_LoF."
                reason = "GTG-to-ATG, not an extra-M prepend, matches the reviewed protein exactly. Resolve the noncanonical-initiation versus LoF annotation conflict before accepting a reference-form rescue."
            record.update(
                {
                    "recommendation": recommendation,
                    "confidence_basis": confidence,
                    "reason": reason,
                    "ortholog_redraw_needed": "not_applicable",
                    "native_full_length_experimentally_validated": False,
                    "initiator_identity_experimentally_supported": symbol == "TRPV6",
                }
            )
        else:
            evidence = adjudication[key]
            record.update(
                {
                    "form_role": "alternative_isoform_screen"
                    if key.startswith("iso:")
                    else "ortholog_counterpart",
                    "selected_transcript": evidence["selected_transcript"],
                    "selected_annotation": evidence["selected_annotation"],
                    "native_cds_sha256": evidence["selected_native_cds_sha256"],
                    "alternative_adjudication": evidence,
                    "better_supported_exact_intended_ATG_transcripts": [],
                    "native_full_length_experimentally_validated": False,
                }
            )
            if symbol == "OPN3":
                model_lengths = sorted(
                    {
                        a["protein_length"]
                        for a in evidence["ATG_alternatives"]
                        if a["annotated_start_and_stop_supported"]
                    }
                )
                if (
                    not model_lengths
                    or max(model_lengths) >= evidence["selected_length"]
                    or evidence["N_terminal_replacement_candidate_transcripts"]
                ):
                    raise ValueError("OPN3 shorter-model review premise changed")
                recommendation = "independent_full_length_curation_or_ortholog_redraw_review"
                reason = f"The annotated ATG models are {model_lengths} amino acids versus selected {evidence['selected_length']} and human {evidence['selected_vs_human_reference']['reference_length']}; their differences extend beyond initiation. They are not equivalent full-length counterpart repairs."
            elif symbol == "ADCY2":
                human_comparison = evidence["selected_vs_human_reference"]
                if (
                    human_comparison["reference_length"] != evidence["selected_length"] + 1
                    or not any(
                        g["location"] == "N_terminal"
                        and g["reference"] == [0, 1]
                        and g["query"] == [0, 0]
                        for g in human_comparison["gaps"]
                    )
                    or any(a["exact_M_plus_selected"] for a in evidence["ATG_alternatives"])
                ):
                    raise ValueError("ADCY2 inferred-initiator review premise changed")
                recommendation = "validate_native_start_or_review_ortholog_form_repick"
                reason = "A single added M is a plausible inference from the human comparison, not a validated start. Available ATG models contain internal sequence losses; none encodes exactly M plus the selected protein."
            elif symbol == "PRRG1":
                candidates = [
                    a
                    for a in evidence["ATG_alternatives"]
                    if a["transcript"] in evidence["N_terminal_replacement_candidate_transcripts"]
                ]
                if not candidates or not all(
                    a["annotated_start_and_stop_supported"] for a in candidates
                ):
                    raise ValueError("PRRG1 annotated-alternative review premise changed")
                descriptions = [
                    f"{a['transcript']} ({a['protein_length']} aa, {'canonical' if 'Ensembl_canonical' in a['annotation']['tags'] else 'alternative'})"
                    for a in candidates
                ]
                cached = [
                    f"{a['transcript']}:{a['exact_topology']['deeptmhmm_label']}"
                    for a in candidates
                    if a["exact_topology"]["found"]
                ]
                recommendation = "review_annotated_ortholog_form_repick_candidates"
                reason = f"ATG models {descriptions} have annotated starts/stops and preserved selected C-terminal sequence. Exact cached v1 topology: {cached}. These are distinct N-terminal forms, not an accepted edit of the {evidence['selected_length']}-aa selection."
            else:
                canonical = [
                    a
                    for a in evidence["ATG_alternatives"]
                    if a["annotation"] and "Ensembl_canonical" in a["annotation"]["tags"]
                ]
                if not {"cds_start_NF", "mRNA_start_NF"} <= set(
                    evidence["selected_annotation"]["tags"]
                ) or not any(a["identical_protein_active_members"] for a in canonical):
                    raise ValueError(
                        "SLC9D1 partial-screen-versus-existing-reference premise changed"
                    )
                recommendation = "retain_separate_isoform_intent_block_native_completeness_claim"
                reason = f"The {evidence['selected_length']}-aa alternative has cds_start_NF and mRNA_start_NF. Canonical model lengths {sorted({a['protein_length'] for a in canonical})} are already represented; replacing this screen arm would erase its intended alternative-form distinction."
            record.update(
                {
                    "recommendation": recommendation,
                    "confidence_basis": "Pinned transcript annotations and exact sequence comparisons; no native-start validation inferred from absence of NF tags.",
                    "reason": reason,
                    "ortholog_redraw_needed": "review_required_if_no_supported_intended_counterpart"
                    if key.startswith("ort:")
                    else "not_applicable",
                }
            )
        output.append(record)
    payload = {
        "inputs": pins,
        "records": output,
        "limits": "Recommendations are not applied decisions. First-codon replacement, choosing a different transcript and redrawing a panel member are distinct actions. Canonical annotation is not permission to replace a distinct alternative isoform.",
    }
    args.out_dir.mkdir(parents=True, exist_ok=False)
    (args.out_dir / "priority_cds_form_recommendations.json").write_text(
        json.dumps(payload, indent=2) + "\n"
    )
    columns = [
        "member_key",
        "form_role",
        "selected_transcript",
        "recommendation",
        "confidence_basis",
        "reason",
        "ortholog_redraw_needed",
        "approved_for_canonical_activation",
    ]
    with (args.out_dir / "priority_cds_form_recommendations.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({c: r[c] for c in columns} for r in output)
    stamp = provenance()
    run = {
        "provenance": stamp,
        "producing_script_commit": stamp.get("git-commit"),
        "producing_script_git_dirty": stamp.get("git-dirty"),
        "producing_script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "helper_script_sha256": {
            Path(module.__file__).name: hashlib.sha256(
                Path(module.__file__).read_bytes()
            ).hexdigest()
            for module in [features, baseline, alternatives, resolver]
        },
        "argv": list(sys.argv),
        "input_manifest_sha256": hashlib.sha256(args.manifest.read_bytes()).hexdigest(),
        "output_sha256": {
            name: hashlib.sha256((args.out_dir / name).read_bytes()).hexdigest()
            for name in [
                "priority_cds_form_recommendations.json",
                "priority_cds_form_recommendations.tsv",
            ]
        },
    }
    (args.out_dir / "priority_cds_form_review_run.json").write_text(
        json.dumps(run, indent=2) + "\n"
    )
    print(
        json.dumps(
            {
                "priority_forms": len(output),
                "recommendations": {r["symbol"]: r["recommendation"] for r in output},
                "sequences_changed": 0,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
