"""Adjudicate annotated ATG alternatives for uncertain non-ATG library forms.

Read-only. A proposed replacement is a distinct protein form requiring review;
neither an ATG nor absent NF tags establishes native full-length validation.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

import audit_library_non_atg as baseline
import pandas as pd
import resolve_compara116_cds as resolver
from Bio import Align
from Bio.Align import substitution_matrices
from Bio.Seq import Seq

from deliverome_analysis import topology


def compare(reference, query):
    """Exact residue/gap accounting; no undocumented homology cutoff."""
    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score, aligner.extend_gap_score = -10, -0.5
    alignment = aligner.align(reference, query)[0]
    coordinates = alignment.coordinates.tolist()
    substitutions, gaps, identical = [], [], 0
    for index in range(len(coordinates[0]) - 1):
        a, b = coordinates[0][index : index + 2]
        c, d = coordinates[1][index : index + 2]
        if a != b and c != d:
            for i, j in zip(range(a, b), range(c, d), strict=True):
                if reference[i] == query[j]:
                    identical += 1
                else:
                    substitutions.append([i, j, reference[i], query[j]])
        elif a != b:
            gaps.append(
                {
                    "side": "reference_only",
                    "reference": [a, b],
                    "query": [c, d],
                    "location": "N_terminal"
                    if c == 0
                    else "C_terminal"
                    if c == len(query)
                    else "internal",
                }
            )
        else:
            gaps.append(
                {
                    "side": "query_only",
                    "reference": [a, b],
                    "query": [c, d],
                    "location": "N_terminal"
                    if a == 0
                    else "C_terminal"
                    if a == len(reference)
                    else "internal",
                }
            )
    suffix = 0
    for a, b in zip(reversed(reference), reversed(query), strict=False):
        if a != b:
            break
        suffix += 1
    relation = (
        "exact"
        if reference == query
        else "exact_N_terminal_extension_or_truncation"
        if suffix == min(len(reference), len(query))
        else "prefix_replacement_with_identical_remainder"
        if suffix >= min(len(reference), len(query)) - 1
        else "internal_or_other_sequence_differences"
    )
    return {
        "reference_length": len(reference),
        "query_length": len(query),
        "length_delta_query_minus_reference": len(query) - len(reference),
        "identical_aligned_residues": identical,
        "substitution_count": len(substitutions),
        "substitutions": substitutions,
        "gaps": gaps,
        "internal_gap_residues": sum(
            max(g["reference"][1] - g["reference"][0], g["query"][1] - g["query"][0])
            for g in gaps
            if g["location"] == "internal"
        ),
        "C_terminal_gap_residues": sum(
            max(g["reference"][1] - g["reference"][0], g["query"][1] - g["query"][0])
            for g in gaps
            if g["location"] == "C_terminal"
        ),
        "exact_common_suffix_length": suffix,
        "unshared_reference_prefix_length": len(reference) - suffix,
        "unshared_query_prefix_length": len(query) - suffix,
        "relationship": relation,
        "reference_first_80": reference[:80],
        "query_first_80": query[:80],
        "coordinates": coordinates,
        "alignment_score": alignment.score,
        "first_alignment_lines": "\n".join(str(alignment).splitlines()[:7]),
    }


def annotated_start_stop(annotation):
    return bool(
        annotation
        and annotation["start_codon_at_CDS_start"]
        and annotation["stop_codon_features"]
        and not annotation["cds_start_NF_tag_present"]
        and not annotation["cds_end_NF_tag_present"]
    )


def annotation_rank(model):
    annotation = model.get("annotation")
    tags = annotation["tags"] if annotation else []
    principal = [t for t in tags if t.lower().startswith("appris_principal")]
    return (
        not annotated_start_stop(annotation),
        not bool(principal),
        "Ensembl_canonical" not in tags,
        model["transcript"],
    )


def read_models(path, wanted):
    out = {}
    for header, dna in resolver.fasta(path):
        transcript = header.split()[0].split(".")[0]
        if transcript in wanted:
            out[transcript] = {
                "cds": dna,
                "gene": resolver.tags(header)["gene"].split(".")[0],
                "transcript_version": header.split()[0],
            }
    return out


def topology_evidence(sequences, species, snapshot):
    result = topology.resolve(
        sequences,
        tool=topology.TOOL_V1,
        cohorts=topology.cohorts_for_species(species),
        snapshot=snapshot,
    )
    out = {}
    for record in result.to_dict("records"):
        fields = [
            "found",
            "sequence_sha256",
            "tool",
            "tool_version",
            "topology_version",
            "cohort",
            "deeptmhmm_label",
            "tm_helix_count",
            "signal_peptide_length",
            "predicted_surface_membrane",
        ]
        out[record["key"]] = {
            k: (None if pd.isna(record.get(k)) else record.get(k)) for k in fields
        }
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--download-inputs", action="store_true")
    args = parser.parse_args()
    pins = json.loads(args.manifest.read_text())
    paths = baseline.load_paths(pins, args.input_dir, args.download_inputs)
    evidence = json.loads(paths["baseline_evidence.json"].read_text())
    cases = json.loads(paths["compara_evidence.json"].read_text())
    native = json.loads(paths["candidate_cds.json"].read_text())
    if "complete_native_cds.json" in paths:
        extra = json.loads(paths["complete_native_cds.json"].read_text())
        for key in native.keys() & extra.keys():
            if native[key]["cds"] != extra[key]["cds"]:
                raise ValueError(f"Complete-native index changes baseline source: {key}")
        native = {**native, **extra}
    selected = json.loads(paths["compara116_candidate_cds.json"].read_text())
    proposed = list(
        csv.DictReader(paths["candidate_manifest_tranches_review.tsv"].open(), delimiter="\t")
    )
    mapping = {
        r["form_id"]: r
        for r in csv.DictReader(paths["compara116_member_map.tsv"].open(), delimiter="\t")
    }
    members = defaultdict(list)
    for row in proposed:
        members[row["key"]].append(row)
    adopted_native = {mapping[form]["member_key"]: value for form, value in selected.items()}
    constructed = (
        json.loads(paths["candidate_constructs_compara116.json"].read_text())
        if "candidate_constructs_compara116.json" in paths
        else {}
    )
    active_protein_members = defaultdict(list)
    unindexed_active = []
    unindexed_species = Counter()
    index_sources = Counter()
    for row in proposed:
        if row.get("retired_reason") or row["library"] not in {"l1", "l2"}:
            continue
        value = adopted_native.get(row["key"], native.get(row["key"]))
        origin = "native_CDS"
        if value is None and row["key"] in constructed:
            value = {"cds": constructed[row["key"]]["construct"]}
            origin = "constructed_ORF"
        if value is None:
            unindexed_active.append(row["candidate_id"])
            unindexed_species[row["species"]] += 1
            continue
        digest = resolver.sha256(str(Seq(value["cds"]).translate()).removesuffix("*"))
        index_sources[origin] += 1
        active_protein_members[(row["species"], digest)].append(
            {
                "candidate_id": row["candidate_id"],
                "member_key": row["key"],
                "library": row["library"],
                "protein_origin": origin,
                "source_record_key": row["key"],
                "source_DNA_sha256": resolver.sha256(value["cds"]),
                "construct_origin": constructed[row["key"]].get("source")
                if origin == "constructed_ORF"
                else None,
            }
        )
    jobs = {}
    for row in evidence["records"]:
        if row["classification"] in {"annotated_non_AUG_start", "documented_non_AUG_initiation"}:
            continue
        jobs[row["member_key"]] = {
            "source_release": 115,
            "native": native[row["member_key"]],
            "record": row,
            "models": row["same_gene_ATG_models"],
        }
    # Prefer release116 evidence for the selected panel; add new PRRG1.
    for case in cases["records"]:
        form = case["form_id"]
        key = mapping[form]["member_key"]
        annotation = case["selected"]["annotation"]
        jobs[key] = {
            "source_release": 116,
            "native": selected[form],
            "record": {
                "member_key": key,
                "symbol": form.split(":")[1],
                "species": "cyno",
                "gene": annotation["transcript"]["attributes"]["gene_id"][0],
                "transcript": selected[form]["transcript"],
                "classification": baseline.classify(annotation),
                "annotation": annotation,
            },
            "models": [
                {"transcript": p["transcript"].split(".")[0], "annotation": p["annotation"]}
                for p in case["other_cyno_release116_models"]
                if p["sequence"].startswith("M")
            ],
        }
    if len(jobs) != 54:
        raise ValueError(f"Expected54 uncertain union members, got{len(jobs)}")
    needed = defaultdict(set)
    for job in jobs.values():
        needed[(job["source_release"], job["record"]["species"])].update(
            m["transcript"] for m in job["models"]
        )
    models = {
        (release, species): read_models(paths[f"cds{release}_{species}.fa.gz"], ids)
        for (release, species), ids in needed.items()
    }
    snapshot = pd.read_parquet(paths["topology.parquet"])
    records = []
    sequences = defaultdict(dict)
    for key, job in jobs.items():
        row = job["record"]
        protein = str(Seq(job["native"]["cds"]).translate()).removesuffix("*")
        ref_key = "lib:" + row["symbol"]
        human = (
            str(Seq(native[ref_key]["cds"]).translate()).removesuffix("*")
            if ref_key in native
            else None
        )
        active = [r for r in members[key] if not r.get("retired_reason")]
        tranche = (
            "l1"
            if any(r["library"] == "l1" for r in active)
            else "l2"
            if any(r["library"] == "l2" for r in active)
            else "retired_baseline"
        )
        priority = (
            ["OPN3", "ADCY2", "PRRG1"].index(row["symbol"])
            if key.startswith("ort:cyno:") and row["symbol"] in {"OPN3", "ADCY2", "PRRG1"}
            else 3
            if tranche == "l1"
            else 4
            if tranche == "l2"
            else 5
        )
        alternatives = []
        for offered in job["models"]:
            source = models[(job["source_release"], row["species"])].get(offered["transcript"])
            if not source or source["gene"] != row["gene"]:
                raise ValueError(f"Alternative transcript/gene identity mismatch: {key}")
            sequence = baseline.usable_atg_model(source["cds"])
            if sequence is None:
                continue
            comparison = compare(sequence, protein)
            alternatives.append(
                {
                    "transcript": offered["transcript"],
                    "transcript_version": source["transcript_version"],
                    "annotation": offered["annotation"],
                    "annotated_start_and_stop_supported": annotated_start_stop(
                        offered["annotation"]
                    ),
                    "native_full_length_experimentally_validated": False,
                    "protein_sha256": resolver.sha256(sequence),
                    "cds_sha256": resolver.sha256(source["cds"]),
                    "protein_length": len(sequence),
                    "vs_selected": comparison,
                    "vs_human_reference": compare(human, sequence) if human else None,
                    "exact_M_plus_selected": sequence == "M" + protein,
                    "candidate_protein_changes": sequence != protein,
                    "identical_protein_active_members": active_protein_members[
                        (row["species"], resolver.sha256(sequence))
                    ],
                }
            )
            sequences[row["species"]][resolver.sha256(sequence)] = sequence
        alternatives.sort(key=annotation_rank)
        eligible = [a for a in alternatives if a["annotated_start_and_stop_supported"]]
        prefix_candidates = [
            a
            for a in eligible
            if a["vs_selected"]["relationship"]
            in {
                "exact_N_terminal_extension_or_truncation",
                "prefix_replacement_with_identical_remainder",
            }
            and not a["vs_selected"]["C_terminal_gap_residues"]
            and a["protein_length"] >= len(protein)
        ]
        recommendation = (
            "review_N_terminal_replacement_candidate"
            if prefix_candidates
            else "retain_blocked_different_protein_models"
            if eligible
            else "retain_blocked_no_start_stop_supported_alternative"
        )
        annotation_ties = (
            [
                a["transcript"]
                for a in eligible
                if annotation_rank(a)[:-1] == annotation_rank(eligible[0])[:-1]
            ]
            if eligible
            else []
        )
        if prefix_candidates and all(
            a["identical_protein_active_members"] for a in prefix_candidates
        ):
            recommendation = "retain_blocked_replacement_would_duplicate_active_form"
        records.append(
            {
                "member_key": key,
                "symbol": row["symbol"],
                "species": row["species"],
                "proposed_library": tranche,
                "priority": priority,
                "active_candidate_ids": [r["candidate_id"] for r in active],
                "candidate_sets": sorted({r.get("set", "") for r in members[key]}),
                "source_release": job["source_release"],
                "selected_transcript": job["native"]["transcript"],
                "selected_native_cds_sha256": resolver.sha256(job["native"]["cds"]),
                "selected_native_protein_sha256": resolver.sha256(protein),
                "selected_length": len(protein),
                "selected_annotation": row["annotation"],
                "selected_classification": row["classification"],
                "human_reference_key": ref_key if human else None,
                "human_reference_starts_ATG": native[ref_key]["cds"].startswith("ATG")
                if human
                else None,
                "selected_vs_human_reference": compare(human, protein) if human else None,
                "ATG_alternatives": alternatives,
                "start_stop_supported_alternative_count": len(eligible),
                "best_annotation_supported_alternative": eligible[0]["transcript"]
                if eligible
                else None,
                "annotation_tied_top_transcripts": annotation_ties,
                "N_terminal_replacement_candidate_transcripts": [
                    a["transcript"] for a in prefix_candidates
                ],
                "recommendation": recommendation,
                "review_status": "curation_required",
                "approved_for_canonical_activation": False,
                "sequences_changed": False,
                "screen_design_limit": "An alternative/shorter isoform may be part of the intended screen. Annotation of a longer canonical model does not authorize replacing that distinct protein form; replacement can erase isoform diversity or duplicate another member.",
                "active_member_duplicate_check_complete_for_species": not unindexed_species[
                    row["species"]
                ],
                "ATG_alternative_transcripts_already_active": [
                    a["transcript"] for a in alternatives if a["identical_protein_active_members"]
                ],
            }
        )
        sequences[row["species"]][resolver.sha256(protein)] = protein
    topo = {
        species: topology_evidence(seq, species, snapshot) for species, seq in sequences.items()
    }
    for row in records:
        row["selected_exact_topology"] = topo[row["species"]][row["selected_native_protein_sha256"]]
        for candidate in row["ATG_alternatives"]:
            candidate["exact_topology"] = topo[row["species"]][candidate["protein_sha256"]]
    records.sort(key=lambda r: (r["priority"], r["member_key"]))
    summary = {
        "uncertain_members": len(records),
        "tranches": dict(Counter(r["proposed_library"] for r in records)),
        "recommendations": dict(Counter(r["recommendation"] for r in records)),
        "members_with_start_stop_supported_alternative": sum(
            bool(r["start_stop_supported_alternative_count"]) for r in records
        ),
        "L2_recommendations": dict(
            Counter(r["recommendation"] for r in records if r["proposed_library"] == "l2")
        ),
        "approved_for_canonical_activation": 0,
        "sequences_changed": 0,
        "active_member_index_missing_raw_CDS_candidate_ids": unindexed_active,
        "active_member_index_missing_raw_CDS_by_species": dict(unindexed_species),
        "active_member_protein_index_sources": dict(index_sources),
    }
    output = {
        "inputs": pins,
        "summary": summary,
        "method": {
            "alignment": "GlobalBLOSUM62, gapopen-10, extension-0.5; firstoptimum, zero-based fullcoordinates and exactdifferencepositions.",
            "annotation_ranking": "Require explicitfirst-codonstart andstopfeatures and noNFtags; preferAPPRISprincipal thenEnsemblcanonical whenpresent. AbsenceofNFaloneisnotcompleteevidence.",
            "prefix_candidate": "AnnotatedATGmodel is at leastaslong asselected; exactcommonC-terminaltail covers allofshorterprotein exceptatmostitsfirstresidue. Thisisareviewcandidate, notnativeidentityorfunctionalproof.",
            "topology": "ExactproteinSHA, speciescohort andactualDeepTMHMMv1provenance only; no othersequence/predictor fallback and no newpredictions.",
            "limits": "Canonical/principal annotation and predictedtopology donotexperimentallyvalidate nativefull-lengthidentity. All54remaincuration_required.",
        },
        "records": records,
    }
    args.out_dir.mkdir(parents=True, exist_ok=True)
    (args.out_dir / "non_atg_alternative_adjudication.json").write_text(
        json.dumps(output, indent=2) + "\n"
    )
    columns = [
        "member_key",
        "proposed_library",
        "priority",
        "symbol",
        "species",
        "source_release",
        "selected_transcript",
        "selected_length",
        "start_stop_supported_alternative_count",
        "best_annotation_supported_alternative",
        "recommendation",
        "review_status",
        "approved_for_canonical_activation",
    ]
    with (args.out_dir / "non_atg_alternative_adjudication.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                *columns,
                "N_terminal_replacement_candidate_transcripts",
                "ATG_alternative_transcripts_already_active",
                "active_member_duplicate_check_complete_for_species",
            ],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in records:
            writer.writerow(
                {k: row[k] for k in columns}
                | {
                    "N_terminal_replacement_candidate_transcripts": ";".join(
                        row["N_terminal_replacement_candidate_transcripts"]
                    ),
                    "ATG_alternative_transcripts_already_active": ";".join(
                        row["ATG_alternative_transcripts_already_active"]
                    ),
                    "active_member_duplicate_check_complete_for_species": row[
                        "active_member_duplicate_check_complete_for_species"
                    ],
                }
            )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
