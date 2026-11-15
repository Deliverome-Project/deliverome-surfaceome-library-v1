"""Classify every baseline L1/L2 non-ATG CDS using pinned Ensembl115 models.

Produces evidence and explicit repair candidates without changing any sequence.
An annotated or experimentally supported non-AUG start is not a missing start.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

import audit_compara116_non_atg as features
import proteoform_decisions as _audit
import resolve_compara116_cds as resolver
from Bio.Seq import Seq

ASSEMBLIES = {"human": "Homo_sapiens.GRCh38", **resolver.ASSEMBLIES}
TRPV6_SUPPORT = {
    "pmid": "23612980",
    "url": "https://pubmed.ncbi.nlm.nih.gov/23612980/",
    "claim": "Experiments establish ACG initiation with methionine for human and mouse TRPV6; ordinary codon-table translation of the first ACG as threonine is not the initiation event.",
}

#: Rows approved for canonical activation, with who decided and who implemented. The audit
#: proposes; only a row named here is written approved, so rerunning the audit keeps it.
#: apply_reference_fixes.py applies them to the designed-CDS record.
APPROVED = {
    key: {
        "decided_by": "Becca Carlson, 2026-11-12 ('Send all fixes'); evidence PR #280 "
        "REFERENCE_FORM_SCREEN.md",
        "implemented_by": "Claude Code session l1-reviewer-package (subagent), PR #283",
        "decision_reason": "Activate the canonical-initiation variant: replace the annotated "
        "non-ATG start codon with ATG; no residue added (apply_reference_fixes.py)",
    }
    for key in ("lib:SLC22A17", "lib:TRPV6")
}
for _key in APPROVED:
    APPROVED[_key]["alternative"] = "not recorded (before the 2026-11-12 alternative column)"
#: OR8G3 (decided by session l1-proteoform-audit, PR #300; proteoform_decisions_applied.tsv):
#: neither Ensembl 115 transcript starts ATG and the GTF does not annotate a start, so
#: the audit proposes an M-plus variant. Under Becca's 2026-11-12 delegation the start codon
#: is replaced instead, the SLC22A17 / TRPV6 remedy: GTG -> ATG gives UniProt P0DMU2 exactly.
APPROVED["lib:OR8G3"] = {
    k: v for k, v in _audit.attribution("lib:OR8G3").items() if k != "sources"
} | {
    "decision_reason": _audit.attribution("lib:OR8G3")["reason"]
    + " Sources: "
    + _audit.attribution("lib:OR8G3")["sources"]
}
APPROVED["lib:OR8G3"].pop("reason")
#: Approved rows whose action differs from the one the audit derives.
ACTION_OVERRIDE = {
    "lib:OR8G3": (
        "replace_start_with_ATG",
        "Explicit canonical-initiation variant for OR8G3: replace the annotated initiating "
        "GTG with ATG; do not insert an extra amino acid. The result equals UniProt P0DMU2.",
    )
}


def classify(annotation, *, documented_non_aug=False):
    if annotation is None:
        return "annotation_unknown"
    if documented_non_aug:
        return "documented_non_AUG_initiation"
    if annotation["cds_start_NF_tag_present"]:
        return "annotated_5prime_incomplete"
    if annotation["start_codon_at_CDS_start"]:
        return "annotated_non_AUG_start"
    if annotation["start_codon_features"]:
        return "start_feature_not_at_CDS_start"
    if annotation["cds_begins_at_transcript_boundary"]:
        return "missing_start_at_transcript_boundary"
    return "start_annotation_unknown"


def protein_relationship(target, candidate):
    if candidate == "M" + target:
        return "exact_M_plus_target"
    if len(candidate) == len(target) and candidate[1:] == target[1:]:
        return "same_tail_first_residue_differs"
    if candidate.endswith(target):
        return "target_is_exact_suffix"
    return "other_protein_difference"


def usable_atg_model(dna):
    if not dna.startswith("ATG") or len(dna) % 3 or set(dna) - set("ACGT"):
        return None
    protein = str(Seq(dna).translate())
    if not protein.endswith("*") or "*" in protein[:-1]:
        return None
    return protein[:-1]


def load_paths(manifest, folder, download):
    output = {}
    for name, item in manifest["files"].items():
        path = folder / name
        if not path.exists() and download:
            from deliverome_analysis import s3

            bucket, key = item["source"].split(":", 1)
            path.parent.mkdir(parents=True, exist_ok=True)
            getattr(s3, bucket).download(key.removeprefix("//"), path)
        if hashlib.sha256(path.read_bytes()).hexdigest() != item["sha256"]:
            raise ValueError(f"Input digest mismatch: {name}")
        output[name] = path
    return output


def inspect_members(members, existing, paths):
    flagged = [
        r
        for r in members
        if r["library"] in {"l1", "l2"}
        and r["key"] in existing
        and not existing[r["key"]]["cds"].upper().startswith("ATG")
    ]
    records = []
    for species, assembly in ASSEMBLIES.items():
        group = [r for r in flagged if r["species"] == species]
        wanted = {existing[r["key"]]["transcript"].split(".")[0] for r in group}
        selected_models, wanted_genes = {}, set()
        cds_path = paths[f"{assembly}.cds.all.fa.gz"]
        for header, dna in resolver.fasta(cds_path):
            transcript = header.split()[0].split(".")[0]
            if transcript in wanted:
                tags = resolver.tags(header)
                gene = tags["gene"].split(".")[0]
                selected_models[transcript] = {
                    "gene": gene,
                    "transcript_version": header.split()[0],
                    "cds": dna,
                }
                wanted_genes.add(gene)
        gtf = features.read_gtf(paths[f"{assembly}.115.gtf.gz"], wanted_genes)
        alternatives = {}
        for header, dna in resolver.fasta(cds_path):
            tags = resolver.tags(header)
            gene = tags["gene"].split(".")[0]
            if gene not in wanted_genes:
                continue
            protein = usable_atg_model(dna)
            if protein is None:
                continue
            transcript = header.split()[0].split(".")[0]
            annotation = (
                features.transcript_evidence(gtf[transcript]) if transcript in gtf else None
            )
            alternatives[transcript] = {
                "gene": gene,
                "transcript": transcript,
                "transcript_version": header.split()[0],
                "protein": protein,
                "cds_sha256": resolver.sha256(dna),
                "protein_sha256": resolver.sha256(protein),
                "annotation": annotation,
            }
        for member in group:
            key = member["key"]
            value = existing[key]
            transcript = value["transcript"].split(".")[0]
            source = selected_models.get(transcript)
            annotation = (
                features.transcript_evidence(gtf[transcript]) if transcript in gtf else None
            )
            target = str(Seq(value["cds"]).translate()).removesuffix("*")
            symbol = key.split(":")[2 if key.startswith("ort:") else 1]
            documented = (
                symbol == "TRPV6" and species in {"human", "mouse"} and value["cds"][:3] == "ACG"
            )
            classification = classify(annotation, documented_non_aug=documented)
            models = []
            for model in alternatives.values():
                if source is None or model["gene"] != source["gene"]:
                    continue
                a = model["annotation"]
                models.append(
                    {k: v for k, v in model.items() if k not in {"gene", "protein"}}
                    | {
                        "protein_length": len(model["protein"]),
                        "first_40_residues": model["protein"][:40],
                        "relationship": protein_relationship(target, model["protein"]),
                        "missing_prefix_residues_if_exact_suffix": len(model["protein"])
                        - len(target)
                        if model["protein"].endswith(target)
                        else None,
                        "has_annotated_start": bool(a and a["start_codon_features"]),
                        "canonical": bool(a and "Ensembl_canonical" in a["tags"]),
                        "five_prime_incomplete": bool(a and a["cds_start_NF_tag_present"]),
                    }
                )

            def order(model, target_length=len(target)):
                return (
                    model["relationship"] != "exact_M_plus_target",
                    model["relationship"] != "same_tail_first_residue_differs",
                    model["five_prime_incomplete"],
                    not model["has_annotated_start"],
                    not model["canonical"],
                    abs(model["protein_length"] - target_length),
                    model["transcript"],
                )

            models.sort(key=order)
            for model in models[:3]:
                model["comparison"] = features.alignment_evidence(
                    alternatives[model["transcript"]]["protein"], target
                )
            exact_m = [
                m
                for m in models
                if m["relationship"] == "exact_M_plus_target"
                and m["has_annotated_start"]
                and not m["five_prime_incomplete"]
            ]
            if classification in {"documented_non_AUG_initiation", "annotated_non_AUG_start"}:
                action = "Preserve native non-AUG record; review an explicit initiation-codon substitution, not an automatic extra-M prepend."
            elif exact_m:
                action = "A same-gene annotated ATG model encodes exactly M plus target; review that pinned model as an explicit protein-changing repair candidate."
            elif models:
                action = "Review same-gene ATG alternative as a distinct protein form; adding M alone is not established to restore the native N terminus."
            else:
                action = "Hold initiation repair until native start evidence or an independently verified full-length CDS is available."
            raw_matches = bool(
                source
                and (
                    value["cds"] == source["cds"]
                    or value["cds"] == source["cds"] + value.get("stop_added", "")
                )
            )
            substitute = documented or bool(annotation and annotation["start_codon_at_CDS_start"])
            proposed_dna = "ATG" + (value["cds"][3:] if substitute else value["cds"])
            proposal = {
                "kind": "replace_first_codon_ATG" if substitute else "prepend_ATG_uncertain_model",
                "native_cds_sha256": resolver.sha256(value["cds"]),
                "proposed_cds_sha256": resolver.sha256(proposed_dna),
                "proposed_ordinary_protein_sha256": resolver.sha256(
                    str(Seq(proposed_dna).translate()).removesuffix("*")
                ),
                "proposed_first_40_residues": str(Seq(proposed_dna[:120]).translate()),
                "native_full_length_restoration_proven": False,
                "applied": False,
                "meaning": "Explicit expression-engineering candidate only; unresolved model defects remain.",
            }
            records.append(
                {
                    "candidate_id": member["candidate_id"],
                    "member_key": key,
                    "library": member["library"],
                    "species": species,
                    "symbol": symbol,
                    "transcript": transcript,
                    "gene": source["gene"] if source else None,
                    "raw_start_codon": value["cds"][:3],
                    "raw_cds_sha256": resolver.sha256(value["cds"]),
                    "ordinary_translation_sha256": resolver.sha256(target),
                    "ordinary_translation_length": len(target),
                    "ordinary_translation_first_80": target[:80],
                    "raw_matches_pinned_transcript_cds": raw_matches,
                    "classification": classification,
                    "annotation": annotation,
                    "same_gene_ATG_models": models,
                    "same_gene_ATG_model_count": len(models),
                    "exact_M_plus_target_with_annotated_start": [m["transcript"] for m in exact_m],
                    "experimental_support": TRPV6_SUPPORT if documented else None,
                    "recommended_action": action,
                    "start_feature_at_CDS_start": annotation["start_codon_at_CDS_start"]
                    if annotation
                    else "unknown",
                    "engineering_proposal": proposal,
                    "sequence_changed": False,
                }
            )
    return records


def boundary_decisions(
    records, *, manifest, compara_evidence=None, selected_cds=None, member_map=None
):
    """Emit checksum-bound expression-engineering proposals without applying them."""
    union = {r["member_key"]: dict(r) for r in records}
    if compara_evidence is not None:
        mapped = {r["form_id"]: r["member_key"] for r in member_map}
        for case in compara_evidence["records"]:
            form = case["form_id"]
            value = selected_cds[form]
            protein = str(Seq(value["cds"]).translate()).removesuffix("*")
            if resolver.sha256(protein) != case["selected_protein_sha256"]:
                raise ValueError(
                    f"Proposed CDS is not the evidence-selected native protein: {form}"
                )
            if value["cds"].startswith("ATG"):
                continue
            key = mapped[form]
            digest = resolver.sha256(value["cds"])
            if key in union:
                if union[key]["raw_cds_sha256"] != digest:
                    raise ValueError(f"Baseline/proposed native CDS collision: {key}")
                union[key]["additional_compara116_evidence"] = case["selected"]["transcript"]
                continue
            annotation = case["selected"]["annotation"]
            union[key] = {
                "member_key": key,
                "raw_cds_sha256": digest,
                "transcript": case["selected"]["transcript"],
                "symbol": form.split(":")[1],
                "raw_start_codon": value["cds"][:3],
                "annotation": annotation,
                "classification": classify(annotation),
                "experimental_support": None,
                "source_release": 116,
            }
    decisions = []
    for key, record in sorted(union.items()):
        annotation = record["annotation"]
        documented = record["classification"] == "documented_non_AUG_initiation"
        substitute = documented or bool(annotation and annotation["start_codon_at_CDS_start"])
        if substitute:
            action = "replace_start_with_ATG"
            completeness = (
                "documented_non_AUG_initiation" if documented else "annotated_non_AUG_start"
            )
            reason = f"Explicit canonical-initiation variant for {record['symbol']}: replace the annotated initiating {record['raw_start_codon']} with ATG; do not insert an extra amino acid."
        else:
            action = "prepend_ATG_for_expression"
            completeness = record["classification"]
            reason = f"Explicit M-plus expression variant for {record['symbol']}; native start/completeness remains unresolved and any existing model defect is retained."
        release = record.get("source_release", 115)
        source_files = compara_evidence["inputs"]["files"] if release == 116 else manifest["files"]
        evidence = {
            "release": release,
            "transcript": record["transcript"],
            "classification": record["classification"],
            "start_codon_at_CDS_start": annotation["start_codon_at_CDS_start"]
            if annotation
            else None,
            "cds_start_NF_tag": annotation["cds_start_NF_tag_present"] if annotation else None,
            "cds_end_NF_tag": annotation["cds_end_NF_tag_present"] if annotation else None,
            "cds_at_transcript_boundary": annotation["cds_begins_at_transcript_boundary"]
            if annotation
            else None,
            "native_full_length_restoration_proven": False,
            "GTF_sources": [v for n, v in source_files.items() if n.endswith(".gtf.gz")],
            "experimental_support": record.get("experimental_support"),
            "additional_release116_transcript": record.get("additional_compara116_evidence"),
        }
        if key in ACTION_OVERRIDE:
            action, reason = ACTION_OVERRIDE[key]
        decisions.append(
            {
                "key": key,
                "original_cds_sha256": record["raw_cds_sha256"],
                "action": action,
                "completeness": completeness,
                "reason": reason,
                "evidence": json.dumps(evidence, sort_keys=True, separators=(",", ":")),
                "review_status": "experimentally_supported_initiator"
                if documented
                else "annotation_supported_pending_review"
                if substitute
                else "curation_required",
                "approved_for_canonical_activation": key in APPROVED,
                **APPROVED.get(
                    key,
                    {
                        "decided_by": "",
                        "implemented_by": "",
                        "decision_reason": "",
                        "alternative": "",
                    },
                ),
            }
        )
    return decisions


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--download-inputs", action="store_true")
    parser.add_argument("--compara-evidence", type=Path)
    parser.add_argument("--selected-cds", type=Path)
    parser.add_argument("--member-map", type=Path)
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text())
    paths = load_paths(manifest, args.input_dir, args.download_inputs)
    existing = json.loads(paths["candidate_cds.json"].read_text())
    members = list(csv.DictReader(paths["candidate_manifest_tranches.tsv"].open(), delimiter="\t"))
    for species, assembly in ASSEMBLIES.items():
        resolver.verify_published_checksum(
            paths[f"{assembly}.115.gtf.gz"], paths[f"{species}_gtf115_CHECKSUMS"]
        )
    rows = inspect_members(members, existing, paths)
    summary = {
        "baseline_non_ATG_members": len(rows),
        "libraries": dict(Counter(r["library"] for r in rows)),
        "species": dict(Counter(r["species"] for r in rows)),
        "classifications": dict(Counter(r["classification"] for r in rows)),
        "with_same_gene_ATG_model": sum(bool(r["same_gene_ATG_models"]) for r in rows),
        "with_exact_M_plus_target_annotated_model": sum(
            bool(r["exact_M_plus_target_with_annotated_start"]) for r in rows
        ),
        "raw_matches_pinned_transcript_cds": sum(
            r["raw_matches_pinned_transcript_cds"] for r in rows
        ),
        "input_manifest_sha256": hashlib.sha256(args.manifest.read_bytes()).hexdigest(),
        "sequences_changed": 0,
    }
    if len(rows) != manifest["expected_baseline_non_ATG"]:
        raise ValueError("Baseline cohort count disagrees with pinned manifest")
    output = {
        "inputs": manifest,
        "summary": summary,
        "method": {
            "classification": "Experimental initiation evidence, NF tag, annotated start codon, transcript boundary, otherwise unknown. Absence of an NF tag does not establish completeness.",
            "model_filter": "Same stable gene; ACGT-only, in-frame CDS; ATG start, terminal stop and no internal stops. No topology or expression suitability inferred.",
            "ranking": "Exact M+target relationship, first-residue-only difference, no start-NF, annotated start, canonical tag, closest length, transcript ID. Alignment is global BLOSUM62 with gap open -10 and extension -0.5 for the first three models.",
            "limits": "Ordinary codon-table translation does not model initiator-M decoding at non-AUG codons. An ATG alternative is a review candidate, not an automatic isoform substitution.",
        },
        "records": rows,
    }
    supplementary = [args.compara_evidence, args.selected_cds, args.member_map]
    if any(supplementary) and not all(supplementary):
        raise ValueError("Union decisions require all three supplementary inputs")
    cases = json.loads(args.compara_evidence.read_text()) if args.compara_evidence else None
    selected_cds = json.loads(args.selected_cds.read_text()) if args.selected_cds else None
    member_map = (
        list(csv.DictReader(args.member_map.open(), delimiter="\t")) if args.member_map else None
    )
    decisions = boundary_decisions(
        rows,
        manifest=manifest,
        compara_evidence=cases,
        selected_cds=selected_cds,
        member_map=member_map,
    )
    output["decision_provenance"] = {
        "applied": False,
        "count": len(decisions),
        "actions": dict(Counter(d["action"] for d in decisions)),
        "supplementary_input_sha256": {
            p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in supplementary if p
        },
        "scope": "Union of baseline non-ATG members and additional proposed Compara116 non-ATG members; native records remain unchanged.",
    }
    args.out_dir.mkdir(parents=True, exist_ok=True)
    (args.out_dir / "library_non_atg_evidence.json").write_text(json.dumps(output, indent=2) + "\n")
    with (args.out_dir / "library_cds_boundary_decisions.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "key",
                "original_cds_sha256",
                "action",
                "completeness",
                "reason",
                "evidence",
                "review_status",
                "approved_for_canonical_activation",
                "decided_by",
                "implemented_by",
                "decision_reason",
                "alternative",
            ],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(decisions)
    columns = [
        "candidate_id",
        "member_key",
        "library",
        "species",
        "symbol",
        "transcript",
        "gene",
        "raw_start_codon",
        "classification",
        "same_gene_ATG_model_count",
        "raw_matches_pinned_transcript_cds",
        "recommended_action",
        "start_feature_at_CDS_start",
    ]
    with (args.out_dir / "library_non_atg_repair_candidates.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                *columns,
                "repair_candidate_transcripts",
                "exact_M_plus_target_with_annotated_start",
            ],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {k: row[k] for k in columns}
                | {
                    "repair_candidate_transcripts": ";".join(
                        m["transcript"] for m in row["same_gene_ATG_models"][:3]
                    ),
                    "exact_M_plus_target_with_annotated_start": ";".join(
                        row["exact_M_plus_target_with_annotated_start"]
                    ),
                }
            )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
