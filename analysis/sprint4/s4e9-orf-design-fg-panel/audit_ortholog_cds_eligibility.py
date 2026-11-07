"""Gate every Compara116 target-gene transcript before protein correspondence scoring.

Read-only and fail closed: uncertain annotation stays in the evidence ledger but
does not enter an automatic draw. NF is an unconfirmed boundary, not proof of
truncation. No native CDS is repaired and no candidate/member record is activated.
The species/protein-hash output preserves a separate record for each stable gene.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

import audit_library_cds as mechanical
import audit_library_completeness as completeness

RULE_VERSION = "ortholog-cds-eligibility-v2"
SPECIES = ("mouse", "cyno")
WAIVABLE_REVIEW_REASONS = {
    "non_ATG_initiation_review",
    "start_feature:no_start_feature",
    "cds_start_NF_not_confirmed",
    "mRNA_start_NF_not_confirmed",
}
RULES = {
    "ORTHCDS001": "Union all CDS, peptide and GTF transcripts in every gene-level eligible target gene.",
    "ORTHCDS002": "Require exact CDS translation/peptide and stable gene/transcript/protein version coherence.",
    "ORTHCDS003": "Require unambiguous, in-frame native ATG CDS, terminal stop, no internal stop.",
    "ORTHCDS004": "Require strand/version/sense-length coherent coordinate-supported boundaries, phase0, no CDS NF tags; mRNA NF remains advisory for CDS-only expression.",
    "ORTHCDS005": "Non-protein_coding biotypes and unresolved evidence require review; no biological inactivity inferred.",
    "ORTHCDS006": "Equivalent encodings rank by version-matched GTF APPRIS or protein-validated external APPRIS, Ensembl canonical, TSL, then stable transcript ID.",
    "ORTHCDS007": "Documented exceptions bind species/gene/transcript version/CDS+protein hashes and explicit review waivers.",
}


def digest(value: str) -> str:
    return hashlib.sha256(value.encode()).hexdigest()


def stable(value: str) -> str:
    return value.split(".", 1)[0]


def read_sequence_models(path: Path, genes: set[str], kind: str) -> dict:
    """Preserve all distinct model records; conflicting identities stay reviewable."""
    result = defaultdict(list)
    header, sequence = None, []

    def store():
        if header is None:
            return
        tags = dict(t.split(":", 1) for t in header.split()[1:] if ":" in t)
        if stable(tags.get("gene", "")) not in genes:
            return
        transcript_version = header.split()[0] if kind == "cds" else tags.get("transcript", "")
        transcript = stable(transcript_version)
        if not transcript:
            raise ValueError(f"Missing transcript identity in {path.name}")
        record = {
            "header": header,
            "tags": tags,
            "transcript_version": transcript_version,
            "sequence": "".join(sequence).upper(),
            "model_id": header.split()[0],
        }
        if record not in result[transcript]:
            result[transcript].append(record)

    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                store()
                header, sequence = line[1:].strip(), []
            elif header is not None:
                sequence.append(line.strip())
    store()
    return {t: sorted(records, key=lambda r: r["header"]) for t, records in result.items()}


def read_target_gtf(path: Path, genes: set[str]) -> dict:
    result = defaultdict(list)
    gene_pattern = re.compile(r'\bgene_id "([^"]+)"')
    transcript_pattern = re.compile(r'\btranscript_id "([^"]+)"')
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] not in completeness.FEATURES:
                continue
            gene = gene_pattern.search(fields[8])
            if not gene or stable(gene[1]) not in genes:
                continue
            transcript = transcript_pattern.search(fields[8])
            if transcript is None:
                continue
            result[stable(transcript[1])].append(
                {
                    "feature": fields[2],
                    "chromosome": fields[0],
                    "start": int(fields[3]),
                    "end": int(fields[4]),
                    "strand": fields[6],
                    "phase": fields[7],
                    "attributes": completeness.attrs(fields[8]),
                }
            )
    return dict(result)


def rank_evidence(annotation: dict, transcript: str, external_support: dict | None = None) -> dict:
    tags = annotation.get("tags", "").split(";")
    principal = [int(m[1]) for tag in tags if (m := re.fullmatch(r"appris_principal_(\d+)", tag))]
    alternative = [
        10 + int(m[1]) for tag in tags if (m := re.fullmatch(r"appris_alternative_(\d+)", tag))
    ]
    gtf_rank = min(principal or alternative or [99])
    external_support = external_support or {"validation_status": "no_external_source_pinned"}
    verified = external_support.get("validation_status") == "verified_identical_protein"
    external_rank = external_support.get("appris_rank", 99) if verified else 99
    appris_rank = gtf_rank if gtf_rank < 99 else external_rank
    source = (
        "Ensembl116_GTF_tags"
        if gtf_rank < 99
        else (
            "APPRIS116_exact_protein_assignment"
            if external_support.get("assignment_source_release") == 116
            else "APPRIS115_exact_protein_transfer"
        )
        if verified
        else "no_verified_APPRIS_support"
    )
    canonical = "Ensembl_canonical" in tags
    tsl = re.match(r"^(\d+)", annotation.get("transcript_support_level", ""))
    tsl_rank = int(tsl[1]) if tsl else 99
    quality = [appris_rank, 0 if canonical else 1, tsl_rank]
    return {
        "appris_source": source,
        "appris_support_source": source,
        "gtf116_appris_rank": gtf_rank,
        "external_appris_support": external_support,
        "appris_rank": appris_rank,
        "ensembl_canonical": canonical,
        "transcript_support_level": annotation.get("transcript_support_level", ""),
        "tsl_rank": tsl_rank,
        "quality_rank_vector": quality,
        "representative_rank_vector": [*quality, transcript],
        "transcript_id_tiebreak_scope": "equivalent encodings of the same species/gene/protein only",
    }


def attribute_values(features: list[dict], key: str) -> set[str]:
    return {value for feature in features for value in feature["attributes"].get(key, [])}


def validate_registry(value: dict) -> list[dict]:
    records = value.get("exceptions", [])
    if not isinstance(records, list) or any(not isinstance(r, dict) for r in records):
        raise ValueError("Exception registry must contain a list of records")
    seen = set()
    for record in records:
        identity = (record.get("species"), record.get("gene"), record.get("transcript_version"))
        required = ("exception_id", "cds_sha256", "protein_sha256", "rationale", "evidence")
        if (
            identity in seen
            or any(not record.get(k) for k in required)
            or record.get("decision") != "approved"
        ):
            raise ValueError(
                "Exception registry requires unique identities, approval, hashes and evidence"
            )
        if any(
            not re.fullmatch(r"[0-9a-f]{64}", str(record[k]))
            for k in ("cds_sha256", "protein_sha256")
        ):
            raise ValueError("Exception registry requires exact SHA256 checksums")
        if (
            identity[0] not in SPECIES
            or not identity[1]
            or not identity[2]
            or "." not in identity[2]
        ):
            raise ValueError("Exception registry requires exact species/gene/versioned transcript")
        if not set(record.get("waive_review_reasons", [])) <= WAIVABLE_REVIEW_REASONS:
            raise ValueError(
                "Identity/integrity or unsupported review waiver in exception registry"
            )
        seen.add(identity)
    return records


def evaluate(
    species: str,
    transcript: str,
    cds_models: list[dict],
    peptide_models: list[dict],
    features: list[dict],
    target_genes: set[str],
    exceptions: list[dict],
    legacy_source: dict | None = None,
) -> tuple[dict, dict | None]:
    reasons, integrity = [], []
    if len(cds_models) != 1:
        reasons.append("missing_CDS_model" if not cds_models else "ambiguous_CDS_models")
    if len(peptide_models) != 1:
        reasons.append(
            "missing_peptide_model" if not peptide_models else "ambiguous_peptide_models"
        )
    cds_model = cds_models[0] if len(cds_models) == 1 else None
    peptide_model = peptide_models[0] if len(peptide_models) == 1 else None
    annotation = completeness.annotation_evidence(features)
    gene_values = {stable(g) for g in attribute_values(features, "gene_id")}
    gene_values |= {stable(m["tags"].get("gene", "")) for m in cds_models + peptide_models}
    gene_values.discard("")
    gene = next(iter(gene_values)) if len(gene_values) == 1 else ""
    if not gene or gene not in target_genes:
        integrity.append("gene_identity_ambiguous_or_outside_target_pool")
    dna = cds_model["sequence"] if cds_model else ""
    peptide = peptide_model["sequence"].removesuffix("*") if peptide_model else ""
    tx_version = cds_model["transcript_version"] if cds_model else ""
    if cds_model:
        synthetic = {
            "candidate_id": transcript,
            "library": "pool",
            "set": "ortholog",
            "species": species,
            "key": f"pool:{species}:{gene}:{transcript}",
        }
        pinned = {"cds": dna, "tags": cds_model["tags"], "transcript_version": tx_version}
        comparison = completeness.audit_row(
            synthetic,
            {"cds": dna, "transcript": transcript},
            None,
            pinned,
            features,
            snapshot="Compara116_target_gene_pool",
            release=116,
            role="documented_Compara116_source",
            priority=0,
        )
        reasons += [r for r in comparison["review_flags"].split(";") if r]
        for issue in mechanical.inspect_sequence(dna)["issues"]:
            if issue == "non_atg_start":
                continue  # neutral initiation review is already explicit above
            integrity.append("CDS_integrity:" + issue)
    else:
        comparison = {"annotation_found": annotation.get("annotation_found", False)}
    if peptide_model and not any(
        key in peptide_model["tags"] for key in ("chromosome", "primary_assembly")
    ):
        reasons.append("outside_frozen_primary_assembly_peptide_universe")
    if cds_model and peptide_model:
        if tx_version != peptide_model["transcript_version"]:
            integrity.append("CDS_peptide_transcript_version_mismatch")
        if cds_model["tags"].get("gene") != peptide_model["tags"].get("gene"):
            integrity.append("CDS_peptide_gene_version_mismatch")
        if len(dna) % 3 == 0 and set(dna) <= set("ACGT"):
            from Bio.Seq import Seq

            translated = str(Seq(dna).translate()).removesuffix("*")
            if translated != peptide:
                integrity.append("CDS_peptide_protein_mismatch")
        protein_ids = attribute_values(features, "protein_id")
        protein_versions = attribute_values(features, "protein_version")
        expected = {pid + "." + version for pid in protein_ids for version in protein_versions}
        if expected != {peptide_model["model_id"]}:
            reasons.append("GTF_peptide_protein_identity_or_version_unconfirmed")
        gene_versions = attribute_values(features, "gene_version")
        expected_genes = {gene + "." + version for version in gene_versions}
        if expected_genes != {cds_model["tags"].get("gene")}:
            reasons.append("GTF_CDS_gene_version_unconfirmed")
    if not annotation.get("annotation_found"):
        reasons.append("unknown_GTF_annotation")
    elif annotation.get("transcript_biotype") != "protein_coding":
        reasons.append("requires_protein_coding_biotype_evidence")
    external_support = verified_legacy_support(species, gene, transcript, peptide, legacy_source)
    rank = rank_evidence(annotation, transcript, external_support)
    # GENCODE distinguishes unconfirmed transcript ends from coding boundaries.
    # UTR/transcript-end uncertainty remains visible without rejecting an otherwise
    # independently supported CDS. Coding NF and every other review rule still block.
    original_reasons = sorted(set(reasons))
    advisory_reasons = sorted(
        set(reasons) & {"mRNA_start_NF_not_confirmed", "mRNA_end_NF_not_confirmed"}
    )
    reasons = [reason for reason in reasons if reason not in advisory_reasons]
    applied, waived = "", []
    for exception in exceptions:
        if (exception["species"], exception["gene"], exception["transcript_version"]) != (
            species,
            gene,
            tx_version,
        ):
            continue
        if exception["cds_sha256"] != digest(dna) or exception["protein_sha256"] != digest(peptide):
            raise ValueError(f"Exception checksum no longer matches {species}/{transcript}")
        applied = exception["exception_id"]
        waived = sorted(set(exception["waive_review_reasons"]) & set(reasons))
        reasons = [r for r in reasons if r not in waived]
    integrity = sorted(set(integrity))
    reasons = sorted(set(reasons))
    eligible = not integrity and not reasons and cds_model is not None and peptide_model is not None
    record = {
        "species": species,
        "gene": gene,
        "transcript": transcript,
        "transcript_version": tx_version,
        "cds": dna,
        "protein": peptide,
        "protein_id_version": peptide_model["model_id"] if peptide_model else "",
        "cds_sha256": digest(dna) if dna else "",
        "protein_sha256": digest(peptide) if peptide else "",
        "eligibility_status": "eligible"
        if eligible
        else "excluded_source_integrity"
        if integrity
        else "review_required",
        "eligibility_basis": (
            "documented_exception" if applied else "canonical_ATG_supported_boundaries"
        )
        if eligible
        else "not_eligible_for_automatic_draw",
        "eligibility_reasons": ";".join([*integrity, *reasons]),
        "review_reasons": reasons,
        "advisory_reasons": advisory_reasons,
        "integrity_reasons": integrity,
        "original_review_reasons": original_reasons,
        "exception_id": applied,
        "waived_review_reasons": waived,
        "rank_evidence": rank,
        "annotation": annotation,
        "comparison_evidence": comparison,
        "chosen_encoding": False,
    }
    ledger = {
        k: v
        for k, v in record.items()
        if k not in {"cds", "protein", "annotation", "comparison_evidence"}
    }
    ledger.update(
        {
            "cds_nt": len(dna),
            "protein_length": len(peptide),
            "start_codon": dna[:3],
            "last_triplet": dna[-3:],
            "cds_model_count": len(cds_models),
            "peptide_model_count": len(peptide_models),
            "gene_candidates": ";".join(sorted(gene_values)),
            **{k: v for k, v in annotation.items()},
            "comparison_evidence": comparison,
        }
    )
    return ledger, record if eligible else None


def choose_encodings(records: list[dict], source_hashes: dict) -> tuple[dict, set[tuple]]:
    grouped = defaultdict(list)
    for record in records:
        grouped[(record["species"], record["gene"], record["protein_sha256"])].append(record)
    output, chosen = {species: {} for species in SPECIES}, set()
    for (species, gene, protein_hash), alternatives in sorted(grouped.items()):
        alternatives.sort(key=lambda r: r["rank_evidence"]["representative_rank_vector"])
        winner = dict(alternatives[0])
        winner["chosen_encoding"] = True
        winner["eligible_alternative_transcripts"] = [
            {
                "transcript": r["transcript"],
                "transcript_version": r["transcript_version"],
                "cds_sha256": r["cds_sha256"],
                "rank_evidence": r["rank_evidence"],
                "representative_reason": "selected_best_annotation_rank"
                if i == 0
                else "equivalent_protein_lower_annotation_rank_or_transcript_tie",
            }
            for i, r in enumerate(alternatives)
        ]
        winner["protein_ids"] = sorted({stable(r["protein_id_version"]) for r in alternatives})
        winner["protein_id_versions"] = sorted({r["protein_id_version"] for r in alternatives})
        winner["source_hashes"] = source_hashes[species]
        output[species].setdefault(protein_hash, {"by_gene": {}})["by_gene"][gene] = winner
        chosen.add((species, gene, winner["transcript"]))
    return output, chosen


def legacy_appris_rows(path: Path) -> list[dict]:
    result = []
    for line in path.read_text().splitlines():
        fields = line.split("\t")
        if len(fields) < 5 or fields[0] == "Gene name":
            continue
        for transcript in fields[2].split(";"):
            result.append(
                {"gene": stable(fields[1]), "transcript": stable(transcript), "label": fields[4]}
            )
    return result


def prepare_legacy_source(
    paths: dict, provenance: dict, targets: set[str], release: int = 115
) -> dict | None:
    if not paths:
        return None
    if set(paths) != {"cds", "gtf", "appris"}:
        raise ValueError("APPRIS support requires version-matched CDS, GTF and APPRIS together")
    assignments = [r for r in legacy_appris_rows(paths["appris"]) if r["gene"] in targets]
    by_identity = defaultdict(list)
    for row in assignments:
        by_identity[(row["gene"], row["transcript"])].append(row["label"])
    wanted = {r["transcript"] for r in assignments}
    cds = read_sequence_models(paths["cds"], targets, "cds")
    return {
        "release": release,
        "records": assignments,
        "assignments": dict(by_identity),
        "cds": {t: models for t, models in cds.items() if t in wanted},
        "gtf": completeness.gtf_selected(paths["gtf"], wanted),
        "provenance": provenance,
    }


def verified_legacy_support(
    species: str, gene: str, transcript: str, target_protein: str, source: dict | None
) -> dict:
    """Use APPRIS ranking support only after pinned CDS/GTF/protein validation."""
    if source is None:
        return {"validation_status": "no_external_source_pinned", "species": species}
    release = source.get("release", 115)
    labels = source["assignments"].get((gene, transcript), [])
    if not labels:
        return {
            "validation_status": "not_assigned_in_pinned_APPRIS_table",
            "assignment_source_release": release,
        }
    evidence = {
        "assignment_source_release": release,
        "assignment_labels": sorted(set(labels)),
        "source_provenance": source["provenance"],
        "gene": gene,
        "transcript": transcript,
        "target_116_protein_sha256": digest(target_protein) if target_protein else "",
        "interpretation": f"Protein-validated APPRIS release{release} support; "
        + (
            "version-matched116 assignment."
            if release == 116
            else "cross-release support, not a116 APPRIS assignment."
        ),
    }
    models = source["cds"].get(transcript, [])
    if len(models) != 1:
        return {**evidence, "validation_status": f"source{release}_CDS_missing_or_ambiguous"}
    model = models[0]
    fs = source["gtf"].get(transcript, [])
    annotation = completeness.annotation_evidence(fs)
    version = (
        model["transcript_version"].rsplit(".", 1)[1] if "." in model["transcript_version"] else ""
    )
    location = next(
        (
            model["tags"][key]
            for key in ("chromosome", "primary_assembly", "scaffold")
            if key in model["tags"]
        ),
        "",
    ).split(":")
    strand = {"1": "+", "-1": "-"}.get(location[4] if len(location) > 4 else "", "")
    genes = attribute_values(fs, "gene_id")
    gene_versions = attribute_values(fs, "gene_version")
    checks = {
        "stable_gene_coherent": genes == {gene} and stable(model["tags"].get("gene", "")) == gene,
        "transcript_version_coherent": bool(
            version and version == annotation.get("transcript_version")
        ),
        "gene_version_coherent": {gene + "." + v for v in gene_versions}
        == {model["tags"].get("gene")},
        "assembly_coherent": bool(location and location[0] == completeness.ASSEMBLIES.get(species)),
        "chromosome_coherent": bool(
            len(location) > 1 and location[1] == annotation.get("chromosome")
        ),
        "strand_coherent": bool(strand and strand == annotation.get("strand")),
    }
    evidence.update(
        source_transcript_version=model["transcript_version"],
        source_cds_sha256=digest(model["sequence"]),
        source_identity_checks=checks,
    )
    if not all(checks.values()):
        return {**evidence, "validation_status": f"source{release}_identity_or_version_unconfirmed"}
    dna = model["sequence"]
    if len(dna) % 3 or set(dna) - set("ACGT") or not dna:
        return {**evidence, "validation_status": f"source{release}_translation_unassessable"}
    from Bio.Seq import Seq

    protein = str(Seq(dna).translate()).removesuffix("*")
    evidence["source_protein_sha256"] = digest(protein)
    evidence[f"source_{release}_protein_sha256"] = digest(protein)
    if "*" in protein or not target_protein:
        return {**evidence, "validation_status": "source_or_target_protein_unassessable"}
    if protein != target_protein:
        return {**evidence, "validation_status": f"source{release}_to_target116_protein_mismatch"}
    ranks = []
    for label in labels:
        match = re.fullmatch(r"(PRINCIPAL|ALTERNATIVE):(\d+)", label)
        if match:
            ranks.append(int(match[2]) + (0 if match[1] == "PRINCIPAL" else 10))
    return {
        **evidence,
        "validation_status": "verified_identical_protein",
        "appris_rank": min(ranks or [99]),
    }


def summarize_appris(ledger: list[dict], targets: dict, legacy: dict) -> tuple[dict, list[dict]]:
    """Separate release116 GTF absence from independently pinned legacy coverage."""
    summary, gene_rows = {}, []
    for species in SPECIES:
        members = [r for r in ledger if r["species"] == species]
        by_identity = {(r["gene"], r["transcript"]): r for r in members}
        sources = [("Ensembl116_GTF_tags", True, None)]
        sources.append(("external_APPRIS", species in legacy, legacy.get(species)))
        summary[species] = {}
        for source, available, legacy_source in sources:
            assignments, principals, available_genes = set(), set(), set()
            source_principals = set()
            if source == "Ensembl116_GTF_tags":
                for row in members:
                    tags = row.get("tags", "").split(";")
                    identity = (row["gene"], row["transcript"])
                    if any(tag.startswith("appris_") for tag in tags):
                        assignments.add(identity)
                        available_genes.add(row["gene"])
                    if any(tag.startswith("appris_principal_") for tag in tags):
                        principals.add(identity)
                source_principals = set(principals)
                interpretation = "Assignments are tags in the pinned release116 GTF. Missing tags do not exclude a model."
            elif legacy_source:
                for row in legacy_source["records"]:
                    if row["gene"] not in targets[species]:
                        continue
                    identity = (row["gene"], row["transcript"])
                    available_genes.add(row["gene"])
                    if identity in by_identity:
                        assignments.add(identity)
                    if row["label"].startswith("PRINCIPAL"):
                        source_principals.add(identity)
                        if identity in by_identity:
                            principals.add(identity)
                stable_overlaps = set(principals)
                principals = {
                    identity
                    for identity in stable_overlaps
                    if by_identity[identity]["rank_evidence"]["external_appris_support"].get(
                        "validation_status"
                    )
                    == "verified_identical_protein"
                }
                interpretation = f"Pinned APPRIS release{legacy_source.get('release', 115)} assignments validated by matching-source CDS/GTF identity/version and exact protein equality to116 peptides; rank carries its exact source release."
            else:
                interpretation = "No external APPRIS table is pinned for this species; this is unavailable evidence, not evidence that APPRIS has no assignment."
            eligible = {
                identity
                for identity in principals
                if by_identity[identity]["eligibility_status"] == "eligible"
            }
            rejected = principals - eligible
            rejected_rows = [
                {
                    "gene": identity[0],
                    "transcript": identity[1],
                    "eligibility_reasons": by_identity[identity]["eligibility_reasons"],
                }
                for identity in sorted(rejected)
            ]
            entry = {
                "source_available": available,
                "interpretation": interpretation,
                "target_genes": len(targets[species]),
                "candidate_transcripts": len(members),
                "genes_with_any_source_assignment": len(available_genes),
                "candidate_transcripts_with_any_assignment": len(assignments),
                "genes_with_candidate_assignment": len({g for g, _ in assignments}),
                "candidate_transcripts_with_principal_assignment": len(principals),
                "genes_with_candidate_principal_assignment": len({g for g, _ in principals}),
                "eligible_principal_transcripts": len(eligible),
                "genes_with_eligible_principal": len({g for g, _ in eligible}),
                "rejected_principal_transcripts": len(rejected),
                "principal_assignments_absent_from_transcript_pool": len(
                    source_principals - set(by_identity)
                ),
                "genes_without_available_source_assignment": len(
                    targets[species] - available_genes
                ),
                "principal_rejection_reason_counts": dict(
                    sorted(
                        Counter(
                            reason
                            for row in rejected_rows
                            for reason in row["eligibility_reasons"].split(";")
                            if reason
                        ).items()
                    )
                ),
                "rejected_principal_records": rejected_rows,
            }
            if legacy_source:
                entry["source_provenance"] = legacy_source["provenance"]
                entry["assignment_source_release"] = legacy_source.get("release", 115)
                overlaps = source_principals & set(by_identity)
                entry["principal_stable_ID_overlaps_before_protein_validation"] = len(overlaps)
                entry["principal_validation_states"] = dict(
                    sorted(
                        Counter(
                            by_identity[identity]["rank_evidence"]["external_appris_support"][
                                "validation_status"
                            ]
                            for identity in overlaps
                        ).items()
                    )
                )
                entry["unverified_principal_overlap_records"] = [
                    {
                        "gene": identity[0],
                        "transcript": identity[1],
                        "validation_evidence": by_identity[identity]["rank_evidence"][
                            "external_appris_support"
                        ],
                    }
                    for identity in sorted(overlaps - principals)
                ]
            summary[species][source] = entry
            for gene in sorted(targets[species]):
                gene_principal = sorted(t for g, t in principals if g == gene)
                gene_eligible = sorted(t for g, t in eligible if g == gene)
                gene_rows.append(
                    {
                        "species": species,
                        "gene": gene,
                        "source": source,
                        "source_available": available,
                        "any_assignment_available": gene in available_genes,
                        "principal_candidate_transcripts": gene_principal,
                        "eligible_principal_transcripts": gene_eligible,
                        "principal_rejection_reasons": {
                            t: by_identity[(gene, t)]["eligibility_reasons"]
                            for t in gene_principal
                            if (gene, t) not in eligible
                        },
                        "interpretation": interpretation,
                    }
                )
    return summary, gene_rows


def write_table(path: Path, rows: list[dict]) -> None:
    columns = list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    k: json.dumps(v, sort_keys=True, separators=(",", ":"))
                    if isinstance(v, (dict, list))
                    else v
                    for k, v in row.items()
                }
            )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pool", type=Path, required=True)
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--exceptions", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    spec = json.loads(args.inputs.read_text())
    if spec.get("rule_version") != RULE_VERSION:
        raise ValueError("Input manifest eligibility rule version differs from this implementation")
    inputs, source_hashes = {}, {species: {} for species in SPECIES}
    legacy_paths, legacy_provenance, native_appris = {}, {}, {}
    for name, item in spec["files"].items():
        native = (
            item.get("release") == 116
            and item.get("species") in SPECIES
            and item.get("kind") in {"cds", "gtf", "pep", "appris"}
        )
        legacy = (
            item.get("release") == 115
            and item.get("species") == "mouse"
            and item.get("kind") in {"cds", "gtf", "appris"}
        )
        if not native and not legacy:
            continue
        path = args.input_dir / name
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
        if actual != item["sha256"]:
            raise ValueError(f"Pinned input checksum differs: {name}")
        if legacy:
            if item["kind"] in legacy_paths:
                raise ValueError("Duplicate source115 APPRIS input role")
            legacy_paths[item["kind"]] = path
            legacy_provenance[item["kind"]] = item
            continue
        if item["kind"] == "appris":
            native_appris[item["species"]] = {"path": path, "provenance": item}
            source_hashes[item["species"]]["appris"] = actual
            continue
        key = (item["species"], item["kind"])
        if key in inputs:
            raise ValueError(f"Duplicate annotation input role: {key}")
        inputs[key] = path
        source_hashes[item["species"]][item["kind"]] = actual
    registry = (
        json.loads(args.exceptions.read_text())
        if args.exceptions
        else spec.get("documented_exception_registry", {"exceptions": []})
    )
    exceptions = validate_registry(registry)
    with args.pool.open() as fh:
        pool_rows = list(csv.DictReader(fh, delimiter="\t"))
    targets = {species: set() for species in SPECIES}
    for row in pool_rows:
        if row["species"] not in SPECIES or str(row["compara_release"]) != "116":
            raise ValueError("Pool species or Compara release differs from this audit")
        if row["status"] == "eligible":
            genes = row["ortholog_ensembl_gene"].split(";")
            if len(genes) != 1 or not genes[0]:
                raise ValueError("Gene-level eligible row has ambiguous/missing target gene")
            targets[row["species"]].add(genes[0])
    legacy_sources = {}
    if legacy_paths:
        legacy_sources["mouse"] = prepare_legacy_source(
            legacy_paths, legacy_provenance, targets["mouse"]
        )
    for species, value in native_appris.items():
        paths = {kind: inputs[(species, kind)] for kind in ("cds", "gtf")}
        paths["appris"] = value["path"]
        native_provenance = {
            item["kind"]: item
            for item in spec["files"].values()
            if item.get("release") == 116
            and item.get("species") == species
            and item.get("kind") in {"cds", "gtf", "appris"}
        }
        legacy_sources[species] = prepare_legacy_source(
            paths, native_provenance, targets[species], release=116
        )
    ledger, eligible = [], []
    for species in SPECIES:
        if any((species, kind) not in inputs for kind in ("cds", "gtf", "pep")):
            raise ValueError(f"CDS/GTF/peptide inputs missing for {species}")
        cds = read_sequence_models(inputs[(species, "cds")], targets[species], "cds")
        peptides = read_sequence_models(inputs[(species, "pep")], targets[species], "pep")
        features = read_target_gtf(inputs[(species, "gtf")], targets[species])
        for transcript in sorted(cds.keys() | peptides.keys() | features.keys()):
            row, record = evaluate(
                species,
                transcript,
                cds.get(transcript, []),
                peptides.get(transcript, []),
                features.get(transcript, []),
                targets[species],
                exceptions,
                legacy_sources.get(species),
            )
            ledger.append(row)
            if record is not None:
                eligible.append(record)
    output, chosen = choose_encodings(eligible, source_hashes)
    for row in ledger:
        row["chosen_encoding"] = (row["species"], row["gene"], row["transcript"]) in chosen
    gene_rows = []
    transcripts_by_gene = defaultdict(list)
    for row in ledger:
        for gene in row["gene_candidates"].split(";"):
            if gene in targets[row["species"]]:
                transcripts_by_gene[(row["species"], gene)].append(row)
    for species in SPECIES:
        for gene in sorted(targets[species]):
            transcripts = transcripts_by_gene[(species, gene)]
            accepted = [r for r in transcripts if r["eligibility_status"] == "eligible"]
            gene_rows.append(
                {
                    "species": species,
                    "gene": gene,
                    "candidate_transcripts": len(transcripts),
                    "eligible_transcripts": len(accepted),
                    "eligible_proteins": len({r["protein_sha256"] for r in accepted}),
                    "automatic_draw_status": "has_eligible_CDS"
                    if accepted
                    else "no_eligible_CDS_requires_review",
                }
            )
    appris_summary, appris_genes = summarize_appris(ledger, targets, legacy_sources)
    args.out.mkdir(parents=True, exist_ok=False)
    for name, value in (
        ("eligible_proteins.json", output),
        ("documented_exception_registry.json", registry),
    ):
        path = args.out / name
        path.write_text(json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n")
        if name == "eligible_proteins.json":
            (args.out / (name + ".gz")).write_bytes(gzip.compress(path.read_bytes(), mtime=0))
    write_table(args.out / "ortholog_transcript_eligibility.tsv", ledger)
    transcript_path = args.out / "ortholog_transcript_eligibility.tsv"
    (args.out / "ortholog_transcript_eligibility.tsv.gz").write_bytes(
        gzip.compress(transcript_path.read_bytes(), mtime=0)
    )
    write_table(args.out / "ortholog_gene_eligibility.tsv", gene_rows)
    write_table(args.out / "input_gene_pool.tsv", pool_rows)
    write_table(args.out / "appris_gene_coverage.tsv", appris_genes)
    (args.out / "appris_coverage.json").write_text(json.dumps(appris_summary, indent=2) + "\n")
    commit = subprocess.run(
        ["git", "-C", str(Path(__file__).parent), "rev-parse", "HEAD"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    summary = {
        "rule_version": RULE_VERSION,
        "rules": RULES,
        "input_manifest": spec,
        "pool_sha256": hashlib.sha256(args.pool.read_bytes()).hexdigest(),
        "input_manifest_sha256": hashlib.sha256(args.inputs.read_bytes()).hexdigest(),
        "exceptions_sha256": hashlib.sha256(args.exceptions.read_bytes()).hexdigest()
        if args.exceptions
        else digest(json.dumps(registry, sort_keys=True)),
        "script_commit": commit,
        "source_dirty": bool(
            subprocess.run(
                ["git", "-C", str(Path(__file__).parent), "status", "--porcelain"],
                check=True,
                capture_output=True,
                text=True,
            ).stdout.strip()
        ),
        "source_hashes": source_hashes,
        "appris_coverage": appris_summary,
        "script_sha256": {
            p.name: hashlib.sha256(p.read_bytes()).hexdigest()
            for p in (Path(__file__), Path(completeness.__file__), Path(mechanical.__file__))
        },
        "command_argv": [Path(__file__).name, *sys.argv[1:]],
        "gene_pool_status_counts": dict(sorted(Counter(r["status"] for r in pool_rows).items())),
        "species": {
            species: {
                "target_genes": len(targets[species]),
                "transcript_status_counts": dict(
                    sorted(
                        Counter(
                            r["eligibility_status"] for r in ledger if r["species"] == species
                        ).items()
                    )
                ),
                "eligible_gene_proteins": sum(len(v["by_gene"]) for v in output[species].values()),
                "genes_without_eligible_CDS": sum(
                    r["species"] == species and not r["eligible_proteins"] for r in gene_rows
                ),
            }
            for species in SPECIES
        },
        "registered_exception_ids": sorted(r["exception_id"] for r in exceptions),
        "exceptions_applied": sorted({r["exception_id"] for r in ledger if r["exception_id"]}),
        "unused_exception_ids": sorted(
            {r["exception_id"] for r in exceptions}
            - {r["exception_id"] for r in ledger if r["exception_id"]}
        ),
        "activation": False,
        "sequences_changed": 0,
        "interpretation": "Eligibility is the recorded automatic-design rule, not evidence of biological functionality or inactivity.",
    }
    summary["output_sha256"] = {
        p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(args.out.iterdir())
    }
    (args.out / "ortholog_cds_eligibility_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps({"species": summary["species"], "transcripts": len(ledger)}, indent=2))


if __name__ == "__main__":
    main()
