"""Score the frozen ortholog pool only after exact gene/protein CDS eligibility.

Retains the correspondence metric and 0.001 tie window. Transcript quality is a
tie preference after reviewed/verified APPRIS and before current/accession. Missing
predictions remain pending candidates, never negative calls or silent exclusions.
No sequence repair, prediction submission, activation or publication occurs here.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import re
import sys
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import pandas as pd
import rebuild_ortholog_pool as pool
import score_ortholog_pool as legacy
import validated_v1_supplements as v1_supplements

from deliverome_analysis.cloning import seqtools
from deliverome_analysis.core.provenance import provenance
from deliverome_analysis.surfaceome import topology

FROZEN_TOPOLOGY_SHA256 = "604521dcbc920f520ea566aee3da638dad6540772c0f38d1792a631ca301e048"
SURFACE_LABELS = {"TM", "SP+TM", "BETA"}


def digest(sequence):
    return hashlib.sha256(sequence.encode()).hexdigest()


def file_digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_scoring_inputs(inputs, manifest):
    """Verify only used frozen inputs; the explicit topology replaces the old snapshot."""
    spec = json.loads(manifest.read_text())
    names = [
        "candidate_cds.json",
        "library_members.json",
        "gpi_set.tsv",
        "uniprot_mouse_reviewed_isoforms.tsv.gz",
        "uniprot_cyno_all_isoforms.tsv.gz",
    ]
    for species in pool.SPECIES:
        names.extend(
            [
                pool.PEPTIDES[species],
                f"{species}_ortholog_one2one_highconf_non_hla_metadata.csv",
                f"{species}_ortholog_one2one_highconf_non_hla.fasta",
            ]
        )
    for name in names:
        if file_digest(inputs / name) != spec["files"][name]["sha256"]:
            raise ValueError(f"Frozen scoring input mismatch: {name}")
    return spec, {name: spec["files"][name]["sha256"] for name in names}


def read_tsv(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, columns=None):
    fields = list(columns or [])
    for row in rows:
        fields.extend(k for k in row if k not in fields)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def eligible_encoding(registry, species, gene, sequence):
    """No accession, length, symbol, cross-gene or cross-species fallback."""
    record = registry.get(species, {}).get(digest(sequence), {}).get("by_gene", {}).get(gene)
    if not record or record.get("eligibility_status") != "eligible":
        return None
    if record["gene"] != gene or record["species"] != species or record["protein"] != sequence:
        raise ValueError("Eligibility registry gene/protein identity mismatch")
    if (
        not record["chosen_encoding"]
        or not record["transcript"]
        or record["transcript_version"].split(".")[0] != record["transcript"]
    ):
        raise ValueError("Eligibility registry representative transcript mismatch")
    cds = record["cds"]
    exception = record.get("eligibility_basis") == "documented_exception"
    if exception and (not record.get("exception_id") or not record.get("source_hashes")):
        raise ValueError("Documented CDS exception lacks explicit ID/source binding")
    if (
        (not cds.startswith("ATG") and not exception)
        or len(cds) % 3
        or set(cds) - set("ACGT")
        or cds[-3:] not in {"TAA", "TAG", "TGA"}
        or seqtools.translate(cds) != sequence + "*"
    ):
        raise ValueError("Eligibility registry encodes an invalid or different protein")
    if digest(cds) != record["cds_sha256"] or digest(sequence) != record["protein_sha256"]:
        raise ValueError("Eligibility registry sequence checksum mismatch")
    vector = record["rank_evidence"]["quality_rank_vector"]
    if len(vector) != 3 or any(not isinstance(v, int) for v in vector):
        raise ValueError("Eligibility record lacks its documented quality rank vector")
    if species == "cyno" and vector[0] != 99:
        raise ValueError("Frozen cyno registry must not claim unsupported APPRIS evidence")
    return record


def rejection_index(path):
    index = defaultdict(set)
    if path is None:
        return index
    for row in read_tsv(path):
        index[(row["species"], row["gene"], row["protein_sha256"])].update(
            r for r in row.get("eligibility_reasons", "").split(";") if r
        )
    return index


def load_anchor_cache(json_paths, text_paths):
    """Use only exact input sequence identity and explicit NetGPI1.1 provenance."""
    records = defaultdict(list)
    inventory = Counter()
    for path in json_paths:
        payload = json.loads(path.read_text())
        inventory["curation_records_not_predictions"] += len(payload.get("curation", []))
        for row in payload.get("predictions", []):
            inventory["prediction_rows"] += 1
            if row.get("tool_version") != "NetGPI-1.1":
                inventory["unsupported_tool_version"] += 1
                continue
            sequence, sha = row.get("sequence"), row.get("sequence_sha256")
            if not sequence or not sha:
                inventory["unresolved_sequence_not_reused"] += 1
                continue
            if row.get("sequence_resolution") != "exact_input":
                inventory["nonexact_resolution_not_reused"] += 1
                continue
            if digest(sequence) != sha or len(sequence) != int(row["sequence_length"]):
                raise ValueError("NetGPI export sequence identity mismatch")
            if not row.get("source_sha256") or row.get("prediction_number") is None:
                inventory["insufficient_provenance_not_reused"] += 1
                continue
            records[sha].append(
                {
                    k: row.get(k)
                    for k in [
                        "run_id",
                        "record_id",
                        "prediction_number",
                        "tool_version",
                        "source_uri",
                        "source_sha256",
                        "sequence_length",
                        "prediction",
                        "likelihood_text",
                        "omega_position",
                    ]
                }
            )
    for path in text_paths:
        text = path.read_text()
        source_digest = file_digest(path)
        if "NetGPI 1.1" not in text:
            raise ValueError("Raw NetGPI cache must declare version1.1")
        for line_number, line in enumerate(text.splitlines(), 1):
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 6 or len(fields[0]) != 64 or set(fields[0]) - set("0123456789abcdef"):
                raise ValueError("Raw NetGPI cache lacks exact submitted sequence digest")
            records[fields[0]].append(
                {
                    "tool_version": "NetGPI-1.1",
                    "source_file": path.name,
                    "source_sha256": source_digest,
                    "source_line_number": line_number,
                    "sequence_length": int(fields[1]),
                    "prediction": fields[2],
                    "likelihood_text": fields[4],
                    "omega_position": None if fields[3] == "-" else int(fields[3]),
                }
            )
    cache = {}
    for sha, entries in records.items():
        unique = {json.dumps(e, sort_keys=True): e for e in entries}
        entries = list(unique.values())
        lengths = {int(e["sequence_length"]) for e in entries}
        if len(lengths) != 1:
            raise ValueError("NetGPI same-digest length contradiction")
        calls = {
            True
            if e["prediction"] == "GPI-Anchored"
            else False
            if float(e["likelihood_text"]) >= 0.55
            else None
            for e in entries
        }
        sites = {e.get("omega_position") for e in entries if e["prediction"] == "GPI-Anchored"}
        conflict = len(calls) > 1 or len(sites) > 1
        call = next(iter(calls)) if not conflict else None
        cache[sha] = {
            "length": next(iter(lengths)),
            "call": call,
            "conflict": conflict,
            "state": "conflicting"
            if conflict
            else "ambiguous"
            if call is None
            else "positive"
            if call
            else "negative",
            "records": entries,
        }
    inventory["exact_sequence_digests"] = len(cache)
    inventory["conflicting_digests"] = sum(c["conflict"] for c in cache.values())
    return cache, dict(inventory)


def attach_predictions(row, candidate, labels, anchors):
    sha = digest(candidate["sequence"])
    candidate["label"] = labels.get(sha, "")
    anchor = anchors.get(sha) if row["human_is_gpi"] else None
    if anchor and anchor["length"] != len(candidate["sequence"]):
        raise ValueError("NetGPI cache length disagrees with exact selected sequence")
    candidate["anchor_evidence"] = anchor
    candidate["anchor_observed"] = anchor is not None
    candidate["anchor_ok"] = anchor["call"] if anchor else None
    # Unknown/contradictory anchor output is not a negative prediction.
    candidate["surface_ok"] = (
        candidate["anchor_ok"]
        if candidate["anchor_ok"] is not None
        else None
        if not candidate["label"]
        else candidate["label"] in SURFACE_LABELS or row["human_is_gpi"]
    )


def verified_principal(candidate):
    evidence = candidate["encoding"]["rank_evidence"]
    rank = evidence["quality_rank_vector"][0]
    source = evidence.get("appris_support_source", evidence.get("appris_source", ""))
    if rank == 99:
        return False
    if source == "APPRIS116_exact_protein_assignment":
        record = candidate["encoding"]
        proof = evidence.get("external_appris_support", {})
        identity_keys = {
            "assembly_coherent",
            "chromosome_coherent",
            "gene_version_coherent",
            "stable_gene_coherent",
            "strand_coherent",
            "transcript_version_coherent",
        }
        valid = (
            record["species"] == "mouse"
            and proof.get("validation_status") == "verified_identical_protein"
            and proof.get("assignment_source_release") == 116
            and proof.get("gene") == record["gene"]
            and proof.get("transcript") == record["transcript"]
            and proof.get("source_transcript_version") == record["transcript_version"]
            and proof.get("appris_rank") == rank == evidence.get("appris_rank")
            and proof.get("source_cds_sha256") == record["cds_sha256"]
            and all(
                proof.get(k) == record["protein_sha256"]
                for k in [
                    "source_protein_sha256",
                    "source_116_protein_sha256",
                    "target_116_protein_sha256",
                ]
            )
            and all(proof.get("source_identity_checks", {}).get(k) is True for k in identity_keys)
        )
        tiers = []
        for label in proof.get("assignment_labels", []):
            match = re.fullmatch(r"(PRINCIPAL|ALTERNATIVE):(\d+)", label)
            if match:
                tiers.append(int(match[2]) + (10 if match[1] == "ALTERNATIVE" else 0))
        valid = valid and min(tiers, default=99) == rank
        for kind in ["appris", "cds", "gtf"]:
            source_record = proof.get("source_provenance", {}).get(kind, {})
            valid = (
                valid
                and source_record.get("release") == 116
                and source_record.get("species") == "mouse"
                and len(source_record.get("sha256", "")) == 64
            )
            if kind in {"cds", "gtf"}:
                valid = valid and source_record.get("sha256") == record["source_hashes"].get(kind)
            elif record["source_hashes"].get("appris"):
                valid = valid and source_record.get("sha256") == record["source_hashes"]["appris"]
        if not valid:
            raise ValueError("APPRIS116 rank/protein/gene/transcript/source proof mismatch")
        return 1 <= rank <= 5
    tags = candidate["encoding"].get("annotation", {}).get("tags", "")
    if source == "Ensembl116_GTF_tags":
        tiers = [
            int(number) + (10 if kind == "alternative" else 0)
            for kind, number in re.findall(
                r"\bappris_(principal|alternative)_(\d+)\b", str(tags).lower()
            )
        ]
        if min(tiers, default=99) == rank:
            return 1 <= rank <= 5
    raise ValueError("Principal rank lacks verified APPRIS support")


def tie_key(alias, candidate, current):
    return (
        not candidate["reviewed"],
        not verified_principal(candidate),
        tuple(candidate["encoding"]["rank_evidence"]["quality_rank_vector"]),
        alias != current,
        alias,
    )


def score_job(job):
    row, reference, candidates, current = job
    scores = {
        sequence: legacy.score(reference, sequence)
        for sequence in sorted({c["sequence"] for c in candidates.values()})
    }
    allowed = {a: c for a, c in candidates.items() if c["surface_ok"] is not False}
    effective = allowed or candidates
    maximum = max(scores[c["sequence"]] for c in effective.values())
    tied = [a for a, c in effective.items() if maximum - scores[c["sequence"]] <= 0.001]
    picked = min(tied, key=lambda a: tie_key(a, candidates[a], current))
    candidate = candidates[picked]
    sequence = candidate["sequence"]
    needs_gpi = row["human_is_gpi"] and not candidate["anchor_observed"]
    uncertain = (
        row["human_is_gpi"] and candidate["anchor_observed"] and candidate["anchor_ok"] is None
    )
    state = (
        "needs_netgpi_and_v1"
        if needs_gpi and not candidate["label"]
        else "needs_netgpi"
        if needs_gpi
        else "needs_netgpi_review_and_v1"
        if uncertain and not candidate["label"]
        else "needs_netgpi_review"
        if uncertain
        else "needs_deeptmhmm_v1"
        if not candidate["label"]
        else "known_non_surface"
        if candidate["surface_ok"] is False
        else "provisional_surface_pick"
    )
    result = {
        **row,
        "picked_acc": picked,
        "picked_source": candidate["source"],
        "picked_reviewed": candidate["reviewed"],
        "picked_score": round(scores[sequence], 4),
        "picked_length": len(sequence),
        "picked_sequence_sha256": digest(sequence),
        "human_reference_sha256": digest(reference),
        "picked_topology_label": candidate["label"],
        "picked_topology_provenance": json.dumps(
            candidate.get("topology_provenance", {}), sort_keys=True
        ),
        "picked_anchor_call": {True: "anchored", False: "not_anchored"}.get(
            candidate["anchor_ok"], ""
        ),
        "topology_tool": topology.TOOL_V1,
        "anchor_uncertain": uncertain,
        "sampling_eligible": not needs_gpi
        and not uncertain
        and candidate["surface_ok"] is not False
        and bool(candidate["label"])
        and (candidate["label"] in SURFACE_LABELS or candidate["anchor_ok"] is True),
        "n_candidates": len(candidates),
        "n_tied": len(tied),
        "n_candidates_unknown_v1": sum(not c["label"] for c in candidates.values()),
        "selection_state": state,
        "sequence": sequence,
        "picked_cds_transcript": candidate["encoding"]["transcript"],
        "picked_cds_sha256": digest(candidate["encoding"]["cds"]),
        "picked_cds_quality_rank": json.dumps(
            candidate["encoding"]["rank_evidence"]["quality_rank_vector"]
        ),
        "picked_cds_rank_evidence": json.dumps(
            candidate["encoding"]["rank_evidence"], sort_keys=True
        ),
    }
    ledger = []
    chosen_key = tie_key(picked, candidate, current)
    names = [
        "reviewed",
        "verified_APPRIS_principal",
        "CDS_transcript_quality",
        "current_alias",
        "lexical_alias",
    ]
    for alias, value in sorted(candidates.items()):
        if alias == picked:
            reason = "provisional_winner"
        elif allowed and value["surface_ok"] is False:
            reason = "known_non_surface_guard"
        elif alias not in tied:
            reason = "lower_human_correspondence_outside_0.001_tie"
        else:
            losing_key = tie_key(alias, value, current)
            difference = next(
                (names[i] for i in range(len(names)) if losing_key[i] != chosen_key[i]),
                "identical_alias",
            )
            reason = "tie_preference_" + difference
        ledger.append(
            {
                "species": row["species"],
                "gene": row["gene"],
                "ortholog_ensembl_gene": row["ortholog_ensembl_gene"],
                "candidate_acc": alias,
                "candidate_source": value["source"],
                "protein_sha256": digest(value["sequence"]),
                "protein_length": len(value["sequence"]),
                "cds_eligibility": "eligible",
                "eligibility_reasons": "",
                "correspondence_score": scores[value["sequence"]],
                "reviewed": value["reviewed"],
                "verified_appris_principal": verified_principal(value),
                "representative_transcript": value["encoding"]["transcript"],
                "quality_rank_vector": json.dumps(
                    value["encoding"]["rank_evidence"]["quality_rank_vector"]
                ),
                "rank_evidence": json.dumps(value["encoding"]["rank_evidence"], sort_keys=True),
                "topology_label": value["label"],
                "topology_provenance": json.dumps(
                    value.get("topology_provenance", {}), sort_keys=True
                ),
                "anchor_state": value["anchor_evidence"]["state"]
                if value["anchor_evidence"]
                else "missing"
                if row["human_is_gpi"]
                else "not_required",
                "anchor_provenance": json.dumps(
                    value["anchor_evidence"]["records"] if value["anchor_evidence"] else [],
                    sort_keys=True,
                ),
                "winner_alias": alias == picked,
                "winner_protein": value["sequence"] == sequence,
                "decision_reason": reason,
            }
        )
    return result, ledger


def collect_jobs(inputs, pairs, registry, rejected, snapshot, anchors, supplements):
    cds = json.loads((inputs / "candidate_cds.json").read_text())
    genes = json.loads((inputs / "library_members.json").read_text())
    references = {
        g: seqtools.translate(cds[f"lib:{g}"]["cds"]).rstrip("*")
        for g in genes
        if f"lib:{g}" in cds
    }
    gpi = {
        r["gene_symbol"]
        for r in legacy.read_tsv(inputs / "gpi_set.tsv")
        if r["gpi_anchored"] == "1"
    }
    jobs, excluded, pair_rows, reconciliation = [], [], [], []
    for species in pool.SPECIES:
        proteins = pool.read_peptides(inputs / pool.PEPTIDES[species])
        targets = {
            r["ortholog_ensembl_gene"]
            for r in pairs
            if r["species"] == species and r["status"] == "eligible"
        }
        by_gene = defaultdict(dict)
        for p in proteins:
            if p["primary_assembly"] and p["gene"] in targets:
                by_gene[p["gene"]][p["protein"]] = {
                    "sequence": p["sequence"],
                    "source": "ensembl",
                    "reviewed": False,
                }
        up = defaultdict(dict)
        name = (
            "uniprot_mouse_reviewed_isoforms.tsv.gz"
            if species == "mouse"
            else "uniprot_cyno_all_isoforms.tsv.gz"
        )
        with gzip.open(inputs / name, "rt") as handle:
            for entry in csv.DictReader(handle, delimiter="\t"):
                symbol = entry.get("Gene Names (primary)", "").strip()
                if symbol:
                    up[symbol][entry["Entry"]] = {
                        "sequence": entry["Sequence"].rstrip("*"),
                        "source": "uniprot",
                        "reviewed": species == "mouse"
                        or entry.get("Reviewed", "").lower().startswith("reviewed"),
                    }
        old_meta = {
            r["query_ensembl_gene_id"]: r
            for r in csv.DictReader(
                (inputs / f"{species}_ortholog_one2one_highconf_non_hla_metadata.csv").open()
            )
            if r.get("status") == "ok"
        }
        cohort = legacy.parse_uniprot_fasta(
            inputs / f"{species}_ortholog_one2one_highconf_non_hla.fasta"
        )
        species_jobs = []
        for original in pairs:
            if original["species"] != species:
                continue
            row = {**original, "human_is_gpi": original["gene"] in gpi}
            if row["status"] != "eligible":
                pair_rows.append(
                    {
                        **row,
                        "pair_state": "upstream_pair_not_eligible",
                        "reason": row["status"],
                        "n_raw_candidates": 0,
                        "n_cds_eligible_aliases": 0,
                    }
                )
                continue
            target = row["ortholog_ensembl_gene"]
            candidates = dict(by_gene[target])
            for symbol in row["ortholog_symbols"].split(";"):
                candidates.update(up.get(symbol, {}))
            old = old_meta.get(row["human_ensembl_gene"], {})
            current = (
                old.get("selected_uniprot_accession", "")
                if old.get("target_ensembl_gene_id") == target
                else ""
            )
            if current in cohort:
                if (
                    current in candidates
                    and candidates[current]["sequence"] != cohort[current]["seq"]
                ):
                    prior_sequence = cohort[current]["seq"]
                    prior_encoding = eligible_encoding(registry, species, target, prior_sequence)
                    excluded.append(
                        {
                            "species": species,
                            "gene": row["gene"],
                            "ortholog_ensembl_gene": target,
                            "candidate_acc": current,
                            "candidate_source": "current_fasta_shadowed",
                            "protein_sha256": digest(prior_sequence),
                            "protein_length": len(prior_sequence),
                            "cds_eligibility": "eligible" if prior_encoding else "excluded",
                            "eligibility_reasons": ""
                            if prior_encoding
                            else "no_eligible_exact_gene_protein_encoding",
                            "decision_reason": "legacy_alias_precedence_retains_frozen_UniProt_sequence",
                        }
                    )
                candidates.setdefault(
                    current,
                    {
                        "sequence": cohort[current]["seq"],
                        "source": "current",
                        "reviewed": old.get("selected_reviewed") == "reviewed",
                    },
                )
            allowed = {}
            for alias, value in sorted(candidates.items()):
                encoding = eligible_encoding(registry, species, target, value["sequence"])
                if encoding:
                    allowed[alias] = {**value, "encoding": encoding}
                else:
                    reasons = rejected.get((species, target, digest(value["sequence"])), set())
                    excluded.append(
                        {
                            "species": species,
                            "gene": row["gene"],
                            "ortholog_ensembl_gene": target,
                            "candidate_acc": alias,
                            "candidate_source": value["source"],
                            "protein_sha256": digest(value["sequence"]),
                            "protein_length": len(value["sequence"]),
                            "cds_eligibility": "excluded",
                            "eligibility_reasons": ";".join(sorted(reasons))
                            or "no_eligible_exact_gene_protein_encoding",
                            "decision_reason": "excluded_before_correspondence_score",
                        }
                    )
            state = (
                "human_designed_cds_missing"
                if row["gene"] not in references
                else "no_cds_eligible_protein"
                if not allowed
                else "scored"
            )
            pair_rows.append(
                {
                    **row,
                    "pair_state": state,
                    "reason": state,
                    "n_raw_candidates": len(candidates),
                    "n_cds_eligible_aliases": len(allowed),
                }
            )
            if state != "scored":
                for alias, value in allowed.items():
                    excluded.append(
                        {
                            "species": species,
                            "gene": row["gene"],
                            "ortholog_ensembl_gene": target,
                            "candidate_acc": alias,
                            "candidate_source": value["source"],
                            "protein_sha256": digest(value["sequence"]),
                            "protein_length": len(value["sequence"]),
                            "cds_eligibility": "eligible",
                            "eligibility_reasons": "",
                            "decision_reason": state,
                        }
                    )
                continue
            species_jobs.append((row, references[row["gene"]], allowed, current))
        seqs = {
            digest(c["sequence"]): c["sequence"]
            for _, _, candidates, _ in species_jobs
            for c in candidates.values()
        }
        report = topology.resolve(
            seqs,
            tool=topology.TOOL_V1,
            cohorts=topology.cohorts_for_species(species),
            snapshot=snapshot,
        )
        labels, topological_evidence, comparisons = v1_supplements.overlay(
            species, report.to_dict("records"), supplements
        )
        reconciliation.extend(comparisons)
        for row, reference, candidates, current in species_jobs:
            for value in candidates.values():
                attach_predictions(row, value, labels, anchors)
                value["topology_provenance"] = topological_evidence[digest(value["sequence"])]
            jobs.append((row, reference, candidates, current))
    return jobs, excluded, pair_rows, reconciliation


def before_after(prior, selected, pairs, registry, candidate_ledger, rejected):
    old = {(r["species"], r["gene"]): r for r in prior}
    new = {(r["species"], r["gene"]): r for r in selected}
    pair_map = {(r["species"], r["gene"]): r for r in pairs}
    current_reasons = defaultdict(set)
    for candidate in candidate_ledger:
        current_reasons[(candidate["species"], candidate["gene"], candidate["protein_sha256"])].add(
            candidate["decision_reason"]
        )
    rows = []
    for key in sorted(old.keys() | new.keys() | pair_map.keys()):
        a, b, pair = old.get(key, {}), new.get(key, {}), pair_map.get(key, {})
        target = b.get(
            "ortholog_ensembl_gene",
            a.get("ortholog_ensembl_gene", pair.get("ortholog_ensembl_gene", "")),
        )
        old_sha = a.get("picked_sequence_sha256", "")
        old_eligible = (
            registry.get(key[0], {})
            .get(old_sha, {})
            .get("by_gene", {})
            .get(target, {})
            .get("eligibility_status")
            == "eligible"
        )
        change = (
            "same_protein"
            if a and b and old_sha == b["picked_sequence_sha256"]
            else "changed_protein"
            if a and b
            else "lost_winner"
            if a
            else "new_winner"
            if b
            else "no_winner"
        )
        reason = (
            "old_winner_failed_exact_CDS_gate"
            if a and not old_eligible
            else "correspondence_surface_or_quality_tie_policy"
            if change == "changed_protein"
            else pair.get("reason", change)
        )
        rows.append(
            {
                "species": key[0],
                "gene": key[1],
                "ortholog_ensembl_gene": target,
                "change": change,
                "reason": reason,
                "old_winner_CDS_eligible": old_eligible,
                "old_winner_CDS_exclusion_reasons": ";".join(
                    sorted(rejected.get((key[0], target, old_sha), set()))
                )
                if not old_eligible
                else "",
                "old_winner_new_ranking_reasons": ";".join(
                    sorted(current_reasons.get((key[0], key[1], old_sha), set()))
                ),
                "old_acc": a.get("picked_acc", ""),
                "old_protein_sha256": old_sha,
                "old_score": a.get("picked_score", ""),
                "old_state": a.get("selection_state", ""),
                "new_acc": b.get("picked_acc", ""),
                "new_protein_sha256": b.get("picked_sequence_sha256", ""),
                "new_score": b.get("picked_score", ""),
                "new_state": b.get("selection_state", pair.get("pair_state", "")),
                "new_cds_transcript": b.get("picked_cds_transcript", ""),
                "human_reference_changed": bool(
                    a and b and a.get("human_reference_sha256") != b.get("human_reference_sha256")
                ),
            }
        )
    return rows


def run(args):
    if args.out.exists():
        raise FileExistsError(f"Output directory already exists: {args.out}")
    spec, used_inputs = verify_scoring_inputs(args.inputs, args.input_manifest)
    built = json.loads((args.pool / "ortholog_pool_summary.json").read_text())
    for name, sha in built["input_sha256"].items():
        if sha != spec["files"][name]["sha256"]:
            raise ValueError(f"Frozen pool input mismatch: {name}")
    if file_digest(args.topology) != FROZEN_TOPOLOGY_SHA256:
        raise ValueError("Topology input is not the frozen71639-row snapshot")
    registry = json.loads(args.eligibility.read_text())
    supplements, v1_manifest = v1_supplements.load(args.v1_supplement, args.v1_run_record)
    rejected = rejection_index(args.eligibility_ledger)
    anchors, anchor_inventory = load_anchor_cache(args.netgpi_json, args.netgpi)
    pairs = legacy.read_tsv(args.pool / "ortholog_pool.tsv")
    jobs, exclusions, pair_rows, reconciliation = collect_jobs(
        args.inputs, pairs, registry, rejected, pd.read_parquet(args.topology), anchors, supplements
    )
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        scored = list(executor.map(score_job, jobs, chunksize=10))
    selected = [r for r, _ in scored]
    ledger = [*exclusions, *(entry for _, rows in scored for entry in rows)]
    raw_count = sum(r["n_raw_candidates"] for r in pair_rows)
    shadow_count = sum(
        r["decision_reason"] == "legacy_alias_precedence_retains_frozen_UniProt_sequence"
        for r in ledger
    )
    if len(ledger) != raw_count + shadow_count:
        raise ValueError("Candidate accounting does not reconcile every scored/excluded alias")
    args.out.mkdir(parents=True, exist_ok=False)
    write_tsv(
        args.out / "ortholog_provisional_picks.tsv",
        [{k: v for k, v in r.items() if k != "sequence"} for r in selected],
    )
    write_tsv(args.out / "ortholog_pair_ledger.tsv", pair_rows)
    write_tsv(
        args.out / "v1_supplement_reconciliation.tsv",
        reconciliation,
        columns=[
            "species",
            "sequence_sha256",
            "state",
            "same_topology",
            "frozen_label",
            "supplement_label",
        ],
    )
    write_tsv(
        args.out / "ortholog_before_after.tsv",
        before_after(
            legacy.read_tsv(args.prior_picks), selected, pair_rows, registry, ledger, rejected
        ),
    )
    columns = []
    for row in ledger:
        columns.extend(k for k in row if k not in columns)
    with (
        (args.out / "ortholog_candidate_ledger.tsv.gz").open("wb") as raw,
        gzip.GzipFile(filename="", fileobj=raw, mode="wb", mtime=0) as gz,
        io.TextIOWrapper(gz, newline="") as handle,
    ):
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(ledger)
    with (args.out / "ortholog_picks.fasta").open("w") as handle:
        for row in selected:
            handle.write(
                f">{row['species']}|{row['gene']}|{row['picked_acc']}\n{row['sequence']}\n"
            )
    pending_v1 = {
        r["picked_sequence_sha256"]: r["sequence"]
        for r in selected
        if not r["picked_topology_label"]
    }
    pending_gpi = {
        r["picked_sequence_sha256"]: r["sequence"]
        for r in selected
        if r["selection_state"] in {"needs_netgpi", "needs_netgpi_and_v1"}
    }
    all_unknown = {
        digest(c["sequence"]): c["sequence"]
        for _, _, candidates, _ in jobs
        for c in candidates.values()
        if not c["label"]
    }
    for name, work in [
        ("missing_deeptmhmm_v1", pending_v1),
        ("pending_netgpi", pending_gpi),
        ("all_unknown_candidate_sequences", all_unknown),
    ]:
        with (args.out / f"{name}.fasta").open("w") as handle:
            for sha, sequence in sorted(work.items()):
                handle.write(f">{sha}\n{sequence}\n")
    uncertain = [r for r in selected if r["anchor_uncertain"]]
    write_tsv(
        args.out / "uncertain_anchor_picks.tsv",
        [{k: v for k, v in r.items() if k != "sequence"} for r in uncertain],
        columns=["species", "gene", "picked_acc", "picked_sequence_sha256", "selection_state"],
    )
    missing_reference = sum(r["pair_state"] == "human_designed_cds_missing" for r in pair_rows)
    appris_coverage = {}
    for species, proteins in registry.items():
        encodings = [
            encoding for protein in proteins.values() for encoding in protein["by_gene"].values()
        ]
        appris_coverage[species] = {
            "eligible_gene_proteins": len(encodings),
            "verified_principal_gene_proteins": sum(
                verified_principal({"encoding": encoding}) for encoding in encodings
            ),
            "without_APPRIS_contribution_not_excluded": sum(
                encoding["rank_evidence"]["quality_rank_vector"][0] == 99 for encoding in encodings
            ),
            "support_sources": dict(
                Counter(
                    encoding["rank_evidence"].get(
                        "appris_support_source",
                        encoding["rank_evidence"].get("appris_source", "none"),
                    )
                    for encoding in encodings
                )
            ),
        }
    summary = {
        "scored_pairs": len(selected),
        "all_input_pairs": len(pair_rows),
        "raw_candidate_aliases": raw_count,
        "additional_shadowed_legacy_aliases": shadow_count,
        "candidate_ledger_rows": len(ledger),
        "pair_states": dict(Counter(r["pair_state"] for r in pair_rows)),
        "candidate_decision_reasons": dict(Counter(r["decision_reason"] for r in ledger)),
        "missing_human_references": missing_reference,
        "netgpi_cache_inventory": anchor_inventory,
        "validated_v1_supplements": v1_manifest,
        "v1_supplement_reconciliation": dict(Counter(r["state"] for r in reconciliation)),
        "APPRIS_coverage": appris_coverage,
        "state_by_species": {
            s: dict(Counter(r["selection_state"] for r in selected if r["species"] == s))
            for s in pool.SPECIES
        },
        "first_batch_deeptmhmm_v1_sequences": len(pending_v1),
        "first_batch_deeptmhmm_v1_residues": sum(map(len, pending_v1.values())),
        "first_batch_netgpi_sequences": len(pending_gpi),
        "uncertain_anchor_picks": len(uncertain),
        "unknown_candidate_sequences": len(all_unknown),
        "eligible_for_sampling": dict(
            Counter(r["species"] for r in selected if r["sampling_eligible"])
        ),
        "panel_ready": not pending_v1 and not pending_gpi,
        "panel_ready_scope": "No missing winner predictor jobs. Known ambiguous/conflicting anchors and missing human references remain explicit sampling exclusions, not biological approvals.",
        "upstream_eligible_pairs": sum(r["status"] == "eligible" for r in pair_rows),
        "excluded_no_eligible_CDS_pairs": sum(
            r["pair_state"] == "no_cds_eligible_protein" for r in pair_rows
        ),
        "unresolved_known_anchor_sampling_exclusions": len(uncertain),
        "known_non_surface_sampling_exclusions": sum(
            r["selection_state"] == "known_non_surface" for r in selected
        ),
        "policy": "Exact species/gene/protein CDS gate; matching blocks divided by longer length; 0.001 tie: reviewed, verified Ensembl116-matched APPRIS principal, quality_rank_vector, current, accession. No stable-ID-only or silent115 APPRIS fallback. Cyno has no APPRIS contribution or eligibility penalty. Unknown predictions remain pending; no activation.",
        "input_sha256": {
            str(p): file_digest(p)
            for p in [
                args.eligibility,
                args.topology,
                args.prior_picks,
                args.input_manifest,
                *args.v1_supplement,
                *args.v1_run_record,
                *(args.netgpi_json),
                *(args.netgpi),
                *([args.eligibility_ledger] if args.eligibility_ledger else []),
            ]
        },
        "frozen_scoring_input_sha256": used_inputs,
    }
    (args.out / "ortholog_scoring_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    stamp = provenance()
    run_record = {
        "provenance": stamp,
        "argv": sys.argv,
        "script_sha256": file_digest(Path(__file__)),
        "outputs": [
            {"name": p.name, "sha256": file_digest(p)}
            for p in sorted(args.out.iterdir())
            if p.is_file()
        ],
    }
    (args.out / "scoring_run.json").write_text(json.dumps(run_record, indent=2) + "\n")
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--input-manifest", type=Path, default=pool.MANIFEST)
    parser.add_argument("--pool", type=Path, required=True)
    parser.add_argument("--eligibility", type=Path, required=True)
    parser.add_argument("--eligibility-ledger", type=Path)
    parser.add_argument("--topology", type=Path, required=True)
    parser.add_argument("--prior-picks", type=Path, required=True)
    parser.add_argument("--netgpi-json", type=Path, action="append", default=[])
    parser.add_argument("--netgpi", type=Path, action="append", default=[])
    parser.add_argument("--v1-supplement", type=Path, action="append", default=[])
    parser.add_argument("--v1-run-record", type=Path, action="append", default=[])
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(run(args), indent=2))


if __name__ == "__main__":
    main()
