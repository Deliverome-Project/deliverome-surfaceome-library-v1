import importlib.util
import json
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
spec = importlib.util.spec_from_file_location(
    "eligible_scorer", HERE / "score_eligible_ortholog_pool.py"
)
scorer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scorer)


def encoding(protein="MA", cds="ATGGCTTAA", gene="G1", species="mouse", quality=None):
    return {
        "gene": gene,
        "species": species,
        "protein": protein,
        "cds": cds,
        "eligibility_status": "eligible",
        "eligibility_basis": "canonical_ATG_supported_boundaries",
        "transcript": "T1",
        "transcript_version": "T1.1",
        "chosen_encoding": True,
        "cds_sha256": scorer.digest(cds),
        "protein_sha256": scorer.digest(protein),
        "source_hashes": {k: scorer.digest("pinned_" + k) for k in ["cds", "gtf", "pep", "appris"]},
        "rank_evidence": {
            "quality_rank_vector": quality or [99, 0, 99],
            "appris_source": "Ensembl116_GTF_tags",
        },
    }


def candidate(protein="MA", cds="ATGGCTTAA", label="TM", quality=None, reviewed=False):
    value = {
        "sequence": protein,
        "source": "ensembl",
        "reviewed": reviewed,
        "encoding": encoding(protein, cds, quality=quality),
    }
    scorer.attach_predictions({"human_is_gpi": False}, value, {scorer.digest(protein): label}, {})
    return value


def add_principal_proof(value):
    record = value["encoding"]
    rank = record["rank_evidence"]["quality_rank_vector"][0]
    record["rank_evidence"].update(
        appris_rank=rank, appris_support_source="APPRIS116_exact_protein_assignment"
    )
    record["rank_evidence"]["external_appris_support"] = {
        "validation_status": "verified_identical_protein",
        "assignment_source_release": 116,
        "gene": record["gene"],
        "transcript": record["transcript"],
        "source_transcript_version": record["transcript_version"],
        "appris_rank": rank,
        "assignment_labels": [f"PRINCIPAL:{rank}"],
        "source_cds_sha256": record["cds_sha256"],
        **{
            k: record["protein_sha256"]
            for k in [
                "source_protein_sha256",
                "source_116_protein_sha256",
                "target_116_protein_sha256",
            ]
        },
        "source_identity_checks": {
            k: True
            for k in [
                "assembly_coherent",
                "chromosome_coherent",
                "gene_version_coherent",
                "stable_gene_coherent",
                "strand_coherent",
                "transcript_version_coherent",
            ]
        },
        "source_provenance": {
            k: {"release": 116, "species": "mouse", "sha256": record["source_hashes"][k]}
            for k in ["appris", "cds", "gtf"]
        },
    }


def row(gpi=False):
    return {"species": "mouse", "gene": "HUMAN", "ortholog_ensembl_gene": "G1", "human_is_gpi": gpi}


def test_gate_does_not_lift_identical_protein_across_gene_or_species():
    registry = {"mouse": {scorer.digest("MA"): {"by_gene": {"G1": encoding()}}}}
    assert scorer.eligible_encoding(registry, "mouse", "G1", "MA") is not None
    assert scorer.eligible_encoding(registry, "mouse", "OTHER", "MA") is None
    assert scorer.eligible_encoding(registry, "cyno", "G1", "MA") is None
    assert scorer.eligible_encoding(registry, "mouse", "G1", "MV") is None


@pytest.mark.parametrize("change", ["cds", "hash", "species", "chosen", "transcript"])
def test_registry_corruption_fails_closed(change):
    value = encoding()
    if change == "cds":
        value["cds"] = "ATGGTTTAA"
    elif change == "hash":
        value["cds_sha256"] = "wrong"
    elif change == "species":
        value["species"] = "cyno"
    elif change == "chosen":
        value["chosen_encoding"] = False
    else:
        value["transcript_version"] = "OTHER.1"
    with pytest.raises(ValueError):
        scorer.eligible_encoding(
            {"mouse": {scorer.digest("MA"): {"by_gene": {"G1": value}}}}, "mouse", "G1", "MA"
        )


def test_documented_checksum_bound_non_atg_exception_is_not_silently_rejected():
    value = encoding("VA", "GTGGCTTAA")
    value.update(eligibility_basis="documented_exception", exception_id="reviewed-case-1")
    registry = {"mouse": {scorer.digest("VA"): {"by_gene": {"G1": value}}}}
    assert scorer.eligible_encoding(registry, "mouse", "G1", "VA") is value
    value["exception_id"] = ""
    with pytest.raises(ValueError, match="explicit ID"):
        scorer.eligible_encoding(registry, "mouse", "G1", "VA")


def test_unknown_topology_winner_remains_pending_instead_of_being_excluded(monkeypatch):
    unknown = candidate("MA", label="")
    known = candidate("MV", "ATGGTTTAA")
    monkeypatch.setattr(
        scorer.legacy, "score", lambda _reference, protein: {"MA": 0.95, "MV": 0.85}[protein]
    )
    selected, ledger = scorer.score_job(
        (row(), "MA", {"UNKNOWN": unknown, "KNOWN": known}, "KNOWN")
    )
    assert selected["picked_acc"] == "UNKNOWN"
    assert selected["selection_state"] == "needs_deeptmhmm_v1"
    assert selected["sampling_eligible"] is False
    assert len(ledger) == 2


def test_quality_tie_precedes_current_alias_but_follows_reviewed(monkeypatch):
    best_quality = candidate("MA", quality=[99, 0, 1])
    current = candidate("MV", "ATGGTTTAA", quality=[99, 1, 99])
    monkeypatch.setattr(
        scorer.legacy, "score", lambda _reference, protein: {"MA": 0.9, "MV": 0.9004}[protein]
    )
    selected, ledger = scorer.score_job(
        (row(), "MA", {"QUALITY": best_quality, "CURRENT": current}, "CURRENT")
    )
    assert selected["picked_acc"] == "QUALITY"
    assert (
        next(r for r in ledger if r["candidate_acc"] == "CURRENT")["decision_reason"]
        == "tie_preference_CDS_transcript_quality"
    )
    current["reviewed"] = True
    assert (
        scorer.score_job((row(), "MA", {"QUALITY": best_quality, "CURRENT": current}, "CURRENT"))[
            0
        ]["picked_acc"]
        == "CURRENT"
    )


def test_only_verified_appris_transfer_is_principal_evidence(monkeypatch):
    principal = candidate("MA", quality=[1, 1, 99])
    add_principal_proof(principal)
    canonical = candidate("MV", "ATGGTTTAA", quality=[99, 0, 1])
    monkeypatch.setattr(scorer.legacy, "score", lambda _reference, _protein: 0.9)
    assert (
        scorer.score_job(
            (row(), "MA", {"PRINCIPAL": principal, "CANONICAL": canonical}, "CANONICAL")
        )[0]["picked_acc"]
        == "PRINCIPAL"
    )
    principal["encoding"]["rank_evidence"]["appris_support_source"] = "stable_id_only_unverified"
    with pytest.raises(ValueError, match="verified APPRIS"):
        scorer.verified_principal(principal)


@pytest.mark.parametrize(
    "field,bad",
    [
        ("validation_status", "unverified"),
        ("assignment_source_release", 115),
        ("gene", "OTHER"),
        ("source_transcript_version", "T1.9"),
        ("target_116_protein_sha256", "wrong"),
        ("assignment_labels", ["PRINCIPAL:5"]),
    ],
)
def test_stale_appris_proof_cannot_contribute_to_ties(field, bad):
    value = candidate("MA", quality=[1, 1, 99])
    add_principal_proof(value)
    value["encoding"]["rank_evidence"]["external_appris_support"][field] = bad
    with pytest.raises(ValueError, match="proof mismatch"):
        scorer.verified_principal(value)


def test_netgpi_reuse_requires_sequence_and_preserves_contradictions(tmp_path):
    base = {
        "run_id": "run",
        "record_id": "one",
        "prediction_number": 1,
        "tool_version": "NetGPI-1.1",
        "sequence": "MA",
        "sequence_sha256": scorer.digest("MA"),
        "sequence_resolution": "exact_input",
        "sequence_length": 2,
        "source_sha256": "source",
        "prediction": "GPI-Anchored",
        "likelihood_text": "0.9",
        "omega_position": 1,
    }
    negative = {
        **base,
        "record_id": "two",
        "prediction": "Not GPI-Anchored",
        "omega_position": None,
    }
    unresolved = {
        **base,
        "record_id": "lost-input",
        "sequence": None,
        "sequence_sha256": None,
        "sequence_resolution": "unresolved_legacy",
    }
    path = tmp_path / "records.json"
    path.write_text(
        json.dumps(
            {"predictions": [base, negative, unresolved], "curation": [{"source_kind": "ad_hoc"}]}
        )
    )
    cache, inventory = scorer.load_anchor_cache([path], [])
    assert cache[scorer.digest("MA")]["state"] == "conflicting"
    assert cache[scorer.digest("MA")]["call"] is None
    assert len(cache[scorer.digest("MA")]["records"]) == 2
    assert inventory["unresolved_sequence_not_reused"] == 1
    assert inventory["curation_records_not_predictions"] == 1


def test_ambiguous_existing_anchor_requires_review_not_false_readiness():
    value = candidate("MA", label="SP")
    anchor = {"length": 2, "call": None, "state": "ambiguous", "records": []}
    scorer.attach_predictions(
        row(True), value, {scorer.digest("MA"): "SP"}, {scorer.digest("MA"): anchor}
    )
    selected, _ = scorer.score_job((row(True), "MA", {"P": value}, ""))
    assert selected["selection_state"] == "needs_netgpi_review"
    assert selected["anchor_uncertain"] is True
    assert selected["sampling_eligible"] is False
