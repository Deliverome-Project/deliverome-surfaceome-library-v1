import importlib.util
from pathlib import Path

import pytest

spec = importlib.util.spec_from_file_location(
    "replay_accessions", Path(__file__).with_name("replay_ortholog_accessions.py")
)
replay = importlib.util.module_from_spec(spec)
spec.loader.exec_module(replay)


def entry(accession, primary="GENE", names="GENE", reviewed="unreviewed", length="100"):
    return {
        "Entry": accession,
        "Gene Names (primary)": primary,
        "Gene Names": names,
        "Reviewed": reviewed,
        "Length": length,
        "Organism (ID)": "10090",
    }


def recorded(accession="A"):
    return {
        "species": "mouse",
        "taxid": "10090",
        "target_gene_symbol": "Gene",
        "query_ensembl_gene_id": "ENSG1",
        "target_ensembl_gene_id": "ENSMUSG1",
        "selected_uniprot_accession": accession,
        "selected_length": "100",
        "n_candidates": "1",
        "status": "ok",
    }


def test_historical_rank_prefers_primary_then_reviewed_then_length_then_accession():
    alias = entry("Z", primary="OTHER", names="OTHER Gene", reviewed="reviewed", length="999")
    primary = entry("A")
    assert replay.rank(primary, " Gene ") > replay.rank(alias, " Gene ")
    reviewed = entry("A", reviewed="reviewed", length="50")
    assert replay.rank(reviewed, "GENE") > replay.rank(primary, "GENE")
    assert replay.rank(entry("A", length="101"), "GENE") > replay.rank(primary, "GENE")
    assert replay.rank(entry("B"), "GENE") > replay.rank(primary, "GENE")


def test_alias_candidates_are_exact_case_insensitive_tokens():
    rows = [entry("A", "OTHER", "OTHER gene"), entry("B", "GENE2", "GENE2", length="999")]
    result = replay.audit([recorded()], rows, "10090")[0]
    assert result["outcome"] == "confirmed_accession"
    assert result["release_candidate_count"] == 1


def test_missing_and_changed_accessions_remain_distinct():
    changed = replay.audit([recorded()], [entry("B")], "10090")[0]
    missing = replay.audit([recorded()], [entry("B", "OTHER", "OTHER")], "10090")[0]
    assert changed["outcome"] == "different_accession"
    assert missing["outcome"] == "no_candidate"
    assert not changed["recorded_accession_present_in_release"]


def test_historical_search_cap_is_not_claimed_as_reproduced():
    rows = [entry(str(i)) for i in range(501)]
    assert (
        replay.audit([recorded("500")], rows, "10090")[0]["outcome"] == "historical_cap_unresolved"
    )


def test_species_contamination_and_duplicate_entries_fail_closed():
    with pytest.raises(ValueError, match="organism"):
        replay.audit([], [entry("A")], "9541")
    with pytest.raises(ValueError, match="Duplicate"):
        replay.audit([], [entry("A"), entry("A")], "10090")


def test_accession_confirmation_does_not_hide_length_change():
    result = replay.audit([recorded()], [entry("A", length="101")], "10090")[0]
    assert result["outcome"] == "confirmed_accession"
    assert not result["same_recorded_length"]
