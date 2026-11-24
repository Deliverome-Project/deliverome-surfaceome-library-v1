import importlib.util
from pathlib import Path

import pytest

spec = importlib.util.spec_from_file_location(
    "exact_compara_cds", Path(__file__).with_name("resolve_compara116_cds.py")
)
resolver = importlib.util.module_from_spec(spec)
spec.loader.exec_module(resolver)


def fixture():
    panel = [
        {
            "form_id": "mouse:GENE:P1",
            "species": "mouse",
            "human_symbol": "GENE",
            "ortholog_ensembl_gene": "ENSG1",
            "ensembl_transcripts": "ENST1",
            "sequence_sha256": resolver.sha256("MAV"),
            "existing_candidate_ids": "old1",
        }
    ]
    members = [{"candidate_id": "old1", "key": "ort:mouse:GENE:OLD"}]
    existing = {
        "ort:mouse:GENE:OLD": {
            "species": "mouse",
            "transcript": "OLDT",
            "cds": "ATGGCTGTTTAA",
            "origin": "ortholog",
            "route": "named transcript",
            "variant_residues": 3,
        }
    }
    peptides = {
        "mouse": {
            "ENST1": {
                "gene": "ENSG1",
                "protein": "MAV",
                "protein_id": "ENSP1.2",
                "transcript_version": "ENST1.2",
            }
        }
    }
    cds = {
        "mouse": {
            "ENST1": {"gene": "ENSG1", "cds": "ATGGCTGTTTAG", "transcript_version": "ENST1.2"}
        }
    }
    return panel, existing, members, peptides, cds


def test_exact_reuse_preserves_existing_record_and_legacy_variant_metadata():
    panel, existing, members, peptides, cds = fixture()
    output, rows = resolver.resolve(panel, existing, members, peptides, cds)
    assert output[panel[0]["form_id"]] == existing["ort:mouse:GENE:OLD"]
    assert rows[0]["translated_protein_sha256"] == panel[0]["sequence_sha256"]
    assert rows[0]["route"] == "exact existing CDS reuse"


def test_equal_length_variant_cannot_reuse_existing_cds():
    panel, existing, members, peptides, cds = fixture()
    existing["ort:mouse:GENE:OLD"]["cds"] = "ATGGCTACTTAA"
    output, rows = resolver.resolve(panel, existing, members, peptides, cds)
    assert output[panel[0]["form_id"]]["cds"] == "ATGGCTGTTTAG"
    assert rows[0]["route"] == "release116 exact CDS"
    assert output[panel[0]["form_id"]]["variant_residues"] == 0


@pytest.mark.parametrize(
    "dna", ["ATGGCTGTTT", "ATGGCNGTTTAA", "ATGTGAGTTTAA", "ATGGCTGTTTAATAA", "ATGGCTACTTAA"]
)
def test_invalid_or_nonmatching_dna_is_rejected(dna):
    assert resolver.verify_cds(dna, "MAV") is None


def test_missing_stop_is_appended_explicitly_without_claiming_native_stop():
    checked = resolver.verify_cds("ATGGCTGTT", "MAV")
    assert checked["cds"] == "ATGGCTGTTTAA"
    assert checked["native_stop"] == ""
    assert checked["stop_added"] == "TAA"
    assert resolver.verify_cds("ATGGCTGTT", "MAV", allow_added_stop=False) is None


def test_exact_non_atg_translation_never_gets_an_invented_methionine():
    checked = resolver.verify_cds("GCTGTTTAA", "AV")
    assert checked["cds"] == "GCTGTTTAA"
    assert resolver.verify_cds("GCTGTTTAA", "MAV") is None


def test_inherited_supplied_stop_is_not_called_native():
    panel, existing, members, peptides, cds = fixture()
    existing["ort:mouse:GENE:OLD"]["stop_added"] = "TAA"
    _, rows = resolver.resolve(panel, existing, members, peptides, cds)
    assert rows[0]["native_terminal_stop"] == ""
    assert rows[0]["stop_added"] == "TAA"
    assert rows[0]["stop_added_inherited"] is True


@pytest.mark.parametrize("field,value", [("gene", "OTHERGENE"), ("transcript_version", "ENST1.3")])
def test_new_cds_requires_selected_gene_and_version(field, value):
    panel, _, members, peptides, cds = fixture()
    cds["mouse"]["ENST1"][field] = value
    with pytest.raises(ValueError, match="No exact"):
        resolver.resolve(panel, {}, members, peptides, cds)


def test_selected_peptide_must_match_panel_hash_before_reuse():
    panel, existing, members, peptides, cds = fixture()
    peptides["mouse"]["ENST1"]["protein"] = "MAT"
    with pytest.raises(ValueError, match="exact target protein/gene"):
        resolver.resolve(panel, existing, members, peptides, cds)


def test_official_bsd_checksum_and_block_count_are_verified(tmp_path):
    data = tmp_path / "input.fa.gz"
    data.write_bytes(b"abc")
    checksum = tmp_path / "CHECKSUMS"
    checksum.write_text("16556 1 input.fa.gz\n")
    resolver.verify_published_checksum(data, checksum)
    data.write_bytes(b"abd")
    with pytest.raises(ValueError, match="checksum mismatch"):
        resolver.verify_published_checksum(data, checksum)
