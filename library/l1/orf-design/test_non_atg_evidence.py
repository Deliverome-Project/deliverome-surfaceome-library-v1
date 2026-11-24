import importlib.util
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).parent


def load(name):
    spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


load("resolve_compara116_cds")
features = load("audit_compara116_non_atg")
audit = load("audit_library_non_atg")
alternatives = load("audit_non_atg_alternatives")
priority = load("review_priority_cds_forms")


def annotation(*, plus=True, start=True, nf=False):
    strand = "+" if plus else "-"
    transcript = {
        "feature": "transcript",
        "start": 100,
        "end": 200,
        "strand": strand,
        "attributes": {"tag": ["cds_start_NF"] if nf else []},
    }
    cds = {"feature": "CDS", "start": 100, "end": 200, "strand": strand, "attributes": {}}
    codon = {
        "feature": "start_codon",
        "start": 100 if plus else 198,
        "end": 102 if plus else 200,
        "strand": strand,
        "attributes": {},
    }
    return features.transcript_evidence([transcript, cds, *([codon] if start else [])])


def test_repeated_gtf_tags_are_retained():
    assert features.attributes('tag "cds_start_NF"; tag "Ensembl_canonical";')["tag"] == [
        "cds_start_NF",
        "Ensembl_canonical",
    ]


@pytest.mark.parametrize("plus,coordinates", [(True, [100, 101, 102]), (False, [200, 199, 198])])
def test_start_feature_must_match_cds_first_codon_in_strand_order(plus, coordinates):
    a = annotation(plus=plus)
    assert a["first_three_CDS_genomic_positions"] == coordinates
    assert a["start_codon_at_CDS_start"] is True
    a["start_codon_at_CDS_start"] = False
    assert audit.classify(a) == "start_feature_not_at_CDS_start"


def test_missing_nf_tag_does_not_establish_completeness():
    a = annotation(start=False)
    assert audit.classify(a) == "missing_start_at_transcript_boundary"
    assert audit.classify(None) == "annotation_unknown"
    assert audit.classify(annotation(nf=True)) == "annotated_5prime_incomplete"


def test_documented_non_aug_is_not_treated_as_an_extra_m_prepend():
    a = annotation(start=False)
    assert audit.classify(a, documented_non_aug=True) == "documented_non_AUG_initiation"


def test_exact_candidate_relationships_do_not_use_length_heuristics():
    assert audit.protein_relationship("WQQA", "MWQQA") == "exact_M_plus_target"
    assert audit.protein_relationship("TWQQA", "MWQQA") == "same_tail_first_residue_differs"
    assert audit.protein_relationship("WQQA", "MAGWQQA") == "target_is_exact_suffix"
    assert audit.protein_relationship("WQQA", "MWQEA") == "other_protein_difference"


def test_same_gene_candidate_dna_must_be_a_usable_atg_orf():
    assert audit.usable_atg_model("ATGGCTTAA") == "MA"
    assert audit.usable_atg_model("ATGGCNTAA") is None
    assert audit.usable_atg_model("ATGTGAGCTTAA") is None
    assert audit.usable_atg_model("ATGGCT") is None


def test_decision_policy_records_uncertainty_and_exact_codon_substitution():
    rows = [
        {
            "member_key": "lib:A",
            "raw_cds_sha256": "a",
            "annotation": annotation(),
            "classification": "annotated_non_AUG_start",
            "symbol": "A",
            "raw_start_codon": "CTG",
            "transcript": "T1",
        },
        {
            "member_key": "lib:B",
            "raw_cds_sha256": "b",
            "annotation": annotation(start=False),
            "classification": "missing_start_at_transcript_boundary",
            "symbol": "B",
            "raw_start_codon": "TGG",
            "transcript": "T2",
        },
    ]
    decisions = audit.boundary_decisions(rows, manifest={"files": {}})
    assert decisions[0]["action"] == "replace_start_with_ATG"
    assert decisions[1]["action"] == "prepend_ATG_for_expression"
    assert decisions[1]["completeness"] == "missing_start_at_transcript_boundary"
    assert decisions[0]["review_status"] == "annotation_supported_pending_review"
    assert decisions[1]["review_status"] == "curation_required"
    assert all(d["approved_for_canonical_activation"] is False for d in decisions)
    assert '"native_full_length_restoration_proven":false' in decisions[1]["evidence"]


def test_union_adds_new_form_once_and_checks_selected_protein():
    form = "cyno:P:NEW"
    dna = "TTAGTTTAA"
    case = {
        "form_id": form,
        "selected_protein_sha256": audit.resolver.sha256("LV"),
        "selected": {"annotation": annotation(start=False), "transcript": "T.1"},
    }
    evidence = {"records": [case], "inputs": {"files": {}}}
    mapping = [{"form_id": form, "member_key": "ort:" + form}]
    args = {
        "manifest": {"files": {}},
        "compara_evidence": evidence,
        "selected_cds": {form: {"cds": dna}},
        "member_map": mapping,
    }
    decisions = audit.boundary_decisions([], **args)
    assert len(decisions) == 1
    assert decisions[0]["key"] == "ort:" + form
    assert decisions[0]["action"] == "prepend_ATG_for_expression"
    args["selected_cds"][form]["cds"] = "TTAGCTTAA"
    with pytest.raises(ValueError, match="evidence-selected"):
        audit.boundary_decisions([], **args)


def test_union_refuses_different_native_cds_for_an_existing_member_key():
    form = "cyno:P:NEW"
    row = {"member_key": "ort:" + form, "raw_cds_sha256": "different"}
    case = {
        "form_id": form,
        "selected_protein_sha256": audit.resolver.sha256("LV"),
        "selected": {"annotation": annotation(start=False), "transcript": "T.1"},
    }
    with pytest.raises(ValueError, match="collision"):
        audit.boundary_decisions(
            [row],
            manifest={"files": {}},
            compara_evidence={"records": [case]},
            selected_cds={form: {"cds": "TTAGTTTAA"}},
            member_map=[{"form_id": form, "member_key": "ort:" + form}],
        )


def test_alternative_difference_accounting_distinguishes_terminal_from_internal():
    terminal = alternatives.compare("MWQQA", "WQQA")
    assert terminal["relationship"] == "exact_N_terminal_extension_or_truncation"
    assert terminal["substitution_count"] == 0
    assert terminal["internal_gap_residues"] == 0
    assert terminal["gaps"][0]["location"] == "N_terminal"
    internal = alternatives.compare("MWQQAGH", "MWQEAGH")
    assert internal["substitution_count"] == 1
    assert internal["relationship"] == "internal_or_other_sequence_differences"


def test_annotated_candidate_requires_start_and_stop_not_just_absent_nf():
    a = annotation(start=False)
    a["stop_codon_features"] = [{"feature": "stop_codon"}]
    assert alternatives.annotated_start_stop(a) is False
    a["start_codon_at_CDS_start"] = True
    assert alternatives.annotated_start_stop(a) is True
    a["cds_end_NF_tag_present"] = True
    assert alternatives.annotated_start_stop(a) is False


def test_alternative_annotation_ranking_preserves_known_principal_support():
    a = annotation()
    a["stop_codon_features"] = [{"feature": "stop_codon"}]
    principal = {"transcript": "T2", "annotation": {**a, "tags": ["appris_principal_1"]}}
    canonical = {"transcript": "T1", "annotation": {**a, "tags": ["Ensembl_canonical"]}}
    assert alternatives.annotation_rank(principal) < alternatives.annotation_rank(canonical)


def test_priority_intended_protein_requires_replacement_not_extra_m():
    evidence = priority.intended_comparison("VDPG", "MDPG")
    assert evidence["first_codon_ATG_replacement_matches"] is True
    assert evidence["extra_M_prepend_matches"] is False
    assert evidence["ordinary_translation_matches"] is False
    assert evidence["aligned_comparison"]["substitution_count"] == 1


def test_priority_missing_m_is_distinct_from_first_residue_substitution():
    evidence = priority.intended_comparison("WQQA", "MWQQA")
    assert evidence["first_codon_ATG_replacement_matches"] is False
    assert evidence["extra_M_prepend_matches"] is True


def test_priority_review_refuses_existing_output_directory(tmp_path, monkeypatch):
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "review_priority_cds_forms.py",
            "--manifest",
            "not-read.json",
            "--input-dir",
            "not-read",
            "--out-dir",
            str(tmp_path),
        ],
    )
    with pytest.raises(FileExistsError, match="already exists"):
        priority.main()
