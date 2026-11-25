"""The synthesis-compliance checker reproduces Twist's hard calls, and the recoder only ever changes synonymous, unprotected codons."""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import pytest

S4E9 = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(S4E9))

import synthesis_compliance as sc  # noqa: E402

from deliverome_analysis import seqtools  # noqa: E402

# Twist's own refusals on the draft (issue_regions), quoted from its report.
TWIST_ME_HITS = [
    "CTGCTCTTCTCACATCT",
    "CTGTGTTTTATACCATCT",
    "AGTGTGTATAAGGAAG",
    "CAGTCATCTTATCACATCT",
]


@pytest.mark.parametrize("hit", TWIST_ME_HITS)
def test_twist_mosaic_end_hits_are_refused(hit):
    assert sc.me_distance("ACGT" * 10 + hit + "ACGT" * 10) <= sc.ME_MAX_REFUSED


def test_unrelated_sequence_is_not_an_me_hit():
    assert sc.me_distance("GATTACAGGCCTAGGCATGCAAGCTTGGCGTAATCATGG") > sc.ME_MAX_REFUSED


def test_gc_profile_counts_barcode_n_as_half():
    assert sc.gc_profile("N" * 50, 50) == [0.5]
    assert sc.gc_profile("GC" * 25, 50) == [1.0]


def test_near_repeat_finds_a_tandem_repeat_with_mismatches():
    unit = "ACGTTGCAAGCTTCGATCGGATCCAAGTTCGA"  # 32 nt
    copy2 = unit[:10] + "T" + unit[11:]  # one mismatch
    s = "GG" + unit + copy2 + unit + "CC"
    n, _i, _j = sc.near_repeat(s)
    assert n >= 64


def _toy_insert(protein_orf: str) -> dict:
    return {
        "rung": "cterm",
        "construct": "toy",
        "design_id": "l1-p0000",
        "gene_symbol": "TOY",
        "five": sc.FIVE,
        "orf": protein_orf,
        "three": "CCAGTCGTTAGCTGCCGTTCC" + "N" * 20 + "GGATCTGGAC",
        "sp_aa": 0,
    }


def _cassette_dna() -> str:
    rows = list(csv.DictReader((S4E9 / "constructs_l1_cterm.tsv").open(), delimiter="\t"))
    fasta = (S4E9 / "constructs_l1_cterm.fasta").read_text().split(">")[1:]
    seq = fasta[0].split("\n", 1)[1].replace("\n", "")
    _f, orf, _t = sc.split_insert(seq)
    assert rows
    return sc.cassette_dna(orf)


def test_recoder_fixes_a_gc_rich_stretch_without_changing_protein_or_cassette():
    gc_rich = "GCCGCCGGCGCCGCGGGCCCCGCCGGCGCCGCGGCCGCCGGCGCCGCGGGCCCCGCC"  # Ala/Gly/Pro, 57 nt
    orf = "ATG" + gc_rich * 2 + "AAAGAAAAA" + _cassette_dna() + "TAA"
    ins = _toy_insert(orf)
    before = sc.metrics(ins["five"] + orf + ins["three"])
    assert before["max_gc50"] > sc.GC50_MAX
    rec = sc.Recoder(ins, seqtools.synonymous_codons(min_usage=sc.MIN_USAGE))
    out = rec.run()
    new_orf = "".join(rec.codons)
    assert seqtools.translate(new_orf) == seqtools.translate(orf)
    assert sc.cassette_dna(new_orf) == sc.cassette_dna(orf)
    assert new_orf[:3] == "ATG"
    assert sc.metrics(out)["max_gc50"] < before["max_gc50"]
    assert rec.edits and all(e["from_codon"] != e["to_codon"] for e in rec.edits)


def test_recoder_rejects_new_polya_signal_and_committed_recodes_use_common_codons():
    # Toy ORF where a candidate edit would create ATTAAA or AATAAA must be rejected by acceptable()
    orf = "ATG" + "GCT" * 5 + "GAA" + "AAG" + _cassette_dna() + "TAA"
    ins = _toy_insert(orf)
    rec = sc.Recoder(ins, seqtools.synonymous_codons(min_usage=sc.MIN_USAGE))
    # Place AAT directly upstream of AAA (codon 6 -> AAT creates AATAAA if codon 7 is AAA)
    rec.codons[6] = "AAC"
    rec.codons[7] = "AAA"
    rec.polya0 = sum(rec.seq().count(p) for p in sc.POLYA)
    assert rec.acceptable(6, "AAT") is False

    # Every committed recode must use a >=10% human codon and include the 14 polya/donor fixes
    syn_map = {
        c: freq
        for aa_list in seqtools.synonymous_codons(min_usage=0.0).values()
        for c, freq in aa_list
    }
    with (S4E9 / "synthesis_recodes.tsv").open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert rows
    assert all(syn_map[r["to_codon"]] >= sc.MIN_USAGE for r in rows)
    reasons = [r["reason"] for r in rows]
    assert reasons.count("domestication_polya") == 4
    assert reasons.count("domestication_splice_donor") == 10
