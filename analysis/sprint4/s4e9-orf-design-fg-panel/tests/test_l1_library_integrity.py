"""L1 library integrity: the tags are intact, and members do not share sequence unexplained.

Two failures these guard against are silent in the lab:

* **A modified tag.** SpyTag003 and ALFA are detected by SpyCatcher and the anti-ALFA
  nanobody; a single changed residue can lose the readout for a member, and it reads as a
  non-expressor. The L1 cassette keeps p003's tag codons (validated in cells) and recodes
  only the GGGGS linkers, so the tag DNA must also be p003's, byte for byte.
* **Shared sequence between members** (homology up to the barcode). A reverse-transcription
  template switch inside a shared stretch moves one member's barcode onto another's ORF.
  The cassette's linker encodings exist to keep C-terminal members under the 100-nt bar;
  members over it must be explained by their own ORFs (the 2026-10-13 decision accepts
  same-gene isoform/ortholog pairs) or be on the accepted list of different-gene families.

The offline tests read committed tables only (``l1_cassette_encodings.tsv``,
``l1_tag_paths.tsv``, ``l1_swap_risk.tsv``), so they run in CI. The live test recomputes
every member's shared tract from the designed CDS on S3 and skips without credentials.
"""

from __future__ import annotations

import csv
import itertools
import sys
from pathlib import Path

import pytest

from deliverome_analysis import seqtools

S4E9 = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(S4E9))

from build_constructs import CASSETTE, CASSETTE_AA, FWD_PRIMER, SPACER, STOP  # noqa: E402

#: The published tag peptides, typed here rather than read from the code under test.
SPYTAG003 = "RGVPHIVMVDAYKRYK"
ALFA = "SRLEEELRRRLTE"
#: ALFA sits between two prolines (P-ALFA-P) that isolate its helix.
P_ALFA_P = "P" + ALFA + "P"
GATE = 100
TAIL = STOP + SPACER + FWD_PRIMER
#: Two C-terminal members may share the constant tail plus at most 14 nt of the terminal
#: GGGGS block (the end rule: the whole 15-nt block is unique per member).
MAX_CASSETTE_SHARED = len(TAIL) + 14

#: Different genes whose ORFs share >100 nt of 3' coding sequence, accepted because the
#: shared stretch lies upstream of the cassette and no encoding can shorten it: clustered
#: protocadherins (shared constant exons) and close paralogs. A new gene here is a decision.
ACCEPTED_DIFFERENT_GENE_HOMOLOGY = frozenset(
    {
        "CNTNAP3",
        "CNTNAP3B",
        "NPY4R",
        "NPY4R2",
        "OPN1LW",
        "OPN1MW",
        "OR2A4",
        "OR2A7",
        "OR2T29",
        "OR2T5",
        "OR4F17",
        "OR4F5",
        "OR8U1",
        "OR8U9",
        "RAET1L",
        "ULBP2",
    }
    | {f"PCDHA{n}" for n in (1, 2, 3, 5, 6, 7, 9, 10, 11, 12, 13)}
    | {"PCDHAC1", "PCDHAC2"}
    | {f"PCDHGA{n}" for n in range(1, 13) if n != 3}
    | {f"PCDHGB{n}" for n in range(1, 8)}
    | {"PCDHGC3", "PCDHGC4", "PCDHGC5"}
)


def _tsv(name: str) -> list[dict]:
    with (S4E9 / name).open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


@pytest.fixture(scope="module")
def encodings() -> list[dict]:
    return [r for r in _tsv("l1_cassette_encodings.tsv") if r["scheme"] == "linkers_only"]


@pytest.fixture(scope="module")
def swap() -> list[dict]:
    return _tsv("l1_swap_risk.tsv")


def _tag_spans() -> dict[str, tuple[int, int]]:
    """Nucleotide span of each tag in the cassette, from its position in the peptide."""
    spans = {}
    for name, pep in (("SpyTag003", SPYTAG003), ("P-ALFA-P", P_ALFA_P)):
        start = CASSETTE_AA.index(pep) * 3
        spans[name] = (start, start + 3 * len(pep))
    return spans


# --------------------------------------------------------------------- tags are intact


def test_the_cassette_peptide_carries_unmodified_spytag003_and_alfa():
    assert CASSETTE_AA == "GGGGS" + SPYTAG003 + "GGGGSGGGGS" + P_ALFA_P + "GGGGS"
    assert seqtools.translate(CASSETTE) == CASSETTE_AA


def test_every_issued_encoding_encodes_the_unmodified_cassette(encodings):
    assert encodings, "no linkers_only encodings found"
    bad = [
        r["candidate_id"]
        for r in encodings
        if seqtools.translate(r["cassette_nt"]) != CASSETTE_AA or r["cassette_aa"] != CASSETTE_AA
    ]
    assert not bad, f"cassette peptide changed for {bad[:5]}"


def test_every_issued_encoding_keeps_p003_tag_codons(encodings):
    spans = _tag_spans()
    bad = [
        (r["candidate_id"], tag)
        for r in encodings
        for tag, (a, b) in spans.items()
        if r["cassette_nt"][a:b] != CASSETTE[a:b]
    ]
    assert not bad, f"tag DNA differs from p003's codons: {bad[:5]}"


def test_every_c_terminal_member_has_exactly_one_issued_encoding(encodings):
    cterm = {r["candidate_id"] for r in _tsv("l1_tag_paths.tsv") if r["placement"] == "C-terminal"}
    issued = [r["candidate_id"] for r in encodings]
    assert len(issued) == len(set(issued)), "a member was issued two encodings"
    assert set(issued) == cterm


def test_no_issued_encoding_carries_a_type_iis_site(encodings):
    sites = seqtools.forbidden_patterns(("BbsI", "Esp3I", "BsaI"))
    bad = [r["candidate_id"] for r in encodings if any(p in r["cassette_nt"] for p in sites)]
    assert not bad, f"BbsI/Esp3I/BsaI site inside a cassette: {bad[:5]}"


# ------------------------------------------------------- homology between members


def test_cassette_encodings_share_at_most_the_tail_plus_14_nt(encodings):
    """Pairwise shared 3' end of cassette + constant tail, the stretch before the barcode."""
    ends = sorted((r["cassette_nt"] + TAIL)[::-1] for r in encodings)
    longest = 0
    for a, b in itertools.pairwise(ends):
        n = 0
        while n < min(len(a), len(b)) and a[n] == b[n]:
            n += 1
        longest = max(longest, n)
    assert longest <= MAX_CASSETTE_SHARED, f"two encodings share {longest} nt before the barcode"


def test_terminal_linker_blocks_are_unique(encodings):
    blocks = [r["cassette_nt"][-15:] for r in encodings]
    assert len(blocks) == len(set(blocks)), "two C-terminal members share a terminal GGGGS block"


def test_no_c_terminal_member_is_over_the_bar(swap):
    over = [
        (r["candidate_id"], r["tract_nt_linkers_only"])
        for r in swap
        if r["placement"] == "C-terminal" and int(r["tract_nt_linkers_only"]) > GATE
    ]
    assert not over, f"C-terminal members over {GATE} nt with the issued encodings: {over[:5]}"


def test_members_over_the_bar_are_same_gene_or_an_accepted_family(swap):
    gene = {r["candidate_id"]: r["gene_symbol"] for r in swap}
    unexplained = sorted(
        {
            (r["gene_symbol"], gene[r["partner_linkers_only"]])
            for r in swap
            if int(r["tract_nt_linkers_only"]) > GATE
            and r["gene_symbol"] != gene[r["partner_linkers_only"]]
            and not {r["gene_symbol"], gene[r["partner_linkers_only"]]}
            <= ACCEPTED_DIFFERENT_GENE_HOMOLOGY
        }
    )
    assert not unexplained, (
        f"different genes share >{GATE} nt before the barcode and are not on the accepted "
        f"list: {unexplained[:5]}"
    )


def test_every_active_member_was_measured(swap):
    placed = {r["candidate_id"] for r in _tsv("l1_tag_paths.tsv")}
    assert {r["candidate_id"] for r in swap} == placed


# ------------------------------------------------- live: recompute from the designed CDS


def _credentials() -> bool:
    from deliverome_analysis import s3

    return s3.credentials_available()


@pytest.mark.skipif(not _credentials(), reason="needs S3 for the designed CDS")
def test_shared_tracts_recomputed_from_designed_cds_match_the_record(swap):
    import l1_swap_risk as sr

    members, unplaced = sr.load_members()
    assert not unplaced
    enc = {r["candidate_id"]: r["cassette_nt"] for r in _tsv("l1_cassette_encodings.tsv")}
    upstream = [m["body"] + (enc[m["id"]] if m["cterm"] else "") + TAIL for m in members]
    worst, _, _ = sr.worst_partner(upstream)
    recorded = {r["candidate_id"]: int(r["tract_nt_linkers_only"]) for r in swap}
    differ = [m["id"] for m, w in zip(members, worst, strict=True) if recorded[m["id"]] != w]
    assert not differ, f"shared tracts differ from l1_swap_risk.tsv for {differ[:5]}"


@pytest.mark.skipif(not _credentials(), reason="needs S3 for the p003 map")
def test_the_cassette_constant_is_p003s_own_dna(tmp_path):
    """The codons every L1 tag keeps are the ones in the built, validated p003 plasmid."""
    from Bio import SeqIO

    from deliverome_analysis.s3 import processed

    path = processed.download("constructs/ordered/p003.gb", tmp_path / "p003.gb")
    p003 = str(SeqIO.read(path, "genbank").seq).upper()
    assert CASSETTE in p003, "the cassette's DNA is not in the p003 map"
