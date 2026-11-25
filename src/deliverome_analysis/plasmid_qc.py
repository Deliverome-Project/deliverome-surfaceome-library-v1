"""Pre-order checks for a Deliverome expression plasmid.

An ordered plasmid is expensive and slow to correct, and the ways it goes wrong are
quiet: a map rebuilt from a superseded design file, a cassette that reverted to the shared
constant, a barcode that drifted from the registry, an ORF that is a different proteoform
from the one the library selected. None of those look wrong in a viewer.

Each check returns a :class:`Finding` rather than raising, so one run reports everything
that is wrong instead of the first thing. ``check_plasmid`` runs them all; the CLI exits
non-zero if any check fails, which is what makes it usable as an order gate.

    uv run --locked python -m deliverome_analysis.plasmid_qc map.gb --spec spec.json

The spec says what this plasmid is supposed to be -- which gene, which designed insert,
which barcode -- because every interesting check is a comparison against an intention, and
a checker that only looks at the file can only find syntax errors.
"""

from __future__ import annotations

import json
from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from . import seqtools

#: The three enzymes every Deliverome ORF is domesticated against. BbsI assembles,
#: Esp3I shuttles, BsaI was never cleared from the backbone.
ASSEMBLY_ENZYMES = ("BbsI", "BsaI")
SHUTTLE_ENZYME = "Esp3I"

KOZAK = "GACCGGTTCTAGAGCCGCCACC"
STOP = "TAA"
SPACER = "CC"
FWD_PRIMER = "AGTCGTTAGCTGCCGTTCC"
REV_PRIMER = "GGATCTGGACCCTTTGGGCAAGCACACGTC"
CS1 = "GCTTTAAGGCCGGTCCTAGCAA"
BARCODE_LEN = 20
#: Bases allowed between the backbone arc and its tail without being the stuffer: the
#: Esp3I site and overhang flanking the insert. Measured at 15 bp on this panel.
MAX_SHUTTLE_SCAR = 40

#: Labels every ordered map must carry. A map whose features are missing or mislabelled
#: misleads whoever opens it, and annotation is the part nobody re-derives.
REQUIRED_FEATURES = (
    "Kozak (AgeI-XbaI)",
    "stop",
    "pad",
    "5' handle",
    "barcode",
    "gap + RT site",
    "10x Capture Sequence 1",
    "EFS_v1",
    "IRES2_v1",
    "PuroR",
    "WPRE",
)

#: The reference cassette encoding carried by p003. A construct that still has it has not
#: been through diversify_cassettes.py, which is the state that puts the same 177 nt in
#: every member and makes the barcode a template-switch substrate for every pair at once.
P003_CASSETTE_PEPTIDE = "GGGGSRGVPHIVMVDAYKRYKGGGGSGGGGSPSRLEEELRRRLTEPGGGGS"
LOOP_CASSETTE_PEPTIDE = "GSPSRLEEELRRRLTEPGS"

#: p003's own DNA for that cassette: the encoding validated in cells, and the one every L1
#: arm started from (``build_constructs.CASSETTE`` in s4e9; its test pins the two equal and
#: checks this DNA is in the p003 map). Two uses: its SpyTag003 and P-ALFA-P codons are the
#: ones an L1 arm must keep, and its GGGGS linker codons are the ones it must NOT keep, since
#: one encoding shared across a pool is an identical 153-nt stretch in every member.
P003_CASSETTE_NT = (
    "GGAGGAGGAGGCAGC"
    "AGAGGCGTGCCACACATCGTGATGGTGGATGCCTACAAGAGATACAAG"
    "GGAGGAGGCGGCTCTGGAGGAGGAGGCAGC"
    "CCC"
    "AGTCGCCTCGAAGAAGAGCTCAGGCGGAGACTGACAGAG"
    "CCT"
    "GGAGGAGGAGGCAGC"
)
#: The two tags inside the cassette. Their DNA is what the readout binds and what was shown to
#: display, so it is held on p003's codons; everything else in the cassette is linker.
TAG_PEPTIDES: dict[str, str] = {"SpyTag003": "RGVPHIVMVDAYKRYK", "P-ALFA-P": "PSRLEEELRRRLTEP"}

#: Murine Ig kappa-chain V-J2-C signal peptide (20 aa / 60 nt ending GGT) plus the +1 Asp
#: every published Igk vector keeps (Thermo pSecTag2/Hygro; Addgene 219810, 259743, 259744).
IGK_LEADER_NT = "ATGGAGACAGACACACTCCTGCTATGGGTACTGCTGCTCTGGGTTCCAGGTTCCACTGGT"
IGK_LEADER_AA = "METDTLLLWVLLLWVPGSTG"
IGK_PLUS1_NT = "GAC"
IGK_PLUS1_AA = "D"

# --- The insert wrapper ------------------------------------------------------------------- #
# Every insert is ``ESP3I_5 | KOZAK | coding (with its stop) | readout | ESP3I_3``. An L1 rung
# insert carries the L1 readout region (:func:`l1_readout`, :func:`wrap_insert`); every plasmid
# ordered before L1, and the L1 rescue arms, carry the p003 readout ``SPACER | FWD_PRIMER |
# barcode | REV_PRIMER | CS1`` (:func:`wrap_p003_insert`, ``P003_STOP_TO_BARCODE``). The rung
# builders, the pool QC and the synthesis checks all assemble or parse these layouts, so they
# are defined once, here.
#
# L1 readout region, adopted 2026-11-22 (Becca: "let's just leave the backbone as is but then
# add the t7 and the feldman flank"). The p003 constants above stay as they are: they describe
# the layout that every plasmid ordered before L1 carries, and the map audits check against
# them. L1 inserts are assembled with ``l1_readout`` instead:
#
#   STOP · SPACER · T7_PROMOTER · FWD_PRIMER (JT832) · barcode · FELDMAN_3PRIME · CS1
#
# * T7_PROMOTER: T7 class III promoter, -17 to +3, forward, so in fixed cells T7 transcribes
#   sense RNA across the handle, barcode and readout region that Feldman's RT primer and padlock
#   read unchanged (in situ T7 transcription before RT, as in PerturbView, Kudo et al. 2025,
#   doi:10.1038/s41587-024-02391-0). It lengthens the sequence every insert shares before its
#   barcode by 20 nt, which the pooled swap model prices (s4e9 swap_t7_scenario.py).
# * FELDMAN_3PRIME: replaces the 10-nt gap + RT site with Feldman's 17-nt padlock construct
#   from lentiGuide-BC (Addgene #127168) followed by his RT site with T->A at position 9, which
#   removes a GTAAG splice donor and a hairpin (reviewer brief, section 9). Same as
#   s4e9 reserve_barcodes.FELDMAN_T9A_3PRIME, which every L1 barcode was screened against.
# The insert still ends at CS1, so the Esp3I junction and the destination backbone are unchanged.
T7_PROMOTER = "TAATACGACTCACTATAGGG"
FELDMAN_3PRIME = "ACTGGCTATTCATTCGC" + "CCTTTGGGAAAGCACACGTC"

#: The Esp3I sites that keep every insert excisable for shuttling, read off p0039_AXL:
#: ``CGTCTC A|GACC...`` at the 5' end (the overhang is the Kozak's own first four bases) and
#: ``...AAAG|A GAGACG`` at the 3' end (bottom strand; the last bases of 10x CS1 + the AG pad).
#: p0039 carries no Esp3I of its own.
ESP3I_5 = "CGTCTCA"
ESP3I_3 = "AGA" + "GAGACG"
#: The barcode position in a designed insert, before a barcode is assigned.
BARCODE_SLOT = "N" * BARCODE_LEN
#: 0-based offset of the coding sequence's first base (the ATG) in a wrapped insert.
CODING_START = len(ESP3I_5) + len(KOZAK)
#: The two constant halves of ``l1_readout`` either side of the barcode, for code that reads an
#: L1 insert's 3' end (seam contexts, swap-risk tails, annotation, QC).
L1_BEFORE_BARCODE = STOP + SPACER + T7_PROMOTER + FWD_PRIMER
L1_AFTER_BARCODE = FELDMAN_3PRIME + CS1
#: The 44 nt before the barcode of every L1 insert that ends on the cassette (a C-terminal
#: arm): the stop codon, the pad, the T7 promoter and the 5' handle. Arms that keep a native
#: C-terminus end on their own stop codon instead, followed by the same pad, T7 and handle.
STOP_TO_BARCODE = L1_BEFORE_BARCODE
#: The 24 nt before the barcode in the p003 layout (stop, pad, 5' handle), for pre-L1 maps.
P003_STOP_TO_BARCODE = STOP + SPACER + FWD_PRIMER
#: Where the wrapped insert replaces p0039's MCS + stuffer: 0-based start, immediately 3' of
#: EFS, and the length replaced (MCS 28 + stuffer 167).
BACKBONE_SLOT_START = 4808
BACKBONE_SLOT_LEN = 195
#: The enzymes every ORF is domesticated against, in the order the builders search them.
DOMESTICATION_ENZYMES: tuple[str, ...] = ("BbsI", SHUTTLE_ENZYME, "BsaI")


def l1_readout(barcode: str) -> str:
    """Everything an L1 insert carries after its stop codon, for one barcode."""
    assert len(barcode) == BARCODE_LEN, "barcodes are 20 nt"
    return STOP + SPACER + T7_PROMOTER + FWD_PRIMER + barcode + FELDMAN_3PRIME + CS1


assert l1_readout(BARCODE_SLOT) == L1_BEFORE_BARCODE + BARCODE_SLOT + L1_AFTER_BARCODE


def wrap_insert(coding: str, barcode: str = BARCODE_SLOT) -> str:
    """The ordered L1 insert around ``coding``, which must already end in its stop codon."""
    return ESP3I_5 + KOZAK + coding + l1_readout(barcode)[len(STOP) :] + ESP3I_3


def wrap_p003_insert(coding: str, barcode: str = BARCODE_SLOT) -> str:
    """The insert around ``coding`` in the p003 readout layout (``coding`` ends in its stop).

    The layout every plasmid ordered before L1 carries; the L1 rescue arms keep it.
    """
    return ESP3I_5 + KOZAK + coding + SPACER + FWD_PRIMER + barcode + REV_PRIMER + CS1 + ESP3I_3


def wrapper_sites_only(insert: str, enzymes: tuple[str, ...] = DOMESTICATION_ENZYMES) -> bool:
    """True when the only Type IIS sites in ``insert`` are the wrapper's two Esp3I."""
    sites = seqtools.find_sites(insert, list(enzymes))
    return len(sites) == 2 and all(s[0] == SHUTTLE_ENZYME for s in sites)


def place_in_backbone(backbone: str, insert: str) -> str:
    """The finished plasmid: ``insert`` in place of the empty p0039's MCS + stuffer slot."""
    return (
        backbone[:BACKBONE_SLOT_START]
        + insert
        + backbone[BACKBONE_SLOT_START + BACKBONE_SLOT_LEN :]
    )


def check_inserts_in_backbone(
    inserts: list[tuple[str, str]],
    backbone: str,
    enzymes: tuple[str, ...] = DOMESTICATION_ENZYMES,
) -> int:
    """Verify each ``(label, insert)`` in ``backbone`` adds only the two wrapper Esp3I sites.

    Compares per-enzyme site counts rather than total count alone, so an extra BsaI/BbsI site
    cannot mask a missing Esp3I site. Returns the backbone's own site count.
    """
    enz_list = list(enzymes)
    bb_sites = seqtools.find_sites(backbone, enz_list)
    expected: dict[str, int] = {e: sum(s[0] == e for s in bb_sites) for e in enzymes}
    expected[SHUTTLE_ENZYME] = expected.get(SHUTTLE_ENZYME, 0) + 2
    for label, insert in inserts:
        vec = place_in_backbone(backbone, insert)
        sites = seqtools.find_sites(vec, enz_list)
        got = {e: sum(s[0] == e for s in sites) for e in enzymes}
        if got != expected:
            raise AssertionError(f"{label}: vector sites {got}, expected {expected}")
    return len(bb_sites)


def translate_insert(insert: str) -> str:
    """Protein of a linear wrapped insert, from the Kozak ATG to the first in-frame stop.

    The stop is not included. A codon outside the genetic code (an ``N``) reads as ``X``, so
    a barcode slot reached before any stop shows up in the protein rather than raising.
    Raises ``ValueError`` when the insert has no Kozak or no ATG right after it. For a
    circular plasmid use :func:`read_orf`.
    """
    k = insert.find(KOZAK)
    if k < 0:
        raise ValueError("insert carries no Kozak")
    start = k + len(KOZAK)
    if insert[start : start + 3] != "ATG":
        raise ValueError("insert does not start ATG after the Kozak")
    aa = []
    for i in range(start, len(insert) - 2, 3):
        c = seqtools.GENETIC_CODE.get(insert[i : i + 3], "X")
        if c == "*":
            break
        aa.append(c)
    return "".join(aa)


def p003_tag_codons() -> dict[int, str]:
    """``{cassette residue index: codon}`` for SpyTag003 and P-ALFA-P on p003's own DNA.

    The ``fixed_codons`` every per-arm cassette encoding is drawn with (Becca, 2026-10-29):
    the tags keep p003's validated codons and only the GGGGS linker blocks are recoded.
    """
    return seqtools.pin_codons(
        P003_CASSETTE_NT, P003_CASSETTE_PEPTIDE, tuple(TAG_PEPTIDES.values())
    )


def keeps_p003_tag_codons(cassette_nt: str) -> bool:
    """True when a 153-nt cassette encoding carries SpyTag003 and P-ALFA-P on p003's codons."""
    return all(cassette_nt[3 * i : 3 * i + 3] == c for i, c in p003_tag_codons().items())


#: Type IIS enzymes the assembly does not use but a downstream cloning step might. Reported as
#: warnings, never gated (recognition sites, top strand).
EXTRA_TYPE_IIS: dict[str, str] = {
    "AarI/PaqCI": "CACCTGC",  # isoschizomers: one site, reported once
    "BtgZI": "GCGATG",
    "SapI": "GCTCTTC",
}
#: Synthesis sanity limits, reported as warnings by :func:`check_synthesis_windows`. The GC
#: ceiling is the one the L1 order was made compliant with (s4e9 ``synthesis_compliance.GC50_MAX``,
#: Twist's 50-nt rule); the floor and the homopolymer limit are Twist's published ones.
SYNTHESIS_GC_WINDOW = 50
SYNTHESIS_GC_MIN = 0.25
SYNTHESIS_GC_MAX = 0.80
SYNTHESIS_MAX_HOMOPOLYMER = 7


@dataclass
class Finding:
    """One check's result. ``ok`` false is a reason not to order."""

    check: str
    ok: bool
    detail: str
    severity: str = "error"  # error | warning


@dataclass
class Spec:
    """What this plasmid is meant to be.

    Every field is optional so a partially specified plasmid still gets the checks its
    spec supports, and the ones it cannot support are reported as skipped rather than
    silently passing.
    """

    plasmid: str
    gene: str | None = None
    designed_insert: str | None = None
    backbone: str | None = None
    barcode: str | None = None
    cassette_nt: str | None = None
    cassette_peptide: str | None = None
    proteoforms: dict[str, str] = field(default_factory=dict)
    #: Per-residue DeepTMHMM **v1** topology of the encoded protein, in the ordered
    #: sequence's own numbering. Passed in rather than fetched: the checks stay pure and
    #: testable, and the caller already has it from
    #: ``deliverome_analysis.topology.for_sequence``. Empty means "not checked".
    topology: str | None = None
    #: True when the gene is GPI-anchored, so a protein with no TM helix is still held in
    #: the membrane. None means "unknown", which is reported rather than assumed either way.
    gpi_anchored: bool | None = None
    #: Native residue the cassette follows (``native[:n] + cassette + native[n:]``), i.e. the
    #: cassette's 0-based start in the encoded protein. None means "not checked".
    cassette_after_residue: int | None = None
    #: The whole tagged protein the ORF must encode, when known. None means "not checked".
    tagged_protein: str | None = None
    #: True only for p003 itself and its own arms: they are allowed to carry p003's shared
    #: cassette DNA. Every other plasmid fails ``linkers_not_shared`` if it does.
    reference_cassette_ok: bool = False
    #: Require SpyTag003 and P-ALFA-P to be on p003's codons (every L1 arm; Becca,
    #: 2026-10-29). Off by default because the 2026-10-07 panel recoded the tags too.
    require_p003_tag_codons: bool = False
    #: ``{label: sequence}`` of backbone features that must each occur exactly once in the
    #: finished plasmid (from the backbone's GenBank). Empty means "not checked".
    backbone_features: dict[str, str] = field(default_factory=dict)
    #: What follows the ORF's stop codon up to the barcode, and what follows the barcode. The
    #: defaults are the legacy readout (pad + 5' handle; gap + RT site). An L1 insert carries the
    #: L1 readout region instead (Becca, 2026-11-22: pad + T7 promoter + 5' handle; Feldman's
    #: padlock + RT site), so its spec names those (``L1_BEFORE_BARCODE``, ``FELDMAN_3PRIME``).
    after_stop: str = SPACER + FWD_PRIMER
    after_barcode: str = REV_PRIMER

    @classmethod
    def from_json(cls, path: Path) -> Spec:
        data = json.loads(Path(path).read_text())
        return cls(**data)


def _circular_find(haystack: str, needle: str) -> int:
    """Index of ``needle`` in a circular sequence, or -1.

    A plasmid has no canonical start, so a map saved from a rotated view holds the same
    molecule with every linear offset changed. Searching the doubled sequence makes the
    checks rotation-invariant instead of failing on a cosmetic difference.
    """
    if not needle:
        return -1
    doubled = haystack + haystack[: len(needle)]
    return doubled.find(needle)


def _common_prefix(a: str, b: str) -> int:
    """Length of the common prefix of two strings (binary search on C-level slice compares)."""
    lo, hi = 0, min(len(a), len(b))
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if a[:mid] == b[:mid]:
            lo = mid
        else:
            hi = mid - 1
    return lo


def _is_rotation(a: str, b: str) -> bool:
    """True when two circular sequences are the same molecule at a different start."""
    return len(a) == len(b) and b in (a + a)


def _circular_count(haystack: str, needle: str) -> int:
    """Occurrences of ``needle`` in a circular sequence (each start position counted once)."""
    if not needle or len(needle) > len(haystack):
        return 0
    doubled = haystack + haystack[: len(needle) - 1]
    n, at = 0, doubled.find(needle)
    while 0 <= at < len(haystack):
        n += 1
        at = doubled.find(needle, at + 1)
    return n


def read_orf(seq: str) -> tuple[int, str] | None:
    """``(start, coding)`` of the ORF from the Kozak ATG to its first in-frame stop.

    ``coding`` includes the stop codon. None when there is no Kozak, no ATG after it, or no
    in-frame stop within one turn of the (circular) sequence.
    """
    k = _circular_find(seq, KOZAK)
    if k < 0:
        return None
    start = (k + len(KOZAK)) % len(seq)
    doubled = seq + seq
    if doubled[start : start + 3] != "ATG":
        return None
    for i in range(start, start + len(seq) - 2, 3):
        if seqtools.translate(doubled[i : i + 3]) == "*":
            return start, doubled[start : i + 3]
    return None


def _tandem_in_orf(seq: str) -> tuple[str, str, int] | None:
    """``(coding, protein, cassette_aa_index)`` when the ORF carries the 51-aa tandem once."""
    orf = read_orf(seq)
    if orf is None:
        return None
    coding = orf[1]
    protein = seqtools.translate(coding[:-3])
    if protein.count(P003_CASSETTE_PEPTIDE) != 1:
        return None
    return coding, protein, protein.find(P003_CASSETTE_PEPTIDE)


def check_backbone(seq: str, spec: Spec) -> Finding:
    """The vector is the destination backbone with the insert in its cloning site.

    The backbone reference is the EMPTY destination vector, which still carries its
    cloning stuffer (the AgeI/XbaI/AsiSI segment the insert replaces). So the finished
    vector is not "backbone plus insert": it is the backbone with one contiguous segment
    swapped out. Testing for the backbone as a substring, or for backbone length plus
    insert length, fails every correct plasmid -- both were tried before this.

    The test used here is the right one: remove the insert, and what is left must be the
    backbone with a single contiguous piece missing. The size of that piece is reported,
    because a sudden change in it means the cloning site moved.
    """
    if not spec.backbone:
        return Finding("backbone", True, "no backbone in spec — not checked", "warning")
    if not spec.designed_insert:
        return Finding(
            "backbone",
            True,
            "no designed insert, so the backbone cannot be isolated — not checked",
            "warning",
        )
    at = _circular_find(seq, spec.designed_insert)
    if at < 0:
        return Finding(
            "backbone",
            True,
            "insert not found, so the backbone cannot be isolated — see insert_matches_design",
            "warning",
        )

    doubled = seq + seq
    remainder = doubled[at + len(spec.designed_insert) : at + len(seq)]
    # Fast path: the remainder starts right after the insert, so in a correct plasmid its
    # opening bases sit in the backbone just after the cloning site. Anchoring there gives the
    # same arc as the general search below in a single comparison, which is what makes this
    # check affordable on a 3,000-insert pool. Any other layout falls back to the search.
    anchor = (spec.backbone + spec.backbone).find(remainder[:50]) if len(remainder) >= 50 else -1
    if anchor >= 0:
        i, j = 0, anchor % len(spec.backbone)
        n = _common_prefix(remainder, spec.backbone[j:] + spec.backbone[:j])
    if anchor < 0 or n < 100:
        n, i, j = seqtools.longest_shared_tract(remainder, spec.backbone)
    if n < 100:
        return Finding(
            "backbone", False, f"shares only {n} bp with the expected backbone — wrong backbone"
        )

    # Rotate both so the shared arc starts at zero, then the rest of the remainder must be
    # the backbone's tail; whatever the backbone has in between is the stuffer.
    rot_r = remainder[i:] + remainder[:i]
    rot_b = spec.backbone[j:] + spec.backbone[:j]
    # The remainder is arc + <shuttle scar> + tail, while the backbone is
    # arc + <stuffer> + tail. The scar is the Esp3I site and overhang flanking the insert;
    # they belong to neither the insert (which starts at the Kozak) nor the stuffer, so a
    # few bases between arc and tail are expected architecture, not a difference.
    tail = rot_r[n:]
    common = 0
    while common < min(len(tail), len(rot_b)) and tail[-1 - common] == rot_b[-1 - common]:
        common += 1
    scar = len(tail) - common
    stuffer = len(rot_b) - n - common
    covered = (n + common) / len(rot_r)
    if scar <= MAX_SHUTTLE_SCAR and covered >= 0.99:
        return Finding(
            "backbone",
            True,
            f"p0039 backbone, insert in the cloning site "
            f"(stuffer replaced: {stuffer} bp; shuttle scar: {scar} bp)",
        )
    return Finding(
        "backbone",
        False,
        f"backbone differs beyond the cloning site: {n} bp arc and {common} bp tail match "
        f"but {scar} bp between them is unaccounted for "
        f"({covered:.1%} of the non-insert sequence explained)",
    )


def check_shuttle_sites(seq: str) -> Finding:
    """Exactly two Esp3I sites, which are the shuttle handles the insert sits between.

    Intact is the requirement here, not absent: Esp3I is how the insert moves to another
    backbone. Any number other than two means the insert cannot be shuttled cleanly.
    """
    sites = _sites_circular(seq, [SHUTTLE_ENZYME])
    n = len({s[2] for s in sites})
    if n == 2:
        return Finding("esp3i_shuttle", True, "2 Esp3I sites, as expected")
    return Finding(
        "esp3i_shuttle",
        False,
        f"{n} Esp3I site(s); the shuttle needs exactly 2 flanking the insert",
    )


def check_domestication(seq: str, spec: Spec) -> Finding:
    """No BbsI or BsaI anywhere: either would cut during Golden Gate assembly."""
    region = spec.designed_insert or seq
    sites = seqtools.find_sites(region, list(ASSEMBLY_ENZYMES))
    if not sites:
        where = "insert" if spec.designed_insert else "whole plasmid"
        return Finding("golden_gate_domestication", True, f"no BbsI or BsaI in the {where}")
    return Finding(
        "golden_gate_domestication",
        False,
        f"{len(sites)} assembly-enzyme site(s) remain: {sites[:4]}",
    )


def check_insert_matches_design(seq: str, spec: Spec) -> Finding:
    """The map carries the CURRENT designed insert, byte for byte.

    This is the check that catches a stale map. A vector rebuilt from a superseded design
    file looks entirely well formed -- right backbone, right features, right length -- and
    differs only in the sequence nobody re-reads.
    """
    if not spec.designed_insert:
        return Finding("insert_matches_design", True, "no designed insert in spec", "warning")
    if _circular_find(seq, spec.designed_insert) >= 0:
        return Finding(
            "insert_matches_design",
            True,
            f"insert matches the design ({len(spec.designed_insert)} nt)",
        )
    shared = seqtools.longest_shared_tract(seq, spec.designed_insert)[0]
    return Finding(
        "insert_matches_design",
        False,
        f"insert differs from the current design; longest shared tract {shared} nt of "
        f"{len(spec.designed_insert)} — the map may have been built from a superseded file",
    )


def check_diversified_linker(seq: str, spec: Spec) -> Finding:
    """The cassette carries this construct's OWN encoding, not the shared reference.

    Hardened 2026-11-21: an issued encoding that IS p003's shared DNA used to pass (the check
    only asked that the spec's encoding be present), and so did an encoding found out of frame
    or twice. It now fails all three.
    """
    if not spec.cassette_nt:
        return Finding("diversified_linker", True, "no cassette in spec", "warning")
    if spec.cassette_nt == P003_CASSETTE_NT and not spec.reference_cassette_ok:
        return Finding(
            "diversified_linker",
            False,
            "the spec's encoding is p003's shared cassette DNA, not an issued per-construct one",
        )
    copies = _circular_count(seq, spec.cassette_nt)
    if copies == 0:
        return Finding(
            "diversified_linker", False, "this construct's issued cassette encoding is not present"
        )
    if copies > 1:
        return Finding(
            "diversified_linker", False, f"the issued encoding occurs {copies} times in the plasmid"
        )
    orf = read_orf(seq)
    if orf is not None:
        at = orf[1].find(spec.cassette_nt)
        if at < 0 or at % 3:
            return Finding(
                "diversified_linker",
                False,
                "the issued encoding is not in the ORF's reading frame",
            )
    peptide = spec.cassette_peptide or seqtools.translate(spec.cassette_nt)
    if seqtools.translate(spec.cassette_nt) != peptide:
        return Finding(
            "diversified_linker",
            False,
            "the issued encoding does not translate to the expected tag peptide",
        )
    return Finding(
        "diversified_linker",
        True,
        f"carries its own {len(spec.cassette_nt)} nt encoding of {peptide[:18]}…",
    )


def check_barcode(seq: str, spec: Spec) -> Finding:
    """The reserved barcode is present, the right length, and between its handles."""
    if not spec.barcode:
        return Finding("barcode", True, "no barcode in spec", "warning")
    if len(spec.barcode) != BARCODE_LEN:
        return Finding(
            "barcode", False, f"barcode is {len(spec.barcode)} nt, expected {BARCODE_LEN}"
        )
    flanked = FWD_PRIMER + spec.barcode + spec.after_barcode
    copies = _circular_count(seq, flanked)
    if copies > 1:
        return Finding("barcode", False, f"barcode slot occurs {copies} times; expected once")
    if copies == 1:
        if _circular_count(seq, FWD_PRIMER) != 1 or _circular_count(seq, spec.after_barcode) != 1:
            return Finding(
                "barcode", False, "a second copy of the 5' handle or RT site is in the plasmid"
            )
        return Finding("barcode", True, "reserved barcode present once, correctly flanked")
    if _circular_find(seq, spec.barcode) >= 0:
        return Finding(
            "barcode", False, "barcode present but not between the 5' handle and the RT site"
        )
    return Finding("barcode", False, "reserved barcode not found")


def check_proteoform(seq: str, spec: Spec) -> Finding:
    """With the tag removed, the ORF is one of the proteoforms we selected.

    Works at the PROTEIN level: the tag is located by its peptide, not by its nucleotide
    encoding, so a map whose cassette carries a different or undiversified encoding is
    still checked for the right proteoform rather than reporting that it cannot be read.
    Ordering an isoform the library did not select is invisible in a map and expensive
    downstream, because the tag site and every structural signal were computed against
    the selected form.
    """
    if not spec.proteoforms:
        return Finding("proteoform", True, "no proteoform in spec", "warning")
    start = _circular_find(seq, KOZAK)
    if start < 0:
        return Finding("proteoform", False, "no Kozak found, so cannot find the start codon")
    doubled = seq + seq
    i = start + len(KOZAK)
    protein = []
    while i < len(doubled) - 2:
        aa = seqtools.translate(doubled[i : i + 3])
        if aa == "*":
            break
        protein.append(aa)
        i += 3
    got = "".join(protein)

    tag = spec.cassette_peptide or (
        seqtools.translate(spec.cassette_nt) if spec.cassette_nt else None
    )
    if tag and tag in got:
        got = got.replace(tag, "", 1)
    elif tag:
        return Finding(
            "proteoform",
            False,
            f"the {len(tag)} aa tag peptide is not in the ORF, so the native sequence "
            "cannot be recovered",
        )

    for name, ref in spec.proteoforms.items():
        if got == ref:
            return Finding("proteoform", True, f"ORF matches the selected proteoform {name}")
    best = max(
        spec.proteoforms,
        key=lambda n: seqtools.longest_shared_tract(got, spec.proteoforms[n])[0],
    )
    return Finding(
        "proteoform",
        False,
        f"untagged ORF is {len(got)} aa and matches no selected proteoform; "
        f"closest is {best} at {len(spec.proteoforms[best])} aa",
    )


def check_membrane_anchor(seq: str, spec: Spec) -> Finding:
    """The encoded protein must have something holding it in the membrane.

    A construct whose protein has neither a transmembrane helix nor a GPI anchor is
    secreted, so no tag site on it can ever be displayed on a live cell. Every other check
    in this module tests DNA-level correctness -- sites, backbone, linkers, barcodes -- or
    proteoform identity. None of them asks whether the protein is a surface protein at all,
    and that is how MICA reached an order list: its L1 reference resolved to
    ``ENST00000449934``, a 332-aa form that has lost the transmembrane domain the 383-aa
    canonical carries, and DeepTMHMM labels it ``SP`` -- secreted. The proteoform checks
    passed, because the insert faithfully encoded the proteoform that was selected. The
    selection ladder answers "which proteoform is the main one", which is a different
    question from "can this proteoform be displayed".

    GPI counts as an anchor: 173 of the 174 zero-TM L1 members are GPI-anchored and are
    perfectly good surface proteins. Only the remaining case is a design error.
    """
    if not spec.topology:
        return Finding("membrane_anchor", True, "no topology in spec — not checked", "warning")
    tm_runs = sum(
        1
        for i, c in enumerate(spec.topology)
        if c == "M" and (i == 0 or spec.topology[i - 1] != "M")
    )
    if tm_runs:
        return Finding("membrane_anchor", True, f"{tm_runs} TM helix/helices")
    if spec.gpi_anchored:
        return Finding("membrane_anchor", True, "no TM helix, but GPI-anchored")
    if spec.gpi_anchored is None:
        return Finding(
            "membrane_anchor",
            False,
            "no TM helix and GPI status unknown — cannot confirm a membrane anchor",
            "warning",
        )
    return Finding(
        "membrane_anchor",
        False,
        "no TM helix and not GPI-anchored: this protein is secreted, so no tag site on "
        "it can be displayed on a live cell",
    )


def check_single_orf(seq: str, spec: Spec) -> Finding:
    """One continuous reading frame from the Kozak ATG to a single stop."""
    start = _circular_find(seq, KOZAK)
    if start < 0:
        return Finding("single_orf", False, "no Kozak found")
    doubled = seq + seq
    i = start + len(KOZAK)
    if doubled[i : i + 3] != "ATG":
        return Finding("single_orf", False, f"no ATG after the Kozak, found {doubled[i : i + 3]}")
    stops = 0
    while i < len(doubled) - 2:
        if seqtools.translate(doubled[i : i + 3]) == "*":
            stops += 1
            break
        i += 3
    if stops != 1:
        return Finding("single_orf", False, "no in-frame stop codon found")
    codons = (i - start - len(KOZAK)) // 3
    # Hardened 2026-11-21: a premature stop or a frameshift used to pass here, because any
    # in-frame stop counted. Every Deliverome insert ends its ORF on a stop followed by the
    # pad and the 5' handle, so when the plasmid has a handle the ORF must end right there.
    if _circular_find(seq, FWD_PRIMER) >= 0:
        after = doubled[i + 3 : i + 3 + len(spec.after_stop)]
        if after != spec.after_stop:
            return Finding(
                "single_orf",
                False,
                f"the ORF stops after {codons} codons but not at the pad + 5' handle: "
                "a premature stop or a frameshift",
            )
    if spec.designed_insert and _circular_find(seq, spec.designed_insert) >= 0:
        ins = _circular_find(seq, spec.designed_insert)
        end = (i + 3 - ins) % len(seq)
        if not 0 < end <= len(spec.designed_insert):
            return Finding("single_orf", False, "the ORF runs out of the designed insert")
    return Finding("single_orf", True, f"single ORF, {codons} codons, ending at the handle")


def check_cassette_placement(seq: str, spec: Spec) -> Finding:
    """The tag peptide occurs once, at the recorded residue, in the recorded protein."""
    peptide = spec.cassette_peptide or (
        seqtools.translate(spec.cassette_nt) if spec.cassette_nt else None
    )
    if not peptide:
        return Finding("cassette_placement", True, "no cassette in spec — not checked", "warning")
    orf = read_orf(seq)
    if orf is None:
        return Finding("cassette_placement", False, "no ORF (Kozak ATG to an in-frame stop)")
    protein = seqtools.translate(orf[1][:-3])
    n = protein.count(peptide)
    if n != 1:
        return Finding("cassette_placement", False, f"the tag peptide occurs {n} times in the ORF")
    at = protein.find(peptide)
    if spec.cassette_after_residue is not None and at != spec.cassette_after_residue:
        return Finding(
            "cassette_placement",
            False,
            f"the tag follows residue {at}, but the design puts it after residue "
            f"{spec.cassette_after_residue}",
        )
    if spec.tagged_protein is not None and protein != spec.tagged_protein:
        n_shared = seqtools.longest_shared_tract(protein, spec.tagged_protein)[0]
        return Finding(
            "cassette_placement",
            False,
            f"the ORF encodes a {len(protein)}-aa protein, not the recorded "
            f"{len(spec.tagged_protein)}-aa tagged protein (longest shared stretch {n_shared} aa)",
        )
    return Finding("cassette_placement", True, f"tag once, after residue {at}")


def check_tag_codons(seq: str, spec: Spec) -> Finding:
    """SpyTag003 and P-ALFA-P are on p003's codons, byte for byte (when the spec requires it)."""
    if not spec.require_p003_tag_codons:
        return Finding("tag_codons", True, "not required by the spec — not checked", "warning")
    found = _tandem_in_orf(seq)
    if found is None:
        return Finding("tag_codons", False, "the ORF does not carry the 51-aa tandem cassette once")
    coding, _, c = found
    cassette = coding[3 * c : 3 * (c + len(P003_CASSETTE_PEPTIDE))]
    bad = []
    for name, tag in TAG_PEPTIDES.items():
        a = P003_CASSETTE_PEPTIDE.index(tag)
        lo, hi = 3 * a, 3 * (a + len(tag))
        if cassette[lo:hi] != P003_CASSETTE_NT[lo:hi]:
            diff = sum(
                x != y for x, y in zip(cassette[lo:hi], P003_CASSETTE_NT[lo:hi], strict=True)
            )
            bad.append(f"{name} differs from p003 at {diff} nt")
    if bad:
        return Finding("tag_codons", False, "; ".join(bad))
    return Finding("tag_codons", True, "SpyTag003 and P-ALFA-P on p003's codons")


def _linker_blocks() -> list[tuple[int, int]]:
    """Nucleotide spans of the cassette's linker blocks (everything that is not a tag)."""
    spans, at = [], 0
    for tag in TAG_PEPTIDES.values():
        a = P003_CASSETTE_PEPTIDE.index(tag)
        if a > at:
            spans.append((3 * at, 3 * a))
        at = a + len(tag)
    if at < len(P003_CASSETTE_PEPTIDE):
        spans.append((3 * at, 3 * len(P003_CASSETTE_PEPTIDE)))
    return spans


def check_linkers_not_shared(seq: str, spec: Spec) -> Finding:
    """The cassette in the ORF is not p003's shared DNA (read from the plasmid, not the spec).

    ``diversified_linker`` compares the plasmid with the spec, so a spec that itself names the
    shared encoding passes it. This one reads the cassette out of the ORF and fails if it is
    p003's DNA, unless the spec marks the plasmid as p003 or one of its own arms. A linker block
    that alone equals p003's is a warning: one block in ~1,500 matches by chance.
    """
    found = _tandem_in_orf(seq)
    if found is None:
        return Finding(
            "linkers_not_shared", True, "no 51-aa tandem in the ORF — not checked", "warning"
        )
    coding, _, c = found
    cassette = coding[3 * c : 3 * (c + len(P003_CASSETTE_PEPTIDE))]
    if cassette == P003_CASSETTE_NT:
        if spec.reference_cassette_ok:
            return Finding("linkers_not_shared", True, "p003's cassette DNA, allowed for p003 arms")
        return Finding(
            "linkers_not_shared",
            False,
            "the cassette is p003's shared DNA: its GGGGS linkers are identical in every "
            "construct that carries it",
        )
    same = [
        f"{lo}-{hi}" for lo, hi in _linker_blocks() if cassette[lo:hi] == P003_CASSETTE_NT[lo:hi]
    ]
    if same and not spec.reference_cassette_ok:
        return Finding(
            "linkers_not_shared",
            False,
            f"linker block(s) at cassette nt {', '.join(same)} still on p003's codons",
            "warning",
        )
    return Finding("linkers_not_shared", True, "linkers recoded away from p003's shared DNA")


def _sites_circular(seq: str, enzymes: list[str]) -> list[tuple[str, str, int]]:
    """Type IIS sites on a circular sequence, each once, positions in ``seq``'s numbering."""
    if not enzymes or len(seq) < min(len(seqtools.ENZYME_SITES[e]) for e in enzymes):
        return []
    pad = max(len(seqtools.ENZYME_SITES[e]) for e in enzymes) - 1
    return [h for h in seqtools.find_sites(seq + seq[:pad], enzymes) if h[2] < len(seq)]


def check_type_iis_junctions(seq: str, spec: Spec) -> Finding:
    """No BbsI/BsaI/Esp3I in the insert beyond its own wrapper, and none made at a junction.

    ``golden_gate_domestication`` looks at the insert alone. A site can also be completed by
    the bases either side of the insert-backbone seam, which no check saw; this one counts the
    plasmid's sites against the insert's plus the backbone remainder's, so a seam site shows as
    the difference. An insert that opens with ``CGTCTC`` and ends with ``GAGACG`` is wrapped,
    and those two Esp3I are its only allowed sites.
    """
    if not spec.designed_insert:
        return Finding("type_iis_junctions", True, "no designed insert — not checked", "warning")
    at = _circular_find(seq, spec.designed_insert)
    if at < 0:
        return Finding(
            "type_iis_junctions", True, "insert not found — see insert_matches_design", "warning"
        )
    enzymes = [*ASSEMBLY_ENZYMES, SHUTTLE_ENZYME]
    insert = spec.designed_insert
    wrapped = insert.startswith("CGTCTC") and insert.endswith("GAGACG")
    in_insert = seqtools.find_sites(insert, enzymes)
    allowed = 2 if wrapped else 0
    if len(in_insert) != allowed or (wrapped and any(h[0] != SHUTTLE_ENZYME for h in in_insert)):
        return Finding(
            "type_iis_junctions",
            False,
            f"insert carries {len(in_insert)} BbsI/BsaI/Esp3I site(s), expected {allowed}: "
            f"{in_insert[:4]}",
        )
    doubled = seq + seq
    remainder = doubled[at + len(insert) : at + len(seq)]
    expected = len(in_insert) + len(seqtools.find_sites(remainder, enzymes))
    total = len(_sites_circular(seq, enzymes))
    if total != expected:
        return Finding(
            "type_iis_junctions",
            False,
            f"{total - expected} BbsI/BsaI/Esp3I site(s) completed across an insert-backbone seam",
        )
    return Finding(
        "type_iis_junctions",
        True,
        f"insert has {len(in_insert)} wrapper site(s) only; no site made at either seam",
    )


def _esp3i_cuts(seq: str) -> list[tuple[int, str]]:
    """``(top-strand cut, 4-nt overhang as top-strand sequence)`` for each Esp3I site.

    Esp3I is CGTCTC(1/5): a forward site at p cuts the top strand at p+7 and leaves the 4 nt
    after it as the overhang; a reverse site (GAGACG at p) cuts the top strand at p-5, and the
    overhang is the 4 nt before p-1.
    """
    out = []
    doubled = seq + seq
    for _, strand, p in _sites_circular(seq, [SHUTTLE_ENZYME]):
        if strand == "+":
            cut = (p + 7) % len(seq)
            out.append((cut, doubled[cut : cut + 4]))
        else:
            cut = (p - 5) % len(seq)
            out.append((cut, doubled[cut : cut + 4]))
    return sorted(out)


def check_shuttle_digest(seq: str, spec: Spec) -> Finding:
    """An Esp3I digest releases the insert with two distinct, non-palindromic overhangs.

    ``esp3i_shuttle`` counts the sites. This digests the plasmid in silico: there must be two
    cuts, the fragment between them must hold the ORF and the barcode slot, and its two 4-nt
    overhangs must be different, each non-palindromic, and not each other's reverse
    complement -- otherwise the fragment can ligate to itself or flip in the next vector.
    """
    cuts = _esp3i_cuts(seq)
    if len(cuts) != 2:
        return Finding("shuttle_digest", False, f"Esp3I cuts the plasmid {len(cuts)} time(s)")
    (c1, _), (c2, _) = cuts
    doubled = seq + seq
    fragments = [doubled[c1 : c2 + 4], doubled[c2 : c1 + len(seq) + 4]]
    held = [f for f in fragments if FWD_PRIMER in f and spec.after_barcode in f]
    if len(held) != 1:
        return Finding(
            "shuttle_digest", False, "the barcode handles are not on exactly one digest fragment"
        )
    frag = held[0]
    if "ATG" not in frag[: frag.find(FWD_PRIMER)]:
        return Finding("shuttle_digest", False, "the released fragment carries no ORF start")
    left, right = frag[:4], frag[-4:]
    rc = seqtools.reverse_complement
    problems = []
    if left == right:
        problems.append("the two overhangs are identical")
    if left == rc(right):
        problems.append("the overhangs are each other's reverse complement")
    problems += [f"{o} is palindromic" for o in (left, right) if o == rc(o)]
    if problems:
        return Finding("shuttle_digest", False, "; ".join(problems))
    return Finding(
        "shuttle_digest",
        True,
        f"Esp3I releases a {len(frag)}-nt insert fragment with overhangs {left} / {right}",
    )


def check_backbone_features(seq: str, spec: Spec) -> Finding:
    """Each backbone feature (promoter, IRES, selection, WPRE...) is present exactly once."""
    if not spec.backbone_features:
        return Finding(
            "backbone_features", True, "no backbone features in spec — not checked", "warning"
        )
    bad = []
    for label, feat in spec.backbone_features.items():
        n = _circular_count(seq, feat) + _circular_count(seq, seqtools.reverse_complement(feat))
        if n != 1:
            bad.append(f"{label} x{n}")
    if bad:
        return Finding("backbone_features", False, "not present exactly once: " + ", ".join(bad))
    return Finding(
        "backbone_features", True, f"{len(spec.backbone_features)} features, each present once"
    )


def check_extra_type_iis(seq: str, spec: Spec) -> Finding:
    """AarI, PaqCI, BtgZI and SapI in the insert: reported as a warning, never gated."""
    region = spec.designed_insert or seq
    # p003's own P-ALFA-P DNA carries a SapI site (GAAGAGC); it is in every construct with the
    # tag and was accepted with the tag codons, so the tags' DNA is masked before counting.
    for tag in TAG_PEPTIDES.values():
        a = P003_CASSETTE_PEPTIDE.index(tag)
        tag_nt = P003_CASSETTE_NT[3 * a : 3 * (a + len(tag))]
        region = region.replace(tag_nt, "N" * len(tag_nt))
    hits = []
    for name, site in EXTRA_TYPE_IIS.items():
        rc = seqtools.reverse_complement(site)
        n = region.count(site) + (region.count(rc) if rc != site else 0)
        if n:
            hits.append(f"{name} x{n}")
    if hits:
        return Finding(
            "extra_type_iis", False, "in the insert outside the tags: " + ", ".join(hits), "warning"
        )
    return Finding(
        "extra_type_iis", True, "no AarI, PaqCI, BtgZI or SapI site in the insert outside the tags"
    )


def check_synthesis_windows(
    seq: str,
    spec: Spec,
    *,
    window: int = SYNTHESIS_GC_WINDOW,
    gc_min: float = SYNTHESIS_GC_MIN,
    gc_max: float = SYNTHESIS_GC_MAX,
    max_homopolymer: int = SYNTHESIS_MAX_HOMOPOLYMER,
) -> Finding:
    """Synthesis sanity on the insert, as a warning: homopolymers and 50-nt GC windows.

    Not in :data:`CHECKS`: a vendor's own screen is the real test, and a run that reads only
    a finished map has no reason to fail on it. ``N`` (a barcode slot) counts as half GC, as in
    ``synthesis_compliance``, and does not form homopolymers.
    """
    region = spec.designed_insert or seq
    problems = []
    longest, run, prev = 0, 0, ""
    for ch in region:
        run = run + 1 if ch == prev and ch != "N" else 1
        prev = ch
        longest = max(longest, run if ch != "N" else 0)
    if longest > max_homopolymer:
        problems.append(f"homopolymer of {longest}")
    if len(region) >= window:
        weight = {"G": 1.0, "C": 1.0, "N": 0.5}
        gc = sum(weight.get(ch, 0.0) for ch in region[:window])
        lo = hi = gc / window
        for i in range(window, len(region)):
            gc += weight.get(region[i], 0.0) - weight.get(region[i - window], 0.0)
            lo, hi = min(lo, gc / window), max(hi, gc / window)
        if lo < gc_min:
            problems.append(f"{window}-nt GC down to {lo:.2f} (< {gc_min})")
        if hi > gc_max:
            problems.append(f"{window}-nt GC up to {hi:.2f} (> {gc_max})")
    if problems:
        return Finding("synthesis_windows", False, "; ".join(problems), "warning")
    return Finding("synthesis_windows", True, "GC windows and homopolymers within limits")


def check_annotations(record: Any, spec: Spec) -> list[Finding]:  # a Bio SeqRecord
    """The GenBank features exist AND describe what is actually at those coordinates.

    Presence alone is weak: coordinates copied from an earlier build still carry the right
    labels while pointing at the wrong bases. So each of the three features whose content
    is known independently -- the Kozak, the tag cassette and the barcode -- is read back
    out of the sequence at its own annotated span and compared.
    """
    seq = str(record.seq).upper()
    labels: dict[str, tuple[int, int]] = {}
    for feat in record.features:
        for label in feat.qualifiers.get("label", []):
            labels.setdefault(label, (int(feat.location.start), int(feat.location.end)))

    out: list[Finding] = []
    missing = [lab for lab in REQUIRED_FEATURES if lab not in labels]
    out.append(
        Finding(
            "annotation_present",
            not missing,
            "all required features annotated"
            if not missing
            else f"{len(missing)} feature(s) missing: {missing}",
        )
    )

    expected: list[tuple[str, str | None]] = [
        ("Kozak (AgeI-XbaI)", KOZAK),
        ("barcode", spec.barcode),
        ("ALFA/SpyTag cassette", spec.cassette_nt),
    ]
    wrong = []
    for label, want in expected:
        if not want or label not in labels:
            continue
        lo, hi = labels[label]
        got = seq[lo:hi]
        if got != want:
            wrong.append(
                f"{label}: annotated {hi - lo} nt do not match the expected {len(want)} nt"
            )
    out.append(
        Finding(
            "annotation_matches_sequence",
            not wrong,
            "annotated spans contain what their labels claim" if not wrong else "; ".join(wrong),
        )
    )
    return out


def check_genbank(path: Path, spec: Spec) -> list[Finding]:
    """Every sequence check, plus the annotation checks a GenBank file makes possible."""
    from Bio import SeqIO

    record = SeqIO.read(str(path), "genbank")
    return (
        check_plasmid(str(record.seq), spec)
        + check_annotations(record, spec)
        + [check_rbp_motifs(record)]
    )


def check_rbp_motifs(record: Any) -> Finding:  # a Bio SeqRecord
    """No RNA-binding-protein motif flagged in designed transcript sequence.

    Pumilio PREs, ARE nonamers, CPEs, cryptic poly(A) signals, near-consensus splice sites,
    TGTA/ATTTA clusters and any motif tiled in identical copies, in a UTR, barcode, handle,
    spacer, tiled element or scar of an expression transcript (see
    :mod:`deliverome_analysis.rbp_motifs`). The miR-L no-target arms reached the bench with
    four PREs tiled in their 3'UTR because nothing looked. Hits inside natural elements
    (CDS, IRES, WPRE, LTR) are not this check's business.
    """
    from . import rbp_motifs

    units = rbp_motifs.screen_record(record)
    if not units:
        return Finding("rbp_motifs", True, "no expression transcript found to screen", "warning")
    bad = [
        f"{h.motif} {h.match} at {u.unit.to_map(h.start)} in {h.region_kind} '{h.region}' "
        f"({h.reason})"
        for u in units
        for h in u.hits
        if h.blocking
    ]
    if bad:
        return Finding("rbp_motifs", False, "; ".join(bad))
    n = sum(len(u.hits) for u in units)
    return Finding(
        "rbp_motifs", True, f"{len(units)} transcript(s), {n} motif hits, none in designed sequence"
    )


# Callable[..., Finding]: check_shuttle_sites takes (seq) and the rest (seq, spec);
# check_plasmid dispatches on identity.
CHECKS: tuple[Callable[..., Finding], ...] = (
    check_backbone,
    check_shuttle_sites,
    check_domestication,
    check_insert_matches_design,
    check_diversified_linker,
    check_barcode,
    check_proteoform,
    check_membrane_anchor,
    check_single_orf,
    check_cassette_placement,
    check_tag_codons,
    check_linkers_not_shared,
    check_type_iis_junctions,
    check_shuttle_digest,
    check_backbone_features,
    check_extra_type_iis,
)


def check_plasmid(seq: str, spec: Spec) -> list[Finding]:
    """Run every check. Order-blocking failures are the ones with ``severity='error'``."""
    seq = seq.upper().replace("\n", "")
    out: list[Finding] = []
    for fn in CHECKS:
        try:
            out.append(fn(seq) if fn is check_shuttle_sites else fn(seq, spec))
        except Exception as exc:  # a check that crashes is a failure, not a pass
            # every entry of CHECKS is a def; Callable[..., Finding] just cannot say so
            name = fn.__name__  # ty: ignore[missing-union-attribute]
            out.append(Finding(name, False, f"check raised {type(exc).__name__}: {exc}"))
    return out


def read_genbank_sequence(path: Path) -> str:
    """Sequence from a GenBank file, without requiring the record to be well annotated."""
    from Bio import SeqIO

    return str(SeqIO.read(str(path), "genbank").seq).upper()


def format_report(plasmid: str, findings: list[Finding]) -> str:
    lines = [f"{plasmid}:"]
    for f in findings:
        mark = "ok  " if f.ok else ("WARN" if f.severity == "warning" else "FAIL")
        lines.append(f"  [{mark}] {f.check:<26} {f.detail}")
    bad = [f for f in findings if not f.ok and f.severity == "error"]
    lines.append(
        f"  {len(findings) - len(bad)}/{len(findings)} checks passed"
        + ("" if not bad else f" — {len(bad)} blocking")
    )
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    import argparse

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("genbank", type=Path)
    ap.add_argument("--spec", type=Path, required=True)
    args = ap.parse_args(argv)

    spec = Spec.from_json(args.spec)
    findings = check_genbank(args.genbank, spec)
    print(format_report(spec.plasmid, findings))
    return 1 if any(not f.ok and f.severity == "error" for f in findings) else 0


if __name__ == "__main__":
    raise SystemExit(main())
