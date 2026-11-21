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


def _is_rotation(a: str, b: str) -> bool:
    """True when two circular sequences are the same molecule at a different start."""
    return len(a) == len(b) and b in (a + a)


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
    sites = seqtools.find_sites(seq + seq[:10], [SHUTTLE_ENZYME])
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
    """The cassette carries this construct's OWN encoding, not the shared reference."""
    if not spec.cassette_nt:
        return Finding("diversified_linker", True, "no cassette in spec", "warning")
    if _circular_find(seq, spec.cassette_nt) < 0:
        return Finding(
            "diversified_linker", False, "this construct's issued cassette encoding is not present"
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
    flanked = FWD_PRIMER + spec.barcode + REV_PRIMER
    if _circular_find(seq, flanked) >= 0:
        return Finding("barcode", True, "reserved barcode present, correctly flanked")
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
    if stops == 1:
        return Finding("single_orf", True, f"single ORF, {(i - start - len(KOZAK)) // 3} codons")
    return Finding("single_orf", False, "no in-frame stop codon found")


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
            name = fn.__name__  # ty: ignore[missing-attribute]
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
