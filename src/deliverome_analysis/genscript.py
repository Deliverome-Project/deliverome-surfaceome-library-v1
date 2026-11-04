"""Assess GenScript deliverables against our own construct maps.

A vendor's **map** is a design; a vendor's **sequence export** is what is in the tube.
They are not the same thing, and p004 is the reason this module exists: its design map
and our S3 copy of it agreed byte-for-byte, while both differed from the delivered
plasmid at 11 bases — enough to turn its PacI and ClaI sites into SbfI and XhoI without
anyone noticing. Agreement with a vendor map is not verification.

What GenScript actually ships, and what each kind is worth:

``.seq`` next to ``.ab1``/``.SQD``
    The synthesised **gene fragment** plus Sanger traces. Verifies the insert only —
    typically 5-12% of the plasmid — which is exactly the region p004's error was *not*
    in. :attr:`SequenceKind.FRAGMENT`.
``* full sequence.seq``
    The **whole plasmid**. This is the one that settles a discrepancy.
    :attr:`SequenceKind.FULL_PLASMID`.
``*.final.fa`` with a NanoPlot report
    A **Nanopore assembly** of the whole plasmid. :attr:`SequenceKind.ASSEMBLY`.

A batch download can also contain no sequence at all — only an invoice and the FAQ.
:func:`scan_batch` reports that as an empty result rather than raising, because "the
vendor never sent sequence for this order" is a finding, not an error.

Everything here treats plasmids as **circular**: matching wraps the origin, and
:func:`circular_diff` finds the rotation before comparing, so two records of the same
plasmid that merely start at different bases are reported as identical.

    from deliverome_analysis import genscript
    found = genscript.scan_batch(Path.home() / "Downloads" / "U367PTSPG0-batchDownload-8ad960")
    report = genscript.audit(found, maps)
"""

from __future__ import annotations

import re
from dataclasses import asdict, dataclass
from enum import StrEnum
from pathlib import Path

from Bio import SeqIO

from deliverome_analysis.seqtools import reverse_complement

__all__ = [
    "MapMatch",
    "SequenceKind",
    "VendorSequence",
    "audit",
    "circular_diff",
    "match_to_maps",
    "read_sequence",
    "scan_batch",
]

_DNA = re.compile(r"[^ACGT]")
_ORDER_ITEM = re.compile(r"(U[0-9A-Z]{7,9}G0)[-_](\d+)")
SEQUENCE_SUFFIXES = (".seq", ".fa", ".fasta", ".gb", ".gbk")


class SequenceKind(StrEnum):
    """How much of a plasmid a vendor file actually covers."""

    FULL_PLASMID = "full_plasmid"
    ASSEMBLY = "assembly"
    FRAGMENT = "fragment"


@dataclass(frozen=True)
class VendorSequence:
    """One sequence file shipped by the vendor."""

    path: Path
    sequence: str
    kind: SequenceKind
    order: str | None = None
    item: str | None = None

    @property
    def length(self) -> int:
        return len(self.sequence)

    @property
    def label(self) -> str:
        return f"{self.order}-{self.item}" if self.order and self.item else self.path.stem


@dataclass(frozen=True)
class MapMatch:
    """The result of checking one vendor sequence against one construct map."""

    vendor: str
    map_name: str
    strand: str
    coverage: float
    diffs: tuple[tuple[int, str, str], ...] = ()

    @property
    def exact(self) -> bool:
        return not self.diffs

    def to_dict(self) -> dict:
        return asdict(self) | {"exact": self.exact, "n_diffs": len(self.diffs)}


def read_sequence(path: Path) -> str:
    """Read a vendor sequence file to a bare uppercase ACGT string.

    Handles the four shapes GenScript ships: a raw ``.seq`` dump (often mixed case and
    with no header at all), FASTA, and GenBank. Non-ACGT characters are dropped, so a
    header line or soft-masking never silently shifts coordinates.
    """
    suffix = path.suffix.lower()
    if suffix in (".gb", ".gbk"):
        return str(SeqIO.read(path, "genbank").seq).upper()
    if suffix in (".fa", ".fasta"):
        return str(SeqIO.read(path, "fasta").seq).upper()
    text = path.read_text()
    if text.lstrip().startswith(">"):
        text = "\n".join(text.splitlines()[1:])
    return _DNA.sub("", text.upper())


def _classify(path: Path) -> SequenceKind:
    name = path.name.lower()
    if "full sequence" in name:
        return SequenceKind.FULL_PLASMID
    if name.endswith((".final.fa", ".final.fasta")):
        return SequenceKind.ASSEMBLY
    return SequenceKind.FRAGMENT


def scan_batch(root: Path) -> list[VendorSequence]:
    """Find every sequence file in a GenScript batch download.

    Returns an empty list when the download carries no sequence — which happens, and is
    worth surfacing rather than treating as a crash: the U4744625G0 batch download
    contains only an invoice and the Gene FAQ.
    """
    if not root.exists():
        raise FileNotFoundError(root)
    out: list[VendorSequence] = []
    for path in sorted(root.rglob("*")):
        if not path.is_file() or path.suffix.lower() not in SEQUENCE_SUFFIXES:
            continue
        seq = read_sequence(path)
        if not seq:
            continue
        m = _ORDER_ITEM.search(path.name) or _ORDER_ITEM.search(path.parent.name)
        out.append(
            VendorSequence(
                path=path,
                sequence=seq,
                kind=_classify(path),
                order=m.group(1) if m else None,
                item=m.group(2) if m else None,
            )
        )
    return out


def circular_diff(a: str, b: str) -> list[tuple[int, str, str]] | None:
    """Per-base differences between two circular sequences of the same length.

    Finds the rotation of ``b`` that aligns with ``a`` before comparing, so records of
    the same plasmid that merely start at different bases compare as identical.
    Positions are 1-based in ``a``. Returns ``None`` when the lengths differ or no
    rotation aligns — meaning these are not the same construct and a base-level diff
    would be meaningless.
    """
    a, b = a.upper(), b.upper()
    if len(a) != len(b) or not a:
        return None
    anchor = a[: min(40, len(a))]
    i = (b + b).find(anchor)
    if i < 0:
        return None
    rotated = (b + b)[i : i + len(b)]
    return [(k + 1, a[k], rotated[k]) for k in range(len(a)) if a[k] != rotated[k]]


def match_to_maps(
    vendor: VendorSequence, maps: dict[str, str], *, max_diffs: int = 0
) -> list[MapMatch]:
    """Every map this vendor sequence is contained in, or differs from only slightly.

    Containment is checked on both strands across the origin. A full-plasmid export is
    additionally compared base by base against same-length maps, so a near-match — the
    p004 case, 11 bases out of 6,201 — is reported with its differences rather than
    silently counting as "no match".
    """
    hits: list[MapMatch] = []
    for name, ref in sorted(maps.items()):
        ref = ref.upper()
        doubled = ref + ref
        for strand, probe in (
            ("fwd", vendor.sequence),
            ("rev", reverse_complement(vendor.sequence)),
        ):
            if probe in doubled:
                hits.append(MapMatch(vendor.label, name, strand, len(vendor.sequence) / len(ref)))
                break
        else:
            if len(ref) == len(vendor.sequence):
                diffs = circular_diff(ref, vendor.sequence)
                if diffs is not None and 0 < len(diffs) <= max_diffs:
                    hits.append(MapMatch(vendor.label, name, "fwd", 1.0, tuple(diffs)))
    return hits


def audit(found: list[VendorSequence], maps: dict[str, str], *, max_diffs: int = 0) -> list[dict]:
    """One flat record per vendor sequence: what it matched, and how much it proves."""
    rows = []
    for vendor in found:
        hits = match_to_maps(vendor, maps, max_diffs=max_diffs)
        rows.append(
            {
                "vendor": vendor.label,
                "file": vendor.path.name,
                "kind": vendor.kind.value,
                "length_nt": vendor.length,
                "matches": ",".join(h.map_name for h in hits),
                "n_matches": len(hits),
                "coverage": round(max((h.coverage for h in hits), default=0.0), 4),
                "exact": all(h.exact for h in hits) if hits else False,
                "n_diffs": sum(len(h.diffs) for h in hits),
                "verdict": _verdict(vendor, hits),
            }
        )
    return rows


def _verdict(vendor: VendorSequence, hits: list[MapMatch]) -> str:
    if not hits:
        return "NO MATCH — investigate"
    if len(hits) > 1:
        return f"region shared by {len(hits)} maps"
    hit = hits[0]
    if hit.diffs:
        return f"{len(hit.diffs)} base difference(s) — the map is stale"
    # Judge by what the sequence actually covers, not by the file's name: a pLann
    # re-annotation of a Nanopore assembly arrives as a .gbk and still covers the
    # whole plasmid, while a "full sequence" export of a sub-region would not.
    if hit.coverage >= 0.99:
        return "whole plasmid verified"
    return f"insert verified ({hit.coverage:.0%} of the plasmid); backbone unverified"
