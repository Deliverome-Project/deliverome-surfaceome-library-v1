"""RNA-binding-protein and RNA-processing motifs on the transcribed strand of a construct.

Construct checks already screen restriction sites, poly(A) signals inside barcodes, splice
sites inside the barcode pool and human miRNA seeds. None of them looks for motifs that a
cellular RNA-binding protein reads to lower expression. The miR-L no-target arms
(P-0017/P-0019/P-0021) show why that matters: their scramble block tiles
``TACTTGTAAATACTATTCTAT`` four times, every copy carries a perfect Pumilio response
element, and so the 3'UTR between barcode and WPRE has four PREs. Those arms express
~2-3x lower than their target-site partners in S6E5. The PREs are a likely cause, not a
proven one. The general lesson is that tiling multiplies a single chance hit, so this
module reports copy number as well as presence.

Motifs (DNA alphabet, sense strand of the transcript):

======================  ==========================  =========================================
name                    pattern                     read by / effect
======================  ==========================  =========================================
``PRE``                 ``TGTA[ACT]ATA``            Pumilio 1/2 (UGUAHAUA); decay, repression
``CFIm_TGTA``           ``TGTA``                    CFIm25/NUDT21; clusters (>=2 within 50 nt)
                                                    recruit alternative 3'-end processing
``ARE_pentamer``        ``ATTTA``                   AU-rich element core; flagged in clusters
``ARE_nonamer``         ``TTATTTA[AT][AT]``         class II ARE (TTP, HuR); decay
``CPE``                 ``TTTTA(A)T``               CPEB; poly(A) length control near a PAS
``PAS``                 ``AATAAA`` / ``ATTAAA``     cryptic cleavage/polyadenylation
``splice_donor``        ``GT[AG]AG``                5' splice site core
``splice_acceptor``     pyrimidine tract + ``YAG``  3' splice site (heuristic, see below)
======================  ==========================  =========================================

The acceptor screen is a heuristic: a ``[CT]AG`` whose 15 preceding nucleotides hold at
least 12 pyrimidines. It ignores the branch point, so it over-calls; treat its hits as
places to look.

:func:`scan` takes a sequence (already oriented as the transcript) plus annotated
:class:`Region` spans and returns every :class:`Hit` with its position, the region it falls
in, and its copy number. Two hits of the same motif are *copies* when they sit in at least
``min_shared`` nt of identical sequence (motif plus flanks), which is what tandem tiling
produces and what independent chance hits almost never do.

:func:`transcription_units` finds the Pol II transcripts in a GenBank plasmid map (internal
promoter to the 3' LTR R end, or to the same-strand poly(A) signal) and
:func:`screen_record` runs :func:`scan` over each, labelling hits with the map's own
features. Whether a flagged hit should block an order is decided by region kind: hits in
designed sequence (UTRs, barcodes, handles, spacers, tiled elements, cloning scars) block,
hits inside natural elements (CDS, IRES, WPRE, LTR, introns) are reported only.
"""

from __future__ import annotations

import re
from collections.abc import Callable, Iterable, Sequence
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, cast

# Biopython SeqRecord / SeqFeature arguments below are typed Any: Bio annotates
# SeqFeature.location as Optional, but every feature parsed from a GenBank map has one, and
# narrowing it at each use would add checks this code has never needed.

# -- who chose what ----------------------------------------------------------------------
#
# Becca, 2026-11-10: "for any manual fixes please note the agent who made the fix /
# adjustment and reasoning". The motif panel was Becca's request; every cutoff and rule below
# that she did not specify was chosen by the agent that implemented the screen. ``DECISIONS``
# carries that record so outputs can quote it (``l1_rbp_motif_summary.json`` does).

IMPLEMENTED_BY = "Claude Code session rbp-motif-screen (worktree blissful-hugle-cab8fa), PR #277"
_BECCA = "Becca, 2026-11-08 (request relayed to the rbp-motif-screen session on 2026-11-10)"
_AGENT = "Claude Code session rbp-motif-screen, 2026-11-10; no separate sign-off, shipped in PR #277 (merged by Becca 2026-11-10)"

#: name -> (decided_by, implemented_by, reason)
DECISIONS: dict[str, tuple[str, str, str]] = {
    "motif_panel": (
        _BECCA,
        IMPLEMENTED_BY,
        "PRE TGTA[ACT]ATA, CFIm25 TGTA clusters, ARE ATTTA and TTATTTA[AT][AT], CPE TTTTAT/"
        "TTTTAAT, PAS AATAAA/ATTAAA, splice donor GT[AG]AG and a simple acceptor screen, as "
        "requested; motivated by the four tiled PREs in the miR-L no-target arms",
    ),
    "cluster_window_50nt": (
        _BECCA,
        IMPLEMENTED_BY,
        "Becca asked to flag >=2 TGTA within ~50 nt; the agent applied the same window to ATTTA",
    ),
    "min_shared_20nt": (
        _AGENT,
        IMPLEMENTED_BY,
        "20 nt of identical context (motif plus flanks) separates tiled copies from chance "
        "repeats: the miR-L scramble unit is 21 nt, and two random 20-mers match at ~1e-12",
    ),
    "acceptor_12_of_15_pyrimidines": (
        _AGENT,
        IMPLEMENTED_BY,
        "Becca asked for 'a simple acceptor screen'; a YAG after >=12 pyrimidines in 15 nt is a "
        "polypyrimidine-tract heuristic with no branch point, chosen to be simple and to over-call",
    ),
    "strong_donor_MAG_GTRAGT": (
        _AGENT,
        IMPLEMENTED_BY,
        "a bare GT[AG]AG occurs about every 500 nt, so only near-consensus donors are flagged",
    ),
    "ltr_pas_to_cleavage_25nt": (
        _AGENT,
        IMPLEMENTED_BY,
        "the transcript ends at the HIV-1 R/U5 boundary, AATAAA + GCTTGCCTTGAGTGCTTCA, read off "
        "the p0039 3' LTR sequence",
    ),
    "block_only_designed_sequence": (
        _AGENT,
        IMPLEMENTED_BY,
        "natural CDSs and backbone elements (IRES2, WPRE, LTR) carry motifs that a gate cannot "
        "act on; Becca's request already separated ORF-intrinsic from shared-element hits",
    ),
    "dropout_never_blocks": (
        _AGENT,
        IMPLEMENTED_BY,
        "the p0039 stuffer is cut out during cloning and never reaches a finished construct",
    ),
}


# -- motifs ------------------------------------------------------------------------------

#: Nucleotides upstream of a ``YAG`` examined for the acceptor's pyrimidine tract, and how
#: many of them must be C/T.
ACCEPTOR_TRACT = 15
ACCEPTOR_MIN_PYRIMIDINES = 12
#: Window for TGTA / ATTTA clusters (start-to-start distance, nt).
CLUSTER_WINDOW = 50
#: Identical sequence two hits must share to count as copies of one unit.
MIN_SHARED = 20
#: Nucleotides from a terminal AATAAA to the cleavage site at the HIV-1 R/U5 boundary
#: (AATAAA + GCTTGCCTTGAGTGCTTCA).
LTR_PAS_TO_CLEAVAGE = 25


@dataclass(frozen=True)
class Motif:
    """One motif. ``pattern`` is a regex on the sense strand; ``finder`` overrides it."""

    name: str
    family: str
    pattern: str
    note: str
    #: Flag every hit (True) or only hits that cluster or tile (False).
    always_flag: bool
    #: Flag when at least two hits start within ``CLUSTER_WINDOW`` nt of each other.
    cluster: bool = False
    finder: Callable[[str], list[tuple[int, int]]] | None = None

    def find(self, seq: str) -> list[tuple[int, int]]:
        if self.finder is not None:
            return self.finder(seq)
        rx = re.compile(f"(?=({self.pattern}))")  # lookahead: overlapping hits
        return [(m.start(), m.start() + len(m.group(1))) for m in rx.finditer(seq)]


def _find_acceptors(seq: str) -> list[tuple[int, int]]:
    out = []
    for m in re.finditer(r"(?=[CT]AG)", seq):
        p = m.start()
        if p < ACCEPTOR_TRACT:
            continue
        tract = seq[p - ACCEPTOR_TRACT : p]
        if sum(c in "CT" for c in tract) >= ACCEPTOR_MIN_PYRIMIDINES:
            out.append((p - ACCEPTOR_TRACT, p + 3))
    return out


#: A donor whose exon side is ``[AC]AG`` and intron side ``GT[AG]AGT`` is close to the
#: consensus; a bare ``GT[AG]AG`` turns up about every 500 nt and is reported unflagged.
STRONG_DONOR = re.compile(r"[AC]AGGT[AG]AGT")

MOTIFS: tuple[Motif, ...] = (
    Motif("PRE", "Pumilio", "TGTA[ACT]ATA", "Pumilio response element UGUAHAUA", True),
    Motif(
        "CFIm_TGTA",
        "CFIm25",
        "TGTA",
        "CFIm25/NUDT21 UGUA; flagged when >=2 lie within 50 nt",
        False,
        cluster=True,
    ),
    Motif(
        "ARE_pentamer",
        "ARE",
        "ATTTA",
        "AU-rich element core AUUUA; flagged when >=2 lie within 50 nt",
        False,
        cluster=True,
    ),
    Motif("ARE_nonamer", "ARE", "TTATTTA[AT][AT]", "class II ARE nonamer UUAUUUAWW", True),
    Motif("CPE_TTTTAT", "CPEB", "TTTTAT", "cytoplasmic polyadenylation element UUUUAU", True),
    Motif("CPE_TTTTAAT", "CPEB", "TTTTAAT", "cytoplasmic polyadenylation element UUUUAAU", True),
    Motif("PAS_AATAAA", "PAS", "AATAAA", "canonical poly(A) signal", True),
    Motif("PAS_ATTAAA", "PAS", "ATTAAA", "second-commonest poly(A) signal", True),
    Motif(
        "splice_donor",
        "splice",
        "GT[AG]AG",
        "5' splice-site core GURAG; flagged when it matches MAG|GURAGU",
        False,
    ),
    Motif(
        "splice_acceptor",
        "splice",
        "",
        f"YAG after >={ACCEPTOR_MIN_PYRIMIDINES}/{ACCEPTOR_TRACT} pyrimidines (heuristic)",
        True,
        finder=_find_acceptors,
    ),
)
MOTIF_NAMES = tuple(m.name for m in MOTIFS)

# -- regions -----------------------------------------------------------------------------

#: Region kinds that are designed sequence. A flagged hit in one of these can be designed
#: away, so it blocks an order unless accepted with a reason.
DESIGNED_KINDS = frozenset(
    {"5'UTR", "3'UTR", "intergenic", "barcode", "handle", "spacer", "tiled", "scar", "kozak"}
)
#: Natural elements: reported, never blocking (a CDS, the IRES, the WPRE carry what they carry).
#: A ``dropout`` stuffer is cut out during cloning and never reaches a finished construct, so
#: it is reported and never blocks either.
NATURAL_KINDS = frozenset(
    {"cds", "ires", "wpre", "ltr", "intron", "exon", "promoter", "polyA", "miRNA", "element"}
)


@dataclass(frozen=True)
class Region:
    """A labelled span of the transcript, 0-based half-open."""

    start: int
    end: int
    label: str
    kind: str


@dataclass
class Hit:
    motif: str
    family: str
    start: int  # 0-based, transcript coordinates
    end: int
    match: str
    region: str = ""
    region_kind: str = ""
    copy_group: int = 0  # 0 = a single copy
    copy_number: int = 1
    cluster_size: int = 1  # same-motif hits within CLUSTER_WINDOW, this one included
    nearest_same_nt: int | None = None  # start-to-start distance to the nearest same motif
    expected: bool = False  # the transcript's own terminal poly(A) signal
    flagged: bool = False
    reason: str = ""

    @property
    def designed(self) -> bool:
        return self.region_kind in DESIGNED_KINDS

    @property
    def blocking(self) -> bool:
        return self.flagged and self.designed

    def as_row(self) -> dict[str, Any]:
        row = asdict(self)
        row["blocking"] = self.blocking
        return row


def _region_for(start: int, end: int, regions: Sequence[Region]) -> Region | None:
    """The most specific (shortest) region overlapping [start, end); None if none does."""
    best = None
    for r in regions:
        if r.start < end and start < r.end:
            if best is None or (r.end - r.start) < (best.end - best.start):
                best = r
    return best


def _shared_context(seq: str, a: int, b: int, length: int) -> int:
    """Length of identical sequence around two hits of the same ``length`` at a < b."""
    left = 0
    while a - left - 1 >= 0 and seq[a - left - 1] == seq[b - left - 1]:
        left += 1
    right = 0
    n = len(seq)
    while b + length + right < n and seq[a + length + right] == seq[b + length + right]:
        right += 1
    return left + length + right


def _copy_groups(seq: str, hits: list[Hit], min_shared: int) -> dict[int, list[int]]:
    """Union hits of one motif that share >= ``min_shared`` nt of identical context."""
    parent = list(range(len(hits)))

    def root(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i in range(len(hits)):
        for j in range(i + 1, len(hits)):
            a, b = hits[i], hits[j]
            if a.match != b.match or a.start == b.start:
                continue
            if _shared_context(seq, a.start, b.start, a.end - a.start) >= min_shared:
                parent[root(j)] = root(i)
    groups: dict[int, list[int]] = {}
    for i in range(len(hits)):
        groups.setdefault(root(i), []).append(i)
    return groups


def scan(
    seq: str,
    regions: Iterable[Region] = (),
    *,
    motifs: Sequence[Motif] = MOTIFS,
    expected_pas: Iterable[int] = (),
    min_shared: int = MIN_SHARED,
    cluster_window: int = CLUSTER_WINDOW,
) -> list[Hit]:
    """Every motif hit in ``seq`` (the transcript, sense strand), sorted by position.

    ``expected_pas`` are start positions of the transcript's intended poly(A) signal(s);
    those hits are reported with ``expected=True`` and never flagged.
    """
    seq = seq.upper().replace("U", "T")
    regions = list(regions)
    expected = set(expected_pas)
    out: list[Hit] = []
    next_group = 1
    for motif in motifs:
        hits = []
        for s, e in motif.find(seq):
            # An acceptor is placed by its YAG (the 3' splice site), not by its tract.
            r = _region_for(e - 3 if motif.finder is _find_acceptors else s, e, regions)
            hits.append(
                Hit(
                    motif.name,
                    motif.family,
                    s,
                    e,
                    seq[s:e],
                    region=r.label if r else "",
                    region_kind=r.kind if r else "",
                )
            )
        if not hits:
            continue
        starts = [h.start for h in hits]
        for h in hits:
            others = [abs(x - h.start) for x in starts if x != h.start]
            h.nearest_same_nt = min(others) if others else None
            h.cluster_size = sum(abs(x - h.start) <= cluster_window for x in starts)
        groups = _copy_groups(seq, hits, min_shared)
        for members in sorted(groups.values()):
            if len(members) < 2:
                continue
            for i in members:
                hits[i].copy_group = next_group
                hits[i].copy_number = len(members)
            next_group += 1
        for h in hits:
            reasons = []
            if motif.family == "PAS" and h.start in expected:
                h.expected = True
            elif motif.always_flag:
                reasons.append(motif.name)
            if motif.cluster and h.cluster_size >= 2:
                reasons.append(f"cluster of {h.cluster_size} within {cluster_window} nt")
            if motif.name == "splice_donor" and STRONG_DONOR.match(seq, max(h.start - 3, 0)):
                if h.start >= 3:
                    reasons.append("near-consensus donor MAG|GTRAGT")
            if h.copy_number >= 2 and not h.expected:
                reasons.append(f"tiled x{h.copy_number}")
            h.flagged = bool(reasons)
            h.reason = "; ".join(reasons)
        out.extend(hits)
    out.sort(key=lambda h: (h.start, h.motif))
    return out


# -- GenBank maps ------------------------------------------------------------------------

#: Promoters that do not make the expression transcript in transduced cells: prokaryotic
#: and phage promoters, and RSV, which drives the genomic RNA in packaging cells.
NON_EXPRESSION_PROMOTERS = ("lac", "t3", "t7", "sp6", "ampr", "bla", "rsv", "cap binding")

#: Feature labels that annotate a sequencing or PCR primer rather than a construct element.
#: Kept out of region labelling so a hit is attributed to the element, not to "WPRE-R".
_ELEMENT_PRIMER = re.compile(r"handle|gap \+ rt|\brt\b|capture|cs1|o-\d+|jt_fw|revprimer", re.I)


def classify(feature_type: str, label: str) -> str | None:
    """Region kind for a GenBank feature, or None to ignore the feature."""
    lab = label.lower()
    if feature_type in ("source", "rep_origin", "protein_bind"):
        return None
    if feature_type in ("primer_bind", "primer"):
        return "handle" if _ELEMENT_PRIMER.search(label) else None
    if feature_type in ("CDS", "sig_peptide", "mat_peptide") or lab == "stop":
        return "cds"
    if "kozak" in lab:
        return "kozak"
    if "stuffer" in lab or "dropout" in lab:
        return "dropout"
    if _SCAR.search(lab):
        return "scar"
    if re.search(r"\bsites? \d+ of \d+|scramble|target array|\bsites \[|\barray\b", lab):
        return "tiled"
    if _ELEMENT_PRIMER.search(label):
        return "handle"
    if "barcode" in lab or re.search(r"\bbc\b|-bc-|_bc\b", lab):
        return "barcode"
    if lab == "pad" or "spacer" in lab:
        return "spacer"
    if "ires" in lab:
        return "ires"
    if "wpre" in lab:
        return "wpre"
    if feature_type == "LTR" or "ltr" in lab:
        return "ltr"
    if feature_type == "intron" or "intron" in lab:
        return "intron"
    if feature_type == "exon":
        return "exon"
    if feature_type == "promoter" or "promoter" in lab or "enhancer" in lab:
        return "promoter"
    if feature_type == "polyA_signal" or "poly(a)" in lab:
        return "polyA"
    if feature_type in ("misc_RNA", "ncRNA") or "mir" in lab:
        return "miRNA"
    if feature_type == "misc_binding":
        return None
    return "element"


#: Cloning scars and restriction-site annotations, matched before "site" or "array" so that
#: "AsiSI site, flanking the array" is a scar and not a tiled element.
_SCAR = re.compile(
    r"\b(esp3i|bsai|bbsi|asisi|pmei|asci|nhei|xhoi|sbfi|agei|xbai|bstbi|pmli|ecori|bamhi)\b"
    r"|\bmcs\b|junction|shuttle"
)


@dataclass
class TranscriptionUnit:
    name: str
    strand: int
    start: int  # plasmid coordinates, 0-based half-open, start < end
    end: int
    seq: str  # sense strand of the transcript
    regions: list[Region] = field(default_factory=list)
    expected_pas: list[int] = field(default_factory=list)

    def to_map(self, tx_pos: int) -> int:
        """1-based plasmid position of transcript position ``tx_pos`` (0-based)."""
        return self.start + tx_pos + 1 if self.strand == 1 else self.end - tx_pos


def _label(feature: Any) -> str:
    q = feature.qualifiers
    return cast(
        str, (q.get("label") or q.get("note") or q.get("gene") or q.get("product") or [""])[0]
    )


def _is_expression_promoter(feature: Any) -> bool:
    if feature.type != "promoter":
        return False
    lab = _label(feature).lower()
    return not any(k in lab for k in NON_EXPRESSION_PROMOTERS)


def transcription_units(record: Any) -> list[TranscriptionUnit]:
    """Pol II expression transcripts in a plasmid map.

    Forward: from the end of the first expression promoter inside the lentiviral cassette
    (or anywhere, without LTRs) to the cleavage site 25 nt past the 3' LTR's AATAAA (the R/U5
    boundary), or to the end of the next forward poly(A) signal. Reverse-strand promoters run
    to the next reverse poly(A) signal. Downstream promoters on the same strand nest inside
    the first transcript and are labelled as regions rather than given their own unit.
    """
    seq = str(record.seq).upper()
    feats = list(record.features)
    ltrs = sorted(
        (f for f in feats if f.type == "LTR" or "ltr" in _label(f).lower()),
        key=lambda f: int(f.location.start),
    )
    five_ltr_end = int(ltrs[0].location.end) if len(ltrs) >= 2 else 0
    three_ltr = ltrs[-1] if len(ltrs) >= 2 else None
    polyas = [f for f in feats if f.type == "polyA_signal"]
    promoters = sorted(
        (f for f in feats if _is_expression_promoter(f)), key=lambda f: int(f.location.start)
    )
    units: list[TranscriptionUnit] = []

    fwd = [
        p
        for p in promoters
        if p.location.strand == 1
        and int(p.location.start) >= five_ltr_end
        and (three_ltr is None or int(p.location.end) <= int(three_ltr.location.start))
    ]
    if fwd:
        p = fwd[0]
        start = int(p.location.end)
        expected: list[int] = []
        if three_ltr is not None:
            ls, le = int(three_ltr.location.start), int(three_ltr.location.end)
            pas = seq.rfind("AATAAA", ls, le)
            end = min(pas + LTR_PAS_TO_CLEAVAGE, le) if pas >= 0 else le
            if pas >= 0:
                expected.append(pas - start)
        else:
            down = sorted(
                (f for f in polyas if f.location.strand == 1 and int(f.location.start) >= start),
                key=lambda f: int(f.location.start),
            )
            end = int(down[0].location.end) if down else start
            pas = seq.find("AATAAA", int(down[0].location.start), end) if down else -1
            if pas >= 0:
                expected.append(pas - start)
        if end > start:
            units.append(_unit(record, seq, feats, f"{_label(p)} (+)", 1, start, end, expected))

    for p in (q for q in promoters if q.location.strand == -1):
        end = int(p.location.start)
        down = [f for f in polyas if f.location.strand == -1 and int(f.location.end) <= end]
        if not down:
            continue
        pa = max(down, key=lambda f: int(f.location.end))
        start = int(pa.location.start)
        # transcript coordinates run from the promoter end down to the poly(A) feature
        pas = seq.rfind("TTTATT", start, int(pa.location.end))
        expected = [end - (pas + 6)] if pas >= 0 else []
        units.append(_unit(record, seq, feats, f"{_label(p)} (-)", -1, start, end, expected))
    return units


def _unit(
    record: Any,
    seq: str,
    feats: list[Any],
    name: str,
    strand: int,
    start: int,
    end: int,
    expected: list[int],
) -> TranscriptionUnit:
    from deliverome_analysis.cloning.seqtools import reverse_complement

    sub = seq[start:end] if strand == 1 else reverse_complement(seq[start:end])
    tu = TranscriptionUnit(name, strand, start, end, sub, expected_pas=expected)
    cds_spans = []
    for f in feats:
        kind = classify(f.type, _label(f))
        if kind is None:
            continue
        fs, fe = int(f.location.start), int(f.location.end)
        if fe <= start or fs >= end:
            continue
        fs, fe = max(fs, start), min(fe, end)
        if strand == 1:
            ts, te = fs - start, fe - start
        else:
            ts, te = end - fe, end - fs
        tu.regions.append(Region(ts, te, _label(f), kind))
        if kind == "cds" and f.location.strand == strand:
            cds_spans.append((ts, te))
    # Unannotated stretches: 5'UTR before the first same-strand CDS, 3'UTR after the last.
    n = len(sub)
    first = min((s for s, _ in cds_spans), default=n)
    last = max((e for _, e in cds_spans), default=n)
    tu.regions.append(Region(0, first, "5'UTR (unannotated)", "5'UTR"))
    if last < n:
        tu.regions.append(Region(last, n, "3'UTR (unannotated)", "3'UTR"))
    if first < last:
        tu.regions.append(Region(first, last, "between CDSs (unannotated)", "intergenic"))
    # The broad unannotated spans must lose to any real feature: make them longest.
    tu.regions.sort(key=lambda r: r.end - r.start)
    return tu


@dataclass
class UnitHits:
    unit: TranscriptionUnit
    hits: list[Hit]


def screen_record(record: Any, **kwargs: Any) -> list[UnitHits]:
    """:func:`scan` every expression transcript in a Biopython GenBank record."""
    out = []
    for tu in transcription_units(record):
        out.append(UnitHits(tu, scan(tu.seq, tu.regions, expected_pas=tu.expected_pas, **kwargs)))
    return out


def screen_genbank(path: str | Path) -> list[UnitHits]:
    from Bio import SeqIO

    return screen_record(SeqIO.read(str(path), "genbank"))


def summarize(hits: Iterable[Hit]) -> dict[str, int]:
    """Hit count per motif (every motif present, zero when absent)."""
    counts = dict.fromkeys(MOTIF_NAMES, 0)
    for h in hits:
        counts[h.motif] += 1
    return counts
