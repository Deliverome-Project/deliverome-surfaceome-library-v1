"""How often the DeepTMHMM v1 call that placed each L1 tag disagrees with experiment.

Every L1 member gets a tag placement from the DeepTMHMM v1 topology of its own designed
protein (``l1_tag_paths.tsv``): N-terminal after the signal peptide, N-terminal at
residue 1 (multipass or type III), C-terminal, or an internal loop. Each of those puts
the cassette on a segment the predictor called extracellular. This script asks, for
every member whose protein can be tied to an experimentally supported topology, which
side of the membrane experiment puts that segment on.

**Evidence.** Two sources, both pinned on ``external://``:

* UniProtKB/Swiss-Prot human and mouse entries with topology features, release 2026_04
  (:data:`UNIPROT_KEY`). A ``TOPO_DOM`` / ``TRANSMEM`` / ``INTRAMEM`` feature counts as
  experimental only if one of its evidence codes is ECO:0000269 (experimental evidence,
  manual assertion; usually a PubMed citation) or ECO:0007744 (structure-derived, manual
  assertion; usually a PDB citation). ECO:0000305 (curator inference), ECO:0000250 (by
  similarity) and ECO:0000255 (sequence analysis) are predictions and are not used.
* TOPDB v3.3 (:data:`TOPDB_KEY`), which includes the PDBTM structure-derived sides. Only
  its ``Experiments`` block is read -- the per-region experimental constraints (structure,
  glycosylation, reporter fusions, chemical labelling, epitope insertion, proteolysis) --
  never its ``Topology`` block, which is a CCTOP prediction constrained by those
  experiments. Experiments of type ``Other`` (annotated linear motifs, sequence motifs,
  tailoring) are not experiments on the protein and are dropped.

**Sequence, not accession.** A member's evidence is matched by protein sequence. An entry
whose sequence is identical to the member's is an ``exact`` match. Otherwise the member is
aligned (global, BLOSUM62) to the same gene's evidence entries, human and mouse, and the
experimental features are carried across aligned residues only (``mapped``, identity at
least :data:`MIN_IDENTITY`). A mapped call is kept only if the member has the same
membrane-crossing context upstream of the site: no annotated TM helix of the entry lost
from the member, and no predicted helix of the member absent from the entry, before the
site. Losing or gaining a helix upstream flips everything after it, so a canonical-entry
side would not describe the isoform there.

**Site and segment.** The tagged site is residue 1 (N at residue 1), the first mature
residue after the cassette (N after SP: ``native_after_residue + 1`` from
``constructs_l1_nterm.tsv``, else the v1 signal peptide end + 1), the last residue
(C-terminal), or ``cassette_after_residue`` for each internal arm in
``constructs_l1_internal.tsv`` (snorkel arms and held members carry no loop site). The
side of the site is the side of the whole membrane-free segment that contains it,
bounded by the member's v1 helices and signal peptide and by every mapped UniProt
TRANSMEM/INTRAMEM feature of any evidence level; any experimental region overlapping that
segment calls its side. A site inside an experimental TM region is called ``membrane``.

**Partial annotations.** UniProt and TOPDB often support only some segments of a protein
experimentally, and give the rest by similarity or prediction. A site is called only from
experimental features that overlap its own segment -- sides are never propagated across
helices from an experimental segment elsewhere -- so a partially annotated entry
contributes evidence where it has it and ``none`` elsewhere. Segments with experimental
regions on both sides are ``conflict`` and are counted separately, not as either verdict.

Verdict per site: ``concordant`` (experiment: outside), ``discordant`` (inside or in the
membrane), ``conflict``, or ``no_evidence``. A member is discordant if any of its tagged
sites is. Also written per member: the experimental side of the N-terminus (mature, after
the v1 signal peptide) and of the C-terminus, against v1's, so a placement chosen on a
wrong N or C call shows even when the tag site itself is fine.

    uv run --locked python topology_vs_experiment.py            # analysis
    uv run --locked python topology_vs_experiment.py --fetch    # re-snapshot evidence
"""

from __future__ import annotations

import argparse
import collections
import csv
import gzip
import hashlib
import json
import re
import shutil
import tempfile
import urllib.request
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd
from Bio import Align
from Bio.Align import substitution_matrices

from deliverome_analysis import library_design as ld
from deliverome_analysis import topology
from deliverome_analysis.s3 import external

HERE = Path(__file__).parent
MANIFEST = HERE / "candidate_manifest_tranches.tsv"
TAG_PATHS = HERE / "l1_tag_paths.tsv"
NTERM = HERE / "constructs_l1_nterm.tsv"
INTERNAL = HERE / "constructs_l1_internal.tsv"
LANE_OVERRIDES = HERE / "lane_overrides.tsv"

OUT_MEMBERS = HERE / "topology_vs_experiment.tsv"
OUT_SITES = HERE / "topology_vs_experiment_sites.tsv"
OUT_SUMMARY = HERE / "topology_vs_experiment_summary.tsv"
OUT_INPUTS = HERE / "topology_vs_experiment_inputs.json"

UNIPROT_RELEASE = "2026_04"
UNIPROT_KEY = (
    f"uniprot/{UNIPROT_RELEASE}/rest-2026-11-20/topology-vs-experiment/"
    "reviewed_human_mouse_topology.tsv.gz"
)
UNIPROT_URL = (
    "https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=tsv"
    "&fields=accession,id,gene_primary,organism_id,length,sequence_version,sequence,"
    "ft_signal,ft_topo_dom,ft_transmem,ft_intramem"
    "&query=reviewed%3Atrue%20AND%20(organism_id%3A9606%20OR%20organism_id%3A10090)"
    "%20AND%20(ft_topo_dom%3A*%20OR%20ft_transmem%3A*%20OR%20ft_intramem%3A*)"
)
TOPDB_VERSION = "v3.3"
TOPDB_KEY = f"topdb/{TOPDB_VERSION}/download-2026-11-20/topdb_all.xml.gz"
TOPDB_URL = f"https://topdb.unitmp.org/data/TOPDB/data/{TOPDB_VERSION}/sets/topdb_all.xml"

#: ECO codes that make a UniProt topology feature experimental.
EXPERIMENTAL_ECO = {"ECO:0000269": "experimental", "ECO:0007744": "structure"}
#: TOPO_DOM notes, by side. "out" is the exoplasmic face (extracellular or the lumen of
#: a secretory-pathway compartment), "in" the cytosolic face. Mitochondrial and
#: peroxisomal compartments are neither and are not used.
OUT_NOTES = (
    "Extracellular",
    "Lumenal",
    "Vesicular",
    "Exoplasmic loop",
    "Intragranular",
    "Vacuolar",
    "Perinuclear space",
)
IN_NOTES = ("Cytoplasmic", "Nuclear", "Mother cell cytoplasmic")
TOPDB_LOC = {
    "Outside": "out",
    "Not Inside": "out",
    "Inside": "in",
    "Not Outside": "in",
    "Membrane": "membrane",
}
TOPDB_EXCLUDED_TYPES = {"Other"}
MIN_IDENTITY = 0.6
#: A helix counts as lost/gained when fewer than half its residues align.
HELIX_ALIGNED_FRACTION = 0.5

PLACEMENTS = [
    "N after SP",
    "N at residue 1, multipass",
    "N at residue 1, single-pass (type III)",
    "C-terminal",
    "internal",
]


# --------------------------------------------------------------------------- evidence


@dataclass
class Region:
    start: int  # 1-based, inclusive, in the evidence entry's numbering
    end: int
    side: str  # out / in / membrane
    kind: str  # experimental / structure / topdb:<Type>/<Subtype>
    cite: str


@dataclass
class Entry:
    source: str  # uniprot / topdb
    accession: str
    mnemonic: str
    gene: str
    sequence: str
    regions: list[Region] = field(default_factory=list)  # experimental only
    helices: list[tuple[int, int]] = field(default_factory=list)  # any-evidence TMs
    signal_end: int = 0  # UniProt SIGNAL end (any evidence), 0 if none


_FEATURE = re.compile(
    r"(SIGNAL|TOPO_DOM|TRANSMEM|INTRAMEM) ([<?]?\d+|\?)(?:\.\.([>?]?\d+|\?))?((?:; /\w+=\"[^\"]*\")*)"
)


def _pos(tok: str | None) -> int | None:
    if tok is None:
        return None
    tok = tok.lstrip("<>")
    return int(tok) if tok.isdigit() else None


def _quals(s: str) -> dict[str, str]:
    return dict(re.findall(r'/(\w+)="([^"]*)"', s))


def _experimental(evidence: str) -> tuple[str, str] | None:
    """(kind, citation) if any evidence code is experimental, else None."""
    kinds, cites = [], []
    for ev in [e.strip() for e in evidence.split(",") if e.strip()]:
        code, _, ref = ev.partition("|")
        if code in EXPERIMENTAL_ECO:
            kinds.append(EXPERIMENTAL_ECO[code])
            if ref:
                cites.append(ref)
    if not kinds:
        return None
    kind = "experimental" if "experimental" in kinds else "structure"
    return kind, ";".join(dict.fromkeys(cites)) or code


def _note_side(note: str) -> str | None:
    if note.startswith(OUT_NOTES):
        return "out"
    if note.startswith(IN_NOTES):
        return "in"
    return None


def load_uniprot(path: Path) -> list[Entry]:
    df = pd.read_csv(path, sep="\t", dtype=str, compression="gzip").fillna("")
    entries = []
    for r in df.to_dict("records"):
        e = Entry(
            "uniprot",
            r["Entry"],
            r["Entry Name"],
            r["Gene Names (primary)"].split(";")[0].strip(),
            r["Sequence"],
        )
        blob = " ".join(
            r[c] for c in ("Signal peptide", "Topological domain", "Transmembrane", "Intramembrane")
        )
        for ftype, a, b, qual in _FEATURE.findall(blob):
            start, end = _pos(a), _pos(b if b else a)
            if start is None or end is None:
                continue
            q = _quals(qual)
            if ftype == "SIGNAL":
                e.signal_end = max(e.signal_end, end)
                continue
            if ftype in ("TRANSMEM", "INTRAMEM"):
                e.helices.append((start, end))
            exp = _experimental(q.get("evidence", ""))
            if exp is None:
                continue
            if ftype == "TOPO_DOM":
                side = _note_side(q.get("note", ""))
                if side is None:
                    continue
            else:
                side = "membrane"
            e.regions.append(Region(start, end, side, exp[0], exp[1]))
        if e.regions:
            entries.append(e)
    return entries


def load_topdb(path: Path, gene_of_mnemonic: dict[str, str]) -> list[Entry]:
    ns = "{https://topdb.unitmp.org}"
    entries = []
    with gzip.open(path, "rb") as fh:
        for _, el in ET.iterparse(fh):
            if el.tag != ns + "TOPDB":
                continue
            mnem = el.get("ID", "")
            if not mnem.endswith(("_HUMAN", "_MOUSE")):
                el.clear()
                continue
            seq = re.sub(r"\s+", "", el.findtext(f"{ns}Sequence/{ns}Seq") or "")
            e = Entry("topdb", mnem, mnem, gene_of_mnemonic.get(mnem, ""), seq)
            for reg in el.iterfind(f"{ns}Topology/{ns}Regions/{ns}Region"):
                if reg.get("Loc") == "Membrane":
                    e.helices.append((int(reg.get("Begin")), int(reg.get("End"))))
            # Structure regions are kept per (PDB entry, chain), and only for a chain that
            # TOPDB also places in the membrane. TOPDB carries soluble-fragment structures
            # too (VEGFR1's kinase domain, 3hng, is labelled "Outside"), and a chain with
            # no membrane-embedded segment has no membrane to have a side of.
            chains: dict[tuple[str, str], list[Region]] = collections.defaultdict(list)
            for reg in el.iterfind(f"{ns}Experiments/{ns}Region"):
                side = TOPDB_LOC.get(reg.get("Loc", ""))
                if side is None:
                    continue
                for exp in reg.iterfind(f"{ns}Exp"):
                    etype = exp.findtext(f"{ns}Type") or ""
                    if etype in TOPDB_EXCLUDED_TYPES:
                        continue
                    sub = exp.findtext(f"{ns}Subtype") or ""
                    ref = exp.get("Ref", "0")
                    pdb = exp.find(f"{ns}pdbRef")
                    cite = (
                        f"PDB:{pdb.get('ID')}"
                        if pdb is not None
                        else (f"PubMed:{ref}" if ref not in ("", "0") else "")
                    )
                    region = Region(
                        int(reg.get("Begin")),
                        int(reg.get("End")),
                        side,
                        f"topdb:{etype}/{sub}",
                        cite,
                    )
                    if pdb is not None:
                        chain = ",".join(
                            c.get("ID", "") for c in pdb.iterfind(f"{ns}Chains/{ns}Chain")
                        )
                        chains[(pdb.get("ID"), chain)].append(region)
                    else:
                        e.regions.append(region)
            for regs in chains.values():
                if any(r.side == "membrane" for r in regs):
                    e.regions.extend(regs)
            if e.regions and seq:
                entries.append(e)
            el.clear()
    return entries


# --------------------------------------------------------------------------- mapping

_ALIGNER = Align.PairwiseAligner()
_ALIGNER.mode = "global"
_ALIGNER.substitution_matrix = substitution_matrices.load("BLOSUM62")
_ALIGNER.open_gap_score = -10
_ALIGNER.extend_gap_score = -0.5
_ALIGNER.end_gap_score = 0


def residue_map(member: str, entry: str) -> tuple[dict[int, int], float]:
    """{entry position: member position} (1-based) over aligned columns, and identity."""
    if member == entry:
        return {i: i for i in range(1, len(member) + 1)}, 1.0
    aln = _ALIGNER.align(member, entry)[0]
    mapping, same, cols = {}, 0, 0
    for (ms, me), (es, _ee) in zip(*aln.aligned, strict=True):
        for k in range(me - ms):
            mapping[es + k + 1] = ms + k + 1
            cols += 1
            same += member[ms + k] == entry[es + k]
    return mapping, (same / cols if cols else 0.0)


def helix_runs(topo: str) -> list[tuple[int, int]]:
    return [(m.start() + 1, m.end()) for m in re.finditer(r"M+", topo)]


@dataclass
class Mapped:
    entry: Entry
    match: str  # exact / mapped
    identity: float
    regions: list[Region]  # in member numbering
    lost_helix_ends: list[int]  # member position after which an entry helix is missing
    gained_helix_starts: list[int]  # member helices with no entry counterpart
    helices: list[tuple[int, int]]  # entry helices in member numbering


def map_entry(member_seq: str, topo: str, e: Entry) -> Mapped | None:
    mapping, ident = residue_map(member_seq, e.sequence)
    exact = member_seq == e.sequence
    if not exact and ident < MIN_IDENTITY:
        return None
    regions = []
    for r in e.regions:
        pos = sorted(mapping[p] for p in range(r.start, r.end + 1) if p in mapping)
        if not pos:
            continue
        run = [pos[0]]
        for p in [*pos[1:], None]:
            if p is not None and p == run[-1] + 1:
                run.append(p)
                continue
            regions.append(Region(run[0], run[-1], r.side, r.kind, r.cite))
            if p is not None:
                run = [p]
    helices, lost = [], []
    for s, t in e.helices:
        pos = [mapping[p] for p in range(s, t + 1) if p in mapping]
        if len(pos) >= HELIX_ALIGNED_FRACTION * (t - s + 1):
            helices.append((min(pos), max(pos)))
        else:
            # where in the member the loss sits: the last aligned residue before it
            before = [mapping[p] for p in range(1, s) if p in mapping]
            lost.append(max(before) if before else 0)
    inv = set(mapping.values())
    gained = []
    for s, t in helix_runs(topo):
        aligned = sum(p in inv for p in range(s, t + 1))
        covered = any(not (t < hs or s > he) for hs, he in helices)
        if aligned < HELIX_ALIGNED_FRACTION * (t - s + 1) or (not covered and not exact):
            gained.append(s)
    if exact:
        lost, gained = [], []
    return Mapped(e, "exact" if exact else "mapped", ident, regions, lost, gained, helices)


# --------------------------------------------------------------------------- calls


def segment(topo: str, site: int, bounds: set[int]) -> tuple[int, int] | None:
    """Membrane-free segment around ``site``; None if the site is in a bound."""
    n = len(topo)
    if site in bounds:
        return None
    lo = site
    while lo - 1 >= 1 and (lo - 1) not in bounds:
        lo -= 1
    hi = site
    while hi + 1 <= n and (hi + 1) not in bounds:
        hi += 1
    return lo, hi


@dataclass
class Call:
    side: str  # out / in / membrane / helix_end / conflict / none / context_changed
    sources: list[str]
    cites: list[str]
    detail: str
    minority: str = ""


#: Share of independent evidence items that must agree for a side to be called.
MAJORITY = 2 / 3
#: A site this close to the end of an experimental helix is a boundary call, not an
#: inversion: the cassette sits at the helix/loop junction and the loop is what we test.
HELIX_END_RESIDUES = 3


def _bounds(topo: str, evidence: list[Mapped] | None) -> set[int]:
    b = {i + 1 for i, c in enumerate(topo) if c in "MS"}
    for m in evidence or []:
        for s, t in m.helices:
            b.update(range(s, t + 1))
    return b


def _vote(topo, site, seg, evidence, *, use_membrane=True):
    """{item: (side, region_kinds, ends)} over independent evidence items.

    An item is one citation (a PDB entry, a PubMed record) within one source, so 40 chains
    of one cryo-EM entry count once and a single mis-placed structure cannot outvote the
    rest of the literature.
    """
    items: dict[tuple[str, str], dict] = {}
    skipped = []
    for m in evidence:
        if m.match == "mapped" and (
            any(p < site for p in m.lost_helix_ends) or any(p < site for p in m.gained_helix_starts)
        ):
            skipped.append(m.entry.accession)
            continue
        for r in m.regions:
            if r.side == "membrane":
                if not use_membrane or not (r.start <= site <= r.end):
                    continue
            elif seg is None or r.end < seg[0] or r.start > seg[1]:
                continue
            for cite in (r.cite or f"{m.entry.accession}:{r.kind}").split(";"):
                it = items.setdefault(
                    (m.entry.source, cite), {"sides": set(), "src": set(), "edge": 99}
                )
                it["sides"].add(r.side)
                it["src"].add(f"{m.entry.source}:{m.entry.accession}({m.match};{r.kind})")
                if r.side == "membrane":
                    it["edge"] = min(it["edge"], site - r.start, r.end - site)
    calls = {}
    for key, it in items.items():
        sides = it["sides"]
        side = (
            "membrane" if "membrane" in sides else (next(iter(sides)) if len(sides) == 1 else None)
        )
        if side is not None:
            calls[key] = (side, it["src"], it["edge"])
    return calls, skipped


def call_site(topo: str, site: int, evidence: list[Mapped]) -> Call:
    seg = segment(topo, site, _bounds(topo, evidence))
    calls, skipped = _vote(topo, site, seg, evidence)
    seg_txt = f"segment {seg[0]}-{seg[1]}" if seg else "site in a helix"
    if not calls:
        side = "context_changed" if skipped else "none"
        det = seg_txt + (f"; mapped context changed: {','.join(skipped)}" if skipped else "")
        return Call(side, [], [], det)
    tally = collections.Counter(c[0] for c in calls.values())
    side, n = tally.most_common(1)[0]
    if n < MAJORITY * len(calls):
        side = "conflict"
    minority = ", ".join(
        f"{k[1]}={v[0]}" for k, v in sorted(calls.items()) if v[0] != side and side != "conflict"
    )
    detail = f"{seg_txt}; items " + " ".join(f"{s}:{c}" for s, c in tally.most_common())
    if side == "membrane":
        edge = min(v[2] for v in calls.values() if v[0] == "membrane")
        detail += f"; {edge} residue(s) from the nearest experimental helix end"
        # what does experiment say about the rest of the loop v1 put the site in?
        loop = segment(topo, site, _bounds(topo, None))
        lcalls, _ = _vote(topo, site, loop, evidence, use_membrane=False)
        lt = collections.Counter(c[0] for c in lcalls.values())
        loop_side = "none"
        if lt:
            top, k = lt.most_common(1)[0]
            loop_side = top if k >= MAJORITY * len(lcalls) else "conflict"
        detail += f"; rest of the v1 loop: {loop_side}"
        if edge <= HELIX_END_RESIDUES and loop_side == "out":
            side = "helix_end"
    used = [k for k, v in calls.items() if side in ("conflict", "helix_end") or v[0] == side]
    return Call(
        side,
        sorted({s for k in used for s in calls[k][1]}),
        sorted({k[1] for k in used if ":" in k[1]}),
        detail,
        minority,
    )


def verdict(side: str) -> str:
    return {
        "out": "concordant",
        "helix_end": "concordant",
        "in": "discordant",
        "membrane": "discordant",
        "conflict": "conflict",
    }.get(side, "no_evidence")


def dtm_side(ch: str) -> str:
    return {"O": "out", "I": "in", "M": "membrane", "S": "signal", "B": "membrane"}.get(ch, ch)


# --------------------------------------------------------------------------- fetch


def fetch() -> None:
    """Snapshot both evidence sources to their dated external:// keys (append-only)."""
    for url, key in ((UNIPROT_URL, UNIPROT_KEY), (TOPDB_URL, TOPDB_KEY)):
        if external.exists(key):
            print(f"exists, not overwriting: external://{key}")
            continue
        with tempfile.TemporaryDirectory() as tmp:
            raw = Path(tmp) / "raw"
            with urllib.request.urlopen(url, timeout=600) as resp, open(raw, "wb") as fh:
                headers = dict(resp.headers.items())
                shutil.copyfileobj(resp, fh)
            data = Path(tmp) / "data.gz"
            if key.endswith(".xml.gz"):
                with open(raw, "rb") as src, gzip.open(data, "wb", compresslevel=9) as dst:
                    shutil.copyfileobj(src, dst)
            else:
                shutil.move(raw, data)
            hdr = Path(tmp) / "headers.json"
            hdr.write_text(json.dumps({"url": url, "headers": headers}, indent=2))
            external.upload(data, key)
            external.upload(hdr, key.rsplit(".", 2)[0] + ".headers.json")
            print(f"uploaded external://{key} (release {headers.get('X-UniProt-Release', '')})")


# --------------------------------------------------------------------------- main


def _wilson(k: int, n: int, z: float = 1.96) -> str:
    """Wilson score interval for k of n, in percent."""
    if n == 0:
        return ""
    p = k / n
    centre = (p + z * z / (2 * n)) / (1 + z * z / n)
    half = z * ((p * (1 - p) / n + z * z / (4 * n * n)) ** 0.5) / (1 + z * z / n)
    return f"{100 * max(0.0, centre - half):.1f}-{100 * min(1.0, centre + half):.1f}"


def _tsv(path: Path) -> list[dict[str, str]]:
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--fetch", action="store_true", help="snapshot the evidence sources")
    args = ap.parse_args(argv)
    if args.fetch:
        fetch()
        return 0

    uni = load_uniprot(external.download(UNIPROT_KEY))
    up_all = pd.read_csv(
        external.download(UNIPROT_KEY), sep="\t", dtype=str, compression="gzip"
    ).fillna("")
    gene_of_mnem = dict(
        zip(
            up_all["Entry Name"],
            up_all["Gene Names (primary)"].str.split(";").str[0].str.strip(),
            strict=True,
        )
    )
    tdb = load_topdb(external.download(TOPDB_KEY), gene_of_mnem)
    entries = uni + tdb
    by_seq: dict[str, list[Entry]] = collections.defaultdict(list)
    by_gene: dict[str, list[Entry]] = collections.defaultdict(list)
    for e in entries:
        by_seq[e.sequence].append(e)
        if e.gene:
            by_gene[e.gene.upper()].append(e)
    print(f"evidence entries: UniProt {len(uni)}, TOPDB {len(tdb)}")

    paths = {r["candidate_id"]: r for r in _tsv(TAG_PATHS)}
    members = ld.designed_proteins(MANIFEST, library=ld.LIBRARY_L1)
    members = {cid: m for cid, m in members.items() if cid in paths}
    report = ld.topology_report(members, tool=topology.TOOL_V1)
    topo_of = {
        r["candidate_id"]: r["per_residue_topology"]
        for r in report.to_dict("records")
        if r["found"]
    }
    design_id = ld.design_ids(MANIFEST)
    nterm_site = {r["candidate_id"]: int(r["native_after_residue"]) + 1 for r in _tsv(NTERM)}
    arms: dict[str, list[dict]] = collections.defaultdict(list)
    for r in _tsv(INTERNAL):
        if r["decision"].startswith("internal_loop") and r["cassette_after_residue"]:
            arms[r["candidate_id"]].append(r)
    overrides = {r["design_id"]: r for r in _tsv(LANE_OVERRIDES)}

    member_rows, site_rows = [], []
    for cid, p in paths.items():
        m = members[cid]
        seq, topo = m["sequence"], topo_of[cid]
        assert len(seq) == len(topo), cid
        assert hashlib.sha256(seq.encode()).hexdigest() == p["sequence_sha256"], cid
        cands = list(by_seq.get(seq, []))
        exact_ids = {id(e) for e in cands}
        if not cands:
            cands = [e for e in by_gene.get(m["gene_symbol"].upper(), []) if id(e) not in exact_ids]
        evidence = [x for x in (map_entry(seq, topo, e) for e in cands) if x is not None]
        exact = any(x.match == "exact" for x in evidence)
        if exact:
            evidence = [x for x in evidence if x.match == "exact"]

        placement = p["placement"]
        sp_end = len(topo) - len(topo.lstrip("S"))
        if placement == "N after SP":
            sites = [("N after SP", nterm_site.get(cid, int(p["signal_peptide_aa"] or 0) + 1))]
        elif placement.startswith("N at residue 1"):
            sites = [("residue 1", 1)]
        elif placement == "C-terminal":
            sites = [("C-terminus", len(seq))]
        else:
            sites = [
                (
                    f"{a['arm_role']} after {a['cassette_after_residue']}",
                    int(a["cassette_after_residue"]),
                )
                for a in arms.get(cid, [])
            ]

        site_verdicts, site_sides = [], []
        for label, site in sites:
            c = call_site(topo, site, evidence)
            v = verdict(c.side)
            site_verdicts.append(v)
            site_sides.append(c.side)
            site_rows.append(
                {
                    "candidate_id": cid,
                    "design_id": design_id.get(cid, cid),
                    "gene_symbol": m["gene_symbol"],
                    "placement": placement,
                    "site": label,
                    "site_residue": site,
                    "deeptmhmm_v1_site": dtm_side(topo[site - 1]),
                    "experimental_site_side": c.side,
                    "verdict": v,
                    "segment": c.detail,
                    "minority_evidence": c.minority,
                    "evidence": "; ".join(c.sources),
                    "citations": "; ".join(c.cites),
                }
            )
        if not sites:
            member_verdict = "no_site"
        elif "discordant" in site_verdicts:
            member_verdict = "discordant"
        elif "concordant" in site_verdicts:
            member_verdict = "concordant"
        elif "conflict" in site_verdicts:
            member_verdict = "conflict"
        else:
            member_verdict = "no_evidence"

        n_call = call_site(topo, sp_end + 1, evidence)
        c_call = call_site(topo, len(seq), evidence)
        uni_sp = [
            x
            for x in evidence
            if x.entry.source == "uniprot" and x.match == "exact" and x.entry.signal_end
        ]
        notes = []
        if placement.startswith("N at residue 1") and uni_sp:
            notes.append(
                f"UniProt SIGNAL 1-{uni_sp[0].entry.signal_end}: residue-1 tag may be cleaved"
            )
        if (
            placement == "N after SP"
            and evidence
            and not any(x.entry.signal_end for x in evidence if x.entry.source == "uniprot")
        ):
            if any(x.entry.source == "uniprot" for x in evidence):
                notes.append("UniProt has no SIGNAL feature")
        member_rows.append(
            {
                "candidate_id": cid,
                "design_id": design_id.get(cid, cid),
                "gene_symbol": m["gene_symbol"],
                "candidate_set": p["candidate_set"],
                "species": p["species"],
                "placement": placement,
                "rung_before_gpi": p["rung_before_gpi"],
                "gpi_rerouted": p["gpi_rerouted"],
                "lane_override": "yes" if design_id.get(cid, cid) in overrides else "",
                "n_tm_v1": p["n_tm"],
                "signal_peptide_aa_v1": sp_end,
                "evidence_match": "exact" if exact else ("mapped" if evidence else "none"),
                "evidence_entries": "; ".join(
                    f"{x.entry.source}:{x.entry.accession}({x.match},{x.identity:.2f})"
                    for x in evidence
                ),
                "n_tagged_sites": len(sites),
                "verdict": member_verdict,
                "discordance": (
                    "wrong side"
                    if "in" in site_sides
                    else ("site inside a helix" if "membrane" in site_sides else "")
                ),
                "site_at_helix_end": "yes" if "helix_end" in site_sides else "",
                "deeptmhmm_v1_n_side": dtm_side(topo[sp_end]) if sp_end < len(topo) else "",
                "experimental_n_side": n_call.side,
                "deeptmhmm_v1_c_side": dtm_side(topo[-1]),
                "experimental_c_side": c_call.side,
                "n_citations": "; ".join(n_call.cites),
                "c_citations": "; ".join(c_call.cites),
                "notes": "; ".join(notes),
                "sequence_sha256": p["sequence_sha256"],
            }
        )

    mem = pd.DataFrame(member_rows)
    sites_df = pd.DataFrame(site_rows)
    mem.to_csv(OUT_MEMBERS, sep="\t", index=False)
    sites_df.to_csv(OUT_SITES, sep="\t", index=False)

    rows = []
    for pl in [*PLACEMENTS, "all"]:
        d = mem if pl == "all" else mem[mem.placement == pl]
        ev = d[d.verdict.isin(["concordant", "discordant"])]
        row = {
            "placement": pl,
            "members": len(d),
            "with_experimental_call": len(ev),
            "via_exact_sequence": int((ev.evidence_match == "exact").sum()),
            "concordant": int((ev.verdict == "concordant").sum()),
            "discordant": int((ev.verdict == "discordant").sum()),
            "discordant_wrong_side": int((ev.discordance == "wrong side").sum()),
            "discordant_site_inside_helix": int((ev.discordance == "site inside a helix").sum()),
            "discordant_pct": round(100 * (ev.verdict == "discordant").mean(), 1)
            if len(ev)
            else "",
            "discordant_pct_95ci": _wilson(int((ev.verdict == "discordant").sum()), len(ev)),
            "concordant_but_site_at_helix_end": int(
                ((ev.verdict == "concordant") & (ev.site_at_helix_end == "yes")).sum()
            ),
            "conflicting_evidence": int((d.verdict == "conflict").sum()),
            "evidence_entry_but_no_call_at_site": int(
                ((d.verdict == "no_evidence") & (d.evidence_match != "none")).sum()
            ),
            "no_evidence_entry": int((d.evidence_match == "none").sum()),
            "no_tagged_site": int((d.verdict == "no_site").sum()),
        }
        # the whole v1 call, not just the tagged site: does v1 put each terminus on the
        # side experiment does?
        for t in ("n", "c"):
            e = d[d[f"experimental_{t}_side"].isin(["out", "in"])]
            row[f"{t}_terminus_with_call"] = len(e)
            row[f"{t}_terminus_v1_disagrees"] = int(
                (e[f"experimental_{t}_side"] != e[f"deeptmhmm_v1_{t}_side"]).sum()
            )
        rows.append(row)
    summary = pd.DataFrame(rows)
    summary.to_csv(OUT_SUMMARY, sep="\t", index=False)
    OUT_INPUTS.write_text(
        json.dumps(
            {
                "uniprot": {
                    "key": f"external://{UNIPROT_KEY}",
                    "release": UNIPROT_RELEASE,
                    "etag": external.etag(UNIPROT_KEY),
                },
                "topdb": {
                    "key": f"external://{TOPDB_KEY}",
                    "version": TOPDB_VERSION,
                    "etag": external.etag(TOPDB_KEY),
                },
                "topology_tool": topology.TOOL_V1,
                "experimental_eco": sorted(EXPERIMENTAL_ECO),
                "topdb_excluded_types": sorted(TOPDB_EXCLUDED_TYPES),
                "min_identity_mapped": MIN_IDENTITY,
                "members": len(mem),
            },
            indent=2,
        )
        + "\n"
    )
    print(summary.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
