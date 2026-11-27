"""Check, member by member, that a C-terminal tag survives on the L1 C-terminal rung.

``l1_tag_paths.py`` puts 218 L1 members on the C-terminal rung because their N-terminus is
cytoplasmic and their C-terminus extracellular by DeepTMHMM v1. That is a topology call. It
does not say whether the C-terminus is still there on the cell surface, or whether a 51-aa
extension is a safe thing to put on it. This script asks those questions for each member and
writes one row per member to ``l1_cterm_review.tsv``. ``build_l1_cterm_arms.py`` builds every
member whose ``build`` is ``build`` and holds the rest.

Every flag is computed from a stated source; none is a judgement typed in by hand.

``cterm_peptide_removed`` (holds the member)
    A UniProt Peptide or Propeptide feature covers the member's last residue, so the
    C-terminus is cut off in the secretory pathway and the tag leaves with it. This is the
    GPI rule (an omega site removes the C-terminus) applied to proteolytic processing.
``gpi`` (holds the member)
    The gene is in ``gpi_set.tsv``. ``l1_tag_paths.py`` already reroutes these, so a C-terminal
    member in the GPI set means the two tables disagree.
``ectodomain_released``
    UniProt records a cleavage that releases the extracellular domain (a sheddase, furin, a
    soluble or secreted form). On a type II protein the C-terminus is on that domain, so the
    tag is released with it. Any extracellular tag on these members has the same exposure;
    there is no other rung to move them to. Recorded alongside the existing ``shed_risk``.
``two_chain``
    UniProt records more than one chain, named as heavy/light or (non-)catalytic chains: an
    activation cleavage splits the ectodomain (the serine proteases, the GGTs). The tag is on
    the C-terminal chain, which stays on the cell only while it stays associated with the
    anchored chain (by a disulfide in the serine proteases); on a gel it can run apart.
``partial_form``
    The member is under half the length of its gene's UniProt canonical protein.
``short_ecto_tail``
    The v1 C-terminal extracellular run is under 5 residues, so the cassette starts at the
    membrane.
``lane_override_topology``
    The member is C-terminal by a curated row in ``lane_overrides.tsv`` (experimental topology
    that reverses v1), so v1's tail length and stalk are not computed for it.
``v2_cterm_cytoplasmic`` / ``v2_no_tm`` / ``no_v2``
    DeepTMHMM v2, on the same sequence, puts the C-terminus in the cytoplasm, calls no TM
    helix at all, or has no row. v1 is the topology of record; these mark weak calls.
``structured_terminus``
    The member's last residue has pLDDT >= 85 in its own structure model, so the terminus
    is inside a confidently predicted fold rather than a free tail.
``cterm_differs_from_uniprot``
    The member's last 15 residues are not in its gene's UniProt canonical sequence, so no
    UniProt feature can be transferred to its C-terminus.
``uniprot_not_checked``
    The member's sequence has no UniProt entry (cyno orthologs from Ensembl), so the
    processing checks above (``cterm_peptide_removed``, ``two_chain``, ``ectodomain_released``)
    were not run on it. The evidence names the human member of the same gene and its flags, as
    context, not as a transferred call. If that human member's C-terminal peptide is removed,
    the review stops: such a member cannot be kept without its own check.

Two measurements are recorded for every member with a model, as evidence for the
structured-terminus calls (``l1_cterm_decisions.py``):

``terminal_ca_neighbors_10a``
    Other CA atoms within 10 A of the last residue's CA (excluding 3 neighbors either side).
    A free tail scores about 2; the M13 metalloproteases' last Trp scores 16-20.
``stalk_residues`` / ``stalk_site_after``
    The stretch from the residue after the last TM helix to the residue before the first run
    of 10 residues at pLDDT >= 70, and its midpoint when it is at least 8 residues long. That
    midpoint is the internal site a member would move to if its C-terminal refold fails.

Sources: the pinned v1/v2 topology snapshot (``deliverome_analysis.surfaceome.topology``), UniProt
entries (REST, saved once to S3 and read back with ``--offline``), the member's structure
model from ``library_structure_payload.tsv`` (the archived AlphaFold DB v6 copy at
``external://alphafold-db/v6/`` or our ColabFold run on S3; never a live download),
``gpi_set.tsv`` and ``l1_tag_paths.tsv``.

    uv run --locked python review_l1_cterm_members.py            # UniProt from the pinned copy
    uv run --locked python review_l1_cterm_members.py --check    # rerun, compare with the record
    uv run --locked python review_l1_cterm_members.py --fetch    # refetch UniProt, then publish

Members are read from ``l1_tag_paths.tsv``, never listed here, so the review and the build
rerun unchanged when the member set changes (for example the Compara 116 ortholog panel).
A changed set needs ``--fetch`` with a new dated ``UNIPROT_PROCESSED_KEY``, its copy to
``external://`` (``publish_upstream_copies.py``) and the new ``UNIPROT_KEY`` / ``UNIPROT_SHA256``
pinned here, and any new C-terminal member needs its cassette encoding issued by ``l1_swap_risk.py`` before ``build_l1_cterm_arms.py``.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
import tempfile
import time
import urllib.parse
import urllib.request
from collections import Counter
from pathlib import Path

from deliverome_analysis.cloning import library_design as ld
from deliverome_analysis.surfaceome import topology

HERE = Path(__file__).parent
MANIFEST = HERE / "candidate_manifest_tranches.tsv"
TAG_PATHS = HERE / "l1_tag_paths.tsv"
GPI_FILE = HERE / "gpi_set.tsv"
PAYLOAD = HERE / "library_structure_payload.tsv"
OUT = HERE / "l1_cterm_review.tsv"
PREFIX = "sprint4/s4e9-orf-design-fg-panel/"
# 2026-11-10 snapshot: the 218 members before PR #283's ortholog adoption and isoform redraw.
# 2026-11-12 snapshot: refetched for the member set at #283's 9928f4b1 (new members joined).
# 2026-11-12-second-pass: refetched for #283's second pass (6bf341f2), which made one more
# member C-terminal.
# 2026-11-13-slc68a1: the second-pass snapshot plus a fetch for SLC68A1 (l1-p1995) only, which
# l1-reviewer-package moved to the C-terminal lane on experimental topology (lane_overrides.tsv).
# Built with --fetch-new, so the other members keep their 2026-11-12 entries.
# 2026-11-21-coreceptor-additions: the slc68a1 snapshot plus KLRK1 (l1-p3164), the one C-terminal
# member among the co-receptor additions (Becca, 2026-11-21), fetched with --fetch-new.
# 2026-11-22-round2-additions: the coreceptor-additions snapshot plus the C-terminal members of the
# round-2 gene additions and the isoform redraw (2026-11-22), fetched with --fetch-new.
# 2026-11-22-mane-gtf: the round2-additions snapshot with the entries of ECE1 (l1-p0537) and GGT5
# (l1-p0696) fetched again, the two C-terminal members whose reference moved to the GTF MANE Select
# protein on 2026-11-22 (--fetch-new --refetch l1-p0537,l1-p0696).
PREVIOUS_UNIPROT_KEY = (
    PREFIX + "library-design/cterm-review/2026-11-22-mane-gtf/uniprot_entries.json.gz"
)
#: Where ``--fetch`` / ``--fetch-new`` write a new snapshot (processed, the only bucket a fetch
#: may write). It is then copied to the append-only external bucket by
#: ``publish_upstream_copies.py`` (mapping in ``upstream_copies.tsv``).
UNIPROT_PROCESSED_KEY = (
    PREFIX + "library-design/cterm-review/2026-11-23-pre-freeze/uniprot_entries.json.gz"
)
#: The snapshot the review reads: the append-only external copy of UNIPROT_PROCESSED_KEY, pinned
#: by SHA-256 (byte-identical to the processed object, checked 2026-11-14), so nothing can change
#: under the review. Move the key and the hash together after a new fetch has been copied.
UNIPROT_KEY = "uniprot/rest-2026-11-23/cterm-review-pre-freeze/uniprot_entries.json.gz"
UNIPROT_SHA256 = "d83dafa75801faa03f347a50c27f8d94e1173cb6c8e0fc904bf47d32f2828823"
COLABFOLD = PREFIX + "predicted-structures/colabfold/"
COLABFOLD_RUNS = ("2026-10-14-mismatched", "2026-10-14", "2026-11-02")
#: candidate_id (l1-mNNNN, plate order) <-> design_id (l1-pNNNN, permanent). The UniProt
#: snapshot and our ColabFold model files are keyed by design_id, which never changes, so a
#: renumbering does not orphan them.
MAPPING = HERE / "l1_member_renumbering.tsv"


def design_ids() -> dict[str, str]:
    """candidate_id -> design_id. Since the per-insert cut (final-e) the two are the same id, so a
    member added after the last renumbering (which runs at the cut) maps to itself."""
    with MAPPING.open() as fh:
        known = {r["candidate_id"]: r["design_id"] for r in csv.DictReader(fh, delimiter="\t")}
    with MANIFEST.open() as fh:
        current = {r["candidate_id"] for r in csv.DictReader(fh, delimiter="\t")}
    return {**{c: c for c in current}, **known}


TAXON = {"human": 9606, "mouse": 10090, "cyno": 9541}

# The thresholds below set flags only; none changes what is built. They were chosen by the
# Claude Code session claude/l1-cterm-build-16eeeb (PR #284) on 2026-11-10 with no person's
# decision behind them (attribution convention, 2026-11-10): change them freely.
STRUCTURED_PLDDT = 85.0
SHORT_TAIL = 5
PARTIAL_FRACTION = 0.5
CTERM_MATCH = 15
RELEASE = re.compile(r"shed|\bsoluble\b|secreted form|released|furin|ADAM\d+", re.I)
TWO_CHAIN = re.compile(r"heavy chain|light chain|non-catalytic chain|catalytic chain", re.I)

#: DeepTMHMM v2 side characters. For every membrane index except the two outer membranes
#: (16, 17) the upper-case character is the cytoplasmic side; for 16 and 17 the sense is
#: inverted but the characters are X/Y for the cytoplasmic side, so upper case = cytoplasm
#: holds across the whole alphabet. Digits: 3/5/7 cytoplasmic, 4/6/8 not. Decoding from
#: accessible-surfaceome ``sources/deeptmhmm2.py`` MEMBRANE_TOPOLOGY_MAP / CYTOPLASMIC_SIDE.
V2_STATES = set("SMBRF>")


def v2_side(ch: str) -> str:
    if not ch or ch in V2_STATES:
        return ""
    if ch in "357":
        return "cytoplasmic"
    if ch in "468":
        return "non-cytoplasmic"
    if ch in "OI":  # v1-alphabet sides on a v2 'Globular' row
        return {"O": "non-cytoplasmic", "I": "cytoplasmic"}[ch]
    return "cytoplasmic" if ch.isupper() else "non-cytoplasmic"


def trailing(s: str, ch: str) -> int:
    n = 0
    while n < len(s) and s[-1 - n] == ch:
        n += 1
    return n


def _get(url: str) -> dict:
    req = urllib.request.Request(url, headers={"User-Agent": "python-urllib"})
    for attempt in range(6):
        try:
            with urllib.request.urlopen(req, timeout=60) as fh:
                return json.load(fh)
        except Exception:
            if attempt == 5:
                raise
            time.sleep(5 * (attempt + 1))
    raise AssertionError


def fetch_uniprot(members: list[dict]) -> dict[str, dict | None]:
    """Each member's gene's reviewed UniProt entry, by gene symbol and taxon."""
    out: dict[str, dict | None] = {}
    for m in members:
        hit = None
        for reviewed in (" AND reviewed:true", ""):
            q = urllib.parse.quote(
                f"gene_exact:{m['gene_symbol']} AND organism_id:{TAXON[m['species']]}{reviewed}"
            )
            res = _get(f"https://rest.uniprot.org/uniprotkb/search?query={q}&format=json&size=1")
            if res["results"]:
                hit = res["results"][0]
                break
        out[m["candidate_id"]] = hit
        time.sleep(0.4)
    return out


def load_uniprot(
    members: list[dict], fetch: bool, fetch_new: bool = False, refetch: frozenset = frozenset()
) -> dict[str, dict | None]:
    """UniProt entries keyed by candidate_id; the snapshot itself is keyed by design_id.

    ``fetch_new`` copies PREVIOUS_UNIPROT_KEY and fetches only the members it lacks, so a lane
    change adds entries without refreshing everyone else's. ``refetch`` (candidate ids) are
    fetched again even when the previous snapshot has them: members whose protein changed.
    """
    from deliverome_analysis.core.s3 import external, processed

    design = design_ids()
    if fetch_new:
        with gzip.open(processed.download(PREVIOUS_UNIPROT_KEY), "rt") as fh:
            prev = json.load(fh)
        new = [
            m
            for m in members
            if design[m["candidate_id"]] not in prev or m["candidate_id"] in refetch
        ]
        print(f"fetching {len(new)} member(s) absent from the previous snapshot")
        cid_of = {d: c for c, d in design.items()}
        entries = {cid_of[d]: v for d, v in prev.items() if d in cid_of}
        entries.update(fetch_uniprot(new))
        fetch = True
    elif fetch:
        entries = fetch_uniprot(members)
    if fetch:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "uniprot_entries.json.gz"
            snap = {design[c]: v for c, v in entries.items()}
            if fetch_new:
                snap = {**prev, **snap}
            with gzip.open(path, "wt") as fh:
                json.dump(snap, fh, sort_keys=True)
            processed.upload(path, UNIPROT_PROCESSED_KEY)
        return entries
    path = Path(external.download(UNIPROT_KEY))
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != UNIPROT_SHA256:
        raise SystemExit(f"external://{UNIPROT_KEY} changed (sha256 {digest}); re-pin it")
    with gzip.open(path, "rt") as fh:
        snap = json.load(fh)
    return {c: snap[d] for c, d in design.items() if d in snap}


#: AlphaFold DB v6 models are read from their archived copies (l1-tag-refold), keyed by the
#: model's file name, so the review never downloads a live model.
AFDB_ARCHIVE = "alphafold-db/v6/"


#: The ``model`` column of a member whose AlphaFold DB model has no archived copy yet (listed in
#: ``tag_refold_afdb_unarchived.tsv``). Its structure-derived columns are blank.
NOT_ARCHIVED = "native model not archived"


def afdb_archive_key(model_url: str) -> str:
    """The external key of an AlphaFold DB model's archived copy."""
    name = model_url.rsplit("/", 1)[1]
    assert name.startswith("AF-") and name.endswith("_v6.pdb"), model_url
    return AFDB_ARCHIVE + name


def terminal_plddt(cid: str, payload: dict[str, dict], cache: Path) -> tuple[float | None, str]:
    """pLDDT of the member's last residue in its own model, and the model's source."""
    model = member_model(cid, payload, cache)
    if model is None:
        row = payload.get(cid)
        if row is None:
            return None, "no model"
        if row["source"] == "colabfold":
            return None, "colabfold model not on S3"
        return None, NOT_ARCHIVED
    plddt, _, name = model
    return plddt[-1], name


def member_model(
    cid: str, payload: dict[str, dict], cache: Path
) -> tuple[list[float | None], list[tuple | None], str] | None:
    """Per-residue pLDDT and CA coordinates in the member's own numbering, and the source."""
    from deliverome_analysis.core.s3 import external, processed

    row = payload.get(cid)
    if row is None:
        return None
    design = design_ids()[cid]
    path = cache / f"{design}.pdb"
    if not path.exists():
        if row["source"] == "colabfold":
            keys = [f"{COLABFOLD}{run}/{design}.pdb" for run in COLABFOLD_RUNS]
            key = next((k for k in keys if processed.exists(k)), None)
            if key is None:
                return None
            processed.download(key, path)
        else:
            # The saved copy, never the live AlphaFold DB: this stage runs offline.
            key = afdb_archive_key(row["model_url"])
            if not external.exists(key):
                # Recorded per member and carried on (the pLDDT, burial and stalk columns stay
                # blank); archive it (tag_refold/archive_afdb_models.py) and rerun to fill them.
                return None
            external.download(key, path)
    atoms = {
        int(line[22:26]): (
            float(line[60:66]),
            tuple(float(line[30 + 8 * i : 38 + 8 * i]) for i in range(3)),
        )
        for line in path.read_text().splitlines()
        if line.startswith("ATOM") and line[12:16].strip() == "CA"
    }
    off, n = int(row["residue_offset"] or 0), int(row["ordered_length_aa"])
    got = [atoms.get(i + off) for i in range(1, n + 1)]
    return (
        [a[0] if a else None for a in got],
        [a[1] if a else None for a in got],
        f"{row['source']} ({row['proteoform_match']})",
    )


#: A folded segment starts at the first run of this many consecutive residues at or above
#: FOLDED_PLDDT; everything between the last TM helix and it is the stalk.
FOLDED_RUN, FOLDED_PLDDT = 10, 70.0
#: A stalk shorter than this cannot take the cassette with a few residues of clearance on
#: each side (membrane and folded domain).
MIN_STALK, STALK_CLEARANCE = 8, 3


def terminal_burial(coords: list[tuple | None]) -> int | None:
    """Other CA atoms within 10 A of the last residue's CA, excluding its 3 neighbors each side."""
    import math

    last = coords[-1]
    if last is None:
        return None
    n = len(coords)
    return sum(
        1
        for i, c in enumerate(coords)
        if c is not None and abs(i - (n - 1)) > 3 and math.dist(c, last) < 10
    )


def stalk(topo: str, plddt: list[float | None]) -> tuple[int, int, int | None]:
    """(first, last) stalk residue (1-based; last < first when there is none) and a site.

    The stalk runs from the residue after the last TM helix to the residue before the first
    folded run. The proposed site is its midpoint, if the stalk leaves STALK_CLEARANCE
    residues on each side; otherwise None.
    """
    first = topo.rfind("M") + 2  # 1-based residue after the last TM residue
    run = 0
    end = len(topo)
    for i in range(first - 1, len(topo)):
        v = plddt[i]
        run = run + 1 if v is not None and v >= FOLDED_PLDDT else 0
        if run == FOLDED_RUN:
            end = i - FOLDED_RUN + 1  # 1-based last stalk residue
            break
    length = end - first + 1
    if length < MIN_STALK:
        return first, end, None
    return first, end, first + length // 2 - 1


#: Flags that come from UniProt processing features (not run without a UniProt entry).
UNIPROT_PROCESSING_FLAGS = ("cterm_peptide_removed", "two_chain", "ectodomain_released")


#: UniProt features that place a residue relative to the membrane.
TOPOLOGY_FEATURES = ("Topological domain", "Transmembrane", "Intramembrane")


def uniprot_cterm_topology(protein: str, entry: dict | None) -> str:
    """UniProt's membrane side at the member's C-terminus, as an independent source.

    Where the member's last residues are in the gene's UniProt sequence, this is the feature
    covering its last residue, with its evidence codes (ECO:0000269 is experimental,
    ECO:0000255 a sequence-model prediction). Where they are not, no feature can be transferred
    to the terminus itself, so the text gives the last residue the member shares with the
    UniProt sequence from the N-terminus, the feature there and the one before it. Blank
    without a UniProt entry.
    """
    if not entry:
        return ""
    seq = entry["sequence"]["value"]

    def at(pos: int) -> dict | None:
        for f in entry.get("features", []):
            s, e = f["location"]["start"].get("value"), f["location"]["end"].get("value")
            if f["type"] in TOPOLOGY_FEATURES and s and e and s <= pos <= e:
                return f
        return None

    def text(f: dict) -> str:
        s, e = f["location"]["start"]["value"], f["location"]["end"]["value"]
        eco = ",".join(sorted({x["evidenceCode"] for x in f.get("evidences", [])}))
        name = f.get("description", "") if f["type"] == "Topological domain" else f["type"]
        return f"{name} {s}-{e}" + (f" ({eco})" if eco else "")

    hit = seq.rfind(protein[-CTERM_MATCH:])
    if hit >= 0:
        f = at(hit + CTERM_MATCH)
        return text(f) if f else "no topology feature at the C-terminus"
    shared = 0
    while shared < min(len(seq), len(protein)) and seq[shared] == protein[shared]:
        shared += 1
    f = at(shared) if shared else None
    if f is None:
        return "C-terminus not in the UniProt sequence; no feature to transfer"
    before = at(f["location"]["start"]["value"] - 1)
    return (
        f"C-terminus not in the UniProt sequence: member 1-{shared} of {len(protein)} is "
        f"UniProt 1-{shared}, which ends in {text(f)}"
        + (f", after {text(before)}" if before else "")
    )


def uniprot_architecture(entry: dict | None) -> dict:
    """UniProt's cleaved signal peptide (as ``start-end``) and its count of TM helices."""
    feats = (entry or {}).get("features", [])
    signal = [
        f"{f['location']['start']['value']}-{f['location']['end']['value']}"
        for f in feats
        if f["type"] == "Signal"
        and f["location"]["start"].get("value")
        and f["location"]["end"].get("value")
    ]
    return {
        "uniprot_signal_peptide": signal[0] if signal else "",
        "uniprot_tm_helices": sum(f["type"] == "Transmembrane" for f in feats) if entry else "",
    }


def uniprot_flags(protein: str, entry: dict | None) -> dict:
    """What UniProt records about the C-terminus this member actually has."""
    if not entry:
        return {"uniprot_acc": "", "matches": "", "flags": [], "evidence": [], "canon_len": 0}
    seq = entry["sequence"]["value"]
    at = seq.rfind(protein[-CTERM_MATCH:])
    end = at + CTERM_MATCH if at >= 0 else None
    flags, evidence = [], []
    feats = entry.get("features", [])
    if end is not None:
        for f in feats:
            s, e = f["location"]["start"].get("value"), f["location"]["end"].get("value")
            if f["type"] in ("Peptide", "Propeptide") and s and e and s <= end <= e:
                flags.append("cterm_peptide_removed")
                evidence.append(f"{f['type']} {s}-{e} {f.get('description', '')}".strip())
        chains = [f for f in feats if f["type"] == "Chain"]
        if len(chains) > 1 and any(TWO_CHAIN.search(f.get("description", "")) for f in chains):
            flags.append("two_chain")
            evidence += [
                f"Chain {f['location']['start'].get('value')}-{f['location']['end'].get('value')} "
                f"{f['description']}"
                for f in feats
                if f["type"] == "Chain" and TWO_CHAIN.search(f.get("description", ""))
            ]
        ptm = [
            t["value"]
            for c in entry.get("comments", [])
            if c["commentType"] == "PTM"
            for t in c.get("texts", [])
        ]
        cuts = [f for f in feats if f["type"] == "Site" and "leav" in f.get("description", "")]
        released = [t for t in ptm if RELEASE.search(t)]
        if released:
            flags.append("ectodomain_released")
            evidence += [t.split(". ")[0][:160] for t in released]
            evidence += [
                f"cleavage {f['location']['start'].get('value')}-"
                f"{f['location']['end'].get('value')} {f['description']}"
                for f in cuts
            ]
    return {
        "uniprot_acc": entry["primaryAccession"],
        "matches": "yes" if end is not None else "no",
        "flags": flags,
        "evidence": evidence,
        "canon_len": len(seq),
    }


def members() -> list[dict]:
    proteins = ld.designed_proteins(MANIFEST, library=ld.LIBRARY_L1)
    with TAG_PATHS.open() as fh:
        rows = [r for r in csv.DictReader(fh, delimiter="\t") if r["placement"] == "C-terminal"]
    for r in rows:
        r["protein"] = proteins[r["candidate_id"]]["sequence"]
    return rows


def _arg_parser() -> argparse.ArgumentParser:
    """The command line ``main`` reads."""
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fetch", action="store_true", help="refetch UniProt and publish it")
    ap.add_argument(
        "--fetch-new",
        action="store_true",
        help="copy the previous snapshot, fetch only members it lacks, and publish",
    )
    ap.add_argument(
        "--refetch",
        default="",
        help="with --fetch-new: candidate ids fetched again although the snapshot has them",
    )
    ap.add_argument("--out", type=Path, default=OUT)
    ap.add_argument("--check", action="store_true", help="rerun and compare with the record")
    return ap


def main(argv=None) -> int:
    args = _arg_parser().parse_args(argv)
    if args.check:
        assert not (args.fetch or args.fetch_new), "--check reads the pinned snapshot"
        args.out = Path(tempfile.mkdtemp()) / OUT.name

    ms = members()
    refetch = frozenset(filter(None, args.refetch.split(",")))
    assert not refetch or args.fetch_new, "--refetch needs --fetch-new"
    uniprot = load_uniprot(ms, args.fetch, args.fetch_new, refetch)
    # Members come from l1_tag_paths.tsv, so a redrawn panel changes them. A member absent
    # from the UniProt snapshot must not read as "no UniProt entry"; refetch instead.
    absent = [m["candidate_id"] for m in ms if m["candidate_id"] not in uniprot]
    if absent:
        raise SystemExit(
            f"{len(absent)} C-terminal member(s) are not in the UniProt snapshot "
            f"({absent[:5]}); the member set changed, so rerun with --fetch to a new dated key"
        )
    with GPI_FILE.open() as fh:
        gpi = {
            r["gene_symbol"] for r in csv.DictReader(fh, delimiter="\t") if r["gpi_anchored"] == "1"
        }
    with PAYLOAD.open() as fh:
        payload = {
            r["candidate_id"]: r
            for r in csv.DictReader(fh, delimiter="\t")
            if r["is_chosen"] == "1"
        }
    cache = Path(tempfile.gettempdir()) / "s4e9-cterm-models-archived"
    cache.mkdir(exist_ok=True)

    # C-terminal by a curated lane override (lane_overrides.tsv: experimental topology that
    # reverses v1), not by v1. v1's sides do not apply to these members, so neither the v1
    # orientation check, the v1 tail length nor the v1-derived stalk is computed for them.
    with (HERE / "lane_overrides.tsv").open() as fh:
        overrides = {
            r["design_id"]: r
            for r in csv.DictReader(fh, delimiter="\t")
            if r["placement"] == "C-terminal"
        }
    out = []
    for m in ms:
        p = m["protein"]
        ov = overrides.get(design_ids()[m["candidate_id"]])
        if ov:
            assert ov["sequence_sha256"] == hashlib.sha256(p.encode()).hexdigest(), (
                f"{m['candidate_id']}: lane override is for a different protein"
            )
        v1 = topology.for_sequence(p, tool=topology.TOOL_V1)
        v2 = topology.for_sequence(p, tool=topology.TOOL_V2)
        assert v1 is not None and v1["tool_version"].startswith(topology.V1_TOOL_PREFIX)
        t1 = v1["per_residue_topology"]
        assert ov or (t1[0] == "I" and t1[-1] == "O"), (
            f"{m['candidate_id']}: not N-in / C-out by v1"
        )
        t2 = v2.get("dtm2_topology_string") if v2 is not None else None
        t2 = t2 if isinstance(t2, str) else ""
        side2 = v2_side(t2[-1:])
        up = uniprot_flags(p, uniprot.get(m["candidate_id"]))
        flags = list(up["flags"])
        if m["gene_symbol"] in gpi:
            flags.append("gpi")
        if not up["uniprot_acc"]:
            flags.append("uniprot_not_checked")
        if up["matches"] == "no":
            flags.append("cterm_differs_from_uniprot")
        if up["canon_len"] and len(p) < PARTIAL_FRACTION * up["canon_len"]:
            flags.append("partial_form")
        tail = trailing(t1, "O")
        if ov:
            flags.append("lane_override_topology")
        elif tail < SHORT_TAIL:
            flags.append("short_ecto_tail")
        if v2 is None:
            flags.append("no_v2")
        elif "M" not in t2:
            flags.append("v2_no_tm")
        elif side2 == "cytoplasmic":
            flags.append("v2_cterm_cytoplasmic")
        plddt, model = terminal_plddt(m["candidate_id"], payload, cache)
        mm = member_model(m["candidate_id"], payload, cache)
        burial = terminal_burial(mm[1]) if mm else None
        s_first, s_last, s_site = stalk(t1, mm[0]) if mm and not ov else (None, None, None)
        if plddt is not None and plddt >= STRUCTURED_PLDDT:
            flags.append("structured_terminus")
        hold = [f for f in flags if f in ("cterm_peptide_removed", "gpi")]
        out.append(
            {
                "candidate_id": m["candidate_id"],
                "design_id": design_ids()[m["candidate_id"]],
                "gene_symbol": m["gene_symbol"],
                "candidate_set": m["candidate_set"],
                "species": m["species"],
                "length_aa": len(p),
                "n_tm_v1": m["n_tm"],
                "cterm_ecto_aa_v1": tail,
                "v2_structural_type": (v2 or {}).get("dtm2_structural_type") or "",
                "v2_cterm_side": side2 or ("no TM" if v2 is not None else "no v2 row"),
                # Evidence for the C-out call (build_l1_cterm_arms.cout_disposition): the length
                # of v1's last TM helix, and UniProt's side at this member's C-terminus.
                "last_tm_aa_v1": len((re.findall("M+", t1) or [""])[-1]),
                "uniprot_cterm_topology": uniprot_cterm_topology(p, uniprot.get(m["candidate_id"])),
                # A cleaved signal peptide and the number of TM helices by UniProt: together
                # with a cytoplasmic C-terminus they describe a type I protein.
                **uniprot_architecture(uniprot.get(m["candidate_id"])),
                "uniprot_acc": up["uniprot_acc"],
                "uniprot_length_aa": up["canon_len"] or "",
                "cterm_in_uniprot": up["matches"],
                "terminal_plddt": "" if plddt is None else f"{plddt:.1f}",
                "terminal_ca_neighbors_10a": "" if burial is None else burial,
                "stalk_residues": (
                    ""
                    if s_first is None
                    else (f"{s_first}-{s_last}" if s_last >= s_first else "none")
                ),
                "stalk_site_after": "" if s_site is None else s_site,
                "model": model,
                "shed_risk": m["shed_risk"],
                "flags": ";".join(flags),
                "cterm_tag": "removed" if hold else "kept",
                "evidence": " | ".join(up["evidence"]),
                "last_10_aa": p[-10:],
                "sequence_sha256": topology.sequence_sha256(p),
            }
        )

    # No UniProt entry: name the human member of the same gene and its flags, so the gap is
    # visible per member. A removed human C-terminal peptide would make keeping it unsafe.
    human = {r["gene_symbol"]: r for r in out if r["species"] == "human" and r["uniprot_acc"]}
    for r in out:
        if "uniprot_not_checked" not in r["flags"].split(";"):
            continue
        h = human.get(r["gene_symbol"])
        if h is None:
            r["evidence"] = "no UniProt entry; processing checks not run; no human member"
            continue
        hflags = [f for f in h["flags"].split(";") if f in UNIPROT_PROCESSING_FLAGS]
        if "cterm_peptide_removed" in hflags:
            raise SystemExit(
                f"{r['design_id']}: no UniProt entry, and human {h['design_id']}'s C-terminal "
                "peptide is removed; check this ortholog's own C-terminus before keeping it"
            )
        r["evidence"] = (
            f"no UniProt entry; processing checks not run; human member {h['design_id']} "
            f"({h['uniprot_acc']}): {', '.join(hflags) or 'no processing flags'}"
        )

    with args.out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out)
    c = Counter(f for r in out for f in r["flags"].split(";") if f)
    print(
        f"{len(out)} C-terminal members; C-terminal tag kept "
        f"{sum(r['cterm_tag'] == 'kept' for r in out)}, removed "
        f"{sum(r['cterm_tag'] == 'removed' for r in out)}"
    )
    for f, n in c.most_common():
        print(f"  {f:<28} {n}")
    for r in out:
        if r["cterm_tag"] == "removed":
            print(f"  C-TAG REMOVED {r['candidate_id']} {r['gene_symbol']}: {r['evidence']}")
    unarchived = [r for r in out if r["model"] == NOT_ARCHIVED]
    for r in unarchived:
        print(f"  NATIVE MODEL NOT ARCHIVED {r['design_id']} {r['gene_symbol']} {r['species']}")
    if args.check:
        assert args.out.read_bytes() == OUT.read_bytes(), f"{OUT.name} does not reproduce"
        print(f"{OUT.name} reproduces byte for byte ({len(out)} members)")
        return 0
    print(f"wrote {args.out.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
