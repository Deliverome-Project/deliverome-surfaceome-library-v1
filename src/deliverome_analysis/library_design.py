"""What we are actually ordering, and which records are about that exact protein.

This module exists because three different identities kept getting conflated, each time
producing an answer that looked right:

* **Which library tranche** -- ``l1`` or ``l2``. Written ``L1``/``L2`` in conversation,
  which is one keystroke from a *tool* version. :data:`LIBRARIES` are the only accepted
  values and they never abbreviate to a bare number.
* **Which topology** -- a per-residue string is produced by a named predictor
  (``topology.TOOL_V1`` / ``TOOL_V2``) under a ``topology_version`` within a ``cohort``.
  All three travel together or the string is unattributable.
* **Which protein** -- a library member's protein is its *designed CDS translated*, which
  is an Ensembl-transcript product. It usually corresponds to a UniProt **isoform**
  record, not the canonical one. So a record fetched by canonical accession -- an
  AlphaFold model, a topology row -- is frequently about a different protein, with a
  plausible length and plausible helix count. Everything here matches on the sequence.

Worked numbers behind those warnings, as of 2026-10-14: of the S4E9 L1 members whose
designed ORF resolved, 119 of 121 matched under ``human_isoforms`` and only 2 under
``human_canonical``; 200 L1 members had an AlphaFold model of a different protein, 23 of
them at identical length. Length is not a proxy for identity -- equal length proves
nothing, only an equal sequence does.

    from deliverome_analysis import library_design as ld
    from deliverome_analysis import topology

    prot = ld.designed_proteins(manifest, library=ld.LIBRARY_L1)
    report = ld.topology_report(prot, tool=topology.TOOL_V1)
    report.loc[~report["found"], "candidate_id"]      # needs a v1 run

The manifest path is a parameter rather than a constant because the tranche assignment
lives with the experiment that made it; the designed CDS is pinned here because it is the
one durable record of what gets synthesised.
"""

from __future__ import annotations

import csv
import json
import re
from collections.abc import Mapping
from pathlib import Path

import pandas as pd

from deliverome_analysis import topology

#: Library tranches. Not tool versions -- see :data:`topology.TOOL_V1`.
LIBRARY_L1 = "l1"
LIBRARY_L2 = "l2"
LIBRARIES: tuple[str, ...] = (LIBRARY_L1, LIBRARY_L2)

#: Designed coding sequences for every candidate, keyed by the manifest's ``key``.
#: The authoritative record of what the library synthesises -- not UniProt, not Ensembl
#: REST, both of which return a *different* proteoform for many members.
#:
#: Append-only by version: each adoption publishes the previous record plus its new
#: entries under a new dated key, with every earlier entry unchanged, and this constant
#: moves to it. 2026-11-10: the ortholog CDS-policy panel (s4e9
#: ``adopt_ortholog_cds_policy.py``) appended 56 ortholog forms to the
#: :data:`BASE_CDS_KEY` record. 2026-11-12: the approved reference fixes (s4e9
#: ``apply_reference_fixes.py``; ledger ``reference_changes.tsv``) replaced seven ``lib:``
#: entries -- five named reference forms and two ATG start replacements -- and changed no
#: other entry; the 2026-11-10 record stays at its own key. 2026-11-12, second round
#: (``reference-fixes-2``): three more named references (LMLN, DISP2, TSPAN10), the OR8G3
#: ATG start and EPHA6's in-frame trim to Q9UF33; the first round's record stays at its key.
#: Then (``ortholog-start-trims``) cyno CCR5 restarted at the native Met aligned to human Met1
#: (s4e9 ``apply_ortholog_start_trims.py``); every other entry unchanged. Then
#: (``ortholog-start-trims-2``) cyno CHRM3, NIPAL4 and OR3A2 restarted the same way, once
#: their DeepTMHMM v1 rows were published; the CCR5-only record stays at its key.
CDS_KEY = (
    "sprint4/s4e9-orf-design-fg-panel/library-candidates/"
    "2026-11-12-ortholog-start-trims-2/candidate_cds.json"
)

#: Domesticated ORF per candidate (``construct``, edit counts), keyed like :data:`CDS_KEY`
#: and versioned with it.
CONSTRUCTS_KEY = (
    "sprint4/s4e9-orf-design-fg-panel/library-candidates/"
    "2026-11-12-ortholog-start-trims-2/candidate_constructs.json"
)

#: The records :data:`CDS_KEY` and :data:`CONSTRUCTS_KEY` extend. They stay the published
#: output of the pre-adoption stages (``library_pipeline.py`` replays those against them).
BASE_CDS_KEY = "sprint4/s4e9-orf-design-fg-panel/library-candidates/candidate_cds.json"
BASE_CONSTRUCTS_KEY = (
    "sprint4/s4e9-orf-design-fg-panel/library-candidates/candidate_constructs.json"
)

#: Columns the MEMBER owns. A topology row carrying any of these has it renamed
#: ``topology_*``, so "which species is this row about" is never one ambiguous column.
_IDENTITY_COLUMNS = frozenset({"species", "gene_symbol", "library", "candidate_set"})

#: Columns :func:`topology_report` always emits, before the topology row's own columns.
REPORT_COLUMNS = (
    "candidate_id",
    "library",
    "candidate_set",
    "species",
    "gene_symbol",
    "found",
    "tool",
    "length",
    "sequence_sha256",
)


def gene_from_key(key: str) -> str:
    """Gene symbol out of a manifest ``key``.

    The field position is NOT constant across candidate sets, which is the trap:

        lib:ADAM12                  -> ADAM12   (field 1)
        iso:CADM1:Q9BY67            -> CADM1    (field 1)
        rescue:ADGRD2               -> ADGRD2   (field 1)
        ort:cyno:ABCB11:D2CFL3      -> ABCB11   (field 2 -- field 1 is the SPECIES)

    Taking field 1 unconditionally returns "cyno"/"mouse" for every ortholog, which then
    matches nothing in any gene-keyed join. That silently dropped GPI and shedding
    annotation for all 196 L1 ortholog members.
    """
    parts = key.split(":")
    if len(parts) < 2:
        return ""
    if parts[0] == "ort":
        return parts[2] if len(parts) > 2 else ""
    return parts[1]


def _check_library(library: str) -> str:
    if library not in LIBRARIES:
        raise ValueError(
            f"unknown library {library!r}; use library_design.LIBRARY_L1 "
            f"({LIBRARY_L1!r}) or LIBRARY_L2 ({LIBRARY_L2!r}). "
            "These are library tranches, not DeepTMHMM tool versions "
            f"(those are topology.TOOL_V1 / TOOL_V2)."
        )
    return library


#: Manifest column naming why a member was retired. A retired member keeps its design id
#: as its ``candidate_id`` and gets no member id, and is no longer built, so every reader
#: here skips it. See s4e9 ``dedupe_library_members.py``.
RETIRED_COLUMN = "retired_reason"

#: Manifest column holding a row's **design id** (``l1-pNNNN`` / ``l2-pNNNN``): the id the
#: row was designed under, permanent and never reused. Design-stage records (panel changes,
#: retirements, adoption ledgers, fold and topology runs) are keyed on it.
DESIGN_ID_COLUMN = "design_id"

#: An active L1 member's ``candidate_id`` is its **member id**, ``l1-mNNNN``, numbered in
#: order of synthesized-insert length (= plate order, 96 per plate, column-major) by s4e9
#: ``renumber_library_members.py``. It is a different spelling from the design id on
#: purpose, so an old ``l1-p`` id can never be read as a member id or vice versa. Member
#: ids are reassigned whenever membership changes, until the order freeze.
MEMBER_ID_RE = re.compile(r"l1-m[0-9]{4}")
DESIGN_ID_RE = re.compile(r"l[12]-p[0-9]{4}")

#: From the first per-insert cut (agreed 2026-11-14/15 by the coordinator, l1-integration,
#: ranker-pu-segment and l1-reviewer-package, after Becca's rule that an internal member may
#: carry 2 or 3 arms) an L1 member can be ordered as several inserts, so three ids replace the
#: member id:
#:
#: * the design id (``l1-pNNNN``) names the member (the protein design); it is also the
#:   member's ``candidate_id``, so a member id no longer moves at a cut;
#: * the **arm id** (``l1-pNNNN-aK``) names one tagged insert design of it, permanently: K is
#:   assigned in creation order and never reused or re-ranked (role and rank are columns);
#: * the **insert id** (``l1-iNNNN``) is the plate position of one ordered insert, numbered
#:   by ordered insert length then arm id at each cut (96 per plate, column-major).
#:
#: ``l1-m`` ids then survive only in dated records.
ARM_ID_RE = re.compile(r"(l[12]-p[0-9]{4})-a([1-9][0-9]*)")
INSERT_ID_RE = re.compile(r"l1-i[0-9]{4}")


def arm_id(design: str, k: int = 1) -> str:
    """The arm id of a design's K-th arm (``l1-p0325``, 2 -> ``l1-p0325-a2``)."""
    if not DESIGN_ID_RE.fullmatch(design) or k < 1:
        raise ValueError(f"not a design id and arm number: {design!r}, {k!r}")
    return f"{design}-a{k}"


def design_of_arm(arm: str) -> str:
    """The design id an arm id belongs to."""
    m = ARM_ID_RE.fullmatch(arm)
    if not m:
        raise ValueError(f"not an arm id: {arm!r}")
    return m.group(1)


def insert_id(n: int) -> str:
    """The insert id at plate position ``n`` (1-based)."""
    if not 1 <= n <= 9999:
        raise ValueError(f"insert position out of range: {n}")
    return f"l1-i{n:04d}"


def design_id(row: Mapping[str, str]) -> str:
    """A manifest row's design id, which design-stage code must key on.

    Falls back to ``candidate_id`` for a manifest written before the renumbering (whose
    ids were all design ids), and refuses a member id with no design id beside it.
    """
    did = (row.get(DESIGN_ID_COLUMN) or "").strip() or row["candidate_id"]
    if not DESIGN_ID_RE.fullmatch(did):
        raise ValueError(f"row has no design id: {row.get('candidate_id')!r} / {did!r}")
    return did


def is_active(row: Mapping[str, str]) -> bool:
    """False for a manifest row that has been retired from its tranche."""
    return not (row.get(RETIRED_COLUMN) or "").strip()


def design_ids(manifest: str | Path) -> dict[str, str]:
    """``{candidate_id: design_id}`` for every manifest row (identity for unnumbered rows).

    Use it wherever a member-keyed record must reach a design-keyed artefact: fold inputs
    and their structure files, topology supplement rows, panel and retirement ledgers.
    """
    with open(manifest) as fh:
        return {r["candidate_id"]: design_id(r) for r in csv.DictReader(fh, delimiter="\t")}


def load_designed_cds() -> dict[str, dict]:
    """The whole designed-CDS record at :data:`CDS_KEY`, keyed by manifest ``key``.

    Includes retired members' entries; the record is append-only, so a retired member's
    CDS stays readable.
    """
    from deliverome_analysis.s3 import processed

    with open(processed.download(CDS_KEY)) as fh:
        return json.load(fh)


def designed_cds(manifest: str | Path, *, library: str) -> dict[str, str]:
    """Designed coding sequence (DNA, as synthesised) per active ``candidate_id``."""
    _check_library(library)
    cds = load_designed_cds()
    with open(manifest) as fh:
        return {
            row["candidate_id"]: cds[row["key"]]["cds"]
            for row in csv.DictReader(fh, delimiter="\t")
            if row.get("library") == library and is_active(row) and row["key"] in cds
        }


def designed_proteins(manifest: str | Path, *, library: str) -> dict[str, dict[str, str]]:
    """Translated designed CDS for one library tranche, keyed by ``candidate_id``.

    Each value carries ``sequence`` plus the identity a caller needs to report what the
    sequence is: ``library``, ``candidate_set``, ``species``, ``gene_symbol``, ``key``,
    ``route`` and ``transcript``. Members with no designed CDS are omitted and listed in
    the returned mapping's absence -- they are not yet designed, and a caller that needs
    to know should compare against the manifest rather than be handed a blank sequence.
    Retired members (:data:`RETIRED_COLUMN` set) are omitted too.
    """
    _check_library(library)
    from deliverome_analysis import seqtools

    cds = load_designed_cds()

    out: dict[str, dict[str, str]] = {}
    with open(manifest) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("library") != library or not is_active(row):
                continue
            entry = cds.get(row["key"])
            if not entry:
                continue
            key = row["key"]
            out[row["candidate_id"]] = {
                "sequence": seqtools.translate(entry["cds"]).rstrip("*"),
                "key": key,
                "library": library,
                "candidate_set": row.get("set", ""),
                "species": row.get("species", ""),
                "gene_symbol": gene_from_key(key),
                "route": row.get("route", ""),
                "transcript": row.get("transcript", ""),
            }
    return out


def topology_report(
    members: Mapping[str, Mapping[str, str]],
    *,
    tool: str | None = topology.PREFERRED_TOOL,
    cohorts: tuple[str, ...] | None = None,
) -> pd.DataFrame:
    """Join designed proteins to their topology, matched on the exact sequence.

    One row per member, carrying the library tranche, the designed sequence's hash, and
    -- when found -- the ``cohort`` / ``topology_version`` / ``tool`` that produced the
    string. Rows with ``found == False`` have had no run of ``tool`` on that exact
    sequence; they keep their identity so the frame doubles as the work list.

    Defaults to :data:`topology.PREFERRED_TOOL`, so an un-parameterised call answers
    "which members lack a v1 string", which is the question that matters for display.
    """
    # Resolve per species. A sequence hash is species-blind and a human protein is often
    # identical to its cyno or mouse ortholog, so an unconstrained match can return the
    # other organism's row -- 66 of 5,888 members were matched across species before this.
    # An explicit `cohorts` argument from the caller wins, so this can still be overridden.
    frames = []
    by_species: dict[str, dict[str, str]] = {}
    for cid, m in members.items():
        by_species.setdefault(m.get("species", ""), {})[cid] = m["sequence"]
    for sp, seqs in by_species.items():
        allowed = cohorts if cohorts is not None else topology.cohorts_for_species(sp)
        frames.append(topology.resolve(seqs, tool=tool, cohorts=allowed))
    res = pd.concat(frames, ignore_index=True).set_index("key")

    rows = []
    for cid, m in members.items():
        r = res.loc[cid].to_dict()
        # The member's ordering identity comes LAST so the topology row cannot overwrite
        # it. It used to come first and be spread over: `topology_public` has its own
        # `species` and `gene_symbol`, so a matched row replaced the member's with the
        # matched protein's -- which is how this frame reported `cynomolgus` and
        # `G7PPA4_MACFA` for a human library reference. Any column the topology row shares
        # with the member identity is renamed `topology_*`, so the collision cannot recur.
        topo_row = {
            (f"topology_{k}" if k in _IDENTITY_COLUMNS else k): v
            for k, v in r.items()
            if k != "key"
        }
        rows.append(
            {
                **topo_row,
                "candidate_id": cid,
                "library": m.get("library", ""),
                "candidate_set": m.get("candidate_set", ""),
                "species": m.get("species", ""),
                "gene_symbol": m.get("gene_symbol", ""),
            }
        )
    out = pd.DataFrame.from_records(rows)
    ordered = [
        *[c for c in REPORT_COLUMNS if c in out.columns],
        *[c for c in out.columns if c not in REPORT_COLUMNS],
    ]
    return out[ordered]


def coverage(report: pd.DataFrame) -> pd.DataFrame:
    """Per-library, per-set found/missing counts from a :func:`topology_report` frame.

    Reports the tranche and the tool on every row so a pasted summary cannot be read as
    being about the other library or the other predictor.
    """
    g = (
        report.groupby(["library", "candidate_set", "tool"], dropna=False)["found"]
        .agg(found="sum", total="count")
        .reset_index()
    )
    g["missing"] = g["total"] - g["found"]
    return g.sort_values(["library", "candidate_set"]).reset_index(drop=True)
