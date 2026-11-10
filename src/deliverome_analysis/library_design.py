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
#: :data:`BASE_CDS_KEY` record.
CDS_KEY = (
    "sprint4/s4e9-orf-design-fg-panel/library-candidates/"
    "2026-11-10-ortholog-cds-policy/candidate_cds.json"
)

#: Domesticated ORF per candidate (``construct``, edit counts), keyed like :data:`CDS_KEY`
#: and versioned with it.
CONSTRUCTS_KEY = (
    "sprint4/s4e9-orf-design-fg-panel/library-candidates/"
    "2026-11-10-ortholog-cds-policy/candidate_constructs.json"
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


#: Manifest column naming why a member was retired. A retired member keeps its
#: ``candidate_id`` -- ids are never reused or renumbered -- but is no longer built, so
#: every reader here skips it. See s4e9 ``dedupe_l1_members.py``.
RETIRED_COLUMN = "retired_reason"


def is_active(row: Mapping[str, str]) -> bool:
    """False for a manifest row that has been retired from its tranche."""
    return not (row.get(RETIRED_COLUMN) or "").strip()


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
