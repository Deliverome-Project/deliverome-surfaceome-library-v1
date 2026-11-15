"""DeepTMHMM topology for a designed protein sequence, matched on the sequence itself.

**Match on the sequence, never on the accession.** An accession identifies a gene; it does
not identify the protein a construct expresses. Our library ORFs are Ensembl transcripts,
and the protein an Ensembl transcript encodes usually corresponds to a UniProt *isoform*
record rather than the canonical one -- of the S4E9 L1 members whose designed ORF resolves
at all, 119 of 121 in one audit matched under ``human_isoforms`` and only 2 under
``human_canonical``. So looking a gene up by its canonical accession returns the topology
of a *different protein*, with a plausible length and plausible helix count, and any tag
site numbered against it is numbered against a sequence the construct does not contain.
The same fact is why AlphaFold DB served the wrong model for 200 L1 members: AFDB serves
the canonical too.

Keying on the sha256 of the sequence makes the answer exact or absent, never wrong. When a
sequence has no record this module returns ``None`` / ``found=False``. It never falls back
to the canonical, the nearest length, or the gene's other isoforms. If you want a
same-gene alternative you must ask for it explicitly and label it as not the ordered form.

Reads a pinned snapshot of ``topology_public`` from ``deliverome-external`` -- see
:data:`SNAPSHOT_KEY`. The snapshot exists because ``sequence`` is not indexed in D1, so a
sequence join there is a full scan of ~71k long TEXT rows: acceptable offline, wrong for a
serve path. Produced by ``scripts/cloud/export_topology_snapshot.py`` in
accessible-surfaceome. Downloads through the cached S3 helper, so repeat calls are free.

Cohorts in the snapshot are proteoform spaces, not datasets: ``human_canonical``,
``human_isoforms``, ``cyno_ortholog``, ``mouse_ortholog``, their ``*_repick`` variants, and
``library_designed_orf`` for the designed ORFs of library members that exist in none of the
others. Rows carry both alphabets -- ``per_residue_topology`` (DeepTMHMM v1, the string the
structure viewer colours by) and ``dtm2_topology_string`` (v2). ``dtm2_v1_alphabet_lossy``
is a v2->v1 translation and is lossy by name: prefer the real v1 string when present.

    from deliverome_analysis import topology

    rec = topology.for_sequence(designed_protein)
    if rec is None:
        ...  # no topology for THIS sequence; do not substitute the canonical
    else:
        rec["per_residue_topology"], rec["cohort"], rec["topology_version"]

    report = topology.resolve({cid: seq for cid, seq in designed.items()})
    gap = report.loc[~report["found"], "key"].tolist()
"""

from __future__ import annotations

import hashlib
from collections.abc import Mapping
from functools import lru_cache

import pandas as pd

# --- pinned topology snapshot. Bump THIS to adopt a newer export. ---------------------
# 2026-11-12: plus 3 topo_2026_11_12_cyno_restarts rows (library_designed_orf,
# deeptmhmm-1.0-academic) for the decided 5' restarts of cyno CHRM3, NIPAL4 and OR3A2
# (s4e9 ortholog_start_trims.tsv; publish_topology_supplement.py, run files under
# processed://sprint4/s4e9-orf-design-fg-panel/topology-runs/2026-11-12_cyno_restarts.*).
# 71,700 rows; every 2026-11-06 row is preserved. Bumped by l1-reviewer-package, PR #283.
SNAPSHOT_KEY = "accessible-surfaceome/topology/2026-11-12-cyno-restarts/topology.parquet"
#: SHA-256 of the object at :data:`SNAPSHOT_KEY`.
SNAPSHOT_SHA256 = "c8fb5874437acd0677da76e2864aa76089ac7fa072f3defd2338f15f3890292e"

# 2026-11-10: the ortholog CDS-policy snapshot, adopted with that panel (s4e9
# adopt_ortholog_cds_policy.py). It is the 2026-11-01 export (topology_public with the
# library_designed_orf cohort, plus three topo_2026_11_01_l1_designed rows) plus the 323
# topo_2026_11_05_l1_compara116 rows and the 58 topo_2026_11_05_ortholog_cds_policy rows
# (deeptmhmm-1.0-academic), 71,697 rows; every 2026-11-01 row is preserved. Adopting it
# changed no existing L1 member's matched row; 166 L2 orthologs gained a v1 match.
ORTHOLOG_CDS_POLICY_SNAPSHOT_KEY = (
    "accessible-surfaceome/topology/2026-11-06-ortholog-cds-policy/topology.parquet"
)
ORTHOLOG_CDS_POLICY_SNAPSHOT_SHA256 = (
    "7683d75a3e04228d95fcfb8fe9b7e9b46b38c0ea254ca42550a666f5d8bbf216"
)

#: The snapshot before 2026-11-10, kept so an analysis run against it can be replayed.
PREVIOUS_SNAPSHOT_KEY = "accessible-surfaceome/topology/2026-11-01/topology.parquet"

# The first Compara 116 replacement panel's 323 v1 predictions (also in SNAPSHOT_KEY).
# Both D1 databases also carry these rows as topo_2026_11_05_l1_compara116.
ORTHOLOG_COMPARA116_SNAPSHOT_KEY = (
    "accessible-surfaceome/topology/2026-11-05-l1-compara116/topology.parquet"
)

# --- DeepTMHMM TOOL versions. NOT library tranches. -----------------------------------
# "L1"/"L2" are *library tranches* (which candidates we order). "v1"/"v2" are *DeepTMHMM
# tool versions* (which predictor produced a topology string). They are unrelated axes,
# and a bare "1" or "2" in a call site cannot say which is meant -- that ambiguity has
# already cost one round of confusion, so these names never abbreviate to a number.
TOOL_V1 = "deeptmhmm_v1"
TOOL_V2 = "deeptmhmm2"

#: What this repo treats as the topology of record, and why.
#:
#: 1. The accessible-surfaceome structure viewer colours by ``per_residue_topology``,
#:    which only the v1 run fills.
#: 2. ``dtm2_v1_alphabet_lossy`` -- v2 projected back onto the v1 alphabet -- is lossy by
#:    name and was flagged lossy for 33 of the 157 designed ORFs published on 2026-10-14.
#:    So v2 is not a drop-in substitute for a real v1 string for a fifth of them.
#:
#: v2 is kept alongside (richer membrane-type and reentrant/interfacial calls) and is
#: worth reporting where the two disagree, but it does not replace v1.
PREFERRED_TOOL = TOOL_V1

#: Column each tool's per-residue string lives in. Different alphabets -- never compare
#: a v1 string to a v2 string character by character.
TOOL_COLUMN = {TOOL_V1: "per_residue_topology", TOOL_V2: "dtm2_topology_string"}

#: Column the structure viewer colours by.
VIEWER_TOPOLOGY_COLUMN = TOOL_COLUMN[TOOL_V1]

#: ``tool_version`` prefix that identifies a row produced by the v1 predictor.
#:
#: **A non-null ``per_residue_topology`` does NOT mean a v1 run.** A dtm2 row carries a
#: v1-alphabet string too -- ``dtm2_v1_alphabet_lossy``, the v2 output projected back --
#: and it differs from the genuine v1 string for 4,644 of the 11,389 human_canonical rows
#: that exist in both versions. Selecting on "has a v1 string" therefore silently accepts
#: the projection for ~41% of rows. The provenance stamp is the only honest discriminator,
#: so that is what :func:`_with_tool` filters on.
V1_TOOL_PREFIX = "deeptmhmm-"

#: ``tool_version`` prefix for the v2 predictor.
V2_TOOL_PREFIX = "deeptmhmm2"

#: Substring marking a ``topology_version`` as provisional rather than a published run.
#: ``topo_test_optA_full`` is a parameter-sweep trial. It is ranked LAST, never excluded:
#: all 280 of its sequences also have a published v1 row, so it is never needed -- but it
#: sorts after every ``topo_2026_*`` name alphabetically, so a newest-version tiebreak
#: picked it for 203 library members before this existed.
PROVISIONAL_VERSION_MARKER = "test"

#: ``tool_version`` values that are not DeepTMHMM output at all -- UniProt feature
#: fallbacks in the mouse and cyno cohorts. They satisfy NEITHER tool: a caller asking
#: for a predictor's topology must not be handed a row no predictor produced.
NON_PREDICTOR_TOOL_VERSIONS = frozenset({"uniprot_features_v1"})

#: Which cohorts may answer for a member of a given species.
#:
#: A sequence hash is species-blind, and a human protein is often byte-identical to its
#: cyno or mouse ortholog over the matched region -- so an unconstrained match can hand a
#: human member its cyno ortholog's topology row. That happened: 66 of 5,888 library
#: members were matched across species, 43 of them human -> cyno. The topology string is
#: usually the same, which is exactly why it went unnoticed, but the row's species,
#: accession and entry name are then another organism's.
#:
#: ``library_designed_orf`` is in every list because its rows are keyed per candidate and
#: hold all three species; a match there is already member-specific.
SPECIES_COHORTS: dict[str, tuple[str, ...]] = {
    "human": ("human_canonical", "human_isoforms", "library_designed_orf"),
    "cyno": (
        "cyno_ortholog",
        "cyno_ortholog_repick",
        "cyno_ortholog_cds_policy",
        "library_designed_orf",
    ),
    "mouse": (
        "mouse_ortholog",
        "mouse_ortholog_repick",
        "mouse_ortholog_cds_policy",
        "library_designed_orf",
    ),
}

#: Spellings seen in the manifest and in topology_public for the same organism.
_SPECIES_ALIASES = {
    "cynomolgus": "cyno",
    "macaca_fascicularis": "cyno",
    "mus_musculus": "mouse",
    "homo_sapiens": "human",
}


def normalize_species(species: str) -> str:
    """Fold the spellings of one organism together (``cynomolgus`` -> ``cyno``).

    Both labels occur: the manifest says ``cyno`` and ``topology_public`` says
    ``cynomolgus``, so a join or a groupby on the raw string silently splits one organism
    in two.
    """
    s = str(species).strip().lower()
    return _SPECIES_ALIASES.get(s, s)


def cohorts_for_species(species: str) -> tuple[str, ...] | None:
    """Allowed cohorts for a species, or None when the species is unknown.

    None means "do not constrain" rather than "match nothing": an unrecognised species
    should degrade to the old species-blind behaviour and be visible, not silently
    resolve to zero rows.
    """
    return SPECIES_COHORTS.get(normalize_species(species))


#: Returned by :func:`resolve` for every input, ahead of the snapshot's own columns.
_KEY_COLUMNS = ("key", "found", "tool", "length", "sequence_sha256")


def _tool_column(tool: str) -> str:
    """Validate a tool name and return its per-residue column.

    Rejects anything that is not one of the two tool constants -- including the bare
    ``1``/``2`` or ``"l1"``/``"v1"`` muddle this module exists to prevent.
    """
    if tool not in TOOL_COLUMN:
        raise ValueError(
            f"unknown topology tool {tool!r}; use topology.TOOL_V1 "
            f"({TOOL_V1!r}) or topology.TOOL_V2 ({TOOL_V2!r}). "
            "Note these are DeepTMHMM tool versions, not library tranches (l1/l2)."
        )
    return TOOL_COLUMN[tool]


def normalize(aa: str) -> str:
    """Canonical form for matching: stripped, upper-case, no trailing stop.

    The snapshot stores translated coding sequences without a stop codon, so a sequence
    that still carries ``*`` would hash differently and silently miss. Internal ``*`` is
    left alone -- it means the translation is wrong and should fail loudly downstream
    rather than be quietly repaired here.
    """
    return str(aa).strip().upper().rstrip("*")


def sequence_sha256(aa: str) -> str:
    """sha256 of :func:`normalize`'d sequence -- the snapshot's join key."""
    return hashlib.sha256(normalize(aa).encode()).hexdigest()


@lru_cache(maxsize=1)
def load_snapshot() -> pd.DataFrame:
    """Load the pinned topology snapshot. Cached in-process.

    One row per (``topology_version``, ``cohort``, ``uniprot_acc_full``), so a sequence
    present in several cohorts or several tool versions has several rows. Use
    :func:`for_sequence` to get the single best row, or pass ``all_rows=True`` to see them.
    """
    from deliverome_analysis.s3 import external

    df = pd.read_parquet(external.download(SNAPSHOT_KEY))
    if "sequence_sha256" not in df.columns:
        raise ValueError(
            f"{SNAPSHOT_KEY} has no sequence_sha256 column -- it predates "
            "sequence-keyed matching and cannot be used here"
        )
    return df


def _with_tool(df: pd.DataFrame, tool: str | None) -> pd.DataFrame:
    """Keep only rows that ``tool`` actually produced.

    Filters on the ``tool_version`` provenance stamp, not on which string column is
    populated. A dtm2 row carries a v1-alphabet string as well -- the lossy back
    projection -- so "has a ``per_residue_topology``" is true for both predictors and
    selecting on it quietly accepts the projection (it differs from the real v1 string for
    4,644 of 11,389 rows present in both versions). Asking for one predictor must never
    return the other's output relabelled.

    ``tool=None`` means "any row", and is the only way to get that laxer behaviour.
    """
    if tool is None:
        return df
    col = _tool_column(tool)
    if col not in df.columns or "tool_version" not in df.columns:
        return df.iloc[0:0]
    tv = df["tool_version"].fillna("")
    produced_by = (
        tv.str.startswith(V1_TOOL_PREFIX) & ~tv.str.startswith(V2_TOOL_PREFIX)
        if tool == TOOL_V1
        else tv.str.startswith(V2_TOOL_PREFIX)
    )
    return df[df[col].notna() & produced_by & ~tv.isin(NON_PREDICTOR_TOOL_VERSIONS)]


def _rank(df: pd.DataFrame, tool: str | None = None) -> pd.DataFrame:
    """Order candidate rows best-first: asked-for tool present, then the other, then newest."""
    primary = _tool_column(tool) if tool else VIEWER_TOPOLOGY_COLUMN
    other = TOOL_COLUMN[TOOL_V2] if primary == TOOL_COLUMN[TOOL_V1] else TOOL_COLUMN[TOOL_V1]
    out = df.copy()
    out["_p"] = out.get(primary, pd.Series(dtype=object)).notna()
    out["_o"] = out.get(other, pd.Series(dtype=object)).notna()
    # Published runs rank above provisional ones regardless of how the names sort.
    out["_published"] = ~out["topology_version"].fillna("").str.contains(PROVISIONAL_VERSION_MARKER)
    return out.sort_values(
        ["_p", "_o", "_published", "topology_version"],
        ascending=[False, False, False, False],
    ).drop(columns=["_p", "_o", "_published"])


def for_sequence(
    aa: str,
    *,
    tool: str | None = PREFERRED_TOOL,
    cohorts: tuple[str, ...] | None = None,
    all_rows: bool = False,
) -> dict | pd.DataFrame | None:
    """Topology for exactly this sequence from ``tool``, or ``None`` if there is none.

    ``None`` means no run of that tool has ever been done on this sequence. It does not
    mean the protein has no topology, and it is never a reason to substitute another
    proteoform's record. Defaults to :data:`PREFERRED_TOOL` (v1), so the default answer is
    "the v1 string or nothing" -- pass ``tool=TOOL_V2`` for v2 or ``tool=None`` for either.
    Restrict to particular proteoform spaces with ``cohorts``. With ``all_rows`` returns
    every matching row, best-first, instead of the best one.
    """
    df = load_snapshot()
    hit = df[df["sequence_sha256"] == sequence_sha256(aa)]
    if cohorts is not None:
        hit = hit[hit["cohort"].isin(cohorts)]
    hit = _with_tool(hit, tool)
    if hit.empty:
        return None if not all_rows else hit
    ranked = _rank(hit, tool)
    return ranked.reset_index(drop=True) if all_rows else ranked.iloc[0].to_dict()


def resolve(
    sequences: Mapping[str, str],
    *,
    tool: str | None = PREFERRED_TOOL,
    cohorts: tuple[str, ...] | None = None,
    snapshot: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Resolve many sequences at once. One row per input, in input order.

    Always returns a row for every key, with ``found`` saying whether this exact sequence
    has topology from ``tool``, and ``tool`` echoed back on every row so a saved report
    can never be read as being about the other predictor. Unmatched inputs keep their key,
    length and hash so the result doubles as the work list for a run -- that is the point
    of not dropping them.

    ``snapshot`` supplies an already checksummed input for replaying an older
    analysis after the shared default snapshot advances. It receives the same
    tool-provenance and cohort filters as the default snapshot.
    """
    df = load_snapshot() if snapshot is None else snapshot
    if cohorts is not None:
        df = df[df["cohort"].isin(cohorts)]
    df = _with_tool(df, tool)
    ranked = _rank(df, tool).drop_duplicates(subset="sequence_sha256", keep="first")
    lookup = ranked.set_index("sequence_sha256")

    rows = []
    for key, aa in sequences.items():
        norm = normalize(aa)
        digest = hashlib.sha256(norm.encode()).hexdigest()
        base = {
            "key": key,
            "found": digest in lookup.index,
            "tool": tool or "any",
            "length": len(norm),
            "sequence_sha256": digest,
        }
        if base["found"]:
            base |= lookup.loc[digest].to_dict()
        rows.append(base)

    out = pd.DataFrame.from_records(rows)
    ordered = [*_KEY_COLUMNS, *[c for c in out.columns if c not in _KEY_COLUMNS]]
    return out[ordered]


def missing(
    sequences: Mapping[str, str],
    *,
    tool: str | None = PREFERRED_TOOL,
    cohorts: tuple[str, ...] | None = None,
) -> list[str]:
    """Keys whose exact sequence has no ``tool`` topology record -- the run work list.

    Defaults to :data:`PREFERRED_TOOL`, so "what still needs running" means "what still
    needs a v1 run" unless you say otherwise.
    """
    report = resolve(sequences, tool=tool, cohorts=cohorts)
    return report.loc[~report["found"], "key"].tolist()


def same_gene_alternatives(aa: str, gene_symbol: str) -> pd.DataFrame:
    """Other proteoforms of the same gene, for when the ordered sequence has no record.

    Explicitly separate from :func:`for_sequence` so that a near-miss can never be
    mistaken for the ordered form. Every row here is a *different protein* from ``aa``;
    label it as such wherever it is shown, and never number a tag site against it.
    """
    df = load_snapshot()
    digest = sequence_sha256(aa)
    alt = df[(df["gene_symbol"] == gene_symbol) & (df["sequence_sha256"] != digest)]
    if alt.empty:
        return alt
    return _rank(alt).reset_index(drop=True)
