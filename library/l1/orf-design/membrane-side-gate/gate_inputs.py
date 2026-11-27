"""Inputs to the membrane-side gate: members, their proteins, topology, models and picks.

Everything is read from pinned records, never from a working tree that can move:

* **Members.** The adopted L1 panel is on branch ``claude/l1-reviewer-package-aa447f``
  (the Compara 116 ortholog panel adopted on 2026-11-10), not yet on ``main``. Its
  manifest, tag-path table and structure payload are read with ``git show`` at
  :data:`PANEL_COMMIT`, together with the designed-CDS record and the topology snapshot
  that commit pins. When that branch merges, point :data:`PANEL_COMMIT` at the merge.
* **v4 picks.** The S5E13 v4 ranker's member picks, ``tag_sites_v4_members.tsv`` at
  :data:`S5E13_COMMIT` on the S5E13 branch (pull request 180, closed; the three files are
  vendored in ``ranker-c7be01a4/`` with their SHA-256), keyed to members through
  that commit's ``l1_member_targets.tsv``.
* **Built loop arms.** ``constructs_loop.tsv`` in this experiment folder.
* **API lanes.** The accessible-surfaceome tag-site API's deterministic internal lanes
  (``disorder``, ``surface_loop``), fetched live and cached with the fetch date. They are
  numbered on the human canonical sequence and carry ``residue_before``/``residue_after``,
  which are checked when a site is carried onto a member.

Models, API responses, UniProt sequences and TmAlphaFold entries are cached under
``~/.cache/deliverome-analysis/membrane-side-gate/``.
"""

from __future__ import annotations

import gzip
import hashlib
import io
import json
import re
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

import pandas as pd

from deliverome_analysis.cloning import seqtools
from deliverome_analysis.surfaceome import topology

HERE = Path(__file__).resolve().parent
EXPERIMENT = HERE.parent
REPO = HERE.parents[3]
CACHE = Path.home() / ".cache" / "deliverome-analysis" / "membrane-side-gate"
UA = "deliverome-analysis membrane-side-gate"

#: ``claude/l1-reviewer-package-aa447f`` at 6bf341f2: the adopted ortholog panel, renumbered
#: to member ids ``l1-mNNNN`` (``design_id`` keeps the old ``l1-pNNNN``), with its second pass
#: (11 more internal-rung members, ortholog start trims). Repointed on 2026-11-12 from
#: 9a1c2c91 by the l1-membrane-side-gate session (PR #291) at the ranker-retrain session's
#: request (PR #245), so every internal member has a mask row.
PANEL_COMMIT = "6bf341f231e0841ffea13e589c4119a909fde93f"
PANEL_DIR = "library/l1/orf-design"
#: The designed-CDS record and topology snapshot :data:`PANEL_COMMIT` pins
#: (``library_design.CDS_KEY`` and ``topology.SNAPSHOT_KEY`` on that branch).
CDS_KEY = (
    "sprint4/s4e9-orf-design-fg-panel/library-candidates/"
    "2026-11-12-ortholog-start-trims/candidate_cds.json"
)
TOPOLOGY_SNAPSHOT_KEY = (
    "accessible-surfaceome/topology/2026-11-06-ortholog-cds-policy/topology.parquet"
)
TOPOLOGY_SNAPSHOT_SHA256 = "7683d75a3e04228d95fcfb8fe9b7e9b46b38c0ea254ca42550a666f5d8bbf216"


def configure(
    *,
    panel_commit: str | None = None,
    cds_key: str | None = None,
    snapshot_key: str | None = None,
    snapshot_sha256: str | None = None,
) -> dict:
    """Override the pinned inputs for one run (a read-only check against other pins).

    Returns the pins in force. The defaults above stay the gate's record; a rerun on other
    pins writes to its own output folder (``membrane_side_gate.py --out-dir``).
    """
    global PANEL_COMMIT, CDS_KEY, TOPOLOGY_SNAPSHOT_KEY, TOPOLOGY_SNAPSHOT_SHA256
    if panel_commit:
        PANEL_COMMIT = panel_commit
    if cds_key:
        CDS_KEY = cds_key
    if snapshot_key:
        if not snapshot_sha256:
            raise SystemExit("--snapshot-key needs --snapshot-sha256")
        TOPOLOGY_SNAPSHOT_KEY, TOPOLOGY_SNAPSHOT_SHA256 = snapshot_key, snapshot_sha256
    return dict(
        panel_commit=PANEL_COMMIT,
        cds_key=CDS_KEY,
        topology_snapshot_key=TOPOLOGY_SNAPSHOT_KEY,
        topology_snapshot_sha256=TOPOLOGY_SNAPSHOT_SHA256,
    )


#: S5E13 v4 picks, after the proteoform-numbering fix (pull request 180).
S5E13_COMMIT = "c7be01a419af8b4cfa3d1570e7a6c9c2b3500346"
S5E13_DIR = "library/l1/tag-site-benchmarking"
#: The three S5E13 files the gate reads, vendored from :data:`S5E13_COMMIT` with their
#: SHA-256 in ``SOURCE.json``. That commit is on PR #180's branch, which closed unmerged, and
#: main's ``tag_sites_v4_members.tsv`` is the older 5ac64d0e table without the
#: proteoform-numbering fix (it differs in 15 picks, including ASIC1 465, SCNN1D 142 and
#: TMEM181 431), so a re-pin to a main commit would change calls. Vendored 2026-11-14 by the
#: l1-membrane-side-gate session at the L1 coordinator session's request.
S5E13_VENDORED = HERE / "ranker-c7be01a4"


def s5e13_tsv(name: str) -> pd.DataFrame:
    """One vendored S5E13 table, refused unless its SHA-256 matches ``SOURCE.json``."""
    src = json.loads((S5E13_VENDORED / "SOURCE.json").read_text())
    path = S5E13_VENDORED / name
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != src["sha256"][name]:
        raise SystemExit(f"{path}: SHA-256 {digest} differs from SOURCE.json")
    return pd.read_csv(path, sep="\t")


COLABFOLD_PREFIX = "sprint4/s4e9-orf-design-fg-panel/predicted-structures/colabfold/"
TAG_API = "https://api.deliverome.org/surfaceome/v1/tag-sites/{}"
UNIPROT_FASTA = "https://rest.uniprot.org/uniprotkb/{}.fasta"


def git_show(commit: str, path: str) -> str:
    from deliverome_analysis.core.provenance import read_pinned_text

    return read_pinned_text(commit, path, repo=REPO)


def git_tsv(commit: str, path: str, **kw) -> pd.DataFrame:
    return pd.read_csv(io.StringIO(git_show(commit, path)), sep="\t", **kw)


# --- network ------------------------------------------------------------------------


def fetch(
    url: str, dest: Path, *, retries: int = 8, missing: tuple[int, ...] = (404,)
) -> Path | None:
    """GET ``url`` into ``dest`` once; a cached miss is an empty file. Honors Retry-After.

    A status in ``missing`` means "no such entry". TmAlphaFold (like PDBTM) answers 500
    rather than 404 for an accession it does not hold, so its callers pass (404, 500) and
    a 500 is accepted as a miss after two attempts.
    """
    if dest.is_file():
        return dest if dest.stat().st_size else None
    dest.parent.mkdir(parents=True, exist_ok=True)
    for attempt in range(retries):
        req = urllib.request.Request(url, headers={"User-Agent": UA})
        try:
            with urllib.request.urlopen(req, timeout=120) as r:
                body = r.read()
            dest.write_bytes(body)
            return dest
        except urllib.error.HTTPError as e:
            if (e.code == 404 and 404 in missing) or (e.code in missing and attempt >= 1):
                dest.write_bytes(b"")
                return None
            print(f"  {url}: HTTP {e.code}, retry {attempt + 1}", file=sys.stderr)
            wait = e.headers.get("Retry-After", "") if e.headers else ""
            time.sleep(float(wait) if wait.replace(".", "", 1).isdigit() else min(2**attempt, 60))
        except (urllib.error.URLError, TimeoutError, ConnectionError) as e:
            print(f"  {url}: {e!r}, retry {attempt + 1}", file=sys.stderr)
            time.sleep(min(2**attempt, 60))
    raise RuntimeError(f"{url} still failing after {retries} attempts")


# --- members ------------------------------------------------------------------------


def use_panel_topology() -> None:
    """Point :mod:`deliverome_analysis.surfaceome.topology` at the panel's snapshot, checking its hash."""
    from deliverome_analysis.core.s3 import external

    path = external.download(TOPOLOGY_SNAPSHOT_KEY)
    digest = hashlib.sha256(Path(path).read_bytes()).hexdigest()
    if digest != TOPOLOGY_SNAPSHOT_SHA256:
        raise SystemExit(f"{TOPOLOGY_SNAPSHOT_KEY}: sha256 {digest} != pinned")
    topology.SNAPSHOT_KEY = TOPOLOGY_SNAPSHOT_KEY
    topology.load_snapshot.cache_clear()


def internal_members(extra: set[str] | frozenset[str] = frozenset()) -> pd.DataFrame:
    """The adopted panel's internal-rung members with designed protein and v1 topology.

    ``extra`` adds members of other rungs by ``candidate_id`` (built loop arms exist for
    two GPCRs whose rung is N-terminal); their ``placement`` column says so.

    One row per member: ``candidate_id``, ``candidate_set``, ``species``, ``gene_symbol``,
    ``protein`` (designed CDS translated), ``topology`` (DeepTMHMM v1 per-residue string),
    ``topology_version``, ``topology_cohort``, ``tool_version``, ``sequence_sha256``.
    The protein's hash must equal the tag-path table's, and the v1 row must be the one
    :func:`topology.for_sequence` returns for it.
    """
    from deliverome_analysis.core.s3 import processed

    use_panel_topology()
    tp = git_tsv(PANEL_COMMIT, f"{PANEL_DIR}/l1_tag_paths.tsv")
    man = git_tsv(PANEL_COMMIT, f"{PANEL_DIR}/candidate_manifest_tranches.tsv", dtype=str)
    man = man.set_index("candidate_id")
    with open(processed.download(CDS_KEY)) as fh:
        cds = json.load(fh)
    ov = pd.read_csv(HERE / "topology_overrides.tsv", sep="\t")
    overrides = dict(zip(ov["design_id"], ov["override"], strict=True))
    rows = []
    sel = tp[(tp["placement"] == "internal") | tp["candidate_id"].isin(extra)]
    for r in sel.itertuples():
        key = man.loc[r.candidate_id, "key"]
        aa = seqtools.translate(cds[key]["cds"]).rstrip("*")
        if topology.sequence_sha256(aa) != r.sequence_sha256:
            raise SystemExit(f"{r.candidate_id}: designed protein hash differs from tag paths")
        # The row the tag-path table classified the member from. A sequence can carry
        # several v1 rows (CACNG7's human canonical is also a cyno repick row); select by
        # (version, cohort), never by newest, and require every v1 row to agree.
        rows_v1 = topology.for_sequence(aa, tool=topology.TOOL_V1, all_rows=True)
        hit = rows_v1[
            (rows_v1["topology_version"] == r.topology_version)
            & (rows_v1["cohort"] == r.topology_cohort)
        ]
        if hit.empty:
            raise SystemExit(f"{r.candidate_id}: tag paths' v1 row not in the snapshot")
        row = hit.iloc[0].to_dict()
        if rows_v1["per_residue_topology"].nunique() != 1:
            raise SystemExit(f"{r.candidate_id}: v1 rows for one sequence disagree")
        v1 = row["per_residue_topology"]
        design_id = man.loc[r.candidate_id, "design_id"]
        # DeepTMHMM v1 stays the topology of record; an attributed override in
        # topology_overrides.tsv is a separate layer, applied here and named in
        # `topology_layer` so no output can mistake the corrected string for v1's own.
        if overrides.get(design_id) == "invert_v1_sides":
            topo = v1.translate(str.maketrans("IO", "OI"))
            layer = "deeptmhmm_v1_inverted_by_override"
        else:
            topo, layer = v1, "deeptmhmm_v1"
        # DeepTMHMM v2's own call on the termini, as an independent check of the (possibly
        # overridden) v1 sides; None when v2 has not been run on this exact protein.
        v2 = topology.for_sequence(aa, tool=topology.TOOL_V2)
        side = {"I": "cytoplasmic", "O": "extracellular"}
        v2_agree = (
            None
            if v2 is None
            else (
                v2.get("n_terminal_orientation") == side.get(topo[0])
                and v2.get("c_terminal_orientation") == side.get(topo[-1])
            )
        )
        rows.append(
            dict(
                candidate_id=r.candidate_id,
                design_id=design_id,
                candidate_set=r.candidate_set,
                placement=r.placement,
                species=r.species,
                gene_symbol=r.gene_symbol,
                protein=aa,
                topology=topo,
                topology_v1=v1,
                topology_layer=layer,
                v2_termini_agree=v2_agree,
                topology_version=row["topology_version"],
                topology_cohort=row["cohort"],
                tool_version=row["tool_version"],
                sequence_sha256=r.sequence_sha256,
            )
        )
    return pd.DataFrame(rows)


def tag_paths() -> pd.DataFrame:
    """The panel's tag-path table: every active L1 member's rung."""
    return git_tsv(PANEL_COMMIT, f"{PANEL_DIR}/l1_tag_paths.tsv")


def chosen_models() -> pd.DataFrame:
    """The ``is_chosen`` structure row per member, from the panel's structure payload."""
    pay = git_tsv(PANEL_COMMIT, f"{PANEL_DIR}/library_structure_payload.tsv")
    ch = pay[pay["is_chosen"] == 1].copy()
    if ch["candidate_id"].duplicated().any():
        raise SystemExit("more than one chosen model for a member")
    return ch.set_index("candidate_id")


def model_path(row: pd.Series) -> Path | None:
    """Local copy of a member's chosen model (AlphaFold DB file, or our ColabFold fold)."""
    from deliverome_analysis.core.s3 import processed

    d = CACHE / "models"
    if row["source"] == "colabfold":
        s3_key = str(row.get("s3_key", "") or "")
        if s3_key.startswith("processed://"):
            # Payloads since PR #283 carry the full key; older ones only the r2 path.
            key = s3_key.removeprefix("processed://")
        else:
            key = COLABFOLD_PREFIX + str(row["r2_key"]).split("library-colabfold/", 1)[1]
        if not processed.list_keys(prefix=key):
            # At 6bf341f2 the payload files the 5 mismatched-ortholog folds under
            # 2026-10-14-mismatched/, but they were published to
            # 2026-11-10-l1-ortholog-mismatched/ under the same name. Read the file where
            # it is; placement still requires its sequence to match the member's protein.
            # Implemented by the l1-membrane-side-gate session (PR #291), 2026-11-12;
            # reported to the reviewer-package session to fix in the payload.
            name = key.rsplit("/", 1)[1]
            alt = [p + name for p in EXTRA_FOLDS_PREFIXES if processed.list_keys(prefix=p + name)]
            if not alt:
                return None
            key = alt[0]
        path = d / ("colabfold_" + key.replace("/", "_"))
        if not path.is_file():
            d.mkdir(parents=True, exist_ok=True)
            processed.download(key, path)
        return path
    url = str(row["model_url"])
    return fetch(url, d / Path(url).name)


#: ColabFold folds of the 61 adopted orthologs that had no AlphaFold DB model, made on
#: 2026-11-10 (Modal spend approved by Becca that day; folded by the reviewer-package
#: session). They are not yet ``is_chosen`` rows of the structure payload, so the gate
#: uses one only for a member with no chosen row, and only when the fold's sequence is
#: identical to the member's designed protein. Decided by the l1-membrane-side-gate session
#: (PR #291), 2026-11-10: the task forbids passing an unmodelled member, and these are the
#: member's own models.
#: The second folder holds folds of the 5 adopted orthologs whose AlphaFold DB model is a
#: different protein (SLC22A8, SLC47A1, MC4R, ICAM1, TM7SF3), named by design id; the
#: reviewer-package session folded them on 2026-11-10 (PR #283) and pointed the gate at them.
EXTRA_FOLDS_PREFIXES = (
    "sprint4/s4e9-orf-design-fg-panel/predicted-structures/colabfold/2026-11-10-l1-orthologs/",
    "sprint4/s4e9-orf-design-fg-panel/predicted-structures/colabfold/"
    "2026-11-10-l1-ortholog-mismatched/",
    # TSPAN10's current (start-trimmed) protein, folded 2026-11-12 by the
    # l1-ortholog-features session (PR #298), named by member id; its AlphaFold DB model is
    # an older sequence and the payload's fold is stale. Added by the l1-membrane-side-gate
    # session (PR #291) at that session's request; used only on an identical sequence.
    "sprint5/s5e13-tag-site-benchmarking/predicted-structures/colabfold/2026-11-12-start-trims/",
)


def extra_folds() -> dict[str, Path]:
    """Local copies of the 2026-11-10 ortholog folds (both folders), keyed by the S3 key."""
    from deliverome_analysis.core.s3 import processed

    out = {}
    keys = [k for pre in EXTRA_FOLDS_PREFIXES for k in processed.list_keys(prefix=pre)]
    for key in keys:
        if key.endswith(".pdb"):
            path = CACHE / "models" / ("colabfold_" + key.replace("/", "_"))
            if not path.is_file():
                path.parent.mkdir(parents=True, exist_ok=True)
                processed.download(key, path)
            out[key] = path
    return out


def afdb_model(acc: str) -> Path | None:
    url = f"https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-model_v6.pdb"
    return fetch(url, CACHE / "models" / Path(url).name)


def uniprot_isoforms(acc: str) -> list[str]:
    """Isoform accessions UniProt lists for ``acc`` (ALTERNATIVE PRODUCTS), cached."""
    p = fetch(
        f"https://rest.uniprot.org/uniprotkb/{acc}.json?fields=cc_alternative_products",
        CACHE / "uniprot" / f"{acc}.isoforms.json",
    )
    if p is None:
        return []
    out = []
    for c in json.loads(p.read_text()).get("comments", []):
        if c.get("commentType") == "ALTERNATIVE PRODUCTS":
            for iso in c.get("isoforms", []):
                out += [i for i in iso.get("isoformIds", []) if i != f"{acc}-1"]
    return out


def uniprot_sequence(acc: str) -> str | None:
    p = fetch(UNIPROT_FASTA.format(acc), CACHE / "uniprot" / f"{acc}.fasta")
    if p is None:
        return None
    return "".join(line.strip() for line in p.read_text().splitlines()[1:])


# --- picks --------------------------------------------------------------------------


def v4_picks() -> pd.DataFrame:
    """S5E13 v4 member picks with their member id. Every v4 pick is a recommended order."""
    v = s5e13_tsv("tag_sites_v4_members.tsv")
    t = s5e13_tsv("l1_member_targets.tsv")
    # S5E13 keys members by their design id (l1-pNNNN); carry it to today's member id.
    man = git_tsv(PANEL_COMMIT, f"{PANEL_DIR}/candidate_manifest_tranches.tsv", dtype=str)
    design2member = dict(zip(man["design_id"], man["candidate_id"], strict=True))
    acc2design = dict(zip(t["uniprot_acc"], t["candidate_id"], strict=True))
    v["design_id"] = v["uniprot_acc"].map(acc2design)
    v["candidate_id"] = v["design_id"].map(design2member)
    # The protein S5E13 designed each member's sites on. A member whose designed protein has
    # changed since (a reference fix or a start trim) has picks numbered on another protein.
    basis = s5e13_tsv("l1_internal_members.tsv")
    v["s5e13_sequence_sha256"] = v["design_id"].map(
        dict(zip(basis["candidate_id"], basis["sequence_sha256"], strict=True))
    )
    if v["candidate_id"].isna().any():
        raise SystemExit("v4 pick with no member")
    if (v["rank"] > v["recommended_sites_to_order"]).any():
        raise SystemExit("v4 table now carries unrecommended ranks; revisit 'shipped'")
    return v


def loop_arms() -> pd.DataFrame:
    """Built panel loop arms (hand-selected sites), on the designed protein of the gene."""
    return pd.read_csv(EXPERIMENT / "constructs_loop.tsv", sep="\t")


def api_sites(gene: str) -> list[dict]:
    """Deterministic internal-lane sites the tag-site API serves for ``gene`` (human)."""
    p = fetch(TAG_API.format(gene), CACHE / "api" / f"{gene}.json")
    if p is None:
        return []
    d = json.loads(p.read_text())
    out = []
    for s in d.get("sites", []):
        # The snorkel lane is a C-terminal intracellular fallback by design.
        if (
            s.get("provenance") == "deterministic_computed"
            and s.get("site_kind") == "internal"
            and s.get("insert_after_residue") is not None
        ):
            out.append(s | {"api_uniprot_acc": d.get("uniprot_acc")})
    return out


def api_fetch_date() -> str:
    """Date the cached API responses were fetched (oldest file), for provenance."""
    files = list((CACHE / "api").glob("*.json"))
    if not files:
        return ""
    return time.strftime("%Y-%m-%d", time.localtime(min(f.stat().st_mtime for f in files)))


# --- TmAlphaFold (TMDET on AlphaFold DB models), the external placement benchmark ------


def tmalphafold(acc: str) -> tuple[dict, Path] | None:
    """TMDET membrane annotation and the TMDET-transformed AlphaFold model for ``acc``.

    ``None`` when TmAlphaFold has no entry. The transformed model has the membrane normal
    along z and the mid-plane at z = 0 (``MEMBRANE/NORMAL/Z`` is the half-thickness).
    """
    j = fetch(
        f"https://tmalphafold.ttk.hu/api/tmdet/{acc}.json",
        CACHE / "tmaf" / f"{acc}.json",
        missing=(404, 500),
    )
    if j is None:
        return None
    meta = json.loads(j.read_text())
    if "MEMBRANE" not in meta:
        return None
    pdb = CACHE / "tmaf" / f"{acc}.tr.pdb"
    if not pdb.is_file():
        page = fetch(
            f"https://tmalphafold.ttk.hu/entry/{acc}",
            CACHE / "tmaf" / f"{acc}.html",
            missing=(404, 500),
        )
        if page is None:
            return None
        m = re.search(rf'href="(/downloads/[0-9a-f/]+/{re.escape(acc)}/trpdb)"', page.read_text())
        if not m:
            return None
        gz = fetch(
            "https://tmalphafold.ttk.hu" + m.group(1),
            CACHE / "tmaf" / f"{acc}.trpdb.gz",
            missing=(404, 500),
        )
        if gz is None:
            return None
        pdb.write_bytes(gzip.decompress(gz.read_bytes()))
    return meta, pdb
