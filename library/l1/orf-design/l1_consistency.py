"""Does every built L1 arm agree with itself, with the manifest, and with what was folded?

One row per built arm across the five rung tables, with three checks:

1. ``dna``: the arm's insert, translated from the Kozak to the first stop, hashes to the row's
   ``protein_sha256``. The DNA we order encodes the protein the row claims.
2. ``layout``: that protein equals the member's current designed protein (the manifest's, by the
   SHA-256 in ``l1_tag_paths.tsv``; sequence from the committed ``tag_refold_natives.tsv``) with the 51-aa cassette at the recorded site,
   ``native[:n] + CASSETTE_AA + native[n:]``. Where a design changes more than an insertion (an Igk
   leader, a snorkel module), the check is structural instead: the cassette occurs once, the
   native C-terminal tail is intact (Igk) or the whole native protein leads (snorkel). ``n/a``
   means the layout cannot be recomputed from the row.
3. ``fold``: the tagged protein was folded by l1-tag-refold (PR #301, ``tag_refold_structures.tsv``)
   on the member's current native sequence. ``not_folded`` is reported, not failed; a fold of a
   *different* tagged protein, or of a stale native, is a mismatch.
4. ``verdict``: l1-tag-refold recorded a verdict for the tagged protein
   (``tag_refold_verdicts.tsv``, matched on ``protein_sha256``, ``mature_sha256`` or
   ``refold_tagged_sha256`` (Igk arms folded as mature chains), any
   design_id). A built arm with no recorded verdict is a mismatch (coordinator, 2026-11-13).
5. ``disposition``: a built arm whose verdict is ``fail`` has a recorded disposition: multipass in
   ``multipass_refold_dispositions.tsv``; type III in its ``refold_disposition`` column; N after SP
   in ``nterm_decisions.tsv``'s ``refold_rule``; a C-terminal-table arm in an explicit "keep"
   decision in ``l1_cterm_decisions.tsv``, or, for a snorkel, a ranker hold in s5e13
   ``l1_internal_decisions.tsv`` (the arm is superseded); an internal arm in that file's
   ``refold_verdict``/``site_history``. A failing arm with none is a mismatch (coordinator,
   2026-11-14).

Members join on ``design_id``; the current native sequence is the one in ``l1_tag_paths.tsv``.
Each mismatch names the rung's owning PR so it can be sent back. Report only: writes
``l1_consistency.tsv`` and ``l1_consistency_summary.json``; ``--check`` exits 1 on any mismatch.

Implemented by Claude Code session l1-integration (draft L1 integration, PR #299), 2026-11-12,
for the coordinator's overnight end-to-end review.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent

RUNGS = {
    # table stem: (owning PR, id column the FASTA header starts with)
    "constructs_l1_cterm": ("#283 (from #284)", "candidate_id"),
    "constructs_l1_nterm": ("#254", "candidate_id"),
    "constructs_l1_multipass": ("#280", "construct_id"),
    "constructs_l1_type3": ("#282", "candidate_id"),
    "constructs_l1_internal": ("#299 (sites #245)", "construct_id"),
}
MAPPING = HERE / "l1_member_renumbering.tsv"
TAG_PATHS = HERE / "l1_tag_paths.tsv"
FOLDS = HERE / "tag_refold_structures.tsv"
VERDICTS = HERE / "tag_refold_verdicts.tsv"
NATIVES = HERE / "tag_refold_natives.tsv"
GPI_SET = HERE / "gpi_set.tsv"
OUT = HERE / "l1_consistency.tsv"
SUMMARY = HERE / "l1_consistency_summary.json"

#: Minimum upstream clearance (aa) from a GPI omega site (positions omega-11..omega-1 form the
#: polar linker / GPI-transamidase recognition motif; matches rescue_round.GPI_OMEGA_UPSTREAM_CLEARANCE).
GPI_OMEGA_MIN_CLEARANCE = 12
VALID_INTERNAL_VETO_CLASSES = frozenset({"none", "soft"})
VALID_REFOLD_DISPOSITIONS = frozenset({"pass", "fail_flagged", "pending"})


def read(path: Path) -> list[dict]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _index_unique(rows: list[dict], key: str, label: str) -> dict[str, dict]:
    """Index ``rows`` by ``key``, failing loudly on missing/empty or duplicate keys."""
    out: dict[str, dict] = {}
    for row in rows:
        k = (row.get(key) or "").strip()
        if not k:
            raise ValueError(f"{label}: missing/empty {key!r}")
        if k in out:
            raise ValueError(f"{label}: duplicate {key} {k!r}")
        out[k] = row
    if len(out) != len(rows):
        raise ValueError(f"{label}: indexed {len(out)} of {len(rows)} rows")
    return out


def _parse_int_strict(raw: object, field: str, ctx: str) -> int:
    """Parse an integer residue/index strictly, rejecting empty, NaN, Inf, or non-integer values."""
    text = ("" if raw is None else str(raw)).strip()
    if not text:
        raise ValueError(f"{ctx}: missing/empty {field}")
    try:
        val = float(text)
    except ValueError as exc:
        raise ValueError(f"{ctx}: non-numeric {field} {raw!r}") from exc
    if not math.isfinite(val) or val != int(val):
        raise ValueError(f"{ctx}: invalid {field} {raw!r}")
    return int(val)


def sha(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


def fasta(path: Path) -> dict[str, str]:
    """Insert DNA keyed by the first token of the header (split on '|' or whitespace)."""
    from deliverome_analysis import seqtools

    out: dict[str, str] = {}
    for header, seq in seqtools.read_fasta(path):
        ident = header.split("|")[0].split()[0]
        if not ident:
            raise ValueError(f"{path.name}: empty FASTA identifier in header {header!r}")
        if ident in out:
            raise ValueError(f"{path.name}: duplicate FASTA identifier {ident!r}")
        out[ident] = seq.upper()
    return out


def translate_insert(insert: str) -> str:
    """Protein encoded from the Kozak ATG to the first in-frame stop."""
    from deliverome_analysis import plasmid_qc as qc

    try:
        return qc.translate_insert(insert)
    except ValueError as exc:
        raise AssertionError(str(exc)) from exc


class NativesUnavailable(RuntimeError):
    """A member that needs its native protein is not in the committed table, and S3 is unreachable."""


def current_natives(needed: set[str] | None = None) -> dict[str, str]:
    """Each L1 member's current designed protein, keyed by ``candidate_id``.

    Read from committed tables: ``l1_tag_paths.tsv`` gives the SHA-256 of each member's current
    designed protein, and l1-tag-refold's ``tag_refold_natives.tsv`` maps that SHA-256 to the
    sequence (each entry is re-hashed, so a sequence is used only if it is exactly that protein).
    Only a member missing from that table falls back to the designed-CDS record on S3, and only if
    it is in ``needed``; with no credentials that raises :class:`NativesUnavailable` rather than
    silently skipping the member.
    """
    natives = {}
    if NATIVES.exists():
        natives = {r["native_sha256"]: r["native_seq"] for r in read(NATIVES)}
    out, gaps = {}, {}
    for tp in read(TAG_PATHS):
        seq = natives.get(tp["sequence_sha256"])
        if seq and sha(seq) == tp["sequence_sha256"]:
            out[tp["candidate_id"]] = seq
        else:
            gaps[tp["candidate_id"]] = tp["sequence_sha256"]
    want = {c: s for c, s in gaps.items() if needed is None or c in needed}
    if want:
        from deliverome_analysis import library_design as ld
        from deliverome_analysis import s3

        if not s3.credentials_available():
            raise NativesUnavailable(
                f"{len(want)} member(s) need a native protein that is not committed: "
                + ", ".join(sorted(want)[:5])
            )
        prot = ld.designed_proteins(HERE / "candidate_manifest_tranches.tsv", library=ld.LIBRARY_L1)
        for cid, s in want.items():
            seq = prot.get(cid, {}).get("sequence", "")
            if seq and sha(seq) == s:
                out[cid] = seq
    return out


S5E13 = HERE.parents[1] / "l1" / "tag-site-benchmarking"
DISPOSITIONS = {
    "multipass": HERE / "multipass_refold_dispositions.tsv",
    "nterm": HERE / "reviewer-package" / "nterm" / "nterm_decisions.tsv",
    "cterm": HERE / "l1_cterm_decisions.tsv",
    "internal": S5E13 / "l1_internal_decisions.tsv",
    #: Per-arm decisions (agreed 2026-11-14), when the ranker has written them.
    "internal_arms": S5E13 / "l1_internal_arm_decisions.tsv",
}


def load_dispositions() -> dict[str, dict]:
    """Where each rung records what was done about a failing refold, keyed for lookup."""
    import re

    def rows(key):
        p = DISPOSITIONS[key]
        return read(p) if p.exists() else []

    cterm_keep: dict[str, list[str]] = {}
    for r in rows("cterm"):
        if (r.get("decision") or "").strip().lower().startswith("keep"):
            for did in re.findall(r"\bl1-p[0-9]{4}\b", r.get("scope", "")):
                cterm_keep.setdefault(did, []).append(r["decision"])
    mp_rows = rows("multipass")
    _index_unique(mp_rows, "construct_id", "multipass_refold_dispositions.tsv")
    nterm_rows = _index_unique(rows("nterm"), "design_id", "nterm_decisions.tsv")
    internal_rows = _index_unique(rows("internal"), "design_id", "l1_internal_decisions.tsv")
    arm_rows = _index_unique(rows("internal_arms"), "arm_id", "l1_internal_arm_decisions.tsv")
    return {
        "multipass": {r["construct_id"]: r["disposition"] for r in mp_rows if r.get("disposition")},
        "nterm": {did: r.get("refold_rule", "") for did, r in nterm_rows.items()},
        "cterm_keep": cterm_keep,
        "internal": internal_rows,
        "internal_arms": arm_rows,
    }


def disposition_of(stem: str, row: dict, dispo: dict) -> str:
    """The recorded disposition of a failing built arm, or "" if the rung records none."""
    import re

    did = row.get("design_id", "")
    internal_arms = dispo.get("internal_arms") or {}
    arm = internal_arms.get(row.get("arm_id") or "")
    if arm is not None:
        # Per-arm rule (agreed 2026-11-14): a failing arm the ranker still builds carries its own
        # ``decision == build`` and a ``flag`` saying why (e.g. so a snorkel is not the member's
        # only arm). The member's row cannot dispose of one arm among several.
        flag = (arm.get("flag") or "").strip()
        return (
            f"internal arm built with flag: {flag[:80]}"
            if arm.get("decision") == "build" and flag
            else ""
        )
    if stem == "constructs_l1_internal" and internal_arms:
        # Fail closed in per-arm mode: an internal arm absent from l1_internal_arm_decisions.tsv
        # must never inherit a member-level disposition from l1_internal_decisions.tsv.
        return ""
    if stem == "constructs_l1_multipass":
        d = dispo["multipass"].get(row.get("construct_id", ""))
        return f"multipass: {d}" if d else ""
    if stem == "constructs_l1_type3":
        d = (row.get("refold_disposition") or "").strip()
        return f"type III: {d}" if d else ""
    if stem == "constructs_l1_nterm":
        d = dispo["nterm"].get(did, "").strip()
        return f"N after SP rule: {d}" if d else ""
    if stem == "constructs_l1_cterm":
        if "snorkel" in (row.get("placement") or ""):
            dec = (dispo["internal"].get(did, {}).get("decision") or "").strip()
            return f"snorkel superseded: member {dec}" if dec.startswith("hold") else ""
        keeps = dispo["cterm_keep"].get(did, [])
        if (row.get("placement") or "") != "C-terminal":
            # A keep decision covers a re-sited arm (stalk, before furin) only if it names that
            # arm's site; "keep at the C-terminus" is about the arm the re-site replaced.
            try:
                site = str(
                    _parse_int_strict(
                        row.get("cassette_after_residue"), "cassette_after_residue", did
                    )
                )
            except ValueError:
                return ""
            keeps = [k for k in keeps if re.search(rf"\b{site}\b", k)]
        return f"C-terminal decision: {keeps[0][:80]}" if keeps else ""
    if stem == "constructs_l1_internal":
        r = dispo["internal"].get(did, {})
        d = (r.get("refold_verdict") or "").strip()
        return f"internal decision: {d}" if d and d != "pass" else ""
    return ""


def armed_members() -> set[str]:
    """candidate_ids (current numbering) of members with at least one built arm in any rung table."""
    cid_of = {
        did: r["candidate_id"]
        for did, r in _index_unique(read(MAPPING), "design_id", MAPPING.name).items()
    }
    out = set()
    for stem in RUNGS:
        tsv = HERE / f"{stem}.tsv"
        if tsv.exists():
            out |= {
                cid_of[r["design_id"]]
                for r in read(tsv)
                if (r.get("insert_nt") or "").strip() and r.get("design_id") in cid_of
            }
    return out


def layout_check(row: dict, protein: str, native: str) -> str:
    from build_constructs import CASSETTE_AA

    design = " ".join(
        (row.get(c) or "") for c in ("placement", "design", "arm", "igk_replaces_residues")
    ).lower()
    if protein.count(CASSETTE_AA) != 1:
        return f"mismatch: cassette occurs {protein.count(CASSETTE_AA)} times"
    if "snorkel" in design:
        return "ok (snorkel: native leads)" if protein.startswith(native) else "mismatch: snorkel"
    if "igk" in design:
        tail = native[-30:]
        return "ok (igk: native tail intact)" if protein.endswith(tail) else "mismatch: igk tail"
    sites = []
    raw_site = (row.get("cassette_after_residue") or "").strip()
    if raw_site:
        try:
            n = _parse_int_strict(raw_site, "cassette_after_residue", row.get("arm_id") or "")
        except ValueError as exc:
            return f"mismatch: {exc}"
        if not (0 <= n <= len(native)):
            return f"mismatch: cassette_after_residue {n} out of range 0..{len(native)}"
        sites.append(n)
    else:
        sites += [1, 0]  # residue-1 tag-alone designs: after Met1 (or before it)
    for n in sites:
        if protein == native[:n] + CASSETTE_AA + native[n:]:
            return "ok"
    if raw_site:
        return f"mismatch: not native[:{sites[0]}] + cassette + native[{sites[0]}:]"
    return "n/a"


def check_internal_arm_decisions(
    internal_rows: list[dict] | None = None,
    arm_decisions: list[dict] | None = None,
) -> dict[str, list[str]]:
    """Fail-closed validation of built internal arms against ``l1_internal_arm_decisions.tsv``.

    Asserts zero silent drops on the ``constructs_l1_internal`` -> ``l1_internal_arm_decisions``
    join, unique ``arm_id``s, no ``veto_class == "hard"`` or missing/NaN ``veto_class``, and a
    well-formed ``refold_disposition`` in :data:`VALID_REFOLD_DISPOSITIONS`.
    """
    if internal_rows is None:
        tsv = HERE / "constructs_l1_internal.tsv"
        internal_rows = read(tsv) if tsv.exists() else []
    built = [r for r in internal_rows if (r.get("insert_nt") or "").strip()]
    if not built:
        return {}
    if arm_decisions is None:
        p = DISPOSITIONS["internal_arms"]
        if not p.exists():
            raise FileNotFoundError(f"missing required internal arm decisions table: {p}")
        arm_decisions = read(p)

    dec_by_id = _index_unique(arm_decisions, "arm_id", "l1_internal_arm_decisions.tsv")
    errors: dict[str, list[str]] = defaultdict(list)
    seen_built: set[str] = set()
    matched = 0
    for row in built:
        arm_id = (row.get("arm_id") or "").strip()
        key = arm_id or row.get("construct_id") or row.get("design_id") or "<missing>"
        if not arm_id:
            errors[key].append("mismatch: built internal arm missing arm_id")
            continue
        if arm_id in seen_built:
            errors[arm_id].append(f"mismatch: duplicate arm_id {arm_id} in constructs_l1_internal")
        seen_built.add(arm_id)
        dec = dec_by_id.get(arm_id)
        if dec is None:
            errors[arm_id].append(
                f"mismatch: silent join drop — {arm_id} missing from l1_internal_arm_decisions.tsv"
            )
            continue
        matched += 1
        if (dec.get("decision") or "").strip() != "build":
            errors[arm_id].append(
                f"mismatch: built internal arm has decision={dec.get('decision')!r}"
            )
        if (dec.get("arm_kind") or "").strip() != "loop":
            errors[arm_id].append(
                f"mismatch: built internal arm has arm_kind={dec.get('arm_kind')!r}"
            )
        veto = (dec.get("veto_class") or "").strip().lower()
        if veto == "hard":
            errors[arm_id].append("mismatch: built internal arm has veto_class='hard'")
        elif veto not in VALID_INTERNAL_VETO_CLASSES:
            errors[arm_id].append(
                f"mismatch: built internal arm has missing/invalid veto_class={dec.get('veto_class')!r}"
            )
        refold_disp = (dec.get("refold_disposition") or "").strip()
        if refold_disp not in VALID_REFOLD_DISPOSITIONS:
            errors[arm_id].append(
                f"mismatch: malformed refold_disposition={dec.get('refold_disposition')!r}"
            )
    if matched != len(built) and not errors:
        raise AssertionError(f"internal join matched {matched} of {len(built)} built arms")
    return dict(errors)


def check_internal_biology(
    internal_rows: list[dict] | None = None,
    gpi_rows: list[dict] | None = None,
    *,
    min_clearance: int = GPI_OMEGA_MIN_CLEARANCE,
) -> dict:
    """Audit built internal arms for GPI omega-site clearance (fail-closed on NaN/missing omega).

    A built internal arm on a GPI-anchored gene must sit at least ``min_clearance`` aa upstream
    of ``omega_site_reference`` (positions ``omega-(min_clearance-1) .. omega-1`` form the
    GPI-transamidase recognition motif, and sites ``>= omega`` are cleaved off with the
    propeptide).
    """
    if internal_rows is None:
        tsv = HERE / "constructs_l1_internal.tsv"
        internal_rows = read(tsv) if tsv.exists() else []
    if gpi_rows is None:
        if not GPI_SET.exists():
            raise FileNotFoundError(f"missing required GPI table: {GPI_SET}")
        gpi_rows = read(GPI_SET)

    anchored = [r for r in gpi_rows if (r.get("gpi_anchored") or "").strip() == "1"]
    gpi_by_gene = _index_unique(anchored, "gene_symbol", "gpi_set.tsv")
    checked = 0
    violations: list[dict] = []
    for row in internal_rows:
        if not (row.get("insert_nt") or "").strip():
            continue
        gene = (row.get("gene_symbol") or "").strip()
        g = gpi_by_gene.get(gene)
        if g is None:
            continue
        checked += 1
        arm_id = (row.get("arm_id") or row.get("construct_id") or "").strip()
        did = (row.get("design_id") or "").strip()
        site = _parse_int_strict(
            row.get("cassette_after_residue"), "cassette_after_residue", arm_id or did
        )
        omega_ref_raw = (g.get("omega_site_reference") or "").strip()
        propep_raw = (g.get("propeptide_aa") or "").strip()
        if omega_ref_raw:
            omega = _parse_int_strict(omega_ref_raw, "omega_site_reference", gene)
        elif propep_raw:
            native_aa = _parse_int_strict(row.get("native_aa"), "native_aa", arm_id or did)
            omega = native_aa - _parse_int_strict(propep_raw, "propeptide_aa", gene)
        else:
            violations.append(
                {
                    "arm_id": arm_id,
                    "design_id": did,
                    "gene_symbol": gene,
                    "site": site,
                    "omega_site_reference": None,
                    "clearance_aa": None,
                    "reason": "GPI gene missing omega_site_reference and propeptide_aa",
                }
            )
            continue
        clearance = omega - site
        if clearance < min_clearance:
            violations.append(
                {
                    "arm_id": arm_id,
                    "design_id": did,
                    "gene_symbol": gene,
                    "site": site,
                    "omega_site_reference": omega,
                    "clearance_aa": clearance,
                    "reason": (
                        f"site {site} is past omega ({omega})"
                        if clearance <= 0
                        else f"site {site} is {clearance} aa upstream of omega ({omega}; min {min_clearance})"
                    ),
                }
            )
    return {
        "min_clearance_aa": min_clearance,
        "gpi_arms_checked": checked,
        "gpi_omega_proximal_count": len(violations),
        "gpi_omega_proximal_arms": violations,
    }


def _check_arm_row(
    stem: str,
    owner: str,
    arm_key: str,
    row: dict,
    m: dict | None,
    dna: dict[str, str],
    folds: dict[str, list[dict]],
    superseded: dict[str, list[dict]],
    by_sha: dict[str, list[dict]],
    verdicts: dict[str, str],
    by_pdb: dict[str, str],
    dispo: dict[str, dict],
    extra_problems: list[str],
) -> dict:
    did = row["design_id"]
    rec = {
        "rung_table": stem,
        "owner": owner,
        "arm_id": arm_key,
        "design_id": did,
        "gene_symbol": row.get("gene_symbol", ""),
        "placement_now": m["placement"] if m else "(not an active L1 member)",
        "dna": "",
        "layout": "",
        "fold": "",
        "verdict": "",
        "disposition": "",
        "mismatch": "",
    }
    problems = list(extra_problems)
    insert = dna.get(arm_key)
    prot = ""
    if insert is None:
        rec["dna"] = "mismatch: no insert in FASTA"
    else:
        try:
            prot = translate_insert(insert)
        except AssertionError as exc:
            prot, rec["dna"] = "", f"mismatch: {exc}"
        if prot:
            rec["dna"] = (
                "ok"
                if sha(prot) == row.get("protein_sha256")
                else "mismatch: translated insert != protein_sha256"
            )
    if m is None:
        problems.append("arm for a member not active in L1")
    elif not m["native"]:
        rec["layout"] = "mismatch: current native protein not found"
    elif insert is not None and rec["dna"] == "ok":
        rec["layout"] = layout_check(row, prot, m["native"])
    shas = {row.get(c) for c in ("protein_sha256", "mature_sha256", "refold_tagged_sha256")} - {
        None,
        "",
    }
    fs = folds.get(did, [])
    hit = [f for f in fs if f.get("tagged_sha256") in shas]
    # The refold tables are append-only, so a protein refolded on its member's new native
    # (a MANE switch that left the tagged protein unchanged, MUC15 2026-11-22) keeps its
    # earlier row: a fold on the current native is enough, and its verdict is the one read.
    current = [
        f
        for f in hit
        if not (f.get("native_sha256") and m and f["native_sha256"] != m["native_sha"])
    ]
    if hit:
        rec["fold"] = "folded" if current else "mismatch: folded on a stale native"
    elif other := [f for s in shas for f in by_sha.get(s, [])]:
        # The same protein folded under another design_id, e.g. #312's igk_tag arms,
        # whose mature protein is their retired rescue twin's (l1-p2361..p2385).
        ids = sorted({f["design_id"] for f in other})
        rec["fold"] = f"folded (same protein folded as {','.join(ids)})"
    elif any(f.get("tagged_sha256") in shas for f in superseded.get(did, [])):
        rec["fold"] = "not_folded (only a superseded fold of this protein)"
    elif fs:
        rec["fold"] = "not_folded (folds exist for other tagged proteins of this member)"
    else:
        rec["fold"] = "not_folded"
    v = next(
        (
            verdicts[s]
            for s in (
                row.get(c) for c in ("protein_sha256", "mature_sha256", "refold_tagged_sha256")
            )
            if s in verdicts
        ),
        None,
    )
    if cv := next((by_pdb[f["pdb_key"]] for f in current if f.get("pdb_key") in by_pdb), None):
        v = cv
    rec["verdict"] = v or "mismatch: no recorded refold verdict"
    if v == "fail":
        rec["disposition"] = (
            disposition_of(stem, row, dispo)
            or "mismatch: failing arm without a recorded disposition"
        )
    problems += [
        val
        for val in (rec["dna"], rec["layout"], rec["verdict"], rec["disposition"])
        if val.startswith("mismatch")
    ]
    if rec["fold"].startswith("mismatch"):
        problems.append(rec["fold"])
    rec["mismatch"] = "; ".join(problems)
    return rec


def build(*, strict_biology: bool = True) -> tuple[list[dict], dict]:
    natives = current_natives(needed=armed_members())
    mapping_rows = _index_unique(read(MAPPING), "candidate_id", MAPPING.name)
    design_of = {cid: r["design_id"] for cid, r in mapping_rows.items()}
    tag_paths = _index_unique(read(TAG_PATHS), "candidate_id", TAG_PATHS.name)
    unmapped_tp = sorted(cid for cid in tag_paths if cid not in design_of)
    if unmapped_tp:
        raise ValueError(
            f"{TAG_PATHS.name}: candidate_ids missing from {MAPPING.name}: {unmapped_tp[:5]}"
        )
    member = {
        design_of[cid]: {
            "candidate_id": cid,
            "placement": tp["placement"],
            "native": natives.get(cid, ""),
            "native_sha": tp["sequence_sha256"],
        }
        for cid, tp in tag_paths.items()
    }
    # Only live tagged folds count. #301 keeps retired folds as kind "superseded" (no pLDDT, no
    # verdict); counting those hid 19 unfolded C-terminal snorkel arms until l1-reviewer-package's
    # fold_order_check (#313) caught it, 2026-11-13.
    folds: dict[str, list[dict]] = defaultdict(list)
    superseded: dict[str, list[dict]] = defaultdict(list)
    by_sha: dict[str, list[dict]] = defaultdict(list)
    for f in read(FOLDS):
        kind = f.get("kind", "tagged")
        if kind == "tagged":
            folds[f["design_id"]].append(f)
            by_sha[f["tagged_sha256"]].append(f)
        elif kind == "superseded":
            superseded[f["design_id"]].append(f)

    verdicts: dict[str, str] = {}
    by_pdb: dict[str, str] = {}
    for v in read(VERDICTS):
        verdicts.setdefault(v["tagged_sha256"], v["verdict"])
        if v.get("pdb_key"):
            by_pdb.setdefault(v["pdb_key"], v["verdict"])

    dispo = load_dispositions()
    internal_tsv = HERE / "constructs_l1_internal.tsv"
    internal_rows = read(internal_tsv) if internal_tsv.exists() else []
    internal_errors = check_internal_arm_decisions(
        internal_rows,
        list(dispo["internal_arms"].values()) if DISPOSITIONS["internal_arms"].exists() else None,
    )
    biology_audit = check_internal_biology(internal_rows)
    bio_by_arm: dict[str, list[str]] = defaultdict(list)
    if strict_biology:
        for hit in biology_audit["gpi_omega_proximal_arms"]:
            bio_by_arm[hit["arm_id"]].append(f"mismatch: GPI omega-proximal ({hit['reason']})")

    out = []
    seen_arms: set[str] = set()
    for stem, (owner, key) in RUNGS.items():
        tsv, fa = HERE / f"{stem}.tsv", HERE / f"{stem}.fasta"
        if not tsv.exists():
            continue
        dna = fasta(fa)
        rows_for_stem = internal_rows if stem == "constructs_l1_internal" else read(tsv)
        for row in rows_for_stem:
            if not (row.get("insert_nt") or "").strip():
                continue  # a decision without DNA (e.g. "moved")
            # Tables that carry an ``arm_id`` key on it (FASTA headers lead with it), so two arms
            # of one member never collide; others keep their legacy key.
            arm_key = row.get("arm_id") or row[key]
            extra: list[str] = []
            if arm_key in seen_arms:
                extra.append(f"mismatch: duplicate arm_id {arm_key}")
            seen_arms.add(arm_key)
            if stem == "constructs_l1_internal":
                extra.extend(internal_errors.get(arm_key, []))
                extra.extend(bio_by_arm.get(arm_key, []))
            rec = _check_arm_row(
                stem,
                owner,
                arm_key,
                row,
                member.get(row["design_id"]),
                dna,
                folds,
                superseded,
                by_sha,
                verdicts,
                by_pdb,
                dispo,
                extra,
            )
            out.append(rec)

    return out, summarize(out, biology_audit=biology_audit)


def summarize(out: list[dict], biology_audit: dict | None = None) -> dict:
    """Counts by check, owner and table for ``l1_consistency_summary.json``."""
    bad = [r for r in out if r["mismatch"]]
    summary = {
        "arms_checked": len(out),
        "mismatches": len(bad),
        "mismatches_by_owner": dict(Counter(r["owner"] for r in bad)),
        "dna": dict(Counter(r["dna"].split(":")[0] for r in out)),
        "layout": dict(Counter(r["layout"].split(" (")[0].split(":")[0] or "skipped" for r in out)),
        "fold": dict(Counter(r["fold"].split(" (")[0].split(":")[0] for r in out)),
        "verdict": dict(Counter(r["verdict"].split(":")[0] for r in out)),
        "failing_arms_disposition": dict(
            Counter(r["disposition"].split(":")[0] for r in out if r["verdict"] == "fail")
        ),
        "by_table": {s: sum(r["rung_table"] == s for r in out) for s in RUNGS},
    }
    if biology_audit is not None:
        summary["biology_audit"] = biology_audit
    return summary


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument(
        "--check", action="store_true", help="check without writing; exit 1 on mismatch"
    )
    ap.add_argument(
        "--strict-biology",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="fail on GPI omega-proximal internal arms in addition to consistency mismatches (default: enabled)",
    )
    ap.add_argument("--no-write", action="store_true")
    args = ap.parse_args(argv)
    rows, summary = build() if args.strict_biology else build(strict_biology=False)
    if not args.check and not args.no_write:
        with OUT.open("w") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(rows)
        SUMMARY.write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps(summary, indent=1))
    if args.check and summary["mismatches"]:
        print(f"{summary['mismatches']} mismatches", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
