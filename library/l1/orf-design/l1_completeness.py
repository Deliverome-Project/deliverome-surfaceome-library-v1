"""Is every active L1 member placed and armed? One row per member, across all five rungs.

The draft L1 is complete when every active member in ``library_manifest_l1.tsv`` has

1. a placement (its rung in ``l1_tag_paths.tsv``), and
2. a built tag arm on that rung, or a recorded decision that it has none.

Each rung's arms come from that rung's own builder and PR, so this script only reads them:

| Placement | Arm table | PR |
|---|---|---|
| C-terminal | ``constructs_l1_cterm.tsv`` | #284 |
| N after SP | ``constructs_l1_nterm.tsv``; held or refused picks in ``reviewer-package/nterm/nterm_site_picks.tsv`` | #254 |
| N at residue 1, multipass | ``constructs_l1_multipass.tsv`` | #280 |
| N at residue 1, single-pass (type III) | ``constructs_l1_type3.tsv`` | #282 |
| internal | ``constructs_l1_internal.tsv`` (built here from s5e13 ``l1_internal_decisions.tsv``); a non-loop decision there (snorkel, reroute) counts as held | #245 / #299 |

Members are joined on ``design_id``, never on ``l1-mNNNN``. #283 renumbers member ids, and
the rung tables were written against different numberings. When an arm table also carries
``candidate_id``, a disagreement with the current mapping is reported as ``stale_member_id``:
the arm was built for this design under an older numbering, and its owner must rekey it.

Status per member:

- ``built``: at least one arm with status decided/built, or no status column (a built table).
- ``proposed``: arms exist but every one is a proposal awaiting a call.
- ``held``: the rung records a hold or refusal for the member, and no arm.
- ``missing``: placed, but no arm and no recorded decision. This is a gap.
- ``unplaced``: active in the manifest, absent from ``l1_tag_paths.tsv``.

An arm whose table records the native sequence (``native_sha256``), or whose tagged protein can be
recomputed as ``native[:n] + cassette + native[n:]`` from ``cassette_after_residue`` and
``protein_sha256``, is also checked against the member's current designed protein. An arm built on
another sequence (``stale_sequence``) does not count, so a start fix or reference swap that
changes the protein reopens the gap until its rung rebuilds.

Report only. It changes no decision and writes ``l1_completeness.tsv`` and
``l1_completeness_summary.json``. ``--check`` exits 1 while any member is ``missing``,
``unplaced`` or ``stale_member_id``.

Implemented by Claude Code session l1-integration (draft L1 integration, PR #299), 2026-11-12,
for the coordinator's target of a complete draft L1 by the morning of 2026-11-12.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
S5E13 = HERE.parents[1] / "l1" / "tag-site-benchmarking"
#: The ranker's per-arm decisions. Once written, the order rule is in force (Becca, 2026-11-14):
#: an internal member has 1-3 loop arms, and a snorkel is never a member's only arm.
ARM_DECISIONS = S5E13 / "l1_internal_arm_decisions.tsv"
#: Rungs that carry exactly one arm per member.
SINGLE_ARM_RUNGS = (
    "N after SP",
    "N at residue 1, multipass",
    "N at residue 1, single-pass (type III)",
)

MANIFEST = HERE / "library_manifest_l1.tsv"
MAPPING = HERE / "l1_member_renumbering.tsv"
TAG_PATHS = HERE / "l1_tag_paths.tsv"
OUT = HERE / "l1_completeness.tsv"
SUMMARY = HERE / "l1_completeness_summary.json"

#: placement -> (arm tables, PR). Missing files are reported, not fatal.
ARMS: dict[str, tuple[list[Path], str]] = {
    "C-terminal": ([HERE / "constructs_l1_cterm.tsv"], "#284"),
    "N after SP": ([HERE / "constructs_l1_nterm.tsv"], "#254"),
    "N at residue 1, multipass": ([HERE / "constructs_l1_multipass.tsv"], "#280"),
    "N at residue 1, single-pass (type III)": ([HERE / "constructs_l1_type3.tsv"], "#282"),
    # Snorkel arms for internal members are C-terminal appends, built by #284's l1_snorkel.py
    # into constructs_l1_cterm.tsv; a C-terminal member is never on the internal rung, so the
    # design_id join cannot cross rungs.
    "internal": (
        [HERE / "constructs_l1_internal.tsv", HERE / "constructs_l1_cterm.tsv"],
        "#299 (loops) / #284 (snorkel)",
    ),
}
#: Records of a decision not to build (yet), by placement.
HOLDS: dict[str, list[Path]] = {
    # #284 records rung calls by a free-text scope ("l1-m1279 (l1-p2487) PLVAP human").
    "C-terminal": [HERE / "l1_cterm_decisions.tsv"],
    "N after SP": [
        HERE / "reviewer-package" / "nterm" / "nterm_decisions.tsv",  # #254's current record
        HERE / "reviewer-package" / "nterm" / "nterm_site_picks.tsv",
    ],
    "internal": [S5E13 / "l1_internal_decisions.tsv"],
}
PROPOSED_SITES: dict[str, list[Path]] = {
    "internal": [S5E13 / "l1_internal_decisions.tsv"],
}

#: Multipass refold verdicts live in their own table (#280), one row per arm that did not pass;
#: every other built multipass arm passed.
MULTIPASS_REFOLD = HERE / "multipass_refold_dispositions.tsv"
#: C-terminal refold verdicts (#284), keyed by design_id and the tagged protein's sha256.
CTERM_REFOLD = HERE / "l1_cterm_refold_verdicts.tsv"
#: l1-tag-refold's verdict per tagged protein (#301); the fallback where a rung records none.
TAG_VERDICTS = HERE / "tag_refold_verdicts.tsv"

BUILT = {"", "decided", "built", "build"}
PROPOSED = {"proposed"}


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
    return out


def by_design(rows: list[dict], *, strict: bool = False) -> dict[str, list[dict]]:
    out: dict[str, list[dict]] = defaultdict(list)
    for r in rows:
        did = (r.get("design_id") or "").strip()
        if not did and r.get("scope"):
            m = re.search(r"\b(l1-p[0-9]{4})\b", r["scope"])
            did = m.group(1) if m else ""
        if did:
            out[did].append(r)
        elif strict:
            raise ValueError(f"row missing design_id in {r.get('_table') or 'arm table'}: {r!r}")
    return out


def _check_arm_table_ids(table_name: str, rows: list[dict]) -> None:
    """Reject duplicate built ``arm_id`` or ``construct_id`` entries in an arm table."""
    seen_arm: set[str] = set()
    seen_cid: set[str] = set()
    for r in rows:
        if "insert_nt" in r and not (r.get("insert_nt") or "").strip():
            continue
        aid = (r.get("arm_id") or "").strip()
        if aid:
            if aid in seen_arm:
                raise ValueError(f"{table_name}: duplicate arm_id {aid!r}")
            seen_arm.add(aid)
        cid = (r.get("construct_id") or "").strip()
        if cid:
            if cid in seen_cid:
                raise ValueError(f"{table_name}: duplicate construct_id {cid!r}")
            seen_cid.add(cid)


def _load_arm_built(path: Path) -> set[str]:
    """Design ids with a ``decision == build`` row, failing closed on duplicate/malformed rows."""
    if not path.exists():
        return set()
    rows = read(path)
    _index_unique(rows, "arm_id", path.name)
    out: set[str] = set()
    for r in rows:
        dec = (r.get("decision") or "").strip()
        if dec not in ("build", "hold"):
            raise ValueError(f"{path.name}: {r['arm_id']} has unexpected decision {dec!r}")
        if dec == "build":
            did = (r.get("design_id") or "").strip()
            if not did:
                raise ValueError(f"{path.name}: {r['arm_id']} missing design_id")
            out.add(did)
    return out


def _multipass_verdict(arm: dict, mp_refold: dict, tag_verdicts: dict[str, dict]) -> str:
    """Resolve a multipass arm's refold verdict without silently defaulting unrecorded fails."""
    m = mp_refold.get(arm.get("construct_id", ""))
    if m:
        verdict = m["verdict"]
        if m.get("disposition", "").startswith("keep"):
            verdict += " (kept as built)"
        return verdict
    if tag_verdicts:
        hit = next(
            (
                tag_verdicts[s]
                for s in (arm.get("protein_sha256"), arm.get("mature_sha256"))
                if s and s in tag_verdicts
            ),
            None,
        )
        if hit is None:
            return "not recorded"
        if hit["verdict"] != "pass":
            return f"{hit['verdict']} (unrecorded in multipass_refold_dispositions)"
    return "pass"


def hold_reason(row: dict) -> str | None:
    """A recorded non-build decision in a hold/pick table row, or None if it is a build."""
    if row.get("arm") == "no_nterm_arm":
        return "arm=no_nterm_arm"
    if row.get("arm") in ("native_site", "igk_swap"):
        return None
    for col in ("status", "decision", "build_status", "verdict", "action"):
        if col == "decision" and "scope" in row and "hold" not in (row.get(col) or "").lower():
            continue  # a #284 decision row that is not a hold
        v = (row.get(col) or "").strip()
        if (
            v
            and v.lower() not in BUILT | {"build_now", "pick"}
            and not v.startswith("internal_loop")
        ):
            return f"{col}={v}"
    return None


def sha(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


def arm_matches_sequence(arm: dict, native: str, native_sha: str) -> bool | None:
    """Was this arm built on the member's current sequence? None when the table cannot say."""
    if "snorkel" in (arm.get("placement") or ""):
        return None  # native + TM module + cassette; not recomputable here
    if arm.get("native_sha256"):
        return arm["native_sha256"] == native_sha
    n, tagged = arm.get("cassette_after_residue"), arm.get("protein_sha256")
    if arm.get("igk_replaces_residues") or arm.get("arm", "native") not in (
        "",
        "native",
        "native_site",
    ):
        return None  # an Igk swap's protein is not native[:n] + cassette + native[n:]
    if n and tagged and native:
        from build_constructs import CASSETTE_AA

        try:
            val = float(n)
        except ValueError:
            return False
        if not math.isfinite(val) or val != int(val):
            return False
        k = int(val)
        if not (0 <= k <= len(native)):
            return False
        return sha(native[:k] + CASSETTE_AA + native[k:]) == tagged
    return None


def arm_shape_problem(placement: str, arms: list[dict]) -> str:
    """'' if a built member's arms fit its rung, else the rule they break.

    - ``snorkel_only``: every arm is a snorkel (a snorkel is never a member's only arm).
    - ``multi_arm_on_single_arm_rung``: more than one arm on a rung that carries exactly one.
    """
    if placement in SINGLE_ARM_RUNGS and len(arms) > 1:
        return "multi_arm_on_single_arm_rung"
    if arms and all("snorkel" in (a.get("placement") or "") for a in arms):
        return "snorkel_only"
    return ""


def orthologs_without_built_human(rows: list[dict]) -> set[str]:
    """design_ids of built mouse or cyno orthologs whose gene has no built human member.

    Becca, 2026-11-14: never substitute a mouse or cyno ortholog for a human member, and do not
    build orthologs for a gene with no good human proteoform. So an ortholog may be built only
    while at least one human member of the same gene is built. Genes are grouped by
    ``gene_symbol``; species by ``species``.
    """
    human_built: set[str] = set()
    for r in rows:
        if r["species"] == "human" and r["status"] == "built":
            human_built.add(r["gene_symbol"])
    return {
        r["design_id"]
        for r in rows
        if r["species"] != "human"
        and r["status"] == "built"
        and r["gene_symbol"] not in human_built
    }


def _load_core_tables() -> tuple[dict[str, str], list[dict], dict[str, dict]]:
    """Load manifest, mapping, and tag paths with strict key uniqueness and join checks."""
    mapping_rows = _index_unique(read(MAPPING), "candidate_id", MAPPING.name)
    mapping = {cid: r["design_id"] for cid, r in mapping_rows.items()}
    manifest = read(MANIFEST)
    _index_unique(manifest, "candidate_id", MANIFEST.name)
    unmapped = sorted(r["candidate_id"] for r in manifest if r["candidate_id"] not in mapping)
    if unmapped:
        raise ValueError(
            f"{MANIFEST.name}: candidate_ids missing from {MAPPING.name}: {unmapped[:5]}"
        )
    paths = _index_unique(read(TAG_PATHS), "candidate_id", TAG_PATHS.name)
    unknown_p = sorted({tp["placement"] for tp in paths.values()} - set(ARMS))
    if unknown_p:
        raise ValueError(f"{TAG_PATHS.name}: unrecognized placement(s): {unknown_p}")
    return mapping, manifest, paths


def build() -> tuple[list[dict], dict]:
    # Current designed proteins from committed tables (S3 only for an armed member they lack).
    sys.path.insert(0, str(HERE))
    import l1_consistency

    natives = l1_consistency.current_natives(needed=l1_consistency.armed_members())
    mapping, manifest, paths = _load_core_tables()

    arm_index: dict[str, dict[str, list[dict]]] = {}
    missing_tables: list[str] = []
    for placement, (tables, _) in ARMS.items():
        rows: list[dict] = []
        for t in tables:
            if t.exists():
                t_rows = read(t)
                _check_arm_table_ids(t.name, t_rows)
                rows += [{**r, "_table": t.name} for r in t_rows]
            else:
                missing_tables.append(str(t.relative_to(HERE.parents[2])))
        arm_index[placement] = by_design(rows, strict=True)
    hold_index = {
        p: by_design([r for t in ts if t.exists() for r in read(t)]) for p, ts in HOLDS.items()
    }
    site_index = {
        p: by_design(
            [
                {**r, "design_id": r.get("design_id") or ""}
                for t in ts
                if t.exists()
                for r in read(t)
            ]
        )
        for p, ts in PROPOSED_SITES.items()
    }

    arm_built = _load_arm_built(ARM_DECISIONS)
    out = []
    shape_problems: dict[str, set[str]] = defaultdict(set)
    for m in manifest:
        cid = m["candidate_id"]
        did = mapping.get(cid, "")
        tp = paths.get(cid)
        row = {
            "design_id": did,
            "candidate_id": cid,
            "gene_symbol": m.get("gene_symbol", ""),
            "candidate_set": m.get("candidate_set", ""),
            "species": m.get("species", ""),
            "placement": tp["placement"] if tp else "",
            "status": "",
            "arms": 0,
            "arm_statuses": "",
            "arm_source": "",
            "note": "",
            "decided_by": "",
        }
        if not tp:
            row["status"] = "unplaced"
            out.append(row)
            continue
        placement = tp["placement"]
        tables, pr = ARMS.get(placement, ([], ""))
        row["arm_source"] = pr
        arms = arm_index.get(placement, {}).get(did, [])
        stale = sorted(
            {a["candidate_id"] for a in arms if a.get("candidate_id") and a["candidate_id"] != cid}
        )
        native = natives.get(cid, "")
        drift = [a for a in arms if arm_matches_sequence(a, native, tp["sequence_sha256"]) is False]
        if drift:
            row["note"] = (row["note"] + "; " if row["note"] else "") + (
                f"stale_sequence: {len(drift)} arm(s) built on another sequence"
            )
            arms = [a for a in arms if a not in drift]
        # A row that records a decision but no DNA (empty insert_nt, e.g. design "moved")
        # is a decided non-arm, not a built arm.
        no_dna = [a for a in arms if "insert_nt" in a and not (a.get("insert_nt") or "").strip()]
        if no_dna and len(no_dna) == len(arms):
            why = no_dna[0].get("design") or no_dna[0].get("decision") or "no insert"
            row["decided_by"] = no_dna[0].get("decided_by", "")
            row["status"] = "held"
            row["note"] = "; ".join(x for x in (row["note"], f"decided without an arm: {why}") if x)
            out.append(row)
            continue
        arms = [a for a in arms if a not in no_dna]
        # The ranker's hold is the decision of record: an arm that a hold supersedes (e.g. a
        # failed snorkel still in constructs_l1_cterm.tsv) does not make the member built.
        ranker_hold = [
            h
            for h in hold_index.get(placement, {}).get(did, [])
            if (h.get("decision") or "").startswith("hold")
        ]
        if did in arm_built:
            # Per-arm mode: the arm rows are the decision of record. A member-level hold that
            # predates them (e.g. SLC35G2, MS4A1: arms built under the control-relative refold
            # bar) does not override arms the per-arm table builds.
            ranker_hold = []
        if arms and ranker_hold:
            row["status"] = "held"
            row["decided_by"] = ranker_hold[0].get("decided_by", "")
            tabs = sorted({a.get("_table", "") for a in arms} - {""}) or ["an arm table"]
            row["note"] = "; ".join(
                x
                for x in (
                    row["note"],
                    f"decision={ranker_hold[0]['decision']}; {len(arms)} superseded arm(s) "
                    f"still in {','.join(tabs)}, awaiting drop by the table's owner",
                )
                if x
            )
            out.append(row)
            continue
        if arms:
            statuses = Counter((a.get("status") or "").strip().lower() for a in arms)
            row["arms"] = len(arms)
            row["arm_statuses"] = ";".join(
                f"{k or 'built'}:{v}" for k, v in sorted(statuses.items())
            )
            if any(s in BUILT for s in statuses):
                row["status"] = "built"
                shape = arm_shape_problem(placement, arms)
                if shape:
                    shape_problems[shape].add(did)
                    row["note"] = "; ".join(x for x in (row["note"], shape) if x)
            elif all(s in PROPOSED for s in statuses):
                row["status"] = "proposed"
            else:
                row["status"] = "held"
        else:
            holds = [h for h in hold_index.get(placement, {}).get(did, []) if hold_reason(h)]
            sites = site_index.get(placement, {}).get(did, [])
            if holds:
                row["status"] = "held"
                row["decided_by"] = holds[0].get("decided_by", "")
                row["note"] = "; ".join(x for x in (row["note"], hold_reason(holds[0])) if x)
            elif sites:
                row["status"] = "missing"
                row["note"] = "; ".join(
                    x for x in (row["note"], f"{len(sites)} proposed site(s), no built arm") if x
                )
            else:
                row["status"] = "missing"
        if stale:
            row["note"] = (row["note"] + "; " if row["note"] else "") + (
                "stale_member_id: arm table has " + ",".join(stale)
            )
        out.append(row)

    status = Counter(r["status"] for r in out)

    # Refold verdicts per rung, over the arms that count as built (current sequence, with DNA).
    ct_refold = {}
    if CTERM_REFOLD.exists():
        for r in read(CTERM_REFOLD):
            ct_refold[(r["design_id"], r.get("tagged_sha256", ""))] = r["verdict"]
    mp_refold = (
        {r["construct_id"]: r for r in read(MULTIPASS_REFOLD)} if MULTIPASS_REFOLD.exists() else {}
    )
    member_of = {
        mapping[cid]: {"candidate_id": cid, "sha": tp["sequence_sha256"]}
        for cid, tp in paths.items()
        if cid in mapping
    }
    tag_verdicts: dict[str, dict] = {}
    if TAG_VERDICTS.exists():
        for r in read(TAG_VERDICTS):
            tag_verdicts.setdefault(r["tagged_sha256"], r)
    snorkel_fail: set[str] = set()
    refold: dict[str, Counter] = {}
    for placement, (tables, _) in ARMS.items():
        active = {
            r["design_id"]
            for r in out
            if r["placement"] == placement and r["status"] in ("built", "proposed")
        }
        c: Counter = Counter()
        for t in tables:
            if not t.exists():
                continue
            for a in read(t):
                if a.get("design_id") not in active or not (a.get("insert_nt") or "").strip():
                    continue
                cid = member_of.get(a["design_id"], {})
                native = natives.get(cid.get("candidate_id", ""), "")
                if arm_matches_sequence(a, native, cid.get("sha", "")) is False:
                    continue
                if (a.get("placement") or placement) not in (
                    placement,
                    "internal, before furin site",
                ):
                    if not (placement == "internal" and t.name == "constructs_l1_cterm.tsv"):
                        continue
                if placement == "N at residue 1, multipass" and mp_refold:
                    verdict = _multipass_verdict(a, mp_refold, tag_verdicts)
                elif t.name == "constructs_l1_cterm.tsv" and ct_refold:
                    verdict = ct_refold.get(
                        (a["design_id"], a.get("protein_sha256", "")), "not recorded"
                    )
                else:
                    verdict = (a.get("refold_verdict") or "").strip() or "not recorded"
                if verdict in ("not recorded", "not_in_table"):
                    # Fall back to l1-tag-refold's own verdict for this tagged protein (#301).
                    # Snorkels had none in their rung's table; a failed snorkel stays built,
                    # pending the ranker's redecision (next v5 site, or hold if none), as decided
                    # by the coordinator, 2026-11-13.
                    hit = next(
                        (
                            tag_verdicts[s]
                            for s in (a.get("protein_sha256"), a.get("mature_sha256"))
                            if s and s in tag_verdicts
                        ),
                        None,
                    )
                    if hit:
                        verdict = f"{hit['verdict']} (#301)"
                        if "snorkel" in (a.get("placement") or "") and hit["verdict"] == "fail":
                            verdict = "fail (#301; snorkel, ranker redecision pending)"
                            snorkel_fail.add(a["design_id"])
                c[verdict] += 1
        refold[placement] = dict(sorted(c.items()))
    for r in out:
        if r["design_id"] in snorkel_fail:
            note = "snorkel refold fail (#301); built, pending the ranker's redecision"
            r["note"] = f"{r['note']}; {note}" if r["note"] else note
    orphans = orthologs_without_built_human(out)
    for r in out:
        if r["design_id"] in orphans:
            note = "ortholog_without_built_human: no human member of this gene is built"
            r["note"] = f"{r['note']}; {note}" if r["note"] else note
    by_placement = defaultdict(Counter)
    for r in out:
        by_placement[r["placement"] or "(unplaced)"][r["status"]] += 1
    summary = {
        "active_members": len(out),
        "status": dict(sorted(status.items())),
        "by_placement": {k: dict(sorted(v.items())) for k, v in sorted(by_placement.items())},
        "stale_member_id": sum("stale_member_id" in r["note"] for r in out),
        "stale_sequence": sum("stale_sequence" in r["note"] for r in out),
        "snorkel_only_members": sorted(shape_problems["snorkel_only"]),
        "multi_arm_on_single_arm_rung": sorted(shape_problems["multi_arm_on_single_arm_rung"]),
        "arm_rule_in_force": ARM_DECISIONS.exists(),
        "orthologs_without_built_human": sorted(orphans),
        "refold_verdicts_by_placement": refold,
        "arm_tables_missing": missing_tables,
        "join_key": "design_id (via l1_member_renumbering.tsv)",
    }
    return out, summary


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--check", action="store_true", help="exit 1 while any member has a gap")
    ap.add_argument("--no-write", action="store_true")
    args = ap.parse_args(argv)
    rows, summary = build()
    if not args.no_write:
        with OUT.open("w") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(rows)
        SUMMARY.write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps(summary, indent=1))
    gaps = summary["status"].get("missing", 0) + summary["status"].get("unplaced", 0)
    shape = summary["multi_arm_on_single_arm_rung"] + (
        summary["snorkel_only_members"] if summary["arm_rule_in_force"] else []
    )
    if shape:
        print(f"arm-shape violations: {shape}", file=sys.stderr)
    orphans = summary["orthologs_without_built_human"]
    if args.check and (
        gaps
        or summary["stale_member_id"]
        or summary["stale_sequence"]
        or summary["arm_tables_missing"]
        or shape
        or orphans
    ):
        print(
            f"incomplete: {gaps} gaps, {summary['stale_member_id']} stale member ids, "
            f"{summary['stale_sequence']} stale sequences, "
            f"{len(summary['arm_tables_missing'])} missing tables, "
            f"{len(orphans)} orthologs built without a built human member {orphans}",
            file=sys.stderr,
        )
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
