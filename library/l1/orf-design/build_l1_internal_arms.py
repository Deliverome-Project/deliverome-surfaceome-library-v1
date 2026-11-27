"""Put the tag cassette at each proposed internal site of every L1 internal-rung member.

This is the internal rung of issue #2 ("assemble the tag cassette into every L1 insert"). Every
input is fixed in the record; nothing is chosen here:

* **which members**: ``placement == "internal"`` in ``l1_tag_paths.tsv`` (DeepTMHMM v1,
  matched by sequence).
* **where**: the loop ranker's decision for each member (s5e13, PR #245,
  ``l1_internal_decisions.tsv``; Becca delegated the call on 2026-11-12). A member whose
  ``decision`` starts with ``internal_loop`` gets one arm, at ``insert_after_residue``. Any
  other decision (``snorkel``, ``route_to_nterm_after_sp``) is listed in the summary with its
  reason and gets no arm here. A decision is used only when its ``sequence_sha256`` equals the
  member's in ``l1_tag_paths.tsv``, so a site called on another proteoform is never applied.
  ``decided_by`` / ``implemented_by`` / ``reason`` are carried through unchanged: the ranker
  owns the site call, and this script owns only the DNA.
* **the ORF**: the member's domesticated ORF from ``library_design.CONSTRUCTS_KEY``. It must
  translate to the member's designed protein plus a stop, and that protein must hash to the
  tag-path ``sequence_sha256``.
* **the cassette**: the 51-aa tandem on this arm's own linker encoding
  (``l1_arm_cassette_encodings.tsv``; SpyTag003 and P-ALFA-P on p003's codons). One shared encoding
  was used until 2026-11-21; Becca, 2026-11-21: "we still don't want shared linkers in the
  middle". An internal cassette shared by every arm lets a template switch inside it join the
  front of one protein to the back and barcode of another.
  The full cassette is kept, leading GGGGS included (Becca, 2026-11-11, "keep full").
* **the layout**: #254's (``build_l1_nterm_arms.py``), which is #284's Esp3I wrapper::

      CGTCTCA | Kozak | ORF[:3n] | cassette | ORF[3n:] (native stop) | CC
              | T7 | 5' handle | barcode (20 N) | Feldman 3' | 10x CS1 | AGA GAGACG

Per arm it asserts that the protein is exactly ``native[:n] + cassette + native[n:]``, that the
junction residues are the site's, and that there is one stop, at the end. A seam that completes
a BbsI/Esp3I/BsaI site is reported and the arm held, never silently recoded. Every built arm is
also checked inside the p0039 backbone.

Members on the rung with no matching site are listed in the summary (``no_site``). They are the
ranker's to fill. ``l1_completeness.py`` reports them as gaps.

**Attribution** (convention of 2026-11-10). Site calls: as carried in each row, from the ranker
sessions (PR #245, PR #298). DNA build: implemented by Claude Code session l1-integration
(draft L1 integration, PR #299), 2026-11-12. The coordinator asked for every member to have a
built arm by morning, and no session owned the internal rung's DNA. Layout and checks are #254's
and #284's, reused.

    uv run --locked python build_l1_internal_arms.py            # build and write
    uv run --locked python build_l1_internal_arms.py --check    # rebuild, compare with the record
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import tempfile
from collections import Counter
from pathlib import Path

import build_l1_nterm_arms as nterm
from build_constructs import CASSETTE, CASSETTE_AA
from diversify_arm_cassettes import PROVISIONAL_ENV, arm_cassette

from deliverome_analysis.cloning import library_design as ld
from deliverome_analysis.cloning import plasmid_qc as qc
from deliverome_analysis.cloning import seqtools

HERE = Path(__file__).resolve().parent
S5E13 = HERE.parents[1] / "l1" / "tag-site-benchmarking"

MANIFEST = HERE / "candidate_manifest_tranches.tsv"
TAG_PATHS = HERE / "l1_tag_paths.tsv"
MAPPING = HERE / "l1_member_renumbering.tsv"
DECISIONS = S5E13 / "l1_internal_decisions.tsv"
#: The ranker's per-arm decisions (agreed schema, 2026-11-14; L1-MULTI-ARM notes): one row per
#: arm, ``arm_id`` = ``l1-pNNNN-aK`` (creation order, never reused), ``arm_kind`` loop|snorkel,
#: ``decision`` build|hold. When present it replaces the one-site-per-member decision: every
#: ``loop`` row with ``decision == build`` becomes one arm here, named by its ``arm_id``.
#: Snorkel arms are built by the C-terminal rung's snorkel builder. Absent, the build is the
#: one-arm-per-member record, byte for byte.
ARM_DECISIONS = S5E13 / "l1_internal_arm_decisions.tsv"
OUT_TSV = HERE / "constructs_l1_internal.tsv"
OUT_FASTA = HERE / "constructs_l1_internal.fasta"
OUT_JSON = HERE / "constructs_l1_internal_summary.json"

BARCODE_SLOT = qc.BARCODE_SLOT
ENZYMES = list(qc.DOMESTICATION_ENZYMES)
BUILT_BY = "Claude Code session l1-integration, PR #299 (DNA build only)"
COLUMNS = [
    "construct_id",
    "design_id",
    "candidate_id",
    "gene_symbol",
    "candidate_set",
    "species",
    "decision",
    "pick_rank",
    "confidence",
    "refold_verdict",
    "status",
    "native_aa",
    "protein_aa",
    "orf_nt",
    "coding_nt",
    "insert_nt",
    "placement",
    "cassette_after_residue",
    "junction",
    "cassette_start_nt",
    "type_iis_sites",
    "snorkel_candidate",
    "site_source",
    "decided_by",
    "implemented_by",
    "reason",
    "built_by",
    "native_sha256",
    "protein_sha256",
]
#: Added in per-arm mode only, so the one-arm record keeps its exact bytes.
ARM_COLUMNS = ["arm_id", "arm_role", "construct_key", "order_verdict", "flag"]


VALID_ARM_DECISIONS = frozenset({"build", "hold"})
VALID_LOOP_VETO_CLASSES = frozenset({"none", "soft"})
VALID_REFOLD_DISPOSITIONS = frozenset({"pass", "fail_flagged", "pending"})


def read_table(path: Path) -> list[dict]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _index_unique(rows: list[dict], key: str, label: str) -> dict[str, dict]:
    out: dict[str, dict] = {}
    for row in rows:
        k = (row.get(key) or "").strip()
        if not k:
            raise ValueError(f"{label}: missing/empty {key!r}")
        if k in out:
            raise ValueError(f"{label}: duplicate {key} {k!r}")
        out[k] = row
    return out


def _parse_int_strict(raw: object, field: str, ctx: str) -> int:
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


def decisions() -> dict[str, dict]:
    """The ranker's one decision per member, keyed by design_id."""
    out = {}
    for r in read_table(DECISIONS):
        did = (r.get("design_id") or "").strip()
        assert did, f"{DECISIONS.name}: row missing design_id"
        assert did not in out, f"{did}: two decisions"
        out[did] = r
    return out


def arm_calls(path: Path | None = None) -> dict[str, list[dict]] | None:
    """The ranker's loop arms to build, per design_id, or None when there is no per-arm table.

    Checks the agreed invariants on the whole table: arm ids are unique and belong to their
    design; decisions are strictly ``build`` or ``hold``; built loop arms carry valid
    ``veto_class`` (never ``hard`` or missing/NaN) and ``refold_disposition`` when those columns
    are present; a member has at most three loop arms to build; a snorkel row is ``build`` only
    if the member also has a loop to build (a snorkel is never a member's only arm).
    """
    path = path or ARM_DECISIONS
    if not path.exists():
        return None
    rows = read_table(path)
    ids = [(r.get("arm_id") or "").strip() for r in rows]
    assert all(ids), f"{path.name}: row missing arm_id"
    assert len(ids) == len(set(ids)), "duplicate arm_id in the per-arm decisions"
    by: dict[str, list[dict]] = {}
    snorkel_build: dict[str, int] = Counter()
    for r in rows:
        did = (r.get("design_id") or "").strip()
        aid = r["arm_id"].strip()
        assert did and aid.startswith(did + "-a"), f"{aid}: not of its design {did!r}"
        dec = (r.get("decision") or "").strip()
        assert dec in VALID_ARM_DECISIONS, f"{aid}: unknown decision {r.get('decision')!r}"
        if dec != "build":
            continue
        if r["arm_kind"] == "loop":
            if "veto_class" in r:
                vc = (r.get("veto_class") or "").strip().lower()
                assert vc in VALID_LOOP_VETO_CLASSES, (
                    f"{aid}: invalid/hard veto_class {r.get('veto_class')!r} on built loop arm"
                )
            if "refold_disposition" in r:
                rd = (r.get("refold_disposition") or "").strip()
                assert rd in VALID_REFOLD_DISPOSITIONS, (
                    f"{aid}: malformed refold_disposition {r.get('refold_disposition')!r}"
                )
            by.setdefault(did, []).append(r)
        elif r["arm_kind"] == "snorkel":
            snorkel_build[did] += 1
        else:
            raise AssertionError(f"{aid}: unknown arm_kind {r['arm_kind']!r}")
    for did, arms in by.items():
        assert len(arms) <= 3, f"{did}: {len(arms)} loop arms to build (at most 3)"
    for did in snorkel_build:
        assert did in by, f"{did}: a snorkel would be the member's only arm"
    return {did: sorted(arms, key=lambda r: r["arm_id"]) for did, arms in by.items()}


def build_one(member: dict, site: dict, orf: str, protein: str, tag_fixed: dict) -> dict | str:
    """The insert for one site, or a reason it cannot be built as proposed."""
    did = member["design_id"]
    n = _parse_int_strict(site.get("insert_after_residue"), "insert_after_residue", did)
    assert seqtools.translate(orf) == protein + "*", f"{did}: ORF does not encode its protein"
    cas = arm_cassette(site.get("arm_id") or f"{did}:internal{n}")
    assert qc.keeps_p003_tag_codons(cas), "not p003's tags"
    assert 0 < n < len(protein), f"{did}: site {n} is not internal"
    junction = f"{protein[n - 1]}{n}|{protein[n]}{n + 1}"
    head, rest = orf[: 3 * n], orf[3 * n :]
    # A native junction that makes a Type IIS site with any encoding of the cassette's ends (MS4A2
    # after 91: ...GAGAC|G) is domesticated once, on p003's encoding, by one synonymous change in
    # the two native codons either side (build_l1_nterm_arms.seam_fix), so the provisional and the
    # real build share one flank and the drawn encoding is screened against it. (2026-11-22)
    seam_edit = ""
    if not qc.wrapper_sites_only(qc.wrap_insert(head + CASSETTE + rest)):
        head, rest, seam_edit = nterm.seam_fix(head, rest, CASSETTE)
    coding = head + cas + rest
    expect = protein[:n] + CASSETTE_AA + protein[n:]
    translated = seqtools.translate(coding)
    assert translated == expect + "*", f"{did}: protein is not native[:n] + cassette + native[n:]"
    assert translated.count("*") == 1, f"{did}: internal stop"
    insert = qc.wrap_insert(coding)
    clean = qc.wrapper_sites_only(insert)
    provisional = os.environ.get(PROVISIONAL_ENV) == "1" and cas == CASSETTE
    if not clean and not provisional:
        return f"Type IIS site at a cassette seam: {seqtools.find_sites(insert, ENZYMES)}"
    # A new arm's provisional build (p003's shared encoding) is written even when that encoding
    # makes a seam site, so diversify_arm_cassettes.py can find its flanks and draw an encoding
    # that avoids it (its draw screens the junctions for the same enzymes); the real build
    # re-checks. Never an order. (2026-11-22, MS4A2 l1-p1253-a1.)
    lo = qc.CODING_START + len(head)
    arm = {k: site.get(k, "") for k in ARM_COLUMNS} if site.get("arm_id") else {}
    if arm:
        arm["construct_key"] = arm["construct_key"] or f"{did}:internal{n}"
    return {
        **arm,
        "provisional_seam": provisional and not clean,
        "construct_id": site.get("arm_id") or f"{did}:internal{n}",
        "seam_edit": seam_edit,
        "design_id": did,
        "candidate_id": member["candidate_id"],
        "gene_symbol": member["gene_symbol"],
        "candidate_set": member["candidate_set"],
        "species": member["species"],
        "decision": site.get("decision_member") or site["decision"],
        "pick_rank": site.get("pick_rank", ""),
        "confidence": site.get("confidence", ""),
        "refold_verdict": site.get("refold_verdict", ""),
        "status": "decided",
        "native_aa": len(protein),
        "protein_aa": len(expect),
        "orf_nt": len(orf),
        "coding_nt": len(coding),
        "insert_nt": len(insert),
        "placement": "internal",
        "cassette_after_residue": n,
        "junction": junction,
        "cassette_start_nt": lo + 1,
        "type_iis_sites": "Esp3I x2 (wrapper)",
        "snorkel_candidate": site.get("snorkel_candidate", ""),
        "site_source": (
            f"{ARM_DECISIONS.name} ({site['arm_role']})"
            if site.get("arm_id")
            else f"{DECISIONS.name} (PR #245)"
        ),
        "decided_by": site["decided_by"],
        "implemented_by": site["implemented_by"],
        "reason": site["reason"] + (f"; seam domesticated: {seam_edit}" if seam_edit else ""),
        "built_by": BUILT_BY,
        "native_sha256": member["sequence_sha256"],
        "protein_sha256": hashlib.sha256(expect.encode()).hexdigest(),
        "insert": insert,
    }


def _load_internal_members() -> tuple[list[dict], dict[str, dict]]:
    """Active internal-rung members and manifest with strict key uniqueness checks."""
    mapping_rows = read_table(MAPPING)
    _index_unique(mapping_rows, "candidate_id", MAPPING.name)
    _index_unique(mapping_rows, "design_id", MAPPING.name)
    tp_rows = read_table(TAG_PATHS)
    _index_unique(tp_rows, "candidate_id", TAG_PATHS.name)
    # candidate_id has equalled design_id since the per-insert cut (final-e); a member added
    # after the last renumbering (which runs at the cut) maps to itself.
    design = {r["candidate_id"]: r["candidate_id"] for r in tp_rows} | {
        r["candidate_id"]: r["design_id"] for r in mapping_rows
    }
    members = [
        {**r, "design_id": design[r["candidate_id"]]}
        for r in tp_rows
        if r["placement"] == "internal"
    ]
    active_manifest = [r for r in read_table(MANIFEST) if ld.is_active(r)]
    manifest = _index_unique(active_manifest, "candidate_id", MANIFEST.name)
    missing_m = sorted(m["candidate_id"] for m in members if m["candidate_id"] not in manifest)
    assert not missing_m, f"internal members missing from active manifest: {missing_m[:5]}"
    return members, manifest


def build() -> tuple[list[dict], dict]:
    from deliverome_analysis.core.s3 import processed

    members, manifest = _load_internal_members()
    proteins = ld.designed_proteins(MANIFEST, library=ld.LIBRARY_L1)
    constructs = json.loads(Path(processed.download(ld.CONSTRUCTS_KEY)).read_text())
    tag_fixed = qc.p003_tag_codons()
    calls = decisions()
    per_arm = arm_calls()
    if per_arm is not None:
        member_dids = {m["design_id"] for m in members}
        orphans = sorted(set(per_arm) - member_dids)
        assert not orphans, f"per-arm build loop decisions not on the internal rung: {orphans[:5]}"

    rows, held, no_site, not_built = [], [], [], []
    for m in sorted(members, key=lambda r: r["design_id"]):
        cid, did = m["candidate_id"], m["design_id"]
        assert cid in proteins, f"{did} ({cid}): missing designed protein"
        protein = proteins[cid]["sequence"]
        assert hashlib.sha256(protein.encode()).hexdigest() == m["sequence_sha256"], (
            f"{did}: designed protein is not the tag-path sequence"
        )
        call = calls.get(did)
        if per_arm is not None:
            raw_loops = per_arm.get(did, [])
            stale_loops = [
                a["arm_id"] for a in raw_loops if a["sequence_sha256"] != m["sequence_sha256"]
            ]
            assert not stale_loops, f"{did}: per-arm loop sequence_sha256 mismatch: {stale_loops}"
            loops = raw_loops
        else:
            loops = None
        if loops:
            # Per-arm mode: the arm rows are the decision of record. A member whose one-site
            # decision was a snorkel or a hold can now carry loop arms (e.g. all_bad: two loops
            # plus a snorkel), so the member row no longer gates the build.
            orf = constructs[manifest[cid]["key"]]["construct"]
            member_decision = call["decision"] if call else ""
            for site in loops:
                out = build_one(
                    m,
                    {**site, "decision_member": member_decision or "internal_loop"},
                    orf,
                    protein,
                    tag_fixed,
                )
                if isinstance(out, str):
                    held.append(
                        {
                            "design_id": did,
                            "arm_id": site["arm_id"],
                            "gene_symbol": m["gene_symbol"],
                            "reason": out,
                        }
                    )
                else:
                    rows.append(out)
            continue
        if call is None or call["sequence_sha256"] != m["sequence_sha256"]:
            no_site.append({"design_id": did, "gene_symbol": m["gene_symbol"]})
            continue
        if not call["decision"].startswith("internal_loop"):
            not_built.append(
                {
                    "design_id": did,
                    "gene_symbol": m["gene_symbol"],
                    "decision": call["decision"],
                    "reason": call["reason"],
                }
            )
            continue
        if per_arm is not None:
            # Per-arm mode and this member has no loop arm to build.
            not_built.append(
                {
                    "design_id": did,
                    "gene_symbol": m["gene_symbol"],
                    "decision": call["decision"],
                    "reason": "no loop arm with decision=build in " + ARM_DECISIONS.name,
                }
            )
            continue
        orf = constructs[manifest[cid]["key"]]["construct"]
        out = build_one(m, call, orf, protein, tag_fixed)
        if isinstance(out, str):
            held.append({"design_id": did, "gene_symbol": m["gene_symbol"], "reason": out})
        else:
            rows.append(out)
    cids = [r["construct_id"] for r in rows]
    assert len(cids) == len(set(cids)), "duplicate construct_id in built internal arms"
    if per_arm is not None:
        expected_arms = sum(len(v) for v in per_arm.values())
        assert len(rows) + len(held) == expected_arms, (
            f"silent arm drop: built {len(rows)} + held {len(held)} != expected {expected_arms}"
        )
    checked = [r for r in rows if not r.get("provisional_seam")]
    n_bb = nterm.vector_check(checked) if checked else 0
    members_built = {r["design_id"] for r in rows}
    summary = {
        "tranche": "l1",
        "rung": "internal",
        "members_on_rung": len(members),
        "members_with_an_arm": len(members_built),
        "arms": len(rows),
        **(
            {
                "arm_decisions": ARM_DECISIONS.name,
                "members_by_arm_count": dict(
                    sorted(Counter(Counter(r["design_id"] for r in rows).values()).items())
                ),
            }
            if per_arm is not None
            else {}
        ),
        "arms_by_decision": dict(sorted(Counter(r["decision"] for r in rows).items())),
        "arms_by_confidence": dict(sorted(Counter(r["confidence"] for r in rows).items())),
        "arms_by_refold_verdict": dict(
            sorted(Counter(r["refold_verdict"] or "none" for r in rows).items())
        ),
        "decided_no_arm_here": not_built,
        "held_at_build": held,
        "no_site": no_site,
        "type_iis": f"every insert: 2 Esp3I (wrapper), 0 BbsI, 0 BsaI; in p0039: backbone's "
        f"{n_bb} + 2",
        "insert_nt_range": (
            [min(r["insert_nt"] for r in rows), max(r["insert_nt"] for r in rows)] if rows else []
        ),
        "cassette": "p003 tandem, this arm's own linker encoding (l1_arm_cassette_encodings.tsv)",
        "constructs_key": ld.CONSTRUCTS_KEY,
        "barcode": "20-N placeholder; barcodes not reserved (issue #1)",
    }
    return rows, summary


def write(rows: list[dict], summary: dict, tsv: Path, fasta: Path, js: Path) -> None:
    columns = COLUMNS + ARM_COLUMNS if any(r.get("arm_id") for r in rows) else COLUMNS
    with tsv.open("w", newline="") as fh:
        w = csv.DictWriter(
            fh, fieldnames=columns, delimiter="\t", extrasaction="ignore", lineterminator="\n"
        )
        w.writeheader()
        w.writerows(rows)
    fasta.write_text(
        "".join(
            f">{r['construct_id']}|{r['gene_symbol']}|after{r['cassette_after_residue']}"
            f"|{r['insert_nt']}nt\n"
            + "\n".join(r["insert"][i : i + 70] for i in range(0, len(r["insert"]), 70))
            + "\n"
            for r in rows
        )
    )
    js.write_text(json.dumps(summary, indent=2) + "\n")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--check", action="store_true", help="rebuild and compare with the record")
    args = ap.parse_args(argv)
    rows, summary = build()
    if args.check:
        tmp = Path(tempfile.mkdtemp())
        write(rows, summary, tmp / "t.tsv", tmp / "f.fasta", tmp / "s.json")
        for new, rec in (
            (tmp / "t.tsv", OUT_TSV),
            (tmp / "f.fasta", OUT_FASTA),
            (tmp / "s.json", OUT_JSON),
        ):
            assert new.read_bytes() == rec.read_bytes(), f"{rec.name} does not reproduce"
        print(f"all {len(rows)} internal arms rebuild byte for byte")
        return 0
    write(rows, summary, OUT_TSV, OUT_FASTA, OUT_JSON)
    short = {
        k: v
        for k, v in summary.items()
        if k not in ("held_at_build", "no_site", "decided_no_arm_here")
    }
    print(json.dumps(short, indent=2))
    for h in summary["held_at_build"]:
        print(f"HELD {h['design_id']} {h['gene_symbol']}: {h['reason']}")
    for h in summary["decided_no_arm_here"]:
        print(f"DECIDED {h['design_id']} {h['gene_symbol']}: {h['decision']}")
    print(f"no site: {len(summary['no_site'])} members")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
