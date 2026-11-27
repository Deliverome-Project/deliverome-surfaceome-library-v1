"""Adopt the CDS-eligible ortholog panel, as reviewed for correspondence, into L1.

Becca's decision of 2026-11-10 (``ortholog_panel_adoption_judgements.tsv``): the L1 ortholog
panel becomes the Compara 116 CDS-gated redraw (``ortholog_panel_cds_eligible.tsv``, PR 262)
with the 11 forms the correspondence review recommends replacing swapped for their proposed
next-in-draw replacement (PR 273), and the 6 forms it asks about kept as drawn.

Attribution (Becca's convention of 2026-11-10): the 17 replace/keep decisions were made by Becca
Carlson on 2026-11-10 ("Keep all 6"; accept the 11 replacements) and implemented by the Claude
Code session l1-reviewer-package (a subagent of it), PR 283. Each decision row in
``ortholog_panel_adoption_judgements.tsv`` carries ``decided_by``, ``implemented_by`` and
``reason``, and every adopted panel row that a decision touched repeats them
(``decided_by`` / ``implemented_by`` / ``decision_reason``); forms taken as drawn say so.
Everything here regenerates from code and pinned inputs.

Every form's CDS is the exact pinned one: the archived final-v2 ``selected_cds.json.gz`` for a
drawn form, and the final-v2 eligibility ledger (``eligible_proteins.json.gz``, the chosen
encoding of the replacement's protein in its ortholog gene) for a replacement. Each must
translate to the panel's protein hash, and its DNA must hash to the recorded CDS hash. No CDS
is repaired here; a form without a pinned eligible CDS stops the run.

Members are reconciled the way ``redraw_panels.py`` and ``dedupe_library_members.py`` do it,
matched on the designed protein, never on accession or length:

| Form's protein is ... | Action |
|---|---|
| an active L1 ortholog member's, with byte-identical CDS | that member keeps its id (``retained_l1``) |
| an active L2 ortholog member's, with byte-identical CDS | next free L1 id; the L2 row is retired ``moved_to:<id>`` (``promoted_l2``) |
| no active member's | next free L1 id, built from its pinned CDS (``new_form``) |

An identical protein whose recorded CDS differs, or that belongs to another set or species,
stops the run rather than silently reusing the old DNA. Every other active L1 ortholog
member is retired ``dropped_from_panel:ortholog_cds_policy``. Ids are never renumbered:
fresh ids start above every live, retired and reserved L1 id (``member_id_allocations``).

The designed-CDS and construct records are append-only: the new records are the previous
records plus the new forms' entries, published under a new dated key, and an existing entry
that would change stops the run. ``library_design.CDS_KEY`` / ``CONSTRUCTS_KEY`` name them.

    uv run --locked python adopt_ortholog_cds_policy.py            # dry run: reproduce, report
    uv run --locked python adopt_ortholog_cds_policy.py --write    # write tables + manifest
    uv run --locked python adopt_ortholog_cds_policy.py --publish  # upload the new S3 keys
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import sys
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
AP = REPO / "library/l1/proteoform-selection"
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(AP))

from member_id_allocations import COMPARA116_RESERVATIONS, next_l1_member_number  # noqa: E402

from deliverome_analysis.cloning import library_design as ld  # noqa: E402
from deliverome_analysis.cloning import seqtools  # noqa: E402

PREFIX = "sprint4/s4e9-orf-design-fg-panel/"
FINAL = PREFIX + "library-design/ortholog-cds-eligibility/2026-11-05/final-v2/"
REBUILD = PREFIX + "library-design/ortholog-rebuild-inputs/2026-11-05/"
#: Where this adoption's records go. New keys only; nothing under them existed before.
#: Pinned: ``library_design.CDS_KEY`` has since moved on (2026-11-12 reference fixes).
RECORDS = PREFIX + "library-candidates/2026-11-10-ortholog-cds-policy/"
POLICY = "ortholog_cds_policy"
RETIRED_REASON = f"dropped_from_panel:{POLICY}"

#: name -> (bucket, key, sha256 of the stored object, sha256 after gunzip or None)
INPUTS = {
    "archive_manifest": (
        "processed",
        FINAL + "archive_manifest.json",
        "8323fec90f4e97339d866802154a4801e723f47bdb9a51dd439a5af3c277ccb0",
        None,
    ),
    "selected_cds": (
        "processed",
        FINAL + "results/panel/selected_cds.json.gz",
        "4f83238205b7f4fa1156d0129d40050300a794cfb9881ace35e1102a2162b222",
        "ccfc8f2d16f8e03e264428159323c1c801c1f2ba5fd7e019d2a483578a5d2a37",
    ),
    "eligible_proteins": (
        "processed",
        FINAL + "results/eligibility/eligible_proteins.json.gz",
        "2c301da225cd910aa74ae1a2aefe42815e4978cea98ca3d60da7c489461ee0d9",
        "1cabf363947d346cda1bae8b8863477bd0c9f5197cbe3efe6a1ea7f93d7fa76d",
    ),
    "provisional_picks": (
        "processed",
        FINAL + "results/scored/ortholog_provisional_picks.tsv.gz",
        "760b3a7d3acd4525df05d88c64849502747a03e5d8ca7538cdab79a484717671",
        "58def62b49dde1973a6721a20957e88371521a747da47f00d07450d48d2322f5",
    ),
    # The library as it stood before this adoption, from the dated snapshots the final-v2
    # draw itself read (ortholog_rebuild_inputs.json; final-v2 archive_manifest.json).
    "base_cds": (
        "processed",
        REBUILD + "candidate_cds.json",
        "03c3413d70e15612868133d15839e092dba20fb05f3e8318d17f48928c91c662",
        None,
    ),
    "base_constructs": (
        "processed",
        FINAL + "inputs/candidate_constructs.json.gz",
        "11d645d255037f0cb2ab855fc438df16178f04c5d90cee6fd56e22cc176e1358",
        "cdc5fa5e119e7590d910b44db03138a829f7a59e6d8ec430d8ebbe9f21008e3f",
    ),
    "pre_manifest": (
        "processed",
        REBUILD + "candidate_manifest_tranches.tsv",
        "7ad0538b831b521a10526f290f3b40eb9ec1f2333f0d674f88a831092c13cdf4",
        None,
    ),
}
PANEL = HERE / "ortholog_panel_cds_eligible.tsv"
REVIEW = HERE / "ortholog_correspondence_review.tsv"
ALTERNATES = HERE / "ortholog_correspondence_alternates.tsv"
JUDGEMENTS = HERE / "ortholog_panel_adoption_judgements.tsv"
MANIFEST = HERE / "candidate_manifest_tranches.tsv"

#: The adoption of 2026-11-10 as this script draws it. ``ortholog_panel_adopted.tsv`` is now the
#: panel of record after the 2026-11-22 narrower feature rule restored four drawn forms
#: (``write_panels_of_record.py``), so this script's own output keeps a dated name.
OUT_PANEL = HERE / "ortholog_panel_adopted_2026-11-10.tsv"
OUT_MAP = HERE / "ortholog_cds_policy_member_map.tsv"
OUT_IDS = HERE / "ortholog_cds_policy_id_allocations.tsv"
OUT_CHANGES = HERE / "panel_changes_ortholog_cds_policy.tsv"
OUT_SUMMARY = HERE / "ortholog_cds_policy_adoption_summary.json"

PANEL_COLUMNS = [
    "species",
    "band",
    "topology_class",
    "human_symbol",
    "form_id",
    "picked_acc",
    "picked_score",
    "sequence_sha256",
    "human_reference_sha256",
    "human_ensembl_gene",
    "ortholog_ensembl_gene",
    "cds_transcript",
    "selected_cds_sha256",
    "cds_source",
    "replaces",
    "decided_by",
    "implemented_by",
    "decision_reason",
]
#: Panel rows no manual decision touched: the seeded draw chose them.
DRAWN = {
    "decided_by": "no decision taken (seeded draw)",
    "implemented_by": "",
    "decision_reason": "",
}
MAP_COLUMNS = [
    "form_id",
    "species",
    "human_symbol",
    "sequence_sha256",
    "selected_cds_sha256",
    "candidate_id",
    "member_key",
    "action",
    "source_candidate_id",
    "record_transcript",
    "selected_transcript",
]


# ------------------------------------------------------------------------- inputs


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def read_tsv(path: Path) -> list[dict]:
    with Path(path).open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], columns: list[str]) -> None:
    with Path(path).open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def fetch(name: str, cache: Path) -> Path:
    """Download one pinned input to ``cache`` and verify it; gunzip if it is compressed."""
    from deliverome_analysis.core import s3

    bucket, key, stored, decompressed = INPUTS[name]
    cache.mkdir(parents=True, exist_ok=True)
    path = cache / Path(key).name
    if not path.exists():
        getattr(s3, bucket).download(key, path)
    if not stored or sha256(path.read_bytes()) != stored:
        raise SystemExit(f"{name}: {bucket}://{key} does not match its pinned SHA-256")
    if decompressed is None:
        return path
    raw = gzip.decompress(path.read_bytes())
    if sha256(raw) != decompressed:
        raise SystemExit(f"{name}: decompressed content does not match its pinned SHA-256")
    out = cache / path.name.removesuffix(".gz")
    out.write_bytes(raw)
    return out


def check_archive(archive: dict) -> None:
    """The committed review panel is the archived final-v2 panel, byte for byte."""
    rel = FINAL + "results/panel/ortholog_panel_cds_eligible.tsv"
    if archive["files"][rel]["sha256"] != sha256(PANEL.read_bytes()):
        raise SystemExit(f"{PANEL.name} differs from the archived final-v2 panel")
    for name in ("selected_cds", "eligible_proteins", "provisional_picks"):
        key = INPUTS[name][1]
        if archive["files"][key]["sha256"] != INPUTS[name][2]:
            raise SystemExit(f"{name}: pin differs from the final-v2 archive manifest")


# ------------------------------------------------------------------ adopted panel


def protein_of(cds: str, expected: str, label: str) -> str:
    aa = seqtools.translate(cds.upper())
    if not aa.endswith("*") or "*" in aa[:-1]:
        raise SystemExit(f"{label}: CDS needs one terminal stop and no internal stop")
    if sha256(aa[:-1].encode()) != expected:
        raise SystemExit(f"{label}: CDS does not translate to the selected protein")
    return aa[:-1]


def adopted_panel(
    panel: list[dict],
    review: list[dict],
    alternates: list[dict],
    judgements: list[dict],
    selected: dict,
    eligible: dict,
    picks: list[dict],
) -> tuple[list[dict], dict]:
    """The 200 adopted forms (replacements in their replaced form's place) and their CDS."""
    rev = {r["form_id"]: r for r in review}
    if [r["form_id"] for r in panel] != [r["form_id"] for r in panel if r["form_id"] in rev]:
        raise SystemExit("the correspondence review does not cover the panel")
    if len(rev) != len(panel):
        raise SystemExit("the correspondence review covers forms outside the panel")
    decided = {j["form_id"]: j for j in judgements}
    flagged = {f for f, r in rev.items() if r["recommendation"] in ("replace", "ask")}
    if set(decided) != flagged:
        raise SystemExit(
            f"every replace/ask form needs exactly one decision: missing "
            f"{sorted(flagged - set(decided))}, extra {sorted(set(decided) - flagged)}"
        )
    for f, j in decided.items():
        r = rev[f]
        if j["review_recommendation"] != r["recommendation"]:
            raise SystemExit(f"{f}: decision written against a different review result")
        if j["decision"] == "replace":
            if j["adopted_form_id"] != r["proposed_replacement"] or not j["adopted_form_id"]:
                raise SystemExit(f"{f}: replacement is not the review's proposed replacement")
        elif j["decision"] == "keep":
            if j["adopted_form_id"] != f:
                raise SystemExit(f"{f}: a kept form must adopt itself")
        else:
            raise SystemExit(f"{f}: unknown decision {j['decision']!r}")
    alt = {a["form_id"]: a for a in alternates}
    pick = {f"{p['species']}:{p['gene']}:{p['picked_acc']}": p for p in picks}

    rows, cds = [], {}
    for r in panel:
        j = decided.get(r["form_id"])
        if j and j["decision"] == "replace":
            fid = j["adopted_form_id"]
            a, p = alt.get(fid), pick.get(fid)
            if a is None or a["recommendation"] != "keep":
                raise SystemExit(f"{fid}: replacement is not a reviewed architecture-keeping form")
            if p is None or p["picked_sequence_sha256"] != a["sequence_sha256"]:
                raise SystemExit(f"{fid}: replacement is not its gene's scored representative")
            if (a["species"], a["band"], a["topology_class"]) != (
                r["species"],
                r["band"],
                r["topology_class"],
            ):
                raise SystemExit(f"{fid}: replacement is outside the replaced form's stratum")
            by_gene = eligible[a["species"]].get(a["sequence_sha256"], {}).get("by_gene", {})
            enc = by_gene.get(p["ortholog_ensembl_gene"])
            if (
                not enc
                or enc["eligibility_status"] != "eligible"
                or not enc.get("chosen_encoding")
                or enc["cds_sha256"] != p["picked_cds_sha256"]
                or enc["transcript"] != p["picked_cds_transcript"]
            ):
                raise SystemExit(f"{fid}: no pinned eligible CDS for the replacement; stopping")
            row = {
                "species": a["species"],
                "band": a["band"],
                "topology_class": a["topology_class"],
                "human_symbol": a["human_symbol"],
                "form_id": fid,
                "picked_acc": a["picked_acc"],
                "picked_score": a["picked_score"],
                "sequence_sha256": a["sequence_sha256"],
                "human_reference_sha256": a["human_reference_sha256"],
                "human_ensembl_gene": p["human_ensembl_gene"],
                "ortholog_ensembl_gene": p["ortholog_ensembl_gene"],
                "cds_transcript": enc["transcript"],
                "selected_cds_sha256": enc["cds_sha256"],
                "cds_source": "final-v2 eligible_proteins (chosen encoding)",
                "replaces": r["form_id"],
            }
        else:
            enc = selected.get(r["form_id"])
            if not enc or enc["cds_sha256"] != r["selected_cds_sha256"]:
                raise SystemExit(f"{r['form_id']}: no pinned selected CDS")
            row = {
                **{c: r[c] for c in PANEL_COLUMNS if c in r},
                "cds_transcript": enc["transcript"],
                "cds_source": "final-v2 selected_cds",
                "replaces": "",
            }
        row.update(
            {
                "decided_by": f"{j['decided_by']}, {j['date']}",
                "implemented_by": j["implemented_by"],
                "decision_reason": j["reason"],
            }
            if j
            else DRAWN
        )
        if sha256(enc["cds"].encode()) != row["selected_cds_sha256"]:
            raise SystemExit(f"{row['form_id']}: CDS bytes differ from the recorded CDS hash")
        if enc["species"] != row["species"]:
            raise SystemExit(f"{row['form_id']}: CDS species differs")
        protein_of(enc["cds"], row["sequence_sha256"], row["form_id"])
        rows.append(row)
        cds[row["form_id"]] = enc
    if len({r["sequence_sha256"] for r in rows}) != len(rows):
        raise SystemExit("two adopted forms encode the same protein")
    if Counter(r["species"] for r in rows) != Counter({"mouse": 100, "cyno": 100}):
        raise SystemExit("the adopted panel must hold 100 mouse and 100 cyno forms")
    return rows, cds


# ------------------------------------------------------------- member reconciliation


def protein_hash(entry: dict) -> str:
    return sha256(seqtools.translate(entry["cds"]).rstrip("*").encode())


def cds_entry(form: dict, enc: dict) -> dict:
    """The designed-CDS record entry for a new form, in the record's own schema."""
    return {
        "cds": enc["cds"],
        "origin": f"ortholog_{form['species']}",
        "route": f"Ensembl116 CDS-eligible ({POLICY}, final-v2)",
        "species": form["species"],
        "transcript": enc["transcript"],
        "variant_residues": "0",
    }


def reconcile(
    panel: list[dict], forms_cds: dict, manifest: list[dict], base_cds: dict
) -> tuple[list[dict], list[dict], dict]:
    """(member map, a copy of the manifest rows, new CDS entries by key)."""
    rows = [dict(r) for r in manifest]
    by_protein: dict[str, list[dict]] = {}
    for r in rows:
        if ld.is_active(r) and r["key"] in base_cds:
            by_protein.setdefault(protein_hash(base_cds[r["key"]]), []).append(r)
    keys = {r["key"] for r in rows}
    # The reservations that existed when this adoption allocated (not its own ledger, which
    # would move every id on a rerun).
    next_n = next_l1_member_number(rows, [COMPARA116_RESERVATIONS])
    mapping, new_cds = [], {}
    for form in panel:
        fid, enc = form["form_id"], forms_cds[form["form_id"]]
        hits = by_protein.get(form["sequence_sha256"], [])
        wrong = [h for h in hits if h["set"] != "ortholog" or h["species"] != form["species"]]
        if wrong:
            raise SystemExit(f"{fid}: protein duplicates {[h['candidate_id'] for h in wrong]}")
        l1 = [h for h in hits if h["library"] == ld.LIBRARY_L1]
        l2 = [h for h in hits if h["library"] == ld.LIBRARY_L2]
        if len(l1) > 1 or (not l1 and len(l2) > 1):
            raise SystemExit(f"{fid}: more than one active member already has this protein")
        src = l1[0] if l1 else (l2[0] if l2 else None)
        if src is not None and base_cds[src["key"]]["cds"] != enc["cds"]:
            # Same protein, different DNA: the pinned CDS is the one to build, and swapping
            # an existing key's CDS would silently change that member's record.
            raise SystemExit(f"{fid}: {src['candidate_id']} encodes it with a different CDS")
        if src is not None and l1:
            cid, key, action = src["candidate_id"], src["key"], "retained_l1"
        else:
            cid = f"l1-p{next_n:04d}"
            next_n += 1
            if src is not None:
                key, action = src["key"], "promoted_l2"
            else:
                key, action = "ort:" + fid, "new_form"
                if key in keys:
                    species, gene, acc = fid.split(":", 2)
                    key = f"ort:{species}:{gene}:{POLICY}:{acc}"
                    if key in keys:
                        raise SystemExit(f"{fid}: versioned key {key} already exists")
                if key in base_cds:
                    raise SystemExit(f"{key}: already in the designed-CDS record")
                new_cds[key] = cds_entry(form, enc)
            keys.add(key)
        mapping.append(
            {
                "form_id": fid,
                "species": form["species"],
                "human_symbol": form["human_symbol"],
                "sequence_sha256": form["sequence_sha256"],
                "selected_cds_sha256": form["selected_cds_sha256"],
                "candidate_id": cid,
                "member_key": key,
                "action": action,
                "source_candidate_id": src["candidate_id"] if src else "",
                "record_transcript": base_cds[src["key"]]["transcript"] if src else "",
                "selected_transcript": enc["transcript"],
            }
        )
    return mapping, rows, new_cds


def build_rows(
    mapping: list[dict], rows: list[dict], cols: list[str], constructs: dict, backbone: str
) -> tuple[list[dict], list[dict], list[str]]:
    """Apply the map to the manifest: retire, move and append. Returns (rows, changes, fasta)."""
    import assemble_candidates

    bb_sites = seqtools.find_sites(backbone, assemble_candidates.ENZYMES)
    by_id = {r["candidate_id"]: r for r in rows}
    keep = {m["candidate_id"] for m in mapping if m["action"] == "retained_l1"}
    changes, added, fasta = [], [], []
    for r in rows:
        if (
            ld.is_active(r)
            and r["library"] == ld.LIBRARY_L1
            and r["set"] == "ortholog"
            and r["candidate_id"] not in keep
        ):
            r[ld.RETIRED_COLUMN] = RETIRED_REASON
            changes.append(
                {
                    "change": "dropped",
                    "candidate_id": r["candidate_id"],
                    "key": r["key"],
                    "detail": "not in the adopted ortholog panel",
                }
            )
    for m in mapping:
        if m["action"] == "retained_l1":
            changes.append(
                {
                    "change": "kept",
                    "candidate_id": m["candidate_id"],
                    "key": m["member_key"],
                    "detail": f"{m['form_id']}: identical protein and CDS",
                }
            )
            continue
        if m["action"] == "promoted_l2":
            src = by_id[m["source_candidate_id"]]
            src[ld.RETIRED_COLUMN] = f"moved_to:{m['candidate_id']}"
            row = {
                **src,
                "candidate_id": m["candidate_id"],
                "library": ld.LIBRARY_L1,
                "panel": "ortholog_panel",
                "prior_candidate_id": src["candidate_id"],
                ld.RETIRED_COLUMN: "",
            }
            detail = f"{m['form_id']}: from {src['candidate_id']} (identical protein and CDS)"
        else:
            con = constructs[m["member_key"]]
            vector, built, clean = assemble_candidates.assemble_one(
                backbone, bb_sites, m["candidate_id"], m["member_key"], con
            )
            if not clean:
                raise SystemExit(f"{m['member_key']}: vector has an unexpected Type IIS pattern")
            row = {
                **dict.fromkeys(cols, ""),
                **{c: str(v) for c, v in built.items() if c in cols},
                "library": ld.LIBRARY_L1,
                "panel": "ortholog_panel",
                "prior_candidate_id": m["candidate_id"],
                ld.RETIRED_COLUMN: "",
            }
            fasta.append(
                assemble_candidates.fasta_record(m["candidate_id"], m["member_key"], vector)
            )
            detail = f"{m['form_id']}: new form, pinned final-v2 CDS"
        added.append(row)
        changes.append(
            {
                "change": "added",
                "candidate_id": m["candidate_id"],
                "key": m["member_key"],
                "detail": detail,
            }
        )
    result = rows + added
    ids = [r["candidate_id"] for r in result]
    if len(ids) != len(set(ids)):
        raise SystemExit("candidate ids are not unique after adoption")
    active = {
        r["candidate_id"]
        for r in result
        if ld.is_active(r) and r["library"] == ld.LIBRARY_L1 and r["set"] == "ortholog"
    }
    if active != {m["candidate_id"] for m in mapping}:
        raise SystemExit("active L1 ortholog members differ from the member map")
    return result, changes, fasta


def domesticate_new(new_cds: dict) -> dict:
    import assemble_candidates
    import domesticate_candidates

    out = {}
    for key, entry in new_cds.items():
        con = domesticate_candidates.domesticate_one(entry)
        if con is None or seqtools.translate(con["construct"]) != seqtools.translate(entry["cds"]):
            raise SystemExit(f"{key}: domestication would change the protein")
        if con["sites_after"] or seqtools.find_sites(con["construct"], assemble_candidates.ENZYMES):
            raise SystemExit(f"{key}: domestication left an internal Type IIS site")
        out[key] = con
    return out


def merge(base: dict, extra: dict, label: str) -> dict:
    clash = sorted(k for k in extra if k in base and base[k] != extra[k])
    if clash:
        raise SystemExit(f"{label}: would change existing entries {clash[:5]}")
    return {**base, **extra}


# ---------------------------------------------------------------------------- run


def reproduce(work: Path, manifest_path: Path | None = None) -> dict:
    """Everything the adoption decides, from pinned inputs only."""
    import candidate_inputs
    from Bio import SeqIO

    cache = work / "inputs"
    archive = json.loads(fetch("archive_manifest", cache).read_text())
    check_archive(archive)
    selected = json.loads(fetch("selected_cds", cache).read_text())
    eligible = json.loads(fetch("eligible_proteins", cache).read_text())
    picks = read_tsv(fetch("provisional_picks", cache))
    base_cds = json.loads(fetch("base_cds", cache).read_text())
    base_con = json.loads(fetch("base_constructs", cache).read_text())
    if manifest_path is None:
        manifest_path = fetch("pre_manifest", cache)
    with Path(manifest_path).open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        cols, manifest = list(reader.fieldnames or []), list(reader)
    if any(r["key"] not in base_cds for r in manifest if r["set"] == "ortholog"):
        raise SystemExit("an ortholog member has no designed CDS in the base record")

    panel, forms_cds = adopted_panel(
        read_tsv(PANEL),
        read_tsv(REVIEW),
        read_tsv(ALTERNATES),
        read_tsv(JUDGEMENTS),
        selected,
        eligible,
        picks,
    )
    mapping, rows, new_cds = reconcile(panel, forms_cds, manifest, base_cds)
    new_con = domesticate_new(new_cds)
    backbone = str(SeqIO.read(candidate_inputs.backbone(), "genbank").seq).upper()
    result, changes, fasta = build_rows(mapping, rows, cols, new_con, backbone)
    cds_record = merge(base_cds, new_cds, "designed CDS")
    con_record = merge(base_con, new_con, "constructs")
    counts = Counter(m["action"] for m in mapping)
    summary = {
        "decision": "Becca, 2026-11-10: adopt the CDS-eligible panel with the 11 correspondence "
        "replacements; keep the 6 ask forms as drawn",
        "panel_forms": len(panel),
        "by_species": dict(sorted(Counter(p["species"] for p in panel).items())),
        "replacements": {p["replaces"]: p["form_id"] for p in panel if p["replaces"]},
        "actions": dict(sorted(counts.items())),
        "retired_old_l1_orthologs": sum(c["change"] == "dropped" for c in changes),
        "l2_rows_moved_to_l1": counts["promoted_l2"],
        "fresh_ids": [m["candidate_id"] for m in mapping if m["action"] != "retained_l1"],
        "transcript_differs_from_record": [
            m["candidate_id"]
            for m in mapping
            if m["record_transcript"] and m["record_transcript"] != m["selected_transcript"]
        ],
        "designed_cds_entries": {"base": len(base_cds), "added": len(new_cds)},
        "construct_entries": {"base": len(base_con), "added": len(new_con)},
        "domestication_edits_in_new_forms": sum(c["edits"] for c in new_con.values()),
        "inputs": {
            name: {"source": f"{b}://{k}", "sha256": h} for name, (b, k, h, _) in INPUTS.items()
        },
        "committed_inputs": {
            p.name: sha256(p.read_bytes()) for p in (PANEL, REVIEW, ALTERNATES, JUDGEMENTS)
        },
        "records": {
            "candidate_cds": f"processed://{RECORDS}candidate_cds.json",
            "candidate_constructs": f"processed://{RECORDS}candidate_constructs.json",
            "added_vectors": f"processed://{RECORDS}added_vectors.fasta",
        },
    }
    return {
        "panel": panel,
        "mapping": mapping,
        "manifest": result,
        "cols": cols,
        "changes": changes,
        "cds_record": cds_record,
        "con_record": con_record,
        "fasta": fasta,
        "summary": summary,
    }


def tsv_text(rows: list[dict], cols: list[str]) -> str:
    lines = ["\t".join(cols)] + ["\t".join(str(r.get(c, "")) for c in cols) for r in rows]
    return "\n".join(lines) + "\n"


def record_bytes(data: dict) -> bytes:
    return json.dumps(data, indent=0, sort_keys=True).encode()


def outputs(r: dict) -> dict[Path, str]:
    """Committed outputs and their exact text."""
    change_cols = ["change", "candidate_id", "key", "detail"]
    return {
        OUT_PANEL: tsv_text(r["panel"], PANEL_COLUMNS),
        OUT_MAP: tsv_text(r["mapping"], MAP_COLUMNS),
        OUT_IDS: tsv_text([m for m in r["mapping"] if m["action"] != "retained_l1"], MAP_COLUMNS),
        OUT_CHANGES: tsv_text(r["changes"], change_cols),
        MANIFEST: tsv_text(r["manifest"], r["cols"]),
        OUT_SUMMARY: json.dumps(
            {
                **r["summary"],
                "record_sha256": {
                    "candidate_cds.json": sha256(record_bytes(r["cds_record"])),
                    "candidate_constructs.json": sha256(record_bytes(r["con_record"])),
                    "added_vectors.fasta": sha256("".join(r["fasta"]).encode()),
                },
            },
            indent=2,
        )
        + "\n",
    }


def unnumbered(text: str) -> str:
    """The manifest as this stage wrote it: design ids as ``candidate_id``, no ``design_id``.

    ``renumber_library_members.py`` runs after every membership change and re-keys active L1
    rows to member ids (``l1-mNNNN``), keeping the design id beside them; it verifies that
    numbering itself. Undoing it here compares the adoption with what adoption decided.
    """
    import csv
    import io

    from deliverome_analysis.cloning import library_design as ld

    rows = list(csv.DictReader(io.StringIO(text), delimiter="\t"))
    if not rows or ld.DESIGN_ID_COLUMN not in rows[0]:
        return text
    cols = [c for c in rows[0] if c != ld.DESIGN_ID_COLUMN]
    return tsv_text([{**r, "candidate_id": ld.design_id(r)} for r in rows], cols)


def ortholog_rows(text: str) -> str:
    """The manifest's header and ortholog rows: the rows this adoption decides.

    Later stages change other rows (2026-11-12: reference fixes, the gated isoform redraw and
    the duplicates they retire), so the adoption is checked on its own rows, which must stay
    byte-identical.
    """
    lines = text.splitlines(keepends=True)
    col = lines[0].rstrip("\n").split("\t").index("set")
    return lines[0] + "".join(x for x in lines[1:] if x.split("\t")[col] == "ortholog")


def masked_trims(text: str) -> str:
    """The manifest text with the design columns of the applied cyno restarts blanked.

    ``ortholog_start_trims.tsv`` (decided 2026-11-12 by l1-proteoform-audit and
    ortholog-correspondence-review, each delegated by Becca; an explicit, recorded exception to
    the policy's ban on start repairs) restarts some adopted cyno forms at the native Met
    aligned to human Met1. ``apply_ortholog_start_trims.py`` owns and verifies those columns;
    every other column of those rows, and every other ortholog row, must still reproduce here.
    """
    import apply_ortholog_start_trims as trims

    keys = trims.applied_keys()
    lines = text.splitlines(keepends=True)
    hdr = lines[0].rstrip("\n").split("\t")
    kcol, blank = hdr.index("key"), {hdr.index(c) for c in trims.DESIGN_COLS if c in hdr}
    out = [lines[0]]
    for line in lines[1:]:
        f = line.rstrip("\n").split("\t")
        if f[kcol] in keys:
            f = ["" if i in blank else x for i, x in enumerate(f)]
            line = "\t".join(f) + "\n"
        out.append(line)
    return "".join(out)


def verify(work: Path) -> str:
    """Rerun the adoption from pinned inputs; every committed output and record must match."""
    from deliverome_analysis.core.s3 import processed

    r = reproduce(work)
    for path, text in outputs(r).items():
        if path == MANIFEST:
            committed = masked_trims(ortholog_rows(unnumbered(path.read_text())))
            text = masked_trims(ortholog_rows(text))
        else:
            committed = path.read_text()
        if committed != text:
            raise AssertionError(f"{path.name} does not reproduce")
    for name, data in (
        ("candidate_cds.json", record_bytes(r["cds_record"])),
        ("candidate_constructs.json", record_bytes(r["con_record"])),
        ("added_vectors.fasta", "".join(r["fasta"]).encode()),
    ):
        got = Path(processed.download(RECORDS + name, work / f"rec_{name}")).read_bytes()
        if got != data:
            raise AssertionError(f"processed://{RECORDS}{name} does not reproduce")
    s = r["summary"]
    return (
        f"{s['panel_forms']} adopted forms ({len(s['replacements'])} replaced by review) "
        f"reproduce from final-v2; {s['actions']}; {s['retired_old_l1_orthologs']} old L1 "
        f"orthologs retired; CDS record = base {s['designed_cds_entries']['base']} + "
        f"{s['designed_cds_entries']['added']} appended"
    )


def main(argv=None) -> int:
    import tempfile

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--write", action="store_true", help="write the committed tables")
    ap.add_argument("--publish", action="store_true", help="upload the new S3 records")
    ap.add_argument(
        "--manifest",
        type=Path,
        help="pre-adoption manifest (default: the pinned 2026-11-05 snapshot)",
    )
    args = ap.parse_args(argv)
    work = Path(tempfile.mkdtemp(prefix="ortholog-cds-policy-"))
    r = reproduce(work, args.manifest)
    print(json.dumps({k: v for k, v in r["summary"].items() if k != "inputs"}, indent=2))
    if args.write:
        for path, text in outputs(r).items():
            path.write_text(text)
            print(f"wrote {path.name}")
    if args.publish:
        from deliverome_analysis.core.s3 import processed

        for path, text in outputs(r).items():
            if path.read_text() != text:
                raise SystemExit(f"{path.name} differs from the rerun; --write and commit first")
        for name, data in (
            ("candidate_cds.json", record_bytes(r["cds_record"])),
            ("candidate_constructs.json", record_bytes(r["con_record"])),
            ("added_vectors.fasta", "".join(r["fasta"]).encode()),
        ):
            key = RECORDS + name
            if processed.exists(key):
                raise SystemExit(f"processed://{key} already exists; records are never overwritten")
            path = work / name
            path.write_bytes(data)
            processed.upload(str(path), key)
            print(f"published processed://{key}")
    if not (args.write or args.publish):
        print("dry run -- --write writes the tables; commit, then --publish uploads the records")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
