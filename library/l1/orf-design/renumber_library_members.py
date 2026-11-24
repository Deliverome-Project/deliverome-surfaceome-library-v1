"""Number L1's ordered inserts in plate order, and key every current record on design ids.

**Per-insert numbering (from the first per-insert cut, 2026-11-15).** An internal member may
now be ordered as 2 or 3 inserts (Becca, 2026-11-14/15), so a member can no longer be a plate
position. Each built insert is numbered instead: inserts are sorted by (ordered insert length,
arm id) and numbered ``l1-i0001`` ...; insert N sits on plate ceil(N/96), wells column-major.
The numbering is written to ``l1_insert_numbering.tsv`` (design id, arm id, insert id, plate,
well, ordered length, rung, key). The member keeps one permanent id, its design id
(``l1-pNNNN``), which is now also its ``candidate_id``; the ``l1-mNNNN`` member id is retired
and survives only in dated records. The re-keying below therefore maps any ``l1-m`` id still
in a current record to its design id, once, and is the identity after that. Schema agreed by
the coordinator, l1-integration, ranker-pu-segment and l1-reviewer-package; implemented by
l1-reviewer-package.

The text below describes the member numbering that held until then (final-d, 2026-11-14).

The L1 library is ordered as purified plasmid in 96-well plates sorted by synthesized-insert
length, so a plate holds members of similar length (reviewer brief, section 8). Becca chose
(2026-11-07) to make the member id follow that order instead of carrying a separate order
number. This stage does that, from the manifest, every time it runs:

* Built active L1 members (``library == l1``, no ``retired_reason``, an insert in the order
  of record) are sorted by (ordered insert length, design id), exactly as
  ``reviewer-package/order_plan.py`` sorts them, and numbered ``l1-m0001`` ... . The ordered
  length is the tagged, synthesis-recoded insert (``synthesis_compliance.tsv``), not the
  untagged draft ``insert_nt``, which is 59-359 nt shorter depending on the rung (decided
  2026-11-14). Member N sits on plate ceil(N/96), wells column-major (A01, B01, ..., H01, A02).
* Held members (decided, but no insert yet) are numbered after every built member, by
  (``insert_nt``, design id), with no plate or well; once built they take their place by
  ordered length.
* The id a row was designed under (``l1-pNNNN``) is kept as ``design_id``, permanent and
  never reused. Every design-stage record (panel changes, retirements, adoption ledgers,
  fold and topology runs) stays keyed on it. The two spellings differ so that an old id can
  never be mistaken for a member id: ``l1-m0026`` and ``l1-p0026`` are different members.
* Retired rows and L2 rows get no member id; their ``candidate_id`` is their design id.
* ``l1_member_renumbering.tsv`` records the mapping, with plate, well, untagged and ordered
  insert length and
  the SHA-256 of each member's designed CDS.
* The L1 records that describe the current library (:data:`RECORDS`) are re-keyed to the
  member id. A record is re-keyed by replacing whole cells that hold an id, so no other
  cell changes; row counts are asserted equal, and every record carrying a
  ``sequence_sha256`` is joined back to the designed protein of the member it now names.

**Membership is not frozen yet** (an ortholog panel redraw awaits adoption; five genes lack
an intact Ensembl copy). Any membership change -- a member added, retired or moved between
tranches -- shifts the numbering, so the stage that makes the change is followed by a rerun
of this one. The numbering that ships is the one taken at the order freeze; set
:data:`ORDER_FREEZE` then, after which ``--write`` refuses to change any assigned member id.

    uv run --locked python renumber_library_members.py --check   # committed state is current
    uv run --locked python renumber_library_members.py --write   # renumber and re-key
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import re
import sys
from collections import Counter
from pathlib import Path

from deliverome_analysis import library_design as ld
from deliverome_analysis import seqtools, topology

HERE = Path(__file__).resolve().parent
AP = HERE.parents[2] / "library/l1/proteoform-selection"
MANIFEST = HERE / "candidate_manifest_tranches.tsv"
MAPPING = HERE / "l1_member_renumbering.tsv"
SUMMARY = HERE / "l1_member_renumbering_summary.json"
RESCUE_CONSTRUCTS = AP / "rescue_constructs.tsv"

#: Date of the order freeze (ISO), or None while membership can still change. Once set,
#: ``--write`` fails if any member id would be reassigned.
ORDER_FREEZE: str | None = None

PLATE_WELLS = 96
WELLS = [f"{r}{c:02d}" for c in range(1, 13) for r in "ABCDEFGH"]  # column-major

#: Records describing the current L1 library, re-keyed to member ids, with the columns that
#: hold a member id (``None``: every string in a JSON record). Everything else that names an
#: L1 id is a design-stage record and stays on design ids (see the README). Columns such as
#: ``topology_acc`` are NOT listed: a supplement topology row is keyed on the design id in
#: the topology snapshot, so that cell is a pointer into another record, not a member id.
RECORDS: dict[str, tuple[str, ...] | None] = {
    "l1_tag_paths.tsv": ("candidate_id",),
    "l1_cassette_encodings.tsv": ("candidate_id",),
    "l1_swap_risk.tsv": (
        "candidate_id",
        "partner_shared",
        "partner_whole_cassette_widened",
        "partner_linkers_only",
    ),
    "l1_swap_risk_summary.json": None,
    "l1_structure_gap.tsv": ("candidate_id",),
    "proteoform_match.tsv": ("candidate_id",),
    "residue_offsets.tsv": ("candidate_id",),
    "library_topology_index.tsv": ("candidate_id",),
    "library_structure_payload.tsv": ("candidate_id",),
    "library_manifest_l1.tsv": ("candidate_id",),
    "l1_internal_members.tsv": ("candidate_id",),
    "l1_deeptmhmm_v1_topology_gap.tsv": ("candidate_id",),
    "structure_coverage_baseline.json": None,
}

#: Records whose writer emits rows sorted by ``candidate_id``; re-sorted after re-keying so
#: a rerun of the writer reproduces them byte for byte.
SORTED_BY_ID = {"l1_cassette_encodings.tsv"}

MAPPING_COLUMNS = (
    "design_id",
    "candidate_id",
    "status",
    "inserts",
    "insert_nt",
    "key",
    "cds_sha256",
)
#: One row per ordered insert, in plate order.
INSERTS = HERE / "l1_insert_numbering.tsv"
INSERT_COLUMNS = (
    "insert_id",
    "plate",
    "well",
    "arm_id",
    "design_id",
    "rung",
    "ordered_insert_nt",
    "key",
)
#: Ordered (tagged, synthesis-recoded) insert length per design id: the length Twist makes and
#: the plates are sorted on. Recoding is synonymous, so the length equals the built insert's.
SYNTHESIS_TABLE = HERE / "synthesis_compliance.tsv"


def ordered_inserts() -> list[dict]:
    """Every insert in the order of record: arm id, design id, rung, ordered length."""
    with SYNTHESIS_TABLE.open() as fh:
        rows = [
            {
                "arm_id": r["arm_id"],
                "design_id": r["design_id"],
                "rung": r["rung"],
                "ordered_insert_nt": int(r["insert_nt"]),
            }
            for r in csv.DictReader(fh, delimiter="\t")
        ]
    arms = [r["arm_id"] for r in rows]
    assert len(set(arms)) == len(arms), "an arm id names two inserts"
    for r in rows:
        assert ld.design_of_arm(r["arm_id"]) == r["design_id"], r["arm_id"]
    return rows


def read_tsv(text: str) -> tuple[list[str], list[dict]]:
    reader = csv.DictReader(io.StringIO(text), delimiter="\t")
    return list(reader.fieldnames or []), list(reader)


def write_tsv(cols: list[str], rows: list[dict]) -> str:
    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=cols, delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(rows)
    return buf.getvalue()


def designed_cds() -> dict[str, str]:
    """Designed CDS per manifest ``key``: the candidate CDS record, plus the 25
    expression-rescue arms, whose designed coding sequence is the leader + receptor
    construct (they are not in ``candidate_cds.json``)."""
    from deliverome_analysis.s3 import processed

    with open(processed.download(ld.CDS_KEY)) as fh:
        out = {k: v["cds"] for k, v in json.load(fh).items()}
    with RESCUE_CONSTRUCTS.open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            key = f"rescue:{r['symbol']}"
            assert key not in out, key
            out[key] = r["construct"]
    return out


def number(rows: list[dict], cds: dict[str, str]) -> tuple[list[dict], list[dict]]:
    """(members, inserts): one row per active L1 member, keyed by design id, and one row per
    ordered insert in plate order."""
    # Plate order is the ordered insert length (decided 2026-11-14 by the coordinator under
    # Becca's 2026-11-12 delegation: the length that matters at the bench is the ordered one),
    # now per insert, ties broken by arm id (per-insert cut, 2026-11-15).
    inserts = sorted(ordered_inserts(), key=lambda r: (r["ordered_insert_nt"], r["arm_id"]))
    active = sorted(
        (r for r in rows if r["library"] == ld.LIBRARY_L1 and ld.is_active(r)), key=ld.design_id
    )
    key_of = {ld.design_id(r): r["key"] for r in active}
    unknown = {i["design_id"] for i in inserts} - set(key_of)
    assert not unknown, (
        f"ordered inserts for designs that are not active L1 members: {sorted(unknown)[:5]}"
    )
    per_member = Counter(i["design_id"] for i in inserts)
    members = []
    for r in active:
        did = ld.design_id(r)
        seq = cds.get(r["key"])
        if seq is None:
            raise SystemExit(f"{did} ({r['key']}): no designed CDS")
        if len(seq) != int(float(r["orf_nt"])):
            raise SystemExit(f"{did}: designed CDS is not orf_nt long")
        members.append(
            {
                "design_id": did,
                "candidate_id": did,
                "status": "built" if per_member[did] else "held",
                "inserts": per_member[did],
                "insert_nt": int(float(r["insert_nt"])),
                "key": r["key"],
                "cds_sha256": hashlib.sha256(seq.encode()).hexdigest(),
            }
        )
    numbered = [
        {
            "insert_id": ld.insert_id(i + 1),
            "plate": i // PLATE_WELLS + 1,
            "well": WELLS[i % PLATE_WELLS],
            "arm_id": ins["arm_id"],
            "design_id": ins["design_id"],
            "rung": ins["rung"],
            "ordered_insert_nt": ins["ordered_insert_nt"],
            "key": key_of[ins["design_id"]],
        }
        for i, ins in enumerate(inserts)
    ]
    return members, numbered


def renumber_manifest(cols: list[str], rows: list[dict], mapping: list[dict]) -> tuple:
    """Manifest with ``design_id`` beside ``candidate_id`` and member ids filled in."""
    to_member = {m["design_id"]: m["candidate_id"] for m in mapping}
    if ld.DESIGN_ID_COLUMN not in cols:
        cols = [cols[0], ld.DESIGN_ID_COLUMN, *cols[1:]]
    out = []
    for r in rows:
        did = ld.design_id(r)
        out.append({**r, ld.DESIGN_ID_COLUMN: did, "candidate_id": to_member.get(did, did)})
    ids = [r["candidate_id"] for r in out]
    assert len(set(ids)) == len(ids) == len(rows), "member ids are not unique"
    assert len({r[ld.DESIGN_ID_COLUMN] for r in out}) == len(rows), "design ids are not unique"
    return cols, out


class Rekey:
    """Map any L1 id (design id, or a member id from an earlier numbering) to its current id."""

    def __init__(self, old_rows: list[dict], mapping: list[dict], prior: list[dict]):
        # A member id from the numbering being replaced -> its design id. Taken from the
        # committed mapping table as well as the manifest, because a design stage may have
        # rewritten the manifest unnumbered while the records still carry member ids.
        self.to_design = {r["candidate_id"]: r["design_id"] for r in prior}
        for r in old_rows:
            if ld.MEMBER_ID_RE.fullmatch(r["candidate_id"]):
                did = ld.design_id(r)
                if self.to_design.setdefault(r["candidate_id"], did) != did:
                    raise SystemExit(
                        f"{r['candidate_id']}: the manifest and {MAPPING.name} disagree on "
                        "its design id"
                    )
        self.to_member = {m["design_id"]: m["candidate_id"] for m in mapping}
        self.changed = 0

    def cell(self, v: str) -> str:
        if ";" in v:
            return ";".join(self.cell(p) for p in v.split(";"))
        if not (ld.MEMBER_ID_RE.fullmatch(v) or ld.DESIGN_ID_RE.fullmatch(v)):
            return v
        if ld.MEMBER_ID_RE.fullmatch(v) and v not in self.to_design:
            raise SystemExit(f"{v}: a member id the manifest does not know")
        new = self.to_member.get(self.to_design.get(v, v), self.to_design.get(v, v))
        self.changed += new != v
        return new

    def value(self, v):
        if isinstance(v, str):
            if ld.MEMBER_ID_RE.fullmatch(v) or ld.DESIGN_ID_RE.fullmatch(v) or ";" in v:
                return self.cell(v)
            # Free text in a JSON record (e.g. a failure message) names members inline.
            return ld.MEMBER_ID_RE.sub(lambda m: self.cell(m.group(0)), v)
        if isinstance(v, list):
            return [self.value(x) for x in v]
        if isinstance(v, dict):
            return {self.cell(k): self.value(x) for k, x in v.items()}
        return v


def rekey_record(
    name: str, text: str, rk: Rekey, proteins: dict[str, str], id_cols: tuple[str, ...] | None
) -> str:
    before = rk.changed
    if name.endswith(".json"):
        out = json.dumps(rk.value(json.loads(text)), indent=2) + "\n"
        # Preserve the writer's own formatting when nothing changed.
        return text if rk.changed == before else out
    cols, rows = read_tsv(text)
    assert id_cols and set(id_cols) <= set(cols), f"{name}: id columns {id_cols} not in {cols}"
    new = [{c: rk.cell(r[c] or "") if c in id_cols else r[c] for c in cols} for r in rows]
    assert len(new) == len(rows), f"{name}: row count changed"
    if name in SORTED_BY_ID:
        new.sort(key=lambda r: r["candidate_id"])
    # Join each re-keyed row back to its member's designed protein.
    if "sequence_sha256" in cols and "candidate_id" in cols:
        bad = [
            r["candidate_id"]
            for r in new
            if r["candidate_id"] in proteins
            and r["sequence_sha256"]
            and r["sequence_sha256"] != proteins[r["candidate_id"]]
        ]
        assert not bad, f"{name}: sequence hash does not match the member's protein: {bad[:5]}"
    return text if rk.changed == before else write_tsv(cols, new)


def build() -> tuple[dict[Path, str], dict]:
    """Every file this stage writes, as text, plus the run summary."""
    old_text = MANIFEST.read_text()
    cols, rows = read_tsv(old_text)
    cds = designed_cds()
    mapping, inserts = number(rows, cds)
    new_cols, new_rows = renumber_manifest(cols, rows, mapping)

    if ORDER_FREEZE and INSERTS.exists():
        _, frozen = read_tsv(INSERTS.read_text())
        was = {r["arm_id"]: r["insert_id"] for r in frozen}
        now = {m["arm_id"]: m["insert_id"] for m in inserts}
        if was != now:
            moved = sorted(a for a in was.keys() | now.keys() if was.get(a) != now.get(a))
            raise SystemExit(
                f"order frozen on {ORDER_FREEZE}: renumbering would move {len(moved)} "
                f"insert(s), e.g. {moved[:5]}"
            )

    # Designed proteins by member id, for the sequence-hash join.
    by_key = {m["key"]: m["candidate_id"] for m in mapping}
    proteins = {
        by_key[k]: topology.sequence_sha256(seqtools.translate(s))
        for k, s in cds.items()
        if k in by_key
    }
    # Each member's CDS hash, joined through the new id, equals its design id's.
    by_design = {ld.design_id(r): r["key"] for r in rows}
    for m in mapping:
        assert by_design[m["design_id"]] == m["key"]

    prior = read_tsv(MAPPING.read_text())[1] if MAPPING.exists() else []
    rk = Rekey(rows, mapping, prior)
    files = {
        MANIFEST: write_tsv(new_cols, new_rows),
        MAPPING: write_tsv(list(MAPPING_COLUMNS), mapping),
        INSERTS: write_tsv(list(INSERT_COLUMNS), inserts),
    }
    rekeyed = {}
    for name, id_cols in RECORDS.items():
        path = HERE / name
        files[path] = rekey_record(name, path.read_text(), rk, proteins, id_cols)
        rekeyed[name] = len(re.findall(r"\bl1-m[0-9]{4}\b", files[path]))

    plates = Counter(i["plate"] for i in inserts)
    arms_per_member = Counter(m["inserts"] for m in mapping)
    summary = {
        "rule": "ordered inserts sorted by (ordered_insert_nt, arm_id), numbered l1-i0001.., "
        "96 per plate, column-major; members keep their design id as candidate_id",
        "order_freeze": ORDER_FREEZE,
        "members": len(mapping),
        "members_built": sum(m["status"] == "built" for m in mapping),
        "members_held": sum(m["status"] == "held" for m in mapping),
        "members_by_inserts": {str(k): v for k, v in sorted(arms_per_member.items())},
        "inserts": len(inserts),
        "inserts_by_rung": dict(Counter(i["rung"] for i in inserts)),
        "plates": len(plates),
        "last_plate_wells": plates[max(plates)] if plates else 0,
        "ordered_insert_nt": {
            "first": inserts[0]["ordered_insert_nt"] if inserts else None,
            "last": inserts[-1]["ordered_insert_nt"] if inserts else None,
            "per_plate": {
                str(p): [
                    min(i["ordered_insert_nt"] for i in inserts if i["plate"] == p),
                    max(i["ordered_insert_nt"] for i in inserts if i["plate"] == p),
                ]
                for p in sorted(plates)
            },
        },
        "members_by_set": dict(
            Counter(r["set"] for r in rows if r["library"] == "l1" and ld.is_active(r))
        ),
        "retired_l1_rows_keeping_design_id": sum(
            r["library"] == "l1" and not ld.is_active(r) for r in rows
        ),
        "member_ids_per_record": rekeyed,
    }
    assert summary["inserts"] == sum(plates.values())
    files[SUMMARY] = json.dumps(summary, indent=1) + "\n"
    return files, summary


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--write", action="store_true", help="renumber and re-key in place")
    g.add_argument(
        "--check", action="store_true", help="exit 1 unless the committed state is current"
    )
    args = ap.parse_args(argv)
    files, summary = build()
    stale = [p.name for p, text in files.items() if not p.exists() or p.read_text() != text]
    if args.check:
        if stale:
            print(f"stale: {stale}; rerun renumber_library_members.py --write", file=sys.stderr)
            return 1
        print(
            f"{summary['inserts']} inserts numbered l1-i0001.. for {summary['members_built']} "
            f"built members; every record current"
        )
        return 0
    for p, text in files.items():
        p.write_text(text)
    print(
        f"{summary['inserts']} inserts of {summary['members_built']} members on "
        f"{summary['plates']} plates; "
        f"wrote {len(stale)} changed file(s): {stale}"
    )
    print(json.dumps(summary["member_ids_per_record"], indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
