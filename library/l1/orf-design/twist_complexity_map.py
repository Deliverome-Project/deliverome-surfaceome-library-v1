"""Map Twist's complexity analysis of the first L1 draft onto the current library.

Twist scored the first L1 draft export, which held untagged ORFs only: the designed CDS
plus 90 nt of cloning flanks, with no cassette and no barcode. Its CSV is vendor output,
stored append-only at :data:`TWIST_KEY` (uploaded 2026-11-12 from the copy Becca shared).
Rows are named ``L1-O####_GENE`` and carry the design id (``l1-pNNNN``) the member had
then. Design ids are permanent, so each row joins to the current manifest on it.

For every row this writes ``twist_complexity_map.tsv``. ``status`` is one of:

* ``in_l1_same_cds``: still an active L1 member, with the designed CDS Twist scored;
* ``in_l1_cds_changed``: still active, but its designed CDS has changed since the draft
  (a reference fix, start trim or restart), so Twist's score is for a different insert;
* ``retired``: no longer an active L1 member (``retired_reason`` says why).

"The CDS Twist scored" is the entry in the pre-adoption record ``BASE_CDS_KEY``, which the
draft was exported from. Its length plus 90 must equal Twist's ``insert_nt``, and this is
asserted for every row. Expression-rescue arms are not in that record; their manifest ORF
length is checked instead, and they count as unchanged.

Requested by the coordinator session under Becca's 2026-11-12 delegation; implemented by
the Claude Code session l1-reviewer-package, PR #283.

    uv run --locked python twist_complexity_map.py
"""

from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path

from deliverome_analysis.cloning import library_design as ld
from deliverome_analysis.core.s3 import external, processed

HERE = Path(__file__).resolve().parent
MANIFEST = HERE / "candidate_manifest_tranches.tsv"
OUT = HERE / "twist_complexity_map.tsv"
SUMMARY = HERE / "twist_complexity_map_summary.json"
TWIST_KEY = "twist/l1-complexity/2026-11-12/l1_insert_complexity_analysis.csv"
TWIST_SHA256 = "b95071e0b255cddd8bd643fe007bc4bcb55e14e5576fc1495e85efa744b4d46b"
FLANK_NT = 90  # the draft export's cloning flanks, measured: insert_nt - len(CDS)


def load_twist() -> list[dict]:
    import hashlib

    path = Path(external.download(TWIST_KEY))
    assert hashlib.sha256(path.read_bytes()).hexdigest() == TWIST_SHA256, "Twist CSV changed"
    return list(csv.DictReader(path.open()))


def main() -> int:
    twist = load_twist()
    man = {r["design_id"]: r for r in csv.DictReader(MANIFEST.open(), delimiter="\t")}
    base = json.loads(Path(processed.download(ld.BASE_CDS_KEY)).read_text())
    cur = json.loads(Path(processed.download(ld.CDS_KEY)).read_text())
    rows = []
    for t in twist:
        m = man[t["design_id"]]
        # Expression-rescue arms are built outside the CDS records (rescue_constructs.tsv);
        # their ORF length in the manifest stands in for the draft's.
        draft = base.get(m["key"], {}).get("cds")
        draft_nt = len(draft) if draft else int(m["orf_nt"])
        assert draft_nt + FLANK_NT == int(t["insert_nt"]), t["name"]
        active = m["library"] == ld.LIBRARY_L1 and not m["retired_reason"]
        if not active:
            status = "retired"
        elif draft is None or cur[m["key"]]["cds"] == draft:
            status = "in_l1_same_cds"
        else:
            status = "in_l1_cds_changed"
        rows.append(
            {
                "twist_name": t["name"],
                "design_id": t["design_id"],
                "gene_symbol": t["name"].split("_", 1)[1],
                "key": m["key"],
                "set": t["set"],
                "twist_score": t["Score"],
                "primary_driver": t["primary_driver"],
                "status": status,
                "member_id": m["candidate_id"] if active else "",
                "retired_reason": m["retired_reason"],
                "draft_insert_nt": t["insert_nt"],
                "current_orf_nt": m["orf_nt"],
                "issues": t["issues"],
                "issue_regions": t["issue_regions"],
            }
        )
    with OUT.open("w", newline="") as fh:
        w = csv.DictWriter(fh, list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    summary = {
        "twist_key": f"external://{TWIST_KEY}",
        "twist_sha256": TWIST_SHA256,
        "rows": len(rows),
        "by_score_and_status": {
            f"{s}|{st}": n
            for (s, st), n in sorted(Counter((r["twist_score"], r["status"]) for r in rows).items())
        },
    }
    SUMMARY.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    for r in rows:
        if r["twist_score"] in ("Not Accepted", "High"):
            print(
                f"  {r['twist_score']:<12} {r['design_id']} {r['gene_symbol']:<10} {r['status']:<18} "
                f"{r['member_id'] or r['retired_reason'][:40]}  {r['primary_driver']}"
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
