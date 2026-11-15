"""Check every ordered L1 member's built ORF against the protein the manifest names.

The ORF that is synthesized is the domesticated ``construct`` in the construct record
(``library_design.CONSTRUCTS_KEY``), not the CDS in the CDS record. For each member of
``reviewer-package/l1_order_plan.tsv`` this checks, from the construct itself:

* it starts ATG, ends in a stop, has no internal stop and is a whole number of codons;
* its translation's SHA-256 and length equal ``library_manifest_l1.tsv``'s
  ``sequence_sha256`` / ``ordered_length_aa`` (expression-rescue arms, which that manifest
  does not carry, are checked against ``rescue_constructs.tsv`` instead);
* the construct record names a source;
* the order plan's ``insert_nt`` equals the ORF plus the fixed 90-nt wrapper
  (``assemble_candidates.assemble_one``: Esp3I sites, Kozak, reverse primer, CS1), so plate
  order and pooling volumes rest on the current ORF.

Writes ``l1_built_orf_check.tsv`` (one row per member) and exits 1 on any failure other
than a stale ``insert_nt``, which is reported for the order-plan owner.

    uv run --locked python analysis/sprint4/s4e9-orf-design-fg-panel/check_l1_built_orfs.py
"""

from __future__ import annotations

import csv
import hashlib
import json
import sys
from collections import Counter
from pathlib import Path

from Bio.Seq import Seq

from deliverome_analysis import library_design as ld
from deliverome_analysis.s3 import processed

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
PLAN = HERE / "reviewer-package" / "l1_order_plan.tsv"
MANIFEST = HERE / "library_manifest_l1.tsv"
RESCUE = REPO / "scratch/becca/appris-principal/rescue_constructs.tsv"
OUT = HERE / "l1_built_orf_check.tsv"
WRAPPER_NT = 90  # Esp3I + Kozak + ... + REV_PRIMER + CS1 + Esp3I around the ORF


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main() -> int:
    constructs = json.loads(Path(processed.download(ld.CONSTRUCTS_KEY)).read_text())
    manifest = {r["manifest_key"]: r for r in read_tsv(MANIFEST)}
    rescue = {r["symbol"]: r["construct"] for r in read_tsv(RESCUE)}
    rows = []
    for p in read_tsv(PLAN):
        key = p["key"]
        rec = constructs.get(key, {})
        orf = rec.get("construct", "")
        prot = str(Seq(orf).translate()) if orf else ""
        body = prot[:-1] if prot.endswith("*") else prot
        sha = hashlib.sha256(body.encode()).hexdigest()
        m = manifest.get(key)
        if m is not None:
            ref, protein_ok = (
                "library_manifest_l1.tsv",
                (m["sequence_sha256"] == sha and int(m["ordered_length_aa"]) == len(body)),
            )
        elif p["set"] == "expression_rescue":
            ref, protein_ok = "rescue_constructs.tsv", orf == rescue.get(key.split(":")[1])
        else:
            ref, protein_ok = "none", False
        checks = {
            "construct_present": bool(orf),
            "atg": orf[:3] == "ATG",
            "stop": prot.endswith("*"),
            "no_internal_stop": "*" not in body,
            "frame": len(orf) % 3 == 0 and bool(orf),
            "protein_matches": protein_ok,
            "source_recorded": bool(rec.get("source")),
        }
        insert_ok = int(p["insert_nt"]) == len(orf) + WRAPPER_NT
        failed = [k for k, ok in checks.items() if not ok]
        rows.append(
            {
                "design_id": p["design_id"],
                "key": key,
                "set": p["set"],
                "orf_nt": len(orf),
                "protein_aa": len(body),
                "protein_sha256": sha,
                "checked_against": ref,
                "source": rec.get("source", ""),
                "order_plan_insert_nt": p["insert_nt"],
                "insert_nt_current": insert_ok,
                "failed": ";".join(failed),
            }
        )
    with OUT.open("w", newline="") as fh:
        w = csv.DictWriter(fh, list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(f"{len(rows)} members; construct record {ld.CONSTRUCTS_KEY}")
    print("failures:", Counter(f for r in rows for f in r["failed"].split(";") if f) or "none")
    stale = [r["design_id"] for r in rows if not r["insert_nt_current"]]
    print(f"order-plan insert_nt not matching the current ORF: {len(stale)} {stale}")
    return 1 if any(r["failed"] for r in rows) else 0


if __name__ == "__main__":
    sys.exit(main())
