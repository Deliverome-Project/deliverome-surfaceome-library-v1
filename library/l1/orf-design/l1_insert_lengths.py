"""Lengths behind section 1's "median CDS length": the CDS alone, with its tag, and the insert.

Three nested lengths, each from its table of record:

* **CDS** -- ``orf_nt`` of the 2,710 active L1 members in ``candidate_manifest_tranches.tsv``
  (the number the brief quotes, ``brief_numbers.py`` ``l1_active_orf_nt_median``). One per CDS.
* **Tagged ORF** -- the coding sequence as ordered, ATG to stop: the CDS with its 51-aa tag
  cassette (153 nt) inserted, and on Igκ rungs the leader in place of the native residues it
  replaces. Translated from the recoded order of record (pinned in
  ``audit_domestication_barcodes``). One per ordered insert.
* **Insert** -- the whole ordered insert: Esp3I wrappers, Kozak, tagged ORF, stop, CC spacer,
  19-nt handle, 20-nt barcode and the 3' readout block. Equal to ``ordered_insert_nt`` in
  ``l1_insert_numbering.tsv``, which the pooling plan uses. One per ordered insert.

Per-insert lengths are also reported per CDS (the median over a CDS's constructs) so all three
can be read on the 2,710.

Written by the Claude Code session l1-swap-risk (2026-11-20) for Becca's comment on the brief's
"of 1,401 nt": "Add the median length with the bc and tag as well for context".

    uv run --locked python l1_insert_lengths.py
"""

from __future__ import annotations

import csv
import json
from collections import defaultdict
from pathlib import Path

import audit_domestication_barcodes as audit
import numpy as np

from deliverome_analysis.cloning import library_design as ld

HERE = Path(__file__).resolve().parent
OUT = HERE / "l1_insert_lengths_summary.json"


def stats(x) -> dict:
    x = np.asarray(x)
    return {
        "n": len(x),
        "median": float(np.median(x)),
        "min": int(x.min()),
        "max": int(x.max()),
    }


def main() -> None:
    manifest = [
        r
        for r in csv.DictReader((HERE / "candidate_manifest_tranches.tsv").open(), delimiter="\t")
        if r["library"] == ld.LIBRARY_L1 and ld.is_active(r)
    ]
    numbered = {
        r["arm_id"]: r
        for r in csv.DictReader((HERE / "l1_insert_numbering.tsv").open(), delimiter="\t")
    }
    rec, digest = audit.recoded_inserts(None)
    assert digest == audit.RECODED_SHA256, "order of record changed"

    tagged, insert = [], []
    per_cds_tagged, per_cds_insert = defaultdict(list), defaultdict(list)
    for r in rec:
        coding, err = audit.coding_of(r["insert"])
        assert not err, r["key"]
        n = numbered[r["key"]]
        assert int(n["ordered_insert_nt"]) == len(r["insert"]), r["key"]
        tagged.append(len(coding))
        insert.append(len(r["insert"]))
        per_cds_tagged[n["design_id"]].append(len(coding))
        per_cds_insert[n["design_id"]].append(len(r["insert"]))

    out = {
        "cds_orf_nt_per_cds": stats([int(r["orf_nt"]) for r in manifest]),
        "tagged_orf_nt_per_insert": stats(tagged),
        "insert_nt_with_barcode_per_insert": stats(insert),
        "tagged_orf_nt_per_cds_median_of_constructs": stats(
            [np.median(v) for v in per_cds_tagged.values()]
        ),
        "insert_nt_per_cds_median_of_constructs": stats(
            [np.median(v) for v in per_cds_insert.values()]
        ),
        "constant_insert_overhead_nt": sorted({i - t for i, t in zip(insert, tagged, strict=True)}),
    }
    OUT.write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
