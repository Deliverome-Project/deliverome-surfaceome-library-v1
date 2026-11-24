"""Attribution table for the 13 designed CDSs whose stop codon was appended by hand-made rule.

The fix lives in ``library/l1/proteoform-selection/resolve_candidate_cds.py`` (``ADDED_STOP``):
an Ensembl 115 CDS record that ends on a sense codon gets TAA appended, and its entry in
``candidate_cds.json`` carries ``stop_added``. This script lists those entries, one row per
member, with who decided the fix, which session implemented it and why, as the 2026-11-10
attribution convention asks. It writes no design record; the CDS JSON keeps its bytes because
both reviewer bundles pin it by hash.

    uv run --locked python stop_codon_fix.py
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

from deliverome_analysis import library_design as ld
from deliverome_analysis.s3 import processed

HERE = Path(__file__).resolve().parent
OUT = HERE / "stop_codon_fix.tsv"
DECIDED_BY = "Becca, 2026-11-04: let's add stop codons then"
IMPLEMENTED_BY = "Claude Code session nterm-reviewer-package, PR #234 (bb393b72)"
REASON = (
    "the Ensembl 115 CDS record ends on a sense codon, so the ORF would read through into the "
    "vector; TAA appended because the native stop cannot be read from the pinned inputs"
)
#: The S3 versions of the records as they were before the fix (versioning is on).
PREVIOUS_VERSIONS = {
    "candidate_cds.json": "ftYKdsVZYDD8M49CEt4rDJMn.x.QmKkn",
    "candidate_constructs.json": "Nkbr70dq6lfjWkb.ASYGWZpU8yRocj_Y",
}


def main() -> int:
    with open(processed.download(ld.CDS_KEY)) as fh:
        cds = json.load(fh)
    fixed = {k for k, v in cds.items() if v.get("stop_added")}
    with (HERE / "candidate_manifest_tranches.tsv").open() as fh:
        rows = [r for r in csv.DictReader(fh, delimiter="\t") if r["key"] in fixed]
    out = [
        {
            "candidate_id": r["candidate_id"],
            "design_id": ld.design_id(r),
            "library": r["library"],
            "key": r["key"],
            "transcript": cds[r["key"]]["transcript"],
            "stop_added": cds[r["key"]]["stop_added"],
            "orf_nt": r["orf_nt"],
            "retired_reason": r["retired_reason"],
            "decided_by": DECIDED_BY,
            "implemented_by": IMPLEMENTED_BY,
            "reason": REASON,
            "previous_s3_versions": "; ".join(f"{k} {v}" for k, v in PREVIOUS_VERSIONS.items()),
        }
        for r in sorted(rows, key=lambda r: r["candidate_id"])
    ]
    keys_in_manifest = {r["key"] for r in rows}
    unlisted = sorted(fixed - keys_in_manifest)
    with OUT.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out)
    print(f"wrote {OUT.name}: {len(out)} manifest rows for {len(fixed)} fixed CDS entries")
    if unlisted:
        print(f"  fixed entries with no manifest row: {unlisted}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
