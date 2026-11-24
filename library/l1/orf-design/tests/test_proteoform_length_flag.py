"""The proteoform audit records a length difference from the canonical in both directions.

Its truncation and retention flags only count canonical residues a member lacks, so the 11
L1 references more than 10% longer or shorter than the UniProt canonical (KCNN2 791 vs 579
aa among them) were never flagged. ``length_vs_canonical`` is the column that closes that;
this checks the helper and that it catches all 11 in the committed audit table.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))

from audit_l1_proteoforms import LENGTH_DIFF_FRACTION, length_vs_canonical  # noqa: E402

#: The references more than 10% longer or shorter than the canonical that trip no audit flag
#: (brief_numbers ``proteoform_audit.unflagged_over_10pct_genes``). Nine were reviewed on
#: 2026-11-22 (PROTEOFORM_REVIEW_2026-11-22.md). Since the MANE-GTF rebuild of 2026-11-22,
#: TMEM145 and LY6H are built within 10% of their canonical and leave the set; SHISA6 and
#: TMEM273, newly MANE references, join it.
UNREVIEWED = {
    "KCNN2", "PLPPR2", "PCDH7", "PCDH1", "HRH2", "KCNC1",
    "GPM6B", "TSPAN17", "GFRA4", "SHISA6", "TMEM273",
}  # fmt: skip


def test_symmetric_and_bounded():
    assert LENGTH_DIFF_FRACTION == 0.10
    assert length_vs_canonical(791, 579) == (0.366, 1)  # longer
    assert length_vs_canonical(90, 100) == (0.1, 0)  # exactly 10%: not more than
    assert length_vs_canonical(269, 299) == (0.1, 1)  # GFRA4: 10.03% shorter
    assert length_vs_canonical(300, 299)[1] == 0
    assert length_vs_canonical(300, 0) == ("", 0)


def test_catches_the_unreviewed_eleven():
    with (HERE / "l1_proteoform_audit.tsv").open() as fh:
        rows = [
            r
            for r in csv.DictReader(fh, delimiter="\t")
            if r["gene_symbol"] in UNREVIEWED and r["candidate_set"] == "library_reference"
        ]
    flagged = {
        r["gene_symbol"]
        for r in rows
        if length_vs_canonical(int(r["protein_aa"]), int(r["canonical_aa"] or 0))[1]
    }
    assert {r["gene_symbol"] for r in rows} == UNREVIEWED
    assert flagged == UNREVIEWED
