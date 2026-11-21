"""How the reference ladder's picks relate to MANE, APPRIS and the UniProt canonical sequence.

Reads the current reference record (the same key ``library_pipeline.REFERENCE_RECORD`` names)
and reports, per library gene:

* how often MANE Select and the APPRIS principal agree, disagree, or only one exists;
* the source finally built (MANE, APPRIS, UniProt canonical, named override);
* how often the chosen protein differs from the gene's UniProt canonical sequence, in which
  direction, and how often the UniProt canonical is not the translation of any Ensembl 115
  transcript of the gene (so it has no native coding sequence to build from).

Writes ``reference_vs_uniprot_summary.json``.

    uv run --locked python analysis/sprint4/s4e9-orf-design-fg-panel/reference_vs_uniprot.py
"""

from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path

from deliverome_analysis.s3 import processed

HERE = Path(__file__).resolve().parent
PREFIX = "sprint4/s4e9-orf-design-fg-panel/"
REFERENCE_RECORD = "library-design/reference-selection/2026-11-12-reference-fixes-2/"
OUT = HERE / "reference_vs_uniprot_summary.json"


def main() -> int:
    key = PREFIX + REFERENCE_RECORD + "reference_selection.tsv"
    with open(processed.download(key)) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    both = [r for r in rows if r["mane_enst"] and r["appris_enst"]]
    disagree = [r for r in both if r["decision"] != "mane_appris_agree"]
    changed = [r for r in rows if r["changes_current_reference"] == "True"]
    diff = [int(r["current_minus_selected"]) for r in changed]
    out = {
        "record": f"processed://{key}",
        "genes": len(rows),
        "ladder_decision": dict(Counter(r["decision"] for r in rows)),
        "source_built": dict(Counter(r["selected_source"] for r in rows)),
        "source_built_by_decision": {
            d: dict(Counter(r["selected_source"] for r in rows if r["decision"] == d))
            for d in sorted({r["decision"] for r in rows})
        },
        "genes_with_both_mane_and_appris": len(both),
        "mane_appris_disagree": len(disagree),
        "mane_appris_disagree_fraction": round(len(disagree) / len(both), 4),
        "disagreements_appris_principal_label": dict(Counter(r["appris_class"] for r in disagree)),
        "disagreements_appris_label_on_mane": dict(
            Counter(r["mane_appris_label"] for r in disagree)
        ),
        "genes_with_uniprot_canonical": sum(1 for r in rows if r["changes_current_reference"]),
        "reference_differs_from_uniprot_canonical": len(changed),
        "differs_by_source": dict(Counter(r["selected_source"] for r in changed)),
        "differs_by_decision": dict(Counter(r["decision"] for r in changed)),
        "canonical_longer_than_reference": sum(1 for d in diff if d > 0),
        "canonical_shorter_than_reference": sum(1 for d in diff if d < 0),
        "same_length_residues_differ": sum(1 for d in diff if d == 0),
        "uniprot_canonical_not_any_ensembl115_translation": sum(
            1 for r in rows if r["current_has_ensembl_transcript"] == "False"
        ),
        "uniprot_canonical_not_any_ensembl115_translation_by_decision": dict(
            Counter(r["decision"] for r in rows if r["current_has_ensembl_transcript"] == "False")
        ),
    }
    OUT.write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
