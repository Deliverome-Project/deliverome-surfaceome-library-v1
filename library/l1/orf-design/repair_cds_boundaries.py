"""Apply explicit, checksum-bound CDS boundary decisions with a reversible edit ledger.

No member is automatically repaired. A decision table must name the original CDS
hash, action, biological rationale and supporting evidence. Adding ATG creates an
engineered protein and invalidates exact-sequence prediction reuse. Adding TAA
preserves the protein. Incomplete or uncertain N termini cannot pass as repaired
full-length proteins: prepend_ATG requires complete_except_initiator evidence.
The separate prepend_ATG_for_expression action permits an explicitly engineered
variant while keeping native completeness unresolved. An annotated non-AUG
start can be replaced by ATG without inserting an extra residue.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

from deliverome_analysis.cloning import seqtools

RULE_VERSION = "cds-boundaries-1"


def digest(value: str) -> str:
    return hashlib.sha256(value.encode()).hexdigest()


def repair(entry: dict, decision: dict) -> tuple[dict, dict]:
    raw = entry["cds"]
    if digest(raw) != decision["original_cds_sha256"]:
        raise ValueError("Original CDS checksum mismatch; refuses reapplication or changed input")
    if not raw or len(raw) % 3 or set(raw) - set("ACGT"):
        raise ValueError("Ambiguous or out-of-frame source requires sequence curation")
    original = seqtools.translate(raw)
    if "*" in original.rstrip("*") or original.endswith("**"):
        raise ValueError("Internal stops require sequence curation")
    if not decision.get("reason") or not decision.get("evidence"):
        raise ValueError("Repair requires an explicit rationale and evidence")
    action = decision["action"]
    if action == "append_TAA":
        if original.endswith("*"):
            raise ValueError("Source already has a stop")
        designed = raw + "TAA"
    elif action in ("prepend_ATG", "prepend_ATG_for_expression"):
        if raw.startswith("ATG"):
            raise ValueError("Source already starts ATG")
        if action == "prepend_ATG" and decision.get("completeness") != "complete_except_initiator":
            raise ValueError("Incomplete or uncertain N terminus requires curation")
        if not original.endswith("*"):
            raise ValueError("Missing stop requires its own explicit decision")
        designed = "ATG" + raw
    elif action == "replace_start_with_ATG":
        if raw.startswith("ATG"):
            raise ValueError("Source already starts ATG")
        if decision.get("completeness") not in (
            "documented_non_AUG_initiation",
            "annotated_non_AUG_start",
        ):
            raise ValueError("Start replacement requires documented non-AUG initiation")
        if not original.endswith("*"):
            raise ValueError("Missing stop requires its own explicit decision")
        designed = "ATG" + raw[3:]
    else:
        raise ValueError(f"Unknown repair action: {action}")
    before, after = original.rstrip("*"), seqtools.translate(designed).rstrip("*")
    expected = {
        "prepend_ATG": "M" + before,
        "prepend_ATG_for_expression": "M" + before,
        "replace_start_with_ATG": "M" + before[1:],
        "append_TAA": before,
    }[action]
    assert after == expected
    ledger = {
        **decision,
        "rule_version": RULE_VERSION,
        "original_cds": raw,
        "designed_cds": designed,
        "designed_cds_sha256": digest(designed),
        "original_protein_sha256": digest(before),
        "designed_protein_sha256": digest(after),
        "protein_changed": before != after,
        "exact_prediction_reuse": before == after,
        "native_completeness_restored": False,
        "engineering_scope": "translation initiation/termination only; structural completeness unresolved",
    }
    value = {**entry, "cds": designed, "boundary_design": ledger}
    return value, ledger


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cds", type=Path, required=True)
    ap.add_argument("--decisions", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    source = json.loads(args.cds.read_text())
    with args.decisions.open() as handle:
        decisions = list(csv.DictReader(handle, delimiter="\t"))
    if len({d["key"] for d in decisions}) != len(decisions):
        raise ValueError("Duplicate decision for a CDS key")
    output, edits = {}, []
    for decision in decisions:
        output[decision["key"]], edit = repair(source[decision["key"]], decision)
        edits.append(edit)
    args.out.mkdir(parents=True, exist_ok=False)
    (args.out / "cds_boundary_overrides.json").write_text(json.dumps(output, indent=2) + "\n")
    (args.out / "cds_boundary_edit_ledger.json").write_text(
        json.dumps(
            {
                "rule_version": RULE_VERSION,
                "input_sha256": {
                    "cds": digest(args.cds.read_text()),
                    "decisions": digest(args.decisions.read_text()),
                },
                "edits": edits,
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
