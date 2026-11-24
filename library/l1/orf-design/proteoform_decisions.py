"""The proteoform decisions this experiment applies, with their attribution.

``proteoform_decisions_applied.tsv`` holds, verbatim, the 14 manifest-changing rows of
``l1_proteoform_decisions.tsv`` (session l1-proteoform-audit, PR #300, commit 8fbcd342;
SHA-256 in each row's ``source_record``): decided by l1-proteoform-audit under Becca's
2026-11-12 delegation ("use your best judgment"), applied here by the session
l1-reviewer-package (subagent), PR #283. The scripts that apply them
(``apply_reference_fixes.py``, ``redraw_panels.py``, ``adopt_ortholog_cds_policy.py``) copy
each row's ``decided_by`` / ``implemented_by`` / ``reason`` / ``alternative_rejected`` /
``sources`` into their own ledgers from here, so the attribution is never retyped.
"""

from __future__ import annotations

import csv
from pathlib import Path

TABLE = Path(__file__).resolve().parent / "proteoform_decisions_applied.tsv"
APPLIER = "applied by l1-reviewer-package (subagent), PR #283"

#: Rows whose decision was reconciled with another session's before it was applied: their
#: attribution and alternatives replace the audit row's, and the evidence is appended.
RECONCILED = {
    "lib:LMLN": {
        "decided_by": "l1-reviewer-package (delegated by Becca 2026-11-12), reconciling "
        "l1-proteoform-audit (canonical) and l1-nterm-type3 (Q96KR4-3, #282 f980ce49)",
        "implemented_by": "Claude Code session l1-reviewer-package (subagent), PR #283",
        "reason_suffix": " Reconciled choice: the UniProt canonical Q96KR4 via ENST00000330198, "
        "whose CDS starts ATG and encodes Q96KR4 exactly (the fallback, Q96KR4-3 already built "
        "in L2 as l2-p0894, was not needed). Both forms keep the same 43-aa signal peptide "
        "(DeepTMHMM v1 SP 1-43, identical N-termini); NetGPI calls both GPI-anchored (canonical "
        "omega 620, isoform-3-like forms omega ~656, likelihood 0.461 each); isoform 3 differs "
        "only by a 37-aa insertion after residue 403 with no evidence for its function, so the "
        "canonical is preferred.",
        "alternative": "Q96KR4-3 (l2-p0894, ENST00000420910, 684 aa); keep the 452-aa "
        "downstream ENST00000482695 form",
    }
}


def by_key() -> dict[str, dict]:
    with TABLE.open() as fh:
        return {r["key"]: r for r in csv.DictReader(fh, delimiter="\t")}


def attribution(key: str) -> dict[str, str]:
    """The ledger fields for one applied decision."""
    r = by_key()[key]
    out = {
        "decided_by": r["decided_by"],
        "implemented_by": f"{r['implemented_by']}; {r['applied_by']}",
        "reason": r["reason"],
        "alternative": r["alternative_rejected"],
        "sources": f"{r['sources']}; {r['source_record']}",
    }
    rec = RECONCILED.get(key)
    if rec:
        out |= {
            "decided_by": rec["decided_by"],
            "implemented_by": rec["implemented_by"],
            "reason": r["reason"] + rec["reason_suffix"],
            "alternative": rec["alternative"],
        }
    return out
