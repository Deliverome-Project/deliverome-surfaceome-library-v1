"""Which library gate each S3E8 surface-proteomics gap candidate fails.

S3E8's cross-method sweep finds surface proteins that its membrane-fractionation
and whole-cell arms both detect; a subset of those are in neither L1 nor L2. A
reviewer asking "why is this protein not in the library?" needs a per-gene answer,
and this script writes it.

The gates are the library's own, from
``analysis/sprint1/s1e3-surface-library-filters/final_library.py`` and described in
section 2 of the reviewer brief. Nothing here re-decides membership: it reports, per
candidate, which gate the surfaceome record fails.

Inputs, both committed so the table regenerates offline:

``s3e8_gap_evidence.json``
    The gate fields of each candidate's surfaceome record, reduced from the public
    API (``/v1/catalog`` and ``/v1/genes/{symbol}``). Every record's
    ``gene.uniprot_acc`` was checked against the accession it was fetched for.

``s3e8_gap_candidates.tsv``
    The candidate list, copied from the S3E8 sweep output named in the brief.

Keyed on ``uniprot_accession``, never the display symbol: the sweep output carries
UniProt *entry-name mnemonics* for three rows (``TMM94``, ``F234B``, ``K1549``,
which truncate the gene name to five characters), so a symbol join would silently
drop them. The ``gene`` column here is the catalog's HGNC symbol.

Run from the repo root:

    uv run --locked python analysis/sprint4/s4e9-orf-design-fg-panel/s3e8_gap_gate_reasons.py
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
EVIDENCE = HERE / "s3e8_gap_evidence.json"
CANDIDATES = HERE / "s3e8_gap_candidates.tsv"
OUT = HERE / "s3e8_gap_gate_reasons.tsv"

#: Tiers the library accepts as surface evidence (canonical union likely).
SURFACE_TIERS = ("canonical", "likely")

#: Localization calls that are not a plasma-membrane protein.
NON_SURFACE_LOCALIZATION = ("secreted_only", "endomembrane_resident")

#: The synthesis cap; longer ORFs are routed to the L2 truncation design.
ORF_NT_CAP = 5000


def orf_nt(protein_length: int) -> int:
    """Coding length including the stop codon."""
    return (protein_length + 1) * 3


def gates_failed(ev: dict) -> tuple[list[str], list[str]]:
    """Return (gate names, human-readable details) for one candidate."""
    rec = ev["deep_dive_record"]
    gates: list[str] = []
    detail: list[str] = []

    if rec is None:
        # No deep-dive record exists. The candidate never passed triage into the
        # deep-dive set, so the triage verdict is the whole story.
        triage = ev.get("catalog_triage") or []
        calls = [t for t in triage if isinstance(t, list) and t]
        verdict = calls[0][0] if calls else None
        reason = calls[0][1] if calls and len(calls[0]) > 1 else None
        gates.append("surface evidence")
        bits = [f"no deep-dive record; triage verdict={verdict or 'none'}"]
        if reason:
            bits.append(f"triage reason={reason}")
        bits.append(
            f"{ev.get('catalog_db_votes', 0)} DB votes, "
            f"{ev.get('catalog_n_sources', 0)} literature sources"
        )
        detail.append("; ".join(bits))
        return gates, detail

    tier = rec.get("deep_dive_tier")
    length = rec.get("protein_length")

    if length is not None and orf_nt(length) > ORF_NT_CAP:
        gates.append("ORF size")
        detail.append(
            f"ORF {orf_nt(length)} nt over the {ORF_NT_CAP:,} nt cap; routed to L2 truncation"
        )
    if tier not in SURFACE_TIERS:
        gates.append("surface evidence")
        detail.append(f"deep_dive_tier={tier}")
    if rec.get("surface_specificity") == "mostly_intracellular":
        gates.append("specificity")
        detail.append("surface_specificity=mostly_intracellular")
    if rec.get("surface_call_reason") in NON_SURFACE_LOCALIZATION:
        gates.append("localization")
        detail.append(f"surface_call_reason={rec['surface_call_reason']}")
    anchored = (rec.get("tm_helix_count") or 0) >= 1
    if not anchored:
        gates.append("anchored")
        detail.append("tm_helix_count=0, not GPI-anchored")
    if rec.get("surface_expression_dependency") == "required" or rec.get(
        "requires_coreceptor_for_expression"
    ):
        gates.append("co-receptor")
        detail.append("obligate partner required")
    if rec.get("low_endogenous_expression"):
        gates.append("expression")
        detail.append("low endogenous expression")

    return gates, detail


def main() -> None:
    evidence = json.loads(EVIDENCE.read_text())
    provenance = evidence.get("_provenance", {})
    candidates = list(csv.DictReader(CANDIDATES.open(), delimiter="\t"))

    rows = []
    for cand in candidates:
        acc = cand["uniprot_accession"]
        ev = evidence[acc]
        rec = ev["deep_dive_record"] or {}
        gates, detail = gates_failed(ev)
        length = rec.get("protein_length")
        rows.append(
            {
                "uniprot_accession": acc,
                "gene": ev["gene_symbol"],
                "surface_tier": rec.get("deep_dive_tier") or "none",
                "length_aa": length if length is not None else "",
                "orf_nt": orf_nt(length) if length is not None else "",
                "gate_failed": "; ".join(gates),
                "detail": "; ".join(detail),
            }
        )

    rows.sort(key=lambda r: r["gene"])
    with OUT.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    unexplained = [r["gene"] for r in rows if not r["gate_failed"]]
    print(f"{OUT.name}: {len(rows)} candidates")
    print(f"  catalog universe: {provenance.get('universe_version')}")
    print(f"  with a stated gate: {len(rows) - len(unexplained)}")
    if unexplained:
        raise SystemExit(f"candidates with no failing gate: {unexplained}")


if __name__ == "__main__":
    main()
