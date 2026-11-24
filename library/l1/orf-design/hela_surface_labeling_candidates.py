"""Where the HeLa surface-labeling off-list candidates fall against the library's gates.

Input: ``hela_surface_labeling_class_a.tsv``, the 16 "Class A" rows of Bobby's off-list
candidate table (proteins that NHS-SS-biotin bead eluates enrich at least 4-fold over the input,
at q < 0.01, in at least 4 of 6 beads, are not canonical or likely surface in the catalog, and
are membrane-anchored with plasma-membrane evidence). Copied, columns reduced, from
``analysis/sprint6/s6e3-streptavidin-bead-ms/s6e3_offlist_surface_candidates.tsv`` on branch
``bobby/s6e3-tic-load`` at commit ``b1615282`` (2026-11-16; SHA-256 of the full file
``eed99fcb…``), which has not merged. The copy keeps the source's columns verbatim.

For each accession (base UniProt accession; the EVA1A;EVA1B protein group is scored as its
lead accession Q9H8M9, EVA1A) it reports, using ``surface_proteomics_gap_candidates`` so the
tests are identical to the surface-proteomics challenge set:

- whether it is already an L1 gene;
- the first S1E3 gate that excludes it, with the catalog fields that decide that gate;
- whether it passes the addition criteria (ORF <= 5,000 nt, the library co-receptor gate,
  DeepTMHMM v1 TM helix or GPI anchor, >= 30 extracellular residues);
- whether it is also in the surface-proteomics challenge set (current pinned snapshot).

    uv run --locked python library/l1/orf-design/hela_surface_labeling_candidates.py
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import surface_proteomics_gap_candidates as G  # noqa: E402

IN = HERE / "hela_surface_labeling_class_a.tsv"
OUT = HERE / "hela_surface_labeling_candidates.tsv"
KEEP = (
    "uniprot_accession",
    "gene",
    "in_l1",
    "s1e3_first_failing_gate",
    "catalog_tier",
    "catalog_surface_specificity",
    "catalog_surface_call_reason",
    "catalog_primary_compartment",
    "verdict",
    "fails",
    "orf_nt",
    "v1_tm_helices",
    "v1_extracellular_res",
    "coreceptor_test",
    "coreceptor_partners",
    "in_surface_proteomics_gap",
    "bobby_rank",
    "bobby_label",
    "beads_enriched",
    "median_log2_enrichment",
    "uniprot_location",
)


def main() -> None:
    with IN.open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    items = {
        r["uniprot_accession"].split(";")[0].split("-")[0]: {"sweep_symbol": r["symbol"], **r}
        for r in rows
    }
    gap = G.gap_accessions()
    scored, _ = G.score_accessions(
        items,
        lambda g: {
            "in_surface_proteomics_gap": g["uniprot_accession"] in gap,
            "bobby_rank": g["rank"],
            "bobby_label": g["symbol"],
            "beads_enriched": g["beads_enriched"],
            "median_log2_enrichment": g["median_log2_enrichment"],
            "uniprot_location": g["uniprot_location"],
        },
    )
    out = [{k: r.get(k, "") for k in KEEP} for r in scored]
    out.sort(key=lambda r: int(r["bobby_rank"]))
    G.write(out, OUT)
    for r in out:
        print(
            f"{r['gene']:<9} in_l1={r['in_l1'] or '-':<6} gate={r['s1e3_first_failing_gate'] or 'no record':<17}"
            f" spec={r['catalog_surface_specificity']:<21} reason={r['catalog_surface_call_reason']:<24}"
            f" {r['verdict']:<36} fails={r['fails']:<12} gap={r['in_surface_proteomics_gap']}"
        )


if __name__ == "__main__":
    main()
