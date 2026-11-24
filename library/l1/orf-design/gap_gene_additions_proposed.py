"""Proposed L1 additions from the surface-proteomics challenge set, in gene_list_additions format.

Becca, on the reviewer brief: add the challenge-set genes that pass the addition criteria
(2026-11-20, "add to the library and update this section"), and FAM234B, MYORG and TMEM94, which
pass every criterion except a co-receptor check that cannot run without a catalog record
(2026-11-22, "in that case we want to add them"). This writes those 18 genes as rows in the
columns of ``gene_list_additions.tsv`` (``added`` left blank: the gene list itself is changed by
the L1 rebuild, not here) to ``gap_gene_additions_proposed.tsv``.

The genes and their values come from ``surface_proteomics_gap_candidates.tsv`` (verdict
``candidate`` or ``candidate pending co-receptor review``). HGNC ids are the HGNC
cross-reference of each accession's UniProt entry, fetched from rest.uniprot.org; the run fails
if an entry has no HGNC id or more than one, or if UniProt's gene name disagrees with the table.

The three genes without a catalog record had their co-receptor requirement reviewed by hand
from UniProt's subunit, interaction and location annotations (release 2026_04, read
2026-11-22). That review is recorded in ``MANUAL_CORECEPTOR`` below with what it rests on, and
copied into each row's evidence.

    uv run --locked python library/l1/orf-design/gap_gene_additions_proposed.py
"""

from __future__ import annotations

import csv
import json
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
CANDIDATES = HERE / "surface_proteomics_gap_candidates.tsv"
OUT = HERE / "gap_gene_additions_proposed.tsv"
COLUMNS = ("gene", "hgnc_id", "uniprot", "added", "decided_by", "source", "evidence")
DECIDED_BY = "Becca Carlson"
REVIEWER = "agent session surface-gate (L1 brief review), 2026-11-22"

#: Hand review of the co-receptor requirement for the genes with no catalog record.
MANUAL_CORECEPTOR = {
    "A2RU67": (
        "no heterologous obligate partner: UniProt 2026_04 gives no subunit annotation and one "
        "binary interactor (ABCD4, 5 experiments), with no requirement for surface delivery; "
        "location single-pass membrane, Golgi outpost in neurons (by similarity to rat)"
    ),
    "Q6NSJ0": (
        "no heterologous obligate partner: UniProt 2026_04 subunit = homodimer (PubMed:36129849), "
        "which the library treats as self-satisfied; it binds IGF2 for IGF2's secretion (IGF2 is "
        "its cargo, not its trafficking partner). Caveat: UniProt places it in the ER and nuclear "
        "membranes, so surface display is uncertain"
    ),
    "Q12767": (
        "no heterologous obligate partner: UniProt 2026_04 subunit = homo-oligomers "
        "(PubMed:38513662), self-satisfied. Caveat: UniProt places it in the ER membrane (ER "
        "Mg2+ uptake), so surface display is uncertain"
    ),
}


def uniprot_hgnc(acc: str) -> tuple[str, str, str]:
    url = f"https://rest.uniprot.org/uniprotkb/{acc}.json"
    req = urllib.request.Request(url, headers={"User-Agent": "deliverome-analysis"})
    with urllib.request.urlopen(req, timeout=60) as resp:
        release = resp.headers.get("X-UniProt-Release", "")
        entry = json.load(resp)
    hgnc = [x["id"] for x in entry.get("uniProtKBCrossReferences", []) if x["database"] == "HGNC"]
    if len(hgnc) != 1:
        raise SystemExit(f"{acc}: {len(hgnc)} HGNC cross-references")
    name = (entry.get("genes") or [{}])[0].get("geneName", {}).get("value", "")
    return hgnc[0], name, release


def main() -> None:
    with CANDIDATES.open() as fh:
        cands = [
            r
            for r in csv.DictReader(fh, delimiter="\t")
            if r["verdict"] in ("candidate", "candidate pending co-receptor review")
        ]
    assert len(cands) == 18, f"{len(cands)} candidates, expected 15 + 3"
    rows, releases = [], set()
    for r in sorted(cands, key=lambda r: r["gene"]):
        acc = r["uniprot_accession"]
        hgnc, name, release = uniprot_hgnc(acc)
        releases.add(release)
        gene = name  # HGNC-approved symbol; the table carries the sweep's mnemonic for 2 rows
        if r["catalog_tier"] != "no deep-dive record" and r["gene"] != name:
            raise SystemExit(f"{acc}: table says {r['gene']}, UniProt says {name}")
        anchor = f"{r['v1_tm_helices']} TM helices" + (
            " + catalog GPI" if r["catalog_gpi"] == "True" else ""
        )
        if acc in MANUAL_CORECEPTOR:
            co = f"co-receptor reviewed by hand ({REVIEWER}): {MANUAL_CORECEPTOR[acc]}"
            waived = "no catalog deep-dive record, so the catalog gates were never evaluated"
        else:
            co = "passes the library co-receptor gate (catalog calls + reviewed overrides)"
            waived = (
                f"catalog tier {r['catalog_tier'] or 'none'}; first S1E3 gate failed: "
                f"{r['s1e3_first_failing_gate']} (waived for challenge-set genes)"
            )
        evidence = (
            f"surface-proteomics challenge set, rule B on sweep-v20 + abundance-v5: "
            f"{r['n_surface_pubs_current']} surface-labeling publications, lod_membrane "
            f"{r['lod_membrane_current']}; ORF {r['orf_nt']} nt (cap 5,000); DeepTMHMM v1 "
            f"({r['topology_version']}) on the UniProt canonical: {anchor}, "
            f"{r['v1_extracellular_res']} extracellular residues (min 30); {co}; {waived}"
        )
        rows.append(
            {
                "gene": gene,
                "hgnc_id": hgnc,
                "uniprot": acc,
                "added": "",
                "decided_by": DECIDED_BY,
                "source": (
                    "reviewer-brief comments 2026-11-20 and 2026-11-22; "
                    "surface_proteomics_gap_candidates.tsv (inputs pinned in "
                    "surface_proteomics_inputs.json)"
                ),
                "evidence": evidence,
            }
        )
    rows.sort(key=lambda r: r["gene"])
    with OUT.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(f"{OUT.name}: {len(rows)} genes; UniProt release {', '.join(sorted(releases))}")
    for r in rows:
        print(f"  {r['gene']:<9} {r['hgnc_id']:<11} {r['uniprot']}")


if __name__ == "__main__":
    main()
