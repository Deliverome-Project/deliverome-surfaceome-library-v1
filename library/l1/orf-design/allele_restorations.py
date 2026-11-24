"""Designed CDSs for library genes whose reference genome copy is broken but a working copy exists.

Six library genes had no designed CDS because every Ensembl 115 transcript of them is
inactivated: GRCh38 carries a stop codon or a frameshift (``genes_without_reference_cds.tsv``).
For four of them a working copy is documented, and Becca decided on 2026-11-22 to build them
("build or and fcgr2c"; OR10AC1 and OR52P1 because most people carry the intact allele):

* ``screened_patch_form``: OR4Q2. The GPCR screen (Tedman et al.) displayed ENST00000710368, the
  patch-scaffold copy (313 aa, one residue short of UniProt's 314), at 1,116 a.u. It is built as
  that transcript, unedited.
* ``common_intact_allele``: a pinned Ensembl 115 transcript with the recorded variants that
  restore the intact allele (each with its rsID and population frequency in the ledger):
  FCGR2C (codon 57 TAG -> CAG; the ORF allele, van der Heijden et al. 2012, PMID 22198951),
  OR10AC1 (one G of a 4-G run removed; intact allele 73%), OR52P1 (one T deleted and GAGC
  inserted, in near-complete linkage; intact allele about 60%).

Each edit is applied to the pinned CDS and the result must translate to the protein the ledger
names (``expected``: a UniProt accession, or ``transcript`` for an unedited form), within the
same 3-residue reference-allele tolerance CDS resolution uses; the residues that still differ
are recorded. OR8J2 (intact allele about 6%) and SIGLEC16 (a restored frame that matches no
UniProt record, so it would need a synthetic CDS) are held, as open decisions in the brief.

The ledger is ``allele_restorations.tsv``; this script writes ``allele_restorations_cds.json``
(``{"lib:<gene>": {cds, transcript, route, edits, ...}}``), which ``add_library_genes.py`` takes
for a gene the resolver leaves unresolved.

Implemented by Claude Code session l1-review-export-scope, 2026-11-22.

    uv run --locked python allele_restorations.py            # build and write
    uv run --locked python allele_restorations.py --check    # rebuild, compare with the record
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path

HERE = Path(__file__).parent
AP = HERE.parents[2] / "library/l1/proteoform-selection"
sys.path.insert(0, str(AP))
import candidate_inputs  # noqa: E402
import resolve_candidate_cds as resolver  # noqa: E402

LEDGER = HERE / "allele_restorations.tsv"
OUT = HERE / "allele_restorations_cds.json"
TOLERANCE = resolver.VARIANT_TOLERANCE


def apply_edits(cds: str, edits: list[dict]) -> str:
    """Apply ``{"pos": 0-based, "ref": str, "alt": str}`` edits, rightmost first."""
    for e in sorted(edits, key=lambda e: -int(e["pos"])):
        pos, ref, alt = int(e["pos"]), e["ref"], e["alt"]
        assert cds[pos : pos + len(ref)] == ref, (
            f"expected {ref!r} at {pos}, found {cds[pos : pos + len(ref)]!r}"
        )
        cds = cds[:pos] + alt + cds[pos + len(ref) :]
    return cds


def build() -> dict[str, dict]:
    cds_all = resolver.read_fasta(candidate_inputs.cds_fasta("human"))
    uniprot = resolver.load_uniprot("human")
    out: dict[str, dict] = {}
    with LEDGER.open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    for r in rows:
        base = cds_all[r["transcript"]][1].upper()
        edits = json.loads(r["edits"])
        cds = apply_edits(base, edits)
        assert len(cds) % 3 == 0, f"{r['gene']}: edited CDS is not a whole number of codons"
        protein = resolver.translate(cds)
        assert protein.endswith("*") and protein.count("*") == 1, (
            f"{r['gene']}: not one stop, at the end"
        )
        protein = protein[:-1]
        if r["expected"] == "transcript":
            want = resolver.translate(base).rstrip("*")
        else:
            want = uniprot[r["expected"]]
        diffs = [
            f"{a}{i + 1}{b}" for i, (a, b) in enumerate(zip(want, protein, strict=False)) if a != b
        ]
        if r["expected"] != "transcript":
            assert len(protein) == len(want), f"{r['gene']}: {len(protein)} aa vs {len(want)}"
        assert len(diffs) <= TOLERANCE + int(r.get("extra_tolerance") or 0), (
            f"{r['gene']}: {len(diffs)} residues differ from {r['expected']}"
        )
        out[f"lib:{r['gene']}"] = {
            "cds": cds,
            "transcript": r["transcript"],
            "route": r["basis"],
            "edits": edits,
            "expected": r["expected"],
            "residues_differing": diffs,
            "protein_sha256": hashlib.sha256(protein.encode()).hexdigest(),
        }
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args(argv)
    built = build()
    for k, v in built.items():
        print(
            f"{k}: {v['transcript']} {v['route']}, {len(v['cds']) // 3 - 1} aa, "
            f"differs from {v['expected']} at {v['residues_differing'] or 'none'}"
        )
    text = json.dumps(built, indent=1, sort_keys=True) + "\n"
    if a.check:
        same = OUT.exists() and OUT.read_text() == text
        print("matches the record" if same else "DIFFERS from the record")
        return 0 if same else 1
    OUT.write_text(text)
    print(f"wrote {OUT.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
