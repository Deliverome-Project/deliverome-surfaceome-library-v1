"""Turn NetGPI's ortholog calls into verdicts, calibrated on our own curated set.

D1's architecture attribute is annotated per human gene, and D1's ortholog topology is
recorded as ``projected_from_human_canonical`` — so an ortholog of a GPI-anchored gene
inherits an anchor nothing has measured. The within-species C-terminal test used for
isoforms cannot settle it either (human/mouse NYX, which differ in anchoring, score 0.37;
human/mouse CD55, which do not, score 0.43). NetGPI 1.1 measures it directly.

A prediction is only worth as much as its measured error rate on the data in hand, so the
same service was run on 148 curated-GPI human references and 100 single-pass references:
sensitivity 89.9%, specificity 99.0%, MCC 0.875 against the paper's 0.895. That run also
showed *where* NetGPI fails, which is what makes the ortholog calls actionable:

* Its confident negatives are nearly all right, and several apparent misses are not misses
  at all — the gene is GPI-anchored in an isoform the canonical reference is not (NCAM1 is
  anchored as NCAM-120 and transmembrane as NCAM-140; ACHE is anchored as the H subunit and
  not as the canonical T subunit; likewise BCAN, TFPI, CD58).
* Its low-confidence negatives are unreliable: every genuine false negative sat below 0.55
  (NYX 0.549, ART4 0.534, FCGR3B 0.445, CD14 0.388), and CD14 and FCGR3B are textbook
  GPI-anchored proteins.

So a negative is acted on only above CONFIDENT, and a gene NetGPI already got wrong on the
human reference is never acted on in the ortholog, whatever the confidence.

    uv run --locked python netgpi_verdicts.py
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path

# Above this, NetGPI's negatives were right on the control set; below it, every genuine
# false negative appeared. Not tuned — read off the control run's failures.
CONFIDENT = 0.55


def read_predictions(path: Path) -> dict[str, dict]:
    out = {}
    for line in path.read_text().splitlines():
        if line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 6:
            continue
        out[f[0]] = {"pred": f[2], "omega": f[3], "lik": float(f[4]), "aa": f[5]}
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--orthologs", default="netgpi_ortholog_predictions.txt")
    ap.add_argument("--control", default="netgpi_control_predictions.txt")
    ap.add_argument("--out", default="netgpi_ortholog_verdicts.tsv")
    args = ap.parse_args()

    control = read_predictions(Path(args.control))
    # Genes where NetGPI called the human reference non-anchored though curation says it is.
    unreliable = {
        k.split("_", 1)[1]
        for k, v in control.items()
        if k.startswith("POS_") and v["pred"] != "GPI-Anchored"
    }
    print(f"genes NetGPI missed on the human reference: {len(unreliable)}")

    rows = []
    for key, v in read_predictions(Path(args.orthologs)).items():
        species, gene, acc = key.split("_", 2)
        if v["pred"] == "GPI-Anchored":
            verdict, why = "anchored", "netgpi_positive"
        elif gene in unreliable:
            verdict, why = "uncertain", "netgpi_missed_this_gene_on_the_human_reference"
        elif v["lik"] < CONFIDENT:
            verdict, why = "uncertain", f"negative_below_{CONFIDENT}"
        else:
            verdict, why = "not_anchored", "netgpi_confident_negative"
        rows.append({"species": species, "gene": gene, "accession": acc,
                     "netgpi": v["pred"], "confidence": v["lik"],
                     "omega_site": v["omega"], "omega_residue": v["aa"],
                     "verdict": verdict, "basis": why})

    rows.sort(key=lambda r: (r["verdict"] != "not_anchored", r["species"], r["gene"]))
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {args.out} ({len(rows)} picks)")
    print(" ", dict(Counter(r["verdict"] for r in rows)))
    print("\nORTHOLOGS TO DROP — the inherited anchor is contradicted:")
    for r in rows:
        if r["verdict"] == "not_anchored":
            print(f"  {r['species']:6} {r['gene']:9} {r['accession']:20} conf {r['confidence']:.3f}")
    print("\nHELD FOR REVIEW — inherited anchor kept, call not trustworthy:")
    for r in rows:
        if r["verdict"] == "uncertain":
            print(f"  {r['species']:6} {r['gene']:9} {r['accession']:20} "
                  f"conf {r['confidence']:.3f}  {r['basis']}")


if __name__ == "__main__":
    main()
