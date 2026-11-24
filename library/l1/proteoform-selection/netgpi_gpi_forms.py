"""Locate the GPI anchor on a specific form of each anchored gene.

The D1 architecture attribute is a claim about the gene: *some* form of it is GPI-anchored.
It does not say the form the library selected is that one, and for a sixth of them it is
not -- NCAM1 is anchored as NCAM-120 and transmembrane as NCAM-140, ACHE as the H subunit
and not the canonical T subunit, and likewise BCAN, TFPI and CD58.

So NetGPI 1.1 is run over every human form of every GPI gene -- the selected reference and
all 770 alternate forms, not just the evidence-filtered cohort, because an anchor sitting on
a form the cohort dropped is exactly what needs finding -- plus the ortholog picks.

Reading the calls uses what the calibration run measured (148 curated-GPI against 100
single-pass human references: sensitivity 89.9%, specificity 99.0%, MCC 0.875 against the
paper's 0.895):

* A negative below CONFIDENT is not acted on. Every genuine false negative on the control
  sat below it (NYX 0.549, ART4 0.534, FCGR3B 0.445, CD14 0.388), and CD14 and FCGR3B are
  textbook GPI-anchored proteins.
* A form shorter than MIN_LENGTH cannot carry a signal that is 18-39 residues by itself, so
  its negative is uninterpretable rather than evidence. These are Ensembl fragments.

A confidence cutoff alone is NOT enough to act on a negative, which this run showed and the
calibration run had not: NetGPI calls every one of CRIPTO's ten forms non-anchored at 0.745,
and Cripto-1 is GPI-anchored (PMID 21055389). CFC1, its EGF-CFC family partner, is missed at
0.93, and GML -- whose HGNC name is literally "GPI anchored molecule like" -- at 0.752. So a
gene with no anchored form is REVIEWED, never acted on: it is equally likely that the gene
attribute is wrong and that NetGPI has a family-shaped blind spot.

A gene whose anchor is found on a different form is the reliable half of the output. All six
are textbook: NCAM-120 against NCAM-140, the ACHE H subunit against the canonical T subunit,
TFPI-beta against TFPI-alpha, and the GPI forms of BCAN, CD55 and CD58.

    uv run --locked python netgpi_gpi_forms.py
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path

from library import library_members

CONFIDENT = 0.55
MIN_LENGTH = 50


def read(path: Path) -> list[dict]:
    out = []
    for line in path.read_text().splitlines():
        if line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 6:
            continue
        out.append({"id": f[0], "length": int(f[1]), "pred": f[2],
                    "omega": f[3], "lik": float(f[4]), "aa": f[5]})
    return out


def classify(r: dict) -> str:
    if r["pred"] == "GPI-Anchored":
        return "anchored"
    if r["length"] < MIN_LENGTH:
        return "too_short"
    if r["lik"] < CONFIDENT:
        return "unresolved"
    return "not_anchored"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--human", default="netgpi_gpi_human_predictions.txt")
    ap.add_argument("--orthologs", default="netgpi_ortholog_predictions.txt")
    ap.add_argument("--out", default="netgpi_gpi_form_calls.tsv")
    ap.add_argument("--summary", default="netgpi_gpi_gene_summary.tsv")
    args = ap.parse_args()

    # The call table is a snapshot; a gene dropped from the library must not reappear here.
    members = library_members()
    forms: dict[str, list[dict]] = defaultdict(list)
    for r in read(Path(args.human)):
        kind, rest = r["id"].split("_", 1)
        gene, tag = rest.rsplit("_", 1)
        if members and gene not in members:
            continue
        forms[gene].append({**r, "kind": kind, "gene": gene, "form": tag,
                            "call": classify(r)})
    for r in read(Path(args.orthologs)):
        species, rest = r["id"].split("_", 1)
        gene, tag = rest.rsplit("_", 1)
        if members and gene not in members:
            continue
        forms[gene].append({**r, "kind": species, "gene": gene, "form": tag,
                            "call": classify(r)})

    rows = [f for g in sorted(forms) for f in forms[g]]
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t",
                           fieldnames=["gene", "kind", "form", "length", "pred",
                                       "lik", "omega", "aa", "call", "id"])
        w.writeheader()
        for r in rows:
            w.writerow({k: r[k] for k in w.fieldnames})
    print(f"wrote {args.out} ({len(rows)} forms over {len(forms)} genes)")
    print("  calls:", dict(Counter(r["call"] for r in rows)))

    summary, verdicts = [], Counter()
    for gene in sorted(forms):
        fs = forms[gene]
        ref = next((f for f in fs if f["kind"] == "ref"), None)
        human = [f for f in fs if f["kind"] in ("ref", "alt")]
        anchored = [f for f in human if f["call"] == "anchored"]
        if ref is None:
            verdict = "no_reference_submitted"
        elif ref["call"] == "anchored":
            verdict = "reference_is_anchored"
        elif anchored:
            verdict = "anchor_on_another_form"
        elif any(f["call"] == "unresolved" for f in human):
            verdict = "unresolved"
        else:
            verdict = "no_anchored_form_found"
        verdicts[verdict] += 1
        summary.append({
            "gene": gene, "verdict": verdict,
            "reference_form": ref["form"] if ref else "",
            "reference_call": ref["call"] if ref else "",
            "reference_conf": ref["lik"] if ref else "",
            "n_human_forms": len(human),
            "n_anchored_human_forms": len(anchored),
            "anchored_forms": ";".join(f["form"] for f in anchored[:6]),
            "mouse_call": next((f["call"] for f in fs if f["kind"] == "mouse"), ""),
            "cyno_call": next((f["call"] for f in fs if f["kind"] == "cyno"), ""),
        })
    with open(args.summary, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(summary[0]), delimiter="\t")
        w.writeheader()
        w.writerows(summary)
    print(f"\nwrote {args.summary} ({len(summary)} genes)")
    for k, v in verdicts.most_common():
        print(f"  {k:26} {v}")

    print("\n=== THE LIBRARY IS CARRYING THE WRONG FORM — anchor is on another form ===")
    for s in summary:
        if s["verdict"] == "anchor_on_another_form":
            print(f"  {s['gene']:10} reference {s['reference_form']:14} "
                  f"({s['reference_call']}, {s['reference_conf']}) -> anchored: "
                  f"{s['anchored_forms']}")

    print("\n=== NO ANCHORED FORM FOUND — the gene-level GPI call looks wrong ===")
    for s in summary:
        if s["verdict"] == "no_anchored_form_found":
            print(f"  {s['gene']:10} {s['n_human_forms']} human forms, reference "
                  f"{s['reference_call']} at {s['reference_conf']}")


if __name__ == "__main__":
    main()
