"""Region-level labels for audited literature sites that the paper places in a loop, not at a residue.

Reads ``audit.tsv`` and, for every row whose verdict is ``region_only`` or
``wrong_species_numbering``, emits the human canonical residue window the paper's
wording supports. The window is a closed range of ``insert_after`` values
(``N`` means between N and N+1). Each window comes from one of three sources, named
in ``basis``:

* ``paper:at-residue``   the paper says "at residue N" without before/after: window N-1..N.
* ``uniprot:<feature>``  the paper names a loop, domain or boundary; the window is that
  UniProt feature on the human entry (``data/uniprot_features/<acc>.json``, EBI Proteins
  API), trimmed to junctions inside the feature.
* ``alignment``          an ortholog position mapped into a poorly conserved stretch; the
  window is the mapped position ±3 from the audit note.

Sites with no defensible window (secreted ZP4; the GABRA2 numbering frame) are
written with an empty window and the reason, so the ranker can skip them explicitly.

Run::

    uv run --locked python library/l1/tag-site-corpus-audit/region_labels.py \
        --audit library/l1/tag-site-corpus-audit/audit.tsv \
        --features library/l1/tag-site-corpus-audit/data/uniprot_features \
        --out library/l1/tag-site-corpus-audit/region_labels.tsv
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

# site_id -> (basis, window or feature spec, reason). A feature spec is
# ("TOPO_DOM"|"DOMAIN"|"PROPEP", begin) naming the UniProt feature that starts at `begin`;
# ("between", end_of_first_domain, begin_of_second_domain) is the linker between two.
RULES: dict[str, tuple] = {
    "ABCB1-internal-None-lit": ("uniprot:TOPO_DOM", ("TOPO_DOM", 68), "HiBiT in 'the 1st extracellular loop' (ECL1)"),
    "AQP1-internal-120-lit": ("paper:at-residue", 120, "'epitope at Thr-120', before/after not stated"),
    "ATP1A3-internal-None-lit": ("uniprot:TOPO_DOM", ("TOPO_DOM", 304), "paralog: rat alpha1 318 is in ECL2; alpha3 ECL2 by UniProt"),
    "ATP1A3-internal-109-lit": ("uniprot:TOPO_DOM", ("TOPO_DOM", 108), "paralog: rat alpha1 119 is in ECL1; alpha3 ECL1 by UniProt"),
    "BACE1-internal-45-lit": ("uniprot:PROPEP", ("PROPEP", 22), "'after the propeptide cleavage site': junction at the propeptide end (45)"),
    "BST2-internal-146-lit": ("paper:at-residue", 146, "'Flag tag at residue 146'"),
    "CD1D-internal-301-lit": ("uniprot:boundary", ("TOPO_DOM", 20), "'between the extracellular and the transmembrane domain': last three junctions of the ectodomain"),
    "CFTR-internal-897-lit": ("paper:at-residue", 898, "'3xHA positioned at position S898' in an isolated TMD2 fragment; no full-length construct"),
    "CLCA4-internal-None-lit": ("uniprot:REGION", ("REGION", 870), "'between the FnIII and GPI anchor domains': disordered stalk before the C-terminal TM"),
    "CLCA4-internal-594-lit": ("paper:at-residue", 595, "nucleotide 1782 after the start codon: codon 594 or 595"),
    "GABRA2-internal-212-lit": ("none", None, "S181TAG in an unstated numbering frame; canonical candidates 209 (mature) or 213 (+32)"),
    "GABRA2-internal-232-lit": ("none", None, "S201TAG in an unstated numbering frame; canonical candidates 229 (mature) or 233 (+32)"),
    "GRIA1-internal-396-lit": ("paper:at-residue", 396, "'at position 396 aa' in a transposon construct"),
    "KCNJ1-internal-113-lit": ("paper:at-residue", 113, "'inserted at Y113' (rat ROMK; same number in human)"),
    "MMP14-internal-111-lit": ("uniprot:PROPEP", ("PROPEP", 21), "'between the propeptide and catalytic domain': junction at the propeptide end (111)"),
    "NRG1-internal-None-lit": ("uniprot:between", ("between", 128, 178), "'between Ig and EGF domains': Ig-like C2 ends 128, EGF-like begins 178"),
    "PMP22-internal-None-lit": ("uniprot:TOPO_DOM", ("TOPO_DOM", 120), "'second extracellular loop'"),
    "SLC46A2-internal-57-lit": ("uniprot:TOPO_DOM", ("TOPO_DOM", 45), "'first predicted extracellular loop'; AfeI site in a codon-optimized CDS"),
    "SLC6A3-internal-None-lit": ("uniprot:TOPO_DOM", ("TOPO_DOM", 172), "EL2 (from the row's second source; the cited paper names no loop)"),
    "SLC6A3-internal-199-lit": ("paper:at-residue", 200, "junction written G199/S201: after 199 or 200"),
    "TMEM123-internal-34-lit": ("paper:at-residue", 34, "'following the signal peptide (Alanine 34)': before or after A34"),
    "TRIL-internal-499-lit": ("alignment", 474, "Xenopus 499|500 maps to human ~473-475"),
    "ADGRG6-internal-550-lit": ("alignment", 567, "zebrafish Y551 maps near human 566-568"),
    "ZP4-internal-468-lit": ("none", None, "porcine secreted ZP4 construct; no window on the human membrane form"),
}
COLS = ["site_id", "gene", "uniprot_acc", "verdict", "window_start", "window_end", "window_width",
        "basis", "reason", "surface_shown", "function_shown", "pmid", "in_training"]


def feature_window(features: list[dict], spec: tuple) -> tuple[int, int]:
    kind = spec[0]
    if kind == "between":
        return spec[1], spec[2] - 1
    hits = [f for f in features if f["type"] == kind and int(f["begin"]) == spec[1]]
    if not hits:
        raise SystemExit(f"no {kind} feature starting at {spec[1]}")
    begin, end = int(hits[0]["begin"]), int(hits[0]["end"])
    if kind == "PROPEP":
        return end, end  # the junction right after the propeptide
    return begin, end - 1  # junctions strictly inside the loop/region


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--audit", required=True, type=Path)
    ap.add_argument("--features", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    out = []
    with args.audit.open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["verdict"] not in ("region_only", "wrong_species_numbering"):
                continue
            rule = RULES.get(r["site_id"])
            if rule is None:
                raise SystemExit(f"no region rule for {r['site_id']}; add one to RULES")
            basis, spec, reason = rule
            if basis == "none":
                start = end = ""
            elif basis == "paper:at-residue":
                start, end = spec - 1, spec
            elif basis == "alignment":
                start, end = spec - 3, spec + 3
            elif basis == "uniprot:boundary":
                feats = json.loads((args.features / f"{r['uniprot_acc']}.json").read_text())["features"]
                _, loop_end = feature_window(feats, spec)
                start, end = loop_end - 1, loop_end + 1  # loop_end is already end-1
            else:
                feats = json.loads((args.features / f"{r['uniprot_acc']}.json").read_text())["features"]
                start, end = feature_window(feats, spec)
            out.append({"site_id": r["site_id"], "gene": r["gene"], "uniprot_acc": r["uniprot_acc"],
                        "verdict": r["verdict"], "window_start": start, "window_end": end,
                        "window_width": (end - start + 1) if start != "" else "",
                        "basis": basis, "reason": reason, "surface_shown": r["surface_shown"],
                        "function_shown": r["function_shown"], "pmid": r["pmid"],
                        "in_training": r["in_training"]})
    with args.out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out)
    usable = [o for o in out if o["window_start"] != ""]
    print(f"{len(out)} region-level sites, {len(usable)} with a window "
          f"({len({o['uniprot_acc'] for o in usable})} proteins); widths: "
          f"{sorted(o['window_width'] for o in usable)}")


if __name__ == "__main__":
    main()
