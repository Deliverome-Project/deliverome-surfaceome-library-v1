"""Per-gene constraints on where a tag may sit, for the whole library.

Tag placement is still an open decision at library scale, and several
independent findings have now stacked constraints onto it. This gathers them
into one table so the decision is made against all of them at once rather than
rediscovering each one gene by gene.

Four rules, each from a different source:

* A GPI gene's C-terminal residues are removed entirely at the omega site, so a
  C-terminal tag is cleaved off the mature protein. This is why FOLR1 takes an
  N-terminal tag on the panel.
* A sheddase substrate loses its distal ectodomain, taking an N-terminal tag
  with it -- the surface stain then reads absent while the protein is present.
  For furin and autolysis substrates the pieces stay associated, so this does
  not apply.
* A class-I PDZ-binding motif at the C-terminus is a sorting signal. Appending
  anything to it reproduces a known loss-of-function mutant. The motif is only
  2.85x enriched over the chance rate for these proteins' own amino-acid
  composition -- 420 hits against 147 expected -- so about a third are
  coincidence. It is graded "at risk", not forbidden, and a gene should be
  checked before the motif alone blocks a tag position.
* A protein whose C-terminus is cytoplasmic cannot display a C-terminal tag at
  all, which is the original panel finding: TFRC was the only member with an
  extracellular C-terminus.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from collections import Counter
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import candidate_inputs as ci

SHEDDASE = ("adam", "secretase", "mmp", "bace", "sheddase", "metalloproteinase")
ASSOCIATES = ("furin", "autolysis", "proprotein convertase", "subtilisin")


def pdz_class_i(seq: str) -> bool:
    """X-S/T-X-phi at positions -3..0; position -2 is third from the end."""
    return len(seq) >= 4 and seq[-3] in "ST" and seq[-1] in "AVILMFWCY"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default="library_tag_constraints.tsv")
    args = ap.parse_args()

    ref = {r["symbol"]: r for r in csv.DictReader(open("reference_selection.tsv"), delimiter="\t")}
    proc = {r["gene"]: r for r in csv.DictReader(open("library_processed_forms.tsv"), delimiter="\t")}

    seqs: dict[str, str] = {}
    for p in ci.uniprot_tables("human"):
        with gzip.open(p, "rt") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                g = (row.get("Gene Names (primary)") or "").strip()
                if g in ref and g not in seqs and row.get("Sequence"):
                    seqs[g] = row["Sequence"].strip()

    rows = []
    for g, r in sorted(ref.items()):
        reasons_no_cterm, reasons_no_nterm, at_risk_cterm = [], [], []
        if r["selected_gene_is_gpi"] == "True":
            reasons_no_cterm.append("GPI: the C-terminus is removed at the omega site")
        label = r["selected_topology_label"]
        if label in ("TM", "SP+TM") and g in seqs and pdz_class_i(seqs[g]):
            at_risk_cterm.append(f"class-I PDZ motif ({seqs[g][-4:]}), 2.85x over chance")
        p = proc.get(g)
        if p:
            by = (p["cleaved_by"] or "").lower()
            if any(x in by for x in SHEDDASE):
                reasons_no_nterm.append(f"sheddase substrate ({p['cleaved_by'][:40]})")
            elif not any(x in by for x in ASSOCIATES) and by.strip() in ("", "unspecified"):
                reasons_no_nterm.append("processed, protease unnamed - shedding unresolved")
        rows.append({
            "gene": g,
            "topology": label,
            "gpi": r["selected_gene_is_gpi"],
            "processed": "yes" if p else "no",
            "c_terminal_tag": ("forbidden" if reasons_no_cterm
                               else "at risk" if at_risk_cterm else "allowed"),
            "n_terminal_tag": ("at risk" if reasons_no_nterm else "allowed"),
            "c_terminal_reason": "; ".join(reasons_no_cterm + at_risk_cterm),
            "n_terminal_reason": "; ".join(reasons_no_nterm),
        })

    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    print(f"  {len(rows)} library members\n")
    for k, n in Counter(r["c_terminal_tag"] for r in rows).most_common():
        print(f"    C-terminal tag {k:10} {n:5}")
    for k, n in Counter(r["n_terminal_tag"] for r in rows).most_common():
        print(f"    N-terminal tag {k:10} {n:5}")
    hard = sum(1 for r in rows if r["c_terminal_tag"] == "forbidden" and r["n_terminal_tag"] == "at risk")
    free = sum(1 for r in rows if r["c_terminal_tag"] == "allowed" and r["n_terminal_tag"] == "allowed")
    print(f"\n    hard-blocked C-term and at-risk N-term : {hard:5}")
    print(f"    neither constraint                     : {free:5}  ({free/len(rows):.0%})")
    kind = lambda r: ("GPI omega site" if "GPI" in r["c_terminal_reason"]
                      else "PDZ motif" if "PDZ" in r["c_terminal_reason"] else "")
    print(f"\n  C-terminal blockers: {dict(Counter(kind(r) for r in rows if r['c_terminal_reason']))}")
    print(f"  N-terminal blockers: "
          f"{dict(Counter('sheddase, named' if 'sheddase substrate' in r['n_terminal_reason'] else 'protease unnamed' for r in rows if r['n_terminal_reason']))}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
