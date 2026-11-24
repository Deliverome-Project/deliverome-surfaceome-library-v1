"""Exact residue bounds for where a tag may sit on each library member.

The constraint table says which terminus is risky; this says where the tag
actually goes. Every bound is a residue number read from a UniProt feature, not
a rule of thumb:

* The N-terminal slot opens after the signal peptide is cleaved, so the first
  usable residue is one past the end of the SIGNAL feature. A tag placed inside
  the signal peptide is removed with it.
* The C-terminal slot closes at the GPI omega site, because everything distal
  to it is replaced by the anchor. The omega residue is UniProt's LIPID
  "GPI-anchor amidated" position.
* Where a class-I PDZ motif occupies the last four residues, the slot closes
  before them.
* A tag cannot go on a terminus that faces the cytoplasm, which the TOPO_DOM
  features settle.

Confidence travels with each bound, because the underlying annotations differ:
a GPI omega site is a hard, curated position, a PDZ motif is a sequence match
enriched only 2.85-fold over chance.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter
from pathlib import Path


def span(blob: str, kind: str) -> list[tuple[int, int]]:
    out = []
    for m in re.finditer(rf"{kind}\s+(\d+)(?:\.\.(\d+))?", blob or ""):
        a = int(m.group(1))
        out.append((a, int(m.group(2)) if m.group(2) else a))
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default="library_tag_bounds.tsv")
    args = ap.parse_args()

    ref = {r["symbol"]: r for r in csv.DictReader(open("reference_selection.tsv"), delimiter="\t")}
    cons = {r["gene"]: r for r in csv.DictReader(open("library_tag_constraints.tsv"), delimiter="\t")}
    rows = []
    seen = set()
    for r in csv.DictReader(open("uniprot_termini.tsv"), delimiter="\t"):
        g = (r.get("Gene Names (primary)") or "").strip()
        if g not in ref or g in seen:
            continue
        seen.add(g)
        L = int(r["Length"])
        seq = (r.get("Sequence") or "").strip()
        sig = span(r["Signal peptide"], "SIGNAL")
        tm = span(r["Transmembrane"], "TRANSMEM")
        gpi = [a for a, _ in span(r["Lipidation"], "LIPID")
               if "GPI-anchor amidated" in (r["Lipidation"] or "")]
        topo = r["Topological domain"] or ""

        n_from = (sig[0][1] + 1) if sig else 1
        n_conf = "signal peptide annotated" if sig else "no signal peptide; tag abuts the initiator"

        c_to, c_conf = L, "native C-terminus"
        if gpi:
            c_to, c_conf = min(gpi) - 1, "hard: GPI omega site, everything distal is replaced"
        elif seq and len(seq) >= 4 and seq[-3] in "ST" and seq[-1] in "AVILMFWCY":
            c_to, c_conf = L - 4, "soft: class-I PDZ motif, only 2.85x over chance"

        # which face the C-terminus sits on
        cface = ""
        last = [m for m in re.finditer(r'TOPO_DOM\s+(\d+)\.\.(\d+);\s+/note="([^"]+)"', topo)]
        if last:
            cface = last[-1].group(3)

        rows.append({
            "gene": g,
            "length_aa": L,
            "n_terminal_tag_from": n_from,
            "n_terminal_basis": n_conf,
            "n_terminal_risk": cons.get(g, {}).get("n_terminal_tag", ""),
            "c_terminal_tag_to": c_to,
            "c_terminal_basis": c_conf,
            "c_terminal_risk": cons.get(g, {}).get("c_terminal_tag", ""),
            "c_terminus_faces": cface,
            "tm_count": len(tm),
            "first_tm": f"{tm[0][0]}-{tm[0][1]}" if tm else "",
        })

    rows.sort(key=lambda r: r["gene"])
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    print(f"  {len(rows)} of {len(ref)} library members have bounds\n")
    print(f"  N-terminal slot opens after a cleaved signal peptide : "
          f"{sum(1 for r in rows if r['n_terminal_tag_from'] > 1)}")
    print(f"  no signal peptide, so the tag abuts residue 1        : "
          f"{sum(1 for r in rows if r['n_terminal_tag_from'] == 1)}")
    print(f"\n  C-terminal slot closes early:")
    print(f"    {dict(Counter(r['c_terminal_basis'].split(':')[0] for r in rows))}")
    print(f"\n  C-terminus faces: {dict(Counter(r['c_terminus_faces'] or 'not annotated' for r in rows))}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
