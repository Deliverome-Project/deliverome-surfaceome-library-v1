"""How much of each oversized protein is actually outside the cell.

KIDINS220 reached our oversized shortlist with nine extracellular residues in
total, which no truncation can fix -- the size cap was never its problem, the
topology was. It entered on the surfaceome catalog's looser "likely" tier rather
than the canonical one, so the obvious question is how many others came in the
same way. This measures the displayable face before any design effort is spent.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from deliverome_analysis import surfaceome

MIN_LOOP = 8        # the tag-site screen needs a contiguous run of at least this
MIN_TOTAL = 50      # below this there is no meaningful antigen to display
# UniProt writes the outside face as "Extracellular" on plasma-membrane proteins
# and "Lumenal" on anything annotated through the secretory pathway. Matching only
# the first silently reported SORL1's 2,056-residue ectodomain as zero.
OUTSIDE = ("extracellular", "lumenal", "lumenal, melanosome", "perinuclear space")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--features", default="giant_features.json")
    ap.add_argument("--out", default="oversized_extracellular_area.tsv")
    args = ap.parse_args()

    feats = json.loads(Path(args.features).read_text())
    canonical = surfaceome.canonical_accessions()
    likely = surfaceome.surface_accessions()

    rows = []
    for gene, rec in sorted(feats.items()):
        topo = [f for f in rec["features"] if f["type"] == "Topological domain"]
        tm = sorted((f for f in rec["features"] if f["type"] == "Transmembrane"),
                    key=lambda f: f["start"])
        sig = next((f for f in rec["features"] if f["type"] == "Signal"), None)
        out = [f for f in topo
               if any(f["description"].lower().startswith(o) for o in OUTSIDE)]
        spans = [f["end"] - f["start"] + 1 for f in out]
        total = sum(spans)
        inferred = False
        if not spans and len(tm) == 1:
            # Several entries annotate only the cytoplasmic side. For a single-pass
            # protein whose C-terminal domain is cytoplasmic, the stretch between
            # the signal peptide and the TM is the ectodomain by construction.
            cyto = [f for f in topo if f["description"].lower().startswith("cytoplasmic")]
            if cyto and min(c["start"] for c in cyto) > tm[0]["end"]:
                start = (sig["end"] + 1) if sig else 1
                length = tm[0]["start"] - start
                if length > 0:
                    spans, total, inferred = [length], length, True
        acc = rec["accession"].split("-")[0]
        rows.append({
            "gene": gene,
            "accession": acc,
            "length_aa": rec["length"],
            "tm_count": len(tm),
            "extracellular_segments": len(out),
            "extracellular_aa": total if spans else "",
            "largest_segment_aa": max(spans) if spans else "",
            "pct_extracellular": f"{total / rec['length']:.1%}" if spans else "",
            "source": "inferred from TM" if inferred else ("annotated" if spans else "none"),
            "surfaceome_tier": ("canonical" if acc in canonical
                                else "likely" if acc in likely else "neither"),
            "displayable": "" if not spans else (
                "no" if (total < MIN_TOTAL or max(spans) < MIN_LOOP) else "yes"),
        })

    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    bad = [r for r in rows if r["displayable"] == "no"]
    unknown = [r for r in rows if r["displayable"] == ""]
    if unknown:
        print("  no topology to judge from, so not cleared either way: "
              + ", ".join(r["gene"] for r in unknown) + "\n")
    print(f"  {'gene':12} {'aa':>6} {'TMs':>4} {'extracellular':>14} {'largest':>8} {'tier':>10}  displayable")
    for r in sorted(rows, key=lambda r: (r["displayable"] != "no", r["gene"])):
        if r["displayable"] == "no" or r["surfaceome_tier"] != "canonical":
            print(f"  {r['gene']:12} {r['length_aa']:6} {r['tm_count']:4} "
                  f"{str(r['extracellular_aa']) + ' aa ' + r['pct_extracellular']:>14} "
                  f"{str(r['largest_segment_aa']):>8} {r['surfaceome_tier']:>10}  {r['displayable'] or '?'}")
    print(f"\n  {len(bad)} with no displayable face; {len(unknown)} with no topology annotation")
    print(f"  of the {len(bad)} undisplayable, surfaceome tier: "
          f"{ {t: sum(1 for r in bad if r['surfaceome_tier'] == t) for t in ('canonical','likely','neither')} }")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
