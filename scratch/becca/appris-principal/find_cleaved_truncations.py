"""Giants that shrink to within budget on their own proteolytic cleavage.

The best-evidenced way to truncate a giant is not to design a cut but to use the
one the cell already makes. Polycystin-1 is the clearest case -- the GPS-cleaved
C-terminal fragment is the only species at the plasma membrane -- and the same
holds for the receptor phosphatase P-subunits, the adhesion-GPCR CTFs and the
retained mucin fragments.

A candidate has to clear three tests, each of which rejected real cases while
this was written:

* The cut must be genuinely internal. A CHAIN feature starting at residue 20 is
  the mature protein after signal-peptide removal, not a cut: LRP1 and SI were
  both reported as shortened on exactly that.
* The retained fragment must still contain the transmembrane span. NOTCH1's
  smallest annotated chain is the intracellular domain, which is soluble.
* The fragment must actually be shorter. A chain covering 98% of the protein
  saves nothing.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path

FIXED, OVERHEAD, CASSETTE, BARCODE = 4087, 131, 153, 20
CEILING = 9000
MAX_FRACTION = 0.9          # the fragment must be a real shortening


def genome(aa: int) -> int:
    return FIXED + aa * 3 + 3 + OVERHEAD + CASSETTE + BARCODE


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default="giant_cleaved_candidates.tsv")
    args = ap.parse_args()

    ov = json.loads(Path("oversized_authoritative.json").read_text())
    proc, term = {}, {}
    for r in csv.DictReader(open("uniprot_processing.tsv"), delimiter="\t"):
        g = (r.get("Gene Names (primary)") or "").strip()
        if g in ov and g not in proc:
            proc[g] = r
    for r in csv.DictReader(open("uniprot_termini.tsv"), delimiter="\t"):
        g = (r.get("Gene Names (primary)") or "").strip()
        if g in ov and g not in term:
            term[g] = r

    site_col = next(c for c in next(iter(proc.values())) if "Site" in c)
    chain_col = next(c for c in next(iter(proc.values())) if "Chain" in c)

    rows, rejected = [], []
    for g, r in sorted(proc.items()):
        t = term.get(g)
        if not t:
            continue
        L = int(t["Length"])
        tms = [(int(a), int(b)) for a, b in
               re.findall(r"TRANSMEM\s+(\d+)\.\.(\d+)", t["Transmembrane"] or "")]
        if not tms:
            continue
        tm_start, tm_end = tms[0]
        sig = re.search(r"SIGNAL\s+\d+\.\.(\d+)", t["Signal peptide"] or "")
        mature = int(sig.group(1)) + 1 if sig else 1

        cands = []
        for m in re.finditer(r'SITE\s+(\d+)(?:\.\.\d+)?;\s+/note="([^"]*)"', r[site_col] or ""):
            if "cleav" in m.group(2).lower():
                cands.append((int(m.group(1)), m.group(2)[:40], "cleavage site"))
        for m in re.finditer(r'CHAIN\s+(\d+)\.\.(\d+);\s+/note="([^"]*)"', r[chain_col] or ""):
            a = int(m.group(1))
            if a > mature:
                cands.append((a - 1, m.group(3)[:40], "chain boundary"))

        for pos, note, kind in cands:
            frag = L - pos
            if pos <= mature:
                rejected.append((g, pos, "cut is the mature-protein boundary, not internal"))
                continue
            if pos >= tm_start:
                rejected.append((g, pos, "fragment would not contain the transmembrane span"))
                continue
            if frag > L * MAX_FRACTION:
                rejected.append((g, pos, "fragment is not meaningfully shorter"))
                continue
            rows.append({
                "gene": g, "full_aa": L, "cut_after": pos, "fragment_aa": frag,
                "shrink": f"{frag/L:.0%}", "ltr_to_ltr_bp": genome(frag),
                "within_9kb": genome(frag) <= CEILING,
                "evidence": kind, "note": note,
                "needs_exogenous_leader": "yes",
            })

    best: dict[str, dict] = {}
    for r in rows:
        if r["gene"] not in best or r["fragment_aa"] < best[r["gene"]]["fragment_aa"]:
            best[r["gene"]] = r
    out = sorted(best.values(), key=lambda r: r["fragment_aa"])

    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out)

    fit = sum(1 for r in out if r["within_9kb"])
    print(f"  {len(out)} oversized genes have a usable internal cleavage; {fit} fit 9 kb\n")
    print(f"  {'gene':10} {'full':>5} {'cut':>6} {'frag':>6} {'shrink':>7} {'genome':>7}  evidence")
    for r in out:
        print(f"  {r['gene']:10} {r['full_aa']:>5} {r['cut_after']:>6} {r['fragment_aa']:>6} "
              f"{r['shrink']:>7} {r['ltr_to_ltr_bp']:>7}  {r['evidence']}")
    from collections import Counter
    print(f"\n  rejected candidates by reason: "
          f"{dict(Counter(x[2] for x in rejected))}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
