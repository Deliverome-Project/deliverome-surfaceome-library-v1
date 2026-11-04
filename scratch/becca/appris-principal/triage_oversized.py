"""Re-triage the oversized set against the authoritative gene list.

The H/M/L judgment in s1e3's oversized_value_annotation.md was written against
a 107-gene list that no longer matches: the set is now 99 genes, and the L-tier
prose is illustrative rather than exhaustive, so a third of the set is not named
anywhere in it. NOTCH4 went untriaged while NOTCH1, NOTCH2 and NOTCH3 were all
graded H, and PTPRZ2 was filed as "structural adhesion / matrix" despite being
a receptor phosphatase and PTPRZ1's paralog.

This regrades all 99 from recorded criteria rather than recollection, so the
reasons travel with the call.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path

# Families where this session has already built a working truncation, so a new
# member costs a build rather than a literature search.
RECIPE = {
    r"^NOTCH\d$": ("NOTCH", "NRR + TM + ICD"),
    r"^PTPR[BQSDFZ]\d?$": ("receptor phosphatase", "Ig block + juxtamembrane stub + TM + PTP domains"),
    r"^LRP(1B|2|4)$": ("LDLR family", "one ligand-binding cluster + TM + tail"),
    r"^FAT\d$": ("FAT cadherin", "EC1-4 + membrane-proximal cassette + TM + tail"),
    r"^DCHS\d$": ("cadherin", "EC1-4 + membrane-proximal cassette + TM + tail"),
    r"^CELSR\d$|^ADGR[A-Z]\d?$|^GPR179$": ("adhesion GPCR", "GPS-cleaved CTF + exogenous leader"),
    r"^MUC\d+[AB]?$": ("membrane mucin", "repeat stub + SEA + TM + tail, exogenous leader"),
    r"^PLXN[AB]\d$|^PLXND1$": ("plexin", "signal + 30 aa stub + TM + cytoplasmic"),
    r"^SDK\d$|^DSCAM": ("Ig adhesion", "N-terminal Ig block + FN3 stub + TM + tail"),
    r"^STAB\d$": ("stabilin", "C-terminal half, exogenous leader"),
}
# Classes the original triage set aside, with the reason restated per class.
SKIP = [
    (r"^CACNA1[A-S]$|^SCN\d+A$|^TRPM\d$|^PKD1L2$|^PKDREJ$", "voltage-gated or TRP channel",
     "multi-pass with little contiguous extracellular surface; the displayable face is short loops"),
    (r"^LCT$|^SI$|^MGAM2?$|^PCSK5$", "brush-border or secreted enzyme",
     "an enzyme rather than a targetable receptor"),
]


# The original triage recorded exactly one piece of gene-specific evidence for a
# skip -- "MUC3B (TPM 0)". A recipe existing does not make a gene worth building
# if it is not expressed, so expression gates the recommendation.
MIN_TPM = 1.0


def expression(gene: str, expr: dict) -> float | None:
    rec = expr.get(gene)
    if not rec:
        return None
    vals = [rec.get(k) for k in ("max_tpm_gtex", "max_ntpm_hpa_tissue", "max_tpm_hpa_celline")]
    vals = [v for v in vals if isinstance(v, (int, float))]
    return max(vals) if vals else None


def classify(gene: str) -> tuple[str, str, str]:
    for pat, (fam, recipe) in RECIPE.items():
        if re.match(pat, gene):
            return "inherits-recipe", fam, recipe
    for pat, cls, why in SKIP:
        if re.match(pat, gene):
            return "skip", cls, why
    return "untriaged", "", ""


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default="oversized_triage.tsv")
    args = ap.parse_args()

    ov = json.loads(Path("oversized_authoritative.json").read_text())
    area = {r["gene"]: r for r in csv.DictReader(open("oversized_extracellular_area.tsv"), delimiter="\t")}
    expr = json.loads(Path("oversized_expression.json").read_text())
    done = {r["gene"].split("_")[0] for r in csv.DictReader(open("giant_truncations.tsv"), delimiter="\t")}
    assessed = {r["gene"].split("_")[0]: r for r in
                csv.DictReader(open("giant_assessment_round2.tsv"), delimiter="\t")}

    rows = []
    for g, v in sorted(ov.items()):
        a = area.get(g, {})
        kind, fam, recipe = classify(g)
        if g in done:
            status, why = "built", "construct built this session"
        elif g in assessed:
            status, why = assessed[g]["verdict"], "assessed this session"
        elif a.get("displayable") == "no":
            status, why = "no-displayable-face", (
                f"{a.get('extracellular_aa')} aa extracellular, largest segment "
                f"{a.get('largest_segment_aa')} aa - nothing to display")
        elif kind == "inherits-recipe":
            tpm = expression(g, expr)
            if tpm is None:
                # Every expression source returns null for these. Absent is not
                # the same as zero, and must not pass silently as expressed --
                # s1e3 read MUC3B's absence as "TPM 0" and skipped it on that.
                status, why = "expression-unknown", (
                    f"a recipe exists ({fam}) but every expression source is null for "
                    f"this gene - confirm it is expressed before spending a build slot")
            elif tpm < MIN_TPM:
                status, why = "skip", (
                    f"a recipe exists ({fam}) but peak expression is {tpm:g} TPM - "
                    f"not expressed, so not worth a build slot")
            else:
                status, why = "ready-to-build", f"{fam}: {recipe}"
        elif kind == "skip":
            status, why = "skip", why
        else:
            status, why = "needs-triage", "no family recipe and no published assessment"
        rows.append({
            "gene": g, "orf_kb": v["orf_kb"], "tier": v["tier"],
            "extracellular_aa": a.get("extracellular_aa", ""),
            "largest_segment_aa": a.get("largest_segment_aa", ""),
            "tm_count": a.get("tm_count", ""),
            "peak_tpm": (lambda t: f"{t:g}" if t is not None else "")(expression(g, expr)),
            "status": status, "family": fam, "reason": why,
        })

    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    from collections import Counter
    c = Counter(r["status"] for r in rows)
    print(f"  {len(rows)} oversized genes re-triaged\n")
    for k, n in c.most_common():
        t1 = sum(1 for r in rows if r["status"] == k and r["tier"] == "T1")
        print(f"    {k:22} {n:4}   ({t1} of them T1)")
    print(f"\n  ready to build now, highest tier first:")
    rb = [r for r in rows if r["status"] == "ready-to-build"]
    for r in sorted(rb, key=lambda r: (r["tier"], -(r["orf_kb"] or 0))):
        print(f"    {r['gene']:9} {r['orf_kb']:5} kb  {r['tier']}  {r['family']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
