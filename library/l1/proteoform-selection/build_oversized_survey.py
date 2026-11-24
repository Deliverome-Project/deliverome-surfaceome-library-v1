"""One row per oversized gene, joining every record this work produced.

The detail lives in six files written at different points and for different
questions. This joins them so a future reader can see, for any of the 99 genes
that clear every library gate except the ORF cap, what was decided, on what
evidence, and what is still open -- without reconstructing the path.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def load(path: str, key: str = "gene") -> dict[str, dict]:
    p = Path(path)
    if not p.exists():
        return {}
    out = {}
    for r in csv.DictReader(p.open(), delimiter="\t"):
        out.setdefault(r[key].split("_")[0], r)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default="OVERSIZED_SURVEY.tsv")
    args = ap.parse_args()

    ov = json.loads(Path("oversized_authoritative.json").read_text())
    triage = load("oversized_triage.tsv")
    area = load("oversized_extracellular_area.tsv")
    built = load("giant_truncations.tsv")
    assess = load("giant_assessment_round2.tsv")
    cleaved = load("giant_cleaved_ranked.tsv")
    clinical = load("giant_clinical_value.tsv")
    spec = load("giant_truncation_spec.tsv")

    rows = []
    for g, v in sorted(ov.items()):
        t, a = triage.get(g, {}), area.get(g, {})
        b, s = built.get(g, {}), assess.get(g, {})
        c, cl = cleaved.get(g, {}), clinical.get(g, {})
        sp = spec.get(g, {})
        rows.append({
            "gene": g,
            "orf_kb": v["orf_kb"],
            "tier": v["tier"],
            "peak_tpm": t.get("peak_tpm", ""),
            "tm_count": a.get("tm_count", ""),
            "extracellular_aa": a.get("extracellular_aa", ""),
            "status": t.get("status", ""),
            "strategy": sp.get("strategy", "") or c.get("verdict", ""),
            "built_aa": b.get("truncated_aa", ""),
            "built_ltr_bp": b.get("ltr_to_ltr_tagged_bp", ""),
            "cleaved_fragment_aa": c.get("fragment_aa", ""),
            "cleaved_ectodomain_aa": c.get("ectodomain_kept", ""),
            "cleavage_evidence": c.get("evidence", ""),
            "assessment": s.get("verdict", ""),
            "confidence": s.get("confidence", "") or sp.get("confidence", ""),
            "pmid": s.get("pmid", "") or sp.get("pmid", ""),
            "clinical_flag": cl.get("clinical_flag", ""),
            "risk": (s.get("risk", "") or "")[:160],
            "detail_in": ";".join(x for x, present in [
                ("giant_truncations.tsv", b), ("giant_assessment_round2.tsv", s),
                ("giant_cleaved_ranked.tsv", c), ("giant_clinical_value.tsv", cl),
                ("oversized_triage.tsv", t)] if present),
        })

    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    from collections import Counter
    print(f"  {len(rows)} oversized genes, one row each\n")
    print(f"  status: {dict(Counter(r['status'] for r in rows))}")
    print(f"\n  how much is recorded per gene:")
    depth = Counter(len(r["detail_in"].split(";")) for r in rows)
    for k in sorted(depth, reverse=True):
        print(f"    {k} source file(s): {depth[k]:3} genes")
    have_design = sum(1 for r in rows if r["built_aa"] or r["cleaved_fragment_aa"])
    print(f"\n  genes with a concrete design (built or cleavage-derived): {have_design}")
    print(f"  genes with only a triage status and topology           : {len(rows) - have_design}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
