"""What the 5,000-nt ORF cap removes from L1, and the packaging budget it was set against.

The cap is ``MAX_ORF_NT`` in the S1E3 filter
(``library/l1/surface-filters/final_library.py``), applied to the gene's
UniProt canonical ORF (residues x 3 + stop). This reports, from pinned inputs:

* every cohort gene over the cap, and which of them the cap alone keeps out: genes that pass
  every other core gate, and genes that meet every other low-literature tier condition;
* screen-validated GPCRs over the cap (that tier does not apply it);
* the largest designed L1 ORFs, since the cap is tested on the canonical, not the ordered form;
* which capped genes have a truncation construct in L2;
* the vector's fixed LTR-to-LTR load, measured from the pinned p0039 backbone map.

Writes ``orf_size_cap_genes.tsv`` (one row per cohort gene over the cap) and
``orf_size_cap_summary.json``. Runs the S1E3 filter offline (``library_pipeline``), or reuses
a staged run with ``--s1e3``.

    uv run --locked python library/l1/orf-design/orf_size_cap.py
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import tempfile
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
AP = REPO / "library/l1/proteoform-selection"
OUT = HERE / "orf_size_cap_genes.tsv"
SUMMARY = HERE / "orf_size_cap_summary.json"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with Path(path).open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def staged_s1e3(arg: Path | None) -> Path:
    if arg:
        return arg.resolve()
    sys.path.insert(0, str(HERE))
    import library_pipeline

    work = Path(tempfile.mkdtemp(prefix="gene-filter-"))
    library_pipeline.run_gene_filter(work)
    return work / "s1e3"


def ltr_load() -> dict:
    """Fixed LTR-to-LTR length of the pinned backbone, with its dropout stuffer removed."""
    sys.path.insert(0, str(AP))
    import candidate_inputs as inputs
    from Bio import SeqIO

    rec = SeqIO.read(inputs.backbone(), "genbank")
    feats = {f.qualifiers.get("label", [""])[0]: f for f in rec.features}
    five = next(f for k, f in feats.items() if k.startswith("5' LTR"))
    three = next(f for k, f in feats.items() if k.startswith("3' LTR"))
    stuffer = next(f for k, f in feats.items() if k.startswith("stuffer"))
    span = int(three.location.end) - int(five.location.start)
    stuffer_bp = int(stuffer.location.end) - int(stuffer.location.start)
    return {
        "backbone": f"processed://{inputs.BACKBONE_KEY}",
        "ltr_to_ltr_bp_with_stuffer": span,
        "stuffer_bp": stuffer_bp,
        "ltr_to_ltr_bp_without_insert": span - stuffer_bp,
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--s1e3", type=Path, help="a staged S1E3 folder from run_gene_filter")
    a = ap.parse_args(argv)
    s1e3 = staged_s1e3(a.s1e3)
    sys.path.insert(0, str(s1e3))
    import final_library as F  # importing runs the filter

    membership = {r["gene"]: r["tier"] for r in read_tsv(s1e3 / "library_membership.tsv")}
    sv = json.loads((s1e3 / "willow_screen_rescue_genes.json").read_text())
    rows_uni = json.loads((s1e3 / "catalog_rows.json").read_text())
    uni = {r["symbol"]: bool(r["db"] & 1) for r in rows_uni}

    giants = read_tsv(AP / "giant_truncations.tsv")
    giant_genes = Counter(g["gene"].split("_")[0] for g in giants)

    rows = []
    for r in F.recs:
        g, f = r["gene"], r["filters"]
        n = F.orf_nt(r)
        if n is None or n <= F.MAX_ORF_NT:
            continue
        core_else = F._passes_except_orf(r)
        ll_else = (
            g not in sv
            and uni.get(g)
            and f.get("n_papers_found") is not None
            and f["n_papers_found"] < 100
            and not F.B.passes_canonical(f)
            and F.anchored(r)
            and (F.mx(g) or 0) >= F.EXPR_MIN_TPM
        )
        rows.append(
            {
                "gene": g,
                "uniprot": r.get("uniprot") or "",
                "canonical_aa": (n - 3) // 3,
                "canonical_orf_nt": n,
                "first_failing_gate": F.first_failing_gate(r) or "",
                "passes_every_other_core_gate": core_else,
                "meets_every_other_low_literature_condition": bool(ll_else) and not core_else,
                "screen_validated_gpcr": g in sv,
                "in_l1_gene_list": g in membership,
                "l2_truncation_constructs": giant_genes.get(g, 0),
            }
        )
    rows.sort(key=lambda d: -d["canonical_orf_nt"])
    with OUT.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    kept_out = [
        d
        for d in rows
        if d["passes_every_other_core_gate"] or d["meets_every_other_low_literature_condition"]
    ]
    manifest = read_tsv(HERE / "library_manifest_l1.tsv")
    orf_nt = sorted(
        ((int(m["ordered_length_aa"]) + 1) * 3, m["candidate_id"], m["gene_symbol"])
        for m in manifest
    )
    summary = {
        "max_orf_nt": F.MAX_ORF_NT,
        "applied_to": "UniProt canonical ORF of the gene record (residues x 3 + stop)",
        "cohort_genes_over_cap": len(rows),
        "first_failing_gate_is_orf_size": sum(
            1 for d in rows if d["first_failing_gate"] == "ORF size"
        ),
        "pass_every_other_core_gate": sum(1 for d in rows if d["passes_every_other_core_gate"]),
        "meet_every_other_low_literature_condition_only": sum(
            1 for d in rows if d["meets_every_other_low_literature_condition"]
        ),
        "kept_out_of_l1_by_the_cap_alone": len(kept_out),
        "screen_validated_gpcrs_over_cap": [d["gene"] for d in rows if d["screen_validated_gpcr"]],
        "over_cap_but_in_l1_gene_list": [d["gene"] for d in rows if d["in_l1_gene_list"]],
        "kept_out_with_l2_truncation_constructs": sorted(
            d["gene"] for d in kept_out if d["l2_truncation_constructs"]
        ),
        "l2_truncation_constructs": len(giants),
        "kept_out_canonical_aa_range": [
            min(d["canonical_aa"] for d in kept_out),
            max(d["canonical_aa"] for d in kept_out),
        ],
        "largest_designed_l1_orfs_nt": [f"{c} {g} {n}" for n, c, g in orf_nt[-5:]],
        "l1_designed_orfs_over_cap": sum(1 for n, _c, _g in orf_nt if n > F.MAX_ORF_NT),
        "vector": ltr_load(),
        "cap_comment_in_filter": (
            "Lentiviral packaging: LTR-LTR budget ~8.5 kb. Current p003 construct -> ORFs up to "
            "~4 kb; with mScarlet removed ~4.7 kb. Provisional cap 5 kb."
        ),
    }
    SUMMARY.write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps(summary, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
