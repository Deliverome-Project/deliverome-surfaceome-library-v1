"""Two questions about the references we build from: how well are they covered, and can
OpenProt adjudicate the MANE/APPRIS disagreements the peptide test could not?

**Coverage.** Reports what fraction of each selected reference's sequence is covered by
observed PeptideAtlas peptides, and what fraction has any OpenProt evidence at all. A
reference nobody has ever seen is not wrong, but it is worth knowing which ones they are.

**Adjudication.** The peptide test failed on 82 of the 111 MANE/APPRIS disagreements, and
mostly for a structural reason: where one candidate is a subsequence of the other there is
no peptide that can tell them apart. **OpenProt does not have that problem.** It scores each
form as its own row, keyed on its transcript, so the MANE form and the APPRIS form get
independent MS and TE scores even when one contains the other. TE in particular is
ribosome-profiling evidence, which reads translation of a transcript rather than presence of
a distinguishing peptide.

Its limits still apply: OpenProt 2.2 is Ensembl GRCh38.106 / UniProt 2022_06, so a
transcript newer than that has no row, and a TE of 0 can mean "not profiled in a tissue that
expresses it" as easily as "not translated".

    uv run --locked python reference_openprot.py
"""

from __future__ import annotations

import argparse
import csv
import json
import statistics
from collections import Counter
from pathlib import Path

from openprot_evidence import OPENPROT_RELEASE, enst_to_ensp, load_openprot

# How much better one candidate's score must be before it counts as a call rather than noise.
TE_MARGIN = 2
MS_MARGIN = 3


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--openprot-dir", default="~/.cache/deliverome/openprot")
    ap.add_argument("--cache-dir", default="~/.cache/deliverome/appris")
    ap.add_argument("--universe", default="isoform_universe_openprot.tsv")
    ap.add_argument("--disagreements", default="reference_disagreements.tsv")
    args = ap.parse_args()

    # ---- 1. coverage of the references we actually build from --------------------
    universe = list(csv.DictReader(open(args.universe), delimiter="\t"))
    refs = [r for r in universe if r["is_reference"] == "True"]
    cov = [float(r["gene_coverage_pct"]) for r in refs]
    te_scored = [r for r in refs if r["openprot_te_score"] not in ("", None)]

    def pct(n: int) -> str:
        return f"{100 * n / len(refs):.1f}%"

    coverage = {
        "references": len(refs),
        "median_sequence_coverage_pct": round(statistics.median(cov), 1),
        "with_any_peptide_coverage": f"{sum(1 for c in cov if c > 0)} ({pct(sum(1 for c in cov if c > 0))})",
        "coverage_over_10pct": f"{sum(1 for c in cov if c >= 10)} ({pct(sum(1 for c in cov if c >= 10))})",
        "coverage_over_50pct": f"{sum(1 for c in cov if c >= 50)} ({pct(sum(1 for c in cov if c >= 50))})",
        "zero_coverage": f"{sum(1 for c in cov if c == 0)} ({pct(sum(1 for c in cov if c == 0))})",
        "with_openprot_row": f"{len(te_scored)} ({pct(len(te_scored))})",
        "openprot_translated": sum(1 for r in te_scored if int(r["openprot_te_score"]) > 0),
    }
    print(json.dumps({"reference_coverage": coverage}, indent=2))

    Path("reference_coverage.tsv").write_text(
        "gene\tlength\tcoverage_pct\tgene_peptides\topenprot_ms\topenprot_te\n"
        + "".join(
            f"{r['gene']}\t{r['length']}\t{r['gene_coverage_pct']}\t{r['gene_peptides']}\t"
            f"{r['openprot_ms_score']}\t{r['openprot_te_score']}\n"
            for r in sorted(refs, key=lambda r: float(r["gene_coverage_pct"]))
        )
    )
    print(f"wrote reference_coverage.tsv ({len(refs)} references)")

    # ---- 2. can OpenProt call the MANE/APPRIS disagreements? ---------------------
    op = load_openprot(sorted(Path(args.openprot_dir).expanduser().glob("*min_2_pep*.tsv"))[0])
    e2p = enst_to_ensp(Path(args.cache_dir).expanduser() / "Homo_sapiens.GRCh38.pep.all.fa.gz")
    dis = list(csv.DictReader(open(args.disagreements), delimiter="\t"))

    def scores(ensts: str) -> tuple[int | None, int | None]:
        best_ms = best_te = None
        for enst in (ensts or "").split(";"):
            ensp = e2p.get(enst.strip())
            rec = op.get(ensp) if ensp else None
            if rec:
                best_ms = max(best_ms or 0, rec["ms"])
                best_te = max(best_te or 0, rec["te"])
        return best_ms, best_te

    out = []
    for d in dis:
        m_ms, m_te = scores(d["mane_enst"])
        a_ms, a_te = scores(d["appris_enst"])
        verdict = "no_openprot_data"
        if m_te is not None and a_te is not None:
            if m_te - a_te >= TE_MARGIN:
                verdict = "openprot_favours_mane"
            elif a_te - m_te >= TE_MARGIN:
                verdict = "openprot_favours_appris"
            elif m_ms is not None and a_ms is not None and m_ms - a_ms >= MS_MARGIN:
                verdict = "openprot_favours_mane_on_ms"
            elif m_ms is not None and a_ms is not None and a_ms - m_ms >= MS_MARGIN:
                verdict = "openprot_favours_appris_on_ms"
            else:
                verdict = "openprot_tied"
        elif m_te is not None:
            verdict = "only_mane_has_a_row"
        elif a_te is not None:
            verdict = "only_appris_has_a_row"
        out.append(
            dict(
                d,
                openprot_mane_ms=m_ms if m_ms is not None else "",
                openprot_mane_te=m_te if m_te is not None else "",
                openprot_appris_ms=a_ms if a_ms is not None else "",
                openprot_appris_te=a_te if a_te is not None else "",
                openprot_verdict=verdict,
            )
        )

    cols = list(out[0])
    with open("reference_disagreements_openprot.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(out)
    print(f"wrote reference_disagreements_openprot.tsv ({len(out)} genes)")

    unresolved = [r for r in out if r["decision"] == "unresolved_default_mane"]
    summary = {
        "openprot_release": OPENPROT_RELEASE,
        "disagreements": len(out),
        "openprot_verdict": dict(Counter(r["openprot_verdict"] for r in out).most_common()),
        "previously_unresolved": {
            "n": len(unresolved),
            "now_called_by_openprot": sum(
                1 for r in unresolved if r["openprot_verdict"].startswith("openprot_favours")
            ),
            "verdict": dict(
                Counter(r["openprot_verdict"] for r in unresolved).most_common()
            ),
        },
        "agreement_where_both_spoke": dict(
            Counter(
                f"peptides={r['decision']} / {r['openprot_verdict']}"
                for r in out
                if r["decision"].startswith(("resolved", "weak"))
                and r["openprot_verdict"].startswith("openprot_favours")
            ).most_common()
        ),
    }
    Path("reference_openprot_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
