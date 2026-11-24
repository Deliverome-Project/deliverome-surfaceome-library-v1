"""Rank forms that no experiment can reach, using APPRIS's per-transcript annotations.

Observed evidence runs out fast. Of 13,687 alternate forms, 5,099 cannot be evidenced by
bottom-up MS even in principle and a further 7,308 simply have none, so roughly 12,400 forms
sit in a bucket no amount of PeptideAtlas or OpenProt will sort. APPRIS annotates **every**
transcript rather than only the principal, which gives a prediction-based ordering for
exactly that bucket.

The columns that matter here, from ``appris_data.appris.txt``:

* **TRIFID score** -- a machine-learning estimate of isoform functionality, and the quantity
  APPRIS itself uses to break ties between principal candidates.
* **APPRIS annotation** -- PRINCIPAL:n, ALTERNATIVE:1 (model conserved in at least three
  non-primate species), ALTERNATIVE:2 (fewer), or **MINOR**, meaning the core modules
  discarded it before tie-breaking began.
* **Transcript support level** -- Ensembl's own confidence in the transcript model, which is
  a different question from whether the protein is functional.
* **CRASH** (signal peptide) and **THUMP** (transmembrane helices) -- carried as columns
  only. **THUMP must not be used to decide whether a form is membrane-anchored**: it scores
  ADRB2, a seven-transmembrane GPCR, as **0** helices, and calls 3 for EGFR, 5 for CD44 and
  5 for TFRC where each has one. It is an APPRIS-internal scoring feature, not a topology
  predictor. The membrane test below uses **DeepTMHMM**, this project's topology of record,
  which returns 7, 1, 1 and 1 for those same proteins.

Prediction is not evidence, so the two are kept in separate columns and the agreement
between them is reported rather than blended into one number.

    uv run --locked python isoform_appris.py
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
import statistics
from collections import Counter, defaultdict
from pathlib import Path

DEEPTMHMM = "data/processed/deeptmhmm"

LABEL_RANK = {"PRINCIPAL": 0, "ALTERNATIVE:1": 1, "ALTERNATIVE:2": 2, "MINOR": 3}


def label_class(label: str) -> str:
    if label.startswith("PRINCIPAL"):
        return "PRINCIPAL"
    return label if label in LABEL_RANK else "MINOR"


def load_appris_full(path: Path) -> dict[str, dict]:
    """{ENST: per-transcript APPRIS annotations}."""
    out = {}
    with path.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            enst = row["Transcript ID"].split(".")[0]
            def num(key):
                v = row.get(key, "").strip()
                try:
                    return float(v)
                except (TypeError, ValueError):
                    return None
            out[enst] = {
                "label": row["APPRIS Annotation"].strip(),
                "trifid": num("Trifid Score"),
                "appris_score": num("APPRIS score"),
                "tsl": row["Transcript support level"].strip(),
                "ccds": "" if row["CCDS ID"].strip() == "-" else row["CCDS ID"].strip(),
                "corsair": num("Conservation (CORSAIR)"),
                "spade": num("Domain Score (SPADE)"),
                "thump": num("Trans-membrane helices (THUMP)"),
                "crash": row["Signal sequence (CRASH)"].strip(),
                "peptides": row["Peptides"].strip(),
            }
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cache-dir", default="~/.cache/deliverome/appris")
    ap.add_argument("--universe", default="isoform_universe_openprot.tsv")
    ap.add_argument("--surfaceome-repo", default="~/Git/accessible-surfaceome")
    args = ap.parse_args()

    cache = Path(args.cache_dir).expanduser()
    appris = load_appris_full(cache / "appris_full.txt")
    print(f"APPRIS transcript annotations: {len(appris)}")

    # DeepTMHMM TM-helix counts, keyed on the full accession (isoform suffix included).
    repo = Path(args.surfaceome_repo).expanduser()
    tm_count: dict[str, int] = {}
    for name in ("deeptmhmm_human_canonical.tsv", "deeptmhmm_human_isoforms.tsv"):
        path = repo / DEEPTMHMM / name
        if not path.exists():
            continue
        with path.open() as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                acc = row.get("uniprot_accession_full") or row.get("uniprot_accession")
                if acc:
                    tm_count[acc] = int(row["tm_helix_count"] or 0)
    print(f"DeepTMHMM TM counts: {len(tm_count)} accessions")

    rows = list(csv.DictReader(open(args.universe), delimiter="\t"))
    # Per gene, the reference's own APPRIS row, so an isoform can be compared to it.
    ref_ann: dict[str, dict] = {}
    ref_tm: dict[str, int] = {}
    for r in rows:
        if r["is_reference"] != "True":
            continue
        if r["uniprot"] in tm_count:
            ref_tm[r["gene"]] = tm_count[r["uniprot"]]
        for enst in (r["ensembl_transcripts"] or "").split(";"):
            a = appris.get(enst.strip())
            if a:
                ref_ann[r["gene"]] = a
                break

    matched = Counter()
    for r in rows:
        best = None
        for enst in (r["ensembl_transcripts"] or "").split(";"):
            a = appris.get(enst.strip())
            if a and (best is None or LABEL_RANK.get(label_class(a["label"]), 3)
                      < LABEL_RANK.get(label_class(best["label"]), 3)):
                best = a
        matched["matched" if best else "no APPRIS row"] += 1
        ref = ref_ann.get(r["gene"])
        r |= {
            "appris_label": best["label"] if best else "",
            "appris_label_class": label_class(best["label"]) if best else "",
            "trifid_score": best["trifid"] if best and best["trifid"] is not None else "",
            "transcript_support_level": best["tsl"] if best else "",
            "appris_ccds": best["ccds"] if best else "",
            "corsair_conservation": best["corsair"] if best and best["corsair"] is not None else "",
            "spade_domains": best["spade"] if best and best["spade"] is not None else "",
            "thump_tm_helices": best["thump"] if best and best["thump"] is not None else "",
            "crash_signal_peptide": best["crash"] if best else "",
            # Membrane anchoring, from DeepTMHMM rather than THUMP. Only a form that loses
            # *every* helix is excluded: losing some changes the protein but can still leave
            # a displayable construct, and that is a judgement for the designer.
            "deeptmhmm_tm_helices": tm_count.get(r["uniprot"], ""),
            "reference_tm_helices": ref_tm.get(r["gene"], ""),
            "loses_all_tm_vs_reference": (
                "True"
                if r["uniprot"] in tm_count
                and r["gene"] in ref_tm
                and ref_tm[r["gene"]] > 0
                and tm_count[r["uniprot"]] == 0
                else "False"
                if r["uniprot"] in tm_count and r["gene"] in ref_tm
                else ""
            ),
            "loses_some_tm_vs_reference": (
                "True"
                if r["uniprot"] in tm_count
                and r["gene"] in ref_tm
                and 0 < tm_count[r["uniprot"]] < ref_tm[r["gene"]]
                else "False"
                if r["uniprot"] in tm_count and r["gene"] in ref_tm
                else ""
            ),
            "trifid_vs_reference": (
                round(best["trifid"] - ref["trifid"], 4)
                if best and ref and best["trifid"] is not None and ref["trifid"] is not None
                else ""
            ),
        }

    out = Path("isoform_universe_scored.tsv")
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {out} ({len(rows)} forms)")

    alt = [r for r in rows if r["is_reference"] == "False"]
    with_tr = [r for r in alt if r["trifid_score"] != ""]

    def has_evidence(r) -> bool:
        return r["tier"] in ("expressed", "likely_expressed", "weak") or (
            r["openprot_te_score"] not in ("", None) and int(r["openprot_te_score"]) > 0
        )

    ev = [r for r in with_tr if has_evidence(r)]
    no_ev = [r for r in with_tr if not has_evidence(r)]
    unreachable = [
        r for r in with_tr
        if r["tier"] in ("undetectable", "not_detected", "not_detected_low_gene_coverage")
        and (r["openprot_te_score"] in ("", None) or int(r["openprot_te_score"]) == 0)
    ]

    summary = {
        "forms": len(rows),
        "alternate_forms": len(alt),
        "appris_match": dict(matched.most_common()),
        "with_trifid": len(with_tr),
        "appris_label_class": dict(
            Counter(r["appris_label_class"] for r in alt if r["appris_label_class"]).most_common()
        ),
        "trifid_median": {
            "forms_with_observed_evidence": round(
                statistics.median(float(r["trifid_score"]) for r in ev), 4
            ) if ev else None,
            "forms_without": round(
                statistics.median(float(r["trifid_score"]) for r in no_ev), 4
            ) if no_ev else None,
            "n_with": len(ev),
            "n_without": len(no_ev),
        },
        "unreachable_by_experiment": {
            "n": len(unreachable),
            "rankable_by_trifid": len(unreachable),
            "appris_label_class": dict(
                Counter(r["appris_label_class"] for r in unreachable).most_common()
            ),
            "trifid_over_0_5": sum(1 for r in unreachable if float(r["trifid_score"]) > 0.5),
        },
        "membrane_test_deeptmhmm": {
            "forms_with_a_call": sum(1 for r in alt if r["loses_all_tm_vs_reference"] != ""),
            "loses_ALL_tm": sum(1 for r in alt if r["loses_all_tm_vs_reference"] == "True"),
            "loses_SOME_tm_kept": sum(
                1 for r in alt if r["loses_some_tm_vs_reference"] == "True"
            ),
            "loses_all_and_in_cohort": sum(
                1
                for r in alt
                if r["loses_all_tm_vs_reference"] == "True" and r["in_current_cohort"] == "True"
            ),
        },
    }
    Path("isoform_appris_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
