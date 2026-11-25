"""Per-gate gene lists for reviewers: who passed and who failed each L1 gene-list gate, and why.

The brief's section 2 gives one count per gate ("Excluded here") and per tier. A reviewer
checking a count, or looking for one gene, needs the genes behind it and the value the gate
actually tested. This writes one table per gate and one per entry tier into ``gate_lists/``,
plus an index (``gate_lists/README.md``) that names the brief number each file backs.

It reruns the cohort filter offline, exactly as the pipeline's ``gene_list`` stage does
(``library_pipeline.run_gene_filter``: accessible-surfaceome release 1.3.0, the pinned
expression tables, the replayed MyGene responses), then imports the filter module that run
produced and reads each gene's tested value from the same functions ``passes`` calls. It
asserts that every pass/fail it writes agrees with the committed
``gene_gate_outcomes.tsv`` and that the tier lists equal the filter's own
``library_membership.tsv``, so these lists cannot drift from the decision they document.

File names carry no experiment IDs: the review repository is read by people outside the lab.

    uv run --locked python gate_lists.py              # stages the filter inputs from S3
    uv run --locked python gate_lists.py --work DIR   # run into DIR, or reuse a finished run there
"""

from __future__ import annotations

import argparse
import csv
import shutil
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
S4E9 = HERE.parent
OUT = HERE / "gate_lists"
OUTCOMES = S4E9 / "gene_gate_outcomes.tsv"
EXCLUSIONS = S4E9 / "gene_list_exclusions.tsv"

#: gate name in final_library.GATES -> output file stem (ordered as the filter applies them)
GATE_FILES = {
    "surface evidence": "01_surface_evidence",
    "accessibility": "02_accessibility",
    "surface fraction": "03_surface_fraction",
    "anchored": "04_anchored",
    "co-receptor": "05_coreceptor",
    "expression": "06_expression",
    "localization": "07_localization",
    "ORF size": "08_orf_size",
    "ectodomain": "09_ectodomain",
}
TIER_FILES = {
    "core": ("10_tier_core", "passes every gate (tiers T1-T3 by architecture risk)"),
    "LL": ("11_tier_low_literature_uniprot", "low-literature UniProt surface tier"),
    "SV": ("12_tier_screen_validated_gpcr", "screen-validated GPCR tier (Tedman et al., 2026)"),
    "OE": ("13_tier_overexpression_surface", "overexpression-surface tier"),
}


def fmt(x) -> str:
    if x is None:
        return ""
    if isinstance(x, float):
        return f"{x:.2f}"
    if isinstance(x, (list, tuple)):
        return ";".join(fmt(v) for v in x)
    return str(x)


def tested_values(F, r) -> dict[str, dict[str, str]]:
    """Per gate, the fields the gate's predicate reads, as the filter module computes them."""
    g, f = r["gene"], r["filters"]
    topo = r.get("canonical_topology") or {}
    partners = (F.part.get(g) or {}).get("partners", [])
    non_self = [p for p in partners if not F._is_self_partner(g, p)]
    meds = {str(p): F.partner_median(p) for p in non_self}
    return {
        "surface evidence": {"deep_dive_tier": fmt(f.get("deep_dive_tier"))},
        "accessibility": {"surface_accessibility": fmt(f.get("surface_accessibility"))},
        "surface fraction": {"surface_specificity": fmt(f.get("surface_specificity"))},
        "anchored": {
            "tm_helix_count_deeptmhmm_v1": fmt(topo.get("tm_helix_count")),
            "subcategory": fmt(f.get("subcategory")),
        },
        "co-receptor": {
            "co_receptor_dependency": fmt(f.get("co_receptor_dependency")),
            "partners": " | ".join(str(p) for p in partners),
            "partner_median_celline_tpm": " | ".join(
                f"{p}={fmt(m) if m is not None else 'unresolved'}" for p, m in meds.items()
            ),
            "allowlisted_alone": "yes" if g in F.CORECEPTOR_ALLOWLIST else "",
        },
        "expression": {
            "low_endogenous_expression": fmt(f.get("low_endogenous_expression")),
            "max_tpm": fmt(F.mx(g)),
        },
        "localization": {"surface_call_reason": fmt(f.get("surface_call_reason"))},
        "ORF size": {"orf_nt": fmt(F.orf_nt(r))},
        "ectodomain": {
            "extracellular_residues_deeptmhmm_v1": fmt(F.ecd_residues(r)),
            "subcategory": fmt(f.get("subcategory")),
        },
    }


def write(path: Path, rows: list[dict], fieldnames: list[str] | None = None) -> None:
    cols = fieldnames if fieldnames is not None else list(rows[0])
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def tier_rows(F, s1: Path, membership: dict, removed: dict) -> dict[str, list[dict]]:
    """Entry tiers: core = every gate passed; the three named tiers as the filter assigned them."""
    by_gene = {r["gene"]: r for r in F.recs}
    # What admitted each named-tier gene: the screen's anti-HA surface immunostaining (a.u.)
    # for SV, the recorded overexpression call for OE, the UniProt surface badge for LL.
    sv = {g: v.get("immunostaining") for g, v in F._sv.items()}
    with (s1 / "oe_calls_136.tsv").open() as fh:
        oe = {r["gene"]: r["reason"] for r in csv.DictReader(fh, delimiter="\t")}
    evidence = {
        "LL": lambda g: "UniProt cell-surface location (catalog low_lit_uniprot badge)",
        "SV": lambda g: f"anti-HA surface immunostaining {fmt(sv.get(g))} a.u.",
        "OE": lambda g: oe.get(g, ""),
    }
    tiers: dict[str, list[dict]] = {k: [] for k in TIER_FILES}
    for g, t in sorted(membership.items()):
        key = t if t in ("LL", "SV", "OE") else "core"
        r = by_gene[g]
        tiers[key].append(
            {
                "gene": g,
                "uniprot_accession": r.get("uniprot") or "",
                "filter_tier": t,
                "first_gate_failed": F.first_failing_gate(r) or "",
                "n_papers_found": fmt(r["filters"].get("n_papers_found")),
                "max_tpm": fmt(F.mx(g)),
                "tier_evidence": evidence[key](g) if key in evidence else "",
                "removed_after_filter": removed.get(g, ""),
            }
        )
    return tiers


def library_gene_rows(committed: dict, membership: dict, removed: dict) -> tuple[list, set]:
    """One table of the library's genes and the set the gene filter contributed.

    The filter's genes less the recorded exclusions, plus the dated additions still active
    (gene_list_additions.tsv; a withdrawn row is not).
    """
    with (HERE.parent / "gene_list_additions.tsv").open() as fh:
        added = [
            r for r in csv.DictReader(fh, delimiter="\t") if not (r.get("withdrawn") or "").strip()
        ]
    lib_rows = [
        {"gene": g, "entered_by": committed[g]["entered_by"], "added": "", "source": "gene filter"}
        for g in sorted(membership)
        if g not in removed
    ]
    in_lib = {r["gene"] for r in lib_rows}
    lib_rows += [
        {
            "gene": r["gene"],
            "entered_by": "dated addition",
            "added": r["added"],
            "source": r["source"],
        }
        for r in added
        if r["gene"] not in in_lib
    ]
    return lib_rows, in_lib


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--work",
        type=Path,
        help="run_gene_filter work directory: reused when finished, else run into (must be empty)",
    )
    a = ap.parse_args(argv)

    work = a.work.resolve() if a.work else Path(tempfile.mkdtemp(prefix="gate-lists-"))
    if not (work / "s1e3" / "library_membership.tsv").exists():
        sys.path.insert(0, str(S4E9))
        import library_pipeline as LP

        if work.exists() and any(work.iterdir()):
            raise SystemExit(f"{work} holds a partial run; pass an empty or finished directory")
        work.mkdir(parents=True, exist_ok=True)
        LP.run_gene_filter(work)
    s1 = work / "s1e3"
    sys.path.insert(0, str(s1))
    import final_library as F  # importing runs the filter

    with OUTCOMES.open() as fh:
        committed = {r["gene"]: r for r in csv.DictReader(fh, delimiter="\t")}
    with (s1 / "library_membership.tsv").open() as fh:
        membership = {r["gene"]: r["tier"] for r in csv.DictReader(fh, delimiter="\t")}
    with EXCLUSIONS.open() as fh:
        removed = {r["gene"]: r["reason"] for r in csv.DictReader(fh, delimiter="\t")}

    if OUT.exists():
        shutil.rmtree(OUT)
    OUT.mkdir()

    gate_names = [name for name, _ in F.GATES]
    assert gate_names == list(GATE_FILES), f"filter gates changed: {gate_names}"
    per_gate: dict[str, list[dict]] = {n: [] for n in gate_names}
    for r in F.recs:
        g = r["gene"]
        first = F.first_failing_gate(r)
        vals = tested_values(F, r)
        for name, ok in F.GATES:
            result = "pass" if ok(r) else "fail"
            assert committed[g][name] == result, f"{g} {name}: {result} vs committed table"
            per_gate[name].append(
                {
                    "gene": g,
                    "uniprot_accession": r.get("uniprot") or "",
                    "result": result,
                    **vals[name],
                    "library_tier": membership.get(g, ""),
                    "removed_after_filter": "yes" if g in removed else "",
                }
            )
            if name == first:
                break  # later gates were never applied to this gene

    counts = {}
    for name, rows in per_gate.items():
        rows.sort(key=lambda d: (d["result"] != "fail", d["gene"]))
        cols = list(rows[0])
        passed = [d for d in rows if d["result"] == "pass"]
        excluded = [d for d in rows if d["result"] == "fail"]
        stem = GATE_FILES[name]
        write(OUT / f"{stem}.tsv", rows, cols)
        write(OUT / f"{stem}_passed.tsv", passed, cols)
        write(OUT / f"{stem}_excluded.tsv", excluded, cols)
        counts[name] = (len(rows), len(passed), len(excluded))

    tiers = tier_rows(F, s1, membership, removed)
    for key, rows in tiers.items():
        write(OUT / f"{TIER_FILES[key][0]}.tsv", rows)
    n_core = len(tiers["core"])
    assert n_core == sum(1 for d in committed.values() if d["entered_by"] == "core")

    total = len(membership)
    removed_from_filter = sorted(g for g in removed if g in membership)
    lib = total - len(removed_from_filter)
    lines = [
        "# Genes at each L1 gene-list gate",
        "",
        "Per-gate and per-tier gene tables written by",
        "[`gate_lists.py`](../gate_lists.py) from an offline rerun of the filter over",
        "accessible-surfaceome release 1.3.0 (the pipeline's `gene_list` stage). For each gate,",
        "`*_passed.tsv` holds only the genes remaining after the gate, `*_excluded.tsv` holds only",
        "the genes excluded at that gate (so each file's row count matches the count shown in the",
        "brief), and the combined `XX_<gate>.tsv` holds all genes tested at that gate (`fail` rows",
        "first). Gates apply in order and stop at a gene's first failure, so a gene appears up to",
        "and including the gate that excluded it, and in no later one. The per-gene summary across",
        "all gates is [`gene_gate_outcomes.tsv`](../../gene_gate_outcomes.tsv); the script asserts",
        "every pass/fail here agrees with it. Topology values are DeepTMHMM v1 (1.0.24) on the",
        "UniProt canonical sequence, as recorded in the release.",
        "",
        "| Gate | Genes tested (combined) | Genes remaining (`*_passed.tsv`) | Genes excluded (`*_excluded.tsv`) |",
        "|---|---:|---:|---:|",
    ]
    for name in gate_names:
        n, npass, nf = counts[name]
        stem = GATE_FILES[name]
        lines.append(
            f"| {name} | [{n:,}]({stem}.tsv) | [{npass:,}]({stem}_passed.tsv) | [{nf:,}]({stem}_excluded.tsv) |"
        )
    lines += ["", "| File | Tier | Genes |", "|---|---|---:|"]
    for key, (stem, label) in TIER_FILES.items():
        lines.append(f"| [`{stem}.tsv`]({stem}.tsv) | {label} | {len(tiers[key]):,} |")
    lines += [
        "",
        f"The filter passes {total:,} genes (core {n_core:,} plus the three tiers); "
        f"{len(removed_from_filter)} removed afterwards ([`gene_list_exclusions.tsv`]"
        f"(../../gene_list_exclusions.tsv): {', '.join(removed_from_filter)}) leaves "
        f"{lib:,}, the library's gene count. Cohort: {len(F.recs):,} genes.",
        "",
    ]
    lib_rows, in_lib = library_gene_rows(committed, membership, removed)
    write(OUT / "14_library_genes.tsv", lib_rows)
    lines.insert(
        len(lines) - 1,
        f"[`14_library_genes.tsv`](14_library_genes.tsv) lists all {len(lib_rows):,} library genes: "
        f"the {len(in_lib):,} above and {len(lib_rows) - len(in_lib)} dated additions "
        "([`gene_list_additions.tsv`](../../gene_list_additions.tsv)).\n",
    )
    (OUT / "README.md").write_text("\n".join(lines))

    print(f"cohort {len(F.recs)}; filter passes {total}; library {lib}")
    for name in gate_names:
        print(
            f"  {name}: tested {counts[name][0]}, passed {counts[name][1]}, failed {counts[name][2]}"
        )
    for key in TIER_FILES:
        print(f"  tier {key}: {len(tiers[key])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
