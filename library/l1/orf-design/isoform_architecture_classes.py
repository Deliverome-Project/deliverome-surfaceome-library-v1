"""Membrane-architecture classes in the library, the isoform cohort and the isoform panel,
and how far the class-stratified draw departs from a uniform one.

Becca, 2026-11-20, on "(type I, type II, type III, GPI, multipass)": "isn't TIII also
multipass? what's the initial distribution of these classes and any other classes to
consider? what is the enrichment vs baseline of uniform?"

**The classes are the isoform sampler's** (``sample_isoform_panel.classify``), assigned per
GENE from the gene's selected reference, never from the isoform's own sequence:

| Class | Rule, in order |
|---|---|
| GPI | the gene is in the GPI set (``gpi_set.tsv``), whatever its helix count |
| multipass N-out | 2 or more TM helices, N-terminus extracellular |
| multipass N-in | 2 or more TM helices, N-terminus cytoplasmic |
| type I | 1 TM helix, a signal peptide, N-terminus extracellular |
| type II | 1 TM helix, no signal peptide, N-terminus cytoplasmic |
| type III | 1 TM helix, no signal peptide, N-terminus extracellular |
| single-pass other | 1 TM helix in any other combination (a signal peptide with a cytoplasmic N-terminus) |
| no TM helix | no TM helix and not GPI; excluded from the panel |

Type III here is UniProt's "single-pass type III" (a signal-anchor protein that puts its
N-terminus outside), so it is single-pass by definition; some textbooks use "type III" for
multipass proteins instead, which is not this usage. The topology is the one the sampler read
when the panel was drawn: DeepTMHMM v1 from the April accessible-surfaceome tables, looked up
by the reference's UniProt accession, or the v1 run on references those tables lacked
(``library_reference_missing_topology.3line``), all from the pinned panel inputs.

Beside it, for the classes worth considering, each gene's BUILT reference protein (the active
L1 ``library_reference`` member) is classified the same way on its own DeepTMHMM v1 string,
matched by sequence, and split further:

* multipass by N-terminal side (outside, e.g. GPCRs, or inside), which decides whether a
  residue-1 tag can be extracellular;
* type II with a C-terminal tail of at most 30 residues after the helix (tail-anchored-like;
  a tail anchor inserts post-translationally, so its tag behaviour can differ from a type II
  signal anchor).

**Since 2026-11-22** (Becca: "yes let's split this up into two"; "yes we should always rebuild
with the latest") the panel is the applied redraw, so the counts follow it: multipass is the two
strata the redraw draws (``redraw_isoform_panel_2026_11_22.split_classify``), the cohort is the
rebuilt cohort (``isoform_cohort_rebuild.tsv``), the gate is the gate that redraw read
(``isoform_gate_rebuild.tsv``, every form judged on its own DeepTMHMM v1 prediction) and the
panel is ``isoform_panel.tsv`` as written by ``write_panels_of_record.py`` (the redraw plus the
six refills). The earlier cohort-of-record and full-v1 gate columns are gone with the old
cohort. Updated by Claude Code session l1-review-export-scope (brief-reconciliation fork).

"Uniform" is the draw the sampler would make without strata: 100 forms uniformly from the
observed tier and 50 from the prediction-only tier, so a class's expected count is
``100 * p_observed + 50 * p_prediction_only``. Enrichment is the panel count over that
expectation, against both the cohort the allocation is computed from and the gate-passing
forms a uniform draw could actually use.

    uv run --locked python isoform_architecture_classes.py
"""

from __future__ import annotations

import csv
import json
import sys
import tempfile
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
AP = REPO / "library/l1/proteoform-selection"
sys.path.insert(0, str(AP))
sys.path.insert(0, str(HERE))

OUT = HERE / "isoform_architecture_classes.tsv"
OUT_GENES = HERE / "isoform_architecture_classes_genes.tsv"
SUMMARY = HERE / "isoform_architecture_classes_summary.json"
FIGURE = HERE / "figures" / "isoform_architecture_classes"
CLASSES = [
    "type I",
    "type II",
    "type III",
    "GPI",
    "multipass N-out",
    "multipass N-in",
    "single-pass other",
    "no TM helix",
    "unclassified",
]
TAIL_ANCHOR_MAX = 30
N_OBSERVED, N_PREDICTED = 100, 50


def read(path: Path) -> list[dict]:
    with Path(path).open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def v1_features(topo: str) -> dict:
    runs, j = [], 0
    while j < len(topo):
        if topo[j] == "M":
            k = j
            while k < len(topo) and topo[k] == "M":
                k += 1
            runs.append((j, k))
            j = k
        else:
            j += 1
    non_signal = [c for c in topo if c != "S"]
    first = next((c for c in non_signal if c in "IO"), None)
    return {
        "tm": len(runs),
        "sp": topo.count("S"),
        "n_term": {"O": "extracellular", "I": "cytoplasmic"}.get(first, "indeterminate"),
        "c_tail": len(topo) - runs[-1][1] if runs else "",
    }


def built_reference_rows(genes, v1: dict, cls: dict, gpi, classify) -> list[dict]:
    """Per gene: the panel class beside the class of its built reference's own v1 string."""
    gene_rows = []
    for g in genes:
        f = v1_features(v1[g]) if g in v1 else None
        own = classify(f["tm"], f["sp"], f["n_term"], g in gpi) if f else "no built reference"
        sub = ""
        if f and own == "multipass":
            sub = f"multipass, N-terminus {f['n_term']}"
        elif f and own == "type II" and f["c_tail"] != "" and f["c_tail"] <= TAIL_ANCHOR_MAX:
            sub = f"type II, C-tail <= {TAIL_ANCHOR_MAX} aa (tail-anchored-like)"
        gene_rows.append(
            {
                "gene": g,
                "panel_class": cls[g],
                "built_reference_class_v1": own,
                "built_reference_subclass": sub,
                "tm_helices_v1": f["tm"] if f else "",
                "signal_peptide_length_v1": f["sp"] if f else "",
                "n_terminus_v1": f["n_term"] if f else "",
                "c_tail_after_last_tm": f["c_tail"] if f else "",
            }
        )
    return gene_rows


def write_tsv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def main() -> int:
    import rebuild_isoform_cohort as rb
    import sample_isoform_panel as s
    from gpi import gpi_genes
    from panel_inputs import stage
    from redraw_isoform_panel_2026_11_22 import split_classify

    from deliverome_analysis import library_design as ld
    from deliverome_analysis import topology

    work = Path(tempfile.mkdtemp(prefix="isoform-classes-"))
    snap = stage(work / "panel")
    gpi = gpi_genes(snap / "gpi_set.tsv")
    by_acc, by_ref_gene = s._topology_sources(
        reference_three_line=snap / "library_reference_missing_topology.3line",
        deeptmhmm_dir=snap,
    )
    base = s.classify
    s.classify = split_classify(base)
    try:
        gene_class = s._gene_class(snap / "reference_selection.tsv", gpi, by_acc, by_ref_gene)
    finally:
        s.classify = base
    genes = rb.library_genes()
    cls = {g: gene_class.get(g, "unclassified") for g in genes}

    # The built reference, classified on its own v1 string (sequence-matched).
    manifest = HERE / "candidate_manifest_tranches.tsv"
    refs = {
        m["gene_symbol"]: m["sequence"]
        for m in ld.designed_proteins(manifest, library=ld.LIBRARY_L1).values()
        if m["candidate_set"] == "library_reference"
    }
    rep = topology.resolve(refs, tool=topology.TOOL_V1)
    v1 = {r["key"]: r["per_residue_topology"] for r in rep.to_dict("records") if r["found"]}
    gene_rows = built_reference_rows(genes, v1, cls, gpi, s.classify)

    cohort = read(rb.OUT_COHORT)
    gate = {r["form_id"]: r["verdict"] for r in read(rb.OUT_GATE)}
    panel = read(HERE / "isoform_panel.tsv")
    fid = lambda r: f"{r['gene']}:{r['uniprot'].strip() or r['ensembl_transcripts']}"  # noqa: E731
    tiers = ("observed", "full_only")

    def counts(rows, key=lambda r: cls.get(r["gene"], "unclassified")) -> Counter:
        return Counter(key(r) for r in rows)

    lib = Counter(cls.values())
    coh = {t: counts([r for r in cohort if r["cohort_tier"] == t]) for t in tiers}
    drawable = {
        t: counts([r for r in cohort if r["cohort_tier"] == t and gate.get(fid(r)) == "pass"])
        for t in tiers
    }
    pan_tier = {
        t: Counter(r["topology_class"] for r in panel if r["panel_tier"] == t)
        for t in ("observed", "prediction_only")
    }
    pan = Counter(r["topology_class"] for r in panel)
    pan_check = counts(panel)
    assert pan == pan_check, "panel classes differ from the gene class map"
    coh_genes = Counter(cls[g] for g in {r["gene"] for r in cohort})
    pan_genes = Counter(cls[g] for g in {r["gene"] for r in panel})

    def expected(pools: dict[str, Counter], c: str) -> float:
        # Uniform within each tier, over the classes the sampler can draw (it drops no-TM).
        out = 0.0
        for t, n in (("observed", N_OBSERVED), ("full_only", N_PREDICTED)):
            usable = sum(v for k, v in pools[t].items() if k not in ("no TM helix", "unclassified"))
            if c not in ("no TM helix", "unclassified") and usable:
                out += n * pools[t][c] / usable
        return out

    rows = []
    for c in CLASSES:
        e_coh, e_draw = expected(coh, c), expected(drawable, c)
        rows.append(
            {
                "class": c,
                "library_genes": lib[c],
                "library_share": round(lib[c] / len(genes), 4),
                "cohort_genes": coh_genes[c],
                "cohort_forms_observed": coh["observed"][c],
                "cohort_forms_prediction_only": coh["full_only"][c],
                "cohort_forms": coh["observed"][c] + coh["full_only"][c],
                "drawable_forms_gate_of_record": drawable["observed"][c] + drawable["full_only"][c],
                "panel_forms": pan[c],
                "panel_forms_observed": pan_tier["observed"][c],
                "panel_forms_prediction_only": pan_tier["prediction_only"][c],
                "panel_genes": pan_genes[c],
                "uniform_expected_from_cohort": round(e_coh, 1),
                "enrichment_vs_uniform_cohort": round(pan[c] / e_coh, 2) if e_coh else "",
                "uniform_expected_from_drawable": round(e_draw, 1),
                "enrichment_vs_uniform_drawable": round(pan[c] / e_draw, 2) if e_draw else "",
            }
        )
    write_tsv(OUT, rows)
    write_tsv(OUT_GENES, gene_rows)

    agree = Counter((r["panel_class"], r["built_reference_class_v1"]) for r in gene_rows)
    summary = {
        "library_genes": len(genes),
        "cohort_forms": len(cohort),
        "panel_forms": len(panel),
        "class_rows": {r["class"]: {k: v for k, v in r.items() if k != "class"} for r in rows},
        "built_reference_subclasses": dict(
            Counter(
                r["built_reference_subclass"] for r in gene_rows if r["built_reference_subclass"]
            )
        ),
        "panel_class_vs_built_reference_class_v1": {
            "same": sum(n for (a, b), n in agree.items() if a == b),
            "different": {f"{a} -> {b}": n for (a, b), n in agree.most_common() if a != b},
        },
        "panel_genes_by_subclass": dict(
            Counter(
                r["built_reference_subclass"]
                for r in gene_rows
                if r["built_reference_subclass"] and r["gene"] in {p["gene"] for p in panel}
            )
        ),
    }
    SUMMARY.write_text(json.dumps(summary, indent=2) + "\n")
    plot(rows)
    for r in rows:
        print(r)
    print(json.dumps(summary["built_reference_subclasses"], indent=1))
    return 0


def plot(rows: list[dict]) -> None:
    import matplotlib.pyplot as plt
    import numpy as np

    from deliverome_analysis.plotting import (
        get_categorical_colors,
        save_figure,
        setup_plotting_style,
    )

    setup_plotting_style()
    drawn = ["type I", "type II", "type III", "GPI", "multipass N-out", "multipass N-in"]
    labels = ["type I", "type II", "type III", "GPI", "multipass\nN-out", "multipass\nN-in"]
    data = {r["class"]: r for r in rows}
    series = [
        ("Library genes", [data[c]["library_genes"] for c in drawn]),
        ("Cohort forms", [data[c]["cohort_forms"] for c in drawn]),
        ("Gate-passing forms", [data[c]["drawable_forms_gate_of_record"] for c in drawn]),
        ("Panel forms", [data[c]["panel_forms"] for c in drawn]),
    ]
    colors = get_categorical_colors(len(series))
    fig, (a, b) = plt.subplots(1, 2, figsize=(16, 6.5), gridspec_kw={"width_ratios": [3, 2]})
    x = np.arange(len(drawn))
    width = 0.2
    for i, (name, vals) in enumerate(series):
        tot = sum(vals)
        share = [100 * v / tot for v in vals]
        bars = a.bar(
            x + (i - 1.5) * width, share, width, label=f"{name} (n = {tot:,})", color=colors[i]
        )
        for bar, v in zip(bars, vals, strict=True):
            a.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + 0.8,
                f"{v:,}",
                ha="center",
                va="bottom",
                fontsize=9,
                rotation=90,
            )
    a.set_xticks(x, labels)
    a.set_ylabel("Share of the six drawn classes (%)")
    a.set_ylim(0, 75)
    a.legend(frameon=False, fontsize=11, loc="upper left")
    a.text(-0.08, 1.02, "a", transform=a.transAxes, fontsize=18, fontweight="bold")

    enr = [data[c]["enrichment_vs_uniform_cohort"] for c in drawn]
    enr_d = [data[c]["enrichment_vs_uniform_drawable"] for c in drawn]
    b.barh(x + 0.18, enr, 0.34, color=colors[1], label="vs uniform over the cohort")
    b.barh(x - 0.18, enr_d, 0.34, color=colors[2], label="vs uniform over gate-passing forms")
    b.axvline(1, color="0.3", lw=1, ls="--")
    for xs, vals in ((x + 0.18, enr), (x - 0.18, enr_d)):
        for yy, v in zip(xs, vals, strict=True):
            b.text(v + 0.04, yy, f"{v:.2f}", va="center", fontsize=10)
    b.set_xlim(0, 3.2)
    b.set_yticks(x, labels)
    b.invert_yaxis()
    b.set_xlabel("Panel count / uniform expectation")
    b.legend(frameon=False, fontsize=11, loc="lower right")
    b.text(-0.2, 1.02, "b", transform=b.transAxes, fontsize=18, fontweight="bold")
    fig.tight_layout()
    save_figure(fig, str(FIGURE))
    plt.close(fig)


if __name__ == "__main__":
    raise SystemExit(main())
