"""s1e2 — AlphaMissense-weighted severity for the 2,243 canonical surfaceome ORFs.

Replaces the flat "32% of missense is damaging" assumption with a **per-gene**
damaging fraction, so the synthesis-error model can be weighted by which genes are
actually intolerant.

⚠️ LICENSE: AlphaMissense predictions are CC BY-NC-SA 4.0 (NON-COMMERCIAL).
   Zenodo record 8208688. Commercial use / publication likely needs legal review.
   Permissive alternatives if that is a blocker: ESM-1v (MIT), SIFT, PolyPhen-2.

Input (1.2 GB source, not committed — regenerate as below):
    /tmp/am_ours2.tsv  rows of AlphaMissense_aa_substitutions.tsv.gz filtered to our
                       2,242 UniProt accessions. Columns:
                       uniprot_id, protein_variant, am_pathogenicity, am_class

    # 1. download to disk with resume, and VERIFY — a truncated stream fails silently
    #    downstream ("gunzip: unexpected end of file" after hours of filtering).
    curl -L --retry 8 -C - -o /tmp/am.tsv.gz \
      "https://zenodo.org/records/8208688/files/AlphaMissense_aa_substitutions.tsv.gz?download=1"
    gzip -t /tmp/am.tsv.gz          # expect silence; 1,207,278,510 bytes

    # 2. filter with awk on an exact first-field key — NOT `grep -F -f`.
    #    grep with 2,242 patterns does unanchored substring matching against every line
    #    of 216M lines: measured ~30 accessions in 30 minutes. awk finishes in ~5.
    gunzip -c /tmp/am.tsv.gz \
      | awk -F'\t' 'NR==FNR{a[$1]=1;next} ($1 in a)' canonical_accessions.txt - \
      > /tmp/am_ours2.tsv          # 28,366,506 rows, 2,235 of 2,243 accessions

Run:
    uv run python library/l1/vector-design/alphamissense_severity.py
"""

import csv
import json
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

from deliverome_analysis.plotting import (
    COLORS,
    get_sequential_palette,
    save_figure,
    setup_plotting_style,
)

HERE = Path(__file__).resolve().parent
AM_TSV = Path("/tmp/am_ours2.tsv")

# error model constants (mirror synthesis_cost_analysis.py)
F_MIS, F_NON = 392 / 549, 23 / 549
INDEL_F = 0.30
P_POOL, P_ABIO, P_FRAG = 0.00647, 0.01000, 0.07
POOL_CAP, ABIO_CAP = 1800, 1200
E_POOL, E_ABIO, E_FRAG = 1 / 3000, 1 / 30000, 1 / 7500


def load_am():
    """uniprot_acc -> (n_variants, n_pathogenic, mean_score)."""
    if not AM_TSV.exists():
        raise SystemExit(
            f"{AM_TSV} not found — run the curl + awk recipe in this module's docstring."
        )
    n = defaultdict(int)
    npath = defaultdict(int)
    score = defaultdict(float)
    classes = defaultdict(int)
    with AM_TSV.open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            acc, _variant, path, cls = parts[0], parts[1], parts[2], parts[3].strip()
            try:
                p = float(path)
            except ValueError:
                continue
            n[acc] += 1
            score[acc] += p
            classes[cls] += 1
            # AlphaMissense am_class values are: pathogenic | benign | ambiguous
            if cls.startswith("pathogenic") or cls == "likely_pathogenic":
                npath[acc] += 1
    print(f"  am_class counts: {dict(classes)}")
    return {a: (n[a], npath[a], score[a] / n[a]) for a in n if n[a]}


def main() -> None:
    setup_plotting_style()
    am = load_am()
    print(f"AlphaMissense: {len(am):,} accessions with variant-level data")

    topo = json.loads((HERE / "canonical_topology.json").read_text())
    accs = json.loads((HERE / "canonical_accessions.json").read_text())  # gene -> acc
    rows = []
    for r in topo:
        acc = accs.get(r["gene"])
        rec = am.get(acc) if acc else None
        if not rec:
            continue
        nv, npath, mean_p = rec
        rows.append(
            dict(
                gene=r["gene"],
                acc=acc,
                cds=r["len"] * 3 + 3,
                label=r["label"],
                tm=r["tm"],
                n_variants=nv,
                frac_pathogenic=npath / nv,
                mean_pathogenicity=mean_p,
            )
        )
    print(f"joined to canonical set: {len(rows):,} of {len(topo):,} genes")
    if not rows:
        raise SystemExit("no overlap — check accession formatting")

    fr = np.array([r["frac_pathogenic"] for r in rows])
    print(f"\nfrac_pathogenic across genes: mean={fr.mean():.3f}  median={np.median(fr):.3f}")
    print(f"  p10={np.percentile(fr, 10):.3f}  p90={np.percentile(fr, 90):.3f}")
    print("  (flat assumption used previously: 0.320)")

    # by topology class
    print("\nby DeepTMHMM class:")
    for lab in sorted({r["label"] for r in rows}):
        sub = [r["frac_pathogenic"] for r in rows if r["label"] == lab]
        print(f"  {lab:<7} n={len(sub):>5}  mean frac_pathogenic={np.mean(sub):.3f}")

    # re-run the severe-error model per gene with its own severity
    def route(L):
        if L <= ABIO_CAP:
            return P_ABIO, E_ABIO
        if L <= POOL_CAP:
            return P_POOL, E_POOL
        return P_FRAG, E_FRAG

    flat_sev = weighted_sev = 0.0
    for r in rows:
        L = r["cds"]
        _, err = route(L)
        for sev, key in ((0.32, "flat"), (r["frac_pathogenic"], "w")):
            p = err * ((1 - INDEL_F) * (F_NON + F_MIS * sev) + INDEL_F)
            v = 1 - (1 - p) ** L
            if key == "flat":
                flat_sev += v
            else:
                weighted_sev += v
    n = len(rows)
    print("\n=== severe-error rate, scenario D routing ===")
    print(f"  flat 32%              : {100 * flat_sev / n:.1f}%")
    print(f"  AlphaMissense-weighted: {100 * weighted_sev / n:.1f}%  (each gene its own fraction)")

    # Sensitivity sweep. The cutoff behind `am_class = pathogenic` was chosen for 90%
    # precision on ClinVar, which is a clinical endpoint, not "does it still reach the
    # surface". Sweep the whole defensible range so the reader can see how little it
    # matters — the indel + nonsense terms dominate and don't depend on this parameter.
    print("\n=== sensitivity to the missense-severity parameter (uniform) ===")
    sweep = [
        (0.070, "Taipale-chained: 0.427 x ~1/6 pathogenic missense mislocalise"),
        (0.245, "VAMP-seq experimental (PTEN + TPMT abundance)"),
        (0.320, "AlphaMissense genome-wide `pathogenic`"),
        (0.427, "AlphaMissense on OUR genes, `pathogenic` class"),
        (0.567, "AlphaMissense counting `ambiguous` as pathogenic too"),
    ]
    for sev, label in sweep:
        tot = 0.0
        for r in rows:
            L = r["cds"]
            _, err = route(L)
            p = err * ((1 - INDEL_F) * (F_NON + F_MIS * sev) + INDEL_F)
            tot += 1 - (1 - p) ** L
        print(f"  {sev:.3f} -> {100 * tot / n:>5.1f}%   {label}")

    out = HERE / "alphamissense_by_gene.tsv"
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(sorted(rows, key=lambda r: -r["frac_pathogenic"]))
    print(f"wrote {out}")

    # ---------- figure ----------
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.2))
    ax = axes[0]
    ax.hist(fr, bins=40, color=COLORS["primary"], alpha=0.85)
    ax.axvline(0.32, ls="--", color="0.35", lw=1.6)
    ax.annotate(
        "flat assumption 0.32",
        xy=(0.32, ax.get_ylim()[1] * 0.92),
        xytext=(6, 0),
        textcoords="offset points",
        fontsize=9,
        color="0.35",
    )
    ax.axvline(fr.mean(), color=COLORS["secondary"], lw=2)
    ax.annotate(
        f"observed mean {fr.mean():.2f}",
        xy=(fr.mean(), ax.get_ylim()[1] * 0.75),
        xytext=(6, 0),
        textcoords="offset points",
        fontsize=9,
        color=COLORS["secondary"],
    )
    ax.set_xlabel("fraction of possible missense classed likely-pathogenic")
    ax.set_ylabel("genes")

    ax = axes[1]
    cds = np.array([r["cds"] for r in rows])
    # hexbin, not scatter: 2,236 overlapping markers both obscure the density and blow
    # the SVG past Notion's 200 KiB inline-attachment limit.
    cmap = LinearSegmentedColormap.from_list(
        "dlv_teal", ["#ffffff", *get_sequential_palette("teal", 6)]
    )
    hb = ax.hexbin(np.log10(cds), fr, gridsize=32, cmap=cmap, mincnt=1, linewidths=0.2)
    fig.colorbar(hb, ax=ax, label="genes per bin", pad=0.02)
    ax.set_xlabel("CDS length (bp)")
    ax.set_ylabel("fraction likely-pathogenic")
    ticks = [300, 1000, 3000, 10000, 40000]
    ax.set_xticks([np.log10(t) for t in ticks])
    ax.set_xticklabels([f"{t:,}" for t in ticks])
    ax.axhline(0.32, ls="--", color="0.35", lw=1.4)

    save_figure(fig, "s1e2_alphamissense_severity", output_dir=str(HERE / "figures"))
    print("wrote figures/s1e2_alphamissense_severity.{svg,png}")


if __name__ == "__main__":
    main()
