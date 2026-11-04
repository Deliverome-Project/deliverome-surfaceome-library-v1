"""Draw a balanced ortholog panel, binned by the human gene's class.

The isoform panel exists because an unstratified draw returns whatever the library is mostly
made of, and multipass proteins dominate the surfaceome by count. The same holds for
orthologs, with one addition: an ortholog's value for tag placement depends on how far it
has drifted from the human reference, and a panel drawn without regard to that would be
mostly near-identical sequences that exercise nothing.

So the draw is stratified twice. **Class** comes from the human gene's selected reference,
never the ortholog's own topology -- the same rule the isoform panel follows, for the same
reason: an ortholog that happens to lose a signal peptide computes as a different class on
its own sequence while still being an ortholog of a type I receptor. **Correspondence** is
split at MEDIAN_SPLIT into a near band and a diverged band, so the panel carries orthologs
that actually differ from the human form as well as ones that barely do.

Only picks that are surface by topology are eligible, and on a GPI gene only picks that
keep the anchor, so the panel cannot propose an ortholog the library would reject.

**Rejected picks are replaced by rule** (``panel_draw.py``): the draw runs unchanged, then
each pick named in ``--reject`` -- or whose protein is identical to the human reference
(correspondence 1.0, so it tests no drift; cyno MMP16 was drawn this way) -- is replaced by
continuing the same seeded draw within its species, band and class. Replacements are logged
in ``ortholog_panel_replacements.tsv``. ``--snapshot`` reads only the pinned inputs in
``panel_inputs.json``; with no rejections it reproduces the published panel byte for byte.

    uv run --locked python sample_ortholog_panel.py --snapshot --reject panel_rejections.tsv
"""

from __future__ import annotations

import argparse
import csv
import json
import random
from collections import Counter, defaultdict
from pathlib import Path

from gpi import gpi_genes as gpi_gene_set
from panel_draw import load_rejections, replace_rejected
from sample_isoform_panel import FLOORS, allocate, gene_class_map

SEED = 20261011
PER_SPECIES = 50          # per correspondence band, so 100 per species
SURFACE_LABELS = {"TM", "SP+TM", "BETA"}
MEDIAN_SPLIT = "median"   # the band boundary is the measured median, not a fixed number


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--selection", default="reference_selection.tsv")
    ap.add_argument("--cache-dir", default="~/.cache/deliverome/appris")
    ap.add_argument("--per-band", type=int, default=PER_SPECIES)
    ap.add_argument("--labels", default="ortholog_topology_labels.tsv",
                    help="topology per ortholog accession; most picks carry their label on "
                         "D1 rather than in any local run, and drawing only from the local "
                         "ones biases the panel to recently-changed picks")
    ap.add_argument("--out", default="ortholog_panel.tsv")
    ap.add_argument("--snapshot", action="store_true",
                    help="read only the pinned inputs in panel_inputs.json, staged from S3")
    ap.add_argument("--reject", default=None, help="panel_rejections.tsv")
    ap.add_argument("--original", action="store_true",
                    help="the draw as first published: no rejections, no replacement")
    args = ap.parse_args()

    snap = None
    if args.snapshot:
        import tempfile

        from panel_inputs import stage

        snap = stage(tempfile.mkdtemp(prefix="panel-inputs-"))
        args.selection = snap / "reference_selection.tsv"
        args.labels = snap / "ortholog_topology_labels.tsv"
        gpi = gpi_gene_set(snap / "gpi_set.tsv")
    else:
        gpi = gpi_gene_set()
    gene_class = gene_class_map(args.selection, gpi, snapshot_dir=snap)
    rejected = {} if args.original else load_rejections(args.reject, "ortholog")
    replaced: list[dict] = []

    def form_id(r: dict) -> str:
        return f"{r['species']}:{r['human_symbol']}:{r['picked_acc']}"

    def invalid(r: dict) -> str:
        if args.original:
            return ""
        return "protein identical to human" if float(r["picked_score"]) >= 1.0 else ""
    label_of: dict[str, str] = {}
    lp = Path(args.labels)
    if lp.exists():
        with lp.open() as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                label_of[row["accession"]] = row["deeptmhmm_label"]

    panel, summary = [], {}
    for species in ("mouse", "cyno"):
        src = (snap or Path(".")) / f"ortholog_repick_{species}.tsv"
        rows = list(csv.DictReader(open(src), delimiter="\t"))
        eligible = []
        for r in rows:
            cls = gene_class.get(r["human_symbol"])
            if cls is None:
                continue
            acc = r["picked_acc"]
            label = (r["picked_topology_label"] or label_of.get(acc)
                     or label_of.get(acc.split("-")[0]) or "")
            anchored = r["picked_anchor_call"] == "anchored"
            # Surface by topology, or anchored on a GPI gene where the helix count is zero
            # by construction.
            if not (label in SURFACE_LABELS or anchored):
                continue
            eligible.append({**r, "topology_class": cls, "species": species,
                             "topology_label": label,
                             "correspondence": float(r["picked_score"])})
        if not eligible:
            continue
        scores = sorted(x["correspondence"] for x in eligible)
        split = scores[len(scores) // 2]
        bands = {"near": [x for x in eligible if x["correspondence"] >= split],
                 "diverged": [x for x in eligible if x["correspondence"] < split]}
        summary[species] = {
            "eligible": len(eligible),
            "correspondence_split": round(split, 4),
            "classes": dict(Counter(x["topology_class"] for x in eligible).most_common()),
            "bands": {},
        }
        rng = random.Random(SEED + (0 if species == "mouse" else 1))
        band_strata: dict[str, tuple] = {}
        for band, pool in bands.items():
            want = min(args.per_band, len(pool))
            alloc = allocate(pool, want)
            by_cls = defaultdict(list)
            for x in pool:
                by_cls[x["topology_class"]].append(x)
            strata = [(cls, rng.sample(by_cls[cls], min(k, len(by_cls[cls]))))
                      for cls, k in alloc.items()]
            band_strata[band] = (pool, by_cls, strata)
        # Replacement runs after both bands are drawn, so the valid picks never move.
        used = {form_id(x) for _, _, strata in band_strata.values()
                for _, ps in strata for x in ps}
        for band, (pool, by_cls, strata) in band_strata.items():
            picked = []
            for cls, ps in strata:
                final, log = replace_rejected(rng, ps, by_cls[cls], form_id, rejected,
                                              invalid, used)
                picked += final
                replaced += [{"species": species, "band": band, "topology_class": cls, **row}
                             for row in log]
            for x in picked:
                panel.append({
                    "species": species, "band": band,
                    "human_symbol": x["human_symbol"],
                    "ortholog_symbol": x["ortholog_symbol"],
                    "picked_acc": x["picked_acc"],
                    "correspondence": x["picked_score"],
                    "topology_label": x["topology_label"],
                    "anchor_call": x["picked_anchor_call"],
                    "topology_class": x["topology_class"],
                })
            summary[species]["bands"][band] = {
                "pool": len(pool), "drawn": len(picked),
                "by_class": dict(Counter(p["topology_class"] for p in picked).most_common()),
            }

    panel.sort(key=lambda r: (r["species"], r["band"], r["human_symbol"]))
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(panel[0]), delimiter="\t")
        w.writeheader()
        w.writerows(panel)
    summary["seed"] = SEED
    summary["per_band"] = args.per_band
    summary["floors"] = FLOORS
    summary["total"] = len(panel)
    summary["replaced"] = len(replaced)
    with open("ortholog_panel_replacements.tsv", "w", newline="") as fh:
        cols = ["species", "band", "topology_class", "rejected", "reason", "replacement",
                "candidates"]
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(replaced)
    for row in replaced:
        print(f"  replaced {row['rejected']} ({row['reason']}) -> {row['replacement'] or 'NONE'}")
    Path("ortholog_panel_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"wrote {args.out} ({len(panel)} orthologs)")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
