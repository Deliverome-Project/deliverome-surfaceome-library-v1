"""Draw a balanced isoform panel: 100 from the observed tier, 50 from prediction-only.

Random sampling would hand back whatever the library is mostly made of. Multipass proteins
dominate the surfaceome by count, so an unstratified draw buries the classes a tag-placement
pipeline most needs to be exercised on — a GPI anchor behaves differently from a type I
receptor, and type II and type III differ from each other in exactly the way that decides
where a tag can go.

**Isoforms are binned by the class of the library gene they belong to, not by their own
topology.** An isoform of a type I receptor that happens to lose its signal peptide computes
as type III on its own sequence, and binning on that produced a panel with 20 type III
members when the library contains **13 type III genes in total**. The gene's class is the
property the panel is meant to span.

Classes, from the gene's selected reference:

* **GPI** — UniProt keyword KW-0336, taken first because a GPI protein's C-terminal signal
  is often called as a transmembrane helix and would otherwise be filed as type I;
* **multipass** — two or more TM helices;
* **type I** — one TM helix, signal peptide, extracellular N-terminus;
* **type II** — one TM helix, no signal peptide, cytoplasmic N-terminus;
* **type III** — one TM helix, no signal peptide, extracellular N-terminus (signal anchor).

Deterministic: the seed is fixed and recorded, so the panel is reproducible and a rerun
does not quietly hand back a different 150. ``--snapshot`` reads only the pinned inputs in
``panel_inputs.json``; with no rejections it reproduces the published panel byte for byte.

**Rejected picks are replaced by rule** (``panel_draw.py``): the draw runs unchanged, then
each pick named in ``--reject`` -- or that is its gene's selected reference form, which the
cohort table could not know once a GPI override moved the reference -- is replaced by
continuing the same seeded draw within its tier and class. Replacements are logged in
``isoform_panel_replacements.tsv``.

**The cohort gate** (Becca, 2026-11-12, "Gate cohort, redraw"; s4e9
``gate_isoform_cohort.py``): with ``--gate isoform_gate.tsv`` every form the gate fails is
invalid, so a failing original pick is replaced and a failing form is never drawn as a
replacement. Under ``--snapshot`` the table must match the hash pinned in
``panel_inputs.json`` (``isoform_gate``); ``--original`` ignores it.

    uv run --locked python sample_isoform_panel.py --snapshot
    uv run --locked python sample_isoform_panel.py --snapshot --reject panel_rejections.tsv \
        --gate isoform_gate.tsv
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

SEED = 20261009
N_OBSERVED = 100
N_PREDICTED = 50
# Floors for the classes too rare to survive a proportional draw. The library holds 33 type
# III, 126 type II and 135 GPI genes against 1,223 multipass and 836 type I, so without
# floors the two big classes take almost everything. Everything else is proportional.
FLOORS = {"type III": 2, "GPI": 5, "type II": 5}

ORIENT = {"O": "extracellular", "I": "cytoplasmic"}


def classify(tm: int, sp: int, n_term: str, is_gpi: bool) -> str:
    if is_gpi:
        return "GPI"
    if tm >= 2:
        return "multipass"
    if tm == 1:
        if sp > 0 and n_term == "extracellular":
            return "type I"
        if sp == 0 and n_term == "cytoplasmic":
            return "type II"
        if sp == 0 and n_term == "extracellular":
            return "type III"
        return "single-pass other"
    return "no TM helix"


def _topology_sources(repo="~/Git/accessible-surfaceome",
                      reference_three_line="library_reference_missing_topology.3line",
                      deeptmhmm_dir=None):
    """Topology by accession and, for references run separately, by gene."""
    repo = Path(repo).expanduser()
    tables = Path(deeptmhmm_dir) if deeptmhmm_dir else repo / "data/processed/deeptmhmm"
    args = type("A", (), {"reference_three_line": reference_three_line})
    # Topology per form: the cohort tables for UniProt-keyed forms, the new run for the rest.
    by_acc: dict[str, dict] = {}
    for name in ("deeptmhmm_human_canonical.tsv", "deeptmhmm_human_isoforms.tsv"):
        path = tables / name
        if not path.exists():
            continue
        with path.open() as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                acc = row.get("uniprot_accession_full") or row.get("uniprot_accession")
                if acc:
                    by_acc[acc] = {
                        "tm": int(row["tm_helix_count"] or 0),
                        "sp": int(row["signal_peptide_length"] or 0),
                        "n_term": "extracellular"
                        if row.get("n_term_extracellular") == "1"
                        else "cytoplasmic",
                    }
    # References run separately because they were missing from the cohort.
    by_ref_gene: dict[str, dict] = {}
    ref3 = Path(args.reference_three_line)
    if ref3.exists():
        lines = ref3.read_text().splitlines()
        for i in range(len(lines) - 2):
            if not lines[i].startswith(">"):
                continue
            name, _label = lines[i][1:].rsplit("|", 1)
            parts = name.split("|")
            gene = parts[1] if len(parts) > 1 else name
            topo = lines[i + 2].strip()
            runs = sum(
                1 for j, c in enumerate(topo) if c == "M" and (j == 0 or topo[j - 1] != "M")
            )
            non_signal = [c for c in topo if c != "S"]
            first = next((c for c in non_signal if c in ORIENT), None)
            by_ref_gene[gene.strip()] = {
                "tm": runs,
                "sp": topo.count("S"),
                "n_term": ORIENT.get(first, "indeterminate"),
            }
    return by_acc, by_ref_gene


def _gene_class(selection_path, gpi_genes, by_acc, by_ref_gene) -> dict[str, str]:
    """Class of each library gene, read off its selected reference."""
    out: dict[str, str] = {}
    with open(selection_path) as fh:
        for srow in csv.DictReader(fh, delimiter="\t"):
            gene = srow["symbol"]
            acc = srow.get("selected_uniprot_equivalent") or srow.get("current_accession") or ""
            topo = by_acc.get(acc) or by_acc.get(acc.split("-")[0]) or by_ref_gene.get(gene)
            if topo is not None:
                out[gene] = classify(
                    topo["tm"], topo["sp"], topo["n_term"], gene in gpi_genes
                )
    return out


def gene_class_map(selection_path, gpi_genes, repo="~/Git/accessible-surfaceome",
                   snapshot_dir=None):
    """Public entry point: build the gene -> class map from the cached topology sources."""
    if snapshot_dir:
        by_acc, by_ref_gene = _topology_sources(
            reference_three_line=Path(snapshot_dir) / "library_reference_missing_topology.3line",
            deeptmhmm_dir=snapshot_dir)
    else:
        by_acc, by_ref_gene = _topology_sources(repo)
    return _gene_class(selection_path, gpi_genes, by_acc, by_ref_gene)


def allocate(pool: list[dict], n: int) -> dict[str, int]:
    """Floors for type III and GPI; the remainder split by largest remainder.

    Largest-remainder rather than rounding each share independently: rounding lets the
    big classes absorb the whole remainder and starve the small ones, which took type II
    from its proportional ~6 slots to 0 despite having 47 candidates.
    """
    by_cls = defaultdict(list)
    for r in pool:
        by_cls[r["topology_class"]].append(r)
    alloc = {c: min(FLOORS.get(c, 0), len(v)) for c, v in by_cls.items()}
    left = n - sum(alloc.values())
    spare = {c: len(v) - alloc[c] for c, v in by_cls.items()}
    total = sum(spare.values())
    if left > 0 and total > 0:
        exact = {c: left * spare[c] / total for c in spare}
        base = {c: min(int(v), spare[c]) for c, v in exact.items()}
        for c, k in base.items():
            alloc[c] += k
        short = left - sum(base.values())
        # hand out the leftover one at a time, largest fractional part first
        order = sorted(
            spare, key=lambda c: (exact[c] - int(exact[c]), spare[c]), reverse=True
        )
        i = 0
        while short > 0 and i < 50 * max(len(order), 1):
            c = order[i % len(order)]
            if alloc[c] - FLOORS.get(c, 0) < spare[c]:
                alloc[c] += 1
                short -= 1
            i += 1
    while sum(alloc.values()) > n:
        c = max(alloc, key=lambda c: alloc[c] - FLOORS.get(c, 0))
        alloc[c] -= 1
    return {c: k for c, k in alloc.items() if k}



def parse_3line(path: Path) -> dict[str, dict]:
    out: dict[str, dict] = {}
    if not path.exists():
        return out
    lines = path.read_text().splitlines()
    for i in range(len(lines) - 2):
        if not lines[i].startswith(">"):
            continue
        name, label = lines[i][1:].rsplit("|", 1)
        topo = lines[i + 2].strip()
        runs = sum(
            1 for j, c in enumerate(topo) if c == "M" and (j == 0 or topo[j - 1] != "M")
        )
        non_signal = [c for c in topo if c != "S"]
        first = next((c for c in non_signal if c in ORIENT), None)
        out[name.strip()] = {
            "tm": runs,
            "sp": topo.count("S"),
            "n_term": ORIENT.get(first, "indeterminate"),
        }
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--surfaceome-repo", default="~/Git/accessible-surfaceome")
    ap.add_argument("--cache-dir", default="~/.cache/deliverome/appris")
    ap.add_argument("--cohorts", default="isoform_cohorts.tsv")
    ap.add_argument("--selection", default="reference_selection.tsv")
    ap.add_argument(
        "--reference-three-line", default="library_reference_missing_topology.3line"
    )
    ap.add_argument("--snapshot", action="store_true",
                    help="read only the pinned inputs in panel_inputs.json, staged from S3")
    ap.add_argument("--reject", default=None, help="panel_rejections.tsv")
    ap.add_argument("--gate", default=None,
                    help="isoform_gate.tsv (gate_isoform_cohort.py): failing forms are invalid")
    ap.add_argument("--original", action="store_true",
                    help="the draw as first published: no rejections, no replacement")
    args = ap.parse_args()

    if args.snapshot:
        import tempfile

        from panel_inputs import stage

        snap = stage(tempfile.mkdtemp(prefix="panel-inputs-"))
        args.cohorts = snap / "isoform_cohorts.tsv"
        args.selection = snap / "reference_selection.tsv"
        gpi_genes = gpi_gene_set(snap / "gpi_set.tsv")
        by_acc, by_ref_gene = _topology_sources(
            reference_three_line=snap / "library_reference_missing_topology.3line",
            deeptmhmm_dir=snap)
    else:
        gpi_genes = gpi_gene_set()
        by_acc, by_ref_gene = _topology_sources(args.surfaceome_repo,
                                                args.reference_three_line)

    # The class of each library GENE, from its selected reference.
    gene_class = _gene_class(args.selection, gpi_genes, by_acc, by_ref_gene)

    rows, unclassified = [], Counter()
    for r in csv.DictReader(open(args.cohorts), delimiter="\t"):
        cls = gene_class.get(r["gene"])
        if cls is None:
            unclassified[r["cohort_tier"]] += 1
            continue
        rows.append({**r, "topology_class": cls})

    # Genes whose own selected reference is secreted or intracellular cannot contribute to a
    # surface panel whatever their isoforms look like. Ten library genes are in this state
    # and they are written out rather than silently dropped: six are reported in the
    # literature as GPI-anchored but carry no GPI keyword and no lipidation feature in
    # UniProt, so classifying them here would be asserting something the annotation does not
    # say. Two others -- PILRA and LMLN -- reached it because the appris_only branch moved
    # their reference onto a soluble isoform.
    non_surface = sorted({
        r["gene"] for r in rows if r["topology_class"] in ("no TM (SP/GLOB)", "no TM helix")
    })
    if non_surface:
        Path("panel_excluded_non_surface_reference.txt").write_text(
            "\n".join(non_surface) + "\n"
        )
        print(
            f"excluded {len(non_surface)} gene(s) whose reference is not surface: "
            f"{', '.join(non_surface)}"
        )
    rows = [
        r for r in rows if r["topology_class"] not in ("no TM (SP/GLOB)", "no TM helix")
    ]

    print(f"classified {len(rows)} forms; no topology for {sum(unclassified.values())}")
    for tier, n in unclassified.most_common():
        print(f"  unclassified in {tier}: {n}")

    pools = {
        "observed": [r for r in rows if r["cohort_tier"] == "observed"],
        "prediction_only": [r for r in rows if r["cohort_tier"] == "full_only"],
    }
    print("\nclass distribution:")
    for tier, pool in pools.items():
        c = Counter(r["topology_class"] for r in pool)
        print(f"  {tier} ({len(pool)}): {dict(c.most_common())}")

    rng = random.Random(SEED)
    plan = {}
    strata = []
    for tier, want in (("observed", N_OBSERVED), ("prediction_only", N_PREDICTED)):
        pool = pools[tier]
        alloc = allocate(pool, want)
        plan[tier] = alloc
        by_cls = defaultdict(list)
        for r in pool:
            by_cls[r["topology_class"]].append(r)
        for cls, k in alloc.items():
            strata.append((tier, cls, by_cls[cls], rng.sample(by_cls[cls], k)))

    # Replacement runs after every original draw, so the valid picks never move.
    def form_id(r: dict) -> str:
        return f"{r['gene']}:{r['uniprot'].strip() or r['ensembl_transcripts']}"

    with open(args.selection) as fh:
        selected = {r["symbol"]: (r.get("selected_uniprot_equivalent") or "").strip()
                    for r in csv.DictReader(fh, delimiter="\t")}

    gate: dict[str, str] = {}
    if args.gate and not args.original:
        import hashlib

        if args.snapshot:
            pin = json.loads((Path(__file__).with_name("panel_inputs.json")).read_text())
            got = hashlib.sha256(Path(args.gate).read_bytes()).hexdigest()
            if got != pin["isoform_gate"]["sha256"]:
                raise SystemExit("--gate table does not match panel_inputs.json isoform_gate")
        with open(args.gate) as fh:
            gate = {g["form_id"]: g["reason"] for g in csv.DictReader(fh, delimiter="\t")
                    if g["verdict"] == "fail"}
            fh.seek(0)
            known = {g["form_id"] for g in csv.DictReader(fh, delimiter="\t")}
        absent = [form_id(r) for r in rows if form_id(r) not in known]
        if absent:
            raise SystemExit(f"{len(absent)} cohort forms have no gate verdict: {absent[:3]}")

    def invalid(r: dict) -> str:
        sel = selected.get(r["gene"], "")
        if sel and r["uniprot"] == sel:
            return "is the gene's selected reference form"
        why = gate.get(form_id(r))
        return f"fails the cohort gate: {why}" if why else ""

    rejected = {} if args.original else load_rejections(args.reject, "isoform")
    if args.original:
        def invalid(r: dict) -> str:
            return ""
    used = {form_id(p) for *_, picked in strata for p in picked}
    panel, replaced = [], []
    for tier, cls, stratum, picked in strata:
        final, log = replace_rejected(rng, picked, stratum, form_id, rejected, invalid, used)
        panel += [{**p, "panel_tier": tier} for p in final]
        replaced += [{"panel_tier": tier, "topology_class": cls, **row} for row in log]

    panel.sort(key=lambda r: (r["panel_tier"], r["topology_class"], r["gene"]))
    with open("isoform_panel.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(panel[0]), delimiter="\t")
        w.writeheader()
        w.writerows(panel)
    print(f"\nwrote isoform_panel.tsv ({len(panel)} forms)")
    with open("isoform_panel_replacements.tsv", "w", newline="") as fh:
        cols = ["panel_tier", "topology_class", "rejected", "reason", "replacement", "candidates"]
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(replaced)
    for row in replaced:
        print(f"  replaced {row['rejected']} ({row['reason']}) -> {row['replacement'] or 'NONE'}")

    summary = {
        "seed": SEED,
        "inputs": "panel_inputs.json snapshot" if args.snapshot else "local files",
        "replaced": len(replaced),
        "panel_size": len(panel),
        "floors": FLOORS,
        "allocation": plan,
        "by_class": dict(Counter(r["topology_class"] for r in panel).most_common()),
        "genes": len({r["gene"] for r in panel}),
        "already_in_cohort": sum(1 for r in panel if r["in_current_cohort"] == "True"),
    }
    Path("isoform_panel_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
