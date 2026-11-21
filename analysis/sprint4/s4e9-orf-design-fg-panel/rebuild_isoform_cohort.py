"""Rebuild the isoform cohort fresh over every library gene, on sequence-matched DeepTMHMM v1,
and draw a PROPOSED isoform panel from it. Nothing of record changes.

Becca, 2026-11-20, on the brief's "Isoforms are drawn from a pool narrower than the library":
"then we need to fix this fresh". Two findings shape what "fresh" means here.

1. **The isoform universe was never limited to the April query set.** ``isoform_universe.py``
   iterates the committed library list (2,369 genes at 88786350), not the April
   ``surfaceome_expressed`` set; the 1,775-of-2,369 figure belongs to the original ortholog pool
   (``check_universe_coverage.py``). The cohort of record already holds 333 forms of 183 genes
   outside the April set, and 26 of the panel's 133 genes are outside it. The one library gene
   the universe lacks is SLC67A1 (formerly SLC22A18), added to the library on 2026-11-02.
2. **What the April files did decide is topology.** The cohort's membrane tests read DeepTMHMM
   labels keyed on UniProt accession from the April accessible-surfaceome tables, falling back
   to a separate October run for forms those tables lacked, and the reference they compared
   against was the universe's MANE/APPRIS pick rather than the library's built reference.

So the rebuild keeps the universe's evidence rules exactly and changes two things: the gene
list is the current library (2,369 genes: the committed list plus SLC67A1), and every
topology decision is DeepTMHMM v1 on the form's own sequence, compared with the library's
built reference protein (the gate's reference).

Steps (``--work`` caches the slow first one):

* ``universe``: the universe, OpenProt and APPRIS steps rerun at their producing commits
  (``reproduce_isoform_universe.py``) on the current library list. PeptideAtlas peptides come
  from the pinned build-607 replay; an accession it lacks is fetched from build 607 live and
  the fetched peptides are written to ``isoform_cohort_rebuild_peptides.json``.
* ``fasta``: every evidence-bar form (observed, or TRIFID > 0.5) and every reference with no v1
  row, for ``isoform_v1_topology.py``'s Modal run and publication.
* ``build``: the fresh cohort, the isoform gate on it (``gate_isoform_cohort.gate``), the
  sampler (``sample_isoform_panel.py``, seed 20261009, its floors and its gene classes from the
  pinned panel inputs, ``panel_rejections.tsv``), and the comparison with ``isoform_panel.tsv``.
  The same sampler is also run on the cohort of record with the full-v1 gate
  (``isoform_gate_full_v1.tsv``), so the effect of the new topology rows alone is separable
  from the effect of the fresh cohort.

The cohort rules, as ``isoform_cohorts.py`` (6fb98222) applies them:

* evidence bar: an observed peptide tier, an OpenProt translation score above 0, or TRIFID
  above 0.5; tier ``observed`` if observed, else ``full_only``;
* lost anchoring: excluded when the reference is a surface label (TM, SP+TM, BETA) and the form
  is not, or the reference has a TM helix and the form none;
* soluble restriction: excluded when the form is not a surface label and does not carry its
  gene's GPI anchor (NetGPI 1.1 on the form where called, else the C-terminal retention test
  against the reference; GPI genes from ``gpi_set.tsv``, the set of record).

A form with no v1 row, or whose gene has no built reference, cannot be judged and is listed
rather than admitted.

    uv run --locked python rebuild_isoform_cohort.py universe --work W
    uv run --locked python rebuild_isoform_cohort.py fasta --work W --base-key K --base-sha256 S \\
        --out topology_runs/<run>.fasta
    uv run --locked python rebuild_isoform_cohort.py build --work W --snapshot-key K \\
        --snapshot-sha256 S --cohort-of-record-snapshot-key K1 --cohort-of-record-snapshot-sha256 S1

Decision: Becca Carlson, 2026-11-20 (review comment 28 on the L1 design brief); the proposed
panel is for review and is not adopted. Implemented by the Claude Code session
l1-isoform-universe (subagent of the brief review).
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import sys
import tempfile
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
AP = REPO / "scratch/becca/appris-principal"
sys.path.insert(0, str(AP))
sys.path.insert(0, str(HERE))

import gate_isoform_cohort as gate_mod  # noqa: E402
import reproduce_isoform_universe as rep  # noqa: E402

LIBRARY_FILE = REPO / "scratch/becca/tag-site-architecture/library_members.json"
PANEL_OF_RECORD = HERE / "isoform_panel.tsv"
REJECTIONS = HERE / "panel_rejections.tsv"
GATE_FULL_V1 = HERE / "isoform_gate_full_v1.tsv"
OUT_COHORT = HERE / "isoform_cohort_rebuild.tsv"
OUT_EXCLUDED = HERE / "isoform_cohort_rebuild_excluded.tsv"
OUT_GATE = HERE / "isoform_gate_rebuild.tsv"
OUT_PANEL = HERE / "isoform_panel_proposed.tsv"
OUT_PANEL_V1 = HERE / "isoform_panel_proposed_cohort_of_record.tsv"
OUT_CHANGES = HERE / "isoform_panel_proposed_changes.tsv"
OUT_GENES = HERE / "isoform_cohort_rebuild_genes.tsv"
OUT_SUMMARY = HERE / "isoform_cohort_rebuild_summary.json"
OUT_PEPTIDES = HERE / "isoform_cohort_rebuild_peptides.json"

OBSERVED_TIERS = {"expressed", "likely_expressed", "weak"}
SURFACE_LABELS = {"TM", "SP+TM", "BETA"}
TRIFID_CUT = 0.5
REFERENCE_FORM = "is the gene's built reference protein"
DECIDED_BY = "Becca Carlson, 2026-11-20 (review comment 28: 'then we need to fix this fresh')"
IMPLEMENTED_BY = "Claude Code session l1-isoform-universe (subagent of the brief review)"
COHORT_COLUMNS = [
    "gene",
    "uniprot",
    "ensembl_transcripts",
    "length",
    "reference_length",
    "source",
    "in_current_cohort",
    "cohort_tier",
    "evidence",
    "peptide_tier",
    "openprot_te_score",
    "trifid_score",
    "appris_label_class",
]


def sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def library_genes() -> list[str]:
    """The current library: the committed list plus active L1 reference genes it lacks."""
    from deliverome_analysis import library_design as ld

    listed = set(json.loads(LIBRARY_FILE.read_text()))
    built = {
        m["gene_symbol"]
        for m in ld.designed_proteins(
            HERE / "candidate_manifest_tranches.tsv", library=ld.LIBRARY_L1
        ).values()
        if m["candidate_set"] == "library_reference"
    }
    return sorted(listed | built)


# ---- step 1: the universe ----------------------------------------------------------------
def build_universe(work: Path) -> Path:
    """Universe, OpenProt and APPRIS steps at their producing commits, current library."""
    out = work / "run"
    if (out / "isoform_universe_scored.tsv").exists():
        return out
    got = rep.stage(work)
    cache, repo = rep.layout(work, got)
    out.mkdir(parents=True, exist_ok=True)
    lib = work / "library_members.json"
    lib.write_text(json.dumps(library_genes()))
    pools = cache / "peptideatlas_build607_pools.json"
    before = json.loads(pools.read_text())["index"].keys()
    code = rep.code_at(rep.COMMITS["universe"], work / "code" / "u")
    rep.run(
        code,
        "isoform_universe.py",
        [
            "--surfaceome-repo",
            str(repo),
            "--cache-dir",
            str(cache),
            "--library",
            str(lib),
            "--selection",
            str(got["reference_selection.tsv"]),
        ],
        out,
    )
    blob = json.loads(pools.read_text())
    new = sorted(set(blob["index"]) - set(before))
    fetched = {
        "atlas_build_id": 607,
        "fetched": "2026-11-20, PeptideAtlas build 607 live (accessions absent from the replay)",
        "index": {a: blob["index"][a] for a in new},
        "pools": {k: blob["pools"][k] for k in sorted({blob["index"][a] for a in new})},
    }
    OUT_PEPTIDES.write_text(json.dumps(fetched, indent=1, sort_keys=True) + "\n")
    rep.run(
        code,
        "openprot_evidence.py",
        ["--openprot-dir", str(work / "openprot"), "--cache-dir", str(cache)],
        out,
    )
    code = rep.code_at(rep.COMMITS["appris"], work / "code" / "a")
    rep.run(
        code,
        "isoform_appris.py",
        ["--cache-dir", str(cache), "--surfaceome-repo", str(repo)],
        out,
    )
    return out


def evidence_rows(work: Path) -> list[dict]:
    """Non-reference universe rows of library genes that clear the evidence bar."""
    genes = set(library_genes())
    with (work / "run" / "isoform_universe_scored.tsv").open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    out = []
    for r in rows:
        if r["is_reference"] != "False" or r["gene"] not in genes:
            continue
        te = int(r["openprot_te_score"] or 0)
        trifid = float(r["trifid_score"] or 0)
        observed = r["tier"] in OBSERVED_TIERS or te > 0
        if observed or trifid > TRIFID_CUT:
            out.append(
                {
                    **r,
                    "cohort_tier": "observed" if observed else "full_only",
                    "evidence": (
                        "protein (deep peptides)"
                        if r["tier"] in ("expressed", "likely_expressed")
                        else "protein (weak)"
                        if r["tier"] == "weak"
                        else "translation (OpenProt TE)"
                        if te > 0
                        else "prediction (TRIFID)"
                    ),
                    "peptide_tier": r["tier"],
                }
            )
    return out


def candidate_table(work: Path) -> Path:
    """The evidence-bar rows in cohort-table shape, for the gate's input loader."""
    path = work / "evidence_bar_forms.tsv"
    rows = evidence_rows(work)
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COHORT_COLUMNS, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    return path


# ---- step 2: the sequences that need a v1 run --------------------------------------------
def write_fasta(args) -> int:
    import isoform_v1_topology as v1

    build_universe(args.work)
    seqs = v1.cohort_sequences(str(candidate_table(args.work)), args.base_key, args.base_sha256)
    if args.out.exists():
        raise SystemExit(f"{args.out} exists; runs use fresh names")
    args.out.write_text("".join(f">{n}\n{s}\n" for n, s in seqs.items()))
    print(f"wrote {args.out} ({len(seqs)} sequences, {sum(map(len, seqs.values()))} residues)")
    return 0


# ---- step 3: cohort, gate, draw ----------------------------------------------------------
def v1_calls(inp: dict, seqs: dict[str, str]) -> dict[str, dict]:
    from deliverome_analysis import topology

    rep_ = topology.resolve(seqs, tool=topology.TOOL_V1, snapshot=inp["snapshot"])
    out = {}
    for row in rep_.to_dict("records"):
        if row["found"]:
            topo = row["per_residue_topology"]
            out[row["key"]] = {
                "label": row["deeptmhmm_label"],
                "tm": sum(
                    1 for j, c in enumerate(topo) if c == "M" and (j == 0 or topo[j - 1] != "M")
                ),
                "tool_version": row["tool_version"],
                "topology_version": row["topology_version"],
            }
    return out


def netgpi_calls() -> dict[tuple[str, str], str]:
    from netgpi_gpi_forms import classify, read

    out = {}
    for r in read(AP / "netgpi_gpi_human_predictions.txt"):
        _kind, rest = r["id"].split("_", 1)
        gene, tag = rest.rsplit("_", 1)
        out[(gene, tag)] = classify(r)
    return out


def fresh_cohort(work: Path, snapshot_key: str, snapshot_sha: str) -> tuple[list, list, dict]:
    from gpi import anchor_retained, gpi_genes

    inp = gate_mod.load_inputs(
        work / "gate-inputs",
        cohort_path=candidate_table(work),
        snapshot_key=snapshot_key,
        snapshot_sha256=snapshot_sha,
    )
    proteins = gate_mod.form_proteins(inp)
    refs = inp["refs"]
    calls = v1_calls(
        inp,
        {f"form:{k}": v for k, v in proteins.items() if v}
        | {f"ref:{g}": s for g, (_, s) in refs.items()},
    )
    gpi = gpi_genes(HERE / "gpi_set.tsv")
    netgpi = netgpi_calls()
    rows = {gate_mod.form_id(r): r for r in evidence_rows(work)}

    def inherits_anchor(r: dict, seq: str) -> tuple[bool, str]:
        if r["gene"] not in gpi:
            return False, ""
        for tag in [r["uniprot"], *(r["ensembl_transcripts"] or "").split(";")]:
            call = netgpi.get((r["gene"], tag.strip()))
            if call in ("anchored", "not_anchored"):
                return call == "anchored", f"NetGPI {call}"
        kept = anchor_retained(seq, refs[r["gene"]][1])
        return kept is not False, "C-terminal retention" if kept is not None else "unknown"

    cohort, excluded = [], []
    for fid, r in rows.items():
        seq = proteins.get(fid, "")
        form, ref = calls.get(f"form:{fid}"), calls.get(f"ref:{r['gene']}")
        base = {
            "form_id": fid,
            "gene": r["gene"],
            "cohort_tier": r["cohort_tier"],
            "evidence": r["evidence"],
            "form_label": form["label"] if form else "",
            "form_tm_helices": form["tm"] if form else "",
            "reference_label": ref["label"] if ref else "",
            "reference_tm_helices": ref["tm"] if ref else "",
        }
        why = ""
        if not seq:
            why = "no protein sequence for the form"
        elif r["gene"] not in refs:
            why = "unjudgeable: gene has no built L1 reference"
        elif seq == refs[r["gene"]][1]:
            # The universe's reference is its MANE/APPRIS pick; where an override built a
            # different form (NCAM1, CD200R1L, ...), that form is the reference, not an isoform.
            why = REFERENCE_FORM
        elif form is None or ref is None:
            why = "unjudgeable: no DeepTMHMM v1 row for the form or its reference"
        elif ref["label"] in SURFACE_LABELS and form["label"] not in SURFACE_LABELS:
            why = "loses surface anchoring: reference surface-labelled, form not"
        elif ref["tm"] > 0 and form["tm"] == 0:
            why = "loses surface anchoring: reference has a TM helix, form none"
        if not why and form and form["label"] not in SURFACE_LABELS:
            anchored, basis = inherits_anchor(r, seq)
            base["gpi_anchor_basis"] = basis
            if not anchored:
                why = (
                    "soluble: GPI gene but the form does not keep the anchor"
                    if r["gene"] in gpi
                    else "soluble: no TM helix and no GPI anchor"
                )
        if why == REFERENCE_FORM:
            # Kept in the pool and rejected at the draw, as the sampler treats its own
            # reference forms, so pool sizes and the seeded picks do not move for it.
            excluded.append({**base, "reason": f"{why} (kept in the pool, rejected at the draw)"})
            cohort.append({k: r.get(k, "") for k in COHORT_COLUMNS})
        elif why:
            excluded.append({**base, "reason": why})
        else:
            cohort.append({k: r.get(k, "") for k in COHORT_COLUMNS})
    return cohort, excluded, inp


def rejections_with_references(workdir: Path, excluded: list[dict]) -> Path:
    """``panel_rejections.tsv`` plus one isoform row per form that is its gene's built reference."""
    with REJECTIONS.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        cols, rows = list(reader.fieldnames or []), list(reader)
    known = {(r["panel"], r["form_id"]) for r in rows}
    for e in excluded:
        if e["reason"].startswith(REFERENCE_FORM) and ("isoform", e["form_id"]) not in known:
            rows.append(
                dict.fromkeys(cols, "")
                | {
                    "panel": "isoform",
                    "form_id": e["form_id"],
                    "reason": REFERENCE_FORM,
                    "decided_by": "rule (sample_isoform_panel.py: a gene's reference is not an isoform)",
                    "implemented_by": IMPLEMENTED_BY,
                }
            )
    path = workdir / "panel_rejections_with_built_references.tsv"
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    return path


def run_sampler(
    workdir: Path, cohorts: Path, gate: Path, rejections: Path = REJECTIONS
) -> tuple[list[dict], list[dict], dict]:
    """The panel sampler with the pinned panel inputs, on ``cohorts`` and ``gate``."""
    import sample_isoform_panel as s
    from gpi import gpi_genes
    from panel_inputs import stage

    snap = stage(workdir / "panel-inputs")
    fake = workdir / "surfaceome"
    (fake / "data/processed/deeptmhmm").mkdir(parents=True, exist_ok=True)
    for name in ("deeptmhmm_human_canonical.tsv", "deeptmhmm_human_isoforms.tsv"):
        dest = fake / "data/processed/deeptmhmm" / name
        if not dest.exists():
            shutil.copy(snap / name, dest)
    pinned_gpi = gpi_genes(snap / "gpi_set.tsv")
    s.gpi_gene_set = lambda *a, **k: pinned_gpi  # the panel of record's GPI set
    run = workdir / f"draw-{cohorts.stem}-{gate.stem}"
    run.mkdir(parents=True, exist_ok=True)
    argv, cwd = sys.argv, os.getcwd()
    sys.argv = [
        "sample_isoform_panel.py",
        "--surfaceome-repo",
        str(fake),
        "--cohorts",
        str(cohorts),
        "--selection",
        str(snap / "reference_selection.tsv"),
        "--reference-three-line",
        str(snap / "library_reference_missing_topology.3line"),
        "--reject",
        str(rejections),
        "--gate",
        str(gate),
    ]
    try:
        os.chdir(run)
        s.main()
    finally:
        sys.argv = argv
        os.chdir(cwd)
    read = lambda p: list(csv.DictReader((run / p).open(), delimiter="\t"))  # noqa: E731
    return (
        read("isoform_panel.tsv"),
        read("isoform_panel_replacements.tsv"),
        json.loads((run / "isoform_panel_summary.json").read_text()),
    )


def write_tsv(path: Path, rows: list[dict], cols: list[str] | None = None) -> None:
    cols = cols or (list(rows[0]) if rows else [])
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, "") for c in cols})


def build(args) -> int:
    work = args.work
    build_universe(work)
    cohort, excluded, inp = fresh_cohort(work, args.snapshot_key, args.snapshot_sha256)
    # The seeded sampler's picks depend on row order, and the record's order follows a hash
    # seed that was never recorded. Forms already in the cohort of record keep its order; new
    # forms follow in universe order, so only a real change in membership moves the draw.
    from panel_inputs import stage as stage_panel

    with (stage_panel(work / "panel-of-record") / "isoform_cohorts.tsv").open() as fh:
        order = {gate_mod.form_id(r): i for i, r in enumerate(csv.DictReader(fh, delimiter="\t"))}
    cohort.sort(key=lambda r: order.get(gate_mod.form_id(r), len(order)))
    write_tsv(OUT_COHORT, cohort, COHORT_COLUMNS)
    write_tsv(
        OUT_EXCLUDED,
        sorted(excluded, key=lambda r: (r["gene"], r["form_id"])),
        [
            "form_id",
            "gene",
            "cohort_tier",
            "evidence",
            "reason",
            "form_label",
            "form_tm_helices",
            "reference_label",
            "reference_tm_helices",
            "gpi_anchor_basis",
        ],
    )
    # The gate, unchanged, on the fresh cohort.
    inp = gate_mod.load_inputs(
        work / "gate-fresh",
        cohort_path=OUT_COHORT,
        snapshot_key=args.snapshot_key,
        snapshot_sha256=args.snapshot_sha256,
    )
    gated = gate_mod.gate(inp)
    for g in gated:
        g |= {
            "decided_by": DECIDED_BY,
            "implemented_by": IMPLEMENTED_BY,
            "decision_reason": gate_mod.DECISION_REASON,
        }
    OUT_GATE.write_text(gate_mod.tsv_text(gated))

    # Reproduce the panel of record first, so the wrapper is proven before it proposes.
    from panel_inputs import stage

    draw_dir = Path(tempfile.mkdtemp(prefix="isoform-draw-"))
    snap = stage(draw_dir / "panel-inputs")
    rec_panel, _, _ = run_sampler(draw_dir, snap / "isoform_cohorts.tsv", HERE / "isoform_gate.tsv")
    with PANEL_OF_RECORD.open() as fh:
        of_record = list(csv.DictReader(fh, delimiter="\t"))
    assert rec_panel == of_record, "the sampler wrapper does not reproduce isoform_panel.tsv"

    v1_panel, v1_repl, v1_sum = run_sampler(
        draw_dir,
        snap / "isoform_cohorts.tsv",
        GATE_FULL_V1,
        rejections_with_references(draw_dir, excluded),
    )
    new_panel, new_repl, new_sum = run_sampler(
        draw_dir, OUT_COHORT, OUT_GATE, rejections_with_references(draw_dir, excluded)
    )
    write_tsv(OUT_PANEL_V1, v1_panel)
    write_tsv(OUT_PANEL, new_panel)

    fid = gate_mod.form_id
    rec_ids = {fid(r): r for r in of_record}
    new_ids = {fid(r): r for r in new_panel}
    v1_ids = {fid(r): r for r in v1_panel}
    gate_new = {g["form_id"]: g for g in gated}
    cohort_ids = {fid(r) for r in cohort}
    excl = {e["form_id"]: e["reason"] for e in excluded}
    with (HERE / "isoform_panel_replacements.tsv").open() as fh:
        record_replacements = {r["replacement"] for r in csv.DictReader(fh, delimiter="\t")}
    changes = []
    for f in sorted(set(rec_ids) | set(new_ids)):
        r = rec_ids.get(f) or new_ids[f]
        status = "kept" if f in rec_ids and f in new_ids else "dropped" if f in rec_ids else "added"
        why = ""
        if status == "dropped":
            why = (
                excl.get(f)
                or (
                    f"fails the gate: {gate_new[f]['reason']}"
                    if f in gate_new and gate_new[f]["verdict"] == "fail"
                    else ""
                )
                or ("not in the rebuilt cohort" if f not in cohort_ids else "")
                or (
                    "a replacement pick of record; the seeded replacement draw now runs over "
                    "the larger gate-passing pool and fills the slot with another form"
                    if f in record_replacements
                    else "not drawn"
                )
            )
        changes.append(
            {
                "form_id": f,
                "gene": r["gene"],
                "status": status,
                "panel_of_record_stratum": (
                    f"{rec_ids[f]['panel_tier']} / {rec_ids[f]['topology_class']}"
                    if f in rec_ids
                    else ""
                ),
                "proposed_stratum": (
                    f"{new_ids[f]['panel_tier']} / {new_ids[f]['topology_class']}"
                    if f in new_ids
                    else ""
                ),
                "in_cohort_of_record_draw_with_full_v1_gate": f in v1_ids,
                "reason_if_dropped": why,
            }
        )
    write_tsv(OUT_CHANGES, changes)

    # Per-gene reach: which library genes can contribute a drawable isoform, before and after.
    with (draw_dir / "panel-inputs/isoform_cohorts.tsv").open() as fh:
        old_cohort = list(csv.DictReader(fh, delimiter="\t"))
    with (HERE / "isoform_gate.tsv").open() as fh:
        old_gate = {g["form_id"]: g["verdict"] for g in csv.DictReader(fh, delimiter="\t")}
    with GATE_FULL_V1.open() as fh:
        v1_gate = {g["form_id"]: g["verdict"] for g in csv.DictReader(fh, delimiter="\t")}
    genes = library_genes()
    per_gene = []
    for g in genes:
        oc = [fid(r) for r in old_cohort if r["gene"] == g]
        nc = [fid(r) for r in cohort if r["gene"] == g]
        per_gene.append(
            {
                "gene": g,
                "cohort_of_record_forms": len(oc),
                "cohort_of_record_drawable_gate_of_record": sum(
                    old_gate.get(f) == "pass" for f in oc
                ),
                "cohort_of_record_drawable_full_v1": sum(v1_gate.get(f) == "pass" for f in oc),
                "rebuilt_cohort_forms": len(nc),
                "rebuilt_drawable": sum(
                    gate_new[f]["verdict"] == "pass" for f in nc if f in gate_new
                ),
                "panel_of_record_forms": sum(1 for f in rec_ids if rec_ids[f]["gene"] == g),
                "proposed_panel_forms": sum(1 for f in new_ids if new_ids[f]["gene"] == g),
            }
        )
    write_tsv(OUT_GENES, per_gene)

    def gsum(rows: list[dict]) -> dict:
        return {"forms": len(rows), "genes": len({r["gene"] for r in rows})}

    reach = lambda k: sum(1 for r in per_gene if r[k] > 0)  # noqa: E731
    summary = {
        "decision": DECIDED_BY,
        "implemented_by": IMPLEMENTED_BY,
        "status": "proposed; the library, panel tables of record, constructs and order unchanged",
        "library_genes": len(genes),
        "topology_snapshot": {"key": args.snapshot_key, "sha256": args.snapshot_sha256},
        "topology_tool": "deeptmhmm_v1",
        "evidence_bar_forms": gsum(evidence_rows(work)),
        "rebuilt_cohort": gsum(cohort)
        | {"by_tier": dict(Counter(r["cohort_tier"] for r in cohort))},
        "excluded": dict(Counter(e["reason"] for e in excluded).most_common()),
        "cohort_of_record": gsum(old_cohort),
        "rebuilt_vs_record_forms": {
            "in_both": len(cohort_ids & {fid(r) for r in old_cohort}),
            "only_rebuilt": len(cohort_ids - {fid(r) for r in old_cohort}),
            "only_record": len({fid(r) for r in old_cohort} - cohort_ids),
        },
        "gate_rebuilt": {
            "pass": sum(g["verdict"] == "pass" for g in gated),
            "fail": sum(g["verdict"] == "fail" for g in gated),
            "unverifiable": sum("unverifiable" in g["reason"] for g in gated),
        },
        "genes_with_a_drawable_form": {
            "cohort_of_record_gate_of_record": reach("cohort_of_record_drawable_gate_of_record"),
            "cohort_of_record_full_v1_gate": reach("cohort_of_record_drawable_full_v1"),
            "rebuilt": reach("rebuilt_drawable"),
        },
        "panel_of_record": gsum(of_record),
        "proposed_cohort_of_record_full_v1_gate": {
            **gsum(v1_panel),
            "by_class": v1_sum["by_class"],
            "replaced": len(v1_repl),
            "same_forms_as_record": len(set(v1_ids) & set(rec_ids)),
        },
        "proposed_rebuilt": {
            **gsum(new_panel),
            "by_class": new_sum["by_class"],
            "allocation": new_sum["allocation"],
            "replaced": len(new_repl),
            "status": dict(Counter(c["status"] for c in changes)),
            "genes_added": sorted(
                {c["gene"] for c in changes if c["status"] == "added"}
                - {c["gene"] for c in changes if c["status"] != "added"}
            ),
            "genes_dropped": sorted(
                {c["gene"] for c in changes if c["status"] == "dropped"}
                - {c["gene"] for c in changes if c["status"] != "dropped"}
            ),
        },
        "outputs_sha256": {
            p.name: sha256(p)
            for p in (
                OUT_COHORT,
                OUT_EXCLUDED,
                OUT_GATE,
                OUT_PANEL,
                OUT_PANEL_V1,
                OUT_CHANGES,
                OUT_GENES,
            )
        },
    }
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({k: v for k, v in summary.items() if k != "outputs_sha256"}, indent=2))
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)
    u = sub.add_parser("universe")
    f = sub.add_parser("fasta")
    f.add_argument("--base-key", required=True)
    f.add_argument("--base-sha256", required=True)
    f.add_argument("--out", type=Path, required=True)
    b = sub.add_parser("build")
    b.add_argument("--snapshot-key", required=True)
    b.add_argument("--snapshot-sha256", required=True)
    for p in (u, f, b):
        p.add_argument("--work", type=Path, required=True)
    args = ap.parse_args(argv)
    if args.cmd == "universe":
        print(build_universe(args.work))
        return 0
    return write_fasta(args) if args.cmd == "fasta" else build(args)


if __name__ == "__main__":
    raise SystemExit(main())
