"""Proposed isoform-panel redraw of 2026-11-22: rebuilt cohort, multipass split by N-side.

**Decided by Becca Carlson, 2026-11-22**: "yes we should always rebuild with the latest";
"yes let's split this up into two" (the multipass stratum); "shouldn't some features be ok
to use as long as they are compatible with display? we are fine with diversity" (the feature
rule). Implemented by Claude Code session l1-review-export-scope (worker fork), 2026-11-22.

What changes against the proposal of 2026-11-20 (``isoform_panel_proposed.tsv``):

* **Strata.** The multipass class is split by the N-terminal side of the gene's selected
  reference (DeepTMHMM v1, read exactly as ``sample_isoform_panel.classify`` reads it):
  ``multipass N-out`` (extracellular) and ``multipass N-in`` (cytoplasmic). Every other
  class, the floors (type III 2, GPI 5, type II 5), the tier split (100 observed, 50
  prediction-only), the seed (20261009) and the largest-remainder allocation are unchanged.
* **Feature rule.** A form is replaced only when it loses a display-deciding feature: the
  signal peptide, the N-terminal side, the first TM helix or the membrane anchor. The isoform
  cohort gate (``gate_isoform_cohort.py``) already tests exactly these (plus the native ATG
  start and an exact Ensembl 115 CDS, which decide whether the form can be built), and the
  rebuilt cohort already excludes forms that lose the membrane or GPI anchor. So the narrower
  rule changes no isoform verdict; it is recorded here so the two panels use one rule.

The draw runs the committed sampler through ``rebuild_isoform_cohort.run_sampler`` (pinned
panel inputs, pinned GPI set) on ``isoform_cohort_rebuild.tsv`` and ``isoform_gate_rebuild.tsv``.
As a guard it first redraws without the split and must reproduce ``isoform_panel_proposed.tsv``.

The diff is against the panel of record (``isoform_panel.tsv``). Kept and dropped picks carry
their design id from the manifest; an added pick carries its existing id when it is already a
manifest member (an L2 member, or a retired L1 row) and otherwise a *proposed* next-free L1 id,
numbered from ``member_id_allocations.next_l1_member_number``. Nothing is allocated: the
manifest, references, CDS records and rung tables are not touched. Each added form is looked
up in the shared topology snapshot by exact sequence (``topology.resolve``, DeepTMHMM v1).

    uv run --locked python redraw_isoform_panel_2026_11_22.py           # write the outputs
    uv run --locked python redraw_isoform_panel_2026_11_22.py --check   # replay, write nothing
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import random
import sys
import tempfile
import types
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import rebuild_isoform_cohort  # noqa: E402, F401  (puts the committed sampler on sys.path)

OUT = HERE / "isoform_panel_redraw_2026-11-22.tsv"
OUT_SUMMARY = HERE / "isoform_panel_redraw_2026-11-22_summary.json"
MANIFEST = HERE / "candidate_manifest_tranches.tsv"
#: The 2026-11-12 panel of record, kept under a dated name. ``isoform_panel.tsv`` now holds
#: the applied 2026-11-22 redraw plus refills (``write_panels_of_record.py``); this redraw
#: is a diff against the 2026-11-12 draw, so it reads the dated copy.
PANEL_OF_RECORD = HERE / "isoform_panel_2026-11-12.tsv"
PROPOSED_2026_11_20 = HERE / "isoform_panel_proposed.tsv"
DECIDED_BY = (
    "Becca Carlson, 2026-11-22 ('yes we should always rebuild with the latest'; "
    "'yes let's split this up into two'; narrower feature rule)"
)
IMPLEMENTED_BY = "Claude Code session l1-review-export-scope (worker fork), 2026-11-22"
N_OUT, N_IN = "multipass N-out", "multipass N-in"
TIER_NAME = {"observed": "observed", "full_only": "prediction_only"}
COLUMNS = [
    "status",
    "form_id",
    "gene",
    "uniprot",
    "ensembl_transcripts",
    "panel_tier",
    "topology_class",
    "panel_of_record_stratum",
    "design_id",
    "design_id_status",
    "manifest_library",
    "manifest_retired_reason",
    "reason",
    "length",
    "v1_topology_in_shared_snapshot",
    "v1_topology_version",
    "v1_topology_in_gate_snapshot",
    "v1_topology_version_gate_snapshot",
    "decided_by",
    "implemented_by",
]


def split_classify(base):
    """``classify`` with the multipass class split by the reference's N-terminal side."""

    def classify(tm: int, sp: int, n_term: str, is_gpi: bool) -> str:
        cls = base(tm, sp, n_term, is_gpi)
        if cls != "multipass":
            return cls
        if n_term == "extracellular":
            return N_OUT
        if n_term == "cytoplasmic":
            return N_IN
        raise ValueError(f"multipass reference with N-terminus {n_term!r}")

    return classify


class PerStratumRandom:
    """A ``random.Random`` stand-in giving every stratum its own seeded stream.

    The sampler draws each stratum with ``rng.sample(stratum, k)`` and each replacement with
    ``rng.choice(candidates)``, where every element shares one (cohort tier, class). The call
    is routed to ``random.Random("<seed>|<tier>|<class>")`` (a string seed, hashed with
    SHA-512 by Python's version-2 seeding, so it is reproducible across runs and platforms).
    A change in one stratum then leaves every other stratum's draw untouched.
    """

    def __init__(self, seed: int):
        self.seed = seed
        self.streams: dict[str, random.Random] = {}

    def _for(self, population) -> random.Random:
        items = list(population)
        names = {f"{TIER_NAME[r['cohort_tier']]}|{r['topology_class']}" for r in items}
        assert len(names) == 1, f"a draw spans strata: {sorted(names)}"
        name = names.pop()
        if name not in self.streams:
            self.streams[name] = random.Random(f"{self.seed}|{name}")
        return self.streams[name]

    def sample(self, population, k):
        return self._for(population).sample(population, k)

    def choice(self, seq):
        return self._for(seq).choice(seq)


def draw(
    work: Path, split: bool, extra: dict[str, str] | None = None, per_stratum: bool = False
) -> tuple[list[dict], list[dict], dict]:
    """The committed sampler on the rebuilt cohort and gate; ``extra`` adds rejections."""
    import rebuild_isoform_cohort as rb
    import sample_isoform_panel as s
    from redraw_panels import REJECTION_COLUMNS, RULE_ROW

    with rb.OUT_EXCLUDED.open() as fh:
        excluded = list(csv.DictReader(fh, delimiter="\t"))
    base, base_random = s.classify, s.random
    if split:
        s.classify = split_classify(base)
    if per_stratum:
        s.random = types.SimpleNamespace(Random=PerStratumRandom)
    try:
        sub = work / (("split" if split else "unsplit") + ("-per-stratum" if per_stratum else ""))
        sub.mkdir(parents=True, exist_ok=True)
        rej = rb.rejections_with_references(sub, excluded)
        if extra:
            with rej.open("a", newline="") as fh:
                w = csv.DictWriter(
                    fh, fieldnames=REJECTION_COLUMNS, delimiter="\t", lineterminator="\n"
                )
                for f, why in sorted(extra.items()):
                    w.writerow(RULE_ROW | {"panel": "isoform", "form_id": f, "reason": why})
        return rb.run_sampler(sub, rb.OUT_COHORT, rb.OUT_GATE, rej)
    finally:
        s.classify, s.random = base, base_random


def draw_buildable(
    work: Path, manifest: list[dict], per_stratum: bool = False
) -> tuple[list[dict], list[dict], dict, dict]:
    """The split draw, iterated as ``redraw_panels.py`` does until every pick is buildable.

    A pick is rejected (and replaced within its stratum by the seeded draw) when no designed
    CDS encodes it, or when its designed protein repeats an active L1 member outside the
    isoform panel or an earlier pick. Returns the panel, its replacements, the sampler's
    summary and the rejections added here.
    """
    import json as _json

    from redraw_panels import iso_keys

    from deliverome_analysis import library_design as ld
    from deliverome_analysis import seqtools
    from deliverome_analysis.s3 import processed

    with open(processed.download(ld.CDS_KEY)) as fh:
        cds = _json.load(fh)
    keys = {m["key"] for m in manifest}
    protein = lambda k: seqtools.translate(cds[k]["cds"]).rstrip("*")  # noqa: E731
    fixed = {
        protein(m["key"]): m["design_id"]
        for m in manifest
        if m["library"] == "l1"
        and not m["retired_reason"]
        and m["set"] != "isoform"
        and m["key"] in cds
    }
    extra: dict[str, str] = {}
    for _ in range(10):
        panel, repl, summary = draw(work, split=True, extra=extra, per_stratum=per_stratum)
        seen, new = dict(fixed), {}
        for r in panel:
            f = f"{r['gene']}:{r['uniprot'].strip() or r['ensembl_transcripts']}"
            hit = sorted(iso_keys(r) & keys)
            if not hit or hit[0] not in cds:
                new[f] = "no designed CDS (no transcript encodes the drawn form)"
                continue
            p = protein(hit[0])
            if p in seen:
                new[f] = f"designed protein duplicates {seen[p]}"
                continue
            seen[p] = hit[0]
        new = {k: v for k, v in new.items() if k not in extra}
        if not new:
            return panel, repl, summary, extra
        extra |= new
    raise SystemExit("no fixed point after 10 passes")


def gene_classes(work: Path) -> Counter:
    """Library genes per (split) class, from the same pinned inputs the draw reads."""
    import sample_isoform_panel as s
    from gpi import gpi_genes
    from panel_inputs import stage

    snap = stage(work / "split" / "panel-inputs")
    base = s.classify
    s.classify = split_classify(base)
    try:
        by_acc, by_ref = s._topology_sources(
            reference_three_line=snap / "library_reference_missing_topology.3line",
            deeptmhmm_dir=snap,
        )
        classes = s._gene_class(
            snap / "reference_selection.tsv", gpi_genes(snap / "gpi_set.tsv"), by_acc, by_ref
        )
    finally:
        s.classify = base
    return Counter(classes.values())


def rebuild_snapshot(work: Path, key: str, sha: str):
    """The topology snapshot the rebuilt gate read (pinned in its summary), checksummed."""
    import hashlib

    import pandas as pd

    from deliverome_analysis.s3 import external

    path = Path(external.download(key, work / "gate-topology.parquet"))
    if hashlib.sha256(path.read_bytes()).hexdigest() != sha:
        raise SystemExit(f"{key} differs from its pinned SHA-256")
    return pd.read_parquet(path)


def read(path: Path) -> list[dict]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def sequences(forms: list[dict]) -> dict[str, str]:
    """form id -> protein, exactly as the gate reads it (UniProt isoform, else Ensembl 115)."""
    import candidate_inputs
    import gate_isoform_cohort as gate_mod
    import resolve_candidate_cds as rcd

    pep_of: dict[str, str] = {}
    for _ensp, (head, seq) in rcd.read_fasta(candidate_inputs.pep_fasta("human")).items():
        tx = next((f.split(":", 1)[1] for f in head.split() if f.startswith("transcript:")), "")
        if tx:
            pep_of[tx.split(".")[0]] = seq.rstrip("*").upper()
    inp = {"cohort": forms, "uniprot": rcd.load_uniprot("human"), "pep_of": pep_of}
    return gate_mod.form_proteins(inp)


def _redraw_summary(
    rows: list[dict],
    summary: dict,
    work,
    gate: dict,
    replacements: list,
    per_stratum: bool,
    buildability,
) -> dict:
    """The redraw's summary: the draw, and its difference from the panel of record."""
    count = Counter(r["status"] for r in rows)
    by_stratum: dict[str, Counter] = {}
    for r in rows:
        stratum = (
            r["panel_of_record_stratum"]
            if r["status"] == "dropped"
            else f"{r['panel_tier']} / {r['topology_class']}"
        )
        by_stratum.setdefault(stratum, Counter())[r["status"]] += 1
    added_rows = [r for r in rows if r["status"] == "added"]
    out = {
        "decided_by": DECIDED_BY,
        "implemented_by": IMPLEMENTED_BY,
        "status": "proposed; manifest, references, CDS records and rung tables unchanged",
        "seed": summary["seed"],
        "panel_size": summary["panel_size"],
        "genes": summary["genes"],
        "floors": summary["floors"],
        "allocation": summary["allocation"],
        "by_class": summary["by_class"],
        "library_genes_by_class": dict(gene_classes(work).most_common()),
        "cohort_forms": len(read(HERE / "isoform_cohort_rebuild.tsv")),
        "gate_pass": sum(g["verdict"] == "pass" for g in gate.values()),
        "replaced_in_draw": len(replacements),
        "random_stream": (
            "one per stratum: random.Random('20261009|<tier>|<class>')"
            if per_stratum
            else "one shared stream: random.Random(20261009), strata in sequence"
        ),
        "diff_vs_panel_of_record": {k: count[k] for k in ("kept", "dropped", "added")},
        "diff_by_stratum": {
            k: {s: v[s] for s in ("kept", "dropped", "added")}
            for k, v in sorted(by_stratum.items())
        },
        "rejected_as_unbuildable_in_redraw": buildability,
        "added_moved_from_l2": sum("moved from L2" in r["design_id_status"] for r in added_rows),
        "added_with_earlier_retired_l1_id": sum(
            "stays retired" in r["design_id_status"] for r in added_rows
        ),
        "added_without_l2_member": sorted(
            r["form_id"] for r in added_rows if "no L2 member" in r["design_id_status"]
        ),
        "new_l1_ids": (
            f"{min(r['design_id'] for r in added_rows)}..{max(r['design_id'] for r in added_rows)}"
            if added_rows
            else ""
        ),
        "added_without_v1_topology_in_shared_snapshot": sorted(
            r["form_id"] for r in added_rows if r["v1_topology_in_shared_snapshot"] != "yes"
        ),
        "added_without_v1_topology_in_gate_snapshot": sorted(
            r["form_id"] for r in added_rows if r["v1_topology_in_gate_snapshot"] != "yes"
        ),
        "would_retire": sorted(r["design_id"] for r in rows if r["status"] == "dropped"),
    }
    return out


def build(work: Path, per_stratum: bool = False) -> tuple[list[dict], dict]:
    import gate_isoform_cohort as gate_mod
    from member_id_allocations import next_l1_member_number
    from redraw_panels import iso_keys

    from deliverome_analysis import topology

    fid = gate_mod.form_id
    unsplit, _, _ = draw(work, split=False)
    assert unsplit == read(PROPOSED_2026_11_20), (
        "the wrapper does not reproduce isoform_panel_proposed.tsv"
    )
    manifest = read(MANIFEST)
    panel, replacements, summary, buildability = draw_buildable(work, manifest, per_stratum)

    record = read(PANEL_OF_RECORD)
    by_key: dict[str, list[dict]] = {}
    for m in manifest:
        if m["set"] == "isoform":
            by_key.setdefault(m["key"], []).append(m)

    def member(r: dict) -> dict | None:
        """The manifest row for a form: the active L1 row first, then L2, then retired."""
        rows = [m for k in iso_keys(r) for m in by_key.get(k, [])]
        rows.sort(key=lambda m: (bool(m["retired_reason"]), m["library"] != "l1", m["design_id"]))
        return rows[0] if rows else None

    rec = {fid(r): r for r in record}
    new = {fid(r): r for r in panel}
    gate = {g["form_id"]: g for g in read(HERE / "isoform_gate_rebuild.tsv")}
    replaced_by = {r["replacement"]: r for r in replacements}
    next_free = next_l1_member_number(manifest)
    added = sorted(set(new) - set(rec))
    seqs = sequences([new[f] for f in added])
    topo = topology.resolve({f: seqs[f] for f in added if seqs.get(f)}, tool=topology.TOOL_V1)
    topo_of = {row["key"]: row for row in topo.to_dict("records")}
    pin = json.loads((HERE / "isoform_cohort_rebuild_summary.json").read_text())[
        "topology_snapshot"
    ]
    gate_topo = topology.resolve(
        {f: seqs[f] for f in added if seqs.get(f)},
        tool=topology.TOOL_V1,
        snapshot=rebuild_snapshot(work, pin["key"], pin["sha256"]),
    )
    gate_topo_of = {row["key"]: row for row in gate_topo.to_dict("records")}

    rows: list[dict] = []
    for f in sorted(set(rec) | set(new), key=lambda f: (f not in new, f not in rec, f)):
        r = new.get(f) or rec[f]
        status = "kept" if f in rec and f in new else "dropped" if f in rec else "added"
        m = member(r)
        row = {
            "status": status,
            "form_id": f,
            "gene": r["gene"],
            "uniprot": r["uniprot"],
            "ensembl_transcripts": r["ensembl_transcripts"],
            "panel_tier": new[f]["panel_tier"] if f in new else "",
            "topology_class": new[f]["topology_class"] if f in new else "",
            "panel_of_record_stratum": (
                f"{rec[f]['panel_tier']} / {rec[f]['topology_class']}" if f in rec else ""
            ),
            "manifest_library": m["library"] if m else "",
            "manifest_retired_reason": m["retired_reason"] if m else "",
            "decided_by": DECIDED_BY,
            "implemented_by": IMPLEMENTED_BY,
        }
        if status != "added":
            assert m and m["library"] == "l1", f"{f}: panel-of-record pick has no L1 row"
            row["design_id"] = m["design_id"]
            row["design_id_status"] = (
                "existing (unchanged)"
                if status == "kept"
                else "existing; would be retired dropped_from_panel"
            )
            if status == "dropped":
                g = gate.get(f)
                row["reason"] = (
                    f"fails the rebuilt gate: {g['reason']}"
                    if g and g["verdict"] == "fail"
                    else "not in the rebuilt cohort"
                    if g is None
                    else "not drawn: the split strata and the rebuilt pool change the seeded draw"
                )
        else:
            # redraw_panels.py gives every pick that is not an active L1 member a new id,
            # moved from its L2 row; an earlier, retired L1 id is never reused.
            l2 = next(
                (x for k in iso_keys(r) for x in by_key.get(k, []) if x["library"] == "l2"), None
            )
            l1 = [x for k in iso_keys(r) for x in by_key.get(k, []) if x["library"] == "l1"]
            now = next((x for x in l1 if not x["retired_reason"]), None)
            l1_old = next((x for x in l1 if x["retired_reason"]), None)
            if now:
                row["design_id"] = now["design_id"]
                status = "allocated by apply_panel_redraw_2026_11_22.py"
            else:
                row["design_id"] = f"l1-p{next_free:04d}"
                next_free += 1
                status = "proposed next free L1 id (not allocated)"
            row["design_id_status"] = (
                status
                + (
                    f"; moved from L2 member {l2['design_id']}"
                    if l2
                    else "; no L2 member (blocks the redraw)"
                )
                + (
                    f"; earlier L1 id {l1_old['design_id']} ({l1_old['retired_reason']}) stays retired"
                    if l1_old
                    else ""
                )
            )
            repl = replaced_by.get(f)
            row["reason"] = (
                f"replacement for {repl['rejected']} ({repl['reason']})" if repl else "drawn"
            )
            seq = seqs.get(f, "")
            row["length"] = len(seq)
            t = topo_of.get(f)
            row["v1_topology_in_shared_snapshot"] = (
                "no sequence"
                if not seq
                else "yes"
                if t and t["found"]
                else "no (not in the shared snapshot)"
            )
            row["v1_topology_version"] = (t or {}).get("topology_version", "") if t else ""
            gt = gate_topo_of.get(f)
            row["v1_topology_in_gate_snapshot"] = "yes" if gt and gt["found"] else "no"
            row["v1_topology_version_gate_snapshot"] = (
                gt.get("topology_version", "") if gt and gt["found"] else ""
            )
        rows.append(row)

    out = _redraw_summary(rows, summary, work, gate, replacements, per_stratum, buildability)
    return rows, out


def tsv_text(rows: list[dict]) -> str:
    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
    w.writeheader()
    for r in rows:
        w.writerow({c: r.get(c, "") for c in COLUMNS})
    return buf.getvalue()


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--check", action="store_true", help="replay; fail if the outputs differ")
    ap.add_argument("--work", type=Path, default=None)
    ap.add_argument(
        "--per-stratum-seed",
        action="store_true",
        help="scenario: each stratum draws from its own seeded stream; outputs *_per_stratum",
    )
    a = ap.parse_args(argv)
    work = a.work or Path(tempfile.mkdtemp(prefix="isoform-redraw-"))
    rows, summary = build(work, a.per_stratum_seed)
    out, out_summary = OUT, OUT_SUMMARY
    if a.per_stratum_seed:
        out = OUT.with_name(OUT.stem + "_per_stratum.tsv")
        out_summary = OUT_SUMMARY.with_name(
            OUT_SUMMARY.name.replace("_summary.json", "_per_stratum_summary.json")
        )
    text, js = tsv_text(rows), json.dumps(summary, indent=2) + "\n"
    print(json.dumps({k: v for k, v in summary.items() if k != "would_retire"}, indent=2))
    if a.check:
        same = out.exists() and out.read_text() == text and out_summary.read_text() == js
        print("matches the record" if same else "DIFFERS from the record")
        return 0 if same else 1
    out.write_text(text)
    out_summary.write_text(js)
    print(f"wrote {out.name}, {out_summary.name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
