"""Give every isoform-cohort form, and its gene's reference, a DeepTMHMM v1 topology.

The isoform gate (``gate_isoform_cohort.py``) judges a form by DeepTMHMM v1 on the form's own
sequence and fails a form with no v1 row as unverifiable. In the panel of record 931 of the
1,894 cohort forms failed for that reason alone (Becca, 2026-11-20: "well then we need to
address this first"). This closes the gap in three rerunnable steps, the pattern of
``publish_topology_supplement.py``:

1. ``fasta`` writes every form protein (and every gene reference) of a cohort table that has
   no v1 row in the base snapshot, one record per distinct sequence, named
   ``<form id>|<gene>|human|<sha256>``. The form protein is the gate's own
   (``gate_isoform_cohort.form_proteins``): the UniProt 2026-10-11 isoform when the row names
   an accession, else the Ensembl 115 translation.
2. ``modal run library/l1/proteoform-selection/modal/deeptmhmm_v1_app.py::sweep`` predicts them
   with the licensed DeepTMHMM 1.0 package on the ``deeptmhmm-v1`` Modal volume. FASTA,
   3-line output and the runtime stamp are committed under ``topology_runs/``.
3. ``publish`` appends one row per predicted sequence to the base snapshot under a fresh
   ``(topology_version, cohort)`` pair -- cohort ``human_isoforms``, version from the run's
   JSON -- asserts the row count on both sides and that no existing row changed, and writes a
   NEW external key. It does not bump ``topology.SNAPSHOT_KEY``: the library's tables of
   record stay on the pinned snapshot until that is decided. Consumers here pass the new key
   explicitly (``regate_isoform_forms.py``).

    uv run --locked python isoform_v1_topology.py fasta --cohort pinned \\
        --out topology_runs/2026-11-20_isoform_cohort.fasta
    uv run --locked python isoform_v1_topology.py publish --run 2026-11-20_isoform_cohort --write

Decision: Becca Carlson, 2026-11-20 (review comment 27 on the L1 design brief). Implemented by
the Claude Code session l1-isoform-universe (subagent of the brief review).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import tempfile
from pathlib import Path

import pandas as pd

from deliverome_analysis.surfaceome import topology

HERE = Path(__file__).resolve().parent
RUNS = HERE / "topology_runs"
RUN_PREFIX = "sprint4/s4e9-orf-design-fg-panel/topology-runs/"
SNAPSHOT_PREFIX = "accessible-surfaceome/topology/"
COHORT = "human_isoforms"
TOOL_VERSION = "deeptmhmm-1.0-academic"
SURFACE_LABELS = {"TM", "SP+TM"}


def sha256_file(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_base(key: str | None, sha: str | None) -> tuple[pd.DataFrame, str, str]:
    from deliverome_analysis.core.s3 import external

    key = key or topology.SNAPSHOT_KEY
    sha = sha or topology.SNAPSHOT_SHA256
    path = Path(external.download(key))
    if sha256_file(path) != sha:
        raise SystemExit(f"{key}: SHA-256 does not match {sha}")
    return pd.read_parquet(path), key, sha


def cohort_sequences(cohort: str, base_key: str | None, base_sha: str | None) -> dict:
    """{name: sequence} for the forms and references of ``cohort`` with no v1 row."""
    import gate_isoform_cohort as g

    work = Path(tempfile.mkdtemp(prefix="isoform-v1-"))
    inp = g.load_inputs(
        work,
        cohort_path=None if cohort == "pinned" else Path(cohort),
        snapshot_key=base_key,
        snapshot_sha256=base_sha,
    )
    proteins = g.form_proteins(inp)
    gene_of = {g.form_id(r): r["gene"] for r in inp["cohort"]}
    wanted: dict[str, tuple[str, str]] = {}
    for fid, seq in proteins.items():
        if seq:
            wanted.setdefault(topology.normalize(seq), (fid, gene_of[fid]))
    genes = set(gene_of.values())
    for gene, (did, seq) in sorted(inp["refs"].items()):
        if gene in genes:
            wanted.setdefault(topology.normalize(seq), (f"reference:{did}", gene))
    report = topology.resolve(
        {s: s for s in wanted}, tool=topology.TOOL_V1, snapshot=inp["snapshot"]
    )
    gap = report.loc[~report["found"], "key"].tolist()
    out = {}
    for seq in gap:
        fid, gene = wanted[seq]
        out[f"{fid}|{gene}|human|{topology.sequence_sha256(seq)}"] = seq
    print(
        f"{len(wanted)} distinct sequences ({len(proteins)} forms + references); "
        f"{len(out)} have no DeepTMHMM v1 row"
    )
    return out


def write_fasta(args) -> int:
    seqs = cohort_sequences(args.cohort, args.base_key, args.base_sha256)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    if args.out.exists():
        raise SystemExit(f"{args.out} exists; runs use fresh names")
    args.out.write_text("".join(f">{n}\n{s}\n" for n, s in seqs.items()))
    print(f"wrote {args.out} ({len(seqs)} sequences, {sum(map(len, seqs.values()))} residues)")
    return 0


def parse_3line(text: str) -> list[dict]:
    lines, out = text.splitlines(), []
    for i in range(len(lines) - 2):
        if lines[i].startswith(">"):
            name, label = lines[i][1:].rsplit("|", 1)
            fid, gene, species, sha = (x.strip() for x in name.split("|"))
            out.append(
                {
                    "form_id": fid,
                    "gene_symbol": gene,
                    "species": species,
                    "sha": sha,
                    "label": label.strip(),
                    "sequence": lines[i + 1].strip(),
                    "topology": lines[i + 2].strip(),
                }
            )
    return out


def publish(args) -> int:
    from publish_topology_supplement import derive

    from deliverome_analysis.core.s3 import external, processed

    meta = json.loads((RUNS / f"{args.run}.json").read_text())
    fasta, three = RUNS / f"{args.run}.fasta", RUNS / f"{args.run}.3line"
    runtime = json.loads((RUNS / f"{args.run}.3line.runtime.json").read_text())
    if runtime.get("failed_shards"):
        raise SystemExit(f"{args.run}: {runtime['failed_shards']} shard(s) failed; rerun them")
    if runtime["submitted_fasta_sha256"] != sha256_file(fasta):
        raise SystemExit(f"{args.run}: the run did not read the committed FASTA")
    submitted = {
        n.split("|")[3]
        for n in fasta.read_text().splitlines()
        if n.startswith(">")
        for n in [n[1:]]
    }
    preds = parse_3line(three.read_text())
    if {p["sha"] for p in preds} != submitted:
        raise SystemExit(f"{args.run}: predictions do not cover exactly the submitted sequences")
    base, base_key, _ = load_base(args.base_key, args.base_sha256)
    version = meta["topology_version"]
    cols = list(base.columns)
    rows = []
    for p in preds:
        seq = topology.normalize(p["sequence"])
        if topology.sequence_sha256(seq) != p["sha"]:
            raise SystemExit(f"{p['form_id']}: sequence does not match its header sha")
        if len(p["topology"]) != len(seq):
            raise SystemExit(f"{p['form_id']}: topology and sequence lengths differ")
        acc = p["form_id"].split(":", 1)[1]
        row = dict.fromkeys(cols)
        row.update(
            {
                "topology_version": version,
                "cohort": COHORT,
                "uniprot_acc": acc.split("-")[0] if not acc.startswith("ENST") else "",
                "uniprot_acc_full": p["form_id"],
                "isoform_id": acc,
                "gene_symbol": p["gene_symbol"],
                "species": "human",
                "is_canonical": 0,
                "protein_length": len(seq),
                "deeptmhmm_label": p["label"],
                "per_residue_topology": p["topology"],
                "predicted_surface_membrane": int(p["label"] in SURFACE_LABELS),
                "predicted_secreted": int(p["label"] == "SP"),
                "tool_version": runtime.get("tool_version", TOOL_VERSION),
                "retrieved_at": meta["run_completed"],
                "sequence_sha256": p["sha"],
                "sequence_length": len(seq),
                **derive(p["topology"]),
            }
        )
        rows.append(row)
    if not all(r["tool_version"].startswith(topology.V1_TOOL_PREFIX) for r in rows):
        raise SystemExit("run is not stamped as DeepTMHMM v1")
    clash = base[(base["topology_version"] == version) & (base["cohort"] == COHORT)]
    if len(clash):
        raise SystemExit(f"({version}, {COHORT}) already has {len(clash)} rows; use a fresh pair")
    add = pd.DataFrame(rows, columns=cols).astype(base.dtypes.to_dict())
    new = pd.concat([base, add], ignore_index=True)
    assert len(new) == len(base) + len(rows), "row count changed unexpectedly"
    assert new.iloc[: len(base)].equals(base), "an existing snapshot row changed"
    v1 = topology._with_tool(base, topology.TOOL_V1)["sequence_sha256"]
    assert not set(v1) & {r["sequence_sha256"] for r in rows}, "a sequence already had v1"
    key = f"{SNAPSHOT_PREFIX}{meta['snapshot']}/topology.parquet"
    labels = pd.Series([r["deeptmhmm_label"] for r in rows]).value_counts().to_dict()
    print(f"{len(rows)} row(s) as ({version}, {COHORT}); {base_key}: {len(base)} -> {len(new)}")
    print(f"labels: {labels}")
    if not args.write:
        print(f"dry run -- would write external://{key}")
        return 0
    if external.exists(key):
        raise SystemExit(f"external://{key} exists; the external bucket is never overwritten")
    out = Path(tempfile.mkdtemp()) / "topology.parquet"
    new.to_parquet(out, index=False)
    digest = sha256_file(out)
    external.upload(str(out), key)
    for path in (fasta, three, RUNS / f"{args.run}.3line.runtime.json", RUNS / f"{args.run}.json"):
        processed.upload(str(path), RUN_PREFIX + path.name)
    meta |= {
        "base_snapshot": base_key,
        "snapshot_key": key,
        "snapshot_sha256": digest,
        "rows_before": len(base),
        "rows_after": len(new),
        "rows_added": len(rows),
        "labels": labels,
    }
    (RUNS / f"{args.run}.published.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(f"wrote external://{key} sha256 {digest}")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)
    f = sub.add_parser("fasta")
    f.add_argument("--cohort", default="pinned", help="'pinned' or a cohort TSV")
    f.add_argument("--out", type=Path, required=True)
    p = sub.add_parser("publish")
    p.add_argument("--run", required=True, help="basename under topology_runs/")
    p.add_argument("--write", action="store_true")
    for s in (f, p):
        s.add_argument("--base-key", default=None, help="snapshot to extend (default: pinned)")
        s.add_argument("--base-sha256", default=None)
    args = ap.parse_args(argv)
    return write_fasta(args) if args.cmd == "fasta" else publish(args)


if __name__ == "__main__":
    raise SystemExit(main())
