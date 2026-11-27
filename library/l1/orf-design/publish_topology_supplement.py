"""Give every active library member a DeepTMHMM v1 topology, reproducibly.

A member whose exact designed protein has no v1 row in the pinned topology snapshot
(``deliverome_analysis.surfaceome.topology.SNAPSHOT_KEY``) cannot be placed: ``l1_tag_paths.py`` skips
it, and nothing may substitute the canonical's topology. This closes that gap in three
committed, rerunnable steps:

1. ``fasta`` writes the designed proteins of every active member with no v1 row, named
   ``<design_id>|<gene>|<species>|<sha256>``. The design id (``l1-pNNNN``), not the
   member id: member ids move until the order freeze, a snapshot row must not.
2. ``modal run library/l1/proteoform-selection/modal/deeptmhmm_v1_app.py::sweep`` predicts
   them with the licensed DeepTMHMM 1.0.24 package on a Modal volume. The run's FASTA and
   3-line output are committed under ``topology_runs/`` and published to S3.
3. ``publish`` derives each snapshot column from the 3-line output (the derivation
   ``upload_topology_ensembl.py`` uses), appends the rows to the current snapshot under a
   fresh ``(topology_version, cohort)`` pair, asserts the row count on both sides, and
   writes a NEW external key. Bump ``topology.SNAPSHOT_KEY`` to adopt it.

Rows go to cohort ``library_designed_orf``, the cohort for designed ORFs that exist in no
other proteoform space, keyed by design id as the 2026-10-14 rows are (those predate the
member renumbering, when every candidate id was a design id). D1
``topology_public`` is not written here: the viewer's production table is loaded from
accessible-surfaceome, so the snapshot leads D1 by these rows until that load is run.

    uv run --locked python publish_topology_supplement.py fasta --out topology_runs/<run>.fasta
    uv run --locked python publish_topology_supplement.py publish --run <run> [--write]

A protein may also be the decided-but-held restart of a member (``ortholog_start_trims.tsv``,
``status = held_pending_v1``): it is accepted only when its sha256 equals that row's
``restarted_protein_sha256``, since the restart cannot be applied before its v1 row exists.
Added 2026-11-12 for the cyno CHRM3, NIPAL4 and OR3A2 restarts (l1-reviewer-package, PR #283).
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
from pathlib import Path

import pandas as pd

from deliverome_analysis.cloning import library_design as ld
from deliverome_analysis.surfaceome import topology

HERE = Path(__file__).parent
MANIFEST = HERE / "candidate_manifest_tranches.tsv"
START_TRIMS = HERE / "ortholog_start_trims.tsv"
RUNS = HERE / "topology_runs"
RUN_PREFIX = "sprint4/s4e9-orf-design-fg-panel/topology-runs/"
SNAPSHOT_PREFIX = "accessible-surfaceome/topology/"
COHORT = "library_designed_orf"
TOOL_VERSION = "deeptmhmm-1.0.24"  # the package on the deeptmhmm-v1 Modal volume
SURFACE_LABELS = {"TM", "SP+TM"}
ORIENT = {"O": "extracellular", "I": "cytoplasmic"}


def derive(topo: str) -> dict:
    """Counts and orientations straight from the per-residue string."""

    def runs(ch: str) -> int:
        return sum(1 for j, c in enumerate(topo) if c == ch and (j == 0 or topo[j - 1] != ch))

    non_signal = [c for c in topo if c != "S"]
    first = next((c for c in non_signal if c in ORIENT), None)
    last = next((c for c in reversed(non_signal) if c in ORIENT), None)
    return {
        "tm_helix_count": runs("M"),
        "beta_strand_count": runs("B"),
        "signal_peptide_length": topo.count("S"),
        "ecd_length_residues": topo.count("O"),
        "icd_length_residues": topo.count("I"),
        "n_terminal_orientation": ORIENT.get(first, "indeterminate"),
        "c_terminal_orientation": ORIENT.get(last, "indeterminate"),
    }


def parse_3line(text: str) -> list[dict]:
    lines, out = text.splitlines(), []
    for i in range(0, len(lines) - 2):
        if lines[i].startswith(">"):
            name, label = lines[i][1:].rsplit("|", 1)
            cid, gene, species, sha = (x.strip() for x in name.split("|"))
            out.append(
                {
                    "candidate_id": cid,
                    "gene_symbol": gene,
                    "species": species,
                    "sha": sha,
                    "label": label.strip(),
                    "sequence": lines[i + 1].strip(),
                    "topology": lines[i + 2].strip(),
                }
            )
    return out


def write_fasta(out: Path) -> int:
    design = ld.design_ids(MANIFEST)
    proteins = {
        design[cid]: v for cid, v in ld.designed_proteins(MANIFEST, library=ld.LIBRARY_L1).items()
    }
    report = topology.resolve({did: v["sequence"] for did, v in proteins.items()})
    gap = sorted(report.loc[~report["found"], "key"])
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as fh:
        for cid in gap:
            p = proteins[cid]
            seq = topology.normalize(p["sequence"])
            fh.write(
                f">{cid}|{p['gene_symbol']}|{p['species']}|{topology.sequence_sha256(seq)}\n{seq}\n"
            )
    print(f"{len(gap)} active L1 member(s) without a v1 topology -> {out}")
    return len(gap)


def held_restarts() -> dict[str, str]:
    """design id -> restarted protein sha256, for restarts decided but held pending v1."""
    import csv

    with START_TRIMS.open() as fh:
        return {
            r["design_id"]: r["restarted_protein_sha256"]
            for r in csv.DictReader(fh, delimiter="\t")
            if r["status"] == "held_pending_v1"
        }


def publish(run: str, write: bool) -> int:
    from deliverome_analysis.core.s3 import external, processed

    meta = json.loads((RUNS / f"{run}.json").read_text())
    fasta, three = RUNS / f"{run}.fasta", RUNS / f"{run}.3line"
    preds = parse_3line(three.read_text())
    # Runs name each protein by its design id; key the members the same way.
    design = ld.design_ids(MANIFEST)
    proteins = {
        design[cid]: v for cid, v in ld.designed_proteins(MANIFEST, library=ld.LIBRARY_L1).items()
    }
    version = meta["topology_version"]
    held = held_restarts()

    base = topology.load_snapshot()
    cols = list(base.columns)
    rows = []
    for p in preds:
        seq = topology.normalize(p["sequence"])
        if topology.sequence_sha256(seq) != p["sha"]:
            raise SystemExit(f"{p['candidate_id']}: sequence does not match its header sha")
        if held.get(p["candidate_id"]) == p["sha"]:
            pass  # the member's decided restart, applied once this row exists
        elif meta.get("pending_additions", {}).get(p["candidate_id"]) == p["sha"]:
            # a dated gene addition (add_library_genes.py), which joins the manifest only once
            # its designed protein has this row; the run records the protein it will build
            pass
        elif topology.normalize(proteins[p["candidate_id"]]["sequence"]) != seq:
            raise SystemExit(f"{p['candidate_id']}: not the member's designed protein")
        if len(p["topology"]) != len(seq):
            raise SystemExit(f"{p['candidate_id']}: topology and sequence lengths differ")
        row = dict.fromkeys(cols)
        row.update(
            {
                "topology_version": version,
                "cohort": COHORT,
                "uniprot_acc": "l1",
                "uniprot_acc_full": p["candidate_id"],
                "isoform_id": p["candidate_id"],
                "gene_symbol": p["gene_symbol"],
                "species": p["species"],
                "is_canonical": 0,
                "protein_length": len(seq),
                "deeptmhmm_label": p["label"],
                "per_residue_topology": p["topology"],
                "predicted_surface_membrane": int(p["label"] in SURFACE_LABELS),
                "predicted_secreted": int(p["label"] == "SP"),
                # the run's own stamp when it recorded one (runtime.json), else the default
                "tool_version": meta.get("tool_version", TOOL_VERSION),
                "retrieved_at": meta["run_completed"],
                "sequence_sha256": p["sha"],
                "sequence_length": len(seq),
                **derive(p["topology"]),
            }
        )
        rows.append(row)

    clash = base[(base["topology_version"] == version) & (base["cohort"] == COHORT)]
    if len(clash):
        raise SystemExit(f"({version}, {COHORT}) already has {len(clash)} rows; use a fresh pair")
    new = pd.concat(
        [base, pd.DataFrame(rows, columns=cols).astype(base.dtypes.to_dict())], ignore_index=True
    )
    assert len(new) == len(base) + len(rows), "row count changed unexpectedly"
    key = f"{SNAPSHOT_PREFIX}{meta['snapshot']}/topology.parquet"
    print(f"{len(rows)} row(s) as ({version}, {COHORT}); snapshot {len(base)} -> {len(new)}")
    for r in rows:
        print(
            f"  {r['uniprot_acc_full']} {r['gene_symbol']:<8} {r['species']:<6} "
            f"{r['deeptmhmm_label']:<6} tm={r['tm_helix_count']} "
            f"N {r['n_terminal_orientation']} / C {r['c_terminal_orientation']}"
        )
    if not write:
        print(f"\ndry run -- would write external://{key}")
        return 0
    if external.exists(key):
        raise SystemExit(f"external://{key} exists; the external bucket is never overwritten")
    out = HERE / f"_{meta['snapshot']}.parquet"
    new.to_parquet(out, index=False)
    external.upload(str(out), key)
    out.unlink()
    for path in (fasta, three, RUNS / f"{run}.json"):
        processed.upload(str(path), RUN_PREFIX + path.name)
    digest = hashlib.sha256(io.BytesIO(three.read_bytes()).getvalue()).hexdigest()
    print(
        f"wrote external://{key}; run files under processed://{RUN_PREFIX} "
        f"(3line sha256 {digest[:12]})"
    )
    print(f"now set topology.SNAPSHOT_KEY = {key!r}")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)
    f = sub.add_parser("fasta")
    f.add_argument("--out", type=Path, required=True)
    p = sub.add_parser("publish")
    p.add_argument("--run", required=True, help="basename under topology_runs/")
    p.add_argument("--write", action="store_true")
    args = ap.parse_args(argv)
    if args.cmd == "fasta":
        write_fasta(args.out)
        return 0
    return publish(args.run, args.write)


if __name__ == "__main__":
    raise SystemExit(main())
