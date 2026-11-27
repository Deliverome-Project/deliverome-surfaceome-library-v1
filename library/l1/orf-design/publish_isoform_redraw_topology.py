"""Carry the redrawn isoform members' DeepTMHMM v1 rows into the shared topology snapshot.

The 2026-11-22 isoform redraw (``apply_panel_redraw_2026_11_22.py``) adds 81 L1 isoform
members. The rebuilt gate judged every one of them on a DeepTMHMM v1 row in its own pinned
snapshot (``isoform_cohort_rebuild_summary.json``: rows ``(topo_2026_11_20_isoform_cohort,
human_isoforms)``, run by ``rebuild_isoform_cohort.py``). For 40 of the members, the designed
protein has no v1 row in the shared snapshot ``topology.SNAPSHOT_KEY``, which was built from
a base that never took those rows. Downstream code (``l1_tag_paths.py``, ``topology.resolve``
on the default snapshot) therefore skips them.

No new prediction is needed. This script copies exactly those rows (matched by the designed
protein's SHA-256, the v1 tool stamp checked), appends them to the current shared snapshot
under a fresh ``(topology_version, cohort)`` pair, ``(topo_2026_11_22_isoform_redraw,
human_isoforms)``, and asserts the row count on both sides. Every other column is the
source row's own, including ``tool_version`` and ``retrieved_at``, so a merged row still
says which run produced it. The result goes to a NEW external key, following the
conventions of ``publish_topology_supplement.py``. Adopt it by bumping
``topology.SNAPSHOT_KEY`` / ``SNAPSHOT_SHA256``.

Implemented by Claude Code session l1-review-export-scope (worker fork), 2026-11-22, for the
redraw Becca Carlson decided on 2026-11-22.

    uv run --locked python publish_isoform_redraw_topology.py --mode refill            # dry run
    uv run --locked python publish_isoform_redraw_topology.py --mode refill --write    # upload
    uv run --locked python publish_isoform_redraw_topology.py --mode redraw --check    # replay

``--mode refill`` (2026-11-22, same day) does the same for PAM and KIT, two members added by
``refill_isoform_panel_2026_11_22.py``, on top of the MANE-GTF snapshot.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
import tempfile
from pathlib import Path

import pandas as pd

from deliverome_analysis.cloning import library_design as ld
from deliverome_analysis.core.s3 import external
from deliverome_analysis.surfaceome import topology

HERE = Path(__file__).resolve().parent
MANIFEST = HERE / "candidate_manifest_tranches.tsv"
GATE_PIN = HERE / "isoform_cohort_rebuild_summary.json"
SOURCE_VERSION = "topo_2026_11_20_isoform_cohort"
COHORT = "human_isoforms"
#: Each publication: the change log naming the added members, the snapshot it builds on
#: (pinned, so ``--check`` replays it after ``topology.SNAPSHOT_KEY`` moves on), its fresh
#: ``topology_version`` and the new key. ``refill`` (2026-11-22) carries PAM and KIT, the two
#: refills of ``refill_isoform_panel_2026_11_22.py`` whose v1 rows were also only in the gate's
#: snapshot, onto the MANE-GTF snapshot.
PREFIX = "accessible-surfaceome/topology/"
MODES = {
    "redraw": {
        "changes": HERE / "panel_redraw_2026-11-22_changes.tsv",
        "base": PREFIX + "2026-11-22-round2-additions/topology.parquet",
        "base_sha256": None,  # the shared snapshot at the time; not pinned when published
        "version": "topo_2026_11_22_isoform_redraw",
        "key": PREFIX + "2026-11-22-isoform-redraw/topology.parquet",
    },
    "refill": {
        "changes": HERE / "panel_refill_2026-11-22_changes.tsv",
        "base": PREFIX + "2026-11-22-mane-gtf/topology.parquet",
        "base_sha256": "5f561945ee6ff47a665d157cd27874170e88e508b1d19c5e41935d6908d80cf7",
        "version": "topo_2026_11_22_isoform_refill",
        "key": PREFIX + "2026-11-22-isoform-refill/topology.parquet",
    },
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build(work: Path, mode: dict) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """(base snapshot, merged rows, new snapshot)."""
    VERSION = mode["version"]
    with mode["changes"].open() as fh:
        added = [
            c["candidate_id"]
            for c in csv.DictReader(fh, delimiter="\t")
            if c.get("panel", "isoform") == "isoform" and c["change"] == "added"
        ]
    design = ld.design_ids(MANIFEST)
    proteins = {
        design[cid]: v["sequence"]
        for cid, v in ld.designed_proteins(MANIFEST, library=ld.LIBRARY_L1).items()
    }
    # Members retired since (e.g. by the MANE-GTF switch) have no designed protein to look up.
    seqs = {d: proteins[d] for d in added if d in proteins}
    base_path = Path(external.download(mode["base"], work / "base.parquet"))
    if mode["base_sha256"] and sha256(base_path) != mode["base_sha256"]:
        raise SystemExit(f"{mode['base']} differs from its pinned SHA-256")
    base = pd.read_parquet(base_path)
    have = topology.resolve(seqs, tool=topology.TOOL_V1, snapshot=base)
    missing = have.loc[~have["found"], ["key", "sequence_sha256"]]

    pin = json.loads(GATE_PIN.read_text())["topology_snapshot"]
    path = Path(external.download(pin["key"], work / "gate.parquet"))
    if sha256(path) != pin["sha256"]:
        raise SystemExit(f"{pin['key']} differs from its pinned SHA-256")
    gate = pd.read_parquet(path)
    src = gate[(gate["topology_version"] == SOURCE_VERSION) & (gate["cohort"] == COHORT)]
    src = topology._with_tool(src, topology.TOOL_V1)
    rows = []
    for key, sha in missing.itertuples(index=False):
        hit = src[src["sequence_sha256"] == sha]
        if len(hit) != 1:
            raise SystemExit(f"{key}: {len(hit)} v1 rows in the gate snapshot")
        rows.append(hit.iloc[0].to_dict() | {"topology_version": VERSION, "cohort": COHORT})
    merged = pd.DataFrame(rows, columns=list(base.columns)).astype(base.dtypes.to_dict())

    clash = base[(base["topology_version"] == VERSION) & (base["cohort"] == COHORT)]
    if len(clash):
        raise SystemExit(f"({VERSION}, {COHORT}) already has {len(clash)} rows")
    new = pd.concat([base, merged], ignore_index=True)
    assert len(new) == len(base) + len(merged), "row count changed unexpectedly"
    after = topology.resolve(seqs, tool=topology.TOOL_V1, snapshot=new)
    assert after["found"].all(), "an added member still has no v1 row"
    return base, merged, new


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--mode", choices=sorted(MODES), default="refill")
    group = ap.add_mutually_exclusive_group()
    group.add_argument("--write", action="store_true")
    group.add_argument("--check", action="store_true", help="compare with the published key")
    a = ap.parse_args(argv)
    mode = MODES[a.mode]
    KEY, VERSION = mode["key"], mode["version"]
    work = Path(tempfile.mkdtemp(prefix="isoform-redraw-topology-"))
    base, merged, new = build(work, mode)
    print(
        f"base {mode['base']}: {len(base)} rows; +{len(merged)} as "
        f"({VERSION}, {COHORT}) -> {len(new)}"
    )
    for r in merged.itertuples():
        print(
            f"  {r.gene_symbol:<10} {r.uniprot_acc_full:<22} {r.deeptmhmm_label:<6} {r.tool_version}"
        )
    if a.check:
        if not external.exists(KEY):
            print(f"external://{KEY} not published")
            return 1
        got = pd.read_parquet(external.download(KEY, work / "published.parquet"))
        same = got.reset_index(drop=True).equals(new.reset_index(drop=True))
        print("matches the published snapshot" if same else "DIFFERS from the published snapshot")
        return 0 if same else 1
    if not a.write:
        print(f"dry run -- would write external://{KEY}")
        return 0
    if external.exists(KEY):
        raise SystemExit(f"external://{KEY} exists; the external bucket is never overwritten")
    out = work / "topology.parquet"
    new.to_parquet(out, index=False)
    external.upload(str(out), KEY)
    print(f"wrote external://{KEY}\nsha256 {sha256(out)}")
    print(f"now set topology.SNAPSHOT_KEY = {KEY!r} and SNAPSHOT_SHA256 to the sha above")
    return 0


if __name__ == "__main__":
    sys.exit(main())
