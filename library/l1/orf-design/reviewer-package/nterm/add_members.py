"""Bring members added to the "N after SP" rung after the 2026-11-12 build through this package.

Gene additions (``gene_list_additions.tsv``) and redraws land members on the rung after its
tables were built. Rerunning the package's scripts from scratch would also rewrite rows of
members that have since left the rung, so the additions go through the same scripts and only
their rows are merged into the committed tables, every other row kept byte for byte.

Subcommands:

``snapshot --members ID,... --date DATE``
    D1 predictor rows (the pinned versions of ``snapshot_inputs.py``) for the human accessions
    whose topology sequence is a new member's designed protein, merged with the current
    snapshot (``nterm_placement.SNAPSHOT``) and published under ``nterm/<DATE>/``. Move
    ``nterm_placement.SNAPSHOT`` to it.
``merge --members ID,... --from DIR``
    Merge the new members' rows of each table a script wrote into ``DIR`` (run the script with
    its output pointed there) into the committed table: rows are matched on ``candidate_id`` (else the
    first column); committed rows of other members are kept verbatim; new rows go at the end.

Written 2026-11-22 for the round-2 additions (l1-p3178..l1-p3214) by Claude Code session
l1-review-export-scope (subagent); the 2026-11-21 co-receptor additions went through the same
steps by hand.

    uv run --locked python add_members.py snapshot --members l1-p3178,... --date 2026-11-22
    uv run --locked python add_members.py merge --members l1-p3178,... --from /tmp/nterm-out
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import nterm_placement  # noqa: E402
import snapshot_inputs as si  # noqa: E402


def snapshot(ids: list[str], date: str, out: Path) -> int:
    from deliverome_analysis.core.s3 import processed
    from deliverome_analysis.surfaceome import topology

    lane = {c: m for c, m in si.members().items() if c in ids}
    missing = sorted(set(ids) - set(lane))
    assert not missing, f"not on the N-after-SP rung: {missing}"
    snap = topology.load_snapshot()
    hit = snap[snap["sequence_sha256"].isin({m["sha"] for m in lane.values()})]
    accs = set(hit[hit["cohort"].isin(si.HUMAN_COHORTS)]["uniprot_acc_full"])
    print(f"{len(lane)} members; {len(accs)} accessions with their sequence")
    new = si.fetch(lane, accs, si.d1_query()) if accs else []
    with gzip.open(processed.download(nterm_placement.SNAPSHOT + si.ROWS_NAME), "rt") as fh:
        old = [json.loads(line) for line in fh]
    seen = {json.dumps(r, sort_keys=True) for r in old}
    added = [r for r in new if json.dumps(r, sort_keys=True) not in seen]
    out.mkdir(parents=True, exist_ok=True)
    si.write_rows(old + added, out / si.ROWS_NAME)
    canary = processed.download(nterm_placement.SNAPSHOT + si.CANARY_NAME)
    (out / si.CANARY_NAME).write_bytes(Path(canary).read_bytes())
    si.publish(out, date)
    print(
        f"{len(old)} rows kept, {len(added)} added; set nterm_placement.SNAPSHOT to nterm/{date}/"
    )
    return 0


def merge(ids: set[str], src: Path) -> int:
    for path in sorted(src.glob("*.tsv")):
        dest = HERE / path.name
        if not dest.exists():
            continue
        with dest.open() as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            cols, old = list(reader.fieldnames or []), list(reader)
        key = "candidate_id" if "candidate_id" in cols else cols[0]
        with path.open() as fh:
            new = [r for r in csv.DictReader(fh, delimiter="\t") if r[key] in ids]
        # every committed row of another member stays where it is; a new member's rows (one or
        # several per member) replace its committed ones and go at the end
        rows = [r for r in old if r[key] not in ids] + new
        before = dest.read_bytes()
        with dest.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n",
                               extrasaction="ignore")  # fmt: skip
            w.writeheader()
            w.writerows(rows)
        old_lines = before.decode().splitlines()
        new_lines = dest.read_text().splitlines()
        col = cols.index(key)
        kept_old = [ln for ln in old_lines[1:] if ln.split("\t")[col] not in ids]
        assert new_lines[1 : 1 + len(kept_old)] == kept_old, f"{dest.name}: another row changed"
        print(f"{dest.name}: {len(new)} row(s) for the new members")
    return 0


#: Per-member tables of this package (``candidate_id``-keyed) that retired members leave.
MEMBER_TABLES = (
    "nterm_members.tsv",
    "nterm_positions.tsv",
    "nterm_build_options.tsv",
    "nterm_site_picks.tsv",
    "overnight_inputs.tsv",
    "overnight_proposals.tsv",
    "overnight_verdicts.tsv",
    "round2_inputs.tsv",
    "nterm_decisions.tsv",
)


def drop_retired() -> int:
    """Remove the rows of members the manifest retires; every other row is kept verbatim.

    Their issued cassette encodings stay issued (``l1_arm_cassette_encodings.tsv`` is not
    touched here)."""
    with (HERE.parents[1] / "candidate_manifest_tranches.tsv").open() as fh:
        retired = {
            r["candidate_id"] for r in csv.DictReader(fh, delimiter="\t") if r["retired_reason"]
        }
    for name in MEMBER_TABLES:
        path = HERE / name
        lines = path.read_text().splitlines()
        col = lines[0].split("\t").index("candidate_id")
        keep = [lines[0]] + [ln for ln in lines[1:] if ln.split("\t")[col] not in retired]
        path.write_text("\n".join(keep) + "\n")
        print(f"{name}: dropped {len(lines) - len(keep)} row(s) of retired members")
    return 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("snapshot")
    s.add_argument("--members", required=True)
    s.add_argument("--date", required=True)
    s.add_argument("--out", type=Path, default=Path.home() / "nterm-snapshot-add")
    m = sub.add_parser("merge")
    m.add_argument("--members", required=True)
    m.add_argument("--from", dest="src", type=Path, required=True)
    sub.add_parser("drop-retired")
    a = ap.parse_args(argv)
    if a.cmd == "drop-retired":
        return drop_retired()
    ids = a.members.split(",")
    if a.cmd == "snapshot":
        return snapshot(ids, a.date, a.out)
    return merge(set(ids), a.src)


if __name__ == "__main__":
    raise SystemExit(main())
