"""Collect ortholog picks that carry no DeepTMHMM call yet.

A pick acquires a topology call either from the published D1 cohorts or from a previous
Modal sweep over the changed picks. Two groups fall through: picks the anchor guard moved
after the last sweep, and unchanged picks whose accession was never in a published cohort.
Both need a call before the surface classification is complete.

    uv run --locked python build_missing_topology_fasta.py --d1-labels <dump dir>
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from build_ortholog_gpi_fasta import sequences

VERSIONS = (
    "topo_2026_10_09_ortholog_repick",
    "topo_2026_09_rescue",
    "topo_2026_05_25",
    "topo_2026_05_16",
)


def known_labels(dump: Path) -> set[str]:
    have: set[str] = set()
    for v in VERSIONS:
        p = dump / f"tp_{v}.json"
        if not p.exists():
            continue
        for r in json.loads(p.read_text())[0]["results"]:
            have.update(x for x in (r["uniprot_acc_full"], r["uniprot_acc"]) if x)
    p = dump / "te_all.json"
    if p.exists():
        for r in json.loads(p.read_text())[0]["results"]:
            have.update(
                x for x in (r["ensembl_protein_id"], r["ensembl_transcript_id"]) if x
            )
    for name in ("ortholog_repick_changed.3line",):
        f = Path(name)
        if not f.exists():
            continue
        lines = f.read_text().splitlines()
        for i in range(len(lines) - 2):
            if lines[i].startswith(">"):
                nm, _lab = lines[i][1:].rsplit("|", 1)
                have.add(nm.split("|")[-1].strip())
    return have


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cache-dir", default="~/.cache/deliverome/appris")
    ap.add_argument("--d1-labels", required=True, help="directory holding the D1 dumps")
    ap.add_argument("--d1-sequences", default="")
    ap.add_argument(
        "--picked", nargs="*", default=[],
        help="FASTAs from `repick_orthologs.py --dump-picked`. These carry the unreviewed "
             "TrEMBL sequences that no cached reviewed table holds.",
    )
    ap.add_argument("--out", default="ortholog_missing_topology.fasta")
    args = ap.parse_args()

    have = known_labels(Path(args.d1_labels))
    seq = sequences(
        Path(args.cache_dir).expanduser(),
        Path(args.d1_sequences).expanduser() if args.d1_sequences else None,
    )
    for f in args.picked:
        name, buf = None, []
        for line in Path(f).read_text().splitlines():
            if line.startswith(">"):
                if name:
                    seq.setdefault(name, "".join(buf))
                name, buf = line[1:].split("|")[-1].strip(), []
            elif line.strip():
                buf.append(line.strip())
        if name:
            seq.setdefault(name, "".join(buf))

    recs, missing = [], []
    for species in ("mouse", "cyno"):
        for r in csv.DictReader(open(f"ortholog_repick_{species}.tsv"), delimiter="\t"):
            acc = r["picked_acc"]
            if r["picked_topology_label"] or acc in have or acc.split("-")[0] in have:
                continue
            s = seq.get(acc) or seq.get(acc.split("-")[0])
            if not s:
                missing.append(f"{species}|{r['human_symbol']}|{acc}")
                continue
            recs.append((f"{species}|{r['human_symbol']}|{acc}", s.rstrip("*")))

    Path(args.out).write_text("".join(f">{n}\n{s}\n" for n, s in recs))
    print(f"wrote {args.out} — {len(recs)} picks needing topology")
    if missing:
        print(f"  no sequence for {len(missing)}: {', '.join(missing)}")
    if recs:
        print(f"  longest {max(len(s) for _, s in recs)} aa, "
              f"{sum(len(s) for _, s in recs):,} residues")


if __name__ == "__main__":
    main()
