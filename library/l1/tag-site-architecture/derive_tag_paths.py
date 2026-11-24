"""Which tag-placement rung each canonical library member lands on, from topology alone.

Reads the committed DeepTMHMM predictions in the surfaceome repo and assigns every
member the rung the deterministic pipeline would reach first:

    ecto-N (after signal-peptide cleavage, or at residue 1) -> ecto-C -> internal

Two members on the same rung can still need different handling, so the per-member table
also carries the TM count: a no-signal-peptide ecto-N member that is multipass is
Tedman's validated GPCR case, while the same rung with a single TM is a type III
protein, where one TM carries the whole topogenic decision and an N-terminal cassette
shifts the charge balance that sets its orientation.

Membership is the s1e3 library file, not the topology export. Those are different sets:
the ``human_canonical_non_hla`` export and the library overlap by only about 75%, so
reading the export as "the library" both analyses ~570 non-members and silently omits
~591 members. The library is the authority for who is in; topology is looked up for them.

Entry names in the 3line headers ("AAAT") are not gene symbols ("SLC1A5"), so symbols
come from the catalog by accession rather than from the FASTA header.

    uv run --locked python library/l1/tag-site-architecture/derive_tag_paths.py \
        --surfaceome-repo ~/Git/accessible-surfaceome/.claude/worktrees/tag-sites-468
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path

HERE = Path(__file__).parent
OUT = HERE / "tag_path_assignments.tsv"
PRED = "data/external/deeptmhmm_surfaceome_predictions/human_canonical_non_hla/predicted_topologies.3line"
CATALOG = "data/processed/catalog/whole_proteome_catalog.tsv"


def parse_3line(path: Path) -> dict[str, tuple[str, str]]:
    """acc -> (sequence, per-residue topology string)."""
    lines = [ln.rstrip("\n") for ln in path.read_text().splitlines()]
    out = {}
    for i in range(0, len(lines) - 2, 3):
        if lines[i].startswith(">"):
            out[lines[i].split("|")[1]] = (lines[i + 1], lines[i + 2])
    return out


def signal_peptide_end(topo: str) -> int:
    """Last residue of a LEADING 'S' run, or 0 when residue 1 is not a signal residue."""
    return 0 if not topo.startswith("S") else len(topo) - len(topo.lstrip("S"))


def classify(topo: str) -> tuple[str, int]:
    """(rung, signal-peptide length). The ladder is ecto-N, then ecto-C, then internal."""
    sp = signal_peptide_end(topo)
    if topo[0] == "O":
        return "N at residue 1", 0
    if sp and sp < len(topo) and topo[sp] == "O":
        return "N after SP", sp
    return ("C-terminal" if topo[-1] == "O" else "internal"), sp


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--surfaceome-repo", type=Path, required=True,
                    help="checkout of accessible-surfaceome holding the DeepTMHMM run")
    ap.add_argument("--extra-3line", action="append", default=[],
                    help="additional 3line files, relative to --surfaceome-repo, holding "
                         "topology for library members absent from the main export")
    args = ap.parse_args()

    recs = parse_3line(args.surfaceome_repo / PRED)
    for extra in args.extra_3line:
        for acc, rec in parse_3line(args.surfaceome_repo / extra).items():
            recs.setdefault(acc, rec)
    symbol = {}
    with (args.surfaceome_repo / CATALOG).open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["uniprot_acc"]:
                symbol[r["uniprot_acc"]] = r["hgnc_symbol"]

    members = set(json.loads((HERE / "library_members.json").read_text()))
    rows, seen = [], set()
    for acc, (seq, topo) in recs.items():
        gene = symbol.get(acc, "")
        if gene not in members or gene in seen:
            continue
        seen.add(gene)
        rung, sp = classify(topo)
        rows.append({
            "gene": gene, "accession": acc, "rung": rung,
            "n_tm": len(re.findall(r"M+", topo)), "signal_peptide_aa": sp,
            "length_aa": len(seq),
        })
    absent = sorted(members - seen)
    (HERE / "library_members_without_topology.txt").write_text("\n".join(absent) + "\n")
    rows.sort(key=lambda r: (r["rung"], r["gene"]))
    with OUT.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {OUT.name}: {len(rows):,} of {len(members):,} library members classified")
    print(f"  {len(absent)} with no topology in any supplied 3line -> "
          f"library_members_without_topology.txt")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
