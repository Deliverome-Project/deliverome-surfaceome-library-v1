"""Collect the ortholog picks whose GPI anchor is inherited rather than measured.

The architecture attribute on D1 is annotated per human gene, and D1's own ortholog topology
records ``topology_projection_source: projected_from_human_canonical`` — so for an ortholog
of a GPI-anchored gene, nothing in hand measures whether that ortholog is anchored too.

The within-species C-terminal retention test used for isoforms cannot answer it either. The
GPI signal is conserved as architecture, not as sequence: human NYX against mouse Nyx, which
differ in anchoring (PMID 15905181), scores 0.37 over the last 30 residues, while human
against mouse CD55, both anchored, scores 0.43. The populations overlap, so no identity
cutoff separates them, and a dedicated GPI-signal predictor is the only instrument left.

This writes the FASTA for that run. All ortholog picks on GPI genes go in, not only the ones
DeepTMHMM calls soluble, so the predictor's agreement can be checked where the topology call
already settles the question.

    uv run --locked python build_ortholog_gpi_fasta.py
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from pathlib import Path

from gpi import curated_gpi_genes

UNIPROT_TABLES = (
    "uniprot_mouse_reviewed_isoforms.tsv.gz",
    "uniprot_cyno_all_isoforms.tsv.gz",
)
ENSEMBL_PEPTIDES = (
    "Mus_musculus.GRCm39.pep.all.fa.gz",
    "Macaca_fascicularis.Macaca_fascicularis_6.0.pep.all.fa.gz",
)


def sequences(cache: Path, extra: Path | None) -> dict[str, str]:
    """Ortholog sequences by accession, from the D1 dump then the local caches."""
    seq: dict[str, str] = {}
    if extra and extra.exists():
        for r in json.loads(extra.read_text()):
            if r.get("sequence"):
                seq[r["uniprot_acc_full"]] = r["sequence"]
    for name in UNIPROT_TABLES:
        path = cache / name
        if not path.exists():
            continue
        with gzip.open(path, "rt") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r.get("Sequence"):
                    seq.setdefault(r["Entry"], r["Sequence"])
    for name in ENSEMBL_PEPTIDES:
        path = cache / name
        if not path.exists():
            continue
        with gzip.open(path, "rt") as fh:
            pid, buf = None, []
            for line in fh:
                if line.startswith(">"):
                    if pid:
                        seq.setdefault(pid, "".join(buf))
                    buf = []
                    pid = line[1:].split()[0].split(".")[0]
                elif pid:
                    buf.append(line.strip())
            if pid:
                seq.setdefault(pid, "".join(buf))
    return seq


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cache-dir", default="~/.cache/deliverome/appris")
    ap.add_argument("--d1-sequences", default="", help="optional JSON dump of topology_public")
    ap.add_argument("--out", default="ortholog_gpi_check.fasta")
    args = ap.parse_args()

    gpi = curated_gpi_genes(args.cache_dir)
    seq = sequences(
        Path(args.cache_dir).expanduser(),
        Path(args.d1_sequences).expanduser() if args.d1_sequences else None,
    )

    recs, missing = [], []
    for species in ("mouse", "cyno"):
        for r in csv.DictReader(open(f"ortholog_repick_{species}.tsv"), delimiter="\t"):
            if r["human_symbol"] not in gpi:
                continue
            acc = r["picked_acc"]
            s = seq.get(acc) or seq.get(acc.split("-")[0])
            if not s:
                missing.append(f"{species}|{r['human_symbol']}|{acc}")
                continue
            # The header carries species and human gene so the result joins straight back.
            recs.append((f"{species}|{r['human_symbol']}|{acc}", s.rstrip("*")))

    Path(args.out).write_text("".join(f">{n}\n{s}\n" for n, s in recs))
    print(f"wrote {args.out} — {len(recs)} sequences, {len(missing)} unresolved")
    print(f"  longest {max(len(s) for _, s in recs)} aa, {sum(len(s) for _, s in recs):,} residues")
    if missing:
        print("  unresolved:", ", ".join(missing[:10]))


if __name__ == "__main__":
    main()
