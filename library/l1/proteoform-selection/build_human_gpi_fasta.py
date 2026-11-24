"""Emit every human form of a GPI-anchored gene, for a GPI-anchor run.

D1's architecture attribute is a statement about the *gene*: that some form of it is
GPI-anchored. It does not say the form this library selected is that one, and often it is
not — NetGPI called 15 of our 148 curated-GPI references non-anchored, and for NCAM1, ACHE,
BCAN, TFPI and CD58 it is right, because the anchor lives on an isoform the canonical
reference is not (NCAM-120 against NCAM-140, the ACHE H subunit against the canonical T).

So the anchor has to be measured per form, and the library then has to carry whichever form
actually bears it. Only GPI genes are emitted -- the question is which of a known anchored
gene's forms carries the anchor, not whether some other gene has one.

Every alternate form goes in, not just the evidence-filtered cohort: if the anchor sits on a
form the cohort dropped, that is exactly what needs to be seen. Orthologs are run separately
by ``build_ortholog_gpi_fasta.py``.

    uv run --locked python build_human_gpi_fasta.py
"""

from __future__ import annotations

import argparse
import csv
import gzip
import re
from pathlib import Path

from gpi import curated_gpi_genes

MAX_RESIDUES = 10_000  # NetGPI rejects anything longer


def uniprot_sequences(cache: Path) -> dict[str, str]:
    seq: dict[str, str] = {}
    for name in ("uniprot_human_reviewed.tsv.gz", "uniprot_human_reviewed_isoforms.tsv.gz"):
        path = cache / name
        if not path.exists():
            continue
        with gzip.open(path, "rt") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r.get("Sequence"):
                    seq.setdefault(r["Entry"], r["Sequence"])
    return seq


def by_gene(cache: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    path = cache / "uniprot_human_reviewed.tsv.gz"
    if path.exists():
        with gzip.open(path, "rt") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                g = (r.get("Gene Names (primary)") or "").strip()
                if g and r.get("Sequence"):
                    out.setdefault(g, r["Sequence"])
    return out


def ensembl_sequences(cache: Path) -> dict[str, str]:
    seq: dict[str, str] = {}
    path = cache / "Homo_sapiens.GRCh38.pep.all.fa.gz"
    if not path.exists():
        return seq
    with gzip.open(path, "rt") as fh:
        enst, buf = None, []
        for line in fh:
            if line.startswith(">"):
                if enst:
                    seq.setdefault(enst, "".join(buf))
                buf = []
                m = re.search(r"transcript:(ENST\d+)", line)
                enst = m.group(1) if m else None
            elif enst:
                buf.append(line.strip())
        if enst:
            seq.setdefault(enst, "".join(buf))
    return seq


def three_line(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    if not path.exists():
        return out
    lines = path.read_text().splitlines()
    for i in range(len(lines) - 2):
        if lines[i].startswith(">"):
            nm, _lab = lines[i][1:].rsplit("|", 1)
            out[nm.strip()] = lines[i + 1].strip()
    return out


def write(path: Path, recs: list[tuple[str, str]]) -> None:
    kept = [(n, s) for n, s in recs if len(s) <= MAX_RESIDUES]
    path.write_text("".join(f">{n}\n{s.rstrip('*')}\n" for n, s in kept))
    skipped = len(recs) - len(kept)
    print(f"wrote {path} — {len(kept)} sequences"
          + (f", {skipped} over {MAX_RESIDUES} aa skipped" if skipped else ""))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cache-dir", default="~/.cache/deliverome/appris")
    ap.add_argument("--manifest", default="topology_shortlist.tsv")
    ap.add_argument("--three-line", default="topology_shortlist.3line")
    ap.add_argument("--universe", default="isoform_universe_scored.tsv")
    args = ap.parse_args()
    cache = Path(args.cache_dir).expanduser()
    gpi = curated_gpi_genes(args.cache_dir)

    acc_seq = uniprot_sequences(cache)
    gene_seq = by_gene(cache)
    enst_seq = ensembl_sequences(cache)
    fid_seq = three_line(Path(args.three_line))
    fid_for = {
        (m["gene"], m["length"]): m["fasta_id"]
        for m in csv.DictReader(open(args.manifest), delimiter="\t")
    }

    # ---- references: the form the precedence actually selected --------------------
    recs, missing = [], []
    for r in csv.DictReader(open("reference_selection.tsv"), delimiter="\t"):
        rung, sym = r["reference_rung"], r["symbol"]
        if sym not in gpi:
            continue
        s = None
        if r["selected_uniprot_equivalent"]:
            s = acc_seq.get(r["selected_uniprot_equivalent"])
        if s is None and rung == "mane" and r["mane_enst"]:
            s = enst_seq.get(r["mane_enst"])
        if s is None and rung == "appris" and r["appris_enst"]:
            s = enst_seq.get(r["appris_enst"])
        if s is None and r["current_accession"]:
            s = acc_seq.get(r["current_accession"])
        if s is None:
            s = gene_seq.get(sym)
        if s is None:
            missing.append(sym)
            continue
        recs.append((f"ref|{sym}|{r['selected_uniprot_equivalent'] or rung}", s))
    write(Path("netgpi_gpi_reference.fasta"), recs)
    if missing:
        print(f"  unresolved references: {len(missing)} — {', '.join(missing[:8])}")

    # ---- every alternate form of those genes ---------------------------------------
    recs, missing = [], []
    seen: set[tuple[str, str]] = set()
    for r in csv.DictReader(open(args.universe), delimiter="\t"):
        if r["gene"] not in gpi or r["is_reference"] == "True":
            continue
        s = acc_seq.get(r["uniprot"]) if r["uniprot"] else None
        if s is None:
            fid = fid_for.get((r["gene"], r["length"]))
            s = fid_seq.get(fid) if fid else None
        if s is None:
            for enst in (r["ensembl_transcripts"] or "").split(";"):
                s = enst_seq.get(enst.strip())
                if s:
                    break
        if s is None:
            missing.append(f"{r['gene']}:{r['uniprot'] or r['ensembl_transcripts']}")
            continue
        tag = r["uniprot"] or (r["ensembl_transcripts"] or "").split(";")[0]
        if (r["gene"], tag) in seen:
            continue
        seen.add((r["gene"], tag))
        recs.append((f"alt|{r['gene']}|{tag}", s))
    write(Path("netgpi_gpi_alt_forms.fasta"), recs)
    if missing:
        print(f"  unresolved alternate forms: {len(missing)} — {', '.join(missing[:8])}")


if __name__ == "__main__":
    main()
