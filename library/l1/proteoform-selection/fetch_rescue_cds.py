"""Fetch the canonical CDS for each expression-rescue receptor.

The rescue arms were chosen by symbol; building the constructs needs their
coding sequence. The record was made through Ensembl REST (release 115 at the
time). ``--offline`` takes the same answer from the pinned release-115 files
instead -- the canonical transcript from the GTF's ``Ensembl_canonical`` tag and
its CDS from the cds FASTA -- so it cannot drift as REST moves to later releases.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
from pathlib import Path

from deliverome_analysis import seqtools


def offline_records(symbols: list[str]) -> dict[str, dict]:
    """Canonical CDS per symbol from the pinned Ensembl 115 GTF and cds FASTA."""
    import candidate_inputs as inputs

    want, canon = set(symbols), {}
    with gzip.open(inputs.gtf("human"), "rt") as fh:
        for line in fh:
            if "\ttranscript\t" not in line or 'tag "Ensembl_canonical"' not in line:
                continue
            attrs = dict(re.findall(r'(\w+) "([^"]*)"', line.split("\t")[8]))
            if attrs.get("gene_name") in want and attrs.get("gene_biotype") == "protein_coding":
                canon[attrs["gene_name"]] = attrs
    cds = {}
    name = None
    with gzip.open(inputs.cds_fasta("human"), "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
                cds[name] = []
            else:
                cds[name].append(line.strip())
    out = {}
    for sym in symbols:
        a = canon.get(sym)
        if a is None:
            continue
        enst = f"{a['transcript_id']}.{a['transcript_version']}"
        seq = "".join(cds[enst]).upper()
        out[sym] = {
            "gene_id": a["gene_id"],
            "assembly": "GRCh38",
            "biotype": a["transcript_biotype"],
            "transcript_id": enst,
            "cds": seq,
            "protein": seqtools.translate(seq),
        }
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--arms", default="expression_rescue_arms.tsv")
    ap.add_argument("--out", default="rescue_cds.json")
    ap.add_argument("--offline", action="store_true", help="pinned Ensembl 115 files, no REST")
    args = ap.parse_args()

    cache = Path(args.out)
    store: dict[str, dict] = json.loads(cache.read_text()) if cache.exists() else {}

    with open(args.arms) as fh:
        symbols = [r["symbol"] for r in csv.DictReader(fh, delimiter="\t")]

    if args.offline:
        store = offline_records(symbols)
        missing = [s for s in symbols if s not in store]
        if missing:
            raise SystemExit(f"no Ensembl_canonical protein-coding transcript for {missing}")
        cache.write_text(json.dumps(store, indent=1, sort_keys=True))
        print(f"{len(store)}/{len(symbols)} receptors from the pinned Ensembl 115 files")
        return 0

    for sym in symbols:
        if sym in store:
            continue
        try:
            rec = seqtools.resolve_gene(sym)
            cds = seqtools.fetch_cds(rec.ensembl_transcript_id)
        except Exception as exc:                       # noqa: BLE001 - report and continue
            print(f"  {sym:9} FAILED: {type(exc).__name__}: {exc}")
            continue
        prot = seqtools.translate(cds)
        store[sym] = {
            "gene_id": rec.ensembl_gene_id,
            "assembly": rec.assembly,
            "biotype": rec.biotype,
            "transcript_id": rec.ensembl_transcript_id,
            "cds": cds,
            "protein": prot,
        }
        flag = "" if cds[:3] == "ATG" and prot.endswith("*") else "  <-- odd ends"
        print(f"  {sym:9} {rec.ensembl_transcript_id:20} {len(cds):5} nt  {len(prot) - 1:4} aa{flag}")
        cache.write_text(json.dumps(store, indent=1, sort_keys=True))

    cache.write_text(json.dumps(store, indent=1, sort_keys=True))
    print(f"\n{len(store)}/{len(symbols)} receptors cached in {cache}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
