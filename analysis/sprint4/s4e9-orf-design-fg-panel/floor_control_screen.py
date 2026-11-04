"""Mechanistic screen for a low-internalization ("floor") control that is easy to tag.

MS4A1/CD20 is the panel's floor control, but it is 4TM with both termini cytoplasmic,
so it needs a loop insertion rather than the terminal cassette the other seven use.
This looks for a replacement or a second floor arm that takes a terminal tag.

Two filters the atlas `intern` SeqGrade alone does not give:

1. **GPI-anchored proteins are excluded outright.** They are the defining cargo of the
   CLIC/GEEC pathway — GPI-AP nanoclusters "act as sorting signals for the selective
   endocytosis of GPI-APs" (Chadda 2007, doi:10.1111/j.1600-0854.2007.00565.x; GRAF1
   marks a pathway internalizing "GPI-linked proteins", Lundmark 2008,
   doi:10.1016/j.cub.2008.10.044). A GPI anchor means clathrin-independent uptake, not
   no uptake — so THY1/CD59/LY6E/NT5E are not floor controls however low the grade says.

2. **Cytoplasmic sorting motifs.** Constitutive clathrin-mediated endocytosis of a
   single-pass type I protein is driven by an adaptor-binding motif in the cytoplasmic
   tail. A short tail with none of them has no CME route — that is the mechanistic
   basis for a floor control, and unlike the SeqGrade it is computable. s2e5 showed the
   SeqGrade does not predict measured constitutive internalization (Mann-Whitney
   P = 0.88 / 0.21), so it is used only to pick the starting pool, never to rank.

Motifs scanned (cytoplasmic tail only):
  YXXO      tyrosine-based, AP-2 mu2         Y-x-x-[FILMV]
  [DE]XXXLL dileucine, AP-2 sigma2/alpha     [DE]-x-x-x-[LI]-[LI]
  NPXY      NPXY, Dab2/ARH                   N-P-x-Y

A hit does not prove internalization and a clean tail does not prove its absence —
this ranks candidates for a literature check, it does not replace one.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path

from deliverome_analysis.s3 import external

CATALOG = "accessible-surfaceome/api-catalog/2026-08-26/catalog.json"
UNIPROT = "https://rest.uniprot.org/uniprotkb/search"
AGENT = "deliverome-analysis floor-control-screen"

MOTIFS = {
    "YXXO": re.compile(r"Y..[FILMV]"),
    "DEXXXLL": re.compile(r"[DE]...[LI][LI]"),
    "NPXY": re.compile(r"NP.Y"),
}


def uniprot_record(acc: str) -> dict | None:
    url = (
        UNIPROT
        + "?"
        + urllib.parse.urlencode(
            {
                "query": f"accession:{acc}",
                "fields": "accession,length,sequence,ft_topo_dom,ft_transmem,ft_signal,ft_lipid",
                "format": "json",
                "size": "1",
            }
        )
    )
    with urllib.request.urlopen(
        urllib.request.Request(url, headers={"User-Agent": AGENT})
    ) as response:
        hits = json.load(response).get("results", [])
    return hits[0] if hits else None


def cytoplasmic_tail(record: dict) -> tuple[str, int, int] | None:
    """(tail sequence, start, end) of the C-terminal cytoplasmic topological domain."""
    seq = record["sequence"]["value"]
    topo = [f for f in record.get("features", []) if f["type"] == "Topological domain"]
    if not topo or "Cytoplasmic" not in topo[-1]["description"]:
        return None
    start = topo[-1]["location"]["start"]["value"]
    end = topo[-1]["location"]["end"]["value"]
    if end != len(seq):
        return None
    return seq[start - 1 : end], start, end


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--out", type=Path, default=Path(__file__).parent / "floor_control_candidates.tsv"
    )
    args = parser.parse_args(argv)

    catalog = json.load(open(external.download(CATALOG)))
    rows = {r["symbol"]: r for r in catalog["rows"]}

    pool = []
    for symbol, row in rows.items():
        facets = row.get("ddf") or {}
        if row.get("intern") not in ("very_low", "low"):
            continue
        if facets.get("subcategory") != "single_pass_T1":  # GPI excluded — GEEC cargo
            continue
        if facets.get("surface_accessibility") != "high":
            continue
        if facets.get("evidence_grade") != "direct_multi_method":
            continue
        if facets.get("confidence") != "high":
            continue
        pool.append((symbol, row))
    # TFRC is the scan's own positive control: its YTRF motif must be found.
    pool.append(("TFRC", rows["TFRC"]))
    print(f"{len(pool) - 1} single-pass T1 candidates + TFRC as the scan's control\n")

    out = []
    for symbol, row in sorted(pool):
        acc = (
            (row.get("uniprot") or [None])[0]
            if isinstance(row.get("uniprot"), list)
            else row.get("uniprot")
        )
        record = uniprot_record(acc) if acc else None
        time.sleep(0.15)
        if not record:
            print(f"{symbol:10s} no UniProt record for {acc}")
            continue
        tail = cytoplasmic_tail(record)
        if not tail:
            print(f"{symbol:10s} no C-terminal cytoplasmic domain — skipped")
            continue
        seq, start, _end = tail
        hits = {name: [m.start() + start for m in rx.finditer(seq)] for name, rx in MOTIFS.items()}
        n_hits = sum(len(v) for v in hits.values())
        facets = row.get("ddf") or {}
        out.append(
            {
                "gene": symbol,
                "uniprot": acc,
                "length_aa": record["sequence"]["length"],
                "orf_nt": record["sequence"]["length"] * 3 + 3,
                "intern_grade": row.get("intern"),
                "tail_start": start,
                "tail_len": len(seq),
                "tail_seq": seq,
                "n_motifs": n_hits,
                "motifs": ";".join(f"{k}@{','.join(map(str, v))}" for k, v in hits.items() if v)
                or "none",
                "expression": f"{facets.get('expression_level')}/{facets.get('expression_breadth')}",
            }
        )

    out.sort(key=lambda r: (r["n_motifs"], r["tail_len"], r["orf_nt"]))
    cols = [
        "gene",
        "uniprot",
        "length_aa",
        "orf_nt",
        "intern_grade",
        "tail_start",
        "tail_len",
        "n_motifs",
        "motifs",
        "expression",
        "tail_seq",
    ]
    with args.out.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=cols, delimiter="\t", extrasaction="ignore", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(out)

    print(
        f"\n{'gene':10s} {'aa':>5s} {'ORF nt':>7s} {'tail':>5s} {'motifs':>6s}  "
        f"{'positions':28s} expression"
    )
    for r in out:
        print(
            f"{r['gene']:10s} {r['length_aa']:5d} {r['orf_nt']:7d} {r['tail_len']:5d} "
            f"{r['n_motifs']:6d}  {r['motifs'][:28]:28s} {r['expression']}"
        )
    print(f"\nWrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
