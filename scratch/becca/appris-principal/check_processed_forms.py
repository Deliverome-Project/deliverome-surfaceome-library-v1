"""Flag library members whose full-length ORF is not the species on the cell.

Some surface proteins are proteolytically processed before or at the membrane:
furin cuts the receptor phosphatases into an E- and a P-subunit, NOTCH is cut at
S1 into a heterodimer, adhesion GPCRs autoproteolyse at the GPS, mucins cut at
the SEA. For these, expressing the full-length ORF does not give you the
full-length protein on the surface -- and in the clearest case, polycystin-1,
the uncleaved form never reaches the surface at all.

That matters twice over for a library. A long ORF may be buying packaging budget
for a stretch of sequence that gets removed anyway, and the construct that looks
like a truncation is often the more faithful one.

Signal comes from UniProt: more than one Chain feature means the entry is
processed into separate chains, and a Site feature describing a cleavage names
the enzyme. Pulled in bulk rather than per gene.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import re
import urllib.parse
import urllib.request
from pathlib import Path

STREAM = "https://rest.uniprot.org/uniprotkb/stream"
FIELDS = "accession,gene_primary,length,ft_chain,ft_site,ft_propep"


def fetch_bulk() -> list[dict]:
    q = "(organism_id:9606) AND (reviewed:true)"
    url = f"{STREAM}?query={urllib.parse.quote(q)}&fields={FIELDS}&format=tsv"
    req = urllib.request.Request(url, headers={"User-Agent": "deliverome-analysis (research)"})
    with urllib.request.urlopen(req, timeout=600) as fh:
        text = fh.read().decode()
    return list(csv.DictReader(io.StringIO(text), delimiter="\t"))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cache", default="uniprot_processing.tsv")
    ap.add_argument("--library", default="../tag-site-architecture/library_members.json")
    ap.add_argument("--out", default="library_processed_forms.tsv")
    args = ap.parse_args()

    cache = Path(args.cache)
    if cache.exists():
        rows = list(csv.DictReader(cache.open(), delimiter="\t"))
        print(f"  {len(rows)} UniProt entries from cache")
    else:
        rows = fetch_bulk()
        with cache.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(rows)
        print(f"  fetched {len(rows)} reviewed human entries")

    members = set(json.loads(Path(args.library).read_text()))
    gene_col = next(c for c in rows[0] if "Gene" in c)
    chain_col = next(c for c in rows[0] if "Chain" in c)
    site_col = next(c for c in rows[0] if "Site" in c)

    out = []
    for r in rows:
        g = (r.get(gene_col) or "").strip()
        if g not in members:
            continue
        chains = r[chain_col].count("CHAIN ")
        sites = re.findall(r'/note="([^"]*[Cc]leav[^"]*)"', r[site_col] or "")
        if chains <= 1 and not sites:
            continue
        enzymes = sorted({s.split("by ")[-1] if "by " in s else "unspecified" for s in sites})
        out.append({
            "gene": g, "accession": r["Entry"], "length_aa": r["Length"],
            "chains": chains, "cleavage_sites": len(sites),
            "cleaved_by": "; ".join(enzymes) or "",
            "full_length_is_surface_species": "no" if (chains > 1 or sites) else "yes",
        })

    out.sort(key=lambda r: -int(r["length_aa"]))
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out)

    from collections import Counter
    print(f"\n  {len(out)} of {len(members)} library members are proteolytically processed")
    print(f"  ({len(out)/len(members):.0%} of the library)")
    big = [r for r in out if int(r["length_aa"]) > 1000]
    print(f"  of those, {len(big)} are over 1,000 aa - where the wasted budget is largest\n")
    print(f"  {'gene':10} {'aa':>5} {'chains':>6}  cleaved by")
    for r in out[:12]:
        print(f"  {r['gene']:10} {r['length_aa']:>5} {r['chains']:>6}  {r['cleaved_by'][:44]}")
    enz = Counter(e for r in out for e in r["cleaved_by"].split("; ") if e)
    print(f"\n  most common annotated protease: {dict(enz.most_common(5))}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
