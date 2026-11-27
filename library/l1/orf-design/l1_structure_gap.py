"""Which Library-1 members have no AlphaFold model, i.e. which ones we would have to fold.

Every tag-site lane and EpicTope read pLDDT, RSA (freesasa) and secondary structure
(pydssp) off an AlphaFold model, so a member without one cannot be scored by any of them.
This enumerates all 2,725 L1 candidates, resolves each to the accession AlphaFold DB is
keyed by, and records whether a model exists.

Resolution differs by set, which is the whole point: reference and rescue candidates are
gene symbols (resolved offline through the surfaceome catalog), isoforms are either a
UniProt isoform id or a bare Ensembl transcript, and orthologs carry their own accession.
An Ensembl-only identifier has to go through UniProt first, and those are the ones that
tend to land on TrEMBL accessions newer than AlphaFold DB.

    uv run --locked python l1_structure_gap.py            # read the cache
    uv run --locked python l1_structure_gap.py --refresh  # re-measure (network, ~minutes)
    uv run --locked python l1_structure_gap.py --library l2 --refresh \
        --out l2_structure_gap.tsv                        # the same question for L2
"""

from __future__ import annotations

import argparse
import csv
import json
import urllib.error
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from deliverome_analysis.surfaceome import surfaceome

HERE = Path(__file__).parent
MANIFEST = HERE / "candidate_manifest_tranches.tsv"
OUT = HERE / "l1_structure_gap.tsv"
WORKERS = 12


def _json(url: str, timeout: int = 25):
    return json.load(urllib.request.urlopen(url, timeout=timeout))


def uniprot_for_transcript(enst: str) -> str | None:
    q = urllib.parse.quote(f"xref:{enst}")
    try:
        hits = _json(
            f"https://rest.uniprot.org/uniprotkb/search?query={q}&fields=accession"
            "&format=json&size=1"
        ).get("results", [])
    except Exception:
        return None
    return hits[0]["primaryAccession"] if hits else None


def has_model(acc: str) -> bool:
    try:
        urllib.request.urlopen(
            f"https://alphafold.ebi.ac.uk/api/prediction/{acc}", timeout=25
        ).read()
        return True
    except urllib.error.HTTPError, urllib.error.URLError, TimeoutError, OSError:
        return False


def members(library: str = "l1") -> list[dict]:
    """One row per candidate in ``library`` with the identifier AlphaFold DB is keyed by."""
    cat = surfaceome.load_catalog()
    sym2acc = dict(zip(cat["symbol"], cat["uniprot_accession"], strict=True))
    rows = []
    for r in csv.DictReader(MANIFEST.open(), delimiter="\t"):
        if r["library"] != library:
            continue
        parts = r["key"].split(":")
        setname = r["set"]
        if setname in ("library_reference", "expression_rescue", "giant_truncation"):
            symbol, ident = parts[-1], sym2acc.get(parts[-1])
        elif setname == "isoform":
            symbol, ident = parts[1], parts[2]
        else:  # ortholog: ort:species:SYMBOL:acc
            symbol, ident = parts[2], (parts[3] if len(parts) > 3 else None)
        rows.append(
            {
                "candidate_id": r["candidate_id"],
                "set": setname,
                "species": r["species"],
                "symbol": symbol,
                "identifier": ident or "",
                "needs_uniprot_lookup": bool(ident and ident.startswith("ENST")),
            }
        )
    return rows


def measure(rows: list[dict]) -> list[dict]:
    # Ensembl-only identifiers need a UniProt hop before AlphaFold DB can be asked.
    todo = sorted({r["identifier"] for r in rows if r["needs_uniprot_lookup"]})
    print(f"resolving {len(todo)} Ensembl identifiers through UniProt...")
    with ThreadPoolExecutor(WORKERS) as pool:
        mapped = dict(zip(todo, pool.map(uniprot_for_transcript, todo), strict=True))
    for r in rows:
        r["accession"] = mapped.get(r["identifier"], r["identifier"]) or ""

    accs = sorted({r["accession"] for r in rows if r["accession"]})
    print(f"checking AlphaFold DB for {len(accs)} distinct accessions...")
    with ThreadPoolExecutor(WORKERS) as pool:
        present = dict(zip(accs, pool.map(has_model, accs), strict=True))
    for r in rows:
        r["has_alphafold"] = bool(present.get(r["accession"], False))
    return rows


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--refresh", action="store_true")
    ap.add_argument("--library", default="l1", help="tranche to measure: l1 or l2")
    ap.add_argument("--out", type=Path, default=None, help="where to cache the result")
    ap.add_argument(
        "--new-only",
        action="store_true",
        help="measure only candidates the cached table lacks and append them; recorded rows "
        "are not re-queried, so a library addition does not silently re-measure the rest",
    )
    args = ap.parse_args(argv)

    out = args.out or (HERE / f"{args.library}_structure_gap.tsv")
    if args.new_only and out.exists():
        have = list(csv.DictReader(out.open(), delimiter="\t"))
        known = {r["candidate_id"] for r in have}
        todo = [r for r in members(args.library) if r["candidate_id"] not in known]
        print(f"{len(todo)} candidate(s) not yet measured")
        rows = have + measure(todo) if todo else have
    elif args.refresh or not out.exists():
        rows = measure(members(args.library))
    else:
        rows = None
    if rows is not None:
        cols = [
            "candidate_id",
            "set",
            "species",
            "symbol",
            "identifier",
            "needs_uniprot_lookup",
            "accession",
            "has_alphafold",
        ]
        with out.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(rows)
    rows = list(csv.DictReader(out.open(), delimiter="\t"))

    import collections

    tally = collections.defaultdict(collections.Counter)
    for r in rows:
        ok = r["has_alphafold"] in ("True", "true", "1")
        tally[r["set"]]["n"] += 1
        tally[r["set"]]["have"] += ok
        if not r["accession"]:
            tally[r["set"]]["unresolved"] += 1
    print(f"\n{'set':<20} {'n':>6} {'has model':>11} {'MISSING':>9} {'unresolved':>11}")
    tot = miss = 0
    for s, c in sorted(tally.items()):
        tot += c["n"]
        miss += c["n"] - c["have"]
        print(f"{s:<20} {c['n']:>6} {c['have']:>11} {c['n'] - c['have']:>9} {c['unresolved']:>11}")
    print(f"{args.library.upper() + ' TOTAL':<20} {tot:>6} {tot - miss:>11} {miss:>9}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
