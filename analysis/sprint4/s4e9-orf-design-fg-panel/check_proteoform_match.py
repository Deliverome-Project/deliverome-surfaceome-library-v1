"""Is the AlphaFold model we would score the SAME PROTEIN as the ORF we would build?

`l1_structure_gap.py` answers "is there a model", which is the weaker question. A model for
the wrong proteoform is worse than no model: every coordinate the lanes emit is numbered
against a sequence the construct does not contain, and nothing downstream flags it.

The comparison is residue-exact. The designed coding sequence comes from
``library-candidates/candidate_cds.json`` on S3 — what the library would actually
synthesise — and is translated and compared against the AlphaFold entry's
``uniprotSequence``. An earlier version of this script compared *lengths*, which is a
necessary condition and not a sufficient one: it cannot see a same-length substitution.

Three outcomes, because the middle one is salvageable rather than broken:

``identical``
    the model is of exactly our protein.
``contained``
    our protein is a contiguous subsequence of the modelled one, so the model is usable
    with a constant residue offset, which is reported.
``different``
    neither; the model cannot be numbered onto the construct and the ORF needs folding.

    uv run --locked python check_proteoform_match.py
    uv run --locked python check_proteoform_match.py --set library_reference
"""

from __future__ import annotations

import argparse
import csv
import json
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from deliverome_analysis import library_design as ld
from deliverome_analysis import seqtools
from deliverome_analysis.s3 import processed

HERE = Path(__file__).parent
GAP = HERE / "l1_structure_gap.tsv"
MANIFEST = HERE / "candidate_manifest_tranches.tsv"
CDS_KEY = ld.CDS_KEY
OUT = HERE / "proteoform_match.tsv"
WORKERS = 12


def af_sequence(acc: str) -> str | None:
    try:
        data = json.load(
            urllib.request.urlopen(f"https://alphafold.ebi.ac.uk/api/prediction/{acc}", timeout=30)
        )
    except Exception:
        return None
    return (data[0].get("uniprotSequence") or None) if data else None


def classify(ours: str, modelled: str) -> tuple[str, int]:
    if ours == modelled:
        return "identical", 0
    at = modelled.find(ours)
    if at >= 0:
        return "contained", at  # construct residue i is modelled residue i + at
    return "different", 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--set", dest="which", default=None)
    ap.add_argument("--out", type=Path, default=OUT)
    ap.add_argument(
        "--new-only",
        action="store_true",
        help="compare only candidates the existing table lacks, and keep its recorded rows",
    )
    args = ap.parse_args(argv)
    kept = []
    if args.new_only and args.out.exists():
        kept = list(csv.DictReader(args.out.open(), delimiter="\t"))
    known = {r["candidate_id"] for r in kept}

    cds = json.load(open(processed.download(CDS_KEY)))
    man = {r["candidate_id"]: r for r in csv.DictReader(MANIFEST.open(), delimiter="\t")}
    rows = [
        r
        for r in csv.DictReader(GAP.open(), delimiter="\t")
        if r["has_alphafold"] in ("True", "true", "1")
        and (args.which is None or r["set"] == args.which)
        and r["candidate_id"] not in known
    ]
    print(f"checking {len(rows)} L1 candidate(s) that have a model, residue-exact")

    with ThreadPoolExecutor(WORKERS) as pool:
        modelled = list(pool.map(af_sequence, [r["accession"] for r in rows]))

    out_rows, skipped = [], 0
    for r, model_seq in zip(rows, modelled, strict=True):
        key = man.get(r["candidate_id"], {}).get("key", "")
        entry = cds.get(key)
        if not entry or not model_seq:
            skipped += 1
            continue
        ours = seqtools.translate(entry["cds"]).rstrip("*")
        verdict, offset = classify(ours, model_seq)
        out_rows.append(
            {
                "candidate_id": r["candidate_id"],
                "set": r["set"],
                "symbol": r["symbol"],
                "accession": r["accession"],
                "orf_aa": len(ours),
                "model_aa": len(model_seq),
                "verdict": verdict,
                "offset": offset,
            }
        )

    out_rows = kept + out_rows
    with args.out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out_rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(sorted(out_rows, key=lambda r: (r["verdict"], r["set"], r["symbol"])))

    n = len(out_rows)
    print(f"\n{n} compared ({skipped} skipped: no designed CDS or no model sequence)")
    for v in ("identical", "contained", "different"):
        sub = [r for r in out_rows if r["verdict"] == v]
        print(f"  {v:<11} {len(sub):>5}  ({len(sub) / n:.1%})")
        by = {}
        for r in sub:
            by[r["set"]] = by.get(r["set"], 0) + 1
        if sub:
            print(f"              by set: {by}")
    worst = [r for r in out_rows if r["verdict"] == "different"]
    for r in sorted(worst, key=lambda r: -abs(int(r["model_aa"]) - int(r["orf_aa"])))[:10]:
        print(
            f"     {r['symbol']:<12} {r['accession']:<8} model {r['model_aa']:>5} "
            f"vs ORF {r['orf_aa']:>5}"
        )
    print(f"\nwrote {args.out.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
