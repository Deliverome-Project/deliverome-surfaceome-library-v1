"""Sequences for the L1 members whose AlphaFold model is a DIFFERENT protein.

`check_proteoform_match.py` finds them: 200 members where the modelled sequence is not the
one the library would synthesise, so any tag site derived from that model is numbered
against a sequence the construct does not contain. Unlike the missing-model set, nothing
needs fetching — the designed coding sequence is already on S3 in
``library-candidates/candidate_cds.json`` and is simply translated here.

**Keyed by candidate_id, not by accession.** The accession's model is precisely the wrong
one, so filing our prediction under it would put the correction in the same slot as the
thing being corrected, and a later reader could not tell them apart.

    uv run --locked python build_mismatch_fold_inputs.py
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

from deliverome_analysis import library_design as ld
from deliverome_analysis import seqtools
from deliverome_analysis.s3 import processed

HERE = Path(__file__).parent
MATCH = HERE / "proteoform_match.tsv"
TRANCHES = HERE / "candidate_manifest_tranches.tsv"
CDS_KEY = ld.CDS_KEY
FASTA = HERE / "fold_inputs_mismatched.fasta"
MANIFEST = HERE / "fold_inputs_mismatched.tsv"


def main() -> int:
    cds = json.load(open(processed.download(CDS_KEY)))
    man = {r["candidate_id"]: r for r in csv.DictReader(TRANCHES.open(), delimiter="\t")}

    rows, no_cds = [], []
    for r in csv.DictReader(MATCH.open(), delimiter="\t"):
        if r["verdict"] != "different":
            continue
        meta = man.get(r["candidate_id"], {})
        entry = cds.get(meta.get("key", ""))
        if not entry:
            no_cds.append(r["candidate_id"])
            continue
        seq = seqtools.translate(entry["cds"]).rstrip("*")
        rows.append(
            {
                # The fold job (and its structure file) is named by the design id, which
                # does not move when members are renumbered.
                "name": meta.get("design_id") or r["candidate_id"],
                "symbol": r["symbol"],
                "species": meta.get("species", ""),
                "set": r["set"],
                "wrong_model_accession": r["accession"],
                "length": len(seq),
                "sequence": seq,
            }
        )
    if no_cds:
        print(f"  no designed CDS for {len(no_cds)}: {no_cds[:5]}")

    rows.sort(key=lambda r: r["name"])
    FASTA.write_text(
        "".join(
            f">{r['name']}|{r['symbol']}|{r['species']}|{r['set']}\n"
            + "\n".join(r["sequence"][i : i + 60] for i in range(0, len(r["sequence"]), 60))
            + "\n"
            for r in rows
        )
    )
    cols = ["name", "symbol", "species", "set", "wrong_model_accession", "length"]
    with MANIFEST.open("w", newline="") as fh:
        w = csv.DictWriter(
            fh, fieldnames=cols, delimiter="\t", lineterminator="\n", extrasaction="ignore"
        )
        w.writeheader()
        w.writerows(rows)

    lens = [r["length"] for r in rows]
    print(
        f"wrote {FASTA.name}: {len(rows)} sequences, {min(lens)}-{max(lens)} aa "
        f"(median {sorted(lens)[len(lens) // 2]}), {sum(lens):,} residues total"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
