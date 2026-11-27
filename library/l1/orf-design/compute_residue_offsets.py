"""Where the ordered protein sits inside the AlphaFold model, for `contained` members.

26 L1 members are ``proteoform_match = contained``: the ordered protein is a contiguous
subsequence of the AFDB model, so the model is usable once residue numbering is shifted.
``residue_offset`` exists for exactly that -- ``construct residue i == model residue
i + offset`` -- but shipped as 0 on every row because nobody computed it. The one member
whose true offset happens to be 0 (HTR7, a clean C-terminal truncation) was therefore the
only one rendered; the other 25 fell through to "no verified structure" despite having a
perfectly good model.

Offsets are measured, not inferred from the length difference: a 47-residue difference
says nothing about which end it is on. Each member's designed protein is located inside
the model's ``uniprotSequence`` with ``str.find``, and anything not found contiguously is
reported rather than given a guessed offset.

    uv run --locked python compute_residue_offsets.py --out residue_offsets.tsv
"""

from __future__ import annotations

import argparse
import csv
import json
import time
import urllib.error
import urllib.request
from pathlib import Path

from deliverome_analysis.cloning import library_design as ld

HERE = Path(__file__).parent
MANIFEST = HERE / "candidate_manifest_tranches.tsv"
PAYLOAD = HERE / "library_structure_payload.tsv"
AFDB = "https://alphafold.ebi.ac.uk/api/prediction/{acc}"


def model_sequence(acc: str, *, attempts: int = 4) -> str | None:
    """The modelled sequence AFDB actually served, not UniProt's current record."""
    for i in range(attempts):
        try:
            with urllib.request.urlopen(AFDB.format(acc=acc), timeout=40) as r:
                return json.load(r)[0]["uniprotSequence"]
        except urllib.error.HTTPError, urllib.error.URLError, KeyError, IndexError:
            if i == attempts - 1:
                return None
            time.sleep(2 * (i + 1))
    return None


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=HERE / "residue_offsets.tsv")
    ap.add_argument("--payload", type=Path, default=PAYLOAD)
    ap.add_argument(
        "--refresh",
        action="store_true",
        help="query the live AlphaFold DB API instead of reading the archived PDBs",
    )
    args = ap.parse_args(argv)

    from check_proteoform_match import _archived_afdb_index, archived_af_sequence

    archive = {} if args.refresh else _archived_afdb_index()

    members: dict[str, dict] = {}
    for library in ld.LIBRARIES:
        members |= ld.designed_proteins(MANIFEST, library=library)

    rows, unresolved = [], []
    with args.payload.open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["proteoform_match"] != "contained":
                continue
            cid = r["candidate_id"]
            ordered = members[cid]["sequence"]
            model = (
                model_sequence(r["uniprot_acc"])
                if args.refresh
                else archived_af_sequence(r["uniprot_acc"], archive)
            )
            if model is None:
                unresolved.append((cid, "AFDB fetch failed"))
                continue
            idx = model.find(ordered)
            if idx < 0:
                # Not contiguous: the `contained` verdict itself is wrong for this one.
                unresolved.append((cid, "not a contiguous subsequence of the model"))
                continue
            rows.append(
                {
                    "candidate_id": cid,
                    "gene_symbol": r["gene_symbol"],
                    "uniprot_acc": r["uniprot_acc"],
                    "ordered_length_aa": len(ordered),
                    "model_length_aa": len(model),
                    "residue_offset": idx,
                    "placement": "n_terminal_start" if idx == 0 else "interior",
                    "display_first_model_residue": idx + 1,
                    "display_last_model_residue": idx + len(ordered),
                }
            )
            time.sleep(0.4)

    with args.out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    nonzero = [r for r in rows if r["residue_offset"]]
    print(f"{len(rows)} contained members located; {len(nonzero)} need a non-zero offset")
    print(
        f"offset range: {min(r['residue_offset'] for r in rows)}"
        f"..{max(r['residue_offset'] for r in rows)}"
    )
    for cid, why in unresolved:
        print(f"  UNRESOLVED {cid}: {why}")
    print(f"wrote {args.out.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
