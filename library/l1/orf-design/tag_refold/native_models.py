"""Index every existing native model by the SHA-256 of the protein it models.

Member ids were renumbered twice on 2026-11-12 (PR #283), so a model is looked up by
sequence, never by candidate_id. Two sources:

* ``library_structure_payload.tsv`` on main: each chosen row's model (AFDB URL or our
  ColabFold key) belongs to the designed protein of that candidate_id under main's
  manifest, provided the lengths agree.
* every ColabFold PDB under ``processed://…/predicted-structures/colabfold/``: its own
  ATOM sequence, read from the file (this includes the 66 ortholog models
  l1-reviewer-package folded on 2026-11-10/11).

A native with no model here is folded untagged by this folder (``n.<sha16>``).

    uv run --locked python library/l1/orf-design/tag_refold/native_models.py
"""

from __future__ import annotations

import csv
import hashlib
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
EXP = HERE.parent
OUT = EXP / "tag_refold_native_models.tsv"
PREFIX = "sprint4/s4e9-orf-design-fg-panel/predicted-structures/colabfold/"
THREE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
    "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
    "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}  # fmt: skip


def sha(s: str) -> str:
    return hashlib.sha256(s.encode()).hexdigest()


def pdb_sequence(path: Path) -> str:
    seq, last = [], None
    for line in path.read_text().splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            key = (line[21], line[22:27])
            if key != last:
                seq.append(THREE.get(line[17:20], "X"))
                last = key
    return "".join(seq)


#: main before #283 merged: the payload's candidate_ids are that manifest's numbering.
MAIN_COMMIT = "5b133d0b"


def _at_main(rel: str) -> Path:
    import tempfile

    from deliverome_analysis.provenance import read_pinned_text

    out = Path(tempfile.mkdtemp()) / Path(rel).name
    out.write_text(
        read_pinned_text(
            MAIN_COMMIT, f"library/l1/orf-design/{rel}", repo=EXP.parents[2]
        )
    )
    return out


def build() -> list[dict[str, str]]:
    from deliverome_analysis import library_design as ld
    from deliverome_analysis.s3 import processed

    rows: dict[str, dict[str, str]] = {}
    designed = ld.designed_proteins(
        _at_main("candidate_manifest_tranches.tsv"), library=ld.LIBRARY_L1
    )
    with _at_main("library_structure_payload.tsv").open() as fh:
        for p in csv.DictReader(fh, delimiter="\t"):
            d = designed.get(p["candidate_id"])
            if p["is_chosen"] != "1" or not d or int(p["ordered_length_aa"]) != len(d["sequence"]):
                continue
            if p["source"] == "alphafold_db":
                ref = p["model_url"]
            else:
                ref = "processed://" + PREFIX + p["r2_key"].removeprefix("library-colabfold/")
            rows.setdefault(
                sha(d["sequence"]),
                {"native_sha256": sha(d["sequence"]), "source": p["source"], "ref": ref,
                 "match": p["proteoform_match"], "via": f"payload {p['candidate_id']} (main)"},
            )  # fmt: skip
    keys = [k for k in processed.list_keys(PREFIX) if k.endswith(".pdb")]
    for k in sorted(keys):
        s = pdb_sequence(processed.download(k))
        rows.setdefault(
            sha(s),
            {"native_sha256": sha(s), "source": "colabfold", "ref": "processed://" + k,
             "match": "identical", "via": "ColabFold model sequence"},
        )  # fmt: skip
    print(
        f"{len(rows)} native sequences with a model ({len(keys)} ColabFold PDBs read)",
        file=sys.stderr,
    )
    return sorted(rows.values(), key=lambda r: r["native_sha256"])


def main() -> None:
    rows = build()
    with OUT.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {OUT.name}: {len(rows)} rows")


if __name__ == "__main__":
    main()
