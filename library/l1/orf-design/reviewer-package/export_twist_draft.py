"""Write the L1 draft library as full p0039 plasmids, in plate order, for a vendor quote.

Each member is the p0039 backbone with its cloning slot replaced by the Esp3I-wrapped
insert the library pipeline assembles (``assemble_candidates.assemble_one``): Kozak, the
designed CDS, the reverse-primer site and 10x CS1. The tag cassette and the 20-nt barcode
are not yet assigned at library scale, so they are absent; a final order will add them.

Order and names come from ``l1_order_plan.tsv`` (plates sorted by insert length). Every
plasmid is checked against the tranche manifest's ``vector_bp`` and for exactly two Esp3I
sites (the insert wrapper) before it is written.

    uv run --locked python export_twist_draft.py --out-dir ~/Downloads/twist-l1-draft
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path

from Bio import SeqIO

from deliverome_analysis import library_design as ld
from deliverome_analysis import seqtools
from deliverome_analysis.s3 import processed

HERE = Path(__file__).parent
S4E9 = HERE.parent
REPO = S4E9.parents[2]
sys.path.insert(0, str(REPO / "library/l1/proteoform-selection"))
import assemble_candidates  # noqa: E402
import candidate_inputs  # noqa: E402

CONSTRUCTS = ld.CONSTRUCTS_KEY
GENE_SYNTHESIS_LIMIT = 5000  # flagged only, for the vendor to comment on


def gene_of(key: str) -> str:
    return key.split(":")[1]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args(argv)
    out = args.out_dir
    out.mkdir(parents=True, exist_ok=True)

    con = json.loads(processed.download(CONSTRUCTS, out / ".candidate_constructs.json").read_text())
    (out / ".candidate_constructs.json").unlink()
    bb = str(SeqIO.read(candidate_inputs.backbone(), "genbank").seq).upper()
    bb_sites = seqtools.find_sites(bb, assemble_candidates.ENZYMES)
    with (S4E9 / "candidate_manifest_tranches.tsv").open() as fh:
        manifest = {r["key"]: r for r in csv.DictReader(fh, delimiter="\t")}
    with (HERE / "l1_order_plan.tsv").open() as fh:
        plan = list(csv.DictReader(fh, delimiter="\t"))

    rows, fasta = [], []
    # The draft is untagged, so a member ordered as several tagged arms is one draft plasmid:
    # keep each member's first insert in plate order and name it by its design id.
    first = {}
    for p in plan:
        first.setdefault(p["design_id"], p)
    for p in first.values():
        key = p["key"]
        vec, row, clean = assemble_candidates.assemble_one(
            bb, bb_sites, p["design_id"], key, con[key]
        )
        assert clean, f"{key}: unexpected Type IIS sites"
        assert len(vec) == int(manifest[key]["vector_bp"]), key
        name = f"{p['design_id']}_{gene_of(key)}"
        fasta.append(f">{name}\n" + "\n".join(vec[i : i + 80] for i in range(0, len(vec), 80)))
        rows.append(
            {
                "name": name,
                "plate": p["plate"],
                "well": p["well"],
                "design_id": p["design_id"],
                "key": key,
                "set": row["set"],
                "insert_nt": row["insert_nt"],
                "plasmid_bp": len(vec),
                "gc_insert": round(
                    sum(
                        b in "GC" for b in vec[assemble_candidates.SLOT_START :][: row["insert_nt"]]
                    )
                    / row["insert_nt"],
                    3,
                ),
            }
        )

    (out / "l1_draft_plasmids.fasta").write_text("\n".join(fasta) + "\n")
    # the manifest carries no sequence: the FASTA is the one copy, so the two cannot disagree
    cols = list(rows[0])
    with (out / "l1_draft_plasmids.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    ins = sorted(r["insert_nt"] for r in rows)
    summary = {
        "members": len(rows),
        "plates": len({r["plate"] for r in rows}),
        "backbone": "p0039 (processed://" + candidate_inputs.BACKBONE_KEY + ")",
        "backbone_bp": len(bb),
        "insert_nt_min_median_max": [ins[0], ins[len(ins) // 2], ins[-1]],
        "plasmid_bp_min_max": [
            min(r["plasmid_bp"] for r in rows),
            max(r["plasmid_bp"] for r in rows),
        ],
        "inserts_over_5kb": sum(i > GENE_SYNTHESIS_LIMIT for i in ins),
        "total_insert_nt": sum(ins),
        "fasta_sha256": hashlib.sha256((out / "l1_draft_plasmids.fasta").read_bytes()).hexdigest(),
    }
    (out / "l1_draft_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
