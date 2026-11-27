"""Drop every domesticated candidate ORF into the p0039 backbone.

This is the pre-tag intermediate. The finished s4e9 insert runs Kozak -> ORF ->
tag cassette -> barcode -> reverse primer site, but tag placement is unresolved
at library scale and barcodes are deliberately not reserved yet, so those two
blocks are left out. What is assembled here is the part that is settled: the
backbone's cloning slot replaced by the Kozak block and the domesticated ORF.

Re-running after tag sites land replaces the insert, not the backbone, so
nothing here is wasted.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from pathlib import Path

from Bio import SeqIO

from deliverome_analysis.cloning import plasmid_qc as qc
from deliverome_analysis.cloning import seqtools

KOZAK = qc.KOZAK
ESP3I_5 = qc.ESP3I_5
REV_PRIMER = qc.REV_PRIMER
CS1 = qc.CS1
ESP3I_3 = qc.ESP3I_3
SLOT_START = qc.BACKBONE_SLOT_START
SLOT_LEN = qc.BACKBONE_SLOT_LEN
ENZYMES = list(qc.DOMESTICATION_ENZYMES)
SET_NAME = {"lib": "library_reference", "iso": "isoform", "ort": "ortholog",
            "rescue": "expression_rescue", "giant": "giant_truncation"}


def assemble_one(bb: str, bb_sites: list, cid: str, key: str, v: dict) -> tuple[str, dict, bool]:
    """(vector, manifest row, matches the expected site pattern) for one candidate."""
    orf = v["construct"]
    # The barcode and tag cassette are deferred; everything structural
    # around them is not, so the wrapper goes on now.
    insert = ESP3I_5 + KOZAK + orf + REV_PRIMER + CS1 + ESP3I_3
    vector = bb[:SLOT_START] + insert + bb[SLOT_START + SLOT_LEN:]
    sites = seqtools.find_sites(vector, ENZYMES)
    esp = [x for x in sites if x[0] == "Esp3I"]
    # two Esp3I and nothing else beyond the backbone's own BsaI
    clean = len(esp) == 2 and len(sites) == len(bb_sites) + 2
    row = {
        "candidate_id": cid,
        "key": key,
        "set": SET_NAME[key.split(":")[0]],
        "species": v["species"],
        "transcript": v["transcript"],
        "route": v.get("route", ""),
        "variant_residues": v.get("variant_residues", ""),
        "orf_nt": len(orf),
        "vector_bp": len(vector),
        "sites_removed": "" if v["sites_before"] is None else v["sites_before"],
        "edits": "" if v["edits"] is None else v["edits"],
        "vector_sites": len(sites),
        "esp3i_sites": len(esp),
        "insert_nt": len(insert),
        "has_igk_leader": key.startswith("rescue"),
        "tag_cassette": "deferred",
        "barcode": "not reserved",
    }
    return vector, row, clean


def fasta_record(cid: str, key: str, vector: str) -> str:
    lines = [vector[j:j + 70] for j in range(0, len(vector), 70)]
    return f">{cid} {key} {len(vector)}bp set={SET_NAME[key.split(':')[0]]}\n" + "\n".join(lines) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--backbone", default=None,
                    help="defaults to the pinned p0039 backbone on S3")
    ap.add_argument("--constructs", default="candidate_constructs.json")
    ap.add_argument("--fasta", default="candidate_vectors.fasta.gz")
    ap.add_argument("--manifest", default="candidate_manifest.tsv")
    args = ap.parse_args()

    import candidate_inputs as inputs
    bb_path = args.backbone or inputs.backbone()
    bb = str(SeqIO.read(bb_path, "genbank").seq).upper()
    if len(bb) != 7672:
        raise SystemExit(f"backbone is {len(bb)} bp, expected 7672")
    bb_sites = seqtools.find_sites(bb, ENZYMES)
    print(f"  backbone {len(bb)} bp, slot {SLOT_LEN} nt at {SLOT_START}, "
          f"backbone sites {len(bb_sites)}")

    d = json.loads(Path(args.constructs).read_text())
    rank = {"lib": 0, "iso": 1, "ort": 2, "rescue": 3, "giant": 4}
    keys = sorted(d, key=lambda k: (rank[k.split(":")[0]], k))

    rows, dirty = [], 0
    with gzip.open(args.fasta, "wt") as fh:
        for i, key in enumerate(keys, 1):
            cid = f"l1-p{i:04d}"
            vector, row, clean = assemble_one(bb, bb_sites, cid, key, d[key])
            dirty += not clean
            fh.write(fasta_record(cid, key, vector))
            rows.append(row)

    with open(args.manifest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    from collections import Counter
    c = Counter(r["set"] for r in rows)
    print(f"\n  assembled {len(rows)} vectors: {rows[0]['candidate_id']} .. {rows[-1]['candidate_id']}")
    for k, n in c.most_common():
        ids = [r["candidate_id"] for r in rows if r["set"] == k]
        print(f"    {k:20} {n:5}   {ids[0]} .. {ids[-1]}")
    print(f"  vector size {min(r['vector_bp'] for r in rows)}-{max(r['vector_bp'] for r in rows)} bp")
    print(f"  vectors NOT matching the expected site pattern: {dirty}")
    ins = [r["insert_nt"] for r in rows]
    import statistics as stat
    print(f"  insert (no barcode/cassette yet): median {stat.median(ins)/1000:.2f} kb, "
          f"{min(ins)}-{max(ins)} nt")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
