"""Build pre-tag p0039 intermediates for the selected Compara116 ortholog forms.

This consumes a form-keyed CDS file; it never resumes from or overwrites historical
candidate records. Tags and barcodes are deferred, as in assemble_candidates.py.
An optional adoption member map supplies stable IDs and keys; without it, IDs are
explicitly pending and the output is for construction review only.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import sys
from pathlib import Path

from Bio import SeqIO

from deliverome_analysis import seqtools

HERE = Path(__file__).resolve().parent
AP = HERE.parents[2] / "library/l1/proteoform-selection"


def construction_tools():
    if str(AP) not in sys.path:
        sys.path.insert(0, str(AP))
    import assemble_candidates
    import domesticate_candidates

    return assemble_candidates, domesticate_candidates


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def build_forms(panel: list[dict], cds: dict, backbone: str, member_map: dict | None = None):
    """Return per-form constructs, assembly rows and vectors after strict checks."""
    assembly, domestication = construction_tools()
    if len(backbone) != 7672:
        raise ValueError(f"p0039 backbone must be 7672 bp, got {len(backbone)}")
    if not panel:
        raise ValueError("selected panel must not be empty")
    form_ids = [r["form_id"] for r in panel]
    if len(set(form_ids)) != len(form_ids):
        raise ValueError("duplicate form_id in selected panel")
    if set(cds) != set(form_ids):
        raise ValueError("CDS forms must exactly match the selected panel")
    if member_map is not None and set(member_map) != set(form_ids):
        raise ValueError("member map must exactly match the selected panel")
    if member_map is not None:
        ids = [r["candidate_id"] for r in member_map.values()]
        if len(set(ids)) != len(ids) or any(not x for x in ids):
            raise ValueError("member IDs must be nonempty and unique")
    bb_sites = seqtools.find_sites(backbone, assembly.ENZYMES)
    constructs, rows, vectors = {}, [], []
    for selected in panel:
        fid = selected["form_id"]
        entry = cds[fid]
        raw = entry["cds"].upper()
        if not raw or len(raw) % 3 or set(raw) - set("ACGT"):
            raise ValueError(f"{fid}: invalid unambiguous complete CDS")
        protein = seqtools.translate(raw)
        if not protein.endswith("*") or "*" in protein[:-1]:
            raise ValueError(f"{fid}: require one terminal stop and no internal stops")
        if sha256(protein[:-1].encode()) != selected["sequence_sha256"]:
            raise ValueError(f"{fid}: CDS does not encode the selected protein")
        if entry["species"] != selected["species"]:
            raise ValueError(f"{fid}: CDS species differs from selected panel")
        # Fresh call for every form: the legacy CLI resumes by key alone and would
        # silently keep an old construct when the CDS behind that key changes.
        made = domestication.domesticate_one({**entry, "cds": raw})
        if made is None or seqtools.translate(made["construct"]) != protein:
            raise ValueError(f"{fid}: domestication changed the selected protein")
        if made["sites_after"] or seqtools.find_sites(made["construct"], assembly.ENZYMES):
            raise ValueError(f"{fid}: domestication left internal Type IIS sites")
        mapping = member_map[fid] if member_map is not None else {}
        cid = mapping.get("candidate_id", f"pending:{fid}")
        key = mapping.get("member_key", f"ort:{fid}")
        if not key.startswith("ort:"):
            raise ValueError(f"{fid}: expected ortholog member key")
        vector, row, clean = assembly.assemble_one(backbone, bb_sites, cid, key, made)
        if not clean:
            raise ValueError(f"{fid}: assembled vector has unexpected Type IIS pattern")
        constructs[fid] = made
        rows.append({**row, "form_id": fid, "sequence_sha256": selected["sequence_sha256"]})
        vectors.append(assembly.fasta_record(cid, key, vector))
    return constructs, rows, vectors


def read_rows(path: Path) -> list[dict]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panel", type=Path, required=True)
    parser.add_argument("--cds", type=Path, required=True)
    parser.add_argument("--backbone", type=Path, required=True, help="pinned p0039 GenBank file")
    parser.add_argument("--backbone-sha256", required=True, help="expected GenBank file checksum")
    parser.add_argument("--member-map", type=Path, help="form_id, member_key, candidate_id TSV")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args(argv)
    if sha256(args.backbone.read_bytes()) != args.backbone_sha256:
        raise ValueError("backbone file checksum differs from the pin")
    panel = read_rows(args.panel)
    mappings = None
    if args.member_map:
        records = read_rows(args.member_map)
        mappings = {r["form_id"]: r for r in records}
        if len(mappings) != len(records):
            raise ValueError("duplicate form_id in member map")
    backbone = str(SeqIO.read(args.backbone, "genbank").seq).upper()
    constructs, rows, vectors = build_forms(
        panel, json.loads(args.cds.read_text()), backbone, mappings
    )
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "compara116_candidate_constructs.json").write_text(
        json.dumps(constructs, indent=2, sort_keys=True) + "\n"
    )
    with (args.out / "compara116_assembly_rows.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
    (args.out / "compara116_pretag_vectors.fasta.gz").write_bytes(
        gzip.compress("".join(vectors).encode(), mtime=0)
    )
    assembly, domestication = construction_tools()
    provenance = {
        "scope": "pre-tag p0039 intermediates; tag cassette deferred; barcode not reserved",
        "construct_count": len(rows),
        "member_ids": "mapped" if mappings is not None else "pending review only",
        "backbone_source": "processed://sprint4/s4e9-orf-design-fg-panel/"
        "p0039_EFS_IRES_puro_backbone.gb",
        "domestication_seed": domestication.SEED,
        "enzymes": domestication.ENZYMES,
        "script_sha256": {
            "builder": sha256(Path(__file__).read_bytes()),
            "domestication": sha256(Path(domestication.__file__).read_bytes()),
            "assembly": sha256(Path(assembly.__file__).read_bytes()),
        },
        "input_sha256": {
            "panel": sha256(args.panel.read_bytes()),
            "cds": sha256(args.cds.read_bytes()),
            "backbone_genbank": args.backbone_sha256,
            **({"member_map": sha256(args.member_map.read_bytes())} if args.member_map else {}),
        },
    }
    (args.out / "compara116_construct_provenance.json").write_text(
        json.dumps(provenance, indent=2) + "\n"
    )
    print(json.dumps(provenance, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
