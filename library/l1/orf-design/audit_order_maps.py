"""Run the pre-order checks over every GenBank map in the p0022-p0038 batch.

Builds each plasmid's spec from the CURRENT design files -- the disorder-anchored
N-terminal rebuild, the loop constructs, the issued cassette encodings, the barcode
registry and the domesticated ORFs -- then checks the map on S3 against it with
``deliverome_analysis.cloning.plasmid_qc``. A map built from a superseded design passes every
structural check and fails only the comparison against the design, which is the point.

    uv run --locked python audit_order_maps.py
"""

from __future__ import annotations

import csv
from pathlib import Path

from deliverome_analysis.cloning import constructs, seqtools
from deliverome_analysis.cloning import plasmid_qc as qc
from deliverome_analysis.core.s3 import processed

HERE = Path(__file__).parent
BACKBONE = HERE / "p0039_EFS_IRES_puro_backbone.fasta"


def read_fasta(path: Path) -> dict[str, str]:
    out, name = {}, None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split("|")[0]
            out[name] = ""
        elif name:
            out[name] += line.strip()
    return out


def rows(name: str) -> list[dict]:
    with (HERE / name).open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def designed_cassette(insert: str) -> str | None:
    """The tandem's DNA in a designed insert, read in frame from its Kozak ATG."""
    start = insert.find("GCCGCCACCATG") + 9
    orf = insert[start:]
    protein = seqtools.translate(orf[: len(orf) // 3 * 3])
    at = protein.find(qc.P003_CASSETTE_PEPTIDE)
    if at < 0:
        return None
    return orf[3 * at : 3 * (at + len(qc.P003_CASSETTE_PEPTIDE))]


def main() -> int:
    backbone = "".join(BACKBONE.read_text().split("\n")[1:]).upper()
    orfs = read_fasta(HERE / "orf_domesticated.fasta")
    nterm = read_fasta(HERE / "constructs_nterm_disorder_anchored.fasta")
    loop = read_fasta(HERE / "constructs_loop.fasta")
    tandem = {r["gene"]: r["cassette_nt"] for r in rows("cassette_encodings.tsv")}
    # The loop arms carry the full 51-aa tandem (constructs_loop.fasta), like every L1 construct.
    # loop_cassette_encodings.tsv holds the 19-aa ALFA-only loop tags they were first designed
    # with; reading it made every loop plasmid fail against a tag it no longer carries.
    loop_cas = {key: designed_cassette(insert) for key, insert in loop.items()}

    plans: list[tuple[str, qc.Spec]] = []
    for r in rows("barcodes.tsv"):
        g = r["gene"]
        plans.append(
            (
                r["plasmid"],
                qc.Spec(
                    plasmid=f"{r['plasmid']} ({g}, N-terminal)",
                    gene=g,
                    designed_insert=nterm.get(g),
                    backbone=backbone,
                    barcode=r["barcode"],
                    cassette_nt=tandem.get(g),
                    cassette_peptide=qc.P003_CASSETTE_PEPTIDE,
                    proteoforms={f"{g} designed ORF": seqtools.translate(orfs[g]).rstrip("*")},
                ),
            )
        )
    for r in rows("barcodes_loop.tsv"):
        g, key = r["gene"], f"{r['gene']}_{r['variant']}"
        plans.append(
            (
                r["plasmid"],
                qc.Spec(
                    plasmid=f"{r['plasmid']} ({key}, loop)",
                    gene=g,
                    designed_insert=loop.get(key),
                    backbone=backbone,
                    barcode=r["barcode"],
                    cassette_nt=loop_cas.get(key),
                    cassette_peptide=qc.P003_CASSETTE_PEPTIDE,
                    proteoforms={f"{g} designed ORF": seqtools.translate(orfs[g]).rstrip("*")},
                ),
            )
        )

    blocked = 0
    for plasmid, spec in plans:
        key = constructs.map_key(plasmid)
        if not processed.exists(key):
            print(
                f"{spec.plasmid}:\n  [FAIL] map_present               "
                f"no processed://{key}\n  0/1 checks passed — 1 blocking"
            )
            blocked += 1
            continue
        findings = qc.check_genbank(Path(processed.download(key)), spec)
        print(qc.format_report(spec.plasmid, findings))
        blocked += any(not f.ok and f.severity == "error" for f in findings)
    print(f"\n{blocked}/{len(plans)} plasmid(s) have at least one blocking failure")
    return 1 if blocked else 0


if __name__ == "__main__":
    raise SystemExit(main())
