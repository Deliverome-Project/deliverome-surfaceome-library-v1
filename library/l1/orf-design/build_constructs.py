"""Assemble the nine non-Tedman ORFs into the p003 cassette format.

Format taken verbatim from Becca's assembled TFRC construct ``p003``
(``processed://constructs/ordered/p003.gb``). Every constant below was read off that map,
not retyped from the s1e2 README, because two README rows disagree with the map (see
FORMAT.md).

    EFS_v1
      → KOZAK              GACCGGTTCTAGAGCCGCCACC      AgeI + XbaI + gccgccacc
      → ORF                ATG … no stop               domesticated, native codons
      → CASSETTE           GGGGS-SpyTag003-GGGGSGGGGS-P-ALFAv1.5-P-GGGGS   (51 aa)
      → STOP               TAA                         +4 = C from the pad → TAAC (see note)
      → PAD                CC
      → FWD primer site    AGTCGTTAGCTGCCGTTCC         19 nt, 10x-compatible handle
      → BARCODE            20 nt, one per ORF
      → REV primer site    GGATCTGGACCCTTTGGGCAAGCACACGTC   10-nt gap + 20-nt RT site
      → 10x CS1            GCTTTAAGGCCGGTCCTAGCAA           22 nt, capture handle
      → IRES2_v1 → mScarlet3 → P2A → PuroR → WPRE

``--placement nterm`` instead puts the cassette inside the mature N-terminus — after the
signal-peptide cleavage site where there is one, after Met1 where there is not. That is
Tedman's own rule for his signal-peptide-bearing receptors, and it is the topologically
correct option for every gene here except the two with no extracellular terminus at all.

Amplicon handle, revised 2026-10-01 (Becca). The 2-nt pad changed ``GG`` -> ``CC`` and
the 5' handle changed ``TGGACCCGGGTCAGGCTCC`` -> ``AGTCGTTAGCTGCCGTTCC``. Barcode, the
10-nt gap, the RT site and 10x CS1 are untouched, and the new handle carries no XmaI, so
the ``JT_fw_XmaI_mod`` workaround is no longer needed. Both changes are confined to the
constants below; every other script imports them.

One consequence worth stating rather than burying: the stop codon's +4 base moves from G
to C, so the termination context goes ``TAAG`` -> ``TAAC``. s1e2 §4.5 picked ``TAAG``
because in mammalian cells a purine at +4 terminates more efficiently than a pyrimidine.
TAA is still the strongest of the three stop codons and readthrough stays low, but these
are C-terminally tagged constructs, so any readthrough product carries a ~30-residue
extension through the barcode. Making the pad ``GC`` instead of ``CC`` would restore
``TAAG`` at no other cost; it has not been changed without asking.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from deliverome_analysis.cloning import plasmid_qc as qc
from deliverome_analysis.cloning import seqtools

HERE = Path(__file__).parent

# --- p003 constants, canonical home in deliverome_analysis.plasmid_qc --------------- #
KOZAK = qc.KOZAK
# GGGGS | SpyTag003 | GGGGSGGGGS | P | ALFA v1.5 | P | GGGGS
CASSETTE = qc.P003_CASSETTE_NT
CASSETTE_AA = qc.P003_CASSETTE_PEPTIDE
STOP = qc.STOP
SPACER = qc.SPACER  # pad between the stop and the 5' handle
FWD_PRIMER = qc.FWD_PRIMER
REV_PRIMER = qc.REV_PRIMER
# 10x Genomics Capture Sequence 1, added 2026-10-06.
CS1 = qc.CS1
TFRC_BARCODE = "TAGCACTGAGTCAGATCTGA"  # in use — never reissue

# L1 readout region, adopted 2026-11-22 (Becca: "let's just leave the backbone as is but then
# add the t7 and the feldman flank"). The constants above stay as they are: they describe the
# p003 layout that every plasmid ordered before L1 carries, and the map audits check against
# them. L1 inserts are assembled with ``l1_readout`` instead:
#
#   STOP · SPACER · T7_PROMOTER · FWD_PRIMER (JT832) · barcode · FELDMAN_3PRIME · CS1
#
# The definitions, and the reasoning for T7 and Feldman's padlock + T9A RT site, live in
# deliverome_analysis.plasmid_qc; they are re-exported here under the names the L1 scripts use.
T7_PROMOTER = qc.T7_PROMOTER
FELDMAN_3PRIME = qc.FELDMAN_3PRIME
l1_readout = qc.l1_readout
#: The two constant halves of ``l1_readout`` either side of the barcode, for code that reads an
#: L1 insert's 3' end (seam contexts, swap-risk tails, annotation, QC).
L1_BEFORE_BARCODE = qc.L1_BEFORE_BARCODE
L1_AFTER_BARCODE = qc.L1_AFTER_BARCODE

ENZYMES = list(qc.DOMESTICATION_ENZYMES)

# UniProt signal-peptide end (last residue of the SP) and whether a terminus is
# extracellular. Both come from design_orfs.py's UniProt pull, not from memory.
TOPOLOGY = {
    "HEG1": {"sp_end": 29, "n_out": True, "c_out": False},
    "ERBB2": {"sp_end": 22, "n_out": True, "c_out": False},
    "AXL": {"sp_end": 25, "n_out": True, "c_out": False},
    "SLC9A3": {"sp_end": 25, "n_out": True, "c_out": False},
    "VSIG10": {"sp_end": 30, "n_out": True, "c_out": False},
    "NECTIN4": {"sp_end": 31, "n_out": True, "c_out": False},
    # omega is the GPI attachment residue: everything after it is the GPI signal, which
    # the transamidase cleaves off and replaces with the anchor.
    "FOLR1": {"sp_end": 24, "n_out": True, "c_out": False, "gpi_omega": 234},
    "SLC2A4": {"sp_end": 0, "n_out": False, "c_out": False},  # no extracellular terminus
    "MS4A1": {"sp_end": 0, "n_out": False, "c_out": False},  # no extracellular terminus
    "PTPRC": {"sp_end": 25, "n_out": True, "c_out": False},  # CD45: ecto 26-577, TM 578-598
}


def read_fasta(path: Path) -> dict[str, str]:
    return {header.split("|")[0]: seq for header, seq in seqtools.read_fasta(path)}


def load_cassettes() -> dict[str, str]:
    """Each gene's own synonymous encoding of the cassette, from diversify_cassettes.py.

    One shared encoding would leave every C-terminal construct identical from the cassette
    through to the barcode -- 177 nt of seamless template-switch substrate, 8.7% by
    Walton's fit. Per-gene encodings move the last point of difference to the stop codon.
    Missing on purpose rather than falling back to CASSETTE: silently reverting one gene
    to the shared encoding would restore the full tract for every pair it belongs to.
    """
    path = HERE / "cassette_encodings.tsv"
    if not path.exists():
        raise SystemExit(f"{path.name} missing — run diversify_cassettes.py first")
    out = {}
    for r in csv.DictReader(path.open(), delimiter="\t"):
        if seqtools.translate(r["cassette_nt"]) != CASSETTE_AA:
            raise SystemExit(f"{r['gene']}: cassette encoding is not synonymous with p003's")
        out[r["gene"]] = r["cassette_nt"]
    return out


def assemble(gene: str, cds: str, barcode: str, placement: str, cassette: str) -> dict:
    """Kozak + ORF-with-cassette + stop + barcode amplicon, as one synthesis insert."""
    if seqtools.translate(cds)[-1] != "*":
        raise RuntimeError(f"{gene}: CDS does not end in a stop codon")
    body = cds[:-3]  # drop the native stop; the cassette supplies the new one
    topo = TOPOLOGY[gene]

    if placement == "cterm":
        # Mirror of the N-terminal guard below, which this path was missing. A cassette on
        # a C-terminus that is not extracellular cannot be reached by anti-ALFA or
        # SpyCatcher on an intact cell, so it is not a surface-display construct whatever
        # else it is good for.
        if topo.get("gpi_omega"):
            return {
                "gene": gene,
                "placement": placement,
                "skipped": (
                    f"GPI-anchored: the transamidase cuts at the omega residue "
                    f"({topo['gpi_omega']}) and replaces everything downstream with the "
                    f"anchor, so a C-terminal cassette is either removed with the GPI "
                    f"signal or blocks anchoring and the protein is secreted. Use nterm."
                ),
            }
        if not topo["c_out"]:
            return {
                "gene": gene,
                "placement": placement,
                "skipped": "C-terminus is not extracellular — the cassette would face the "
                "cytosol and cannot be stained on an intact cell; use nterm or a loop",
            }
        coding = body + cassette
        note = "cassette at the C-terminus, p003 layout, per-gene codons"
    elif placement == "nterm":
        if not topo["n_out"]:
            return {
                "gene": gene,
                "placement": placement,
                "skipped": "no extracellular terminus — cassette must go in an "
                "extracellular loop; terminal placement is not possible",
            }
        cut = topo["sp_end"] * 3  # nt offset just past the signal peptide (0 if none)
        coding = body[:cut] + cassette + body[cut:]
        note = (
            f"cassette after the signal peptide (residue {topo['sp_end']})"
            if topo["sp_end"]
            else "cassette after Met1"
        )
    else:
        raise ValueError(placement)

    insert = KOZAK + coding + STOP + SPACER + FWD_PRIMER + barcode + REV_PRIMER + CS1
    protein = seqtools.translate(coding + STOP)
    if protein[-1] != "*" or "*" in protein[:-1]:
        raise RuntimeError(f"{gene}: internal stop or missing stop after assembly")
    if CASSETTE_AA not in protein:
        raise RuntimeError(f"{gene}: cassette peptide not intact after assembly")
    sites = seqtools.find_sites(insert, ENZYMES)
    return {
        "gene": gene,
        "placement": placement,
        "note": note,
        "insert_nt": len(insert),
        "cds_nt": len(coding) + 3,
        "protein_aa": len(protein) - 1,
        "barcode": barcode,
        "sites_in_insert": ";".join(f"{e}@{p}{s}" for e, s, p in sites) or "none",
        "sites_n": len(sites),
        "insert": insert,
        "protein": protein.rstrip("*"),
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--placement", choices=["cterm", "nterm"], default="nterm")
    parser.add_argument("--out-dir", type=Path, default=HERE)
    args = parser.parse_args(argv)

    cds = read_fasta(HERE / "orf_domesticated.fasta")
    cassettes = load_cassettes()
    barcodes = {
        r["gene"]: r["barcode"]
        for r in csv.DictReader((HERE / "barcodes.tsv").open(), delimiter="\t")
    }
    for gene, bc in barcodes.items():
        if bc == TFRC_BARCODE:
            raise RuntimeError(f"{gene} was issued TFRC's in-use barcode")

    rows, skipped = [], []
    for gene, seq in cds.items():
        if gene not in barcodes:
            if gene in ("SLC2A4", "MS4A1"):
                why = "no extracellular terminus — built by build_loop_inserts.py"
            else:
                # The Tedman five are domesticated here but not yet assembled: their tag
                # placement is still open (Tedman's 9-aa HA at position 2 versus our
                # 51-aa cassette, and the Igk leader question for CNR1 and CCR5).
                why = "ORF designed; tag placement still open, no barcode issued"
            print(f"{gene:8s} not built — {why}")
            continue
        if gene not in cassettes:
            raise SystemExit(f"{gene}: no cassette encoding issued — rerun diversify_cassettes.py")
        record = assemble(gene, seq, barcodes[gene], args.placement, cassettes[gene])
        if "skipped" in record:
            skipped.append(record)
            print(f"{gene:8s} SKIPPED — {record['skipped']}")
            continue
        rows.append(record)
        print(
            f"{gene:8s} insert {record['insert_nt']:5d} nt  "
            f"protein {record['protein_aa']:5d} aa  "
            f"TypeIIS in insert: {record['sites_in_insert']}  | {record['note']}"
        )

    suffix = args.placement
    table = args.out_dir / f"constructs_{suffix}.tsv"
    cols = [
        "gene",
        "placement",
        "note",
        "insert_nt",
        "cds_nt",
        "protein_aa",
        "barcode",
        "sites_n",
        "sites_in_insert",
    ]
    with table.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=cols, delimiter="\t", extrasaction="ignore", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)

    fasta = args.out_dir / f"constructs_{suffix}.fasta"
    fasta.write_text(
        "".join(
            f">{r['gene']}|insert|{suffix}|bc={r['barcode']}|{r['insert_nt']}nt\n"
            + "\n".join(r["insert"][i : i + 60] for i in range(0, len(r["insert"]), 60))
            + "\n"
            for r in rows
        )
    )

    if skipped:
        (args.out_dir / f"skipped_{suffix}.json").write_text(json.dumps(skipped, indent=2) + "\n")
    print(f"\nWrote {table} and {fasta}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
