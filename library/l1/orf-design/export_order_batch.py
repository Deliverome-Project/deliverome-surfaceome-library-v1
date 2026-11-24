"""Export the order batch as named GenBank files, for cross-checking in Benchling.

The plasmid number alone does not say what a construct IS, so a folder of ``p0046.gb``
is unreadable next to a registry. Each file is named for its architecture instead:

    p0022_HEG1_nativeSP_Nterm.gb          tag after the native signal peptide
    p0046_ACKR3_noSP_Nterm_leaderfree.gb  no signal peptide, nothing added
    p0043_CCR5_noSP_Nterm_igkLeader.gb    no signal peptide, Igk leader added
    p0048_VSIG10_swapSP30_Nterm_igkLeader.gb  native signal peptide REPLACED by the leader
    p0032_MS4A1_noSP_loop146.gb           tag inside an extracellular loop
    p0037_TFRC_noSP_Cterm.gb              tag at the extracellular C-terminus

The signal-peptide term names the boundary the CONSTRUCT WAS BUILT AGAINST -- the
``sp_end`` rule in ``rebuild_disorder_anchored.JUNCTIONS``, which is SignalP 6 -- and not
the L1 manifest's ``signal_peptide_aa``, which is DeepTMHMM v1. The two predictors
disagree for two genes in this batch, and reading the manifest put the wrong number in
both filenames: ERBB2 exported as ``nativeSP21`` against a build boundary of 22, and
VSIG10 as ``nativeSP27`` against 30. VSIG10 is the case that matters, because the 3-residue
gap is on the one construct whose signal peptide actually fails -- SignalP gives its
cleavage site only 0.348 confidence -- so a filename asserting a crisp 1-27 boundary
claimed precision the data does not support, on the construct least entitled to it.

For a gene with no entry in ``JUNCTIONS`` the manifest value still decides, which is how
the GPCR and loop arms get "noSP". All five panel GPCRs carry none, and UniProt annotates
no SIGNAL feature for any of them either, so "noSP" on a GPCR is the native state rather
than something the design removed.

The GenBank DEFINITION line is rewritten to the same description, because that is what
Benchling shows in a list.

    uv run --locked python export_order_batch.py --out-dir ~/Downloads/s4e9-order
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import shutil
from pathlib import Path

from deliverome_analysis import constructs

HERE = Path(__file__).parent
#: Where rebuild_order_maps.py leaves the assembled maps.
CACHE = Path.home() / ".cache" / "deliverome-analysis" / "rebuilt_maps"
TABLES = (
    ("barcodes.tsv", None),
    ("barcodes_loop.tsv", None),
    ("barcodes_tfrc.tsv", "in_order_batch"),
    ("barcodes_gpcr.tsv", None),
    ("barcodes_sp_swap.tsv", None),
)


def _rows() -> list[dict]:
    out: list[dict] = []
    for table, gate in TABLES:
        path = HERE / table
        if not path.exists():
            continue
        with path.open() as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if gate and r.get(gate) != "1":
                    continue
                out.append(r)
    return out


def _locus(plasmid: str, gene: str, suffix: str) -> str:
    """A LOCUS name carrying the architecture, inside GenBank's 16-character cap.

    Benchling lists by LOCUS, so this is the string that becomes the construct name on
    import. 16 characters is the hard limit -- exceeding it is what produced
    ``p0039_ACKR3_plai`` -- so the separator is dropped, then the gene truncated, rather
    than letting the suffix fall off the end where the architecture lives.
    """
    for candidate in (f"{plasmid}_{gene}_{suffix}", f"{plasmid}_{gene}{suffix}"):
        if len(candidate) <= 16:
            return candidate
    room = 16 - len(plasmid) - 1 - len(suffix)
    return f"{plasmid}_{gene[:room]}{suffix}"


def _describe(row: dict, sp_aa: int, site: str | None) -> tuple[str, str, str]:
    """(filename stem, LOCUS name, human description) for one construct."""
    gene = row["gene"]
    construct = row.get("construct") or ""
    sp = f"nativeSP{sp_aa}" if sp_aa else "noSP"
    sp_text = f"native signal peptide 1-{sp_aa}" if sp_aa else "no signal peptide"
    if construct.endswith("_igk") and row.get("arm") == "igk_leader_swap":
        # A swap arm REMOVED its signal peptide. Calling it "nativeSP30" would say the
        # construct keeps a signal peptide it does not have, and "noSP" would say the
        # protein never had one. Neither is true, so this arm gets its own term.
        return (
            f"{row['plasmid']}_{gene}_swapSP{sp_aa}_Nterm_igkLeader",
            _locus(row["plasmid"], gene, "igk"),
            f"native signal peptide 1-{sp_aa} REPLACED by the Igk leader + Asp; "
            f"51-aa cassette after the first mature residue",
        )
    if construct.endswith("_igk"):
        return (
            f"{row['plasmid']}_{gene}_{sp}_Nterm_igkLeader",
            _locus(row["plasmid"], gene, "igk"),
            f"{sp_text}; Igk leader + Asp ADDED upstream; 51-aa cassette after the leader",
        )
    if construct.endswith("_plain"):
        return (
            f"{row['plasmid']}_{gene}_{sp}_Nterm_leaderfree",
            _locus(row["plasmid"], gene, "lf"),
            f"{sp_text}; nothing added upstream; 51-aa cassette after Met1",
        )
    if site:
        return (
            f"{row['plasmid']}_{gene}_{sp}_loop{site}",
            _locus(row["plasmid"], gene, f"L{site}"),
            f"{sp_text}; 51-aa cassette inside an extracellular loop, after residue {site}",
        )
    if gene == "TFRC":
        return (
            f"{row['plasmid']}_{gene}_{sp}_Cterm",
            _locus(row["plasmid"], gene, "C"),
            f"{sp_text}; type II, 51-aa cassette at the extracellular C-terminus",
        )
    return (
        f"{row['plasmid']}_{gene}_{sp}_Nterm",
        _locus(row["plasmid"], gene, "sp"),
        f"{sp_text}; 51-aa cassette immediately after it, at residue {sp_aa + 1}",
    )


def _retitle(text: str, locus: str, description: str) -> str:
    """Set LOCUS, ACCESSION and DEFINITION -- the fields a registry names the file by.

    Asserts the sequence and feature count survive, because a rename that silently
    dropped annotation would be worse than the unhelpful name it replaced.
    """
    from io import StringIO

    from Bio import SeqIO

    rec = SeqIO.read(StringIO(text), "genbank")
    n_features, n_bases = len(rec.features), len(rec.seq)
    rec.id = rec.name = locus
    rec.description = description
    rec.annotations["accession"] = locus
    out = StringIO()
    SeqIO.write(rec, out, "genbank")
    result = out.getvalue()
    back = SeqIO.read(StringIO(result), "genbank")
    if len(back.features) != n_features or len(back.seq) != n_bases:
        raise RuntimeError(f"{locus}: renaming lost content")
    return result


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args(argv)
    dest = args.out_dir.expanduser()
    (dest / "genbank").mkdir(parents=True, exist_ok=True)

    with (HERE / "library_manifest_l1.tsv").open() as fh:
        sp = {
            r["gene_symbol"]: int(r["signal_peptide_aa"])
            for r in csv.DictReader(fh, delimiter="\t")
            if r["candidate_set"] == "library_reference"
        }
    # The build rule wins where it exists: a filename has to describe the construct, and
    # the construct was cut at sp_end, not at the manifest's DeepTMHMM boundary.
    from rebuild_disorder_anchored import JUNCTIONS as _JUNCTIONS

    sp.update({g: spec["sp_end"] for g, spec in _JUNCTIONS.items()})
    with (HERE / "constructs_loop.tsv").open() as fh:
        sites = {
            (r["gene"], r["variant"]): r["insert_after_residue"]
            for r in csv.DictReader(fh, delimiter="\t")
        }

    manifest = []
    for row in sorted(_rows(), key=lambda r: r["plasmid"]):
        src = CACHE / f"{row['plasmid']}.gb"
        if not src.exists():
            raise SystemExit(f"no rebuilt map for {row['plasmid']}; run rebuild_order_maps.py")
        site = sites.get((row["gene"], row.get("variant") or ""))
        stem, locus, description = _describe(row, sp[row["gene"]], site)
        if len(locus) > 16:
            raise SystemExit(f"{locus!r} exceeds the GenBank 16-char LOCUS cap")
        out = dest / "genbank" / f"{stem}.gb"
        out.write_text(_retitle(src.read_text(), locus, description))
        manifest.append(
            {
                "plasmid": row["plasmid"],
                "gene": row["gene"],
                "file": out.name,
                "locus": locus,
                "architecture": description,
                "barcode": row.get("barcode", ""),
                "sha256_gb": hashlib.sha256(out.read_bytes()).hexdigest()[:16],
                "s3_key": "processed://" + constructs.map_key(row["plasmid"]),
            }
        )

    with (dest / "MANIFEST.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(manifest[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(manifest)

    fig = HERE / "figures" / "s4e9_order_batch_map.png"
    if fig.exists():
        shutil.copy2(fig, dest / fig.name)
    print(f"{len(manifest)} GenBank files -> {dest / 'genbank'}")
    for m in manifest:
        print(f"   {m['file']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
