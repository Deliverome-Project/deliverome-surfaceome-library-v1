"""Rebuild every p0022-p0038 GenBank map from the CURRENT design, and publish what passes.

The maps on S3 were assembled from ``constructs_nterm.fasta``, which is superseded: it
places AXL's cassette inside AXL's own signal peptide, and seven of the eight N-terminal
junctions have since moved. ``annotate_insert.SOURCES`` now points at the disorder-anchored
rebuild, so re-running the assembler produces the corrected maps.

Every map is checked with ``deliverome_analysis.plasmid_qc`` against a spec built from the
current design files **before** it is uploaded. A map that fails is written locally and
reported, never published: publishing an unchecked map is how the stale ones got there.

    uv run --locked python rebuild_order_maps.py --dry-run
    uv run --locked python rebuild_order_maps.py
"""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from pathlib import Path

from deliverome_analysis import constructs, seqtools
from deliverome_analysis import plasmid_qc as qc
from deliverome_analysis.s3 import processed

HERE = Path(__file__).parent
#: Outside the experiment folder: upload-experiment-outputs.sh walks every file under
#: an experiment, and these maps already have a home at processed://constructs/designed/.
OUT = Path.home() / ".cache" / "deliverome-analysis" / "rebuilt_maps"
PREFIX = "constructs/designed"


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


def plans() -> list[tuple[str, str, qc.Spec]]:
    """``(plasmid, construct name in the design FASTA, spec)`` for the whole batch."""
    backbone = "".join(
        (HERE / "p0039_EFS_IRES_puro_backbone.fasta").read_text().split("\n")[1:]
    ).upper()
    orfs = read_fasta(HERE / "orf_domesticated.fasta")
    nterm = read_fasta(HERE / "constructs_nterm_disorder_anchored.fasta")
    loop = read_fasta(HERE / "constructs_loop.fasta")
    arms = read_fasta(HERE / "constructs_tag_arms.fasta")
    tandem = {r["gene"]: r["cassette_nt"] for r in rows("cassette_encodings.tsv")}
    # The loop arms draw their 51-aa encoding from cassette_encodings.tsv now, keyed
    # per CONSTRUCT (SLC2A4_E67), not from the retired 19-aa loop_tag_nt column.
    loop_cas = {r["gene"]: r["cassette_nt"] for r in rows("cassette_encodings.tsv")}

    out = []
    for r in rows("barcodes.tsv"):
        g = r["gene"]
        out.append(
            (
                r["plasmid"],
                g,
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
        out.append(
            (
                r["plasmid"],
                key,
                qc.Spec(
                    plasmid=f"{r['plasmid']} ({key}, loop)",
                    gene=g,
                    designed_insert=loop.get(key),
                    backbone=backbone,
                    barcode=r["barcode"],
                    # The loop arms carry the SAME 51-aa cassette as the terminal
                    # constructs as of 2026-10-18, so they declare the p003 peptide.
                    # Declaring LOOP_CASSETTE_PEPTIDE here made the annotation check
                    # compare 153 annotated nt against an expected 57.
                    cassette_nt=loop_cas.get(key),
                    cassette_peptide=qc.P003_CASSETTE_PEPTIDE,
                    proteoforms={f"{g} designed ORF": seqtools.translate(orfs[g]).rstrip("*")},
                ),
            )
        )
    # The recycled GPCR slots. These had no map and had never been through plasmid_qc:
    # plans() read the three original barcode tables only, so the arms that took over
    # p0034/p0037/p0038 were built and barcoded but never checked. An ordered construct
    # that no checker has seen is the gap this whole module exists to close.

    from build_gpcr_rescue_arms import ARMS as _ARMS
    from build_gpcr_rescue_arms import EXTRA_ARMS as _EXTRA
    from build_gpcr_rescue_arms import LEADER_AA as _LEADER_AA
    from build_gpcr_rescue_arms import PLUS1_AA as _PLUS1_AA

    gpcr = read_fasta(HERE / "constructs_gpcr_rescue.fasta")
    for r in rows("barcodes_gpcr.tsv"):
        name = r["construct"]
        if name not in gpcr:
            print(f"  {r['plasmid']}: {name} not in constructs_gpcr_rescue.fasta — skipped")
            continue
        g = r["gene"]
        # The cassette is a property of the ARM, not only of the gene. AVPR2 carries
        # the 51-aa tandem by default, but AVPR2_plain_alfa is the 19-aa variant that
        # exists precisely to differ on cassette size -- reading ARMS[gene] alone
        # declared the wrong peptide, and the proteoform check then could not recover
        # the native chain.
        kind = _EXTRA.get(name, {}).get("cassette") or _ARMS[g]["cassette"]
        out.append(
            (
                r["plasmid"],
                name,
                qc.Spec(
                    plasmid=f"{r['plasmid']} ({name}, N-terminal leader-free)",
                    gene=g,
                    designed_insert=gpcr[name],
                    backbone=backbone,
                    barcode=r["barcode"],
                    # The cassette has to be declared or the proteoform check compares the
                    # TAGGED protein to the untagged reference and fails by exactly the
                    # cassette length (422 vs 371 for AVPR2 = 51 aa; 371 vs 352 for CCR5
                    # = 19 aa). Which peptide depends on the arm: the GPCRs with a free
                    # N-terminal ectodomain carry the 51-aa tandem, the disulfide-tethered
                    # ones the 19-aa ALFA, and each is issued from its own encoding file.
                    # Keyed by CONSTRUCT, not gene: diversify_cassettes.py issues a GPCR
                    # arm its own encoding under `CNR1_plain`, while a plain terminal
                    # construct is keyed by the bare gene. And which peptide an arm
                    # carries is declared in ARMS, so read it there rather than inferring
                    # it from which dict happens to hold the key.
                    cassette_nt=(
                        tandem.get(name) or tandem.get(g)
                        if kind == "tandem"
                        else loop_cas.get(name)
                    ),
                    cassette_peptide=(
                        qc.P003_CASSETTE_PEPTIDE if kind == "tandem" else qc.LOOP_CASSETTE_PEPTIDE
                    ),
                    # An igk arm encodes LEADER + Asp + the receptor, so the untagged
                    # chain is 21 residues longer than the native ORF and the check fails
                    # by exactly that much against the bare proteoform. Declare both: the
                    # native form for the leader-free arms and the leader-bearing form for
                    # the igk ones, since the construct really does encode it before
                    # signal-peptidase cleavage.
                    proteoforms=(
                        {
                            f"{g} L1 designed ORF": seqtools.translate(orfs[g]).rstrip("*"),
                            f"{g} + Igk leader": _LEADER_AA
                            + _PLUS1_AA
                            + seqtools.translate(orfs[g]).rstrip("*"),
                        }
                        if g in orfs
                        else {}
                    ),
                ),
            )
        )

    # The signal-peptide SWAP arms. Unlike an igk GPCR arm, which PREPENDS the leader to a
    # receptor that has none, these REPLACE the member's own signal peptide -- so the
    # untagged chain is neither the native ORF nor the native ORF plus 21 residues, and
    # declaring either would fail the proteoform check by the length of the removed signal
    # peptide. Both forms are declared: the designed ORF for reference, and the chain this
    # construct actually encodes before signal-peptidase cleavage.
    from build_sp_swap_arms import SWAPS as _SWAPS

    swap = read_fasta(HERE / "constructs_sp_swap.fasta")
    for r in rows("barcodes_sp_swap.tsv"):
        name = r["construct"]
        if name not in swap:
            print(f"  {r['plasmid']}: {name} not in constructs_sp_swap.fasta — skipped")
            continue
        g = r["gene"]
        sp_end = _SWAPS[g]["sp_end"]
        native = seqtools.translate(orfs[g]).rstrip("*") if g in orfs else ""
        out.append(
            (
                r["plasmid"],
                name,
                qc.Spec(
                    plasmid=f"{r['plasmid']} ({name}, Igk leader replacing the native SP)",
                    gene=g,
                    designed_insert=swap[name],
                    backbone=backbone,
                    barcode=r["barcode"],
                    cassette_nt=tandem.get(name) or tandem.get(g),
                    cassette_peptide=qc.P003_CASSETTE_PEPTIDE,
                    proteoforms=(
                        {
                            f"{g} L1 designed ORF": native,
                            f"{g} mature + Igk leader (native SP 1-{sp_end} removed)": (
                                _LEADER_AA + _PLUS1_AA + native[sp_end:]
                            ),
                        }
                        if native
                        else {}
                    ),
                ),
            )
        )

    for r in rows("barcodes_tfrc.tsv"):
        key = f"TFRC_{r['arm']}"
        if key not in arms:
            continue
        # Skip arms cut from the batch. Without this the linker 3x/4x arms were rebuilt
        # and published to constructs/designed/pNNNN.gb on every run, which is what kept
        # recreating them under a P-number after they were parked -- and CLAUDE.md is
        # explicit that a parked map is "named for what it is, never with a P-number", so
        # a parked map cannot be mistaken for a registered construct and the number stays
        # free to allocate. Their map of record is constructs/scratch/tfrc-linker-{3x,4x}.gb.
        if r.get("in_order_batch", "1") == "0":
            print(
                f"  skipping {r['plasmid']} (TFRC {r['arm']}): {r.get('batch_note', 'not in the order batch')}"
            )
            continue
        out.append(
            (
                r["plasmid"],
                key,
                qc.Spec(
                    plasmid=f"{r['plasmid']} ({key}, TFRC arm)",
                    gene="TFRC",
                    designed_insert=arms.get(key),
                    backbone=backbone,
                    barcode=r["barcode"],
                    proteoforms={"TFRC designed ORF": seqtools.translate(orfs["TFRC"]).rstrip("*")}
                    if "TFRC" in orfs
                    else {},
                ),
            )
        )
    return out


def _stamp_identity(text: str, plasmid: str, construct: str, sp_by_gene: dict[str, int]) -> str:
    """Rewrite the record's own name to the plasmid it IS, not the backbone it came from.

    ``assemble_vector.py`` builds onto the p0039 backbone and does not know the
    P-number, so it names every record ``p0039_<insert>``. The GenBank LOCUS field is
    capped at 16 characters, so what reached a registry on import was
    ``p0039_ACKR3_plai`` -- the wrong plasmid number AND a truncated construct name, on
    all 23 maps. Benchling lists by LOCUS, so every construct imported under the
    backbone's identity.

    Round-tripping through SeqIO is safe here because the file was written by SeqIO in
    the first place; the caller asserts the sequence and feature count survive.
    """
    from io import StringIO

    from Bio import SeqIO
    from Bio.SeqFeature import FeatureLocation, SeqFeature
    from build_gpcr_rescue_arms import LEADER
    from build_sp_swap_arms import SWAPS as _SWAPS

    rec = SeqIO.read(StringIO(text), "genbank")
    n_features, n_bases = len(rec.features), len(rec.seq)
    rec.id = rec.name = plasmid
    rec.description = f"{plasmid}: {construct} in the p0039 backbone"
    rec.annotations["accession"] = plasmid

    # "ORF, N-terminal" / "ORF, remainder" do not say WHICH protein, which is useless in
    # a registry listing a hundred maps. Name them, and split out what the N-terminal
    # block actually contains: for an igk arm the added Igk leader, for a panel arm the
    # native signal peptide. Both are overlays rather than replacements, so the blocks
    # still tile the insert exactly and annotate_insert's assertion still holds.
    gene = construct.split("_")[0]
    extra: list[SeqFeature] = []
    for feat in rec.features:
        label = feat.qualifiers.get("label", [""])[0]
        if label not in ("ORF, N-terminal", "ORF, remainder"):
            continue
        feat.qualifiers["label"] = [f"{gene} ORF ({label.split(', ')[1]})"]
        feat.qualifiers["gene"] = [gene]
        if label != "ORF, N-terminal":
            continue
        start = int(feat.location.start)
        if construct.endswith("_igk"):
            extra.append(
                SeqFeature(
                    FeatureLocation(start, start + len(LEADER)),
                    type="sig_peptide",
                    qualifiers={
                        # ADDED and REPLACES are different constructs, and a map that
                        # says the wrong one misstates whether the member still has its
                        # own signal peptide. A GPCR arm has none to replace; a swap arm
                        # had one and lost it.
                        "label": (
                            [
                                f"Igk leader (REPLACES native SP 1-{_SWAPS[gene]['sp_end']}, "
                                f"{len(LEADER) // 3} aa)"
                            ]
                            if gene in _SWAPS
                            else [f"Igk leader (ADDED, {len(LEADER) // 3} aa)"]
                        ),
                        "note": [
                            "mouse Ig kappa chain V-III signal peptide, not native to "
                            + gene
                            + (
                                f"; substituted for the native signal peptide 1-"
                                f"{_SWAPS[gene]['sp_end']}, which this cassette abolishes "
                                f"(tagged p(SP) 0.000137)"
                                if gene in _SWAPS
                                else ""
                            )
                        ],
                    },
                )
            )
            extra.append(
                SeqFeature(
                    FeatureLocation(start + len(LEADER), start + len(LEADER) + 3),
                    type="misc_feature",
                    qualifiers={"label": ["Asp spacer"]},
                )
            )
        elif sp_by_gene.get(gene, 0):
            sp_aa = sp_by_gene[gene]
            extra.append(
                SeqFeature(
                    FeatureLocation(start, start + sp_aa * 3),
                    type="sig_peptide",
                    qualifiers={
                        "label": [f"{gene} signal peptide (native, 1-{sp_aa})"],
                        "note": [f"DeepTMHMM v1 on the designed proteoform; cleaved after {sp_aa}"],
                    },
                )
            )
    rec.features.extend(extra)
    rec.features.sort(key=lambda f: (int(f.location.start), -len(f)))
    n_features += len(extra)

    out = StringIO()
    SeqIO.write(rec, out, "genbank")
    result = out.getvalue()
    back = SeqIO.read(StringIO(result), "genbank")
    if len(back.features) != n_features or len(back.seq) != n_bases:
        raise RuntimeError(
            f"{plasmid}: renaming lost content -- {len(back.features)}/{n_features} "
            f"features, {len(back.seq)}/{n_bases} bases"
        )
    return result


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dry-run", action="store_true", help="build and check, upload nothing")
    args = ap.parse_args(argv)

    OUT.mkdir(parents=True, exist_ok=True)
    with (HERE / "library_manifest_l1.tsv").open() as fh:
        sp_by_gene = {
            r["gene_symbol"]: int(r["signal_peptide_aa"])
            for r in csv.DictReader(fh, delimiter="\t")
            if r["candidate_set"] == "library_reference"
        }
    published, held = [], []
    for plasmid, name, spec in plans():
        built = subprocess.run(
            [sys.executable, str(HERE / "assemble_vector.py"), name, "--out-dir", str(OUT)],
            capture_output=True,
            text=True,
            cwd=HERE,
            check=False,
        )
        gb = OUT / f"p0039_{name}_vector.gb"
        if built.returncode != 0 or not gb.exists():
            print(f"{spec.plasmid}: assembly failed\n{built.stderr.strip()[-400:]}")
            held.append(plasmid)
            continue
        findings = qc.check_genbank(gb, spec)
        print(qc.format_report(spec.plasmid, findings))
        if any(not f.ok and f.severity == "error" for f in findings):
            held.append(plasmid)
            continue
        target = OUT / f"{plasmid}.gb"
        target.write_text(_stamp_identity(gb.read_text(), plasmid, name, sp_by_gene))
        if not args.dry_run:
            processed.upload(target, constructs.map_key(plasmid))  # wherever the map is filed now
        published.append(plasmid)

    print(
        f"\n{len(published)} map(s) pass and were "
        f"{'checked only (dry run)' if args.dry_run else "published to each map's current processed://constructs/ folder"}"
        f": {', '.join(published)}"
    )
    if held:
        print(f"{len(held)} held back, not published: {', '.join(held)}")
    return 1 if held else 0


if __name__ == "__main__":
    raise SystemExit(main())
