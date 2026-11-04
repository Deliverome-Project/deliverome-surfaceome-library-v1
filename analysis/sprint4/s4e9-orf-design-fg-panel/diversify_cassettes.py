"""Issue one distinct DNA encoding of the tag cassette per gene.

Every C-terminal construct ends ORF | CASSETTE | STOP | pad | handle | barcode. The
cassette and everything after it are the same constant in all of them, so walking back
from the barcode the constructs stay identical for 177 nt. Reverse transcription switches
between the two genomes a virion carries, and a switch anywhere in an identical stretch is
seamless, so that 177 nt is a barcode/ORF swap substrate for every pair at once -- 8.7% by
Walton's fit, against the <=100 bp (5.0%) bar check_swap_risk.py enforces.

The cassette peptide has to stay exactly what p003 encodes; its codons do not. Giving each
gene a different synonymous encoding moves the last point of difference from the start of
the cassette to the stop codon, leaving only STOP + pad + handle shared -- 24 nt, 1.3%.

This writes the encodings to ``cassette_encodings.tsv`` rather than generating them inside
build_constructs.py, so the sequence that goes to synthesis is a reviewed artefact and not
the output of an RNG that could drift between runs. build_constructs.py reads that file.

The p003 ``CASSETTE`` constant stays as the reference encoding and is not reissued: p003
is a built plasmid, and the TFRC linker-arm constructs are its arms.

**Encodings already in the file are pinned, not redrawn.** Members are served in descending
swap-risk order and each draw avoids the k-mers already issued, so adding a member reorders
the queue and shifts every later draw. Issuing for ten keys instead of eight moved seven of
the eight built as p0022-p0029 when this was tried without pinning; only AXL, which sorts
first, held. Those eight are passed back through
``seqtools.diversify_library_linkers(existing=...)``. ``--reissue-all`` redraws everything,
which is only correct before anything is built.

**Only the GGGGS linkers are recoded; SpyTag003 and P-ALFA-P keep p003's codons**
(Becca, 2026-10-29). The 2026-10-07 scheme recoded all 51 codons, so the 25-construct
order carries 22-23 distinct encodings of each tag. The tags are what the readout binds
(SpyCatcher, anti-ALFA nanobody), and p003's encoding of them is the one validated in
cells; the swap tract is set by the 3'-most GGGGS, which still varies. The 23 encodings
already built are pinned as issued and are labelled ``whole_cassette`` in the ``scheme``
column; every new draw is ``linkers_only``.

The end-uniqueness rule widens with the scheme. Every C-terminal member shares the 3'
flank, so each needs its own 3' end. The 2026-10-07 rule asked for distinct 9-nt ends,
and only 60 of the 80 GGS encodings are free of a 5-base run -- which is why the encoder
ran out after ~60 members when the library's 216 C-terminal members were drawn at once.
12-nt ends give 180, still short; the whole 15-nt terminal GGGGS gives 540.
``LINKER_ONLY_JUNCTION_IDENTITY`` = 14 therefore makes that block unique per member, which
caps the tract any two C-terminal members share before the barcode at 24 + 14 = 38 nt.
See ``l1_swap_risk.py``.

Rows are keyed by gene for the one-construct-per-gene panel members and by construct name
for the GPCR rescue arms, which need two encodings for one gene.

    uv run --locked python diversify_cassettes.py
    uv run --locked python diversify_cassettes.py --check   # re-audit, write nothing
"""

from __future__ import annotations

import argparse
import csv
import random
from pathlib import Path

from build_constructs import (
    CASSETTE,
    CASSETTE_AA,
    ENZYMES,
    FWD_PRIMER,
    SPACER,
    STOP,
    TOPOLOGY,
    read_fasta,
)
from build_gpcr_rescue_arms import arm_flanks
from build_sp_swap_arms import sp_swap_flanks

from deliverome_analysis import seqtools

HERE = Path(__file__).parent
OUT = HERE / "cassette_encodings.tsv"
SEED = 2026_10_07  # fixed so the issued encodings are reproducible
FLANK = 150  # how much context the junction checks get; 150 nt is ample for k=20

#: Held on p003's codons in every new encoding: SpyTag003, and ALFA v1.5 with its flanking
#: prolines. Everything else in CASSETTE_AA is a GGGGS block and is recoded.
TAG_SEGMENTS = ("RGVPHIVMVDAYKRYK", "PSRLEEELRRRLTEP")
#: How far a cassette may stay identical to a same-flank partner at each junction.
#: Every C-terminal member shares the 3' flank (stop, pad, handle), so this bounds the
#: shared tract before the barcode at 24 + this many nt.
WHOLE_CASSETTE_JUNCTION_IDENTITY = seqtools.DEFAULT_MAX_JUNCTION_IDENTITY  # 2026-10-07, 9-nt ends
LINKER_ONLY_JUNCTION_IDENTITY = 14  # 2026-10-29, the whole 15-nt terminal GGGGS
SCHEME_WHOLE, SCHEME_LINKERS = "whole_cassette", "linkers_only"


def tag_codons() -> dict[int, str]:
    """The fixed tag codons, read off p003's validated cassette."""
    return seqtools.pin_codons(CASSETTE, CASSETTE_AA, TAG_SEGMENTS)


def scheme_of(dna: str) -> str:
    """Which scheme an encoding follows: tags on p003's codons, or recoded."""
    fixed = tag_codons()
    return (
        SCHEME_LINKERS
        if all(dna[3 * i : 3 * i + 3] == c for i, c in fixed.items())
        else SCHEME_WHOLE
    )


def loop_flanks(cds: dict[str, str], flank: int = 150) -> dict[str, tuple[str, str]]:
    """Each loop construct's sequence either side of its insertion point."""
    from build_loop_inserts import VARIANTS

    out: dict[str, tuple[str, str]] = {}
    for gene, variants in VARIANTS.items():
        if gene not in cds:
            continue
        body = cds[gene][:-3]
        for variant, (after, _why) in variants.items():
            cut = after * 3
            out[f"{gene}_{variant}"] = (body[:cut][-flank:], body[cut:][:flank])
    return out


def member_flanks(cds: dict[str, str], genes: list[str]) -> dict[str, tuple[str, str]]:
    """Each gene's sequence either side of the cassette in its C-terminal construct.

    The C-terminal placement is the one that sets the risk -- there the cassette sits
    between the ORF and the barcode with nothing gene-specific after it -- so the
    encodings are chosen against those flanks and reused for the N-terminal build.
    """
    out = {}
    for gene in genes:
        body = cds[gene][:-3]  # native stop dropped; the cassette supplies the new one
        out[gene] = (body[-FLANK:], (STOP + SPACER + FWD_PRIMER)[:FLANK])
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--check", action="store_true", help="audit the existing file only")
    ap.add_argument(
        "--reissue-all",
        action="store_true",
        help="redraw every encoding instead of pinning the ones already issued",
    )
    args = ap.parse_args(argv)

    cds = read_fasta(HERE / "orf_domesticated.fasta")
    barcodes = {r["gene"] for r in csv.DictReader((HERE / "barcodes.tsv").open(), delimiter="\t")}
    genes = sorted(g for g in cds if g in barcodes and g in TOPOLOGY)
    flanks = member_flanks(cds, genes)
    flanks.update(arm_flanks(cds, flank=FLANK, cassette="tandem"))
    # The loop arms carry the 51-aa cassette as of 2026-10-18, so they are drawn here too.
    # Per CONSTRUCT, not per gene: SLC2A4 has two loop sites and MS4A1 two, and issuing
    # one encoding per gene would put identical codons in both arms of a pair -- which is
    # exactly the shared tract the diversifier exists to shorten, and on the pair whose
    # difference is the experiment.
    flanks.update(loop_flanks(cds, flank=FLANK))
    # The signal-peptide swap arms. Their 5' flank carries the Igk leader in place of
    # the native signal peptide, so an encoding drawn against the native arm's flank
    # would be screened against sequence this construct does not contain.
    flanks.update(sp_swap_flanks(cds, flank=FLANK))

    on_file = {}
    if OUT.exists():
        on_file = {r["gene"]: r["cassette_nt"] for r in csv.DictReader(OUT.open(), delimiter="\t")}
    fixed = tag_codons()

    if args.check:
        issued = on_file
    else:
        pinned = {} if args.reissue_all else {k: v for k, v in on_file.items() if k in flanks}
        if pinned:
            print(
                f"  pinning {len(pinned)} already-issued encoding(s): {', '.join(sorted(pinned))}"
            )
        dropped = sorted(set(on_file) - set(flanks))
        if dropped:
            print(f"  dropped (no longer a member): {', '.join(dropped)}")
        plan = seqtools.diversify_library_linkers(
            flanks,
            CASSETTE_AA,
            kmer_policy="minimize",
            enzymes=ENZYMES,
            existing=pinned or None,
            fixed_codons=fixed,
            max_junction_identity=LINKER_ONLY_JUNCTION_IDENTITY,
            rng=random.Random(SEED),
        )
        issued = {g: v.sequence for g, v in plan.variants.items()}
        for name, value in sorted(plan.audit().items()):
            print(f"  {name:<26} {value}")

    # Every encoding must be synonymous with p003's, distinct from it and from each other,
    # and must not create a Type IIS site against its own flanks.
    for gene, dna in issued.items():
        if seqtools.translate(dna) != CASSETTE_AA:
            raise SystemExit(f"{gene}: encoding does not translate to the p003 cassette")
        if dna == CASSETTE:
            raise SystemExit(f"{gene}: encoding is the p003 reference, not a new one")
        if gene not in on_file and scheme_of(dna) != SCHEME_LINKERS:
            raise SystemExit(f"{gene}: a new encoding recodes a tag; only linkers may vary")
        f5, f3 = flanks[gene]
        sites = seqtools.find_sites(f5 + dna + f3, ENZYMES)
        if sites:
            raise SystemExit(f"{gene}: Type IIS site in or across the cassette: {sites}")
    if len(set(issued.values())) != len(issued):
        raise SystemExit("two genes were issued the same encoding")

    worst = max(
        (seqtools.longest_shared_tract(a, b)[0], x, y)
        for i, (x, a) in enumerate(sorted(issued.items()))
        for y, b in sorted(issued.items())[i + 1 :]
    )
    tail = len(STOP + SPACER + FWD_PRIMER)
    print(
        f"\n{len(issued)} encodings, all synonymous with p003's {len(CASSETTE)} nt cassette"
        f"\n  longest tract shared by any two: {worst[0]} nt ({worst[1]} vs {worst[2]})"
        f"\n  constant tail left before the barcode: {tail} nt = "
        f"{seqtools.barcode_swap_rate(tail):.1%} predicted swapping"
    )

    if not args.check:
        with OUT.open("w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(
                ["gene", "cassette_aa", "cassette_nt", "gc", "differs_from_p003_nt", "scheme"]
            )
            for gene in sorted(issued):
                dna = issued[gene]
                diff = sum(1 for a, b in zip(dna, CASSETTE, strict=True) if a != b)
                gc = f"{seqtools.gc_content(dna):.3f}"
                w.writerow([gene, CASSETTE_AA, dna, gc, diff, scheme_of(dna)])
        print(f"\nwrote {OUT.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
