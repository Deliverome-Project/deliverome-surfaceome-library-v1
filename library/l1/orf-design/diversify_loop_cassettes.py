"""Issue one distinct DNA encoding of the loop tag per loop construct.

``diversify_cassettes.py`` does this for the 51-aa terminal cassette and writes
``cassette_encodings.tsv``. The loop arms were never run through it: every one of them
carried the single hard-coded ``LOOP_TAG`` constant, so all four shared the same 57 nt.

A shared constant is the dangerous shape here, and it is worse than a long pairwise
homology even when the pairwise number is bigger. Reverse transcription switches between
the two genomes a virion carries, and a switch anywhere in an identical stretch is
seamless. A 972-nt tract shared by the two SLC2A4 variants is a substrate for that one
pair, and a pair is a small fraction of a pool; a tag every member carries is a substrate
for *every* pair at once. So the long same-gene tracts are accepted (Becca, 2026-10-12 —
two tag positions in one ORF cannot be recoded apart) and the universal one is removed.

The tag peptide stays exactly ``GSPSRLEEELRRRLTEPGS``; only its codons change. Encodings
are written here rather than generated inside ``build_loop_inserts.py`` so the sequence
that goes to synthesis is a reviewed artefact and not the output of an RNG that could
drift between runs, matching how the terminal cassette is handled.

**Encodings already in the file are pinned, not redrawn.** Members are served in
descending swap-risk order and each draw avoids the k-mers already issued, so changing
the construct set reorders the queue and shifts every later draw -- rewriting encodings
for arms that did not change. A variant label carries its own residue and position
(``SLC2A4_P66``), so a key that survives a site change is by construction a construct
whose flanks are untouched, and its issued DNA is passed back through
``seqtools.diversify_library_linkers(existing=...)``. ``--reissue-all`` redraws
everything, which is only correct before anything is built.

    uv run --locked python diversify_loop_cassettes.py
    uv run --locked python diversify_loop_cassettes.py --check   # re-audit, write nothing
"""

from __future__ import annotations

import argparse
import csv
import random
from pathlib import Path

from build_constructs import ENZYMES, read_fasta
from build_gpcr_rescue_arms import arm_flanks
from build_loop_inserts import LOOP_TAG, LOOP_TAG_AA, VARIANTS

from deliverome_analysis.cloning import seqtools

HERE = Path(__file__).parent
OUT = HERE / "loop_cassette_encodings.tsv"
SEED = 2026_10_12  # fixed so the issued encodings are reproducible
FLANK = 150  # context for the junction checks; 150 nt is ample for k=20


def member_flanks(cds: dict[str, str]) -> dict[str, tuple[str, str]]:
    """Every construct carrying the cassette peptide being drawn, keyed by construct name.

    Two kinds share this peptide and therefore have to be drawn against each other, or
    two constructs in one pool could be issued encodings sharing a long tract:

    * loop insertions, where the tag sits mid-ORF with real coding sequence on both sides
      rather than the shared constant tail the terminal cassette gets;
    * the GPCR rescue arms, where it sits after Met1 with the Kozak (and, on the leader
      arm, the Igk leader) upstream.
    """
    out: dict[str, tuple[str, str]] = {}
    for gene, variants in VARIANTS.items():
        body = cds[gene][:-3]  # drop the native stop
        for variant, (after, _) in variants.items():
            cut = after * 3
            out[f"{gene}_{variant}"] = (body[:cut][-FLANK:], body[cut:][:FLANK])
    out.update(arm_flanks(cds, flank=FLANK))
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
    flanks = member_flanks(cds)

    on_file = {}
    if OUT.exists():
        on_file = {
            r["construct"]: r["loop_tag_nt"] for r in csv.DictReader(OUT.open(), delimiter="\t")
        }

    if args.check:
        issued = on_file
    else:
        pinned = {} if args.reissue_all else {k: v for k, v in on_file.items() if k in flanks}
        dropped = sorted(set(on_file) - set(flanks))
        if pinned:
            print(
                f"  pinning {len(pinned)} already-issued encoding(s): {', '.join(sorted(pinned))}"
            )
        if dropped:
            print(f"  dropped (site changed or removed): {', '.join(dropped)}")
        plan = seqtools.diversify_library_linkers(
            flanks,
            LOOP_TAG_AA,
            kmer_policy="minimize",
            enzymes=ENZYMES,
            existing=pinned or None,
            rng=random.Random(SEED),
        )
        issued = {name: v.sequence for name, v in plan.variants.items()}
        for name, value in sorted(plan.audit().items()):
            print(f"  {name:<26} {value}")

    # Same contract as the terminal cassette: synonymous with the reference, distinct
    # from it and from each other, and no Type IIS site created against its own flanks.
    for name, dna in issued.items():
        if seqtools.translate(dna) != LOOP_TAG_AA:
            raise SystemExit(f"{name}: encoding does not translate to the loop tag")
        if dna == LOOP_TAG:
            raise SystemExit(f"{name}: encoding is the reference LOOP_TAG, not a new one")
        f5, f3 = flanks[name]
        sites = seqtools.find_sites(f5 + dna + f3, ENZYMES)
        if sites:
            raise SystemExit(f"{name}: Type IIS site in or across the loop tag: {sites}")
    if len(set(issued.values())) != len(issued):
        raise SystemExit("two loop constructs were issued the same encoding")

    worst = max(
        (seqtools.longest_shared_tract(a, b)[0], x, y)
        for i, (x, a) in enumerate(sorted(issued.items()))
        for y, b in sorted(issued.items())[i + 1 :]
    )
    print(
        f"\n{len(issued)} encodings, all synonymous with the {len(LOOP_TAG)} nt reference"
        f"\n  longest tract shared by any two: {worst[0]} nt ({worst[1]} vs {worst[2]})"
    )

    if not args.check:
        with OUT.open("w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(["construct", "loop_tag_aa", "loop_tag_nt", "gc", "differs_from_ref_nt"])
            for name in sorted(issued):
                dna = issued[name]
                diff = sum(1 for a, b in zip(dna, LOOP_TAG, strict=True) if a != b)
                w.writerow([name, LOOP_TAG_AA, dna, f"{seqtools.gc_content(dna):.3f}", diff])
        print(f"\nwrote {OUT.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
