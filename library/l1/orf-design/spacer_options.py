"""Screen the alternatives to Du's TACT spacer in the miR-L target array.

Each 21-nt repeat is 4 nt that must NOT pair with the miR-L guide, followed by the
17 nt that must. Only position 4 does any pairing work, so it simply must not be A.
T, C and G all satisfy that -- and Du picked T, which is the one that also completes a
conserved miR-374-5p 8mer (TATTATAA) across the junction.

The question this answers is whether C or G is safe to substitute. It is not enough
that they drop miR-374: a different filler can introduce a different seed match, a
restriction or polyadenylation motif, or extra pairing back to the guide, and it
changes the AU context flanking the site, which is one of the known determinants of
miRNA site efficacy.

Every candidate is scored in the real flanked context (AsiSI ... array ... PmeI), not
on the bare 84-mer, because the junction with the cloning sites is itself capable of
creating a site -- that is exactly how the miR-374 one arises.

`--series` instead asks why DIMMER worked in Du's hands despite carrying the site.
Their complementarity series answers it: the miR-374 8mer needs the filler to end in T
immediately before `ATTATA`, which is only reachable at 17 nt. At 18 and 19 nt the site
cannot form for *any* filler, and those designs compensate just as well as 4x17 -- so
compensation tracks site count and complementarity, not miR-374.

    uv run --locked python spacer_options.py
    uv run --locked python spacer_options.py --tiers all_human
    uv run --locked python spacer_options.py --series
"""

from __future__ import annotations

import argparse
from pathlib import Path

from deliverome_analysis.cloning.primers import load_mirna_seed_sites
from deliverome_analysis.cloning.seqtools import gc_content, reverse_complement

REPO = Path(__file__).resolve().parents[3]
SEED_TSV = REPO / "library/l1/vector-design/mirna_seed_sites.tsv"
TIER_SETS = {
    "conserved+": {"broadly_conserved", "conserved"},
    "all_human": {"broadly_conserved", "conserved", "poorly_conserved", "unassigned"},
}

SITE = "ATTATAATGCTTATCTA"  # the 17 nt complementary to the miR-L guide
# Du's fully complementary 1x21nt site, read off the hairpin passenger arm of
# Addgene 218262. Shorter variants pair over fewer nt from the guide's 3' end, which
# eats into this site's 5' end; the eaten bases become non-pairing filler.
FULL21 = "AGGAATTATAATGCTTATCTA"
M374 = "TATTATA"  # revcomp of miR-374-5p positions 2-8
GUIDE = reverse_complement(SITE)
FLANK5, FLANK3 = "GGCGCGCCGCGATCGC", "GTTTAAACGTCC"  # AscI+AsiSI ... PmeI, as built
FORBIDDEN = {
    "CGTCTC": "Esp3I",
    "GAGACG": "Esp3I (rev)",
    "GGTCTC": "BsaI",
    "GAGACC": "BsaI (rev)",
    "GAAGAC": "BbsI",
    "GTCTTC": "BbsI (rev)",
    "AATAAA": "poly(A) signal",
    "ATTAAA": "poly(A) signal (variant)",
    "TTAATTAA": "PacI",
    "ATCGAT": "ClaI",
    "CCCGGG": "XmaI",
    "GGTAAG": "splice donor",
    "TTTTT": "Pol III terminator",
}


def block(spacer: str) -> str:
    return (spacer + SITE) * 4


def context(spacer: str) -> str:
    return FLANK5 + block(spacer) + FLANK3


def seed_hits(seq: str, sites: dict) -> list[tuple[str, str, str, int]]:
    return [
        (seq[i : i + 7], fam, tier, i + 1)
        for i in range(len(seq) - 6)
        for fam, tier in sites.get(seq[i : i + 7], ())
    ]


def extra_guide_pairing(spacer: str) -> int:
    """Longest run of the guide's 3' end that the spacer re-pairs with.

    0 is the design intent: the spacer exists to end complementarity. Anything above 0
    lengthens the duplex beyond the 17 nt Du characterised.
    """
    best = 0
    for n in range(1, len(spacer) + 1):
        tail = spacer[-n:]
        if reverse_complement(tail) == GUIDE[len(GUIDE) - 17 - 0 :][:0] or False:
            pass
        # the base 5' of the site pairs with the guide base just 3' of the 17-mer
        if reverse_complement(spacer[-n:]) == GUIDE[:0]:
            pass
    # direct test: does the extended site extend the perfect duplex?
    ext = spacer + SITE
    for n in range(1, len(spacer) + 1):
        sub = ext[len(spacer) - n :]
        if reverse_complement(sub) in (GUIDE + "N" * 8):
            best = n
    return best


def au_context(spacer: str) -> float:
    """AU fraction of the 8 nt flanking each site on the 5' side.

    Local AU richness is one of the determinants of miRNA site efficacy, so a filler
    change is not purely neutral for miR-L even though it does no base pairing.
    """
    unit = spacer + SITE
    win = (unit * 2)[len(unit) - 8 + len(spacer) : len(unit) + len(spacer)]
    return sum(c in "AT" for c in win) / len(win)


def series_report() -> None:
    """Which of Du's complementarity variants can carry the miR-374 site at all."""
    from itertools import product

    note = {17: "yes (Fig 1H)", 18: "yes (Fig 3B)", 19: "yes (Fig 3B)", 21: "no (Fig 1F)"}
    print("complementarity  filler  fillers that create a miR-374 8mer   compensates")
    for n in (21, 19, 18, 17):
        comp, flen = FULL21[21 - n :], 21 - n
        bad = ["".join(f) for f in product("ACGT", repeat=flen) if M374 in ("".join(f) + comp) * 4]
        frac = f"{len(bad)} of {4**flen}"
        print(f"  {n} nt{'':<11}{flen} nt   {frac:<36} {note[n]}")
    print(
        "\nThe site is only reachable at 17 nt, and only for a filler ending in T.\n"
        "Du's TACT is one of those; 4x18 and 4x19 cannot carry it under any filler,\n"
        "and they compensate as well as 4x17. So the site is not what makes or breaks\n"
        "dosage compensation -- which is why it is worth testing rather than assuming."
    )


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tiers", choices=list(TIER_SETS), default="conserved+")
    ap.add_argument("--series", action="store_true", help=series_report.__doc__)
    args = ap.parse_args(argv)
    if args.series:
        series_report()
        return 0
    sites = load_mirna_seed_sites(frozenset(TIER_SETS[args.tiers]), path=SEED_TSV)
    print(f"seed table: {len(sites):,} target sites at tier {args.tiers}\n")

    base = context("TACT")
    for spacer in ("TACT", "TACC", "TACG", "TACA"):
        ctx = context(spacer)
        hits = seed_hits(ctx, sites)
        motifs = [f"{m} ({n})" for m, n in FORBIDDEN.items() if m in ctx]
        pairs = spacer[-1] == "A"
        tag = "  <- Du's choice, as built" if spacer == "TACT" else ""
        print(f"{spacer}{tag}")
        print(f"   re-pairs with the guide: {'YES - unusable' if pairs else 'no'}")
        print(
            f"   GC of the 84-mer: {gc_content(block(spacer)) * 100:.1f}%"
            f"   5' AU context of each site: {au_context(spacer) * 100:.0f}%"
        )
        print(f"   forbidden motifs: {motifs or 'none'}")
        print(f"   {args.tiers} seed sites in the flanked block: {len(hits)}")
        for site, fam, tier, pos in hits:
            print(f"      {site}  seed {fam}  ({tier})  at {pos}")
        if spacer != "TACT":
            changed = sum(1 for a, b in zip(base, ctx, strict=True) if a != b)
            print(f"   differs from the built array at {changed} of {len(base)} positions")
        print()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
