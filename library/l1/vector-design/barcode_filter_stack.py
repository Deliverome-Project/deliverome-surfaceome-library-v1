"""s1e2 — the published barcode filter stack, applied to our DETERMINISTIC barcodes.

Each filter below is one published rule, tagged with where it comes from and whether we
verified it against a primary source or are carrying it from a secondhand summary. Applies
only to the designed (deterministic) barcode — a degenerate window cannot be filtered at all,
which is handled separately in barcode_design.py.

SOURCES

  Slutskin  Vainberg Slutskin, Weingarten-Gabbay, Nir, Weinberger & Segal, Nat Commun 2018.
            PMID 29410437. The closest biological precedent: a 12-nt barcode inside a
            mammalian transcript, so they worried the barcode itself might alter expression.
            Reported rules: 12 nt; every barcode >=3 nt from every other; exclude seed regions
            for HIGHLY EXPRESSED miRNAs; exclude AATAAA and its point mutants; exclude cloning
            restriction sites.
            STATUS: VERIFIED verbatim against the open-access copy, PMC5802814. Methods:
            a "12 nt barcode sequence"; barcodes "designed to differ from other barcodes in
            the library by 3 nt or more"; designed to "avoid the introduction of known
            regulatory elements such as seed regions for highly expressed miRNAs, the
            canonical polyadenylation signal (AATAAA), and its point mutants." The seed is
            defined as "the MRE bases complementary to bases 2-7 of the miRNA", i.e. a 6-mer.
            NOTE "3 nt or more" does not say whether that is Hamming or edit distance. Our
            min-Levenshtein-3 is at least as strong either way, so we are not under-shooting.

  FREE      Hawkins, Jones, Finkelstein & Press, PNAS 2018. PMID 29925596.
            STATUS: constraints VERIFIED against the reference implementation source,
            finkelsteinlab/freebarcodes `generate.py`, `possible_barcode_iterator`.

            ⚠️ One correction to the secondhand summary. It reported "no GGC". The code
            forbids a G followed by a C at all -- every GC dinucleotide:

                if prev_seq[-1] == 'G':   # Illumina has higher errors with GGC
                    bad_bases += 'C'

            The comment says GGC; the code implements GC. That is far more restrictive, and
            combined with a 40-60% GC window it is close to self-contradictory, so we treat it
            as optional and measure its cost rather than adopting it silently.

  Shiroguchi Shiroguchi et al., PNAS 2012 -- GC 40-60%, homopolymers <=4. Strictly weaker than
            FREE on homopolymers, so FREE governs.

  Lyons     Lyons et al., Sci Rep 2017 -- min Hamming distance, homopolymer cap, GC range,
            arbitrary blacklisted subsequences. Our min-Levenshtein-3 is stronger than Hamming.

  O'Brien   medRxiv 2026 preprint: FREE barcodes -> sequence-quality filter -> exclude
            mammalian miRNA seeds via miRBase, 35 barcodes per construct.
            STATUS: NOT verified; preprint, cited via ResearchGate in the summary we received.
            Recorded as a recent implementation precedent only, not as a basis for a rule.

SEED DEFINITION. Slutskin defines the seed-matching region as complementarity to miRNA bases
2-7, i.e. a 6-mer. Our existing filter uses the 7mer-m8 target (bases 2-8). The 6-mer is
strictly more aggressive: a 16-mer holds 11 six-mer windows against 10 seven-mer windows, and
the 6-mer space is 4x smaller, so cost rises roughly 4x. Both are measured below.

For a transcribed barcode the sequence to exclude is the REVERSE COMPLEMENT of the miRNA seed,
not the miRNA's own 5' sequence -- that is what an mRNA would present. mirna_seed_sites.tsv
already stores reverse-complemented DNA target sites.

Run:
    uv run python library/l1/vector-design/barcode_filter_stack.py
"""

import csv
import itertools
import random
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
SEED_TSV = HERE / "mirna_seed_sites.tsv"

# Type IIS sites to exclude, both strands. BsaI added per design review: it was never cleared
# from any construct (see README 4.10), so a barcode carrying one is a live hazard if a
# BsaI-based toolkit is ever used, and it costs almost nothing to exclude now.
ENZYME_SITES = {
    "BbsI": ("GAAGAC", "GTCTTC"),
    "Esp3I/BsmBI": ("CGTCTC", "GAGACG"),
    "BsaI": ("GGTCTC", "GAGACC"),
}
ALL_SITES = tuple(s for pair in ENZYME_SITES.values() for s in pair)

GC_MIN, GC_MAX = 0.40, 0.60
TARGET_N = 30_000  # up to 30,000 distinct ORFs
LENGTHS = [12, 13, 14, 15, 16, 17, 18, 20]
SEED = 20260812


def revcomp(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def polyadenylation_hexamers():
    """AATAAA and its single-base point mutants (Slutskin)."""
    out = {"AATAAA"}
    canonical = "AATAAA"
    for i in range(6):
        for b in "ACGT":
            if b != canonical[i]:
                out.add(canonical[:i] + b + canonical[i + 1 :])
    return out


def load_seeds(tiers, width):
    """DNA target sites to exclude. width=7 -> 7mer-m8; width=6 -> Slutskin's bases 2-7.

    The 6-mer target is the 7-mer target minus its first base: target7 = revcomp(seed+m8), and
    dropping seed position 8 removes the 5'-most base of the reverse complement.
    """
    out = set()
    with SEED_TSV.open() as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("target_site"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) >= 3 and p[2] in tiers:
                site = p[0]
                out.add(site if width == 7 else site[1:])
    return out


# ---- individual filters, each returning True when the barcode PASSES ----
def f_gc(bc):
    return GC_MIN <= (bc.count("G") + bc.count("C")) / len(bc) <= GC_MAX


def f_no_triples(bc):
    """FREE: no homopolymer triple, i.e. max run 2. Stricter than our previous run<=3."""
    return all(bc[i] != bc[i + 1] or bc[i + 1] != bc[i + 2] for i in range(len(bc) - 2))


def f_no_gc_dinuc(bc):
    """FREE as implemented: no G followed by C anywhere. Optional -- see docstring."""
    return "GC" not in bc


def f_no_ggc(bc):
    """FREE as described in the secondhand summary. Kept so the two can be compared."""
    return "GGC" not in bc


def f_no_selfcomp3(bc):
    """FREE: no 3-bp self-complementarity between two NON-OVERLAPPING arms.

    The arms must not overlap. A naive "does any 3-mer's reverse complement appear anywhere"
    test rejects `ACGT` outright, because revcomp(ACG) = CGT overlaps it at offset 1 -- that is
    a reading artifact, not a hairpin. FREE's implementation enforces separation between the
    arms (`range(len(prev_seq) - 4)` with the `prev_seq[i+1:i+3]` offset), so this is a
    hairpin/self-annealing rule, not a k-mer blacklist.
    """
    n = len(bc)
    for i in range(n - 2):
        target = revcomp(bc[i : i + 3])
        for j in range(i + 3, n - 2):
            if bc[j : j + 3] == target:
                return False
    return True


def f_no_enzyme(bc):
    return not any(s in bc for s in ALL_SITES)


def f_no_polya(bc, hexes=None):
    hexes = hexes or polyadenylation_hexamers()
    return not any(h in bc for h in hexes)


def f_no_atg(bc):
    return "ATG" not in bc


# 5' splice-donor core. The consensus 5' splice site is MAG|GTRAGT; GTRAG is its
# invariant core (the GT dinucleotide plus the strongly preferred RAG). A barcode sits in
# the transcribed 3'UTR of a lentiviral construct, so a cryptic donor there can splice out
# everything downstream of it -- in the p003 architecture that is the IRES and the
# selection marker. Added 2026-09-26: this criterion was missing from the original stack,
# and pool index 7 (GGTAAGACGCGCCTGAGTAG, opening with the textbook GGTAAG) reached the
# 30k pool because of it and had to be skipped by hand in the P-series registry.
SPLICE_DONOR = re.compile(r"GT[AG]AG")


def f_no_splice_donor(bc):
    """No GT[AG]AG anywhere in the barcode.

    Intrinsic only -- a donor formed across a junction with the flanking primer sites
    depends on flanks the generator cannot know, and has to be screened per construct at
    reservation time.
    """
    return SPLICE_DONOR.search(bc) is None


def make_seed_filter(seeds, width):
    def f(bc):
        return not any(bc[i : i + width] in seeds for i in range(len(bc) - width + 1))

    return f


def _pack(s):
    v = 0
    for ch in s:
        v = (v << 2) | "ACGT".index(ch)
    return v


def build_code(length, min_distance, ok, rng, target, max_tries):
    """Greedy set with guaranteed min Levenshtein distance (deletion-neighbourhood filter)."""
    k = min_distance - 1
    seen, accepted, tries = set(), 0, 0
    while accepted < target and tries < max_tries:
        tries += 1
        bc = "".join(rng.choice("ACGT") for _ in range(length))
        if not ok(bc):
            continue
        nb = (
            [_pack(bc)]
            if k == 0
            else [
                _pack("".join(bc[i] for i in range(length) if i not in d))
                for d in itertools.combinations(range(length), k)
            ]
        )
        if seen.isdisjoint(nb):
            seen.update(nb)
            accepted += 1
    return accepted, tries


def main():
    rng = random.Random(SEED)
    polya = polyadenylation_hexamers()
    seed7 = load_seeds({"broadly_conserved", "conserved"}, 7)
    seed6 = load_seeds({"broadly_conserved", "conserved"}, 6)
    seed7_all = load_seeds({"broadly_conserved", "conserved", "poorly_conserved", "unassigned"}, 7)
    seed6_all = load_seeds({"broadly_conserved", "conserved", "poorly_conserved", "unassigned"}, 6)
    print(f"polyA hexamers (AATAAA + point mutants): {len(polya)}")
    print(f"seed target sites -- conserved+: {len(seed7)} 7-mers, {len(seed6)} 6-mers")
    print(f"                    all human : {len(seed7_all)} 7-mers, {len(seed6_all)} 6-mers")

    # ---- cumulative survival, at the length we care about ----
    stack = [
        ("GC 40-60%", f_gc),
        ("no homopolymer triple (FREE)", f_no_triples),
        ("no ATG anywhere", f_no_atg),
        ("no AATAAA + point mutants (Slutskin)", lambda bc: f_no_polya(bc, polya)),
        ("no BbsI / Esp3I / BsaI, both strands", f_no_enzyme),
        ("no 3-bp self-complementarity (FREE)", f_no_selfcomp3),
        ("no conserved+ 7mer-m8 seed", make_seed_filter(seed7, 7)),
    ]
    optional = [
        ("+ no conserved+ 6mer seed (Slutskin defn)", make_seed_filter(seed6, 6)),
        ("+ no GGC (summary's reading of FREE)", f_no_ggc),
        ("+ no GC dinucleotide (FREE as coded)", f_no_gc_dinuc),
        ("+ no ALL-human 7mer-m8 seed", make_seed_filter(seed7_all, 7)),
    ]

    for L in (16, 20):
        print(f"\n=== cumulative survival, {L}-nt deterministic barcode (MC 120,000) ===")
        n_mc = 120_000
        cands = ["".join(rng.choice("ACGT") for _ in range(L)) for _ in range(n_mc)]
        alive = cands
        print(f"{'filter':44s} {'kept':>9s} {'cum. frac':>10s} {'eff. space':>12s}")
        for name, fn in stack:
            alive = [b for b in alive if fn(b)]
            frac = len(alive) / n_mc
            print(f"{name:44s} {len(alive):>9,} {frac:>10.4f} {frac * 4.0**L:>12.3e}")
        base_alive, base_frac = alive, len(alive) / n_mc
        for name, fn in optional:
            sub = [b for b in base_alive if fn(b)]
            frac = len(sub) / n_mc
            print(
                f"{name:44s} {len(sub):>9,} {frac:>10.4f} {frac * 4.0**L:>12.3e}"
                f"   ({100 * (1 - frac / base_frac):.0f}% of survivors lost)"
            )

    # ---- does the full stack still reach 30,000 at min Levenshtein 3? ----
    def full_ok(bc):
        return all(fn(bc) for _, fn in stack)

    def full_ok_6mer(bc):
        return full_ok(bc) and make_seed_filter(seed6, 6)(bc)

    rows = []
    for label, pred in (("stack (7mer seed)", full_ok), ("stack + 6mer seed", full_ok_6mer)):
        print(f"\n=== min length for {TARGET_N:,} barcodes at min Levenshtein 3 — {label} ===")
        print(f"{'nt':>4} {'built':>9} {'drawn':>10} {'reaches target':>15}")
        for L in LENGTHS:
            got, tries = build_code(L, 3, pred, rng, TARGET_N, 60 * TARGET_N)
            print(f"{L:>4} {got:>9,} {tries:>10,} {'YES' if got >= TARGET_N else 'no':>15}")
            rows.append(
                dict(
                    filter_set=label,
                    barcode_length=L,
                    built=got,
                    target=TARGET_N,
                    reaches_target=int(got >= TARGET_N),
                    candidates_drawn=tries,
                )
            )
            if got >= TARGET_N and L >= 16:
                break

    out = HERE / "barcode_filter_stack.tsv"
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
