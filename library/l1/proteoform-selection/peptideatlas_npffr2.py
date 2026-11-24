"""Does observed MS evidence support our NPFFR2 start codon, or the shorter one?

NPFFR2 is the one non-olfactory gene where Tedman and APPRIS both disagree with our
reference and agree with each other: we build a 522-residue protein, they both build 420,
differing by a 102-residue N-terminal extension. Sequence annotation cannot settle which
start codon is used; observed peptides can speak to it.

PeptideAtlas build **Human 2026-01** (atlas_build_id 607, 5,309,343 distinct peptides) is
the orthogonal source: uniformly reprocessed public MS data, not an annotation. Its entry
for Q9Y5X5 is built on the 522-residue UniProt sequence, so its coordinates are ours and
the two are directly comparable -- this is asserted by the test below, not assumed.

The question is not "are peptides seen in residues 1-102" on its own, because absence is
weak evidence for a low-abundance GPCR. It is whether that region is *less* covered than
the rest of the same protein in the same experiments, given how many detectable tryptic
peptides each region should yield.

PeptideAtlas asks automated clients not to query its CGI and to use bulk downloads
instead, so the eight observed peptides below were read once from the protein page rather
than scraped; the page URL and build are recorded so the read is checkable.

    uv run --locked python peptideatlas_npffr2.py
"""

from __future__ import annotations

import csv
import gzip
import json
from pathlib import Path

from scipy.stats import binomtest, fisher_exact

ACCESSION = "Q9Y5X5"
ATLAS_BUILD = "Human 2026-01"
ATLAS_BUILD_ID = 607
ATLAS_URL = (
    "https://db.systemsbiology.net/sbeams/cgi/PeptideAtlas/GetProtein"
    f"?atlas_build_id={ATLAS_BUILD_ID}&protein_name={ACCESSION}&action=QUERY"
)
RETRIEVED = "2026-10-08"

# The disputed region: our canonical's residues 1..DISPUTED_END are absent from the
# 420-residue form that both Tedman and APPRIS build.
DISPUTED_END = 102

# PeptideAtlas reports 8 distinct peptides / 15 total observations for this protein, at
# ProteinProphet 0.999 and identification status "canonical". Read from ATLAS_URL.
OBSERVED = [
    ("PAp01529466", 500, "SAEKPQQELVMEELK", 3),
    ("PAp16492675", 471, "SHVLINTSNQLVQESTFQNPHGETLLYR", 3),
    ("PAp17296865", 469, "AKSHVLINTSNQLVQESTFQNPHGETLLYR", 2),
    ("PAp19127498", 487, "FQNPHGETL", 1),
    ("PAp01665959", 310, "EDWPNQEMR", 1),
    ("PAp12792205", 434, "IYGFFNENF", 1),
    ("PAp16487557", 360, "NQEQWHVVSR", 1),
    ("PAp08112799", 504, "PQQELVMEELK", 3),
]

# Length window within which a tryptic peptide is routinely identified in a shotgun run.
MIN_LEN, MAX_LEN = 7, 30


def tryptic(seq: str) -> list[tuple[int, str]]:
    """Fully-tryptic peptides as (1-based start, sequence). Cleaves after K/R, not before P."""
    cuts = [0]
    cuts += [i + 1 for i, c in enumerate(seq) if c in "KR" and i + 1 < len(seq) and seq[i + 1] != "P"]
    cuts += [len(seq)]
    return [(cuts[i] + 1, seq[cuts[i] : cuts[i + 1]]) for i in range(len(cuts) - 1)]


def main() -> None:
    cache = Path("~/.cache/deliverome/appris").expanduser()
    seq = None
    with gzip.open(cache / "uniprot_human_reviewed_isoforms.tsv.gz", "rt") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["Entry"] == ACCESSION:
                seq = row["Sequence"]
                break
    assert seq, f"{ACCESSION} not found in the cached UniProt stream"
    print(f"{ACCESSION} UniProt canonical: {len(seq)} aa")

    # PeptideAtlas coordinates must be on this sequence, or the comparison is meaningless.
    for pid, start, pep, _ in OBSERVED:
        assert seq[start - 1 : start - 1 + len(pep)] == pep, (
            f"{pid} does not match {ACCESSION} at residue {start}; PeptideAtlas is not "
            "using our numbering and this analysis does not apply"
        )
    print(f"all {len(OBSERVED)} PeptideAtlas peptides sit on our numbering")

    detectable = [(s, p) for s, p in tryptic(seq) if MIN_LEN <= len(p) <= MAX_LEN]
    first_observed = min(s for _, s, _, _ in OBSERVED)

    disputed = [(s, p) for s, p in detectable if s + len(p) - 1 <= DISPUTED_END]
    covered_region = [(s, p) for s, p in detectable if s >= first_observed]
    observed_starts = {s for _, s, _, _ in OBSERVED}

    hit_disputed = sum(1 for s, _ in disputed if s in observed_starts)
    hit_covered = sum(1 for s, _ in covered_region if s in observed_starts)

    table = [
        [hit_disputed, len(disputed) - hit_disputed],
        [hit_covered, len(covered_region) - hit_covered],
    ]
    _, p_value = fisher_exact(table, alternative="less")

    print(f"\ndetectable fully-tryptic peptides ({MIN_LEN}-{MAX_LEN} aa): {len(detectable)}")
    print(
        f"  in the disputed extension (1-{DISPUTED_END}): "
        f"{hit_disputed} of {len(disputed)} observed"
    )
    for s, p in disputed:
        print(f"    {s:4d}-{s + len(p) - 1:4d}  {p}")
    print(
        f"  from residue {first_observed} on: {hit_covered} of {len(covered_region)} observed"
    )
    print(f"\nFisher exact, one-sided (less coverage in the extension): p = {p_value:.4f}")

    # A second test that does not depend on the digest model at all: if observed peptides
    # fell anywhere in the protein with equal chance, the extension would take its share of
    # them in proportion to its length.
    share = DISPUTED_END / len(seq)
    binom = binomtest(0, len(OBSERVED), share, alternative="less")
    print(
        f"binomial, 0 of {len(OBSERVED)} peptides in a region that is "
        f"{share:.1%} of the sequence: p = {binom.pvalue:.4f}"
    )
    verdict = (
        "suggestive" if min(p_value, binom.pvalue) < 0.05 else "not statistically supported"
    )
    print(
        f"\nVerdict: {verdict}. Zero observations in the disputed extension is consistent "
        "with the shorter form, but the coverage is too thin to separate that from the "
        "region simply not having been sampled."
    )

    Path("peptideatlas_npffr2.json").write_text(
        json.dumps(
            {
                "accession": ACCESSION,
                "atlas_build": ATLAS_BUILD,
                "atlas_build_id": ATLAS_BUILD_ID,
                "atlas_url": ATLAS_URL,
                "retrieved": RETRIEVED,
                "canonical_length": len(seq),
                "disputed_extension": [1, DISPUTED_END],
                "distinct_peptides_observed": len(OBSERVED),
                "total_observations": sum(n for *_, n in OBSERVED),
                "detectable_tryptic_in_extension": len(disputed),
                "observed_in_extension": hit_disputed,
                "detectable_tryptic_from_first_observation": len(covered_region),
                "observed_from_first_observation": hit_covered,
                "fisher_one_sided_p": round(p_value, 4),
                "binomial_one_sided_p": round(binom.pvalue, 4),
                "extension_fraction_of_sequence": round(share, 4),
                "verdict": verdict,
            },
            indent=2,
        )
        + "\n"
    )
    print("wrote peptideatlas_npffr2.json")


if __name__ == "__main__":
    main()
