"""Design the barcode-amplicon reverse primer for the P-0050-backbone set (P-0051..P-0075).

On the P-0039 constructs the barcode amplicon is O-0013 (the 5' handle, forward) with O-0012
(reverse, in IRES2). The P-0050 set has no IRES2 -- AscI is followed directly by the P-0008
SV40 enhancer/promoter -- so O-0012 has no site and the only inventory reverse that binds,
O-0014 in WPRE, gives a 1,145-bp product that spans the whole SV40-PuroR unit.

This picks a reverse primer inside the SV40 enhancer/promoter, which is byte-identical in all
25 members, in p0050 itself and in P-0008, so one pair reads every construct and the amplicon
differs only by barcode. Every candidate window in the first 200 nt after AscI is scored
with ``deliverome_analysis.primers.check_pair`` against O-0013 on the p0053 (AXL) map: it must
trip no threshold on its own, keep the pair's Tm difference inside NEB's 5 degC limit, give one
product on every one of the 25 maps, and the product must be short enough for a 150-cycle
read (150-250 bp preferred). Ties go to the shortest amplicon.

O-0013 is included in every check because it is the fixed partner; its own flags (if any) are
reported but cannot be designed away here.

    uv run --locked python design_p0050_bc_reverse.py            # reads maps from S3
    uv run --locked python design_p0050_bc_reverse.py --maps DIR  # or a local folder of p00NN.gb
"""

from __future__ import annotations

import argparse
import json
import tempfile
from pathlib import Path

from Bio import SeqIO

from deliverome_analysis import constructs, primers
from deliverome_analysis.seqtools import reverse_complement

HERE = Path(__file__).parent
OUT = HERE / "p0050_bc_reverse_primer.json"
FWD = "AGTCGTTAGCTGCCGTTCC"  # O-0013, TFRC_bc_fwd_JT832
ASCI = "GGCGCGCC"
MEMBERS = [f"p{n:04d}" for n in range(51, 76)]
REFERENCE = "p0053"
WINDOW_NT = 200
LENGTHS = range(18, 26)
AMPLICON = (150, 250)


def _load(name: str, maps: Path | None, work: Path) -> str:
    if maps and (maps / f"{name}.gb").exists():
        path = maps / f"{name}.gb"
    else:
        from deliverome_analysis import s3

        path = work / f"{name}.gb"
        s3.processed.download(constructs.map_key(name), path)
    return str(SeqIO.read(path, "genbank").seq).upper()


def candidates(template: str) -> list[dict]:
    """Every reverse-primer window in the SV40 region after AscI, scored as a pair with O-0013."""
    start = template.index(ASCI) + len(ASCI)
    region = template[start : start + WINDOW_NT]
    out = []
    for length in LENGTHS:
        for i in range(0, len(region) - length + 1):
            rev = reverse_complement(region[i : i + length])
            rep = primers.check_pair(FWD, rev, template, fwd_name="O-0013", rev_name="candidate")
            rev_flags = rep.reverse.concerns()
            pair_flags = [c for c in rep.concerns() if not c.startswith(("forward:", "reverse:"))]
            product = rep.products[0][2] if len(rep.products) == 1 else None
            out.append(
                {
                    "sequence": rev,
                    "offset_after_asci": i,
                    "length": length,
                    "tm_c": rep.reverse.tm_c,
                    "tm_difference": rep.tm_difference,
                    "anneal_c": rep.anneal_c,
                    "amplicon_bp": product,
                    "reverse_flags": rev_flags,
                    "pair_flags": pair_flags,
                    "cross_dimer_dg": rep.cross_dimer.dg,
                }
            )
    return out


def _disqualifying(flags: list[str]) -> list[str]:
    """Pair flags that rule a candidate out.

    The conserved-miRNA-seed note is informational: it reports sites in the amplicon span,
    which is construct sequence that exists whether or not this primer is made, so it is
    carried into the report but does not reject a primer.
    """
    return [f for f in flags if "miRNA seed" not in f]


def choose(cands: list[dict]) -> list[dict]:
    """Candidates with a clean reverse primer, a clean pair, one product of usable size.

    Ranking: within NEB's 1 degC of O-0013 (so a single annealing temperature serves both
    backbones), the shortest amplicon wins -- the barcode read is the point, and a shorter
    product sequences and pools more evenly.
    """
    ok = [
        c
        for c in cands
        if not c["reverse_flags"]
        and not _disqualifying(c["pair_flags"])
        and c["amplicon_bp"]
        and AMPLICON[0] <= c["amplicon_bp"] <= AMPLICON[1]
        and c["tm_difference"] <= 1.0
    ]
    return sorted(ok, key=lambda c: (c["amplicon_bp"], c["tm_difference"]))


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument(
        "--maps", type=Path, help="folder of p00NN.gb; default: processed://constructs/designed/"
    )
    args = ap.parse_args(argv)
    work = Path(tempfile.mkdtemp(prefix="p0050rev-"))

    ref = _load(REFERENCE, args.maps, work)
    cands = candidates(ref)
    ranked = choose(cands)
    if not ranked:
        raise SystemExit("no candidate passes every check; loosen AMPLICON or LENGTHS")
    best = ranked[0]

    # The pair must give exactly one product, of the same size, on every member and on p0050.
    sizes: dict[str, int | None] = {}
    for name in [*MEMBERS, "p0050"]:
        prods = primers.amplicons(FWD, best["sequence"], _load(name, args.maps, work))
        sizes[name] = prods[0][2] if len(prods) == 1 else None
    member_sizes = {sizes[m] for m in MEMBERS}
    assert member_sizes == {best["amplicon_bp"]}, sizes
    assert sizes["p0050"] is None, (
        "O-0013 has no site on the empty backbone, so no product expected"
    )

    full = primers.check_pair(
        FWD, best["sequence"], ref, fwd_name="O-0013", rev_name="p0050_bc_rev_SV40"
    )
    print(primers._format_pair(full))
    print(f"\nproduct on all 25 members: {best['amplicon_bp']} bp; on p0050: none (no O-0013 site)")
    print("\nrunners-up:")
    for c in ranked[1:4]:
        print(
            f"  {c['sequence']}  Tm {c['tm_c']:.1f}  dTm {c['tm_difference']:.1f}  {c['amplicon_bp']} bp"
        )
    OUT.write_text(
        json.dumps(
            {
                "forward": {"id": "O-0013", "sequence": FWD},
                "reverse": best,
                "runners_up": ranked[1:6],
                "n_candidates": len(cands),
                "n_passing": len(ranked),
                "products_by_member": sizes,
                "pair_report": full.to_dict(),
            },
            indent=2,
            default=str,
        )
    )
    print(f"\nwrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
