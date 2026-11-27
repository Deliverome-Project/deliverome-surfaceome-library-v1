"""Esp3I/BsmBI shuttle sites flanking the insert, for one-pot transfer to another backbone.

s1e2 §4.11 holds the analysis; this is the build. Two sites, pointing INWARD so both
recognition sequences stay behind in the donor and the released fragment carries no scar:

    …EFS promoter… CGTCTC N [OH5] Kozak-ORF-cassette-stop-barcode [OH3] N GAGACG …IRES2…
                   <--- cuts right --->                          <--- cuts left ---

Esp3I is CGTCTC(1/5): one spacer base, then a 4-nt 5' overhang. `GAGACG` is the reverse
complement, so the downstream site reads CGTCTC on the bottom strand and cuts leftward.

**Overhangs are chosen on NEB ligation-fidelity data**, dataset "Esp3I, 1x T4 DNA Ligase
buffer, 37-16 cycling" (NEBridge Ligase Fidelity Viewer, read 2026-09-22):

    pair            fidelity   correct-pair counts   mismatch counts
    GACC + CGTC       95%        217 / 232            8 (GACC:CGTC), 4 (GGTC:GACG)
    GACC + TCAG      100%        217 / 355            0          <-- chosen
    GACC + TCCA      100%        217 / 317            0
    GACC + TCAC      100%        217 / 239            0
    GACC + CACG      100%        217 / 177            0

The obvious choice, OH3 = CGTC (the natural last four bases of the reverse primer site),
loses 5% because GACC differs from GACG -- CGTC's ligation partner -- at a single
position, so the two ends cross-ligate. TCAG is clean and is also the strongest ligator
in the set. It costs a 2-nt pad (AG) after the reverse primer site, inside the
transferred fragment but in the non-coding gap before the IRES, so nothing coding,
no primer site and no barcode is touched.

OH5 = GACC is the first four bases of the p003 Kozak block -- existing sequence, no scar
-- and is retained: it pairs cleanly with TCAG and moving it would mean either splitting
the Kozak across the junction or padding the 5'UTR.

Why Esp3I and not BbsI or BsaI (s1e2 §4.10/§4.11): BbsI is the assembly enzyme and
cannot double as the shuttle; BsaI has 2-4 sites in every backbone scanned and would cut
it. Esp3I is at 0 sites in p003, p004 and pRDA_722, and is what Tedman used for the same
purpose. It is nominally reserved for guide cloning (s1e2 §1) — that reservation is the
only objection and is contested.

The ORFs here are domesticated against BbsI + Esp3I + BsaI, so these two are the only
Esp3I sites in the insert. This module proves that, and simulates the digest.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from build_constructs import (
    CS1,
    ENZYMES,
    KOZAK,
    read_fasta,
)

from deliverome_analysis.cloning import seqtools

HERE = Path(__file__).parent

ESP3I = "CGTCTC"
ESP3I_RC = "GAGACG"
SPACER_BASE = "A"  # the single N between recognition site and overhang
OH5 = KOZAK[:4]  # GACC — first four bases of the Kozak block
PAD = "AG"  # 2-nt pad after the insert's 3' end, to reach the OH3 overhang
OH3 = (CS1 + PAD)[-4:]  # AAAG — the insert now ends at 10x CS1, not at the RT site

UPSTREAM = ESP3I + SPACER_BASE  # sits between the promoter and OH5
DOWNSTREAM = PAD + SPACER_BASE + ESP3I_RC  # pad, spacer, then the reverse-orientation site


def check_overhangs() -> list[str]:
    """Structural requirements for a 2-part Golden Gate: distinct, non-palindromic,
    and not each other's complement (which would let the fragment circularise or
    insert backwards)."""
    problems = []
    if OH5 == OH3:
        problems.append(f"overhangs identical ({OH5})")
    for name, oh in (("OH5", OH5), ("OH3", OH3)):
        if oh == seqtools.reverse_complement(oh):
            problems.append(f"{name} {oh} is palindromic — self-ligates")
    if OH5 == seqtools.reverse_complement(OH3):
        problems.append(f"OH5 {OH5} is the complement of OH3 {OH3} — fragment can flip")
    return problems


def add_shuttle(insert: str) -> str:
    """Wrap a p003-format insert in the two inward-pointing Esp3I sites."""
    if not insert.startswith(KOZAK):
        raise RuntimeError("insert does not start with the p003 Kozak block")
    if not insert.endswith(CS1):
        raise RuntimeError("insert does not end with the 10x capture sequence")
    return UPSTREAM + insert + DOWNSTREAM


def digest(shuttled: str) -> str:
    """Simulate Esp3I digestion and return the released fragment (top strand, from the
    5' overhang through the 3' overhang inclusive) — what lands in the new backbone."""
    sites = seqtools.find_sites(shuttled, ["Esp3I"])
    if len(sites) != 2:
        raise RuntimeError(f"expected exactly 2 Esp3I sites, found {len(sites)}: {sites}")
    (_, strand5, pos5), (_, strand3, pos3) = sites
    if strand5 != "+" or strand3 != "-":
        raise RuntimeError(
            f"sites must point inward: upstream on +, downstream on -; got {strand5}/{strand3}"
        )
    # + site: CGTCTC at pos5, spacer at pos5+6, overhang starts pos5+7
    start = pos5 + len(ESP3I) + 1
    # - site: GAGACG at pos3; spacer at pos3-1; overhang ends at pos3-2
    end = pos3 - 1
    return shuttled[start:end]


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", default="constructs_nterm.fasta")
    parser.add_argument("--out-dir", type=Path, default=HERE)
    args = parser.parse_args(argv)

    problems = check_overhangs()
    print(f"OH5 = {OH5} (Kozak[:4])   OH3 = {OH3} (10x CS1 + {PAD!r} pad)")
    print("NEB Esp3I 37-16 cycling: 100% predicted fidelity, 0 mismatch counts")
    print("overhang checks:", "; ".join(problems) if problems else "pass")
    if problems:
        return 1

    inserts = read_fasta(HERE / args.source)
    rows, out = [], []
    for gene, insert in inserts.items():
        shuttled = add_shuttle(insert)
        released = digest(shuttled)
        if released != insert + PAD:
            raise RuntimeError(f"{gene}: digest does not regenerate insert+pad exactly")
        sites = {e: seqtools.find_sites(shuttled, [e]) for e in ENZYMES}
        if len(sites["Esp3I"]) != 2 or sites["BbsI"] or sites["BsaI"]:
            raise RuntimeError(f"{gene}: unexpected Type IIS content {sites}")
        rows.append(
            {
                "gene": gene,
                "insert_nt": len(insert),
                "shuttled_nt": len(shuttled),
                "added_nt": len(shuttled) - len(insert),
                "esp3i_sites": len(sites["Esp3I"]),
                "bbsi_sites": len(sites["BbsI"]),
                "bsai_sites": len(sites["BsaI"]),
                "oh5": OH5,
                "oh3": OH3,
                "digest_regenerates_insert": True,
            }
        )
        out.append((gene, shuttled))
        print(
            f"{gene:9s} {len(insert):5d} -> {len(shuttled):5d} nt (+{len(shuttled) - len(insert)})  "
            f"Esp3I={len(sites['Esp3I'])} BbsI={len(sites['BbsI'])} BsaI={len(sites['BsaI'])}  "
            f"digest==insert ✓"
        )

    stem = Path(args.source).stem
    table = args.out_dir / f"{stem}_shuttle.tsv"
    with table.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
    fasta = args.out_dir / f"{stem}_shuttle.fasta"
    fasta.write_text(
        "".join(
            f">{g}|shuttle|Esp3I-flanked|{len(s)}nt\n"
            + "\n".join(s[i : i + 60] for i in range(0, len(s), 60))
            + "\n"
            for g, s in out
        )
    )
    print(f"\nWrote {table} and {fasta}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
