"""The puro-only intermediate: p003 with mScarlet3 and P2A deleted.

**No P-number.** This held p0036 until 2026-10-18. It was never ordered and never
registered (the Reagent Inventory has no P-0036), and ``empty_backbone.py`` superseded it
by deleting the insert as well, so one synthesised plasmid serves every ORF. Per the repo
convention a designed-but-parked construct is named for what it is, so the number went
back to the pool. Kept because p0039 is derived from this deletion and the reasoning below
is why the IRES coupling survives it.

All of this experiment's constructs go into a puro-only backbone to mimic the library.
mScarlet3 sits in the backbone, not in any insert, so this is a destination-vector build
and none of the insert designs change.

**IRES2 → PuroR, not SV40 → PuroR.** P008 already provides a puro-only vector, but it
drives PuroR from its own SV40 promoter as a separate transcription unit, which
decouples selection from ORF expression — a cell can survive puromycin while
transcribing the ORF poorly. p003 couples them through the IRES, and this keeps that.

The edit is a single clean deletion::

    p003        …IRES2 · linker·ATG | mScarlet3 · linker · P2A | PuroR · TAA · WPRE…
                                     └───── delete 7969..8721 ─────┘
    (this)      …IRES2 · linker·ATG | PuroR · TAA · WPRE…

**No initiator needs adding.** PuroR in p003 begins ``ACC GAG TAC`` (T-E-Y) — it has no
ATG of its own, because it is translated as a P2A downstream product. But the IRES
linker already carries an in-frame ``ATG`` at 7966..7968 (it was mScarlet3's upstream
Met). Deleting exactly 7969..8721 joins that ATG straight onto PuroR's second codon and
regenerates the **native** PuroR protein, ``MTEYKPTVRLATRDDVPRAVR…TWCMTRKPGA``.

Because nothing 5' of 7969 is touched, the IRES2-to-initiator spacing — the parameter
that sets IRES initiation efficiency — is bit-for-bit what p003 already validated.

⚠️ Input is the local-only p003 GenBank (same gap flagged in s1e2/enzyme_scan.py).
"""

from __future__ import annotations

import argparse
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

from deliverome_analysis import seqtools

HERE = Path(__file__).parent
P003 = (
    Path.home()
    / "Downloads"
    / "p003_pag-lenti_236083_esp3i_domesticated_efs_tfrc_bc_alfa_iresfragmid_wpre_puro_assembled.gb"
)

DELETE = (7969, 8721)  # mScarlet3 + linker + P2A, 1-based inclusive
ATG_AT = 7966  # the in-frame initiator already present in the IRES linker
NATIVE_PURO_START = "MTEYKPTVRLATRDDVPRAVR"
NATIVE_PURO_END = "TWCMTRKPGA"


def build(p003_path: Path = P003) -> tuple[str, str]:
    if not p003_path.exists():
        raise SystemExit(f"p003 GenBank not found: {p003_path}")
    parent = str(SeqIO.read(p003_path, "genbank").seq).upper()
    backbone = parent[: DELETE[0] - 1] + parent[DELETE[1] :]

    orf = backbone[ATG_AT - 1 : ATG_AT - 1 + 600]
    protein = str(Seq(orf).translate())
    if orf[:3] != "ATG":
        raise RuntimeError("the IRES linker ATG is not where expected")
    if not protein.startswith(NATIVE_PURO_START):
        raise RuntimeError(f"PuroR N-terminus is not native: {protein[:24]}")
    if not protein.rstrip("*").endswith(NATIVE_PURO_END):
        raise RuntimeError(f"PuroR C-terminus is not native: {protein[-14:]}")
    if protein.count("*") != 1 or not protein.endswith("*"):
        raise RuntimeError("PuroR ORF is not a single clean reading frame")
    sites = seqtools.find_sites(
        backbone[DELETE[0] - 70 : DELETE[0] + 90], ["BbsI", "Esp3I", "BsaI"]
    )
    if sites:
        raise RuntimeError(f"deletion created Type IIS site(s): {sites}")
    # everything 5' of the deletion is untouched -> IRES spacing preserved
    if backbone[: DELETE[0] - 1] != parent[: DELETE[0] - 1]:
        raise RuntimeError("sequence 5' of the deletion changed")
    return parent, backbone


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=HERE)
    args = ap.parse_args(argv)
    parent, backbone = build()

    print(f"p003    {len(parent):>7,} bp")
    print(
        f"delete  {DELETE[0]}..{DELETE[1]}  ({DELETE[1] - DELETE[0] + 1} nt: "
        f"mScarlet3 + linker + P2A)"
    )
    print(f"puro-only  {len(backbone):>7,} bp  (-{len(parent) - len(backbone)})\n")
    print("new junction (IRES2 linker -> PuroR):")
    print(f"  …{backbone[7932:7968]}|{backbone[7968:7995]}…")
    orf = backbone[ATG_AT - 1 : ATG_AT - 1 + 600]
    print(
        f"\nPuroR: {len(orf)} nt, native protein, "
        f"{str(Seq(orf).translate())[:16]}…{str(Seq(orf).translate())[-11:]}"
    )
    print(
        "all checks pass: native N- and C-termini, single reading frame, "
        "no Type IIS created, IRES spacing untouched"
    )

    out = args.out_dir / "puro_only_backbone.fasta"
    out.write_text(
        f">puro_only_backbone|puro-only intermediate|p003 minus mScarlet3+P2A "
        f"(del {DELETE[0]}..{DELETE[1]})|{len(backbone)}bp\n"
        + "\n".join(backbone[i : i + 60] for i in range(0, len(backbone), 60))
        + "\n"
    )
    print(f"\nWrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
