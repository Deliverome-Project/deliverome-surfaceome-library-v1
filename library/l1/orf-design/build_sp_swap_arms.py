"""Signal-peptide SWAP arms: the Igk leader replacing a member's own weak signal peptide.

Distinct from both neighbouring builders, and the distinction is the point:

* ``rebuild_disorder_anchored.py`` keeps the native signal peptide and puts the cassette
  a few residues downstream of its cleavage site.
* ``build_gpcr_rescue_arms.py`` **prepends** the leader to a receptor that has no signal
  peptide at all, so the receptor keeps its own Met1 and nothing is removed.
* this module **replaces** residues ``1..sp_end`` with the leader plus its Asp, so the
  member loses its own signal peptide and gains a defined one.

Only a member whose native signal peptide fails re-prediction once tagged earns an arm
here. The swap is not an improvement offered to a working design; it is the repair for a
design that the predictors say does not express.

**VSIG10 is the only such member in the batch.** Re-predicting all eight native-SP arms
with their cassettes in place (``canary_rp_1792371023``, SignalP 6.0+h) preserved the
cleavage site on six, moved HEG1's by one residue, and ABOLISHED VSIG10's: the tagged
design is called OTHER with p(SP) 0.000137 and no cleavage site at all, against 0.999832
for the untagged protein. Junction 32 -- the only other residue in its 30..32 admissible
interval -- also fails when MODELLED, at p(SP) 0.0 (``canary_rp_1792371276``).

Those two results are not equally sourced and the distinction matters when quoting them.
The junction-31 result is on the REAL construct: its tagged sequence comes from
``processed://constructs/designed/p0026.gb``. The junction-32 sequence was assembled for
that test and exists nowhere else -- verified residue-exact here (sha256 0df6fab8bb253e86
equals native-with-cassette-after-32, reconstructed independently), but never designed,
built or registered. So the honest claim is "the built construct fails, and the one
alternative position its interval allows also fails when modelled" -- which is still
enough to say not to expect a one-residue fix, but is not "both positions tested".

The submitted sequences are archived at
``processed://sprint4/s4e9-orf-design-fg-panel/reprediction/`` -- one
``<run_id>/submission.json`` per run plus a ``provenance.json``. Cite those keys rather
than the Modal run ids: a Modal Volume is working storage and the published runs recorded
the shas without the sequences, so for a while it was possible to confirm THAT sixteen
sequences were predicted without being able to recover WHICH. Verified here on
2026-10-19: 18 records, zero sha mismatches; all eight tagged sequences reconstruct
exactly as native-with-cassette-at-the-stated-residue; and the junction-31 VSIG10 record
is byte-identical to p0026's own protein (sha256 1bfd4c87bedd29ce), which is the link
that makes this a prediction about the construct rather than about something like it.

The cause is visible in the native prediction rather than the tagged one. VSIG10's own
signal peptide is the weakest-DEFINED in the set: SignalP is certain one exists
(p(SP) 0.999832) but not where it ends, giving the cleavage site only **0.348** -- the
11th lowest of 1,074 L1 members with a signal peptide, first percentile. DeepTMHMM v1 and
SignalP 6 disagree about it too (1-27 against 1-30). A marginal c-region plus 39 residues
of glycine-rich cassette inside SignalP's 70-residue window is what tips the call.

**The leader fixes it, measured rather than assumed** (``canary_rp_1792372578``):

    arm                         call   cleavage   p(SP)      cleavage site p
    VSIG10 native, untagged     SP           30   0.999832             0.348
    VSIG10 native, tagged       OTHER      none   0.000137              none
    VSIG10 igk,    untagged     SP           20   0.999844             0.980
    VSIG10 igk,    tagged       SP           20   0.999833             0.982

So the swap does not merely restore the call, it produces a BETTER-defined cleavage site
than VSIG10 ever had natively (0.982 against 0.348), and the cassette costs it nothing.
Both arms are built: the native arm is the comparator, and it ships as a PREDICTED
FAILURE rather than as an oversight -- if it expresses after all, SignalP was wrong about
it, which is worth one well to learn. Becca, 2026-10-18.

The mature chains differ by exactly one residue, so the pair is a clean leader contrast:

    native arm   ...EA | V-GGGGS-SpyTag003-...      cleaves after 30, mature starts V31
    igk arm      ...TG | D-V-GGGGS-SpyTag003-...    cleaves after 20, mature starts D

    uv run --locked python build_sp_swap_arms.py
    uv run --locked python build_sp_swap_arms.py --write
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import build_constructs as bc
from build_gpcr_rescue_arms import LEADER, LEADER_AA, PLUS1, PLUS1_AA

from deliverome_analysis import seqtools

HERE = Path(__file__).parent
SEED = 20261018

#: gene -> the swap. ``sp_end`` is the rule this folder already uses for the native arm
#: (``rebuild_disorder_anchored.JUNCTIONS[gene]["sp_end"]``, SignalP 6); ``junction`` is
#: that arm's cassette position, reused unchanged so the two arms differ ONLY in the
#: leader. Keeping the junction fixed is what makes the pair readable.
SWAPS = {
    "VSIG10": {
        "sp_end": 30,
        "junction": 31,
        "native_cleavage_p": 0.3477,
        "igk_cleavage_p": 0.9825,
        "canary": "canary_rp_1792372578",
        "why": (
            "its native signal peptide is abolished by the cassette (tagged p(SP) "
            "0.000137, called OTHER, no cleavage site; canary_rp_1792371023) at both "
            "residues of its admissible interval. The native cleavage site was only "
            "0.348-confident even untagged -- first percentile of the 1,074 L1 members "
            "with a signal peptide -- and the Igk leader replaces it with a 0.982 site "
            "that survives the cassette intact."
        ),
    },
}


def arm_names() -> list[str]:
    return [f"{gene}_igk" for gene in SWAPS]


def sp_swap_flanks(cds: dict[str, str], flank: int = 150) -> dict[str, tuple[str, str]]:
    """Each swap arm's sequence either side of its cassette, for the encoding draw.

    Same contract as ``build_gpcr_rescue_arms.arm_flanks``: the diversifier chooses a
    synonymous cassette encoding against these, so the 5' flank must already carry the
    leader rather than the signal peptide it replaces.
    """
    out: dict[str, tuple[str, str]] = {}
    for gene, spec in SWAPS.items():
        if gene not in cds:
            continue
        body = cds[gene][:-3]
        head = body[spec["sp_end"] * 3 : spec["junction"] * 3]
        upstream = bc.KOZAK + LEADER + PLUS1 + head
        out[f"{gene}_igk"] = (upstream[-flank:], body[spec["junction"] * 3 :][:flank])
    return out


def build(gene: str, cds: str, cassette_nt: str, barcode: str) -> dict:
    """Kozak + leader + Asp + mature-with-cassette + stop + barcode amplicon."""
    spec = SWAPS[gene]
    sp_end, junction = spec["sp_end"], spec["junction"]
    body = cds[:-3]  # drop the native stop; STOP is appended below
    native = seqtools.translate(cds).rstrip("*")
    if junction <= sp_end:
        raise SystemExit(f"{gene}: junction {junction} is inside the signal peptide")
    if junction * 3 > len(body):
        raise SystemExit(f"{gene}: junction {junction} is past the ORF")

    head = body[sp_end * 3 : junction * 3]  # mature residues kept ahead of the cassette
    coding = LEADER + PLUS1 + head + cassette_nt + body[junction * 3 :]

    tamed, edits = seqtools.domesticate(coding, bc.ENZYMES, seed=SEED)
    expect = LEADER_AA + PLUS1_AA + native[sp_end:junction] + bc.CASSETTE_AA + native[junction:]
    got = seqtools.translate(tamed)
    if got != expect:
        raise RuntimeError(f"{gene}_igk: protein changed\n  want {expect[:70]}\n  got  {got[:70]}")
    # The published leader and the drawn cassette encoding are both fixed inputs: an
    # edit inside either is a silent substitution of something a record points at.
    if tamed[: len(LEADER)] != LEADER:
        raise RuntimeError(f"{gene}_igk: domestication edited the published Igk leader")
    if cassette_nt not in tamed:
        raise RuntimeError(f"{gene}_igk: domestication edited the drawn cassette encoding")

    insert = (
        bc.KOZAK + tamed + bc.STOP + bc.SPACER + bc.FWD_PRIMER + barcode + bc.REV_PRIMER + bc.CS1
    )
    protein = seqtools.translate(tamed + bc.STOP)
    if protein[-1] != "*" or "*" in protein[:-1]:
        raise SystemExit(f"{gene}_igk: internal stop or missing stop")

    sites = seqtools.find_sites(insert, bc.ENZYMES)
    cut = len(LEADER) + len(PLUS1) + len(head)
    seams = (cut, cut + len(cassette_nt))
    at_seam = [s for s in sites if any(c - 7 < s[2] - len(bc.KOZAK) < c + 1 for c in seams)]

    # Located in the insert THIS call assembled, never restated from a table or read back
    # from a FASTA the same run rewrites.
    i = insert.find(cassette_nt)
    if i < 0:
        raise SystemExit(f"{gene}_igk: cassette not found in the assembled insert")
    return {
        "construct": f"{gene}_igk",
        "gene": gene,
        "arm": "igk_leader_swap",
        "insert_after_residue": (i - len(bc.KOZAK)) // 3,
        "junction_source": "located in the insert assembled by this call",
        "native_sp_removed": f"1-{sp_end}",
        "native_junction_residue": junction,
        "leader_aa": len(LEADER_AA),
        "mature_starts": PLUS1_AA + native[sp_end:junction],
        "native_cleavage_p": spec["native_cleavage_p"],
        "igk_cleavage_p": spec["igk_cleavage_p"],
        "reprediction_run": spec["canary"],
        "insert_nt": len(insert),
        "protein_aa": len(protein) - 1,
        "domestication_edits": len(edits),
        "type_iis_sites": len(sites),
        "sites_at_a_seam": len(at_seam),
        "barcode": barcode,
        "why": spec["why"],
        "insert": insert,
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--write", action="store_true")
    args = ap.parse_args(argv)

    cds = bc.read_fasta(HERE / "orf_domesticated.fasta")
    cassettes = bc.load_cassettes()
    table = HERE / "barcodes_sp_swap.tsv"
    if not table.exists():
        raise SystemExit(f"{table.name} missing — reserve a barcode first (reserve_barcodes.py)")
    barcodes = {r["construct"]: r["barcode"] for r in csv.DictReader(table.open(), delimiter="\t")}

    rows = []
    for gene in sorted(SWAPS):
        name = f"{gene}_igk"
        if name not in cassettes:
            raise SystemExit(f"{name}: no cassette encoding — run diversify_cassettes.py")
        if name not in barcodes:
            raise SystemExit(f"{name}: no barcode in {table.name}")
        rows.append(build(gene, cds[gene], cassettes[name], barcodes[name]))

    for r in rows:
        print(
            f"{r['construct']:14} removes native SP {r['native_sp_removed']}, "
            f"cassette after residue {r['insert_after_residue']} of the leader-bearing ORF, "
            f"mature starts {r['mature_starts']!r}, {r['insert_nt']} nt, "
            f"{r['protein_aa']} aa, {r['domestication_edits']} edit(s), "
            f"{r['sites_at_a_seam']} seam site(s)"
        )

    if not args.write:
        print("\n(dry run) nothing written; pass --write")
        return 0

    cols = [k for k in rows[0] if k != "insert"]
    with (HERE / "constructs_sp_swap.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({k: r[k] for k in cols})
    with (HERE / "constructs_sp_swap.fasta").open("w") as fh:
        for r in rows:
            fh.write(f">{r['construct']}|bc={r['barcode']}|{r['insert_nt']}nt\n{r['insert']}\n")
    print(f"\nWrote constructs_sp_swap.tsv and constructs_sp_swap.fasta ({len(rows)} arm(s))")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
