"""Reserve barcodes for this experiment's constructs from the running P-series registry.

The registry is the single source of truth for which pool barcodes are spoken for:

    processed://sprint1/s1e2-fragmid-vector-design/barcode_assignments.tsv

Its own header says to reserve from the pool and append there *before* assigning a
barcode to a new construct. This script does that, and refuses to hand out anything
already listed.

Two screens beyond "not already used", both inferred from why indices 2 and 7 were
skipped ("boundary-enzyme or splice-donor motif"):

* **boundary enzyme** — no BbsI / Esp3I / BsaI site anywhere in the barcode read in its
  real junction context (forward primer site + barcode + reverse primer site + the
  Esp3I shuttle that follows). A barcode clean in isolation can still create a site
  across a junction.

  These flanks are **fixed across the whole P-series**, not specific to this experiment:
  the barcode sits between the same forward and reverse primer sites in P003, P008-P013
  and p0016-p0021, verified against each construct's map. So one junction screen covers
  every construct in the series and no per-construct flank configuration is needed. The
  Esp3I shuttle is new here and absent from the earlier constructs, but it changes
  nothing — the pool yields the same 94 junction failures with or without it.
* **splice donor** — no ``GT[AG]AG`` in that same context. Pool index 7 is
  ``GGTAAGACGCGCCTGAGTAG``, which opens with ``GGTAAG``, the canonical 5' splice-site
  consensus; in a lentiviral transcript a cryptic donor here would splice out
  everything downstream. The s1e2 generator's filter stack does **not** screen for this
  (it covers GC, homopolymer, ATG, polyA, hairpin, miRNA seed, Type IIS and distance),
  so it has to be applied at reservation time.

Nothing is written to S3 by this script — it prints the rows to append and writes them
beside the local barcode tables. Appending to the shared registry is a separate,
deliberate step.
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

from build_constructs import CS1, ENZYMES, FWD_PRIMER, REV_PRIMER

from deliverome_analysis import seqtools
from deliverome_analysis.s3 import processed

HERE = Path(__file__).parent
REGISTRY_KEY = "sprint1/s1e2-fragmid-vector-design/barcode_assignments.tsv"
POOL_KEY = "sprint1/s1e2-fragmid-vector-design/barcodes_20nt_30k.tsv"

SPLICE_DONOR = re.compile(r"GT[AG]AG")
# Trailing context a barcode actually sits in: reverse primer site, then the AG pad and
# the inward-pointing Esp3I shuttle site (see shuttle.py).
TRAILING = REV_PRIMER + CS1 + "AG" + "A" + "GAGACG"

# Constructs needing a barcode, in the order they get P-numbers.
CONSTRUCTS = [
    ("HEG1", "nterm", "barcodes.tsv", "HEG1"),
    ("ERBB2", "nterm", "barcodes.tsv", "ERBB2"),
    ("AXL", "nterm", "barcodes.tsv", "AXL"),
    ("SLC9A3", "nterm", "barcodes.tsv", "SLC9A3"),
    ("VSIG10", "nterm", "barcodes.tsv", "VSIG10"),
    ("NECTIN4", "nterm", "barcodes.tsv", "NECTIN4"),
    ("FOLR1", "nterm", "barcodes.tsv", "FOLR1"),
    ("PTPRC", "nterm", "barcodes.tsv", "PTPRC"),
    # Loop sites, keyed to match ``barcodes_loop.tsv``'s own (gene, variant) exactly.
    #
    # A stale entry here does not fail loudly: the key simply matches no row, the arm looks
    # unreserved, and the allocator draws it a barcode and a P-number that a BUILT
    # construct already holds. That is not hypothetical. This list still said
    # ``MS4A1|P79``, ``MS4A1|I162`` and ``CCR5|S271`` on 2026-10-18 -- names retired when
    # the sites were repointed to F146, Y161 and F263 -- and a dry run accordingly offered
    # to issue p0048-p0050 to three constructs that do not exist, handing one of them
    # ``GCATTCTAGCATCTTCGTCC`` at index 27, the barcode CCR5 F263 is built with. Nothing
    # shipped wrong because the allocator had not been re-run since the repointing; the
    # collision was latent, and the dry-run flag is what surfaced it.
    #
    # The retired names are SLC2A4 P66/E193, MS4A1 Q187/N171, CCR5 S271 and MS4A1 P79/I162.
    # Repoint this list in the same commit that repoints a site.
    ("SLC2A4", "loop E67", "barcodes_loop.tsv", "SLC2A4|E67"),
    ("SLC2A4", "loop V319", "barcodes_loop.tsv", "SLC2A4|V319"),
    ("MS4A1", "loop F146", "barcodes_loop.tsv", "MS4A1|F146"),
    ("MS4A1", "loop Y161", "barcodes_loop.tsv", "MS4A1|Y161"),
    ("CCR5", "loop F263", "barcodes_loop.tsv", "CCR5|F263"),
    ("AVPR2", "loop L302", "barcodes_loop.tsv", "AVPR2|L302"),
    ("TFRC", "cterm", "barcodes_tfrc.tsv", "TFRC|cterm"),
    # The GPCR rescue arms, in the order barcodes_gpcr.tsv holds them. Three genes carry
    # a leader-free/igk pair; ACKR3 and AVPR2 are leader-free only, so a failure on
    # either cannot be split between cassette load and the receptor -- the one design
    # asymmetry left in the batch, kept because the batch is capped at 23.
    ("CNR1", "leader_free", "barcodes_gpcr.tsv", "CNR1_plain"),
    ("CNR1", "igk_leader", "barcodes_gpcr.tsv", "CNR1_igk"),
    ("CCR5", "leader_free", "barcodes_gpcr.tsv", "CCR5_plain"),
    ("CCR5", "igk_leader", "barcodes_gpcr.tsv", "CCR5_igk"),
    ("OR51E2", "leader_free", "barcodes_gpcr.tsv", "OR51E2_plain"),
    ("OR51E2", "igk_leader", "barcodes_gpcr.tsv", "OR51E2_igk"),
    ("ACKR3", "leader_free", "barcodes_gpcr.tsv", "ACKR3_plain"),
    ("AVPR2", "leader_free", "barcodes_gpcr.tsv", "AVPR2_plain"),
    # The signal-peptide SWAP arm. VSIG10's native arm is a predicted non-expressor
    # once tagged -- SignalP 6 calls the tagged design OTHER with no cleavage site at
    # either residue of its admissible interval -- so it gains an Igk twin in which the
    # leader REPLACES the failing signal peptide rather than preceding it. Measured to
    # restore a 0.982 cleavage site with the cassette in place (canary_rp_1792372578).
    # Architecturally its own layout, hence its own table; see build_sp_swap_arms.py.
    ("VSIG10", "igk_leader_swap", "barcodes_sp_swap.tsv", "VSIG10_igk"),
    # AVPR2_igk, ACKR3_igk and the 19-aa AVPR2_plain_alfa arm are NOT here: all three
    # were cut on 2026-10-18 to hold the batch at 23. Every remaining arm carries the
    # full 51-aa cassette, so no arm varies cassette length against position.
    # TFRC linker_3x/4x are NOT here: they were cut from the batch, parked under
    # descriptive names in constructs/scratch/, and released their barcodes to the GPCR
    # arms. Reviving one needs a fresh draw and a new row here.
]
# P-numbers this script cannot learn from the barcode registry, and must never issue.
#
# The registration authority for a P-number is the **Reagent Inventory** sheet (which
# syncs to Notion), not this registry: the registry only knows constructs that drew a
# barcode from the s1e2 pool. So a vector can be fully registered and still be invisible
# here. That is how the 2026-10-18 renumbering put two ordered constructs onto P-0039 and
# P-0040 -- both registered, neither barcoded from the pool -- and a third onto p0036, a
# backbone this folder builds. Numbers are cheap; a number meaning two plasmids is not.
#
# Re-read the sheet's Plasmid tab before any renumbering and update this set.
# Reagent Inventory, Plasmid tab, read 2026-10-18. P-0005 and P-0010 are gaps in that
# sheet; left alone rather than recycled, since nothing records why they were skipped.
REGISTERED_PLASMIDS = (
    {1, 2, 3, 4, 6, 7, 8, 9}
    | set(range(11, 22))  # the miR-L iFFL series p0016-p0021 among them
    | {
        39,  # pAG Lenti insert-agnostic library backbone (empty_backbone.py)
        40,  # miR-L iFFL PGK intron-split PuroR, TACC target (assign_p0040_barcode.py)
    }
)
#: Built in this folder and barcode-free, so the registry cannot report them either.
#
# 36 is NOT here. It was listed until 2026-10-18 for the puro-only intermediate, but that
# plasmid was never ordered and never registered -- the Reagent Inventory has no P-0036 --
# and p0039 superseded it by removing the insert as well. Per the repo convention a
# designed-but-parked construct is named for what it is, never with a P-number, so the
# number is released. The artifact is now puro_only_backbone.fasta.
NON_BARCODED_PLASMIDS = {
    39,  # p0039, the insert-agnostic destination vector (empty_backbone.py)
}
#: Everything the allocator must skip over.
UNAVAILABLE_PLASMIDS = REGISTERED_PLASMIDS | NON_BARCODED_PLASMIDS


def read_registry() -> tuple[set[str], set[int], list[str], set[tuple[str, str]], set[int]]:
    """Barcodes, pool indices, skipped indices, already-reserved (gene, arm), used P-numbers.

    The last two are what make this script safe to re-run. Without them it
    re-draws every entry in CONSTRUCTS from a fixed starting P-number, which
    silently reassigns barcodes that are already registered — and, worse, looks
    like a normal successful run while doing it.
    """
    text = Path(processed.download(REGISTRY_KEY)).read_text()
    lines = text.splitlines()
    body = [ln for ln in lines if ln and not ln.startswith("#")]
    rows = list(csv.DictReader(body, delimiter="\t"))
    used_bc = {r["barcode"] for r in rows if r["barcode"]}
    used_ix = {int(r["pool_index"]) for r in rows if r["pool_index"]}
    done = {(r["orf"], r["notes"].split(";")[0].strip()) for r in rows if r["orf"]}
    used_p = set(NON_BARCODED_PLASMIDS)
    for r in rows:
        m = re.fullmatch(r"p(\d+)", r["construct"].strip())
        if m:
            used_p.add(int(m.group(1)))
    skipped = set()
    for ln in lines:
        if ln.startswith("#") and "skipped" in ln:
            skipped |= {int(n) for n in re.findall(r"\b(\d+)\b", ln)}
    return used_bc, used_ix, sorted(skipped), done, used_p


def screens_pass(bc: str) -> tuple[bool, str]:
    context = FWD_PRIMER + bc + TRAILING
    sites = seqtools.find_sites(context, ENZYMES)
    # the shuttle's own Esp3I site is expected; anything else is a boundary-enzyme hit
    extra = [s for s in sites if s[0] != "Esp3I"] + [
        s for s in sites if s[0] == "Esp3I" and s[2] < len(FWD_PRIMER) + len(bc)
    ]
    if extra:
        return False, f"boundary enzyme {extra}"
    if SPLICE_DONOR.search(context):
        return False, f"splice donor {SPLICE_DONOR.search(context).group(0)}"
    return True, ""


def _locally_used_plasmids() -> set[int]:
    """P-numbers already spoken for in this folder's barcode tables."""
    out: set[int] = set()
    for table in {t for _g, _a, t, _k in CONSTRUCTS}:
        path = HERE / table
        if not path.exists():
            continue
        with path.open() as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                m = re.fullmatch(r"p(\d{3,4})", (row.get("plasmid") or "").strip())
                if m:
                    out.add(int(m.group(1)))
    return out


def _locally_assigned() -> tuple[set[tuple[str, str]], set[str]]:
    """((gene, arm) already assigned here, the barcodes they hold).

    Read from the per-construct tables rather than the registry, because those are what
    the built maps were made from.
    """
    done: set[tuple[str, str]] = set()
    bcs: set[str] = set()
    for gene, arm, table, key in CONSTRUCTS:
        path = HERE / table
        if not path.exists():
            continue
        with path.open() as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                bc = (row.get("barcode") or "").strip().upper()
                if not bc or set(bc) - set("ACGT"):
                    continue
                ident = row.get("construct") or "|".join(
                    x for x in (row.get("gene"), row.get("variant") or row.get("arm")) if x
                )
                if ident == key or row.get("gene") == key:
                    done.add((gene, arm))
                    bcs.add(bc)
    return done, bcs


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=HERE)
    ap.add_argument(
        "--dry-run",
        action="store_true",
        help="draw and screen but write nothing. This allocator was one-shot with no way "
        "to look first, which is why it was never safe to run against a stale target list.",
    )
    args = ap.parse_args(argv)

    used_bc, used_ix, skipped, done, used_p = read_registry()
    # A construct already assigned in one of this folder's barcode tables is DONE, even
    # when the shared registry has never heard of it. The registry is 16 reservations
    # behind reality -- they were never appended -- so trusting it alone makes every
    # existing construct look unreserved and draws it a second, different barcode. The
    # built GenBank maps already contain the first one, so that would silently split each
    # construct's identity between its map and its table.
    local, local_bc = _locally_assigned()
    done |= local
    used_bc |= local_bc
    print(
        f"registry: {len(used_bc)} barcode(s) in use, pool indices used {sorted(used_ix)}, "
        f"skipped {skipped}, P-numbers used up to p{max(used_p):04d}"
    )

    wanted = [c for c in CONSTRUCTS if (c[0], c[1]) not in done]
    for gene, arm, _t, _k in CONSTRUCTS:
        if (gene, arm) in done:
            print(f"  already reserved, leaving alone: {gene} {arm}")
    if not wanted:
        print("\nevery construct already has a reservation — nothing to draw")
        return 0

    # Key by the file's own `index` column, never by row position. The pool is a
    # downselection: it holds 28,506 rows whose indices run to 30,000, so position and
    # index diverge at the first gap. Walking positions while vetoing against `used_ix`
    # -- which holds real indices from the registry -- mixes the two namespaces, and
    # records an index that points at a different barcode than the one handed out.
    # That is how p0037 and p0038 came to claim indices 23 and 24 on 2026-10-06 when
    # their barcodes are at 25 and 26; index 23 does not exist in the pool at all.
    pool = {
        int(r["index"]): r["barcode"].strip().upper()
        for r in csv.DictReader(Path(processed.download(POOL_KEY)).open(), delimiter="\t")
        if r.get("index")
    }

    reserved, rejected = [], []
    for idx in sorted(pool):
        if len(reserved) == len(wanted):
            break
        if idx in used_ix or idx in skipped:
            continue
        bc = pool[idx]
        if bc in used_bc:
            continue
        ok, why = screens_pass(bc)
        if not ok:
            rejected.append((idx, bc, why))
            continue
        reserved.append((idx, bc))
    if len(reserved) < len(wanted):
        raise RuntimeError("pool exhausted before every construct had a barcode")
    for idx, bc in reserved:
        if pool[idx] != bc:
            raise RuntimeError(f"index {idx} does not hold {bc}; namespaces have diverged")

    if rejected:
        print("\nrejected by screen:")
        for idx, bc, why in rejected:
            print(f"  index {idx:3d} {bc}  — {why}")

    print(f"\n{'plasmid':9s} {'gene':9s} {'arm':11s} {'idx':>4s}  barcode")
    rows = []
    # The P-number must clear LOCAL assignments too, not just the registry's. The
    # registry knows nothing of p0041-p0049 -- those reservations were never appended --
    # so deriving the next number from it alone handed out p0041 a second time, to a
    # different construct. Same root cause as the barcode side above, and the same fix:
    # the registry is a lagging view, the local tables are the current one.
    local_p = _locally_used_plasmids()
    taken = used_p | local_p | UNAVAILABLE_PLASMIDS
    next_p = max(taken) + 1
    for (gene, arm, _table, _key), (idx, bc) in zip(wanted, reserved, strict=True):
        # Step over anything already spoken for rather than trusting max()+1 to clear it.
        # max()+1 only works while every taken number sits below the high-water mark; a
        # number registered in the sheet *above* it would be issued silently.
        while next_p in taken:
            next_p += 1
        plasmid = f"p{next_p:04d}"
        taken.add(next_p)
        next_p += 1
        print(f"{plasmid:9s} {gene:9s} {arm:11s} {idx:4d}  {bc}")
        rows.append(
            {
                "barcode": bc,
                "pool_index": idx,
                "construct": plasmid,
                "orf": gene,
                "position": "",
                "assigned_by": "s1e2 pool",
                "status": "designed",
                "notes": f"{arm}; analysis/sprint4/s4e9-orf-design-fg-panel; "
                f"provisional P-number, not yet built",
            }
        )

    out = args.out_dir / "barcode_reservations.tsv"
    if args.dry_run:
        print(f"\n(dry run) {len(rows)} row(s) drawn and screened; {out.name} NOT written")
        return 0
    with out.open("w", newline="") as h:
        w = csv.DictWriter(h, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(f"\nWrote {out} ({len(rows)} new row(s)) — register them with append_reservations.py,")
    print("which is the step that writes to the shared S3 registry.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
