"""Build truncated, surface-displayable constructs for the oversized genes.

The library caps ORFs for lentiviral packaging, which drops 107 genes including
most of the delivery-relevant giants. s1e3 gathered per-gene literature
precedents for expressing a shortened form of each; this turns those precedents
into actual sequences.

Every cut is anchored on UniProt coordinates rather than on a residue number
copied from a paper, and each construct is checked for the thing that makes
these work at all: the product must still carry a membrane anchor, and must
start with a signal sequence -- either the gene's own, or the Igk leader where
the fragment begins in the middle of the protein and has none.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

IGK = "METDTLLLWVLLLWVPGSTG"          # see igk_leader.tsv
MIN_ECTO_SPACER = 30                   # residues required between cleavage and the TM
FIXED_PAYLOAD = 4087                   # p0039 LTR-to-LTR with the slot removed
INSERT_OVERHEAD = 131                  # Kozak + primer sites + CS1 + Esp3I wrapper
TAG_CASSETTE = 153                     # the 51-residue ALFA/SpyTag cassette
BARCODE = 20
# 9 kb LTR-to-LTR, counted with the tag cassette and barcode in place rather than
# against the pre-tag intermediate -- the intermediate is not what gets packaged.
GENOME_CEILING = 9000


def segments(keep: str, rec: dict) -> list[tuple[int, int]]:
    """Resolve a spec's keep-rule into 1-based inclusive residue ranges."""
    feats = rec["features"]
    tm = next(f for f in feats if f["type"] == "Transmembrane")
    sig = next((f for f in feats if f["type"] == "Signal"), None)
    out: list[tuple[int, int]] = []
    for part in keep.split("+"):
        part = part.strip()
        if part == "signal":
            if sig:
                out.append((sig["start"], sig["end"]))
        elif part.startswith("last:"):
            n = int(part.split(":")[1])
            out.append((rec["length"] - n + 1, rec["length"]))
        elif part.startswith("tm:"):
            out.append((tm["start"], rec["length"]))
        elif ":" in part:
            a, b = part.split(":")
            out.append((int(a), rec["length"] if b == "end" else int(b)))
        else:
            raise ValueError(f"unparsed keep-rule part: {part!r}")
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--spec", default="giant_truncation_spec.tsv")
    ap.add_argument("--features", default="giant_features.json")
    ap.add_argument("--out", default="giant_truncations.tsv")
    args = ap.parse_args()

    feats = json.loads(Path(args.features).read_text())
    rows, problems = [], []
    for spec in csv.DictReader(open(args.spec), delimiter="\t"):
        gene = spec["gene"]
        base = gene.split("_")[0]
        rec = feats.get(base)
        if rec is None:
            problems.append(f"{gene}: no UniProt record")
            continue
        segs = segments(spec["keep"], rec)
        prot = "".join(rec["sequence"][a - 1:b] for a, b in segs)
        if spec["leader"] == "igk":
            prot = IGK + prot

        tm = next(f for f in rec["features"] if f["type"] == "Transmembrane")
        has_anchor = any(a <= tm["start"] and b >= tm["end"] for a, b in segs)
        # Keeping residue 1 is not the same as keeping a signal peptide. SDK1 has
        # no annotated signal at all -- it opens with a 98-residue Gly-rich
        # extension and the mature protein starts at Q99 -- so a construct built
        # from residue 1 would have no leader while looking like it did.
        sig_feat = next((f for f in rec["features"] if f["type"] == "Signal"), None)
        starts_clean = spec["leader"] == "igk" or (
            sig_feat is not None
            and any(a <= sig_feat["start"] and b >= sig_feat["end"] for a, b in segs))
        if not has_anchor:
            problems.append(f"{gene}: the kept segments do not contain the transmembrane span")
        # An exogenous leader REPLACES the native signal; it does not stack in
        # front of it. Two tandem signal sequences would be cleaved after the
        # first, leaving the second as a hydrophobic mature N-terminus that can
        # act as an unintended anchor. PLXNB3 was built this way before this
        # check existed.
        if spec["leader"] == "igk" and sig_feat is not None and any(
                a <= sig_feat["start"] and b >= sig_feat["end"] for a, b in segs):
            problems.append(
                f"{gene}: carries the Igk leader AND keeps its native signal peptide - "
                f"two tandem signal sequences. Drop the signal segment from the keep rule.")

        if not starts_clean:
            problems.append(
                f"{gene}: no leader - it neither carries the Igk leader nor keeps an "
                f"annotated signal peptide"
                + ("" if sig_feat else " (and UniProt annotates none for this protein)"))

        # Carrying a TM is not the same as inserting into the membrane. The
        # signal sequence has to be cleaved before the TM can act as a
        # stop-transfer anchor, and a leader butted straight against the TM
        # leaves no cleavage site -- the TM gets read as a signal anchor and the
        # topology inverts. Require a real juxtamembrane stub between them.
        offset, tm_at = (len(IGK) if spec["leader"] == "igk" else 0), None
        for a, b in segs:
            if a <= tm["start"] <= b:
                tm_at = offset + (tm["start"] - a) + 1
                break
            offset += b - a + 1
        sig_len = len(IGK) if spec["leader"] == "igk" else next(
            (f["end"] for f in rec["features"] if f["type"] == "Signal"), 0)
        ecto = (tm_at - sig_len - 1) if tm_at else -1
        cyto = (len(prot) - tm_at - (tm["end"] - tm["start"]) ) if tm_at else -1
        if ecto < MIN_ECTO_SPACER:
            problems.append(
                f"{gene}: only {ecto} aa between signal cleavage and the TM; "
                f"needs at least {MIN_ECTO_SPACER}")

        orf_nt = len(prot) * 3 + 3
        genome = FIXED_PAYLOAD + orf_nt + INSERT_OVERHEAD
        genome_tagged = genome + TAG_CASSETTE + BARCODE
        if genome_tagged > GENOME_CEILING:
            problems.append(
                f"{gene}: {genome_tagged} bp LTR-to-LTR with the tag cassette and barcode, "
                f"over the {GENOME_CEILING} bp ceiling by {genome_tagged - GENOME_CEILING}")
        rows.append({
            "gene": gene,
            "accession": rec["accession"],
            "full_aa": rec["length"],
            "truncated_aa": len(prot),
            "shrink": f"{len(prot) / rec['length']:.1%}",
            "orf_nt": orf_nt,
            "ltr_to_ltr_pretag_bp": genome,
            "ltr_to_ltr_tagged_bp": genome_tagged,
            "within_9kb": genome_tagged <= GENOME_CEILING,
            "segments": ";".join(f"{a}-{b}" for a, b in segs),
            "ecto_spacer_aa": ecto,
            "cyto_tail_aa": cyto,
            "leader": spec["leader"],
            "confidence": spec["confidence"],
            "pmid": spec["pmid"],
            "strategy": spec["strategy"],
            "protein": prot,
        })

    rows.sort(key=lambda r: r["truncated_aa"])
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    print(f"  {'gene':12} {'full':>6} {'trunc':>6} {'shrink':>7} {'LTR-LTR':>7}  9kb   confidence")
    for r in rows:
        print(f"  {r['gene']:12} {r['full_aa']:6} {r['truncated_aa']:6} {r['shrink']:>7} "
              f"{r['ltr_to_ltr_tagged_bp']:7}  {'yes' if r['within_9kb'] else ' NO':4}  {r['confidence']}")
    print(f"\n  {sum(1 for r in rows if r['within_9kb'])}/{len(rows)} within "
          f"{GENOME_CEILING} bp LTR-to-LTR, tag cassette and barcode included")
    if problems:
        print("\n  problems:")
        for p in problems:
            print(f"    {p}")
    else:
        print("  every construct keeps its transmembrane span and starts with a leader")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
