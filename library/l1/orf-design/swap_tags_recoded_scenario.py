"""Barcode-swap model for a counterfactual order whose tags are recoded too, not only the linkers.

The order of record (``audit_domestication_barcodes.RECODED_KEY``) gives every construct its own
encoding of the four GGGGS linkers, but SpyTag003 and P-ALFA-P keep p003's codons in all of
them (Becca, 2026-10-29), so every pair of constructs still shares the two tags (48 and 45 nt).
The brief's bottom-line swap table needs the third row, "tags recoded too", measured rather
than estimated (Becca's 2026-11-22 comment: "we need to know the bottom line for the lib with or
without codon diversification").

This builds that counterfactual and runs the same pooled model on it:

1. For every insert of the order of record, the 51-aa cassette is redrawn whole, one synonymous
   encoding per construct, with ``seqtools.diversify_library_linkers`` and no fixed codons. The
   existing linker encodings cannot be pinned while the tags change (a pin fixes the whole
   cassette), so every cassette is redrawn. Constraints as in ``diversify_arm_cassettes.py``:
   no BbsI/Esp3I/BsaI in the cassette or at either seam, the 15-nt junction rule
   (``LINKER_ONLY_JUNCTION_IDENTITY``), variable 20-mers issued elsewhere avoided where possible,
   and runs of up to 5 identical bases for any group of more than 400 constructs sharing the
   sequence before their cassette (the Igκ group), 4 otherwise. Seeded, so it reproduces.
2. Each new cassette is spliced into its insert in place of the old one. The protein must be
   unchanged and the insert must carry only its two wrapper Esp3I sites and no BbsI/BsaI.
3. ``l1_swap_pool_model`` runs unchanged on the counterfactual inserts. Its committed outputs are
   not touched: the run writes only ``swap_tags_recoded_scenario.json`` (summary) here and its
   per-insert tables to a scratch folder given by ``--detail-dir``.

Nothing here changes the order of record, an encodings table or a builder.

    uv run --locked python swap_tags_recoded_scenario.py
    uv run --locked python swap_tags_recoded_scenario.py --check   # redraw; encodings must match
"""

from __future__ import annotations

import argparse
import hashlib
import json
import random
import sys
import tempfile
from collections import Counter
from pathlib import Path

import audit_domestication_barcodes as audit
import l1_swap_pool_model as model
from build_constructs import CASSETTE_AA, ENZYMES
from diversify_cassettes import LINKER_ONLY_JUNCTION_IDENTITY

from deliverome_analysis.cloning import seqtools

HERE = Path(__file__).resolve().parent
OUT = HERE / "swap_tags_recoded_scenario.json"
SEED = 2026_11_22
FLANK = 150
CROWDED_GROUP = 400
CROWDED_MAX_HOMOPOLYMER = 5
HEAD = audit.WRAP_5 + audit.KOZAK


def cassette_span(insert: str) -> tuple[int, int]:
    """(start, end) of the cassette DNA in the insert, read in frame from the ATG."""
    coding, err = audit.coding_of(insert)
    assert not err, err
    at = seqtools.translate(coding).find(CASSETTE_AA)
    assert at >= 0, "no cassette"
    start = len(HEAD) + 3 * at
    return start, start + 3 * len(CASSETTE_AA)


def redraw(inserts: dict[str, str]) -> dict[str, str]:
    """One whole-cassette encoding per insert, tags included."""
    members = {}
    for key, ins in inserts.items():
        s, e = cassette_span(ins)
        members[key] = (ins[max(0, s - FLANK) : s], ins[e : e + FLANK])
    groups = Counter(f5[-60:] for f5, _ in members.values())
    # sorted, not a set: set order follows the string hash seed, and draw order sets the encodings
    crowded = sorted(k for k, (f5, _) in members.items() if groups[f5[-60:]] > CROWDED_GROUP)
    common = dict(
        kmer_policy="minimize",
        enzymes=ENZYMES,
        max_junction_identity=LINKER_ONLY_JUNCTION_IDENTITY,
    )
    first = seqtools.diversify_library_linkers(
        {k: members[k] for k in crowded},
        CASSETTE_AA,
        rng=random.Random(SEED),
        max_homopolymer=CROWDED_MAX_HOMOPOLYMER,
        **common,
    )
    pinned = {k: first.variants[k].sequence for k in crowded}
    plan = seqtools.diversify_library_linkers(
        members, CASSETTE_AA, existing=pinned, rng=random.Random(SEED), **common
    )
    return {k: plan.variants[k].sequence for k in members}


def splice(inserts: dict[str, str], cassettes: dict[str, str]) -> dict[str, str]:
    out = {}
    for key, ins in inserts.items():
        s, e = cassette_span(ins)
        new = ins[:s] + cassettes[key] + ins[e:]
        old_c, _ = audit.coding_of(ins)
        new_c, err = audit.coding_of(new)
        assert not err and seqtools.translate(new_c) == seqtools.translate(old_c), key
        sites = seqtools.find_sites(new, ENZYMES)
        assert len(sites) == 2 and all(x[0] == "Esp3I" for x in sites), (key, sites)
        out[key] = new
    return out


def longest_shared(seqs: list[str]) -> int:
    """Longest stretch any two of ``seqs`` share (binary search on an L-mer index)."""

    def collides(n: int) -> bool:
        seen: dict[str, int] = {}
        for i, s in enumerate(seqs):
            for j in range(len(s) - n + 1):
                k = s[j : j + n]
                if seen.setdefault(k, i) != i:
                    return True
        return False

    lo, hi = 0, max(len(s) for s in seqs)
    while lo < hi:
        mid = (lo + hi + 1) // 2
        lo, hi = (mid, hi) if collides(mid) else (lo, mid - 1)
    return lo


def encodings_sha(cassettes: dict[str, str]) -> str:
    blob = "".join(f"{k}\t{v}\n" for k, v in sorted(cassettes.items()))
    return hashlib.sha256(blob.encode()).hexdigest()


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--check", action="store_true", help="redraw; fail if the encodings differ")
    ap.add_argument("--detail-dir", type=Path, default=None, help="per-insert model tables")
    a = ap.parse_args(argv)

    rec, digest = audit.recoded_inserts(None)
    assert digest == audit.RECODED_SHA256, "order of record changed"
    inserts = {r["key"]: r["insert"] for r in rec}
    cassettes = redraw(inserts)
    sha = encodings_sha(cassettes)
    if a.check:
        recorded = json.loads(OUT.read_text())["encodings_sha256"]
        print("encodings reproduce" if sha == recorded else "encodings DIFFER")
        return 0 if sha == recorded else 1
    assert len(set(cassettes.values())) == len(cassettes), "two inserts share a cassette"
    new = splice(inserts, cassettes)
    current = [inserts[k][slice(*cassette_span(inserts[k]))] for k in inserts]

    # Run the pooled model on the counterfactual inserts, its outputs redirected.
    detail = a.detail_dir or Path(tempfile.mkdtemp(prefix="swap-tags-recoded-"))
    detail.mkdir(parents=True, exist_ok=True)
    real_recoded = audit.recoded_inserts
    audit.recoded_inserts = lambda _p: (
        [{**r, "insert": new[r["key"]]} for r in rec],
        audit.RECODED_SHA256,
    )
    model.OUT_TSV = detail / "l1_swap_pool_model.tsv"
    model.OUT_PAIRS = detail / "l1_swap_pool_pairs.tsv"
    model.OUT_JSON = detail / "l1_swap_pool_model_summary.json"
    model.plot = lambda *args, **kwargs: None
    try:
        model.main()
    finally:
        audit.recoded_inserts = real_recoded
    summary = json.loads(model.OUT_JSON.read_text())

    out = {
        "scenario": "SpyTag003 and P-ALFA-P recoded per construct as well as the linkers",
        "order_of_record": audit.RECODED_KEY,
        "order_of_record_sha256": audit.RECODED_SHA256,
        "seed": SEED,
        "inserts": len(new),
        "encodings_sha256": sha,
        "pinning": "none: every cassette redrawn whole (a pin fixes the whole cassette)",
        "longest_shared_cassette_nt": {
            "current_order": longest_shared(current),
            "tags_recoded": longest_shared(list(cassettes.values())),
        },
        "library_fraction_wrong_protein": summary["library_fraction_wrong_protein"],
        "wrong_protein_by_block_region": summary["wrong_protein_by_block_region"],
        "per_barcode_impurity": {
            k: v for k, v in summary["per_barcode_impurity"].items() if k != "by_rung_median"
        },
        "model": "l1_swap_pool_model.main() unchanged, on the spliced inserts",
    }
    OUT.write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out["library_fraction_wrong_protein"], indent=2))
    print("longest shared cassette nt:", out["longest_shared_cassette_nt"])
    print(f"wrote {OUT.name}; per-insert tables in {detail}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
