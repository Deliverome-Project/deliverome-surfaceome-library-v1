"""Barcode-swap model for a counterfactual order with a T7 promoter just 5' of the amplicon handle.

Becca, 2026-11-22: "maybe we just add the t7 promoter site just 5' of jt832 and recalculate
swapping as long as it doesn't hurt too much." The promoter (``T7``, forward, ending in the GGG
transcription start) is inserted between the 2-nt pad after the stop codon and the 5' amplicon
handle (JT832, ``build_constructs.FWD_PRIMER``). In fixed cells T7 then transcribes sense RNA
across the handle, the barcode and the readout region, so the Feldman RT primer and padlock read
it unchanged (PerturbView, Kudo et al. 2025; CROPseq-multi, Walton et al.).

The cost is that the promoter sits between the ORF and the barcode, so every construct shares
20 more nt there, and a template switch inside it joins one ORF to another's barcode. This
inserts it into every insert of the order of record (protein unchanged; checked for new
BbsI/Esp3I/BsaI sites and for GT[AG]AG splice-donor motifs across both seams) and runs
``l1_swap_pool_model`` unchanged on the result, its outputs redirected. Nothing in the order of
record changes.

Implemented by Claude Code session l1-review-export-scope, 2026-11-22.

    uv run --locked python swap_t7_scenario.py
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import tempfile
from pathlib import Path

import audit_domestication_barcodes as audit
import l1_swap_pool_model as model
from build_constructs import ENZYMES, FWD_PRIMER

from deliverome_analysis import seqtools

HERE = Path(__file__).parent
OUT = HERE / "swap_t7_scenario.json"
#: T7 class III promoter, -17 to +3 (the +1 GGG start).
T7 = "TAATACGACTCACTATAGGG"
SPLICE_DONOR = re.compile(r"GT[AG]AG")


def toggle_t7(insert: str) -> tuple[str, int]:
    """If T7 is already present before FWD_PRIMER (as in the adopted order of record), remove it
    to model the pre-T7 baseline; otherwise insert T7 to model adding it."""
    at = insert.find(FWD_PRIMER)
    assert at > 0 and insert.find(FWD_PRIMER, at + 1) < 0, "handle not found exactly once"
    if at >= len(T7) and insert[at - len(T7) : at] == T7:
        return insert[: at - len(T7)] + insert[at:], -len(T7)
    return insert[:at] + T7 + insert[at:], len(T7)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--detail-dir", type=Path, default=None)
    a = ap.parse_args(argv)
    rec, digest = audit.recoded_inserts(None)
    assert digest == audit.RECODED_SHA256, "order of record changed"
    toggled = {r["key"]: toggle_t7(r["insert"]) for r in rec}
    new = {k: seq for k, (seq, _) in toggled.items()}
    deltas = {d for _, d in toggled.values()}
    assert len(deltas) == 1, f"mixed T7 presence across order: {deltas}"
    delta = next(iter(deltas))
    for r in rec:
        old, ins = r["insert"], new[r["key"]]
        sites = lambda s: sorted((e, st) for e, st, _ in seqtools.find_sites(s, ENZYMES))  # noqa: E731
        assert sites(ins) == sites(old), f"{r['key']}: Type IIS sites changed"
    t7_seq = rec[0]["insert"] if delta < 0 else new[rec[0]["key"]]
    at = t7_seq.find(FWD_PRIMER)
    seam = t7_seq[at - len(T7) - 12 : at + 12]
    donors = [m.start() for m in SPLICE_DONOR.finditer(seam)]
    seams = {
        (r["insert"] if delta < 0 else new[r["key"]])[
            (r["insert"] if delta < 0 else new[r["key"]]).find(FWD_PRIMER) - len(T7) - 12 : (
                r["insert"] if delta < 0 else new[r["key"]]
            ).find(FWD_PRIMER)
            + 12
        ]
        for r in rec
    }

    detail = a.detail_dir or Path(tempfile.mkdtemp(prefix="swap-t7-"))
    detail.mkdir(parents=True, exist_ok=True)
    real = audit.recoded_inserts
    audit.recoded_inserts = lambda _p: ([{**r, "insert": new[r["key"]]} for r in rec], digest)
    model.OUT_TSV = detail / "l1_swap_pool_model.tsv"
    model.OUT_PAIRS = detail / "l1_swap_pool_pairs.tsv"
    model.OUT_JSON = detail / "l1_swap_pool_model_summary.json"
    model.plot = lambda *args, **kwargs: None
    try:
        model.main()
    finally:
        audit.recoded_inserts = real
    summary = json.loads(model.OUT_JSON.read_text())
    scenario_desc = (
        "counterfactual without T7 promoter (removed 20 nt 5' of JT832; order of record carries T7)"
        if delta < 0
        else "forward T7 promoter (TAATACGACTCACTATAGGG) just 5' of the JT832 handle"
    )
    out = {
        "scenario": scenario_desc,
        "order_of_record": audit.RECODED_KEY,
        "order_of_record_sha256": digest,
        "inserts": len(new),
        "added_nt_per_insert": delta,
        "distinct_seams": len(seams),
        "splice_donor_motifs_in_seam": donors,
        "library_fraction_wrong_protein": summary["library_fraction_wrong_protein"],
        "wrong_protein_by_block_region": summary["wrong_protein_by_block_region"],
        "model": "l1_swap_pool_model.main() unchanged, on the counterfactual inserts",
    }
    OUT.write_text(json.dumps(out, indent=2) + "\n")
    print(
        json.dumps(
            {
                k: out[k]
                for k in (
                    "splice_donor_motifs_in_seam",
                    "distinct_seams",
                    "library_fraction_wrong_protein",
                )
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
