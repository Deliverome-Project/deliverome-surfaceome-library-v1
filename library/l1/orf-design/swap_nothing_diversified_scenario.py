"""Barcode-swap model for a counterfactual order with no codon diversification at all.

The brief's bottom-line swap figure (Becca's 2026-11-22 comment: "we need to know the bottom line
for the lib with or without codon diversification") compares three scenarios on one order and
one shared tail. ``swap_tags_recoded_scenario.py`` gives "tags recoded too" and the model on the
order of record gives "linkers diversified". This gives the third, "nothing diversified", on the
same order of record (``audit_domestication_barcodes.RECODED_KEY``, T7 promoter included):

1. In every insert, the 51-aa cassette is replaced by p003's single encoding
   (``build_constructs.CASSETTE``), so every arm carries the same cassette DNA.
2. In every snorkel insert, the 55-aa snorkel module (``l1_snorkel.PREFIX_AA``) is replaced by
   the original shared encoding (``l1_snorkel.PREFIX_DNA``), undoing
   ``diversify_snorkel_modules.py``.
3. The protein must be unchanged and the insert must keep its Type IIS sites exactly (the two
   wrapper Esp3I sites, no new BbsI/Esp3I/BsaI at a seam).
4. ``l1_swap_pool_model`` runs unchanged on the counterfactual inserts. Its committed outputs are
   not touched: the run writes only ``swap_nothing_diversified_scenario.json`` (summary) here and
   its per-insert tables to a scratch folder given by ``--detail-dir``.

Nothing here changes the order of record, an encodings table or a builder. Everything else in
each insert (native CDS as recoded for synthesis, the shared 3' tail with T7, the barcode) is the
order of record's, so the only difference from the order-of-record run is the cassette and
snorkel-module DNA.

Implemented by Claude Code session l1-review-export-scope, 2026-11-22.

    uv run --locked python swap_nothing_diversified_scenario.py
"""

from __future__ import annotations

import argparse
import json
import sys
import tempfile
from pathlib import Path

import audit_domestication_barcodes as audit
import l1_snorkel
import l1_swap_pool_model as model
from build_constructs import CASSETTE, CASSETTE_AA, ENZYMES
from swap_tags_recoded_scenario import HEAD, longest_shared

from deliverome_analysis.cloning import seqtools

HERE = Path(__file__).resolve().parent
OUT = HERE / "swap_nothing_diversified_scenario.json"


def span_of(insert: str, peptide: str) -> tuple[int, int] | None:
    """(start, end) of ``peptide``'s DNA in the insert, read in frame from the ATG; None if absent."""
    coding, err = audit.coding_of(insert)
    assert not err, err
    protein = seqtools.translate(coding)
    at = protein.find(peptide)
    if at < 0:
        return None
    assert protein.find(peptide, at + 1) < 0, f"{peptide[:12]}... found twice"
    start = len(HEAD) + 3 * at
    return start, start + 3 * len(peptide)


def unify(insert: str) -> tuple[str, bool]:
    """The insert with p003's cassette DNA and, for a snorkel, the shared module DNA."""
    new = insert
    span = span_of(new, CASSETTE_AA)
    assert span is not None, "no cassette"
    new = new[: span[0]] + CASSETTE + new[span[1] :]
    module = span_of(new, l1_snorkel.PREFIX_AA)
    if module is not None:
        new = new[: module[0]] + l1_snorkel.PREFIX_DNA + new[module[1] :]
    old_c, _ = audit.coding_of(insert)
    new_c, err = audit.coding_of(new)
    assert not err and seqtools.translate(new_c) == seqtools.translate(old_c), "protein changed"

    def sites(s: str) -> list[tuple[str, int]]:
        return sorted((e, st) for e, st, _ in seqtools.find_sites(s, ENZYMES))

    assert sites(new) == sites(insert), "Type IIS sites changed"
    return new, module is not None


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--detail-dir", type=Path, default=None, help="per-insert model tables")
    a = ap.parse_args(argv)

    rec, digest = audit.recoded_inserts(None)
    assert digest == audit.RECODED_SHA256, "order of record changed"
    inserts = {r["key"]: r["insert"] for r in rec}
    new, snorkels = {}, 0
    for key, ins in inserts.items():
        try:
            new[key], is_snorkel = unify(ins)
        except AssertionError as e:
            raise SystemExit(f"{key}: {e}") from e
        snorkels += is_snorkel
    changed = sum(new[k] != inserts[k] for k in inserts)

    # Run the pooled model on the counterfactual inserts, its outputs redirected.
    detail = a.detail_dir or Path(tempfile.mkdtemp(prefix="swap-nothing-diversified-"))
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
        "scenario": (
            "no codon diversification: every cassette on p003's single encoding "
            "(build_constructs.CASSETTE), every snorkel module on l1_snorkel.PREFIX_DNA"
        ),
        "order_of_record": audit.RECODED_KEY,
        "order_of_record_sha256": audit.RECODED_SHA256,
        "inserts": len(new),
        "inserts_changed": changed,
        "snorkel_modules_replaced": snorkels,
        "protein_unchanged": True,
        "type_iis_sites_unchanged": True,
        "longest_shared_cassette_nt": {
            "current_order": longest_shared(
                [inserts[k][slice(*span_of(inserts[k], CASSETTE_AA))] for k in inserts]
            ),
            "nothing_diversified": len(CASSETTE),
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
    print(f"changed {changed} of {len(new)} inserts; {snorkels} snorkel modules replaced")
    print(f"wrote {OUT.name}; per-insert tables in {detail}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
