"""One arm per gene for the 25 L1 expression-rescue genes.

Becca, 2026-11-13: "we are pooling the whole L1 library together. why resuce pairs no i think
we should make a single choice for each L1". The library is one pool, so a gene keeps one member
and one insert; the L1 goes from 2,735 active members to 2,710.

**Rule.** Each gene keeps its library reference, built with its lane's default design, and its
expression-rescue member is retired. For the 23 genes on the "N at residue 1, multipass" rung the
lane default for a GPCR is Igk + tag (Becca, 2026-11-10), which is exactly the rescue arm's
protein, so the leader the rescue set was meant to test is kept and only the duplicate goes. For
the two genes on the "N after SP" lane (ADGRD2, GLP2R) the reference's design handles the native
signal peptide; the rescue arm does not improve on it (ADGRD2) or would put the tag in front of a
signal peptide predicted to be cleaved (GLP2R).

Refold verdicts never separate the two arms: for the 23 rung genes the kept and retired proteins
are identical once the leader is cleaved, so the retired twin's verdict is the kept arm's. The
call rests on the leader and lane evidence in each row.

Attribution: the one-arm rule is Becca's (2026-11-13); the per-gene choice is
l1-nterm-res1-multipass's under her 2026-11-12 delegation ("use your best judgment"); implemented
by the Claude Code session l1-nterm-res1-multipass. Retirement in the manifest is applied by the
writer (l1-reviewer-package).

    uv run --locked python rescue_one_arm.py
"""

from __future__ import annotations

import csv
from pathlib import Path

import pandas as pd

from deliverome_analysis import library_design as ld
from deliverome_analysis.s3 import processed

HERE = Path(__file__).parent
MANIFEST = HERE / "candidate_manifest_tranches.tsv"
TAG_PATHS = HERE / "l1_tag_paths.tsv"
MULTIPASS = HERE / "constructs_l1_multipass.tsv"
OUT = HERE / "rescue_one_arm.tsv"
TAGGED = "sprint4/s4e9-orf-design-fg-panel/predicted-structures/tagged/2026-11-12/"
RESCUE_VERDICTS = TAGGED + "verdicts_rescue.tsv"
NTERM_VERDICTS = TAGGED + "verdicts_nterm.tsv"
DATE = "2026-11-13"
DECIDED_BY = (
    'Becca, 2026-11-13 (one arm per gene: "make a single choice for each L1"); arm chosen by '
    'l1-nterm-res1-multipass (delegated by Becca 2026-11-12, "use your best judgment")'
)
IMPLEMENTED_BY = "Claude Code session l1-nterm-res1-multipass"

RUNG_REASON = (
    "GPCR on the multipass residue-1 rung: the lane default is Igk + tag (Becca 2026-11-10, on "
    "Aashish Manglik's suggestion), which is the rescue arm's protein, so the reference is built "
    "Igk + tag and the rescue member is retired as its duplicate. The leader is kept for a "
    "receptor chosen for poor display"
)
RUNG_REJECTED = (
    "keep the rescue member and retire the reference (same protein; the reference is the gene's "
    "library member every other table joins on); or keep the reference tag-alone (drops the "
    "leader from a poor expressor, against the GPCR default)"
)
RUNG_CHANGE = (
    "evidence that the Igk leader lowers display of a poorly displayed GPCR; none in Tedman or "
    "the literature reviewed"
)

#: Genes whose reference is on another lane: (reason, rejected, what would change the call).
OFF_RUNG = {
    "ADGRD2": (
        "reference is on the N-after-SP lane with the tag after its native 18-residue signal "
        "peptide (DeepTMHMM v1); both arms pass refold. The rescue arm swaps the native signal "
        "peptide for Igk, and nothing shows the native one underperforms, so the lane's design "
        "for a native signal peptide stands",
        "retire the reference and keep the Igk signal-peptide swap (rescue member l1-p2361)",
        "evidence that ADGRD2's own signal peptide targets poorly compared with Igk",
    ),
    "GLP2R": (
        "reference is on the N-after-SP lane, built with Igk replacing the signal peptide that "
        "DeepTMHMM v2 calls at 1-66 (mouse Glp2r has a classical one); passes refold. The rescue "
        "arm keeps residues 1-66 behind Igk + tag, so a cleaved signal peptide would take the tag "
        "with it, the likely reason Tedman's HA-after-Met1 construct read about 0.12x the median",
        "retire the reference and keep the rescue arm (Igk + tag at Met1, in front of the "
        "predicted signal peptide)",
        "evidence that GLP2R's N-terminal segment is not cleaved (mature N-terminus by Edman or "
        "mass spectrometry)",
    ),
}


def read_table(path: Path) -> list[dict]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def verdict_text(r: pd.Series | None) -> str:
    if r is None:
        return "no verdict"
    why = r["verdict_reasons"] if isinstance(r["verdict_reasons"], str) else ""
    return r["verdict"] + (f" ({why})" if why else "")


def main() -> int:
    manifest = read_table(MANIFEST)
    # All 25 are retired (one arm per gene); the record keeps them.
    rescue = [r for r in manifest if r["set"] == "expression_rescue"]
    refs = {
        r["key"].split(":")[1]: r
        for r in manifest
        if r["set"] == "library_reference" and ld.is_active(r) and r["library"] == "l1"
    }
    lanes = {r["candidate_id"]: r["placement"] for r in read_table(TAG_PATHS)}
    built = pd.read_csv(MULTIPASS, sep="\t").set_index("design_id")
    rv = pd.read_csv(processed.download(RESCUE_VERDICTS), sep="\t").set_index("design_id")
    nv = pd.read_csv(processed.download(NTERM_VERDICTS), sep="\t").drop_duplicates("design_id")
    nv = nv.set_index("design_id")
    assert len(rescue) == 25, f"expected 25 rescue members, found {len(rescue)}"

    rows = []
    for m in sorted(rescue, key=lambda r: r["key"]):
        gene = m["key"].split(":")[1]
        ref = refs[gene]
        d_ref, d_res = ref["design_id"], m["design_id"]
        lane = lanes[ref["candidate_id"]]
        retired_v = rv.loc[d_res] if d_res in rv.index else None
        if gene in OFF_RUNG:
            reason, rejected, change = OFF_RUNG[gene]
            kept_design = nv.loc[d_ref, "arm"] if d_ref in nv.index else ""
            kept_v = nv.loc[d_ref] if d_ref in nv.index else None
            kept_src, kind = f"processed://{NTERM_VERDICTS}", "exception: reference off the rung"
        else:
            assert lane == "N at residue 1, multipass", f"{gene}: reference lane {lane!r}"
            reason, rejected, change = RUNG_REASON, RUNG_REJECTED, RUNG_CHANGE
            kept_design = built.loc[d_ref, "design"]
            # Same mature protein as the retired twin, so the twin's verdict is the kept arm's.
            kept_v, kept_src, kind = retired_v, f"processed://{RESCUE_VERDICTS}", "rule"
            if kept_design == "igk_tag":
                assert built.loc[d_ref, "mature_sha256"] == retired_v["tagged_sha256"], gene
        rows.append(
            {
                "gene_symbol": gene,
                "kept_design_id": d_ref,
                "kept_candidate_id": ref["candidate_id"],
                "kept_lane": lane,
                "kept_design": kept_design,
                "kept_refold": verdict_text(kept_v),
                "kept_refold_source": kept_src,
                "retired_design_id": d_res,
                "retired_candidate_id": m["candidate_id"],
                "retired_design": retired_v["design"] if retired_v is not None else "",
                "retired_refold": verdict_text(retired_v),
                "rule_or_exception": kind,
                "reason": reason,
                "rejected_alternative": rejected,
                "would_change_the_call": change,
                "decided_by": DECIDED_BY,
                "implemented_by": IMPLEMENTED_BY,
                "decision_date": DATE,
            }
        )
    res = pd.DataFrame(rows)
    res.to_csv(OUT, sep="\t", index=False, lineterminator="\n")
    print(
        f"wrote {OUT.name}: {len(res)} genes; {res['rule_or_exception'].value_counts().to_dict()}"
    )
    print(res[["gene_symbol", "kept_design_id", "kept_design", "kept_refold"]].to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
