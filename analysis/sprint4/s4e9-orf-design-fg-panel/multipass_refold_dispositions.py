"""Per-member dispositions for the multipass residue-1 arms that failed the tagged refold.

Reads l1-tag-refold's verdicts for the rung (`VERDICTS`, ColabFold model_1_ptm, 3 recycles,
native MSA reused; thresholds in `score_tagged_folds.py`) and writes
`multipass_refold_dispositions.tsv`, one row per arm that did not pass, plus the clustering
summary printed below. Pass/fail thresholds are the refold session's; this file only decides
what each failure means for the build.

Failure classes, assigned in order:

* ``epitope_packed`` -- at least ``CONFIDENT_EPITOPE`` SpyTag003/ALFA residues confidently
  packed on the native domain. The only class where the model is confident the tag is
  occluded.
* ``fold_disturbed`` -- the native region moved against the native model (native-region
  lDDT < 0.85).
* ``tail_local`` -- only the region within 15 residues of the cassette moved (lDDT < 0.7), with
  the native region intact (lDDT >= 0.85). That region is the member's own N-tail, which is
  disordered in the native model too, so its lDDT compares two disordered conformations.
* ``low_confidence_burial`` -- SpyTag003 buried (RSA below the refold threshold) in a stretch
  the model is not confident about (SpyTag003 pLDDT < 50), with fewer than
  ``CONFIDENT_EPITOPE`` confidently packed epitope residues: AF2 compacting a disordered tag.

**Family-level rule (2026-11-14): every failing arm is kept as built and flagged.** The options
are keep-as-built-flagged or hold, since the alternatives tested do not help:

* The cassette size is fixed ("keep full", Becca 2026-11-11).
* A 10-residue spacer between the cassette and native residue 2 rescued none of the 12
  fold-disturbed arms (`verdicts_multipass_spacer.tsv`). A 30-residue spacer made the
  snorkel arms worse.
* No internal site can be vetted for these members (no membrane-side mask or per-residue
  features off the internal rung).

No arm has a specific reason to hold, because each verdict is a single-model prediction and
the pooled screen measures display directly. Each row records the pattern, the metric behind
it, the rejected alternative, and the bench result that would change it. The
"buried SpyTag003 at the mature N-terminus" pattern (53 arms) is an open design question for
Becca, written up in `NTERM_RES1_MULTIPASS.md`.

Attribution: decided_by l1-nterm-res1-multipass (delegated by Becca 2026-11-12, "use your best
judgment"); implemented by the same session, PR #280. Verdicts by l1-tag-refold.

    uv run --locked python multipass_refold_dispositions.py
"""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd
from scipy.stats import fisher_exact

from deliverome_analysis.s3 import processed

HERE = Path(__file__).parent
TAGGED = "sprint4/s4e9-orf-design-fg-panel/predicted-structures/tagged/2026-11-12/"
VERDICTS = TAGGED + "verdicts_multipass.tsv"
#: The expression-rescue arms' verdicts. Since one arm per gene (Becca, 2026-11-13), 23 library
#: references on this rung are built Igk + tag, the same mature protein as their retired rescue
#: twin, so the twin's verdict is the reference's, matched on the mature protein's sha256.
RESCUE_VERDICTS = TAGGED + "verdicts_rescue.tsv"
PLAN = HERE / "nterm_res1_multipass.tsv"
CONSTRUCTS = HERE / "constructs_l1_multipass.tsv"
OUT = HERE / "multipass_refold_dispositions.tsv"
#: Confidently packed epitope residues that count as a real occlusion call. Five is the
#: smallest count at which the packed residues span a turn of SpyTag003 or ALFA.
CONFIDENT_EPITOPE = 5
#: SpyTag003 pLDDT below which a burial call is treated as low confidence (AF2's own
#: "low confidence" band starts at 50).
LOW_PLDDT = 50
DELEGATED = 'l1-nterm-res1-multipass (delegated by Becca 2026-11-12, "use your best judgment")'
IMPLEMENTED_BY = (
    "Claude Code session l1-nterm-res1-multipass (PR #280; per-member record 2026-11-14)"
)

DECISION_DATE = "2026-11-14"

#: failure class -> the pattern name used in the write-up.
PATTERN = {
    "epitope_packed": "epitope packed",
    "fold_disturbed": "fold disturbed",
    "tail_local": "tail only",
    "low_confidence_burial": "buried SpyTag003",
}

#: failure class -> the bench result that would change the disposition.
BENCH = {
    "epitope_packed": (
        "on non-permeabilized cells, ALFA stains but SpyCatcher003 does not: SpyTag003 is "
        "occluded on the surface, so move it (C-terminal or loop) for L2"
    ),
    "fold_disturbed": (
        "neither ALFA nor SpyCatcher003 stains on non-permeabilized cells while the permeabilized "
        "or total signal shows expression: the arm misfolds or is retained, so drop or redesign "
        "for L2"
    ),
    "tail_local": (
        "neither ALFA nor SpyCatcher003 stains on non-permeabilized cells while the permeabilized "
        "or total signal shows expression"
    ),
    "low_confidence_burial": (
        "on non-permeabilized cells, ALFA stains but SpyCatcher003 does not: the burial is real, "
        "so move SpyTag003 for L2"
    ),
}

DISPOSITION = {
    "epitope_packed": (
        "keep as built, flagged; display risk: the model packs SpyTag003 onto the receptor's "
        "extracellular surface. ALFA stays exposed, so read ALFA as the display control",
        "move to an internal loop site (none assessed; these are class A receptors with short "
        "loops) or shorten the cassette (ruled out: keep full)",
    ),
    "fold_disturbed": (
        "keep as built, flagged; structural risk: the cassette perturbs the region next to TM1 "
        "in a single-model fold, and a 10-residue spacer rescued none of the 12 tested",
        "add a 10- or 30-residue spacer (0 of 12 rescued; 30 made snorkels worse), move to an "
        "internal site (none can be vetted), or hold the member",
    ),
    "tail_local": (
        "keep as built, flagged: the native fold is intact (native-region lDDT >= 0.85); only the "
        "member's own disordered N-tail next to the cassette moved",
        "treat as a failure and redesign",
    ),
    "low_confidence_burial": (
        "keep as built, flagged: SpyTag003 is called buried with no confident epitope packing "
        "(pLDDT < 50 for most), typical of AF2 compacting a disordered tag. ALFA stays exposed, "
        "so read ALFA as the display control",
        "hold the arm, or move SpyTag003 (no site assessed; the cassette is fixed)",
    ),
}


def klass(r: pd.Series) -> str:
    if r["epitope_confident_contact_res"] >= CONFIDENT_EPITOPE:
        return "epitope_packed"
    reasons = str(r["verdict_reasons"])
    if "native-region lDDT" in reasons:
        return "fold_disturbed"
    if "within 15 of the site" in reasons:
        return "tail_local"
    if "SpyTag003 buried" in str(r["verdict_reasons"]) and r["spytag003_plddt"] < LOW_PLDDT:
        return "low_confidence_burial"
    # Burial with SpyTag003 pLDDT >= 50 but < CONFIDENT_EPITOPE packed residues.
    return "low_confidence_burial"


def metric(r) -> str:
    """The numbers the class rests on, for one failing arm."""
    if r.failure_class == "fold_disturbed":
        return (
            f"native-region lDDT {r.lddt_confident:.3f} (pass >= 0.85); TM-score {r.tm_native:.3f}"
        )
    if r.failure_class == "tail_local":
        return (
            f"site lDDT below 0.7 ({r.verdict_reasons}); native-region lDDT {r.lddt_confident:.3f}"
        )
    return (
        f"SpyTag003 RSA {r.spytag003_rsa:.3f} (buried < 0.25), pLDDT {r.spytag003_plddt:.1f}; "
        f"{int(r.epitope_confident_contact_res)} epitope residues confidently packed; "
        f"ALFA RSA {r.alfa_rsa:.3f}"
    )


def family(r: pd.Series) -> str:
    if r["gpcr"] != 1:
        return "non-GPCR"
    if r["tedman_class"] == "A (Olfactory)" or re.match(r"OR\d", r["gene_symbol"]):
        return "olfactory GPCR"
    if r["gene_symbol"].startswith("TAS2R"):
        return "taste GPCR"
    return "other GPCR"


def main(argv: list[str] | None = None) -> int:
    import argparse

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--check", action="store_true", help="compare with the record; write nothing")
    args = ap.parse_args(argv)
    v = pd.read_csv(processed.download(VERDICTS), sep="\t")
    etag = processed.list_etags(VERDICTS).get(VERDICTS, "")
    plan = pd.read_csv(PLAN, sep="\t").set_index("design_id")
    built = pd.read_csv(CONSTRUCTS, sep="\t")
    # The file also carries l1-tag-refold's folds of the order of record (no `design`); their
    # verdicts match this rung's 561 of 561, so only this build's rows are read.
    v = v[v["design"].notna()].copy()
    v["construct_id"] = v["design_id"] + ":" + v["design"]
    built_mature = built.set_index("construct_id")["mature_sha256"]
    # A construct refolded after a rebuild has a row per sequence: keep the current one.
    v = v[
        ~v["construct_id"].isin(built_mature.index)
        | (v["tagged_sha256"] == v["construct_id"].map(built_mature))
    ].drop_duplicates("construct_id", keep="last")
    built_ids = set(built.loc[built["insert_nt"].notna(), "construct_id"])
    extra = set(v["construct_id"]) - built_ids
    # The rescue references' tag-alone verdicts are for an arm replaced by Igk + tag.
    rescue_refs = set(built.loc[built["reason"].str.startswith("one arm per gene"), "design_id"])
    dropped = {c for c in extra if c.split(":")[0] in rescue_refs and c.endswith(":tag_alone")}
    assert not extra - dropped, f"verdicts for arms no longer built: {sorted(extra - dropped)}"
    v = v[~v["construct_id"].isin(dropped)]
    if dropped:
        print(f"{len(dropped)} tag-alone verdicts dropped: those references are now Igk + tag")
    mature = built.set_index("construct_id")["mature_sha256"]
    rv = pd.read_csv(processed.download(RESCUE_VERDICTS), sep="\t")
    rv_by_sha = rv.drop_duplicates("tagged_sha256").set_index("tagged_sha256")
    take = [c for c in built_ids - set(v["construct_id"]) if mature[c] in rv_by_sha.index]
    if take:
        add = rv_by_sha.loc[[mature[c] for c in take]].reset_index()
        add["construct_id"] = take
        add["design_id"] = [c.split(":")[0] for c in take]
        add["design"] = [c.split(":")[1] for c in take]
        add["gene_symbol"] = add["design_id"].map(built.set_index("design_id")["gene_symbol"])
        v = pd.concat([v, add], ignore_index=True)
        print(f"{len(take)} verdicts taken from the retired rescue twins (same mature protein)")
    # Arms built after the verdict file was written have no verdict yet: listed, not dropped.
    missing = sorted(built_ids - set(v["construct_id"]))
    if missing:
        add = built[built["construct_id"].isin(missing)][
            ["construct_id", "design_id", "gene_symbol", "design"]
        ]
        add = add.assign(
            verdict="not_folded",
            tagged_sha256=add["construct_id"].map(built.set_index("construct_id")["mature_sha256"]),
        )
        v = pd.concat([v, add], ignore_index=True)
    assert (v["tagged_sha256"] == v["construct_id"].map(mature)).all(), "verdicts are stale"
    v["from_rescue"] = v["construct_id"].isin(set(take))
    for col in ("gpcr", "tedman_class", "n_tail_aa", "candidate_set", "species"):
        v[col] = v["design_id"].map(plan[col])
    v["tedman_class"] = v["tedman_class"].fillna("")
    v["family"] = v.apply(family, axis=1)

    scored = v[v["verdict"].isin(["pass", "fail"])]
    print(f"verdicts {VERDICTS} (etag {etag}): {v['verdict'].value_counts().to_dict()}")
    for name, col in (("design", "design"), ("family", "family")):
        t = pd.crosstab(scored[col], scored["verdict"])
        t["fail_rate"] = (t["fail"] / t.sum(axis=1)).round(3)
        print(f"\nby {name}:\n{t.to_string()}")
    t = pd.crosstab(scored["design"], scored["verdict"])
    print(f"\ndesign, Fisher exact: p = {fisher_exact(t.values).pvalue:.2g}")
    for fam in ("other GPCR", "taste GPCR", "olfactory GPCR"):
        t = pd.crosstab(scored["family"] == fam, scored["verdict"])
        print(f"{fam} vs the rest, Fisher exact: p = {fisher_exact(t.values).pvalue:.2g}")
    tail = pd.cut(scored["n_tail_aa"], [0, 20, 30, 45, 60, 10_000])
    t = pd.crosstab(tail, scored["verdict"])
    t["fail_rate"] = (t["fail"] / t.sum(axis=1)).round(3)
    print(f"\nby native N-tail (aa):\n{t.to_string()}")

    out = v[v["verdict"] != "pass"].copy()
    out["failure_class"] = out.apply(
        lambda r: klass(r) if r["verdict"] == "fail" else "not_folded", axis=1
    )
    rows = []
    for r in out.itertuples():
        if r.failure_class == "not_folded":
            disp, rejected = (
                "keep as built; unscored (not folded in the overnight run)",
                "hold the arm until folded",
            )
        else:
            disp, rejected = DISPOSITION[r.failure_class]
        rows.append(
            {
                "construct_id": r.construct_id,
                "design_id": r.design_id,
                "gene_symbol": r.gene_symbol,
                "species": r.species,
                "candidate_set": r.candidate_set,
                "design": r.design,
                "family": r.family,
                "n_tail_aa": r.n_tail_aa,
                "verdict": r.verdict,
                "failure_class": r.failure_class,
                "verdict_reasons": r.verdict_reasons if isinstance(r.verdict_reasons, str) else "",
                "spytag003_plddt": r.spytag003_plddt,
                "spytag003_rsa": r.spytag003_rsa,
                "alfa_rsa": r.alfa_rsa,
                "epitope_confident_contact_res": r.epitope_confident_contact_res,
                "lddt_confident": r.lddt_confident,
                "tm_native": r.tm_native,
                "pattern": PATTERN.get(r.failure_class, "not folded"),
                "metric": metric(r) if r.failure_class in PATTERN else "",
                "disposition": disp,
                "build_changed": 0,
                "rejected_alternative": rejected,
                "bench_result_that_would_change_it": BENCH.get(r.failure_class, "a fold"),
                "decided_by": DELEGATED,
                "implemented_by": IMPLEMENTED_BY,
                "decision_date": DECISION_DATE,
                "verdict_source": (
                    f"processed://{RESCUE_VERDICTS if r.from_rescue else VERDICTS}; l1-tag-refold"
                    + ("" if r.from_rescue else f" (etag {etag})")
                ),
                "tagged_sha256": r.tagged_sha256,
            }
        )
    res = pd.DataFrame(rows).sort_values(["failure_class", "design_id"])
    # An unfolded arm has no count; keep the column integer rather than let NaN make it float.
    res["epitope_confident_contact_res"] = res["epitope_confident_contact_res"].astype("Int64")
    text = res.to_csv(sep="\t", index=False, lineterminator="\n")
    if args.check:
        same = OUT.read_text() == text
        print(
            f"{OUT.name}: " + ("record matches rebuild" if same else "record DIFFERS from rebuild")
        )
        return 0 if same else 1
    OUT.write_text(text)
    print(f"\nwrote {OUT.name}: {len(res)} rows; {res['failure_class'].value_counts().to_dict()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
