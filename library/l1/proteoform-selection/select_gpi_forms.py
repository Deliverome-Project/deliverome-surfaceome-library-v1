"""Choose the best GPI-anchored form in each slot: reference, isoform, mouse, cyno.

The architecture attribute says a gene is anchored; NetGPI says which of its forms carries
the anchor. This fills the four slots the library needs with forms that actually do.

Anchoring is a **veto, not a sort key**: a form is excluded only when NetGPI confidently
says it is not anchored. Among everything else the project's own source precedence decides —
MANE > APPRIS > targeted experimental evidence > UniProt canonical, the same order
``select_reference.py`` enforces — so a marginally-better-anchored form does not displace a
better-sourced one. Within a precedence tier a confirmed anchored form is preferred over an
unresolved one.

**NetGPI's likelihood is not a ranking signal.** For an anchored call it is the confidence in
the omega-site *position*, not in the anchoring, so ranking by it would sort by how
crisply the attachment point is located rather than by how good the form is. It is used
here only as the binary filter it is. Ranking comes from the evidence this project already
has: observed peptides first, then APPRIS principal status, then translation evidence, then
TRIFID, preferring a reviewed UniProt form over an Ensembl-only model.

The reference slot keeps the current selection whenever it is anchored, so a gene only moves
when it has to. Ortholog slots are already filled by the anchor guard in
``repick_orthologs.py``, which picks the best-corresponding candidate that keeps the anchor.

    uv run --locked python select_gpi_forms.py
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path

from library import library_members
from gpi import gpi_genes

# Mirrors select_reference.PRECEDENCE.
PRECEDENCE = {"mane": 0, "appris": 1, "evidence": 2, "uniprot_canonical": 3}

TIER_RANK = {
    "expressed": 0, "likely_expressed": 1, "weak": 2,
    "not_detected_low_gene_coverage": 3, "not_detected": 4, "undetectable": 5,
}
APPRIS_RANK = {"PRINCIPAL": 0, "ALTERNATIVE": 1}


def num(v: str, default: float = 0.0) -> float:
    try:
        return float(v)
    except (TypeError, ValueError):
        return default


def source_rung(r: dict, mane_enst: str, appris_enst: str) -> str:
    """Which precedence rung this form sits on."""
    transcripts = {t.strip() for t in (r.get("ensembl_transcripts") or "").split(";")}
    if mane_enst and mane_enst in transcripts:
        return "mane"
    if (appris_enst and appris_enst in transcripts) or (
        r.get("appris_label_class") or ""
    ).startswith("PRINCIPAL"):
        return "appris"
    # Observed peptides only. `select_reference.py` sets this rung at three supporting
    # peptides across ten samples, so the `weak` tier is below it by definition and a lone
    # Ribo-seq event is not evidence of this kind at all. Both still rank, lower down, as
    # tie-breaks. Counting them here promoted an undetectable Ensembl model over NCAM-120,
    # and a weak-tier fragment over CD55's own canonical.
    if TIER_RANK.get(r.get("tier", ""), 9) <= 1:
        return "evidence"
    return "uniprot_canonical"


def rank(r: dict, mane_enst: str, appris_enst: str, anchored: bool) -> tuple:
    """Lower is better. Source precedence leads; anchoring only breaks ties within a rung."""
    acc = r.get("uniprot") or ""
    return (
        PRECEDENCE[source_rung(r, mane_enst, appris_enst)],
        0 if anchored else 1,
        # Within the UniProt rung the canonical accession is that rung's own representative:
        # an unsuffixed entry outranks its isoforms and any Ensembl-only model.
        0 if (acc and "-" not in acc) else 1,
        TIER_RANK.get(r.get("tier", ""), 9),
        APPRIS_RANK.get((r.get("appris_label_class") or "").split(":")[0], 2),
        -num(r.get("openprot_te_score", "")),
        -num(r.get("trifid_score", "")),
        0 if r.get("uniprot") else 1,              # reviewed form over Ensembl-only model
        0 if r.get("in_current_cohort") == "True" else 1,
        -int(r.get("length") or 0),
    )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cache-dir", default="~/.cache/deliverome/appris")
    ap.add_argument("--calls", default="netgpi_gpi_form_calls.tsv")
    ap.add_argument("--universe", default="isoform_universe_scored.tsv")
    ap.add_argument("--reference", default="reference_selection.tsv")
    ap.add_argument("--overrides", default="gpi_form_overrides.tsv")
    ap.add_argument("--out", default="gpi_selected_forms.tsv")
    args = ap.parse_args()

    gpi = gpi_genes()
    anchored: dict[str, set[str]] = defaultdict(set)
    call_of: dict[tuple[str, str], str] = {}
    ref_form: dict[str, str] = {}
    for r in csv.DictReader(open(args.calls), delimiter="\t"):
        call_of[(r["gene"], r["form"])] = r["call"]
        # Human slots draw only on human forms. The call table also holds the mouse and
        # cyno picks, and letting those into the pool put mouse Cd14 (P10810) forward as a
        # replacement for the human CD14 reference.
        if r["call"] == "anchored" and r["kind"] in ("ref", "alt"):
            anchored[r["gene"]].add(r["form"])
        if r["kind"] == "ref":
            ref_form[r["gene"]] = r["form"]

    meta: dict[tuple[str, str], dict] = {}
    for r in csv.DictReader(open(args.universe), delimiter="\t"):
        if r["gene"] not in gpi:
            continue
        tag = r["uniprot"] or (r["ensembl_transcripts"] or "").split(";")[0]
        meta[(r["gene"], tag)] = r

    override: dict[str, str] = {}
    op = Path(args.overrides)
    if op.exists():
        for r in csv.DictReader(op.open(), delimiter="\t"):
            override[r["gene"]] = r["form"]

    enst: dict[str, tuple[str, str]] = {}
    for r in csv.DictReader(open(args.reference), delimiter="\t"):
        enst[r["symbol"]] = (r["mane_enst"], r["appris_enst"])

    ortho: dict[str, dict[str, dict]] = defaultdict(dict)
    for species in ("mouse", "cyno"):
        for r in csv.DictReader(open(f"ortholog_repick_{species}.tsv"), delimiter="\t"):
            if r["human_symbol"] in gpi:
                ortho[r["human_symbol"]][species] = r

    rows = []
    members = library_members()
    for gene in sorted(gpi & set(ref_form) & (members or set(ref_form))):
        anc = anchored.get(gene, set())
        cur = ref_form[gene]
        mane_e, appris_e = enst.get(gene, ("", ""))
        # Eligible: every human form of this gene that is not a confident negative.
        eligible = [
            f for (g, f), c in call_of.items()
            if g == gene and f != cur and c != "not_anchored"
            and (gene, f) in meta
        ]
        alts = sorted(
            eligible,
            key=lambda f: rank(meta[(gene, f)], mane_e, appris_e, f in anc),
        )
        cur_call = call_of.get((gene, cur), "")
        if gene in override and override[gene] in anc:
            # A named form, for the cases where the ranking is arbitrary rather than wrong.
            ref_pick, ref_status = override[gene], "moved_override"
        elif cur in anc:
            ref_pick, ref_status = cur, "kept"
        elif cur_call != "not_anchored":
            # Only a confident negative justifies moving. An `unresolved` call is the one
            # the calibration said not to trust — CD14's reference is called non-anchored at
            # 0.388 and CD14 is a textbook GPI-anchored protein.
            ref_pick, ref_status = cur, f"kept_reference_{cur_call or 'no_call'}"
        elif alts:
            ref_pick, ref_status = alts[0], "moved"
        else:
            ref_pick, ref_status = cur, "no_anchored_form"
        # Isoform slot: the next best form that is not the reference pick.
        iso = [f for f in alts if f != ref_pick]
        iso_pick = iso[0] if iso else ""

        row = {"gene": gene,
               "reference_form": ref_pick, "reference_status": ref_status,
               "previous_reference": cur if ref_status == "moved" else "",
               "n_anchored_human_forms": len(anc),
               "isoform_form": iso_pick,
               "reference_rung": source_rung(meta[(gene, ref_pick)], mane_e, appris_e)
               if (gene, ref_pick) in meta else "current",
               "reference_anchor_call": call_of.get((gene, ref_pick), ""),
               "isoform_evidence": (meta.get((gene, iso_pick), {}) or {}).get("tier", ""),
               "isoform_anchor_call": call_of.get((gene, iso_pick), ""),
               "n_other_anchored_isoforms": max(len(iso) - 1, 0)}
        for species in ("mouse", "cyno"):
            o = ortho[gene].get(species)
            row[f"{species}_form"] = o["picked_acc"] if o else ""
            row[f"{species}_anchored"] = (o["picked_anchor_call"] or "no_call") if o else "no_ortholog"
            row[f"{species}_correspondence"] = o["picked_score"] if o else ""
        rows.append(row)

    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {args.out} ({len(rows)} genes)\n")
    print("reference slot:", dict(Counter(r["reference_status"] for r in rows)))
    print("isoform slot filled with an anchored form:",
          sum(1 for r in rows if r["isoform_form"]))
    for species in ("mouse", "cyno"):
        print(f"{species} slot:",
              dict(Counter(r[f"{species}_anchored"] for r in rows)))

    print("\n=== reference moved onto an anchored form ===")
    for r in rows:
        if r["reference_status"] == "moved":
            print(f"  {r['gene']:9} {r['previous_reference']:16} -> {r['reference_form']:18}"
                  f" rung={r['reference_rung']:18} {r['reference_anchor_call']}")
    print("\n=== no anchored form anywhere — review, do not drop ===")
    for r in rows:
        if r["reference_status"].startswith("no_anchored_form"):
            print(f"  {r['gene']:9} reference {r['reference_form']}")


if __name__ == "__main__":
    main()
