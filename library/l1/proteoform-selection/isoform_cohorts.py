"""Emit the isoform cohort at two evidence tiers, so a smaller order is a filter not a rebuild.

**Tier 1 — full (1,928 forms, 962 genes).** Everything that clears the candidate bar and
keeps its membrane anchoring. 1,104 of these rest on APPRIS TRIFID alone, which is a
prediction that tracks the observed ordering well (AUC 0.731) but is still a prediction.

**Tier 2 — observed (824 forms, 564 genes).** The subset where something was actually seen:
a diagnostic peptide in PeptideAtlas, or an OpenProt translation score above zero.

The line between the tiers is deliberately *observed vs predicted* rather than a threshold.
The tighter alternatives available in the data -- TE>=3 giving 280 forms, TE>=5 giving 199 --
are arbitrary cuts on a count of studies, and a count of studies is a measure of how often a
protein has been looked at as much as of whether the isoform is real. The observed/predicted
boundary is the one that means something.

Both tiers are written so the full set stays usable; tier 2 is a column, not a separate
pipeline.

    uv run --locked python isoform_cohorts.py
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from collections import Counter
from pathlib import Path
from library import library_members
from gpi import anchor_retained
from gpi import gpi_genes as gpi_gene_set

OBSERVED_TIERS = {"expressed", "likely_expressed", "weak"}
# s1e3's topology_status, applied as an absolute test rather than only relative to the
# reference. The relative test ("does this form lose the TM helices its reference has?")
# cannot fire when the reference has no topology call, and 12 forms reached the cohort that
# way -- real soluble isoforms such as ADAM12-S, soluble PVR and shed OSMR. They are
# genuine proteins and not surface-displayable, which is what the library is for.
SURFACE_LABELS = {"TM", "SP+TM", "BETA"}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--universe", default="isoform_universe_scored.tsv")
    ap.add_argument("--verdicts", default="topology_verdicts.tsv")
    ap.add_argument("--manifest", default="topology_shortlist.tsv")
    ap.add_argument("--three-line", default="topology_shortlist.3line")
    ap.add_argument("--surfaceome-repo", default="~/Git/accessible-surfaceome")
    ap.add_argument("--cache-dir", default="~/.cache/deliverome/appris")
    ap.add_argument("--netgpi-calls", default="netgpi_gpi_form_calls.tsv")
    args = ap.parse_args()

    verdicts = {
        (r["gene"], r["length"], r["uniprot"]): r
        for r in csv.DictReader(open(args.verdicts), delimiter="\t")
    }
    # The universe is a snapshot and keeps genes the library has since dropped.
    members = library_members()
    rows = [
        r for r in csv.DictReader(open(args.universe), delimiter="\t")
        if r["is_reference"] == "False" and (not members or r["gene"] in members)
    ]

    # DeepTMHMM label per form, from the cohort tables or the new run.
    repo = Path(args.surfaceome_repo).expanduser()
    label: dict[str, str] = {}
    for name in ("deeptmhmm_human_canonical.tsv", "deeptmhmm_human_isoforms.tsv"):
        path = repo / "data/processed/deeptmhmm" / name
        if path.exists():
            with path.open() as fh:
                for row in csv.DictReader(fh, delimiter="\t"):
                    acc = row.get("uniprot_accession_full") or row.get("uniprot_accession")
                    if acc:
                        label[acc] = row["deeptmhmm_label"]
    fid_for = {
        (m["gene"], m["length"]): m["fasta_id"]
        for m in csv.DictReader(open(args.manifest), delimiter="\t")
        if m["is_reference"] == "False"
    }
    label_for_fid: dict[str, str] = {}
    three = Path(args.three_line)
    if three.exists():
        lines = three.read_text().splitlines()
        for i in range(len(lines) - 2):
            if lines[i].startswith(">"):
                nm, lb = lines[i][1:].rsplit("|", 1)
                label_for_fid[nm.strip()] = lb.strip()
    # Reference and isoform sequences, so a GPI gene's anchor is inherited by an isoform
    # only when that isoform keeps the C-terminus carrying the signal.
    cache = Path(args.cache_dir).expanduser()
    ref_seq: dict[str, str] = {}
    hu = cache / "uniprot_human_reviewed.tsv.gz"
    if hu.exists():
        with gzip.open(hu, "rt") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                sym = (row.get("Gene Names (primary)") or "").strip()
                if sym and row.get("Sequence"):
                    ref_seq.setdefault(sym, row["Sequence"])
    iso_seq: dict[str, str] = {}
    hi = cache / "uniprot_human_reviewed_isoforms.tsv.gz"
    if hi.exists():
        with gzip.open(hi, "rt") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                if row.get("Sequence"):
                    iso_seq[row["Entry"]] = row["Sequence"]
    seq_for_fid: dict[str, str] = {}
    if three.exists():
        lines = three.read_text().splitlines()
        for i in range(len(lines) - 2):
            if lines[i].startswith(">"):
                nm, _lb = lines[i][1:].rsplit("|", 1)
                seq_for_fid[nm.strip()] = lines[i + 1].strip()

    def form_seq(r) -> str | None:
        if r["uniprot"] and r["uniprot"] in iso_seq:
            return iso_seq[r["uniprot"]]
        fid = fid_for.get((r["gene"], r["length"]))
        return seq_for_fid.get(fid) if fid else None

    gpi_genes = gpi_gene_set()

    # NetGPI measures the anchor on the form itself. The C-terminal retention test is only a
    # proxy for it, and a proxy that inverts when the reference is not the anchored form:
    # it asks whether the isoform keeps the REFERENCE's C-terminus, so for ACHE, BCAN and
    # TFPI -- whose references carry no anchor -- it rejected precisely the anchored isoform
    # the reference was being moved onto. Where NetGPI has a call it decides.
    netgpi_call: dict[tuple[str, str], str] = {}
    calls = Path(args.netgpi_calls)
    if calls.exists():
        with calls.open() as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                netgpi_call[(row["gene"], row["form"])] = row["call"]

    def inherits_anchor(r) -> bool:
        """Does this isoform carry the gene's GPI anchor?"""
        if r["gene"] not in gpi_genes:
            return False
        for tag in [r["uniprot"]] + (r["ensembl_transcripts"] or "").split(";"):
            call = netgpi_call.get((r["gene"], tag.strip()))
            if call in ("anchored", "not_anchored"):
                return call == "anchored"
        return anchor_retained(form_seq(r), ref_seq.get(r["gene"])) is not False

    def topology_label(r) -> str | None:
        lab = label.get(r["uniprot"]) if r["uniprot"] else None
        if lab is None:
            fid = fid_for.get((r["gene"], r["length"]))
            lab = label_for_fid.get(fid) if fid else None
        return lab

    def is_surface(r) -> bool | None:
        lab = topology_label(r)
        if lab is None:
            return None
        return lab in SURFACE_LABELS or inherits_anchor(r)

    def te(r) -> int:
        v = r["openprot_te_score"]
        return int(v) if v not in ("", None) else 0

    def trifid(r) -> float:
        v = r["trifid_score"]
        return float(v) if v not in ("", None) else 0.0

    def excluded(r) -> bool:
        v = verdicts.get((r["gene"], r["length"], r["uniprot"]))
        return r["loses_all_tm_vs_reference"] == "True" or (
            v is not None and v["verdict"] == "exclude_lost_surface"
        )

    out = []
    restricted: list[tuple[str, str, str | None, bool]] = []
    for r in rows:
        observed = r["tier"] in OBSERVED_TIERS or te(r) > 0
        passes = observed or trifid(r) > 0.5
        if not passes or excluded(r):
            continue
        if is_surface(r) is False:
            restricted.append((r["gene"], r["uniprot"] or r["ensembl_transcripts"],
                               topology_label(r), observed,
                               "anchor_not_retained" if r["gene"] in gpi_genes
                               else "no_anchor"))
            continue
        out.append(
            {
                "gene": r["gene"],
                "uniprot": r["uniprot"],
                "ensembl_transcripts": r["ensembl_transcripts"],
                "length": r["length"],
                "reference_length": r["reference_length"],
                "source": r["source"],
                "in_current_cohort": r["in_current_cohort"],
                "cohort_tier": "observed" if observed else "full_only",
                "evidence": (
                    "protein (deep peptides)" if r["tier"] in ("expressed", "likely_expressed")
                    else "protein (weak)" if r["tier"] == "weak"
                    else "translation (OpenProt TE)" if te(r) > 0
                    else "prediction (TRIFID)"
                ),
                "peptide_tier": r["tier"],
                "openprot_te_score": r["openprot_te_score"],
                "trifid_score": r["trifid_score"],
                "appris_label_class": r["appris_label_class"],
            }
        )

    with open("isoform_cohorts.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]), delimiter="\t")
        w.writeheader()
        w.writerows(out)
    print(f"wrote isoform_cohorts.tsv ({len(out)} forms)")
    # The soluble restriction: evidence-passing forms dropped for being explicitly SP or
    # GLOB with no GPI anchor on the parent gene.
    print(
        f"soluble restriction removed {len(restricted)} evidence-passing forms "
        f"({sum(1 for x in restricted if x[3])} of them observed-tier) over "
        f"{len({x[0] for x in restricted})} genes: "
        + dict(Counter(x[2] for x in restricted).most_common()).__repr__()
    )
    with open("isoform_soluble_restricted.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["gene", "form", "deeptmhmm_label", "observed_tier", "reason"])
        w.writerows(sorted(restricted))

    obs = [r for r in out if r["cohort_tier"] == "observed"]
    summary = {
        "tier_1_full": {
            "forms": len(out),
            "genes": len({r["gene"] for r in out}),
            "already_in_cohort": sum(1 for r in out if r["in_current_cohort"] == "True"),
            "evidence": dict(Counter(r["evidence"] for r in out).most_common()),
        },
        "tier_2_observed": {
            "forms": len(obs),
            "genes": len({r["gene"] for r in obs}),
            "already_in_cohort": sum(1 for r in obs if r["in_current_cohort"] == "True"),
            "evidence": dict(Counter(r["evidence"] for r in obs).most_common()),
            "reduction_vs_full": f"{100 * (1 - len(obs) / len(out)):.0f}%",
        },
    }
    Path("isoform_cohorts_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
