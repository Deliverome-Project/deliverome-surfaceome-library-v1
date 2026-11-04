"""Apply s1e3's topology rule to the newly predicted forms and decide which to exclude.

The rule is taken from `analysis/sprint1/s1e3-surface-library-filters/build_library.py` on
purpose, so the 1,527 new forms are judged exactly as the existing 2,238 were: the DeepTMHMM
label decides, with TM / SP+TM / BETA meaning surface, SP meaning secreted and GLOB meaning
intracellular. Inventing a second rule here would leave the library with two inconsistent
filters, which is the failure that started this thread when APPRIS THUMP was used instead of
a topology predictor.

A form is excluded when it is **not surface** while its reference **is** -- it has lost the
membrane anchoring the construct depends on. Losing only some helices is kept, per the
design call: that changes the protein but can still leave something displayable.

    uv run --locked python apply_topology.py
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path

SURFACE_LABELS = {"TM", "SP+TM", "BETA"}


def parse_3line(path: Path) -> dict[str, dict]:
    """{fasta id: {label, topology string, tm_helices}} from DeepTMHMM's 3-line output."""
    out: dict[str, dict] = {}
    lines = path.read_text().splitlines()
    i = 0
    while i < len(lines):
        if not lines[i].startswith(">"):
            i += 1
            continue
        header = lines[i][1:]
        name, _, label = header.partition("|")
        # The id itself contains "|", so the label is whatever follows the LAST separator.
        name = header.rsplit("|", 1)[0].strip()
        label = header.rsplit("|", 1)[1].strip()
        topo = lines[i + 2].strip() if i + 2 < len(lines) else ""
        # A helix is a maximal run of M; count transitions into M rather than M characters.
        helices = sum(
            1 for j, c in enumerate(topo) if c == "M" and (j == 0 or topo[j - 1] != "M")
        )
        out[name] = {"label": label, "topology": topo, "tm_helices": helices}
        i += 3
    return out


def status(label: str) -> str:
    """s1e3's topology_status, label-first."""
    lab = (label or "").upper()
    if lab in SURFACE_LABELS:
        return "surface"
    if lab == "SP":
        return "secreted"
    return "intracellular"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--three-line", default="topology_shortlist.3line")
    ap.add_argument("--manifest", default="topology_shortlist.tsv")
    ap.add_argument("--surfaceome-repo", default="~/Git/accessible-surfaceome")
    args = ap.parse_args()

    pred = parse_3line(Path(args.three_line))
    manifest = list(csv.DictReader(open(args.manifest), delimiter="\t"))
    print(f"predictions: {len(pred)}; manifest rows: {len(manifest)}")

    # Existing DeepTMHMM calls, for references that already had one.
    repo = Path(args.surfaceome_repo).expanduser()
    existing: dict[str, dict] = {}
    for name in ("deeptmhmm_human_canonical.tsv", "deeptmhmm_human_isoforms.tsv"):
        path = repo / "data/processed/deeptmhmm" / name
        if path.exists():
            with path.open() as fh:
                for row in csv.DictReader(fh, delimiter="\t"):
                    acc = row.get("uniprot_accession_full") or row.get("uniprot_accession")
                    if acc:
                        existing[acc] = {
                            "label": row["deeptmhmm_label"],
                            "tm_helices": int(row["tm_helix_count"] or 0),
                        }

    # Reference status per gene: from this run where it was queued, else from the cohort.
    ref_status: dict[str, dict] = {}
    for m in manifest:
        if m["is_reference"] == "True" and m["fasta_id"] in pred:
            p = pred[m["fasta_id"]]
            ref_status[m["gene"]] = {"label": p["label"], "tm_helices": p["tm_helices"]}
    universe = list(csv.DictReader(open("isoform_universe_scored.tsv"), delimiter="\t"))
    for r in universe:
        if r["is_reference"] == "True" and r["gene"] not in ref_status:
            e = existing.get(r["uniprot"])
            if e:
                ref_status[r["gene"]] = e

    rows = []
    for m in manifest:
        if m["is_reference"] == "True":
            continue
        p = pred.get(m["fasta_id"])
        if not p:
            continue
        ref = ref_status.get(m["gene"])
        st = status(p["label"])
        ref_st = status(ref["label"]) if ref else ""
        verdict = (
            "exclude_lost_surface"
            if ref_st == "surface" and st != "surface"
            else "keep_lost_some_tm"
            if ref and st == "surface" and 0 < p["tm_helices"] < ref["tm_helices"]
            else "keep"
            if ref
            else "undecidable_no_reference_call"
        )
        rows.append(
            {
                **{k: m[k] for k in ("gene", "uniprot", "ensembl_transcripts", "length",
                                     "source", "tier", "trifid_score", "reason")},
                "deeptmhmm_label": p["label"],
                "topology_status": st,
                "tm_helices": p["tm_helices"],
                "reference_label": ref["label"] if ref else "",
                "reference_tm_helices": ref["tm_helices"] if ref else "",
                "verdict": verdict,
            }
        )

    with open("topology_verdicts.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote topology_verdicts.tsv ({len(rows)} alternate forms)")

    summary = {
        "sequences_predicted": len(pred),
        "alternate_forms_judged": len(rows),
        "references_predicted": sum(1 for m in manifest if m["is_reference"] == "True"),
        "verdict": dict(Counter(r["verdict"] for r in rows).most_common()),
        "topology_status": dict(Counter(r["topology_status"] for r in rows).most_common()),
        "deeptmhmm_label": dict(Counter(r["deeptmhmm_label"] for r in rows).most_common()),
        "excluded_by_source": dict(
            Counter(
                r["source"] for r in rows if r["verdict"] == "exclude_lost_surface"
            ).most_common()
        ),
        "excluded_with_evidence": sum(
            1
            for r in rows
            if r["verdict"] == "exclude_lost_surface"
            and r["tier"] in ("expressed", "likely_expressed", "weak")
        ),
    }
    Path("topology_verdicts_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
