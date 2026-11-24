"""Rank the library's alternate isoforms by evidence that they are actually expressed.

Isoform membership today rests on UniProt annotating the isoform, which says a splice form
has been *described*, not that it is *made*. This scores each alternate isoform on observed
protein evidence instead: peptides in PeptideAtlas that fall in residues unique to that
isoform among all proteoforms of its gene, weighted by how many distinct samples saw them.

The scoring mirrors the reference-selection logic. A peptide counts for an isoform only if:

* it lies in a stretch present in that isoform and **absent from every other proteoform of
  the same gene** -- otherwise it says the gene is expressed, not which form;
* it was observed in PeptideAtlas build Human 2026-01;
* ``n_samples`` is carried through, because a peptide seen once and a peptide seen in 1,852
  samples are not the same evidence.

The denominator matters as much as the hits. An isoform differing by a few residues has
nowhere for a diagnostic peptide to sit, so ``detectable`` records whether the isoform had
at least one 7-residue unique stretch at all. Without it, "no evidence" and "cannot be
evidenced" look the same, and they are not.

    uv run --locked python isoform_expression.py --surfaceome-repo ~/Git/accessible-surfaceome
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

from select_reference import (
    ATLAS_BUILD,
    MAX_PEP,
    MIN_PEP,
    MIN_SAMPLES,
    fetch_peptides,
)
from select_principal import UNIPROT_HUMAN, UNIPROT_RELEASE, load_uniprot_proteoforms, parse_uniprot_fasta

def diagnostic_windows(seq: str, others: list[str]) -> bool:
    """Does any peptide-length window of `seq` occur in no other proteoform of the gene?

    This is the question that decides whether an isoform is *evidenceable* at all, and it is
    not the same as having unique residues. Most alternate isoforms are cassette-exon skips
    and so are strict subsequences of the canonical -- every residue they contain is present
    in the longer form. What is novel is the **junction** the skip creates: the sequence
    either side exists, the join does not. A window spanning that join is diagnostic even
    though none of its residues is.
    """
    for i in range(len(seq) - MIN_PEP + 1):
        win = seq[i : i + MIN_PEP]
        if not any(win in o for o in others):
            return True
    return False


def diagnostic_peptides(
    seq: str, others: list[str], observed: dict[str, tuple[int, int]]
) -> list[tuple[str, int, int]]:
    """Observed peptides present in this isoform and in no other proteoform of the gene.

    Driven from the observed peptides rather than from the sequence, which is both far
    cheaper and the right semantics: only a peptide somebody actually saw can be evidence.
    """
    hits = []
    for pep, (n_obs, n_samples) in observed.items():
        if pep in seq and not any(pep in o for o in others):
            hits.append((pep, n_obs, n_samples))
    return sorted(hits, key=lambda x: (x[2], x[1], len(x[0])), reverse=True)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--surfaceome-repo", default="~/Git/accessible-surfaceome")
    ap.add_argument("--cache-dir", default="~/.cache/deliverome/appris")
    ap.add_argument("--library", default="../tag-site-architecture/library_members.json")
    ap.add_argument("--selection", default="reference_selection.tsv")
    ap.add_argument("--offline", action="store_true")
    args = ap.parse_args()

    repo = Path(args.surfaceome_repo).expanduser()
    cache = Path(args.cache_dir).expanduser()
    seqdir = repo / "data/external/deeptmhmm_surfaceome_inputs"

    library = set(json.loads(Path(args.library).read_text()))
    uniprot = load_uniprot_proteoforms(cache / UNIPROT_HUMAN)
    cohort_iso = parse_uniprot_fasta(seqdir / "human_isoforms_from_afdb_non_hla.fasta")
    in_cohort = set(cohort_iso)
    selection = {
        r["symbol"]: r for r in csv.DictReader(open(args.selection), delimiter="\t")
    }

    genes = sorted(
        g for g in library if any(p["is_isoform"] for p in uniprot.get(g, []))
    )
    print(f"library genes with at least one UniProt isoform: {len(genes)}")

    # Query every proteoform, not just the canonical. PeptideAtlas indexes isoform
    # biosequences separately (P16070-3 returns a different peptide set from P16070), so
    # asking only for the canonical hides exactly the isoform-specific peptides this is
    # looking for -- which is what an earlier version of this script did, and it found 3
    # expressed isoforms in 2,802.
    accessions = sorted({p["accession"] for g in genes for p in uniprot[g]})
    observed = (
        {a: {} for a in accessions}
        if args.offline
        else fetch_peptides(accessions, cache)
    )

    rows = []
    for gene in genes:
        forms = uniprot[gene]
        canonical = next((p for p in forms if not p["is_isoform"]), None)
        if not canonical:
            continue
        # Pool every proteoform's peptides for the gene; the containment test below is what
        # assigns a peptide to one isoform, not which query returned it.
        obs: dict[str, tuple[int, int]] = {}
        for form in forms:
            obs.update(observed.get(form["accession"], {}))
        sel = selection.get(gene, {})
        for iso in forms:
            if not iso["is_isoform"]:
                continue
            others = [p["seq"] for p in forms if p["accession"] != iso["accession"]]
            ranked = diagnostic_peptides(iso["seq"], others, obs)
            best_samples = max((h[2] for h in ranked), default=0)
            detectable = diagnostic_windows(iso["seq"], others)
            rows.append(
                {
                    "gene": gene,
                    "isoform": iso["accession"],
                    "isoform_length": len(iso["seq"]),
                    "canonical_length": len(canonical["seq"]),
                    "in_current_cohort": iso["accession"] in in_cohort,
                    "is_selected_reference": sel.get("selected_uniprot_equivalent")
                    == iso["accession"],
                    "detectable": detectable,
                    "diagnostic_peptides": len(ranked),
                    "best_n_samples": best_samples,
                    "best_n_observations": max((h[1] for h in ranked), default=0),
                    "peptide_examples": ";".join(h[0] for h in ranked[:3]),
                    "n_peptides_for_gene": len(obs),
                    "tier": (
                        "expressed"
                        if len(ranked) >= 2 and best_samples >= MIN_SAMPLES
                        else "likely_expressed"
                        if ranked and best_samples >= MIN_SAMPLES
                        else "weak"
                        if ranked
                        else "not_detected"
                        if detectable
                        else "undetectable"
                    ),
                }
            )

    columns = list(rows[0])
    rows.sort(key=lambda r: (-r["best_n_samples"], -r["diagnostic_peptides"], r["gene"]))
    with open("isoform_priority.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=columns, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote isoform_priority.tsv ({len(rows)} isoforms)")

    tiers = Counter(r["tier"] for r in rows)
    cohort_rows = [r for r in rows if r["in_current_cohort"]]
    summary = {
        "uniprot_release": UNIPROT_RELEASE,
        "peptideatlas_build": ATLAS_BUILD,
        "genes_with_isoforms": len(genes),
        "isoforms_scored": len(rows),
        "tier": dict(tiers.most_common()),
        "undetectable_share": round(100 * tiers["undetectable"] / len(rows), 1),
        "current_cohort": {
            "scored": len(cohort_rows),
            "tier": dict(Counter(r["tier"] for r in cohort_rows).most_common()),
        },
        "expressed_but_not_in_cohort": sum(
            1
            for r in rows
            if r["tier"] in ("expressed", "likely_expressed") and not r["in_current_cohort"]
        ),
        "in_cohort_but_not_detected": sum(
            1 for r in cohort_rows if r["tier"] == "not_detected"
        ),
    }
    Path("isoform_priority_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
