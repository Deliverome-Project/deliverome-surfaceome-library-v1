"""Rank alternate forms of every library gene, over UniProt *and* Ensembl, with the
gene's own MS coverage carried alongside so "not detected" can be read honestly.

Two things this adds over `isoform_expression.py`:

**The universe is UniProt plus Ensembl.** Tedman's isoform arm showed the limit of a
UniProt-only list: 85 of his 177 deliberately non-canonical constructs match no UniProt
proteoform at all. An experimentalist working from Ensembl considers forms UniProt has not
annotated, so both are pooled and deduplicated by protein sequence. Each form records
whether UniProt knows it, Ensembl knows it, or both.

**"Not detected" is qualified by the gene's coverage.** A form with no diagnostic peptide in
a protein that is itself barely seen says nothing. Sequence coverage of the gene's reference
form is computed from the same observed peptides, and a form that is undetected in a
poorly-covered gene is tiered separately from one undetected in a gene that is covered well.
The first is uninformative and should stay in the library; only the second is a real
negative.

    uv run --locked python isoform_universe.py --surfaceome-repo ~/Git/accessible-surfaceome
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import statistics
from collections import Counter, defaultdict
from pathlib import Path

from select_principal import (
    APPRIS_RELEASE,
    UNIPROT_HUMAN,
    UNIPROT_RELEASE,
    load_uniprot_proteoforms,
    parse_appris,
    parse_ensembl_pep,
    parse_uniprot_fasta,
    select_principal,
)
from select_reference import (
    ATLAS_BUILD,
    MIN_PEP,
    MIN_SAMPLES,
    fetch_peptides,
    translations_by_gene,
)

TAKES_APPRIS = {"appris_only", "resolved_appris"}

# A form this close to the reference -- same length, at most this many residues apart -- is
# the reference under a different allele or a cDNA sequence conflict, not an alternate form.
# Matches resolve_candidate_cds.VARIANT_TOLERANCE, which would otherwise build it from the
# reference's own transcript (CD177: Q8N6Q3 is 3 conflicts from GRCh38; 2026-10-29).
ALLELE_TOLERANCE = 3


def same_protein_modulo_allele(seq: str, ref: str | None) -> bool:
    if ref is None or len(seq) != len(ref):
        return False
    return sum(a != b for a, b in zip(seq, ref, strict=True)) <= ALLELE_TOLERANCE


# Below this fraction of the reference sequence covered by observed peptides, the gene is
# too thinly seen for "no diagnostic peptide" to mean anything. Chosen from the data -- see
# the distribution printed at the end of a run.
LOW_COVERAGE_PCT = 10.0


def coverage_pct(seq: str, peptides: list[str]) -> float:
    """Percent of `seq` covered by any observed peptide."""
    covered = bytearray(len(seq))
    for pep in peptides:
        start = seq.find(pep)
        while start != -1:
            for i in range(start, start + len(pep)):
                covered[i] = 1
            start = seq.find(pep, start + 1)
    return 100.0 * sum(covered) / len(seq) if seq else 0.0


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

    library = sorted(set(json.loads(Path(args.library).read_text())))
    uniprot = load_uniprot_proteoforms(cache / UNIPROT_HUMAN)
    cohort_iso = set(parse_uniprot_fasta(seqdir / "human_isoforms_from_afdb_non_hla.fasta"))
    selection = {r["symbol"]: r for r in csv.DictReader(open(args.selection), delimiter="\t")}
    appris = parse_appris(cache / f"appris_homo_sapiens_{APPRIS_RELEASE}.txt")
    pep = parse_ensembl_pep(cache / "Homo_sapiens.GRCh38.pep.all.fa.gz")
    ensembl_by_gene = translations_by_gene(cache / "Homo_sapiens.GRCh38.pep.all.fa.gz")

    # ENST per sequence, so an Ensembl-only form can be named.
    enst_for_seq: dict[str, list[str]] = defaultdict(list)
    for enst, seq in pep.items():
        enst_for_seq[seq].append(enst)

    appris_by_symbol: dict[str, list[dict]] = defaultdict(list)
    for rows in appris.values():
        for r in rows:
            appris_by_symbol[r["symbol"]].append(r)
    ensg_for_acc: dict[str, set[str]] = defaultdict(set)
    with (repo / "data/external/uniprot_ensembl_xrefs/ensg_to_uniprot.tsv").open() as fh:
        next(fh)
        for line in fh:
            g, a = line.split("\t")[:2]
            ensg_for_acc[a].add(g)

    # ---- build the per-gene universe -------------------------------------------
    universe: dict[str, list[dict]] = {}
    for gene in library:
        forms: dict[str, dict] = {}
        for p in uniprot.get(gene, []):
            forms[p["seq"]] = {
                "uniprot": p["accession"],
                "is_uniprot_canonical": not p["is_isoform"],
                "ensts": [],
            }
        accs = [p["accession"] for p in uniprot.get(gene, []) if not p["is_isoform"]]
        ensgs = sorted({g for a in accs for g in ensg_for_acc.get(a, ())}) or sorted(
            {r["ensg"] for r in appris_by_symbol.get(gene, [])}
        )
        for g in ensgs:
            for seq in ensembl_by_gene.get(g, ()):
                rec = forms.setdefault(
                    seq, {"uniprot": "", "is_uniprot_canonical": False, "ensts": []}
                )
                rec["ensts"] = sorted(enst_for_seq.get(seq, []))
        if forms:
            universe[gene] = [dict(v, seq=k) for k, v in forms.items()]
    print(f"library genes with any form: {len(universe)}")
    print(f"distinct forms: {sum(len(v) for v in universe.values())}")

    # ---- which form is the reference -------------------------------------------
    reference: dict[str, str] = {}
    for gene, forms in universe.items():
        rows = [
            r for g in sorted(ensg_for_acc.get(
                next((f["uniprot"] for f in forms if f["is_uniprot_canonical"]), ""), ()))
            for r in appris.get(g, [])
        ] or appris_by_symbol.get(gene, [])
        mane = [pep[r["enst"]] for r in rows if r["mane"] and r["enst"] in pep]
        sel = select_principal(rows, pep) if rows else None
        dec = selection.get(gene, {}).get("decision", "")
        pick = (sel["seq"] if sel else None) if dec in TAKES_APPRIS else (
            mane[0] if mane else (sel["seq"] if sel else None)
        )
        # The library's selected form wins over the MANE/APPRIS pick: reference_overrides.tsv
        # moves six GPI genes onto the anchored form, and reference_selection.tsv carries
        # the result. Ignoring it offered ACHE's selected P22303-2 back as an "alternate"
        # and hid the T subunit (2026-10-29).
        chosen = selection.get(gene, {}).get("selected_uniprot_equivalent", "").strip()
        pick = next((f["seq"] for f in forms if chosen and f["uniprot"] == chosen), pick)
        if pick is None:
            pick = next((f["seq"] for f in forms if f["is_uniprot_canonical"]), None)
        if pick:
            reference[gene] = pick

    accessions = sorted({f["uniprot"] for v in universe.values() for f in v if f["uniprot"]})
    print(f"accessions to query: {len(accessions)}")
    observed = (
        fetch_peptides(accessions, cache, offline=args.offline)
    )

    rows_out = []
    for gene, forms in sorted(universe.items()):
        ref = reference.get(gene)
        obs: dict[str, tuple[int, int]] = {}
        for f in forms:
            if f["uniprot"]:
                obs.update(observed.get(f["uniprot"], {}))
        # Peptides actually belonging to this gene: those found in one of its forms.
        gene_peps = [p for p in obs if any(p in f["seq"] for f in forms)]
        gene_cov = coverage_pct(ref, gene_peps) if ref else 0.0
        gene_best = max((obs[p][1] for p in gene_peps), default=0)

        for f in forms:
            is_ref = ref is not None and (
                f["seq"] == ref or same_protein_modulo_allele(f["seq"], ref)
            )
            others = [o["seq"] for o in forms if o["seq"] != f["seq"]]
            hits = sorted(
                (
                    (p, obs[p][0], obs[p][1])
                    for p in gene_peps
                    if p in f["seq"] and not any(p in o for o in others)
                ),
                key=lambda x: (x[2], x[1], len(x[0])),
                reverse=True,
            )
            detectable = any(
                not any(f["seq"][i : i + MIN_PEP] in o for o in others)
                for i in range(len(f["seq"]) - MIN_PEP + 1)
            )
            best = max((h[2] for h in hits), default=0)
            low_cov = gene_cov < LOW_COVERAGE_PCT
            tier = (
                "expressed" if len(hits) >= 2 and best >= MIN_SAMPLES
                else "likely_expressed" if hits and best >= MIN_SAMPLES
                else "weak" if hits
                else "undetectable" if not detectable
                else "not_detected_low_gene_coverage" if low_cov
                else "not_detected"
            )
            rows_out.append(
                {
                    "gene": gene,
                    "is_reference": is_ref,
                    "uniprot": f["uniprot"],
                    "ensembl_transcripts": ";".join(f["ensts"][:4]),
                    "source": (
                        "both" if f["uniprot"] and f["ensts"]
                        else "uniprot_only" if f["uniprot"]
                        else "ensembl_only"
                    ),
                    "length": len(f["seq"]),
                    "reference_length": len(ref) if ref else "",
                    "in_current_cohort": f["uniprot"] in cohort_iso,
                    "tier": "reference" if is_ref else tier,
                    "diagnostic_peptides": len(hits),
                    "best_n_samples": best,
                    "peptide_examples": ";".join(h[0] for h in hits[:3]),
                    "detectable": detectable,
                    "gene_coverage_pct": round(gene_cov, 1),
                    "gene_peptides": len(gene_peps),
                    "gene_best_n_samples": gene_best,
                }
            )

    rows_out.sort(key=lambda r: (-r["best_n_samples"], -r["diagnostic_peptides"], r["gene"]))
    with open("isoform_universe.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows_out[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows_out)
    print(f"wrote isoform_universe.tsv ({len(rows_out)} alternate forms)")

    covs = sorted({(r["gene"], r["gene_coverage_pct"]) for r in rows_out})
    cv = [c for _g, c in covs]
    summary = {
        "uniprot_release": UNIPROT_RELEASE,
        "peptideatlas_build": ATLAS_BUILD,
        "genes": len(universe),
        "forms_total": len(rows_out),
        "reference_forms": sum(1 for r in rows_out if r["is_reference"]),
        "alternate_forms": sum(1 for r in rows_out if not r["is_reference"]),
        "source": dict(Counter(r["source"] for r in rows_out).most_common()),
        "tier": dict(Counter(r["tier"] for r in rows_out).most_common()),
        "low_coverage_threshold_pct": LOW_COVERAGE_PCT,
        "gene_coverage_pct": {
            "median": round(statistics.median(cv), 1) if cv else None,
            "genes_under_threshold": sum(1 for c in cv if c < LOW_COVERAGE_PCT),
            "genes_total": len(cv),
        },
        "ensembl_only_forms": sum(1 for r in rows_out if r["source"] == "ensembl_only"),
        "ensembl_only_expressed": sum(
            1
            for r in rows_out
            if r["source"] == "ensembl_only" and r["tier"].startswith(("expressed", "likely"))
        ),
    }
    Path("isoform_universe_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
