"""Emit exactly the sequences that need a DeepTMHMM call, and nothing else.

Scoped against what s1e3 already decided. `library/l1/surface-filters`
classifies every isoform in the UniProt cohort by `topology_status` and already drops the
ones that are secreted or intracellular, so those are settled and must not be re-run. What
s1e3 never saw are the Ensembl-sourced forms, which is the gap its own
`isoform_expression_scoping.md` identified.

A form earns a topology call only if it would change a decision: **observed evidence of
expression** (a diagnostic peptide, or OpenProt translation) **or a high TRIFID rank**.

Two corrections over the first version of this shortlist:

* it treated "no TM comparison available" as "the form needs a call", but that value is also
  empty when it is the **reference** that lacks a call. 61 forms were queued that already
  had one;
* and it omitted those references. A TM count is only meaningful against another count, so
  a gene whose reference has never been run cannot be judged however many isoforms are run.
  **233 references are queued alongside the alternates.**

752 candidates need no run at all: both they and their reference already have calls, so they
can be decided from data in hand.

    uv run --locked python shortlist_for_topology.py
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

from select_principal import (
    APPRIS_RELEASE,
    UNIPROT_HUMAN,
    load_uniprot_proteoforms,
    parse_appris,
    parse_ensembl_pep,
)
from select_reference import translations_by_gene

TRIFID_CUT = 0.5


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--surfaceome-repo", default="~/Git/accessible-surfaceome")
    ap.add_argument("--cache-dir", default="~/.cache/deliverome/appris")
    ap.add_argument("--universe", default="isoform_universe_scored.tsv")
    args = ap.parse_args()

    repo = Path(args.surfaceome_repo).expanduser()
    cache = Path(args.cache_dir).expanduser()
    rows = [
        r
        for r in csv.DictReader(open(args.universe), delimiter="\t")
        if r["is_reference"] == "False"
    ]

    def te(r) -> int:
        v = r["openprot_te_score"]
        return int(v) if v not in ("", None) else 0

    def trifid(r) -> float | None:
        v = r["trifid_score"]
        return float(v) if v not in ("", None) else None

    def evidenced(r) -> bool:
        return r["tier"] in ("expressed", "likely_expressed", "weak") or te(r) > 0

    # Accessions DeepTMHMM has already been run on.
    called: set[str] = set()
    for name in ("deeptmhmm_human_canonical.tsv", "deeptmhmm_human_isoforms.tsv"):
        path = repo / "data/processed/deeptmhmm" / name
        if path.exists():
            with path.open() as fh:
                for row in csv.DictReader(fh, delimiter="\t"):
                    acc = row.get("uniprot_accession_full") or row.get("uniprot_accession")
                    if acc:
                        called.add(acc)

    all_rows = list(csv.DictReader(open(args.universe), delimiter="\t"))
    refs = {r["gene"]: r for r in all_rows if r["is_reference"] == "True"}

    candidates = [r for r in rows if evidenced(r) or (trifid(r) or 0) > TRIFID_CUT]
    short = [r for r in candidates if r["uniprot"] not in called]
    decidable = [
        r
        for r in candidates
        if r["uniprot"] in called and (refs.get(r["gene"], {}).get("uniprot") or "_") in called
    ]
    # A reference with no call makes its whole gene unjudgeable, so it joins the run.
    ref_needed = [
        refs[g]
        for g in sorted({r["gene"] for r in candidates})
        if g in refs and (refs[g].get("uniprot") or "_") not in called
    ]
    print(f"candidates meeting the bar: {len(candidates)}")
    print(f"  already decidable with no run: {len(decidable)}")
    print(f"  alternates needing a call:     {len(short)}")
    print(f"  references needing a call:     {len(ref_needed)}")
    short = short + ref_needed
    print(f"  -> run size: {len(short)} sequences over {len({r['gene'] for r in short})} genes")

    # Rebuild the sequences, the same way the universe was built.
    uniprot = load_uniprot_proteoforms(cache / UNIPROT_HUMAN)
    pep = parse_ensembl_pep(cache / "Homo_sapiens.GRCh38.pep.all.fa.gz")
    ensembl_by_gene = translations_by_gene(cache / "Homo_sapiens.GRCh38.pep.all.fa.gz")
    appris = parse_appris(cache / f"appris_homo_sapiens_{APPRIS_RELEASE}.txt")
    appris_by_symbol: dict[str, list[dict]] = defaultdict(list)
    for rs in appris.values():
        for r in rs:
            appris_by_symbol[r["symbol"]].append(r)
    ensg_for_acc: dict[str, set[str]] = defaultdict(set)
    with (repo / "data/external/uniprot_ensembl_xrefs/ensg_to_uniprot.tsv").open() as fh:
        next(fh)
        for line in fh:
            g, a = line.split("\t")[:2]
            ensg_for_acc[a].add(g)

    seq_for: dict[tuple[str, str, int], str] = {}
    for gene in {r["gene"] for r in short}:
        for p in uniprot.get(gene, []):
            seq_for[(gene, p["accession"], len(p["seq"]))] = p["seq"]
        accs = [p["accession"] for p in uniprot.get(gene, []) if not p["is_isoform"]]
        ensgs = sorted({g for a in accs for g in ensg_for_acc.get(a, ())}) or sorted(
            {r["ensg"] for r in appris_by_symbol.get(gene, [])}
        )
        for g in ensgs:
            for s in ensembl_by_gene.get(g, ()):
                seq_for.setdefault((gene, "", len(s)), s)

    out_fa, manifest, missing = [], [], 0
    for i, r in enumerate(short, 1):
        seq = seq_for.get((r["gene"], r["uniprot"], int(r["length"])))
        if seq is None and r["ensembl_transcripts"]:
            seq = next(
                (pep[e.strip()] for e in r["ensembl_transcripts"].split(";") if e.strip() in pep),
                None,
            )
        if seq is None or len(seq) != int(r["length"]):
            missing += 1
            continue
        name = f"form{i:05d}|{r['gene']}|{r['uniprot'] or (r['ensembl_transcripts'].split(';')[0])}"
        out_fa.append(f">{name}\n{seq}\n")
        manifest.append(
            {
                "fasta_id": name,
                "gene": r["gene"],
                "uniprot": r["uniprot"],
                "ensembl_transcripts": r["ensembl_transcripts"],
                "length": r["length"],
                "source": r["source"],
                "tier": r["tier"],
                "openprot_te_score": r["openprot_te_score"],
                "trifid_score": r["trifid_score"],
                "appris_label_class": r["appris_label_class"],
                "is_reference": r["is_reference"],
                "reason": (
                    "reference_needed_for_comparison"
                    if r["is_reference"] == "True"
                    else "evidence+trifid"
                    if evidenced(r) and (trifid(r) or 0) > TRIFID_CUT
                    else "evidence"
                    if evidenced(r)
                    else "trifid"
                ),
            }
        )

    Path("topology_shortlist.fasta").write_text("".join(out_fa))
    with open("topology_shortlist.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(manifest[0]), delimiter="\t")
        w.writeheader()
        w.writerows(manifest)
    print(f"wrote topology_shortlist.fasta ({len(out_fa)} sequences; {missing} unresolved)")
    print(
        json.dumps(
            {
                "shortlisted": len(manifest),
                "genes": len({m["gene"] for m in manifest}),
                "reason": dict(Counter(m["reason"] for m in manifest).most_common()),
                "source": dict(Counter(m["source"] for m in manifest).most_common()),
                "residues_total": sum(int(m["length"]) for m in manifest),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
