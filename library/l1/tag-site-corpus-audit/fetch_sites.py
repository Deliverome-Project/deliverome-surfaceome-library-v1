"""Enumerate the literature-retrieved internal tag sites behind the S5E13 ranker.

Inputs
------
* ``site_features.tsv`` from branch ``claude/protein-tag-location-cv-83085b``
  (PR 180; gitignored on main). Pass its path with ``--site-features``;
  the 176 ``label == 1`` rows are the training positives.
* The public API ``https://api.deliverome.org/surfaceome/v1/tag-sites/{SYMBOL}``
  (no credentials) for every gene in those rows.
* UniProt REST for the human canonical sequence of every accession, to check
  each row's ``residue_before`` / ``residue_after``.

Outputs (in ``--out``)
----------------------
* ``api_rows.json``     every live ``literature_retrieved`` + ``internal`` row, verbatim.
                        Reused on rerun when present (delete it to refetch).
* ``sites.tsv``         one row per live site with the fields the audit needs, an
                        ``in_training`` flag, the UniProt flanking check, the first
                        PMID and citation.
* ``sequences.fasta``   the canonical sequences used.

Run::

    uv run --locked python library/l1/tag-site-corpus-audit/fetch_sites.py \
        --site-features /path/to/site_features.tsv \
        --out library/l1/tag-site-corpus-audit/data
"""

from __future__ import annotations

import argparse
import csv
import json
import time
import urllib.request
from pathlib import Path

API = "https://api.deliverome.org/surfaceome/v1/tag-sites/{}"
UNIPROT = "https://rest.uniprot.org/uniprotkb/{}.fasta"
# UniProt REST returned HTTP 500 for a subset of accessions on 2026-10-29; the EBI
# Proteins API serves the same canonical entry and is tried first, REST as fallback.
EBI_PROTEINS = "https://www.ebi.ac.uk/proteins/api/proteins/{}"


def get(url: str, retries: int = 6) -> bytes:
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "deliverome-analysis tag-site audit"})
            with urllib.request.urlopen(req, timeout=60) as r:
                return r.read()
        except Exception:  # noqa: BLE001 - retry any transport error
            if attempt == retries - 1:
                raise
            time.sleep(2 * (attempt + 1))
    raise RuntimeError(url)


def canonical_sequence(acc: str) -> str:
    try:
        entry = json.loads(get(EBI_PROTEINS.format(acc), retries=2))
        assert entry["accession"] == acc, (acc, entry.get("accession"))
        return entry["sequence"]["sequence"]
    except Exception:  # noqa: BLE001
        text = get(UNIPROT.format(acc)).decode()
        return "".join(line.strip() for line in text.splitlines() if not line.startswith(">"))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--site-features", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    with args.site_features.open() as fh:
        positives = [r for r in csv.DictReader(fh, delimiter="\t") if r["label"] == "1"]
    training_ids = {r["site_id"] for r in positives}
    genes = sorted({r["gene_symbol"] for r in positives})
    print(f"{len(positives)} training positives across {len(genes)} genes")

    api_cache = args.out / "api_rows.json"
    api_rows = json.loads(api_cache.read_text()) if api_cache.exists() else []
    for g in ([] if api_rows else genes):
        payload = json.loads(get(API.format(g)))
        for s in payload.get("sites", []):
            if s.get("provenance") == "literature_retrieved" and s.get("site_kind") == "internal":
                api_rows.append(s)
        time.sleep(0.2)
    api_cache.write_text(json.dumps(api_rows, indent=1))
    live_ids = {s["site_id"] for s in api_rows}
    print(f"{len(api_rows)} live literature internal rows; "
          f"{len(training_ids & live_ids)} of the 176 training ids found live; "
          f"missing: {sorted(training_ids - live_ids)}")

    fasta = args.out / "sequences.fasta"
    seqs: dict[str, str] = {}
    if fasta.exists():
        for block in fasta.read_text().split(">")[1:]:
            head, *lines = block.splitlines()
            seqs[head.strip()] = "".join(lines)
    for acc in sorted({s["uniprot_acc"] for s in api_rows} - set(seqs)):
        seqs[acc] = canonical_sequence(acc)
        time.sleep(0.5)
    with (args.out / "sequences.fasta").open("w") as fh:
        for acc, seq in seqs.items():
            fh.write(f">{acc}\n{seq}\n")

    cols = ["site_id", "gene", "uniprot_acc", "in_training", "insert_after_residue",
            "residue_before", "residue_after", "uniprot_before", "uniprot_after",
            "flank_matches_uniprot", "protein_length", "tag_type", "evidence_type",
            "position_evidence", "validation_level", "confidence", "entailment_verified",
            "position_frame", "cited_tag_residue", "n_sources", "pmid", "citation",
            "all_pmids", "rationale"]
    with (args.out / "sites.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for s in sorted(api_rows, key=lambda x: (x["gene_symbol"], x["insert_after_residue"] or -1)):
            seq = seqs[s["uniprot_acc"]]
            n = s["insert_after_residue"]  # None for region-only rows (never in training)
            ub = seq[n - 1] if n is not None and 1 <= n <= len(seq) else ""
            ua = seq[n] if n is not None and 0 <= n < len(seq) else ""
            src = s.get("sources") or []
            w.writerow({
                "site_id": s["site_id"], "gene": s["gene_symbol"], "uniprot_acc": s["uniprot_acc"],
                "in_training": int(s["site_id"] in training_ids),
                "insert_after_residue": n, "residue_before": s.get("residue_before"),
                "residue_after": s.get("residue_after"), "uniprot_before": ub, "uniprot_after": ua,
                "flank_matches_uniprot": "" if n is None else int(ub == s.get("residue_before") and ua == s.get("residue_after")),
                "protein_length": len(seq), "tag_type": s.get("tag_type"),
                "evidence_type": s.get("evidence_type"), "position_evidence": s.get("position_evidence"),
                "validation_level": s.get("validation_level"), "confidence": s.get("confidence"),
                "entailment_verified": s.get("entailment_verified"), "position_frame": s.get("position_frame"),
                "cited_tag_residue": s.get("cited_tag_residue"), "n_sources": len(src),
                "pmid": src[0].get("pmid") if src else "", "citation": src[0].get("citation") if src else "",
                "all_pmids": ";".join(str(x.get("pmid")) for x in src),
                "rationale": (s.get("rationale") or "").replace("\t", " ").replace("\n", " "),
            })


if __name__ == "__main__":
    main()
