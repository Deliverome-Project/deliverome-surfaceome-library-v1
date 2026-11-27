"""Replay the historical symbol-to-accession rule against immutable UniProt tables.

This audit never chooses library replacements or writes library/construct records.
The original API responses were not saved; confirmations mean the same accession
wins in the pinned release, not that an original response has been recovered.
"""

from __future__ import annotations

import argparse
import collections
import csv
import gzip
import hashlib
import json
from pathlib import Path


def normalize(value: str) -> str:
    return (value or "").strip().upper()


def rank(row: dict[str, str], symbol: str) -> tuple:
    """Exact f9c73691 historical choose_best_candidate sort key."""
    target = normalize(symbol)
    try:
        length = int((row.get("Length") or "").strip())
    except ValueError:
        length = 0
    reviewed = (row.get("Reviewed") or "").strip().lower()
    return (
        normalize(row.get("Gene Names (primary)")) == target,
        target in {normalize(t) for t in (row.get("Gene Names") or "").split()},
        reviewed in {"1", "true", "t", "yes", "y", "reviewed"},
        length,
        (row.get("Entry") or "").strip(),
    )


def read_rows(path: Path, delimiter: str) -> list[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def audit(metadata: list[dict[str, str]], entries: list[dict[str, str]], taxid: str):
    index: dict[str, list] = collections.defaultdict(list)
    accessions = {}
    for entry in entries:
        if entry["Organism (ID)"] != taxid:
            raise ValueError("Taxonomy table includes an unexpected organism")
        if entry["Entry"] in accessions:
            raise ValueError("Duplicate accession in taxonomy table")
        accessions[entry["Entry"]] = entry
        for token in {normalize(t) for t in entry["Gene Names"].split()}:
            index[token].append(entry)
    results = []
    for original in metadata:
        if original["taxid"] != taxid:
            raise ValueError("Metadata taxonomy mismatch")
        symbol = original["target_gene_symbol"]
        candidates = index.get(normalize(symbol), [])
        # The historical request used size=500 without pagination. If exceeded,
        # its undocumented server ordering cannot be recreated from this table.
        capped = len(candidates) > 500
        best = max(candidates, key=lambda row: rank(row, symbol), default={})
        old = original["selected_uniprot_accession"]
        old_entry = accessions.get(old, {})
        winner = best.get("Entry", "")
        outcome = (
            "not_recorded_ok"
            if original["status"] != "ok"
            else "historical_cap_unresolved"
            if capped
            else "confirmed_accession"
            if old == winner
            else "different_accession"
            if winner
            else "no_candidate"
        )
        results.append(
            {
                "species": original["species"],
                "query_ensembl_gene_id": original["query_ensembl_gene_id"],
                "target_ensembl_gene_id": original["target_ensembl_gene_id"],
                "target_gene_symbol": symbol,
                "recorded_status": original["status"],
                "recorded_accession": old,
                "recorded_candidate_count": original["n_candidates"],
                "release_candidate_count": len(candidates),
                "release_winner_accession": winner,
                "release_winner_reviewed": best.get("Reviewed", ""),
                "recorded_length": original["selected_length"],
                "release_winner_length": best.get("Length", ""),
                "recorded_accession_present_in_release": old in accessions,
                "recorded_release_primary_gene": old_entry.get("Gene Names (primary)", ""),
                "recorded_release_gene_names": old_entry.get("Gene Names", ""),
                "recorded_release_length": old_entry.get("Length", ""),
                "recorded_accession_in_symbol_candidates": any(
                    c["Entry"] == old for c in candidates
                ),
                "same_recorded_length": bool(
                    winner and old == winner and original["selected_length"] == best.get("Length")
                ),
                "outcome": outcome,
            }
        )
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--download-inputs", action="store_true")
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text())
    paths = {}
    for item in manifest["inputs"]:
        path = args.input_dir / item["name"]
        if not path.exists() and args.download_inputs:
            from deliverome_analysis.core.s3 import external

            path.parent.mkdir(parents=True, exist_ok=True)
            external.download(item["source"].removeprefix("external://"), path)
        if hashlib.sha256(path.read_bytes()).hexdigest() != item["sha256"]:
            raise ValueError(f"Checksum mismatch: {item['name']}")
        paths[item["name"]] = path
        if path.suffix == ".headers":
            headers = dict(
                line.split(":", 1) for line in path.read_text().lower().splitlines() if ":" in line
            )
            if headers.get("x-uniprot-release", "").strip() != manifest["uniprot_release"].lower():
                raise ValueError(f"Unexpected UniProt release: {path.name}")
    rows = []
    summary = {"release": manifest["uniprot_release"], "species": {}}
    for species, taxid in [("mouse", "10090"), ("cyno", "9541")]:
        entries = read_rows(paths[f"{species}.tsv.gz"], "\t")
        count_headers = dict(
            line.split(":", 1)
            for line in paths[f"{species}_count.headers"].read_text().lower().splitlines()
            if ":" in line
        )
        if len(entries) != int(count_headers["x-total-results"]):
            raise ValueError(f"Incomplete taxonomy snapshot: {species}")
        metadata = read_rows(paths[f"{species}_metadata.csv"], ",")
        result = audit(metadata, entries, taxid)
        rows.extend(result)
        summary["species"][species] = {
            "release_entries": len(entries),
            "recorded_rows": len(metadata),
            "outcomes": dict(sorted(collections.Counter(r["outcome"] for r in result).items())),
            "confirmed_same_length": sum(r["same_recorded_length"] for r in result),
            "max_candidates": max(r["release_candidate_count"] for r in result),
        }
        for validation in manifest["search_validations"]:
            if validation["species"] != species:
                continue
            symbol = normalize(validation["symbol"])
            offline = {
                r["Entry"]
                for r in entries
                if symbol in {normalize(t) for t in r["Gene Names"].split()}
            }
            live = {r["Entry"] for r in read_rows(paths[validation["name"]], "\t")}
            if offline != live:
                raise ValueError(f"Exact-search validation failed: {species} {symbol}")
    summary["exact_search_validations_passed"] = len(manifest["search_validations"])
    summary["total_outcomes"] = dict(
        sorted(collections.Counter(r["outcome"] for r in rows).items())
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    with (args.output_dir / "ortholog_accession_replay.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
    (args.output_dir / "ortholog_accession_replay_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
