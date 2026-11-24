"""Score all release-specific ortholog proteoforms against the actual human designs.

This prepares the expanded pool for the seeded panel draw. It uses the existing
SequenceMatcher correspondence metric (matching blocks / longer sequence length),
the 0.001 tie tolerance, and reviewed / APPRIS / current-pick tie preferences.
The final lexical accession tiebreak makes otherwise equal new candidates stable.

Unknown topology is a pending prediction, never a reason to exclude a gene from
the rebuilt pool. This script writes an explicit prediction worklist and refuses
to present the pool as ready for sampling while that work remains. GPI candidates
also need sequence-verified NetGPI calls before the anchor guard can be finalized.

UniProt candidates retain the established symbol-based within-gene lookup; the
Ensembl candidates come directly from Compara's target stable gene ID. UniProt
gene-to-accession provenance is the separate review issue #19.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import sys
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import pandas as pd
import rebuild_ortholog_pool as pool

HERE = Path(__file__).parent
AP = HERE.parents[2] / "library/l1/proteoform-selection"
sys.path.insert(0, str(AP))

from repick_orthologs import score  # noqa: E402
from select_principal import parse_appris, parse_uniprot_fasta, select_principal  # noqa: E402

from deliverome_analysis import seqtools, topology  # noqa: E402


def read_tsv(path):
    return list(csv.DictReader(path.open(), delimiter="\t"))


def best_candidate(job):
    """Score once per distinct sequence; retain all provenance aliases in ties."""
    row, reference, candidates, principal, current = job
    by_sequence = {c["sequence"] for c in candidates.values()}
    scores = {s: score(reference, s) for s in sorted(by_sequence)}
    # Preserve the existing soluble-form guard: an unknown is not a negative call.
    allowed = {a: c for a, c in candidates.items() if c["surface_ok"] is not False}
    effective = allowed or candidates
    best = max(scores[c["sequence"]] for c in effective.values())
    tied = [a for a, c in effective.items() if best - scores[c["sequence"]] <= 0.001]
    picked = min(
        tied,
        key=lambda a: (
            not candidates[a]["reviewed"],
            not principal or candidates[a]["sequence"] != principal,
            a != current,
            a,
        ),
    )
    candidate = candidates[picked]
    sequence = candidate["sequence"]
    needs_gpi = row["human_is_gpi"] and not candidate.get("anchor_observed", False)
    uncertain_gpi = (
        row["human_is_gpi"]
        and candidate.get("anchor_observed", False)
        and candidate.get("anchor_ok") is None
    )
    return {
        **row,
        "picked_acc": picked,
        "picked_source": candidate["source"],
        "picked_reviewed": candidate["reviewed"],
        "picked_score": round(scores[sequence], 4),
        "picked_length": len(sequence),
        "picked_sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
        "human_reference_sha256": hashlib.sha256(reference.encode()).hexdigest(),
        "picked_topology_label": candidate["label"],
        "picked_anchor_call": {True: "anchored", False: "not_anchored"}.get(
            candidate.get("anchor_ok"), ""
        ),
        "topology_tool": topology.TOOL_V1,
        "anchor_uncertain": uncertain_gpi,
        "sampling_eligible": (
            not needs_gpi
            and candidate["surface_ok"] is not False
            and bool(candidate["label"])
            and (
                candidate["label"] in {"TM", "SP+TM", "BETA"} or candidate.get("anchor_ok") is True
            )
        ),
        "n_candidates": len(candidates),
        "n_tied": len(tied),
        "n_candidates_unknown_v1": sum(not c["label"] for c in candidates.values()),
        "selection_state": (
            "needs_netgpi_and_v1"
            if needs_gpi and not candidate["label"]
            else "needs_netgpi"
            if needs_gpi
            else "needs_deeptmhmm_v1"
            if not candidate["label"]
            else "known_non_surface"
            if candidate["surface_ok"] is False
            else "provisional_surface_pick"
        ),
        "sequence": sequence,
    }


def read_v1_supplement(path):
    """Validate the submitted sequence digest and every predicted residue."""
    if path is None:
        return {}
    lines = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    if len(lines) % 3:
        raise ValueError("v1 supplement must contain complete three-line records")
    result = {}
    for i in range(0, len(lines), 3):
        header, sequence, topo = lines[i : i + 3]
        if not header.startswith(">") or "|" not in header:
            raise ValueError("v1 supplement header lacks its label")
        digest, label = header[1:].rsplit("|", 1)
        digest = digest.strip()
        if digest != hashlib.sha256(sequence.encode()).hexdigest():
            raise ValueError("v1 supplement sequence does not match submitted digest")
        if len(topo) != len(sequence) or set(topo) - set("IOMSB"):
            raise ValueError("v1 supplement topology does not match the sequence")
        if digest in result:
            raise ValueError("duplicate v1 supplement sequence")
        result[digest] = label.strip()
    return result


def read_netgpi(path):
    if path is None:
        return {}
    result = {}
    for line in path.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 6 or len(fields[0]) != 64:
            raise ValueError("NetGPI result lacks a full submitted sequence digest")
        if fields[0] in result:
            raise ValueError("duplicate NetGPI result")
        positive = fields[2] == "GPI-Anchored"
        call = positive if positive or float(fields[4]) >= 0.55 else None
        result[fields[0]] = (int(fields[1]), call)
    return result


def run(inputs, pool_dir, out, workers, v1_supplement=None, netgpi=None):
    spec = pool.verify_inputs(inputs)
    pool_summary = json.loads((pool_dir / "ortholog_pool_summary.json").read_text())
    for name, digest in pool_summary["input_sha256"].items():
        if digest != spec["files"][name]["sha256"]:
            raise ValueError(f"pool was built from different input: {name}")
    cds = json.loads((inputs / "candidate_cds.json").read_text())
    supplemental_labels = read_v1_supplement(v1_supplement)
    anchor_calls = read_netgpi(netgpi)
    genes = json.loads((inputs / "library_members.json").read_text())
    references = {
        g: seqtools.translate(cds[f"lib:{g}"]["cds"]).rstrip("*")
        for g in genes
        if f"lib:{g}" in cds
    }
    gpi_rows = read_tsv(inputs / "gpi_set.tsv")
    gpi = {r["gene_symbol"] for r in gpi_rows if r["gpi_anchored"] == "1"}
    pairs = read_tsv(pool_dir / "ortholog_pool.tsv")
    jobs, missing = [], []
    # Resolve exact sequences through the shared v1 API; it validates predictor provenance
    # and keeps species cohorts separate. No accession/length fallback is permitted.
    for species in pool.SPECIES:
        proteins = pool.read_peptides(inputs / pool.PEPTIDES[species])
        by_gene = defaultdict(dict)
        for r in proteins:
            if r["primary_assembly"]:
                by_gene[r["gene"]][r["protein"]] = {
                    "sequence": r["sequence"],
                    "source": "ensembl",
                    "reviewed": False,
                }
        uniprot_name = (
            "uniprot_mouse_reviewed_isoforms.tsv.gz"
            if species == "mouse"
            else "uniprot_cyno_all_isoforms.tsv.gz"
        )
        up = defaultdict(dict)
        with gzip.open(inputs / uniprot_name, "rt") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                symbol = r.get("Gene Names (primary)", "").strip()
                if symbol:
                    up[symbol][r["Entry"]] = {
                        "sequence": r["Sequence"].rstrip("*"),
                        "source": "uniprot",
                        "reviewed": species == "mouse"
                        or r.get("Reviewed", "").lower().startswith("reviewed"),
                    }
        old_meta = {
            r["query_ensembl_gene_id"]: r
            for r in csv.DictReader(
                (inputs / f"{species}_ortholog_one2one_highconf_non_hla_metadata.csv").open()
            )
            if r.get("status") == "ok"
        }
        cohort = parse_uniprot_fasta(inputs / f"{species}_ortholog_one2one_highconf_non_hla.fasta")
        principal = {}
        if species == "mouse":
            appris = parse_appris(inputs / "appris_mus_musculus_e115v50.txt")
            pep = {r["transcript"]: r["sequence"] for r in proteins}
            for gene, rows in appris.items():
                selected = select_principal(rows, pep)
                if selected and selected["seq"]:
                    principal[gene] = selected["seq"]
        species_jobs = []
        for row in pairs:
            if row["species"] != species or row["status"] != "eligible":
                continue
            gene = row["gene"]
            if gene not in references:
                missing.append({**row, "selection_state": "human_designed_cds_missing"})
                continue
            target = row["ortholog_ensembl_gene"]
            candidates = dict(by_gene[target])
            for symbol in row["ortholog_symbols"].split(";"):
                candidates.update(up.get(symbol, {}))
            old = old_meta.get(row["human_ensembl_gene"], {})
            current = ""
            if old.get("target_ensembl_gene_id") == target:
                current = old.get("selected_uniprot_accession", "")
                if current in cohort:
                    candidates.setdefault(
                        current,
                        {
                            "sequence": cohort[current]["seq"],
                            "source": "current",
                            "reviewed": old.get("selected_reviewed") == "reviewed",
                        },
                    )
            species_jobs.append(
                (
                    {**row, "human_is_gpi": gene in gpi},
                    references[gene],
                    candidates,
                    principal.get(target),
                    current,
                )
            )
        seqs = {
            hashlib.sha256(c["sequence"].encode()).hexdigest(): c["sequence"]
            for j in species_jobs
            for c in j[2].values()
        }
        report = topology.resolve(
            seqs,
            tool=topology.TOOL_V1,
            cohorts=topology.cohorts_for_species(species),
            snapshot=pd.read_parquet(inputs / "topology.parquet"),
        )
        labels = {
            r["key"]: r.get("deeptmhmm_label", "") if r["found"] else ""
            for r in report.to_dict("records")
        }
        for row, ref, candidates, principal_seq, current in species_jobs:
            for candidate in candidates.values():
                digest = hashlib.sha256(candidate["sequence"].encode()).hexdigest()
                label = supplemental_labels.get(digest, labels[digest])
                candidate["label"] = label
                anchor = None
                if row["human_is_gpi"] and digest in anchor_calls:
                    length, anchor = anchor_calls[digest]
                    if length != len(candidate["sequence"]):
                        raise ValueError("NetGPI sequence length does not match submitted digest")
                candidate["anchor_ok"] = anchor
                candidate["anchor_observed"] = row["human_is_gpi"] and digest in anchor_calls
                candidate["surface_ok"] = (
                    anchor
                    if anchor is not None
                    else (
                        None
                        if not label
                        else label in {"TM", "SP+TM", "BETA"} or row["human_is_gpi"]
                    )
                )
            jobs.append((row, ref, candidates, principal_seq, current))
    out.mkdir(parents=True, exist_ok=True)
    with ProcessPoolExecutor(max_workers=workers) as executor:
        selected = list(executor.map(best_candidate, jobs, chunksize=10))
    # Retain the full unknown inventory, but predict winners first and rescore after
    # each batch. Predicting all alternatives would spend compute on losing forms.
    unknown = {}
    for row, _ref, candidates, _principal, _current in jobs:
        for acc, c in candidates.items():
            if not c["label"]:
                digest = hashlib.sha256(c["sequence"].encode()).hexdigest()
                unknown.setdefault(digest, (f"{row['species']}|{row['gene']}|{acc}", c["sequence"]))
    with (out / "all_unknown_candidate_sequences.fasta").open("w") as fh:
        for digest, (name, seq) in sorted(unknown.items()):
            fh.write(f">{name}|sha256={digest}\n{seq}\n")
    winner_v1 = {r["picked_sequence_sha256"]: r for r in selected if not r["picked_topology_label"]}
    winner_gpi = {
        r["picked_sequence_sha256"]: r for r in selected if "netgpi" in r["selection_state"]
    }
    for name, worklist in (("missing_deeptmhmm_v1", winner_v1), ("pending_netgpi", winner_gpi)):
        with (out / f"{name}.fasta").open("w") as fh:
            for digest, r in sorted(worklist.items()):
                fh.write(f">{digest}\n{r['sequence']}\n")
    with (out / "ortholog_picks.fasta").open("w") as fh:
        for r in selected:
            fh.write(f">{r['species']}|{r['gene']}|{r['picked_acc']}\n{r['sequence']}\n")
    pool.write_tsv(
        out / "ortholog_provisional_picks.tsv",
        [{k: v for k, v in r.items() if k != "sequence"} for r in selected],
    )
    if missing:
        pool.write_tsv(out / "unresolved_references.tsv", missing)
    summary = {
        "input_sha256": {n: e["sha256"] for n, e in spec["files"].items()},
        "supplement_sha256": {
            p.name: hashlib.sha256(p.read_bytes()).hexdigest()
            for p in (v1_supplement, netgpi)
            if p is not None
        },
        "scored_pairs": len(selected),
        "missing_human_references": len(missing),
        "state_by_species": {
            s: dict(Counter(r["selection_state"] for r in selected if r["species"] == s))
            for s in pool.SPECIES
        },
        "unknown_candidate_sequences": len(unknown),
        "unknown_candidate_residues": sum(len(v[1]) for v in unknown.values()),
        "first_batch_deeptmhmm_v1_sequences": len(winner_v1),
        "first_batch_deeptmhmm_v1_residues": sum(len(r["sequence"]) for r in winner_v1.values()),
        "first_batch_netgpi_sequences": len(winner_gpi),
        "panel_ready": not winner_v1 and not winner_gpi,
        "eligible_for_sampling": dict(
            Counter(r["species"] for r in selected if r["sampling_eligible"])
        ),
        "uncertain_anchor_picks": sum(r["anchor_uncertain"] for r in selected),
        "reason": "complete missing v1 and GPI anchor calls before sampling; every excluded or unresolved pair remains in the accounting tables",
    }
    (out / "ortholog_scoring_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--pool", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--v1-supplement", type=Path)
    parser.add_argument("--netgpi", type=Path)
    args = parser.parse_args()
    print(
        json.dumps(
            run(args.inputs, args.pool, args.out, args.workers, args.v1_supplement, args.netgpi),
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
