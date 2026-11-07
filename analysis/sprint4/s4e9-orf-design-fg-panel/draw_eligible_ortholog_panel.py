"""Seeded review-only ortholog draw after full-pool CDS and prediction checks.

Freeze the published panel's species/band/class quotas and measured band cutoffs.
Sort candidates explicitly before shuffling; this is a new deterministic full draw,
not a claim to reproduce the old iteration order. Pending winning predictions block
selection. No quota borrowing, boundary repairs, ID allocation or activation occurs.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import random
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd
import rebuild_ortholog_pool as pool
from draw_rebuilt_ortholog_panel import classify
from export_ortholog_predictions import fasta

from deliverome_analysis import library_design as ld
from deliverome_analysis import seqtools, topology

RULE_VERSION = "ortholog-cds-draw-v1"


def read(path):
    with Path(path).open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def truth(value):
    return value is True or value == "True"


def write(path, rows, columns):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def frozen_policy(old_panel, old_summary):
    quotas = Counter((r["species"], r["band"], r["topology_class"]) for r in old_panel)
    if len(old_panel) != 200 or any(
        sum(n for (s, b, _), n in quotas.items() if (s, b) == (species, band)) != 50
        for species in ("mouse", "cyno")
        for band in ("near", "diverged")
    ):
        raise ValueError("Published quota source must contain four 50-member bands")
    return {
        "seed": old_summary["seed"],
        "cutoffs": {
            s: old_summary["species"][s]["correspondence_split"] for s in ("mouse", "cyno")
        },
        "quotas": {"|".join(k): n for k, n in sorted(quotas.items())},
        "order": "species mouse,cyno; band near,diverged; class lexical; candidates stable human gene,target gene,protein hash,accession; per-stratum shuffle; consume first globally unique proteins",
        "band_score_precision": "published picked_score (four decimal places)",
    }


def draw(rows, fixed, policy):
    """Return a tentative full draw and explicit shortfalls, never relax quotas."""
    groups = defaultdict(list)
    identities = set()
    for row in rows:
        identity = (row["species"], row["human_ensembl_gene"])
        if identity in identities:
            raise ValueError(f"Multiple representatives for species/human gene: {identity}")
        identities.add(identity)
        groups[(row["species"], row["band"], row["topology_class"])].append(row)
    selected, rejected, deficits = [], [], []
    used = dict(fixed)
    for species in ("mouse", "cyno"):
        rng = random.Random(policy["seed"] + (species == "cyno"))
        for band in ("near", "diverged"):
            classes = sorted(
                k.split("|")[2] for k in policy["quotas"] if k.startswith(f"{species}|{band}|")
            )
            for cls in classes:
                key = (species, band, cls)
                target = policy["quotas"]["|".join(key)]
                candidates = sorted(
                    groups[key],
                    key=lambda r: (
                        r["human_ensembl_gene"],
                        r["ortholog_ensembl_gene"],
                        r["sequence_sha256"],
                        r["picked_acc"],
                    ),
                )
                rng.shuffle(candidates)
                picked = 0
                for row in candidates:
                    h = row["sequence_sha256"]
                    if h == row["human_reference_sha256"]:
                        rejected.append(
                            {"form_id": row["form_id"], "reason": "exact_human_reference_protein"}
                        )
                        continue
                    if h in used:
                        rejected.append(
                            {"form_id": row["form_id"], "reason": "protein_duplicates:" + used[h]}
                        )
                        continue
                    if picked < target:
                        selected.append(row)
                        used[h] = row["form_id"]
                        picked += 1
                if picked != target:
                    deficits.append(
                        {
                            "species": species,
                            "band": band,
                            "topology_class": cls,
                            "required": target,
                            "available_after_uniqueness": picked,
                            "shortfall": target - picked,
                        }
                    )
    return selected, rejected, deficits


def validate_encoding(encoding, sequence, species, gene):
    """Recheck the exact material to be emitted, not only its descriptive fields."""
    dna = encoding["cds"]
    if (encoding["species"], encoding["gene"]) != (species, gene):
        raise ValueError("Representative encoding species/gene mismatch")
    if encoding["eligibility_status"] != "eligible" or not encoding["chosen_encoding"]:
        raise ValueError("Representative encoding is not an eligible selected CDS")
    if (
        not encoding["transcript"]
        or encoding["transcript_version"].split(".")[0] != encoding["transcript"]
    ):
        raise ValueError("Representative transcript identity mismatch")
    if not dna or set(dna) - set("ACGT") or len(dna) % 3 or dna[-3:] not in {"TAA", "TAG", "TGA"}:
        raise ValueError("Representative CDS fails mechanical integrity")
    translated = seqtools.translate(dna)
    if "*" in translated[:-1] or translated[:-1] != sequence or encoding["protein"] != sequence:
        raise ValueError("Representative CDS does not encode scored protein")
    if (
        hashlib.sha256(dna.encode()).hexdigest() != encoding["cds_sha256"]
        or hashlib.sha256(sequence.encode()).hexdigest() != encoding["protein_sha256"]
    ):
        raise ValueError("Representative sequence checksum mismatch")


def classify_human(cds, gpi, snapshot):
    human = {
        k[4:]: seqtools.translate(v["cds"]).rstrip("*")
        for k, v in cds.items()
        if k.startswith("lib:")
    }
    report = topology.resolve(
        human,
        tool=topology.TOOL_V1,
        cohorts=topology.cohorts_for_species("human"),
        snapshot=pd.read_parquet(snapshot),
    )
    result = {}
    for row in report.to_dict("records"):
        if not row["found"]:
            continue
        first = next((c for c in row["per_residue_topology"] if c in "IO"), "")
        result[row["key"]] = classify(
            int(row["tm_helix_count"]),
            int(row["signal_peptide_length"]),
            {"I": "cytoplasmic", "O": "extracellular"}.get(first, "indeterminate"),
            row["key"] in gpi,
        )
    return result


def run(args):
    if args.out.exists():
        raise FileExistsError(args.out)
    inputs = {
        name: getattr(args, name)
        for name in (
            "old_panel",
            "old_summary",
            "cds",
            "constructs",
            "manifest",
            "gpi",
            "topology",
            "eligible",
        )
    }
    inputs.update(
        {
            "scores": args.scored / "ortholog_provisional_picks.tsv",
            "score_summary": args.scored / "ortholog_scoring_summary.json",
            "proteins": args.scored / "ortholog_picks.fasta",
        }
    )
    old_panel = read(args.old_panel)
    policy = frozen_policy(old_panel, json.loads(args.old_summary.read_text()))
    native = json.loads(args.cds.read_text())
    constructs = json.loads(args.constructs.read_text())
    manifest = read(args.manifest)
    encodings = json.loads(args.eligible.read_text())
    gpi = {r["gene_symbol"] for r in read(args.gpi) if r["gpi_anchored"] == "1"}
    classes = classify_human(native, gpi, args.topology)
    fixed, existing = {}, defaultdict(list)
    for row in manifest:
        key = row["key"]
        dna = native.get(key, {}).get("cds") or constructs.get(key, {}).get("construct")
        if not dna:
            if row["library"] == "l1" and ld.is_active(row):
                raise ValueError(f"Active L1 member missing protein: {key}")
            continue
        h = hashlib.sha256(seqtools.translate(dna).rstrip("*").encode()).hexdigest()
        if row["set"] == "ortholog":
            species, gene = key.split(":")[1:3]
            existing[(species, gene, h)].append(row["candidate_id"])
        elif row["library"] == "l1" and ld.is_active(row):
            fixed[h] = row["candidate_id"]
    sequences = fasta(inputs["proteins"])
    scores = read(inputs["scores"])
    score_summary = json.loads(inputs["score_summary"].read_text())
    candidates, pending, excluded = [], [], []
    for row in scores:
        species, gene, acc = row["species"], row["gene"], row["picked_acc"]
        state = row["selection_state"]
        if state in {
            "needs_netgpi_and_v1",
            "needs_netgpi",
            "needs_deeptmhmm_v1",
            "needs_netgpi_review_and_v1",
        }:
            pending.append(
                {
                    "species": species,
                    "gene": gene,
                    "selection_state": state,
                    "sequence_sha256": row["picked_sequence_sha256"],
                }
            )
        if not truth(row["sampling_eligible"]):
            excluded.append({"species": species, "gene": gene, "reason": state})
            continue
        sequence = sequences[f"{species}|{gene}|{acc}"]
        h = hashlib.sha256(sequence.encode()).hexdigest()
        if h != row["picked_sequence_sha256"]:
            raise ValueError("Scored sequence checksum mismatch")
        encoding = encodings[species][h]["by_gene"][row["ortholog_ensembl_gene"]]
        validate_encoding(encoding, sequence, species, row["ortholog_ensembl_gene"])
        human = seqtools.translate(native[f"lib:{gene}"]["cds"]).rstrip("*")
        if hashlib.sha256(human.encode()).hexdigest() != row["human_reference_sha256"]:
            raise ValueError("Scored human reference differs from supplied native CDS")
        if gene not in classes:
            raise ValueError(f"Missing exact human reference topology: {gene}")
        band = "near" if float(row["picked_score"]) >= policy["cutoffs"][species] else "diverged"
        candidates.append(
            {
                "species": species,
                "human_symbol": gene,
                "ortholog_symbol": row["ortholog_symbols"],
                "picked_acc": acc,
                "picked_score": row["picked_score"],
                "topology_label": row["picked_topology_label"],
                "topology_tool": topology.TOOL_V1,
                "anchor_call": row["picked_anchor_call"],
                "topology_class": classes[gene],
                "form_id": f"{species}:{gene}:{acc}",
                "human_ensembl_gene": row["human_ensembl_gene"],
                "ortholog_ensembl_gene": row["ortholog_ensembl_gene"],
                "sequence_sha256": h,
                "human_reference_sha256": row["human_reference_sha256"],
                "ensembl_transcripts": encoding["transcript"],
                "existing_candidate_ids": ";".join(sorted(existing[(species, gene, h)])),
                "encoding_source": "eligibility_ranked_ensembl116_CDS",
                "in_april_ortholog_list": row["in_april_ortholog_list"],
                "band": band,
                "selected_cds_sha256": hashlib.sha256(encoding["cds"].encode()).hexdigest(),
                "cds_eligibility_status": encoding["eligibility_status"],
            }
        )
    args.out.mkdir(parents=True)
    ready = bool(score_summary["panel_ready"]) and not pending
    selected, rejected, deficits = draw(candidates, fixed, policy) if ready else ([], [], [])
    status = (
        "blocked_pending_predictions"
        if not ready
        else "blocked_stratum_shortfall"
        if deficits
        else "review_panel_ready"
    )
    if status == "review_panel_ready":
        if len(selected) != 200 or len({r["sequence_sha256"] for r in selected}) != 200:
            raise ValueError("Panel failed size/uniqueness check")
        selected.sort(key=lambda r: (r["species"], r["band"], r["human_symbol"]))
        pool.write_tsv(args.out / "ortholog_panel_cds_eligible.tsv", selected)
        selected_cds = {
            r["form_id"]: encodings[r["species"]][r["sequence_sha256"]]["by_gene"][
                r["ortholog_ensembl_gene"]
            ]
            for r in selected
        }
        (args.out / "selected_cds.json").write_text(json.dumps(selected_cds, indent=2) + "\n")
    write(
        args.out / "pending_predictions.tsv",
        pending,
        ["species", "gene", "selection_state", "sequence_sha256"],
    )
    write(args.out / "sampling_exclusions.tsv", excluded, ["species", "gene", "reason"])
    write(args.out / "draw_rejections.tsv", rejected, ["form_id", "reason"])
    write(
        args.out / "stratum_shortfalls.tsv",
        deficits,
        [
            "species",
            "band",
            "topology_class",
            "required",
            "available_after_uniqueness",
            "shortfall",
        ],
    )
    old = {(r["species"], r["human_ensembl_gene"]): r for r in old_panel}
    new = (
        {(r["species"], r["human_ensembl_gene"]): r for r in selected}
        if status == "review_panel_ready"
        else {}
    )
    changes = []
    if new:
        for key in sorted(old.keys() | new.keys()):
            before, after = old.get(key, {}), new.get(key, {})
            same = before.get("sequence_sha256") == after.get("sequence_sha256")
            changes.append(
                {
                    "species": key[0],
                    "human_ensembl_gene": key[1],
                    "human_symbol": (after or before)["human_symbol"],
                    "change": "retained_protein"
                    if before and after and same
                    else "changed_form"
                    if before and after
                    else "added_gene"
                    if after
                    else "removed_gene",
                    "old_protein_sha256": before.get("sequence_sha256", ""),
                    "new_protein_sha256": after.get("sequence_sha256", ""),
                    "old_transcripts": before.get("ensembl_transcripts", ""),
                    "new_transcript": after.get("ensembl_transcripts", ""),
                    "decision": "eligible_pool_full_seeded_draw_with_frozen_quotas",
                }
            )
    write(
        args.out / "panel_changes.tsv",
        changes,
        [
            "species",
            "human_ensembl_gene",
            "human_symbol",
            "change",
            "old_protein_sha256",
            "new_protein_sha256",
            "old_transcripts",
            "new_transcript",
            "decision",
        ],
    )
    result = {
        "rule_version": RULE_VERSION,
        "status": status,
        "policy": policy,
        "sampleable_candidates": len(candidates),
        "pending_winners": len(pending),
        "selected": len(selected) if status == "review_panel_ready" else 0,
        "shortfalls": deficits,
        "changes": dict(Counter(r["change"] for r in changes)),
        "canonical_activation": False,
        "input_sha256": {name: digest(path) for name, path in inputs.items()},
        "command_argv": sys.argv,
        "script_sha256": digest(__file__),
        "source_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
        "source_dirty": bool(
            subprocess.check_output(["git", "status", "--porcelain"], text=True).strip()
        ),
    }
    result["output_sha256"] = {p.name: digest(p) for p in sorted(args.out.iterdir())}
    (args.out / "ortholog_eligible_panel_summary.json").write_text(
        json.dumps(result, indent=2) + "\n"
    )
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in (
        "old-panel",
        "old-summary",
        "cds",
        "constructs",
        "manifest",
        "gpi",
        "topology",
        "eligible",
        "scored",
        "out",
    ):
        parser.add_argument("--" + name, type=Path, required=True)
    print(json.dumps(run(parser.parse_args()), indent=2))


if __name__ == "__main__":
    main()
