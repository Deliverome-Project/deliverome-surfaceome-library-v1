"""Draw 200 orthologs from the rebuilt pool with the established seed and strata.

New genes do not need an old manifest entry: an exact Ensembl translation of the
Compara target gene is a native coding-sequence source. An existing designed
protein may also supply the exact form. Unbuildable forms, identical-to-human
forms and duplicate designed proteins are replaced within their original
species/band/class stratum using the established replacement rule.

Outputs are the replacement panel and its explicit changes. This does not
allocate member IDs, retire constructs or alter ordered material.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import random
import sys
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd
import rebuild_ortholog_pool as pool
from export_ortholog_predictions import fasta

from deliverome_analysis.cloning import library_design as ld
from deliverome_analysis.cloning import seqtools
from deliverome_analysis.surfaceome import topology

AP = Path(__file__).resolve().parents[3] / "library/l1/proteoform-selection"
sys.path.insert(0, str(AP))
from panel_draw import replace_rejected  # noqa: E402
from sample_isoform_panel import allocate, classify  # noqa: E402
from sample_ortholog_panel import PER_SPECIES, SEED  # noqa: E402


def read(path):
    return list(csv.DictReader(path.open(), delimiter="\t"))


def draw(rows, rejected):
    output, replacements, summary = [], [], {}
    for species in ("mouse", "cyno"):
        eligible = [r for r in rows if r["species"] == species]
        scores = sorted(float(r["picked_score"]) for r in eligible)
        split = scores[len(scores) // 2]
        bands = {
            "near": [r for r in eligible if float(r["picked_score"]) >= split],
            "diverged": [r for r in eligible if float(r["picked_score"]) < split],
        }
        rng = random.Random(SEED + (0 if species == "mouse" else 1))
        strata = {}
        for band, band_rows in bands.items():
            by_class = defaultdict(list)
            for row in band_rows:
                by_class[row["topology_class"]].append(row)
            allocation = allocate(band_rows, PER_SPECIES)
            strata[band] = (
                by_class,
                [(c, rng.sample(by_class[c], k)) for c, k in allocation.items()],
            )
        used = {
            r["form_id"]
            for _, groups in strata.values()
            for _, selected in groups
            for r in selected
        }
        for band, (by_class, groups) in strata.items():
            for cls, selected in groups:
                picked, log = replace_rejected(
                    rng,
                    selected,
                    by_class[cls],
                    lambda r: r["form_id"],
                    rejected,
                    lambda r: "protein identical to human" if float(r["picked_score"]) >= 1 else "",
                    used,
                )
                if len(picked) != len(selected):
                    raise ValueError(
                        f"cannot fill {species}/{band}/{cls}; no valid replacement remains"
                    )
                output.extend({**r, "band": band} for r in picked)
                replacements.extend(
                    {"species": species, "band": band, "topology_class": cls, **r} for r in log
                )
        summary[species] = {"eligible": len(eligible), "correspondence_split": split}
    return output, replacements, summary


def run(inputs, scored, out):
    pool.verify_inputs(inputs)
    summary = json.loads((scored / "ortholog_scoring_summary.json").read_text())
    if not summary["panel_ready"]:
        raise ValueError("complete the pending predictions before the full-pool draw")
    cds = json.loads((inputs / "candidate_cds.json").read_text())
    manifest = read(inputs / "candidate_manifest_tranches.tsv")
    gpi = {r["gene_symbol"] for r in read(inputs / "gpi_set.tsv") if r["gpi_anchored"] == "1"}
    proteins = {key: seqtools.translate(v["cds"]).rstrip("*") for key, v in cds.items()}
    human = {key[4:]: value for key, value in proteins.items() if key.startswith("lib:")}
    report = topology.resolve(
        human,
        tool=topology.TOOL_V1,
        cohorts=topology.cohorts_for_species("human"),
        snapshot=pd.read_parquet(inputs / "topology.parquet"),
    )
    classes = {}
    for r in report.to_dict("records"):
        if not r["found"]:
            continue
        topo = r["per_residue_topology"]
        first = next((c for c in topo if c in "IO"), "")
        classes[r["key"]] = classify(
            int(r["tm_helix_count"]),
            int(r["signal_peptide_length"]),
            {"I": "cytoplasmic", "O": "extracellular"}.get(first, "indeterminate"),
            r["key"] in gpi,
        )
    sequences = {name: seq for name, seq in fasta(scored / "ortholog_picks.fasta").items()}
    encodings = {s: defaultdict(set) for s in pool.SPECIES}
    for species in pool.SPECIES:
        for r in pool.read_peptides(inputs / pool.PEPTIDES[species]):
            if r["primary_assembly"]:
                encodings[species][(r["gene"], r["sequence"])].add(r["transcript"])
    existing = defaultdict(list)
    fixed, old = {}, set()
    for r in manifest:
        if r["key"] not in proteins:
            continue
        sequence = proteins[r["key"]]
        digest = hashlib.sha256(sequence.encode()).hexdigest()
        if r["set"] == "ortholog":
            species, gene = r["key"].split(":")[1:3]
            existing[(species, gene, digest)].append(r["candidate_id"])
            if r["library"] == "l1" and ld.is_active(r):
                old.add((species, gene, digest))
        elif r["library"] == "l1" and ld.is_active(r):
            fixed[sequence] = r["candidate_id"]
    eligible, rejected = [], {}
    for row in read(scored / "ortholog_provisional_picks.tsv"):
        if row["sampling_eligible"] != "True":
            continue
        species, gene, acc = row["species"], row["gene"], row["picked_acc"]
        if gene not in classes:
            raise ValueError(f"{gene}: exact human reference lacks a v1 topology class")
        seq = sequences[f"{species}|{gene}|{acc}"]
        native = sorted(encodings[species].get((row["ortholog_ensembl_gene"], seq), ()))
        old_ids = sorted(existing.get((species, gene, row["picked_sequence_sha256"]), ()))
        fid = f"{species}:{gene}:{acc}"
        if not native and not old_ids:
            rejected[fid] = "no exact Ensembl translation or existing designed CDS"
        if seq in fixed:
            rejected[fid] = "protein duplicates active non-ortholog " + fixed[seq]
        eligible.append(
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
                "form_id": fid,
                "human_ensembl_gene": row["human_ensembl_gene"],
                "ortholog_ensembl_gene": row["ortholog_ensembl_gene"],
                "sequence_sha256": row["picked_sequence_sha256"],
                "ensembl_transcripts": ";".join(native),
                "existing_candidate_ids": ";".join(old_ids),
                "encoding_source": "ensembl116" if native else "existing_designed_cds",
                "in_april_ortholog_list": row["in_april_ortholog_list"],
                "sequence": seq,
            }
        )
    for _ in range(20):
        selected, replacements, species_summary = draw(eligible, rejected)
        seen = dict(fixed)
        extra = {}
        for r in selected:
            if r["sequence"] in seen:
                extra[r["form_id"]] = "protein duplicates " + seen[r["sequence"]]
            else:
                seen[r["sequence"]] = r["form_id"]
        if not extra:
            break
        rejected.update(extra)
    else:
        raise ValueError("duplicate replacement did not converge")
    selected.sort(key=lambda r: (r["species"], r["band"], r["human_symbol"]))
    assert len(selected) == 200
    assert len({r["sequence"] for r in selected}) == 200
    new = {(r["species"], r["human_symbol"], r["sequence_sha256"]) for r in selected}
    out.mkdir(parents=True, exist_ok=True)
    pool.write_tsv(
        out / "ortholog_panel_compara116.tsv",
        [{k: v for k, v in r.items() if k != "sequence"} for r in selected],
    )
    pool.write_tsv(
        out / "panel_rejections_compara116.tsv",
        [{"form_id": f, "reason": reason} for f, reason in sorted(rejected.items())],
    )
    if replacements:
        pool.write_tsv(out / "panel_replacements_compara116.tsv", replacements)
    changes = [
        {"change": "removed", "species": s, "gene": g, "sequence_sha256": h}
        for s, g, h in sorted(old - new)
    ]
    changes += [
        {"change": "added", "species": s, "gene": g, "sequence_sha256": h}
        for s, g, h in sorted(new - old)
    ]
    pool.write_tsv(out / "panel_changes_compara116.tsv", changes)
    result = {
        "seed": SEED,
        "per_band": PER_SPECIES,
        "total": len(selected),
        "species": species_summary,
        "retained_proteins": len(old & new),
        "removed_proteins": len(old - new),
        "added_proteins": len(new - old),
        "drawn_from_genes_absent_from_april_list": sum(
            r["in_april_ortholog_list"] == "False" for r in selected
        ),
        "new_forms_without_existing_member_id": sum(
            not r["existing_candidate_ids"] for r in selected
        ),
        "classes": dict(Counter(r["topology_class"] for r in selected)),
        "status": "replacement panel selected; member-ID allocation and construct regeneration are downstream",
    }
    (out / "ortholog_panel_compara116_summary.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--scored", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(run(args.inputs, args.scored, args.out), indent=2))


if __name__ == "__main__":
    main()
