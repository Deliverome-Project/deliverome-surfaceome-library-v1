"""Read-only transcript-completeness evidence for every active L1/L2 member.

NF means an annotated boundary is not confirmed, not proven truncation. ATG and
absence of NF tags do not establish full length. No model is repicked or repaired.
Baseline annotations are pinned115 comparators, not an assertion of original
source release; newly resolved and selected Compara116 forms also
receive separately labelled release116 observations, without transferring evidence
across a CDS mismatch. Construct-only rescue/giant arms remain explicit.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

RULE_VERSION = "library-completeness-v1"
RULES = {
    "COMP001": "Audit all active L1/L2 rows; count retired exclusions explicitly.",
    "COMP002": "Use the documented comparison release and exact transcript identity, never infer original release or use gene-symbol fallback.",
    "COMP003": "Require exact pinned CDS or a checksum-bound recorded stop/boundary exception.",
    "COMP004": "NF tags report unconfirmed start/end; no inference of proven truncation.",
    "COMP005": "Compare complete strand-oriented start/stop triplets at spliced CDS boundaries.",
    "COMP006": "Retain biotype, TSL, phase and gene/assembly mapping evidence without automatic exclusion.",
    "COMP007": "Separately sourced rescue/giant constructed arms are not missing transcript CDSs.",
}
ASSEMBLIES = {"human": "GRCh38", "mouse": "GRCm39", "cyno": "Macaca_fascicularis_6.0"}
FEATURES = {"transcript", "exon", "CDS", "start_codon", "stop_codon"}


def sha(value: str) -> str:
    return hashlib.sha256(value.encode()).hexdigest()


def attrs(text: str) -> dict:
    result = defaultdict(list)
    for key, value in re.findall(r'(\S+) "([^"]*)"', text):
        result[key].append(value)
    return dict(result)


def fasta_selected(path: Path, wanted: set[str]) -> dict:
    """Index exact named transcripts only; retain FASTA assembly/gene/version tags."""
    result, header, sequence = {}, None, []

    def store():
        if header is None:
            return
        transcript = header.split()[0].split(".")[0]
        if transcript not in wanted:
            return
        record = {
            "header": header,
            "transcript_version": header.split()[0],
            "cds": "".join(sequence).upper(),
            "tags": dict(t.split(":", 1) for t in header.split()[1:] if ":" in t),
        }
        if transcript in result and result[transcript] != record:
            raise ValueError(f"Conflicting FASTA transcript identity: {transcript}")
        result[transcript] = record

    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                store()
                header, sequence = line[1:].strip(), []
            elif header is not None:
                sequence.append(line.strip())
    store()
    return result


def gtf_selected(path: Path, wanted: set[str]) -> dict:
    """One streaming pass, O(1) transcript lookup; no large gene-name regex."""
    result = defaultdict(list)
    transcript_pattern = re.compile(r'\btranscript_id "([^"]+)"')
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] not in FEATURES:
                continue
            match = transcript_pattern.search(fields[8])
            if not match or match[1].split(".")[0] not in wanted:
                continue
            result[match[1].split(".")[0]].append(
                {
                    "feature": fields[2],
                    "chromosome": fields[0],
                    "start": int(fields[3]),
                    "end": int(fields[4]),
                    "strand": fields[6],
                    "phase": fields[7],
                    "attributes": attrs(fields[8]),
                }
            )
    return dict(result)


def oriented_positions(features: list[dict], plus: bool) -> list[int]:
    result = []
    for f in sorted(features, key=lambda f: f["start"], reverse=not plus):
        result.extend(
            range(f["start"], f["end"] + 1) if plus else range(f["end"], f["start"] - 1, -1)
        )
    return result


def annotation_evidence(features: list[dict] | None) -> dict:
    if not features:
        return {"annotation_found": False}
    loci = {(f["chromosome"], f["strand"]) for f in features}
    if len(loci) != 1:
        return {"annotation_found": False, "mapping_error": "multiple_chromosome_or_strand_loci"}
    chromosome, strand = next(iter(loci))
    if strand not in {"+", "-"}:
        return {"annotation_found": False, "mapping_error": "unknown_strand"}
    plus = strand == "+"
    by_kind = {kind: [f for f in features if f["feature"] == kind] for kind in FEATURES}
    transcript = by_kind["transcript"][0] if len(by_kind["transcript"]) == 1 else None
    coding = sorted(by_kind["CDS"], key=lambda f: f["start"], reverse=not plus)
    tags = sorted({t for f in features for t in f["attributes"].get("tag", [])})
    merged = defaultdict(set)
    for feature in features:
        for key, values in feature["attributes"].items():
            merged[key].update(values)
    first = oriented_positions(coding[:1], plus)[:3] if coding else []
    # A codon can be split across exon/CDS segments. Include the next segment
    # when the first coding segment contributes fewer than three positions.
    if coding and len(first) < 3:
        first = oriented_positions(coding[:3], plus)[:3]
    last = oriented_positions(coding[-1:], plus)[-3:] if coding else []
    if coding and len(last) < 3:
        last = oriented_positions(coding[-3:], plus)[-3:]
    start_positions = oriented_positions(by_kind["start_codon"], plus)
    stop_positions = oriented_positions(by_kind["stop_codon"], plus)
    exons = sorted(by_kind["exon"], key=lambda f: f["start"], reverse=not plus)

    def transcript_position(position):
        offset = 0
        for exon in exons:
            if exon["start"] <= position <= exon["end"]:
                return offset + (position - exon["start"] if plus else exon["end"] - position)
            offset += exon["end"] - exon["start"] + 1
        return None

    stop_relation = "no_stop_feature"
    if stop_positions:
        if len(stop_positions) == 3 and stop_positions == last:
            stop_relation = "overlaps_last_three_CDS_bases"
        elif coding and exons:
            end_position = transcript_position(last[-1]) if last else None
            mapped = [transcript_position(position) for position in stop_positions]
            stop_relation = (
                "immediately_after_CDS_in_spliced_transcript"
                if end_position is not None and mapped == [end_position + i for i in (1, 2, 3)]
                else "not_at_CDS_boundary"
            )
        else:
            stop_relation = "unassessable_without_exons_or_CDS"
    start_relation = (
        "first_three_CDS_bases"
        if len(first) == 3 and start_positions == first
        else "not_at_CDS_boundary"
        if start_positions
        else "no_start_feature"
    )

    def attr(key):
        return ";".join(sorted(merged.get(key, set())))

    five_prime = (transcript["start"] if plus else transcript["end"]) if transcript else None
    return {
        "annotation_found": True,
        "transcript_feature_found": transcript is not None,
        "chromosome": chromosome,
        "strand": strand,
        "gene_id": attr("gene_id"),
        "gene_name": attr("gene_name"),
        "gene_biotype": attr("gene_biotype"),
        "transcript_biotype": attr("transcript_biotype"),
        "transcript_version": attr("transcript_version"),
        "transcript_support_level": attr("transcript_support_level"),
        "tags": ";".join(tags),
        "cds_start_NF": "cds_start_NF" in tags,
        "cds_end_NF": "cds_end_NF" in tags,
        "mRNA_start_NF": "mRNA_start_NF" in tags,
        "mRNA_end_NF": "mRNA_end_NF" in tags,
        "first_CDS_phase": coding[0]["phase"] if coding else "",
        "CDS_feature_nt": sum(f["end"] - f["start"] + 1 for f in coding),
        "first_CDS_positions": ";".join(map(str, first)),
        "last_CDS_positions": ";".join(map(str, last)),
        "annotated_start_positions": ";".join(map(str, start_positions)),
        "annotated_stop_positions": ";".join(map(str, stop_positions)),
        "start_feature_relation": start_relation,
        "stop_feature_relation": stop_relation,
        "cds_begins_at_transcript_boundary": bool(first and five_prime == first[0]),
        "exon_segments": ";".join(f"{f['start']}-{f['end']}" for f in exons),
    }


def match_cds(entry: dict, pinned: dict | None) -> dict:
    source = entry.get("cds", "")
    boundary = entry.get("boundary_design", {})
    exception = ""
    if boundary:
        original = boundary.get("original_cds", "")
        if sha(source) != boundary.get("designed_cds_sha256") or sha(original) != boundary.get(
            "original_cds_sha256"
        ):
            return {"cds_match": "invalid_engineering_checksum", "known_exception": ""}
        source = original
        exception = "checksum_bound_boundary_design:" + str(boundary.get("action", ""))
    if pinned is None:
        return {"cds_match": "transcript_missing_from_pinned_CDS", "known_exception": exception}
    target = pinned["cds"]
    stop_added = entry.get("stop_added", "")
    if source == target:
        match = "exact_pinned_CDS"
    elif stop_added in {"TAA", "TAG", "TGA"} and source == target + stop_added:
        match = "exact_pinned_CDS_plus_recorded_stop"
        exception = ";".join(filter(None, [exception, "recorded_stop_added:" + stop_added]))
    else:
        match = "differs_from_pinned_CDS"
    return {
        "cds_match": match,
        "known_exception": exception,
        "pinned_cds_sha256": sha(target),
        "pinned_cds_nt": len(target),
        "comparison_cds_sha256": sha(source),
        "length_difference_nt": len(source) - len(target),
        "nucleotide_differences_if_equal_length": sum(
            a != b for a, b in zip(source, target, strict=True)
        )
        if len(source) == len(target)
        else "",
    }


def audit_row(
    row: dict,
    entry: dict | None,
    construct: dict | None,
    pinned: dict | None,
    features: list[dict] | None,
    *,
    snapshot: str,
    release: int,
    role: str,
    priority: int,
) -> dict:
    entry, construct = entry or {}, construct or {}
    raw = entry.get("cds", "")
    separate = (
        not raw
        and bool(construct.get("construct"))
        and (
            (
                row["set"] == "expression_rescue"
                and construct.get("source", "").startswith("rescue_arm")
            )
            or (
                row["set"] == "giant_truncation"
                and construct.get("source", "").startswith("giant_truncation")
            )
        )
    )
    annotation = annotation_evidence(features)
    match = (
        match_cds(entry, pinned)
        if raw
        else {"cds_match": "construct_only_arm" if separate else "missing_raw_CDS"}
    )
    flags = []
    if not separate:
        if not raw:
            flags.append("missing_raw_CDS")
        if not str(entry.get("transcript", "")):
            flags.append("missing_transcript_identity")
        if not annotation.get("annotation_found"):
            flags.append("unknown_transcript_annotation")
        if match["cds_match"] not in {"exact_pinned_CDS", "exact_pinned_CDS_plus_recorded_stop"}:
            flags.append(match["cds_match"])
        if annotation.get("annotation_found"):
            if not annotation["transcript_feature_found"]:
                flags.append("missing_transcript_GTF_feature")
            for tag in ("cds_start_NF", "cds_end_NF", "mRNA_start_NF", "mRNA_end_NF"):
                if annotation[tag]:
                    flags.append(tag + "_not_confirmed")
            if annotation["start_feature_relation"] != "first_three_CDS_bases":
                flags.append("start_feature:" + annotation["start_feature_relation"])
            if annotation["stop_feature_relation"] not in {
                "overlaps_last_three_CDS_bases",
                "immediately_after_CDS_in_spliced_transcript",
            }:
                flags.append("stop_feature:" + annotation["stop_feature_relation"])
            if annotation["first_CDS_phase"] != "0":
                flags.append("first_CDS_phase:" + annotation["first_CDS_phase"])
            if annotation["transcript_biotype"] not in {"", "protein_coding"}:
                flags.append("transcript_biotype_review:" + annotation["transcript_biotype"])
    fasta_tags = pinned.get("tags", {}) if pinned else {}
    fasta_gene = fasta_tags.get("gene", "").split(".")[0]
    gtf_gene = annotation.get("gene_id", "")
    location = next(
        (
            fasta_tags[name]
            for name in ("chromosome", "primary_assembly", "scaffold")
            if name in fasta_tags
        ),
        "",
    )
    location_parts = location.split(":")
    assembly_match = bool(location_parts and location_parts[0] == ASSEMBLIES.get(row["species"]))
    gene_match = bool(fasta_gene and fasta_gene == gtf_gene)
    chromosome_match = bool(
        len(location_parts) > 1 and location_parts[1] == annotation.get("chromosome")
    )
    fasta_strand = {"1": "+", "-1": "-"}.get(
        location_parts[4] if len(location_parts) > 4 else "", ""
    )
    strand_match = bool(fasta_strand and fasta_strand == annotation.get("strand"))
    fasta_version = (
        pinned.get("transcript_version", "").rsplit(".", 1)[1]
        if pinned and "." in pinned.get("transcript_version", "")
        else ""
    )
    gtf_version = annotation.get("transcript_version", "")
    version_match = bool(fasta_version and gtf_version and fasta_version == gtf_version)
    length_relation, pinned_sense, gtf_sense = "unassessable_missing_annotation_or_CDS", "", ""
    if pinned and annotation.get("annotation_found"):
        target = pinned["cds"]
        pinned_sense = len(target) - (3 if target[-3:] in {"TAA", "TAG", "TGA"} else 0)
        stop_relation = annotation["stop_feature_relation"]
        if stop_relation in {
            "overlaps_last_three_CDS_bases",
            "immediately_after_CDS_in_spliced_transcript",
        }:
            gtf_sense = annotation["CDS_feature_nt"] - (
                3 if stop_relation == "overlaps_last_three_CDS_bases" else 0
            )
            length_relation = (
                "sense_coding_lengths_match"
                if pinned_sense == gtf_sense
                else "sense_coding_lengths_differ"
            )
        else:
            length_relation = "unassessable_stop_annotation"
    if not separate and raw and not raw.startswith("ATG"):
        flags.append("non_ATG_initiation_review")
    if not separate and pinned and annotation.get("annotation_found"):
        if not gene_match:
            flags.append("stable_gene_mapping_mismatch")
        if not assembly_match or not chromosome_match:
            flags.append("assembly_or_chromosome_mapping_unconfirmed")
        if not strand_match:
            flags.append("FASTA_GTF_strand_mismatch_or_unconfirmed")
        if not version_match:
            flags.append("FASTA_GTF_transcript_version_mismatch_or_unconfirmed")
        if length_relation != "sense_coding_lengths_match":
            flags.append("FASTA_GTF_CDS_length:" + length_relation)
    start_status = (
        "ATG"
        if raw.startswith("ATG")
        else "non_ATG"
        if raw
        else "construct_only_ATG"
        if separate and construct["construct"].startswith("ATG")
        else "raw_CDS_unavailable"
    )
    return {
        "snapshot": snapshot,
        "priority": priority,
        "candidate_id": row["candidate_id"],
        "library": row["library"],
        "candidate_set": row["set"],
        "species": row["species"],
        "member_key": row["key"],
        "annotation_release": release,
        "annotation_role": role,
        "source_transcript": entry.get("transcript", ""),
        "source_origin": entry.get("origin", ""),
        "source_route": entry.get("route", ""),
        "source_cds_sha256": sha(raw) if raw else "",
        "source_start_codon": raw[:3],
        "source_last_triplet": raw[-3:],
        "start_status": start_status,
        "source_stop_added": entry.get("stop_added", ""),
        "construct_source": construct.get("source", ""),
        "separate_construct_only_arm": separate,
        "rule_version": RULE_VERSION,
        "annotation_status": "construct_only_not_transcript_CDS"
        if separate
        else "annotation_review"
        if flags
        else "boundary_features_supported",
        "review_flags": ";".join(sorted(set(flags))),
        "full_length_proven_by_audit": False,
        "fasta_transcript_version": pinned.get("transcript_version", "") if pinned else "",
        "fasta_gene": fasta_gene,
        "fasta_assembly_location": location,
        "stable_gene_match": gene_match,
        "assembly_match": assembly_match,
        "chromosome_match": chromosome_match,
        "fasta_strand": fasta_strand,
        "strand_match": strand_match,
        "fasta_transcript_version_number": fasta_version,
        "transcript_version_match": version_match,
        "pinned_sense_coding_nt": pinned_sense,
        "gtf_sense_coding_nt": gtf_sense,
        "FASTA_GTF_CDS_length_relation": length_relation,
        **match,
        **annotation,
    }


def table(path: Path) -> list[dict]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write(path: Path, rows: list[dict]) -> None:
    columns = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--snapshot",
        nargs=4,
        action="append",
        required=True,
        metavar=("LABEL", "MANIFEST", "CDS", "CONSTRUCTS"),
    )
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--baseline-cds", type=Path, required=True)
    parser.add_argument("--member-map", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    specification = json.loads(args.inputs.read_text())
    inputs = {}
    for name, item in specification["files"].items():
        path = args.input_dir / name
        if hashlib.sha256(path.read_bytes()).hexdigest() != item["sha256"]:
            raise ValueError(f"Pinned input checksum mismatch: {name}")
        inputs[(item["release"], item["species"], item["kind"])] = path
    baseline = json.loads(args.baseline_cds.read_text())
    selected = {r["member_key"] for r in table(args.member_map)}
    snapshots, observations, wanted, provenance = {}, [], defaultdict(set), {}
    for label, manifest, cds, constructs in args.snapshot:
        if label in snapshots:
            raise ValueError("Duplicate snapshot label")
        manifest_rows, native, built = (
            table(Path(manifest)),
            json.loads(Path(cds).read_text()),
            json.loads(Path(constructs).read_text()),
        )
        active = [
            r
            for r in manifest_rows
            if r["library"] in {"l1", "l2"} and not r.get("retired_reason", "").strip()
        ]
        if len({r["candidate_id"] for r in active}) != len(active):
            raise ValueError("Duplicate active member identity")
        snapshots[label] = {
            "manifest_rows": len(manifest_rows),
            "active_rows": len(active),
            "retired_rows_excluded": sum(
                bool(r.get("retired_reason", "").strip()) for r in manifest_rows
            ),
        }
        provenance[label] = {
            name: {
                "file": Path(path).name,
                "sha256": hashlib.sha256(Path(path).read_bytes()).hexdigest(),
            }
            for name, path in (("manifest", manifest), ("cds", cds), ("constructs", constructs))
        }
        for row in active:
            key, entry = row["key"], native.get(row["key"], {})
            release = 115 if key in baseline or not entry else 116 if key in selected else 0
            priority = (
                0
                if row["library"] == "l1" and row["set"] == "library_reference"
                else 1
                if row["library"] == "l1"
                and key in selected
                and not entry.get("cds", "").startswith("ATG")
                else 2
                if row["library"] == "l1"
                else 3
            )
            roles = [
                (
                    release,
                    "pinned_baseline_comparison_release"
                    if release == 115
                    else "documented_Compara116_resolution_release"
                    if release == 116
                    else "unknown_comparison_release",
                )
            ]
            if key in selected and release != 116:
                roles.append((116, "supplementary_Compara116_same_transcript"))
            for annotation_release, role in roles:
                observations.append(
                    (label, row, entry, built.get(key), annotation_release, role, priority)
                )
                transcript = entry.get("transcript", "").split(".")[0]
                if transcript and annotation_release:
                    wanted[(annotation_release, row["species"])].add(transcript)
    indexes = {}
    for key, transcripts in sorted(wanted.items()):
        release, species = key
        cds_path, gtf_path = (
            inputs.get((release, species, "cds")),
            inputs.get((release, species, "gtf")),
        )
        if cds_path is None or gtf_path is None:
            raise ValueError(f"Annotation inputs absent for release/species: {key}")
        indexes[key] = (fasta_selected(cds_path, transcripts), gtf_selected(gtf_path, transcripts))
    rows = []
    for label, row, entry, construct, release, role, priority in observations:
        transcript = entry.get("transcript", "").split(".")[0]
        fastas, gtfs = indexes.get((release, row["species"]), ({}, {}))
        rows.append(
            audit_row(
                row,
                entry,
                construct,
                fastas.get(transcript),
                gtfs.get(transcript),
                snapshot=label,
                release=release,
                role=role,
                priority=priority,
            )
        )
    rows.sort(
        key=lambda r: (r["priority"], r["snapshot"], r["candidate_id"], r["annotation_release"])
    )
    grouped_rows = []
    for label, summary in snapshots.items():
        primary = [
            r
            for r in rows
            if r["snapshot"] == label
            and r["annotation_role"] != "supplementary_Compara116_same_transcript"
        ]
        if len(primary) != summary["active_rows"]:
            raise ValueError("An active member disappeared from the primary audit")
        summary["counts_by_library_set_species_start"] = dict(
            sorted(
                Counter(
                    "|".join(
                        str(r[k]) for k in ("library", "candidate_set", "species", "start_status")
                    )
                    for r in primary
                ).items()
            )
        )
        summary["primary_review_flag_counts"] = dict(
            sorted(
                Counter(
                    flag for r in primary for flag in r["review_flags"].split(";") if flag
                ).items()
            )
        )
        summary["primary_annotation_status"] = dict(
            sorted(Counter(r["annotation_status"] for r in primary).items())
        )
        groups = defaultdict(list)
        for row in primary:
            groups[
                tuple(row[k] for k in ("library", "candidate_set", "species", "start_status"))
            ].append(row)
        for group, members in sorted(groups.items()):
            counts = {
                "snapshot": label,
                **dict(
                    zip(("library", "candidate_set", "species", "start_status"), group, strict=True)
                ),
                "members": len(members),
                "members_with_review_flags": sum(bool(r["review_flags"]) for r in members),
            }
            for flag in ("cds_start_NF", "cds_end_NF", "mRNA_start_NF", "mRNA_end_NF"):
                counts[flag + "_members"] = sum(r.get(flag) is True for r in members)
            counts["unknown_annotation_members"] = sum(
                not r.get("annotation_found")
                for r in members
                if not r["separate_construct_only_arm"]
            )
            counts["exact_pinned_CDS_or_recorded_stop_members"] = sum(
                r["cds_match"] in {"exact_pinned_CDS", "exact_pinned_CDS_plus_recorded_stop"}
                for r in members
            )
            counts["flag_counts_json"] = json.dumps(
                dict(
                    sorted(
                        Counter(
                            f for r in members for f in r["review_flags"].split(";") if f
                        ).items()
                    )
                ),
                sort_keys=True,
            )
            grouped_rows.append(counts)
    args.out.mkdir(parents=True, exist_ok=False)
    ledger = args.out / "library_completeness.tsv"
    write(ledger, rows)
    write(args.out / "library_completeness_groups.tsv", grouped_rows)
    (args.out / "library_completeness.tsv.gz").write_bytes(
        gzip.compress(ledger.read_bytes(), mtime=0)
    )
    write(args.out / "library_completeness_review.tsv", [r for r in rows if r["review_flags"]])
    commit = subprocess.run(
        ["git", "-C", str(Path(__file__).parent), "rev-parse", "HEAD"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    summary = {
        "rule_version": RULE_VERSION,
        "rules": RULES,
        "snapshots": snapshots,
        "annotation_inputs": specification,
        "snapshot_provenance": provenance,
        "script_commit": commit,
        "script_has_uncommitted_changes": bool(
            subprocess.run(
                [
                    "git",
                    "-C",
                    str(Path(__file__).parent),
                    "status",
                    "--porcelain",
                    "--",
                    Path(__file__).name,
                ],
                check=True,
                capture_output=True,
                text=True,
            ).stdout.strip()
        ),
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "command_argv": [Path(__file__).name, *sys.argv[1:]],
        "baseline_cds_sha256": hashlib.sha256(args.baseline_cds.read_bytes()).hexdigest(),
        "member_map_sha256": hashlib.sha256(args.member_map.read_bytes()).hexdigest(),
        "sequences_changed": 0,
        "automatic_model_replacements": 0,
        "interpretation": "NF means a boundary is not confirmed; neither ATG nor absence of NF proves full length.",
    }
    summary["output_sha256"] = {
        p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(args.out.iterdir())
    }
    (args.out / "library_completeness_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps({"snapshots": snapshots, "rows": len(rows)}, indent=2))


if __name__ == "__main__":
    main()
