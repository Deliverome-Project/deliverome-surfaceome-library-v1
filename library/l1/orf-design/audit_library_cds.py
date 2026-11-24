"""Audit every active L1/L2 raw CDS and domesticated ORF, without changing either.

Expression-rescue and giant-truncation arms have their own constructed ORF record;
absence from the raw transcript-CDS table is recorded separately from a missing ORF.
All other missing inputs, start/stop/frame/ambiguity defects and protein differences
remain explicit. This is a sequence-integrity audit, not an expression experiment.

Each --snapshot argument is LABEL MANIFEST_TSV CDS_JSON CONSTRUCTS_JSON. Multiple
snapshots give a per-member transition table as well as the complete evidence ledger.
Output goes to a fresh directory; no canonical records or sequences are modified.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import subprocess
import sys
from collections import Counter
from pathlib import Path

from Bio.Seq import Seq

RULE_VERSION = "library-cds-audit-v1"
RULES = {
    "CDS001": "Every active member has a domesticated constructed ORF.",
    "CDS002": "Raw transcript CDS exists unless a separately sourced rescue/giant arm applies.",
    "CDS101": "First triplet is ATG; non-ATG is flagged for initiation review.",
    "CDS102": "Nucleotide length is divisible by three.",
    "CDS103": "Every nucleotide is A, C, G or T (case-insensitive).",
    "CDS104": "Final in-frame codon is TAA, TAG or TGA.",
    "CDS105": "No standard-code stop occurs before the final in-frame codon.",
    "CDS106": "Complete diagnostic codons translate with standard nuclear table 1.",
    "CDS201": "Raw CDS and domesticated ORF encode identical proteins, excluding one terminal stop.",
    "CDS301": "Rescue/giant arm source explicitly supports constructed-ORF-only classification.",
}
ISSUE_RULE = {
    "non_atg_start": "CDS101",
    "length_not_multiple_of_3": "CDS102",
    "ambiguous_or_invalid_bases": "CDS103",
    "missing_terminal_stop": "CDS104",
    "internal_stop": "CDS105",
    "translation_error": "CDS106",
}
STOPS = {"TAA", "TAG", "TGA"}
CONSTRUCT_ONLY_SOURCES = {
    "expression_rescue": "rescue_arm",
    "giant_truncation": "giant_truncation",
}


def digest(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


def inspect_sequence(value: str | None) -> dict:
    """Use the standard nuclear code; retain uncertainty instead of normalizing it."""
    if not value:
        return {
            "available": False,
            "length_nt": 0,
            "sequence_sha256": "",
            "start_codon": "",
            "last_triplet": "",
            "length_mod3": "",
            "ambiguous_positions_1based": "",
            "internal_stop_positions_1based": "",
            "protein_sha256": "",
            "protein_length": "",
            "issues": ["missing_sequence"],
        }
    if not isinstance(value, str):
        raise ValueError("CDS/construct sequence must be a string")
    dna = value.upper()
    ambiguous = [str(i + 1) for i, base in enumerate(dna) if base not in "ACGT"]
    issues = []
    if dna[:3] != "ATG":
        issues.append("non_atg_start")
    if len(dna) % 3:
        issues.append("length_not_multiple_of_3")
    if ambiguous:
        issues.append("ambiguous_or_invalid_bases")
    terminal_stop = len(dna) % 3 == 0 and dna[-3:] in STOPS
    if not terminal_stop:
        issues.append("missing_terminal_stop")
    # Translate complete codons only for diagnostics, while the frame defect stays
    # blocking. Ambiguous sequences never receive a definitive protein hash.
    complete = dna[: len(dna) // 3 * 3]
    try:
        protein = str(Seq(complete).translate(table=1)) if complete else ""
    except Exception as exc:
        protein = ""
        issues.append("translation_error:" + type(exc).__name__)
    internal = [
        str(i + 1)
        for i, aa in enumerate(protein)
        if aa == "*" and not (terminal_stop and i == len(protein) - 1)
    ]
    if internal:
        issues.append("internal_stop")
    encoded = protein[:-1] if terminal_stop else protein
    definitive = not ambiguous and len(dna) % 3 == 0 and bool(encoded) and "*" not in encoded
    return {
        "available": True,
        "length_nt": len(dna),
        "sequence_sha256": digest(value),
        "start_codon": dna[:3],
        "last_triplet": dna[-3:],
        "length_mod3": len(dna) % 3,
        "ambiguous_positions_1based": ";".join(ambiguous),
        "internal_stop_positions_1based": ";".join(internal),
        "protein_sha256": digest(encoded) if definitive else "",
        "protein_length": len(encoded) if protein else "",
        "issues": issues,
    }


def audit_member(row: dict, raw: dict | None, constructed: dict | None, snapshot: str) -> dict:
    raw, constructed = raw or {}, constructed or {}
    raw_audit = inspect_sequence(raw.get("cds"))
    construct_audit = inspect_sequence(constructed.get("construct"))
    expected_source = CONSTRUCT_ONLY_SOURCES.get(row["set"])
    construct_only = (
        not raw_audit["available"]
        and construct_audit["available"]
        and expected_source is not None
        and str(constructed.get("source", "")).startswith(expected_source)
    )
    if raw_audit["available"]:
        evidence = "raw_cds_and_construct" if construct_audit["available"] else "raw_cds_only"
    elif construct_only:
        evidence = "construct_only_" + row["set"]
    elif construct_audit["available"]:
        evidence = "construct_only_missing_raw_cds"
    else:
        evidence = "missing_both_sequences"
    issues = [
        "raw:" + issue
        for issue in raw_audit["issues"]
        if not (construct_only and issue == "missing_sequence")
    ]
    issues += ["construct:" + issue for issue in construct_audit["issues"]]
    raw_hash, constructed_hash = raw_audit["protein_sha256"], construct_audit["protein_sha256"]
    if raw_hash and constructed_hash:
        protein_match = "identical" if raw_hash == constructed_hash else "different"
        if protein_match == "different":
            issues.append("cds_vs_construct_protein_mismatch")
    elif construct_only:
        protein_match = "not_applicable_separate_construct_arm"
    else:
        protein_match = "unassessable_missing_or_invalid_sequence"
    rule_hits = []
    for issue in issues:
        if issue == "cds_vs_construct_protein_mismatch":
            rule_hits.append("CDS201")
        elif issue == "raw:missing_sequence":
            rule_hits.append("CDS002")
        elif issue == "construct:missing_sequence":
            rule_hits.append("CDS001")
        else:
            rule_hits.append(ISSUE_RULE[issue.split(":")[1]])

    def engineering_keys(record):
        return {
            k: v
            for k, v in record.items()
            if any(term in k.lower() for term in ("start", "stop", "added", "engineer"))
        }

    result = {
        "snapshot": snapshot,
        "candidate_id": row["candidate_id"],
        "library": row["library"],
        "candidate_set": row["set"],
        "species": row.get("species", ""),
        "member_key": row["key"],
        "raw_transcript": raw.get("transcript", ""),
        "manifest_transcript": row.get("transcript", ""),
        "construct_transcript": constructed.get("transcript", ""),
        "raw_origin": raw.get("origin", ""),
        "raw_route": raw.get("route", ""),
        "construct_source": constructed.get("source", ""),
        "construct_route": constructed.get("route", ""),
        "evidence_class": evidence,
        "cds_vs_construct_protein": protein_match,
        "status": "flagged" if issues else "sequence_checks_pass",
        "rule_version": RULE_VERSION,
        "flagged_rule_ids": ";".join(sorted(set(rule_hits))),
        "construct_only_rule_id": "CDS301" if construct_only else "",
        "raw_engineering_metadata": json.dumps(engineering_keys(raw), sort_keys=True),
        "construct_engineering_metadata": json.dumps(engineering_keys(constructed), sort_keys=True),
        "construct_edits_sha256": digest(json.dumps(constructed.get("edits", []), sort_keys=True)),
        "issues": ";".join(issues),
    }
    for label, details in (("raw", raw_audit), ("construct", construct_audit)):
        result.update(
            {label + "_" + key: value for key, value in details.items() if key != "issues"}
        )
        result[label + "_issues"] = ";".join(details["issues"])
    return result


def audit_snapshot(
    label: str, manifest: list[dict], cds: dict, constructs: dict
) -> tuple[list, dict]:
    active = [
        r
        for r in manifest
        if r["library"] in {"l1", "l2"} and not r.get("retired_reason", "").strip()
    ]
    ids = [r["candidate_id"] for r in active]
    if len(ids) != len(set(ids)):
        raise ValueError("Active member IDs are not unique")
    rows = [audit_member(r, cds.get(r["key"]), constructs.get(r["key"]), label) for r in active]
    summary = {
        "manifest_rows": len(manifest),
        "active_members": len(rows),
        "retired_rows_excluded": sum(bool(r.get("retired_reason", "").strip()) for r in manifest),
        "non_l1_l2_rows_excluded": sum(r["library"] not in {"l1", "l2"} for r in manifest),
        "by_library": {},
    }
    for library in ("l1", "l2"):
        members = [r for r in rows if r["library"] == library]
        summary["by_library"][library] = {
            "active_members": len(members),
            "flagged_members": sum(r["status"] == "flagged" for r in members),
            "evidence_classes": dict(sorted(Counter(r["evidence_class"] for r in members).items())),
            "issue_member_counts": dict(
                sorted(
                    Counter(
                        issue for r in members for issue in r["issues"].split(";") if issue
                    ).items()
                )
            ),
        }
    return rows, summary


def transitions(before: list[dict], after: list[dict]) -> list[dict]:
    old = {r["candidate_id"]: r for r in before}
    new = {r["candidate_id"]: r for r in after}
    result = []
    for cid in sorted(old.keys() | new.keys()):
        a, b = old.get(cid, {}), new.get(cid, {})
        result.append(
            {
                "candidate_id": cid,
                "library": (a or b)["library"],
                "member_key_before": a.get("member_key", ""),
                "member_key_after": b.get("member_key", ""),
                "active_before": bool(a),
                "active_after": bool(b),
                "issues_before": a.get("issues", ""),
                "issues_after": b.get("issues", ""),
                "sequence_changed": bool(
                    a
                    and b
                    and (
                        a["raw_sequence_sha256"] != b["raw_sequence_sha256"]
                        or a["construct_sequence_sha256"] != b["construct_sequence_sha256"]
                    )
                ),
            }
        )
    return result


def write(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--snapshot",
        nargs=4,
        metavar=("LABEL", "MANIFEST", "CDS", "CONSTRUCTS"),
        action="append",
        required=True,
    )
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    labels = [values[0] for values in args.snapshot]
    if len(labels) != len(set(labels)):
        raise ValueError("Snapshot labels must be unique")
    all_rows, snapshots, provenance = [], {}, {}
    by_snapshot = []
    for label, manifest, cds, constructs in args.snapshot:
        with Path(manifest).open() as fh:
            manifest_rows = list(csv.DictReader(fh, delimiter="\t"))
        rows, summary = audit_snapshot(
            label,
            manifest_rows,
            json.loads(Path(cds).read_text()),
            json.loads(Path(constructs).read_text()),
        )
        all_rows.extend(rows)
        by_snapshot.append(rows)
        snapshots[label] = summary
        provenance[label] = {
            name: {
                "file": Path(path).name,
                "sha256": hashlib.sha256(Path(path).read_bytes()).hexdigest(),
            }
            for name, path in (("manifest", manifest), ("cds", cds), ("constructs", constructs))
        }
    args.out.mkdir(parents=True, exist_ok=False)
    write(args.out / "library_cds_audit.tsv", all_rows)
    ledger = args.out / "library_cds_audit.tsv"
    (args.out / "library_cds_audit.tsv.gz").write_bytes(gzip.compress(ledger.read_bytes(), mtime=0))
    flagged = [r for r in all_rows if r["status"] == "flagged"]
    if flagged:
        write(args.out / "library_cds_issues.tsv", flagged)
    if len(by_snapshot) == 2:
        write(args.out / "library_cds_transitions.tsv", transitions(*by_snapshot))
    commit = subprocess.run(
        ["git", "-C", str(Path(__file__).parent), "rev-parse", "HEAD"],
        capture_output=True,
        text=True,
        check=True,
    ).stdout.strip()
    script_status = subprocess.run(
        [
            "git",
            "-C",
            str(Path(__file__).parent),
            "status",
            "--porcelain",
            "--",
            Path(__file__).name,
        ],
        capture_output=True,
        text=True,
        check=True,
    ).stdout.strip()
    summary = {
        "snapshots": snapshots,
        "rule_version": RULE_VERSION,
        "rule_definitions": RULES,
        "script_commit": commit,
        "script_has_uncommitted_changes": bool(script_status),
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "command_argv": [Path(__file__).name, *sys.argv[1:]],
        "input_provenance": provenance,
        "translation_table": 1,
        "sequences_changed": False,
        "scope": "all active L1/L2 records; raw transcript CDS and domesticated pre-tag ORFs",
        "limitations": [
            "Sequence checks do not establish initiation or expression in cells.",
            "Separate rescue/giant constructed arms lack a directly comparable raw transcript CDS.",
        ],
    }
    summary["output_sha256"] = {
        p.name: hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted(args.out.iterdir())
        if p.suffix != ".json"
    }
    (args.out / "library_cds_audit_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
