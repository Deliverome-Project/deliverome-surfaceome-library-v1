"""Append the sequence-verified ortholog predictions to a new topology snapshot.

The old snapshot and every prior row are preserved. The same JSONL is loaded
into both D1 databases by accessible-surfaceome's topology importer.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

from deliverome_analysis.core.s3 import external

TOOL_VERSIONS = ("deeptmhmm-1.0.24", "deeptmhmm-1.0-academic")


def append_records(base: pd.DataFrame, records: list[dict], tool_version: str) -> pd.DataFrame:
    """Append already derived importer rows without relabelling historical predictions."""
    if tool_version not in TOOL_VERSIONS:
        raise ValueError("unsupported DeepTMHMM tool version")
    if not records:
        raise ValueError("no topology records to append")
    rows = []
    keys = set()
    for record in records:
        sequence = record["sequence"]
        if record["tool_version"] != tool_version:
            raise ValueError("record tool version differs from the explicitly selected version")
        if not sequence or sequence != sequence.strip().upper() or "*" in sequence:
            raise ValueError("sequence must be normalized, nonempty and stop-free")
        if len(sequence) != len(record["per_residue_topology"]):
            raise ValueError("topology differs from sequence length")
        if set(record["per_residue_topology"]) - set("SIOMB"):
            raise ValueError("invalid DeepTMHMM v1 topology alphabet")
        if record.get("protein_length") != len(sequence):
            raise ValueError("protein_length differs from sequence length")
        row = {col: record.get(col) for col in base.columns}
        row.update(
            {
                "uniprot_acc": record["uniprot_accession"],
                "uniprot_acc_full": record["uniprot_accession_full"],
                "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
                "sequence_length": len(sequence),
            }
        )
        key = (row["topology_version"], row["cohort"], row["uniprot_acc_full"])
        if not all(isinstance(part, str) and part.strip() for part in key):
            raise ValueError("missing topology natural key")
        if key in keys:
            raise ValueError("duplicate topology record")
        keys.add(key)
        clash = base[(base["topology_version"] == key[0]) & (base["cohort"] == key[1])]
        if len(clash):
            raise ValueError("the new version/cohort already exists in the base snapshot")
        rows.append(row)
    result = pd.concat(
        [
            base,
            pd.DataFrame(rows, columns=base.columns).astype(base.dtypes.to_dict()),
        ],
        ignore_index=True,
    )
    assert len(result) == len(base) + len(records)
    assert result.iloc[: len(base)].reset_index(drop=True).equals(base.reset_index(drop=True))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--jsonl", type=Path, required=True)
    parser.add_argument("--base", required=True)
    parser.add_argument("--key", required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--tool-version", choices=TOOL_VERSIONS, default=TOOL_VERSIONS[0])
    parser.add_argument("--base-sha256")
    parser.add_argument("--expected-base-rows", type=int)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    if args.base == args.key:
        raise ValueError("the new snapshot key must differ from the frozen base")
    if args.tool_version == "deeptmhmm-1.0-academic" and (
        not args.base_sha256 or args.expected_base_rows is None
    ):
        raise ValueError("academic model append requires explicit base hash and row count")
    base_path = Path(external.download(args.base))
    base_sha256 = hashlib.sha256(base_path.read_bytes()).hexdigest()
    if args.base_sha256 and base_sha256 != args.base_sha256:
        raise ValueError("base snapshot SHA256 mismatch")
    base = pd.read_parquet(base_path)
    if args.expected_base_rows is not None and len(base) != args.expected_base_rows:
        raise ValueError("base snapshot row count mismatch")
    records = [json.loads(line) for line in args.jsonl.read_text().splitlines() if line.strip()]
    result = append_records(base, records, args.tool_version)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    result.to_parquet(args.out, index=False)
    meta = {
        "base_key": args.base,
        "base_sha256": base_sha256,
        "new_key": args.key,
        "tool_version": args.tool_version,
        "rows_before": len(base),
        "rows_added": len(records),
        "rows_after": len(result),
        "source_jsonl_sha256": hashlib.sha256(args.jsonl.read_bytes()).hexdigest(),
        "snapshot_sha256": hashlib.sha256(args.out.read_bytes()).hexdigest(),
    }
    args.out.with_suffix(".json").write_text(json.dumps(meta, indent=2) + "\n")
    print(json.dumps(meta, indent=2))
    if args.write:
        if external.exists(args.key):
            raise ValueError("the new external snapshot key already exists; never overwrite it")
        external.upload(str(args.out), args.key)
        print("published external://" + args.key)


if __name__ == "__main__":
    main()
