"""Verify the published ColabFold models are the protein their payload row claims.

The payload records ``proteoform_match = ours_by_construction`` for every ColabFold model:
the model is asserted to be our designed CDS translated, not a record fetched by
accession. That assertion is what makes the model usable for tag-site geometry, and it has
never been checked against the deposited coordinates.

This checks it the only way that settles it -- on the sequence. For each model it reads the
chain out of the PDB, translates the member's designed CDS from
``library-candidates/candidate_cds.json``, and requires exact equality. Equal length is not
evidence: 16 L1 members have matched a topology row of identical length but 1-3 different
residues, and 23 of 200 wrong AlphaFold models were the same length as the right one.

It also recomputes mean pLDDT from the C-alpha B-factors, where ColabFold writes it, and
compares that against the ``fold_summary.tsv`` published beside the models -- so a score
quoted downstream is checked against the coordinates it came from rather than trusted.

    uv run --locked python verify_colabfold_models.py --folder 2026-10-14-mismatched
    uv run --locked python verify_colabfold_models.py --folder 2026-10-14-mismatched --limit 5
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import os
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from deliverome_analysis.cloning import library_design as ld
from deliverome_analysis.core.s3 import client
from deliverome_analysis.surfaceome.membrane_side import read_model

HERE = Path(__file__).parent
EXP = "sprint4/s4e9-orf-design-fg-panel"
MODEL_PREFIX = f"{EXP}/predicted-structures/colabfold"
PAYLOAD = HERE / "library_structure_payload.tsv"
MANIFEST = HERE / "candidate_manifest_tranches.tsv"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--folder", required=True, help="date folder under predicted-structures/colabfold/"
    )
    ap.add_argument("--limit", type=int, default=0, help="check only the first N models")
    ap.add_argument("--out", type=Path, default=None, help="write per-model results as JSON")
    args = ap.parse_args(argv)

    c = client()
    bucket = os.environ["S3_BUCKET_PROCESSED"]
    base = f"{MODEL_PREFIX}/{args.folder}"

    # The shared helper owns the candidate_id -> designed-CDS mapping (candidate_cds.json
    # is keyed by design id, not candidate id) and translates it the same way the payload
    # builder does, so this checks the model against the protein of record rather than a
    # second opinion about what that protein is.
    designed = ld.designed_proteins(MANIFEST, library=ld.LIBRARY_L1)

    summary = {
        r["name"]: r
        for r in csv.DictReader(
            io.StringIO(
                c.get_object(Bucket=bucket, Key=f"{base}/fold_summary.tsv")["Body"].read().decode()
            ),
            delimiter="\t",
        )
    }

    payload = {
        r["candidate_id"]: r
        for r in csv.DictReader(PAYLOAD.open(), delimiter="\t")
        # Match the folder as a path segment: a bare substring test makes "2026-10-14"
        # also select the 182 models under "2026-10-14-mismatched/".
        if r["source"] == "colabfold"
        and f"/{args.folder}/" in r["s3_key"]
        and r["is_chosen"] == "1"
    }
    ids = sorted(payload)
    if args.limit:
        ids = ids[: args.limit]
    print(
        f"{len(ids)} chosen models in {args.folder}; designed proteins available: {len(designed)}"
    )

    def check(cid: str) -> dict:
        row = payload[cid]
        out: dict = {"candidate_id": cid, "gene_symbol": row["gene_symbol"]}
        entry = designed.get(cid)
        if not entry:
            return {**out, "verdict": "no_designed_protein"}
        want = entry["sequence"]
        # The payload records the exact key; do not rebuild it. Some folds are named by
        # design id rather than candidate id (the 2026-11-10 ortholog round), so a
        # reconstructed "{candidate_id}.pdb" 404s on those.
        key = row["s3_key"].replace("processed://", "")
        out["model_key"] = key
        stem = Path(key).stem
        with tempfile.TemporaryDirectory() as td:
            local = Path(td) / f"{stem}.pdb"
            c.download_file(bucket, key, str(local))
            model = read_model(local)
        got, plddt = model.sequence, model.plddt
        out |= {
            "model_len": len(got),
            "designed_len": len(want),
            "payload_ordered_len": int(row["ordered_length_aa"] or 0),
            "sequence_identical": got == want,
            "mean_plddt_from_model": round(float(plddt.mean()), 2),
            "mean_plddt_in_summary": (summary.get(stem, {}).get("mean_plddt") or "").strip()
            or None,
            "mean_plddt_in_payload": (row["mean_plddt"] or "").strip() or None,
        }
        if got != want:
            out["first_diff"] = next(
                (i for i, (a, b) in enumerate(zip(got, want, strict=False)) if a != b),
                min(len(got), len(want)),
            )
        out["verdict"] = "ok" if out["sequence_identical"] else "SEQUENCE_MISMATCH"
        return out

    with ThreadPoolExecutor(max_workers=8) as pool:
        res = list(pool.map(check, ids))

    bad = [r for r in res if r["verdict"] != "ok"]
    checked = [r for r in res if "sequence_identical" in r]
    unchecked = [r for r in res if "sequence_identical" not in r]
    print(
        f"\nmodels compared                        : {len(checked)}/{len(res)}"
        + (f"  ({len(unchecked)} had no designed protein)" if unchecked else "")
    )
    print(
        f"sequence identical to the designed CDS : {sum(1 for r in checked if r['sequence_identical'])}/{len(checked)}"
    )
    print(
        f"length agrees with payload ordered_len : {sum(1 for r in checked if r['model_len'] == r['payload_ordered_len'])}/{len(checked)}"
    )
    scored = [r for r in checked if r.get("mean_plddt_in_summary")]
    drift = [
        r
        for r in scored
        if abs(float(r["mean_plddt_in_summary"]) - r["mean_plddt_from_model"]) > 0.5
    ]
    print(f"published mean pLDDT within 0.5 of coords : {len(scored) - len(drift)}/{len(scored)}")
    missing_in_payload = [r for r in scored if not r.get("mean_plddt_in_payload")]
    print(f"score in S3 but absent from the payload   : {len(missing_in_payload)}/{len(scored)}")
    if bad:
        print("\nMISMATCHES:")
        for r in bad[:20]:
            print(
                f"  {r['candidate_id']} {r['gene_symbol']}: {r['verdict']} {r.get('first_diff', '')}"
            )
    if drift:
        print("\npLDDT drift (published vs coordinates):")
        for r in drift[:10]:
            print(
                f"  {r['candidate_id']}: summary {r['mean_plddt_in_summary']} vs model {r['mean_plddt_from_model']}"
            )
    if args.out:
        args.out.write_text(json.dumps(res, indent=1))
        print(f"\nwrote {args.out}")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
