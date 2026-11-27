"""Publish untagged native folds and append them to the native-model index and structure payload.

A member whose designed protein has no sequence-identical model (no AlphaFold DB entry, no
earlier ColabFold fold) is folded untagged by ``scripts/fold_tagged_modal.py`` from a queue of
``kind == native`` rows. This publishes those folds under a fresh, write-once processed prefix,
appends each to ``tag_refold_native_models.tsv`` (source ``colabfold``), and appends one chosen
row per member to ``library_structure_payload.tsv`` for every member named, whose model is then
found BY SEQUENCE in the native-model index (an AlphaFold DB model identical to the designed
protein first, else our ColabFold fold). Every existing row of both tables is kept byte for
byte (asserted); a member's earlier payload rows describe an earlier protein and stay as the
trail, so a reader takes the last chosen row of a candidate (the 2026-11-22 rows here).

Written 2026-11-22 for the round-2 additions, the isoform redraw and refills, the ortholog
restorations and the GTF MANE Select switches, by Claude Code session l1-review-export-scope
(subagent), under Becca's 2026-11-21 "Yes, rebuild now" and 2026-11-22 "Yes, in this build".

    uv run --locked python library/l1/orf-design/tag_refold/publish_native_folds.py \\
        --folds ~/.cache/deliverome-analysis/tagged-folds/natives-2026-11-22 \\
        --queues tag_refold_queue_natives_2026_11_22.tsv,tag_refold_queue_natives_2026_11_22_mane.tsv,tag_refold_queue_natives_2026_11_22_refill.tsv \\
        --prefix sprint4/s4e9-orf-design-fg-panel/predicted-structures/tagged/2026-11-22-natives \\
        --release natives-2026-11-22 --members l1-p3178,...
"""

from __future__ import annotations

import argparse
import hashlib
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
EXP = HERE.parent
sys.path.insert(0, str(HERE))

from coreceptor_additions_publish import append, tsv  # noqa: E402
from native_models import pdb_sequence  # noqa: E402

INDEX = EXP / "tag_refold_native_models.tsv"
PAYLOAD = EXP / "library_structure_payload.tsv"
IMPLEMENTED = "Claude Code session l1-review-export-scope (subagent), 2026-11-22"


def sha(s: str) -> str:
    return hashlib.sha256(s.encode()).hexdigest()


def mean_plddt(pdb: Path) -> float:
    b = [
        float(ln[60:66])
        for ln in pdb.read_text().splitlines()
        if ln.startswith("ATOM") and ln[12:16].strip() == "CA"
    ]
    return round(sum(b) / len(b), 2)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--folds", required=True)
    ap.add_argument("--queues", required=True, help="comma-separated queue files (in s4e9)")
    ap.add_argument("--prefix", required=True, help="fresh processed prefix, no trailing slash")
    ap.add_argument("--release", required=True, help="payload release label")
    ap.add_argument("--members", required=True, help="candidate ids to give a payload row")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    from deliverome_analysis.cloning import library_design as ld
    from deliverome_analysis.core.s3 import processed
    from deliverome_analysis.surfaceome import topology

    folds = Path(a.folds).expanduser()
    natives = {}
    for q in a.queues.split(","):
        for r in tsv(EXP / q)[1]:
            if r["kind"] == "native":
                natives[r["fold_id"]] = r
    existing = set(processed.list_keys(a.prefix + "/"))
    if existing:
        raise SystemExit(f"processed://{a.prefix}/ is in use; publish under a fresh prefix")
    index_cols, index_rows = tsv(INDEX)
    have = {r["native_sha256"] for r in index_rows}
    new_index, new_keys = [], []
    for fid, q in sorted(natives.items()):
        pdb, scores = folds / f"{fid}.pdb", folds / f"{fid}.scores.json.gz"
        assert pdb.exists() and scores.exists(), f"{fid}: fold missing from {folds}"
        assert sha(q["native_seq"]) == q["native_sha256"], fid
        assert pdb_sequence(pdb) == q["native_seq"], f"{fid}: model is not the queued protein"
        for f in (pdb, scores):
            key = f"{a.prefix}/{f.name}"
            assert key not in existing, f"processed://{key} exists; publish under a fresh prefix"
            new_keys.append((f, key))
        if q["native_sha256"] not in have:
            have.add(q["native_sha256"])
            new_index.append(
                {
                    "native_sha256": q["native_sha256"],
                    "source": "colabfold",
                    "ref": f"processed://{a.prefix}/{pdb.name}",
                    "match": "identical",
                    "via": f"untagged native fold {a.release} of {q['design_id']} "
                    f"({q['gene_symbol']}); {IMPLEMENTED}",
                }
            )

    # payload rows, one chosen model per member, by sequence
    model = {}
    for r in index_rows + new_index:
        h = r["native_sha256"]
        rank = 0 if r["source"] == "alphafold_db" and r["match"] == "identical" else 1
        if r["match"] == "identical" and (h not in model or rank < model[h][0]):
            model[h] = (rank, r)
    local = {r["ref"]: folds / Path(r["ref"]).name for r in new_index if r["source"] == "colabfold"}
    prot = ld.designed_proteins(EXP / "candidate_manifest_tranches.tsv", library=ld.LIBRARY_L1)
    paths = {r["candidate_id"]: r for r in tsv(EXP / "l1_tag_paths.tsv")[1]}
    pay_cols, pay_rows = tsv(PAYLOAD)
    prev_acc = {r["candidate_id"]: r["uniprot_acc"] for r in pay_rows if r["uniprot_acc"]}
    snap = topology.load_snapshot()
    acc_by_sha = {
        h: a_
        for h, a_ in zip(snap["sequence_sha256"], snap["uniprot_acc_full"], strict=True)
        if isinstance(a_, str) and not a_.startswith("ENS")
    }
    payload, missing = [], []
    for cid in a.members.split(","):
        seq = prot[cid]["sequence"]
        h = sha(seq)
        assert paths[cid]["sequence_sha256"] == h, f"{cid}: tag paths describe another protein"
        if h not in model:
            missing.append(cid)
            continue
        r = model[h][1]
        t = paths[cid]
        base = {
            "release": a.release,
            "candidate_id": cid,
            "gene_symbol": t["gene_symbol"],
            "species": t["species"],
            "candidate_set": t["candidate_set"],
            "r2_key": "",
            "proteoform_match": "identical",
            "residue_offset": "0",
            "ordered_length_aa": str(len(seq)),
            "model_length_aa": str(len(seq)),
            "is_chosen": "1",
        }
        if r["source"] == "alphafold_db":
            acc = r["ref"].rsplit("AF-", 1)[1].split("-F1-")[0]
            payload.append(
                base
                | {
                    "uniprot_acc": acc,
                    "source": "alphafold_db",
                    "s3_key": "",
                    "model_url": r["ref"],
                    "mean_plddt": "",
                    "note": f"AlphaFold DB v6 model identical to the designed protein ({a.release})",
                }
            )
        else:
            pdb = local.get(r["ref"])
            if pdb is None:  # an earlier ColabFold fold of the same protein
                pdb = Path(processed.download(r["ref"].removeprefix("processed://")))
            assert pdb_sequence(pdb) == seq, f"{cid}: model is not the designed protein"
            payload.append(
                base
                | {
                    "uniprot_acc": prev_acc.get(cid) or acc_by_sha.get(h, ""),
                    "source": "colabfold",
                    "s3_key": r["ref"],
                    "model_url": "",
                    "mean_plddt": f"{mean_plddt(pdb):.2f}",
                    "proteoform_match": "ours_by_construction",
                    "note": f"folded from the designed protein (untagged native fold, {a.release})",
                }
            )
    if missing:
        raise SystemExit(f"no sequence-identical model for {missing}")
    print(f"{len(new_keys)} objects, {len(new_index)} index rows, {len(payload)} payload rows")
    if a.dry_run:
        return
    for f, key in new_keys:
        processed.upload(f, key)
    print("index", append(INDEX, new_index, index_cols))
    print("payload", append(PAYLOAD, payload, pay_cols))
    print(f"published under processed://{a.prefix}/")


if __name__ == "__main__":
    main()
