#!/usr/bin/env python3
"""Backfill protrusion_14 and d_other_domain for all L1 internal CDSs, and build the ProDomino queue.

Reuses the exact geometry and InterPro domain-seam functions from
`library/l1/tag-site-benchmarking/structure_global_features.py` against each
L1 internal CDS's native structure model (AlphaFold DB v6 or ColabFold native PDB).

Outputs (under tag-site-models/data/, gitignored):
    l1_structure_global_features.tsv   per-residue protrusion_14 and d_other_domain
    l1_prodomino_queue.json            L1 internal sequences missing ProDomino (+ CLDN17 parity check)
    l1_backfill_summary.json           coverage summary

    uv run --locked python library/l1/orf-design/backfill_l1_structure_global.py
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO / "library/l1/tag-site-benchmarking"))
import structure_global_features as SGF  # noqa: E402

DATA = HERE / "tag-site-models" / "data"
DEFAULT_BUNDLE = Path.home() / "l1-bundle-2026-11-21b"
INTERPRO_URL = (
    "https://www.ebi.ac.uk/interpro/api/entry/all/protein/uniprot/{acc}/?page_size=200&format=json"
)


def resolve_pdb(source: str, ref: str, bundle: Path) -> Path | None:
    if pd.isna(source) or pd.isna(ref):
        return None
    if source == "alphafold_db":
        fname = ref.rsplit("/", 1)[-1]
        p = bundle / "s3/external/alphafold-db/v6" / fname
        return p if p.is_file() else None
    if source == "colabfold":
        rel = ref.removeprefix("processed://")
        p = bundle / "s3/processed" / rel
        return p if p.is_file() else None
    return None


def fetch_interpro(acc: str, cache_dir: Path) -> dict | None:
    cache_file = cache_dir / f"{acc}.json"
    if cache_file.is_file():
        try:
            return json.loads(cache_file.read_text())
        except json.JSONDecodeError:
            pass
    url = INTERPRO_URL.format(acc=acc)
    req = urllib.request.Request(url, headers={"User-Agent": "deliverome-analysis/1.0"})
    try:
        with urllib.request.urlopen(req, timeout=20) as resp:
            raw = resp.read().decode("utf-8")
            data = json.loads(raw) if raw.strip() else {"results": []}
    except urllib.error.HTTPError, urllib.error.URLError, TimeoutError, json.JSONDecodeError:
        data = {"results": []}
    cache_file.write_text(json.dumps(data))
    return data


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bundle", type=Path, default=DEFAULT_BUNDLE)
    ap.add_argument("--skip-interpro", action="store_true")
    args = ap.parse_args()

    DATA.mkdir(parents=True, exist_ok=True)
    ip_cache = DATA / "interpro_cache"
    ip_cache.mkdir(exist_ok=True)

    tp = pd.read_csv(HERE / "l1_tag_paths.tsv", sep="\t", low_memory=False)
    internal = tp[tp["placement"] == "internal"].copy()
    nm = pd.read_csv(HERE / "tag_refold_native_models.tsv", sep="\t", low_memory=False)
    tgt = pd.read_csv(REPO / "library/l1/tag-site-benchmarking/l1_adopted_targets.tsv", sep="\t")
    sha_to_acc = dict(zip(tgt["sequence_sha256"], tgt["uniprot_acc"], strict=False))

    m = internal.merge(nm, left_on="sequence_sha256", right_on="native_sha256", how="left")
    m["uniprot_acc"] = m["sequence_sha256"].map(sha_to_acc)
    # Also parse UniProt acc from AFDB filename when not in l1_adopted_targets
    for i, r in m[m["uniprot_acc"].isna() & (m["source"] == "alphafold_db")].iterrows():
        fname = str(r["ref"]).rsplit("/", 1)[-1]
        if fname.startswith("AF-") and "-F1-" in fname:
            m.at[i, "uniprot_acc"] = fname.split("-")[1]

    # Check which accessions already have ProDomino in all_extracellular_targets.tsv
    allx = pd.read_csv(DATA / "all_extracellular_targets.tsv", sep="\t", low_memory=False)
    has_prodom_accs = set(allx.loc[allx["prodomino"].notna(), "uniprot_acc"].dropna().unique())

    # 1. Fetch InterPro JSONs concurrently for unique UniProt accessions
    unique_accs = sorted(
        {
            a
            for a in m["uniprot_acc"].dropna().unique()
            if isinstance(a, str) and len(a) >= 6 and not a.startswith("l1-")
        }
    )
    interpro_by_acc: dict[str, dict] = {}
    if not args.skip_interpro:
        print(f"Fetching/loading InterPro domains for {len(unique_accs)} UniProt accessions...")
        with ThreadPoolExecutor(max_workers=12) as pool:
            futs = {pool.submit(fetch_interpro, a, ip_cache): a for a in unique_accs}
            done = 0
            for fut in as_completed(futs):
                a = futs[fut]
                interpro_by_acc[a] = fut.result() or {"results": []}
                done += 1
                if done % 100 == 0 or done == len(unique_accs):
                    print(f"  InterPro {done}/{len(unique_accs)}")

    # 2. Compute protrusion_14 + d_other_domain from local PDBs, and extract sequences
    out_rows = []
    queue_items = []
    seen_queue_sha: set[str] = set()
    n_pdb_ok = 0
    n_seam_ok = 0
    n_sha_match = 0

    for r in m.itertuples(index=False):
        pdb = resolve_pdb(r.source, r.ref, args.bundle)
        if pdb is None:
            continue
        chains = SGF.chain_of(pdb)
        if not chains:
            continue
        res = SGF.residues_with_ca(chains[0])
        if not res:
            continue
        n_pdb_ok += 1
        seq = SGF.sequence(res)
        seq_sha = hashlib.sha256(seq.encode()).hexdigest()
        if seq_sha == r.sequence_sha256:
            n_sha_match += 1

        atoms = SGF.heavy_atoms(res)
        prot14 = SGF.protrusion(res, atoms, SGF.PROTRUSION_RADIUS)

        acc = r.uniprot_acc if pd.notna(r.uniprot_acc) else r.sequence_sha256[:12]
        d_dom: dict[int, float] = {}
        if acc in interpro_by_acc:
            _src, instances = SGF.domain_instances(interpro_by_acc[acc])
            if len(instances) >= 2:
                max_res = max(int(x.id[1]) for x in res)
                own, ann = SGF.assign_domains(instances, max_res)
                d_dom = SGF.d_other_domain(res, own, ann)
                if d_dom:
                    n_seam_ok += 1

        for x in res:
            rn = int(x.id[1])
            out_rows.append(
                {
                    "sequence_sha256": r.sequence_sha256,
                    "uniprot_acc": acc,
                    "gene_symbol": r.gene_symbol,
                    "residue": rn,
                    "protrusion_14": prot14.get(rn),
                    "d_other_domain": round(d_dom[rn], 3) if rn in d_dom else None,
                }
            )

        # Queue for ProDomino if missing (or if CLDN17 P56750 parity check)
        if (
            acc not in has_prodom_accs or acc == "P56750"
        ) and r.sequence_sha256 not in seen_queue_sha:
            seen_queue_sha.add(r.sequence_sha256)
            queue_items.append(
                {
                    "sequence_sha256": r.sequence_sha256,
                    "uniprot_acc": acc,
                    "gene_symbol": r.gene_symbol,
                    "sequence": seq,
                    "is_parity_check": acc == "P56750",
                }
            )

    out_df = pd.DataFrame(out_rows)
    out_df.to_csv(DATA / "l1_structure_global_features.tsv", sep="\t", index=False)
    (DATA / "l1_prodomino_queue.json").write_text(json.dumps(queue_items, indent=2) + "\n")

    summary = {
        "n_internal_l1": len(m),
        "n_pdb_found": n_pdb_ok,
        "n_pdb_seq_sha_exact_match": n_sha_match,
        "n_proteins_with_protrusion_14": int(out_df["sequence_sha256"].nunique()),
        "n_proteins_with_d_other_domain": n_seam_ok,
        "n_residues_scored": len(out_df),
        "n_already_have_prodomino": int(m["uniprot_acc"].isin(has_prodom_accs).sum()),
        "n_queued_for_prodomino": len(queue_items),
    }
    (DATA / "l1_backfill_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
