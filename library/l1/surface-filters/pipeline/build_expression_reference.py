#!/usr/bin/env python3
"""Re-derive the deterministic expression reference from raw sources.

This is the *port* of deliverome-internal's
`data/analysis/genomewide_expression_filtering/build_genomewide_expression_list.py`
into this repo — re-derived from raw (NOT reused), with a download-traceability
manifest and an added **cell-line median** column (for co-receptor ubiquity).

Sources (versions pinned by URL; recorded to download_traceability.json):
  * GTEx v10 gene median TPM  (max across 68 tissues)          -> max_tpm_gtex
  * HPA RNA tissue consensus nTPM (max across 51 tissues)      -> max_ntpm_hpa_tissue
  * HPA RNA cell line TPM (max AND median across 1,206 lines)  -> max_tpm_hpa_celline,
                                                                  median_tpm_hpa_celline
  * TCGA (UCSC Xena TcgaTargetGtex RSEM gene TPM, log2(TPM+1e-3);
    TCGA-* samples only; max)                                  -> max_tpm_tcga
  * HGNC complete set (protein-coding filter + Ensembl collapse)

Derived: max_expression_any = max of the four source maxima; passes_expression_
threshold = max_expression_any >= 10. See external/expression_tpm/PROVENANCE.md
for units (nTPM vs TPM) and caveats.

Promotion: intended to move to `scripts/build_expression_reference.py` and the
outputs to `s3://deliverome-external/expression/<version>/` (admin-write; the
scientist role is read-only on external). Run: `uv run python build_expression_reference.py`
"""

from __future__ import annotations

import gzip
import hashlib
import io
import json
import statistics
import urllib.request
import zipfile
from datetime import UTC, datetime
from pathlib import Path

HERE = Path(__file__).resolve().parent
RAW = HERE / "raw"
OUT = HERE.parent / "external" / "expression_tpm"
TPM_THRESHOLD = 10.0

# NOTE (2026-08-04): all four source URLs verified 200. GTEx v10 lives under
# `adult-gtex/bulk-gex/v10/rna-seq/` (HYPHENS, not underscores — that was the 404).
# The download step skips any file already present in pipeline/raw/.
SOURCES = {
    "gtex_v10_gene_median_tpm": "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz",
    "hpa_rna_tissue_consensus": "https://www.proteinatlas.org/download/tsv/rna_tissue_consensus.tsv.zip",
    "hpa_rna_celline": "https://www.proteinatlas.org/download/tsv/rna_celline.tsv.zip",
    "tcga_target_gtex_rsem_gene_tpm": "https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/TcgaTargetGtex_rsem_gene_tpm.gz",
    "hgnc_complete_set": "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt",
}


def download(name: str, url: str, trace: dict) -> Path:
    RAW.mkdir(parents=True, exist_ok=True)
    dest = RAW / url.split("/")[-1]
    if not dest.exists():
        print(f"  downloading {name} <- {url}")
        urllib.request.urlretrieve(url, dest)
    sha = hashlib.sha256(dest.read_bytes()).hexdigest()
    trace[name] = {
        "url": url,
        "file": dest.name,
        "bytes": dest.stat().st_size,
        "sha256": sha,
        "retrieved_at": datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ"),
    }
    return dest


def _open_maybe_zip(path: Path, inner_suffix: str = ".tsv"):
    if path.suffix == ".zip":
        zf = zipfile.ZipFile(path)
        member = next(n for n in zf.namelist() if n.endswith(inner_suffix))
        return io.TextIOWrapper(zf.open(member), encoding="utf-8")
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return open(path, encoding="utf-8")


def gtex_max(path):
    import pandas as pd

    df = pd.read_csv(path, sep="\t", skiprows=2, compression="gzip")
    df["ensembl_base"] = df["Name"].str.split(".").str[0]
    cols = [c for c in df.columns if c not in ("Name", "Description", "ensembl_base")]
    df[cols] = df[cols].apply(pd.to_numeric, errors="coerce")
    return df.assign(max_tpm_gtex=df[cols].max(axis=1))[["ensembl_base", "max_tpm_gtex"]]


def hpa_tissue_max(path):
    import pandas as pd

    df = pd.read_csv(_open_maybe_zip(path), sep="\t")
    df["ensembl_base"] = df["Gene"].str.split(".").str[0]
    return (
        df.groupby("ensembl_base", as_index=False)["nTPM"]
        .max()
        .rename(columns={"nTPM": "max_ntpm_hpa_tissue"})
    )


def hpa_celline_max_median(path):
    """Stream the (gene-grouped) cell-line matrix -> per-gene max AND median TPM."""
    import pandas as pd

    rows, cur, vals = [], None, []
    with _open_maybe_zip(path) as fh:
        fh.readline()
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 4:
                continue
            g = p[0].split(".")[0]
            if g != cur:
                if cur and vals:
                    rows.append((cur, max(vals), statistics.median(vals)))
                cur, vals = g, []
            try:
                vals.append(float(p[3]))
            except ValueError:
                pass
    if cur and vals:
        rows.append((cur, max(vals), statistics.median(vals)))
    return pd.DataFrame(
        rows, columns=["ensembl_base", "max_tpm_hpa_celline", "median_tpm_hpa_celline"]
    )


def tcga_max(path):
    import pandas as pd

    best = {}
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = [i for i, s in enumerate(header[1:], 1) if s.startswith("TCGA-")]
        for line in fh:
            p = line.rstrip("\n").split("\t")
            base = p[0].split(".")[0]
            m = (2 ** max(float(p[i]) for i in idx)) - 0.001
            if base not in best or m > best[base]:
                best[base] = m
    return pd.DataFrame([{"ensembl_base": b, "max_tpm_tcga": v} for b, v in best.items()])


def main():
    import pandas as pd

    trace = {}
    paths = {n: download(n, u, trace) for n, u in SOURCES.items()}
    print("processing…")
    df = gtex_max(paths["gtex_v10_gene_median_tpm"])
    for other in (
        hpa_tissue_max(paths["hpa_rna_tissue_consensus"]),
        hpa_celline_max_median(paths["hpa_rna_celline"]),
        tcga_max(paths["tcga_target_gtex_rsem_gene_tpm"]),
    ):
        df = df.merge(other, on="ensembl_base", how="outer")
    val_cols = ["max_tpm_gtex", "max_ntpm_hpa_tissue", "max_tpm_hpa_celline", "max_tpm_tcga"]
    df["max_expression_any"] = df[val_cols].max(axis=1)
    df["passes_expression_threshold"] = df["max_expression_any"] >= TPM_THRESHOLD
    # HGNC protein-coding collapse
    hg = pd.read_csv(paths["hgnc_complete_set"], sep="\t", dtype=str, low_memory=False)
    hg = hg[(hg["status"] == "Approved") & (hg["locus_group"] == "protein-coding gene")].copy()
    hg["ensembl_base"] = hg["ensembl_gene_id"].str.split(".").str[0]
    df = df.merge(
        hg[["ensembl_base", "symbol"]].rename(columns={"symbol": "gene_symbol_canonical"}),
        on="ensembl_base",
        how="left",
    )
    OUT.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT / "genomewide_expression_tpm_rederived.csv", index=False)
    (OUT / "download_traceability.json").write_text(json.dumps(trace, indent=2))
    print(
        f"wrote {OUT / 'genomewide_expression_tpm_rederived.csv'} ({len(df)} genes) + download_traceability.json"
    )


if __name__ == "__main__":
    main()
