#!/usr/bin/env python3
"""Build an authoritative ENST -> UniProt-isoform-ID map for the genes with
isoforms in Willow's GPCR screen (Tedman 2025), so screen transcripts can be
matched to our catalog's UniProt-isoform-keyed topology records.

Source: UniProt REST `xref_ensembl` (each Ensembl transcript cross-ref carries an
`isoformId` like P35372-3). Cached per-accession under uniprot_cache/.

Run: uv run --locked python build_isoform_map.py
Out: willow_enst_uniprot_isoform_map.json  {ENST(no version): {gene, acc, isoform_id}}
"""
# ruff: noqa: E741

import json
import os
import time
import urllib.request

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
XLSX = os.path.join(HERE, "external/willow_gpcr_screen/tedman2025_gpcr_screen_media-2.xlsx")
CACHE = os.path.join(HERE, "external/willow_gpcr_screen/uniprot_cache")
OUT = os.path.join(HERE, "willow_enst_uniprot_isoform_map.json")
os.makedirs(CACHE, exist_ok=True)

recs = {json.loads(l)["gene"]: json.loads(l) for l in open(os.path.join(HERE, "cohort.jsonl"))}

iso = pd.read_excel(XLSX, sheet_name="Isoforms")
iso_genes = sorted(set(iso["Isoform Transcript"].str.split("_").str[0]))
# accession per gene from our cohort
acc_of = {g: recs[g]["uniprot"] for g in iso_genes if g in recs}
print(f"{len(iso_genes)} isoform-screened genes; {len(acc_of)} have a cohort UniProt accession")


def fetch(acc):
    cp = os.path.join(CACHE, f"{acc}.json")
    if os.path.exists(cp):
        return json.load(open(cp))
    url = f"https://rest.uniprot.org/uniprotkb/{acc}.json?fields=xref_ensembl"
    req = urllib.request.Request(
        url, headers={"User-Agent": "deliverome-analysis/1.0 (becca@deliverome.org)"}
    )
    d = json.load(urllib.request.urlopen(req, timeout=30))
    json.dump(d, open(cp, "w"))
    time.sleep(0.15)
    return d


emap = {}
miss = []
for g, acc in acc_of.items():
    try:
        d = fetch(acc)
    except Exception as e:
        miss.append((g, acc, str(e)[:60]))
        continue
    for x in d.get("uniProtKBCrossReferences", []):
        if x.get("database") != "Ensembl":
            continue
        enst = x.get("id", "").split(".")[0]
        iso_id = x.get("isoformId")  # e.g. P35372-3 ; None if entry has no isoforms
        if enst:
            emap[enst] = {"gene": g, "acc": acc, "isoform_id": iso_id}

json.dump(emap, open(OUT, "w"), indent=0)
n_iso = sum(1 for v in emap.values() if v["isoform_id"])
print(f"mapped {len(emap)} ENST ({n_iso} carry an isoformId); fetch failures: {len(miss)}")
if miss:
    print("  misses:", miss[:10])

# coverage against the screen's isoform transcripts
iso["enst"] = iso["Isoform Transcript"].str.extract(r"(ENST\d+)")
covered = iso["enst"].isin(emap).sum()
print(f"screen isoform transcripts covered by map: {covered}/{len(iso)}")
