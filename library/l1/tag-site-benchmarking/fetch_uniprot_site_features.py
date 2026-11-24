#!/usr/bin/env python3
"""UniProt residue annotations for the corpus, the L1 members and the screened proteins.

Pulls disulfide bonds, glycosylation, binding / active / metal / functional
sites and the sequence from the EBI Proteins API (UniProt REST was returning
500s when this was written; the EBI API serves the same UniProtKB features),
one accession at a time with a 0.3 s pause. Writes one row per feature:

    uniprot_site_features.tsv   uniprot_acc, type, description, begin, end

A DISULFID feature spans its two cysteines (begin and end are the partners);
an interchain bond has begin == end. Accessions with no entry are listed in
uniprot_site_features.json rather than silently dropped.

    uv run --locked python fetch_uniprot_site_features.py
"""

from __future__ import annotations

import json
import re
import time
import urllib.error
import urllib.request
from datetime import UTC, datetime
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
OUT = HERE / "uniprot_site_features.tsv"
META = HERE / "uniprot_site_features.json"
TYPES = ["DISULFID", "CARBOHYD", "BINDING", "ACT_SITE", "METAL", "SITE", "LIPID"]
API = "https://www.ebi.ac.uk/proteins/api/features/{acc}?types=" + ",".join(TYPES)
SCREENED = ["P78348", "P22001", "P56373", "P63252", "P35561"]  # P35561: mouse Kir2.1
ENDONB = ["P02786", "P05556", "P18084", "P30530", "Q8N131"]
# UniProt accession grammar; designed ORFs (l1-pNNNN) have no entry to fetch
ACCESSION = re.compile(r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$")
CACHE = HERE / "data" / "uniprot_features"  # gitignored; one JSON per accession, so reruns resume


def accessions() -> list[str]:
    sf = pd.read_csv(HERE / "site_features.tsv", sep="\t", usecols=["uniprot_acc", "label"])
    corpus = set(sf.loc[sf["label"] == 1, "uniprot_acc"])
    members = set(pd.read_csv(HERE / "l1_member_targets.tsv", sep="\t")["uniprot_acc"])
    return sorted(a for a in corpus | members | set(SCREENED) | set(ENDONB) if ACCESSION.match(a))


def fetch(acc: str) -> dict | None:
    hit = CACHE / f"{acc}.json"
    if hit.is_file():
        return json.loads(hit.read_text()) or None
    d = _get(acc)
    CACHE.mkdir(parents=True, exist_ok=True)
    hit.write_text(json.dumps(d or {}))
    time.sleep(0.3)
    return d


def _get(acc: str) -> dict | None:
    req = urllib.request.Request(
        API.format(acc=acc),
        headers={"Accept": "application/json", "User-Agent": "deliverome-analysis s5e13"},
    )
    for attempt in range(4):
        try:
            with urllib.request.urlopen(req, timeout=60) as r:
                return json.loads(r.read())
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None
            time.sleep(2 + 3 * attempt)
        except urllib.error.URLError, TimeoutError:
            time.sleep(2 + 3 * attempt)
    raise RuntimeError(f"EBI Proteins API failed for {acc}")


def main() -> int:
    accs = accessions()
    rows, missing, seqs = [], [], {}
    for i, acc in enumerate(accs):
        d = fetch(acc)
        if d is None:
            missing.append(acc)
            continue
        seqs[acc] = d.get("sequence", "")
        for f in d.get("features", []):
            rows.append(
                {
                    "uniprot_acc": acc,
                    "type": f.get("type"),
                    # description before the positions: it is often empty, and an
                    # empty last field is stripped by the trailing-whitespace hook
                    "description": (f.get("description") or "").replace("\t", " ").strip(),
                    "begin": f.get("begin"),
                    "end": f.get("end"),
                }
            )
        if i % 50 == 0:
            print(f"{i}/{len(accs)} {acc}", flush=True)
    pd.DataFrame(rows).to_csv(OUT, sep="\t", index=False)
    META.write_text(
        json.dumps(
            {
                "fetched_at": datetime.now(UTC).isoformat(),
                "source": API.split("{")[0],
                "types": TYPES,
                "n_accessions": len(accs),
                "no_entry": missing,
                "sequences": seqs,
            },
            indent=1,
        )
        + "\n"
    )
    print(f"wrote {OUT.name}: {len(rows)} features on {len(accs) - len(missing)} accessions")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
