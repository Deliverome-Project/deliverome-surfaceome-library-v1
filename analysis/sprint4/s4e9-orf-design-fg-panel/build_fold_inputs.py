"""Sequences for the Library-1 members AlphaFold DB has no model for.

`l1_structure_gap.py` finds them; this fetches what they actually are so they can be
folded. Identifiers come in two flavours and neither can be assumed:

* UniProt accessions (the human isoforms resolve to recently-added TrEMBL entries, which
  have sequences even though AlphaFold DB predates them)
* Ensembl protein ids (`ENSMFAP…` cyno, `ENSMUSP…` mouse), which UniProt does not serve

Each resolved accession is re-checked against AlphaFold DB before being written out, so a
candidate that was only missing because a gene-symbol lookup failed — TPO is the one here —
drops out instead of being folded needlessly.

    uv run --locked python build_fold_inputs.py
"""

from __future__ import annotations

import csv
import json
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).parent
GAP = HERE / "l1_structure_gap.tsv"
FASTA = HERE / "fold_inputs.fasta"
MANIFEST = HERE / "fold_inputs.tsv"
WORKERS = 8


def uniprot_seq(acc: str) -> str | None:
    try:
        d = json.load(
            urllib.request.urlopen(
                f"https://rest.uniprot.org/uniprotkb/{acc}.json?fields=sequence", timeout=30
            )
        )
        return d["sequence"]["value"]
    except Exception:
        return None


def ensembl_seq(pid: str, *, attempts: int = 6) -> str | None:
    """Ensembl REST rate-limits hard, so this retries on 429 and is called serially.

    Running it through the thread pool silently returned nothing for all 30 Ensembl ids;
    the requests were being throttled, not the ids being wrong.
    """
    url = f"https://rest.ensembl.org/sequence/id/{pid}?type=protein"
    req = urllib.request.Request(url, headers={"Content-Type": "application/json"})
    for attempt in range(attempts):
        try:
            return json.load(urllib.request.urlopen(req, timeout=30))["seq"]
        except urllib.error.HTTPError as exc:
            # 429 is the documented rate limit; 5xx are transient server hiccups under
            # load and succeed on a retry, so neither is treated as a real failure.
            if exc.code == 429 or exc.code >= 500:
                wait = float(exc.headers.get("Retry-After", 1)) + 1.5 * (attempt + 1)
                time.sleep(wait)
                continue
            print(f"    {pid}: HTTP {exc.code}")
            return None
        except Exception as exc:  # report rather than swallow
            print(f"    {pid}: {type(exc).__name__}")
            time.sleep(1.0)
    print(f"    {pid}: gave up after {attempts} attempts")
    return None


def has_model(acc: str) -> bool:
    try:
        urllib.request.urlopen(
            f"https://alphafold.ebi.ac.uk/api/prediction/{acc}", timeout=25
        ).read()
        return True
    except Exception:
        return False


def fetch(row: dict) -> dict:
    acc = row["accession"]
    seq = ensembl_seq(acc) if acc.startswith("ENS") else uniprot_seq(acc)
    row["sequence"] = seq or ""
    row["length"] = len(seq or "")
    return row


def main() -> int:
    rows = [
        r
        for r in csv.DictReader(GAP.open(), delimiter="\t")
        if r["has_alphafold"] not in ("True", "true", "1") and r["accession"]
    ]
    print(f"{len(rows)} candidates with a resolvable identifier and no model")

    uni = [r for r in rows if not r["accession"].startswith("ENS")]
    ens = [r for r in rows if r["accession"].startswith("ENS")]
    with ThreadPoolExecutor(WORKERS) as pool:
        uni = list(pool.map(fetch, uni))
    print(f"  fetching {len(ens)} Ensembl sequences serially (REST rate limit)...")
    ens = [fetch(r) for r in ens]
    rows = uni + ens
    missing_seq = [r for r in rows if not r["sequence"]]
    rows = [r for r in rows if r["sequence"]]
    if missing_seq:
        print(
            f"  no sequence retrieved for {len(missing_seq)}: "
            f"{[r['accession'] for r in missing_seq]}"
        )

    # A second AlphaFold check on the RESOLVED accession: the first pass keyed some
    # candidates off a gene symbol, and a lookup miss there looks identical to a real gap.
    with ThreadPoolExecutor(WORKERS) as pool:
        present = list(pool.map(has_model, [r["accession"] for r in rows]))
    drop = [r["accession"] for r, p in zip(rows, present, strict=True) if p]
    if drop:
        print(f"  dropping {len(drop)} that do have a model after all: {drop}")
    rows = [r for r, p in zip(rows, present, strict=True) if not p]

    # One sequence per accession; the same ortholog can appear under several candidates.
    seen: dict[str, dict] = {}
    for r in rows:
        seen.setdefault(r["accession"], r)
    uniq = list(seen.values())

    FASTA.write_text(
        "".join(
            f">{r['accession']}|{r['symbol']}|{r['species']}|{r['set']}\n"
            + "\n".join(r["sequence"][i : i + 60] for i in range(0, len(r["sequence"]), 60))
            + "\n"
            for r in sorted(uniq, key=lambda r: r["accession"])
        )
    )
    cols = ["accession", "symbol", "species", "set", "identifier", "length"]
    with MANIFEST.open("w", newline="") as fh:
        w = csv.DictWriter(
            fh, fieldnames=cols, delimiter="\t", lineterminator="\n", extrasaction="ignore"
        )
        w.writeheader()
        w.writerows(sorted(uniq, key=lambda r: r["accession"]))

    lens = [r["length"] for r in uniq]
    print(
        f"\nwrote {FASTA.name}: {len(uniq)} sequences, "
        f"{min(lens)}-{max(lens)} aa (median {sorted(lens)[len(lens) // 2]}), "
        f"{sum(lens):,} residues total"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
