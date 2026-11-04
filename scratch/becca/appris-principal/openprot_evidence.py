"""Attach OpenProt's MS and Ribo-seq evidence to each alternate form in the universe.

OpenProt re-annotates the proteome from Ensembl, RefSeq and UniProt together and scores
every protein on two axes we do not otherwise have in one place:

* **MS score** -- unique peptides summed per study, so it aggregates across studies rather
  than reporting one build the way PeptideAtlas does;
* **TE (Translation Event) score** -- the number of studies with a significant detection of
  **translation**, which is ribosome-profiling evidence. This is the gap the GENCODE
  Ribo-seq ORF catalogue could not fill, because that catalogue holds only ORFs *outside*
  annotated CDSs while OpenProt scores the annotated forms themselves.

`NPFFR2` is the worked example: the 522-residue UniProt canonical we used to build from
scores TE=0 and is not even typed as a RefProt, while the 420-residue form that MANE, CCDS,
APPRIS and Tedman all use scores TE=3. That is independent ribosome-profiling support for
the call the reference policy already made.

**Caveat on currency.** OpenProt 2.2 is built on Ensembl GRCh38.106 and UniProt 2022_06,
several years behind the Ensembl 115 / UniProt 2026_03 used everywhere else here. A form
missing from OpenProt may simply postdate it, so absence is not evidence.

    uv run --locked python openprot_evidence.py
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

OPENPROT_RELEASE = "2.2 (Ensembl GRCh38.106, UniProt 2022_06)"


def load_openprot(path: Path) -> dict[str, dict]:
    """{accession: best scores}. One row per transcript, so scores are taken at their max."""
    out: dict[str, dict] = {}
    with path.open() as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for row in rd:
            acc = row["protein accession"]
            if not acc:
                continue
            try:
                ms = int(row["MS score"] or 0)
                te = int(row["TE Translation Event score"] or 0)
            except ValueError:
                continue
            # OpenProt carries versioned accessions (ENSP00000307822.7, NP_004876.3) while
            # our side is unversioned. Indexing both spellings is what makes the Ensembl and
            # RefSeq joins work at all -- keyed on the versioned form alone, they matched 0.
            for key in {acc, acc.split(".")[0]}:
                rec = out.setdefault(
                    key,
                    {
                        "type": row["protein type"],
                        "length": row["protein length (a.a.)"],
                        "gene": row["gene symbol"],
                        "ms": 0,
                        "te": 0,
                    },
                )
                rec["ms"] = max(rec["ms"], ms)
                rec["te"] = max(rec["te"], te)
    return out


def enst_to_ensp(pep_gz: Path) -> dict[str, str]:
    out = {}
    with gzip.open(pep_gz, "rt") as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            ensp = line[1:].split()[0].split(".")[0]
            m = re.search(r"\btranscript:(\S+)", line)
            if m:
                out[m.group(1).split(".")[0]] = ensp
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--openprot-dir", default="~/.cache/deliverome/openprot")
    ap.add_argument("--cache-dir", default="~/.cache/deliverome/appris")
    ap.add_argument("--universe", default="isoform_universe.tsv")
    args = ap.parse_args()

    opdir = Path(args.openprot_dir).expanduser()
    tsvs = sorted(opdir.glob("*min_2_pep*.tsv"))
    if not tsvs:
        sys.exit(f"no OpenProt TSV in {opdir}")
    print(f"reading {tsvs[0].name}")
    op = load_openprot(tsvs[0])
    print(f"OpenProt accessions: {len(op)}")

    e2p = enst_to_ensp(Path(args.cache_dir).expanduser() / "Homo_sapiens.GRCh38.pep.all.fa.gz")
    rows = list(csv.DictReader(open(args.universe), delimiter="\t"))

    matched = Counter()
    for r in rows:
        hit = None
        if r["uniprot"] and r["uniprot"] in op:
            hit = op[r["uniprot"]]
            matched["by uniprot accession"] += 1
        if hit is None:
            for enst in (r["ensembl_transcripts"] or "").split(";"):
                ensp = e2p.get(enst)
                if ensp and ensp in op:
                    hit = op[ensp]
                    matched["by ensembl protein"] += 1
                    break
        if hit is None:
            matched["no OpenProt row"] += 1
        r["openprot_type"] = hit["type"] if hit else ""
        r["openprot_ms_score"] = hit["ms"] if hit else ""
        r["openprot_te_score"] = hit["te"] if hit else ""
        r["openprot_translated"] = (
            "" if hit is None else str(hit["te"] > 0)
        )

    out = Path("isoform_universe_openprot.tsv")
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {out} ({len(rows)} forms)")

    scored = [r for r in rows if r["openprot_te_score"] != ""]
    translated = [r for r in scored if int(r["openprot_te_score"]) > 0]
    summary = {
        "openprot_release": OPENPROT_RELEASE,
        "forms": len(rows),
        "match": dict(matched.most_common()),
        "with_openprot_scores": len(scored),
        "ribo_seq_translated": len(translated),
        "rescued_by_riboseq": {
            # forms our peptide tiers could not call, that OpenProt says are translated
            tier: sum(1 for r in translated if r["tier"] == tier)
            for tier in sorted({r["tier"] for r in translated})
        },
        "te_score_by_tier": {
            tier: sum(1 for r in scored if r["tier"] == tier and int(r["openprot_te_score"]) > 0)
            for tier in sorted({r["tier"] for r in scored})
        },
    }
    Path("openprot_evidence_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
