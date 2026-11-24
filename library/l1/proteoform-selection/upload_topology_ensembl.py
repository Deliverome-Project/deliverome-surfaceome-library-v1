"""Load the DeepTMHMM v1 calls for Ensembl-only forms into `topology_ensembl_public`.

Derives every column from the 3-line output rather than carrying numbers across from
another file, so the row and the prediction cannot drift apart. `hgnc_id` is resolved
through D1's own `gene_identifier_public` keyed on the **Ensembl gene ID**, never on the
gene symbol -- the `topology_release` notes record a past collision where symbol-keyed
resolution filed MYO18A's row under TIAF1, a different protein that shares the primary gene
name.

Two forms in the verdict file are marked `undecidable_no_reference_call` because their
gene's reference was never matched to a form. They are not undecidable: `TMC4` 706 aa and
`KIR3DL3` 410 aa **are** their gene's selected reference, so they are skipped as references
rather than loaded as alternates. `TMC4` 712 aa is a genuine alternate and compares against
706.

    uv run --locked python upload_topology_ensembl.py              # dry run
    uv run --locked python upload_topology_ensembl.py --apply
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import subprocess
import sys
import tempfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

TOPOLOGY_VERSION = "topo_2026_10_09_ensembl_forms"
TOOL_VERSION = "deeptmhmm-1.0.24"
COHORT = "human_ensembl_forms"
SPECIES = "human"
DATABASE = "surfaceome_public"
# Batch by statement size, not row count: each row carries a full protein sequence AND a
# full per-residue topology string, so a fixed row count overruns SQLITE_TOOBIG on
# long proteins and wastes round trips on short ones.
MAX_STATEMENT_BYTES = 90_000

# wrangler is pinned in accessible-surfaceome's package.json, and npx resolves it from the
# directory it runs in -- outside that repo it fails with an empty stderr, which is a silent
# non-zero exit rather than a message.
WRANGLER_CWD = Path("~/Git/accessible-surfaceome").expanduser()

SURFACE_LABELS = {"TM", "SP+TM"}
ORIENT = {"O": "extracellular", "I": "cytoplasmic"}


def d1(command: str, database: str = DATABASE) -> list[dict]:
    proc = subprocess.run(
        ["npx", "--yes", "wrangler", "d1", "execute", database, "--remote",
         "--command", command, "--json"],
        capture_output=True, text=True, cwd=WRANGLER_CWD,
    )
    if proc.returncode != 0:
        raise SystemExit(
            f"wrangler exited {proc.returncode}\nstderr: {proc.stderr[-1500:]}\n"
            f"stdout: {proc.stdout[-1500:]}"
        )
    payload = json.loads(proc.stdout)
    blocks = payload if isinstance(payload, list) else payload.get("result", [])
    return [r for b in blocks for r in (b.get("results") or [])]


def enst_to_ensp(pep_gz: Path) -> dict[str, str]:
    """{ENST: ENSP} from the Ensembl peptide FASTA headers.

    The table is keyed on the **protein** id, so the transcript id cannot stand in for it.
    """
    import gzip

    out: dict[str, str] = {}
    with gzip.open(pep_gz, "rt") as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            ensp = line[1:].split()[0].split(".")[0]
            m = re.search(r"\btranscript:(\S+)", line)
            if m:
                out[m.group(1).split(".")[0]] = ensp
    return out


def parse_3line(path: Path) -> dict[str, dict]:
    out: dict[str, dict] = {}
    lines = path.read_text().splitlines()
    for i in range(0, len(lines) - 2):
        if not lines[i].startswith(">"):
            continue
        header = lines[i][1:]
        name, label = header.rsplit("|", 1)
        out[name.strip()] = {
            "label": label.strip(),
            "sequence": lines[i + 1].strip(),
            "topology": lines[i + 2].strip(),
        }
    return out


def derive(topology: str) -> dict:
    """Counts and orientations straight from the per-residue string."""
    runs = lambda ch: sum(
        1 for j, c in enumerate(topology) if c == ch and (j == 0 or topology[j - 1] != ch)
    )
    non_signal = [c for c in topology if c != "S"]
    first = next((c for c in non_signal if c in ORIENT), None)
    last = next((c for c in reversed(non_signal) if c in ORIENT), None)
    return {
        "tm_helix_count": runs("M"),
        "beta_strand_count": runs("B"),
        "signal_peptide_length": topology.count("S"),
        "ecd_length_residues": topology.count("O"),
        "icd_length_residues": topology.count("I"),
        "n_terminal_orientation": ORIENT.get(first, "indeterminate"),
        "c_terminal_orientation": ORIENT.get(last, "indeterminate"),
    }


def sql_str(v) -> str:
    if v is None or v == "":
        return "NULL"
    return "'" + str(v).replace("'", "''") + "'"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--apply", action="store_true")
    ap.add_argument("--three-line", default="topology_shortlist.3line")
    ap.add_argument("--manifest", default="topology_shortlist.tsv")
    ap.add_argument("--verdicts", default="topology_verdicts.tsv")
    ap.add_argument("--universe", default="isoform_universe_scored.tsv")
    ap.add_argument("--cache-dir", default="~/.cache/deliverome/appris")
    args = ap.parse_args()

    e2p = enst_to_ensp(
        Path(args.cache_dir).expanduser() / "Homo_sapiens.GRCh38.pep.all.fa.gz"
    )

    pred = parse_3line(Path(args.three_line))
    manifest = {m["fasta_id"]: m for m in csv.DictReader(open(args.manifest), delimiter="\t")}
    verdicts = list(csv.DictReader(open(args.verdicts), delimiter="\t"))
    universe = list(csv.DictReader(open(args.universe), delimiter="\t"))

    refs = {r["gene"]: r for r in universe if r["is_reference"] == "True"}
    sel_len = {r["gene"]: r.get("length") for r in universe if r["is_reference"] == "True"}
    ensg_for_gene = {
        r["gene"]: (r.get("ensembl_gene_id") or "")
        for r in universe
        if r["is_reference"] == "True"
    }

    # ENSG for each form, taken from the universe row that produced it.
    form_ensg: dict[tuple[str, str], str] = {}
    for r in universe:
        form_ensg[(r["gene"], r["length"])] = r.get("ensembl_gene_id", "")

    # Only forms with NO UniProt accession belong in this table; the rest go to
    # topology_public through the existing uploader.
    # Driven from the manifest, keyed on fasta_id. An earlier version matched predictions
    # to manifest rows on (gene, length), but 63 genes carry two or more *different*
    # sequences of the same length -- those all resolved to one manifest row, got one ENSP,
    # and collapsed on the primary key, leaving the wrong sequence stored under it.
    verdict_by_id = {}
    for v in verdicts:
        verdict_by_id.setdefault((v["gene"], v["length"]), []).append(v)

    rows, skipped = [], Counter()
    for fasta_id, m in manifest.items():
        if m["is_reference"] == "True":
            skipped["reference, not an alternate"] += 1
            continue
        if m["uniprot"]:
            skipped["has a UniProt accession (belongs in topology_public)"] += 1
            continue
        p = pred.get(fasta_id)
        if not p:
            skipped["no prediction"] += 1
            continue
        if m["length"] == sel_len.get(m["gene"]):
            skipped["is the gene's reference, not an alternate"] += 1
            continue
        enst = (m["ensembl_transcripts"] or "").split(";")[0]
        ensp = e2p.get(enst, "")
        if not ensp:
            skipped["no ENSP for the transcript"] += 1
            continue
        d = derive(p["topology"])
        rows.append(
            {
                "topology_version": TOPOLOGY_VERSION,
                "cohort": COHORT,
                "ensembl_protein_id": ensp,
                "ensembl_transcript_id": enst,
                "ensembl_gene_id": "",
                "hgnc_id": "",
                "gene_symbol": m["gene"],
                "species": SPECIES,
                "reference_uniprot_acc_full": refs.get(m["gene"], {}).get("uniprot", ""),
                "sequence": p["sequence"],
                "protein_length": len(p["sequence"]),
                "deeptmhmm_label": p["label"],
                **d,
                "per_residue_topology": p["topology"],
                "predicted_surface_membrane": 1 if p["label"] in SURFACE_LABELS else 0,
                "predicted_secreted": 1 if p["label"] == "SP" else 0,
                "tool_version": TOOL_VERSION,
                "retrieved_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            }
        )

    # Every ENSP must be unique, or rows silently overwrite each other on the primary key.
    dupes = [k for k, n in Counter(r["ensembl_protein_id"] for r in rows).items() if n > 1]
    if dupes:
        raise SystemExit(
            f"{len(dupes)} ENSPs appear more than once -- loading would overwrite rows. "
            f"First few: {dupes[:5]}"
        )

    print(f"rows to load: {len(rows)}")
    for k, n in skipped.most_common():
        print(f"  skipped — {k}: {n}")

    # hgnc_id + ensembl_gene from D1, joined on the **reference's UniProt accession**.
    # Every one of these forms belongs to a gene whose selected reference is a UniProt
    # entry, so the accession is available and is a stable key. Joining on the gene symbol
    # instead is what filed MYO18A's topology row under TIAF1 once already.
    accs = sorted({(r["reference_uniprot_acc_full"] or "").split("-")[0] for r in rows} - {""})
    lookup: dict[str, dict] = {}
    for i in range(0, len(accs), 300):
        inlist = ",".join(sql_str(a) for a in accs[i : i + 300])
        for rec in d1(
            f"SELECT uniprot_acc, hgnc_id, ensembl_gene, hgnc_symbol "
            f"FROM gene_identifier_public WHERE uniprot_acc IN ({inlist});"
        ):
            lookup[rec["uniprot_acc"]] = rec
    hit = 0
    for r in rows:
        rec = lookup.get((r["reference_uniprot_acc_full"] or "").split("-")[0])
        if rec:
            r["hgnc_id"] = rec.get("hgnc_id") or ""
            r["ensembl_gene_id"] = rec.get("ensembl_gene") or ""
            hit += 1
    print(
        f"  resolved hgnc_id for {hit}/{len(rows)} rows "
        f"({100 * hit / max(len(rows), 1):.0f}%) via the reference accession"
    )

    cols = list(rows[0])
    values = []
    for r in rows:
        values.append("(" + ",".join(
            str(r[c]) if isinstance(r[c], int) else sql_str(r[c]) for c in cols
        ) + ")")

    # Pack values into statements under the size limit.
    batches, cur, size = [], [], 0
    for v in values:
        if cur and size + len(v) > MAX_STATEMENT_BYTES:
            batches.append(cur)
            cur, size = [], 0
        cur.append(v)
        size += len(v) + 1
    if cur:
        batches.append(cur)

    if not args.apply:
        print(f"\nDRY RUN — {len(values)} rows in {len(batches)} statements "
              f"(max {MAX_STATEMENT_BYTES:,} bytes each).")
        print("first row:\n  " + values[0][:400])
        print(f"\ntopology_version = {TOPOLOGY_VERSION}")
        return 0

    done = 0
    for batch in batches:
        stmt = (
            f"INSERT OR REPLACE INTO topology_ensembl_public ({','.join(cols)}) "
            f"VALUES {','.join(batch)};"
        )
        with tempfile.NamedTemporaryFile("w", suffix=".sql", delete=False) as fh:
            fh.write(stmt)
            path = fh.name
        proc = subprocess.run(
            ["npx", "--yes", "wrangler", "d1", "execute", DATABASE, "--remote", f"--file={path}"],
            capture_output=True, text=True, cwd=WRANGLER_CWD,
        )
        Path(path).unlink(missing_ok=True)
        if proc.returncode != 0:
            sys.stderr.write(proc.stderr[-1500:])
            return proc.returncode
        done += len(batch)
        print(f"  loaded {done}/{len(values)}")
    print(f"done — topology_version {TOPOLOGY_VERSION}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
