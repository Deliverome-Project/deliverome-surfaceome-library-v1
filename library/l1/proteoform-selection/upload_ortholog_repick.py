"""Load DeepTMHMM v1 topology for the re-picked orthologs, without disturbing the live rows.

The current ortholog picks stay exactly as they are in `topology_public` — the main
surfaceome database continues to serve them. These rows land under their **own
topology_version**, which is part of the primary key, so the two coexist and the re-pick can
be compared against the live cohort or abandoned without a restore.

Rows route by how the picked proteoform is identified, which is also what decides the table:

* a **UniProt accession** → `topology_public`, the existing table, new version;
* an **Ensembl protein id** → `topology_ensembl_public`, since `topology_public` requires a
  non-null UniProt accession in its key.

Species prefixes differ — mouse is `ENSMUSP`, cynomolgus `ENSMFAP`, only human is the bare
`ENSP` — so routing tests the shape rather than a literal prefix. Checking for `ENSP` alone
reported all 843 as UniProt-keyed when 442 are not.

    uv run --locked python upload_ortholog_repick.py            # dry run
    uv run --locked python upload_ortholog_repick.py --apply
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

DEFAULT_VERSION = "topo_2026_10_09_ortholog_repick"
TOOL_VERSION = "deeptmhmm-1.0.24"
DATABASE = "surfaceome_public"
WRANGLER_CWD = Path("~/Git/accessible-surfaceome").expanduser()
MAX_STATEMENT_BYTES = 90_000

SURFACE_LABELS = {"TM", "SP+TM"}
ORIENT = {"O": "extracellular", "I": "cytoplasmic"}
SPECIES_NAME = {"mouse": "mouse", "cyno": "cynomolgus"}
# Ensembl protein ids are ENS + an optional species code + P + digits; UniProt accessions
# never match this shape.
ENSEMBL_PROTEIN = re.compile(r"^ENS[A-Z]*P\d+$")


def parse_3line(path: Path) -> dict[str, dict]:
    out: dict[str, dict] = {}
    lines = path.read_text().splitlines()
    for i in range(len(lines) - 2):
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


def execute(statements: list[str], apply: bool) -> int:
    if not apply:
        return 0
    for stmt in statements:
        with tempfile.NamedTemporaryFile("w", suffix=".sql", delete=False) as fh:
            fh.write(stmt)
            path = fh.name
        proc = subprocess.run(
            ["npx", "--yes", "wrangler", "d1", "execute", DATABASE, "--remote", f"--file={path}"],
            capture_output=True, text=True, cwd=WRANGLER_CWD,
        )
        Path(path).unlink(missing_ok=True)
        if proc.returncode != 0:
            sys.stderr.write((proc.stderr or proc.stdout)[-1500:])
            return proc.returncode
    return 0


def batch(values: list[str], cols: list[str], table: str) -> list[str]:
    out, cur, size = [], [], 0
    for v in values:
        if cur and size + len(v) > MAX_STATEMENT_BYTES:
            out.append(
                f"INSERT OR REPLACE INTO {table} ({','.join(cols)}) VALUES {','.join(cur)};"
            )
            cur, size = [], 0
        cur.append(v)
        size += len(v) + 1
    if cur:
        out.append(
            f"INSERT OR REPLACE INTO {table} ({','.join(cols)}) VALUES {','.join(cur)};"
        )
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--apply", action="store_true")
    ap.add_argument("--three-line", default="ortholog_repick_changed.3line")
    ap.add_argument(
        "--version", default=DEFAULT_VERSION,
        help="topology_version to write under. Always a new value for a new run: the "
             "existing picks stay on their own version for the main database.",
    )
    ap.add_argument(
        "--changed-only", action="store_true",
        help="restrict to picks the re-pick moved; a gap-fill run needs every pick",
    )
    args = ap.parse_args()

    pred = parse_3line(Path(args.three_line))
    # Key by accession as well as by the composite id: a gap-fill run names its records by
    # the HUMAN symbol while the re-pick names them by the ortholog's, and the accession is
    # the field both agree on.
    by_acc = {k.rsplit("|", 1)[-1]: v for k, v in pred.items()}
    print(f"predictions: {len(pred)}")

    uniprot_rows, ensembl_rows, skipped = [], [], Counter()
    now = datetime.now(timezone.utc).isoformat(timespec="seconds")
    for species in ("mouse", "cyno"):
        for r in csv.DictReader(open(f"ortholog_repick_{species}.tsv"), delimiter="\t"):
            if args.changed_only and r["changed"] != "True":
                continue
            fid = f"{species}|{r['ortholog_symbol']}|{r['picked_acc']}"
            p = pred.get(fid) or by_acc.get(r["picked_acc"])
            if not p:
                skipped["no prediction"] += 1
                continue
            d = derive(p["topology"])
            common = {
                "topology_version": args.version,
                "species": SPECIES_NAME[species],
                "gene_symbol": r["ortholog_symbol"],
                "sequence": p["sequence"],
                "protein_length": len(p["sequence"]),
                "deeptmhmm_label": p["label"],
                **d,
                "per_residue_topology": p["topology"],
                "predicted_surface_membrane": 1 if p["label"] in SURFACE_LABELS else 0,
                "predicted_secreted": 1 if p["label"] == "SP" else 0,
                "tool_version": TOOL_VERSION,
                "retrieved_at": now,
            }
            acc = r["picked_acc"]
            if ENSEMBL_PROTEIN.match(acc):
                ensembl_rows.append(
                    {
                        **common,
                        "cohort": f"{species}_ortholog_repick",
                        "ensembl_protein_id": acc,
                        "ensembl_transcript_id": "",
                        "ensembl_gene_id": "",
                        "hgnc_id": "",
                        "reference_uniprot_acc_full": "",
                    }
                )
            else:
                base = acc.split("-")[0]
                uniprot_rows.append(
                    {
                        **common,
                        "cohort": f"{species}_ortholog_repick",
                        "hgnc_id": "",
                        "uniprot_acc": base,
                        "uniprot_acc_full": acc,
                        "isoform_id": acc,
                        "is_canonical": 0 if "-" in acc else 1,
                    }
                )
    print(f"  -> topology_public (UniProt keyed):        {len(uniprot_rows)}")
    print(f"  -> topology_ensembl_public (Ensembl keyed): {len(ensembl_rows)}")
    for k, n in skipped.most_common():
        print(f"  skipped — {k}: {n}")

    for label, rows in (("uniprot", uniprot_rows), ("ensembl", ensembl_rows)):
        key = (
            (r["topology_version"], r["cohort"], r.get("uniprot_acc_full") or r["ensembl_protein_id"])
            for r in rows
        )
        dupes = [k for k, n in Counter(key).items() if n > 1]
        if dupes:
            raise SystemExit(f"{label}: {len(dupes)} duplicate keys — would overwrite: {dupes[:3]}")

    statements = []
    if uniprot_rows:
        cols = list(uniprot_rows[0])
        vals = ["(" + ",".join(
            str(r[c]) if isinstance(r[c], int) else sql_str(r[c]) for c in cols
        ) + ")" for r in uniprot_rows]
        statements += batch(vals, cols, "topology_public")
    if ensembl_rows:
        cols = list(ensembl_rows[0])
        vals = ["(" + ",".join(
            str(r[c]) if isinstance(r[c], int) else sql_str(r[c]) for c in cols
        ) + ")" for r in ensembl_rows]
        statements += batch(vals, cols, "topology_ensembl_public")

    if not args.apply:
        print(f"\nDRY RUN — {len(statements)} statements, version {args.version}")
        print("  existing rows are untouched: the version is part of the primary key")
        return 0

    rc = execute(statements, True)
    if rc == 0:
        print(f"done — {args.version}")
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
