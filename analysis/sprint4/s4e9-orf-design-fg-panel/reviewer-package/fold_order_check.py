"""Independent check that every L1 fold is the protein in the order of record.

l1-tag-refold proved its folds against the draft integration build (#299) and the #303 recodes
(``tag_refold_construct_check.tsv``). This check starts from the other end, the file Twist
receives: the dated, recoded order-of-record FASTA. It shares no code with the refold session.

For every insert in the FASTA it

1. translates from the wrapper's Kozak (GCCGCCACC) through the first stop, with its own
   standard codon table;
2. looks up a tagged fold whose ``tagged_sha256`` equals the protein, or else the protein minus
   an N-terminal peptide (a leader folded as the cleaved mature chain), and records that
   peptide so a reader can see it is a leader and nothing else;
3. downloads that fold's PDB, checks the file's SHA-256 against ``pdb_sha256``, and reads the
   residues actually present in the model (CA atoms, chain A), which must equal the folded
   sequence exactly.

Step 3 is what makes the check independent of the refold tables: the residues come from the
coordinates, not from a hash someone wrote down.

Output: ``fold_order_check.tsv`` (one row per ordered insert) and
``fold_order_check_summary.json``. Requested by the coordinator (2026-11-13, "the
fold-equals-ordered-insert proof"), under Becca's overnight delegation of 2026-11-12;
implemented by the Claude Code session l1-reviewer-package.

    uv run --locked python analysis/sprint4/s4e9-orf-design-fg-panel/reviewer-package/fold_order_check.py
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import re
import subprocess
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
EXP = HERE.parent
REPO = EXP.parents[2]
EXP_REL = "analysis/sprint4/s4e9-orf-design-fg-panel"

#: The order of record (synthesis_compliance.FASTA_KEY), pinned by SHA-256.
FASTA_KEY = (
    "sprint4/s4e9-orf-design-fg-panel/synthesis/2026-11-15-final-e/l1_order_inserts_recoded.fasta"
)
FASTA_SHA256 = "20e8f275"  # prefix; superseded: 2026-11-13/ (d1257f0a), -one-arm/ (1f29936b), -one-arm-renamed/ (b426e433)
#: l1-tag-refold's structure table (PR #301) at the commit checked.
REFOLD_COMMIT = "541c6c18"  # main at the final regeneration (2026-11-14); was 3a7dfbd9 (#301)
STRUCTURES = f"{EXP_REL}/tag_refold_structures.tsv"
#: A mature chain may drop at most this many N-terminal residues (leader + SP); longer
#: differences are reported as unmatched, never silently accepted.
MAX_LEADER = 60
WORKERS = 2  # Becca's cap on local parallelism after the 2026-11-13 quiet period

OUT = EXP / "fold_order_check.tsv"
SUMMARY = EXP / "fold_order_check_summary.json"

#: An arm id, the first header token of a per-insert order (``l1-pNNNN-aK``).
ARM_ID = re.compile(r"l[12]-p[0-9]{4}-a[1-9][0-9]*")
BASES = "TCAG"
AMINO = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON = {a + b + c: AMINO[16 * i + 4 * j + k]
         for i, a in enumerate(BASES) for j, b in enumerate(BASES) for k, c in enumerate(BASES)}  # fmt: skip
THREE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
         "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
         "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}  # fmt: skip


def sha(s: str | bytes) -> str:
    return hashlib.sha256(s.encode() if isinstance(s, str) else s).hexdigest()


def translate(nt: str) -> str:
    start = nt.find("GCCGCCACC")
    if start < 0 or nt[start + 9 : start + 12] != "ATG":
        raise ValueError("no GCCGCCACC-ATG start")
    prot = []
    for i in range(start + 9, len(nt) - 2, 3):
        aa = CODON[nt[i : i + 3]]
        if aa == "*":
            return "".join(prot)
        prot.append(aa)
    raise ValueError("no stop codon")


def pdb_sequence(data: bytes) -> str:
    if data[:2] == b"\x1f\x8b":
        data = gzip.decompress(data)
    seq, seen = [], set()
    for line in data.decode().splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA" and line[21] in ("A", " "):
            key = (line[21], line[22:27])
            if key not in seen:
                seen.add(key)
                seq.append(THREE.get(line[17:20], "X"))
    return "".join(seq)


def main() -> int:
    from deliverome_analysis.s3 import processed

    fasta_bytes = Path(processed.download(FASTA_KEY)).read_bytes()
    fasta_sha = sha(fasta_bytes)
    assert fasta_sha.startswith(FASTA_SHA256), f"order of record changed: {fasta_sha}"
    inserts, name = [], None
    for line in fasta_bytes.decode().splitlines():
        if line.startswith(">"):
            name = line[1:]
            inserts.append([name, ""])
        elif line:
            inserts[-1][1] += line.strip().upper()

    table = subprocess.run(["git", "show", f"{REFOLD_COMMIT}:{STRUCTURES}"], cwd=REPO,
                           check=True, capture_output=True, text=True).stdout  # fmt: skip
    folds: dict[str, dict[str, str]] = {}
    for r in csv.DictReader(table.splitlines(), delimiter="\t"):
        if r["kind"] == "tagged" and r.get("pdb_key"):
            folds.setdefault(r["tagged_sha256"], r)

    rows = []
    for header, nt in inserts:
        tok = header.split()
        if ARM_ID.fullmatch(tok[0]):
            # Per-insert order (from the first per-insert cut): arm id, insert id, design, rung.
            arm, insert, design, rung = tok[:4]
        else:  # an order of record cut before then: member-named construct, design, rung
            arm, insert = "", tok[0].split("|")[0].split(":")[0]
            design = tok[1] if len(tok) > 1 else ""
            rung = tok[2] if len(tok) > 2 else ""
        prot = translate(nt)
        fold, dropped = None, ""
        for k in range(0, min(MAX_LEADER, len(prot) - 1) + 1):
            if (fold := folds.get(sha(prot[k:]))) is not None:
                dropped = prot[:k]
                break
        rows.append({"insert": insert, "arm_id": arm, "design_id": design, "rung": rung, "protein_aa": len(prot), "protein_sha256": sha(prot),
                     "fold_id": fold["fold_id"] if fold else "",
                     "folded_form": ("" if not fold else "precursor" if not dropped else "mature chain"),
                     "n_terminal_peptide_not_folded": dropped,
                     "pdb_key": fold["pdb_key"] if fold else "",
                     "pdb_sha256_expected": fold["pdb_sha256"] if fold else "",
                     "_seq": prot[len(dropped):]})  # fmt: skip

    def check(row: dict) -> dict:
        if not row["pdb_key"]:
            row.update(pdb_sha256_ok="", model_residues_equal="", status="no fold")
            return row
        key = row["pdb_key"].removeprefix("processed://")
        data = Path(processed.download(key)).read_bytes()
        row["pdb_sha256_ok"] = "yes" if sha(data) == row["pdb_sha256_expected"] else "no"
        row["model_residues_equal"] = "yes" if pdb_sequence(data) == row["_seq"] else "no"
        ok = row["pdb_sha256_ok"] == "yes" and row["model_residues_equal"] == "yes"
        row["status"] = "ok" if ok else "MISMATCH"
        return row

    with ThreadPoolExecutor(WORKERS) as pool:
        rows = list(pool.map(check, rows))

    cols = [c for c in rows[0] if not c.startswith("_")]
    with OUT.open("w", newline="") as fh:
        w = csv.DictWriter(fh, cols, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    summary = {
        "order_of_record": f"processed://{FASTA_KEY}",
        "order_of_record_sha256": fasta_sha,
        "refold_table": f"{STRUCTURES} at {REFOLD_COMMIT} (PR #301)",
        "inserts": len(rows),
        "status": dict(Counter(r["status"] for r in rows)),
        "folded_form": dict(Counter(r["folded_form"] for r in rows if r["fold_id"])),
        "n_terminal_peptides_not_folded": dict(
            Counter(
                r["n_terminal_peptide_not_folded"]
                for r in rows
                if r["n_terminal_peptide_not_folded"]
            )
        ),
        "by_rung": {k: dict(v) for k, v in _by_rung(rows).items()},
        "decided_by": "coordinator under Becca's 2026-11-12 overnight delegation",
        "implemented_by": "Claude Code session l1-reviewer-package, 2026-11-13",
    }
    SUMMARY.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({k: summary[k] for k in ("inserts", "status", "folded_form")}, indent=2))
    return 0 if summary["status"].get("ok") == len(rows) else 1


def _by_rung(rows: list[dict]) -> dict[str, Counter]:
    out: dict[str, Counter] = {}
    for r in rows:
        out.setdefault(r["rung"], Counter())[r["status"]] += 1
    return out


if __name__ == "__main__":
    raise SystemExit(main())
