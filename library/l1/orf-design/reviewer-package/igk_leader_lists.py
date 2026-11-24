"""Gene lists for the Igκ-leader table of the brief ("Signal peptides and the Igκ leader").

Becca, on the 2026-11-22 brief: "we will also need lists of genes here for reviewers to
evaluate". Each count in that table is one file in ``igk_lists/``, written here from the decision
tables the builds read, so a list cannot drift from the construct it describes:

* N after SP: ``nterm/nterm_decisions.tsv``. ``arm`` is ``igk_swap`` (the leader replaces the
  native signal peptide) or ``native_site`` (the CDS keeps its own). An Igκ swap is split by
  what triggered it: a weak predicted cleavage site (``weak_sp_flag``), a signal peptide only
  DeepTMHMM v1 calls (``sp_review_verdict``, one file per verdict), or neither (no candidate
  site keeps the native cleavage with the tag in place).
* Residue 1, multipass: ``constructs_l1_multipass.tsv``, ``design`` ``igk_tag`` or ``tag_alone``.
* Residue 1, type III: ``constructs_l1_type3.tsv``, the same, with its ``subgroup``.

Each row carries the design id, gene, candidate set, species and the column that put it in the
group, verbatim. ``igk_lists/README.md`` gives each file's count. Rerun after any rung build.

Implemented by Claude Code session l1-review-export-scope (fork), 2026-11-22.

    uv run --locked python igk_leader_lists.py
"""

from __future__ import annotations

import csv
import shutil
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
S4E9 = HERE.parent
OUT = HERE / "igk_lists"
FIELDS = ["design_id", "gene_symbol", "candidate_set", "species", "group", "basis"]


def _tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def nterm_group(r: dict[str, str]) -> tuple[str, str] | None:
    if r["arm"] == "native_site":
        return "nterm_native", r.get("reason", "")
    if r["arm"] != "igk_swap":
        return None
    if (r.get("weak_sp_flag") or "").strip():
        return "nterm_igk_weak_signal_peptide", r["weak_sp_flag"]
    verdict = (r.get("sp_review_verdict") or "").strip()
    if verdict:
        return f"nterm_igk_deeptmhmm_only_{verdict}", r.get("sp_sources_checked", "")
    return "nterm_igk_native_cleavage_lost", r.get("reason", "")


def groups() -> dict[str, list[dict[str, str]]]:
    out: dict[str, list[dict[str, str]]] = defaultdict(list)

    def add(group: str, r: dict[str, str], basis: str) -> None:
        out[group].append(
            {
                "design_id": r["design_id"],
                "gene_symbol": r["gene_symbol"],
                "candidate_set": r.get("candidate_set", ""),
                "species": r.get("species", ""),
                "group": group,
                "basis": " ".join(basis.split()),
            }
        )

    for r in _tsv(HERE / "nterm" / "nterm_decisions.tsv"):
        g = nterm_group(r)
        if g:
            add(g[0], r, g[1])
    for r in _tsv(S4E9 / "constructs_l1_multipass.tsv"):
        add(f"multipass_{r['design']}", r, r.get("reason", ""))
    for r in _tsv(S4E9 / "constructs_l1_type3.tsv"):
        add(f"type3_{r['design']}", r, f"{r.get('subgroup', '')}: {r.get('reason', '')}")
    return out


def main() -> int:
    data = groups()
    if OUT.exists():
        shutil.rmtree(OUT)
    OUT.mkdir()
    lines = [
        "# Igκ-leader groups",
        "",
        "One file per group of the brief's Igκ-leader table, written by",
        "[`igk_leader_lists.py`](../igk_leader_lists.py) from the decision tables.",
        "",
        "| File | CDSs |",
        "|---|---:|",
    ]
    for g in sorted(data):
        rows = sorted(data[g], key=lambda r: (r["gene_symbol"], r["design_id"]))
        with (OUT / f"{g}.tsv").open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(rows)
        lines.append(f"| [`{g}.tsv`]({g}.tsv) | {len(rows):,} |")
        print(f"{g}: {len(rows)}")
    (OUT / "README.md").write_text("\n".join(lines) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
