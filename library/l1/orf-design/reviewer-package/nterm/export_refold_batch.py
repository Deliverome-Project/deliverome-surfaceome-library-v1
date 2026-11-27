"""Tagged N-after-SP sequences for l1-tag-refold to fold, one FASTA and one table per batch.

Batch ``nterm_held_2026_11_12`` holds the two held groups whose refold is the evidence Becca
needs:
- the 354 members whose only preserved sites lie inside the predicted fold, each at its
  re-prediction-gated site;
- the 53 members whose cleavage shifts upstream of the cassette, at the site proposed in
  ``nterm_site_picks.tsv``.

Batch ``nterm_igk_2026_11_12`` holds the Igκ-swapped tagged forms, one per member (written by
``--igk`` once the Igκ calls exist).

Every tagged sequence is ``native[:n] + p003 tandem + native[n:]``, or the Igκ form, and is
named by its SHA-256. The native protein's SHA-256 is given too, so the folder can reuse the
native MSA, with gap columns for the 51 cassette residues starting after residue ``n``.

**Attribution** (convention of 2026-11-10). Sites come from ``nterm_site_picks.tsv``:
re-prediction-gated (Becca, 2026-11-11), or proposed for the 53 by the Claude Code session
``nterm-reviewer-package`` (PR #254). The batch is assembled by that session for the
2026-11-12 overnight refold (coordinator plan ``OVERNIGHT-L1-2026-11-12.md``).

    uv run --locked python export_refold_batch.py
"""

from __future__ import annotations

import csv
import hashlib
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
S4E9 = HERE.parents[1]
sys.path.insert(0, str(S4E9))
sys.path.insert(0, str(HERE))

import build_constructs as bc  # noqa: E402
import nterm_placement as npl  # noqa: E402

OUT = HERE / "refold_batches"
GROUPS = ("held_inside_fold", "held_moved_upstream")


def sha(s: str) -> str:
    return hashlib.sha256(s.encode()).hexdigest()


def write_batch(name: str, rows: list[dict]) -> None:
    OUT.mkdir(exist_ok=True)
    cols = [k for k in rows[0] if k != "tagged_sequence"]
    with (OUT / f"{name}.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(
            fh, fieldnames=cols, delimiter="\t", lineterminator="\n", extrasaction="ignore"
        )
        w.writeheader()
        w.writerows(rows)
    with (OUT / f"{name}.fasta").open("w") as fh:
        for r in rows:
            fh.write(
                f">{r['tagged_sha256']} {r['candidate_id']} {r['gene_symbol']} "
                f"{r['arm']} after{r['insert_after_residue']}\n{r['tagged_sequence']}\n"
            )
    print(f"{name}: {len(rows)} sequences")


def main() -> int:
    from deliverome_analysis.core.s3 import processed

    members, _ = npl.load_members(processed)
    with (HERE / "nterm_site_picks.tsv").open() as fh:
        picks = [r for r in csv.DictReader(fh, delimiter="\t") if r["status"] in GROUPS]
    rows = []
    for g in GROUPS:
        for p in (r for r in picks if r["status"] == g):
            native = members[p["candidate_id"]]["sequence"]
            n = int(p["insert_after_residue"])
            tagged = native[:n] + bc.CASSETTE_AA + native[n:]
            assert sha(tagged) == p["tagged_protein_sha256"], p["candidate_id"]
            rows.append(
                {
                    "order": len(rows) + 1,
                    "candidate_id": p["candidate_id"],
                    "design_id": p["design_id"],
                    "gene_symbol": p["gene_symbol"],
                    "group": g,
                    "arm": "native_sp_tagged",
                    "insert_after_residue": n,
                    "cassette_residues": f"{n + 1}-{n + len(bc.CASSETTE_AA)}",
                    "cleavage_site": p["cleavage_site"],
                    "inside_predicted_fold": p["inside_predicted_fold"],
                    "pick_rule": p["pick_rule"],
                    "native_length": len(native),
                    "native_sha256": sha(native),
                    "tagged_length": len(tagged),
                    "tagged_sha256": sha(tagged),
                    "tagged_sequence": tagged,
                }
            )
    write_batch("nterm_held_2026_11_12", rows)
    if "--batch2" in sys.argv:
        write_batch("nterm_batch2_2026_11_12", batch2(members))
    return 0


def batch2(members: dict) -> list[dict]:
    """One proposed tagged form per member from ``overnight_proposals.tsv`` (Igκ, ortholog, new).

    Each sequence is rebuilt from the member's designed protein and must hash to the form
    SignalP scored. For Igκ forms the leader plus Asp replaces residues ``1..c``.
    """
    from build_gpcr_rescue_arms import LEADER_AA, PLUS1_AA

    with (HERE / "overnight_proposals.tsv").open() as fh:
        props = [r for r in csv.DictReader(fh, delimiter="\t") if r["proposal"] != "hold"]
    with (HERE / "overnight_inputs.tsv").open() as fh:
        parent_c = {
            (r["candidate_id"], r["arm"]): int(r["cleavage_site"])
            for r in csv.DictReader(fh, delimiter="\t")
            if r["arm"] in ("native", "igk_native")
        }
    order = {"new_member": 0, "igk_failing_native": 1, "igk_disagree": 2, "ortholog": 3}
    rows = []
    for p in sorted(props, key=lambda r: (order[r["group"]], r["candidate_id"])):
        native = members[p["candidate_id"]]["sequence"]
        n = int(p["insert_after_residue"])
        if p["proposal"] == "igk_swap":
            c = parent_c[(p["candidate_id"], "igk_native")]
            parent = LEADER_AA + PLUS1_AA + native[c:]
        else:
            parent = native
        tagged = parent[:n] + bc.CASSETTE_AA + parent[n:]
        assert sha(tagged) == p["tagged_sha256"], p["candidate_id"]
        rows.append(
            {
                "order": len(rows) + 1,
                "candidate_id": p["candidate_id"],
                "design_id": p["design_id"],
                "gene_symbol": p["gene_symbol"],
                "group": p["group"],
                "arm": "igk_tagged" if p["proposal"] == "igk_swap" else "native_sp_tagged",
                "insert_after_residue": n,
                "cassette_residues": f"{n + 1}-{n + len(bc.CASSETTE_AA)}",
                "cleavage_site": p["parent_cleavage"],
                "inside_predicted_fold": p["inside_predicted_fold"],
                "pick_rule": "reprediction_gated (proposed)",
                "native_length": len(native),
                "native_sha256": sha(native),
                "parent_sha256": sha(parent),
                "tagged_length": len(tagged),
                "tagged_sha256": sha(tagged),
                "tagged_sequence": tagged,
            }
        )
    return rows


if __name__ == "__main__":
    raise SystemExit(main())
