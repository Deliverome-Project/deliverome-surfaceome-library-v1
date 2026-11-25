"""Index the v8_prodomino refolds' verdicts: source rows, then verdict and viewer rows.

``coreceptor_additions_publish.py`` uploaded the v8 batches' models and appended their metrics and
structures (``…/tagged/2026-11-22-v8-internal`` and ``…-v8-snorkel``), but it writes a verdict only
for a tagged protein that already has a row in ``tag_refold_sources.tsv``, and the v8 arms had
none. This adds one source row per built or held v8 arm (keyed by its permanent arm id in S5E13
``l1_internal_arm_decisions.tsv``, matched on ``refold_tagged_sha256``), then appends the verdict
and viewer rows for them from the scored metrics, with the same columns and rules as the publish
script. Every existing row is kept byte for byte (asserted by ``append``). A snorkel whose exact
protein was folded before reuses that fold's metrics row.

Decided by Becca Carlson, 2026-11-22 ("i want to stay with the v8 picks exactly"); implemented by
Claude Code session l1-review-export-scope, 2026-11-22.

    uv run --locked python library/l1/orf-design/tag_refold/v8_sources_verdicts.py
"""

from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
EXP = HERE.parent
REPO = EXP.parents[2]
sys.path.insert(0, str(HERE))

from coreceptor_additions_publish import append, tsv  # noqa: E402
from publish_tagged_structures import VERDICT_COLS  # noqa: E402

ARMS = REPO / "library/l1/tag-site-benchmarking" / "l1_internal_arm_decisions.tsv"
PREFIX = "sprint4/s4e9-orf-design-fg-panel/predicted-structures/tagged/2026-11-22-v8-"
BATCHES = {"internal": "internal", "snorkel": "cterm"}
DECIDED = 'Becca Carlson, 2026-11-22 ("i want to stay with the v8 picks exactly")'


def main() -> int:
    arms = {r["refold_tagged_sha256"]: r for r in tsv(ARMS)[1] if r["refold_tagged_sha256"]}
    metrics: dict[str, dict] = {}
    for r in tsv(EXP / "tag_refold_metrics.tsv")[1]:
        metrics[r["tagged_sha256"]] = r  # the latest row for a protein wins
    structures = {}
    for r in tsv(EXP / "tag_refold_structures.tsv")[1]:
        if r["pdb_key"].startswith("processed://" + PREFIX):
            structures[r["tagged_sha256"]] = r
    have_src = {
        (r["tagged_sha256"], r["source_id"]) for r in tsv(EXP / "tag_refold_sources.tsv")[1]
    }
    have_v = {(r["tagged_sha256"], r["source_id"]) for r in tsv(EXP / "tag_refold_verdicts.tsv")[1]}
    sources, verdicts, viewer = [], [], []
    for batch, rung in BATCHES.items():
        queue = tsv(EXP / f"tag_refold_queue_2026_11_22_v8_{batch}.tsv")[1]
        for q in queue:
            if q["kind"] != "tagged":
                continue
            h = q["tagged_sha256"]
            arm = arms.get(h)
            if arm is None:  # a v8 site that the delta did not keep
                continue
            src = {
                "tagged_sha256": h, "rung": rung, "source_id": arm["arm_id"],
                "design_id": arm["design_id"], "gene_symbol": arm["gene"], "tier": q["tier"],
                "reason": f"v8_prodomino {arm['arm_kind']} site (2026-11-22)",
                "arm": arm["arm_id"], "design": "v8", "rank": arm.get("design_rank", ""),
                "loop_member": "", "review_flags": "",
                "source": f"designed protein + cassette (v8 {batch} batch, 2026-11-22)",
            }  # fmt: skip
            if (h, arm["arm_id"]) not in have_src:
                sources.append(src)
            m, s = metrics.get(h), structures.get(h)
            assert m is not None and s is not None, f"{arm['arm_id']}: no metrics or structure"
            if (h, arm["arm_id"]) in have_v:
                continue
            v = {**src, "pdb_key": s["pdb_key"], "decided_by": DECIDED}
            v.update({c: m.get(c, "") for c in VERDICT_COLS})
            verdicts.append(v)
            n = q["native_sha256"]
            viewer.append({**v, "native_sha256": n, "ordered_protein_sha256": "",
                           "native_model_path": "",
                           "msa_key": f"processed://{PREFIX}{batch}/msa/{n}.a3m.gz"})  # fmt: skip
    print("sources", append(EXP / "tag_refold_sources.tsv", sources))
    print("verdicts", append(EXP / "tag_refold_verdicts.tsv", verdicts))
    print("viewer", append(EXP / "tag_refold_viewer_manifest.tsv", viewer))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
