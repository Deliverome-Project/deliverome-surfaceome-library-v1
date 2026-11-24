"""The membrane-side gate's verdict at every internal tag site that is actually ordered.

The designs of record are S5E13's final internal decisions (``l1_internal_decisions.tsv``:
the v5 ranker under the gate's mask, with a ligand-contact veto). This script reads them at
a pinned commit, checks the file's SHA-256, and gives each ordered site the gate's verdict on
the member's own model, at the pins of record on ``main``:

* ``internal_loop`` (and ``internal_loop_orientation_unconfirmed``): the junction's verdict,
  computed by the gate's own ``score_picks`` (thresholds from ``calibration.json``, the manual
  partner-interface exclusion, the topology overrides).
* ``snorkel``: a snorkel adds a helix and the tag at the C-terminus, so it is displayed only
  if the C-terminus is cytoplasmic. The verdict is the side of the C-terminal flank on the
  model (``snorkel_candidates.c_terminal_side``): ``pass_snorkel`` when it lies below the
  mid-plane, ``fail_snorkel_c_terminus_outside`` when above, ``unscorable`` when unmodelled.

Members are joined on ``sequence_sha256`` of the designed protein, never on member id.
Read-only on the shipped gate tables: writes ``gate_on_ordered_sites.tsv`` and
``gate_on_ordered_sites_summary.json`` into a new dated folder.

    uv run --locked python gate_on_ordered_sites.py --out-dir ordered-sites-2026-11-14
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import gate_inputs as gi  # noqa: E402
import membrane_side_gate as g  # noqa: E402
from placements import member_placements  # noqa: E402
from snorkel_candidates import c_terminal_side  # noqa: E402

#: The final internal decisions, as merged on main (S5E13, PR #245).
DECISIONS_COMMIT = "debbc418"
DECISIONS_PATH = "library/l1/tag-site-benchmarking/l1_internal_decisions.tsv"
DECISIONS_SHA256 = "db1dbf7da96f8aa17d47cb6340552f58e455c0013cf8ec83684df8cb6c9eaa5a"
#: main's pins of record on 2026-11-14 (unchanged since 2cc9432d).
MAIN_PINS = dict(
    panel_commit="34910ec68217eff67e10a8299524ae93bb210b8a",
    cds_key="sprint4/s4e9-orf-design-fg-panel/library-candidates/"
    "2026-11-12-ortholog-start-trims-2/candidate_cds.json",
    snapshot_key="accessible-surfaceome/topology/2026-11-12-cyno-restarts/topology.parquet",
    snapshot_sha256="c8fb5874437acd0677da76e2864aa76089ac7fa072f3defd2338f15f3890292e",
)
LOOP_DECISIONS = ("internal_loop", "internal_loop_orientation_unconfirmed")


def load_decisions(commit: str = DECISIONS_COMMIT, sha: str = DECISIONS_SHA256) -> pd.DataFrame:
    from deliverome_analysis.provenance import read_pinned_bytes

    raw = read_pinned_bytes(commit, DECISIONS_PATH, repo=gi.REPO)
    digest = hashlib.sha256(raw).hexdigest()
    if digest != sha:
        raise SystemExit(f"decisions SHA-256 {digest} differs from the pin")
    return pd.read_csv(io.BytesIO(raw), sep="\t")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--decisions-commit", default=DECISIONS_COMMIT)
    ap.add_argument("--decisions-sha256", default=DECISIONS_SHA256)
    args = ap.parse_args(argv)
    out = HERE / args.out_dir if not args.out_dir.is_absolute() else args.out_dir
    out.mkdir(parents=True, exist_ok=True)

    pins = gi.configure(**MAIN_PINS)
    dec = load_decisions(args.decisions_commit, args.decisions_sha256)
    ordered = dec[dec["decision"].isin((*LOOP_DECISIONS, "snorkel"))].copy()

    tp = gi.tag_paths()
    by_sha = tp.set_index("sequence_sha256")["candidate_id"]
    extra = {by_sha[s] for s in ordered["sequence_sha256"] if s in by_sha.index}
    members = gi.internal_members(extra)
    placed = member_placements(members)
    g.V2_AGREES.update(dict(zip(members["candidate_id"], members["v2_termini_agree"], strict=True)))
    g.CID_OF_DESIGN.update(dict(zip(members["design_id"], members["candidate_id"], strict=True)))
    th = g.thresholds()
    cid_of_sha = dict(zip(members["sequence_sha256"], members["candidate_id"], strict=True))
    mem = members.set_index("candidate_id")

    rows = []
    loops = []
    for r in ordered.itertuples():
        cid = cid_of_sha.get(r.sequence_sha256)
        base = dict(
            design_id=r.design_id,
            decisions_member_id=r.member_id,
            candidate_id=cid,
            gene_symbol=r.gene,
            candidate_set=r.candidate_set,
            decision=r.decision,
            insert_after_residue=r.insert_after_residue,
            sequence_sha256=r.sequence_sha256,
        )
        if cid is None:
            rows.append(
                base
                | dict(verdict="not_evaluated", reason="protein not in the panel at main's pins")
            )
            continue
        if r.decision == "snorkel":
            p = placed[cid]
            m = mem.loc[cid]
            if p["status"] != "placed":
                rows.append(base | dict(verdict="unscorable", reason=p["status"]))
                continue
            cz = c_terminal_side(m.topology, p["z"], p["plddt"])
            v = (
                "unscorable"
                if np.isnan(cz)
                else ("pass_snorkel" if cz < 0 else "fail_snorkel_c_terminus_outside")
            )
            rows.append(
                base
                | dict(
                    verdict=v,
                    height_A=round(cz, 1) if not np.isnan(cz) else np.nan,
                    topology_tool=m.topology_layer,
                    reason=f"C-terminal flank (15 residues after the last helix) at {cz:+.1f} A",
                )
            )
            continue
        loops.append(
            dict(
                candidate_id=cid,
                source="ordered_internal",
                site_id=f"{r.design_id}-ordered",
                rank=r.pick_rank,
                shipped=1,
                insert_after_residue=int(r.insert_after_residue),
                numbering_check="v1_O_check",
            )
        )
        rows.append(base | dict(_loop_key=f"{r.design_id}-ordered"))

    scored = g.score_picks(pd.DataFrame(loops), members, placed, th).set_index("site_id")
    final = []
    for row in rows:
        key = row.pop("_loop_key", None)
        if key is not None:
            s = scored.loc[key]
            row |= dict(
                verdict=s["verdict"],
                membrane_verdict=s.get("membrane_verdict", s["verdict"]),
                height_A=s.get("height_A"),
                plddt=s.get("plddt"),
                v1_state=s.get("v1_state"),
                residue_before=s.get("residue_before"),
                residue_after=s.get("residue_after"),
                numbering_check=s.get("numbering_check"),
                topology_tool=s.get("topology_tool"),
                reason=s.get("reason"),
            )
        final.append(row)
    table = pd.DataFrame(final).sort_values(["verdict", "design_id"])
    table["implemented_by"] = "l1-membrane-side-gate session (Claude Code)"
    table.to_csv(out / "gate_on_ordered_sites.tsv", sep="\t", index=False)

    flagged = table[~table["verdict"].isin(["pass", "pass_snorkel"])]
    summary = dict(
        pins=pins,
        decisions=dict(
            commit=args.decisions_commit, path=DECISIONS_PATH, sha256=args.decisions_sha256
        ),
        thresholds=th.__dict__,
        ordered_sites=len(table),
        by_decision=table["decision"].value_counts().to_dict(),
        by_verdict=table["verdict"].value_counts().to_dict(),
        by_decision_and_verdict={
            f"{d}|{v}": int(n)
            for (d, v), n in table.groupby(["decision", "verdict"]).size().items()
        },
        not_passing=flagged[
            ["design_id", "gene_symbol", "decision", "insert_after_residue", "verdict", "height_A"]
        ].to_dict("records"),
    )
    (out / "gate_on_ordered_sites_summary.json").write_text(
        json.dumps(summary, indent=1, default=str) + "\n"
    )
    print(
        json.dumps({k: summary[k] for k in ("ordered_sites", "by_decision_and_verdict")}, indent=1)
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
