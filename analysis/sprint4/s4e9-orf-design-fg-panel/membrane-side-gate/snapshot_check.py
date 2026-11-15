"""Diff two membrane-side gate runs, member by member and site by site, keyed by design id.

Built for the 2026-11-14 check that the gate's calls do not depend on which DeepTMHMM
snapshot it read (the gate had pinned the 2026-11-06 snapshot; the snapshot of record became
2026-11-12). Read-only: it compares two output folders written by
``membrane_side_gate.py --out-dir`` and writes the differences to a third.

A *call* is, per member: its rung (`placement`), placement status, DeepTMHMM v1 row
(`topology_version`, cohort, layer), whether its extracellular face is supported, and its
per-residue and per-junction masks; per site: its gate verdict.

    uv run --locked python snapshot_check.py BASE_DIR NEW_DIR OUT_DIR
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

MEMBER_FIELDS = [
    "placement",
    "status",
    "topology_version",
    "topology_cohort",
    "topology_tool",
    "orientation_supported",
]
MASK_FIELDS = ["usable", "sequence_sha256", "mask", "junction_mask"]


def _site_key(df: pd.DataFrame) -> pd.Series:
    # API site ids carry the member id after "@"; member ids are renumbered between panels,
    # so key every site on the design id plus the id's member-independent part.
    return df["design_id"] + "|" + df["source"] + "|" + df["site_id"].str.split("@").str[0]


def compare(base: Path, new: Path) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    def members(d: Path) -> pd.DataFrame:
        m = pd.read_csv(d / "l1_internal_member_gate.tsv", sep="\t", keep_default_na=False)
        k = pd.read_csv(d / "l1_extracellular_mask_residues.tsv", sep="\t", keep_default_na=False)
        return m.merge(k[["design_id", *MASK_FIELDS]], on="design_id", how="left").set_index(
            "design_id"
        )

    mb, mn = members(base), members(new)
    rows = []
    for d in sorted(set(mb.index) | set(mn.index)):
        if d not in mn.index or d not in mb.index:
            rows.append(dict(design_id=d, change="only_in_" + ("base" if d in mb.index else "new")))
            continue
        b, n = mb.loc[d], mn.loc[d]
        changed = [f for f in MEMBER_FIELDS + MASK_FIELDS if str(b[f]) != str(n[f])]
        if changed:
            rows.append(
                dict(
                    design_id=d,
                    gene_symbol=n["gene_symbol"],
                    change=";".join(changed),
                    **{f"base_{f}": b[f] for f in changed if f not in ("mask", "junction_mask")},
                    **{f"new_{f}": n[f] for f in changed if f not in ("mask", "junction_mask")},
                )
            )
    member_diff = pd.DataFrame(rows)

    gb = pd.read_csv(base / "l1_internal_pick_gate.tsv", sep="\t")
    gn = pd.read_csv(new / "l1_internal_pick_gate.tsv", sep="\t")
    gb["key"], gn["key"] = _site_key(gb), _site_key(gn)
    j = gb[["key", "design_id", "gene_symbol", "shipped", "insert_after_residue", "verdict"]].merge(
        gn[["key", "insert_after_residue", "verdict"]],
        on="key",
        how="outer",
        suffixes=("_base", "_new"),
        indicator=True,
    )
    # Two missing values are equal here (an unscorable site has no residue on either side).
    same_res = (j["insert_after_residue_base"] == j["insert_after_residue_new"]) | (
        j["insert_after_residue_base"].isna() & j["insert_after_residue_new"].isna()
    )
    site_diff = j[
        (j["_merge"] != "both") | (j["verdict_base"] != j["verdict_new"]) | ~same_res
    ].copy()
    site_diff["_merge"] = site_diff["_merge"].astype(str)

    sb = json.loads((base / "gate_summary.json").read_text())
    sn = json.loads((new / "gate_summary.json").read_text())
    common = len(set(mb.index) & set(mn.index))
    summary = dict(
        base=str(base),
        new=str(new),
        base_pins={k: sb.get(k) for k in ("panel_commit", "cds_key", "topology_snapshot_key")},
        new_pins={k: sn.get(k) for k in ("panel_commit", "cds_key", "topology_snapshot_key")},
        members_base=len(mb),
        members_new=len(mn),
        members_in_both=common,
        members_with_any_change=int((member_diff["change"].str.len() > 0).sum())
        if len(member_diff)
        else 0,
        member_changes_by_field={
            f: int(member_diff["change"].str.contains(f).sum()) if len(member_diff) else 0
            for f in MEMBER_FIELDS + MASK_FIELDS + ["only_in_base", "only_in_new"]
        },
        sites_compared=int((j["_merge"] == "both").sum()),
        site_verdict_changes=int(
            ((j["_merge"] == "both") & (j["verdict_base"] != j["verdict_new"])).sum()
        ),
        shipped_site_verdict_changes=int(
            (
                (j["_merge"] == "both")
                & (j["verdict_base"] != j["verdict_new"])
                & (j["shipped"] == 1)
            ).sum()
        ),
        sites_only_in_one_run=int((j["_merge"] != "both").sum()),
        base_shipped_by_verdict=sb.get("shipped_by_verdict"),
        new_shipped_by_verdict=sn.get("shipped_by_verdict"),
    )
    return member_diff, site_diff, summary


def main(argv: list[str] | None = None) -> int:
    base, new, out = (Path(a) for a in (argv or sys.argv[1:]))
    out.mkdir(parents=True, exist_ok=True)
    member_diff, site_diff, summary = compare(base, new)
    member_diff.to_csv(out / "member_diff.tsv", sep="\t", index=False)
    site_diff.to_csv(out / "site_diff.tsv", sep="\t", index=False)
    (out / "diff_summary.json").write_text(json.dumps(summary, indent=1, default=str) + "\n")
    print(json.dumps(summary, indent=1, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
