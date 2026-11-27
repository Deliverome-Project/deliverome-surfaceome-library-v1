"""Which internal members are left with only tiny extracellular loops: candidates for a snorkel.

A snorkel (Brown et al. 2013, PLoS One 8:e73255, PMID 24023844) appends a transmembrane helix
and then the tag to a protein's cytoplasmic C-terminus, so the tag is displayed outside the
cell without being inserted into any loop of the protein. Becca asked on 2026-11-11 that
members with only tiny loops be considered for it seriously.

Two questions, both answered from the membrane-side gate's own placements:

1. **What is "tiny"?** For every published internal tag with surface display shown, in a
   multipass protein (the calibration set of ``calibrate.py``), the run of consecutive
   junctions that pass the gate around the published junction is measured on the protein's
   own AlphaFold DB model. The shortest such run is the smallest loop a working tag is
   known to sit in. A member whose longest run of passing junctions is shorter has no loop
   of a size any published tag has used.
2. **Would a snorkel work for the member?** Only if the C-terminus is cytoplasmic, so that
   one more helix carries the tag outside. Read on the model: the mean height of the
   confidently modelled residues in the 15 after the last DeepTMHMM v1 helix.

Writes ``l1_snorkel_candidates.tsv`` and ``snorkel_summary.json``. Proposes; changes nothing.

    uv run --locked python snorkel_candidates.py
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import gate_inputs as gi  # noqa: E402
from membrane_side_gate import junction, thresholds  # noqa: E402
from placements import member_placements, place_protein  # noqa: E402

from deliverome_analysis.surfaceome import membrane_side as ms  # noqa: E402
from deliverome_analysis.surfaceome import topology  # noqa: E402

GOLD = "scratch/becca/tag-site-corpus-audit/gold_standard.tsv"
#: Residues after the last v1 helix read for the C-terminal side.
C_FLANK = 15


def passing_junctions(top: str, z: np.ndarray, th) -> np.ndarray:
    """Boolean per junction n (index n-1): v1 O at n and junction height >= outside."""
    return np.array(
        [top[n - 1] == "O" and junction(z, n) >= th.outside for n in range(1, len(top) + 1)]
    )


def run_containing(ok: np.ndarray, n: int) -> int:
    """Length of the run of passing junctions that contains junction n (0 if n fails)."""
    i = n - 1
    if not ok[i]:
        return 0
    a = i
    while a > 0 and ok[a - 1]:
        a -= 1
    b = i
    while b + 1 < len(ok) and ok[b + 1]:
        b += 1
    return b - a + 1


def longest_run(ok: np.ndarray) -> int:
    return max(
        (len(m.group()) for m in re.finditer("1+", "".join("1" if x else "0" for x in ok))),
        default=0,
    )


def c_terminal_side(top: str, z: np.ndarray, plddt: np.ndarray) -> float:
    """Mean height of confident residues in the C_FLANK after the last v1 helix (NaN if none)."""
    last = top.rfind("M")
    if last < 0:
        return float("nan")
    idx = [
        i
        for i in range(last + 1, min(last + 1 + C_FLANK, len(top)))
        if plddt[i] >= ms.FIT_MIN_PLDDT and not np.isnan(z[i])
    ]
    return float(np.mean(z[idx])) if idx else float("nan")


def published_loop_runs(th) -> pd.DataFrame:
    from deliverome_analysis.core.s3 import processed

    gold = pd.read_csv(processed.download(GOLD), sep="\t")
    gold = gold[gold["surface_shown"] == "yes"]
    rows, cache = [], {}
    for s in gold.itertuples():
        acc = s.uniprot_acc
        if acc not in cache:
            path = gi.afdb_model(acc)
            cache[acc] = None
            if path is not None:
                model = ms.read_model(path)
                row = topology.for_sequence(model.sequence, tool=topology.TOOL_V1)
                if row is not None:
                    top = row["per_residue_topology"]
                    p = place_protein(model.sequence, top, path, acc)
                    if p["status"] == "placed":
                        cache[acc] = (top, p)
        if cache[acc] is None:
            continue
        top, p = cache[acc]
        n = int(s.insert_after_residue)
        if not 1 <= n <= len(top):
            continue
        ok = passing_junctions(top, p["z"], th)
        rows.append(
            dict(
                gene=s.gene,
                uniprot_acc=acc,
                insert_after_residue=n,
                v1_tm_helices=len(re.findall("M+", top)),
                junction_height_A=round(junction(p["z"], n), 1),
                passing_run_at_site=run_containing(ok, n),
                v1_O_run_at_site=next(
                    (b - a + 1 for a, b in ms.segments(top, "O") if a <= n <= b), 0
                ),
            )
        )
    return pd.DataFrame(rows)


def main() -> int:
    th = thresholds()
    gi.use_panel_topology()
    members = gi.internal_members()
    placed = member_placements(members)

    pub = published_loop_runs(th)
    multi = pub[(pub["v1_tm_helices"] >= 2) & (pub["passing_run_at_site"] > 0)]
    tiny = int(multi["passing_run_at_site"].min())

    gate = pd.read_csv(HERE / "l1_internal_pick_gate.tsv", sep="\t")
    shipped_pass = gate[(gate["shipped"] == 1) & (gate["verdict"] == "pass")]
    has_pass = set(shipped_pass["candidate_id"])
    rows = []
    for m in members.itertuples():
        p = placed[m.candidate_id]
        if p["status"] != "placed":
            continue
        ok = passing_junctions(m.topology, p["z"], th)
        longest = longest_run(ok)
        cz = c_terminal_side(m.topology, p["z"], p["plddt"])
        rows.append(
            dict(
                candidate_id=m.candidate_id,
                design_id=m.design_id,
                gene_symbol=m.gene_symbol,
                candidate_set=m.candidate_set,
                species=m.species,
                v1_tm_helices=len(re.findall("M+", m.topology)),
                longest_passing_run=longest,
                below_smallest_published_loop=int(longest < tiny),
                has_passing_shipped_site=int(m.candidate_id in has_pass),
                v1_c_terminus=m.topology[-1],
                c_terminal_flank_height_A=round(cz, 1),
                c_terminus_cytoplasmic_on_model=int(cz < 0) if not np.isnan(cz) else -1,
                topology_tool=m.topology_layer,
                topology_version=m.topology_version,
                topology_cohort=m.topology_cohort,
                decided_by="Becca, 2026-11-11 (consider a snorkel lane for tiny-loop members); "
                "cut-off derived from published sites",
                implemented_by="l1-membrane-side-gate session (Claude Code), PR #291",
            )
        )
    cand = pd.DataFrame(rows)
    cand.to_csv(HERE / "l1_snorkel_candidates.tsv", sep="\t", index=False)
    pub.to_csv(HERE / "snorkel_published_loop_runs.tsv", sep="\t", index=False)

    flagged = cand[cand["below_smallest_published_loop"] == 1]
    summary = dict(
        smallest_published_passing_run=tiny,
        published_multipass_sites=len(multi),
        published_passing_run_percentiles={
            q: float(multi["passing_run_at_site"].quantile(q / 100)) for q in (0, 5, 10, 25, 50)
        },
        published_smallest=multi.nsmallest(5, "passing_run_at_site")[
            ["gene", "insert_after_residue", "passing_run_at_site", "v1_O_run_at_site"]
        ].to_dict("records"),
        placed_internal_members=len(cand),
        members_below_smallest_published_loop=len(flagged),
        of_which_c_terminus_cytoplasmic=int(
            (flagged["c_terminus_cytoplasmic_on_model"] == 1).sum()
        ),
        of_which_without_passing_shipped_site=int((flagged["has_passing_shipped_site"] == 0).sum()),
        members=flagged.sort_values("longest_passing_run")[
            [
                "design_id",
                "gene_symbol",
                "longest_passing_run",
                "c_terminal_flank_height_A",
                "has_passing_shipped_site",
            ]
        ].to_dict("records"),
    )
    (HERE / "snorkel_summary.json").write_text(json.dumps(summary, indent=1, default=str) + "\n")
    print(json.dumps({k: v for k, v in summary.items() if k != "members"}, indent=1, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
