"""The membrane-side gate's verdict at each rescue-round site (S5E13, 2026-11-14).

The rescue round proposes tag sites on held members that are not on the internal rung:
mostly single-pass (type I/II) and GPI-anchored ectodomains, plus six multipass members with
a signal peptide. None has a DeepTMHMM v1 row, so the topology here is UniProt's: the
member's extracellular residues as the ranker published them, and the transmembrane and
intramembrane spans of the member's span reference (UniProt), mapped onto the member by
sequence alignment (identical residues only).

Placement, per member, on the member's own model (identical sequence, checked by SHA-256):

* TMDET, when TmAlphaFold holds the reference accession and its model superposes on ours
  within 2 A (the gate's rule), for any number of helices;
* otherwise the gate's helix-axis fit, which needs at least two confident TM helices.

What a verdict means depends on the protein, so each row says its ``basis``:

* ``multipass``: the gate's calibrated verdict, as for internal loops.
* ``single_pass_tmdet``: the height is real, but the gate's thresholds were calibrated on
  multipass proteins, and a single-pass ectodomain's angle to the membrane is not fixed by
  the model (a stalk or hinge lets it swing). A height at or above the outside threshold is
  a pass; anything lower is ``review_single_pass_geometry``, never a fail. Decided
  2026-11-14 by the l1-membrane-side-gate session (Claude Code) for the rescue round under
  Becca's 2026-11-14 request: the model cannot show that a low single-pass site is in the
  membrane, only that this model's pose puts it low.
* ``not_placeable``: a GPI anchor (a lipid, absent from the model) or one TM helix without a
  TMDET placement. The site rests on its sequence distance from the TM or omega site and on
  the refold.

Every row carries ``seq_distance_to_membrane``: residues from the junction to the nearest TM
or intramembrane residue, or to the omega site for GPI members.

    uv run --locked python gate_rescue_round.py --out-dir rescue-round-2026-11-14b
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import gate_inputs as gi  # noqa: E402
import membrane_side_gate as g  # noqa: E402
from placements import place_protein  # noqa: E402

from deliverome_analysis import membrane_side as ms  # noqa: E402

PREFIX = "sprint5/s5e13-tag-site-benchmarking/rescue_round/2026-11-14b/"
INPUTS = (
    "rescue_round_sites.tsv",
    "l1_rescue_member_features.tsv",
    "l1_rescue_member_extracellular.tsv",
)
DECIDED_BY = (
    "l1-membrane-side-gate session (Claude Code), rescue round under Becca's 2026-11-14 request"
)
IMPLEMENTED_BY = "l1-membrane-side-gate session (Claude Code)"


def read_inputs() -> tuple[dict[str, pd.DataFrame], dict[str, str]]:
    from deliverome_analysis.s3 import processed

    out, shas = {}, {}
    for name in INPUTS:
        p = Path(processed.download(PREFIX + name))
        shas[name] = hashlib.sha256(p.read_bytes()).hexdigest()
        out[name] = pd.read_csv(p, sep="\t", low_memory=False)
    return out, shas


def uniprot_membrane_spans(acc: str) -> list[tuple[int, int]]:
    """1-based Transmembrane and Intramembrane spans UniProt lists for ``acc``."""
    p = gi.fetch(
        f"https://rest.uniprot.org/uniprotkb/{acc}.json?fields=ft_transmem,ft_intramem",
        gi.CACHE / "uniprot" / f"{acc}.membrane.json",
    )
    spans = []
    for f in json.loads(p.read_text()).get("features", []) if p else []:
        if f["type"] in ("Transmembrane", "Intramembrane"):
            s, e = f["location"]["start"]["value"], f["location"]["end"]["value"]
            if s is not None and e is not None:
                spans.append((int(s), int(e)))
    return spans


def model_for(feat: pd.Series) -> Path:
    ref = str(feat["model_ref"])
    if feat["model_source"] == "alphafold_db":
        return gi.fetch(ref, gi.CACHE / "models" / Path(ref).name)
    from deliverome_analysis.s3 import processed

    key = ref.removeprefix("processed://")
    dest = gi.CACHE / "models" / ("colabfold_" + key.replace("/", "_"))
    if not dest.is_file():
        dest.write_bytes(Path(processed.download(key)).read_bytes())
    return dest


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args(argv)
    out = HERE / args.out_dir if not args.out_dir.is_absolute() else args.out_dir
    out.mkdir(parents=True, exist_ok=True)

    tabs, shas = read_inputs()
    sites = tabs["rescue_round_sites.tsv"]
    feats = (
        tabs["l1_rescue_member_features.tsv"].drop_duplicates("design_id").set_index("design_id")
    )
    ec = tabs["l1_rescue_member_extracellular.tsv"]
    ec["member_id"] = ec["uniprot_acc"].str.split("@").str[0]
    th = g.thresholds()

    members = {}
    for d, s in sites.groupby("design_id"):
        s0 = s.iloc[0]
        feat = feats.loc[d]
        path = model_for(feat)
        model = ms.read_model(path)
        seq_sha = hashlib.sha256(model.sequence.encode()).hexdigest()
        if seq_sha != s0["sequence_sha256"]:
            members[d] = dict(status="model_sequence_differs_from_member")
            continue
        seq = model.sequence
        L = len(seq)
        topo = ["I"] * L
        # Input join only: both input tables carry the ranker snapshot's member ids (a22130a2).
        for r in ec.loc[ec["member_id"] == s0["member_id"], "residue"]:
            topo[int(r) - 1] = "O"
        ref = str(s0["span_reference"])
        ref_seq = gi.uniprot_sequence(ref)
        rmap, _ = ms.map_to_model(seq, ref_seq) if ref_seq else ([None] * L, "none")
        ref_to_member = {j: i for i, j in enumerate(rmap) if j is not None}
        tm = set()
        for a, b in uniprot_membrane_spans(ref):
            for k in range(a - 1, b):
                if k in ref_to_member:
                    tm.add(ref_to_member[k])
        for i in tm:
            topo[i] = "M"
        topo = "".join(topo)
        n_tm_model = len(ms.segments(topo, "M"))
        omega = s0["omega_site_member"]
        info = dict(
            seq=seq,
            topo=topo,
            tm=sorted(tm),
            omega=None if pd.isna(omega) else int(omega),
            n_tm=int(s0["n_tm"]),
            n_tm_mapped=n_tm_model,
            model=f"{feat['model_source']}:{feat['model_ref']}",
        )
        if not tm:
            members[d] = info | dict(
                status="not_placeable",
                why="GPI-anchored: no transmembrane helix; the anchor is a lipid the model does not contain"
                if info["omega"]
                else "no transmembrane span maps onto the member",
            )
            continue
        afdb = ref if feat["model_source"] != "alphafold_db" else str(feat["sequence_accession"])
        p = place_protein(seq, topo, path, afdb_acc=afdb)
        # A single-pass member is placeable only by TMDET: the helix-axis fit needs two
        # helices, and an alignment gap can split one UniProt TM into two segments
        # (BRICD5), which must not count as two.
        single = info["n_tm"] < 2
        if p["status"] != "placed" or (single and p["placement_source"] != "tmdet"):
            rmsd = p.get("tmdet_rmsd", np.nan)
            tm_why = (
                "no TmAlphaFold entry"
                if rmsd is None or np.isnan(rmsd)
                else f"TmAlphaFold model {rmsd:.2f} A RMSD from ours, above the 2 A rule"
            )
            why = (
                f"one TM helix and no TMDET placement ({tm_why}); "
                "the helix-axis fit needs two helices"
                if single
                else p["status"]
            )
            members[d] = info | dict(status="not_placeable", why=why, placement=p)
            continue
        members[d] = info | dict(status="placed", placement=p)

    rows = []
    for r in sites.itertuples():
        m = members[r.design_id]
        n = int(r.insert_after_residue)
        base = dict(
            design_id=r.design_id,
            gene=r.gene,
            species=r.species,
            rescue_rank=r.rescue_rank,
            insert_after_residue=n,
            junction=r.junction,
            sequence_sha256=r.sequence_sha256,
            n_tm=r.n_tm,
        )
        if m["status"] == "model_sequence_differs_from_member":
            rows.append(base | dict(verdict="unscorable", basis=m["status"]))
            continue
        assert m["seq"][n - 1] + "|" + m["seq"][n] == r.junction, (r.design_id, n)
        dists = [abs(n - (i + 1)) for i in m["tm"]] + ([abs(n - m["omega"])] if m["omega"] else [])
        base |= dict(
            seq_distance_to_membrane=min(dists) if dists else np.nan,
            junction_topology=m["topo"][n - 1] + m["topo"][n],
            model=m["model"],
        )
        if m["status"] == "not_placeable":
            rows.append(
                base | dict(verdict="not_placeable", basis="not_placeable", reason=m["why"])
            )
            continue
        p = m["placement"]
        z = g.junction(p["z"], n)
        pl = float(np.nanmean(p["plddt"][n - 1 : n + 1]))
        ok = ms.orientation_supported(p["topo_for"], p["topo_against"], p["kr_out"], p["kr_in"], th)
        v = ms.site_verdict(z, pl, ok, th)
        multipass = m["n_tm"] >= 2
        basis = "multipass" if multipass else "single_pass_tmdet"
        if not multipass and v not in ("pass", "review_orientation", "unscorable"):
            v = "review_single_pass_geometry"
        rows.append(
            base
            | dict(
                verdict=v,
                height_A=round(z, 1),
                plddt=round(pl, 1),
                basis=basis,
                placement_source=p["placement_source"],
                tmdet_rmsd=p.get("tmdet_rmsd"),
                orientation_supported=ok,
                topo_for=p["topo_for"],
                topo_against=p["topo_against"],
                kr_out=p["kr_out"],
                kr_in=p["kr_in"],
                reason=f"junction at {z:+.1f} A, pLDDT {pl:.0f}; thresholds in calibration.json",
            )
        )
    t = pd.DataFrame(rows)
    t["topology_source"] = (
        "UniProt (span reference, aligned) + ranker's extracellular residues; not DeepTMHMM v1"
    )
    t["decided_by"], t["implemented_by"] = DECIDED_BY, IMPLEMENTED_BY
    t = t.sort_values(["design_id", "rescue_rank"])
    t.to_csv(out / "rescue_round_gate.tsv", sep="\t", index=False)
    summary = dict(
        inputs={PREFIX + k: v for k, v in shas.items()},
        thresholds=th.__dict__,
        sites=len(t),
        members=int(t["design_id"].nunique()),
        by_basis_and_verdict={
            f"{b}|{v}": int(k) for (b, v), k in t.groupby(["basis", "verdict"]).size().items()
        },
        not_passing=t[t["verdict"] != "pass"][
            ["design_id", "gene", "insert_after_residue", "verdict", "height_A", "basis"]
        ].to_dict("records"),
    )
    (out / "rescue_round_gate_summary.json").write_text(
        json.dumps(summary, indent=1, default=str) + "\n"
    )
    print(json.dumps(summary["by_basis_and_verdict"], indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
