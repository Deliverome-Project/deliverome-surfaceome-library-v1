"""Calibrate the membrane-side gate against independent evidence; derive its thresholds.

Three questions, three evidence sets, none of which the placement itself uses:

1. **Is the slab where TMDET puts it?** TMDET (the algorithm behind PDBTM) has been run
   on AlphaFold DB models by TmAlphaFold. For every AlphaFold-DB-sourced member with a
   TmAlphaFold entry, our model is superposed on TmAlphaFold's and each residue's height
   is read in both frames. This measures the placement error, as a function of pLDDT.
2. **Is the extracellular face the right one?** UniProt topological domains backed by
   experimental evidence (ECO:0000269) say which side a stretch is on. Read on the
   members whose designed protein is the UniProt canonical sequence.
3. **How deep do tags that work sit?** The audited literature gold standard
   (``processed://scratch/becca/tag-site-corpus-audit/gold_standard.tsv``): internal tag
   insertions confirmed at the cited residue, 120 of them with surface display shown.
   Each is placed on its protein's AlphaFold DB model with the same code as the members.

Writes ``calibration.json`` (the thresholds the gate reads, and every number they rest
on), ``calibration_tmdet.tsv``, ``calibration_uniprot_side.tsv`` and
``calibration_published_sites.tsv``, and the figure ``figures/membrane_side_calibration.png``.

    uv run --locked python calibrate.py
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
from placements import member_placements, place_protein  # noqa: E402

from deliverome_analysis.surfaceome import membrane_side as ms  # noqa: E402
from deliverome_analysis.surfaceome import topology  # noqa: E402

GOLD_PREFIX = "scratch/becca/tag-site-corpus-audit/"
EBI_FEATURES = "https://www.ebi.ac.uk/proteins/api/features/{}?categories=TOPOLOGY"
PLDDT_BINS = [0, 50, 70, 90, 101]
PLDDT_LABELS = ["<50", "50-70", "70-90", ">=90"]


def junction(z: np.ndarray, n: int) -> float:
    """Height of an insert-after-``n`` junction (1-based): mean of residues n and n+1."""
    vals = [z[n - 1]] + ([z[n]] if n < len(z) else [])
    vals = [v for v in vals if not np.isnan(v)]
    return float(np.mean(vals)) if vals else float("nan")


# --- 1. TMDET -----------------------------------------------------------------------


def tmdet_comparison(members: pd.DataFrame, placed: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Helix-axis placement against TMDET, on the members TMDET placed.

    Per member: TMDET status, superposition RMSD, normal angle. Per residue: both heights,
    each oriented by the same topology vote, so a side flip shows as a sign difference.
    """
    prot_rows, res_rows = [], []
    for r in members.itertuples():
        p = placed[r.candidate_id]
        if p.get("model_source") != "alphafold_db" or p.get("status") != "placed":
            continue
        base = dict(candidate_id=r.candidate_id, gene_symbol=r.gene_symbol)
        if "tmdet_half_thickness" not in p:
            prot_rows.append(base | dict(tmdet="no_entry"))
            continue
        used = p["placement_source"] == "tmdet"
        prot_rows.append(
            base
            | dict(
                tmdet="used" if used else "different_model",
                tmdet_rmsd=p["tmdet_rmsd"],
                tmdet_half_thickness=p["tmdet_half_thickness"],
                helix_axis_vs_tmdet_deg=p.get("helix_axis_vs_tmdet_deg", np.nan),
            )
        )
        if not used or "z_helix_axis" not in p:
            continue
        for i in range(len(r.protein)):
            zt, zh = p["z"][i], p["z_helix_axis"][i]
            if np.isnan(zt) or np.isnan(zh):
                continue
            res_rows.append((r.candidate_id, i + 1, r.topology[i], float(p["plddt"][i]), zh, zt))
    res = pd.DataFrame(res_rows, columns=["candidate_id", "residue", "v1", "plddt", "z", "z_tmdet"])
    return pd.DataFrame(prot_rows), res


# --- 2. UniProt experimental sides --------------------------------------------------


def uniprot_topo(acc: str) -> list[tuple[int, int, str, bool]]:
    """(start, end, 'O'/'I', experimental?) per TOPO_DOM of the canonical entry."""
    p = gi.fetch(EBI_FEATURES.format(acc), gi.CACHE / "ebi" / f"{acc}.json")
    if p is None:
        return []
    out = []
    for f in json.loads(p.read_text()).get("features", []):
        if f.get("type") != "TOPO_DOM":
            continue
        d = (f.get("description") or "").lower()
        side = (
            "O"
            if d.startswith(("extracellular", "lumenal", "vesicular"))
            else ("I" if d.startswith("cytoplasmic") else None)
        )
        if side is None:
            continue
        exp = any(e.get("code") == "ECO:0000269" for e in f.get("evidences", []))
        try:
            out.append((int(f["begin"]), int(f["end"]), side, exp))
        except KeyError, ValueError:
            continue
    return out


def uniprot_side(members: pd.DataFrame, placed: dict, z_side: float) -> pd.DataFrame:
    chosen = gi.chosen_models()
    rows = []
    for r in members.itertuples():
        p = placed[r.candidate_id]
        if p["status"] != "placed" or r.candidate_set != "library_reference":
            continue
        if p["model_source"] != "alphafold_db" or p["model_match"] != "identical":
            continue
        acc = re.search(r"AF-(.+?)-F1", str(chosen.loc[r.candidate_id, "model_url"])).group(1)
        doms = [d for d in uniprot_topo(acc) if d[3]]
        if not doms:
            continue
        agree = disagree = 0
        for a, b, side, _ in doms:
            for i in range(a - 1, min(b, len(r.protein))):
                zi, pl = p["z"][i], p["plddt"][i]
                if np.isnan(zi) or pl < ms.FIT_MIN_PLDDT or abs(zi) < z_side:
                    continue
                if (zi > 0) == (side == "O"):
                    agree += 1
                else:
                    disagree += 1
        if agree + disagree == 0:
            continue
        rows.append(
            dict(
                candidate_id=r.candidate_id,
                design_id=r.design_id,
                gene_symbol=r.gene_symbol,
                uniprot_acc=acc,
                residues_agree=agree,
                residues_disagree=disagree,
                protein_side="agrees" if agree > disagree else "inverted",
                side_source=p["side_source"],
                positive_inside_agrees=p["positive_inside_agrees"],
            )
        )
    return pd.DataFrame(rows)


# --- 3. Published sites -------------------------------------------------------------


def published_sites() -> pd.DataFrame:
    from deliverome_analysis.core.s3 import processed

    gold = pd.read_csv(processed.download(GOLD_PREFIX + "gold_standard.tsv"), sep="\t")
    failed = pd.read_csv(processed.download(GOLD_PREFIX + "failed_sites.tsv"), sep="\t")
    gold["set"] = "gold_standard"
    failed["set"] = "failed_site"
    failed["surface_shown"] = "construct_failed"
    sites = pd.concat([gold, failed], ignore_index=True)
    rows = []
    cache: dict[str, dict | None] = {}
    for s in sites.itertuples():
        acc = s.uniprot_acc
        if acc not in cache:
            path = gi.afdb_model(acc)
            if path is None:
                cache[acc] = None
            else:
                model = ms.read_model(path)
                row = topology.for_sequence(model.sequence, tool=topology.TOOL_V1)
                cache[acc] = (
                    None
                    if row is None
                    else place_protein(model.sequence, row["per_residue_topology"], path, acc)
                    | dict(topology=row["per_residue_topology"])
                )
        p = cache[acc]
        base = dict(
            set=s.set,
            uniprot_acc=acc,
            gene=s.gene,
            insert_after_residue=int(s.insert_after_residue),
            surface_shown=s.surface_shown,
        )
        if p is None or p.get("status") != "placed":
            rows.append(base | dict(status="no_model_or_v1" if p is None else p["status"]))
            continue
        n = int(s.insert_after_residue)
        if not 1 <= n <= len(p["z"]):
            rows.append(base | dict(status="residue_out_of_range"))
            continue
        rows.append(
            base
            | dict(
                status="placed",
                v1_state=p["topology"][n - 1],
                z=round(junction(p["z"], n), 2),
                plddt=round(float(p["plddt"][n - 1]), 1),
                placement_source=p["placement_source"],
                v1_tm_helices=len(re.findall("M+", p["topology"])),
                topo_for=p["topo_for"],
                topo_against=p["topo_against"],
                kr_out=p["kr_out"],
                kr_in=p["kr_in"],
            )
        )
    return pd.DataFrame(rows)


def derive(
    prot: pd.DataFrame,
    res: pd.DataFrame,
    pub: pd.DataFrame,
    side: pd.DataFrame,
    members: pd.DataFrame,
    placed: dict,
) -> dict:
    """Every threshold the gate uses, each from one stated piece of evidence."""
    ev: dict = {}
    ev["tmdet_status_counts"] = prot["tmdet"].value_counts().to_dict()
    used = prot[prot["tmdet"] == "used"]
    ev["tmdet_rmsd_median"] = float(prot["tmdet_rmsd"].median())
    ev["helix_axis_vs_tmdet_deg_median"] = float(used["helix_axis_vs_tmdet_deg"].median())
    ev["helix_axis_vs_tmdet_deg_p90"] = float(used["helix_axis_vs_tmdet_deg"].quantile(0.9))
    ev["tmdet_half_thickness_median"] = float(used["tmdet_half_thickness"].median())
    ev["tmdet_half_thickness_p10_p90"] = [
        float(used["tmdet_half_thickness"].quantile(q)) for q in (0.1, 0.9)
    ]
    ev["side_flips_helix_axis_vs_tmdet_members"] = int(
        res.assign(flip=np.sign(res["z"]) != np.sign(res["z_tmdet"]))
        .loc[lambda d: (d["z"].abs() > 15) & (d["z_tmdet"].abs() > 15)]
        .groupby("candidate_id")["flip"]
        .mean()
        .gt(0.5)
        .sum()
    )

    # Placement error against TMDET, on loop residues (v1 not M) within 40 A of the
    # mid-plane -- the residues a tag site can be -- by pLDDT of the residue.
    loops = res[(res["v1"] != "M") & (res["z_tmdet"].abs() <= 40)].copy()
    loops["abs_dz"] = (loops["z"] - loops["z_tmdet"]).abs()
    edges = [0, 30, 40, 50, 60, 70, 80, 90, 101]
    loops["plddt_bin"] = pd.cut(loops["plddt"], edges, right=False)
    by_bin = loops.groupby("plddt_bin", observed=True)["abs_dz"].describe(percentiles=[0.5, 0.9])
    ev["tmdet_abs_dz_by_plddt"] = {
        str(k): dict(n=int(v["count"]), median=round(v["50%"], 2), p90=round(v["90%"], 2))
        for k, v in by_bin.iterrows()
    }
    # Side agreement with TMDET-in-sign-aligned frame: residues our placement calls
    # beyond the core on one side that TMDET puts beyond its own boundary on the other.
    ev["tmdet_loop_residues"] = len(loops)
    conf_loops = loops[loops["plddt"] >= 50]

    # 1. core: TMDET's own hydrophobic half-thickness, the median over the members it
    # placed. A junction this close to the mid-plane is in the membrane by TMDET's measure.
    core = round(ev["tmdet_half_thickness_median"], 1)

    # 2. outside: the shallowest published, surface-displayed tag site in a protein whose
    # placement is defined the same way the members' are (>= 2 v1 helices, side backed by
    # the topology vote). No working published site would be refused at this height.
    pub_ok = pub[
        (pub["status"] == "placed")
        & (pub["set"] == "gold_standard")
        & (pub["surface_shown"] == "yes")
        & (pub["v1_tm_helices"] >= 2)
    ]
    ev["published_displayed_multipass_sites"] = len(pub_ok)
    ev["published_displayed_multipass_z_percentiles"] = {
        q: round(float(pub_ok["z"].quantile(q / 100)), 1) for q in (0, 1, 5, 10, 25, 50)
    }
    shallow = pub_ok.nsmallest(5, "z")[["gene", "insert_after_residue", "z", "plddt"]]
    ev["published_displayed_shallowest"] = shallow.to_dict("records")
    outside = float(np.floor(pub_ok["z"].min() * 2) / 2)

    # margin: the median helix-axis error against TMDET at pLDDT >= 50. A junction this
    # close below the pass line is within the typical error of passing: review, not fail.
    margin = round(float(conf_loops["abs_dz"].median()), 1)
    ev["helix_axis_abs_dz_median_plddt50"] = margin

    # 3. min_plddt: the lowest pLDDT bin, among an unbroken run from the top, whose
    # 90th-percentile error against TMDET is reliably smaller than `core`: the upper end of
    # a 95% bootstrap interval for that percentile, resampling whole proteins (residues of
    # one protein are not independent), must be below `core`. Below it, an error of that
    # size can carry a junction from the cytoplasmic side across the mid-plane, so a
    # failing call is only a hold. A point estimate alone flipped this cut-off between 40
    # and 0 when 15 members were added on 2026-11-12, because the [30, 40) bin sits at
    # 13.0 +/- 0.3 A; the interval makes the rule stable (l1-membrane-side-gate, PR #291).
    rng = np.random.default_rng(20261112)
    min_plddt = 101.0
    ev["tmdet_p90_upper95_by_plddt"] = {}
    for interval, g in sorted(
        loops.groupby("plddt_bin", observed=True), key=lambda kv: -kv[0].left
    ):
        prots = g["candidate_id"].unique()
        by_prot = {c: d["abs_dz"].to_numpy() for c, d in g.groupby("candidate_id")}
        boots = []
        for _ in range(200):
            pick = rng.choice(prots, size=len(prots), replace=True)
            boots.append(np.percentile(np.concatenate([by_prot[c] for c in pick]), 90))
        upper = float(np.percentile(boots, 97.5))
        ev["tmdet_p90_upper95_by_plddt"][str(interval)] = round(upper, 2)
        if upper < core and len(g) >= 50:
            min_plddt = float(interval.left)
        else:
            break

    # 4/5. Orientation: the weakest topology vote, and the largest flank K/R excess on the
    # extracellular face, seen in members whose face UniProt experimental topology confirms.
    conf = side[side["protein_side"] == "agrees"]
    ratios, excess = [], []
    for cid in conf["candidate_id"]:
        p = placed[cid]
        ratios.append(p["topo_for"] / max(p["topo_against"], 1))
        excess.append(p["kr_out"] - p["kr_in"])
    ev["uniprot_confirmed_members"] = len(conf)
    ev["uniprot_inverted_members"] = side.loc[
        side["protein_side"] != "agrees", "gene_symbol"
    ].tolist()
    ev["uniprot_residues_agree_disagree"] = [
        int(side["residues_agree"].sum()),
        int(side["residues_disagree"].sum()),
    ]
    ev["confirmed_min_topology_ratio"] = round(min(ratios), 2)
    ev["confirmed_max_kr_excess_outside"] = int(max(excess))
    th = dict(
        outside=outside,
        margin=margin,
        core=core,
        min_plddt=min_plddt,
        min_topology_ratio=float(np.floor(min(ratios) * 10) / 10),
        max_kr_excess_outside=int(max(excess)),
    )
    attribution = {
        k: dict(
            decided_by="derived by calibrate.py from the evidence above (rule written by the "
            "l1-membrane-side-gate session, PR #291, 2026-11-10, to Becca's requirement that "
            "every cutoff be justified from evidence)",
            implemented_by="l1-membrane-side-gate session (Claude Code), PR #291",
        )
        for k in th
    }
    return dict(
        thresholds=th,
        attribution=attribution,
        evidence=ev,
        # By design id: member ids are renumbered between panels (2026-11-14: l1-m1545 was
        # SLC68A1 at 6bf341f2 and SLC7A10 on main), design ids are not.
        uniprot_inverted_design_ids=side.loc[
            side["protein_side"] != "agrees", "design_id"
        ].tolist(),
    )


def figure(res: pd.DataFrame, pub: pd.DataFrame, cal: dict) -> None:
    import matplotlib.pyplot as plt

    from deliverome_analysis.core.plotting import save_figure, setup_plotting_style

    setup_plotting_style()
    th = cal["thresholds"]
    fig, (a, b) = plt.subplots(1, 2, figsize=(13, 4.2))
    loops = res[(res["v1"] != "M") & (res["z_tmdet"].abs() <= 40)]
    for lo, hi, lab in ((70, 101, "pLDDT ≥ 70"), (50, 70, "pLDDT 50-70"), (0, 50, "pLDDT < 50")):
        d = loops[(loops["plddt"] >= lo) & (loops["plddt"] < hi)]
        dens, edges = np.histogram((d["z"] - d["z_tmdet"]).clip(-20, 20), bins=80, density=True)
        a.step(edges[:-1], dens, where="post", label=f"{lab} (n = {len(d):,})")
    a.set_xlabel("Height, helix-axis fit minus TMDET (Å)")
    a.set_ylabel("Density")
    a.set_ylim(0, a.get_ylim()[1] * 1.35)
    a.legend(frameon=False, fontsize=7, loc="upper left")
    ok = pub[
        (pub["status"] == "placed") & (pub["set"] == "gold_standard") & (pub["v1_tm_helices"] >= 2)
    ]
    groups = [
        ("yes", "Surface display shown"),
        ("not_reported", "Display not reported"),
        ("no", "Shown not displayed"),
    ]
    for i, (k, lab) in enumerate(groups):
        d = ok[ok["surface_shown"] == k]
        b.scatter(
            d["z"],
            np.full(len(d), i) + np.random.default_rng(i).uniform(-0.25, 0.25, len(d)),
            s=12,
            alpha=0.7,
            label=f"{lab} (n = {len(d)})",
        )
    o, r = th["outside"], th["outside"] - th["margin"]
    b.axvline(o, ls="--", lw=1, color="0.3", label=f"pass at {o:g} Å")
    b.axvline(r, ls=":", lw=1, color="0.3", label=f"review from {r:g} Å")
    b.axvline(-th["core"], ls="-.", lw=1, color="0.3", label=f"cytoplasmic below {-th['core']:g} Å")
    b.legend(frameon=False, fontsize=7, loc="center left", bbox_to_anchor=(1.01, 0.5))
    b.set_yticks(range(len(groups)), [g[1] for g in groups])
    b.set_xlabel("Junction height above mid-plane (Å)")
    fig.tight_layout()
    save_figure(fig, HERE / "figures" / "membrane_side_calibration")


def main() -> int:
    gi.use_panel_topology()
    members = gi.internal_members()
    placed = member_placements(members)

    prot, res = tmdet_comparison(members, placed)
    prot.to_csv(HERE / "calibration_tmdet.tsv", sep="\t", index=False)
    pub = published_sites()
    pub.to_csv(HERE / "calibration_published_sites.tsv", sep="\t", index=False)
    side = uniprot_side(members, placed, z_side=ms.HALF_THICKNESS)
    side.to_csv(HERE / "calibration_uniprot_side.tsv", sep="\t", index=False)
    cal = derive(prot, res, pub, side, members, placed)
    cal["inputs"] = dict(
        panel_commit=gi.PANEL_COMMIT,
        topology_snapshot=gi.TOPOLOGY_SNAPSHOT_KEY,
        topology_tool=topology.TOOL_V1,
        gold_standard="processed://" + GOLD_PREFIX + "gold_standard.tsv",
        tmalphafold="https://tmalphafold.ttk.hu/api/tmdet/{acc}.json (TMDET on AlphaFold DB models)",
        uniprot_features="EBI Proteins API, TOPO_DOM with ECO:0000269",
    )
    (HERE / "calibration.json").write_text(json.dumps(cal, indent=1, default=str) + "\n")
    figure(res, pub, cal)
    print(json.dumps(cal["thresholds"], indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
