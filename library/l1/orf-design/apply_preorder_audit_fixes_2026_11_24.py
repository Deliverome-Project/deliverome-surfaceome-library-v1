#!/usr/bin/env python3
"""Apply the 2026-11-24 pre-order biological audit fixes to the L1 internal arm tables.

Fixes all 14 internal-tagged CDSs caught by the pre-order audit (0 CDSs lost; -9 defective or
redundant arms; 8 one-for-one swaps to passing refolded sites):

* Cytoplasmic loop arm (1 dropped): ANO4 l1-p0113-a6 (@685, z=-37.9 A TMDET).
* Membrane-core jmask=m arms (4 dropped): SLC44A1 l1-p1961-a3/a4 (@198/@103) and isoform
  l1-p3281-a4/a5 (@198/@103).
* Cleaved N-terminal propeptide arms (3 swapped, zero-GPU): ATRN @79 -> @1038 (l1-p0163-a1),
  CDH3 @90 -> @645 (l1-p0347-a4), CDH5 @27 -> @586 (l1-p0349-a1).
* GPI omega-proximal arms (<12 aa before omega; 1 dropped + 3 swapped): LYPD6 drop a1@146,
  MDGA2 @998 -> @151 (l1-p1204-a2, zero-GPU), BST2 @155 -> @149 (l1-p0188-a3, Modal pass),
  PRND @149 fail -> @97 (l1-p1675-a3, Modal pass).
* Same-loop <25 aa close pairs (3 dropped + 2 swapped, zero-GPU): TRPV3 @473 -> @615
  (l1-p2278-a2), XPR1 @426 -> @503 (l1-p2354-a3), HCN2 drop a4@401, AQP4 drop a4@140,
  SLC7A11 drop a3@312.
"""

from __future__ import annotations

import csv
import hashlib
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
S5E13 = HERE.parents[1] / "l1/tag-site-benchmarking"
ARM_DECISIONS = S5E13 / "l1_internal_arm_decisions.tsv"
ARM_IDS = S5E13 / "l1_internal_arm_ids.tsv"
VERDICTS_TSV = HERE / "tag_refold_verdicts.tsv"
STRUCTURES_TSV = HERE / "tag_refold_structures.tsv"

sys.path.insert(0, str(HERE / "tag_refold"))
import score_tagged_folds as stf  # noqa: E402
from build_refold_inputs import CASSETTE  # noqa: E402

from deliverome_analysis.cloning import library_design as ld  # noqa: E402

DECIDED_BY = "Becca Carlson, 2026-11-24 (pre-order audit)"
IMPLEMENTED_BY = "pre-order biological audit (2026-11-24)"
DECIDED_ON = "2026-11-24"
S3_PREFIX = (
    "processed://sprint4/s4e9-orf-design-fg-panel/"
    "predicted-structures/tagged/2026-11-24-pre-order-audit/"
)

NEW_MODAL_ARMS = [
    ("l1-p0188", "BST2", "l1-p0188-a3", 149, "t.5549ef4a27f2485b"),
    ("l1-p1675", "PRND", "l1-p1675-a3", 97, "t.95d8abd64639a995"),
]

DROP_ARMS = {
    "l1-p0113-a6",  # ANO4 @685 cytoplasmic (z=-37.9 A)
    "l1-p1961-a3",  # SLC44A1 @198 membrane core
    "l1-p1961-a4",  # SLC44A1 @103 membrane core
    "l1-p3281-a4",  # SLC44A1 iso @198 membrane core
    "l1-p3281-a5",  # SLC44A1 iso @103 membrane core
    "l1-p1181-a1",  # LYPD6 @146 GPI omega-proximal + N147 glycan
    "l1-p0846-a4",  # HCN2 @401 close pair (5 aa from @396)
    "l1-p0131-a4",  # AQP4 @140 close pair (12 aa from @152)
    "l1-p2015-a3",  # SLC7A11 @312 close pair (13 aa from @299)
}

# old_arm_id -> (new_arm_id, new_site, new_flag, reason)
SWAP_ARMS: dict[str, tuple[str, int, str, str]] = {
    "l1-p0163-a2": (
        "l1-p0163-a1",
        1038,
        "audit_flag: membrane-side gate usable=0 (single-pass; site in mature ectodomain 84-1279)",
        "pre-order audit swap: move out of cleaved propeptide 29-83 to passing ectodomain @1038",
    ),
    "l1-p0347-a5": (
        "l1-p0347-a4",
        645,
        "",
        "pre-order audit swap: move out of cleaved furin propeptide 25-107 to passing stalk @645",
    ),
    "l1-p0349-a2": (
        "l1-p0349-a1",
        586,
        "",
        "pre-order audit swap: move out of cleaved furin propeptide 26-47 to passing stalk @586",
    ),
    "l1-p1204-a3": (
        "l1-p1204-a2",
        151,
        "",
        "pre-order audit swap: move out of GPI omega motif (@998 vs omega 931/1000) to passing @151",
    ),
    "l1-p0188-a2": (
        "l1-p0188-a3",
        149,
        "",
        "pre-order audit swap: move >=12 aa upstream of GPI omega=161 to Modal-passing @149",
    ),
    "l1-p1675-a2": (
        "l1-p1675-a3",
        97,
        "",
        "pre-order audit swap: move 55 aa upstream of GPI omega=152 to Modal-passing @97",
    ),
    "l1-p2278-a3": (
        "l1-p2278-a2",
        615,
        "",
        "pre-order audit swap: replace close pair @473 (3 aa from @470) with distinct-loop @615",
    ),
    "l1-p2354-a5": (
        "l1-p2354-a3",
        503,
        "",
        "pre-order audit swap: replace close pair @426 (8 aa from @434) with distinct-loop @503",
    ),
}

# Clear stale close-pair audit_flags on the surviving partners of the 5 resolved pairs
CLEAR_FLAGS = {
    "l1-p2278-a1",
    "l1-p0846-a3",
    "l1-p2354-a4",
    "l1-p0131-a5",
    "l1-p2015-a4",
}


def register_new_arm_ids() -> None:
    rows = list(csv.DictReader(ARM_IDS.open(), delimiter="\t"))
    existing = {r["arm_id"] for r in rows}
    added = 0
    for did, _gene, aid, pos, _fid in NEW_MODAL_ARMS:
        if aid not in existing:
            rows.append({"arm_id": aid, "design_id": did, "construct_key": f"{did}:internal{pos}"})
            added += 1
    with ARM_IDS.open("w", newline="") as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=["arm_id", "design_id", "construct_key"],
            delimiter="\t",
            lineterminator="\n",
        )
        w.writeheader()
        w.writerows(rows)
    print(f"arm_ids: +{added} registered ({len(rows)} total)")


def register_modal_refolds(fold_dir: Path) -> None:
    q = {
        r["fold_id"]: r
        for r in csv.DictReader(Path("/tmp/gpi_omega_rescue_queue_v2.tsv").open(), delimiter="\t")
    }
    ni = {
        r["native_sha256"]: r
        for r in csv.DictReader((HERE / "tag_refold_native_models.tsv").open(), delimiter="\t")
    }
    v_rows = list(csv.DictReader(VERDICTS_TSV.open(), delimiter="\t"))
    v_cols = list(v_rows[0].keys())
    v_seen = {r["tagged_sha256"] for r in v_rows}

    s_rows = list(csv.DictReader(STRUCTURES_TSV.open(), delimiter="\t"))
    s_cols = list(s_rows[0].keys())
    s_seen = {r["tagged_sha256"] for r in s_rows}

    for did, gene, aid, pos, fid in NEW_MODAL_ARMS:
        item = q[fid]
        res = stf.score(item, fold_dir, ni)
        pdb_path = fold_dir / f"{fid}.pdb"
        pdb_sha = hashlib.sha256(pdb_path.read_bytes()).hexdigest()
        tsha = item["tagged_sha256"]
        nsha = item["native_sha256"]
        pdb_key = f"{S3_PREFIX}{fid}.pdb"
        scores_key = f"{S3_PREFIX}{fid}.scores.json.gz"

        if tsha not in v_seen:
            v_rec = {c: "" for c in v_cols}
            v_rec.update(
                {
                    "tagged_sha256": tsha,
                    "rung": "internal",
                    "source_id": aid,
                    "design_id": did,
                    "gene_symbol": gene,
                    "tier": "1",
                    "reason": "GPI omega-clearance pre-order audit fix (2026-11-24)",
                    "arm": aid,
                    "design": "v8",
                    "rank": "1",
                    "source": "designed protein + cassette (pre-order audit batch, 2026-11-24)",
                    "pdb_key": pdb_key,
                    "verdict": res["verdict"],
                    "verdict_reasons": res["verdict_reasons"],
                    "lddt_confident": res["lddt_confident"],
                    "lddt_site": res["lddt_site"],
                    "tm_native": res["tm_native"],
                    "tm_confident": res["tm_confident"],
                    "rmsd_confident": res["rmsd_confident"],
                    "dplddt_native": res["dplddt_native"],
                    "dplddt_site": res["dplddt_site"],
                    "mean_plddt": res["mean_plddt"],
                    "cassette_plddt": res["cassette_plddt"],
                    "spytag003_plddt": res["spytag003_plddt"],
                    "alfa_plddt": res["alfa_plddt"],
                    "spytag003_rsa": res["spytag003_rsa"],
                    "alfa_rsa": res["alfa_rsa"],
                    "cassette_contact_res": res["cassette_contact_res"],
                    "epitope_confident_contact_res": res["epitope_confident_contact_res"],
                    "contacted_native_residues": res.get("contacted_native_residues", ""),
                    "cassette_longest_helix_excl_alfa": res["cassette_longest_helix_excl_alfa"],
                    "spytag003_longest_helix": res["spytag003_longest_helix"],
                    "native_model": res["native_model"],
                    "decided_by": DECIDED_BY,
                }
            )
            v_rows.append(v_rec)
            v_seen.add(tsha)

        if tsha not in s_seen:
            s_rows.append(
                {
                    "fold_id": fid,
                    "kind": "tagged",
                    "rung": "loop",
                    "design_id": did,
                    "gene_symbol": gene,
                    "tagged_sha256": tsha,
                    "native_sha256": nsha,
                    "cassette_start_1based": str(pos + 1),
                    "cassette_end_1based": str(pos + len(CASSETTE)),
                    "insert_after_native_residue": str(pos),
                    "mean_plddt": str(res["mean_plddt"]),
                    "verdict": res["verdict"],
                    "pdb_key": pdb_key,
                    "scores_key": scores_key,
                    "pdb_sha256": pdb_sha,
                    "colabfold_version": "1.6.2 (ghcr.io/sokrypton/colabfold:1.6.2-cuda12)",
                    "settings": "alphafold2_ptm model_1, num_recycle=3, no relax, no templates",
                    "msa_source": (
                        "native protein's MSA from api.colabfold.com (mmseqs2_uniref_env, filtered),"
                        " cached per native SHA-256; non-native residues as gap columns"
                    ),
                    "implemented_by": IMPLEMENTED_BY,
                }
            )
            s_seen.add(tsha)

    with VERDICTS_TSV.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=v_cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(v_rows)
    with STRUCTURES_TSV.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=s_cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(s_rows)
    print(f"verdicts: {len(v_rows)} rows; structures: {len(s_rows)} rows")


def update_arm_decisions() -> None:
    prots = ld.designed_proteins(HERE / "candidate_manifest_tranches.tsv", library=ld.LIBRARY_L1)
    verdicts = {}
    for r in csv.DictReader(VERDICTS_TSV.open(), delimiter="\t"):
        verdicts[r["tagged_sha256"]] = r

    rows = list(csv.DictReader(ARM_DECISIONS.open(), delimiter="\t"))
    cols = list(rows[0].keys())
    out = []
    dropped = 0
    swapped = 0

    for r in rows:
        aid = r["arm_id"]
        if aid in DROP_ARMS:
            dropped += 1
            continue
        if aid in SWAP_ARMS:
            new_aid, new_pos, new_flag, reason = SWAP_ARMS[aid]
            did = r["design_id"]
            nat = prots[did]["sequence"]
            tag_sha = hashlib.sha256(
                (nat[:new_pos] + CASSETTE + nat[new_pos:]).encode()
            ).hexdigest()
            v = verdicts[tag_sha]
            assert v["verdict"] == "pass", (
                f"expected passing refold for {new_aid}, got {v['verdict']}"
            )
            r["arm_id"] = new_aid
            r["construct_key"] = f"{did}:internal{new_pos}"
            r["insert_after_residue"] = str(new_pos)
            r["refold_verdict"] = v["verdict"]
            r["refold_reasons"] = v["verdict_reasons"]
            r["refold_source"] = v["source"]
            r["refold_disposition"] = "pass"
            r["refold_tagged_sha256"] = tag_sha
            r["veto_class"] = "none"
            r["flag"] = new_flag
            r["reason"] = reason
            r["decided_by"] = DECIDED_BY
            r["implemented_by"] = IMPLEMENTED_BY
            r["decided_on"] = DECIDED_ON
            swapped += 1
        elif aid in CLEAR_FLAGS:
            r["flag"] = ""
        out.append(r)

    # Re-normalize loop_best / loop_2 / loop_3 roles per design_id for affected CDSs
    by_did: dict[str, list[dict]] = {}
    for r in out:
        if r["decision"] == "build" and r["arm_role"].startswith("loop"):
            by_did.setdefault(r["design_id"], []).append(r)
    role_names = ("loop_best", "loop_2", "loop_3")
    for did, loop_rows in by_did.items():
        if any(
            (r["arm_id"] in {v[0] for v in SWAP_ARMS.values()})
            or any(d.startswith(f"{did}-") for d in DROP_ARMS)
            for r in loop_rows
        ):
            for idx, r in enumerate(loop_rows):
                r["arm_role"] = role_names[idx]

    with ARM_DECISIONS.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out)
    n_build = sum(1 for r in out if r["decision"] == "build")
    n_loop = sum(1 for r in out if r["decision"] == "build" and r["arm_role"].startswith("loop"))
    print(
        f"arm_decisions: -{dropped} dropped, {swapped} swapped -> {n_build} built ({n_loop} loop)"
    )


def record_arm_reissues() -> None:
    reissues_tsv = HERE / "l1_arm_reissues.tsv"
    prots = ld.designed_proteins(HERE / "candidate_manifest_tranches.tsv", library=ld.LIBRARY_L1)
    rows = list(csv.DictReader(reissues_tsv.open(), delimiter="\t"))
    cols = list(rows[0].keys())
    have = {r["old_key"] for r in rows}
    added = 0
    for aid in sorted(DROP_ARMS | set(SWAP_ARMS)):
        if aid in have:
            continue
        did = aid.rsplit("-a", 1)[0]
        p_sha = hashlib.sha256(prots[did]["sequence"].encode()).hexdigest()
        reason = (
            SWAP_ARMS[aid][3]
            if aid in SWAP_ARMS
            else "pre-order audit drop: arm removed from built order; encoding stays issued, unused"
        )
        rows.append(
            {
                "old_key": aid,
                "new_key": "",
                "design_id": did,
                "new_protein_sha256": p_sha,
                "reason": reason,
                "decided_by": DECIDED_BY,
                "implemented_by": IMPLEMENTED_BY,
            }
        )
        added += 1
    with reissues_tsv.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(f"arm_reissues: +{added} recorded ({len(rows)} total)")


def main() -> int:
    fold_dir = Path("/tmp/gpi-omega-fix-2026-11-24b/gpi-omega-fix-2026-11-24b")
    register_new_arm_ids()
    register_modal_refolds(fold_dir)
    update_arm_decisions()
    record_arm_reissues()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
