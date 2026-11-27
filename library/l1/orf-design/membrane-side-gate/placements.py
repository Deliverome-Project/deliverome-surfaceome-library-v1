"""Place each member's own model in the membrane; heights on the member's residues.

A member's model is its ``is_chosen`` structure: the AlphaFold DB model when that model is
of exactly the designed protein (or contains it), else our ColabFold fold of the designed
CDS. Model residues are matched to the designed protein by sequence
(:func:`membrane_side.map_to_model`), never by accession or by length.

Placement of record, best first:

1. **TMDET**, when TmAlphaFold has placed the same AlphaFold DB entry: our copy of the
   model is superposed on TmAlphaFold's on residues modelled at pLDDT >= 70, and used in
   TMDET's frame if the superposition RMSD is at most :data:`SAME_MODEL_RMSD` (the
   entries are then the same prediction; TmAlphaFold's date from 2022 and AlphaFold DB
   has since re-issued some models).
2. **Helix axis** (:func:`membrane_side.fit_helix_axis`) otherwise -- every ColabFold
   fold, and AlphaFold DB entries TmAlphaFold does not hold.

The helix-axis heights are computed for every member either way, so the calibration can
measure them against TMDET where both exist.
"""

from __future__ import annotations

import hashlib
import pickle
import re
from pathlib import Path

import numpy as np
import pandas as pd
from gate_inputs import CACHE, chosen_models, model_path, tmalphafold

from deliverome_analysis.surfaceome import membrane_side as ms

PLACEMENT_CACHE = CACHE / "placements_v9.pkl"
#: Largest C-alpha RMSD (A) at which our AlphaFold DB model and TmAlphaFold's count as
#: the same prediction. Same-version pairs superpose at a median 0.26 A (calibration).
#: Hand-set, not derived: decided and implemented by the l1-membrane-side-gate session
#: (PR #291), 2026-11-10, as a conventional same-conformation cutoff. The RMSDs are not
#: bimodal (471 members: 389 <= 1 A, 54 in 1-3 A, 28 > 3 A), and the choice matters for
#: 9 of the 106 shipped sites on members in the 1-3 A range, which take a different
#: verdict on the helix-axis fit (CNGA1 362 among them); the README lists them.
SAME_MODEL_RMSD = 2.0


def _afdb_accession(model_url: str) -> str | None:
    m = re.search(r"AF-(.+?)-F1", str(model_url))
    return m.group(1) if m else None


def tmdet_frame(model: ms.Model, acc: str) -> tuple[np.ndarray, np.ndarray, float, float] | None:
    """(rotation, translation, rmsd, TMDET half-thickness) into TMDET's frame, or None."""
    tm = tmalphafold(acc)
    if tm is None:
        return None
    meta, path = tm
    theirs = ms.read_model(path)
    to_theirs, _ = ms.map_to_model(model.sequence, theirs.sequence)
    pairs = [
        (i, j)
        for i, j in enumerate(to_theirs)
        if j is not None and model.plddt[i] >= ms.FIT_MIN_PLDDT
    ]
    if len(pairs) < 30:
        return None
    rot, t, rmsd = ms.superpose(model.ca[[i for i, _ in pairs]], theirs.ca[[j for _, j in pairs]])
    half = float(meta["MEMBRANE"]["NORMAL"]["_attributes"]["Z"])
    return rot, t, rmsd, half


def _on_protein(values: np.ndarray, rmap: list[int | None]) -> np.ndarray:
    return np.array([np.nan if j is None else values[j] for j in rmap])


def place_protein(protein: str, topology: str, path: Path, afdb_acc: str | None = None) -> dict:
    """Heights (z, A) and pLDDT on ``protein``'s residues, or the reason there are none."""
    model = ms.read_model(path)
    rmap, match = ms.map_to_model(protein, model.sequence)
    mapped = sum(j is not None for j in rmap)
    out = dict(model_match=match, model_mapped=mapped, model_length=len(model.sequence))
    if mapped < 0.9 * len(protein):
        return out | dict(status="model_does_not_cover_protein")
    topo_m = ms.topology_on_model(topology, rmap, len(model.sequence))

    helix = ms.place(model, topo_m)
    if helix is not None:
        out |= dict(z_helix_axis=_on_protein(helix.z, rmap), n_helices=helix.n_helices)

    record = None
    if afdb_acc:
        frame = tmdet_frame(model, afdb_acc)
        if frame is not None:
            rot, t, rmsd, half = frame
            out |= dict(tmdet_rmsd=round(rmsd, 2), tmdet_half_thickness=half)
            if rmsd <= SAME_MODEL_RMSD:
                record = ms.orient(model, ms.placement_from_frame(model, rot, t, "tmdet"), topo_m)
                if helix is not None:
                    cos = abs(float(record.normal @ helix.normal))
                    out["helix_axis_vs_tmdet_deg"] = round(
                        float(np.degrees(np.arccos(min(cos, 1.0)))), 1
                    )
        else:
            out |= dict(tmdet_rmsd=np.nan)
    if record is None:
        record = helix
    if record is None:
        return out | dict(status="fewer_than_two_confident_helices")
    return out | dict(
        status="placed",
        placement_source=record.source,
        z=_on_protein(record.z, rmap),
        plddt=_on_protein(model.plddt, rmap),
        side_source=record.side_source,
        topo_for=record.topo_for,
        topo_against=record.topo_against,
        kr_out=record.kr_out,
        kr_in=record.kr_in,
        positive_inside_agrees=record.positive_inside_agrees,
    )


_EXTRA: dict[str, tuple[Path, str]] | None = None
_PM: dict[str, str] | None = None


def _place_extra_fold(r) -> dict:
    """A member with no chosen model: an AlphaFold DB or 2026-11-10 ColabFold model of
    exactly its protein, if one exists.

    The structure payload is not yet regenerated for the 182 adopted orthologs, so their
    models have no ``is_chosen`` row. Tried in order: the AlphaFold DB model of the
    accession ``proteoform_match.tsv`` records for the member, then the 2026-11-10 folds.
    Either is used only when its sequence is identical to the member's protein. Decided by
    the l1-membrane-side-gate session (PR #291), 2026-11-10, after the reviewer-package
    session pointed out that 13 such members have an identical AlphaFold DB model.
    """
    global _EXTRA, _PM
    from gate_inputs import PANEL_COMMIT, PANEL_DIR, afdb_model, extra_folds, git_tsv

    if _PM is None:
        pm = git_tsv(PANEL_COMMIT, f"{PANEL_DIR}/proteoform_match.tsv", dtype=str)
        _PM = dict(zip(pm["candidate_id"], pm["accession"], strict=False))
    acc = _PM.get(r.candidate_id)
    # The accession's own model, then those of its UniProt isoforms: a member whose protein
    # is an isoform (SLC22A17 is Q8WUG5-2 since the reference fixes) has the canonical's
    # accession recorded. Isoforms added 2026-11-12 by the l1-membrane-side-gate session
    # (PR #291), after the l1-ortholog-features session (PR #298) found SLC22A17's model.
    from gate_inputs import uniprot_isoforms

    tried = [acc, *uniprot_isoforms(acc)] if isinstance(acc, str) and acc else []
    for a in tried:
        path = afdb_model(a)
        if path is None or ms.read_model(path).sequence != r.protein:
            continue
        return place_protein(r.protein, r.topology, path, a) | dict(
            model_source="alphafold_db_not_in_payload",
            model_ref=f"https://alphafold.ebi.ac.uk/files/AF-{a}-F1-model_v6.pdb",
            payload_match="identical_by_sequence",
        )
    if _EXTRA is None:
        _EXTRA = {ms.read_model(path).sequence: (path, key) for key, path in extra_folds().items()}
    hit = _EXTRA.get(r.protein)
    if hit is None:
        return dict(status="no_model")
    path, key = hit
    return place_protein(r.protein, r.topology, path) | dict(
        model_source="colabfold_fold_not_in_payload",
        model_ref=key,
        payload_match="identical_by_sequence",
    )


def member_placements(members: pd.DataFrame, *, refresh: bool = False) -> dict[str, dict]:
    """Placement per ``candidate_id`` (cached on disk; ``refresh`` recomputes)."""
    cache = {}
    if PLACEMENT_CACHE.is_file() and not refresh:
        cache = pickle.loads(PLACEMENT_CACHE.read_bytes())
    chosen = chosen_models()
    out = {}
    for r in members.itertuples():
        # The topology string itself is in the key, so a rerun on another snapshot (or an
        # override) can never be served a placement oriented by a different string.
        key = (
            r.candidate_id,
            r.sequence_sha256,
            r.topology_version,
            hashlib.sha256(r.topology.encode()).hexdigest(),
        )
        if key in cache:
            out[r.candidate_id] = cache[key]
            continue
        if r.candidate_id not in chosen.index:
            res = _place_extra_fold(r)
        else:
            row = chosen.loc[r.candidate_id]
            path = model_path(row)
            acc = _afdb_accession(row["model_url"]) if row["source"] == "alphafold_db" else None
            if path is None:
                res = dict(status="model_not_downloadable")
            else:
                res = place_protein(r.protein, r.topology, path, acc)
            res |= dict(
                model_source=row["source"],
                model_ref=row["model_url"] if row["source"] != "colabfold" else row["r2_key"],
                payload_match=row["proteoform_match"],
            )
        cache[key] = out[r.candidate_id] = res
    PLACEMENT_CACHE.parent.mkdir(parents=True, exist_ok=True)
    PLACEMENT_CACHE.write_bytes(pickle.dumps(cache))
    return out
