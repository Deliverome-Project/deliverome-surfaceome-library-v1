"""Which side of the membrane a residue of a structure model sits on, and how deep.

A tag cassette inserted into a loop is displayed only if that loop is outside the cell.
Topology predictors say which loops are outside by labelling residues, and they can be
wrong about a whole loop: DeepTMHMM v1 merges two helices of every anoctamin, and the
cytoplasmic loop between them, into one long "extracellular" segment. A structure model
of the protein carries independent evidence, because a loop that is really cytoplasmic
sits on the far side of the helix bundle from the real extracellular loops.

This module places a single-chain model (AlphaFold DB or ColabFold) in a planar membrane
and reports, for every model residue, its C-alpha height ``z`` above the mid-plane along
the membrane normal, in angstroms, positive on the extracellular face.

**Placement.** Two sources, best first:

* :func:`placement_from_frame` -- an external placement of the same model, given as the
  rigid transform into a frame whose z axis is the membrane normal and whose origin is
  on the mid-plane. TMDET's placements of AlphaFold DB models (TmAlphaFold) come this way;
  the caller superposes its copy of the model on TMDET's and passes the transform.
* :func:`fit_helix_axis` -- the normal is the mean axis of the predicted TM helices
  (C-alpha atoms modelled at pLDDT >= :data:`FIT_MIN_PLDDT`, helices of at least
  :data:`MIN_HELIX_RESIDUES` such residues, each sign-aligned to the first) and the
  mid-plane passes through their C-alpha centroid. Only helix *positions* are read from
  the topology string, never its sides. Its error against TMDET is measured by the s4e9
  ``membrane-side-gate`` calibration, which also records why a hydrophobic-slab search
  was tried and dropped (it tilted the normal a median 17 degrees from TMDET's).

**Side.** A slab has two faces and its geometry cannot say which is outside. Two
independent signals are computed and both are reported (:func:`orient`):

* the majority of a topology string's ``O`` against ``I`` residues away from the
  mid-plane. A vote over every confidently modelled loop residue, so one wrong loop -- the
  error this module exists to catch -- cannot flip it; a fully inverted prediction can,
  which is what the second signal is for.
* the positive-inside rule (von Heijne): Lys and Arg are enriched on the cytoplasmic
  face, counted at the helix flanks only (:data:`FLANK_RESIDUES`), not in large domains.

Everything is pure computation on coordinates; fetching models is the caller's job.
"""

from __future__ import annotations

import re
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np

#: Half-thickness, in angstroms, of the band either side of the mid-plane that the side
#: signals ignore (the topology vote reads residues beyond a third of it; the K/R count
#: reads flank residues from 5 A inside it to :data:`FLANK_BAND` beyond it). The
#: plasma-membrane hydrocarbon core is ~30 A thick. It shapes the side vote only; the
#: gate's cut-offs are calibrated separately (:class:`GateThresholds`).
HALF_THICKNESS = 15.0

#: Residues below this pLDDT do not define the helix axis or vote on the side.
FIT_MIN_PLDDT = 70.0

#: A helix contributes to the axis only with at least this many confident residues.
MIN_HELIX_RESIDUES = 10

#: Height band beyond each face that the positive-inside count reads, in angstroms.
FLANK_BAND = 15.0

#: Sequence distance from a membrane-core residue within which a Lys/Arg counts as a
#: helix flank for the positive-inside rule. The rule is about the ends of the loops
#: that connect helices, not about large domains on either side (von Heijne 1989 counts
#: loops shorter than ~60 residues); a K/R-rich ectodomain otherwise outvotes the flanks.
FLANK_RESIDUES = 10

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
    "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
    "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V", "MSE": "M",
}  # fmt: skip


@dataclass(frozen=True)
class Model:
    """One chain of a structure model: sequence, C-alpha coordinates, pLDDT."""

    sequence: str
    ca: np.ndarray
    plddt: np.ndarray

    def __post_init__(self) -> None:
        if not (len(self.ca) == len(self.plddt) == len(self.sequence)):
            raise ValueError("sequence, coordinates and pLDDT must have one entry per residue")


def read_model(path: str | Path, chain: str | None = None) -> Model:
    """Read the first (or named) chain of a PDB or mmCIF model.

    pLDDT is read from the C-alpha B-factor, which is where AlphaFold DB and ColabFold
    write it. Residues without a C-alpha are skipped.
    """
    from Bio.PDB import MMCIFParser, PDBParser

    path = Path(path)
    cif = path.suffix in (".cif", ".mmcif")
    parser = MMCIFParser(QUIET=True) if cif else PDBParser(QUIET=True)
    structure = parser.get_structure("m", str(path))
    chains = list(structure[0])
    ch = next(c for c in chains if c.id == chain) if chain else chains[0]
    res = [r for r in ch if r.id[0] == " " and "CA" in r]
    seq = "".join(THREE_TO_ONE.get(r.get_resname(), "X") for r in res)
    ca = np.array([r["CA"].coord for r in res], dtype=float)
    plddt = np.array([r["CA"].get_bfactor() for r in res], dtype=float)
    return Model(seq, ca, plddt)


def map_to_model(query: str, model_seq: str) -> tuple[list[int | None], str]:
    """Index of each ``query`` residue in the model (0-based), and how it was matched.

    ``identical`` and ``contained`` (the query is an exact substring of the model) map
    every residue. Anything else is a global alignment (BLOSUM62, gap -10/-0.5) in which
    only residues aligned to the *same* amino acid are mapped: a residue of a different
    proteoform is not evidence about ours.
    """
    if query == model_seq:
        return list(range(len(query))), "identical"
    at = model_seq.find(query)
    if at >= 0:
        return list(range(at, at + len(query))), "contained"
    from Bio import Align
    from Bio.Align import substitution_matrices

    al = Align.PairwiseAligner()
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    al.open_gap_score, al.extend_gap_score = -10.0, -0.5
    al.mode = "global"
    best = al.align(query, model_seq)[0]
    out: list[int | None] = [None] * len(query)
    for (qa, qb), (ma, _mb) in zip(*best.aligned, strict=True):
        for k in range(qb - qa):
            if query[qa + k] == model_seq[ma + k]:
                out[qa + k] = ma + k
    return out, "aligned"


def topology_on_model(
    topology: str, residue_map: Sequence[int | None], n_model: int
) -> list[str | None]:
    """Project a per-residue topology string of the query onto model residues."""
    out: list[str | None] = [None] * n_model
    for i, j in enumerate(residue_map):
        if j is not None and i < len(topology):
            out[j] = topology[i]
    return out


def segments(topology: str, state: str) -> list[tuple[int, int]]:
    """1-based inclusive (start, end) runs of ``state`` in a topology string."""
    return [(m.start() + 1, m.end()) for m in re.finditer(f"{re.escape(state)}+", topology)]


def superpose(a: np.ndarray, b: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    """Rotation ``r`` and translation ``t`` minimising |a @ r.T + t - b|, and the RMSD."""
    ca, cb = a.mean(0), b.mean(0)
    u, _, vt = np.linalg.svd((a - ca).T @ (b - cb))
    d = np.sign(np.linalg.det(u @ vt))
    r = (u @ np.diag([1.0, 1.0, d]) @ vt).T
    t = cb - ca @ r.T
    return r, t, float(np.sqrt((((a @ r.T + t) - b) ** 2).sum(1).mean()))


@dataclass
class Placement:
    """A membrane placement of one model.

    ``normal`` (model frame) points from the cytoplasmic to the extracellular face once
    :func:`orient` has run; ``z`` is each model residue's C-alpha height above the
    mid-plane along it, in angstroms. ``source`` names the placement method.
    """

    normal: np.ndarray
    z: np.ndarray
    source: str
    n_helices: int = 0
    #: Topology-string residues (confident, away from the mid-plane) agreeing /
    #: disagreeing with the chosen face. Zero when no topology was given.
    topo_for: int = 0
    topo_against: int = 0
    #: Flank Lys+Arg on the face called extracellular / cytoplasmic.
    kr_out: int = 0
    kr_in: int = 0
    #: Which signal set the face: "topology", "positive_inside", or "none".
    side_source: str = "none"

    @property
    def positive_inside_agrees(self) -> bool | None:
        """Whether the K/R count agrees with the face chosen (None on a tie)."""
        if self.kr_in == self.kr_out:
            return None
        return self.kr_in > self.kr_out


def placement_from_frame(
    model: Model, rotation: np.ndarray, translation: np.ndarray, source: str
) -> Placement:
    """Placement from a rigid transform into a membrane frame (normal = z, mid-plane z = 0).

    The sign is not yet chosen; :func:`orient` chooses it.
    """
    z = (model.ca @ rotation.T + translation)[:, 2]
    return Placement(normal=rotation.T @ np.array([0.0, 0.0, 1.0]), z=z, source=source)


def fit_helix_axis(
    model: Model,
    topology_on_model: Sequence[str | None],
    *,
    min_plddt: float = FIT_MIN_PLDDT,
    min_helix_residues: int = MIN_HELIX_RESIDUES,
    min_helices: int = 2,
) -> Placement | None:
    """Membrane normal from the mean TM-helix axis; mid-plane through the helix C-alphas.

    ``None`` when fewer than ``min_helices`` helices have ``min_helix_residues`` residues
    modelled at ``min_plddt``: one helix (a single-pass or GPI protein) does not define a
    plane, and the side vote for such a protein rests on a handful of residues.
    """
    lab = "".join(t or "-" for t in topology_on_model)
    helices, idx_all = [], []
    for m in re.finditer("M+", lab):
        idx = [i for i in range(m.start(), m.end()) if model.plddt[i] >= min_plddt]
        if len(idx) >= min_helix_residues:
            helices.append(model.ca[idx])
            idx_all += idx
    if len(helices) < min_helices:
        return None
    axes = []
    for pts in helices:
        a = pts[-4:].mean(0) - pts[:4].mean(0)
        axes.append(a / np.linalg.norm(a))
    normal = np.mean([a if a @ axes[0] >= 0 else -a for a in axes], axis=0)
    normal /= np.linalg.norm(normal)
    centre = model.ca[idx_all].mean(0)
    z = (model.ca - centre) @ normal
    return Placement(normal=normal, z=z, source="helix_axis", n_helices=len(helices))


def _flank_kr(
    model: Model, z: np.ndarray, half_thickness: float, min_plddt: float = FIT_MIN_PLDDT
) -> tuple[int, int]:
    """Lys/Arg just outside each face (above, below) and near a membrane-core residue.

    Only residues modelled at ``min_plddt`` count, as for the topology vote. A disordered
    tail's coordinates are not a position: DISP2's cytoplasmic termini (pLDDT 26-32) lie
    across the slab in its AlphaFold model and once put 4 of its 10 "outside" K/R there.
    """
    idx = np.flatnonzero(np.abs(z) <= half_thickness - 5.0)
    if not len(idx):
        return 0, 0
    pos = np.arange(len(z))
    near = np.abs(pos[:, None] - idx[None, :]).min(1) <= FLANK_RESIDUES
    kr = np.array([a in "KR" for a in model.sequence])
    band = (np.abs(z) >= half_thickness - 5.0) & (np.abs(z) <= half_thickness + FLANK_BAND)
    sel = kr & near & band & (model.plddt >= min_plddt)
    return int((sel & (z > 0)).sum()), int((sel & (z < 0)).sum())


def orient(
    model: Model,
    placement: Placement,
    topology_on_model: Sequence[str | None] | None = None,
    *,
    half_thickness: float = HALF_THICKNESS,
    min_plddt: float = FIT_MIN_PLDDT,
) -> Placement:
    """Choose which face is extracellular and flip ``z`` so it is positive there.

    The topology vote decides when there is one and it is not tied; the positive-inside
    count is recorded either way, and decides only when there is no topology vote.
    Callers that need both to agree check :func:`orientation_supported`.
    """
    z = placement.z
    conf = model.plddt >= min_plddt
    kr_up, kr_dn = _flank_kr(model, z, half_thickness, min_plddt)
    t_for = t_against = 0
    flip = None
    source = "none"
    if topology_on_model is not None:
        lab = np.array([t if t else "-" for t in topology_on_model])
        away = conf & (np.abs(z) > half_thickness / 3)
        up = ((lab == "O") & (z > 0)) | ((lab == "I") & (z < 0))
        dn = ((lab == "O") & (z < 0)) | ((lab == "I") & (z > 0))
        up_votes, dn_votes = int((up & away).sum()), int((dn & away).sum())
        if up_votes != dn_votes:
            flip = dn_votes > up_votes
            t_for, t_against = max(up_votes, dn_votes), min(up_votes, dn_votes)
            source = "topology"
    if flip is None and kr_up != kr_dn:
        flip = kr_up > kr_dn  # more K/R on the + face means + is cytoplasmic
        source = "positive_inside"
    sign = -1.0 if flip else 1.0
    placement.normal = placement.normal * sign
    placement.z = placement.z * sign
    placement.kr_out, placement.kr_in = (kr_dn, kr_up) if flip else (kr_up, kr_dn)
    placement.topo_for, placement.topo_against = t_for, t_against
    placement.side_source = source
    return placement


def place(model: Model, topology_on_model: Sequence[str | None], **kw) -> Placement | None:
    """:func:`fit_helix_axis` then :func:`orient`."""
    p = fit_helix_axis(model, topology_on_model, **kw)
    return None if p is None else orient(model, p, topology_on_model)


# --- the gate -----------------------------------------------------------------------


@dataclass(frozen=True)
class GateThresholds:
    """Cut-offs for :func:`site_verdict`. No defaults: a caller derives them from evidence.

    ``outside`` -- junction height (A) at or above which a site passes.
    ``margin`` -- below ``outside`` by less than this, the placement error could still put
    the junction outside: reported for review, not failed.
    ``core`` -- hydrophobic half-thickness (A). A failing junction below ``-core`` is
    labelled cytoplasmic; one above it, in the membrane. The label only describes the
    failure.
    ``min_plddt`` -- below this pLDDT a failing verdict is reported as a hold, because the
    model's coordinates for the junction are not evidence enough to drop a pick.
    ``min_topology_ratio`` -- the topology vote for the chosen face must beat the vote
    against it by this factor.
    ``max_kr_excess_outside`` -- more flank Lys/Arg on the chosen extracellular face than
    on the cytoplasmic one, by more than this, contradicts the face.
    """

    outside: float
    margin: float
    core: float
    min_plddt: float
    min_topology_ratio: float
    max_kr_excess_outside: int


def orientation_supported(
    topo_for: int, topo_against: int, kr_out: int, kr_in: int, th: GateThresholds
) -> bool:
    """Whether the face called extracellular is backed by both side signals."""
    if topo_for < th.min_topology_ratio * max(topo_against, 1):
        return False
    return (kr_out - kr_in) <= th.max_kr_excess_outside


#: Verdicts :func:`site_verdict` returns.
VERDICTS = (
    "pass",
    "review_interface",
    "review_orientation",
    "hold_low_plddt",
    "fail_membrane",
    "fail_cytoplasmic",
    "unscorable",
)
#: A caller's manual exclusion (for example a partner-assembly interface) is reported as
#: ``fail_partner_interface``; :func:`site_verdict` never returns it.


def site_verdict(z: float, plddt: float, orientation_ok: bool, th: GateThresholds) -> str:
    """Gate verdict for one insertion junction at height ``z`` (A) with pLDDT ``plddt``."""
    if np.isnan(z):
        return "unscorable"
    if not orientation_ok:
        return "review_orientation"
    if z >= th.outside:
        return "pass"
    if z >= th.outside - th.margin:
        return "review_interface"
    if plddt < th.min_plddt:
        return "hold_low_plddt"
    return "fail_cytoplasmic" if z < -th.core else "fail_membrane"
