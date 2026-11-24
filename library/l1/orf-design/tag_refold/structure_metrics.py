"""Structure comparisons for tagged-vs-native models, with the residue correspondence known.

Every comparison here is between two models of the same native residues (a tagged
construct's native region and the native model), so no sequence or structural alignment
is needed: residue i of one is residue i of the other. That makes three standard scores
cheap and exact:

* :func:`tm_score` -- TM-score (Zhang & Skolnick 2004) under the fixed correspondence,
  maximised over superpositions by the TM-score program's fragment-seed search.
* :func:`lddt_ca` -- C-alpha lDDT (Mariani et al. 2013): superposition-free, so a domain
  that moves on a hinge costs only the contacts across the hinge, not every residue.
* :func:`kabsch_rmsd` -- RMSD after optimal superposition.

Plus per-residue relative solvent accessibility (Shrake-Rupley via Biopython, divided by
Tien et al. 2013 theoretical maxima), heavy-atom contacts, and a φ/ψ helix call.
"""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import numpy as np

#: Tien et al. 2013 (PLoS ONE 8:e80635), theoretical maximum ASA in Å^2.
MAX_ASA = {
    "ALA": 129.0, "ARG": 274.0, "ASN": 195.0, "ASP": 193.0, "CYS": 167.0, "GLN": 225.0,
    "GLU": 223.0, "GLY": 104.0, "HIS": 224.0, "ILE": 197.0, "LEU": 201.0, "LYS": 236.0,
    "MET": 224.0, "PHE": 240.0, "PRO": 159.0, "SER": 155.0, "THR": 172.0, "TRP": 285.0,
    "TYR": 263.0, "VAL": 174.0,
}  # fmt: skip
THREE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
    "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
    "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}  # fmt: skip


class Model:
    """One single-chain model: sequence, C-alpha coordinates, pLDDT (B-factor) and residues."""

    def __init__(self, path: Path):
        from Bio.PDB import PDBParser

        opener = gzip.open if str(path).endswith(".gz") else open
        with opener(path, "rt") as fh:
            st = PDBParser(QUIET=True).get_structure("m", fh)
        chain = next(st[0].get_chains())
        self.structure = st
        self.residues = [r for r in chain if r.id[0] == " " and "CA" in r]
        self.seq = "".join(THREE.get(r.get_resname(), "X") for r in self.residues)
        self.ca = np.array([r["CA"].coord for r in self.residues], dtype=float)
        self.plddt = np.array([r["CA"].bfactor for r in self.residues], dtype=float)

    def rsa(self) -> np.ndarray:
        """Relative solvent accessibility per residue (Shrake-Rupley, residue level)."""
        from Bio.PDB.SASA import ShrakeRupley

        ShrakeRupley().compute(self.structure[0], level="R")
        return np.array(
            [r.sasa / MAX_ASA.get(r.get_resname(), 200.0) for r in self.residues], dtype=float
        )

    def heavy_atoms(self) -> tuple[np.ndarray, np.ndarray]:
        """Coordinates of every heavy atom and the residue index each belongs to."""
        xyz, idx = [], []
        for i, r in enumerate(self.residues):
            for a in r:
                if a.element != "H":
                    xyz.append(a.coord)
                    idx.append(i)
        return np.array(xyz, dtype=float), np.array(idx)

    def helix(self) -> np.ndarray:
        """True where φ is -100..-30 and ψ is -80..-10: the core right-handed alpha basin."""
        from Bio.PDB.Polypeptide import PPBuilder

        out = np.zeros(len(self.residues), dtype=bool)
        pos = {id(r): i for i, r in enumerate(self.residues)}
        for pp in PPBuilder().build_peptides(self.structure[0]):
            for r, (phi, psi) in zip(pp, pp.get_phi_psi_list(), strict=True):
                if phi is None or psi is None or id(r) not in pos:
                    continue
                phi, psi = np.degrees(phi), np.degrees(psi)
                out[pos[id(r)]] = -100 <= phi <= -30 and -80 <= psi <= -10
        return out


def read_scores(path: Path) -> dict:
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        return json.load(fh)


def _kabsch(p: np.ndarray, q: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Rotation and translation that best map p onto q (least squares)."""
    pc, qc = p.mean(0), q.mean(0)
    h = (p - pc).T @ (q - qc)
    u, _, vt = np.linalg.svd(h)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    rot = vt.T @ np.diag([1, 1, d]) @ u.T
    return rot, qc - pc @ rot.T


def kabsch_rmsd(p: np.ndarray, q: np.ndarray) -> float:
    rot, t = _kabsch(p, q)
    return float(np.sqrt(((p @ rot.T + t - q) ** 2).sum(1).mean()))


def tm_score(p: np.ndarray, q: np.ndarray, norm_len: int | None = None) -> float:
    """TM-score of p against reference q under the fixed residue correspondence.

    The superposition is searched as the TM-score program does: seed on contiguous
    fragments of length L, L/2, L/4 ... (>= 4), then iterate, superposing on the residues
    within a growing distance cutoff, keeping the best score seen.
    """
    n = len(p)
    lt = norm_len or n
    if n < 5:
        return float("nan")
    d0 = max(0.5, 1.24 * np.cbrt(max(lt, 19) - 15) - 1.8)
    best = 0.0

    def score(rot, t):
        d = np.sqrt(((p @ rot.T + t - q) ** 2).sum(1))
        return (1 / (1 + (d / d0) ** 2)).sum() / lt, d

    frag = n
    while frag >= 4:
        step = max(1, frag // 2)
        for start in range(0, n - frag + 1, step):
            sel = np.arange(start, start + frag)
            for _ in range(20):
                rot, t = _kabsch(p[sel], q[sel])
                s, d = score(rot, t)
                best = max(best, s)
                for cut in (d0, d0 + 1, d0 + 2, d0 + 4):
                    new = np.where(d < cut)[0]
                    if len(new) >= 3:
                        break
                if len(new) < 3 or np.array_equal(new, sel):
                    break
                sel = new
        frag //= 2
    return float(best)


def lddt_ca(model: np.ndarray, ref: np.ndarray, subset: np.ndarray | None = None) -> np.ndarray:
    """Per-residue C-alpha lDDT of ``model`` against ``ref`` (15 Å inclusion, 0.5/1/2/4 Å).

    ``subset`` (bool mask) restricts both the scored residues and their partners, so
    disordered stretches in the reference neither score nor count as partners.
    Residues with no partner get NaN.
    """
    if subset is None:
        subset = np.ones(len(ref), dtype=bool)
    idx = np.where(subset)[0]
    r, m = ref[idx], model[idx]
    dr = np.sqrt(((r[:, None] - r[None]) ** 2).sum(-1))
    dm = np.sqrt(((m[:, None] - m[None]) ** 2).sum(-1))
    incl = (dr < 15.0) & ~np.eye(len(idx), dtype=bool)
    diff = np.abs(dr - dm)
    kept = sum((diff < t) & incl for t in (0.5, 1.0, 2.0, 4.0)) / 4.0
    n = incl.sum(1)
    out = np.full(len(ref), np.nan)
    with np.errstate(invalid="ignore", divide="ignore"):
        out[idx] = np.where(n > 0, kept.sum(1) / n, np.nan)
    return out


def residue_contacts(
    xyz: np.ndarray, res: np.ndarray, a: np.ndarray, b: np.ndarray, cutoff: float = 4.5
) -> set[tuple[int, int]]:
    """Residue pairs (i in a, j in b) with any heavy atoms within ``cutoff`` Å."""
    from scipy.spatial import cKDTree

    ia = np.isin(res, a)
    ib = np.isin(res, b)
    if not ia.any() or not ib.any():
        return set()
    ta, tb = cKDTree(xyz[ia]), cKDTree(xyz[ib])
    ra, rb = res[ia], res[ib]
    pairs = set()
    for i, js in enumerate(ta.query_ball_tree(tb, cutoff)):
        for j in js:
            pairs.add((int(ra[i]), int(rb[j])))
    return pairs
