"""Recompute the tag-site signals from a DEPOSITED structure and see what changes.

The `surface_loop` lane and EpicTope both read pLDDT, RSA (freesasa) and secondary
structure (pydssp) off an AlphaFold model. `plot_loop_sites_vs_experiment.py` showed the
deposited structures disagree about secondary structure at every site we picked. This runs
the same two functions over the experimental coordinates instead, so the comparison is
method-matched: same pydssp, same freesasa, same Tien-2013 normalisation, different input.

Three things this has to get right, each of which silently corrupts the answer otherwise:

* **Strip the Fabs.** 8vgn is CD20 bound to a therapeutic Fab that covers EL2. Computing
  RSA on the complex buries exactly the surface we are asking about. Only chains matching
  the target accession are kept.
* **Numbering.** SIFTS returns no UniProt mapping for 8vgn, so the offset is derived by
  aligning the modelled chain's own sequence to the UniProt sequence. No offset is assumed.
* **Oligomer state.** CD20 is a constitutive dimer and the AlphaFold model is a monomer, so
  RSA is reported both ways: the partner chain buries real surface that the monomer model
  cannot know about. The difference is the point, not a nuisance.

Needs pydssp and freesasa, which are accessible-surfaceome dependencies, and imports that
package so the methodology is identical rather than merely similar. Run it with that
repository's interpreter, and point --signals-dir at the cached AlphaFold signals the
tag-site pipeline wrote (``data/external/tag_site_tuning/signals`` in that checkout):

    <accessible-surfaceome>/.venv/bin/python experimental_signals.py \
        --signals-dir <accessible-surfaceome>/data/external/tag_site_tuning/signals
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import urllib.request
from pathlib import Path

import freesasa
from accessible_surfaceome.tag_sites.signals import _MAX_ASA, per_residue_ss
from accessible_surfaceome.tag_sites.surface_loop import (
    LOOP_SS,
    MIN_LOOP_LEN,
    RSA_MIN,
)
from Bio.PDB import PDBIO, MMCIFParser, Select
from Bio.PDB.Polypeptide import protein_letters_3to1

HERE = Path(__file__).parent
#: Downloaded coordinates are a cache, so they live outside the experiment folder:
#: upload-experiment-outputs.sh walks every file under it and would publish them into
#: processed:// alongside the real outputs.
WORK = (
    Path(os.environ.get("XDG_CACHE_HOME", Path.home() / ".cache")) / "deliverome-analysis" / "pdb"
)
TARGETS = {
    "MS4A1": {
        "acc": "P11836",
        "pdb": "8vgn",
        "extracellular": [(79, 84), (142, 188)],
        "sites": {79: "P79 (p0032)", 161: "Y161 (p0033)", 171: "N171 (retired)"},
    },
    "SLC2A4": {
        "acc": "P14672",
        "pdb": "7wsm",
        "extracellular": [(46, 81), (133, 142), (193, 201), (309, 323), (375, 384), (439, 445)],
        "sites": {66: "P66 (p0030)", 193: "E193 (retired)", 320: "G320 (p0031)"},
    },
}


class ChainSelect(Select):
    def __init__(self, keep: set[str]):
        self.keep = keep

    def accept_chain(self, chain):
        return chain.id in self.keep

    def accept_residue(self, residue):
        return residue.id[0] == " "  # drop waters/ligands


def uniprot_sequence(acc: str) -> str:
    url = f"https://rest.uniprot.org/uniprotkb/{acc}.json?fields=sequence"
    return json.load(urllib.request.urlopen(url))["sequence"]["value"]


def chain_sequence(chain) -> tuple[str, list[int]]:
    seq, nums = [], []
    for res in chain:
        if res.id[0] != " " or "CA" not in res:
            continue
        try:
            seq.append(protein_letters_3to1[res.get_resname()])
        except KeyError:
            seq.append("X")
        nums.append(res.id[1])
    return "".join(seq), nums


def offset_to_uniprot(seq: str, nums: list[int], ref: str) -> int | None:
    """auth_number -> UniProt number offset, from a unique sequence anchor.

    Derived, never assumed: a wrong offset would silently move every site.
    """
    for start in range(0, max(1, len(seq) - 30)):
        probe = seq[start : start + 30]
        if "X" in probe:
            continue
        hit = ref.find(probe)
        if hit >= 0 and ref.find(probe, hit + 1) == -1:
            return (hit + 1) - nums[start]
    return None


def rsa_for_chain(pdb_path: Path, chain_id: str) -> dict[int, float]:
    """freesasa RSA restricted to one chain, so a complex can be scored chain-aware."""
    st = freesasa.Structure(str(pdb_path))
    areas = freesasa.calc(st).residueAreas()
    out: dict[int, float] = {}
    for ch in areas:
        if ch != chain_id:
            continue
        for resnum, ra in areas[ch].items():
            try:
                r = int(resnum)
            except ValueError:
                continue
            mx = _MAX_ASA.get(getattr(ra, "residueType", ""), None)
            out[r] = (ra.total / mx) if mx else 0.0
    return out


def window_max(rsa: dict[int, float], res: int, half: int = 2) -> float:
    return max((rsa.get(r, 0.0) for r in range(res - half, res + half + 1)), default=0.0)


def loop_len(ss: dict[int, str], res: int) -> int:
    if ss.get(res, "?") not in LOOP_SS:
        return 0
    lo = hi = res
    while ss.get(lo - 1, "?") in LOOP_SS:
        lo -= 1
    while ss.get(hi + 1, "?") in LOOP_SS:
        hi += 1
    return hi - lo + 1


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--signals-dir",
        required=True,
        type=Path,
        help="directory of cached AlphaFold signal JSONs ({SYMBOL}.json.gz)",
    )
    args = ap.parse_args(argv)
    af_signals: Path = args.signals_dir

    WORK.mkdir(parents=True, exist_ok=True)
    records: list[dict] = []
    for gene, spec in TARGETS.items():
        pdb_id, acc = spec["pdb"], spec["acc"]
        cif = WORK / f"{pdb_id}.cif"
        if not cif.exists():
            urllib.request.urlretrieve(f"https://files.rcsb.org/download/{pdb_id.upper()}.cif", cif)
        ref = uniprot_sequence(acc)
        model = MMCIFParser(QUIET=True).get_structure(pdb_id, str(cif))[0]

        hits = {}
        for chain in model:
            seq, nums = chain_sequence(chain)
            if len(seq) < 40:
                continue
            off = offset_to_uniprot(seq, nums, ref)
            if off is not None:
                hits[chain.id] = off
        if not hits:
            print(f"{gene}: no chain matched {acc}")
            continue
        target = sorted(hits)[0]
        off = hits[target]
        print(
            f"\n=== {gene} / {pdb_id} === chains matching {acc}: "
            f"{sorted(hits)} (using {target}); auth->UniProt offset {off:+d}"
        )

        io = PDBIO()
        io.set_structure(model)
        mono = WORK / f"{pdb_id}_{target}_monomer.pdb"
        io.save(str(mono), ChainSelect({target}))
        oligo = WORK / f"{pdb_id}_{acc}_oligomer.pdb"
        io.save(str(oligo), ChainSelect(set(hits)))

        ss_auth = per_residue_ss(str(mono))
        rsa_mono_auth = rsa_for_chain(mono, target)
        rsa_oligo_auth = rsa_for_chain(oligo, target)
        ss = {k + off: v for k, v in ss_auth.items()}
        rsa_mono = {k + off: v for k, v in rsa_mono_auth.items()}
        rsa_oligo = {k + off: v for k, v in rsa_oligo_auth.items()}

        af = json.load(gzip.open(af_signals / f"{gene}.json.gz"))
        af_ss, af_rsa = af["ss"], af["rsa"]

        af_ss_map = {i + 1: c for i, c in enumerate(af_ss)}
        af_rsa_map = {i + 1: v for i, v in enumerate(af_rsa)}

        print(
            f"  {'site':<22} {'AF ss':>6} {'exp ss':>7} {'AF loop':>8} {'exp loop':>9} "
            f"{'AF wRSA':>8} {'exp wRSA':>9}  verdict"
        )
        for pos, label in spec["sites"].items():
            a_loop, e_loop = loop_len(af_ss_map, pos), loop_len(ss, pos)
            a_rsa, e_rsa = window_max(af_rsa_map, pos), window_max(rsa_mono, pos)
            a_ok = af_ss_map.get(pos) in LOOP_SS and a_loop >= MIN_LOOP_LEN and a_rsa >= RSA_MIN
            e_ok = ss.get(pos) in LOOP_SS and e_loop >= MIN_LOOP_LEN and e_rsa >= RSA_MIN
            print(
                f"  {label:<22} {af_ss_map.get(pos, '-'):>6} {ss.get(pos, '-'):>7} "
                f"{a_loop:>8} {e_loop:>9} {a_rsa:>8.2f} {e_rsa:>9.2f}  "
                f"AF {'PASS' if a_ok else 'FAIL'} / EXP {'PASS' if e_ok else 'FAIL'}"
            )
            records.append(
                {
                    "gene": gene,
                    "pdb": pdb_id,
                    "residue": pos,
                    "site": label,
                    "af_ss": af_ss_map.get(pos, ""),
                    "exp_ss": ss.get(pos, ""),
                    "af_loop_len": a_loop,
                    "exp_loop_len": e_loop,
                    "af_window_rsa": f"{a_rsa:.3f}",
                    "exp_window_rsa": f"{e_rsa:.3f}",
                    "exp_oligomer_rsa": f"{window_max(rsa_oligo, pos):.3f}",
                    "af_passes": a_ok,
                    "exp_passes": e_ok,
                }
            )

        # What would the lane propose from the experimental structure?
        for name, S, R in (
            ("AlphaFold", af_ss_map, af_rsa_map),
            ("experimental (monomer)", ss, rsa_mono),
            ("experimental (oligomer)", ss, rsa_oligo),
        ):
            passed = []
            for lo, hi in spec["extracellular"]:
                for p in range(lo, hi + 1):
                    if S.get(p, "?") not in LOOP_SS:
                        continue
                    if loop_len(S, p) < MIN_LOOP_LEN:
                        continue
                    if window_max(R, p) < RSA_MIN:
                        continue
                    passed.append(p)
            print(
                f"  sites passing LOOP_SS + MIN_LOOP_LEN{MIN_LOOP_LEN} + "
                f"RSA{RSA_MIN} from {name}: {passed if passed else 'NONE'}"
            )
            records.append(
                {
                    "gene": gene,
                    "pdb": pdb_id,
                    "residue": "",
                    "site": f"lane:{name}",
                    "af_ss": "",
                    "exp_ss": "",
                    "af_loop_len": "",
                    "exp_loop_len": "",
                    "af_window_rsa": "",
                    "exp_window_rsa": "",
                    "exp_oligomer_rsa": "",
                    "af_passes": "",
                    "exp_passes": ";".join(str(p) for p in passed),
                }
            )

    out = HERE / "experimental_vs_alphafold.tsv"
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(records[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(records)
    print(f"\nwrote {out.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
