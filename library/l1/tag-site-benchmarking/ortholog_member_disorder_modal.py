"""AlphaFold-disorder (pLDDT, windowed RSA) and metapredict for the 27 adopted orthologs.

The recipe is accessible-surfaceome's `modal/disorder_app.py` (`alphafold_disorder_shard`,
`metapredict_shard`), which wrote every `disorder_public` row the tag-site model trained on,
with one change: the model is passed in, not fetched from AlphaFold DB by accession,
because four of these members have no AlphaFold DB model and are our ColabFold folds.
Everything else is the same: Debian `dssp` through Biopython, pLDDT from the CA B-factor
(byte = round(pLDDT * 2.55)), RSA as DSSP relative ASA meaned over +/-12 residues
(byte = round(RSA * 255)), metapredict 3.0.2 disorder (byte = round(score * 255)).

Inputs come from `ortholog_member_features.py prepare`; results come back as bytes and are
written locally. Nothing is uploaded from the container.

    uv run --locked --group modal modal run ortholog_member_disorder_modal.py
"""

from __future__ import annotations

import json
from pathlib import Path

import modal

HERE = Path(__file__).resolve().parent
WORK = HERE / "data" / "ortholog_members"

metapredict_image = (
    modal.Image.debian_slim(python_version="3.11")
    .pip_install("metapredict==3.0.2", "numpy<2")
    .run_commands(
        "python -c \"import metapredict as m; m.predict_disorder('MKTIIALSYIFCLVFADYKDDDDK')\""
    )
)
afd_image = (
    modal.Image.debian_slim(python_version="3.11")
    .apt_install("dssp")
    .pip_install("biopython", "numpy<2")
)
if modal.is_local():  # the image is built from the laptop; containers only re-import this file
    from deliverome_analysis.core.trust_store import harden

    metapredict_image = harden(metapredict_image)  # no compromised eMudhra roots
    afd_image = harden(afd_image)

app = modal.App("s5e13-ortholog-member-disorder")

SUBSTITUTIONS = {"U": "C", "O": "K", "B": "N", "Z": "Q", "J": "L", "X": "A"}


@app.function(image=metapredict_image, cpu=2.0, timeout=1800)
def metapredict_one(rec: dict) -> dict:
    import metapredict as meta

    seq = "".join(SUBSTITUTIONS.get(c, c) for c in rec["sequence"])
    scores = meta.predict_disorder(seq)
    assert len(scores) == len(seq)
    return {
        "key": rec["key"],
        "version": getattr(meta, "__version__", "3.0.2"),
        "disorder_hex": bytes(min(255, max(0, round(s * 255))) for s in scores).hex(),
    }


@app.function(image=afd_image, cpu=2.0, timeout=1800)
def alphafold_disorder_one(rec: dict) -> dict:
    import subprocess

    from Bio.PDB import PDBParser
    from Bio.PDB.DSSP import DSSP

    pdb = Path("/tmp") / f"{rec['key']}.pdb"
    # Same DBREF strip as the reference: a 10-character accession overflows the field.
    lines = [ln for ln in rec["pdb"].splitlines() if not ln.startswith("DBREF")]
    # ColabFold writes no HEADER or CRYST1 record, and mkdssp 4.2 reads a file that does
    # not start with HEADER as mmCIF and fails. AlphaFold DB files carry both (the P 1 unit
    # cell below); adding the same two lines to a fold changes no coordinate, it only lets
    # DSSP read the file.
    if not any(ln.startswith("CRYST1") for ln in lines):
        cryst = "CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1           1          "
        lines.insert(0, cryst)
    if not lines[0].startswith("HEADER"):
        lines.insert(0, "HEADER".ljust(80))
    pdb.write_text("\n".join(lines) + "\n")
    model = next(iter(PDBParser(QUIET=True).get_structure(rec["key"], str(pdb))))
    try:
        dssp = DSSP(model, str(pdb), dssp="mkdssp")
    except Exception as exc:
        probe = subprocess.run(
            ["mkdssp", str(pdb), "/tmp/probe.dssp"], capture_output=True, text=True
        )
        return {"key": rec["key"], "error": f"{exc}; mkdssp: {probe.stderr[-600:]}"}
    plddt, rsa = {}, {}
    for key in dssp.keys():
        rsa[key[1][1]] = float(dssp[key][3])
    for res in model.get_residues():
        if "CA" in res:
            plddt[res.id[1]] = float(res["CA"].get_bfactor())
    n = max(plddt) if plddt else 0
    win = 12
    smoothed = []
    for i in range(1, n + 1):
        vals = [rsa[j] for j in range(max(1, i - win), min(n, i + win) + 1) if j in rsa]
        smoothed.append(sum(vals) / len(vals) if vals else 0.0)
    ver = subprocess.run(["mkdssp", "--version"], capture_output=True, text=True)
    return {
        "key": rec["key"],
        "model_residues": n,
        "dssp_version": (ver.stdout or ver.stderr).strip().splitlines()[0] if ver else "",
        "plddt_hex": bytes(
            min(255, max(0, round(plddt.get(i, 0.0) * 2.55))) for i in range(1, n + 1)
        ).hex(),
        "rsa_hex": bytes(min(255, max(0, round(v * 255))) for v in smoothed).hex(),
    }


@app.local_entrypoint()
def main(inputs: str = str(WORK / "modal_inputs.json"), out: str = str(WORK / "modal_out.json")):
    recs = json.loads(Path(inputs).read_text())
    with_model = [r for r in recs if r["pdb"]]
    afd = list(alphafold_disorder_one.map(with_model))
    mp = list(metapredict_one.map(recs))
    Path(out).write_text(json.dumps({"alphafold-disorder": afd, "metapredict": mp}))
    print(f"alphafold-disorder {len(afd)}, metapredict {len(mp)} -> {out}")
