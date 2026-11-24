"""Run DeepTMHMM **v1** on Modal, for the shortlisted alternate forms.

v1 rather than v2 deliberately: every reference topology call in this project is v1, and a
TM-helix count is only meaningful against another count from the same predictor. v2 would
also have to re-run the references to be comparable, and its signal-peptide recall is worse
(92.9% against v1's 99.5%).

The package is license-gated (CC BY-NC-SA 4.0) and is **not** committed. It is uploaded once
to a Modal volume from the copy in ``deliverome-internal``::

    modal volume create deeptmhmm-v1
    modal volume put deeptmhmm-v1 \\
      ~/Git/deliverome-internal/analyses/surface-proteome/data/external/deeptmhmm/DeepTMHMM-Academic-License-v1.0 \\
      /DeepTMHMM-Academic-License-v1.0

Dependencies mirror ``deeptmhmm_uv_requirements.txt``, which is the pinned modern-stack set
this project already uses locally; the upstream README's torch 1.5.0+cu92 is a legacy wheel
that no longer installs.

    modal run modal/deeptmhmm_v1_app.py::canary --n 8
    modal run modal/deeptmhmm_v1_app.py::sweep --fasta ../topology_shortlist.fasta
"""

from __future__ import annotations

import hashlib
import json
import re
from functools import lru_cache
import time
from pathlib import Path

import modal

# A path the base image definitely does not already populate: Modal refuses to
# mount a volume over a non-empty directory, and /pkg collided.
MOUNT = "/deeptmhmm"
PACKAGE = f"{MOUNT}/DeepTMHMM-Academic-License-v1.0"
# DeepTMHMM holds every sequence's ESM embedding in memory, so a shard is kept small enough
# that a long protein cannot push the worker into swap.
SHARD_RESIDUES = 60_000
LONG_SEQUENCE = 3_000

image = (
    modal.Image.debian_slim(python_version="3.11")
    .apt_install("libhdf5-dev", "build-essential")
    .pip_install(
        "biopython==1.79",
        "fair-esm==0.4.0",
        "h5py==3.12.1",
        "PeptideBuilder==1.1.0",
        "matplotlib==3.9.2",
        "numpy==1.26.4",
        "requests==2.32.5",
        "tqdm==4.67.1",
        "torch==2.5.1",
    )
)
if modal.is_local():  # the image is built from the laptop; containers only re-import this file
    from deliverome_analysis.trust_store import harden

    image = harden(image)  # no compromised eMudhra roots: docs/security/emudhra-trust-store-audit.md

app = modal.App("surfaceome-deeptmhmm-v1")
package_volume = modal.Volume.from_name("deeptmhmm-v1", create_if_missing=True)
out_volume = modal.Volume.from_name("deeptmhmm-v1-runs", create_if_missing=True)


def _read_fasta(text: str) -> list[tuple[str, str]]:
    out, name, buf = [], None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if name:
                out.append((name, "".join(buf)))
            name, buf = line[1:].strip(), []
        elif name:
            buf.append(line.strip())
    if name:
        out.append((name, "".join(buf)))
    return out


def _shard(records: list[tuple[str, str]], budget: int) -> list[list[tuple[str, str]]]:
    shards, cur, total = [], [], 0
    for rec in records:
        if cur and total + len(rec[1]) > budget:
            shards.append(cur)
            cur, total = [], 0
        cur.append(rec)
        total += len(rec[1])
    if cur:
        shards.append(cur)
    return shards


def _sha_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


@lru_cache(maxsize=1)
def _runtime_identity() -> dict:
    """Observe the mounted academic package; never infer a patch release."""
    package = Path(PACKAGE)
    predict = package / "predict.py"
    version = re.search(r'ArgumentParser\(description=[\'"]DeepTMHMM v([^\'"]+)[\'"]', predict.read_text())
    if not version or version.group(1) != "1.0":
        raise ValueError("Mounted predict.py does not declare the expected DeepTMHMM v1.0")
    models = [package / f"deeptmhmm_cv_{i}.model" for i in range(5)]
    checkpoints = sorted(package.rglob("*.pt"))
    if not checkpoints:
        raise ValueError("Mounted package has no ESM checkpoint to fingerprint")
    files = sorted(set(package.rglob("*.py")) | set(models) | set(checkpoints))
    fingerprints = {str(path.relative_to(package)): _sha_file(path) for path in files}
    return {
        "tool_version": "deeptmhmm-1.0-academic",
        "model_version": version.group(1),
        "version_source": 'predict.py ArgumentParser description="DeepTMHMM v1.0"',
        "runner_sha256": _sha_file(Path(__file__)),
        "package_predict_sha256": fingerprints["predict.py"],
        "package_files_sha256": fingerprints,
        "model_files": [str(path.relative_to(package)) for path in models],
        "esm_checkpoint_files": [str(path.relative_to(package)) for path in checkpoints],
    }


def _run_predict(records: list[tuple[str, str]], tag: str) -> dict:
    """Invoke the vendored predict.py and return its 3-line topology output."""
    import subprocess
    import sys
    import tempfile

    started = time.time()
    identity = _runtime_identity()
    with tempfile.TemporaryDirectory() as tmp:
        fasta = Path(tmp) / "in.fasta"
        fasta.write_text("".join(f">{n}\n{s}\n" for n, s in records))
        outdir = Path(tmp) / "out"  # predict.py refuses a directory that already exists
        proc = subprocess.run(
            [sys.executable, "predict.py", "--fasta", str(fasta), "--output-dir", str(outdir)],
            cwd=PACKAGE,
            capture_output=True,
            text=True,
        )
        three = outdir / "predicted_topologies.3line"
        if proc.returncode != 0 or not three.exists():
            return {
                "tag": tag,
                "n": len(records),
                "ok": False,
                "error": (proc.stderr or proc.stdout)[-2000:],
            }
        return {
            "tag": tag,
            "n": len(records),
            "ok": True,
            "seconds": round(time.time() - started, 1),
            "three_line": three.read_text(),
            "runtime_identity": identity,
            "submitted_fasta_sha256": _sha_file(fasta),
            "output_sha256": _sha_file(three),
            "sequence_sha256s": [hashlib.sha256(seq.encode()).hexdigest() for _, seq in records],
        }


@app.function(
    image=image,
    gpu="T4",
    timeout=60 * 60,
    volumes={MOUNT: package_volume},
)
def predict_shard(payload: dict) -> dict:
    return _run_predict([tuple(r) for r in payload["records"]], payload["tag"])


@app.function(
    image=image,
    gpu="A10G",
    timeout=120 * 60,
    volumes={MOUNT: package_volume},
)
def predict_long(payload: dict) -> dict:
    """Long sequences get more memory; v1's ESM step is what OOMs, not the CRF."""
    return _run_predict([tuple(r) for r in payload["records"]], payload["tag"])


@app.local_entrypoint()
def canary(n: int = 8, fasta: str = "../topology_shortlist.fasta") -> None:
    """Prove the image and the volume before paying for the full sweep."""
    records = _read_fasta(Path(fasta).read_text())[:n]
    res = predict_shard.remote({"records": records, "tag": "canary"})
    print(json.dumps({k: v for k, v in res.items() if k != "three_line"}, indent=2))
    if res.get("ok"):
        print(res["three_line"][:600])


@app.local_entrypoint()
def sweep(fasta: str = "../topology_shortlist.fasta", out: str = "../topology_shortlist.3line") -> None:
    if any(Path(out + suffix).exists() for suffix in ("", ".runtime.json", ".failures.json")):
        raise FileExistsError("Prediction outputs must use fresh paths")
    records = _read_fasta(Path(fasta).read_text())
    short = [r for r in records if len(r[1]) <= LONG_SEQUENCE]
    long_ = [r for r in records if len(r[1]) > LONG_SEQUENCE]
    shards = [{"records": s, "tag": f"s{i}"} for i, s in enumerate(_shard(short, SHARD_RESIDUES))]
    longs = [{"records": [r], "tag": f"L{i}"} for i, r in enumerate(long_)]
    print(f"{len(records)} sequences -> {len(shards)} shards + {len(longs)} long singletons")

    results = list(predict_shard.map(shards)) + (
        list(predict_long.map(longs)) if longs else []
    )
    ok = [r for r in results if r.get("ok")]
    bad = [r for r in results if not r.get("ok")]
    Path(out).write_text("".join(r["three_line"] for r in ok))
    Path(out + ".failures.json").write_text(json.dumps(bad, indent=2))
    identities = [r["runtime_identity"] for r in ok]
    if identities and any(identity != identities[0] for identity in identities[1:]):
        raise ValueError("Prediction workers used different package identities")
    runtime = dict(identities[0]) if identities else {}
    runtime.update({
        "schema_version": "deeptmhmm-academic-sweep-1",
        "submitted_fasta_sha256": _sha_file(Path(fasta)),
        "output_sha256": _sha_file(Path(out)),
        "shards": [{k: v for k, v in r.items() if k != "three_line"} for r in results],
        "failed_shards": len(bad),
    })
    Path(out + ".runtime.json").write_text(json.dumps(runtime, indent=2, sort_keys=True) + "\n")
    print(f"wrote {out}: {len(ok)} shards ok, {len(bad)} failed")
