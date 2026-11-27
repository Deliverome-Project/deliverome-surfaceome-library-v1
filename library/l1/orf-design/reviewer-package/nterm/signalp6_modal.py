"""Run SignalP 6.0 (slow-sequential) on Modal with every version and file pinned and recorded.

Made so that a SignalP call this package depends on can be reproduced, to the same standard as
``library/l1/proteoform-selection/modal/deeptmhmm_v1_app.py``: the licensed package sits on a
Modal volume and is fingerprinted file by file at run time, and every run writes a
``provenance.json`` naming exactly what produced it.

**What is pinned, and where.**

| Piece | Where | Check made on every run |
|---|---|---|
| SignalP 6.0 package code, version ``6.0+h`` | volume ``signalp6-package``: ``/signalp-6-package/`` | SHA-256 of every file must equal the DTU tarball's member (``/tar_member_sha256.json``) |
| The DTU tarball it came from | ``signalp-6.0i.slow_sequential.tar.gz``, SHA-256 ``3193855d…ba024`` (``/TARBALL.json``) | recorded |
| Model checkpoints, as shipped (CPU) | volume ``signalp6-models``: ``/cpu/sequential_models_signalp6/`` | each SHA-256 must equal the tarball's |
| Checkpoints converted for GPU (used to predict) | volume ``signalp6-models``: ``/gpu/sequential_models_signalp6/`` | SHA-256 recorded (conversion is deterministic but not re-checked against the tarball) |
| Container | debian-slim, Python 3.10, the pins in ``PIP_PINS`` | ``pip freeze``, torch and CUDA versions and the GPU name recorded |
| Invocation | ``--mode slow-sequential --organism eukarya --bsize 64 --format none`` on a T4 | recorded, with the input FASTA SHA-256 and the exact sequence order of every shard |

SignalP's numbers depend on batch composition (the record measured shifts up to 1.0e-4 in
partial batches), so shards keep the input order and are a multiple of 64 long, and the order
is recorded. ``validate`` reruns one shard of the 2026-10-19 sweep ``sweep_rp_1792381414`` in
its original order and requires identical calls.

**Attribution** (convention of 2026-11-10). Decided by Becca, 2026-11-10 ("make sure exact
version sha etc are clearly saved and noted as we did for deeptmhmm"). Implemented by the
Claude Code session ``nterm-reviewer-package`` (PR #254). The package volume was uploaded from
Becca's licensed tarball on 2026-11-10; it cannot be rebuilt from this repo, so its hashes are
the record.

    uv run --locked --group modal modal run signalp6_modal.py::fingerprint
    uv run --locked --group modal modal run signalp6_modal.py::validate
    uv run --locked --group modal modal run signalp6_modal.py::run --fasta in.fasta --run-id ID --out DIR
"""

from __future__ import annotations

import hashlib
import json
import subprocess
import time
from pathlib import Path

import modal

PKG_MOUNT = "/sp6pkg"
MODELS_MOUNT = "/sp6models"
RUNS_MOUNT = "/sp6runs"
PACKAGE_DIR = f"{PKG_MOUNT}/signalp-6-package"
MODE, ORGANISM, BATCH_SIZE = "slow-sequential", "eukarya", 64
SHARD_SIZE = 1024  # a multiple of BATCH_SIZE, as the sweep used
GPU = "T4"
PIP_PINS = ("numpy==1.26.4", "matplotlib==3.8.4", "tqdm==4.66.5")
TORCH_PIN = "torch==1.13.1+cu117"
TORCH_INDEX = "https://download.pytorch.org/whl/cu117"

image = (
    modal.Image.debian_slim(python_version="3.10")
    .pip_install(TORCH_PIN, extra_index_url=TORCH_INDEX)
    .pip_install(*PIP_PINS)
)
if modal.is_local():  # the image is built from the laptop; containers only re-import this file
    from deliverome_analysis.core.trust_store import harden

    image = harden(image)

app = modal.App("deliverome-signalp6-pinned")
package_volume = modal.Volume.from_name("signalp6-package").read_only()
models_volume = modal.Volume.from_name("signalp6-models").read_only()
runs_volume = modal.Volume.from_name("signalp6-pinned-runs", create_if_missing=True)
VOLUMES = {PKG_MOUNT: package_volume, MODELS_MOUNT: models_volume, RUNS_MOUNT: runs_volume}


def _sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1 << 23), b""):
            h.update(block)
    return h.hexdigest()


def _identity() -> dict:
    """Fingerprint the mounted package and checkpoints; refuse anything not from the tarball."""
    tar = json.loads(Path(PKG_MOUNT, "tar_member_sha256.json").read_text())
    by_tail = {k.split("signalp-6-package/", 1)[1]: v["sha256"] for k, v in tar.items()}
    pkg = Path(PACKAGE_DIR)
    package_files = {
        str(p.relative_to(pkg)): _sha(p)
        for p in sorted(pkg.rglob("*"))
        if p.is_file() and "__pycache__" not in p.parts
    }
    bad = {k: v for k, v in package_files.items() if by_tail.get(k) != v}
    if bad:
        raise RuntimeError(f"package files differ from the DTU tarball: {sorted(bad)}")
    ckpt = {}
    for kind in ("cpu", "gpu"):
        d = Path(MODELS_MOUNT, kind, "sequential_models_signalp6")
        ckpt[kind] = {p.name: _sha(p) for p in sorted(d.glob("*.pt"))}
    want = {k.rsplit("/", 1)[1]: v for k, v in by_tail.items() if k.endswith(".pt")}
    if ckpt["cpu"] != want:
        raise RuntimeError("CPU checkpoints on signalp6-models differ from the DTU tarball")
    return {
        "tool": "SignalP 6.0",
        "package_version": "6.0+h",
        "tarball": json.loads(Path(PKG_MOUNT, "TARBALL.json").read_text()),
        "package_files_sha256": package_files,
        "package_files_match_tarball": True,
        "checkpoints_cpu_sha256": ckpt["cpu"],
        "checkpoints_cpu_match_tarball": True,
        "checkpoints_gpu_sha256": ckpt["gpu"],
        "volumes": {"package": "signalp6-package", "models": "signalp6-models"},
    }


def _environment() -> dict:
    import platform

    import torch

    gpu = subprocess.run(
        ["nvidia-smi", "--query-gpu=name,driver_version", "--format=csv,noheader"],
        capture_output=True,
        text=True,
        check=False,
    ).stdout.strip()
    freeze = subprocess.run(["pip", "freeze"], capture_output=True, text=True, check=False).stdout
    return {
        "python": platform.python_version(),
        "torch": str(torch.__version__),
        "cuda": str(torch.version.cuda),
        "cudnn": int(torch.backends.cudnn.version() or 0),
        "gpu": gpu,
        "pip_freeze": sorted(freeze.split()),
    }


@app.function(image=image, volumes=VOLUMES, timeout=1800)
def fingerprint_remote() -> dict:
    return _identity()


@app.local_entrypoint()
def fingerprint(out: str = "signalp6_identity.json"):
    """Check the package and checkpoints against the DTU tarball and save their hashes."""
    ident = fingerprint_remote.remote()
    Path(out).write_text(json.dumps(ident, indent=1, sort_keys=True) + "\n")
    print(json.dumps({k: v for k, v in ident.items() if "sha256" not in k}, indent=1))
    print(
        f"package files: {len(ident['package_files_sha256'])}; "
        f"CPU checkpoints match tarball: {ident['checkpoints_cpu_match_tarball']}; written {out}"
    )


@app.function(image=image, gpu=GPU, volumes=VOLUMES, timeout=3 * 3600, retries=0)
def predict_shard(records: list[tuple[str, str]], run_id: str, shard_id: str) -> dict:
    """SignalP 6.0 on one shard, in the order given. Returns the raw outputs and the provenance."""
    import os
    import tempfile

    started = time.time()
    identity = _identity()
    env = dict(os.environ, PYTHONPATH=PACKAGE_DIR, PYTHONDONTWRITEBYTECODE="1")
    with tempfile.TemporaryDirectory() as tmp:
        fasta = Path(tmp) / "in.fasta"
        fasta.write_text(
            "".join(
                f">{n}\n" + "\n".join(s[i : i + 60] for i in range(0, len(s), 60)) + "\n"
                for n, s in records
            )
        )
        out = Path(tmp) / "out"
        cmd = [
            "python",
            "-c",
            # the package's own console script, signalp6 = signalp:predict
            "import sys, signalp; sys.argv[0] = 'signalp6'; signalp.predict()",
            "--fastafile",
            str(fasta),
            "--output_dir",
            str(out),
            "--model_dir",
            f"{MODELS_MOUNT}/gpu",
            "--mode",
            MODE,
            "--organism",
            ORGANISM,
            "--format",
            "none",
            "--bsize",
            str(BATCH_SIZE),
            "--write_procs",
            "1",
        ]
        proc = subprocess.run(cmd, capture_output=True, text=True, env=env, check=False)
        if proc.returncode:
            raise RuntimeError(proc.stderr[-3000:])
        files = {p.name: p.read_text() for p in out.iterdir() if p.is_file()}
        fasta_sha = _sha(fasta)
    dest = Path(RUNS_MOUNT, run_id, shard_id)
    dest.mkdir(parents=True, exist_ok=True)
    for name, text in files.items():
        (dest / name).write_text(text)
    prov = {
        "run_id": run_id,
        "shard": shard_id,
        "n": len(records),
        "ids_in_order": [n for n, _ in records],
        "input_fasta_sha256": fasta_sha,
        "command": [*cmd[:3], "--fastafile", "in.fasta", "--output_dir", "out", *cmd[7:]],
        "mode": MODE,
        "organism": ORGANISM,
        "batch_size": BATCH_SIZE,
        "gpu_type": GPU,
        "seconds": round(time.time() - started, 1),
        "identity": identity,
        "environment": _environment(),
    }
    (dest / "provenance.json").write_text(json.dumps(prov, indent=1, sort_keys=True))
    runs_volume.commit()
    return {"files": files, "provenance": prov}


def parse_results(text: str) -> dict[str, dict]:
    """``prediction_results.txt`` -> {id: call}; cleavage_site is the last SP residue."""
    import re

    out = {}
    for line in text.splitlines():
        if line.startswith("#") or not line.strip():
            continue
        c = line.split("\t")
        m = re.search(r"CS pos: (\d+)-\d+\. Pr: ([\d.]+)", line)
        out[c[0].split()[0]] = {
            "prediction": c[1],
            "other_probability": float(c[2]),
            "sp_probability": float(c[3]),
            "cleavage_site": int(m.group(1)) if m else None,
            "cleavage_probability": float(m.group(2)) if m else None,
        }
    return out


def _read_fasta(path: str) -> list[tuple[str, str]]:
    recs, name, buf = [], None, []
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            if name:
                recs.append((name, "".join(buf)))
            name, buf = line[1:].split()[0], []
        elif line.strip():
            buf.append(line.strip())
    if name:
        recs.append((name, "".join(buf)))
    return recs


def _local_provenance() -> dict:
    here = Path(__file__).resolve()

    def git(*a: str) -> str:
        return subprocess.run(
            ["git", *a], capture_output=True, text=True, cwd=here.parent
        ).stdout.strip()

    return {
        "runner": str(here.relative_to(Path(git("rev-parse", "--show-toplevel")))),
        "runner_sha256": _sha(here),
        "git_commit": git("rev-parse", "HEAD"),
        "git_dirty": bool(git("status", "--porcelain", "--", str(here))),
        "pip_pins": [TORCH_PIN, *PIP_PINS],
        "torch_index": TORCH_INDEX,
    }


def _run(records: list[tuple[str, str]], run_id: str, out: Path) -> dict[str, dict]:
    out.mkdir(parents=True, exist_ok=True)
    shards = [records[i : i + SHARD_SIZE] for i in range(0, len(records), SHARD_SIZE)]
    calls, provs = {}, []
    args = [(s, run_id, f"s{i:04d}") for i, s in enumerate(shards)]
    for i, res in enumerate(predict_shard.starmap(args)):
        sd = out / f"s{i:04d}"
        sd.mkdir(exist_ok=True)
        for name, text in res["files"].items():
            (sd / name).write_text(text)
        (sd / "provenance.json").write_text(json.dumps(res["provenance"], indent=1, sort_keys=True))
        calls |= parse_results(res["files"]["prediction_results.txt"])
        provs.append(res["provenance"])
    summary = {
        "run_id": run_id,
        "n": len(records),
        "shards": len(shards),
        "local": _local_provenance(),
        "identity": provs[0]["identity"],
        "environment": provs[0]["environment"],
        "seconds_per_shard": [p["seconds"] for p in provs],
    }
    (out / "provenance.json").write_text(json.dumps(summary, indent=1, sort_keys=True))
    (out / "calls.json").write_text(json.dumps(calls, indent=1, sort_keys=True))
    return calls


@app.local_entrypoint()
def run(fasta: str, run_id: str, out: str):
    calls = _run(_read_fasta(fasta), run_id, Path(out))
    print(f"{len(calls)} calls written to {out}/calls.json")


@app.local_entrypoint()
def validate(out: str = "signalp6_validate", shard: str = "s0000"):
    """Rerun one shard of sweep_rp_1792381414 in its original order; every call must match."""
    from deliverome_analysis.core.s3 import processed

    pre = "sprint4/s4e9-orf-design-fg-panel/reprediction/sweep_rp_1792381414/"
    sub = {
        r["sequence_sha256"]: r["sequence"]
        for r in json.loads(Path(processed.download(pre + "submission.json")).read_text())
    }
    pinned = {
        r["sequence_sha256"]: r
        for r in json.loads(Path(processed.download(pre + "reprediction.json")).read_text())
    }
    order = json.loads(
        Path(__file__).with_name(f"sweep_rp_1792381414_{shard}_ids.json").read_text()
    )
    calls = _run([(s, sub[s]) for s in order], f"validate_sweep_rp_1792381414_{shard}", Path(out))
    diff = [
        s
        for s in order
        if (calls[s]["prediction"], calls[s]["cleavage_site"])
        != (pinned[s]["prediction"], pinned[s]["cleavage_site"])
    ]
    dp = max(abs(calls[s]["sp_probability"] - pinned[s]["sp_probability"]) for s in order)
    dc = max(
        abs((calls[s]["cleavage_probability"] or 0) - (pinned[s]["cleavage_probability"] or 0))
        for s in order
    )
    report = {
        "shard": shard,
        "n": len(order),
        "calls_differ": diff,
        "max_abs_sp_probability_difference": dp,
        "max_abs_cleavage_probability_difference": dc,
    }
    (Path(out) / "validation.json").write_text(json.dumps(report, indent=1))
    print(json.dumps(report, indent=1))
