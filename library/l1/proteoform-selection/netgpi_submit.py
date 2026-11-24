"""Submit a FASTA to the NetGPI 1.1 service and collect the tabular predictions.

NetGPI's trained weights are not distributed — the GitHub repository is MIT-licensed
training code whose `.gitignore` excludes the `*.pt` ensemble, and DTU publishes only the
dataset — so the service is the only way to run the published tool. Sequences submitted
here are public UniProt and Ensembl entries.

The service runs one job per site at a time and the long output format times out above 100
entries, so this submits `short` output in batches and waits for each before the next.

    uv run --locked python netgpi_submit.py netgpi_reference.fasta --out netgpi_reference_predictions.txt
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
import time
from pathlib import Path

CGI = "https://services.healthtech.dtu.dk/cgi-bin/webface2.cgi"
CONFIG = "/var/www/services/services/NetGPI-1.1/webface.cf"
RESULT = "https://services.healthtech.dtu.dk/services/NetGPI-1.1/tmp/{job}/output_protein_type.txt"
JOBID = re.compile(r"jobid=([A-F0-9]+)")


def records(path: Path) -> list[tuple[str, str]]:
    out, name, buf = [], None, []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name:
                out.append((name, "".join(buf)))
            name, buf = line[1:], []
        elif line.strip():
            buf.append(line.strip())
    if name:
        out.append((name, "".join(buf)))
    return out


def curl(args: list[str]) -> tuple[int, str, str]:
    p = subprocess.run(["curl", "-sS", *args], capture_output=True, text=True)
    return p.returncode, p.stdout, p.stderr


def submit(fasta: Path, retries: int = 40) -> str:
    """Post the batch; returns the job id. Waits out the one-job-per-site limit."""
    for attempt in range(1, retries + 1):
        rc, body, _ = curl(["-D", "-", "-o", "/dev/null", "-X", "POST", CGI,
                            "-F", f"configfile={CONFIG}",
                            "-F", f"uploadfile=@{fasta};type=text/plain",
                            "-F", "format=short"])
        m = JOBID.search(body)
        if m:
            return m.group(1)
        if attempt == 1:
            print("    queue busy, waiting for the running job to finish", flush=True)
        time.sleep(30)
    raise SystemExit(f"could not submit {fasta} after {retries} attempts")


def collect(job: str, expected: int, timeout: int = 1800) -> list[str]:
    url = RESULT.format(job=job)
    deadline = time.time() + timeout
    while time.time() < deadline:
        rc, body, _ = curl(["-w", "%{http_code}", url])
        if body.endswith("200"):
            text = body[:-3]
            rows = [ln for ln in text.splitlines() if ln and not ln.startswith("#")]
            if len(rows) >= expected:
                return rows
        time.sleep(15)
    raise SystemExit(f"job {job} did not produce {expected} rows in {timeout}s")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("fasta")
    ap.add_argument("--out", required=True)
    ap.add_argument("--batch", type=int, default=700)
    ap.add_argument("--work", default="", help="scratch dir for batch files")
    args = ap.parse_args()

    recs = records(Path(args.fasta))
    work = Path(args.work) if args.work else Path(args.fasta).parent
    batches = [recs[i:i + args.batch] for i in range(0, len(recs), args.batch)]
    print(f"{args.fasta}: {len(recs)} sequences in {len(batches)} batch(es)", flush=True)

    rows: list[str] = []
    for i, batch in enumerate(batches, 1):
        bf = work / f".netgpi_batch_{i}.fasta"
        bf.write_text("".join(f">{n}\n{s}\n" for n, s in batch))
        job = submit(bf)
        print(f"  batch {i}/{len(batches)}: {len(batch)} seqs, job {job}", flush=True)
        rows += collect(job, len(batch))
        print(f"  batch {i}/{len(batches)}: done ({len(rows)} rows so far)", flush=True)
        bf.unlink(missing_ok=True)

    Path(args.out).write_text(
        "# NetGPI 1.1 short output\n"
        "# ID\tSeq-length\tPred. GPI-Anchored\tOmega-site pos.\tLikelihood\tAmino-acid\n"
        + "\n".join(rows) + "\n"
    )
    print(f"wrote {args.out} ({len(rows)} predictions)")


if __name__ == "__main__":
    main()
