"""Rebuild the isoform universe, the isoform cohorts and the ortholog repicks as they were made.

The panels draw from two tables made on 2026-10-08 to 10-11: ``isoform_cohorts.tsv`` (every
alternate form worth considering, tiered by evidence) and ``ortholog_repick_{mouse,cyno}.tsv``
(each gene's mouse and cyno proteoform). Their scripts in ``library/l1/proteoform-selection``
have changed since, so this reruns **each step at the commit that produced it**, offline,
from the inputs pinned in ``reproduce_isoform_universe_inputs.json``:

| Step | Script | Commit | Checked against |
|---|---|---|---|
| universe: UniProt + Ensembl forms, PeptideAtlas tiers | isoform_universe.py | 88786350 | its committed summary, byte for byte |
| OpenProt 2.2 MS / translation evidence | openprot_evidence.py | 88786350 | committed summary, byte for byte |
| APPRIS TRIFID, membrane-loss filter | isoform_appris.py | a9646c44 | committed summary, byte for byte |
| topology shortlist | shortlist_for_topology.py | 39e19e2c | committed FASTA (see below) |
| DeepTMHMM v1 verdicts, from the committed 3-line output | apply_topology.py | e944c759 | committed summary, byte for byte |
| NetGPI form calls, then the cohorts | netgpi_gpi_forms.py, isoform_cohorts.py | 6fb98222 | the record (see below) |
| mouse and cyno ortholog repick | repick_orthologs.py | c4dbe752 | the records, byte for byte |

**What cannot be reproduced byte for byte, and why.** ``translations_by_gene`` returns Python
sets, whose iteration order follows the interpreter's per-process string-hash seed, and that
seed was not recorded. Two outputs inherit it:

* the order of a gene's Ensembl-only forms. ``isoform_cohorts.tsv`` is written in that order,
  so a few hundred of its 1,894 rows (all Ensembl-only; 426 to 428 in test runs) sit in a
  different position, and how many varies with each run's seed.
  Every row, the gene order and the summary reproduce exactly; the check requires that.
* which translation the shortlist writes for an Ensembl-only form whose gene has another
  translation of the same length (keyed ``(gene, "", length)``, first one wins). Against the
  1,527 committed sequences, a run differs on 78 to 95 such same-length siblings. DeepTMHMM ran on the committed
  sequences, so the committed FASTA and 3-line file are used as the record, and the check
  requires every difference to be exactly this case.

Two inputs are reconstructions rather than copies; each is the only variant that reproduces
the records:

* PeptideAtlas returns one pooled peptide table per 20-accession request, so a form's tier
  depends on which accessions shared its request. The pinned table replays the cache as it
  filled on 2026-10-08: the reference disagreements, the canonicals of genes with isoforms,
  all their proteoforms, then the rest of the universe (stage lists pinned alongside).
* ``surface_subcategory.tsv`` is the 2026-11-01 D1 architecture snapshot without EPGN's
  ``GPI_anchored`` row, which D1 gained after 10-11 (``gpi_exclusions.tsv`` now lists it as
  an annotation error). The surfaceome ortholog DeepTMHMM tables are deliberately absent: the
  repick ran without them and took its labels from the committed 3-line files.

    uv run --locked python reproduce_isoform_universe.py                # both checks
    uv run --locked python reproduce_isoform_universe.py --step orthologs
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tarfile
import tempfile
import zipfile
from io import BytesIO
from pathlib import Path

HERE = Path(__file__).parent
REPO = HERE.parents[2]
INPUTS = HERE / "reproduce_isoform_universe_inputs.json"
AP = "library/l1/proteoform-selection"
COMMITS = {
    "universe": "88786350",
    "openprot": "88786350",
    "appris": "a9646c44",
    "shortlist": "39e19e2c",
    "apply": "e944c759",
    "cohorts": "6fb98222",
    "orthologs": "c4dbe752",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def stage(work: Path) -> dict[str, Path]:
    """Download every pinned input into ``work/in`` and check its SHA-256."""
    from deliverome_analysis.core import s3

    spec = json.loads(INPUTS.read_text())
    got = {}
    for name, entry in spec["files"].items():
        bucket, _, key = entry["source"].partition(":")
        dest = work / "in" / name
        dest.parent.mkdir(parents=True, exist_ok=True)
        if not dest.exists():
            getattr(s3, bucket).download(key, dest)
        if sha256(dest) != entry["sha256"]:
            raise SystemExit(f"{name}: SHA-256 does not match {INPUTS.name}")
        got[name] = dest
    return got


def code_at(commit: str, dest: Path) -> Path:
    """The appris-principal folder (and the library list beside it) as of ``commit``."""
    if not dest.exists():
        dest.mkdir(parents=True)
        blob = subprocess.run(
            [
                "git",
                "-C",
                str(REPO),
                "archive",
                commit,
                AP,
                "library/l1/tag-site-architecture/library_members.json",
            ],
            check=True,
            capture_output=True,
        ).stdout
        tarfile.open(fileobj=BytesIO(blob)).extractall(dest, filter="data")
    return dest / AP


def run(code: Path, script: str, args: list[str], cwd: Path) -> None:
    env = {**os.environ, "PYTHONPATH": str(code)}
    proc = subprocess.run(
        [sys.executable, str(code / script), *args],
        cwd=cwd,
        env=env,
        capture_output=True,
        text=True,
    )
    if proc.returncode:
        raise SystemExit(f"{script} failed:\n{proc.stderr[-3000:]}")


def layout(work: Path, got: dict[str, Path]) -> tuple[Path, Path]:
    """The cache directory and accessible-surfaceome checkout the scripts expect."""
    spec = json.loads(INPUTS.read_text())["files"]
    cache, repo = work / "cache", work / "surfaceome"
    for name, entry in spec.items():
        where = entry.get("cache_as")
        if where:
            (cache / where).parent.mkdir(parents=True, exist_ok=True)
            shutil.copy(got[name], cache / where)
        where = entry.get("repo_path")
        if where:
            (repo / where).parent.mkdir(parents=True, exist_ok=True)
            shutil.copy(got[name], repo / where)
    op = work / "openprot"
    op.mkdir(exist_ok=True)
    with zipfile.ZipFile(got["openprot_min2pep.tsv.zip"]) as z:
        for member in z.namelist():
            if member.endswith(".tsv"):
                (op / Path(member).name).write_bytes(z.read(member))
    return cache, repo


def record(key: str, dest: Path) -> Path:
    from deliverome_analysis.core.s3 import processed

    return Path(processed.download(key, dest))


PANEL = "sprint4/s4e9-orf-design-fg-panel/library-design/panel-inputs/2026-11-01/"


def build_universe(work: Path) -> tuple[Path, dict[str, Path], Path, Path]:
    """Steps 1-3: universe, OpenProt evidence, APPRIS scoring. Returns the run folder."""
    got = stage(work)
    cache, repo = layout(work, got)
    out = work / "run"
    out.mkdir(exist_ok=True)
    lib88 = code_at(COMMITS["universe"], work / "code" / "u").parent / (
        "tag-site-architecture/library_members.json"
    )
    code = code_at(COMMITS["universe"], work / "code" / "u")
    # reference_selection.tsv is read only for its decision column, and only to send the
    # "appris_only" / "resolved_appris" genes to APPRIS; none were resolved_appris, and an
    # appris_only gene has no MANE, so the 2026-11-01 record decides identically.
    run(
        code,
        "isoform_universe.py",
        [
            "--surfaceome-repo",
            str(repo),
            "--cache-dir",
            str(cache),
            "--library",
            str(lib88),
            "--selection",
            str(got["reference_selection.tsv"]),
        ],
        out,
    )
    run(
        code,
        "openprot_evidence.py",
        [
            "--openprot-dir",
            str(work / "openprot"),
            "--cache-dir",
            str(cache),
        ],
        out,
    )
    code = code_at(COMMITS["appris"], work / "code" / "a")
    run(
        code,
        "isoform_appris.py",
        [
            "--cache-dir",
            str(cache),
            "--surfaceome-repo",
            str(repo),
        ],
        out,
    )
    return out, got, cache, repo


def read_fasta(path: Path) -> dict[str, str]:
    out, name = {}, None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].strip()
            out[name] = ""
        elif name:
            out[name] += line.strip()
    return out


def same_json(a: Path, b: Path) -> bool:
    return json.loads(a.read_text()) == json.loads(b.read_text())


def verify_isoform_universe(work: Path) -> str:
    out, got, cache, repo = build_universe(work)
    for step, name in (
        ("universe", "isoform_universe_summary.json"),
        ("openprot", "openprot_evidence_summary.json"),
        ("appris", "isoform_appris_summary.json"),
    ):
        code = code_at(COMMITS[step], work / "code" / step)
        assert same_json(out / name, code / name), f"{name} does not reproduce"

    code = code_at(COMMITS["shortlist"], work / "code" / "s")
    run(
        code,
        "shortlist_for_topology.py",
        [
            "--surfaceome-repo",
            str(repo),
            "--cache-dir",
            str(cache),
        ],
        out,
    )
    ours, rec = (
        read_fasta(out / "topology_shortlist.fasta"),
        read_fasta(code / "topology_shortlist.fasta"),
    )
    ours_by = {n.split("|", 1)[1]: s for n, s in ours.items()}
    rec_by = {n.split("|", 1)[1]: (n, s) for n, s in rec.items()}
    assert set(ours_by) == set(rec_by), "the shortlist names a different set of forms"
    sys.path.insert(0, str(code))
    from select_reference import translations_by_gene

    by_gene = translations_by_gene(cache / "Homo_sapiens.GRCh38.pep.all.fa.gz")
    # Forms whose recorded sequence has a same-length sibling (another translation of the same
    # Ensembl gene), i.e. the only forms a hash-ordered tie can change. Counted from the record
    # and the release in sorted order, so the number reported does not depend on this run's
    # hash seed.
    tied = sorted(
        tag
        for tag, (_, seq) in rec_by.items()
        if any(
            seq in t and any(len(o) == len(seq) and o != seq for o in t) for t in by_gene.values()
        )
    )
    for tag in sorted(rec_by):
        _, seq = rec_by[tag]
        if ours_by[tag] == seq:
            continue
        # Only an Ensembl-only form whose same-length sibling won the hash-ordered tie.
        assert (
            tag in tied
            and len(seq) == len(ours_by[tag])
            and any(seq in t and ours_by[tag] in t for t in by_gene.values())
        ), f"{tag}: unexplained shortlist difference"
    # The DeepTMHMM run is keyed by the committed ids, so the manifest follows them.
    with (out / "topology_shortlist.tsv").open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        cols, rows = reader.fieldnames, {r["fasta_id"].split("|", 1)[1]: r for r in reader}
    with (out / "topology_shortlist.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for name in rec:
            w.writerow({**rows[name.split("|", 1)[1]], "fasta_id": name})
    shutil.copy(code / "topology_shortlist.fasta", out / "topology_shortlist.fasta")

    code = code_at(COMMITS["apply"], work / "code" / "t")
    shutil.copy(code / "topology_shortlist.3line", out / "topology_shortlist.3line")
    run(code, "apply_topology.py", ["--surfaceome-repo", str(repo)], out)
    assert same_json(
        out / "topology_verdicts_summary.json", code / "topology_verdicts_summary.json"
    ), "verdicts do not reproduce"

    code = code_at(COMMITS["cohorts"], work / "code" / "c")
    for name in ("surface_subcategory.tsv", "gpi_supplement.tsv", "gpi_exclusions.tsv"):
        shutil.copy(got[name], code / name)
    run(
        code,
        "netgpi_gpi_forms.py",
        [
            "--human",
            str(code / "netgpi_gpi_human_predictions.txt"),
            "--orthologs",
            str(code / "netgpi_ortholog_predictions.txt"),
        ],
        out,
    )
    run(
        code,
        "isoform_cohorts.py",
        [
            "--surfaceome-repo",
            str(repo),
            "--cache-dir",
            str(cache),
        ],
        out,
    )
    assert same_json(out / "isoform_cohorts_summary.json", code / "isoform_cohorts_summary.json"), (
        "cohort summary does not reproduce"
    )
    rec_rows = record(PANEL + "isoform_cohorts.tsv", work / "rec.tsv").read_text().splitlines()
    new_rows = (out / "isoform_cohorts.tsv").read_text().splitlines()
    assert sorted(new_rows) == sorted(rec_rows), "isoform_cohorts.tsv rows do not reproduce"
    assert [r.split("\t")[0] for r in new_rows] == [r.split("\t")[0] for r in rec_rows], (
        "isoform_cohorts.tsv gene order does not reproduce"
    )
    moved = [n for n, r in zip(new_rows, rec_rows, strict=True) if n != r]
    assert all(r.split("\t")[5] == "ensembl_only" for r in moved), "a non-hash-ordered row moved"
    # The rows a hash-ordered tie can permute: Ensembl-only rows of genes with more than one,
    # counted from the record so the report is the same on every run.
    per_gene: dict[str, int] = {}
    for r in rec_rows[1:]:
        f = r.split("\t")
        if f[5] == "ensembl_only":
            per_gene[f[0]] = per_gene.get(f[0], 0) + 1
    permutable = sum(n for n in per_gene.values() if n > 1)
    return (
        f"isoform universe chain reruns at its producing commits: summaries byte for byte; "
        f"all {len(rec_rows) - 1} cohort rows and the gene order reproduce; only the "
        f"{permutable} Ensembl-only rows of genes with several may change order within their gene, "
        f"and only the {len(tied)} shortlist forms with a same-length sibling translation may "
        "change sequence (both hash-seed ties, checked row by row)"
    )


def verify_ortholog_repick(work: Path) -> str:
    out, got, cache, repo = build_universe(work)
    code = code_at(COMMITS["orthologs"], work / "code" / "o")
    for name in ("surface_subcategory.tsv", "gpi_supplement.tsv"):
        shutil.copy(got[name], code / name)
    for name in ("ortholog_repick_changed.3line", "ortholog_missing_topology.3line"):
        shutil.copy(code / name, out / name)
    shutil.copy(code / "netgpi_ortho_candidate_predictions.txt", out)
    for sp in ("mouse", "cyno"):
        run(
            code,
            "repick_orthologs.py",
            [
                "--species",
                sp,
                "--surfaceome-repo",
                str(repo),
                "--cache-dir",
                str(cache),
            ],
            out,
        )
        rec = record(PANEL + f"ortholog_repick_{sp}.tsv", work / f"rec_{sp}.tsv")
        assert sha256(out / f"ortholog_repick_{sp}.tsv") == sha256(rec), (
            f"ortholog_repick_{sp}.tsv does not reproduce its record"
        )
    return "mouse and cyno ortholog repicks rerun at c4dbe752 and reproduce byte for byte"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--step", choices=["isoforms", "orthologs"])
    ap.add_argument("--work", type=Path, default=None)
    args = ap.parse_args(argv)
    work = args.work or Path(tempfile.mkdtemp(prefix="isoform-universe-"))
    if args.step in (None, "isoforms"):
        print(verify_isoform_universe(work / "iso"))
    if args.step in (None, "orthologs"):
        print(verify_ortholog_repick(work / "orth"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
