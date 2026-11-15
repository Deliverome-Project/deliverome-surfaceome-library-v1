"""Every AlphaFold DB model the L1 design read is in the one archive, with its SHA-256.

The design read AlphaFold DB (AFDB) model files by URL: the structure payload
(``library_structure_payload.tsv``, read by the C-terminal review, the lane check, the
membrane-side gate and the S5E13 feature scripts) and the refold natives
(``tag_refold_native_models.tsv``, read by ``tag_refold/score_tagged_folds.py``). A URL is not a
record, so every such model is archived, write-once, under ``external://alphafold-db/v6/`` by
l1-tag-refold's ``tag_refold/archive_afdb_models.py``, which records each model's archive key and
SHA-256 in ``tag_refold_afdb_archive.tsv``. That archive is the single source of truth for these
models (coordinator decision, 2026-11-14).

This script only checks it (pipeline stage ``structure_models``):

* every AFDB URL the two tables name has a row in ``tag_refold_afdb_archive.tsv`` with an
  ``external://`` key and a SHA-256. A model the archive's owner lists as not yet archived
  (``tag_refold_afdb_unarchived.tsv``) is reported as a counted exception; any other model the
  design read that is not archived fails;
* each archived object exists at its key with the recorded size (``--deep`` also downloads it
  and checks its SHA-256).

Requested by the coordinator (2026-11-14); implemented by the Claude Code session
l1-reviewer-package.

    uv run --locked python pin_structure_models.py check
    uv run --locked python pin_structure_models.py check --deep
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PAYLOAD = HERE / "library_structure_payload.tsv"
NATIVES = HERE / "tag_refold_native_models.tsv"
ARCHIVE = HERE / "tag_refold_afdb_archive.tsv"
UNARCHIVED = HERE / "tag_refold_afdb_unarchived.tsv"


def _tsv(path: Path) -> list[dict]:
    if not path.exists():
        return []
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def afdb_urls() -> list[str]:
    """Every AFDB model file the design read: payload models and refold natives."""
    urls = {
        r["model_url"] for r in _tsv(PAYLOAD) if r["source"] == "alphafold_db" and r["model_url"]
    }
    urls |= {r["ref"] for r in _tsv(NATIVES) if r["source"] == "alphafold_db"}
    # A payload row with no accession (TPO, "no_model") names no file; it is not a model read.
    return sorted(u for u in urls if "/AF--F1-" not in u)


def problems(deep: bool = False) -> tuple[list[str], int]:
    if not ARCHIVE.exists():
        return [f"{ARCHIVE.name} is missing (l1-tag-refold's archive_afdb_models.py writes it)"], 0
    archive = {r["afdb_url"]: r for r in _tsv(ARCHIVE)}
    out = []
    missing = [u for u in afdb_urls() if u not in archive]
    # Models the archive itself lists as not yet archived (AFDB refuses this machine) are a
    # named, counted exception; any other missing model fails.
    known = {r["afdb_url"] for r in _tsv(UNARCHIVED)}
    unexpected = [u for u in missing if u not in known]
    if unexpected:
        out.append(
            f"{len(unexpected)} AFDB models the design read are neither archived nor listed in "
            f"{UNARCHIVED.name}: {unexpected[:5]}"
        )
    bad = [u for u, r in archive.items() if not r["archive_key"].startswith("external://")
           or len(r["sha256"]) != 64]  # fmt: skip
    if bad:
        out.append(f"archive rows without an external:// key or a SHA-256: {bad[:5]}")
    if out:
        return out, len(archive)
    from deliverome_analysis.s3 import external

    for url in afdb_urls():
        if url not in archive:
            continue
        r = archive[url]
        key = r["archive_key"].removeprefix("external://")
        if deep:
            data = Path(external.download(key)).read_bytes()
            if hashlib.sha256(data).hexdigest() != r["sha256"]:
                out.append(f"{r['archive_key']}: SHA-256 differs")
        elif int(external._head(key)["ContentLength"]) != int(r["bytes"]):
            out.append(f"{r['archive_key']}: size differs from the archive record")
    return out, len(archive)


def check(deep: bool = False) -> str:
    found, n = problems(deep)
    if found:
        raise SystemExit("structure models: " + "; ".join(found))
    urls = afdb_urls()
    archive = {r["afdb_url"] for r in _tsv(ARCHIVE)}
    pending = [u for u in urls if u not in archive]
    return (
        f"{len(urls) - len(pending)} of the {len(urls)} AFDB models the design read are in the "
        f"external://alphafold-db/v6/ archive ({n} archived), each present with its recorded "
        + ("SHA-256" if deep else "size (SHA-256 with --deep)")
        + (
            f"; {len(pending)} not yet archived, each listed in {UNARCHIVED.name} (exception)"
            if pending
            else ""
        )
    )


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("command", choices=("check",))
    ap.add_argument("--deep", action="store_true")
    args = ap.parse_args(argv)
    print(check(args.deep))
    return 0


if __name__ == "__main__":
    sys.exit(main())
