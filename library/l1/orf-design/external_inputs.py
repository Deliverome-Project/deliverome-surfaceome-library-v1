"""Every external source the L1 design reads, and the saved copy that replaces it.

``external_inputs.tsv`` has one row per script and external source (UniProt, Ensembl, MyGene,
AlphaFold DB, TmAlphaFold, NetGPI, D1, the licensed predictors and the Modal images). For each
row it records the version or release used, where the saved copy is, its SHA-256, and a status:

=================  ===========================================================================
``saved``          the script reads the saved copy by default; a live refresh is opt-in
``saved_not_read`` a saved copy exists, but rerunning the script queries the live source; the
                   committed output is what the build reads
``service_record`` a web service that cannot itself be saved; its responses are saved
``not_saveable``   a licensed tool or a container runtime; its version and file hashes are
                   recorded instead (``saved_copy`` names that record)
``output_only``    no raw response was saved; only the script's committed output exists
=================  ===========================================================================

``feeds_order`` says how the row reaches the L1 order: ``verify`` (a ``library_pipeline.py
verify`` stage runs the script), ``upstream`` (its committed output feeds a build) or ``no``.

``saved_copy`` is one of ``git:<path>`` (a committed file, SHA-256 in the ``sha256`` column),
``s3:<bucket>://<key>`` (SHA-256 in the column), or ``manifest:<path>`` (a committed manifest
that carries one SHA-256 per object: ``library_pipeline_inputs_saved.tsv``,
``structure_models_saved.tsv``, ``nterm/nterm_bundle_manifest.tsv``).

Checks:

* :func:`scan` (``tests/test_external_inputs.py``, offline): every s4e9 script that makes a
  network call is listed; a ``verify`` row is ``saved``; an ``upstream`` row has a saved copy
  or a recorded tool identity, or is an ``output_only`` row with a note saying what stands in
  for the response; every ``git:`` copy exists and matches its SHA-256 (a record the renumber
  stage re-keys is pinned by git and checked by that stage instead).
* :func:`check` (pipeline stage ``external_inputs``): the scan, plus every ``s3:`` copy and
  every object of every ``manifest:`` copy is on S3 with the recorded SHA-256.

Requested by the coordinator (2026-11-14, Becca's bar: "fully reproducible in 5 years or more
with all inputs pinned and saved so we aren't relying on external"); implemented by the Claude
Code session l1-reviewer-package.

    uv run --locked python external_inputs.py pin     # fill the sha256 column of git: rows
    uv run --locked python external_inputs.py scan    # offline
    uv run --locked python external_inputs.py check   # with S3
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

from deliverome_analysis.provenance import sha256_file, verify_s3_objects

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
REGISTRY = HERE / "external_inputs.tsv"
STATUSES = {"saved", "saved_not_read", "service_record", "not_saveable", "output_only"}
FEEDS = {"verify", "upstream", "no"}
#: A call that reaches the network: urllib, requests/httpx, curl/wget, a container image
#: pulled from a registry or built with pip/apt, and the shared helpers that query MyGene and
#: Ensembl. A URL in a docstring is not a call.
NETWORK = re.compile(
    r"urlopen\(|urlretrieve\(|requests\.(?:get|post|Session)\(|httpx\.|from_registry\(|"
    r"\.pip_install\(|\.apt_install\(|[\"']curl[\"']|[\"']wget[\"']|"
    r"gene_resolver\.(?:resolve|to_ensembl|to_canonical)\(|seqtools\.(?:resolve_gene|fetch_cds)\("
)
SKIP = ("tests/", "review-repo/")
#: Verify-stage scripts whose saved copy exists but whose reader on ``main`` still fetches live,
#: with the session switching it (2026-11-14). Remove each entry when its change lands; the
#: scan then requires the row to be ``saved``.
VERIFY_PENDING: dict[str, str] = {}  # emptied 2026-11-14: the C-terminal review reads saved copies


def rows() -> list[dict]:
    with REGISTRY.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def fetching_scripts() -> list[str]:
    out = []
    for path in sorted(HERE.rglob("*.py")):
        rel = path.relative_to(HERE).as_posix()
        if rel.startswith(SKIP) or "/review-repo/" in rel:
            continue
        if NETWORK.search(path.read_text()):
            out.append(rel)
    return out


def _git_path(copy: str) -> Path:
    return HERE / copy.removeprefix("git:")


def _manifest_objects(copy: str) -> list[tuple[str, str, str]]:
    """(bucket, key, sha256) for every S3 object a manifest lists."""
    path = HERE / copy.removeprefix("manifest:")
    with path.open() as fh:
        recs = list(csv.DictReader(fh, delimiter="\t"))
    out = []
    for r in recs:
        if "bucket" in r:
            out.append((r["bucket"], r["key"], r["sha256"]))
        elif r.get("ref", "").startswith("processed://"):
            out.append(("processed", r["ref"].removeprefix("processed://"), r["sha256"]))
        elif r.get("external_key", "").startswith("external://"):
            out.append(("external", r["external_key"].removeprefix("external://"), r["sha256"]))
        elif r.get("archive_key", "").startswith("external://"):
            out.append(("external", r["archive_key"].removeprefix("external://"), r["sha256"]))
    return out


def renumbered() -> set[str]:
    """Committed records the renumber stage re-keys (their bytes change with member ids)."""
    sys.path.insert(0, str(HERE))
    import renumber_library_members

    return set(renumber_library_members.RECORDS)


def scan() -> list[str]:
    """Offline problems with the registry (empty when it is complete)."""
    problems = []
    recs = rows()
    listed = {r["script"] for r in recs}
    for script in fetching_scripts():
        if script not in listed:
            problems.append(f"{script} makes a network call but has no row in {REGISTRY.name}")
    for r in recs:
        where = f"{r['script']} / {r['source']}"
        if r["status"] not in STATUSES or r["feeds_order"] not in FEEDS:
            problems.append(f"{where}: unknown status or feeds_order")
            continue
        if (r["feeds_order"] == "verify" and r["status"] != "saved"
                and r["script"] not in VERIFY_PENDING):  # fmt: skip
            problems.append(f"{where}: a verify stage must read a saved copy, not {r['status']}")
        if r["feeds_order"] == "upstream" and r["status"] == "output_only" and not r["note"]:
            problems.append(f"{where}: output_only upstream row needs a note")
        if r["status"] in {"saved", "saved_not_read", "service_record", "not_saveable"}:
            if r["feeds_order"] != "no" and not r["saved_copy"]:
                problems.append(f"{where}: {r['status']} needs a saved_copy")
        copy = r["saved_copy"]
        if copy.startswith("git:"):
            path = _git_path(copy)
            if not path.exists():
                problems.append(f"{where}: {path.relative_to(HERE)} is missing")
            elif copy.removeprefix("git:") in renumbered():
                continue  # re-keyed at every renumber; git pins it and the renumber stage checks it
            elif sha256_file(path) != r["sha256"]:
                problems.append(f"{where}: {path.relative_to(HERE)} SHA-256 changed")
        elif copy.startswith("manifest:"):
            if not (HERE / copy.removeprefix("manifest:")).exists():
                problems.append(f"{where}: manifest {copy} is missing")
        elif copy.startswith("s3:"):
            if not re.fullmatch(r"[0-9a-f]{64}", r["sha256"]):
                problems.append(f"{where}: s3 copy without a SHA-256")
        elif copy:
            problems.append(f"{where}: saved_copy {copy!r} is not git:, s3: or manifest:")
    return problems


def check() -> str:
    problems = scan()
    objects: dict[tuple[str, str], str] = {}
    for r in rows():
        copy = r["saved_copy"]
        if copy.startswith("s3:"):
            bucket, key = copy.removeprefix("s3:").split("://", 1)
            objects[(bucket, key)] = r["sha256"]
        elif copy.startswith("manifest:"):
            for bucket, key, sha in _manifest_objects(copy):
                objects[(bucket, key)] = sha

    items = [(bucket, key, sha) for (bucket, key), sha in sorted(objects.items())]
    missing, changed = verify_s3_objects(items, max_workers=12)
    problems.extend(missing)
    problems.extend(f"{label}: SHA-256 changed" for label in changed)
    if problems:
        raise SystemExit("external inputs:\n  " + "\n  ".join(problems))
    by = {}
    for r in rows():
        by[r["status"]] = by.get(r["status"], 0) + 1
    return (
        f"{len(rows())} external sources listed ({by}); {len(objects)} saved S3 objects "
        "unchanged; every network-calling script is listed"
    )


def pin() -> int:
    recs = rows()
    for r in recs:
        if r["saved_copy"].startswith("git:"):
            r["sha256"] = sha256_file(_git_path(r["saved_copy"]))
    with REGISTRY.open("w", newline="") as fh:
        w = csv.DictWriter(fh, list(recs[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(recs)
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("command", choices=("pin", "scan", "check"))
    args = ap.parse_args(argv)
    if args.command == "pin":
        return pin()
    if args.command == "scan":
        problems = scan()
        print("\n".join(problems) or "registry complete")
        return 1 if problems else 0
    print(check())
    return 0


if __name__ == "__main__":
    sys.exit(main())
