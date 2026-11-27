"""Archive every AlphaFold DB model an s4e9 or s5e13 analysis reads, write-once.

Several analyses read AlphaFold DB models by URL only. AFDB can replace a file behind the same
URL (it is on its sixth version) and, from this machine, began refusing file downloads on
2026-11-13, so a URL is not a record. This copies each model to the external bucket under its
AFDB version:

    external://alphafold-db/v<N>/AF-<accession>-F1-model_v<N>.pdb

(``tag_refold_afdb_sources.tsv`` names bytes another session downloaded, with their SHA-256)
and writes ``tag_refold_afdb_archive.tsv`` (one row per model: AFDB URL, version, archive key,
SHA-256, bytes, where the bytes came from, and which readers use it). Readers:

| tag | models | read by |
|---|---|---|
| ``refold`` | native model of every live tagged-refold verdict | ``tag_refold/score_tagged_folds.py`` (reads only the archive) |
| ``native_index`` | every AFDB model in ``tag_refold_native_models.tsv`` (natives with or without a live verdict) | ``tag_refold/`` scripts; l1-reviewer-package's manifest |
| ``payload`` | every ``model_url`` in ``library_structure_payload.tsv`` (L1 and L2 members, orthologs included) | ``review_l1_cterm_members.py``, ``nterm_res1_lane_check.py``, ``membrane-side-gate``, s5e13 ``ortholog_member_features.py`` / ``designed_orf_feature_plan.py`` |
| ``gate`` | the AFDB files pinned in ``membrane-side-gate/input_manifest.tsv`` | ``membrane-side-gate`` (its SHA-256 must match) |
| ``s5e13`` | the v6 model of every accession in s5e13 ``site_features.tsv`` | s5e13 ``fetch_corpus_models.py``, ``panel_loop_arms.py``, ``tag_site_model_v5_geometry.py`` |
| ``request`` | models other sessions asked for (``tag_refold_afdb_requests.tsv``) | as listed there |

Bytes come from the local caches the readers wrote (``~/.cache/deliverome-analysis/**``), which
must agree with each other. A model no cache holds is not fetched (AFDB refuses this machine,
and routing around that is not this script's call); it is listed in
``tag_refold_afdb_unarchived.tsv`` with its readers. Refold natives use the exact copies they were
scored against and must contain their native protein; gate models must match the gate's recorded
SHA-256.

Keys are new and never overwritten (the external bucket is append-only by policy): a key that
already exists is compared with its recorded SHA-256 and left alone.

Requested by the coordinator after an independent reproducibility audit (2026-11-14), widened
the same day to every s4e9/s5e13 reader; implemented by Claude Code session l1-tag-refold.

    uv run --locked python library/l1/orf-design/tag_refold/archive_afdb_models.py
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import re
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
EXP = HERE.parent
REPO = EXP.parents[2]
S5E13 = REPO / "analysis" / "l1" / "tag-site-benchmarking"
sys.path.insert(0, str(HERE))

import structure_metrics as sm  # noqa: E402

CACHES = Path.home() / ".cache" / "deliverome-analysis"
NATIVE_CACHE = CACHES / "native-models"
OUT = EXP / "tag_refold_afdb_archive.tsv"
UNARCHIVED = EXP / "tag_refold_afdb_unarchived.tsv"
PREFIX = "alphafold-db"
FILES = "https://alphafold.ebi.ac.uk/files/"
NAME = re.compile(r"(AF-[^/_]+-F1-model_v(\d+)\.pdb)$")
SCORED = "alphafold.ebi.ac.uk, 2026-11-12/13 (the copies the refold verdicts were scored against)"
LOCAL = "alphafold.ebi.ac.uk, downloaded earlier by the reader (local cache)"


def tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def refold_natives() -> dict[str, str]:
    """AFDB URL -> native sequence, for every native a live refold verdict compares against."""
    index = {r["native_sha256"]: r for r in tsv(EXP / "tag_refold_native_models.tsv")}
    natives = {r["native_sha256"]: r["native_seq"] for r in tsv(EXP / "tag_refold_natives.tsv")}
    out: dict[str, str] = {}
    for r in tsv(EXP / "tag_refold_metrics.tsv"):
        if r.get("verdict") in ("superseded", "not_folded", "", None):
            continue
        p = index.get(r["native_sha256"])
        if p and p["ref"].startswith("https://"):
            out.setdefault(p["ref"], natives[r["native_sha256"]])
    return out


def wanted() -> tuple[dict[str, set[str]], dict[str, str], dict[str, str]]:
    """URL -> reader tags; refold URL -> native sequence; gate URL -> pinned SHA-256."""
    readers: dict[str, set[str]] = defaultdict(set)
    natives = refold_natives()
    for url in natives:
        readers[url].add("refold")
    for r in tsv(EXP / "tag_refold_native_models.tsv"):
        if r["ref"].startswith("https://"):
            readers[r["ref"]].add("native_index")
    for r in tsv(EXP / "library_structure_payload.tsv"):
        if r["model_url"] and not NAME.search(r["model_url"]):
            # no accession in the URL (TPO l1-m2323 on 2026-11-14): nothing to archive; reported
            print(f"  payload {r['candidate_id']}: malformed model_url {r['model_url']}, skipped")
        elif r["model_url"]:
            readers[r["model_url"]].add("payload")
    gate: dict[str, str] = {}
    for r in tsv(EXP / "membrane-side-gate" / "input_manifest.tsv"):
        m = NAME.search(r["file"])
        if m:
            readers[FILES + m.group(1)].add("gate")
            gate[FILES + m.group(1)] = r["sha256"]
    for acc in {r["uniprot_acc"] for r in tsv(S5E13 / "site_features.tsv")} - {""}:
        readers[f"{FILES}AF-{acc}-F1-model_v6.pdb"].add("s5e13")
    for r in tsv(EXP / "tag_refold_afdb_requests.tsv"):
        readers[r["model_url"]].add("request")
    for url in readers:
        if not url.startswith(FILES) or not NAME.search(url):
            raise SystemExit(f"not an AFDB file URL: {url}")
    return readers, natives, gate


def local_copies() -> dict[str, set[Path]]:
    """AFDB file name -> local copies, from every reader cache (``<acc>__`` prefixes stripped)."""
    out: dict[str, set[Path]] = defaultdict(set)
    for p in CACHES.rglob("*AF-*-F1-model_v*.pdb"):
        m = NAME.search(p.name.split("__")[-1])
        if m and p.is_file() and p.stat().st_size > 0:  # empty = a failed download
            out[m.group(1)].add(p)
    return out


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--missing", help="write the URLs no local cache holds to this file and stop")
    args = ap.parse_args()
    readers, natives, gate = wanted()
    named = {r["afdb_url"]: r for r in tsv(EXP / "tag_refold_afdb_sources.tsv")}
    copies = local_copies()
    if args.missing:
        miss = sorted(u for u in readers if not copies.get(u.rsplit("/", 1)[1]))
        Path(args.missing).write_text("".join(u + "\n" for u in miss))
        print(f"{len(miss)} of {len(readers)} models have no local copy -> {args.missing}")
        return

    from deliverome_analysis.core.s3 import external

    existing = set(external.list_keys(PREFIX + "/"))
    previous = {r["afdb_url"]: r for r in tsv(OUT)} if OUT.exists() else {}
    rows, uploaded, problems, unarchived = [], 0, [], []
    for url in sorted(readers):
        name, version = NAME.search(url).groups()
        key = f"{PREFIX}/v{version}/{name}"
        if url in natives:
            paths, source = {NATIVE_CACHE / name} & copies.get(name, set()), SCORED
        else:
            paths, source = copies.get(name, set()), LOCAL
            if url in named:  # bytes another session downloaded, named with their SHA-256
                source = named[url]["source"]
        if not paths:
            unarchived.append({"afdb_url": url, "afdb_version": f"v{version}",
                               "readers": ";".join(sorted(readers[url])),
                               "reason": "no local copy; AFDB refuses this machine"})  # fmt: skip
            continue
        digests = {sha(p) for p in paths}
        if url in natives and {sha(p) for p in copies[name]} != digests:
            print(
                f"  note {name}: other caches differ from the scored copy; the scored copy is archived"
            )
        elif len(digests) > 1:
            problems.append(f"{name}: local copies differ ({', '.join(str(p) for p in paths)})")
            continue
        path = sorted(paths)[0]
        digest = sha(path)
        if url in named and named[url]["sha256"] != digest:
            problems.append(f"{name}: differs from the SHA-256 its downloader recorded")
            continue
        if url in gate and gate[url] != digest:
            problems.append(f"{name}: differs from the SHA-256 the membrane-side gate pinned")
            continue
        if url in natives:
            model = sm.Model(path).seq
            if model.find(natives[url]) < 0 and model.find(natives[url][1:]) < 1:
                problems.append(f"{name}: does not contain its native protein")
                continue
        if not path.read_bytes().lstrip().startswith((b"HEADER", b"ATOM", b"REMARK", b"MODEL")):
            problems.append(f"{name}: not a PDB file")
            continue
        if key in existing:
            if previous.get(url, {}).get("sha256", digest) != digest:
                problems.append(f"{name}: archived copy differs from the local one")
                continue
            source = previous.get(url, {}).get("source", source)
        else:
            external.upload(path, key, content_type="chemical/x-pdb")
            uploaded += 1
        rows.append({
            "afdb_url": url, "afdb_version": f"v{version}", "archive_key": f"external://{key}",
            "sha256": digest, "bytes": path.stat().st_size, "source": source,
            "readers": ";".join(sorted(readers[url])),
        })  # fmt: skip
    if problems:
        print("\n".join(problems), file=sys.stderr)
        raise SystemExit(f"{len(problems)} models not archived")
    lost = set(previous) - {r["afdb_url"] for r in rows}
    if lost:
        raise SystemExit(
            f"{len(lost)} previously archived models would leave the table, e.g. {sorted(lost)[0]}"
        )
    with OUT.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    with UNARCHIVED.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["afdb_url", "afdb_version", "readers", "reason"],
                           delimiter="\t", lineterminator="\n")  # fmt: skip
        w.writeheader()
        w.writerows(unarchived)
    if any("refold" in r["readers"] for r in unarchived):
        raise SystemExit("a refold native is not archived; the refold scorer would stop")
    print(f"{len(unarchived)} models not archived (no local copy) -> {UNARCHIVED.name}")
    print(
        f"{len(rows)} AFDB models archived under external://{PREFIX}/ ({uploaded} uploaded now); wrote {OUT.name}"
    )


if __name__ == "__main__":
    main()
