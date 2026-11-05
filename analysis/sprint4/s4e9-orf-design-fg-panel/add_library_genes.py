"""Add genes to the L1 library after its gene list was frozen, without moving anything else.

The 2,368-gene root is frozen (``library_pipeline.py``, stage ``gene_list``). A gene admitted
later is listed in ``gene_list_additions.tsv`` with the decision, the date and the surfaceome
release that admitted it, and enters as a library reference through the same stages every
other reference went through, run from the same pinned inputs:

1. **reference** -- ``select_reference.py --offline`` over root + additions. Every recorded
   row must come back byte for byte; only the added genes' rows are new.
2. **candidate_cds** -- ``resolve_candidate_cds.py`` against that table. Every recorded CDS
   must come back unchanged; the added gene's ``lib:<gene>`` entry is taken.
3. **constructs** -- domesticated with the same enzymes and seed
   (``domesticate_candidates.domesticate_one``) and dropped into the pinned p0039 backbone
   (``assemble_candidates.assemble_one``). The vector must show the expected site pattern.
4. **tranches** -- appended to ``candidate_manifest_tranches.tsv`` as the next free L1 id.
   Nothing is renumbered. Its designed protein may not repeat an active member's.
5. **topology** -- its designed protein must already have a DeepTMHMM v1 row in the pinned
   snapshot; if not, the script stops and names ``publish_topology_supplement.py``.

An added gene is not eligible for the isoform or ortholog panels: those draw from the pinned
isoform universe and ortholog repick (``panel_inputs.json``), which predate it.

**Rebuilt genes** are root genes that had no member because GRCh38 carries a
loss-of-function copy (Ensembl ``protein_coding_LoF``) and the ladder's pick was that
transcript. A named override with basis ``reference_allele_lof`` in
``reference_overrides.tsv`` points the reference at the intact UniProt protein; the gene then
goes through the same stages as an addition. Only its own reference row may change, and every
other recorded row and CDS must still reproduce.

The pre-tranche ``candidate_manifest.tsv`` and ``candidate_vectors.fasta.gz`` on S3 number
members by sorted position, a namespace the tranche split abandoned, so they are left as the
historical intermediate. Added genes' vectors go to ``library-candidates/added_vectors.fasta``
under their final L1 id.

    uv run --locked python add_library_genes.py              # dry run: reproduce and report
    uv run --locked python add_library_genes.py --write      # append to the committed records
    uv run --locked python add_library_genes.py --publish    # after committing: upload to S3
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import library_pipeline as lp  # noqa: E402

sys.path.insert(0, str(lp.AP))
import assemble_candidates  # noqa: E402
import candidate_inputs  # noqa: E402
import domesticate_candidates  # noqa: E402
from Bio import SeqIO  # noqa: E402

from deliverome_analysis import library_design as ld  # noqa: E402
from deliverome_analysis import seqtools, topology  # noqa: E402

MANIFEST = HERE / "candidate_manifest_tranches.tsv"
COVERAGE = lp.AP / "candidate_cds_coverage.tsv"
SUMMARY = lp.AP / "reference_selection_summary.json"
DESIGN = "library-design/"
CANDIDATES = "library-candidates/"
VECTORS = CANDIDATES + "added_vectors.fasta"
OVERRIDES = lp.AP / "reference_overrides.tsv"
REBUILD_BASIS = "reference_allele_lof"


def rebuilt_genes() -> list[str]:
    """Root genes whose reference a ``reference_allele_lof`` override repoints."""
    with OVERRIDES.open() as fh:
        return [
            r["gene"] for r in csv.DictReader(fh, delimiter="\t") if r["basis"] == REBUILD_BASIS
        ]


def unix(path: Path) -> list[str]:
    return path.read_bytes().replace(b"\r\n", b"\n").decode().splitlines()


def reproduce(work: Path) -> dict:
    """Rerun every stage for root + additions and check nothing recorded moved."""
    with lp.ADDITIONS.open() as fh:
        additions = [r["gene"] for r in csv.DictReader(fh, delimiter="\t")]
    rebuilt = rebuilt_genes()
    additions += [g for g in rebuilt if g not in additions]
    with MANIFEST.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        cols, rows = list(reader.fieldnames or []), list(reader)
    pending = [g for g in additions if not any(r["key"] == f"lib:{g}" for r in rows)]

    # 1. reference
    out = lp.run_reference(work / "ref")
    rec = unix(lp.record(DESIGN + "reference_selection.tsv", work / "rec_ref.tsv"))
    new = unix(out / "reference_selection.tsv")
    assert new[0] == rec[0], "reference_selection.tsv header changed"
    # A rebuilt gene's row is replaced by its override; no other recorded row may move.
    gone = sorted(line for line in set(rec) - set(new) if line.split("\t")[0] not in rebuilt)
    assert not gone, f"{len(gone)} recorded reference row(s) changed, e.g. {gone[0][:80]}"
    added_rows = {line.split("\t")[0] for line in set(new) - set(rec)}
    assert added_rows <= set(additions), f"reference rows for unlisted genes: {added_rows}"
    assert set(pending) <= added_rows, f"pending {pending} not new to the reference record"
    dis = lp.record(DESIGN + "reference_disagreements.tsv", work / "rec_dis.tsv")
    assert lp.same_bytes(out / "reference_disagreements.tsv", dis), "disagreements changed"

    # 2. designed CDS
    cds = lp.run_candidate_cds(work / "cds", out / "reference_selection.tsv")
    rec_cds = json.loads(lp.record(CANDIDATES + "candidate_cds.json", work / "rc.json").read_text())
    moved = sorted(k for k in rec_cds if k not in cds or cds[k]["cds"] != rec_cds[k]["cds"])
    assert not moved, f"recorded CDS changed for {moved[:5]}"
    keys = [f"lib:{g}" for g in additions]
    lost = [k for k in keys if k not in cds]
    if lost:
        raise SystemExit(f"no transcript encodes the selected reference for {lost}")
    coverage = {
        line.split("\t")[0]: line
        for line in unix(work / "cds" / "candidate_cds_coverage.tsv")
        if line.split("\t")[0] in keys
    }

    # 3. constructs
    rec_con = json.loads(
        lp.record(CANDIDATES + "candidate_constructs.json", work / "rk.json").read_text()
    )
    bb = str(SeqIO.read(candidate_inputs.backbone(), "genbank").seq).upper()
    bb_sites = seqtools.find_sites(bb, assemble_candidates.ENZYMES)
    active = {
        seqtools.translate(rec_cds[r["key"]]["cds"]).rstrip("*"): r["candidate_id"]
        for r in rows
        if r["library"] == ld.LIBRARY_L1 and ld.is_active(r) and r["key"] in rec_cds
    }
    next_n = max(int(r["candidate_id"].split("-p")[1]) for r in rows if r["library"] == "l1") + 1

    ids = {r["key"]: r["candidate_id"] for r in rows if r["library"] == ld.LIBRARY_L1}
    added, constructs, fasta = [], {}, []
    for key in keys:
        protein = seqtools.translate(cds[key]["cds"]).rstrip("*")
        if key not in ids and protein in active:
            raise SystemExit(f"{key}: designed protein duplicates {active[protein]}")
        hit = topology.for_sequence(protein)
        if hit is None:
            raise SystemExit(
                f"{key}: no DeepTMHMM v1 row for its designed protein; run "
                "publish_topology_supplement.py fasta / modal sweep / publish first"
            )
        con = domesticate_candidates.domesticate_one(cds[key])
        if con is None:
            raise SystemExit(f"{key}: domestication would change the protein")
        cid = ids.get(key)
        if cid is None:
            cid = f"l1-p{next_n:04d}"
            next_n += 1
        vector, row, clean = assemble_candidates.assemble_one(bb, bb_sites, cid, key, con)
        if not clean:
            raise SystemExit(f"{key}: vector does not show the expected restriction-site pattern")
        constructs[key] = con
        fasta.append(assemble_candidates.fasta_record(cid, key, vector))
        entry = {
            **{c: "" for c in cols},
            **{c: str(v) for c, v in row.items() if c in cols},
            "library": ld.LIBRARY_L1,
            "prior_candidate_id": cid,
        }
        if key in ids:
            recorded = next(r for r in rows if r["candidate_id"] == cid)
            assert recorded == entry, f"{cid}: committed manifest row does not reproduce"
        else:
            added.append(entry)
        print(
            f"  {cid} {key}: {cds[key]['transcript']} ({cds[key]['route']}), "
            f"{len(protein)} aa, {con['edits']} domestication edit(s), {len(vector)} bp vector; "
            f"v1 {hit['deeptmhmm_label']} tm={hit['tm_helix_count']} "
            f"({hit['topology_version']}/{hit['cohort']})"
        )
    return {
        "pending": pending,
        "cols": cols,
        "rows": rows,
        "added": added,
        "reference_dir": out,
        "cds": {**rec_cds, **{k: cds[k] for k in keys}},
        "constructs": {**rec_con, **constructs},
        "coverage": coverage,
        "fasta": fasta,
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--write", action="store_true", help="append to the committed records")
    ap.add_argument("--publish", action="store_true", help="upload the S3 records")
    args = ap.parse_args(argv)

    work = Path(tempfile.mkdtemp(prefix="add-genes-"))
    r = reproduce(work)
    print(f"\nreference and CDS records reproduce; pending additions: {r['pending'] or 'none'}")

    if args.write and r["added"]:
        with MANIFEST.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=r["cols"], delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(r["rows"] + r["added"])
        # Replace a rebuilt gene's "unresolved" line rather than adding a second one.
        lines = [ln for ln in unix(COVERAGE) if ln.split("\t")[0] not in r["coverage"]]
        COVERAGE.write_text("\n".join(lines + list(r["coverage"].values())) + "\n")
        SUMMARY.write_text((r["reference_dir"] / "reference_selection_summary.json").read_text())
        print(f"appended {len(r['added'])} row(s) to {MANIFEST.name}; updated {COVERAGE.name}")

    if args.publish:
        if r["added"]:
            raise SystemExit("run --write and commit before --publish")
        from deliverome_analysis.s3 import processed

        for name in ("reference_selection.tsv", "reference_selection_summary.json"):
            processed.upload(str(r["reference_dir"] / name), lp.PREFIX + DESIGN + name)
        for name, data in (
            ("candidate_cds.json", r["cds"]),
            ("candidate_constructs.json", r["constructs"]),
        ):
            path = work / name
            path.write_text(json.dumps(data, indent=0, sort_keys=True))
            processed.upload(str(path), lp.PREFIX + CANDIDATES + name)
        path = work / "added_vectors.fasta"
        path.write_text("".join(r["fasta"]))
        processed.upload(str(path), lp.PREFIX + VECTORS)
        processed.upload(str(MANIFEST), lp.PREFIX + CANDIDATES + MANIFEST.name)
        print(f"published reference, CDS, constructs, {len(r['fasta'])} added vector(s), manifest")
    elif not args.write:
        print("dry run -- --write appends locally; commit, then --publish uploads")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
