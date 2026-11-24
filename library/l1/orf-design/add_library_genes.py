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
from member_id_allocations import (  # noqa: E402
    RENUMBER_HINT,
    next_l1_member_number,
    to_design_space,
)

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
RESTORED = HERE / "allele_restorations_cds.json"
REBUILD_BASIS = "reference_allele_lof"


def rebuilt_genes() -> list[str]:
    """Root genes whose reference a ``reference_allele_lof`` override repoints."""
    with OVERRIDES.open() as fh:
        return [
            r["gene"] for r in csv.DictReader(fh, delimiter="\t") if r["basis"] == REBUILD_BASIS
        ]


def unix(path: Path) -> list[str]:
    return path.read_bytes().replace(b"\r\n", b"\n").decode().splitlines()


def _rerun_reference(
    work: Path,
    additions: list[str],
    withdrawn: set,
    rebuilt,
    pending: list[str],
    restored_genes: list[str],
) -> Path:
    """Rerun the reference stage and check every recorded row comes back unchanged."""
    # 1. reference
    out = lp.run_reference(work / "ref")
    rec = unix(lp.record(lp.REFERENCE_RECORD + "reference_selection.tsv", work / "rec_ref.tsv"))
    new = unix(out / "reference_selection.tsv")
    assert new[0] == rec[0], "reference_selection.tsv header changed"
    # A rebuilt gene's row is replaced by its override; no other recorded row may move.
    # A withdrawn addition's row leaves with it (SCNN1B, 2026-11-22); its record stays dated.
    gone = sorted(
        line
        for line in set(rec) - set(new)
        if line.split("\t")[0] not in rebuilt and line.split("\t")[0] not in withdrawn
    )
    assert not gone, f"{len(gone)} recorded reference row(s) changed, e.g. {gone[0][:80]}"
    added_rows = {line.split("\t")[0] for line in set(new) - set(rec)}
    assert added_rows <= set(additions), f"reference rows for unlisted genes: {added_rows}"
    # Root genes being restored, and an addition already in the record but left unresolved
    # (SIGLEC16), have reference rows; every other pending gene must get a new one.
    recorded = {line.split("\t")[0] for line in rec}
    pending_new = [g for g in pending if g not in restored_genes and g not in recorded]
    assert set(pending_new) <= added_rows, f"pending {pending_new} not new to the reference record"
    dis = lp.record(lp.REFERENCE_RECORD + "reference_disagreements.tsv", work / "rec_dis.tsv")
    # An added gene whose MANE and APPRIS disagree adds its own disagreement row (CD79B and
    # SIRPB1, 2026-11-21); every recorded row must still come back unchanged.
    rec_dis, new_dis = unix(dis), unix(out / "reference_disagreements.tsv")
    assert new_dis[0] == rec_dis[0], "reference_disagreements.tsv header changed"
    lost_dis = sorted(set(rec_dis) - set(new_dis))
    assert not lost_dis, f"{len(lost_dis)} recorded disagreement row(s) changed: {lost_dis[0][:80]}"
    extra_dis = {line.split("\t")[0] for line in set(new_dis) - set(rec_dis)}
    assert extra_dis <= set(additions), f"disagreement rows for unlisted genes: {extra_dis}"
    return out


def _restore_and_cover(
    restored: dict, keys: list[str], cds: dict, work: Path
) -> tuple[list[str], list[str], dict]:
    """Add restored alleles to ``cds``; the keys with a CDS, the unresolved and coverage lines."""
    for k, entry in restored.items():
        if k in keys and k not in cds:
            cds[k] = {
                "cds": entry["cds"],
                "origin": "library",
                "route": entry["route"],
                "species": "human",
                "transcript": entry["transcript"],
                "variant_residues": len(entry["residues_differing"]),
                "restoration": {e: entry[e] for e in ("edits", "expected", "residues_differing")},
            }
    unresolved = [k for k in keys if k not in cds]
    if unresolved:
        print(f"unresolved (no member): {unresolved}")
    all_keys, keys = keys, [k for k in keys if k in cds]
    coverage = {
        line.split("\t")[0]: line
        for line in unix(work / "cds" / "candidate_cds_coverage.tsv")
        if line.split("\t")[0] in all_keys  # an unresolved addition keeps its "unresolved" line
    }
    for k, entry in restored.items():  # a restored gene's coverage line says how it was built
        if k in cds and cds[k].get("restoration"):
            coverage[k] = "\t".join(
                [
                    k,
                    "library",
                    entry["transcript"],
                    "ok",
                    entry["route"],
                    f"{len(entry['cds'])} nt; allele_restorations.tsv",
                ]
            )
    return keys, unresolved, coverage


def reproduce(work: Path) -> dict:
    """Rerun every stage for root + additions and check nothing recorded moved."""
    additions = lp.additions()
    withdrawn = set(lp.additions(withdrawn=True))
    rebuilt = rebuilt_genes()
    additions += [g for g in rebuilt if g not in additions]
    # Root genes whose reference genome copy is broken but a working copy is recorded
    # (allele_restorations.tsv: a screened patch form or a common intact allele, 2026-11-22).
    restored = json.loads(RESTORED.read_text()) if RESTORED.exists() else {}
    restored_genes = [k.split(":", 1)[1] for k in restored]
    additions += [g for g in restored_genes if g not in additions]
    with MANIFEST.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        cols, rows = list(reader.fieldnames or []), to_design_space(list(reader))
    pending = [g for g in additions if not any(r["key"] == f"lib:{g}" for r in rows)]

    # 1. reference
    out = _rerun_reference(work, additions, withdrawn, rebuilt, pending, restored_genes)

    # 2. designed CDS
    cds = lp.run_candidate_cds(work / "cds", out / "reference_selection.tsv")
    rec_cds = json.loads(lp.record(CANDIDATES + "candidate_cds.json", work / "rc.json").read_text())
    # The resolver's record predates the named references of 2026-11-12 (apply_reference_fixes),
    # so exactly those come back different, as in library_pipeline.verify_candidate_cds.
    import apply_reference_fixes as fixes

    named = {f"lib:{g}" for g in fixes.NAMED}
    moved = sorted(
        k
        for k in rec_cds
        if (k not in cds or cds[k]["cds"] != rec_cds[k]["cds"])
        and k not in named
        and not (k not in cds and k.split(":", 1)[1] in withdrawn)
    )
    assert not moved, f"recorded CDS changed for {moved[:5]}"
    keys = [f"lib:{g}" for g in additions]
    # A gene no transcript encodes gets no member, as for the root genes the resolver leaves
    # "unresolved" (FCGR2C, OR10AC1, ...): it stays on the gene list and is reported, never
    # built from a guessed sequence (SIGLEC16, 2026-11-21).
    keys, unresolved, coverage = _restore_and_cover(restored, keys, cds, work)

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
    next_n = next_l1_member_number(rows)

    ids = {r["key"]: r["candidate_id"] for r in rows if r["library"] == ld.LIBRARY_L1}
    added, constructs, fasta = [], {}, []
    no_topology: list[str] = []
    for key in keys:
        protein = seqtools.translate(cds[key]["cds"]).rstrip("*")
        if key not in ids and protein in active:
            raise SystemExit(f"{key}: designed protein duplicates {active[protein]}")
        hit = topology.for_sequence(protein)
        cid = ids.get(key)
        if cid is None:
            cid = f"l1-p{next_n:04d}"
            next_n += 1
        if hit is None:
            # Collected, not fatal one at a time, so one DeepTMHMM v1 run covers every addition.
            gene = key.split(":", 1)[1]
            sha = __import__("hashlib").sha256(protein.encode()).hexdigest()
            no_topology.append(f">{cid}|{gene}|human|{sha}\n{protein}\n")
            continue
        con = domesticate_candidates.domesticate_one(cds[key])
        if con is None:
            raise SystemExit(f"{key}: domestication would change the protein")
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
        if "design_id" in cols:  # the permanent design id, added to the manifest 2026-11-13
            entry["design_id"] = cid
        if key in ids:
            recorded = next(r for r in rows if r["candidate_id"] == cid)
            diff = {
                c: (recorded.get(c), entry.get(c)) for c in cols if recorded.get(c) != entry.get(c)
            }
            assert recorded == entry, f"{cid}: committed manifest row does not reproduce: {diff}"
        else:
            added.append(entry)
        print(
            f"  {cid} {key}: {cds[key]['transcript']} ({cds[key]['route']}), "
            f"{len(protein)} aa, {con['edits']} domestication edit(s), {len(vector)} bp vector; "
            f"v1 {hit['deeptmhmm_label']} tm={hit['tm_helix_count']} "
            f"({hit['topology_version']}/{hit['cohort']})"
        )
    if no_topology:
        out_fa = HERE / "topology_runs" / "additions_missing_v1.fasta"
        out_fa.write_text("".join(no_topology))
        raise SystemExit(
            f"{len(no_topology)} addition(s) have no DeepTMHMM v1 row for their designed protein; "
            f"wrote {out_fa.relative_to(HERE)}: run the modal sweep on it, then "
            "publish_topology_supplement.py publish"
        )
    # A withdrawn addition's member is retired, never deleted or renumbered.
    retire = {
        r["candidate_id"]: f"withdrawn_from_gene_list:{r['key'].split(':', 1)[1]} "
        "(Becca 2026-11-22; gene_list_additions.tsv)"
        for r in rows
        if r["library"] == ld.LIBRARY_L1
        and ld.is_active(r)
        and r["key"].startswith("lib:")
        and r["key"].split(":", 1)[1] in withdrawn
    }
    for r in rows:
        if r["candidate_id"] in retire:
            r[ld.RETIRED_COLUMN] = retire[r["candidate_id"]]
    return {
        "pending": pending,
        "retired": sorted(retire),
        "unresolved": unresolved,
        "cols": cols,
        "rows": rows,
        "added": added,
        "reference_dir": out,
        "cds": {**rec_cds, **{k: cds[k] for k in keys}},
        "new_cds": {k: cds[k] for k in keys},
        "new_constructs": constructs,
        "constructs": {**rec_con, **constructs},
        "coverage": coverage,
        "fasta": fasta,
    }


PEPTIDES = "peptideatlas_build607_pools.json"


def fetch_peptides(work: Path, snapshot: str) -> int:
    """Extend the pinned PeptideAtlas cache with the added genes' accessions (append-only)."""
    from deliverome_analysis.s3 import external

    out = lp.run_reference(work / "ref", offline=False)
    print(f"reference ladder ran online: {out}")
    old = json.loads(lp.stage_inputs([PEPTIDES], work / "pinned")[PEPTIDES].read_text())
    path = work / "ref" / "cache" / PEPTIDES
    new = json.loads(path.read_text())
    assert all(new["index"].get(a) == k for a, k in old["index"].items()), "an indexed entry moved"
    assert all(new["pools"].get(k) == v for k, v in old["pools"].items()), "a pooled table moved"
    added = sorted(set(new["index"]) - set(old["index"]))
    if not added:
        print("the pinned cache already holds every accession; nothing to publish")
        return 0
    key = f"peptideatlas/build-607/{snapshot}/{PEPTIDES}"
    if external.exists(key):
        raise SystemExit(f"external://{key} already exists; pick another snapshot name")
    external.upload(path, key)
    spec = json.loads(lp.INPUTS.read_text())
    entry = spec["files"][PEPTIDES]
    prior = entry["source"].split("/")[-2]
    entry["origin"] = (
        f"PeptideAtlas GetPeptides CGI, atlas_build_id 607 (Human 2026-01): the {prior} queries "
        f"(sha256 {entry['sha256'][:8]}…) unchanged, plus {len(added)} accession(s) queried "
        f"{snapshot.removeprefix('queries-')[:10]} by select_reference.fetch_peptides for the "
        f"genes added to the list ({', '.join(added)})"
    )
    entry["source"] = f"external:{key}"
    entry["sha256"] = lp.sha256(path)
    lp.INPUTS.write_text(json.dumps(spec, indent=1, ensure_ascii=False) + "\n")
    print(f"published external://{key} (+{len(added)} accessions) and repinned {lp.INPUTS.name}")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--write", action="store_true", help="append to the committed records")
    ap.add_argument("--publish", action="store_true", help="upload the S3 records")
    ap.add_argument(
        "--fetch-peptides",
        metavar="SNAPSHOT",
        help="query PeptideAtlas for accessions the pinned cache lacks, publish the extended "
        "cache as external peptideatlas/build-607/<SNAPSHOT>/ and repin it; then rerun",
    )
    ap.add_argument("--tag", default="2026-11-22-round2-additions", help="dated record key")
    args = ap.parse_args(argv)

    work = Path(tempfile.mkdtemp(prefix="add-genes-"))
    if args.fetch_peptides:
        return fetch_peptides(work, args.fetch_peptides)
    r = reproduce(work)
    print(f"\nreference and CDS records reproduce; pending additions: {r['pending'] or 'none'}")

    if r["retired"]:
        print(f"retired (withdrawn from the gene list): {r['retired']}")
    if args.write and (r["added"] or r["unresolved"] or r["retired"]):
        with MANIFEST.open("w", newline="") as fh:  # unchanged when nothing new was added
            w = csv.DictWriter(fh, fieldnames=r["cols"], delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(r["rows"] + r["added"])
        print(RENUMBER_HINT)
        # Replace a rebuilt gene's "unresolved" line rather than adding a second one.
        lines = [ln for ln in unix(COVERAGE) if ln.split("\t")[0] not in r["coverage"]]
        COVERAGE.write_text("\n".join(lines + list(r["coverage"].values())) + "\n")
        SUMMARY.write_text((r["reference_dir"] / "reference_selection_summary.json").read_text())
        (lp.AP / "reference_selection_sequences.tsv").write_text(
            (r["reference_dir"] / "reference_selection_sequences.tsv").read_text()
        )
        print(f"appended {len(r['added'])} row(s) to {MANIFEST.name}; updated {COVERAGE.name}")

    if args.publish:
        if r["added"] or r["retired"]:
            raise SystemExit("run --write and commit before --publish")
        # Records are versioned under dated keys (library_pipeline.REFERENCE_RECORD,
        # library_design.CDS_KEY / CONSTRUCTS_KEY): publish the next version, every earlier
        # entry unchanged plus the additions, then move the three constants to it.
        from deliverome_analysis.s3 import processed

        tag = args.tag
        ref_prefix = f"{lp.PREFIX}{DESIGN}reference-selection/{tag}/"
        cand_prefix = f"{lp.PREFIX}{CANDIDATES}{tag}/"
        for key in (ref_prefix, cand_prefix):
            if processed.list_keys(prefix=key):
                raise SystemExit(f"processed://{key} exists; choose a fresh --tag")
        for name in (
            "reference_selection.tsv",
            "reference_disagreements.tsv",
            "reference_selection_summary.json",
        ):
            processed.upload(str(r["reference_dir"] / name), ref_prefix + name)
        added_keys = [k for k in r["new_cds"] if k.startswith("lib:")]
        for name, current_key in (
            ("candidate_cds.json", ld.CDS_KEY),
            ("candidate_constructs.json", ld.CONSTRUCTS_KEY),
        ):
            data = json.loads(
                Path(processed.download(current_key, work / f"cur_{name}")).read_text()
            )
            source = r["new_cds"] if name == "candidate_cds.json" else r["new_constructs"]
            for k in added_keys:
                assert k not in data or data[k] == source[k], (
                    f"{k}: differs from the current record"
                )
                data[k] = source[k]
            path = work / name
            path.write_text(json.dumps(data, indent=0, sort_keys=True))
            processed.upload(str(path), cand_prefix + name)
        prev = processed.download(lp.PREFIX + VECTORS, work / "prev_vectors.fasta")
        path = work / "added_vectors.fasta"
        # Every addition's vector is rebuilt in r["fasta"]; the earlier file's must be among them.
        text = "".join(r["fasta"])
        old = Path(prev).read_text().split(">")[1:]
        assert all(">" + rec in text for rec in old), "an earlier added vector changed"
        path.write_text(text)
        processed.upload(str(path), cand_prefix + "added_vectors.fasta")
        processed.upload(str(MANIFEST), lp.PREFIX + CANDIDATES + MANIFEST.name)
        print(
            f"published {ref_prefix} and {cand_prefix} ({len(added_keys)} addition(s)); now set "
            f'library_pipeline.REFERENCE_RECORD = "{DESIGN}reference-selection/{tag}/" and '
            f"library_design.CDS_KEY / CONSTRUCTS_KEY to {cand_prefix}"
        )
    elif not args.write:
        print("dry run -- --write appends locally; commit, then --publish uploads")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
