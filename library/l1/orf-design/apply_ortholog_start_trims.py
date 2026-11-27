"""Restart four cyno L1 orthologs at the native Met aligned to human Met1 (5' trim only).

Decided by the sessions l1-proteoform-audit (``l1_proteoform_decisions.tsv``, PR #300) and
ortholog-correspondence-review (``ortholog_restart_decisions.tsv``, PR 297), each under
Becca's 2026-11-12 delegation ("use your best judgment"), independently and with the same
restart hashes; applied by the Claude Code session l1-reviewer-package (subagent), PR #283.
``ortholog_start_trims.tsv`` is the ledger, one row per member with attribution, the rejected
alternative and the restart's expected CDS and protein SHA-256.

Each cyno Ensembl 116 model starts upstream of the residue aligned to human Met1, and Ensembl
116 has no other model for the gene. The restart drops the codons before that native in-frame
ATG and edits no codon. It is an explicit, recorded exception to ``ORTHOLOG_CDS_POLICY.md``'s
"no implicit start repairs" (column ``policy_exception``).

A restart is applied only when its protein has a DeepTMHMM v1 row in the pinned topology
snapshot (``status = applied``); the topology and tag-path stages need one, and no prediction
is run here. The other rows are ``held_pending_v1``: decided, applied once the v1 rows are
published. Every row is checked either way: the current CDS hashes to ``original_cds_sha256``,
the trimmed CDS and protein hash to the recorded values, the protein starts M, has no internal
stop and has the stated length.

The CDS and construct records extend the 2026-11-12 reference-fixes-2 records: only the
applied members' entries change, under a new dated key. Every other ortholog entry and every
other active L1 ortholog manifest row is asserted byte-identical. The adopted panel itself
(``adopt_ortholog_cds_policy.py``) is unchanged; its verify masks the applied members' design
columns, which this script verifies.

    uv run --locked python apply_ortholog_start_trims.py            # dry run
    uv run --locked python apply_ortholog_start_trims.py --write    # manifest + summary
    uv run --locked python apply_ortholog_start_trims.py --publish  # after committing: S3
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import library_pipeline as lp  # noqa: E402
from member_id_allocations import RENUMBER_HINT, to_design_space  # noqa: E402

sys.path.insert(0, str(lp.AP))
import assemble_candidates  # noqa: E402
import candidate_inputs  # noqa: E402
import domesticate_candidates  # noqa: E402
from Bio import SeqIO  # noqa: E402

from deliverome_analysis.cloning import library_design as ld  # noqa: E402
from deliverome_analysis.cloning import seqtools  # noqa: E402
from deliverome_analysis.surfaceome import topology  # noqa: E402

MANIFEST = HERE / "candidate_manifest_tranches.tsv"
LEDGER = HERE / "ortholog_start_trims.tsv"
SUMMARY = HERE / "ortholog_start_trims_summary.json"
PREV = "sprint4/s4e9-orf-design-fg-panel/library-candidates/2026-11-12-reference-fixes-2/"
# CCR5 alone was applied under 2026-11-12-ortholog-start-trims/ (kept); all four since the
# v1 rows for CHRM3, NIPAL4 and OR3A2 were published (topo_2026_11_12_cyno_restarts).
RECORDS = "sprint4/s4e9-orf-design-fg-panel/library-candidates/2026-11-12-ortholog-start-trims-2/"
DESIGN_COLS = (
    "transcript",
    "route",
    "variant_residues",
    "orf_nt",
    "vector_bp",
    "sites_removed",
    "edits",
    "vector_sites",
    "esp3i_sites",
    "insert_nt",
    "has_igk_leader",
)


def sha(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


def ledger() -> list[dict]:
    with LEDGER.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def applied_keys() -> set[str]:
    return {r["key"] for r in ledger() if r["status"] == "applied"}


def restarted(entry: dict, row: dict) -> dict:
    """The entry restarted at ``restart_codon``, after every recorded check."""
    key = row["key"]
    if sha(entry["cds"]) != row["original_cds_sha256"]:
        raise SystemExit(f"{key}: CDS differs from the reviewed original")
    if entry["transcript"] != row["transcript"]:
        raise SystemExit(f"{key}: transcript is {entry['transcript']}")
    cds = entry["cds"][3 * (int(row["restart_codon"]) - 1) :]
    prot = seqtools.translate(cds)
    checks = {
        "starts ATG": cds.startswith("ATG"),
        "protein starts M": prot.startswith("M"),
        "ends in a stop": prot.endswith("*"),
        "no internal stop": "*" not in prot[:-1],
        f"{row['protein_aa']} aa": len(prot) - 1 == int(row["protein_aa"]),
        "CDS hash": sha(cds) == row["restarted_cds_sha256"],
        "protein hash": sha(prot.rstrip("*")) == row["restarted_protein_sha256"],
    }
    bad = [k for k, ok in checks.items() if not ok]
    if bad:
        raise SystemExit(f"{key}: restart fails {bad}")
    return {**entry, "cds": cds, "trimmed_upstream_codons": int(row["restart_codon"]) - 1}


def reproduce(work: Path) -> dict:
    from deliverome_analysis.core.s3 import processed

    prev_cds = json.loads(
        Path(processed.download(PREV + "candidate_cds.json", work / "c")).read_text()
    )
    prev_con = json.loads(
        Path(processed.download(PREV + "candidate_constructs.json", work / "k")).read_text()
    )
    with MANIFEST.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        cols, rows = list(reader.fieldnames or []), to_design_space(list(reader))
    excl = lp.exclusions()

    def _adopted_member(r: dict) -> bool:
        if r["library"] != ld.LIBRARY_L1 or r["set"] != "ortholog":
            return False
        if ld.is_active(r):
            return True
        ret = (r.get(ld.RETIRED_COLUMN) or "").strip()
        if ret.startswith("parent_gene_dropped_ectodomain_below_30aa:"):
            gene = ret.split(":", 1)[1].split(" ", 1)[0]
            assert gene in excl, f"{r['candidate_id']}: unexpected parent-gene retirement ({gene})"
            return True
        return False

    active = {r["key"]: r for r in rows if _adopted_member(r)}
    assert len(active) == 200, f"{len(active)} adopted L1 orthologs, not 200"
    before = {k: "\t".join(r[c] for c in cols) for k, r in active.items()}

    bb = str(SeqIO.read(candidate_inputs.backbone(), "genbank").seq).upper()
    bb_sites = seqtools.find_sites(bb, assemble_candidates.ENZYMES)
    changed, new_con, fasta, report = {}, {}, [], []
    for row in ledger():
        key = row["key"]
        assert (
            key in active
            and ld.is_active(active[key])
            and active[key]["candidate_id"] == row["design_id"]
        ), key
        e = restarted(prev_cds[key], row)  # every row is checked, applied or held
        hit = topology.for_sequence(seqtools.translate(e["cds"]).rstrip("*"), tool=topology.TOOL_V1)
        v1 = hit["topology_version"] if hit else ""
        report.append({"key": key, "status": row["status"], "v1_topology_version": v1})
        if row["status"] != "applied":
            continue
        if not hit:
            raise SystemExit(f"{key}: applied but its restarted protein has no v1 row")
        con = domesticate_candidates.domesticate_one(e)
        if con is None or seqtools.translate(con["construct"]) != seqtools.translate(e["cds"]):
            raise SystemExit(f"{key}: domestication would change the protein")
        r = active[key]
        vector, built, clean = assemble_candidates.assemble_one(
            bb, bb_sites, r["candidate_id"], key, con
        )
        if not clean:
            raise SystemExit(f"{key}: vector does not show the expected site pattern")
        r.update({c: str(built[c]) for c in DESIGN_COLS if c in built})
        changed[key], new_con[key] = e, con
        fasta.append(assemble_candidates.fasta_record(r["candidate_id"], key, vector))

    after = {k: "\t".join(r[c] for c in cols) for k, r in active.items()}
    same = [k for k in active if k not in changed and after[k] == before[k]]
    assert len(same) == 200 - len(changed), "an ortholog row outside the applied trims changed"
    cds_record = {**prev_cds, **changed}
    con_record = {**prev_con, **new_con}
    moved = sorted(k for k in prev_cds if cds_record[k] != prev_cds[k])
    assert moved == sorted(changed), f"unexpected CDS changes {moved}"
    summary = {
        "applied": sorted(changed),
        "held_pending_v1": sorted(x["key"] for x in report if x["status"] != "applied"),
        "v1": {x["key"]: x["v1_topology_version"] or "none" for x in report},
        "other_active_orthologs_unchanged": len(same),
        "records": {n: f"processed://{RECORDS}{n}" for n in _names()},
        "record_sha256": {
            "candidate_cds.json": sha(record_text(cds_record)),
            "candidate_constructs.json": sha(record_text(con_record)),
            "changed_vectors.fasta": sha("".join(fasta)),
        },
    }
    return {
        "rows": rows,
        "cols": cols,
        "cds_record": cds_record,
        "con_record": con_record,
        "fasta": fasta,
        "summary": summary,
    }


def _names() -> tuple[str, ...]:
    return ("candidate_cds.json", "candidate_constructs.json", "changed_vectors.fasta")


def record_text(data: dict) -> str:
    return json.dumps(data, indent=0, sort_keys=True)


def published(r: dict) -> dict[str, bytes]:
    return {
        "candidate_cds.json": record_text(r["cds_record"]).encode(),
        "candidate_constructs.json": record_text(r["con_record"]).encode(),
        "changed_vectors.fasta": "".join(r["fasta"]).encode(),
    }


def verify(work: Path) -> str:
    from deliverome_analysis.core.s3 import processed

    r = reproduce(work)
    assert SUMMARY.read_text() == json.dumps(r["summary"], indent=2) + "\n", "summary differs"
    for name, data in published(r).items():
        got = Path(processed.download(RECORDS + name, work / f"rec_{name}")).read_bytes()
        assert got == data, f"processed://{RECORDS}{name} does not reproduce"
    with MANIFEST.open() as fh:
        committed = {x["key"]: x for x in to_design_space(list(csv.DictReader(fh, delimiter="\t")))}
    built = {x["key"]: x for x in r["rows"]}
    for key in r["summary"]["applied"]:
        assert all(committed[key][c] == built[key][c] for c in DESIGN_COLS), key
    s = r["summary"]
    return (
        f"{len(s['applied'])} cyno restart applied ({', '.join(s['applied'])}), "
        f"{len(s['held_pending_v1'])} held pending a v1 row; every restart matches its recorded "
        f"CDS and protein hash; {s['other_active_orthologs_unchanged']} other active L1 "
        "orthologs unchanged; records byte for byte"
    )


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--write", action="store_true")
    ap.add_argument("--publish", action="store_true")
    args = ap.parse_args(argv)
    work = Path(tempfile.mkdtemp(prefix="ortholog-trims-"))
    r = reproduce(work)
    print(json.dumps(r["summary"], indent=2))
    if args.write:
        SUMMARY.write_text(json.dumps(r["summary"], indent=2) + "\n")
        with MANIFEST.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=r["cols"], delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(r["rows"])
        print(f"wrote {SUMMARY.name}, {MANIFEST.name}")
        print(RENUMBER_HINT)
    if args.publish:
        from deliverome_analysis.core.s3 import processed

        if SUMMARY.read_text() != json.dumps(r["summary"], indent=2) + "\n":
            raise SystemExit("summary differs from the rerun; --write and commit first")
        for name, data in published(r).items():
            key = RECORDS + name
            if processed.exists(key):
                raise SystemExit(f"processed://{key} exists; records are never overwritten")
            path = work / f"pub_{name}"
            path.write_bytes(data)
            processed.upload(str(path), key)
            print(f"published processed://{key}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
