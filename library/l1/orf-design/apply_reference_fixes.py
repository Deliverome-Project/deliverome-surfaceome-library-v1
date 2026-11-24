"""Apply the approved L1 reference-form fixes: named references and start-codon replacements.

Two decisions, both implemented by the Claude Code session l1-reviewer-package (subagent),
PR #283:

* **Signal-peptide references** (Becca Carlson, 2026-11-11, "Apply all"; "tIII only for
  now"; reviewed by session type3-risk-review, recorded in PR #282 by session
  l1-nterm-type3): CD200R1L, LRRC37B and CACHD1 move to the UniProt canonical form that
  keeps the signal peptide, already built in L2.
* **Reference-form screen** (Becca Carlson, 2026-11-12, "Send all fixes"; evidence PR #280
  ``REFERENCE_FORM_SCREEN.md``): ADGRF3 moves to its signal-peptide form Q8IZF5; RGMB to
  Q6NW40, the form that starts at the signal peptide instead of 41 residues upstream of it;
  and the approved ``library_cds_boundary_decisions.tsv`` rows replace the annotated
  non-ATG start of SLC22A17 (CTG) and TRPV6 (ACG) with ATG, adding no residue. EPHA6 was
  approved to move to Q9UF33, but no Ensembl 115 transcript encodes Q9UF33 exactly, so it is
  not applied (the decision said not to choose another form).
* **Second round** (decided by session l1-proteoform-audit under Becca's 2026-11-12
  delegation, "use your best judgment"; ``l1_proteoform_decisions.tsv``, PR #300; applied
  here by l1-reviewer-package (subagent), PR #283): LMLN, DISP2 and TSPAN10 move to their UniProt canonical form (the
  APPRIS-only forms are N- or C-terminally truncated Ensembl models); OR8G3's GTG start is
  replaced with ATG (``library_cds_boundary_decisions.tsv``), giving UniProt P0DMU2; and
  EPHA6's CDS is trimmed in frame to the Met that starts Q9UF33's signal peptide, removing
  the 94 upstream codons of the MANE transcript, which gives Q9UF33 exactly. Their ledger
  rows copy the audit row's ``decided_by`` / ``implemented_by`` / ``reason`` /
  ``alternative_rejected`` / ``sources`` verbatim from ``proteoform_decisions_applied.tsv``
  (``proteoform_decisions.py``).

The named references go through the ladder's own mechanism: rows in
``library/l1/proteoform-selection/reference_overrides.tsv`` (pinned as a pipeline input),
then ``select_reference.py --offline`` and ``resolve_candidate_cds.py`` rerun from pinned
inputs. Only the named genes' reference rows and ``lib:`` CDS may change; each new CDS must
encode the named UniProt protein exactly, start ATG and come from the expected transcript.
A start replacement is applied only to an active L1 ``lib:`` row approved in the decisions
table whose current CDS hashes to the row's ``original_cds_sha256``.

The CDS and construct records are versioned: the previous records with the changed ``lib:``
entries (and ADGRF3's, which is replaced in place like MICA's on 2026-11-01), published under
a new dated key. Every other entry is unchanged and the earlier records stay where they are.
``reference_changes.tsv`` is the ledger, one row per change with ``decided_by``,
``implemented_by`` and ``reason``; the manifest rows of the changed members take the new
design columns (transcript, route, lengths, sites), and their identity is kept.

    uv run --locked python apply_reference_fixes.py            # dry run: reproduce, report
    uv run --locked python apply_reference_fixes.py --write    # ledger, summary, manifest
    uv run --locked python apply_reference_fixes.py --publish  # after committing: S3 records
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
import proteoform_decisions as pd_  # noqa: E402
from member_id_allocations import RENUMBER_HINT, to_design_space  # noqa: E402

sys.path.insert(0, str(lp.AP))
import assemble_candidates  # noqa: E402
import candidate_inputs  # noqa: E402
import domesticate_candidates  # noqa: E402
import resolve_candidate_cds as rcd  # noqa: E402
from Bio import SeqIO  # noqa: E402

from deliverome_analysis import library_design as ld  # noqa: E402
from deliverome_analysis import seqtools, topology  # noqa: E402

MANIFEST = HERE / "candidate_manifest_tranches.tsv"
DECISIONS = HERE / "library_cds_boundary_decisions.tsv"
LEDGER = HERE / "reference_changes.tsv"
SUMMARY = HERE / "reference_changes_summary.json"

#: The records this change extends (the 2026-11-10 ortholog adoption's).
PREV_CDS_KEY = (
    "sprint4/s4e9-orf-design-fg-panel/library-candidates/"
    "2026-11-10-ortholog-cds-policy/candidate_cds.json"
)
PREV_CONSTRUCTS_KEY = PREV_CDS_KEY.replace("candidate_cds.json", "candidate_constructs.json")
#: The reference record before this change.
PREV_REFERENCE = "library-design/"
#: The first round's records (2026-11-12): kept, and checked to differ from this record only
#: in the second round's entries.
ROUND1_RECORDS = "sprint4/s4e9-orf-design-fg-panel/library-candidates/2026-11-12-reference-fixes/"
#: Where this change's records go (both rounds, cumulative over the 2026-11-10 records);
#: nothing existed under these prefixes before.
RECORDS = "sprint4/s4e9-orf-design-fg-panel/library-candidates/2026-11-12-reference-fixes-2/"
REFERENCE_RECORDS = (
    "sprint4/s4e9-orf-design-fg-panel/library-design/reference-selection/"
    "2026-11-12-reference-fixes-2/"
)

SP_SWAP = {
    "decided_by": "Becca Carlson, 2026-11-11 (Apply all; tIII only for now); reviewed by "
    "session type3-risk-review; recorded in PR #282 by session l1-nterm-type3",
    "implemented_by": "Claude Code session l1-reviewer-package (subagent), PR #283",
}
SCREEN = {
    "decided_by": "Becca Carlson, 2026-11-12 ('Send all fixes'); evidence PR #280 "
    "REFERENCE_FORM_SCREEN.md",
    "implemented_by": "Claude Code session l1-reviewer-package (subagent), PR #283",
}
#: Second-round rows: attribution and reason come from the audit table at ledger time.
AUDITED = {"decided_by": "see proteoform_decisions_applied.tsv", "implemented_by": ""}
#: gene -> (form, expected Ensembl 115 transcript, attribution, reason). The override rows in
#: reference_overrides.tsv carry the same form, basis and attribution.
NAMED = {
    "CD200R1L": (
        "Q6Q8B3",
        "ENST00000398214",
        SP_SWAP,
        "type III reference lacks the signal peptide of the UniProt canonical; take "
        "the signal-peptide form (N after SP)",
    ),
    "LRRC37B": (
        "Q96QE4",
        "ENST00000341671",
        SP_SWAP,
        "type III reference lacks the signal peptide of the UniProt canonical; take "
        "the signal-peptide form (N after SP)",
    ),
    "CACHD1": (
        "Q5VU97",
        "ENST00000651257",
        SP_SWAP,
        "type III reference lacks the signal peptide of the UniProt canonical; take "
        "the signal-peptide form (N after SP)",
    ),
    "ADGRF3": (
        "Q8IZF5",
        "ENST00000311519",
        SCREEN,
        "reference lacks the signal peptide of the UniProt canonical; take Q8IZF5-1",
    ),
    "RGMB": (
        "Q6NW40",
        None,
        SCREEN,
        "reference RGMB-201 adds 41 residues upstream of the signal peptide; take the "
        "form that starts at it",
    ),
    "LMLN": (
        "Q96KR4",
        "ENST00000330198",
        AUDITED,
        "see proteoform_decisions_applied.tsv",
    ),
    "DISP2": (
        "A7MBM2",
        "ENST00000267889",
        AUDITED,
        "see proteoform_decisions_applied.tsv",
    ),
    "TSPAN10": (
        "Q9H1Z9",
        None,
        AUDITED,
        "see proteoform_decisions_applied.tsv",
    ),
}
NO_ALTERNATIVE = "not recorded (before the 2026-11-12 alternative column)"
#: In-frame N-terminal trims of the existing CDS: key -> spec. The trimmed protein must equal
#: the UniProt form exactly; no sequence is added.
TRIMS = {
    "lib:EPHA6": {
        "form": "Q9UF33",
        "transcript": "ENST00000389672",
        "codons": 94,
        "original_cds_sha256": "297836165d7b5f243d5a36783e85ebb3f86ad3685de0187b1161f1ca4cf70bd4",
        **AUDITED,
        "reason": "see proteoform_decisions_applied.tsv",
        "alternative": "see proteoform_decisions_applied.tsv",
    }
}
#: Start replacements whose result must equal a UniProt protein exactly.
START_EXPECT = {"lib:OR8G3": "P0DMU2"}
#: Approved but not applied, with why (reported, never forced).
NOT_APPLIED: dict[str, str] = {}
START_FIX = {
    **SCREEN,
    "reason": "approved canonical-initiation variant: replace the annotated non-ATG start "
    "codon with ATG; no residue added",
}
#: The entries the second round changed.
ROUND2_KEYS = ("lib:DISP2", "lib:EPHA6", "lib:LMLN", "lib:OR8G3", "lib:TSPAN10")
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
LEDGER_COLS = [
    "candidate_id",
    "key",
    "change",
    "form",
    "old_transcript",
    "new_transcript",
    "old_protein_aa",
    "new_protein_aa",
    "old_cds_sha256",
    "new_cds_sha256",
    "v1_topology_version",
    "v1_label",
    "decided_by",
    "implemented_by",
    "reason",
    "alternative",
    "sources",
]


def sha(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


def unix(path: Path) -> list[str]:
    return Path(path).read_bytes().replace(b"\r\n", b"\n").decode().splitlines()


def protein(entry: dict) -> str:
    return seqtools.translate(entry["cds"]).rstrip("*")


def approved_starts() -> dict[str, dict]:
    with DECISIONS.open() as fh:
        return {
            r["key"]: r
            for r in csv.DictReader(fh, delimiter="\t")
            if r["approved_for_canonical_activation"] == "True"
        }


def start_replaced(entry: dict, decision: dict) -> dict:
    """The entry with its first codon replaced by ATG, after checking it is the reviewed CDS."""
    if decision["action"] != "replace_start_with_ATG":
        raise SystemExit(f"{decision['key']}: approved action is {decision['action']}")
    if sha(entry["cds"]) != decision["original_cds_sha256"]:
        raise SystemExit(f"{decision['key']}: CDS differs from the reviewed original")
    if entry["cds"][:3] == "ATG":
        raise SystemExit(f"{decision['key']}: already starts ATG")
    return {**entry, "cds": "ATG" + entry["cds"][3:], "start_codon_replaced": entry["cds"][:3]}


def trimmed(entry: dict, spec: dict) -> dict:
    """The entry with its first ``codons`` codons removed, after checking it is the reviewed
    CDS and that the cut lands on an ATG."""
    if sha(entry["cds"]) != spec["original_cds_sha256"]:
        raise SystemExit(f"{spec['form']}: CDS differs from the reviewed original")
    if entry["transcript"] != spec["transcript"]:
        raise SystemExit(f"{spec['form']}: transcript is {entry['transcript']}")
    cds = entry["cds"][3 * spec["codons"] :]
    if not cds.startswith("ATG"):
        raise SystemExit(f"{spec['form']}: the trim does not start at an ATG")
    return {**entry, "cds": cds, "trimmed_upstream_codons": spec["codons"]}


def reproduce(work: Path) -> dict:
    from deliverome_analysis.s3 import processed

    uniprot = rcd.load_uniprot("human")
    prev_cds = json.loads(Path(processed.download(PREV_CDS_KEY, work / "pc.json")).read_text())
    prev_con = json.loads(
        Path(processed.download(PREV_CONSTRUCTS_KEY, work / "pk.json")).read_text()
    )
    base = json.loads(Path(processed.download(ld.BASE_CDS_KEY, work / "bc.json")).read_text())

    # 1. reference ladder with the named rows; only the named genes may move.
    out = lp.run_reference(work / "ref")
    rec = unix(lp.record(PREV_REFERENCE + "reference_selection.tsv", work / "rec_ref.tsv"))
    new = unix(out / "reference_selection.tsv")
    assert new[0] == rec[0], "reference_selection.tsv header changed"
    moved = {line.split("\t")[0] for line in set(rec) ^ set(new)}
    assert moved == set(NAMED), f"reference rows changed for {sorted(moved ^ set(NAMED))}"
    rows_by_gene = {line.split("\t")[0]: line.split("\t") for line in new}
    hdr = new[0].split("\t")
    for gene, (form, *_rest) in NAMED.items():
        r = dict(zip(hdr, rows_by_gene[gene], strict=True))
        assert r["reference_rung"] == "named" and r["selected_uniprot_equivalent"] == form, (
            f"{gene}: ladder did not name {form}"
        )

    # 2. designed CDS; every entry not named reproduces the resolver's record.
    cds = lp.run_candidate_cds(work / "cds", out / "reference_selection.tsv")
    named_keys = {f"lib:{g}" for g in NAMED}
    differ = sorted(k for k in base if k not in cds or cds[k]["cds"] != base[k]["cds"])
    assert set(differ) == named_keys, f"CDS changed for {sorted(set(differ) ^ named_keys)}"
    changed: dict[str, dict] = {}
    for gene, (form, want_tx, _att, _why) in NAMED.items():
        key = f"lib:{gene}"
        e = cds[key]
        assert e["variant_residues"] == 0, f"{key}: not an exact encoding"
        assert protein(e) == uniprot[form], f"{key}: does not encode {form}"
        assert e["cds"].startswith("ATG"), f"{key}: CDS does not start ATG"
        if want_tx:
            assert e["transcript"] == want_tx, f"{key}: {e['transcript']} is not {want_tx}"
        changed[key] = e
    # 3. approved start-codon replacements on active L1 references.
    with MANIFEST.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        cols, rows = list(reader.fieldnames or []), to_design_space(list(reader))
    active_l1 = {r["key"]: r for r in rows if r["library"] == ld.LIBRARY_L1 and ld.is_active(r)}
    for key, dec in sorted(approved_starts().items()):
        if key not in active_l1:
            continue
        changed[key] = start_replaced(prev_cds[key], dec)
        if key in START_EXPECT:
            want = START_EXPECT[key]
            assert protein(changed[key]) == uniprot[want], f"{key}: does not encode {want}"
    # 3b. in-frame N-terminal trims to the UniProt start.
    for key, spec in sorted(TRIMS.items()):
        assert key in active_l1, f"{key}: not an active L1 member"
        changed[key] = trimmed(prev_cds[key], spec)
        got, want = protein(changed[key]), uniprot[spec["form"]]
        assert got == want, f"{key}: trimmed protein is not {spec['form']}"

    # 4. constructs and vectors for the changed entries.
    bb = str(SeqIO.read(candidate_inputs.backbone(), "genbank").seq).upper()
    bb_sites = seqtools.find_sites(bb, assemble_candidates.ENZYMES)
    new_con, fasta, ledger, design_rows = {}, [], [], {}
    audited = {k for k in pd_.by_key() if k.startswith("lib:")}
    for key in sorted(changed):
        e = changed[key]
        con = domesticate_candidates.domesticate_one(e)
        if con is None or seqtools.translate(con["construct"]) != seqtools.translate(e["cds"]):
            raise SystemExit(f"{key}: domestication would change the protein")
        r = active_l1[key]
        cid = r["candidate_id"]
        vector, built, clean = assemble_candidates.assemble_one(bb, bb_sites, cid, key, con)
        if not clean:
            raise SystemExit(f"{key}: vector does not show the expected site pattern")
        new_con[key] = con
        fasta.append(assemble_candidates.fasta_record(cid, key, vector))
        design_rows[key] = {c: str(built[c]) for c in DESIGN_COLS if c in built}
        hit = topology.for_sequence(protein(e), tool=topology.TOOL_V1)
        gene = key.split(":", 1)[1]
        alt = NO_ALTERNATIVE
        if gene in NAMED:
            form, _tx, att, why = NAMED[gene]
            kind = "named_reference"
        elif key in TRIMS:
            spec = TRIMS[key]
            form, att, why, alt = spec["form"], spec, spec["reason"], spec["alternative"]
            kind = f"trim_{spec['codons']}_upstream_codons"
        else:
            dec = approved_starts()[key]
            if dec["decided_by"].startswith("Becca"):
                form, att, why = "", START_FIX, START_FIX["reason"]
            else:
                form, att, why = START_EXPECT.get(key, ""), dec, dec["decision_reason"]
                alt = dec["alternative"]
            kind = f"start_{e['start_codon_replaced']}_to_ATG"
        ledger.append(
            {
                "candidate_id": cid,
                "key": key,
                "change": kind,
                "form": form,
                "old_transcript": prev_cds[key]["transcript"],
                "new_transcript": e["transcript"],
                "old_protein_aa": len(protein(prev_cds[key])),
                "new_protein_aa": len(protein(e)),
                "old_cds_sha256": sha(prev_cds[key]["cds"]),
                "new_cds_sha256": sha(e["cds"]),
                "v1_topology_version": hit["topology_version"] if hit else "MISSING",
                "v1_label": hit["deeptmhmm_label"] if hit else "",
                "decided_by": att["decided_by"],
                "implemented_by": att["implemented_by"],
                "reason": why,
                "alternative": alt,
                "sources": "",
                # Second-round rows carry the audit's own attribution, verbatim.
                **(pd_.attribution(key) if key in audited else {}),
            }
        )
    assert set(audited) <= set(changed), f"audited but not applied: {set(audited) - set(changed)}"
    for r in rows:
        if r["key"] in design_rows:
            r.update(design_rows[r["key"]])
    cds_record = {**prev_cds, **changed}
    con_record = {**prev_con, **new_con}
    assert set(cds_record) == set(prev_cds) and set(con_record) == set(prev_con)
    summary = {
        "named_references": {g: NAMED[g][0] for g in NAMED},
        "not_applied": NOT_APPLIED,
        "start_codon_replacements": sorted(
            k for k in changed if k.split(":")[1] not in NAMED and k not in TRIMS
        ),
        "n_terminal_trims": {k: TRIMS[k]["form"] for k in sorted(TRIMS)},
        "changed_entries": sorted(changed),
        "v1_missing": [x["key"] for x in ledger if x["v1_topology_version"] == "MISSING"],
        "records": {
            "candidate_cds": f"processed://{RECORDS}candidate_cds.json",
            "candidate_constructs": f"processed://{RECORDS}candidate_constructs.json",
            "changed_vectors": f"processed://{RECORDS}changed_vectors.fasta",
            "reference_selection": f"processed://{REFERENCE_RECORDS}reference_selection.tsv",
        },
        "record_sha256": {
            "candidate_cds.json": sha(record_text(cds_record)),
            "candidate_constructs.json": sha(record_text(con_record)),
            "changed_vectors.fasta": sha("".join(fasta)),
            "reference_selection.tsv": hashlib.sha256(
                (out / "reference_selection.tsv").read_bytes()
            ).hexdigest(),
        },
    }
    return {
        "reference_dir": out,
        "cds_record": cds_record,
        "con_record": con_record,
        "fasta": fasta,
        "ledger": ledger,
        "rows": rows,
        "cols": cols,
        "summary": summary,
    }


def record_text(data: dict) -> str:
    return json.dumps(data, indent=0, sort_keys=True)


def tsv_text(rows: list[dict], cols: list[str]) -> str:
    lines = ["\t".join(cols)] + ["\t".join(str(r.get(c, "")) for c in cols) for r in rows]
    return "\n".join(lines) + "\n"


def published(r: dict) -> dict[str, tuple[str, bytes]]:
    out = r["reference_dir"]
    return {
        "candidate_cds.json": (RECORDS, record_text(r["cds_record"]).encode()),
        "candidate_constructs.json": (RECORDS, record_text(r["con_record"]).encode()),
        "changed_vectors.fasta": (RECORDS, "".join(r["fasta"]).encode()),
        **{
            name: (REFERENCE_RECORDS, (out / name).read_bytes())
            for name in (
                "reference_selection.tsv",
                "reference_disagreements.tsv",
                "reference_selection_summary.json",
            )
        },
    }


def verify(work: Path) -> str:
    """Rerun from pinned inputs; the ledger, summary and every published record must match,
    and each changed member's manifest design columns must be the rebuilt ones."""
    from deliverome_analysis.s3 import processed

    r = reproduce(work)
    assert LEDGER.read_text() == tsv_text(r["ledger"], LEDGER_COLS), f"{LEDGER.name} differs"
    assert SUMMARY.read_text() == json.dumps(r["summary"], indent=2) + "\n", "summary differs"
    for name, (prefix, data) in published(r).items():
        got = Path(processed.download(prefix + name, work / f"rec_{name}")).read_bytes()
        assert got == data, f"processed://{prefix}{name} does not reproduce"
    with MANIFEST.open() as fh:
        committed = {x["key"]: x for x in to_design_space(list(csv.DictReader(fh, delimiter="\t")))}
    built = {x["key"]: x for x in r["rows"]}
    for x in r["ledger"]:
        assert all(committed[x["key"]][c] == built[x["key"]][c] for c in DESIGN_COLS), x["key"]
    # The first round's record differs from this one only in this round's entries.
    round1 = json.loads(
        Path(
            processed.download(ROUND1_RECORDS + "candidate_cds.json", work / "r1.json")
        ).read_text()
    )
    second = sorted(k for k in r["cds_record"] if round1.get(k) != r["cds_record"][k])
    assert set(second) == set(ROUND2_KEYS), f"round-two entries are {second}"
    s = r["summary"]
    return (
        f"{len(s['named_references'])} named references and "
        f"{len(s['start_codon_replacements'])} start-codon replacements and "
        f"{len(s['n_terminal_trims'])} N-terminal trim reproduce (round two: "
        f"{', '.join(k.split(':')[1] for k in ROUND2_KEYS)}); records and "
        f"ledger byte for byte; not applied: {', '.join(s['not_applied']) or 'none'}"
    )


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--write", action="store_true")
    ap.add_argument("--publish", action="store_true")
    args = ap.parse_args(argv)
    work = Path(tempfile.mkdtemp(prefix="reference-fixes-"))
    r = reproduce(work)
    print(json.dumps(r["summary"], indent=2))
    for x in r["ledger"]:
        print(
            f"  {x['candidate_id']} {x['key']}: {x['old_transcript']} {x['old_protein_aa']} aa -> "
            f"{x['new_transcript']} {x['new_protein_aa']} aa ({x['change']}; v1 "
            f"{x['v1_label']} {x['v1_topology_version']})"
        )
    if args.write:
        LEDGER.write_text(tsv_text(r["ledger"], LEDGER_COLS))
        SUMMARY.write_text(json.dumps(r["summary"], indent=2) + "\n")
        (lp.AP / "reference_selection_summary.json").write_bytes(
            (r["reference_dir"] / "reference_selection_summary.json").read_bytes()
        )
        with MANIFEST.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=r["cols"], delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(r["rows"])
        print(f"wrote {LEDGER.name}, {SUMMARY.name}, {MANIFEST.name}")
        print(RENUMBER_HINT)
    if args.publish:
        from deliverome_analysis.s3 import processed

        if LEDGER.read_text() != tsv_text(r["ledger"], LEDGER_COLS):
            raise SystemExit("ledger differs from the rerun; --write and commit first")
        for name, (prefix, data) in published(r).items():
            key = prefix + name
            if processed.exists(key):
                raise SystemExit(f"processed://{key} exists; records are never overwritten")
            path = work / f"pub_{name}"
            path.write_bytes(data)
            processed.upload(str(path), key)
            print(f"published processed://{key}")
    if not (args.write or args.publish):
        print("dry run -- --write writes the ledger and manifest; commit, then --publish")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
