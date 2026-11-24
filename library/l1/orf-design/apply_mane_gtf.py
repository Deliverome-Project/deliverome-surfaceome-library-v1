"""Rerun the L1 reference ladder with MANE Select read from the GTF, and carry the switches through.

The ladder took MANE from the APPRIS e115v50 file's MANE column, which omits 775 of Ensembl 115's
MANE Select transcripts. 66 library genes therefore looked "APPRIS only", and in 58 of them the
library builds a protein that is not the MANE Select one (``PROTEOFORM_REVIEW_2026-11-22.md``).
**Decided by Becca Carlson, 2026-11-22** ("Yes, in this build", to rerunning the ladder with GTF
MANE). Implemented by Claude Code session l1-review-export-scope (subagent).

Steps, each from pinned inputs:

1. **evidence** -- for every gene whose reference the GTF MANE would move, a three-way look at
   ours (the current designed protein), MANE and the UniProt canonical. PeptideAtlas build 607
   peptides observed for *any* reviewed accession of the gene (``REVIEW_CACHE``), counted for a
   form when contained in it and not the other and proteotypic (``proteoform_review.py``). The
   call (``three_way_call``) is the corrected policy: MANE, unless observed evidence clearly
   favors ours; where both forms are observed, or ours only weakly, the current reference is
   kept and the gene listed as an open decision. Written to ``mane_gtf_decisions_2026-11-22.tsv``.
2. **overrides** -- named references for the genes kept on their current form, and for the two
   non-MANE switches the review recommended (PTPRM, LY6H), appended to
   ``library/l1/proteoform-selection/reference_overrides.tsv`` and pinned as a new key.
3. **build** -- ``select_reference.py --offline`` over the gene list: only the 66 genes' rows may
   change, and only the switched genes' proteins. ``resolve_candidate_cds.py``: only their
   ``lib:`` CDS may change. Each new CDS is domesticated and dropped into the p0039 backbone as
   ``apply_reference_fixes.py`` does; the manifest rows keep their identity and take the new
   design columns (the RGMB precedent). Every new protein must have a DeepTMHMM v1 row.
4. **publish** -- new dated records under ``TAG``: the reference record, and the CDS and
   construct records (the current ones with the switched ``lib:`` entries replaced).

    uv run --locked python apply_mane_gtf.py evidence            # decisions table
    uv run --locked python apply_mane_gtf.py overrides --publish  # pin the named rows
    uv run --locked python apply_mane_gtf.py build [--write] [--publish]
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
import sys
import tempfile
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import library_pipeline as lp  # noqa: E402
import proteoform_review as pr  # noqa: E402

sys.path.insert(0, str(lp.AP))
import assemble_candidates  # noqa: E402
import candidate_inputs  # noqa: E402
import domesticate_candidates  # noqa: E402
from Bio import SeqIO  # noqa: E402
from select_principal import load_uniprot_proteoforms  # noqa: E402

from deliverome_analysis import library_design as ld  # noqa: E402
from deliverome_analysis import seqtools, topology  # noqa: E402

TAG = "2026-11-22-mane-gtf"
DECISIONS = HERE / "mane_gtf_decisions_2026-11-22.tsv"
LEDGER = HERE / "mane_gtf_changes_2026-11-22.tsv"
SUMMARY = HERE / "mane_gtf_changes_2026-11-22_summary.json"
MANIFEST = HERE / "candidate_manifest_tranches.tsv"
OVERRIDES = lp.AP / "reference_overrides.tsv"
OVERRIDES_KEY = f"{lp.PREFIX}library-design/pipeline-inputs/{TAG}/reference_overrides.tsv"
REF_PREFIX = f"{lp.PREFIX}library-design/reference-selection/{TAG}/"
CAND_PREFIX = f"{lp.PREFIX}library-candidates/{TAG}/"
#: The review's PeptideAtlas cache: every reviewed accession of the genes in question.
REVIEW_CACHE = (
    "peptideatlas/build-607/queries-2026-11-22-mane-review/peptideatlas_build607_pools.json"
)
REVIEW_CACHE_SHA = "60ce32deb51168fe94c732f9a73b203befb5bf29793d0f4a33ec9b8a741823cc"
#: A peptide "seen" is in at least this many samples; "strong" is at least this many such
#: peptides (the ladder's MIN_SAMPLES and STRONG_EVIDENCE).
MIN_SAMPLES, STRONG = 10, 3

DECIDED = (
    "Becca Carlson, 2026-11-22 ('Yes, in this build': rerun the reference ladder with GTF MANE; "
    "three-way calls delegated: MANE unless evidence is strong for the other form, ambiguous "
    "ones kept and listed)"
)
IMPLEMENTED = "Claude Code session l1-review-export-scope (subagent)"
#: Non-MANE switches the proteoform review recommended (PROTEOFORM_REVIEW_2026-11-22.md).
NONMANE = {
    "PTPRM": (
        "P28827",
        "ENST00000332175",
        "evidence_for_appris",
        "MANE (P28827-2) and the APPRIS principal (P28827) disagree; the APPRIS form is a pure "
        "deletion, so the ladder's region test saw no evidence. Counting junction-spanning "
        "peptides the canonical has 7 distinguishing peptides (4 in >= 10 samples, best 28) and "
        "MANE none, which meets the ladder's bar for overruling MANE",
    ),
    "LY6H": (
        "O94772",
        "ENST00000430474",
        "upstream_extension",
        "MANE (O94772-2, ENST00000342752) starts 21 aa upstream of the signal peptide; the "
        "canonical, the APPRIS principal, starts at it and gives the same mature GPI protein. "
        "Same case as RGMB (2026-11-12)",
    ),
}


def sha(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


def protein(entry: dict) -> str:
    return seqtools.translate(entry["cds"]).rstrip("*")


def unix(path: Path) -> list[str]:
    return path.read_bytes().replace(b"\r\n", b"\n").decode().splitlines()


def three_way_call(ours: tuple[int, int, int], mane: tuple[int, int, int]) -> str:
    """``switch_to_mane`` / ``keep_current_evidence`` / ``keep_current_open``.

    Each argument is (distinguishing peptides, of them seen in >= MIN_SAMPLES samples, best
    sample count) for that form against the other. MANE is the default. Ours is kept when it
    has at least STRONG peptides seen in >= MIN_SAMPLES samples and MANE has none so seen (the
    ladder's bar); kept and left open when ours has any such peptide but not that clear a case
    (both forms observed, or ours only one or two peptides deep).
    """
    if ours[1] == 0:
        return "switch_to_mane"
    if ours[1] >= STRONG and mane[1] == 0:
        return "keep_current_evidence"
    return "keep_current_open"


def stage(work: Path) -> dict:
    names = [
        "Homo_sapiens.GRCh38.pep.all.fa.gz",
        "uniprot_human_reviewed_isoforms.tsv.gz",
        "mane_select_ensembl115.tsv",
    ]
    return lp.stage_inputs(names, work / "in")


def evidence(work: Path) -> list[dict]:
    from deliverome_analysis.s3 import external, processed

    got = stage(work)
    pep = {}
    import gzip

    cur = ""
    with gzip.open(got["Homo_sapiens.GRCh38.pep.all.fa.gz"], "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                cur = next(f for f in line.split() if f.startswith("transcript:"))
                cur = cur.split(":")[1].split(".")[0]
                pep[cur] = []
            else:
                pep[cur].append(line.strip())
    pep = {k: "".join(v).rstrip("*") for k, v in pep.items()}
    with got["mane_select_ensembl115.tsv"].open() as fh:
        mane_by_name = {
            r["gene_name"]: r["transcript_id"] for r in csv.DictReader(fh, delimiter="\t")
        }
    uniprot = load_uniprot_proteoforms(got["uniprot_human_reviewed_isoforms.tsv.gz"])
    cds = ld.load_designed_cds()
    refsel = list(
        csv.DictReader(
            open(processed.download(lp.PREFIX + lp.REFERENCE_RECORD + "reference_selection.tsv")),
            delimiter="\t",
        )
    )
    cache = Path(external.download(REVIEW_CACHE, work / "atlas.json"))
    assert lp.sha256(cache) == REVIEW_CACHE_SHA, "review cache differs from its pin"
    blob = json.loads(cache.read_text())
    allp = [(s, p["seq"]) for s, recs in uniprot.items() for p in recs]
    rows = []
    for r in refsel:
        g = r["symbol"]
        key = f"lib:{g}"
        if r["mane_enst"] or g not in mane_by_name or key not in cds:
            continue
        m_enst = mane_by_name[g]
        mane, ours = pep.get(m_enst, ""), protein(cds[key])
        if not mane or mane == ours or r["reference_rung"] == "named":
            continue
        can = next(
            (p["seq"] for p in uniprot.get(g, []) if p["accession"] == r["current_accession"]), ""
        )
        accs = [p["accession"] for p in uniprot.get(g, [])]
        missing = [a for a in accs if a not in blob["index"]]
        assert not missing, f"{g}: review cache lacks {missing}"
        obs: dict[str, tuple[int, int]] = {}
        for a in accs:
            for p_, v in blob["pools"][blob["index"][a]].items():
                if p_ not in obs or v[1] > obs[p_][1]:
                    obs[p_] = (v[0], v[1])

        def proteotypic(sub: str, g=g) -> bool:
            return not any(sub in s for sy, s in allp if sy != g)

        def stats(hits) -> tuple[int, int, int]:
            return (
                len(hits),
                sum(1 for h in hits if h[2] >= MIN_SAMPLES),
                max((h[2] for h in hits), default=0),
            )

        so = pr.distinguishing(ours, mane, obs, proteotypic)
        sm = pr.distinguishing(mane, ours, obs, proteotypic)
        o, m = stats(so), stats(sm)
        acc_of = {p["seq"]: p["accession"] for p in uniprot.get(g, [])}
        rows.append(
            {
                "gene_symbol": g,
                "current_source": r["selected_source"],
                "current_transcript": cds[key]["transcript"],
                "current_form": acc_of.get(ours, ""),
                "current_aa": len(ours),
                "mane_select": m_enst,
                "mane_form": acc_of.get(mane, ""),
                "mane_aa": len(mane),
                "canonical_acc": r["current_accession"],
                "canonical_aa": len(can),
                "current_is_canonical": int(ours == can),
                "mane_is_canonical": int(mane == can),
                "ours_peptides": o[0],
                "ours_ge_10_samples": o[1],
                "ours_best_samples": o[2],
                "ours_examples": ";".join(h[0] for h in so[:3]),
                "mane_peptides": m[0],
                "mane_ge_10_samples": m[1],
                "mane_best_samples": m[2],
                "mane_examples": ";".join(h[0] for h in sm[:3]),
                "call": three_way_call(o, m),
            }
        )
    return rows


def write_tsv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def read_tsv(path: Path) -> list[dict]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def override_rows() -> list[dict]:
    """The named rows this change adds: kept forms (open or evidence) and the non-MANE switches."""
    out = []
    for d in read_tsv(DECISIONS):
        if d["call"] == "switch_to_mane":
            continue
        form = d["current_form"] or d["current_transcript"]
        why = (
            f"GTF MANE Select {d['mane_select']} ({d['mane_form'] or 'no UniProt entry'}, "
            f"{d['mane_aa']} aa) would replace the current {form} ({d['current_aa']} aa). "
            f"Distinguishing PeptideAtlas peptides: current {d['ours_peptides']} "
            f"({d['ours_ge_10_samples']} in >= 10 samples, best {d['ours_best_samples']}), "
            f"MANE {d['mane_peptides']} ({d['mane_ge_10_samples']}, best {d['mane_best_samples']})"
        )
        if d["call"] == "keep_current_evidence":
            basis, why = "evidence_for_current", why + "; the evidence favors the current form"
        else:
            basis, why = (
                "open_decision_keep_current",
                why
                + ("; both forms observed or the current one only weakly: kept pending a decision"),
            )
        out.append(
            {
                "gene": d["gene_symbol"],
                "form": form,
                "basis": basis,
                "note": why,
                "decided_by": DECIDED,
                "implemented_by": IMPLEMENTED,
                "alternative": f"MANE Select {d['mane_select']}",
            }
        )
    for g, (form, tx, basis, why) in NONMANE.items():
        out.append(
            {
                "gene": g,
                "form": form,
                "basis": basis,
                "note": f"{why}; expected transcript {tx}",
                "decided_by": "Becca Carlson, 2026-11-22 (approved the proteoform review's "
                "switches with the GTF MANE rerun)",
                "implemented_by": IMPLEMENTED,
                "alternative": "the ladder's MANE pick",
            }
        )
    return out


def overrides(publish: bool) -> int:
    rows = read_tsv(OVERRIDES)
    cols = list(rows[0])
    have = {r["gene"] for r in rows}
    new = [r for r in override_rows() if r["gene"] not in have]
    if new:
        with OVERRIDES.open("a", newline="") as fh:
            csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n").writerows(new)
        print(f"appended {len(new)} row(s) to {OVERRIDES.name}: {[r['gene'] for r in new]}")
    want = {r["gene"]: r for r in override_rows()}
    now = {r["gene"]: r for r in read_tsv(OVERRIDES)}
    assert all(now[g]["form"] == r["form"] for g, r in want.items()), "an override row differs"
    if publish:
        from deliverome_analysis.s3 import processed

        spec = json.loads(lp.INPUTS.read_text())
        entry = spec["files"]["reference_overrides.tsv"]
        if entry["sha256"] == lp.sha256(OVERRIDES):
            print("already pinned")
            return 0
        if processed.exists(OVERRIDES_KEY):
            raise SystemExit(f"processed://{OVERRIDES_KEY} exists")
        processed.upload(str(OVERRIDES), OVERRIDES_KEY)
        entry.update(
            {
                "source": f"processed:{OVERRIDES_KEY}",
                "sha256": lp.sha256(OVERRIDES),
                "origin": "library/l1/proteoform-selection/reference_overrides.tsv as of "
                f"{TAG}: the 2026-11-12 round-two rows unchanged, plus the GTF-MANE rerun's "
                "named rows (apply_mane_gtf.py; mane_gtf_decisions_2026-11-22.tsv)",
            }
        )
        lp.INPUTS.write_text(json.dumps(spec, indent=1, ensure_ascii=False) + "\n")
        print(f"published processed://{OVERRIDES_KEY} and repinned")
    return 0


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


def _ledger_row(
    k: str, g: str, e: dict, r: dict, p: str, cur: dict, hit: dict, decisions: dict, ref_rows: dict
) -> dict:
    """One switched reference's ledger row: old and new CDS, v1 topology, evidence, who and why."""
    d = decisions.get(g, {})
    return {
        "candidate_id": r["candidate_id"],
        "design_id": r.get("design_id", ""),
        "key": k,
        "change": "named_reference" if g in NONMANE else "mane_select_from_gtf",
        "form": ref_rows[g]["selected_uniprot_equivalent"],
        "old_transcript": cur[k]["transcript"],
        "new_transcript": e["transcript"],
        "old_protein_aa": len(protein(cur[k])),
        "new_protein_aa": len(p),
        "old_cds_sha256": sha(cur[k]["cds"]),
        "new_cds_sha256": sha(e["cds"]),
        "new_protein_sha256": topology.sequence_sha256(p),
        "v1_topology_version": hit["topology_version"],
        "v1_label": hit["deeptmhmm_label"],
        "v1_tm": int(hit["tm_helix_count"] or 0),
        "evidence_current": (
            f"{d['ours_peptides']}/{d['ours_ge_10_samples']}/{d['ours_best_samples']}" if d else ""
        ),
        "evidence_mane": (
            f"{d['mane_peptides']}/{d['mane_ge_10_samples']}/{d['mane_best_samples']}" if d else ""
        ),
        "decided_by": DECIDED if g not in NONMANE else "Becca Carlson, 2026-11-22",
        "implemented_by": IMPLEMENTED,
        "reason": NONMANE[g][3]
        if g in NONMANE
        else "MANE/APPRIS disagreement hidden by the APPRIS file's MANE column; the "
        "policy default is MANE and the current form has no peptide seen in >= 10 samples",
    }


def _retired_duplicates(dups: list[tuple[str, str]], active: dict) -> dict:
    """Isoform members whose protein a switched reference now is: retired, never deleted."""
    retire = {}
    for k, dup in dups:
        row = next(x for x in active.values() if x["candidate_id"] == dup)
        if row["set"] != "isoform" or row["key"].split(":")[1] != k.split(":", 1)[1]:
            raise SystemExit(f"{k}: protein duplicates {dup} ({row['key']})")
        retire[dup] = (
            f"duplicates_reference:{k} (its reference moved to this protein, GTF MANE rerun "
            "2026-11-22; isoform slot to be refilled by the panel redraw)"
        )
    return retire


def build(work: Path, snapshot_key: str | None) -> dict:
    from deliverome_analysis.s3 import processed

    if snapshot_key:
        topology.SNAPSHOT_KEY = snapshot_key  # a supplement not yet adopted in topology.py
    decisions = {d["gene_symbol"]: d for d in read_tsv(DECISIONS)}
    switch = {g for g, d in decisions.items() if d["call"] == "switch_to_mane"} | set(NONMANE)
    gap = set(decisions) | {
        r["gene_symbol"] for r in read_tsv(HERE / "proteoform_review_2026-11-22_mane_gap.tsv")
    }

    # 1. reference ladder
    out = lp.run_reference(work / "ref")
    rec = unix(lp.record(lp.REFERENCE_RECORD + "reference_selection.tsv", work / "rec.tsv"))
    new = unix(out / "reference_selection.tsv")
    assert new[0] == rec[0], "reference_selection.tsv header changed"
    moved = {ln.split("\t")[0] for ln in set(rec) ^ set(new)}
    assert moved <= gap | set(NONMANE), f"rows changed outside the MANE gap: {moved - gap}"
    # By protein hash: a gene kept on its form changes source (named) but not protein.
    seq_old = {
        ln.split("\t")[0]: ln.split("\t")[3]
        for ln in unix(lp.AP / "reference_selection_sequences.tsv")
    }
    seq_new = {
        ln.split("\t")[0]: ln.split("\t")[3]
        for ln in unix(out / "reference_selection_sequences.tsv")
    }
    assert set(seq_old) == set(seq_new), "gene list differs from the committed sequences file"
    prot_moved = {g for g in seq_new if seq_new[g] != seq_old[g]}
    assert prot_moved == switch, (
        f"proteins changed {sorted(prot_moved ^ switch)} differently from the decisions"
    )
    hdr = new[0].split("\t")
    ref_rows = {ln.split("\t")[0]: dict(zip(hdr, ln.split("\t"), strict=True)) for ln in new}
    for g, (form, *_r) in NONMANE.items():
        assert ref_rows[g]["reference_rung"] == "named", f"{g}: not named"
        assert ref_rows[g]["selected_uniprot_equivalent"] == form, f"{g}: not {form}"

    # 2. designed CDS: only the switched genes' lib: entries move
    cds_new = lp.run_candidate_cds(work / "cds_new", out / "reference_selection.tsv")
    cds_old = lp.run_candidate_cds(work / "cds_old")
    keys = {f"lib:{g}" for g in switch}
    differ = {k for k in set(cds_old) | set(cds_new) if cds_old.get(k) != cds_new.get(k)}
    # A gene whose row now carries a MANE transcript resolves through it ("named transcript")
    # instead of by protein match: same CDS, only the route moves (7 genes, among them the
    # 2026-11-12 named references CACHD1, DISP2, RGMB and TSPAN10).
    route_only = {
        k
        for k in differ - keys
        if k in cds_old
        and k in cds_new
        and {f: v for f, v in cds_old[k].items() if f != "route"}
        == {f: v for f, v in cds_new[k].items() if f != "route"}
    }
    assert differ - route_only == keys, (
        f"resolver output moved for {sorted(differ - route_only ^ keys)}"
    )
    import apply_reference_fixes as fixes

    # Those named references' record entries equal the resolver's (verify_candidate_cds), so
    # they take the new route; every other route-only entry keeps its recorded bytes.
    reroute = {k for k in route_only if k.split(":", 1)[1] in fixes.NAMED}
    cur = json.loads(Path(processed.download(ld.CDS_KEY, work / "cur_cds.json")).read_text())
    con_cur = json.loads(
        Path(processed.download(ld.CONSTRUCTS_KEY, work / "cur_con.json")).read_text()
    )
    # A switched gene must not carry a later fix on its CDS (a start replacement, trim or
    # restoration): the new entry would silently drop it.
    for k in keys:
        assert cur[k] == cds_old[k], f"{k}: current record carries a fix on top of the resolver"
    uniprot = load_uniprot_proteoforms(stage(work)["uniprot_human_reviewed_isoforms.tsv.gz"])
    for g in switch:
        e = cds_new[f"lib:{g}"]
        assert e["variant_residues"] == 0, f"{g}: not an exact encoding"
        assert e["cds"].startswith("ATG"), f"{g}: does not start ATG"
        want = ref_rows[g]["selected_uniprot_equivalent"]
        if want:
            seq = next(p["seq"] for p in uniprot[g] if p["accession"] == want)
            assert protein(e) == seq, f"{g}: does not encode {want}"
        if g in NONMANE:
            assert e["transcript"] == NONMANE[g][1], f"{g}: {e['transcript']}"
        else:
            assert e["transcript"] == decisions[g]["mane_select"], f"{g}: not its MANE transcript"

    # 3. topology, constructs, manifest rows
    with MANIFEST.open() as fh:
        mrows = list(csv.DictReader(fh, delimiter="\t"))
    active = {r["key"]: r for r in mrows if r["library"] == ld.LIBRARY_L1 and ld.is_active(r)}
    proteins_active = {
        protein(cur[r["key"]]): r["candidate_id"] for r in active.values() if r["key"] in cur
    }
    bb = str(SeqIO.read(candidate_inputs.backbone(), "genbank").seq).upper()
    bb_sites = seqtools.find_sites(bb, assemble_candidates.ENZYMES)
    missing_v1, ledger, design_rows, new_con, fasta = [], [], {}, {}, []
    dups: list[tuple[str, str]] = []
    for k in sorted(keys):
        g = k.split(":", 1)[1]
        e, r = cds_new[k], active[k]
        p = protein(e)
        dup = proteins_active.get(p)
        if dup is not None and dup != r["candidate_id"]:
            dups.append((k, dup))
        hit = topology.for_sequence(p, tool=topology.TOOL_V1)
        if hit is None:
            did = r.get("design_id") or r["candidate_id"]
            missing_v1.append(f">{did}|{g}|human|{topology.sequence_sha256(p)}\n{p}\n")
            continue
        con = domesticate_candidates.domesticate_one(e)
        if con is None or seqtools.translate(con["construct"]) != seqtools.translate(e["cds"]):
            raise SystemExit(f"{k}: domestication would change the protein")
        vector, built, clean = assemble_candidates.assemble_one(
            bb, bb_sites, r["candidate_id"], k, con
        )
        if not clean:
            raise SystemExit(f"{k}: vector does not show the expected site pattern")
        new_con[k] = con
        fasta.append(assemble_candidates.fasta_record(r["candidate_id"], k, vector))
        design_rows[k] = {c: str(built[c]) for c in DESIGN_COLS if c in built}
        ledger.append(_ledger_row(k, g, e, r, p, cur, hit, decisions, ref_rows))
    # A switched reference whose new protein is already an active isoform-panel member: the
    # protein stays in the library as the reference, so the isoform member is retired (never
    # deleted or renumbered) and its panel slot is left for the isoform redraw to refill, as the
    # redraw's own duplicate rule would. Any other duplicate stops the build.
    retire = _retired_duplicates(dups, active)
    if missing_v1:
        path = HERE / "topology_runs" / f"{TAG}_missing_v1.fasta"
        path.write_text("".join(missing_v1))
        raise SystemExit(
            f"{len(missing_v1)} new protein(s) have no DeepTMHMM v1 row; wrote {path.name}: run "
            "the modal sweep, publish_topology_supplement.py publish, then rerun with "
            "--snapshot-key"
        )
    assert all(cur[k] == cds_old[k] for k in reroute), "a named reference's entry moved"
    cds_record = {**cur, **{k: cds_new[k] for k in keys | reroute}}
    con_record = {**con_cur, **new_con}
    assert set(cds_record) == set(cur) and set(con_record) == set(con_cur)
    return {
        "reference_dir": out,
        "switch": sorted(switch),
        "ledger": ledger,
        "design_rows": design_rows,
        "cds_record": cds_record,
        "reroute": sorted(reroute),
        "retire": retire,
        "route_only": sorted(route_only),
        "con_record": con_record,
        "fasta": fasta,
        "ref_rows": ref_rows,
    }


def record_text(data: dict) -> str:
    return json.dumps(data, indent=0, sort_keys=True)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("evidence")
    o = sub.add_parser("overrides")
    o.add_argument("--publish", action="store_true")
    b = sub.add_parser("build")
    b.add_argument("--write", action="store_true")
    b.add_argument("--publish", action="store_true")
    b.add_argument("--snapshot-key", help="topology snapshot to read (a new supplement)")
    a = ap.parse_args(argv)
    work = Path(tempfile.mkdtemp(prefix="mane-gtf-"))
    if a.cmd == "evidence":
        rows = evidence(work)
        write_tsv(DECISIONS, rows)
        print(f"wrote {DECISIONS.name}: {dict(Counter(r['call'] for r in rows))}")
        return 0
    if a.cmd == "overrides":
        return overrides(a.publish)
    r = build(work, a.snapshot_key)
    for x in r["ledger"]:
        print(
            f"  {x['candidate_id']} {x['key']}: {x['old_transcript']} {x['old_protein_aa']} aa -> "
            f"{x['new_transcript']} {x['new_protein_aa']} aa ({x['form'] or 'no UniProt'}; v1 "
            f"{x['v1_label']} tm={x['v1_tm']})"
        )
    summary = {
        "tag": TAG,
        "switched": len(r["switch"]),
        "switched_genes": r["switch"],
        "isoform_members_retired": r["retire"],
        "kept_open": sorted(
            d["gene_symbol"] for d in read_tsv(DECISIONS) if d["call"] == "keep_current_open"
        ),
        "kept_on_evidence": sorted(
            d["gene_symbol"] for d in read_tsv(DECISIONS) if d["call"] == "keep_current_evidence"
        ),
        "selected_source": dict(
            Counter(x["selected_source"] for x in r["ref_rows"].values() if x["symbol"] != "symbol")
        ),
        "decision": dict(
            Counter(x["decision"] for x in r["ref_rows"].values() if x["symbol"] != "symbol")
        ),
        "changes_current_reference": sum(
            1 for x in r["ref_rows"].values() if x["changes_current_reference"] == "True"
        ),
        "records": {
            "reference_selection": f"processed://{REF_PREFIX}reference_selection.tsv",
            "candidate_cds": f"processed://{CAND_PREFIX}candidate_cds.json",
            "candidate_constructs": f"processed://{CAND_PREFIX}candidate_constructs.json",
        },
        "record_sha256": {
            "candidate_cds.json": sha(record_text(r["cds_record"])),
            "candidate_constructs.json": sha(record_text(r["con_record"])),
            "switched_vectors.fasta": sha("".join(r["fasta"])),
            "reference_selection.tsv": lp.sha256(r["reference_dir"] / "reference_selection.tsv"),
        },
    }
    print(json.dumps({k: v for k, v in summary.items() if k != "switched_genes"}, indent=1))
    if a.write:
        write_tsv(LEDGER, r["ledger"])
        SUMMARY.write_text(json.dumps(summary, indent=1) + "\n")
        for name in ("reference_selection_sequences.tsv", "reference_selection_summary.json"):
            shutil.copy(r["reference_dir"] / name, lp.AP / name)
        # The resolver's coverage lines for the switched keys.
        cov = {
            ln.split("\t")[0]: ln
            for ln in unix(work / "cds_new" / "candidate_cds_coverage.tsv")
            if ln.split("\t")[0] in {x["key"] for x in r["ledger"]}
        }
        path = lp.AP / "candidate_cds_coverage.tsv"
        lines = [cov.get(ln.split("\t")[0], ln) for ln in unix(path)]
        path.write_text("\n".join(lines) + "\n")
        # Manifest: read fresh (other agents edit it), update only the switched rows' design
        # columns in place, and write once, checking the row count right before.
        with MANIFEST.open() as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            cols, rows = list(reader.fieldnames or []), list(reader)
        n = len(rows)
        hit = 0
        retired = 0
        for row in rows:
            if row["library"] != ld.LIBRARY_L1 or not ld.is_active(row):
                continue
            if row["key"] in r["design_rows"]:
                row.update(r["design_rows"][row["key"]])
                hit += 1
            if row["candidate_id"] in r["retire"]:
                row[ld.RETIRED_COLUMN] = r["retire"][row["candidate_id"]]
                retired += 1
        assert hit == len(r["design_rows"]), f"updated {hit} of {len(r['design_rows'])} rows"
        assert retired == len(r["retire"]), f"retired {retired} of {len(r['retire'])}"
        with MANIFEST.open() as fh:
            assert sum(1 for _ in fh) - 1 == n, "manifest changed under us; rerun --write"
        with MANIFEST.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(rows)
        print(f"wrote {LEDGER.name}, {SUMMARY.name}; {hit} manifest row(s) updated in place")
    if a.publish:
        from deliverome_analysis.s3 import processed

        for prefix in (REF_PREFIX, CAND_PREFIX):
            if processed.list_keys(prefix=prefix):
                raise SystemExit(f"processed://{prefix} exists; records are never overwritten")
        for name in (
            "reference_selection.tsv",
            "reference_disagreements.tsv",
            "reference_selection_summary.json",
        ):
            processed.upload(str(r["reference_dir"] / name), REF_PREFIX + name)
        for name, data in (
            ("candidate_cds.json", record_text(r["cds_record"])),
            ("candidate_constructs.json", record_text(r["con_record"])),
            ("switched_vectors.fasta", "".join(r["fasta"])),
        ):
            path = work / name
            path.write_text(data)
            processed.upload(str(path), CAND_PREFIX + name)
        print(f"published {REF_PREFIX} and {CAND_PREFIX}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
