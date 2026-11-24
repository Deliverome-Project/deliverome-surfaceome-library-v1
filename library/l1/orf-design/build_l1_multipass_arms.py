"""Build the tagged inserts for the L1 "N at residue 1, multipass" rung.

Every member of the rung (`l1_tag_paths.tsv`) gets one row and one design, from the plan in
`nterm_res1_multipass.tsv` (see NTERM_RES1_MULTIPASS.md):

* ``igk_tag`` -- Igk leader + Asp + the member's Met1 + terminal cassette + native 2..end.
* ``tag_alone`` -- Met1 + terminal cassette + native 2..end.
* ``moved`` -- not built here; the member belongs in another lane, recorded in the row.
* ``replace`` -- not built; the member itself should be replaced, recorded in the row.

Every row records ``decided_by``, ``implemented_by``, ``reason`` and ``rejected_alternative``.

| Decision | Decided by |
|---|---|
| Igk + tag by default for the rung's GPCRs | Becca, 2026-11-10 (on Aashish Manglik's suggestion) |
| Full 51-aa cassette at residue 1, leading GGGGS kept | Becca, 2026-11-11 ("keep full"); l1-nterm-type3, PR #282 |
| Isoform gate and redraw; ADGRF3 to its signal-peptide form | Becca, 2026-11-12; applied by l1-reviewer-package, PR #283 |
| One arm per gene for the expression-rescue genes: the 23 references on this rung built Igk + tag, their rescue twins retired (`rescue_one_arm.tsv`) | Becca, 2026-11-13; per-gene choice by this session, delegated |
| Per-member calls: Igk for the non-GPCRs; leader-free for CRHR2; lane moves for MFSD2A x2, AQP12B and GLP2R (since applied by the writer) replacing cyno NIPAL4 and the OPN3 isoform (since restarted and dropped); recoding the OR8G3 start to ATG (since applied) | this session, delegated by Becca 2026-11-12 ("use your best judgment") |

**Layout**, identical to `build_l1_type3_arms.py` (the L1 Esp3I wrapper)::

    CGTCTCA | Kozak | [Igk leader (60 nt) | GAC] | ATG | cassette (51 aa, p003 tag codons, own linkers) |
            ORF codons 2..end incl. native stop | CC | T7 | 5' handle | barcode (20 N) | Feldman 3' |
            10x CS1 | AGA GAGACG

For each built row the build checks: the exact encoded protein; a single stop, at the end; and
only the two wrapper Esp3I sites in the insert (no other BbsI/Esp3I/BsaI), both alone and in
the p0039 backbone. `multipass_arms_mature.fasta` holds each built arm's mature protein
(the leader removed) for refolding.

    uv run --locked python build_l1_multipass_arms.py
    uv run --locked python build_l1_multipass_arms.py --check
"""

from __future__ import annotations

import argparse
import csv
import functools
import hashlib
import json
import sys
from pathlib import Path

from build_constructs import CASSETTE_AA
from diversify_arm_cassettes import arm_cassette

from deliverome_analysis import library_design as ld
from deliverome_analysis import plasmid_qc as qc
from deliverome_analysis import seqtools, topology

LEADER = qc.IGK_LEADER_NT
LEADER_AA = qc.IGK_LEADER_AA
PLUS1 = qc.IGK_PLUS1_NT
PLUS1_AA = qc.IGK_PLUS1_AA

HERE = Path(__file__).parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO / "library/l1/proteoform-selection"))

MANIFEST = HERE / "candidate_manifest_tranches.tsv"
TAG_PATHS = HERE / "l1_tag_paths.tsv"
PLAN = HERE / "nterm_res1_multipass.tsv"
#: Applied start trims (#283): carried into the arm rows of the members they changed.
START_TRIMS = HERE / "ortholog_start_trims.tsv"
OUT_TSV = HERE / "constructs_l1_multipass.tsv"
OUT_FASTA = HERE / "constructs_l1_multipass.fasta"
OUT_MATURE = HERE / "multipass_arms_mature.fasta"
#: The pinned constructs record, the same one `library_design` resolves CDS from.
CONSTRUCTS_KEY = ld.CONSTRUCTS_KEY
RUNG = "N at residue 1, multipass"
BARCODE_SLOT = qc.BARCODE_SLOT
ENZYMES = list(qc.DOMESTICATION_ENZYMES)

IMPLEMENTED_BY = "Claude Code session l1-nterm-res1-multipass (PR #280)"
GPCR_DECIDED = "Becca, 2026-11-10 (Igk by default for GPCRs, on Aashish Manglik's suggestion)"
CASSETTE_DECIDED = 'Becca, 2026-11-11 ("keep full"); implemented by l1-nterm-type3, PR #282'
REDRAW_DECIDED = (
    'Becca, 2026-11-12 ("Gate cohort, redraw"); applied by l1-reviewer-package, PR #283'
)

#: Becca delegated protein-level calls on 2026-11-12: "try to not queue protein specific things
#: for management decisions ... use your best judgment to complete a new and improved evidence
#: based l1 completely by the morning." Every call below is this session's, made on the
#: evidence cited, and records the rejected alternative so it can be overturned.
DELEGATED = 'l1-nterm-res1-multipass (delegated by Becca 2026-11-12, "use your best judgment")'

#: design_id -> (lane, reason, rejected alternative) for members moved off this rung. Sent
#: to l1-reviewer-package as `apply`. Not built here: their arms belong to the new lane.
MOVED = {
    "l1-p1213": ("internal (TM5-TM6 extracellular domain)", "MFSD2A N-terminus is cytoplasmic: two cryo-EM structures (PMIDs 34349262, 34135507) place the glycosylated ECD in the TM5-TM6 loop and the post-TM12 helix inside; a residue-1 tag would be cytoplasmic and a leader would invert the transporter", "keep on this rung with Igk"),
    "l1-p2466": ("internal (TM5-TM6 extracellular domain)", "MFSD2A isoform with the same N-terminus as l1-p1213", "keep on this rung with Igk"),
    "l1-p0128": ("internal", "AQP12B: v1's extra N-terminal helix conflicts with the aquaporin fold (N- and C-termini cytoplasmic); AQP12 is intracellular when expressed (PMID 19726746)", "keep on this rung with Igk"),
    "l1-p0724": ("N after SP", "GLP2R: class B GPCR; v2 calls a signal peptide (cleavage near G64/S66, lane check medium confidence); a residue-1 tag would sit on a cleaved segment; Tedman's HA-after-Met1 construct read 121 a.u.", "keep on this rung with Igk"),
}  # fmt: skip

#: design_id -> (reason, rejected) for members not built because the member itself should be
#: replaced. Recorded here; the replacement belongs to the manifest's writer. Each decision is
#: tied to the sequence it was made on (REPLACED_SHA): once the writer swaps the sequence, the
#: new form is built normally and this row stays as the trail.
REPLACE = {
    "l1-p2962": ("cyno NIPAL4 ortholog carries a 214-residue Arg/Lys/Gly-rich extension before Met215, which aligns to human NIPAL4 Met1; DeepTMHMM v2 calls that N-terminus cytoplasmic, so a leader would sit in front of a likely cytoplasmic segment; replace with a proteoform starting at Met215 (lane check, PR #296)", "build with Igk + tag"),
    "l1-p2477": ("OPN3 isoform is a 189-aa fragment with helices 1-4 only (no helix 7, no retinal-binding Lys); it passes the isoform gate but cannot fold as a GPCR, and the canonical OPN3 (l1-p1353) is already on this rung; drop or replace (lane check, PR #296)", "build with Igk + tag as the gate allows"),
}  # fmt: skip

#: Sequence each REPLACE decision was made on (sha256 of the designed protein).
REPLACED_SHA = {
    "l1-p2962": "141a857c93a6d082ab04b1ab6e5be8e8c0293ad6b3e16cc1d89f3ab095282538",
    "l1-p2477": "",  # retired by the proteoform audit before a sequence swap
}

#: SLC61A1 stays on this rung (decided by the manifest's writer, 2026-11-13). This session
#: had moved it to internal on the major facilitator fold argument; that view is kept as the
#: rejected alternative.
SLC61A1 = {"l1-p1994", "l1-p3105"}
SLC61A1_DECIDED = (
    'l1-reviewer-package (delegated by Becca 2026-11-12, "use your best judgment"), 2026-11-13'
)
SLC61A1_REASON = (
    "stays on N at residue 1: DeepTMHMM v1 and v2, SignalP 6 and the AlphaFold model (13 "
    "crossings, pLDDT 89) all put the N-terminus outside, and the lane check (PR #296) rates "
    "the lane as holding; Igk + tag as for the rung's other non-GPCRs"
)

#: TSPO2, the lane check's one unresolved member (PR #296), decided 2026-11-13.
TSPO2 = "l1-p2302"
TSPO2_DECIDED = (
    'l1-nterm-res1-multipass (delegated by Becca 2026-11-12, "use your best judgment"), 2026-11-13'
)
TSPO2_REASON = (
    "kept on N at residue 1 with Igk: with 5 helices exactly one terminus is extracellular. "
    "DeepTMHMM v1 (the topology of record) and TSPO1 homology (N-terminus non-cytosolic; PMIDs "
    "9442055, 24653034) put the N-terminus outside, the 4-residue tail MRLQ carries the weak "
    "positive-inside bias for N-out (delta +1), and the Igk arm passes the tagged refold. Only "
    "DeepTMHMM v2 calls N-in. Low confidence: an N-glycosylation reporter or surface staining "
    "of an N- vs C-tagged pair would settle it"
)
TSPO2_REJECTED = "move to the C-terminal lane (right only if v2's N-in call is correct)"

#: design_id -> (codon, reason, rejected) for members whose start codon is recoded to ATG here.
RECODE_START = {
    "l1-p1525": ("GTG", "OR8G3 opens on GTG, which initiates poorly in a construct; recoded to ATG with no added residue, as decided for SLC22A17 and TRPV6", "keep in curation, unbuilt"),
}  # fmt: skip

CRHR2_REASON = (
    "CRHR2 is built leader-free: its N-terminal pseudo signal peptide is uncleaved and cannot "
    "target the ER on its own (Rutz 2006, PMID 16766521). Behind a leader it would become an "
    "untested internal hydrophobic segment. Tedman displayed this transcript leader-free at "
    "2,385 a.u. (screen median about 1,000)"
)
#: One arm per gene for the 25 expression-rescue genes. Becca, 2026-11-13: "we are pooling the
#: whole L1 library together. why resuce pairs no i think we should make a single choice for
#: each L1". The per-gene choice is this session's, under her 2026-11-12 delegation, recorded
#: gene by gene in `rescue_one_arm.tsv` (`rescue_one_arm.py`). For the 23 genes on this rung the
#: library reference is built Igk + tag, the GPCR default, which is the same protein as its
#: rescue twin; the twin is retired.
ONE_ARM_DECIDED = (
    'Becca, 2026-11-13 ("make a single choice for each L1"); arm chosen by '
    'l1-nterm-res1-multipass (delegated by Becca 2026-11-12, "use your best judgment")'
)
ONE_ARM_REASON = (
    "one arm per gene: the expression-rescue twin is the same Igk + tag protein, so it is "
    "retired and this reference carries the gene (rescue_one_arm.tsv)"
)


def read_table(path: Path) -> list[dict]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


@functools.cache
def _added_2026_11_21_22() -> set[str]:
    """Members added to L1 on 2026-11-21/22: dated gene additions (gene_list_additions.tsv) and
    the 2026-11-22 isoform/ortholog redraw (panel_redraw_2026-11-22_changes.tsv) and its
    refills (panel_refill_2026-11-22_changes.tsv)."""
    out: set[str] = set()
    genes = {
        r["gene"]
        for r in read_table(HERE / "gene_list_additions.tsv")
        if r["added"] >= "2026-11-21"
    }
    for r in read_table(MANIFEST):
        if r["key"].startswith("lib:") and r["key"].split(":", 1)[1] in genes:
            out.add(r["design_id"])
    for name in ("panel_redraw_2026-11-22_changes.tsv", "panel_refill_2026-11-22_changes.tsv"):
        redraw = HERE / name
        if redraw.exists():
            out |= {r["candidate_id"] for r in read_table(redraw) if r["change"].startswith("add")}
    return out


DELEGATED_2026_11_21 = (
    'Claude Code session l1-review-export-scope, under Becca\'s 2026-11-21 "Yes, rebuild now"'
)


def design_for(row: dict, rescue_genes: set[str]) -> tuple[str, str, str, str, str]:
    """(design, status, decided_by, reason, rejected_alternative) for one member."""
    d = row["design_id"]
    if d in MOVED:
        lane, reason, rejected = MOVED[d]
        return "moved", "decided", DELEGATED, f"moved to {lane}: {reason}", rejected
    if d in REPLACE and row["sequence_sha256"] == REPLACED_SHA[d]:
        reason, rejected = REPLACE[d]
        return "replace", "decided", DELEGATED, reason, rejected
    if d in SLC61A1:
        return (
            "igk_tag",
            "decided",
            SLC61A1_DECIDED,
            SLC61A1_REASON,
            "move to internal (MFS N-in fold argument and a 13-residue TM1; this session's view)",
        )
    if row["gpcr"] == "1":
        if row["gene_symbol"] == "CRHR2":
            return "tag_alone", "decided", DELEGATED, CRHR2_REASON, "Igk + tag (GPCR default)"
        if row["candidate_set"] == "library_reference" and row["gene_symbol"] in rescue_genes:
            return (
                "igk_tag",
                "decided",
                f"{GPCR_DECIDED}; {ONE_ARM_DECIDED}",
                ONE_ARM_REASON,
                "tag alone, keeping the rescue twin as a within-gene leader contrast (no pairs "
                "in the pooled L1, Becca 2026-11-13)",
            )
        reason = "GPCR on the multipass residue-1 rung"
        decided_by = GPCR_DECIDED
        if d in RECODE_START:
            reason += f"; {RECODE_START[d][1]}"
            decided_by += f"; start recode by {DELEGATED}"
        return "igk_tag", "decided", decided_by, reason, "tag alone"
    verdict = row["planned_leader"]
    if d == TSPO2:
        return "igk_tag", "decided", TSPO2_DECIDED, TSPO2_REASON, TSPO2_REJECTED
    if verdict == "igk":
        why = f"literature review verdict igk ({row['planned_leader_reason']})"
        return "igk_tag", "decided", DELEGATED, why, "tag alone"
    if not verdict and d in _added_2026_11_21_22():
        why = (
            "added 2026-11-21/22 (gene addition or panel redraw); not separately reviewed, so built as the "
            "rung's non-GPCR default: Igk + tag, which removes the leaderless translocation and "
            "TM1-flip risks at no identified cost (DeepTMHMM v1 places its N-terminus outside)"
        )
        return "igk_tag", "decided", DELEGATED_2026_11_21, why, "tag alone"
    if not verdict:
        why = (
            "entered with the #283 redraw, not separately reviewed; built as its family's "
            "reviewed members (all igk: N-out established, or the leader removes the leaderless "
            "translocation and TM1-flip risks at no identified cost)"
        )
        return "igk_tag", "decided", DELEGATED, why, "tag alone"
    raise ValueError(f"{d}: unhandled non-GPCR verdict {verdict!r}")


def build() -> list[dict]:
    from deliverome_analysis.s3 import processed

    paths = {r["candidate_id"]: r for r in read_table(TAG_PATHS)}
    plan = read_table(PLAN)
    rung = sorted(c for c, r in paths.items() if r["placement"] == RUNG)
    assert sorted(r["candidate_id"] for r in plan) == rung, "plan table is not the current rung"
    manifest_rows = read_table(MANIFEST)
    manifest = {r["candidate_id"]: r for r in manifest_rows if ld.is_active(r)}
    rescue_genes = {
        r["key"].split(":")[1] for r in manifest_rows if r["set"] == "expression_rescue"
    }
    proteins = ld.designed_proteins(MANIFEST, library=ld.LIBRARY_L1)
    constructs = json.loads(Path(processed.download(CONSTRUCTS_KEY)).read_text())
    assert seqtools.translate(qc.P003_CASSETTE_NT) == CASSETTE_AA and len(CASSETTE_AA) == 51
    assert seqtools.translate(LEADER) == LEADER_AA and seqtools.translate(PLUS1) == PLUS1_AA
    # A decision whose member has already left the rung (the writer applied the move, or the
    # member was retired) is reported, not an error: the record in this file stays the trail.
    applied = sorted(
        (set(MOVED) | set(REPLACE) | set(RECODE_START)) - {r["design_id"] for r in plan}
    )
    if applied:
        print(f"decisions already applied off this rung: {applied}")

    # Sequence changes another session applied to a member (e.g. an ortholog restarted at its
    # native Met). Matched on the restarted protein's hash, so a row only claims a change
    # that is actually in the sequence being built.
    trims = {}
    if START_TRIMS.exists():
        for t in read_table(START_TRIMS):
            if t["status"] == "applied":
                trims[(t["design_id"], t["restarted_protein_sha256"])] = t

    rows = []
    for p in sorted(plan, key=lambda r: (r["design_id"], r["candidate_id"])):
        cid = p["candidate_id"]
        native = proteins[cid]["sequence"].rstrip("*")
        assert paths[cid]["sequence_sha256"] == topology.sequence_sha256(native), f"{cid} moved"
        for design, status, decided_by, reason, rejected in [design_for(p, rescue_genes)]:
            row = {
                "construct_id": f"{p['design_id']}:{design}",
                "design_id": p["design_id"],
                "candidate_id": cid,
                "gene_symbol": p["gene_symbol"],
                "candidate_set": p["candidate_set"],
                "species": p["species"],
                "gpcr": p["gpcr"],
                "design": design,
                "status": status,
                "native_aa": len(native),
                "protein_aa": "",
                "mature_aa": "",
                "coding_nt": "",
                "insert_nt": "",
                "cassette": "",
                "type_iis_sites": "",
                "decided_by": decided_by,
                "cassette_decided_by": "",
                "implemented_by": IMPLEMENTED_BY,
                "reason": reason,
                "rejected_alternative": rejected,
                "sequence_change": "",
                "sequence_decided_by": "",
                "sequence_implemented_by": "",
                "protein_sha256": "",
                "mature_sha256": "",
                "native_sha256": paths[cid]["sequence_sha256"],
                "insert": "",
                "mature": "",
            }
            trim = trims.get((p["design_id"], paths[cid]["sequence_sha256"]))
            if trim:
                row.update(
                    sequence_change=f"restarted at native Met{trim['restart_codon']} ({trim['reason']})",
                    sequence_decided_by=trim["decided_by"],
                    sequence_implemented_by=trim["implemented_by"],
                )
            if design not in {"moved", "replace"}:
                orf = constructs[manifest[cid]["key"]]["construct"]
                assert seqtools.translate(orf) == native + "*", (
                    f"{cid}: ORF does not encode protein"
                )
                # Recode here only while the CDS record still carries the native start; once
                # the writer applies the ATG fix to the record, the ORF already starts ATG.
                if p["design_id"] in RECODE_START and orf[:3] != "ATG":
                    assert orf[:3] == RECODE_START[p["design_id"]][0], f"{cid}: unexpected start"
                    orf = "ATG" + orf[3:]
                    native = "M" + native[1:]
                assert orf[:3] == "ATG", f"{cid}: ORF does not start with ATG"
                head = (LEADER + PLUS1) if design == "igk_tag" else ""
                coding = head + orf[:3] + arm_cassette(f"{p['design_id']}:{design}") + orf[3:]
                mature = native[0] + CASSETTE_AA + native[1:]
                if design == "igk_tag":
                    mature = PLUS1_AA + mature
                protein = (LEADER_AA if design == "igk_tag" else "") + mature
                got = seqtools.translate(coding)
                assert got == protein + "*" and got.count("*") == 1, f"{cid}: protein mismatch"
                insert = qc.wrap_insert(coding)
                assert qc.wrapper_sites_only(insert), (
                    f"{cid}: {seqtools.find_sites(insert, ENZYMES)}"
                )
                row.update(
                    protein_aa=len(protein),
                    mature_aa=len(mature),
                    coding_nt=len(coding),
                    insert_nt=len(insert),
                    cassette="full 51 aa; p003 tag codons, own linker encoding",
                    type_iis_sites="Esp3I x2 (wrapper)",
                    cassette_decided_by=CASSETTE_DECIDED,
                    protein_sha256=hashlib.sha256(protein.encode()).hexdigest(),
                    mature_sha256=hashlib.sha256(mature.encode()).hexdigest(),
                    insert=insert,
                    mature=mature,
                )
            rows.append(row)
    vector_check([r for r in rows if r["insert"]])
    return rows


def vector_check(rows: list[dict]) -> None:
    """Every insert in the p0039 backbone: backbone sites plus the two wrapper Esp3I, no more."""
    import candidate_inputs
    from Bio import SeqIO

    bb = str(SeqIO.read(candidate_inputs.backbone(), "genbank").seq).upper()
    qc.check_inserts_in_backbone([(r["construct_id"], r["insert"]) for r in rows], bb)


def render(rows: list[dict]) -> tuple[str, str, str]:
    cols = [c for c in rows[0] if c not in {"insert", "mature"}]
    tsv = "\t".join(cols) + "\n"
    tsv += "".join("\t".join(str(r[c]) for c in cols) + "\n" for r in rows)
    built = [r for r in rows if r["insert"]]
    fasta = "".join(
        f">{r['construct_id']} {r['candidate_id']} {r['gene_symbol']} {r['species']} {r['status']}\n"
        f"{r['insert']}\n"
        for r in built
    )
    mature = "".join(
        f">{r['construct_id']} {r['candidate_id']} {r['gene_symbol']} native_sha256={r['native_sha256']}\n"
        f"{r['mature']}\n"
        for r in built
    )
    return tsv, fasta, mature


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--check", action="store_true", help="rebuild and compare with the record")
    args = ap.parse_args(argv)
    rows = build()
    tsv, fasta, mature = render(rows)
    outs = {OUT_TSV: tsv, OUT_FASTA: fasta, OUT_MATURE: mature}
    if args.check:
        same = all(p.read_text() == text for p, text in outs.items())
        print("record matches rebuild" if same else "record DIFFERS from rebuild")
        return 0 if same else 1
    for p, text in outs.items():
        p.write_text(text)
    counts: dict[tuple[str, str], int] = {}
    for r in rows:
        counts[(r["design"], r["status"])] = counts.get((r["design"], r["status"]), 0) + 1
    members = len({r["candidate_id"] for r in rows})
    print(f"wrote {OUT_TSV.name}, {OUT_FASTA.name}, {OUT_MATURE.name}: {members} members")
    for k, n in sorted(counts.items()):
        print(f"  {k[0]:10s} {k[1]:9s} {n}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
