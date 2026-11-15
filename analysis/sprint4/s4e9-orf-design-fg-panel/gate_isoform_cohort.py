"""Gate every isoform-cohort form before the panel draw: native start, and the canonical's
N-terminal architecture kept.

Decision: Becca Carlson, 2026-11-12 ("Gate cohort, redraw"). Implemented by the Claude Code
session l1-reviewer-package (subagent), PR #283. The gate is applied to every candidate form
of ``isoform_cohorts.tsv`` (the pinned panel snapshot) before the draw. The isoform sampler
reads the table this script writes through its validity hook (``invalid()`` in
``sample_isoform_panel.py``, applied by ``panel_draw.replace_rejected``), so a failing
original pick is replaced and a failing form is never drawn as a replacement.

A form passes only if both hold.

1. **Native ATG start and a supported CDS start.** Its protein starts with M, and an Ensembl
   115 transcript in the pinned CDS file encodes exactly that protein with a CDS that starts
   ATG. The transcript is found by exact peptide match against the pinned Ensembl 115 pep
   file, and its CDS record must translate to the same protein (no allele tolerance).
2. **The canonical's N-terminal orientation and first TM helix kept**, by DeepTMHMM v1 (the
   topology of record), against the gene's selected reference: the active L1
   ``library_reference`` member's designed protein (``library_design.designed_proteins``).
   Both proteins are looked up by sequence in the pinned topology snapshot. The form fails if
   the reference has a signal peptide and the form has none, or if the N-terminal orientation
   after any signal peptide differs. For a non-GPI reference it also fails if the reference's
   first TM helix sequence is absent from the form, or is called membrane (``M``) over fewer
   than half its residues in the form's v1 string. A form or reference with no v1 row fails as
   unverifiable; no prediction is run here.

A form's protein is its UniProt isoform sequence (UniProt 2026-10-11, the resolver's source)
when the cohort row names an accession, otherwise the Ensembl 115 translation of its first
transcript. The GPI set is ``gpi_set.tsv`` (the GPI set of record).

    uv run --locked python gate_isoform_cohort.py            # write isoform_gate.tsv + summary
    uv run --locked python gate_isoform_cohort.py --check    # fail unless the table reproduces
    uv run --locked python gate_isoform_cohort.py --publish  # upload the table's pinned key
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
AP = REPO / "scratch/becca/appris-principal"
sys.path.insert(0, str(AP))

MANIFEST = HERE / "candidate_manifest_tranches.tsv"
GPI_SET = HERE / "gpi_set.tsv"
OUT = HERE / "isoform_gate.tsv"
SUMMARY = HERE / "isoform_gate_summary.json"
#: Where the table is pinned for ``--snapshot`` replay; ``panel_inputs.json`` records its hash.
#: ``-2``: regenerated after the second-round reference fixes (DISP2, TSPAN10, LMLN, EPHA6
#: references changed length; no verdict changed). The first table stays at its own key.
GATE_KEY = (
    "sprint4/s4e9-orf-design-fg-panel/library-design/panel-inputs/"
    "2026-11-12-isoform-gate-2/isoform_gate.tsv"
)

DECIDED_BY = "Becca Carlson, 2026-11-12 (Gate cohort, redraw)"
IMPLEMENTED_BY = "Claude Code session l1-reviewer-package (subagent), PR #283"
DECISION_REASON = (
    "a panel form must start at a native ATG that an Ensembl 115 CDS supports, and keep its "
    "canonical's signal peptide, N-terminal orientation and first TM helix (DeepTMHMM v1)"
)
COLUMNS = [
    "form_id",
    "gene",
    "uniprot",
    "ensembl_transcripts",
    "cohort_tier",
    "verdict",
    "reason",
    "transcript",
    "form_length",
    "reference_design_id",
    "reference_length",
    "reference_gpi",
    "reference_sp",
    "form_sp",
    "reference_n_term",
    "form_n_term",
    "first_tm_helix",
    "first_tm_helix_membrane_fraction",
    "decided_by",
    "implemented_by",
    "decision_reason",
]
FLIP = {"I": "O", "O": "I"}


def form_id(r: dict) -> str:
    """The sampler's form id (``sample_isoform_panel.form_id``)."""
    return f"{r['gene']}:{r['uniprot'].strip() or r['ensembl_transcripts']}"


def sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def n_term(topo: str) -> tuple[int, str]:
    """(signal-peptide length, N-terminal side after it: I, O or ?)."""
    sp = len(topo) - len(topo.lstrip("S"))
    rest = topo[sp:]
    if not rest:
        return sp, "?"
    if rest[0] in "IO":
        return sp, rest[0]
    # Starts in a helix: the N-terminus is on the side opposite the first loop after it.
    after = next((c for c in rest.lstrip(rest[0]) if c in "IO"), None)
    return sp, FLIP.get(after, "?")


def first_helix(topo: str) -> tuple[int, int] | None:
    start = topo.find("M")
    if start < 0:
        return None
    end = start
    while end < len(topo) and topo[end] == "M":
        end += 1
    return start, end


def load_inputs(work: Path) -> dict:
    import candidate_inputs
    import resolve_candidate_cds as rcd
    from gpi import gpi_genes
    from panel_inputs import stage

    from deliverome_analysis import library_design as ld
    from deliverome_analysis import topology
    from deliverome_analysis.s3 import external

    snap = stage(work / "panel")
    with (snap / "isoform_cohorts.tsv").open() as fh:
        cohort = list(csv.DictReader(fh, delimiter="\t"))
    pep = rcd.read_fasta(candidate_inputs.pep_fasta("human"))
    cds = rcd.read_fasta(candidate_inputs.cds_fasta("human"))
    by_protein: dict[str, list[str]] = defaultdict(list)
    pep_of: dict[str, str] = {}
    for _ensp, (head, seq) in pep.items():
        tx = next((f.split(":", 1)[1] for f in head.split() if f.startswith("transcript:")), "")
        tx = tx.split(".")[0]
        if tx:
            prot = seq.rstrip("*").upper()
            by_protein[prot].append(tx)
            pep_of[tx] = prot
    # Keyed on the design id, which renumbering never changes.
    did = ld.design_ids(MANIFEST)
    refs = {
        m["gene_symbol"]: (did[cid], m["sequence"])
        for cid, m in ld.designed_proteins(MANIFEST, library=ld.LIBRARY_L1).items()
        if m["candidate_set"] == "library_reference"
    }
    snap_path = Path(external.download(topology.SNAPSHOT_KEY, work / "topology.parquet"))
    if sha256(snap_path) != topology.SNAPSHOT_SHA256:
        raise SystemExit("topology snapshot differs from topology.SNAPSHOT_SHA256")
    import pandas as pd

    return {
        "cohort": cohort,
        "uniprot": rcd.load_uniprot("human"),
        "pep_of": pep_of,
        "by_protein": by_protein,
        "cds": cds,
        "translate": rcd.translate,
        "refs": refs,
        "gpi": gpi_genes(GPI_SET),
        "snapshot": pd.read_parquet(snap_path),
    }


def gate(inp: dict) -> list[dict]:
    from deliverome_analysis import topology

    cohort, refs = inp["cohort"], inp["refs"]
    proteins: dict[str, str] = {}
    for r in cohort:
        acc = r["uniprot"].strip()
        txs = [t.split(".")[0] for t in r["ensembl_transcripts"].split(";") if t.strip()]
        seq = inp["uniprot"].get(acc) if acc else None
        if seq is None:
            seq = next((inp["pep_of"][t] for t in txs if t in inp["pep_of"]), "")
        proteins[form_id(r)] = seq
    seqs = {f"form:{k}": v for k, v in proteins.items() if v}
    seqs |= {f"ref:{g}": s for g, (_, s) in refs.items()}
    topo = topology.resolve(seqs, tool=topology.TOOL_V1, snapshot=inp["snapshot"])
    v1 = {
        row["key"]: row["per_residue_topology"] for row in topo.to_dict("records") if row["found"]
    }

    out = []
    for r in cohort:
        fid = form_id(r)
        prot = proteins[fid]
        reasons: list[str] = []
        row = {c: "" for c in COLUMNS} | {
            "form_id": fid,
            "gene": r["gene"],
            "uniprot": r["uniprot"],
            "ensembl_transcripts": r["ensembl_transcripts"],
            "cohort_tier": r["cohort_tier"],
            "form_length": len(prot),
            "decided_by": DECIDED_BY,
            "implemented_by": IMPLEMENTED_BY,
            "decision_reason": DECISION_REASON,
        }
        # ---- 1. native ATG start, supported by an Ensembl 115 CDS -------------------------
        if not prot:
            reasons.append("no protein sequence for the form")
        else:
            if not prot.startswith("M"):
                reasons.append("protein does not start with M")
            named = [t.split(".")[0] for t in r["ensembl_transcripts"].split(";") if t.strip()]
            hits = inp["by_protein"].get(prot, [])
            hits = sorted(hits, key=lambda t: (t not in named, t))
            good, non_atg = [], []
            for t in hits:
                rec = inp["cds"].get(t)
                if rec is None:
                    continue
                dna = rec[1].upper()
                if inp["translate"](dna).rstrip("*") != prot:
                    continue
                (good if dna.startswith("ATG") else non_atg).append(t)
            if good:
                row["transcript"] = good[0]
            elif non_atg:
                row["transcript"] = non_atg[0]
                reasons.append(f"Ensembl 115 CDS starts {inp['cds'][non_atg[0]][1][:3]}, not ATG")
            elif hits:
                row["transcript"] = hits[0]
                reasons.append("no Ensembl 115 CDS record translates exactly to the protein")
            else:
                reasons.append("no Ensembl 115 transcript encodes exactly this protein")
        # ---- 2. the reference's N-terminal architecture ----------------------------------
        ref = refs.get(r["gene"])
        if ref is None:
            reasons.append("unverifiable: gene has no active L1 library_reference member")
        else:
            row["reference_design_id"] = ref[0]
            row["reference_length"] = len(ref[1])
            gpi = r["gene"] in inp["gpi"]
            row["reference_gpi"] = str(gpi)
            rt, ft = v1.get(f"ref:{r['gene']}"), v1.get(f"form:{fid}")
            if rt is None:
                reasons.append("unverifiable: reference has no DeepTMHMM v1 row")
            if ft is None:
                reasons.append("unverifiable: form has no DeepTMHMM v1 row")
            if rt is not None and ft is not None:
                rsp, rn = n_term(rt)
                fsp, fn = n_term(ft)
                row |= {
                    "reference_sp": rsp,
                    "form_sp": fsp,
                    "reference_n_term": rn,
                    "form_n_term": fn,
                }
                if rsp and not fsp:
                    reasons.append("loses the reference's signal peptide")
                if rn != fn:
                    reasons.append(f"N-terminus {fn} after any signal peptide, reference {rn}")
                helix = None if gpi else first_helix(rt)
                if helix is not None:
                    hseq = ref[1][helix[0] : helix[1]]
                    row["first_tm_helix"] = f"{helix[0] + 1}-{helix[1]}"
                    at = prot.find(hseq)
                    if at < 0:
                        reasons.append("reference's first TM helix sequence absent")
                    else:
                        frac = ft[at : at + len(hseq)].count("M") / len(hseq)
                        row["first_tm_helix_membrane_fraction"] = f"{frac:.2f}"
                        if frac < 0.5:
                            reasons.append("reference's first TM helix not called membrane")
        row["verdict"] = "fail" if reasons else "pass"
        row["reason"] = "; ".join(reasons)
        out.append(row)
    return out


def tsv_text(rows: list[dict]) -> str:
    lines = ["\t".join(COLUMNS)] + ["\t".join(str(r[c]) for c in COLUMNS) for r in rows]
    return "\n".join(lines) + "\n"


def summarise(rows: list[dict]) -> dict:
    fail = [r for r in rows if r["verdict"] == "fail"]
    reasons = Counter(
        x.split(":")[0] if x.startswith("unverifiable") else x
        for r in fail
        for x in r["reason"].split("; ")
    )
    first = Counter(r["reason"].split("; ")[0] for r in fail)
    return {
        "decision": DECIDED_BY,
        "implemented_by": IMPLEMENTED_BY,
        "forms": len(rows),
        "pass": len(rows) - len(fail),
        "fail": len(fail),
        "by_tier": {
            t: dict(Counter(r["verdict"] for r in rows if r["cohort_tier"] == t))
            for t in sorted({r["cohort_tier"] for r in rows})
        },
        "fail_reasons_any": dict(sorted(reasons.items(), key=lambda kv: -kv[1])),
        "fail_reason_first": dict(sorted(first.items(), key=lambda kv: -kv[1])),
        "unverifiable": sum("unverifiable" in r["reason"] for r in fail),
        "unverifiable_only": sum(
            all(x.startswith("unverifiable") for x in r["reason"].split("; ")) for r in fail
        ),
        "topology_tool": "deeptmhmm_v1",
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--publish", action="store_true")
    args = ap.parse_args(argv)
    work = Path(tempfile.mkdtemp(prefix="isoform-gate-"))
    rows = gate(load_inputs(work))
    text = tsv_text(rows)
    summary = json.dumps(summarise(rows), indent=2) + "\n"
    if args.check:
        assert OUT.read_text() == text, f"{OUT.name} does not reproduce"
        assert SUMMARY.read_text() == summary, f"{SUMMARY.name} does not reproduce"
        print(f"{OUT.name} reproduces ({len(rows)} forms)")
        return 0
    if args.publish:
        from deliverome_analysis.s3 import processed

        if OUT.read_text() != text:
            raise SystemExit(f"{OUT.name} differs from the rerun; write and commit first")
        if processed.exists(GATE_KEY):
            raise SystemExit(f"processed://{GATE_KEY} exists; pinned keys are never overwritten")
        processed.upload(str(OUT), GATE_KEY)
        print(f"published processed://{GATE_KEY} sha256 {sha256(OUT)}")
        return 0
    OUT.write_text(text)
    SUMMARY.write_text(summary)
    print(summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
