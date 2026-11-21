"""Which surface-proteomics gap genes could be added to L1 under Becca's 2026-11-20 criteria.

The challenge set is every protein that passes the cross-method survey's rule B and is in
neither L1 nor L2 (``surface_proteomics_logodds_compare.py`` explains rule B). Becca, 2026-11-20:
add such a gene "as long as they have below the largest ectodomain limit and no co receptor
issues and reasonable anchoring". This script turns that into four tests and reports, per gene,
the values each test reads. It changes nothing: the library, ``gene_list_additions.tsv``,
constructs, barcodes and the order are untouched, because adding a CDS changes the order of
record, which is Becca's call.

======================= ===================================================================
test                    how it is decided
======================= ===================================================================
size                    ORF (canonical length + stop, x3) <= ``final_library.MAX_ORF_NT``
                        (5,000 nt), the library's synthesis/packaging cap. The filter has no
                        upper ectodomain limit; ``ecto_vs_l1_max`` compares each gene with
                        the largest ectodomain any L1 gene carries, for reference only.
co-receptor             ``final_library.coreceptor_ok``, i.e. the library's own gate with the
                        reviewed override layer (coreceptor_calls.tsv, coreceptor_rationales
                        .json, the INCLUDE_ALONE allowlist). ``not evaluable`` when the gene has
                        no catalog deep-dive record, because the partner fields live there.
anchored                DeepTMHMM v1 (``deliverome_analysis.topology``, matched on the exact
                        canonical sequence) calls >= 1 TM helix, or the catalog architecture
                        annotation calls it GPI-anchored.
ectodomain              >= ``final_library.MIN_ECD_RES`` (30) extracellular residues in the v1
                        string; GPI-anchored genes are exempt, as in the library filter.
======================= ===================================================================

The surface-evidence, accessibility, surface-fraction, expression and localization gates are
**not** applied (Becca's criteria waive them; the survey is the surface evidence), but each
gene's outcome at them is reported so the waiver is visible.

Inputs: the challenge set is the rule-B not-in-library list of the snapshot that
``surface_proteomics_inputs.json`` names as ``current``, computed from
``surface_proteomics_counts.tsv`` (written from the pinned S3 keys by
``surface_proteomics_snapshot.py``; the ``snapshots`` column says which other pinned snapshots
also list the gene); the S1E3 filter run
offline over surfaceome release 1.3.0 (``library_pipeline.run_gene_filter``, which checks every
input's SHA-256 and reproduces the 2,370 filter passes), and UniProt canonical sequences from
the key pinned as ``canonical_sequences`` in ``surface_proteomics_inputs.json``
(UniProt 2026_03 reviewed human, SHA-256 checked). Where the catalog carries a canonical topology, its length is compared with the
sequence's and a mismatch is reported in ``canonical_length_check``. Needs S3 read access.

    uv run --locked python analysis/sprint4/s4e9-orf-design-fg-panel/surface_proteomics_gap_candidates.py
"""

from __future__ import annotations

import csv
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
INPUTS = HERE / "surface_proteomics_inputs.json"
OUT = HERE / "surface_proteomics_gap_candidates.tsv"

sys.path.insert(0, str(HERE))
import library_pipeline as L  # noqa: E402

from deliverome_analysis import topology  # noqa: E402


def gap_accessions() -> dict[str, dict]:
    """Rule-B gap of the current snapshot, with the other pinned snapshots that list each gene."""
    import json

    import pandas as pd
    from surface_proteomics_logodds_compare import COUNTS, score

    current = json.loads(INPUTS.read_text())["current"]
    df = score(pd.read_csv(COUNTS, sep="	", keep_default_na=False))
    gap = df[df.pass_membrane & ~df.in_library]
    out: dict[str, dict] = {}
    for r in gap[gap.snapshot.eq(current)].itertuples():
        out[r.uniprot_accession] = {
            "sweep_symbol": r.gene_symbol,
            "snapshots": list(gap.snapshot[gap.uniprot_accession.eq(r.uniprot_accession)]),
            "n_surface_pubs": r.n_surface_pubs,
            "lod_membrane": r.lod_membrane,
        }
    return out


def canonical_sequences(dest: Path) -> dict[str, str]:
    """UniProt 2026_03 canonical sequences (isoform entries carry a ``-N`` suffix; skipped)."""
    import hashlib
    import json

    from deliverome_analysis import s3

    pin = json.loads(INPUTS.read_text())["canonical_sequences"]
    s3.external.download(pin["key"].split("://", 1)[1], dest / "human.fasta")
    got = hashlib.sha256((dest / "human.fasta").read_bytes()).hexdigest()
    assert got == pin["sha256"], f"{pin['key']}: SHA-256 is not the pinned one"
    seqs, acc = {}, None
    for line in (dest / "human.fasta").read_text().splitlines():
        if line.startswith(">"):
            acc = line[1:].split()[0].split("|")[1]
        elif acc and "-" not in acc:
            seqs[acc] = seqs.get(acc, "") + line.strip()
    return seqs


def segments(topo: str, ch: str) -> int:
    n, prev = 0, ""
    for c in topo:
        if c == ch and prev != ch:
            n += 1
        prev = c
    return n


def main() -> None:
    gap = gap_accessions()
    work = Path(tempfile.mkdtemp(prefix="gap_candidates_"))
    passed = L.run_gene_filter(work)  # stages the S1E3 inputs and runs the filter offline
    assert len(passed) == 2370, f"S1E3 filter gave {len(passed)} passes, not 2,370"
    import final_library as F  # the copy run_gene_filter staged and ran

    by_acc = {r.get("uniprot"): r for r in F.recs}
    lib = set(passed)
    l1_max_ecto = max(F.ecd_residues(r) for r in F.recs if r["gene"] in lib)
    (work / "seq").mkdir()
    seqs = canonical_sequences(work / "seq")
    topo = topology.resolve(
        {a: seqs[a] for a in gap if a in seqs}, tool=topology.TOOL_V1
    ).set_index("key")

    rows = []
    for acc, g in sorted(gap.items(), key=lambda kv: kv[1]["sweep_symbol"]):
        rec = by_acc.get(acc)
        seq = seqs.get(acc)
        t = topo.loc[acc] if acc in topo.index else None
        v1 = t["per_residue_topology"] if t is not None and t["found"] else None
        gpi = bool(rec and rec["filters"].get("subcategory") == "GPI_anchored")
        tm = segments(v1, "M") if v1 else None
        ecto = v1.count("O") if v1 else None
        orf = (len(seq) + 1) * 3 if seq else None
        gates = {name: ("pass" if ok(rec) else "fail") for name, ok in F.GATES} if rec else {}

        size_ok = orf is not None and orf <= F.MAX_ORF_NT
        cat_len = len((rec or {}).get("canonical_topology", {}).get("per_residue_topology") or "")
        if not cat_len or not seq:
            length_check = "no catalog topology" if seq else "no sequence"
        else:
            length_check = "same length" if cat_len == len(seq) else f"catalog {cat_len} aa"
        if rec is None:
            co = "not evaluable"
        else:
            co = "pass" if F.coreceptor_ok(rec["gene"], rec["filters"]) else "fail"
        anch = "pass" if gpi or (tm or 0) >= 1 else ("fail" if v1 else "no v1 topology")
        ecd = "pass" if gpi or (ecto or 0) >= F.MIN_ECD_RES else ("fail" if v1 else "no v1")
        tests = (size_ok, co == "pass", anch == "pass", ecd == "pass")
        if all(tests):
            verdict = "candidate"
        elif co == "not evaluable" and size_ok and anch == "pass" and ecd == "pass":
            verdict = "candidate pending co-receptor review"
        else:
            verdict = "not a candidate"
        failed = [
            n
            for n, ok in zip(("size", "co-receptor", "anchored", "ectodomain"), tests, strict=True)
            if not ok
        ]
        rows.append(
            {
                "uniprot_accession": acc,
                "gene": rec["gene"] if rec else g["sweep_symbol"],
                "verdict": verdict,
                "fails": "; ".join(failed),
                "snapshots": ", ".join(g["snapshots"]),
                "n_surface_pubs_current": g["n_surface_pubs"],
                "lod_membrane_current": g["lod_membrane"],
                "length_aa": len(seq) if seq else "",
                "canonical_length_check": length_check,
                "orf_nt": orf or "",
                "size_test": "pass" if size_ok else "fail",
                "coreceptor_test": co,
                "coreceptor_dependency": (rec or {})
                .get("filters", {})
                .get("co_receptor_dependency", "")
                if rec
                else "",
                "coreceptor_partners": "; ".join(
                    map(str, (rec or {}).get("co_receptor_partners") or [])
                ),
                "topology_tool": "deeptmhmm_v1" if v1 else "",
                "topology_version": t["topology_version"] if v1 else "",
                "v1_tm_helices": tm if v1 else "",
                "v1_signal_peptide": ("S" in v1) if v1 else "",
                "catalog_gpi": gpi,
                "anchored_test": anch,
                "v1_extracellular_res": ecto if v1 else "",
                "ectodomain_test": ecd,
                "ecto_vs_l1_max": f"{ecto}/{l1_max_ecto}" if v1 else "",
                "catalog_tier": (rec["filters"].get("deep_dive_tier") or "")
                if rec
                else "no deep-dive record",
                "s1e3_first_failing_gate": (F.first_failing_gate(rec) or "passes") if rec else "",
                **{f"gate:{k}": v for k, v in gates.items()},
            }
        )

    cols = list(dict.fromkeys(k for r in rows for k in r))
    with OUT.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n", restval="")
        w.writeheader()
        w.writerows(rows)
    from collections import Counter

    print(f"{OUT.name}: {len(rows)} gap accessions; L1 largest ectodomain {l1_max_ecto} aa")
    print(Counter(r["verdict"] for r in rows))
    for r in rows:
        if r["verdict"] != "not a candidate":
            print(
                f"  {r['gene']:<10} {r['verdict']:<38} orf={r['orf_nt']} tm={r['v1_tm_helices']} "
                f"ecto={r['v1_extracellular_res']} gpi={r['catalog_gpi']} co={r['coreceptor_test']}"
            )


if __name__ == "__main__":
    main()
