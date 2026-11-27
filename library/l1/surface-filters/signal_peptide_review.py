"""Reproduce the frozen S1E3 SP prioritization; no live catalog or model calls.

Run with --fetch to download checksum-pinned inputs from processed S3.
Default execution reads only the explicitly supplied input directory.
"""

import argparse
import collections
import csv
import gzip
import hashlib
import json
import pathlib
import re

MANIFEST = json.loads(pathlib.Path(__file__).with_name("signal_peptide_snapshot.json").read_text())


def read_tedman(path):
    """Read original workbook, retaining 1-based Excel rows and nulls."""
    from openpyxl import load_workbook

    wb = load_workbook(path, read_only=True, data_only=True)
    try:
        sheet = wb["Canonical"]
        values = sheet.iter_rows(values_only=True)
        headers = next(values)
        for index, values_row in enumerate(values, start=2):
            row = dict(zip(headers, values_row, strict=True))
            if row.get("Receptor Name"):
                row["excel_row"] = index
                yield row
    finally:
        wb.close()


def prepare_inputs(directory, fetch=False):
    """Refuse changed inputs; download only on explicit --fetch."""
    if fetch:
        from deliverome_analysis.core.s3 import processed

        directory.mkdir(parents=True, exist_ok=True)
        for name, metadata in MANIFEST["inputs"].items():
            processed.download(metadata["s3_key"], directory / name)
    for name, metadata in MANIFEST["inputs"].items():
        observed = hashlib.sha256((directory / name).read_bytes()).hexdigest()
        if observed != metadata["sha256"]:
            raise ValueError(f"Input checksum mismatch: {name}")


def analyze(SRC, OUT):
    members = list(csv.DictReader((SRC / "library_membership.tsv").open(), delimiter="\t"))
    cohort = {r["gene"]: r for r in map(json.loads, gzip.open(SRC / "cohort.jsonl.gz", "rt"))}
    tedman = {}
    for row in read_tedman(SRC / "tedman-v1.xlsx"):
        gene = row["gene_id"] or row["Receptor Name"].split("_ENST")[0]
        gene = {"GPR182": "ACKR5"}.get(gene, gene)
        assert gene not in tedman
        tedman[gene] = row

    direct = {
        "ADRB2": (
            "Human receptor, insect-cell expression: improved functional expression; host transfer remains inference.",
            "https://pubmed.ncbi.nlm.nih.gov/1331042/",
        ),
        "CNR1": (
            "Added SP improved stability and surface expression; direct mechanistic rescue of long N-tail translocation.",
            "https://pubmed.ncbi.nlm.nih.gov/12920192/",
        ),
        "GIPR": (
            "Native-SP replacement plus FLAG tag improved membrane colocalization; tag is a co-modification.",
            "https://www.frontiersin.org/journals/pharmacology/articles/10.3389/fphar.2020.01271/full",
        ),
        "GRIK1": (
            "SP replacement improved surface/synaptic trafficking in rodent neurons; human transfer remains inference.",
            "https://www.nature.com/articles/s41467-018-07403-7",
        ),
    }
    feasible = {
        "TLR4": (
            "Functional surface expression with RAGE SP; superiority to native SP not established.",
            "https://pmc.ncbi.nlm.nih.gov/articles/PMC2774436/",
        ),
        "KCNE3": (
            "Human KCNE3 with IL2-SP and mCherry formed functional KCNQ1 complex; SP-only gain not established.",
            "https://pmc.ncbi.nlm.nih.gov/articles/PMC7083075/",
        ),
        "PRNP": (
            "SP swaps improve translocation; steady-state surface expression unchanged in unstressed cells.",
            "https://pmc.ncbi.nlm.nih.gov/articles/PMC533048/",
        ),
    }
    negative = {"DRD1", "DRD2", "HTR2A", "AVPR2", "GCGR", "GLP1R"}
    rows = []
    for member in members:
        g = member["gene"]
        r = cohort[g]
        t = r["canonical_topology"]
        sub = r["filters"]["subcategory"]
        fam = r["uniprot_family"] or ""
        nt = t["n_terminal_orientation"]
        sp = t["signal_peptide_length"]
        tr = tedman.get(g)
        bona = sub == "GPCR" and fam.startswith(("G protein-coupled", "chemokine-like"))
        sensory = g.startswith(("TAS", "VN1R")) or "Opsin subfamily" in fam
        source = ""
        note = ""
        review = ""
        if g in {"MRAP", "MRAP2"}:
            tier = "D"
            cls = "Dual topology accessory proteins"
            note = "Do not force a single N-out topology with blanket SP addition."
            source = "https://pmc.ncbi.nlm.nih.gov/articles/PMC2992810/"
        elif nt != "extracellular":
            tier = "D"
            cls = "N-in or conflicting topology"
            note = "No blanket cleavable-SP addition; native N-terminal cytoplasmic orientation may be required."
            if g == "FZD1":
                review = "GPCR family versus N-in annotation conflict; manual sequence review before eligibility decision."
        elif g in direct:
            tier = "A"
            cls = "Direct improvement precedent"
            note, source = direct[g]
        elif re.match(r"^OR\d", g):
            tier = "B"
            cls = "Olfactory receptors"
            note = "Lucy cleavable-SP family precedent, mainly mouse receptors; benefit not demonstrated for every human member or without accessory proteins."
            source = "https://pmc.ncbi.nlm.nih.gov/articles/PMC3698168/"
        elif g in feasible:
            tier = "B"
            cls = "Specific feasibility or translocation precedent"
            note, source = feasible[g]
        elif bona and not sensory:
            tier = "B"
            cls = "Other non-sensory GPCRs"
            note = "Broad GPCR library engineering precedent; confidence in worth testing, not in an expression gain."
            source = "https://pmc.ncbi.nlm.nih.gov/articles/PMC11290489/"
        else:
            tier = "C"
            if bona and sensory:
                cls = "Taste, vomeronasal and opsin GPCRs"
                note = "Compatible N-out topology; no clean general cleavable-SP rescue precedent established in this review."
            elif sub == "GPI_anchored":
                cls = "Other GPI anchored proteins"
                note = "Plausible SP replacement; preserve C-terminal GPI processing signal. ER entry may not be limiting."
            elif sub == "single_pass_T1" and sp:
                cls = "Other type I proteins with native SP"
                note = "Plausible SP replacement; native-SP annotation establishes compatibility, not a defective native signal."
            elif sub == "single_pass_T1":
                cls = "Other N-out single-pass proteins without SP"
                note = "Exploratory addition; verify native topology and mature N terminus before design."
            else:
                cls = "Other N-out membrane proteins"
                note = "Topology-only hypothesis; multipass insertion/assembly or accessory proteins may dominate trafficking."
        if g in negative:
            note += " Comparative HA-SP study did not establish general functional improvement; do not treat as a rescue positive."
            source += " https://pmc.ncbi.nlm.nih.gov/articles/PMC9855867/"
        if g == "DRD2":
            note += " Lucy rescue reported for a deletion mutant, not intact full-length receptor."
        if g == "CNR2":
            note += " Direct human preprotrypsin-SP expression precedent, with weaker surface-specific evidence."
            source += " https://pmc.ncbi.nlm.nih.gov/articles/PMC3703952/"
        if g == "OR51E2":
            note += " Directly tested with Lucy but already surface-positive without SP; no gene-specific surface-rescue gain established."
        if g == "MC4R":
            note += " Functional prolactin-SP construct precedent, not proof of increased surface abundance."
        if g in {"FZD4", "FZD5", "FZD6", "FZD10"}:
            note += " DuET used a replacement 5HT3A signal peptide."
        if g == "AMHR2":
            review = "Type II category conflicts with N-out/native-SP topology; manual annotation review."
        if sub == "GPCR" and not bona:
            review = "Catalog GPCR label lacks corroborating GPCR family annotation; not included in GPCR-precedent class."
        val = tr["Immunostaining Intensity"] if tr else None
        rows.append(
            dict(
                gene=g,
                library_tier=member["tier"],
                confidence_tier=tier,
                candidate_class=cls,
                design=(
                    "hold for topology review"
                    if tier == "D"
                    else "replace native SP"
                    if sp
                    else "add cleavable SP"
                ),
                native_sp_length=sp,
                n_orientation=nt,
                tm_count=t["tm_helix_count"],
                tedman_intensity=val,
                tedman_sd=tr["Immunostaining Standard Deviation"] if tr else None,
                tedman_below_1000=None if val is None else val < 1000,
                tedman_construct=tr["Receptor Name"] if tr else None,
                tedman_excel_row=tr["excel_row"] if tr else None,
                rationale=note,
                review_note=review,
                evidence_url=source,
            )
        )
    assert len(rows) == len({r["gene"] for r in rows}) == 2344
    assert sum(r["tedman_intensity"] is not None for r in rows) == 528
    assert sum(r["tedman_below_1000"] is True for r in rows) == 253
    summary = []
    for (tier, cls), count in sorted(
        collections.Counter((r["confidence_tier"], r["candidate_class"]) for r in rows).items()
    ):
        rr = [r for r in rows if (r["confidence_tier"], r["candidate_class"]) == (tier, cls)]
        summary.append(
            dict(
                tier=tier,
                classification=cls,
                n=count,
                tedman_measured=sum(r["tedman_intensity"] is not None for r in rr),
                tedman_low=sum(r["tedman_below_1000"] is True for r in rr),
                addition=sum(r["design"] == "add cleavable SP" for r in rr),
                replacement=sum(r["design"] == "replace native SP" for r in rr),
            )
        )
    OUT.mkdir(parents=True, exist_ok=True)
    payload = {
        "analysis_date": "2026-09-11",
        "membership_source": MANIFEST["inputs"]["library_membership.tsv"]["s3_key"],
        "membership_sha256": hashlib.sha256(
            (SRC / "library_membership.tsv").read_bytes()
        ).hexdigest(),
        "cohort_sha256": hashlib.sha256((SRC / "cohort.jsonl.gz").read_bytes()).hexdigest(),
        "tedman_version": "2025-09-19 v1 Doc S1 Canonical",
        "summary": summary,
        "proteins": rows,
    }
    (OUT / "s1e3-signal-peptide-candidates.json").write_text(json.dumps(payload, indent=2) + "\n")

    expected = MANIFEST["expected"]
    selections = {
        "all": rows,
        "precedent": [r for r in rows if r["confidence_tier"] in {"A", "B"}],
        "broad": [r for r in rows if r["confidence_tier"] in {"A", "B", "C"}],
        "low_tedman": [r for r in rows if r["tedman_below_1000"] is True],
    }
    for name, subset in selections.items():
        if len(subset) != expected[name]:
            raise ValueError(f"Unexpected {name} count: {len(subset)}")
        with (OUT / f"signal_peptide_{name}.tsv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(sorted(subset, key=lambda r: r["gene"]))
    observed = {
        "measured": sum(r["tedman_intensity"] is not None for r in rows),
        "native_sp": sum(r["design"] == "replace native SP" for r in rows),
        "addition": sum(r["design"] == "add cleavable SP" for r in rows),
    }
    for name, value in observed.items():
        if value != expected[name]:
            raise ValueError(f"Unexpected {name} count: {value}")
    (OUT / "signal_peptide_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({name: len(subset) for name, subset in selections.items()}))
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=pathlib.Path, required=True)
    parser.add_argument("--output-dir", type=pathlib.Path, required=True)
    parser.add_argument("--fetch", action="store_true")
    args = parser.parse_args()
    prepare_inputs(args.input_dir, fetch=args.fetch)
    analyze(args.input_dir, args.output_dir)


if __name__ == "__main__":
    main()
