"""How each library gene without a buildable GRCh38 reference copy entered the gene list, and why.

Two groups are covered, both from pinned inputs only:

* the 9 library genes with no L1 member of their own: 4 paralogs whose identical protein is
  built under another symbol, and 5 whose only Ensembl 115 copies are loss-of-function;
* the 4 genes whose GRCh38 primary-assembly copy is loss-of-function but which are built from
  an intact copy on an alternative-haplotype or patch scaffold (KIR2DS4, OR12D1, OR5AL1,
  OR11H7). ``check_library_transcripts.py`` finds these library-wide; this explains them.

For each gene: the route into the list (the core gates, or which named tier and on what
evidence), the reference the ladder chose, and what blocks or replaces the build.

The entry route is recomputed by running the S1E3 filter offline over accessible-surfaceome
release 1.3.0 (``library_pipeline.run_gene_filter``), so this needs S3 read access and a few
minutes; pass ``--s1e3`` to reuse a staged run.

    uv run --locked python analysis/sprint4/s4e9-orf-design-fg-panel/genes_without_reference_cds.py
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
AP = REPO / "scratch/becca/appris-principal"
OUT = HERE / "genes_without_reference_cds.tsv"
REFERENCE_RECORD = (
    "sprint4/s4e9-orf-design-fg-panel/library-design/reference-selection/"
    "2026-11-12-reference-fixes-2/reference_selection.tsv"
)

PARALOGS = {"OPN1MW2": "OPN1MW", "OPN1MW3": "OPN1MW", "OR2A42": "OR2A1", "UPK3BL2": "UPK3BL1"}
NO_INTACT_COPY = ["FCGR2C", "OR8J2", "OR10AC1", "OR52P1", "OR4Q2"]
BUILT_OFF_PRIMARY = ["KIR2DS4", "OR12D1", "OR5AL1", "OR11H7"]

TIER_NAMES = {
    "T1": "core gates (architecture tier T1)",
    "T2": "core gates (architecture tier T2)",
    "T3": "core gates (architecture tier T3)",
    "LL": "low-literature UniProt surface tier",
    "SV": "screen-validated GPCR tier",
    "OE": "overexpression surface tier",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with Path(path).open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def ensembl_copies(genes: list[str]) -> dict[str, list[dict]]:
    """Every pinned Ensembl 115 CDS of each gene, translated, with where it lies."""
    sys.path.insert(0, str(AP))
    import candidate_inputs as inputs
    from Bio.Seq import Seq

    out: dict[str, list[dict]] = {g: [] for g in genes}
    head, buf = None, []

    def flush():
        if not head:
            return
        sym = re.search(r"gene_symbol:(\S+)", head)
        if not sym or sym.group(1) not in out:
            return
        cds = "".join(buf).upper()
        prot = str(Seq(cds[: len(cds) // 3 * 3]).translate()).rstrip("*")
        region = head.split()[2]
        out[sym.group(1)].append(
            {
                "transcript": head.split()[0].split(".")[0],
                "biotype": re.search(r"transcript_biotype:(\S+)", head).group(1),
                "primary": region.startswith("chromosome:"),
                "scaffold": region.split(":")[2],
                "protein": prot,
                "stops": [i + 1 for i, c in enumerate(prot) if c == "*"],
            }
        )

    with gzip.open(inputs.cds_fasta("human"), "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                flush()
                head, buf = line[1:].strip(), []
            else:
                buf.append(line.strip())
    flush()
    return out


def uniprot_sequences() -> dict[str, str]:
    sys.path.insert(0, str(AP))
    import candidate_inputs as inputs

    out = {}
    for path in inputs.uniprot_tables("human"):
        with gzip.open(path, "rt") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                out.setdefault(row["Entry"], row["Sequence"].strip().upper())
    return out


def describe(copies: list[dict], target: str | None) -> str:
    """One clause per copy: where it lies, its biotype, and how its translation fails."""
    parts = []
    for c in copies:
        where = "GRCh38 primary" if c["primary"] else c["scaffold"]
        if c["protein"] == target:
            state = "identical to UniProt"
        elif c["stops"]:
            first = c["stops"][0]
            ref = f" (UniProt {target[first - 1]})" if target and len(target) >= first else ""
            state = f"{len(c['stops'])} internal stop(s), first at codon {first}{ref}"
        elif target and len(c["protein"]) == len(target):
            n = sum(1 for x, y in zip(c["protein"], target, strict=True) if x != y)
            state = f"no stop; {n} residue(s) differ from UniProt"
        else:
            state = f"no stop; {len(c['protein'])} aa vs UniProt {len(target or '')} aa"
        parts.append(f"{c['transcript']} ({where}, {c['biotype']}): {state}")
    return "; ".join(parts)


def staged_s1e3(arg: Path | None) -> Path:
    if arg:
        return arg.resolve()
    sys.path.insert(0, str(HERE))
    import library_pipeline

    work = Path(tempfile.mkdtemp(prefix="gene-filter-"))
    library_pipeline.run_gene_filter(work)
    return work / "s1e3"


def entry_basis(F, r: dict, tier: str, sv: dict) -> str:
    f, g = r["filters"], r["gene"]
    catalog = (
        "catalog canonical tier"
        if F.B.passes_canonical(f)
        else "catalog likely tier"
        if F.B.passes_likely(f)
        else f"not in the catalog's canonical or likely tier (deep_dive_tier={f.get('deep_dive_tier')})"
    )
    tpm = F.mx(g)
    tpm_s = f"max TPM {tpm:.1f}" if tpm is not None else "max TPM not recorded"
    orf = F.orf_nt(r)
    if tier in ("T1", "T2", "T3"):
        expr = (
            "expression gate passed on the catalog's low_endogenous_expression=False"
            if f.get("low_endogenous_expression") is False
            else "expression gate passed on max TPM >= 5"
        )
        return f"{catalog}; passes every gate ({expr}); {tpm_s}; canonical ORF {orf} nt"
    if tier == "LL":
        return (
            f"{catalog}; UniProt places it at the cell surface; {f.get('n_papers_found')} papers "
            f"(< 100); anchored ({(r.get('canonical_topology') or {}).get('tm_helix_count')} "
            f"TM helices); {tpm_s} (>= 5); canonical ORF {orf} nt (<= 5,000)"
        )
    if tier == "SV":
        s = sv[g]
        return (
            f"{catalog}; anti-HA surface immunostaining {s['immunostaining']:,.1f} a.u. "
            f"(>= 1,000) in the GPCRome screen; failed core gates {', '.join(s['prior_fail_gates'])}"
        )
    return catalog


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--s1e3", type=Path, help="a staged S1E3 folder from run_gene_filter")
    a = ap.parse_args(argv)
    s1e3 = staged_s1e3(a.s1e3)
    sys.path.insert(0, str(s1e3))
    import final_library as F  # importing runs the filter

    tiers = {r["gene"]: r["tier"] for r in read_tsv(s1e3 / "library_membership.tsv")}
    sv = json.loads((s1e3 / "willow_screen_rescue_genes.json").read_text())
    recs = {r["gene"]: r for r in F.recs}

    from deliverome_analysis.s3 import processed

    refsel = {r["symbol"]: r for r in read_tsv(processed.download(REFERENCE_RECORD))}
    manifest = {r["key"]: r for r in read_tsv(HERE / "candidate_manifest_tranches.tsv")}
    active = {
        r["gene_symbol"]: r
        for r in read_tsv(HERE / "library_manifest_l1.tsv")
        if r["candidate_set"] == "library_reference"
    }
    copies = ensembl_copies([*NO_INTACT_COPY, *BUILT_OFF_PRIMARY])
    uniprot = uniprot_sequences()
    tx_check = {
        r["gene"]: r
        for r in read_tsv(HERE / "library_transcript_check.tsv")
        if r["set"] == "library_reference"
    }

    rows = []
    for gene in [*PARALOGS, *NO_INTACT_COPY, *BUILT_OFF_PRIMARY]:
        r, ref, tier = recs[gene], refsel[gene], tiers[gene]
        f = r["filters"]
        notes = " ".join(str(v) for k, v in f.items() if k.endswith("_rationale") and v).lower()
        if gene in PARALOGS:
            m = manifest[f"lib:{gene}"]
            group = "identical paralog"
            why = (
                f"its reference ({ref['selected_source']}, {ref['mane_enst']}) resolves to a CDS "
                f"whose protein is identical to {PARALOGS[gene]}'s; design {m['design_id']} is "
                f"retired ({m['retired_reason']}) and the protein is built once, as "
                f"{active[PARALOGS[gene]]['candidate_id']}"
            )
            built = active[PARALOGS[gene]]["candidate_id"] + f" ({PARALOGS[gene]})"
        elif gene in NO_INTACT_COPY:
            group = "no intact copy in Ensembl 115"
            why = describe(copies[gene], uniprot.get(r.get("uniprot") or ""))
            built = "none"
        else:
            t = tx_check[gene]
            group = "primary copy loss-of-function; built from a non-primary scaffold"
            why = (
                describe(
                    [c for c in copies[gene] if c["primary"]], uniprot.get(r.get("uniprot") or "")
                )
                + f"; built from {t['transcript']} on {t['region'].split(':')[2]}"
            )
            built = f"{t['design_id']} (= {t['uniprot_identical'] or 'no UniProt match'})"
        rows.append(
            {
                "gene": gene,
                "uniprot": r.get("uniprot") or "",
                "group": group,
                "entered_by": TIER_NAMES[tier],
                "entry_evidence": entry_basis(F, r, tier, sv),
                "catalog_notes_pseudogene": "pseudogene" in notes,
                "ladder_decision": ref["decision"],
                "reference_source": ref["selected_source"],
                "reference_transcript": ref["mane_enst"] or ref["appris_enst"],
                "why": why,
                "built_as": built,
            }
        )

    with OUT.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    for row in rows:
        print(f"{row['gene']:8} {row['entered_by']:38} {row['group']}")
    print(f"wrote {OUT.name}: {len(rows)} genes")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
