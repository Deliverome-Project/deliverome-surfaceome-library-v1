"""Every L1 member's coding sequence, checked against the Ensembl transcript it was taken from.

The reference-allele check (``scratch/becca/appris-principal/check_reference_allele_lof.py``)
covers only the three genes built from an alternative-haplotype or patch transcript. This
runs the same questions over the whole of L1, from pinned files only:

Per active L1 member (``library_manifest_l1.tsv``, 2,710 rows; the 25 expression-rescue arms
are constructed ORFs with no source transcript and are checked by ``check_l1_built_orfs.py``):

* the designed CDS (``library_design.CDS_KEY``) is whole codons, starts ATG, ends in a stop
  and has no internal stop;
* its translation's SHA-256 equals the manifest's ``sequence_sha256``, the protein the
  built-ORF check then follows into the synthesized construct;
* the named transcript is found in the pinned Ensembl FASTA (115; 116 for orthologs adopted
  under the Compara 116 policy), with its biotype and whether it lies on the primary
  assembly or on an alternative-haplotype / patch scaffold;
* the designed CDS equals that transcript's CDS base for base, and if not, why
  (a recorded start-codon replacement, upstream trim or added stop);
* for human members, which UniProt 2026_03 reviewed sequence the protein equals, if any.

Per library gene (human, the 2,369 genes): whether any primary-assembly Ensembl 115
transcript of the gene is ``protein_coding`` and translates without an internal stop. A gene
with none is one whose GRCh38 copy cannot be built.

Writes ``library_transcript_check.tsv`` (one row per member),
``library_transcript_genes.tsv`` (one row per library gene) and
``library_transcript_check_summary.json``; exits 1 if any member fails a CDS-integrity or
identity test.

    uv run --locked python analysis/sprint4/s4e9-orf-design-fg-panel/check_library_transcripts.py
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

from Bio.Seq import Seq

from deliverome_analysis import library_design as ld
from deliverome_analysis.s3 import external, processed

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
AP = REPO / "scratch/becca/appris-principal"
sys.path.insert(0, str(AP))

import candidate_inputs as inputs  # noqa: E402

MANIFEST = HERE / "library_manifest_l1.tsv"
REFERENCE_RECORD = (
    "sprint4/s4e9-orf-design-fg-panel/library-design/reference-selection/"
    "2026-11-12-reference-fixes-2/reference_selection.tsv"
)
OUT = HERE / "library_transcript_check.tsv"
OUT_GENES = HERE / "library_transcript_genes.tsv"
SUMMARY = HERE / "library_transcript_check_summary.json"
STOPS = {"TAA", "TAG", "TGA"}


def sha(s: str) -> str:
    return hashlib.sha256(s.encode()).hexdigest()


def read_fasta(path: Path) -> dict[str, tuple[str, str]]:
    out, name, head, buf = {}, None, None, []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    out[name] = (head, "".join(buf).upper())
                head = line[1:].strip()
                name = head.split()[0].split(".")[0]
                buf = []
            else:
                buf.append(line.strip())
    if name:
        out[name] = (head, "".join(buf).upper())
    return out


def field(head: str, key: str) -> str:
    m = re.search(rf"\b{key}:(\S+)", head)
    return m.group(1) if m else ""


def translate(cds: str) -> str:
    return str(Seq(cds[: len(cds) // 3 * 3]).translate())


def load_uniprot() -> tuple[dict[str, list[str]], dict[str, str]]:
    """(protein sequence -> accessions, accession -> sequence) from the pinned 2026_03 tables."""
    by_seq: dict[str, list[str]] = defaultdict(list)
    by_acc: dict[str, str] = {}
    for path in inputs.uniprot_tables("human"):
        with gzip.open(path, "rt") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                acc, seq = row["Entry"], row["Sequence"].strip().upper()
                if acc in by_acc:
                    continue
                by_acc[acc] = seq
                by_seq[seq].append(acc)
    return by_seq, by_acc


def cds_fasta_116(species: str) -> Path:
    stem = inputs.ASSEMBLIES[species]
    key = f"ensembl/release-116/cds/{stem}.cds.all.fa.gz"
    dest = inputs.CACHE / key.replace("/", "__")
    if not dest.exists():
        dest.parent.mkdir(parents=True, exist_ok=True)
        external.download(key, dest)
    return dest


def main() -> int:
    with MANIFEST.open() as fh:
        members = list(csv.DictReader(fh, delimiter="\t"))
    designed = ld.load_designed_cds()
    with open(processed.download(REFERENCE_RECORD)) as fh:
        refsel = {r["symbol"]: r for r in csv.DictReader(fh, delimiter="\t")}
    uni_by_seq, uni_by_acc = load_uniprot()

    fasta = {
        ("human", "115"): read_fasta(inputs.cds_fasta("human")),
        ("mouse", "115"): read_fasta(inputs.cds_fasta("mouse")),
        ("cyno", "115"): read_fasta(inputs.cds_fasta("cyno")),
        ("mouse", "116"): read_fasta(cds_fasta_116("mouse")),
        ("cyno", "116"): read_fasta(cds_fasta_116("cyno")),
    }

    rows, failures, flags = [], [], []
    for m in members:
        rec = designed[m["manifest_key"]]
        cds = rec["cds"].upper()
        prot = translate(cds)
        body = prot[:-1] if prot.endswith("*") else prot
        enst = m["transcript"].split(".")[0]
        release = ""
        src = None
        for rel in ("115", "116"):
            hit = fasta.get((m["species"], rel), {}).get(enst)
            if hit:
                release, src = rel, hit
                break
        head, src_cds = src if src else ("", "")
        src_prot = translate(src_cds) if src_cds else ""
        src_body = src_prot[:-1] if src_prot.endswith("*") else src_prot
        region = head.split()[2] if head else ""
        edits = [
            k
            for k in ("start_codon_replaced", "trimmed_upstream_codons", "stop_added")
            if rec.get(k)
        ]
        if not src_cds:
            vs_source = "transcript not in pinned FASTA"
        elif src_cds == cds:
            vs_source = "identical"
        elif edits:
            vs_source = "differs by recorded edit: " + ",".join(edits)
        elif src_body == body:
            vs_source = "same protein, different codons"
        else:
            vs_source = "differs"
        accs = uni_by_seq.get(body, []) if m["species"] == "human" else []
        canonical_acc = (refsel.get(m["gene_symbol"]) or {}).get("current_accession", "")
        row = {
            "design_id": m["candidate_id"],
            "key": m["manifest_key"],
            "set": m["candidate_set"],
            "species": m["species"],
            "gene": m["gene_symbol"],
            "transcript": enst,
            "route": m["route"],
            "ensembl_release": release,
            "transcript_biotype": field(head, "transcript_biotype"),
            "region": region,
            "primary_assembly": region.startswith("chromosome:") if region else "",
            "cds_nt": len(cds),
            "whole_codons": len(cds) % 3 == 0,
            "starts_atg": cds[:3] == "ATG",
            "ends_in_stop": cds[-3:] in STOPS,
            "internal_stops": body.count("*"),
            "protein_equals_manifest": sha(body) == m["sequence_sha256"],
            "cds_vs_source_transcript": vs_source,
            "source_internal_stops": src_body.count("*") if src_cds else "",
            "variant_residues_vs_uniprot": rec.get("variant_residues", ""),
            "uniprot_identical": ";".join(accs),
            "equals_gene_uniprot_canonical": (
                bool(canonical_acc) and uni_by_acc.get(canonical_acc) == body
                if m["species"] == "human"
                else ""
            ),
        }
        rows.append(row)
        problems = []
        if not row["whole_codons"]:
            problems.append("not whole codons")
        if not row["ends_in_stop"]:
            problems.append("no terminal stop")
        if row["internal_stops"]:
            problems.append("internal stop")
        if not row["protein_equals_manifest"]:
            problems.append("protein differs from manifest")
        if vs_source in ("differs", "transcript not in pinned FASTA"):
            problems.append(f"CDS vs transcript: {vs_source}")
        if row["transcript_biotype"] not in ("protein_coding", ""):
            # Not a sequence defect (the CDS is tested above); listed for review.
            flags.append(
                f"{row['design_id']} {row['gene']}: source biotype {row['transcript_biotype']}"
            )
        if problems:
            failures.append(f"{row['design_id']} {row['gene']}: {'; '.join(problems)}")

    # Gene level: can the GRCh38 primary-assembly copy of each library gene be built? A gene
    # is found by every Ensembl gene id its members' transcripts and its MANE / APPRIS picks
    # carry, and by its symbol, because the symbols of recently renamed genes (the new SLC
    # names, ACKR5) differ between HGNC and Ensembl 115.
    human = fasta[("human", "115")]
    by_symbol: dict[str, set[str]] = defaultdict(set)
    by_ensg: dict[str, set[str]] = defaultdict(set)
    for enst, (head, _seq) in human.items():
        by_symbol[field(head, "gene_symbol")].add(enst)
        by_ensg[field(head, "gene").split(".")[0]].add(enst)

    def ensg(enst: str) -> str:
        hit = human.get(enst.split(".")[0])
        return field(hit[0], "gene").split(".")[0] if hit else ""

    gene_rows = []
    for gene in sorted(refsel):
        own = [r["transcript"] for r in rows if r["gene"] == gene and r["species"] == "human"]
        picks = [t for k in ("mane_enst", "appris_enst") for t in refsel[gene][k].split(";") if t]
        ids = {ensg(t) for t in own + picks} - {""}
        tx_ids = set(by_symbol.get(gene, set()))
        for i in ids:
            tx_ids |= by_ensg[i]
        tx = [(t, *human[t]) for t in sorted(tx_ids)]
        prim = [t for t in tx if t[1].split()[2].startswith("chromosome:")]
        translates = [
            t for t in prim if len(t[2]) % 3 == 0 and "*" not in translate(t[2]).rstrip("*")
        ]
        coding = [t for t in translates if field(t[1], "transcript_biotype") == "protein_coding"]
        if not prim:
            status = "absent from primary assembly"
        elif coding:
            status = "intact protein_coding copy"
        elif translates:
            status = "translates intact, but labeled " + ",".join(
                sorted({field(t[1], "transcript_biotype") for t in translates})
            )
        else:
            status = "no intact copy"
        built = sorted(
            {r["design_id"] for r in rows if r["gene"] == gene and r["species"] == "human"}
        )
        gene_rows.append(
            {
                "gene": gene,
                "ensembl_genes": ";".join(sorted(ids)),
                "ensembl_transcripts": len(tx),
                "primary_assembly_transcripts": len(prim),
                "primary_protein_coding_LoF": sum(
                    1 for t in prim if field(t[1], "transcript_biotype") == "protein_coding_LoF"
                ),
                "primary_translating_intact": len(translates),
                "primary_status": status,
                "l1_human_members": ";".join(built),
                "built_from_non_primary": ";".join(
                    sorted(
                        r["design_id"]
                        for r in rows
                        if r["gene"] == gene
                        and r["species"] == "human"
                        and r["primary_assembly"] is False
                    )
                ),
            }
        )

    for path, data in ((OUT, rows), (OUT_GENES, gene_rows)):
        with path.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(data[0]), delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(data)

    human_rows = [r for r in rows if r["species"] == "human"]
    ref_rows = [r for r in rows if r["set"] == "library_reference"]
    summary = {
        "members_checked": len(rows),
        "by_set_species": dict(Counter(f"{r['set']}:{r['species']}" for r in rows)),
        "cds_integrity_failures": len(failures),
        "failures": failures,
        "flagged_source_biotype": flags,
        "ensembl_release_of_source": dict(Counter(r["ensembl_release"] for r in rows)),
        "source_transcript_biotype": dict(Counter(r["transcript_biotype"] for r in rows)),
        "cds_vs_source_transcript": dict(Counter(r["cds_vs_source_transcript"] for r in rows)),
        "human_members_off_primary_assembly": sorted(
            f"{r['design_id']} {r['gene']} {r['region'].split(':')[2]}"
            for r in human_rows
            if r["primary_assembly"] is False
        ),
        "human_members": len(human_rows),
        "human_members_identical_to_a_uniprot_reviewed_sequence": sum(
            1 for r in human_rows if r["uniprot_identical"]
        ),
        "human_members_variant_residues_vs_uniprot": dict(
            Counter(r["variant_residues_vs_uniprot"] for r in human_rows)
        ),
        "references": len(ref_rows),
        "references_equal_gene_uniprot_canonical": sum(
            1 for r in ref_rows if r["equals_gene_uniprot_canonical"] is True
        ),
        "library_genes": len(gene_rows),
        "genes_with_primary_LoF_transcript": sum(
            1 for g in gene_rows if g["primary_protein_coding_LoF"]
        ),
        "genes_by_primary_status": dict(Counter(g["primary_status"] for g in gene_rows)),
        "genes_without_intact_protein_coding_primary_copy": [
            f"{g['gene']}: {g['primary_status']} (members: {g['l1_human_members'] or 'none'})"
            for g in gene_rows
            if g["primary_status"] != "intact protein_coding copy"
        ],
        "inputs": {
            "manifest": MANIFEST.name,
            "designed_cds": f"processed://{ld.CDS_KEY}",
            "reference_selection": f"processed://{REFERENCE_RECORD}",
            **inputs.provenance(),
            "cds_mouse_116": "external://ensembl/release-116/cds/Mus_musculus.GRCm39.cds.all.fa.gz",
            "cds_cyno_116": (
                "external://ensembl/release-116/cds/"
                "Macaca_fascicularis.Macaca_fascicularis_6.0.cds.all.fa.gz"
            ),
        },
    }
    SUMMARY.write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps({k: v for k, v in summary.items() if k != "inputs"}, indent=1))
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
