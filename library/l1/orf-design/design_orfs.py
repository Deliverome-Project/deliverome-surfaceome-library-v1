"""Golden Gate domestication of the nine non-Tedman ORFs on Brenna's 15-receptor panel.

Scope: the panel minus TFRC (already ordered) and minus the five receptors that
exist as Tedman constructs (CNR1, AVPR2, ACKR3, CCR5, OR51E2 — reuse his designs
verbatim; his ORFs are nucleotide-identical to today's Ensembl canonical CDS).

What this does, per gene:

1. Canonicalize the symbol through ``deliverome_analysis.surfaceome.gene_resolver`` (MyGene),
   giving the Ensembl **gene** id. Symbols are never the join key.
2. Resolve that gene to its Ensembl **canonical transcript** and pull the native CDS
   (``deliverome_analysis.cloning.seqtools``).
3. Cross-check the translated protein against the reviewed UniProt entry reached from
   the same ENSG — identity of the *sequence*, not agreement of two labels. A mismatch
   is a hard stop, not a warning.
4. Remove every Type IIS recognition site with the **fewest possible nucleotide edits**
   (``seqtools.domesticate`` → DNA Chisel ``AvoidChanges``). No codon optimization:
   native codon usage is load-bearing for co-translational folding in multi-pass
   membrane proteins (s1e2 §4.2), and it is also what Tedman kept.
5. Gate on the protein being byte-for-byte identical before and after.

The enzyme set is a real decision and is reported per enzyme so the cost of each is
visible: BbsI is the Fragmid assembly enzyme (s1e2 §1), Esp3I/BsmBI is what Tedman used
as his shuttle enzyme and is contested for guide cloning (s1e2 §4.11), BsaI was never
cleared from the backbone. Scrubbing all three keeps every downstream option open;
``--enzymes`` narrows it.

Output is the domesticated CDS only — the tag cassette is NOT added here. Tag placement
is topology-dependent and is a separate decision (TFRC is the only panel member with an
extracellular C-terminus).
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
import time
import urllib.parse
import urllib.request
from datetime import UTC, datetime
from importlib.metadata import version
from pathlib import Path

from deliverome_analysis.cloning import seqtools
from deliverome_analysis.surfaceome import gene_resolver

# Brenna's panel minus TFRC (already ordered as p003).
PANEL = [
    "HEG1",
    "ERBB2",
    "AXL",
    "SLC9A3",
    "VSIG10",
    "NECTIN4",
    "SLC2A4",
    "MS4A1",
    "FOLR1",
    "PTPRC",
    # The five Tedman receptors. These were initially listed as "reuse verbatim", but
    # that was wrong for our backbone: Tedman silenced **BsmBI/Esp3I only**, so every
    # one of them still carries BbsI and/or BsaI sites (12 across the five; ACKR3 alone
    # has 6). BbsI is our assembly enzyme, so they have to be domesticated like any other
    # ORF. The protein and Tedman's tag placement are unaffected -- domestication is
    # synonymous-only, and he kept native codons outside his own BsmBI edits too.
    "CNR1",
    "AVPR2",
    "ACKR3",
    "CCR5",
    "OR51E2",
]

DEFAULT_ENZYMES = ["BbsI", "Esp3I", "BsaI"]

UNIPROT = "https://rest.uniprot.org/uniprotkb/search"


def uniprot_by_ensg(ensg: str) -> dict:
    """Reviewed human UniProt entry cross-referenced to this Ensembl gene id.

    Joins on the stable id, never the display symbol. Returns accession, length,
    sequence and the signal-peptide / topology features that decide tag placement.
    """
    query = f"(xref:ensembl-{ensg}) AND (reviewed:true) AND (organism_id:9606)"
    url = (
        UNIPROT
        + "?"
        + urllib.parse.urlencode(
            {
                "query": query,
                "fields": "accession,length,sequence,ft_signal,ft_topo_dom,ft_lipid",
                "format": "json",
                "size": "5",
            }
        )
    )
    request = urllib.request.Request(url, headers={"User-Agent": seqtools.USER_AGENT})
    with urllib.request.urlopen(request) as response:
        payload = json.load(response)
    hits = payload.get("results", [])
    if not hits:
        raise RuntimeError(f"no reviewed UniProt entry cross-referenced to {ensg}")
    hit = hits[0]
    features = hit.get("features", [])

    def spans(kind):
        return [f for f in features if f["type"] == kind]

    signal = spans("Signal")
    topo = spans("Topological domain")
    lipid = spans("Lipidation")
    return {
        "accession": hit["primaryAccession"],
        "length": hit["sequence"]["length"],
        "sequence": hit["sequence"]["value"],
        "signal_peptide": f"1-{signal[0]['location']['end']['value']}" if signal else "",
        "n_terminus": topo[0]["description"] if topo else "",
        "c_terminus": topo[-1]["description"] if topo else "",
        "gpi": ";".join(f["description"] for f in lipid if "GPI" in f.get("description", "")),
    }


def edit_positions(before: str, after: str) -> list[int]:
    """0-based positions where domestication changed a base."""
    if len(before) != len(after):
        raise RuntimeError("domestication changed sequence length")
    return [i for i, (a, b) in enumerate(zip(before, after, strict=True)) if a != b]


def read_fasta_records(path: Path) -> dict[str, str]:
    """header -> sequence, keeping the full pipe-delimited header."""
    out, key = {}, None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            key = line[1:]
            out[key] = ""
        elif key:
            out[key] += line.strip()
    return out


def verify(out_dir: Path, enzymes: list[str], seed: int) -> int:
    """Re-domesticate every committed native CDS and demand a byte-identical result.

    This is the reproducibility claim as a runnable check rather than an assertion in a
    README. It needs no network: the native CDS is the committed input, so a failure here
    means the *solver* drifted (a dnachisel or numpy upgrade), which is exactly the change
    a stable set of versioned Ensembl transcript ids cannot catch. Re-fetching instead
    would conflate that with Ensembl having revised the transcript.

    Writes reproducibility.json recording the versions this ran under and a checksum per
    gene, so a future mismatch can be attributed instead of merely detected.
    """
    native = read_fasta_records(out_dir / "orf_native.fasta")
    built = read_fasta_records(out_dir / "orf_domesticated.fasta")
    by_gene = {h.split("|")[0]: s for h, s in built.items()}

    rows, bad = [], []
    for header, seq in sorted(native.items()):
        gene = header.split("|")[0]
        if gene not in by_gene:
            bad.append(f"{gene}: in orf_native.fasta but not orf_domesticated.fasta")
            continue
        redone, _warn = seqtools.domesticate(seq, enzymes, preserve_translation=True, seed=seed)
        ok = redone == by_gene[gene]
        if not ok:
            n = sum(1 for a, b in zip(redone, by_gene[gene], strict=False) if a != b)
            bad.append(f"{gene}: {n} nt differ from the committed design")
        rows.append(
            {
                "gene": gene,
                "transcript": header.split("|")[2] if header.count("|") >= 2 else "",
                "native_sha256": hashlib.sha256(seq.encode()).hexdigest(),
                "domesticated_sha256": hashlib.sha256(by_gene[gene].encode()).hexdigest(),
                "reproduced": ok,
            }
        )
        print(f"  {gene:12s} {'ok' if ok else 'MISMATCH'}")

    record = {
        "checked_utc": datetime.now(UTC).isoformat(timespec="seconds"),
        "enzymes": enzymes,
        "seed": seed,
        "versions": {m: version(m) for m in ("dnachisel", "biopython", "numpy")},
        "genes": rows,
    }
    (out_dir / "reproducibility.json").write_text(json.dumps(record, indent=2) + "\n")
    print(
        f"\n{sum(r['reproduced'] for r in rows)}/{len(rows)} reproduced byte-for-byte "
        f"-> reproducibility.json"
    )
    for b in bad:
        print(f"  FAIL {b}")
    return 1 if bad else 0


def design(symbol: str, enzymes: list[str], seed: int) -> dict:
    resolved = gene_resolver.resolve([symbol])[symbol]
    if resolved["notfound"] or not resolved["symbol"]:
        raise RuntimeError(f"{symbol}: gene_resolver could not canonicalize")
    approved = resolved["symbol"]

    gene = seqtools.resolve_gene(approved)
    if gene.ensembl_gene_id not in resolved["ensembl"]:
        raise RuntimeError(
            f"{approved}: Ensembl canonical gene {gene.ensembl_gene_id} is not among "
            f"the resolver's ids {resolved['ensembl']}"
        )

    native = seqtools.fetch_cds(gene.ensembl_transcript_id)
    protein = seqtools.translate(native).rstrip("*")

    reference = uniprot_by_ensg(gene.ensembl_gene_id)
    if protein != reference["sequence"]:
        raise RuntimeError(
            f"{approved}: Ensembl canonical protein ({len(protein)} aa) does not match "
            f"UniProt {reference['accession']} ({reference['length']} aa)"
        )

    before = seqtools.find_sites(native, enzymes)
    domesticated, warnings = seqtools.domesticate(
        native, enzymes, preserve_translation=True, seed=seed
    )
    if seqtools.translate(domesticated).rstrip("*") != protein:
        raise RuntimeError(f"{approved}: BUG — domestication changed the protein")
    after = seqtools.find_sites(domesticated, enzymes)
    if after:
        raise RuntimeError(f"{approved}: {len(after)} site(s) survive domestication: {after}")

    edits = edit_positions(native, domesticated)
    return {
        "gene": approved,
        "queried_as": symbol,
        "uniprot": reference["accession"],
        "ensembl_gene": gene.ensembl_gene_id_versioned,
        "ensembl_transcript": gene.ensembl_transcript_id,
        "ensembl_protein": gene.ensembl_protein_id,
        "aa": len(protein),
        "cds_nt": len(native),
        "signal_peptide": reference["signal_peptide"] or "none",
        "n_terminus": reference["n_terminus"] or "?",
        "c_terminus": reference["gpi"] or reference["c_terminus"] or "?",
        "sites_before": ";".join(f"{e}@{p}{s}" for e, s, p in before) or "none",
        "sites_before_n": len(before),
        "per_enzyme_before": {e: sum(1 for hit in before if hit[0] == e) for e in enzymes},
        "sites_after_n": len(after),
        "edits": len(edits),
        "edit_positions": ";".join(str(p) for p in edits) or "",
        "gc_before": round(seqtools.gc_content(native), 4),
        "gc_after": round(seqtools.gc_content(domesticated), 4),
        "warnings": " | ".join(warnings) or "none",
        "native_cds": native,
        "domesticated_cds": domesticated,
        "protein": protein,
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--enzymes", nargs="+", default=DEFAULT_ENZYMES)
    parser.add_argument("--genes", nargs="+", default=PANEL)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--out-dir", type=Path, default=Path(__file__).parent)
    parser.add_argument(
        "--on-error",
        choices=["stop", "skip"],
        default="stop",
        help="stop (default) halts on the first gene that fails a check, which is what a "
        "hand-checked panel wants. skip records the gene in skipped_orfs.tsv and "
        "continues, which is what a library-scale run wants: an Ensembl-canonical / "
        "UniProt-canonical disagreement is a per-gene scientific call about which "
        "sequence is authoritative, and halting the batch on the first one hides "
        "how many others there are. A skipped gene is never designed from the "
        "mismatched sequence.",
    )
    parser.add_argument(
        "--cache-dir",
        default=None,
        help="where per-gene records are cached for resume (default: <out-dir>/.orf_cache).",
    )
    parser.add_argument(
        "--verify",
        action="store_true",
        help="re-domesticate every CDS in the out-dir's orf_native.fasta and require a "
        "byte-identical match to orf_domesticated.fasta. Offline; exits non-zero on "
        "any drift.",
    )
    args = parser.parse_args(argv)

    if args.verify:
        return verify(args.out_dir, args.enzymes, args.seed)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    # Per-gene cache, written the moment a gene completes. Each record costs a MyGene
    # call, an Ensembl CDS fetch, a UniProt fetch and a constraint solve, and the batch
    # is otherwise only written at the end -- so a network timeout at gene 49 of 145
    # discarded 49 genes of work and the rerun paid for all of them again. A cached gene
    # is reused as-is; delete the directory to force a refetch.
    cache = Path(args.cache_dir) if args.cache_dir else args.out_dir / ".orf_cache"
    cache.mkdir(parents=True, exist_ok=True)

    results, skipped = [], []
    for symbol in args.genes:
        started = time.time()
        hit = cache / f"{symbol}.json"
        if hit.exists():
            record = json.loads(hit.read_text())
            results.append(record)
            print(f"{record['gene']:8s} cached", flush=True)
            continue
        try:
            record = design(symbol, args.enzymes, args.seed)
        except Exception as exc:
            # Deliberately broad: the first batch run died on a socket TimeoutError, which
            # a (RuntimeError, ValueError) clause does not catch, so one flaky fetch took
            # down the whole run. Every reason is recorded with its exception type.
            if args.on_error == "stop":
                raise
            skipped.append({"gene": symbol, "reason": f"{type(exc).__name__}: {exc}"})
            print(f"{symbol:8s} SKIPPED — {type(exc).__name__}: {exc}", flush=True)
            continue
        record["seconds"] = round(time.time() - started, 1)
        hit.write_text(json.dumps(record))
        per = record["per_enzyme_before"]
        print(
            f"{record['gene']:8s} {record['uniprot']:7s} {record['aa']:5d} aa  "
            f"{record['cds_nt']:5d} nt  sites "
            + " ".join(f"{e}={per[e]}" for e in args.enzymes)
            + f"  total={record['sites_before_n']:2d} -> 0  "
            f"edits={record['edits']:2d}  ({record['seconds']}s)",
            flush=True,
        )
        results.append(record)

    columns = [
        "gene",
        "uniprot",
        "ensembl_gene",
        "ensembl_transcript",
        "ensembl_protein",
        "aa",
        "cds_nt",
        "signal_peptide",
        "n_terminus",
        "c_terminus",
        "sites_before_n",
        "sites_before",
        "sites_after_n",
        "edits",
        "edit_positions",
        "gc_before",
        "gc_after",
        "warnings",
    ]
    if skipped:
        skip_table = args.out_dir / "skipped_orfs.tsv"
        with skip_table.open("w", newline="") as handle:
            w = csv.DictWriter(
                handle, fieldnames=["gene", "reason"], delimiter="\t", lineterminator="\n"
            )
            w.writeheader()
            w.writerows(skipped)
        print(f"\n{len(skipped)} gene(s) skipped -> {skip_table.name}")

    table = args.out_dir / "orf_design.tsv"
    with table.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=columns, delimiter="\t", extrasaction="ignore", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(results)

    enzymes = "+".join(args.enzymes)
    fasta = args.out_dir / "orf_domesticated.fasta"
    fasta.write_text(
        "".join(
            f">{r['gene']}|{r['uniprot']}|{r['ensembl_transcript']}|"
            f"domesticated:{enzymes}|{r['edits']}edits\n"
            + "\n".join(
                r["domesticated_cds"][i : i + 60] for i in range(0, len(r["domesticated_cds"]), 60)
            )
            + "\n"
            for r in results
        )
    )
    native = args.out_dir / "orf_native.fasta"
    native.write_text(
        "".join(
            f">{r['gene']}|{r['uniprot']}|{r['ensembl_transcript']}|native\n"
            + "\n".join(r["native_cds"][i : i + 60] for i in range(0, len(r["native_cds"]), 60))
            + "\n"
            for r in results
        )
    )

    print(f"\nEnzymes scrubbed: {enzymes}")
    print(f"Wrote {table}, {fasta}, {native}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
