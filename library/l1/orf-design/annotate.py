"""Per-gene annotation: the single source of topology for construct design.

Everything a placement decision needs, fetched rather than typed. Built to run over the
whole surfaceome, not a hand-listed panel — so nothing here is keyed to a gene name.

Per gene it resolves the symbol through ``deliverome_analysis.surfaceome.gene_resolver``, takes the
Ensembl canonical transcript, and reads the reviewed UniProt entry reached from the same
stable id. It never joins on a display symbol.

**GPI proteins are the case that forces this to be careful.** UniProt annotates no
``Topological domain`` features for a GPI-anchored protein, so a naive read of
``features[0]`` returns nothing and a placement rule consuming it would mis-place the
tag. FOLR1 is the worked example, and there are ~156 GPI proteins in the catalog. When a
GPI anchor is annotated, the mature N-terminus (after signal-peptide cleavage) is
extracellular and the C-terminus is **removed** at the omega site, so a C-terminal tag is
physically cleaved off. Both facts are recorded explicitly rather than inferred later.
"""

from __future__ import annotations

import argparse
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import asdict, dataclass
from pathlib import Path

from deliverome_analysis.cloning import seqtools
from deliverome_analysis.surfaceome import gene_resolver

HERE = Path(__file__).parent
UNIPROT = "https://rest.uniprot.org/uniprotkb/search"


@dataclass
class Annotation:
    """Everything the placement rule is allowed to look at."""

    gene: str
    queried_as: str
    uniprot: str
    ensembl_gene: str
    ensembl_transcript: str
    ensembl_protein: str
    aa: int
    cds_nt: int
    signal_peptide_end: int  # 0 when there is none
    n_term_extracellular: bool
    c_term_extracellular: bool
    gpi_anchored: bool
    gpi_omega_site: int | None
    n_term_ecto_len: int  # residues of extracellular N-terminus before the first TM
    c_term_ecto_len: int
    n_tm: int
    notes: str


def _retry(fn, *, tries: int = 4, base: float = 2.0, what: str = ""):
    """Retry a network call with backoff. Surfaceome-scale runs make thousands of
    requests; a transient timeout must not abort the batch."""
    for attempt in range(1, tries + 1):
        try:
            return fn()
        except (TimeoutError, urllib.error.URLError, ConnectionError) as err:
            if attempt == tries:
                raise RuntimeError(f"{what}: gave up after {tries} attempts ({err})") from err
            time.sleep(base**attempt)
    raise AssertionError("unreachable")


def _uniprot(ensg: str) -> dict:
    query = f"(xref:ensembl-{ensg}) AND (reviewed:true) AND (organism_id:9606)"
    url = (
        UNIPROT
        + "?"
        + urllib.parse.urlencode(
            {
                "query": query,
                "format": "json",
                "size": "1",
                "fields": "accession,length,sequence,ft_signal,ft_topo_dom,ft_transmem,ft_lipid",
            }
        )
    )

    def once():
        request = urllib.request.Request(url, headers={"User-Agent": seqtools.USER_AGENT})
        with urllib.request.urlopen(request, timeout=30) as response:
            return json.load(response).get("results", [])

    hits = _retry(once, what=f"UniProt {ensg}")
    if not hits:
        raise RuntimeError(f"no reviewed UniProt entry cross-referenced to {ensg}")
    return hits[0]


def annotate(symbol: str) -> Annotation:
    resolved = gene_resolver.resolve([symbol])[symbol]
    if resolved["notfound"] or not resolved["symbol"]:
        raise RuntimeError(f"{symbol}: gene_resolver could not canonicalize")
    approved = resolved["symbol"]
    gene = _retry(lambda: seqtools.resolve_gene(approved), what=f"Ensembl {approved}")
    if gene.ensembl_gene_id not in resolved["ensembl"]:
        raise RuntimeError(
            f"{approved}: canonical gene {gene.ensembl_gene_id} not among {resolved['ensembl']}"
        )

    cds = _retry(
        lambda: seqtools.fetch_cds(gene.ensembl_transcript_id),
        what=f"Ensembl CDS {gene.ensembl_transcript_id}",
    )
    protein = seqtools.translate(cds).rstrip("*")
    entry = _uniprot(gene.ensembl_gene_id)
    if protein != entry["sequence"]["value"]:
        raise RuntimeError(
            f"{approved}: Ensembl canonical protein ({len(protein)} aa) != UniProt "
            f"{entry['primaryAccession']} ({entry['sequence']['length']} aa)"
        )

    features = entry.get("features", [])

    def spans(kind):
        return [f for f in features if f["type"] == kind]

    signal = spans("Signal")
    topo = spans("Topological domain")
    tms = spans("Transmembrane")
    gpi = [f for f in spans("Lipidation") if "GPI" in f.get("description", "")]

    sp_end = signal[0]["location"]["end"]["value"] if signal else 0
    omega = gpi[0]["location"]["start"]["value"] if gpi else None
    notes = []

    if gpi:
        # No topological domains are annotated; the anchor defines the topology.
        n_out, c_out = True, False
        n_ecto = (omega - 1) - sp_end if omega else len(protein) - sp_end
        c_ecto = 0
        notes.append(f"GPI-anchored at {omega}; C-terminus cleaved at the omega site")
    elif topo:
        first, last = topo[0], topo[-1]
        n_out = first["description"].startswith("Extra")
        c_out = last["description"].startswith("Extra") and last["location"]["end"]["value"] == len(
            protein
        )
        n_ecto = first["location"]["end"]["value"] - sp_end if n_out else 0
        c_ecto = (
            last["location"]["end"]["value"] - last["location"]["start"]["value"] + 1
            if c_out
            else 0
        )
    else:
        raise RuntimeError(
            f"{approved}: UniProt annotates neither topological domains nor a GPI anchor "
            f"— topology cannot be determined, refusing to guess"
        )

    if n_out and not signal and n_ecto and n_ecto < 10:
        notes.append(
            f"no signal peptide and only {n_ecto} extracellular residues before TM1, which "
            f"is therefore the signal-anchor — a large N-terminal insert risks blocking "
            f"translocation (positive-inside rule)"
        )

    return Annotation(
        gene=approved,
        queried_as=symbol,
        uniprot=entry["primaryAccession"],
        ensembl_gene=gene.ensembl_gene_id_versioned,
        ensembl_transcript=gene.ensembl_transcript_id,
        ensembl_protein=gene.ensembl_protein_id,
        aa=len(protein),
        cds_nt=len(cds),
        signal_peptide_end=sp_end,
        n_term_extracellular=n_out,
        c_term_extracellular=c_out,
        gpi_anchored=bool(gpi),
        gpi_omega_site=omega,
        n_term_ecto_len=max(n_ecto, 0),
        c_term_ecto_len=max(c_ecto, 0),
        n_tm=len(tms),
        notes="; ".join(notes),
    )


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("genes", nargs="*", help="symbols; default reads genes.tsv")
    ap.add_argument("--out", type=Path, default=HERE / "annotations.json")
    args = ap.parse_args(argv)

    symbols = args.genes or [
        line.split("\t")[0].strip()
        for line in (HERE / "genes.tsv").read_text().splitlines()[1:]
        if line.strip()
    ]

    out = {}
    for symbol in symbols:
        record = annotate(symbol)
        out[record.gene] = asdict(record)
        term = (
            "C-term ecto"
            if record.c_term_extracellular
            else "N-term ecto"
            if record.n_term_extracellular
            else "no ecto terminus"
        )
        print(
            f"{record.gene:9s} {record.uniprot:7s} {record.aa:5d} aa  SP={record.signal_peptide_end:3d}  "
            f"TM={record.n_tm:2d}  {'GPI ' if record.gpi_anchored else '    '}{term:16s}"
            f"{'  ⚠ ' + record.notes if record.notes else ''}"
        )
    args.out.write_text(json.dumps(out, indent=2) + "\n")
    print(f"\nWrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
