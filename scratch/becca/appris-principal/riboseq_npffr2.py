"""Does the GENCODE Ribo-seq ORF catalogue settle NPFFR2's start codon? And what does
Ensembl itself say?

Ribosome profiling reports translation directly, so a Ribo-seq ORF catalogue is the right
instinct for a start-codon question that sequence annotation cannot settle. The GENCODE
catalogue (Phase 1, GENCODE v35; Phase 2 primary and comprehensive, GENCODE v45) is the
community consensus set.

Two things come out, and the second is the one that matters:

1. **NPFFR2 has no Ribo-seq ORF of any type**, in any of the three sets. The catalogue also
   cannot answer this class of question by construction: it catalogues ORFs "not currently
   represented as protein-coding CDSs in reference gene annotation", and an in-frame
   N-terminal extension of a CDS is a change *to* a CDS, so it is out of scope. The type
   vocabulary is checked here rather than assumed.
2. **No Ensembl 115 transcript encodes our 522-residue reference at all.** Every translation
   at ENSG00000056291 is 420, 423 or 304 residues, and the MANE Select transcript, carrying
   CCDS3552.1, encodes the 420-residue form that both Tedman and APPRIS chose.

    uv run --locked python riboseq_npffr2.py
"""

from __future__ import annotations

import csv
import gzip
import json
import re
from collections import Counter
from pathlib import Path

import requests

ACCESSION = "Q9Y5X5"
GENE = "NPFFR2"
ALIASES = ["NPFFR2", "NPFF2", "GPR74", "NPGPR"]
ENSG = "ENSG00000056291"

RIBOSEQ_BASE = "https://ftp.ebi.ac.uk/pub/databases/gencode/riboseq_orfs"
RIBOSEQ_SETS = {
    "phase1_all": ("phase1/Ribo-seq_ORFs.all.bed", "v35"),
    "phase2_primary": ("phase2/Ribo-seq_ORFs.primary.bed", "v45"),
    "phase2_comprehensive": ("phase2/Ribo-seq_ORFs.comprehensive.bed", "v45"),
}


def fetch(cache: Path) -> dict[str, Path]:
    cache.mkdir(parents=True, exist_ok=True)
    out = {}
    for name, (rel, _) in RIBOSEQ_SETS.items():
        dest = cache / Path(rel).name
        if not dest.exists():
            resp = requests.get(f"{RIBOSEQ_BASE}/{rel}", timeout=600)
            resp.raise_for_status()
            dest.write_bytes(resp.content)
            print(f"  fetched {dest.name}")
        out[name] = dest
    return out


def main() -> None:
    cache = Path("~/.cache/deliverome/riboseq").expanduser()
    appris_cache = Path("~/.cache/deliverome/appris").expanduser()
    beds = fetch(cache)

    result: dict = {"gene": GENE, "accession": ACCESSION, "ensembl_gene": ENSG}

    # ---- 1. is the gene in any Ribo-seq ORF set, under any alias? ----------------
    for name, path in beds.items():
        names = [line.split("\t")[3] for line in path.read_text().splitlines() if line]
        # The name field is <orf_id>:<gene_name>:<last>. Per the Phase 1 README that last
        # field is the *transcript biotype*; in the Phase 2 files it is the *ORF type*.
        # Only Phase 2 therefore carries the ORF-type vocabulary.
        genes = {n.split(":")[1] for n in names if n.count(":") >= 2}
        types = Counter(n.split(":")[-1] for n in names)
        hits = {a: sum(1 for n in names if f":{a}:" in n) for a in ALIASES}
        result[name] = {
            "n_orfs": len(names),
            "n_genes": len(genes),
            "orf_types" if name.startswith("phase2") else "transcript_biotypes": dict(
                types.most_common()
            ),
            "hits_for_aliases": hits,
        }
        print(f"{name}: {len(names)} ORFs over {len(genes)} genes; {GENE} hits {sum(hits.values())}")
        print(
            f"   {'ORF types' if name.startswith('phase2') else 'transcript biotypes'}: "
            f"{dict(types.most_common())}"
        )

    assert all(
        sum(result[n]["hits_for_aliases"].values()) == 0 for n in RIBOSEQ_SETS
    ), "an alias now has a Ribo-seq ORF; the conclusion below needs revisiting"

    # No type in the vocabulary is an in-frame N-terminal extension of a CDS.
    all_types = set()
    for n in RIBOSEQ_SETS:
        if n.startswith("phase2"):
            all_types |= set(result[n]["orf_types"])
    result["orf_type_vocabulary"] = sorted(all_types)
    print(f"\nPhase 2 ORF-type vocabulary: {sorted(all_types)}")

    # ---- 2. what does Ensembl itself encode at this locus? -----------------------
    uniprot = {}
    with gzip.open(appris_cache / "uniprot_human_reviewed_isoforms.tsv.gz", "rt") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["Entry"].split("-")[0] == ACCESSION:
                uniprot[row["Entry"]] = row["Sequence"]

    translations: list[tuple[str, str]] = []
    cur, buf, hdr = None, [], ""
    with gzip.open(appris_cache / "Homo_sapiens.GRCh38.pep.all.fa.gz", "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if cur and ENSG in hdr:
                    translations.append((cur, "".join(buf).rstrip("*")))
                hdr = line
                m = re.search(r"\btranscript:(\S+)", line)
                cur = m.group(1).split(".")[0] if m else None
                buf = []
            elif cur:
                buf.append(line.strip())
    if cur and ENSG in hdr:
        translations.append((cur, "".join(buf).rstrip("*")))

    mane, ccds = {}, {}
    with (appris_cache / "appris_homo_sapiens_e115v50.txt").open() as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) > 5 and f[1].split(".")[0] == ENSG:
                enst = f[2].split(".")[0]
                ccds[enst] = None if f[3] in ("-", "") else f[3]
                mane[enst] = "MANE_Select" in f[5]

    canonical = uniprot[ACCESSION]
    print(f"\nUniProt proteoforms: {({k: len(v) for k, v in uniprot.items()})}")
    print(f"\n{'transcript':18s} {'len':>5s} {'CCDS':13s} {'MANE':6s}  UniProt equivalent")
    rows = []
    for enst, seq in sorted(translations, key=lambda x: -len(x[1])):
        eq = next((k for k, v in uniprot.items() if v == seq), "")
        rows.append(
            {
                "transcript": enst,
                "length": len(seq),
                "ccds": ccds.get(enst) or "",
                "mane_select": bool(mane.get(enst)),
                "uniprot_equivalent": eq,
            }
        )
        print(
            f"{enst:18s} {len(seq):>5d} {ccds.get(enst) or '-':13s} "
            f"{'MANE' if mane.get(enst) else '':6s}  {eq or '(none)'}"
        )

    encodes_canonical = any(seq == canonical for _, seq in translations)
    short = next((s for _, s in translations if len(s) == 420), None)
    result |= {
        "uniprot_proteoform_lengths": {k: len(v) for k, v in uniprot.items()},
        "ensembl_115_translations": rows,
        "ensembl_encodes_uniprot_canonical": encodes_canonical,
        "short_form_is_suffix_of_canonical": bool(short and canonical.endswith(short)),
        "n_terminal_extension_length": len(canonical) - len(short) if short else None,
    }
    print(f"\nany Ensembl 115 transcript encoding the {len(canonical)}-aa UniProt canonical: "
          f"{encodes_canonical}")
    print(f"the 420-aa form is an exact suffix of ours: {result['short_form_is_suffix_of_canonical']} "
          f"(difference is exactly the first {result['n_terminal_extension_length']} residues)")

    Path("riboseq_npffr2.json").write_text(json.dumps(result, indent=2) + "\n")
    print("wrote riboseq_npffr2.json")


if __name__ == "__main__":
    main()
