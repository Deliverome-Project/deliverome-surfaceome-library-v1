"""Resolve every candidate construct to a coding sequence that is verified to
encode the form the library actually selected.

Candidates are the library's human references, the full isoform cohort, both
full ortholog sets, and the expression-rescue arms.

The trap this guards against: a gene's MANE or APPRIS transcript is not always
the form the library picked. Six GPI references were named outright and others
come off the UniProt-canonical rung, so their gene's canonical transcript
encodes a different protein. The test is therefore translation equality against
the selected protein wherever that protein is known, not a start-codon or
length heuristic -- Ensembl annotates some CDS entries on a non-ATG start, and
rejecting those loses a sequence that is in fact correct.

Sequences come from the pinned Ensembl 115 FTP release, archived on S3, rather
than from thousands of live REST calls.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

import candidate_inputs as inputs

ENSEMBL = inputs.ASSEMBLIES
# How far an Ensembl translation may sit from the UniProt sequence of the same
# length and still be accepted. These are reference-allele differences between
# GRCh38 and UniProt, not different isoforms, and refusing them throws away a
# transcript that is otherwise exactly right. Anything accepted this way is
# recorded with its distance so the call stays visible.
VARIANT_TOLERANCE = 3
#: Appended when an Ensembl CDS record ends on a sense codon. Ensembl leaves the stop out
#: of some records (four KIR transcripts on alternate-haplotype scaffolds and nine cyno
#: models in release 115) although the protein is complete and equals UniProt. An ORF
#: without a stop would read through into the vector, so one is supplied: TAA, the stop
#: the library's C-terminal constructs already use. The native codon cannot be read from
#: the pinned inputs (the primary-assembly GTF does not cover those scaffolds, and no
#: genome sequence is pinned). Entries that get one carry ``stop_added``.
ADDED_STOP = "TAA"
STOPS = {"TAA", "TAG", "TGA"}
CODONS = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L",
    "CTG": "L", "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M", "GTT": "V", "GTC": "V",
    "GTA": "V", "GTG": "V", "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "CCT": "P",
    "CCC": "P", "CCA": "P", "CCG": "P", "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A", "TAT": "Y", "TAC": "Y", "CAT": "H",
    "CAC": "H", "CAA": "Q", "CAG": "Q", "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E", "TGT": "C", "TGC": "C", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGT": "S", "AGC": "S", "AGA": "R",
    "AGG": "R", "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    "TAA": "*", "TAG": "*", "TGA": "*",
}


def translate(cds: str) -> str:
    return "".join(CODONS.get(cds[i:i + 3], "X") for i in range(0, len(cds) - 2, 3))


def read_fasta(path: Path) -> dict[str, tuple[str, str]]:
    out, name, head, buf = {}, None, None, []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    out[name] = (head, "".join(buf))
                head = line[1:].strip()
                name = head.split()[0].split(".")[0]
                buf = []
            else:
                buf.append(line.strip())
    if name:
        out[name] = (head, "".join(buf))
    return out


def load_uniprot(species: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for path in inputs.uniprot_tables(species):
        with gzip.open(path, "rt") if path.suffix == ".gz" else path.open() as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                acc, seq = row.get("Entry"), row.get("Sequence")
                if acc and seq:
                    out[acc] = seq.strip().upper()
    return out


class Resolver:
    """Finds, for one species, a transcript whose CDS encodes a given protein."""

    def __init__(self, species: str):
        self.cds = read_fasta(inputs.cds_fasta(species))
        pep = read_fasta(inputs.pep_fasta(species))
        self.uniprot = load_uniprot(species)
        # every transcript a protein sequence maps to, not just the first: the
        # first may carry a CDS record that cannot be used.
        self.by_protein: dict[str, list[str]] = defaultdict(list)
        for _ensp, (head, seq) in pep.items():
            m = re.search(r"transcript:(\S+)", head)
            if m:
                self.by_protein[seq.rstrip("*").upper()].append(m.group(1).split(".")[0])
        self.by_symbol: dict[str, list[str]] = defaultdict(list)
        for enst, (head, _seq) in self.cds.items():
            m = re.search(r"gene_symbol:(\S+)", head)
            if m:
                self.by_symbol[m.group(1).upper()].append(enst)

    def usable(self, enst: str, want_protein: str | None,
               want_len: int | None) -> tuple[str, int] | None:
        """(CDS, residues differing from the wanted protein) or None."""
        rec = self.cds.get(enst.split(".")[0])
        if rec is None:
            return None
        seq = rec[1].upper()
        if len(seq) % 3:
            return None
        if set(seq) - set("ACGT"):
            # An ambiguity code cannot be codon-optimised or ordered as DNA.
            return None
        prot = translate(seq).rstrip("*")
        if "*" in prot:                      # a premature stop is never right
            return None
        if seq[-3:] not in STOPS:
            seq += ADDED_STOP
        if want_protein is not None:
            if len(prot) != len(want_protein):
                return None
            d = sum(a != b for a, b in zip(prot, want_protein))
            return (seq, d) if d <= VARIANT_TOLERANCE else None
        if want_len is not None:
            return (seq, 0) if len(prot) == want_len else None
        return (seq, 0) if seq.startswith("ATG") else None

    def find(self, transcripts: list[str], acc: str | None, symbol: str | None,
             want_len: int | None) -> tuple[str, str, str, int] | None:
        """(transcript, cds, route, variant distance), best match first.

        Every route is tried and the closest hit wins, so an exact match is
        never passed over for a tolerated variant found earlier.
        """
        want_protein = self.uniprot.get(acc) if acc else None
        best: tuple[str, str, str, int] | None = None

        def offer(enst: str, hit, route: str):
            nonlocal best
            if hit and (best is None or hit[1] < best[3]):
                best = (enst.split(".")[0], hit[0], route, hit[1])

        for enst in transcripts:
            offer(enst, self.usable(enst, want_protein, want_len), "named transcript")
        if best and best[3] == 0:
            return best
        if want_protein:
            for enst in self.by_protein.get(want_protein, []):
                offer(enst, self.usable(enst, want_protein, None), "protein match")
            if best and best[3] == 0:
                return best
            # same length, a few residues apart: the reference-allele case
            for enst in self.by_symbol.get((symbol or "").upper(), []):
                offer(enst, self.usable(enst, want_protein, None), "allele variant")
            if best and best[3] == 0:
                return best
        if symbol:
            for enst in self.by_symbol.get(symbol.upper(), []):
                hit = self.usable(enst, None, want_len)
                if hit and want_protein:
                    # Matched on length, but a UniProt sequence exists, so report
                    # how far this actually is rather than implying a clean hit.
                    prot = translate(hit[0]).rstrip("*")
                    d = (sum(a != b for a, b in zip(prot, want_protein))
                         if len(prot) == len(want_protein) else -1)
                    hit = (hit[0], d if d >= 0 else 0)
                offer(enst, hit, "gene symbol")
        return best


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)

    ap.add_argument("--out", default="candidate_cds.json")
    ap.add_argument("--report", default="candidate_cds_coverage.tsv")
    args = ap.parse_args()

    res = {}
    for sp in ENSEMBL:
        res[sp] = Resolver(sp)
        print(f"  {sp:6} CDS {len(res[sp].cds):7}  proteins {len(res[sp].by_protein):7}  "
              f"symbols {len(res[sp].by_symbol):7}  uniprot {len(res[sp].uniprot):7}")

    out: dict[str, dict] = {}
    rows: list[dict] = []

    def take(key, sp, transcripts, acc, symbol, want_len, group):
        hit = res[sp].find(transcripts, acc, symbol, want_len)
        if hit:
            enst, seq, route, dist = hit
            out[key] = {"species": sp, "transcript": enst, "cds": seq,
                        "origin": group, "route": route, "variant_residues": dist}
            if res[sp].cds[enst][1].upper()[-3:] not in STOPS:
                out[key]["stop_added"] = ADDED_STOP
            rows.append({"key": key, "group": group, "transcript": enst,
                         "status": "ok", "route": route,
                         "detail": f"{len(seq)} nt" + (f"; {dist} residue(s) from the "
                                                      f"UniProt sequence" if dist else "")
                         + ("; Ensembl record has no stop, TAA appended"
                            if "stop_added" in out[key] else "")})
        else:
            rows.append({"key": key, "group": group, "transcript": ";".join(transcripts) or "-",
                         "status": "unresolved", "route": "-",
                         "detail": f"no transcript encodes {acc or symbol or 'the selected form'}"})

    def ints(v):
        return int(v) if v.strip().isdigit() else None

    for r in csv.DictReader(open("reference_selection.tsv"), delimiter="\t"):
        rung = r["reference_rung"]
        ts = [r["mane_enst"]] if rung == "mane" and r["mane_enst"].strip() else \
             [r["appris_enst"]] if rung == "appris" and r["appris_enst"].strip() else \
             [x for x in (r["mane_enst"], r["appris_enst"]) if x.strip()]
        take(f"lib:{r['symbol']}", "human", ts,
             r["selected_uniprot_equivalent"].strip() or None, r["symbol"],
             ints(r["selected_length"]), "library")

    for r in csv.DictReader(open("isoform_cohorts.tsv"), delimiter="\t"):
        ts = [x for x in r["ensembl_transcripts"].split(";") if x]
        # Ensembl-only forms carry no accession, so the accession alone would
        # collapse every such form of a gene onto one key and silently drop all
        # but the last. The first transcript discriminates them.
        tag = r["uniprot"].strip() or (ts[0] if ts else f"len{r['length']}")
        take(f"iso:{r['gene']}:{tag}", "human", ts,
             r["uniprot"].strip() or None, r["gene"], ints(r["length"]), "isoform")

    for sp, path in (("mouse", "ortholog_repick_mouse.tsv"), ("cyno", "ortholog_repick_cyno.tsv")):
        for r in csv.DictReader(open(path), delimiter="\t"):
            acc = r["picked_acc"].strip()
            take(f"ort:{sp}:{r['human_symbol']}:{acc}", sp, [acc] if acc.startswith("ENS") else [],
                 acc if not acc.startswith("ENS") else None,
                 r["ortholog_symbol"], ints(r["picked_length"]), f"ortholog_{sp}")

    if len(out) + sum(1 for r in rows if r["status"] != "ok") != len(rows):
        seen = Counter(r["key"] for r in rows)
        clash = {k: n for k, n in seen.items() if n > 1}
        raise SystemExit(f"key collision would drop candidates: {len(clash)} keys, e.g. "
                         f"{list(clash)[:3]}")

    Path(args.out).write_text(json.dumps(out, indent=0, sort_keys=True))
    with open(args.report, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["key", "group", "transcript", "status", "route", "detail"],
                           delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    by = Counter((r["group"], r["status"]) for r in rows)
    route = Counter(r["route"] for r in rows if r["status"] == "ok")
    variant = sum(1 for v in out.values() if v.get("variant_residues"))
    print(f"\n  {'group':16} {'ok':>6} {'unresolved':>11}")
    for g in sorted({k[0] for k in by}):
        print(f"  {g:16} {by[(g,'ok')]:6} {by[(g,'unresolved')]:11}")
    print(f"\n  routes used: {dict(route)}")
    print(f"  accepted with 1-{VARIANT_TOLERANCE} residue(s) differing: {variant}")
    print(f"  stop codon appended (Ensembl record ends on a sense codon): "
          f"{sorted(k for k, v in out.items() if 'stop_added' in v)}")
    print(f"  resolved {len(out)} of {len(rows)} ({len(out)/len(rows):.1%})")
    Path("candidate_inputs_provenance.json").write_text(
        json.dumps(inputs.provenance(), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
