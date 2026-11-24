"""Select a reference isoform per library gene with APPRIS, and audit it against the
UniProt sequences we actually store.

APPRIS annotates every Ensembl transcript with a principal-isoform class. PRINCIPAL:1 is
the strongest: a single CCDS-supported transcript whose structural and functional evidence
agrees. Weaker classes (PRINCIPAL:2..5, PRINCIPAL:M) are assigned when that evidence does
not single one transcript out. APPRIS makes this call *within each species independently*,
which is what makes the cross-species check at the end necessary.

Three things come out of this:

1. **Selection.** The APPRIS principal transcript per library gene, for human and mouse,
   with the class it was assigned and whether it is also MANE Select.
2. **Agreement with what we store.** The principal transcript's protein sequence compared
   to the UniProt canonical we currently build ORFs from, and to every UniProt isoform we
   hold for that gene. A principal that matches neither is reported in full: it means
   APPRIS is pointing at a proteoform our pipeline has never seen.
3. **Human/mouse pair consistency.** APPRIS picks per species, so a human principal and a
   mouse principal need not be the same splice form. Where they differ substantially in
   length the ortholog construct is not the counterpart of the human one, and a measured
   cross-species difference could be splice form rather than species.

Cynomolgus is absent by construction: APPRIS e115v50 covers Macaca mulatta, our orthologs
are Macaca fascicularis, and the two macaques differ at ~9.7M SNVs, so rhesus is not used
as a stand-in. The cyno cohort is written out as an explicit gap list.

    uv run --locked python select_principal.py \
        --surfaceome-repo ~/Git/accessible-surfaceome \
        --cache-dir ~/.cache/deliverome/appris
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
import statistics
from collections import defaultdict
from pathlib import Path

import requests

APPRIS_RELEASE = "e115v50"
ENSEMBL_RELEASE = "115"
UNIPROT_RELEASE = "2026_03"
UNIPROT_HUMAN = "uniprot_human_reviewed_isoforms.tsv.gz"
UNIPROT_MOUSE = "uniprot_mouse_reviewed_isoforms.tsv.gz"

# Best (most confident) first. APPRIS assigns one of these per transcript; a gene's
# selection uses the best class present for that gene, never a mix of classes.
PRINCIPAL_ORDER = [
    "PRINCIPAL:1",
    "PRINCIPAL:2",
    "PRINCIPAL:3",
    "PRINCIPAL:4",
    "PRINCIPAL:5",
    "PRINCIPAL:M",
]
PRINCIPAL_RANK = {lab: i for i, lab in enumerate(PRINCIPAL_ORDER)}

# Relative length difference above which a human/mouse principal pair is called
# inconsistent. Reported alongside the full distribution so it can be re-cut without
# rerunning the sequence work.
PAIR_LENGTH_TOLERANCE = 0.10


# --------------------------------------------------------------------------- parsing


def read_fasta(path: Path) -> dict[str, str]:
    """Return {header: sequence}. Handles plain and gzipped files."""
    opener = gzip.open if path.suffix == ".gz" else open
    out: dict[str, str] = {}
    header: str | None = None
    chunks: list[str] = []
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if header is not None:
                    out[header] = "".join(chunks)
                header = line[1:].rstrip("\n")
                chunks = []
            else:
                chunks.append(line.strip())
    if header is not None:
        out[header] = "".join(chunks)
    return out


def parse_uniprot_fasta(path: Path) -> dict[str, dict]:
    """{accession: {symbol, seq, is_isoform}} from a UniProt-style FASTA."""
    out = {}
    for header, seq in read_fasta(path).items():
        parts = header.split("|")
        if len(parts) < 3:
            continue
        acc = parts[1]
        gn = re.search(r"\bGN=(\S+)", header)
        out[acc] = {
            "accession": acc,
            "symbol": gn.group(1) if gn else None,
            "seq": seq,
            "is_isoform": "-" in acc,
        }
    return out


def parse_ensembl_pep(path: Path) -> dict[str, str]:
    """{unversioned ENST: protein sequence} from an Ensembl pep.all FASTA."""
    out = {}
    for header, seq in read_fasta(path).items():
        m = re.search(r"\btranscript:(\S+)", header)
        if not m:
            continue
        enst = m.group(1).split(".")[0]
        out[enst] = seq.rstrip("*")
    return out


def parse_appris(path: Path) -> dict[str, list[dict]]:
    """{ENSG: [row, ...]}.

    The human table carries a header row and a sixth MANE column; the mouse table has
    neither. Both are handled rather than assumed.
    """
    by_gene: dict[str, list[dict]] = defaultdict(list)
    with path.open() as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 5 or f[0] == "Gene name":
                continue
            by_gene[f[1].split(".")[0]].append(
                {
                    "symbol": f[0],
                    "ensg": f[1].split(".")[0],
                    "enst": f[2].split(".")[0],
                    "ccds": None if f[3] in ("-", "") else f[3],
                    "label": f[4],
                    "mane": "MANE_Select" in f[5] if len(f) > 5 else False,
                }
            )
    return by_gene


# --------------------------------------------------------------------------- selection


def select_principal(rows: list[dict], pep: dict[str, str]) -> dict | None:
    """Pick one principal transcript for a gene.

    Takes the best principal class present, then collapses the transcripts in that class
    by *protein sequence* — several transcripts commonly share one protein and differ only
    in UTR, which is not an ambiguity. A genuine ambiguity (several distinct proteins in
    the same class) is broken by MANE Select, then by length, and is flagged either way.
    """
    principal = [r for r in rows if r["label"] in PRINCIPAL_RANK]
    if not principal:
        return None
    best = min(PRINCIPAL_RANK[r["label"]] for r in principal)
    tier = [r for r in principal if PRINCIPAL_RANK[r["label"]] == best]

    by_seq: dict[str, list[dict]] = defaultdict(list)
    for r in tier:
        seq = pep.get(r["enst"])
        if seq:
            by_seq[seq].append(r)
    if not by_seq:
        return {
            "label": PRINCIPAL_ORDER[best],
            "enst": ";".join(sorted(r["enst"] for r in tier)),
            "seq": None,
            "n_distinct_proteins": 0,
            "ambiguous": False,
            "mane": any(r["mane"] for r in tier),
            "ccds": None,
            "note": "no Ensembl 115 peptide for any transcript in the principal class",
        }

    if len(by_seq) == 1:
        seq = next(iter(by_seq))
        chosen = by_seq[seq]
    else:
        mane_seqs = [s for s, rs in by_seq.items() if any(r["mane"] for r in rs)]
        seq = mane_seqs[0] if len(mane_seqs) == 1 else max(by_seq, key=len)
        chosen = by_seq[seq]

    return {
        "label": PRINCIPAL_ORDER[best],
        "enst": ";".join(sorted(r["enst"] for r in chosen)),
        "seq": seq,
        "n_distinct_proteins": len(by_seq),
        "ambiguous": len(by_seq) > 1,
        "mane": any(r["mane"] for r in chosen),
        "ccds": next((r["ccds"] for r in chosen if r["ccds"]), None),
        "note": "",
    }


def classify_difference(seq: str, ref: str) -> dict:
    """Describe how a principal differs from the reference, in terms that matter here.

    Which end a difference falls on decides whether it changes a design: an N-terminal
    difference moves the signal peptide and every post-SP placement derived from it, a
    C-terminal one moves the GPI omega site and the C-terminal tag arm, and a pure
    substitution moves neither.
    """
    out = {
        "len_principal": len(seq),
        "len_reference": len(ref),
        "len_delta": len(seq) - len(ref),
    }
    if len(seq) == len(ref):
        out["difference_kind"] = "substitution_only"
        changed = [(i + 1, a, b) for i, (a, b) in enumerate(zip(ref, seq)) if a != b]
        out["substitutions_vs_canonical"] = len(changed)
        out["_substitutions"] = changed
        return out
    short, long_ = (seq, ref) if len(seq) < len(ref) else (ref, seq)
    if long_.endswith(short):
        out["difference_kind"] = "n_terminal"
    elif long_.startswith(short):
        out["difference_kind"] = "c_terminal"
    else:
        out["difference_kind"] = "internal"
    return out


def compare_to_stored(
    seq: str,
    canonical: dict | None,
    isoforms: list[dict],
    uniprot_proteoforms: list[dict],
) -> dict:
    """Classify a principal protein against the UniProt sequences available for the gene.

    Four outcomes, in decreasing comfort: it is the canonical we already build from; it is
    an isoform we already hold; it is a UniProt proteoform that exists but is not in our
    cohort; or it is a sequence no reviewed UniProt entry for that gene carries, which
    means APPRIS and UniProt disagree about what the protein is.
    """
    if canonical and seq == canonical["seq"]:
        return {"match": "canonical", "matched_accession": canonical["accession"]}
    for iso in isoforms:
        if seq == iso["seq"]:
            return {"match": "stored_isoform", "matched_accession": iso["accession"]}
    for pf in uniprot_proteoforms:
        if seq == pf["seq"]:
            return {
                "match": "uniprot_proteoform_not_stored",
                "matched_accession": pf["accession"],
            }

    detail = {"match": "no_uniprot_match", "matched_accession": ""}
    if canonical:
        detail |= classify_difference(seq, canonical["seq"])
    return detail


# --------------------------------------------------------------------------- loading


def annotate_allelic_variants(rows: list[dict], cache: Path) -> None:
    """Mark substitution-only differences that are documented UniProt natural variants.

    A principal that differs from the UniProt sequence at a handful of positions is usually
    not a different protein: UniProt's reference entry and the GRCh38/GRCm39 reference
    genome carry different alleles at a common polymorphism. UniProt annotates those
    positions as natural variants, so the question is answerable rather than a judgement
    call -- and the answer decides whether the difference touches the design (it does not)
    or signals a genuinely different proteoform (it does).
    """
    store = cache / "uniprot_variants.json"
    seen: dict[str, dict[str, list[str]]] = (
        json.loads(store.read_text()) if store.exists() else {}
    )
    todo = sorted(
        {
            r["accession"]
            for r in rows
            if r.get("difference_kind") == "substitution_only"
            and r["accession"] not in seen
        }
    )
    for n, acc in enumerate(todo, 1):
        try:
            resp = requests.get(
                f"https://rest.uniprot.org/uniprotkb/{acc}.json",
                params={"fields": "ft_variant"},
                timeout=60,
            )
            resp.raise_for_status()
            feats = resp.json().get("features", [])
        except Exception as exc:  # network trouble must not silently become "not a variant"
            seen[acc] = {"__error__": [str(exc)]}
            continue
        at: dict[str, list[str]] = {}
        for f in feats:
            if f.get("type") != "Natural variant":
                continue
            pos = str(f["location"]["start"]["value"])
            alts = f.get("alternativeSequence", {}).get("alternativeSequences", [])
            at.setdefault(pos, []).extend(alts)
        seen[acc] = at
        if n % 25 == 0:
            print(f"  UniProt variants: {n}/{len(todo)}")
    store.write_text(json.dumps(seen, indent=0, sort_keys=True))

    for r in rows:
        if r.get("difference_kind") != "substitution_only":
            continue
        at = seen.get(r["accession"], {})
        if "__error__" in at:
            r["substitutions_are_annotated_variants"] = "unknown"
            continue
        changed = r.get("_substitutions") or []
        annotated = sum(1 for pos, _ref, alt in changed if alt in at.get(str(pos), []))
        at_variant_position = sum(1 for pos, _ref, _alt in changed if str(pos) in at)
        r["substitutions_are_annotated_variants"] = (
            "all" if annotated == len(changed) and changed
            else "positions_only" if at_variant_position == len(changed) and changed
            else "partial" if at_variant_position
            else "none"
        )


def load_uniprot_proteoforms(path: Path) -> dict[str, list[dict]]:
    """{primary gene symbol: [proteoform, ...]} from a reviewed-proteome TSV stream.

    The stream is requested with ``includeIsoform=true``, so each gene carries its
    canonical entry and every annotated isoform. The canonical is the entry whose
    accession has no ``-N`` suffix, and is placed first.
    """
    out: dict[str, list[dict]] = defaultdict(list)
    with gzip.open(path, "rt") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            sym = row["Gene Names (primary)"].strip()
            if not sym:
                continue
            acc = row["Entry"]
            rec = {
                "accession": acc,
                "symbol": sym,
                "seq": row["Sequence"],
                "is_isoform": "-" in acc,
            }
            if rec["is_isoform"]:
                out[sym].append(rec)
            else:
                out[sym].insert(0, rec)
    return dict(out)


def load_ortholog_metadata(path: Path) -> list[dict]:
    """Rows of a *_metadata.csv, which carry the human ENSG the ortholog was found from."""
    with path.open() as fh:
        return [r for r in csv.DictReader(fh) if r.get("status") == "ok"]


# --------------------------------------------------------------------------- main


def audit_gene(
    stored: dict,
    isoforms: list[dict],
    uniprot_proteoforms: list[dict],
    ensgs: list[str],
    appris: dict[str, list[dict]],
    pep: dict[str, str],
    by_symbol: dict[str, list[dict]] | None = None,
) -> dict:
    """Select the principal for one gene and compare it to the stored sequences."""
    rec = {
        "symbol": stored["symbol"],
        "accession": stored["accession"],
        "ensembl_gene": ";".join(ensgs),
        "len_stored": len(stored["seq"]),
    }
    rows = [r for g in ensgs for r in appris.get(g, [])]
    if not rows and by_symbol and stored["symbol"]:
        # Fall back to the APPRIS gene symbol only when no cross-reference exists, and
        # say so, so a symbol-based join is never mistaken for an ID-based one.
        rows = by_symbol.get(stored["symbol"], [])
        if rows:
            rec["ensembl_gene"] = ";".join(sorted({r["ensg"] for r in rows}))
            rec["ensg_source"] = "appris_symbol"
    rec.setdefault("ensg_source", "xref" if ensgs else "")

    if not rows:
        rec["status"] = "no_appris_entry"
        return rec

    sel = select_principal(rows, pep)
    if sel is None:
        rec["status"] = "no_principal_class"
        return rec

    rec |= {
        "appris_class": sel["label"],
        "appris_enst": sel["enst"],
        "appris_ccds": sel["ccds"] or "",
        "mane_select": sel["mane"],
        "n_distinct_principal_proteins": sel["n_distinct_proteins"],
        "ambiguous_principal": sel["ambiguous"],
    }
    if sel["seq"] is None:
        rec |= {"status": "no_peptide", "note": sel["note"]}
        return rec

    rec["len_principal"] = len(sel["seq"])
    cmp = compare_to_stored(sel["seq"], stored, isoforms, uniprot_proteoforms)
    rec |= cmp
    rec["status"] = cmp["match"]
    return rec


def write_tsv(path: Path, rows: list[dict], columns: list[str]) -> None:
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"wrote {path}  ({len(rows)} rows)")


COLUMNS = [
    "symbol",
    "accession",
    "sequence_source",
    "ensembl_gene",
    "ensg_source",
    "status",
    "appris_class",
    "appris_enst",
    "appris_ccds",
    "mane_select",
    "ambiguous_principal",
    "n_distinct_principal_proteins",
    "len_stored",
    "len_principal",
    "matched_accession",
    "difference_kind",
    "len_reference",
    "len_delta",
    "substitutions_vs_canonical",
    "substitutions_are_annotated_variants",
    "note",
]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--surfaceome-repo", default="~/Git/accessible-surfaceome")
    ap.add_argument("--cache-dir", default="~/.cache/deliverome/appris")
    ap.add_argument(
        "--library",
        default="../tag-site-architecture/library_members.json",
        help="JSON list of library gene symbols",
    )
    ap.add_argument("--out-dir", default=".")
    args = ap.parse_args()

    repo = Path(args.surfaceome_repo).expanduser()
    cache = Path(args.cache_dir).expanduser()
    out_dir = Path(args.out_dir)
    seqdir = repo / "data/external/deeptmhmm_surfaceome_inputs"

    library = sorted(set(json.loads(Path(args.library).read_text())))
    print(f"library members: {len(library)}")

    cohort = parse_uniprot_fasta(seqdir / "human_canonical_non_hla.fasta")
    human_iso = parse_uniprot_fasta(seqdir / "human_isoforms_from_afdb_non_hla.fasta")
    mouse = parse_uniprot_fasta(seqdir / "mouse_ortholog_one2one_highconf_non_hla.fasta")
    cyno = parse_uniprot_fasta(seqdir / "cyno_ortholog_one2one_highconf_non_hla.fasta")
    uniprot_human = load_uniprot_proteoforms(cache / UNIPROT_HUMAN)
    uniprot_mouse = load_uniprot_proteoforms(cache / UNIPROT_MOUSE)
    print(
        f"stored: cohort {len(cohort)}, isoforms {len(human_iso)}, mouse {len(mouse)}, "
        f"cyno {len(cyno)}"
    )
    print(
        f"UniProt {UNIPROT_RELEASE} proteoforms: human "
        f"{sum(len(v) for v in uniprot_human.values())} over {len(uniprot_human)} genes, "
        f"mouse {sum(len(v) for v in uniprot_mouse.values())} over {len(uniprot_mouse)} genes"
    )

    cohort_by_symbol = {r["symbol"]: r for r in cohort.values() if r["symbol"]}
    isoforms_by_symbol: dict[str, list[dict]] = defaultdict(list)
    for iso in human_iso.values():
        if iso["symbol"]:
            isoforms_by_symbol[iso["symbol"]].append(iso)

    # Human accession -> ENSG. Joining on the accession keeps the symbol out of the key
    # path; it is only ever the library's own identifier.
    ensg_human: dict[str, set[str]] = defaultdict(set)
    with (repo / "data/external/uniprot_ensembl_xrefs/ensg_to_uniprot.tsv").open() as fh:
        next(fh)
        for line in fh:
            ensg, acc = line.split("\t")[:2]
            ensg_human[acc].add(ensg)

    # The ortholog cohorts were built by Compara one-to-one high-confidence orthology and
    # their metadata carries both ends of each pair, so the ortholog ENSG and the
    # human/mouse pairing both come from there rather than from a symbol match.
    mouse_meta = load_ortholog_metadata(
        seqdir / "mouse_ortholog_one2one_highconf_non_hla_metadata.csv"
    )
    cyno_meta = load_ortholog_metadata(
        seqdir / "cyno_ortholog_one2one_highconf_non_hla_metadata.csv"
    )
    mouse_ensg = {
        r["selected_uniprot_accession"]: {r["target_ensembl_gene_id"]}
        for r in mouse_meta
        if r["target_ensembl_gene_id"]
    }
    human_ensg_for_mouse_acc = {
        r["selected_uniprot_accession"]: r["query_ensembl_gene_id"] for r in mouse_meta
    }

    appris_h = parse_appris(cache / f"appris_homo_sapiens_{APPRIS_RELEASE}.txt")
    appris_m = parse_appris(cache / f"appris_mus_musculus_{APPRIS_RELEASE}.txt")
    pep_h = parse_ensembl_pep(cache / "Homo_sapiens.GRCh38.pep.all.fa.gz")
    pep_m = parse_ensembl_pep(cache / "Mus_musculus.GRCm39.pep.all.fa.gz")
    print(f"APPRIS genes: human {len(appris_h)}, mouse {len(appris_m)}")
    print(f"Ensembl {ENSEMBL_RELEASE} peptides: human {len(pep_h)}, mouse {len(pep_m)}")

    appris_h_by_symbol: dict[str, list[dict]] = defaultdict(list)
    for rows in appris_h.values():
        for r in rows:
            appris_h_by_symbol[r["symbol"]].append(r)
    appris_m_by_symbol: dict[str, list[dict]] = defaultdict(list)
    for rows in appris_m.values():
        for r in rows:
            appris_m_by_symbol[r["symbol"]].append(r)

    # ---- human, over the whole library -------------------------------------------
    res_h: list[dict] = []
    for symbol in library:
        stored = cohort_by_symbol.get(symbol)
        source = "deeptmhmm_cohort"
        if stored is None:
            pf = uniprot_human.get(symbol)
            stored = pf[0] if pf and not pf[0]["is_isoform"] else None
            source = f"uniprot_reviewed_{UNIPROT_RELEASE}"
        if stored is None:
            res_h.append(
                {"symbol": symbol, "status": "no_uniprot_canonical", "sequence_source": ""}
            )
            continue
        rec = audit_gene(
            stored,
            isoforms_by_symbol.get(symbol, []),
            uniprot_human.get(symbol, []),
            sorted(ensg_human.get(stored["accession"], ())),
            appris_h,
            pep_h,
            appris_h_by_symbol,
        )
        rec["sequence_source"] = source
        res_h.append(rec)

    # ---- mouse, over the ortholog cohort -----------------------------------------
    res_m = []
    for acc, stored in sorted(mouse.items(), key=lambda kv: kv[1]["symbol"] or ""):
        rec = audit_gene(
            stored,
            [],
            uniprot_mouse.get(stored["symbol"] or "", []),
            sorted(mouse_ensg.get(acc, ())),
            appris_m,
            pep_m,
            appris_m_by_symbol,
        )
        rec["sequence_source"] = "deeptmhmm_cohort"
        res_m.append(rec)

    annotate_allelic_variants(res_h + res_m, cache)

    write_tsv(out_dir / "appris_principal_human.tsv", res_h, COLUMNS)
    write_tsv(out_dir / "appris_principal_mouse.tsv", res_m, COLUMNS)

    # The specific ask: principals that are not one of the sequences we store. Split by
    # whether UniProt carries the sequence at all, because the two need different fixes --
    # one is a cohort gap we can fill, the other is APPRIS and UniProt disagreeing.
    not_stored_status = {"uniprot_proteoform_not_stored", "no_uniprot_match"}
    unmatched = [
        dict(r, species=sp)
        for sp, res in (("human", res_h), ("mouse", res_m))
        for r in res
        if r.get("status") in not_stored_status
    ]
    unmatched.sort(
        key=lambda r: (r["species"], r["status"], -abs(r.get("len_delta") or 0))
    )
    write_tsv(out_dir / "appris_principal_not_stored.tsv", unmatched, ["species"] + COLUMNS)

    # ---- human/mouse principal pair consistency ----------------------------------
    human_by_ensg = {
        g: r for r in res_h for g in str(r.get("ensembl_gene", "")).split(";") if g
    }
    pairs = []
    for m in res_m:
        hum_ensg = human_ensg_for_mouse_acc.get(m["accession"])
        h = human_by_ensg.get(hum_ensg) if hum_ensg else None
        if not h or not h.get("len_principal") or not m.get("len_principal"):
            continue
        lh, lm = h["len_principal"], m["len_principal"]
        rel = abs(lh - lm) / max(lh, lm)
        pairs.append(
            {
                "human_symbol": h["symbol"],
                "mouse_symbol": m["symbol"],
                "human_accession": h["accession"],
                "mouse_accession": m["accession"],
                "human_class": h.get("appris_class", ""),
                "mouse_class": m.get("appris_class", ""),
                "human_enst": h.get("appris_enst", ""),
                "mouse_enst": m.get("appris_enst", ""),
                "len_human": lh,
                "len_mouse": lm,
                "rel_length_diff": round(rel, 4),
                "consistent": rel <= PAIR_LENGTH_TOLERANCE,
            }
        )
    pairs.sort(key=lambda r: -r["rel_length_diff"])
    write_tsv(out_dir / "appris_pair_consistency.tsv", pairs, list(pairs[0]))

    # ---- cyno: no APPRIS coverage exists for M. fascicularis ---------------------
    human_symbol_for_ensg = {
        g: r["symbol"] for r in res_h for g in str(r.get("ensembl_gene", "")).split(";") if g
    }
    cyno_rows = []
    for r in cyno_meta:
        acc = r["selected_uniprot_accession"]
        if acc not in cyno:
            continue
        cyno_rows.append(
            {
                "symbol": cyno[acc]["symbol"],
                "accession": acc,
                "len_stored": len(cyno[acc]["seq"]),
                "cyno_ensembl_gene": r["target_ensembl_gene_id"],
                "human_ensembl_gene": r["query_ensembl_gene_id"],
                "human_symbol": human_symbol_for_ensg.get(r["query_ensembl_gene_id"], ""),
                "in_library": r["query_ensembl_gene_id"] in human_symbol_for_ensg,
            }
        )
    cyno_rows.sort(key=lambda r: r["symbol"] or "")
    write_tsv(
        out_dir / "cyno_no_appris.tsv",
        cyno_rows,
        [
            "symbol",
            "accession",
            "len_stored",
            "cyno_ensembl_gene",
            "human_ensembl_gene",
            "human_symbol",
            "in_library",
        ],
    )

    def variant_verdicts(res):
        sub = [r for r in res if r.get("difference_kind") == "substitution_only"]
        return dict(
            sorted(
                (
                    (v, sum(1 for r in sub if r.get("substitutions_are_annotated_variants") == v))
                    for v in {
                        r.get("substitutions_are_annotated_variants", "") for r in sub
                    }
                ),
                key=lambda kv: -kv[1],
            )
        )

    def difference_kinds(res):
        unmatched = [r for r in res if r.get("status") == "no_uniprot_match"]
        return dict(
            sorted(
                (
                    (k, sum(1 for r in unmatched if r.get("difference_kind") == k))
                    for k in {r.get("difference_kind", "") for r in unmatched}
                ),
                key=lambda kv: -kv[1],
            )
        )

    def tally(res):
        t = defaultdict(int)
        for r in res:
            t[r.get("status", "?")] += 1
        return dict(sorted(t.items(), key=lambda kv: -kv[1]))

    rels = [p["rel_length_diff"] for p in pairs]
    summary = {
        "appris_release": APPRIS_RELEASE,
        "ensembl_release": ENSEMBL_RELEASE,
        "uniprot_release": UNIPROT_RELEASE,
        "library_members": len(library),
        "human": {
            "status": tally(res_h),
            "sequence_source": dict(
                sorted(
                    ((s, sum(1 for r in res_h if r.get("sequence_source") == s))
                     for s in {r.get("sequence_source", "") for r in res_h}),
                    key=lambda kv: -kv[1],
                )
            ),
            "class": dict(
                sorted(
                    ((lab, sum(1 for r in res_h if r.get("appris_class") == lab))
                     for lab in PRINCIPAL_ORDER),
                    key=lambda kv: -kv[1],
                )
            ),
            "ambiguous_principal": sum(1 for r in res_h if r.get("ambiguous_principal")),
            "principal_is_mane_select": sum(1 for r in res_h if r.get("mane_select")),
            "ensg_from_symbol_fallback": sum(
                1 for r in res_h if r.get("ensg_source") == "appris_symbol"
            ),
            "no_uniprot_match_by_difference": difference_kinds(res_h),
            "substitution_only_variant_verdict": variant_verdicts(res_h),
        },
        "mouse": {
            "cohort": len(mouse),
            "status": tally(res_m),
            "no_uniprot_match_by_difference": difference_kinds(res_m),
            "substitution_only_variant_verdict": variant_verdicts(res_m),
        },
        "pairs": {
            "n": len(pairs),
            "consistent": sum(1 for p in pairs if p["consistent"]),
            "inconsistent": sum(1 for p in pairs if not p["consistent"]),
            "tolerance": PAIR_LENGTH_TOLERANCE,
            "rel_length_diff_median": round(statistics.median(rels), 4) if rels else None,
            "rel_length_diff_p90": (
                round(sorted(rels)[int(0.9 * len(rels))], 4) if rels else None
            ),
        },
        "cyno": {
            "cohort": len(cyno),
            "paired_to_a_library_gene": sum(1 for r in cyno_rows if r["in_library"]),
            "appris_coverage": (
                "none - APPRIS e115v50 covers Macaca mulatta, not M. fascicularis (taxon 9541)"
            ),
        },
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
