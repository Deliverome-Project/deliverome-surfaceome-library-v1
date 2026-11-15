"""Topogenic features of the L1 members on the "N at residue 1, single-pass (type III)" rung.

Question: if the 51-aa tag cassette goes after Met1 of a member with no signal peptide
and one TM helix, does the TM still translocate the (now longer, more positive) N-terminus,
or should these members get an Igk leader ahead of the cassette?

For each member on the rung in ``l1_tag_paths.tsv`` this reads the member's OWN ordered
protein (designed CDS, matched by sequence hash) and its DeepTMHMM v1 string -- the
topology of record -- and reports:

* the native N-tail (residues before TM1, Met1 included) and its length once the cassette
  is inserted after Met1 (native + 51);
* the charge (K+R minus D+E) of the 15 residues on each side of TM1 and their difference,
  ``delta_charge = C-flank - N-flank`` (Hartmann, Rapoport & Lodish 1989: positive favours
  an extracellular N-terminus), for the native and the tagged protein;
* TM1 length and mean Kyte-Doolittle hydropathy;
* the v2 call on the same sequence (kept alongside, never substituted);
* whether the member's N-terminal region looks like an uncalled signal peptide (an
  N-terminal hydrophobic stretch), and whether the gene's UniProt record annotates one.

Only His is left uncharged; at cytosolic pH it is mostly neutral.

    uv run --locked python type3_nterm.py            # features -> type3_nterm_features.tsv
    uv run --locked python type3_nterm.py --refresh-uniprot
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import urllib.parse
import urllib.request
from pathlib import Path

from deliverome_analysis import library_design as ld
from deliverome_analysis import topology

HERE = Path(__file__).parent
TAG_PATHS = HERE / "l1_tag_paths.tsv"
MANIFEST = HERE / "candidate_manifest_tranches.tsv"
OUT = HERE / "type3_nterm_features.tsv"
UNIPROT_CACHE = HERE / "type3_uniprot_annotations.json"
ARMS_FASTA = HERE / "type3_arms.fasta"
RUNG = "N at residue 1, single-pass (type III)"

#: GGGGS | SpyTag003 | GGGGSGGGGS | P-ALFA v1.5-P | GGGGS, the fixed L1 cassette
#: (README, "Cassette encoding: only the linkers vary"). Asserted against the issued
#: encodings below so a cassette change cannot pass silently.
FULL_CASSETTE = "GGGGSRGVPHIVMVDAYKRYKGGGGSGGGGSPSRLEEELRRRLTEPGGGGS"
#: The full 51-aa cassette is kept at residue 1, leading GGGGS included (Becca, 2026-11-11,
#: "keep full"; implemented by session l1-nterm-type3). Dropping the GGGGS put SpyTag003's
#: VPHIVMV right after Met1, and DeepTMHMM v1 then called a cleavable signal peptide ending
#: inside SpyTag003 for 12 of the 62 tag-alone arms (0 of 62 with the full cassette).
CASSETTE = FULL_CASSETTE
#: Murine Igk V-J2-C leader + the +1 Asp every Igk vector keeps (build_gpcr_rescue_arms).
#: The Igk arm is leader + Asp + Met1 + cassette + native 2..end.
IGK = "METDTLLLWVLLLWVPGSTG"
IGK_PLUS1 = "D"
FLANK = 15

KD = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5, "G": -0.4,
    "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6, "S": -0.8,
    "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}  # fmt: skip


def arms(aa: str) -> dict[str, str]:
    return {
        "native": aa,
        "tag": aa[0] + CASSETTE + aa[1:],
        "igk_tag": IGK + IGK_PLUS1 + aa[0] + CASSETTE + aa[1:],
    }


def charge(seq: str) -> int:
    return sum(seq.count(a) for a in "KR") - sum(seq.count(a) for a in "DE")


def kd_mean(seq: str) -> float:
    return sum(KD.get(a, 0.0) for a in seq) / len(seq) if seq else 0.0


def segments(topo: str, state: str) -> list[tuple[int, int]]:
    """0-based half-open runs of ``state``."""
    return [(m.start(), m.end()) for m in re.finditer(f"{state}+", topo)]


def max_kd_window(seq: str, width: int = 8) -> float:
    if len(seq) < width:
        return kd_mean(seq)
    return max(kd_mean(seq[i : i + width]) for i in range(len(seq) - width + 1))


def features(aa: str, topo: str) -> dict:
    tms = segments(topo, "M")
    if not tms:
        return {"n_tm_v1": 0}
    s, e = tms[0]
    tm = aa[s:e]
    n_tail = aa[:s]
    n15, c15 = aa[max(0, s - FLANK) : s], aa[e : e + FLANK]
    tagged = aa[0] + CASSETTE + aa[1:]
    ts = s + len(CASSETTE)
    tn15 = tagged[max(0, ts - FLANK) : ts]
    return {
        "n_tm_v1": len(tms),
        "tm1_start": s + 1,
        "tm1_end": e,
        "tm1_len": e - s,
        "tm1_kd_mean": round(kd_mean(tm), 2),
        "n_tail_native_aa": s,
        "n_tail_tagged_aa": s + len(CASSETTE),
        "n_tail_charge_native": charge(n_tail),
        "n_tail_charge_tagged": charge(aa[0] + CASSETTE + n_tail[1:]),
        "n_flank15_charge_native": charge(n15),
        "n_flank15_charge_tagged": charge(tn15),
        "c_flank15_charge": charge(c15),
        "delta_charge_native": charge(c15) - charge(n15),
        "delta_charge_tagged": charge(c15) - charge(tn15),
        "c_tail_aa": len(aa) - e,
        # An uncalled signal peptide shows as a hydrophobic stretch in the first ~35
        # residues that is not TM1 itself.
        "nterm35_max_kd8": round(max_kd_window(aa[: min(35, s)]), 2) if s >= 8 else "",
        "n_cys_n_tail": n_tail.count("C"),
        "n_glyc_sequons_n_tail": len(re.findall(r"N[^P][ST]", n_tail)),
        "n_tail_head20": n_tail[:20],
    }


def canonical_offset(member: str, canonical: str) -> tuple[int | None, float]:
    """Canonical residue (1-based) that the member's residue 1 aligns to, and identity.

    Local alignment of the member's first 60 residues against the canonical's first 400,
    so an N-terminal truncation (alternative start, ortholog CDS) shows as an offset > 1
    and the signal peptide it would have skipped can be read off the UniProt record.
    """
    from Bio import Align
    from Bio.Align import substitution_matrices

    al = Align.PairwiseAligner(mode="local")
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    al.open_gap_score, al.extend_gap_score = -10, -0.5
    q, t = member[:60], canonical[:400]
    a = al.align(q, t)[0]
    (qs, qe), (ts, te) = a.aligned[0][0], a.aligned[1][0]
    ident = sum(x == y for x, y in zip(q[qs:qe], t[ts:te], strict=True)) / max(1, qe - qs)
    # project member residue 1 onto the canonical
    return ts - qs + 1, round(ident, 2)


def v2_summary(row: dict) -> dict:
    if not row.get("found"):
        return {"v2_found": 0}
    t = row.get("per_residue_topology") or ""
    return {
        "v2_found": 1,
        "v2_label": row.get("deeptmhmm_label", ""),
        "v2_sp_aa": int(row.get("signal_peptide_length") or 0),
        "v2_n_tm": int(row.get("tm_helix_count") or 0),
        "v2_n_term": row.get("n_terminal_orientation", ""),
        "v2_tm1_start": (t.find("M") + 1) if "M" in t else "",
    }


def uniprot_annotations(symbols: list[str], refresh: bool) -> dict:
    """Reviewed human UniProt record per HGNC symbol, cached in UNIPROT_CACHE.

    Looked up by exact gene name among reviewed human entries, and the returned record's
    primary gene name is checked against the query, so an alias hit cannot slip in.
    """
    cache = json.loads(UNIPROT_CACHE.read_text()) if UNIPROT_CACHE.exists() else {}
    for sym in symbols:
        if sym in cache and not refresh:
            continue
        q = urllib.parse.quote(f"gene_exact:{sym} AND organism_id:9606 AND reviewed:true")
        url = f"https://rest.uniprot.org/uniprotkb/search?query={q}&format=json&size=5"
        with urllib.request.urlopen(url, timeout=60) as fh:
            hits = json.load(fh)["results"]
        hits = [
            h
            for h in hits
            if any(g.get("geneName", {}).get("value") == sym for g in h.get("genes", []))
        ]
        if len(hits) != 1:
            cache[sym] = {"accession": None, "n_hits": len(hits)}
            continue
        d = hits[0]
        feats = d.get("features", [])

        def span(f):
            return [f["location"]["start"]["value"], f["location"]["end"]["value"]]

        topo_notes = []
        for c in d.get("comments", []):
            for loc in c.get("subcellularLocations", []):
                if "topology" in loc:
                    topo_notes.append(loc["topology"]["value"])
        cache[sym] = {
            "accession": d["primaryAccession"],
            "entry_version": d.get("entryAudit", {}).get("entryVersion"),
            "length": d.get("sequence", {}).get("length"),
            "sequence": d.get("sequence", {}).get("value"),
            "signal": [span(f) for f in feats if f["type"] == "Signal"],
            "transmembrane": [span(f) for f in feats if f["type"] == "Transmembrane"],
            "topo_dom": [
                [*span(f), f.get("description", "")]
                for f in feats
                if f["type"] == "Topological domain"
            ],
            "topology_notes": sorted(set(topo_notes)),
        }
    UNIPROT_CACHE.write_text(json.dumps(cache, indent=1, sort_keys=True) + "\n")
    return cache


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--refresh-uniprot", action="store_true")
    args = ap.parse_args(argv)

    with TAG_PATHS.open() as fh:
        rung = [r for r in csv.DictReader(fh, delimiter="\t") if r["placement"] == RUNG]
    print(f"{len(rung)} members on the {RUNG!r} rung")

    with (HERE / "l1_cassette_encodings.tsv").open() as fh:
        issued = {r["cassette_aa"] for r in csv.DictReader(fh, delimiter="\t")}
    assert issued == {FULL_CASSETTE}, f"issued cassettes differ from the fixed one: {issued}"

    members = ld.designed_proteins(MANIFEST, library=ld.LIBRARY_L1)
    # design_id is stable across member renumbering (candidate_id is plate order), so carry
    # it for rekeying.
    with MANIFEST.open() as fh:
        design_ids = {r["candidate_id"]: r["design_id"] for r in csv.DictReader(fh, delimiter="\t")}
    sub = {r["candidate_id"]: members[r["candidate_id"]] for r in rung}
    v1 = ld.topology_report(sub, tool=topology.TOOL_V1).set_index("candidate_id")
    v2 = ld.topology_report(sub, tool=topology.TOOL_V2).set_index("candidate_id")

    rows = []
    for r in rung:
        cid = r["candidate_id"]
        aa = sub[cid]["sequence"]
        t1 = v1.loc[cid]
        assert t1["found"] and topology.sequence_sha256(aa) == r["sequence_sha256"], cid
        assert str(t1["tool_version"]).startswith("deeptmhmm-"), (cid, t1["tool_version"])
        topo = t1["per_residue_topology"]
        assert topo[0] == "O" and "S" not in topo, cid
        rows.append(
            {
                "candidate_id": cid,
                "design_id": design_ids[cid],
                "candidate_set": r["candidate_set"],
                "species": r["species"],
                "gene_symbol": r["gene_symbol"],
                "length_aa": len(aa),
                "v1_topology_version": t1["topology_version"],
                "v1_tool_version": t1["tool_version"],
                "v1_label": t1.get("deeptmhmm_label", ""),
                "v1_uniprot_acc_full": t1.get("uniprot_acc_full", ""),
                **features(aa, topo),
                **v2_summary(v2.loc[cid].to_dict()),
                "shed_risk": r["shed_risk"],
            }
        )

    # Gene-level UniProt annotation for the human canonical, to separate members whose
    # gene is annotated type III from those whose gene carries a signal peptide the
    # ordered form lacks (or v1 missed).
    ann = uniprot_annotations(sorted({r["gene_symbol"] for r in rows}), args.refresh_uniprot)
    for row in rows:
        a = ann.get(row["gene_symbol"], {})
        row["uniprot_canonical"] = a.get("accession") or ""
        row["uniprot_signal"] = ";".join(f"{s}-{e}" for s, e in a.get("signal", []))
        row["uniprot_tm"] = ";".join(f"{s}-{e}" for s, e in a.get("transmembrane", []))
        row["uniprot_topology"] = ";".join(a.get("topology_notes", []))
        row["uniprot_length"] = a.get("length", "")
        if a.get("sequence"):
            off, ident = canonical_offset(sub[row["candidate_id"]]["sequence"], a["sequence"])
            row["member_res1_at_canonical"] = off
            row["nterm_identity"] = ident
            sp = a.get("signal") or []
            row["member_has_canonical_sp"] = "" if not sp else int(off == 1 and ident >= 0.5)

    # Three arms per member for re-prediction with DeepTMHMM v1: native (reproduces the
    # stored string), tag alone (cassette after Met1, the current rung design), and Igk
    # leader + Asp + the tag-alone protein (Met1, cassette, native 2..end) -- the
    # build_gpcr_rescue_arms layout, so the two arms differ by the leader alone.
    with ARMS_FASTA.open("w") as fh:
        for row in rows:
            aa = sub[row["candidate_id"]]["sequence"]
            for arm, seq in arms(aa).items():
                fh.write(f">sp|{row['candidate_id']}__{arm}|{row['gene_symbol']}\n{seq}\n")
    print(f"wrote {ARMS_FASTA.name}: {3 * len(rows)} sequences")

    cols: list[str] = []
    for row in rows:
        cols += [k for k in row if k not in cols]
    with OUT.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n", restval="")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {OUT.name}: {len(rows)} rows")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
