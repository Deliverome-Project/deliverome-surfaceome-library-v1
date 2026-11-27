"""Second-source review of members whose signal peptide only DeepTMHMM v1 calls.

Class ``refused_sp_disagree_deeptmhmm_v1_only`` (58 members): DeepTMHMM v1 calls a signal
peptide 1..c and SignalP 6.0 does not. The Igκ swap deletes residues 1..c, so it is only right
if 1..c is a cleaved signal peptide. If the segment is a signal anchor or the first
transmembrane helix, deleting it changes the protein's topology.

Rule (decided by the coordinator under Becca's 2026-11-12 delegation, 2026-11-14; implemented
by the Claude Code session nterm-reviewer-package):

| Evidence on the designed protein | Verdict |
|---|---|
| UniProt SIGNAL ends within ``SITE_TOLERANCE`` of c, or DeepTMHMM v2 calls a signal peptide ending there | ``cleaved_sp_confirmed``: keep the Igκ swap |
| UniProt TRANSMEM, or a DeepTMHMM v2 TM helix, starts within the first ``ANCHOR_WINDOW`` residues | ``signal_anchor_or_tm``: do not delete; route by the N-terminus side |
| neither | ``unconfirmed``: hold |

Sources are matched on sequence, never accession: a UniProt entry counts only when its first
``NTERM_MATCH`` residues equal the designed protein's, and DeepTMHMM v2 is read for exactly the
designed sequence (``topology.for_sequence(..., tool=TOOL_V2)``). A routed member's rung is
read from DeepTMHMM v2 (UniProt topology when v2 has no row): N-terminus outside with one TM
helix → type III residue 1; outside with several → residue-1 multipass; N-terminus inside →
C-terminal when the C-terminus is outside, else internal.

    uv run --locked python sp_disagree_review.py snapshot   # UniProt rows, frozen to S3
    uv run --locked python sp_disagree_review.py review
"""

from __future__ import annotations

import csv
import hashlib
import io
import json
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import nterm_placement as npl  # noqa: E402

CLASS = "refused_sp_disagree_deeptmhmm_v1_only"
SITE_TOLERANCE = 3
ANCHOR_WINDOW = 60
NTERM_MATCH = 70
PREFIX = "sprint4/s4e9-orf-design-fg-panel/reviewer-package/nterm/uniprot-sp-disagree/"
FIELDS = "accession,gene_primary,sequence,ft_signal,ft_transmem,ft_topo_dom"
API = "https://rest.uniprot.org/uniprotkb/stream"
SNAP_META = HERE / "uniprot_sp_disagree_snapshot.json"
SNAP_LOCAL = HERE / "uniprot_sp_disagree.tsv"
OUT = HERE / "sp_disagree_review.tsv"
DECIDED_BY = "coordinator under Becca's 2026-11-12 delegation, 2026-11-14"
IMPLEMENTED_BY = "Claude Code session nterm-reviewer-package"


def _tsv(path: Path) -> list[dict]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


#: Members with no census SignalP row whose own native SignalP 6.0 call (pinned runs, including
#: ``nterm_sp6_2026_11_14h`` on the six orthologs round 2 never scored natively) is not a
#: signal peptide. Their Igκ swaps rest on DeepTMHMM v1 alone too (coordinator, 2026-11-14).
EXTRA_CLASS = "refused_no_signalp_row"


def _native_not_sp() -> set[str]:
    import hashlib

    import nterm_overnight as no

    from deliverome_analysis.core.s3 import processed

    calls = no.prior_calls(processed)
    seqs, _ = npl.load_members(processed)
    out = set()
    for r in _tsv(HERE / "nterm_members.tsv"):
        if r["class"] != EXTRA_CLASS or r["candidate_id"] not in seqs:
            continue
        h = hashlib.sha256(seqs[r["candidate_id"]]["sequence"].encode()).hexdigest()
        call = calls.get(h)
        if call and call["prediction"] != "SP":
            out.add(r["candidate_id"])
    return out


def members() -> list[dict]:
    extra = _native_not_sp()
    return [
        r
        for r in _tsv(HERE / "nterm_members.tsv")
        if r["class"] == CLASS or r["candidate_id"] in extra
    ]


def snapshot() -> int:
    from deliverome_analysis.core.s3 import processed

    genes = sorted({r["gene_symbol"] for r in members()})
    q = " OR ".join(f"gene_exact:{g}" for g in genes)
    url = (
        API
        + "?"
        + urllib.parse.urlencode(
            {
                "query": f"({q}) AND organism_id:9606 AND reviewed:true",
                "fields": FIELDS,
                "format": "tsv",
                "includeIsoform": "true",
            }
        )
    )
    with urllib.request.urlopen(url, timeout=180) as r:
        release = r.headers.get("X-UniProt-Release")
        data = r.read()
    SNAP_LOCAL.write_bytes(data)
    key = PREFIX + f"{release}/uniprot_features.tsv"
    processed.upload(str(SNAP_LOCAL), key)
    meta = {
        "release": release,
        "key": key,
        "sha256": hashlib.sha256(data).hexdigest(),
        "genes": len(genes),
        "entries": data.decode().count("\n") - 1,
        "fields": FIELDS,
        "query": "gene_exact AND organism_id:9606 AND reviewed:true, isoforms included",
        "fetched_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }
    SNAP_META.write_text(json.dumps(meta, indent=1) + "\n")
    print(json.dumps(meta, indent=1))
    return 0


#: UniProt evidence codes seen on SIGNAL features (Evidence & Conclusion Ontology).
ECO = {
    "ECO:0000269": "experimental evidence",
    "ECO:0000250": "by similarity",
    "ECO:0000255": "sequence-analysis prediction",
    "ECO:0000305": "curator inference",
}
#: A primary paper reporting cleavage at a site only DeepTMHMM v2 calls. Verified against the
#: abstract on 2026-11-14 (nterm-reviewer-package): signal peptidase "cleaves off a large
#: N-terminal signal peptide of about 64 amino acids" from CNNM2.
PAPERS = {
    "CNNM2": "de Baaij et al. 2012, J Biol Chem 287:13644, PMID 22399287 "
    "(doi:10.1074/jbc.M112.342204): signal peptidase cleaves an N-terminal signal peptide "
    "of about 64 residues",
}


def _eco(field: str) -> str:
    m = re.search(r"SIGNAL [^;]*;[^S]*?/evidence=\"+(ECO:\d+)", field or "")
    return m.group(1) if m else ""


def _spans(text: str, kind: str) -> list[tuple[int, int, str]]:
    out = []
    for m in re.finditer(kind + r" (\d+)\.\.(\d+)([^;]*;\s*/note=\"([^\"]*)\")?", text or ""):
        out.append((int(m.group(1)), int(m.group(2)), m.group(4) or ""))
    return out


def _load_uniprot() -> dict[str, list[dict]]:
    from deliverome_analysis.core.s3 import processed

    meta = json.loads(SNAP_META.read_text())
    if not SNAP_LOCAL.exists():
        processed.download(meta["key"], str(SNAP_LOCAL))
    data = SNAP_LOCAL.read_bytes()
    assert hashlib.sha256(data).hexdigest() == meta["sha256"], "UniProt snapshot hash mismatch"
    out: dict[str, list[dict]] = {}
    for e in csv.DictReader(io.StringIO(data.decode()), delimiter="\t"):
        for g in (e["Gene Names (primary)"] or "").split(";"):
            out.setdefault(g.strip(), []).append(e)
    return out


def dtm2(seq: str) -> dict | None:
    """DeepTMHMM v2 on exactly this sequence: SP end, TM helices, N- and C-terminus side."""
    from deliverome_analysis.surfaceome import topology as topo

    row = topo.for_sequence(seq, tool=topo.TOOL_V2)
    if not row:
        return None
    raw = row.get("dtm2_topology_string")
    if not isinstance(raw, str) or not raw:
        return None
    sp_end = len(raw) - len(raw.lstrip("S")) if raw.startswith("S") else 0
    tms = [(m.start() + 1, m.end()) for m in re.finditer(r"M+", raw)]
    # Membrane sides need the v1-alphabet projection (I/O); the raw v2 alphabet encodes
    # compartments with other letters. Without the projection the side is left blank.
    lossy = row.get("dtm2_v1_alphabet_lossy")
    lossy = lossy if isinstance(lossy, str) else ""
    nside = next((ch for ch in lossy[sp_end:] if ch in "IO"), "")
    cside = next((ch for ch in reversed(lossy) if ch in "IO"), "")
    return {
        "sp_end": sp_end,
        "tms": tms,
        "n_side": nside,
        "c_side": cside,
        "version": row.get("topology_version", ""),
    }


def _side(note: str) -> str:
    if "Extracellular" in note or "Lumenal" in note:
        return "O"
    return "I" if "Cytoplasmic" in note else ""


ORTHO_SITE_IDENTITY = 0.80
ORTHO_WINDOW = 90


def human_ortholog_signal(seq: str, entries: list[dict]) -> dict | None:
    """The human canonical UniProt SIGNAL aligned onto an ortholog's N-terminus."""
    from Bio import Align
    from Bio.Align import substitution_matrices

    canon = next((e for e in entries if "-" not in e["Entry"]), None)
    if not canon:
        return None
    sp = _spans(canon["Signal peptide"], "SIGNAL")
    if not sp:
        return None
    h = sp[0][1]
    eco = _eco(canon["Signal peptide"])
    eco_txt = f"{eco} ({ECO.get(eco, 'see UniProt')})" if eco else "no evidence code"
    al = Align.PairwiseAligner(mode="global", open_gap_score=-10, extend_gap_score=-0.5)
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    a = al.align(canon["Sequence"][:ORTHO_WINDOW], seq[:ORTHO_WINDOW])[0]
    # human residue i (1-based) -> ortholog residue j (1-based) where both are aligned
    pos = {}
    for (hs, he), (os_, _oe) in zip(*a.aligned, strict=True):
        for k in range(he - hs):
            pos[hs + k + 1] = os_ + k + 1
    hseq, oseq = canon["Sequence"], seq

    def ident(lo: int, hi: int) -> float:
        cols = range(max(1, lo), hi + 1)
        same = sum(1 for i in cols if i in pos and hseq[i - 1] == oseq[pos[i] - 1])
        return same / len(cols)

    sp_id, site_id = ident(1, h), ident(h - 4, h + 5)
    gap = not (h in pos and h + 1 in pos and pos[h + 1] == pos[h] + 1)
    site = pos.get(h, 0)
    passes = site_id >= ORTHO_SITE_IDENTITY and not gap and site > 0
    text = (
        f"human ortholog UniProt {canon['Entry']} SIGNAL 1-{h} [{eco_txt}] aligns to ortholog 1-{site}: "
        f"{sp_id:.0%} identity over the signal peptide, {site_id:.0%} over +/-5 of the site, "
        f"{'a gap' if gap else 'no gap'} at the site"
    )
    return {
        "text": text,
        "passes": passes,
        "site": site,
        "eco": eco_txt,
        "source": f"human ortholog UniProt {canon['Entry']} SIGNAL 1-{h}, aligned ({site_id:.0%} over +/-5)",
    }


def route(n_side: str, c_side: str, n_tm: int) -> str:
    if n_side == "O":
        return "type III residue 1" if n_tm == 1 else "residue-1 multipass"
    if n_side == "I":
        return "C-terminal" if c_side == "O" else "internal"
    return ""


def review() -> int:
    from deliverome_analysis.core.s3 import processed

    seqs, _ = npl.load_members(processed)
    uni = _load_uniprot()
    rows = []
    for m in members():
        cid = m["candidate_id"]
        seq = seqs[cid]["sequence"]
        c = int(m["deeptmhmm_v1_cleavage"])
        checked, sig, anchor = [], [], []
        topo_u = None
        hits = [
            e
            for e in uni.get(m["gene_symbol"], [])
            if e["Sequence"][:NTERM_MATCH] == seq[:NTERM_MATCH]
        ]
        if hits:
            e = hits[0]
            checked.append(f"UniProt {e['Entry']} (first {NTERM_MATCH} aa identical)")
            for _a, b, _ in _spans(e["Signal peptide"], "SIGNAL"):
                if abs(b - c) <= SITE_TOLERANCE:
                    sig.append(f"UniProt SIGNAL 1-{b}")
                else:
                    checked.append(f"UniProt SIGNAL 1-{b} (not within {SITE_TOLERANCE} of {c})")
            tm = [(a, b) for a, b, _ in _spans(e["Transmembrane"], "TRANSMEM")]
            if tm and tm[0][0] <= ANCHOR_WINDOW:
                note = _spans(e["Transmembrane"], "TRANSMEM")[0][2]
                anchor.append(
                    f"UniProt TRANSMEM {tm[0][0]}-{tm[0][1]}" + (f" ({note})" if note else "")
                )
            doms = _spans(e["Topological domain"], "TOPO_DOM")
            if tm:
                first = next((n for a, b, n in doms if b < tm[0][0]), "")
                last = next((n for a, b, n in reversed(doms) if a > tm[-1][1]), "")
                topo_u = (_side(first), _side(last), len(tm))
        else:
            checked.append(
                f"UniProt: no reviewed entry or isoform with the designed first {NTERM_MATCH} aa"
            )
        v2 = dtm2(seq)
        if v2:
            checked.append(f"DeepTMHMM v2 ({v2['version']})")
            if v2["sp_end"] and abs(v2["sp_end"] - c) <= SITE_TOLERANCE:
                sig.append(f"DeepTMHMM v2 SP 1-{v2['sp_end']}")
            elif v2["sp_end"]:
                checked.append(
                    f"DeepTMHMM v2 SP 1-{v2['sp_end']} (not within {SITE_TOLERANCE} of {c})"
                )
            if v2["tms"] and v2["tms"][0][0] <= ANCHOR_WINDOW and not v2["sp_end"]:
                anchor.append(f"DeepTMHMM v2 TM {v2['tms'][0][0]}-{v2['tms'][0][1]}")
        else:
            checked.append("DeepTMHMM v2: no row for the designed sequence")
        # A cleaved signal peptide at another site: UniProt's curated end first, else v2's.
        uni_end = next(
            (b for _a, b, _ in (_spans(hits[0]["Signal peptide"], "SIGNAL") if hits else [])), 0
        )
        reanchor, reanchor_src, reanchor_ev, hold = 0, "", "", ""
        if uni_end:
            eco = _eco(hits[0]["Signal peptide"])
            reanchor, reanchor_src = uni_end, f"UniProt {hits[0]['Entry']} SIGNAL 1-{uni_end}"
            reanchor_ev = (
                f"{eco} ({ECO.get(eco, 'see UniProt')})"
                if eco
                else "no evidence code on the UniProt feature"
            )
        elif v2 and v2["sp_end"] and m["gene_symbol"] in PAPERS:
            reanchor, reanchor_src = v2["sp_end"], f"DeepTMHMM v2 SP 1-{v2['sp_end']}"
            reanchor_ev = PAPERS[m["gene_symbol"]]
        elif v2 and v2["sp_end"]:
            hold = (
                f"only DeepTMHMM v2 calls a cleaved signal peptide (1-{v2['sp_end']}); no curated "
                "UniProt SIGNAL on this sequence, the residues around the v2 site do not align "
                "with a canonical UniProt SIGNAL, and no verified paper reports cleavage there"
            )
        else:
            hold = "no independent source calls a cleaved signal peptide"
        # Ortholog with no same-sequence UniProt entry: the human canonical entry's SIGNAL,
        # carried over by alignment (coordinator, 2026-11-14). Supporting evidence for a kept
        # swap; for a held member it re-anchors at the aligned site when the site region is
        # conserved (>= ORTHO_SITE_IDENTITY over +/-5, no gap at the site).
        ortho = ""
        if not hits and m["species"] != "human":
            al = human_ortholog_signal(seq, uni.get(m["gene_symbol"], []))
            if al:
                ortho = al["text"]
                if not sig and not anchor and not reanchor and al["passes"]:
                    reanchor, reanchor_src = al["site"], al["source"]
                    reanchor_ev, hold = al["eco"], ""
        if sig:
            verdict, to = "cleaved_sp_confirmed", ""
        elif anchor:
            verdict = "signal_anchor_or_tm"
            if v2 and v2["tms"] and v2["n_side"]:
                to = route(v2["n_side"], v2["c_side"], len(v2["tms"]))
            elif topo_u and topo_u[0]:
                to = route(*topo_u)
            else:
                to = ""
        else:
            verdict, to = "unconfirmed", ""
        rows.append(
            {
                "candidate_id": cid,
                "design_id": m["design_id"],
                "gene_symbol": m["gene_symbol"],
                "candidate_set": m["candidate_set"],
                "deeptmhmm_v1_cleavage": c,
                "sources_checked": "; ".join(checked),
                "cleaved_sp_evidence": "; ".join(sig),
                "anchor_or_tm_evidence": "; ".join(anchor),
                "verdict": verdict,
                "route_to": to,
                "human_ortholog_signal": ortho,
                "reanchor_cleavage": reanchor if verdict == "unconfirmed" and reanchor else "",
                "reanchor_source": reanchor_src if verdict == "unconfirmed" and reanchor else "",
                "reanchor_evidence": reanchor_ev if verdict == "unconfirmed" and reanchor else "",
                "hold_reason": hold if verdict == "unconfirmed" and not reanchor else "",
                "decided_by": DECIDED_BY,
                "implemented_by": IMPLEMENTED_BY,
            }
        )
    with OUT.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    from collections import Counter

    print(Counter((r["verdict"], r["route_to"]) for r in rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(snapshot() if sys.argv[1:2] == ["snapshot"] else review())
