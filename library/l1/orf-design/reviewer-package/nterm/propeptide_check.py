"""Does any N-after-SP cassette sit on sequence the mature protein loses with a pro-peptide?

A cassette placed after the signal-peptide cleavage site still lies upstream of the mature
N-terminus when the protein carries an N-terminal pro-peptide: the convertase that removes the
pro-peptide (furin or a related PC, at R-X-[KR]-R) removes the cassette with it. The site rule
(``admissible_interval.py``) knows the signal-peptide call and disorder only, so this checks
every decided arm against two independent sources:

- **UniProt PROPEP features** on the human reviewed entry for the gene, mapped onto the
  designed protein by sequence: the feature's last 8 residues plus the next 8 (the mature
  N-terminus) must occur exactly once in the designed protein, otherwise the row says
  ``unmapped``;
- **the furin motif** R-X-[KR]-R in the designed protein between the cleavage site and the
  first 60 mature residues, flagged when the cassette lies upstream of the motif's cleavage.

``snapshot`` freezes the UniProt rows (with the release and the SHA-256 of the file) to S3;
``check`` reads that snapshot only, so the check is reproducible without the network.

**Attribution** (convention of 2026-11-10): written by the Claude Code session
``nterm-reviewer-package`` (PR #254), 2026-11-12, for the overnight design review the
coordinator relayed from Becca ("review every step of the in-silico design").

    uv run --locked python propeptide_check.py snapshot
    uv run --locked python propeptide_check.py check
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

PREFIX = "sprint4/s4e9-orf-design-fg-panel/reviewer-package/nterm/uniprot-propeptides/"
FIELDS = "accession,gene_primary,sequence,ft_signal,ft_propep,ft_chain"
API = "https://rest.uniprot.org/uniprotkb/stream"
OUT = HERE / "propeptide_check.tsv"
SUMMARY = HERE / "propeptide_check_summary.json"
FURIN = re.compile(r"(?=(R.[KR]R))")
FURIN_WINDOW = 60  # mature residues past the cleavage site scanned for a furin motif
FLANK = 8
#: An N-terminal pro-peptide starts at the residue after the signal peptide; UniProt and the
#: designed cleavage call can differ by a few residues.
START_SLACK = 5


def _spans(text: str, kind: str) -> list[tuple[int, int]]:
    return [(int(a), int(b)) for a, b in re.findall(kind + r" (\d+)\.\.(\d+)", text or "")]


def snapshot() -> int:
    """Query UniProt (human, reviewed) for every gene on the rung and publish the rows."""
    from deliverome_analysis.s3 import processed

    with (HERE / "nterm_decisions.tsv").open() as fh:
        genes = sorted({r["gene_symbol"] for r in csv.DictReader(fh, delimiter="\t")})
    rows, release = [], None
    for i in range(0, len(genes), 80):
        q = " OR ".join(f"gene_exact:{g}" for g in genes[i : i + 80])
        url = (
            API
            + "?"
            + urllib.parse.urlencode(
                {
                    "query": f"({q}) AND organism_id:9606 AND reviewed:true",
                    "fields": FIELDS,
                    "format": "tsv",
                }
            )
        )
        for attempt in range(5):
            try:
                with urllib.request.urlopen(url, timeout=120) as r:
                    release = r.headers.get("X-UniProt-Release") or release
                    body = r.read().decode()
                break
            except Exception:  # retry transient HTTP errors
                time.sleep(2**attempt)
        else:
            raise SystemExit(f"UniProt query failed for genes {i}..{i + 80}")
        part = list(csv.DictReader(io.StringIO(body), delimiter="\t"))
        rows += part
        print(f"genes {i}-{i + 80}: {len(part)} entries", flush=True)
    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(rows)
    data = buf.getvalue().encode()
    sha = hashlib.sha256(data).hexdigest()
    key = PREFIX + f"{release}/uniprot_features.tsv"
    local = HERE / "uniprot_features.tsv"
    local.write_bytes(data)
    processed.upload(str(local), key)
    meta = {
        "release": release,
        "key": key,
        "sha256": sha,
        "entries": len(rows),
        "genes_queried": len(genes),
        "query": "gene_exact AND organism_id:9606 AND reviewed:true",
        "fields": FIELDS,
        "fetched_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }
    (HERE / "uniprot_features_snapshot.json").write_text(json.dumps(meta, indent=1) + "\n")
    print(json.dumps(meta, indent=1))
    return 0


def _download(key: str, dest: str) -> str:
    """``external://…`` keys read the append-only external bucket; bare keys read processed."""
    from deliverome_analysis.s3 import external, processed

    if key.startswith("external://"):
        return external.download(key.removeprefix("external://"), dest)
    return processed.download(key, dest)


def _load_snapshot() -> list[dict]:

    meta = json.loads((HERE / "uniprot_features_snapshot.json").read_text())
    local = HERE / "uniprot_features.tsv"
    if not local.exists():
        _download(meta["key"], str(local))
    data = local.read_bytes()
    assert hashlib.sha256(data).hexdigest() == meta["sha256"], "UniProt snapshot hash mismatch"
    rows = list(csv.DictReader(io.StringIO(data.decode()), delimiter="\t"))
    # Genes that joined the rung after the snapshot, fetched with the same query (2026-11-21).
    for sup in meta.get("supplements", []):
        local = HERE / f"uniprot_features.{sup['sha256'][:12]}.tsv"
        if not local.exists():
            _download(sup["key"], str(local))
        extra = local.read_bytes()
        assert hashlib.sha256(extra).hexdigest() == sup["sha256"], "UniProt supplement mismatch"
        rows += list(csv.DictReader(io.StringIO(extra.decode()), delimiter="\t"))
    return rows


#: Pro-peptides removed only when an extracellular protease activates the protein, not during
#: secretion: the PAR activation peptides (the tethered ligand starts at the new N-terminus) and
#: trypsin-fold zymogen activation peptides (the new N-terminus inserts into the activation
#: pocket). A cassette there is present on the resting surface protein; moving it past the
#: pro-peptide would block activation, so these arms keep their site. Decided by
#: nterm-reviewer-package (delegated by Becca 2026-11-12), from the UniProt notes ("Removed for
#: receptor activation", "Activation peptide") and the protease family.
ACTIVATION_ONLY = {
    "F2R": "PAR1: thrombin removes it to expose the tethered ligand",
    "F2RL1": "PAR2: trypsin removes it to expose the tethered ligand",
    "F2RL2": "PAR3: thrombin removes it to expose the tethered ligand",
    "F2RL3": "PAR4: thrombin removes it to expose the tethered ligand",
    "PRSS8": "prostasin zymogen activation peptide; the mature N-terminus forms the active site",
    "PRSS21": "testisin zymogen activation peptide; the mature N-terminus forms the active site",
    "PRSS41": "trypsin-fold zymogen activation peptide; the mature N-terminus forms the active site",
    "MEP1B": "meprin beta pro-peptide, removed by trypsin outside the cell",
}

#: Pro-peptides cut but not released: the fragments stay disulfide-linked on the surface, so a
#: cassette in the N-terminal fragment is kept.
STAYS_ATTACHED = {
    "RGMA": "autocatalytic cut at 168/169; by similarity to RGMB/RGMC the N- and C-terminal "
    "fragments stay disulfide-linked",
}


def entries_by_gene() -> dict[str, list[dict]]:
    out: dict[str, list[dict]] = {}
    for e in _load_snapshot():
        for g in (e["Gene Names (primary)"] or "").split(";"):
            out.setdefault(g.strip(), []).append(e)
    return out


def _find(motif: str, seq: str, max_mismatch: int = 2) -> list[int]:
    hits = [m.start() for m in re.finditer(re.escape(motif), seq)]
    if hits:
        return hits
    k = len(motif)
    best = [
        i
        for i in range(len(seq) - k + 1)
        if sum(a != b for a, b in zip(motif, seq[i : i + k], strict=True)) <= max_mismatch
    ]
    return best


def nterm_propeptides(
    seq: str, entries: list[dict], native_cleavage: int
) -> tuple[str, list[tuple]]:
    """N-terminal UniProt pro-peptides mapped onto ``seq``: (start, end, how) in its numbering.

    N-terminal means the pro-peptide starts within START_SLACK of the end of UniProt's own
    signal peptide (or of ``native_cleavage`` when UniProt has none). It is mapped by the 16
    residues around its last residue, exactly or with at most 2 mismatches (orthologs), and
    left ``unmapped`` when that does not give one place.
    """
    status, props = ("no_entry" if not entries else "no_propeptide"), []
    for e in entries:
        useq = e["Sequence"]
        sig = _spans(e["Signal peptide"], "SIGNAL")
        sp_end = sig[0][1] if sig else native_cleavage
        for a, b in _spans(e["Propeptide"], "PROPEP"):
            if (
                a > sp_end + 1 + START_SLACK
            ):  # not N-terminal: C-terminal GPI-attachment peptides etc.
                continue
            if useq == seq:
                props.append((a, b, "identical"))
                continue
            lo = max(0, b - FLANK)
            hits = _find(useq[lo : b + FLANK], seq)
            if len(hits) == 1:
                end = hits[0] + (b - lo)
                props.append((max(1, end - (b - a)), end, "mapped"))
            else:
                props.append((a, b, "unmapped"))
    return ("propeptide" if props else status), props


def native_cleavage(r: dict) -> int:
    if r["arm"] == "igk_swap":
        return int(r["igk_replaces_residues"].split("-")[1])
    return int(r["cleavage_site"] or 0)


def check() -> int:
    from deliverome_analysis.s3 import processed

    members, _ = npl.load_members(processed)
    entries = entries_by_gene()
    with (HERE / "nterm_decisions.tsv").open() as fh:
        dec = [r for r in csv.DictReader(fh, delimiter="\t")]
    out = []
    for r in dec:
        if r["candidate_id"] not in members:  # routed in from another rung (ATRN, NRG2); no arm
            assert r["arm"] in ("no_nterm_arm", "routed_internal"), r["candidate_id"]
            continue
        seq = members[r["candidate_id"]]["sequence"]
        built = r["arm"] in ("native_site", "igk_swap")
        c = native_cleavage(r) if built else 0
        n = int(r["native_after_residue"]) if built else None
        status, props = nterm_propeptides(seq, entries.get(r["gene_symbol"], []), c)
        furin = [
            m.start() + 4 for m in FURIN.finditer(seq) if c < m.start() + 4 <= c + FURIN_WINDOW
        ]
        lost_prop = [
            p for p in props if n is not None and p[2] != "unmapped" and p[0] - 1 <= n < p[1]
        ]
        lost_furin = [p for p in furin if n is not None and n < p]
        out.append(
            {
                "candidate_id": r["candidate_id"],
                "design_id": r["design_id"],
                "gene_symbol": r["gene_symbol"],
                "arm": r["arm"],
                "native_cleavage": c if built else "",
                "cassette_after_native_residue": "" if n is None else n,
                "uniprot": status,
                "n_terminal_propeptides": ";".join(f"{a}-{b}:{how}" for a, b, how in props),
                "propeptide_removal": (
                    ""
                    if not props
                    else "activation_only"
                    if r["gene_symbol"] in ACTIVATION_ONLY
                    else "stays_attached"
                    if r["gene_symbol"] in STAYS_ATTACHED
                    else "constitutive"
                ),
                "furin_motifs_cleave_after": ";".join(map(str, furin)),
                "cassette_inside_propeptide": int(bool(lost_prop)),
                "cassette_upstream_of_furin_motif": int(bool(lost_furin)),
            }
        )
    with OUT.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out)
    built = [o for o in out if o["arm"] in ("native_site", "igk_swap")]
    summary = {
        "uniprot_snapshot": json.loads((HERE / "uniprot_features_snapshot.json").read_text()),
        "members": len(out),
        "uniprot_status": {
            k: sum(o["uniprot"] == k for o in out)
            for k in ("propeptide", "no_propeptide", "no_entry")
        },
        "built_arms": len(built),
        "cassette_inside_constitutive_propeptide": [
            o["gene_symbol"]
            for o in built
            if o["cassette_inside_propeptide"] and o["propeptide_removal"] == "constitutive"
        ],
        "cassette_inside_activation_only_propeptide": [
            o["gene_symbol"]
            for o in built
            if o["cassette_inside_propeptide"] and o["propeptide_removal"] == "activation_only"
        ],
        "cassette_upstream_of_furin_motif_without_annotation": [
            o["gene_symbol"]
            for o in built
            if o["cassette_upstream_of_furin_motif"] and not o["cassette_inside_propeptide"]
        ],
        "unmapped_propeptides": [
            o["gene_symbol"] for o in out if "unmapped" in o["n_terminal_propeptides"]
        ],
    }
    SUMMARY.write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps(summary, indent=1))
    return 0


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "check"
    raise SystemExit(snapshot() if cmd == "snapshot" else check())
