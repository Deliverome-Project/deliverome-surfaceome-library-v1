"""Plasmid QC over every insert of the L1 pool, plus the checks only a pool can fail.

``deliverome_analysis.plasmid_qc`` was written for single ordered plasmids (the P-series maps)
and had never run on the L1 library. This runs it on all of it (Becca, 2026-11-21: "we want the
plasmid qc to run on the pool and consider hardening it as well"; and "we still don't want
shared linkers in the middle").

**Per insert.** Each insert is placed into the p0039 backbone the way the pool is built
(``assemble_candidates.py``: the Esp3I-wrapped insert replaces the 195-nt MCS + stuffer slot
after EFS, keeping its two Esp3I sites so it can be shuttled later; p0039 has no Esp3I of its
own). Every check in ``plasmid_qc.CHECKS`` then runs on the finished plasmid, with a spec built
from the arm tables:

* ``designed_insert`` the insert itself, ``backbone`` the empty p0039, ``barcode`` the 20-N slot;
* ``cassette_nt`` the arm's issued encoding: ``l1_arm_cassette_encodings.tsv`` by ``arm_id`` when
  it exists, else ``l1_cassette_encodings.tsv`` for a C-terminal-table member; with neither, p0003's
  shared DNA, which ``diversified_linker`` and ``linkers_not_shared`` then fail;
* ``require_p003_tag_codons`` on (SpyTag003 and P-ALFA-P keep p003's codons, Becca 2026-10-29);
* for a plain insertion (not an Igκ leader arm, not a snorkel): the member's native protein as
  its proteoform, the cassette's residue, and the whole tagged protein;
* the p0039 features (EFS, IRES2, PuroR, WPRE, LTRs, origin...) that occur once in the empty
  backbone, which must still occur once.

On top of those, per insert: the encoded protein hashes to the arm table's ``protein_sha256``;
the layout matches the native (``l1_consistency.layout_check``); the insert opens with the
Esp3I site, Kozak and ATG and ends with the pool's single readout block (pad, 5' handle, 20-N
slot, RT site, 10x CS1, Esp3I site); and synthesis windows (warning only).

**Across the pool.**

* ``cassette_not_shared``: no two inserts carry the same cassette DNA. Two that do share the
  whole 153-nt cassette, an identical stretch a reverse-transcription template switch can cross
  without a trace. Gated.
* ``cassette_region_tract``: the longest stretch of identical sequence each insert shares with
  any insert of another member, within its cassette ± 150 nt (cut at the barcode slot).
  Reported per insert with the partner; a warning above the 100-nt swap gate when the partner is
  another gene. Computed exactly, by a per-insert binary search over the length with a
  rolling-hash index of every window.
* ``completeness``: every member ``l1_completeness.tsv`` calls built has exactly its arms in the
  pool, and every insert belongs to a built member; arm ids are unique.

Reads committed files only (the five rung tables and FASTAs, ``l1_completeness.tsv``, the
committed natives), so it runs offline. ``--order FASTA`` checks a different insert file, such
as the Twist-recoded order of record, matched to the tables by arm id.

Writes ``l1_pool_qc.tsv`` (one row per insert, one column per check: ``pass``, ``FAIL``,
``warn`` or ``skip``, with the details of everything not passing) and
``l1_pool_qc_summary.json``. Exits 1 if any error-level check fails.

Implemented by Claude Code session l1-review-export-scope (fork for pool QC), 2026-11-21, at
Becca's request; decided by Becca, 2026-11-21 (quoted above).

    uv run --locked python l1_pool_qc.py
    uv run --locked python l1_pool_qc.py --order l1_order_inserts_recoded.fasta
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(REPO / "scratch/becca/appris-principal"))

import assemble_candidates as ac  # noqa: E402
import l1_arm_ids as arm_ids  # noqa: E402
import l1_consistency as cons  # noqa: E402

from deliverome_analysis import plasmid_qc as qc  # noqa: E402
from deliverome_analysis import seqtools  # noqa: E402

#: Short rung name -> arm table stem (the five rungs of the order; the retired rescue rung is
#: kept as a record and not ordered).
RUNGS = {
    "cterm": "constructs_l1_cterm",
    "nterm": "constructs_l1_nterm",
    "multipass": "constructs_l1_multipass",
    "type3": "constructs_l1_type3",
    "internal": "constructs_l1_internal",
}
BACKBONE_GB = HERE / "p0039_EFS_IRES_puro_backbone.gb"
ARM_ENCODINGS = HERE / "l1_arm_cassette_encodings.tsv"
CTERM_ENCODINGS = HERE / "l1_cassette_encodings.tsv"
COMPLETENESS = HERE / "l1_completeness.tsv"
OUT_TSV = HERE / "l1_pool_qc.tsv"
OUT_JSON = HERE / "l1_pool_qc_summary.json"

#: Context either side of the cassette for the shared-tract measure.
WINDOW_FLANK = 150
#: The swap gate (``l1_swap_risk.GATE``, Walton's fit): a tract above it is a warning.
SWAP_GATE = 100
BARCODE_SLOT = "N" * qc.BARCODE_LEN
READOUT_BLOCK = qc.SPACER + qc.FWD_PRIMER + BARCODE_SLOT + ac.REV_PRIMER + ac.CS1 + ac.ESP3I_3
OPENING = ac.ESP3I_5 + qc.KOZAK + "ATG"
#: Backbone feature types worth checking for presence once (primer sites repeat by design).
FEATURE_TYPES = {"promoter", "CDS", "rep_origin", "LTR", "misc_feature"}
MIN_FEATURE_NT = 100

PASS, FAIL, WARN, SKIP = "pass", "FAIL", "warn", "skip"


# --------------------------------------------------------------------------- inputs


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def sha(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


def read_fasta(path: Path) -> list[tuple[str, str]]:
    out: list[tuple[str, str]] = []
    name, buf = None, []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                out.append((name, "".join(buf).upper()))
            name, buf = line[1:], []
        else:
            buf.append(line.strip())
    if name is not None:
        out.append((name, "".join(buf).upper()))
    return out


def load_inserts(order: Path | None = None) -> tuple[list[dict], list[str]]:
    """Every insert of the five rungs, with its table row and arm id.

    With ``order``, each insert's DNA is taken from that file instead, matched by arm id (the
    first header token); the second value lists arm ids present on one side only.
    """
    decisions = arm_ids.load_decisions()
    out = []
    for rung, stem in RUNGS.items():
        rows = read_tsv(HERE / f"{stem}.tsv")
        key = "construct_id" if rows and "construct_id" in rows[0] else "candidate_id"
        by_id = {r[key]: r for r in rows}
        by_id.update({r["arm_id"]: r for r in rows if r.get("arm_id")})
        for header, seq in read_fasta(HERE / f"{stem}.fasta"):
            ident = header.split("|")[0].split()[0]
            row = by_id[ident]
            design = row.get("design_id", "")
            out.append(
                {
                    "rung": rung,
                    "arm_id": arm_ids.arm_id_for(rung, row, design, decisions),
                    "design_id": design,
                    "candidate_id": row.get("candidate_id", ""),
                    "gene": row.get("gene_symbol", ""),
                    "row": row,
                    "insert": seq,
                }
            )
    arm_ids.check_unique(i["arm_id"] for i in out)
    unmatched: list[str] = []
    if order is not None:
        ordered = {h.split()[0]: s for h, s in read_fasta(order)}
        for i in out:
            if i["arm_id"] in ordered:
                i["insert"] = ordered[i["arm_id"]]
            else:
                unmatched.append(f"{i['arm_id']} not in the order")
        known = {i["arm_id"] for i in out}
        unmatched += [f"{a} in the order but in no arm table" for a in ordered if a not in known]
        out = [i for i in out if i["arm_id"] in ordered]
    return out, unmatched


def issued_encodings() -> tuple[dict[str, str], dict[str, str], str]:
    """``(by arm_id, by C-terminal candidate_id, source)`` of issued cassette encodings."""
    by_arm: dict[str, str] = {}
    if ARM_ENCODINGS.exists():
        by_arm = {r["arm_id"]: r["cassette_nt"] for r in read_tsv(ARM_ENCODINGS)}
    by_cid = {r["candidate_id"]: r["cassette_nt"] for r in read_tsv(CTERM_ENCODINGS)}
    source = ARM_ENCODINGS.name if by_arm else f"{CTERM_ENCODINGS.name} only (no per-arm table)"
    return by_arm, by_cid, source


def backbone() -> tuple[str, dict[str, str]]:
    """The empty p0039 and the features that occur exactly once in it."""
    from Bio import SeqIO

    rec = SeqIO.read(str(BACKBONE_GB), "genbank")
    bb = str(rec.seq).upper()
    feats: dict[str, str] = {}
    for f in rec.features:
        label = (f.qualifiers.get("label") or [""])[0]
        if f.type not in FEATURE_TYPES or not label or label in feats:
            continue
        lo, hi = int(f.location.start), int(f.location.end)
        if lo < ac.SLOT_START + ac.SLOT_LEN and hi > ac.SLOT_START:
            continue  # the MCS + stuffer the insert replaces
        s = str(f.extract(rec.seq)).upper()
        if len(s) < MIN_FEATURE_NT:
            continue
        n = qc._circular_count(bb, s) + qc._circular_count(bb, seqtools.reverse_complement(s))
        if n == 1:
            feats[label] = s
    return bb, feats


def assemble(bb: str, insert: str) -> str:
    """The finished plasmid: the wrapped insert in place of p0039's MCS + stuffer slot."""
    return bb[: ac.SLOT_START] + insert + bb[ac.SLOT_START + ac.SLOT_LEN :]


# --------------------------------------------------------------------------- per insert


def is_plain(ins: dict) -> bool:
    """A plain insertion: ``native[:n] + cassette + native[n:]`` (no leader swap, no snorkel)."""
    row = ins["row"]
    if (row.get("igk_replaces_residues") or "").strip() or arm_ids.is_snorkel(row):
        return False
    design = " ".join((row.get(c) or "") for c in ("placement", "design", "arm")).lower()
    return "igk" not in design and "snorkel" not in design


def cassette_site(ins: dict) -> int | None:
    """The cassette's 0-based residue in a plain insertion, from the arm table."""
    row = ins["row"]
    if (row.get("cassette_after_residue") or "").strip():
        return int(float(row["cassette_after_residue"]))
    if ins["rung"] in ("multipass", "type3"):
        return 1  # after Met1
    return None


def spec_for(
    ins: dict,
    bb: str,
    feats: dict[str, str],
    natives: dict[str, str],
    by_arm: dict[str, str],
    by_cid: dict[str, str],
) -> qc.Spec:
    # The per-arm table keys a single-arm rung (nterm, multipass, type3) by its design id and
    # an internal arm by its arm id; try both, then the row's own id.
    cassette = next(
        (
            by_arm[k]
            for k in (ins["arm_id"], ins["design_id"], ins["row"].get("construct_id", ""))
            if k in by_arm
        ),
        None,
    )
    if cassette is None and ins["rung"] == "cterm":
        cassette = by_cid.get(ins["candidate_id"])
    spec = qc.Spec(
        plasmid=ins["arm_id"],
        gene=ins["gene"],
        designed_insert=ins["insert"],
        backbone=bb,
        barcode=BARCODE_SLOT,
        cassette_nt=cassette or qc.P003_CASSETTE_NT,
        cassette_peptide=qc.P003_CASSETTE_PEPTIDE,
        require_p003_tag_codons=True,
        backbone_features=feats,
    )
    native = natives.get(ins["candidate_id"], "")
    n = cassette_site(ins)
    if native and is_plain(ins) and n is not None:
        spec.proteoforms = {ins["design_id"]: native}
        spec.cassette_after_residue = n
        spec.tagged_protein = native[:n] + qc.P003_CASSETTE_PEPTIDE + native[n:]
    return spec


def extra_findings(ins: dict, native: str) -> list[qc.Finding]:
    """Checks that need the arm table or the pool's fixed layout."""
    insert, row = ins["insert"], ins["row"]
    out = []
    orf = qc.read_orf(insert + "A" * 3)  # linear insert; padding keeps read_orf from wrapping
    protein = seqtools.translate(orf[1][:-3]) if orf else ""
    want = {row.get(c, "") for c in ("protein_sha256", "mature_sha256")} - {""}
    if not want:
        out.append(qc.Finding("protein_matches_table", True, "no protein hash in table", "warning"))
    else:
        out.append(
            qc.Finding(
                "protein_matches_table",
                sha(protein) in want,
                "encoded protein hashes to the arm table's"
                if sha(protein) in want
                else "encoded protein differs from the arm table's protein_sha256",
            )
        )
    if not native:
        out.append(qc.Finding("layout", True, "no committed native — not checked", "warning"))
    else:
        verdict = cons.layout_check(row, protein, native)
        out.append(
            qc.Finding(
                "layout",
                verdict.startswith("ok"),
                verdict,
                "warning" if verdict == "n/a" else "error",
            )
        )
    ok_open = insert.startswith(OPENING)
    tail_ok = False
    if orf is not None:
        end = insert.find(orf[1]) + len(orf[1])
        tail_ok = insert[end:] == READOUT_BLOCK
    detail = []
    if not ok_open:
        detail.append("does not open with Esp3I + Kozak + ATG")
    if not tail_ok:
        detail.append("after the stop is not the pool's readout block (pad, handle, 20 N, RT, CS1)")
    out.append(qc.Finding("insert_ends", ok_open and tail_ok, "; ".join(detail) or "as the pool's"))
    out.append(qc.check_synthesis_windows(insert, qc.Spec(plasmid="", designed_insert=insert)))
    return out


def cassette_window(insert: str) -> tuple[str, str] | None:
    """``(cassette DNA, cassette ± WINDOW_FLANK nt cut at the barcode slot)``."""
    orf = qc.read_orf(insert + "AAA")
    if orf is None:
        return None
    protein = seqtools.translate(orf[1][:-3])
    if protein.count(qc.P003_CASSETTE_PEPTIDE) != 1:
        return None
    s = orf[0] + 3 * protein.find(qc.P003_CASSETTE_PEPTIDE)
    e = s + len(qc.P003_CASSETTE_NT)
    window = insert[max(0, s - WINDOW_FLANK) : e + WINDOW_FLANK]
    return insert[s:e], window.split("N")[0]


# --------------------------------------------------------------------------- across the pool


def cassette_sharing(cassette: str, holders: Counter) -> qc.Finding:
    """Gated: an insert whose cassette DNA another insert also carries shares the whole cassette."""
    if not cassette:
        return qc.Finding("cassette_not_shared", False, "no 51-aa cassette found")
    n = holders[cassette] - 1
    if n == 0:
        return qc.Finding("cassette_not_shared", True, "own cassette DNA")
    return qc.Finding(
        "cassette_not_shared",
        False,
        f"shares its whole {len(cassette)}-nt cassette DNA with {n} other insert(s)"
        + (" (p003's shared encoding)" if cassette == qc.P003_CASSETTE_NT else ""),
    )


_BASE = np.uint64(0x9E3779B97F4A7C15)  # odd 64-bit multiplier; arithmetic wraps mod 2**64
_CODE = {"A": 1, "C": 2, "G": 3, "T": 4}


class TractIndex:
    """Polynomial hashes of every substring of every window, fully vectorized.

    With ``y[k] = x[k] * B**k`` and ``S`` its running sum, ``S[p+L] - S[p]`` is a substring's
    hash times ``B**p``; multiplying by ``B**(N-p)`` puts every substring on the same footing,
    so equal substrings get equal values wherever they sit. All of it wraps mod 2**64.
    """

    def __init__(self, windows: list[str], owners: list[int]) -> None:
        self.windows = windows
        lens = np.array([len(w) for w in windows], dtype=np.int64)
        self.lens = lens
        self.starts = np.concatenate([[0], np.cumsum(lens + 1)[:-1]]).astype(np.int64)
        text = "\x00".join(windows) + "\x00"  # a separator no window contains
        codes = np.frombuffer(text.encode(), dtype=np.uint8)
        x = np.zeros(len(codes), dtype=np.uint64)
        for ch, v in _CODE.items():
            x[codes == ord(ch)] = v
        n = len(codes)
        with np.errstate(over="ignore"):
            self.pw = np.cumprod(np.full(n + 1, _BASE, dtype=np.uint64))
            self.pw = np.concatenate([[np.uint64(1)], self.pw[:-1]])  # pw[k] = B**k
            self.cum = np.concatenate([[np.uint64(0)], np.cumsum(x * self.pw[:n])])
        self.n = n
        self.window_at = np.repeat(np.arange(len(windows), dtype=np.int64), lens + 1)
        self.owner_at = np.repeat(np.asarray(owners, dtype=np.int64), lens + 1)
        self.pos = np.arange(n, dtype=np.int64) - np.repeat(self.starts, lens + 1)

    def table(self, length: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """``(hash, window, owner, start)`` of every in-window L-mer, sorted by hash, owner."""
        p = np.arange(self.n - length + 1, dtype=np.int64)
        p = p[self.pos[p] + length <= self.lens[self.window_at[p]]]
        with np.errstate(over="ignore"):
            h = (self.cum[p + length] - self.cum[p]) * self.pw[self.n - p]
        own = self.owner_at[p]
        order = np.lexsort((own, h))
        return h[order], self.window_at[p][order], own[order], p[order]

    @staticmethod
    def groups(h: np.ndarray, own: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Per element: its hash group's first and last index (groups sorted by owner)."""
        edge = np.concatenate([[True], h[1:] != h[:-1]])
        gid = np.cumsum(edge) - 1
        first = np.nonzero(edge)[0]
        last = np.concatenate([first[1:], [len(h)]]) - 1
        return first[gid], last[gid]

    def shared(self, length: int) -> np.ndarray:
        """Per window: a window of a different owner sharing one of its L-mers, or -1."""
        h, win, own, _ = self.table(length)
        out = np.full(len(self.windows), -1, dtype=np.int64)
        if len(h):
            f, last = self.groups(h, own)
            multi = own[f] != own[last]
            other = np.where(own[f] != own, win[f], win[last])
            out[win[multi]] = other[multi]
        return out

    def confirm(self, i: int, j: int, length: int) -> bool:
        """Window j really holds a length-L substring of window i (no hash collision)."""
        wi, wj = self.windows[i], self.windows[j]
        return any(wi[q : q + length] in wj for q in range(len(wi) - length + 1))


#: Coarse ladder for the first pass: every window is bracketed to within this many nt by one
#: index build per rung, then binary-searched inside its bracket.
LADDER_STEP = 16


def longest_shared(windows: list[str], owners: list[int]) -> tuple[list[int], list[int | None]]:
    """Per window, the longest substring it shares with a window of another owner, and that
    window. A ladder pass brackets every window at once, then a parallel binary search finishes
    inside the brackets: each round groups windows by trial length, so a length's index is
    built once for every window trying it (about 100 builds for the L1 pool, not one per
    window)."""
    idx = TractIndex(windows, owners)
    n = len(windows)
    lo = np.zeros(n, dtype=np.int64)  # shared at lo (0 trivially)
    hi = np.array([len(w) + 1 for w in windows], dtype=np.int64)  # not shared at hi
    partner = np.full(n, -1, dtype=np.int64)

    def trial(length: int, who: np.ndarray) -> None:
        got = idx.shared(length)[who]
        ok = got >= 0
        lo[who[ok]] = length
        partner[who[ok]] = got[ok]
        hi[who[~ok]] = np.minimum(hi[who[~ok]], length)

    for length in range(LADDER_STEP, int(hi.max()), LADDER_STEP):
        who = np.nonzero((lo < length) & (length < hi))[0]
        if len(who):
            trial(length, who)
    while True:
        active = np.nonzero(hi - lo > 1)[0]
        if len(active) == 0:
            break
        mid = (lo[active] + hi[active]) // 2
        for m in np.unique(mid):
            trial(int(m), active[mid == m])
    out: list[int | None] = []
    for i in range(n):
        j = int(partner[i])
        if lo[i] > 0 and j >= 0 and not idx.confirm(i, j, int(lo[i])):
            raise AssertionError(f"hash collision between windows {i} and {j}")
        out.append(j if lo[i] > 0 and j >= 0 else None)
    return lo.tolist(), out


# --------------------------------------------------------------------------- run


def run(order: Path | None = None) -> tuple[list[dict], dict]:
    t0 = time.time()
    inserts, unmatched = load_inserts(order)
    bb, feats = backbone()
    by_arm, by_cid, enc_source = issued_encodings()
    natives = cons.current_natives(needed={i["candidate_id"] for i in inserts})

    rows: list[dict] = []
    findings_of: dict[str, list[qc.Finding]] = {}
    for ins in inserts:
        spec = spec_for(ins, bb, feats, natives, by_arm, by_cid)
        found = qc.check_plasmid(assemble(bb, ins["insert"]), spec)
        found += extra_findings(ins, natives.get(ins["candidate_id"], ""))
        findings_of[ins["arm_id"]] = found

    # cassette_not_shared
    cas_of: dict[str, str] = {}
    windows: list[str] = []
    for ins in inserts:
        cw = cassette_window(ins["insert"])
        cas_of[ins["arm_id"]] = cw[0] if cw else ""
        windows.append(cw[1] if cw else "")
    holders = Counter(c for c in cas_of.values() if c)
    for ins in inserts:
        findings_of[ins["arm_id"]].append(cassette_sharing(cas_of[ins["arm_id"]], holders))

    # cassette_region_tract
    designs = sorted({i["design_id"] for i in inserts})
    owner = {d: k for k, d in enumerate(designs)}
    t1 = time.time()
    tract, partner = longest_shared(windows, [owner[i["design_id"]] for i in inserts])
    tract_s = time.time() - t1
    for k, ins in enumerate(inserts):
        p = partner[k]
        other = inserts[p] if p is not None else None
        same_gene = bool(other) and other["gene"] == ins["gene"]
        over = tract[k] > SWAP_GATE and not same_gene
        ins["tract"], ins["partner"] = tract[k], (other["arm_id"] if other else "")
        findings_of[ins["arm_id"]].append(
            qc.Finding(
                "cassette_region_tract",
                not over,
                f"{tract[k]} nt shared with {ins['partner'] or 'no other insert'}"
                + (" (same gene)" if same_gene else ""),
                "warning",
            )
        )

    # pool-level: completeness
    pool: list[qc.Finding] = []
    built = {r["design_id"]: r for r in read_tsv(COMPLETENESS) if r["status"] == "built"}
    arms_of = Counter(i["design_id"] for i in inserts)
    missing = sorted(d for d in built if d not in arms_of)
    stray = sorted(d for d in arms_of if d not in built)
    wrong_n = sorted(
        d
        for d, r in built.items()
        if d in arms_of and (r.get("arms") or "").isdigit() and int(r["arms"]) != arms_of[d]
    )
    pool.append(
        qc.Finding(
            "completeness",
            not (missing or stray or wrong_n or unmatched),
            f"{len(built)} built members, {len(inserts)} inserts; missing {missing[:5]} "
            f"({len(missing)}), not a built member {stray[:5]} ({len(stray)}), arm count off "
            f"{wrong_n[:5]} ({len(wrong_n)}), order vs tables {unmatched[:5]} ({len(unmatched)})",
        )
    )

    checks: list[str] = []
    for found in findings_of.values():
        for f in found:
            if f.check not in checks:
                checks.append(f.check)
    counts: dict[str, Counter] = defaultdict(Counter)
    for ins in inserts:
        found = {f.check: f for f in findings_of[ins["arm_id"]]}
        row = {
            "arm_id": ins["arm_id"],
            "design_id": ins["design_id"],
            "gene_symbol": ins["gene"],
            "rung": ins["rung"],
            "worst_partner_tract_nt": ins["tract"],
            "worst_partner": ins["partner"],
        }
        problems = []
        for c in checks:
            f = found.get(c)
            if f is None:
                status = SKIP
            elif f.ok and f.severity == "warning" and "not checked" in f.detail:
                status = SKIP
            elif f.ok:
                status = PASS
            elif f.severity == "warning":
                status = WARN
            else:
                status = FAIL
            row[c] = status
            counts[c][status] += 1
            if status in (FAIL, WARN) and f is not None:
                problems.append(f"{c}: {f.detail}")
        row["problems"] = " | ".join(problems)
        rows.append(row)

    t = np.array(tract)
    failing = sorted(c for c in checks if counts[c][FAIL]) + [
        f.check for f in pool if not f.ok and f.severity == "error"
    ]
    summary = {
        "inserts": len(inserts),
        "inserts_source": str(order) if order else "committed constructs_l1_*.fasta",
        "inserts_sha256": sha("".join(i["insert"] for i in inserts)),
        "cassette_encodings": enc_source,
        "checks": {c: dict(counts[c]) for c in checks},
        "pool_checks": {f.check: {"ok": f.ok, "detail": f.detail} for f in pool},
        "failing_checks": failing,
        "cassette_dna_variants": len(holders),
        "largest_shared_cassette_group": max(holders.values()) if holders else 0,
        "cassette_region_tract_nt": {
            "median": float(np.median(t)) if len(t) else None,
            "p90": float(np.percentile(t, 90)) if len(t) else None,
            "max": int(t.max()) if len(t) else None,
            "max_insert": inserts[int(t.argmax())]["arm_id"] if len(t) else None,
            "over_swap_gate_other_gene": counts["cassette_region_tract"][WARN],
            "window": f"cassette ± {WINDOW_FLANK} nt, cut at the barcode slot",
        },
        "backbone_features_checked": sorted(feats),
        "runtime_s": {"total": round(time.time() - t0, 1), "tract": round(tract_s, 1)},
        "decided_by": 'Becca, 2026-11-21 ("we want the plasmid qc to run on the pool")',
        "implemented_by": "Claude Code session l1-review-export-scope (pool QC fork)",
    }
    return rows, summary


def write(rows: list[dict], summary: dict) -> None:
    with OUT_TSV.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    OUT_JSON.write_text(json.dumps(summary, indent=2) + "\n")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--order", type=Path, help="insert FASTA to check instead of the arm tables'")
    ap.add_argument("--no-write", action="store_true", help="report only")
    a = ap.parse_args(argv)
    rows, summary = run(a.order.expanduser() if a.order else None)
    if not a.no_write:
        write(rows, summary)
    print(json.dumps({k: summary[k] for k in ("inserts", "checks", "failing_checks")}, indent=1))
    print(json.dumps(summary["cassette_region_tract_nt"]), summary["runtime_s"])
    return 1 if summary["failing_checks"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
