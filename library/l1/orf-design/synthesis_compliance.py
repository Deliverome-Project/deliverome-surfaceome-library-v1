"""Synthesis compliance for the final tagged L1 inserts: a local check and minimal synonymous recodes.

Twist scored the first L1 draft (untagged ORFs; ``twist_complexity_map.py``) and refused 25
inserts. This stage checks the inserts we will actually order -- each rung's tagged arm as
built, wrapped as the vendor receives it (Esp3I, Kozak, tagged ORF, readout block with the
20-nt barcode as ``N``, Esp3I) -- and makes the fewest synonymous codon changes that clear
each problem.

**The checker approximates Twist; it is not Twist's algorithm.** Its rules were calibrated
against Twist's calls on the 2,735 draft inserts, rebuilt here exactly (``--calibrate``):

| Rule | Local check | Agreement with Twist on the draft |
|---|---|---|
| Tn5 mosaic end (impermissible) | edit distance <= 3 between the 19-nt ME (either strand) and any substring | 13 of 13 refused at 3; every accepted insert at >= 4 |
| 100-bp windowed GC (refused) | any 100-bp window above 90% GC | Twist refused 11 with windows of 91-98%, but accepted some at 91-92%, so this errs strict |
| 50-bp windowed GC (flag) | any 50-bp window above 80% GC | all 568 Twist flags, plus 134 Twist did not flag (stricter) |
| Direct repeat (refused) | a direct repeat of >= 200 nt at >= 90% identity, tandem copies included | the only refusal is CLEC4M (204 nt by Twist, 433 nt here). Twist's rule is not reproducible: CD209 (same 69-nt offset, 453 nt) and MUC21 (984 nt) scored High, so this errs strict |

Targets after recoding, per insert: every 50-bp window at or below 80% GC (which keeps
every 100-bp window inside spec) and at or below 75% where the window spans the cassette's
last base (margin: the cassette encodings are fixed, so a vendor flag there would mean
reissuing them; reported, not required), ME edit distance >= 4, no exact direct repeat of 40 nt or
longer (the length that drives Twist's High repeat scores), and no >= 90%-identity repeat
of 100 nt or longer. ``N`` barcode bases count as 0.5 GC and take part in no repeat or ME match.

**What may change:** only codons of the member's own protein, each to a synonymous codon
with human usage >= 10% (``seqtools.synonymous_codons``). **What may not:** Met1; the
51-aa cassette (its linker codons are the per-member encodings the swap-risk check
governs) and one codon either side of it; an Igκ leader; and the native signal peptide
(``l1_tag_paths.tsv``'s DeepTMHMM v1 ``signal_peptide_aa``) plus 3 codons after its
cleavage site. A candidate edit is rejected if it creates a BbsI, Esp3I or BsaI site, a
consensus splice donor (``GT[AG]AGT``) or a new ME match, or if the protein changes; the
whole insert is re-checked after recoding, and the cassette DNA, wrapper and readout block
are asserted byte-identical. Edits are greedy and local: the worst window first, the
codon change that helps most, ties to the more frequent codon.

Inputs are each rung's ``constructs_l1_<rung>.fasta`` and ``.tsv`` at the commit in
:data:`RUNG_REFS`: the integration commit :data:`INPUT_REF` for the five tag-placement rungs.
(The 25 expression-rescue arms were retired on 2026-11-13, one design per gene in one pool, so
no rescue rung is read; the 2026-11-13 order of record that held them is kept, superseded.)
They are read with ``git show`` and
checksummed. The recoder is deterministic (no random choice, ties broken by codon usage and
then position), so ``--check`` reproduces every table byte for byte.

Outputs: ``synthesis_compliance.tsv`` (one row per insert: before/after metrics,
Twist-like category, status), ``synthesis_recodes.tsv`` (one row per codon changed, with
``decided_by``/``implemented_by``), ``synthesis_compliance_summary.json``, and the recoded
inserts as FASTA on S3 (too large for git; its SHA-256 is in the summary).

Requested by the coordinator session under Becca's 2026-11-12 delegation, after Twist's
complexity report; implemented by the Claude Code session l1-reviewer-package, PR #283.

    uv run --locked python synthesis_compliance.py              # check + recode, write tables
    uv run --locked python synthesis_compliance.py --calibrate  # rules vs Twist on the draft
    uv run --locked python synthesis_compliance.py --publish    # recoded FASTA to S3
    uv run --locked python synthesis_compliance.py --check      # committed tables reproduce
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import sys
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(HERE))
import l1_arm_ids as arm_ids  # noqa: E402

from deliverome_analysis import plasmid_qc as qc  # noqa: E402
from deliverome_analysis import seqtools  # noqa: E402
from deliverome_analysis.plasmid_qc import P003_CASSETTE_PEPTIDE  # noqa: E402

#: Integration commit (PR #299) whose rung builds this stage checks. Bump to re-run on a newer build.
#: 496d3482 (2026-11-13) carries #312: the 23 rescue genes' GPCR references built Igκ + tag and
#: SLC61A1 x2 built; the 25 rescue members retired (Becca, 2026-11-13: one design per gene, one
#: pool). Bumped from 355b524d by l1-reviewer-package at the coordinator's request.
#: 2026-11-14, final regeneration: main at cdb78068 (#336 C-terminal stalk holds, snorkel list
#: 2026-11-14b and the SLC68A1 arm; #331/#333 internal redesign; #324 disposition gate). A
#: commit on main, so main alone rebuilds the order (pinned_refs stage). 496d3482 was the
#: draft integration branch and is superseded.
INPUT_REF = "e3336caf6b0e4c084cbfeec9e5d189015794284e"
RUNG_REFS = {
    "nterm": INPUT_REF,
    "multipass": INPUT_REF,
    "type3": INPUT_REF,
    "internal": INPUT_REF,
    "cterm": INPUT_REF,
}
RUNGS = tuple(RUNG_REFS)
REL = "library/l1/orf-design"
TAG_PATHS = HERE / "l1_tag_paths.tsv"
OUT = HERE / "synthesis_compliance.tsv"
INSERT_NUMBERING = HERE / "l1_insert_numbering.tsv"
RECODES = HERE / "synthesis_recodes.tsv"
SUMMARY = HERE / "synthesis_compliance_summary.json"
#: The order of record. Versioned by date: a new run writes a new key, never over an old one.
FASTA_KEY = "sprint4/s4e9-orf-design-fg-panel/synthesis/2026-11-23-pre-freeze/l1_order_inserts_recoded.fasta"
#: Superseded orders of record, kept unchanged: 2026-11-21-arm-linkers (sha256 8b1c4ed3…, the
#: same inserts with the 109 snorkel arms on one shared module encoding); 2026-11-15-final-e (sha256 20e8f275…, the same
#: 3,255 inserts with the 2,942 upstream arms on p003's one shared cassette encoding); 2026-11-13 (sha256 d1257f0a…, 2,716 inserts
#: including the 25 retired rescue arms and 23 tag-alone GPCR references); 2026-11-13-one-arm
#: (sha256 1f29936b…, 2,678 inserts, the current builds but named by pre-renumber member ids,
#: and uploaded from an uncommitted tree); 2026-11-13-one-arm-renamed (sha256 b426e433…, 2,678
#: inserts from the draft integration build 496d3482, which is not on main).
BARCODE_BLOCK = HERE / "l1_barcode_block_provisional.tsv"
BARCODE_EXTENSION = HERE / "l1_barcode_block_extension_2026-11-15.tsv"

ME = "CTGTCTCTTATACACATCT"  # Tn5 mosaic end
ME_MAX_REFUSED = 3
GC50_MAX = 0.80  # Medium-complexity threshold (reported; not recoded since 2026-11-23)
GC50_HIGH = 0.90  # High-complexity threshold: only High and Not Accepted inserts are recoded (Becca, 2026-11-23)
CASSETTE_END_GC50 = 0.75  # reported metric for windows spanning the cassette's last base
GC100_REFUSED = 0.90  # Twist refused 0.91-0.98; set below its refusals, so it errs strict
REPEAT_TARGET = 40  # no exact direct repeat this long after recoding
REPEAT_REFUSED = 200  # near-repeat length; Twist's own repeat rule could not be reproduced
NEAR_TARGET = 100  # no >=90%-identity direct repeat this long after recoding
MIN_USAGE = 0.10
IGK_LEADER = qc.IGK_LEADER_AA + qc.IGK_PLUS1_AA
SP_MARGIN = 3
DONOR = re.compile("GT[AG]AGT")
FIVE = qc.ESP3I_5 + qc.KOZAK
DECIDED_BY = 'l1-reviewer-package (delegated by Becca 2026-11-12, "use your best judgment")'
IMPLEMENTED_BY = "Claude Code session l1-reviewer-package, PR #283 (synthesis_compliance.py)"


rc = seqtools.reverse_complement


def me_distance(s: str) -> int:
    """Smallest edit distance between the ME (either strand) and any substring of ``s``."""
    best = len(ME)
    for t in (s, rc(s)):
        prev = list(range(len(ME) + 1))
        for c in t:
            cur = [0]
            for i, pc in enumerate(ME, 1):
                cur.append(min(prev[i] + 1, cur[i - 1] + 1, prev[i - 1] + (pc != c or c == "N")))
            best = min(best, cur[-1])
            prev = cur
    return best


def gc_profile(s: str, w: int) -> list[float]:
    val = [1.0 if c in "GC" else 0.5 if c == "N" else 0.0 for c in s]
    if len(s) < w:
        return [sum(val) / max(len(s), 1)]
    cur = sum(val[:w])
    out = [cur / w]
    for i in range(w, len(s)):
        cur += val[i] - val[i - w]
        out.append(cur / w)
    return out


def longest_repeat(s: str, k: int = 12) -> tuple[int, int, int]:
    """(length, i, j) of the longest exact direct repeat, k-mer seeded; N never matches."""
    seen: dict[str, list[int]] = {}
    for i in range(len(s) - k + 1):
        km = s[i : i + k]
        if "N" not in km:
            seen.setdefault(km, []).append(i)
    best = (0, 0, 0)
    for pos in seen.values():
        if len(pos) < 2:
            continue
        for a in range(len(pos)):
            for b in range(a + 1, len(pos)):
                i, j = pos[a], pos[b]
                if i > 0 and s[i - 1] == s[j - 1] and s[i - 1] != "N":
                    continue  # not a maximal start
                n = k
                while j + n < len(s) and i + n < j and s[i + n] == s[j + n] and s[j + n] != "N":
                    n += 1
                if n > best[0]:
                    best = (n, i, j)
    return best


def near_repeat(s: str, k: int = 12, max_mismatch: float = 0.10) -> tuple[int, int, int]:
    """(length, i, j): the longest direct repeat at >= 90% identity (ungapped), 12-mer seeded.

    Twist refused CLEC4M for a 204-nt repeat that is only 135 nt as an exact match: its
    neck repeats differ at a few bases. Copies may overlap (tandem repeats). Each diagonal (offset ``d = j - i``) that shares a
    12-mer is scanned for the longest segment whose mismatches are at most 10% of its length.
    """
    seen: dict[str, list[int]] = {}
    for i in range(len(s) - k + 1):
        km = s[i : i + k]
        if "N" not in km:
            seen.setdefault(km, []).append(i)
    offsets = {b - a for pos in seen.values() for a in pos for b in pos if b > a}
    best = (0, 0, 0)
    for d in offsets:
        m = [s[i] == s[i + d] and s[i] != "N" for i in range(len(s) - d)]
        lo, mism = 0, 0
        for hi in range(len(m)):
            mism += not m[hi]
            while lo < hi and (mism * 1.0 > max_mismatch * (hi - lo + 1) or not m[lo]):
                mism -= not m[lo]
                lo += 1
            # tandem copies may overlap (CLEC4M's pair is offset by 69 nt over 204 nt)
            if m[hi] and hi - lo + 1 > best[0]:
                best = (hi - lo + 1, lo, lo + d)
    return best


def metrics(s: str) -> dict:
    rep = longest_repeat(s)
    return {
        "max_gc50": round(max(gc_profile(s, 50)), 4),
        "max_gc100": round(max(gc_profile(s, 100)), 4),
        "me_distance": me_distance(s),
        "longest_repeat_nt": rep[0],
        "near_repeat_nt": near_repeat(s)[0],
    }


def category(m: dict) -> str:
    """A Twist-like score from the local metrics (an approximation)."""
    if (
        m["me_distance"] <= ME_MAX_REFUSED
        or m["max_gc100"] > GC100_REFUSED
        or m["near_repeat_nt"] >= REPEAT_REFUSED
    ):
        return "Not Accepted"
    if (
        m["longest_repeat_nt"] >= REPEAT_TARGET
        or m["near_repeat_nt"] >= NEAR_TARGET
        or m["max_gc50"] >= 0.90
    ):
        return "High"
    if m["max_gc50"] > GC50_MAX or m["longest_repeat_nt"] >= 18:
        return "Medium"
    return "Low"


def compliant(m: dict) -> bool:
    """True unless the insert scores High or Not Accepted (Becca, 2026-11-23: leave Low/Medium as-is)."""
    return (
        m["max_gc50"] < GC50_HIGH
        and m["max_gc100"] <= GC100_REFUSED
        and m["me_distance"] > ME_MAX_REFUSED
        and m["longest_repeat_nt"] < REPEAT_TARGET
        and m["near_repeat_nt"] < NEAR_TARGET
    )


# --------------------------------------------------------------------------- inputs


def git_show(path: str, ref: str = INPUT_REF) -> str:
    from deliverome_analysis.provenance import read_pinned_text

    return read_pinned_text(ref, path, repo=REPO)


def split_insert(s: str) -> tuple[str, str, str]:
    """(5' wrapper, ORF incl. stop, 3' readout + wrapper)."""
    assert s.startswith(FIVE), "unexpected 5' wrapper"
    body = s[len(FIVE) :]
    i = 0
    while i + 3 <= len(body) and seqtools.GENETIC_CODE.get(body[i : i + 3]) != "*":
        i += 3
    return FIVE, body[: i + 3], body[i + 3 :]


def current_name(name: str, design_id: str, member_of_design: dict[str, str]) -> str:
    """``name`` with a leading ``l1-mNNNN`` member id replaced by the design's current one."""
    if re.match(r"l1-m\d{4}", name):
        return member_of_design[design_id] + name[len("l1-m0000") :]
    return name


def load_inserts() -> tuple[list[dict], dict[str, str]]:
    sp = {}
    with TAG_PATHS.open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            sp[r.get("design_id") or r["candidate_id"]] = int(r.get("signal_peptide_aa") or 0)
    design_of_member = {}
    with (HERE / "candidate_manifest_tranches.tsv").open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            design_of_member[r["candidate_id"]] = r["design_id"]
    # Rung builds name inserts by the member id current when they were built; member ids follow
    # plate order and move when membership changes, so each insert is renamed to its member id
    # now (l1_member_renumbering.tsv) and its design id is kept beside it. Without this, 2,318
    # of 2,678 insert names disagreed with the plate map on 2026-11-13 (l1-reviewer-package).
    member_of_design = {}
    with (HERE / "l1_member_renumbering.tsv").open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            member_of_design[r["design_id"]] = r["candidate_id"]
    # One row per insert, keyed by its permanent arm id (l1_arm_ids.py); the insert id is its
    # plate position from the last numbering (renumber_library_members.py), blank for an arm
    # not yet numbered (the first pass of a cut, before the numbering it feeds).
    decisions = arm_ids.load_decisions()
    insert_of_arm = {}
    if INSERT_NUMBERING.exists():
        with INSERT_NUMBERING.open() as fh:
            insert_of_arm = {
                r["arm_id"]: r["insert_id"] for r in csv.DictReader(fh, delimiter="\t")
            }
    out, shas = [], {}
    for rung in RUNGS:
        fa_text = git_show(f"{REL}/constructs_l1_{rung}.fasta", RUNG_REFS[rung])
        tsv_text = git_show(f"{REL}/constructs_l1_{rung}.tsv", RUNG_REFS[rung])
        shas[f"constructs_l1_{rung}.fasta"] = hashlib.sha256(fa_text.encode()).hexdigest()
        shas[f"constructs_l1_{rung}.tsv"] = hashlib.sha256(tsv_text.encode()).hexdigest()
        rows = list(csv.DictReader(tsv_text.splitlines(), delimiter="\t"))
        by_sha: dict[str, list[dict]] = {}
        for r in rows:
            by_sha.setdefault(r["protein_sha256"], []).append(r)
        for rec in fa_text.split(">")[1:]:
            head, seq = rec.split("\n", 1)
            seq = seq.replace("\n", "").upper()
            five, orf, three = split_insert(seq)
            prot = seqtools.translate(orf).rstrip("*")
            psha = hashlib.sha256(prot.encode()).hexdigest()
            row = by_sha[psha][0]
            did = row.get("design_id") or design_of_member.get(row["candidate_id"], "")
            if not did.startswith("l1-p"):
                did = design_of_member.get(row["candidate_id"], did)
            arm = arm_ids.arm_id_for(rung, row, did, decisions)
            out.append(
                {
                    "rung": rung,
                    "arm_id": arm,
                    "insert_id": insert_of_arm.get(arm, ""),
                    "construct": current_name(head.split()[0], did, member_of_design),
                    "design_id": did,
                    "gene_symbol": row["gene_symbol"],
                    "five": five,
                    "orf": orf,
                    "three": three,
                    "sp_aa": sp.get(did, 0),
                }
            )
    arm_ids.check_unique(r["arm_id"] for r in out)
    # Order of record: plate order once every insert is numbered, else arm id.
    out.sort(key=lambda r: (r["insert_id"] or "~", r["arm_id"]))
    return out, shas


# --------------------------------------------------------------------------- recoding


def protected_codons(protein: str, sp_aa: int) -> set[int]:
    """0-based codon indices that may not change."""
    keep = {0}
    c = protein.find(P003_CASSETTE_PEPTIDE)
    assert c >= 0 and protein.count(P003_CASSETTE_PEPTIDE) == 1, "cassette not found once"
    keep |= set(range(max(c - 1, 0), c + len(P003_CASSETTE_PEPTIDE) + 1))
    if protein.startswith(IGK_LEADER):
        keep |= set(range(len(IGK_LEADER) + SP_MARGIN))
    if sp_aa and not protein.startswith(IGK_LEADER):
        keep |= set(range(sp_aa + SP_MARGIN))
    return keep


class Recoder:
    def __init__(self, ins: dict, syn: dict[str, list[tuple[str, float]]]):
        self.ins = ins
        self.codons = seqtools.split_codons(ins["orf"])
        self.protein = seqtools.translate(ins["orf"])
        self.keep = protected_codons(self.protein.rstrip("*"), ins["sp_aa"])
        self.syn = syn
        self.off = len(ins["five"])
        self.edits: list[dict] = []
        self.sites0 = len(seqtools.find_sites(self.seq(), list(qc.DOMESTICATION_ENZYMES)))
        self.donors0 = len(DONOR.findall(self.seq()))
        c = self.protein.find(P003_CASSETTE_PEPTIDE)
        self.cas_end = self.off + 3 * (c + len(P003_CASSETTE_PEPTIDE))

    def seq(self) -> str:
        return self.ins["five"] + "".join(self.codons) + self.ins["three"]

    def options(self, ci: int) -> list[tuple[str, float]]:
        aa = seqtools.GENETIC_CODE[self.codons[ci]]
        if ci in self.keep or aa == "*":
            return []
        return [(c, f) for c, f in self.syn[aa] if c != self.codons[ci]]

    def acceptable(self, ci: int, new: str) -> bool:
        old = self.codons[ci]
        self.codons[ci] = new
        s = self.seq()
        ok = (
            len(seqtools.find_sites(s, list(qc.DOMESTICATION_ENZYMES))) <= self.sites0
            and len(DONOR.findall(s)) <= self.donors0
        )
        if ok:
            p = self.off + 3 * ci
            ok = me_distance(s[max(p - 25, 0) : p + 28]) > ME_MAX_REFUSED
        self.codons[ci] = old
        return ok

    def apply(self, ci: int, new: str, reason: str) -> None:
        self.edits.append(
            {
                "residue": ci + 1,
                "amino_acid": seqtools.GENETIC_CODE[new],
                "from_codon": self.codons[ci],
                "to_codon": new,
                "reason": reason,
            }
        )
        self.codons[ci] = new

    def window_codons(self, start: int, end: int) -> list[int]:
        """Codon indices overlapping insert bases [start, end)."""
        lo = max((start - self.off) // 3, 0)
        hi = min((end - self.off + 2) // 3, len(self.codons))
        return list(range(lo, hi))

    def gc_target(self, w: int) -> float:
        """Strictly below 90% (High threshold), so High-GC stretches clear with minimal synonymous edits."""
        return GC50_HIGH - 0.02  # 88% (44/50), strictly < 90% High cutoff

    def fix_gc(self, limit: int = 600) -> None:
        stuck: set[int] = set()
        for _ in range(limit):
            s = self.seq()
            prof = gc_profile(s, 50)
            over = [
                (prof[w] - self.gc_target(w), w)
                for w in range(len(prof))
                if prof[w] >= GC50_HIGH and w not in stuck
            ]
            if not over:
                return
            _excess, worst = max(over)
            best = None
            for ci in self.window_codons(worst, worst + 50):
                cur = self.codons[ci]
                g0 = sum(b in "GC" for b in cur)
                for c, f in self.options(ci):
                    gain = g0 - sum(b in "GC" for b in c)
                    if gain > 0 and (best is None or (gain, f) > (best[0], best[1])):
                        if self.acceptable(ci, c):
                            best = (gain, f, ci, c)
            if best is None:
                stuck.add(worst)
                continue
            self.apply(best[2], best[3], "gc50_high")

    def fix_me(self) -> None:
        for _ in range(10):
            s = self.seq()
            if me_distance(s) > ME_MAX_REFUSED:
                return
            # locate: slide a 40-nt window, find where the distance is minimal
            hit = min(range(0, len(s) - 30, 5), key=lambda i: me_distance(s[i : i + 40]))
            lo = max(hit - 20, 0)
            cands = []
            for ci in self.window_codons(hit, hit + 40):
                for c, f in self.options(ci):
                    old = self.codons[ci]
                    self.codons[ci] = c
                    d = me_distance(self.seq()[lo : hit + 60])  # the match region itself
                    self.codons[ci] = old
                    if d > ME_MAX_REFUSED and self.acceptable(ci, c):
                        cands.append((d, f, ci, c))
            if not cands:
                return
            _d, _f, ci, c = max(cands)
            self.apply(ci, c, "tn5_me")

    def fix_repeats(self, limit: int = 600) -> None:
        """Split the longest repeat until none is over target.

        A repeat is a diagonal: base ``p`` pairs with ``p - d``. Each edit lands at the middle
        of the longest perfectly matching run on that diagonal, in whichever copy has an
        editable codon there, and must make at least one paired base differ (checked against
        both neighbours ``p - d`` and ``p + d``, so a tandem copy is not matched instead). A
        codon is edited at most once, so edits cannot oscillate.
        """
        edited: set[int] = set()
        for _ in range(limit):
            s = self.seq()
            n, i, j = longest_repeat(s)
            if n < REPEAT_TARGET:
                n, i, j = near_repeat(s)
                if n < NEAR_TARGET:
                    return
            d = j - i
            # longest perfect run on the diagonal within the second copy's span
            best_run, cur = (0, j), 0
            for p in range(j, min(j + n, len(s))):
                if s[p] == s[p - d]:
                    cur += 1
                    if cur > best_run[0]:
                        best_run = (cur, p - cur + 1)
                else:
                    cur = 0
            mid = best_run[1] + best_run[0] // 2
            choice = None
            for centre in (mid, mid - d):  # second copy, else the paired base in the first
                cand = sorted(
                    (ci for ci in self.window_codons(centre - 60, centre + 60) if ci not in edited),
                    key=lambda ci: abs(self.off + 3 * ci + 1 - centre),
                )
                for ci in cand:
                    p0 = self.off + 3 * ci
                    for c, _f in sorted(self.options(ci), key=lambda x: -x[1]):
                        breaks = any(
                            c[k] != self.codons[ci][k]
                            and all(
                                not (0 <= q < len(s)) or c[k] != s[q]
                                for q in (p0 + k - d, p0 + k + d)
                            )
                            for k in range(3)
                        )
                        if breaks and self.acceptable(ci, c):
                            choice = (ci, c)
                            break
                    if choice:
                        break
                if choice:
                    break
            if not choice:
                return
            edited.add(choice[0])
            self.apply(choice[0], choice[1], "direct_repeat")

    def run(self) -> str:
        self.fix_me()
        self.fix_repeats()
        self.fix_gc()
        self.fix_me()
        new = self.seq()
        # record the net change per codon (its last reason), not the search history
        orig = seqtools.split_codons(self.ins["orf"])
        last = {e["residue"]: e["reason"] for e in self.edits}
        self.edits = [
            {
                "residue": ci + 1,
                "amino_acid": seqtools.GENETIC_CODE[c],
                "from_codon": orig[ci],
                "to_codon": c,
                "reason": last[ci + 1],
            }
            for ci, c in enumerate(self.codons)
            if c != orig[ci]
        ]
        assert seqtools.translate("".join(self.codons)) == self.protein, "protein changed"
        assert new.startswith(self.ins["five"]) and new.endswith(self.ins["three"])
        return new


def cassette_end_gc(s: str, cas_end: int) -> float:
    """Highest GC of the 50-bp windows spanning the cassette's last base."""
    prof = gc_profile(s, 50)
    lo, hi = max(cas_end - 50, 0), min(cas_end, len(prof))
    return max(prof[lo:hi]) if hi > lo else 0.0


def cassette_dna(orf: str) -> str:
    p = seqtools.translate(orf)
    c = p.find(P003_CASSETTE_PEPTIDE)
    return orf[3 * c : 3 * (c + len(P003_CASSETTE_PEPTIDE))]


# --------------------------------------------------------------------------- main


def unfixable_reason(before: dict, after: dict, rec: Recoder) -> str:
    bits = []
    if after["max_gc50"] >= GC50_HIGH:
        bits.append(
            f"50-bp GC {after['max_gc50']:.0%} remains: the window is protected codons or has no lower-GC synonym"
        )
    if after["max_gc100"] > GC100_REFUSED:
        bits.append(
            f"100-bp GC {after['max_gc100']:.0%} remains: the window is protected codons or has no lower-GC synonym"
        )
    if after["me_distance"] <= ME_MAX_REFUSED:
        bits.append("Tn5 ME-like match remains in protected sequence")
    if after["longest_repeat_nt"] >= REPEAT_TARGET:
        bits.append(
            f"{after['longest_repeat_nt']}-nt exact repeat remains in protected or invariant codons"
        )
    if after["near_repeat_nt"] >= NEAR_TARGET:
        bits.append(f"{after['near_repeat_nt']}-nt near-identical repeat remains")
    return "; ".join(bits)


def tsv_text(data: list[dict]) -> str:
    import io

    buf = io.StringIO()
    w = csv.DictWriter(buf, list(data[0]), delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(data)
    return buf.getvalue()


def run(check: bool = False) -> int:
    inserts, shas = load_inserts()
    syn = seqtools.synonymous_codons(min_usage=MIN_USAGE)
    rows, recodes, fasta = [], [], []
    for ins in inserts:
        s0 = ins["five"] + ins["orf"] + ins["three"]
        before = metrics(s0)
        rec = Recoder(ins, syn)
        cas0 = cassette_end_gc(s0, rec.cas_end)
        s1 = rec.run() if not compliant(before) else s0
        after = metrics(s1) if rec.edits else before
        cas1 = cassette_end_gc(s1, rec.cas_end)
        assert cassette_dna("".join(rec.codons)) == cassette_dna(ins["orf"]), "cassette changed"
        status = (
            "pass"
            if not rec.edits and compliant(before)
            else "fixed"
            if compliant(after)
            else "not_fully_fixable"
        )
        rows.append(
            {
                "rung": ins["rung"],
                "arm_id": ins["arm_id"],
                "insert_id": ins["insert_id"],
                "construct": ins["construct"],
                "design_id": ins["design_id"],
                "gene_symbol": ins["gene_symbol"],
                "insert_nt": len(s0),
                **{f"before_{k}": v for k, v in before.items()},
                "before_category": category(before),
                **{f"after_{k}": v for k, v in after.items()},
                "after_category": category(after),
                "before_cassette_end_gc50": round(cas0, 4),
                "after_cassette_end_gc50": round(cas1, 4),
                "codons_recoded": len(rec.edits),
                "status": status,
                "remaining": unfixable_reason(before, after, rec)
                if status == "not_fully_fixable"
                else "",
                "recoded_insert_sha256": hashlib.sha256(s1.encode()).hexdigest(),
                "decided_by": DECIDED_BY if rec.edits else "",
                "implemented_by": IMPLEMENTED_BY if rec.edits else "",
            }
        )
        for e in rec.edits:
            recodes.append(
                {
                    "rung": ins["rung"],
                    "arm_id": ins["arm_id"],
                    "construct": ins["construct"],
                    "design_id": ins["design_id"],
                    "gene_symbol": ins["gene_symbol"],
                    **e,
                    "decided_by": DECIDED_BY,
                    "implemented_by": IMPLEMENTED_BY,
                }
            )
        # Header agreed with l1-integration (2026-11-15): arm id first, then the plate
        # position the vendor's file names, the member's design id and the rung.
        fasta.append(
            f">{ins['arm_id']} {ins['insert_id'] or '-'} {ins['design_id']} {ins['rung']}\n{s1}\n"
        )
    fasta_text = "".join(fasta)
    # No provisional L1 barcode may occur as a whole 20-mer (either strand) inside any insert.
    barcodes = set()
    for table in (BARCODE_BLOCK, BARCODE_EXTENSION):
        if table.exists():
            with table.open() as fh:
                barcodes |= {r["barcode"].upper() for r in csv.DictReader(fh, delimiter="\t")}
    seen = set()
    for rec in fasta:
        seq = rec.split("\n")[1]
        for t in (seq, rc(seq)):
            seen.update(t[i : i + 20] for i in range(len(t) - 19))
    barcode_hits = sorted(barcodes & seen)
    assert not barcode_hits, f"provisional barcodes inside inserts: {barcode_hits[:5]}"
    summary = {
        "input_refs": RUNG_REFS,
        "input_sha256": shas,
        "inserts": len(rows),
        "rules": {
            "gc50_max": GC50_MAX,
            "cassette_end_gc50_margin": CASSETTE_END_GC50,
            "gc100_refused_above": GC100_REFUSED,
            "me_edit_distance_refused_at_or_below": ME_MAX_REFUSED,
            "repeat_target_below_nt": REPEAT_TARGET,
            "near_repeat_target_below_nt": NEAR_TARGET,
            "repeat_refused_at_nt": REPEAT_REFUSED,
            "min_codon_usage": MIN_USAGE,
        },
        "status": dict(Counter(r["status"] for r in rows)),
        "category_before": dict(Counter(r["before_category"] for r in rows)),
        "category_after": dict(Counter(r["after_category"] for r in rows)),
        "cassette_end_windows_over_margin": {
            "before": sum(r["before_cassette_end_gc50"] > CASSETTE_END_GC50 for r in rows),
            "after": sum(r["after_cassette_end_gc50"] > CASSETTE_END_GC50 for r in rows),
        },
        "codons_recoded": len(recodes),
        "inserts_recoded": len({r["arm_id"] for r in recodes}),
        "members_recoded": len({r["design_id"] for r in recodes}),
        "members": len({r["design_id"] for r in rows}),
        "recodes_by_reason": dict(Counter(r["reason"] for r in recodes)),
        "provisional_barcodes_checked": len(barcodes),
        "provisional_barcodes_inside_an_insert": len(barcode_hits),
        "recode_log_sha256": hashlib.sha256(tsv_text(recodes).encode()).hexdigest(),
        "recoded_fasta": f"processed://{FASTA_KEY}",
        "recoded_fasta_sha256": hashlib.sha256(fasta_text.encode()).hexdigest(),
    }
    outputs = {
        OUT: tsv_text(rows),
        RECODES: tsv_text(recodes),
        SUMMARY: json.dumps(summary, indent=2) + "\n",
    }
    if check:
        stale = [p.name for p, text in outputs.items() if p.read_text() != text]
        if stale:
            raise SystemExit(f"stale: {stale}; rerun synthesis_compliance.py")
        print(f"{len(rows)} inserts; {summary['codons_recoded']} codons recoded; tables reproduce")
        return 0
    for p, text in outputs.items():
        p.write_text(text)
    (HERE / "_synthesis_recoded.fasta").write_text(fasta_text)
    print(json.dumps(summary, indent=2))
    return 0


def calibrate() -> int:
    """Local rules against Twist's calls on the rebuilt draft inserts."""
    import twist_complexity_map as tcm

    from deliverome_analysis import library_design as ld
    from deliverome_analysis.s3 import processed

    twist = tcm.load_twist()
    man = {}
    with (HERE / "candidate_manifest_tranches.tsv").open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            man[r["design_id"]] = r
    recs = [
        json.loads(Path(processed.download(k)).read_text())
        for k in (ld.BASE_CONSTRUCTS_KEY, ld.CONSTRUCTS_KEY)
    ]
    agree = Counter()
    for t in twist:
        key = man[t["design_id"]]["key"]
        s = next(
            ins
            for rec in recs
            if key in rec
            for ins in [
                qc.ESP3I_5 + qc.KOZAK + rec[key]["construct"] + qc.REV_PRIMER + qc.CS1 + qc.ESP3I_3
            ]
            if len(ins) == int(t["insert_nt"])
        )
        agree[(t["Score"], category(metrics(s)))] += 1
    for (tw, mine), n in sorted(agree.items()):
        print(f"  Twist {tw:<13} local {mine:<13} {n}")
    return 0


def publish() -> int:
    from deliverome_analysis.s3 import processed

    local = HERE / "_synthesis_recoded.fasta"
    summ = json.loads(SUMMARY.read_text())
    assert hashlib.sha256(local.read_bytes()).hexdigest() == summ["recoded_fasta_sha256"]
    if processed.exists(FASTA_KEY):
        raise SystemExit(f"processed://{FASTA_KEY} exists; bump FASTA_KEY's date for a new run")
    processed.upload(str(local), FASTA_KEY)
    local.unlink()
    print(f"published processed://{FASTA_KEY}")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--calibrate", action="store_true")
    ap.add_argument("--publish", action="store_true")
    ap.add_argument("--check", action="store_true", help="committed tables are current")
    args = ap.parse_args(argv)
    if args.calibrate:
        return calibrate()
    if args.publish:
        return publish()
    return run(check=args.check)


if __name__ == "__main__":
    raise SystemExit(main())
