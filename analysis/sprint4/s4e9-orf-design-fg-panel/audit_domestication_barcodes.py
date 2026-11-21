"""Audit two design steps across the whole L1 library: domestication and barcodes.

The question for each step is whether what the pipeline claims is true of the final inserts,
and which claims nothing tests. Every number below is recomputed from the inserts themselves.
None is copied from a summary. The final inserts are the five tag-rung arm tables on the draft
L1 integration branch (PR #299), read at a pinned commit straight from git. The multipass rung is
read after the one-arm-per-gene rebuild (PR #312), which retired the 25 expression-rescue members:

| Rung | Arm table |
|---|---|
| C-terminal and snorkel | ``constructs_l1_cterm.fasta`` |
| N after signal peptide | ``constructs_l1_nterm.fasta`` |
| N at residue 1, multipass | ``constructs_l1_multipass.fasta`` |
| N at residue 1, type III | ``constructs_l1_type3.fasta`` |
| internal | ``constructs_l1_internal.fasta`` |

**Domestication**, per insert:

* Type IIS: exactly the two wrapper Esp3I sites (``CGTCTC`` at the 5' end, ``GAGACG`` at the 3'
  end), and no BbsI, BsaI or other Esp3I site on either strand. The audit also screens enzymes
  the s4e9 code never lists (AarI, PaqCI, BtgZI, SapI). They are reported, not gated.
* Kozak: the insert opens with the Esp3I site, the p003 Kozak block and ATG.
* Protein: the coding sequence (Kozak to the first in-frame stop) translates to the protein its
  arm table records (``protein_sha256``, or ``mature_sha256`` where the table has one), with
  exactly one stop, at the end.
* Edits: every difference between a member's designed CDS (``CDS_KEY``) and its domesticated
  ORF (``CONSTRUCTS_KEY``) is logged per codon in ``l1_domestication_edits.tsv``, with its
  position, old codon and new codon, and is checked to be synonymous. The record's ``edits``
  count is checked against the logged nucleotide differences.
* Unlogged edits: with its cassette cut out, an insert's coding sequence should equal the
  domesticated ORF. Where an Igκ leader replaces the first residues (``igk_replaces_residues``), it
  should equal the leader followed by the domesticated ORF from the first kept residue. Native
  nucleotides that are not verbatim are counted per insert (``native_nt_not_verbatim``; -1 means a
  leader arm breaks that rule).
* Motifs: domestication only claims to clear the three enzymes. Its edits are still screened for a
  new ATG (any frame), splice donor (``GT[AG]AG``) or poly(A) signal (``AATAAA``, ``ATTAAA``).
  These counts are informational, because nothing in the pipeline claims them.

**Barcodes.** No member has a barcode yet: every insert carries a 20-N slot. What exists is the
block of 3,200 pool barcodes reserved for L1 (``l1_barcode_block_provisional.tsv``,
``reserve_l1_block.py``). The audit re-derives that block's properties:

* length, alphabet, GC and longest homopolymer;
* pairwise Hamming distance over all 20 nt and over the 10-nt in situ prefix;
* every screen the brief states, re-run (Type IIS, splice donor and junction miRNA seeds under
  both 3' flanks, and the RBP motifs);
* no barcode, on either strand, occurs inside any insert or the p0039 backbone;
* uniqueness against the shared registry;
* that the draw is deterministic.

**Swap tract.** The longest stretch each insert shares with any other member's insert, measured
on the final DNA up to the barcode slot. Arms of the same member (alternative designs) are not
compared with each other.

**Twist recodes.** The order of record is the Twist-compliant recoded set
(``synthesis_compliance.py``, PR #303), pinned here by SHA-256. Every invariant is re-run on it, and
each recoded insert is compared with its original:

* protein, table SHA-256, wrapper, Kozak, readout block and barcode slot unchanged;
* every codon change synonymous, and the count equal to ``synthesis_recodes.tsv``;
* cassette DNA unchanged, and on the C-terminal rung still the member's issued encoding with a
  unique terminal linker block;
* no new Type IIS or consensus splice-donor site;
* the swap tract recomputed and compared with the original's.

The pinned order predates that rebuild. Its retired rescue inserts are counted
(``retired_members_still_in_order``), and inserts rebuilt since it was made are listed
(``rebuilt_since_order_not_compared``) instead of being compared.

``--recoded FASTA`` checks a different file instead.

Writes ``l1_insert_audit.tsv`` (one row per insert), ``l1_domestication_edits.tsv`` (one row per
edited codon) and ``l1_domestication_barcode_audit.json``.

Implemented by the Claude Code session l1-cterm-build (branch
claude/l1-domestication-barcode-audit) on 2026-11-12, at the coordinator's request for Becca's
overnight review of every in-silico design step. Report only: it changes no design.

    uv run --locked python audit_domestication_barcodes.py
    uv run --locked python audit_domestication_barcodes.py --recoded twist_recoded.fasta
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import re
import subprocess
import sys
from collections import Counter
from pathlib import Path

import l1_arm_ids as arm_ids
import numpy as np
from build_constructs import CASSETTE_AA, CS1, FWD_PRIMER, KOZAK, REV_PRIMER, SPACER

from deliverome_analysis import library_design as ld
from deliverome_analysis import seqtools

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
REL = HERE.relative_to(REPO)
sys.path.insert(0, str(REPO / "scratch/becca/appris-principal"))

import assemble_candidates  # noqa: E402

#: The commit whose arm tables are audited: main at cdb78068 (2026-11-14), the build the
#: 2026-11-14-final order of record was made from (synthesis_compliance.INPUT_REF). Earlier
#: audits read the draft integration branch (496d3482, then 5cfd91c0, PR #299), not on main.
INTEGRATION_REF = "8432c0104c8fdaa57e78173d18fc5bb7064d91bd"
RUNGS = {
    "C-terminal/snorkel": "constructs_l1_cterm",
    "N after SP": "constructs_l1_nterm",
    "N at residue 1, multipass": "constructs_l1_multipass",
    "N at residue 1, type III": "constructs_l1_type3",
    "internal": "constructs_l1_internal",
}
#: The 25 expression-rescue members (l1-p2361..l1-p2385) are retired, one arm per gene (PR #312,
#: Becca 2026-11-13: "make a single choice for each L1"); ``constructs_l1_rescue`` is kept as a
#: record and is not ordered, so the audit does not read it.
#: A rung read at a later commit than INTEGRATION_REF, compared with the order as rebuilt since
#: it (``rebuilt_since_order``). Empty: every rung is read at the commit the order was built from.
RUNG_REF: dict[str, str] = {}
#: The arm tables' short rung names (l1_arm_ids, synthesis_compliance).
SHORT_RUNG = {
    "C-terminal/snorkel": "cterm",
    "N after SP": "nterm",
    "N at residue 1, multipass": "multipass",
    "N at residue 1, type III": "type3",
    "internal": "internal",
}
#: Design ids of the retired rescue members; the superseded 2026-11-13 order still held them.
RETIRED_RESCUE = {f"l1-p{n:04d}" for n in range(2361, 2386)}
GATED_ENZYMES = ["BbsI", "Esp3I", "BsaI"]
EXTRA_TYPE_IIS = {  # recognition sites, top strand; reported only
    "AarI": "CACCTGC",
    "PaqCI": "CACCTGC",
    "BtgZI": "GCGATG",
    "SapI": "GCTCTTC",
}
WRAP_5 = assemble_candidates.ESP3I_5
WRAP_3 = assemble_candidates.ESP3I_3
BARCODE_SLOT = "N" * 20
SPLICE_DONOR = re.compile(r"GT[AG]AG")
POLYA = ("AATAAA", "ATTAAA")

OUT_INSERTS = HERE / "l1_insert_audit.tsv"
OUT_EDITS = HERE / "l1_domestication_edits.tsv"
OUT_JSON = HERE / "l1_domestication_barcode_audit.json"


# --------------------------------------------------------------------------- inputs


def git_show(path: str, ref: str = INTEGRATION_REF) -> str:
    return subprocess.run(
        ["git", "show", f"{ref}:{REL}/{path}"],
        cwd=REPO,
        check=True,
        capture_output=True,
        text=True,
    ).stdout


def read_fasta(text: str) -> list[tuple[str, str]]:
    out, name, seq = [], None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if name is not None:
                out.append((name, "".join(seq)))
            name, seq = line[1:], []
        else:
            seq.append(line.strip())
    if name is not None:
        out.append((name, "".join(seq)))
    return out


def arm_rows(stem: str, ref: str) -> dict[str, dict]:
    rows = list(csv.DictReader(io.StringIO(git_show(f"{stem}.tsv", ref)), delimiter="\t"))
    key = "construct_id" if rows and "construct_id" in rows[0] else "candidate_id"
    out = {r[key]: r for r in rows}
    # A snorkel row's FASTA header leads with its arm id (constructs_l1_cterm.arm_id).
    out.update({r["arm_id"]: r for r in rows if r.get("arm_id")})
    return out


def insert_key(ident: str, design_id: str) -> str:
    """An insert's identity that survives a member renumber: member ids (``l1-m…``) renumber,
    so an insert named by one is keyed by its permanent design id; arm ids built on the design
    id (``l1-pNNNN:…``) already are."""
    return design_id if ident.startswith("l1-m") and design_id else ident


def load_inserts(ref: str) -> list[dict]:
    """Every built insert, keyed by its permanent arm id (``l1_arm_ids``); ``legacy_key`` is
    the key an order of record cut before the per-insert numbering used, so its recoded file
    can still be matched."""
    decisions = arm_ids.load_decisions()
    out = []
    for rung, stem in RUNGS.items():
        rung_ref = RUNG_REF.get(rung, ref)
        rows = arm_rows(stem, rung_ref)
        for header, seq in read_fasta(git_show(f"{stem}.fasta", rung_ref)):
            ident = re.split(r"[|\s]", header)[0]
            row = rows.get(ident)
            assert row is not None, f"{stem}: no table row for {ident}"
            design = row.get("design_id", "")
            out.append(
                {
                    "rung": rung,
                    "id": ident,
                    "arm_id": arm_ids.arm_id_for(SHORT_RUNG[rung], row, design, decisions),
                    "legacy_key": insert_key(ident, design),
                    "header": header,
                    "insert": seq.upper(),
                    "row": row,
                }
            )
    for i in out:
        i["key"] = i["arm_id"]
    arm_ids.check_unique(i["arm_id"] for i in out)
    return out


# --------------------------------------------------------------------------- inserts


def type_iis(insert: str) -> dict:
    sites = seqtools.find_sites(insert, GATED_ENZYMES)
    esp = [s for s in sites if s[0] == "Esp3I"]
    wrapper = (
        insert.startswith(WRAP_5)
        and insert.endswith(WRAP_3)
        and len(esp) == 2
        and len(sites) == 2
        and {s[2] for s in esp} == {0, len(insert) - 6}
    )
    rc = seqtools.reverse_complement(insert)
    extra = sorted({name for name, site in EXTRA_TYPE_IIS.items() if site in insert or site in rc})
    return {"sites": sites, "wrapper_only": wrapper, "extra_type_iis": extra}


def coding_of(insert: str) -> tuple[str, str]:
    """(coding sequence from ATG to the first in-frame stop inclusive, a problem or '')."""
    head = WRAP_5 + KOZAK
    if not insert.startswith(head + "ATG"):
        return "", "Kozak block or start codon not as designed"
    body = insert[len(head) :]
    for i in range(0, len(body) - 2, 3):
        if seqtools.translate(body[i : i + 3]) == "*":
            return body[: i + 3], ""
    return body, "no in-frame stop"


def protein_matches(protein: str, row: dict) -> bool:
    import hashlib

    h = hashlib.sha256(protein.encode()).hexdigest()
    return h in {row.get("protein_sha256"), row.get("mature_sha256")}


def verbatim_native(coding: str, orf: str) -> int:
    """Native nucleotides of the domesticated ORF (stop excluded) not present verbatim in the
    insert's coding sequence as one shared prefix plus one shared suffix. The cassette (or a
    leader) sits between them, so anything else is a change the arm made to native codons."""
    a, b = coding[:-3], orf[:-3]
    p = 0
    while p < min(len(a), len(b)) and a[p] == b[p]:
        p += 1
    s = 0
    while s < min(len(a), len(b)) - p and a[-1 - s] == b[-1 - s]:
        s += 1
    return max(0, len(b) - p - s)


# --------------------------------------------------------------------------- domestication edits


def edit_log(manifest: dict, cds: dict, con: dict) -> tuple[list[dict], dict]:
    rows, stats = [], Counter()
    for cid, m in manifest.items():
        if m["key"].startswith("rescue:") or m["design_id"] in RETIRED_RESCUE:
            stats["retired_rescue_members_skipped"] += 1
            continue
        stats["members"] += 1
        if m["key"] not in cds or m["key"] not in con:
            stats["not_in_cds_record"] += 1
            continue
        else:
            raw = cds[m["key"]]["cds"]
        dom = con[m["key"]]["construct"]
        if len(raw) != len(dom):
            stats["length_changed"] += 1
            continue
        nt_diff = sum(1 for x, y in zip(raw, dom, strict=True) if x != y)
        recorded = con[m["key"]]["edits"]
        if recorded not in (None, "None") and nt_diff != int(recorded or 0):
            stats["edit_count_mismatch"] += 1
        if seqtools.translate(raw) != seqtools.translate(dom):
            stats["translation_changed"] += 1
        stats["members_edited"] += nt_diff > 0
        for i in range(0, len(raw), 3):
            old, new = raw[i : i + 3], dom[i : i + 3]
            if old == new:
                continue
            win_old, win_new = raw[max(0, i - 5) : i + 8], dom[max(0, i - 5) : i + 8]
            row = {
                "design_id": m["design_id"],
                "candidate_id": cid,
                "gene_symbol": ld.gene_from_key(m["key"]),
                "codon": i // 3 + 1,
                "nt_position": i + 1,
                "old_codon": old,
                "new_codon": new,
                "aa": seqtools.translate(old),
                "synonymous": seqtools.translate(old) == seqtools.translate(new),
                "new_atg": win_new.count("ATG") > win_old.count("ATG"),
                "new_splice_donor": len(SPLICE_DONOR.findall(win_new))
                > len(SPLICE_DONOR.findall(win_old)),
                "new_polya": sum(win_new.count(p) for p in POLYA)
                > sum(win_old.count(p) for p in POLYA),
            }
            rows.append(row)
            stats["codons_edited"] += 1
            stats["nonsynonymous"] += not row["synonymous"]
            stats["edits_new_atg"] += row["new_atg"]
            stats["edits_new_splice_donor"] += row["new_splice_donor"]
            stats["edits_new_polya"] += row["new_polya"]
    return rows, dict(stats)


# --------------------------------------------------------------------------- swap tract


def upstream_of_barcode(insert: str) -> str | None:
    i = insert.find(BARCODE_SLOT)
    if i < 0:
        return None
    return insert[len(WRAP_5) : i]


def same_member_swap(rows: list[dict], ups: list[str | None], have: list[int]) -> dict:
    """The shared tract up to the barcode between two inserts of one member, and the swap rate
    it predicts (``seqtools.barcode_swap_rate``, Walton). Different-member pairs are measured by
    ``worst_partner_excluding_self``; two arms of one member share the native sequence after
    the later site, so their pairs are reported here, apart (Becca, 2026-11-15: one pool; a
    swap between them misattributes the site, not the protein)."""
    by_member: dict[str, list[int]] = {}
    for k in have:
        by_member.setdefault(rows[k]["design_id"], []).append(k)
    pairs = []
    for ks in by_member.values():
        for a in range(len(ks)):
            for b in range(a + 1, len(ks)):
                x, y = ups[ks[a]], ups[ks[b]]
                n = 0
                while n < min(len(x), len(y)) and x[-1 - n] == y[-1 - n]:
                    n += 1
                pairs.append(n)
    multi = [ks for ks in by_member.values() if len(ks) > 1]
    for r in rows:
        r["same_member_inserts"] = 1
    for ks in multi:
        for k in ks:
            rows[k]["same_member_inserts"] = len(ks)
    gate = 100
    if not pairs:
        return {"pairs": 0, "members": 0, "inserts_in_multi_arm_members": 0}
    arr = np.array(pairs)
    rate = np.array([seqtools.barcode_swap_rate(int(n)) for n in arr])
    return {
        "pairs": len(pairs),
        "members": len(multi),
        "inserts_in_multi_arm_members": sum(len(ks) for ks in multi),
        "fraction_of_inserts": round(sum(len(ks) for ks in multi) / len(have), 4),
        "tract_nt_median": float(np.median(arr)),
        "tract_nt_min": int(arr.min()),
        "tract_nt_max": int(arr.max()),
        "pairs_over_100": int((arr > gate).sum()),
        "modelled_swap_median": round(float(np.median(rate)), 4),
        "modelled_swap_max": round(float(rate.max()), 4),
        "modelled_swap_at_100_nt_gate": round(seqtools.barcode_swap_rate(gate), 4),
    }


def worst_partner_excluding_self(upstream: list[str], owner: list[str]) -> list[tuple[int, int]]:
    """For each sequence, (longest common suffix with a different owner's sequence, partner)."""

    def lcs(a: str, b: str) -> int:
        n = 0
        while n < min(len(a), len(b)) and a[-1 - n] == b[-1 - n]:
            n += 1
        return n

    order = sorted(range(len(upstream)), key=lambda i: upstream[i][::-1])
    adj = [lcs(upstream[order[k]], upstream[order[k + 1]]) for k in range(len(order) - 1)]
    out = [(0, -1)] * len(upstream)
    for pos, i in enumerate(order):
        best = (0, -1)
        for step in (-1, 1):
            j, run = pos, 10**9
            while 0 <= j + step < len(order):
                run = min(run, adj[min(j, j + step)])
                j += step
                if run <= best[0]:
                    break
                if owner[order[j]] != owner[i]:
                    best = (run, order[j])
                    break
        out[i] = best
    return out


# --------------------------------------------------------------------------- barcodes


def hamming_min(codes: list[str], width: int) -> tuple[int, int]:
    """(minimum pairwise Hamming distance over the first ``width`` nt, pairs at that minimum)."""
    lut = {c: i for i, c in enumerate("ACGT")}
    arr = np.array([[lut[c] for c in s[:width]] for s in codes], dtype=np.int8)
    best, count = width + 1, 0
    for i in range(len(arr) - 1):
        d = (arr[i + 1 :] != arr[i]).sum(axis=1)
        m = int(d.min())
        if m < best:
            best, count = m, int((d == m).sum())
        elif m == best:
            count += int((d == m).sum())
    return best, count


def levenshtein_min(codes: list[str], below: int = 3) -> tuple[int, list[tuple[str, str]]]:
    """Pairs closer than ``below`` in edit distance, found exactly.

    Two 20-mers within edit distance 2 share an identical 6-mer in each of three blocks' worth
    of positions (pigeonhole over blocks 0-5, 6-11, 12-19, shifted by at most 2), so every such
    pair shares some 6-mer. Candidates come from a 6-mer index, and each is scored with a full
    dynamic program, so the result is exact, not a bound.
    """
    from collections import defaultdict

    k = 6
    index = defaultdict(set)
    for i, c in enumerate(codes):
        for j in range(len(c) - k + 1):
            index[c[j : j + k]].add(i)

    def lev(a: str, b: str) -> int:
        prev = list(range(len(b) + 1))
        for i, x in enumerate(a, 1):
            cur = [i]
            for j, y in enumerate(b, 1):
                cur.append(min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (x != y)))
            prev = cur
        return prev[-1]

    seen, close, best = set(), [], 99
    for ids in index.values():
        for i in ids:
            for j in ids:
                if i < j and (i, j) not in seen:
                    seen.add((i, j))
                    d = lev(codes[i], codes[j])
                    best = min(best, d)
                    if d < below:
                        close.append((codes[i], codes[j]))
    return best, close


def longest_homopolymer(s: str) -> int:
    return max(len(m.group(0)) for m in re.finditer(r"(.)\1*", s))


def barcode_audit(inserts: list[dict], backbone: str) -> dict:
    import reserve_barcodes as rb
    import reserve_l1_block as rlb
    from barcode_prefix_screen import PREFIX_BASES, PREFIX_MIN_DISTANCE

    with (HERE / "l1_barcode_block_provisional.tsv").open() as fh:
        block = list(csv.DictReader(fh, delimiter="\t"))
    # The 2026-11-15 extension (L1-res-3201 on, Becca) is audited with the block, as one set.
    extension = []
    if rlb.EXT_OUT.exists():
        with rlb.EXT_OUT.open() as fh:
            extension = list(csv.DictReader(fh, delimiter="\t"))
    codes = [r["barcode"].strip().upper() for r in block + extension]
    gc = [seqtools.gc_content(c) for c in codes]
    hp = [longest_homopolymer(c) for c in codes]
    screen_fail = Counter()
    for c in codes:
        ok, why = rb.l1_screens_pass(c)
        if not ok:
            screen_fail[why.split(" [")[0]] += 1
        ok, why = rlb.rbp_pass(c)
        if not ok:
            screen_fail["RBP motif"] += 1
    rc = [seqtools.reverse_complement(c) for c in codes]
    hay = [backbone, seqtools.reverse_complement(backbone)] + [i["insert"] for i in inserts]
    hay += [seqtools.reverse_complement(i["insert"]) for i in inserts]
    joined = "|".join(hay)
    in_dna = [c for c, r in zip(codes, rc, strict=True) if c in joined or r in joined]

    from deliverome_analysis.s3 import processed

    text = Path(processed.download(rb.REGISTRY_KEY)).read_text()
    body = [ln for ln in text.splitlines() if ln and not ln.startswith("#")]
    reg = list(csv.DictReader(body, delimiter="\t"))
    reg_count = Counter(r["barcode"].strip().upper() for r in reg if r["barcode"])
    block_rows = [r for r in reg if r["construct"].startswith("L1-res-")]
    elsewhere = [c for c in codes if reg_count[c] != 1]

    # Determinism: the draw takes pool barcodes in index order, no random number generator.
    # Replaying it against the registry as it stood before the block (rows before the first
    # L1-res row) must give back the same 3,200, in the same order.
    first = next(i for i, r in enumerate(reg) if r["construct"].startswith("L1-res-"))
    before = reg[:first]
    used_bc = {r["barcode"].strip().upper() for r in before if r["barcode"]}
    used_ix = {int(r["pool_index"]) for r in before if r["pool_index"]}
    real_read = rlb.read_registry

    def replay_registry():
        _bc, _ix, skipped, done, used_p = real_read()
        return used_bc, used_ix, skipped, done, used_p

    rlb.read_registry = replay_registry
    try:
        replay, _why = rlb.draw(len(block))
    finally:
        rlb.read_registry = real_read
    replay_codes = [bc for _ix, bc in replay]
    if extension:
        # The extension replays against the registry as it stood before it (rows before
        # L1-res-3201): every passing barcode after the block's last pool index.
        first_ext = next(
            i for i, r in enumerate(reg) if r["construct"] == extension[0]["construct"]
        )
        before = reg[:first_ext]
        used_bc = {r["barcode"].strip().upper() for r in before if r["barcode"]}
        used_ix = {int(r["pool_index"]) for r in before if r["pool_index"]}
        rlb.read_registry = replay_registry
        try:
            last = max(int(r["pool_index"]) for r in block)
            ext_replay, _why = rlb.draw(None, after=last)
        finally:
            rlb.read_registry = real_read
        replay_codes += [bc for _ix, bc in ext_replay]

    full_min, full_pairs = hamming_min(codes, 20)
    lev_min, lev_close = levenshtein_min(codes)
    pre_min, pre_pairs = hamming_min(codes, PREFIX_BASES)
    return {
        "block_size": len(codes),
        "block_rows": len(block),
        "extension_rows": len(extension),
        "unique": len(set(codes)) == len(codes),
        "length_20_acgt": all(len(c) == 20 and not set(c) - set("ACGT") for c in codes),
        "gc_min": round(min(gc), 3),
        "gc_max": round(max(gc), 3),
        "longest_homopolymer_max": max(hp),
        "contains_atg": sum("ATG" in c for c in codes),
        "hamming_min_full_20nt": full_min,
        "pairs_at_full_min": full_pairs,
        "levenshtein_min_over_kmer_candidates": lev_min,
        "pairs_levenshtein_below_3": lev_close,
        "pool_claim": "s1e2 pool: min Levenshtein 3 and min Hamming 4 (s1e2 README)",
        "hamming_min_prefix": pre_min,
        "prefix_rule": f"{PREFIX_MIN_DISTANCE} over the first {PREFIX_BASES} nt",
        "pairs_at_prefix_min": pre_pairs,
        "screen_failures_both_flanks": dict(screen_fail),
        "found_in_an_insert_or_backbone": in_dna,
        "registry_block_rows": len(block_rows),
        "not_exactly_once_in_registry": elsewhere,
        "draw_replays_identically": replay_codes == codes,
        "draw_replay_note": "replayed against the registry rows that precede the block, with "
        "today's local barcode tables",
    }


# --------------------------------------------------------------------------- main


def audit(inserts: list[dict], cds: dict, con: dict, manifest_by_design: dict) -> list[dict]:
    rows = []
    for ins in inserts:
        r = ins["row"]
        t = type_iis(ins["insert"])
        coding, problem = coding_of(ins["insert"])
        protein = seqtools.translate(coding)
        m = manifest_by_design.get(r.get("design_id", ""))
        native_off = ""
        if m is not None and coding and m["key"] in con:
            # Cut the cassette out, then compare with the domesticated ORF. Residues a leader
            # replaces (``igk_replaces_residues``) are expected to differ; nothing else is.
            k = protein.find(CASSETTE_AA)
            core = coding if k < 0 else coding[: 3 * k] + coding[3 * (k + len(CASSETTE_AA)) :]
            replaced = r.get("igk_replaces_residues") or ""
            n_rep = 0
            if replaced:
                lo, hi = (int(x) for x in replaced.split("-"))
                n_rep = hi - lo + 1
            orf = con[m["key"]]["construct"]
            if n_rep:
                native_off = 0 if core.endswith(orf[3 * n_rep :]) else -1
            else:
                native_off = verbatim_native(core, orf)
        rows.append(
            {
                "rung": ins["rung"],
                "id": ins["id"],
                "arm_id": ins["arm_id"],
                "design_id": r.get("design_id", ""),
                "gene_symbol": r.get("gene_symbol", ""),
                "insert_nt": len(ins["insert"]),
                "type_iis_sites": len(t["sites"]),
                "wrapper_only": t["wrapper_only"],
                "extra_type_iis_reported": ";".join(t["extra_type_iis"]),
                "kozak_and_start": not problem.startswith("Kozak"),
                "single_stop_at_end": bool(coding)
                and protein.count("*") == 1
                and protein.endswith("*"),
                "protein_matches_table": bool(coding) and protein_matches(protein.rstrip("*"), r),
                "barcode_slots": ins["insert"].count(BARCODE_SLOT),
                "native_nt_not_verbatim": native_off,  # -1: leader arm not leader + ORF
                "problem": problem,
            }
        )
    return rows


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ref", default=INTEGRATION_REF, help="git commit holding the arm tables")
    ap.add_argument(
        "--recoded",
        type=Path,
        help="a recoded FASTA to check instead of the pinned order of record",
    )
    args = ap.parse_args(argv)

    import candidate_inputs
    from Bio import SeqIO

    from deliverome_analysis.s3 import processed

    inserts = load_inserts(args.ref)
    cds = json.loads(Path(processed.download(ld.CDS_KEY)).read_text())
    con = json.loads(Path(processed.download(ld.CONSTRUCTS_KEY)).read_text())
    manifest = {
        r["candidate_id"]: {**r, "design_id": ld.design_id(r)}
        for r in csv.DictReader(
            io.StringIO(git_show("candidate_manifest_tranches.tsv", args.ref)), delimiter="\t"
        )
        if r["library"] == ld.LIBRARY_L1 and ld.is_active(r)
    }
    by_design = {m["design_id"]: m for m in manifest.values()}
    backbone = str(SeqIO.read(candidate_inputs.backbone(), "genbank").seq).upper()

    rows = audit(inserts, cds, con, by_design)
    edits, edit_stats = edit_log(manifest, cds, con)

    ups = [upstream_of_barcode(i["insert"]) for i in inserts]
    have = [k for k, u in enumerate(ups) if u is not None]
    owner = [rows[k]["design_id"] or rows[k]["id"] for k in have]
    worst = worst_partner_excluding_self([ups[k] for k in have], owner)
    for (tract, partner), k in zip(worst, have, strict=True):
        rows[k]["worst_partner_tract_nt"] = tract
        rows[k]["worst_partner"] = rows[have[partner]]["id"] if partner >= 0 else ""
    tracts = np.array([w[0] for w in worst])
    same_member = same_member_swap(rows, ups, have)

    summary = {
        "integration_ref": args.ref,
        "inserts": len(rows),
        "by_rung": dict(Counter(r["rung"] for r in rows)),
        "insert_checks_failed": {
            k: [f"{r['rung']}:{r['id']}" for r in rows if not r[k]]
            for k in (
                "wrapper_only",
                "kozak_and_start",
                "single_stop_at_end",
                "protein_matches_table",
            )
        },
        "extra_type_iis_reported": dict(
            Counter(e for r in rows for e in r["extra_type_iis_reported"].split(";") if e)
        ),
        "barcode_slots_not_one": [
            f"{r['rung']}:{r['id']}" for r in rows if r["barcode_slots"] != 1
        ],
        "native_nt_not_verbatim": {
            "inserts_with_any": sum(1 for r in rows if r["native_nt_not_verbatim"] not in ("", 0)),
            "by_rung": dict(
                Counter(r["rung"] for r in rows if r["native_nt_not_verbatim"] not in ("", 0))
            ),
        },
        "domestication": edit_stats,
        "swap_tract_final_dna": {
            "inserts_measured": len(tracts),
            "median_nt": float(np.median(tracts)),
            "max_nt": int(tracts.max()),
            "over_80": int((tracts > 80).sum()),
            "over_100": int((tracts > 100).sum()),
        },
        "swap_tract_same_member": same_member,
        "barcodes": barcode_audit(inserts, backbone),
    }
    recoded, digest = recoded_inserts(args.recoded)
    # An order cut before the per-insert numbering names inserts the old way; key it by arm.
    arm_of_legacy = {(i["rung"], i["legacy_key"]): i["arm_id"] for i in inserts}
    for r in recoded:
        if "key" not in r:
            r["key"] = arm_of_legacy.get((r["rung"], r["legacy_key"]), r["legacy_key"])
    if args.recoded is None and digest != RECODED_SHA256:
        raise SystemExit(f"{RECODED_KEY} changed (sha256 {digest}); re-pin it deliberately")
    # The encodings as of the audited commit, so their candidate ids match its arm tables
    # (a later renumber re-keys the working-tree copy).
    by_design_cand = {
        r["candidate_id"]: r["cassette_nt"]
        for r in csv.DictReader(
            io.StringIO(git_show("l1_cassette_encodings.tsv", args.ref)), delimiter="\t"
        )
    }
    stale = rebuilt_since_order(inserts)
    rec_summary, facts = recoded_check(recoded, inserts, by_design_cand, stale)
    # The swap tract again, on the recoded DNA, per insert against its original. Retired members
    # are left out; a rebuilt insert has no recoded version yet, so its built DNA stands in.
    rec_map = {
        (r["rung"], r["key"]): r["insert"] for r in recoded if r["rung"] != "expression rescue"
    }
    rec_map.update(
        {(i["rung"], i["key"]): i["insert"] for i in inserts if (i["rung"], i["key"]) in stale}
    )
    rups = [
        upstream_of_barcode(rec_map.get((inserts[k]["rung"], inserts[k]["key"]), "")) for k in have
    ]
    rworst = worst_partner_excluding_self([u or "" for u in rups], owner)
    rt = np.array([w[0] for w in rworst])
    for (tract, _p), k in zip(rworst, have, strict=True):
        rows[k]["recoded_codons"] = facts.get((inserts[k]["rung"], inserts[k]["key"]), {}).get(
            "recoded_codons", ""
        )
        rows[k]["recoded_worst_partner_tract_nt"] = tract
    rec_summary["swap_tract_recoded_dna"] = {
        "median_nt": float(np.median(rt)),
        "max_nt": int(rt.max()),
        "over_100": int((rt > 100).sum()),
        "inserts_whose_tract_grew": int((rt > tracts).sum()),
        "largest_growth_nt": int((rt - tracts).max()),
    }
    rec_summary["recoded_fasta_sha256"] = digest
    summary["recoded"] = rec_summary

    for path, table in ((OUT_INSERTS, rows), (OUT_EDITS, edits)):
        with path.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(table[0]), delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(table)
    OUT_JSON.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({k: v for k, v in summary.items() if k != "barcodes"}, indent=2)[:6000])
    print(json.dumps(summary["barcodes"], indent=2))
    return 0


RECODED_KEY = "sprint4/s4e9-orf-design-fg-panel/synthesis/2026-11-21-snorkel-modules/l1_order_inserts_recoded.fasta"
#: The Twist-compliant recoded inserts, the order of record (synthesis_compliance.py), regenerated
#: by l1-reviewer-package from main at cdb78068 on 2026-11-14: 2,629 inserts on the five rungs,
#: named by member id (ordered-length numbering). Superseded pins: 2026-11-13-one-arm-renamed/
#: b426e433 (2,678 inserts, draft integration 496d3482), 2026-11-13/ d1257f0a (2,716 inserts,
#: with the 25 retired rescue arms) and the undated #303 slice 9b1c9f55.
RECODED_SHA256 = "66332be47f53968fadf1eb3953bae584d49393028d6d1c9c64eb0baa6e3623d5"
#: The synthesis_recodes.tsv this order was made with (l1-reviewer-package, final regeneration).
RECODES_REF = "3b00f73ea842a69150738bdabd9c7cdbd79b9436"
RUNG_OF_TAG = {
    "cterm": "C-terminal/snorkel",
    "nterm": "N after SP",
    "multipass": "N at residue 1, multipass",
    "type3": "N at residue 1, type III",
    "internal": "internal",
    "rescue": "expression rescue",
}


def recoded_inserts(path: Path | None) -> tuple[list[dict], str]:
    """The recoded set, keyed like the originals (rung, first header token)."""
    import hashlib

    if path is None:
        from deliverome_analysis.s3 import processed

        path = Path(processed.download(RECODED_KEY))
    raw = path.read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    out = []
    for header, seq in read_fasta(raw.decode()):
        tok = header.split()
        rung = RUNG_OF_TAG[tok[-1]]
        if ld.ARM_ID_RE.fullmatch(tok[0]):  # per-insert order: arm id, insert id, design, rung
            out.append({"rung": rung, "id": tok[0], "key": tok[0], "insert": seq.upper()})
            continue
        ident, design_id = re.split(r"[|\s]", header)[0], tok[1]
        out.append(
            {
                "rung": rung,
                "id": ident,
                "legacy_key": insert_key(ident, design_id),
                "insert": seq.upper(),
            }
        )
    return out, digest


def rebuilt_since_order(inserts: list[dict]) -> set[tuple[str, str]]:
    """Inserts whose DNA changed after the order of record was built: the rungs read at a later
    ref (RUNG_REF) compared with the same arm at INTEGRATION_REF."""
    out = set()
    for rung in RUNG_REF:
        before = {
            re.split(r"[|\s]", h)[0]: q.upper()
            for h, q in read_fasta(git_show(f"{RUNGS[rung]}.fasta", INTEGRATION_REF))
        }
        out |= {
            (i["rung"], i["key"])
            for i in inserts
            if i["rung"] == rung and before.get(i["id"]) != i["insert"]
        }
    return out


def recoded_check(
    recoded: list[dict], inserts: list[dict], encodings: dict, stale: set | frozenset = frozenset()
) -> tuple[dict, dict]:
    """Every audit invariant on the recoded set, compared insert by insert with the original.

    Inserts in ``stale`` were rebuilt after the order of record, so the order holds an older
    design for them; they are counted, not compared. Recoded inserts of retired members (the
    rescue rung) are counted as still in the order.

    Returns (summary, per-insert recode facts keyed by (rung, id)).
    """
    orig = {(i["rung"], i["key"]): i for i in inserts if (i["rung"], i["key"]) not in stale}
    retired = {(r["rung"], r["key"]) for r in recoded if r["rung"] == "expression rescue"}
    rec = {(r["rung"], r["key"]): r for r in recoded if (r["rung"], r["key"]) not in retired}
    recodes = Counter()
    text = git_show("synthesis_recodes.tsv", RECODES_REF)
    for r in csv.DictReader(io.StringIO(text), delimiter="\t"):
        ident = re.split(r"[|\s]", r["construct"])[0]
        key = r.get("arm_id") or insert_key(ident, r["design_id"])  # per-insert cuts log arm ids
        recodes[(RUNG_OF_TAG[r["rung"]], key)] += 1
    problems = Counter()
    facts = {}
    terminal_blocks = Counter()
    for key, o in orig.items():
        r = rec.get(key)
        if r is None:
            problems["missing from the recoded set"] += 1
            continue
        a, b = r["insert"], o["insert"]
        t = type_iis(a)
        if not t["wrapper_only"]:
            problems["Type IIS site beyond the wrapper"] += 1
        if set(t["extra_type_iis"]) - set(type_iis(b)["extra_type_iis"]):
            problems["new unscreened Type IIS site (AarI/PaqCI/BtgZI/SapI)"] += 1
        ca, pa = coding_of(a)
        cb, _ = coding_of(b)
        if pa:
            problems[pa] += 1
            continue
        prot_a, prot_b = seqtools.translate(ca), seqtools.translate(cb)
        if prot_a != prot_b:
            problems["protein changed"] += 1
        if not protein_matches(prot_a.rstrip("*"), o["row"]):
            problems["protein does not match the arm table"] += 1
        if (
            len(ca) != len(cb)
            or a[len(WRAP_5 + KOZAK) + len(ca) :] != b[len(WRAP_5 + KOZAK) + len(cb) :]
        ):
            problems["readout block, barcode slot or wrapper changed"] += 1
        diffs = [k for k in range(0, len(ca), 3) if ca[k : k + 3] != cb[k : k + 3]]
        nonsyn = [
            k
            for k in diffs
            if seqtools.translate(ca[k : k + 3]) != seqtools.translate(cb[k : k + 3])
        ]
        if nonsyn:
            problems["nonsynonymous codon change"] += 1
        if len(diffs) != recodes[key]:
            problems["codon changes differ from synthesis_recodes.tsv"] += 1
        k = prot_b.find(CASSETTE_AA)
        if k >= 0:
            span = slice(3 * k, 3 * (k + len(CASSETTE_AA)))
            if ca[span] != cb[span]:
                problems["cassette DNA changed"] += 1
            enc = encodings.get(o["row"].get("candidate_id", ""))
            if o["rung"] == "C-terminal/snorkel" and enc and ca[span] != enc:
                problems["C-terminal cassette is not the issued encoding"] += 1
            if o["rung"] == "C-terminal/snorkel":
                terminal_blocks[ca[span][-15:]] += 1
        if SPLICE_DONOR.pattern and len(re.findall("GT[AG]AGT", a)) > len(
            re.findall("GT[AG]AGT", b)
        ):
            problems["new consensus splice donor"] += 1
        facts[key] = {"recoded_codons": len(diffs)}
    dup_ends = sum(n for n in terminal_blocks.values() if n > 1)
    extra = set(rec) - set(orig) - set(stale)
    return (
        {
            "recoded_fasta": f"processed://{RECODED_KEY}",
            "recoded_inserts": len(rec) + len(retired),
            "retired_members_still_in_order": len(retired),
            "rebuilt_since_order_not_compared": sorted(f"{a}:{b}" for a, b in stale),
            "not_in_original_set": len(extra),
            "inserts_recoded": sum(1 for f in facts.values() if f["recoded_codons"]),
            "codons_recoded": sum(f["recoded_codons"] for f in facts.values()),
            "problems": dict(problems),
            "cterm_terminal_linker_blocks_shared": dup_ends,
        },
        facts,
    )


TAIL_AFTER_CODING = SPACER + FWD_PRIMER
_ = (CS1, REV_PRIMER)  # the 3' constants are part of the wrapper checked in type_iis

if __name__ == "__main__":
    raise SystemExit(main())
