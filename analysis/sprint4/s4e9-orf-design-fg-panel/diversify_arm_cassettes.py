"""Issue one linker encoding of the tag cassette per N-terminal, internal, multipass and type III arm.

Until 2026-11-21 only the C-terminal rung had per-member encodings (``l1_cassette_encodings.tsv``,
from ``l1_swap_risk.py``): the swap model measures the tract a member shares up to its barcode,
and a cassette upstream of the ORF's own C-terminus never enters it. Every other arm carried
p003's one encoding (``build_constructs.CASSETTE``), so 2,942 inserts shared the same 153 nt.
A template switch inside that stretch joins the front of one ORF to the back and barcode of
another, and for an internal arm that is a chimera of two proteins under one barcode.

**Decided by Becca, 2026-11-21:** "we still don't want shared linkers in the middle." So every
arm gets its own encoding under the same rules as the C-terminal one (Becca, 2026-10-29):

* SpyTag003 and P-ALFA-P keep p003's codons (``TAG_SEGMENTS``); only the four GGGGS blocks vary.
* No two arms that share flanking sequence may share the first or last 15 nt of the cassette
  (``LINKER_ONLY_JUNCTION_IDENTITY`` = 14). This is what separates the 609 arms whose cassette
  follows Met1 (multipass and type III), which all share the Kozak + ATG before it.
* The 574 arms whose cassette follows the Igκ leader + Asp + Met1 share the 95 nt before it, so
  each needs its own first GGGGS block, and only 540 of those blocks have no run of 5 bases. They
  are drawn first with runs of up to 5 allowed (``CROWDED_MAX_HOMOPOLYMER``; GGGGG at most, which
  the synthesis rules in ``synthesis_compliance.py`` accept); every other arm keeps the default 4.
* Variable 20-mers already issued anywhere in the pool are avoided where possible
  (``kmer_policy="minimize"``), and no BbsI / Esp3I / BsaI site may appear in the cassette or at
  either seam with the arm's own flanks.

The 313 C-terminal encodings are passed in pinned, with their own flanks, so the new draws avoid
them and they never change. Encodings issued here are pinned too: rerunning draws only arms that
have none. ``--reissue-all`` redraws this file's arms, which is correct only before any is built.

Implemented by Claude Code session l1-review-export-scope, 2026-11-21.

    uv run --locked python diversify_arm_cassettes.py            # issue and write
    uv run --locked python diversify_arm_cassettes.py --check    # replay, write nothing
"""

from __future__ import annotations

import argparse
import csv
import random
import sys
import time
from collections import Counter
from pathlib import Path

from build_constructs import CASSETTE, CASSETTE_AA, ENZYMES
from diversify_cassettes import LINKER_ONLY_JUNCTION_IDENTITY, TAG_SEGMENTS

from deliverome_analysis import seqtools

HERE = Path(__file__).parent
OUT = HERE / "l1_arm_cassette_encodings.tsv"
CTERM_ENCODINGS = HERE / "l1_cassette_encodings.tsv"
#: Rungs whose cassette sits upstream of the ORF's own C-terminus, and their construct FASTAs.
RUNGS = ("nterm", "internal", "multipass", "type3")
SEED = 2026_11_21
FLANK = 150
#: A group of arms sharing their 5' flank larger than this cannot each get a distinct first
#: GGGGS block under the default run limit (540 blocks), so it is drawn with runs up to 5.
CROWDED_GROUP = 400
CROWDED_MAX_HOMOPOLYMER = 5
FIELDS = ["arm_id", "rung", "scheme", "cassette_aa", "cassette_nt", "gc", "issued_batch"]


def read_fasta(path: Path) -> dict[str, str]:
    """``{arm_id: insert}``; the arm id is the header up to the first ``|`` or space."""
    out: dict[str, str] = {}
    name = ""
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split("|")[0].split()[0]
            assert name not in out, f"{path.name}: duplicate {name}"
            out[name] = ""
        else:
            out[name] += line.strip()
    return out


def cassette_span(insert: str, cassette: str) -> int:
    """Start of ``cassette`` in ``insert``; it must occur exactly once."""
    at = insert.find(cassette)
    assert at >= 0 and insert.find(cassette, at + 1) < 0, "cassette not found exactly once"
    return at


def upstream_arms(encodings: dict[str, str] | None = None) -> dict[str, tuple[str, str, str]]:
    """``{arm_id: (rung, flank5, flank3)}`` for every upstream arm.

    The arm's cassette is located by its issued encoding when it has one (after a rebuild) and
    by p003's shared encoding otherwise (before one), so the flanks are the same either way.
    """
    encodings = encodings or {}
    out: dict[str, tuple[str, str, str]] = {}
    for rung in RUNGS:
        for arm, insert in read_fasta(HERE / f"constructs_l1_{rung}.fasta").items():
            cas = encodings.get(arm, CASSETTE)
            if cas not in insert:
                cas = CASSETTE
            at = cassette_span(insert, cas)
            assert arm not in out, f"{arm} is on two rungs"
            out[arm] = (rung, insert[max(0, at - FLANK) : at], insert[at + len(cas) :][:FLANK])
    return out


def cterm_members() -> tuple[dict[str, tuple[str, str]], dict[str, str]]:
    """Flanks and pinned encodings of the built C-terminal arms, keyed ``cterm:<arm>``."""
    issued = {
        r["candidate_id"]: r["cassette_nt"]
        for r in csv.DictReader(CTERM_ENCODINGS.open(), delimiter="\t")
    }
    flanks: dict[str, tuple[str, str]] = {}
    pinned: dict[str, str] = {}
    table = {
        r.get("arm_id") or r["candidate_id"]: r["candidate_id"]
        for r in csv.DictReader((HERE / "constructs_l1_cterm.tsv").open(), delimiter="\t")
    }
    for arm, insert in read_fasta(HERE / "constructs_l1_cterm.fasta").items():
        nt = issued[table[arm]]
        at = cassette_span(insert, nt)
        key = f"cterm:{arm}"
        flanks[key] = (insert[max(0, at - FLANK) : at], insert[at + len(nt) :][:FLANK])
        pinned[key] = nt
    return flanks, pinned


_ISSUED: dict[str, str] = {}


def arm_cassette(arm_id: str) -> str:
    """The cassette DNA issued to one upstream arm (the builders' only source for it).

    Raises ``KeyError`` for an arm with no encoding: run this script first. There is no fallback
    to p003's shared encoding, which is what let 2,942 arms share it.
    """
    if not _ISSUED:
        _ISSUED.update({a: r["cassette_nt"] for a, r in issued().items()})
    nt = _ISSUED[arm_id]
    fixed = seqtools.pin_codons(CASSETTE, CASSETTE_AA, TAG_SEGMENTS)
    assert seqtools.translate(nt) == CASSETTE_AA, f"{arm_id}: not the 51-aa cassette"
    assert all(nt[3 * i : 3 * i + 3] == c for i, c in fixed.items()), f"{arm_id}: not p003's tags"
    assert nt != CASSETTE, f"{arm_id}: carries p003's shared encoding"
    return nt


def issued() -> dict[str, dict[str, str]]:
    if not OUT.exists():
        return {}
    return {r["arm_id"]: r for r in csv.DictReader(OUT.open(), delimiter="\t")}


def issue(
    reissue_all: bool = False,
) -> tuple[dict[str, dict[str, str]], seqtools.LibraryLinkerPlan]:
    previous = {} if reissue_all else issued()
    arms = upstream_arms({a: r["cassette_nt"] for a, r in previous.items()})
    stale = sorted(set(previous) - set(arms))
    assert not stale, f"issued encodings for arms no longer built: {stale[:5]}"
    cflanks, cpinned = cterm_members()
    members = {a: (f5, f3) for a, (_, f5, f3) in arms.items()} | cflanks
    existing = {a: r["cassette_nt"] for a, r in previous.items()} | cpinned
    fixed = seqtools.pin_codons(CASSETTE, CASSETTE_AA, TAG_SEGMENTS)
    groups = Counter(f5[-60:] for f5, _ in members.values())
    crowded = {a for a, (_, f5, _) in arms.items() if groups[f5[-60:]] > CROWDED_GROUP}
    t0 = time.time()
    common = dict(
        kmer_policy="minimize",
        enzymes=ENZYMES,
        fixed_codons=fixed,
        max_junction_identity=LINKER_ONLY_JUNCTION_IDENTITY,
    )
    # Crowded arms first, with C-terminal and earlier issues pinned; then everyone else.
    first = seqtools.diversify_library_linkers(
        {a: members[a] for a in crowded | set(existing)},
        CASSETTE_AA,
        existing=existing,
        rng=random.Random(SEED),
        max_homopolymer=CROWDED_MAX_HOMOPOLYMER,
        **common,
    )
    existing = existing | {a: first.variants[a].sequence for a in crowded}
    plan = seqtools.diversify_library_linkers(
        members, CASSETTE_AA, existing=existing, rng=random.Random(SEED), **common
    )
    batch = max((int(r["issued_batch"]) for r in previous.values()), default=0) + 1
    rows: dict[str, dict[str, str]] = {}
    for arm, (rung, _, _) in sorted(arms.items()):
        nt = plan.variants[arm].sequence
        assert seqtools.translate(nt) == CASSETTE_AA, arm
        assert all(nt[3 * i : 3 * i + 3] == c for i, c in fixed.items()), f"{arm}: tag codons"
        rows[arm] = previous.get(arm) or {
            "arm_id": arm,
            "rung": rung,
            "scheme": "linkers_only",
            "cassette_aa": CASSETTE_AA,
            "cassette_nt": nt,
            "gc": f"{seqtools.gc_content(nt):.3f}",
            "issued_batch": str(batch),
        }
    print(
        f"{len(rows)} upstream arms: {len(rows) - len(previous)} issued in batch {batch}, "
        f"{len(previous)} pinned, {len(cpinned)} C-terminal pinned; {len(crowded)} crowded "
        f"(runs up to {CROWDED_MAX_HOMOPOLYMER}); {time.time() - t0:.0f} s"
    )
    return rows, plan


def write(rows: dict[str, dict[str, str]]) -> None:
    with OUT.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows.values())


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--check", action="store_true", help="replay; fail if the file would change")
    ap.add_argument("--reissue-all", action="store_true")
    a = ap.parse_args(argv)
    rows, plan = issue(a.reissue_all)
    audit = plan.audit()
    print("audit:", audit)
    distinct = len({r["cassette_nt"] for r in rows.values()})
    assert distinct == len(rows), f"only {distinct} distinct encodings for {len(rows)} arms"
    if a.check:
        same = rows == issued()
        print("matches the record" if same else "DIFFERS from the record")
        return 0 if same else 1
    write(rows)
    print(f"wrote {OUT.name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
