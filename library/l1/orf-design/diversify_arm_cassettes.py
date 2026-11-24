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

**Adding arms.** A new arm has no encoding yet, and this script finds an arm's flanks in its
built insert. So: build the rung once with ``L1_PROVISIONAL_CASSETTE=1`` (new arms get p003's
encoding, existing arms their issued one), run this script (it issues only the new arms'
encodings, every issued one pinned), then rebuild the rung without the variable.

**A member whose protein changed** (a reference moved to the GTF MANE Select protein, 2026-11-22)
keeps its old arm's encoding issued but unused, and its new arm gets a new one. Its upstream arm
key (the insert's FASTA id) is mapped to ``<old key>@<new protein sha256[:12]>`` by
``l1_arm_reissues.tsv`` (``--supersede FILE`` adds the rows for a table with ``candidate_id`` and
``new_protein_sha256``); ``arm_cassette`` and the flank search read the mapped key, and the old
key's row is written back unchanged, like a retired member's.
"""

from __future__ import annotations

import argparse
import csv
import os
import random
import re
import sys
import time
from collections import Counter
from pathlib import Path

from build_constructs import CASSETTE, CASSETTE_AA, ENZYMES
from diversify_cassettes import LINKER_ONLY_JUNCTION_IDENTITY

from deliverome_analysis import plasmid_qc as qc
from deliverome_analysis import seqtools

HERE = Path(__file__).parent
OUT = HERE / "l1_arm_cassette_encodings.tsv"
REISSUES = HERE / "l1_arm_reissues.tsv"
REISSUE_FIELDS = (
    "old_key", "new_key", "design_id", "new_protein_sha256", "reason", "decided_by",
    "implemented_by",
)  # fmt: skip


def reissues() -> dict[str, str]:
    """``{old arm key: new arm key}`` for arms whose member's protein changed."""
    if not REISSUES.exists():
        return {}
    with REISSUES.open() as fh:
        return {r["old_key"]: r["new_key"] for r in csv.DictReader(fh, delimiter="\t")}


def held_arms() -> set[str]:
    """Upstream arm keys issued an encoding and then held by their rung (a refold failure the
    rung's rule does not build, 2026-11-22): the encoding stays issued, unused."""
    held = set()
    with (HERE / "reviewer-package/nterm/nterm_decisions.tsv").open() as fh:
        held |= {
            r["candidate_id"] for r in csv.DictReader(fh, delimiter="\t")
            if r["arm"] == "no_nterm_arm"
        }  # fmt: skip
    arms = HERE.parents[1] / "sprint5/s5e13-tag-site-benchmarking/l1_internal_arm_decisions.tsv"
    with arms.open() as fh:
        held |= {r["arm_id"] for r in csv.DictReader(fh, delimiter="\t") if r["decision"] == "hold"}
    return held


def encoding_key(arm: str) -> str:
    # an internal arm is superseded with no new key: its new arms carry new arm ids instead
    return reissues().get(arm) or arm


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
    return seqtools.fasta_by_id(path)


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
            arm = encoding_key(arm)
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
#: Set to "1" only for the provisional build of new arms (see arm_cassette).
PROVISIONAL_ENV = "L1_PROVISIONAL_CASSETTE"


def arm_cassette(arm_id: str) -> str:
    """The cassette DNA issued to one upstream arm (the builders' only source for it).

    Raises ``KeyError`` for an arm with no encoding: run this script first. There is no fallback
    to p003's shared encoding, which is what let 2,942 arms share it.
    """
    if not _ISSUED:
        _ISSUED.update({a: r["cassette_nt"] for a, r in issued().items()})
    arm_id = encoding_key(arm_id)
    if arm_id not in _ISSUED and os.environ.get(PROVISIONAL_ENV) == "1":
        # A new arm: build it once on p003's encoding so this script can find its flanks,
        # issue its own encoding, then rebuild without the variable. Never for an order.
        return CASSETTE
    nt = _ISSUED[arm_id]
    assert seqtools.translate(nt) == CASSETTE_AA, f"{arm_id}: not the 51-aa cassette"
    assert qc.keeps_p003_tag_codons(nt), f"{arm_id}: not p003's tags"
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
    # An arm of a retired member keeps its issued encoding (written back unchanged, never handed
    # on); any other issued encoding with no arm is an error.
    with (HERE / "candidate_manifest_tranches.tsv").open() as fh:
        retired = {
            r["candidate_id"] for r in csv.DictReader(fh, delimiter="\t") if r["retired_reason"]
        }
    gone = {a for a in set(previous) - set(arms)}
    member = {a: a.split(":")[0].rsplit("-a", 1)[0] for a in gone}
    superseded = set(reissues())  # an old protein's arm: kept issued, unused
    stale = sorted(
        a for a in gone if member[a] not in retired and a not in superseded and a not in held_arms()
    )
    assert not stale, f"issued encodings for arms no longer built: {stale[:5]}"
    kept_retired = {a: previous[a] for a in gone}
    cflanks, cpinned = cterm_members()
    members = {a: (f5, f3) for a, (_, f5, f3) in arms.items()} | cflanks
    existing = {a: r["cassette_nt"] for a, r in previous.items() if a in arms} | cpinned
    # An arm no longer built (a retired member's, or superseded: its member's protein changed)
    # keeps its issued encoding: pinned with no flank, so its k-mers are in the pool and no new
    # arm can draw it (2026-11-22: four new arms had drawn retired members' encodings).
    for a in sorted(gone):
        members[a], existing[a] = ("N" * FLANK, "N" * FLANK), previous[a]["cassette_nt"]
    fixed = qc.p003_tag_codons()
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
    rows: dict[str, dict[str, str]] = dict(kept_retired)
    for arm, (rung, _, _) in sorted(arms.items()):
        nt = plan.variants[arm].sequence
        assert seqtools.translate(nt) == CASSETTE_AA, arm
        assert qc.keeps_p003_tag_codons(nt), f"{arm}: tag codons"
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
        w.writerows(r for _, r in sorted(rows.items()))


def supersede(path: Path, reason: str, decided_by: str, implemented_by: str) -> int:
    """Add reissue rows for every issued upstream arm of the members in ``path`` whose protein
    is now ``new_protein_sha256`` (the arm keys are found in the issued encodings)."""
    have = reissues()
    rows = list(csv.DictReader(REISSUES.open(), delimiter="\t")) if REISSUES.exists() else []
    old = issued()
    with path.open() as fh:
        changes = {
            r["candidate_id"]: r["new_protein_sha256"] for r in csv.DictReader(fh, delimiter="\t")
        }
    added = 0
    for key in sorted(old):
        if "@" in key:
            continue
        cid = key.split(":")[0]
        internal = bool(ARM_RE.fullmatch(cid))  # new arms get new arm ids (per-arm table)
        if internal:
            cid = cid.rsplit("-a", 1)[0]
        if cid not in changes or key in have:
            continue
        rows.append({
            "old_key": key, "new_key": "" if internal else f"{key}@{changes[cid][:12]}",
            "design_id": cid,
            "new_protein_sha256": changes[cid], "reason": reason, "decided_by": decided_by,
            "implemented_by": implemented_by,
        })  # fmt: skip
        added += 1
    with REISSUES.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=REISSUE_FIELDS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(f"{REISSUES.name}: {added} arm(s) superseded ({len(rows)} in all)")
    return 0


ARM_RE = re.compile(r"l[12]-p[0-9]{4}-a[0-9]+")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--check", action="store_true", help="replay; fail if the file would change")
    ap.add_argument("--reissue-all", action="store_true")
    ap.add_argument("--supersede", type=Path, help="table of members whose protein changed")
    a = ap.parse_args(argv)
    if a.supersede:
        return supersede(
            a.supersede,
            "the member's reference moved to the GTF MANE Select protein (2026-11-22); the old "
            "arm's encoding stays issued, unused, and the new arm is issued its own",
            "Becca Carlson, 2026-11-22 ('Yes, in this build')",
            "Claude Code session l1-review-export-scope (subagent), 2026-11-22",
        )
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
