"""Issue one DNA encoding of the snorkel module per snorkel arm.

A snorkel arm ends ``native ORF | module | cassette | TAA``, where the module is 55 aa:
``GGSGGSGGSGG`` | PDGFRB 556-566 | PDGFRB 530-555 (the TM helix) | ``GSGGGGS``
(``l1_snorkel.PREFIX_AA``). Until 2026-11-21 every snorkel arm carried one encoding of it,
``l1_snorkel.PREFIX_DNA``, so the 109 snorkel inserts shared up to 176 nt in the middle of the
protein, linkers included. **Decided by Becca, 2026-11-21** ("we still don't want shared linkers
in the middle"; "Yes, recode them" for the snorkel modules).

Every codon of the module is drawn per arm, not only the linkers'. ``PREFIX_DNA`` was a solver
encoding (``seqtools.dnachisel_optimize``, 2026-11-12) that no construct has been shown to
display, unlike p003's tag codons, so holding its TM-helix codons would keep 111 nt of shared
sequence for no tested benefit. The protein is unchanged, so every snorkel refold still applies.
The rules are the cassette's: synonymous codons only, no BbsI/Esp3I/BsaI site in the module or at
either seam, variable 20-mers issued elsewhere in the pool avoided where possible, and arms that
share flanking sequence differ within 15 nt of the junction. How much any module still shares
with any other insert is measured by the pool QC (``l1_pool_qc.py``).

Encodings already in the file are pinned; ``--reissue-all`` redraws them, which is correct only
before any is built. ``l1_swap_risk.py`` keeps ``PREFIX_DNA`` as the flank its cassette batches
were issued against, so replaying those batches is unaffected.

Implemented by Claude Code session l1-review-export-scope, 2026-11-21.

    uv run --locked python diversify_snorkel_modules.py            # issue and write
    uv run --locked python diversify_snorkel_modules.py --check    # replay, write nothing
"""

from __future__ import annotations

import argparse
import csv
import random
import sys
import time
from pathlib import Path

import l1_snorkel
from build_constructs import ENZYMES
from diversify_arm_cassettes import read_fasta
from diversify_cassettes import LINKER_ONLY_JUNCTION_IDENTITY

from deliverome_analysis.cloning import seqtools

HERE = Path(__file__).parent
OUT = HERE / "l1_snorkel_module_encodings.tsv"
SEED = 2026_11_2101
FLANK = 150
FIELDS = ["candidate_id", "arm_id", "module_aa", "module_nt", "gc", "issued_batch"]


def snorkel_arms() -> dict[str, tuple[str, str, str]]:
    """``{candidate_id: (arm_id, flank5, flank3)}`` for every built snorkel arm.

    The module is located by its issued encoding when it has one and by ``PREFIX_DNA``
    otherwise, so the flanks (the native ORF end before it, the cassette after) are the same
    before and after a rebuild.
    """
    rows = {
        r.get("arm_id") or r["candidate_id"]: r
        for r in csv.DictReader((HERE / "constructs_l1_cterm.tsv").open(), delimiter="\t")
        if r["placement"].startswith("snorkel")
    }
    previous = {c: r["module_nt"] for c, r in issued().items()}
    out: dict[str, tuple[str, str, str]] = {}
    for name, insert in read_fasta(HERE / "constructs_l1_cterm.fasta").items():
        if name not in rows:
            continue
        cid = rows[name]["candidate_id"]
        module = previous.get(cid, l1_snorkel.PREFIX_DNA)
        if module not in insert:
            module = l1_snorkel.PREFIX_DNA
        at = insert.find(module)
        assert at >= 0 and insert.find(module, at + 1) < 0, f"{cid}: module not found once"
        end = at + len(module)
        out[cid] = (
            rows[name].get("arm_id", ""),
            insert[max(0, at - FLANK) : at],
            insert[end:][:FLANK],
        )
    assert len(out) == len(rows), "a snorkel arm in the table has no insert"
    return out


def retired_members() -> set[str]:
    with (HERE / "candidate_manifest_tranches.tsv").open() as fh:
        return {
            r["candidate_id"] for r in csv.DictReader(fh, delimiter="\t") if r["retired_reason"]
        }


def issued() -> dict[str, dict[str, str]]:
    if not OUT.exists():
        return {}
    return {r["candidate_id"]: r for r in csv.DictReader(OUT.open(), delimiter="\t")}


def module_dna(candidate_id: str) -> str:
    """The module DNA issued to one snorkel arm (the builder's only source for it).

    A new snorkel arm has none yet: under ``L1_PROVISIONAL_CASSETTE=1`` it is built once on the
    shared ``PREFIX_DNA`` so ``snorkel_arms`` can find its flanks, this script issues its own,
    and the rung is rebuilt without the variable (as ``diversify_arm_cassettes.arm_cassette``).
    """
    import os

    from diversify_arm_cassettes import PROVISIONAL_ENV

    if candidate_id not in issued() and os.environ.get(PROVISIONAL_ENV) == "1":
        return l1_snorkel.PREFIX_DNA
    nt = issued()[candidate_id]["module_nt"]
    assert seqtools.translate(nt) == l1_snorkel.PREFIX_AA, f"{candidate_id}: module protein changed"
    assert nt != l1_snorkel.PREFIX_DNA, f"{candidate_id}: carries the shared module encoding"
    return nt


def issue(
    reissue_all: bool = False,
) -> tuple[dict[str, dict[str, str]], seqtools.LibraryLinkerPlan]:
    previous = {} if reissue_all else issued()
    arms = snorkel_arms()
    # A retired snorkel arm keeps its issued module (written back unchanged); any other
    # issued module with no arm is an error.
    retired = retired_members()
    # A member whose protein changed (l1_arm_reissues.tsv, 2026-11-22) keeps its old snorkel
    # module issued and unused, like a retired member's.
    from diversify_arm_cassettes import REISSUES

    if REISSUES.exists():
        with REISSUES.open() as fh:
            retired |= {r["design_id"] for r in csv.DictReader(fh, delimiter="\t")}
    stale = sorted(set(previous) - set(arms) - retired)
    assert not stale, f"issued modules for arms no longer built: {stale[:5]}"
    members = {c: (f5, f3) for c, (_, f5, f3) in arms.items()}
    existing = {c: r["module_nt"] for c, r in previous.items() if c in arms}
    # A module issued to an arm no longer built (retired, superseded or dropped) stays pinned
    # with no flank, so no new arm can draw it (2026-11-22: four new v8 snorkels had drawn
    # unbuilt members' modules; diversify_arm_cassettes had the same fix that morning).
    for c in sorted(set(previous) - set(arms)):
        members[c], existing[c] = ("N" * FLANK, "N" * FLANK), previous[c]["module_nt"]
    t0 = time.time()
    plan = seqtools.diversify_library_linkers(
        members,
        l1_snorkel.PREFIX_AA,
        kmer_policy="minimize",
        enzymes=ENZYMES,
        max_junction_identity=LINKER_ONLY_JUNCTION_IDENTITY,
        existing=existing or None,
        rng=random.Random(SEED),
    )
    batch = max((int(r["issued_batch"]) for r in previous.values()), default=0) + 1
    rows: dict[str, dict[str, str]] = {
        c: r for c, r in previous.items() if c in retired and c not in arms
    }
    for cid, (arm, _, _) in sorted(arms.items()):
        nt = plan.variants[cid].sequence
        assert seqtools.translate(nt) == l1_snorkel.PREFIX_AA, cid
        rows[cid] = previous.get(cid) or {
            "candidate_id": cid,
            "arm_id": arm,
            "module_aa": l1_snorkel.PREFIX_AA,
            "module_nt": nt,
            "gc": f"{seqtools.gc_content(nt):.3f}",
            "issued_batch": str(batch),
        }
    print(
        f"{len(rows)} snorkel arms: {len(rows) - len(previous)} issued in batch {batch}, "
        f"{len(previous)} pinned; {time.time() - t0:.0f} s"
    )
    return rows, plan


def write(rows: dict[str, dict[str, str]]) -> None:
    with OUT.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(r for _, r in sorted(rows.items()))


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--check", action="store_true", help="replay; fail if the file would change")
    ap.add_argument("--reissue-all", action="store_true")
    a = ap.parse_args(argv)
    rows, plan = issue(a.reissue_all)
    print("audit:", plan.audit())
    distinct = len({r["module_nt"] for r in rows.values()})
    assert distinct == len(rows), f"only {distinct} distinct modules for {len(rows)} arms"
    if a.check:
        same = rows == issued()
        print("matches the record" if same else "DIFFERS from the record")
        return 0 if same else 1
    write(rows)
    print(f"wrote {OUT.name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
