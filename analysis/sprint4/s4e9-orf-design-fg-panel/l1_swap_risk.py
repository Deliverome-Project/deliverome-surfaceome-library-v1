"""Barcode-swap substrate across every L1 member, under each cassette-encoding scheme.

``check_swap_risk.py`` gates the panel constructs that have been built. This is the same
question for the library: for each L1 member, how long is the tract it shares, up to the
barcode, with the member it shares most with? A reverse-transcription template switch
inside that tract moves a barcode onto another member's ORF unseen (Walton's fit,
``seqtools.barcode_swap_rate``). Gate <= 100 bp, <= 80 preferred.

Every L1 insert ends ``... | TAA | CC | 19-nt handle | barcode``, so 24 nt is the floor.
What sits before the stop depends on placement:

* **C-terminal** members end ``ORF | 51-aa cassette | TAA``; the cassette is the last
  coding sequence and its encoding decides the tract.
* **N-terminal and internal** members end on the ORF's own C-terminus; the cassette is
  upstream and does not enter the tract. Isoforms that keep their gene's C-terminus share
  that tail with it, which the 2026-10-13 decision accepts.

Three schemes for the C-terminal cassette are measured on the same members:

``shared``
    p003's encoding in every member (the starting point).
``whole_cassette``
    every codon recoded per member, as ``diversify_cassettes.py`` drew the 25-construct
    order (2026-10-07), with its 9-nt end-uniqueness rule.
``whole_cassette_widened``
    every codon recoded, but under the 15-nt end rule -- the control that says whether
    recoding the tags buys anything once the end rule is wide enough to finish.
``linkers_only``
    SpyTag003 and P-ALFA-P held on p003's codons; only the four GGGGS blocks vary
    (Becca, 2026-10-29). The issued encodings are written to ``l1_cassette_encodings.tsv``.

Members come from ``candidate_manifest_tranches.tsv`` (retired duplicates excluded, see
``dedupe_l1_members.py``), placement from ``l1_tag_paths.tsv``, sequence from the designed
CDS on S3.

    uv run --locked python l1_swap_risk.py
"""

from __future__ import annotations

import argparse
import csv
import json
import random
import time
from collections import Counter
from pathlib import Path

import numpy as np
from build_constructs import CASSETTE, CASSETTE_AA, ENZYMES, FWD_PRIMER, SPACER, STOP
from diversify_cassettes import (
    LINKER_ONLY_JUNCTION_IDENTITY,
    TAG_SEGMENTS,
    WHOLE_CASSETTE_JUNCTION_IDENTITY,
)

from deliverome_analysis import library_design as ld
from deliverome_analysis import seqtools

HERE = Path(__file__).parent
MANIFEST = HERE / "candidate_manifest_tranches.tsv"
TAG_PATHS = HERE / "l1_tag_paths.tsv"
ENCODINGS = HERE / "l1_cassette_encodings.tsv"
OUT_TSV = HERE / "l1_swap_risk.tsv"
OUT_JSON = HERE / "l1_swap_risk_summary.json"

TAIL = STOP + SPACER + FWD_PRIMER  # the constant 24 nt before every barcode
GATE, PREFERRED = 100, 80
SEED = 2026_10_29
FLANK = 150


def common_suffix(a: str, b: str) -> int:
    n = 0
    while n < min(len(a), len(b)) and a[-1 - n] == b[-1 - n]:
        n += 1
    return n


def worst_partner(upstream: list[str]) -> tuple[np.ndarray, np.ndarray, int]:
    """Each member's longest common suffix with any other, that partner, and pairs over GATE.

    Sorting the reversed sequences puts every member's best suffix partner next to it, and
    the suffix two members share is the minimum adjacent LCP between them.
    """
    order = sorted(range(len(upstream)), key=lambda i: upstream[i][::-1])
    adj = [common_suffix(upstream[order[k]], upstream[order[k + 1]]) for k in range(len(order) - 1)]
    worst = np.zeros(len(upstream), dtype=int)
    partner = np.full(len(upstream), -1)
    for k, lcp in enumerate(adj):
        a, b = order[k], order[k + 1]
        if lcp > worst[a]:
            worst[a], partner[a] = lcp, b
        if lcp > worst[b]:
            worst[b], partner[b] = lcp, a
    pairs_over = run = 0
    for lcp in adj:
        run = run + 1 if lcp > GATE else 0
        pairs_over += run
    return worst, partner, pairs_over


def domesticate_junction(body: str) -> tuple[str, str]:
    """Recode an ORF's last codon synonymously if it completes a Type IIS site with the cassette.

    The cassette opens on Gly, which is GGN, so an ORF ending in the five bases that
    complete a site (SLC41A2 ends ``...GGA GAC``, and ``GAC GG`` closes Esp3I) cannot be
    rescued by any cassette codon choice. Returns the body and a description of the edit,
    or an empty string when no edit was needed.
    """
    sites = seqtools.forbidden_patterns(ENZYMES)

    def clashes(seq: str) -> bool:
        return any(p in seq[-5:] + g for p in sites for g in ("GGA", "GGC", "GGG", "GGT"))

    if not clashes(body):
        return body, ""
    last = body[-3:]
    for codon, _ in seqtools.synonymous_codons()[seqtools.translate(last)]:
        if codon != last and not clashes(body[:-3] + codon):
            return body[:-3] + codon, f"last codon {last}->{codon}"
    raise SystemExit(f"no synonymous last codon clears the cassette junction for ...{body[-9:]}")


def load_members() -> tuple[list[dict], list[str]]:
    """Members with a placement, and the ids of active members that have none yet.

    A member with no DeepTMHMM v1 topology has no tag path, so whether its cassette sits
    before the barcode is unknown; it is reported, not guessed.
    """
    proteins = ld.designed_proteins(MANIFEST, library=ld.LIBRARY_L1)
    keys = set(proteins)
    cds = ld.designed_cds(MANIFEST, library=ld.LIBRARY_L1)
    out = []
    for r in csv.DictReader(TAG_PATHS.open(), delimiter="\t"):
        cid = r["candidate_id"]
        if cid not in keys:
            continue
        seq = cds[cid]
        body = seq[:-3] if seqtools.translate(seq[-3:]) == "*" else seq
        edit = ""
        if r["placement"] == "C-terminal":
            body, edit = domesticate_junction(body)
        out.append(
            {
                "id": cid,
                "gene": r["gene_symbol"],
                "set": r["candidate_set"],
                "placement": r["placement"],
                "cterm": r["placement"] == "C-terminal",
                "body": body,
                "junction_edit": edit,
            }
        )
    return out, sorted(keys - {m["id"] for m in out})


def issued_batches() -> dict[str, tuple[str, int]]:
    """``{candidate_id: (cassette_nt, issued_batch)}`` from ``l1_cassette_encodings.tsv``.

    Encodings are pinned once issued, exactly as ``diversify_cassettes.py`` pins the panel's:
    members are served in risk order, so without pinning, adding two members re-drew 156 of
    216 encodings. Each records the batch it was issued in, so ``--verify`` can replay the
    issuance: batch k is drawn over every member of batches 1..k with batches 1..k-1 pinned.
    """
    if not ENCODINGS.exists():
        return {}
    with ENCODINGS.open() as fh:
        return {
            r["candidate_id"]: (r["cassette_nt"], int(r["issued_batch"]))
            for r in csv.DictReader(fh, delimiter="\t")
        }


def issued_encodings(flanks: dict) -> dict[str, str]:
    """Issued encodings for members still drawn."""
    return {cid: nt for cid, (nt, _) in issued_batches().items() if cid in flanks}


def draw(
    members: list[dict],
    *,
    fixed: dict[int, str],
    identity: int,
    pin: bool = False,
    existing: dict[str, str] | None = None,
) -> tuple[dict, str]:
    flanks = {m["id"]: (m["body"][-FLANK:], TAIL) for m in members if m["cterm"]}
    if existing is None:
        existing = issued_encodings(flanks) if pin else {}
    t0 = time.time()
    try:
        plan = seqtools.diversify_library_linkers(
            flanks,
            CASSETTE_AA,
            kmer_policy="minimize",
            enzymes=ENZYMES,
            fixed_codons=fixed,
            max_junction_identity=identity,
            existing=existing or None,
            rng=random.Random(SEED),
        )
    except RuntimeError as exc:
        issued = str(exc).split(" ran out of encoding space after ")[-1].split(" of ")[0]
        return {}, f"FAILED after {issued} of {len(flanks)} C-terminal members: {exc}"
    return plan, (
        f"issued {len(plan.variants) - len(existing)} new, pinned {len(existing)}, "
        f"in {time.time() - t0:.0f} s"
    )


def summarize(members: list[dict], upstream: list[str]) -> tuple[dict, np.ndarray, np.ndarray]:
    worst, partner, pairs_over = worst_partner(upstream)
    cterm = np.array([m["cterm"] for m in members])
    q = np.percentile(worst, [50, 90, 99])
    imax = int(np.argmax(worst))
    over = [i for i in range(len(members)) if worst[i] > GATE]
    return (
        {
            "members": len(members),
            "worst_partner_tract_median": float(q[0]),
            "p90": float(q[1]),
            "p99": float(q[2]),
            "max": int(worst[imax]),
            "max_pair": [members[imax]["id"], members[partner[imax]]["id"]],
            "members_at_floor_24nt": int((worst == len(TAIL)).sum()),
            "members_over_preferred_80": int((worst > PREFERRED).sum()),
            "members_over_gate_100": len(over),
            "pairs_over_gate_100": pairs_over,
            "members_over_gate_by_set": dict(Counter(members[i]["set"] for i in over)),
            "cterm_members": int(cterm.sum()),
            "cterm_worst_partner_median": float(np.median(worst[cterm])) if cterm.any() else None,
            "cterm_worst_partner_max": int(worst[cterm].max()) if cterm.any() else None,
            "swap_rate_at_median": round(seqtools.barcode_swap_rate(int(q[0])), 4),
            "swap_rate_at_max": round(seqtools.barcode_swap_rate(int(worst[imax])), 4),
        },
        worst,
        partner,
    )


def verify_issuance() -> str:
    """Replay every issued batch from scratch and require the recorded encodings back."""
    members, _ = load_members()
    tag_fixed = seqtools.pin_codons(CASSETTE, CASSETTE_AA, TAG_SEGMENTS)
    batches = issued_batches()
    cterm = {m["id"] for m in members if m["cterm"]}
    missing = sorted(cterm - set(batches))
    assert not missing, f"C-terminal members with no issued encoding: {missing[:5]}"
    for k in sorted({b for _, b in batches.values()}):
        upto = {cid for cid, (_, b) in batches.items() if b <= k}
        pins = {cid: nt for cid, (nt, b) in batches.items() if b < k}
        drawn = [m for m in members if not m["cterm"] or m["id"] in upto]
        plan, note = draw(
            drawn,
            fixed=tag_fixed,
            identity=LINKER_ONLY_JUNCTION_IDENTITY,
            existing=pins,
        )
        assert plan, f"batch {k}: {note}"
        for cid in upto - set(pins):
            assert plan.variants[cid].sequence == batches[cid][0], f"batch {k}: {cid} differs"
    n = len({b for _, b in batches.values()})
    return f"all {len(batches)} linker encodings replay from scratch in {n} issued batch(es)"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--no-write", action="store_true", help="print the summary only")
    ap.add_argument("--verify", action="store_true", help="replay the issued batches and stop")
    args = ap.parse_args(argv)
    if args.verify:
        print(verify_issuance())
        return 0

    members, unplaced = load_members()
    if unplaced:
        print(f"{len(unplaced)} active member(s) have no tag path and are not measured: {unplaced}")
    tag_fixed = seqtools.pin_codons(CASSETTE, CASSETTE_AA, TAG_SEGMENTS)
    report: dict = {
        "members": len(members),
        "unplaced_not_measured": unplaced,
        "placement": dict(Counter(m["placement"] for m in members)),
        "gate_bp": GATE,
        "preferred_bp": PREFERRED,
        "seed": SEED,
        "cterm_junction_edits": {
            m["id"]: f"{m['gene']}: {m['junction_edit']}" for m in members if m["junction_edit"]
        },
    }
    scenarios: dict[str, dict[str, str]] = {"shared": {m["id"]: CASSETTE for m in members}}
    for name, identity in (
        ("whole_cassette", WHOLE_CASSETTE_JUNCTION_IDENTITY),
        ("whole_cassette_widened", LINKER_ONLY_JUNCTION_IDENTITY),
    ):
        whole, note = draw(members, fixed={}, identity=identity)
        report[f"{name}_draw"] = note
        if whole:
            scenarios[name] = {k: v.sequence for k, v in whole.variants.items()}
    linkers, note = draw(members, fixed=tag_fixed, identity=LINKER_ONLY_JUNCTION_IDENTITY, pin=True)
    report["linkers_only_draw"] = note
    if not linkers:
        raise SystemExit(f"linkers-only draw failed: {note}")
    scenarios["linkers_only"] = {k: v.sequence for k, v in linkers.variants.items()}
    report["linkers_only_audit"] = linkers.audit()

    per_member: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for name, enc in scenarios.items():
        upstream = [m["body"] + (enc[m["id"]] if m["cterm"] else "") + TAIL for m in members]
        report[name], worst, partner = summarize(members, upstream)
        per_member[name] = (worst, partner)
        print(f"\n{name}:")
        for key in (
            "worst_partner_tract_median",
            "max",
            "max_pair",
            "members_at_floor_24nt",
            "members_over_gate_100",
            "cterm_worst_partner_median",
            "cterm_worst_partner_max",
        ):
            print(f"  {key:<30} {report[name][key]}")
    for name in ("whole_cassette", "whole_cassette_widened", "linkers_only"):
        print(f"{name} draw: {report[f'{name}_draw']}")

    if args.no_write:
        return 0
    batches = issued_batches()
    new_batch = max((b for _, b in batches.values()), default=0) + 1
    with OUT_TSV.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        cols = ["candidate_id", "gene_symbol", "candidate_set", "placement"]
        for name in per_member:
            cols += [f"tract_nt_{name}", f"partner_{name}"]
        w.writerow(cols)
        for i, m in enumerate(members):
            row = [m["id"], m["gene"], m["set"], m["placement"]]
            for worst, partner in per_member.values():
                row += [int(worst[i]), members[partner[i]]["id"] if partner[i] >= 0 else ""]
            w.writerow(row)
    with ENCODINGS.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["candidate_id", "scheme", "cassette_aa", "cassette_nt", "gc", "issued_batch"])
        for cid in sorted(scenarios["linkers_only"]):
            dna = scenarios["linkers_only"][cid]
            assert seqtools.translate(dna) == CASSETTE_AA
            assert all(dna[3 * i : 3 * i + 3] == c for i, c in tag_fixed.items())
            batch = batches[cid][1] if cid in batches else new_batch
            w.writerow(
                [cid, "linkers_only", CASSETTE_AA, dna, f"{seqtools.gc_content(dna):.3f}", batch]
            )
    OUT_JSON.write_text(json.dumps(report, indent=2, default=str) + "\n")
    print(f"\nwrote {OUT_TSV.name}, {ENCODINGS.name}, {OUT_JSON.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
