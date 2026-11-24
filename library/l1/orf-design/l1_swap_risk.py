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

Members come from ``candidate_manifest_tranches.tsv`` (retired members excluded, see
``dedupe_l1_members.py``), placement from ``l1_tag_paths.tsv``, sequence from the designed
CDS on S3. A retired member keeps the encoding it was issued: it stays in
``l1_cassette_encodings.tsv``, new batches are drawn around it, and it is not measured.

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

import l1_snorkel
import numpy as np
from build_constructs import CASSETTE, CASSETTE_AA, ENZYMES
from diversify_cassettes import (
    LINKER_ONLY_JUNCTION_IDENTITY,
    WHOLE_CASSETTE_JUNCTION_IDENTITY,
)

from deliverome_analysis import library_design as ld
from deliverome_analysis import plasmid_qc as qc
from deliverome_analysis import seqtools

HERE = Path(__file__).parent
MANIFEST = HERE / "candidate_manifest_tranches.tsv"
TAG_PATHS = HERE / "l1_tag_paths.tsv"
ENCODINGS = HERE / "l1_cassette_encodings.tsv"
OUT_TSV = HERE / "l1_swap_risk.tsv"
OUT_JSON = HERE / "l1_swap_risk_summary.json"

#: The constant run before every barcode: stop, pad, T7 promoter, 5' handle (44 nt). It was
#: stop, pad, handle (24 nt) until the L1 readout region (Becca, 2026-11-22: "add the t7 and the
#: feldman flank"); plasmid_qc.STOP_TO_BARCODE == plasmid_qc.L1_BEFORE_BARCODE.
TAIL = qc.STOP_TO_BARCODE
GATE, PREFERRED = 100, 80
SEED = 2026_10_29
FLANK = 150


common_suffix = seqtools.common_suffix


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


#: Active members that hold an issued encoding but are no longer C-terminal, keyed by design
#: id, with the designed-CDS record their encoding was drawn against. Their encoding stays
#: issued (it is never redrawn or handed on), so later batches are drawn around it and a
#: replay of its batch needs the flank it was issued with. TSPAN10 (l1-p2284): issued in
#: batch 1 on the 235-aa APPRIS model; the 2026-11-12 reference swap to Q9H1Z9 (decided by
#: l1-proteoform-audit under Becca's delegation, applied by l1-reviewer-package, PR #283)
#: made it internal. SLC41A2 (l1-p1957): issued in batch 1, then moved to the residue-1
#: multipass rung on 2026-11-14 because experimental topology puts its C-terminus inside
#: (lane_overrides.tsv; l1-nterm-res1-multipass under Becca's delegation). Its CDS is unchanged
#: since batch 1, so the base record is the flank it was drawn on.
RELEASED = {
    "l1-p2284": "sprint4/s4e9-orf-design-fg-panel/library-candidates/"
    "2026-11-12-reference-fixes/candidate_cds.json",
    "l1-p1957": "sprint4/s4e9-orf-design-fg-panel/library-candidates/candidate_cds.json",
    # Rescue round, 2026-11-15 (coordinator under Becca's delegation; lane_overrides.tsv): held
    # C-terminal members moved to the internal rung with loop arms. Each keeps its batch-1
    # encoding as a holder, drawn on the base record's CDS (--verify replays batch 1 on it).
    "l1-p2342": "sprint4/s4e9-orf-design-fg-panel/library-candidates/candidate_cds.json",  # VTCN1
    "l1-p2628": "sprint4/s4e9-orf-design-fg-panel/library-candidates/candidate_cds.json",  # VTCN1
    "l1-p2487": "sprint4/s4e9-orf-design-fg-panel/library-candidates/candidate_cds.json",  # PLVAP
    "l1-p2528": "sprint4/s4e9-orf-design-fg-panel/library-candidates/candidate_cds.json",  # TRPV1
}


def load_members(with_retired_holders: bool = False) -> tuple[list[dict], list[str]]:
    """Members with a placement, and the ids of active members that have none yet.

    A member with no DeepTMHMM v1 topology has no tag path, so whether its cassette sits
    before the barcode is unknown; it is reported, not guessed.

    ``with_retired_holders`` adds every retired member that holds an issued encoding, as a
    C-terminal member (``retired``: True), in its manifest position. Those encodings stay
    issued -- an issued encoding is never redrawn or handed to another member -- so a new
    batch is drawn around them and a replay of an earlier batch needs their flanks.
    Measurements (``l1_swap_risk.tsv``) are over active members only.
    """
    proteins = ld.designed_proteins(MANIFEST, library=ld.LIBRARY_L1)
    keys = set(proteins)
    cds = ld.designed_cds(MANIFEST, library=ld.LIBRARY_L1)
    paths = {r["candidate_id"]: r for r in csv.DictReader(TAG_PATHS.open(), delimiter="\t")}
    holders: dict[str, str] = {}
    released: set[str] = set()
    issued = issued_batches()
    for r in csv.DictReader(MANIFEST.open(), delimiter="\t"):
        did = ld.design_id(r)
        if r["library"] == ld.LIBRARY_L1 and ld.is_active(r) and did in RELEASED:
            assert r["candidate_id"] in issued, f"{did}: released but holds no encoding"
            assert paths[r["candidate_id"]]["placement"] != "C-terminal", f"{did} is C-terminal"
            released.add(r["candidate_id"])
            if with_retired_holders:
                from deliverome_analysis.s3 import processed

                with open(processed.download(RELEASED[did])) as fh:
                    holders[r["candidate_id"]] = json.load(fh)[r["key"]]["cds"]
    # Snorkel members dropped from the built list (l1_snorkel.reserved_designs) keep their
    # issued encoding; it was drawn on the snorkel flank, so it is replayed on that flank.
    reserved_snorkels: set[str] = set()
    if with_retired_holders:
        # from the full manifest: a reserved snorkel's member may since be retired (CLRN3's
        # isoform, 2026-11-22 redraw) and so absent from l1_member_renumbering.tsv
        design_to_cid = {
            ld.design_id(r): r["candidate_id"]
            for r in csv.DictReader(MANIFEST.open(), delimiter="\t")
        }
        for did in l1_snorkel.reserved_designs():
            cid = design_to_cid[did]
            assert cid in issued, f"{did}: reserved snorkel holds no encoding"
            if cid in cds:
                holders[cid] = cds[cid]
            else:  # since retired (2026-11-22 isoform redraw): its CDS from the full record
                key = next(
                    r["key"]
                    for r in csv.DictReader(MANIFEST.open(), delimiter="\t")
                    if r["candidate_id"] == cid
                )
                holders[cid] = ld.load_designed_cds()[key]["cds"]
            reserved_snorkels.add(cid)
    if with_retired_holders:
        record = ld.load_designed_cds()
        for r in csv.DictReader(MANIFEST.open(), delimiter="\t"):
            if (
                r["library"] == ld.LIBRARY_L1
                and not ld.is_active(r)
                and r["candidate_id"] in issued
            ):
                holders[r["candidate_id"]] = record[r["key"]]["cds"]
    order = [
        r["candidate_id"]
        for r in csv.DictReader(MANIFEST.open(), delimiter="\t")
        if r["library"] == ld.LIBRARY_L1
    ]
    # Snorkel members (l1_snorkel.py) end on the snorkel module and then the cassette, so like
    # the C-terminal rung their cassette sits right before the barcode and needs an issued,
    # unique encoding. Keyed by design_id in the design table, mapped to candidate_id here.
    design_of = {
        r["candidate_id"]: r["design_id"]
        for r in csv.DictReader((HERE / "l1_member_renumbering.tsv").open(), delimiter="\t")
    }
    snorkels = set(l1_snorkel.built_designs())
    # Every arm ever issued a snorkel module (l1_snorkel_module_encodings.tsv) is a snorkel,
    # including ones since retired, whose rows leave the arm decisions.
    with (HERE / "l1_snorkel_module_encodings.tsv").open() as fh:
        snorkels |= {r["candidate_id"] for r in csv.DictReader(fh, delimiter="\t")}
    out = []
    for cid in order:
        if cid in reserved_snorkels:
            gene = paths[cid]["gene_symbol"] if cid in paths else ""  # retired: off the tag paths
            seq, cset, placement = holders[cid], "reserved", "snorkel"
        elif cid in holders:
            # A retired member keeps the flank its encoding was issued on: a snorkel member's
            # (retired by the 2026-11-22 isoform redraw, e.g. SLC22A24 l1-p3149) is the snorkel's.
            was_snorkel = design_of.get(cid) in snorkels
            seq, gene, cset = holders[cid], "", "retired"
            placement = "snorkel" if was_snorkel else "C-terminal"
        elif cid in keys and cid in paths:
            r = paths[cid]
            seq, gene, cset, placement = (
                cds[cid],
                r["gene_symbol"],
                r["candidate_set"],
                r["placement"],
            )
        else:
            continue
        body = seq[:-3] if seqtools.translate(seq[-3:]) == "*" else seq
        edit = ""
        snorkel = design_of.get(cid) in snorkels or cid in reserved_snorkels
        if snorkel:
            # The flank its cassette batches were issued against. Built snorkel arms carry a
            # per-arm module since 2026-11-21 (l1_snorkel_module_encodings.tsv); the cassette sits
            # after it either way, and replaying the batches needs the issuance-time flank.
            body = body + l1_snorkel.PREFIX_DNA
        if placement == "C-terminal":
            body, edit = domesticate_junction(body)
        out.append(
            {
                "id": cid,
                "gene": gene,
                "set": cset,
                "placement": placement,
                # True when the cassette is the last coding sequence before the barcode.
                "cterm": placement == "C-terminal" or snorkel,
                "body": body,
                "junction_edit": edit,
                "retired": cid in holders,
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
            "members_at_floor_tail": int((worst == len(TAIL)).sum()),
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
    members, _ = load_members(with_retired_holders=True)
    tag_fixed = qc.p003_tag_codons()
    batches = issued_batches()
    cterm = {m["id"] for m in members if m["cterm"]}
    missing = sorted(cterm - set(batches))
    stale = sorted(set(batches) - cterm)
    assert not stale, (
        f"issued encodings held by neither a C-terminal nor a retired member: {stale[:5]}"
    )
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
        for cid in sorted(upto - set(pins)):
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

    everyone, _ = load_members(with_retired_holders=True)
    # Measured: every active member at its current placement (a released holder, such as
    # TSPAN10, is measured as the internal member it now is, not by its old flank).
    members, unplaced = load_members()
    if unplaced:
        print(f"{len(unplaced)} active member(s) have no tag path and are not measured: {unplaced}")
    tag_fixed = qc.p003_tag_codons()
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
    # Drawn around every issued encoding, a retired member's included, and pinned: only
    # members without one get a new encoding, in a new batch.
    linkers, note = draw(
        everyone, fixed=tag_fixed, identity=LINKER_ONLY_JUNCTION_IDENTITY, pin=True
    )
    report["linkers_only_draw"] = note
    if not linkers:
        raise SystemExit(f"linkers-only draw failed: {note}")
    issued_all = {k: v.sequence for k, v in linkers.variants.items()}
    scenarios["linkers_only"] = {m["id"]: issued_all[m["id"]] for m in members if m["cterm"]}
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
            "members_at_floor_tail",
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
        # design_id is the permanent key (member ids are regenerated by renumbering); every
        # reader joins on it, and candidate_id is kept as the plate-order label.
        did_of = ld.design_ids(MANIFEST)
        w.writerow(
            [
                "candidate_id",
                "design_id",
                "scheme",
                "cassette_aa",
                "cassette_nt",
                "gc",
                "issued_batch",
            ]
        )
        for cid in sorted(issued_all):
            dna = issued_all[cid]
            assert seqtools.translate(dna) == CASSETTE_AA
            assert all(dna[3 * i : 3 * i + 3] == c for i, c in tag_fixed.items())
            batch = batches[cid][1] if cid in batches else new_batch
            w.writerow(
                [
                    cid,
                    did_of[cid],
                    "linkers_only",
                    CASSETTE_AA,
                    dna,
                    f"{seqtools.gc_content(dna):.3f}",
                    batch,
                ]
            )
    OUT_JSON.write_text(json.dumps(report, indent=2, default=str) + "\n")
    print(f"\nwrote {OUT_TSV.name}, {ENCODINGS.name}, {OUT_JSON.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
