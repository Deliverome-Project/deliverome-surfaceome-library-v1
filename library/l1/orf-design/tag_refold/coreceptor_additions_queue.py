"""Queue the refolds of the 25 co-receptor additions (Becca, 2026-11-21: "Yes, rebuild now").

The full queue builder (``build_refold_inputs.py``) re-reads every rung's snapshotted export, so
this appends only the new constructs, built with the same ``construct()`` (residue map, cassette
span, SHA-256s), with the same fold settings and verdict tests downstream. Every existing row of
the queue, natives, tagged and sources tables is kept byte for byte (asserted).

Constructs, one fold per tagged SHA-256, each from its rung's own record:

* N after SP, 22 arms: ``reviewer-package/nterm/refold_requests_coreceptor_additions.tsv``; the
  tagged precursor is the built insert translated (``constructs_l1_nterm.fasta``), and must hash
  to the request's ``tagged_sha256``. An Igk-swap arm is folded as its precursor, as before.
* internal, SCNN1B: s5e13 ``refold_requests_coreceptor_additions.tsv`` (``tagged_protein``).
* residue-1 multipass, SLC51A: the Igk mature chain (``multipass_arms_mature.fasta``).
* C-terminal, KLRK1: the designed protein plus the cassette at its C-terminus.
* natives with no model of exactly their sequence (AFDB v6 differs for GRIN1, HLA-DQA1 and
  HLA-DQB1): an untagged fold of the designed protein, ``native.<sha12>``.

Writes ``tag_refold_queue_coreceptor_additions.tsv`` (the fold script's queue, sequences inline)
and appends to ``tag_refold_queue.tsv``, ``tag_refold_natives.tsv``, ``tag_refold_tagged.tsv``
and ``tag_refold_sources.tsv``. Implemented by Claude Code session l1-review-export-scope
(subagent), 2026-11-21.

    uv run --locked python library/l1/orf-design/tag_refold/coreceptor_additions_queue.py
    uv run --locked python library/l1/orf-design/tag_refold/coreceptor_additions_queue.py --round2
    uv run --locked python library/l1/orf-design/tag_refold/coreceptor_additions_queue.py \
        --natives l1-p3210,... --out tag_refold_queue_natives_2026_11_22.tsv
    uv run --locked python library/l1/orf-design/tag_refold/coreceptor_additions_queue.py \
        --all-rungs --out tag_refold_queue_2026_11_22_all_rungs.tsv

``--all-rungs`` (2026-11-22): one batch for every built or requested arm, on every rung, whose
tagged protein has never been folded (not in ``tag_refold_queue.tsv``): each rung's request
files (N after SP ``refold_requests.tsv`` and ``refold_requests_alternatives.tsv``; s5e13's
``refold_requests_*.tsv`` for the internal rung) plus the rung's built arms (``constructs_l1_*``:
the N-terminal and C-terminal precursors, the type III and residue-1 multipass mature chains with
the Igk leader removed, as ``build_refold_inputs.py`` folds them). Natives without a model of
their exact sequence are queued too. Membership is by sequence, so it covers the round-2
additions, the redraw and refills, the MANE switches, re-sited arms and new flag-not-hold arms
alike, and an arm already folded is never folded twice.

``--all-rungs --members ID,...`` (2026-11-22 cut fill) sweeps only the named members' requests and
built arms, by the same rules. It queued the four C-terminal arms of the 2026-11-22-rebuild order
of record whose only fold of the exact protein was superseded (CD70, BRICD5, CD209, SCARA5: built
again under Becca's 2026-11-22 flag-not-hold rule), so ``fold_order_check.py`` finds a live fold
for every ordered insert:

    uv run --locked python library/l1/orf-design/tag_refold/coreceptor_additions_queue.py \
        --all-rungs --members l1-p0314,l1-p0184,l1-p0274,l1-p1763 \
        --out tag_refold_queue_2026_11_22_cut_fill.tsv

``--native-switch ID,... --out NAME`` (2026-11-22): a member whose native protein changed (a MANE
switch) while its built tagged protein did not keeps a live fold made on the old native's MSA and
scored against the old native's model (MUC15, l1-p1266: the Igk leader replaces the residues the
switch added). Each such tagged protein is folded again on the current native, under the fold id
``t.<tagged sha16>.n<native sha12>`` so the earlier fold's rows stay as they are:

    uv run --locked python library/l1/orf-design/tag_refold/coreceptor_additions_queue.py \
        --native-switch l1-p1266 --out tag_refold_queue_2026_11_22_native_switch.tsv
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
EXP = HERE.parent
REPO = EXP.parents[2]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(EXP))

from build_refold_inputs import CASSETTE, construct, fasta, sha, translate_insert  # noqa: E402

from deliverome_analysis import library_design as ld  # noqa: E402

QUEUE = EXP / "tag_refold_queue.tsv"
NATIVES = EXP / "tag_refold_natives.tsv"
TAGGED = EXP / "tag_refold_tagged.tsv"
SOURCES = EXP / "tag_refold_sources.tsv"
NATIVE_INDEX = EXP / "tag_refold_native_models.tsv"
OUT = EXP / "tag_refold_queue_coreceptor_additions.tsv"
NTERM_REQ = EXP / "reviewer-package/nterm/refold_requests_coreceptor_additions.tsv"
LOOP_REQ = (
    REPO / "library/l1/tag-site-benchmarking/refold_requests_coreceptor_additions.tsv"
)
DECIDED = 'l1-review-export-scope (subagent), under Becca\'s 2026-11-21 "Yes, rebuild now"'
IMPLEMENTED = "Claude Code session l1-review-export-scope (subagent), 2026-11-21"
REASON = "co-receptor addition (2026-11-21): the refold can change the call"
REASON_ALL_RUNGS = (
    "2026-11-22 addition, MANE switch or re-sited arm: the refold can change the call"
)


def tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open() as fh:
        r = csv.DictReader(fh, delimiter="\t")
        return list(r.fieldnames or []), list(r)


def constructs() -> list[dict]:
    proteins = ld.designed_proteins(EXP / "candidate_manifest_tranches.tsv", library=ld.LIBRARY_L1)
    gene = {c: p["gene_symbol"] for c, p in proteins.items()}
    out = []
    inserts = {h.split("|")[0]: s for h, s in fasta(EXP / "constructs_l1_nterm.fasta")}
    for r in tsv(NTERM_REQ)[1]:
        cid = r["candidate_id"]
        tagged = translate_insert(inserts[cid])
        assert sha(tagged) == r["tagged_sha256"], (
            f"{cid}: built insert is not the requested protein"
        )
        out.append(construct(
            "nterm", cid, r["design_id"], gene[cid], proteins[cid]["sequence"], tagged,
            source=f"{NTERM_REQ.name} (nterm rung, 2026-11-21), arm {r['arm']}",
            tier=1, order=0, reason=REASON, extra={"arm": r["arm"]},
        ))  # fmt: skip
    for r in tsv(LOOP_REQ)[1]:
        d = r["design_id"]
        out.append(construct(
            "loop", r["arm_id"], d, gene[d], proteins[d]["sequence"], r["tagged_protein"],
            source=f"{LOOP_REQ.name} (s5e13 internal rung, 2026-11-21)", tier=1, order=0,
            reason=REASON, extra={"rank": r["rank"], "loop_member": d},
            published_sha=r["tagged_protein_sha256"],
        ))  # fmt: skip
    for h, s in fasta(EXP / "multipass_arms_mature.fasta"):
        arm, d = h.split()[0], h.split()[1]
        if d != "l1-p3165":
            continue
        out.append(construct(
            "multipass", arm, d, gene[d], proteins[d]["sequence"], s,
            source="multipass_arms_mature.fasta (2026-11-21), mature chain", tier=1, order=0,
            reason=REASON, extra={"design": "igk_tag"},
        ))  # fmt: skip
    d = "l1-p3164"
    native = proteins[d]["sequence"]
    out.append(construct(
        "cterm", d, d, gene[d], native, native + CASSETTE,
        source="designed protein + cassette (C-terminal rung, 2026-11-21)", tier=1, order=0,
        reason=REASON,
    ))  # fmt: skip
    return out


def round2() -> list[dict]:
    """The nterm rung's replacement and alternative sites for members whose first arm failed.

    ``nterm_decide.py`` sends every remaining SignalP-passing alternative at once
    (``refold_requests.tsv`` and ``refold_requests_alternatives.tsv``). The requests name the
    tagged protein by SHA-256 only; each is rebuilt here (native site: native[:n] + cassette +
    native[n:]; Igk swap: leader + Asp + native[c:n] + cassette + native[n:], c searched) and
    must hash to it.
    """
    from build_constructs import CASSETTE_AA
    from build_gpcr_rescue_arms import LEADER_AA, PLUS1_AA

    proteins = ld.designed_proteins(EXP / "candidate_manifest_tranches.tsv", library=ld.LIBRARY_L1)
    members = {f"l1-p{n}" for n in range(3153, 3178)}
    out = []
    for name in ("refold_requests.tsv", "refold_requests_alternatives.tsv"):
        for r in tsv(EXP / "reviewer-package/nterm" / name)[1]:
            cid = r["candidate_id"]
            if cid not in members:
                continue
            nat, n = proteins[cid]["sequence"], int(r["insert_after_residue"])
            if r["arm"] == "native_site":
                cands = [nat[:n] + CASSETTE_AA + nat[n:]]
            else:
                # the site counts precursor residues (leader + Asp + native[c:])
                k = n - len(LEADER_AA) - len(PLUS1_AA)
                cands = [
                    LEADER_AA + PLUS1_AA + nat[c:][:k] + CASSETTE_AA + nat[c:][k:]
                    for c in range(5, 80)
                ]
            hit = [t for t in cands if sha(t) == r["tagged_sha256"]]
            assert hit, f"{cid} {r['arm']} {n}: no construct hashes to the request"
            out.append(construct(
                "nterm", cid, r["design_id"], proteins[cid]["gene_symbol"], nat, hit[0],
                source=f"{name} (nterm rung, 2026-11-21 round 2), arm {r['arm']}", tier=1,
                order=0, reason="co-receptor addition, round 2: replaces a failed arm",
                extra={"arm": r["arm"]},
            ))  # fmt: skip
    return out


def all_rungs(members: set[str] | None = None) -> list[dict]:
    """Every requested or built arm, on every rung, whose tagged protein was never folded.

    ``members`` replaces the swept set of members (default: the 2026-11-22 additions, redraw,
    refills, restorations, MANE switches and the two re-sited co-receptor members)."""
    from build_constructs import CASSETTE_AA
    from build_gpcr_rescue_arms import LEADER_AA, PLUS1_AA

    proteins = ld.designed_proteins(EXP / "candidate_manifest_tranches.tsv", library=ld.LIBRARY_L1)
    design = ld.design_ids(EXP / "candidate_manifest_tranches.tsv")
    # folded = queued, except a fold whose only metric row is superseded (no live verdict)
    superseded = {
        r["tagged_sha256"]
        for r in tsv(EXP / "tag_refold_metrics.tsv")[1]
        if r["verdict"] == "superseded"
    } - {r["tagged_sha256"] for r in tsv(EXP / "tag_refold_verdicts.tsv")[1]}
    folded = {r["tagged_sha256"] for r in tsv(QUEUE)[1]} - superseded
    out: dict[str, dict] = {}
    src = "2026-11-22 all-rungs batch" if members is None else "2026-11-22 cut fill"
    # Built arms are swept only for the members whose arms are new or changed on 2026-11-22: the
    # additions, redraw, refills and restorations (l1-p3178 on), the MANE switches, and the two
    # co-receptor members re-sited on the N-terminal rung (CHRNB4 at 28, HLA-DQA1 at 31).
    mane = {r["candidate_id"] for r in tsv(EXP / "mane_gtf_changes_2026-11-22.tsv")[1]}
    target = {c for c in proteins if c >= "l1-p3178"} | mane | {"l1-p3160", "l1-p3169"}
    if members is not None:
        assert members <= set(proteins), f"not L1 members: {sorted(members - set(proteins))}"
        target = set(members)

    def add(rung, source_id, cid, tagged, source, extra=None, published=None):
        h = sha(tagged)
        if h in folded or h in out:
            return
        nat = proteins[cid]["sequence"]
        out[h] = construct(
            rung, source_id, design.get(cid, cid), proteins[cid]["gene_symbol"], nat, tagged,
            source=source, tier=1, order=0,
            reason=REASON_ALL_RUNGS,
            extra=extra or {}, published_sha=published,
        )  # fmt: skip

    # N after SP: requests, then built arms
    for name in ("refold_requests.tsv", "refold_requests_alternatives.tsv"):
        for r in tsv(EXP / "reviewer-package/nterm" / name)[1]:
            cid = r["candidate_id"]
            if cid not in proteins or r["tagged_sha256"] in folded:
                continue
            nat, n = proteins[cid]["sequence"], int(r["insert_after_residue"])
            if r["arm"] == "native_site":
                cands = [nat[:n] + CASSETTE_AA + nat[n:]]
            else:
                k = n - len(LEADER_AA) - len(PLUS1_AA)
                cands = [
                    LEADER_AA + PLUS1_AA + nat[c:][:k] + CASSETTE_AA + nat[c:][k:]
                    for c in range(5, 80)
                ]
            hit = [t for t in cands if sha(t) == r["tagged_sha256"]]
            assert hit, f"{cid} {r['arm']} {n}: no construct hashes to the request"
            add("nterm", cid, cid, hit[0], f"{name} (nterm rung, {src}), arm {r['arm']}",
                {"arm": r["arm"]})  # fmt: skip
    for h, s in fasta(EXP / "constructs_l1_nterm.fasta"):
        cid = h.split("|")[0].split()[0]
        if cid in target:
            add("nterm", cid, cid, translate_insert(s), f"constructs_l1_nterm.fasta ({src})")
    # internal: every s5e13 request file, by its tagged protein
    s5 = REPO / "library/l1/tag-site-benchmarking"
    for path in sorted(s5.glob("refold_requests_*.tsv")):
        cols, rows = tsv(path)
        if "tagged_protein" not in cols:
            continue
        for r in rows:
            d = r["design_id"]
            if d in proteins and r["tagged_protein"]:
                if sha(proteins[d]["sequence"]) != r["sequence_sha256"]:
                    continue  # a request made on an earlier protein of the member
                add("loop", r.get("arm_id") or f"{d}:after{r['insert_after_residue']}", d,
                    r["tagged_protein"], f"{path.name} (s5e13 internal rung, {src})",
                    {"rank": r.get("rank", ""), "loop_member": d},
                    r.get("tagged_protein_sha256"))  # fmt: skip
    # residue-1 multipass and type III: the mature chain (Igk leader removed)
    for h, m in fasta(EXP / "multipass_arms_mature.fasta"):
        arm, d = h.split()[0], h.split()[1]
        if d in target:
            add("multipass", arm, d, m, f"multipass_arms_mature.fasta ({src}), mature chain")
    for h, s in fasta(EXP / "constructs_l1_type3.fasta"):
        cid = h.split("|")[0].split()[0].split(":")[0]
        if cid not in target:
            continue
        t = translate_insert(s)
        if t.startswith(LEADER_AA):
            t = t[len(LEADER_AA) :]
        add("type3", cid, cid, t, f"constructs_l1_type3.fasta ({src})")
    # C-terminal: the designed protein plus the cassette (a member held pending its refold has
    # no built arm yet), then every built arm and snorkel arm as its precursor
    paths = {r["candidate_id"]: r for r in tsv(EXP / "l1_tag_paths.tsv")[1]}
    for cid in sorted(target):
        if paths.get(cid, {}).get("placement") == "C-terminal":
            add("cterm", cid, cid, proteins[cid]["sequence"] + CASSETTE_AA,
                f"designed protein + cassette (C-terminal rung, {src})")  # fmt: skip
    for h, s in fasta(EXP / "constructs_l1_cterm.fasta"):
        arm = h.split("|")[0].split()[0]
        cid = ld.design_of_arm(arm) if ld.ARM_ID_RE.fullmatch(arm) else arm.split(":")[0]
        if cid not in target:
            continue
        add("cterm", arm, cid, translate_insert(s), f"constructs_l1_cterm.fasta ({src})")
    return list(out.values())


def native_switch(ids: list[str]) -> list[dict]:
    """Live tagged folds of these members' built arms made on a native that is no longer the
    member's, rebuilt on the current native (same tagged protein, new residue map and MSA).
    Unbuilt alternatives folded on the old native are left alone."""
    built = {
        r[c]
        for path in EXP.glob("constructs_l1_*.tsv")
        for r in tsv(path)[1]
        for c in ("protein_sha256", "mature_sha256", "refold_tagged_sha256")
        if r.get(c)
    }
    proteins = ld.designed_proteins(EXP / "candidate_manifest_tranches.tsv", library=ld.LIBRARY_L1)
    tagged = {r["tagged_sha256"]: r["tagged_seq"] for r in tsv(TAGGED)[1]}
    live = [r for r in tsv(EXP / "tag_refold_structures.tsv")[1] if r["kind"] == "tagged"]
    out: dict[str, dict] = {}
    for d in ids:
        nat = proteins[d]["sequence"]
        h = sha(nat)
        current = {
            r["tagged_sha256"] for r in live if r["design_id"] == d and r["native_sha256"] == h
        }
        for r in live:
            t = r["tagged_sha256"]
            if r["design_id"] != d or r["native_sha256"] == h or t in current or t in out:
                continue
            if t not in built:
                continue
            c = construct(
                r["rung"], d, d, proteins[d]["gene_symbol"], nat, tagged[t],
                source=f"{r['fold_id']} refolded on the current native ({h[:12]}), 2026-11-22",
                tier=1, order=0,
                reason=f"native switch: {r['fold_id']} was folded on {r['native_sha256'][:12]}",
            )  # fmt: skip
            out[t] = {**c, "fold_id": f"t.{t[:16]}.n{h[:12]}"}
    return list(out.values())


def native_folds(ids: list[str]) -> list[dict]:
    """Untagged folds of members' designed proteins (no identical AlphaFold DB model): the
    comparison models their tagged refolds and the s5e13 features need. One per sequence."""
    proteins = ld.designed_proteins(EXP / "candidate_manifest_tranches.tsv", library=ld.LIBRARY_L1)
    out = {}
    for d in ids:
        s = proteins[d]["sequence"]
        h = sha(s)
        out.setdefault(h, {
            "fold_id": f"native.{h[:12]}", "kind": "native", "tier": 1, "order": -1,
            "rung": "native", "design_id": d, "source_id": "",
            "gene_symbol": proteins[d]["gene_symbol"], "native_aa": len(s), "tagged_aa": len(s),
            "native_sha256": h, "tagged_sha256": h, "blocks": f"[[0, 0, {len(s)}]]",
            "native_seq": s, "tagged_seq": s,
            "source": "native of a member with no model of this exact sequence",
            "reason": "the comparison model for its tagged forms",
            "decided_by": DECIDED, "implemented_by": IMPLEMENTED,
        })  # fmt: skip
    return list(out.values())


def main() -> None:
    global OUT, DECIDED, IMPLEMENTED, REASON_ALL_RUNGS
    if "--natives" in sys.argv:
        # --natives ID,ID,... --out NAME: one batch of native folds (no rung constructs)
        ids = sys.argv[sys.argv.index("--natives") + 1].split(",")
        OUT = EXP / sys.argv[sys.argv.index("--out") + 1]
        new = native_folds(ids)
        qcols = tsv(QUEUE)[0]
        with OUT.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=[*qcols, "native_seq", "tagged_seq"],
                               delimiter="\t", lineterminator="\n", extrasaction="ignore")  # fmt: skip
            w.writeheader()
            w.writerows(new)
        print(f"{len(new)} native folds -> {OUT.name}")
        return
    if "--native-switch" in sys.argv:
        OUT = EXP / sys.argv[sys.argv.index("--out") + 1]
        DECIDED = (
            "coordinator, 2026-11-22: a tagged fold must be on the member's current native "
            "(l1_consistency); under Becca's 2026-11-22 rebuild"
        )
        IMPLEMENTED = "Claude Code session l1-review-export-scope (subagent), 2026-11-22"
        new = native_switch(sys.argv[sys.argv.index("--native-switch") + 1].split(","))
    elif "--all-rungs" in sys.argv:
        OUT = EXP / sys.argv[sys.argv.index("--out") + 1]
        members = None
        if "--members" in sys.argv:
            members = set(sys.argv[sys.argv.index("--members") + 1].split(","))
            DECIDED = (
                "Becca, 2026-11-22: every ordered insert needs a fold of its exact protein; "
                "flag-not-hold builds the refold-held C-terminal arms"
            )
            IMPLEMENTED = "Claude Code session l1-review-export-scope (subagent), 2026-11-22"
            REASON_ALL_RUNGS = (
                "2026-11-22 cut fill: an ordered insert whose only fold of its exact protein "
                "was superseded"
            )
        new = all_rungs(members)
    elif "--round2" in sys.argv:
        OUT = EXP / "tag_refold_queue_coreceptor_additions_r2.tsv"
        new = round2()
    else:
        new = constructs()
    have_model = {r["native_sha256"] for r in tsv(NATIVE_INDEX)[1]}
    qcols, qrows = tsv(QUEUE)
    folded = {r["tagged_sha256"] for r in qrows}
    rows, natives = [], {}
    folded = {r["fold_id"] for r in qrows}
    for c in new:
        fid = c.get("fold_id") or f"t.{c['tagged_sha256'][:16]}"
        assert fid not in folded, f"{c['source_id']}: already queued"
        natives[c["native_sha256"]] = c["native_seq"]
        rows.append({**c, "fold_id": fid, "kind": "tagged",
                     "decided_by": DECIDED, "implemented_by": IMPLEMENTED})  # fmt: skip
    for h, s in sorted(natives.items()):
        if h in have_model or h in {r["native_sha256"] for r in qrows if r["kind"] == "native"}:
            continue
        rows.append({
            "fold_id": f"native.{h[:12]}", "kind": "native", "tier": 1, "order": -1,
            "rung": "native", "design_id": next(c["design_id"] for c in new if c["native_sha256"] == h),
            "source_id": "", "gene_symbol": next(c["gene_symbol"] for c in new if c["native_sha256"] == h),
            "native_aa": len(s), "tagged_aa": len(s), "native_mapped_aa": "", "cassette_start": "",
            "cassette_end": "", "insert_after_native_residue": "", "leader_end": "", "sp_end": "",
            "tm1_span_folded": "", "native_sha256": h, "tagged_sha256": h,
            "blocks": f"[[0, 0, {len(s)}]]", "native_seq": s, "tagged_seq": s,
            "source": "native of a queued construct with no model of this exact sequence",
            "reason": "the comparison model for its tagged forms",
            "decided_by": DECIDED, "implemented_by": IMPLEMENTED,
        })  # fmt: skip
    # the fold script's own queue: every column plus the sequences inline
    with OUT.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=[*qcols, "native_seq", "tagged_seq"], delimiter="\t",
                           lineterminator="\n", extrasaction="ignore")  # fmt: skip
        w.writeheader()
        w.writerows(rows)

    def append(path: Path, new_rows: list[dict], key: str) -> int:
        cols, old = tsv(path)
        before = path.read_bytes()
        seen = {r[key] for r in old}
        add = [r for r in new_rows if r[key] not in seen]
        with path.open("a", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n",
                               extrasaction="ignore")  # fmt: skip
            w.writerows(add)
        assert path.read_bytes().startswith(before), f"{path.name}: an existing row changed"
        return len(add)

    n_q = append(QUEUE, rows, "fold_id")
    n_n = append(NATIVES, [{"native_sha256": h, "native_seq": s} for h, s in natives.items()],
                 "native_sha256")  # fmt: skip
    n_t = append(TAGGED, [{"tagged_sha256": r["tagged_sha256"], "tagged_seq": r["tagged_seq"]}
                          for r in rows], "tagged_sha256")  # fmt: skip
    src = [{**r, "tier": r["tier"]} for r in rows if r["kind"] == "tagged"]
    scols, sold = tsv(SOURCES)
    have = {(r["tagged_sha256"], r["source_id"]) for r in sold}
    before = SOURCES.read_bytes()
    with SOURCES.open("a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=scols, delimiter="\t", lineterminator="\n",
                           extrasaction="ignore")  # fmt: skip
        w.writerows(r for r in src if (r["tagged_sha256"], r["source_id"]) not in have)
    assert SOURCES.read_bytes().startswith(before)
    print(f"{len(rows)} folds ({sum(r['kind'] == 'native' for r in rows)} native) -> {OUT.name}; "
          f"appended queue {n_q}, natives {n_n}, tagged {n_t}, sources {len(src)}")  # fmt: skip


if __name__ == "__main__":
    main()
