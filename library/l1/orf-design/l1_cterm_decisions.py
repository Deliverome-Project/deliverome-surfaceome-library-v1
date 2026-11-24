"""The C-terminal rung's protein-level calls, one reversible row per decision.

Becca delegated these on 2026-11-12: "try to not queue protein specific things for
management decisions or anything that would block you from fully completing things. use
your best judgment to complete a new and improved evidence based l1 completely by the
morning." Rows from that work read ``decided_by = l1-cterm-build (delegated by Becca
2026-11-12, "use your best judgment")`` and name the alternative they rejected, so any row
can be overturned by flipping it. The ITM2B site is the exception: Becca chose it herself
on 2026-11-10. The rows of 2026-11-14 (KLRF1, the PLVAP isoform and the C-out calls) were decided by the coordinator and say so.

No row is typed per member. Which members move, stay or are held comes from
``build_l1_cterm_arms.dispositions()`` (``STALK_RULE``, ``COUT_RULE``, ``RULE_EXCEPTIONS``);
this script writes the row that explains each. Evidence comes from ``l1_cterm_review.tsv``
(UniProt processing and topology, DeepTMHMM v1 and v2, terminal pLDDT, terminal burial,
stalk) and the two refold-verdict tables, and is copied into each row. Writes
``l1_cterm_decisions.tsv``.

    uv run --locked python l1_cterm_decisions.py
    uv run --locked python l1_cterm_decisions.py --check    # compare with the record
"""

from __future__ import annotations

import csv
import io
import sys
from pathlib import Path

HERE = Path(__file__).parent
REVIEW = HERE / "l1_cterm_review.tsv"
VERDICTS = HERE / "l1_cterm_refold_verdicts.tsv"
OUT = HERE / "l1_cterm_decisions.tsv"

IMPLEMENTED_FLAG = "Claude Code session l1-review-export-scope (subagent), 2026-11-22"
DELEGATED = 'l1-cterm-build (delegated by Becca 2026-11-12, "use your best judgment")'
IMPLEMENTED = "Claude Code session l1-cterm-build (claude/l1-cterm-build-16eeeb), PR #284"
#: The stalk-site refold dispositions (PR #301 verdicts) came later, on their own branch.
IMPLEMENTED_STALK = (
    "Claude Code session l1-cterm-build (claude/l1-cterm-stalk-refold-holds), 2026-11-13"
)
IMPLEMENTED_AUDIT = (
    "Claude Code session l1-cterm-build (claude/l1-cterm-audit-followups), 2026-11-14"
)

#: The five genes recorded by hand on 2026-11-12 (PR #284); their rows keep that attribution.
#: Which genes get a row is derived (``build_l1_cterm_arms.released_genes``), not this list.
FURIN_ROWS_2026_11_12 = ("COL23A1", "COL25A1", "EDA", "TNFSF12", "TNFSF13B")

RELEASE_REASON = (
    "The processing is native, so the arm keeps native sequence: a mutated site would put a "
    "construct the cell does not make on the surface. A tag before a releasing cut would report "
    "the membrane stub, not the protein. At the C-terminus the tag reports the uncut fraction "
    "still on the surface; any extracellular tag on these members shares that exposure. ITM2B "
    "differs: furin removes only a 23-aa C-terminal peptide and the folded domain stays, so its "
    "tag moves before the site"
)
FOLD_REASON = (
    "No refold verdict yet. The N-terminus is cytoplasmic, so the only alternative is an "
    "internal stalk site, and the most buried termini (M13 family) have no stalk long enough "
    "to take the cassette. The cassette carries its own GGGGS linker"
)
TAAC_REASON = (
    "TAA is the strongest stop codon either way. CC is the pad of every construct designed and "
    "ordered since the 2026-10-01 handle revision, so keeping it keeps L1 comparable with them. "
    "The pad is shared by every L1 member, so changing it would change every insert on every "
    "rung for a +4 purine effect on termination that we have not measured. The readthrough "
    "product would be a short extension through the barcode handle"
)


def scope_of(r: dict) -> str:
    return f"{r['candidate_id']} ({r['design_id']}) {r['gene_symbol']} {r['species']}"


def evidence_of(r: dict) -> str:
    return (
        f"terminal pLDDT {r['terminal_plddt'] or 'n/a'}; terminal CA neighbors "
        f"{r['terminal_ca_neighbors_10a'] or 'n/a'}; stalk {r['stalk_residues'] or 'n/a'}"
    )


def cout_evidence(r: dict) -> str:
    return (
        f"DeepTMHMM v1: {r['n_tm_v1']} TM, last helix {r['last_tm_aa_v1']} residues, "
        f"extracellular tail {r['cterm_ecto_aa_v1']}; DeepTMHMM v2: "
        f"{r['v2_structural_type'] or 'no row'}, C-terminus {r['v2_cterm_side']}; UniProt "
        f"{r['uniprot_acc'] or 'no entry'}: {r['uniprot_cterm_topology'] or 'n/a'}"
        + (
            f"; UniProt signal peptide {r['uniprot_signal_peptide'] or 'none'}, "
            f"{r['uniprot_tm_helices']} TM"
            if r["uniprot_acc"]
            else ""
        )
    )


def structured_outcome(d) -> str:
    """What the rule made of a structured-terminus member's C-terminal tag (d: Disposition)."""
    if d.clause == "refold_pass":
        return "keep the cassette at the C-terminus: its refold passes"
    if d.clause == "stalk_refold_fail":
        return (
            f"superseded: held (not built), as the C-terminal and stalk-site (residue "
            f"{d.site}) refolds both fail; see its 'tagged refold failed' row"
        )
    if d.clause in ("stalk_refold_fail_flagged", "fold_fail_no_stalk_flagged"):
        return (
            "keep the cassette at the C-terminus, flagged: no site passes its refold and refold "
            "failures no longer hold (Becca, 2026-11-22); see its 'tagged refold failed' row"
        )
    if d.clause == "stalk_refold_pass":
        return (
            f"superseded: cassette after stalk residue {d.site}, whose refold "
            "passes; see its 'tagged refold failed' row"
        )
    if d.clause == "stalk_refold_pending":
        return f"superseded: stalk residue {d.site}, pending that design's refold"
    if d.clause in ("exception", "fold_fail_no_stalk"):
        return "superseded: held (not built); see its 'tagged refold failed' row"
    if d.action == "hold":
        return "superseded: held (not built); see its C-out call row"
    return (
        "keep the cassette at the C-terminus: the refold fails but the fold stays intact; see "
        "its 'tagged refold failed' row"
    )


def refold_row(r: dict, v: dict, sv: dict | None, d, arms) -> dict:
    """The 'tagged refold failed' row: what STALK_RULE (or its exception) does, and why."""
    by, impl, reason = DELEGATED, IMPLEMENTED, arms.STALK_RULE
    why = (
        "fold disturbed"
        if "lDDT" in v["verdict_reasons"]
        else "fold intact; a helix forms inside the cassette linker"
        if "helix inside the cassette" in v["verdict_reasons"]
        else "fold intact; epitope partly packed"
    )
    if d.clause == "stalk_refold_fail":
        decision = (
            f"hold (not built): the stalk-site design (residue {d.site}) "
            f"also fails its refold ({sv['verdict_reasons']}; site lDDT "
            f"{sv['lddt_site'] or 'n/a'}; PR #301)"
        )
        rejected = "keep the stalk-site arm flagged; keep the C-terminal tag"
        change = "a passing extracellular site from the ranker"
        impl = IMPLEMENTED_STALK
    elif d.clause == "stalk_refold_pass":
        decision = (
            f"move the cassette to stalk residue {d.site} (its refold "
            f"passes: native-region lDDT {sv['lddt_confident']}; PR #301)"
        )
        rejected = "keep the C-terminal tag; hold the member"
        change = "a failed refold of the stalk-site design: then hold"
        impl = IMPLEMENTED_STALK
    elif d.clause == "stalk_refold_pending":
        decision = f"move the cassette to stalk residue {d.site}, pending that design's refold"
        rejected = "keep the C-terminal tag; hold the member"
        change = "a failed refold of the stalk-site design: then hold"
    elif d.clause == "exception":
        decision = (
            "hold (not built): a recorded exception to the stalk rule, which alone gives "
            f"'{d.rule_action}'; left as on main"
        )
        rejected = "build the C-terminal arm, as the uniform rule gives"
        change = "the coordinator's decision on this exception; a native model and a passing site"
        by, impl = arms.COORDINATOR, arms.IMPLEMENTED_FOLLOWUP
        reason = arms.RULE_EXCEPTIONS[r["design_id"]][2]
    elif d.clause == "released_ectodomain":
        decision = f"keep the C-terminal tag ({why})"
        rejected = "stalk site: it lies before the ectodomain cut and would tag the stub"
        change = "a site after the furin cut that passes its refold"
        if "uniprot_not_checked" in r["flags"].split(";"):
            # An ortholog with no UniProt entry: say what the call rests on.
            decision += (
                "; the release call is transferred from the human gene, as this ortholog has "
                "no UniProt entry"
                + (
                    f", and its stalk site (residue {r['stalk_site_after']}) was never refolded"
                    if r["stalk_site_after"]
                    else ""
                )
            )
            change += "; this ortholog's own processing record; a refold of its stalk site"
    elif d.clause == "decided_site":
        decision = f"keep the cassette after residue {d.site} (Becca's site)"
        rejected = f"stalk residue {r['stalk_site_after']}"
        change = "Becca revisiting the site"
    elif d.clause in ("stalk_refold_fail_flagged", "fold_fail_no_stalk_flagged"):
        decision = "keep the cassette at the C-terminus, flagged (refold_flag=fail): " + (
            f"the stalk-site design (residue {r['stalk_site_after']}) also fails "
            f"({sv['verdict_reasons']})"
            if d.clause == "stalk_refold_fail_flagged" and sv
            else "the fold fails and the member has no usable stalk site"
        )
        rejected = "hold the member (the rule before 2026-11-22)"
        change = "a passing site (stalk or loop) from a later refold"
        by, impl = arms.FLAG_RULE_BY, IMPLEMENTED_FLAG
    elif d.clause == "fold_fail_no_stalk":
        decision = "hold (not built): the fold fails and the member has no usable stalk site"
        rejected = "keep the C-terminal tag, as for a packed epitope with no stalk site"
        change = "a native model and a passing site; a passing alternative site"
        by, impl = arms.COORDINATOR, arms.IMPLEMENTED_FOLLOWUP
    else:
        assert d.clause == "no_stalk_site", (r["design_id"], d)
        decision = f"keep the C-terminal tag, flagged ({why})"
        rejected = "none available: no stalk of 8+ residues"
        change = "a passing internal loop site from the ranker"
    return {
        "scope": scope_of(r),
        "question": "tagged refold failed",
        "decision": decision,
        "rejected_alternatives": rejected,
        "evidence": f"{v['verdict_reasons']}; native-region lDDT "
        f"{v['lddt_confident']}; site lDDT {v['lddt_site']}; " + evidence_of(r),
        "decided_by": by,
        "implemented_by": impl,
        "reason": reason,
        "would_change_it": change,
    }


def packing_row(r: dict, v: dict, sv: dict, d, arms) -> dict:
    """A member moved to its stalk site although its fold is intact (audit Q5: KLRF1)."""
    return {
        "scope": scope_of(r),
        "question": "moved to the stalk site on epitope packing alone",
        "decision": f"keep the cassette after stalk residue {d.site}",
        "rejected_alternatives": "return the cassette to the C-terminus, where a member with a "
        "packed epitope and no stalk site stays",
        "evidence": f"C-terminal refold: {v['verdict_reasons']}; native-region lDDT "
        f"{v['lddt_confident']}. Stalk-site refold: {sv['verdict']}, native-region lDDT "
        f"{sv['lddt_confident']}, {sv['epitope_confident_contact_res']} epitope residues packed",
        "decided_by": arms.COORDINATOR,
        "implemented_by": arms.IMPLEMENTED_FOLLOWUP,
        "reason": "The stalk-site arm passes its refold and the C-terminal form packs epitope "
        "residues on the domain, so the stalk site is the arm more likely to be detected. The "
        "rule is the same for every member: " + arms.STALK_RULE,
        "would_change_it": "evidence that the C-terminal epitopes are accessible, or that the "
        "stalk insertion disturbs the protein",
    }


def cout_rows(r: dict, d, cout, arms) -> list[dict]:
    """The row on whether the C-terminus is outside the cell (COUT_RULE), for one member."""
    flags = set(r["flags"].split(";"))
    weak = [f for f in arms.WEAK_COUT_FLAGS if f in flags]
    override = "lane_override_topology" in flags
    flagged = cout.clause == "predicted_cytoplasmic_kept"
    if not (weak or override or flagged or cout.action == arms.HOLD):
        return []
    note = arms.COUT_NOTES.get(r["design_id"], "")
    if weak or override:
        question = "weak C-out call (" + ", ".join(weak or ["lane_override_topology"]) + ")"
    elif cout.clause == "ortholog_of_held_human":
        question = "C-out call: the gene's human member is not built on UniProt's topology"
    else:
        question = (
            "C-out call: UniProt annotates the C-terminus cytoplasmic; DeepTMHMM v1 and v2 "
            "both call it extracellular"
        )
    if cout.action == arms.HOLD:
        decision = "hold (not built): " + (note.split(". ")[0] if note else cout.reason)
        rejected = "build the C-terminal arm on the DeepTMHMM call alone"
        change = (
            "experimental evidence that this form's C-terminus is extracellular; or a passing "
            "site on the rung its topology points to, which would move it there"
        )
    elif flagged:
        decision = (
            "keep the cassette at the C-terminus (kept; DeepTMHMM v1 and v2 agree C-out; "
            "UniProt's cytoplasmic annotation is a prediction; flagged for Becca)"
        )
        rejected = "take the member out of the build on a predicted annotation"
        change = "experimental evidence that the C-terminus is cytoplasmic; Becca's review"
    elif d.action == arms.HOLD:
        decision = (
            "C-out call stands on DeepTMHMM v1; the member has no arm for another reason "
            "(see its 'tagged refold failed' row)"
        )
        rejected = "drop the member for orientation as well"
        change = "an independent source placing the C-terminus in the cytoplasm"
    else:
        decision = "keep the cassette at the C-terminus"
        rejected = "take the member out of the build until its orientation is measured"
        change = "an independent source placing the C-terminus in the cytoplasm"
    return [
        {
            "scope": scope_of(r),
            "question": question,
            "decision": decision,
            "rejected_alternatives": rejected,
            "evidence": cout_evidence(r),
            "decided_by": arms.COORDINATOR,
            "implemented_by": arms.IMPLEMENTED_FOLLOWUP,
            "reason": ". ".join(
                x
                for x in (
                    cout.reason,
                    note,
                    arms.COUT_REROUTE.get(r["design_id"], "") + " (a note, not a build)"
                    if r["design_id"] in arms.COUT_REROUTE
                    else "",
                    arms.COUT_RULE,
                )
                if x
            ),
            "would_change_it": change,
        }
    ]


def rows() -> list[dict]:
    import build_l1_cterm_arms as arms

    with REVIEW.open() as fh:
        review = list(csv.DictReader(fh, delimiter="\t"))
    with VERDICTS.open() as fh:
        verdicts = {r["design_id"]: r for r in csv.DictReader(fh, delimiter="\t")}
    with arms.STALK_REFOLD_VERDICTS.open() as fh:
        stalk = {r["design_id"]: r for r in csv.DictReader(fh, delimiter="\t")}
    released = arms.released_genes(review)
    disposition = arms.dispositions(review, verdicts, stalk)
    couts = arms.cout_dispositions(review)
    out = []
    for r in review:
        design = r["design_id"]
        d = disposition[design]
        cout = couts[design]
        if r["gene_symbol"] in released and design not in arms.RELOCATED:
            furin = "furin" in (r["evidence"] or "").lower()
            out.append(
                {
                    "scope": scope_of(r),
                    "question": "UniProt records processing that cuts or releases the ectodomain"
                    + (" (furin)" if furin else ""),
                    "decision": "keep the cassette at the C-terminus, native sequence",
                    "rejected_alternatives": "mutate the cleavage site; tag before the site",
                    "evidence": (r["evidence"] or "UniProt annotation on the human gene")
                    + f"; {evidence_of(r)}",
                    "decided_by": DELEGATED,
                    "implemented_by": IMPLEMENTED
                    if r["gene_symbol"] in FURIN_ROWS_2026_11_12
                    else IMPLEMENTED_AUDIT,
                    "reason": RELEASE_REASON,
                    "would_change_it": "evidence that the surface fraction is undetectable",
                }
            )
        # a member added after the last refold batch has no verdict yet; its refold row comes
        # with the verdict, while its C-out rows are written now
        v = verdicts.get(design)
        if v and v["verdict"] == "fail" and cout.action != arms.HOLD:
            out.append(refold_row(r, v, stalk.get(design), d, arms))
            if d.clause == "stalk_refold_pass" and "lDDT" not in v["verdict_reasons"]:
                out.append(packing_row(r, v, stalk[design], d, arms))
        out += cout_rows(r, d, cout, arms)
        if "structured_terminus" in r["flags"].split(";") and r["cterm_tag"] == "kept":
            site = r["stalk_site_after"]
            out.append(
                {
                    "scope": scope_of(r),
                    "question": "terminus ends inside a confidently predicted fold",
                    "decision": structured_outcome(d),
                    "rejected_alternatives": (
                        f"cassette after stalk residue {site}"
                        if site
                        else "none available: no stalk of 8+ residues, N-terminus cytoplasmic"
                    ),
                    "evidence": evidence_of(r),
                    "decided_by": DELEGATED,
                    "implemented_by": IMPLEMENTED,
                    "reason": FOLD_REASON,
                    "would_change_it": (
                        f"a refold that disrupts the native fold: then move to stalk residue {site}"
                        if site
                        else "a refold that disrupts the native fold: then hold the member"
                    ),
                }
            )
    out.append(
        {
            "scope": "every L1 insert (shared 2-nt pad after the stop)",
            "question": "stop context TAAC (pad CC) or TAAG (pad GC)",
            "decision": "keep TAAC (pad CC)",
            "rejected_alternatives": "TAAG (pad GC)",
            "evidence": "build_constructs.py: pad CC since the 2026-10-01 handle revision",
            "decided_by": DELEGATED,
            "implemented_by": IMPLEMENTED,
            "reason": TAAC_REASON,
            "would_change_it": "a readthrough measurement on a tagged construct",
        }
    )
    return out


def main(argv=None) -> int:
    out = rows()
    buf = io.StringIO(newline="")
    w = csv.DictWriter(buf, fieldnames=list(out[0]), delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(out)
    if "--check" in (sys.argv[1:] if argv is None else argv):
        assert OUT.read_text() == buf.getvalue(), f"{OUT.name} does not reproduce"
        print(f"{OUT.name} reproduces ({len(out)} decisions)")
        return 0
    OUT.write_text(buf.getvalue())
    print(f"wrote {len(out)} decisions to {OUT.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
