"""Put the tag cassette on the C-terminus of every L1 C-terminal member, as synthesis inserts.

This is the C-terminal rung of issue #2 ("assemble the tag cassette into every L1 insert").
Every input it uses is already fixed in the record; nothing is chosen here:

* **which members**: ``placement == "C-terminal"`` in ``l1_tag_paths.tsv`` (DeepTMHMM v1,
  matched by sequence). A member whose C-terminus is cut off (``cterm_tag == "removed"`` in
  ``l1_cterm_review.tsv``) is built only at a decided internal site (``RELOCATED``: ITM2B,
  after residue 235); any other such member is held and listed in the summary.
* **where the cassette goes, and who is held**: ``dispositions()``, which applies two stated
  rules to committed tables and nothing typed per member: ``COUT_RULE`` (a weak C-out call
  that an independent source contradicts is held) and ``STALK_RULE`` (a failing C-terminal
  refold moves the cassette to the stalk site if that design passes its own refold, holds
  the member if it fails too, and otherwise keeps the C-terminus). A member left as it was
  where the rule says otherwise is a row in ``RULE_EXCEPTIONS`` with its reason.
* **the ORF**: the member's domesticated ORF from ``candidate_constructs.json`` (stage 10),
  which must translate to the member's designed protein (``candidate_cds.json``) plus a stop.
* **the cassette**: the member's own linkers-only encoding from ``l1_cassette_encodings.tsv``
  (stage 14). SpyTag003 and P-ALFA-P on p003's codons, GGGGS linkers recoded.
* **the layout**: p003's, from ``build_constructs.py``, in the Esp3I wrapper
  ``assemble_candidates.py`` puts every L1 insert in::

      CGTCTCA | Kozak | ORF (native stop removed) | cassette (51 aa) | TAA | CC
              | T7 | 5' handle | barcode (20 nt) | Feldman 3' | 10x CS1 | AGA GAGACG

Barcodes are not reserved yet (issue #1), so the barcode is a 20-N placeholder. Everything
upstream of it (the whole coding sequence and the swap-relevant stretch) is final.

Per member the build asserts that the protein is exactly ``native + cassette``; that the
coding sequence has one stop, at its end; that the insert has the two wrapper Esp3I sites and
no other BbsI/Esp3I/BsaI site, alone and inside the p0039 backbone; and that the tag codons
are p003's. One junction can need an edit: an ORF whose last bases complete a Type IIS site
with the cassette's first Gly codon gets a synonymous last-codon change
(``l1_swap_risk.domesticate_junction``; SLC41A2 is the known case). Across the library it then
measures, on the built sequences, the longest stretch each member shares with any other up to
its barcode (``l1_swap_risk.worst_partner``), and scans each insert for the vendor
complexity features ``synthesis_risk_scan.py`` reports, keeping only those that touch the
cassette or its junctions: the ORF's own features exist with or without the tag.

    uv run --locked python build_l1_cterm_arms.py            # build and write
    uv run --locked python build_l1_cterm_arms.py --check    # rebuild, compare with the record
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
import tempfile
from collections import Counter
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

from build_constructs import CASSETTE_AA, STOP
from diversify_snorkel_modules import module_dna
from l1_swap_risk import GATE, PREFERRED, domesticate_junction, worst_partner
from synthesis_risk_scan import (
    GC_WINDOW,
    MAX_LOOP,
    MIN_DIRECT_REPEAT,
    MIN_HOMOPOLYMER,
    MIN_STEM,
    direct_repeats,
    gc_windows,
    hairpins,
    homopolymers,
)

from deliverome_analysis.cloning import library_design as ld
from deliverome_analysis.cloning import plasmid_qc as qc
from deliverome_analysis.cloning import seqtools

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
AP = REPO / "library/l1/proteoform-selection"
sys.path.insert(0, str(AP))

import l1_snorkel  # noqa: E402

MANIFEST = HERE / "candidate_manifest_tranches.tsv"
TAG_PATHS = HERE / "l1_tag_paths.tsv"
MAPPING = HERE / "l1_member_renumbering.tsv"
REVIEW = HERE / "l1_cterm_review.tsv"
ENCODINGS = HERE / "l1_cassette_encodings.tsv"
OUT_TSV = HERE / "constructs_l1_cterm.tsv"
OUT_FASTA = HERE / "constructs_l1_cterm.fasta"
OUT_JSON = HERE / "constructs_l1_cterm_summary.json"
# The pinned domesticated-ORF record (extended by #283 for new members); never a hard-coded key.
CONSTRUCTS_KEY = ld.CONSTRUCTS_KEY
PUBLISH_PREFIX = "sprint4/s4e9-orf-design-fg-panel/library-design/cterm-arms/2026-11-12/"

BARCODE_SLOT = qc.BARCODE_SLOT
ENZYMES = list(qc.DOMESTICATION_ENZYMES)
TAIL = qc.STOP_TO_BARCODE  # constant between the cassette and the barcode

#: Held members that get the cassette at a decided internal position instead, keyed by
#: candidate_id: (insert after this residue, the native residues that must follow it).
#: ITM2B (Becca, 2026-11-10): furin cuts after R243 and releases the C-terminal Bri23
#: peptide, so a C-terminal tag is lost, and a tag placed straight after R243 would leave the
#: whole cut motif in front of it. The cassette goes after E235 instead: past the BRICHOS
#: domain (137-231, its 164-223 disulfide untouched) and the confidently predicted residues
#: (pLDDT >= 86 through 234), in the low-confidence linker (pLDDT 76 at 236, 62 at 240), and
#: far enough upstream that the eight residues before the cut (P8-P1, 236-243) stay
#: contiguous. The tag stays on the membrane form after furin; ADAM10 later sheds that
#: ectodomain, as for the other ectodomain-released members.
#:
#: Attribution (convention of 2026-11-10: every manual fix names who decided it and who
#: implemented it). Decided by Becca Carlson on 2026-11-10, choosing "Insert before furin site"
#: over mutating the furin site. Residue 235 was proposed by, and implemented in, the Claude
#: Code session on branch claude/l1-cterm-build-16eeeb (PR #284). That session also found that
#: the first suggestion, a tag straight after R243, would still be cut off by furin. The fix
#: regenerates from this constant; the inputs it read (UniProt features, the member's model)
#: are in the UniProt snapshot and library_structure_payload.tsv.
ITM2B_REASON = (
    "furin releases the C-terminal Bri23 peptide (244-266), taking a C-terminal tag with it; "
    "cassette after E235: past BRICHOS (137-231), in the linker, P8-P1 (236-243) kept intact"
)
DECIDED_BY = "Becca Carlson, 2026-11-10"
IMPLEMENTED_BY = "Claude Code session claude/l1-cterm-build-16eeeb, PR #284"
#: The junction edit (SLC41A2's last codon GAC->GAT) is computed, not hand-picked: the rule
#: recodes an ORF's last codon synonymously when it completes a Type IIS site with the cassette's
#: first Gly codon. The rule entered in commit e545db2e (2026-10-29, l1_swap_risk.py) and is
#: reapplied here to the domesticated ORF; its rows are attributed so the edit is traceable.
JUNCTION_DECIDED_BY = "rule (l1_swap_risk.domesticate_junction, commit e545db2e, 2026-10-29)"
JUNCTION_REASON = (
    "synonymous last-codon change: the native codon plus the cassette's GGN made Esp3I"
)


DELEGATED_BY = 'l1-cterm-build (delegated by Becca 2026-11-12, "use your best judgment")'
#: Decisions of 2026-11-14 (the audit follow-up): the uniform stalk rule and the C-out holds.
COORDINATOR = (
    "coordinator under Becca's 2026-11-12 delegation and her 2026-11-14 instruction to manage "
    "the C-terminal rung"
)
IMPLEMENTED_FOLLOWUP = "coordinator subagent, branch claude/l1-cterm-audit-followups-2, 2026-11-14"


def attribution(at_cterm: bool, junction_edit: str, refold: str = "") -> dict[str, str]:
    """decided_by / implemented_by / reason for a row; blank only when nothing was adjusted."""
    if refold:
        return {"decided_by": DELEGATED_BY, "implemented_by": IMPLEMENTED_BY, "reason": refold}
    if not at_cterm:
        return {"decided_by": DECIDED_BY, "implemented_by": IMPLEMENTED_BY, "reason": ITM2B_REASON}
    if junction_edit:
        return {
            "decided_by": JUNCTION_DECIDED_BY,
            "implemented_by": IMPLEMENTED_BY,
            "reason": f"{JUNCTION_REASON} ({junction_edit})",
        }
    return {"decided_by": "", "implemented_by": "", "reason": ""}


RELOCATED = {
    # design_id (permanent, survives renumbering): (insert after residue, residues that follow)
    "l1-p0976": (235, "TIKGIQKR"),  # ITM2B human
    "l1-p2672": (235, "TIRGIQKR"),  # ITM2B mouse
}

# --- Where a C-terminal member's cassette goes: one rule, read from the verdict tables --------
#
# Evidence. l1_cterm_refold_verdicts.tsv is a snapshot of l1-tag-refold's verdict on each
# member's C-terminally tagged protein (processed://.../predicted-structures/tagged/2026-11-12/
# verdicts_cterm.tsv); l1_cterm_stalk_refold_verdicts.tsv is its verdict on each stalk-site
# design (PR #301, 8af3b01d). Both join on design_id and are checked against the tagged
# protein's SHA-256. The stalk site is the review's ``stalk_site_after``: the midpoint of a
# stalk of 8+ residues between the last TM helix and the folded domain.
REFOLD_VERDICTS = HERE / "l1_cterm_refold_verdicts.tsv"
STALK_REFOLD_VERDICTS = HERE / "l1_cterm_stalk_refold_verdicts.tsv"

CTERM, STALK, HOLD, RELOCATE = "C-terminal", "stalk", "hold", "relocated"

FLAG_RULE_BY = "Becca Carlson, 2026-11-22 ('Yes, build them'): refold failures flag, never hold"

STALK_RULE = (
    "A C-terminal arm moves to its stalk site when the C-terminally tagged protein fails its "
    "refold verdict (the fold is disturbed, the epitope is buried or packed on the native "
    "domain past the verdict's recorded limit, or a helix forms inside the cassette) and the "
    "stalk-site design passes its own refold. When the stalk-site design fails too, or there is "
    "no usable stalk site, and the failure is in the fold (native-region lDDT or lDDT near the "
    "site under the verdict's bar), the C-terminal arm is built and flagged (Becca, 2026-11-22: "
    "refold failures flag, never hold); a mechanistic failure at both sites is held. It stays at the C-terminus when the refold "
    "passes; when it has no usable stalk site and the fold is intact (epitope packing or a "
    "cassette helix alone), flagged; or when the ectodomain is natively released (a stalk site "
    "lies before the cut and would tag the membrane stub)"
)
REFOLD_REASON = (
    "C-terminal refold failed ({why}); cassette moved to stalk residue {site}, between the "
    "membrane and the folded domain (rule recorded before the verdicts in l1_cterm_decisions.py)"
)


@dataclass(frozen=True)
class Disposition:
    """What is done with one C-terminal member, and the clause of the rule that says so."""

    action: str  # CTERM, STALK, HOLD or RELOCATE
    clause: str
    site: int | None = None  # the residue the cassette follows (stalk site, or decided site)
    reason: str = ""
    #: What the uniform rule alone gives, where a recorded exception overrides it.
    rule_action: str = ""


def released_genes(review: Iterable[dict]) -> set[str]:
    """Genes whose ectodomain UniProt records as cut or released (any member flagged).

    Derived from the saved UniProt snapshot through ``l1_cterm_review.tsv``'s
    ``ectodomain_released`` flag, not typed by hand. Gene level, so an ortholog with no
    UniProt entry (``uniprot_not_checked``) is covered by its human member's record.
    """
    return {r["gene_symbol"] for r in review if "ectodomain_released" in r["flags"].split(";")}


def refold_disposition(
    review: dict, verdict: dict, stalk_verdict: dict | None, released: bool
) -> Disposition:
    """:data:`STALK_RULE` for one member whose C-terminus survives (``cterm_tag == kept``).

    ``verdict`` is the refold verdict of the C-terminally tagged protein, ``stalk_verdict`` the
    verdict of the stalk-site design (None before it has been refolded), ``released`` whether
    the gene's ectodomain is natively released.
    """
    if verdict["verdict"] == "pass":
        return Disposition(CTERM, "refold_pass")
    assert verdict["verdict"] == "fail", f"{review['design_id']}: verdict {verdict['verdict']!r}"
    why = verdict["verdict_reasons"]
    if released:
        return Disposition(CTERM, "released_ectodomain", reason=why)
    if not review["stalk_site_after"]:
        # The verdict names an lDDT only when the fold itself is under the bar
        # (tag_refold/score_tagged_folds.py: LDDT_PASS, LOCAL_LDDT_PASS).
        if "lDDT" in why:
            # Becca, 2026-11-22: a refold failure flags, never holds; the C-terminus is the
            # rung's own first pick.
            return Disposition(
                CTERM,
                "fold_fail_no_stalk_flagged",
                reason=f"refold_flag=fail: {why} ({review['gene_symbol']}: no usable stalk site)",
            )
        return Disposition(CTERM, "no_stalk_site", reason=why)
    site = int(review["stalk_site_after"])
    if stalk_verdict is None:
        return Disposition(STALK, "stalk_refold_pending", site, why)
    if stalk_verdict["verdict"] == "pass":
        return Disposition(STALK, "stalk_refold_pass", site, why)
    assert stalk_verdict["verdict"] == "fail", review["design_id"]
    if "lDDT" in why:
        # Becca, 2026-11-22: neither site passes, so the C-terminus (the rung's first pick) is
        # built, flagged. A mechanistic C-terminal failure (packing, a cassette helix) keeps
        # the hold below.
        return Disposition(
            CTERM,
            "stalk_refold_fail_flagged",
            reason=f"refold_flag=fail: {why} | stalk residue {site}: "
            f"{stalk_verdict['verdict_reasons']}",
        )
    return Disposition(
        HOLD,
        "stalk_refold_fail",
        site,
        f"{review['gene_symbol']}: C-terminal and stalk-site (residue {site}) refolds both fail "
        f"({why} | {stalk_verdict['verdict_reasons']})",
    )


#: Members left as they are on main where :func:`refold_disposition` gives something else.
#: Each is a recorded exception, not a silent difference; :func:`dispositions` fails if the
#: rule comes to agree with one (then the row must go). design_id: (action, flag text, reason).
#: None at present (2026-11-14: the PLVAP isoform's hold is a clause of the rule, not an exception).
RULE_EXCEPTIONS: dict[str, tuple[str, str, str]] = {}

# --- Is the C-terminus outside the cell? -------------------------------------------------------
#: Review flags that mark a weak C-out call (``review_l1_cterm_members.py``).
WEAK_COUT_FLAGS = ("short_ecto_tail", "v2_cterm_cytoplasmic", "v2_no_tm", "no_v2")
COUT_RULE = (
    "DeepTMHMM v1 is the topology of record, so a C-out call is kept unless an independent "
    "source contradicts it. The member is held when (1) UniProt annotates its C-terminus "
    "cytoplasmic on experimental evidence (ECO:0000269); (2) UniProt describes a type I "
    "protein, with a cleaved signal peptide, one TM helix and a cytoplasmic C-terminus, so "
    "that v1's C-out call rests on reading the signal peptide as the only TM helix; (3) the "
    "call is weak (a v1 extracellular tail under 5 residues, or DeepTMHMM v2 disagreeing, "
    "calling no TM helix or having no row) and UniProt annotates the C-terminus cytoplasmic on "
    "any evidence; or (4) v2 puts the C-terminus in the cytoplasm and neither UniProt (an "
    "extracellular annotation at this member's C-terminus) nor a curated lane override "
    "(lane_overrides.tsv) places it outside. An ortholog with no UniProt entry follows its "
    "gene's human member under (1) and (2). A cytoplasmic annotation that is only a "
    "sequence-model prediction (ECO:0000255), on a call where v1 and v2 agree and the protein "
    "is not type I, is kept and flagged"
)
#: Clauses of :func:`cout_disposition` that hold a member.
COUT_HOLD_CLAUSES = (
    "uniprot_experimental_cytoplasmic",
    "uniprot_type_i",
    "uniprot_cytoplasmic",
    "v2_cytoplasmic_unsupported",
    "ortholog_of_held_human",
)
#: Member-specific evidence quoted in the decision row, beyond the review's columns. A paper is
#: named here only after its PubMed record was read (2026-11-14) and the claim found in it.
COUT_NOTES = {
    "l1-p2528": (
        "C-terminal orientation unconfirmed; canonical TRPV1 C-terminus is cytoplasmic; "
        "candidate for the internal rung. The isoform (521 aa) is canonical TRPV1 (UniProt "
        "Q8NER1, 839 aa) residues 1-515 plus six residues of its own, so it ends five residues "
        "into canonical helix S3 (511-532), which follows the cytoplasmic S2-S3 loop (498-510); "
        "it lacks S4-S6, the pore and the canonical cytoplasmic C-terminus (684-839, "
        "ECO:0000269). v1's C-out call rests on a 10-residue third helix (511-520) and one "
        "extracellular residue; v2 calls two helices and a cytoplasmic C-terminus (498-521)"
    ),
    "l1-p1957": (
        "C-terminus cytoplasmic by experiment; candidate for the N at residue 1 rung. The member "
        "is identical to UniProt Q96JW4, whose topological domains all cite PMID 16984228 "
        "(Sahni, Nelson and Scharenberg 2007, Biochem J 401:505, doi:10.1042/BJ20060673): "
        "extracellular 1-162 through cytoplasmic 565-573, the last feature. The abstract states "
        "that heterologously expressed SLC41A2 'is a plasma-membrane protein with an "
        "N-terminus-outside/C-terminus-inside 11-TM (transmembrane)-span topology'. DeepTMHMM "
        "v1 and v2 both call 11 helices in the reverse orientation (N-terminus in, C-terminus "
        "out). Only the abstract was read, not the experiments behind it"
    ),
    "l1-p2342": (
        "C-terminus cytoplasmic if UniProt's type I topology is right; candidate for the N "
        "after SP rung. UniProt Q7Z7D3 (identical to the member) describes a single-pass type "
        "I protein: "
        "signal peptide 1-24, extracellular 25-259, TM 260-280, cytoplasmic 281-282, all "
        "sequence-model annotations (ECO:0000255). DeepTMHMM v1 instead calls residues 7-29 "
        "the only TM helix and everything after it extracellular, so it reads the signal "
        "peptide as a signal anchor and calls no helix in the hydrophobic C-terminal segment "
        "(260-280); v2 does the same. If UniProt is right, a C-terminal cassette sits in the "
        "cytoplasm, two residues after the TM helix. VTCN1 is not in gpi_set.tsv"
    ),
    "l1-p2628": (
        "Held with its human member (coordinator, 2026-11-14). The cyno protein is 98.6% "
        "identical to human VTCN1, with the same C-terminal segment and the same v1 call"
    ),
}
#: Where a held member could go instead. A note for the next cut, not a build: nothing reads it.
COUT_REROUTE = {
    "l1-p1957": "Candidate for re-routing: SLC41A2 to N at residue 1 (multipass, N-out)",
    "l1-p2342": "Candidate for re-routing: VTCN1 to N after signal peptide (type I, SP 1-24)",
    "l1-p2628": "Candidate for re-routing: VTCN1 to N after signal peptide (type I, SP 1-24)",
    "l1-p2528": "Candidate for re-routing: TRPV1 isoform to the internal rung",
}


def cout_disposition(review: dict, human: dict | None = None) -> Disposition:
    """:data:`COUT_RULE` for one member, from its ``l1_cterm_review.tsv`` row.

    ``human`` is the review row of the gene's human member, used only for an ortholog that has
    no UniProt entry of its own (``uniprot_not_checked``).
    """
    flags = set(review["flags"].split(";"))
    weak = [f for f in WEAK_COUT_FLAGS if f in flags]
    uniprot = review["uniprot_cterm_topology"]
    if "lane_override_topology" in flags:
        return Disposition(CTERM, "lane_override", reason="curated row in lane_overrides.tsv")
    if "uniprot_not_checked" in flags and human is not None:
        h = cout_disposition(human)
        if h.clause in ("uniprot_experimental_cytoplasmic", "uniprot_type_i"):
            return Disposition(
                HOLD,
                "ortholog_of_held_human",
                reason=f"no UniProt entry of its own; its human member {human['design_id']} is "
                f"not built ({h.reason})",
            )
    if uniprot.startswith("Cytoplasmic"):
        if "ECO:0000269" in uniprot:
            return Disposition(
                HOLD,
                "uniprot_experimental_cytoplasmic",
                reason=f"UniProt annotates the C-terminus cytoplasmic on experimental "
                f"evidence: {uniprot}",
            )
        if review["uniprot_signal_peptide"] and review["uniprot_tm_helices"] == "1":
            return Disposition(
                HOLD,
                "uniprot_type_i",
                reason=f"UniProt describes a type I protein (signal peptide "
                f"{review['uniprot_signal_peptide']}, one TM helix, C-terminus {uniprot})",
            )
        if weak:
            return Disposition(HOLD, "uniprot_cytoplasmic", reason=f"UniProt: {uniprot}")
        return Disposition(
            CTERM,
            "predicted_cytoplasmic_kept",
            reason="DeepTMHMM v1 and v2 agree C-out; UniProt's cytoplasmic annotation is a "
            "prediction",
        )
    if not weak:
        return Disposition(CTERM, "not_weak")
    if "v2_cterm_cytoplasmic" in flags and not uniprot.startswith("Extracellular"):
        return Disposition(
            HOLD,
            "v2_cytoplasmic_unsupported",
            reason="C-terminal orientation unconfirmed: DeepTMHMM v2 puts the C-terminus in "
            "the cytoplasm and no independent source places it outside",
        )
    support = []
    if review["v2_cterm_side"] == "non-cytoplasmic":
        support.append("DeepTMHMM v2 agrees (C-terminus non-cytoplasmic)")
    if uniprot.startswith("Extracellular"):
        support.append(f"UniProt: {uniprot}")
    return Disposition(
        CTERM,
        "weak_kept",
        reason="; ".join(support) or "DeepTMHMM v1 alone; no source contradicts it",
    )


def cout_dispositions(review: Iterable[dict]) -> dict[str, Disposition]:
    """:func:`cout_disposition` for every member, each ortholog with its gene's human member."""
    review = list(review)
    human = {r["gene_symbol"]: r for r in review if r["species"] == "human" and r["uniprot_acc"]}
    return {r["design_id"]: cout_disposition(r, human.get(r["gene_symbol"])) for r in review}


def _index_unique(rows: Iterable[dict], key: str, label: str) -> dict[str, dict]:
    out: dict[str, dict] = {}
    for row in rows:
        k = (row.get(key) or "").strip()
        if not k:
            raise ValueError(f"{label}: missing/empty {key!r}")
        if k in out:
            raise ValueError(f"{label}: duplicate {key} {k!r}")
        out[k] = row
    return out


def dispositions(
    review: Iterable[dict] | None = None,
    verdicts: dict[str, dict] | None = None,
    stalk_verdicts: dict[str, dict] | None = None,
) -> dict[str, Disposition]:
    """Every C-terminal member's disposition, keyed by design_id, from committed tables.

    Order: a member whose C-terminus is cut off (``cterm_tag == removed``) is built only at a
    decided site (:data:`RELOCATED`) and is otherwise held; then :func:`cout_disposition`;
    then :func:`refold_disposition`, overridden only by :data:`RULE_EXCEPTIONS`.
    """
    review = read_table(REVIEW) if review is None else list(review)
    _index_unique(review, "design_id", REVIEW.name)
    if verdicts is None:
        verdicts = _index_unique(read_table(REFOLD_VERDICTS), "design_id", REFOLD_VERDICTS.name)
    if stalk_verdicts is None:
        stalk_verdicts = _index_unique(
            read_table(STALK_REFOLD_VERDICTS), "design_id", STALK_REFOLD_VERDICTS.name
        )
    released = released_genes(review)
    couts = cout_dispositions(review)
    out, consulted = {}, set()
    for r in review:
        d = r["design_id"]
        assert r["cterm_tag"] in ("kept", "removed"), (
            f"{d}: unexpected cterm_tag {r['cterm_tag']!r}"
        )
        if r["cterm_tag"] == "removed":
            out[d] = (
                Disposition(RELOCATE, "decided_site", RELOCATED[d][0], ITM2B_REASON)
                if d in RELOCATED
                else Disposition(HOLD, "cterm_removed", reason=r["flags"])
            )
            continue
        cout = couts[d]
        if cout.action == HOLD:
            out[d] = cout
            continue
        if d not in verdicts:
            # A member added after the last refold batch: not built until its verdict exists.
            out[d] = Disposition(HOLD, "refold_pending", reason="C-terminal refold not yet run")
            continue
        rule = refold_disposition(
            r, verdicts[d], stalk_verdicts.get(d), r["gene_symbol"] in released
        )
        if rule.site is not None or rule.clause == "stalk_refold_fail_flagged":
            consulted.add(d)
        if d in RULE_EXCEPTIONS:
            action, text, _ = RULE_EXCEPTIONS[d]
            assert action != rule.action, f"{d}: the rule now gives {action}; drop its exception"
            rule = Disposition(action, "exception", reason=text, rule_action=rule.action)
        out[d] = rule
    assert set(RULE_EXCEPTIONS) <= set(out), "an exception names a member not on the rung"
    unused = set(stalk_verdicts) - consulted
    assert not unused, f"stalk-site verdicts the rule never read: {sorted(unused)}"
    return out


def read_table(path: Path) -> list[dict]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def load_inputs() -> tuple[list[dict], dict, dict, dict, dict]:
    """The C-terminal members and their review, encoding, manifest and construct records.

    Every table is joined on design_id, the permanent key; member ids (candidate_id) are
    regenerated by each renumbering and are kept only as plate-order labels. The tag-path
    table has no design_id, so it reaches its design through the manifest of the same commit
    and is checked against the review's protein hash. The returned dicts are keyed by the
    current member id for the builder below.
    """
    from deliverome_analysis.core.s3 import processed

    manifest_rows = [r for r in read_table(MANIFEST) if ld.is_active(r)]
    _index_unique(manifest_rows, "candidate_id", MANIFEST.name)
    did_of = {r["candidate_id"]: ld.design_id(r) for r in manifest_rows}
    cid_of = {d: c for c, d in did_of.items()}
    assert len(cid_of) == len(did_of), "two active members share a design id"
    tp_rows = read_table(TAG_PATHS)
    _index_unique(tp_rows, "candidate_id", TAG_PATHS.name)
    unmapped_tp = sorted(r["candidate_id"] for r in tp_rows if r["candidate_id"] not in did_of)
    assert not unmapped_tp, (
        f"{TAG_PATHS.name}: inactive/missing manifest members: {unmapped_tp[:5]}"
    )
    paths = {did_of[r["candidate_id"]]: r for r in tp_rows}
    review = _index_unique(read_table(REVIEW), "design_id", REVIEW.name)
    cterm = sorted(d for d, r in paths.items() if r["placement"] == "C-terminal")
    assert sorted(review) == cterm, "l1_cterm_review.tsv is not the current C-terminal set"
    for d in cterm:
        assert review[d]["sequence_sha256"] == paths[d]["sequence_sha256"], d
        assert review[d]["candidate_id"] == cid_of[d], f"{d}: review carries a stale member id"
    encodings = _index_unique(read_table(ENCODINGS), "design_id", ENCODINGS.name)
    missing = [d for d in cterm if d not in encodings]
    assert not missing, f"C-terminal members with no issued cassette encoding: {missing}"
    manifest = {ld.design_id(r): r for r in manifest_rows}
    constructs = json.loads(Path(processed.download(CONSTRUCTS_KEY)).read_text())
    by_member = {cid_of[d]: r for d, r in encodings.items() if d in cid_of}
    return (
        [paths[d] for d in sorted(cterm, key=lambda d: cid_of[d])],
        {cid_of[d]: r for d, r in review.items()},
        by_member,
        {cid_of[d]: r for d, r in manifest.items()},
        constructs,
    )


def junction_features(insert: str, lo: int, hi: int) -> list[str]:
    """Vendor complexity features that overlap the cassette span [lo, hi) of the insert."""
    out = []
    for s, e, gc in gc_windows(insert, GC_WINDOW):
        if s < hi and e > lo:
            out.append(f"gc{gc:.2f}@{s}-{e}")
    for j, i, n in direct_repeats(insert, MIN_DIRECT_REPEAT):
        if (j < hi and j + n > lo) or (i < hi and i + n > lo):
            out.append(f"repeat{n}@{j},{i}")
    for s, base, n in homopolymers(insert, MIN_HOMOPOLYMER):
        if s < hi and s + n > lo:
            out.append(f"homopolymer{base}{n}@{s}")
    for s, e, n in hairpins(insert, MIN_STEM, MAX_LOOP):
        if s < hi and e > lo:
            out.append(f"hairpin{n}@{s}-{e}")
    return out


def build_one(
    member: dict,
    orf: str,
    protein: str,
    cassette: str,
    tag_fixed: dict,
    after: int | None = None,
    module: str = "",
) -> dict:
    """Cassette after residue ``after`` (default: the last residue), then TAA.

    ``module`` is DNA appended to the native C-terminus before the cassette (the snorkel's
    linker and TM helix, issued per arm in ``l1_snorkel_module_encodings.tsv``); only with the
    default ``after``.
    """
    cid = member["candidate_id"]
    assert not (module and after is not None), f"{cid}: a module goes only at the C-terminus"
    after = len(protein) if after is None else after
    module_aa = seqtools.translate(module) if module else ""
    assert seqtools.translate(orf) == protein + "*", f"{cid}: ORF does not encode its protein"
    assert seqtools.translate(cassette) == CASSETTE_AA, f"{cid}: cassette is not the 51-aa tandem"
    assert qc.keeps_p003_tag_codons(cassette), f"{cid}: tag codons are not p003's"
    if after == len(protein):
        body, edit = domesticate_junction(orf[:-3])
        head, rest = body + module, ""
    else:
        head, rest, edit = orf[: 3 * after], orf[3 * after : -3], ""
    coding = head + cassette + rest + STOP
    translated = seqtools.translate(coding)
    expect = protein[:after] + module_aa + CASSETTE_AA + protein[after:]
    assert translated == expect + "*", f"{cid}: protein is not native[:n] + cassette + native[n:]"
    assert translated.count("*") == 1, f"{cid}: internal stop"
    insert = qc.wrap_insert(coding)
    assert qc.wrapper_sites_only(insert), (
        f"{cid}: Type IIS sites {seqtools.find_sites(insert, ENZYMES)}"
    )
    lo = qc.CODING_START + len(head)
    hi = lo + len(cassette) + (len(TAIL) if not rest else 0)
    return {
        "candidate_id": cid,
        "design_id": member["design_id"],
        "gene_symbol": member["gene_symbol"],
        "candidate_set": member["candidate_set"],
        "species": member["species"],
        "native_aa": len(protein),
        "protein_aa": len(expect),
        "orf_nt": len(orf),
        "coding_nt": len(coding),
        "insert_nt": len(insert),
        "cassette_after_residue": after + len(module_aa),
        "placement": (
            "snorkel (TM module, then cassette)"
            if module
            else "C-terminal"
            if after == len(protein)
            else "internal, stalk (refold)"
            if member.get("refold_reason")
            else "internal, before furin site"
        ),
        # Whether the insert ends ORF-cassette-stop (C-terminal and snorkel arms) or carries
        # native sequence after the cassette (stalk sites, ITM2B): the brief counts these apart.
        "cassette_ends_orf": "no" if rest else "yes",
        # Blank on rule-built rows: only a manual override carries an attribution.
        **attribution(after == len(protein), edit, member.get("refold_reason", "")),
        "cassette_start_nt": lo + 1,
        "junction_edit": edit,
        "type_iis_sites": "Esp3I x2 (wrapper)",
        "cassette_features": ";".join(junction_features(insert, lo - 15, hi)),
        "protein_sha256": hashlib.sha256(expect.encode()).hexdigest(),
        "insert": insert,
        "protein": expect,
        "upstream_of_barcode": head + cassette + rest + TAIL,
    }


def build_snorkels(
    manifest: dict, constructs: dict, encodings: dict, proteins: dict, tag_fixed: dict
) -> list[dict]:
    """Snorkel arms for the internal-rung members in ``l1_snorkel.built_designs()``.

    The designs are the per-arm decisions' snorkel rows when that table exists (each keyed on
    its permanent ``arm_id``), else the pinned snorkel list. The translated protein must hash to
    the design's ``snorkel_protein_sha256`` (for per-arm designs, the protein that was folded)
    and equal its ``snorkel_protein`` where the list gives one; the member's designed protein
    must equal the design's native sequence (by SHA-256); the member must be on the internal rung.
    """
    designs = l1_snorkel.built_designs()
    paths = _index_unique(read_table(TAG_PATHS), "candidate_id", TAG_PATHS.name)
    cid_of = {
        did: r["candidate_id"]
        for did, r in _index_unique(read_table(MAPPING), "design_id", MAPPING.name).items()
    }
    rows = []
    seen_arm_ids: set[str] = set()
    for design, d in sorted(designs.items()):
        assert design in cid_of, f"{design}: snorkel design missing from {MAPPING.name}"
        cid = cid_of[design]
        assert cid in paths, f"{cid}: snorkel member missing from {TAG_PATHS.name}"
        tp = paths[cid]
        assert tp["placement"] == "internal", f"{cid}: snorkel member not on the internal rung"
        assert tp["sequence_sha256"] == d["sequence_sha256"], f"{cid}: native protein changed"
        assert cid in encodings, f"{cid}: no issued cassette encoding (run l1_swap_risk.py)"
        member = {**tp, "design_id": design}
        row = build_one(
            member,
            constructs[manifest[cid]["key"]]["construct"],
            proteins[cid]["sequence"],
            encodings[cid]["cassette_nt"],
            tag_fixed,
            module=module_dna(cid),
        )
        if "snorkel_protein" in d:
            assert row["protein"] == d["snorkel_protein"], f"{cid}: protein differs from the design"
        assert row["protein_sha256"] == d["snorkel_protein_sha256"], (
            f"{cid}: translated snorkel is not the protein that was folded"
        )
        arm_id = (d.get("arm_id") or "").strip()
        if arm_id:
            assert arm_id not in seen_arm_ids, f"duplicate snorkel arm_id {arm_id!r}"
            seen_arm_ids.add(arm_id)
        row.update(
            {
                "arm_id": arm_id,
                "decided_by": d["decided_by"],
                "implemented_by": IMPLEMENTED_BY
                + " (module DNA: l1_snorkel_module_encodings.tsv, issued per arm 2026-11-21)",
                "reason": "snorkel: " + d["reason"],
                "review_flags": "",
            }
        )
        rows.append(row)
    assert len(rows) == len(designs), f"built {len(rows)} of {len(designs)} snorkel designs"
    return rows


def vector_check(rows: list[dict]) -> int:
    """Every arm in the p0039 backbone: backbone sites plus the two wrapper Esp3I, no more."""
    import candidate_inputs
    from Bio import SeqIO

    bb = str(SeqIO.read(candidate_inputs.backbone(), "genbank").seq).upper()
    return qc.check_inserts_in_backbone([(r["candidate_id"], r["insert"]) for r in rows], bb)


def swap_check(rows: list[dict], manifest: dict, constructs: dict) -> dict:
    """Worst-partner tract up to the barcode, over every active L1 member, on built DNA.

    C-terminal members carry their arm; every other member ends on its own domesticated ORF
    (its cassette is upstream), exactly as ``l1_swap_risk.py`` models them.
    """
    built = {r["candidate_id"]: r["upstream_of_barcode"] for r in rows}
    assert len(built) == len(rows), "duplicate candidate_id in C-terminal/snorkel rows"
    ids, upstream = [], []
    for cid, m in sorted(manifest.items()):
        if m["library"] != ld.LIBRARY_L1:
            continue
        orf = constructs[m["key"]]["construct"]
        ids.append(cid)
        upstream.append(built.get(cid, orf[:-3] + TAIL))
    worst, partner, _ = worst_partner(upstream)
    cterm = [i for i, c in enumerate(ids) if c in built]
    for i in cterm:
        built_row = next(r for r in rows if r["candidate_id"] == ids[i])
        built_row["worst_partner_tract_nt"] = int(worst[i])
        built_row["worst_partner"] = ids[partner[i]] if partner[i] >= 0 else ""
    tracts = [int(worst[i]) for i in cterm]
    return {
        "members_measured": len(ids),
        "cterm_median_nt": sorted(tracts)[len(tracts) // 2],
        "cterm_max_nt": max(tracts),
        "cterm_over_preferred_80": sum(t > PREFERRED for t in tracts),
        "cterm_over_gate_100": sum(t > GATE for t in tracts),
    }


COLUMNS = [
    "candidate_id",
    "design_id",
    "gene_symbol",
    "candidate_set",
    "species",
    "native_aa",
    "protein_aa",
    "orf_nt",
    "coding_nt",
    "insert_nt",
    "placement",
    "cassette_ends_orf",
    "cassette_after_residue",
    "cassette_start_nt",
    "junction_edit",
    "type_iis_sites",
    "cassette_features",
    "worst_partner_tract_nt",
    "worst_partner",
    "review_flags",
    "decided_by",
    "implemented_by",
    "reason",
    "protein_sha256",
    # The permanent arm id (l1-pNNNN-aK) of a snorkel arm built from the per-arm decisions; blank
    # on the C-terminal rung's own arms. The FASTA header leads with it where present.
    "arm_id",
]


def build(pending_ok: bool = False) -> tuple[list[dict], dict]:
    """Every arm and the summary. ``pending_ok`` builds a stalk-site design that has no refold
    verdict yet, only so that it can be exported for refolding."""
    members, review, encodings, manifest, constructs = load_inputs()
    proteins = ld.designed_proteins(MANIFEST, library=ld.LIBRARY_L1)
    tag_fixed = qc.p003_tag_codons()
    rows, held = [], []
    snorkel_rows = build_snorkels(manifest, constructs, encodings, proteins, tag_fixed)
    verdicts = _index_unique(read_table(REFOLD_VERDICTS), "design_id", REFOLD_VERDICTS.name)
    stalk_verdicts = _index_unique(
        read_table(STALK_REFOLD_VERDICTS), "design_id", STALK_REFOLD_VERDICTS.name
    )
    disposition = dispositions(review.values(), verdicts, stalk_verdicts)

    def hold(m: dict, design: str, flags: str) -> None:
        held.append(
            {
                "candidate_id": m["candidate_id"],
                "design_id": design,
                "gene_symbol": m["gene_symbol"],
                "flags": flags,
            }
        )

    for m in members:
        cid = m["candidate_id"]
        design = review[cid]["design_id"]
        native = proteins[cid]["sequence"]
        d = disposition[design]
        if d.clause == "refold_pending":
            hold(m, design, "refold_pending: " + d.reason)
            continue
        # Every verdict the rule read is for this member's own C-terminally tagged protein
        # (a relocated member's verdict is for its built, relocated protein instead).
        # (a member held before any refold, e.g. on its C-out call, may have none yet)
        v = verdicts.get(design)
        tagged = hashlib.sha256((native + CASSETTE_AA).encode()).hexdigest()
        if v is not None and review[cid]["cterm_tag"] != "removed":
            assert v["tagged_sha256"] == tagged, (
                f"{cid}: the C-terminal refold verdict is for a different or missing tagged protein"
            )
        verdict_reasons = v["verdict_reasons"] if v is not None else ""
        after, refold = None, ""
        if d.action == RELOCATE:
            after, motif = RELOCATED[design]
            assert native[after : after + len(motif)] == motif, f"{cid}: furin context moved"
        elif d.clause == "cterm_removed":
            hold(m, design, d.reason)
            continue
        elif d.clause in COUT_HOLD_CLAUSES:
            hold(m, design, "cout_unconfirmed: " + d.reason)
            continue
        elif d.clause in ("exception", "fold_fail_no_stalk"):
            assert d.action == HOLD, f"{cid}: only hold exceptions are implemented"
            hold(m, design, "refold_fail: " + d.reason)
            continue
        elif d.clause in ("fold_fail_no_stalk_flagged", "stalk_refold_fail_flagged"):
            refold = d.reason  # built at the C-terminus, the refold verdict carried as a flag
        elif d.site is not None:  # the stalk site: built, then kept or held on its own verdict
            assert pending_ok or d.clause != "stalk_refold_pending", (
                f"{cid}: its stalk-site design (after residue {d.site}) has no refold verdict; "
                "export it with --export-proteins, refold it, and add its verdict row"
            )
            after = d.site
            refold = REFOLD_REASON.format(why=verdict_reasons, site=after)
        orf = constructs[manifest[cid]["key"]]["construct"]
        m = {**m, "design_id": design, "refold_reason": refold}
        row = build_one(m, orf, native, encodings[cid]["cassette_nt"], tag_fixed, after)
        row["review_flags"] = review[cid]["flags"]
        if d.site is not None and d.action != RELOCATE and design in stalk_verdicts:
            assert stalk_verdicts[design]["tagged_sha256"] == row["protein_sha256"], (
                f"{cid}: the stalk-site verdict is for a different tagged protein"
            )
        if d.action == HOLD:
            hold(m, design, "refold_fail: " + d.reason)
            continue
        rows.append(row)
    rows += snorkel_rows
    assert len({r["candidate_id"] for r in rows}) == len(rows), "duplicate candidate_id in rows"
    assert len({r["design_id"] for r in rows}) == len(rows), "duplicate design_id in rows"
    non_empty_arms = [r["arm_id"] for r in rows if r.get("arm_id")]
    assert len(non_empty_arms) == len(set(non_empty_arms)), (
        "duplicate arm_id in C-term/snorkel rows"
    )
    assert len(rows) + len(held) == len(members) + len(snorkel_rows), (
        f"silent row drop: built {len(rows)} + held {len(held)} != {len(members) + len(snorkel_rows)}"
    )
    n_bb = vector_check(rows)
    swap = swap_check(rows, manifest, constructs)
    summary = {
        "cterm_members": len(members),
        "snorkel_members": len(snorkel_rows),
        "built": len(rows),
        "held": held,
        "relocated": {
            r["candidate_id"]: f"{r['gene_symbol']}: cassette after residue "
            f"{r['cassette_after_residue']} (before the furin site)"
            for r in rows
            if r["placement"] == "internal, before furin site"
        },
        "refold_moves": {
            r["candidate_id"]: f"{r['gene_symbol']}: cassette after stalk residue "
            f"{r['cassette_after_residue']}"
            for r in rows
            if r["placement"] == "internal, stalk (refold)"
        },
        "by_set": dict(Counter(r["candidate_set"] for r in rows)),
        "junction_edits": {
            r["candidate_id"]: f"{r['gene_symbol']}: {r['junction_edit']}"
            for r in rows
            if r["junction_edit"]
        },
        "type_iis": f"every insert: 2 Esp3I (wrapper), 0 BbsI, 0 BsaI; in p0039: backbone's "
        f"{n_bb} + 2",
        "inserts_with_cassette_features": sum(bool(r["cassette_features"]) for r in rows),
        "insert_nt_range": [min(r["insert_nt"] for r in rows), max(r["insert_nt"] for r in rows)],
        "barcode": "20-N placeholder; barcodes not reserved (issue #1)",
        "swap": swap,
    }
    return rows, summary


def write(rows: list[dict], summary: dict, tsv: Path, fasta: Path, js: Path) -> None:
    with tsv.open("w", newline="") as fh:
        w = csv.DictWriter(
            fh, fieldnames=COLUMNS, delimiter="\t", extrasaction="ignore", lineterminator="\n"
        )
        w.writeheader()
        w.writerows(rows)
    fasta.write_text(
        "".join(
            f">{r.get('arm_id') or r['candidate_id']}|{r['gene_symbol']}|{r['species']}|"
            f"{'cterm' if r['placement'] == 'C-terminal' else 'after' + str(r['cassette_after_residue'])}"
            f"|{r['insert_nt']}nt\n"
            + "\n".join(r["insert"][i : i + 70] for i in range(0, len(r["insert"]), 70))
            + "\n"
            for r in rows
        )
    )
    js.write_text(json.dumps(summary, indent=2) + "\n")


def export_proteins(rows: list[dict], out: Path) -> None:
    """Tagged proteins for a structure check: does the cassette leave the native fold intact?"""
    out.mkdir(parents=True, exist_ok=True)
    review = {r["candidate_id"]: r for r in read_table(REVIEW)}
    (out / "l1_cterm_tagged_proteins.fasta").write_text(
        "".join(
            f">{r['design_id']}|{r['candidate_id']}|{r['gene_symbol']}|{r['species']}\n"
            + "\n".join(r["protein"][i : i + 70] for i in range(0, len(r["protein"]), 70))
            + "\n"
            for r in rows
        )
    )
    cols = [
        "design_id",
        "candidate_id",
        "gene_symbol",
        "candidate_set",
        "species",
        "placement",
        "cassette_after_residue",
        "native_aa",
        "protein_aa",
        "cassette_start_aa",
        "cassette_end_aa",
        "terminal_plddt_native",
        "review_flags",
        "protein_sha256",
    ]
    with (out / "l1_cterm_tagged_proteins.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for r in rows:
            a = int(r["cassette_after_residue"])
            w.writerow(
                [
                    r["design_id"],
                    r["candidate_id"],
                    r["gene_symbol"],
                    r["candidate_set"],
                    r["species"],
                    r["placement"],
                    a,
                    r["native_aa"],
                    r["protein_aa"],
                    a + 1,
                    a + len(CASSETTE_AA),
                    review.get(r["candidate_id"], {}).get("terminal_plddt", ""),
                    r["review_flags"],
                    r["protein_sha256"],
                ]
            )
    print(f"wrote {len(rows)} tagged proteins to {out}")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--check", action="store_true", help="rebuild and compare with the record")
    ap.add_argument("--publish", action="store_true", help="upload the FASTA and table to S3")
    ap.add_argument(
        "--export-proteins",
        type=Path,
        metavar="DIR",
        help="write the tagged proteins (FASTA + table keyed by design_id) for refolding",
    )
    args = ap.parse_args(argv)
    rows, summary = build(pending_ok=bool(args.export_proteins))
    if args.export_proteins:
        export_proteins(rows, args.export_proteins)
        return 0
    if args.check:
        tmp = Path(tempfile.mkdtemp())
        write(rows, summary, tmp / "t.tsv", tmp / "f.fasta", tmp / "s.json")
        for new, rec in (
            (tmp / "t.tsv", OUT_TSV),
            (tmp / "f.fasta", OUT_FASTA),
            (tmp / "s.json", OUT_JSON),
        ):
            assert new.read_bytes() == rec.read_bytes(), f"{rec.name} does not reproduce"
        print(
            f"all {len(rows)} arms rebuild byte for byte "
            f"({len(summary['relocated'])} relocated, {len(summary['held'])} held)"
        )
        return 0
    write(rows, summary, OUT_TSV, OUT_FASTA, OUT_JSON)
    print(json.dumps({k: v for k, v in summary.items() if k != "held"}, indent=2))
    for h in summary["held"]:
        print(f"HELD {h['candidate_id']} {h['gene_symbol']}: {h['flags']}")
    if args.publish:
        from deliverome_analysis.core.s3 import processed

        for path in (OUT_TSV, OUT_FASTA, OUT_JSON, REVIEW):
            processed.upload(path, PUBLISH_PREFIX + path.name)
            print(f"uploaded processed://{PUBLISH_PREFIX}{path.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
