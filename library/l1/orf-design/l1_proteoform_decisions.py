"""Per-member decisions on every L1 member the proteoform audit flags (``l1_proteoform_audit.tsv``).

Requested by Becca on 2026-11-11 ("review each individual sequence and record overrides or
modifications to mane or appris"). Reviewed by hand and implemented by the Claude Code session
l1-proteoform-audit (branch ``claude/l1-proteoform-audit-c4b011``, PR 300). This script writes
the decision table only; it changes no manifest, panel or CDS record. The single writer for
those (the l1-reviewer-package session, then l1-integration) applies every row.

Every row is ``apply``. On 2026-11-12 Becca delegated protein-level calls ("try to not queue
protein specific things for management decisions ... use your best judgment to complete a new
and improved evidence based l1 completely by the morning"; OVERNIGHT-L1-2026-11-12.md, "RULE
CHANGE"). So each judgement is decided here on the evidence and recorded as ``decided_by =
l1-proteoform-audit (delegated by Becca 2026-11-12, "use your best judgment")``, with the
alternative that was rejected in ``alternative_rejected`` so the call can be reversed. The one
objective defect (OR8G3's non-ATG start) is attributed to the overnight hard rule, and the four
GPI-form keeps to Becca's own earlier overrides. A "keep" changes nothing in the manifest.

Each decision cites the record it rests on: the UniProt accession or isoform the member's
protein equals (exact sequence match against the pinned UniProt 2026-10-11 isoform table), the
Ensembl 115 transcript that encodes a replacement, the existing L2 member that already
carries it, and the audit's DeepTMHMM v1 / Pfam 39.0 calls. No literature is cited: no claim
below needed more than the sequence records.

    uv run --locked python library/l1/orf-design/l1_proteoform_decisions.py
"""

from __future__ import annotations

import csv
from pathlib import Path

HERE = Path(__file__).resolve().parent
AUDIT = HERE / "l1_proteoform_audit.tsv"
OUT = HERE / "l1_proteoform_decisions.tsv"
#: Post-apply check: what each decision did to the current manifest (rerun after every apply).
CHECK = HERE / "l1_proteoform_decisions_check.tsv"
TRIM_LEDGER = HERE / "ortholog_start_trims.tsv"
#: The member columns describe the member *as decided* (audit at commit 8fbcd342, before the
#: writer applied anything); once a row is in OUT, those columns are read back from it rather
#: than from the current audit, where an applied change has already cleared the flags.
SNAPSHOT_COLUMNS = (
    "candidate_id",
    "member_id",
    "key",
    "candidate_set",
    "species",
    "gene_symbol",
    "transcript",
    "flags",
    "flag_class",
)

IMPLEMENTED_BY = (
    "Claude Code session l1-proteoform-audit, PR #300 (branch claude/l1-proteoform-audit-c4b011)"
)
RULE = "rule (objective defect; overnight hard rules, Becca 2026-11-11)"
DELEGATED = 'l1-proteoform-audit (delegated by Becca 2026-11-12, "use your best judgment")'
GPI_OVERRIDE = "Becca Carlson (GPI-form override, reference_overrides.tsv, before 2026-11-10)"

# Reasons shared by many members.
LADDER = (
    "MANE Select and the APPRIS principal agree on this protein, so the reference ladder "
    "(select_reference.py rung 1) chose it; the UniProt canonical is the odd one out. The flag "
    "is a difference from UniProt, not a defect."
)
MANE_DEFAULT = (
    "MANE Select kept by the ladder: MANE and APPRIS disagree and observed evidence did not "
    "overturn MANE (select_reference.py rung 3)."
)
APPRIS_ONLY = (
    "The gene has no MANE Select; the APPRIS principal is the ladder's pick (rung 2). "
    "It keeps its signal peptide, TM count and orientation per DeepTMHMM v1."
)
SP_VARIANT = (
    "The member's signal peptide differs in sequence from the UniProt canonical's (alternative "
    "first exon or N-terminal residues), but DeepTMHMM v1 still calls a cleaved signal peptide "
    "on the member and the mature chain is unchanged in orientation. Not a lost signal peptide."
)
PFAM_THRESHOLD = (
    "Coverage of the human reference is >= 0.97 and the family is still present at a lower "
    "copy number, below the gathering threshold, or as a clan-mate family on the same segment; "
    "the domain is not missing from the sequence."
)
ORTHO_DIVERGENT = (
    "Compara 116 one-to-one ortholog covering >= 97% of the human reference with no N- or "
    "C-terminal gap; 54-60% identity is ordinary mouse-human divergence for this gene, not a "
    "gene-model error. Topology and signal peptide match human (DeepTMHMM v1)."
)


def ortho_ext(met: int) -> str:
    return (
        f"The cyno Ensembl 116 model starts {met - 1} residues upstream of the residue aligned to "
        f"human Met1, which is an in-frame, native ATG (ortholog Met{met}). The extension is not "
        "in the human protein, carries no signal peptide (DeepTMHMM v1) and would sit on the "
        "extracellular N-terminus, where the tag goes. Ensembl 116 has no other model for this "
        f"cyno gene. Starting at Met{met} makes the protein exactly as long as the human "
        "reference, which marks the extension as a gene-model error. The repair is a 5' trim "
        "to a native codon; no codon is edited. It is an explicit, recorded exception to the "
        "ortholog CDS policy's ban on implicit start repairs."
    )


#: design_id -> (decision, target, status, decided_by, reason, sources)
DECISIONS: dict[str, tuple[str, str, str, str, str, str]] = {
    # ---- apply: objective defect --------------------------------------------------------
    "l1-p1525": (
        "repair",
        "codon 1 GTG -> ATG (protein becomes UniProt P0DMU2 exactly)",
        "apply",
        RULE,
        "Both Ensembl 115 OR8G3 transcripts (ENST00000641572, ENST00000641848) start with GTG, "
        "so the ordered protein begins Val. Replacing codon 1 with ATG gives the reviewed "
        "UniProt sequence P0DMU2 residue for residue. Same remedy Becca approved for SLC22A17 "
        "and TRPV6 on 2026-11-12 (reference_changes.tsv); the non-ATG adjudication found no "
        "start-supported alternative transcript.",
        "UniProt P0DMU2; Ensembl 115 ENST00000641572, ENST00000641848; "
        "non_atg_alternative_adjudication.tsv",
    ),
    # ---- changes: references ------------------------------------------------
    "l1-p2284": (
        "swap",
        "l2-p1634 (iso:TSPAN10:Q9H1Z9; Ensembl 115 ENST00000611590 and four others encode it)",
        "apply",
        DELEGATED,
        "The APPRIS-only reference ENST00000696047 is a 235-aa Ensembl model with no UniProt "
        "isoform. It keeps 3 of the canonical's 4 TM helices (DeepTMHMM v1) and stops 123 "
        "residues early, so it loses TM4 and the C-terminal part of the large extracellular loop. The "
        "355-aa UniProt canonical is encoded by five Ensembl 115 transcripts and already sits in L2.",
        "UniProt Q9H1Z9; Ensembl 115 ENST00000696047, ENST00000611590",
    ),
    "l1-p0513": (
        "swap",
        "l2-p0455 (iso:DISP2:A7MBM2; Ensembl 115 ENST00000267889)",
        "apply",
        DELEGATED,
        "The APPRIS-only reference ENST00000949524 is an Ensembl-only 1,184-aa model with no "
        "UniProt isoform. 217 canonical residues are unaligned internally and 2 of its 12 TM "
        "helices (DeepTMHMM v1: 10 vs 12). The 1,401-aa canonical has its own Ensembl 115 "
        "transcript and already sits in L2.",
        "UniProt A7MBM2; Ensembl 115 ENST00000949524, ENST00000267889",
    ),
    "l1-p1106": (
        "swap",
        "UniProt canonical Q96KR4 via Ensembl 115 ENST00000330198 (no L2 member carries it)",
        "apply",
        DELEGATED,
        "The APPRIS-only reference ENST00000482695 (452 aa, no UniProt isoform) starts at "
        "canonical residue 233. It lacks the signal peptide DeepTMHMM v1 calls on the canonical "
        "(43 aa) and canonical residues 1-232, so as ordered it has no route to the secretory "
        "pathway (v1: no SP, no TM). The 647-aa canonical has its own Ensembl 115 transcript. "
        "A new CDS has to be designed for it.",
        "UniProt Q96KR4; Ensembl 115 ENST00000482695, ENST00000330198",
    ),
    "l1-p0578": (
        "repair",
        "start at native Met95 (5' trim; protein becomes UniProt Q9UF33 exactly)",
        "apply",
        DELEGATED,
        "MANE and APPRIS agree on a 1,130-aa protein that adds 94 residues before the signal "
        "peptide of the 1,036-aa UniProt canonical. With the extension DeepTMHMM v1 no longer "
        "calls a signal peptide: it reads the signal peptide as a second TM helix and puts the "
        "N-terminus inside (2 TM, N-in; canonical: SP + 1 TM, N-out). The mouse and cyno models "
        "in L2 are also 1,130 aa, so the extension may be translated, but nothing in this audit shows "
        "the extended protein reaching the surface. The 1,036-aa reviewed UniProt form has a "
        "clean, cleavable signal peptide at residue 1 and the canonical type I topology, which is "
        "the form a surface-display library needs. Met95 is a native ATG, so the repair is a "
        "5' trim with no codon edit.",
        "UniProt Q9UF33; Ensembl 115 ENST00000389672 (MANE Select)",
    ),
    "l1-p2038": (
        "keep",
        "",
        "apply",
        DELEGATED,
        "Flagged for 14 vs 12 TM helices against F5H094, the readthrough locus's canonical. The "
        "ordered protein is UniProt Q9NPD5-2, an SLCO1B3 isoform, while lib:SLCO1B3 is a "
        "separate L1 reference and F5H094 itself is an active L1 isoform (iso:SLCO1B3-SLCO1B7:"
        "F5H094). Q9NPD5-2 is a reviewed UniProt isoform and differs from the SLCO1B3 reference, "
        "so the slot orders a real, distinct OATP1B3 form rather than a duplicate. Kept, "
        "consistent with the ladder (APPRIS PRINCIPAL:5, the only candidate).",
        "UniProt Q9NPD5-2, F5H094; Ensembl 115 ENST00000540229",
    ),
    # ---- changes: isoforms ----------------------------------------------------
    "l1-p2420": (
        "redraw",
        "replace within the isoform stratum (as the cohort-gate redraw does)",
        "apply",
        DELEGATED,
        "UniProt Q9BXN2-9 is 77 aa: the cytoplasmic tail, TM and a stalk fragment, with no "
        "C-type lectin domain (Pfam PF00059 absent). Becca named CLEC7A as a fragment isoform "
        "on 2026-11-11. It passes the cohort gate because the gate tests only the N-terminus "
        "and first TM.",
        "UniProt Q9BXN2-9",
    ),
    "l1-p2477": (
        "redraw",
        "replace within the isoform stratum",
        "apply",
        DELEGATED,
        "ENST00000469376 is a 189-aa OPN3 fragment with no UniProt isoform: it keeps 4 of the "
        "reference's 7 TM helices (DeepTMHMM v1) and stops 207 residues early. A half GPCR is "
        "not a surface isoform worth a library slot.",
        "Ensembl 115 ENST00000469376; UniProt Q9H1Y3 (reference)",
    ),
    "l1-p2511": (
        "redraw",
        "replace within the isoform stratum",
        "apply",
        DELEGATED,
        "UniProt P48066-2 is 208 aa, 3 of the transporter's 12 TM helices (DeepTMHMM v1), "
        "33% of the reference. Becca named SLC6A11 as a fragment isoform on 2026-11-11.",
        "UniProt P48066-2",
    ),
    "l1-p3119": (
        "redraw",
        "replace within the isoform stratum",
        "apply",
        DELEGATED,
        "UniProt O15438-5 (design l1-p3119, formerly L2 l2-p0011, moved into L1 by the 2026-11-12 redraw) is 572 aa "
        "of a 1,527-aa ABC transporter: 10 of 17 TM helices and neither nucleotide-binding "
        "domain (Pfam PF00005 x2 and one PF00664 absent). It cannot form a functional "
        "transporter.",
        "UniProt O15438-5",
    ),
    "l1-p2470": (
        "redraw",
        "replace within the isoform stratum",
        "apply",
        DELEGATED,
        "ENST00000917943 (537 aa, no UniProt isoform) differs from the L1 reference by one "
        "residue deleted inside the signal peptide (DeepTMHMM v1 SP 30 vs 31). The mature "
        "protein is therefore identical to the reference's, so the slot adds no isoform "
        "information.",
        "Ensembl 115 ENST00000917943; UniProt Q92692 (reference)",
    ),
    # ---- changes: orthologs ---------------------------------------------------
    "l1-p2938": (
        "repair",
        "start at native Met32 (5' trim; protein length then equals human)",
        "apply",
        DELEGATED,
        ortho_ext(32),
        "Ensembl 116 ENSMFAT00000030386; human L1 reference l1-p0375",
    ),
    "l1-p2959": (
        "repair",
        "start at native Met28 (5' trim; protein length then equals human)",
        "apply",
        DELEGATED,
        ortho_ext(28),
        "Ensembl 116 ENSMFAT00000006686; human L1 reference l1-p0248",
    ),
    "l1-p2963": (
        "repair",
        "start at native Met33 (5' trim; protein length then equals human)",
        "apply",
        DELEGATED,
        ortho_ext(33),
        "Ensembl 116 ENSMFAT00000018821; human L1 reference l1-p1440",
    ),
    "l1-p2962": (
        "repair",
        "start at native Met215 (5' trim; protein length then equals human)",
        "apply",
        DELEGATED,
        ortho_ext(215)
        + " The 214-residue extension is Lys/Arg/Gly-rich and is a gene-model error.",
        "Ensembl 116 ENSMFAT00000064353; human L1 reference l1-p1292",
    ),
}

# ---- keeps ---------------------------------------------------------------------------
for did in (
    "l1-p0034 l1-p0245 l1-p0259 l1-p0516 l1-p0582 l1-p0688 l1-p0752 l1-p0766 l1-p0964 "
    "l1-p1003 l1-p1307 l1-p1386 l1-p1495 l1-p1526 l1-p1529 l1-p1545 l1-p1660 l1-p1764 "
    "l1-p1918 l1-p2011 l1-p2214 l1-p0184 l1-p0162 l1-p1799 l1-p1861 l1-p1862"
).split():
    DECISIONS[did] = (
        "keep",
        "",
        "apply",
        DELEGATED,
        LADDER,
        "reference_selection.tsv; Ensembl 115 pep",
    )
DECISIONS["l1-p0128"] = (
    "keep",
    "",
    "apply",
    DELEGATED,
    LADDER + " DeepTMHMM v1 calls the MANE form (UniProt A6NM10-2, 12 extra residues) 7 TM, "
    "N-out, against 6 TM, N-in for the canonical; aquaporins are 6-TM, N-in proteins, so the "
    "v1 call on this member is suspect. Any tag placement for AQP12B should not rest on it.",
    "UniProt A6NM10, A6NM10-2; Ensembl 115 ENST00000407834 (MANE), ENST00000621682",
)
for did in ("l1-p2363", "l1-p2367"):
    DECISIONS[did] = (
        "keep",
        "",
        "apply",
        DELEGATED,
        "Expression-rescue arm: Igk leader plus the same receptor as the MANE / APPRIS "
        "reference, which the ladder chose (MANE and APPRIS agree). The flag compares it with "
        "the longer UniProt canonical.",
        "rescue_constructs.tsv; reference_selection.tsv",
    )
for did in ("l1-p0115", "l1-p1835", "l1-p1720", "l1-p2344"):
    DECISIONS[did] = ("keep", "", "apply", DELEGATED, MANE_DEFAULT, "reference_selection.tsv")
for did in ("l1-p0882", "l1-p0524", "l1-p0329", "l1-p1289", "l1-p1976"):
    DECISIONS[did] = ("keep", "", "apply", DELEGATED, APPRIS_ONLY, "reference_selection.tsv")
DECISIONS["l1-p1264"] = (
    "keep",
    "",
    "apply",
    DELEGATED,
    APPRIS_ONLY + " The 273-aa protein is UniProt P15941-8, MUC1 without the variable "
    "tandem-repeat region; the canonical's repeat array is not synthesizable as a gene "
    "fragment, so the full canonical is not an orderable alternative.",
    "UniProt P15941-8; reference_selection.tsv",
)
DECISIONS["l1-p1160"] = (
    "keep",
    "",
    "apply",
    DELEGATED,
    "APPRIS-only reference (no MANE). DeepTMHMM v1 calls a 41-aa signal peptide and one TM, "
    "N-out, on the member and 0 SP, 2 TM, N-in on the 649-aa canonical, so the member's "
    "predicted architecture is the more plausible surface form of the two. The canonical stays in L2 (l2-p0915).",
    "UniProt Q86X29; Ensembl 115 ENST00000354900, ENST00000621372",
)
for did in ("l1-p0058", "l1-p0227", "l1-p0347", "l1-p1093", "l1-p1096", "l1-p1100", "l1-p1182"):
    DECISIONS[did] = (
        "keep",
        "",
        "apply",
        DELEGATED,
        SP_VARIANT,
        "l1_reference_signal_peptide_audit.tsv; reference_signal_peptides_uniprot.tsv",
    )
for did, form in (
    ("l1-p0017", "P22303-2 (H subunit)"),
    ("l1-p0173", "Q96GW7-2"),
    ("l1-p1275", "P13591-3 (NCAM-120)"),
    ("l1-p2129", "P10646-2 (TFPI-beta)"),
):
    DECISIONS[did] = (
        "keep",
        "",
        "apply",
        GPI_OVERRIDE,
        f"Becca's GPI-form override ({form}); the domains and TM the flag reports as lost are "
        "the ones the soluble or TM canonical carries and the GPI form does not. No change.",
        "library/l1/proteoform-selection/reference_overrides.tsv",
    )
DECISIONS["l1-p2448"] = (
    "keep",
    "",
    "apply",
    DELEGATED,
    "UniProt P17693-3, an HLA-G membrane isoform with the signal peptide, alpha-1 domain, TM "
    "and tail but no alpha-2 / alpha-3 domains. A named, membrane-anchored splice form, not a "
    "fragment; it passes the cohort gate.",
    "UniProt P17693-3",
)
DECISIONS["l1-p2445"] = (
    "keep",
    "",
    "apply",
    DELEGATED,
    "UniProt P42261-6: an alternative N-terminus with its own signal peptide (DeepTMHMM v1, "
    "16 aa) on otherwise reference sequence; passes the cohort gate. Flagged only because the "
    "canonical signal-peptide sequence differs.",
    "UniProt P42261-6",
)
DECISIONS["l1-p2496"] = (
    "keep",
    "",
    "apply",
    DELEGATED,
    "UniProt Q6AZY7-2 (466 of 606 aa) keeps the TM and the region N-terminal to the collagen segments "
    "and lacks the two collagen-like segments, the only Pfam families on the reference. A "
    "named UniProt isoform of reasonable length, not a fragment.",
    "UniProt Q6AZY7-2",
)
for did in ("l1-p2942", "l1-p2956", "l1-p2951"):
    DECISIONS[did] = (
        "keep",
        "",
        "apply",
        DELEGATED,
        "The cyno model extends 10-21 residues beyond human Met1, but the extension precedes a "
        "signal peptide DeepTMHMM v1 still calls, so it is removed with the signal peptide and "
        "the mature protein starts where human does.",
        "audit columns ortholog_n_extension, member_v1_sp",
    )
for did in ("l1-p3022", "l1-p3036", "l1-p3039", "l1-p3040", "l1-p3060"):
    DECISIONS[did] = (
        "keep",
        "",
        "apply",
        DELEGATED,
        ORTHO_DIVERGENT + (" " + PFAM_THRESHOLD if did == "l1-p3036" else ""),
        "audit columns ortholog_identity_to_human, ortholog_human_coverage",
    )
DECISIONS["l1-p3043"] = (
    "keep",
    "",
    "apply",
    DELEGATED,
    "The human OR10G6 reference has 21 N-terminal residues (MLEGVEHLLLLLLLTDVNSKE) that the "
    "mouse protein lacks, and DeepTMHMM v1 calls them a signal peptide. Olfactory receptors "
    "have no signal peptide; the human N-terminus is the anomaly, not the mouse one.",
    "audit columns comparator_n_missing, comparator_v1_sp",
)
for did in ("l1-p2945", "l1-p2957", "l1-p3063", "l1-p3068", "l1-p3084", "l1-p3113"):
    DECISIONS[did] = (
        "keep",
        "",
        "apply",
        DELEGATED,
        PFAM_THRESHOLD,
        "audit columns member_pfam, comparator_pfam",
    )


#: The rejected alternative for each change decision, so Becca can reverse it.
# 2026-11-12, second pass (provenance review for the morning brief): references whose built
# protein is not their ladder rung's and that had no decision row. The audit's
# ladder_source_unrecorded flag finds them.
ALT_ALLELE = (
    "The gene's Ensembl 115 model sits on an alternative-haplotype / novel-patch scaffold "
    "(the ladder's 'no reference, use UniProt canonical' rung). The build is that model's "
    "own CDS, which differs from UniProt canonical {acc} at {diffs}; the residues are "
    "declared in the manifest's variant_residues. Keep: it is an assembled human allele with "
    "a real CDS, while the UniProt sequence has no Ensembl CDS and would have to be edited "
    "in. The audit finds no internal stop and the same length as the canonical."
)
DECISIONS["l1-p1382"] = (
    "keep",
    "",
    "apply",
    DELEGATED,
    "Reference-allele loss of function, same pattern as KIR2DS4 / OR12D1 / OR5AL1. The "
    "ladder's APPRIS pick ENST00000553765 (chr14 primary assembly) is protein_coding_LoF: "
    "Q227 is a stop. The build is ENST00000708727 on patch scaffold HG2526_HG2573_PATCH, "
    "the intact 314-aa allele, identical to UniProt canonical Q8NGC8. Keep the intact allele; "
    "the reference_overrides.tsv row (basis reference_allele_lof) is sent to the writer so "
    "the ladder record names it.",
    "Ensembl 115 pep (ENSP00000451021 vs ENSP00000517307); UniProt Q8NGC8",
)
for did, acc, diffs in (
    ("l1-p1532", "P0C7N1", "3 residues (H20R, V109I, R165C)"),
    ("l1-p1537", "P0C7N8", "1 residue (C193G)"),
    ("l1-p2115", "P59539", "1 residue (C76Y)"),
):
    DECISIONS[did] = (
        "keep",
        "",
        "apply",
        DELEGATED,
        ALT_ALLELE.format(acc=acc, diffs=diffs),
        f"Ensembl 115 pep; UniProt {acc}; candidate_manifest_tranches.tsv variant_residues",
    )

ALTERNATIVE: dict[str, str] = {
    "l1-p1525": "drop OR8G3 (no start-supported alternative transcript exists)",
    "l1-p2284": "keep the 235-aa APPRIS model",
    "l1-p0513": "keep the 1,184-aa APPRIS model",
    "l1-p1106": "keep the 452-aa APPRIS model; or swap to L2 isoform Q96KR4-3 (684 aa)",
    "l1-p0578": "keep the 1,130-aa MANE / APPRIS form with its 94-residue extension",
    "l1-p2038": "drop the readthrough gene's reference",
    "l1-p2420": "keep the 77-aa isoform",
    "l1-p2477": "keep the 189-aa fragment",
    "l1-p2511": "keep the 208-aa fragment",
    "l1-p3119": "keep the 572-aa fragment",
    "l1-p2470": "keep the near-duplicate isoform",
    "l1-p2938": "redraw from the Compara 116 pool (loses the cyno CHRM3 arm)",
    "l1-p2959": "redraw from the Compara 116 pool (loses the cyno CCR5 arm)",
    "l1-p2963": "redraw from the Compara 116 pool (loses the cyno OR3A2 arm)",
    "l1-p2962": "redraw from the Compara 116 pool (loses the cyno NIPAL4 arm)",
    "l1-p0128": "swap to the 295-aa UniProt canonical (ENST00000621682)",
    "l1-p1160": "swap to the 649-aa canonical, l2-p0915",
    "l1-p1382": "build the ladder's APPRIS transcript ENST00000553765 (stop at codon 227)",
    "l1-p1532": "edit the CDS to the UniProt P0C7N1 residues",
    "l1-p1537": "edit the CDS to the UniProt P0C7N8 residue",
    "l1-p2115": "edit the CDS to the UniProt P59539 residue",
    "l1-p1264": "none orderable (the canonical's tandem-repeat array cannot be synthesized)",
}
KEEP_ALTERNATIVE = "replace with the comparator form (UniProt canonical or human reference)"


def _read(path: Path) -> dict[str, dict[str, str]]:
    if not path.exists():
        return {}
    with path.open() as fh:
        return {r["design_id"]: r for r in csv.DictReader(fh, delimiter="\t")}


def check(
    audit: dict[str, dict[str, str]], record: dict[str, dict[str, str]]
) -> list[dict[str, str]]:
    """Each decision against the current audit: did the change land, is a keep still a keep."""
    held = {d: r["status"] for d, r in _read(TRIM_LEDGER).items()}
    rows = []
    for did in sorted(DECISIONS):
        decision = DECISIONS[did][0]
        a = audit.get(did)
        if a is None:
            outcome = "dropped by redraw" if decision == "redraw" else "MISSING"
        elif decision == "keep":
            outcome = "kept"
        elif a["flag_class"] in ("apply", "proposed"):
            outcome = (
                f"pending ({held[did]})" if held.get(did, "").startswith("held") else "NOT APPLIED"
            )
        else:
            outcome = "applied"
        rows.append(
            {
                "design_id": did,
                "gene_symbol": (a or record.get(did, {})).get("gene_symbol", ""),
                "decision": decision,
                "outcome": outcome,
                "flags_now": (a or {}).get("flags", ""),
                "flag_class_now": (a or {}).get("flag_class", ""),
                "translation_check_now": (a or {}).get("translation_check", ""),
            }
        )
    return rows


def _write(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def main() -> int:
    from collections import Counter

    audit = _read(AUDIT)
    record = _read(OUT)
    flagged = {d for d, r in audit.items() if r["flag_class"] in ("apply", "proposed")}
    missing = flagged - set(DECISIONS)
    if missing:
        raise SystemExit(f"flagged members without a decision: {sorted(missing)}")
    gone = {d for d in DECISIONS if d not in audit and d not in record}
    if gone:
        raise SystemExit(f"decided members neither in the audit nor in the record: {sorted(gone)}")
    info = {
        "no_v1_topology",
        "structure_not_ours",
        "duplicate_protein_l2",
        "reference_not_mane_appris_canonical",
    }
    rows = []
    for did in sorted(DECISIONS):
        if did in record:
            member = {c: record[did][c] for c in SNAPSHOT_COLUMNS}
        else:
            a = audit[did]
            member = {c: a[c] for c in SNAPSHOT_COLUMNS}
            member["flags"] = ";".join(f for f in a["flags"].split(";") if f and f not in info)
        decision, target, status, decided_by, reason, sources = DECISIONS[did]
        rows.append(
            {
                "design_id": did,
                **member,
                "decision": decision,
                "target": target,
                "status": status,
                "decided_by": decided_by,
                "implemented_by": IMPLEMENTED_BY,
                "reason": reason,
                "alternative_rejected": ALTERNATIVE.get(did, KEEP_ALTERNATIVE),
                "sources": sources,
            }
        )
    _write(OUT, rows)
    checked = check(audit, _read(OUT))
    _write(CHECK, checked)
    print(len(rows), Counter((r["status"], r["decision"]) for r in rows))
    print(Counter((r["decision"], r["outcome"]) for r in checked))
    bad = [r for r in checked if r["outcome"] in ("MISSING", "NOT APPLIED")]
    if bad:
        raise SystemExit(f"decisions not reflected in the manifest: {bad}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
