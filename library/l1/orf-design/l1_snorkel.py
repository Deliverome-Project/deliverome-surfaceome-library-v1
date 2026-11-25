"""Snorkel constructs for internal-rung L1 members whose C-terminus is cytoplasmic.

A snorkel appends a single TM helix to the protein's own cytoplasmic C-terminus, so the
cassette that follows it is extracellular
([Brown et al., 2013](https://doi.org/10.1371/journal.pone.0073255)). The members and the
protein design are not chosen here. They come from the s5e13 ranker work: PR #245,
``processed://sprint5/s5e13-tag-site-benchmarking/l1_snorkel_designs.tsv``, pinned below by
SHA-256. That table is decided_by ``ranker-pu-segment (delegated by Becca 2026-11-12)``, and
its rows carry their own reason and rejected alternative.

The module after the native C-terminus is 106 aa: GGSGGSGGSGG, then PDGFRB 556-566
(QKKPRYEIRWK), then PDGFRB 530-555 (the TM helix), then GSGGGGS, then the library's 51-aa
cassette. This file supplies the DNA for the 55 residues before the cassette
(:data:`PREFIX_DNA`). The cassette itself is each member's own issued linkers-only encoding
(``l1_cassette_encodings.tsv``), exactly as on the C-terminal rung, so the stretch before the
barcode stays unique per member.

:data:`PREFIX_DNA` was encoded once on 2026-11-12 with ``seqtools.dnachisel_optimize``
(human codon usage, no BbsI/Esp3I/BsaI, GC window, homopolymer limit, seed 20261112) and is
pinned here, so a solver change cannot alter it. It is the encoding the cassette batches were
issued against. Built arms no longer carry it: since 2026-11-21 each snorkel arm has its own
module encoding (``diversify_snorkel_modules.py``, ``l1_snorkel_module_encodings.tsv``; Becca,
2026-11-21: no shared linkers in the middle). Attribution: encoded
by the Claude Code session l1-cterm-build (claude/l1-cterm-build-16eeeb, PR #284). Becca
delegated the call on 2026-11-12, and no person decided this encoding.
"""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path

from deliverome_analysis import seqtools

#: The ranker publishes each snorkel list at its own dated key, which never changes. To adopt a
#: new list, move DESIGNS_KEY and DESIGNS_SHA256 together to the ranker's new dated key.
DESIGNS_KEY = (
    "sprint5/s5e13-tag-site-benchmarking/l1_snorkel_designs/2026-11-14b/l1_snorkel_designs.tsv"
)
# Pinned 2026-11-13 to the ranker's list after the spacer test (PR #245 bf831be8): 26 members.
# The 45-member list below minus 21 whose snorkel arm fails its own refold (a 30-residue spacer
# rescued 0 of 5 tested, and both controls broke), plus P2RX5 (l1-p1552; every designed internal
# site lies within 2 residues of the ATP-pocket ligand contact) and P2RX2 (l1-p1550; its
# contact-free internal site after residue 89 fails its refold, and no other clean site exists).
# Decided by ranker-pu-segment under the coordinator's hold rule (delegated by Becca 2026-11-12).
# Re-pinned 2026-11-14 (PR #245 322a2bc0, 29 members): the same 26 plus GJA8 (l1-p0708), OTOP1
# (l1-p1543) and SLC23A2 (l1-p1894), whose internal sites are all vetoed by verified literature
# ligand contacts; cassette batch 9. No removals, and no protein or design_id changed.
# P2RX2's and the three new snorkel arms have not been refolded yet.
# Re-pinned 2026-11-14 again (PR #331, 28 members): SLC22A24 l1-p1884 leaves (the proteoform-fix
# redesign gives it an internal site after residue 59); its batch-6 encoding stays reserved.
# The other SLC22A24 member, l1-p3149, stays a snorkel. No additions.
DESIGNS_SHA256 = "2b1bb0772611cc553a574564845cffb357f02d425d356ca248481ea9a9933405"

#: The list cassette batch 6 was issued against (PR #245 bb9cd99d, 45 members; the undated key,
#: which the ranker leaves at this table). Kept because an issued encoding is never redrawn or
#: handed on: a member dropped from DESIGNS_KEY keeps its batch-6 encoding as a reserved holder,
#: and replaying batch 6 needs its snorkel flank.
ISSUED_KEY = "sprint5/s5e13-tag-site-benchmarking/l1_snorkel_designs.tsv"
ISSUED_SHA256 = "3ddc9336a3c792e07287424afe18e7590558f7acaba2fdf3037bcbf5ccde9d93"

PREFIX_AA = "GGSGGSGGSGG" + "QKKPRYEIRWK" + "KVVVISAILALVVLTVISLIILIMLW" + "GSGGGGS"
PREFIX_DNA = (
    "GGTGGCTCCGGCGGAAGCGGGGGATCTGGGGGTCAAAAGAAACCCAGGTACGAGATCAGATGGAAGAAAGTCGTGGTAATTTCC"
    "GCCATCCTGGCTCTGGTTGTGCTTACCGTGATCAGTCTCATTATATTGATTATGTTATGGGGCAGCGGAGGGGGCGGCTCA"
)
assert seqtools.translate(PREFIX_DNA) == PREFIX_AA


#: Committed review copies of the two pinned lists (S3 stays their durable home). The issued
#: list is committed in S5 at its own key; the 2026-11-14b list is copied to this experiment's
#: root, the only place git keeps a TSV. Either source must match the pinned SHA-256.
COMMITTED_COPIES = {
    DESIGNS_KEY: Path(__file__).resolve().parent / "l1_snorkel_designs_2026-11-14b.tsv",
    ISSUED_KEY: (
        Path(__file__).resolve().parents[2] / "l1/tag-site-benchmarking" / "l1_snorkel_designs.tsv"
    ),
}


def committed_copy(key: str) -> Path:
    """The committed review copy of a pinned key (it may be absent)."""
    return COMMITTED_COPIES[key]


def _read(key: str, expected_sha256: str) -> dict[str, dict]:
    path = committed_copy(key)
    if not path.exists():
        from deliverome_analysis.s3 import processed

        path = Path(processed.download(key))
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != expected_sha256:
        raise SystemExit(f"{key} changed (sha256 {digest}); re-pin it deliberately")
    with path.open() as fh:
        rows = {r["design_id"]: r for r in csv.DictReader(fh, delimiter="\t")}
    for d, r in rows.items():
        assert r["status"] == "decided", f"{d}: snorkel design not decided"
        assert r["module_aa"][: len(PREFIX_AA)] == PREFIX_AA, f"{d}: module differs from PREFIX_AA"
    return rows


def pinned_designs(expected_sha256: str = DESIGNS_SHA256) -> dict[str, dict]:
    """The snorkel designs built now, keyed by design_id, after checking the file's SHA-256."""
    return _read(DESIGNS_KEY, expected_sha256)


#: The ranker's per-arm decisions (agreed schema, 2026-11-14), final for L1 on 2026-11-15:
#: ``processed://sprint5/s5e13-tag-site-benchmarking/l1_internal_arms_final/2026-11-15e/``
#: (PR #377 9fa02863), committed in S5. Once present, its ``snorkel`` rows with
#: ``decision == build`` are the snorkel arms built (one per member, the snorkel_companion of
#: two loop arms), replacing the pinned list above. Re-pin deliberately with the ranker.
ARM_DECISIONS = (
    Path(__file__).resolve().parents[2]
    / "l1/tag-site-benchmarking"
    / "l1_internal_arm_decisions.tsv"
)
# 2026-11-21: re-pinned after s5e13 l1_internal_additions.py appended SCNN1B's arm (l1-p3174-a1),
# the internal co-receptor addition, then its refold pass; 2026-11-22: rows of retired members
# (the isoform redraw, SCNN1B) dropped, their encodings kept as issued; 16 refold holds built at
# rank 1, flagged (l1_internal_flag_holds.py; Becca 2026-11-22). Later on 2026-11-22: the internal
# members of the 2026-11-22 additions and the GTF MANE Select switches decided (old proteins' rows
# replaced, l1_internal_additions.py --set additions_2026_11_22), EMP1's flagged site corrected
# from 0 to 124, then their refold verdicts applied.
# Later still on 2026-11-22: the v8_prodomino internal sites (Becca: "i want to stay with the v8
# picks exactly") applied with l1_internal_additions.py --apply-delta (1,267 arm rows: 1,119 loop,
# 148 snorkel; refold failures flagged, the 24 mechanistic failures held), then every arm id
# re-assigned through the l1_internal_arm_ids.tsv registry so an arm moved to a new site gets a new K.
# 2026-11-23: the interface/sequon veto rerun on the v8 arms (interface_veto_v8.py, branch
# claude/l1-interface-veto-v8 3b252ac6; Becca: "we should do v8") drops AQP5's second loop arm
# (l1-p0132-a2, after 36) and flags its first; no other arm changes.
# 2026-11-23-pre-freeze: retired 7 dropped internal members, applied 21 zero-GPU clash/fail &
# lit-52 upgrades + 4 Modal clash/propeptide swaps (AREG @193, FZD1 @588, OTOP3 @244, SLC6A20 @536),
# dropped 3 lit-52 hard-vetoed multi-arm sites (SLC1A3 a4, STEAP4 a5, XKR8 a4), and added the 5 new
# internal CDSs (ANO9, TMEM117, CLDN18.2, KCNK2, SCNN1A).
# 2026-11-24: table hygiene & audit flags (PLLP refold_disposition -> fail_flagged, 6 v5 loops
# veto_class -> none, transparent audit_flag entries on SLC44A1/ATRN, 4 GPI omega-proximal arms,
# and 5 close-pair loop pairs; zero decision or construct sequence changes).
ARM_DECISIONS_SHA256 = "928e99acac407ff16f35d59818e4ad2679357326937ec8694306bce050403eab"


def arm_designs(expected_sha256: str = ARM_DECISIONS_SHA256) -> dict[str, dict] | None:
    """Snorkel arms from the per-arm decisions, keyed by design_id, or None without the table.

    Every snorkel uses the same module (``PREFIX_AA`` + the 51-aa cassette), so its protein is
    the member's native protein plus that module; ``snorkel_protein_sha256`` is the SHA-256 the
    ranker folded (``refold_tagged_sha256``), which the builder checks the translated insert
    against. A member carries at most one snorkel arm.
    """
    if not ARM_DECISIONS.exists():
        return None
    digest = hashlib.sha256(ARM_DECISIONS.read_bytes()).hexdigest()
    if digest != expected_sha256:
        raise SystemExit(f"{ARM_DECISIONS.name} changed (sha256 {digest}); re-pin it deliberately")
    out: dict[str, dict] = {}
    with ARM_DECISIONS.open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["arm_kind"] != "snorkel" or r["decision"] != "build":
                continue
            assert r["design_id"] not in out, f"{r['design_id']}: two snorkel arms"
            assert r["refold_tagged_sha256"], f"{r['arm_id']}: no folded protein recorded"
            out[r["design_id"]] = {
                "design_id": r["design_id"],
                "arm_id": r["arm_id"],
                "arm_role": r["arm_role"],
                "gene": r["gene"],
                "sequence_sha256": r["sequence_sha256"],
                "snorkel_protein_sha256": r["refold_tagged_sha256"],
                "refold_verdict": r["refold_verdict"],
                "flag": r.get("flag", ""),
                "status": "decided",
                "reason": r["reason"],
                "decided_by": r["decided_by"],
                "implemented_by": r["implemented_by"],
            }
    return out


def built_designs() -> dict[str, dict]:
    """The snorkel arms built now: the per-arm decisions when present, else the pinned list."""
    arms = arm_designs()
    return arms if arms is not None else pinned_designs()


def reserved_designs() -> dict[str, dict]:
    """Snorkel designs that hold an issued encoding but are no longer built.

    Their encoding stays issued (never redrawn or handed on), so later batches are drawn around
    their snorkel flank. Every member that was ever a built snorkel holds one: the issued list
    (batch 6) and the pinned list (later batches); those not built now are reserved.
    """
    built = built_designs()
    ever = {**_read(ISSUED_KEY, ISSUED_SHA256), **pinned_designs()}
    return {d: r for d, r in ever.items() if d not in built}
