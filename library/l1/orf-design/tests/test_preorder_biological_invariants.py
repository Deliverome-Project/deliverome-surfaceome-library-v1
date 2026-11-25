"""Offline biological and fail-closed pre-order invariants for L1 internal arms.

These checks run strictly offline (zero S3 or network calls) in both the source
analysis repository and the clean-path exported review repository. They guard
against the pre-order failure classes uncovered during the 2026-11-24 audit:

1. GPI omega-proximal internal loop insertions (< 12 aa upstream of omega).
2. Internal loop insertions inside cleaved N-terminal propeptides.
3. Membrane-core ('m'), cytoplasmic ('c'), or unflagged usable==0 internal loops,
   plus cytoplasmic/membrane-core verdicts in the unmasked multipass audit.
4. Same-loop internal arm pairs closer than 25 aa (with review-band 'r' dips kept
   continuous via membrane_mask.masked_regions).
5. Decision-table hygiene (veto_class, refold_disposition enums), silent arm drops
   in topology_vs_experiment_sites.tsv, and rung-filtered snorkel refold counts.
"""

from __future__ import annotations

import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

S4E9 = Path(__file__).resolve().parents[1]
S5E13 = S4E9.parents[1] / "l1/tag-site-benchmarking"
REVIEWER = S4E9 / "reviewer-package"

for _path in (S4E9, S5E13, REVIEWER):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

DECISIONS_TSV = S5E13 / "l1_internal_arm_decisions.tsv"
GPI_SET_TSV = S4E9 / "gpi_set.tsv"
GPI_ANCHORED_GENES_TSV = S4E9 / "gpi_anchored_genes.tsv"
TAG_PATHS_TSV = S4E9 / "l1_tag_paths.tsv"
CTERM_REVIEW_TSV = S4E9 / "l1_cterm_review.tsv"
CONSTRUCTS_CTERM_TSV = S4E9 / "constructs_l1_cterm.tsv"
CONSTRUCTS_INTERNAL_TSV = S4E9 / "constructs_l1_internal.tsv"
CONSTRUCT_RUNGS = ("cterm", "nterm", "multipass", "type3", "internal")
PROPEPTIDE_CHECK_TSV = REVIEWER / "nterm" / "propeptide_check.tsv"
UNIPROT_FEATURES_TSV = REVIEWER / "nterm" / "uniprot_features.tsv"
MASK_RESIDUES_TSV = S4E9 / "membrane-side-gate" / "l1_extracellular_mask_residues.tsv"
UNMASKED_AUDIT_TSV = S4E9 / "membrane-side-gate" / "unmasked_multipass_audit_2026_11_24.tsv"
TOPO_VS_EXP_SITES_TSV = S4E9 / "topology_vs_experiment_sites.tsv"
FOLD_ORDER_CHECK_TSV = S4E9 / "fold_order_check.tsv"
TAG_REFOLD_STRUCTURES_TSV = S4E9 / "tag_refold_structures.tsv"

GPI_OMEGA_CLEARANCE_AA = 12
SAME_LOOP_MIN_SPACING_AA = 25

#: Known UniProt canonical GPI omega sites where NetGPI reference omega differs.
KNOWN_UNIPROT_GPI_OMEGA: dict[str, int] = {
    "MDGA2": 931,
    "LYPD6": 147,
}

#: Known constitutively cleaved N-terminal propeptides of internal-tagged 0/1-TM proteins.
PINNED_CLEAVED_NTERM_PROPEPTIDES: dict[str, tuple[int, int]] = {
    "ATRN": (29, 83),
    "CDH3": (25, 107),
    "CDH5": (26, 47),
    "BACE1": (22, 45),
}

#: Propeptides removed only on extracellular activation or retained via disulfides.
NON_CONSTITUTIVE_PROPEPTIDE_GENES = frozenset(
    {"F2R", "F2RL1", "F2RL2", "F2RL3", "PRSS8", "PRSS21", "PRSS41", "MEP1B", "RGMA"}
)

#: Intentional post-ADAM17 juxtamembrane stalk site (not inside the N-terminal propeptide).
ALLOWED_PROPEPTIDE_EXEMPT_ARMS: frozenset[tuple[str, int]] = frozenset({("AREG", 193)})

#: Known failing cytoplasmic/membrane-core sites from the 2026-11-24 unmasked multipass audit.
KNOWN_UNMASKED_AUDIT_FAILS: dict[tuple[str, int], str] = {
    ("l1-p0113", 685): "ANO4 fail_cytoplasmic (z=-37.89 A)",
}

VALID_REFOLD_DISPOSITIONS = frozenset({"pass", "fail_flagged", "pending"})
FORBIDDEN_MASK_STATES = frozenset({"m", "c"})
LOOP_BOUNDARY_STATES = frozenset({"M", ".", "c"})


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _built_internal_loop_arms() -> list[dict[str, str]]:
    return [
        row
        for row in _read_tsv(DECISIONS_TSV)
        if row["decision"] == "build" and row["arm_role"].startswith("loop")
    ]


def _pos_int(value: str) -> int:
    return int(float(value))


def _load_mask_by_sha() -> dict[str, dict[str, str]]:
    by_sha: dict[str, dict[str, str]] = {}
    for row in _read_tsv(MASK_RESIDUES_TSV):
        sha = row["sequence_sha256"]
        if sha and sha not in by_sha:
            by_sha[sha] = row
    return by_sha


def _parse_propeptide_check_spans() -> dict[str, list[tuple[int, int]]]:
    spans: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for row in _read_tsv(PROPEPTIDE_CHECK_TSV):
        if row["propeptide_removal"] != "constitutive" or not row["n_terminal_propeptides"]:
            continue
        for token in row["n_terminal_propeptides"].split(";"):
            match = re.match(r"^(\d+)-(\d+):(\w+)$", token.strip())
            if match and match.group(3) != "unmapped":
                spans[row["design_id"]].append((int(match.group(1)), int(match.group(2))))
    return spans


def _parse_local_uniprot_nterm_propeptides() -> dict[str, list[tuple[int, int, int]]]:
    if not UNIPROT_FEATURES_TSV.is_file():
        return {}
    spans_by_gene: dict[str, list[tuple[int, int, int]]] = defaultdict(list)
    for row in _read_tsv(UNIPROT_FEATURES_TSV):
        propep_text = row.get("Propeptide", "")
        if not propep_text:
            continue
        sig_spans = [
            (int(a), int(b))
            for a, b in re.findall(r"SIGNAL (\d+)\.\.(\d+)", row.get("Signal peptide", ""))
        ]
        if not sig_spans:
            continue
        sp_end = sig_spans[0][1]
        seq_len = len(row.get("Sequence", ""))
        genes = [g.strip() for g in row.get("Gene Names (primary)", "").split(";") if g.strip()]
        for a_str, b_str in re.findall(r"PROPEP (\d+)\.\.(\d+)", propep_text):
            start, end = int(a_str), int(b_str)
            if start > sp_end + 6:
                continue
            for gene in genes:
                if gene not in NON_CONSTITUTIVE_PROPEPTIDE_GENES:
                    spans_by_gene[gene].append((start, end, seq_len))
    return spans_by_gene


def test_zero_gpi_omega_proximal_internal_loop_insertions() -> None:
    """Every built internal loop arm on a GPI gene stays >= 12 aa upstream of omega."""
    gpi_ref_omega: dict[str, int] = {}
    for row in _read_tsv(GPI_SET_TSV):
        if row["gpi_anchored"] == "1" and row["omega_site_reference"]:
            gpi_ref_omega[row["gene_symbol"]] = _pos_int(row["omega_site_reference"])

    checked = 0
    violations: list[str] = []
    for row in _built_internal_loop_arms():
        gene = row["gene"]
        omegas: list[tuple[str, int]] = []
        if gene in gpi_ref_omega:
            omegas.append(("gpi_set_ref", gpi_ref_omega[gene]))
        if gene in KNOWN_UNIPROT_GPI_OMEGA:
            omegas.append(("uniprot_pinned", KNOWN_UNIPROT_GPI_OMEGA[gene]))
        if not omegas:
            continue
        checked += 1
        pos = _pos_int(row["insert_after_residue"])
        min_omega = min(w for _, w in omegas)
        max_allowed = min_omega - GPI_OMEGA_CLEARANCE_AA
        if pos > max_allowed:
            violations.append(
                f"{row['arm_id']} ({gene} @{pos}): min_omega={min_omega} {omegas}, "
                f"clearance={min_omega - pos} aa < {GPI_OMEGA_CLEARANCE_AA} aa"
            )

    assert checked > 0, "expected at least one built internal loop arm on a GPI gene"
    assert not violations, "GPI omega-proximal internal loop arms found:\n  " + "\n  ".join(
        violations
    )


def test_zero_cleaved_propeptide_internal_loop_insertions() -> None:
    """No built internal loop arm sits inside a constitutively cleaved N-terminal propeptide."""
    pc_spans = _parse_propeptide_check_spans()
    uf_spans = _parse_local_uniprot_nterm_propeptides()
    native_len_by_arm = {
        row["construct_id"]: int(row["native_aa"]) for row in _read_tsv(CONSTRUCTS_INTERNAL_TSV)
    }

    violations: list[str] = []
    for row in _built_internal_loop_arms():
        arm_id = row["arm_id"]
        design_id = row["design_id"]
        gene = row["gene"]
        pos = _pos_int(row["insert_after_residue"])
        if (gene, pos) in ALLOWED_PROPEPTIDE_EXEMPT_ARMS:
            continue

        candidate_spans: list[tuple[str, int, int]] = []
        for start, end in pc_spans.get(design_id, []):
            candidate_spans.append(("propeptide_check.tsv", start, end))
        if gene in PINNED_CLEAVED_NTERM_PROPEPTIDES:
            start, end = PINNED_CLEAVED_NTERM_PROPEPTIDES[gene]
            candidate_spans.append(("pinned", start, end))
        arm_len = native_len_by_arm.get(arm_id)
        for start, end, seq_len in uf_spans.get(gene, []):
            if arm_len is None or arm_len == seq_len:
                candidate_spans.append(("uniprot_features.tsv", start, end))

        for source, start, end in candidate_spans:
            if start - 1 <= pos < end:
                violations.append(
                    f"{arm_id} ({gene} @{pos}) inside cleaved N-term propeptide "
                    f"{start}..{end} ({source})"
                )
                break

    assert not violations, (
        "Built internal loop arms inside cleaved N-terminal propeptides:\n  "
        + "\n  ".join(violations)
    )


def test_zero_membrane_core_cytoplasmic_or_unflagged_unusable_mask_loop_arms() -> None:
    """Built loop arms in the membrane-side mask avoid 'm'/'c' and flag usable==0 rows."""
    mm_by_sha = _load_mask_by_sha()
    checked = 0
    violations: list[str] = []

    for row in _built_internal_loop_arms():
        sha = row["sequence_sha256"]
        if sha not in mm_by_sha:
            continue
        checked += 1
        mr = mm_by_sha[sha]
        usable = int(float(mr["usable"]))
        pos = _pos_int(row["insert_after_residue"])
        mask = mr["mask"]
        jmask = mr["junction_mask"]
        arm_id = row["arm_id"]
        gene = row["gene"]

        if usable == 1:
            if not (1 <= pos <= len(mask) and 1 <= pos <= len(jmask)):
                violations.append(f"{arm_id} ({gene} @{pos}): out of mask bounds ({len(jmask)})")
                continue
            mc, jc = mask[pos - 1], jmask[pos - 1]
            if jc in FORBIDDEN_MASK_STATES or mc in FORBIDDEN_MASK_STATES:
                violations.append(
                    f"{arm_id} ({gene} @{pos}, usable=1): forbidden mask='{mc}' jmask='{jc}'"
                )
        else:
            if not row["flag"].strip():
                violations.append(
                    f"{arm_id} ({gene} @{pos}, usable=0): missing flag on unusable mask row"
                )
            if jmask and 1 <= pos <= len(jmask) and jmask[pos - 1] in FORBIDDEN_MASK_STATES:
                violations.append(
                    f"{arm_id} ({gene} @{pos}, usable=0): forbidden jmask='{jmask[pos - 1]}'"
                )

    assert checked > 0, "expected built internal loop arms in l1_extracellular_mask_residues.tsv"
    assert not violations, "Membrane-side mask invariant violations:\n  " + "\n  ".join(violations)


def test_zero_unmasked_multipass_audit_cytoplasmic_or_membrane_core_loop_arms() -> None:
    """No built loop arm fails the unmasked multipass cytoplasmic/membrane-core audit."""
    audit_by_arm: dict[str, str] = {}
    audit_by_site: dict[tuple[str, int], str] = {}
    if UNMASKED_AUDIT_TSV.is_file():
        for row in _read_tsv(UNMASKED_AUDIT_TSV):
            verdict = row["verdict"]
            audit_by_arm[row["arm_id"]] = verdict
            audit_by_site[(row["design_id"], _pos_int(row["insert_after_residue"]))] = verdict

    bad_verdicts = {"fail_cytoplasmic", "fail_membrane_core"}
    violations: list[str] = []
    for row in _built_internal_loop_arms():
        arm_id = row["arm_id"]
        design_id = row["design_id"]
        gene = row["gene"]
        pos = _pos_int(row["insert_after_residue"])

        verdict = audit_by_arm.get(arm_id) or audit_by_site.get((design_id, pos))
        if verdict in bad_verdicts:
            violations.append(f"{arm_id} ({gene} @{pos}): unmasked audit verdict={verdict}")
        elif (design_id, pos) in KNOWN_UNMASKED_AUDIT_FAILS:
            violations.append(
                f"{arm_id} ({gene} @{pos}): {KNOWN_UNMASKED_AUDIT_FAILS[(design_id, pos)]}"
            )

    assert not violations, "Unmasked multipass audit failures in built loops:\n  " + "\n  ".join(
        violations
    )


def test_zero_same_loop_close_pairs_under_25aa() -> None:
    """Built loop arms on the same continuous extracellular loop are >= 25 aa apart."""
    import membrane_mask as mm_mod

    mm_by_sha = _load_mask_by_sha()
    arms_by_sha: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in _built_internal_loop_arms():
        arms_by_sha[row["sequence_sha256"]].append(row)

    checked_cds = 0
    violations: list[str] = []
    for sha, rows in sorted(arms_by_sha.items()):
        if len(rows) < 2 or sha not in mm_by_sha:
            continue
        mr = mm_by_sha[sha]
        if int(float(mr["usable"])) != 1:
            continue
        checked_cds += 1
        mask, jmask = mr["mask"], mr["junction_mask"]
        ordered = sorted(
            ((_pos_int(r["insert_after_residue"]), r["arm_id"], r["gene"]) for r in rows),
            key=lambda item: item[0],
        )
        positions = [p for p, _, _ in ordered]
        regions = mm_mod.masked_regions(positions, (mask, jmask))

        for i in range(len(ordered)):
            p1, a1, gene = ordered[i]
            for j in range(i + 1, len(ordered)):
                p2, a2, _ = ordered[j]
                if abs(p2 - p1) >= SAME_LOOP_MIN_SPACING_AA:
                    continue
                between = mask[p1 - 1 : p2]
                no_boundary = not any(ch in LOOP_BOUNDARY_STATES for ch in between)
                if regions[i] == regions[j] or no_boundary:
                    violations.append(
                        f"{gene}: {a1}@{p1} and {a2}@{p2} (sep={p2 - p1} aa < "
                        f"{SAME_LOOP_MIN_SPACING_AA} aa, same_region={regions[i] == regions[j]}, "
                        f"mask_between={between!r})"
                    )

    assert checked_cds > 0, "expected multi-arm usable CDS entries in membrane-side mask"
    assert not violations, "Same-loop close pairs (<25 aa) found:\n  " + "\n  ".join(violations)


def test_internal_arm_decisions_veto_class_and_refold_disposition_hygiene() -> None:
    """Every built loop arm has valid non-hard veto_class and enum refold_disposition."""
    built = _built_internal_loop_arms()
    assert built, "no built internal loop arms found"

    bad_veto: list[str] = []
    bad_disp: list[str] = []
    for row in built:
        veto = row["veto_class"].strip()
        if not veto or veto.lower() == "nan" or veto == "hard":
            bad_veto.append(f"{row['arm_id']} ({row['gene']}): veto_class={row['veto_class']!r}")
        disp = row["refold_disposition"].strip()
        if disp not in VALID_REFOLD_DISPOSITIONS:
            bad_disp.append(
                f"{row['arm_id']} ({row['gene']}): refold_disposition={row['refold_disposition']!r}"
            )

    assert not bad_veto, "Invalid veto_class on built loop arms:\n  " + "\n  ".join(bad_veto)
    assert not bad_disp, "Invalid refold_disposition on built loop arms:\n  " + "\n  ".join(
        bad_disp
    )


def test_topology_vs_experiment_internal_site_count_matches_built_loops() -> None:
    """topology_vs_experiment_sites.tsv internal count equals constructs_l1_internal.tsv count."""
    tve_internal = [
        row for row in _read_tsv(TOPO_VS_EXP_SITES_TSV) if row["placement"] == "internal"
    ]
    constructs_internal = _read_tsv(CONSTRUCTS_INTERNAL_TSV)
    assert len(tve_internal) == len(constructs_internal), (
        f"topology_vs_experiment_sites.tsv internal count ({len(tve_internal)}) != "
        f"constructs_l1_internal.tsv count ({len(constructs_internal)})"
    )


def test_rung_numbers_snorkel_refold_matches_built_arms() -> None:
    """rung_numbers.snorkel() pass + fail equals built_arms (guards rung == 'cterm' filter)."""
    import rung_numbers as rn

    sn = rn.snorkel()
    refold = sn["snorkel_refold"]
    built_arms = sn["built_arms"]
    assert built_arms > 0
    assert refold.get("pass", 0) + refold.get("fail", 0) == built_arms, (
        f"snorkel_refold pass+fail ({refold}) != built_arms ({built_arms})"
    )
    assert sum(refold.values()) == built_arms, (
        f"snorkel_refold total ({refold}) != built_arms ({built_arms})"
    )


def test_cterm_rung_zero_gpi_and_extracellular_topology() -> None:
    """Zero GPI-anchored proteins sit in constructs_l1_cterm.tsv, and every non-snorkel
    C-terminal arm has an extracellular C-terminus (no cytosolic-C-tail member built
    without a snorkel)."""
    import build_l1_cterm_arms as b_cterm

    cterm_rows = _read_tsv(CONSTRUCTS_CTERM_TSV)
    assert cterm_rows, "constructs_l1_cterm.tsv is empty"

    gpi_genes = {r["gene_symbol"] for r in _read_tsv(GPI_SET_TSV) if r["gpi_anchored"] == "1"}
    if GPI_ANCHORED_GENES_TSV.is_file():
        gpi_genes |= {r["gene_symbol"] for r in _read_tsv(GPI_ANCHORED_GENES_TSV)}
    assert gpi_genes, "expected GPI-anchored genes in gpi_set.tsv"

    tag_paths = {r["candidate_id"]: r for r in _read_tsv(TAG_PATHS_TSV)}
    review_rows = _read_tsv(CTERM_REVIEW_TSV)
    review_by_id = {r["design_id"]: r for r in review_rows}
    human_review_by_gene = {r["gene_symbol"]: r for r in review_rows if r["species"] == "human"}

    gpi_violations: list[str] = []
    topology_violations: list[str] = []
    non_snorkel_checked = 0
    snorkel_checked = 0

    for row in cterm_rows:
        did = row["design_id"]
        gene = row["gene_symbol"]
        placement = row["placement"]
        tp = tag_paths.get(did)

        if gene in gpi_genes:
            gpi_violations.append(f"{did} ({gene}, placement={placement!r}): in GPI-anchored set")
        if tp is None:
            topology_violations.append(f"{did} ({gene}): missing from l1_tag_paths.tsv")
            continue
        if tp["gpi_rerouted"] != "0":
            gpi_violations.append(
                f"{did} ({gene}): l1_tag_paths.tsv gpi_rerouted={tp['gpi_rerouted']!r}"
            )

        if "snorkel" in placement.lower():
            snorkel_checked += 1
            if tp["placement"] != "internal":
                topology_violations.append(
                    f"{did} ({gene}): snorkel arm has l1_tag_paths placement={tp['placement']!r}"
                )
            continue

        non_snorkel_checked += 1
        if tp["placement"] != "C-terminal":
            topology_violations.append(
                f"{did} ({gene}, placement={placement!r}): non-snorkel arm has "
                f"l1_tag_paths placement={tp['placement']!r} (expected 'C-terminal')"
            )
        rev = review_by_id.get(did)
        if rev is None:
            topology_violations.append(
                f"{did} ({gene}): non-snorkel arm missing from l1_cterm_review.tsv"
            )
            continue
        cout = b_cterm.cout_disposition(rev, human_review_by_gene.get(gene))
        if cout.action == b_cterm.HOLD:
            topology_violations.append(
                f"{did} ({gene}): non-snorkel arm fails COUT_RULE ({cout.clause}: {cout.reason})"
            )

    assert non_snorkel_checked > 0 and snorkel_checked > 0
    assert not gpi_violations, "GPI-anchored proteins found in constructs_l1_cterm.tsv:\n  " + (
        "\n  ".join(gpi_violations)
    )
    assert not topology_violations, (
        "C-terminal rung topology invariant violations:\n  " + "\n  ".join(topology_violations)
    )


def test_unmasked_multipass_internal_arms_have_complete_non_failing_audit() -> None:
    """100% of unmasked multipass internal arms are covered by unmasked_multipass_audit
    with a non-failing verdict, and usable==0 multipass arms are either audited or
    explicitly flagged."""
    tag_paths = {r["candidate_id"]: r for r in _read_tsv(TAG_PATHS_TSV)}
    mask_by_sha = _load_mask_by_sha()
    assert UNMASKED_AUDIT_TSV.is_file(), f"missing {UNMASKED_AUDIT_TSV.name}"

    audit_by_arm: dict[str, str] = {}
    audit_by_site: dict[tuple[str, int], str] = {}
    for row in _read_tsv(UNMASKED_AUDIT_TSV):
        verdict = row["verdict"].strip()
        audit_by_arm[row["arm_id"]] = verdict
        audit_by_site[(row["design_id"], _pos_int(row["insert_after_residue"]))] = verdict

    unmasked_mp_checked = 0
    unusable_mp_checked = 0
    missing_audit: list[str] = []
    failing_audit: list[str] = []

    for row in _read_tsv(CONSTRUCTS_INTERNAL_TSV):
        did = row["design_id"]
        arm_id = row["arm_id"]
        gene = row["gene_symbol"]
        pos = _pos_int(row["cassette_after_residue"])
        n_tm = int(float(tag_paths.get(did, {}).get("n_tm", 0) or 0))
        verdict = audit_by_arm.get(arm_id) or audit_by_site.get((did, pos))

        if verdict is not None and verdict.startswith("fail"):
            failing_audit.append(f"{arm_id} ({gene} @{pos}): audit verdict={verdict!r}")

        if n_tm < 2:
            continue

        mr = mask_by_sha.get(row["native_sha256"])
        if mr is None:
            unmasked_mp_checked += 1
            if verdict is None:
                missing_audit.append(
                    f"{arm_id} ({did} {gene} @{pos}, n_tm={n_tm}): absent from "
                    f"l1_extracellular_mask_residues.tsv and unmasked_multipass_audit"
                )
        elif int(float(mr["usable"])) == 0:
            unusable_mp_checked += 1
            if verdict is None and "usable=0" not in row["flag"]:
                missing_audit.append(
                    f"{arm_id} ({did} {gene} @{pos}, n_tm={n_tm}): usable==0 in mask and "
                    f"neither audited nor explicitly flagged (flag={row['flag']!r})"
                )

    assert unmasked_mp_checked > 0, "expected unmasked multipass built internal arms"
    assert unusable_mp_checked >= 0
    assert not missing_audit, "Unaudited unmasked/unusable multipass internal arms:\n  " + (
        "\n  ".join(missing_audit)
    )
    assert not failing_audit, "Built internal arms failing unmasked multipass audit:\n  " + (
        "\n  ".join(failing_audit)
    )


def test_fold_to_order_parity_across_all_five_rungs() -> None:
    """Every built arm across all 5 constructs_l1_*.tsv files appears in fold_order_check.tsv
    with status=='ok' and has its tagged fold recorded in tag_refold_structures.tsv
    (kind=='tagged')."""
    import l1_arm_ids as arm_ids

    decisions = arm_ids.load_decisions()
    foc_rows = _read_tsv(FOLD_ORDER_CHECK_TSV)
    foc_by_arm = {r["arm_id"]: r for r in foc_rows}
    assert len(foc_by_arm) == len(foc_rows), "duplicate arm_id in fold_order_check.tsv"

    trs_tagged_by_fold = {
        r["fold_id"]: r for r in _read_tsv(TAG_REFOLD_STRUCTURES_TSV) if r["kind"] == "tagged"
    }
    trs_tagged_shas = {
        r["tagged_sha256"] for r in trs_tagged_by_fold.values() if r["tagged_sha256"]
    }
    assert trs_tagged_by_fold and trs_tagged_shas

    built_arms: set[str] = set()
    missing_foc: list[str] = []
    bad_foc: list[str] = []
    missing_trs: list[str] = []
    sha_mismatch: list[str] = []

    for rung in CONSTRUCT_RUNGS:
        rows = _read_tsv(S4E9 / f"constructs_l1_{rung}.tsv")
        assert rows, f"constructs_l1_{rung}.tsv is empty"
        for row in rows:
            arm_id = arm_ids.arm_id_for(rung, row, row["design_id"], decisions)
            assert arm_id not in built_arms, f"duplicate built arm_id {arm_id!r}"
            built_arms.add(arm_id)

            fr = foc_by_arm.get(arm_id)
            if fr is None:
                missing_foc.append(f"{rung}:{arm_id}")
                continue
            if (
                fr["status"] != "ok"
                or fr["pdb_sha256_ok"] != "yes"
                or fr["model_residues_equal"] != "yes"
            ):
                bad_foc.append(
                    f"{rung}:{arm_id} (status={fr['status']!r}, pdb_ok={fr['pdb_sha256_ok']!r}, "
                    f"res_eq={fr['model_residues_equal']!r})"
                )

            tr = trs_tagged_by_fold.get(fr["fold_id"])
            if tr is None or tr["tagged_sha256"] not in trs_tagged_shas:
                missing_trs.append(f"{rung}:{arm_id} (fold_id={fr['fold_id']!r})")
                continue

            row_shas = {
                s
                for s in (
                    row.get("protein_sha256", ""),
                    row.get("mature_sha256", ""),
                    row.get("refold_tagged_sha256", ""),
                )
                if s
            }
            if fr["protein_sha256"] not in row_shas:
                sha_mismatch.append(
                    f"{rung}:{arm_id} fold_order_check protein_sha256={fr['protein_sha256'][:12]} "
                    f"not in construct row shas"
                )
            elif tr["tagged_sha256"] not in row_shas and not (
                fr["folded_form"] == "mature chain" and fr["n_terminal_peptide_not_folded"].strip()
            ):
                sha_mismatch.append(
                    f"{rung}:{arm_id} tagged_sha256={tr['tagged_sha256'][:12]} not in construct "
                    f"row shas (folded_form={fr['folded_form']!r})"
                )

    extra_foc = sorted(set(foc_by_arm) - built_arms)
    assert not missing_foc, f"Built arms missing from fold_order_check.tsv: {missing_foc[:10]}"
    assert not extra_foc, f"Extra arms in fold_order_check.tsv not in builds: {extra_foc[:10]}"
    assert not bad_foc, f"Non-ok rows in fold_order_check.tsv: {bad_foc[:10]}"
    assert not missing_trs, (
        f"Built arms missing kind=='tagged' row in tag_refold_structures.tsv: {missing_trs[:10]}"
    )
    assert not sha_mismatch, f"Fold-to-order SHA-256 mismatches: {sha_mismatch[:10]}"


def test_no_minimal_ectodomain_non_gpi_reference_genes() -> None:
    """Every active non-GPI L1 reference gene must clear the >=30 aa ectodomain gate.

    Guards against the 2026-11-25 regression where tier LL in final_library.py omitted the
    ``ecd_ok(r)`` gate and admitted 20 minimal-ectodomain (<30 aa) non-GPI genes.
    Runs strictly offline from committed tables so it executes in both repos.
    """
    gpi_genes = {
        r["gene_symbol"]
        for r in _read_tsv(GPI_SET_TSV)
        if r.get("gene_symbol") and r.get("gpi_anchored") == "True"
    }
    ecto_excluded = {
        r["gene"]
        for r in _read_tsv(REVIEWER / "gate_lists" / "09_ectodomain_excluded.tsv")
        if r.get("gene")
    }
    gate_outcomes = {
        r["gene"]: r for r in _read_tsv(S4E9 / "gene_gate_outcomes.tsv") if r.get("gene")
    }
    violations: list[str] = []
    for r in _read_tsv(TAG_PATHS_TSV):
        if r["candidate_set"] != "library_reference":
            continue
        gene = r["gene_symbol"]
        if gene in gpi_genes:
            continue
        cid = r["candidate_id"]
        if gene in ecto_excluded:
            violations.append(f"{cid} ({gene}): in 09_ectodomain_excluded.tsv")
        outcome = gate_outcomes.get(gene)
        if outcome and outcome.get("ectodomain") == "fail":
            violations.append(f"{cid} ({gene}): gene_gate_outcomes ectodomain=fail")
        if int(r["n_tm"]) == 1 and int(r["longest_ecto_run_aa"]) < 30:
            violations.append(
                f"{cid} ({gene}): 1-TM non-GPI longest_ecto_run_aa={r['longest_ecto_run_aa']} < 30"
            )
    assert not violations, (
        f"Active non-GPI L1 reference genes below 30 aa ectodomain gate ({len(violations)}): "
        f"{violations[:10]}"
    )
