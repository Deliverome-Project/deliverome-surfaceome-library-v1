"""The membrane-side mask: which L1 junctions sit on the extracellular face.

DeepTMHMM v1 stays the topology of record, and this is a separate layer on top
of it. The s4e9 membrane-side gate (PR #291,
`library/l1/orf-design/membrane-side-gate/`) writes one row
per L1 internal member, keyed on the ORDERED protein's `sequence_sha256`. Its
`mask` string has one character per residue; 'E' means v1 says O and the
member's own structure model puts the residue at least 14 Å above the membrane
mid-plane, the shallowest published displayed site.

A tag inserted after residue n is eligible only if `junction_mask[n-1]` is 'E':
the gate's own per-junction verdict (mean height of residues n and n+1 at least
14 Å, residue n v1 O), reproduced by a committed test on PR #291. The per-residue
`mask` only sets region boundaries. (The first version required mask[n-1] and
mask[n] both 'E', which is stricter than the gate: on the adopted panel it
refused 40 junctions the gate passes at 14.0-21.5 Å.)

Without the mask the ranker designed 45 of 1,074 shipped internal sites on the
cytoplasmic side or inside the bilayer. v1 merged an anoctamin helix pair and its
cytoplasmic loop into one O segment (7 sites), or placed whole short loops on the
wrong side (7). The other 31 sit in the membrane core: 10 in v1 O loops that
barely leave the membrane (at most a quarter of their residues reach 14 Å; ATP1A1-4,
AQP10, MIP, KCNJ14, TRPM2, SLC12A8, SLC26A11), and 21 at the base of a loop that
does reach outside, 2-8 residues short of its nearest outside residue. The gate's
`cause` column (PR #291, l1_internal_pick_changes.tsv) carries this per site.

Manual exclusion (r4): the H7-H8 extracellular loop holding the SYGQ alpha-beta
assembly motif is 'x' in both masks for ATP1A1-4 and ATP12A. Na,K-ATPase alpha needs
beta to leave the ER, D1 structures show beta contacts across that loop, and
alanine scanning of SYGQ breaks assembly (Colonna et al. 1997, PMID 9139681).
Decided by Becca 2026-11-11 as the narrowest exclusion that covers it; per-gene
attribution is in the gate's manual_exclusions.tsv. Partner interfaces in general
stay evidence-only. The filter needs no special case: anything not 'E' is dropped.

Decided by Becca, 2026-11-10 (requirement) and 2026-11-11 ("talk to the ranker
so it knows to not run without the mask"). Implemented in the ranker by the
ranker-pu-segment session (Claude Code), PR #245. The mask itself was built by
the l1-membrane-side-gate session, PR #291.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
MASK_KEY = (
    "sprint4/s4e9-orf-design-fg-panel/membrane-side-gate/2026-11-12-r8/"
    "l1_extracellular_mask_residues.tsv"
)
LOCAL = REPO / "library/l1/orf-design/membrane-side-gate" / "l1_extracellular_mask_residues.tsv"
ELIGIBLE = "E"
#: Mask characters that remain within one extracellular loop (eligible 'E' and shallow
#: review-band 'r' dips); regions split only at transmembrane/non-O ('.', 'M', 'i'),
#: membrane-core ('m'), cytoplasmic ('c'), unmodelled ('u'), or excluded ('x') boundaries.
LOOP_CONTINUOUS = frozenset(("E", "r"))


class MaskMissing(RuntimeError):
    """A target the ranker was asked to design has no usable membrane-side mask."""


def load() -> pd.DataFrame:
    if not LOCAL.is_file():
        from deliverome_analysis.core import s3

        LOCAL.parent.mkdir(parents=True, exist_ok=True)
        s3.processed.download(MASK_KEY, str(LOCAL))
    m = pd.read_csv(LOCAL, sep="\t")
    u = m["usable"] == 1  # an unusable row (single-pass, unconfirmed face) may carry no mask
    for col in ("mask", "junction_mask"):
        assert (m.loc[u, col].str.len() == m.loc[u, "length_aa"]).all(), (
            f"a usable {col} is not one character per residue"
        )
    return m


def accession_to_sha() -> dict[str, str]:
    """The per-residue frame accession of each L1 member -> its ordered protein's sha256."""
    t = pd.read_csv(HERE / "l1_member_targets.tsv", sep="\t")
    r = pd.read_csv(HERE / "l1_member_reconciliation.tsv", sep="\t")
    j = t.merge(r[["candidate_id", "sequence_sha256"]], on="candidate_id", how="left")
    out = dict(zip(j["uniprot_acc"], j["sequence_sha256"], strict=False))
    # the adopted panel (tag_site_model_v5.py) also maps members exported after the 566
    # and the 27 Compara 116 orthologs and the 11 start-trimmed members, keyed on member id
    # (tag_sites_v5_orthologs.py, l1_ortholog_member_targets.tsv)
    for name in (
        "l1_adopted_targets.tsv",
        "l1_ortholog_member_targets.tsv",
        "l1_start_trim_member_targets.tsv",
    ):
        adopted = HERE / name
        if adopted.is_file():
            a = pd.read_csv(adopted, sep="\t")
            out.update(zip(a["uniprot_acc"], a["sequence_sha256"], strict=False))
    return out


def masks_for(accessions, mask: pd.DataFrame | None = None) -> tuple[dict[str, str], pd.DataFrame]:
    """Mask string per accession, and the targets that have none or an unusable one.

    Never falls back to v1 O: a target listed in the second return value must not
    be designed.
    """
    mask = load() if mask is None else mask
    by_sha = mask.drop_duplicates("sequence_sha256").set_index("sequence_sha256")
    sha = accession_to_sha()
    ok, refused = {}, []
    for acc in accessions:
        s = sha.get(acc)
        if s is None or s not in by_sha.index:
            refused.append({"uniprot_acc": acc, "reason": "no mask row for this protein"})
            continue
        row = by_sha.loc[s]
        if int(row["usable"]) != 1:
            refused.append(
                {
                    "uniprot_acc": acc,
                    "gene": row["gene_symbol"],
                    "reason": f"mask unusable ({row['status']}; orientation_supported="
                    f"{row['orientation_supported']})",
                }
            )
            continue
        ok[acc] = (row["mask"], row["junction_mask"])
    return ok, pd.DataFrame(refused)


def eligible_after(masks: tuple[str, str], residues) -> np.ndarray:
    """True where the gate passes "insert after residue n": junction_mask[n-1] == 'E'."""
    _mask, junction = masks
    res = np.asarray(residues, int)
    return np.array([1 <= r <= len(junction) and junction[r - 1] == ELIGIBLE for r in res], bool)


def masked_regions(residues, masks: tuple[str, str] | None = None) -> np.ndarray:
    """Region index for eligible residues, with boundaries set by the per-residue mask.

    Runs are contiguous residues staying on the extracellular face (mask in {'E', 'r'} or
    an eligible junction), so a short review-band ('r') dip does not split a single
    extracellular loop into separate regions; regions split only at transmembrane ('.'/'M'),
    membrane-core ('m'), cytoplasmic ('c'), or excluded/unmodelled boundaries.
    """
    res = np.asarray(residues, int)
    if masks is not None:
        mask, junction = masks
        run_of, run = {}, -1
        for i in range(len(mask)):
            if mask[i] in LOOP_CONTINUOUS or junction[i] == ELIGIBLE:
                if i == 0 or (i not in run_of):  # position i (1-based) not in the previous run
                    run += 1
                run_of[i + 1] = run
        missing = set(res.tolist()) - set(run_of)
        assert not missing, f"eligible residues outside every region: {sorted(missing)[:5]}"
        ids = np.array([run_of[r] for r in res])
        # renumber 0..k-1 in sequence order over the runs that hold an eligible residue
        uniq = {v: k for k, v in enumerate(sorted(set(ids.tolist())))}
        return np.array([uniq[v] for v in ids], int)
    order = np.argsort(res)
    out = np.empty(len(res), int)
    idx, prev = -1, None
    for i in order:
        if prev is None or res[i] != prev + 1:
            idx += 1
        out[i] = idx
        prev = res[i]
    return out


def apply(
    scored: pd.DataFrame, mask: pd.DataFrame | None = None
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Keep only eligible junctions, recompute regions on them, and list refused targets.

    `scored` holds one row per extracellular (v1 O) residue with `uniprot_acc`,
    `residue` and `region_index`. The v1 region is kept as `v1_region_index`.
    """
    ok, refused = masks_for(sorted(set(scored["uniprot_acc"])), mask)
    parts = []
    for acc, g in scored[scored["uniprot_acc"].isin(ok)].groupby("uniprot_acc"):
        g = g[eligible_after(ok[acc], g["residue"])].copy()
        if g.empty:
            refused = pd.concat(
                [
                    refused,
                    pd.DataFrame(
                        [{"uniprot_acc": acc, "reason": "no eligible junction under the mask"}]
                    ),
                ],
                ignore_index=True,
            )
            continue
        g["v1_region_index"] = g["region_index"]
        g["region_index"] = masked_regions(g["residue"], ok[acc])
        parts.append(g)
    out = pd.concat(parts, ignore_index=True) if parts else scored.iloc[0:0].copy()
    return out, refused
