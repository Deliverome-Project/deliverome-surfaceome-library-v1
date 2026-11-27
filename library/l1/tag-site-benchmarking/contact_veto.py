"""Ligand contacts on each L1 member's current ordered protein, for the internal-site veto.

Becca, 2026-11-12: internal sites should combine the random-forest prediction with "hard
negatives like the literature agent ligand sites". `compare_internal_policies.py` chose v5
plus a pick-time veto: no insert after n with a ligand contact at n-2 ... n+3 (VETO_HALF
residues either side of the junction; see below). That kept every tolerated-site metric and
removed every contact hit.

Three sources, each moved onto the member's CURRENT ordered protein (by sequence_sha256):

    d1_ligand              D1 structure-derived extracellular ligand contacts
                           (`member_interface_contacts.tsv`, class ligand)
    transferred_ligand     contacts transferred from a relative at >= 50% alignment
                           identity (`member_transferred_ligand_contacts.tsv`,
                           passes_threshold); precision 0.68 on held-out proteins
    literature             contacts the literature agents found and verify_scaleup checked
                           against the cited paper (PR #255), already numbered on the
                           current ordered protein

The first two are numbered on the 2026-11-07 ordered proteins. Where a member's sequence
has changed since, each contact is carried over by global alignment (BLOSUM62) and kept
only where the aligned residue is identical.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
#: The veto window for a cassette inserted after residue n: any contact at residues n-2 ... n+3,
#: i.e. within 2 residues of either residue flanking the insertion (n and n+1). Intended: the
#: insert disturbs both sides of the junction. Note that the benchmark that chose the veto
#: (compare_internal_policies.py via ranker_hard_negative_layer.near) used |residue - contact|
#: <= 2, i.e. n-2 ... n+2, one residue shorter on the C-terminal side; on the 2026-11-15
#: decisions that difference moved one site (HAS3 41).
VETO_HALF = 2
LIGAND_CLASSES = {
    "small_ligand",
    "peptide_ligand",
    "protein_ligand",
}  # scores are smoothed over +/-2, as in ranker_hard_negative_layer
#: Built 2026-11-14 with PR #255's own verify_scaleup.py and build_panel_table.py: the 394-gene
#: verified set plus the 26 outputs added since, verified the same way. Coverage of the 581
#: internal members: 437 searched (154 with verified contacts, 283 with none found) and 144 not
#: searched. A member with no literature contact is "none found" or "not searched", never "none".
#: The build code is scratch/becca/lit-contact-agents at LIT_CODE_COMMIT, on main since PR #255
#: merged (main a22130a2, 2026-11-14), so pinned_refs checks it like any other pin. LIT_CODE_KEY
#: is a saved copy of the same tree (SHA-256 LIT_CODE_SHA256), checked file by file against the
#: commit on 2026-11-14 (1,046 files, all identical) by l1-reviewer-package; it was the record
#: while the commit was off main. The table itself is the dated, hashed LITERATURE_KEY below.
LIT_CODE_COMMIT = "b63244fa54ab1315dd80f27396ae1503a5d4467b"
LIT_CODE_KEY = (
    "sprint5/s5e13-tag-site-benchmarking/lit_ligand_contacts/2026-11-14/"
    "build_code/lit_contact_agents_b63244fa.tar.gz"
)
LIT_CODE_SHA256 = "e1c6a21429e516a59b8b685e375d70de1ef0657a0412e45b9578c74fae8a5134"
LITERATURE_KEY = (
    "sprint5/s5e13-tag-site-benchmarking/lit_ligand_contacts/2026-11-14/lit_ligand_contacts_l1.tsv"
)
LITERATURE_LOCAL = HERE / "data" / "lit_ligand_contacts_2026-11-14.tsv"
COVERAGE_KEY = (
    "sprint5/s5e13-tag-site-benchmarking/lit_ligand_contacts/2026-11-14/lit_ligand_coverage_l1.tsv"
)
COVERAGE_LOCAL = HERE / "data" / "lit_ligand_coverage_2026-11-14.tsv"


def _residue_map(src: str, dst: str) -> dict[int, int]:
    """1-based residue of `src` -> residue of `dst`, identical aligned residues only."""
    if src == dst:
        return {i + 1: i + 1 for i in range(len(src))}
    from Bio import Align
    from Bio.Align import substitution_matrices

    al = Align.PairwiseAligner(mode="global")
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    al.open_gap_score, al.extend_gap_score = -10, -0.5
    a = al.align(src, dst)[0]
    out = {}
    for (s0, s1), (d0, _d1) in zip(*a.aligned, strict=False):
        for k in range(s1 - s0):
            if src[s0 + k] == dst[d0 + k]:
                out[s0 + k + 1] = d0 + k + 1
    return out


def literature_table(
    download: bool = True, classes: set[str] | frozenset[str] = frozenset(LIGAND_CLASSES)
) -> pd.DataFrame | None:
    if not LITERATURE_LOCAL.is_file() and download:
        from deliverome_analysis.core import s3

        if not s3.processed.exists(LITERATURE_KEY):
            return None
        LITERATURE_LOCAL.parent.mkdir(parents=True, exist_ok=True)
        s3.processed.download(LITERATURE_KEY, str(LITERATURE_LOCAL))
    if not LITERATURE_LOCAL.is_file():
        return None
    lit = pd.read_csv(LITERATURE_LOCAL, sep="\t")
    lit = lit[lit["citation_verified"].astype(str).str.lower() == "true"]
    # the veto uses ligands only, as the chosen policy: partner and homomer interfaces were
    # benchmarked as hard negatives and rejected (README, "Interfaces as hard negatives");
    # contact_lookup(evidence=True) reads the other classes as evidence, never for the veto
    return lit[lit["contact_class"].isin(classes)]


#: how a member's D1 and transferred ligand contacts were looked up (`contact_lookup`)
STRUCTURE_BASIS = ("own", "routed_in", "via_human_member", "via_same_gene_member", "none")
#: sources the veto (n-2 ... n+3) acts on; everything else contact_lookup returns is evidence only
VETO_SOURCES = frozenset({"d1_ligand", "transferred_ligand", "literature"})
#: D1 interface classes carried as evidence (evidence=True), labelled by class
EVIDENCE_D1_CLASSES = {"protein_partner": "d1_protein_partner"}
CONTACT_COLUMNS = ["sequence_sha256", "residue", "source", "via_design_id"]


def contact_lookup(
    current: pd.DataFrame,
    prot: dict[str, str],
    download: bool = True,
    human_of: dict[str, list[str]] | None = None,
    routed_ids: set[str] | None = None,
    gene_of: dict[str, str] | None = None,
    evidence: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Ligand contacts on the current ordered proteins, and what was looked up per member.

    `current` has design_id and sequence_sha256 for every member to cover; `prot` maps
    sequence_sha256 to protein sequence and must include the 2026-11-07 proteins.
    `human_of` maps a non-human member's design_id to the design_ids of the human members
    of its gene. `routed_ids` are the design_ids routed_in_contacts.py covered
    (l1_internal_decisions.ROUTED_IN), so one with no contacts still counts as looked up.
    `gene_of` maps a member's design_id to its gene symbol, for the same-gene carry.
    `evidence=True` also returns D1 protein-partner and literature partner/homomer contacts,
    labelled by class (`d1_protein_partner`, `literature_<class>`); they are never vetoed
    (`VETO_SOURCES`), because the interface benchmark rejected them as negatives.

    D1 and transferred contacts come from one of three places, recorded in the lookup's
    `structure_basis` column:

        own               the member is in the 2026-11-07 tables (`l1_internal_members.tsv`)
        routed_in         routed_in_ligand_contacts.tsv (members routed in later)
        via_human_member  an ortholog adopted after 2026-11-07: each human member's contacts
                          are carried onto the ortholog's ordered protein by global alignment,
                          identical aligned residues only (`_residue_map`), as for a changed
                          human sequence. The human member that supplied each contact is
                          recorded in `via_design_id`.
        via_same_gene_member
                          a member whose design_id is not in the 2026-11-07 tables (a start
                          trim, a new design_id): contacts of every 2026-11-07 member of the
                          same gene, carried the same way
        none              no lookup was possible; the member must not be labelled as having
                          D1 or transferred contacts checked

    Before 2026-11-15 the adopted orthologs fell through to `none` silently, while their
    decision rows still said "d1_ligand, transferred_ligand".
    """
    old = pd.read_csv(HERE / "l1_internal_members.tsv", sep="\t")
    old_sha = dict(zip(old["candidate_id"], old["sequence_sha256"], strict=False))
    d1 = pd.read_csv(HERE / "member_interface_contacts.tsv", sep="\t")
    classes = {"ligand": "d1_ligand", **(EVIDENCE_D1_CLASSES if evidence else {})}
    d1 = d1[d1["class"].isin(classes)].assign(source=lambda x: x["class"].map(classes))
    tr = pd.read_csv(HERE / "member_transferred_ligand_contacts.tsv", sep="\t")
    tr = tr[tr["passes_threshold"].astype(bool)].assign(source="transferred_ligand")
    both = pd.concat(
        [
            d1[["candidate_id", "residue", "source"]],
            tr[["candidate_id", "residue", "source"]],
        ]
    )
    by_member = {did: g for did, g in both.groupby("candidate_id")}
    routed = HERE / "routed_in_ligand_contacts.tsv"
    rc = pd.read_csv(routed, sep="\t") if routed.is_file() else None
    routed_shas = set(rc["sequence_sha256"]) if rc is not None else set()
    human_of, gene_of = human_of or {}, gene_of or {}
    old_by_gene = old.groupby("gene_symbol")["candidate_id"].apply(sorted).to_dict()

    def carry(src_did: str, sha_now: str, via: str | None) -> list[dict]:
        g = by_member.get(src_did)
        if g is None:
            return []
        m = _residue_map(prot[old_sha[src_did]], prot[sha_now])
        return [
            {"sequence_sha256": sha_now, "residue": m[r], "source": s, "via_design_id": via}
            for r, s in zip(g["residue"].astype(int), g["source"], strict=False)
            if r in m
        ]

    rows, lookup = [], []
    for did, sha_now in zip(current["design_id"], current["sequence_sha256"], strict=False):
        basis, via = "none", []
        if sha_now in prot and old_sha.get(did) in prot:
            basis = "own"
            rows += carry(did, sha_now, None)
        elif sha_now in routed_shas or did in (routed_ids or set()):
            basis = "routed_in"
        elif sha_now in prot:
            via = [h for h in human_of.get(did, []) if old_sha.get(h) in prot]
            if via:
                basis = "via_human_member"
            else:
                via = [
                    h
                    for h in old_by_gene.get(gene_of.get(did), [])
                    if h != did and old_sha.get(h) in prot
                ]
                basis = "via_same_gene_member" if via else "none"
            for h in via:
                rows += carry(h, sha_now, h)
        lookup.append(
            {"design_id": did, "sequence_sha256": sha_now, "structure_basis": basis,
             "via_design_id": ";".join(via)}
        )  # fmt: skip
    if rc is not None:  # members routed into the internal rung later (routed_in_contacts.py)
        rc = rc[rc["sequence_sha256"].isin(set(current["sequence_sha256"]))]
        rows += [
            {"sequence_sha256": s, "residue": int(r), "source": src, "via_design_id": None}
            for s, r, src in zip(rc["sequence_sha256"], rc["residue"], rc["source"], strict=False)
        ]
    lit = literature_table(download)
    if lit is not None:
        keep = lit[lit["sequence_sha256"].isin(set(current["sequence_sha256"]))]
        rows += [
            {"sequence_sha256": s, "residue": int(r), "source": "literature", "via_design_id": None}
            for s, r in zip(keep["sequence_sha256"], keep["residue"], strict=False)
        ]
    ev = literature_table(download, classes={"homomer", "partner"}) if evidence else None
    if ev is not None:
        keep = ev[ev["sequence_sha256"].isin(set(current["sequence_sha256"]))]
        rows += [
            {"sequence_sha256": s, "residue": int(r), "source": f"literature_{c}",
             "via_design_id": None}
            for s, r, c in zip(keep["sequence_sha256"], keep["residue"], keep["contact_class"],
                               strict=False)
        ]  # fmt: skip
    searched = literature_searched(download)
    lk = pd.DataFrame(lookup)
    lk["literature_searched"] = (
        lk["sequence_sha256"].isin(searched) if searched is not None else False
    )
    contacts = pd.DataFrame(rows, columns=CONTACT_COLUMNS).drop_duplicates(
        ["sequence_sha256", "residue", "source"]
    )
    return contacts, lk


def sources_available(lookup_row) -> str:
    """The `contact_sources_available` label, from what `contact_lookup` actually did."""
    basis = lookup_row["structure_basis"]
    if basis == "none":
        s = "no d1_ligand or transferred_ligand lookup (not in the contact tables)"
    elif basis == "via_human_member":
        s = f"d1_ligand, transferred_ligand via human member {lookup_row['via_design_id']}"
    elif basis == "via_same_gene_member":
        s = f"d1_ligand, transferred_ligand via same-gene member {lookup_row['via_design_id']}"
    else:
        s = "d1_ligand, transferred_ligand"
    lit = bool(lookup_row["literature_searched"])
    return s + (", literature" if lit else " (literature not run for this gene)")


def member_contacts(
    current: pd.DataFrame,
    prot: dict[str, str],
    download: bool = True,
    human_of: dict[str, list[str]] | None = None,
    routed_ids: set[str] | None = None,
    gene_of: dict[str, str] | None = None,
) -> pd.DataFrame:
    """One row per (sequence_sha256, residue, source) on the current ordered proteins."""
    return contact_lookup(current, prot, download, human_of, routed_ids, gene_of)[0]


def literature_searched(download: bool = True) -> set[str] | None:
    """sequence_sha256 of members whose gene had a literature run (None: no coverage table)."""
    if not COVERAGE_LOCAL.is_file() and download:
        from deliverome_analysis.core import s3

        if not s3.processed.exists(COVERAGE_KEY):
            return None
        COVERAGE_LOCAL.parent.mkdir(parents=True, exist_ok=True)
        s3.processed.download(COVERAGE_KEY, str(COVERAGE_LOCAL))
    if not COVERAGE_LOCAL.is_file():
        return None
    c = pd.read_csv(COVERAGE_LOCAL, sep="\t")
    return set(c.loc[c["searched"].astype(str).str.lower() == "true", "sequence_sha256"])


def vetoed(contacts: pd.DataFrame, sha: str, n: int) -> list[str]:
    """Sources of any contact at residues n-VETO_HALF ... n+1+VETO_HALF (insertion after n)."""
    c = contacts[contacts["sequence_sha256"] == sha]
    near = c[(c["residue"] >= n - VETO_HALF) & (c["residue"] <= n + 1 + VETO_HALF)]
    return sorted(set(near["source"]))


#: contact groups shown per arm (evidence, not all of them veto): which D1 classes and literature
#: contact_class values fall in each
GROUPS = {
    "ligand": ({"ligand"}, LIGAND_CLASSES),
    "partner": ({"protein_partner"}, {"partner"}),
    "homomer": ({"homomer"}, {"homomer"}),
}


def classified_contacts(current: pd.DataFrame, prot: dict[str, str]) -> pd.DataFrame:
    """Every contact on the current ordered proteins, by group, with source and evidence.

    One row per (sequence_sha256, residue, group, source): D1 structure contacts (moved by
    identical alignment, as in member_contacts), transferred ligand contacts, literature rows
    (verified; partner and homomer classes are the agents' calls and are labelled so), and the
    routed-in / rescue contact tables. D1 holds no homomer contacts at all (no observation pairs
    a protein with its own accession), so homomer rows come only from the literature.
    """
    old = pd.read_csv(HERE / "l1_internal_members.tsv", sep="\t")
    old_sha = dict(zip(old["candidate_id"], old["sequence_sha256"], strict=False))
    cur = dict(zip(current["design_id"], current["sequence_sha256"], strict=False))
    d1 = pd.read_csv(HERE / "member_interface_contacts.tsv", sep="\t")
    tr = pd.read_csv(HERE / "member_transferred_ligand_contacts.tsv", sep="\t")
    tr = tr[tr["passes_threshold"].astype(bool)].assign(**{"class": "ligand"})
    rows = []
    for name, df, ev in (("d1", d1, "structure"), ("transferred", tr, "transferred")):
        for did, g in df.groupby("candidate_id"):
            sha_now, sha_then = cur.get(did), old_sha.get(did)
            if sha_now is None or sha_then not in prot or sha_now not in prot:
                continue
            m = _residue_map(prot[sha_then], prot[sha_now])
            for r, cls in zip(g["residue"].astype(int), g["class"], strict=False):
                grp = next((k for k, (d1c, _l) in GROUPS.items() if cls in d1c), None)
                if grp and r in m:
                    rows.append((sha_now, m[r], grp, f"{name}_{cls}", ev))
    lit = literature_table(download=False)
    if LITERATURE_LOCAL.is_file():
        full = pd.read_csv(LITERATURE_LOCAL, sep="\t")
        full = full[full["citation_verified"].astype(str).str.lower() == "true"]
        full = full[full["sequence_sha256"].isin(set(cur.values()))]
        for s, r, cls, ev in zip(
            full["sequence_sha256"],
            full["residue"],
            full["contact_class"],
            full["evidence"],
            strict=False,
        ):
            grp = next((k for k, (_d, lc) in GROUPS.items() if cls in lc), None)
            if grp:
                rows.append((s, int(r), grp, f"literature_{cls}", str(ev)))
    del lit
    for f in ("routed_in_ligand_contacts.tsv", "rescue_round_ligand_contacts.tsv"):
        p = HERE / f
        if p.is_file():
            t = pd.read_csv(p, sep="\t")
            for s, r, src in zip(t["sequence_sha256"], t["residue"], t["source"], strict=False):
                rows.append((s, int(r), "ligand", src, "structure" if src == "d1_ligand" else src))
    return pd.DataFrame(
        rows, columns=["sequence_sha256", "residue", "group", "source", "evidence"]
    ).drop_duplicates()


def nearest(contacts: pd.DataFrame, sha: str, n: int, group: str) -> dict:
    """Nearest contact of a group to the insertion after n.

    `nearest_<group>_distance` is residues from the junction (0 on either flanking residue);
    `nearest_<group>_residue` is that contact's residue, numbered on the member's ordered protein.
    (Before 2026-11-15 the distance was written as `nearest_<group>_residues`, which readers
    took both ways.)
    """
    c = contacts[(contacts["sequence_sha256"] == sha) & (contacts["group"] == group)]
    if c.empty:
        return {f"nearest_{group}_distance": None, f"nearest_{group}_residue": None,
                f"nearest_{group}_source": "", f"nearest_{group}_evidence": ""}  # fmt: skip
    d = (c["residue"] - n).where(c["residue"] <= n, c["residue"] - (n + 1)).abs()
    i = d.idxmin()
    return {
        f"nearest_{group}_distance": int(d[i]),
        f"nearest_{group}_residue": int(c.loc[i, "residue"]),
        f"nearest_{group}_source": c.loc[i, "source"],
        f"nearest_{group}_evidence": c.loc[i, "evidence"],
    }
