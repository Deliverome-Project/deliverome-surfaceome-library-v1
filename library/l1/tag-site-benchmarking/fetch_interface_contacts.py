#!/usr/bin/env python3
"""Extracellular interface residues from D1, split by what binds there.

Reads the active release of the surfaceome D1 contact tables
(`contact_observation`, `contact_active_release`) for every accession in the
published corpus, the L1 members and the screened proteins, drops the
contexts that cannot be extracellular (`non_extracellular:*`), and classifies
each observation:

    antibody        AACDB, IEDB, SAbDab, Thera-SAbDab, any `therapeutic` / `tool`
                    category, and PDB partners named as Fab / scFv / VHH /
                    nanobody / antibody. An epitope is exposed surface the
                    protein does not need for function, so this class is NOT
                    used as a hard negative; it is written so it can be checked.
    protein_partner the partner is another UniProt accession (PDB/PDBe), or the
                    category is `receptor_partner`: co-receptor and other
                    hetero-protein interfaces
    homomer         the partner is the protein's own accession
    ligand          everything else: small molecules and peptides (BioLiP, PDBe
                    chemical components, GPCRdb, IUPHAR) and `endogenous_large`
                    ligands
    nonhuman_partner / mhc_presentation  protein partners that are not an endogenous
                    interface (see refine_partner); not used as hard negatives

    interface_contacts.tsv   uniprot_acc, residue, class, source, partner, pdb
    interface_contacts.json  release, counts, classification rules

    uv run --locked python fetch_interface_contacts.py
"""

from __future__ import annotations

import json
import re
import sys
from collections import Counter
from datetime import UTC, datetime
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[3] / "src"))

from deliverome_analysis.core.env import load_env  # noqa: E402

load_env()

from score_l1_library import d1  # noqa: E402

OUT = HERE / "interface_contacts.tsv"
META = HERE / "interface_contacts.json"
ANTIBODY_SOURCES = {"AACDB", "IEDB", "SAbDab", "Thera-SAbDab"}
ANTIBODY_NAME = re.compile(r"\b(fab|scfv|vhh|nanobody|antibody|igg|mab|darpin)\b", re.I)
ACCESSION = re.compile(r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$")
# Contexts that cannot touch an extracellular residue. Everything else is kept --
# GPCR and transporter pockets are labelled `membrane_spanning_site` although
# they include extracellular-loop residues -- and which positions count is
# decided later, by intersecting with the DeepTMHMM v1 extracellular residues.
NOT_EXTRACELLULAR = ("non_extracellular",)


def classify(acc: str, source: str, o: dict) -> str:
    partner = str(o.get("partner") or "")
    cat = str(o.get("category") or "")
    if (
        source in ANTIBODY_SOURCES
        or cat in ("therapeutic", "tool")
        or ANTIBODY_NAME.search(partner)
    ):
        return "antibody"
    base = partner.split("-")[0]
    if ACCESSION.match(base):
        return "homomer" if base == acc.split("-")[0] else "protein_partner"
    if cat == "receptor_partner":
        return "protein_partner"
    return "ligand"


# MHC and T-cell receptor chains: the "contact" is a peptide of the protein presented
# by MHC, not an interface of the folded protein.
MHC = re.compile(
    r"^(1[ABC]\d+|DR[AB]\d*|DQ[AB]\d*|DP[AB]\d*|DM[AB]|DO[AB]|2B[A-Z0-9]+|B2MG|HLA\w*|TR[ABDG]\w*)_HUMAN$"
)


def refine_partner(cls: str, partner: str) -> str:
    """Split protein partners that are not an endogenous interface of the folded protein.

    nonhuman_partner   a viral, bacterial or other non-human protein (a pathogen
                       receptor contact, a toxin), not an interface the cell uses
    mhc_presentation   HLA / beta2-microglobulin / TCR chains holding a peptide of the
                       protein
    Neither is used as a hard negative.
    """
    if cls != "protein_partner":
        return cls
    from fetch_uniprot_site_features import fetch

    base = partner.split("-")[0]
    d = fetch(base) if ACCESSION.match(base) else None
    if not d:
        return cls
    if d.get("taxid") not in (9606, None):
        return "nonhuman_partner"
    if MHC.match(str(d.get("entryName", ""))):
        return "mhc_presentation"
    return cls


def accessions() -> list[str]:
    sf = pd.read_csv(HERE / "site_features.tsv", sep="\t", usecols=["uniprot_acc", "label"])
    accs = set(sf.loc[sf["label"] == 1, "uniprot_acc"])
    accs |= set(pd.read_csv(HERE / "l1_member_targets.tsv", sep="\t")["uniprot_acc"])
    accs |= {"P78348", "P22001", "P56373", "P63252", "P02786", "P05556", "P18084", "P30530"}
    accs |= {"Q8N131"}
    return sorted(a for a in accs if ACCESSION.match(a.split("-")[0]))


def main() -> int:
    rel = d1("SELECT release_id FROM contact_active_release")[0]["release_id"]
    accs = accessions()
    rows, n_obs = [], Counter()
    for i in range(0, len(accs), 40):
        batch = accs[i : i + 40]
        ph = ",".join("?" * len(batch))
        for r in d1(
            "SELECT uniprot_acc, source, context, observation_json FROM contact_observation "
            f"WHERE release_id = ? AND uniprot_acc IN ({ph});",
            [rel, *batch],
        ):
            if str(r["context"]).startswith(NOT_EXTRACELLULAR):
                continue
            o = json.loads(r["observation_json"] or "{}")
            cls = refine_partner(
                classify(r["uniprot_acc"], r["source"], o), str(o.get("partner") or "")
            )
            n_obs[cls] += 1
            for p in o.get("positions", []):
                if isinstance(p, (int, float)):
                    rows.append(
                        {
                            "uniprot_acc": r["uniprot_acc"],
                            "residue": int(p),
                            "class": cls,
                            "source": r["source"],
                            "context": r["context"],
                            "partner": str(o.get("partner") or "").replace("\t", " "),
                            "pdb": o.get("pdb") or "",
                        }
                    )
    df = pd.DataFrame(rows).drop_duplicates()
    df.to_csv(OUT, sep="\t", index=False)
    by = df.drop_duplicates(["uniprot_acc", "residue", "class"]).groupby("class")
    META.write_text(
        json.dumps(
            {
                "fetched_at": datetime.now(UTC).isoformat(),
                "d1_release": rel,
                "accessions_queried": len(accs),
                "contexts": "all except non_extracellular:*; positions are UniProt numbering of uniprot_acc",
                "observations_by_class": dict(n_obs),
                "residues_by_class": by.size().to_dict(),
                "proteins_by_class": by["uniprot_acc"].nunique().to_dict(),
                "antibody_rule": "AACDB/IEDB/SAbDab/Thera-SAbDab, category therapeutic or tool, "
                "or partner named Fab/scFv/VHH/nanobody/antibody/IgG/mAb/DARPin",
            },
            indent=2,
        )
        + "\n"
    )
    print(json.dumps(json.loads(META.read_text()), indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
