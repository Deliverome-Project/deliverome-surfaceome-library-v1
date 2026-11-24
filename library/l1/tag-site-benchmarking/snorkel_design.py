#!/usr/bin/env python3
"""Snorkel constructs for the internal-rung members decided as `snorkel`, at protein level.

A snorkel tag (Brown et al. 2013, PLoS One 8:e73255, PMID 24023844) displays a tag outside
the cell without touching the protein's own loops. A linker, a transmembrane domain and the
tag are appended to the protein's cytoplasmic C-terminus. From the paper (pSNKL-Q):

    L1   22 residues: an 11-residue glycine/serine-rich sequence, then the 11 residues that
         sit next to the membrane on the cytoplasmic side of mouse PDGFRB (556-566), moved
         to the TM's N-terminal side
    TMD  mouse PDGFRB 530-555, in reversed topology (the paper counts it as 25 residues)
    L2   GS + FLAG + GMQ; the paper reports a tag 13 residues from the TM stains better than
         one 2 residues away, and extracellular regions under 40 residues all performed alike
    tag  HA

The paper does not print the 11-residue glycine/serine sequence or say how it altered the
flanking residues. So this is a Brown-style design, not pSNKL-Q itself:

    L1_GS         GGSGGSGGSGG                     (ours: 11 G/S residues)
    L1_JUXTA      QKKPRYEIRWK                     (mouse PDGFRB 556-566, UniProt P05622)
    TMD           KVVVISAILALVVLTVISLIILIMLW      (mouse PDGFRB 530-555, UniProt P05622)
    L2            GSGGGGS                         (ours: the cassette's first Arg sits 13
                                                   residues past the TM, as the paper's
                                                   better-staining FLAG placement did)
    TAG           the library's full 51-aa cassette (s4e9 CASSETTE_AA) in place of HA

The positive-inside rule holds as built. The cytoplasmic flank (QKKPRYEIRWK, plus K530 at
the TM's inner end) is Lys/Arg-rich. The first positive residue on the outside is the
cassette's Arg, 13 residues from the TM.

Decided by the ranker-pu-segment session, delegated by Becca 2026-11-12 ("use your best
judgment"). The rejected alternative is the paper's HA tag and unpublished linker. The DNA
arm (codons, seams, checks) belongs with the C-terminal arm builder: a snorkel arm is a
C-terminal append.

    uv run --locked python snorkel_design.py       # writes l1_snorkel_designs/VERSION/, offline
    uv run --locked python snorkel_design.py --live-uniprot   # also re-checks P05622 live
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[3] / "src"))

from build_internal_tagged import CASSETTE_AA  # noqa: E402
from l1_internal_decisions import DECIDED_BY, IMPLEMENTED_BY, head_tables  # noqa: E402

PDGFRB_MOUSE = "P05622"
L1_GS = "GGSGGSGGSGG"
L1_JUXTA = "QKKPRYEIRWK"  # P05622 556-566
TMD = "KVVVISAILALVVLTVISLIILIMLW"  # P05622 530-555
L2 = "GSGGGGS"
MODULE = L1_GS + L1_JUXTA + TMD + L2 + CASSETTE_AA
#: Dated version. l1_snorkel_designs.tsv at the experiment root is pinned by main's s4e9
#: l1_snorkel.py (sha256 3ddc9336, 45 members) and is never rewritten; each changed list is
#: published under a new dated key and the C-terminal builder moves its pin
#: (s5e13_publish.PINNED refuses an in-place overwrite).
# 2026-11-14: member ids re-keyed to main; +GJA8, OTOP1, SLC23A2 (literature veto)
# 2026-11-14b: SLC22A24 l1-p1884 leaves (isoform-fix redesign gives it an internal site)
VERSION = "2026-11-14b"
OUT = HERE / "l1_snorkel_designs" / VERSION / "l1_snorkel_designs.tsv"
OUT_JSON = HERE / "l1_snorkel_designs" / VERSION / "l1_snorkel_designs.json"


#: P05622 as UniProt held it on 2026-11-14 (sequence version 2), committed so the check runs offline
PDGFRB_FASTA = HERE / "pdgfrb_mouse_P05622.fasta"


def check_against_uniprot(live: bool = False) -> str:
    """The PDGFRB segments must match the UniProt record they are cited to.

    Reads the committed copy of P05622 by default; `--live-uniprot` fetches the current record
    and also requires it to equal the committed copy.
    """
    saved = "".join(ln.strip() for ln in PDGFRB_FASTA.read_text().splitlines()[1:])
    seq = saved
    if live:
        import urllib.request

        req = urllib.request.Request(
            f"https://www.ebi.ac.uk/proteins/api/proteins/{PDGFRB_MOUSE}",
            headers={"Accept": "application/json", "User-Agent": "deliverome-analysis s5e13"},
        )
        seq = json.loads(urllib.request.urlopen(req, timeout=60).read())["sequence"]["sequence"]
        assert seq == saved, "UniProt P05622 has changed since the committed copy"
    assert seq[529:555] == TMD, "TMD is not P05622 530-555"
    assert seq[555:566] == L1_JUXTA, "juxtamembrane segment is not P05622 556-566"
    return seq


def main() -> int:
    check_against_uniprot(live="--live-uniprot" in sys.argv[1:])
    assert len(L1_GS + L1_JUXTA) == 22 and len(TMD) == 26 and len(CASSETTE_AA) == 51
    first_pos_out = next(i for i, c in enumerate(L2 + CASSETTE_AA, 1) if c in "KR")
    assert first_pos_out >= 13, first_pos_out
    dec = pd.read_csv(HERE / "l1_internal_decisions.tsv", sep="\t")
    sn = dec[dec["decision"] == "snorkel"]
    _internal, _ids, prot = head_tables()
    rows = []
    for r in sn.itertuples(index=False):
        native = prot[r.sequence_sha256]
        full = native + MODULE
        rows.append(
            {
                "member_id": r.member_id,
                "design_id": r.design_id,
                "gene": r.gene,
                "sequence_sha256": r.sequence_sha256,
                "native_length": len(native),
                "native_c_terminus": native[-10:],
                "module_aa": MODULE,
                "module_length": len(MODULE),
                "snorkel_protein": full,
                "snorkel_protein_sha256": hashlib.sha256(full.encode()).hexdigest(),
                "components": "L1_GS|L1_JUXTA(P05622 556-566)|TMD(P05622 530-555)|L2|CASSETTE_AA",
                "status": "decided",
                "reason": "Brown-style snorkel (PMID 24023844), adapted: the paper's TMD and "
                "juxtamembrane segment from mouse PDGFRB; our G/S linker and L2, since the paper "
                "does not print them; the library's 51-aa cassette in place of HA",
                "rejected_alternative": "pSNKL-Q as published (unpublished linker sequence; HA tag)",
                "decided_by": DECIDED_BY,
                "implemented_by": IMPLEMENTED_BY,
            }
        )
    out = pd.DataFrame(rows)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT, sep="\t", index=False)
    OUT_JSON.write_text(
        json.dumps(
            {
                "citation": "Brown et al. 2013, PLoS One 8:e73255, PMID 24023844, "
                "doi:10.1371/journal.pone.0073255",
                "module": {
                    "L1_GS": L1_GS,
                    "L1_JUXTA": L1_JUXTA,
                    "TMD": TMD,
                    "L2": L2,
                    "TAG": CASSETTE_AA,
                },
                "module_length": len(MODULE),
                "first_positive_residue_outside_tm": first_pos_out,
                "members": out["member_id"].tolist(),
            },
            indent=2,
        )
        + "\n"
    )
    print(out[["member_id", "design_id", "gene", "native_length", "module_length"]].to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
