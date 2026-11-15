"""Build the tagged inserts for the L1 expression-rescue members (``rescue:*`` in the manifest).

**Not ordered.** One arm per gene (Becca, 2026-11-13: "make a single choice for each L1"): all 25
rescue members are retired and each gene is carried by its library reference
(`rescue_one_arm.tsv`, `rescue_one_arm.py`). For the 23 genes on the multipass rung the reference
is now built Igk + tag in `constructs_l1_multipass.tsv`, the same protein as the arm below. This
builder and its tables stay as the record of what the retired members were, until the order
assembly and the insert audits stop reading them.

Each rescue member is the murine Igk leader + Asp + a GPCR's own ORF, paired with that gene's
library reference, so the pair measures the leader within one gene. Their cassette was deferred.
This builds each as Igk + tag, the multipass rung's design (`build_l1_multipass_arms.py`):

    Igk leader + Asp | the receptor's Met1 | full 51-aa cassette | receptor residues 2..end

The references are built tag-alone on the multipass rung, so each pair differs by the leader
only.

**Decisions**

| Decision | Decided by | Implemented by |
|---|---|---|
| Build the rescue members as Igk + tag | l1-reviewer-package (delegated by Becca 2026-11-12), 2026-11-13 | session l1-nterm-res1-multipass |
| Full 51-aa cassette, leading GGGGS kept | Becca, 2026-11-11 ("keep full") | session l1-nterm-type3, PR #282 |
| Each arm gets its own linker-only cassette encoding; the receptor DNA stays as designed (reverting a receptor recode) | Becca, 2026-11-13 | session l1-nterm-res1-multipass |
| ADGRD2: Igk replaces its 18-residue native signal peptide (Igk + Asp + cassette + mature receptor from residue 19) | l1-reviewer-package (delegated by Becca 2026-11-12), 2026-11-13, on this session's recommendation | session l1-nterm-res1-multipass |

The layout, the protein check, the single-stop check and the two-wrapper-Esp3I check (alone and
in the p0039 backbone) are the multipass builder's.

    uv run --locked python build_l1_rescue_arms.py
    uv run --locked python build_l1_rescue_arms.py --check
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from build_constructs import CASSETTE, CASSETTE_AA, CS1, FWD_PRIMER, KOZAK, REV_PRIMER, SPACER
from build_gpcr_rescue_arms import LEADER, LEADER_AA, PLUS1, PLUS1_AA
from build_l1_multipass_arms import (
    BARCODE_SLOT,
    CASSETTE_DECIDED,
    ENZYMES,
    assemble_candidates,
    read_table,
    vector_check,
)

from deliverome_analysis import library_design as ld
from deliverome_analysis import seqtools

HERE = Path(__file__).parent
MANIFEST = HERE / "candidate_manifest_tranches.tsv"
OUT_TSV = HERE / "constructs_l1_rescue.tsv"
OUT_FASTA = HERE / "constructs_l1_rescue.fasta"
OUT_MATURE = HERE / "rescue_arms_mature.fasta"

DECIDED = (
    'l1-reviewer-package (delegated by Becca 2026-11-12, "use your best judgment"), 2026-11-13'
)
IMPLEMENTED_BY = "Claude Code session l1-nterm-res1-multipass"
REASON = (
    "rescue member is Igk leader + Asp + the receptor ORF with its cassette deferred; built as "
    "Igk + tag so it pairs with its tag-alone reference by the leader only"
)
REJECTED = "ship untagged (invisible to SpyCatcher/ALFA readout)"

#: design_id -> (native signal-peptide length, decided_by, reason, rejected) for rescue members
#: whose receptor carries its own signal peptide: Igk replaces it, then the cassette, then the
#: mature receptor. Keeping the native signal behind the cassette would let its cleavage take
#: the tag.
SP_SWAP = {
    "l1-p2361": (
        18,
        'l1-reviewer-package (delegated by Becca 2026-11-12, "use your best judgment"), 2026-11-13',
        "ADGRD2's receptor carries an 18-residue native signal peptide (DeepTMHMM v1 on the "
        "reference l1-p0061); Igk replaces it, as in the established signal-peptide swap "
        "(p0048 VSIG10), so the cassette sits on the mature N-terminus. The reference keeps "
        "its native signal peptide on the N-after-SP rung, so the pair still contrasts the leader",
        "Igk + tag at Met1 (cassette in front of an internal signal sequence that could be "
        "cleaved with it); tag after the native signal peptide (the reference's own design)",
    ),
}


#: Linker diversification. Each rescue arm gets its own linker-only cassette encoding: the four
#: GGGGS blocks recoded, SpyTag003 and P-ALFA-P held on p003's validated codons, drawn by the
#: library's own diversifier with the settings of `diversify_cassettes.py`. The receptor DNA
#: stays exactly as designed.
#:
#: Decided by Becca, 2026-11-13: "sharing bases if that's the prot seq is fine I jsut want to
#: diversify the linker since these proteins are at low abundance in the pool". That reverses a
#: synonymous recode of the receptor (l1-reviewer-package's rule, built at 9f9599d3), which is
#: the rejected alternative. Consequence, recorded rather than hidden: the linker sits upstream
#: of the receptor, so each rescue arm still shares its whole receptor-to-barcode stretch with
#: its tag-alone reference, and the predicted barcode swap within each pair is unchanged.
LINKER_DECIDED = (
    "Becca, 2026-11-13 (\"sharing bases if that's the prot seq is fine I jsut want to diversify "
    'the linker since these proteins are at low abundance in the pool")'
)
LINKER_REJECTED = "synonymous recode of the receptor CDS (built at 9f9599d3, reverted)"
REGISTRY = HERE / "l1_cassette_encodings.tsv"
LINKER_SEED = 2026_11_13
FLANK = 150


def issue_linker_encodings(flanks: dict[str, tuple[str, str]]) -> dict[str, str]:
    """One linker-only cassette encoding per rescue arm, distinct from the library registry."""
    import random

    from diversify_cassettes import LINKER_ONLY_JUNCTION_IDENTITY, tag_codons

    plan = seqtools.diversify_library_linkers(
        flanks,
        CASSETTE_AA,
        kmer_policy="minimize",
        enzymes=ENZYMES,
        fixed_codons=tag_codons(),
        max_junction_identity=LINKER_ONLY_JUNCTION_IDENTITY,
        rng=random.Random(LINKER_SEED),
    )
    issued = {name: v.sequence for name, v in plan.variants.items()}
    registry = {r["cassette_nt"] for r in read_table(REGISTRY)} if REGISTRY.exists() else set()
    fixed = tag_codons()
    for name, dna in issued.items():
        assert seqtools.translate(dna) == CASSETTE_AA, f"{name}: encoding changes the cassette"
        assert all(dna[3 * i : 3 * i + 3] == c for i, c in fixed.items()), f"{name}: tag recoded"
        assert dna != CASSETTE and dna not in registry, f"{name}: encoding already issued"
    assert len(set(issued.values())) == len(issued), "two rescue arms share an encoding"
    return issued


def build() -> list[dict]:
    from deliverome_analysis.s3 import processed

    # All 25 are retired (one arm per gene, Becca 2026-11-13); the record still builds them.
    rescue = [r for r in read_table(MANIFEST) if r["set"] == "expression_rescue"]
    constructs = json.loads(Path(processed.download(ld.CONSTRUCTS_KEY)).read_text())
    head = LEADER + PLUS1
    assert seqtools.translate(CASSETTE) == CASSETTE_AA and len(CASSETTE_AA) == 51
    flanks: dict[str, tuple[str, str]] = {}
    for m in rescue:
        d, receptor = m["design_id"], constructs[m["key"]]["construct"][len(head) :]
        if d in SP_SWAP:
            f5, f3 = head, receptor[3 * SP_SWAP[d][0] :]
        else:
            f5, f3 = head + receptor[:3], receptor[3:]
        flanks[d] = (f5[-FLANK:], f3[:FLANK])
    encodings = issue_linker_encodings(flanks)

    rows = []
    for m in sorted(rescue, key=lambda r: r["design_id"]):
        d, key = m["design_id"], m["key"]
        orf = constructs[key]["construct"]
        assert orf.startswith(head), f"{d}: rescue ORF does not open with Igk leader + GAC"
        receptor = orf[len(head) :]
        assert receptor[:3] == "ATG", f"{d}: receptor ORF does not start with ATG"
        cassette = encodings[d]
        native = seqtools.translate(orf).rstrip("*")
        row = {
            "construct_id": f"{d}:igk_tag",
            "design_id": d,
            "candidate_id": m["candidate_id"],
            "key": key,
            "gene_symbol": key.split(":")[1],
            "candidate_set": "expression_rescue",
            "design": "igk_tag",
            "status": "decided",
            "native_aa": len(native),
            "protein_aa": "",
            "mature_aa": "",
            "coding_nt": "",
            "insert_nt": "",
            "cassette": "",
            "type_iis_sites": "",
            "decided_by": DECIDED,
            "cassette_decided_by": "",
            "implemented_by": IMPLEMENTED_BY,
            "reason": REASON,
            "rejected_alternative": REJECTED,
            "protein_sha256": "",
            "mature_sha256": "",
            # The rescue member's own designed protein (Igk leader + Asp + receptor), and the
            # library reference protein its receptor part is identical to.
            "native_sha256": hashlib.sha256(native.encode()).hexdigest(),
            "reference_sha256": hashlib.sha256(
                seqtools.translate(constructs["lib:" + key.split(":")[1]]["construct"])
                .rstrip("*")
                .encode()
            ).hexdigest(),
            "cassette_encoding": cassette,
            "cassette_encoding_decided_by": LINKER_DECIDED,
            "cassette_encoding_rejected": LINKER_REJECTED,
            "insert": "",
            "mature": "",
        }
        receptor_aa = native[len(LEADER_AA) + 1 :]
        if d in SP_SWAP:
            sp, decided_by, reason, rejected = SP_SWAP[d]
            row.update(
                construct_id=f"{d}:igk_sp_swap",
                design="igk_sp_swap",
                decided_by=decided_by,
                reason=reason,
                rejected_alternative=rejected,
            )
            coding = head + cassette + receptor[3 * sp :]
            mature = PLUS1_AA + CASSETTE_AA + receptor_aa[sp:]
        else:
            coding = head + receptor[:3] + cassette + receptor[3:]
            mature = PLUS1_AA + receptor_aa[0] + CASSETTE_AA + receptor_aa[1:]
        protein = LEADER_AA + mature
        got = seqtools.translate(coding)
        assert got == protein + "*" and got.count("*") == 1, f"{d}: protein mismatch"
        insert = (
            assemble_candidates.ESP3I_5
            + KOZAK
            + coding
            + SPACER
            + FWD_PRIMER
            + BARCODE_SLOT
            + REV_PRIMER
            + CS1
            + assemble_candidates.ESP3I_3
        )
        sites = seqtools.find_sites(insert, ENZYMES)
        assert len(sites) == 2 and all(s[0] == "Esp3I" for s in sites), f"{d}: {sites}"
        row.update(
            protein_aa=len(protein),
            mature_aa=len(mature),
            coding_nt=len(coding),
            insert_nt=len(insert),
            cassette="full 51 aa, linkers diversified, tags on p003 codons",
            type_iis_sites="Esp3I x2 (wrapper)",
            cassette_decided_by=CASSETTE_DECIDED,
            protein_sha256=hashlib.sha256(protein.encode()).hexdigest(),
            mature_sha256=hashlib.sha256(mature.encode()).hexdigest(),
            insert=insert,
            mature=mature,
        )
        rows.append(row)
    vector_check([r for r in rows if r["insert"]])
    return rows


def render(rows: list[dict]) -> tuple[str, str, str]:
    cols = [c for c in rows[0] if c not in {"insert", "mature"}]
    tsv = "\t".join(cols) + "\n"
    tsv += "".join("\t".join(str(r[c]) for c in cols) + "\n" for r in rows)
    built = [r for r in rows if r["insert"]]
    fasta = "".join(
        f">{r['construct_id']} {r['candidate_id']} {r['gene_symbol']} {r['status']}\n{r['insert']}\n"
        for r in built
    )
    mature = "".join(
        f">{r['construct_id']} {r['candidate_id']} {r['gene_symbol']} "
        f"rescue_protein_sha256={r['native_sha256']} reference_sha256={r['reference_sha256']}\n"
        f"{r['mature']}\n"
        for r in built
    )
    return tsv, fasta, mature


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--check", action="store_true", help="rebuild and compare with the record")
    args = ap.parse_args(argv)
    rows = build()
    outs = dict(zip((OUT_TSV, OUT_FASTA, OUT_MATURE), render(rows), strict=True))
    if args.check:
        same = all(p.read_text() == text for p, text in outs.items())
        print("record matches rebuild" if same else "record DIFFERS from rebuild")
        return 0 if same else 1
    for p, text in outs.items():
        p.write_text(text)
    built = sum(1 for r in rows if r["insert"])
    print(f"wrote {OUT_TSV.name}, {OUT_FASTA.name}, {OUT_MATURE.name}: {built} built, "
          f"{len(rows) - built} held")  # fmt: skip
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
