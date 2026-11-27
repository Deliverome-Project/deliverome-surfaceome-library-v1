"""Paired N-terminal arms for the panel's weak GPCR expressors: leader-free vs Igk leader.

Each gene gets two constructs that differ by the Igk leader and nothing else, so the
comparison reads as "does a cleavable leader rescue surface display". The leader-free arm
is not merely a control: it **is** the form the 7,066-member library would build for a
GPCR, so it has to be validated whatever the leader does.

**Why these receptors.** Tedman's immunostaining screen (external://tedman-gpcr-surface-
screen/2025-09-19/, Canonical sheet, `Immunostaining Intensity`) ranks the five panel
GPCRs CNR1 129.6 < CCR5 302.2 < OR51E2 923.1 < ACKR3 2199.5 < AVPR2 6366.6, against a
screen median of 873.2. CCR5 is the second weakest. OR51E2 is the olfactory arm, the
library's otherwise unsampled class; it sits mid-screen, and because Tedman ran without
olfactory RTP chaperones that means it reaches the surface autonomously.

**Not a signal-peptide substitution.** None of the five carries a signal peptide --
UniProt annotates none, and independently none appears among the 68 constructs Tedman
suffixed `_SP_HA`, which are the receptors that do have one (PAR1-4, FSHR/LHCGR/TSHR,
LGR4/5/6, the adhesion class-B set) and whose HA therefore had to go after it. So the
leader is **prepended** here; there is nothing to replace. It supplies the initiator and
the receptor keeps its own start codon, leaving the mature chain as the native receptor
preceded by a single Asp.

**Every arm carries the 51-aa tandem** (Becca, 2026-10-18). The 19-aa ALFA-only cassette
is retired from the batch, so cassette size is constant across all 23 constructs and no
comparison in the panel confounds it with position or with a leader.

**Tedman already validated the POSITION on all five, so the open variable is SIZE, not
placement.** Tedman 2025 tagged every class A and olfactory receptor with a single HA
immediately after the start codon -- the same junction these arms use -- keeping the native
sequence otherwise, and read surface display by anti-HA on unfixed cells. All five panel
GPCRs are class A plain-HA constructs, so none carried an exogenous leader, and each still
reported a surface value: AVPR2 6,367, ACKR3 2,200, OR51E2 923, CCR5 302, CNR1 130 a.u.
ACKR3 in particular reached 2,200 with an N-terminal insertion present DESPITE its
C34-C287 tether. The `_SP_HA` suffix in that dataset marks the class B/C/F placement --
HA after the native signal-peptide cleavage site -- not an added SP; it tracks class
(26/30 of class B, 0/393 of olfactory).

**So the tether argument is weaker than it reads, and what remains is cassette size.**
CCR5 C20-C269, OR51E2 C4-C168 and ACKR3 C34-C287 each bond the short N-terminal ectodomain
to a distal partner, and a leaderless signal-anchor protein must post-translationally
translocate everything upstream of TM1 against the positive-inside rule -- the SLC20A2
hazard. But a tethered leaderless N-terminus demonstrably already accepted an HA tag
(roughly 9-11 residues with the start context; the paper does not state the exact length).
Our cassette is 51 aa, about five times that. The untested proposition is the LENGTH, not
whether the position tolerates an insertion at all.

**What absorbs the risk is the leader, which only two of the three have.** CCR5 and
OR51E2 each carry a leader-free/igk pair, so if the leader-free arm fails, its igk twin
says whether co-translational routing rescues it. ACKR3 is leader-free only, so a failure
there cannot be split between the tether and the cassette -- the one design asymmetry
left after capping the batch at 23. AVPR2 is also leader-free only but carries no
N-terminal cysteine at all, so it has no tether to confound.

**An internal ACKR3 arm was evaluated and rejected, on size rather than on gates**
(tag-sites session, 2026-10-18). Two findings, and the second is the one that decides it.

No ACKR3 ECL has a gate-passing site: ECL2's distal segment is strand 196-198, a 3-residue
coil 199-201, then helix 202-209 at pLDDT 84-92 -- the Q187 failure mode again, and the
~202 position this module briefly considered is DSSP H. ECL3 280-286 is the best stretch
in any ACKR3 loop (genuine coil, 18-25 A clearance, RSA up to 0.67) and fails only a
DSSP-loop floor of 8 by one residue.

But an internal arm is not a size hedge, it is a SMALLER HOST, and that is independent of
every gate verdict. DSSP-irregular run length by region: N-terminal domain O1-45 has 39
residues, ECL1 7, ECL2 6 then 3, ECL3 7. The N-terminus is the only region on this protein
with room for a 51-aa insert by a wide margin, so moving the cassette inward concentrates
the same length risk into a seventh of the space. If the open question is whether 51 aa is
too much for a tethered leaderless N-terminus -- which is what Tedman's HA result narrows
it to -- then a shorter cassette at the validated N-terminal junction is the hedge, not a
full-length cassette somewhere smaller. That is a wet-lab call.

**Topological run length is not usable space, and this matters beyond ACKR3.** Our
``l1_tag_paths`` records ACKR3's longest inter-TM ectodomain run as 29 aa (ECL2 181-209),
while that region's irregular run is 6 residues then 3. MIN_LOOP_LEN measures the
topological run; a DSSP-loop floor measures the coil. For multipass proteins the two can
differ several-fold, so the 29 aa figure overstates insertable space by roughly 5x here.
Any library-wide internal-site count resting on topological runs inherits that gap.

**Upstream retirement to carry:** that pipeline's rank-1 ACKR3 internal site, P35, is
retired. It reports 12.1 A clear because the gate's hazard set contains only the ANNOTATED
C117-C196; residue 35 sits 3.81 A from C287's SG and 4.54 A from C34's. Unlike the CD20
enclosure bug this is a missing INPUT rather than a missing concept -- the existing 10 A
veto would have rejected it had UniProt annotated C34-C287. Our tag path for ACKR3 is
``n_terminal`` at residue 1 and never consumed P35, so nothing here depended on it.

**The igk contrast is the experiment Tedman did not run.** Of 766 canonical constructs,
ZERO genes appear with both a plain and an `_SP_HA` form, so that screen contains no
within-gene test of adding an N-terminal leader, and it never added an exogenous leader to
a class A receptor at all. Comparing its 68 `_SP_HA` rows against the 698 plain ones is
confounded by class, since the suffix tracks class. Each plain/igk pair here is therefore
a controlled contrast that does not exist upstream.

Barcodes are deliberately not reserved. ``reserve_barcodes.py`` is a one-shot allocator
against a shared registry, so these arms carry a placeholder until the order is settled.

    uv run --locked python build_gpcr_rescue_arms.py
    uv run --locked python build_gpcr_rescue_arms.py --write
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import build_constructs as bc
from build_loop_inserts import LOOP_TAG_AA

from deliverome_analysis.cloning import plasmid_qc as qc
from deliverome_analysis.cloning import seqtools

HERE = Path(__file__).parent
SEED = 20261011

#: Murine Ig kappa-chain V-J2-C signal peptide, 20 aa / 60 nt ending GGT, plus the +1 Asp
#: every published Igk vector keeps. Canonical source: ``deliverome_analysis.plasmid_qc``.
LEADER = qc.IGK_LEADER_NT
LEADER_AA = qc.IGK_LEADER_AA
PLUS1 = qc.IGK_PLUS1_NT
PLUS1_AA = qc.IGK_PLUS1_AA

ARMS = {
    "CNR1": {
        "tedman_au": 129.6,
        "n_ecto_aa": 116,
        "n_disulfide": "none",
        "cassette": "tandem",
        "why": (
            "weakest expressor on the panel; the only one whose ectodomain takes the 51-aa "
            "library cassette; and the library already plans both arms for it "
            "(l1-p0448 library_reference, l1-p2364 expression_rescue)"
        ),
    },
    "CCR5": {
        "tedman_au": 302.2,
        "n_ecto_aa": 30,
        "n_disulfide": "C20-C269",
        "cassette": "tandem",
        "why": "second weakest expressor on the panel in Tedman's screen",
    },
    "OR51E2": {
        "tedman_au": 923.1,
        "n_ecto_aa": 23,
        "n_disulfide": "C4-C168",
        "cassette": "tandem",
        "why": "olfactory arm: the library's otherwise unsampled GPCR class",
    },
    "ACKR3": {
        "tedman_au": 2199.5,
        "n_ecto_aa": 45,
        "n_disulfide": "3 N-terminal Cys (21, 26, 34); pairing not resolved here",
        "cassette": "tandem",
        "why": (
            "Brenna's constitutive-internalisation case: ACKR3 recycles without an agonist, "
            "which no other arm in the batch covers, and it is the only one of her five "
            "GPCRs that had no construct at all. Carries the 51-aa tandem like every "
            "other arm, over a MEASURED tether argument that is worth keeping on the "
            "record because it predicts where this arm fails. UniProt annotates "
            "ACKR3 only C117-C196, and that from a PROSITE rule (ECO:0000255); C21, C26 "
            "and C34 carry no annotation at all, which reads naively as a free N-terminus. "
            "SG-SG distances in AF-P25106-F1-model_v6 say otherwise: **C34-C287 is 2.03 A, "
            "a formed disulfide** (pLDDT 80/89), the same distance the model gives the "
            "PROSITE-inferred C117-C196 bond, so the N-terminus IS tethered to a distal "
            "partner exactly as CCR5's C20-C269 and OR51E2's C4-C168 are. C21 and C26 sit "
            "6.08 A apart and 25 A from C34, so they are free, but at pLDDT 51 -- the "
            "disordered N-terminal tip -- that half is weakly supported. Measured "
            "2026-10-18 after the tag-location session flagged that absence of annotation "
            "here is uninvestigated rather than negative. The tether at C34 argued for "
            "capping this arm at 19 aa; that cap was dropped the same day to hold the "
            "cassette constant across the batch. So this is the arm most likely to fail "
            "for translocation load, and -- being leader-free only -- the one whose "
            "failure the batch cannot attribute. An igk twin is the follow-up that would."
        ),
    },
    "AVPR2": {
        "tedman_au": 6366.8,
        "n_ecto_aa": 38,
        "n_disulfide": "none",
        "cassette": "tandem",
        "why": (
            "the STRONGEST expressor on the panel and Brenna's high-expression/high-rate "
            "corner, which nothing else in the batch occupies. Takes the 51-aa library "
            "cassette for the same reason CNR1 does: its 38-aa N-terminal ectodomain "
            "carries NO cysteine, so it is a free extension rather than a chain tethered "
            "to a distal partner (contrast CCR5 C20-C269 and OR51E2 C4-C168, the tether "
            "that once argued for capping those two at 19 aa). Pairing it with CNR1 "
            "brackets the expression range on the one cassette the library will actually "
            "build. AVPR2 and CNR1 are the two arms carrying it over no tether at all, "
            "so they read the cassette itself rather than the cassette-plus-tether. "
            "It also carries LESS post-translational translocation load than CNR1 "
            "(38+51=89 residues upstream of TM1 against 116+51=167), so if the tandem "
            "fails on AVPR2 the cassette is the problem, not the burden."
        ),
    },
}

#: construct suffix -> whether it carries the leader
LAYOUTS = {"plain": False, "igk": True}

#: Which layouts each gene actually builds. Not every receptor earns an Igκ arm.
#:
#: Brenna's stated rationale for the leader is "signal peptide for low expressors", and
#: against her screen median of 873.2 a.u. only two of the five qualify: CNR1 at 130
#: (0.15x) and CCR5 at 302 (0.35x). AVPR2 at 6,367 is 7.3x the median -- the STRONGEST
#: expressor on the panel -- and ACKR3 at 2,200 is 2.5x, so a rescue arm on either tests
#: rescue on a receptor that does not need rescuing. Both are dropped (Becca, 2026-10-18).
#:
#: OR51E2 at 923 is at the median and is kept "just to see" (Becca) despite having no
#: rescue rationale either, and despite Lv 2016's olfactory rescue being an N-terminal
#: BRIL fusion rather than a signal peptide -- so this arm tests a different mechanism
#: from the one shown to work for ORs. Recorded so the result is read for what it is.
#:
#: The leader-free arm is built for every gene regardless: it is the form the 7,066-member
#: library itself would construct, so it has to be validated whatever the leader does.
BUILD_LAYOUTS = {
    "CNR1": ("plain", "igk"),
    "CCR5": ("plain", "igk"),
    "OR51E2": ("plain", "igk"),
    "AVPR2": ("plain",),
    "ACKR3": ("plain",),
}

#: One-off arms outside the gene x layout grid. Empty since 2026-10-18:
#: AVPR2_plain_alfa existed only to vary cassette size, and with every arm on the
#: 51-aa cassette there is no size to vary.
EXTRA_ARMS: dict[str, dict] = {}


def _l1_cds() -> dict[str, str]:
    """Gene -> its L1 designed CDS, from the library's own record.

    ``candidate_cds.json`` is keyed by the manifest key, so a library reference is
    ``lib:<GENE>``. That file is the authoritative record of what the library synthesises;
    this folder's ``orf_domesticated.fasta`` is a panel-local copy.
    """
    import json

    from deliverome_analysis.cloning import library_design as ld
    from deliverome_analysis.core.s3 import processed

    with open(processed.download(ld.CDS_KEY)) as fh:
        raw = json.load(fh)
    return {g: raw[f"lib:{g}"]["cds"] for g in ARMS if f"lib:{g}" in raw}


def arm_names() -> list[str]:
    return [f"{gene}_{arm}" for gene in ARMS for arm in BUILD_LAYOUTS[gene]] + list(EXTRA_ARMS)


def arm_flanks(
    cds: dict[str, str], flank: int = 150, cassette: str = "alfa_only"
) -> dict[str, tuple[str, str]]:
    """Each arm's sequence either side of its cassette, for the encoding draw.

    The cassette goes after the receptor's own Met1, so the 3' flank is the rest of the
    ORF and is shared by both arms of a gene; the 5' flank differs because only one of
    them has the leader in front.

    Filtered by ``cassette`` so that only arms sharing a peptide are drawn against each
    other. Every arm in the batch is now ``tandem``, so the filter is a single group
    today; it stays because the 19-aa ALFA peptide is still buildable and a revived arm
    must not be diversified against a different peptide's encodings.
    """
    out: dict[str, tuple[str, str]] = {}
    for gene, spec in ARMS.items():
        if spec["cassette"] != cassette:
            continue
        body = cds[gene][:-3]  # drop the native stop; the construct supplies its own
        for arm in BUILD_LAYOUTS[gene]:
            leader = LAYOUTS[arm]
            upstream = bc.KOZAK + (LEADER + PLUS1 if leader else "") + body[:3]
            out[f"{gene}_{arm}"] = (upstream[-flank:], body[3:][:flank])
    for name, spec in EXTRA_ARMS.items():
        if spec["cassette"] != cassette:
            continue
        body = cds[spec["gene"]][:-3]
        upstream = bc.KOZAK + (LEADER + PLUS1 if spec["leader"] else "") + body[:3]
        out[name] = (upstream[-flank:], body[3:][:flank])
    return out


def _reserved_barcodes() -> dict[str, dict[str, str]]:
    """construct -> its reserved slot, from ``barcodes_gpcr.tsv``.

    Only the arms selected for the order carry one. The rest keep the ``N`` placeholder,
    because ``reserve_barcodes.py`` is a one-shot allocator against a shared registry and
    an unordered arm must not consume an index.
    """
    import csv

    path = HERE / "barcodes_gpcr.tsv"
    if not path.exists():
        return {}
    with path.open() as fh:
        return {r["construct"]: r for r in csv.DictReader(fh, delimiter="\t")}


def build(
    gene: str,
    arm: str,
    cds: str,
    cassette_nt: str,
    cassette_aa: str,
    barcode: str | None = None,
) -> dict:
    """Leader + Asp + receptor-with-its-ATG + cassette + rest of receptor.

    ``build_constructs.assemble`` would place a cassette at offset 0 for a gene with no
    signal peptide -- before the start codon -- while its own note says "after Met1". That
    path has never run, because no signal-peptide-free gene is in its topology table. The
    offset is explicit here so this code does not inherit the discrepancy.
    """
    leader = LAYOUTS[arm]
    body = cds[:-3]
    native = seqtools.translate(cds).rstrip("*")
    prefix = (LEADER + PLUS1) if leader else ""
    coding = prefix + body[:3] + cassette_nt + body[3:]

    tamed, edits = seqtools.domesticate(coding, bc.ENZYMES, seed=SEED)
    expect = (LEADER_AA + PLUS1_AA if leader else "") + native[0] + cassette_aa + native[1:]
    got = seqtools.translate(tamed)
    if got != expect:
        raise RuntimeError(
            f"{gene}_{arm}: protein changed\n  want {expect[:70]}\n  got  {got[:70]}"
        )
    if leader and tamed[: len(LEADER)] != LEADER:
        raise RuntimeError(f"{gene}_{arm}: domestication edited the published leader")

    # A reserved arm gets its real barcode so the junction context is screened as built;
    # the placeholder would hide a site created by the barcode/primer seam.
    bc_nt = barcode or "N" * 20
    insert = bc.KOZAK + tamed + bc.STOP + bc.SPACER + bc.FWD_PRIMER + bc_nt + bc.REV_PRIMER + bc.CS1
    scan = insert if barcode else insert.replace("N" * 20, "")
    sites = seqtools.find_sites(scan, bc.ENZYMES)
    if sites:
        raise RuntimeError(f"{gene}_{arm}: {len(sites)} Type IIS site(s) remain: {sites}")

    spec = ARMS[gene]
    return {
        "construct": f"{gene}_{arm}",
        "gene": gene,
        "arm": "igk_leader" if leader else "leader_free",
        "leader_aa": LEADER_AA if leader else "",
        "plus1_aa": PLUS1_AA if leader else "",
        "tedman_au": spec["tedman_au"],
        "n_term_ecto_aa": spec["n_ecto_aa"],
        "n_term_disulfide": spec["n_disulfide"],
        "cassette": spec["cassette"],
        "cassette_aa_len": len(cassette_aa),
        "receptor_met1_kept": True,
        "signal_peptide_replaced": False,
        "coding_nt": len(tamed),
        "protein_aa": len(got),
        "domestication_edits": len(edits),
        "insert_nt_with_barcode_placeholder": len(insert),
        "barcode": barcode or "not reserved",
        "insert": insert,
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--write", action="store_true")
    args = ap.parse_args(argv)

    # The ORF comes from the L1 library's own designed CDS, not this folder's
    # orf_domesticated.fasta. A rescue arm has to encode the protein the library actually
    # builds, or the arm answers a question about a different proteoform -- the same
    # distinction that put 131 of 326 L1 structure tabs on a UniProt isoform rather than
    # the canonical. Asserted rather than assumed: every arm's L1 CDS must translate to
    # exactly what the panel FASTA holds, so switching the source cannot move a sequence.
    cds = _l1_cds()
    panel = bc.read_fasta(HERE / "orf_domesticated.fasta")
    missing = sorted(set(ARMS) - set(cds))
    if missing:
        raise SystemExit(f"no L1 designed CDS for: {missing}")
    reserved = _reserved_barcodes()
    for gene in sorted(ARMS):
        if gene in panel:
            a = seqtools.translate(cds[gene]).rstrip("*")
            b = seqtools.translate(panel[gene]).rstrip("*")
            if a != b:
                raise SystemExit(
                    f"{gene}: the L1 CDS and the panel ORF encode different proteins "
                    f"({len(a)} aa vs {len(b)} aa) — resolve which proteoform the arm is "
                    f"for before building it"
                )

    # Each cassette size is issued by the diversifier that owns its peptide, so an arm
    # reads its encoding from that file and nowhere else. Falling back to the other
    # peptide's file, or to the reference constant, would silently reintroduce a tract
    # every member shares.
    sources = {
        "alfa_only": ("loop_cassette_encodings.tsv", "construct", "loop_tag_nt", LOOP_TAG_AA),
        "tandem": ("cassette_encodings.tsv", "gene", "cassette_nt", bc.CASSETTE_AA),
    }
    encodings: dict[str, dict[str, str]] = {}
    for kind, (fname, key, col, _) in sources.items():
        path = HERE / fname
        encodings[kind] = (
            {r[key]: r[col] for r in csv.DictReader(path.open(), delimiter="\t")}
            if path.exists()
            else {}
        )

    rows = []
    for gene, spec in ARMS.items():
        kind = spec["cassette"]
        if kind not in sources:
            raise SystemExit(f"{gene}: unknown cassette {kind!r}")
        fname, _, _, cassette_aa = sources[kind]
        for arm in BUILD_LAYOUTS[gene]:
            name = f"{gene}_{arm}"
            if name not in encodings[kind]:
                raise SystemExit(
                    f"{name}: no {kind} cassette encoding issued — rerun the diversifier "
                    f"that writes {fname}"
                )
            slot = reserved.get(name)
            rows.append(
                build(
                    gene,
                    arm,
                    cds[gene],
                    encodings[kind][name],
                    cassette_aa,
                    barcode=slot["barcode"] if slot else None,
                )
            )

    for name, spec in EXTRA_ARMS.items():
        kind = spec["cassette"]
        if name not in encodings[kind]:
            raise SystemExit(f"{name}: no {kind} cassette encoding issued — rerun the diversifier")
        _, _, _, cassette_aa = sources[kind]
        slot = reserved.get(name)
        rows.append(
            build(
                spec["gene"],
                "plain" if not spec["leader"] else "igk",
                cds[spec["gene"]],
                encodings[kind][name],
                cassette_aa,
                barcode=slot["barcode"] if slot else None,
            )
            | {"construct": name}
        )

    for r in rows:
        print(
            f"{r['construct']:<14} {r['arm']:<12} au={r['tedman_au']:<7} "
            f"ecto={r['n_term_ecto_aa']:<4} cassette={r['cassette_aa_len']} aa  "
            f"protein {r['protein_aa']:>4} aa  insert {r['insert_nt_with_barcode_placeholder']:>5} nt  "
            f"edits={r['domestication_edits']}"
        )

    if args.write:
        cols = [k for k in rows[0] if k != "insert"]
        with (HERE / "constructs_gpcr_rescue.tsv").open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows([{k: r[k] for k in cols} for r in rows])
        (HERE / "constructs_gpcr_rescue.fasta").write_text(
            "".join(
                # `bc=` is what check_swap_risk.py keys on. Without it a reserved arm is
                # reported as NOT COVERED and its ORF-to-barcode distance is never
                # measured -- the gate would pass while saying nothing about it.
                f">{r['construct']}|{r['arm']}|{r['cassette']}|"
                f"{r['insert_nt_with_barcode_placeholder']}nt"
                + (f"|bc={r['barcode']}" if r["barcode"] != "not reserved" else "")
                + "\n"
                + "\n".join(r["insert"][i : i + 60] for i in range(0, len(r["insert"]), 60))
                + "\n"
                for r in rows
            )
        )
        print("\nwrote constructs_gpcr_rescue.tsv and constructs_gpcr_rescue.fasta")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
