"""Loop-insertion constructs for the panel members with no extracellular terminus.

These arms carry the SAME p003 cassette as the terminal constructs
(GGGGS-SpyTag003-GGGGSGGGGS-P-ALFA-P-GGGGS, 51 aa), so every construct in the batch is
SpyCatcher-conjugatable and cassette size is constant across the panel. The reasoning,
and what it leaves unsettled for the library, is recorded at ``LOOP_TAG`` below.

The accessible-surfaceome tag_sites pipeline proposes ALFA or DogTag at 15 aa with "GS
both sides" for every surface_loop site it emits, and that 19-aa P-flanked ALFA v1.5
peptide is still built here as ``ALFA_P_FLANKED``; it is simply not what any arm in this
batch uses.

Each gene carries one site the deterministic pipeline proposed and one Becca selected to
probe a cutoff, so the pair reads as "the lane's pick versus the thing the lane refuses".
Pipeline sites are recorded in tag_sites_loop_insertions.json; the selected ones in
loop_site_overrides.json, which is a separate file because that JSON is a verbatim record
of pipeline output and hand-picked sites must not be attributed to it.

**SLC2A4 P66 is the lead site.** The pipeline ranks surface_loop candidates by
(-own RSA, conservation, -gap_freq) — pLDDT is only a >=70 pass/fail reliability gate,
never a ranking key. P66 wins on that key: own RSA 0.929, conservation 0.0. The rank does
not survive into the emitted JSON, which is re-sorted by residue and then
lexicographically by site_id, so file order is not a ranking.

The literature agrees. Two independent lineages put the tag in GLUT4's first exofacial
loop, where P66 sits: Kanai 1993 J Biol Chem (PMID 7686158), c-myc in the "first
ectodomain", which became the field-standard GLUT4myc reporter; and Quon 1994 PNAS
(PMID 8202531, doi:10.1073/pnas.91.12.5587), HA1 "in the first exofacial loop".

**The two comparison arms moved on 2026-10-13** (Becca), from E193 to G320 and from N171
to P79. Each replacement probes a different gate:

* SLC2A4 **G320** sits in the 309-323 extracellular region, 254 residues from P66, so the
  pair separates "this loop" from "this protein". E193 was the first residue of a 9-aa
  region, 1 residue from TM5, and its arm confounded the two.
* MS4A1 **P79** sits in the 6-residue EL1 (79-84), 1 residue from TM1 — a 51-aa cassette
  into a 6-aa loop, the sharpest test in the batch of the ``MIN_LOOP_LEN`` floor. If it
  fails it indicts the floor, not CD20.
* Dropping N171 also removes a hazard nobody had flagged: CD20 carries a **disulfide bond
  167-183** (UniProt P11836), so inserting after 171 put the cassette *inside* a
  disulfide-closed loop. Y161 (dist 20 to the nearest TM) and P79 are both outside that bracket.

Each variant is its own construct and gets its own barcode.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from build_constructs import (
    CASSETTE,
    CASSETTE_AA,
    CS1,
    ENZYMES,
    FWD_PRIMER,
    KOZAK,
    REV_PRIMER,
    SPACER,
    STOP,
    read_fasta,
)

from deliverome_analysis.cloning import seqtools

HERE = Path(__file__).parent

# P-flanked ALFA v1.5, p003 codons, plus a GS linker each side -> 19 aa.
ALFA_P_FLANKED = "CCCAGTCGCCTCGAAGAAGAGCTCAGGCGGAGACTGACAGAGCCT"  # P-SRLEEELRRRLTE-P
GS = "GGCAGC"  # G-S
# The loop arms now carry the SAME 51-aa cassette as the terminal constructs
# (Becca, 2026-10-18). They were 19-aa ALFA-only on the argument that a 51-aa insert does
# not belong in a short extracellular loop -- but that argument conflates "fits inside the
# loop" with "is tolerated at the insertion point". An insertion ADDS length; it does not
# have to be accommodated within the existing span, and the cassette's GGGGS linkers are
# flexible, so loop length is not obviously the binding constraint.
#
# What this buys: all 23 constructs now carry SpyTag003, so every one is SpyCatcher-
# conjugatable. The 19-aa cassette has ALFA only, which silently excluded 6 internal arms
# from that chemistry -- and the README had "decide whether internal-rung members carry
# SpyTag003" as an OPEN question, so 19 aa was a placeholder nobody had ratified.
#
# The library-wide consequence is real and is NOT settled by this: MIN_LOOP_LEN was
# calibrated for a 19-aa cargo. At 51 aa the floor has to be recalibrated, and the README
# notes that the 51-aa tandem does not fit SLC16A9/12/14's 16-residue runs -- so a whole
# paralog block's internal-rung membership turns on it. This change is for the panel.
LOOP_TAG = CASSETTE
LOOP_TAG_AA = CASSETTE_AA

# gene -> variant -> (insert_after_residue, provenance)
VARIANTS = {
    "SLC2A4": {
        "V319": (
            319,
            "selected (Becca, 2026-10-15); tag-site model rank 2, within-protein 0.570 / "
            "absolute 0.562; extracellular 309-322 (14 aa, v1), 10/3 residues of room, "
            "4 from TM8 — separate-region arm, moved 1 residue off G320",
        ),
        "E67": (
            67,
            "selected (Becca, 2026-10-15); tag-site model rank 1, within-protein 0.748 / "
            "absolute 0.553; exofacial loop 1 (46-81, 36 aa), 21/14 residues of room, "
            "15 from TM2 — the GLUT4myc literature loop (Kanai 1993, PMID 7686158), "
            "moved 1 residue off P66",
        ),
    },
    "AVPR2": {
        "L302": (
            302,
            "selected (Becca, 2026-10-18) to answer what a LOOP insertion costs in display "
            "relative to a free terminus. AVPR2 is the panel member that can: it is the "
            "strongest expressor (6,367 a.u. in Tedman, against a screen median of 873), so "
            "a loss is measurable against a high baseline, and it already carries a "
            "terminal arm (p0047) on the same receptor. ECL3 (293-306, 14 aa, v1) between "
            "TM6 (272-292) and TM7 (307-328); 9/4 residues of room, 5 from TM7, and no "
            "cysteine within 10 residues, so no disulfide is at risk. Terminal p0047 and "
            "this loop arm carry the SAME 51-aa cassette, so the pair varies position "
            "alone. That is what putting every arm on the full cassette bought: the "
            "earlier design paired a 51-aa terminus against a 19-aa loop and so varied "
            "position and cassette size together, which is why a third AVPR2 arm was "
            "drafted to collapse the confound. With the cassette held constant the "
            "confound is gone and that arm is unnecessary. CCR5 repeats the same "
            "position contrast on a weak expressor (302 a.u.; terminal p0042 against its "
            "S271 loop arm), so the two bracket the question at both ends of the range",
        ),
    },
    "MS4A1": {
        "F146": (
            146,
            "selected 2026-10-18 after the tag-sites session showed Q187 and 182 are both "
            "inside a confidently modelled HELIX (AF-P11836-F1-model_v6 residues 178-190 "
            "are a contiguous alpha-helical dihedral run at pLDDT 84-93, confirmed here "
            "independently of their pydssp), and 5.6 A from C167. Topology alone could not "
            "see that: 187 is outside the 167-183 disulfide SPAN yet 5.6 A from its "
            "cysteine in space. 146 is the best remaining CD20 junction -- a genuine "
            "7-residue coil in the N-terminal half of the 142-188 run, RSA 0.73, 16.0 A "
            "clear of both cysteines, the largest clearance of any candidate measured. "
            "It fails ONE criterion, and softly: the tag-sites pipeline floors a DSSP loop "
            "at 8 residues and this is 7. Note that is a different gate from our "
            "MIN_LOOP_LEN, which measures the TOPOLOGICAL run (47 aa here) and passes it "
            "easily -- the two floors are not interchangeable.\n\n"
            "CHALLENGED AND UPHELD 2026-10-18, then the challenge was WITHDRAWN at source. "
            "A third read called 142-150 a helix, which would have given this arm Q187's "
            "defect. Its provenance turned out to be the answer: that column is `ss3`, "
            "netsurfp-3.0 -- a SEQUENCE predictor that never sees the structure -- so "
            "against structure-based reads it had no standing. A fourth read, "
            "Kabsch-Sander hydrogen-bond DSSP run on the model, segments it "
            "140-142 H / 143-152 C / 153-156 H / 157-158 G / 159-171 C / 172-174 G / "
            "175-177 C / 178-192 H, reproducing pydssp's 159-171 boundary exactly. 146 is "
            "coil on three of four reads and the dissenter is the one that never saw a "
            "model. It is also the best CD20 junction on both axes: absolute 0.6322, which "
            "is the region ceiling, and 14.7 A heavy / 16.0 A CA clear of both cysteines.\n\n"
            "The two structure-based methods that upheld it before that withdrawal: "
            "per-residue dihedral basins across 143-149 are other/beta/ALPHA/beta/beta/"
            "beta -- mixed irregular, mostly beta and other, not an alpha run -- and "
            "pydssp finds NO regular hydrogen-bond pattern anywhere in 143-149, which is "
            "what actually makes a helix. A separate smoothed read (alphafold-disorder DSSP, "
            "+/-12) dissented the other way, and smoothing explains it: a "
            "25-residue window pulls the adjacent 13-residue 159-171 coil across 156-158 "
            "and washes this 7-residue coil into its flanking helices, which is why that "
            "read's two disagreements point in opposite directions -- a real structural "
            "difference would not.\n\n"
            "One honest caveat on 146 itself: it is the single alpha-basin residue in that "
            "otherwise beta/other run, i.e. the LEAST irregular position in it. 145 and "
            "143 are both 'other' basin, better founded on the basin axis alone. But 145 "
            "FAILS the clearance that matters, against BOTH cysteines: 8.5 A heavy-to-SG "
            "and 6.7 A heavy-to-heavy from C167, and 9.8 A from C183. That is a "
            "Y161-class marginal -- passing on CA-CA and CA-SG, failing on both "
            "heavy-atom conventions. 146 clears every convention against both cysteines "
            "(14.7 / 13.2 from C167, 15.6 from C183) and is the only residue in 143-149 "
            "that does, so it is better defended than when it was attacked. "
            "Not moved -- across 140-192 the "
            "per-residue basin and the DSSP call disagree at 7 of 53 scattered positions, "
            "so neither method survives being read one residue at a time and THE RUN IS "
            "THE UNIT OF EVIDENCE. 146 sits in a run with no regular element, which is the "
            "question that matters for threading an insert.\n\n"
            "After-157 was proposed and REJECTED: 153-158 is a clean six-residue alpha run "
            "by both dihedrals and pydssp, so that move would have taken a shipped "
            "construct out of a coil and put it in a helix -- Q187 in reverse",
        ),
        "Y161": (
            161,
            "the pipeline's rank-1 CD20 site, corrected 2026-10-18 from after-162 back to "
            "after-161. This was an OFF-BY-ONE, not a repoint: the endorsed site is "
            "insert-after-161 (between Y161 and I162) and is LABELLED by the residue "
            'before the junction, so reading "Y161" as "move one residue to I162" '
            "stepped off it. The residue matters -- measured on "
            "AF-P11836-F1-model_v6 here: after-161 clears C167 by 11.9 A CA-SG against "
            'after-162\'s 8.3 A. The earlier "20/26 residues of room" measured the C side '
            "to TM4 rather than to a cysteine.\n\n"
            "CORRECTION 2026-10-18, and it narrows this arm's margin: C167 is NOT the "
            "constraining partner from 161 -- C183 is. The bond is C167-C183 and C183 is "
            "the nearer of the pair from this junction. All three conventions, after-161 "
            "to C183: 11.3 A CA-CA, 10.0 A CA-SG, 8.0 A nearest-heavy-atom. So it PASSES "
            "the gate's own 10 A CA-CA floor at 11.3, sits exactly on it in CA-SG, and "
            "FAILS it in heavy-atom -- the rule that retired after-162 and Q187.\n\n"
            "RESOLVED 2026-10-18: this arm SHIPS AT 161. The clearance calibration "
            "found 8.0 A heavy-heavy unremarkable among published working sites -- 26 "
            "of 130 (20.0%) sit below it, and 11 of 49 functionally-validated sites "
            "(22.4%) sit below it, so it is around the 20th percentile rather than "
            "near the tightest. The question dissolved WITHOUT anyone settling which "
            "convention is physically right, which is the only honest way it could "
            "have closed. The study went further: the gate's own 10 A CA-CA floor "
            "REJECTS 22 of 130 published working sites (16.9%) and 11 of 49 "
            "functionally-validated ones (22.4%), so the floor is too strict rather "
            "than 161 being marginal. ADCYAP1R1 D24 is the load-bearing case: "
            "surface_and_function, position validated, inserted ONE residue from "
            "bonded C25 at 3.84 A CA-CA, so any floor at or above 4 A CA-CA rejects a "
            "site whose tagged protein both reached the surface and still worked. "
            "Note the corpus is positive-only: 'too strict' is established, 'permissive "
            "enough' is not obtainable from it. T159 remains the convention-proof "
            "alternative: inside the 159-171 DSSP loop that both independent reads agree "
            "on, 159 and 160 are the ONLY positions clearing 10 A on all three "
            "conventions (159 is 11.7 / 11.4 / 10.4) and 160 is a proline. NOT moved: the "
            "arm passes the adopted convention, so this is a judgement call rather than a "
            "defect fix",
        ),
    },
    "CCR5": {
        "F263": (
            263,
            "ECL3 (258-277, 20 aa, v1), moved off S271 on 2026-10-18 and then off F264 to "
            "the deterministic pipeline's own published site F263 (det_path surface_loop, "
            "residue_range F260-N268). S271 failed three ways on AF-P51681-F1-model_v6: "
            "pydssp calls it H in a helix contiguous from 269 into TM7 with no break, so "
            "the distal half of this loop is TM7's N-terminal extension; 5.6 A CA-CA from "
            "C269, under the 10 A floor; and buried at RSA 0.06. The earlier note called "
            "271 safe for sitting two residues PAST C269 and outside the disulfide "
            "bracket, which is the Q187 error -- outside a span is not clearance. "
            "263 is the one arm in this batch that needs no caveat: coil in a NINE-residue "
            "run 260-268, which clears the DSSP loop floor of 8 outright, RSA 0.70, "
            "15.1 A CA-CA, pLDDT 91.6. It beats 264 on both exposure (0.70 against 0.39) "
            "and clearance (15.1 against 13.7). Same loop throughout, so the position "
            "contrast against the p0042 N-terminal arm is unchanged.\n\n"
            "Clearances here are CA-CA, the convention features.py calibrates its 10 A "
            "floor on. CA-SG and nearest-heavy-atom are NOT a fixed offset from it and can "
            "fall either side: after-161 on MS4A1 is 11.3 A CA-CA but 10.0 CA-SG, while "
            "after-263 here is 15.1 CA-CA and 17.8 CA-SG. Mixing conventions at the "
            "boundary is a category error.\n\n"
            "THE SCORE WOULD HAVE PICKED THE HELIX. The v4 tag-site model's absolute "
            "scores across this loop run 271 = 0.3756, 264 = 0.3496, 263 = 0.3352 -- "
            "exactly inverse to the structural ordering, preferring the helical residue "
            "inside C269's reach and ranking the one good coil residue 34th of 70. The "
            "session that owns that model says it would not have found 263 and would have "
            "argued against it. Recorded here because it is the clearest evidence in this "
            "panel that the score must not be the tiebreaker on a GPCR loop: geometry and "
            "secondary structure decide, and the score ranks within what they allow.\n\n"
            "The cause is a GAP IN WHAT THE FEATURE SET REPRESENTS, not a bad feature. An "
            "earlier version of this note blamed `ss3`, the sequence-predicted secondary "
            "structure; an ablation refutes that. ss3 reads H at 271, agreeing with the "
            "structure at the position in question, and forcing it to coil LOWERS the "
            "score by ~0.007 everywhere, so the model mildly penalises ss3-coil rather "
            "than rewarding it -- a seventh of the 0.046 gap, with the sign reversed. What "
            "ranks 271 top is the continuous block: rsa 0.373, disorder 0.165, "
            "kd_hydropathy_w7 -0.673, against 263's 0.337 / 0.129 / +0.193. With no "
            "cysteine geometry and no structure-derived secondary structure among the 13 "
            "features, the model orders a region by exposure, local disorder and "
            "hydrophilicity -- and in a GPCR extracellular loop the most exposed, most "
            "hydrophilic, most disorder-like stretch is the one abutting TM7, which is the "
            "helix extension. It is working as designed on axes that do not encode the "
            "constraint that matters, so the failure is structural to the feature set and "
            "not fixable by cleaning one column. ss3 is separately unreliable -- it "
            "disagrees with model DSSP on 18% of extracellular residues across these five "
            "proteins, 36% on MS4A1 -- so every v4 score quoted in this folder inherits "
            "that, but it is not what produced the inversion",
        ),
    },
}


def _l1_cds() -> dict[str, str]:
    """Gene -> its L1 designed CDS, keyed ``lib:<GENE>`` in the library's own record."""
    import json

    from deliverome_analysis.cloning import library_design as ld
    from deliverome_analysis.core.s3 import processed

    with open(processed.download(ld.CDS_KEY)) as fh:
        raw = json.load(fh)
    return {g: raw[f"lib:{g}"]["cds"] for g in VARIANTS if f"lib:{g}" in raw}


def build(gene: str, variant: str, cds: str, after: int, barcode: str, tag: str) -> dict:
    """``tag`` is this construct's own encoding of LOOP_TAG_AA, from
    ``cassette_encodings.tsv`` (issued by ``diversify_cassettes.py``). Every loop arm used
    to carry the single LOOP_TAG constant, which put the same 57 nt in all of them.

    Not ``loop_cassette_encodings.tsv``: that file holds the retired 19-aa ALFA-only
    peptide and is vestigial now that every arm carries the 51-aa cassette. Reading it
    here would hand a loop arm a 57-nt encoding of the wrong peptide.
    """
    if seqtools.translate(tag) != LOOP_TAG_AA:
        raise RuntimeError(f"{gene} {variant}: issued tag does not encode the loop peptide")
    body = cds[:-3]  # drop the native stop
    cut = after * 3
    coding = body[:cut] + tag + body[cut:]
    protein = seqtools.translate(coding + STOP)

    native = seqtools.translate(cds).rstrip("*")
    if protein.rstrip("*") != native[:after] + LOOP_TAG_AA + native[after:]:
        raise RuntimeError(f"{gene} {variant}: insertion did not land where intended")
    if "*" in protein[:-1]:
        raise RuntimeError(f"{gene} {variant}: internal stop")

    insert = KOZAK + coding + STOP + SPACER + FWD_PRIMER + barcode + REV_PRIMER + CS1
    sites = seqtools.find_sites(insert, ENZYMES)
    if sites:
        raise RuntimeError(f"{gene} {variant}: insertion created Type IIS site(s): {sites}")
    return {
        "gene": gene,
        "variant": variant,
        "insert_after_residue": after,
        "residue_before": native[after - 1],
        "residue_after": native[after],
        "tag": "SpyTag003/ALFA tandem, P-flanked ALFA, GGGGS linkers",
        "tag_aa": len(LOOP_TAG_AA),
        "protein_aa": len(protein) - 1,
        "insert_nt": len(insert),
        "barcode": barcode,
        "provenance": VARIANTS[gene][variant][1],
        "insert": insert,
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, default=HERE)
    args = parser.parse_args(argv)

    # Same rule as the GPCR rescue arms: the ORF is the L1 library's designed CDS, not
    # this folder's panel copy, so a loop construct encodes the protein the library
    # actually builds. Asserted, not assumed -- a mismatch stops the build.
    # Build from the DOMESTICATED copy, and check it against the L1 library's CDS.
    #
    # Both are needed and neither substitutes for the other. ``candidate_cds.json`` is the
    # library's record of what it synthesises but it is NOT Type IIS-cleared: building the
    # loop arms straight off it reintroduced five Esp3I/BbsI/BsaI sites into SLC2A4, because
    # ``build()`` here inserts a cassette and screens, it does not domesticate the way the
    # GPCR arms do. ``orf_domesticated.fasta`` is the same protein with those codons
    # cleared, so it is the buildable form.
    #
    # What the L1 comparison buys is the guarantee the source switch was really after: that
    # the panel's domesticated copy still encodes the library's protein. Domestication is
    # synonymous, so the proteins must match exactly; if they ever diverge, this folder's
    # copy has drifted from the library and the arm would answer for a different proteoform.
    cds = read_fasta(HERE / "orf_domesticated.fasta")
    l1 = _l1_cds()
    for gene in sorted(VARIANTS):
        if gene not in l1:
            raise SystemExit(f"{gene} is not an L1 library reference — no CDS to check against")
        a = seqtools.translate(l1[gene]).rstrip("*")
        b = seqtools.translate(cds[gene]).rstrip("*")
        if a != b:
            raise SystemExit(
                f"{gene}: the panel's domesticated ORF and the L1 designed CDS encode "
                f"different proteins ({len(b)} aa vs {len(a)} aa) — the panel copy has "
                f"drifted from the library; resolve the proteoform before building"
            )
    print(f"L1 proteoform check: {len(VARIANTS)} gene(s) match the library's designed CDS")
    barcodes = {
        (r["gene"], r["variant"]): r["barcode"]
        for r in csv.DictReader((HERE / "barcodes_loop.tsv").open(), delimiter="\t")
    }
    tags = {
        r["construct"] if "construct" in r else r["gene"]: r["cassette_nt"]
        for r in csv.DictReader((HERE / "cassette_encodings.tsv").open(), delimiter="\t")
    }

    rows = []
    for gene, variants in VARIANTS.items():
        for variant, (after, _) in variants.items():
            key = f"{gene}_{variant}"
            if key not in tags:
                raise SystemExit(
                    f"{key}: no loop-tag encoding issued — rerun diversify_cassettes.py"
                )
            record = build(gene, variant, cds[gene], after, barcodes[(gene, variant)], tags[key])
            rows.append(record)
            print(
                f"{gene:8s} {variant:6s} after {record['insert_after_residue']:4d} "
                f"({record['residue_before']}|{record['residue_after']})  "
                f"insert {record['insert_nt']:5d} nt  protein {record['protein_aa']:4d} aa  "
                f"bc={record['barcode']}"
            )

    cols = [
        "gene",
        "variant",
        "insert_after_residue",
        "residue_before",
        "residue_after",
        "tag",
        "tag_aa",
        "protein_aa",
        "insert_nt",
        "barcode",
        "provenance",
    ]
    table = args.out_dir / "constructs_loop.tsv"
    with table.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=cols, delimiter="\t", extrasaction="ignore", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)

    fasta = args.out_dir / "constructs_loop.fasta"
    fasta.write_text(
        "".join(
            f">{r['gene']}_{r['variant']}|loop-insert|bc={r['barcode']}|{r['insert_nt']}nt\n"
            + "\n".join(r["insert"][i : i + 60] for i in range(0, len(r["insert"]), 60))
            + "\n"
            for r in rows
        )
    )
    print(f"\nWrote {table} and {fasta}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
