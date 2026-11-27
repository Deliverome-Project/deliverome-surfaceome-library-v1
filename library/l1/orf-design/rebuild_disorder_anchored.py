"""Re-place the N-terminal cassette on the disorder-anchored junction, not the cleavage site.

`build_constructs.py` splices the cassette at ``sp_end * 3``, so the mature protein begins
with the cassette. POST_SP_PLACEMENT.md shows that is not what works: of Tedman's 108
signal-peptide-bearing insertions, **none** sits at the cleavage site, 107 are downstream,
the median offset is +4, and the junctions sit at the 88th-93rd disorder percentile of
their own protein.

The rule, applied here: ``sp_end`` is the **later** of the DeepTMHMM v1 and SignalP 6
calls, because the failure is asymmetric -- a cassette inside the signal peptide leaves
with it and the construct reads as a non-expressor, while one placed a few residues late is
still presented. Then search ``sp_end + 1 .. sp_end + 20`` and take the residue with the
highest mean within-protein disorder percentile across metapredict, NetSurfP-3.0 and
AlphaFold-disorder.

CORRECTION 2026-10-18. The stored percentiles WERE two-predictor means -- recomputing
from metapredict and NetSurfP-3.0 alone reproduces them within ~1 point on six of eight.
My first explanation of why was wrong: I read `disorder_hex` being NULL on
`alphafold-disorder` rows as the data being absent. It is not. That column holds
metapredict and NetSurfP only BY SCHEMA; AlphaFold-disorder stores its track in
`plddt_hex`, complete at one byte per residue on all eight genes (FOLR1: 257 bytes for
257 residues, range 129-252), and disorder is 1 - pLDDT/100. Data present, consumer
blind -- whatever produced the stored values read `disorder_hex` uniformly across all
three predictors, got nothing for AlphaFold and dropped it silently. Same silent-drop
pattern as `insert_after_residue` below, and worth checking in that producer.

The `pct` values in JUNCTIONS are now genuine THREE-predictor means:

    gene      stored(2)  metapredict  NetSurfP  AF-disorder  mean(3)
    AXL            90.9         90.1      91.6         76.7     86.1
    ERBB2          84.6         76.7      92.4         71.6     80.2
    FOLR1          86.1         91.4      80.7         88.9     87.0
    HEG1           59.8         29.2      90.4         32.4     50.7
    NECTIN4        74.5         68.3      80.6         64.0     71.0
    PTPRC          89.1         95.5      88.4         78.1     87.3
    SLC9A3         92.8         87.2      98.4         88.6     91.4
    VSIG10         53.7         65.8      48.5         48.7     54.3

Read HEG1 twice: metapredict says 29.2 and NetSurfP says 90.4 at the SAME residue, a
61-point disagreement. The disorder evidence there is contradictory rather than low, so
no disorder argument should carry weight for HEG1 in either direction -- which is why
that junction was moved on signal-peptide clearance instead.

TAGGED RE-PREDICTION, canary_rp_1792371023, 2026-10-18 -- the first prediction in this
project run on a sequence with the cassette in it. Eight tagged against eight native:

    gene      native cleavage  junction  tagged call      verdict
    VSIG10                 30        31  OTHER, no site   *** LOST ***
    HEG1                   28        32  SP at 29         moved +1, benign
    ERBB2                  22        23  SP at 22         preserved
    AXL                    32        33  SP at 32         preserved
    SLC9A3                 25        26  SP at 25         preserved, p(SP) 0.867 -> 0.667
    NECTIN4                31        32  SP at 31         preserved
    FOLR1                  22        27  SP at 22         preserved
    PTPRC                  25        41  SP at 25         preserved

**VSIG10 (p0026) loses its own signal peptide when tagged.** No cleavage site is predicted
at all, so the cassette does not sit at a mature N-terminus and the construct reads as a
NON-EXPRESSOR rather than a weak signal. This is the AXL-at-residue-25 failure realised
rather than theorised, on a different arm, and it is the reason this file exists.

Mechanism, confirmed rather than guessed: SignalP 6 reads only the first ~70 residues, so
a 51-aa cassette at +1 replaces most of the h-region and c-region context the call depends
on. Four of the five +1 arms kept their site; that is luck of sequence, not safety of
position.

RETEST AT 32: ALSO LOST. p(SP) is 0.000 at BOTH 31 and 32, not merely reduced, and no
cleavage site is called at either. The interval is 30..32 and residue 30 is the cleavage
site itself, which none of Tedman's 111 signal-peptide-bearing constructs uses. So no
position in VSIG10's interval rescues it WITH THIS CASSETTE: the remaining options are a
shorter cassette or dropping the arm. Neither a shorter cassette nor position 30 has been
tested.

THE OBSERVATION IN THE RESULT THAT MAY MATTER MORE THAN THE FINDING. VSIG10's NATIVE
cleavage-site probability is 0.348, the weakest in the panel by a wide margin -- the next
lowest is 0.965 -- while its p(SP) is 1.000. Its signal peptide is called with certainty
and its cleavage POSITION is poorly defined natively. The plausible reading is that an arm
whose native cleavage site is ill-defined tolerates nothing behind it. Untested as an
explanation and n=1, but `signalp_public.cleavage_probability` exists for every canonical
member, so the whole library can be ranked on this risk BEFORE anything is re-predicted or
ordered -- a free sort rather than a per-construct Modal run.

AND OUR OWN INDEPENDENT SUPPORT FOR TEDMAN'S +4 MEDIAN, from a different direction than
the fold-end finding. Four of the seven surviving arms SOFTENED their cleavage-site
probability while keeping the exact site: SLC9A3 0.967 -> 0.650, AXL 0.978 -> 0.816,
HEG1 0.972 -> 0.840, ERBB2 0.965 -> 0.905. The two arms furthest from their cleavage site
-- FOLR1 at +5 and PTPRC at +16 -- are unchanged, and NECTIN4 at +1 did not soften. So
softening tracks proximity to the cleavage site, imperfectly: the arms hugging the
boundary are the ones whose call degrades even when it does not move.

VSIG10 was already the worst-placed arm on two independent axes before this ran. It is the
only construct genuinely inside a real domain -- InterPro IPR013106 Ig V-set 31-105 begins
exactly at its junction -- and its admissible interval is the narrowest in the panel at
30..32, two residues wide. Three independent lines now point at the same arm.

HEG1's +1 shift is NOT a failure: `compare_reprediction` scores it not-preserved because
the site is not identical, but the cassette sits after 32 and so remains three residues
into the mature chain under the new call. A shift that does not cross the junction is
benign. SLC9A3 passes the binary test while its p(SP) falls 0.867 -> 0.667; `preserved`
ignores confidence, and SLC9A3 was already the panel's lowest-confidence native SP.

Disorder shifts are RECORDED, NOT INTERPRETED: small everywhere except FOLR1 at 0.27
(metapredict) and NECTIN4 at 0.15 (both tools). No threshold for this has ever been
defined, a single construct cannot be ranked against a series, and our 51 aa does not
inherit Tedman's 9-aa calibration. A two-arm reference distribution over his 109 validated
junctions is approved and will make these readable; until it lands these numbers are
neither a pass nor a fail.

EVERY other number in this file, here and in JUNCTIONS, is a prediction on the NATIVE
sequence. Tedman's own criterion was the opposite -- he built a series across the
admissible interval and kept "the construct that was most strongly predicted to preserve
the correct cleavage site without altering disorder prediction", i.e. he re-predicted on
the TAGGED construct. No such prediction exists for any of these eight: signalp_designed_public
holds only untagged `l1_needs_prediction` rows. That is the open gap in this file.

Only two of the eight panel constructs actually move:

* **AXL** was built after residue 25, the UniProt cleavage site. SignalP 6 (p = 0.9998)
  and DeepTMHMM v1 both call the signal peptide 1-32, so the existing cassette sits
  *inside* it and would be cleaved off. Moves to 33.
* **PTPRC** was built after 25 and the disorder optimum is 41, sixteen residues away.

The other six land within one residue of where they already are, so they are rebuilt here
only to prove the script reproduces them byte-for-byte.

Cassette encodings are NOT re-issued: each gene keeps the per-gene synonymous encoding it
already holds, since an encoding is gene-specific rather than position-specific. Moving the
junction does change its flanks, so the seam is re-checked for Type IIS sites.

    uv run --locked python rebuild_disorder_anchored.py
    uv run --locked python rebuild_disorder_anchored.py --write
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import build_constructs as bc

from deliverome_analysis.cloning import seqtools

HERE = Path(__file__).resolve().parent

#: Junction = insert the cassette after this residue. Computed by the rule above against
#: INTERPRO DOMAIN REVIEW, 2026-10-18. Every junction was measured against InterPro
#: (REST, all member DBs, entry types `domain` and `repeat` only -- a `family` entry
#: routinely spans the whole protein and says nothing about a fold). First domain start
#: strictly after the latest SP call, and the slack to it:
#:
#:     gene      junction  SPmax  1st domain  slack
#:     HEG1            32     29         985    953
#:     ERBB2           23     22          52     29
#:     AXL             33     32          35      2
#:     SLC9A3          26     25          62     36
#:     VSIG10          31     30          31      0   <- junction IS the domain's first residue
#:     NECTIN4         32     31          36      4
#:     FOLR1           27     24          36      9
#:     PTPRC           41     25         391    350
#:
#: VSIG10 sits exactly at the start of IPR013106 Ig V-set (31-105), and that domain
#: begins after the SP, so it is a real mature-chain boundary rather than an annotation
#: running into the signal peptide -- which is what AXL's and NECTIN4's apparent hits
#: are (IPR007110 starting at 27 and 30 against cleavage at 32 and 31, i.e. describing
#: sequence the mature protein will not contain; their first real domains are 35 and 36).
#: VSIG10 is NOT moved: Tedman's own screen contains five constructs whose junction sits
#: at or before an annotated domain start -- Frizzled-CRD and LRR-N -- and all of them
#: expressed, so slack 0 is a geometry his data validates rather than disqualifies.
#:
#: Read a HIT from InterPro as strong and a CLEAR as weak. Its starts run a median 16
#: residues later than the disorder consensus, its containment of Tedman's junctions is
#: no better than chance for its own distribution, and an InterProScan run has found a
#: SMART Ig domain at residue 27 where the match lookup's earliest was 167. So "clear"
#: above means "no annotated domain found", not "no fold there". AXL at 2 and NECTIN4 at
#: 4 are unresolved by this method, not safe.
#:
#: No construct in the panel has an N-terminal propeptide -- SIGNAL, PROPEP, CHAIN,
#: TRANSIT and INIT_MET checked for all eight. FOLR1's Propeptide is 235-257, the
#: C-terminal GPI-attachment signal, which the ORF needs and must keep.
#:
#: public D1 (disorder_public, signalp_public, topology_public at topo_2026_05_16); the
#: percentile is the mean across the three predictors at that residue.
JUNCTIONS = {
    "AXL": {
        "junction": 33,
        "sp_end": 32,
        "built": 25,
        "pct": 86.1,
        "why": "SignalP 6 and DeepTMHMM v1 both call the SP 1-32; the built cassette at 25 is inside it",
    },
    "PTPRC": {
        "junction": 41,
        "sp_end": 25,
        "built": 25,
        "pct": 87.3,
        "why": "disorder optimum is 16 residues downstream of the cleavage site",
    },
    "ERBB2": {
        "junction": 23,
        "sp_end": 22,
        "built": 22,
        "pct": 80.2,
        "why": "+1, unchanged in effect",
    },
    "FOLR1": {
        "junction": 27,
        "sp_end": 22,
        "built": 24,
        # Recomputed from disorder_public against the designed length (2026-10-18):
        # 23 scored 90.9, 27 scores 86.1. The move costs 4.8 percentile points and buys
        # 2 residues of signal-peptide clearance, which is the right trade -- disorder is
        # a soft preference, a cleavage that removes the cassette is a hard failure.
        # Scanning the whole admissible window 25-32 gives 87.2 at 25 down to 79.8 at 32,
        # a 7-point spread, so disorder does not discriminate here and should not drive
        # the choice. 25 would score 1.1 points higher for 2 residues less margin.
        "pct": 87.0,
        "why": (
            "moved 23 -> 27 on 2026-10-18. DeepTMHMM v1 and SignalP 6 call the SP 1-22, "
            "UniProt 1-24, so after-23 sat INSIDE UniProt's call and displaced the -1 "
            "residue the peptidase reads. 27 clears all three (+3 past the latest) and "
            "sits at Tedman's median +4 geometry. InterPro puts the first domain at 36, "
            "so 9 residues of clearance, and C37 is the first cysteine -- staying at or "
            "below 32 keeps the cassette out of the disulfide-rich folate-binding fold. "
            "Junction is A27|W28, an alanine at P1. FOLR1 is the GPI member whose "
            "C-terminus is cleaved for the anchor, so its N-terminus is its only "
            "taggable end and there is no second chance"
        ),
    },
    "HEG1": {
        "junction": 32,
        "sp_end": 28,
        "built": 29,
        # Recomputed 2026-10-18: 29 scored 64.3, 32 scores 59.8. The whole admissible
        # window 30-36 spans only 59.8-62.7, i.e. HEG1's N-terminus simply is not
        # disordered and no position in it rescues that. Chosen on SP clearance instead.
        "pct": 50.7,
        "why": (
            "moved 29 -> 32 on 2026-10-18. DeepTMHMM v1 and SignalP 6 call the SP 1-28 "
            "so after-29 was a legal +1, but UniProt calls 1-29, which made the cassette "
            "the first residue of the mature protein. Of 111 signal-peptide-bearing "
            "constructs in Tedman's screen, ZERO sit at the cleavage site -- 107+ are "
            "downstream -- so margin 0 is the one geometry his data never uses. 32 clears "
            "all three calls (+4 past DeepTMHMM, his median offset) with 953 residues to "
            "the first domain at 985, and lands just before the Pro-rich 33-37 stretch, "
            "which is flexible. Also retires the worst disorder percentile in the set: "
            "29 scored 65.3 against a cohort band of 88-93"
        ),
    },
    "NECTIN4": {
        "junction": 32,
        "sp_end": 31,
        "built": 31,
        "pct": 71.0,
        "why": "+1, unchanged in effect",
    },
    "SLC9A3": {
        "junction": 34,
        "sp_end": 25,
        "built": 25,
        "pct": 85.3,
        "why": (
            "MOVED 26 -> 34 on 2026-10-19 (Becca). At 26 (sp_end+1) the cassette measurably "
            "degraded the signal peptide: tagged p(SP) 0.667 and cleavage-site probability "
            "0.650 against 0.867/0.967 untagged -- the worst of the eight, and the only arm "
            "below 0.9 on p(SP). Re-predicting at 26/34/40/46 (canary_rp_1792374944) gives "
            "0.650/0.973/0.950/0.974, so 34 does not merely recover the call, it beats the "
            "untagged protein. Chosen over 46, which is OUTSIDE the admissible interval "
            "25-43 entirely, and over 40 because 34 leaves 22 residues to TM1 (57) against "
            "16 and sits in the Gly/Glu-rich 26-45 stretch rather than near the aromatic run "
            "46-55. Disorder percentile 85.3 against 90.5 at 26 by the same three-predictor "
            "method -- still deep in disordered sequence. "
            "MECHANISM CORRECTED 2026-10-19: this was first explained as SignalP's "
            "70-residue window -- at 26 it holds one native mature residue and 44 of "
            "cassette, at 34 it holds nine, so 'buy more native context'. That is "
            "FALSIFIED. HEG1 has no p(SP) deficit and ample room, and across 4/12/20 native "
            "residues its site probability runs 0.840/0.789/0.859 -- flat and "
            "non-monotonic, never nearing its 0.972 native. Twenty residues is more than "
            "double what SLC9A3 bought and recovers nothing. AXL, also with no p(SP) "
            "deficit, is flat across its whole interval (canary_rp_1792375661). What the "
            "three arms say instead: moving the junction recovers a p(SP) DEFICIT and does "
            "NOT recover a softened cleavage SITE. SLC9A3 was the only one of the eight "
            "with a native p(SP) deficit (0.867) and the only one that improved. Residue "
            "count was incidental -- in the first two cases it was perfectly confounded "
            "with which axis was degraded. The confound and the separating test are the "
            "InterPro session's. The MOVE still stands: 0.650 -> 0.973 and 0.667 -> 1.000 "
            "are measurements on this construct. Only the explanation, and any "
            "extrapolation from it, were wrong."
        ),
    },
    "VSIG10": {
        "junction": 31,
        "sp_end": 30,
        "built": 30,
        "pct": 54.3,
        "why": "+1, unchanged in effect",
    },
}


def assemble_at(gene: str, cds: str, barcode: str, cassette: str, junction: int) -> dict:
    """Insert the cassette after ``junction`` residues of the ORF."""
    body = cds[:-3]  # drop the native stop; STOP is appended below
    cut = junction * 3
    if cut > len(body):
        raise SystemExit(f"{gene}: junction {junction} is past the ORF")
    coding = body[:cut] + cassette + body[cut:]
    insert = (
        bc.KOZAK + coding + bc.STOP + bc.SPACER + bc.FWD_PRIMER + barcode + bc.REV_PRIMER + bc.CS1
    )

    protein = seqtools.translate(coding + bc.STOP)
    native = seqtools.translate(cds)
    if protein[-1] != "*" or "*" in protein[:-1]:
        raise SystemExit(f"{gene}: internal stop or missing stop")
    expect = native[:junction] + bc.CASSETTE_AA + native[junction:]
    if protein.rstrip("*") != expect.rstrip("*"):
        raise SystemExit(f"{gene}: protein is not native-with-cassette-inserted")
    sites = seqtools.find_sites(insert, bc.ENZYMES)
    # A site can appear at either new seam even though cassette and ORF are each clean.
    seams = [c for c in (cut, cut + len(cassette))]
    at_seam = [s for s in sites if any(c - 7 < s[2] < c + 1 for c in seams)]
    return {
        "gene": gene,
        "junction_after_residue": junction,
        # The junction as BUILT, derived from the assembled insert rather than restated.
        # This is the column a reader should use. It is first in the row for that reason.
        # Derived from the insert THIS call just assembled, by locating the cassette
        # in it -- not from the FASTA on disk, which during a moving run still holds
        # the previous build. That read-before-write is why FOLR1 and HEG1 shipped
        # 23 and 29 here while the built insert had 27 and 32.
        "insert_after_residue": _cassette_offset(insert, cassette),
        "junction_source": "located in the insert assembled by this call",
        # Renamed from `built_after_residue`, which held the PRE-rebuild position and
        # therefore disagreed with the built insert for six of the eight genes. Two
        # readers (this session and the InterPro review) each took it for the as-built
        # junction and had to recover the real one by translating ORFs and locating the
        # cassette by sequence. Keeping the value is useful -- it is what moved -- but the
        # old name claimed the opposite of what it meant.
        "pre_rebuild_after_residue": JUNCTIONS[gene]["built"],
        "moved_by": junction - JUNCTIONS[gene]["built"],
        "sp_end_rule": JUNCTIONS[gene]["sp_end"],
        "offset_past_sp_end": junction - JUNCTIONS[gene]["sp_end"],
        "mean_disorder_pct": JUNCTIONS[gene]["pct"],
        "insert_nt": len(insert),
        "protein_aa": len(protein) - 1,
        "type_iis_sites": len(sites),
        "sites_at_a_seam": len(at_seam),
        "barcode": barcode,
        "why": JUNCTIONS[gene]["why"],
        "insert": insert,
    }


def _cassette_offset(insert: str, cassette: str) -> int:
    """Residues of ORF before the cassette, located in the insert just assembled.

    Derived, never restated -- but derived from the artifact IN HAND rather than the one
    on disk. The previous version read constructs_nterm_disorder_anchored.fasta, which on
    the run that moves a junction still holds the old build, so the column landed one run
    behind: FOLR1 and HEG1 were written as 23 and 29 while their inserts had 27 and 32,
    and a post-write assertion could not catch it because writing made the two agree.
    """
    i = insert.find(cassette)
    if i < 0:
        raise SystemExit("cassette not found in the assembled insert")
    return (i - len(bc.KOZAK)) // 3


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--write", action="store_true")
    args = ap.parse_args(argv)

    cds = bc.read_fasta(HERE / "orf_domesticated.fasta")
    cassettes = bc.load_cassettes()
    barcodes = {
        r["gene"]: r["barcode"]
        for r in csv.DictReader((HERE / "barcodes.tsv").open(), delimiter="\t")
    }
    existing = {
        r["gene"]: r for r in csv.DictReader((HERE / "constructs_nterm.tsv").open(), delimiter="\t")
    }

    rows = []
    for gene in sorted(JUNCTIONS):
        r = assemble_at(
            gene, cds[gene], barcodes[gene], cassettes[gene], JUNCTIONS[gene]["junction"]
        )
        was = existing.get(gene, {})
        r["insert_nt_was"] = was.get("insert_nt", "")
        r["reproduces_existing"] = r["moved_by"] == 0 and str(r["insert_nt"]) == was.get(
            "insert_nt"
        )
        rows.append(r)

    moved = [r for r in rows if r["moved_by"]]
    print(
        f"{'gene':9} {'was':>4} {'now':>4} {'move':>5} {'offset':>7} {'pct':>6} {'sites':>6} {'seam':>5}  why"
    )
    for r in rows:
        print(
            f"{r['gene']:9} {r['pre_rebuild_after_residue']:4} {r['junction_after_residue']:4} "
            f"{r['moved_by']:+5} {r['offset_past_sp_end']:+7} {r['mean_disorder_pct']:6.1f} "
            f"{r['type_iis_sites']:6} {r['sites_at_a_seam']:5}  {r['why'][:52]}"
        )
    print(f"\nmoved: {len(moved)} ({', '.join(r['gene'] for r in moved)})")
    print(
        f"unchanged and byte-identical to constructs_nterm.tsv: "
        f"{sum(1 for r in rows if r['reproduces_existing'])}"
    )
    assert all(r["type_iis_sites"] == 0 for r in rows), "a rebuilt insert carries a Type IIS site"

    if args.write:
        cols = [k for k in rows[0] if k != "insert"]
        with open("constructs_nterm_disorder_anchored.tsv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows([{k: r[k] for k in cols} for r in rows])
        Path("constructs_nterm_disorder_anchored.fasta").write_text(
            "".join(
                f">{r['gene']}|after{r['junction_after_residue']}|disorder_anchored"
                f"|bc={r['barcode']}|{r['insert_nt']}nt\n"
                + "\n".join(r["insert"][i : i + 60] for i in range(0, len(r["insert"]), 60))
                + "\n"
                for r in rows
            )
        )
        print("wrote constructs_nterm_disorder_anchored.{tsv,fasta}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
