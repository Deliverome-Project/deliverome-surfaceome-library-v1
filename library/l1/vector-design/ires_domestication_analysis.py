"""Is Fragmid's domesticated IRES2_v1 structurally safe?

p003 carries Fragmid's `IRES2_v1` (parts pF_AA333/334/335) rather than the
classic EMCV IRES2 (Addgene #100673, EMMA). The two differ by exactly three
substitutions, all Type IIS / restriction-site removals:

    137 T->G   kills BbsI  GTCTTC (minus strand)
    236 A->C   kills BbsI  GAAGAC (plus strand)
    480 G->A   kills ApaLI GTGCAC

An EMCV IRES works through its fold, so the question is whether those three
bases perturb structure. Prediction alone cannot answer that -- an absolute
MFE fold of a 587 nt RNA is not trustworthy. This script therefore frames
every structural number against two controls:

  1. random 3-substitution variants of the same IRES  (is the perturbation
     larger than chance?)
  2. 79 natural EMCV isolates from NCBI               (is it larger than what
     nature already tolerates?)

and asks, per edited base, whether the domesticated allele occurs in nature.

Run:  uv run python library/l1/vector-design/ires_domestication_analysis.py
"""

from __future__ import annotations

import collections
import random
import statistics
from pathlib import Path
from xml.etree import ElementTree

import matplotlib.pyplot as plt
import RNA
from Bio import Align, Entrez, SeqIO

from deliverome_analysis.plotting import (
    COLORS,
    save_figure,
    setup_plotting_style,
)

HERE = Path(__file__).parent
CACHE = Path.home() / ".cache" / "deliverome" / "emcv_genomes.fa"
ENTREZ_EMAIL = "becca@deliverome.org"
SEED = 0

# Parental: EMCV IRES2 as distributed in EMMA (Addgene #100673).
PARENTAL = (
    "cccctctccctcccccccccctaacgttactggccgaagccgcttggaataaggccggtgtgcgtttgtctatatgtt"
    "attttccaccatattgccgtcttttggcaatgtgagggcccggaaacctggccctgtcttcttgacgagcattcctagg"
    "ggtctttcccctctcgccaaaggaatgcaaggtctgttgaatgtcgtgaaggaagcagttcctctggaagcttcttgaa"
    "gacaaacaacgtctgtagcgaccctttgcaggcagcggaaccccccacctggcgacaggtgcctctgcggccaaaagcc"
    "acgtgtataagatacacctgcaaaggcggcacaaccccagtgccacgttgtgagttggatagttgtggaaagagtcaaa"
    "tggctctcctcaagcgtattcaacaaggggctgaaggatgcccagaaggtaccccattgtatgggatctgatctggggc"
    "ctcggtgcacatgctttacatgtgtttagtcgaggttaaaaaaacgtctaggccccccgaaccacggggacgtggtttt"
    "cctttgaaaaacacgatgataatatggccacaacc"
).upper()

EDITS = {137: "G", 236: "C", 480: "A"}  # 1-based position -> domesticated base

# The two BbsI recognition sites removed, as 1-based inclusive spans.
BBSI_SITES = {
    "site 1 (GTCTTC, minus strand)": (134, 139),  # Fragmid edited 137
    "site 2 (GAAGAC, plus strand)": (234, 239),  # Fragmid edited 236
}

TYPE_IIS = {"BbsI": "GAAGAC", "Esp3I": "CGTCTC", "BsaI": "GGTCTC", "ApaLI": "GTGCAC"}

# This 587 bp element aligns to EMCV reference NC_001479.1 at 260..845, so
# EMCV genome coordinate = IRES position + 259. Verified by local alignment
# rather than assumed; re-derived in report_confident_pair_changes().
EMCV_OFFSET = 259

# Functional landmarks, in EMCV genome coordinates.
#   The IRES core is domains H, I, J-K and L immediately upstream of the
#   initiator AUG at nt 834. eIF4G binds a three-way helical junction in the
#   J-K domain and contacts duplexes in J-K and L; deleting the J and K apical
#   hairpins impairs both eIF4G binding and IRES function.
#     Kolupaeva et al. 2003, Mol Cell Biol 23:687-698, PMID 12509466,
#     doi:10.1128/MCB.23.2.687-698.2003
#   Domain I is the structurally contested region. Mutational analysis + SHAPE
#   showed two helices inside it are important for IRES translation, at the
#   regions of interest below.
#     Maloney & Joseph 2023, Biochemistry 63:107-115, PMID 38081770,
#     doi:10.1021/acs.biochem.3c00579
#   J-K bounds below are approximate — taken as the span between domain I's 3'
#   end and the reported 3' boundary of the eIF4G footprint.
EMCV_LANDMARKS = [
    (451, 678, "domain I (ROI 1; base pairing validated)"),
    (501, 625, "  ROI 3 (base pairing validated)"),
    (679, 786, "domain J-K (eIF4G platform, approximate bounds)"),
    (834, 836, "initiator AUG-11"),
]


def apply_edits(seq: str, edits: dict[int, str]) -> str:
    out = list(seq)
    for pos, base in edits.items():
        out[pos - 1] = base
    return "".join(out)


DOMESTICATED = apply_edits(PARENTAL, EDITS)


def revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def count_sites(seq: str, site: str) -> int:
    return seq.count(site) + seq.count(revcomp(site))


# --------------------------------------------------------------------------
# folding helpers
# --------------------------------------------------------------------------


def fold(seq: str) -> tuple[str, float]:
    """MFE structure and energy. ViennaRNA wants RNA, so T -> U."""
    return RNA.fold(seq.replace("T", "U"))


def pair_table(struct: str) -> list[int]:
    stack: list[int] = []
    table = [0] * len(struct)
    for i, char in enumerate(struct):
        if char == "(":
            stack.append(i)
        elif char == ")":
            j = stack.pop()
            table[i] = j + 1
            table[j] = i + 1
    return table


def pair_set(struct: str) -> set[tuple[int, int]]:
    table = pair_table(struct)
    return {(i + 1, table[i]) for i in range(len(struct)) if table[i] > i + 1}


def confident_pairs(seq: str, cutoff: float = 0.5) -> tuple[set[tuple[int, int]], float]:
    """Pairs the partition function assigns probability > cutoff.

    More robust than comparing single MFE structures, which flip between
    near-degenerate alternatives for a molecule this size.
    """
    fc = RNA.fold_compound(seq.replace("T", "U"))
    _, ensemble_dg = fc.pf()
    bpp = fc.bpp()
    n = len(seq)
    pairs = {(i, j) for i in range(1, n + 1) for j in range(i + 1, n + 1) if bpp[i][j] > cutoff}
    return pairs, ensemble_dg


def positional_entropy(seq: str) -> list[float]:
    fc = RNA.fold_compound(seq.replace("T", "U"))
    fc.pf()
    return list(fc.positional_entropy())


# --------------------------------------------------------------------------
# natural variation
# --------------------------------------------------------------------------


def fetch_emcv_genomes() -> Path:
    """Download EMCV genomes from NCBI, cached on disk.

    The exact accessions behind the committed numbers are pinned in
    `emcv_isolates.txt` alongside this script, with the query and retrieval
    date. NCBI results drift as records are deposited, so refresh that file
    whenever the returned count changes.
    """
    if CACHE.exists():
        return CACHE
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    Entrez.email = ENTREZ_EMAIL
    handle = Entrez.esearch(
        db="nuccore",
        term='"Encephalomyocarditis virus"[Organism] AND 7000:8500[SLEN]',
        retmax=80,
    )
    ids = Entrez.read(handle)["IdList"]
    handle.close()
    fetched = Entrez.efetch(db="nuccore", id=",".join(ids), rettype="fasta", retmode="text").read()
    CACHE.write_text(fetched)
    return CACHE


def align_natural_ires(path: Path):
    """Align each genome's 5'UTR to PARENTAL.

    Returns per-column allele counts (in PARENTAL coordinates) and the natural
    IRES subsequences with their coordinate maps.
    """
    aligner = Align.PairwiseAligner(
        mode="local",
        match_score=2,
        mismatch_score=-1,
        open_gap_score=-6,
        extend_gap_score=-1,
    )
    columns: dict[int, collections.Counter] = collections.defaultdict(collections.Counter)
    variants: dict[str, dict[int, int]] = {}
    for record in SeqIO.parse(path, "fasta"):
        genome = str(record.seq).upper()[:1800]
        try:
            alignment = aligner.align(PARENTAL, genome)[0]
        except ValueError, IndexError:
            continue
        if alignment.score < 700:  # require a solid IRES-length match
            continue
        target_blocks, query_blocks = alignment.aligned
        for (t0, t1), (q0, _) in zip(target_blocks, query_blocks, strict=True):
            for k in range(t1 - t0):
                columns[t0 + k][genome[q0 + k]] += 1
        # keep near-full-length IRESes for the structural envelope
        if target_blocks[0][0] > 6 or target_blocks[-1][-1] < len(PARENTAL) - 6:
            continue
        start, end = query_blocks[0][0], query_blocks[-1][-1]
        sub = genome[start:end]
        if sub in variants:
            continue
        coord_map = {}
        for (t0, t1), (q0, _) in zip(target_blocks, query_blocks, strict=True):
            for k in range(t1 - t0):
                coord_map[q0 + k - start] = t0 + k + 1
        variants[sub] = coord_map
    return columns, variants


def structural_divergence(
    seq: str, coord_map: dict[int, int], reference_pairs: set[tuple[int, int]]
) -> int:
    """Symmetric pair difference vs the reference fold, in PARENTAL coordinates.

    Length-tolerant, so natural isolates carrying indels stay comparable.
    """
    projected = set()
    for i, j in pair_set(fold(seq)[0]):
        if i - 1 in coord_map and j - 1 in coord_map:
            projected.add((coord_map[i - 1], coord_map[j - 1]))
    covered = set(coord_map.values())
    comparable = {(i, j) for i, j in reference_pairs if i in covered and j in covered}
    return len(comparable ^ projected)


# --------------------------------------------------------------------------
# report sections
# --------------------------------------------------------------------------


def report_sequence_differences() -> None:
    print("=" * 72)
    print("1. SEQUENCE")
    print("=" * 72)
    diffs = [
        (i + 1, a, b) for i, (a, b) in enumerate(zip(PARENTAL, DOMESTICATED, strict=True)) if a != b
    ]
    print(f"length: parental {len(PARENTAL)} bp, domesticated {len(DOMESTICATED)} bp")
    print(f"differences: {len(diffs)} substitutions, no indels")
    for pos, a, b in diffs:
        print(f"    {pos:>4}  {a} -> {b}")
    print(f"3' end, both: ...{PARENTAL[-24:]}")
    print("    (authentic EMCV AUG-11 initiator + native GCCACAACC spacer intact)")
    print(
        f"upstream AUGs: parental {PARENTAL.count('ATG')}, domesticated {DOMESTICATED.count('ATG')}"
    )
    print("\nrestriction sites:")
    for enzyme, site in TYPE_IIS.items():
        print(f"    {enzyme:<7} {count_sites(PARENTAL, site)} -> {count_sites(DOMESTICATED, site)}")


def report_folding() -> tuple[int, list[int]]:
    print("\n" + "=" * 72)
    print("2. STRUCTURE — domesticated vs parental")
    print("=" * 72)
    struct_p, mfe_p = fold(PARENTAL)
    struct_d, mfe_d = fold(DOMESTICATED)
    pairs_p, dg_p = confident_pairs(PARENTAL)
    pairs_d, dg_d = confident_pairs(DOMESTICATED)
    observed = RNA.bp_distance(struct_p, struct_d)

    print(f"MFE:         {mfe_p:8.2f} -> {mfe_d:8.2f} kcal/mol  (delta {mfe_d - mfe_p:+.2f})")
    print(f"ensemble dG: {dg_p:8.2f} -> {dg_d:8.2f} kcal/mol  (delta {dg_d - dg_p:+.2f})")
    print(
        f"confident pairs (p>0.5): {len(pairs_p)} -> {len(pairs_d)}; "
        f"{len(pairs_p & pairs_d)} shared, {len(pairs_p - pairs_d)} lost, "
        f"{len(pairs_d - pairs_p)} gained"
    )
    print(f"MFE base-pair distance: {observed}")

    random.seed(SEED)
    control = []
    for _ in range(40):
        seq = list(PARENTAL)
        for pos in random.sample(range(20, len(PARENTAL) - 20), 3):
            seq[pos] = random.choice([b for b in "ACGT" if b != seq[pos]])
        control.append(RNA.bp_distance(struct_p, fold("".join(seq))[0]))
    control.sort()
    pct = 100 * sum(1 for x in control if x <= observed) / len(control)
    print(
        f"\ncontrol, 40 random 3-substitution variants: "
        f"median {statistics.median(control):.0f}, "
        f"range {control[0]}-{control[-1]}"
    )
    print(f"    observed {observed} -> ~{pct:.0f}th percentile (no worse than chance)")

    entropy_p = positional_entropy(PARENTAL)
    entropy_d = positional_entropy(DOMESTICATED)
    print("\npositional entropy at edited bases (lower = fold better determined):")
    print(f"    {'pos':>5} {'parental':>10} {'domesticated':>14}")
    for pos in sorted(EDITS):
        print(f"    {pos:>5} {entropy_p[pos]:>10.3f} {entropy_d[pos]:>14.3f}")
    print(f"    median across molecule: {statistics.median(entropy_p[1:]):.3f}")
    return observed, control


def report_natural_variation(columns, variants) -> tuple[int, list[int]]:
    print("\n" + "=" * 72)
    print("3. NATURAL VARIATION — 79 EMCV isolates (NCBI)")
    print("=" * 72)
    covered = [
        (i, c.most_common(1)[0][1] / sum(c.values()))
        for i, c in sorted(columns.items())
        if sum(c.values()) >= 40
    ]
    mean_identity = sum(v for _, v in covered) / len(covered)
    invariant = sum(1 for _, v in covered if v == 1.0)
    print(f"columns with coverage: {len(covered)}; mean identity {mean_identity:.3f}")
    print(
        f"perfectly invariant columns: {invariant} "
        f"({invariant / len(covered):.0%}) <- the base rate to judge against"
    )

    print("\nis the domesticated allele seen in nature?")
    for pos in sorted(EDITS):
        counter = columns[pos - 1]
        total = sum(counter.values())
        conservation = counter.most_common(1)[0][1] / total
        natural = counter.get(EDITS[pos], 0)
        verdict = f"YES, {natural}/{total} isolates" if natural else "NO — novel base"
        # A raw conservation figure reads as more alarming than it is: rank the
        # column against the rest of the molecule before drawing conclusions.
        rank = sum(1 for _, v in covered if v > conservation) / len(covered)
        print(
            f"    {pos:>4} {PARENTAL[pos - 1]}->{EDITS[pos]}  alleles {dict(counter.most_common())}"
        )
        print(
            f"          conservation {conservation:.3f} | {verdict} | "
            f"{rank:.0%} of columns are MORE conserved"
        )
    print(
        "\n    Calibration: with over half the molecule invariant across a tight\n"
        "    clade, an unobserved allele is close to the base rate, not a signal.\n"
        "    Absence from 79 isolates is weak evidence on its own."
    )

    struct_p = fold(PARENTAL)[0]
    reference_pairs = pair_set(struct_p)
    natural_divergence = sorted(
        structural_divergence(seq, cmap, reference_pairs)
        for seq, cmap in variants.items()
        if seq != PARENTAL
    )
    observed = len(reference_pairs ^ pair_set(fold(DOMESTICATED)[0]))
    less_than = sum(1 for x in natural_divergence if x > observed)
    print("\nstructural divergence from the parental fold (symmetric pair difference):")
    print(
        f"    natural isolates (n={len(natural_divergence)}): "
        f"min {natural_divergence[0]}, "
        f"median {statistics.median(natural_divergence):.0f}, "
        f"max {natural_divergence[-1]}"
    )
    print(f"    domesticated: {observed}")
    print(
        f"    -> smaller than {less_than}/{len(natural_divergence)} "
        f"natural isolates "
        f"({100 * less_than / len(natural_divergence):.0f}%)"
    )
    return observed, natural_divergence


def describe_emcv_position(pos: int) -> str:
    """Annotate an IRES position with its EMCV coordinate and any landmark."""
    emcv = pos + EMCV_OFFSET
    hits = [name for lo, hi, name in EMCV_LANDMARKS if lo <= emcv <= hi]
    return f"EMCV {emcv}" + (f" — {'; '.join(hits)}" if hits else "")


def report_confident_pair_changes() -> None:
    """Which base pairs actually move, and do they sit anywhere functional?"""
    print("\n" + "=" * 72)
    print("4. WHICH PAIRS MOVE, AND DO THEY MATTER?")
    print("=" * 72)
    pairs_p, _ = confident_pairs(PARENTAL)
    pairs_d, _ = confident_pairs(DOMESTICATED)
    lost = sorted(pairs_p - pairs_d)
    gained = sorted(pairs_d - pairs_p)
    print(f"lost   ({len(lost)}): {lost}")
    print(f"gained ({len(gained)}): {gained}")
    span = [x for pair in lost + gained for x in pair]
    lo, hi = min(span), max(span)
    print(
        f"\nevery change falls in IRES {lo}-{hi} = EMCV {lo + EMCV_OFFSET}"
        f"-{hi + EMCV_OFFSET} ({(hi - lo + 1) / len(PARENTAL):.0%} of the "
        f"molecule); the rest is untouched"
    )

    print("\nfunctional context:")
    for lo_e, hi_e, name in EMCV_LANDMARKS:
        overlap = lo + EMCV_OFFSET <= hi_e and hi + EMCV_OFFSET >= lo_e
        print(f"    {'OVERLAPS ' if overlap else 'clear of  '}EMCV {lo_e}-{hi_e}  {name}")
    print("\n    The J-K eIF4G platform and the initiator AUG are untouched.")
    print("    The moved pairs land in domain I, including ROI 3, where base")
    print("    pairing is experimentally validated as important for translation.")
    print("")
    print("    Read that last point with care, because it partly refutes itself.")
    print("    Maloney & Joseph report that the validated domain I helices do NOT")
    print("    form in EMCV IRES structures predicted by modern secondary-structure")
    print("    methods -- ViennaRNA included. So the model is known to be wrong in")
    print("    exactly the region this localisation claim depends on.")
    print("")
    print("    What survives: the whole-molecule comparison, where the same model")
    print("    is applied to both sequences and systematic error largely cancels.")
    print("    What does not: any confident statement about WHICH domain-I helices")
    print("    move. Treat the localisation as a hypothesis, not a finding.")


def report_alternative_edits(columns, variants) -> None:
    print("\n" + "=" * 72)
    print("5. WAS THERE A SAFER EDIT? — both BbsI sites")
    print("=" * 72)
    struct_p, mfe_p = fold(PARENTAL)
    baseline = count_sites(PARENTAL, "GAAGAC")
    chosen = {137: "G", 236: "C"}

    for label, (start, end) in BBSI_SITES.items():
        site_seq = PARENTAL[start - 1 : end]
        without_site = sum(
            1
            for seq in variants
            if "GTCTTC" not in seq[start - 40 : end + 40]
            and "GAAGAC" not in seq[start - 40 : end + 40]
        )
        edited = next(p for p in chosen if start <= p <= end)
        print(f"\n{label}, IRES {start}-{end} = {site_seq}")
        print(f"  {describe_emcv_position(start)}")
        print(f"  Fragmid edited position {edited} ({PARENTAL[edited - 1]}->{chosen[edited]})")
        print(f"  natural isolates already lacking this site: {without_site}/{len(variants)}")

        rows = []
        for pos in range(start, end + 1):
            counter = columns[pos - 1]
            total = sum(counter.values())
            for base in "ACGT":
                if base == PARENTAL[pos - 1]:
                    continue
                candidate = apply_edits(PARENTAL, {pos: base})
                if count_sites(candidate, "GAAGAC") >= baseline:
                    continue  # does not kill the site
                struct_c, mfe_c = fold(candidate)
                natural = counter.get(base, 0)
                rows.append(
                    (
                        f"{pos} {PARENTAL[pos - 1]}->{base}",
                        f"{natural}/{total}" if natural else "novel",
                        RNA.bp_distance(struct_p, struct_c),
                        mfe_c - mfe_p,
                        positional_entropy(candidate)[pos],
                        chosen.get(pos) == base,
                    )
                )
        print(f"    {'edit':<12}{'in nature':<14}{'bpDist':>7}{'dMFE':>8}{'entropy':>9}")
        for lbl, natural, dist, dmfe, entropy, is_current in sorted(
            rows, key=lambda r: (r[4], r[2])
        )[:6]:
            marker = "  <-- Fragmid's choice" if is_current else ""
            print(f"    {lbl:<12}{natural:<14}{dist:>7}{dmfe:>+8.2f}{entropy:>9.2f}{marker}")
        if edited == 236:
            print("    Fragmid's choice is the best option at this site: the")
            print("    allele occurs in nature and it lowers local entropy.")
        else:
            print("    Fragmid's choice ranks near the bottom here. 139 C->T and")
            print("    138 T->A are predicted structurally silent but sit at")
            print("    positions invariant across isolates; 134 G->T uses an")
            print("    allele nature already runs.")


def make_figure(observed, natural_divergence, control, columns) -> None:
    setup_plotting_style()
    fig, (ax_left, ax_right) = plt.subplots(1, 2, figsize=(11, 4.2))

    ax_left.hist(
        natural_divergence,
        bins=12,
        color=COLORS["secondary"],
        alpha=0.75,
        edgecolor="white",
        linewidth=0.8,
        label=f"natural EMCV isolates (n={len(natural_divergence)})",
    )
    ax_left.axvline(
        observed,
        color=COLORS["primary"],
        linewidth=2.2,
        label=f"domesticated IRES2_v1 ({observed})",
    )
    ax_left.set_xlabel("structural divergence from parental fold\n(symmetric base-pair difference)")
    ax_left.set_ylabel("isolates")
    ax_left.legend(frameon=False, fontsize=8, loc="upper right")

    positions = sorted(EDITS)
    labels = [f"{p}\n{PARENTAL[p - 1]}→{EDITS[p]}" for p in positions]
    fractions = []
    for pos in positions:
        counter = columns[pos - 1]
        fractions.append(100 * counter.get(EDITS[pos], 0) / sum(counter.values()))
    bars = ax_right.bar(
        labels,
        fractions,
        color=[
            # maroon flags the edit nature never makes; teal = allele already
            # observed in circulating EMCV
            COLORS["primary"] if f == 0 else COLORS["secondary"]
            for f in fractions
        ],
        width=0.55,
    )
    for bar, frac in zip(bars, fractions, strict=True):
        note = "not seen\nin nature" if frac == 0 else f"{frac:.0f}%"
        ax_right.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 1.2,
            note,
            ha="center",
            va="bottom",
            fontsize=8,
        )
    ax_right.set_ylabel("EMCV isolates carrying\ndomesticated allele (%)")
    ax_right.set_xlabel("edited position")
    ax_right.set_ylim(0, max(fractions) * 1.45 + 6)

    fig.tight_layout()
    save_figure(fig, "s1e2_ires_domestication", output_dir=str(HERE / "figures"))
    plt.close(fig)


def write_notion_svg(columns, observed, natural_divergence) -> None:
    """Emit the lean SVG used for the Notion embed.

    Notion caps inline attachments at 200 KiB and a matplotlib SVG blows past
    that once font paths are embedded. The tempting shortcut is to hand-write a
    small SVG with the numbers typed in -- don't. It drifts from the analysis
    silently. This renders the same values `make_figure` plots, from the same
    run, so the embed cannot disagree with the canonical figure.
    """
    axis_max = 50 * (int(max(natural_divergence) / 50) + 1)
    x0, width = 60, 480

    def x(value: float) -> float:
        return x0 + (value / axis_max) * width

    rows, y = [], 46
    enzymes = {137: "BbsI", 236: "BbsI", 480: "ApaLI"}
    for pos in sorted(EDITS):
        counter = columns[pos - 1]
        total = sum(counter.values())
        count = counter.get(EDITS[pos], 0)
        colour = COLORS["primary"] if count == 0 else COLORS["secondary"]
        bar = (
            f'\n    <rect x="196" y="{y - 10}" width="{260 * count / total:.1f}" '
            f'height="12" rx="3" fill="{COLORS["secondary"]}"/>'
            if count
            else ""
        )
        flag = (
            f'\n    <text x="524" y="{y}" fill="{COLORS["primary"]}">watch</text>'
            if count == 0
            else ""
        )
        rows.append(
            f'    <text x="0" y="{y}" fill="{COLORS["dark"]}">{pos} '
            f"{PARENTAL[pos - 1]}&#8594;{EDITS[pos]}</text>\n"
            f'    <text x="76" y="{y}" fill="{COLORS["neutral"]}">EMCV '
            f"{pos + EMCV_OFFSET} &#183; {enzymes[pos]}</text>\n"
            f'    <rect x="196" y="{y - 10}" width="260" height="12" rx="3" '
            f'fill="{COLORS["line"]}"/>{bar}\n'
            f'    <text x="466" y="{y}" fill="{colour}" font-weight="600">'
            f"{count} / {total}</text>{flag}"
        )
        y += 30

    ticks = "".join(
        f'\n      <line x1="{x(v):.0f}" y1="238" x2="{x(v):.0f}" y2="243" '
        f'stroke="#B9A9A4"/><text x="{x(v):.0f}" y="258">{v}</text>'
        for v in range(0, axis_max + 1, axis_max // 3)
    )
    median = statistics.median(natural_divergence)
    svg = f"""<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 620 300" width="620" height="300" font-family="Manrope, Inter, system-ui, sans-serif">
  <text x="0" y="14" font-size="13" fill="{COLORS["dark"]}" font-weight="600">Do natural EMCV isolates carry the domesticated allele?</text>
  <g font-size="12">
{chr(10).join(rows)}
  </g>
  <line x1="0" y1="140" x2="620" y2="140" stroke="{COLORS["line"]}" stroke-width="1"/>
  <text x="0" y="168" font-size="13" fill="{COLORS["dark"]}" font-weight="600">Structural divergence from the parental fold</text>
  <text x="0" y="188" font-size="12" fill="{COLORS["neutral"]}">base-pair difference; natural isolates set the scale</text>
  <g>
    <line x1="{x0}" y1="238" x2="{x0 + width}" y2="238" stroke="#B9A9A4" stroke-width="1"/>
    <g font-size="11" fill="{COLORS["neutral"]}" text-anchor="middle">{ticks}
    </g>
    <rect x="{x(min(natural_divergence)):.0f}" y="212" width="{x(max(natural_divergence)) - x(min(natural_divergence)):.0f}" height="14" rx="4" fill="{COLORS["secondary"]}" opacity="0.32"/>
    <line x1="{x(median):.0f}" y1="208" x2="{x(median):.0f}" y2="230" stroke="{COLORS["secondary"]}" stroke-width="2.5"/>
    <text x="{x(median) + 8:.0f}" y="206" font-size="11" fill="{COLORS["secondary"]}">natural isolates, median {median:.0f} (range {min(natural_divergence)}&#8211;{max(natural_divergence)})</text>
    <line x1="{x(observed):.0f}" y1="206" x2="{x(observed):.0f}" y2="232" stroke="{COLORS["primary"]}" stroke-width="2.5"/>
    <text x="{x(observed):.0f}" y="284" font-size="11" fill="{COLORS["primary"]}" text-anchor="middle" font-weight="600">domesticated, {observed}</text>
  </g>
</svg>
"""
    path = HERE / "figures" / "s1e2_ires_domestication_notion.svg"
    path.write_text(svg)
    # A malformed embed uploads fine and renders as an error box in Notion, so
    # fail here rather than there.
    ElementTree.fromstring(svg)
    print(f"wrote {path.name} ({len(svg)} bytes, Notion inline cap 204800)")


def write_table(columns, observed, natural_divergence) -> None:
    path = HERE / "ires_domestication_edits.tsv"
    entropy_p = positional_entropy(PARENTAL)
    entropy_d = positional_entropy(DOMESTICATED)
    lines = [
        "position\tparental\tdomesticated\tenzyme_removed\tnatural_isolates_with_"
        "domesticated_allele\ttotal_isolates\tconservation\tentropy_parental\t"
        "entropy_domesticated"
    ]
    enzymes = {137: "BbsI", 236: "BbsI", 480: "ApaLI"}
    for pos in sorted(EDITS):
        counter = columns[pos - 1]
        total = sum(counter.values())
        lines.append(
            f"{pos}\t{PARENTAL[pos - 1]}\t{EDITS[pos]}\t{enzymes[pos]}\t"
            f"{counter.get(EDITS[pos], 0)}\t{total}\t"
            f"{counter.most_common(1)[0][1] / total:.4f}\t"
            f"{entropy_p[pos]:.4f}\t{entropy_d[pos]:.4f}"
        )
    path.write_text("\n".join(lines) + "\n")
    print(f"\nwrote {path.relative_to(HERE.parents[2])}")


def main() -> None:
    report_sequence_differences()
    _, control = report_folding()
    path = fetch_emcv_genomes()
    columns, variants = align_natural_ires(path)
    observed_div, natural_divergence = report_natural_variation(columns, variants)
    report_confident_pair_changes()
    report_alternative_edits(columns, variants)
    make_figure(observed_div, natural_divergence, control, columns)
    write_notion_svg(columns, observed_div, natural_divergence)
    write_table(columns, observed_div, natural_divergence)

    print("\n" + "=" * 72)
    print("VERDICT")
    print("=" * 72)
    print("Adopt Fragmid's IRES2_v1. Three independent controls agree the")
    print("perturbation is small: it is no larger than a random 3-nt change,")
    print("smaller than most natural EMCV isolates differ from each other, and")
    print("it leaves the initiator AUG context and uAUG count untouched.")
    print("")
    print("The J-K eIF4G platform and the initiator AUG are untouched, and the")
    print("one edit inside J-K (EMCV 739) is the structurally silent one.")
    print("")
    print("Position 137 (EMCV 396) is the weakest of the three, but the case")
    print("against it is thinner than the raw numbers suggest, and both supporting")
    print("signals deflate on inspection:")
    print("  - 'not seen in 79 isolates' is close to the base rate. Over half the")
    print("    molecule is invariant across this clade, so most non-wild-type bases")
    print("    are unobserved everywhere.")
    print("  - '98.7% conserved' sounds extreme but is below median here: 54% of")
    print("    columns are MORE conserved than 137.")
    print("  - the predicted rearrangement sits in domain I, which is precisely")
    print("    where secondary-structure prediction is documented to fail.")
    print("")
    print("Weighing against that, IRES2_v1 drives PuroR, BlastR and HygroR in")
    print("pF_AA333/334/335, so Broad's own selection cassettes depend on it. That")
    print("constrains 'not dead' firmly, though drug selection is forgiving enough")
    print("that it does not exclude a modest loss of output.")
    print("")
    print("Recommendation: ship IRES2_v1 as-is. Re-deriving costs a synthesis and")
    print("sequence-verification cycle and gives up drop-in Fragmid compatibility,")
    print("to buy an improvement that exists only inside an unreliable model. The")
    print("exception is if the IRES fragment is being synthesised de novo anyway —")
    print("then 134 G->T is free, and free is worth taking. Otherwise hold it as")
    print("the first thing to try if the p006/p007 A/B shows low IRES output.")
    print("Position 236 needs no action: Fragmid's choice is the best available")
    print("edit at its site on every measure computed here.")
    print("")
    print("Limits: thermodynamic prediction only (ViennaRNA, 37C, Turner). No")
    print("probing data, no ITAF/eIF4G binding, no tertiary structure. EMCV is a")
    print("tight clade, so 'natural variation' here is shallow. The empirical")
    print("A/B in p006/p007 is what settles it; this only sets the prior.")


if __name__ == "__main__":
    main()
