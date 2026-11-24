"""What the per-pair barcode-swap rate means for L1 as one pooled lentiviral library.

``seqtools.barcode_swap_rate`` (Walton et al., bioRxiv 2024.03.17.585235, PMID 38558968,
PMC10979941) gives the fraction of integrants whose two markers are discordant when the
two genomes in a virion differ and the sequence between the markers is identical between
them: ``Y(x) = 1/2 (1 - exp(-2kx))``, ``k = 5.3e-4`` per bp. Walton calibrated it on
designs where the intervening sequence is constant vector, shared by *every* library
member, so any heterozygous virion can swap. In L1 a long identical tract is shared only by
a specific partner (two constructs of one CDS, an isoform and its reference), and a
template switch needs that partner in the same virion. How often that happens is the
partner's share of the pool, about 1 in 3,255, not 1.

The model, per ordered pair (upstream template ``a``, barcode template ``b``):

* Virions pair genomes at random from the pool (Walton: "essentially randomly paired");
  the partner of ``b`` is ``a`` with probability ``p_a``, its molar share of the pool.
  ``p`` follows the length-weighted pool of record (3.5-fold per kb, ``l1_order_plan.tsv``):
  pmol is proportional to ng / vector length.
* Every maximal identical block of >= 20 nt between ``a`` and ``b`` upstream of the barcode
  is a place a switch can land seamlessly; a switch there happens with probability
  ``Y(L)`` (first order: blocks add). The product is ``a`` up to the switch and ``b`` from it,
  so it carries ``b``'s barcode.
* The product is classified by translating it: ``correct`` (``b``'s protein), ``mislabel``
  (``a``'s intact protein under ``b``'s barcode) or ``chimera`` (neither: a fusion, a
  truncation, a frameshift, or a tag lost or doubled). Where the block lies in ``b`` is
  recorded: 5' constant head, tag cassette, native CDS, or the 3' tail (stop, CC, handle).

``impurity[b] = sum_a p_a * sum_blocks Y(L) * [product != b's protein]`` is the expected
fraction of cells carrying ``b``'s barcode that express something else, with every virion
heterozygous (the upper bound; diluting the library plasmid with carrier lowers it).

Specific pairs (sharing a 24-mer seen in at most ``RARE`` inserts) are computed exactly;
for each barcode the rest of the pool is sampled (``PARTNERS`` partners drawn by pool share)
and scaled by the share it represents.

Inputs: the recoded order of record (pinned by SHA-256 in ``audit_domestication_barcodes``),
``reviewer-package/l1_order_plan.tsv``. Writes ``l1_swap_pool_model.tsv`` (per insert),
``l1_swap_pool_pairs.tsv`` (same-CDS pairs, and other pairs sharing > 50 nt before the
barcode), ``l1_swap_pool_model_summary.json`` and ``reviewer-package/figures/l1_swap_pool.png``.

Written by the Claude Code session l1-swap-risk (2026-11-20) to answer Becca's 2026-11-20
comment on section 1 of the brief ("doesn't make sense if we are pooling at low per
construct rates"). Report only: it changes no design.

    uv run --locked python l1_swap_pool_model.py
"""

from __future__ import annotations

import csv
import json
import os
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import audit_domestication_barcodes as audit
import numpy as np
from build_constructs import CASSETTE_AA

from deliverome_analysis import seqtools

HERE = Path(__file__).resolve().parent
PLAN = HERE / "reviewer-package" / "l1_order_plan.tsv"
FIGURES = HERE / "reviewer-package" / "figures"
OUT_TSV = HERE / "l1_swap_pool_model.tsv"
OUT_PAIRS = HERE / "l1_swap_pool_pairs.tsv"
OUT_JSON = HERE / "l1_swap_pool_model_summary.json"

WEIGHT = "x3.5"  # the pool of record: 3.5-fold per kb (BRIEF.md section 8)
MIN_BLOCK = 20  # shortest identical block counted as switch substrate (Walton's shortest: 19 bp)
K = 24  # k-mer for finding specific pairs
RARE = 20  # a 24-mer in at most this many inserts marks a specific pair
PARTNERS = 96  # sampled background partners per barcode
SEED = 2026_11_20
GATE = 100
PAIR_ROW_NT = 50
LEADER_AA = (
    80  # a tag at or before this residue is an N-terminal tag (SP and Igκ leaders <= ~60 aa)
)
HEAD = audit.WRAP_5  # every upstream region starts after the 5' Esp3I wrapper

# ----------------------------------------------------------------------------- inputs


def load() -> list[dict]:
    rec, digest = audit.recoded_inserts(None)
    if digest != audit.RECODED_SHA256:
        raise SystemExit("order of record changed; re-pin it in audit_domestication_barcodes")
    plan = {r["arm_id"]: r for r in csv.DictReader(PLAN.open(), delimiter="\t")}
    out = []
    for r in rec:
        p = plan[r["key"]]
        up = audit.upstream_of_barcode(r["insert"])
        coding, err = audit.coding_of(r["insert"])
        assert up is not None and not err, r["key"]
        protein = seqtools.translate(coding)
        cas = protein.find(CASSETTE_AA)
        start = len(audit.KOZAK)  # ATG offset in the upstream region
        out.append(
            {
                "arm_id": r["key"],
                "insert_id": p["insert_id"],
                "design_id": p["design_id"],
                "rung": p["rung"],
                "set": p["set"],
                "up": up,
                "protein": protein,
                # region boundaries in upstream coordinates
                "atg": start,
                "cas0": start + 3 * cas,
                "cas1": start + 3 * (cas + len(CASSETTE_AA)),
                "cas_aa": cas,
                "stop_end": start + len(coding),
                "pmol": float(p[f"pool_ng_per_pool_{WEIGHT}"]) / float(p["vector_bp"]),
                "insert_nt": len(r["insert"]),
            }
        )
    total = sum(i["pmol"] for i in out)
    for i in out:
        i["p"] = i["pmol"] / total
    return out


# ----------------------------------------------------------------------------- pair model


def blocks(a: str, b: str, kmers_a: dict[str, list[int]]) -> list[tuple[int, int, int]]:
    """Maximal identical blocks (i in a, j in b, length) of at least MIN_BLOCK nt."""
    out = []
    for j in range(len(b) - MIN_BLOCK + 1):
        for i in kmers_a.get(b[j : j + MIN_BLOCK], ()):
            if i > 0 and j > 0 and a[i - 1] == b[j - 1]:
                continue  # inside a block already found from its left end
            n = MIN_BLOCK
            while i + n < len(a) and j + n < len(b) and a[i + n] == b[j + n]:
                n += 1
            out.append((i, j, n))
    return out


def kmer_index(s: str) -> dict[str, list[int]]:
    idx: dict[str, list[int]] = defaultdict(list)
    for i in range(len(s) - MIN_BLOCK + 1):
        idx[s[i : i + MIN_BLOCK]].append(i)
    return idx


def region(b: dict, j: int) -> str:
    if j < b["atg"]:
        return "5' head"
    if b["cas0"] <= j < b["cas1"]:
        return "tag cassette"
    if j < b["stop_end"] - 3:
        return "native CDS"
    return "3' tail"


def product_protein(a: dict, b: dict, i: int, j: int, n: int) -> str:
    m = n // 2
    prod = a["up"][: i + m] + b["up"][j + m :]
    body = prod[len(audit.KOZAK) :]
    if not body.startswith("ATG"):
        return ""  # the start codon itself was lost
    # A frameshifted product can run into the barcode slot without a stop; translate leniently
    # (codons with N read as X) rather than with the audit's strict check.
    aa = []
    for k in range(0, len(body) - 2, 3):
        c = seqtools.GENETIC_CODE.get(body[k : k + 3], "X")
        aa.append(c)
        if c == "*":
            break
    return "".join(aa)


def leader_swap(prot: str, b: dict) -> bool:
    """b's protein from its tag onward, behind a different short leader (<= LEADER_AA residues).

    The product of a switch inside the shared tag cassette when both members carry the tag near
    their N-terminus: b's mature protein, tagged, with a's signal peptide or Igκ leader. The
    leader is cleaved, so the displayed protein is b's; only the leader's efficiency differs.
    """
    tail = b["protein"][b["cas_aa"] :]
    return b["cas_aa"] <= LEADER_AA and prot.endswith(tail) and len(prot) - len(tail) <= LEADER_AA


def pair_terms(a: dict, b: dict, idx_a: dict | None = None) -> list[tuple[str, str, int, float]]:
    """(outcome, region in b, block length, Y) for every block, product = a upstream, b barcode."""
    idx_a = idx_a if idx_a is not None else kmer_index(a["up"])
    out = []
    for i, j, n in blocks(a["up"], b["up"], idx_a):
        prot = product_protein(a, b, i, j, n)
        if prot == b["protein"]:
            outcome = "correct"
        elif prot == a["protein"]:
            outcome = "mislabel"
        elif leader_swap(prot, b):
            outcome = "leader swap"
        else:
            outcome = "chimera"
        out.append((outcome, region(b, j + n // 2), n, seqtools.barcode_swap_rate(n)))
    return out


common_suffix = seqtools.common_suffix


# ----------------------------------------------------------------------------- workers

_INS: list[dict] = []
_IDX: list[dict] = []


def _init(ins):
    global _INS, _IDX
    _INS = ins
    _IDX = [None] * len(ins)


def _idx(k: int):
    if _IDX[k] is None:
        _IDX[k] = kmer_index(_INS[k]["up"])
    return _IDX[k]


def _background(job):
    """Sampled background for barcode b: per-outcome and per-region expected impurity."""
    b, partners = job
    acc = Counter()
    for a in partners:
        for outcome, reg, _n, y in pair_terms(_INS[a], _INS[b], _idx(a)):
            acc[(outcome, reg)] += y
    return b, {f"{o}|{r}": v / len(partners) for (o, r), v in acc.items()}


def _specific(job):
    a, b = job
    return a, b, pair_terms(_INS[a], _INS[b], _idx(a))


# ----------------------------------------------------------------------------- main


def specific_pairs(ins: list[dict]) -> set[tuple[int, int]]:
    owners: dict[str, set[int]] = defaultdict(set)
    for k, i in enumerate(ins):
        u = i["up"]
        for s in range(len(u) - K + 1):
            owners[u[s : s + K]].add(k)
    pairs = set()
    for ks in owners.values():
        if 1 < len(ks) <= RARE:
            ks = sorted(ks)
            pairs.update((x, y) for n, x in enumerate(ks) for y in ks[n + 1 :])
    return pairs


def _write_tsv(path: Path, rows: list[dict]) -> None:
    """``rows`` as a tab-separated table with the first row's keys as the header."""
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, rows[0].keys(), delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def main() -> None:
    ins = load()
    n = len(ins)
    p = np.array([i["p"] for i in ins])
    spec = specific_pairs(ins)
    partners_of: dict[int, set[int]] = defaultdict(set)
    for x, y in spec:
        partners_of[x].add(y)
        partners_of[y].add(x)

    rng = np.random.default_rng(SEED)
    jobs = []
    for b in range(n):
        excl = partners_of[b] | {b}
        q = p.copy()
        q[list(excl)] = 0
        jobs.append((b, rng.choice(n, size=PARTNERS, p=q / q.sum()).tolist()))
    bg_share = np.array([1 - p[list(partners_of[b] | {b})].sum() for b in range(n)])

    workers = max(1, (os.cpu_count() or 2) - 1)
    with ProcessPoolExecutor(workers, initializer=_init, initargs=(ins,)) as ex:
        background = dict(ex.map(_background, jobs, chunksize=16))
        ordered = [(x, y) for x, y in spec] + [(y, x) for x, y in spec]
        specific = list(ex.map(_specific, ordered, chunksize=64))

    # per barcode
    rows = []
    spec_imp = np.zeros(n)
    spec_mis = np.zeros(n)
    spec_chi = np.zeros(n)
    spec_lead = np.zeros(n)
    region_total = Counter()
    pair_rows = []
    for a, b, terms in specific:
        mis = sum(y for o, _r, _n, y in terms if o == "mislabel")
        chi = sum(y for o, _r, _n, y in terms if o == "chimera")
        spec_mis[b] += p[a] * mis
        spec_chi[b] += p[a] * chi
        spec_lead[b] += p[a] * sum(y for o, _r, _n, y in terms if o == "leader swap")
        for o, r, _n, y in terms:
            if o != "correct":
                region_total[("specific", o, r)] += p[b] * p[a] * y
        A, B = ins[a], ins[b]
        tract = common_suffix(A["up"], B["up"])
        # keep the table under the repo's 5 MB limit: every same-CDS pair, and other pairs only
        # where they share more than PAIR_ROW_NT before the barcode (the rest share a linker)
        if a < b and (A["design_id"] == B["design_id"] or tract > PAIR_ROW_NT):
            pair_rows.append(
                {
                    "arm_a": A["arm_id"],
                    "arm_b": B["arm_id"],
                    "relation": "same CDS" if A["design_id"] == B["design_id"] else "different CDS",
                    "set_a": A["set"],
                    "set_b": B["set"],
                    "suffix_tract_nt": tract,
                    # the suffix runs back from the barcode: 3' tail (stop, CC, handle) first,
                    # then native CDS until the two constructs diverge
                    "suffix_in_native_cds_nt": max(0, tract - (len(B["up"]) - (B["stop_end"] - 3))),
                    "per_pair_rate_walton": round(seqtools.barcode_swap_rate(tract), 4),
                    "co_packaging_probability": round(2 * p[a] * p[b], 9),
                    # of the cells carrying b's barcode, the fraction carrying a's construct
                    "mislabel_of_b_cells_from_a": round(
                        p[a] * seqtools.barcode_swap_rate(tract), 9
                    ),
                    "pool_expected_swaps_per_integrant": round(
                        2 * p[a] * p[b] * seqtools.barcode_swap_rate(tract), 11
                    ),
                }
            )
    spec_imp = spec_mis + spec_chi

    bg_mis = np.zeros(n)
    bg_chi = np.zeros(n)
    bg_lead = np.zeros(n)
    for b in range(n):
        for key, v in background[b].items():
            o, r = key.split("|")
            if o == "mislabel":
                bg_mis[b] += bg_share[b] * v
            elif o == "chimera":
                bg_chi[b] += bg_share[b] * v
            elif o == "leader swap":
                bg_lead[b] += bg_share[b] * v
            if o != "correct":
                region_total[("background", o, r)] += p[b] * bg_share[b] * v

    for k, i in enumerate(ins):
        rows.append(
            {
                "insert_id": i["insert_id"],
                "arm_id": i["arm_id"],
                "design_id": i["design_id"],
                "rung": i["rung"],
                "set": i["set"],
                "pool_share": f"{p[k]:.3e}",
                "specific_partners": len(partners_of[k]),
                "specific_mislabel": f"{spec_mis[k]:.3e}",
                "specific_chimera": f"{spec_chi[k]:.3e}",
                "background_mislabel": f"{bg_mis[k]:.4f}",
                "background_chimera": f"{bg_chi[k]:.4f}",
                "leader_swap": f"{spec_lead[k] + bg_lead[k]:.4f}",
                "impurity": f"{spec_imp[k] + bg_mis[k] + bg_chi[k]:.4f}",
            }
        )
    _write_tsv(OUT_TSV, rows)
    _write_tsv(OUT_PAIRS, pair_rows)

    same = [r for r in pair_rows if r["relation"] == "same CDS"]
    same_rate = np.array([r["per_pair_rate_walton"] for r in same])
    same_pool = np.array([r["pool_expected_swaps_per_integrant"] for r in same])
    imp = spec_imp + bg_mis + bg_chi
    by_rung = defaultdict(list)
    for k, i in enumerate(ins):
        by_rung[i["rung"]].append(k)

    def wq(x):  # pool-weighted library fraction
        return float((p * x).sum())

    summary = {
        "source": "Walton, Qin & Blainey, bioRxiv 2024.03.17.585235 (PMID 38558968): k=5.3e-4/bp",
        "inserts": n,
        "pool_weighting": WEIGHT,
        "pool_share_median": float(np.median(p)),
        "pool_share_1_over": round(1 / float(np.median(p))),
        "min_block_nt": MIN_BLOCK,
        "background_partners_sampled_per_barcode": PARTNERS,
        "specific_pairs": len(spec),
        "same_cds_pairs": {
            "pairs": len(same),
            "suffix_tract_median_nt": float(np.median([r["suffix_tract_nt"] for r in same])),
            "over_gate": int(sum(r["suffix_tract_nt"] > GATE for r in same)),
            "per_pair_rate_median": float(np.median(same_rate)),
            "per_pair_rate_max": float(same_rate.max()),
            "co_packaging_probability_median": float(
                np.median([r["co_packaging_probability"] for r in same])
            ),
            "library_fraction_all_791_pairs": float(same_pool.sum()),
            "suffix_tract_nt_in_native_cds_median_fraction": float(
                np.median(
                    [
                        r["suffix_in_native_cds_nt"] / r["suffix_tract_nt"]
                        for r in same
                        if r["suffix_tract_nt"]
                    ]
                )
            ),
        },
        "different_cds_specific_pairs": {
            "pairs_over_50_nt_written": len(pair_rows) - len(same),
            "over_gate": int(
                sum(r["suffix_tract_nt"] > GATE for r in pair_rows if r["relation"] != "same CDS")
            ),
        },
        "library_fraction_wrong_protein": {
            "specific_mislabel": wq(spec_mis),
            "specific_chimera": wq(spec_chi),
            "background_mislabel": wq(bg_mis),
            "background_chimera": wq(bg_chi),
            "total": wq(imp),
            "not_counted_leader_swap": wq(spec_lead + bg_lead),
        },
        "per_barcode_impurity": {
            "median": float(np.median(imp)),
            "max": float(imp.max()),
            "specific_part_median": float(np.median(spec_imp)),
            "specific_part_max": float(spec_imp.max()),
            "by_rung_median": {
                r: {
                    "inserts": len(ks),
                    "impurity": float(np.median(imp[ks])),
                    "mislabel": float(np.median(bg_mis[ks] + spec_mis[ks])),
                    "chimera": float(np.median(bg_chi[ks] + spec_chi[ks])),
                    "leader_swap": float(np.median(bg_lead[ks] + spec_lead[ks])),
                }
                for r, ks in by_rung.items()
            },
        },
        # Long constructs carry more specific-partner impurity (BRIEF.md section 2): their
        # siblings and orthologs are long too, and the length-weighted pool gives them more
        # molecules. Quartiles are by insert length; the rank correlation is Spearman's.
        "specific_impurity_by_insert_length": length_summary(ins, spec_imp),
        "wrong_protein_by_block_region": {
            f"{src}|{o}|{r}": v for (src, o, r), v in sorted(region_total.items())
        },
        "rate_at_gate": seqtools.barcode_swap_rate(GATE),
        "rate_at_cassette_153": seqtools.barcode_swap_rate(153),
    }
    OUT_JSON.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    plot(ins, p, same, spec_imp, bg_mis, bg_chi, bg_lead + spec_lead, summary)


def _ranks(x: np.ndarray) -> np.ndarray:
    """Average ranks (ties share the mean rank), for a Spearman correlation without scipy."""
    order = np.argsort(x, kind="mergesort")
    ranks = np.empty(len(x))
    ranks[order] = np.arange(len(x))
    for v in np.unique(x):
        tied = x == v
        ranks[tied] = ranks[tied].mean()
    return ranks


def length_summary(ins: list[dict], spec_imp: np.ndarray, cut: float = 0.01) -> dict:
    """Specific-partner impurity by insert-length quartile, its Spearman correlation with
    length, and the barcodes above ``cut`` (with their length range)."""
    length = np.array([i["insert_nt"] for i in ins], dtype=float)
    q = np.quantile(length, [0.25, 0.5, 0.75])
    quart = np.searchsorted(q, length, side="right")  # 0 = shortest quarter
    over = np.flatnonzero(spec_imp > cut)
    return {
        "median_by_length_quartile": {
            str(k + 1): float(np.median(spec_imp[quart == k])) for k in range(4)
        },
        "quartile_upper_bounds_nt": [float(v) for v in q],
        "spearman_vs_insert_length": float(np.corrcoef(_ranks(length), _ranks(spec_imp))[0, 1]),
        "barcodes_over_cut": {
            "cut": cut,
            "count": len(over),
            "max": float(spec_imp.max()),
            "insert_nt_range": [int(length[over].min()), int(length[over].max())]
            if len(over)
            else [],
            "design_ids": sorted({ins[k]["design_id"] for k in over}),
        },
    }


def plot(ins, p, same, spec_imp, bg_mis, bg_chi, lead, summary) -> None:
    import matplotlib.pyplot as plt

    from deliverome_analysis.plotting import get_color, save_figure, setup_plotting_style

    setup_plotting_style()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 8))

    # a: the per-pair rate across the 791 same-CDS pairs, and the same pairs weighted by pool
    rate = np.array([r["per_pair_rate_walton"] for r in same])
    pool = np.array([r["mislabel_of_b_cells_from_a"] for r in same])
    bins = np.logspace(-9, 0, 46)
    ax1.hist(
        rate, bins=bins, color=get_color("primary"), alpha=0.85, label="per pair, both co-packaged"
    )
    ax1.hist(
        pool,
        bins=bins,
        color=get_color("secondary"),
        alpha=0.85,
        label="share of one barcode's cells, in the pool",
    )
    gate = summary["rate_at_gate"]
    ax1.axvline(gate, color=get_color("neutral"), ls="--", lw=1.5, label="100-nt gate")
    ax1.set_xscale("log")
    ax1.set_xlabel("modelled swap rate (Walton)")
    ax1.set_ylabel("same-CDS construct pairs")
    ax1.legend(frameon=False, loc="upper left")
    ax1.text(-0.12, 1.02, "a", transform=ax1.transAxes, fontweight="bold")

    # b: per-barcode impurity, by source, every insert
    series = [
        ("specific partner (same CDS, isoform)", spec_imp, get_color("primary")),
        ("any partner: shared 3' tail (mislabel)", bg_mis, get_color("accent")),
        ("any partner: shared tag cassette, fusion or truncation", bg_chi, get_color("secondary")),
        (
            "any partner: shared tag cassette, leader swap only (N-terminal tags)",
            lead[lead > 0],
            get_color("quaternary"),
        ),
    ]
    for label, x, c in series:
        xs = np.sort(np.clip(x, 1e-6, None))
        ax2.plot(xs, np.arange(1, len(xs) + 1) / len(xs), color=c, lw=2.2, label=label)
    ax2.axvline(gate, color=get_color("neutral"), ls="--", lw=1.5, label="100-nt gate")
    ax2.set_xscale("log")
    ax2.set_xlim(1e-6, 1)
    ax2.set_xlabel("expected fraction of a barcode's cells carrying another construct")
    ax2.set_ylabel("cumulative fraction of inserts")
    ax2.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.2), fontsize="x-small")
    ax2.text(-0.12, 1.02, "b", transform=ax2.transAxes, fontweight="bold")
    fig.tight_layout()
    save_figure(fig, "l1_swap_pool", output_dir=FIGURES)


def replot() -> None:
    """Redraw the figure from the written tables, without rerunning the model."""
    summary = json.loads(OUT_JSON.read_text())
    same = [
        {
            k: float(v) if k in ("per_pair_rate_walton", "mislabel_of_b_cells_from_a") else v
            for k, v in r.items()
        }
        for r in csv.DictReader(OUT_PAIRS.open(), delimiter="\t")
        if r["relation"] == "same CDS"
    ]
    t = list(csv.DictReader(OUT_TSV.open(), delimiter="\t"))
    col = lambda k: np.array([float(r[k]) for r in t])  # noqa: E731
    plot(
        None,
        None,
        same,
        col("specific_mislabel") + col("specific_chimera"),
        col("background_mislabel"),
        col("background_chimera"),
        col("leader_swap"),
        summary,
    )


if __name__ == "__main__":
    import sys

    replot() if "--replot" in sys.argv else main()
