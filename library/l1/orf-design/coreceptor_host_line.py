"""Co-receptor gate against the screening line: what a K562 (or lower-median) rule would change.

The S1E3 co-receptor gate (``final_library.coreceptor_ok``) keeps a gene whose co-receptor
dependency is ``required`` only if every listed obligate partner has a median TPM >= 5
across the 1,206 HPA v23 cell lines. That asks whether the partner is ubiquitous. The
question raised in review is a different one: is the partner present in the line we screen
in first (K562)? This script answers it without changing the gate.

It reruns the S1E3 filter offline (``library_pipeline.run_gene_filter``, the same rerun the
pipeline verifies against the 2,370-gene list), reads the pinned HPA v23 cell-line TPM
matrix, and for every gene the partner-TPM rule decides it reports each partner's median,
K562 and best T-, B- and myeloid-lineage TPM. Partners are resolved through
``gene_resolver`` (MyGene, the catalog's own config) and joined to HPA on Ensembl gene ID;
aliases of one gene (``CLR`` / ``CALCRL``) collapse to one partner, and free-text partners
that do not resolve are flagged rather than guessed.

Partners are also scored line by line over the S2E18 screening shortlist (its 42 candidate lines,
29 of them in HPA v23) and the 10 of those in house or on order.

Outputs (experiment root, mirrored to processed://sprint4/s4e9-orf-design-fg-panel/):

- ``coreceptor_host_line_genes.tsv``: one row per rule-decided gene.
- ``coreceptor_host_line_partners.tsv``: one row per (gene, listed partner string).
- ``coreceptor_host_line_options.tsv``: each candidate rule, the genes it adds and drops
  relative to the current gate, among genes that pass every other S1E3 gate.
- ``figures/coreceptor_host_line.png``.

    uv run --locked python coreceptor_host_line.py
"""

from __future__ import annotations

import hashlib
import io
import json
import os
import re
import sys
import tempfile
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

from deliverome_analysis import gene_resolver
from deliverome_analysis.plotting import create_figure, get_color, save_figure, setup_plotting_style
from deliverome_analysis.s3 import external

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))

# HPA 23.0 (2023-06-19, Ensembl 109). The same bytes are also at external://hpa/v23/, pinned
# again on 2026-11-20 before this earlier pin was noticed.
HPA_KEY = "human-protein-atlas/2023-05/rna_celline.tsv.zip"
HPA_SHA256 = "e5b11779308c2720cfe89d9cb19932ad4c0c3fea0fa37e8c61a14d501995852a"
DESC_KEY = "hpa/v23/rna_celline_description.tsv.zip"
DESC_SHA256 = "93696a4e703b8da4bf4f27699ded3ddb170787549c1b29b80a9bc543f3162389"
MYGENE_CACHE = HERE / "coreceptor_host_line_mygene"  # committed: the resolver's answer
CACHE = HERE / "data"  # gitignored local cache of the parsed TPM matrix
EXPECTED_GENES = 2370  # library_pipeline.verify_gene_list: the filter's output

SCREEN_LINE = "K-562"
GATE_TPM = 5.0  # final_library.coreceptor_ok: median cell-line TPM >= 5

# Lineage from HPA's own "Disease subtype". Lines HPA leaves blank are assigned from
# Cellosaurus; myeloma (plasma cell), HAP1 and K562 itself are kept out of all three.
T_SUBTYPE = re.compile(r"T-cell|ATL")
B_SUBTYPE = re.compile(r"B-cell|DLBCL|Burkitt")
MYELOID_SUBTYPE = re.compile(r"Acute Myelogenous")  # CML blast-crisis lines can be lymphoid
BLANK_SUBTYPE = {
    "JURKAT": "T",
    "Jurkat E6.1": "T",
    "MOLT-4": "T",
    "BJAB": "B",
    "CCRF-SB": "B",
    "Daudi": "B",
    "REH": "B",
    "TMD8": "B",
    "U-2932": "B",
    "U-698-M": "B",
    "WSU-FSCCL": "B",
    "HEL": "myeloid",
    "HL-60": "myeloid",
    "HMC-1": "myeloid",
    "MOLM-13": "myeloid",
    "NB4": "myeloid",
    "THP-1": "myeloid",
    "U-937": "myeloid",
    "UKE-1": "myeloid",
}

# Screening shortlist: the candidate lines of the cell-line panel experiment (S2E18,
# panel.candidates()). HPA names are matched after normalization; these four are the same
# Cellosaurus line under another HPA name (HEK293T is read from parental HEK293, as S2E18 does).
PANEL_DIR = HERE.parents[1] / "sprint2/s2e18-cell-line-panel-selection"
HPA_ALIASES = {
    "HEK293T": "HEK293",  # parental line; HPA has no HEK293T
    "BJ-5ta": "BJ1-hTERT",  # CVCL_6573
    "HEP3B": "Hep 3B2.1-7",  # CVCL_0326
    "SUM159": "SUM159PT",  # CVCL_5423
}
# Lines in house or on order, per the S2E18 README "Result" (2026-09-08): 3 in house plus the
# 9-line order. HK-2 and GM12878 are not in HPA v23.
IN_HAND = ["K562", "HEK293T", "HeLa", "U937", "NALM6", "SH-SY5Y", "A549", "HuH-7", "hTERT-RPE1",
           "U2OS", "HK-2", "GM12878"]  # fmt: skip
SHORTLIST_MIN_LINES = 2
# The adopted rule: every required partner (one per group of alternatives) at the ARCHS4 CPM
# equivalent of >= 10 TPM in >= 2 of the 42 shortlist lines, ARCHS4 scoring every line; no K562
# clause. Becca chose 10 TPM / 2 lines on 2026-11-21 (then HPA for 29 lines + ARCHS4 for 13) and
# "ARCHS4 for all" on 2026-11-22 (brief comment). Implemented by the co-receptor host-line agent
# session.
ADOPTED_TPM = 10
ADOPTED = "ADOPTED: ARCHS4 for all 42 lines, >= 10 CPM in >= 2"
# Becca, 2026-11-22: "let's do 10 cpm". The ARCHS4 cut for the 10 TPM rule is set to 10 CPM rather
# than the kappa-optimal 15, because kappa is flat from 10 to 20 CPM (0.732 / 0.758 / 0.755) and
# the 15 CPM optimum dropped the TYROBP-paired receptors on U937's 14.4 CPM. The calibration table
# is still written as the evidence for that choice.
ADOPTED_CPM = 10.0
NEAR_CUT_FACTOR = 2  # flag an addition whose second-best line is under 2x the CPM cut
# ARCHS4 v2.5 per-line consensus built by the cell-line panel experiment (log2(CPM + 1), gene
# symbols x panel lines): covers all 42 shortlist lines, including the 13 HPA lacks. CPM is not
# TPM, so the CPM cut matching each TPM cut is calibrated on the 29 lines both sources share.
# Pinned append-only copy (PROVENANCE.md beside it) of processed://scratch/becca/cell-line-selection/.
ARCHS4_KEY = "archs4/v2.5/cell-line-panel-consensus/2026-09-04/rnaseq_consensus_archs4.parquet"
ARCHS4_SHA256 = "a3a789acdd51c440d536c52e61ec440ddd5d3b5e8f24f3a185d7703b78e56ee1"
CPM_GRID = [
    1,
    2,
    3,
    5,
    7,
    10,
    15,
    20,
    30,
    40,
    50,
    75,
    100,
]  # "even 2 of the shortlist lines" (Becca, 2026-11-21)

# Genes whose listed partners are alternatives, not all required: the gene needs one partner
# from each group. Proposed from the partner biology for the genes where "all" and "any"
# disagree; every other gene keeps the gate's all-partners reading. Pending review.
ALTERNATIVES = {
    "SLC16A1": [["BSG", "EMB"]],  # MCT1 pairs with basigin or embigin
    "SLC16A8": [["BSG", "EMB"]],
    "ST14": [["SPINT1", "SPINT2"]],  # matriptase with HAI-1 or HAI-2
    "GRIN1": [["GRIN2A", "GRIN2B", "GRIN2C", "GRIN2D", "GRIN3A", "GRIN3B"]],  # any GluN2/3
    "KCND2": [["KCNIP1", "KCNIP2", "DPP6", "DPP10"]],  # Kv4.2 with a KChIP or a DPP
    "GABRA5": [["GABRB1", "GABRB2", "GABRB3"], ["GABRG1", "GABRG2", "GABRG3"]],
    "GABRB1": [
        ["GABRA1", "GABRA2", "GABRA3", "GABRA4", "GABRA5", "GABRA6"],
        ["GABRG1", "GABRG2", "GABRG3", "GABRD"],
    ],
}

# Genes the review asked about by name; shown in panel b whether or not they pass the
# other gates.
NAMED = ["CD3D", "CD3E", "CD3G", "CD247", "PTCRA", "CD79A", "CD79B", "ITGAL", "ITGAM", "ITGAX"]


def sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def fetch(key: str, digest: str, dest: Path) -> Path:
    path = dest / Path(key).name
    if not path.exists() or sha256(path) != digest:
        external.download(key, path)
    assert sha256(path) == digest, f"{key}: SHA-256 does not match the pinned file"
    return path


def read_zip_tsv(path: Path, **kw) -> pd.DataFrame:
    with zipfile.ZipFile(path) as z:
        (name,) = z.namelist()
        with z.open(name) as fh:
            return pd.read_csv(io.TextIOWrapper(fh), sep="\t", **kw)


def hpa_tpm() -> pd.DataFrame:
    """Gene (Ensembl, unversioned) x cell line TPM, from the pinned v23 file."""
    CACHE.mkdir(exist_ok=True)
    pq = CACHE / f"hpa_rna_celline_{HPA_SHA256[:12]}.parquet"
    if not pq.exists():
        z = fetch(HPA_KEY, HPA_SHA256, CACHE)
        long = read_zip_tsv(z, usecols=["Gene", "Cell line", "TPM"])
        long.pivot(index="Gene", columns="Cell line", values="TPM").to_parquet(pq)
    return pd.read_parquet(pq)


def lineages(lines: list[str]) -> dict[str, list[str]]:
    desc = read_zip_tsv(fetch(DESC_KEY, DESC_SHA256, CACHE), dtype=str).fillna("")
    assert set(desc["Cell line"]) == set(lines), "description and TPM files list different lines"
    out = {"T": [], "B": [], "myeloid": []}
    for _, r in desc.iterrows():
        line, sub = r["Cell line"], r["Disease subtype"]
        if line == SCREEN_LINE or r["Disease"] not in ("Leukemia", "Lymphoma"):
            continue
        lin = BLANK_SUBTYPE.get(line) if not sub else None
        lin = lin or ("T" if T_SUBTYPE.search(sub) else None)
        lin = lin or ("B" if B_SUBTYPE.search(sub) else None)
        lin = lin or ("myeloid" if MYELOID_SUBTYPE.search(sub) else None)
        if lin:
            out[lin].append(line)
    return out


def shortlist(lines: list[str], names: list[str] | None = None) -> tuple[list[str], list[str]]:
    """Panel lines present in HPA v23 (as HPA names), and the panel lines it lacks."""
    sys.path.insert(0, str(PANEL_DIR))
    import panel

    norm = {re.sub(r"[^a-z0-9]", "", x.lower()): x for x in lines}
    have, missing = [], []
    for name in names or panel.candidates():
        hit = HPA_ALIASES.get(name) or norm.get(re.sub(r"[^a-z0-9]", "", name.lower()))
        (have if hit in lines else missing).append(hit if hit in lines else name)
    return have, missing


def archs4_cpm() -> pd.DataFrame:
    path = CACHE / f"archs4_consensus_{ARCHS4_SHA256[:12]}.parquet"
    if not path.exists() or sha256(path) != ARCHS4_SHA256:
        external.download(ARCHS4_KEY, path)
    assert sha256(path) == ARCHS4_SHA256, "ARCHS4 consensus differs from the pinned file"
    return 2 ** pd.read_parquet(path) - 1


def hpa_symbols() -> pd.Series:
    """HPA's own Ensembl -> gene name pairing (unique names only), to join symbol-keyed ARCHS4."""
    path = CACHE / f"hpa_ensembl_symbol_{HPA_SHA256[:12]}.tsv"
    if not path.exists():
        z = fetch(HPA_KEY, HPA_SHA256, CACHE)
        with zipfile.ZipFile(z) as zf, zf.open(zf.namelist()[0]) as fh:
            chunks = pd.read_csv(
                io.TextIOWrapper(fh), sep="\t", usecols=["Gene", "Gene name"], chunksize=2_000_000
            )
            pd.concat(c.drop_duplicates() for c in chunks).drop_duplicates().to_csv(
                path, sep="\t", index=False
            )
    m = pd.read_csv(path, sep="\t")
    return m[~m["Gene name"].duplicated(keep=False)].set_index("Gene")["Gene name"]


def calibrate(tpm, cpm, sym, pairs) -> pd.DataFrame:
    """Agreement of HPA TPM >= t with ARCHS4 CPM >= c over the shared lines, all genes."""
    genes = [e for e in tpm.index if e in sym.index and sym[e] in cpm.index]
    x = tpm.loc[genes, [h for h, _ in pairs]].to_numpy()
    y = cpm.loc[[sym[e] for e in genes], [a for _, a in pairs]].to_numpy()
    rows = []
    for t in (GATE_TPM, 10.0):
        hx = x >= t
        for c in CPM_GRID:
            hy = y >= c
            tp, fp = (hx & hy).sum(), (~hx & hy).sum()
            fn, tn = (hx & ~hy).sum(), (~hx & ~hy).sum()
            po = (tp + tn) / hx.size
            pe = ((tp + fp) * (tp + fn) + (fn + tn) * (fp + tn)) / hx.size**2
            rows.append({"tpm_cut": t, "cpm_cut": c, "agreement": po, "kappa": (po - pe) / (1 - pe),
                         "sensitivity": tp / (tp + fn), "specificity": tn / (tn + fp),
                         "genes": len(genes), "lines": len(pairs)})  # fmt: skip
    return pd.DataFrame(rows)


def lines_ok_sym(cpm, groups: list[list[str]], lines: list[str], thr: float):
    """lines_ok on ARCHS4 symbols; None when a requirement has no partner ARCHS4 measures."""
    groups = [[x for x in grp if x in cpm.index] for grp in groups]
    if not all(groups):
        return None
    ok = pd.Series(True, index=lines)
    for grp in groups:
        ok &= (cpm.loc[grp, lines] >= thr).any(axis=0)
    return ok


def lines_ok(tpm: pd.DataFrame, reqs: list[list[str]], lines: list[str], thr: float) -> pd.Series:
    """Per line: every requirement group has at least one partner at >= thr TPM."""
    ok = pd.Series(True, index=lines)
    for grp in reqs:
        ok &= (tpm.loc[grp, lines] >= thr).any(axis=0)
    return ok


def rerun_gene_filter():
    """The S1E3 filter, offline; returns the imported final_library module."""
    import library_pipeline

    work = Path(tempfile.mkdtemp(prefix="coreceptor-host-line-"))
    passed = library_pipeline.run_gene_filter(work)
    assert len(passed) == EXPECTED_GENES, f"filter rerun gave {len(passed)} genes"
    os.chdir(work / "s1e3")
    import final_library as fl

    return fl, passed


def others_pass(fl, r) -> bool:
    """Every S1E3 gate except the co-receptor one (as in final_library's membership)."""
    g, f = r["gene"], r["filters"]
    return (
        (fl.B.passes_canonical(f) or fl.B.passes_likely(f))
        and f.get("surface_accessibility") != "no"
        and f.get("surface_specificity") in ("surface_dominant", "mixed")
        and fl.anchored(r)
        and fl.expr_ok(g, f)
        and fl.localization_ok(r)
        and fl.orf_ok(r)
        and fl.ecd_ok(r)
    )


def rule_decided(fl, r) -> bool:
    """The gene's co-receptor outcome comes from the partner-TPM rule, not a review."""
    g, f = r["gene"], r["filters"]
    partners = fl.part.get(g, {}).get("partners", [])
    reviewed = fl.reviewed_partner_ok(
        fl.cohort_accessions.get(g), partners, fl.reviewed_partner_median
    )
    return (
        reviewed is None
        and g not in fl.CORECEPTOR_ALLOWLIST
        and f.get("co_receptor_dependency") == "required"
    )


def _rule_options(genes: pd.DataFrame) -> dict:
    """Each candidate co-receptor rule, by label, as a per-gene pass mask."""

    def rule(col, thr):
        return genes["self_only"] | (genes[col] >= thr)

    options = {
        "current gate (all partners, median >= 5, gate's own partner lookup)": genes[
            "coreceptor_ok_now"
        ],
        "all partners, median >= 5, corrected partner lookup": rule("all_median", 5),
        "all partners, median >= 2": rule("all_median", 2),
        "all partners, median >= 1": rule("all_median", 1),
        "all partners, median >= 10": rule("all_median", 10),
        "all partners, K562 >= 5": rule("all_k562", 5),
        "all partners, K562 >= 1": rule("all_k562", 1),
        "all partners, median >= 5 OR K562 >= 5": rule("all_median", 5) | rule("all_k562", 5),
        "curated alternatives, median >= 5": rule("curated_median", 5),
        "curated alternatives, K562 >= 5": rule("curated_k562", 5),
        "curated alternatives, median >= 5 OR K562 >= 5": rule("curated_median", 5)
        | rule("curated_k562", 5),
        "all partners, >= 5 TPM in >= 2 shortlist lines": genes["all_shortlist_lines_ge5"]
        >= SHORTLIST_MIN_LINES,
        "curated alternatives, >= 5 TPM in >= 1 shortlist line": genes[
            "curated_shortlist_lines_ge5"
        ]
        >= 1,
        "curated alternatives, >= 5 TPM in >= 2 shortlist lines": genes[
            "curated_shortlist_lines_ge5"
        ]
        >= SHORTLIST_MIN_LINES,
        "curated alternatives, >= 10 TPM in >= 2 shortlist lines": genes[
            "curated_shortlist_lines_ge10"
        ]
        >= SHORTLIST_MIN_LINES,
        "curated alternatives, >= 5 TPM in >= 2 shortlist lines OR in K562": (
            genes["curated_shortlist_lines_ge5"] >= SHORTLIST_MIN_LINES
        )
        | rule("curated_k562", 5),
        "HPA 29 + ARCHS4 13, >= 5 TPM (or CPM equivalent) in >= 2 of 42 lines OR in K562": (
            genes["combined_lines_ge5"] >= SHORTLIST_MIN_LINES
        )
        | rule("curated_k562", 5),
        ADOPTED: genes["self_only"]
        | (genes[f"archs4_lines_ge{ADOPTED_TPM}"] >= SHORTLIST_MIN_LINES),
        "HPA 29 + ARCHS4 13, >= 10 TPM (ARCHS4 lines >= 10 CPM) in >= 2 of 42 lines "
        "(adopted 2026-11-21, replaced 2026-11-22)": genes["self_only"]
        | (genes[f"combined_lines_ge{ADOPTED_TPM}"] >= SHORTLIST_MIN_LINES),
        "HPA 29 + ARCHS4 13, >= 10 TPM (ARCHS4 lines >= 10 CPM) in >= 2 of 42 lines OR in K562": (
            genes["combined_lines_ge10"] >= SHORTLIST_MIN_LINES
        )
        | rule("curated_k562", 10),
        "ARCHS4 alone, CPM equivalent of >= 5 TPM in >= 2 of 42 lines OR in K562": (
            genes["archs4_lines_ge5"] >= SHORTLIST_MIN_LINES
        )
        | (genes["archs4_k562_ok"] == True),  # noqa: E712
        "ARCHS4 alone, >= 15 CPM (kappa-optimal) in >= 2 of 42 lines "
        "(adopted 2026-11-22, replaced the same day by 10 CPM)": genes["self_only"]
        | (genes["archs4_lines_ge15cpm"] >= SHORTLIST_MIN_LINES),
        "curated alternatives, >= 10 TPM in >= 2 shortlist lines OR in K562": (
            genes["curated_shortlist_lines_ge10"] >= SHORTLIST_MIN_LINES
        )
        | rule("curated_k562", 10),
        "curated alternatives, >= 5 TPM in >= 2 in-house/ordered lines": genes[
            "curated_in_hand_lines_ge5"
        ]
        >= SHORTLIST_MIN_LINES,
        "curated alternatives, median >= 5 OR >= 5 TPM in >= 2 shortlist lines": rule(
            "curated_median", 5
        )
        | (genes["curated_shortlist_lines_ge5"] >= SHORTLIST_MIN_LINES),
        "any partner, median >= 5": rule("any_median", 5),
        "any partner, K562 >= 5": rule("any_k562", 5),
        "any partner, median >= 5 OR K562 >= 5": rule("any_median", 5) | rule("any_k562", 5),
    }
    return options


def main() -> int:
    setup_plotting_style()
    fl, passed = rerun_gene_filter()
    tpm = hpa_tpm()
    assert tpm.shape[1] == 1206, tpm.shape
    med = tpm.median(axis=1)
    # The gate reads a re-derived per-gene median; it must be this file's median.
    gate_med = pd.Series(fl.med_ens)
    shared = gate_med.index.intersection(med.index)
    assert len(shared) == len(gate_med), "gate medians list genes the v23 file lacks"
    assert np.allclose(gate_med[shared], med[shared].round(1), atol=0.051), "medians differ"
    lin = lineages(list(tpm.columns))

    decided = [r for r in fl.recs if rule_decided(fl, r)]
    tokens = sorted(
        {
            fl._partner_sym(p)
            for r in decided
            for p in fl.part.get(r["gene"], {}).get("partners", [])
        }
    )
    # The gate's own map of misresolved free-text partners (final_library.PARTNER_SYMBOL_FIXES);
    # each target is resolved here too, so its Ensembl ID still comes from the resolver.
    fixes = fl.PARTNER_SYMBOL_FIXES
    res = gene_resolver.resolve(set(tokens) | set(fixes.values()), cache_dir=MYGENE_CACHE)
    for tok, target in fixes.items():
        assert res[target]["symbol"] == target, f"{target} does not resolve to itself"
        res[tok] = {**res[target], "override": target}
    short, short_missing = shortlist(list(tpm.columns))
    assert SCREEN_LINE in short and len(short) == 29, (len(short), short_missing)
    sys.path.insert(0, str(PANEL_DIR))
    import panel

    cpm = archs4_cpm()
    assert set(panel.candidates()) <= set(cpm.columns), "ARCHS4 consensus lacks a shortlist line"
    hpa_of = {}
    for name in panel.candidates():
        got, _ = shortlist(list(tpm.columns), [name])
        hpa_of[name] = got[0] if got else None
    pairs = [(h, a) for a, h in hpa_of.items() if h]  # (HPA name, ARCHS4 name), 29 lines
    archs4_only = [a for a, h in hpa_of.items() if not h]  # the 13 HPA lacks
    assert len(pairs) == 29 and len(archs4_only) == 13
    hsym = hpa_symbols()
    calib = calibrate(tpm, cpm, hsym, pairs)
    best = calib.loc[calib.groupby("tpm_cut")["kappa"].idxmax()].set_index("tpm_cut")["cpm_cut"]
    cpm_cut_calibrated = {5: float(best[GATE_TPM]), 10: float(best[10.0])}
    cpm_cut = {5: cpm_cut_calibrated[5], 10: ADOPTED_CPM}
    hand, hand_missing = shortlist(list(tpm.columns), IN_HAND)
    assert set(hand) <= set(short) and len(hand) == 10, (hand, hand_missing)

    prow, grow = [], []
    line_pass: dict[int, dict[str, pd.Series]] = {5: {}, 10: {}}  # HPA 29 + ARCHS4 13
    line_pass_a4: dict[int, dict[str, pd.Series]] = {5: {}, 10: {}}  # ARCHS4, all 42 lines
    a4_score: dict[str, pd.Series] = {}  # ARCHS4 CPM of the limiting requirement, per line
    panel_of_hpa = dict(pairs)
    for r in decided:
        g, f = r["gene"], r["filters"]
        listed = fl.part.get(g, {}).get("partners", [])
        by_ens: dict[str, dict] = {}
        for p in listed:
            tok = fl._partner_sym(p)
            hit = res.get(tok, {})
            ens = [e for e in hit.get("ensembl", []) if e in tpm.index]
            self_p = fl._is_self_partner(g, p)
            status = (
                "self"
                if self_p
                else "unresolved"
                if hit.get("notfound") or not hit.get("ensembl")
                else "ambiguous"
                if hit.get("ambiguous") or len(ens) > 1
                else "not_in_hpa"
                if not ens
                else "resolved"
            )
            e = ens[0] if status == "resolved" else None
            row = {
                "gene": g,
                "listed_partner": p,
                "token": tok,
                "resolved_symbol": hit.get("symbol"),
                "ensembl": e,
                "status": status,
                "symbol_fixed_to": fixes.get(tok),
                "gate_median_tpm": fl.partner_median(p),  # what the gate itself saw
                "median_tpm": None if e is None else round(float(med[e]), 1),
                "k562_tpm": None if e is None else float(tpm.loc[e, SCREEN_LINE]),
                "lines_ge5_tpm": None if e is None else int((tpm.loc[e] >= GATE_TPM).sum()),
                "shortlist_lines_ge5_tpm": None
                if e is None
                else int((tpm.loc[e, short] >= GATE_TPM).sum()),
            }
            for k, cols in lin.items():
                if e is not None:
                    best = tpm.loc[e, cols]
                    row[f"best_{k}_tpm"] = float(best.max())
                    row[f"best_{k}_line"] = best.idxmax()
                    row[f"{k}_lines_ge5_tpm"] = int((best >= GATE_TPM).sum())
            gm = row["gate_median_tpm"]
            row["gate_lookup_differs"] = (gm is not None or row["median_tpm"] is not None) and (
                gm is None or row["median_tpm"] is None or abs(gm - row["median_tpm"]) > 0.051
            )
            prow.append(row)
            if e is not None:
                by_ens.setdefault(e, row)
        non_self = [p for p in listed if not fl._is_self_partner(g, p)]
        dist = pd.DataFrame(by_ens.values())
        unresolved = sorted(
            {
                x["listed_partner"]
                for x in prow
                if x["gene"] == g and x["status"] not in ("resolved", "self")
            }
        )

        self_only = bool(listed) and not non_self  # homo-oligomer: the gate passes it

        def lim(col, how, dist=dist, self_only=self_only):
            if self_only or not len(dist):
                return np.nan
            return float(getattr(dist[col], how)())

        out = {
            "gene": g,
            "uniprot": fl.cohort_accessions.get(g),
            "passes_other_gates": others_pass(fl, r),
            "in_library_now": g in passed,
            "coreceptor_ok_now": fl.coreceptor_ok(g, f),
            "self_only": self_only,
            "partners_listed": "; ".join(map(str, listed)),
            "partners_resolved": "; ".join(dist["resolved_symbol"]) if len(dist) else "",
            "partners_unresolved": "; ".join(map(str, unresolved)),
        }
        for col, name in [
            ("median_tpm", "median"),
            ("k562_tpm", "k562"),
            ("best_T_tpm", "best_T"),
            ("best_B_tpm", "best_B"),
            ("best_myeloid_tpm", "best_myeloid"),
        ]:
            out[f"all_{name}"] = lim(col, "min")  # the limiting partner
            out[f"any_{name}"] = lim(col, "max")  # the best partner
        for col, name in [("median_tpm", "median"), ("k562_tpm", "k562")]:
            if g in ALTERNATIVES and len(dist):
                v = dist.set_index("resolved_symbol")[col]
                groups = [[x for x in grp if x in v.index] for grp in ALTERNATIVES[g]]
                assert all(groups), f"{g}: an alternatives group has no resolved partner"
                out[f"curated_{name}"] = float(min(v[grp].max() for grp in groups))
            else:
                out[f"curated_{name}"] = out[f"all_{name}"]
        # Line-by-line view over the shortlist: the gene passes in a line when each requirement
        # (a partner, or one member of a curated alternatives group) clears the threshold there.
        ens_of = (
            dict(zip(dist["resolved_symbol"], dist["ensembl"], strict=True)) if len(dist) else {}
        )
        reqs_all = [[e] for e in ens_of.values()]
        reqs_cur = (
            [[ens_of[x] for x in grp if x in ens_of] for grp in ALTERNATIVES[g]]
            if g in ALTERNATIVES
            else reqs_all
        )
        for tag, reqs in (("all", reqs_all), ("curated", reqs_cur)):
            for thr in (GATE_TPM, 10.0):
                n = (
                    len(short)
                    if self_only
                    else int(lines_ok(tpm, reqs, short, thr).sum())
                    if reqs
                    else 0
                )
                out[f"{tag}_shortlist_lines_ge{int(thr)}"] = n
        out["curated_in_hand_lines_ge5"] = (
            len(hand)
            if self_only
            else int(lines_ok(tpm, reqs_cur, hand, GATE_TPM).sum())
            if reqs_cur
            else 0
        )

        # HPA for its 29 lines plus ARCHS4 (calibrated CPM cut) for the 13 it lacks; and ARCHS4
        # alone over all 42.
        # ARCHS4 rows are symbols: reach them through the partner's Ensembl ID and HPA v23's name
        # for it (CDC50A is TMEM30A there), falling back to the resolver's symbol.
        def a4(x, ens_of=ens_of):
            name = hsym.get(ens_of[x])
            return name if name in cpm.index else x

        sym_groups = (
            [[a4(x) for x in grp if x in ens_of] for grp in ALTERNATIVES[g]]
            if g in ALTERNATIVES
            else [[a4(x)] for x in ens_of]
        )
        out["archs4_missing_partners"] = "; ".join(
            x for grp in sym_groups for x in grp if x not in cpm.index
        )
        for thr in (5, 10):
            if self_only:
                out[f"combined_lines_ge{thr}"] = len(pairs) + len(archs4_only)
                out[f"archs4_lines_ge{thr}"] = len(hpa_of)
                out["archs4_lines_ge15cpm"] = len(hpa_of)
                continue
            a13 = lines_ok_sym(cpm, sym_groups, archs4_only, cpm_cut[thr]) if sym_groups else None
            h29 = (
                lines_ok(tpm, reqs_cur, [h for h, _ in pairs], float(thr))
                if reqs_cur
                else pd.Series(False, index=[h for h, _ in pairs])
            )
            out[f"combined_lines_ge{thr}"] = int(h29.sum()) + (0 if a13 is None else int(a13.sum()))
            s13 = a13 if a13 is not None else pd.Series(False, index=archs4_only)
            line_pass[thr][g] = pd.concat([h29.rename(index=panel_of_hpa), s13])[list(hpa_of)]
            if thr == 5:
                out["archs4_only_lines_ge5_names"] = (
                    "" if a13 is None else ", ".join(a13[a13].index)
                )
            a42 = lines_ok_sym(cpm, sym_groups, list(hpa_of), cpm_cut[thr]) if sym_groups else None
            if thr == 10:  # the kappa-optimal cut, kept so the 15 CPM decision reproduces
                a15 = (
                    lines_ok_sym(cpm, sym_groups, list(hpa_of), cpm_cut_calibrated[10])
                    if sym_groups
                    else None
                )
                out["archs4_lines_ge15cpm"] = np.nan if a15 is None else int(a15.sum())
            line_pass_a4[thr][g] = a42 if a42 is not None else pd.Series(False, index=list(hpa_of))
            if sym_groups and all(x in cpm.index for grp in sym_groups for x in grp):
                a4_score[g] = pd.concat(
                    [cpm.loc[grp, list(hpa_of)].max(axis=0) for grp in sym_groups], axis=1
                ).min(axis=1)
            out[f"archs4_lines_ge{thr}"] = np.nan if a42 is None else int(a42.sum())
            if thr == 5:
                out["archs4_k562_ok"] = None if a42 is None else bool(a42["K562"])
        out["curated_shortlist_lines_ge5_names"] = (
            ", ".join(lines_ok(tpm, reqs_cur, short, GATE_TPM).pipe(lambda x: x[x].index))
            if reqs_cur and not self_only
            else ""
        )
        if len(dist):
            limiting = dist.loc[dist["k562_tpm"].idxmin()]
            out["k562_limiting_partner"] = limiting["resolved_symbol"]
        grow.append(out)

    partners = pd.DataFrame(prow)
    genes = pd.DataFrame(grow).sort_values("gene").reset_index(drop=True)

    options = _rule_options(genes)
    cand = genes["passes_other_gates"].to_numpy()
    now = genes["coreceptor_ok_now"].to_numpy()
    # A gene can be in the library through another route while failing this gate (CALHM3,
    # CATSPERE, CATSPERG): an addition is a gene not yet in the library; a drop is a gene the
    # gate keeps today and the rule would not.
    in_lib = genes["in_library_now"].to_numpy()
    names = genes["gene"].to_numpy()
    orows = []
    for label, ok in options.items():
        ok = ok.to_numpy()
        genes[f"ok: {label}"] = ok
        orows.append(
            {
                "rule": label,
                "genes_kept": int((ok & cand).sum()),
                "added": int((ok & ~in_lib & cand).sum()),
                "dropped": int((~ok & now & cand).sum()),
                "added_genes": ", ".join(names[ok & ~in_lib & cand]),
                "dropped_genes": ", ".join(names[~ok & now & cand]),
            }
        )
    opts = pd.DataFrame(orows)

    genes.to_csv(
        HERE / "coreceptor_host_line_genes.tsv", sep="\t", index=False, float_format="%.1f"
    )
    partners.to_csv(
        HERE / "coreceptor_host_line_partners.tsv", sep="\t", index=False, float_format="%.1f"
    )
    opts.to_csv(HERE / "coreceptor_host_line_options.tsv", sep="\t", index=False)
    # Which K562 + 2-line set covers the most genes: a gene is covered when every partner it needs
    # is present (>= thr TPM, or the CPM equivalent) in at least one of the three lines. Genes
    # that pass every other gate and have a partner other than themselves.
    cand_genes = set(genes.loc[genes["passes_other_gates"] & ~genes["self_only"], "gene"])

    def as_mats(lp):
        return {thr: pd.DataFrame(v).T.loc[sorted(cand_genes & set(v))] for thr, v in lp.items()}

    mats = as_mats(line_pass_a4)  # the adopted source; trios and additions use it
    long = pd.concat(
        [m.astype(int).assign(source=src, threshold_tpm=thr).rename_axis("gene").reset_index()
         for src, ms in (("archs4_all_42", mats), ("hpa29_archs4_13", as_mats(line_pass)))
         for thr, m in ms.items()]
    )  # fmt: skip
    long.to_csv(HERE / "coreceptor_host_line_line_pass.tsv", sep="\t", index=False)
    trows = []
    others = [x for x in hpa_of if x != "K562"]
    for i, l2 in enumerate(others):
        for l3 in others[i + 1 :]:
            row = {
                "line2": l2,
                "line3": l3,
                "both_in_house_or_ordered": l2 in IN_HAND and l3 in IN_HAND,
            }
            for thr, m in mats.items():
                k = m["K562"]
                cov = k | m[l2] | m[l3]
                row[f"covered_ge{thr}"] = int(cov.sum())
                row[f"beyond_k562_ge{thr}"] = int((cov & ~k).sum())
                row[f"added_genes_ge{thr}"] = ", ".join(m.index[cov & ~k])
            trows.append(row)
    trios = pd.DataFrame(trows).sort_values(
        ["covered_ge10", "covered_ge5", "line2", "line3"], ascending=[False, False, True, True]
    )
    trios.insert(0, "genes_scored", len(mats[10]))
    trios.insert(1, "k562_alone_ge10", int(mats[10]["K562"].sum()))
    trios.to_csv(HERE / "coreceptor_host_line_trios.tsv", sep="\t", index=False)

    # Genes the adopted rule adds to the library, in the gene_list_additions.tsv shape plus the
    # evidence columns. Proposed only: the gene list changes when these rows are appended there.
    m10 = mats[ADOPTED_TPM]
    adds = genes[genes["passes_other_gates"] & genes[f"ok: {ADOPTED}"] & ~genes["in_library_now"]]
    arows = []
    for _, r in adds.iterrows():
        lines = list(m10.columns[m10.loc[r["gene"]]])
        arows.append({
            "gene": r["gene"], "uniprot": r["uniprot"], "added": "", "decided_by": "Becca Carlson",
            "source": "co-receptor gate restated 2026-11-22 (coreceptor_host_line.py; ARCHS4 v2.5 "
                      "cell-line panel consensus, all 42 shortlist lines)",
            "evidence": f"passes every other S1E3 gate; partners {r['partners_resolved']} are >= "
                        f"{cpm_cut[ADOPTED_TPM]:.0f} CPM in {len(lines)} of 42 "
                        f"shortlist lines: {', '.join(lines)}; median-rule limiting partner "
                        f"{r['curated_median']:.1f} TPM",
            "n_lines": len(lines), "in_house_or_ordered": ", ".join(x for x in lines if x in IN_HAND),
            "displayable_in_k562": "K562" in lines,
            "second_line_cpm": round(float(a4_score[r["gene"]].nlargest(2).iloc[-1]), 1),
            "near_cut": bool(a4_score[r["gene"]].nlargest(2).iloc[-1]
                             < NEAR_CUT_FACTOR * cpm_cut[ADOPTED_TPM]),
            "hla_locus": r["gene"].startswith("HLA-"),
        })  # fmt: skip
    pd.DataFrame(arows).sort_values(["n_lines", "gene"], ascending=[False, True]).to_csv(
        HERE / "coreceptor_rule_additions_proposed.tsv", sep="\t", index=False
    )
    calib.to_csv(HERE / "coreceptor_host_line_archs4_calibration.tsv", sep="\t", index=False,
                 float_format="%.3f")  # fmt: skip
    meta = {
        "archs4": f"external://{ARCHS4_KEY}",
        "archs4_sha256": ARCHS4_SHA256,
        "archs4_cpm_cut_calibrated_for_tpm": cpm_cut_calibrated,
        "archs4_cpm_cut_used_for_tpm": cpm_cut,
        "archs4_only_lines": archs4_only,
        "hpa": f"external://{HPA_KEY}",
        "hpa_sha256": HPA_SHA256,
        "cell_lines": int(tpm.shape[1]),
        "lineage_lines": {k: len(v) for k, v in lin.items()},
        "shortlist_lines_in_hpa": short,
        "shortlist_lines_not_in_hpa": short_missing,
        "in_hand_lines_in_hpa": hand,
        "in_hand_lines_not_in_hpa": hand_missing,
        "rule_decided_genes": len(genes),
        "rule_decided_passing_other_gates": int(cand.sum()),
        "of_those_in_library_now": int((cand & now).sum()),
        "partner_strings": len(partners),
        "partner_status": partners["status"].value_counts().to_dict(),
    }
    (HERE / "coreceptor_host_line_summary.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(json.dumps(meta, indent=2))
    print(opts[["rule", "genes_kept", "added", "dropped"]].to_string(index=False))
    plot(genes, mats[ADOPTED_TPM], a4_score, ADOPTED)
    return 0


def plot(genes: pd.DataFrame, m10: pd.DataFrame, a4_score: dict, adopted: str) -> None:
    """Figure on the adopted rule: ARCHS4 CPM over the 42 shortlist lines, 10 CPM in >= 2 lines."""
    cut = ADOPTED_CPM

    def lg(v):
        return np.log10(np.asarray(v, dtype=float) + 0.1)

    fig, (a, b) = create_figure(1, 2, figsize=(30, 16), gridspec_kw={"width_ratios": [1, 1.35]})
    c = genes[genes["passes_other_gates"] & ~genes["self_only"]].copy()
    c = c[c["gene"].isin(list(a4_score))]
    # The rule's test value: the limiting requirement's CPM in the gene's second-best line. A gene
    # passes exactly when this is >= the cut.
    c["second"] = [float(a4_score[g].nlargest(2).iloc[-1]) for g in c["gene"]]
    c["k562"] = [float(a4_score[g]["K562"]) for g in c["gene"]]
    ok = c[f"ok: {adopted}"].to_numpy()
    lib = c["in_library_now"].to_numpy()
    status = np.where(lib, "in", np.where(ok, "add", "out"))
    style = {
        "in": ("in L1 now", get_color("secondary")),
        "add": ("added by the rule", get_color("accent")),
        "out": ("excluded", get_color("neutral")),
    }
    jit = np.random.default_rng(0).uniform(-0.3, 0.3, len(c))
    x, y = lg(c["second"]), c["archs4_lines_ge10"].to_numpy() + jit
    for key, (lab, col) in style.items():
        m = status == key
        a.scatter(
            x[m],
            y[m],
            s=80,
            color=col,
            alpha=0.85,
            edgecolor="white",
            linewidth=0.6,
            label=f"{lab} (n = {int(m.sum())})",
            zorder=3,
        )
    k = c["k562"].to_numpy() >= cut
    a.scatter(
        x[k],
        y[k],
        s=200,
        facecolor="none",
        edgecolor=get_color("dark"),
        linewidth=1.3,
        label=f"partners >= {cut:.0f} CPM in K562 (n = {int(k.sum())})",
        zorder=4,
    )
    a.axvline(lg(cut), color=get_color("dark"), lw=1.5)
    a.text(
        lg(cut) + 0.04,
        41,
        f"rule: every partner >= {cut:.0f} CPM\nin >= {SHORTLIST_MIN_LINES} lines",
        fontsize=14,
        va="top",
    )
    # In L1 through another route although this gate fails them, labeled one by one; the
    # antigen-receptor subunits (partners only in Jurkat) share one label.
    for n, (_, r) in enumerate(c[lib & ~c["coreceptor_ok_now"].to_numpy()].iterrows()):
        a.annotate(f"{r['gene']} (in L1)", (lg(r["second"]), r["archs4_lines_ge10"]),
                   fontsize=12, xytext=(-15 if r["gene"] == "CATSPERE" else 30, 70 + 30 * n), textcoords="offset points",
                   arrowprops={"arrowstyle": "-", "color": get_color("secondary"), "lw": 0.8})  # fmt: skip
    tcr = c[c["gene"].isin(["CD3D", "CD3E", "PTCRA", "TRB"]).to_numpy() & (status == "out")]
    if len(tcr):
        a.annotate(", ".join(sorted(tcr["gene"])) + " (Jurkat only)",
                   (lg(tcr["second"]).mean(), tcr["archs4_lines_ge10"].iloc[0]), fontsize=12,
                   xytext=(120, 190), textcoords="offset points", ha="center",
                   arrowprops={"arrowstyle": "-", "color": get_color("dark"), "lw": 0.8})  # fmt: skip
    ticks = [0, 1, 10, 100, 1000]
    a.set_xticks(lg(ticks))
    a.set_xticklabels([str(t) for t in ticks])
    a.set_ylim(-1, 43)
    a.set_xlabel("limiting partner in the gene's second-best line (ARCHS4 CPM)")
    a.set_ylabel(f"shortlist lines (of 42) where every partner is >= {cut:.0f} CPM")
    a.legend(loc="upper center", bbox_to_anchor=(0.5, -0.1), ncol=2, fontsize=13, frameon=False)
    a.text(-0.1, 1.02, "a", transform=a.transAxes, fontsize=26, fontweight="bold")

    # b: which lines display each added gene (and the antigen-receptor subunits it leaves out).
    adds = c[status == "add"].sort_values(["archs4_lines_ge10", "gene"], ascending=[False, True])
    named = c[c["gene"].isin(["CD3D", "CD3E", "PTCRA", "TRB"]) & (status == "out")]
    rows = list(adds["gene"]) + list(named["gene"])
    hand = [x for x in IN_HAND if x in m10.columns and x != "K562"]
    rest = sorted(
        [x for x in m10.columns if x not in hand and x != "K562" and m10.loc[rows, x].any()]
    )
    cols = ["K562", *hand, *rest]
    mat = m10.loc[rows, cols]
    b.axvspan(-0.5, 0.5, color=get_color("primary"), alpha=0.12, lw=0)
    b.axvspan(0.5, len(hand) + 0.5, color=get_color("secondary"), alpha=0.07, lw=0)
    yy, xx = np.nonzero(mat.to_numpy())
    colr = [
        get_color("accent") if rows[r] in set(adds["gene"]) else get_color("neutral") for r in yy
    ]
    b.scatter(xx, yy, s=110, marker="s", c=colr, edgecolor="white", linewidth=0.5)
    b.set_xticks(range(len(cols)))
    b.set_xticklabels(cols, rotation=90, fontsize=12)
    b.set_yticks(range(len(rows)))
    b.set_yticklabels(rows, fontsize=12)
    b.set_ylim(len(rows) - 0.5, -0.5)
    b.set_xlim(-0.5, len(cols) - 0.5)
    b.axhline(len(adds) - 0.5, color=get_color("dark"), lw=1)
    b.text(0, -1.2, "K562", ha="center", fontsize=12, color=get_color("primary"))
    b.text(
        (len(hand) + 1) / 2,
        -1.2,
        "in house or ordered",
        ha="center",
        fontsize=12,
        color=get_color("secondary"),
    )
    b.text(-0.12, 1.02, "b", transform=b.transAxes, fontsize=26, fontweight="bold")
    b.grid(False)
    fig.tight_layout()
    save_figure(fig, "coreceptor_host_line", output_dir=HERE / "figures")


if __name__ == "__main__":
    raise SystemExit(main())
