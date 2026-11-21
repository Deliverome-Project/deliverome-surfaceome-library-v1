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
# ARCHS4 v2.5 per-line consensus built by the cell-line panel experiment (log2(CPM + 1), gene
# symbols x panel lines): covers all 42 shortlist lines, including the 13 HPA lacks. CPM is not
# TPM, so the CPM cut matching each TPM cut is calibrated on the 29 lines both sources share.
ARCHS4_KEY = "scratch/becca/cell-line-selection/rnaseq_consensus_archs4.parquet"
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
    from deliverome_analysis.s3 import processed

    path = CACHE / f"archs4_consensus_{ARCHS4_SHA256[:12]}.parquet"
    if not path.exists() or sha256(path) != ARCHS4_SHA256:
        processed.download(ARCHS4_KEY, path)
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
    cpm_cut = {5: float(best[GATE_TPM]), 10: float(best[10.0])}
    hand, hand_missing = shortlist(list(tpm.columns), IN_HAND)
    assert set(hand) <= set(short) and len(hand) == 10, (hand, hand_missing)

    prow, grow = [], []
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
                continue
            a13 = lines_ok_sym(cpm, sym_groups, archs4_only, cpm_cut[thr]) if sym_groups else None
            h29 = (
                lines_ok(tpm, reqs_cur, [h for h, _ in pairs], float(thr))
                if reqs_cur
                else pd.Series(False, index=[h for h, _ in pairs])
            )
            out[f"combined_lines_ge{thr}"] = int(h29.sum()) + (0 if a13 is None else int(a13.sum()))
            if thr == 5:
                out["archs4_only_lines_ge5_names"] = (
                    "" if a13 is None else ", ".join(a13[a13].index)
                )
            a42 = lines_ok_sym(cpm, sym_groups, list(hpa_of), cpm_cut[thr]) if sym_groups else None
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
        "HPA 29 + ARCHS4 13, >= 10 TPM (or CPM equivalent) in >= 2 of 42 lines OR in K562": (
            genes["combined_lines_ge10"] >= SHORTLIST_MIN_LINES
        )
        | rule("curated_k562", 10),
        "ARCHS4 alone, CPM equivalent of >= 5 TPM in >= 2 of 42 lines OR in K562": (
            genes["archs4_lines_ge5"] >= SHORTLIST_MIN_LINES
        )
        | (genes["archs4_k562_ok"] == True),  # noqa: E712
        "ARCHS4 alone, CPM equivalent of >= 10 TPM in >= 2 of 42 lines": genes["archs4_lines_ge10"]
        >= SHORTLIST_MIN_LINES,
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
    cand = genes["passes_other_gates"].to_numpy()
    now = genes["coreceptor_ok_now"].to_numpy()
    names = genes["gene"].to_numpy()
    orows = []
    for label, ok in options.items():
        ok = ok.to_numpy()
        genes[f"ok: {label}"] = ok
        orows.append(
            {
                "rule": label,
                "genes_kept": int((ok & cand).sum()),
                "added": int((ok & ~now & cand).sum()),
                "dropped": int((~ok & now & cand).sum()),
                "added_genes": ", ".join(names[ok & ~now & cand]),
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
    calib.to_csv(HERE / "coreceptor_host_line_archs4_calibration.tsv", sep="\t", index=False,
                 float_format="%.3f")  # fmt: skip
    meta = {
        "archs4": f"processed://{ARCHS4_KEY}",
        "archs4_sha256": ARCHS4_SHA256,
        "archs4_cpm_cut_for_tpm": cpm_cut,
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
    plot(genes, partners, lin, len(short))
    return 0


def plot(genes: pd.DataFrame, partners: pd.DataFrame, lin: dict, n_short: int) -> None:
    def lg(v):
        return np.log10(np.asarray(v, dtype=float) + 0.1)

    fig, (a, b) = create_figure(1, 2, figsize=(22, 12), gridspec_kw={"width_ratios": [1.15, 1]})
    c = genes[genes["passes_other_gates"] & ~genes["self_only"]].dropna(subset=["curated_median"])
    med_ok = c["curated_median"] >= GATE_TPM
    sh_ok = c["curated_shortlist_lines_ge5"] >= SHORTLIST_MIN_LINES
    groups = [
        (med_ok & sh_ok, "passes both", get_color("secondary")),
        (med_ok & ~sh_ok, "median only (a shortlist rule drops)", get_color("primary")),
        (~med_ok & sh_ok, "shortlist only (a shortlist rule adds)", get_color("accent")),
        (~med_ok & ~sh_ok, "passes neither", get_color("neutral")),
    ]
    # Integer counts overlap; a fixed small vertical jitter separates the points.
    jit = np.random.default_rng(0).uniform(-0.3, 0.3, len(c))
    y = c["curated_shortlist_lines_ge5"].to_numpy() + jit
    x = lg(c["curated_median"])
    for m, lab, col in groups:
        m = m.to_numpy()
        a.scatter(
            x[m],
            y[m],
            s=70,
            color=col,
            alpha=0.85,
            edgecolor="white",
            linewidth=0.6,
            label=f"{lab} (n = {int(m.sum())})",
        )
    k = (c["curated_k562"] >= GATE_TPM).to_numpy()
    a.scatter(
        x[k],
        y[k],
        s=170,
        facecolor="none",
        edgecolor=get_color("dark"),
        linewidth=1.2,
        label=f"limiting partner >= 5 TPM in K562 (n = {int(k.sum())})",
    )
    t = lg(GATE_TPM)
    a.axvline(t, color=get_color("dark"), lw=1.2)
    a.axhline(SHORTLIST_MIN_LINES - 0.5, color=get_color("dark"), lw=1.2, ls="--")
    a.text(t - 0.03, 24, "current gate:\nmedian 5 TPM", fontsize=13, va="center", ha="right")
    a.text(
        lg(200), SHORTLIST_MIN_LINES - 0.4, f">= {SHORTLIST_MIN_LINES} shortlist lines", fontsize=13
    )
    # Genes below the median cut that clear at least one shortlist line, listed by line count
    # (the points overlap too much to label individually).
    low = c[~med_ok & (c["curated_shortlist_lines_ge5"] >= 1)]
    rows_txt = []
    for n, grp in low.groupby("curated_shortlist_lines_ge5"):
        names = sorted(grp["gene"])
        txt = ", ".join(names) if len(names) <= 9 else ", ".join(names[:9]) + f" +{len(names) - 9}"
        rows_txt.append(f"{int(n)}: {txt}")
    a.text(
        lg(6),
        17.5,
        "below the median cut, by shortlist lines\n" + "\n".join(rows_txt[::-1]),
        fontsize=11.5,
        va="top",
        linespacing=1.5,
    )
    ticks = [0, 1, 5, 10, 100, 1000]
    a.set_xticks(lg(ticks))
    a.set_xticklabels([str(v) for v in ticks])
    a.set_ylim(-0.8, n_short + 0.8)
    a.set_xlabel("limiting partner, median TPM over 1,206 lines")
    a.set_ylabel(f"shortlist lines (of {n_short}) where every partner is >= 5 TPM")
    a.legend(loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=2, fontsize=12, frameon=False)
    a.text(-0.08, 1.02, "a", transform=a.transAxes, fontsize=24, fontweight="bold")

    # b: the named immune complexes; one row per gene, its limiting partner by K562.
    rows = genes.set_index("gene").reindex(NAMED).dropna(subset=["all_k562"])
    cols = [
        ("all_k562", "K562", get_color("primary"), "o"),
        ("all_median", "median, 1,206 lines", get_color("dark"), "|"),
        ("all_best_T", f"best T line (of {len(lin['T'])})", get_color("secondary"), "s"),
        ("all_best_B", f"best B line (of {len(lin['B'])})", get_color("quaternary"), "D"),
        (
            "all_best_myeloid",
            f"best myeloid line (of {len(lin['myeloid'])})",
            get_color("accent"),
            "^",
        ),
    ]
    y = np.arange(len(rows))[::-1]
    for col, lab, colr, mk in cols:
        b.scatter(
            lg(rows[col]),
            y,
            color=colr,
            marker=mk,
            s=140 if mk != "|" else 400,
            label=lab,
            zorder=3,
            linewidths=2.5 if mk == "|" else 0.6,
            edgecolor="white" if mk != "|" else None,
        )
    b.axvline(t, color=get_color("dark"), lw=1.2)
    b.set_yticks(y)
    b.set_yticklabels(
        [
            f"{g}  ({rows.loc[g, 'k562_limiting_partner']}; "
            f"{int(rows.loc[g, 'curated_shortlist_lines_ge5'])}/{n_short})"
            for g in rows.index
        ],
        fontsize=14,
    )
    b.set_xticks(lg(ticks))
    b.set_xticklabels([str(x) for x in ticks])
    b.set_xlabel("limiting partner TPM")
    b.legend(loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=3, fontsize=12, frameon=False)
    b.text(-0.32, 1.02, "b", transform=b.transAxes, fontsize=24, fontweight="bold")
    fig.tight_layout()
    save_figure(fig, "coreceptor_host_line", output_dir=HERE / "figures")


if __name__ == "__main__":
    raise SystemExit(main())
