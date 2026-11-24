#!/usr/bin/env python3
"""Does the model find the EndoNB knock-in sites? Leave-the-gene-out.

Lenaerts et al. 2025 put an ALFA tag into five endogenous loci and showed the
tagged protein still reaches the surface. That is the only set of internal
insertion sites in the literature that was (a) chosen by somebody else, (b)
knocked into the genome rather than overexpressed from a plasmid, and (c)
reported as working. It is the closest thing we have to a prospective test.

`endonb_sodn.py` puts the five sites in our numbering by translating the donor
oligos. This scores them. Three of the five genes carry a published site in our
training corpus, so every number here comes from a model refitted with that
GENE held out entirely -- otherwise the model would be recognising its own
training data and the answer would mean nothing.

ITGB1 is the strongest of the three: its corpus site IS the EndoNB site, so
holding the gene out removes the answer. TFRC is the most interesting: its
corpus site is residue 111, nowhere near 290, so under +/-10 masking residue
290 was an explicit NEGATIVE in training. The model is being asked to rank
highly a residue it was taught was not a site.

    uv run --locked python endonb_scored.py
"""

from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[3] / "src"))

from independent_model_search import make, smooth, within_protein  # noqa: E402
from tag_site_model_v4 import FEATS, KIND, SPEC, training_frame  # noqa: E402

TOL = 2  # the scores are smoothed over +/-2, so credit a call within 2 residues


def leave_gene_out(train, allx, acc, feats):
    """Rank and absolute scores for every extracellular residue of `acc`.

    Keyed on ACCESSION, not symbol: the cohort table falls back to the
    accession wherever the symbol was missing upstream, so AXL's rows are
    filed under "P30530" and a symbol filter silently returns nothing.
    """
    sub = allx[allx["uniprot_acc"] == acc].copy().reset_index(drop=True)
    held = set(allx.loc[allx["uniprot_acc"] == acc, "gene_symbol"]) | {acc}
    tr = train[~train["gene_symbol"].isin(held) & (train["uniprot_acc"] != acc)]
    tr = tr.reset_index(drop=True)
    y = tr["label"].to_numpy()
    mz, mraw = make(KIND, feats), make(KIND, feats)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        mz.fit(within_protein(tr, feats, "z"), y)
        mraw.fit(tr, y)
        sub["rank_score"] = smooth(sub, mz.predict_proba(within_protein(sub, feats, "z"))[:, 1], 2)
        sub["absolute_score"] = smooth(sub, mraw.predict_proba(sub)[:, 1], 2)
    return sub.sort_values("rank_score", ascending=False).reset_index(drop=True)


def rank_of(sub, residue, col):
    """Best rank achieved within TOL residues of the target, on `col`."""
    s = sub.sort_values(col, ascending=False).reset_index(drop=True)
    near = s.index[(s["residue"] - residue).abs() <= TOL]
    if len(near) == 0:
        return None, None
    i = int(near.min())
    return i + 1, float(s.loc[i, col])


def main() -> int:
    sites = json.loads((HERE / "endonb_sites.json").read_text())
    spec = json.loads(SPEC.read_text())
    ceilings = np.array(spec["published_protein_ceilings_oof"])
    med = float(np.quantile(ceilings, 0.5))

    pos, allx, train = training_frame()
    # ITGB5 is in neither the library nor the published corpus, so no features
    # were ever exported for it. `export_site_features.py --extra-accessions
    # P18084` fills that in; the rows are kept beside this script rather than
    # folded into the cohort table, which other results already depend on.
    extra = HERE / "itgb5_extracellular.tsv"
    if extra.is_file():
        e = pd.read_csv(extra, sep="\t", low_memory=False)
        for c in ("ss3", "topology_class", "region_class"):
            e[c] = e[c].fillna("").astype(str)
        e["gene_symbol"] = e["gene_symbol"].fillna(e["uniprot_acc"]).astype(str)
        allx = pd.concat([allx, e], ignore_index=True)
    corpus = {a: sorted(v) for a, v in pos.groupby("uniprot_acc")["residue"].apply(list).items()}

    # one call per gene: the per-isoform rows all name the same site
    seen: set[str] = set()
    wanted = []
    for s in sites:
        if s["proteoform"] in set(allx["uniprot_acc"]) and s["gene"] not in seen:
            seen.add(s["gene"])
            wanted.append(s)
    rows = []
    for s in wanted:
        gene, acc, res = s["gene"], s["proteoform"], s["insert_after_residue"]
        sub = leave_gene_out(train, allx, acc, FEATS)
        held_out = acc not in set(train.loc[train["label"] == 1, "uniprot_acc"])
        n = len(sub)
        r_rank, r_score = rank_of(sub, res, "rank_score")
        _a_rank, a_score = rank_of(sub, res, "absolute_score")
        top = sub.iloc[0]
        ceiling = float(sub["absolute_score"].max())
        row = {
            "gene": gene,
            "uniprot_acc": acc,
            "paper_label": s["paper_label"],
            "endonb_residue_ours": res,
            "extracellular_residues": n,
            "rank_of_endonb_site": r_rank,
            "percentile": round(100 * (1 - (r_rank - 1) / n), 1) if r_rank else None,
            "rank_score_there": round(r_score, 4) if r_score else None,
            "absolute_score_there": round(a_score, 4) if a_score else None,
            "model_top_residue": int(top["residue"]),
            "model_top_rank_score": round(float(top["rank_score"]), 4),
            "protein_ceiling": round(ceiling, 4),
            "ceiling_clears_median_published_protein": bool(ceiling >= med),
            "same_region_as_model_top": bool(
                sub.loc[(sub["residue"] - res).abs() <= TOL, "region_index"]
                .isin([top["region_index"]])
                .any()
            ),
            "corpus_sites_for_this_gene": corpus.get(acc, []),
            "never_had_a_positive_anywhere": held_out,
            "endonb_site_was_a_training_negative": bool(
                acc in corpus and all(abs(p - res) > 10 for p in corpus[acc])
            ),
        }
        for other in ("prodomino", "epictope"):
            if other in sub.columns and sub[other].notna().any():
                orank, _ = rank_of(sub, res, other)
                row[f"{other}_rank"] = orank
                row[f"{other}_percentile"] = (
                    round(100 * (1 - (orank - 1) / n), 1) if orank else None
                )
            else:
                row[f"{other}_rank"] = None
                row[f"{other}_percentile"] = None
        rows.append(row)
        print(
            f"{gene:8s} {s['paper_label']:5s} -> ours after {res:<4d}  "
            f"rank {r_rank:>4}/{n:<4} ({row['percentile']:>5.1f}th pct)  "
            f"model's own top {row['model_top_residue']:<4d}  "
            f"ProDomino {row['prodomino_rank']}  EpicTope {row['epictope_rank']}",
            flush=True,
        )

    out = pd.DataFrame(rows)
    out.to_csv(HERE / "endonb_scored.tsv", sep="\t", index=False)
    skipped = sorted({s["gene"] for s in sites} - {r["gene"] for r in rows})
    print(f"\nwrote endonb_scored.tsv ({len(out)} sites)")
    if skipped:
        print(
            f"not scored, no feature rows in this cohort: {', '.join(skipped)} "
            "-- neither is in the 515 nor in the published corpus, so no features were "
            "ever exported for them"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
