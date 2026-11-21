#!/usr/bin/env python3
"""Does the L1 refold check tell tag insertions that work from ones that fail?

The refold check (``tag_refold/score_tagged_folds.py``) folds each tagged L1 construct with
ColabFold and compares it with a model of the untagged protein. Its thresholds were set on 17
untagged controls only, so nothing has yet shown that a construct it passes works better
than one it fails. This script puts the check on constructs whose outcome is known.

The L1 cassette (51 aa, ``GGGGS``-SpyTag003-``GGGGSGGGGS``-``P``-ALFA-``P``-``GGGGS``, the
amino acids every ``l1_cassette_encodings.tsv`` row translates to) is inserted at:

    positive          a confirmed literature internal tag site (corpus audit gold standard)
                      where the paper shows surface display or function
    known_failure     a literature site where the tagged construct failed: the 5 rows of
                      ``failed_sites.tsv`` plus the gold-standard sites shown NOT to reach
                      the surface
    function_lost     gold-standard sites shown to reach the surface (or not tested) but to
                      lose function; reported, not pooled with the failures
    decoy_core        one buried, confidently placed helix/strand residue of the same
                      protein's extracellular region (relative SASA < 0.05 in its AFDB model,
                      pLDDT >= 80), more than 10 residues from any literature-cited site
    decoy_hardneg     one hard negative of the shipped site ranker's training frame
                      (buried core, disulfide Cys, N-glycan, junction or ligand contact) in
                      the same protein, more than 10 residues from any cited site

Every tagged protein and every untagged native is folded with the L1 refold settings
(``scripts/fold_tagged_modal.py``: ColabFold 1.6.2, alphafold2_ptm model 1, 3 recycles, no
relax, no templates, native MSA with the cassette as gap columns) and scored with
``score_tagged_folds.score()`` unchanged, as an internal-loop construct, twice:

    afdb_ref   against the protein's archived AFDB model, as production scores most members
    own_ref    against our own untagged refold with the same settings

Implemented by Claude Code (subagent of the l1-brief-review session, 2026-11-20) for Becca's
question whether the refold check discriminates; classes and decoy rules are this session's.

    uv run --locked python analysis/sprint4/s4e9-orf-design-fg-panel/refold_check_validation.py build
    uv run --locked --group modal modal run --detach scripts/fold_tagged_modal.py \\
        --queue analysis/sprint4/s4e9-orf-design-fg-panel/refold_check_validation_queue.tsv \\
        --run refold-validation-2026-11-20 --tiers 1 --spawn
    uv run --locked --group modal modal run scripts/fold_tagged_modal.py \\
        --run refold-validation-2026-11-20 --collect-only \\
        --out-dir ~/.cache/deliverome-analysis/tagged-folds/refold-validation-2026-11-20
    uv run --locked python analysis/sprint4/s4e9-orf-design-fg-panel/refold_check_validation.py score
    uv run --locked python analysis/sprint4/s4e9-orf-design-fg-panel/refold_check_validation.py report
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(HERE / "tag_refold"))

import score_tagged_folds as stf  # noqa: E402
import structure_metrics as sm  # noqa: E402
from build_refold_inputs import CASSETTE  # noqa: E402

AUDIT = REPO / "scratch/becca/tag-site-corpus-audit"
OUT = HERE / "refold-check-validation"
FIG = HERE / "figures"
FOLDS = Path.home() / ".cache/deliverome-analysis/tagged-folds/refold-validation-2026-11-20"
RUN = "refold-validation-2026-11-20"
S3_OUT = "sprint4/s4e9-orf-design-fg-panel/refold-check-validation/2026-11-20/"
#: The shipped site ranker's pinned feature and training tables (tag_site_model_benchmark.py).
S3_INPUTS = "sprint5/s5e13-tag-site-benchmarking/v5/2026-11-14/"
INPUT_SHA256 = {
    "v5_training_table.tsv": "e7d948688cf47f1cbf1ad3a6ff39a4e256f1822a6091c9b3174122568e10d23c",
    "all_extracellular_targets.tsv": "d3dcdac545b1bfe322df2965f5510c45c68eee780fb4b8c08628ea21462ffefe",
}
#: KCNH7 has no archived AFDB model; its canonical sequence is pinned here by SHA-256
#: (UniProt Q9NS40 v2, fetched 2026-11-20) and it is scored against its own refold only.
UNIPROT_ONLY = {"Q9NS40": "https://rest.uniprot.org/uniprotkb/Q9NS40.fasta"}
UNIPROT_SHA256 = {"Q9NS40": "4b4bf5412fee1e9e82191d444e0788e1f5ca809cee5a9a72392e055a2ea0ece1"}
SEED = 20261120


def T(name: str) -> Path:
    """A committed table at the experiment root (the only place .gitignore lets a TSV live)."""
    return HERE / f"refold_check_validation_{name}"


FLANK = 10  # a decoy lies more than this many residues from any literature-cited site
CORE_RSA = 0.05
CORE_PLDDT = 80.0
#: Small-protein stratum: 22 of the 32 distinct genes among L1's 37 held CDSs (EMP1-3, PMP22,
#: MS4A, TM4SF, CD52, CD70, ...) are 300 residues or shorter, where a 51-residue cassette is a
#: sixth or more of the chain.
SMALL_AA = 300
CLASSES = ["positive", "known_failure", "function_lost", "decoy_core", "decoy_hardneg"]
TIER = {"known_failure": 1, "positive": 2, "function_lost": 3, "decoy_core": 3, "decoy_hardneg": 4}


def sha(s: str) -> str:
    return hashlib.sha256(s.encode()).hexdigest()


def pinned_inputs() -> tuple[pd.DataFrame, pd.DataFrame]:
    data = OUT / "data"
    data.mkdir(parents=True, exist_ok=True)
    for name, want in INPUT_SHA256.items():
        p = data / name
        if not p.is_file():
            from deliverome_analysis.s3 import processed

            processed.download(S3_INPUTS + name, p)
        got = hashlib.sha256(p.read_bytes()).hexdigest()
        assert got == want, f"{name}: SHA-256 {got} != pinned {want}"
    allx = pd.read_csv(data / "all_extracellular_targets.tsv", sep="\t", low_memory=False)
    tr = pd.read_csv(data / "v5_training_table.tsv", sep="\t")
    return allx, tr


def afdb_rows() -> dict[str, dict[str, str]]:
    with (HERE / "tag_refold_afdb_archive.tsv").open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return {r["afdb_url"].split("/AF-")[1].split("-F1")[0]: r for r in rows}


def native_of(acc: str, arch: dict) -> tuple[str, str, sm.Model | None]:
    """(sequence, reference label, AFDB model): the human canonical protein the sites number."""
    if acc in arch:
        m = sm.Model(stf.archived_afdb(arch[acc]["afdb_url"]))
        return m.seq, arch[acc]["afdb_url"], m
    import urllib.request

    with urllib.request.urlopen(UNIPROT_ONLY[acc], timeout=60) as r:
        seq = "".join(r.read().decode().splitlines()[1:])
    assert sha(seq) == UNIPROT_SHA256[acc], f"{acc}: UniProt sequence changed since 2026-11-20"
    return seq, "", None


def tagged_row(cls: str, acc: str, gene: str, site: int, native: str, note: str, order: int):
    tagged = native[:site] + CASSETTE + native[site:]
    nsha, tsha = sha(native), sha(tagged)
    return {
        "fold_id": f"v.{cls}.{acc}.{site}",
        "kind": "tagged",
        "tier": TIER[cls],
        "order": order,
        "rung": "loop",
        "design_id": f"{acc}:{site}",
        "source_id": note,
        "gene_symbol": gene,
        "class": cls,
        "uniprot_acc": acc,
        "native_aa": len(native),
        "tagged_aa": len(tagged),
        "cassette_start": site,
        "cassette_end": site + len(CASSETTE),
        "insert_after_native_residue": site,
        "leader_end": "",
        "sp_end": "",
        "tm1_span_folded": "",
        "native_sha256": nsha,
        "tagged_sha256": tsha,
        "blocks": json.dumps([[0, 0, site], [site + len(CASSETTE), site, len(native) - site]]),
        "native_seq": native,
        "tagged_seq": tagged,
    }


def build() -> None:
    gold = pd.read_csv(AUDIT / "gold_standard.tsv", sep="\t")
    failed = pd.read_csv(AUDIT / "failed_sites.tsv", sep="\t")
    audit = pd.read_csv(AUDIT / "audit.tsv", sep="\t")
    allx, tr = pinned_inputs()
    arch = afdb_rows()

    sites: list[tuple[str, str, str, int, str]] = []  # class, acc, gene, site, note
    for g in gold.itertuples():
        if g.surface_shown == "no":
            cls = "known_failure"  # incl. PTH1R 93: function shown, surface shown absent
        elif g.surface_shown == "yes" or g.function_shown == "yes":
            cls = "positive"
        elif g.function_shown == "no":
            cls = "function_lost"
        else:
            continue  # neither display nor function reported
        note = (
            f"gold_standard PMID {g.pmid}; surface {g.surface_shown}; function {g.function_shown}"
        )
        sites.append((cls, g.uniprot_acc, g.gene, int(g.insert_after_residue), note))
    for f in failed.itertuples():
        sites.append(
            ("known_failure", f.uniprot_acc, f.gene, int(f.insert_after_residue),
             f"failed_sites PMID {f.pmid}: {f.why_failed}")
        )  # fmt: skip

    # residues any literature row cites, for the decoy exclusion zone
    cited: dict[str, set[int]] = {}
    for acc, r in zip(audit["uniprot_acc"], audit["insert_after_residue"], strict=True):
        if pd.notna(r):
            cited.setdefault(acc, set()).add(int(r))
    for _, acc, _, s, _ in sites:
        cited.setdefault(acc, set()).add(s)

    natives: dict[str, tuple[str, str, sm.Model | None]] = {}
    for acc in sorted({s[1] for s in sites}):
        natives[acc] = native_of(acc, arch)

    rng = np.random.default_rng(SEED)
    genes = {acc: gene for _, acc, gene, _, _ in sites}
    for acc in sorted({s[1] for s in sites if s[0] == "positive"}):
        seq, _, model = natives[acc]
        far = np.array(sorted(cited[acc]))
        ecd = allx[allx["uniprot_acc"] == acc]
        if model is None or ecd.empty:
            continue
        ecd_res = set(ecd["residue"].astype(int))
        rsa = model.rsa()
        ss = ecd.set_index("residue")["ss3"].to_dict()
        cand = [
            r
            for r in range(FLANK + 1, len(seq) - FLANK)
            if r in ecd_res
            and np.abs(far - r).min() > FLANK
            and rsa[r - 1] < CORE_RSA
            and model.plddt[r - 1] >= CORE_PLDDT
            and ss.get(r) in ("H", "E")
        ]
        if cand:
            # the most buried 5-residue window among the candidates
            w = {r: float(np.mean(rsa[r - 3 : r + 2])) for r in cand}
            r = min(cand, key=lambda k: (w[k], k))
            sites.append(("decoy_core", acc, genes[acc], r,
                          f"buried core: AFDB RSA {rsa[r - 1]:.3f}, pLDDT {model.plddt[r - 1]:.0f}, ss {ss[r]}"))  # fmt: skip
        hn = tr[(tr["uniprot_acc"] == acc) & (tr["weight"] == 5) & (tr["label"] == 0)]
        hn = [int(r) for r in hn["residue"] if np.abs(far - int(r)).min() > FLANK]
        if hn:
            r = int(rng.choice(sorted(hn)))
            sites.append(("decoy_hardneg", acc, genes[acc], r,
                          "site-ranker hard negative (weight 5)"))  # fmt: skip

    rows = []
    seen = set()
    for order, (cls, acc, gene, site, note) in enumerate(sites):
        seq = natives[acc][0]
        assert 0 < site < len(seq), (acc, site, len(seq))
        row = tagged_row(cls, acc, gene, site, seq, note, order)
        if row["fold_id"] in seen:
            continue
        seen.add(row["fold_id"])
        rows.append(row)
    # one untagged native per protein, at the tier of its first tagged form
    for acc, (seq, ref, _) in natives.items():
        tier = min(r["tier"] for r in rows if r["uniprot_acc"] == acc)
        rows.append({
            "fold_id": f"n.{sha(seq)[:16]}", "kind": "native", "tier": tier, "order": -1,
            "rung": "native", "design_id": acc, "source_id": ref or UNIPROT_ONLY[acc],
            "gene_symbol": genes[acc], "class": "native", "uniprot_acc": acc,
            "native_aa": len(seq), "tagged_aa": len(seq), "cassette_start": "",
            "cassette_end": "", "insert_after_native_residue": "", "leader_end": "", "sp_end": "",
            "tm1_span_folded": "", "native_sha256": sha(seq), "tagged_sha256": sha(seq),
            "blocks": json.dumps([[0, 0, len(seq)]]), "native_seq": seq, "tagged_seq": seq,
        })  # fmt: skip
    OUT.mkdir(parents=True, exist_ok=True)
    q = pd.DataFrame(rows)
    q.to_csv(T("queue.tsv"), sep="\t", index=False)
    print(q.groupby(["tier", "class"]).agg(n=("fold_id", "size"), residues=("tagged_aa", "sum")))
    print(f"total {len(q)} folds, {q['tagged_aa'].sum():,} residues")


def _index_afdb(q: pd.DataFrame) -> dict:
    arch = afdb_rows()
    idx = {}
    for acc, nsha in q[["uniprot_acc", "native_sha256"]].drop_duplicates().itertuples(index=False):
        if acc in arch:
            url = arch[acc]["afdb_url"]
            idx[nsha] = {"native_sha256": nsha, "ref": url, "local": str(stf.archived_afdb(url))}
    return idx


def score() -> None:
    q = pd.read_csv(T("queue.tsv"), sep="\t", keep_default_na=False, dtype=str)
    rows = q.to_dict("records")
    idx = _index_afdb(q)
    from concurrent.futures import ProcessPoolExecutor

    for label, index in (("afdb_ref", idx), ("own_ref", {})):
        tasks = []
        for r in rows:
            if r["kind"] == "native":
                if label == "own_ref":
                    continue
                # the untagged refold against AFDB: what the settings alone cost, per protein
                r = {**r, "rung": "control", "kind": "control"}
            tasks.append((r, FOLDS, index))
        with ProcessPoolExecutor(8) as ex:
            scored = list(ex.map(stf._score_one, tasks, chunksize=2))
        df = pd.DataFrame(scored)
        meta = q.set_index("fold_id")[["class", "uniprot_acc", "source_id"]]
        df = df.drop(columns=[c for c in ("source_id",) if c in df]).join(meta, on="fold_id")
        df.insert(0, "reference", label)
        df.to_csv(T(f"metrics_{label}.tsv"), sep="\t", index=False)
        print(label, df.groupby(["class", "verdict"]).size().to_string(), sep="\n")


# ------------------------------------------------------------------ report


def auroc(pos: np.ndarray, neg: np.ndarray) -> float:
    """P(score of a positive > score of a negative), ties half: the Mann-Whitney AUROC."""
    pos, neg = pos[np.isfinite(pos)], neg[np.isfinite(neg)]
    if not len(pos) or not len(neg):
        return float("nan")
    d = pos[:, None] - neg[None, :]
    return float(((d > 0).sum() + 0.5 * (d == 0).sum()) / d.size)


def boot_auroc(df: pd.DataFrame, col: str, neg: str, n: int = 2000) -> tuple[float, float]:
    """95% CI resampling proteins (sites on one protein are not independent)."""
    rng = np.random.default_rng(SEED)
    accs = df["uniprot_acc"].unique()
    by = {a: g for a, g in df.groupby("uniprot_acc")}
    vals = []
    for _ in range(n):
        s = pd.concat([by[a] for a in rng.choice(accs, len(accs))])
        vals.append(
            auroc(
                s.loc[s["class"] == "positive", col].to_numpy(float),
                s.loc[s["class"] == neg, col].to_numpy(float),
            )
        )
    v = np.array(vals)
    v = v[np.isfinite(v)]
    return float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))


def wilson(k: int, n: int) -> tuple[float, float]:
    if n == 0:
        return float("nan"), float("nan")
    z, p = 1.96, k / n
    c = (p + z * z / (2 * n)) / (1 + z * z / n)
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return c - h, c + h


#: Continuous scores, oriented so that higher = "looks fine" to the check.
SCORES = {
    "lddt_confident": "native-region lDDT",
    "lddt_site": "site lDDT (+/-15)",
    "epitope_rsa_min": "epitope exposure (min RSA)",
    "neg_epitope_contacts": "epitope contacts (negated)",
}


def load(label: str) -> pd.DataFrame:
    df = pd.read_csv(T(f"metrics_{label}.tsv"), sep="\t")
    df = df[df["class"].isin(CLASSES) & df["verdict"].isin(["pass", "fail"])].copy()
    df["epitope_rsa_min"] = df[["spytag003_rsa", "alfa_rsa"]].min(axis=1)
    df["neg_epitope_contacts"] = -df["epitope_confident_contact_res"]
    return df


def report() -> None:
    out_rates, out_auc, out_pairs = [], [], []
    for label in ("afdb_ref", "own_ref"):
        df = load(label)
        for size, keep in (
            ("all", df["native_aa"] > 0),
            (f"<= {SMALL_AA} aa", df["native_aa"] <= SMALL_AA),
            (f"> {SMALL_AA} aa", df["native_aa"] > SMALL_AA),
        ):
            for cls in CLASSES:
                s = df[keep & (df["class"] == cls)]
                k = int((s["verdict"] == "pass").sum())
                lo, hi = wilson(k, len(s))
                out_rates.append({"reference": label, "size_bin": size, "class": cls,
                                  "n_sites": len(s), "n_proteins": s["uniprot_acc"].nunique(),
                                  "pass": k, "pass_rate": round(k / len(s), 3) if len(s) else "",
                                  "wilson95_lo": round(lo, 3), "wilson95_hi": round(hi, 3)})  # fmt: skip
        for neg in ("known_failure", "decoy_core", "decoy_hardneg", "function_lost"):
            for col, name in [*SCORES.items(), ("pass01", "verdict (pass=1)")]:
                d = df.assign(pass01=(df["verdict"] == "pass").astype(float))
                sub = d[d["class"].isin(["positive", neg])]
                a = auroc(
                    sub.loc[sub["class"] == "positive", col].to_numpy(float),
                    sub.loc[sub["class"] == neg, col].to_numpy(float),
                )
                lo, hi = boot_auroc(sub, col, neg) if np.isfinite(a) else (np.nan, np.nan)
                out_auc.append({"reference": label, "negative_class": neg, "score": name,
                                "n_positive": int((sub["class"] == "positive").sum()),
                                "n_negative": int((sub["class"] == neg).sum()),
                                "auroc": round(a, 3), "boot95_lo": round(lo, 3),
                                "boot95_hi": round(hi, 3)})  # fmt: skip
        # within-protein: does the confirmed site score better than the same protein's decoy?
        for neg in ("decoy_core", "decoy_hardneg", "known_failure"):
            for col, name in SCORES.items():
                wins = ties = n = 0
                for _acc, g in df.groupby("uniprot_acc"):
                    p = g.loc[g["class"] == "positive", col].dropna()
                    q = g.loc[g["class"] == neg, col].dropna()
                    if len(p) and len(q):
                        n += 1
                        dv = p.mean() - q.mean()
                        wins += dv > 0
                        ties += dv == 0
                out_pairs.append({"reference": label, "negative_class": neg, "score": name,
                                  "n_proteins": n, "positive_higher": wins, "tied": ties})  # fmt: skip
    rates = pd.DataFrame(out_rates)
    aucs = pd.DataFrame(out_auc)
    pairs = pd.DataFrame(out_pairs)
    rates.to_csv(T("pass_rates.tsv"), sep="\t", index=False)
    aucs.to_csv(T("auroc.tsv"), sep="\t", index=False)
    pairs.to_csv(T("within_protein.tsv"), sep="\t", index=False)
    pd.set_option("display.width", 200)
    print(rates.to_string(index=False))
    print(aucs[aucs["negative_class"] != "function_lost"].to_string(index=False))
    print(pairs.to_string(index=False))
    figure()


def figure() -> None:
    import matplotlib.pyplot as plt

    from deliverome_analysis.plotting import (
        CATEGORICAL_PALETTE,
        COLORS,
        save_figure,
        setup_plotting_style,
    )

    df = load("afdb_ref")
    rates = pd.read_csv(T("pass_rates.tsv"), sep="\t")
    rates = rates[rates["reference"] == "afdb_ref"]
    setup_plotting_style(font_scale=1.0)
    shown = [c for c in ("positive", "known_failure", "decoy_core", "decoy_hardneg")
             if (df["class"] == c).any()]  # fmt: skip
    labels = {
        "positive": "literature\nsite works",
        "known_failure": "literature\nsite failed",
        "decoy_core": "decoy:\nburied core",
        "decoy_hardneg": "decoy: ranker\nhard negative",
    }
    color = {c: CATEGORICAL_PALETTE[i] for i, c in enumerate(shown)}
    panels = [
        ("lddt_confident", r"native-region C$\alpha$ lDDT", stf.LDDT_PASS),
        ("lddt_site", r"site C$\alpha$ lDDT ($\pm$15 residues)", stf.LOCAL_LDDT_PASS),
        ("epitope_rsa_min", "epitope RSA (lower of the two)", stf.EPITOPE_RSA_PASS),
    ]
    fig, axes = plt.subplots(1, 4, figsize=(18, 5.0), gridspec_kw={"width_ratios": [1, 1, 1, 0.9]})
    rng = np.random.default_rng(0)
    for ax, (col, ylab, thr) in zip(axes[:3], panels, strict=True):
        for x, c in enumerate(shown):
            v = df.loc[df["class"] == c, col].dropna().to_numpy(float)
            ax.scatter(x + rng.uniform(-0.25, 0.25, len(v)), v, s=12, alpha=0.7,
                       color=color[c], linewidths=0)  # fmt: skip
            if len(v):
                ax.hlines(np.median(v), x - 0.33, x + 0.33, color=COLORS["dark"], lw=1.5)
            ax.text(x, 1.01, f"n={len(v)}", ha="center", va="bottom", fontsize=8,
                    transform=ax.get_xaxis_transform())  # fmt: skip
        ax.axhline(thr, color=COLORS["primary"], ls="--", lw=1)
        ax.set_xticks(range(len(shown)), [labels[c] for c in shown], fontsize=8)
        ax.set_ylabel(ylab)
    ax = axes[3]
    bins = [("all", 0.0, "all proteins", 1.0), (f"<= {SMALL_AA} aa", 0.36, f"proteins <= {SMALL_AA} aa", 0.45)]  # fmt: skip
    for size, dx, _lab, alpha in bins:
        r_all = rates[rates["size_bin"] == size].set_index("class")
        for x, c in enumerate(shown):
            r = r_all.loc[c]
            if not r["n_sites"]:
                continue
            ax.bar(x - 0.18 + dx, r["pass_rate"], color=color[c], width=0.34, alpha=alpha)
            ax.errorbar(x - 0.18 + dx, r["pass_rate"],
                        yerr=[[r["pass_rate"] - r["wilson95_lo"]], [r["wilson95_hi"] - r["pass_rate"]]],
                        color=COLORS["dark"], capsize=2, lw=1)  # fmt: skip
            ax.text(x - 0.18 + dx, 0.03, f"{int(r['pass'])}/{int(r['n_sites'])}", ha="center", rotation=90,
                    va="bottom", fontsize=7, color=COLORS["dark"])  # fmt: skip
    from matplotlib.patches import Patch

    ax.legend(
        handles=[Patch(color=COLORS["neutral"], alpha=a, label=lab) for _, _, lab, a in bins],
        fontsize=7, frameon=False, loc="upper center", bbox_to_anchor=(0.5, 1.13), ncol=2,
    )  # fmt: skip
    ax.set_ylim(0, 1)
    ax.set_xticks(range(len(shown)), [labels[c] for c in shown], fontsize=8)
    ax.set_ylabel("fraction passing the refold check")
    fig.tight_layout()
    save_figure(fig, "refold_check_validation", output_dir=FIG)


def upload() -> None:
    """Publish each output under its own key (never a folder-wide upload)."""
    from deliverome_analysis.s3 import processed

    names = ["queue.tsv", "metrics_afdb_ref.tsv", "metrics_own_ref.tsv", "pass_rates.tsv",
             "auroc.tsv", "within_protein.tsv"]  # fmt: skip
    for n in names:
        processed.upload(T(n), S3_OUT + n)
        print("uploaded", S3_OUT + n)
    png = FIG / "refold_check_validation.png"
    processed.upload(png, S3_OUT + "figures/" + png.name)
    print("uploaded", S3_OUT + "figures/" + png.name)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("step", choices=["build", "score", "report", "figure", "upload"])
    a = ap.parse_args()
    {"build": build, "score": score, "report": report, "figure": figure, "upload": upload}[a.step]()


if __name__ == "__main__":
    main()
