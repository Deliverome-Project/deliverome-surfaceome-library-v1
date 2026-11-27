"""Which L1 CDSs need a host accessory factor to reach the surface, and does K562 have it?

The library is screened first in K562. Some receptors reach the plasma membrane of a
heterologous cell only when an accessory protein that their native cell type expresses is
co-expressed: odorant receptors need RTP1/RTP2 (Saito et al., 2004), and some bitter taste
receptors respond better with RTP3/RTP4 (Behrens et al., 2006). The S1E3 co-receptor gate
does not catch these: the catalog lists an accessory factor for only two of the ORs in L1, and
for the rest it lists nothing, so the gate has nothing to test. A CDS whose factor K562 lacks
will read as "fails to display" for a reason unrelated to the CDS.

This lists every active L1 CDS (``l1_tag_paths.tsv``) in a family with a documented
accessory-factor or obligate-partner requirement, with:

- the factor(s), the evidence and the PMIDs (each resolved and read at the abstract level;
  see ``EVIDENCE``), and whether the evidence is for this receptor or for its family;
- the factor's TPM in K562 and in HEK293 (HPA v23; HPA has no HEK293T line);
- the entry route of the CDS's gene (core, low-literature UniProt tier, screen-validated GPCR
  tier, overexpression tier, dated addition), from the reviewer-package gate lists;
- the gene's surface immunostaining in the Tedman et al. (2026) HEK293T GPCR screen, which
  co-expressed no accessory factor (no RTP, REEP or RIC8B), and whether the screen's ORF is
  the same Ensembl transcript as the L1 CDS.

Genes are joined on Ensembl gene ID from ``gene_resolver`` (MyGene; the response is committed
in ``accessory_factor_dependence_mygene/``), never on symbol. Orthologs and isoforms inherit
their human gene's route and factor.

Nothing in the library, constructs or order changes.

    uv run --locked python accessory_factor_dependence.py [--upload]
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import re
import zipfile
from pathlib import Path

import pandas as pd

from deliverome_analysis.core.s3 import external, processed
from deliverome_analysis.surfaceome import gene_resolver

HERE = Path(__file__).parent
CACHE = HERE / "data"  # gitignored
MYGENE = HERE / "accessory_factor_dependence_mygene"  # committed resolver answer
PREFIX = "sprint4/s4e9-orf-design-fg-panel/"
OUT_PREFIX = PREFIX + "accessory-factors/"

HPA_KEY = "hpa/v23/rna_celline.tsv.zip"  # HPA 23.0 (2023-06-19, Ensembl 109)
HPA_SHA256 = "e5b11779308c2720cfe89d9cb19932ad4c0c3fea0fa37e8c61a14d501995852a"
TEDMAN_KEY = "tedman-gpcr-surface-screen/2025-09-19/tedman2025_gpcr_screen_media-2.xlsx"
TEDMAN_SHA256 = "c0fbcbd1df5ab488e026239fb55aa92a7a8d9171154bec20bc837f6fb9685286"
GATE_LISTS = {  # reviewer-package gate lists: one table per entry tier
    "core": "reviewer-package/gate_lists/10_tier_core.tsv",
    "low-literature UniProt tier": "reviewer-package/gate_lists/11_tier_low_literature_uniprot.tsv",
    "screen-validated GPCR tier": "reviewer-package/gate_lists/12_tier_screen_validated_gpcr.tsv",
    "overexpression-surface tier": "reviewer-package/gate_lists/13_tier_overexpression_surface.tsv",
}
SCREEN_LINE, HEK_LINE = "K-562", "HEK293"
PRESENT_TPM = 5.0  # same bar as the co-receptor gate
SCREEN_POSITIVE = 1000.0  # anti-HA a.u.; the screen-validated tier's threshold (willow_rescue.py)

# Evidence, each PMID resolved in PubMed on 2026-11-20 and its claim checked against the
# abstract (and, for Tedman, the full text). "scope" says whether the evidence is for the
# receptor itself or for its family.
EVIDENCE = {
    "OR": dict(
        factors=[["RTP1", "RTP2"], ["REEP1"]],
        need="surface delivery",
        claim=(
            "RTP1 and RTP2 promote cell-surface expression of ORs in HEK293T, REEP1 more weakly "
            "(Saito 2004); the short RTP1 isoform RTP1S gives robust surface expression "
            "(Zhuang & Matsunami 2007); standard OR surface-expression protocol co-expresses "
            "RTPs (Zhuang & Matsunami 2008); a subset of ORs reaches the surface without RTPs "
            "(Sharma 2017)"
        ),
        pmids="15550249;17387175;18772867;28262096",
        scope="family",
    ),
    "TAS2R": dict(
        factors=[["RTP3", "RTP4"], ["REEP2"]],
        need="surface delivery (some receptors)",
        claim=(
            "TAS2Rs show marginal surface expression in HEK293T; RTP/REEP co-expression changed "
            "responses of some hTAS2Rs, with better surface localization by biotinylation, and "
            "TAS2R16 co-immunoprecipitates with RTP3 and RTP4 (Behrens 2006); REEP2 enhances "
            "bitter and sweet receptor responses (Ilegems 2010)"
        ),
        pmids="16720576;20943918",
        scope="family; receptor-specific effect",
    ),
    "TAS1R": dict(
        factors=[["REEP2"]],
        need="signalling, not surface delivery",
        claim=(
            "REEP2 associates with TAS1R2 and TAS1R3 and enhances sweet-receptor responses, "
            "but does not increase their cell-surface expression (Ilegems 2010)"
        ),
        pmids="20943918",
        scope="family",
    ),
    "CHRNA7": dict(
        factors=[["TMEM35A"]],
        need="surface delivery",
        claim=(
            "alpha7 does not assemble into functional pentamers in cell lines; NACHO (TMEM35A) "
            "promotes alpha7 folding, Golgi maturation and cell-surface expression (Gu 2016)"
        ),
        pmids="26875622",
        scope="receptor",
    ),
    "GABRA": dict(
        factors=[["GABRB1", "GABRB2", "GABRB3"]],
        need="surface delivery (obligate partner)",
        claim=(
            "single GABA-A subunits are retained in the ER; alpha1 reaches the surface only "
            "with beta2 (Connolly 1996, alpha1/beta2/gamma2L in A293 cells). Applied here to "
            "the other alpha subunits by family, not tested per subunit"
        ),
        pmids="8550630",
        scope="alpha1 tested; family extrapolation",
    ),
    "ATP1A": dict(
        factors=[["ATP1B1", "ATP1B2", "ATP1B3"]],
        need="surface delivery (obligate partner)",
        claim=(
            "the beta subunit is required for structural and functional maturation of the "
            "Na,K-ATPase alpha subunit (Geering 2008, review)"
        ),
        pmids="18695395",
        scope="family (review)",
    ),
}
# Families the screen measured or the brief asked about, with no accessory factor found in
# the literature checked here. Listed so their screen values sit beside the others.
NO_FACTOR = {"TAAR", "VN1R"}
# Heteromer subunits whose partner the catalog does not list, so the co-receptor gate never
# tested it. Literature not checked here: listed as candidates, not counted.
CANDIDATES = {
    "CHRNA4": "CHRNB2;CHRNB4",
    "CHRNA6": "CHRNB2;CHRNB3;CHRNB4",
    "CHRNB3": "CHRNA4;CHRNA6",
    "GABRE": "GABRA1;GABRB1;GABRB2;GABRB3",
    "GABRP": "GABRA1;GABRB1;GABRB2;GABRB3",
    "GRIK4": "GRIK1;GRIK2;GRIK3",
    "SCNN1A": "SCNN1B;SCNN1G",
    "SCNN1D": "SCNN1B;SCNN1G",
}
# The two L1 ORs for which the catalog does list an accessory factor (co-receptor host-line
# partner table, branch claude/l1-coreceptor-host-line-e60c8b): the gate dropped RTP1S as
# unresolvable for OR4S2 and failed OR5B21 on RTP1/RTP2, yet both entered by a named tier.
CATALOG_LISTED = {"OR4S2": "RTP1S (unresolved)", "OR5B21": "RTP1; RTP2"}


def family(symbol: str) -> str | None:
    if re.match(r"^OR\d", symbol):
        return "OR"
    for fam in ("TAS2R", "TAS1R", "TAAR", "VN1R", "GABRA", "ATP1A"):
        if symbol.startswith(fam):
            return fam
    if symbol == "CHRNA7":
        return "CHRNA7"
    if symbol in CANDIDATES:
        return "candidate"
    return None


def sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def fetch(bucket, key: str, digest: str | None) -> Path:
    CACHE.mkdir(exist_ok=True)
    path = CACHE / Path(key).name
    if not path.exists() or (digest and sha256(path) != digest):
        bucket.download(key, path)
    if digest:
        assert sha256(path) == digest, f"{key}: SHA-256 does not match the pinned file"
    return path


def hpa(lines: list[str]) -> pd.DataFrame:
    z = fetch(external, HPA_KEY, HPA_SHA256)
    with zipfile.ZipFile(z) as zf, zf.open(zf.namelist()[0]) as fh:
        long = pd.read_csv(io.TextIOWrapper(fh), sep="\t", usecols=["Gene", "Cell line", "TPM"])
    long = long[long["Cell line"].isin(lines)]
    return long.pivot(index="Gene", columns="Cell line", values="TPM")


def _upload_outputs() -> None:
    """Publish the three outputs to processed://."""
    for name in (
        "accessory_factor_dependence.tsv",
        "accessory_factor_dependence_factors.tsv",
        "accessory_factor_dependence_summary.json",
    ):
        print("uploaded", processed.upload(HERE / name, OUT_PREFIX + name))


def main(upload: bool) -> int:
    tp = pd.read_csv(HERE / "l1_tag_paths.tsv", sep="\t", dtype=str)
    tp["family"] = tp.gene_symbol.map(family)
    rows = tp[tp.family.notna()].copy()

    gate = {}
    for route, key in GATE_LISTS.items():
        path = fetch(processed, PREFIX + key, None)
        gate[route] = (pd.read_csv(path, sep="\t", dtype=str), sha256(path))
    additions = pd.read_csv(HERE / "gene_list_additions.tsv", sep="\t", dtype=str)
    if "withdrawn" in additions:  # a withdrawn addition is no longer on the gene list
        additions = additions[additions["withdrawn"].fillna("").str.strip() == ""]
    tx = fetch(external, TEDMAN_KEY, TEDMAN_SHA256)
    screen = pd.read_excel(tx, sheet_name="Canonical")
    screen = screen[screen.gene_id.notna()]
    screen["enst"] = screen["Receptor Name"].str.extract(r"(ENST\d+)")

    factors = sorted(
        {f for ev in EVIDENCE.values() for grp in ev["factors"] for f in grp}
        | {f for v in CANDIDATES.values() for f in v.split(";")}
    )
    symbols = (
        set(tp.gene_symbol)
        | set(screen.gene_id)
        | set(factors)
        | set(additions.gene)
        | {g for df, _ in gate.values() for g in df.gene}
    )
    res = gene_resolver.resolve(symbols, cache_dir=MYGENE)

    tpm = hpa([SCREEN_LINE, HEK_LINE])

    def ens(sym: str) -> str | None:
        """One Ensembl gene per symbol. MyGene also returns alternate-haplotype loci (MHC
        region ORs, the TAS2R cluster); keep the IDs HPA (primary assembly) measures."""
        ids = sorted(res.get(sym, {}).get("ensembl") or [])
        if len(ids) > 1:
            ids = [i for i in ids if i in tpm.index]
        return ids[0] if len(ids) == 1 else None

    unresolved = sorted({s for s in rows.gene_symbol if ens(s) is None})
    assert not unresolved, f"L1 symbols without one Ensembl gene: {unresolved}"

    route_of: dict[str, str] = {}
    gate_unresolved: list[str] = []
    for route, (df, _) in gate.items():
        for g in df.gene:
            e = ens(g)
            if e is None:
                gate_unresolved.append(g)  # reported; every listed CDS must still get a route
                continue
            route_of.setdefault(e, route)
    for g in additions.gene:
        route_of.setdefault(ens(g), "dated addition")

    def tpm_of(sym: str, line: str) -> float | None:
        e = ens(sym)
        return None if e is None or e not in tpm.index else round(float(tpm.at[e, line]), 1)

    # Screen: one value per gene (max over the gene's canonical rows, as willow_rescue.py).
    screen["ensembl"] = screen.gene_id.map(ens)
    by_gene = screen.dropna(subset=["ensembl"]).sort_values(
        "Immunostaining Intensity", ascending=False
    )
    by_gene = by_gene.drop_duplicates("ensembl").set_index("ensembl")

    manifest = pd.read_csv(HERE / "candidate_manifest_tranches.tsv", sep="\t", dtype=str)
    transcript = manifest.set_index("candidate_id").transcript.to_dict()

    out = []
    for r in rows.itertuples():
        e = ens(r.gene_symbol)
        fam = r.family
        ev = EVIDENCE.get(fam)
        if ev:
            groups = ev["factors"]
            # one factor from each group suffices; the group's best K562 value counts
            k562 = [max((tpm_of(f, SCREEN_LINE) or 0.0) for f in grp) for grp in groups]
            hek = [max((tpm_of(f, HEK_LINE) or 0.0) for f in grp) for grp in groups]
            primary_k562, primary_hek = k562[0], hek[0]
            category = (
                "accessory factor, not for surface"
                if ev["need"].startswith("signalling")
                else "accessory factor or obligate partner for surface"
            )
            factor_txt = " + optional ".join("/".join(g) for g in groups)
        elif fam == "candidate":
            groups = [CANDIDATES[r.gene_symbol].split(";")]
            primary_k562 = max((tpm_of(f, SCREEN_LINE) or 0.0) for f in groups[0])
            primary_hek = max((tpm_of(f, HEK_LINE) or 0.0) for f in groups[0])
            category = "candidate obligate partner, literature not checked"
            factor_txt = "/".join(groups[0])
        else:
            primary_k562 = primary_hek = None
            category = "no accessory factor found"
            factor_txt = ""
        s = by_gene.loc[e] if e in by_gene.index else None
        l1_tx = transcript.get(r.candidate_id) or ""
        out.append(
            dict(
                candidate_id=r.candidate_id,
                candidate_set=r.candidate_set,
                species=r.species,
                gene_symbol=r.gene_symbol,
                ensembl_gene=e,
                family=fam,
                category=category,
                entry_route=route_of.get(e, "not on a gate list"),
                catalog_listed_partner=CATALOG_LISTED.get(r.gene_symbol, ""),
                factor=factor_txt,
                need=ev["need"] if ev else "",
                evidence=ev["claim"] if ev else "",
                pmids=ev["pmids"] if ev else "",
                evidence_scope=ev["scope"] if ev else "",
                factor_k562_tpm=primary_k562,
                factor_hek293_tpm=primary_hek,
                factor_present_k562=(None if primary_k562 is None else primary_k562 >= PRESENT_TPM),
                tedman_measured=s is not None,
                tedman_class=s["GPCR Class"] if s is not None else "",
                tedman_construct=s["Receptor Name"] if s is not None else "",
                tedman_same_transcript=(
                    (s["enst"] == l1_tx) if s is not None and r.species == "human" else None
                ),
                tedman_immunostaining_au=(
                    round(float(s["Immunostaining Intensity"]), 1) if s is not None else None
                ),
                tedman_surface_positive=(
                    float(s["Immunostaining Intensity"]) >= SCREEN_POSITIVE
                    if s is not None
                    else None
                ),
            )
        )
    table = pd.DataFrame(out)
    assert (table.entry_route != "not on a gate list").all(), table[
        table.entry_route == "not on a gate list"
    ][["candidate_id", "gene_symbol"]]
    table = table.sort_values(["category", "family", "gene_symbol", "candidate_id"])
    table.to_csv(HERE / "accessory_factor_dependence.tsv", sep="\t", index=False)

    ftab = pd.DataFrame(
        [
            dict(
                factor=f,
                ensembl_gene=ens(f),
                k562_tpm=tpm_of(f, SCREEN_LINE),
                hek293_tpm=tpm_of(f, HEK_LINE),
            )
            for f in factors
        ]
    )
    ftab.to_csv(HERE / "accessory_factor_dependence_factors.tsv", sep="\t", index=False)

    surf = table[table.category == "accessory factor or obligate partner for surface"]
    lacking = surf[~surf.factor_present_k562.astype(bool)]
    summary = {
        "inputs": {
            "l1_tag_paths.tsv": sha256(HERE / "l1_tag_paths.tsv"),
            "hpa": f"external://{HPA_KEY} {HPA_SHA256}",
            "tedman": f"external://{TEDMAN_KEY} {TEDMAN_SHA256}",
            **{f"processed://{PREFIX}{GATE_LISTS[k]}": d for k, (_, d) in gate.items()},
        },
        "l1_active_cdss": len(tp),
        "rows": len(table),
        "gate_list_genes_without_one_ensembl_id": sorted(gate_unresolved),
        "by_category_family": {
            f"{c} | {f}": int(n)
            for (c, f), n in table.groupby(["category", "family"]).size().items()
        },
        "surface_dependent_lacking_factor_in_k562": {
            "cdss": len(lacking),
            "by_family": lacking.family.value_counts().to_dict(),
            "by_route": lacking.entry_route.value_counts().to_dict(),
            "by_set": lacking.candidate_set.value_counts().to_dict(),
            "tedman_measured": int(lacking.tedman_measured.sum()),
            "tedman_surface_positive": int(lacking.tedman_surface_positive.fillna(False).sum()),
            "tedman_same_transcript": int(lacking.tedman_same_transcript.fillna(False).sum()),
            "by_family_route_screen": {
                f"{fam} | {route} | {call}": int(n)
                for (fam, route, call), n in lacking.assign(
                    call=lacking.tedman_surface_positive.map(
                        {
                            True: f">= {SCREEN_POSITIVE:.0f} a.u.",
                            False: f"< {SCREEN_POSITIVE:.0f} a.u.",
                        }
                    ).fillna("not in screen")
                )
                .groupby(["family", "entry_route", "call"])
                .size()
                .items()
            },
        },
        "tedman_class_context": {
            cls: {
                "receptors": len(v),
                "median_au": round(float(v.median()), 0),
                "fraction_ge_threshold": round(float((v >= SCREEN_POSITIVE).mean()), 2),
            }
            for cls, v in screen.groupby("GPCR Class")["Immunostaining Intensity"]
        },
        "or_and_taste_by_route": {
            f"{fam} | {route}": int(n)
            for (fam, route), n in table[table.family.isin(["OR", "TAS1R", "TAS2R"])]
            .groupby(["family", "entry_route"])
            .size()
            .items()
        },
    }
    (HERE / "accessory_factor_dependence_summary.json").write_text(
        json.dumps(summary, indent=1, default=str) + "\n"
    )
    print(json.dumps(summary, indent=1, default=str))
    print(ftab.to_string(index=False))

    if upload:
        _upload_outputs()
    return 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--upload", action="store_true", help="publish outputs to processed://")
    raise SystemExit(main(ap.parse_args().upload))
