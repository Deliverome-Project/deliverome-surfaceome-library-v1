"""Does every library proteoform actually have a structure to score tag sites from?

The tag-site lanes and EpicTope both need an AlphaFold model: pLDDT, RSA (freesasa) and
secondary structure (pydssp) are all computed off it. A candidate with no model cannot be
scored at all, by either. This measures that, and separately how many have a *deposited*
structure we could cross-check against -- see plot_loop_sites_vs_experiment.py for why
that cross-check matters.

Deposited coverage is counted **per gene**, which is an upper bound on what is usable: a
gene can have many structures and none covering the tag site. AXL is the worked example --
10 PDB chains, all of them the kinase domain (514-818), none reaching its N-terminal
cassette junction at 33.

Sampling is seeded, and results are cached in ``structure_coverage.tsv`` so the figure
redraws without re-querying. Pass --refresh to re-measure.

    uv run --locked python structure_coverage.py
    uv run --locked python structure_coverage.py --refresh
"""

from __future__ import annotations

import argparse
import csv
import json
import random
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

from deliverome_analysis.plotting import (
    COLORS,
    create_figure,
    get_color,
    save_figure,
    setup_plotting_style,
)

HERE = Path(__file__).parent
MANIFEST = HERE / "candidate_manifest_tranches.tsv"
CACHE = HERE / "structure_coverage.tsv"
SAMPLE_N = 35
SEED = 7


def _get(url: str, timeout: int = 25):
    return json.load(urllib.request.urlopen(url, timeout=timeout))


def has_alphafold(acc: str) -> bool:
    try:
        urllib.request.urlopen(
            f"https://alphafold.ebi.ac.uk/api/prediction/{acc}", timeout=20
        ).read()
        return True
    except urllib.error.HTTPError, urllib.error.URLError, TimeoutError:
        return False


def accession_for(key_part: str, *, symbol: bool) -> str | None:
    """UniProt accession for a gene symbol or an Ensembl transcript."""
    if symbol:
        q = f"gene_exact:{key_part} AND organism_id:9606 AND reviewed:true"
    else:
        q = f"xref:{key_part}"
    try:
        hits = _get(
            "https://rest.uniprot.org/uniprotkb/search?query="
            + urllib.parse.quote(q)
            + "&fields=accession&format=json&size=1"
        ).get("results", [])
    except Exception:
        return None
    return hits[0]["primaryAccession"] if hits else None


def measure() -> list[dict]:
    rows = list(csv.DictReader(MANIFEST.open(), delimiter="\t"))
    out = []
    for setname in (
        "library_reference",
        "expression_rescue",
        "isoform",
        "ortholog",
        "giant_truncation",
    ):
        members = [r for r in rows if r["set"] == setname]
        keys = [r["key"].split(":") for r in members]
        # how the candidate is identified decides whether a model can even be looked up
        if setname in ("library_reference", "expression_rescue", "giant_truncation"):
            ids = [(k[-1], True) for k in keys]
        elif setname == "isoform":
            ids = [(k[2], not k[2].startswith("ENST")) for k in keys]
        else:
            ids = [(k[3], False) for k in keys if len(k) > 3]
        rng = random.Random(SEED)
        sample = rng.sample(ids, min(SAMPLE_N, len(ids)))
        have = 0
        checked = 0
        for ident, is_symbol in sample:
            acc = ident
            if is_symbol or ident.startswith("ENST"):
                acc = accession_for(ident, symbol=is_symbol)
            if not acc:
                continue
            checked += 1
            have += has_alphafold(acc)
        out.append(
            {
                "set": setname,
                "candidates": len(members),
                "sampled": checked,
                "with_alphafold": have,
                "alphafold_frac": f"{have / checked:.3f}" if checked else "",
            }
        )
        print(f"  {setname:<20} n={len(members):>5}  AlphaFold {have}/{checked}")
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--refresh", action="store_true", help="re-measure instead of using the cache")
    args = ap.parse_args(argv)

    if args.refresh or not CACHE.exists():
        data = measure()
        with CACHE.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(data[0]), delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(data)
    data = list(csv.DictReader(CACHE.open(), delimiter="\t"))

    setup_plotting_style()
    fig, ax = create_figure(figsize=(9.8, 4.8))
    data = sorted(data, key=lambda r: float(r["alphafold_frac"] or 0))
    ys = range(len(data))
    for y, r in zip(ys, data, strict=True):
        frac = float(r["alphafold_frac"] or 0)
        low = frac < 0.9
        ax.barh(
            y,
            frac * 100,
            height=0.52,
            zorder=2,
            color=get_color("primary") if low else COLORS["secondary"],
        )
        ax.annotate(
            f"{frac:.0%}  ({r['with_alphafold']}/{r['sampled']} sampled)",
            (frac * 100 + 1.4, y),
            va="center",
            fontsize=8.6,
            color=COLORS["dark"],
        )
    ax.axvline(100, color=COLORS["line"], lw=1.0, zorder=1)
    ax.set_yticks(list(ys))
    ax.set_yticklabels([f"{r['set']}\nn={int(r['candidates']):,}" for r in data], fontsize=8.6)
    ax.set_xlim(0, 132)
    ax.set_xlabel("candidates with an AlphaFold model (%), seeded sample")
    for side in ("right", "top"):
        ax.spines[side].set_visible(False)
    fig.tight_layout()
    save_figure(fig, "s4e9_structure_coverage", output_dir=HERE / "figures")
    print("wrote figures/s4e9_structure_coverage.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
