"""Does the new MANE+APPRIS reference agree with Tedman's constructs better than ours did?

Tedman's 766 screen-validated GPCRs are the only gene set with an independent,
experimentally exercised choice of reference proteoform, so they are the test of whether the
new selection policy moved toward or away from what somebody actually built and screened.

Four references per gene: **ours** (the UniProt canonical we build from today), **Tedman's**
(the Ensembl transcript he flagged canonical), **MANE Select**, and **the new policy's pick**
(MANE and APPRIS where they agree; the evidence-adjudicated choice where they do not, read
from ``reference_selection.tsv``).

    uv run --locked python compare_to_tedman.py
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

from select_principal import (
    APPRIS_RELEASE,
    UNIPROT_HUMAN,
    load_uniprot_proteoforms,
    parse_appris,
    parse_ensembl_pep,
    parse_uniprot_fasta,
    select_principal,
)

# Which candidate the policy's decision resolves to. Mirrors select_reference.py.
TAKES_APPRIS = {"appris_only", "resolved_appris"}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--surfaceome-repo", default="~/Git/accessible-surfaceome")
    ap.add_argument("--cache-dir", default="~/.cache/deliverome/appris")
    ap.add_argument(
        "--tedman",
        default=(
            "~/Git/accessible-surfaceome/.claude/worktrees/tag-sites-468"
            "/data/tag_sites/tedman_gpcr_controls.tsv"
        ),
    )
    ap.add_argument("--selection", default="reference_selection.tsv")
    args = ap.parse_args()

    repo = Path(args.surfaceome_repo).expanduser()
    cache = Path(args.cache_dir).expanduser()
    seqdir = repo / "data/external/deeptmhmm_surfaceome_inputs"

    tedman = list(csv.DictReader(Path(args.tedman).expanduser().open(), delimiter="\t"))
    selection = {
        r["symbol"]: r for r in csv.DictReader(open(args.selection), delimiter="\t")
    }
    cohort = parse_uniprot_fasta(seqdir / "human_canonical_non_hla.fasta")
    cohort_by_symbol = {r["symbol"]: r for r in cohort.values() if r["symbol"]}
    uniprot = load_uniprot_proteoforms(cache / UNIPROT_HUMAN)
    pep = parse_ensembl_pep(cache / "Homo_sapiens.GRCh38.pep.all.fa.gz")
    appris = parse_appris(cache / f"appris_homo_sapiens_{APPRIS_RELEASE}.txt")

    appris_by_symbol: dict[str, list[dict]] = defaultdict(list)
    for rows in appris.values():
        for r in rows:
            appris_by_symbol[r["symbol"]].append(r)

    ensg_for_acc: dict[str, set[str]] = defaultdict(set)
    with (repo / "data/external/uniprot_ensembl_xrefs/ensg_to_uniprot.tsv").open() as fh:
        next(fh)
        for line in fh:
            g, a = line.split("\t")[:2]
            ensg_for_acc[a].add(g)

    # Tedman's reference choice per gene: the transcript he flagged canonical.
    his = {}
    for r in tedman:
        if r["is_canonical"] == "true":
            his.setdefault(r["gene_symbol"], r)

    rows = []
    for symbol, t in sorted(his.items()):
        tedman_seq = pep.get(t["ensembl_transcript_id"])
        stored = cohort_by_symbol.get(symbol)
        pf = uniprot.get(symbol)
        ours = (stored or (pf[0] if pf and not pf[0]["is_isoform"] else None))
        ours_seq = ours["seq"] if ours else None

        arows = [
            x for g in sorted(ensg_for_acc.get(t["uniprot_acc"], ())) for x in appris.get(g, [])
        ] or appris_by_symbol.get(symbol, [])
        mane_rows = [x for x in arows if x["mane"] and x["enst"] in pep]
        mane_seq = pep[mane_rows[0]["enst"]] if mane_rows else None
        sel = select_principal(arows, pep) if arows else None
        appris_seq = sel["seq"] if sel else None

        decision = selection.get(symbol, {}).get("decision", "")
        if symbol not in selection:
            new_seq = None  # not a library member; the policy was never applied to it
        elif decision in TAKES_APPRIS:
            new_seq = appris_seq
        else:
            new_seq = mane_seq or appris_seq

        rows.append(
            {
                "gene_symbol": symbol,
                "in_library": symbol in selection,
                "decision": decision,
                "tedman_len": len(tedman_seq) if tedman_seq else "",
                "ours_len": len(ours_seq) if ours_seq else "",
                "mane_len": len(mane_seq) if mane_seq else "",
                "new_len": len(new_seq) if new_seq else "",
                "tedman_vs_ours": _cmp(tedman_seq, ours_seq),
                "tedman_vs_mane": _cmp(tedman_seq, mane_seq),
                "tedman_vs_new": _cmp(tedman_seq, new_seq),
                "new_vs_ours": _cmp(new_seq, ours_seq),
                "new_minus_tedman": (
                    len(new_seq) - len(tedman_seq) if new_seq and tedman_seq else ""
                ),
            }
        )

    with open("tedman_vs_new_selection.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote tedman_vs_new_selection.tsv ({len(rows)} genes)")

    def rate(key: str, subset=None) -> dict:
        rs = [r for r in (subset or rows) if r[key] in ("same", "differs")]
        c = Counter(r[key] for r in rs)
        return {
            "compared": len(rs),
            "same": c["same"],
            "differs": c["differs"],
            "agreement": round(100 * c["same"] / len(rs), 1) if rs else None,
        }

    lib = [r for r in rows if r["in_library"]]
    summary = {
        "tedman_genes": len(rows),
        "of_which_library_members": len(lib),
        "all_tedman_genes": {
            k: rate(k) for k in ("tedman_vs_ours", "tedman_vs_mane", "tedman_vs_new")
        },
        "library_members_only": {
            k: rate(k, lib) for k in ("tedman_vs_ours", "tedman_vs_mane", "tedman_vs_new")
        },
        "new_vs_ours": rate("new_vs_ours", lib),
        "moved_toward_tedman": sum(
            1 for r in lib if r["tedman_vs_ours"] == "differs" and r["tedman_vs_new"] == "same"
        ),
        "moved_away_from_tedman": sum(
            1 for r in lib if r["tedman_vs_ours"] == "same" and r["tedman_vs_new"] == "differs"
        ),
    }
    Path("tedman_vs_new_selection.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


def _cmp(a: str | None, b: str | None) -> str:
    if not a or not b:
        return ""
    return "same" if a == b else "differs"


if __name__ == "__main__":
    main()
