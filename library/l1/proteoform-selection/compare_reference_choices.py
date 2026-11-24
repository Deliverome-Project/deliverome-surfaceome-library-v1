"""Three reference proteoforms for the same gene: ours, Tedman's, and APPRIS's.

Every tagged construct starts from a choice of which protein sequence the gene *is*. Three
choices are on the table and they are made on different grounds:

* **ours** - the UniProt canonical, which is what ``build_constructs.py`` splices into and
  what every topology, cleavage and disorder track was computed on.
* **Tedman's** - the Ensembl transcript each of his 946 screen-validated GPCR constructs was
  built from, read from ``tedman_gpcr_controls.tsv``. 177 of his constructs are deliberately
  non-canonical isoforms, so his reference choice is taken as the transcript he flagged
  canonical for that gene.
* **APPRIS** - the principal isoform selected in ``select_principal.py``.

Disagreement only matters in so far as it moves residue numbering, so the length difference
is reported directly rather than a similarity score: a tag placed at residue N of the wrong
proteoform lands somewhere else entirely.

    uv run --locked python compare_reference_choices.py \
        --surfaceome-repo ~/Git/accessible-surfaceome \
        --tedman ~/Git/accessible-surfaceome/.claude/worktrees/tag-sites-468/data/tag_sites/tedman_gpcr_controls.tsv
"""

from __future__ import annotations

import argparse
import csv
import json
import statistics
from collections import Counter, defaultdict
from pathlib import Path

from select_principal import (
    APPRIS_RELEASE,
    classify_difference,
    ENSEMBL_RELEASE,
    UNIPROT_HUMAN,
    load_uniprot_proteoforms,
    parse_appris,
    parse_ensembl_pep,
    parse_uniprot_fasta,
    select_principal,
)


def histogram(deltas: list[int]) -> dict[str, int]:
    """Signed length differences, bucketed. Sign is kept: a longer reference and a shorter
    one fail differently, and the N-terminal cases are the ones that move a tag site."""
    buckets = [
        ("<= -100", lambda d: d <= -100),
        ("-99..-20", lambda d: -99 <= d <= -20),
        ("-19..-5", lambda d: -19 <= d <= -5),
        ("-4..-1", lambda d: -4 <= d <= -1),
        ("0", lambda d: d == 0),
        ("1..4", lambda d: 1 <= d <= 4),
        ("5..19", lambda d: 5 <= d <= 19),
        ("20..99", lambda d: 20 <= d <= 99),
        (">= 100", lambda d: d >= 100),
    ]
    return {name: sum(1 for d in deltas if test(d)) for name, test in buckets}


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
    ap.add_argument("--out-dir", default=".")
    args = ap.parse_args()

    repo = Path(args.surfaceome_repo).expanduser()
    cache = Path(args.cache_dir).expanduser()
    out_dir = Path(args.out_dir)
    seqdir = repo / "data/external/deeptmhmm_surfaceome_inputs"

    tedman = list(csv.DictReader(Path(args.tedman).expanduser().open(), delimiter="\t"))
    print(f"Tedman constructs: {len(tedman)}")

    cohort = parse_uniprot_fasta(seqdir / "human_canonical_non_hla.fasta")
    cohort_by_symbol = {r["symbol"]: r for r in cohort.values() if r["symbol"]}
    uniprot = load_uniprot_proteoforms(cache / UNIPROT_HUMAN)
    by_accession = {
        r["accession"]: r for recs in uniprot.values() for r in recs
    }
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
            ensg, acc = line.split("\t")[:2]
            ensg_for_acc[acc].add(ensg)

    # Tedman's reference choice per gene: the transcript he flagged canonical. Genes where
    # he only ever built non-canonical isoforms have no reference choice to compare.
    his_reference: dict[str, dict] = {}
    isoform_constructs: dict[str, list[dict]] = defaultdict(list)
    for r in tedman:
        if r["is_canonical"] == "true":
            his_reference.setdefault(r["gene_symbol"], r)
        else:
            isoform_constructs[r["gene_symbol"]].append(r)
    print(
        f"genes with a Tedman canonical transcript: {len(his_reference)}; "
        f"genes with only non-canonical constructs: "
        f"{len(set(isoform_constructs) - set(his_reference))}"
    )

    rows = []
    for symbol, t in sorted(his_reference.items()):
        acc = t["uniprot_acc"]
        mine = cohort_by_symbol.get(symbol) or by_accession.get(acc)
        mine_source = (
            "deeptmhmm_cohort" if symbol in cohort_by_symbol else "uniprot_reviewed"
        )
        his_seq = pep.get(t["ensembl_transcript_id"])

        ensgs = sorted(ensg_for_acc.get(acc, ()))
        arows = [x for g in ensgs for x in appris.get(g, [])] or appris_by_symbol.get(
            symbol, []
        )
        sel = select_principal(arows, pep) if arows else None

        rec = {
            "gene_symbol": symbol,
            "uniprot_acc": acc,
            "mine_source": mine_source if mine else "",
            "mine_len": len(mine["seq"]) if mine else "",
            "tedman_enst": t["ensembl_transcript_id"],
            "tedman_len": len(his_seq) if his_seq else "",
            "tedman_in_ensembl_115": bool(his_seq),
            "appris_class": sel["label"] if sel else "",
            "appris_enst": sel["enst"] if sel else "",
            "appris_len": len(sel["seq"]) if sel and sel["seq"] else "",
            "n_tedman_isoform_constructs": len(isoform_constructs.get(symbol, [])),
        }
        for label, a, b in (
            ("tedman_vs_mine", his_seq, mine["seq"] if mine else None),
            ("appris_vs_mine", sel["seq"] if sel else None, mine["seq"] if mine else None),
            ("appris_vs_tedman", sel["seq"] if sel else None, his_seq),
        ):
            if not a or not b:
                continue
            rec[label] = "same" if a == b else "differs"
            rec[label.replace("_vs_", "_minus_")] = len(a) - len(b)
            if a != b:
                rec[label + "_kind"] = classify_difference(a, b)["difference_kind"]
        rows.append(rec)

    columns = [c for r in rows for c in r]
    columns = list(dict.fromkeys(columns))
    with (out_dir / "reference_choice_three_way.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote reference_choice_three_way.tsv  ({len(rows)} genes)")

    def agreement(key: str) -> dict:
        vals = [r[key] for r in rows if key in r]
        return {"compared": len(vals), **dict(Counter(vals))}

    def kinds(key: str) -> dict:
        k = [r[key + "_kind"] for r in rows if key + "_kind" in r]
        return dict(sorted(Counter(k).items(), key=lambda kv: -kv[1]))

    def spread(key: str) -> dict:
        d = [r[key] for r in rows if key in r and r[key] != 0]
        vs = key.replace("_minus_", "_vs_")
        same_length = sum(
            1 for r in rows if r.get(vs) == "differs" and r.get(key) == 0
        )
        if not d:
            return {"n_differing_in_length": 0, "differ_at_same_length": same_length}
        a = sorted(abs(x) for x in d)
        return {
            "n_differing_in_length": len(d),
            "differ_at_same_length": same_length,
            "shorter_than_reference": sum(1 for x in d if x < 0),
            "median_abs": statistics.median(a),
            "p90_abs": a[int(0.9 * len(a))],
            "max_abs": a[-1],
            "within_5_aa": sum(1 for x in a if x <= 5),
            "histogram": histogram(d),
            "difference_kind": kinds(vs),
        }

    summary = {
        "appris_release": APPRIS_RELEASE,
        "ensembl_release": ENSEMBL_RELEASE,
        "tedman_constructs": len(tedman),
        "genes_compared": len(rows),
        "tedman_transcripts_retired_from_ensembl_115": sum(
            1 for r in rows if not r["tedman_in_ensembl_115"]
        ),
        "agreement": {
            k: agreement(k)
            for k in ("tedman_vs_mine", "appris_vs_mine", "appris_vs_tedman")
        },
        "length_difference_where_they_differ": {
            k: spread(k)
            for k in ("tedman_minus_mine", "appris_minus_mine", "appris_minus_tedman")
        },
    }
    (out_dir / "reference_choice_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
