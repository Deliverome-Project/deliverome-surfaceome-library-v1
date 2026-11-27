"""Split the 7,066 library candidates into an L1 synthesis tranche and an L2 remainder.

L1 is what gets synthesised first: every human reference, the 25 expression-rescue arms,
and the two prioritisation panels (150 isoforms, 200 orthologs). L2 is the rest -- the
isoforms and orthologs outside the panels, plus the giant truncations.

Why the manifest is re-cut rather than annotated in place: every row currently carries an
``l1-pNNNN`` id, which asserts first-tranche membership for all 7,066. Renumbering is safe
*now* and only now -- ``tag_cassette`` is ``deferred`` and ``barcode`` is ``not reserved``
on every row, so no id is bound to a physical reagent yet. Once barcodes are reserved the
id becomes the stable key a sequencing readout joins on and must never move.

Determinism: ids are assigned by (set precedence, key) with no randomness, so a re-run on
the same inputs produces byte-identical output. ``prior_candidate_id`` is carried so the
previous numbering stays recoverable.

    uv run --locked python split_library_tranches.py            # dry run, prints the split
    uv run --locked python split_library_tranches.py --write    # write locally
    uv run --locked python split_library_tranches.py --write --publish
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from deliverome_analysis.core.s3 import processed

PREFIX = "sprint4/s4e9-orf-design-fg-panel/"
MANIFEST = PREFIX + "library-candidates/candidate_manifest.tsv"
ISOFORM_PANEL = PREFIX + "library-design/isoform_panel.tsv"
ORTHOLOG_PANEL = PREFIX + "library-design/ortholog_panel.tsv"

#: L1 before L2, and within a tranche this order, so ids group by what they are.
SET_ORDER = {
    "library_reference": 0,
    "expression_rescue": 1,
    "isoform": 2,
    "ortholog": 3,
    "giant_truncation": 4,
}
OUT_COLS = ["candidate_id", "library", "panel", "prior_candidate_id"]


def read_tsv(key: str) -> list[dict[str, str]]:
    return list(csv.DictReader(open(processed.download(key)), delimiter="\t"))


def panel_keys() -> tuple[set[str], set[str], list[tuple[str, str]]]:
    """Manifest keys named by each panel, plus the panel rows that match nothing.

    The ortholog panel reconstructs its key exactly. The isoform panel names a gene and a
    semicolon-separated transcript list, so a row is matched on gene AND transcript; a
    gene-only fallback would silently pull in a sibling isoform the panel did not draw,
    which is the quiet wrong-construct failure this join exists to avoid.
    """
    manifest = read_tsv(MANIFEST)
    iso_keys = {r["key"] for r in manifest if r["set"] == "isoform"}
    ort_keys = {r["key"] for r in manifest if r["set"] == "ortholog"}

    iso_hit: set[str] = set()
    unmatched: list[tuple[str, str]] = []
    for row in read_tsv(ISOFORM_PANEL):
        gene, accession = row["gene"], row["uniprot"]
        # The manifest keys an isoform by WHICHEVER identifier resolved it -- a UniProt
        # accession for some (iso:ACHE:P22303-2), an Ensembl transcript for others
        # (iso:CACNA2D1:ENST00000957014). The panel carries both, so accept either.
        # Matching on transcript alone recovers 62 of 150; the accession carries the rest.
        wanted = {f"iso:{gene}:{accession}"} if accession else set()
        wanted |= {
            f"iso:{gene}:{t.strip()}"
            for t in (row.get("ensembl_transcripts") or "").split(";")
            if t.strip()
        }
        found = wanted & iso_keys
        if found:
            iso_hit |= found
        else:
            unmatched.append((f"isoform:{gene}", accession or "<no accession>"))

    ort_hit: set[str] = set()
    for row in read_tsv(ORTHOLOG_PANEL):
        key = f"ort:{row['species']}:{row['human_symbol']}:{row['picked_acc']}"
        if key in ort_keys:
            ort_hit.add(key)
            continue
        # Do NOT fall back to (species, gene): the manifest may hold several proteoforms
        # of one ortholog gene and taking the wrong one substitutes a construct the panel
        # did not draw. Report instead.
        siblings = sorted(
            k for k in ort_keys if k.startswith(f"ort:{row['species']}:{row['human_symbol']}:")
        )
        unmatched.append(
            (
                f"ortholog:{row['species']}:{row['human_symbol']}",
                f"panel {row['picked_acc']} vs manifest {siblings or '<none>'}",
            )
        )
    return iso_hit, ort_hit, unmatched


def split(
    manifest: list[dict[str, str]], iso_hit: set[str], ort_hit: set[str]
) -> list[dict[str, str]]:
    rows = []
    for r in manifest:
        s, key = r["set"], r["key"]
        panel = ""
        if s in ("library_reference", "expression_rescue"):
            tranche = "l1"
        elif key in iso_hit:
            tranche, panel = "l1", "isoform_panel"
        elif key in ort_hit:
            tranche, panel = "l1", "ortholog_panel"
        else:
            tranche = "l2"
        rows.append({**r, "library": tranche, "panel": panel})
    out = []
    for tranche in ("l1", "l2"):
        group = sorted(
            (r for r in rows if r["library"] == tranche),
            key=lambda r: (SET_ORDER.get(r["set"], 99), r["key"]),
        )
        for n, r in enumerate(group, 1):
            out.append(
                {
                    **r,
                    "prior_candidate_id": r["candidate_id"],
                    "candidate_id": f"{tranche}-p{n:04d}",
                }
            )
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--write", action="store_true")
    ap.add_argument("--publish", action="store_true")
    args = ap.parse_args(argv)

    manifest = read_tsv(MANIFEST)
    iso_hit, ort_hit, unmatched = panel_keys()
    rows = split(manifest, iso_hit, ort_hit)

    counts: dict[str, dict[str, int]] = {}
    for r in rows:
        counts.setdefault(r["library"], {}).setdefault(r["set"], 0)
        counts[r["library"]][r["set"]] += 1
    summary = {
        "source_manifest": f"processed://{MANIFEST}",
        "isoform_panel_matched": len(iso_hit),
        "ortholog_panel_matched": len(ort_hit),
        "unmatched_panel_rows": [{"panel_row": a, "accession": b} for a, b in unmatched],
        "counts": counts,
        "totals": {k: sum(v.values()) for k, v in counts.items()},
    }
    print(json.dumps(summary, indent=1))
    assert len(rows) == len(manifest), "row count changed"
    assert len({r["candidate_id"] for r in rows}) == len(rows), "duplicate id"

    if not args.write:
        print("\ndry run -- re-run with --write")
        return 0
    cols = list(manifest[0]) + [c for c in OUT_COLS if c not in manifest[0]]
    out = Path("candidate_manifest_tranches.tsv")
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    Path("library_tranche_summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(f"\nwrote {out} and library_tranche_summary.json")
    if args.publish:
        for p, key in (
            (out, "candidate_manifest_tranches.tsv"),
            (Path("library_tranche_summary.json"), "library_tranche_summary.json"),
        ):
            processed.upload(
                str(p),
                PREFIX + "library-candidates/" + key,
                content_type="text/tab-separated-values"
                if p.suffix == ".tsv"
                else "application/json",
            )
            print("published", PREFIX + "library-candidates/" + key)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
