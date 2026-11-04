"""The authoritative GPI-anchored gene set, from NetGPI 1.1 plus curation.

A GPI anchor matters for tag placement because the omega site **deletes the C-terminus**:
an ecto-C tag on a GPI-anchored protein is removed with the propeptide, so the member needs
an internal site instead. Getting the set wrong therefore mis-plans real constructs.

Two sources, and neither alone is right:

* **NetGPI 1.1** (`netgpi_gpi_human_predictions.txt`, 918 forms over 143 reference and 770
  alternative forms). Measured on 148 curated-GPI references + 100 single-pass references:
  **sensitivity 89.9%, specificity 99.0%, MCC 0.875** against the paper's 0.895.
* **Curated annotation**, reached here through the ``gpi_anchored`` triage reason in
  ``triage_run_public`` (130 genes).

They disagree in both directions, and the disagreements are informative rather than noise:

* NetGPI finds **27 genes curation missed**, because a gene can be anchored in a form that
  is not the canonical reference — NCAM1 is anchored as NCAM-120 and transmembrane as
  NCAM-140, ACHE as the H subunit and not the canonical T subunit, likewise BCAN, TFPI and
  CD58. Scanning alternative forms is what surfaces these.
* Curation holds **genes NetGPI calls negative with low confidence**, and those negatives
  are not trustworthy: every genuine false negative on the control run sat below
  **0.55** — CD14 0.388, FCGR3B 0.445, ART4 0.534, NYX 0.549 — and CD14 and FCGR3B are
  textbook GPI-anchored proteins.

So the rule, taken from ``netgpi_verdicts.py``'s calibration rather than invented here:

    anchored = NetGPI says anchored in ANY form
               OR curation says anchored and NetGPI's negative is below CONFIDENT (0.55)

A NetGPI negative is only acted on **above** CONFIDENT. Below it the prediction is treated
as uninformative and curation stands.

Emits every gene with its evidence and score, so a later reader can see which source
carried it and at what confidence rather than having to re-derive this.

    uv run --locked python build_gpi_set.py
"""

from __future__ import annotations

import argparse
import collections
import csv
import re
from pathlib import Path

HERE = Path(__file__).parent
SCRATCH = HERE.parents[2] / "scratch/becca/appris-principal"
NETGPI_HUMAN = SCRATCH / "netgpi_gpi_human_predictions.txt"
CURATED = HERE / "gpi_anchored_genes.tsv"

#: Above this, NetGPI's negatives were right on the control set; below it, every genuine
#: false negative appeared. Read off the control run's failures, not tuned. Same constant
#: as ``netgpi_verdicts.CONFIDENT`` -- keep them equal.
CONFIDENT = 0.55

_ID = re.compile(r"^(ref|alt)_(.+?)_(ENST\d+|[A-Z0-9]+(?:-\d+)?)$")


def read_netgpi(path: Path) -> list[dict]:
    """NetGPI short-output rows. ``id`` is ``{ref|alt}_{GENE}_{acc-or-transcript}``."""
    out = []
    for line in path.read_text().splitlines():
        if line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 6:
            continue
        m = _ID.match(f[0])
        if not m:
            continue
        out.append(
            {
                "form": m.group(1),
                "gene": m.group(2),
                "ident": m.group(3),
                "anchored": f[2].strip().startswith("GPI"),
                "omega": f[3].strip(),
                "likelihood": float(f[4]),
                "omega_aa": f[5].strip(),
            }
        )
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--netgpi", type=Path, default=NETGPI_HUMAN)
    ap.add_argument("--curated", type=Path, default=CURATED)
    ap.add_argument("--out", type=Path, default=HERE / "gpi_set.tsv")
    args = ap.parse_args(argv)

    preds = read_netgpi(args.netgpi)
    by_gene: dict[str, list[dict]] = collections.defaultdict(list)
    for p in preds:
        by_gene[p["gene"]].append(p)
    with args.curated.open() as fh:
        curated = {r["gene_symbol"].strip() for r in csv.DictReader(fh, delimiter="\t")}

    print(
        f"NetGPI forms: {len(preds)} over {len(by_gene)} genes "
        f"({sum(1 for p in preds if p['form'] == 'ref')} reference, "
        f"{sum(1 for p in preds if p['form'] == 'alt')} alternative)"
    )
    print(f"curated (triage gpi_anchored): {len(curated)}")

    rows = []
    for gene in sorted(set(by_gene) | curated):
        forms = by_gene.get(gene, [])
        anchored = [f for f in forms if f["anchored"]]
        negatives = [f for f in forms if not f["anchored"]]
        # The CONFIDENT threshold was measured on 148 curated-GPI human **reference**
        # sequences, so it is the REFERENCE form's negative that it applies to. Taking the
        # max across all forms instead reads an alternative form's confident negative as if
        # it were the reference's, and alternative forms are different proteins that are
        # legitimately not anchored: ART4's reference scores 0.534 (a known false negative)
        # while its alternatives score 0.991/0.983, and FCGR3B's reference 0.445 against an
        # alternative's 0.994. Using the max flipped exactly those genes out of the set.
        ref_neg = next((f["likelihood"] for f in negatives if f["form"] == "ref"), None)
        alt_neg_max = max((f["likelihood"] for f in negatives if f["form"] == "alt"), default=None)
        best_neg = ref_neg
        in_curated = gene in curated

        if anchored:
            call, basis = True, "netgpi_anchored"
        elif in_curated and (best_neg is None or best_neg < CONFIDENT):
            # Negative below CONFIDENT is uninformative; curation stands.
            call, basis = True, "curated_netgpi_negative_low_confidence"
        elif in_curated:
            # NetGPI is confidently negative AND curation says anchored. Flagged, not
            # silently resolved -- this is the set a human should look at.
            call, basis = True, "conflict_curated_vs_confident_netgpi_negative"
        else:
            call, basis = False, "netgpi_not_anchored"

        best_anch = max(anchored, key=lambda f: f["likelihood"], default=None)
        rows.append(
            {
                "gene_symbol": gene,
                "gpi_anchored": int(call),
                "basis": basis,
                "netgpi_anchored_forms": len(anchored),
                "netgpi_total_forms": len(forms),
                "anchored_in_reference_form": int(any(f["form"] == "ref" for f in anchored)),
                "anchored_only_in_alternative_form": int(
                    bool(anchored) and not any(f["form"] == "ref" for f in anchored)
                ),
                "best_anchored_likelihood": f"{best_anch['likelihood']:.3f}" if best_anch else "",
                "omega_site": best_anch["omega"] if best_anch else "",
                "omega_aa": best_anch["omega_aa"] if best_anch else "",
                # The reference form's negative -- the one CONFIDENT is calibrated on.
                "reference_negative_likelihood": f"{ref_neg:.3f}" if ref_neg is not None else "",
                # Informational only: never decides the call, for the reason above.
                "strongest_alt_negative_likelihood": (
                    f"{alt_neg_max:.3f}" if alt_neg_max is not None else ""
                ),
                "in_curated_set": int(in_curated),
            }
        )

    with args.out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    anchored_rows = [r for r in rows if r["gpi_anchored"]]
    print(f"\nwrote {args.out.name}: {len(rows)} genes scored, {len(anchored_rows)} GPI-anchored")
    print(f"\n{'basis':52s} {'genes':>6}")
    for b, n in collections.Counter(r["basis"] for r in rows).most_common():
        print(f"{b:52s} {n:>6}")

    alt_only = [r for r in anchored_rows if r["anchored_only_in_alternative_form"]]
    print(
        f"\nanchored ONLY in a non-reference form ({len(alt_only)}) -- the cases a "
        f"canonical-only annotation cannot see:"
    )
    print("  " + ", ".join(r["gene_symbol"] for r in alt_only[:20]))
    conflicts = [r for r in rows if r["basis"].startswith("conflict")]
    if conflicts:
        print(
            f"\nCONFLICTS for human review ({len(conflicts)}): curation says anchored, "
            f"NetGPI confidently disagrees"
        )
        for r in conflicts:
            print(
                f"  {r['gene_symbol']:12s} reference-form negative {r['reference_negative_likelihood']}"
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
